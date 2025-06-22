import tkinter as tk
from tkinter import filedialog, ttk, messagebox
from PIL import Image, ImageTk, ImageDraw, ImageFilter
from pathlib import Path
from collections import defaultdict
import json, itertools, math
import numpy as np

# --- Constants ------------------------------------------------------------
MAX_EDGE, MAX_PIXELS = 8000, 50e6
LINE_COLOR = "#00ffff"
FILL_STIPPLE = "gray50"
POINT_COLOR = "#ff3333"
CTRL_R = 3

MASK_SUFFIX = ".mask.png"
ANNO_SUFFIX = ".anno.json"
RAW_EXT = ".raw"
BIN_EXT = ".bin"
RAW_WIDTH_DEF = 1000

class ImageViewer(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Image Annotation Tool")

        # ---------- TOP TOOLBAR ----------------------------------------
        toolbar = ttk.Frame(self, padding=2)
        toolbar.grid(row=0, column=0, columnspan=2, sticky="ew")

        # Left Block: Navigation --------------------------------------
        nav = ttk.Frame(toolbar)
        nav.pack(side="left", padx=(0, 8))
        for text, cmd in (("Open", self.open_image),
                          ("Prev", lambda: self.show_image(-1)),
                          ("Next", lambda: self.show_image(+1))):
            ttk.Button(nav, text=text, command=cmd, width=5).pack(side="left")

        self.index_var = tk.StringVar()
        ttk.Entry(nav, width=5, textvariable=self.index_var).pack(side="left")
        ttk.Button(nav, text="Go", command=self.go_index).pack(side="left")

        self.path_var = tk.StringVar(value="No image loaded")
        ttk.Label(nav, textvariable=self.path_var, width=40, anchor="w").pack(side="left", padx=10)
                

        # ---- RAW width selector ------------------------------------------
        raw_frame = ttk.Frame(toolbar)
        raw_frame.pack(side="left", padx=(10,0))

        ttk.Label(raw_frame, text="RAW W").pack(side="left")
        self.raw_width_var = tk.IntVar(value=RAW_WIDTH_DEF)
        self.raw_width_cb = ttk.Combobox(
            raw_frame, width=8, textvariable=self.raw_width_var,
            state='disabled', values=(RAW_WIDTH_DEF,)
        )
        self.raw_width_cb.pack(side="left")
        self.raw_width_var.trace_add("write", lambda *_: self._on_raw_width_change())
        
        # RAW bit深度選択（uint8 / uint16）
        ttk.Label(raw_frame, text="Type").pack(side="left", padx=(8, 0))
        self.raw_type_var = tk.StringVar(value="uint8")
        self.raw_type_cb = ttk.Combobox(
            raw_frame, width=6, textvariable=self.raw_type_var,
            state="readonly", values=("uint8", "uint16")
        )
        self.raw_type_cb.pack(side="left")

        # コールバックが必要なら以下も追加（例）
        self.raw_type_var.trace_add("write", lambda *_: self._on_raw_width_change())   

        ttk.Separator(nav, orient="vertical").pack(side="left", fill="y", pady=2)

        # Center Block: Annotation tools ----------------------
        anno = ttk.Frame(toolbar)             
        anno.pack(side="left", padx=(8, 8))
        for text, cmd in (("Del", self._delete_annotations),
                          ("Copy", self._copy_annotations)):
            ttk.Button(toolbar, text=text, command=cmd, width=4).pack(side="left")
        
        self.lasso_var = tk.BooleanVar(value=False)
        lasso = ttk.Frame(toolbar)
        lasso.pack(side="left", padx=(8, 0))
        self.lasso_btn = ttk.Button(
            lasso, text="Lasso OFF", width=8, command=self._toggle_lasso)
        self.lasso_btn.pack(side="left")


  
        # --- Lasso parameters -----------------------------------------
        self.tol_var = tk.IntVar(value=12)
        self.eps_var = tk.DoubleVar(value=2.5)
        parm = ttk.Frame(lasso)
        parm.pack(side="left", padx=(6, 0))
        # Tol
        ttk.Label(parm, text="T").pack(side="left")
        tol_scale = ttk.Scale(parm, from_=1, to=50, orient="horizontal",
                              variable=self.tol_var, length=60)
        tol_scale.pack(side="left")
        tol_val = ttk.Label(parm, width=2, anchor="e")
        tol_val.pack(side="left")
        self.tol_var.trace_add("write",
            lambda *_: tol_val.config(text=str(self.tol_var.get())))
        # Eps
        ttk.Label(parm, text="Eps").pack(side="left", padx=(6,0))
        eps_scale = ttk.Scale(parm, from_=0.5, to=10,
                              orient="horizontal", variable=self.eps_var,
                              length=60)
        eps_scale.pack(side="left")
        eps_val = ttk.Label(parm, width=4, anchor="e")
        eps_val.pack(side="left")
        self.eps_var.trace_add("write",
            lambda *_: eps_val.config(text=f"{self.eps_var.get():.1f}"))
        
        # ---- Scrollable canvas -------------------------------------------
        self.canvas = tk.Canvas(self, bg="grey80")
        self.hsb = tk.Scrollbar(self, orient="horizontal", command=self.canvas.xview)
        self.vsb = tk.Scrollbar(self, orient="vertical", command=self.canvas.yview)
        self.canvas.configure(xscrollcommand=self.hsb.set, yscrollcommand=self.vsb.set)

        self.canvas.grid(row=1, column=0, sticky="nsew")
        self.vsb.grid(row=1, column=1, sticky="ns")
        self.hsb.grid(row=2, column=0, sticky="ew")

        self.rowconfigure(1, weight=1)
        self.columnconfigure(0, weight=1)

        # ---- Mouse bindings ----------------------------------------------
        self.canvas.bind("<MouseWheel>", self._on_wheel_vertical)
        self.canvas.bind("<Shift-MouseWheel>", self._on_wheel_horizontal)        
        self.canvas.bind("<Control-MouseWheel>", self._on_ctrl_wheel)
        self.canvas.bind("<Motion>", self._on_motion)
        
        # center button drag
        self.canvas.bind("<ButtonPress-2>", self._on_mouse_press)
        self.canvas.bind("<B2-Motion>", self._on_mouse_drag)

        # Left click = add vertex
        # * plain : add vertex
        # * Shift : edit (drag existing vertex)
        self.canvas.bind("<ButtonPress-1>", self._on_left_down)
        self.canvas.bind("<B1-Motion>", self._on_left_drag)
        self.canvas.bind("<ButtonRelease-1>", self._on_left_up)
        self.canvas.bind("<Double-Button-1>", self._on_double_click)

        # Esc or Right click = Cancel
        self.bind_all("<Escape>", self._on_cancel_poly)
        #   通常:   右クリック = 頂点取消
        #   Lasso: 右クリック = 抽出マスク破棄
        self.canvas.bind("<ButtonPress-3>", self._on_right_click)

        # ========================  State  ==================================
        self.images: list[Path] = []
        self.idx: int = -1
        self.tk_img = None
        self.img_id = None
        self.orig_img = None
        self.orig_w = self.orig_h = 0
        self.current_scale = 1.0
        self._drag_mode = False

        # polygon drawing
        self.polygons = []
        self.ctrl_ids   = [] # oval
        self.poly_items = [] # canvas item ID
        self.drag_idx = None # edditing vertex index
        self.active_id = None
        
        self.annotations = defaultdict(list)
        self.curr_poly  = [] # image coordinate
        self._in_lasso_build = False  # lasso mode

        # ---- Status bar ---------------------------------------------------
        status_frame = ttk.Frame(self)
        status_frame.grid(row=3, column=0, columnspan=2, sticky="ew")

        self.status_var = tk.StringVar()
        ttk.Label(status_frame, textvariable=self.status_var, anchor="w").pack(side="left", fill="x", expand=True)

        self.count_var = tk.StringVar(value="Pt 0  Ed 0  Pg 0")
        ttk.Label(status_frame, textvariable=self.count_var, anchor="w").pack(side="left", padx=(10, 10))
        
        # ---- Zoom slider --------------------------------------------------
        self.scale_var = tk.DoubleVar(value=1.0)
        scale_frame = ttk.Frame(status_frame)
        scale_frame.pack(side="right", padx=(10, 10))
        ttk.Label(scale_frame, text="Scale:").pack(side="left")
        self.scale_slider = ttk.Scale(
            scale_frame, from_=0.25, to=4.0,
            orient="horizontal", variable=self.scale_var,
            command=self._on_scale_change, length=120
        )
        self.scale_slider.pack(side="left")
        self.scale_label = ttk.Label(scale_frame, text="100 %")
        self.scale_label.pack(side="left", padx=(2,0))


    # --- File handling ---
    def open_image(self):
        path = filedialog.askopenfilename(
            filetypes=[("Image", "*.png;*.jpg;*.jpeg;*.tif;*.bmp;*.raw;*.bin")]
            )
        if not path:
            return
        # mask / anno are not selectable
        if path.endswith((MASK_SUFFIX, ANNO_SUFFIX)):
            tk.messagebox.showwarning("Skip", "Mask/anno file can not be directly opened")
            return

        self.build_list(Path(path))
        self.idx = self.images.index(Path(path))
        self.display()

    def build_list(self, first_path: Path):
        folder = first_path.parent
        exts = {".png", ".jpg", ".jpeg", ".tif", ".bmp"}
        self.images = []
        for p in sorted(folder.iterdir()):
            if p.suffix.lower() in exts | {RAW_EXT, BIN_EXT}:
                if not (p.name.endswith(MASK_SUFFIX) or p.name.endswith(ANNO_SUFFIX)):
                    self.images.append(p)

    # ------------------------------------------------------------------
    #  Navigation
    # ------------------------------------------------------------------
    def show_image(self, delta: int):
        if not self.images: 
            return
        self.idx = (self.idx + delta) % len(self.images)
        self.display()
        
    def go_index(self):
        if not self.images:
            return
        try:
            i = int(self.index_var.get())
            if 0 <= i < len(self.images):
                self.idx = i
                self.display()
        except ValueError:
            pass # ignore invalid input

    def read_raw_image(self, fsize, new_w, path):
        depth = self.raw_type_var.get()
        mode = "L" if depth == "uint8" else "I;16L"  # Little-endian 16bit
        bytes_per_pixel = 1 if depth == "uint8" else 2

        # 再ロード
        H = fsize // new_w // bytes_per_pixel or 1 
        with open(path, "rb") as fp:
            buf = fp.read()

        expected_bytes = new_w * H * bytes_per_pixel
        if len(buf) < expected_bytes:
            return  # 安全対策：データ不足

        raw = Image.frombytes(mode, (new_w, H), buf[:expected_bytes])
        if depth =="uint16":
            arr = np.array(raw, dtype=np.uint16)  # PIL → NumPy
            arr = (arr / 256).clip(0, 255).astype(np.uint8)  # 16bit → 8bit
            raw = Image.fromarray(arr, mode="L")  # NumPy → PIL（8bit画像に変換）
        return raw
        
    def open_raw_image(self, path):
        fsize = path.stat().st_size

        # width candidates
        widths = self._possible_raw_widths(fsize)
        self.raw_width_cb["values"] = widths
        self.raw_width_cb.configure(state="readonly")
        
        # decide width
        def _best_width(cands):
            s = math.isqrt(fsize)
            return min(cands, key=lambda w: abs(w - s))
        
        W = (
            self.raw_width_var.get()
            if self.raw_width_var.get() in widths
            else _best_width(widths)
        )
        self.raw_width_var.set(W)

        raw_image = self.read_raw_image(fsize, W, path)
        return raw_image
    
    def display_main(self):        
        self.orig_w, self.orig_h = self.orig_img.size

        # ==== Canvas image item ==========================================
        self.tk_img = ImageTk.PhotoImage(self.orig_img)
        if self.img_id is None:
            self.img_id = self.canvas.create_image(0,0,anchor="nw", image=self.tk_img)
        else:
            self.canvas.itemconfig(self.img_id, image=self.tk_img)

        self.canvas.config(scrollregion=self.canvas.bbox(self.img_id))

    # ------------------------------------------------------------------
    #  Display / redraw helpers
    # ------------------------------------------------------------------
    def display(self):
        path = self.images[self.idx]
        if path.suffix.lower() in (RAW_EXT, BIN_EXT):
            
            self.orig_img = self.open_raw_image(path)
        else:
            self.raw_width_cb.configure(state="disabled")
            self.orig_img = Image.open(path)

        self.display_main()

        self._apply_zoom(self.scale_var.get())

        # ==== Side info ===================================================
        self.path_var.set(f"{path.name} ({self.idx}/{len(self.images)-1})")
        self.index_var.set(self.idx)
    
        # ==== Annotation ==================================================
        self._clear_polygons()
        self._load_annotation_files(path)

        cx, cy = self.canvas.coords(self.img_id)
        self._refresh_ctrl_sizes(cx, cy)

    # ------------------------------------------------------------------
    #  RAW width helpers
    # ------------------------------------------------------------------
    def _possible_raw_widths(self, fsize: int) -> list[int]:
        base = {w for w in range(2, math.isqrt(fsize) + 1) if fsize % w == 0}
        cand  = {w for w in base if w <= MAX_EDGE}
        cand |= {fsize // w for w in base if fsize // w <= MAX_EDGE}
        if not cand:
            cand = {fsize} if fsize <= MAX_EDGE else {RAW_WIDTH_DEF}
        return sorted(cand)
    
    def _on_raw_width_change(self):
        if self.idx < 0 or not self.images:
            return
        path = self.images[self.idx]
        if path.suffix.lower() not in (RAW_EXT, BIN_EXT):
            return  # only for RAW/BIN
        # 幅を取得して再読込
        try:
            new_w = int(self.raw_width_var.get())
        except ValueError:
            return        
        
        # 保存していたズーム・パンを保持
        scale = self.current_scale
        x0, y0 = self.canvas.xview()[0], self.canvas.yview()[0]

        fsize = path.stat().st_size

        if fsize % new_w:
            return  # invalid width
        depth = self.raw_type_var.get()
        mode = "L" if depth == "uint8" else "I;16L"  # Little-endian 16bit
        bytes_per_pixel = 1 if depth == "uint8" else 2

        # 再ロード
        H = fsize // new_w // bytes_per_pixel or 1 
        with open(path, "rb") as fp:
            buf = fp.read()

        expected_bytes = new_w * H * bytes_per_pixel
        if len(buf) < expected_bytes:
            return  # 安全対策：データ不足

        self.orig_img = self.read_raw_image(fsize, new_w, path)

        self.display_main()
        # re‑scale & pan
        self.current_scale = 1.0
        self._apply_zoom(scale)

        self.canvas.xview_moveto(x0)
        self.canvas.yview_moveto(y0)
        # refresh
        cx, cy = self.canvas.coords(self.img_id)
        self._refresh_ctrl_sizes(cx, cy)

    def _on_raw_type_change(self):
        pass

    def _load_annotation_files(self, img_path: Path):        
        anno_path = img_path.with_suffix(ANNO_SUFFIX)
        if not anno_path.exists(): return
        with open(anno_path) as f:
            data = json.load(f)
        img_x0, img_y0 = self.canvas.coords(self.img_id)
        for v in data.get("polygons", []):
            ctrl, edge, fill = [], [], None
            # 頂点と隣接辺
            for i, (x_img, y_img) in enumerate(v):
                x_canvas = img_x0 + x_img * self.current_scale
                y_canvas = img_y0 + y_img * self.current_scale
                cid = self.canvas.create_oval(x_canvas-CTRL_R, 
                                              y_canvas-CTRL_R, 
                                              x_canvas+CTRL_R, 
                                              y_canvas+CTRL_R,
                                              fill=POINT_COLOR, outline="")
                ctrl.append(cid)
                if i:  # 0 番目以外は直前とのエッジ
                    x_prev, y_prev = v[i-1]
                    eid = self.canvas.create_line(img_x0 + x_prev*self.current_scale,
                                                  img_y0 + y_prev*self.current_scale,
                                                  x_canvas, 
                                                  y_canvas, 
                                                  fill=LINE_COLOR, width=2)
                    edge.append(eid)
            # 閉ループおよび塗りつぶしは 1 回だけ作成
            if len(v) >= 3:
                eid = self.canvas.create_line(img_x0 + v[-1][0]*self.current_scale,
                                              img_y0 + v[-1][1]*self.current_scale,
                                              img_x0 + v[0][0]*self.current_scale,
                                              img_y0 + v[0][1]*self.current_scale,
                                              fill=LINE_COLOR, width=2)
                edge.append(eid)
                scaled = [coord*self.current_scale for pt in v for coord in pt]
                fill = self.canvas.create_polygon(*scaled, fill=LINE_COLOR, outline=LINE_COLOR,
                                                  width=2, stipple=FILL_STIPPLE)
            # 登録
            self.polygons.append({'v': v, 'ctrl': ctrl, 'edge': edge, 'fill': fill})
            self.ctrl_ids.extend(ctrl); self.poly_items.extend(edge + ([fill] if fill else []))
        self._update_counts()           # ← ここでカウント表示更新


    # ------------------------------------------------------------------
    #  Zoom / pan ---------------------------------------------------------
    def _on_scale_change(self, _):
        self._apply_zoom(float(self.scale_var.get()))

    def _on_ctrl_wheel(self, event):
        factor = 1.1 if event.delta > 0 else 0.9
        new_scale = min(4.0, max(0.25, self.scale_var.get() * factor))
        self.scale_var.set(round(new_scale, 2))
        self._apply_zoom(float(self.scale_var.get()))

    def _apply_zoom(self, scale: float, cursor: tuple[int, int] | None = None):
        if self.orig_img is None:
            return
        # ---- clamp scale to limits --------------------------------------
        w = int(self.orig_w * scale)
        h = int(self.orig_h * scale)
        if w > MAX_EDGE or h > MAX_EDGE or w * h > MAX_PIXELS:
            scale = min(
                MAX_EDGE / self.orig_w,
                MAX_EDGE / self.orig_h,
                (MAX_PIXELS / (self.orig_w * self.orig_h)) ** 0.5,
            )
            w, h = int(self.orig_w * scale), int(self.orig_h * scale)
            self.scale_var.set(round(scale, 2))

        # ---- update image -----------------------------------------------
        old_scale = self.current_scale
        img_resized = self.orig_img.resize((w, h), Image.LANCZOS)
        self.tk_img = ImageTk.PhotoImage(img_resized)
        self.canvas.itemconfig(self.img_id, image=self.tk_img)
        self.canvas.config(scrollregion=self.canvas.bbox(self.img_id))

        # ---- scale existing drawings ------------------------------------
        factor = scale / old_scale if old_scale else 1.0
        cx, cy = self.canvas.coords(self.img_id)  # 画像左上を基準
        if abs(factor-1.0) > 1e-6:      ### スライダー操作でも図形を追従させる
            if cursor:
                cx = self.canvas.canvasx(cursor[0])
                cy = self.canvas.canvasy(cursor[1])
            else:
                cx, cy = self.canvas.coords(self.img_id)  # 画像左上を基準
            # 画像 (self.img_id) は除外し、注釈アイテムのみ拡大
            for item in itertools.chain(self.ctrl_ids, self.poly_items):
                self.canvas.scale(item, cx, cy, factor, factor)
            self.canvas.update_idletasks()
        self.current_scale = scale
        self._refresh_ctrl_sizes(cx,cy)
        self._update_scale_label()

    def _update_scale_label(self):
        self.scale_label["text"] = f"{int(self.scale_var.get()*100)} %"

    # ---- pan ------------------------------------------------------------
    def _on_mouse_press(self, event):
        self.canvas.scan_mark(event.x, event.y)

    def _on_mouse_drag(self, event):
        self.canvas.scan_dragto(event.x, event.y, gain=1)

    def _on_double_click(self, event):
        """Double-click action depends on Lasso mode."""
        # --- Auto-lasso enabled ----------------------------------------
        if self.lasso_var.get():
            self._auto_lasso(event)        # 既存メソッドを呼び出し
            return
        """Finish current polygon when double-left-clicked anywhere."""
        if self.active_id is None:
            return  # 何も描画していない
        poly = self.polygons[self.active_id]
        if len(poly['v']) < 3 or poly['fill'] is not None:
            return  # 頂点不足 or 既に閉じている
        img_x0, img_y0 = self.canvas.coords(self.img_id)
        self._close_active_polygon(img_x0, img_y0)
        self._update_counts()

    def _on_wheel_vertical(self, event):
        self.canvas.yview_scroll(int(-event.delta / 120), "units")

    def _on_wheel_horizontal(self, event):
        self.canvas.xview_scroll(int(-event.delta / 120), "units")
        
    # ------------------------------------------------------------------
    #  Cursor / status bar
    # ------------------------------------------------------------------
    def _on_motion(self, event):
        if self.orig_img is None:
            return
        wx = self.canvas.canvasx(event.x)
        wy = self.canvas.canvasy(event.y)
        img_x0, img_y0 = self.canvas.coords(self.img_id)
        ix = int((wx - img_x0) / self.current_scale)
        iy = int((wy - img_y0) / self.current_scale)
        if 0 <= ix < self.orig_w and 0 <= iy < self.orig_h:
            pix = self.orig_img.getpixel((ix, iy))
            self.status_var.set(f"({ix}, {iy}) {pix}")
        else:
            self.status_var.set("")

    def _on_left_down(self, event):
        if self.lasso_var.get():
            return
        img_x0, img_y0 = self.canvas.coords(self.img_id)
        # canvas coordinate (cursor position)
        x_canvas = self.canvas.canvasx(event.x)
        y_canvas = self.canvas.canvasy(event.y)
        shift_pressed = bool(event.state & 0x0001)

        # 1) hit-test: did the cursor land on an existing control point?
        if shift_pressed:
            nearest = self.canvas.find_overlapping(
                x_canvas - CTRL_R, y_canvas - CTRL_R,
                x_canvas + CTRL_R, y_canvas + CTRL_R
                )
            for pid, poly in enumerate(self.polygons):
                for idx, cid in enumerate(poly['ctrl']):
                    if cid in nearest:
                        # start drag mode
                        self.active_id = pid
                        self.drag_idx = idx
                        self._drag_mode = True
                        return # continue _on_l-eft_drag
                
        # 2) add vertex to new/conventional polygon
        if self.active_id is None or self.polygons[self.active_id]['fill'] is not None:
            self.polygons.append({'v':[], 'ctrl':[], 'edge':[], 'fill':None})           
            self.active_id = len(self.polygons) - 1
            
        self._add_vertex_to_active(event, img_x0, img_y0)

    def _add_vertex_to_active(self, event, img_x0, img_y0):
        ctrl_pressed = bool(event.state & 0x0004)
        x_canvas = self.canvas.canvasx(event.x) 
        y_canvas = self.canvas.canvasy(event.y) 
        x_img = (x_canvas - img_x0)/self.current_scale
        y_img = (y_canvas - img_y0)/self.current_scale
        # image coordinate
        if ctrl_pressed:
            x_img, y_img = self._insert_edge_snapped_vertex(x_img, y_img)
            x_canvas = x_img*self.current_scale + img_x0
            y_canvas = y_img*self.current_scale + img_y0

        poly = self.polygons[self.active_id]
        # click head point -> auto close
        if not self._in_lasso_build and len(poly['v']) >= 3:
            first_x_canvas = img_x0 + poly['v'][0][0]*self.current_scale
            first_y_canvas = img_y0 + poly['v'][0][1]*self.current_scale
            dist = math.hypot(x_canvas - first_x_canvas,
                              y_canvas - first_y_canvas)
            if dist < CTRL_R*2.5 / self.current_scale:
                self._close_active_polygon(img_x0, img_y0)
                self._update_counts()
                return


        # normal proc
        # draw edge in canvas
        if len(poly['v']) >= 1:
            x1, y1 = poly['v'][-1]
            lid = self.canvas.create_line(img_x0 + x1 * self.current_scale,
                                          img_y0 + y1 * self.current_scale,
                                          x_canvas,
                                          y_canvas,
                                          fill=LINE_COLOR, width=2)
            poly['edge'].append(lid)
            self.poly_items.append(lid)

        # draw ctrl in canvas
        r = CTRL_R
        cid = self.canvas.create_oval(x_canvas - r, 
                                      y_canvas - r, 
                                      x_canvas + r, 
                                      y_canvas + r, 
                                      fill=POINT_COLOR, outline="")
        poly['v'].append((x_img, y_img))
        poly['ctrl'].append(cid)
        self.ctrl_ids.append(cid) 

        self._autosave_json()
        self._update_counts()

    def _close_active_polygon(self, img_x0, img_y0):
        poly = self.polygons[self.active_id]

        x1, y1 = poly['v'][-1]
        x0, y0 = poly['v'][0]
        poly['edge'].append(
            self.canvas.create_line(
                img_x0 + x1*self.current_scale,
                img_y0 + y1*self.current_scale,
                img_x0 + x0*self.current_scale,
                img_y0 + y0*self.current_scale,
                fill=LINE_COLOR,width=2
            )
        )
        self.poly_items.append(poly['edge'][-1]) ### 追加
        # fill
        scaled = []
        for x, y in poly['v']:
            scaled.extend([img_x0 + x*self.current_scale,
                           img_y0 + y*self.current_scale])
        poly['fill'] = self.canvas.create_polygon(
            *scaled, fill=LINE_COLOR, outline=LINE_COLOR, width=2,
            stipple=FILL_STIPPLE
        )
        self.poly_items.append(poly['fill'])
        self.active_id = None

        self._autosave_json()
        self._save_mask()

    def _insert_edge_snapped_vertex(self, x_img, y_img, win=5):
        
        # crop local patch & detect edge magnitude
        x0 = max(x_img - win, 0); x1 = min(x_img + win + 1, self.orig_w)
        y0 = max(y_img - win, 0); y1 = min(y_img + win + 1, self.orig_h)
        patch_gray = self.orig_img.crop((x0, y0, x1, y1))\
                                .convert("L")\
                                .filter(ImageFilter.FIND_EDGES)

        arr = np.asarray(patch_gray, dtype=np.uint8)
        inner = arr[1:-1, 1:-1]          # shape = (h, w)

        dy, dx = np.unravel_index(inner.argmax(), inner.shape)
        dx += 1; dy += 1  # adjust to local patch

        x_img = x0 + dx
        y_img = y0 + dy

        print(x_img, y_img)


        return x_img, y_img

    # ------------------------------------------------------------------
    #  Right-click dispatcher
    # ------------------------------------------------------------------
    def _on_right_click(self, event):
        if self.lasso_var.get():
            # simply clear last generated lasso polygon if存在
            if self.active_id is not None and self.polygons[self.active_id]['fill'] is not None:
                self._on_cancel_poly()        # 既存ロジックで削除
            self.status_var.set("Lasso result cancelled")
        else:
            self._on_cancel_poly()

    def _on_left_drag(self, event):
        if not self._drag_mode:
            return
        if self.active_id is None or self.drag_idx is None:
            return
        poly = self.polygons[self.active_id]
        img_x0, img_y0 = self.canvas.coords(self.img_id)
        x_canvas = self.canvas.canvasx(event.x) 
        y_canvas = self.canvas.canvasy(event.y) 
        x_img = (x_canvas - img_x0)/self.current_scale
        y_img = (y_canvas - img_y0)/self.current_scale
        poly['v'][self.drag_idx] = (x_img, y_img)
        # 2) update marker position
        r = CTRL_R
        cid = poly['ctrl'][self.drag_idx]
        self.canvas.coords(cid, 
                           x_canvas - r, 
                           y_canvas - r, 
                           x_canvas + r, 
                           y_canvas + r)
        # 3) re-draw
        for lid in poly['edge']:
            self.canvas.delete(lid)
            if lid in self.poly_items:
                self.poly_items.remove(lid)
        poly['edge'].clear()
        for i in range(1, len(poly['v'])):
            x1, y1 = poly['v'][i-1]
            x2, y2 = poly['v'][i]
            poly['edge'].append(
                self.canvas.create_line(
                    img_x0+x1*self.current_scale,
                    img_y0+y1*self.current_scale,
                    img_x0+x2*self.current_scale,
                    img_y0+y2*self.current_scale,
                    fill=LINE_COLOR, width=2
                )
            )
            self.poly_items.append(poly['edge'][-1])

        # closed loop
        if poly['fill']:
            self.canvas.delete(poly['fill'])
            if poly['fill'] in self.poly_items:
                self.poly_items.remove(poly['fill'])

            x1, y1 = poly['v'][-1]
            x0, y0 = poly['v'][0]
            eid = self.canvas.create_line(
                img_x0 + x1*self.current_scale,
                img_y0 + y1*self.current_scale,
                img_x0 + x0*self.current_scale,
                img_y0 + y0*self.current_scale,
                fill=LINE_COLOR, width=2
            )
            poly['edge'].append(eid)
            self.poly_items.append(eid)
        scaled = []
        for x, y in poly['v']:
            scaled.extend([img_x0 + x*self.current_scale,
                           img_y0 + y*self.current_scale])
        fid = self.canvas.create_polygon(
            *scaled, fill=LINE_COLOR, outline=LINE_COLOR, 
            width=2, stipple=FILL_STIPPLE
        )
        poly['fill'] = fid
        self.poly_items.append(fid)

    def _on_left_up(self, event):
        if self._drag_mode:
            self._drag_mode = False
            self.drag_idx = None
            self._autosave_json()
            self._save_mask()
            return

    def _on_finish_poly(self, event=None):
        if len(self.curr_poly) < 3:
            return
        # delete unused vertex
        if event is not None and self.curr_poly[-1] == self.curr_poly[-2]:
            self.curr_poly.pop()
        poly = [(int(p[0]), int(p[1])) for p in self.curr_poly]
        self.annotations[self.images[self.idx]].append(poly)
        # display fill region
        scaled = [c*self.current_scale for pt in poly for c in pt]
        pid = self.canvas.create_polygon(*scaled,
                                         fill="#00ffff55", outline=LINE_COLOR, width=2)
        self.poly_items.append(pid)
        self.curr_poly.clear(); self.poly_items.clear()
    
    def _on_cancel_poly(self, event=None):
        # no editting polygon
        if self.active_id is None:
            if len(self.polygons) == 0:
                return
            else:
                self.active_id = len(self.polygons) - 1
        poly = self.polygons[self.active_id]
        # 1) delete last ctrl point and marker
        if poly['ctrl']:
            self.canvas.delete(poly['ctrl'].pop())
        if poly['v']:
            poly['v'].pop()
        # 2) delete last edge
        if poly['edge']:
            self.canvas.delete(poly['edge'].pop())
        if poly['fill']:
            self.canvas.delete(poly['fill'])
            poly['fill'] = None
            if len(poly['v']) < 2:
                self.active_id = None
        if not poly['v'] and not poly['ctrl'] and not poly['edge']:
            del self.polygons[self.active_id]
            self.active_id = None  # これで次回は新規ポリゴンが自動生成される
        self._autosave_json()        
        self._update_counts()
    
    def _on_mode_change(self, *_):
        if self.mode_var.get() == "Annotate":
            self._bind_annotate_mode()
        else:
            self._bind_pan_mode()

    def _refresh_ctrl_sizes(self, img_x0, img_y0):
        #img_x0, img_y0 = self.canvas.coords(self.img_id)
        for poly in self.polygons:
            for (x, y), cid in zip(poly['v'], poly['ctrl']):
                cx = img_x0 + x * self.current_scale
                cy = img_y0 + y * self.current_scale
                self.canvas.coords(cid, cx-CTRL_R, cy-CTRL_R,
                                        cx+CTRL_R, cy+CTRL_R)
                
    def _autosave_json(self):
        if not self.images: return
        img_path = self.images[self.idx]
        out = {'polygons':[pg['v'] for pg in self.polygons]}
        with open(img_path.with_suffix(ANNO_SUFFIX), 'w') as f:
            json.dump(out,f,indent=2)

    def _save_mask(self):
        if not self.images: return
        img_path = self.images[self.idx]
        if not self.polygons: return
        m = Image.new('1', (self.orig_w, self.orig_h))
        draw = ImageDraw.Draw(m)
        for pg in self.polygons:
            if len(pg['v'])>=3:
                draw.polygon(pg['v'],1)
        m.save(img_path.with_suffix(MASK_SUFFIX))

    # ------------------------------------------------------------------
    #  Utility
    # ------------------------------------------------------------------
    def _clear_polygons(self):
        for item in itertools.chain(self.ctrl_ids, self.poly_items):
            self.canvas.delete(item)
        self.polygons.clear()
        self.ctrl_ids.clear()
        self.poly_items.clear()
        self.active_id = self.drag_idx = None
        self._update_counts()

    def _delete_annotations(self):
        if not self.images:
            return
        img_path = self.images[self.idx]
        # clear canvas
        self._clear_polygons()
        # delete existed file
        for suf in (ANNO_SUFFIX, MASK_SUFFIX):
            p = img_path.with_suffix(suf)
            if p.exists():
                p.unlink()
        self._autosave_json() # save empty
        self.status_var.set("Annotations deleted")
        self._update_counts()
        if self.lasso_var.get():
            self._toggle_lasso()

    def _copy_annotations(self):
        if not self.images:
            return
        src_idx = None
        # top priority : next
        if self.idx > 0 and self.images[self.idx-1].with_suffix(ANNO_SUFFIX).exists():
            src_idx = self.idx - 1
        elif self.idx < len(self.images) - 1 and self.images[self.idx + 1].with_suffix(ANNO_SUFFIX).exists():
            src_idx = self.idx + 1
        if src_idx is None:
            messagebox.showinfo("CopyAnno", "No adjacent annotations found.")
            return
        src_path = self.images[src_idx]
        with open(src_path.with_suffix(ANNO_SUFFIX)) as f:
            data = json.load(f)
        polygons = data.get("polygons", [])
        # draw
        self._clear_polygons()
        self._render_polygons(polygons)
        self._autosave_json(); self._save_mask()
        self.status_var.set(f"Copied annotations from {'prev' if src_idx < self.idx else 'next'} image")
        self._update_counts()

    # assist drawing
    def _render_polygons(self, poly_list):
        img_x0, img_y0 = self.canvas.coords(self.img_id)
        for v in poly_list:
            ctrl, edge = [], []
            # 頂点 & 辺
            for i, (x, y) in enumerate(v):
                cx = img_x0 + x * self.current_scale
                cy = img_y0 + y * self.current_scale
                cid = self.canvas.create_oval(cx-CTRL_R, cy-CTRL_R, cx+CTRL_R, cy+CTRL_R,
                                               fill=POINT_COLOR, outline="")
                ctrl.append(cid)
                if i:
                    x_prev, y_prev = v[i-1]
                    eid = self.canvas.create_line(img_x0 + x_prev*self.current_scale,
                                                   img_y0 + y_prev*self.current_scale,
                                                   cx, cy, fill=LINE_COLOR, width=2)
                    edge.append(eid)
            fill = None
            if len(v) >= 3:
                eid = self.canvas.create_line(img_x0 + v[-1][0]*self.current_scale,
                                               img_y0 + v[-1][1]*self.current_scale,
                                               img_x0 + v[0][0]*self.current_scale,
                                               img_y0 + v[0][1]*self.current_scale,
                                               fill=LINE_COLOR, width=2)
                edge.append(eid)
                scaled = [coord*self.current_scale for pt in v for coord in pt]
                fill = self.canvas.create_polygon(*scaled, fill=LINE_COLOR, outline=LINE_COLOR,
                                                   width=2, stipple=FILL_STIPPLE)
            self.polygons.append({'v': v, 'ctrl': ctrl, 'edge': edge, 'fill': fill})
            self.ctrl_ids.extend(ctrl); self.poly_items.extend(edge + ([fill] if fill else []))
        self._update_counts()

    def _update_counts(self):
        """Update count_var label with current Pt / Edge / Poly numbers."""
        n_poly = sum(pg['fill'] is not None for pg in self.polygons)
        n_pts  = sum(len(pg['v'])   for pg in self.polygons)
        n_edg  = sum(len(pg['edge']) for pg in self.polygons)
        self.count_var.set(f"Pt {n_pts}  Ed {n_edg}  Pg {n_poly}")

    # ------------------------------------------------------------------
    #  Toggle Lasso Mode
    # ------------------------------------------------------------------
    def _toggle_lasso(self):
        """Enable / disable auto-lasso mode."""
        state = not self.lasso_var.get()
        self.lasso_var.set(state)
        self.lasso_btn.config(text=f"Lasso {'ON' if state else 'OFF'}")
        self.status_var.set("Lasso mode ON" if state else "Lasso mode OFF")

    def _auto_lasso(self, event):
        """
        Flood-fill from clicked pixel, extract external contour,
        simplify with RDP, then create a draggable polygon.
        """
        tol          = int(self.tol_var.get())
        simplify_eps = float(self.eps_var.get())
        if self.orig_img.mode != "L":
            img_gray = self.orig_img.convert("L")
        else:
            img_gray = self.orig_img

        img_x0, img_y0 = self.canvas.coords(self.img_id)
        cx = int((self.canvas.canvasx(event.x) - img_x0) / self.current_scale)
        cy = int((self.canvas.canvasy(event.y) - img_y0) / self.current_scale)
        if not (0 <= cx < self.orig_w and 0 <= cy < self.orig_h):
            return

        import numpy as np
        import cv2   # needs opencv-python

        arr  = np.array(img_gray, dtype=np.uint8)
        # mask must be (H+2, W+2)
        mask = np.zeros((arr.shape[0] + 2, arr.shape[1] + 2), np.uint8)
        retval, _, _, _ = cv2.floodFill(
            arr.copy(), mask, (cx, cy), 255,
            loDiff=tol, upDiff=tol,
            flags=cv2.FLOODFILL_MASK_ONLY | 8
        )
        # Remove the 1-pixel border that floodFill added
        mask = mask[1:-1, 1:-1]

        if retval <= 20:     # too small
            return
        # find external contour
        contours, _ = cv2.findContours(mask, cv2.RETR_EXTERNAL,
                                       cv2.CHAIN_APPROX_NONE)
        if not contours:
            return
        cnt = max(contours, key=cv2.contourArea).squeeze()
        # Douglas-Peucker simplify
        cnt = cv2.approxPolyDP(cnt, simplify_eps, True).squeeze()
        if cnt.ndim == 1:    # ≒ line
            return
        # register as new polygon
        img_x0, img_y0 = self.canvas.coords(self.img_id)
        self.polygons.append({'v':[], 'ctrl':[], 'edge':[], 'fill':None})
        self.active_id = len(self.polygons) - 1
        self._in_lasso_build = True
        for x, y in cnt:
            evt = tk.Event()
            canvas_x = img_x0 + x * self.current_scale
            canvas_y = img_y0 + y * self.current_scale
            evt.x = canvas_x - self.canvas.canvasx(0)
            evt.y = canvas_y - self.canvas.canvasy(0)
            evt.state = 0          # ← 追加：修飾キー無しを明示
            self._add_vertex_to_active(evt, img_x0, img_y0)
        # 完成扱いにする
        self._close_active_polygon(img_x0, img_y0)
        self._in_lasso_build = False
        self._update_counts()        

if __name__ == "__main__":
    ImageViewer().mainloop()