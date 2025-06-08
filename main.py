import tkinter as tk
from tkinter import filedialog, ttk, messagebox
from PIL import Image, ImageTk, ImageDraw
from pathlib import Path
from collections import defaultdict
import json, itertools, math

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

        # ========================  UI  ====================================
        toolbar = ttk.Frame(self)
        toolbar.grid(row=0, column=0, columnspan=2, sticky="ew")

        ttk.Button(toolbar, text="Open", command=self.open_image).pack(side="left")
        ttk.Button(toolbar, text="Prev", command=lambda: self.show_image(-1)).pack(side="left")
        ttk.Button(toolbar, text="Next", command=lambda: self.show_image(+1)).pack(side="left")
        ttk.Button(toolbar, text="DelAnno", command=self._delete_annotations).pack(side="left")
        ttk.Button(toolbar, text="CopyAnno", command=self._copy_annotations).pack(side="left")
        

        self.index_var = tk.StringVar()
        ttk.Entry(toolbar, width=5, textvariable=self.index_var).pack(side="left")
        ttk.Button(toolbar, text="Go", command=self.go_index).pack(side="left")

        self.path_var = tk.StringVar(value="No image loaded")
        ttk.Label(toolbar, textvariable=self.path_var, width=40, anchor="w").pack(side="left", padx=10)
        self.count_var = tk.StringVar(value="Pt 0  Ed 0  Pg 0")
        ttk.Label(toolbar, textvariable=self.count_var).pack(side="left", padx=(10, 0))
        # ---- Zoom slider --------------------------------------------------
        self.scale_var = tk.DoubleVar(value=1.0)

        scale_frame = ttk.Frame(toolbar)
        scale_frame.pack(side="left", padx=(10,0))
        ttk.Label(scale_frame, text="Scale").pack(side="left")
        self.scale_slider = ttk.Scale(
            scale_frame, from_=0.25, to=4.0,
            orient="horizontal", variable=self.scale_var,
            command=self._on_scale_change, length=120
        )
        self.scale_slider.pack(side="left")
        self.scale_label = ttk.Label(scale_frame, text="100 %")
        self.scale_label.pack(side="left", padx=(2,0))

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
        self.canvas.bind("<ButtonPress-1>", self._on_left_down)
        #self.canvas.bind("<B1-Motion>", self._on_left_drag)
        #self.canvas.bind("<ButtonRelease-1>", self._on_left_up)

        # Esc or Right click = Cancel
        self.bind_all("<Escape>", self._on_cancel_poly)
        self.canvas.bind("<ButtonPress-3>", self._on_cancel_poly)

        # ========================  State  ==================================
        self.images: list[Path] = []
        self.idx: int = -1
        self.tk_img = None
        self.img_id = None
        self.orig_img = None
        self.orig_w = self.orig_h = 0
        self.current_scale = 1.0

        # polygon drawing
        self.polygons = []
        self.ctrl_ids   = [] # oval
        self.poly_items = [] # canvas item ID
        self.drag_idx = None # edditing vertex index
        self.active_id = None
        
        self.annotations = defaultdict(list)
        self.curr_poly  = [] # image coordinate

        # ---- Status bar ---------------------------------------------------
        self.status_var = tk.StringVar()
        ttk.Label(self, textvariable=self.status_var, anchor="w").grid(row=3, column=0, columnspan=2, sticky="ew")

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

    # ------------------------------------------------------------------
    #  Display / redraw helpers
    # ------------------------------------------------------------------
    def display(self):
        path = self.images[self.idx]
        if path.suffix.lower() in (RAW_EXT, BIN_EXT):
            fsize = path.stat().st_size
            widths = self._possible_raw_widths(fsize)
            self.raw_width_cb["values"] = widths
            self.raw_width_cb.configure(state="readonly")
            
            def _best_width(cands):
                s = math.isqrt(fsize)
                return min(cands, key=lambda w: abs(w - s))
            W = (
                self.raw_width_var.get()
                if self.raw_width_var.get() in widths
                else _best_width(widths)
            )
            
            self.raw_width_var.set(W)
            H = fsize // W or 1
            with open(path, 'rb') as fp:
                buf = fp.read()
            self.orig_img = Image.frombytes('L', (W, H), buf)
        else:
            self.raw_width_cb.configure(state="disabled")
            self.orig_img = Image.open(path)
        self.orig_w, self.orig_h = self.orig_img.size

        # ==== Canvas image item ==========================================
        self.tk_img = ImageTk.PhotoImage(self.orig_img)
        if self.img_id is None:
            self.img_id = self.canvas.create_image(0,0,anchor="nw", image=self.tk_img)
        else:
            self.canvas.itemconfig(self.img_id, image=self.tk_img)

        self.canvas.config(scrollregion=self.canvas.bbox(self.img_id))
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
        fsize = path.stat().st_size
        if fsize % new_w:
            return  # invalid width
        # 保存していたズーム・パンを保持
        scale = self.current_scale
        x0, y0 = self.canvas.xview()[0], self.canvas.yview()[0]
        # 再ロード
        H = fsize // new_w or 1 
        with open(path, "rb") as fp:
            buf = fp.read()
        self.orig_img = Image.frombytes("L", (new_w, H), buf)
        self.orig_w, self.orig_h = self.orig_img.size
        self.tk_img = ImageTk.PhotoImage(self.orig_img)
        self.canvas.itemconfig(self.img_id, image=self.tk_img)
        self.canvas.config(scrollregion=self.canvas.bbox(self.img_id))
        # re‑scale & pan
        self.current_scale = 1.0
        self._apply_zoom(scale)
        self.canvas.xview_moveto(x0)
        self.canvas.yview_moveto(y0)
        cx, cy = self.canvas.coords(self.img_id)
        self._refresh_ctrl_sizes(cx, cy)


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
        img_x0, img_y0 = self.canvas.coords(self.img_id)
        self._apply_zoom(new_scale, cursor=(img_x0, img_y0))
        #self._apply_zoom(new_scale, cursor=(event.x, event.y))

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
        img_x0, img_y0 = self.canvas.coords(self.img_id)
        # canvas coordinate
        x_canvas = self.canvas.canvasx(event.x)
        y_canvas = self.canvas.canvasy(event.y)
        
        # 1) inputted ctrl point hit test (all polygon)
        # for pid, poly in enumerate(self.polygons):
        #     if poly['ctrl']:
        #         nearest = self.canvas.find_closest(x_canvas, y_canvas)[0]
        #         if nearest in poly['ctrl']:
        #             self.active_id = pid
        #             self.drag_idx = poly['ctrl'].index(nearest)
        #             return # continue _on_l-eft_drag
                
        # 2) add vertex to new/conventional polygon
        if self.active_id is None:
            self.polygons.append({'v':[], 'ctrl':[], 'edge':[], 'fill':None})           
            self.active_id = len(self.polygons) - 1
            
        self._add_vertex_to_active(event, img_x0, img_y0)

    def _add_vertex_to_active(self, event, img_x0, img_y0):
        poly = self.polygons[self.active_id]
        # image coordinate
        x_canvas = self.canvas.canvasx(event.x) 
        y_canvas = self.canvas.canvasy(event.y) 
        x_img = (x_canvas - img_x0)/self.current_scale
        y_img = (y_canvas - img_y0)/self.current_scale

        # click head point -> auto close
        if len(poly['v']) >= 3:
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

    def _on_left_drag(self, event):
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
                self.canvas.create_line(img_x0+x1*self.current_scale,
                                        img_y0+y1*self.current_scale,
                                        img_x0+x2*self.current_scale,
                                        img_y0+y2*self.current_scale,
                                        fill=LINE_COLOR, width=2)
            )
            self.poly_items.append(poly['edge'][-1])

        # closed loop
        if poly['fill']:
            self.canvas.delete(poly['fill'])
            if poly['fill'] in self.poly_items:
                self.poly_items.remove(poly['fill'])
            x1, y1 = poly['v'][-1]
            x0, y0 = poly['v'][0]
            poly['edge'].append(
                self.canvas.create_line(
                    img_x0+x1*self.current_scale,
                    img_y0+y1*self.current_scale,
                    img_x0+x0*self.current_scale,
                    img_y0+y0*self.current_scale,
                    fill=LINE_COLOR, width=2
                )
            )
        scaled = []
        for x, y in poly['v']:
            scaled.extend([img_x0 + x*self.current_scale,
                           img_y0 + y*self.current_scale])
        poly['fill'] = self.canvas.create_polygon(
            *scaled, fill=LINE_COLOR, outline=LINE_COLOR, 
            width=2, stipple=FILL_STIPPLE
        )
        self.poly_items.append(poly['fill'])

    def _on_left_up(self, event):
        if self.drag_idx == 0:
            poly = self.polygons[self.active_id]
            if len(poly['v']) >= 3:
                img_x0, img_y0 = self.canvas.coords(self.img_id)
                self._close_active_polygon(img_x0, img_y0)
        else:
            self.drag_idx = None
            self._autosave_json()
            self._save_mask()

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
        if poly['v']==[] and poly['ctrl']==[] and poly['edge']==[]:
            del self.polygons[self.active_id]
            self.active_id -= 1
            if self.active_id < 0:
                self.active_id = None
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

if __name__ == "__main__":
    ImageViewer().mainloop()