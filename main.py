import tkinter as tk
from tkinter import filedialog, ttk
from PIL import Image, ImageTk
from pathlib import Path
from collections import defaultdict


MAX_EDGE, MAX_PIXELS = 8000, 50e6
LINE_COLOR = "#00ffff"
FILL_STIPPLE = "gray50"
POINT_COLOR = "#ff3333"
CTRL_R = 3

class ImageViewer(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Image Annotation Tool")

        # --- UI ---
        toolbar = ttk.Frame(self)
        toolbar.grid(row=0, column=0, columnspan=2, sticky="ew")

        ttk.Button(toolbar, text="Open", command=self.open_image).pack(side="left")
        ttk.Button(toolbar, text="Prev", command=lambda: self.show_image(-1)).pack(side="left")
        ttk.Button(toolbar, text="Next", command=lambda: self.show_image(+1)).pack(side="left")

        self.index_var = tk.StringVar()
        ttk.Entry(toolbar, width=5, textvariable=self.index_var).pack(side="left")
        ttk.Button(toolbar, text="Go", command=self.go_index).pack(side="left")

        self.path_var = tk.StringVar(value="No image loaded")
        ttk.Label(toolbar, textvariable=self.path_var).pack(side="left", padx=10)

        # --- Zoom slider ---
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

        # Scrollable canvas
        self.canvas = tk.Canvas(self, bg="grey80")
        self.hsb = tk.Scrollbar(self, orient="horizontal", command=self.canvas.xview)
        self.vsb = tk.Scrollbar(self, orient="vertical", command=self.canvas.yview)
        self.canvas.configure(xscrollcommand=self.hsb.set, yscrollcommand=self.vsb.set)

        self.canvas.grid(row=1, column=0, sticky="nsew")
        self.vsb.grid(row=1, column=1, sticky="ns")
        self.hsb.grid(row=2, column=0, sticky="ew")

        self.rowconfigure(1, weight=1)
        self.columnconfigure(0, weight=1)

        # --- Mouse bindings ---
        self.canvas.bind("<MouseWheel>", self._on_wheel_vertical)
        self.canvas.bind("<Shift-MouseWheel>", self._on_wheel_horizontal)
        
        self.canvas.bind("<Control-MouseWheel>", self._on_ctrl_wheel)
        self.canvas.bind("<Motion>", self._on_motion)
        
        # center button drag
        self.canvas.bind("<ButtonPress-2>", self._on_mouse_press)
        self.canvas.bind("<B2-Motion>", self._on_mouse_drag)

        # Left click = add vertex
        self.canvas.bind("<ButtonPress-1>", self._on_left_down)
        self.canvas.bind("<B1-Motion>", self._on_left_drag)
        self.canvas.bind("<ButtonRelease-1>", self._on_left_up)
        # Esc or Right click = Cancel
        self.bind_all("<Escape>", self._on_cancel_poly)
        self.canvas.bind("<ButtonPress-3>", self._on_cancel_poly)

        # --- State ---
        self.images, self.idx = [], -1
        self.tk_img = None
        self.img_id = None
        self.orig_img = None
        self.orig_w = self.orig_h = 0
        self.current_scale = 1.0
        self.annotations = defaultdict(list)
        self.curr_poly  = [] # image coordinate
        self.poly_items = [] # canvas item ID
        self.ctrl_ids   = [] # oval
        self.drag_idx = None # edditing vertex index
        self.polygons = []
        self.active_id = None
        

        # --- Status bar ---
        self.status_var = tk.StringVar()
        ttk.Label(self, textvariable=self.status_var, anchor="w").grid(row=3, column=0, columnspan=2, sticky="ew")

    # --- File handling ---
    def open_image(self):
        path = filedialog.askopenfilename(
            filetypes=[("Image", "*.png;*.jpg;*.jpeg;*.tif;*.bmp")]
            )
        if not path:
            return
        self.build_list(Path(path))
        self.idx = self.images.index(Path(path))
        self.display()

    def build_list(self, first_path: Path):
        folder = first_path.parent
        exts = {".png", ".jpg", ".jpeg", ".tif", ".bmp"}
        self.images = [p for p in sorted(folder.iterdir()) if p.suffix.lower() in exts]

    # --- Navigation --- 
    def show_image(self, delta):
        if not self.images: return
        self.idx = (self.idx + delta) % len(self.images)
        self.display()
        
    def go_index(self):
        if not self.images: return
        try:
            i = int(self.index_var.get())
            if 0 <= i < len(self.images):
                self.idx = i
                self.display()
        except ValueError:
            pass # ignore invalid input

    # --- Display ---
    def display(self):
        path = self.images[self.idx]
        self.orig_img = Image.open(path)
        self.orig_w, self.orig_h = self.orig_img.size
        self.current_scale = 1.0
        self.scale_var.set(1.0)

        self.tk_img = ImageTk.PhotoImage(self.orig_img)
        if self.img_id is None:
            self.img_id = self.canvas.create_image(0,0,anchor="nw", image=self.tk_img)
        else:
            self.canvas.itemconfig(self.img_id, image=self.tk_img)

        self.canvas.config(scrollregion=self.canvas.bbox(self.img_id))
        self.path_var.set(f"{path.name} ({self.idx}/{len(self.images)-1})")
        self.index_var.set(self.idx)
    

    def _on_scale_change(self, _):
        self._apply_zoom(float(self.scale_var.get()))

    def _on_ctrl_wheel(self, event):
        factor = 1.1 if event.delta > 0 else 0.9
        new_scale = min(4.0, max(0.25, self.scale_var.get() * factor))
        self.scale_var.set(round(new_scale, 2))
        self._apply_zoom(new_scale, cursor=(event.x, event.y))
    

    def _apply_zoom(self, scale, cursor=None):
        if self.orig_img is None:
            return
        
        # 1. check limited upper
        w, h = int(self.orig_w * scale), int(self.orig_h * scale)
        if w > MAX_EDGE or h > MAX_EDGE or w * h > MAX_PIXELS:
            scale = min(MAX_EDGE / self.orig_w,
                        MAX_EDGE / self.orig_h,
                        (MAX_PIXELS / (self.orig_w * self.orig_h)) ** 0.5)
            w, h = int(self.orig_w * scale), int(self.orig_h * scale)
            self.scale_var.set(round(scale, 2))
        
        # 2. image generation
        old_scale = self.current_scale
        
        img = self.orig_img.resize((w, h), Image.LANCZOS)
        self.tk_img = ImageTk.PhotoImage(img)
        self.canvas.itemconfig(self.img_id, image=self.tk_img)
        self.canvas.config(scrollregion=self.canvas.bbox(self.img_id))

        # 3. pan based cursor position
        if cursor:
            factor = scale / old_scale
            cx = self.canvas.canvasx(cursor[0])
            cy = self.canvas.canvasy(cursor[1])
            self.canvas.scale("all", cx, cy, factor, factor)
            self.canvas.update_idletasks()
    
        self.current_scale = scale
        self._refresh_ctrl_sizes()
        self._update_scale_label()

    def _update_scale_label(self):
        self.scale_label["text"] = f"{int(self.scale_var.get()*100)} %"
    # --- Pan / Scroll ---
    def _on_mouse_press(self, event):
        self.canvas.scan_mark(event.x, event.y)

    def _on_mouse_drag(self, event):
        self.canvas.scan_dragto(event.x, event.y, gain=1)

    def _on_wheel_vertical(self, event):
        self.canvas.yview_scroll(int(-event.delta / 120), "units")
        
    def _on_wheel_horizontal(self, event):
        self.canvas.xview_scroll(int(-event.delta / 120), "units")
        
    def _on_motion(self, event):
        if self.orig_img is None:
            return
        # world
        wx = self.canvas.canvasx(event.x)
        wy = self.canvas.canvasy(event.y)
        # image anchor
        img_x0, img_y0 = self.canvas.coords(self.img_id)
        # cordinate in image
        ix = int((wx - img_x0) / self.current_scale)
        iy = int((wy - img_y0) / self.current_scale)
        if 0 <= ix < self.orig_w and 0 <= iy < self.orig_h:
            pix = self.orig_img.getpixel((ix, iy))
            self.status_var.set(f"({ix}, {iy}) {pix}")
        else:
            self.status_var.set("")



    def _on_add_vertex(self, event):
        # 1) image coordinate system
        img_x0, img_y0 = self.canvas.coords(self.img_id)
        x_img = (self.canvas.canvasx(event.x) - img_x0) / self.current_scale
        y_img = (self.canvas.canvasy(event.y) - img_y0) / self.current_scale
        self.curr_poly.append((x_img, y_img))

        # 2) canvas coordinate system
        if len(self.curr_poly) > 1:
            x1, y1 = self.curr_poly[-2]
            x2, y2 = self.curr_poly[-1]
            x1c = img_x0 + x1 * self.current_scale
            y1c = img_y0 + y1 * self.current_scale
            x2c = img_x0 + x2 * self.current_scale
            y2c = img_y0 + y2 * self.current_scale,
            line_id = self.canvas.create_line(x1c, y1c, x2c, y2c,
                                              fill=LINE_COLOR, width=2)
            self.poly_items.append(line_id)
        
        # ctrl marker (radius = 3px)
        r = CTRL_R
        cx = img_x0 + x_img * self.current_scale
        cy = img_y0 + y_img * self.current_scale
        cid = self.canvas.create_oval(cx-r, cy-r, cx+r, cy+r,
                                      fill=POINT_COLOR, outline="" )
        self.ctrl_ids.append(cid)
    
    def _on_left_down(self, event):
        img_x0, img_y0 = self.canvas.coords(self.img_id)
        # canvas coordinate
        cx, cy = self.canvas.canvasx(event.x), self.canvas.canvasy(event.y)
        
        # 1) inputted ctrl point hit test (all polygon)
        for pid, poly in enumerate(self.polygons):
            if poly['ctrl']:
                nearest = self.canvas.find_closest(cx, cy)[0]
                if nearest in poly['ctrl']:
                    self.active_id = pid
                    self.drag_idx = poly['ctrl'].index(nearest)
                    return # continue _on_left_drag
        # 2) add vertex to new/conventional polygon
        if self.active_id is None:
            self.polygons.append({'v':[], 'ctrl':[], 'edge':[], 'fill':None})           
            self.active_id = len(self.polygons) - 1
            
        self._add_vertex_to_active(event, img_x0, img_y0)

    def _add_vertex_to_active(self, event, img_x0, img_y0):
        poly = self.polygons[self.active_id]
        # image coordinate
        x_img = (self.canvas.canvasx(event.x) - img_x0)/self.current_scale
        y_img = (self.canvas.canvasy(event.y) - img_y0)/self.current_scale
        poly['v'].append((x_img, y_img))

        # canvas edge
        if len(poly['v']) > 1:
            x1, y1 = poly['v'][-2]; x2, y2 = poly['v'][-1]
            lid = self.canvas.create_line(img_x0 + x1 * self.current_scale,
                                          img_y0 + y1 * self.current_scale,
                                          img_x0 + x2 * self.current_scale,
                                          img_y0 + y2 * self.current_scale,
                                          fill=LINE_COLOR, width=2)
            poly['edge'].append(lid)
        r, cx, cy = CTRL_R, img_x0+x_img*self.current_scale, img_y0+y_img*self.current_scale
        cid = self.canvas.create_oval(cx-r, cy-r, cx+r, cy+r, fill=POINT_COLOR, outline="")
        poly['ctrl'].append(cid)

        # click new a head point auto close
        if len(poly['v']) >= 3:
            x0, y0 = poly['v'][0]
            if (abs(x_img - x0) + abs(y_img - y0)) * self.current_scale < 10:
                self._close_active_polygon(img_x0, img_y0)

    def _close_active_polygon(self, img_x0, img_y0):
        poly = self.polygons[self.active_id]

        x1, y1 = poly['v'][-1]
        x0, y0 = poly['v'][0]
        poly['edge'].append(
            self.canvas.create_line(
                img_x0+x1*self.current_scale,
                img_y0+y1*self.current_scale,
                img_x0+x0*self.current_scale,
                img_y0+y0*self.current_scale,
                fill=LINE_COLOR,width=2
            )
        )
        # fill
        scaled = []
        for x, y in poly['v']:
            scaled.extend([img_x0 + x*self.current_scale,
                           img_y0 + y*self.current_scale])
        poly['fill'] = self.canvas.create_polygon(
            *scaled, fill=LINE_COLOR, outline=LINE_COLOR, width=2,
            stipple=FILL_STIPPLE
        )

        self.active_id = None

    def _on_left_drag(self, event):
        if self.active_id is None or self.drag_idx is None:
            return

        poly = self.polygons[self.active_id]
        img_x0, img_y0 = self.canvas.coords(self.img_id)
        x_img = (self.canvas.canvasx(event.x) - img_x0) / self.current_scale
        y_img = (self.canvas.canvasy(event.y) - img_y0) / self.current_scale
        poly['v'][self.drag_idx] = (x_img, y_img)
        # 2) update marker position
        r = CTRL_R
        cx = img_x0 + x_img * self.current_scale
        cy = img_y0 + y_img * self.current_scale
        cid = poly['ctrl'][self.drag_idx]
        self.canvas.coords(cid, cx-r, cy-r, cx+r, cy+r)
        # 3) re-draw
        for lid in poly['edge']:
            self.canvas.delete(lid)
        poly['edge'].clear()
        for i in range(1, len(poly['v'])):
            x1, y1 = poly['v'][i-1]; x2, y2 = poly['v'][i]
            poly['edge'].append(
                self.canvas.create_line(img_x0+x1*self.current_scale,
                                        img_y0+y1*self.current_scale,
                                        img_x0+x2*self.current_scale,
                                        img_y0+y2*self.current_scale,
                                        fill=LINE_COLOR, width=2)
            )
        # closed loop
        if poly['fill']:
            self.canvas.delete(poly['fill'])
            poly['edge'].append(self.canvas.create_line(
                img_x0+poly['v'][-1][0]*self.current_scale,
                img_y0+poly['v'][-1][1]*self.current_scale,
                img_x0+poly['v'][0][0]*self.current_scale,
                img_y0+poly['v'][0][1]*self.current_scale,
                fill=LINE_COLOR, width=2
            ))
        scaled = []
        for x, y in poly['v']:
            scaled.extend([img_x0 + x*self.current_scale,
                           img_y0 + y*self.current_scale])
            poly['fill'] = self.canvas.create_polygon(
                *scaled, fill=LINE_COLOR, outline=LINE_COLOR, 
                width=2, stipple=FILL_STIPPLE
            )

    def _on_left_up(self, event):
        self.drag_idx = None

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
        if self.active_id is None: return
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
    
    def _on_mode_change(self, *_):
        if self.mode_var.get() == "Annotate":
            self._bind_annotate_mode()
        else:
            self._bind_pan_mode()

    def _refresh_ctrl_sizes(self):
        img_x0, img_y0 = self.canvas.coords(self.img_id)
        for poly in self.polygons:
            for (x, y), cid in zip(poly['v'], poly['ctrl']):
                cx = img_x0 + x * self.current_scale
                cy = img_y0 + y * self.current_scale
                self.canvas.coords(cid, cx-CTRL_R, cy-CTRL_R,
                                        cx+CTRL_R, cy+CTRL_R)
if __name__ == "__main__":
    ImageViewer().mainloop()