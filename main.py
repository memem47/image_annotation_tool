import tkinter as tk
from tkinter import filedialog, ttk
from PIL import Image, ImageTk
from pathlib import Path


MAX_EDGE, MAX_PIXELS = 8000, 50e6

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
        self.canvas.bind("<ButtonPress-1>", self._on_mouse_press)
        self.canvas.bind("<B1-Motion>", self._on_mouse_drag)
        self.canvas.bind("<MouseWheel>", self._on_wheel_vertical)
        self.canvas.bind("<Shift-MouseWheel>", self._on_wheel_horizontal)
        self.canvas.bind("<Control-MouseWheel>", self._on_ctrl_wheel)

        # --- State ---
        self.images, self.idx = [], -1
        self.tk_img = None
        self.img_id = None
        self.orig_img = None
        self.orig_w = self.orig_h = 0
        self.current_scale = 1.0

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
                        MAX_PIXELS / (self.orig_w * self.orig_h) ** 0.5)
            w, h = int(self.orig_w * scale), int(self.orig_h * scale)
            self.scale_var.set(round(scale, 2))
        
        # 2. image generation
        self.current_scale = scale
        img = self.orig_img.resize((w, h), Image.LANCZOS)
        self.tk_img = ImageTk.PhotoImage(img)
        self.canvas.itemconfig(self.img_id, image=self.tk_img)
        self.canvas.config(scrollregion=self.canvas.bbox(self.img_id))

        # 3. pan based cursor position
        if cursor:
            cx = self.canvas.canvasx(cursor[0])
            cy = self.canvas.canvasy(cursor[1])
            self.canvas.scale("all", cx, cy, factor:=scale / self.current_scale, factor)
    
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
        
if __name__ == "__main__":
    ImageViewer().mainloop()