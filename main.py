import tkinter as tk
from tkinter import filedialog, ttk
from PIL import Image, ImageTk
from pathlib import Path

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

        self.canvas.bind("<MouseWheel>", self._on_mousewheel)

        # --- State ---
        self.images = []
        self.idx = -1
        self.ti_img = None


    # --- Functions ---
    def open_image(self):
        path = filedialog.askopenfilename(
            filetypes=[("Image", "*.png;*.jpg;*.jpeg;*.tif;*.bmp")]
            )
        if not path:
            return
        self.build_list(Path(path))
        self.idex = self.images.index(Path(path))
        self.display()

    def build_list(self, first_path: Path):
        folder = first_path.parent
        exts = {".png", ".jpg", ".jpeg", ".tif", ".bmp"}
        self.images = [p for p in sorted(folder.iterdir()) if p.suffix.lower() in exts]

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

    def display(self):
        path = self.images[self.idx]
        img = Image.open(path)
        self.tk_img = ImageTk.PhotoImage(img)
        self.canvas.delete("all")
        self.canvas.create_image(0, 0, anchor="nw", image=self.tk_img)
        self.canvas.config(scrollregion=self.canvas.bbox("all"))
        self.path_var.set(f"{path.name} ({self.idx}/{len(self.images)-1})")
        self.index_var.set(self.idx)
    
    def _on_mousewheel(self, event):
        self.canvas.yview_scroll(int(-1 * (event.delta / 120)), "units")

if __name__ == "__main__":
    ImageViewer().mainloop()