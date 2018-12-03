"""View and modify datafile"""
import h5py

from tkinter import *
from tkinter import ttk

from ssat import *


class Application(ttk.Frame):
    MARGIN = "8"

    def __init__(self, parent):
        super().__init__(parent, padding=self.MARGIN)
        self.parent = parent

    def run(self) -> int:
        self.parent.mainloop()
        return 1


def main() -> int:
    with h5py.File(cfg['DATAFILE'], 'a') as f:
        pass
    root = Tk()
    root.title(cfg['DATAFILE'])
    app = Application(root)
    return app.run()
