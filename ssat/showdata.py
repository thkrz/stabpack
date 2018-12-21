"""Show and alter data parameters"""
import h5py
import numpy as np

import gi
gi.require_version('Gtk', '3.0')
from gi.repository import Gtk
from matplotlib.backends.backend_gtk3agg import (FigureCanvasGTK3Agg as
                                                 FigureCanvas)
from matplotlib.figure import Figure

from ssat import *
from ssat.ssatlib.gui import Dialog, Separator


class Window(Dialog):
    def __init__(self):
        super().__init__("Data")
        self.notebook = Gtk.Notebook()
        self.add_viewport(self.notebook)
        self.notebook.append_page(self.init_design(), Gtk.Label('Design'))
        self.notebook.append_page(self.init_model(), Gtk.Label('Model'))

    def init_design(self):
        box = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=self.spacing)
        box.pack_start(Separator("Test"), False, True, 0)
        box.pack_start(Gtk.Box(), True, True, 0)
        return box

    def init_model(self):
        f = Figure(figsize=(5, 4), dpi=100)
        a = f.add_subplot(111)
        t = np.arange(0.0, 3.0, 0.01)
        s = np.sin(2 * np.pi * t)
        a.plot(t, s)

        sw = Gtk.ScrolledWindow()
        sw.set_border_width(10)
        canvas = FigureCanvas(f)  # a Gtk.DrawingArea
        canvas.set_size_request(800, 600)
        sw.add_with_viewport(canvas)
        return sw

    def on_commit(self, widget):
        print("Commit")


def main():
    win = Window()
    win.connect("destroy", Gtk.main_quit)
    win.show_all()
    Gtk.main()
    return 1
