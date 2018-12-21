import gi
gi.require_version('Gtk', '3.0')
from gi.repository import Gtk


class Dialog(Gtk.Window):
    def __init__(self, title, spacing=8):
        super().__init__(title=title)
        self.spacing = spacing
        self.pane = Gtk.Box(
            orientation=Gtk.Orientation.VERTICAL,
            spacing=spacing,
            margin=spacing)
        self.add(self.pane)

        buttons = Gtk.Box(spacing=spacing)

        button = Gtk.Button(label="Cancel")
        button.connect("clicked", self.on_cancel)
        buttons.pack_end(button, False, False, 0)

        button = Gtk.Button(label="Apply")
        button.connect("clicked", self.on_commit)
        buttons.pack_end(button, False, False, 0)

        button = Gtk.Button(label="OK")
        button.connect("clicked", self.on_ok)
        buttons.pack_end(button, False, False, 0)

        self.pane.pack_end(buttons, False, False, 0)
        self.pane.pack_end(
            Gtk.Separator(orientation=Gtk.Orientation.HORIZONTAL),
            False,
            True,
            0)

    def add_viewport(self, view):
        self.pane.pack_start(view, True, True, 0)

    def on_cancel(self, widget):
        self.destroy()

    def on_commit(self, widget):
        pass

    def on_ok(self, widget):
        self.on_commit(widget)
        self.on_cancel(widget)


class Separator(Gtk.Box):
    def __init__(self, title=None):
        super().__init__()
        if title:
            self.pack_start(Gtk.Label(title), False, False, 0)
        self.pack_start(
            Gtk.Separator(orientation=Gtk.Orientation.HORIZONTAL),
            True,
            True,
            0)
