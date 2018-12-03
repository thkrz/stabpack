import h5py
import numpy as np


class SPI:
    def __init__(self, name, mode='a'):
        self.f = h5py.File(name, mode)

    def __getitem__(self, key):
        pass
