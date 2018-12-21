import numpy as np

from typing import Tuple

from ssat.ssatlib.constants import G


@dataclass
class Layer:
    z: float
    c: float
    phi: float
    rho_b: float
    a: Tuple[float, float]
    theta: float
    Ks: float
    E: float
    pr: float
    rho_p: float
    n: float
    weights: np.array

    def gamma(self):
        return self.rho_b * G


class Slope:
    def __getitem__(self, i):
        if type(i) == int:
            return self.__layers[i]
        if type(i) == str:
            return self.__layers[self.__names.index(i)]
        if type(i) == Point:
            y = self.boundary(i.x)
            for l in self.__layers:
                if i.y > y:
                    return l
                y -= l.z / np.cos(self.dip)

    def __init__(self, relief, dip=0, outcrop=0):
        self.__x = relief[:, 0]
        self.__y = relief[:, 1]
        self.__layers = []
        self.dip = dip
        if outcrop < self.__x[0]:
            self.b = self.__y[0] # Needed? maybe interp takes care of it
        else:
            self.b = self.height(outcrop)

    def boundary(self, x):
        return self.b + np.tan(np.radians(self.dip)) * x

    def height(self, x):
        return np.interp(x, self.__x, self.__y)


def load(name):
    s = Slope()
