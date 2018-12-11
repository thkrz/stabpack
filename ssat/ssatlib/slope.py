import numpy as np

from ssat.ssatlib.constants import G


@dataclass
class Facies(MultiPolygon):
    c: float
    phi: float
    rho_b: float
    a: (float, float)
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
    def __getitem__(self, key):
        if type(key) != int:
            key = self.__names.index(key)
        return self.__facieses[key]

    def __init__(self, relief):
        self.__relief = relief
        self.__facieses = []
        self.__names = []

    def __iter__(self):
        self.__curr = 0
        return self

    def __len__(self):
        last = self.__relief.shape[0] - 1
        return self.__relief[last][0] - self.__relief[0][0]

    def __next__(self):
        keys = self.__facieses.keys()
        i = self.__curr
        if i == len(keys):
            raise StopIteration
        self.__curr += 1
        return self.__facieses[keys[i]]

    def __setitem__(self, key, value):
        self.__names.append(key)
        self.__facieses.append(value)

    def height(self, x):
        return np.interp(x, self.__relief[:, 0], self.__relief[:, 1])


def load(name):
    s = Slope()
