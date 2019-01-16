import numpy as np

from itertools import chain

from ssat.ssatlib import shape

KEYWORDS = [{
    "BEGINSLOPE": None,
    "ENDSLOPE": None,
    "BEGINLAYER": None,
    "ENDLAYER": None,
    "COMMENT": None,
}, {
    "CALC_EXT": "span",
    "CRACKS": "cracks",
    "DIP": "alpha",
    "FAULTS": "faults",
    "OUTCROP": "o",
    "PROFILE": "ridge"
}, {
    "BULK_DENSITY": "rho_b",
    "COHESION": "c",
    "FRICTION_ANGLE": "phi",
    "GRAIN_SIZES": "pdf",
    "HYDRAULIC_CONDUCTIVITY": "K",
    "NAME": "name",
    "PARTICLE_DENSITY": "rho_p",
    "POISSONS_RATIO": "pr",
    "POROSITY": "n",
    "RETENTION_MODEL": "beta",
    "THICKNESS": "h",
    "WATER_CONTENT": "theta",
    "YOUNGS_MODULUS": "E"
}]
WORDS = list(chain.from_iterable([w.keys() for w in KEYWORDS]))


class Stratum(shape.Edge):
    def __init__(self):
        super().__init__(0, 0)

    def fissure_depth(self):
        return 2. * self.c / self.gamma() * np.tan(
            np.radians(45. + .5 * self.phi))

    def isdebris(self):
        return self.name == 'DEBRIS'


class Slope:
    def __init__(self, name, mu, dim=2):
        self.dim = dim

        self.__strata = []
        self.__parse(name)

        self.ridge = shape.Boundary(self.ridge[:, 0], self.ridge[:, 1])
        if not hasattr(self.span):
            self.span = self.ridge.span
        if hasattr(self.cracks):
            self.cracks = np.atleast_1d(self.cracks)

        top = self.ridge
        i = 0
        s = self.__strata[i]
        while s.isdebris():
            v = top.vdist(mu)
            top.mirror()
            v += top.vdist(mu)
            top.mirror()
            vmax = np.amax(v)
            y = (s.h[1] - s.h[0]) * (v / vmax) + s.h[0]
            s.bottom = shape.Boundary(top.x, top.y - y)
            top = s.bottom
            i += 1
            s = self.__strata[i]
        b = self.top(self.o)
        m = np.arctan(np.radians(self.alpha))
        for s in self.__strata[i:]:
            s.m = m
            s.b = b - s.h
            b = s.b

    def __iter__(self):
        self.__curr = 0
        return self

    def __next__(self):
        i = self.__curr
        if i == len(self.__strata):
            raise StopIteration
        self.__curr += 1
        return self.__strata[i]

    def __parse(self, name):
        c = False
        eof = False
        k, o = None, None
        v = []
        w = None
        with open(name) as dat:
            ln = next(dat).strip()
            while True:
                try:
                    peek = next(dat).strip()
                except StopIteration:
                    if eof:
                        break
                    eof = True
                if ln == "BEGINSLOPE":
                    o = self
                    w = KEYWORDS[1]
                elif ln == "ENDSLOPE":
                    break
                elif ln == "BEGINLAYER":
                    o = Stratum()
                    self.__strata.append(o)
                    w = KEYWORDS[2]
                elif ln == "ENDLAYER":
                    o = self
                    w = KEYWORDS[1]
                elif ln == "COMMENT":
                    c = True
                elif ln in w.keys():
                    k = ln
                elif ln:
                    if not c:
                        v.append(ln.split())
                    if peek in WORDS:
                        if c:
                            c = False
                        else:
                            if k == "NAME":
                                v = v[0][0]
                            else:
                                v = np.squeeze(np.array(v).astype(np.float))
                            setattr(o, w[k], v)
                            v = []
                ln = peek

    def top(self, x):
        return np.interp(x, self.ridge.x, self.ridge.y)
