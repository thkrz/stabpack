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


class Stratum:
    def __init__(self):
        self.bottom = None

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

        self.ridge = shape.LineString.fromarray(self.ridge)
        if self.cracks is not None:
            self.cracks = np.atleast_1d(self.cracks)

        alpha = np.radians(self.alpha)
        a = shape.Point(self.o, self.ridge.interp(self.o))
        for s in self.__strata:
            a[1] -= s.h / np.cos(alpha)
            s.bottom = shape.Line.fromrotation(a, alpha)

    def __iter__(self):
        self.__iter = iter(self.__strata)
        return self

    def __next__(self):
        return next(self.__iter)

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
                    map(lambda n: setattr(o, n, None), w.values())
                elif ln == "ENDSLOPE":
                    break
                elif ln == "BEGINLAYER":
                    o = Stratum()
                    w = KEYWORDS[2]
                    map(lambda n: setattr(o, n, None), w.values())
                    self.__strata.append(o)
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

    # def kinetic_energy(self, mu, forward=True, m=1.):
    #     threshold = .5 * np.pi
    #     v = np.zeros(self.n)
    #     for i in range(self.n - 1):
    #         j = i + 1
    #         x0 = self[i if forward else self.n - 1 - i]
    #         x1 = self[j if forward else self.n - 1 - j]
    #         dy = x0[1] - x1[1]
    #         if not forward:
    #             dy = -dy
    #         alpha = np.arctan(dy / np.abs(x1[0] - x0[0]))
    #         a = G * (np.sin(alpha) - np.cos(alpha) * mu)
    #         if alpha > threshold or a / G > .6:
    #             v[j] = np.sqrt(2. * G * np.abs(dy)) + v[i] * np.sin(alpha)
    #             continue
    #         if a == 0:
    #             v[j] = v[i]
    #             continue
    #         s = np.linalg.norm(x1 - x0)
    #         p = v[i] / a
    #         q = -2. * s / a
    #         det = p**2 - q
    #         if det >= 0:
    #             t1 = -p + np.sqrt(det)
    #             t2 = -p - np.sqrt(det)
    #             if t1 < 0 and t2 > 0:
    #                 t = t2
    #             elif t2 < 0 and t1 > 0:
    #                 t = t1
    #             elif t1 > 0 and t2 > 0:
    #                 t = np.minimum(t1, t2)
    #             else:
    #                 t = 0.0
    #             v[j] = t * a + v[i]
    #     return .5 * m * np.pow(v, 2)
