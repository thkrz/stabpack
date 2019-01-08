import numpy as np

from itertools import chain

from ssat.ssatlib.constants import G

KEYWORDS = [{
    "BEGINSLOPE": None,
    "ENDSLOPE": None,
    "BEGINLAYER": None,
    "ENDLAYER": None,
    "COMMENT": None,
}, {
    "CRACKS": "cracks",
    "DIP": "alpha",
    "OUTCROP": "a",
    "PROFILE": "omega",
    "UNCONFORMITIES": "shifts"
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
    def depth(self):
        return 2. * self.c / self.gamma() * np.tan(
            np.radians(45. + .5 * self.phi))

    def gamma(self):
        return self.rho_b * G


class Slope:
    def __init__(self, name, dim=2, cr=.4):
        self.dim = dim

        self.__strata = []
        self.__parse__(name)

        self.b = self.top(self.a)
        self.m = np.arctan(self.alpha)

        n = self.omega.shape[0] - 1
        self.span = self.omega[n, 0]
        if self.omega[n, 1] > np.mean(self.omega[:, 1]):
            self.omega = self.omega[::-1, :]
            self.omega[:, 0] = self.span - self.omega[:, 0]

        for i, S in enumerate(self.__strata):
            if S.name != 'DEBRIS':
                break
            if i == 0:
                omega = self.omega
            else:
                omega = self.__strata[i - 1].omega
            S.omega = self.__runoff__(omega, S.h, cr)

    def __parse__(self, name):
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

    def __runoff__(self, omega, x, cr):
        def velocity(arr):
            n = arr.shape[0]
            v = np.zeros(n)
            for i in range(n - 1):
                j = i + 1
                x0 = arr[i]
                x1 = arr[j]
                alpha = np.arctan((x0[1] - x1[1]) / (x1[0] - x0[0]))
                a = G * (np.sin(alpha) - np.cos(alpha) * cr)
                s = np.linalg.norm(x1 - x0)
                p = 2. * v[i] / a
                q = -2. * s / a
                det = (p / 2.)**2 - q
                if det < 0:
                    v[j] = 0.0
                else:
                    t1 = -p / 2. + np.sqrt(det)
                    t2 = -p / 2. - np.sqrt(det)
                    t = np.maximum(t1, t2)
                    v[j] = t * a + v[i]
            return v

        v = velocity(omega) + velocity(omega[::-1])
        v /= np.amax(v)
        d = np.zeros_like(omega)
        d[:, 1] = (x[0] - x[1]) * v + x[1]
        return omega - d

    def bottoms(self, x):
        b = self.b
        for s in self.__strata:
            b -= s.h
            yield x * self.m + self.b

    def iscrack(self, x, eps):
        for crack in self.cracks:
            if np.abs(x - crack) < eps:
                return True
        return False

    def stratum(self, x):
        if type(x) == str:
            for s in self.__strata:
                if s.name == x:
                    return s
        elif type(x) == int:
            return self.__strata[x]
        for y in self.bottoms(x[0]):
            if x[1] > y:
                return s
        return None

    def top(self, x):
        return np.interp(x, self.omega[:, 0], self.omega[:, 1])
