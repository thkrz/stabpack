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


class Profile:
    def __init__(self, x, y):
        self.x = x
        self.y = y
        if x.size != y.size:
            raise ValueError
        self.n = x.size
        if self.y[self.n - 1] > np.mean(self.y):
            self.mirror()

    def __velocity(self, cr):
        v = np.zeros(self.n)
        for i in range(self.n - 1):
            j = i + 1
            x0 = np.asarray((self.x[i], self.y[i]))
            x1 = np.asarray((self.x[j], self.y[j]))
            alpha = np.arctan((x0[1] - x1[1]) / (x1[0] - x0[0]))
            a = G * (np.sin(alpha) - np.cos(alpha) * cr)
            s = np.linalg.norm(x1 - x0)
            p = 2. * v[i] / a
            q = -2. * s / a
            det = (p / 2.)**2 - q
            if det >= 0:
                t1 = -p / 2. + np.sqrt(det)
                t2 = -p / 2. - np.sqrt(det)
                if t1 < 0:
                    t = t2
                elif t2 < 0:
                    t = t1
                elif t1 > 0 and t2 > 0:
                    t = np.minimum(t1, t2)
                else:
                    t = 0.0
                v[j] = t * a + v[i]
        return v

    def mirror(self):
        xmax = self.x[self.n - 1]
        self.x = xmax - self.x[::-1]
        self.y = self.y[::-1]

    def span(self):
        return self.x[self.n - 1]

    def top(self, x):
        return np.interp(x, self.x, self.y)

    def velocity(self, cr):
        v = self.__velocity(cr)
        self.mirror()
        v += self.__velocity(cr)
        self.mirror()
        return v


class Stratum:
    def depth(self):
        return 2. * self.c / self.gamma() * np.tan(
            np.radians(45. + .5 * self.phi))

    def gamma(self):
        return self.rho_b * G

    def isdebris(self):
        return self.name == 'DEBRIS'


class Slope:
    def __init__(self, name, dim=2, cr=.2):
        self.dim = dim

        self.__strata = []
        self.__parse(name)

        self.omega = Profile(self.omega[:, 0], self.omega[:, 1])
        self.span = self.omega.span()

        self.b = self.top(self.a)
        self.m = np.arctan(self.alpha)

        top = self.omega
        for s in self.__strata:
            if not s.isdebris():
                break
            v = top.velocity(cr)
            d = (s.h[1] - s.h[0]) * (v / np.amax(v)) + s.h[0]
            s.omega = Profile(self.omega.x, top.y - d)
            top = s.omega

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

    def bottoms(self, x):
        b = self.b
        for s in self.__strata:
            if s.isdebris():
                y = s.omega.top(x)
            else:
                b -= s.h
                y = x * self.m + self.b
            yield y

    def iscrack(self, x, eps):
        for crack in np.asarray(self.cracks):
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

    def top(self, x, relief=True):
        if relief:
            return self.omega.top(x)
        stratum = None
        for s in self.__strata:
            if s.isdebris():
                stratum = s
            else:
                break
        if stratum:
            return stratum.omega.top(x)
        return None
