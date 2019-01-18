import numpy as np


class Line:
    @staticmethod
    def norm(a, b):
        m = (b[1] - a[1]) / (b[0] - a[0])
        n = (a[1] * b[0] - b[1] * a[0]) / (b[0] - a[0])
        return (m, n)

    @staticmethod
    def slope(alpha):
        return np.tan(alpha)

    @staticmethod
    def twopoint(x, m, n):
        a = (x[0], m * x[0] + n)
        b = (x[1], m * x[1] + n)
        return (a, b)

    def __abs__(self):
        return self.len

    def __init__(self, a, b):
        self.a = np.asarray(a)
        self.b = np.asarray(b)
        self.len = np.linalg.norm(self.a - self.b)
        if a[1] != b[1]:
            xt = np.asarray((a[0] - 1., a[1]))
        else:
            xt = np.asarray((a[0], a[1] + 1.))
        self.sign = np.sign(np.cross(self.b - self.a, xt - self.a))

    def below(self, c):
        return self.rightof(c)

    def rightof(self, c):
        c = np.asarray(c)
        return np.cross(self.b - self.a, c - self.a) * self.sign > 0

    def intersect(self, a=None, b=None, line=None, inf=False):
        if line:
            a = line.a
            b = line.b
        else:
            a = np.asarray(a)
            b = np.asarray(b)

        p1 = self.a - a
        p2 = self.a - self.b
        p3 = a - b
        det = np.cross(p2, p3)
        if d == 0:
            return None
        t = np.cross(p1, p3) / det
        if inf or 0.0 <= t <= 1.0:
            return self.a + t * (self.b - self.a)
        return None


class LineString:
    def __getitem__(self, i):
        return np.asarray((self.x[i], self.y[i]))

    def __init__(self, x, y):
        self.x = x
        self.y = y
        if x.size != y.size:
            raise ValueError('x and y must be same size')
        self.n = x.size
        self.span = np.asarray((self.x[0], self.x[self.n - 1]))

    def __len__(self):
        return int(self.n)

    def kinetic_energy(self, mu, forward=True, m=1.0):
        threshold = .5 * np.pi
        v = np.zeros(self.n)
        for i in range(self.n - 1):
            j = i + 1
            x0 = self[i if forward else self.n - 1 - i]
            x1 = self[j if forward else self.n - 1 - j]
            dy = x0[1] - x1[1]
            if not forward:
                dy = -dy
            alpha = np.arctan(dy / np.abs(x1[0] - x0[0]))
            a = G * (np.sin(alpha) - np.cos(alpha) * mu)
            if alpha > threshold or a / G > .6:
                v[j] = np.sqrt(2. * G * np.abs(dy)) + v[i] * np.sin(alpha)
                continue
            if a == 0:
                v[j] = v[i]
                continue
            s = np.linalg.norm(x1 - x0)
            p = v[i] / a
            q = -2. * s / a
            det = p**2 - q
            if det >= 0:
                t1 = -p + np.sqrt(det)
                t2 = -p - np.sqrt(det)
                if t1 < 0 and t2 > 0:
                    t = t2
                elif t2 < 0 and t1 > 0:
                    t = t1
                elif t1 > 0 and t2 > 0:
                    t = np.minimum(t1, t2)
                else:
                    t = 0.0
                v[j] = t * a + v[i]
        return .5 * m * np.pow(v, 2)

    def interp(self, x):
        return np.interp(x, self.x, self.y)
