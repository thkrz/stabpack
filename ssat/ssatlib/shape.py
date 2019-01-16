import numpy as np


class Boundary:
    def __init__(self, x, y):
        self.x = x
        self.y = y
        if x.size != y.size:
            raise ValueError('x and y must be same size')
        self.n = x.size
        self.span = np.asarray((self.x[0], self.x[self.n - 1]))

    def __len__(self):
        return int(np.round(self.span[1]))

    def __vdist(self, mu, forward):
        threshold = .5 * np.pi
        v = np.zeros(self.n)
        for i in range(self.n - 1):
            j = i + 1
            k = i if forward else self.n - 1 - i
            l = j if forward else k - 1
            x0 = np.asarray((self.x[k], self.y[l]))
            x1 = np.asarray((self.x[k], self.y[l]))
            dy = x0[1] - x1[1]
            if forward:
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
        return v

    def interp(self, x):
        return np.interp(x, self.x, self.y)

    def vdist(self, mu):
        return self.__vdist(mu, True) + self.__vdist(mu, False)


class Edge:
    def __init__(self, m, b):
        self.m = m
        self.b = b

    def above(self, x):
        return x[1] <= self.m * x[0] + self.b

    def below(self, x):
        return x[1] >= self.m * x[0] + self.b

    def leftof(self, x):
        return self.above(x)

    def rightof(self, x):
        return self.below(x)
