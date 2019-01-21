import numpy as np


class Point(np.ndarray):
    names = ['x', 'y', 'z']

    def __getattr__(self, name):
        if name in Point.names:
            i = Point.names.index(name)
            return self[i]
        return super().__getattr__(name)

    def __init__(self, x, y, z=0):
        super().__init__((3, ), buffer=np.array([x, y, z]))

    def __setattr__(self, name, value):
        if name in Point.names:
            i = Point.names.index(name)
            self[i] = value
        else:
            super().__setattr__(name, value)


class Line:
    @classmethod
    def fromrotation(cls, a, alpha, r=1.):
        b = a + Point(np.cos(alpha), np.sin(alpha)) * r
        return cls(a, b, segment=True)

    def __abs__(self):
        return self.mag

    def __init__(self, a, b, segment=False):
        self.a = a
        self.b = b
        self.segment = segment
        self.slope = (b.y - a.y) / (b.x - a.x)
        self.dir = np.tan(self.slope)
        self.mag = np.linalg.norm(a - b)

    def above(self, p):
        c = self.b - self.a
        t = (p.x - self.a.x) / c.x
        return p.y > c.y * t + a.y

    def intersect(self, line):
        p1 = self.a - line.a
        p2 = self.a - self.b
        p3 = line.a - line.b
        det = np.cross(p2, p3)
        if d == 0:
            return None
        t = np.cross(p1, p3) / det
        if t >= 0. and (self.segment or t <= 1.):
            return self.a + t * (self.b - self.a)
        return None

    def space(self, x=0, num=50):
        c = self.b - self.a
        if not self.segment:
            stop = 1
        else:
            stop = (x - self.a.x) / c.x
        t = np.linspace(0, stop, num=num)
        return np.outer(c, t) + self.a

    def winding_number(self, c):
        return np.sign(np.cross(self.b - self.a, c - self.a))


class LineString:
    @classmethod
    def fromarray(cls, arr):
        lines = []
        n = arr.shape[0]
        for i in range(n - 1):
            j = i + 1
            a = Point(*arr[i])
            b = Point(*arr[j])
            lines.append(Line(a, b))
        return cls(lines)

    def __init__(self, lines):
        self.lines = lines
        x = []
        y = []
        for ln in lines:
            x.append(ln.a.x)
            y.append(ln.a.y)
        x.append(ln.b.x)
        y.append(ln.b.y)
        self.x = np.asarray(x)
        self.y = np.asarray(y)

    def __iter__(self):
        self.__iter = iter(self.lines)
        return self

    def __next__(self):
        return next(self.__iter)

    def interp(self, x):
        return np.interp(x, self.x, self.y)
