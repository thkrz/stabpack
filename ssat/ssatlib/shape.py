import numba
import numpy as np

from collections import namedtuple


def __cross(a, b):
    return a.x * b.y - a.y * b.x


class Bounds:
    def __init__(self, points):
        self.xmin = points[0].x
        self.ymin = points[0].y
        self.xmax = points[0].x
        self.ymax = points[0].y
        for x, y in points:
            if self.xmin > x:
                self.xmin = x
            if self.xmax < x:
                self.xmax = x
            if self.ymin > y:
                self.ymin = y
            if self.ymax < y:
                self.ymax = y

    def contains(self, p):
        return self.xmin <= p.x <= self.xmax and self.ymin <= p.y <= self.ymax


class Point:
    def __add__(self, p):
        return Point((self.x + p.x, self.y + p.y))

    def __eq__(self, p):
        return self.x == p.x and self.y == p.y

    def __init__(self, x):
        self.x = x[0]
        self.y = x[1]

    def __sub__(self, p):
        return Point((self.x - p.x, self.y - p.y))


class Edge:
    def __init__(self, a, b):
        self.a = a
        self.b = b

    def inside(self, p):
        return __cross(p - self.a, self.b - self.a) <= 0

    def intersect(self, a, b):
        dx = Point((a.x - b.x, self.a.x - self.b.x))
        dy = Point((a.y - b.y, self.a.y - self.b.y))
        div = __cross(dx, dy)
        if div == 0:
            return None
        d = Point((__cross(a, b), __cross(self.a, self.b)))
        x = __cross(d, dx) / div
        y = __cross(d, dy) / div
        return Point((x, y))


class Polygon:
    def __getitem__(self, i):
        if i == self.__len:
            i = 0
        return self.__points[i]

    def __init__(self, points, closed=False):
        if closed:
            points = points[:-1]
        self.__points = [Point(p) for p in points]
        self.__len = len(self.__points)
        self.__bounds = Bounds(self.__points)

    def __iter__(self):
        self.__curr = 0
        return self

    def __len__(self):
        return self.__len + 1

    def __next__(self):
        i = self.__curr
        if i > self.__len or self.__len == 0:
            raise StopIteration
        if i == self.__len:
            i = 0
        self.__curr += 1
        return self.__points[i]

    def area(self):
        a = 0
        for e in self.edges():
            a += __cross(e.a, e.b)
        return .5 * a

    def bounds(self):
        return self.__bounds

    def clip(self, p):
        points = self.__points.copy()
        for e in p.edges():
            points_c = points.copy()
            points.clear()
            s = points_c[-1]
            for p in points_c:
                if e.inside(p):
                    if not e.inside(s):
                        points.append(e.intersect(s, p))
                    points.append(p)
                elif e.inside(s):
                    points.append(e.intersect(s, p))
                s = p
        if len(points) < 3:
            return None
        return Polygon(points)

    @numba.jit
    def contains(self, p):
        def rot(a, b, c):
            if a.y == b.y == c.y:
                if b.x <= a.x <= c.x or c.x <= a.x <= b.x:
                    return 0
                return 1
            if a == b:
                return 0
            if b.y > c.y:
                swap = b.y
                b.y = c.y
                c.y = swap
            if a.y <= b.y or a.y > c.y:
                return 1
            return np.sign(__cross(b - a, c - a))

        if not self.__bounds.contains(p):
            return False
        t = -1
        for i in range(self.__len):
            j = i + 1
            if j == self.__len:
                j = 0
            t *= rot(p, self.__points[i], self.__points[j])
            if t == 0:
                break
        return t == -1 or t == 0

    def edges(self):
        e = []
        for i in range(self.__len):
            j = i + 1
            if j == self.__len:
                j = 0
            e.append(Edge(self.__points[i], self.__points[j]))
        return e

    @numba.jit
    def split(self, a, alpha):
        x = self.__bounds.xmin if alpha > 0 else self.__bounds.xmax
        b = Point((x, x * np.tan(alpha) + a.y))
        sec = []
        for e in self.edges():
            p = e.intersect(a, b)
            if p:
                sec.append(p)
        n = len(sec) / 2 + 1
        P = [] * n
        return (Polygon(A), Polygon(B))


class MultiPolygon:
    def __init__(self, polygons):
        self.__polygons = polygons
        self.features = {}

    def __iter__(self):
        self.__curr = 0
        return self

    def __next__(self):
        i = self.__curr
        if i == len(self.__polygons):
            raise StopIteration
        self.__curr += 1
        return self.__polygons[i]

    def contains(self, p):
        for poly in self.__polygons:
            if poly.contains(p):
                return True
        return False
