import numba
import numpy as np

from collections import namedtuple

Bounds = namedtuple('Bounds', ['xmin', 'ymin', 'xmax', 'ymax'])
Point = namedtuple('Point', ['x', 'y'])


class Edge:
    def __init__(self, a, b):
        self.a = a
        self.b = b

    def includes(self, p):
        d = (p.x - self.a.x) * (self.b.y - self.a.y) - (p.y - self.a.y) * (
            self.b.x - self.a.x)
        return d <= 0

    def intersect(self, a, b):
        def det(a, b):
            return a.x * b.y - a.y * b.x

        dx = Point((a.x - b.x, self.a.x - self.b.x))
        dy = Point((a.y - b.y, self.a.y - self.b.y))
        div = det(dx, dy)
        if div == 0:
            return None
        d = Point((det(a, b), det(self.a, self.b)))
        x = det(d, dx) / div
        y = det(d, dy) / div
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

        xmin = self.__points[0].x
        ymin = self.__points[0].y
        xmax = self.__points[0].x
        ymax = self.__points[0].y
        for x, y in self.__points:
            if xmin > x:
                xmin = x
            if xmax < x:
                xmax = x
            if ymin > y:
                ymin = y
            if ymax < y:
                ymax = y
        self.__bounds = Bounds(xmin, ymin, xmax, ymax)

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
            a += e.a.x * e.b.y - e.b.x * e.a.y
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
                if e.includes(p):
                    if not e.includes(s):
                        points.append(e.intersect(s, p))
                    points.append(p)
                elif e.includes(s):
                    points.append(e.intersect(s, p))
                s = p
        if len(points) < 3:
            return None
        return Polygon(points)

    @numba.jit
    def contains(self, p):
        def wind(a, b, c):
            if a.y == b.y == c.y:
                if b.x <= a.x <= c.x or c.x <= a.x <= b.x:
                    return 0
                return 1
            if a.y == b.y and a.x == b.x:
                return 0
            if b.y > c.y:
                swap = b.y
                b.y = c.y
                c.y = swap
            if a.y <= b.y or a.y > c.y:
                return 1
            det = (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x)
            return np.sign(det)

        xmin, ymin, xmax, ymax = self.__bounds
        if xmin > p.x or xmax < p.x or ymin < p.y or ymax > p.y:
            return False
        t = -1
        for i in range(self.__len):
            j = i + 1
            if j == self.__len:
                j = 0
            t *= wind(p, self.__points[i], self.__points[j])
            if t == 0:
                break
        return t == 1 or t == 0

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
        points = []
        for i, e in enumerate(self.edges()):
            p = e.intersect(a, b)
            if p:
                points.append(p)
            else:
                points.append(e.a)


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
