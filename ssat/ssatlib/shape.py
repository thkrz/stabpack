import numpy as np

from collections import namedtuple

Bounds = namedtuple('Bounds', ['xmin', 'ymin', 'xmax', 'ymax'])
Point = namedtuple('Point', ['x', 'y'])


@dataclass
class Edge:
    a: Point
    b: Point
    visited: bool


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

        xmin = points[0].x
        ymin = points[0].y
        xmax = points[0].x
        ymax = points[0].y
        for x, y in points:
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
        self.curr = 0
        return self

    def __len__(self):
        return self.__len + 1

    def __next__(self):
        i = self.curr
        if i > self.__len or self.__len == 0:
            raise StopIteration
        if i == self.__len:
            i = 0
        self.curr += 1
        return self.__points[i]

    def bounds(self):
        return self.__bounds

    def contains(self, p):
        def wind(a, b, c) -> int:
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
            if det > 0:
                return -1
            elif det < 0:
                return 1
            return 0

        p = Point(p)
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
            e.append(Edge(self.__points[i], self.__points[j], False))
        return e

    def split(self, a, b):
        pass


class MultiPolygon:
    def __init__(self, polygons):
        self.__polygons = polygons
        self.features = {}

    def __iter__(self):
        self.curr = 0
        return self

    def __next__(self):
        i = self.curr
        if i == len(self.__polygons):
            raise StopIteration
        self.curr += 1
        return self.__polygons[i]

    def contains(self, p):
        for poly in self.__polygons:
            if poly.contains(p):
                return True
        return False
