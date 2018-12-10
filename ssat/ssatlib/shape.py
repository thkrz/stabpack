import numpy as np


@dataclass
class Edge:
    a: (float, float)
    b: (float, float)
    visited: bool


class Polygon:
    def __init__(self, points):
        self.__points = points
        self.__len = len(points)

    def __iter__(self):
        self.curr = 0
        return self

    def __len__(self):
        return self.__len

    def __next__(self):
        i = self.curr
        if i > self.__len or self.__len == 0:
            raise StopIteration
        if i == self.__len:
            i = 0
        self.curr += 1
        return self.__points[i]

    def bounds(self):
        xmin = self.__points[0][0]
        ymin = self.__points[0][1]
        xmax = self.__points[0][0]
        ymax = self.__points[0][1]
        for x, y in self.__points:
            if xmin > x:
                xmin = x
            if xmax < x:
                xmax = x
            if ymin > y:
                ymin = y
            if ymax < y:
                ymax = y
        return (xmin, ymin, xmax, ymax)

    def edges(self):
        e = []
        for i in range(self.__len):
            j = i + 1
            if j == self.__len:
                j = 0
            e.append(Edge(self.__points[i], self.__points[j], False))
        return e

    def contains(self, p):
        pass

    def split(self, a, b):
        pass

    def toraster(self, cell_size):
        pass
