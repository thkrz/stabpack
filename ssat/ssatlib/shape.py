import numpy as np


@dataclass
class Edge:
    a: (float, float)
    b: (float, float)
    visited: bool


class Polygon:
    def __init__(self, points=[]):
        self.points = points

    def __iter__(self):
        self.curr = 0
        return self

    def __len__(self):
        return len(self.points) + 1

    def __next__(self):
        n = len(self.points)
        i = self.curr
        if i > n or n == 0:
            raise StopIteration
        if i == n:
            i = 0
        self.curr += 1
        return self.points[i]

    def bounds(self):
        pass

    def edges(self):
        n = len(self.points)
        e = []
        for i in range(n):
            j = i + 1
            if j == n:
                j = 0
            e.append(Edge(self.points[i], self.points[j], False))
        return e

    def contains(self, p):
        pass

    def split(self, a, b):
        pass

    def toraster(self, cell_size):
        pass
