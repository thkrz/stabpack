import numpy as np
import re

from ssat.ssatlib.shape import Polygon


@dataclass
class Layer:
    c: float
    phi: float
    weights: np.array


class Slope:
    def __init__(self, name):
        comments = re.compile('#.*$')
        lines = []
        with open(name) as f:
            for s in f:
                s = comments.sub('', s).strip()
                if s:
                    lines.append(s)
        lines = iter(lines)
        n = int(next(lines))
        points = []
        for i in range(n):
            s = next(lines).split()
            points.append([float(d) for d in s])
        self.relief = np.array(points)
        points.append((points[0][0], points[-1][1]))
        p = Polygon(points)
