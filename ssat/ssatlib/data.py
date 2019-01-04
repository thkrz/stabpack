import numpy as np

from ssat.ssatlib.constants import G

# def zc(c, phi, gamma):
#     2. * c / gamma * np.tan(np.radians(45. + .5 * phi))


def asarray(obj, rank, dtype):
    if len(obj) == 0:
        return None
    if rank < 2:
        obj = obj[0]
    if rank < 1:
        obj = obj[0]
    if dtype == str:
        return obj
    return np.asarray(obj).astype(dtype)


class Layer:
    def unit_weight(self):
        return self.bulk_density * G


class Slope:
    def __init__(self, name, dim=2):
        self.__dim = dim
        self.__layers = []

        s = [{
            'cracks': ([], 2, np.float),
            'dip': ([], 1, np.float),
            'outcrop': ([], 1, np.float),
            'profile': ([], 2, np.float)
        }, {
            'bulk_density': ([], 0, np.float),
            'cohesion': ([], 0, np.float),
            'friction_angle': ([], 0, np.float),
            'grain_sizes': ([], 2, np.float),
            'hydraulic_conductivity': ([], 0, np.float),
            'name': ([], 0, str),
            'particle_density': ([], 0, np.float),
            'poissons_ratio': ([], 0, np.float),
            'porosity': ([], 0, np.float),
            'retention_model': ([], 1, np.float),
            'thickness': ([], 0, np.float),
            'water_content': ([], 0, np.float),
            'youngs_modulus': ([], 0, np.float)
        }]
        d, k = None, None
        with open(name) as dat:
            for ln in dat:
                ln = ln.strip()
                if not ln:
                    continue
                if ln == 'BEGINSLOPE':
                    d = s[0]
                elif ln == 'ENDSLOPE':
                    break
                elif ln == 'BEGINLAYER':
                    d = s[-1]
                    s.append(d.copy())
                elif ln == 'ENDLAYER':
                    d = s[0]
                elif ln == 'COMMENT':
                    k = ln
                elif ln.lower() in d.keys():
                    k = ln.lower()
                elif k != 'COMMENT':
                    d[k][0].append(ln.split())

        for k, v in s[0].items():
            setattr(self, k, asarray(*v))

        for d in s[1:-1]:
            obj = Layer()
            self.__layers.append(obj)
            for k, v in d.items():
                setattr(obj, k, asarray(*v))

        if self.cracks:
            r = []
            for xy in self.cracks:
                lyr = self[xy]
                z = self.top(xy)
                h = 2. * lyr.cohesion / lyr.unit_weight() * np.tan(
                    np.radians(45. + .5 * lyr.friction_angle))
                r.append(np.array(xy.tolist() + [z - h]))
            self.cracks = np.asarray(r)
        if self.outcrop:
            self.outcrop = np.array(self.outcrop.tolist() +
                                    [self.top(outcrop)])

        for layer in self.__layers:
            if layer.particle_density is None:
                layer.particle_density = 2.65
            if layer.grain_sizes:
                if layer.porosity is None:
                    pass
                if layer.retention_model is None:
                    pass
