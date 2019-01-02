import numpy as np


class Layer:
    __words__ = {
        'bulk_density': False,
        'cohesion': False,
        'friction_angle': False,
        'grain_size_distribution': True,
        'hydraulic_conductivity': False,
        'name': False,
        'particle_density': False,
        'poissons_ratio': False,
        'porosity': False,
        'retention_model': False,
        'thickness': False,
        'water_content': False,
        'youngs_modulus': False
    }


class Slope:
    __words__ = {
        'cracks': True,
        'dip': False,
        'outcrop': False,
        'profile': True
    }

    def __init__(self, dim=2):
        self.__dim = dim
        self.__layers = []

    def dump(self):
        print('BEGINSLOPE')
        print('ENDSLOPE')

    def top(self, x):
        if self.__dim == 2:
            return np.interp(x, self.profile[:, 0], self.profile[:, 1])
        return None

    @staticmethod
    def load(name):
        def toarray(o):
            for k, _ in o.__words__.items():
                a = getattr(o, k, None)
                if type(a) == list:
                    setattr(o, k, np.array(a))
                elif a is None:
                    setattr(o, k, None)

        s = Slope()
        c = False
        o, k = None, None
        with open(name) as dat:
            for ln in dat:
                ln = ln.strip()
                if not ln:
                    continue
                if ln == 'BEGINSLOPE':
                    c = False
                    o = s
                elif ln == 'ENDSLOPE':
                    break
                elif ln == 'BEGINLAYER':
                    c = False
                    o = Layer()
                    s.__layers.append(o)
                elif ln == 'ENDLAYER':
                    c = False
                    o = s
                elif ln == 'COMMENT':
                    c = True
                elif ln.lower() in o.__words__.keys():
                    c = False
                    k = ln.lower()
                elif not c:
                    v = ln.split()
                    for i, e in enumerate(v):
                        if e.isdigit():
                            v[i] = float(e)
                    if len(v) == 1:
                        v = v[0]
                    a = getattr(o, k, None)
                    if o.__words__[k]:
                        if a is None:
                            a = []
                            setattr(o, k, a)
                        a.append(v)
                    else:
                        setattr(o, k, v)

        toarray(s)
        for o in s.__layers:
            toarray(o)
        s.profile[:, 0] -= s.profile[0, 0]
        n = s.profile.shape[0] - 1
        s.span = s.profile[n, 0]
        return s
