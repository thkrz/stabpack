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

    def __getitem__(self, key):
        c = type(key)
        if c == str:
            for s in self.strata:
                if key == s.name:
                    return s
        elif c == int:
            return self.strata[key]
        return None

    def __init__(self):
        self.strata = []

    def dump(self):
        print('BEGINSLOPE')
        print('ENDSLOPE')

    def top(self, x):
        return np.interp(x, self.profile[:, 0], self.profile[:, 1])

    @staticmethod
    def load(name):
        def __toarray(o):
            for k, b in o.__words__.items():
                a = getattr(o, k, None)
                if type(a) == list:
                    setattr(o, k, np.array(a))

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
                    __toarray(s)
                    break
                elif ln == 'BEGINLAYER':
                    c = False
                    o = Layer()
                    s.strata.append(o)
                elif ln == 'ENDLAYER':
                    c = False
                    __toarray(o)
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
        return s
