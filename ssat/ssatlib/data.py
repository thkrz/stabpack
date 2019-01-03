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
        s = {(k, []) for k in ['cracks', 'dip', 'outcrop', 'profile']}
        l = {(k, [])
             for k in [
                 'bulk_density', 'cohesion', 'friction_angle',
                 'grain_size_distribution', 'hydraulic_conductivity', 'name',
                 'particle_density', 'poissons_ratio', 'porosity',
                 'retention_model', 'thickness', 'water_content',
                 'youngs_modulus'
             ]}
        lyr = []
        c = False
        d, k = None, None
        with open(name) as dat:
            for ln in dat:
                ln = ln.strip()
                if ln == 'BEGINSLOPE':
                    c = False
                    d = s
                elif ln == 'ENDSLOPE':
                    break
                elif ln == 'BEGINLAYER':
                    c = False
                    d = l.copy()
                    lyr.append(d)
                elif ln == 'ENDLAYER':
                    c = False
                    d = s
                elif ln == 'COMMENT':
                    c = True
                elif ln.lower() in d.keys():
                    c = False
                    k = ln.lower()
                elif ln == 'ASSUME':
                    print(k)
                elif not c:
                    v = ln.split()
                    for i, e in enumerate(v):
                        if e.isdigit():
                            v[i] = float(e)
                    d[k].append(v)
