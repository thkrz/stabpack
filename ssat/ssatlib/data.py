import numpy as np

from ssat.ssatlib.constants import G

KEYWORDS = {
    'BEGINSLOPE': '',
    'ENDSLOPE': '',
    'BEGINLAYER': '',
    'ENDLAYER': '',
    'COMMENT': '',
    'CRACKS': 'cracks',
    'DIP': 'm',
    'OUTCROP': 'b',
    'PROFILE': 'omega',
    'BULK_DENSITY': 'rho_b',
    'COHESION': 'c',
    'FRICTION_ANGLE': 'phi',
    'GRAIN_SIZES': 'pdf',
    'HYDRAULIC_CONDUCTIVITY': 'K',
    'NAME': 'name',
    'PARTICLE_DENSITY': 'rho_p',
    'POISSONS_RATIO': 'pr',
    'POROSITY': 'n',
    'RETENTION_MODEL': 'beta',
    'THICKNESS': 'h',
    'WATER_CONTENT': 'theta',
    'YOUNGS_MODULUS': 'E'
}


class Layer:
    def depth(self):
        return 2. * self.c / self.gamma() * np.tan(
            np.radians(45. + .5 * self.phi))

    def gamma(self):
        return self.rho_b * G


class Slope:
    def __init__(self, name, dim=2):
        self.__dim = dim
        self.__layers = []

        c = False
        quit = False
        k, o = None, None
        val = []
        with open(name) as dat:
            ln = next(dat).strip()
            while True:
                try:
                    peek = next(dat).strip()
                except StopIteration:
                    if quit:
                        break
                    quit = True
                if ln == "BEGINSLOPE":
                    o = self
                elif ln == "ENDSLOPE":
                    break
                elif ln == "BEGINLAYER":
                    o = Layer()
                    self.__layers.append(o)
                elif ln == "ENDLAYER":
                    o = self
                elif ln == "COMMENT":
                    c = True
                elif ln in KEYWORDS.keys():
                    k = ln
                elif ln:
                    if not c:
                        val.append(ln.split())
                    if peek in KEYWORDS.keys():
                        if c:
                            c = False
                        else:
                            if k == "NAME":
                                val = val[0][0]
                            else:
                                val = np.squeeze(
                                    np.array(val).astype(np.float))
                            setattr(o, KEYWORDS[k], val)
                            val = []
                ln = peek
