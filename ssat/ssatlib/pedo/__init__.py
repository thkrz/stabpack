import numpy as np

from collections import OrderedDict

ppf = None
diameter = OrderedDict({
    'Cl': 0.002,
    'FSi': 0.0063,
    'MSi': 0.02,
    'CSi': 0.063,
    'FSa': 0.2,
    'MSa': 0.63,
    'CSa': 2.0,
    'FGr': 6.3,
    'MGr': 20.0,
    'CGr': 63.0,
    'Co': 200.0,
    'Bo': 630.0,
    'LBo': 2000.0
})


def diaspace():
    return np.array(diameter.values())


def phi(d):
    return -np.log2(d)


def classify(w):
    def less(o):
        q, p = o
        q = q.lower()
        if p > 1.:
            if p <= 5.:
                q = '(v{})'.format(c)
            elif p <= 20.:
                q = '({})'.format(c)
            return q
        return ''

    clay = w[0]
    silt = w[1:4].sum()
    sand = w[4:7].sum()
    mud = silt + clay
    gravel = w[7:10].sum()
    if gravel > 1.:
        classes = [('G', gravel), ('S', sand), ('M', mud)]
    else:
        classes = [('S', sand), ('SI', silt), ('C', clay)]
    rank = sorted(classes, key=lambda x: x[1])
    return less(rank[2]) + less(rank[1]) + rank[0][0]
