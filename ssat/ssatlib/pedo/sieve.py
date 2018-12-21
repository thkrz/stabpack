import numpy as np

from collections import OrderedDict

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
    return np.array(list(diameter.values()))


def grade(w):
    def fmt(o):
        q, p = o
        q = q.lower()
        if p > .01:
            if p <= .05:
                q = '(v{})'.format(q)
            elif p <= .2:
                q = '({})'.format(q)
            return q
        return ''

    c = w[0]
    si = w[1:4].sum()
    s = w[4:7].sum()
    m = si + c
    g = w[7:10].sum()
    if g > .01:
        d = [('G', g), ('S', s), ('M', m)]
    else:
        d = [('S', s), ('SI', si), ('C', c)]
    rank = sorted(d, key=lambda x: -x[1])
    return fmt(rank[2]) + fmt(rank[1]) + rank[0][0]


def sort(ppf):
    d95 = phi(ppf(.95))
    d84 = phi(ppf(.84))
    d16 = phi(ppf(.16))
    d5 = phi(ppf(.05))
    si = (d84 - d16) / 4. + (d95 - d5) / 6.6
    if si < .35:
        return 'very well sorted'
    if si < .5:
        return 'well sorted'
    if si < 1.:
        return 'moderately sorted'
    if si < 2.:
        return 'poorly sorted'
    if si < 4.:
        return 'very poorly sorted'
    return 'extremely poorly sorted'


def phi(d):
    return -np.log2(d)
