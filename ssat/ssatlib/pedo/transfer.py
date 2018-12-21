n = 0
U = 0
d = None

def meter(x):
    return x * .001

def hazen():
    N = 6e-4
    phi = 1. + 10. * (n - 0.26)
    de = meter(d(.1))
    return N * phi * de

def slichter():
    N = 1e-2
    phi = n**3.287
    de = meter(d(.1))
    return N * phi * de

