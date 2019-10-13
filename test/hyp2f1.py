#!/usr/bin/python3
import sys

from scipy.special import hyp2f1

params = [float(n) for n in sys.argv[1:]]
print(hyp2f1(*params))
