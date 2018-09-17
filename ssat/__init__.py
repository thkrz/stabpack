import sys

cfg = {'DIMENSION': 2, 'DATAFILE': ''}


def die(s):
    log(s)
    sys.exit(1)


def log(s, end='\n'):
    sys.stderr.write(s + end)
    sys.stderr.flush()
