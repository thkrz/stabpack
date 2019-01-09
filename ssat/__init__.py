import sys

cfg = {
    'DIMENSION': 2,
    'DATAFILE': '',
    'CRR': 0.2,
    'FOS': 1.3,
    'MAXMU': 2.0,
    'PRECIPITATION': '',
    'PARADX': 1.0,
    'NPARAM': 100,
    'MAXDEPTH': 0.3,
    'MINDESC': 10.0,
    'MINLEN': 2.0,
    'DEGREE': 3
}


def die(s, exit_status=1):
    log(s)
    sys.exit(exit_status)


def log(s, end='\n'):
    sys.stderr.write(s + end)
    sys.stderr.flush()
