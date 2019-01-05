import sys

cfg = {
    'DIMENSION': 2,
    'DATAFILE': '',
    'PRECIPITATION': '',
    'PARADX': 1.0,
    'MINLEN': 2.0,
    'DEGREE': 3
}


def die(s, exit_status=1):
    log(s)
    sys.exit(exit_status)


def log(s, end='\n'):
    sys.stderr.write(s + end)
    sys.stderr.flush()
