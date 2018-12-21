import sys

cfg = {
    'DIMENSION': 2,
    'DATAFILE': '',
    'LOGFILE': 'ssat.log',
    'HYDRO_MODEL': 'auto',
    'STABILITY_MODEL': 'auto',
    'SUCTION_MODEL': 'auto'
}


def die(s, exit_status=1):
    log(s)
    sys.exit(exit_status)


def log(s, end='\n'):
    sys.stderr.write(s + end)
    sys.stderr.flush()
