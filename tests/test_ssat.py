#!/usr/bin/env python3
import os
import shlex
import sys
sys.path.insert(0,
                os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from ssat import cfg, die
from ssat import mesh

rc = {'mesh': mesh.main}


def usage():
    die('usage: ssat cfgfile command...')


if __name__ == '__main__':
    if len(sys.argv) < 3:
        usage()

    with open(sys.argv.pop(1)) as f:
        for tok in shlex.split(f.read(), comments=True):
            key, val = tuple(tok.split('=', 2))
            if key not in cfg.keys():
                die('invalid config file.')
            try:
                cfg[key] = type(cfg[key])(val)
            except ValueError:
                die('invalid config file.')
    while len(sys.argv) > 1:
        cmd = sys.argv.pop(1)
        sys.argv[0] = cmd
        if cmd not in rc.keys():
            usage()
        if rc[cmd]():
            break
    sys.exit(0)
