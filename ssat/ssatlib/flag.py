import sys

flags = {}
__argv0 = sys.argv.pop(0)


def parse() -> bool:
    argv = iter(sys.argv)
    p = []
    for a in argv:
        if a.startswith('-'):
            v = None
            if '=' in a:
                k, v = tuple(a.split('=', 2))
            else:
                k = a
            k = k[1:]
            if k not in flags.keys():
                return True
            c = type(flags[k])
            if c == bool:
                flags[k] = True
            else:
                if v is None:
                    try:
                        v = next(argv)
                    except StopIteration:
                        return True
                    p.append(v)
                flags[k] = c(v)
            p.append(a)
    for e in p:
        sys.argv.remove(e)
    return False
