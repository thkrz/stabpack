def phi(d):
    return -np.log2(d)

class Distribution:
    def __init__(self, d, arr):
        self.d = {5: d(.05)}

    def issorted(self):
        si = (self.d[84] - self.d[16]) / 4. + (self.d[95] - self.d[5]) / 6.6
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
