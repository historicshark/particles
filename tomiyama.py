import numpy as np


def drag_coefficient_tomiyama_pure(re, eo):
    cd = np.maximum(np.minimum(16 / re * (1 + .15 * re**.687), 48 / re), 8 / 3 * eo / (eo + 4))
    cd[cd > 200] = 200
    return cd


def drag_coefficient_tomiyama_low(re, eo):
    cd = np.maximum(np.minimum(24 / re * (1 + .15 * re ** .687), 72 / re), 8 / 3 * eo / (eo + 4))
    cd[cd > 200] = 200
    return cd


def drag_coefficient_tomiyama_high(re, eo):
    cd = np.maximum(24 / re * (1 + .15 * re) ** .687, 8 / 3 * eo / (eo + 4))
    cd[cd > 200] = 200
    return cd


def lift_coefficient_tomiyama(re, eo, radius):
    d_h = 2 * radius * (1 + .1163 * eo**.757)
    eo_d = eo * d_h**2 / (2 * radius)**2
    f = 0.00105 * eo_d**3 - .0159 * eo_d**2 - .0204 * eo_d + .474
    ind = eo_d < 4
    cl = np.zeros_like(re)
    cl[ind] = np.minimum(0.288 * np.tanh(0.121 * re[ind]), f[ind])
    cl[~ind] = f[~ind]
    return cl
