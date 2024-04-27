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
