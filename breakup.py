import numpy as np
from scipy.special import gamma

rng = np.random.default_rng()
u_distribution_factor = gamma(1) / gamma(0.5)**2


def detect_breakup(we, eo):
    e_b = 0.00105 * eo**3 - .0159 * eo**2 - .0204 * eo + .474
    zeta = ((1 + 2 * e_b**1.6075) / (3 * e_b**(2*1.6075/3)))**(-1/1.6075)
    return we > 12 * zeta


def radius_fraction():
    x = rng.random()
    return (u_distribution_factor * x**(-.5) * (1 - x)**(-.5))**(1/3)


def new_bubble(r, x):
    # not vectorized
    ratio = radius_fraction()
    r1 = ratio * r
    r2 = (1 - ratio) * r

    new_position_magnitude = 1.2 * (r1 + r2)
    vec = rng.random((2,))
    new_position = new_position_magnitude * vec / np.sqrt(vec @ vec)
    return r1, r2, new_position
