import numpy as np


def drag_coefficient_mei(re, eo):
    cd = 24.0 / re * (2.0 / 3.0 + (12.0 / re + 3.0 / 4.0 * (1 + 3.315 / np.sqrt(re)))**(-1))
    cd[cd > 200] = 200
    return cd
