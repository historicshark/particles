#pragma once

#include <cmath>

double drag_coefficient_mei(double re, double eo)
{
    double cd = 24.0 / re * (2.0 / 3.0 + std::pow(12.0 / re + 3.0 / 4.0 * (1 + 3.315 / np.sqrt(re)), -1));
    if (cd > 200)
    {
        cd = 200;
    }
    return cd;
}
