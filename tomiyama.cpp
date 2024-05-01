#include "tomiyama.hpp"

double drag_coefficient_tomiyama_pure(double re, double eo)
{
    double cd = std::max(std::min(16. / re * (1. + .15 * std::pow(re,0.687)), 48. / re), 8. / 3. * eo / (eo + 4.));
    if (cd > 200)
    {
        cd = 200;
    }
    return cd;
}

double drag_coefficient_tomiyama_low(double re, double eo)
{
    double cd = std::max(std::min(24. / re * (1. + .15 * std::pow(re,0.687)), 72. / re), 8. / 3. * eo / (eo + 4.));
    if (cd > 200)
    {
        cd = 200;
    }
    return cd;
}

double drag_coefficient_tomiyama_high(double re, double eo)
{
    double cd = std::max(24. / re * (1. + .15 * std::pow(re,0.687)), 8. / 3. * eo / (eo + 4.));
    if (cd > 200)
    {
        cd = 200;
    }
    return cd;
}

double lift_coefficient_tomiyama(double re, double eo_d, double radius)
{
    double f = 0.00105 * std::pow(eo_d,3) - .0159 * std::pow(eo_d,2) - .0204 * eo_d + .474;
    double cl;
    if (eo_d > 4)
    {
        cl = f;
    }
    else
    {
        cl = std::max(0.288 * std::tanh(0.121 * re), f);
    }
    return cl;
}
