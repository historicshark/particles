#include "coalescence_breakup.hpp"

std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<> dis(0.0, 1.0);

bool detect_breakup(double we, double eo)
{
    double e_b = 0.00105 * std::pow(eo,3) - .0159 * std::pow(eo,2) - .0204 * eo + .474;
    double zeta = std::pow((1. + 2. * std::pow(e_b, 1.6075)) / (3. * std::pow(e_b,2.*1.6075/3.)), -1./1.6075);
    return we > 12. * zeta;
}

double radius_fraction()
{
    double x = dis(gen);
    return std::pow(u_distribution_factor * std::pow(x, -0.5) * std::pow(1.0 - x, -0.5), 1./3.);
}

std::tuple<double, double, Vector> breakup_radii_and_new_position(double r, Vector x)
{
    double ratio = radius_fraction();
    double r1 = ratio * r;
    double r2 = (1.0 - ratio) * r;
    double new_position_magnitude = 1.2 * (r1 + r2);
    Vector vec = {dis(gen), dis(gen)};
    Vector new_position = new_position_magnitude * vec / vec.norm();
    return {r1, r2, new_position};
}

double equivalent_diameter(double r1, double r2)
{
    return 2. / (1. / (2. * r1) + 1. / (2. * r2));
}

double drainage_time(double deq, double rho_l, double sigma)
{
    return std::sqrt(std::pow(deq, 3) * rho_l / (128. * sigma)) * log_factor;
}

std::tuple<double, Vector> coalesced_bubble_radius_and_velocity(double r1, double r2, Vector u1, Vector u2)
{
    double d13 = 8. * std::pow(r1, 3);
    double d23 = 8. * std::pow(r2, 3);
    
    double r_new = std::pow(d13 + d23, 1./3.) / 2.;
    Vector u_new = (d13 * u1 + d23 * u2) / (d13 + d23);
    return {r_new, u_new};
}
