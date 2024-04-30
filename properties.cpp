#include "properties.hpp"

double calculate_area_xs(double radius)
{
    return std::numbers::pi * std::pow(radius, 2);
}

double calculate_volume(double radius)
{
    return 4.0 / 3.0 * std::numbers::pi * std::pow(radius, 3);
}

double calculate_mass(double radius, double rho_p, double rho_l)
{
    return (rho_p + 0.5 * rho_l) * calculate_volume(radius);
}

double reynolds_number(Vector u_rel, double radius, double rho_l, double mu_l)
{
    return 2 * rho_l * u_rel.norm() * radius / mu_l;
}

double eotvos_number(Vector g, double rho_p, double rho_l, double radius, double sigma)
{
    return -g[1] * std::abs(rho_l - rho_p) * std::pow(2 * radius, 2) / sigma;
}

double weber_number(double rho_l, Vector u_rel, double radius, double sigma)
{
    return rho_l * u_rel.norm_squared() * radius / sigma;
}

double kinetic_energy(Vector u, double m)
{
    return 0.5 * m * u.norm_squared();
}

double potential_energy(double m, Vector g, Vector x, double ymin)
{
    return m * g[1] * (ymin - x[1]);
}

double mechanical_energy(Vector x, Vector u, double m, Vector g, double ymin)
{
    return kinetic_energy(u, m) + potential_energy(m, g, x, ymin);
}
