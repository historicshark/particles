#include "simple.hpp"

double fluid_film_thickness(double mu_l, double mu_p, double radius, double sigma, bool is_particle_collision)
{
    double factor = std::sqrt(3. / 8.);
    if (is_particle_collision)
    {
        factor = std::sqrt(3. / 16.);
    }
    
    return factor * std::sqrt(mu_l * mu_p * std::pow(radius, 2) / sigma);
}

double contact_radius(double radius, double deformation)
{
    return radius * (-30. + std::sqrt(900. + .424 * deformation / radius));
}

Vector force_elastic_simple(double radius, double deformation, double sigma, Vector n)
{
    return -n * (radius * sigma * (18.5 * std::pow(deformation / radius, 2) + 2. * deformation / radius));
}

Vector force_viscous_simple(double deformation, double h0, double radius, Vector velocity, double mu_l, Vector n, bool is_particle_collision)
{
    double r_a = contact_radius(radius, deformation);
    
    double cbc = 1;
    if (is_particle_collision)
    {
        cbc = .25;
    }
    
    return -n * (velocity.dot(n) * cbc * 6 * mu_l / std::numbers::pi 
                 * .34 * std::pow(deformation / radius + .0002, -0.5) 
                 * (4 * std::sqrt(std::pow(radius, 3) / h0) + 3. * r_a * radius / h0));
}

Vector acceleration_particle_collision_simple(Vector position, Vector other_position, Vector velocity, double radius, double other_radius, double mass, double mu_l, double mu_p, double sigma)
{
    double h0 = fluid_film_thickness(mu_l, mu_p, radius, sigma, true);
    double deformation = radius - 0.5 * (other_position - position).norm() + h0 / 2.;
    if (deformation < 0)
    {
        deformation = 0;
    }
    Vector distance_vector = other_position - position;
    Vector n = distance_vector / distance_vector.norm();
    Vector f_e = force_elastic_simple(radius, deformation, sigma);
    Vector f_v = force_viscous_simple(deformation, h0, radius, velocity, mu_l, mu_p, n, true);
    Vector a = (f_e + f_v) / mass;
    return a;
}

Vector acceleration_wall_collision_simple(std::array<bool, 4> wall_collision, std::array<double, 4> wall_overlap, Vector position, Vector velocity, double radius, double mass, double mu_p, double mu_l, double sigma)
{
    Vector a;
    double h0 = fluid_film_thickness(mu_l, mu_p, radius, sigma, false);
    for (int i = 0; i != 4; i++)
    {
        if (wall_collision[i])
        {
            double deformation = wall_overlap[i] + h0;
            Vector f_e = force_elastic_simple(radius, deformation, sigma);
            Vector f_v;
            if (wall_overlap[i] > h0)
            {
                f_v = force_viscous_simple(deformation, h0, radius, velocity, mu_l, wall_normals[i], false);
            }
            a += (f_e + f_v) / m;
        }
    }
    return a;
}
