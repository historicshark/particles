#include "integration.hpp"


Walls wall_overlap(Vector x, double radius, Walls walls)
{
    Walls overlap = {radius - (x[0] - walls[0]),
                     radius - (walls[1] - x[0]),
                     radius - (x[1] - walls[2]),
                     radius - (walls[3] - x[1])};
    
    return overlap;
}

WallContact detect_wall_contact(Walls overlap)
{
    WallContact result;
    std::transform(overlap.cbegin(), overlap.cend(), result.begin(), [](const double& x){return x > 0;});
    return result;
}

bool detect_particle_contact(Vector position1, Vector position2, double radius1, double radius2)
{
    double distance = (position1 - position2).norm();
    return distance < radius1 + radius2;
}

std::vector<std::tuple<size_t, size_t>> new_particle_contact_list(const std::vector<Vector>& position, const std::vector<double>& radius)
{
    std::vector<std::tuple<size_t, size_t>> particle_contact;
    for (size_t i = 0; i != position.size(); i++)
    {
        for (size_t j = i + 1; j != position.size(); j++)
        {
            if (detect_particle_contact(position[i], position[j], radius[i], radius[j]))
            {
                particle_contact.push_back({i,j});
            }
        }
    }
    return particle_contact;
}

Vector acceleration_drag(Vector particle_velocity,
                         Vector flow_velocity,
                         double radius,
                         double rho_l,
                         double mu_l,
                         double mass,
                         Vector g,
                         double rho_p,
                         double sigma)
{
    Vector u_rel = flow_velocity - particle_velocity;
    double re = reynolds_number(u_rel, radius, rho_l, mu_l);
    double eo = eotvos_number(g, rho_p, rho_l, radius, sigma);
    double cd = drag_coefficient_tomiyama_pure(re, eo);
    return 0.5 * rho_l * u_rel.norm() * u_rel * calculate_area_xs(radius) * cd / mass;
}


Vector acceleration_gravity(double rho_p, double rho_l, Vector g)
{
    return (rho_p - rho_l) / rho_p * g;
}

Vector acceleration_lift(double radius, double rho_l, Vector position, Vector velocity, Vector flow_velocity, std::string flow_type, std::vector<double> parameters)
{
    Vector u_rel = flow_velocity - velocity;
    double curl = curl_of_velocity(position, flow_type, parameters);
    Vector vec = {-u_rel[1] * curl, u_rel[0] * curl};
    return calculate_volume(radius) * rho_l * vec / 4.;
}


std::vector<Vector> apply_accelerations(std::vector<Vector>& position,
                                        std::vector<Vector>& velocity,
                                        std::vector<double>& radius,
                                        std::vector<double>& mass,
                                        double rho_p,
                                        double rho_l,
                                        double mu_l,
                                        double mu_p,
                                        Vector g,
                                        double sigma,
                                        Walls walls,
                                        double dt,
                                        std::string flow_type,
                                        std::vector<double> parameters,
                                        bool drag,
                                        bool gravity,
                                        bool particle_collisions,
                                        bool wall_collisions,
                                        bool lift)
{
    std::vector<Vector> acceleration(position.size());
    
    std::vector<Vector> flow_velocity(position.size());
    std::transform(position.begin(), position.end(), flow_velocity.begin(), [flow_type, parameters](Vector& v){
        return interpolate_flow_properties_fast(v, flow_type, parameters);
    });
    
    for (auto [a, x, u, r, m, u_l] : ranges::views::zip(acceleration, position, velocity, radius, mass, flow_velocity))
    {
        if (drag)
        {
            a += acceleration_drag(u, u_l, r, rho_l, mu_l, m, g, rho_p, sigma);
        }
        if (gravity)
        {
            a += acceleration_gravity(rho_p, rho_l, g);
        }
        if (wall_collisions)
        {
            Walls overlap = wall_overlap(x, r, walls);
            WallContact wall_collision = detect_wall_contact(overlap);
            a += acceleration_wall_collision_simple(wall_collision, overlap, x, u, r, m, mu_p, mu_l, sigma);
        }
        if (lift)
        {
            a += acceleration_lift(r, rho_l, x, u, u_l, flow_type, parameters);
        }
    }
    if (particle_collisions)
    {
        auto contact_list = new_particle_contact_list(position, radius);
        for (auto [i,j] : contact_list)
        {
            //        std::cout << i << " " << j << std::endl;
            acceleration[i] += acceleration_particle_collision_simple(position[i], position[j], velocity[i], radius[i], radius[j], mass[i], mu_l, mu_p, sigma);
            acceleration[j] += acceleration_particle_collision_simple(position[j], position[i], velocity[j], radius[j], radius[i], mass[j], mu_l, mu_p, sigma);
        }
    }
    return acceleration;
}

std::vector<Vector> apply_velocities(std::vector<Vector>& velocity, double dt)
{
    std::vector<Vector> kx;
    std::transform(velocity.cbegin(), velocity.cend(), std::back_inserter(kx), [&dt](const Vector& u){return u * dt;});
    return kx;
}


void integrate_rk4(std::vector<Vector>& x_n,
                   std::vector<Vector>& u_n,
                   std::vector<Vector>& x,
                   std::vector<Vector>& u,
                   double dt,
                   std::vector<double>& radius,
                   std::vector<double>& mass,
                   double rho_p,
                   double rho_l,
                   double mu_l,
                   double mu_p,
                   Vector g,
                   double sigma,
                   Walls walls,
                   std::string flow_type,
                   std::vector<double> parameters,
                   bool drag,
                   bool gravity,
                   bool particle_collisions,
                   bool wall_collisions,
                   bool lift)
{
    auto k1u = apply_accelerations(x_n, u_n,
                                   radius, mass, rho_p, rho_l, mu_l, mu_p, g, sigma, walls, dt, flow_type, parameters, drag, gravity, particle_collisions, wall_collisions, lift);
    
    std::vector<Vector> x_n2;
    std::vector<Vector> u_n2;
    std::transform(x_n.begin(), x_n.end(), u_n.begin(), std::back_inserter(x_n2), [dt](Vector& v, Vector& k){return v + dt * k / 2;});
    std::transform(u_n.begin(), u_n.end(), k1u.begin(), std::back_inserter(u_n2), [dt](Vector& v, Vector& k){return v + dt * k / 2;});
    auto k2u = apply_accelerations(x_n2, u_n2,
                                   radius, mass, rho_p, rho_l, mu_l, mu_p, g, sigma, walls, dt, flow_type, parameters,  drag, gravity, particle_collisions, wall_collisions, lift);
    
    std::vector<Vector> x_n3;
    std::vector<Vector> u_n3;
    std::transform(x_n.begin(), x_n.end(), u_n2.begin(), std::back_inserter(x_n3), [dt](Vector& v, Vector& k){return v + dt * k / 2;});
    std::transform(u_n.begin(), u_n.end(), k2u.begin(), std::back_inserter(u_n3), [dt](Vector& v, Vector& k){return v + dt * k / 2;});
    auto k3u = apply_accelerations(x_n3, u_n3,
                                   radius, mass, rho_p, rho_l, mu_l, mu_p, g, sigma, walls, dt, flow_type, parameters,  drag, gravity, particle_collisions, wall_collisions, lift);
    
    std::vector<Vector> x_n4;
    std::vector<Vector> u_n4;
    std::transform(x_n.begin(), x_n.end(), u_n3.begin(), std::back_inserter(x_n4), [dt](Vector& v, Vector& k){return v + dt * k;});
    std::transform(u_n.begin(), u_n.end(), k3u.begin(), std::back_inserter(u_n4), [dt](Vector& v, Vector& k){return v + dt * k;});
    auto k4u = apply_accelerations(x_n4, u_n4,
                                   radius, mass, rho_p, rho_l, mu_l, mu_p, g, sigma, walls, dt, flow_type, parameters,  drag, gravity, particle_collisions, wall_collisions, lift);
    
    for (size_t i = 0; i != u.size(); i++)
    {
        u[i] = u_n[i] + dt * (k1u[i] + 2 * (k2u[i] + k3u[i]) + k4u[i]) / 6;
    }
    for (size_t i = 0; i != x.size(); i++)
    {
        x[i] = x_n[i] + dt * (u_n[i] + 2 * (u_n2[i] + u_n3[i]) + u_n4[i]) / 6;
    }
}

void integrate_euler(std::vector<Vector>& x_n,
                     std::vector<Vector>& u_n,
                     std::vector<Vector>& x,
                     std::vector<Vector>& u,
                     double dt,
                     std::vector<double>& radius,
                     std::vector<double>& mass,
                     double rho_p,
                     double rho_l,
                     double mu_l,
                     double mu_p,
                     Vector g,
                     double sigma,
                     Walls walls,
                     std::string flow_type,
                     std::vector<double> parameters,
                     bool drag,
                     bool gravity,
                     bool particle_collisions,
                     bool wall_collisions,
                     bool lift)
{
    auto k1u = apply_accelerations(x_n, u_n,
                                   radius, mass, rho_p, rho_l, mu_l, mu_p, g, sigma, walls, dt, flow_type, parameters,  drag, gravity, particle_collisions, wall_collisions, lift);
    
    for (auto [x_n_, u_n_, x_, u_, k1u_] : ranges::views::zip(x_n, u_n, x, u, k1u))
    {
        u_ = u_n_ + dt * k1u_;
        x_ = 0.5 * k1u_ * dt * dt + u_ * dt + x_n_;
    }
}
