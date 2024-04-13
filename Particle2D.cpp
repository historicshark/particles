auto Particle::distance_vector(Particle& other) { return other.position() - position(); }
auto Particle::distance_vector(Vector other) { return other - position(); }

auto Particle::distance(Particle& other) { return distance_vector(other).norm(); }
auto Particle::distance(Vector other) { return distance_vector(other).norm(); }

void Particle::set_fluid_properties(double rho, double mu)
{
    rho_l = rho;
    mu_l = mu;
}

void Particle::update_reynolds_number(Vector u_p)
{
    re = rho_l * (u_p - v_l).norm() * diameter() / mu_l;
}

auto Particle::kinetic_energy()
{
    return 0.5 * mass() * velocity().norm_squared();
}

auto Particle::drag_coefficient()
{
    if (re > std::numeric_limits<double>::epsilon())
    {
        return 24.0 / re * (2.0 / 3.0 + std::pow(12.0 / re + 3.0 / 4.0 * (1 + 3.315 / std::sqrt(re)), -1));
    }
    else
    {
        return 0.0;
    }
};

auto Particle::drag_acceleration(Vector u_p)
{
    update_reynolds_number(u_p);
    return 1.0 / 2.0 * rho_l / rho * (v_l - u_p).norm() * (v_l - u_p) * area_xs() / volume() * drag_coefficient();
};

auto Particle::buoyancy_acceleration(Vector g = {0, -9.81})
{
    return (rho - rho_l) / rho * g;
}

auto Particle::wall_contact_acceleration(const double dt, double epsilon, std::vector<double> walls)
{
    double fx = 0, fy = 0;
    
    if(r > distance({walls[0],x[1]}) || r > distance({walls[1],x[1]}))
    {
        fx += -(1 + epsilon) * u[0] / dt;
    }
    if(r > distance({x[0],walls[2]}) || r > distance({x[0],walls[3]}))
    {
        fy += -(1 + epsilon) * u[1] / dt;
    }
    
    return Vector{fx, fy};
}

auto Particle::time_to_collision(Particle& other)
{
    auto a = (velocity_prev() - other.velocity_prev()).norm_squared();
    auto b = 2 * ((other.velocity_prev(0) - velocity_prev(0)) * (other.position_prev(0) - position_prev(0)) + (other.velocity_prev(1) - velocity_prev(1)) * (other.position_prev(1) - position_prev(1)));
    auto c = (position_prev() - other.position_prev()).norm_squared() - std::pow(radius() + other.radius(), 2);
    
    auto t1 = (-b + std::sqrt(b*b - 4 * a * c)) / (2 * a);
    auto t2 = (-b - std::sqrt(b*b - 4 * a * c)) / (2 * a);
    
    if (t1 > 0 && t2 > 0)
    {
        if (t1 > t2)
        {
            return t2;
        }
        return t1;
    }
    else if (t1 > 0)
    {
        return t1;
    }
    return t2;
}

void Particle::collision_acceleration(Particle& other, double dt, double epsilon)
{
    auto n = (other.position() - position()) / (other.position() - position()).norm();
    auto dv = dot(other.velocity() - velocity(), n);
    auto m_ratio = other.mass() / (other.mass() + mass());
    auto u_prime = u + n * dot(other.velocity() - velocity(), n) * (1 + epsilon) * other.mass() / (other.mass() + mass());
    a_collision = (u_prime - u) / dt;
}

auto Particle::apply_accelerations(Vector u_p,
                            double dt,
                            std::vector<double> walls,
                            double epsilon,
                            Vector g,
                            bool drag,
                            bool grav,
                            bool wall,
                            bool collision)
{
    Vector acceleration;
    if (drag) { acceleration += drag_acceleration(u_p); }
    if (grav) { acceleration += buoyancy_acceleration(g); }
    if (wall) { acceleration += wall_contact_acceleration(dt, epsilon, walls); }
    if (collision) { acceleration += a_collision; }
    return acceleration;
}

void Particle::integrate(double dt,
                         std::vector<double> walls,
                         double epsilon,
                         Vector g,
                         bool drag,
                         bool grav,
                         bool wall,
                         bool collision /* = false */)
{
    // fix apply u and x simultaneously
    auto a1 = apply_accelerations(u_n,                 dt, walls, epsilon, g, drag, grav, wall, collision);
//    auto a2 = apply_accelerations(u_n + dt / 2.0 * a1, dt, walls, epsilon, g, drag, grav, wall, collision);
//    auto a3 = apply_accelerations(u_n + dt / 2.0 * a2, dt, walls, epsilon, g, drag, grav, wall, collision);
//    auto a4 = apply_accelerations(u_n + dt * a3,       dt, walls, epsilon, g, drag, grav, wall, collision);
    
//    u = u_n + dt / 6.0 * (a1 + 2*a2 + 2*a3 + a4);
//
//    x = x_n + (dt + std::pow(dt,2) / 2 + std::pow(dt,3) / 6 + std::pow(dt,4) / 24) * u;
    
    u = u_n + dt * a1;
    x = x_n + dt * u;
}

void Particle::update()
{
    x_n = x;
    u_n = u;
};

void Particle::integrate_update(double dt,
                                std::vector<double> walls,
                                double epsilon,
                                Vector g,
                                bool drag,
                                bool grav,
                                bool wall,
                                bool collision /* = false */)
{
    integrate(dt, walls, epsilon, g, drag, grav, wall, collision);
    update();
}
