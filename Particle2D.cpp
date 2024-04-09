auto Particle::distance_vector(Particle other) { return other.x - x; }
auto Particle::distance_vector(Vector other) { return other - x; }

auto Particle::distance(Particle other) { return distance_vector(other).norm(); }
auto Particle::distance(Vector other) { return distance_vector(other).norm(); }

void Particle::set_fluid_properties(double rho, double mu)
{
    rho_l = rho;
    mu_l = mu;
}

void Particle::update_reynolds_number(Vector u_p, double rho_l, double mu_l, Vector v_l = {0,0})
{
    re = rho_l * (u_p - v_l).norm() * diameter() / mu_l;
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

auto Particle::drag_force(Vector u_p, double rho_l, double mu_l, Vector v_l = {0,0})
{
    update_reynolds_number(u_p, rho_l, mu_l, v_l);
    return 1.0 / 2.0 * rho_l / rho * (v_l - u_p).norm() * (v_l - u_p) * area_xs() / volume() * drag_coefficient();
};

auto Particle::buoyancy_force(double rho_l, Vector g = {0, -9.81})
{
    return (rho - rho_l) / rho * g;
}

auto Particle::wall_contact_force(const double dt, double epsilon, std::vector<double> walls)
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

auto Particle::apply_forces(Vector u_p, 
                            double dt, 
                            double rho_l, 
                            double mu_l, 
                            std::vector<double> walls, 
                            double epsilon, 
                            Vector v_l,
                            Vector g,
                            bool drag,
                            bool grav,
                            bool wall)
{
    Vector force;
    if (drag) { force += drag_force(u_p, rho_l, mu_l, v_l); }
    if (grav) { force += buoyancy_force(rho_l, g); }
    if (wall) { force += wall_contact_force(dt, epsilon, walls); }
    return force;
}

void Particle::update(double dt, 
                      double rho_l, 
                      double mu_l, 
                      std::vector<double> walls, 
                      double epsilon,
                      Vector v_l, 
                      Vector g,
                      bool drag,
                      bool grav,
                      bool wall)
{
    auto k1 = apply_forces(u_n,                 dt, rho_l, mu_l, walls, epsilon, v_l, g, drag, grav, wall);
    auto k2 = apply_forces(u_n + dt / 2.0 * k1, dt, rho_l, mu_l, walls, epsilon, v_l, g, drag, grav, wall);
    auto k3 = apply_forces(u_n + dt / 2.0 * k2, dt, rho_l, mu_l, walls, epsilon, v_l, g, drag, grav, wall);
    auto k4 = apply_forces(u_n + dt * k3,       dt, rho_l, mu_l, walls, epsilon, v_l, g, drag, grav, wall);
    
    u = u_n + dt / 6.0 * (k1 + 2*k2 + 2*k3 + k4);
    
    x = x_n + (dt + std::pow(dt,2) / 2 + std::pow(dt,3) / 6 + std::pow(dt,4) / 24) * u;
    
    x_n = x;
    u_n = u;
};
