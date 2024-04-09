#pragma once

#include <array>
#include <cmath>
#include <stdexcept>
#include <functional>
#include <numbers>

#include "define.h"

class Particle
{
    int id;
    
    // Current timestep
    Vector x;
    Vector u;
    
    // Previous timestep
    Vector x_n;
    Vector u_n;
    
    double r;
    double rho;
    
    double rho_l;
    double mu_l;
    
    double re;
        
    public:
    Particle(int id, double r, double rho, double x, double y, double u, double v, double rho_l, double mu_l)
    : id(id)
    , r(r)
    , rho(rho)
    , x({x,y})
    , u({u,v})
    , x_n({x,y})
    , u_n({u,v})
    , rho_l(rho_l)
    , mu_l(mu_l)
    {};
    
    Particle(int id, double r, double rho, Vector x, Vector u, double rho_l, double mu_l)
    : id(id)
    , r(r)
    , rho(rho)
    , x(x)
    , u(u)
    , x_n(x)
    , u_n(u)
    , rho_l(rho_l)
    , mu_l(mu_l)
    {};
    
    Particle(std::array<std::string,7> line)
    : id(std::stoi(line[0]))
    , r(std::stod(line[1]))
    , rho(std::stod(line[2]))
    , x({std::stod(line[3]),std::stod(line[4])})
    , u({std::stod(line[5]),std::stod(line[6])})
    {
        x_n = x;
        u_n = u;
    };
    
    auto get_id() { return id; };
    auto diameter() { return 2 * r; };
    auto radius() { return r; };
    auto position() { return x; };
    auto position(const int dim)
    {
        if (dim == 0) { return x[0]; }
        else if (dim == 1) { return x[1]; }
        else { throw std::invalid_argument("ndim"); };
    }
    auto velocity() { return u; };
    auto velocity(const int dim)
    {
        if (dim == 0) { return u[0]; }
        else if (dim == 1) { return u[1]; }
        else { throw std::invalid_argument("ndim"); };
    }
    auto volume() { return 4./3. * std::numbers::pi * pow(r, 3); };
    auto area_xs() { return std::numbers::pi * pow(r, 2); };
    
    auto distance_vector(Particle other);
    auto distance_vector(Vector other);
    auto distance(Particle other);
    auto distance(Vector other);
    void set_fluid_properties(double rho, double mu);
    void update_reynolds_number(Vector u_p, double rho_l, double mu_l, Vector v_l);
    auto reynolds_number() { return re; };
    auto drag_coefficient();
    auto drag_force(Vector u_p, double rho_l, double mu_l, Vector v_l);
    auto buoyancy_force(double rho_l, Vector g);
    auto wall_contact_force(const double dt, double epsilon, std::vector<double> walls);
    // auto collision_force();
    auto apply_forces(Vector u_p, 
                      double dt, 
                      double rho_l, 
                      double mu_l, 
                      std::vector<double> walls, 
                      double epsilon, 
                      Vector v_l,
                      Vector g,
                      bool drag,
                      bool grav,
                      bool wall);
    void update(double dt, 
                double rho_l, 
                double mu_l, 
                std::vector<double> walls, 
                double epsilon,
                Vector v_l, 
                Vector g,
                bool drag,
                bool grav,
                bool wall);
    auto position_string() { return x.string(); };
    auto velocity_string() { return u.string(); };
};

#include "Particle2D.cpp"
