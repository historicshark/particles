#pragma once

#include <array>
#include <cmath>
#include <stdexcept>
#include <functional>
#include <numbers>

#include "define.h"
#include "Vector2D.h"

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
    Vector v_l;
    
    double re;
    
    Vector a_collision;
    bool collision_done;
    
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
    auto position_prev() { return x_n; };
    auto position(const int dim)
    {
        if (dim == 0) { return x[0]; }
        else if (dim == 1) { return x[1]; }
        else { throw std::invalid_argument("ndim"); };
    }
    auto position_prev(const int dim)
    {
        if (dim == 0) { return x_n[0]; }
        else if (dim == 1) { return x_n[1]; }
        else { throw std::invalid_argument("ndim"); };
    }
    auto velocity() { return u; };
    auto velocity_prev() { return u_n; };
    auto velocity(const int dim)
    {
        if (dim == 0) { return u[0]; }
        else if (dim == 1) { return u[1]; }
        else { throw std::invalid_argument("ndim"); };
    }
    auto velocity_prev(const int dim)
    {
        if (dim == 0) { return u_n[0]; }
        else if (dim == 1) { return u_n[1]; }
        else { throw std::invalid_argument("ndim"); };
    }
    auto volume() { return 4./3. * std::numbers::pi * pow(r, 3); };
    auto area_xs() { return std::numbers::pi * pow(r, 2); };
    auto mass() { return rho * volume(); };
    auto collided() { return collision_done; };
    void set_colllided(bool a) { collision_done = a; };
    
    bool operator==(Particle& other) { return get_id() == other.get_id(); };
    bool operator!=(Particle& other) { return get_id() != other.get_id(); };
    
    auto distance_vector(Particle& other);
    auto distance_vector(Vector other);
    auto distance(Particle& other);
    auto distance(Vector other);
    void set_fluid_properties(double rho, double mu);
    void set_fluid_velocity(Vector v) { v_l = v; };
    void update_reynolds_number(Vector u_p);
    auto reynolds_number() { return re; };
    auto kinetic_energy();
    auto drag_coefficient();
    auto drag_acceleration(Vector u_p);
    auto buoyancy_acceleration(Vector g);
    auto wall_contact_acceleration(const double dt, double epsilon, std::vector<double> walls);
    auto time_to_collision(Particle& other);
    void collision_acceleration(Particle& other, double dt, double epsilon);
    auto apply_accelerations(Vector u_p, 
                      double dt,
                      std::vector<double> walls,
                      double epsilon,
                      Vector g,
                      bool drag,
                      bool grav,
                      bool wall,
                      bool collision);
    void integrate(double dt,
                   std::vector<double> walls,
                   double epsilon,
                   Vector g,
                   bool drag,
                   bool grav,
                   bool wall,
                   bool collision=false);
    void update();
    void integrate_update(double dt,
                          std::vector<double> walls,
                          double epsilon,
                          Vector g,
                          bool drag,
                          bool grav,
                          bool wall,
                          bool collision=false);
    auto position_string() { return x.string(); };
    auto velocity_string() { return u.string(); };
};

#include "Particle2D.cpp"
