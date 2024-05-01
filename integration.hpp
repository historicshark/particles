#pragma once

#include <string>
#include <range/v3/view/zip.hpp>

#include "Vector2D.hpp"
#include "util.hpp"
#include "tomiyama.hpp"
#include "simple.hpp"
#include "properties.hpp"
#include "flow_properties.hpp"

Walls wall_overlap(Vector x, double radius, Walls walls);

WallContact detect_wall_contact(Walls overlap);

bool detect_particle_contact(Vector position1, Vector position2, double radius1, double radius2);

std::vector<std::tuple<size_t, size_t>> new_particle_contact_list(const std::vector<Vector>& position, const std::vector<double>& radius);

Vector acceleration_drag(Vector particle_velocity,
                         Vector flow_velocity,
                         double radius,
                         double rho_l,
                         double mu_l,
                         double mass,
                         Vector g,
                         double rho_p,
                         double sigma);

Vector acceleration_gravity(double rho_p, double rho_l, Vector g);

//Vector acceleration_lift();

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
                                        bool wall_collisions);

std::vector<Vector> apply_velocities(std::vector<Vector>& velocity, double dt);

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
                   bool wall_collisions);

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
                     bool wall_collisions);
