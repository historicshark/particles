#pragma once

#include <cmath>

#include "Vector2D.hpp"
#include "properties.hpp"
#include "simple.hpp"

Vector acceleration_drag(Vector particle_velocity, Vector flow_velocity, double radius, double rho_l, double mu_l, double mass, Vector g, double rho_p, double sigma);

Vector acceleration_gravity(double rho_p, double rho_l, Vector g);

Vector apply_accelerations(Vector position,
                           Vector velocity,
                           Vector flow_velocity,
                           double radius,
                           double rho_p,
                           double mass,
                           double rho_l,
                           double mu_l,
                           double mu_p,
                           Vector g,
                           double sigma,
                           std::array<double,4> walls,
                           double dt);

