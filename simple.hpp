#pragma once

#include <cmath>
#include <array>
#include <numbers>

#include <range/v3/view/zip.hpp>

#include "Vector2D.hpp"
#include "util.hpp"

double fluid_film_thickness(double mu_l, double mu_p, double radius, double sigma, bool is_particle_collision);
double contact_radius(double radius, double deformation);
Vector force_elastic_simple(double radius, double deformation, double sigma, Vector n);
Vector force_viscous_simple(double deformation, double h0, double radius, Vector velocity, double mu_l, Vector n, bool is_particle_collision);
Vector acceleration_particle_collision_simple(Vector position, Vector other_position, Vector velocity, double radius, double other_radius, double mass, double mu_l, double mu_p, double sigma);
Vector acceleration_wall_collision_simple(std::array<bool, 4> wall_collision, Walls wall_overlap, Vector position, Vector velocity, double radius, double mass, double mu_p, double mu_l, double sigma);
