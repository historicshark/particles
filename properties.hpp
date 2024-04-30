#pragma once

#include <vector>
#include <cmath>
#include <array>
#include <numbers>

#include "Vector2D.hpp"
#include "util.hpp"

double calculate_area_xs(double radius);
double calculate_volume(double radius);
double calculate_mass(double radius, double rho_p, double rho_l);
double reynolds_number(Vector u_rel, double radius, double rho_l, double mu_l);
double eotvos_number(Vector g, double rho_p, double rho_l, double radius, double sigma);
double weber_number(double rho_l, Vector u_rel, double radius, double sigma);
double kinetic_energy(Vector u, double m);
double potential_energy(double m, Vector g, Vector x, double ymin);
double mechanical_energy(Vector x, Vector u, double m, Vector g, double ymin);
