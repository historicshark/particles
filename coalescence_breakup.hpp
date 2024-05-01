#pragma once

#include <cmath>
#include <random>
#include <tuple>

#include "Vector2D.hpp"

const double u_distribution_factor = std::tgamma(1.) / (std::tgamma(0.5) * std::tgamma(0.5));

const double log_factor = std::log(10000);

bool detect_breakup(double we, double eo);

double radius_fraction();

std::tuple<double, double, Vector> breakup_radii_and_new_position(double r, Vector x);

double equivalent_diameter(double r1, double r2);

double drainage_time(double deq, double rho_l, double sigma);

std::tuple<double, Vector> coalesced_bubble_radius_and_velocity(double r1, double r2, Vector u1, Vector u2);

bool detect_coalescence(double r1, double r2, double rho_l, double sigma, double t_contact);
