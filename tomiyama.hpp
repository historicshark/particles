#pragma once


#include <algorithm>
#include <cmath>

double drag_coefficient_tomiyama_pure(double re, double eo);
double drag_coefficient_tomiyama_low(double re, double eo);
double drag_coefficient_tomiyama_high(double re, double eo);
double lift_coefficient_tomiyama(double re, double eo_d, double radius);
