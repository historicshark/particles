#pragma once

#include <cmath>
#include <string>
#include <vector>

#include "Vector2D.hpp"

Vector interpolate_flow_properties_fast(Vector x, std::string flow_type, std::vector<double> parameters);

double curl_of_velocity(Vector x, std::string flow_type, std::vector<double> parameters);
