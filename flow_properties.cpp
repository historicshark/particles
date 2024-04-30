#include "flow_properties.hpp"

Vector interpolate_flow_properties_fast(Vector x, std::string flow_type, std::vector<double> parameters)
{
    Vector flow_velocity;
    if (flow_type == "uniform")
    {
        // parameters = [velocity-x]
        flow_velocity = {parameters[0], 0};
    }
    else if (flow_type == "vortex")
    {
        // parameters = [core-diameter, max-tangential-velocity]
        double r = x.norm();
        double theta = std::atan2(x[1], x[0]);
        double u_theta;
        if (r < parameters[0])
        {
            u_theta = r * parameters[1] / parameters[0];
        }
        else
        {
            u_theta = parameters[0] * parameters[1] / r;
        }
        flow_velocity = {-u_theta * std::sin(theta), u_theta * std::cos(theta)};
    }
    else if (flow_type == "test")
    {
        if (x[0] < 0)
        {
            flow_velocity = {parameters[0], 0};
        }
        else
        {
            flow_velocity = {-parameters[0], 0};
        }
    }
    return flow_velocity;
}
