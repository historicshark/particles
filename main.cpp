#include <iostream>
#include <cmath>
#include <random>
#include <vector>
#include <array>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>

#include "Vector2D.hpp"
#include "properties.hpp"
#include "util.hpp"

#include "simple.hpp"


void calculate(double dt, 
               double end_time, 
               int n_frames, 
               double rho_l, 
               double mu_l, 
               double mu_p, 
               double epsilon, 
               Vector g, 
               double sigma, 
               std::array<double, 4> walls,
               std::string flow_type, 
               std::vector<double> parameters)
{
    std::vector<double> radius;
    std::vector<double> rho_p;
    std::vector<Vector> x;
    std::vector<Vector> u;
    
    auto read_1d_csv = [](std::string filename, std::vector<double>& list)    
    {
        std::ifstream radius_file(filename);
        if (radius_file.is_open())
        {
            std::string line;
            while (std::getline(radius_file, line))
            {
                list.push_back(std::stod(line));
            }
        }
    };
    
    auto read_2d_csv = [](std::string filename, std::vector<Vector>& list)
    {
        std::ifstream x_file(filename);
        if (x_file.is_open())
        {
            std::string line;
            while (std::getline(x_file, line))
            {
                std::stringstream ss(line);
                std::string entry;
                std::getline(ss, entry, ',');
                double vx = std::stod(entry);
                std::getline(ss, entry, ',');
                double vy = std::stod(entry);
                list.push_back({vx, vy});
            }
        }
    };
    
    read_1d_csv("radius.csv", radius);
    read_1d_csv("rho_p.csv", rho_p);
    read_2d_csv("x.csv", x);
    read_2d_csv("u.csv", u);
    
    auto x_n = x;
    auto u_n = u;
    
    auto t_save = linspace(0., end_time, n_frames);
    int i_save = 0;
    for (double t = 0; t != end_time; t += dt)
    {
        std::vector<double> mass;
        std::transform(radius.cbegin(), radius.cend(), rho_p.cbegin(), std::back_inserter(mass), [&rho_l](const double& r, const double& rho){return calculate_mass(r, rho, rho_l);});
        
    }
}

int main(int argc, char* argv[])
{
    double dt, end_time, gx, gy, mu_l, rho_l, sigma, mu_p, xmin, xmax, ymin, ymax;
    int n_frames;
    bool drag, gravity, particle_collisions;
    std::string flow_type;
    std::vector<double> parameters;
    
    std::ifstream options_file("options.txt");
    if (options_file.is_open())
    {
        std::string line;
        std::getline(options_file, line);
        dt = std::stod(line);
        std::getline(options_file, line);
        end_time = std::stod(line);
        std::getline(options_file, line);
        n_frames = std::stoi(line);
        std::getline(options_file, line);
        gx = std::stod(line);
        std::getline(options_file, line);
        gy = std::stod(line);
        std::getline(options_file, line);
        mu_l = std::stod(line);
        std::getline(options_file, line);
        rho_l = std::stod(line);
        std::getline(options_file, line);
        sigma = std::stod(line);
        std::getline(options_file, line);
        mu_p = std::stod(line);
        std::getline(options_file, line);
        xmin = std::stod(line);
        std::getline(options_file, line);
        xmax = std::stod(line);
        std::getline(options_file, line);
        ymin = std::stod(line);
        std::getline(options_file, line);
        ymax = std::stod(line);
        
        std::getline(options_file, line);
        drag = line != "0";
        std::getline(options_file, line);
        gravity = line != "0";
        std::getline(options_file, line);
        particle_collisions = line != "0";
        
        std::getline(options_file, flow_type);
        
        std::getline(options_file, line);
        std::stringstream ss(line);
        std::string entry;
        while (std::getline(ss, entry))
        {
            parameters.push_back(std::stod(entry));
        }
    }
    
    Vector g = {gx, gy};
    std::array<double,4> walls = {xmin, xmax, ymin, ymax};
    
    calculate(dt, end_time, n_frames, rho_l, mu_l, mu_p, 1, g, sigma, walls, flow_type, parameters);
}