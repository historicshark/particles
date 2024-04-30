#include <iostream>
#include <cmath>
#include <random>
#include <vector>
#include <array>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>

#include <range/v3/view/zip.hpp>

#include "Vector2D.hpp"
#include "integration.hpp"
#include "properties.hpp"
#include "util.hpp"
#include "tomiyama.hpp"
#include "flow_properties.hpp"
#include "simple.hpp"


void update_particle_contact(const std::vector<Vector>& position, const std::vector<double>& radius, std::vector<std::array<size_t, 2>>& contact_list, std::vector<double>& contact_time_list, double dt)
{
    for (size_t i = 0; i != position.size(); i++)
    {
        for (size_t j = i + 1; j != position.size(); j++)
        {
            std::array<size_t, 2> elem = {i,j};
            auto in_list_it = std::find(contact_list.begin(), contact_list.end(), elem);
            bool in_list = in_list_it != contact_list.end();
            bool are_colliding = detect_particle_contact(position[i], position[j], radius[i], radius[j]);
            if (in_list && are_colliding)
            {
                size_t loc = std::distance(contact_list.begin(), in_list_it);
                if (are_colliding)
                {
                    contact_time_list[loc] += dt;
                }
                else
                {
                    contact_list[loc] = contact_list.back();
                    contact_time_list[loc] = contact_time_list.back();
                    contact_list.pop_back();
                    contact_time_list.pop_back();
                }
            }
            else if (!in_list && are_colliding)
            {
                contact_list.push_back(elem);
                contact_time_list.push_back(0);
            }
        }
    }
}

void calculate(double dt, 
               double end_time, 
               int n_frames,
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
               bool wall_collisions,
               bool particle_collisions)
{
    std::cout << "_________________________\n";
    std::vector<double> radius;
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
    
    read_1d_csv("r.csv", radius);
    read_2d_csv("x.csv", x);
    read_2d_csv("u.csv", u);
    
    auto x_n = x;
    auto u_n = u;
    
    auto t_save = linspace(0., end_time, n_frames);
    auto t_save_it = t_save.cbegin();
    
    std::ofstream position_file("position.csv");
    std::ofstream velocity_file("velocity.csv");
    std::ofstream results_file("results.csv");
    std::ofstream radius_file("radius.csv");
    std::ofstream time_file("t.csv");
    
    std::vector<std::array<size_t, 2>> contact_list;
    std::vector<double> contact_time_list;
    
    std::vector<Vector> flow_velocity(u.size());
    std::transform(x.cbegin(), x.cend(), flow_velocity.begin(), [flow_type, parameters](const Vector& v){return interpolate_flow_properties_fast(v, flow_type, parameters);});
    
    std::vector<double> mass(u.size());
    std::transform(radius.cbegin(), radius.cend(), mass.begin(), [rho_p, rho_l](const double& r){return calculate_mass(r, rho_p, rho_l);});
    
    for (double t = 0; t < end_time; t += dt)
    {
        std::transform(x.begin(), x.end(), flow_velocity.begin(), [flow_type, parameters](Vector& v){return interpolate_flow_properties_fast(v, flow_type, parameters);});
        
        update_particle_contact(x, radius, contact_list, contact_time_list, dt);
        
        integrate_rk4(x_n,
                      u_n,
                      x,
                      u,
                      dt,
                      flow_velocity,
                      radius,
                      mass,
                      rho_p,
                      rho_l,
                      mu_l,
                      mu_p,
                      g,
                      sigma,
                      walls,
                      drag,
                      gravity,
                      particle_collisions,
                      wall_collisions);

        if (t > *t_save_it)
        {
            t_save_it++;
            std::cout << "t = " << t << "\t" << contact_list.size() << " " << x[0].string() << std::endl;
            
            std::for_each(x.cbegin(), x.cend(), [&position_file](const Vector& v){position_file << v.string() << ",";});
            position_file << "\n";
            
            std::for_each(u.cbegin(), u.cend(), [&velocity_file](const Vector& v){velocity_file << v.string() << ",";});
            velocity_file << "\n";
            
            std::for_each(radius.cbegin(), radius.cend(), [&radius_file](const double& r){radius_file << r << ",";});
            radius_file << "\n";
            
            time_file << t << "\n";
            
            double ke = 0;
            double pe = 0;
            
            for (auto [pos, vel, m] : ranges::views::zip(x, u, mass))
            {
                ke += kinetic_energy(vel, m);
                pe += potential_energy(m, g, pos, walls[2]);
            }
            
            results_file << ke << "," << pe << "\n";
            
        }
        
        x_n = x;
        u_n = u;
    }
    std::cout << "end" << std::endl;
//    std::ranges::views::z
//    for (auto [pos, vel] : ranges::views::zip(x, u))
//    {
//        std::cout << pos.string() << "\t" << vel.string() << std::endl;
//        pos = {0,0};
//    }
//    for (auto [pos, vel] : ranges::views::zip(x, u))
//    {
//        std::cout << pos.string() << "\t" << vel.string() << std::endl;
//    }
}

int main(int argc, char* argv[])
{
    std::cout << "main" << std::endl;
    double dt, end_time, gx, gy, mu_l, rho_l, sigma, rho_p, mu_p, xmin, xmax, ymin, ymax;
    int n_frames;
    bool drag, gravity, wall_collisions, particle_collisions;
    std::string flow_type;
    std::vector<double> parameters;
    
    std::ifstream options_file("options.txt");
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
    rho_p = std::stod(line);
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
    wall_collisions = line != "0";
    std::getline(options_file, line);
    particle_collisions = line != "0";
    
    std::getline(options_file, flow_type);
    
    std::getline(options_file, line);
    std::stringstream ss(line);
    std::string entry;
    while (std::getline(ss, entry, ','))
    {
        parameters.push_back(std::stod(entry));
    }
    
    Vector g = {gx, gy};
    Walls walls = {xmin, xmax, ymin, ymax};
    
    calculate(dt,
              end_time,
              n_frames,
              rho_p,
              rho_l,
              mu_l,
              mu_p,
              g,
              sigma,
              walls,
              flow_type,
              parameters,
              drag,
              gravity,
              wall_collisions,
              particle_collisions);
}
