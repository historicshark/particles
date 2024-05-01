#include <iostream>
#include <cmath>
#include <random>
#include <vector>
#include <array>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <tuple>

#include <range/v3/view/zip.hpp>

#include "Vector2D.hpp"
#include "integration.hpp"
#include "properties.hpp"
#include "util.hpp"
#include "tomiyama.hpp"
#include "flow_properties.hpp"
#include "simple.hpp"
#include "coalescence_breakup.hpp"


void update_particle_contact(const std::vector<Vector>& position,
                             const std::vector<double>& radius,
                             std::vector<std::tuple<size_t, size_t>>& contact_list,
                             std::vector<double>& contact_time_list,
                             double dt)
{
    for (size_t i = 0; i != position.size(); i++)
    {
        for (size_t j = i + 1; j != position.size(); j++)
        {
            auto elem = std::make_tuple(i,j);
            auto in_list_it = std::find(contact_list.begin(), contact_list.end(), elem);
            bool in_list = in_list_it != contact_list.end();
            bool are_colliding = detect_particle_contact(position[i], position[j], radius[i], radius[j]);
            if (in_list)
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
            else if (are_colliding)
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
               bool particle_collisions,
               bool lift)
{
    std::cout << "-------------------------\n";
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
    
    auto t_loadbar = linspace(0., end_time, 25);
    auto t_loadbar_it = t_loadbar.cbegin();
    
    std::ofstream position_file("position.csv");
    std::ofstream velocity_file("velocity.csv");
    std::ofstream results_file("results.csv");
    std::ofstream radius_file("radius.csv");
    std::ofstream time_file("time.csv");
    
    std::vector<std::tuple<size_t, size_t>> contact_list;
    std::vector<double> contact_time_list;
    
    std::vector<double> mass(u.size());
    std::transform(radius.cbegin(), radius.cend(), mass.begin(), [rho_p, rho_l](const double& r){
        return calculate_mass(r, rho_p, rho_l);
    });
    
    for (double t = 0; t < end_time; t += dt)
    {
        std::vector<Vector> flow_velocity(x.size());
        std::transform(x.begin(), x.end(), flow_velocity.begin(), [flow_type, parameters](Vector& v){
            return interpolate_flow_properties_fast(v, flow_type, parameters);
        });
        
        // Evaluate breakup
        std::vector<double> weber(x.size());
        std::vector<double> eotvos(x.size());
        std::vector<double> new_radius;
        std::vector<double> new_mass;
        std::vector<Vector> new_position;
        std::vector<Vector> new_velocity;
        for (auto [we, eo, pos, vel, r, u_l] : ranges::views::zip(weber, eotvos, x, u, radius, flow_velocity))
        {
            we = weber_number(rho_l, vel - u_l, r, sigma);
            eo = eotvos_number(g, rho_p, rho_l, r, sigma);
            if (detect_breakup(we, eo))
            {
                auto [new_r1, new_r2, new_pos] = breakup_radii_and_new_position(r, pos);
                r = new_r1;
                new_radius.push_back(new_r2);
                new_position.push_back(new_pos);
                new_velocity.push_back(vel);
                new_mass.push_back(calculate_mass(new_r2, rho_p, rho_l));
            }
        }
        if (!new_radius.empty())
        {
            for (const auto [new_r, new_pos, new_vel, new_m] : ranges::views::zip(new_radius,
                                                                                  new_position,
                                                                                  new_velocity,
                                                                                  new_mass))
            {
                radius.push_back(new_r);
                x_n.push_back(new_pos);
                x.push_back(new_pos);
                u_n.push_back(new_vel);
                u.push_back(new_vel);
                mass.push_back(new_m);
            }
            
            update_particle_contact(x, radius, contact_list, contact_time_list, 0);
        }
        
        // Evaluate coalescence
        std::vector<size_t> particles_to_remove;
        for (auto [contact_pair, t_contact] : ranges::views::zip(contact_list, contact_time_list))
        {
            auto [i,j] = contact_pair;
            if (detect_coalescence(radius[i], radius[j], rho_l, sigma, t_contact))
            {
                auto [r_new, u_new] = coalesced_bubble_radius_and_velocity(radius[i], radius[j], u[i], u[j]);
                radius[i] = r_new;
                u[i] = u_new;
                u_n[i] = u_new;
                mass[i] = calculate_mass(r_new, rho_p, rho_l);
                particles_to_remove.push_back(j);
            }
        }
        if (!particles_to_remove.empty())
        {
            std::sort(particles_to_remove.begin(), particles_to_remove.end());
            for (auto it = particles_to_remove.rbegin(); it != particles_to_remove.rend(); it++)
            {
                x[*it] = x.back();
                x_n[*it] = x_n.back();
                u[*it] = u.back();
                u_n[*it] = u_n.back();
                radius[*it] = radius.back();
                mass[*it] = mass.back();
                
                x.pop_back();
                x_n.pop_back();
                u.pop_back();
                u_n.pop_back();
                radius.pop_back();
                mass.pop_back();
            }
        }
        
        update_particle_contact(x, radius, contact_list, contact_time_list, dt);
        
        integrate_rk4(x_n,
                      u_n,
                      x,
                      u,
                      dt,
                      radius,
                      mass,
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
                      particle_collisions,
                      wall_collisions,
                      lift);

        if (t > *t_save_it)
        {
            if (t > *t_loadbar_it)
            {
                std::cout << "#" << std::flush;
                t_loadbar_it++;
            }
            
            t_save_it++;
            
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
    std::for_each(x.cbegin(), x.cend(), [&position_file](const Vector& v){position_file << v.string() << ",";});
    position_file << "\n";
    
    std::for_each(u.cbegin(), u.cend(), [&velocity_file](const Vector& v){velocity_file << v.string() << ",";});
    velocity_file << "\n";
    
    std::for_each(radius.cbegin(), radius.cend(), [&radius_file](const double& r){radius_file << r << ",";});
    radius_file << "\n";
    
    time_file << end_time << "\n";
    
    double ke = 0;
    double pe = 0;
    
    for (auto [pos, vel, m] : ranges::views::zip(x, u, mass))
    {
        ke += kinetic_energy(vel, m);
        pe += potential_energy(m, g, pos, walls[2]);
    }
    
    results_file << ke << "," << pe << "\n";
    
    std::cout << "#\n-------------------------\n\n";
}

int main(int argc, char* argv[])
{
    double dt, end_time, gx, gy, mu_l, rho_l, sigma, rho_p, mu_p, xmin, xmax, ymin, ymax;
    int n_frames;
    bool drag, gravity, wall_collisions, particle_collisions, lift;
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
    std::getline(options_file, line);
    lift = line != "0";
    
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
              particle_collisions,
              lift);
}
