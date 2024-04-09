#include <cmath>
#include <iostream>
#include <fstream>
#include <cstring>
#include <string>
#include <random>
#include <numeric>
#include <vector>

#include "define.h"
#include "Vector2D.h"
#include "Particle2D.h"

int main(int argc, char* argv[])
{
    int i_arg = 1;
    std::string domain_filepath = argv[i_arg++];
    std::string particle_filepath = argv[i_arg++];
    std::string time_filepath = argv[i_arg++];
    std::string position_filepath = argv[i_arg++];
    std::string velocity_filepath = argv[i_arg++];
    
    double gx = std::stod(argv[i_arg++]);
    double gy = std::stod(argv[i_arg++]);
    double rho_l = std::stod(argv[i_arg++]);
    double mu_l = std::stod(argv[i_arg++]);
    double epsilon = std::stod(argv[i_arg++]);
    double xmin = std::stod(argv[i_arg++]);
    double xmax = std::stod(argv[i_arg++]);
    double ymin = std::stod(argv[i_arg++]);
    double ymax = std::stod(argv[i_arg++]);
    double dt = std::stod(argv[i_arg++]);
    double end_time = std::stod(argv[i_arg++]);
    
    bool drag = false;
    if (strcmp(argv[i_arg++], "True") == 0) { drag = true; }
    bool grav = false;
    if (strcmp(argv[i_arg++], "True") == 0) { grav = true; }
    bool wall = false;
    if (strcmp(argv[i_arg++], "True") == 0) { wall = true; }
    
    Vector g{gx, gy};
    Vector v_l;
    
    std::vector<double> walls {xmax, xmin, ymax, ymin};
    
    std::vector<Particle> particles;
    {
        std::ifstream particle_file(particle_filepath);
        if (particle_file.is_open())
        {
            std::string line;
            while (std::getline(particle_file, line))
            {
                std::array<std::string, 7> in;
                std::stringstream ss(line);
                std::string entry;
                for (int i = 0; i != 7; i++)
                {
                    std::getline(ss, entry, ',');
                    in[i] = entry;
                }
                Particle p = {in};
                p.set_fluid_properties(rho_l, mu_l);
                p.set_fluid_velocity(v_l);
                particles.push_back(p);
            }
        }
    }
    
    std::ofstream time_file(time_filepath);
    std::ofstream position_file(position_filepath);
    std::ofstream velocity_file(velocity_filepath);
    
    for (double t = 0; t < end_time; t += dt)
    {
        time_file << t << "\n";
        for (auto& particle : particles)
        {
            particle.update(dt, walls, epsilon, g, drag, grav, wall);
            position_file << particle.position_string() << ",";
            velocity_file << particle.velocity_string() << ",";
        }
        position_file << "\n";
        velocity_file << "\n";
    }
    time_file.close();
    position_file.close();
    velocity_file.close();
}
