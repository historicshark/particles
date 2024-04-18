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
    std::string results_filepath = argv[i_arg++];
    
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
                p.set_gravity(g);
                particles.push_back(p);
            }
        }
    }
    
    std::ofstream time_file(time_filepath);
    std::ofstream position_file(position_filepath);
    std::ofstream velocity_file(velocity_filepath);
    std::ofstream results_file(results_filepath);
    
    results_file << "ke,pe,me\n";
    
    double dt_substeps_before_collision = dt / 10.;
//    int n_steps = end_time / dt;
//    int step_print = n_steps / 100;
//    int i = 0;
//    double time = 0;
    
    for (double t = 0; t < end_time; t += dt)
    {
//        std::cout << "step" << std::endl;
//        if (i == step_print)
//        {
//            std::cout << "time = " << time << std::endl;
//            i = 0;
//        }
        time_file << t << "\n";
        
        // Reset collision flag
        std::for_each(particles.begin(), particles.end(), [](Particle& p){p.set_colllided(false);});
        
        for (auto& particle : particles)
        {
            if (!particle.collided())
            {
                // Calculate next position
                particle.integrate(dt, walls, epsilon, g, drag, grav, wall);
                
                // Check to see if particle has collided with any other particle
                for (auto& other : particles)
                {
                    if (particle != other && !other.collided())
                    {
                        other.integrate(dt, walls, epsilon, g, drag, grav, wall);
                        
                        auto distance = particle.distance(other);
                        if (particle.distance(other) < particle.radius() + other.radius())
                        {
                            particle.integrate_update(dt_substeps_before_collision, walls, epsilon, g, drag, grav, wall);
                            other.integrate_update(dt_substeps_before_collision, walls, epsilon, g, drag, grav, wall);
                            
                            auto min_distance = particle.radius() + other.radius();
                            double time_before_collision = 0;
                            
                            distance = particle.distance(other);
                            while (particle.distance(other) > min_distance)
                            {
                                particle.integrate_update(dt_substeps_before_collision, walls, epsilon, g, drag, grav, wall);
                                other.integrate_update(dt_substeps_before_collision, walls, epsilon, g, drag, grav, wall);
                                time_before_collision += dt_substeps_before_collision;
                                distance = particle.distance(other);
                            }
                            
                            double dt_substep_after_collision = (dt - time_before_collision);
                            
                            particle.collision_acceleration(other, dt_substep_after_collision, epsilon);
                            other.collision_acceleration(particle, dt_substep_after_collision, epsilon);
                            
                            particle.integrate_update(dt_substep_after_collision, walls, epsilon, g, drag, grav, wall, /* collision = */ true);
                            other.integrate_update(dt_substep_after_collision, walls, epsilon, g, drag, grav, wall, /* collision = */ true);
                            
                            // These two particles have integrated through dt so do not integrate them again
                            particle.set_colllided(true);
                            other.set_colllided(true);
                        }
                        else
                        {
                            
                        }
                    }
                }
                
//                if (!particle.collided())
//                {
//                    particle.update();
//                }
                
                particle.update();
            }
            
            position_file << particle.position_string() << ",";
            velocity_file << particle.velocity_string() << ",";
        }
        position_file << "\n";
        velocity_file << "\n";
        
        auto ke = std::accumulate(particles.begin(), particles.end(), 0.0, [] (double ke, Particle& particle) { return ke + particle.kinetic_energy(); });
        auto pe = std::accumulate(particles.begin(), particles.end(), 0.0, [ymin] (double pe, Particle& particle) { return pe + particle.potential_energy(ymin); });
        auto me = std::accumulate(particles.begin(), particles.end(), 0.0, [ymin] (double me, Particle& particle) { return me + particle.mechanical_energy(ymin); });
        results_file << ke << "," << pe << "," << me << "\n";
//        i++;
//        time += dt;
    }
    time_file.close();
    position_file.close();
    velocity_file.close();
}
