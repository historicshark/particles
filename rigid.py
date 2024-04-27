import numpy as np


def acceleration_particle_collision_rigid(particle_collision,
                                          particle_collision_loc,
                                          position,
                                          velocity,
                                          radius,
                                          mass,
                                          mu_l,
                                          mu_p,
                                          sigma,
                                          epsilon,
                                          dt):
    a = np.zeros((np.sum(particle_collision_loc), 2))
    for i, (x, u, m, loc) in enumerate(zip(position[particle_collision_loc], velocity[particle_collision_loc], mass[particle_collision_loc], particle_collision)):
        distance_vector = position[loc] - x
        n = distance_vector / np.sqrt(np.sum(distance_vector * distance_vector, axis=1))[:,np.newaxis]
        u_prime = u + n * (np.sum((velocity[loc] - u) * n, axis=1) * (1 + epsilon) * mass[loc] / (mass[loc] + m))[:,np.newaxis]
        a[i] = np.sum((u_prime - u) / dt, axis=0)
    return a


def acceleration_wall_collision_rigid(wall_collision,
                                      wall_collision_loc,
                                      wall_overlap,
                                      position,
                                      velocity,
                                      mass,
                                      epsilon,
                                      dt):
    wall_collision_rigid_loc = np.any(np.hsplit(wall_collision, 2), axis=2).transpose()
    a = -(1 + epsilon) * velocity[wall_collision_rigid_loc] / dt
    return a
