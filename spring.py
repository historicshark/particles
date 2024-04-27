import numpy as np

from other import magnitude


def acceleration_particle_collision_spring(particle_collision,
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
    k = 1e-4
    c = 1e-9
    for i, (x, u, r, m, loc) in enumerate(zip(position[particle_collision_loc],
                                              velocity[particle_collision_loc],
                                              radius[particle_collision_loc],
                                              mass[particle_collision_loc],
                                              particle_collision)):
        distance_vector = position[loc] - x
        distance = magnitude(distance_vector)
        n = distance_vector / np.sqrt(np.sum(distance_vector * distance_vector, axis=1))[:, np.newaxis]
        overlap = (r + radius[loc] - distance) / 2
        a_k = (k * overlap) / m
        a_c = (c * np.sum(u * n, axis=1)) / m
        a[i] = np.sum(-n * (a_k + a_c)[:,np.newaxis], axis=0)
    return a


def acceleration_wall_collision_spring(wall_collision,
                                       wall_collision_loc,
                                       wall_overlap,
                                       position,
                                       velocity,
                                       radius,
                                       mass,
                                       mu_p,
                                       mu_l,
                                       sigma,
                                       epsilon,
                                       dt):
    """
    wall_collision_loc: (n_contacting, 4)
    walls: (xmin, xmax, ymin, ymax)
    wall_overlap: (n_contacting, 4)
    """
    a = np.zeros((np.sum(wall_collision_loc), 2))
    k = 1e-3
    c = 1e-8
    n = np.array([(1,0), (-1,0), (0,1), (0,-1)])
    for i, (wc, loc, xov, u, m) in enumerate(zip(wall_collision,
                                                np.any(np.hsplit(wall_collision,2), axis=2).transpose(),
                                                wall_overlap,
                                                velocity[wall_collision_loc],
                                                mass[wall_collision_loc])):
        f_k = (k * xov[wc]) / m
        f_c = -(c * np.sum(np.broadcast_to(u, n[wc].shape) * n[wc], axis=1)) / m
        a[i] = np.sum(n[wc] * (f_k + f_c), axis=0)
    return a
