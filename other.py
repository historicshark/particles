import numpy as np


def magnitude(a):
    return np.sqrt(np.sum(a*a, axis=1))


def drag_coefficient_none(re, eo):
    return np.zeros_like(re)


def acceleration_particle_collision_none(particle_collision,
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
    return np.zeros((np.sum(particle_collision_loc), 2))
