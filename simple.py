import numpy as np

from other import magnitude


def fluid_film_thickness(mu_l, mu_p, radius, sigma, is_particle_collision):
    # assuming radius has already been indexed
    if is_particle_collision:
        factor = np.sqrt(3 / 16)
    else:
        factor = np.sqrt(3 / 8)

    return factor * np.sqrt(mu_l * mu_p * radius ** 2 / sigma)


def contact_radius(radius, deformation):
    return radius * (-30 + np.sqrt(900 + .424 * deformation / radius))


def force_elastic_simple(radius, deformation, sigma):
    return radius * sigma * (18.5 * (deformation / radius) ** 2 + 2 * deformation / radius)


def force_viscous_simple(deformation, h0, r, u, n, mu_l, is_particle_collision):
    r_a = contact_radius(r, deformation)
    if is_particle_collision:
        cbc = .25
    else:
        cbc = 1.
    return (np.sum(u * n, axis=1) * cbc * 6 * mu_l / np.pi
            * .34 * (deformation / r + .0002) ** (-.5)
            * (4 * np.sqrt(r ** 3 / h0) + 3 * r_a * r / h0))


def force_viscous_simple_alt(r, u, xov, n, mu_l):
    return -6 * np.pi * mu_l * np.sum(u * n, axis=1) * r * r / xov


def force_tangential_simple():
    pass


def acceleration_particle_collision_simple(particle_collision,
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
    for i, (x, u, r, m, loc) in enumerate(zip(position[particle_collision_loc],
                                              velocity[particle_collision_loc],
                                              radius[particle_collision_loc],
                                              mass[particle_collision_loc],
                                              particle_collision)):
        h0 = fluid_film_thickness(mu_l, mu_p, r, sigma, True)
        deformation = r - .5 * magnitude(position[loc] - x) + h0 / 2
        deformation[deformation < 0] = 0.0
        distance_vector = position[loc] - x
        n = distance_vector / np.sqrt(np.sum(distance_vector * distance_vector, axis=1))[:, np.newaxis]
        f_e = force_elastic_simple(r, deformation, sigma)
        f_v = force_viscous_simple(deformation, h0, r, u, n, mu_l, True)
        a[i] = np.sum(-n * (f_e + f_v)[:, np.newaxis], axis=0) / m
    return a


def acceleration_wall_collision_simple(wall_collision,
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
    a = np.zeros((np.sum(wall_collision_loc), 2))
    n = np.array([(1, 0), (-1, 0), (0, 1), (0, -1)])
    for i, (x, u, xov, r, m, loc) in enumerate(zip(position[wall_collision_loc],
                                                   velocity[wall_collision_loc],
                                                   wall_overlap,
                                                   radius[wall_collision_loc],
                                                   mass[wall_collision_loc],
                                                   wall_collision)):
        h0 = fluid_film_thickness(mu_l, mu_p, r, sigma, False)
        deformation = xov[loc] + h0
        f_e = force_elastic_simple(r, deformation, sigma)
        close = xov[loc] > h0
        f_v = np.zeros_like(f_e)
        f_v[close] = force_viscous_simple(deformation, h0, r, u, -n[loc], mu_l, False)
        a[i] = np.sum(n[loc] * (f_e + f_v), axis=0) / m
    return a
