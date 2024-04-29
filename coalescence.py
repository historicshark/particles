import numpy as np


def update_contact_time(contact_time, particle_collision, dt):
    contact_time[particle_collision] += dt
    contact_time[~particle_collision] = 0
    return contact_time


def drainage_time(radius, rho_l, sigma):
    r1, r2 = np.meshgrid(radius, radius)
    deq = 2 * (1 / (2 * r1) + 1 / (2 * r2))**(-1)
    return np.sqrt(deq**3 * rho_l / (128 * sigma)) * np.log(1e4)


def check_coalescence(contact_time, t_drainage):
    return contact_time > t_drainage


def merge_particles(particles_to_merge_matrix, radius, velocity):
    merged = []
    for i, (loc, r, u) in enumerate(zip(particles_to_merge_matrix,
                                        radius,
                                        velocity)):
        if np.any(loc) and i not in merged:
            r_other = radius[loc][0]
            u_other = velocity[loc][0]
            merged.append(np.nonzero(loc)[0])
            da3 = (2 * r)**3
            db3 = (2 * r_other)**3
            new_radius = (da3 + db3)**(1/3)
            new_velocity = (da3 * np.sqrt(np.sum(u * u)) + db3 * np.sqrt(np.sum(u_other * u_other))) / (da3 + db3) * u / np.sqrt(np.sum(u * u))
            radius[i] = new_radius
            velocity[i] = new_velocity
    return radius, velocity, merged
