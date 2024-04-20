import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import pandas as pd
from pathlib import Path
import itertools
import argparse

drag = False
gravity = True
wall_collisions = True
particle_collisions = False

def main():

    dt = 0.01
    end_time = 100
    n_frames = 200
    save_animation = False

    g = np.array([0, -9.81])
    mu_l = 1e-3
    rho_l = 0

    # Particles
    n_particles = 2

    rho_p = 1.2

    diameter = 10
    diameter_stddev = 0
    diameter_min = 1

    velocity = 20
    velocity_stddev = 5

    particle_density = 1.2

    epsilon = .8

    # Domain
    nx = 20
    ny = 20

    # xmin = -.002
    # xmax = .002
    # ymin = -.002
    # ymax = .002

    xmin = -100
    xmax = 100
    ymin = -100
    ymax = 100
    walls = np.array([xmin, xmax, ymin, ymax])

    flow_type = 'uniform'
    parameters = [0]

    # Initialization
    radius, rho_p = particle_properties(n_particles, diameter, diameter_stddev, diameter_min, rho_p)
    x, u = initial_conditions(n_particles, xmin, xmax, ymin, ymax, velocity, velocity_stddev)

    # remove overlapping particles
    n_particles, particles_to_keep = get_overlapping_index(n_particles, radius, x)

    radius = radius[particles_to_keep]
    rho_p = rho_p[particles_to_keep]
    x = x[particles_to_keep]
    u = u[particles_to_keep]

    # Grid
    coordinates, connectivity, background_flow = calculate_background_flow(xmin, xmax, ymin, ymax, nx, ny, flow_type, parameters)

    np.savetxt('radius.csv', radius)
    np.savetxt('rho_p.csv', rho_p)
    np.savetxt('x.csv', x)
    np.savetxt('u.csv', u)
    np.savetxt('coor.csv', coordinates)
    np.savetxt('con.csv', connectivity)
    np.savetxt('background_flow.csv', background_flow)

    calculate(dt, end_time, n_frames, n_particles, rho_l, mu_l, epsilon, g, walls, flow_type, parameters)

    animate(dt, end_time, n_frames, 'ke', xmin, xmax, ymin, ymax)

    return


def calculate(dt, end_time, n_frames, n_particles, rho_l, mu_l, epsilon, g, walls, flow_type, parameters):
    radius = np.genfromtxt('radius.csv')
    rho_p = np.genfromtxt('rho_p.csv')
    x = np.genfromtxt('x.csv')
    x_n = x.copy()
    u = np.genfromtxt('u.csv')
    u_n = u.copy()
    coordinates = np.genfromtxt('coor.csv')
    connectivity = np.genfromtxt('con.csv', dtype=int)
    background_flow = np.genfromtxt('background_flow.csv')

    m = mass(radius, rho_p)

    t = np.arange(0, end_time, dt)
    n_t = len(t)
    frames = np.linspace(0, n_t - 1, n_frames).astype(int)
    t_save = t[frames]

    n_substeps = 10
    dt_substeps = dt / n_substeps
    
    position = []
    results = []

    position_file = open('position.csv', 'w')
    results_file = open('results.csv', 'w')
    radius_file = open('radius.csv', 'w')
    for t_step in t:
        # if np.any(np.abs(x) > 95):
        #     print('wall')

        # flow_velocity = interpolate_flow_properties(x, coordinates, flow_velocity, connectivity)
        flow_velocity = interpolate_flow_properties_fast(x, flow_type, parameters)

        # x, u = integrate_euler(dt, x_n, u_n, flow_velocity, radius, rho_p, m, epsilon, rho_l, mu_l, g, walls)
        x, u = integrate_rk4(dt, x_n, u_n, flow_velocity, radius, rho_p, m, epsilon, rho_l, mu_l, g, walls)
        x_n = x
        u_n = u

        ke = np.sum(kinetic_energy(u, m))
        pe = np.sum(kinetic_energy(u, m))

        if t_step in t_save:
            print(f'time = {t_step}')

            for x_ in x.flatten():
                position_file.write(f'{x_},')
            position_file.write('\n')

            for r_ in radius:
                radius_file.write(f'{r_},')
            radius_file.write('\n')

            results_file.write(f'{ke},{pe}\n')

    np.savetxt('radius.csv', [radius])
    np.savetxt('t.csv', t_save)
    position_file.close()
    results_file.close()
    radius_file.close()
    return


def magnitude(a):
    return np.sqrt(np.sum(a*a, axis=1))


def area_xs(radius):
    return np.pi * radius**2


def volume(radius):
    return 4.0 / 3.0 * np.pi * radius**3


def mass(radius, rho_p):
    return rho_p * volume(radius)


def reynolds_number(u_rel, radius, rho_l, mu_l):
    return 2 * rho_l * magnitude(u_rel) * radius / mu_l


def drag_coefficient(re):
    cd = 24.0 / re * (2.0 / 3.0 + 1.0 / (12.0 / re + 3.0 / 4.0 * (1 + 3.315 / np.sqrt(re))))
    cd[re < 1e-8] = 0.0
    return cd


def drag_acceleration(particle_velocity, flow_velocity, radius, rho_l, mu_l, m):
    u_rel = particle_velocity - flow_velocity
    re = reynolds_number(u_rel, radius, rho_l, mu_l)
    cd = drag_coefficient(re)
    return 0.5 * rho_l * magnitude(u_rel) * u_rel * area_xs(radius) * cd / m


def buoyancy_acceleration(rho_p, rho_l, g):
    return (rho_p - rho_l) / rho_p * g


def detect_wall_contact(x, radius, walls):
    collision_xmin = radius > x[:,0] - walls[0]
    collision_xmax = radius > walls[1] - x[:,0]
    collision_ymin = radius > x[:,1] - walls[2]
    collision_ymax = radius > walls[3] - x[:,1]
    return np.stack((np.logical_or(collision_xmin, collision_xmax), np.logical_or(collision_ymin, collision_ymax)), axis=1)


def wall_collision_acceleration(wall_collision, particle_velocity, epsilon, dt):
    wall_collision_loc = np.any(wall_collision, axis=1)
    # a = np.zeros(particle_velocity.shape)
    a = -(1 + epsilon) * particle_velocity[wall_collision] / dt
    return a


def in_cell(x, vertices):
    vectors = vertices - x
    test = np.sign(np.cross(vectors, vectors[np.r_[1:4, 0]]))
    return np.all(test > 0) or np.all(test < 0)


def interpolate(x, vertices, velocity):
    A = np.array([[vertices[0, 0] * vertices[0, 1], vertices[0, 0], vertices[0, 1], 1],
                  [vertices[1, 0] * vertices[1, 1], vertices[1, 0], vertices[1, 1], 1],
                  [vertices[2, 0] * vertices[2, 1], vertices[2, 0], vertices[2, 1], 1],
                  [vertices[3, 0] * vertices[3, 1], vertices[3, 0], vertices[3, 1], 1],
                  ])
    coeff = np.linalg.solve(A, velocity)
    return coeff[0] * x[0] * x[1] + coeff[1] * x[0] + coeff[2] * x[1] + coeff[3]


def interpolate_flow_properties(x, coordinates,  flow_velocity, connectivity):
    u_l = np.zeros(x.shape)
    for i, particle in enumerate(x):
        for cell in connectivity:
            vertices = coordinates[cell]
            if in_cell(particle, vertices):
                u_l[i] = interpolate(particle, vertices, flow_velocity[cell])
                break
    return u_l


def interpolate_flow_properties_fast(x, flow_type, parameters):
    if flow_type == 'uniform':
        u_l = np.full(x[:,0].shape, parameters[0])
        v_l = np.zeros(x[:,1].shape)
    else:
        u_l = np.zeros(x[:,0].shape)
        v_l = np.zeros(x[:,1].shape)

    return np.stack((u_l, v_l), axis=1)


def apply_accelerations(particle_position, particle_velocity, flow_velocity, radius, rho_p, m, epsilon, rho_l, mu_l, g, wall_collision, dt):
    a = np.zeros(particle_velocity.shape)
    if drag:
        a += drag_acceleration(particle_velocity, flow_velocity, radius, rho_l, mu_l, m)
    if gravity:
        a += buoyancy_acceleration(rho_p, rho_l, g)
    if wall_collisions:
        a[wall_collision] = wall_collision_acceleration(wall_collision, particle_velocity, epsilon, dt)
    return a


def integrate_euler(dt, x_n, u_n, flow_velocity, radius, rho_p, m, epsilon, rho_l, mu_l, g, walls):
    wall_collision = detect_wall_contact(x_n, radius, walls)
    a1 = apply_accelerations(x_n, u_n, flow_velocity, radius, rho_p, m, epsilon, rho_l, mu_l, g, wall_collision, dt)
    u = u_n + a1 * dt
    x = 0.5 * a1 * dt**2 + u * dt + x_n
    return x, u


def integrate_rk4(dt, x_n, u_n, flow_velocity, radius, rho_p, m, epsilon, rho_l, mu_l, g, walls):
    wall_collision = detect_wall_contact(x_n, radius, walls)

    if not np.any(wall_collision):
        k1u = apply_accelerations(x_n,           u_n,           flow_velocity, radius, rho_p, m, epsilon, rho_l, mu_l, g, wall_collision, dt) * dt
        k1x = u_n * dt

        k2u = apply_accelerations(x_n + k1x / 2, u_n + k1u / 2, flow_velocity, radius, rho_p, m, epsilon, rho_l, mu_l, g, wall_collision, dt) * dt
        k2x = (u_n + k1u / 2) * dt

        k3u = apply_accelerations(x_n + k2x / 2, u_n + k2u / 2, flow_velocity, radius, rho_p, m, epsilon, rho_l, mu_l, g, wall_collision, dt) * dt
        k3x = (u_n + k2u / 2) * dt

        k4u = apply_accelerations(x_n + k3x,     u_n + k3u,     flow_velocity, radius, rho_p, m, epsilon, rho_l, mu_l, g, wall_collision, dt) * dt
        k4x = (u_n + k3u) * dt

        u = u_n + (k1u + 2 * (k2u + k3u) + k4u) / 6
        x = x_n + (k1x + 2 * (k2x + k3x) + k4x) / 6

    else:
        x, u = integrate_euler(dt, x_n, u_n, flow_velocity, radius, rho_p, m, epsilon, rho_l, mu_l, g, walls)

    return x, u


def kinetic_energy(u, m):
    return 0.5 * m * magnitude(u)**2


def particle_properties(n_particles, diameter, diameter_stddev, diameter_min, rho_p):
    rng = np.random.default_rng()
    rp = rng.normal(diameter / 2, diameter_stddev / 2, (n_particles,))
    rp[rp < diameter_min / 2] = diameter_min / 2

    rho_p = np.full((n_particles,), rho_p)

    return rp, rho_p


def initial_conditions(n_particles, xmin, xmax, ymin, ymax, velocity, velocity_stddev):
    rng = np.random.default_rng()
    middle = 0.9
    x0 = middle * ((xmax - xmin) * rng.random((n_particles,)) - (xmax - xmin) / 2)
    y0 = middle * ((ymax - ymin) * rng.random((n_particles,)) - (ymax - ymin) / 2)

    u0 = rng.normal(velocity, velocity_stddev, (n_particles,))
    v0 = rng.normal(velocity, velocity_stddev, (n_particles,))

    return np.stack((x0,y0), axis=1), np.stack((u0,v0), axis=1)


def get_overlapping_index(n_particles, radius, x):
    particles_to_keep = np.full((n_particles,), True)
    if n_particles > 1:
        for i in range(n_particles):
            for j in range(n_particles):
                if i != j:
                    if particles_to_keep[i] and particles_to_keep[j] and np.sqrt((x[i,0] - x[j,0])**2 + (x[i,1] - x[j,1])**2) < radius[i] + radius[j]:
                        particles_to_keep[j] = False
                        n_particles -= 1
    return n_particles, particles_to_keep


def calculate_background_flow(xmin, xmax, ymin, ymax, nx, ny, flow_type, parameters):
    x_ = np.linspace(xmin, xmax, nx)
    y_ = np.linspace(ymin, ymax, ny)

    x, y = np.meshgrid(x_, y_)

    x = x.flatten()
    y = y.flatten()

    coor = np.stack((x, y), axis=1)
    con = np.zeros(((nx - 1) * (ny - 1), 4), dtype=int)
    for i in range(ny - 1):
        for j in range(nx - 1):
            ielem = (nx - 1) * i + j
            idx1 = nx * i + j
            idx2 = nx * i + 1 + j
            idx3 = nx * (i + 1) + 1 + j
            idx4 = nx * (i + 1) + j
            con[ielem] = [idx1, idx2, idx3, idx4]

    if flow_type == 'uniform':
        u_l = np.full(x.shape, parameters[0])
        v_l = np.zeros(x.shape)
    else:
        u_l = np.zeros(x.shape)
        v_l = np.zeros(x.shape)

    return coor, con, np.stack((u_l, v_l), axis=1)


def animate(dt, end_time, n_frames, tracked_var: str, xmin, xmax, ymin, ymax, save=False) -> None:
    time = np.loadtxt('t.csv')
    position = np.genfromtxt('position.csv', delimiter=',')[:, :-1]
    results = np.genfromtxt('results.csv', delimiter=',')
    rp = np.genfromtxt('radius.csv', delimiter=',')[:, :-1]

    fig, (ax, ax_bar) = plt.subplots(2, 1, height_ratios=[10, 1])
    wall_x = (xmin, xmax, xmax, xmin, xmin)
    wall_y = (ymin, ymin, ymax, ymax, ymin)
    walls = ax.plot(wall_x, wall_y, '')

    theta = np.linspace(0, 2 * np.pi)
    cos_theta = np.cos(theta)
    sin_theta = np.sin(theta)

    def particle_xy_plot(i_particle, r, time_step):
        return r * cos_theta + position[time_step, 2 * i_particle], r * sin_theta + position[time_step, 2 * i_particle + 1]

    particles = []
    for i, radius in enumerate(rp[0]):
        x_particle, y_particle = particle_xy_plot(i, radius, 0)
        particles.append(ax.plot(x_particle, y_particle, '-k')[0])
    ax.axis('equal')
    ax.set(xlim=(xmin, xmax), ylim=(xmin, ymax))

    tracked_var_i = 0
    if tracked_var == 'ke':
        tracked_var_i = 0

    bar, = ax_bar.barh(tracked_var.upper(), results[0, tracked_var_i])
    ax_bar.set_xlim(results[:, tracked_var_i].min(), results[:, tracked_var_i].max())

    def update(frame):
        for i_particle, (r, particle) in enumerate(zip(rp[frame], particles)):
            x_particle, y_particle = particle_xy_plot(i_particle, r, frame)
            particle.set_xdata(x_particle)
            particle.set_ydata(y_particle)
        bar.set_width(results[frame, tracked_var_i])
        return particles

    n_time = len(time)
    frames = np.linspace(0, n_time - 1, n_frames).astype(int)
    ani = animation.FuncAnimation(fig=fig, func=update, frames=frames, interval=1)

    fig, ax = plt.subplots()
    ax.plot(time, results[:, tracked_var_i])
    ax.set_xlabel('t')
    ax.set_ylabel(tracked_var.upper())
    ax.grid()

    if save:
        print('Saving animation...')
        ani.save(filename='animation.gif', writer="pillow")
        fig.savefig('tracked_var.png', dpi=400)
        print(f"Animation saved")
    else:
        plt.show()

    return

if __name__ == '__main__':
    main()
