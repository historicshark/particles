import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import pandas as pd
from pathlib import Path
import itertools

from particle import Particle


def main():
    # --------------------------------------------
    # Inputs
    # --------------------------------------------

    # Fluid
    gx = 0
    gy = -9.81
    mu_l = 1
    rho_l = 0

    # Particles
    n_particles = 30

    diameter = 10
    diameter_stddev = 3
    diameter_min = 1

    velocity = 0
    velocity_stddev = 20

    rho_p = 1

    epsilon = .5

    drag = False
    gravity = False
    wall_collisions = True
    particle_collisions = True

    # Domain
    nx = 100
    ny = 100

    xmin = -100
    xmax = 100
    ymin = -100
    ymax = 100

    # Time
    dt = 0.001
    end_time = 150
    n_frames = 800
    save_animation = False

    # --------------------------------------------
    # Particles
    # --------------------------------------------
    # rp = np.full((n_particles,), 10)
    # x0 = np.array([50, -50])
    # y0 = np.array([50, -50])
    # u0 = np.array([-15, 15])
    # v0 = np.array([-15, 15])

    rng = np.random.default_rng()
    rp = rng.normal(diameter / 2, diameter_stddev / 2, (n_particles,))
    rp[rp < diameter_min / 2] = diameter_min / 2

    middle = 0.9
    x0 = middle * ((xmax - xmin) * rng.random((n_particles,)) - (xmax - xmin) / 2)
    y0 = middle * ((ymax - ymin) * rng.random((n_particles,)) - (ymax - ymin) / 2)

    u0 = rng.normal(velocity, velocity_stddev, (n_particles,))
    v0 = rng.normal(velocity, velocity_stddev, (n_particles,))

    rho = np.full((n_particles,), rho_p)
    # rho = np.array([rho_p, 10*rho_p])

    def distance(x1, y1, x2, y2):
        return np.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)

    # Check initially overlapping particles. Just delete one of them for now
    if n_particles > 1:
        particles_to_keep = np.full(rp.shape, True)
        n_particles_removed = 0
        for i in range(n_particles):
            for j in range(n_particles):
                if i != j:
                    if particles_to_keep[i] and particles_to_keep[j] and distance(x0[i], y0[i], x0[j], y0[j]) < rp[i] + rp[j]:
                        particles_to_keep[j] = False
                        n_particles_removed += 1
        n_particles -= n_particles_removed

        print(f'Removed {n_particles_removed} initially overlapping particles')

        rp = rp[particles_to_keep]
        rho = rho[particles_to_keep]
        x0 = x0[particles_to_keep]
        y0 = y0[particles_to_keep]
        u0 = u0[particles_to_keep]
        v0 = v0[particles_to_keep]
    else:
        x0 = np.zeros((n_particles,))
        y0 = np.zeros((n_particles,))

    position = np.stack((x0, y0), axis=1)
    velocity = np.stack((u0, v0), axis=1)
    v_l = np.zeros((2,))
    g = np.array([gx, gy])

    particles = []
    for i in range(n_particles):
        particles.append(
            Particle(i, rp[i], epsilon, rho[i], position[i], velocity[i], rho_l, mu_l, v_l, g, drag, gravity, wall_collisions, particle_collisions))

    particle_combinations = []
    if particle_collisions:
        for permutation in itertools.permutations(particles, 2):
            if (permutation[1], permutation[0]) not in particle_combinations:
                particle_combinations.append(permutation)

    # --------------------------------------------
    # Domain
    # --------------------------------------------
    walls = [xmax, xmin, ymax, ymin]

    x = np.linspace(xmin, xmax, nx)
    y = np.linspace(ymin, ymax, ny)

    X, Y = np.meshgrid(x, y)

    # --------------------------------------------
    # Run calculations
    # --------------------------------------------
    cwd = Path.cwd()

    domain_filepath = cwd.joinpath('domain.csv')
    particle_filepath = cwd.joinpath('particles.csv')
    time_filepath = cwd.joinpath('time.csv')
    position_filepath = cwd.joinpath('position.csv')
    velocity_filepath = cwd.joinpath('velocity.csv')
    results_filepath = cwd.joinpath('results.csv')

    print('Running calculation...')

    calculate(particles, particle_combinations, walls, dt, end_time, time_filepath, position_filepath, velocity_filepath, results_filepath)

    print('Finished calculations')

    # --------------------------------------------
    # Animation
    # --------------------------------------------
    if gravity:
        tracked_var = 'me'
    else:
        tracked_var = 'ke'

    animate(time_filepath, position_filepath, velocity_filepath, results_filepath, tracked_var, rp, xmin, xmax, ymin, ymax, n_frames, save_animation)


def calculate(particles, particle_combinations, walls, dt, end_time, time_filepath, position_filepath, velocity_filepath, results_filepath):
    time = np.arange(0, end_time, dt)

    # position_file = open(position_filepath, 'w')
    # velocity_file = open(velocity_filepath, 'w')
    # results_file = open(results_filepath, 'w')

    pd.DataFrame(time).to_csv(time_filepath, header=None, index=False)

    n_substeps = 10
    dt_substeps = dt / n_substeps

    with open(position_filepath, 'w') as position_file, open(velocity_filepath, 'w') as velocity_file, open(results_filepath, 'w') as results_file:
        results_file.write('ke,pe,me\n')
        for t in time:
            ke = 0.0
            pe = 0.0
            me = 0.0

            for particle in particles:
                particle.integrate(dt)

                # If particle hits a wall, reduce time step to integrate up to collision
                if particle.detect_wall_contact(walls):
                    i_substep = 0
                    keepgoing = True
                    while keepgoing:
                        particle.integrate(dt_substeps)
                        i_substep += 1
                        if particle.detect_wall_contact(walls) or i_substep > n_substeps:
                            keepgoing = False
                        else:
                            particle.update()

                    dt_remaining = dt - dt_substeps * i_substep
                    particle.integrate_update(dt_remaining)

            # particle collision detection
            for particle, other in particle_combinations:
                if not particle.finished_update and not other.finished_update and particle.distance(other) < particle.r + other.r:
                    i_substep = 0
                    keepgoing = True
                    while keepgoing:
                        particle.integrate(dt_substeps)
                        other.integrate(dt_substeps)
                        if particle.distance(other) < particle.r + other.r or i_substep > n_substeps:
                            keepgoing = False
                        else:
                            particle.update()
                            other.update()
                            i_substep += 1

                    dt_remaining = dt - dt_substeps * i_substep

                    particle.collision_acceleration(other, dt_remaining)
                    other.collision_acceleration(particle, dt_remaining)

                    particle.integrate_update(dt_remaining, True)
                    other.integrate_update(dt_remaining, True)

            for particle in particles:
                if not particle.finished_update:
                    particle.update()

                ke += particle.kinetic_energy()
                pe += particle.potential_energy(walls[3])
                me += particle.mechanical_energy(walls[3])

                position_file.write(particle.position_string() + ',')
                velocity_file.write(particle.velocity_string() + ',')

            position_file.write('\n')
            velocity_file.write('\n')
            results_file.write(f'{ke},{pe},{me}\n')

    return


def animate(time_file, position_file, velocity_file, results_file, tracked_var: str, rp, xmin, xmax, ymin, ymax, n_frames, save) -> None:
    cwd = Path.cwd()

    time = pd.read_csv(time_file, header=None)
    position = pd.read_csv(position_file, header=None)
    # velocity = pd.read_csv(velocity_file, header=None)
    results = pd.read_csv(results_file)

    fig, (ax, ax_bar) = plt.subplots(2, 1, height_ratios=[10, 1])
    wall_x = (xmin, xmax, xmax, xmin, xmin)
    wall_y = (ymin, ymin, ymax, ymax, ymin)
    walls = ax.plot(wall_x, wall_y, '')

    theta = np.linspace(0, 2 * np.pi)
    cos_theta = np.cos(theta)
    sin_theta = np.sin(theta)

    particles = []
    particle_xy_plot = lambda i, r, time_step: (r * cos_theta + position[2 * i][time_step], r * sin_theta + position[2 * i + 1][time_step])
    for i, r in enumerate(rp):
        x_particle, y_particle = particle_xy_plot(i, r, 0)
        particles.append(ax.plot(x_particle, y_particle, '-k')[0])
    ax.axis('equal')
    ax.set(xlim=(xmin, xmax), ylim=(xmin, ymax))

    bar, = ax_bar.barh(tracked_var.upper(), results[tracked_var][0])
    ax_bar.set_xlim(results[tracked_var].min(), results[tracked_var].max())

    def update(frame):
        for i, (r, particle) in enumerate(zip(rp, particles)):
            x_particle, y_particle = particle_xy_plot(i, r, frame)
            particle.set_xdata(x_particle)
            particle.set_ydata(y_particle)
        bar.set_width(results[tracked_var][frame])
        return particles

    n_time = len(time)
    frames = np.linspace(0, n_time - 1, n_frames).astype(int)
    ani = animation.FuncAnimation(fig=fig, func=update, frames=frames, interval=1)

    fig, ax = plt.subplots()
    ax.plot(time, results[tracked_var])
    ax.set_xlabel('t')
    ax.set_ylabel(tracked_var.upper())
    ax.grid()

    if save:
        print('Saving animation...')
        ani.save(filename=cwd.joinpath('animation.gif'), writer="pillow")
        fig.savefig('tracked_var.png', dpi=400)
        print(f"Animation saved to '{cwd.joinpath('animation.gif')}'")
    else:
        plt.show()

    return


if __name__ == '__main__':
    main()
