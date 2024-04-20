import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import pandas as pd
from pathlib import Path
import itertools
import argparse

from cell import Cell
from particle import Particle

# --------------------------------------------
# Global options
# --------------------------------------------
drag = True
gravity = False
wall_collisions = True
particle_collisions = False

dt = 0.0001
end_time = 1
n_frames = 100
save_animation = False

cwd = Path.cwd()
particle_filepath = cwd.joinpath('particles.csv')
domain_filepath = cwd.joinpath('domain.csv')
con_filepath = cwd.joinpath('con.csv')
time_filepath = cwd.joinpath('time.csv')
position_filepath = cwd.joinpath('position.csv')
velocity_filepath = cwd.joinpath('velocity.csv')
results_filepath = cwd.joinpath('results.csv')

def initialize():
    # --------------------------------------------
    # Inputs
    # --------------------------------------------
    # Fluid
    gx = 0
    gy = -9.81
    mu_l = 1e-3
    rho_l = 1000

    # Particles
    n_particles = 1

    diameter = .00004
    diameter_stddev = .00001
    diameter_min = .00002

    velocity = 0
    velocity_stddev = 0

    rho_p = 1.2

    epsilon = 1

    # Domain
    nx = 30
    ny = 30

    xmin = -.002
    xmax = .002
    ymin = -.002
    ymax = .002

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
                    if particles_to_keep[i] and particles_to_keep[j] and distance(x0[i], y0[i], x0[j], y0[j]) < rp[i] + \
                            rp[j]:
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

    v_l = np.zeros((2,))
    g = np.array([gx, gy])

    with open(particle_filepath, 'w') as f:
        for i in range(n_particles):
            f.write(f'{i},{rp[i]},{epsilon},{rho[i]},{x0[i]},{y0[i]},{u0[i]},{v0[i]},{rho_l},{mu_l},{v_l[0]},{v_l[1]},{g[0]},{g[1]},{xmax},{xmin},{ymax},{ymin}\n')

    # --------------------------------------------
    # Domain
    # --------------------------------------------
    walls = [xmax, xmin, ymax, ymin]

    x = np.linspace(xmin, xmax, nx)
    y = np.linspace(ymin, ymax, ny)

    X_, Y_ = np.meshgrid(x, y)

    X = X_.flatten()
    Y = Y_.flatten()

    coor = np.stack((X, Y), axis=1)
    con = np.zeros(((nx - 1) * (ny - 1), 4), dtype=int)
    for i in range(ny - 1):
        for j in range(nx - 1):
            ielem = (nx - 1) * i + j
            idx1 = nx * i + j
            idx2 = nx * i + 1 + j
            idx3 = nx * (i + 1) + 1 + j
            idx4 = nx * (i + 1) + j
            con[ielem] = [idx1, idx2, idx3, idx4]

    # Free vortex
    # u_theta = Gamma / (2 pi r)
    # scale = .00002
    # r = np.sqrt(X**2 + Y**2)
    # theta = np.arctan2(Y,X)
    # u_theta = scale / r
    # u_l = -u_theta * np.sin(theta)
    # v_l = u_theta * np.cos(theta)

    # uniform
    u_l = np.full(X.shape, (xmax-xmin)/20.)
    v_l = np.zeros(X.shape)

    np.savetxt(con_filepath, con)
    np.savetxt(domain_filepath, np.concatenate((coor, u_l[:, np.newaxis], v_l[:, np.newaxis]), axis=1))

    # fig, ax = plt.subplots()
    # c = ax.contourf(X_, Y_, scale / np.sqrt(X_**2 + Y_**2))
    # ax.axis('equal')
    # fig.colorbar(c)
    # plt.show()
    return


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--from_file', action='store_true', help='Run from particle file')
    args = parser.parse_args()
    from_file = args.from_file

    if not from_file:
        initialize()

    # --------------------------------------------
    # Run calculations
    # --------------------------------------------

    print('Running calculation...')

    calculate(dt, end_time)

    print('Finished calculations')

    # --------------------------------------------
    # Animation
    # --------------------------------------------
    if gravity:
        tracked_var = 'me'
    else:
        tracked_var = 'ke'

    animate(tracked_var, n_frames, save_animation)


def calculate(dt, end_time):
    time = np.arange(0, end_time, dt)
    time_progress = np.linspace(0,end_time, 100)
    i_progress = 0

    pd.DataFrame(time).to_csv(time_filepath, header=None, index=False)

    particles = []
    with open(particle_filepath) as f:
        for line in f:
            particles.append(Particle.from_csv(line, drag, gravity, wall_collisions, particle_collisions))

    particle_combinations = []
    if particles[0].collisions_on:
        for permutation in itertools.permutations(particles, 2):
            if (permutation[1], permutation[0]) not in particle_combinations:
                particle_combinations.append(permutation)

    con = np.genfromtxt(con_filepath, dtype=int)
    domain = np.genfromtxt(domain_filepath)
    coor = domain[:, :2]
    u_l = domain[:, 2:]

    cells = []
    for connectivity in con:
        cells.append(Cell(coor[connectivity], u_l[connectivity]))

    n_substeps = 10
    dt_substeps = dt / n_substeps

    with open(position_filepath, 'w') as position_file, open(velocity_filepath, 'w') as velocity_file, open(results_filepath, 'w') as results_file:
        results_file.write('ke,pe,me\n')
        for t in time:
            if t > time_progress[i_progress]:
                i_progress += 1
                print(f'{t = }')
            ke = 0.0
            pe = 0.0
            me = 0.0

            for particle in particles:
                for cell in itertools.takewhile(lambda cell_: not particle.in_cell_and_interpolate(cell_), cells):
                    pass

                particle.integrate(dt)

                # If particle hits a wall, reduce time step to integrate up to collision
                if particle.detect_wall_contact():
                    i_substep = 0
                    keepgoing = True
                    while keepgoing:
                        particle.integrate(dt_substeps)
                        i_substep += 1
                        if particle.detect_wall_contact() or i_substep > n_substeps:
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
                pe += particle.potential_energy()
                me += particle.mechanical_energy()

                position_file.write(particle.position_string() + ',')
                velocity_file.write(particle.velocity_string() + ',')

            position_file.write('\n')
            velocity_file.write('\n')
            results_file.write(f'{ke},{pe},{me}\n')

    return


def animate(tracked_var: str, n_frames, save) -> None:
    particles = pd.read_csv(particle_filepath, header=None)
    time = pd.read_csv(time_filepath, header=None)
    position = pd.read_csv(position_filepath, header=None)
    # velocity = pd.read_csv(velocity_filepath, header=None)
    results = pd.read_csv(results_filepath)

    rp = particles[1]
    xmax = particles[14][0]
    xmin = particles[15][0]
    ymax = particles[16][0]
    ymin = particles[17][0]

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
