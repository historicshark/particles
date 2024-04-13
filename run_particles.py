import csv
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import pandas as pd
from pathlib import Path
import platform
import subprocess

def main():
    # --------------------------------------------
    # Inputs
    # --------------------------------------------
    
    # Fluid
    gx = 0
    gy = -9.81
    mu_l = 1
    rho_l = 1
    
    # Particles
    n_particles = 3
    
    diameter = 5
    diameter_stddev = .5
    diameter_min = 1
    
    velocity = 0
    velocity_stddev = 15
    
    rho_p = 1
    
    epsilon = 1
    
    drag = False
    gravity = False
    walls = True
    
    # Domain
    nx = 100
    ny = 100
    
    xmin = -100
    xmax = 100
    ymin = -100
    ymax = 100
    
    # Time
    dt = 0.001
    end_time = 20
    n_frames = 800
    save_animation = False
    
    # --------------------------------------------
    # Particles
    # --------------------------------------------
    # rp = np.full((n_particles,), 10)
    # x0 = np.array([50, -50])
    # y0 = np.array([0, 0])
    # u0 = np.array([-15,15])
    # v0 = np.array([0,0])
    
    rng = np.random.default_rng()
    rp = rng.normal(diameter/2, diameter_stddev/2, (n_particles,))
    rp[rp < diameter_min/2] = diameter_min/2
    
    rho = np.full((n_particles,), rho_p)

    middle = 0.9
    x0 = middle * ((xmax - xmin) * rng.random((n_particles,)) - (xmax - xmin) / 2)
    y0 = middle * ((ymax - ymin) * rng.random((n_particles,)) - (ymax - ymin) / 2)

    u0 = rng.normal(velocity, velocity_stddev, (n_particles,))
    v0 = rng.normal(velocity, velocity_stddev, (n_particles,))
    
    def distance(x1, y1, x2, y2):
        return np.sqrt((x1 - x2)**2 + (y1 - y2)**2)
    
#   Check initially overlapping particles. Just delete one of them for now
    particles_to_keep = np.full(rp.shape, True)
    n_particles_removed = 0
    for i in range(n_particles):
        for j in range(n_particles):
            if i != j:
                if particles_to_keep[i] and particles_to_keep[j] and distance(x0[i], y0[i], x0[j], y0[j]) < rp[i] + rp[j]:
                    particles_to_keep[j] = False
                    print(f'{i = }, {j = }')
                    n_particles_removed += 1
    n_particles -= n_particles_removed
    
    print(f'Removed {n_particles_removed} initially overlapping particles')
    print(f'{particles_to_keep = }')
    
    rp = rp[particles_to_keep]
    rho = rho[particles_to_keep]
    x0 = x0[particles_to_keep]
    y0 = y0[particles_to_keep]
    u0 = u0[particles_to_keep]
    v0 = v0[particles_to_keep]
    
    # --------------------------------------------
    # Domain
    # --------------------------------------------
    x = np.linspace(xmin, xmax, nx)
    y = np.linspace(ymin, ymax, ny)
    
    X, Y = np.meshgrid(x, y)
    
    # --------------------------------------------
    # Run calculations
    # --------------------------------------------
    
    cwd = Path.cwd()
    
    domain_file = cwd.joinpath('domain.csv')
    particle_file = cwd.joinpath('particles.csv')
    time_file = cwd.joinpath('time.csv')
    position_file = cwd.joinpath('position.csv')
    velocity_file = cwd.joinpath('velocity.csv')
    results_file = cwd.joinpath('results.csv')
    
    time_file.unlink(missing_ok=True)
    position_file.unlink(missing_ok=True)
    velocity_file.unlink(missing_ok=True)
    results_file.unlink(missing_ok=True)
    
    command = []
    if platform.system() == 'Darwin':
        command.append('./main')
        executable = cwd.joinpath('main')
        executable.unlink(missing_ok=True)
        executable.symlink_to('/Users/mschmidt/Library/Developer/Xcode/DerivedData/particles_xcode-ewvadkdtekwkuzgyselujxvcctby/Build/Products/Debug/particles_xcode')
    else:
        command.append('main.exe')
    for x in [domain_file, particle_file, time_file, position_file, velocity_file, results_file]:
        command.append(x)
    
    parameters = [gx, gy, mu_l, rho_l, epsilon, xmin, xmax, ymin, ymax, dt, end_time, str(drag), str(gravity), str(walls)]
    
    for parameter in parameters: command.append(str(parameter))
    
    with open(domain_file, 'w') as f:
        writer = csv.writer(f)
        writer.writerow([gx, gy, mu_l, rho_l, epsilon, xmin, xmax, ymin, ymax, dt, end_time])
    
    pd.DataFrame({'r': rp, 'rho': rho, 'x0': x0, 'y0': y0, 'u0': u0, 'v0': v0}).to_csv(particle_file, header=None)
    
    print('Running calculation...')
    print(command)
    result = subprocess.run(command)
    print(f'\n{result = }\n')
    print('Finished calculations')
    
    # --------------------------------------------
    # Animation
    # --------------------------------------------
    animate(time_file, position_file, velocity_file, results_file, rp, xmin, xmax, ymin, ymax, n_frames, save_animation)
    
    return result

def animate(time_file, position_file, velocity_file, results_file, rp, xmin, xmax, ymin, ymax, n_frames, save):
    cwd = Path.cwd()
    
    time = pd.read_csv(time_file, header=None)
    position = pd.read_csv(position_file, header=None)
    # velocity = pd.read_csv(velocity_file, header=None)
    results = pd.read_csv(results_file)
    
    fig, (ax, ax_bar) = plt.subplots(2,1, height_ratios=[10,1])
    wall_x = (xmin,xmax,xmax,xmin,xmin)
    wall_y = (ymin,ymin,ymax,ymax,ymin)
    walls = ax.plot(wall_x, wall_y, '')
    
    theta = np.linspace(0, 2 * np.pi)
    cos_theta = np.cos(theta)
    sin_theta = np.sin(theta)
    
    particles = []
    particle_xy_plot = lambda i, r, time_step: (r * cos_theta + position[2*i][time_step], r * sin_theta + position[2*i+1][time_step])
    for i, r in enumerate(rp):
        x_particle, y_particle = particle_xy_plot(i, r, 0)
        particles.append(ax.plot(x_particle, y_particle, '-k')[0])
    ax.axis('equal')
    ax.set(xlim=(xmin, xmax), ylim=(xmin,ymax))
    
    bar, = ax_bar.barh('KE', results['ke'][0])
    ax_bar.set_xlim(results['ke'].min(), results['ke'].max())
    
    def update(frame):
        for i, (r, particle) in enumerate(zip(rp, particles)):
            x_particle, y_particle = particle_xy_plot(i, r, frame)
            particle.set_xdata(x_particle)
            particle.set_ydata(y_particle)
        bar.set_width(results['ke'][frame])
        return particles
    
    n_time = len(time)
    frames = np.linspace(0, n_time - 1, n_frames).astype(int)
    ani = animation.FuncAnimation(fig=fig, func=update, frames=frames, interval=1)
    
    if save:
        print('Saving animation...')
        ani.save(filename=cwd.joinpath('animation.gif'), writer="pillow")
        print(f"Animation saved to '{cwd.joinpath('animation.gif')}'")
    else:
        plt.show()

if __name__ == "__main__":
    main()
