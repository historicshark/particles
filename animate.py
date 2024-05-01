import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import pandas as pd
import argparse

def main():
    tracked_var = 'ke'
    save_animation = False
    
    with open('options.txt') as f:
        for i in range(10):
            f.readline()
        xmin = float(f.readline())
        xmax = float(f.readline())
        ymin = float(f.readline())
        ymax = float(f.readline())
    time = np.loadtxt('time.csv')
    position = pd.read_csv('position.csv').to_numpy()[:, :-1]
    results = np.genfromtxt('results.csv', delimiter=',')
    rp = pd.read_csv('radius.csv').to_numpy()[:, :-1]

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
        if np.isnan(radius):
            particles.append(ax.plot(x_particle, y_particle, '-w')[0])
        else:
            particles.append(ax.plot(x_particle, y_particle, '-k')[0])
    ax.axis('equal')
    ax.set(xlim=(xmin, xmax), ylim=(xmin, ymax))

    tracked_var_i = 0
    if tracked_var == 'ke':
        tracked_var_i = 0
    elif tracked_var == 'pe':
        tracked_var_i = 1
    elif tracked_var == 'me':
        tracked_var_i = 2

    bar, = ax_bar.barh(tracked_var.upper(), results[0, tracked_var_i])
    # ax_bar.set_xlim(results[:, tracked_var_i].min(), results[:, tracked_var_i].max())

    def update(frame):
        for i_particle, (r, particle) in enumerate(zip(rp[frame], particles)):
            x_particle, y_particle = particle_xy_plot(i_particle, r, frame)
            particle.set_xdata(x_particle)
            particle.set_ydata(y_particle)
        bar.set_width(results[frame, tracked_var_i])
        return particles

    n_time = len(time)
    ani = animation.FuncAnimation(fig=fig, func=update, frames=n_time-1, interval=1)

    fig, ax = plt.subplots()
    ax.plot(time, results[:, tracked_var_i])
    ax.set_xlabel('t')
    ax.set_ylabel(tracked_var.upper())
    ax.grid()

    if save_animation:
        print('Saving animation...')
        ani.save(filename='animation.gif', writer="pillow")
        fig.savefig('tracked_var.png', dpi=400)
        print(f"Animation saved")
    else:
        plt.show()

    return


if __name__ == '__main__':
    main()
