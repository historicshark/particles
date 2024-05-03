import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import pandas as pd
import argparse
import csv

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-save', action='store_true', help='Save animation to file')
    args = parser.parse_args()
    
    save_animation = args.save
    tracked_var = ''
    animate(tracked_var, save_animation)
    return


def animate(tracked_var, save_animation):
    with open('options.txt') as f:
        for i in range(10):
            f.readline()
        xmin = float(f.readline())
        xmax = float(f.readline())
        ymin = float(f.readline())
        ymax = float(f.readline())
    
    wall_x = (xmin, xmax, xmax, xmin, xmin)
    wall_y = (ymin, ymin, ymax, ymax, ymin)
    
    theta = np.linspace(0, 2 * np.pi)
    cos_theta = np.cos(theta)
    sin_theta = np.sin(theta)
    
    fig, ax = plt.subplots(dpi=300)
    ax.set_ylim(ymin,ymax)
    ax.axis('equal')
    
    artists = []
    with open('position.csv') as x_file, open('radius.csv') as r_file, open('time.csv') as t_file:
        x_reader = csv.reader(x_file)
        r_reader = csv.reader(r_file)
        t_reader = csv.reader(t_file)
        
        for x_row_, r_row_ in zip(x_reader, r_reader):
            artists_step = []
            artists_step.append(ax.plot(wall_x, wall_y, 'k-')[0])
            x_row = [float(x) for x in x_row_[:-1]]
            r_row = [float(x) for x in r_row_[:-1]]
            # t_row = [float(x) for x in t_row_[:-1]]
            for x, y, r in zip(x_row[::2], x_row[1::2], r_row):
                if np.abs(x) < 10 and np.abs(y) < 10:
                    x_plot = r * cos_theta + x
                    y_plot = r * sin_theta + y
                    artists_step.append(ax.plot(x_plot, y_plot, 'k-')[0])
            artists.append(artists_step)
    
    ani = animation.ArtistAnimation(fig=fig, artists=artists, interval=1)
    
    if save_animation:
        print('Saving animation...')
        ani.save(filename='animation.gif', writer="pillow")
        if tracked_var:
            fig.savefig('tracked_var.png', dpi=400)
        print("Animation saved")
    else:
        plt.show()

    return


if __name__ == '__main__':
    main()
