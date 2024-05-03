import numpy as np
import matplotlib.pyplot as plt
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('test_case', help='Test case to run', choices=['test','uniform','uniform-vertical','vortex','couette','poiseuille','breakup','terminal','custom'])
    args = parser.parse_args()
    test_case = args.test_case
    setup(test_case)
    return


def setup(test_case):
    dt = 1e-5
    end_time = .015
    n_frames = 200
    n_particles = 10
    
    drag = True
    particle_collisions = True
    gravity = False
    wall_collisions = True
    lift = True
    breakup = False

    g = np.array([0, -9.81])
    mu_l = 1e-3
    rho_l = 1000
    sigma = .072
    
    rho_p = 1.2
    mu_p = 1e-5

    diameter = 5e-3
    diameter_stddev = 1e-3
    diameter_min = 50e-6

    velocity = 0
    velocity_stddev = 0
    
    # Domain
    nx = 20
    ny = 20
        
    match test_case:
        case 'test':
            xmin = -.0025
            xmax = .0025
            ymin = -.0025
            ymax = .0025
            n_particles = 2
            flow_type = 'test'
            parameters = [.08]
        case 'uniform':
            xmin = -.0025
            xmax = .0025
            ymin = -.0025
            ymax = .0025
            flow_type = 'uniform'
            parameters = [.02]
        case 'uniform-vertical':
            gravity = True
            end_time = 1.5
            xmin = -.08
            xmax = .08
            ymin = 0
            ymax = 1
            flow_type = 'uniform-vertical'
            parameters = [0]
        case 'vortex':
            xmin = -.005
            xmax = .005
            ymin = -.005
            ymax = .005
            flow_type = 'vortex'
            parameters = [.002, .005]
            # parameters = [core-radius, max-tangential-velocity]
        case 'couette':
            gravity = True
            xmin = -.08
            xmax = .08
            ymin = 0
            ymax = 1
            flow_type = 'couette'
            parameters = [0.01, 0, mu_l, xmin, xmax]
            # u = U y/b + 1/(2 mu) (dpdx) (y^2 - by)
            # parameters = [U, dpdx, mu, ymin, ymax]
        case 'poiseuille':
            gravity = True
            dt = 1e-5
            end_time = 5
            xmin = -.08
            xmax = .08
            ymin = 0
            ymax = 1
            diameter = 10e-5
            diameter_stddev = 0
            flow_type = 'poiseuille'
            parameters = [0.1, (xmax+xmin)/2, xmax-xmin]
            # parameters = [U_max, center, R]
        case 'breakup':
            breakup = True
            gravity = True
            n_particles = 1
            dt = 1e-7
            end_time = 0.0004
            xmin = -.005
            xmax = .005
            ymin = 0
            ymax = .02
            diameter = 0.005
            flow_type = 'uniform-vertical'
            parameters = [0.1]
        case 'terminal':
            gravity = True
            wall_collisions = False
            diameter = 200e-6
            diameter_stddev = 0
            n_particles = 1
            flow_type = 'uniform-vertical'
            parameters = [0]
            xmin = -0.01
            xmax = 0.01
            ymin = 0
            ymax = 5
            diameter_stddev = 0
    
    walls = np.array([xmin, xmax, ymin, ymax])
    
    fig, ax = plt.subplots()
    ax.set_xlabel('X [m]')
    ax.set_ylabel('Y [m]')
    ax.axis('equal')
    
    # Initialization
    radius, rho_p = particle_properties(n_particles, diameter, diameter_stddev, diameter_min, rho_p)
    x, u = initial_conditions(n_particles, xmin, xmax, ymin, ymax, velocity, velocity_stddev, flow_type, parameters, ax, test_case)
    coordinates, connectivity, background_flow = calculate_background_flow(xmin, xmax, ymin, ymax, nx, ny, flow_type, parameters, fig, ax)
    
    np.savetxt('r.csv', radius, delimiter=',')
    np.savetxt('x.csv', x, delimiter=',')
    np.savetxt('u.csv', u, delimiter=',')
    np.savetxt('coor.csv', coordinates, delimiter=',')
    np.savetxt('con.csv', connectivity, delimiter=',')
    np.savetxt('background_flow.csv', background_flow, delimiter=',')
    
    fig.savefig('initial.png', dpi=300)

    with open('options.txt', 'w') as f:
        f.write(f'''{dt}
{end_time}
{n_frames}
{g[0]}
{g[1]}
{mu_l}
{rho_l}
{sigma}
{mu_p}
{rho_p}
{xmin}
{xmax}
{ymin}
{ymax}
{int(drag)}
{int(gravity)}
{int(wall_collisions)}
{int(particle_collisions)}
{int(lift)}
{int(breakup)}
{flow_type}
{','.join(str(p) for p in parameters)}''')
    plt.show()
    return


def particle_properties(n_particles, diameter, diameter_stddev, diameter_min, rho_p):
    rng = np.random.default_rng()
    rp = rng.normal(diameter / 2, diameter_stddev / 2, (n_particles,))
    rp[rp < diameter_min / 2] = diameter_min / 2

    return rp, rho_p


def initial_conditions(n_particles, xmin, xmax, ymin, ymax, velocity, velocity_stddev, flow_type, parameters, ax, test_case):
    rng = np.random.default_rng()
    middle = 0.9
    x0 = (xmax - xmin) * ((2 * middle - 1) * rng.random((n_particles,)) + 1 - middle) + xmin
    y0 = (ymax - ymin) * ((2 * middle - 1) * rng.random((n_particles,)) + 1 - middle) + ymin
    # x0 = np.array((-20,20))
    # y0 = np.zeros((n_particles,))
    
    if flow_type == 'test':
        midx = (xmax - xmin) / 2
        midy = (ymax + ymin) / 2
        x0 = np.array((xmin/2,xmax/2))
        y0 = np.full_like(x0, midy)
        u0 = np.array(((parameters[0],0),(-parameters[0],0)))
    elif flow_type in ['couette', 'uniform-vertical', 'poiseuille']:
        y0 = (ymax/3 - ymin) * ((2 * middle - 1) * rng.random((n_particles,)) + 1 - middle) + ymin
    
    if test_case == 'breakup':
        x0 = np.zeros((n_particles,))
        y0 = np.full_like(x0, (ymax+ymin)/6)

    # u0 = rng.normal(velocity, velocity_stddev, (n_particles,2))
    # v0 = rng.normal(velocity, velocity_stddev, (n_particles,))

    # u0 = np.zeros((n_particles,))
    # v0 = np.zeros((n_particles,))
    
    ax.plot(x0,y0, 'ko')
    
    x0 = np.stack((x0,y0), axis=1)

    u0 = interpolate_flow_properties_fast(x0, flow_type, parameters)
    return x0, u0


def calculate_background_flow(xmin, xmax, ymin, ymax, nx, ny, flow_type, parameters, fig, ax):
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
    
    u_l = interpolate_flow_properties_fast(np.stack((x,y), axis=1), flow_type, parameters)
    
    c = ax.contourf(x_, y_, np.reshape(np.sqrt(np.sum(u_l * u_l, axis=1)), (len(y_), len(x_))), levels=41)
    cbar = fig.colorbar(c)
    cbar.ax.set_ylabel('Velocity Magnitude [m/s]')

    return coor, con, u_l


def interpolate_flow_properties_fast(x, flow_type, parameters):
    if flow_type == 'uniform':
        # parameters = [velocity-x]
        u_l = np.full(x[:,0].shape, parameters[0])
        v_l = np.zeros(x[:,1].shape)
    elif flow_type == 'uniform-vertical':
        # parameters = [velocity-y]
        v_l = np.full(x[:,0].shape, parameters[0])
        u_l = np.zeros_like(v_l)
    elif flow_type == 'vortex':
        # parameters = [core-diameter, max-tangential-velocity]
        r = np.sqrt(np.sum(x**2, axis=1))
        theta = np.arctan2(x[:,1], x[:,0])
        u_theta = np.zeros_like(r)
        core = r < parameters[0]
        u_theta[core] = r[core] * parameters[1] / parameters[0]
        u_theta[~core] = parameters[0] * parameters[1] / r[~core]
        u_l = -u_theta * np.sin(theta)
        v_l = u_theta * np.cos(theta)
    elif flow_type == 'test':
        u_l = np.zeros_like(x[:,0])
        u_l[x[:,0] < 0] = parameters[0]
        u_l[x[:,0] > 0] = -parameters[0]
        v_l = np.zeros_like(u_l)
    elif flow_type == 'couette':
        # parameters = [U, dpdx, mu, ymin, ymax]
        v_l = parameters[0] * (x[:,0] - parameters[3]) / (parameters[4] - parameters[3]) + 1. / (2. * parameters[2]) * parameters[1] * ((x[:,0] - parameters[3])**2 - (parameters[4] - parameters[3]) * (x[:,0] - parameters[3]))
        u_l = np.zeros_like(v_l)
    elif flow_type == 'poiseuille':
        v_l = parameters[0] * (1 - ((x[:,0] - parameters[1])/parameters[2])**2)
        u_l = np.zeros_like(v_l)
    else:
        u_l = np.zeros(x[:,0].shape)
        v_l = np.zeros(x[:,1].shape)

    return np.stack((u_l, v_l), axis=1)


if __name__ == '__main__':
    main()
