import numpy as np
import argparse

def main():
    drag = True
    particle_collisions = True
    gravity = False
    wall_collisions = True
#    particle_breakup_coalescence = False
    
    dt = 1e-5
    end_time = .5
    n_frames = 200

    g = np.array([0, -9.81])
    mu_l = 1e-3
    rho_l = 1000
    sigma = .072

    # Particles
    n_particles = 50

    rho_p = 1.2
    mu_p = 1e-5

    diameter = 200e-6
    diameter_stddev = 0
    diameter_min = 50e-6

    velocity = 0
    velocity_stddev = 0

    epsilon = .8

    # Domain
    nx = 20
    ny = 20

    xmin = -.0025
    xmax = .0025
    ymin = -.0025
    ymax = .0025
    walls = np.array([xmin, xmax, ymin, ymax])

#    flow_type = 'uniform'
#    parameters = [.02]

    flow_type = 'vortex'
    parameters = [.001, .005]

#    flow_type = 'test'
#    parameters = [.05]

    # Initialization
    radius, rho_p = particle_properties(n_particles, diameter, diameter_stddev, diameter_min, rho_p)
    x, u = initial_conditions(n_particles, xmin, xmax, ymin, ymax, velocity, velocity_stddev, flow_type, parameters)
    coordinates, connectivity, background_flow = calculate_background_flow(xmin, xmax, ymin, ymax, nx, ny, flow_type, parameters)
    
    np.savetxt('r.csv', radius, delimiter=',')
    np.savetxt('x.csv', x, delimiter=',')
    np.savetxt('u.csv', u, delimiter=',')
    np.savetxt('coor.csv', coordinates, delimiter=',')
    np.savetxt('con.csv', connectivity, delimiter=',')
    np.savetxt('background_flow.csv', background_flow, delimiter=',')

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
{flow_type}
{','.join(str(p) for p in parameters)}''')
    
    return


def particle_properties(n_particles, diameter, diameter_stddev, diameter_min, rho_p):
    rng = np.random.default_rng()
    rp = rng.normal(diameter / 2, diameter_stddev / 2, (n_particles,))
    rp[rp < diameter_min / 2] = diameter_min / 2

    return rp, rho_p


def initial_conditions(n_particles, xmin, xmax, ymin, ymax, velocity, velocity_stddev, flow_type, parameters):
    rng = np.random.default_rng()
    middle = 0.9
    x0 = (xmax - xmin) * ((2 * middle - 1) * rng.random((n_particles,)) + 1 - middle) + xmin
    y0 = (ymax - ymin) * ((2 * middle - 1) * rng.random((n_particles,)) + 1 - middle) + ymin
    # x0 = np.array((-20,20))
    # y0 = np.zeros((n_particles,))

    x0 = np.stack((x0,y0), axis=1)

    u0 = interpolate_flow_properties_fast(x0, flow_type, parameters)

    if flow_type == 'test':
        midx = (xmax - xmin) / 2
        midy = (ymax + ymin) / 2
        x0 = np.array((xmin/2,xmax/2))
        y0 = np.full_like(x0, midy)
        x0 = np.stack((x0, y0), axis=1)
        u0 = np.array(((parameters[0],0),(-parameters[0],0)))

    # u0 = rng.normal(velocity, velocity_stddev, (n_particles,2))
    # v0 = rng.normal(velocity, velocity_stddev, (n_particles,))

    # u0 = np.zeros((n_particles,))
    # v0 = np.zeros((n_particles,))

    return x0, u0


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

    u_l = interpolate_flow_properties_fast(np.stack((x,y), axis=1), flow_type, parameters)

    return coor, con, u_l


def interpolate_flow_properties_fast(x, flow_type, parameters):
    if flow_type == 'uniform':
        # parameters = [velocity-x]
        u_l = np.full(x[:,0].shape, parameters[0])
        v_l = np.zeros(x[:,1].shape)
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
        u_l = np.array([parameters[0], -parameters[0]])
        v_l = np.zeros_like(u_l)
    else:
        u_l = np.zeros(x[:,0].shape)
        v_l = np.zeros(x[:,1].shape)

    return np.stack((u_l, v_l), axis=1)


if __name__ == '__main__':
    main()
