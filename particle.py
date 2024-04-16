import numpy as np


class Particle:
    def __init__(self, id_, r, epsilon, rho, x, u, rho_l, mu_l, v_l, g, drag, grav, wall):
        self.id = id_
        self.r = r
        self.epsilon = epsilon
        self.rho = rho
        self.x = x
        self.u = u
        self.x_n = x
        self.u_n = u
        self.rho_l = rho_l
        self.mu_l = mu_l
        self.v_l = v_l
        self.g = g
        self.drag = drag
        self.grav = grav
        self.wall = wall
        self.collided = False
        self.reynolds_number(self.u)

    def diameter(self):
        return 2 * self.r

    def volume(self):
        return 4. / 3. * np.pi * self.r ** 3

    def area_xs(self):
        return np.pi * self.r ** 2

    def mass(self):
        return self.rho * self.volume()

    def distance_vector(self, other):
        return other.x - self.x

    def distance(self, other):
        return np.sqrt(self.distance_vector(other) @ self.distance_vector(other))

    def reynolds_number(self, u_p):
        v_rel = u_p - self.v_l
        self.re = self.rho_l * np.sqrt(v_rel @ v_rel) * self.diameter() / self.mu_l
        return

    def kinetic_energy(self):
        return 0.5 * self.mass() * (self.u @ self.u)

    def drag_coefficient(self):
        if self.re > 1e-10:
            return 24.0 / self.re * (2.0 / 3.0 + (12.0 / self.re + 3.0 / 4.0 * (1 + 3.315 / np.sqrt(self.re)) ** (-1)))
        else:
            return 0.0

    def drag_acceleration(self, u_p):
        self.reynolds_number(u_p)
        v_rel = self.v_l - u_p
        return 1.0 / 2.0 * self.rho_l / self.rho * np.sqrt(v_rel @ v_rel) * v_rel * self.area_xs() / self.volume() * self.drag_coefficient()

    def buoyancy_acceleration(self):
        return (self.rho - self.rho_l) / self.rho * self.g

    def wall_contact_acceleration(self, dt, walls):
        a = np.zeros((2,))
        if self.r > walls[0] - self.x[0] or self.r > self.x[0] - walls[1]:
            a[0] += -(1 + self.epsilon) * self.u[0] / dt
        if self.r > walls[2] - self.x[1] or self.r > self.x[1] - walls[3]:
            a[1] += -(1 + self.epsilon) * self.u[1] / dt

        return a

    def collision_acceleration(self, other, dt):
        dr = self.distance_vector(other)
        n = dr / np.sqrt(dr @ dr)
        # dv = np.dot(other.u - self.u, n)
        # m_ratio = other.mass() / (other.mass() + self.mass())
        u_prime = self.u + n * np.dot(other.u - self.u, n) * (1 + self.epsilon) * other.mass() / (other.mass() + self.mass())
        self.a_collision = (u_prime - self.u) / dt
        return

    def apply_accelerations(self, u_p, dt, walls, collision: bool):
        a = np.zeros((2,))
        if self.drag: a += self.drag_acceleration(u_p)
        if self.grav: a += self.buoyancy_acceleration()
        if self.wall: a += self.wall_contact_acceleration(dt, walls)
        if collision: a += self.a_collision
        return a

    def integrate(self, dt, walls, collision=False):
        a1 = self.apply_accelerations(self.u_n, dt, walls, collision)

        self.u = self.u_n + dt * a1
        self.x = self.x_n + dt * self.u
        return

    def update(self):
        self.x_n = self.x
        self.u_n = self.u
        return

    def integrate_update(self, dt, walls, collision=False):
        self.integrate(dt, walls, collision)
        self.update()
        return

    def position_string(self):
        return f'{self.x[0]},{self.x[1]}'

    def velocity_string(self):
        return f'{self.u[0]},{self.u[1]}'
