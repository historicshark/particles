import numpy as np

from cell import Cell

class Particle:
    def __init__(self, id_, r, epsilon, rho, x, u, rho_l, mu_l, v_l, g, xmax, xmin, ymax, ymin, drag, grav, wall, collisions):
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
        self.drag_on = drag
        self.gravity_on = grav
        self.walls_on = wall
        self.collisions_on = collisions
        self.re = 0.0
        self.walls = np.array([xmax, xmin, ymax, ymin])
        self.wall_collision = np.full((4,), False)
        self.wall_contact = np.full((4,), False)
        self.a_collision = 0.0
        self.finished_update = False

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

    def potential_energy(self):
        return self.mass() * -self.g[1] * (self.x[1] - self.walls[3])

    def mechanical_energy(self):
        return self.kinetic_energy() + self.potential_energy()

    def drag_coefficient(self):
        if self.re > 1e-10:
            return 24.0 / self.re * (2.0 / 3.0 + (12.0 / self.re + 3.0 / 4.0 * (1 + 3.315 / np.sqrt(self.re)))**(-1))
        else:
            return 0.0

    def drag_acceleration(self, u_p):
        self.reynolds_number(u_p)
        v_rel = self.v_l - u_p
        return 0.5 * self.rho_l * np.sqrt(v_rel @ v_rel) * v_rel * self.area_xs() * self.drag_coefficient() / self.mass()

    def buoyancy_acceleration(self):
        return (self.rho - self.rho_l) / self.rho * self.g

    def detect_wall_contact(self):
        self.wall_collision[:] = False
        if self.walls_on:
            if self.r > self.walls[0] - self.x[0]:
                self.wall_collision[0] = True
            if self.r > self.x[0] - self.walls[1]:
                self.wall_collision[1] = True
            if self.r > self.walls[2] - self.x[1]:
                self.wall_collision[2] = True
            if self.r > self.x[1] - self.walls[3]:
                self.wall_collision[3] = True
        return np.any(self.wall_collision)

    def wall_contact_acceleration(self, dt):
        a = np.zeros((2,))
        if self.walls_on:
            if self.wall_collision[0] or self.wall_collision[1]:
                a[0] += -(1 + self.epsilon) * self.u[0] / dt
            if self.wall_collision[2] or self.wall_collision[3]:
                a[1] += -(1 + self.epsilon) * self.u[1] / dt
        return a

    def collision_acceleration(self, other, dt):
        dr = self.distance_vector(other)
        n = dr / np.sqrt(dr @ dr)
        # dv = np.dot(other.u - self.u, n)
        # m_ratio = other.mass() / (other.mass() + self.mass())
        u_prime = self.u + n * np.dot(other.u - self.u, n) * (1 + self.epsilon) * other.mass() / (
                other.mass() + self.mass())
        self.a_collision = (u_prime - self.u) / dt
        return

    def apply_accelerations(self, u_p, dt, collision: bool):
        a = np.zeros((2,))
        if self.drag_on:
            a += self.drag_acceleration(u_p)
        if self.gravity_on:
            a += self.buoyancy_acceleration()
        if self.walls_on:
            a += self.wall_contact_acceleration(dt)
        if collision:
            a += self.a_collision
        return a

    def integrate(self, dt, collision=False):
        a1 = self.apply_accelerations(self.u_n,               dt, collision)
        v1 = self.u_n
        a2 = self.apply_accelerations(self.u_n + dt * a1 / 2, dt, collision)
        v2 = self.u_n + dt * v1 / 2
        a3 = self.apply_accelerations(self.u_n + dt * a2 / 2, dt, collision)
        v3 = self.x_n
        a4 = self.apply_accelerations(self.u_n + dt * a3,     dt, collision)

        self.u = self.u_n + dt * a1
        # self.x = 0.5 *

        self.finished_update = False
        return

    def update(self):
        self.x_n = self.x
        self.u_n = self.u
        self.finished_update = True
        return

    def integrate_update(self, dt, collision=False):
        self.integrate(dt, collision)
        self.update()
        return

    def interpolate(self, cell):
        A = np.array([[cell.vertices[0, 0] * cell.vertices[0, 1], cell.vertices[0, 0], cell.vertices[0, 1], 1],
                      [cell.vertices[1, 0] * cell.vertices[1, 1], cell.vertices[1, 0], cell.vertices[1, 1], 1],
                      [cell.vertices[2, 0] * cell.vertices[2, 1], cell.vertices[2, 0], cell.vertices[2, 1], 1],
                      [cell.vertices[3, 0] * cell.vertices[3, 1], cell.vertices[3, 0], cell.vertices[3, 1], 1],
                      ])
        coeff = np.linalg.solve(A, cell.velocity)
        return coeff[0] * self.x[0] * self.x[1] + coeff[1] * self.x[0] + coeff[2] * self.x[1] + coeff[3]

    def in_cell(self, cell: Cell):
        vectors = cell.vertices - self.x
        test = np.sign(np.cross(vectors, vectors[np.r_[1:4,0]]))
        return np.all(test > 0) or np.all(test < 0)

    def in_cell_and_interpolate(self, cell: Cell):
        if self.in_cell(cell):
            self.v_l = self.interpolate(cell)
            return True
        return False

    def position_string(self):
        return f'{self.x[0]},{self.x[1]}'

    def velocity_string(self):
        return f'{self.u[0]},{self.u[1]}'

    @classmethod
    def from_csv(cls, line: str, drag, gravity, wall_collisions, particle_collisions):
        p = line[:-1].split(',')
        id_ = int(p[0])
        r, epsilon, rho, x0, y0, u0, v0, rho_l, mu_l, vlx, vly, gx, gy, xmax, xmin, ymax, ymin = map(float, p[1:18])

        particle = cls(id_, r, epsilon, rho, np.array([x0, y0]), np.array([u0, v0]), rho_l, mu_l, np.array([vlx, vly]), np.array([gx, gy]), xmax, xmin, ymax, ymin, drag, gravity, wall_collisions, particle_collisions)
        return particle
