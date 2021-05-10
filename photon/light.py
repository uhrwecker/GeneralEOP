import numpy as np


class PhotonProperties:
    def __init__(self, emitter, chi=5, alpha=0, beta=np.pi / 2):
        self.emitter = emitter

        self.chi = chi
        self.alpha = alpha
        self.beta = beta

        self.E = None
        self.L = None
        self.Q = None

        self.dt = None
        self.dr = None
        self.dtheta = None
        self.dphi = None

        self.setup()

    def setup(self):
        self.E, self.L, self.Q = self.calculate_com()

        self.dt, self.dr, self.dtheta, self.dphi = self.calculate_ic()

    def calculate_com(self):
        vr, vphi, gamma = self.emitter.get_velocities()
        r0 = self.emitter.get_r0()

        E = self._E(self.chi, gamma, r0, vr, vphi, self.alpha, self.beta)
        L = self._L(self.chi, gamma, r0, vr, vphi, self.alpha, self.beta, np.pi / 2)
        Q = self._Q(r0, self.beta, L, np.pi / 2)

        return E, L, Q

    def calculate_ic(self):
        r = self.emitter.get_r0()

        dt = self.E / (1 - 2 / r)
        omega = self.Q - self.L ** 2 * (np.cos(np.pi / 2) / np.sin(np.pi / 2)) ** 2
        if omega < 0:
            omega = np.abs(omega)
        dtheta = self.chi / r * np.cos(self.beta)#np.sqrt(omega) / r**2
        dphi = self.L / (r * np.sin(np.pi / 2)) ** 2

        dr = np.sqrt(self.E ** 2 - (1 - 2 / r) * (self.Q + self.L ** 2) / r ** 2)
        dr = self._check_dr_sign(dr)

        return dt, dr, dtheta, dphi

    def get_ic(self):
        return self.dt, self.dr, self.dtheta, self.dphi

    def get_com(self):
        return self.E, self.L, self.Q

    def get_angles(self):
        return self.chi, self.alpha, self.beta

    def set_alpha(self, alpha, recalc=True):
        self.alpha = alpha
        if recalc:
            self.setup()

    def set_beta(self, beta, recalc=True):
        self.beta = beta
        if recalc:
            self.setup()

    def _E(self, chi, gamma, r, vr, vphi, alpha, beta):
        factor = chi * gamma * np.sqrt(1 - 2 / r)

        return factor * (1 + vr * np.cos(alpha) * np.sin(beta) + vphi * np.sin(alpha) * np.sin(beta))

    def _L(self, chi, gamma, r, vr, vphi, alpha, beta, theta):
        f1 = gamma * vr * vphi / (1 + gamma)
        f2 = 1 / gamma + gamma * vphi ** 2 / (1 + gamma)

        factor = gamma * chi * r * np.sin(theta)

        return factor * (vphi + f1 * np.cos(alpha) * np.sin(beta) + f2 * np.sin(alpha) * np.sin(beta))

    def _Q(self, r, beta, L, theta):
        q = self.chi**2 * r ** 2 * np.cos(beta) ** 2 + L ** 2 / np.tan(theta) ** 2
        if q < 0 and np.abs(q) < 1e-10:
            return np.abs(q)
        else:
            return q

    def _check_dr_sign(self, dr):
        vr, vphi, gamma = self.emitter.get_velocities()
        r0 = self.emitter.get_r0()

        dr = np.sqrt(1 - 2 / r0) * self.chi * gamma * (
                vr + 1 / gamma * (1 + gamma ** 2 * vr ** 2 / (1 + gamma)) * np.cos(self.alpha) * np.sin(self.beta) \
                + gamma * vr * vphi * np.sin(self.alpha) * np.sin(self.beta) / (1 + gamma))

        return dr