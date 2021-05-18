import numpy as np


class PhotonProperties:
    """
        Class for encapsuling the photon properties.
    """
    def __init__(self, emitter, chi=1, alpha=0, beta=np.pi / 2):
        """
            :param emitter: object; em.emitter.EmitterProperties object
            :param chi: float; scaling parameter for the photon momentum, invariant for the application
            :param alpha: float; phi-like emission angle
            :param beta: float; theta-like emission angle
        """
        self.emitter = emitter

        self.chi = chi
        self.alpha = alpha
        self.beta = beta

        # COM:
        self.E = None
        self.L = None
        self.Q = None

        # initial velocities:
        self.dt = None
        self.dr = None
        self.dtheta = None
        self.dphi = None

        self.setup()

    def setup(self):
        """
            Simple setup routine.
        """
        self.E, self.L, self.Q = self.calculate_com()

        self.dt, self.dr, self.dtheta, self.dphi = self.calculate_ic()

    def calculate_com(self):
        """
            Routine to calculate the constants of motion for the photon geodesic.
            :return: iterable; [E, L, Q]
        """
        vr, vphi, gamma = self.emitter.get_velocities()
        r0 = self.emitter.get_r0()

        E = self._E(self.chi, gamma, r0, vr, vphi, self.alpha, self.beta)
        L = self._L(self.chi, gamma, r0, vr, vphi, self.alpha, self.beta, np.pi / 2)
        Q = self._Q(r0, self.beta, L, np.pi / 2)

        return E, L, Q

    def calculate_ic(self):
        """
            Routine to calculate the initial velocities for the photon geodesic.
            :return: iterable; [t', r', theta', phi']
        """
        r = self.emitter.get_r0()

        dt = self.E / (1 - 2 / r)
        omega = self.Q - self.L ** 2 * (np.cos(np.pi / 2) / np.sin(np.pi / 2)) ** 2
        if omega < 0:
            omega = np.abs(omega)
        dtheta = self.chi / r * np.cos(self.beta)
        dphi = self.L / (r * np.sin(np.pi / 2)) ** 2

        dr = np.sqrt(self.E ** 2 - (1 - 2 / r) * (self.Q + self.L ** 2) / r ** 2)
        dr = self._check_dr_sign(dr)

        return dt, dr, dtheta, dphi

    def get_ic(self):
        """
            Getter for the initial velocities.
            :return: iterable; [t', r', theta', phi']
        """
        return self.dt, self.dr, self.dtheta, self.dphi

    def get_com(self):
        """
            Getter for the constants of motion.
            :return: iterable; [E, L, Q]
        """
        return self.E, self.L, self.Q

    def get_angles(self):
        """
            Getter for the emission angles.
            :return: iterable; [chi, alpha, beta]
        """
        return self.chi, self.alpha, self.beta

    def set_alpha(self, alpha, recalc=True):
        """
            Set (and possibly recalc) the emission angle alpha.
            :param alpha: phi-like emission angle
            :param recalc: bool; determines if the COM and physical velocities are recalculated for the given s.
        """
        self.alpha = alpha
        if recalc:
            self.setup()

    def set_beta(self, beta, recalc=True):
        """
            Set (and possibly recalc) the emission angle beta.
            :param beta: theta-like emission angle
            :param recalc: bool; determines if the COM and physical velocities are recalculated for the given s.
        """
        self.beta = beta
        if recalc:
            self.setup()

    def _E(self, chi, gamma, r, vr, vphi, alpha, beta):
        """
            (private)
            Calculate the energy of the photon. This subroutine includes all equations necessary.
            :param chi: float; scaling parameter for the photon momentum, invariant for the application
            :param gamma: float; common gamma
            :param r: float; initial orbit radius parameter
            :param vr: float; radial velocity of the emitter
            :param vphi: float; orbital velocity of the emitter
            :param alpha: float; phi-like emission angle
            :param beta: float; theta-like emisson angle
            :return: float; energy of the photon geodesic
        """
        factor = chi * gamma * np.sqrt(1 - 2 / r)

        return factor * (1 + vr * np.cos(alpha) * np.sin(beta) + vphi * np.sin(alpha) * np.sin(beta))

    def _L(self, chi, gamma, r, vr, vphi, alpha, beta, theta):
        """
            (private)
            Calculate the angular momentum of the photon. This subroutine includes all equations necessary.
            :param chi: float; scaling parameter for the photon momentum, invariant for the application
            :param gamma: float; common gamma
            :param r: float; initial orbit radius parameter
            :param vr: float; radial velocity of the emitter
            :param vphi: float; orbital velocity of the emitter
            :param alpha: float; phi-like emission angle
            :param beta: float; theta-like emisson angle
            :param theta: float; initial theta parameter
            :return: float; angular momentum of the photon geodesic
        """
        f1 = gamma * vr * vphi / (1 + gamma)
        f2 = 1 / gamma + gamma * vphi ** 2 / (1 + gamma)

        factor = gamma * chi * r * np.sin(theta)

        return factor * (vphi + f1 * np.cos(alpha) * np.sin(beta) + f2 * np.sin(alpha) * np.sin(beta))

    def _Q(self, r, beta, L, theta):
        """
            (private)
            Calculate the Carter constant of the photon. This subroutine includes all equations necessary.
            :param r: float; initial orbit radius parameter
            :param beta: float; theta-like emission angle
            :param L: float; angular momentum constant of motion
            :param theta: float; initial theta parameter
            :return: float; Carter constant of motion
        """
        q = self.chi**2 * r ** 2 * np.cos(beta) ** 2 + L ** 2 / np.tan(theta) ** 2
        if q < 0 and np.abs(q) < 1e-10:
            return np.abs(q)
        else:
            return q

    def _check_dr_sign(self, dr):
        """
            (private)
            Check the sign of r'.
            :param dr: float; derivative of r
            :return: float; r' with proper sign
        """
        vr, vphi, gamma = self.emitter.get_velocities()
        r0 = self.emitter.get_r0()

        dr = np.sqrt(1 - 2 / r0) * self.chi * gamma * (
                vr + 1 / gamma * (1 + gamma ** 2 * vr ** 2 / (1 + gamma)) * np.cos(self.alpha) * np.sin(self.beta) \
                + gamma * vr * vphi * np.sin(self.alpha) * np.sin(self.beta) / (1 + gamma))
        return dr
