import numpy as np

from em import circular, isco


class EmitterProperties:
    """
        Class containing all emitter properties.
    """
    def __init__(self, s, r0=None, rotation='positive'):
        """
            :param s: float; spin of the timelike object
            :param r0: float; orbit of the timelike object
            :param rotation: rotation: ['positive', 'negative']; determining equations are oblivious to the sign of L;
            this is determined by the rotation parameter.
        """
        # spin of the emitter:
        self.s = s

        # sense of rotation:
        self.rotation = rotation

        self.r0 = r0
        self.E = None
        self.L = None

        # orbit velocities:
        self.vr = None
        self.vphi = None
        self.gamma = None

        self.setup()

    def setup(self):
        """
            Basic setup function.
        """
        self.E, self.L, self.r0 = self.calculate_com()

        self.vr, self.vphi, self.gamma = self.calculate_vel()

    def calculate_com(self):
        """
            Routine to calculate the constants of motion.
            Right now, only calculates the circular orbit at specified r0.
            :return: [E, L, r0]; orbit parameters.
        """

        # always calculate the isco, to compare:
        E, r0, L = isco.get_com(self.s, self.rotation)

        # if a distance is given, calculate the COM for the emitter at this distance
        if self.r0:
            E, r1, L = circular.get_com(self.s, self.r0, self.rotation)

            # if the distance is smaller than the ISCO, print warning
            if r0 > r1:
                print('Warning: your given radius is smaller than the corresponding ISCO.')

            r0 = r1

        return E, L, r0

    def calculate_vel(self):
        """
            Routine to encapsule the calculation of the physical velocities.
            :return: [vr, vphi, gamma]; return the radial and orbit velocities, as well as common gamma.
        """
        vphi = self._vphi()
        vr = 0  # self._vr(); only circular orbits considered, thus vr = 0

        return vr, vphi, 1/np.sqrt(1 - vr**2 - vphi**2)

    def set_s(self, s, recalc=True):
        """
            Set (and possibly recalc) the spin parameter
            :param s: float; spin of the timelike object
            :param recalc: bool; determines if the COM and physical velocities are recalculated for the given s.
        """
        self.s = s
        if recalc:
            self.setup()

    def set_r0(self, r0, recalc=True):
        """
            Set (and possibly recalc) the orbit radius parameter.
            :param r0: float; orbit of the timelike object
            :param recalc: bool; determines if the COM and physical velocities are recalculated for the given s.
        """
        self.r0 = r0
        if recalc:
            self.setup()

    def get_com(self):
        """
            Getter for the constants of motion.
            :return: [E, L]; energy and angular momentum of the orbit.
        """
        return self.E, self.L

    def get_r0(self):
        """
            Getter for the orbit radius parameter.
            :return: float; orbit of the timelike object.
        """
        return self.r0

    def get_velocities(self):
        """
            Getter for the physical velocities of the timelike object.
            :return: [vr, vphi, gamma]; physical velocities of the emitter, as well as well known gamma.
        """
        return self.vr, self.vphi, self.gamma

    def get_sense_of_rotation(self):
        """
            Getter for the sense of rotation.
            :return: ['positive', 'negative]; determines the sign of the angular momentum of the orbit.
        """
        return self.rotation

    def get_s(self):
        """
            Getter for the spin of the timelike object
            :return: float; spin of the timelike object.
        """
        return self.s

    def _vphi(self):
        """
            (private)
            Calculate the orbital velocity. This subroutine includes all equations necessary.
            :return: float; orbital velocity of the timelike object.
        """
        upper = (self.r0 ** 3 + 2 * self.s ** 2) * (self.L - self.s * self.E)
        lower = (self.r0 ** 3 * self.E - self.s * self.L) * (self.r0 ** 3 - self.s ** 2)

        return np.sqrt(1 - 2 / self.r0) * self.r0 ** 2 * upper / lower

    def _vr(self):
        """
            (private)
            Calculate the radial velocity. This subroutine includes all equations necessary.
            Note that this is 0 (most of the times), as only circular motion is considered.
            :return: float; radial velocity of the timelike object.
        """
        root = (self.r0 ** 2 * self.E - self.s / self.r0 * self.L) ** 2 - \
               self.r0 * (self.r0 - 2) * ((self.r0 ** 3 - self.s ** 2) ** 2 / self.r0 ** 4 +
                                          (self.L - self.s * self.E) ** 2)

        if root < 0:
            print(f'warning: root smaller 0; might be ok, root is {root}')  #

        return self.r0 / (self.r0 ** 3 * self.E - self.s * self.L) * np.sqrt(np.abs(root))
