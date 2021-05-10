import numpy as np

from em import circular, isco


class EmitterProperties:
    def __init__(self, s, r0=None, rotation='positive'):
        # spin of the emitter
        self.s = s

        # sense of rotation
        self.rotation = rotation

        self.r0 = r0
        self.E = None
        self.L = None

        self.vr = None
        self.vphi = None
        self.gamma = None

        self.setup()

    def setup(self):
        self.E, self.L, self.r0 = self.calculate_com()

        self.vr, self.vphi, self.gamma = self.calculate_vel()

    def calculate_com(self):
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
        vphi = self._vphi()
        vr = 0  # self._vr(); only circular orbits considered, thus vr = 0

        return vr, vphi, 1/np.sqrt(1 - vr**2 - vphi**2)

    def set_s(self, s, recalc=True):
        self.s = s
        if recalc:
            self.setup()

    def set_r0(self, r0, recalc=True):
        self.r0 = r0
        if recalc:
            self.setup()

    def get_com(self):
        return self.E, self.L

    def get_r0(self):
        return self.r0

    def get_velocities(self):
        return self.vr, self.vphi, self.gamma

    def get_sense_of_rotation(self):
        return self.rotation

    def get_s(self):
        return self.s

    def _vphi(self):
        upper = (self.r0 ** 3 + 2 * self.s ** 2) * (self.L - self.s * self.E)
        lower = (self.r0 ** 3 * self.E - self.s * self.L) * (self.r0 ** 3 - self.s ** 2)

        return np.sqrt(1 - 2 / self.r0) * self.r0 ** 2 * upper / lower

    def _vr(self):
        root = (self.r0 ** 2 * self.E - self.s / self.r0 * self.L) ** 2 - \
               self.r0 * (self.r0 - 2) * ((self.r0 ** 3 - self.s ** 2) ** 2 / self.r0 ** 4 +
                                          (self.L - self.s * self.E) ** 2)

        if root < 0:
            print(f'warning: root smaller 0; might be ok, root is {root}')  #

        return self.r0 / (self.r0 ** 3 * self.E - self.s * self.L) * np.sqrt(np.abs(root))
