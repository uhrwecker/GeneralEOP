import numpy as np
from scipy.integrate import odeint

from em import emitter
from photon import light, geodesics


class ODESolver:
    """
        Wrapper class for solving the differential equations for the motion of timelike geodesics.
    """
    def __init__(self, s, theta0, phi0, chi, alpha, beta, r=None, rotation='positive',
                 start=0, stop=7, num=5000, abserr=1e-7, relerr=1e-7):
        """
            :param s: float; spin of the timelike object
            :param theta0: float; initial angle theta of the lightlike particle
            :param phi0: float; initial angle phi of the lightlike particle
            :param chi: float; scaling parameter for the photon momentum, invariant for the application
            :param alpha: float; phi-like emission angle
            :param beta: float; theta-like emission angle
            :param r: float; initial radius parameter of the circular orbit of the timelike object
            :param rotation: rotation: ['positive', 'negative']; determining equations are oblivious to the sign of L;
                   this is determined by the rotation parameter.
            :param start: float; start of the affine parameter range for integration
            :param stop: float; stop of the affine parameter range for integration
            :param num: int; number of steps in integration
            :param abserr: float; least absolute error for integration
            :param relerr: float; least relative error for integration
        """

        self.emitter = emitter.EmitterProperties(s, r, rotation)
        self.photon = light.PhotonProperties(self.emitter, chi, alpha, beta)

        self.start = start
        self.stop = stop
        self.num = num

        self.sigma = np.linspace(start, stop, num=num)

        self.abserr = abserr
        self.relerr = relerr

        self.t0 = 0
        self.r0 = self.emitter.r0
        self.theta0 = theta0
        self.phi0 = phi0

        self.dt, self.dr, self.dtheta, self.dphi = self.photon.get_ic()

    def solve(self):
        """
            Main routine for solving with the previously specified initial conditions
            :return: iter; [sigma, result] where sigma is the array of affine parameter, and result includes all [x, x']
        """
        psi = np.array([self.t0, self.dt, self.r0, self.dr, self.theta0, self.dtheta, self.phi0, self.dphi])

        result = odeint(geodesics.geod, psi, self.sigma, atol=self.abserr, rtol=self.relerr)
        if (result[:, 2] < 2.5).any():
            result[:, 2] = np.array([np.nan for n in range(len(result[:, 2]))])

        return self.sigma, result

    def get_data(self):
        """
            Routine to get all important data in a dictionary, for saving purposes.
            :return: dict; dictionary that includes all important hyperparameters.
        """
        s = self.emitter.get_s()
        vr, vphi, gamma = self.emitter.get_velocities()
        Eph, Lph, Qph = self.photon.get_com()
        Eem, Lem = self.emitter.get_com()
        rotation = self.emitter.get_sense_of_rotation()

        data = {}
        data['Emitter'] = {}
        data['Photon'] = {}
        data['Observer'] = {}
        data['Numerics'] = {}

        data['Emitter']['s'] = s
        data['Emitter']['v_r'] = vr
        data['Emitter']['v_phi'] = vphi
        data['Emitter']['gamma'] = gamma
        data['Emitter']['E_emitter'] = Eem
        data['Emitter']['J_emitter'] = Lem
        data['Emitter']['emitter_sense_of_rotation'] = rotation

        data['Photon']['E_photon'] = Eph
        data['Photon']['L_photon'] = Lph
        data['Photon']['Q_photon'] = Qph
        data['Photon']['t0'] = self.t0
        data['Photon']['r0'] = self.r0
        data['Photon']['theta0'] = self.theta0
        data['Photon']['phi0'] = self.phi0

        data['Numerics']['solver_abs_err'] = self.abserr
        data['Numerics']['solver_rel_err'] = self.relerr
        data['Numerics']['setpoints'] = self.num

        return data

    def set_alpha(self, alpha, recalc=True):
        """
            Set (and possibly recalc) the emission angle alpha.
            :param alpha: phi-like emission angle
            :param recalc: bool; determines if the COM and physical velocities are recalculated for the given s.
        """
        self.photon.set_alpha(alpha, recalc)
        if recalc:
            self.dt, self.dr, self.dtheta, self.dphi = self.photon.get_ic()

    def set_beta(self, beta, recalc=True):
        """
            Set (and possibly recalc) the emission angle beta.
            :param beta: theta-like emission angle
            :param recalc: bool; determines if the COM and physical velocities are recalculated for the given s.
        """
        self.photon.set_beta(beta, recalc)
        if recalc:
            self.dt, self.dr, self.dtheta, self.dphi = self.photon.get_ic()