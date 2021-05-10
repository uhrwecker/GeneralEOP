import numpy as np
from scipy.integrate import odeint
import os
import sys
import contextlib

from em import emitter
from photon import light, geodesics


def fileno(file_or_fd):
    fd = getattr(file_or_fd, 'fileno', lambda: file_or_fd)()
    if not isinstance(fd, int):
        raise ValueError("Expected a file (`.fileno()`) or a file descriptor")
    return fd


@contextlib.contextmanager
def stdout_redirected(to=os.devnull, stdout=None):
    """
    https://stackoverflow.com/a/22434262/190597 (J.F. Sebastian)
    """
    if stdout is None:
       stdout = sys.stdout

    stdout_fd = fileno(stdout)
    # copy stdout_fd before it is overwritten
    # NOTE: `copied` is inheritable on Windows when duplicating a standard stream
    with os.fdopen(os.dup(stdout_fd), 'wb') as copied:
        stdout.flush()  # flush library buffers that dup2 knows nothing about
        try:
            os.dup2(fileno(to), stdout_fd)  # $ exec >&to
        except ValueError:  # filename
            with open(to, 'wb') as to_file:
                os.dup2(to_file.fileno(), stdout_fd)  # $ exec > to
        try:
            yield stdout # allow code to be run with the redirected stdout
        finally:
            # restore stdout to its previous value
            # NOTE: dup2 makes stdout_fd inheritable unconditionally
            stdout.flush()
            os.dup2(copied.fileno(), stdout_fd)  # $ exec >&copied


class ODESolver:
    def __init__(self, s, theta0, phi0, chi, alpha, beta, r=None, rotation='positive',
                 start=0, stop=7, num=5000, abserr=1e-7, relerr=1e-7):

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
        psi = np.array([self.t0, self.dt, self.r0, self.dr, self.theta0, self.dtheta, self.phi0, self.dphi])

        #with stdout_redirected():
        result = odeint(geodesics.geod, psi, self.sigma, atol=self.abserr, rtol=self.relerr)
        if (result[:, 2] < 2.5).any():
            result[:, 2] = np.array([np.nan for n in range(len(result[:, 2]))])

        return self.sigma, result

    def get_data(self):
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
        self.photon.set_alpha(alpha, recalc)
        if recalc:
            self.dt, self.dr, self.dtheta, self.dphi = self.photon.get_ic()

    def set_beta(self, beta, recalc=True):
        self.photon.set_beta(beta, recalc)
        if recalc:
            self.dt, self.dr, self.dtheta, self.dphi = self.photon.get_ic()