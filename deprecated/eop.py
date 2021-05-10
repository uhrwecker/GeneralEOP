import numpy as np
import scipy.optimize as opt

from ode import algorithms


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx


class EmitterObserverProblem:
    def __init__(self, solver, r_obs, theta_obs, phi_obs):
        self.solver = solver
        self.robs = r_obs
        self.thetaobs = theta_obs
        self.phiobs = phi_obs

    def find_critical_angles(self, amin, amax, bmin, bmax):
        converged = False
        step = 0

        #alpha = np.nanmean([amin, amax])
        #beta = #np.nanmean([bmin, bmax])
        #self.solver.set_beta(beta)

        while not converged:
            alpha = self.find_critical_alpha(amin, amax, 0)

            self.solver.set_alpha(alpha)
            if self._check_if_at_observer():
                converged = True
                print(f'Converged at step {step}.')
                break

            beta = self.find_critical_beta(bmin, bmax, 0)
            self.solver.set_beta(beta)
            if self._check_if_at_observer():
                converged = True
                print(f'Converged at step {step}.')
                break

            if step == 100:
                print('Broke because to many steps')
                break

            step += 1

        return alpha, beta

    def _check_if_at_observer(self):
        _, data = self.solver.solve()

        tflag = False
        pflag = False

        r = data[:, 2]
        t = data[:, 4]
        p = data[:, 6]

        idx_r = find_nearest(r, self.robs)

        if np.isclose(t[idx_r] % (np.pi), self.thetaobs, rtol=1e-4):
            tflag = True

        if np.isclose(p[idx_r] % (2 * np.pi), self.phiobs, rtol=1e-4):
            pflag = True

        if tflag and pflag:
            return True
        else:
            return False

    def find_critical_alpha(self, mini, maxi, beta):
        #try:
        return algorithms.find_critical_angle(self._alpha_root, self.phiobs, mini, maxi, (self.robs, ), 'alpha')
        #except:
        #    raise ValueError(f'Alpha is not in range {mini} to {maxi} at beta = {beta}.')

    def find_critical_beta(self, mini, maxi, alpha):
        #try:
        return algorithms.find_critical_angle(self._beta_root, self.phiobs, mini, maxi, (self.robs,), 'beta')
        #except:
        #    raise ValueError(f'Beta is not in range {mini} to {maxi} for alpha = {alpha}.')

    def vary_beta(self, _):
        def _betaa(b, angle):
            self.solver.set_beta(b)
            _, data = self.solver.solve()

            r = data[:, 2]
            t = data[:, 4]
            p = data[:, 6]

            #x = r * np.cos(p) * np.sin(t)
            y = r * np.sin(p) * np.sin(t)
            z = r * np.cos(t)

            zn = y * np.sin(angle) + z * np.cos(angle)

            return np.mean(zn)

        _, data = self.solver.solve()
        yobs = self.robs * np.sin(self.phiobs) * np.sin(self.thetaobs)
        zobs = self.robs * np.cos(self.thetaobs)

        turn_angle = np.arctan(yobs/zobs)-np.pi/2

        beta = opt.root_scalar(_betaa, args=(turn_angle,), bracket=(0.1, np.pi), x0=3*np.pi/4).root

        return beta


    def vary_alpha(self, crit):
        # left
        idx_r = np.nan
        critical = np.nan
        alphas = np.linspace(np.pi, crit, num=50)
        for alpha in alphas:
            self.solver.set_alpha(alpha)

            _, data = self.solver.solve()

            if data[:, 2][-1] < 2:
                pass
            else:
                idx_r = find_nearest(data[:, 2], self.robs)
                critical = alpha
                break

        if idx_r == np.nan or critical == np.nan:
            raise ValueError()

        return algorithms.find_critical_angle(self._alpha_root, self.phiobs, critical, crit, (self.robs, ), 'alpha')



    def _beta_root(self, beta, robs):
        #self.solver.set_alpha(np.pi, False)
        self.solver.set_beta(beta)

        _, data = self.solver.solve()
        idx_r = find_nearest(data[:, 2], robs)

        if data[:, 2][-1] < 2:
            return 0

        return data[:, -4][idx_r]

    def _alpha_root(self, alpha, robs):
        self.solver.set_alpha(alpha)#, False)
        #beta = self.vary_beta(0)
        #self.solver.set_beta(beta)

        _, data = self.solver.solve()
        idx_r = find_nearest(data[:, 2], robs)

        if data[:, 2][-1] < 2:
            return 0

        return data[:, -2][idx_r]

