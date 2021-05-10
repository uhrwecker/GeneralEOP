import numpy as np
import time
import pandas as pd
import configparser as cp

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def get_indices_of_k_smallest(arr, k):
    idx = np.argpartition(arr.ravel(), k)
    return tuple(np.array(np.unravel_index(idx, arr.shape))[:, range(min(k, 0), max(k, 0))])


class EmitterObserverProblem:
    def __init__(self, solver, r_obs, theta_obs, phi_obs):
        self.solver = solver
        self.robs = r_obs
        self.thetaobs = theta_obs
        self.phiobs = phi_obs

    def find_critical_angles(self, half='right', amin=None, amax=None, bmin=None, bmax=None):
        if not amin:
            amin = 0
        if not amax:
            amax = np.pi
        if not bmin:
            bmin = 0
        if not bmax:
            bmax = np.pi

        if half == 'left':
            amin += np.pi
            amax += np.pi

        start = time.time()

        alpha, beta, flag = self.bad_montecarlo(amin, amax, bmin, bmax, n=17, max_step=40)

        now = time.time()
        print(f'Finding the critical angle took {now-start}s.')

        return alpha, beta, flag


    def bad_montecarlo(self, amin, amax, bmin, bmax, n=3, max_step=50):
        converged = False
        once = False
        step = 0
        inc = 0.05

        last_ab = [0, 0, 0, 0]

        while not converged and step < max_step:
            parameters = [[(a, b) for a in np.linspace(amin, amax, endpoint=True, num=n)] for b in
                          np.linspace(bmin, bmax, endpoint=True, num=n)]
            amin, amax, bmin, bmax, flag = self._generate_solutions(parameters)

            if flag:
                print(f'Converged at step {step} / {max_step}!')
                print(f'- Alpha = {amin}, Beta = {bmin}.')
                converged = True
                break

            if amin in last_ab and amax in last_ab and bmin in last_ab and bmax in last_ab:
                guess_alpha = (amax + amin) / 2
                guess_beta = (bmax + bmin) / 2

                amin = guess_alpha - inc
                amax = guess_alpha + inc
                bmin = guess_beta - inc / 2
                bmax = guess_beta + inc / 2

                print('!!! Same alphas and betas as before.')

            if np.abs(amin - amax) < 1e-6 and np.abs(bmin - bmax) < 1e-6:
                amin -= inc
                amax += inc
                bmin -= inc/2
                bmax += inc/2

                inc /= 2
                if n < 3:
                    #amin = np.nanmean([amin, amax])
                    #bmin = np.nanmean([bmin, bmax])
                    #converged = True
                    n = 3
                #    print('Warning: n = 0 (?)')
                step -= 2
                if step < 0:
                    step = 0
                if inc < 1e-5 and not once:
                    inc = 0.1
                    n = 25
                    print('Retrying ...')
                    once = True
                    step -= 5
                elif inc < 1e-5 and once:
                    step = max_step

            last_ab[0] = amin
            last_ab[1] = amax
            last_ab[2] = bmin
            last_ab[3] = bmax
            step += 1

            if step % 5 == 0:
                print(f'Algo for phi0 = {self.solver.phi0} at step {step} / {max_step}.')
                print(f'- alpha between {amin} and {amax}.')
                print(f'- beta  between {bmin} and {bmax} \n')

            if step == max_step:
                print('The algorithm did not converge fast enough.\n')

        return amin, bmin, converged

    def save_to_csv(self, alpha, beta, folder='./'):
        start = time.time()

        self.solver.set_alpha(alpha, False)
        self.solver.set_beta(beta)

        sigma, data = self.solver.solve()

        df = pd.DataFrame({
            'sigma': sigma,
            't': data[:, 0],
            'dt': data[:, 1],
            'r': data[:, 2],
            'dr': data[:, 3],
            'theta': data[:, 4],
            'dtheta': data[:, 5],
            'phi': data[:, 6],
            'dphi': data[:, 7],
        })

        if beta > np.pi / 2:
            tag = folder + 'up_%.8f' % data[:, 6][0]
        else:
            tag = folder + 'down_%.8f' % data[:, 6][0]

        self._write_config(tag, alpha, beta)
        df.to_csv(tag+'.csv', index=False)

        print(f'Saving took {time.time() - start}s.')

    def _write_config(self, tag, alpha, beta):
        config = cp.ConfigParser()

        info = self.solver.get_data()
        info['Observer']['r_obs'] = self.robs
        info['Observer']['phi_obs'] = self.phiobs
        info['Observer']['theta_obs'] = self.thetaobs

        info['Emitter']['alpha'] = alpha
        info['Emitter']['beta'] = beta

        config['Emitter'] = info['Emitter']
        config['Observer'] = info['Observer']
        config['Photon'] = info['Photon']
        config['Numerics'] = info['Numerics']

        fp = tag + '.ini'
        with open(fp, 'w') as cf:
            config.write(cf)

    def _generate_solutions(self, params):
        m = []
        alpha = None
        beta = None

        for fixed_b in params:
            row = []
            for a, b in fixed_b:
                self.solver.set_alpha(a, False)
                self.solver.set_beta(b)

                _, data = self.solver.solve()

                r = data[:, 2]
                t = data[:, 4]
                p = data[:, 6]

                if self._check_if_at_observer(r, t, p):
                    alpha = a
                    beta = b
                    break

                idx_r = find_nearest(r, self.robs)

                dist = np.sqrt(((t[idx_r] % np.pi) - self.thetaobs)**2 + ((p[idx_r] % (2*np.pi)) - self.phiobs)**2)
                row.append(dist)
            m.append(row)

        if alpha and beta:
            return alpha, None, beta, None, True

        m = np.array(m)
        try:
            idx = get_indices_of_k_smallest(m, 4)
        except:
            print(m)
            raise ValueError
        para_matrix = np.array(params)[idx[:3]]
        #print(para_matrix)

        return np.amin(para_matrix[:, 0]), np.amax(para_matrix[:, 0]), np.amin(para_matrix[:, 1]), np.amax(para_matrix[:, 1]), False


    def _check_if_at_observer(self, r, t, p):

        tflag = False
        pflag = False

        idx_r = find_nearest(r, self.robs)

        #print(t[idx_r], np.pi/3, p[idx_r], np.pi/2)

        if np.isclose(t[idx_r] % (np.pi), self.thetaobs, atol=1e-4, rtol=1e-3):
            tflag = True

        if np.isclose(p[idx_r] % (2 * np.pi), self.phiobs, atol=1e-4, rtol=1e-3):
            pflag = True

        if tflag and pflag:
            return True
        else:
            return False
