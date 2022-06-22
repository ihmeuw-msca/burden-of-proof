# nonlinear mixed effects model
import numpy as np
import ipopt
from copy import deepcopy
from limetr import utils


class LimeTr:
    def __init__(self, n, k_beta, k_gamma, Y, F, JF, Z,
                 S=None, share_obs_std=False,
                 C=None, JC=None, c=None,
                 H=None, JH=None, h=None,
                 uprior=None, gprior=None, lprior=None,
                 certain_inlier_id=None,
                 inlier_percentage=1.0):
        """
        Create LimeTr object, for general mixed effects model

        Parameters
        ----------
        n : ndarray
            study sizes, n[i] is the number of observation for ith study.
        k_beta : int
            dimension of beta
        k_gamma : int
            dimension of gamma
        Y : ndarray
            study observations
        F : function
            return the predict observations given beta
        JF : function
            return the jacobian function of F
        Z : ndarray
            covariates matrix for the random effect
        S : optional, ndarray
            observation standard deviation
        """
        # pass in the dimension
        self.n = np.array(n)
        self.m = len(n)
        self.N = sum(n)
        self.k_beta = k_beta
        self.k_gamma = k_gamma
        # if include measurement error also as variable
        if S is not None:
            self.std_flag = 0
            self.k_delta = 0
        elif share_obs_std:
            self.std_flag = 1
            self.k_delta = 1
        else:
            self.std_flag = 2
            self.k_delta = self.m

        self.k = self.k_beta + self.k_gamma + self.k_delta
        self.k_total = self.k

        self.idx_beta = slice(0, self.k_beta)
        self.idx_gamma = slice(self.k_beta, self.k_beta + self.k_gamma)
        self.idx_delta = slice(self.k_beta + self.k_gamma, self.k)
        self.idx_split = np.cumsum(np.insert(n, 0, 0))[:-1]

        # pass in the data
        self.Y = Y
        self.F = F
        self.JF = JF
        self.Z = Z
        self.S = S
        if self.std_flag == 0:
            self.V = S**2

        # pass in the priors
        self.use_constraints = (C is not None)
        self.use_regularizer = (H is not None)
        self.use_uprior = (uprior is not None)
        self.use_gprior = (gprior is not None)
        self.use_lprior = (lprior is not None)

        self.C = C
        self.JC = JC
        self.c = c
        if self.use_constraints:
            self.constraints = C
            self.jacobian = JC
            self.num_constraints = C(np.zeros(self.k)).size
            self.cl = c[0]
            self.cu = c[1]
        else:
            self.num_constraints = 0
            self.cl = []
            self.cu = []

        self.H = H
        self.JH = JH
        self.h = h
        if self.use_regularizer:
            self.num_regularizer = H(np.zeros(self.k)).size
            self.hm = self.h[0]
            self.hw = 1.0/self.h[1]**2
        else:
            self.num_regularizer = 0

        if self.use_uprior:
            self.uprior = uprior
        else:
            self.uprior = np.array([
                [-np.inf]*self.k_beta + [0.0]*self.k_gamma +
                [1e-7]*self.k_delta,
                [np.inf]*self.k
            ])
            self.use_uprior = True

        self.lb = self.uprior[0]
        self.ub = self.uprior[1]

        if self.use_gprior:
            self.gprior = gprior
            self.gm = gprior[0]
            self.gw = 1.0/gprior[1]**2

        if self.use_lprior:
            self.lprior = lprior
            self.lm = lprior[0]
            self.lw = np.sqrt(2.0)/lprior[1]

            # double dimension pass into ipopt
            self.k_total += self.k

            # extend the constraints matrix
            if self.use_constraints:
                def constraints(x):
                    v = x[:self.k]
                    v_abs = x[self.k:]

                    vec1 = C(v)
                    vec2 = np.hstack((v_abs - (v - self.lm),
                                      v_abs + (v - self.lm)))

                    return np.hstack((vec1, vec2))

                def jacobian(x):
                    v = x[:self.k]
                    v_abs = x[self.k:]
                    Id = np.eye(self.k)

                    mat1 = JC(v)
                    mat2 = np.block([[-Id, Id], [Id, Id]])

                    return np.vstack((mat1, mat2))
            else:
                def constraints(x):
                    v = x[:self.k]
                    v_abs = x[self.k:]

                    vec = np.hstack((v_abs - v, v_abs + v))

                    return vec

                def jacobian(x):
                    v = x[:self.k]
                    v_abs = x[self.k:]
                    Id = np.eye(self.k)

                    mat = np.block([[-Id, Id], [Id, Id]])

                    return mat

            self.num_constraints += 2*self.k
            self.constraints = constraints
            self.jacobian = jacobian
            self.cl = np.hstack((self.cl, np.zeros(2*self.k)))
            self.cu = np.hstack((self.cu, np.repeat(np.inf, 2*self.k)))

            # extend the regularizer matrix
            if self.use_regularizer:
                def H_new(x):
                    v = x[:self.k]

                    return H(v)

                def JH_new(x):
                    v = x[:self.k]

                    return np.hstack((JH(v),
                                      np.zeros((self.num_regularizer,
                                                self.k))))

                self.H = H_new
                self.JH = JH_new

            # extend Uniform priors
            if self.use_uprior:
                uprior_abs = np.array([[0.0]*self.k, [np.inf]*self.k])
                self.uprior = np.hstack((self.uprior, uprior_abs))
                self.lb = self.uprior[0]
                self.ub = self.uprior[1]

        # trimming option
        self.use_trimming = (0.0 < inlier_percentage < 1.0)
        self.certain_inlier_id = certain_inlier_id
        self.inlier_percentage = inlier_percentage
        self.num_inliers = np.floor(inlier_percentage*self.N)
        self.num_outliers = self.N - self.num_inliers
        self.w = np.repeat(self.num_inliers/self.N, self.N)

        if self.certain_inlier_id is not None:
            self.certain_inlier_id = np.unique(self.certain_inlier_id)
            self.active_trimming_id = np.array(
                [i
                 for i in range(self.N)
                 if i not in self.certain_inlier_id])
        else:
            self.active_trimming_id = None

        # specify solution to be None
        self.soln = None
        self.info = None
        self.beta = np.zeros(self.k_beta)
        self.gamma = np.repeat(0.01, self.k_gamma)
        self.delta = np.repeat(0.01, self.k_delta)

        # check the input
        self.check()

    def check(self):
        assert self.Y.shape == (self.N,)
        assert self.Z.shape == (self.N, self.k_gamma)
        if self.S is not None:
            assert self.S.shape == (self.N,)
            assert np.all(self.S > 0.0)

        if self.use_constraints:
            assert self.c.shape == (2, self.num_constraints)
            assert np.all(self.cl <= self.cu)

        if self.use_regularizer:
            assert self.h.shape == (2, self.num_regularizer)
            assert np.all(self.h[1] > 0.0)

        if self.use_uprior:
            assert np.all(self.lb <= self.ub)

        if self.use_gprior:
            assert np.all(self.gprior[1] > 0.0)

        assert 0.0 < self.inlier_percentage <= 1.0
        if self.use_trimming and self.certain_inlier_id is not None:
            assert isinstance(self.certain_inlier_id, np.ndarray)
            assert self.certain_inlier_id.dtype == int
            assert self.certain_inlier_id.ndim == 1
            assert self.certain_inlier_id.size <= self.num_inliers
            assert np.min(self.certain_inlier_id) >= 0
            assert np.max(self.certain_inlier_id) < self.N

        if self.k > self.N:
            print('Warning: information insufficient!')

    def objective(self, x, use_ad=False):
        # unpack variable
        beta = x[self.idx_beta]
        gamma = x[self.idx_gamma]
        delta = x[self.idx_delta]

        gamma[gamma <= 0.0] = 0.0

        # trimming option
        if self.use_trimming:
            sqrt_w = np.sqrt(self.w)
            sqrt_W = sqrt_w.reshape(self.N, 1)
            F_beta = self.F(beta)*sqrt_w
            Y = self.Y*sqrt_w
            Z = self.Z*sqrt_W
            if self.std_flag == 0:
                V = self.V**self.w
            elif self.std_flag == 1:
                V = np.repeat(delta[0], self.N)**self.w
            elif self.std_flag == 2:
                V = np.repeat(delta, self.n)**self.w
        else:
            F_beta = self.F(beta)
            Y = self.Y
            Z = self.Z
            if self.std_flag == 0:
                V = self.V
            elif self.std_flag == 1:
                V = np.repeat(delta[0], self.N)
            elif self.std_flag == 2:
                V = np.repeat(delta, self.n)

        # residual and variance
        R = Y - F_beta
        D = utils.VarMat(V, Z, gamma, self.n)

        val = 0.5*self.N*np.log(2.0*np.pi)

        if use_ad:
            # should only use when testing
            varmat = D
            D = varmat.varMat()
            inv_D = varmat.invVarMat()

            val += 0.5*np.log(np.linalg.det(D))
            val += 0.5*R.dot(inv_D.dot(R))
        else:
            val += 0.5*D.logDet()
            val += 0.5*R.dot(D.invDot(R))

        # add gpriors
        if self.use_regularizer:
            val += 0.5*self.hw.dot((self.H(x) - self.hm)**2)

        if self.use_gprior:
            val += 0.5*self.gw.dot((x[:self.k] - self.gm)**2)

        if self.use_lprior:
            val += self.lw.dot(x[self.k:])

        return val

    def gradient(self, x, use_ad=False, eps=1e-12):
        if use_ad:
            # should only use when testing
            g = np.zeros(self.k_total)
            z = x + 0j
            for i in range(self.k_total):
                z[i] += eps*1j
                g[i] = self.objective(z, use_ad=use_ad).imag/eps
                z[i] -= eps*1j

            return g

        # unpack variable
        beta = x[self.idx_beta]
        gamma = x[self.idx_gamma]
        delta = x[self.idx_delta]

        gamma[gamma <= 0.0] = 0.0

        # trimming option
        if self.use_trimming:
            sqrt_w = np.sqrt(self.w)
            sqrt_W = sqrt_w.reshape(self.N, 1)
            F_beta = self.F(beta)*sqrt_w
            JF_beta = self.JF(beta)*sqrt_W
            Y = self.Y*sqrt_w
            Z = self.Z*sqrt_W
            if self.std_flag == 0:
                V = self.V**self.w
            elif self.std_flag == 1:
                V = np.repeat(delta[0], self.N)**self.w
            elif self.std_flag == 2:
                V = np.repeat(delta, self.n)**self.w
        else:
            F_beta = self.F(beta)
            JF_beta = self.JF(beta)
            Y = self.Y
            Z = self.Z
            if self.std_flag == 0:
                V = self.V
            elif self.std_flag == 1:
                V = np.repeat(delta[0], self.N)
            elif self.std_flag == 2:
                V = np.repeat(delta, self.n)

        # residual and variance
        R = Y - F_beta
        D = utils.VarMat(V, Z, gamma, self.n)

        # gradient for beta
        DR = D.invDot(R)
        g_beta = -JF_beta.T.dot(DR)

        # gradient for gamma
        DZ = D.invDot(Z)
        g_gamma = 0.5*np.sum(Z*DZ, axis=0) -\
            0.5*np.sum(
                np.add.reduceat(DZ.T*R, self.idx_split, axis=1)**2,
                axis=1)

        # gradient for delta
        if self.std_flag == 0:
            g_delta = np.array([])
        elif self.std_flag == 1:
            d = -DR**2 + D.invDiag()
            if self.use_trimming:
                v = np.repeat(delta[0], self.N)
                d *= self.w*(v**(self.w - 1.0))
            g_delta = 0.5*np.array([np.sum(d)])
        elif self.std_flag == 2:
            d = -DR**2 + D.invDiag()
            if self.use_trimming:
                v = np.repeat(delta, self.n)
                d *= self.w*(v**(self.w - 1.0))
            g_delta = 0.5*(np.add.reduceat(d, self.idx_split))

        g = np.hstack((g_beta, g_gamma, g_delta))

        # add gradient from the regularizer
        if self.use_regularizer:
            g += self.JH(x).T.dot((self.H(x) - self.hm)*self.hw)

        # add gradient from the gprior
        if self.use_gprior:
            g += (x[:self.k] - self.gm)*self.gw

        # add gradient from the lprior
        if self.use_lprior:
            g = np.hstack((g, self.lw))

        return g

    def objectiveTrimming(self, w):
        t = (self.Z**2).dot(self.gamma)
        r = self.Y - self.F(self.beta)
        if self.std_flag == 0:
            v = self.V
        elif self.std_flag == 1:
            v = np.repeat(self.delta[0], self.N)
        elif self.std_flag == 2:
            v = np.repeat(self.delta, self.n)
        d = v + t

        val = 0.5*np.sum(r**2*w/d)
        val += 0.5*self.N*np.log(2.0*np.pi) + 0.5*w.dot(np.log(d))

        return val

    def gradientTrimming(self, w, use_ad=False, eps=1e-10):
        if use_ad:
            # only use when testing
            g = np.zeros(self.N)
            z = w + 0j
            for i in range(self.N):
                z[i] += eps*1j
                g[i] = self.objectiveTrimming(z).imag/eps
                z[i] -= eps*1j

            return g

        t = (self.Z**2).dot(self.gamma)
        r = (self.Y - self.F(self.beta))**2
        if self.std_flag == 0:
            v = self.V
        elif self.std_flag == 1:
            v = np.repeat(self.delta[0], self.N)
        elif self.std_flag == 2:
            v = np.repeat(self.delta, self.n)
        d = v + t

        g = 0.5*r/d
        g += 0.5*np.log(d)

        return g

    def optimize(self, x0=None, print_level=0, max_iter=100, tol=1e-8,
                 acceptable_tol=1e-6,
                 nlp_scaling_method=None,
                 nlp_scaling_min_value=None):
        if x0 is None:
            x0 = np.hstack((self.beta, self.gamma, self.delta))
            if self.use_lprior:
                x0 = np.hstack((x0, np.zeros(self.k)))

        assert x0.size == self.k_total

        opt_problem = ipopt.problem(
            n=int(self.k_total),
            m=int(self.num_constraints),
            problem_obj=self,
            lb=self.uprior[0],
            ub=self.uprior[1],
            cl=self.cl,
            cu=self.cu
        )

        opt_problem.addOption('print_level', print_level)
        opt_problem.addOption('max_iter', max_iter)
        opt_problem.addOption('tol', tol)
        opt_problem.addOption('acceptable_tol', acceptable_tol)
        if nlp_scaling_method is not None:
            opt_problem.addOption('nlp_scaling_method', nlp_scaling_method)
        if nlp_scaling_min_value is not None:
            opt_problem.addOption('nlp_scaling_min_value', nlp_scaling_min_value)

        soln, info = opt_problem.solve(x0)

        self.soln = soln
        self.info = info
        self.beta = soln[self.idx_beta]
        self.gamma = soln[self.idx_gamma]
        self.delta = soln[self.idx_delta]

    def fitModel(self, x0=None,
                 inner_print_level=0,
                 inner_max_iter=20,
                 inner_tol=1e-8,
                 inner_acceptable_tol=1e-6,
                 inner_nlp_scaling_method=None,
                 inner_nlp_scaling_min_value=None,
                 outer_verbose=False,
                 outer_max_iter=100,
                 outer_step_size=1.0,
                 outer_tol=1e-6,
                 normalize_trimming_grad=False):

        if not self.use_trimming:
            self.optimize(x0=x0,
                          print_level=inner_print_level,
                          max_iter=inner_max_iter,
                          acceptable_tol=inner_acceptable_tol,
                          nlp_scaling_method=inner_nlp_scaling_method,
                          nlp_scaling_min_value=inner_nlp_scaling_min_value)

            return self.beta, self.gamma, self.w

        self.soln = x0

        num_iter = 0
        err = outer_tol + 1.0

        while err >= outer_tol:
            self.optimize(x0=self.soln,
                          print_level=inner_print_level,
                          max_iter=inner_max_iter,
                          tol=inner_tol,
                          acceptable_tol=inner_acceptable_tol,
                          nlp_scaling_method=inner_nlp_scaling_method,
                          nlp_scaling_min_value=inner_nlp_scaling_min_value)

            w_grad = self.gradientTrimming(self.w)
            if normalize_trimming_grad:
                w_grad /= np.linalg.norm(w_grad)
            w_new = utils.projCappedSimplex(
                self.w - outer_step_size*w_grad,
                self.num_inliers,
                active_id=self.active_trimming_id)

            err = np.linalg.norm(w_new - self.w)/outer_step_size
            np.copyto(self.w, w_new)

            num_iter += 1

            if outer_verbose:
                obj = self.objectiveTrimming(self.w)
                print('iter %4d, obj %8.2e, err %8.2e' % (num_iter, obj, err))

            if num_iter >= outer_max_iter:
                print('reach max outer iter')
                break

        return self.beta, self.gamma, self.w

    def estimateRE(self):
        """
        estimate random effect after fitModel
        """
        if self.soln is None:
            print('Please fit the model first.')
            return None

        if self.std_flag == 0:
            S = self.S
        elif self.std_flag == 1:
            S = np.sqrt(np.repeat(self.delta[0], self.N))
        elif self.std_flag == 2:
            S = np.sqrt(np.repeat(self.delta, self.n))

        if self.use_trimming:
            R = (self.Y - self.F(self.beta))*np.sqrt(self.w)
            Z = self.Z*(np.sqrt(self.w).reshape(self.N, 1))
        else:
            R = self.Y - self.F(self.beta)
            Z = self.Z

        iV = 1.0/S**2
        iVZ = Z*iV.reshape(iV.size, 1)

        r = np.split(R, np.cumsum(self.n)[:-1])
        v = np.split(S**2, np.cumsum(self.n)[:-1])
        z = np.split(Z, np.cumsum(self.n)[:-1], axis=0)
        ivz = np.split(iVZ, np.cumsum(self.n)[:-1], axis=0)

        u = []
        for i in range(self.m):
            rhs = ivz[i].T.dot(r[i])
            tmp = z[i]*self.gamma
            mat = np.diag(v[i]) + tmp.dot(z[i].T)
            vec = self.gamma*rhs - tmp.T.dot(np.linalg.solve(mat, tmp.dot(rhs)))
            u.append(vec)

        self.u = np.vstack(u)

        return self.u

    def get_re_vcov(self):
        if self.soln is None:
            print('Please fit the model first.')
            return None

        if self.std_flag == 0:
            S = self.S
        elif self.std_flag == 1:
            S = np.sqrt(np.repeat(self.delta[0], self.N))
        elif self.std_flag == 2:
            S = np.sqrt(np.repeat(self.delta, self.n))

        if self.use_trimming:
            R = (self.Y - self.F(self.beta))*np.sqrt(self.w)
            Z = self.Z*(np.sqrt(self.w).reshape(self.N, 1))
        else:
            R = self.Y - self.F(self.beta)
            Z = self.Z

        r = np.split(R, np.cumsum(self.n)[:-1])
        v = np.split(S**2, np.cumsum(self.n)[:-1])
        z = np.split(Z, np.cumsum(self.n)[:-1], axis=0)

        vcov = []
        for i in range(self.m):
            hessian = (z[i].T/v[i]).dot(z[i]) + np.diag(1.0/self.gamma)
            vcov.append(np.linalg.pinv(hessian))

        return vcov

    def estimate_re(self,
                    beta: np.ndarray = None,
                    gamma: np.ndarray = None,
                    use_gamma: bool = True) -> np.ndarray:
        beta = self.beta if beta is None else beta
        gamma = self.gamma if gamma is None else gamma
        r = np.split(self.Y - self.F(beta), np.cumsum(self.n)[:-1])
        z = np.split(self.Z, np.cumsum(self.n)[:-1], axis=0)
        v = np.split(self.S**2, np.cumsum(self.n)[:-1])

        u = []
        for i in range(self.m):
            rhs = (z[i].T/v[i]).dot(r[i])
            if use_gamma:
                q = (z[i].T/v[i]).dot(z[i])*gamma + np.identity(self.k_gamma)
                u.append(gamma[:, None]*np.linalg.inv(q).dot(rhs))
            else:
                q = (z[i].T/v[i]).dot(z[i])
                u.append(np.linalg.inv(q).dot(rhs))

        return np.vstack(u)

    def get_gamma_fisher(self, gamma: np.ndarray) -> np.ndarray:
        z = np.split(self.Z, np.cumsum(self.n)[:-1], axis=0)
        v = np.split(self.S**2, np.cumsum(self.n)[:-1])
        H = np.zeros((self.k_gamma, self.k_gamma))
        for i in range(self.m):
            q = np.diag(v[i]) + (z[i]*gamma).dot(z[i].T)
            q = z[i].T.dot(np.linalg.inv(q).dot(z[i]))
            H += 0.5*(q**2)
        return H

    def simulateData(self, beta_t, gamma_t, sim_prior=True, sim_re=True):
        # sample random effects and measurement error
        if sim_re:
            u = np.random.randn(self.m, self.k_gamma)*np.sqrt(gamma_t)
        else:
            if not hasattr(self, 'u'):
                self.estimateRE()
            u = self.u

        U = np.repeat(u, self.n, axis=0)
        ZU = np.sum(self.Z*U, axis=1)

        if self.std_flag == 0:
            S = self.S
        elif self.std_flag == 1:
            S = np.sqrt(np.repeat(self.delta[0], self.N))
        elif self.std_flag == 2:
            S = np.sqrt(np.repeat(self.delta, self.n))

        E = np.random.randn(self.N)*S

        self.Y = self.F(beta_t) + ZU + E

        if sim_prior:
            if self.use_gprior:
                valid_id = ~np.isinf(self.gprior[1])
                valid_num = np.sum(valid_id)
                self.gm = np.zeros(self.k)
                self.gm[valid_id] = self.gprior[0][valid_id] +\
                    np.random.randn(valid_num)*self.gprior[1][valid_id]

            if self.use_regularizer:
                valid_id = ~np.isinf(self.h[1])
                valid_num = np.sum(valid_id)
                self.hm = np.zeros(self.num_regularizer)
                self.hm[valid_id] = self.h[0][valid_id] +\
                    np.random.randn(valid_num)*self.h[1][valid_id]

    @classmethod
    def testProblem(cls,
                    use_trimming=False,
                    use_constraints=False,
                    use_regularizer=False,
                    use_uprior=False,
                    use_gprior=False,
                    know_obs_std=True,
                    share_obs_std=False):
        m = 10
        n = [5]*m
        N = sum(n)
        k_beta = 3
        k_gamma = 2
        if know_obs_std:
            k_delta = 0
        elif share_obs_std:
            k_delta = 1
        else:
            k_delta = m
        k = k_beta + k_gamma + k_delta

        beta_t = np.random.randn(k_beta)
        gamma_t = np.random.rand(k_gamma)*0.09 + 0.01

        X = np.random.randn(N, k_beta)
        Z = np.random.randn(N, k_gamma)

        S = np.random.rand(N)*0.09 + 0.01
        V = S**2
        D = np.diag(V) + (Z*gamma_t).dot(Z.T)

        U = np.random.multivariate_normal(np.zeros(N), D)
        E = np.random.randn(N)*S

        Y = X.dot(beta_t) + U + E

        def F(beta, X=X):
            return X.dot(beta)

        def JF(beta, X=X):
            return X

        # constraints, regularizer and priors
        if use_constraints:
            M = np.ones((1, k))

            def C(x, M=M):
                return M.dot(x)

            def JC(x, M=M):
                return M

            c = np.array([[0.0], [1.0]])
        else:
            C, JC, c = None, None, None

        if use_regularizer:
            M = np.ones((1, k))

            def H(x, M=M):
                return M.dot(x)

            def JH(x, M=M):
                return M

            h = np.array([[0.0], [2.0]])
        else:
            H, JH, h = None, None, None

        if use_uprior:
            uprior = np.array([[0.0]*k, [np.inf]*k])
        else:
            uprior = None

        if use_gprior:
            gprior = np.array([[0.0]*k, [2.0]*k])
        else:
            gprior = None

        if use_trimming:
            inlier_percentage = 0.9
        else:
            inlier_percentage = 1.0

        if not know_obs_std:
            S = None

        return cls(n, k_beta, k_gamma, Y, F, JF, Z, S=S,
                   C=C, JC=JC, c=c,
                   H=H, JH=JH, h=h,
                   uprior=uprior, gprior=gprior,
                   inlier_percentage=inlier_percentage,
                   share_obs_std=share_obs_std)

    @classmethod
    def testProblemLasso(cls):
        m = 100
        n = [1]*m
        N = sum(n)
        k_beta = 150
        k_gamma = 1
        k = k_beta + k_gamma

        beta_t = np.zeros(k_beta)
        beta_t[np.random.choice(k_beta, 5)] = np.sign(np.random.randn(5))
        gamma_t = np.zeros(k_gamma)

        X = np.random.randn(N, k_beta)
        Z = np.ones((N, k_gamma))
        Y = X.dot(beta_t)
        S = np.repeat(1.0, N)

        weight = 0.1*np.linalg.norm(X.T.dot(Y), np.inf)

        def F(beta):
            return X.dot(beta)

        def JF(beta):
            return X

        uprior = np.array([[-np.inf]*k_beta + [0.0], [np.inf]*k_beta + [0.0]])
        lprior = np.array([[0.0]*k, [np.sqrt(2.0)/weight]*k])

        return cls(n, k_beta, k_gamma, Y, F, JF, Z, S=S,
                   uprior=uprior, lprior=lprior)

    @staticmethod
    def sampleSoln(lt, sample_size=1, print_level=0, max_iter=100,
                   sim_prior=True, sim_re=True):
        beta_samples = np.zeros((sample_size, lt.k_beta))
        gamma_samples = np.zeros((sample_size, lt.k_gamma))

        beta_t = lt.beta.copy()
        gamma_t = lt.gamma.copy()

        lt_copy = deepcopy(lt)
        lt_copy.uprior[:, lt.k_beta:] = np.vstack((gamma_t, gamma_t))

        for i in range(sample_size):
            lt_copy.simulateData(beta_t, gamma_t,
                                 sim_prior=sim_prior, sim_re=sim_re)
            lt_copy.optimize(x0=np.hstack((beta_t, gamma_t)),
                             print_level=print_level,
                             max_iter=max_iter)

            u_samples = lt_copy.estimateRE()

            beta_samples[i] = lt_copy.beta.copy()
            gamma_samples[i] = np.maximum(
                lt.uprior[0, lt.k_beta:],
                np.minimum(lt.uprior[1, lt.k_beta:], np.var(u_samples, axis=0))
            )

            print('sampling solution progress %0.2f' % ((i + 1)/sample_size),
                  end='\r')

        return beta_samples, gamma_samples


def get_baseline_model(model: LimeTr):
    n = np.copy(model.n)
    k_beta = 1
    k_gamma = 1
    Y = model.Y.copy()
    intercept = np.ones((Y.size, 1))
    def F(beta): return intercept.dot(beta)
    def JF(beta): return intercept
    Z = intercept
    S = model.S.copy()
    w = model.w.copy()

    baseline_model = LimeTr(n, k_beta, k_gamma, Y, F, JF, Z, S=S)
    baseline_model.optimize()
    return baseline_model


def get_fe_pred(model: LimeTr):
    return model.F(model.beta)


def get_re_pred(model: LimeTr):
    re = model.estimateRE()
    return np.sum(model.Z*np.repeat(re, model.n, axis=0), axis=1)


def get_varmat(model: LimeTr):
    S = model.S**model.w
    Z = model.Z*np.sqrt(model.w)[:, None]
    n = model.n
    gamma = model.gamma
    return utils.VarMat(S**2, Z, gamma, n)


def get_marginal_rvar(model: LimeTr):
    residual = (model.Y - get_fe_pred(model))*np.sqrt(model.w)
    varmat = get_varmat(model)
    return residual.dot(varmat.invDot(residual))/model.w.sum()


def get_conditional_rvar(model: LimeTr):
    residual = (model.Y - get_fe_pred(model) - get_re_pred(model))*np.sqrt(model.w)
    varvec = (model.S**model.w)**2
    return (residual/varvec).dot(residual)/model.w.sum()


def get_marginal_R2(model: LimeTr,
                    baseline_model: LimeTr = None) -> float:
    if baseline_model is None:
        baseline_model = get_baseline_model(model)

    return 1 - get_marginal_rvar(model)/get_marginal_rvar(baseline_model)


def get_conditional_R2(model: LimeTr,
                       baseline_model: LimeTr = None):
    if baseline_model is None:
        baseline_model = get_baseline_model(model)

    return 1 - get_conditional_rvar(model)/get_conditional_rvar(baseline_model)


def get_R2(model: LimeTr):
    baseline_model = get_baseline_model(model)
    return {
        "conditional_R2": get_conditional_R2(model, baseline_model),
        "marginal_R2": get_marginal_R2(model, baseline_model)
    }


def get_rmse(model: LimeTr):
    return np.sqrt(get_marginal_rvar(model))
