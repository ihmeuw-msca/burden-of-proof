# test function lprior


def limetr_lasso():
    import numpy as np
    from limetr.__init__ import LimeTr

    ok = True
    # setup test problem
    # -------------------------------------------------------------------------
    model = LimeTr.testProblemLasso()

    tol = 1e-6

    # test lasso
    # -------------------------------------------------------------------------
    model.optimize()
    beta = model.beta
    zero_idx = np.abs(beta) <= 1e-8
    beta[zero_idx] = 0.0

    # calculate the gradient
    g_beta = -model.JF(beta).T.dot(model.Y - model.F(beta))
    for i in range(model.k_beta):
        if beta[i] == 0.0 and np.abs(g_beta[i]) < model.lw[i]:
            g_beta[i] = 0.0
        else:
            g_beta[i] += np.sign(beta[i])*model.lw[i]

    err = np.linalg.norm(g_beta)
    ok = ok and err < tol

    if not ok:
        print('err', err)

    return ok
