# test function gradient


def limetr_gradient():
    import numpy as np
    from limetr.__init__ import LimeTr

    ok = True
    # setup test problem
    # -------------------------------------------------------------------------
    model = LimeTr.testProblem(use_trimming=True,
                               use_constraints=True,
                               use_regularizer=True,
                               use_uprior=True,
                               use_gprior=True,
                               know_obs_std=False,
                               share_obs_std=True)

    tol = 1e-6

    # test the gradient
    # -------------------------------------------------------------------------
    x = np.random.randn(model.k)
    x[model.idx_gamma] = 0.1
    x[model.idx_delta] = 0.1

    tr_grad = model.gradient(x, use_ad=True)
    my_grad = model.gradient(x)

    err = np.linalg.norm(tr_grad - my_grad)
    ok = ok and err < tol

    if not ok:
        print('err', err)
        print('tr_grad', tr_grad)
        print('my_grad', my_grad)

    return ok
