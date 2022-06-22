# test function objective


def limetr_objective():
    import numpy as np
    from limetr.__init__ import LimeTr

    ok = True
    # setup test problem
    # -------------------------------------------------------------------------
    model = LimeTr.testProblem(use_constraints=True,
                               use_regularizer=True,
                               use_uprior=True,
                               use_gprior=True,
                               know_obs_std=False)

    tol = 1e-8

    # test objective
    # -------------------------------------------------------------------------
    x = np.random.randn(model.k)
    x[model.idx_gamma] = 0.1
    x[model.idx_delta] = 0.1

    tr_obj = model.objective(x, use_ad=True)
    my_obj = model.objective(x)

    err = np.abs(tr_obj - my_obj)
    ok = ok and err < tol

    if not ok:
        print('err', err)
        print('tr_obj', tr_obj)
        print('my_obj', my_obj)

    return ok
