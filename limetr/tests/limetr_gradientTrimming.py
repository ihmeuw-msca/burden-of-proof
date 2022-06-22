# test function gradientTrimming


def limetr_gradientTrimming():
    import numpy as np
    from limetr.__init__ import LimeTr

    ok = True
    # setup test problem
    # -------------------------------------------------------------------------
    model = LimeTr.testProblem(use_trimming=True)

    # decouple all the studies
    model.n = np.array([1]*model.N)

    tol = 1e-8

    # test gradientTrimming
    # -------------------------------------------------------------------------
    x = np.hstack((model.beta, model.gamma))
    w = model.w

    tr_grad = model.gradientTrimming(w, use_ad=True)
    my_grad = model.gradientTrimming(w)

    err = np.linalg.norm(tr_grad - my_grad)
    ok = ok and err < tol

    if not ok:
        print('err', err)
        print('tr_grad', tr_grad)
        print('my_grad', my_grad)

    return ok
