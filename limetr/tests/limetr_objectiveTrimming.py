# test function objectiveTrimming


def limetr_objectiveTrimming():
    import numpy as np
    from limetr.__init__ import LimeTr

    ok = True
    # setup test problem
    # -------------------------------------------------------------------------
    model = LimeTr.testProblem(use_trimming=True)

    tol = 1e-8

    # test objectiveTrimming
    # -------------------------------------------------------------------------
    x = np.hstack((model.beta, model.gamma))
    w = model.w

    r = model.Y - model.F(model.beta)
    t = (model.Z**2).dot(model.gamma)
    d = model.V + t

    tr_obj = 0.5*np.sum(r**2*w/d) + 0.5*model.N*np.log(2.0*np.pi)\
        + 0.5*w.dot(np.log(d))
    my_obj = model.objectiveTrimming(w)

    err = np.abs(tr_obj - my_obj)
    ok = ok and err < tol

    if not ok:
        print('err', err)
        print('tr_obj', tr_obj)
        print('my_obj', my_obj)

    return ok
