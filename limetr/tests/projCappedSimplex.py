# test function dot


def projCappedSimplex():
    import numpy as np
    from limetr.utils import projCappedSimplex

    ok = True
    # setup test problem
    # -------------------------------------------------------------------------
    w = np.ones(10)
    sum_w = 9.0

    tr_w = np.repeat(0.9, 10)
    my_w = projCappedSimplex(w, sum_w)

    tol = 1e-10
    err = np.linalg.norm(tr_w - my_w)

    ok = ok and err < tol

    if not ok:
        print('tr_w', tr_w)
        print('my_w', my_w)

    return ok
