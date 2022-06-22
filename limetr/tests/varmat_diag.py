# check utils diag


def varmat_diag():
    import numpy as np
    from limetr.utils import VarMat

    ok = True
    tol = 1e-10
    # setup problem
    # -------------------------------------------------------------------------
    mat = VarMat.testProblem()
    D = mat.varMat()

    tr_y = np.diag(D)

    my_y = mat.diag()

    err = np.linalg.norm(tr_y - my_y)
    ok = ok and err < tol

    if not ok:
        print('err in diag')
        print('err:', err)

    return ok
