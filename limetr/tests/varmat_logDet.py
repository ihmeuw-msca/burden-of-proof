# check utils logDet


def varmat_logDet():
    import numpy as np
    from limetr.utils import VarMat

    ok = True
    tol = 1e-10
    # setup problem
    # -------------------------------------------------------------------------
    mat = VarMat.testProblem()
    D = mat.varMat()

    tr_y = np.log(np.linalg.det(D))
    my_y = mat.logDet()

    err = np.linalg.norm(tr_y - my_y)
    ok = ok and err < tol

    if not ok:
        print('err in logDet')
        print('err:', err)

    return ok
