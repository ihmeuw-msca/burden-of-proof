# check utils invDiag


def varmat_invDiag():
    import numpy as np
    from limetr.utils import VarMat

    ok = True
    tol = 1e-10
    # setup problem
    # -------------------------------------------------------------------------
    mat = VarMat.testProblem()
    inv_D = mat.invVarMat()

    tr_y = np.diag(inv_D)

    my_y = mat.invDiag()

    err = np.linalg.norm(tr_y - my_y)
    ok = ok and err < tol

    if not ok:
        print('err in invDiag')
        print('err:', err)

    return ok
