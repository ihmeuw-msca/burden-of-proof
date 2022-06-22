# check utils dot


def varmat_invDot():
    import numpy as np
    from limetr.utils import VarMat

    ok = True
    tol = 1e-10
    # setup problem
    # -------------------------------------------------------------------------
    mat = VarMat.testProblem()
    inv_D = mat.invVarMat()
    x = np.random.randn(mat.N)
    X = np.random.randn(mat.N, 5)

    tr_y = inv_D.dot(x)
    tr_Y = inv_D.dot(X)

    my_y = mat.invDot(x)
    my_Y = mat.invDot(X)

    err = np.linalg.norm(tr_y - my_y) + np.linalg.norm(tr_Y - my_Y)
    ok = ok and err < tol

    if not ok:
        print('err in invDot')
        print('err:', err)

    return ok
