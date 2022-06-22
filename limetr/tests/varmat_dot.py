# check utils dot


def varmat_dot():
    import numpy as np
    from limetr.utils import VarMat

    ok = True
    tol = 1e-10
    # setup problem
    # -------------------------------------------------------------------------
    mat = VarMat.testProblem()
    D = mat.varMat()
    x = np.random.randn(mat.N)
    X = np.random.randn(mat.N, 5)

    tr_y = D.dot(x)
    tr_Y = D.dot(X)

    my_y = mat.dot(x)
    my_Y = mat.dot(X)

    err = np.linalg.norm(tr_y - my_y) + np.linalg.norm(tr_Y - my_Y)
    ok = ok and err < tol

    if not ok:
        print('err in dot')
        print('err:', err)

    return ok
