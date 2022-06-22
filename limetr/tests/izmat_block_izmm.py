# check utils block_izmm


def izmat_block_izmm():
    import numpy as np
    from limetr.special_mat import izmat

    ok = True
    tol = 1e-10
    # problem 1, tall matrix
    # -------------------------------------------------------------------------
    n, k = 6, 3
    l = min(n, k)

    z = np.random.randn(n, k)
    x = np.random.randn(n, 5)
    
    my_u = np.zeros(n*l)
    my_s = np.zeros(l)
    izmat.lsvd(z, my_u, my_s)

    tr_y = x + z.dot(z.T.dot(x))
    my_y = np.zeros((n, 5), order='F')
    izmat.block_izmm(my_u, my_s**2, x, my_y)

    err = np.linalg.norm(tr_y - my_y)
    ok = ok and err < tol

    if not ok:
        print('err in block_izmm tall matrix')
        print('err:', err)

    # problem 2, fat matrix
    # -------------------------------------------------------------------------
    n, k = 3, 6
    l = min(n, k)

    z = np.random.randn(n, k)
    x = np.random.randn(n, 5)
    
    my_u = np.zeros(n*l)
    my_s = np.zeros(l)
    izmat.lsvd(z, my_u, my_s)

    tr_y = x + z.dot(z.T.dot(x))
    my_y = np.zeros((n, 5), order='F')
    izmat.block_izmm(my_u, my_s**2, x, my_y)

    err = np.linalg.norm(tr_y - my_y)
    ok = ok and err < tol

    if not ok:
        print('err in block_izmm fat matrix')
        print('err:', err)

    return ok
