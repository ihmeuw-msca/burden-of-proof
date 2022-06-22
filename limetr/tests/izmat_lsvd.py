# check utils lsvd


def izmat_lsvd():
    import numpy as np
    from limetr.special_mat import izmat

    ok = True
    tol = 1e-10
    # problem 1, tall matrix
    # -------------------------------------------------------------------------
    n, k = 6, 3
    z = np.random.randn(n, k)
    tr_u, tr_s, tr_vt = np.linalg.svd(z, full_matrices=False)
    my_u = np.zeros(tr_u.size)
    my_s = np.zeros(tr_s.size)
    izmat.lsvd(z, my_u, my_s)

    err = np.linalg.norm(my_u.reshape(k, n).T - tr_u)
    ok = ok and err < tol

    if not ok:
        print('err in lsvd tall matrix')
        print('err:', err)

    # problem 2, fat matrix
    # -------------------------------------------------------------------------
    n, k = 3, 6
    z = np.random.randn(n, k)
    tr_u, tr_s, tr_vt = np.linalg.svd(z, full_matrices=False)
    my_u = np.zeros(tr_u.size)
    my_s = np.zeros(tr_s.size)
    izmat.lsvd(z, my_u, my_s)

    err = np.linalg.norm(np.abs(my_u.reshape(n, n).T) - np.abs(tr_u))
    ok = ok and err < tol

    if not ok:
        print('err in lsvd fat matrix')
        print('err:', err)

    return ok
