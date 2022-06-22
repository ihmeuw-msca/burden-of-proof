# check utils izdiag


def izmat_izdiag():
    import numpy as np
    from limetr.special_mat import izmat

    ok = True
    tol = 1e-10
    # setup problem
    # -------------------------------------------------------------------------
    k = 3
    n = np.array([5, 2, 4])
    m = n.size

    z_list = [np.random.randn(n[i], k) for i in range(m)]

    z = np.vstack(z_list)
    
    ns = np.minimum(n, k)
    nu = ns*n
    nx = n
    nz = n

    u = np.zeros(nu.sum())
    s = np.zeros(ns.sum())

    izmat.zdecomp(nz, nu, ns, z, u, s)
    my_y = izmat.izdiag(n.sum(), nu, ns, nx, u, s**2)

    y_list = [np.diag(np.eye(n[i]) + z_list[i].dot(z_list[i].T))
              for i in range(m)]

    tr_y = np.hstack(y_list)

    err = np.linalg.norm(tr_y - my_y)

    if not ok:
        print('err in izdiag')
        print('err:', err)

    return ok
