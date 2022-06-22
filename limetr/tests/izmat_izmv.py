# check utils izmv


def izmat_izmv():
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
    x_list = [np.random.randn(n[i]) for i in range(m)]

    z = np.vstack(z_list)
    x = np.hstack(x_list)
    
    ns = np.minimum(n, k)
    nu = ns*n
    nx = n
    nz = n

    u = np.zeros(nu.sum())
    s = np.zeros(ns.sum())

    izmat.zdecomp(nz, nu, ns, z, u, s)
    my_y = izmat.izmv(nu, ns, nx, u, s**2, x)

    y_list = [x_list[i] + z_list[i].dot(z_list[i].T.dot(x_list[i]))
              for i in range(m)]

    tr_y = np.hstack(y_list)

    err = np.linalg.norm(tr_y - my_y)

    if not ok:
        print('err in izmv')
        print('err:', err)

    return ok
