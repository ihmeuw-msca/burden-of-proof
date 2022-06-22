# check utils izeig


def izmat_izeig():
    import numpy as np
    from scipy.linalg import block_diag
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
    
    my_eig = izmat.izeig(sum(n), n, ns, s**2)
    tr_eig, vec = np.linalg.eig(block_diag(*[
            np.eye(n[i]) + z_list[i].dot(z_list[i].T)
            for i in range(len(n))
        ]))

    err = np.linalg.norm(tr_eig - my_eig)

    if not ok:
        print('err in izeig')
        print('err:', err)

    return ok
