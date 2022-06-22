# check utils zdecomp


def izmat_zdecomp():
    import numpy as np
    from limetr.special_mat import izmat

    ok = True
    tol = 1e-10
    # setup problem
    # -------------------------------------------------------------------------
    k = 3
    n = [5, 2, 4]

    z_list = []
    tr_u_list = []
    tr_s_list = []
    for i in range(len(n)):
        z_list.append(np.random.randn(n[i], k))
        u, s, vt = np.linalg.svd(z_list[-1], full_matrices=False)
        tr_u_list.append(u)
        tr_s_list.append(s)

    z = np.vstack(z_list)
    tr_u = np.hstack([u.reshape(u.size, order='F') for u in tr_u_list])
    tr_s = np.hstack(tr_s_list)

    my_u = np.zeros(tr_u.size)
    my_s = np.zeros(tr_s.size)

    nz = [z_sub.shape[0] for z_sub in z_list]
    nu = [u_sub.size for u_sub in tr_u_list]
    ns = [s_sub.size for s_sub in tr_s_list]

    izmat.zdecomp(nz, nu, ns, z, my_u, my_s)


    if not ok:
        print('err in zdecomp')
        print('err:', err)

    return ok
