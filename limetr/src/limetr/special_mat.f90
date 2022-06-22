! fortran helper functions for special mat calculation


MODULE IZMAT


CONTAINS


! LSVD compute the svd of a matrix and output the singular value and the left
! singular vector
! -----------------------------------------------------------------------------
SUBROUTINE LSVD(sz1, sz2, su, ss, z, u, s)
    IMPLICIT NONE
    INTEGER(8) :: sz1, sz2, su, ss
    REAL(8) :: z(sz1, sz2), u(su), s(ss)
    INTEGER(8) :: lwork, info
    REAL(8), allocatable :: work(:)

    !f2py INTENT(IN) :: z
    !f2py INTENT(HIDE), DEPEND(z) :: sz1 = SIZE(z, 1), sz2 = SIZE(z, 2)
    !f2py INTENT(HIDE), DEPEND(u) :: su = SIZE(u)
    !f2py INTENT(HIDE), DEPEND(s) :: sp = SIZE(s)
    !f2py INTENT(INOUT) :: u, s

    lwork = 2*MAX(1, 3*ss + MAX(sz1, sz2), 5*ss)

    ALLOCATE(work(lwork))
    CALL DGESVD('S', 'N', sz1, sz2, z, sz1, s, u, sz1, z, ss, &
                work, lwork, info)
    DEALLOCATE(work)

END SUBROUTINE LSVD


! ZDECOMP compute the svd of the block of Z matrix
! -----------------------------------------------------------------------------
SUBROUTINE ZDECOMP(sz1, sz2, su, ss, m, nz, nu, ns, z, u, s)
    IMPLICIT NONE
    INTEGER(8) :: sz1, sz2, su, ss, m
    INTEGER(8) :: nz(m), nu(m), ns(m)
    REAL(8) :: z(sz1, sz2), u(su), s(ss)
    INTEGER(8) :: az, bz, au, bu, as, bs, i

    !f2py INTENT(IN) :: nz, nu, ns, z
    !f2py INTENT(HIDE), DEPEND(nz) :: m = SIZE(nz)
    !f2py INTENT(HIDE), DEPEND(z) :: sz1 = SIZE(z, 1), sz2 = SIZE(z, 2)
    !f2py INTENT(HIDE), DEPEND(u) :: su = SIZE(u)
    !f2py INTENT(HIDE), DEPEND(s) :: ss = SIZE(s)
    !f2py INTENT(INOUT) :: u, s

    az = 1
    au = 1
    as = 1
    bz = 0
    bu = 0
    bs = 0

    DO i = 1, m
        bz = bz + nz(i)
        bu = bu + nu(i)
        bs = bs + ns(i)

        CALL LSVD(nz(i), sz2, nu(i), ns(i), z(az:bz, :), u(au:bu), s(as:bs))

        az = az + nz(i)
        au = au + nu(i)
        as = as + ns(i)
    END DO

END SUBROUTINE ZDECOMP


! BLOCK_IZMV calculate matrix vector multiplication of block izmat
! -----------------------------------------------------------------------------
SUBROUTINE BLOCK_IZMV(su, sd, sx, u, d, x, y)
    IMPLICIT NONE
    INTEGER(8) :: su, sd, sx
    REAL(8) :: u(su), d(sd), x(sx), y(sx)
    REAL(8) :: c(sd)

    !f2py INTENT(IN) :: u, d, x
    !f2py INTENT(HIDE), DEPEND(u) :: su = SIZE(u)
    !f2py INTENT(HIDE), DEPEND(d) :: sd = SIZE(d)
    !f2py INTENT(HIDE), DEPEND(x) :: sx = SIZE(x)
    !f2py INTENT(INOUT) :: y

    CALL DGEMV('T', sx, sd, 1.d0, u, sx, x, 1, 0.d0, c, 1)
    c = c * d
    CALL DGEMV('N', sx, sd, 1.d0, u, sx, c, 1, 0.d0, y, 1)
    y = y + x

END SUBROUTINE BLOCK_IZMV


! IZMV calculate matrix vector multiplication of izmat
! -----------------------------------------------------------------------------
FUNCTION IZMV(su, sd, sx, m, nu, nd, nx, u, d, x) RESULT(y)
    IMPLICIT NONE
    INTEGER(8) :: su, sd, sx, m
    INTEGER(8) :: nu(m), nd(m), nx(m)
    REAL(8) :: u(su), d(sd), x(sx)
    REAL(8) :: y(sx)
    INTEGER(8) :: ax, bx, au, bu, ad, bd, i

    !f2py INTENT(IN) :: u, d, x
    !f2py INTENT(HIDE), DEPEND(u) :: su = SIZE(u)
    !f2py INTENT(HIDE), DEPEND(d) :: sd = SIZE(d)
    !f2py INTENT(HIDE), DEPEND(x) :: sx = SIZE(x)
    !f2py INTENT(HIDE), DEPEND(nx) :: m = SIZE(nx)

    ax = 1
    au = 1
    ad = 1
    bx = 0
    bu = 0
    bd = 0

    DO i = 1, m
        bx = bx + nx(i)
        bu = bu + nu(i)
        bd = bd + nd(i)

        CALL BLOCK_IZMV(nu(i), nd(i), nx(i), u(au:bu), d(ad:bd), x(ax:bx), &
                        y(ax:bx))

        ax = ax + nx(i)
        au = au + nu(i)
        ad = ad + nd(i)
    END DO

END FUNCTION IZMV


! BLOCK_IZMM calculate matrix matrix multiplication of block izmat
! -----------------------------------------------------------------------------
SUBROUTINE BLOCK_IZMM(su, sd, sx1, sx2, u, d, x, y)
    IMPLICIT NONE
    INTEGER(8) :: su, sd, sx1, sx2
    REAL(8) :: u(su), d(sd), x(sx1, sx2), y(sx1, sx2)
    INTEGER(8) :: i
    REAL(8) :: c(sd, sx2)

    !f2py INTENT(IN) :: u, d, x
    !f2py INTENT(HIDE), DEPEND(u) :: su = SIZE(u)
    !f2py INTENT(HIDE), DEPEND(d) :: sd = SIZE(d)
    !f2py INTENT(HIDE), DEPEND(x) :: sx1 = SIZE(x, 1)
    !f2py INTENT(HIDE), DEPEND(x) :: sx2 = SIZE(x, 2)
    !f2py INTENT(INOUT) :: y

    CALL DGEMM('T', 'N', sd, sx2, sx1, 1.d0, u, sx1, x, sx1, 0.d0, c, sd)
    DO i = 1, sx2
        c(:, i) = c(:, i) * d
    END DO
    CALL DGEMM('N', 'N', sx1, sx2, sd, 1.d0, u, sx1, c, sd, 0.d0, y, sx1)
    y = y + x

END SUBROUTINE BLOCK_IZMM


! IZMM calculate matrix matrix multiplication of izmat
! -----------------------------------------------------------------------------
FUNCTION IZMM(su, sd, sx1, sx2, m, nu, nd, nx, u, d, x) RESULT(y)
    IMPLICIT NONE
    INTEGER(8) :: su, sd, sx1, sx2, m
    INTEGER(8) :: nu(m), nd(m), nx(m)
    REAL(8) :: u(su), d(sd), x(sx1, sx2)
    REAL(8) :: y(sx1, sx2)
    INTEGER(8) :: ax, bx, au, bu, ad, bd, i

    !f2py INTENT(IN) :: u, d, x
    !f2py INTENT(HIDE), DEPEND(u) :: su = SIZE(u)
    !f2py INTENT(HIDE), DEPEND(d) :: sd = SIZE(d)
    !f2py INTENT(HIDE), DEPEND(x) :: sx1 = SIZE(x, 1)
    !f2py INTENT(HIDE), DEPEND(x) :: sx2 = SIZE(x, 2)
    !f2py INTENT(HIDE), DEPEND(nx) :: m = SIZE(nx)

    ax = 1
    au = 1
    ad = 1
    bx = 0
    bu = 0
    bd = 0

    DO i = 1, m
        bx = bx + nx(i)
        bu = bu + nu(i)
        bd = bd + nd(i)

        CALL BLOCK_IZMM(nu(i), nd(i), nx(i), sx2, u(au:bu), d(ad:bd), &
                        x(ax:bx, :), y(ax:bx, :))

        ax = ax + nx(i)
        au = au + nu(i)
        ad = ad + nd(i)
    END DO

END FUNCTION IZMM


! IZEIG calculate the eigenvalues of the matrix in the iz form
! -----------------------------------------------------------------------------
FUNCTION IZEIG(sd, m, sz, nz, nd, d) RESULT(v)
    IMPLICIT NONE
    INTEGER(8) :: sd, m, sz
    INTEGER(8) :: nz(m), nd(m)
    INTEGER(8) :: i, ad, bd, av, bv
    REAL(8) :: d(sd), v(sz)

    !f2py INTENT(IN) :: sz, nz, nd, d
    !f2py INTENT(HIDE), DEPEND(d) :: sd = SIZE(d)
    !f2py INTENT(HIDE), DEPEND(nd) :: m = SIZE(nd)

    ad = 1
    bd = 0
    av = 1
    bv = 0

    v = 1.d0

    DO i = 1, m
        bd = bd + nd(i)
        bv = av + nd(i) - 1

        v(av:bv) = v(av:bv) + d(ad:bd)

        ad = ad + nd(i)
        av = av + nz(i)
    END DO

END FUNCTION IZEIG


! BLOCK_IZDIAG return the diagonal of the izmat block
! -----------------------------------------------------------------------------
SUBROUTINE BLOCK_IZDIAG(su, sd, sx, u, d, x)
    IMPLICIT NONE
    INTEGER(8) :: su, sd, sx
    REAL(8) :: u(su), d(sd), x(sx)

    !f2py INTENT(IN) :: u, d
    !f2py INTENT(HIDE), DEPEND(u) :: su = SIZE(u)
    !f2py INTENT(HIDE), DEPEND(d) :: sd = SIZE(d)
    !f2py INTENT(HIDE), DEPEND(x) :: sx = SIZE(x)
    !f2py INTENT(INOUT) :: x

    CALL DGEMV('N', sx, sd, 1.d0, u**2, sx, d, 1, 0.d0, x, 1)

    x = x + 1.d0

END SUBROUTINE BLOCK_IZDIAG


! IZDIAG return the diagonal of the izmat
! -----------------------------------------------------------------------------
FUNCTION IZDIAG(su, sd, sx, m, nu, nd, nx, u, d) RESULT(x)
    IMPLICIT NONE
    INTEGER(8) :: su, sd, sx, m
    INTEGER(8) :: nu(m), nd(m), nx(m)
    REAL(8) :: u(su), d(sd), x(sx)
    INTEGER(8) :: au, bu, ad, bd, ax, bx, i

    !f2py INTENT(IN) :: sx, nu, nd, nx, u, d
    !f2py INTENT(HIDE), DEPEND(u) :: su = SIZE(u)
    !f2py INTENT(HIDE), DEPEND(d) :: sd = SIZE(d)
    !f2py INTENT(HIDE), DEPEND(nd) :: m = SIZE(nd)

    au = 1
    bu = 0
    ad = 1
    bd = 0
    ax = 1
    bx = 0

    DO i = 1, m
        bu = bu + nu(i)
        bd = bd + nd(i)
        bx = bx + nx(i)

        CALL BLOCK_IZDIAG(nu(i), nd(i), nx(i), u(au:bu), d(ad:bd), x(ax:bx))

        au = au + nu(i)
        ad = ad + nd(i)
        ax = ax + nx(i)
    END DO

END FUNCTION IZDIAG


END MODULE IZMAT
