module source_terms_module
    implicit none
    contains

!---------------------------------------------------------------------
!Source steps for EDC (mbHSGN) and BBBD (HSGN).
!The Riemann solver (rp1_dispersive_SWEs_fwave_splitting.f90) only contains the hydrostatic
!bathymetry term and -c^2*ubar*\partial_x(eta). The non-hydrostatic bathymetry terms are integrated here,
!implicitly and together with the relaxation terms:
!one SDIRK(2,2) step of a local linear system per cell, with h frozen and b_x = cell-centred
!bottom slope (source_terms.py: cell_centred_bx). The local systems are neutral for any b_x, so
!the step is stable for any dt and any (under-resolved) bathymetry jump.
!---------------------------------------------------------------------

subroutine src_EDC_1D(q,bx,qout,dt,csq,gamma,mx)
    !Integrates, cell by cell with h frozen, the linear system for (u,w,p)
    !     u_t = -gamma b_x p / h
    !     w_t =  gamma     p / h
    !     p_t =  c^2 b_x u / h - 2 c^2 w / h
    !i.e. the relaxation source plus the two non-hydrostatic bathymetry terms of mbHSGN with SDIRK(2,2).
    implicit none
    integer, intent(in) :: mx
    double precision, intent(in) :: dt
    double precision, intent(in) :: q(4,mx), bx(mx)
    double precision, intent(in) :: gamma, csq
    double precision, intent(out) :: qout(4,mx)

    !Locals
    integer :: i
    double precision :: h,hu,hw,hp

    do i=1,mx
        !Copy input to output
        qout(:,i) = q(:,i)
        h = q(1,i)
        hu = q(2,i)
        hw = q(3,i)
        hp = q(4,i)
        call one_step_SDIRK22_EDC(h,hu,hw,hp,bx(i),csq,gamma,dt)
        !Update output
        qout(2,i) = hu
        qout(3,i) = hw
        qout(4,i) = hp
    end do

end subroutine src_EDC_1D

subroutine src_EDC_1D_radial(q,r,bx,qout,dt,csq,gamma,mx)
    !As src_EDC_1D (bx = d(b)/dr), followed by the explicit radial terms
    implicit none
    integer, intent(in) :: mx
    double precision, intent(in) :: dt
    double precision, intent(in) :: q(4,mx)
    double precision, intent(in) :: r(mx), bx(mx)
    double precision, intent(in) :: gamma, csq
    double precision, intent(out) :: qout(4,mx)

    !Locals
    integer :: i
    double precision :: h,hu,hw,hp
    double precision ::ratio,u

    do i=1,mx
        !Copy input to output
        qout(:,i) = q(:,i)
        h = q(1,i)
        hu = q(2,i)
        hw = q(3,i)
        hp = q(4,i)
        call one_step_SDIRK22_EDC(h,hu,hw,hp,bx(i),csq,gamma,dt)
        !Radial terms
        ratio = dt/r(i)
        qout(1,i) = h-ratio*q(2,i)
        u = q(2,i)/h
        qout(2,i) = hu-ratio*(h*u**2)
        qout(3,i) = hw - ratio*u*hw
        qout(4,i) = hp - ratio*(u*hp+h*u*csq)
    end do
end subroutine src_EDC_1D_radial

subroutine src_BBBD_1D(q,bx,qout,dt,csq,mx)
    !Integrates, cell by cell with h frozen, the linear system for y = (u,w,sg,p,pb)
    !     u_t  = -b_x pb / h
    !     w_t  =  pb / h
    !     sg_t = (12 p - 6 pb) / h
    !     p_t  = (-c^2 sg - c^2 b_x u) / h
    !     pb_t = (6 c^2 b_x u - 6 c^2 w + 3 c^2 sg) / h
    !i.e. the relaxation source plus the bathymetry terms that are not in the Riemann solver  with SDIRK(2,2).
    implicit none
    integer, intent(in) :: mx
    double precision, intent(in) :: dt
    double precision, intent(in) :: q(6,mx), bx(mx)
    double precision, intent(in) :: csq
    double precision, intent(out) :: qout(6,mx)

    !Locals
    integer :: i
    double precision :: h,hu,hw,hsg,hp,hpb

    do i=1,mx
        !Copy input to output
        qout(:,i) = q(:,i)
        h = q(1,i)
        hu = q(2,i)
        hw = q(3,i)
        hsg = q(4,i)
        hp = q(5,i)
        hpb = q(6,i)
        call one_step_SDIRK22_BBBD(h,hu,hw,hsg,hp,hpb,bx(i),csq,dt)
        !Update output
        qout(2,i) = hu
        qout(3,i) = hw
        qout(4,i) = hsg
        qout(5,i) = hp
        qout(6,i) = hpb
    end do

end subroutine src_BBBD_1D

subroutine src_BBBD_1D_radial(q,r,bx,qout,dt,csq,mx)
    !As src_BBBD_1D (bx = d(b)/dr), followed by the explicit radial terms
    implicit none
    integer, intent(in) :: mx
    double precision, intent(in) :: dt
    double precision, intent(in) :: q(6,mx)
    double precision, intent(in) :: r(mx), bx(mx)
    double precision, intent(in) :: csq
    double precision, intent(out) :: qout(6,mx)

    !Locals
    integer :: i
    double precision :: h,hu,hw,hsg,hp,hpb
    double precision :: ratio,u

    do i=1,mx
        !Copy input to output
        qout(:,i) = q(:,i)
        h = q(1,i)
        hu = q(2,i)
        hw = q(3,i)
        hsg = q(4,i)
        hp = q(5,i)
        hpb = q(6,i)
        call one_step_SDIRK22_BBBD(h,hu,hw,hsg,hp,hpb,bx(i),csq,dt)
        !Radial terms
        ratio = dt/r(i)
        qout(1,i) = h-ratio*q(2,i)
        u = q(2,i)/h
        qout(2,i) = hu-ratio*(h*u**2)
        qout(3,i) = hw - ratio*u*hw
        qout(4,i) = hsg - ratio*u*hsg
        qout(5,i) = hp - ratio*(u*hp+h*u*csq)
        qout(6,i) = hpb - ratio*u*hpb
    end do

end subroutine src_BBBD_1D_radial

subroutine src_DR_1D(q,qout,dt,csq,mx)
    implicit none
    integer, intent(in) :: mx
    double precision, intent(in) :: dt
    double precision, intent(in) :: q(5,mx)
    double precision, intent(in) :: csq
    double precision, intent(out) :: qout(5,mx)

    !Locals
    integer :: i
    double precision :: h,hw,hp

    do i=1,mx
        !Copy input to output
        qout(:,i) = q(:,i)
        h = q(1,i)
        hw = q(3,i)
        hp = q(4,i)
        !The linear source term is the same as EDC's system
        call one_exact_step_EDC(h,hw,hp,csq,1.5d0,dt)
        !Update output
        qout(3,i) = hw
        qout(4,i) = hp
    end do

end subroutine src_DR_1D

!---------------------------------------------------------------------
!Exact integration of the relaxation terms (flat bottom; used by the DR system only)
!---------------------------------------------------------------------
subroutine one_exact_step_EDC(h, hw, hp, csq_space, gammaEDC,dt)
    implicit none
    double precision, intent(in) :: h, csq_space, gammaEDC, dt
    double precision, intent(inout) :: hw, hp
    !Locals
    double precision :: p,w,expdtA(2,2),result(2,1)
    double precision :: b, posc
    w = hw/h
    p = hp/h
    posc = dt*2.d0*csq_space/h
    b = dt*gammaEDC/h
    !Calculate the exact solution for the half problem:
    !(w,p)_t = h^{-1}(gamma*p,-2*csq*w) = A(p,w),
    !where A= [[0,gamma],[-2*csq,0]]/h.
    !The exact solution is given by:
    !(w,p)(dt) = exp(dt*A)*(w,p)(0).
    expdtA(1,1) = cos(dsqrt(b*posc)) !cosh(i*a) = cos(a)
    expdtA(1,2) = dsqrt(b/posc)*sin(dsqrt(b*posc)) !sinh(i*a) = i*sin(a)
    expdta(2,1) = -dsqrt(posc/b)*sin(dsqrt(b*posc)) !sin(i*a) = i*sinh(a)
    expdtA(2,2) = cos(dsqrt(b*posc)) !cosh(i*a) = cos(a)
    result(:,1) = (/w,p/)
    result = matmul(expdtA,result)
    !Compute w,p
    ! w = expdtA(1,1)*w + expdtA(1,2)*p
    ! p = expdtA(2,1)*w + expdtA(2,2)*p
    !Output
    hw = result(1,1)*h
    hp = result(2,1)*h
end subroutine one_exact_step_EDC

subroutine one_step_SDIRK22_EDC(h, hu, hw, hp, bx, csq_space, gammaEDC, dt)
    !SDIRK(2,2) step for y=(u,w,p), y' = A y with A = [[0,0,a],[0,0,b],[c,d,0]],
    !a = -gamma bx/h, b = gamma/h, c = c^2 bx/h, d = -2 c^2/h  (see src_EDC_1D).
    !(I - G dt A) y = r is solved in closed form (det = 1 - (G dt)^2 (ac+bd) >= 1).
    implicit none
    double precision, intent(in) :: h, bx, csq_space, gammaEDC, dt
    double precision, intent(inout) :: hu, hw, hp
    double precision :: a, b, c, d, G, k, den
    double precision :: u, w, p, u1, w1, p1, u2, w2, p2

    if (h<1.e-3 .or. csq_space < 1.e-3) then
        !Do nothing
        return
    end if
    u = hu/h
    w = hw/h
    p = hp/h
    a = -gammaEDC*bx/h
    b = gammaEDC/h
    c = csq_space*bx/h
    d = -2.d0*csq_space/h
    G = (1.d0 - 1.d0/dsqrt(2.d0))*dt
    k = dt - G
    den = 1.d0 - G*G*(a*c + b*d)
    !Stage 1: (I - G A) y1 = y
    p1 = (p + G*c*u + G*d*w)/den
    u1 = u + G*a*p1
    w1 = w + G*b*p1
    !Stage 2: (I - G A) y2 = y + (dt-G) A y1
    u = u + k*a*p1
    w = w + k*b*p1
    p = p + k*(c*u1 + d*w1)
    p2 = (p + G*c*u + G*d*w)/den
    u2 = u + G*a*p2
    w2 = w + G*b*p2
    hu = u2*h
    hw = w2*h
    hp = p2*h
end subroutine one_step_SDIRK22_EDC

subroutine one_step_SDIRK22_BBBD(h, hu, hw, hsg, hp, hpb, bx, csq_space, dt)
    !SDIRK(2,2) step for y = (u,w,sg,p,pb), y' = A y (see src_BBBD_1D).
    !With bx = 0 it is the SDIRK(2,2) step of the relaxation terms alone.
    implicit none
    double precision, intent(in) :: h, bx, csq_space, dt
    double precision, intent(inout) :: hu, hw, hsg, hp, hpb
    !Matrix handling
    double precision, dimension(5,5) :: A
    double precision, dimension(5) :: qn, qnp1

    if (h<1.e-3 .or. csq_space < 1.e-3) then
        !Do nothing
        return
    end if
    qn = (/hu/h, hw/h, hsg/h, hp/h, hpb/h/)
    A = 0.d0
    A(1,5) = -bx/h
    A(2,5) = 1.d0/h
    A(3,4) = 12.d0/h
    A(3,5) = -6.d0/h
    A(4,1) = -csq_space*bx/h
    A(4,3) = -csq_space/h
    A(5,1) = 6.d0*csq_space*bx/h
    A(5,2) = -6.d0*csq_space/h
    A(5,3) = 3.d0*csq_space/h

    call one_general_sdirk2_step(dt,A,qn,qnp1)
    hu  = qnp1(1)*h
    hw  = qnp1(2)*h
    hsg = qnp1(3)*h
    hp  = qnp1(4)*h
    hpb = qnp1(5)*h
end subroutine one_step_SDIRK22_BBBD

subroutine one_general_sdirk2_step(dt, A, qn, qnp1)
    implicit none
    ! Single-call SDIRK2 step for dq/dt = A*q
    ! All arrays are DOUBLE PRECISION (fortran default double)
    !
    ! Inputs:
    !   dt      - DOUBLE PRECISION scalar timestep
    !   A(n,n)  - DOUBLE PRECISION system matrix
    !   qn(n)  - DOUBLE PRECISION current state
    !
    ! Output:
    !   qnp1(n) - DOUBLE PRECISION next state
    !
    ! Method:
    !   gamma = 1 - 1/sqrt(2)
    !   M = I - gamma*dt*A
    !   factorize M (LU)
    !   solve M * q1 = qn
    !   form rhs2 = qn + dt*(1-gamma)*A*q1
    !   solve M * q2 = rhs2
    !   qnp1 = q2
    !
    ! Requires LAPACK: dgetrf, dgetrs
    !
    double precision, intent(in)  :: dt
    double precision, intent(in)  :: A(:,:)
    double precision, intent(in)  :: qn(:)
    double precision, intent(out) :: qnp1(:)

    integer :: n, info, lda, ldb, nrhs
    integer, allocatable :: ipiv(:)
    double precision, allocatable :: M(:,:), rhs(:), q1(:), q2(:)
    double precision :: gamma
    integer :: i

    ! For LAPACK external declarations 
    external dgetrf, dgetrs

    n = size(qn)
    if (size(A,1) /= n .or. size(A,2) /= n .or. size(qnp1) /= n) then
        stop "Dimension mismatch in sdirk2_step"
    end if

    gamma = 1.0d0 - 1.0d0/sqrt(2.0d0)

    lda = n
    ldb = n
    nrhs = 1

    allocate(M(n,n), rhs(n), q1(n), q2(n), ipiv(n))

    ! Build M = I - gamma*dt*A
    M = -gamma * dt * A
    do i = 1, size(M,1)
        M(i,i) = M(i,i) + 1.0d0
    end do

    ! LU factorization of M (overwrites M with LU factors)
    call dgetrf(n, n, M, lda, ipiv, info)
    if (info /= 0) then
        stop "dgetrf failed in one_general_sdirk2_step (singular or illegal argument)"
    end if

    ! Stage 1: solve M * q1 = qn  (dgetrs uses LU factors in M and ipiv)
    rhs = qn
    call dgetrs('N', n, nrhs, M, lda, ipiv, rhs, ldb, info)
    if (info /= 0) then
        stop "dgetrs failed at stage 1 in one_general_sdirk2_step"
    end if
    q1 = rhs

    ! Build RHS for stage 2: qn + dt*(1-gamma)*A*q1
    rhs = qn + dt * (1.0d0 - gamma) * matmul(A, q1)

    ! Stage 2: solve M * q2 = rhs  (reuse LU factors in M and ipiv)
    call dgetrs('N', n, nrhs, M, lda, ipiv, rhs, ldb, info)
    if (info /= 0) then
        stop "dgetrs failed at stage 2 in one_general_sdirk2_step"
    end if
    q2 = rhs

    qnp1 = q2

    deallocate(M, rhs, q1, q2, ipiv)
end subroutine one_general_sdirk2_step

end module source_terms_module