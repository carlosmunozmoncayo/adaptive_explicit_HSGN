module source_terms_module
    implicit none
    contains

subroutine src_EDC_1D(q,qout,dt,csq,gamma,mx)
    implicit none
    integer, intent(in) :: mx
    double precision, intent(in) :: dt
    double precision, intent(in) :: q(4,mx)
    double precision, intent(in) :: gamma, csq
    double precision, intent(out) :: qout(4,mx)

    !Locals
    integer :: i
    double precision :: h,hw,hp

    do i=1,mx
        !Copy input to output
        qout(:,i) = q(:,i)
        h = q(1,i)
        hw = q(3,i)
        hp = q(4,i)
        !Exact integration of Source terms for EDC and BBBD systems
        ! call one_exact_step_EDC(h,hw,hp,csq,gamma,dt)
        ! call one_step_BE_EDC(h,hw,hp,csq,gamma,dt)
        call one_step_SDIRK22_EDC(h,hw,hp,csq,gamma,dt)
        ! call one_step_RKGL4_EDC(h,hw,hp,csq,gamma,dt)
        !Update output
        qout(3,i) = hw
        qout(4,i) = hp
    end do

end subroutine src_EDC_1D

subroutine src_EDC_1D_radial(q,r,qout,dt,csq,gamma,mx)
    implicit none
    integer, intent(in) :: mx
    double precision, intent(in) :: dt
    double precision, intent(in) :: q(4,mx)
    double precision, intent(in) :: r(mx)
    double precision, intent(in) :: gamma, csq
    double precision, intent(out) :: qout(4,mx)

    !Locals
    integer :: i
    double precision :: h,hw,hp
    double precision ::ratio,u

    do i=1,mx
        !Copy input to output
        qout(:,i) = q(:,i)
        h = q(1,i)
        hw = q(3,i)
        hp = q(4,i)
        !Exact integration of Source terms for EDC and BBBD systems
        ! call one_exact_step_EDC(h,hw,hp,csq,gamma,dt)
        call one_step_SDIRK22_EDC(h,hw,hp,csq,gamma,dt)
        !Update output
        qout(3,i) = hw
        qout(4,i) = hp
        !Radial terms
        ratio = dt/r(i)
        qout(1,i) = h-ratio*q(2,i)
        u = q(2,i)/h
        qout(2,i) = q(2,i)-ratio*(h*u**2)
        qout(3,i) = hw - ratio*u*hw
        qout(4,i) = hp - ratio*(u*hp+h*u*csq)
    end do
end subroutine src_EDC_1D_radial

subroutine src_BBBD_1D(q,qout,dt,csq,mx)
    implicit none
    integer, intent(in) :: mx
    double precision, intent(in) :: dt
    double precision, intent(in) :: q(6,mx)
    double precision, intent(in) :: csq
    double precision, intent(out) :: qout(6,mx)

    !Locals
    integer :: i
    double precision :: h,hw,hsg,hp,hpb

    do i=1,mx
        !Copy input to output
        qout(:,i) = q(:,i)
        h = q(1,i)
        hw = q(3,i)
        hsg = q(4,i)
        hp = q(5,i)
        hpb = q(6,i)
        !Exact integration of Source terms for EDC and BBBD systems
        ! call one_exact_step_BBBD(h,hw,hsg,hp,hpb,csq,dt)
        call one_step_SDIRK22_BBBD(h,hw,hsg,hp,hpb,csq,dt)
        !Update output
        qout(3,i) = hw
        qout(4,i) = hsg
        qout(5,i) = hp
        qout(6,i) = hpb
    end do

end subroutine src_BBBD_1D

subroutine src_BBBD_1D_radial(q,r,qout,dt,csq,mx)
    implicit none
    integer, intent(in) :: mx
    double precision, intent(in) :: dt
    double precision, intent(in) :: q(6,mx)
    double precision, intent(in) :: r(mx)
    double precision, intent(in) :: csq
    double precision, intent(out) :: qout(6,mx)

    !Locals
    integer :: i
    double precision :: h,hw,hsg,hp,hpb
    double precision :: ratio,u

    do i=1,mx
        !Copy input to output
        qout(:,i) = q(:,i)
        h = q(1,i)
        hw = q(3,i)
        hsg = q(4,i)
        hp = q(5,i)
        hpb = q(6,i)
        !Exact integration of Source terms for EDC and BBBD systems
        call one_exact_step_BBBD(h,hw,hsg,hp,hpb,csq,dt)
        !Update output
        qout(3,i) = hw
        qout(4,i) = hsg
        qout(5,i) = hp
        qout(6,i) = hpb
        !Radial terms
        ratio = dt/r(i)
        qout(1,i) = h-ratio*q(2,i)
        u = q(2,i)/h
        qout(2,i) = q(2,i)-ratio*(h*u**2)
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
!Exact integration of Source terms for EDC and BBBD systems
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

subroutine one_step_BE_EDC(h, hw, hp, csq_space,gammaEDC, dt)
            !A backward Euler step for the EDC equations
            implicit none
            double precision, intent(in) :: h, csq_space,gammaEDC, dt
            double precision, intent(inout) :: hw, hp
            !Locals
            double precision :: p,w
            double precision ::denom,b,c
            w = hw/h
            p = hp/h
            denom = 1.d0 + (dt*dt*2.d0*csq_space*gammaEDC)/(h*h)
            b = dt*gammaEDC/h
            c = -2.d0*dt*csq_space/h
            !Output (Inverting the matrix explicitly and multiplying back h)
            hw = h*(w + b*p)/denom 
            hp = h*(c*w + p)/denom
end subroutine one_step_BE_EDC


subroutine one_step_SDIRK22_EDC(h, hw, hp, csq_space, gammaEDC, dt)
    !Take one step using the SDIRK2(2) method
    !Intended to solve q_t = A(q) with A(q) = A*q
    !And q = (w,p), 
    implicit none
    double precision, intent(in) :: h, csq_space, gammaEDC, dt
    double precision, intent(inout) :: hw, hp
    !Matrix handling
    double precision, dimension(2,2) :: A
    double precision, dimension(2) :: qn, qnp1

    if (h<1.e-3 .or. csq_space < 1.e-3) then
        !Do nothing
        return
    end if
    qn = (/hw/h,hp/h/)
    A = 0.d0
    A(1,2) = gammaEDC/h
    A(2,1) = -2*csq_space/h

    call one_general_sdirk2_step(dt,A,qn,qnp1)
    hw = qnp1(1)*h
    hp = qnp1(2)*h
end subroutine one_step_SDIRK22_EDC

subroutine one_exact_step_BBBD(h, hw, hsg, hp, hpb, csq_space, dt)
    implicit none

    double precision, intent(in) :: h, csq_space, dt
    double precision, intent(inout) :: hw, hsg, hp, hpb
    !Locals
    complex(kind=8) :: q(4,1)
    complex(kind=8) :: R(4,4), Diag(4), Rinv(4,4)
    complex(kind=8) :: z,sq7,a1,a2,a3,a4,a5,a6,a7,a8,c
    complex(kind=8) :: img = cmplx(0.d0,1.d0)
    integer :: i,j
    integer :: ipiv(4), info


    c = sqrt(csq_space)
    sq7 = sqrt(7.d0)
    a1 = 3.d0 + sq7
    a2 = -3.d0 + sq7
    a3 = 2.d0 + sq7
    a4 = -2.d0 + sq7
    a5 = 1.d0 - sq7
    a6 = 1.d0 + sq7
    a7 = sqrt(2.d0 - 2.d0*sq7/3.d0)
    a8 = sqrt(1.d0 - sq7/3.d0)

    R(1,:) = (1.d0/sqrt(6.d0))*(/img/sqrt(a1),-img/sqrt(a1),img/sqrt(-a2),1.d0/sqrt(a2)/)
    R(2,:) = (/ -img*a3*a8, img*a3*a8, -img*a4*a7/a2, img*a4*a7/a2 /)
    R(3,:) = (c/6.d0)*(/ a5, a5, a6, a6 /)
    R(4,:) = (/ c,c,c,c/)

    Diag = (dt*sqrt(6.d0)*c/h)*(/ -img*sqrt(a1),img*sqrt(a1),-sqrt(a2),sqrt(a2) /)

    !Initial conditions
    q(:,1) = (/hw,hsg,hp,hpb/)/h
    !More efficient matmul
    !!!!!!!!!!!!
    !Old approach (numerically unstable)
    ! Rinv = complmatinv4(R) 
    !q = matmul(Rinv,q)
    !!!!!!!!!!!!
    ! Now use LAPACK to solve the linear system
    call zgesv(4, 1, R, 4, ipiv, q, 4, info)
    if (info /= 0) then
        print *, "ZGESV failed in one_exact_step_BBBD, info =", info
        stop
    end if
    do i=1,4
        q(i,1)=exp(Diag(i))*q(i,1)
    end do
    ! Multiply by R again (using matmul)
    ! since we overwrote R with LU factors, we rebuild it:
    R(1,:) = (1.d0/sqrt(6.d0))*(/img/sqrt(a1),-img/sqrt(a1),img/sqrt(-a2),1.d0/sqrt(a2)/)
    R(2,:) = (/ -img*a3*a8, img*a3*a8, -img*a4*a7/a2, img*a4*a7/a2 /)
    R(3,:) = (c/6.d0)*(/ a5, a5, a6, a6 /)
    R(4,:) = (/ c,c,c,c/)

    q = matmul(R,q)
    hw  = RealPart(q(1,1))*h
    hsg = RealPart(q(2,1))*h
    hp  = RealPart(q(3,1))*h
    hpb = RealPart(q(4,1))*h
end subroutine one_exact_step_BBBD

subroutine one_step_SDIRK22_BBBD(h, hw, hsg, hp, hpb, csq_space, dt)
    !Take one step using the SDIRK2(2) method
    !Intended to solve q_t = A(q) with A(q) = A*q
    !And q = (w,sg,p,pb), 
    implicit none
    double precision, intent(in) :: h, csq_space, dt
    double precision, intent(inout) :: hw, hsg, hp, hpb
    !Matrix handling
    double precision, dimension(4,4) :: A
    double precision, dimension(4) :: qn, qnp1

    if (h<1.e-3 .or. csq_space < 1.e-3) then
        !Do nothing
        return
    end if
    qn = (/hw/h,hsg/h,hp/h,hpb/h/)
    A = 0.d0
    A(1,4) = 1.d0/h
    A(2,3) = 12.d0/h
    A(2,4) = -6.d0/h
    A(3,2) = -csq_space/h
    A(4,1) = -6.d0*csq_space/h
    A(4,2) = 3.d0*csq_space/h

    call one_general_sdirk2_step(dt,A,qn,qnp1)
    hw = qnp1(1)*h
    hsg = qnp1(2)*h
    hp = qnp1(3)*h
    hpb = qnp1(4)*h
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