module bouss_module
    !Parameters for Boussinesq equations
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    use cparam_mod, only: boussEquations
    use cparam_mod, only: c_sq => boussEDCcsq 
    !Reference hyperbolic relaxation parameter (Just used for EDC)
    !Approximate SGN (3/2) or Sainte-Marie equations (2)
    use cparam_mod, only: gammaEDC => boussEDCgamma 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Mapped grid and geophysical parameters
    use cparam_mod, only:  g => grav, dry_tolerance, sea_level
    use cparam_mod, only: inclusion_idx, inclusion_type, srcintegrator
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    implicit none
    contains
    !!!! Timestepping routines for source terms
        subroutine one_step_BE_EDC(h, hw, hp, csq_space, dt)
            !A backward Euler step for the EDC equations
            implicit none
            double precision, intent(in) :: h, csq_space, dt
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

        subroutine one_step_SDIRK22_EDC(h, hw, hp, csq_space, dt)
            !Take one step using the SDIRK2(2) method
            !Intended to solve q_t = A(q) with A(q) = A*q
            !And q = (w,p), 
            implicit none
            double precision, intent(in) :: h, csq_space, dt
            double precision, intent(inout) :: hw, hp
            !Matrix handling
            double precision, dimension(2,2) :: A
            double precision, dimension(2) :: qn, qnp1

            if (h<dry_tolerance) then
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

        subroutine one_exact_step_EDC(h, hw, hp, csq_space, dt)
            implicit none
            double precision, intent(in) :: h, csq_space, dt
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
            expdtA(2,1) = -dsqrt(posc/b)*sin(dsqrt(b*posc)) !sin(i*a) = i*sinh(a)
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
    !!!!Riemann solvers    

end module bouss_module