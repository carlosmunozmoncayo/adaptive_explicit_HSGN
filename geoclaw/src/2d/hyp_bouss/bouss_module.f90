module bouss_module
   !  use amr_module, only: NEEDS_TO_BE_SET
   !  use geoclaw_module, only: g => grav
    use geoclaw_module, only: dry_tolerance

    implicit none

    !Additional parameters, not used now:
    integer ::  boussSolver
    ! parameters set in setrun:
    integer :: boussEquations       ! Which equations? 0=SWE, 1=Madsen, 2=SGN, -1=EDC, -2=BBBD
    integer:: boussMinLevel, boussMaxLevel
    double precision :: boussTstart, boussTfinal
    double precision :: boussEDCcsq        !Reference EDC (and BBBD) Hyperbolic relaxation parameter
    double precision :: boussEDCgamma       !Use EDC to approximate SGN (3/2) or Sainte-Marie equations (2)
    integer :: bouss_csq_index     ! Index of hyperbolic relaxation parameter in aux array
    integer :: boussDecayRate_index ! Index of decay rate in aux array
    integer :: boussTransitionType ! Transition method from EDCEs (and BBBDEs) to SWEs
    integer :: bouss_transition_type_fun ! Transition method (see transition_function below)
    double precision :: bouss_transition_low !Below this depth/distance from region use SWEs
    double precision :: bouss_transition_up !Above this depth/distance from region use BTEs
    real(kind=8), dimension(2) :: projection_center !(useful if coordinate_system == 2)

    !Global variables to turn on/off BTEs
    logical :: use_bouss = .FALSE.
    double precision :: small_tol_bouss = 1.d-13
    double precision :: large_tol_bouss = 1.d-6

    
    save


contains

!======================================================================

subroutine set_bouss
    
    ! Set things up for Boussinesq solver, in particular
    ! create and factor tridiagonal matrix for implicit solves.

    use geoclaw_module, only: earth_radius,deg2rad,coordinate_system,sea_level

    implicit none
    ! integer, intent(in) :: mx,my mbc, mthbc(2)
    integer :: i, iunit
    character*25 fname

    iunit = 7
    fname = 'bouss.data'
 !  # open the unit with new routine from Clawpack 4.4 to skip over
 !  # comment lines starting with #:
    call opendatafile(iunit, fname)

    read(iunit,*) boussEquations    ! which equations
    read(iunit,*) boussMinLevel     ! minimum level to use BTEs
    read(iunit,*) boussMaxLevel     ! maximum level to use BTEs
    read(iunit,*) boussSolver !Not used here
    read(iunit,*) boussTstart !Time to start using BTEs
    read(iunit,*) boussTfinal !Time to stop using BTEs
    read(iunit,*) boussEDCcsq       ! Reference EDC Hyperbolic relaxation parameter
    read(iunit,*) boussEDCgamma     ! Use EDC to approximate SGN (3/2) or Sainte-Marie equations (2)
    read(iunit,*) bouss_csq_index   ! Index of hyperbolic relaxation parameter in aux array
    read(iunit,*) boussDecayRate_index ! Index of decay rate in aux array
    read(iunit,*) boussTransitionType   ! Transition method from EDCEs to SWEs
    read(iunit,*) bouss_transition_type_fun ! Transition method (see transition_function below)
    read(iunit,*) bouss_transition_low !Below this depth/distance from region use SWEs
    read(iunit,*) bouss_transition_up !Above this depth/distance from region use BTEs
    read(iunit,*) projection_center ! useful if coordinate_system == 2

    print *, "########################"
    print *, "Reading Bouss data from file bouss.data, inside set_bouss"
    print *, "boussEquations = ", boussEquations
    print *, "boussMinLevel = ", boussMinLevel
    print *, "boussMaxLevel = ", boussMaxLevel
    print *, "boussSolver = ", boussSolver
    print *, "boussTmin = ", boussTstart
    print *, "boussTmax = ", boussTfinal
    print *, "boussEDCcsq = ", boussEDCcsq
    print *, "boussEDCgamma = ", boussEDCgamma
    print *, "bouss_csq_index = ", bouss_csq_index
    print *, "boussDecayRate_index = ", boussDecayRate_index
    print *, "boussTransitionType = ", boussTransitionType
    print *, "bouss_transition_type_fun = ", bouss_transition_type_fun
    print *, "bouss_transition_low = ", bouss_transition_low
    print *, "bouss_transition_up = ", bouss_transition_up
    print *, "projection_center = ", projection_center
    print *, "########################"
end subroutine set_bouss

!======================================================================
! Routines for Hyperbolic relaxations of dispersive water wave models
!======================================================================
!Escalante, Dumbser and Castro's system, 2019 (EDC)
subroutine riemann_normal_HypRel_EDC(mu, nv, grav, gamma, c_sq, &
                               hl, ul, vl, wl, pl, bl, csq_l,&
                               hr, ur, vr, wr, pr, br, csq_r, &
                               fw, sw)

   !Normal 2D Riemann solver for Escalante, Dumbser, and Castro's
   !hyperbolic relaxation of Sainte-Marie's dispersive water wave system
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !When used in 2D:
   !  Input: mu and nv according to the normal direction (analogous to Geoclaw's solver for SWEs)
   !  Output: 5 f-waves with 5 components each (fw) and 5 wave speeds (sw)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !When used in 1D :
   !  Input: mu=2, nv =3, vl=vr=0
   !  Output: The 3rd component of all the f-waves must be discarded
   !          The 3rd f-wave and the 3rd wave speed must be discarded
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !${CLAW}/riemann/src/geoclaw_riemann_utils.f is probably a better home for this routine
   !I'll leave it here for now
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

   !Input
   integer :: mu, nv
   double precision :: grav, gamma, c_sq
   double precision :: hl, ul, vl, wl, pl, bl, csq_l
   double precision :: hr, ur, vr, wr, pr, br, csq_r

   !Output
   double precision :: fw(5,5), sw(5)

   !Local
   double precision :: hlsq, hrsq
   double precision :: hbar, ubar, vbar, wbar, pbar, csqbar, CEbar
   double precision :: lamb1bar, lamb5bar
   double precision :: lamb1tilde, lamb5tilde
   double precision :: R(5,5), R_inv(5,5)
   double precision :: alpha(5,1), total_fluct(5,1)
   double precision :: etaL, etaR
   double precision :: use_full_bath, bath_correction_term(2)

   integer :: mw

   etaL = hl+bl!+sea_level not necessary since these two are subtracted below
   etaR = hr+br!+sea_level not necessary since these two are subtracted below
   !Roe averages
   hbar = 0.5*(hl+hr)
   hlsq = dsqrt(hl)
   hrsq = dsqrt(hr)
   ubar = (hlsq*ul+hrsq*ur)/(hlsq+hrsq)
   vbar = (hlsq*vl+hrsq*vr)/(hlsq+hrsq)
   wbar = (hlsq*wl+hrsq*wr)/(hlsq+hrsq)
   pbar = (hlsq*pl+hrsq*pr)/(hlsq+hrsq)

   !Still need to verify if there is difference between taking
   !the average, min, or max
   csqbar = 0.5d0*(csq_l+csq_r)

   CEbar=dsqrt(grav*hbar+csqbar+pbar)

   !Roe speeds
   lamb1tilde=ubar-CEbar
   lamb5tilde=ubar+CEbar

   !Einfeldt speeds
   lamb1bar=min(lamb1tilde, ul-dsqrt(grav*hl+csq_l+pl))
   lamb5bar=max(lamb5tilde, ur+dsqrt(grav*hr+csq_r+pr))
   
   !Matrix of right eigenvectors
   R(1,:)= [1.d0,       1.d0,      0.d0, 0.d0, 1.d0]
   R(mu,:)=[lamb1bar,   ubar,      0.d0, 0.d0, lamb5bar]
   R(nv,:)=[vbar,       0.d0,      0.d0, 1.d0, vbar]
   R(4,:)= [wbar,       0.d0,      1.d0, 0.d0, wbar]
   R(5,:)= [pbar+csqbar,-grav*hbar,0.d0, 0.d0, pbar+csqbar]

   use_full_bath = bath_sharpness_sensor(bl,br)

   !Computing total fluctuation to split in f-waves
   total_fluct(1,1)=hr*ur-hl*ul
   total_fluct(mu,1)=hr*ur**2+hr*pr+0.5d0*grav*hr**2-&
                  hl*ul**2-hl*pl-0.5d0*grav*hl**2+&
                  (grav*hbar)*(br-bl)
    bath_correction_term(1) = use_full_bath*gamma*pbar*(br-bl)
   total_fluct(nv,1)=ur*hr*vr-ul*hl*vl
   total_fluct(4,1)=ur*hr*wr-ul*hl*wl!-gamma*pbar !Handled with operator splitting
   total_fluct(5,1)=ur*hr*(pr+csq_r)-ul*hl*(pl+csq_l)!
   bath_correction_term(2) = -csqbar*ubar*((etaR-etaL)+use_full_bath*(br-bl)) 
   !The last term in the last expression is a consistency/stability term,
   !it is suppressed when there are large relative bathymetry jumps
   total_fluct(mu,1) = total_fluct(mu,1) + bath_correction_term(1)
   total_fluct(5,1) = total_fluct(5,1) + bath_correction_term(2)

   !Solving system explicitly for the 5x5 case
   R_inv(:,1)=[ubar*CEbar+grav*hbar,2*(csqbar+pbar),&
      -2*grav*hbar*wbar, -2*grav*hbar*vbar, -ubar*CEbar+grav*hbar]
   R_inv(:,mu)=[-CEbar, 0.d0, 0.d0, 0.d0, CEbar]
   R_inv(:,nv)=[0.d0, 0.d0, 0.d0, 2.d0*CEbar**2, 0.d0]
   R_inv(:,4)=[0.d0, 0.d0, 2*CEbar**2, 0.d0, 0.d0]
   R_inv(:,5)=[1.d0,-2.d0,-2*wbar,-2*vbar,1.d0]
   R_inv=(0.5d0/CEbar**2)*R_inv
   
   alpha = matmul(R_inv,total_fluct)
   
   ! Computing f-waves
   do mw = 1, 5 !Replace by mwaves if necessary
      fw(:, mw) = alpha(mw,1) * R(:, mw)
   END do
   
   !Computing wave speeds
   sw(1) = min(ul-sqrt(grav*hl+csqbar+pl),lamb1bar)
   sw(2) = ubar
   sw(3) = ubar
   sw(4) = ubar
   sw(5) = max(ur+sqrt(grav*hr+csqbar+pr),lamb5bar)

end subroutine riemann_normal_HypRel_EDC


subroutine riemann_normal_HypRel_DR(mu,nv, g, alphaLDR, c_sq, &
   hl, ul, vl, wl, pl, sgl, bl, csq_l, &
   hr, ur, vr, wr, pr, sgr, br, csq_r, &
   fw, sw)
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !Normal Riemann solver for Duran and Richard's system
   !Since the matrix of right eigenvectors is invertible, but its
   !analytical inverse is very ugly, we invert it numerically.
   !before I was using LU decomposition from lu_module, but now I
   !am switching lo LAPACK's dgesv routine for better everything
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   implicit none

   !Input
   integer :: mu, nv
   double precision :: g, alphaLDR, c_sq
   double precision :: hl, ul, vl, wl, pl, sgl, bl, csq_l
   double precision :: hr, ur, vr, wr, pr, sgr, br, csq_r

   !Output
   double precision :: fw(6,6), sw(6)

   !Local
   double precision :: hlsq, hrsq
   double precision :: h, u, v, w, p, sigma
   double precision :: saux(7) !Auxiliary variables
   double precision :: lamb_bar(6), lamb_tilde(6)
   double precision :: R(6,6), Rsolve(6,6) !R would get overriden during the linear solve
   double precision :: Rinv(6,6), alpha(6,1), denom_Rinv
   double precision :: total_fluct(6,1), total_fluct_solve(6,1)
   double precision :: bath_correction_term(3)
   ! double precision :: bath_sharpness_sensor !For stability term !No need to declare again
   integer :: j
   integer :: mw, j_debug

   !Variables for LAPACK dgesv
   integer :: n, info, lda, ldb, nrhs
   integer, allocatable :: ipiv(:)
   !external declarations 
   external dgetrf, dgetrs

   integer :: M, RANK, LWORK
   double precision :: RCOND
   double precision, allocatable :: A(:,:), B(:,:), singular(:),  WORK(:)
   double precision :: WORKQ(1)
   integer :: linear_solve_flag

   linear_solve_flag = 0 !Set to 1 if some part of the linear solve fails

   ! print*, "alphaLDR = ", alphaLDR
   !Roe averages
   h = 0.5d0*(hl+hr)
   hlsq = dsqrt(hl)
   hrsq = dsqrt(hr)
   u = (hlsq*ul+hrsq*ur)/(hlsq+hrsq)
   v = (hlsq*vl+hrsq*vr)/(hlsq+hrsq)
   w = (hlsq*wl+hrsq*wr)/(hlsq+hrsq)
   p = (hlsq*pl+hrsq*pr)/(hlsq+hrsq)
   sigma = (hlsq*sgl+hrsq*sgr)/(hlsq+hrsq)

   !Still need to verify if there is difference between taking
   !the average, min, or max
   c_sq = 0.5d0*(csq_l+csq_r)

   saux(:) = [dsqrt(g*h+p+c_sq*alphaLDR), &
               g*h+alphaLDR*(c_sq*alphaLDR+p), &
               dsqrt(h)*w, &
               alphaLDR*(c_sq*alphaLDR+p),&
               0.d0,&
               0.d0,&
               sqrt(g*(alphaLDR-1.d0)/alphaLDR)]
   saux(5) = 3.d0*saux(3)*(g*h*alphaLDR+saux(4))
   saux(6) = saux(1)*sigma*(g*h+saux(4))

   !Roe speeds
   lamb_tilde(1) = u-saux(1)
   lamb_tilde(2) = u-sqrt(h)*saux(7)
   lamb_tilde(3) = u
   lamb_tilde(4) = u
   lamb_tilde(5) = u+sqrt(h)*saux(7)
   lamb_tilde(6) = u+saux(1)

   !Einfeldt speeds
   lamb_bar(1) = min(lamb_tilde(1), ul-dsqrt(g*hl+csq_l*alphaLDR+pl))
   lamb_bar(2) = min(lamb_tilde(2), ul-dsqrt(hl)*saux(7))
   lamb_bar(3) = u
   lamb_bar(4) = u
   lamb_bar(5) = max(lamb_tilde(4), ur+dsqrt(hr)*saux(7))
   lamb_bar(6) = max(lamb_tilde(5), ur+dsqrt(g*hr+csq_r*alphaLDR+pr))

   !Roe matrix of right eigenvectors
   R(1,:) = [saux(2)*(saux(5) + saux(6)),0.d0,0.d0,1.d0,0.d0,saux(2)*(saux(5) - saux(6))]
   R(mu,:) = [(u - saux(1))*saux(2)*(saux(5) + saux(6)),0.d0,0.d0,u,0.d0,(u + saux(1))*saux(2)*(saux(5) - saux(6))]
   R(nv,:) = [v*saux(2)*(saux(5) + saux(6)),0.d0,1.d0,0.d0,0.d0,v*saux(2)*(saux(5) - saux(6))]
   R(4,:) = [(w*((-1 + 3*alphaLDR)*g*h + 2*saux(4))*(saux(5) + saux(6)))/2.,saux(7)/2.,0.d0,-0.5*w,&
               -0.5*saux(7),(w*((-1 + 3*alphaLDR)*g*h + 2*saux(4))*(saux(5) - saux(6)))/2.]
   R(5,:) = [(alphaLDR*c_sq + p)*saux(2)*(saux(5) + saux(6)),0.d0,0.d0,-(g*h),0.d0,(alphaLDR*c_sq + p)*saux(2)*(saux(5) - saux(6))]
   R(6,:) = [(sigma*saux(2) + 3*alphaLDR*saux(1)*saux(3))*(saux(5) + saux(6)),1.d0,0.d0,sigma,1.d0,&
               (-3*alphaLDR*saux(1)*saux(3) + sigma*(g*h + saux(4)))*(saux(5) - saux(6))]
   Rsolve = R !Copying R to Rsolve since R would get overriden during the linear solve

   
   !Computing total fluctuation to split in f-waves
   total_fluct(1,1)=hr*ur-hl*ul
   total_fluct(mu,1)=hr*ur**2+hr*pr+0.5d0*g*hr**2-&
                  hl*ul**2-hl*pl-0.5d0*g*hl**2+&
                  g*h*(br-bl)
   total_fluct(nv,1)=hr*ur*vr-hl*ul*vl
   total_fluct(4,1)=ur*hr*wr-ul*hl*wl+&
                     0.5d0*((1.d0-alphaLDR)/alphaLDR)*((hr**1.5)*sgr-(hl**1.5)*sgl)
   bath_correction_term(1) = 0.75d0*((alphaLDR-1.d0)/alphaLDR)*g*sigma*dsqrt(h)*(hr-hl) !Non-conservative term
   total_fluct(5,1)= hr*ur*(pr+c_sq*alphaLDR)-hl*ul*(pl+c_sq*alphaLDR) 
   bath_correction_term(2) = -c_sq*alphaLDR*u*(hr-hl) !Non-conservative term
   total_fluct(6,1)=hr*ur*sgr-2.d0*dsqrt(hr)**3.d0*wr-&
                     (hl*ul*sgl-2.d0*dsqrt(hl)**3.d0*wl)

   !Call sharpness sensor to avoid instability due to bathymetry
   total_fluct(4,1) = total_fluct(4,1) + bath_correction_term(1)
   total_fluct(5,1) = total_fluct(5,1) + bath_sharpness_sensor(bl,br)*bath_correction_term(2)

   !To invert R, first check its conditioning, we just need to check
   !that the denominator of the analytical inverse is not too small
   denom_Rinv = (saux(5)-saux(6))*(saux(5)+ saux(6))*saux(1)*saux(2)*saux(7)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   if(maxval(abs(total_fluct)) < 1.d-30) then
   !If nothing is happpening at this interface, move on
         total_fluct(:,1) = 0.d0
   else if (abs(denom_Rinv) > 1.e-16) then
      !R well-conditioned, build Rinv analytically and solve directly
      Rinv(:,1) = [(g*h + u*saux(1))/(2.*saux(1)**2*saux(2)*(saux(5) + saux(6))),&
         (g*h*w*((1 - 3*alphaLDR)*g*h - 2*saux(4)) - 3*alphaLDR*u*saux(1)**2*saux(3)*saux(7) +&
         saux(2)*((alphaLDR*c_sq + p)*w - (alphaLDR*c_sq + g*h + p)*sigma*saux(7)))/(2.*saux(1)**2*saux(2)*saux(7)),&
         -((g*h*v)/saux(1)**2),1 - (g*h)/saux(1)**2,(g*h*w*((-1 + 3*alphaLDR)*g*h + 2*saux(4)) -&
         3*alphaLDR*u*saux(1)**2*saux(3)*saux(7) + saux(2)*(-((alphaLDR*c_sq + p)*w) - (alphaLDR*c_sq + g*h +&
         p)*sigma*saux(7)))/(2.*saux(1)**2*saux(2)*saux(7)),(g*h - u*saux(1))/(2.*saux(1)**2*saux(2)*(saux(5) - saux(6)))]
      Rinv(:,mu) = [-0.5*1/(saux(1)*saux(2)*(saux(5) + saux(6))),(3*alphaLDR*saux(3))/(2.*saux(2)),0.d0,0.d0,&
         (3*alphaLDR*saux(3))/(2.*saux(2)),1/(2*saux(1)*saux(2)*saux(5) - 2*saux(1)*saux(2)*saux(6))]
      Rinv(:,nv) = [0.d0,0.d0,1.d0,0.d0,0.d0,0.d0]
      Rinv(:,4) = [0.d0,1/saux(7),0.d0,0.d0,-(1/saux(7)),0.d0]
      Rinv(:,5) = [1/(2.*saux(1)**2*saux(2)*(saux(5) + saux(6))),-0.5*(w*((-1 + 3*alphaLDR)*g*h + saux(2) +&
         2*saux(4)))/(saux(1)**2*saux(2)*saux(7)),-(v/saux(1)**2),-saux(1)**(-2),(w*((-1 + 3*alphaLDR)*g*h + saux(2) +&
         2*saux(4)))/(2.*saux(1)**2*saux(2)*saux(7)),1/(2.*saux(1)**2*saux(2)*(saux(5) - saux(6)))]
      Rinv(:,6) = [0d0,0.5d0,0d0,0d0,0.5d0,0d0]

      total_fluct_solve = matmul(Rinv,total_fluct)
      total_fluct = total_fluct_solve
   else
      ! print *, "Warning: Matrix R is ill-conditioned in riemann_1D_HypRel_DR, attempting robust linear solve."
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Some columns of R are nearly linearly dependent, so Rinv is ill-conditioned
      !Try using LU (fast but less robust)
      !If that fails, use least-squares with SVD (slower but more robust)
      !If that fails, set total_fluct to zero and warn user
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      total_fluct_solve = total_fluct !Backup copy in case we need it for SVD
      !Proceed to linear solve
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Using LU decomposition from LAPACK dgesv
      n = 6
      nrhs = 1
      lda = 6
      ldb = 6
      ! allocate(ipiv(n))
      ! call dgetrf(n, n, Rsolve, lda, ipiv, info)
      ! if (info /= 0) then
      !    !If factorization failed, set total_fluct to zero, but warn user
      !    ! print *, "Warning: LAPACK dgesv failed to factorize matrix in riemann_1D_HypRel_DR."
      !    ! total_fluct(:) = 0.d0
         linear_solve_flag = 1
      ! else 
      !    !Solve the system
      !    call dgetrs('N', n, nrhs, Rsolve, lda, ipiv, total_fluct, ldb, info)
      !    if (info /= 0) then
      !       !If solve failed, set total_fluct to zero, but warn user
      !       ! total_fluct(:) = 0.d0
      !       ! print *, "Warning: LAPACK dgesv failed to solve for wave strengths in riemann_1D_HypRel_DR."
      !       linear_solve_flag = 1
      !    end if
      ! end if
      ! deallocate(ipiv)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !If LU failed, use SVD-based least squares
      if (linear_solve_flag == 1) then
         total_fluct = total_fluct_solve !Restore backup copy
         M=6
         RCOND = 1.e-12
         allocate(A(5,5), B(5,1), singular(5))
         A = R
         B(:,1) = total_fluct(:,1)     
         ! --- Workspace query ---
         LWORK = -1

         call dgelss(M, N, NRHS, A, LDA, B, LDB, singular, RCOND, RANK, WORKQ, LWORK, INFO)

         LWORK = int(WORKQ(1))
         allocate(WORK(LWORK))

         ! Restore RHS because DGELSS overwrote B during the query
         B(:,1) = total_fluct(:,1) 

         ! --- Actual solve ---
         call dgelss(M, N, NRHS, A, LDA, B, LDB, singular, RCOND, RANK, WORK, LWORK, INFO)

         if (INFO /= 0) then
            print *, "DGELSS failed, info=", INFO
            total_fluct(:,1)  = 0d0
         else
            total_fluct(:,1)  = B(:,1)
         end if
         ! print *, "Singular values of R:", singular

         !Check if total_fluct is greater than a small threshold
         ! if (maxval(abs(total_fluct)) > 1.d-14) then
         !    print*, 'total fluct not zeroed'
         ! end if
         deallocate(A, B, singular, WORK)
      end if
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   end if
   ! Computing f-waves
   do mw = 1, 6 !Replace by mwaves if necessary
      fw(:, mw) = total_fluct(mw,1) * R(:, mw)
   END do
   
   !Computing wave speeds
   do j = 1, 6
      sw(j) = lamb_bar(j)
   end do

end subroutine riemann_normal_HypRel_DR

subroutine riemann_transverse_HypRel_DR(ixy,u_input,v_input,h,w,p,sigma,g,alphaLDR,c_sq,localasdq, &
                                        waves,s)
   implicit none
   !Input
   integer :: ixy, k
   double precision :: u_input, v_input, h, w, p,sigma, c_sq, localasdq(6,1)
   double precision :: g, alphaLDR

   !Output
   double precision ::  waves(6,6), s(6)

   !Local
   double precision :: CE, lamb(6), R(6,6),  Rsolve(6,6) 
   double precision :: lolcalasdq_solve(6,1), saux(7)
   double precision :: alpha(6,1), Rinv(6,6)
   double precision :: denom_Rinv
   double precision :: u,v !Local velocities in normal and transverse directions
   integer :: mu, nv

   !Variables for LAPACK dgesv
   integer :: n, info, lda, ldb, nrhs
   integer, allocatable :: ipiv(:)
   !external declarations 
   external dgetrf, dgetrs

   integer :: M, RANK, LWORK
   double precision :: RCOND
   double precision, allocatable :: A(:,:), B(:,:), singular(:),  WORK(:)
   double precision :: WORKQ(1)
   integer :: linear_solve_flag

   linear_solve_flag = 0 !Set to 1 if some part of the linear solve fails
   if (ixy.eq.1) then
      !Horizontal normal swipe (Construct R^y(u,v))
      u = v_input
      v = u_input
      mu = 3
      nv = 2
   else
      !Vertical normal swipe (Construct R^x(u,v))
      u = u_input
      v = v_input
      mu = 2
      nv = 3
   end if

   saux(:) = [dsqrt(g*h+p+c_sq*alphaLDR), &
               g*h+alphaLDR*(c_sq*alphaLDR+p), &
               dsqrt(h)*w, &
               alphaLDR*(c_sq*alphaLDR+p),&
               0.d0,&
               0.d0,&
               sqrt(g*(alphaLDR-1.d0)/alphaLDR)]
   saux(5) = 3.d0*saux(3)*(g*h*alphaLDR+saux(4))
   saux(6) = saux(1)*sigma*(g*h+saux(4))

        
   R(1,:) = [saux(2)*(saux(5) + saux(6)),0.d0,0.d0,1.d0,0.d0,saux(2)*(saux(5) - saux(6))]
   R(mu,:) = [(u - saux(1))*saux(2)*(saux(5) + saux(6)),0.d0,0.d0,u,0.d0,(u + saux(1))*saux(2)*(saux(5) - saux(6))]
   R(nv,:) = [v*saux(2)*(saux(5) + saux(6)),0.d0,1.d0,0.d0,0.d0,v*saux(2)*(saux(5) - saux(6))]
   R(4,:) = [(w*((-1 + 3*alphaLDR)*g*h + 2*saux(4))*(saux(5) + saux(6)))/2.,saux(7)/2.,0.d0,-0.5*w,&
               -0.5*saux(7),(w*((-1 + 3*alphaLDR)*g*h + 2*saux(4))*(saux(5) - saux(6)))/2.]
   R(5,:) = [(alphaLDR*c_sq + p)*saux(2)*(saux(5) + saux(6)),0.d0,0.d0,-(g*h),0.d0,(alphaLDR*c_sq + p)*saux(2)*(saux(5) - saux(6))]
   R(6,:) = [(sigma*saux(2) + 3*alphaLDR*saux(1)*saux(3))*(saux(5) + saux(6)),1.d0,0.d0,sigma,1.d0,&
               (-3*alphaLDR*saux(1)*saux(3) + sigma*(g*h + saux(4)))*(saux(5) - saux(6))]

   s(1) = u - saux(1)
   s(2) = u - sqrt(h)*saux(7)
   s(3) = u
   s(4) = u
   s(5) = u + sqrt(h)*saux(7)
   s(6) = u + saux(1)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !To invert R, first check its conditioning, we just need to check
   !that the denominator of the analytical inverse is not too small
   denom_Rinv = (saux(5)-saux(6))*(saux(5)+ saux(6))*saux(1)*saux(2)*saux(7)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   if(maxval(abs(localasdq)) < 1.d-30) then
   !If nothing is happpening at this interface, move on
         lolcalasdq_solve(:,1) = 0.d0
   else if (abs(denom_Rinv) > 1.e-16) then
      !Building Rinv analytically
      Rinv(:,1) = [(g*h + u*saux(1))/(2.*saux(1)**2*saux(2)*(saux(5) + saux(6))),&
         (g*h*w*((1 - 3*alphaLDR)*g*h - 2*saux(4)) - 3*alphaLDR*u*saux(1)**2*saux(3)*saux(7) +&
         saux(2)*((alphaLDR*c_sq + p)*w - (alphaLDR*c_sq + g*h + p)*sigma*saux(7)))/(2.*saux(1)**2*saux(2)*saux(7)),&
         -((g*h*v)/saux(1)**2),1 - (g*h)/saux(1)**2,(g*h*w*((-1 + 3*alphaLDR)*g*h + 2*saux(4)) -&
         3*alphaLDR*u*saux(1)**2*saux(3)*saux(7) + saux(2)*(-((alphaLDR*c_sq + p)*w) - (alphaLDR*c_sq + g*h +&
         p)*sigma*saux(7)))/(2.*saux(1)**2*saux(2)*saux(7)),(g*h - u*saux(1))/(2.*saux(1)**2*saux(2)*(saux(5) - saux(6)))]
      Rinv(:,mu) = [-0.5*1/(saux(1)*saux(2)*(saux(5) + saux(6))),(3*alphaLDR*saux(3))/(2.*saux(2)),0.d0,0.d0,&
         (3*alphaLDR*saux(3))/(2.*saux(2)),1/(2*saux(1)*saux(2)*saux(5) - 2*saux(1)*saux(2)*saux(6))]
      Rinv(:,nv) = [0.d0,0.d0,1.d0,0.d0,0.d0,0.d0]
      Rinv(:,4) = [0.d0,1/saux(7),0.d0,0.d0,-(1/saux(7)),0.d0]
      Rinv(:,5) = [1/(2.*saux(1)**2*saux(2)*(saux(5) + saux(6))),-0.5*(w*((-1 + 3*alphaLDR)*g*h + saux(2) +&
         2*saux(4)))/(saux(1)**2*saux(2)*saux(7)),-(v/saux(1)**2),-saux(1)**(-2),(w*((-1 + 3*alphaLDR)*g*h + saux(2) +&
         2*saux(4)))/(2.*saux(1)**2*saux(2)*saux(7)),1/(2.*saux(1)**2*saux(2)*(saux(5) - saux(6)))]
      Rinv(:,6) = [0d0,0.5d0,0d0,0d0,0.5d0,0d0]

      lolcalasdq_solve = matmul(Rinv,localasdq)
   else
      Rsolve = R !Copying R to Rsolve since R would get overriden during the linear solve
      lolcalasdq_solve(:,1) = localasdq(:,1) !Backup copy in case we need it for SVD
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Then try using LU (fast but less robust)
      !If that fails, use least-squares with SVD (slower but more robust)
      !If that fails, set localasdq to zero and warn user
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Using LU decomposition from LAPACK dgesv
      n = 6
      nrhs = 1
      lda = 6
      ldb = 6
      ! allocate(ipiv(n))
      ! call dgetrf(n, n, Rsolve, lda, ipiv, info)
      ! if (info /= 0) then
      !    linear_solve_flag = 1
      ! else
      !    !Solve the system
      !    call dgetrs('N', n, nrhs, Rsolve, lda, ipiv, lolcalasdq_solve, ldb, info)
      !    if (info /= 0) then
            linear_solve_flag = 1
      !    end if
      ! end if
      ! deallocate(ipiv)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !If LU failed, use SVD-based least squares
      if (linear_solve_flag == 1) then
         lolcalasdq_solve(:,1) = localasdq(:,1) !Restore backup copy
         M=6
         RCOND = 1.e-12
         allocate(A(6,6), B(6,1), singular(6))
         A = R
         B(:,1) = lolcalasdq_solve(:,1)    
         ! --- Workspace query ---
         LWORK = -1
         call dgelss(M, N, NRHS, A, LDA, B, LDB, singular, RCOND, RANK, WORKQ, LWORK, INFO)
         LWORK = int(WORKQ(1))
         allocate(WORK(LWORK))
         ! Restore RHS because DGELSS overwrote B during the query
         B(:,1) = localasdq(:,1)
         ! --- Actual solve ---
         call dgelss(M, N, NRHS, A, LDA, B, LDB, singular, RCOND, RANK, WORK, LWORK, INFO)
         if (INFO /= 0) then
            print *, "DGELSS failed, info=", INFO
            lolcalasdq_solve(:,1) = 0.d0
         else
            lolcalasdq_solve(:,1) = B(:,1)
         end if
         deallocate(A, B, singular, WORK)
      end if
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   end if
   do k = 1, 6  !Compute waves
      waves(:, k) = lolcalasdq_solve(k,1) * R(:, k)
   END do
   
end subroutine riemann_transverse_HypRel_DR


function bath_sharpness_sensor(bL,bR)
   !Computes a sensor to detect sharp bathymetry changes
   !relative to the water depth. Returns a value between 0 and 1
   !Cubic polynomial transition

   implicit none

   !Input
   double precision,intent(in) :: bL, bR
   !Output
   double precision :: bath_sharpness_sensor

   double precision :: delta_b, delta_h, ratio, ratio_tolerance, eps
   double precision :: t
   delta_b = abs(bR - bL)
   delta_h = max(0.5d0*abs(bL + bR), 1.e-6) !Safeguard against divisions by zero

   if (delta_h < 1.e-6) then
      bath_sharpness_sensor = 1.d0 !Full correction if very shallow water
      return
   end if

   ratio = delta_b / delta_h
   ratio_tolerance = 0.01d0
   eps = 0.05d0
   !Smooth cubic transition between 1 and 0
   if (ratio <= ratio_tolerance) then
      bath_sharpness_sensor = 1.d0
   else if (ratio > ratio_tolerance .and. ratio < ratio_tolerance + eps) then
      t = (ratio - ratio_tolerance)/eps
      bath_sharpness_sensor = -3.d0*t**2 + 2.d0*t**3 + 1.d0
   else
      bath_sharpness_sensor = 0.d0
   end if

end function bath_sharpness_sensor

subroutine riemann_transverse_HypRel_EDC(ixy,u,v,h,w,p,grav,c_sq,localasdq, &
                                        waves,s)
    implicit none
    !Input
    integer :: ixy, k
    double precision :: u, v, h, w, p, c_sq, localasdq(5,1)
    double precision :: grav

    !Output
    double precision ::  waves(5,5), s(5)

    !Local
    double precision :: CE, lamb1, lamb5, R(5,5),  R_inv(5,5)
    double precision :: alpha(5,1)

    if (ixy.eq.1) then
        !Horizontal normal swipe (Construct R^y(u,v))
        CE=dsqrt(grav*h+c_sq+p)
        lamb1 = v-CE
        lamb5 = v+CE
        R(1,:)=[1.d0,   1.d0,   0.d0,   0.d0,   1.d0]
        R(2,:)=[u,      0.d0,   0.d0,   1.d0,   u]
        R(3,:)=[lamb1,  v,      0.d0,   0.d0,   lamb5]
        R(4,:)=[w,      0.d0,   1.d0,   0.d0,   w]
        R(5,:)=[p+c_sq, -grav*h,0.d0,   0.d0,   p+c_sq]

        R_inv(:,1)=[v*CE+grav*h,2*(c_sq+p),-2*grav*h*w, -2*grav*h*u, -v*CE+grav*h]
        R_inv(:,2)=[0.d0, 0.d0, 0.d0, 2.d0*CE**2, 0.d0]
        R_inv(:,3)=[-CE, 0.d0, 0.d0, 0.d0, CE]
        R_inv(:,4)=[0.d0, 0.d0, 2*CE**2, 0.d0, 0.d0]
        R_inv(:,5)=[1.d0,-2.d0,-2*w,-2*u,1.d0]
        R_inv=(0.5d0/CE**2)*R_inv

        s(1) = lamb1
        s(2) = v
        s(3) = v
        s(4) = v
        s(5) = lamb5
    else
        !Vertical normal swipe (Construct R^x(u,v))
        CE=dsqrt(grav*h+c_sq+p)
        lamb1 = u-CE
        lamb5 = u+CE
        R(1,:)=[1.d0,   1.d0,   0.d0,   0.d0,   1.d0]
        R(2,:)=[lamb1,  u,      0.d0,   0.d0,   lamb5]
        R(3,:)=[v,      0.d0,   0.d0,   1.d0,   v]
        R(4,:)=[w,      0.d0,   1.d0,   0.d0,   w]
        R(5,:)=[p+c_sq, -grav*h,0.d0,   0.d0,   p+c_sq]

        R_inv(:,1)=[u*CE+grav*h,2*(c_sq+p),-2*grav*h*w, -2*grav*h*v, -u*CE+grav*h]
        R_inv(:,2)=[-CE, 0.d0, 0.d0, 0.d0, CE]
        R_inv(:,3)=[0.d0, 0.d0, 0.d0, 2.d0*CE**2, 0.d0]
        R_inv(:,4)=[0.d0, 0.d0, 2*CE**2, 0.d0, 0.d0]
        R_inv(:,5)=[1.d0,-2.d0,-2*w,-2*v,1.d0]
        R_inv=(0.5d0/CE**2)*R_inv

        s(1) = lamb1
        s(2) = u
        s(3) = u
        s(4) = u
        s(5) = lamb5
    end if

    !Compute alphas
    alpha = matmul(R_inv,localasdq)
    do k=1,5
        waves(:,k) = alpha(k,1)*R(:,k)
    end do
end subroutine riemann_transverse_HypRel_EDC


function Source_HypRel_EDC(hw, hp, h, csq_space)
    implicit none
    double precision :: Source_HypRel_EDC(2), hw, hp, h, csq_space
    
    if (abs(h) > dry_tolerance) then
        Source_HypRel_EDC(1) = boussEDCgamma*hp/h
        Source_HypRel_EDC(2) = -2.d0*csq_space*hw/h
    end if
end function 

!Bassi, Bonaventura, Busto, Dumbser, 2020 (BBBD)
subroutine riemann_normal_HypRel_BBBD(mu, nv, grav, c_sq, &
    hl, ul, vl, wl, sgl, pl, pbl, bl, csq_l,&
    hr, ur, vr, wr, sgr, pr, pbr, br, csq_r, &
    fw, sw)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !When used in 2D:
    !  Input: mu and nv according to the normal direction (analogous to Geoclaw's solver for SWEs)
    !  Output: 7 f-waves with 7 components each (fw) and 7 wave speeds (sw)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !When used in 1D :
    !  Input: mu=2, nv =3, vl=vr=0
    !  Output: The 3rd component of all the f-waves must be discarded
    !          The 3rd f-wave and the 3rd wave speed must be discarded
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !${CLAW}/riemann/src/geoclaw_riemann_utils.f is probably a better home for this routine
    !I'll leave it here for now
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    implicit none

    !Input
    integer, intent(in) :: mu, nv
    double precision, intent(in) :: grav, c_sq
    double precision, intent(in) :: hl, ul, vl, wl, pl, bl, csq_l
    double precision, intent(in) :: hr, ur, vr, wr, pr, br, csq_r
    double precision, intent(in) :: sgl, sgr, pbl, pbr

    !Output
    double precision, intent(out) :: fw(7,7), sw(7)

    !Local
    double precision :: hlsq, hrsq, etaL, etaR
    double precision :: hbar, ubar, vbar, wbar, pbar, csqbar, CEbar
    double precision :: sgbar, pbbar
    double precision :: delb
    double precision :: lamb1bar, lamb7bar
    double precision :: R(7,7), R_inv(7,7)
    double precision :: alpha(7,1), total_fluct(7,1)
    double precision :: use_full_bath, bath_correction_term(3)

    integer :: mw

    etaL = hL + bL
    etaR = hR + bR

    !Roe averages
    hbar = 0.5*(hl+hr)
    hlsq = dsqrt(hl)
    hrsq = dsqrt(hr)
    ubar = (hlsq*ul+hrsq*ur)/(hlsq+hrsq)
    vbar = (hlsq*vl+hrsq*vr)/(hlsq+hrsq)
    wbar = (hlsq*wl+hrsq*wr)/(hlsq+hrsq)
    sgbar = (hlsq*sgl+hrsq*sgr)/(hlsq+hrsq)
    pbar = (hlsq*pl+hrsq*pr)/(hlsq+hrsq)
    pbbar = (hlsq*pbl+hrsq*pbr)/(hlsq+hrsq)

    !Still need to verify if there is difference between taking
    !the average, min, or max
    csqbar = 0.5d0*(csq_l+csq_r)

    CEbar=dsqrt(grav*hbar+csqbar+pbar)

    !Roe speeds
    lamb1bar=ubar-CEbar
    lamb7bar=ubar+CEbar

    !Matrix of right eigenvectors
    R(1,:)= [1.d0,       1.d0,      0.d0, 0.d0, 0.d0, 0.d0, 1.d0]
    R(mu,:)=[lamb1bar,   ubar,      0.d0, 0.d0, 0.d0, 0.d0, lamb7bar]
    R(nv,:)=[vbar,       0.d0,      0.d0, 0.d0, 1.d0, 0.d0, vbar]
    R(4,:)= [wbar,       0.d0,      0.d0, 1.d0, 0.d0, 0.d0, wbar]
    R(5,:)= [sgbar,      0.d0,      1.d0, 0.d0, 0.d0, 0.d0,  sgbar]
    R(6,:)= [pbar+csqbar,-grav*hbar,0.d0, 0.d0, 0.d0, 0.d0, pbar+csqbar]
    R(7,:)= [pbbar,      0.d0,      0.d0, 0.d0, 0.d0, 1.d0,  pbbar]
    
    delb = br-bl
    !Computing total fluctuation to split in f-waves
    total_fluct(1,1)=hr*ur-hl*ul
    bath_correction_term(1) = pbbar*delb
    total_fluct(mu,1) = hr*ur**2.d0 +hr*pr+0.5d0*grav*hr**2.d0 &
                    - hl*ul**2.d0 -hl*pl-0.5d0*grav*hl**2.d0 &
                    + grav*hbar*delb
    total_fluct(nv,1)=hr*ur*vr-hl*ul*vl
    total_fluct(4,1)=hr*ur*wr-hl*ul*wl
    total_fluct(5,1)=hr*ur*sgr-hl*ul*sgl
    bath_correction_term(2) = csqbar*ubar*delb
    total_fluct(6,1)=hr*ur*(csq_r+pr)-hl*ul*(csq_l+pl)-csqbar*ubar*(etaR-etaL)
    bath_correction_term(3) = -6.d0*csqbar*ubar*delb
    total_fluct(7,1)=hr*ur*pbr-hl*ul*pbl

    use_full_bath = bath_sharpness_sensor(bl,br)
    total_fluct(mu,1) = total_fluct(mu,1) + use_full_bath*bath_correction_term(1)
    total_fluct(6,1) = total_fluct(6,1) + use_full_bath*bath_correction_term(2)
    total_fluct(7,1) = total_fluct(7,1) + use_full_bath*bath_correction_term(3)

    
    R_inv(:,1)=(1.d0/CEbar**2.d0)*[0.5d0*(grav*hbar+CEbar*ubar),&
    csqbar+pbar, -grav*hbar*sgbar, -grav*hbar*wbar, -grav*hbar*vbar,&
    -grav*hbar*pbbar,0.5d0*(grav*hbar-CEbar*ubar)]
    R_inv(:,mu)=[-0.5d0/CEbar, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.5d0/CEbar]
    R_inv(:,nv)=[0.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 0.d0]
    R_inv(:,4)= [0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0]
    R_inv(:,5)= [0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0,  0.d0]
    R_inv(:,6)=(1.d0/CEbar**2.d0)*[0.5d0,-1.d0,-sgbar,&
        -wbar,-vbar, -pbbar, 0.5d0]
    R_inv(:,7)=[0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0]

    alpha = matmul(R_inv,total_fluct)

    ! Computing f-waves
    do mw = 1, 7 !Replace by mwaves if necessary
        fw(:, mw) = alpha(mw,1) * R(:, mw)
    END do

    !Computing wave speeds
    !Einfeldt speeds for outer waves
    sw(1) = min(ul-sqrt(grav*hl+csqbar+pl),lamb1bar)
    sw(2:6) = ubar
    sw(7) = max(ur+sqrt(grav*hr+csqbar+pr),lamb7bar)

end subroutine riemann_normal_HypRel_BBBD

subroutine riemann_transverse_HypRel_BBBD(ixy,u,v,h,w,sg,p,pb, &
    grav,csq,localasdq,waves,s)
    
    implicit none
    !Input
    integer, intent(in) :: ixy
    double precision, intent(in) :: u, v, h, w, p, csq, localasdq(7,1)
    double precision, intent(in) :: sg, pb
    double precision, intent(in) :: grav

    !Output
    double precision, intent(out) ::  waves(7,7), s(7)

    !Local
    integer :: k
    double precision :: CE, lamb1, lamb7, R(7,7),  R_inv(7,7)
    double precision :: alpha(7,1)

    if (ixy.eq.1) then
        !Horizontal normal swipe (Construct R^y(u,v))
        CE=dsqrt(grav*h+csq+p)
        lamb1 = v-CE
        lamb7 = v+CE

        !Matrix of right eigenvectors
        R(1,:)= [1.d0,   1.d0,  0.d0, 0.d0, 0.d0, 0.d0, 1.d0]
        R(2,:)=[u,       0.d0,  0.d0, 0.d0, 1.d0, 0.d0, u]
        R(3,:)=[lamb1,   v,     0.d0, 0.d0, 0.d0, 0.d0, lamb7]
        R(4,:)= [w,      0.d0,  0.d0, 1.d0, 0.d0, 0.d0, w]
        R(5,:)= [sg,     0.d0,  1.d0, 0.d0, 0.d0, 0.d0, sg]
        R(6,:)= [p+csq,-grav*h, 0.d0, 0.d0, 0.d0, 0.d0, p+csq]
        R(7,:)= [pb,     0.d0,  0.d0, 0.d0, 0.d0, 1.d0, pb]

        R_inv(:,1)=(1.d0/CE**2.d0)*[0.5d0*(grav*h+CE*v),&
        csq+p, -grav*h*sg, -grav*h*w, -grav*h*u,&
        -grav*h*pb,0.5d0*(grav*h-CE*v)]
        R_inv(:,3)=[-0.5d0/CE, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.5d0/CE]
        R_inv(:,2)=[0.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 0.d0]
        R_inv(:,4)= [0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0]
        R_inv(:,5)= [0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0,  0.d0]
        R_inv(:,6)=(1.d0/CE**2.d0)*[0.5d0,-1.d0,-sg,&
            -w,-u, -pb, 0.5d0]
        R_inv(:,7)=[0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0]

        s(1) = lamb1
        s(2:6) = v
        s(7) = lamb7
    else
        !Vertical normal swipe (Construct R^x(u,v))
        CE=dsqrt(grav*h+csq+p)
        lamb1 = u-CE
        lamb7 = u+CE

        R(1,:)= [1.d0,   1.d0,   0.d0, 0.d0, 0.d0, 0.d0, 1.d0]
        R(2,:)=[lamb1,   u,      0.d0, 0.d0, 0.d0, 0.d0, lamb7]
        R(3,:)=[v,       0.d0,   0.d0, 0.d0, 1.d0, 0.d0, v]
        R(4,:)= [w,      0.d0,   0.d0, 1.d0, 0.d0, 0.d0, w]
        R(5,:)= [sg,     0.d0,   1.d0, 0.d0, 0.d0, 0.d0,  sg]
        R(6,:)= [pb+csq,-grav*h, 0.d0, 0.d0, 0.d0, 0.d0, p+csq]
        R(7,:)= [pb,     0.d0,   0.d0, 0.d0, 0.d0, 1.d0,  pb]


        R_inv(:,1)=(1.d0/CE**2.d0)*[0.5d0*(grav*h+CE*u),&
        csq+p, -grav*h*sg, -grav*h*w, -grav*h*v,&
        -grav*h*pb,0.5d0*(grav*h-CE*u)]
        R_inv(:,2)=[-0.5d0/CE, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.5d0/CE]
        R_inv(:,3)=[0.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 0.d0]
        R_inv(:,4)= [0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0]
        R_inv(:,5)= [0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0,  0.d0]
        R_inv(:,6)=(1.d0/CE**2.d0)*[0.5d0,-1.d0,-sg,&
            -w,-v, -pb, 0.5d0]
        R_inv(:,7)=[0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0]

        s(1) = lamb1
        s(2:6) = u
        s(7) = lamb7
    end if

    !Compute alphas
    alpha = matmul(R_inv,localasdq)
    do k=1,7
        waves(:,k) = alpha(k,1)*R(:,k)
    end do
end subroutine riemann_transverse_HypRel_BBBD

function Source_HypRel_BBBD(h,hw,hsg,hp,hpb,csq_space)
    implicit none
    double precision :: Source_HypRel_BBBD(4)
    double precision :: h, hw, hsg, hp, hpb, csq_space
    double precision :: w, sg, p, pb
    
    if (abs(h) > dry_tolerance) then
        w = hw/h
        sg = hsg/h
        p = hp/h
        pb = hpb/h
        Source_HypRel_BBBD(1) = pb
        Source_HypRel_BBBD(2) = -6.d0*pb+12.d0*p
        Source_HypRel_BBBD(3) = -2.d0*csq_space*sg
        Source_HypRel_BBBD(4) = -6.d0*csq_space*(w-0.5d0*sg)
    end if
end function 

function transition_function(s,s0,s1,optional_transition_type_fun,optional_shape) result(fs)
    !Returns f(s), where f(s0)=0 and f(s1)=1
    implicit none
    double precision, intent(in) :: s, s0, s1
    integer, optional, intent(in) :: optional_transition_type_fun
    double precision, optional, intent(in) :: optional_shape
    !Output
    double precision :: fs
    !Locals
    double precision :: shape = 1.d0 !Different effects depending on transition_type
    integer :: transition_type_fun = 1
    double precision :: sn !Normalized s (0<=sn<=1)
    !Edge cases
    if (s < s0) then
        fs = 0.d0
        return
    else if (s >= s1) then
        fs = 1.d0
        return
    end if
    !Normalize and get linear transition
    sn = (s-s0)/(s1-s0)
    if (present(optional_shape)) shape = optional_shape
    if (present(optional_transition_type_fun)) transition_type_fun = optional_transition_type_fun
    !Transition function
    select case(transition_type_fun)
        case(1) !Linear (can be made a Heaviside if s0==s1)
            fs = sn
        case(2) !Hyperbolic tangent
            !Larger shape -> sharper transition, better match at 0 and 1
            !Smaller shape -> smoother transition, does not match 0 and 1
            fs = 0.5*(tanh(10.d0*shape*(sn-0.5))+1)
        case(3) !Engsig-Karup (2007)
            !Matches 0 and 1, smooth in the middle, sharper at the ends
            fs = -2.d0*sn**3+3.d0*sn**2
        case(4) !Mayer et al. (1998)
            !Very smooth at 0, does not match 1
            !Smaller shape -> smoother transition at 0, better match at 1
            !and faster decrease close to 1
            fs = (0.01d0*shape*sn)**3+(1.d0-0.01d0*shape)*sn**6
        case(5) !Modified Mayer et al. (1998)
            !Very smooth at 1, does not match 0
            !Smaller shape -> smoother transition at 1, better match at 0
            !and faster decrease close to 0
            fs = 1.d0-((0.01d0*shape*(1.d0-sn))**3+(1.d0-0.01d0*shape)*(1.d0-sn)**6)
        case default
            fs = sn !Linear
    end select
end function transition_function

end module bouss_module
