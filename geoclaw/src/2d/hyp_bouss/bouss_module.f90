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
    double precision :: boussEDCalpha      !Hyperbolic relaxation parameter alpha (EDC and BBBD):
                                           !c^2 = alpha^2 * g * H, H = local depth at rest (set in setaux)
                                           !EDC/BBBD: p-equation in non-conservative form, transmission
                                           !Riemann solvers (riemann_normal_HypRel_EDC/BBBD, riemann_transverse_HypRel_EDC_BBBD)
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
    read(iunit,*) boussEDCalpha     ! Hyperbolic relaxation parameter alpha (c^2 = alpha^2 g H)
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
    print *, "boussEDCalpha = ", boussEDCalpha
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
!======================================================================
! Normal Riemann solvers of EDC (mbHSGN) and BBBD (HSGN) with the p-equation in
! non-conservative form (c^2 outside the divergence), e.g. for EDC
!     (hp)_t + div(hu p) + c^2 [div(hu) - u.grad(eta)] = c^2 u.grad(b) - 2 c^2 w,
! with c^2 = c^2(x) taken from the aux array (alpha^2 g H, and 0 in cells where the SWEs are solved).
! Transmission-type solver: every wave carries the c^2 of the cell it enters (left-going acoustic wave
! c^2_L, right-going acoustic wave c^2_R, contact wave the c^2 of its upwind side), with speeds
! ubar -/+ sqrt(g hbar + pbar + c^2_side), and the total fluctuation is decomposed exactly into these
! f-waves. The p-update of a cell is then c^2_cell times the mass-flux jump entering it (the product
! c^2 div(hu) is unambiguous even where c^2 jumps). If the whole fan of the problem linearised with the
! c^2 of one of the cells lies inside that cell (ubar - C_R > 0, or ubar + C_L < 0), every wave carries
! that c^2, and the slow acoustic wave keeps its linearised speed (no Einfeldt bound), so that it cannot
! change cell.
! Bathymetry: only the hydrostatic term g*hbar*db (momentum) and the non-conservative product
! -c^2*ubar*d(eta) (p-equation) are in the Riemann solver; the non-hydrostatic bathymetry terms are
! integrated implicitly, together with the relaxation terms, in the source step (src2 / src1).
! Roe-type averages: hbar arithmetic; ubar, vbar, wbar, pbar (and sgbar, pbbar) sqrt(h)-weighted.
!======================================================================

!Escalante, Dumbser and Castro's system, 2019 (EDC)
subroutine riemann_normal_HypRel_EDC(mu, nv, grav, gamma, c_sq, &
                               hl, ul, vl, wl, pl, bl, csq_l,&
                               hr, ur, vr, wr, pr, br, csq_r, &
                               fw, sw)

   !Normal 2D Riemann solver for Escalante, Dumbser, and Castro's
   !hyperbolic relaxation of Sainte-Marie's dispersive water wave system
   !(transmission type, non-conservative p-equation: see the comment above)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !When used in 2D:
   !  Input: mu and nv according to the normal direction (analogous to Geoclaw's solver for SWEs)
   !  Output: 5 f-waves with 5 components each (fw) and 5 wave speeds (sw)
   !  Waves: 1 and 5 acoustic, 2 contact, 3 passive (w), 4 passive (transverse velocity)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !When used in 1D :
   !  Input: mu=2, nv =3, vl=vr=0
   !  Output: The 3rd component of all the f-waves must be discarded
   !          The 4th f-wave and the 4th wave speed must be discarded
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !gamma and c_sq are not used (kept for the interface)

   implicit none

   !Input
   integer :: mu, nv
   double precision :: grav, gamma, c_sq
   double precision :: hl, ul, vl, wl, pl, bl, csq_l
   double precision :: hr, ur, vr, wr, pr, br, csq_r

   !Output
   double precision :: fw(5,5), sw(5)

   !Local
   double precision :: hlsq, hrsq, etaL, etaR
   double precision :: hbar, ubar, vbar, wbar, pbar
   double precision :: cl, cr, lamb1, lamb5, csql, csqr, csqc, m, n
   double precision :: tf(5), beta(5)
   integer :: fan

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

   !Acoustic speeds with the c^2 of each side
   cl = dsqrt(grav*hbar + pbar + csq_l)
   cr = dsqrt(grav*hbar + pbar + csq_r)
   !Whole fan inside one cell: every wave carries the c^2 of that cell
   fan = 0
   csql = csq_l
   csqr = csq_r
   if (ubar - cr > 0.d0) then
      fan = 1
      csql = csq_r
      cl = cr
   else if (ubar + cl < 0.d0) then
      fan = -1
      csqr = csq_l
      cr = cl
   end if
   lamb1 = ubar - cl
   lamb5 = ubar + cr
   !Contact wave: c^2 of its upwind side
   csqc = csql
   if (ubar >= 0.d0) csqc = csqr

   !c^2-free parts of the total fluctuation
   tf(1) = hr*ur - hl*ul
   tf(2) = hr*ur**2 + hr*pr + 0.5d0*grav*hr**2 - hl*ul**2 - hl*pl - 0.5d0*grav*hl**2 + grav*hbar*(br-bl)
   tf(3) = hr*ur*vr - hl*ul*vl
   tf(4) = hr*ur*wr - hl*ul*wl
   tf(5) = hr*ur*pr - hl*ul*pl
   !Strengths: contact wave (it carries -c^2_c ubar d(eta)), acoustic waves, passive waves
   beta(2) = (pbar*tf(1) - tf(5) + csqc*ubar*(etaR-etaL))/(grav*hbar + pbar + csqc)
   m = tf(1) - beta(2)
   n = tf(2) - ubar*beta(2)
   beta(1) = (lamb5*m - n)/(lamb5 - lamb1)
   beta(5) = (n - lamb1*m)/(lamb5 - lamb1)
   beta(3) = tf(4) - wbar*(beta(1) + beta(5))
   beta(4) = tf(3) - vbar*(beta(1) + beta(5))

   ! Computing f-waves (they add up to the total fluctuation; the hp-row is
   ! d(hup) + c^2_L beta_1 + c^2_c beta_2 + c^2_R beta_5 - c^2_c ubar d(eta))
   fw = 0.d0
   fw(1,1) = beta(1)
   fw(mu,1) = beta(1)*lamb1
   fw(nv,1) = beta(1)*vbar
   fw(4,1) = beta(1)*wbar
   fw(5,1) = beta(1)*(pbar + csql)
   fw(1,2) = beta(2)
   fw(mu,2) = beta(2)*ubar
   fw(5,2) = -beta(2)*grav*hbar
   fw(4,3) = beta(3)
   fw(nv,4) = beta(4)
   fw(1,5) = beta(5)
   fw(mu,5) = beta(5)*lamb5
   fw(nv,5) = beta(5)*vbar
   fw(4,5) = beta(5)*wbar
   fw(5,5) = beta(5)*(pbar + csqr)

   !Computing wave speeds (Einfeldt-type bounds for the outer waves)
   sw(1) = min(ul - dsqrt(grav*hl + csq_l + pl), lamb1)
   if (fan == 1) sw(1) = lamb1
   sw(2) = ubar
   sw(3) = ubar
   sw(4) = ubar
   sw(5) = max(ur + dsqrt(grav*hr + csq_r + pr), lamb5)
   if (fan == -1) sw(5) = lamb5

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
   !Only used by the DR solver (riemann_normal_HypRel_DR). EDC and BBBD do not use a sensor:
   !their non-hydrostatic bathymetry terms are integrated implicitly in src2.

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

    !Normal 2D Riemann solver for Bassi, Bonaventura, Busto and Dumbser's hyperbolic relaxation
    !(transmission type, non-conservative p-equation: see the comment above riemann_normal_HypRel_EDC)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !When used in 2D:
    !  Input: mu and nv according to the normal direction (analogous to Geoclaw's solver for SWEs)
    !  Output: 7 f-waves with 7 components each (fw) and 7 wave speeds (sw)
    !  Waves: 1 and 7 acoustic, 2 contact, 3 passive (sigma), 4 passive (w),
    !         5 passive (transverse velocity), 6 passive (pb)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !When used in 1D :
    !  Input: mu=2, nv =3, vl=vr=0
    !  Output: The 3rd component of all the f-waves must be discarded
    !          The 5th f-wave and the 5th wave speed must be discarded
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !c_sq is not used (kept for the interface)

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
    double precision :: hbar, ubar, vbar, wbar, pbar, sgbar, pbbar
    double precision :: cl, cr, lamb1, lamb7, csql, csqr, csqc, m, n
    double precision :: tf(7), beta(7)
    integer :: fan

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

    !Acoustic speeds with the c^2 of each side
    cl = dsqrt(grav*hbar + pbar + csq_l)
    cr = dsqrt(grav*hbar + pbar + csq_r)
    !Whole fan inside one cell: every wave carries the c^2 of that cell
    fan = 0
    csql = csq_l
    csqr = csq_r
    if (ubar - cr > 0.d0) then
       fan = 1
       csql = csq_r
       cl = cr
    else if (ubar + cl < 0.d0) then
       fan = -1
       csqr = csq_l
       cr = cl
    end if
    lamb1 = ubar - cl
    lamb7 = ubar + cr
    !Contact wave: c^2 of its upwind side
    csqc = csql
    if (ubar >= 0.d0) csqc = csqr

    !c^2-free parts of the total fluctuation
    tf(1) = hr*ur - hl*ul
    tf(2) = hr*ur**2 + hr*pr + 0.5d0*grav*hr**2 - hl*ul**2 - hl*pl - 0.5d0*grav*hl**2 + grav*hbar*(br-bl)
    tf(3) = hr*ur*vr - hl*ul*vl
    tf(4) = hr*ur*wr - hl*ul*wl
    tf(5) = hr*ur*sgr - hl*ul*sgl
    tf(6) = hr*ur*pr - hl*ul*pl
    tf(7) = hr*ur*pbr - hl*ul*pbl
    !Strengths: contact wave (it carries -c^2_c ubar d(eta)), acoustic waves, passive waves
    beta(2) = (pbar*tf(1) - tf(6) + csqc*ubar*(etaR-etaL))/(grav*hbar + pbar + csqc)
    m = tf(1) - beta(2)
    n = tf(2) - ubar*beta(2)
    beta(1) = (lamb7*m - n)/(lamb7 - lamb1)
    beta(7) = (n - lamb1*m)/(lamb7 - lamb1)
    beta(3) = tf(5) - sgbar*(beta(1) + beta(7))
    beta(4) = tf(4) - wbar*(beta(1) + beta(7))
    beta(5) = tf(3) - vbar*(beta(1) + beta(7))
    beta(6) = tf(7) - pbbar*(beta(1) + beta(7))

    ! Computing f-waves (they add up to the total fluctuation; the hp-row is
    ! d(hup) + c^2_L beta_1 + c^2_c beta_2 + c^2_R beta_7 - c^2_c ubar d(eta))
    fw = 0.d0
    fw(1,1) = beta(1)
    fw(mu,1) = beta(1)*lamb1
    fw(nv,1) = beta(1)*vbar
    fw(4,1) = beta(1)*wbar
    fw(5,1) = beta(1)*sgbar
    fw(6,1) = beta(1)*(pbar + csql)
    fw(7,1) = beta(1)*pbbar
    fw(1,2) = beta(2)
    fw(mu,2) = beta(2)*ubar
    fw(6,2) = -beta(2)*grav*hbar
    fw(5,3) = beta(3)
    fw(4,4) = beta(4)
    fw(nv,5) = beta(5)
    fw(7,6) = beta(6)
    fw(1,7) = beta(7)
    fw(mu,7) = beta(7)*lamb7
    fw(nv,7) = beta(7)*vbar
    fw(4,7) = beta(7)*wbar
    fw(5,7) = beta(7)*sgbar
    fw(6,7) = beta(7)*(pbar + csqr)
    fw(7,7) = beta(7)*pbbar

    !Computing wave speeds (Einfeldt-type bounds for the outer waves)
    sw(1) = min(ul - dsqrt(grav*hl + csq_l + pl), lamb1)
    if (fan == 1) sw(1) = lamb1
    sw(2:6) = ubar
    sw(7) = max(ur + dsqrt(grav*hr + csq_r + pr), lamb7)
    if (fan == -1) sw(7) = lamb7

end subroutine riemann_normal_HypRel_BBBD

!======================================================================
! Transverse Riemann solver of EDC (mbHSGN, meqn = 5) and BBBD (HSGN, meqn = 7), transmission type,
! for the non-conservative p-equation (see riemann_normal_HypRel_EDC). The normal fluctuation asdq
! entering a cell (the middle cell, state q) is split, as in Clawpack's rpt2 for variable-coefficient
! acoustics, with one transmission problem per face of that cell:
!   bmasdq: what moves down in the problem at the bottom face (waves of the cell below | of the middle cell),
!   bpasdq: what moves up in the problem at the top face (waves of the middle cell | of the cell above);
! the waves that stay in the middle cell (reflected at the face) are dropped. Each acoustic wave carries
! P = p + c^2 of the cell it lives in (P = 0 in an SWE cell, dr = 0), with speed vt -/+ sqrt(g h + P),
! vt the velocity across the sweep; the contact and passive waves (speed vt) are lumped. If the whole fan
! of a face lies in one cell (as in the normal solver), all its waves use that cell's P. If the middle
! cell is SWE, the non-hydrostatic rows of asdq are dropped.
!======================================================================
subroutine transverse_split_HypRel(it, ip, meqn, q, gh, vt, plo, clo, phi, chi, a, dr_md, W, sp)
    ! Split of a (rows h, momentum across, hp) into a down-going acoustic wave carrying P_lo
    ! (speed vt - C_lo, W(:,1)), an up-going one carrying P_hi (vt + C_hi, W(:,3)) and the
    ! contact/passive remainder (vt, W(:,2))
    implicit none
    integer, intent(in) :: it, ip, meqn
    double precision, intent(in) :: q(meqn), gh, vt, plo, clo, phi, chi, a(meqn), dr_md
    double precision, intent(out) :: W(meqn,3), sp(3)
    double precision :: m, n, den, blo, bhi
    m = a(it) - vt*a(1)
    n = a(ip) + gh*a(1)
    den = clo*(gh + phi) + chi*(gh + plo)
    blo = (chi*n - m*(gh + phi))/den
    bhi = (clo*n + m*(gh + plo))/den
    W(:,1) = blo*q/q(1)
    W(it,1) = blo*(vt - clo)
    W(:,3) = bhi*q/q(1)
    W(it,3) = bhi*(vt + chi)
    if (dr_md <= small_tol_bouss) then
        W(4:meqn,1) = 0.d0
        W(4:meqn,3) = 0.d0
    end if
    W(ip,1) = blo*plo
    W(ip,3) = bhi*phi
    W(:,2) = a - W(:,1) - W(:,3)
    sp(1) = vt - clo
    sp(2) = vt
    sp(3) = vt + chi
end subroutine transverse_split_HypRel

subroutine riemann_transverse_HypRel_EDC_BBBD(ixy, meqn, q, grav, csq_dn, csq_md, csq_up, &
                                              dr_dn, dr_md, dr_up, asdq, bmasdq, bpasdq)
    ! Input: q = state of the middle cell; csq_* and dr_* = c^2 and SWE switch (aux(dr_index)) of the
    ! cells below (dn), in the middle (md) and above (up); asdq = fluctuation entering the middle cell.
    ! Output: the down/up-going flux differences (without the lat-lon factors dxdc).
    implicit none
    integer, intent(in) :: ixy, meqn
    double precision, intent(in) :: q(meqn), grav, csq_dn, csq_md, csq_up, dr_dn, dr_md, dr_up, asdq(meqn)
    double precision, intent(out) :: bmasdq(meqn), bpasdq(meqn)
    integer :: it, ip, k
    double precision :: a(meqn), gh, vt, p, pdn, pmd, pup, cdn, cmd, cup, W(meqn,3), sp(3)

    
    if (ixy == 2) then 
      it = 2
    else
      it = 3
    end if
    if (meqn == 7) then 
      ip = 6
    else
      ip = 5
    end if
    gh = grav*q(1)
    vt = q(it)/q(1)
    a = asdq
    p = q(ip)/q(1)
    if (dr_md <= small_tol_bouss) then
        a(4:meqn) = 0.d0
        p = 0.d0
    end if
    pdn = 0.d0
    pmd = 0.d0
    pup = 0.d0
    if (dr_dn > small_tol_bouss) pdn = p + csq_dn
    if (dr_md > small_tol_bouss) pmd = p + csq_md
    if (dr_up > small_tol_bouss) pup = p + csq_up
    cdn = dsqrt(max(gh + pdn, 1.d-2*gh))
    cmd = dsqrt(max(gh + pmd, 1.d-2*gh))
    cup = dsqrt(max(gh + pup, 1.d-2*gh))
    bmasdq = 0.d0
    bpasdq = 0.d0
    ! bottom face (below | middle); nothing goes down if the whole fan lies in the middle cell
    if (vt - cmd <= 0.d0) then
        if (vt + cdn < 0.d0) then     ! whole fan in the cell below
            call transverse_split_HypRel(it, ip, meqn, q, gh, vt, pdn, cdn, pdn, cdn, a, dr_md, W, sp)
        else
            call transverse_split_HypRel(it, ip, meqn, q, gh, vt, pdn, cdn, pmd, cmd, a, dr_md, W, sp)
        end if
        do k = 1, 3
            if (sp(k) < 0.d0) bmasdq = bmasdq + sp(k)*W(:,k)
        end do
    end if
    ! top face (middle | above); nothing goes up if the whole fan lies in the middle cell
    if (vt + cmd >= 0.d0) then
        if (vt - cup > 0.d0) then     ! whole fan in the cell above
            call transverse_split_HypRel(it, ip, meqn, q, gh, vt, pup, cup, pup, cup, a, dr_md, W, sp)
        else
            call transverse_split_HypRel(it, ip, meqn, q, gh, vt, pmd, cmd, pup, cup, a, dr_md, W, sp)
        end if
        do k = 1, 3
            if (sp(k) > 0.d0) bpasdq = bpasdq + sp(k)*W(:,k)
        end do
    end if
end subroutine riemann_transverse_HypRel_EDC_BBBD

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

!======================================================================
! local source systems: relaxation terms + non-hydrostatic bathymetry terms,
! one SDIRK(2,2) step in closed form, h frozen (called from src2 for every Boussinesq cell).
! EDC / mbHSGN, y = (u, v, w, p):
!     u_t = -gamma bx p/h,  v_t = -gamma by p/h,  w_t = gamma p/h,
!     p_t = (c^2 (u bx + v by) - 2 c^2 w)/h
! BBBD / HSGN, y = (u, v, w, sg, p, pb):
!     u_t = -bx pb/h,  v_t = -by pb/h,  w_t = pb/h,  sg_t = (12 p - 6 pb)/h,
!     p_t = (-c^2 sg - c^2 (u bx + v by))/h,
!     pb_t = (6 c^2 (u bx + v by) - 6 c^2 w + 3 c^2 sg)/h
!======================================================================
subroutine sdirk22_edc(h, u, v, w, p, bx, by, csq, gam, dt)
    implicit none
    double precision, intent(in) :: h, bx, by, csq, gam, dt
    double precision, intent(inout) :: u, v, w, p
    double precision :: a1, a2, bw, c1, c2, d, G, k, den
    double precision :: u1, v1, w1, p1, ru, rv, rw, rp
    a1 = -gam*bx/h
    a2 = -gam*by/h
    bw = gam/h
    c1 = csq*bx/h
    c2 = csq*by/h
    d = -2.d0*csq/h
    G = (1.d0 - 1.d0/dsqrt(2.d0))*dt
    k = dt - G
    den = 1.d0 - G*G*(c1*a1 + c2*a2 + d*bw)      ! = 1 - G^2 K >= 1 when K <= 0
    ! stage 1: (I - G A) y1 = y
    p1 = (p + G*(c1*u + c2*v + d*w))/den
    u1 = u + G*a1*p1
    v1 = v + G*a2*p1
    w1 = w + G*bw*p1
    ! stage 2: (I - G A) y2 = y + k A y1
    ru = u + k*a1*p1
    rv = v + k*a2*p1
    rw = w + k*bw*p1
    rp = p + k*(c1*u1 + c2*v1 + d*w1)
    p = (rp + G*(c1*ru + c2*rv + d*rw))/den
    u = ru + G*a1*p
    v = rv + G*a2*p
    w = rw + G*bw*p
end subroutine sdirk22_edc

subroutine solve_stage_bbbd(G, ax, ay, e, s1, s2, q1, qx, qy, mx, my, m2, m3, &
                            ru, rv, rw, rs, rp, rb, u, v, w, sg, p, pb)
    ! (I - G A) y = r for the BBBD matrix, by elimination
    implicit none
    double precision, intent(in) :: G, ax, ay, e, s1, s2, q1, qx, qy, mx, my, m2, m3
    double precision, intent(in) :: ru, rv, rw, rs, rp, rb
    double precision, intent(out) :: u, v, w, sg, p, pb
    double precision :: D1, P0, P1, S0, S1c, D2
    D1 = 1.d0 - G*G*q1*s1
    P0 = (rp + G*(qx*ru + qy*rv) + G*q1*rs)/D1
    P1 = G*G*(q1*s2 + qx*ax + qy*ay)/D1
    S0 = rs + G*s1*P0
    S1c = G*s1*P1 + G*s2
    D2 = 1.d0 - G*G*(mx*ax + my*ay + m2*e) - G*m3*S1c
    pb = (rb + G*(mx*ru + my*rv + m2*rw) + G*m3*S0)/D2
    p = P0 + P1*pb
    sg = S0 + S1c*pb
    u = ru + G*ax*pb
    v = rv + G*ay*pb
    w = rw + G*e*pb
end subroutine solve_stage_bbbd

subroutine sdirk22_bbbd(h, u, v, w, sg, p, pb, bx, by, csq, dt)
    implicit none
    double precision, intent(in) :: h, bx, by, csq, dt
    double precision, intent(inout) :: u, v, w, sg, p, pb
    double precision :: ax, ay, e, s1, s2, q1, qx, qy, mx, my, m2, m3, G, k
    double precision :: u1, v1, w1, sg1, p1, pb1
    ax = -bx/h
    ay = -by/h
    e = 1.d0/h
    s1 = 12.d0/h
    s2 = -6.d0/h
    q1 = -csq/h
    qx = -csq*bx/h
    qy = -csq*by/h
    mx = 6.d0*csq*bx/h
    my = 6.d0*csq*by/h
    m2 = -6.d0*csq/h
    m3 = 3.d0*csq/h
    G = (1.d0 - 1.d0/dsqrt(2.d0))*dt
    k = dt - G
    call solve_stage_bbbd(G, ax, ay, e, s1, s2, q1, qx, qy, mx, my, m2, m3, &
                          u, v, w, sg, p, pb, u1, v1, w1, sg1, p1, pb1)
    call solve_stage_bbbd(G, ax, ay, e, s1, s2, q1, qx, qy, mx, my, m2, m3, &
                          u + k*ax*pb1, v + k*ay*pb1, w + k*e*pb1, &
                          sg + k*(s1*p1 + s2*pb1), &
                          p + k*(q1*sg1 + qx*u1 + qy*v1), &
                          pb + k*(mx*u1 + my*v1 + m2*w1 + m3*sg1), &
                          u, v, w, sg, p, pb)
end subroutine sdirk22_bbbd

end module bouss_module
