module bouss_module
    use geoclaw_module, only:  dry_tolerance
    implicit none

    
    ! parameters set in setrun:
    integer :: boussEquations       ! Which equations? 0=SWE, 1=Madsen, 2=SGN, -1=EDC
    real(kind=8) :: boussMinDepth   ! depth below which to switch to SWE
    double precision :: boussEDCalpha      !Hyperbolic relaxation parameter alpha (EDC and BBBD):
                                           !c^2 = alpha^2 * g * H, H = local depth at rest (set in setaux)
                                           !EDC/BBBD: p-equation in non-conservative form, transmission
                                           !Riemann solvers (riemann_normal_HypRel_EDC/BBBD)
    double precision :: boussEDCgamma       !Use EDC to approximate SGN (3/2) or Sainte-Marie equations (2)

    ! Hyperbolic relaxation systems (boussEquations < 0): fields set by setaux in the aux array
    ! (aux(1) = topography, aux(2) = capacity), as in the 2D code:
    integer, parameter :: bouss_csq_index = 3       ! c^2 (0 where the SWEs are solved)
    integer, parameter :: boussDecayRate_index = 4  ! 1 where the dispersive system is solved, 0 for SWEs
    double precision, parameter :: small_tol_bouss = 1.d-13

    integer :: bc_xlo, bc_xhi
    real(kind=8) :: B_param, alpha

    ! for higher-order derivatives:
    real(kind=8), allocatable, dimension(:) :: rinv, dxm, dxc, cm2r, cp2r, c02r
    real(kind=8), allocatable, dimension(:) :: cm2, cp2, c02, cm1, cp1, c01

    ! for initial water depth and switching to SWE:
    real(kind=8) :: hmin
    real(kind=8), allocatable, dimension(:) :: h0
    logical, allocatable, dimension(:) :: useBouss
    
    ! for tridiagonal solver:
    integer, allocatable, dimension(:) :: IPIV
    real(kind=8), allocatable, dimension(:) :: D, DL, DU, DU2    

    save


contains

!======================================================================

    subroutine set_bouss(mx,mbc,mthbc)
    
    ! Set things up for Boussinesq solver, in particular
    ! create and factor tridiagonal matrix for implicit solves.

    use geoclaw_module, only: earth_radius,deg2rad,coordinate_system,sea_level
    use grid_module, only: xcell, grid_type
    use topo_module, only: zcell
    
    implicit none
    integer, intent(in) :: mx, mbc, mthbc(2)
    real(kind=8) :: DLi, r, denom
    integer :: i, iunit
    character*25 fname
    real(kind=8), dimension(1-mbc:mx+mbc) :: h02, h03
    integer(kind=4) :: INFO

    iunit = 7
    fname = 'bouss.data'
 !  # open the unit with new routine from Clawpack 4.4 to skip over
 !  # comment lines starting with #:
    call opendatafile(iunit, fname)

    read(7,*) boussEquations    ! which equations
    read(7,*) boussMinDepth     ! depth below which to switch to SWE
    read(7,*) boussEDCalpha     ! Hyperbolic relaxation parameter alpha (c^2 = alpha^2 g H)
    read(7,*) boussEDCgamma     ! Use EDC to approximate SGN (3/2) or Sainte-Marie equations (2)
    
    if (boussEquations==1) B_param = 1.d0/15.d0  ! parameter in MS equations
    if (boussEquations==2) then
        alpha = 1.d0!1.153  ! parameter in SGN equations
    end if
    
    allocate(rinv(0:mx+1),dxm(0:mx+2),dxc(0:mx+1))
    allocate(cm1(0:mx+1),cp1(0:mx+1),c01(0:mx+1))
    allocate(cm2(0:mx+1),cp2(0:mx+1),c02(0:mx+1))
    allocate(cm2r(0:mx+1),cp2r(0:mx+1),c02r(0:mx+1))
    allocate(D(mx+2), DL(mx+1), DU(mx+2), DU2(mx+2), IPIV(mx+2))
    allocate(h0(-1:mx+2), useBouss(-1:mx+2))


    ! Boundary conditions to impose in computing Boussinesq update:
    
    if ((mthbc(1)==2) .or. (mthbc(2)==2)) then
        write(6,*) '*** Periodic BCs not supported in bouss_module'
        stop
    endif

    ! Use Dirichlet BC with value 0 for correction in ghost cells by default:
    bc_xlo = 0
    bc_xhi = 0
    
    ! Check if wall BCs specified (also be used for radial symmetry at r=0):
    
    if (mthbc(1)==3) then
        ! For wall boundary conditions at left boundary:
        bc_xlo = 3
    endif

    if (mthbc(2)==3) then
        ! For wall boundary conditions at right boundary:
        bc_xhi = 3
    endif

    if (mthbc(1)==0) then
        ! For wavemaker BC at left
        bc_xlo = 3
    endif

    ! To try out Neumann BCs:
    ! This doesn't seem to work well, so not a general option
    !bc_xlo = 1
    !bc_xhi = 1
    
    
    !------------------------------------
    ! coefficients needed for second-order derivatives:

    do i=0,mx+2
        if (coordinate_system .eq. 2) then
            ! x is latitude, convert dx to meters:
            dxm(i) = (xcell(i) - xcell(i-1)) * earth_radius * deg2rad
        else
            ! x is meters
            dxm(i) = xcell(i) - xcell(i-1)
        endif
    enddo

    do i=0,mx+1
        if (grid_type == 0) then
            ! uniform grid
            dxc(i) = 2*dxm(i)
            cm2(i) = 1.d0/dxm(i)**2
            c02(i) = -2.d0/dxm(i)**2
            cp2(i) = 1.d0/dxm(i)**2
            ! first derivative coefficients:
            cm1(i) = -1.d0/dxc(i)
            c01(i) = 0.d0
            cp1(i) = 1.d0/dxc(i)
        else
            dxc(i) = (xcell(i+1) - xcell(i-1))  ! = 2*dx for uniform grid!
            cm2(i) = 2.d0 / (dxm(i)*dxc(i))
            cp2(i) = 2.d0 / (dxm(i+1)*dxc(i))
            c02(i) = -(cm2(i) + cp2(i))
            ! first derivative coefficients:
            denom = (dxm(i)*dxm(i+1)*(dxm(i)+dxm(i+1)))
            ! CORRECT??
            cm1(i) = -dxm(i+1)**2 / denom
            c01(i) = (dxm(i+1)**2 - dxm(i)**2) / denom
            cp1(i) = +dxm(i)**2 / denom
        endif
        
        if (coordinate_system == 1) then
            cm2r(i) = cm2(i)
            cp2r(i) = cp2(i)
            c02r(i) = c02(i)
            rinv(i) = 0.d0
        else if (coordinate_system == -1) then
            ! x = radial coordinate r >= 0
            ! include factors for ((1/r)*(r*q)_r)_r 
            ! using form q_{rr} + (1/r)q_r - (1/r**2)q
            r = xcell(i)
            rinv(i) = 1.d0/r
            cm2r(i) = cm2(i) - 1.d0/(r*dxc(i))
            cp2r(i) = cp2(i) + 1.d0/(r*dxc(i))
            c02r(i) = c02(i) - 1.d0/(r**2)
            ! first derivative coefficients for d/dr are unchanged 
            ! note: *not* including 1/r term in this derivative
        else if (coordinate_system == 2) then
            write(6,*) '*** latitude coordinates correctly implemented in Bouss??'
            stop
            
            ! PROBABLY NOT CORRECT -- r derivative??
            ! x = latitude coordinate -90 < x < 90 degrees
            ! r = sin(x) goes to zero at both poles
            ! include factors for ((1/r)*(r*q)_r)_r 
            ! using form q_{rr} + (1/r)q_r - (1/r**2)q
            r = sin(deg2rad*xcell(i))
            rinv(i) = 1.d0/r
            cm2r(i) = cm2(i) - 1.d0/(r*dxc(i))
            cp2r(i) = cp2(i) + 1.d0/(r*dxc(i))
            c02r(i) = c02(i) - 1.d0/(r**2)
            ! first derivative coefficients for d/dr are unchanged 
            ! note: *not* including 1/r term in this derivative        
        endif
    enddo

    ! initial depth and where to switch to SWE
    ! (useBouss is used by the MS and SGN solvers; the hyperbolic relaxation systems use the
    ! same criterion through aux(boussDecayRate_index,:), set in setaux):
    
    do i=1-mbc,mx+mbc
        h0(i) = max(sea_level - zcell(i), 0.d0)
    enddo
    
    do i=2-mbc,mx+mbc-1
        hmin = min(h0(i-1), h0(i), h0(i+1))
        useBouss(i) = (hmin > boussMinDepth)
        !write(6,*) '+++ useBouss: ',i,xcell(i),h0(i),useBouss(i)
    enddo        
    !Initialize first and last ghost cells
    useBouss(-1) = useBouss(0)
    useBouss(mx+2) = useBouss(mx+1)
    !------------------------------------
    
    if (boussEquations==1) then
    
        ! Madsen-Sorensen
        ! Form tridiagonal matrix and factor
        ! Done once at start and then used in each time step
        ! (for SGN, we need to build in each time step)

        call build_tridiag_ms(mx,mbc,mthbc)

    endif

    end subroutine set_bouss

!======================================================================
! Routines for Madsen-Sorensen equations
!======================================================================

subroutine build_tridiag_ms(mx,mbc,mthbc)

    use geoclaw_module, only: sea_level, coordinate_system
    use grid_module, only: xcell, grid_type
    use topo_module, only: zcell
    
    implicit none
    integer, intent(in) :: mx, mbc, mthbc(2)
    
    integer :: i
    real(kind=8) :: DLi, r
    real(kind=8), dimension(1-mbc:mx+mbc) :: h02, h03
    integer(kind=4) :: INFO
    
    !h0 = sea_level - aux(1,:)  # resting depth

    do i=1-mbc,mx+mbc
        !h0(i) = sea_level - zcell(i) ! now set in set_bouss
  666   format(i3,3e16.6)
        h02(i)=(max(0.,h0(i)))**2
        h03(i)=(max(0.,h0(i)))**3
    enddo

    ! initialize to identity matrix:
    D = 1.d0
    DU= 0.d0
    DL= 0.d0
    
    ! First and last rows (rows 1 and mx+2) of matrix corresponds to BCs
    ! Row i+1 of matrix corresponds to equation for cell i

    do i=1,mx
    
        ! Modify row i+1 (equation for i'th grid cell) to form (I-d^2) operator,
        ! unless this is a cell where SWE are used, then leave as row of I
        ! and set RHS to 0, so no modification to SWE result.
        
        if ((h0(i) > boussMinDepth)) then
            
          ! Replace this row of identity matrix with dispersion terms
          ! Note that cm2r(i)*w(i-1) + c02r(i)*w(i) + cp2r(i)*w(i+1) gives
          ! approximation to d^2 w = w_{xx} in plane wave case, or
          !                  d^2 w = (1/r * (r*w)_r)_r  in radial case
          ! Also h02 is h0**2 and h03 is h0**3:

          D(i+1) = 1.d0 - c02r(i)*((B_param+.5d0)*h02(i) &
                - 1.d0/6.d0*h03(i)/h0(i))

          DU(i+1)= -cp2r(i)*((B_param+.5d0)*h02(i) &
                - 1.d0/6.d0*h03(i)/h0(i+1))

          DL(i)= -cm2r(i)*((B_param+.5d0)*h02(i) &
                - 1.d0/6.d0*h03(i)/h0(i-1))

        !else
        !  write(66,*) 'i, h0(i): ',i,h0(i)
        endif
        
    enddo

    ! left boundary
    ! No change ==> Dirichlet value 0 in ghost cell, Q_0 = 0
    ! For other BCs, change superdiagonal in row 1:

    if (bc_xlo==1) then
        ! for Neumann BC at left: impose Q_0 = Q_1
        DU(1) = -1.d0     
    endif
    if (bc_xlo==3) then
        ! for wall-reflecting BC at left: impose Q_0 = -Q_1
        DU(1) = 1.d0     
    endif


    ! right boundary
    ! No change ==> Dirichlet value 0 in ghost cell, Q_{mx+1} = 0
    ! For other BCs, change subdiagonal in row mx+2:

    if (bc_xhi==1) then
        ! for Neumann at right: impose Q_mx = Q_{mx+1}
        DL(mx+1) = -1.d0
    endif
    if (bc_xhi==3) then
        ! for wall at right: impose Q_mx = - Q_{mx+1}
        DL(mx+1) = 1.d0
    endif

    
    if (.false.) then
        ! DEBUG:
        write(66,*) 'D matrix:'
        do i=1,mx+2
            if (i>1) then
                DLi = DL(i-1)
            else
                DLi = -9999.  ! since no subdiag in first row
            endif
            write(66,661) DLi, D(i), DU(i)
  661       format(3e16.6)
        enddo
        write(66,*) '===== end of D'
    endif

    ! factor the tridiagonal matrix:    
    call DGTTRF( mx+2, DL, D, DU, DU2, IPIV, INFO )

end subroutine build_tridiag_ms

!------------------------------------------------------------------

subroutine solve_tridiag_ms(mx,meqn,mbc,dx,q,maux,aux,psi)

    ! Set up right hand side for tridagonal system and then solve it
    ! psi first holds RHS and then the solution to be passed back.
    
    ! Note that matrix has already been factored in set_bouss, so the
    ! arrays D, DL, DU, DU2, IPIV are already set up for the solver.
    
    use geoclaw_module, only: g => grav, sea_level
    use geoclaw_module, only: coordinate_system
    use grid_module, only: xcell

    implicit none

    integer, intent(in) :: mx,meqn,mbc,maux
    real(kind=8), intent(in) :: dx

    real(kind=8), intent(in) ::  q(meqn,1-mbc:mx+mbc)
    real(kind=8), intent(in) :: aux(maux,1-mbc:mx+mbc)
    real(kind=8), intent(out) :: psi(mx+2)

    real(kind=8)  depth
    real(kind=8), dimension(1-mbc:mx+mbc) :: h0, h02, h03, eta, hu2
    real(kind=8), dimension(1-mbc:mx+mbc) :: h_eta_x,h0_eta_x,s1,s1_h0
    real(kind=8)  h0_eta_xxx,s1_xx,s1_h0_xx
    integer :: i,j,k,kk,iL,iR
    integer(kind=4) :: INFO
    integer :: LDB

    h0 = sea_level - aux(1,:)

    do i=1-mbc,mx+mbc
        ! compute hu2 = h*u**2
        if (q(1,i).gt.dry_tolerance) then
            hu2(i)= q(2,i)**2/q(1,i)
        else
            hu2(i)= 0.d0
        endif
        
        !hu2(i) = 0.d0  ! Debug: turning off nonlinear terms
        
        ! h0 is resting depth, h02 = h0**2 and h03 = h0**3:
        eta(i)= q(1,i) - h0(i)
        h02(i)=(max(0.d0,h0(i)))**2
        h03(i)=(max(0.d0,h0(i)))**3
        
    enddo

    ! compute these quantities, extending to one set of ghost cells:
    ! h_eta_x  = h * eta_x
    ! h0_eta_x = h0 * eta_x
    
    ! defaults to 0 if h0==0 at one of the stencil points:
    h_eta_x = 0.d0
    h0_eta_x = 0.d0
        
    do i= 0,mx+1
        if (minval(h0(i-1:i+1)) > 0.d0) then
            ! if h0 > 0 at all stencil points.  Note dxc(i) = 2*dx if uniform:
            h_eta_x(i)= q(1,i)*(eta(i+1)-eta(i-1))/dxc(i)
            h0_eta_x(i)= h0(i)*(eta(i+1)-eta(i-1))/dxc(i)
        endif
    enddo
     

    
    do i=0,mx+1
    
        ! compute approximation to
        !   s1 = (h*u**2)_x + g*h*eta_x   (plane wave) or
        !   s1 = (1/r)*(r*h*u**2)_r + r*g*h*eta_r   (radial)
        !      = (h*u**2)_r + (1/r)*h*u**2 + r*g*h*eta_r
    
        s1(i) =  (hu2(i+1) - hu2(i-1)) / dxc(i) + g*h_eta_x(i)
        !s1(i) =  g*h_eta_x(i)  ! debug
        !s1(i) =  (hu2(i+1) - hu2(i-1)) / dxc(i) 
                 
        if (coordinate_system == -1) then
            ! radial: add term for (1/r)*h*u**2:
            s1(i) =  s1(i) + hu2(i) / xcell(i) 
        else if (coordinate_system == 2) then
            ! latitude
            write(6,*) '*** latitude not yet implemented in bouss_module'
            stop
        endif
        
        ! compute s1_h0 = s1 / h0:
        
        if (h0(i) > dry_tolerance) then
            s1_h0(i)=s1(i)/h0(i)
        else
            s1_h0(i)=0.d0
        endif
    enddo
     

    ! Compute right hand side psi:
    
    psi = 0.d0

    do i=1,mx

       k = i+1  ! i'th equation corresponds to row k=i+1 of RHS

       if (h0(i) <= boussMinDepth) then
            ! no dispersive term:
            psi(k) = 0.d0

        else
            ! Note that cm2r(i)*w(i-1) + c02r(i)*w(i) + cp2r(i)*w(i+1) gives
            ! approximation to d^2 w = w_{xx} 
            ! or in radial case to d^2 w = (1/r * (rw)_r)_r

            ! compute s1_xx, approximates d^2 s1:

            s1_xx = cp2r(i)*s1(i+1) + c02r(i)*s1(i) + cm2r(i)*s1(i-1)

            ! compute s1_h0_xx, approximates d^2 s1_h0 = d^2(s1 / h0):
            
            s1_h0_xx = cp2r(i)*s1_h0(i+1) + c02r(i)*s1_h0(i) + cm2r(i)*s1_h0(i-1)
        
            ! compute h0_eta_xxx, approximates d^2 h0_eta_x:
            
            h0_eta_xxx = cp2r(i)*h0_eta_x(i+1) + c02r(i)*h0_eta_x(i) + &
                    cm2r(i)*h0_eta_x(i-1)
                    

            ! Right-hand side:
            
            psi(k) = -(B_param+.5d0) * h02(i) * s1_xx &
                        + h03(i)/6.d0 * s1_h0_xx &
                        + B_param * g * h02(i) * h0_eta_xxx
            if (.false.) then
                write(66,*) '+++i,rhs(i+1): ',i, psi(k)
                write(66,*) (B_param+.5d0) * h02(i), s1_xx
                write(66,*) h03(i)/6.d0,  s1_h0_xx 
                write(66,*) B_param * g * h02(i) * h0_eta_xxx
                write(66,*) B_param,  g,  h02(i),  h0_eta_xxx
                write(66,*) s1(i)
                write(66,*) h0_eta_x(i)
                write(66,*) eta(i)
            endif
              
                               
        endif


        ! OLD NOTATION, NOT USED NOW:
        !if ((h0(i) > sw_depth0) .and. (h0(i) < sw_depth1)) then
        !    ! reduce psi linearly between sw_depth0 and sw_depth1:
        !    psi(k) = psi(k) * (h0(i)-sw_depth0)/(sw_depth1-sw_depth0)
        !endif
            
    enddo
    
    if (.false.) then
        write(66,*) 'RHS:'
        do i=1,mx+2
            write(66,*) '+++i,rhs: ',i, psi(i)
        enddo
    endif
    
    ! -------------------------
    ! solve tridiagonal system.  
    ! psi is RHS on entry, solution on exit:
    
    LDB = mx+2  ! leading dimension of right-hand side matrix psi
    
    call DGTTRS( 'N', mx+2, 1, DL, D, DU, DU2, IPIV, psi, LDB, INFO )
                                
    if (.false.) then
        write(66,*) 'Solution psi: '
        do i=1,mx+2
            write(66,*) '+++ i,psi: ',i,psi(i)
        enddo
    endif

    return
end subroutine solve_tridiag_ms


!======================================================================
! Routines for SGN equations
!======================================================================

subroutine build_tridiag_sgn(meqn,mbc,mx,xlower,dx,q,maux,aux)

    use geoclaw_module, only: sea_level, coordinate_system
    use grid_module, only: xcell, grid_type
    use topo_module, only: zcell
    
    implicit none
    integer, intent(in) :: mbc,mx,meqn,maux
    real(kind=8), intent(in) :: xlower,dx
    real(kind=8), intent(in) ::  aux(maux,1-mbc:mx+mbc)
    real(kind=8), intent(inout) ::  q(meqn,1-mbc:mx+mbc)
    
    ! Locals
    integer :: i
    real(kind=8) :: DLi, r, Brr, denom, hmin
    real(kind=8), dimension(1-mbc:mx+mbc) :: h, B, Br, hr
    real(kind=8), dimension(1-mbc:mx+mbc) :: etar, eta, phi, w

    integer(kind=4) :: INFO

    do i=1-mbc,mx+mbc
        h(i) = q(1,i)  ! depth
        B(i) = aux(1,i)  ! topo
        eta(i) = B(i) + h(i)
    enddo
    
    ! reset useBouss(:) if it should be based on h rather than h0 or slope, etc:
    if (.false.) then
        do i=2-mbc,mx+mbc-1
            hmin = min(h(i-1), h(i), h(i+1))
            useBouss(i) = (h(i) > boussMinDepth)
        enddo
    endif
    
    do i=2-mbc,mx+mbc-1
        Br(i) = cm1(i)*B(i-1) + c01(i)*B(i) + cp1(i)*B(i+1)
        hr(i) = cm1(i)*h(i-1) + c01(i)*h(i) + cp1(i)*h(i+1)
        etar(i) = cm1(i)*eta(i-1) + c01(i)*eta(i) + cp1(i)*eta(i+1)
    enddo

    ! initialize to identity matrix:
    D = 1.d0
    DU= 0.d0
    DL= 0.d0
    
    ! First and last rows (rows 1 and mx+2) of matrix corresponds to BCs
    ! Row i+1 of matrix corresponds to equation for cell i

    do i=1,mx
    
        ! Modify row i+1 (equation for i'th grid cell) to form (I-d^2) operator,
        ! unless this is a cell where SWE are used, then leave as row of I
        ! and set RHS to 0, so no modification to SWE result.
        
        !if ((h(i) > boussMinDepth)) then
        if (useBouss(i)) then
            
          ! Replace this row of identity matrix with dispersion terms
          ! Note that cm2r(i)*w(i-1) + c02r(i)*w(i) + cp2r(i)*w(i+1) gives
          ! approximation to d^2 w = w_{xx} in plane wave case, or
          !                  d^2 w = (1/r * (r*w)_r)_r  in radial case
          ! while cm2, c02, cp2 are coefficients for pure second derivative
          ! even in radial case.
          
          Brr = cm2(i)*B(i-1) + c02(i)*B(i) + cp2(i)*B(i+1)
          
          D(i+1) = 1.d0 + alpha*(-h(i)**2/3.d0 *c02r(i) &
                    -h(i)*hr(i)*(c01(i) + rinv(i)) &
                    + 0.5d0*h(i)*(Brr - rinv(i)*Br(i)) + etar(i)*Br(i)) !CORRECT??
          
          DU(i+1) = alpha*(-h(i)**2/3.d0 *cp2r(i) - h(i)*hr(i)*cp1(i))

          DL(i) = alpha*(-h(i)**2/3.d0 *cm2r(i) - h(i)*hr(i)*cm1(i))
        endif
    enddo

    ! left boundary
    ! No change ==> Dirichlet value 0 in ghost cell, Q_0 = 0
    ! For other BCs, change superdiagonal in row 1:

    if (bc_xlo==1) then
        ! for Neumann BC at left: impose Q_0 = Q_1
        DU(1) = -1.d0     
    endif
    if (bc_xlo==3) then
        ! for wall-reflecting BC at left: impose Q_0 = -Q_1
        DU(1) = 1.d0     
    endif


    ! right boundary
    ! No change ==> Dirichlet value 0 in ghost cell, Q_{mx+1} = 0
    ! For other BCs, change subdiagonal in row mx+2:

    if (bc_xhi==1) then
        ! for Neumann at right: impose Q_mx = Q_{mx+1}
        DL(mx+1) = -1.d0
    endif
    if (bc_xhi==3) then
        ! for wall at right: impose Q_mx = - Q_{mx+1}
        DL(mx+1) = 1.d0
    endif

    
    if (.false.) then
        ! DEBUG:
        write(66,*) 'D matrix:'
        do i=1,mx+2
            if (i>1) then
                DLi = DL(i-1)
            else
                DLi = -9999.  ! since no subdiag in first row
            endif
            write(66,661) DLi, D(i), DU(i)
  661       format(3e16.6)
        enddo
        write(66,*) '===== end of D'
    endif

    ! factor the tridiagonal matrix:    
    call DGTTRF( mx+2, DL, D, DU, DU2, IPIV, INFO )

end subroutine build_tridiag_sgn

!------------------------------------------------------------------    

subroutine solve_tridiag_sgn(mx,meqn,mbc,dx,q,maux,aux,psi)

    ! For the SGN equations
    
    ! Set up right hand side for tridagonal system and then solve it
    ! psi first holds RHS and then the solution to be passed back.
    
    ! Note that matrix has already been factored in src1, so the
    ! arrays D, DL, DU, DU2, IPIV are already set up for the solver.
    
    use geoclaw_module, only: g => grav, sea_level
    use geoclaw_module, only: coordinate_system
    use grid_module, only: xcell

    implicit none

    integer, intent(in) :: mx,meqn,mbc,maux
    real(kind=8), intent(in) :: dx

    real(kind=8), intent(in) ::  q(meqn,1-mbc:mx+mbc)
    real(kind=8), intent(in) :: aux(maux,1-mbc:mx+mbc)
    real(kind=8), intent(out) :: psi(mx+2)

    real(kind=8), dimension(1-mbc:mx+mbc) :: h, eta, u, B, phi, w
    real(kind=8) :: Br, hr, etar, Brr, wr, phir, ur
    integer :: i,j,k,kk,iL,iR
    integer(kind=4) :: INFO
    integer :: LDB

    h = q(1,:) ! depth

    do i=1-mbc,mx+mbc
        if (q(1,i).gt.dry_tolerance) then
            u(i)= q(2,i)/q(1,i)
        else
            u(i)= 0.d0
        endif            

        B(i) = aux(1,i)  ! topo
        eta(i) = B(i) + h(i)
    enddo

        
    do i=0,mx+1
        Brr = cm2(i)*B(i-1) + c02(i)*B(i) + cp2(i)*B(i+1)  ! pure 2nd deriv
        w(i) = u(i)**2 * Brr
        ur = cm1(i)*u(i-1) + c01(i)*u(i) + cp1(i)*u(i+1)
        phi(i) = ur**2 + rinv(i)*ur*u(i) + (rinv(i)*ur)**2  ! CORRECT??
    enddo        
    

    ! Compute right hand side psi:
    
    psi = 0.d0

    do i=1,mx

       k = i+1  ! i'th equation corresponds to row k=i+1 of RHS

       !if (h(i) <= boussMinDepth) then
       if (.not. useBouss(i)) then
            ! no dispersive term:
            psi(k) = 0.d0

        else
            Br = cm1(i)*B(i-1) + c01(i)*B(i) + cp1(i)*B(i+1)   ! pure 1st deriv
            hr = cm1(i)*h(i-1) + c01(i)*h(i) + cp1(i)*h(i+1)
            etar = cm1(i)*eta(i-1) + c01(i)*eta(i) + cp1(i)*eta(i+1)
            wr = cm1(i)*w(i-1) + c01(i)*w(i) + cp1(i)*w(i+1)
            phir = cm1(i)*phi(i-1) + c01(i)*phi(i) + cp1(i)*phi(i+1)
                    

            ! Right-hand side:
            
            psi(k) = g/alpha * etar &
                    + 2*h(i)*(h(i)/3.d0 * phir + phi(i)*(hr + 0.5d0*Br)) &
                    + 0.5d0*h(i)*wr + w(i)*etar
            
            if (.false.) then
                write(66,*) '+++i,rhs(i+1): ',i, psi(k)
            endif                       
        endif


        ! OLD NOTATION, NOT USED NOW:
        !if ((h(i) > sw_depth0) .and. (h(i) < sw_depth1)) then
        !    ! reduce psi linearly between sw_depth0 and sw_depth1:
        !    psi(k) = psi(k) * (h(i)-sw_depth0)/(sw_depth1-sw_depth0)
        !endif

    enddo
    
    if (.false.) then
        write(66,*) 'RHS:'
        do i=1,mx+2
            write(66,*) '+++i,rhs: ',i, psi(i)
        enddo
    endif
    
    ! -------------------------
    ! solve tridiagonal system.  
    ! psi is RHS on entry, solution on exit:
    
    LDB = mx+2  ! leading dimension of right-hand side matrix psi
    
    call DGTTRS( 'N', mx+2, 1, DL, D, DU, DU2, IPIV, psi, LDB, INFO )
                                
    if (.false.) then
        write(66,*) 'Solution psi: '
        do i=1,mx+2
            write(66,*) '+++ i,psi: ',i,psi(i)
        enddo
    endif

    return
end subroutine solve_tridiag_sgn

!======================================================================
! Routines for Hyperbolic relaxations of dispersive water wave models
!======================================================================
!======================================================================
! Normal Riemann solvers of EDC (mbHSGN) and BBBD (HSGN) with the p-equation in
! NON-CONSERVATIVE form (c^2 outside the divergence), e.g. for EDC
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
   
   if (ubar >= 0.d0) then
      csqc = csqr
   else
      csqc = csql
   end if 

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
   if (fan == 1) then
      sw(1) = lamb1
   else
      sw(1) = min(ul - dsqrt(grav*hl + csq_l + pl), lamb1)
   end if
   sw(2:4) = ubar
   if (fan == -1) then
      sw(5) = lamb5
   else
      sw(5) = max(ur + dsqrt(grav*hr + csq_r + pr), lamb5)
   end if

end subroutine riemann_normal_HypRel_EDC

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
    if (ubar >= 0.d0) then
       csqc = csqr
    else
        csqc = csql
    end if

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
    if (fan == 1) then
        sw(1) = lamb1
    else
        sw(1) = min(ul - dsqrt(grav*hl + csq_l + pl), lamb1)
    end if
    sw(2:6) = ubar
    if (fan == -1) then
        sw(7) = lamb7
    else
        sw(7) = max(ur + dsqrt(grav*hr + csq_r + pr), lamb7)
    end if
end subroutine riemann_normal_HypRel_BBBD


!======================================================================
! Local source systems: relaxation terms + non-hydrostatic bathymetry terms,
! one SDIRK(2,2) step in closed form, h frozen (called from src1 for every Boussinesq cell).
! 2D form; in 1D call with by = 0, v = 0.
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
    ! (I - G A) y = r for the BBBD src3 matrix, by elimination
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
