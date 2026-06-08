! =====================================================
subroutine rpn2(ixy,maxm,meqn,mwaves,maux,mbc,mx,&
                      ql,qr,auxl,auxr,fwave,s,amdq,apdq)
! =====================================================
! Riemann solver for HSGN in 2D.
! Modified to handle mapped grids and inclusions with Neumann BCs.

! waves: 5 
! equations: 5

! Conserved quantities:
!       1   h (water height)
!       2   hu (x-momentum)
!       3   hv (y-momentum)
!       4   hw (related to vertical velocity)
!       5   hp (related to pressure)

! Auxiliary quantities:
!       1 a_x
!       2 a_y   where (a_x,a_y) is unit normal to left face
!       3 length_ratio_left ratio of length of left face to dyc
!       4 b_x
!       5 b_y   where (b_x,b_y) is unit normal to bottom face
!       6 length_ratio_bottom   ratio of length of bottom face to dxc
!       7 cell_area   ratio of cell area to dxc*dyc
!         (approximately Jacobian of mapping function)
!       inclusion_idx (if any) index for inclusion (0 if none)

! solve Riemann problems along one slice of data.

! On input, ql contains the state vector at the left edge of each cell
!           qr contains the state vector at the right edge of each cell

! This data is along a slice in the x-direction if ixy=1
!                            or the y-direction if ixy=2.
! On output, wave contains the waves,
!            s the speeds,
!            amdq the  left-going flux difference  A^- \Delta q
!            apdq the right-going flux difference  A^+ \Delta q


! Note that the i'th Riemann problem has left state qr(i-1,:)
!                                    and right state ql(i,:)
! From the basic clawpack routines, this routine is called with ql = qr
    !Parameters for Boussinesq equations
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !For the moment, the only system that requires changes in the RS
    !is EDC HypRel (boussEquations = -1)
    use cparam_mod, only: boussEquations
    use cparam_mod, only: c_sq => boussEDCcsq 
    !Reference hyperbolic relaxation parameter (Just used for EDC)
    !Approximate SGN (3/2) or Sainte-Marie equations (2)
    use cparam_mod, only: gamma => boussEDCgamma 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Mapped grid and geophysical parameters
    use cparam_mod, only:  g => grav!, sea_level
    use cparam_mod, only: inclusion_idx, inclusion_type
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    implicit none

    integer, intent(in) :: ixy, maxm, meqn, mwaves, maux, mbc, mx
    real(kind=8), intent(in) :: ql(meqn,1-mbc:maxm+mbc)
    real(kind=8), intent(in) :: qr(meqn,1-mbc:maxm+mbc)
    real(kind=8), intent(in) :: auxl(maux,1-mbc:maxm+mbc)
    real(kind=8), intent(in) :: auxr(maux,1-mbc:maxm+mbc)

    real(kind=8), intent(out) :: s(mwaves,1-mbc:maxm+mbc)
    real(kind=8), intent(out) :: fwave(meqn, mwaves,1-mbc:maxm+mbc)
    real(kind=8), intent(out) :: amdq(meqn,1-mbc:maxm+mbc)
    real(kind=8), intent(out) :: apdq(meqn,1-mbc:maxm+mbc)

    !Locals
    integer :: i, mw,m ,mu,nv

    ! double precision :: inside_inclusionL, inside_inclusionR

    !For mapped grid
    integer :: inx, iny, ilenrat
    double precision :: st, ct ! sin and cos of angle between grid lines and x-axis

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    double precision :: hL, hR, huL, huR, hvL, hvR
    ! double precision :: hwL, hwR, hpL, hpR
    double precision :: uL, uR, vL,vR, pL, pR
    double precision :: wL, wR !For EDC and BBBD HypRel
    ! double precision :: sgL, sgR, pbL, pbR !For BBBD HypRel
    double precision :: fwHypRel(meqn,mwaves), swHypRel(mwaves)
    ! double precision :: csq_l, csq_r, drL, drR
    double precision :: bL, bR
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    double precision :: inside_inclusionL, inside_inclusionR

    !No bathymetry for the moment
    bL = 0.d0
    bR = 0.d0


    if (ixy == 1) then
        !Riemann problem in the x-direction
        inx = 1 
        iny = 2 
        ilenrat = 3
        ! mu = 2 !Index for normal momentum component
        ! nv = 3 !Index for tangential momentum component
    else
        !Riemann problem in the y-direction
        inx = 4 
        iny = 5 
        ilenrat = 6
        ! mu = 3 !Index for normal momentum component
        ! nv = 2 !Index for tangential momentum component
    endif

    amdq(:,:)  = 0.d0
    apdq(:,:)  = 0.d0
    fwave(:,:,:) = 0.d0

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if ((trim(inclusion_type) == "neumann") .and. (inclusion_idx .gt. 0) ) then
        !Account for inclusion. I am splitting into two loops to avoid the overhead, but maybe there's
        !a better way to do this without repeating code
        do i = 2-mbc, mx+mbc
            !inform of a bad riemann problem from the start
            if((qr(1,i-1).lt.0.d0).or.(ql(1,i) .lt. 0.d0)) then
                write(*,*) 'Negative input: hl,hr,i=',qr(1,i-1),ql(1,i),i
            endif
            !Initialize Riemann problem for grid interface
            do mw=1,mwaves
                s(mw,i)=0.d0
                do m=1,meqn
                    fwave(m,mw,i)=0.d0
                enddo
            enddo

            ! Rotation matrix:
            !               [ cos(\theta) sin(\theta) ]
            !               [-sin(\theta) cos(\theta) ]
            ! Note that this reduces to the identity on a standard cartesian grid.
            !Get normal vector to the interface from aux array
            ct  = auxl(inx,i)
            st  = auxl(iny,i)

            inside_inclusionR = auxl(inclusion_idx,i)
            inside_inclusionL = auxr(inclusion_idx,i-1)

            if ( (inside_inclusionL .gt. 0.5d0) .and. (inside_inclusionR .gt. 0.5d0) ) then
                !Both cell centers inside inclusion: Do nothing at this interface
                cycle
            else if ( (inside_inclusionL .lt. 0.5d0) .and. (inside_inclusionR .lt. 0.5d0) ) then
                !Both cell centers outside inclusion: Proceed as usual
                hL = qr(1,i-1) 
                hR = ql(1,i)
                huL = ct*qr(2,i-1) + st*qr(3,i-1) !Normal velocity component
                huR = ct*ql(2,i) + st*ql(3,i) !Normal velocity component
                hvL = -st*qr(2,i-1) + ct*qr(3,i-1) !Tangential velocity component
                hvR = -st*ql(2,i) + ct*ql(3,i) !Tangential velocity component
                wL = qr(4,i-1)/hL
                wR = ql(4,i)/hR
                pL = qr(5,i-1)/hL
                pR = ql(5,i)/hR
            else if ( (inside_inclusionL .gt. 0.5d0) .and. (inside_inclusionR .lt. 0.5d0) ) then
                !Left cell center inside inlcusion, right cell center outside inclusion
                !constant extrapolation and reflect normal velocity component
                hR = ql(1,i)
                huR = ct*ql(2,i) + st*ql(3,i) !Normal velocity component
                wR = ql(4,i)/hR
                pR = ql(5,i)/hR
                hL = hR
                huL = -huR !Reflect normal velocity component
                wL = wR
                pL = pR
            else 
                !Left cell center outside inlcusion, right cell center inside inclusion
                !constant extrapolation and reflect normal velocity component
                hL = qr(1,i-1) 
                huL = ct*qr(2,i-1) + st*qr(3,i-1) !Normal velocity component
                wL = qr(4,i-1)/hL
                pL = qr(5,i-1)/hL
                hR = hL
                huR = -huL !Reflect normal velocity component
                wR = wL
                pR = pL
            end if

            if (boussEquations==-1 .or. boussEquations==0) then
                !EDC HypRel Boussinesq-type equations
                uL = huL/hL
                vL = hvL/hL
                uR=huR/hR
                vR=hvR/hR
                
                !Solve 2D Riemann problem in normal direction
                !We'll drop the tangential velocity component, and recover x,y components later
                call riemann_normal_HypRel_EDC(2,3,g,gamma, c_sq,&
                                    hL, uL, vL, wL, pL, bL, c_sq, &
                                    hR, uR, vR, wR, pR, bR, c_sq, &
                                    fwHypRel, swHypRel)     
                
                !Reconstruct f-waves and speeds in original 2D variables
                !rotating back to x,y components
                fwave(1,:,i) = fwHypRel(1,:)
                fwave(2,:,i) = ct*fwHypRel(2,:) - st*fwHypRel(3,:)
                fwave(3,:,i) = st*fwHypRel(2,:) + ct*fwHypRel(3,:)
                ! fwave(2,:,i) = fwHypRel(2,:)
                ! fwave(3,:,i) = fwHypRel(3,:)
                fwave(4,:,i) = fwHypRel(4,:)
                fwave(5,:,i) = fwHypRel(5,:)
                
                do mw=1,mwaves
                    fwave(:,mw,i) = auxl(ilenrat,i)*fwave(:,mw,i) !Scale by length ratio
                    s(mw,i) = auxl(ilenrat,i)*swHypRel(mw) !Scale by length ratio
                    ! s(mw,i) = swHypRel(mw) !Do not scale wave speeds
                end do
                do  mw=1,mwaves
                    if (s(mw,i) <-1.e-16) then
                        amdq(:,i) = amdq(:,i) + fwave(:,mw,i)
                    else if (s(mw,i)>1.e-16) then
                        apdq(:,i)  = apdq(:,i) + fwave(:,mw,i)
                    else
                        amdq(:,i) = amdq(:,i) + 0.5d0 * fwave(:,mw,i)
                        apdq(:,i) = apdq(:,i) + 0.5d0 * fwave(:,mw,i)
                    endif
                enddo
            else
                !Other Boussinesq-type equations (BBBD, DR)
                !To be implemented later
                write(*,*) 'Riemann solver for this Boussinesq system not yet implemented'
                stop
            end if
        end do
    else
        do i = 2-mbc, mx+mbc
            !inform of a bad riemann problem from the start
            if((qr(1,i-1).lt.0.d0).or.(ql(1,i) .lt. 0.d0)) then
                write(*,*) 'Negative input: hl,hr,i=',qr(1,i-1),ql(1,i),i
            endif
            !Initialize Riemann problem for grid interface
            do mw=1,mwaves
                s(mw,i)=0.d0
                do m=1,meqn
                    fwave(m,mw,i)=0.d0
                enddo
            enddo

            ! Rotation matrix:
            !               [ cos(\theta) sin(\theta) ]
            !               [-sin(\theta) cos(\theta) ]
            ! Note that this reduces to the identity on a standard cartesian grid.
            !Get normal vector to the interface from aux array
            ct  = auxl(inx,i)
            st  = auxl(iny,i)

            !Get data for Riemann Problem
            hL = qr(1,i-1) 
            hR = ql(1,i)
            huL = ct*qr(2,i-1) + st*qr(3,i-1) !Normal velocity component
            huR = ct*ql(2,i) + st*ql(3,i) !Normal velocity component
            hvL = -st*qr(2,i-1) + ct*qr(3,i-1) !Tangential velocity component
            hvR = -st*ql(2,i) + ct*ql(3,i) !Tangential velocity component

            if (boussEquations==-1 .or. boussEquations==0) then
                !EDC HypRel Boussinesq-type equations
                uL = huL/hL
                vL = hvL/hL
                uR=huR/hR
                vR=hvR/hR
                wL = qr(4,i-1)/hL
                wR = ql(4,i)/hR
                pL = qr(5,i-1)/hL
                pR = ql(5,i)/hR
                
                !Solve 2D Riemann problem in normal direction
                !We'll drop the tangential velocity component, and recover x,y components later
                call riemann_normal_HypRel_EDC(2,3,g,gamma, c_sq,&
                                    hL, uL, vL, wL, pL, bL, c_sq, &
                                    hR, uR, vR, wR, pR, bR, c_sq, &
                                    fwHypRel, swHypRel)     
                
                !Reconstruct f-waves and speeds in original 2D variables
                !rotating back to x,y components
                fwave(1,:,i) = fwHypRel(1,:)
                fwave(2,:,i) = ct*fwHypRel(2,:) - st*fwHypRel(3,:)
                fwave(3,:,i) = st*fwHypRel(2,:) + ct*fwHypRel(3,:)
                fwave(4,:,i) = fwHypRel(4,:)
                fwave(5,:,i) = fwHypRel(5,:)
                
                do mw=1,mwaves
                    fwave(:,mw,i) = auxl(ilenrat,i)*fwave(:,mw,i) !Scale by length ratio
                    s(mw,i) = auxl(ilenrat,i)*swHypRel(mw) !Scale by length ratio
                end do
                do  mw=1,mwaves
                    if (s(mw,i) <-1.e-16) then
                        amdq(:,i) = amdq(:,i) + fwave(:,mw,i)
                    else if (s(mw,i)>1.e-16) then
                        apdq(:,i)  = apdq(:,i) + fwave(:,mw,i)
                    else
                        amdq(:,i) = amdq(:,i) + 0.5d0 * fwave(:,mw,i)
                        apdq(:,i) = apdq(:,i) + 0.5d0 * fwave(:,mw,i)
                    endif
                enddo
            else
                !Other Boussinesq-type equations (BBBD, DR)
                !To be implemented later
                write(*,*) 'Riemann solver for this Boussinesq system not yet implemented'
                stop
            end if
        end do
    end if

    contains

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


end subroutine rpn2