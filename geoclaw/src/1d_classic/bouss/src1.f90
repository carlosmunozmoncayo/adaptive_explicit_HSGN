subroutine src1(meqn,mbc,mx,xlower,dx,q,maux,aux,t,dt)
    ! Called to update q by solving source term equation
    ! $q_t = \psi(q)$ over time dt starting at time t.
    !
    ! This default version integrates manning friction or other friction terms if present
    
    ! Version from 1D branch of GeoClaw originally from Dave George

    ! Also handles radial source term, if desired.
    ! Note: assumes radial about lower boundary, wall bc should be imposed

    
    use geoclaw_module, only: dry_tolerance, grav, DEG2RAD
    use geoclaw_module, only: friction_forcing
    use geoclaw_module, only: frictioncoeff => friction_coefficient
    use geoclaw_module, only: earth_radius, coordinate_system
    use grid_module, only: xcell
    use bouss_module, only: boussEquations, alpha, cm1, c01, cp1, &
                      solve_tridiag_ms, build_tridiag_sgn, &
                      solve_tridiag_sgn, useBouss
    !Parameters for Hyperbolic Relaxation EDC
    !Local relaxation parameter c^2 and SWE/Boussinesq indicator (aux, set in setaux)
    use bouss_module, only: csq_index => bouss_csq_index, dr_index => boussDecayRate_index
    use bouss_module, only: small_tol_bouss
    use bouss_module, only: gammaEDC => boussEDCgamma !Approximate SGN (3/2) or Sainte-Marie equations (2)
    use bouss_module, only: sdirk22_edc, sdirk22_bbbd !Local implicit source steps

    implicit none
    integer, intent(in) :: mbc,mx,meqn,maux
    real(kind=8), intent(in) :: xlower,dx,t,dt
    real(kind=8), intent(in) ::  aux(maux,1-mbc:mx+mbc)
    real(kind=8), intent(inout) ::  q(meqn,1-mbc:mx+mbc)

    ! Locals
    real(kind=8) :: eta(0:mx+1)
    real(kind=8) :: gamma, rcell, tanxR, etax
    real(kind=8) :: rk_stage(1:mx,4), delt
    integer ::  i,k,ii,rk_order
    real(kind=8)  q0(meqn,1-mbc:mx+mbc)
    real(kind=8) psi(mx+2)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Local EDC Hyperbolic Relaxation
    !The 2 above should be passed by bouss_module
    double precision, dimension(4,2) :: rk_EDC
    double precision :: bath !Local bathymetry
    double precision :: csq_space !Local hyperbolic relaxation parameter
    double precision :: h, hp, hu, hw, u, p, w
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Local BBBD
    double precision :: hwout, hsgout, hpout, hpbout
    !Gradients of the bathymetry and of c^2 at cell centres
    double precision :: bx

    if (frictioncoeff.gt.0.d0 .and. friction_forcing) then
          ! integrate source term based on Manning formula
            do i=1,mx
               if (q(1,i)<=dry_tolerance) then
                  q(2,i) = 0.0
               else
                  gamma= dsqrt(q(2,i)**2)*(grav*frictioncoeff**2)/(q(1,i)**(7.0/3.0))
                  q(2,i)= q(2,i)/(1.d0 + dt*gamma)
              endif
            enddo
    endif

!      ----------------------------------------------------------------
!Forcing terms for radially symmetric equations
    if (coordinate_system == -1 .and. boussEquations>=0) then
        ! radial source term for SWE:
        do i=1,mx
            if (q(1,i) .gt. dry_tolerance) then
                ! x is radial coordinate in meters, x>=0, 
                ! u = radial velocity
                q(1,i) = q(1,i) - dt/xcell(i) * q(2,i)
                u = q(2,i)/q(1,i)
                q(2,i) = q(2,i) - dt/xcell(i) * q(1,i)*u**2
            end if
        end do
    end if

    if (coordinate_system == -1 .and. boussEquations==-1) then
        do i=1,mx
            if (q(1,i) .gt. dry_tolerance) then
                ! x is radial coordinate in meters, x>=0, 
                ! u = radial velocity
                h = q(1,i)
                q(1,i) = h - dt/xcell(i) * q(2,i)
                u = q(2,i)/h
                q(2,i) = q(2,i) - dt/xcell(i) * h*u**2
                if (aux(dr_index,i) > small_tol_bouss) then
                    q(3,i) = q(3,i) - dt/xcell(i) * u*q(3,i)
                    q(4,i) = q(4,i) - dt/xcell(i) * (u*q(4,i)+h*u*aux(csq_index,i))
                end if
            end if
        enddo
    end if

    if (coordinate_system == -1 .and. boussEquations==-2) then
        do i=1,mx
            if (q(1,i) .gt. dry_tolerance) then
                ! x is radial coordinate in meters, x>=0, 
                ! u = radial velocity
                h = q(1,i)
                q(1,i) = h - dt/xcell(i) * q(2,i)
                u = q(2,i)/h
                q(2,i) = q(2,i) - dt/xcell(i) * h*u**2
                if (aux(dr_index,i) > small_tol_bouss) then
                    q(3,i) = q(3,i) - dt/xcell(i) * u*q(3,i)
                    q(4,i) = q(4,i) - dt/xcell(i) * u*q(4,i)
                    q(5,i) = q(5,i) - dt/xcell(i) * (u*q(5,i)+h*u*aux(csq_index,i))
                    q(6,i) = q(6,i) - dt/xcell(i) * u*q(6,i)
                end if
            end if
        enddo
    end if

!      ----------------------------------------------------------------

    if (coordinate_system == 2) then
        ! source term for x = latitude in degrees -90 <= x <= 90,
        ! u = velocity in latitude direction (m/s) on sphere:
        do i=1,mx
            if (q(1,i) .gt. dry_tolerance) then
                tanxR = tan(xcell(i)*DEG2RAD) / earth_radius
                q(1,i) = q(1,i) + dt * tanxR * q(2,i)
                u = q(2,i)/q(1,i)
                q(2,i) = q(2,i) + dt * tanxR * q(1,i)*u**2
            endif
         enddo
     endif

    ! -------------------------------------------------
    if (boussEquations > 0) then
    ! Boussinesq terms for MS and SGN 

        rk_order = 1  ! 1 for Forward Euler, 2 for second-order RK
        delt = dt / rk_order

        if (boussEquations == 2) then
            ! for SGN, need to factor matrix each step
            call build_tridiag_sgn(meqn,mbc,mx,xlower,dx,q,maux,aux)
        endif
                
        q0  = q
          
        do i=0,mx+1
           if (q(1,i) > dry_tolerance) then
               eta(i) = q(1,i)+aux(1,i)
           else
               eta(i) = 0.d0
           endif
        enddo

        !-----------------------
        ! First stage (only stage for rk_order == 1):
        
        if (boussEquations == 1) then
          call solve_tridiag_ms(mx,meqn,mbc,dx,q0,maux,aux,psi)
          ! returns solution psi = source term for Madsen
        else if (boussEquations == 2) then
          call solve_tridiag_sgn(mx,meqn,mbc,dx,q0,maux,aux,psi)
          ! modify solution psi for source term of SGN:
          do i=1,mx
              if (useBouss(i)) then
                  etax = cm1(i)*eta(i-1) + c01(i)*eta(i) + cp1(i)*eta(i+1)
                  psi(i+1) = q(1,i) * (grav/alpha * etax - psi(i+1))
              else
                  psi(i+1) = 0.d0
              endif
           enddo
        endif
              
        ! Forward Euler update to momentum q(2,:):
        ! Note psi(1) used for BC, so psi(i+1) updates q(2,i):
        q0(2,1:mx) = q0(2,1:mx) + delt*psi(2:mx+1)
          

        !-----------------------
        if (rk_order == 1) then
            q = q0 ! and we are done
              
        else if (rk_order == 2) then

            ! Second stage for 2-stage R-K, solve for psi at midpoint
            ! in time based on q0 computed in first stage:
            
            if (boussEquations == 1) then
                call solve_tridiag_ms(mx,meqn,mbc,dx,q0,maux,aux,psi)
                ! returns solution psi = source term for Madsen
            else if (boussEquations == 2) then
                call solve_tridiag_sgn(mx,meqn,mbc,dx,q0,maux,aux,psi)
                ! modify solution psi for source term of SGN:
                do i=1,mx
                    ! note that h,eta were not changed by first stage
                    if (useBouss(i)) then
                        etax = cm1(i)*eta(i-1) + c01(i)*eta(i) + cp1(i)*eta(i+1)
                        psi(i+1) = q(1,i) * (grav/alpha * etax - psi(i+1))
                    else
                        psi(i+1) = 0.d0
                    endif
                enddo
            endif
              
            ! Second stage is midpoint method dt=delt*2:
            ! Note psi(1) used for BC, so psi(i+1) updates q(2,i):
            q(2,1:mx) = q(2,1:mx) + 2.d0*delt*psi(2:mx+1)

        endif ! rk_order == 2

    endif ! end of Bouss terms

    ! -------------------------------------------------
    !Hyperbolic relaxation systems (EDC = mbHSGN, BBBD = HSGN): relaxation terms AND the
    !non-hydrostatic bathymetry terms, which are not in the Riemann solver,
    !integrated implicitly with one SDIRK(2,2) step of a local linear system
    !per cell (closed form, h frozen). Stable for any dt and any bathymetry jump.
    !b_x at cell centres from the physical cell centres xcell (non-uniform grids): the 3-point
    !formula of bouss_module (cm1, c01, cp1; = centred difference on uniform grids).
    !c^2 = aux(csq_index), the cell value (alpha^2 g H, set in setaux; 0 where the SWEs are solved).
    !The p-equation is in non-conservative form (c^2 outside the divergence), so there is no c^2_x
    !term: c^2 (hu)_x is in the (transmission-type) Riemann solver, which also handles the jump of
    !c^2 to 0 at the SWE/Boussinesq switch.
    !Only where the dispersive system is solved (aux(dr_index) = 1, as useBouss).
    if (boussEquations == -1 .or. boussEquations == -2) then
        do i=max(2-mbc,0),min(mx+mbc-1,mx+1) !i=1,mx+mbc
            if (aux(dr_index,i) > small_tol_bouss) then
                bx = cm1(i)*aux(1,i-1) + c01(i)*aux(1,i) + cp1(i)*aux(1,i+1)
                if (boussEquations == -1) then
                    call one_step_SDIRK22_EDC_1D(q(1,i), q(2,i), q(3,i), q(4,i), bx, &
                                              aux(csq_index,i), dt)
                else
                    call one_step_SDIRK22_BBBD_1D(q(1,i), q(2,i), q(3,i), q(4,i), q(5,i), q(6,i),&
                                               bx, aux(csq_index,i), dt)
                endif
            endif
        end do
    end if

    contains

        subroutine one_step_SDIRK22_EDC_1D(h, hu, hw, hp, bx, csq_space, dt)
            !One SDIRK(2,2) step for y = (u,w,p) (h frozen):
            !  u_t = -gamma b_x p/h,  w_t = gamma p/h,  p_t = (c^2 b_x u - 2 c^2 w)/h
            implicit none
            double precision, intent(in) :: h, bx, csq_space, dt
            double precision, intent(inout) :: hu, hw, hp
            double precision :: u, v, w, p

            if (h<dry_tolerance .or. csq_space < 1.e-6) then
                !Do nothing
                return
            end if
            u = hu/h
            v = 0.d0
            w = hw/h
            p = hp/h
            call sdirk22_edc(h, u, v, w, p, bx, 0.d0, csq_space, gammaEDC, dt)
            hu = h*u
            hw = h*w
            hp = h*p
        end subroutine one_step_SDIRK22_EDC_1D

        subroutine one_step_SDIRK22_BBBD_1D(h, hu, hw, hsg, hp, hpb, bx, csq_space, dt)
            !One SDIRK(2,2) step for y = (u,w,sg,p,pb) (h frozen):
            !  u_t = -b_x pb/h,  w_t = pb/h,  sg_t = (12 p - 6 pb)/h,
            !  p_t = (-c^2 sg - c^2 b_x u)/h,  pb_t = (6 c^2 b_x u - 6 c^2 w + 3 c^2 sg)/h
            implicit none
            double precision, intent(in) :: h, bx, csq_space, dt
            double precision, intent(inout) :: hu, hw, hsg, hp, hpb
            double precision :: u, v, w, sg, p, pb

            if (h<dry_tolerance .or. csq_space < 1.e-6) then
                !Do nothing
                return
            end if
            u = hu/h
            v = 0.d0
            w = hw/h
            sg = hsg/h
            p = hp/h
            pb = hpb/h
            call sdirk22_bbbd(h, u, v, w, sg, p, pb, bx, 0.d0, csq_space, dt)
            hu = h*u
            hw = h*w
            hsg = h*sg
            hp = h*p
            hpb = h*pb
        end subroutine one_step_SDIRK22_BBBD_1D

end subroutine src1

