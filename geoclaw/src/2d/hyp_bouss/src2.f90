subroutine src2(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux,t,dt)
      
    use geoclaw_module, only: g => grav, coriolis_forcing, coriolis
    use geoclaw_module, only: sea_level
    use geoclaw_module, only: friction_forcing, friction_depth
    use geoclaw_module, only: manning_coefficient
    use geoclaw_module, only: manning_break, num_manning
    use geoclaw_module, only: spherical_distance, coordinate_system
    use geoclaw_module, only: RAD2DEG, pi, dry_tolerance, DEG2RAD
    use geoclaw_module, only: rho_air, rho
    use geoclaw_module, only: earth_radius, sphere_source

    use geoclaw_module, only: speed_limit  ! used if no friction
      
    use storm_module, only: wind_forcing, pressure_forcing, wind_drag
    use storm_module, only: wind_index, pressure_index
    use storm_module, only: storm_direction, storm_location

    use friction_module, only: variable_friction, friction_index

    !Parameters for Hyperbolic Relaxation 
    use bouss_module, only: boussEquations
    use bouss_module, only: sdirk22_edc, sdirk22_bbbd !Local implicit source steps
    use bouss_module, only: gammaEDC => boussEDCgamma !Approximate SGN (3/2) or Sainte-Marie equations (2)
    use bouss_module, only: csq_index => bouss_csq_index !Index of the hyperbolic relaxation parameter
    use bouss_module, only: dr_index => boussDecayRate_index !Index of the decay rate
    use bouss_module, only: small_tol=> small_tol_bouss
    use bouss_module, only: large_tol=> large_tol_bouss

    implicit none
    
    ! Input parameters
    integer, intent(in) :: meqn,mbc,mx,my,maux
    double precision, intent(in) :: xlower,ylower,dx,dy,t,dt
    
    ! Output
    double precision, intent(inout) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
    double precision, intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

    ! Locals
    integer :: i, j, nman
    real(kind=8) :: h, hu, hv, gamma, dgamma, y, fdt, a(2,2), coeff
    real(kind=8) :: xm, xc, xp, ym, yc, yp, dx_meters, dy_meters
    real(kind=8) :: u, v, hu0, hv0
    real(kind=8) :: tau, wind_speed, theta, phi, psi, P_gradient(2), S(2)
    real(kind=8) :: Ddt, sloc(2)
    real(kind=8) :: tanyR, huv, huu, hvv
    real(kind=8) :: speed, sratio, xs, ys
    real(kind=8) :: eta

    ! Algorithm parameters

    ! Parameter controls when to zero out the momentum at a depth in the
    ! friction source term
    real(kind=8), parameter :: depth_tolerance = 1.0d-30

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Local EDC Hyperbolic Relaxation
    double precision :: bath !Local bathymetry
    double precision :: csq_space !Local hyperbolic relaxation parameter
    double precision :: hp, hw
    ! Local BBBD Hyperbolic Relaxation
    double precision :: hsg, hpb
    ! Source term substepping
    double precision :: dt_LSC !Linear stability condition for the source term
    double precision :: t_substep, dt_remainder
    ! Gradient of the bathymetry at cell centres (in metres)
    double precision :: bx, by, dx_m, dy_m, ycell
    double precision :: hu_dum, hv_dum
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! ----------------------------------------------------------------
    ! Spherical geometry source term(s)
    !
    ! These should be included for shallow water on the sphere, 
    ! at least in the mass term, but were only added in v5.9.0 as an option.
    ! rundata.geo_data.sphere_source can now be set in setrun.py
    ! Set sphere_source = 0 to omit source terms for backward compatibility
    ! sphere_source = 1 should become the default?

    if ((coordinate_system == 2) .and. (sphere_source > 0)) then
        ! add in spherical source term in mass equation 
        ! if sphere_source in [1,2],
        ! and also in momentum equations if sphere_source == 2
        do j=1,my
            y = ylower + (j - 0.5d0) * dy
            tanyR = tan(y*DEG2RAD) / earth_radius
            do i=1,mx
                if (q(1,i,j) > dry_tolerance) then
                    ! source term in mass equation:
                    q(1,i,j) = q(1,i,j) + dt * tanyR * q(3,i,j)

                    if (sphere_source == 2) then
                        ! Momentum source terms that drop out if linearized:
                        ! These seem to have very little effect for
                        ! practical problems
                        huv = q(2,i,j)*q(3,i,j)/q(1,i,j)
                        huu = q(2,i,j)*q(2,i,j)/q(1,i,j)
                        hvv = q(3,i,j)*q(3,i,j)/q(1,i,j)
                        q(2,i,j) = q(2,i,j) + dt * tanyR * 2.d0*huv
                        q(3,i,j) = q(3,i,j) + dt * tanyR * (hvv - huu)
                    endif
                endif
            enddo
        enddo
    endif
                
    ! ----------------------------------------------------------------                
    ! Friction source term
    if (friction_forcing) then
        do j=1,my
            do i=1,mx
                ! Extract appropriate momentum
                if (q(1,i,j) < depth_tolerance) then
                    q(2:3,i,j) = 0.d0
                else
                    ! Apply friction source term only if in shallower water
                    if (q(1,i,j) <= friction_depth) then
                        if (.not.variable_friction) then
                            do nman = num_manning, 1, -1
                                if (aux(1,i,j) .lt. manning_break(nman)) then
                                    coeff = manning_coefficient(nman)
                                endif
                            enddo
                        else
                            coeff = aux(friction_index,i,j)
                        endif
                        
                        ! Calculate source term
                        gamma = sqrt(q(2,i,j)**2 + q(3,i,j)**2) * g     &   
                              * coeff**2 / (q(1,i,j)**(7.d0/3.d0))
                        dgamma = 1.d0 + dt * gamma
                        q(2, i, j) = q(2, i, j) / dgamma
                        q(3, i, j) = q(3, i, j) / dgamma
                    endif
                endif
            enddo
        enddo
    ! End of friction source term

    else  

      ! if not friction_forcing: 

      ! Use speed_limit to scale back unphysical speeds that might
      ! arise from computing hu / h when h is very small.  This should
      ! only be needed in this routine if friction is not being used or
      ! the on-shore manning_coefficient is very small. 
      ! (In the latter case, the user might need to tweak this by hand.)

      do j=1-mbc,my+mbc
        do i=1-mbc,mx+mbc
            if (q(1,i,j) > 0.d0) then
                speed = sqrt((q(2,i,j)**2 + q(3,i,j)**2)) / q(1,i,j)
                if (speed > speed_limit) then
                    sratio = speed_limit / speed
                    q(2,i,j) = q(2,i,j) * sratio
                    q(3,i,j) = q(3,i,j) * sratio

                    if (.false.) then
                        ! write out info useful for investigating topo:
                        ! might give LOTS of output if no friction!
                        xs = xlower + (i-0.5d0)*dx
                        ys = ylower + (j-0.5d0)*dy
                        write(6,604) t, i,j, mx,my, dx
                        write(6,603) speed,q(1,i,j),aux(1,i,j),xs,ys

 604                    format('src2 at t =',f10.2, '  i,j,mx,my:',4i4, &
                               '  dx = ',f10.7)
 603                    format('     reset s =',f9.2,'  h=',e11.3, ' B=',f8.2,&
                               '  x,y = ', f11.6,',',f10.6)
                    endif
               endif
            endif
        enddo
      enddo

    ! end of speed limiting
    endif  ! end of else (i.e., if not friction_forcing)


    ! Coriolis source term
    ! TODO: May want to remove the internal calls to coriolis as this could 
    !       lead to slow downs.
    if (coriolis_forcing) then
        do j=1,my
            y = ylower + (j - 0.5d0) * dy
            fdt = coriolis(y) * dt ! Calculate f dependent on coordinate system

            ! Calculate matrix components
            a(1,1) = 1.d0 - 0.5d0 * fdt**2 + fdt**4 / 24.d0
            a(1,2) =  fdt - fdt**3 / 6.d0
            a(2,1) = -fdt + fdt**3 / 6.d0
            a(2,2) = a(1,1)

            do i=1,mx
                q(2,i,j) = q(2, i, j) * a(1,1) + q(3, i, j) * a(1,2)
                q(3,i,j) = q(2, i, j) * a(2,1) + q(3, i, j) * a(2,2)
            enddo
        enddo
    endif
    ! End of coriolis source term

    ! wind -----------------------------------------------------------
    if (wind_forcing) then
        ! Need storm location and direction for sector based wind drag
        sloc = storm_location(t)
        theta = storm_direction(t)
        do j=1,my
            yc = ylower + (j - 0.5d0) * dy
            do i=1,mx
                xc = xlower + (i - 0.5d0) * dx
                if (q(1,i,j) > dry_tolerance) then
                    psi = atan2(yc - sloc(2), xc - sloc(1))
                    if (theta > psi) then
                        phi = (2.d0 * pi - theta + psi) * RAD2DEG
                    else
                        phi = (psi - theta) * RAD2DEG 
                    endif
                    wind_speed = sqrt(aux(wind_index,i,j)**2        &
                                    + aux(wind_index+1,i,j)**2)
                    tau = wind_drag(wind_speed, phi) * rho_air * wind_speed / rho(1)
                    q(2,i,j) = q(2,i,j) + dt * tau * aux(wind_index,i,j)
                    q(3,i,j) = q(3,i,j) + dt * tau * aux(wind_index+1,i,j)
                endif
            enddo
        enddo
    endif
    ! ----------------------------------------------------------------

    ! Atmosphere Pressure --------------------------------------------
    ! Handled in Riemann solver
    ! if (pressure_forcing) then
    !     do j=1,my  
    !         ym = ylower + (j - 1.d0) * dy
    !         yc = ylower + (j - 0.5d0) * dy
    !         yp = ylower + j * dy
    !         do i=1,mx  
    !             xm = xlower + (i - 1.d0) * dx
    !             xc = xlower + (i - 0.5d0) * dx
    !             xp = xlower + i * dx
                
    !             if (coordinate_system == 2) then
    !                 ! Convert distance in lat-long to meters
    !                 dx_meters = spherical_distance(xp,yc,xm,yc)
    !                 dy_meters = spherical_distance(xc,yp,xc,ym)
    !             else
    !                 dx_meters = dx
    !                 dy_meters = dy
    !             endif

    !             ! Extract depths
    !             h = q(1,i,j)

    !             ! Calculate gradient of Pressure
    !             P_gradient(1) = (aux(pressure_index,i+1,j) &
    !                            - aux(pressure_index,i-1,j)) / (2.d0 * dx_meters)
    !             P_gradient(2) = (aux(pressure_index,i,j+1) &
    !                            - aux(pressure_index,i,j-1)) / (2.d0 * dy_meters)

    !                 ! Modify momentum in each layer
    !             if (h > dry_tolerance) then
    !                 q(2, i, j) = q(2, i, j) - dt * h * P_gradient(1) / rho(1)
    !                 q(3, i, j) = q(3, i, j) - dt * h * P_gradient(2) / rho(1)
    !             end if
    !         enddo
    !     enddo
    ! endif

    ! ----------------------------------------------------------------
    ! Hyperbolic relaxation systems (EDC = mbHSGN, BBBD = HSGN): relaxation terms and the
    ! non-hydrostatic bathymetry terms, which are not in the Riemann solver (no
    ! bathymetry sensor), integrated implicitly with one SDIRK(2,2) step of a local linear system
    ! per cell (closed form, h frozen). Stable for any dt and any bathymetry jump.
    !  - grad(b) at cell centres by centred differences of aux(1,:,:) (ghost cells are set by
    !    setaux, also on AMR patches), in metres: on lat-lon grids (coordinate_system == 2)
    !    b_x = db/(2 R dlambda cos(phi_j)), b_y = db/(2 R dphi), phi_j = cell-centre latitude.
    if (boussEquations == -1 .or. boussEquations == -2) then
        dy_m = dy
        if (coordinate_system == 2) dy_m = earth_radius*DEG2RAD*dy
        do j=1,my
            ycell = ylower + (j-0.5d0)*dy
            dx_m = dx
            if (coordinate_system == 2) dx_m = earth_radius*DEG2RAD*dx*cos(DEG2RAD*ycell)
            do i=1,mx
                if (aux(dr_index,i,j)>small_tol) then
                    h = q(1,i,j)
                    csq_space = aux(csq_index,i,j)
                    bx = (aux(1,i+1,j) - aux(1,i-1,j))/(2.d0*dx_m)
                    by = (aux(1,i,j+1) - aux(1,i,j-1))/(2.d0*dy_m)
                    if (boussEquations == -1) then
                        call one_step_SDIRK22_EDC(h, q(2,i,j), q(3,i,j), q(4,i,j), q(5,i,j), &
                                                  bx, by, csq_space, dt)
                    else
                        call one_step_SDIRK22_BBBD(h, q(2,i,j), q(3,i,j), q(4,i,j), q(5,i,j), &
                                                   q(6,i,j), q(7,i,j), bx, by, csq_space, dt)
                    endif
                end if
            enddo
        enddo
    endif

    if (boussEquations == -3) then
        ! DR: relaxation terms only (its bathymetry terms are in the Riemann solver);
        ! SDIRK(2,2) on (w,p), velocities untouched
        do j=1,my
            do i=1,mx
                if (aux(dr_index,i,j)>small_tol) then
                    h = q(1,i,j)
                    csq_space = aux(csq_index,i,j) 
                    hu_dum = q(2,i,j)
                    hv_dum = q(3,i,j)
                    call one_step_SDIRK22_EDC(h, hu_dum, hv_dum, q(4,i,j), q(5,i,j), &
                                              0.d0, 0.d0, csq_space, dt)
                end if
            enddo
        enddo
    endif


    contains

        subroutine one_step_SDIRK22_EDC(h, hu, hv, hw, hp, bx, by, csq_space, dt)
            !One SDIRK(2,2) step for y = (u,v,w,p) (h frozen):
            !  (u,v)_t = -gamma grad(b) p/h,  w_t = gamma p/h,
            !  p_t = (c^2 u.grad(b) - 2 c^2 w)/h
            !With grad(b) = 0 it is the SDIRK(2,2) step of the relaxation terms alone.
            implicit none
            double precision, intent(in) :: h, bx, by, csq_space, dt
            double precision, intent(inout) :: hu, hv, hw, hp
            double precision :: u, v, w, p

            if (h<dry_tolerance .or. csq_space < small_tol) then
                !Do nothing
                return
            end if
            u = hu/h
            v = hv/h
            w = hw/h
            p = hp/h
            call sdirk22_edc(h, u, v, w, p, bx, by, csq_space, gammaEDC, dt)
            hu = h*u
            hv = h*v
            hw = h*w
            hp = h*p
        end subroutine one_step_SDIRK22_EDC

        subroutine one_step_SDIRK22_BBBD(h, hu, hv, hw, hsg, hp, hpb, bx, by, csq_space, dt)
            !One SDIRK(2,2) step for y = (u,v,w,sg,p,pb) (h frozen):
            !  (u,v)_t = -grad(b) pb/h,  w_t = pb/h,  sg_t = (12 p - 6 pb)/h,
            !  p_t = (-c^2 sg - c^2 u.grad(b))/h,
            !  pb_t = (6 c^2 u.grad(b) - 6 c^2 w + 3 c^2 sg)/h
            implicit none
            double precision, intent(in) :: h, bx, by, csq_space, dt
            double precision, intent(inout) :: hu, hv, hw, hsg, hp, hpb
            double precision :: u, v, w, sg, p, pb

            if (h<dry_tolerance .or. csq_space < small_tol) then
                !Do nothing
                return
            end if
            u = hu/h
            v = hv/h
            w = hw/h
            sg = hsg/h
            p = hp/h
            pb = hpb/h
            call sdirk22_bbbd(h, u, v, w, sg, p, pb, bx, by, csq_space, dt)
            hu = h*u
            hv = h*v
            hw = h*w
            hsg = h*sg
            hp = h*p
            hpb = h*pb
        end subroutine one_step_SDIRK22_BBBD

end subroutine src2
