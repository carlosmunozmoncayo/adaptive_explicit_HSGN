! ============================================
subroutine b4step2(mbc, mx, my, meqn, q, xlower, ylower, dx, dy, t, dt,        &
                   maux, aux, actualstep)
! ============================================
!
! # called before each call to step
! # use to set time-dependent aux arrays or perform other tasks.
!
! This particular routine sets negative values of q(1,i,j) to zero,
! as well as the corresponding q(m,i,j) for m=1,meqn.
! This is for problems where q(1,i,j) is a depth.
! This should occur only because of rounding error.
!
! Also calls movetopo if topography might be moving.

    use geoclaw_module, only: dry_tolerance
    use geoclaw_module, only: g => grav
    use geoclaw_module, only: speed_limit, coordinate_system
    use topo_module, only: num_dtopo,topotime
    use topo_module, only: aux_finalized
    use topo_module, only: xlowdtopo,xhidtopo,ylowdtopo,yhidtopo

    use amr_module, only: xlowdomain => xlower
    use amr_module, only: ylowdomain => ylower
    use amr_module, only: xhidomain => xupper
    use amr_module, only: yhidomain => yupper
    use amr_module, only: xperdom,yperdom,spheredom,NEEDS_TO_BE_SET
    use amr_module, only: outunit

    use storm_module, only: set_storm_fields

    use bouss_module, only: c_sq_index => bouss_csq_index
    use bouss_module, only: dr_index => boussDecayRate_index
    use bouss_module, only: boussTransitionType
    use bouss_module, only: boussEquations
    use bouss_module, only: boussTstart, boussTfinal
    use bouss_module, only: boussTransitionType
    use bouss_module, only: use_bouss
    use bouss_module, only: small_tol_bouss, large_tol_bouss



    implicit none

    ! Subroutine arguments
    integer, intent(in) :: meqn
    integer, intent(inout) :: mbc,mx,my,maux
    real(kind=8), intent(inout) :: xlower, ylower, dx, dy, t, dt
    real(kind=8), intent(inout) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
    real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
    logical, intent(in) :: actualstep

    ! Local storage
    integer :: index,i,j,k,dummy
    real(kind=8) :: h,u,v,sratio, xs,ys, s
    integer :: neighborhood_idx(8,2)
    integer :: ii,jj,kk,iidx,jidx,indicator_neighborhood
    real(kind=8) :: detadx_limited, detady_limited
    real(kind=8) :: etaC, etaL, etaR, etaB, etaT


    ! Check for NaNs in the solution
    call check4nans(meqn,mbc,mx,my,q,t,1)

    ! check for h < 0 and reset to zero
    ! check for h < drytolerance
    ! set hu = hv = 0 in all these cells
    forall(i=1-mbc:mx+mbc, j=1-mbc:my+mbc,q(1,i,j) < dry_tolerance)
        q(1,i,j) = max(q(1,i,j),0.d0)
        q(1:meqn,i,j) = 0.d0
    end forall

    ! Check for fluid speed sqrt(u**2 + v**2) > speed_limit
    ! and reset by scaling (u,v) down to this value (preserving direction)
    ! Note: similar check is done in getmaxspeed
    ! This helps avoid too many dt reductions when flow off very steep topo
    ! with delta B larger than fluid depth gives big speeds in Riemann solution
    ! (shallow water equations aren't valid for flow off a cliff)
    do j=1-mbc,my+mbc
        do i=1-mbc,mx+mbc
            if (q(1,i,j) > 0.d0) then
                s = sqrt((q(2,i,j)**2 + q(3,i,j)**2)) / q(1,i,j)
                if (s > speed_limit) then
                    sratio = speed_limit / s
                    q(2,i,j) = q(2,i,j) * sratio
                    q(3,i,j) = q(3,i,j) * sratio

                    if (.true.) then
                        ! write out info useful for investigating topo:
                        xs = xlower + (i-0.5d0)*dx
                        ys = ylower + (j-0.5d0)*dy
                        write(6,604) t, i,j, mx,my, dx
                        write(outunit,604) t, i,j, mx,my, dx
                        write(6,603) s,q(1,i,j),aux(1,i,j),xs,ys
                        write(outunit,603) s,q(1,i,j),aux(1,i,j),xs,ys

 604                    format('b4step2 at t =',f10.2, '  i,j,mx,my:',4i4, &
                               '  dx = ',f10.7)
 603                    format('     reset s =',f9.2,'  h=',e11.3, ' B=',f8.2,&
                               '  x,y = ', f11.6,',',f10.6)
                    endif
               endif
            endif
        enddo
    enddo

    !-------------------------------------------------------
    !Turning on/off dispersive terms
    !Uses the fact that we are hitting 
    !BoussTstart and BoussTfinal exactly (up to machine precision)
    !in tick.f.
    if ((abs(boussTstart-t)<=small_tol_bouss .or. abs(boussTfinal-t)<=small_tol_bouss) &
        .and. actualstep) then
        !Reset aux array to turn on/off dispersive terms
        ! aux(dr_index,:,:) = NEEDS_TO_BE_SET !Always recompute this array to avoid undefined values
        call setaux(mbc,mx,my,xlower,ylower,dx,dy,maux,aux)
        if (abs(boussTstart-t)<=small_tol_bouss) then
            print *, '------------------------------------------------'
            print *, 'Turning on dispersive terms at t=', t
            print *, "USEBOUSS = ", use_bouss
            print *, '------------------------------------------------'
        else
            print *, '------------------------------------------------'
            print *, 'Turning off dispersive terms at t=', t
            print *, "USEBOUSS = ", use_bouss
            print *, '------------------------------------------------'
        end if
    end if
    !-------------------------------------------------------

    if (aux_finalized < 2 .and. actualstep) then
        ! topo arrays might have been updated by dtopo more recently than
        ! aux arrays were set unless at least 1 step taken on all levels
        aux(1,:,:) = NEEDS_TO_BE_SET ! new system checks this val before setting
        !Always recompute this array to avoid undefined values
        ! if (0==boussTransitionType) then
        !     !Depth-based transition, reset Bouss component of  aux array
        !     aux(dr_index,:,:) = NEEDS_TO_BE_SET
        ! end if
        call setaux(mbc,mx,my,xlower,ylower,dx,dy,maux,aux)
    endif

    !Set dispersive variables according to aux array
    ! call setaux(mbc,mx,my,xlower,ylower,dx,dy,maux,aux)
    

    ! if (boussEquations==-1 .or. boussEquations==-2) then
    if (boussEquations<0) then
        !EDC or BBBD (mbHSGN or HSGN)
        forall(i=1-mbc:mx+mbc, j=1-mbc:my+mbc)
            q(4:meqn,i,j) = aux(dr_index,i,j) * q(4:meqn,i,j)
        end forall
    end if

    ! if (boussEquations==-3) then
    if (.false.) then
        !DR (HiSGN)
        !Just loop inside grid because we need the (potentially ghost) 
        !neighboring values for sigma
        neighborhood_idx(:,1) = [-1,-1,-1,0 ,0, 1,1,1]
        neighborhood_idx(:,2) = [-1, 0, 1,-1,1,-1,0,1]
        do i=1,mx 
            do j=1,my
                q(4:meqn,i,j) = aux(dr_index,i,j) * q(4:meqn,i,j)
                !For this one, it does not suffice to zero out all the variables,
                !since the last variable (sigma for the moment), approximates 
                !sqrt(h)*1.159*(detadx + detady), which is used in the dispersive terms.
                !So, in cells where dispersive terms are turned off but are neighbors
                !to cells where they are on, we need to do something for sigma other than zeroing it out.
                !What I am doing now is setting it to the relaxation limit, computing limited
                !approximations of detadx and detady.
                if (aux(dr_index,i,j) < small_tol_bouss) Then
                    indicator_neighborhood = 0 !Tells us if any neighbor has dispersive terms on
                    ! Go around neighborhood of this cell
                    do kk=1,8
                        iidx = ii + neighborhood_idx(kk,1)
                        jidx = jj + neighborhood_idx(kk,2)
                        if (aux(dr_index,iidx,jidx) >= small_tol_bouss) then
                            indicator_neighborhood = 1
                        end if
                    end do
                    if (indicator_neighborhood == 1) then
                        !Compute limited detadx and detady
                        !This will not work at all when the coordinate system is not Cartesian
                        if (coordinate_system /= 1) then
                            stop "b4step2 HiSGN: Coordinate system not Cartesian. Cannot compute limited detadx and detady yet!"
                        end if
                        etaC = q(1,i,j) + aux(1,i,j)
                        etaL = q(1,i-1,j) + aux(1,i-1,j)
                        etaR = q(1,i+1,j) + aux(1,i+1,j)
                        etaB = q(1,i,j-1) + aux(1,i,j-1)
                        etaT = q(1,i,j+1) + aux(1,i,j+1)
                        detadx_limited = minmod( (etaR - etaC)/dx , (etaC - etaL)/dx )
                        detady_limited = minmod( (etaT - etaC)/dy , (etaC - etaB)/dy )
                        q(6,i,j) = q(1,i,j)*1.159d0 * sqrt(q(1,i,j)) * (detadx_limited + detady_limited)    
                    end if
                end if
            end do
        end do
    end if

    ! Set wind and pressure aux variables for this grid
    if (actualstep) then
        call set_storm_fields(maux,mbc,mx,my,xlower,ylower,dx,dy,t,aux)
    end if

    contains
    real(kind=8) function minmod(a,b)
        implicit none
        real(kind=8), intent(in) :: a,b
        if (a*b <= 0.d0) then
            minmod = 0.d0
        else
            minmod = sign(1.d0,a) * min(abs(a),abs(b))
        endif
    end function minmod

end subroutine b4step2
