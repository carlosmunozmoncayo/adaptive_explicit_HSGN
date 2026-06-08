! ============================================================================
subroutine tfluct1(maxnx,meqn,mwaves,maux,mbc,mx,ql,qr,auxl,auxr,adq)
    ! ============================================================================
    !   "Internal" Riemann solver for EDC and BBBD BTEs in 1D.
    !   The riemann problem is solved by assuming a discontinuity at the
    !   center of the i'th cell.
    !
    !   On input, ql contains the state vector at the left edge of each cell
    !             qr contains the state vector at the right edge of each cell
    !             auxl contains the auxiliary vector at the left edge of each cell
    !             auxr contains the state vector at the right edge of each cell
    !             maxnx is the number of physical points (without ghost cells)
    !             mbc is the number of ghost cells
    !             meqn is the number of equations
    !             ixyz is the dimension index
    !             mx is the size of the patch for the dimension corresponding
    !               to the value of ixyz
    !
    !   On output, adq contains the decomposition of the flux difference
    !              f(qr(i)) - f(ql(i)).

    
    
        implicit none
        integer,          intent(in)  :: maxnx, mx, meqn, mwaves, maux, mbc
        double precision, intent(in)  :: auxl(maux,1-mbc:maxnx+mbc)
        double precision, intent(in)  :: auxr(maux,1-mbc:maxnx+mbc)
        double precision, intent(in)  ::   ql(meqn,1-mbc:maxnx+mbc)
        double precision, intent(in)  ::   qr(meqn,1-mbc:maxnx+mbc)
    
        double precision, intent(out) :: adq(meqn,1-mbc:maxnx+mbc)
    
        integer :: i
    
    !   local storage
    !   ---------------
        double precision :: hr, hl, ur, ul, wl, wr, pl, pr
        double precision :: csq_r, csq_l, br, bl
        double precision :: sgr, sgl, pbl, pbr
        double precision :: hlsq, hrsq
        double precision :: hbar, ubar, wbar, sgbar, pbar, pbbar, csqbar, delb
        double precision :: grav, c_sq, gama, eta_sq, dry_tolerance
        integer :: bouss_equations
        common /cparam/  grav, c_sq, gama, eta_sq, dry_tolerance,bouss_equations
   
        do i = 1,mx
            hr = qr(1,i)
            hl = ql(1,i)
            ur = qr(2,i)/hr
            ul = ql(2,i)/hl
            if (bouss_equations == 1) then !EDC
                wl = qr(3,i)/hr;    wr = ql(3,i)/hl
                pr = qr(4,i)/hr;    pl = ql(4,i)/hl
                csq_r = auxr(2,i);  csq_l = auxl(2,i)
                br = auxr(1,i);     bl = auxl(1,i)
                delb = br-bl
                hbar = 0.5*(hl+hr)
                hlsq = dsqrt(hl);   hrsq = dsqrt(hr)
                ubar = (hlsq*ul+hrsq*ur)/(hlsq+hrsq)
                wbar = (hlsq*wl+hrsq*wr)/(hlsq+hrsq)
                pbar = (hlsq*pl+hrsq*pr)/(hlsq+hrsq)
                csqbar = (hlsq*csq_l+hrsq*csq_l)/(hlsq+hrsq)
                pbar = (hlsq*pl+hrsq*pr)/(hlsq+hrsq)

                adq(1,i) = hr*ur-hl*ul
                adq(2,i) = hr*ur**2+hr*pr+0.5d0*grav*hr**2-&
                            hl*ul**2-hl*pl-0.5d0*grav*hl**2+&
                            (grav*hbar+gama*pbar)*(br-bl)
                adq(3,i) = hr*ur*wr-hl*ul*wl
                adq(4,i) = hr*ur*(csq_r+pr)-hl*ul*(csq_l+pl) &
                            -csqbar*ubar*(hr-hl)

            else if (bouss_equations == 2) then !BBBD
                wl = qr(3,i)/hr;    wr = ql(3,i)/hl
                sgr = qr(4,i)/hr;   sgl = ql(4,i)/hl
                pr = qr(5,i)/hr;    pl = ql(5,i)/hl
                pbr = qr(6,i)/hr;   pbl = ql(6,i)/hl
                br = auxr(1,i);     bl = auxl(1,i)
                csq_r = auxr(2,i);  csq_l = auxl(2,i)
                hbar = 0.5*(hl+hr)
                hlsq = dsqrt(hl)
                hrsq = dsqrt(hr)
                ubar = (hlsq*ul+hrsq*ur)/(hlsq+hrsq)
                wbar = (hlsq*wl+hrsq*wr)/(hlsq+hrsq)
                sgbar = (hlsq*sgl+hrsq*sgr)/(hlsq+hrsq)
                pbar = (hlsq*pl+hrsq*pr)/(hlsq+hrsq)
                pbbar =(hlsq*pbl+hrsq*pbr)/(hlsq+hrsq)
                csqbar = (hlsq*csq_l+hrsq*csq_l)/(hlsq+hrsq)
                delb = br-bl
                adq(1,i) = hr*ur-hl*ul
                adq(2,i) = hr*ur**2.d0 +hr*pr+0.5d0*grav*hr**2.d0 &
                        - hl*ul**2.d0 -hl*pl-0.5d0*grav*hl**2.d0 &
                        + (grav*hbar+pbbar)*delb
                adq(3,i) = hr*ur*wr-hl*ul*wl
                adq(4,i) = hr*ur*sgr-hl*ul*sgl
                adq(5,i) = hr*ur*(csq_r+pr)-hl*ul*(csq_l+pl) &
                        -csqbar*ubar*(hr-hl)
                adq(6,i) = hr*ur*pbr-hl*ul*pbl &
                        -6.d0*csqbar*ubar*delb
            end if
        enddo
        return
    end subroutine tfluct1

    