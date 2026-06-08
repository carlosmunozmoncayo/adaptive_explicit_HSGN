! ==============================================================
subroutine evec(maxnx,num_eqn,num_ghost,mx,q,auxl,auxr,evl,evr)
! ==============================================================
!
!	Calculation of left and right eigenvectors
!
!   On input, q contains the cell average state vector
!             maxnx contains the number of physical points 
!               (without ghost cells)
!             num_ghost is the number of ghost cells
!             num_eqn the number of equations 
!
!   On output, evl(i) and evr(i) contain left/right eigenvectors
!   at interface i-1/2


    implicit none
    integer,          intent(in) :: maxnx, num_eqn, num_ghost, mx
    double precision, intent(in) :: q(num_eqn,maxnx+2*num_ghost)
    double precision, intent(in) :: auxl(num_eqn,maxnx+2*num_ghost)
    double precision, intent(in) :: auxr(num_eqn,maxnx+2*num_ghost)
    double precision, intent(out) :: evl(num_eqn,num_eqn,maxnx+2*num_ghost)
    double precision, intent(out) :: evr(num_eqn,num_eqn,maxnx+2*num_ghost)

!   local storage
!   ---------------
    ! common /cparam/ gamma
    integer :: i
    integer :: mx2

    double precision :: grav, c_sq, gama, eta_sq, dry_tolerance
    integer :: bouss_equations
    common /cparam/  grav, c_sq, gama, eta_sq, dry_tolerance,bouss_equations

    double precision :: hr, hl, ur, ul, wl, wr, pl, pr
    double precision :: sgr, sgl, pbl, pbr
    double precision :: csq_r, csq_l
    double precision :: hlsq, hrsq
    double precision :: hbar, ubar, wbar, sgbar, pbar, pbbar, csqbar
    double precision :: CEbar, lamb1stbar, lamblastbar

    mx2 = size(q,2)

    ! gamma1 = gamma - 1.d0

    do i=2,mx2
        hr = q(1,i)
        hl = q(1,i-1)
        ur = q(2,i)/hr
        ul = q(2,i-1)/hl
        if (bouss_equations == 1) then !EDC
            wr = q(3,i)/hr;    wl = q(3,i-1)/hl
            pr = q(4,i)/hr;    pl = q(4,i-1)/hl
            csq_r = auxr(2,i);  csq_l = auxl(2,i)
            !Roe averages
            hbar = 0.5d0*(hl+hr)
            hlsq = dsqrt(hl)
            hrsq = dsqrt(hr)
            ubar = (hlsq*ul+hrsq*ur)/(hlsq+hrsq)
            wbar = (hlsq*wl+hrsq*wr)/(hlsq+hrsq)
            pbar = (hlsq*pl+hrsq*pr)/(hlsq+hrsq)
            csqbar = 0.5d0*(csq_l+csq_r)
            CEbar=dsqrt(grav*hbar+csqbar+pbar)
            lamb1stbar=ubar-CEbar
            lamblastbar=ubar+CEbar

            ! Construct matrix of right eigenvectors
            evr(1,:,i)=[1.d0,           1.d0,       0.d0,      1.d0]
            evr(2,:,i)=[lamb1stbar,       ubar,       0.d0,      lamblastbar]
            evr(3,:,i)=[wbar,           0.d0,       1.d0,      wbar]
            evr(4,:,i)=[pbar+csqbar, -grav*hbar,    0.d0,      pbar+csqbar]

            ! Construct matrix of left eigenvectors
            evl(:,1,i)=[ubar*CEbar+grav*hbar,2*(csqbar+pbar),&
                -2*grav*hbar*wbar, -ubar*CEbar+grav*hbar]
            evl(:,2,i)=[-CEbar, 0.d0, 0.d0, CEbar]
            evl(:,3,i)=[0.d0, 0.d0, 2*CEbar**2, 0.d0]
            evl(:,4,i)=[1.d0,-2.d0,-2*wbar,1.d0]
            evl(:,:,i)=(0.5d0/CEbar**2)*evl(:,:,i)

            
        else if (bouss_equations==2) then !BBBD
            wr = q(3,i)/hr;    wl = q(3,i-1)/hl
            sgr = q(4,i)/hr;    sgl = q(4,i-1)/hl
            pr = q(5,i)/hr;    pl = q(5,i-1)/hl
            pbr = q(6,i)/hr;    pbl = q(6,i-1)/hl
            csq_r = auxr(2,i);  csq_l = auxl(2,i)

            !Roe averages
            hbar = 0.5*(hl+hr)
            hlsq = dsqrt(hl)
            hrsq = dsqrt(hr)
            ubar = (hlsq*ul+hrsq*ur)/(hlsq+hrsq)
            wbar = (hlsq*wl+hrsq*wr)/(hlsq+hrsq)
            sgbar = (hlsq*sgl+hrsq*sgr)/(hlsq+hrsq)
            pbar = (hlsq*pl+hrsq*pr)/(hlsq+hrsq)
            pbbar =(hlsq*pbl+hrsq*pbr)/(hlsq+hrsq)
            csqbar = 0.5d0*(csq_l+csq_r)
            CEbar=dsqrt(grav*hbar+csqbar+pbar)

            !Roe speeds
            lamb1stbar=ubar-CEbar
            lamblastbar=ubar+CEbar

            ! Construct matrix of right eigenvectors
            evr(1,:,i)= [1.d0,       1.d0,      0.d0, 0.d0, 0.d0, 1.d0]
            evr(2,:,i)=[lamb1stbar,    ubar,      0.d0, 0.d0, 0.d0, lamblastbar]
            evr(3,:,i)= [wbar,       0.d0,      0.d0, 1.d0, 0.d0, wbar]
            evr(4,:,i)= [sgbar,      0.d0,      1.d0, 0.d0, 0.d0,  sgbar]
            evr(5,:,i)= [pbar+csqbar,-grav*hbar,0.d0, 0.d0, 0.d0, pbar+csqbar]
            evr(6,:,i)= [pbbar,      0.d0,      0.d0, 0.d0, 1.d0,  pbbar]

            ! Construct matrix of left eigenvectors
            evl(:,:,i) = 0.d0
            evl(:,1,i) = [0.5d0*(grav*hbar+CEbar*ubar), csqbar+pbar, -grav*hbar*sgbar, &
                -grav*hbar*wbar, -grav*hbar*pbbar, 0.5d0*(grav*hbar-CEbar*ubar)]
            evl(:,1,i) = evl(:,1,i)/CEbar**2.d0
            evl(1,2,i) = -0.5d0/CEbar
            evl(6,2,i) = 0.5d0/CEbar
            evl(4,3,i) = 1.d0
            evl(3,4,i) = 1.d0
            evl(:,5,i)=[0.5d0,-1.d0,-sgbar,-wbar, -pbbar, 0.5d0]/CEbar**2.d0
            evl(5,6,i) = 1.d0
        endif
    enddo

    return
end subroutine evec
