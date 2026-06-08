! =====================================================
subroutine rpt2(ixy,imp,maxm,meqn,mwaves,maux,mbc,mx,ql,qr,aux1,aux2,aux3,asdq,bmasdq,bpasdq)
! =====================================================

! # Riemann solver in the transverse direction for NLSH in 2D
! # Split asdq into down-going flux bmasdq and up-going flux bpasdq.
! # imp=1  means  asdq=amdq,    imp=2 means asdq=apdq
! # Not sure if this is well coded yet. We'll have to try.  
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
    use cparam_mod, only:  g => grav, drytol => dry_tolerance, sea_level
    use cparam_mod, only: inclusion_idx, inclusion_type
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    implicit none
    !Input
    integer, intent(in)  :: ixy,imp,maxm,meqn,mwaves,maux,mbc,mx
    double precision, dimension(meqn,1-mbc:maxm+mbc), intent(in) :: ql,qr
    double precision, dimension(maux,1-mbc:maxm+mbc), intent(in) :: aux1,aux2,aux3
    double precision, dimension(meqn,1-mbc:maxm+mbc), intent(in) :: asdq
    !Output
    double precision, dimension(meqn,1-mbc:maxm+mbc), intent(out) :: bmasdq,bpasdq


    !Locals
    integer :: i, mw, ix1, ixm1
    real(kind=8) :: delta(5,1), s(5)
    real(kind=8) :: waves(5,5)

    real(kind=8) :: h,u,v,w,p


    !For mapped grid
    integer :: ixtran, iytran, ilenrat, i1
    double precision :: st, ct

    if (ixy.eq.1) then
        ixtran = 4
        iytran = 5
        ilenrat = 6
    else
        ixtran = 1
        iytran = 2
        ilenrat = 3
    endif

    bmasdq(:,:)  = 0.d0
    bpasdq(:,:)  = 0.d0

    if (imp.eq.1) then
        ! asdq = amdq, moving to left
        ix1 = 2-mbc
        ixm1 = mx+mbc
    else
        !asdq = apdq, moving to right
        ix1 = 1-mbc
        ixm1 = mx+mbc
    endif

    do i=ix1,ixm1
        i1 = i + imp - 2
        !Get rotation angle from aux array
        !up-going:
        ct = aux3(ixtran,i1)
        st = aux3(iytran,i1)
        !Rotate forward asdq
        delta(:,1) = asdq(:,i)
        delta(2,1) = ct*asdq(2,i)+st*asdq(3,i)
        delta(3,1) = -st*asdq(2,i)+ct*asdq(3,i)
        !Get primitive variables
        h = ql(1,i1)
        u = (ct*ql(2,i1) + st*ql(3,i1)) / h
        v = (-st*ql(2,i1) + ct*ql(3,i1)) / h
        w = ql(4,i1) / h
        p = ql(5,i1)
        !Get waves from Riemann solver (already rotated back inside)
        call mapped_transverse_HypRel_EDC(ct,st,u,v,h,w,p,g,c_sq,delta,waves,s)
        !Split fluctuations
        do mw=1,mwaves
            if (s(mw) > 1.e-16) then
                bpasdq(:,i) = bpasdq(:,i) + aux3(ilenrat,i1)*s(mw)*waves(:,mw)
            endif
        end do

        !down-going:
        ct = aux2(ixtran,i1)
        st = aux2(iytran,i1)
        !Rotate forward asdq
        delta(:,1) = asdq(:,i)
        delta(2,1) = ct*asdq(2,i)+st*asdq(3,i)
        delta(3,1) = -st*asdq(2,i)+ct*asdq(3,i)
        !Get primitive variables
        u = (ct*ql(2,i1) + st*ql(3,i1)) / h
        v = (-st*ql(2,i1) + ct*ql(3,i1)) / h
        !Get waves from Riemann solver (already rotated back inside)
        call mapped_transverse_HypRel_EDC(ct,st,h,u,v,w,p,g,c_sq,delta,waves,s)
        !Get coefficients for the 1d problem in x direction
        !Split fluctuations
        do mw=1,mwaves
            if (s(mw) < -1.e-16) then
                bmasdq(:,i) = bmasdq(:,i) + aux2(ilenrat,i1)*s(mw)*waves(:,mw)
            endif
        end do      
    end do

    contains
        subroutine mapped_transverse_HypRel_EDC(ct,st,h,u,v,w,p,grav,c_sq,localasdq,waves, s)
            implicit none
            !Input
            double precision, intent(in) :: ct,st !Rotation cos and sin
            double precision, intent(in) :: u, v, h, w, p, c_sq, localasdq(5,1)
            double precision, intent(in) :: grav
            !u should be normal velocity component
            !v should be tangential velocity component

            !Output
            double precision, intent(out) ::  waves(5,5),s(5) 

            !Local
            double precision :: CE, lamb1, lamb5, R(5,5),  R_inv(5,5), alpha(5,1)
            integer :: k

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

            alpha = matmul(R_inv, localasdq)

            !Rotate back waves
            do k=1,5
                waves(1,k) = alpha(k,1)*R(1,k)
                waves(2,k) = alpha(k,1)*(ct*R(2,k)-st*R(3,k))
                waves(3,k) = alpha(k,1)*(st*R(2,k)+ct*R(3,k))
                waves(4,k) = alpha(k,1)*R(4,k)
                waves(5,k) = alpha(k,1)*R(5,k)
            end do
        end subroutine mapped_transverse_HypRel_EDC
end subroutine rpt2
