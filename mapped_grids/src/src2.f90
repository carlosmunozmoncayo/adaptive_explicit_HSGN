subroutine src2(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux,t,dt)

!> Called to update q by solving source term equation 
!! $q_t = S(q)$ over time dt starting at time t.
!!
!! Since we can show that mass is conserved for this substep, the source
!! is linear and we can solve exactly by
!! q(dt)= exp(dt S) q(0)
    !Parameters for Boussinesq equations
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !For the moment, the only system that requires changes in the RS
    !is EDC HypRel (boussEquations = -1)
    use cparam_mod, only: boussEquations
    use cparam_mod, only: c_sq => boussEDCcsq 
    !Reference hyperbolic relaxation parameter (Just used for EDC)
    !Approximate SGN (3/2) or Sainte-Marie equations (2)
    use cparam_mod, only: gammaEDC => boussEDCgamma 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Mapped grid and geophysical parameters
    use cparam_mod, only:  g => grav, dry_tolerance, sea_level
    use cparam_mod, only: inclusion_idx, inclusion_type, srcintegrator
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Import integrators from bouss_module
    use bouss_module, only: one_exact_step_EDC
    use bouss_module, only: one_step_SDIRK22_EDC
    use bouss_module, only: one_step_BE_EDC
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    implicit none
    integer, intent(in) :: mbc,mx,my,meqn,maux
    real(kind=8), intent(in) :: xlower,ylower,dx,dy,t,dt
    real(kind=8), intent(in) ::  aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
    real(kind=8), intent(inout) ::  q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)

    !Locals
    real(kind=8) :: h
    ! Conserved quantities:
    !       1   h (water height)
    !       2   hu (x-momentum)
    !       3   hv (y-momentum)
    !       4   hw (vertical momentum)
    !       5   hp (non-hydrostatic pressure term)

    integer :: i,j

    if (boussEquations <0 .or. c_sq >1.e-3) then
        !Checking integrator type
        ! print *, 'srcintegrator =', srcintegrator    
        if (trim(srcintegrator)=='exact') then
            ! print *, 'Using exact integrator'
            ! loop over all grid points
            do i = 1, mx
                do j = 1, my
                        h = q(1,i,j)
                        call one_exact_step_EDC(h, q(4,i,j), q(5,i,j), c_sq, dt)
                end do
            end do

        else if (trim(srcintegrator)=='dirk22') then   
            ! print *, 'Using DIRK22 integrator'
            ! loop over all grid points
            do i = 1, mx
                do j = 1, my
                        h = q(1,i,j)
                        call one_step_SDIRK22_EDC(h, q(4,i,j), q(5,i,j), c_sq, dt)
                end do
            end do

        else if (trim(srcintegrator)=='be') then
            ! print *, 'Using Backward Euler integrator'
            ! loop over all grid points
            do i = 1, mx
                do j = 1, my
                        h = q(1,i,j)
                        call one_step_BE_EDC(h, q(4,i,j), q(5,i,j), c_sq, dt)
                end do
            end do

        else
            print *, 'Error: srcintegrator not recognized. Use exact, dirk22 or BE.'
            print *, 'You provided: ', trim(srcintegrator)
            stop
        end if
    end if


end subroutine src2