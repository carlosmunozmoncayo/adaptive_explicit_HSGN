subroutine src1d(meqn,mbc,mx1d,q1d,maux,aux1d,t,dt)

    ! This routine should be a simplified version of src2 
    ! which applies source terms for a 1-d slice of data along the
    ! edge of a grid.  This is called only from qad where the conservative
    ! fix-up is applied and is used to apply source terms over partial
    ! time steps to the coarse grid cell values used in solving Riemann 
    ! problems at the interface between coarse and fine grids.
 
    ! If the source terms depend only on q, it should be easy to 
    ! adapt src2 to create this routine, just loop over 1:mx1d.
    ! If the source terms are more complicated, it may not be easy.
 
    ! The code may work fine without applying source terms in this
    ! context, so using this dummy routine might be successful even when
    ! source terms are present. 
 
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
    integer, intent(in) :: meqn,mbc,mx1d,maux
    real(kind=8), intent(in) :: t, dt
    real(kind=8), intent(in) ::  aux1d(maux,mx1d)
    real(kind=8), intent(inout) ::  q1d(meqn,mx1d)


     !Locals
    real(kind=8) :: h
    ! Conserved quantities:
    !       1   h (water height)
    !       2   hu (x-momentum)
    !       3   hv (y-momentum)
    !       4   hw (vertical momentum)
    !       5   hp (non-hydrostatic pressure term)

    integer :: i

    if (boussEquations <0 .or. c_sq >1.e-3) then
        !Checking integrator type
        ! print *, 'srcintegrator =', srcintegrator    
        if (trim(srcintegrator)=='exact') then
            ! print *, 'Using exact integrator'
            ! loop over all grid points
            do i=1,mx1d
                ! if (aux1d(dr_index,i)>small_tol) then
                h = q1d(1,i)
                call one_exact_step_EDC(h, q1d(4,i), q1d(5,i), c_sq, dt)
                ! end if
            enddo
        else if (trim(srcintegrator)=='dirk22') then   
            ! print *, 'Using DIRK22 integrator'
            ! loop over all grid points
            do i = 1, mx1d
                h = q1d(1,i)
                call one_step_SDIRK22_EDC(h, q1d(4,i), q1d(5,i), c_sq, dt)
            end do
        else if (trim(srcintegrator)=='be') then
            ! print *, 'Using Backward Euler integrator'
            ! loop over all grid points
            do i = 1, mx1d
                h = q1d(1,i)
                call one_step_BE_EDC(h, q1d(4,i), q1d(5,i), c_sq, dt)
            end do
        else
            print *, 'Error: srcintegrator not recognized. Use exact, dirk22 or BE.'
            print *, 'You provided: ', trim(srcintegrator)
            stop
        end if
    endif

end subroutine src1d
