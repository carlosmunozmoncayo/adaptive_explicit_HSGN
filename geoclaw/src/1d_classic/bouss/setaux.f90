subroutine setaux(mbc,mx,xlower,dx,maux,aux)

    ! Called at start of computation before calling qinit, and
    ! when AMR is used, also called every time a new grid patch is created.
    ! Use to set auxiliary arrays aux(1:maux, 1-mbc:mx+mbc, 1-mbc:my+mbc).
    ! Note that ghost cell values may need to be set if the aux arrays
    ! are used by the Riemann solver(s).
    !
    ! This version sets aux(1,:) to b(x) for shallow flow
    ! and aux(2,:) to the ratio of cell width to dxc (capacity function).
    !
    ! For the hyperbolic relaxation systems (boussEquations < 0: EDC = mbHSGN, BBBD = HSGN)
    ! it also sets, as in the 2D code (geoclaw/src/2d/hyp_bouss/setaux.f90):
    !   aux(bouss_csq_index,:)      = c^2 = alpha^2 * g * H * f   (relaxation parameter)
    !   aux(boussDecayRate_index,:) = f, with f = 1 where the dispersive system is solved and
    !                                 f = 0 where the SWEs are solved,
    ! where alpha = boussEDCalpha (bouss_EDC_alpha in setrun.py), H = sea_level - b is the depth
    ! at rest, minimum over the cell and its two neighbours (clipped at 0), and f = 1 iff
    ! H > boussMinDepth (bouss_min_depth in setrun.py; the same criterion as useBouss).
    ! This requires num_aux >= 4 in setrun.py.

    !use geoclaw_module, only: dry_tolerance !uncomment if needed
    use geoclaw_module, only: grav, sea_level
    use geoclaw_module, only: earth_radius, coordinate_system, DEG2RAD
    use grid_module, only: xp_edge,mx_edge
    use topo_module, only: zcell
    use bouss_module, only: boussEquations, boussMinDepth
    use bouss_module, only: alpha_relax => boussEDCalpha
    use bouss_module, only: csq_index => bouss_csq_index
    use bouss_module, only: dr_index => boussDecayRate_index

    implicit none
    integer, intent(in) :: mbc,mx,maux
    real(kind=8), intent(in) :: xlower,dx
    real(kind=8), intent(inout) ::  aux(maux,1-mbc:mx+mbc)

    !locals
    integer :: i,i0,i1,j
    real(kind=8) :: xcell,a
    real(kind=8) :: h_rest(1-mbc:mx+mbc), h_at_rest

    if (mx+1 .ne. mx_edge) then
        write(6,*) 'mx_edge from celledges.data must agree with mx+1'
        write(6,*) 'mx_edge = ',mx_edge
        write(6,*) 'mx = ',mx
        stop
        endif

    ! store topo values in aux(1,:):
    ! computation of cell averages is done in read_topo_file or topo_update

    !write(6,*) '+++ Before resetting aux, aux(1,200) = ',aux(1,200)

    do i=1,mx
        aux(1,i) = zcell(i)
        !write(6,*) '+++ i,zcell: ',i,zcell(i)
        ! capacity function for nonuniform grid:
        aux(2,i) = (xp_edge(i+1) - xp_edge(i))/dx
        if (coordinate_system == 2) then
            ! additional factor for latitude coords on sphere:
            ! dx in meters = (dx in degrees) * (pi/180) * earth_radius
            aux(2,i) = aux(2,i)* DEG2RAD * earth_radius
        endif
        if (aux(2,i) <= 0.d0) then
            write(6,*) '+++ i,xp_edge(i),xp_edge(i+1): ',i,xp_edge(i),xp_edge(i+1)
            endif
    enddo


    aux(:,0) = aux(:,1)
    aux(:,-1) = aux(:,1)
    aux(:,mx+1) = aux(:,mx)
    aux(:,mx+2) = aux(:,mx)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Hyperbolic relaxation parameter and SWE/Boussinesq indicator for EDC and BBBD
    if (boussEquations < 0) then
        if (maux < max(csq_index, dr_index)) then
            write(6,*) '*** setaux: the hyperbolic relaxation systems (bouss_equations < 0) need'
            write(6,*) '*** num_aux >= ', max(csq_index, dr_index), ' in setrun.py, got ', maux
            stop
        endif
        do i=1-mbc,mx+mbc
            h_rest(i) = max(sea_level - zcell(i), 0.d0)
        enddo
        do i=1-mbc,mx+mbc
            ! depth at rest: minimum over the cell and its neighbours (inside the array)
            h_at_rest = minval(h_rest(max(i-1,1-mbc):min(i+1,mx+mbc)))
            if (h_at_rest > boussMinDepth) then
                aux(dr_index,i) = 1.d0
            else
                aux(dr_index,i) = 0.d0
            endif
            aux(csq_index,i) = (alpha_relax**2)*grav*h_at_rest * aux(dr_index,i)
        enddo
    endif

end subroutine setaux
