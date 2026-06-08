subroutine qinit(meqn,mbc,mx,xlower,dx,q,maux,aux)

    ! Set initial conditions for the q array.

    use geoclaw_module, only: sea_level
    use geoclaw_module, only: grav
    use bouss_module, only : gamma => boussEDCgamma, boussEquations
    use grid_module, only: xcell
    use aux_module, only: read_data, bilinear_interp_1d, wellprepared

    ! uncomment if any of these needed...
    !use geoclaw_module, only: dry_tolerance, grav

    implicit none

    integer, intent(in) :: meqn,mbc,mx,maux
    real(kind=8), intent(in) :: xlower,dx
    real(kind=8), intent(in) :: aux(maux,1-mbc:mx+mbc)
    real(kind=8), intent(inout) :: q(meqn,1-mbc:mx+mbc)

    !locals
    integer :: i
    real(kind=8) :: h_at_rest, ampl, width, x0

    !Reading data from file
    real(kind=8), allocatable :: xdata(:), wdata(:), pdata(:)
    integer :: mxdata
    
    real(kind=8) :: w(mx), p(mx)


    h_at_rest = 100.d0
    ampl = 5.0d0  ! amplitude
    width = 100.    ! controls width of Gaussian
    x0 = 0.d0   ! initial location of Gaussian

    call read_data("w0_p0.txt",xdata, wdata, pdata, mxdata)

    do i=1,mx
        q(2:meqn,i) = 0.d0
        q(1,i) = h_at_rest + ampl * exp(-((xcell(i)-x0)/width)**2)
        w(i) = bilinear_interp_1d(xdata, wdata, mxdata, xcell(i))
        p(i) = bilinear_interp_1d(xdata, pdata, mxdata, xcell(i))
        if (wellprepared==1) then
            q(3,i) = q(1,i)*w(i)  ! momentum
            q(4,i) = q(1,i)*p(i)  ! auxiliary variable
        else
            q(3,i) = 0.d0
            q(4,i) = 0.d0
        end if
   enddo

end subroutine qinit
