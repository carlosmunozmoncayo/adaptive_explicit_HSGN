module aux_module
    implicit none
    integer :: wellprepared
    contains
subroutine read_data(fname_tail, x,w,p, mx)
    !Reads data from a text file.
    character(len=*), intent(in) :: fname_tail
    real(8), allocatable, intent(out) :: x(:), w(:), p(:)
    integer, intent(out) :: mx

    integer :: iounit, i, j, ios
    character(len=512)  :: env_path
    character(len=1024) :: fullpath
    integer :: length

    ! Get environment variable LOCALGEOCLAW (set by geoclaw/local_CLAW_env.sh)
    call get_environment_variable("LOCALGEOCLAW", env_path, length, status=ios)
    if (ios /= 0) then
      print *, "Error: environment variable LOCALGEOCLAW not set (source geoclaw/local_CLAW_env.sh)"
      stop
    end if
    env_path = env_path(:length)

    ! Build full path
    fullpath = trim(env_path)//"/examples/1d_classic/well_prepared_test/"//trim(fname_tail)
    print *, "Reading data from file: ", trim(fullpath)

    ! Open text file
    open(newunit=iounit, file=trim(fullpath), status='old', action='read', iostat=ios)
    if (ios /= 0) then
      print *, "Error opening file: ", trim(fullpath)
      stop
    end if

    ! Read dimensions
    read(iounit, *) mx
    allocate(x(mx), w(mx), p(mx))

    ! Read x coordinates
    read(iounit, *) x
    ! Read w values    
    read(iounit, *) w
    ! Read p values    
    read(iounit, *) p

    close(iounit)
end subroutine read_data


function bilinear_interp_1d(x, u_data, mx, xi) result(val)
    ! Performs bilinear interpolation on a 2D grid.
    ! Inputs:
    !   x      : 1D arrays of grid coordinates in x direction
    !   u_data     : 1D array of data values at grid points
    !   mx         : number of grid points in x direction
    !   xi         : coordinate where interpolation is desired
    ! Output:
    !   val        : interpolated value at xi
    !How to use:
    !   val = bilinear_interp_1d(x, u_data, mx, xi)
    real(8), intent(in) :: x(:), u_data(:)
    integer, intent(in) :: mx
    real(8), intent(in) :: xi
    real(8) :: val
    integer :: i
    real(8) :: x1,x2, q1,q2,tx

    ! Find index i such that x(i)<=xi<=x(i+1)
    i = max(1, min(mx-1, count(x <= xi)))

    x1 = x(i);   x2 = x(i+1)
    q1 = u_data(i);     q2 = u_data(i+1)

    tx = (xi-x1)/(x2-x1)

    val = (1-tx)*q1 + tx*q2
  end function bilinear_interp_1d

end module aux_module
