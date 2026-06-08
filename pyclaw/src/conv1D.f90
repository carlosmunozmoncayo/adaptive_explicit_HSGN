subroutine conv1d(qin,qout,filter_window,T,T0,mx,mw,meqn,dx)
    implicit none
    integer, intent(in) :: mx,mw,meqn
    double precision, intent(in) :: qin(meqn,mx), dx
    double precision, intent(in) :: filter_window(mw)
    double precision, intent(in) :: T0
    double precision, intent(in) :: T(mx)
    double precision, intent(out) :: qout(meqn,mx)

    !For all the points i where T(i)>T0, replace q(meqn,i) by
    !the convolution of q(meqn,i) with the filter_window
    !Locals
    integer :: i,j

    do i=1,mx
        !Copy input to output
        qout(:,i) = qin(:,i)
        if (T(i)>T0) then
            !Convolution
            do j=1,meqn
                qout(j,i) = sum(filter_window*qin(j,i-mw: i+mw))*dx
            end do
        end if
    end do
end subroutine conv1d