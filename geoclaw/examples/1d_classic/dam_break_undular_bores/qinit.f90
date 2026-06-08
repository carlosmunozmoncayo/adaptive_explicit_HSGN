subroutine qinit(meqn,mbc,mx,xlower,dx,q,maux,aux)

    ! Set initial conditions for the q array.

    use geoclaw_module, only: sea_level
    use geoclaw_module, only: grav
    use bouss_module, only : gamma => boussEDCgamma, boussEquations
    use grid_module, only: xcell

    ! uncomment if any of these needed...
    !use geoclaw_module, only: dry_tolerance, grav

    implicit none

    integer, intent(in) :: meqn,mbc,mx,maux
    real(kind=8), intent(in) :: xlower,dx
    real(kind=8), intent(in) :: aux(maux,1-mbc:mx+mbc)
    real(kind=8), intent(inout) :: q(meqn,1-mbc:mx+mbc)

    !locals
    integer :: i
    real(kind=8) :: H1, H2, Froude, effective_Froude
    real(kind=8) :: alpha

    Froude = 1.1702d0 !flag_froude_number PLEASE DO NOT DELETE THIS COMMENT
    H1 = 1.d0
    !Fr = sqrt(H2/(2*H1**2)*(H2+H1))

    if (abs(Froude-1.1702d0) < 1.d-10) then
        H1 = 0.2477
        H2 = 0.3685
        ! effective_Froude = sqrt(H2/(2*H1**2)*(H2+H1))
        ! print *, "Using H1 and H2 from digitized data for Froude 1.1702"
        ! print *, "Effective Froude number: ", effective_Froude
    else if (abs(Froude-1.104d0) < 1.d-10) then
        H1 = 0.2510
        H2 = 0.3230
        ! effective_Froude = sqrt(H2/(2*H1**2)*(H2+H1))
        ! print *, "Using H1 and H2 from digitized data for Froude 1.104"
        ! print *, "Effective Froude number: ", effective_Froude
    end if

    ! H2 = sqrt(1+8*Froude**2)-1.d0
    ! H2 = 0.5d0*H2 !The only positive solution
    ! H2 = 0.5d0*H1*(sqrt(1+8*Froude**2)-1.d0) !Compute H2 from H1 and Froude number
    Froude = sqrt(0.5d0*(H2/H1)*(1+H2/H1)) !Compute Froude number from H1 and H2 to check consistency
    print *, "H1, H2, Froude: ", H1, H2, Froude


    do i=1,mx
        q(2:meqn,i) = 0.d0
        ! if (xcell(i) < 0.0d0) then
        !     q(1,i) = H2
        ! else
        !     q(1,i) = H1
        ! end if
        ! alpha = transition_function(xcell(i),-0.2d0,0.2d0,3) !Smooth out dam break initial condition
        if (xcell(i) <0.d0) then
            alpha = 0.d0
        else
            alpha = 1.d0
        end if
        q(1,i) = H1+(1-alpha)*(H2-H1)
   enddo

   contains
    function transition_function(s,s0,s1,optional_transition_type_fun,optional_shape) result(fs)
        !Returns f(s), where f(s0)=0 and f(s1)=1
        implicit none
        double precision, intent(in) :: s, s0, s1
        integer, optional, intent(in) :: optional_transition_type_fun
        double precision, optional, intent(in) :: optional_shape
        !Output
        double precision :: fs
        !Locals
        double precision :: shape = 1.d0 !Different effects depending on transition_type
        integer :: transition_type_fun = 1
        double precision :: sn !Normalized s (0<=sn<=1)
        !Edge cases
        if (s < s0) then
            fs = 0.d0
            return
        else if (s >= s1) then
            fs = 1.d0
            return
        end if
        !Normalize and get linear transition
        sn = (s-s0)/(s1-s0)
        if (present(optional_shape)) shape = optional_shape
        if (present(optional_transition_type_fun)) transition_type_fun = optional_transition_type_fun
        !Transition function
        select case(transition_type_fun)
            case(1) !Linear (can be made a Heaviside if s0==s1)
                fs = sn
            case(2) !Hyperbolic tangent
                !Larger shape -> sharper transition, better match at 0 and 1
                !Smaller shape -> smoother transition, does not match 0 and 1
                fs = 0.5*(tanh(10.d0*shape*(sn-0.5))+1)
            case(3) !Engsig-Karup (2007)
                !Matches 0 and 1, smooth in the middle, sharper at the ends
                fs = -2.d0*sn**3+3.d0*sn**2
            case(4) !Mayer et al. (1998)
                !Very smooth at 0, does not match 1
                !Smaller shape -> smoother transition at 0, better match at 1
                !and faster decrease close to 1
                fs = (0.01d0*shape*sn)**3+(1.d0-0.01d0*shape)*sn**6
            case(5) !Modified Mayer et al. (1998)
                !Very smooth at 1, does not match 0
                !Smaller shape -> smoother transition at 1, better match at 0
                !and faster decrease close to 0
                fs = 1.d0-(0.01d0*shape*(1.d0-sn))**3+(1.d0-0.01d0*shape)*(1.d0-sn)**6
            case default
                fs = sn !Linear
        end select
    end function transition_function

end subroutine qinit
