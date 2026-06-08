subroutine qinit(meqn,mbc,mx,xlower,dx,q,maux,aux)

    ! Set initial conditions for the q array.

    use geoclaw_module, only: sea_level
    use geoclaw_module, only: grav
    use bouss_module, only : gamma => boussEDCgamma, boussEquations
    use grid_module, only: xcell
    use solitary_module, only: solitary_wave_Bristeau

    ! uncomment if any of these needed...
    !use geoclaw_module, only: dry_tolerance, grav

    implicit none

    integer, intent(in) :: meqn,mbc,mx,maux
    real(kind=8), intent(in) :: xlower,dx
    real(kind=8), intent(in) :: aux(maux,1-mbc:mx+mbc)
    real(kind=8), intent(inout) :: q(meqn,1-mbc:mx+mbc)

    !locals
    integer :: i
    real(kind=8) :: xi, h0, ampl, L, d1, direction, eval_point
    real(kind=8) :: x, y
    real(kind=8) :: h, u, v, p, w !To receive the 1D SGN soliton

    L = 800.d0 !Lenght of the domain in the x and y direction (the domain is square)
    d1 = sqrt(2.d0)*L/4.d0 !Half the distance from the center to the corner of the domain

    h0 = 1.0d0 !Reference depth
    ampl = 0.32d0 !Soliton amplitude
    eval_point = 0.d0 
    direction = 1.0d0 !Direction of the soliton (1 for up-going, -1 for down-going)


    do i=1,mx
        xi = abs(xcell(i)) !Distance to the line y=-x
         if (xcell(i)<0.d0) then !Set bottom left 
          if (xi<d1) then !Set up-going soliton 
            eval_point = xi-0.5d0*d1
            direction = 1.0d0
          else !Set down-going soliton
            eval_point = xi-1.5d0*d1
            direction = -1.0d0
          end if
          call solitary_wave_Bristeau(eval_point, h0, ampl, h0, 1.5d0, grav, 0.d0, h, u, p, w, 0.d0, direction)
          !Set the state variables
          q(1,i) = h
          q(2,i) = u*q(1,i)
          if (boussEquations <0) then !HSGN
            q(3,i) = -q(1,i)*w
            q(4,i) = q(1,i)*p
          end if
        else !Set top right
          if (xi<d1) then !Set down-going soliton
            eval_point = xi-0.5d0*d1
            direction = -1.0d0
          else !Set up-going soliton
            eval_point = xi-1.5d0*d1
            direction = 1.0d0
          end if
          call solitary_wave_Bristeau(eval_point, h0, ampl, h0, 1.5d0, grav, 0.d0, h, u, p, w, 0.d0, direction)
          !Set the state variables
          q(1,i) = h
          q(2,i) = u*q(1,i)
          if (boussEquations <0) then !HSGN
            q(3,i) = q(1,i)*w
            q(4,i) = q(1,i)*p
          end if
        end if
        if (xcell(i)<0.d0) then !Set bottom left 
          if (xi<d1) then !Set up-going soliton 
            eval_point = xi-0.5d0*d1
            direction = 1.0d0
          else !Set down-going soliton
            eval_point = xi-1.5d0*d1
            direction = -1.0d0
          end if
          call solitary_wave_Bristeau(eval_point, h0, ampl, h0, 1.5d0, grav, 0.d0, h, u, p, w, 0.d0, direction)
          !Set the state variables
          q(1,i) = h
          q(2,i) = u*q(1,i)
          if (boussEquations <0) then !HSGN
            q(3,i) = -q(1,i)*w
            q(4,i) = q(1,i)*p
          end if
        else !Set top right
          if (xi<d1) then !Set down-going soliton
            eval_point = xi-0.5d0*d1
            direction = -1.0d0
          else !Set up-going soliton
            eval_point = xi-1.5d0*d1
            direction = 1.0d0
          end if
          call solitary_wave_Bristeau(eval_point, h0, ampl, h0, 1.5d0, grav, 0.d0, h, u, p, w, 0.d0, direction)
          !Set the state variables
          q(1,i) = h
          q(2,i) = u*q(1,i)
          if (boussEquations <0) then !HSGN
            q(3,i) = q(1,i)*w
            q(4,i) = q(1,i)*p
          end if
        end if
   enddo
end subroutine qinit
