
subroutine qinit(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux)
    
    use geoclaw_module, only: sea_level, g => grav
    use bouss_module, only: boussEquations
    use solitary_module, only: solitary_wave_Bristeau
    
    implicit none
    
    ! Subroutine arguments
    integer, intent(in) :: meqn,mbc,mx,my,maux
    real(kind=8), intent(in) :: xlower,ylower,dx,dy
    real(kind=8), intent(inout) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
    real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
    
    ! Locals
    integer :: i,j
    real(kind=8) :: xi, h0, ampl, L, d1, direction, eval_point
    real(kind=8) :: x, y
    real(kind=8) :: h, u, v, p, w !To receive the 1D SGN soliton

    L = 800.d0 !Lenght of the domain in the x and y direction (the domain is square)
    d1 = sqrt(2.d0)*L/4.d0 !Half the distance from the center to the corner of the domain

    h0 = 1.0d0 !Reference depth
    ampl = 0.32d0 !Soliton amplitude
    eval_point = 0.d0 
    direction = 1.0d0 !Direction of the soliton (1 for up-going, -1 for down-going)


    ! Set flat state based on sea_level
    q = 0.d0
    forall(i=1-mbc:mx+mbc, j=1-mbc:my+mbc)
        q(1,i,j) = max(0.d0, sea_level - aux(1,i,j))
    end forall


    do i=1-mbc,mx+mbc
      x = (xlower + (i-0.5d0)*dx)
      do j=1-mbc,my+mbc
        y = (ylower + (j-0.5d0)*dy)
        xi = (y+x)/sqrt(2.d0)  !Distance to the line y=-x
        if (xi .lt. 0) then !Set bottom left 
          if (abs(xi)<d1) then !Set up-going soliton 
            eval_point = xi+0.5d0*d1
            direction = 1.0d0
          else !Set down-going soliton
            eval_point = xi+1.5d0*d1
            direction = -1.0d0
          end if
          call solitary_wave_Bristeau(eval_point, h0, ampl, h0, 1.5d0, g, 0.d0, h, u, p, w, 0.d0, direction)
          !Set the state variables
          q(1,i,j) = h
          q(2,i,j) = u*q(1,i,j)*sqrt(2.d0)/2.d0 !Rotate velocity
          q(3,i,j) = u*q(1,i,j)*sqrt(2.d0)/2.d0 !Rotate velocity
          q(4,i,j) = q(1,i,j)*w
          q(5,i,j) = q(1,i,j)*p
        else !Set top right
          if (abs(xi)<d1) then !Set down-going soliton
            eval_point = xi-0.5d0*d1
            direction = -1.0d0
          else !Set up-going soliton
            eval_point = xi-1.5d0*d1
            direction = 1.0d0
          end if
          call solitary_wave_Bristeau(eval_point, h0, ampl, h0, 1.5d0, g, 0.d0, h, u, p, w, 0.d0, direction)
          !Set the state variables
          q(1,i,j) = h
          q(2,i,j) = u*q(1,i,j)*sqrt(2.d0)/2.d0 !Rotate velocity
          q(3,i,j) = u*q(1,i,j)*sqrt(2.d0)/2.d0 !Rotate velocity
          q(4,i,j) = q(1,i,j)*w
          q(5,i,j) = q(1,i,j)*p
        end if
        
      enddo
  enddo

    
end subroutine qinit
