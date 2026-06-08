!     =====================================================
subroutine qinit(meqn,mbc,mx,my,xlower,ylower,&
                  dx,dy,q,maux,aux)
!     =====================================================
!
!     # Set initial conditions for q.
!     # All the variables are 0 except for the real part of z0, which is a radial smoothed front (Busto et al. 2021).
!
   use map_module, only: mapc2p
   use cparam_mod, only: grav
   use cparam_mod, only: gamma => boussEDCgamma
   use cparam_mod, only: csq => boussEDCcsq
   use cparam_mod, only: boussEquations
   use solitary_module, only: solitary_wave_Bristeau
   implicit none 
   ! intrinsic :: sqrt, tanh
   double precision :: q(meqn, 1-mbc:mx+mbc, 1-mbc:my+mbc)
   double precision :: aux(maux, 1-mbc:mx+mbc, 1-mbc:my+mbc)

   integer :: meqn, mbc, mx, my, maux
   double precision :: xlower, ylower, dx, dy
   integer :: i, j

!     #Local variables
!Some variables for solitary wave
   double precision :: Hcap, A, d, x0, direction
   double precision :: h, u, v, p, w, etax
   double precision :: xcell, ycell
   double precision :: xp, yp
    

   !Setting parameters for solitary wave of Bristeau's system
   Hcap = 1.d0 !Still water depth
   A = 0.4d0 !Amplitude of solitary wave
   d = Hcap
   x0 = 0.d0
   direction = 1.d0

   q(:,:,:) = 0.d0
   q(1,:,:) = 1.d0

!     # set concentration profile
!     ---------------------------
   do i=1-mbc,mx+mbc
      xcell = xlower + (i-0.5d0)*dx
      do j=1-mbc,my+mbc
         ycell = ylower + (j-0.5d0)*dy 
         call mapc2p(xcell, ycell, xp, yp)
         ! q(1,i,j) = 1.d0 + 0.1*exp(-((xp-20.d0)**2.d0)-yp**2.d0)
         call solitary_wave_Bristeau(xp, Hcap, A, d,&
            gamma, grav, 0.d0, h, u, p, w, etax, x0, direction)  
      
         v = 0.d0
         q(1,i,j) = h
         q(2,i,j) = h*u
         q(3,i,j) = h*v
         if (csq > 1.e-3 .or. boussEquations == 0) then
            q(4,i,j) = h*w
            q(5,i,j) = h*p
         else
            q(4,i,j) = 0.d0
            q(5,i,j) = 0.d0
         end if
      enddo
   enddo
end subroutine qinit

      