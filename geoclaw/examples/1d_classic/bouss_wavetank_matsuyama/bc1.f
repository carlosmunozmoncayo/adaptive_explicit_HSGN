
c
c
c     =================================================================
      subroutine bc1(meqn,mbc,mx,xlower,dx,q,maux,aux,t,dt,mthbc)
c     =================================================================
c
c     # Standard boundary condition choices for claw2
c
c     # Modified for 1d GeoClaw to extend topo aux(1,:) to ghost cells.
c
c     # At each boundary  k = 1 (left),  2 (right):
c     #   mthbc(k) =  0  for user-supplied BC's (must be inserted!)
c     #            =  1  for zero-order extrapolation
c     #            =  2  for periodic boundary coniditions
c     #            =  3  for solid walls, assuming this can be implemented
c     #                  by reflecting the data about the boundary and then
c     #                  negating the 2'nd component of q.
c     ------------------------------------------------
c
c     # Extend the data from the computational region
c     #      i = 1, 2, ..., mx2
c     # to the virtual cells outside the region, with
c     #      i = 1-ibc  and   i = mx+ibc   for ibc=1,...,mbc
c
      use geoclaw_module, only: grav
      use bouss_module, only: c_sq => boussEDCcsq !Reference hyperbolic relaxation parameter
      use bouss_module, only: gamma => boussEDCgamma !Approximate SGN (3/2) or Sainte-Marie equations (2)

      implicit double precision (a-h,o-z)
      dimension q(meqn,1-mbc:mx+mbc)
      dimension aux(maux,1-mbc:mx+mbc)

      dimension mthbc(2)

      pi = 4.d0*atan(1.d0)
      h0 = 4.d0
      
      ampl_eta = 0.03
      ampl_u = ampl_eta * sqrt(grav/h0)
      !write(6,*) '+++ ampl_eta, ampl_u: ',ampl_eta, ampl_u
      !ampl = 0.047  #0.15d0
      tperiod = 20.d0

c
c
c-------------------------------------------------------
c     # left boundary:
c-------------------------------------------------------
      do ibc=1,mbc
         aux(1,1-ibc)=aux(1,1)
      enddo

      go to (100,110,120,130) mthbc(1)+1
c
  100 continue
c     # wave maker
      if (t < tperiod) then
          s = ampl_u * sin(2.d0*pi*t/tperiod)
      else
          s = 0.d0
      endif
      do ibc=1,mbc
          aux(1,1-ibc) = aux(1,ibc)
          do m=1,meqn
              q(m,1-ibc) = q(m,ibc)
          enddo
          u = q(2,ibc)/q(1,ibc)
          q(2:meqn,1-ibc) = 0.d0
          q(2,1-ibc) = (2.d0*s - u) * q(1,1-ibc)
      enddo

      !Setting the last 2 components of q
      ! do ibc=1,mbc !Using FD approximations of Bristeau's definition of p and w
      !    h0 = q(1,1-ibc)
      !    u0 = q(2,1-ibc)/h0
      !    h1 = q(1,ibc)
      !    u1 = q(2,ibc)/h1
      !    w = -(h0+h1)*(u1-u0)/(8.d0*dx)
      !    !One sided derivatives
      !    h0u0 = q(2,1-ibc)
      !    h1u1 = q(2,2-ibc)
      !    h2u2 = q(2,3-ibc)
      !    dxh0u0 = ((-3.d0/2.d0)*h0u0 + 2.d0*h1u1 - (1.d0/2.d0)*h2u2)/dx
      !    h1u1 = q(2,ibc)
      !    h0u0 = q(2,ibc-1)
      !    hm1um1 = q(2,ibc-2)
      !    dxh1u1 = ((1.d0/2.d0)*hm1um1 - 2.d0*h0u0 + (3.d0/2.d0)*h1u1)/dx
      !    p = -gamma*w*(h1*u1-h0*u0)/(2.d0*dx)
      !    p = p - gamma*((h1+h0)/2.d0)*(u1-u0)*(dxh0u0+dxh1u1)/(8.d0*dx)

      !    q(3,1-ibc) = (h0+h1)*w - q(3,ibc)
      !    q(4,1-ibc) = (h0+h1)*p - q(4,ibc)
      ! enddo

      ! do ibc=1,mbc !Using multiples of the right going eigenvector
      !    h0 = q(1,1-ibc)
      !    hu0 = q(2,1-ibc)
      !    u0 = hu0/h0
      !    w0 = q(3,1-ibc)/h0
      !    p0 = q(4,1-ibc)/h0
      !    eigenval4 = u0+sqrt(grav*h0+c_sq+p0)
      !    q(1,1-ibc) = hu0/eigenval4
      !    q(3,1-ibc) = hu0*w0/eigenval4
      !    q(4,1-ibc) = hu0*(p0+c_sq)/eigenval4
      ! enddo

      
      go to 199
c
  110 continue
c     # zero-order extrapolation:
      do ibc=1,mbc
         do m=1,meqn
            q(m,1-ibc) = q(m,1)
         end do
      end do
      go to 199

  120 continue
c     # periodic:
      do ibc=1,mbc
         do m=1,meqn
            q(m,1-ibc) = q(m,mx+1-ibc)
         end do
      end do
      go to 199

  130 continue
c     # solid wall (assumes 2'nd component is velocity or momentum in x):
      do ibc=1,mbc
         do m=1,meqn
            q(m,1-ibc) = q(m,ibc)
         end do
c        # negate the normal velocity:
         q(2,1-ibc) = -q(2,ibc)
      end do
      go to 199

  199 continue

c
c-------------------------------------------------------
c     # right boundary:
c-------------------------------------------------------
      do ibc=1,mbc
         aux(1,mx+ibc)=aux(1,mx)
      enddo

      go to (200,210,220,230) mthbc(2)+1
c
  200 continue
c     # user-specified boundary conditions go here in place of error output
      write(6,*) '*** ERROR *** mthbc(2)=0 and no BCs specified in bc2'
      stop
      go to 299

  210 continue
c     # zero-order extrapolation:
      do ibc=1,mbc
         do m=1,meqn
            q(m,mx+ibc) = q(m,mx)
         end do
      end do
      go to 299

  220 continue
c     # periodic:
      do ibc=1,mbc
         do m=1,meqn
            q(m,mx+ibc) = q(m,ibc)
         end do
      end do
      go to 299

  230 continue
c     # solid wall (assumes 2'nd component is velocity or momentum in x):
      do ibc=1,mbc
         do m=1,meqn
            q(m,mx+ibc) = q(m,mx+1-ibc)
         end do
c        # negate the normal velocity:
         q(2,mx+ibc) = -q(2,mx+1-ibc)
      end do
      go to 299

  299 continue
c
      return
      end
