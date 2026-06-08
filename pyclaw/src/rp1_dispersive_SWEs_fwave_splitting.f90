! =========================================================
subroutine rp1(maxmx,meqn,mwaves,maux,mbc,mx,ql,qr,auxl,auxr,fwave,s,amdq,apdq)
! =========================================================

! solve Riemann problem for the 1D shallow water equations using David George's
! augmented Riemann solver or Escalante's system using a Roe-type Riemann solver.
! waves: 4
! equations: 4
! Conserved quantities:
!       1 depth (h)
!       2 momentum (hu)
!       3 depth times vertical velocity (hw)
!       4 depth times non-hydrostatic pressure (hp)
!The first component of the auxiliary array contains the bathymetry
!The second component of the auxiliary array contains a flag that
!determines wether SWEs or Escalante's hyperbolic relaxation of Saint Marie's system is solved
!on the corresponding cell
    implicit none !(type, external)
    !external :: sgesv  
    
    integer, intent(in) :: maxmx, meqn, mwaves, mbc, mx, maux
    double precision, dimension(meqn,1-mbc:maxmx+mbc), intent(in) :: ql, qr
    double precision, dimension(maux,1-mbc:maxmx+mbc), intent(in) :: auxl, auxr
    double precision, dimension(meqn, 1-mbc:maxmx+mbc), intent(out) :: amdq, apdq
    double precision, dimension(mwaves, 1-mbc:maxmx+mbc), intent(out) :: s
    double precision, dimension(meqn, mwaves, 1-mbc:maxmx+mbc), intent(out) :: fwave

    double precision :: u_l, u_r, h_l, h_r
    double precision :: w_l, w_r, p_l, p_r, b_l, b_r, csq_l, csq_r 
    double precision :: hbar
    double precision :: sg_l, sg_r, pb_l, pb_r
    integer :: i, k
    double precision :: grav, c_sq, gama, eta_sq

    !EDC and BBBD RSs output
    double precision :: fwEDC(5,5), swEDC(5)
    double precision :: fwBBBD(7,7), swBBBD(7)
    integer :: bouss_equations
    integer :: skip_fourth_idx(6)
    integer :: skip_fifth_idx(6)

    !Additional variables for Geoclaw's solver
    integer maxiter
    real(kind=8) fw(2,3) !2 equations, 3 waves for 1D solver
    real(kind=8) sw(3)  
    real(kind=8) hu_l,hu_r,se1,se2
    real(kind=8) hvl,hvr,vl,vr
    real(kind=8) phil, phir
    double precision :: dry_tolerance, chat, uhat
    double precision :: alphaLDR
    !PROCEDURE(riemann_aug_JCP), POINTER :: PTR_TO_SUB

    common /cparam/  grav, c_sq, gama, eta_sq, dry_tolerance,bouss_equations, alphaLDR
    skip_fourth_idx = (/ 1, 2, 3, 5, 6, 7 /)
    skip_fifth_idx = (/ 1, 2, 3, 4, 6, 7 /)
    maxiter=1

    do i=2-mbc,mx+mbc
        h_l = qr(1,i-1)
        h_r = ql(1,i)
        u_l = qr(2,i-1) / qr(1,i-1)
        u_r = ql(2,i  ) / ql(1,i  )
        if (bouss_equations==1) then !EDC
            w_l = qr(3,i-1) / qr(1,i-1)
            w_r = ql(3,i  ) / ql(1,i  )
            p_l = qr(4,i-1) / qr(1,i-1)
            p_r = ql(4,i  ) / ql(1,i  )
        else if (bouss_equations==2) then !BBBD
            w_l = qr(3,i-1) / qr(1,i-1)
            w_r = ql(3,i  ) / ql(1,i  )
            sg_l = qr(4,i-1) / qr(1,i-1)
            sg_r = ql(4,i  ) / ql(1,i  )
            p_l = qr(5,i-1) / qr(1,i-1)
            p_r = ql(5,i  ) / ql(1,i  )
            pb_l = qr(6,i-1) / qr(1,i-1)
            pb_r = ql(6,i  ) / ql(1,i  )
        else if (bouss_equations==3) then  !DR
            w_l = qr(3,i-1) / qr(1,i-1)
            w_r = ql(3,i  ) / ql(1,i  )
            p_l = qr(4,i-1) / qr(1,i-1)
            p_r = ql(4,i  ) / ql(1,i  )
            sg_l = qr(5,i-1) / qr(1,i-1)
            sg_r = ql(5,i  ) / ql(1,i  )
        end if
        b_l = auxr(1, i - 1)
        b_r = auxl(1,i)
        csq_l = auxr(2, i - 1)
        csq_r = auxl(2,i)

        c_sq = 0.5d0*(csq_l+csq_r)
 
         if ((max(csq_l,csq_r) < 1.e-10) .or. (bouss_equations==0) &! if ((max(csq_l,csq_r) < 1.e-10) .or. (bouss_equations==0) &
            .or. (h_l < dry_tolerance) .or. (h_r < dry_tolerance)) then
         !Solve SWEs
            ! # Left states
            h_l = qr(1, i - 1)
            if (h_l > dry_tolerance) then
                u_l = qr(2, i - 1) / h_l
                phil = qr(1, i - 1) * u_l**2 + 0.5d0 * grav * qr(1, i - 1)**2
            else
                u_l = 0.d0
                phil = 0.d0
            end if
            b_l = auxr(1, i - 1)
            
            ! # Right states
            h_r = ql(1, i)
            if (h_r > dry_tolerance) then
                u_r = ql(2, i)/h_r
                phir = ql(1, i) * u_r**2 + 0.5d0 * grav * ql(1, i)**2
            else
                u_r = 0.d0
                phir = 0.d0
            end if
            b_r = auxl(1, i)

            ! # Roe average states
            hbar = 0.5 * (h_r + h_l)
            uhat = (sqrt(h_l) * u_l + sqrt(h_r) * u_r) / (sqrt(h_l) + sqrt(h_r))
            chat = sqrt(grav * hbar)
            
            ! # Wave speeds
            sE1 = min(u_l - sqrt(grav * h_l), uhat - chat)
            sE2 = max(u_r + sqrt(grav * h_r), uhat + chat)
            !sE2 = 0.5d0*(s(1,i)+s(3,i))


            hu_L=h_l*u_l
            hu_R=h_r*u_r

            hvL=0.d0
            hvR=0.d0
            vl=0.d0
            vr=0.d0
            
            sw=0.d0
            fw=0.d0 
            
            call riemann_aug_JCP(maxiter,2,3,h_L,h_R,hu_L,hu_R, &
            hvL,hvR,b_l,b_r,u_L,u_R,vL,vR,phiL,phiR,sE1,sE2,dry_tolerance,grav,sw,fw)
            
            !Speeds
            s(1,i) = sw(1)
            s(2:mwaves-1,i) = sw(2)
            s(mwaves,i) = sw(3)
            
            !fwaves
            do k=1,mwaves
                fwave(:,k,i)=0.d0
            end do  
            fwave(1:2,1,i) = fw(:,1)
            fwave(1:2,2,i) = fw(:,2)
            fwave(1:2,3:mwaves-1,i) = 0.d0
            fwave(1:2,mwaves,i) = fw(:,3)
            
         else if (bouss_equations==1) then
         !Solve EDC
            ! print *, "SOLVING EDC RS!!!"
            call riemann_normal_HypRel_EDC(2, 3, grav, gama, c_sq, &
                               h_l, u_l, 0.d0, w_l, p_l, b_l, csq_l,&
                               h_r, u_r, 0.d0, w_r, p_r, b_r, csq_r, &
                               fwEDC, swEDC)
            do k=1,mwaves
               s(k,i) = swEDC(skip_fourth_idx(k))
               fwave(1:2,k,i) = fwEDC(1:2,skip_fourth_idx(k)) 
               !fwEDC(3,:) has no use in 1D
               fwave(3:4,k,i) = fwEDC(4:5,skip_fourth_idx(k))
            end do
            ! call riemman_1D_HypRel_EDC(gama ,grav, h_l, u_l, w_l, p_l, b_l, csq_l, &
            ! h_r, u_r, w_r, p_r, b_r, csq_r, fwEDC1D, swEDC1D)
            ! do k=1,mwaves
            !    s(k,i) = swEDC1D(k)
            !    fwave(:,k,i) = fwEDC1D(:,k) 
            ! end do
            
         else if (bouss_equations==2) then
         !Solve BBBD
            ! print *, "SOLVING BBBD RS!!!"
            call riemann_normal_HypRel_BBBD(2, 3, grav, c_sq, &
                               h_l, u_l, 0.d0, w_l, sg_l, p_l, pb_l, b_l, csq_l,&
                               h_r, u_r, 0.d0, w_r, sg_r, p_r, pb_r, b_r, csq_r, &
                               fwBBBD, swBBBD)
            do k=1,mwaves
               s(k,i) = swBBBD(skip_fifth_idx(k))
               fwave(1:2,k,i) = fwBBBD(1:2,skip_fifth_idx(k)) 
               !fwBBBD(3,:) has no use in 1D
               fwave(3:6,k,i) = fwBBBD(4:7,skip_fifth_idx(k))
            end do

         else if (bouss_equations==3) then
         !Solve DR
            call riemann_1D_HypRel_DR(grav, alphaLDR, c_sq, &
                              h_l, u_l, w_l, p_l, sg_l, b_l, csq_l, &
                              h_r, u_r, w_r, p_r, sg_r, b_r, csq_r, &
                              fwEDC, swEDC) !The same number of entries
            do k=1,mwaves
               s(k,i) = swEDC(k)
               fwave(:,k,i) = fwEDC(:,k) 
            end do
         end if

         ! Computing left and right going fluctuations
         amdq(:,i)=0.d0
         apdq(:,i)=0.d0
         do k=1, mwaves
             if (s(k, i) < 0.d0) then
                 amdq(:, i) = amdq(:, i) + fwave(:, k, i)
             else if (s(k, i) > 0.d0) then
                 apdq(:, i) = apdq(:, i) + fwave(:, k, i)
             else
                  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  ! Debugging
                  ! if (max(h_l+b_l, h_r+b_r)>1.e-5) then
                  !    print *, "Warning: Riemann solver returned zero speed at i=", i, " k=", k
                  !    print *, "s(k,i)=", s(k, i), "h_l=", h_l, "h_r=", h_r, "u_l=", u_l, "u_r=", u_r
                  !    print *, "Location (i)=", i
                  ! endif
                  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 amdq(:, i) = amdq(:, i) + 0.5d0 * fwave(:, k, i)
                 apdq(:, i) = apdq(:, i) + 0.5d0 * fwave(:, k, i)
             end if
         enddo
         if (h_l<dry_tolerance.AND. h_r<dry_tolerance) then
            amdq(:,i)=0.d0
            apdq(:,i)=0.d0
         end if
    end do 

end subroutine rp1

!=============================================================================
!Subroutines for SWEs
!=============================================================================
!-----------------------------------------------------------------------
subroutine riemann_aug_JCP(maxiter,num_eqn,num_waves,hL,hR,huL,huR, &
         hvL,hvR,bL,bR,uL,uR,vL,vR,phiL,phiR,sE1,sE2,drytol,g,sw,fw)

         ! solve shallow water equations given single left and right states
         ! This solver is described in J. Comput. Phys. (6): 3089-3113, March 2008
         ! Augmented Riemann Solvers for the Shallow Equations with Steady States and Inundation

         ! To use the original solver call with maxiter=1.

         ! This solver allows iteration when maxiter > 1. The iteration seems to help with
         ! instabilities that arise (with any solver) as flow becomes transcritical over variable topo
         ! due to loss of hyperbolicity.



         implicit none

         !input
         integer num_eqn,num_waves,maxiter
         real(kind=8) fw(num_eqn,num_waves)
         real(kind=8) sw(num_waves)
         real(kind=8) hL,hR,huL,huR,bL,bR,uL,uR,phiL,phiR,sE1,sE2
         real(kind=8) hvL,hvR,vL,vR
         real(kind=8) drytol,g


         !local
         integer m,mw,k,iter
         real(kind=8) A(3,3)
         real(kind=8) r(3,3)
         real(kind=8) lambda(3)
         real(kind=8) del(3)
         real(kind=8) beta(3)

         real(kind=8) delh,delhu,delphi,delb,delnorm
         real(kind=8) rare1st,rare2st,sdelta,raremin,raremax
         real(kind=8) criticaltol,convergencetol,raretol
         real(kind=8) s1s2bar,s1s2tilde,hbar,hLstar,hRstar,hustar
         real(kind=8) huRstar,huLstar,uRstar,uLstar,hstarHLL
         real(kind=8) deldelh,deldelphi
         real(kind=8) s1m,s2m,hm
         real(kind=8) det1,det2,det3,determinant

         logical rare1,rare2,rarecorrector,rarecorrectortest,sonic

         !determine del vectors
         delh = hR-hL
         delhu = huR-huL
         delphi = phiR-phiL
         delb = bR-bL
         delnorm = delh**2 + delphi**2

         call riemanntype(hL,hR,uL,uR,hm,s1m,s2m,rare1,rare2,1,drytol,g)


         lambda(1)= min(sE1,s2m) !Modified Einfeldt speed
         lambda(3)= max(sE2,s1m) !Modified Eindfeldt speed
         sE1=lambda(1)
         sE2=lambda(3)
         lambda(2) = 0.d0  ! ### Fix to avoid uninitialized value in loop on mw -- Correct?? ###

         
         hstarHLL = max((huL-huR+sE2*hR-sE1*hL)/(sE2-sE1),0.d0) ! middle state in an HLL solve

   !     !determine the middle entropy corrector wave------------------------
         rarecorrectortest=.false.
         rarecorrector=.false.
         if (rarecorrectortest) then
            sdelta=lambda(3)-lambda(1)
            raremin = 0.5d0
            raremax = 0.9d0
            if (rare1.and.sE1*s1m.lt.0.d0) raremin=0.2d0
            if (rare2.and.sE2*s2m.lt.0.d0) raremin=0.2d0
            if (rare1.or.rare2) then
               !see which rarefaction is larger
               rare1st=3.d0*(sqrt(g*hL)-sqrt(g*hm))
               rare2st=3.d0*(sqrt(g*hR)-sqrt(g*hm))
               if (max(rare1st,rare2st).gt.raremin*sdelta.and. & 
                  max(rare1st,rare2st).lt.raremax*sdelta) then
                     rarecorrector=.true.
                  if (rare1st.gt.rare2st) then
                     lambda(2)=s1m
                  elseif (rare2st.gt.rare1st) then
                     lambda(2)=s2m
                  else
                     lambda(2)=0.5d0*(s1m+s2m)
                  endif
               endif
            endif
            if (hstarHLL.lt.min(hL,hR)/5.d0) rarecorrector=.false.
         endif

   !     ## Is this correct 2-wave when rarecorrector == .true. ??
         do mw=1,num_waves
            r(1,mw)=1.d0
            r(2,mw)=lambda(mw)
            r(3,mw)=(lambda(mw))**2
         enddo
         if (.not.rarecorrector) then
            lambda(2) = 0.5d0*(lambda(1)+lambda(3))
   !         lambda(2) = max(min(0.5d0*(s1m+s2m),sE2),sE1)
            r(1,2)=0.d0
            r(2,2)=0.d0
            r(3,2)=1.d0
         endif
   !     !---------------------------------------------------
   !     !determine the steady state wave -------------------
         criticaltol = 1.d-6
         deldelh = -delb
         deldelphi = -g*0.5d0*(hR+hL)*delb

   !     !determine a few quanitites needed for steady state wave if iterated
         hLstar=hL
         hRstar=hR
         uLstar=uL
         uRstar=uR
         huLstar=uLstar*hLstar
         huRstar=uRstar*hRstar

         !iterate to better determine the steady state wave
         convergencetol=1.d-6
         do iter=1,maxiter
            !determine steady state wave (this will be subtracted from the delta vectors)
            if (min(hLstar,hRstar).lt.drytol.and.rarecorrector) then
               rarecorrector=.false.
               hLstar=hL
               hRstar=hR
               uLstar=uL
               uRstar=uR
               huLstar=uLstar*hLstar
               huRstar=uRstar*hRstar
               lambda(2) = 0.5d0*(lambda(1)+lambda(3))
   !           lambda(2) = max(min(0.5d0*(s1m+s2m),sE2),sE1)
               r(1,2)=0.d0
               r(2,2)=0.d0
               r(3,2)=1.d0
            endif

            hbar =  max(0.5d0*(hLstar+hRstar),0.d0)
            s1s2bar = 0.25d0*(uLstar+uRstar)**2 - g*hbar
            s1s2tilde= max(0.d0,uLstar*uRstar) - g*hbar

   !        !find if sonic problem
            sonic=.false.
            if (abs(s1s2bar).le.criticaltol) sonic=.true.
            if (s1s2bar*s1s2tilde.le.criticaltol) sonic=.true.
            if (s1s2bar*sE1*sE2.le.criticaltol) sonic = .true.
            if (min(abs(sE1),abs(sE2)).lt.criticaltol) sonic=.true.
            if (sE1.lt.0.d0.and.s1m.gt.0.d0) sonic = .true.
            if (sE2.gt.0.d0.and.s2m.lt.0.d0) sonic = .true.
            if ((uL+sqrt(g*hL))*(uR+sqrt(g*hR)).lt.0.d0) sonic=.true.
            if ((uL-sqrt(g*hL))*(uR-sqrt(g*hR)).lt.0.d0) sonic=.true.

   !        !find jump in h, deldelh
            if (sonic) then
               deldelh =  -delb
            else
               deldelh = delb*g*hbar/s1s2bar
            endif
   !        !find bounds in case of critical state resonance, or negative states
            if (sE1.lt.-criticaltol.and.sE2.gt.criticaltol) then
               deldelh = min(deldelh,hstarHLL*(sE2-sE1)/sE2)
               deldelh = max(deldelh,hstarHLL*(sE2-sE1)/sE1)
            elseif (sE1.ge.criticaltol) then
               deldelh = min(deldelh,hstarHLL*(sE2-sE1)/sE1)
               deldelh = max(deldelh,-hL)
            elseif (sE2.le.-criticaltol) then
               deldelh = min(deldelh,hR)
               deldelh = max(deldelh,hstarHLL*(sE2-sE1)/sE2)
            endif

   !        !find jump in phi, deldelphi
            if (sonic) then
               deldelphi = -g*hbar*delb
            else
               deldelphi = -delb*g*hbar*s1s2tilde/s1s2bar
            endif
   !        !find bounds in case of critical state resonance, or negative states
            deldelphi=min(deldelphi,g*max(-hLstar*delb,-hRstar*delb))
            deldelphi=max(deldelphi,g*min(-hLstar*delb,-hRstar*delb))

            del(1)=delh-deldelh
            del(2)=delhu
            del(3)=delphi-deldelphi

   !        !Determine determinant of eigenvector matrix========
            det1=r(1,1)*(r(2,2)*r(3,3)-r(2,3)*r(3,2))
            det2=r(1,2)*(r(2,1)*r(3,3)-r(2,3)*r(3,1))
            det3=r(1,3)*(r(2,1)*r(3,2)-r(2,2)*r(3,1))
            determinant=det1-det2+det3

   !        !solve for beta(k) using Cramers Rule=================
            do k=1,3
               do mw=1,3
                     A(1,mw)=r(1,mw)
                     A(2,mw)=r(2,mw)
                     A(3,mw)=r(3,mw)
               enddo
               A(1,k)=del(1)
               A(2,k)=del(2)
               A(3,k)=del(3)
               det1=A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))
               det2=A(1,2)*(A(2,1)*A(3,3)-A(2,3)*A(3,1))
               det3=A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1))
               beta(k)=(det1-det2+det3)/determinant
            enddo

            !exit if things aren't changing
            if (abs(del(1)**2+del(3)**2-delnorm).lt.convergencetol) exit
            delnorm = del(1)**2+del(3)**2
            !find new states qLstar and qRstar on either side of interface
            hLstar=hL
            hRstar=hR
            uLstar=uL
            uRstar=uR
            huLstar=uLstar*hLstar
            huRstar=uRstar*hRstar
            do mw=1,num_waves
               if (lambda(mw).lt.0.d0) then
                  hLstar= hLstar + beta(mw)*r(1,mw)
                  huLstar= huLstar + beta(mw)*r(2,mw)
               endif
            enddo
            do mw=num_waves,1,-1
               if (lambda(mw).gt.0.d0) then
                  hRstar= hRstar - beta(mw)*r(1,mw)
                  huRstar= huRstar - beta(mw)*r(2,mw)
               endif
            enddo

            if (hLstar.gt.drytol) then
               uLstar=huLstar/hLstar
            else
               hLstar=max(hLstar,0.d0)
               uLstar=0.d0
            endif
            if (hRstar.gt.drytol) then
               uRstar=huRstar/hRstar
            else
               hRstar=max(hRstar,0.d0)
               uRstar=0.d0
            endif

         enddo ! end iteration on Riemann problem

         do mw=1,num_waves
            sw(mw)=lambda(mw)
            fw(1,mw)=beta(mw)*r(2,mw)
            fw(2,mw)=beta(mw)*r(3,mw)  !This is equation (21) of the paper
            !fw(3,mw)=beta(mw)*r(2,mw)
         enddo
      

         return

   end subroutine riemann_aug_JCP!-------------------------------------------------

   subroutine riemanntype(hL,hR,uL,uR,hm,s1m,s2m,rare1,rare2, &
                     maxiter,drytol,g)

         !determine the Riemann structure (wave-type in each family)


         implicit none

         !input
         real(kind=8) hL,hR,uL,uR,drytol,g
         integer maxiter

         !output
         real(kind=8) s1m,s2m
         logical rare1,rare2

         !local
         real(kind=8) hm,u1m,u2m,um,delu
         real(kind=8) h_max,h_min,h0,F_max,F_min,dfdh,F0,slope,gL,gR
         integer iter



   !     !Test for Riemann structure

         h_min=min(hR,hL)
         h_max=max(hR,hL)
         delu=uR-uL

         if (h_min.le.drytol) then
            hm=0.d0
            um=0.d0
            s1m=uR+uL-2.d0*sqrt(g*hR)+2.d0*sqrt(g*hL)
            s2m=uR+uL-2.d0*sqrt(g*hR)+2.d0*sqrt(g*hL)
            if (hL.le.0.d0) then
               rare2=.true.
               rare1=.false.
            else
               rare1=.true.
               rare2=.false.
            endif

         else
            F_min= delu+2.d0*(sqrt(g*h_min)-sqrt(g*h_max))
            F_max= delu + (h_max-h_min)*(sqrt(.5d0*g*(h_max+h_min)/(h_max*h_min)))

            if (F_min.gt.0.d0) then !2-rarefactions

               hm=(1.d0/(16.d0*g))* &
                        max(0.d0,-delu+2.d0*(sqrt(g*hL)+sqrt(g*hR)))**2
               um=sign(1.d0,hm)*(uL+2.d0*(sqrt(g*hL)-sqrt(g*hm)))

               s1m=uL+2.d0*sqrt(g*hL)-3.d0*sqrt(g*hm)
               s2m=uR-2.d0*sqrt(g*hR)+3.d0*sqrt(g*hm)

               rare1=.true.
               rare2=.true.

            elseif (F_max.le.0.d0) then !2 shocks

   !           !root finding using a Newton iteration on sqrt(h)===
               h0=h_max
               do iter=1,maxiter
                  gL=sqrt(.5d0*g*(1/h0 + 1/hL))
                  gR=sqrt(.5d0*g*(1/h0 + 1/hR))
                  F0=delu+(h0-hL)*gL + (h0-hR)*gR
                  dfdh=gL-g*(h0-hL)/(4.d0*(h0**2)*gL)+gR-g*(h0-hR)/(4.d0*(h0**2)*gR)
                  slope=2.d0*sqrt(h0)*dfdh
                  h0=(sqrt(h0)-F0/slope)**2
               enddo
                  hm=h0
                  u1m=uL-(hm-hL)*sqrt((.5d0*g)*(1/hm + 1/hL))
                  u2m=uR+(hm-hR)*sqrt((.5d0*g)*(1/hm + 1/hR))
                  um=.5d0*(u1m+u2m)

                  s1m=u1m-sqrt(g*hm)
                  s2m=u2m+sqrt(g*hm)
                  rare1=.false.
                  rare2=.false.

            else !one shock one rarefaction
               h0=h_min

               do iter=1,maxiter
                  F0=delu + 2.d0*(sqrt(g*h0)-sqrt(g*h_max)) &
                        + (h0-h_min)*sqrt(.5d0*g*(1/h0+1/h_min))
                  slope=(F_max-F0)/(h_max-h_min)
                  h0=h0-F0/slope
               enddo

               hm=h0
               if (hL.gt.hR) then
                  um=uL+2.d0*sqrt(g*hL)-2.d0*sqrt(g*hm)
                  s1m=uL+2.d0*sqrt(g*hL)-3.d0*sqrt(g*hm)
                  s2m=uL+2.d0*sqrt(g*hL)-sqrt(g*hm)
                  rare1=.true.
                  rare2=.false.
               else
                  s2m=uR-2.d0*sqrt(g*hR)+3.d0*sqrt(g*hm)
                  s1m=uR-2.d0*sqrt(g*hR)+sqrt(g*hm)
                  um=uR-2.d0*sqrt(g*hR)+2.d0*sqrt(g*hm)
                  rare2=.true.
                  rare1=.false.
               endif
            endif
         endif

         return

   end  subroutine riemanntype!----------------------------------------------------------------

!=============================================================================
!Subroutines for EDC system
!=============================================================================
subroutine riemann_normal_HypRel_EDC(mu, nv, grav, gamma, c_sq, &
                               hl, ul, vl, wl, pl, bl, csq_l,&
                               hr, ur, vr, wr, pr, br, csq_r, &
                               fw, sw)

   !Normal 2D Riemann solver for Escalante, Dumbser, and Castro's
   !hyperbolic relaxation of Sainte-Marie's dispersive water wave system
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !When used in 2D:
   !  Input: mu and nv according to the normal direction (analogous to Geoclaw's solver for SWEs)
   !  Output: 5 f-waves with 5 components each (fw) and 5 wave speeds (sw)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !When used in 1D :
   !  Input: mu=2, nv =3, vl=vr=0
   !  Output: The 3rd component of all the f-waves must be discarded
   !          The 3rd f-wave and the 3rd wave speed must be discarded
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !${CLAW}/riemann/src/geoclaw_riemann_utils.f is probably a better home for this routine
   !I'll leave it here for now
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !Using hydrostatic reconstruction for the bathymetry source term

   implicit none

   !Input
   integer :: mu, nv
   double precision :: grav, gamma, c_sq
   double precision :: hl, ul, vl, wl, pl, bl, csq_l
   double precision :: hr, ur, vr, wr, pr, br, csq_r

   !Output
   double precision :: fw(5,5), sw(5)

   !Local
   double precision :: hlsq, hrsq, etaL, etaR
   double precision :: hbar, ubar, vbar, wbar, pbar, csqbar, CEbar
   double precision :: lamb1bar, lamb5bar
   double precision :: lamb1tilde, lamb5tilde
   double precision :: R(5,5), R_inv(5,5)
   double precision :: alpha(5,1), total_fluct(5,1)
   double precision :: bath_correction_term

   double precision :: bath_sharpness_sensor !For stability term

   integer :: mw

   etaL = hL + bL
   etaR = hR + bR
   !Note that now, etaL = etaR
   !Roe averages
   hbar = 0.5d0*(hL+hR)
   hlsq = dsqrt(hL)
   hrsq = dsqrt(hR)
   ubar = (hlsq*ul+hrsq*ur)/(hlsq+hrsq)
   vbar = 0.d0!(hlsq*vl+hrsq*vr)/(hlsq+hrsq)
   wbar = (hlsq*wl+hrsq*wr)/(hlsq+hrsq)
   pbar = (hlsq*pl+hrsq*pr)/(hlsq+hrsq)

   !Still need to verify if there is difference between taking
   !the average, min, or max
   csqbar = 0.5d0*(csq_l+csq_r)

   CEbar=dsqrt(grav*hbar+csqbar+pbar)

   !Roe speeds
   lamb1tilde=ubar-CEbar
   lamb5tilde=ubar+CEbar

   !Einfeldt speeds
   lamb1bar=min(lamb1tilde, ul-dsqrt(grav*hl+csq_l+pl))
   lamb5bar=max(lamb5tilde, ur+dsqrt(grav*hr+csq_r+pr))
   
   !Matrix of right eigenvectors
   R(1,:)=[1.d0,     1.d0,       0.d0,    0.d0,    1.d0]
   R(mu,:)=[lamb1bar, ubar,      0.d0,    0.d0,    lamb5bar]
   R(nv,:)=[vbar,    0.d0,       0.d0,    1.d0,    vbar]
   R(4,:)=[wbar,     0.d0,       1.d0,    0.d0,    wbar]
   R(5,:)=[pbar+csqbar, -grav*hbar,0.d0,  0.d0,   pbar+csqbar]

   !Computing total fluctuation to split in f-waves
   total_fluct(1,1)=hr*ur-hl*ul
   total_fluct(mu,1)=hr*ur**2+hr*pr+0.5d0*grav*hr**2-&
                  hl*ul**2-hl*pl-0.5d0*grav*hl**2+&
                  (grav*hbar+gamma*pbar)*(br-bl)
   total_fluct(nv,1)=ur*hr*vr-ul*hl*vl
   total_fluct(4,1)=ur*hr*wr-ul*hl*wl

   total_fluct(5,1)=ur*hr*(pr+csq_r)-ul*hl*(pl+csq_l)!
   bath_correction_term = -csqbar*ubar*((etaR-etaL)+bath_sharpness_sensor(bl,br)*(br-bl)) 
   !The last term in the last expression is a consistency/stability term,
   !it is suppressed when there are large relative bathymetry jumps
   total_fluct(5,1) = total_fluct(5,1) + bath_correction_term 

   !Solving system explicitly for the 5x5 case
   R_inv(:,1)=[ubar*CEbar+grav*hbar,2*(csqbar+pbar),&
      -2*grav*hbar*wbar, 0.d0, -ubar*CEbar+grav*hbar]
   R_inv(:,mu)=[-CEbar, 0.d0, 0.d0, 0.d0, CEbar]
   R_inv(:,nv)=[0.d0, 0.d0, 0.d0, 2.d0*CEbar**2, 0.d0] 
   R_inv(:,4)=[0.d0, 0.d0, 2*CEbar**2, 0.d0, 0.d0]
   R_inv(:,5)=[1.d0,-2.d0,-2*wbar,0.d0,1.d0]
   R_inv=(0.5d0/CEbar**2)*R_inv

   alpha = matmul(R_inv,total_fluct)
   
   ! Computing f-waves
   do mw = 1, 5 !Replace by mwaves if necessary
      fw(:, mw) = alpha(mw,1) * R(:, mw)
   END do
   
   !Computing wave speeds
   sw(1) = min(ul-sqrt(grav*hl+csqbar+pl),lamb1bar)
   sw(2) = ubar
   sw(3) = ubar
   sw(4) = ubar
   sw(5) = max(ur+sqrt(grav*hr+csqbar+pr),lamb5bar)

end subroutine riemann_normal_HypRel_EDC


function bath_sharpness_sensor(bL,bR)
   !Computes a sensor to detect sharp bathymetry changes
   !relative to the water depth. Returns a value between 0 and 1
   !Cubic polynomial transition

   implicit none

   !Input
   double precision,intent(in) :: bL, bR
   !Output
   double precision :: bath_sharpness_sensor

   double precision :: delta_b, delta_h, ratio, ratio_tolerance, eps
   double precision :: t
   delta_b = abs(bR - bL)
   delta_h = 0.5d0*abs(bL + bR)

   if (delta_h < 1.e-6) then
      bath_sharpness_sensor = 1.d0 !Full correction if very shallow water
      return
   end if

   ratio = delta_b / delta_h
   ratio_tolerance = 0.01d0
   eps = 0.05d0
   if (ratio <= ratio_tolerance) then
      bath_sharpness_sensor = 1.d0
   else if (ratio > ratio_tolerance .and. ratio < ratio_tolerance + eps) then
      t = (ratio - ratio_tolerance)/eps
      bath_sharpness_sensor = -3.d0*t**2 + 2.d0*t**3 + 1.d0
   else
      bath_sharpness_sensor = 0.d0
   end if

end function bath_sharpness_sensor


!=============================================================================
!Subroutines for BBBD system
!=============================================================================
subroutine riemann_normal_HypRel_BBBD(mu, nv, grav, c_sq, &
    hl, ul, vl, wl, sgl, pl, pbl, bl, csq_l,&
    hr, ur, vr, wr, sgr, pr, pbr, br, csq_r, &
    fw, sw)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !When used in 2D:
    !  Input: mu and nv according to the normal direction (analogous to Geoclaw's solver for SWEs)
    !  Output: 7 f-waves with 7 components each (fw) and 7 wave speeds (sw)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !When used in 1D :
    !  Input: mu=2, nv =3, vl=vr=0
    !  Output: The 3rd component of all the f-waves must be discarded
    !          The 3rd f-wave and the 3rd wave speed must be discarded
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !${CLAW}/riemann/src/geoclaw_riemann_utils.f is probably a better home for this routine
    !I'll leave it here for now
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    implicit none

    !Input
    integer, intent(in) :: mu, nv
    double precision, intent(in) :: grav, c_sq
    double precision, intent(in) :: hl, ul, vl, wl, pl, bl, csq_l
    double precision, intent(in) :: hr, ur, vr, wr, pr, br, csq_r
    double precision, intent(in) :: sgl, sgr, pbl, pbr

    !Output
    double precision, intent(out) :: fw(7,7), sw(7)

    !Local
    double precision :: hlsq, hrsq, etaL, etaR
    double precision :: hbar, ubar, vbar, wbar, pbar, csqbar, CEbar
    double precision :: sgbar, pbbar
    double precision :: delb
    double precision :: lamb1bar, lamb7bar
    double precision :: R(7,7), R_inv(7,7)
    double precision :: alpha(7,1), total_fluct(7,1)
    double precision :: use_full_bath, bath_correction_term(3)

    double precision :: bath_sharpness_sensor !For stability term

    integer :: mw

    etaL = hL + bL
    etaR = hR + bR

    !Roe averages
    hbar = 0.5*(hl+hr)
    hlsq = dsqrt(hl)
    hrsq = dsqrt(hr)
    ubar = (hlsq*ul+hrsq*ur)/(hlsq+hrsq)
    vbar = (hlsq*vl+hrsq*vr)/(hlsq+hrsq)
    wbar = (hlsq*wl+hrsq*wr)/(hlsq+hrsq)
    sgbar = (hlsq*sgl+hrsq*sgr)/(hlsq+hrsq)
    pbar = (hlsq*pl+hrsq*pr)/(hlsq+hrsq)
    pbbar = (hlsq*pbl+hrsq*pbr)/(hlsq+hrsq)

    !Still need to verify if there is difference between taking
    !the average, min, or max
    csqbar = 0.5d0*(csq_l+csq_r)

    CEbar=dsqrt(grav*hbar+csqbar+pbar)

    !Roe speeds
    lamb1bar=ubar-CEbar
    lamb7bar=ubar+CEbar

    !Matrix of right eigenvectors
    R(1,:)= [1.d0,       1.d0,      0.d0, 0.d0, 0.d0, 0.d0, 1.d0]
    R(mu,:)=[lamb1bar,   ubar,      0.d0, 0.d0, 0.d0, 0.d0, lamb7bar]
    R(nv,:)=[vbar,       0.d0,      0.d0, 0.d0, 1.d0, 0.d0, vbar]
    R(4,:)= [wbar,       0.d0,      0.d0, 1.d0, 0.d0, 0.d0, wbar]
    R(5,:)= [sgbar,      0.d0,      1.d0, 0.d0, 0.d0, 0.d0,  sgbar]
    R(6,:)= [pbar+csqbar,-grav*hbar,0.d0, 0.d0, 0.d0, 0.d0, pbar+csqbar]
    R(7,:)= [pbbar,      0.d0,      0.d0, 0.d0, 0.d0, 1.d0,  pbbar]
    
    delb = br-bl
    !Computing total fluctuation to split in f-waves
    total_fluct(1,1)=hr*ur-hl*ul
    bath_correction_term(1) = pbbar*delb
    total_fluct(mu,1) = hr*ur**2.d0 +hr*pr+0.5d0*grav*hr**2.d0 &
                    - hl*ul**2.d0 -hl*pl-0.5d0*grav*hl**2.d0 &
                    + grav*hbar*delb
    total_fluct(nv,1)=hr*ur*vr-hl*ul*vl
    total_fluct(4,1)=hr*ur*wr-hl*ul*wl
    total_fluct(5,1)=hr*ur*sgr-hl*ul*sgl
    bath_correction_term(2) = csqbar*ubar*delb
    total_fluct(6,1)=hr*ur*(csq_r+pr)-hl*ul*(csq_l+pl)-csqbar*ubar*(etaR-etaL)
    bath_correction_term(3) = -6.d0*csqbar*ubar*delb
    total_fluct(7,1)=hr*ur*pbr-hl*ul*pbl

    use_full_bath = 1.d0!bath_sharpness_sensor(bl,br)
    total_fluct(mu,1) = total_fluct(mu,1) + use_full_bath*bath_correction_term(1)
    total_fluct(6,1) = total_fluct(6,1) + use_full_bath*bath_correction_term(2)
    total_fluct(7,1) = total_fluct(7,1) + use_full_bath*bath_correction_term(3)

    
    R_inv(:,1)=(1.d0/CEbar**2.d0)*[0.5d0*(grav*hbar+CEbar*ubar),&
    csqbar+pbar, -grav*hbar*sgbar, -grav*hbar*wbar, -grav*hbar*vbar,&
    -grav*hbar*pbbar,0.5d0*(grav*hbar-CEbar*ubar)]
    R_inv(:,mu)=[-0.5d0/CEbar, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.5d0/CEbar]
    R_inv(:,nv)=[0.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 0.d0]
    R_inv(:,4)= [0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0]
    R_inv(:,5)= [0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0,  0.d0]
    R_inv(:,6)=(1.d0/CEbar**2.d0)*[0.5d0,-1.d0,-sgbar,&
        -wbar,-vbar, -pbbar, 0.5d0]
    R_inv(:,7)=[0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0]

    alpha = matmul(R_inv,total_fluct)

    ! Computing f-waves
    do mw = 1, 7 !Replace by mwaves if necessary
        fw(:, mw) = alpha(mw,1) * R(:, mw)
    END do

    !Computing wave speeds
    !Einfeldt speeds for outer waves
    sw(1) = min(ul-sqrt(grav*hl+csqbar+pl),lamb1bar)
    sw(2:6) = ubar
    sw(7) = max(ur+sqrt(grav*hr+csqbar+pr),lamb7bar)

end subroutine riemann_normal_HypRel_BBBD

!=============================================================================
!Subroutines for DR system Duran and Richard 2025
!=============================================================================

subroutine riemann_1D_HypRel_DR(g, alphaLDR, c_sq, &
   hl, ul, wl, pl, sgl, bl, csq_l, &
   hr, ur, wr, pr, sgr, br, csq_r, &
   fw, sw)
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !Normal Riemann solver for Duran and Richard's system
   !Since the matrix of right eigenvectors is invertible, but its
   !analytical inverse is very ugly, we invert it numerically.
   !before I was using LU decomposition from lu_module, but now I
   !am switching lo LAPACK's dgesv routine for better everything
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   implicit none

   !Input
   ! integer :: mu, nv
   double precision :: g, alphaLDR, c_sq
   double precision :: hl, ul, wl, pl, sgl, bl, csq_l
   double precision :: hr, ur, wr, pr, sgr, br, csq_r

   !Output
   double precision :: fw(5,5), sw(5)

   !Local
   double precision :: hlsq, hrsq
   double precision :: h, u, w, p, sigma
   double precision :: saux(7) !Auxiliary variables
   double precision :: lamb_bar(5), lamb_tilde(5)
   double precision :: R(5,5), Rsolve(5,5) !R would get overriden during the linear solve
   double precision :: Rinv(5,5), alpha(5,1)
   double precision :: total_fluct(5,1), total_fluct_solve(5,1)
   double precision :: bath_correction_term(3)
   double precision :: denom_Rinv(2)
   double precision :: bath_sharpness_sensor !For stability term
   integer :: j


   integer :: mw, j_debug

   !Variables for LAPACK dgesv
   integer :: n, info, lda, ldb, nrhs
   integer, allocatable :: ipiv(:)
   !external declarations 
   external dgetrf, dgetrs

   integer :: M, RANK, LWORK
   double precision :: RCOND
   double precision, allocatable :: A(:,:), B(:,:), S(:),  WORK(:)
   double precision :: WORKQ(1)
   integer :: linear_solve_flag

   linear_solve_flag = 0 !Set to 1 if some part of the linear solve fails

   ! print*, "alphaLDR = ", alphaLDR
   !Roe averages
   h = 0.5d0*(hl+hr)
   hlsq = dsqrt(hl)
   hrsq = dsqrt(hr)
   u = (hlsq*ul+hrsq*ur)/(hlsq+hrsq)
   w = (hlsq*wl+hrsq*wr)/(hlsq+hrsq)
   p = (hlsq*pl+hrsq*pr)/(hlsq+hrsq)
   sigma = (hlsq*sgl+hrsq*sgr)/(hlsq+hrsq)

   !Still need to verify if there is difference between taking
   !the average, min, or max
   c_sq = 0.5d0*(csq_l+csq_r)

   saux(:) = [dsqrt(g*h+p+c_sq*alphaLDR), &
               g*h+alphaLDR*(c_sq*alphaLDR+p), &
               dsqrt(h)*w, &
               alphaLDR*(c_sq*alphaLDR+p),&
               0.d0,&
               0.d0,&
               sqrt(g*(alphaLDR-1.d0)/alphaLDR)]
   saux(5) = 3.d0*saux(3)*(g*h*alphaLDR+saux(4))
   saux(6) = saux(1)*sigma*(g*h+saux(4))

   !Roe speeds
   lamb_tilde(1) = u-saux(1)
   lamb_tilde(2) = u-sqrt(g*h*(alphaLDR-1.d0)/alphaLDR)
   lamb_tilde(3) = u
   lamb_tilde(4) = u+sqrt(g*h*(alphaLDR-1.d0)/alphaLDR)
   lamb_tilde(5) = u+saux(1)

   !Einfeldt speeds
   lamb_bar(1) = min(lamb_tilde(1), ul-dsqrt(g*hl+csq_l*alphaLDR+pl))
   lamb_bar(2) = min(lamb_tilde(2), ul-dsqrt(g*hl*(alphaLDR-1.d0)/alphaLDR))
   lamb_bar(3) = u
   lamb_bar(4) = max(lamb_tilde(4), ur+dsqrt(g*hr*(alphaLDR-1.d0)/alphaLDR))
   lamb_bar(5) = max(lamb_tilde(5), ur+dsqrt(g*hr+csq_r*alphaLDR+pr))
   !Roe matrix of right eigenvectors
   R(1,:) = [saux(2)*(saux(5) + saux(6)),0.d0,1.d0,0.d0,saux(2)*(saux(5) - saux(6))]
   R(2,:) = [(u - saux(1))*saux(2)*(saux(5) + saux(6)),0.d0,u,0.d0,(u + saux(1))*saux(2)*(saux(5) - saux(6))]
   R(3,:) = [(w*((-1 + 3*alphaLDR)*g*h + 2*saux(4))*(saux(5) + saux(6)))/2.,saux(7)/2.,-0.5*w,&
               -0.5*saux(7),(w*((-1 + 3*alphaLDR)*g*h + 2*saux(4))*(saux(5) - saux(6)))/2.]
   R(4,:) = [(alphaLDR*c_sq + p)*saux(2)*(saux(5) + saux(6)),0.d0,-(g*h),0.d0,(alphaLDR*c_sq + p)*saux(2)*(saux(5) - saux(6))]
   R(5,:) = [(sigma*saux(2) + 3*alphaLDR*saux(1)*saux(3))*(saux(5) + saux(6)),1.d0,sigma,1.d0,&
               (-3*alphaLDR*saux(1)*saux(3) + sigma*(g*h + saux(4)))*(saux(5) - saux(6))]

   ! do j = 1, 5 !Normalizing columns of R for better conditioning
   !    colmax = maxval(abs(R(:,j)))
   !    if (colmax > 0d0) then
   !       R(:,j) = R(:,j) / colmax
   !    endif
   ! end do

   Rsolve = R !Copying R to Rsolve since R would get overriden during the linear solve

   
   !Computing total fluctuation to split in f-waves
   total_fluct(1,1)=hr*ur-hl*ul
   total_fluct(2,1)=hr*ur**2+hr*pr+0.5d0*g*hr**2-&
                  hl*ul**2-hl*pl-0.5d0*g*hl**2+&
                  g*h*(br-bl)
   total_fluct(3,1)=ur*hr*wr-ul*hl*wl+&
                     0.5d0*((1.d0-alphaLDR)/alphaLDR)*((hr**1.5)*sgr-(hl**1.5)*sgl)
   bath_correction_term(1) = 0.75d0*((alphaLDR-1.d0)/alphaLDR)*g*sigma*dsqrt(h)*(hr-hl) !Non-conservative term
   total_fluct(4,1)= hr*ur*(pr+c_sq*alphaLDR)-hl*ul*(pl+c_sq*alphaLDR) 
   bath_correction_term(2) = -c_sq*alphaLDR*u*(hr-hl) !Non-conservative term
   total_fluct(5,1)=hr*ur*sgr-2.d0*dsqrt(hr)**3.d0*wr-&
                     (hl*ul*sgl-2.d0*dsqrt(hl)**3.d0*wl)

   !Call sharpness sensor to avoid instability due to bathymetry
   total_fluct(3,1) = total_fluct(3,1) + bath_correction_term(1)
   total_fluct(4,1) = total_fluct(4,1) + bath_sharpness_sensor(bl,br)*bath_correction_term(2)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !Try building Rinv, and check if we have good conditioning (really looking for a zero denominator)
   denom_Rinv(1) = saux(5)+ saux(6)
   denom_Rinv(2) = saux(5)- saux(6)
   ! if (abs(denom_Rinv(1))>1.e-15 .and. abs(denom_Rinv(2))>1.e-15) then
   if (.false.) then
      !Building Rinv analytically
      Rinv(:,1) = [(g*h + u*saux(1))/(2.*saux(1)**2*saux(2)*(saux(5) + saux(6))),&
         (g*h*w*((1 - 3*alphaLDR)*g*h - 2*saux(4)) - 3*alphaLDR*u*saux(1)**2*saux(3)*saux(7) +&
         saux(2)*((alphaLDR*c_sq + p)*w - (alphaLDR*c_sq + g*h + p)*sigma*saux(7)))/(2.*saux(1)**2*saux(2)*saux(7)),&
         1 - (g*h)/saux(1)**2,&
         (g*h*w*((-1 + 3*alphaLDR)*g*h + 2*saux(4)) - 3*alphaLDR*u*saux(1)**2*saux(3)*saux(7) +&
         saux(2)*(-((alphaLDR*c_sq + p)*w) - (alphaLDR*c_sq + g*h + p)*sigma*saux(7)))/(2.*saux(1)**2*saux(2)*saux(7)),&
         (g*h - u*saux(1))/(2.*saux(1)**2*saux(2)*(saux(5) - saux(6)))]
      Rinv(:,2) = [-0.5*1/(saux(1)*saux(2)*(saux(5) + saux(6))),(3*alphaLDR*saux(3))/(2.*saux(2)),0.d0,&
         (3*alphaLDR*saux(3))/(2.*saux(2)),1/(2*saux(1)*saux(2)*saux(5) - 2*saux(1)*saux(2)*saux(6))]
      Rinv(:,3) = [0.d0,1/saux(7),0.d0,-(1/saux(7)),0.d0]
      Rinv(:,4) = [1/(2.*saux(1)**2*saux(2)*(saux(5) + saux(6))),-0.5*(w*((-1 + 3*alphaLDR)*g*h + saux(2) +&
         2*saux(4)))/(saux(1)**2*saux(2)*saux(7)),&
         -saux(1)**(-2),(w*((-1 + 3*alphaLDR)*g*h + saux(2) + 2*saux(4)))/(2.*saux(1)**2*saux(2)*saux(7)),&
         1/(2.*saux(1)**2*saux(2)*(saux(5) - saux(6)))]
      Rinv(:,5) = [0d0,0.5d0,0d0,0.5d0,0d0]

      alpha = matmul(Rinv, total_fluct)
      total_fluct(:,1) = alpha(:,1)
   else
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Try to solve numerically (should skip directly to least-squares, but try LU first anyway)
      !Linear solve R*alpha = total_fluct
      !First, skip linear solve if total_fluct is close to zero
      !Then try using LU (fast but less robust)
      !If that fails, use least-squares with SVD (slower but more robust)
      !If that fails, set total_fluct to zero and warn user
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! print *, "Warning: R matrix poorly conditioned in riemann_1D_HypRel_DR. Attempting numerical solve."
      ! total_fluct_solve(:,1) = total_fluct(:,1) !Backing up total_fluct since it gets overwritten
      if(maxval(abs(total_fluct)) < 1.d-15) then
         total_fluct(:,1) = 0.d0
      else
         !Proceed to linear solve
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !Using LU decomposition from LAPACK dgesv
         n = 5
         nrhs = 1
         lda = 5
         ldb = 5
         allocate(ipiv(n))
         call dgetrf(n, n, Rsolve, lda, ipiv, info)
         if (info /= 0) then
            !If factorization failed, set total_fluct to zero, but warn user
            ! print *, "Warning: LAPACK dgesv failed to factorize matrix in riemann_1D_HypRel_DR."
            ! total_fluct(:) = 0.d0
            linear_solve_flag = 1
         else 
            !Solve the system
            call dgetrs('N', n, nrhs, Rsolve, lda, ipiv, total_fluct, ldb, info)
            if (info /= 0) then
               !If solve failed, set total_fluct to zero, but warn user
               ! total_fluct(:) = 0.d0
               ! print *, "Warning: LAPACK dgesv failed to solve for wave strengths in riemann_1D_HypRel_DR."
               linear_solve_flag = 1
            end if
         end if
         deallocate(ipiv)
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !If LU failed, use SVD-based least squares
         if (linear_solve_flag == 1) then
            total_fluct = total_fluct_solve !Restore total_fluct before attempting SVD solve
            M=5
            RCOND = 1.e-12
            allocate(A(5,5), B(5,1), S(5))
            A = R
            B(:,1) = total_fluct(:,1)    
            ! --- Workspace query ---
            LWORK = -1

            call dgelss(M, N, NRHS, A, LDA, B, LDB, S, RCOND, RANK, WORKQ, LWORK, INFO)

            LWORK = int(WORKQ(1))
            allocate(WORK(LWORK))

            ! Restore RHS because DGELSS overwrote B during the query
            B(:,1) = total_fluct(:,1)

            ! --- Actual solve ---
            call dgelss(M, N, NRHS, A, LDA, B, LDB, S, RCOND, RANK, WORK, LWORK, INFO)

            if (INFO /= 0) then
               print *, "DGELSS failed, info=", INFO
               total_fluct(:,1) = 0d0
            else
               total_fluct(:,1) = B(:,1)
            end if
            ! print *, "Singular values of R:", S

            !Check if total_fluct is greater than a small threshold
            ! if (maxval(abs(total_fluct)) > 1.d-14) then
            !    print*, 'total fluct not zeroed'
            ! end if
            deallocate(A, B, S, WORK)
         end if
      end if
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   end if
   ! Computing f-waves
   do mw = 1, 5 !Replace by mwaves if necessary
      fw(:, mw) = total_fluct(mw,1) * R(:, mw)
   END do
   
   !Computing wave speeds
   sw(1) = lamb_bar(1)
   sw(2) = lamb_bar(2)
   sw(3) = lamb_bar(3)
   sw(4) = lamb_bar(4)
   sw(5) = lamb_bar(5)

end subroutine riemann_1D_HypRel_DR
