c     Wave limiter for the hyperbolic relaxation systems: copy of amrclaw/src/2d/inlinelimiter.f with
c     'strength' limiting for EDC (mbHSGN) and BBBD (HSGN): for a wave with a nonzero h-component
c     (acoustic and contact waves, whose eigenvectors have 1 in the h-row) the smoothness ratio is the
c     ratio of the h-components, i.e. of the wave strengths beta, instead of the dot product over all
c     components. The dot product is dominated by the hp-row, (p+c^2)*beta, which mixes the
c     cell-to-cell variation of c^2 = alpha^2 g H into the ratio. Waves without an h-component
c     (single-component passive waves) and all other systems use the usual dot product.
      subroutine limiter(maxm,meqn,mwaves,mbc,mx,wave,s,mthlim)
      use bouss_module, only: boussEquations
      implicit double precision (a-h,o-z)
      dimension mthlim(mwaves)
      dimension wave(meqn, mwaves, 1-mbc:maxm+mbc)
      dimension    s(mwaves, 1-mbc:maxm+mbc)
      logical strength
      strength = (boussEquations == -1 .or. boussEquations == -2)
      do 200 mw=1,mwaves
         if (mthlim(mw) .eq. 0) go to 200
         dotr = 0.d0
         hprev = wave(1,mw,0)
         hcur = wave(1,mw,0)
         do 190 i = 0, mx+1
            wnorm2 = 0.d0
            dotl = dotr
            dotr = 0.d0
            do 5 m=1,meqn
               wnorm2 = wnorm2 + wave(m,mw,i)**2
               dotr = dotr + wave(m,mw,i)*wave(m,mw,i+1)
    5          continue
            hprev = hcur
            hcur = wave(1,mw,i)
            if (i.eq.0) go to 190
            if (wnorm2.eq.0.d0) go to 190
            if (strength .and. hcur .ne. 0.d0) then
               if (s(mw,i) .gt. 0.d0) then
                  r = hprev / hcur
               else
                  r = wave(1,mw,i+1) / hcur
               endif
            else
               if (s(mw,i) .gt. 0.d0) then
                  r = dotl / wnorm2
               else
                  r = dotr / wnorm2
               endif
            endif
            go to (10,20,30,40,50) mthlim(mw)
   10       continue
            wlimitr = dmax1(0.d0, dmin1(1.d0, r))
            go to 170
   20       continue
            wlimitr = dmax1(0.d0, dmin1(1.d0, 2.d0*r), dmin1(2.d0, r))
            go to 170
   30       continue
            wlimitr = (r + dabs(r)) / (1.d0 + dabs(r))
            go to 170
   40       continue
            c = (1.d0 + r)/2.d0
            wlimitr = dmax1(0.d0, dmin1(c, 2.d0, 2.d0*r))
            go to 170
   50       continue
            wlimitr = r
            go to 170
  170       continue
            do 180 m=1,meqn
               wave(m,mw,i) = wlimitr * wave(m,mw,i)
  180          continue
  190       continue
  200    continue
      return
      end
