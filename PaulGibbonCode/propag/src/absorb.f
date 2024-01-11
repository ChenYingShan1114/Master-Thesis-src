

c =============================================================


      subroutine absorb
      include 'dynafoc.h'
c  Absorber on RH boundary:  nu=1 on boundary

c  renorm absorption rates: ensure xnumax*dt<<1

      xnumax = xnumax/dt
      xnumin = xnumin/dt

      if (mabs.eq.1) then

c  exponential with flat top
        sabs = -nabs*dr/log(xnumin)
        nflat = nr-nabs/2
        rfl = rl - dr*nabs/2
        do i=0,nr+1
          if (i.lt.nflat-nabs/2) then 
            xnu(i) = 0.
          else if (i.ge.nflat) then
            xnu(i) = xnumax
          else
            xnu(i) = amax1(xnumax*exp((r(i)-rfl)/sabs),0.)
          endif
        end do
      else if (mabs.eq.2) then
c  quadratic
        ra = r(nr-nabs)
        s = rl-ra
        do i=0,nr
          if (i.lt.nr-nabs) then
            xnu(i) = 0.
          else
             xnu(i) = amax1((r(i)-ra)**2/s**2*xnumax,0.)
c             xnu(i) = xnumax
          endif
        end do

      else
c  none
        do i=1,nr
          xnu(i)=0.
        end do
 
      endif


      end
