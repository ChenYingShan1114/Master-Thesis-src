       

c =============================================================

      subroutine launch
      include 'dynafoc.h'


c  ilas = 1-3:  launch pulse from LH boundary - fixed shape

      t = dt*itime

      do i=1,nr-1
        rc = r(i)
        if (ilas.eq.1) then
c        a(i,1) = a0*exp(-rc**2/2/sigma**2 - (t-t0)**2/2/tau**2)
c  Penetrante and Bardsley pulse
          argr = amin1( rc**2/sigma**2,20.)
          argt = amin1( 2*(t-t0)/tau,20.)
          a(i,1) = a0*exp(-argr)*1./cosh(argt)

        else if (ilas.eq.2) then
c  Gaussian x, constant with linear rise-time
          argr = amin1( rc**2/sigma**2,20.)
          if (t.lt.tau) then
            a(i,1) = a0*sqrt(t/tau)*exp(-argr)
          else 
            a(i,1) = a0*exp(-argr)
          endif

         else if (ilas.eq.3) then
c  Gaussian x and t
           argr = amin1( rc**2/2./sigma**2,20.)
           argt = amin1( (t-t0)**2/2./tau**2,20.)
           a(i,1) = a0*exp(-argr-argt)

        endif

      end do
      a(0,1) = a(1,1)
      do i=1,nr
        epond(i,1) = 0.25*((abs(a(i+1,1)))**2
     :                     - (abs(a(i-1,1)))**2)/2./dr
      end do

      epond(0,1) = 0.

      end
