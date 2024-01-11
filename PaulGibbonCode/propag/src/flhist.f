
  
c  ================================================================== 
c
c     Fluid histories
c
c  ================================================================== 

      subroutine flhist
      include 'dynafoc.h'
      real work1(0:nrmax),rc(0:nrmax),work2(0:nrmax)
      data itc/1/
      save itc

c  Average electron density and ES energy
      emass=0.
      ese = 0.
      km=max(nz/2,1)
      do i=1,nr-1
        emass = emass + rhoe(i,km)*r(i)*dr
        ese = ese + 0.5*Er(i,km)**2*r(i)*dr
      end do
      edena(itc) = emass*2/rl**2
      Ues(itc) = ese

c  Electrostatic energy
      ese = 0.
    
      itc = itc + 1

      end


















