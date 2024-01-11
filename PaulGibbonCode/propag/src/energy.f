
c     ================================================================
c     
c     Particle energies
c     
c     ================================================================


      subroutine energy(i1,n)
      include 'dynafoc.h'
      data ntc/1/
      save ntc
      real mass
      real*8 demax

      complex yes
      if (i1.eq.1) then
         mass=me
         uth(ntc)=0.
         udr(ntc)=0.
      else
         mass=mi
      endif

c     potential energy

      ese=0.
      do k=1,nzp
         do i=0,nr
            ese = ese + 0.5*er(i,k)**2*dr
         end do
      end do

      ues(ntc)=ese

c     total
      utot(ntc)=ues(ntc)


c     Max electron density
      demax=0.
      do k=1,nzp
         do i=0,nr
            demax=max(demax,abs(rhoe(i,k)))
         end do
      end do
      denem(ntc) = demax
c     Density on axis - ave. of 1st 3 cells
      km = max(t0/dzp,1.)
      dene0(ntc) = abs(rhoe(0,km)+rhoe(1,km)+rhoe(2,km))/3.
      ntc=ntc+1

      end
