
c  ==================================================================
c
c      Initial profile graphs
c
c  ================================================================== 

      subroutine grinp
      include 'dynafoc.h'
      real work1(0:nrmax),rc(0:nrmax),reno(0:nrmax),work2(0:nrmax)
      rlno = rl*rfac

c  laser amplitude
      do i=0,nr
        reno(i) = r(i)*rfac
      end do

c  Fixed grid quantities

      call grxy(reno,xnu,nr,10,1,1
     :,'      r        ','     xnu       ','xnua            ')
c     :,0.,rlno,rlno/5.,0.,2*a0,.25*a0)
      end
