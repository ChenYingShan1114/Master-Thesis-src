
c  ================================================================
c
c     Rezone particles
c
c  ================================================================


      subroutine prezone
      include 'dynafoc.h'

      do k=1,nzp-1
c  re-label particles
        ifirst = nse*(k-1)
        do l=1,nse
          ro(ifirst+l) = ro(ifirst+nse+l)
          rn(ifirst+l) = rn(ifirst+nse+l)
          qc(ifirst+l) = qc(ifirst+nse+l)
          vr(ifirst+l) = vr(ifirst+nse+l)
          ur(ifirst+l) = ur(ifirst+nse+l)
        end do

      end do
c  insert fresh slice at k=nzp
      ifirst = nse*(nzp-1)+1
      call denprof(ifirst,nse)
      call v1q(ifirst,nse,vte)
      call ixq(ifirst,nse)
      call calcu(ne)
c  recalculate density
      call eden
c  neutralise new slice
      if (ni.eq.0) call addneut(nzp)

      end
