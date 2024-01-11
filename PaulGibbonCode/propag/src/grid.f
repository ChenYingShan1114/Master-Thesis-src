

c =============================================================
 

      subroutine grid
      include 'dynafoc.h'
      
      if (igrid.eq.2) then
c  uniform action-weighting - gaussian spacing
        alpha = rl/sigma0/log(2.*nr)
        do j=1,nr+1
          r(j) = -alpha*sigma0*log(1.- (j-0.5)/nr)
        end do
        
      else
c  uniform grid
        do j=0,nr+1
          r(j) = dr*j
        end do
      endif

      end
