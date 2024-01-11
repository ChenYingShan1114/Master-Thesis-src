

      
c     =============================================================
c     
c     Compute electron density due to ponderomotive force:
c     adiabatic (static or neutralising) electrons
c     
c     =============================================================

      subroutine density(tim)
      include 'dynafoc.h'

      real rhetot

      call cputime(tim1,v)




      if (melec.eq.1) then

c     adiabatic electrons

         do k=1,nz
            do i=1,nr
               dri = 0.5*(r(i+1)-r(i-1))
c  total unperturbed electron density summed over species
               rhetot=0.
               do jspec = 1,nion
                  rhetot = rhetot + jspec*rhoi(i,k,jspec)
               end do
               rhe = rhetot  
     :       + 2 * ( ( gamma(i+1,k) - 2*gamma(i,k) + gamma(i-1,k) )/dri**2
     :       + ( gamma(i+1,k) - gamma(i-1,k) )/(2*r(i)*dri) )

c  prevent -ve electron densities
               rhoe(i,k) = max(0.,rhe)            
            end do

c     special treatment for r=0
            dri = r(1)
c  total unperturbed electron density summed over species
               rhetot=0.
               do jspec = 1,nion
                  rhetot = rhetot + 1.*jspec*rhoi(0,k,jspec)
               end do
            rhe = rhetot 
     :           + 2 * 2*( gamma(1,k) - gamma(0,k) )/dri**2
            rhoe(0,k) = max(0.,rhe)
         end do


      else 

c     melec=0 - uniform density (hydrogen): no ponderomotive effects

         if (gas .eq. 'He') then
            do k=1,nz
               do i=0,nr
                  rhoe(i,k) = rhoi(i,k,1) + 2 * rhoi(i,k,2)
               end do
            end do
         else
            do k=1,nz
               do i=0,nr
                  rhoe(i,k) = rhoi(i,k,1)
               end do
            end do
         endif
         
      endif

      do i=0,nr
         rhoe(i,0) = rhoe(i,1)
      end do


      call cputime(tim2,v)
      tim = tim + tim2-tim1

      end










