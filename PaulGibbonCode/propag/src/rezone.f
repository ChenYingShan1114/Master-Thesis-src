

c     ================================================================
c     
c     Rezone plasma grid quantities
c     
c     ================================================================


      subroutine rezone(dte)
      include 'dynafoc.h'

c     rezone ion densities: 
c     move plasma back through grid with Dz=c.dt

c     ion density
      do k=1,nz-1
         do i=0,nr+1
            do jspec=0,nion
               rhoi(i,k,jspec) = rhoi(i,k+1,jspec)
            end do
            rhoe(i,k) = rhoe(i,k+1)  
            rhoatom(i,k) = rhoatom(i,k+1)
         end do
      end do
      

c     em spectral probe position

      if (dt*itime.gt.tprobe) then
         kprobe = kprobe-1
      endif

      
      if (mioniz.eq.0) then
         do i=0,nr
c     fresh preionized plasma at RHB
            rhoi(i,nz,1) = rho0
         end do
      else         
         do i=0,nr
            do jspec=1,nion
               rhoi(i,nz,jspec) = 0.0

             if (zedge .lt. zl+3*zbord .and. zend .gt. zl-3*zbord) then
c     fresh vacuum or gas at RHB
c     longitudinal gas density profile:    tanh leading and trailing edges
               zmid = zedge+zplas/2.
               z1 = zl - (zedge-zbord)
               rhos = rho0*0.5*(1. + tanh(z1/zbord))
               z2 = zl - (zend+zbord)
               rhof = rho0 - rho0*0.5*(1. + tanh(z2/zbord))

c     clamp at min density (rho0/1e3)
               rhomin = rho0/1.e8
               if (zl.lt.zmid) then
                  rhoatom(i,nz) = max(rhos,rhomin)
               else
                  rhoatom(i,nz) = max(rhof,rhomin)
               endif

               rhoi(i,nz,0) = 0.
            else
               rhoatom(i,nz) = 0.
               rhoi(i,nz,0) = 0.
            endif

            end do
         end do
      endif
      

c     update plasma boundary
      zedge = zedge - dte
c     shift RH boundary if ionization on
      if (mioniz.ge.1) zend = zend-dte

      end



