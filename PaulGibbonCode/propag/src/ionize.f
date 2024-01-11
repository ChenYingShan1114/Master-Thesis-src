
c     ================================================================
c     
c     Ionization model
c     
c     ================================================================

      subroutine ioniz(dte,tim)
      include 'dynafoc.h'
      complex eem
      
      parameter (eryd=13.61)
      parameter (niter=20)    ! max # iterations
      parameter (err_max=0.02)   ! max fractional density error allowed
      real eion(maxspc)
      real*8 alpha(maxspc), rhonew(0:maxspc), rhohalf(0:maxspc),
     : rhoold(0:maxspc), delta_n(maxspc)
      integer noconv
      data noconv /0/
      save noconv

      call cputime(tim1,1)

      
c     initial  gas density is rho0

      alpham=1.e-20

      if (mioniz.eq.1) then


c     Ionisation energies in Rydbergs (13.61 eV)
         if (gas .eq. 'He') then
c     helium
            eion(1)= 24.59/eryd
            eion(2)= 54.418/eryd
            alpha(1) = 0.
            alpha(2) = 0.
         else
c     hydrogen
            eion(1)=1.
            alpha(1) = 0.
            alpha(2) = 0.
         endif

         
c     compute ion density from Keldysh ionization model 

c     wa0 is the atomic frequency normalised to w0
c     Ea0 is the atomic E-field, normalised to m.w0.c/e
c     rhoi(i,k,0) is the local density of neutrals (unionized gas)
         
         do k=1,nz
            z = k*dz
            do i=0,nr
               elas = abs(a(i,k)) ! ignore EM phase


c     ionisation rates
               do jspec=1,nion
                  alpha(jspec) = 4.*wa0*eion(jspec)**2.5
     :*Ea0/elas * exp( -2./3.*eion(jspec)**1.5*Ea0/elas )

                  alpham = max(alpha(jspec),alpham)
               end do

c     1st guess for ion densities

               do jspec=0,nion
                  rhonew(jspec) = rhoi(i,k,jspec)
               end do

               it = 1
               err_n = 1.0

               do while (it .lt. niter .and. err_n .gt. err_max) 

c     1/2 steps   
                  do jspec=0,nion
                     rhohalf(jspec)=0.5*(rhonew(jspec)+rhoi(i,k,jspec))
                     rhoold(jspec) = rhonew(jspec)
                  end do

c     1st ion
                  rhonew(1) = rhoi(i,k,1)
     :                 + dte*alpha(1)*rhohalf(0)
     :                 - dte*alpha(2)*rhohalf(1)

c     if (rhoi(i,k,1).lt.0) then
c     write (*,*) '***** -ve n1 ****** ',i*dr,z
c     :,' n1: ',rhoi(i,k,1),' n2 ',rhoi(i,k,2),' n0',rhoi(i,k,0)
c     :,' a1: ',alpha(1),' a2: ',alpha(2)
c     endif

c     prevent going -ve and from exceeding gas density
                  rhonew(1) = min(rhoatom(i,k), max(0.,rhonew(1)))

                  if (nion .gt. 1) then
c     2nd ion - turn this into a loop for higher Z gases
                     rhonew(2) = rhoi(i,k,2)
     :                    + dte*alpha(2)*rhohalf(1)

c     clamp to local gas density
                     rhonew(2) = min(rhoatom(i,k),rhonew(2))
                  else
                     rhonew(2) = 0.
                  endif

c     new density of remaining neutrals
                  rhonew(0) = rhoatom(i,k) - rhonew(1) - rhonew(2)
                  rhonew(0) = max(rhonew(0),0.)

c  find change in iterations
                  err_n = 0.
                  do js=1,nion
                     delta_n(js) = (rhonew(js) - rhoold(js))/rho0
                     err_n = max(err_n,delta_n(js))
                  end do
                  it = it + 1

               end do

               if (err_n .gt. err_max) then
                  noconv = noconv+1
                  write (11,*) 
     : 'Rate equation not converged - count:',noconv,' error:',err_n
c                  write (6,*) ' z',z,' r',dr*i,' n0:'
c     :,rhonew(0),' n1:',rhonew(1),' n2:',rhonew(2),' error:',err_n 
               endif

c    update densities with final iteration
               do jspec=0,nion
                  rhoi(i,k,jspec) = rhonew(jspec)
               end do
              
            end do           
         end do



      else if (mioniz.eq.2) then

c     preplasma with tanh borders
         alpham=1.
         zmid = zedge+zplas/2.
c     zbord = zplas/10.
         z1 = zdens(nz) - (zedge-zbord)
         rhos = rho0*0.5*(1. + tanh(z1/zbord))
         z2 = zdens(nz) - (zend+zbord)
         rhof = rho0 - rho0*0.5*(1. + tanh(z2/zbord))

c     clamp at min density (rho0/1e3)
         rhomin = rho0/1.e3
         do i=0,nr
            if (zdens(nz).lt.zmid) then
               rhoi(i,nz,1) = max(rhos,rhomin)
            else
               rhoi(i,nz,1) = max(rhof,rhomin)
            endif
         end do

      else if (mioniz.eq.3) then

c     preplasma:   tanh borders in z
c     density channel created by ion motion
         alpham=1.
         zmid = zedge+zplas/2.
c     zbord = zplas/10.
         z1 = zdens(nz) - (zedge-zbord)
         rhos = rho0*0.5*(1. + tanh(z1/zbord))
         z2 = zdens(nz) - (zend+zbord)
         rhof = rho0 - rho0*0.5*(1. + tanh(z2/zbord))

         do i=0,nr
            if (zdens(nz).lt.zmid) then
               rhoi(i,nz,1) = max(rhos*frion(i),1.e-20)
            else
               rhoi(i,nz,1) = max(rhof*frion(i),1.e-20)
            endif
         end do

      else 
         do i=0,nr
            if (zdens(nz).ge.zedge .and. zdens(nz).lt.zend) then
               rhoi(i,nz,1) = rho0
            else
               rhoi(i,nz,1)=0.
            endif
         end do
      endif 
      
c     ionization timestep
      dtioniz = 0.2/alpham

      call cputime(tim2,1)
      tim = tim + tim2-tim1
      
      end







