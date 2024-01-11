
c     ================================================================
c     
c     Ionization model
c     
c     ================================================================

      subroutine ioniz(dte,tim)
      include 'dynafoc.h'
      complex eem
      
      parameter (eryd=13.61)
      real eion(maxspc)
      real alpha(maxspc)


      call cputime(tim1,1)

      
c     initial  gas density is rho0

      alpham=1.e-20

      if (mioniz.eq.1) then
c     Ionisation energies in Rydbergs (13.61 eV)
         if (gas .eq. 'He') then
c     helium
            eion(1)= 24.59/eryd
            eion(2)= 54.418/eryd
         else
c     hydrogen
            eion(1)=1.
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
     :*   Ea0/elas * exp( -2./3.*eion(jspec)**1.5*Ea0/elas )

                  alpham = max(alpha(jspec),alpham)
               end do


c     1st ion
               rhoi(i,k,1) = rhoi(i,k,1)
     :              + dte*alpha(1)*rhoi(i,k,0)
     :              - dte*alpha(2)*rhoi(i,k,1)

c               if (rhoi(i,k,1).lt.0) then
c                  write (*,*) '***** -ve n1 ****** ',i*dr,z
c     :,' n1: ',rhoi(i,k,1),' n2 ',rhoi(i,k,2),' n0',rhoi(i,k,0)
c     :,' a1: ',alpha(1),' a2: ',alpha(2)
c               endif

c     prevent going -ve from exceeding gas density
               rhoi(i,k,1) = min(rhoatom(i,k), max(0.,rhoi(i,k,1)) )

               if (nion .gt. 1) then
c     2nd ion - turn this into a loop for higher Z gases
                  rhoi(i,k,2) = rhoi(i,k,2)
     :                 + dte*alpha(2)*rhoi(i,k,1)

c     clamp to local gas density
                  rhoi(i,k,2) = min(rhoatom(i,k),rhoi(i,k,2))

               endif

c  new density of remaining neutrals
               rhoi(i,k,0) = rhoatom(i,k) - rhoi(i,k,1) - rhoi(i,k,2)
               rhoi(i,k,0) = max(rhoi(i,k,0),0.)


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







