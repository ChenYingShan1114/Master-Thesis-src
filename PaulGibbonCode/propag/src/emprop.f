      
c     =============================================================
c     
c     EMPROP - code kernel
c     
c     EM wave propagation using nonlinear Klein-Gordon equation (cf: NLSE)
c     2nd order implicit scheme (3 iterations - nitem - is usually enough)
c     - includes: relativistic self-focussing (saturated nonlinearity)
c     diffraction
c     refraction due to electron density modification via
c     ionisation and ponderomotive expulsion (cavitation)
c     - conserves wave action
c     
c     =============================================================
      
      subroutine emprop(dte,tim)
      include 'dynafoc.h'

      complex an,aip1,ai,aim1,aic,a2,d2a2,dadr,acurr,anold
     :     ,yr,yc1,yc2,yc3
      dimension an(0:nrmax,0:nzmax)
      real acmax, drhomax
      parameter (err_max = 0.01)
      integer ncount, nit
      data nit /4/
      save nit
c     data ncount /0/
c     save noconv

      call cputime(tim1,v)
    
      ncount = 0

c     first guess:  an(i,k)=a(i,k)
      do k=1,nz
         do i=0,nr+1
            an(i,k) = a(i,k)
         end do
      end do
      acmax = (0.,0.)
      drhomax = 0.
      if (nr.eq.1) then
         yc1 = (0.,0.)
         yc2 = (0.,0.)
      else
         yc1 = yi*0.5*dte
         yc2 = yi*0.25*dte
      endif
      
c     index of start of plasma
      kedge = zedge/dz
      kend = zend/dz
!      write(*,*) zedge,zend,kedge,kend
      err_n=0.

      do k=1,nz

c     plasma-vacuum boundary
         if (k.ge.kedge .and. k.le.kend) then
            yc3 = yi*0.5*dte
         else 
            yc3 = (0.,0.)
         endif


c     iteration loop is ALL r because of 2nd derivatives
         
         do it=1, nit 

            do i=1,nr

c     intermediates
               aip1 = 0.5*( an(i+1,k) + a(i+1,k) )
               ai = 0.5*( an(i,k) + a(i,k) )
               aim1 = 0.5*( an(i-1,k) + a(i-1,k) )
               aic = conjg(ai)
               a2  = ai*aic

               anold = an(i,k)
               if (mrel.eq.1) then  
                  gam = sqrt( 1.+fpol*abs(a2) )
               else
                  gam=1.
               endif

c     total unperturbed electron density summed over species
               rhetot=0.
               do jspec = 1,nion
                  rhetot = rhetot + jspec*rhoi(i,k,jspec)
               end do
c               drhomax = max(drhomax,rhetot-rhoe(i,k)/gam)
c     refractive index depends on vacuum/plasma dispersion relation 
               acurr = (rhetot - rhoe(i,k)/gam)*ai
c               acurr =  - rhoe(i,k)/gam*ai
c     grid spacing
               dri = 0.5*( r(i+1) - r(i-1) )
               

               an(i,k) = a(i,k)
c     Diffraction
     :              - yc1*( aip1 - 2*ai + aim1 )/dri**2
     :              - yc2*( aip1 - aim1 )/r(i)/dri
c     Self focussing and density
     :              - yc3*acurr
c     RH boundary absorption
     :              - dte*xnu(i)*ai
               err1 = (an(i,k)-abs(anold))/a0
               err_n= max(err_n,err1)
            end do

            
c     update for a(0) at r=0

c     intermediates
            aip1 = 0.5*( an(1,k) + a(1,k) )
            ai = 0.5*( an(0,k) + a(0,k) )
            aic = conjg(ai)
            a2  = ai*aic
            anold = an(0,k)

            if (mrel.eq.1) then  
               gam = sqrt( 1.+fpol*abs(a2) )
            else
               gam=1.
            endif

c     total unperturbed electron density summed over species
            rhetot=0.
            do jspec = 1,nion
               rhetot = rhetot + jspec*rhoi(0,k,jspec)
            end do

c     refractive index depends on vacuum/plasma dispersion relation 
            acurr = (rhetot - rhoe(0,k)/gam)*ai
c            acurr = -rhoe(0,k)/gam*ai

c     grid spacing
            dri = r(1) 
c            write (*,*) rhetot,rhoe(0,k),acurr
            an(0,k) = a(0,k)
c     Diffraction
     :           - 2*yc1*( aip1 - ai )/dri**2
c     Self focussing and density
     :           - yc3*acurr

            err1 = (an(0,k)-abs(anold))/a0
            err_n = max(err1,err_n)

         end do
         if (err_n .gt. err_max) then
            ncount = ncount+1
         endif
      end do


      if (ncount .gt. 1) then
         write(11,*) 'EM amplitude not converged at timestep ',itime
     :        ,'- count: ',ncount,
     :        ' error: ',err_n
         if (ncount.gt.nz/4) then
            nit = min(nit + 1,nitem)
            write (11,*) 'Increasing # iterations to ',nit
         endif
      endif

c     reset a; compute gamma at n+1
      do k=1,nz
         do i=0,nr
            a(i,k) = an(i,k)   
         end do

c     continuous at r=rl
         a(nr+1,k) = a(nr,k)

      end do

c     define gamma and fpond on em grid
c     -  these values used if nz=1

      do k=1,nz
         do i=1,nr+1
            gamma(i,k) = sqrt( 1.+ fpol*abs(a(i,k))**2)
            epond(i,k) = 0.25*(abs(a(i+1,k))**2-abs(a(i-1,k))**2)/2/dr
         end do
c     r=0
         gamma(0,k) = sqrt( 1.+ fpol*abs(a(0,k))**2)
         epond(0,k)=0.
      end do

      call cputime(tim2,v)
      tim = tim + tim2-tim1
      end





