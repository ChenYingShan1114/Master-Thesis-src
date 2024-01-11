 


      
c  ================================================================== 
c
c    em pulse histories
c
c  ================================================================== 


      subroutine emhist
      include 'dynafoc.h'
      real work1(0:nrmax),rc(0:nrmax),work2(0:nrmax)
      complex dadr
      data itc/1/
      save itc

      Asum = 0.
      rsum = 0.
      hsum1 = 0.
      hsum2 = 0.
      hsum3 = 0.
      km = max(t0/dz,1.0)
      kmp = max(nz/2,1)
      a(nr+1,km) = (0.,0.)   

      do i=1,nr
        dri = r(i)-r(i-1)
        rmid = 0.5*(r(i)+r(i-1))
        aden = 0.5*(   a(i,km)*conjg( a(i,km)   ) 
     :             + a(i-1,km)*conjg( a(i-1,km) ) )
        gam = sqrt(1.+aden*fpol)
        gam2 = sqrt(1.+a(i,km)*conjg(a(i,km))*fpol)
        gam1 = sqrt(1.+a(i-1,km)*conjg(a(i-1,km))*fpol)
c  Wave action
        Asum = Asum + 2*pi*aden*rmid*dri
c  Beam radius
        rsum = rsum + 2*pi*aden*rmid**3*dri
c  Hamiltonian 
        dadr = (a(i,km)-a(i-1,km))/dri 
        dgdr = (gam2-gam1)/dri
        hsum1 = hsum1 + 0.5*dadr*conjg(dadr)*rmid*dri
        hsum2 = hsum2 - mrel*(1.5-fpol)*wp0**2*(gam-1.)**2*rmid*dri  
        hsum3 = hsum3 - dgdr**2*rmid*dri    
      end do

      action(itc) = Asum
      rbeam(itc) = sqrt(rsum/Asum)*rfac
      H(itc) = 2*pi*(hsum1+hsum2+hsum3)



c  normalised em intensity on axis
      a2r0(itc) = a(0,km)*conjg(a(0,km))*sfac

c  actual intensity on axis in Wcm^-2
      H2(itc) = a2r0(itc)*xi0

c  em phase on axis
      aphi(itc) = pha(a(0,km))

c  critical density
      bden = 1.e21/xlam**2
c  electron density on axis
      rhoer0(itc) = rhoe(0,km)
c  gas density on axis
      rhoar0(itc) = rhoatom(0,km)
c  ion density on axis
      rho0r0(itc) = rhoi(0,km,0)
      rho1r0(itc) = rhoi(0,km,1)
      rho2r0(itc) = rhoi(0,km,2)

c  Thomson scatter integration: I*ne
c  (cf Abel inversion formula)
c  Absolute intensity I = (a/a0)**2*I0

      do i=0,nr
        work1(i) = 0.
        xsc = dr*i
        ny = (rl-xsc)/dr

        do j=0,ny
          rsc = xsc+(j+0.5)*dr
          jsc = rsc/dr
          ysc = rsc/sqrt(rsc**2-xsc**2)
          xi = a(jsc,km)*conjg(a(jsc,km))
          gam2 = 1.+0.5*xi
          work1(i) = work1(i) + xi/gam2*rhoe(jsc,km)*ysc*dr
        end do
      end do

c  Thomson image on axis - average over r=rthom
      iav = rthom/dr
      iav = max(iav,1)
      a2ne(itc)=0.
      do j=0,iav
        a2ne(itc) = a2ne(itc)+work1(j)/iav
      end do

c  Power contours 
      pown = Asum

      do icon=1,ncon
        psum=0.
        i=0
        do while(psum.lt.pcon(icon).and.i.lt.nr)
          i=i+1
          dri = r(i)-r(i-1)
          rmid = 0.5*(r(i)+r(i-1))
          aden = 0.5*(   a(i,km)*conjg( a(i,km)   ) 
     :               + a(i-1,km)*conjg( a(i-1,km) ) )
          psum = psum + 2*pi*rmid*aden*dri/pown
        end do
        rcon(icon,itc) = sqrt( abs( rmid**2 
     : - (psum - pcon(icon))/pi/aden*pown) )
      end do
         

      itc = itc + 1

      end






