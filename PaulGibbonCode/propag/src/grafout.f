      
      
c     ==================================================================
c     
c     Snapshots
c     
c     ================================================================== 


      subroutine grafout
      include 'dynafoc.h'
      real work1(0:nrmax),rc(0:nrmax),reno(0:nrmax),work2(0:nrmax)
      real zem(0:nzmax),work3(0:nrmax),work4(0:nrmax), work5(0:nrmax)
      complex w
      data iplot/0/
      save iplot

      rlno = rl*rfac
      zlno = zl*tfac
      zmin = 10*int(-t0*tfac/10)
      zmax = 10*int((zmin+zlno)/10)
      km = max(t0/dz,1.)

c     laser amplitude
      do i=0,nr
         work1(i) = abs(a(i,km))**2
         work2(i) = rhoe(i,km)/rhoe0
c     cell centres
         reno(i) = r(i)*rfac
      end do

c     em z-grid
      do i=1,nz
         zem(i) = (i*dz-t0)*tfac
      end do


c     Fixed grid quantities

c     em intensity vs r
      call grxy(reno,work1,nr,500+iplot,0,1
     :     ,'      r        ','     |a|^2     ','emar'//ctime)
c     :,0.,rplmax*rfac,rplmax/5.,0.,2*a0**2,.2*a0**2)

c     em intensity vs z
      do k=1,nz
         work1(k) = abs(a(0,k))**2
      end do
      call grxy(zem,work1,nz,1300+iplot,1,1
     :     ,'      t/fs     ','       |a|^2   ','emaz'//ctime)
c     :     ,zmin,zmax,zlno/5.,0.,1.2*a0**2,.2*a0**2)


c     em phase vs z
      do k=1,nz
         work1(k) = pha(a(0,k))
      end do
      call grxy2(zem,work1,nz,1350+iplot,1,-1
     :     ,'      t/fs     ','     arg(a)    ','emph'//ctime
     :     ,zmin,zmax,zlno/5.,-4.,4.,2.)

      

c     Adiabatic plasma

c     electron density -  vs r at z=z0

      call grxy2(reno,work2,nr,200+iplot,0,-1
     :     ,'      r        ','     ne        ','rher'//ctime
     :     ,0.,rplmax*rfac,rlno/5.,0.,1.5,.5)


c     electron and ion densities -  vs z at r=0
      do k=1,nz
         work1(k) = rhoe(0,k)/rhoe0
         work2(k) =  rhoi(0,k,1)/rho0
         work3(k) =  rhoi(0,k,2)/rho0
         work4(k) = rhoatom(0,k)/rho0
         work5(k) = rhoi(0,k,0)/rho0
      end do
      call grxy2(zem,work4,nz,1000+iplot,1,-1
     :     ,'      t/fs     ','     na        ','natz'//ctime
     :     ,zmin,zmax,zlno/5.,0.,1.5,.5)

      call grxy2(zem,work1,nz,1100+iplot,1,-1
     :     ,'      t/fs     ','     ne        ','rhez'//ctime
     :     ,zmin,zmax,zlno/5.,0.,1.5,.5)

      call grxy2(zem,work2,nz,1050+iplot,1,-1
     :     ,'      t/fs     ','     ni(1+)    ','ni1z'//ctime
     :     ,zmin,zmax,zlno/5.,0.,1.5,.5)

      call grxy2(zem,work3,nz,1150+iplot,1,-1
     :     ,'      t/fs     ','     ni(2+)    ','ni2z'//ctime
     :     ,zmin,zmax,zlno/5.,0.,1.5,.5)
      
      call grxy2(zem,work5,nz,1200+iplot,1,-1
     :     ,'      t/fs     ','     n0        ','ni0z'//ctime
     :     ,zmin,zmax,zlno/5.,0.,1.5,.5)
      

      iplot=iplot + 1

      end















