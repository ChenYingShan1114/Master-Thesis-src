
c     ==================================================================
c     
c     2D Time-history plots
c     
c     ================================================================== 

      subroutine grzsurf
      include 'dynafoc.h'
      real work1(0:nrmax),rc(0:nrmax),reno(0:nrmax),work2(0:nrmax)
     :     ,work3(0:nrmax) 
     :     ,test1(0:nrmax), test2(0:nrmax), test3(0:nrmax)
      character chars(1)*8,chx*15,chy*15,ctitle*16,cfile*12,cid*1
      character cax*15
      data itc/1/
      save itc

c     offset for focussing geometry
      if (mfoc.eq.1) then
         zoff = Rfoc
      else
         zoff = 0.
      endif

      rlno = rl*rfac
      tnorm = (dt*itime-zoff)*rfac/1000. 
      znorm = (dt*nt-zoff)*rfac

      km=max(t0/dz,1.)
      nrplot = rplmax/dr


      if(mod(nr,isurf).eq.0) then
         npr=nrplot/isurf
      else
         npr=nrplot/isurf+1
      endif

      
c     work1 = laser intensity 
c     work2 = plasma density 
c     work3 = Thomson scatter diagnostic: Int(ne * I)dy

      do i=2*isurf+nrplot,1,-isurf
         ig=(i+isurf-1)/isurf
         work1(ig) = a(i,km)*conjg(a(i,km))/a0**2
         work2(ig) = rhoe(i,km)/rhoe0
         test1(ig) = rhoi(i,km,0)/rho0
         test2(ig) = rhoi(i,km,1)/rho0
         test3(ig) = rhoi(i,km,2)/rho0
         reno(ig) = -r(i)*rfac
      end do

      do i=0,nrplot+2*isurf,isurf
         ig=(i+isurf-1)/isurf
         work1(ig) = a(i,km)*conjg(a(i,km))*sfac
         work2(ig) = rhoe(i,km)/rhoe0
         test1(ig) = rhoi(i,km,0)/rho0
         test2(ig) = rhoi(i,km,1)/rho0
         test3(ig) = rhoi(i,km,2)/rho0
         reno(ig) = r(i)*rfac
      end do

c     Thomson scatter integration: I*ne
c     Absolute intensity I = (a/a0)**2*I0

      do i=0,nrplot+2*isurf,isurf
         ig=(i+isurf-1)/isurf
         work3(ig) = 0.
         xsc = dr*i

         do j=0,nrplot
            ysc = j*dr
            rsc = sqrt(xsc**2 + ysc**2)
            isc = rsc/dr
            work3(ig) = work3(ig) + work1(isc)*work2(isc)*rho0*xi0*dr
         end do
      end do

      
c     data
      write (70,104) (tnorm,-reno(i),work1(i)*xi0
     : ,work2(i)*rhoe0*den0/rho0, test1(i)*rho0*den0/rho0
     : , test2(i)*rho0*den0/rho0
     : , test3(i)*rho0*den0/rho0,i=npr,1,-1)
      write (70,104) (tnorm,reno(i),work1(i)*xi0
     : ,work2(i)*rhoe0*den0/rho0, test1(i)*rho0*den0/rho0
     : , test2(i)*rho0*den0/rho0
     : , test3(i)*rho0*den0/rho0, i=0,npr)
c     write (71,103) (tnorm,reno(i),work2(i),i=0,npr)
c     write (72,103) (tnorm,-reno(i),work3(i),i=npr,1,-1)
c     write (72,103) (tnorm,reno(i),work3(i),i=0,npr)

c     intensity
c     if (itc.eq.1) write (70,104) '! nx ',2*npr+1,' ny',nt/izsurf
c     :,' xmin ',-rlno,' xmax ',rlno,' ymin ',-zoff,' ymax ',znorm

c      write (70,105) (work1(i),i=npr,0,-1)
c      write (70,105) (work1(i),i=1,npr)
c      write (70,*) ' '
c     write (70,103) (work1(i),i=npr,1,-1)
c     write (70,103) (work1(i),i=0,npr)
      
c     electron density
c     if (itc.eq.1) write (71,104) '! nx ',npr+1,' ny',nt/izsurf
c     :,' xmin ',-rlno,' xmax ',0,' ymin ',-zoff,' ymax ',znorm
      write (71,105) (work2(i),i=npr,0,-1)
c      write (71,105) (work2(i),i=0,npr)

      write (71,*) ' '

c     Thomson signal
c     if (itc.eq.1) write (72,104) '! nx ',2*npr+1,' ny',nt/izsurf
c     :,' xmin ',-rlno,' xmax ',rlno,' ymin ',-zoff,' ymax ',znorm
      write (72,105) (work3(i),i=npr,1,-1)
      write (72,105) (work3(i),i=0,npr)
      write (72,*) ' '


c     on-axis line-out: average over 3*dr
c     a2ne(itc) = work3(0)+work3(1)+work3(2)
      itc = itc+1



 103  format (3(1pe12.4))
 104  format (7(1pe12.4))
 105  format (1pe12.4)
      end







