      

c     ======================================================
      
      subroutine ghist
      include 'dynafoc.h'
      real xt(0:ntmax),work1(0:ntmax),work2(0:ntmax),work3(0:ntmax)
      character cax*15

      nhist = nt/ihist

c     time axis -> z if steady-state
      if (inorm.eq.2) then
         hfac = rfac
         cax = '  z (microns)  '
      else
         cax = '     t         '
         hfac = tfac
      endif        

c     offset for focussing geometry
      if (mfoc.eq.1) then
         zoff = Rfoc
      else
         zoff = 0.
      endif

      do i=1,nhist
         xt(i) = (i*dt*ihist-zoff)*hfac
      end do

      if (melec.eq.2) then

c     Lagrangian fluid histories

c     ES energy
         call grxy(xt,Ues,nhist,6100,1,1
     :        ,cax,'   Ues         ','Uese            ')

c     Wave-breaking condition
         call grxy(xt,wbk,nhist,6200,1,1
     :        ,cax,'   wbr flag    ','wbkf            ')

      endif

c     EM pulse histories

c     action
      call grxy(xt,action,nhist,5000,1,1
     :     ,cax,'   Action      ','emac            ')

c     beam radius
      call grxy(xt,rbeam,nhist,5100,1,1
     :     ,cax,'r_{1/e} (mu)   ','rmsb            ')

c     selected power contour
      do i=1,nhist
         work1(i) = rfac*rcon(ipcon,i)
      end do

      call grxy(xt,work1,nhist,5150,1,1
     :     ,cax,' r(P_0/2) (mu) ','rpo2            ')

c     Hamiltonian
      call grxy(xt,H,nhist,5200,1,1
     :     ,cax,'     H         ','hami            ')

c     On-axis intensity in Wcm-2
      call grxy(xt,H2,nhist,5350,1,1
     :     ,cax,'  I0/Wcm^{-2}  ','i0r0            ')

c     On-axis intensity
      call grxy(xt,a2r0,nhist,5300,1,1
     :     ,cax,'    I/I0(r=0)  ','a2r0            ')

c     On-axis phase
      call grxy(xt,aphi,nhist,5600,1,1
     :     ,cax,'    Phi(0,0)   ','aphi            ')

c     electron and gas density on axis
      bden = den0/rho0
      do i=1,nhist
         work1(i) = bden*rhoer0(i)
         work2(i) = bden*rhoar0(i)
      end do

      call grxy(xt,work1,nhist,5400,1,1
     :     ,cax,'   n_e(r=0)    ','rer0            ')
c     :,0.,trun*hfac,dt*100*hfac,0.,1.5,0.5)

      call grxy(xt,work2,nhist,5420,1,1
     :     ,cax,' n_{gas}(r=0)  ','nar0            ')
c     :,0.,trun*hfac,dt*100*hfac,0.,1.5,0.5)

c     ion density on axis
      bden = den0/rho0
      do i=1,nhist
         work1(i) = bden*rho1r0(i)
         work2(i) = bden*rho2r0(i)
         work3(i) = bden*rho0r0(i)
      end do
      call grxy(xt,work3,nhist,5440,1,1
     :     ,cax,'   n_0         ','n0r0            ')
      call grxy(xt,work1,nhist,5450,1,1
     :     ,cax,'   n_i(1+)     ','n1r0            ')

      call grxy(xt,work2,nhist,5460,1,1
     :     ,cax,'   n_i(2+)     ','n2r0            ')


c     Thomson source:  integrated ne(0)*I(0)
      call grxy(xt,a2ne,nhist,5500,1,1
     :     ,cax,' Thomson (w0)  ','a2ne            ')

c     mass cons
      call grxy(xt,edena,nhist,6000,1,1
     :     ,cax,'   Qtot        ','rhet            ')



      end



