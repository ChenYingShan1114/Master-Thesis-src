

c     =============================================================
      

      subroutine outini
      include 'dynafoc.h'
      
      drmin=dr
      do i=1,nr
         drmin=amin1(r(i)-r(i-1),drmin)
      end do

      dtmax = drmin**2/4.
c     Debye length
      xld = vte/wp0

c     convert c/w0 to microns
      fmic = xlam/2./pi
c     convert c/w0 to fs
      ffs = 1./w0
c     power
      power = pi*a0**2*(sigma*wp0)**2
      pcrit = 16*pi
      pc = 0.0175/wp0**2
      popc = power/pcrit
      bden = den0
      spot = sigma0*fmic
      tpulse = tau*tfac/0.6
c     local spot size and intensity at Rfoc
      sig = sigma0*sqrt(1.+Rfoc**2/Zray**2)
      af = a0*sigma0/sig
      xIf = xI0*(sigma0/sig)**2
c  cavitation fraction
      sigwp = sigma0*wp0  
      dncav = min(1.,a0**2/2/sigwp**2)

      do istr = 6,15,9
         write (istr,101) tpulse,tpulse/tfac,spot,xI0,xI0*xlam**2
     :, den0, 1.e21/xlam**2


 101     format (/'Laser-plasma parameters'
     :        /'======================='/
     :        /'Pulse length (fs)    ',f12.2,' = ',f12.2,' w0**-1'
     :        /'Spot size (mu)       ',f12.2
     :        /'Intensity  (Wcm-2)   ',1pe12.3
     :        /'Irradiance (Wcm-2mu2)',1pe12.3
     :        /'Background density:  ',1pe12.2
     :        /'Critical density n_c:',1pe12.2)

         write (istr,'(a20,f12.2/a20,f12.2/a20,a12)') 
     : 'Plasma length (mu):  ',zplen*rfac,
     : 'Leading edge (mu):   ',zedge*rfac,
     : 'Target gas:          ',gas
     
        
         write (istr,103) popc*pc,pc,popc,rhoe0,dncav


 103     format (
     :        'Laser power (TW)     ',f12.4
     :        /'Pcrit (TW)           ',f12.2
     :        /'P/Pc:                ',f12.4
     :        /'n_e/n_c:             ',1pe12.2
     :        /'Cavitation Dn_e/n_e: ',1pe12.2)


         write (istr,102) Zray*rfac,fno,fmic*Rfoc,fmic*Dlens
     :        ,sig*fmic,xIf

 102     format (/'Focussing:'
     :        /'=========='/
     :        /'Rayleigh length (mu):',f13.3
     :        /'f-number             ',f13.3
     :        /'Focal length: (mu)   ',f13.3
     :        /'Lens width:   (mu)   ',f13.3
     :        /'Initial spot (mu):   ',f13.2
     :        /'Initial Intensity:   ',1pe12.2
     :        /)

c     Rel. focussing length
         if (popc.gt.1) then
            Zc = Zray*(popc-1.)**(-0.5)
            write (15,105) Zc,Zc*xlam/2/pi
 105        format ('Focussing length: ',f13.3,' = ',f12.2,' microns')
         else
            Zc=0.
         endif

         write (istr,'(/a)') 'Power contours at:'
         write (istr,'(f15.3)')  (pcon(i),i=1,ncon)
c     surface 'movie' plots
         nrplot = rplmax/dr

         if(mod(nr,isurf).eq.0) then
            npr=nrplot/isurf
         else
            npr=nrplot/isurf+1
         endif
         write (istr,'(/a)') 'Mesh for (r,z) plots:'
         write (istr,'(2i6)') nt/izsurf,2*npr+1
         write (istr,'(a20,f12.4,a3,f12.4)') 'Plot region (mu): '
     :,rplmax*rfac,' x ',trun*rfac
         write(istr,'(/a/)') 'Other parameters (code units):'
      end do
      
c  timestep for NLSE stability
      dtem = dr**2/4.
      
      call i0prnt('    nt',nt)
      call i0prnt('    nr',nr)
      call i0prnt('    nz',nz)
      call r0form('    dt',dt,'f13.3')
      call r0form('  dtem',dtem,'f13.3')
c      call r0form('dtioni',dtioniz,'f13.3')
      call i0prnt('  izsurf',izsurf)
      call i0prnt('icycem',icycem)
      call r0form('    rl',rl,'f12.2')
      call r0form('    zl',zl,'f12.2')
      call r0form('    dr',dr,'f14.4')
      call r0form('    dz',dz,'f14.4')
      call r0form('Zrmass',Zrmass,'f14.4')
      call r0form('  wp0 ',wp0,'f13.3')
      call r0form(' a0(f)',a0,'f13.3')
      call r0form('Io(18)',xi0/1.e18,'f13.3')
      call r0form('  tau ',tau,'f12.2')
      call r0form('  t0  ',t0,'f12.2')
      call r0form(' sigma',sigma,'f12.2')
      call r0form('sigma0',sigma0,'f12.2')
      call r0form(' f #  ',fno,'f11.1')
      call r0form(' Rfoc ',Rfoc,'f11.1')
      call r0form('  D   ',Dlens,'f11.1')
      call r0form('  Zray',Zray,'f11.1')
      call r0form('min dr',drmin,'f12.2')
      call r0form('max dt',dtmax,'f12.2')
      call i0prnt(' melec',melec)
      call i0prnt(' igrid',igrid)
      call r0form(' wa0  ',wa0,'f12.2')
      call r0form(' Ea0  ',Ea0,'f12.2')

      end






