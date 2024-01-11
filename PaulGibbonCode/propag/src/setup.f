c ===========================================
c
c     SETUP
c
c  convert inputs to code norms and setup simulation
c  according to laser, target configuration
c
c =========================================== 
      
      subroutine setup
      include 'dynafoc.h'
      real work1(0:nrmax)

c     laser frequency in 10^15 s**-1
      w0 = 2*pi*c/(xlam*1.e-6)/1.e15

c     laser peak irradiance
      if (ptw.ne.0) then
c     get intensity from I/P power (TW) and focal spot size (microns)
         xi0 = 1.e20*ptw/pi/sigma**2
         a0 = sqrt(xi0*xlam**2/1.38e18)
         
      else if (xi0.ne.0) then
c     I/P intensity
         a0 = sqrt(xi0*xlam**2/1.38e18)

      else
         xi0 = 1.38e18*a0**2
      endif

      if (gas .eq. 'He') then
         nion=2
      else 
         nion=1
      endif


c     i/p timestep normalised to plasma period: convert to laser periods
      
      if (inorm.eq.3) then 
c     i/p in ES plasma units wp**-1, c/wp
         rl = rl/wp0
         zl = zl/wp0
         rplmax = rplmax/wp0
         sigma = sigma/wp0
         tau = tau/wp0
         t0 = t0/wp0
         trun = trun/wp0
c     time in wp**-1
         tfac = wp0
         rfac = wp0
         sfac = 1./a0**2

      else if (inorm.eq.2) then

c     i/p in experimental units:  lengths in microns, times in fs
c     rfac is c/w0; gas density den0 in cm**-3

         wpe = 5.64e4*sqrt(nion*den0)
c    plasma frequency normalised to w0
         wp0 = 1.e-15*wpe/w0
         rfac = xlam/2/pi
         sigma = sigma/rfac
         rl = rl/rfac
         Dlens = Dlens/rfac
         rplmax = rplmax/rfac
         rthom = rthom/rfac
         zplas = zplas/rfac
         zbord = zbord/rfac
         zfoff = zfoff/rfac
         trun = trun/rfac
c     time in fs
         tfac = 1./w0
         tau = tau/tfac
         zl = zl/rfac
         t0 = t0/rfac
         tprobe = tprobe/tfac
c     surface plot norm
         sfac = 1./a0**2
         
      else 
         tfac = 1.
         sfac = 1.
         rfac = 1.
         sfac = 1./a0**2
      endif 

      dr = rl/nr

c     convert pulse length FWHM to 1/e
      tau = tau*0.6

c     em uniform grid spacing
      dr = rl/nr
      dz = zl/nz
c     em advection condition (implicit in plasma rezone)
c     dt = dz
c  timestep for NLSE stability
      dtem = dr**2/4.
      if (nz.eq.1) then
         dt = dtem/2.
      else
         dt=dz
      endif

c     Focussing lens
      if (mfoc.eq.1) then
         Rfoc = Dlens*fno
      endif

c     Determine initial position of plasma edge

      if (mplas.eq.1) then
c     semi-infinite plasma, starting at rear of box
         zedge = zl
         zplen = zplas

      else if (mplas.eq.2 .and. mfoc.eq.1) then
c     focus at centre of plasma
         zedge = Rfoc-zplas/2.
         zplen = zplas       

      else if (mplas.eq.3 .and. mfoc.eq.1) then
c     focus at plasma leading edge
         zedge = Rfoc-zfoff
         zplen = zplas       

      else
         zedge = 0.
         zplen = 2*zplas
      endif
      
      zend = zedge+zplen
      if (rplmax.eq.0) rplmax=rl
      rplmax=amin1(rplmax,rl)

c     em probe starts at RHB
      kprobe = nz



c     beam power (gaussian) at peak of pulse
      pcrit = 16*pi
      power = pi*a0**2*sigma**2
      popc = power/pcrit

c     Rayleigh length (norm)
      Zray = sigma**2
      sigma0 = sigma
c     gas density normalised to nc
      rho0 = wp0*wp0/nion
c  max electron density
      rhoe0 = rho0*nion

c     normalised atomic field
      Ea0 = 0.16*xlam

c     normalised atomic frequency
      wa0 = 22.13*xlam 

c     max ionization rate
      alpha = 4.*wa0*Ea0/a0*exp(-2./3.*Ea0/a0)

      dtioniz=1./alpha
c     initial sub timestep
      dts = dt/icycem


      if (trun.ne.0) then
         nt = trun/dt
         igr = max(1.0,igr/rfac/dt)
         iout = max(1.0,iout/rfac/dt)
         imovie = max(1.0,max(imovie/rfac/dt,nt/80.))
         izsurf = max(izsurf/rfac/dt,nt/300.0)
c         izsurf = max(0.0, nt/60000.0)
c         izsurf = max(0.0, nt/6000.0)
         ihist = max(ihist/rfac/dt,nt/200.0)
      endif    
      end












