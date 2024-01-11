

c =============================================================
 
 

      subroutine init
      include 'dynafoc.h'
      
      wp0 = 0.1
      sigma = 10.
      tau = 10./wp0
      t0 = 5./wp0
      rl = 50.
      zl = 100.
      Zrmass =1./100.
      a0=0.
      den0 = 1.e18
      trun = 0.
      ptw = 0.

      nr = 50
      nz = 10
      nt = 10    
      ihist = 1
      ilas = 1
      wbcm = 1.
      rfac = 1.
      tfac = 1.

c  time label in w0**-1
      mtimlbl = 1
c  adiabatic electrons
      melec = 1
c  self-foc on
      mrel = 1
c  initialise pulse at focus
      mfoc = 0
c  pre-ionized plasma
      mioniz = 0
      mplas = 0
      iout = 1
      imovie=1
      izsurf = 10
      icycem = 1
      gas='H'
      nion = 1
c  em absorber on RH boundary
      mabs = 2
      nabs = 10
      isubm=2
      igrid = 1
      xnumin = 1.e-6
      xnumax = 0.5
      zplen = 10000.
      zend = 10000.
      omemin = 0.
      omemax = 0.
      fpol = 0.5
      ipcon = 1
      rthom = 1.

      timec=0.
      zbord=1.

      read (10,idynf)
      open (11,file='error.log')
      end
