c open inpu/output files for unix 
c	open(unit=5,file='alum.d',status='old')
c	open(unit=6,file='oalum.o',status='unknown')
c	open(unit=13,file='falum.d',status='unknown')
c	open(unit=16,file='wch1.o',status='unknown',form='unformatted')
c
c                           med103
c                           =======
c                  a one-dimensional lagrangian
c                  ----------------------------
c                  laser-plasma interaction code
c                  -----------------------------
c
c-------------------------------------------------------------------------
c                    Changes made by Paul Gibbon
c 		     ===========================
c
c       3/8/98          instantaneous intensity and timestep written out
c			to stream 90 (subroutine LASER)
c
c	28/11/95	nt1,nt2 mesh min and max for excited pop. O/P
c			included in input deck.  
c			Defaults: nt1=1, nt2=mesh
c
c-------------------------------------------------------------------------
c
c                      a. djaoui ' s version
c     I have attempted to modify the code such that it is compatible 
c     with the old versions (ie med101) and will run a med101 input data
c     without the need for any changes. The only exception is 
c     the variable SAHA. In order to run the NLTE average atom
c     ionisation model SAHA should be set to 2.
c     For very short pulses the user can now specify  DTMAX in the input 
c     nlhf=t  include high field correction to ib 
c     nltnl=t  include tunnel ionisation in the rates
c
c     this program must be compiled in double precision on the ibm
c     using the autodbl option of the fortvs2 compiler
c
c     x-ray laser post-processor package (radran and subs by sjr/ad1)
c     and neon-like  for si(14) ar(18) ti(22) fe(26)
c     ge(32) kr(36) (ad1)
c     are included and may be turned on via the namelist file.
c     mprint has been tidied up and supports an output file to
c     stream 15 for use with flipper the graphics post-processor
c     (by pr).
c     xrl code has also been modified slightly (amr) to give
c     effective gain as observed by experimentalists.
c     stream 11 has been added to give output to interface
c     with flipper carrying gain data both spatially dependent
c     and effective gain (amr).
c     code has also been altered so that it will compile and
c     run identically on the cray and the ibm. the only loss
c     has been that of the timing subroutines daytim, runtim
c     and jobtim (amr).
c                 index of subroutines
c                 prologue                                          clas
c     labrun      label the run
c     clear       clear variables and arrays
c     preset      set default values
c     data        define   data specific to run
c     auxval      set auxiliary values
c     inital      define physical initial conditions
c     start       prepare to start the calculation
c     source      change standard initial conditions
c     chinic      check initial conditions
c                 calculation                                       clas
c     stepon      step on the calculation
c     motion      controls the motion
c     laser(1)    produces the laser power
c     gie         the terms of the energy equation at level n-1
c     stit        starts the iteration by assuming level 1 = 3
c     statei(1)   describes the state of ions
c     statee(1)   describes the state of electrons
c     coulog      evaluates the coulomb logarithm
c     absorb      evaluates the absorption of laser light
c     brems       the loss rate through bremsstrahlung
c     xchang(1)   the rate of exchange of energy
c     fusion      the energy released by thermonuclear reactions
c     hcduct      thermal conductivities
c     abcd        calculates the terms a, b, c and d
c     findt       the temperature equation (gauss-elimination)
c     timstp      most restrictive value of timestep
c     speed       the hydrodynamic velocities
c     neuman      the artificial viscous pressure
c     moveon      advances the lagrangian coordinates
c     volume      volumes and densities
c     cverge      examines the convergence of ti, te and u
c     energy      kinetic, thermal energy. energy input and gain.
c     shift       shifts the levels 1 - 5 back
c     formp       forms the hydrodynamic pressure
c     boundy(1)   sets the boundary values
c     exam        examines present state of system
c     revers      revers the calculation if break occurs
c     burnup      burn-up of deuterium and tritium
c                 output                                            clas
c     output(1)   control the output
c     mprint(1)   main printer
c     select      select the output
c                 epilogue                                          clas
c     tesend      test for completion of run
c     endrun      terminate the run
c                 diagnostics                                       clas
c     report(3)   reports diagnostic information
c     clist(2)    print common variables
c     arrays(2)   print common arrays
c                 dummy subroutines                                  cla
c     hdcopy(1)   production of hardcopies
c                 index of common blocks
c                 general olympus data
c     combas      basic system parameters
c     comddp      diagnostics and program development
c                 physics
c     comhyd      hydrodynamics
c     comth       thermodynamics
c     comie       ions and electrons
c     comlas      laser variables
c     comfus      thermonuclear reactions
c     comen       energies
c     comcon      physics control
c                 numerical scheme
c     comnc       numerical control parameters
c     comnum      mesh and numerical methods
c                 house-keeping
c     comadm      administrative variables
c                 input-output
c     comout      input-output control variables
c                 index of common variables
c                 physics
c     comhyd      hydrodynamics
c     dm          cell masses                                          r
c     p3          hydrodynamic pressure at level 3                     r
c     pini        initial hydrodynamic pressure                        r
c     r1          coordinates of cell boundaries at level 1            r
c     r3          coordinates of cell boundaries at level 3            r
c     r5          coordinates of cell boundaries at level 5            r
c     rho1        physical density at level 1                          r
c     rho3        physical density at level 3                          r
c     rho5        physical density at level 5                          r
c     rhoini      initial value of rho                                 r
c     rhor        integral of rho * r                                  r
c     rini        initial dimension of system                          r
c     time        the real time in the calculation                     r
c     u2          hydrodynamic velocities at level 2                   r
c     u4          hydrodynamic velocities at level 4                   r
c     uedge       boundary velocity at level 3                         r
c     v1          specific volumes at level 1                          r
c     v3          specific volumes at level 3                          r
c     v5          specific volumes at level 5                          r
c     comth       thermodynamics
c     ddroe1      partial derivative of electron energy                r
c     ddroe3      partial derivative of electron energy                r
c     ddroi1      partial derivative of ion energy                     r
c     ddroi3      partial derivative of ion energy                     r
c     ddte1       partial derivative of electron energy                r
c     ddte3       partial derivative of electron energy                r
c     ddti1       partial derivative of ion energy                     r
c     ddti3       partial derivative of ion energy                     r
c     gammae      ratio of specific heats (electrons)                  r
c     gammai      ratio of specific heats (ions)                       r
c     xappae      electron thermal conductivity                        r
c     xappai      ion thermal conductivity                             r
c     pe1         electron pressure at level 1                         r
c     pe3         electron pressure at level 3                         r
c     pi1         ion pressure at level 1                              r
c     pi3         ion pressure at level 3                              r
c     te1         electron temperature at level 1                      r
c     te3         electron temperature at level 3                      r
c     teini       initial electron temperature                         r
c     ti1         ion temperature at level 1                           r
c     ti3         ion temperature at level 3                           r
c     tiini       initial ion temperature                              r
c     comie       ions and electrons
c     brems1      rate of bremsstrahlung at level 1                    r
c     brems3      rate of bremsstrahlung at level 3                    r
c     degen       degree of electron degeneracy                        r
c     degmax      upper limit of partial degeneracy                    r
c     degmin      lower limit of partial degeneracy                    r
c     effz        average charge number                                r
c     fz1         average charge at level 1
c     fz3         average charge at level 3
c     fzsq1       average <z*z> at level 1
c     fzsq3       average <z*z> at level 3                             r
c     eixch2      rate of ion-electron energy exchange                 r
c     xlc          the coulomb logarithm
c     xmieff       average ion mass number
c     ne          electron number density                              r
c     ni          number density of ions and neutrals                  r
c     omega1      time rate of ion-electron energy exchange            r
c     pmass       proton mass                                          r
c     comlas      laser variables
c     alpha1      absorption coefficient                               r
c     elas1       current laser energy deposited                       r
c     xamda1      wave-length of laser light                           r
c     xaser1      description of laser pulse                           r
c     plas1       laser power as a function of time                    r
c     rabs1       coordinate with critical density                     r
c     rocri1      critical density value for absorption                r
c     xecri1      critical electron density                            r
c     xl1         rate of energy absorbed at level 1                   r
c     xl3         rate of energy absorbed at level 3                   r
c     nabs1       number of cell with critical density                 i
c     comfus      thermonuclear reactions
c     deuter      initial deuterium fraction                           r
c     f1d         fraction of deuterium at level 1                     r
c     f3d         fraction of deuterium at level 3                     r
c     f1h         fraction of hydrogen at level 1                      r
c     f3h         fraction of hydrogen at level 3                      r
c     f1he3       fraction of helium3 at level 1                       r
c     f3he3       fraction of helium3 at level 3                       r
c     f1he4       fraction of helium4 at level 1                       r
c     f3he4       fraction of helium4 at level 3                       r
c     f1neu       fraction of neutrons released at level 1             r
c     f3neu       fraction of neutrons released at level 3             r
c     f1ntrl      fraction of neutrals at level 1                      r
c     f3ntrl      fraction of neutrals at level 3                      r
c     f1t         fraction of tritium at level 1                       r
c     f3t         fraction of tritium at level 3                       r
c     f1x         fraction of an extra element at level 1              r
c     f3x         fraction of an extra element at level 3              r
c     heliu3      initial fraction of helium3                          r
c     heliu4      initial fraction of helium4                          r
c     hydrog      initial fraction of hydrogen                         r
c     xetral      initial fraction of neutrals                         r
c     xtrlms      mass number of neutral element                       r
c     pneut1      neutron power at level 1                             r
c     pneut3      neutron power at level 3                             r
c     r1dd        number of d-d reactions at level 1                   r
c     r3dd        number of d-d reactions at level 3                   r
c     r1dhe3      number of d-he3 reactions at level 1                 r
c     r3dhe3      number of d-he3 reactions at level 3                 r
c     r1dt        number of d-t reactions at level 1                   r
c     r3dt        number of d-t reactions at level 3                   r
c     rneut1      neutron flux at level 1                              r
c     rneut3      neutron flux at level 3                              r
c     tinucl      lower threshold for thermonuclear reactions          r
c     totneu      total number of neutrons released                    r
c     tritiu      initial fraction of tritium                          r
c     xmass       mass number of extra element                         r
c     xtra        initial fraction of extra element                    r
c     xz          charge number of extra element                       r
c     ye1         fusion energy apportioned to electrons               r
c     ye3         fusion energy apportioned to electrons               r
c     yi1         fusion energy apportioned to ions                    r
c     yi3         fusion energy apportioned to ions                    r
c     yield       ratio of thermonuclear/laser energies                r
c     comen       energies
c     eefuse      thermonuclear energy given to electrons              r
c     eerror      error in the energy calculation                      r
c     eifuse      thermonuclear energy given to ions                   r
c     eindt1      energy input at level 1                              r
c     eindt3      energy input at level 3                              r
c     einput      total energy input                                   r
c     eloss       total energy loss through ff radiation               r
c     en          total energy of system                               r
c     eneutr      total energy of escaping neutrons                    r
c     pv          total thermal energy of system                       r
c     usqm        total kinetic energy of system                       r
c     comcon      physics control
c     pmin        permissible minimum of pressure                      r
c     rhomin      permissible minimum of density                       r
c     temin       permissible minimum of electron temperature          r
c     timin       permissible minimum of ion temperature               r
c     umin        smallest velocity allowed for safety                 r
c     mstep       timestep counter                                     i
c     ncase       selects boundary conditions                          i
c     ngeom       selects the geometry                                 i
c     mlbrms      program switch for bremsstrahlung                    l
c     mlecon      program switch for electron heat conduction          l
c     mlfuse      program switch for thermonuclear reactions           l
c     mlicon      program switch for ion heat conduction               l
c     mlx         program switch for i-e exchange                      l
c     nlabs       switch : absorption by inverse bremsstrahl.          l
c     nlbrms      switch : bremsstrahlung                              l
c     nlburn      switch : burn-up of d and t                          l
c     nlcri1      switch : dump laser power at rocri1                  l
c     nldepo      switch : deposit thermonuclear energies              l
c     nlecon      switch : electron thermal conduction                 l
c     nlfuse      switch : thermonuclear reactions                     l
c     nlicon      switch : ion thermal conduction                      l
c     nlmove      switch : hydrodynamic motion                         l
c     nlpfe       switch : perfect electron gas laws                   l
c     nlpfi       switch : perfect ion gas laws                        l
c     nlx         switch : ion-electron energy exchange                l
c                 numerical scheme
c     comnc       numerical control parameters
c     ak0         control the time centering                           r
c     ak1         control dt by soundspeed                             r
c     ak2         control dt by variations in volume                   r
c     ak3         control dt by variations in ti                       r
c     ak4         control dt by variations in te                       r
c     ak5         spare dt control parameter                           r
c     ak6         spare dt control parameter                           r
c     ak7         control rate of lase energy input                   r
c     ak8         controls rate of change of energies in each zone AA
c     ak9         control rate of change of z*  (max allowed ak9)AA    r
c     bneum       determines the rate of shock heating                 r
c     deltat      user specified value of dt                           r
c     dt2         timestep (t3-t1)                                     r
c     dt3         timestep (t4-t2)                                     r
c     dt4         timestep (t5-t3)                                     r
c     dtemax      controls convergence of te                           r
c     dtimax      controls convergence of ti                           r
c     dtmax       maximum value of dt                                  r
c     dtmin       minimum value of dt                                  r
c     dumax       controls convergence of u                            r
c     rdt2        = 1.0 / dt2                                          r
c     rdt3        = 1.0 / dt3                                          r
c     rdt4        = 1.0 / dt4                                          r
c     nceldt      the cell number that sets the limit of dt            i
c     ncondt      condition determining dt                             i
c     nit         current number of iterations                         i
c     nitmax      maximum number of iterations allowed                 i
c     break       break the calculation                                l
c     nlgoon      iterations converge                                  l
c     nlite       switch : solve equations by iterations               l
c     comnum      mesh and numerical methods
c     ae          coefficient in electron energy equation              r
c     ai          coefficient in ion energy equation                   r
c     be          coefficient in electron energy equation              r
c     bi          coefficient in ion energy equation                   r
c     ce          coefficient in electron energy equation              r
c     ci          coefficient in ion energy equation                   r
c     de          coefficient in electron energy equation              r
c     di          coefficient in ion energy equation                   r
c     e           auxiliary array used for gauss elimination           r
c     f           auxiliary array used for gauss elimination           r
c     ge          coefficient in electron energy equation              r
c     gi          coefficient in ion energy equation                   r
c     q2          viscous pressure at level 2                          r
c     q4          viscous pressure at level 4                          r
c     teite       level 3 value of te from previous iteration          r
c     tiite       level 3 value of ti from previous iteration          r
c     uite        level 4 value of u from previous iteration           r
c     mesh        number of cells in mesh                              i
c     nj          auxiliary variable (=njp1-1=mesh)                    i
c     njm1        auxiliary variable (=nj-1)                           i
c     njp1        number of cell boundaries (=mesh+1)                  i
c     nl          number of cells in mesh (=mesh)                      i
c     nlm1        auxiliary variable (=nl-1)                           i
c     nlp1        auxiliary variable (=nl+1)                           i
c                 house-keeping
c     comadm      administrative variables
c     piq         general purpose array                                r
c     tstop       maximum value of time permitted                      r
c     maxdim      maximum permissible mesh size                        i
c     maxrun      maximum number of timesteps permitted                i
c     ndump       frequency of dumping common areas                    i
c     nrep        repetition of film frames                            i
c     nldump      switch : dumping of common areas                     l
c     nlemp       switch : emergency printing                          l
c                 input-output
c     comout      input-output control variables
c     buf1        output buffer (r)                                    r
c     buf2        output buffer (u)                                    r
c     buf3        output buffer (density)                              r
c     buf4        output buffer (pressure)                             r
c     buf5        output buffer (ion temperature)                      r
c     buf6        output buffer (electron temperature)                 r
c     scp         input-output scaling factor for pressure             r
c     scr         input-output scaling factor for coordinates          r
c     scrho       input-output scaling factor for density              r
c     scte        input-output scaling factor for te                   r
c     scti        input-output scaling factor for ti                   r
c     sctime      input-output scaling factor for time                 r
c     nfilm       frequency of film frame production                   i
c     nhdcpy      frequency of hardcopy production                     i
c     np1         array printing selector : start                      i
c     np2         array printing selector : end                        i
c     np3         array printing selector : increment                  i
c     nprnt       frequency of printing                                i
c     nlfilm      switch : film production                             l
c     nlhcpy      switch : hardcopy production                         l
c     nlprnt      switch : printing                                    l
c     ztnext      time for next printout
c     -----------------------------------------------------------
c     it is possible
c     to  alter the basic performance of the medusa code. for
c     all 4 cases of boundary conditions it is possible  to  alter
c     pressure, temperature, velocity etc. at the boundary cell in
c     subprogram "boundy" (2.25).  also the  shape  of  the  laser
c     pulse  for  calculations on case 1 can be changed in subpro-
c     gram  "laser"  (2.3)  provided  the  appropriate  quantities
c     (laser  power,  energy,  pulse  length,  switch on-off times
c     etc.) are stored in array "xaser1" as described in (2.3) and
c     also  in  the  list  presented below. use can be made of the
c     all-purpose buffer "piq" as long as those elements of  "piq"
c     which  appear  in the list given below are not affected.  if
c     particular modifications to  the  computational  scheme  are
c     required  these  should  be  coded up in subprogram "expert"
c     (0.4) as explained in the "computer physics  communications"
c     paper  which  describes  the olympus system.  finally, it is
c     possible to change many of the coefficients which  determine
c     the physical processes involved by setting selected elements
c     of array "piq" to values different from  zero  as  explained
c     below.   in  the  list below we present all variables of the
c     namelist "newrun" which is  read  in  by  subprogram  "data"
c     (1.4).  the  list  contains  examples of the range of values
c     possible and also refers to the default value of a variable.
c     'newrun'  which is read in by subprogram  'data' (1.4). the list
c     contains examples of the range of values possible and also refers
c     to the default value of a variable.
c
c******** c o n t e n t s  o f  n e w r u n  *********
c
cname  value       purpose of variable                     default value
c
cak0=50.,          control the time centering                  50.0
cak1=0.25,          control dt by soundspeed                    0.25
cak2=0.25,          control dt by variations in volume          0.25
cak3=0.25,         control dt by variations in ti              0.25
cak4=0.25,         control dt by variations in te              0.25
cak5=0.0,          spare dt control parameter                  0.0
cak6=0.0,          spare dt control parameter                  0.0
cak7=0.0,          spare dt control parameter                  0.0
cak8=0.0,          control rate of input energy                0.0
cak9=0.0,          control rate of ionisation energy           0.0
cbneum=1.2,        determines the rate of shock heating        1.0
cdeltat=1.0e-20,   user specified value of dt                  0.0
cdeuter=0.5,       initial fraction of deuterium               0.5
cdtemax=0.2,       controls convergence of te                  0.1
cdtimax=0.2,       controls convergence of ti                  0.1
cdumax=0.2,        controls convergence of u                   0.1
cgammae=1.667,     ratio of specific heats (electrons)         1.667
cgammai=1.667,     ratio of specific heats (ions)              1.667
cheliu3=0.0,       initial fraction of helium3                 0.0
cheliu4=0.0,       initial fraction of helium4                 0.0
chydrog=0.0,       initial fraction of hydrogen                0.0
cxamda1=1.0e-6,    wavelength of laser light                   1.0e-5
cnetral=0.0,       initial fraction of neutrals                0.0
cntrlms=0.0,       atomic mass number of neutrals              0.0
crhoini=100.0,     initial constant density                    124.0
crini=0.001,       initial dimension of system                 0.00048
cscp=1.0,          scale input-output of pressure              1.0
cscr=100.0,        scale input-output of coordinates           1.0
cscrho=1000.0,     scale input-output of density               1.0
cscte=1.0,         scale input-output of electron temperature  1.0
cscti=1.0,         scale input-output of ion temperatureature  1.0
csctime=1.0e6,     scale input-output of time                  1.0
cteine=5.0e.,       initial electron temperature                5.0e3
ctiini=5.0e3,      i5itial ion temperature                     5.0e3
ctinucl=1.0e7,     ignition temperature (fusion)               1.0e7
ctritiu=0.5,       initial fraction of tritium                 0.5
ctstop=0.1,        maximum value of time permitted             1.0e-6
cxmass=0.0,        mass number of extra element                0.0
cxtra=0.0,         initial fraction of extra element           0.0
cxz=0.0,           charge number of extra element              0.0
cmesh=40,          number of cells in mesh                     40
cncase=1,          boundary conditions : laser                 1
cncase=2,          boundary conditions : applied wall temperature
cncase=3,          boundary conditions : applied wall pressure
cncase=4,          boundary conditions : applied wall velocity
cndump=1,          frequency of dumping common areas           10
cnfilm=1,          frequency of film frame production          1
cngeom=1,          plane geometry
cngeom=2,          cylindrical geometry
cngeom=3,          spherical geometry                          3
cnhdcpy=500,       frequency of hardcopy production            100
cnitmax=4,         maximum number of iterations allowed        5
cnp1=1,            array printing selector : start             1
cnp2=40,           array printing selector : end               40
cnp3=2,            array printing selector : increment         2
cnprnt=100,        frequency of printing                       100
cnrep=1,           repetition of film frames                   0
cnlabs=f,          switch : absorption by inverse bremsstrahl. t
cnlbrms=f,         switch : bremsstrahlung                     t
cnlburn=t,         switch : burn-up of d and t                 t
cnlcri1=f,         switch : dump laser power at rocri1         t
cnldepo=t,         switch : deposit thermonuclear energies     t
cnldump=f,         switch : dumping of common areas            f
cnlecon=t,         switch : electron thermal conduction        t
cnlemp=t,          switch : emergency printing                 t
cnlfilm=t,         switch : film production                    t
cnlfuse=t,         switch : thermonuclear reactions            t
cnlhcpy=t,         switch : hardcopy production                t
cnlicon=f,         switch : ion thermal conduction             t
cnlite=t,          switch : solve equations by iterations      t
cnlmove=t,         switch : hydrodynamic motion                t
cnlpfe=f,          switch : perfect electron gas laws          t
cnlpfi=t,          switch : perfect ion gas laws               t
cnlprnt=t,         switch : printing                           t
cnlx=f,            switch : ion-electron energy exchange       t
c
cxaser1(1)=0.0,    time to switch on the laser                 0.0
cxaser1(2)=1.0e9,  initial power level                         0.0
cxaser1(3)=2.0e-8, pulse duration                              0.0
cxaser1(4)=2.0,    pulse shape factor                          0.0
cxaser1(5)=2.0e-8, time to switch off the laser                0.0
cxaser1(6)=1.0e14, maximum permissible power level             0.0
cxaser1(7)=2.0e4,  maximum energy available                    0.0
cc
cpiq(3)=1.0,       determines boundary value of pressure or te 0.0
cpiq(4)=1.25,      determines boundary value of pressure       0.0
cpiq(5)=1.0,       determines boundary value of te             0.0
cpiq(6)=1.0e4,     determines boundary value of pressure       0.0
cpiq(8)=1.0,       modify boundary pressure in case 1          0.0
cpiq(11)=1.0,      change initial coordinates                  0.0
cpiq(12)=1.0,      change initial density distribution         0.0
cpiq(13)=1.0,      change initial composition of pellet        0.0
cpiq(20)=0.1,      modify coulomb logarithm coefficient        0.0
cpiq(21)=5.0,      modify absorption coefficient               0.0
cpiq(24)=2.0,      modify bremsstrahlung coefficient           0.0
cpiq(25)=1.0,      modify exchange rate coefficient            0.0
cpiq(26)=-0.999,   modify ion conduction coefficient           0.0
cpiq(27)=-0.5,     modify electron conduction coefficient      0.0
cpiq(31)=1.0,      print all of buf1                           0.0
cpiq(32)=-1.0,     no printing of buf2                         0.0
cpiq(33)=0.0,      normal printing of buf3 from np1 to np2/np3 0.0
cpiq(34)=-1.0,     no printing of buf4                         0.0
cpiq(35)=0.0,      normal printing of buf5 from np1 to np2/np3 0.0
cpiq(36)=3.0,      print all of buf6                           0.0
c
c                   c1.1     basic system parameters
cnrun=500,         number of timesteps to be run            999999
c                  additional pulse shapes
c
c up to five gaussion pulses if gauss=1.0,   anpuls=1.,2.,3.,4.,5.,
c 1st : xaser1(6)=pmax :    xaser1(4)plenth, :  xaser1(2)=pmult,
c     :fraction of pmax:    plenth for nth   :  pmult for nth pulse
c 2nd : piq(38)= 1.00, :    piq(42)=6.0e-12, :  piq(46)=168.334,
c 3rd : piq(39)= 0.0,  :    piq(43)=0.00000, :  piq(47)=0.0,
c 4th : piq(40)= 0.0,  :    piq(44)=0.00000, :  piq(48)=0.0,
c 5th : piq(41)= 0.0,  :    piq(45)=0.00000, :  piq(49)=0.0,
c
c       triangular pulse
c piq(37) =gauss < 0 for triangular pulse
c ad1 gauss = -1 (triangular) -2 (triangular with flat top)
c triangular pulse with prepulse and flat top can be set as
c prepulse power from 0 at time 0 to xaser1(21) w/m**2 at time xaser1(20
c linear incresase of power to pmax=xaser1(6) for plenth=xaser1(4)
c flat top at xaser1(6) for xaser1(22)
c then linear decrease to zero at toff=xaser1(5)
c        gauss must be = -2.
c        time dependent excitation ionisation + additional variables
c
c nlres=.true. restart case from previous run
c saha=2.0, time dependent excitation and ionisation
c bbtrap trapping factor for bb radiation 0 no trap 1 full trap
c fbtrap trapping factor for bound-free radiation
c
c idrflg=1, switch dielectronic recombination on
c drmul =   multiplying factor for dielectronic recombination
c fwide =   width of focus in x-ray laser experiment
c
c                  x-ray laser calculations
c     icxrl sets whether xrl calculations are carried out
c     icxrl=1 true       icxrl=0 false
c choose ionisation stage
c istage = 1   h-like
c istage = 2   li-lik
c istage = 3   na-like
c istage = 4   ne-like
c istage = 5   ni-like
c set optical depth flag
c idflg=0  - no optical depth correction
c idflg=1  - include optical depth correction
c ropmul is the factor multiplying the optical depth
c     ropmul=0.0
c set thermal band flag
c itbflg=0 - no thermal band
c itbflg=1 - force population of highest level to be in lte
c     itbflg=1
c set motional stark broadening flag
c istflg=0 - do not include motional stark broadening
c istflg=1 - include motional stark broadening on the balmer alpha line
c     istflg=0
c set photopumping flag and values
c ipuflg=0 - do not include photopumping
c ipuflg=/0 - include photopumping
c tpon - time pump on
c tpoff - time pump off
c nlp - lower level for the pump
c nup - upper level for the pump
c     ipuflg=1
c     tpon=2.15e-10
c     tpoff=2.85e-10
c     nlp=1
c     nup=3
c     rmpd=0.002
c setup wavelength of the lasing line(in cm)
c     dlamda=81.0e-08
c setup lengths used in the axial line intensity calculation
c     nlmax=2
c     fl(1)=0.0001
c     fl(2)=0.55
c tape15 is used for dumping the simulation state in order to restart
c tape16 is used for readind restart information when nlres=.true.
c 0.0 fortran main program
c c1.1 basic system parameters
c end of basic system parameters
c ad1 * suppress error messages
c     overflow  (gt 16**63 = 7.2e75)
c      CALL errset(207,256,-1,-1)
c     underflow (lt 16**-65 = 5.4e-79)
c      CALL errset(208,256,-1,-1)
c time allocated to job
c     call jobtim(altime)
c set up the basic control data
      CALL basic
c print date and time
c     call daytim
c control the run
      CALL cotrol
      STOP 1000
      END


      SUBROUTINE abcd
      PARAMETER(kk=1001,underf=1.E-30)
c 2.14 calculates the terms a, b, c and d
c c1.9 development and diagnostic parameters
      DIMENSION nadump(20),npdump(20),nvdump(20)
      LOGICAL nlched,nlhead(9),nlomt1(50),nlomt2(50),nlomt3(50),nlrept
      COMMON /comddp/ maxdum,mxdump,nadump,nclass,npdump,npoint,nsub,
     &                nvdump,nlched,nlhead,nlomt1,nlomt2,nlomt3,nlrept
c end of development and diagnostic parameters
      LOGICAL lphi(kk)
      DIMENSION deltae(50),ebdy(51),emid(50),ha(kk),hb(kk),hc(kk),hd(kk)
     &          ,hdiff(kk),he(kk),heta(kk),hf(kk),hg(kk),hgrp(kk),hj(kk)
     &          ,hlas(kk),hloss(50),hncold(kk),hp(kk),hphi(kk),htaug(kk)
     &          ,htherm(kk),htr(kk),xhot1(kk,50),xhot3(kk,50)
      COMMON /comhot/ igrp,ngroup,lphi,deltae,ebdy,ehot,emid,ha,hb,hc,
     &                hd,hdiff,he,heta,hf,hg,hgrp,hj,hlas,hloss,hncold,
     &                hp,hphi,htaug,htherm,htr,xhot1,xhot3,thot
c c2.1 hydrodynamics
      DIMENSION dm(kk),p3(kk),r1(kk),r3(kk),r5(kk),rho1(kk),rho3(kk),
     &          rho5(kk),u2(kk),u4(kk),v1(kk),v3(kk),v5(kk),lstat(kk)
      COMMON /comhyd/ dm,p3,pini,r1,r3,r5,rho1,rho3,rho5,rhoini,rhor,
     &                rini,time,u2,u4,uedge,v1,v3,v5
c end of hydrodynamics
c c2.2 thermodynamics
      DIMENSION ddroe1(kk),ddroe3(kk),ddroi1(kk),ddroi3(kk),ddte1(kk),
     &          ddte3(kk),ddti1(kk),ddti3(kk),xappae(kk),xappai(kk),
     &          pe1(kk),pe3(kk),pi1(kk),pi3(kk),te1(kk),te3(kk),ti1(kk),
     &          ti3(kk)
      COMMON /comth / ddroe1,ddroe3,ddroi1,ddroi3,ddte1,ddte3,ddti1,
     &                ddti3,gammae,gammai,xappae,xappai,pe1,pe3,pi1,pi3,
     &                te1,te3,teini,ti1,ti3,tiini
c end of thermodynamics
c c2.3 ions and electrons
      DIMENSION brems1(kk),brems3(kk),ddz1(kk),ddz3(kk),degen(kk),
     &          effz(kk),eixch2(kk),fidash(kk),fz1(kk),fz3(kk),fzsq1(kk)
     &          ,fzsq3(kk),xlc(kk),xmieff(kk),xne(kk),xni(kk),omega1(kk)
      COMMON /comie / brems1,brems3,ddz1,ddz3,degen,degmax,degmin,effz,
     &                eixch2,fidash,fz1,fz3,fzsq1,fzsq3,xlc,xmieff,xne,
     &                xni,omega1,pmass
c end of ions and electrons
c c2.4 laser variables
      DIMENSION alpha1(kk),xaser1(kk),xl1(kk),xl3(kk)
      COMMON /comlas/ alpha1,elas1,xamda1,xaser1,xecri1,plas1,rabs1,
     &                rocri1,xl1,xl3,nabs1
c end of laser variables
c c2.5 thermonuclear reactions
      DIMENSION f1d(kk),f1h(kk),f1he3(kk),f1he4(kk),f1neu(kk),f1ntrl(kk)
     &          ,f1t(kk),f1x(kk),f1x1(kk),f1x2(kk),f1x3(kk),f1x4(kk),
     &          f3d(kk),f3h(kk),f3he3(kk),f3he4(kk),f3neu(kk),f3ntrl(kk)
     &          ,f3t(kk),f3x(kk),f3x1(kk),f3x2(kk),f3x3(kk),f3x4(kk),
     &          r1dd(kk),r1dhe3(kk),r1dt(kk),r3dd(kk),r3dhe3(kk),
     &          r3dt(kk),ye1(kk),ye3(kk),yi1(kk),yi3(kk)
      COMMON /comfus/ deuter,f1d,f1h,f1he3,f1he4,f1neu,f1ntrl,f1t,f1x,
     &                f1x1,f1x2,f1x3,f1x4,f3d,f3h,f3he3,f3he4,f3neu,
     &                f3ntrl,f3t,f3x,f3x1,f3x2,f3x3,f3x4,heliu3,heliu4,
     &                hydrog,xetral,xtrlms,pneut1,pneut3,r1dd,r1dhe3,
     &                r1dt,r3dd,r3dhe3,r3dt,rneut1,rneut3,tinucl,totneu,
     &                tritiu,xmass,xmass1,xmass2,xmass3,xmass4,xtra,
     &                xtra1,xtra2,xtra3,xtra4,xz,xz1,xz2,xz3,xz4,ye1,
     &                ye3,yi1,yi3,yield
c end of thermonuclear reactions
c c2.7 physics control
      LOGICAL mlbrms,mlecon,mlfuse,mlicon,mlx,nlabs,nlbrms,nlburn,
     &        nlcri1,nldepo,nlecon,nlfuse,nlicon,nlmove,nlpfe,nlpfi,nlx,
     &        nlhf,nltnl
      COMMON /comcon/ mstep,ncase,ngeom,mlbrms,mlecon,mlfuse,mlicon,mlx,
     &                nlabs,nlbrms,nlburn,nlcri1,nldepo,nlecon,nlfuse,
     &                nlicon,nlmove,nlpfe,nlpfi,nlx,nlhf,nltnl,pmin,
     &                rhomin,temin,timin,umin
c end of physics control
c c3.2 mesh and numerical methods
      DIMENSION ae(kk),ai(kk),be(kk),bi(kk),ce(kk),ci(kk),de(kk),di(kk),
     &          e(kk),f(kk),ge(kk),gi(kk),q2(kk),q4(kk),teite(kk),
     &          tiite(kk),uite(kk)
      COMMON /comnum/ mesh,nj,njm1,njp1,nl,nlm1,nlp1,ae,ai,be,bi,ce,ci,
     &                de,di,e,f,ge,gi,q2,q4,teite,tiite,uite
c end of mesh and numerical methods
c c3.1 numerical control parameters
      LOGICAL break,nlgoon,nlite
      COMMON /comnc / nceldt,ncondt,nit,nitmax,break,nlgoon,nlite,ak0,
     &                ak1,ak2,ak3,ak4,ak5,ak6,ak7,ak8,ak9,bneum,deltat,
     &                dt2,dt3,dt4,dtemax,dtimax,dtmax,dtmin,dumax,rdt2,
     &                rdt3,rdt4
c end of numerical control parameters
c end of basic system parameters
c c4.1 - administrative variables
      LOGICAL nldump,nlemp
      DIMENSION piq(302)
      COMMON /comadm/ maxdim,maxrun,ndump,nrep,nldump,nlemp,piq,tstop
c c4.1 - end administrative variables
      COMMON /dpi   / ecbb(kk),ecbf(kk),eraddr(kk),erbf(kk),eati(kk),
     &                dcvz(kk),efbx2(kk),efrx2(kk),erbb(kk),erfb(kk),
     &                acfc(10,kk),acfr(10,kk),bcfc(10,kk),bcfr(10,kk),
     &                edfb(kk),ter(kk),dqdtim(kk),nmaxr
      DATA iclass,isub/2,14/
      IF(nlomt2(isub))RETURN
c ad1 debug output
c     write(6,*)mstep,xl3(nl),pe3(nl),te3(nl),u4(nlp1)
c 1. coefficients for perfect gases
      igeom=ngeom-1
c zkm is the ratio  k/m(proton)
      zkm=8255.17
c the terms at l=nl are included in the loop, but will be modified
      DO 100 l=2,nl
       j=l
c treat factors common to ions and electrons
       zpi=zkm/xmieff(l)
       zpe=zkm*fz3(l)/xmieff(l)
       zmra=(r3(j)**igeom)*dt2/(dm(l)*(r3(j+1)-r3(j-1)))
       zmrc=(r3(j+1)**igeom)*dt2/(dm(l)*(r3(j+2)-r3(j)))
       zdrho=rho3(l)-rho1(l)
       zdv=v3(l)-v1(l)
       zpi3dv=0.5*zdv*zpi*rho3(l)
       zpe3dv=0.5*zdv*zpe*rho3(l)
       zdz=fz3(l)-fz1(l)
c evaluate the coefficients of equation 62
c eq.63
       ai(l)=-zmra*xappai(j)
       ae(l)=-zmra*xappae(j)
c eq.64
       ci(l)=-zmrc*xappai(j+1)
       ce(l)=-zmrc*xappae(j+1)
c eq.65 transfer pressure term p*dv/dt to left hand side of eq.62
       bi(l)=0.5*(ddti1(l)+ddti3(l))-ai(l)-ci(l)+zpi3dv
       be(l)=0.5*(ddte1(l)+ddte3(l)+dcvz(l))-ae(l)-ce(l)+zpe3dv
c eq.66
       di(l)=0.5*dt2*(-eixch2(l)+yi3(l))
     &       +0.5*(ti1(l)*(ddti1(l)+ddti3(l))-zdv*pi1(l))
     &       +0.5*zdrho*(ddroi1(l)+ddroi3(l))
       de(l)=0.5*dt2*(eixch2(l)+ye3(l)+xl3(l)+brems3(l))
     &       +0.5*(te1(l)*(ddte1(l)+ddte3(l)-dcvz(l))-zdv*pe1(l))
     &       +0.5*zdrho*(ddroe1(l)+ddroe3(l)) + htherm(l)
     &       +dt2*(efbx2(l)+efrx2(l))
       if(abs(piq(58)-1.).lt.underf) then
         de(l)=de(l)  - 0.5 * zdz * (ddz1(l) + ddz3(l))
       endif        
100   CONTINUE
c     2. set coefficients at boundary points
c     first define quantities common to all cases
      zpi=zkm/xmieff(1)
      zpe=zkm*fz3(1)/xmieff(1)
c     ad1
      zmrc=(r3(2)**igeom)*dt2/(dm(1)*(r3(3)-r3(1)))
      zdrho=rho3(1)-rho1(1)
      zdv=v3(1)-v1(1)
      zdz=fz3(1)-fz1(1)
      zpi3dv=0.5*zdv*zpi*rho3(1)
      zpe3dv=0.5*zdv*zpe*rho3(1)
c change the terms at the boundary (prescription after eq.67)
      ai(1)=0.0
      ae(1)=0.0
      ci(1)=-zmrc*xappai(2)
      ce(1)=-zmrc*xappae(2)
      bi(1)=0.5*(ddti1(1)+ddti3(1))-ci(1)+zpi3dv
      be(1)=0.5*(ddte1(1)+ddte3(1)+dcvz(1))-ce(1)+zpe3dv
      di(1)=0.5*dt2*(-eixch2(1)+yi3(1))
     &      +0.5*(ti1(1)*(ddti1(1)+ddti3(1))-zdv*pi1(1))
     &      +0.5*zdrho*(ddroi1(1)+ddroi3(1))
      de(1)=0.5*dt2*(eixch2(1)+ye3(1)+xl3(1)+brems3(1))
     &      +0.5*(te1(1)*(ddte1(1)+ddte3(1)-dcvz(1))-zdv*pe1(1))
     &      +0.5*zdrho*(ddroe1(1)+ddroe3(1))+htherm(1)
     &      +dt2*(efbx2(1)+efrx2(1))
       if(abs(piq(58)-1.).lt.underf) then
         de(1)=de(1)  - 0.5 * zdz * (ddz1(1) + ddz3(1))
       endif
c     xxxxxxxxxxxxxxxxxxxxxx
      IF(ncase.eq.1)THEN
      ELSEIF(ncase.eq.2)THEN
c 2.2 te or ti determined at the boundary
c to work out the incoming heat flux we assume the point nj+2 to be
c placed at the same distance from nj+1 as the point nj is.
       zmrc=(r3(njp1)**igeom)*dt2/(dm(nl)*(r3(njp1)-r3(nj))*2.0)
       zpe=zkm*fz3(nl)/xmieff(nl)
       zpi=zkm/xmieff(nl)
       zdv=v3(nl)-v1(nl)
       zpi3dv=0.5*zdv*zpi*rho3(nl)
       zpe3dv=0.5*zdv*zpe*rho3(nl)
c        the conductivities at the wall point
       ci(nl)=-zmrc*xappai(njp1)
       ce(nl)=-zmrc*xappae(njp1)
c        now use the value of c to determine b
       bi(nl)=0.5*(ddti1(nl)+ddti3(nl))-ai(nl)-ci(nl)+zpi3dv
       be(nl)=0.5*(ddte1(nl)+ddte3(nl)+dcvz(nl))-ae(nl)-ce(nl)+zpe3dv
c        subtract the term c*t from the right hand side of eq.62
       di(nl)=di(nl)-ti3(nlp1)*ci(nl)
       de(nl)=de(nl)-te3(nlp1)*ce(nl)
c        set the c-coefficients to zero after the subtraction
       ci(nl)=0.0
       ce(nl)=0.0
c        applied pressure
c        in this case the heat flux entering from the outside is zero.
c        this case is therefore exactly like case 1
c        2.4 boundary velocity determined
c        in case 4 no heat enters the system, thus
       GOTO 200
      ELSEIF(ncase.eq.3)THEN
      ELSEIF(ncase.ne.4)THEN
       CALL gotoer(' MEDUSA: abcd: goto error, ncase is',ncase)
      ENDIF
c     xxxxxxxxxxxxxxxxxx
c 2.1 thermally insulated boundary (laser case)
      zpi=zkm/xmieff(nl)
      zpe=zkm*fz3(nl)/xmieff(nl)
      zdv=v3(nl)-v1(nl)
      zpi3dv=0.5*zdv*zpi*rho3(nl)
      zpe3dv=0.5*zdv*zpe*rho3(nl)
      ci(nl)=0.0
      ce(nl)=0.0
      bi(nl)=0.5*(ddti1(nl)+ddti3(nl))-ai(nl)+zpi3dv
      be(nl)=0.5*(ddte1(nl)+ddte3(nl)+dcvz(nl))-ae(nl)+zpe3dv
c 3. move pressure term for non-perfect ions
 200  IF(.not.(nlpfi))THEN
       DO 250 l=1,nl
        zpi=zkm/xmieff(l)
        zdv=v3(l)-v1(l)
        zpi3dv=0.5*zdv*zpi*rho3(l)
c           subtract the pi3/ti3 * dv term, because pi3 is not nkti3
        bi(l)=bi(l)-zpi3dv
c           we must also subtract this term from di but now we use pi3
        di(l)=di(l)-zdv*pi3(l)*0.5
 250   CONTINUE
      ENDIF
c     4. move pressure term for non-perfect electrons
      IF(nlpfe)RETURN
      DO 300 l=1,nl
       zpe=zkm*fz3(l)/xmieff(l)
       zdv=v3(l)-v1(l)
       zdz=fz3(l)-fz1(l)
       zpe3dv=0.5*zdv*zpe*rho3(l)
c        subtract the pe3/te3 * dv term, because pe3 is not nkte3
       be(l)=be(l)-zpe3dv
c in the degenerate (and tf) case we have combined the terms pe3*dv and
c ddroe*drho, the latter being written as -ddroe/(v*v) * dv
       de(l)=0.5*dt2*(eixch2(l)+ye3(l)+xl3(l)+brems3(l))
     &       +0.5*(te1(l)*(ddte1(l)+ddte3(l)-dcvz(l)))
     &       -0.5*zdv*(ddroe1(l)+ddroe3(l))+ htherm(l)
     &       +dt2*(efbx2(l)+efrx2(l))
       if(abs(piq(58)-1.).lt.underf) then
         de(l)=de(l)  - 0.5 * zdz * (ddz1(l) + ddz3(l))
       endif
 300  CONTINUE
c remember to redefine the coefficient de at the outer boundary in
c case 2 and only if degeneracy is present
      IF(.not.(ncase.eq.2.and.(.not.nlpfe)))RETURN
      zmrc=(r3(njp1)**igeom)*dt2/(dm(nl)*(r3(njp1)-r3(nj))*2.0)
      ce(nl)=-zmrc*xappae(njp1)
c now subtract the term c*t from the rhs of eq.62
      de(nl)=de(nl)-te3(nlp1)*ce(nl)
c reset ce once again
      ce(nl)=0.0
      RETURN
      END


      SUBROUTINE arrays(kgroup,kblock)
      PARAMETER(kk=1001)
c 5.3 print common arrays
c c1.1 basic system parameters
      CHARACTER label1*80,label2*80,label3*80,label4*80,label5*80,
     &          label6*80,label7*80,label8*80
      COMMON /comlab/ label1,label2,label3,label4,label5,label6,label7,
     &                label8
c end of basic system parameters
c c1.9 development and diagnostic parameters
      DIMENSION nadump(20),npdump(20),nvdump(20)
      LOGICAL nlched,nlhead(9),nlomt1(50),nlomt2(50),nlomt3(50),nlrept
      COMMON /comddp/ maxdum,mxdump,nadump,nclass,npdump,npoint,nsub,
     &                nvdump,nlched,nlhead,nlomt1,nlomt2,nlomt3,nlrept
c end of development and diagnostic parameters
c c2.1 hydrodynamics
      DIMENSION dm(kk),p3(kk),r1(kk),r3(kk),r5(kk),rho1(kk),rho3(kk),
     &          rho5(kk),u2(kk),u4(kk),v1(kk),v3(kk),v5(kk),lstat(kk)
      COMMON /comhyd/ dm,p3,pini,r1,r3,r5,rho1,rho3,rho5,rhoini,rhor,
     &                rini,time,u2,u4,uedge,v1,v3,v5
c end of hydrodynamics
c c2.2 thermodynamics
      DIMENSION ddroe1(kk),ddroe3(kk),ddroi1(kk),ddroi3(kk),ddte1(kk),
     &          ddte3(kk),ddti1(kk),ddti3(kk),xappae(kk),xappai(kk),
     &          pe1(kk),pe3(kk),pi1(kk),pi3(kk),te1(kk),te3(kk),ti1(kk),
     &          ti3(kk)
      COMMON /comth / ddroe1,ddroe3,ddroi1,ddroi3,ddte1,ddte3,ddti1,
     &                ddti3,gammae,gammai,xappae,xappai,pe1,pe3,pi1,pi3,
     &                te1,te3,teini,ti1,ti3,tiini
c end of thermodynamics
c c2.3 ions and electrons
      DIMENSION brems1(kk),brems3(kk),ddz1(kk),ddz3(kk),degen(kk),
     &          effz(kk),eixch2(kk),fidash(kk),fz1(kk),fz3(kk),fzsq1(kk)
     &          ,fzsq3(kk),xlc(kk),xmieff(kk),xne(kk),xni(kk),omega1(kk)
      COMMON /comie / brems1,brems3,ddz1,ddz3,degen,degmax,degmin,effz,
     &                eixch2,fidash,fz1,fz3,fzsq1,fzsq3,xlc,xmieff,xne,
     &                xni,omega1,pmass
c end of ions and electrons
c c2.4 laser variables
      DIMENSION alpha1(kk),xaser1(kk),xl1(kk),xl3(kk)
      COMMON /comlas/ alpha1,elas1,xamda1,xaser1,xecri1,plas1,rabs1,
     &                rocri1,xl1,xl3,nabs1
c end of laser variables
c c2.5 thermonuclear reactions
      DIMENSION f1d(kk),f1h(kk),f1he3(kk),f1he4(kk),f1neu(kk),f1ntrl(kk)
     &          ,f1t(kk),f1x(kk),f1x1(kk),f1x2(kk),f1x3(kk),f1x4(kk),
     &          f3d(kk),f3h(kk),f3he3(kk),f3he4(kk),f3neu(kk),f3ntrl(kk)
     &          ,f3t(kk),f3x(kk),f3x1(kk),f3x2(kk),f3x3(kk),f3x4(kk),
     &          r1dd(kk),r1dhe3(kk),r1dt(kk),r3dd(kk),r3dhe3(kk),
     &          r3dt(kk),ye1(kk),ye3(kk),yi1(kk),yi3(kk)
      COMMON /comfus/ deuter,f1d,f1h,f1he3,f1he4,f1neu,f1ntrl,f1t,f1x,
     &                f1x1,f1x2,f1x3,f1x4,f3d,f3h,f3he3,f3he4,f3neu,
     &                f3ntrl,f3t,f3x,f3x1,f3x2,f3x3,f3x4,heliu3,heliu4,
     &                hydrog,xetral,xtrlms,pneut1,pneut3,r1dd,r1dhe3,
     &                r1dt,r3dd,r3dhe3,r3dt,rneut1,rneut3,tinucl,totneu,
     &                tritiu,xmass,xmass1,xmass2,xmass3,xmass4,xtra,
     &                xtra1,xtra2,xtra3,xtra4,xz,xz1,xz2,xz3,xz4,ye1,
     &                ye3,yi1,yi3,yield
c end of thermonuclear reactions
c c3.2 mesh and numerical methods
      DIMENSION ae(kk),ai(kk),be(kk),bi(kk),ce(kk),ci(kk),de(kk),di(kk),
     &          e(kk),f(kk),ge(kk),gi(kk),q2(kk),q4(kk),teite(kk),
     &          tiite(kk),uite(kk)
      COMMON /comnum/ mesh,nj,njm1,njp1,nl,nlm1,nlp1,ae,ai,be,bi,ce,ci,
     &                de,di,e,f,ge,gi,q2,q4,teite,tiite,uite
c end of mesh and numerical methods
c c4.1 - administrative variables
      LOGICAL nldump,nlemp
      DIMENSION piq(302)
      COMMON /comadm/ maxdim,maxrun,ndump,nrep,nldump,nlemp,piq,tstop
c c4.1 - end administrative variables
c c5.1 input-output control variables
      LOGICAL nlfilm,nlhcpy,nlprnt
      DIMENSION buf1(kk),buf2(kk),buf3(kk),buf4(kk),buf5(kk),buf6(kk)
      REAL nprnt,nproc
      COMMON /comout/ nfilm,nhdcpy,np1,np2,np3,nlfilm,nlhcpy,nlprnt,
     &        ztnext,nprnt,nproc,buf1,buf2,buf3,buf4,buf5,buf6,scp,
     &                scr,scrho,scte,scti,sctime,nt1,nt2
c end of input-output control variables
c 1. general olympus data
      IF(.not.nlemp)RETURN
c hot electron dump
      CALL hdump(2)
c select group to be printed
      IF(kgroup.ne.0)THEN
       IF(kgroup.eq.1)THEN
       ELSEIF(kgroup.eq.2)THEN
        GOTO 300
       ELSEIF(kgroup.eq.3)THEN
        GOTO 1000
       ELSEIF(kgroup.eq.4)THEN
        GOTO 1200
       ELSEIF(kgroup.eq.5)THEN
        GOTO 1300
       ELSE
        CALL gotoer(' MEDUSA: arrays: goto error, kgroup is',kgroup)
       ENDIF
      ENDIF
      IF(kblock.ne.0)THEN
       IF(kblock.eq.1)THEN
       ELSEIF(kblock.eq.2.or.kblock.eq.3.or.kblock.eq.4.or.
     &        kblock.eq.5.or.kblock.eq.6.or.kblock.eq.7.or.kblock.eq.8)
     &        THEN
        GOTO 200
       ELSEIF(kblock.eq.9)THEN
        GOTO 100
       ELSE
        CALL gotoer(' MEDUSA: arrays: goto error, kblock is',kblock)
       ENDIF
      ENDIF
c 1.1 basic system parameters
      CALL page
      CALL carray('label1',label1)
      CALL carray('label2',label2)
      CALL carray('label3',label3)
      CALL carray('label4',label4)
      CALL page
      CALL carray('label5',label5)
      CALL carray('label6',label6)
      CALL carray('label7',label7)
      CALL carray('label8',label8)
      IF(kblock.ne.0)RETURN
c 1.9 diagnostics and program development
 100  CALL page
      CALL iarray('nadump  ',nadump,20)
      CALL iarray('npdump  ',npdump,20)
      CALL iarray('nvdump  ',nvdump,20)
      CALL larray('nlhead  ',nlhead,9)
      CALL page
      CALL larray('nlomt1  ',nlomt1,50)
      CALL larray('nlomt2  ',nlomt2,99)
      CALL larray('nlomt3  ',nlomt3,50)
 200  IF(kgroup.ne.0)RETURN
c 2. physics
c select block to be printed
 300  IF(kblock.ne.0)THEN
       IF(kblock.eq.1)THEN
       ELSEIF(kblock.eq.2)THEN
        GOTO 400
       ELSEIF(kblock.eq.3)THEN
        GOTO 500
       ELSEIF(kblock.eq.4)THEN
        GOTO 600
       ELSEIF(kblock.eq.5)THEN
        GOTO 700
       ELSEIF(kblock.eq.6)THEN
        GOTO 800
       ELSEIF(kblock.eq.7)THEN
        GOTO 900
       ELSE
        CALL gotoer(' MEDUSA: arrays: goto error, kblock is',kblock)
       ENDIF
      ENDIF
c 2.1 hydrodynamics
      CALL page
      CALL rarray('dm      ',dm,mesh)
      CALL rarray('p3      ',p3,mesh)
      CALL rarray('r1      ',r1,mesh)
      CALL rarray('r3      ',r3,mesh)
      CALL page
      CALL rarray('r5      ',r5,mesh)
      CALL rarray('rho1    ',rho1,mesh)
      CALL rarray('rho3    ',rho3,mesh)
      CALL rarray('rho5    ',rho5,mesh)
      CALL page
      CALL rarray('u2      ',u2,mesh)
      CALL rarray('u4      ',u4,mesh)
      CALL rarray('v1      ',v1,mesh)
      CALL rarray('v3      ',v3,mesh)
      CALL page
      CALL rarray('v5      ',v5,mesh)
      IF(kblock.ne.0)RETURN
c 2.2 thermodynamics
 400  CALL page
      CALL rarray('ddroe1  ',ddroe1,mesh)
      CALL rarray('ddroe3  ',ddroe3,mesh)
      CALL rarray('ddroi1  ',ddroi1,mesh)
      CALL rarray('ddroi3  ',ddroi3,mesh)
      CALL page
      CALL rarray('ddte1   ',ddte1,mesh)
      CALL rarray('ddte3   ',ddte3,mesh)
      CALL rarray('ddti1   ',ddti1,mesh)
      CALL rarray('ddti3   ',ddti3,mesh)
      CALL page
      CALL rarray('xappae  ',xappae,mesh)
      CALL rarray('xappai  ',xappai,mesh)
      CALL rarray('pe1     ',pe1,mesh)
      CALL rarray('pe3     ',pe3,mesh)
      CALL page
      CALL rarray('pi1     ',pi1,mesh)
      CALL rarray('pi3     ',pi3,mesh)
      CALL rarray('te1     ',te1,mesh)
      CALL rarray('te3     ',te3,mesh)
      CALL page
      CALL rarray('ti1     ',ti1,mesh)
      CALL rarray('ti3     ',ti3,mesh)
      IF(kblock.ne.0)RETURN
c 2.3 ions and electrons
 500  CALL page
      CALL rarray('brems1  ',brems1,mesh)
      CALL rarray('brems3  ',brems3,mesh)
      CALL rarray('degen   ',degen,mesh)
      CALL rarray('effz    ',effz,mesh)
      CALL rarray('fz1     ',fz1,mesh)
      CALL rarray('fz3     ',fz3,mesh)
      CALL rarray('fzsq1   ',fzsq1,mesh)
      CALL rarray('fzsq3   ',fzsq3,mesh)
      CALL page
      CALL rarray('eixch2  ',eixch2,mesh)
      CALL rarray('xlc      ',xlc,mesh)
      CALL rarray('xmieff   ',xmieff,mesh)
      CALL rarray('xne      ',xne,mesh)
      CALL page
      CALL rarray('xni      ',xni,mesh)
      CALL rarray('omega1  ',omega1,mesh)
      IF(kblock.ne.0)RETURN
c 2.4 laser variables
 600  CALL page
      CALL rarray('alpha1  ',alpha1,mesh)
      CALL rarray('xaser1  ',xaser1,15)
      CALL rarray('xl1     ',xl1,mesh)
      CALL rarray('xl3     ',xl3,mesh)
      IF(kblock.ne.0)RETURN
c 2.5 thermonuclear reactions
 700  CALL page
      CALL rarray('f1d     ',f1d,mesh)
      CALL rarray('f3d     ',f3d,mesh)
      CALL rarray('f1h     ',f1h,mesh)
      CALL rarray('f3h     ',f3h,mesh)
      CALL page
      CALL rarray('f1he3   ',f1he3,mesh)
      CALL rarray('f3he3   ',f3he3,mesh)
      CALL rarray('f1he4   ',f1he4,mesh)
      CALL rarray('f3he4   ',f3he4,mesh)
      CALL page
      CALL rarray('f1neu   ',f1neu,mesh)
      CALL rarray('f3neu   ',f3neu,mesh)
      CALL rarray('f1ntrl  ',f1ntrl,mesh)
      CALL rarray('f3ntrl  ',f3ntrl,mesh)
      CALL page
      CALL rarray('f1t     ',f1t,mesh)
      CALL rarray('f3t     ',f3t,mesh)
      CALL rarray('f1x     ',f1x,mesh)
      CALL rarray('f3x     ',f3x,mesh)
      CALL page
      CALL rarray('f1x1    ',f1x1,mesh)
      CALL rarray('f3x1    ',f3x1,mesh)
      CALL rarray('f1x2    ',f1x2,mesh)
      CALL rarray('f3x2    ',f3x2,mesh)
      CALL page
      CALL rarray('f1x3    ',f1x3,mesh)
      CALL rarray('f3x3    ',f3x3,mesh)
      CALL rarray('f1x4    ',f1x4,mesh)
      CALL rarray('f3x4    ',f3x4,mesh)
      CALL page
      CALL rarray('r1dd    ',r1dd,mesh)
      CALL rarray('r3dd    ',r3dd,mesh)
      CALL rarray('r1dhe3  ',r1dhe3,mesh)
      CALL rarray('r3dhe3  ',r3dhe3,mesh)
      CALL page
      CALL rarray('r1dt    ',r1dt,mesh)
      CALL rarray('r3dt    ',r3dt,mesh)
      CALL rarray('ye1     ',ye1,mesh)
      CALL rarray('ye3     ',ye3,mesh)
      CALL page
      CALL rarray('yi1     ',yi1,mesh)
      CALL rarray('yi3     ',yi3,mesh)
      IF(kblock.ne.0)RETURN
c 2.6 energies
 800  IF(kblock.ne.0)RETURN
c 2.7 physics control
 900  IF(kgroup.ne.0)RETURN
c 3. numerical scheme
c select block to be printed
 1000 IF(kblock.ne.0)THEN
       IF(kblock.eq.1)THEN
       ELSEIF(kblock.eq.2)THEN
        GOTO 1100
       ELSE
        CALL gotoer(' MEDUSA: arrays: goto error, kblock is',kblock)
       ENDIF
      ENDIF
c 3.1 numerical control parameters
      IF(kblock.ne.0)RETURN
c 3.2 mesh and numerical methods
 1100 CALL page
      CALL rarray('ae      ',ae,mesh)
      CALL rarray('ai      ',ai,mesh)
      CALL rarray('be      ',be,mesh)
      CALL rarray('bi      ',bi,mesh)
      CALL page
      CALL rarray('ce      ',ce,mesh)
      CALL rarray('ci      ',ci,mesh)
      CALL rarray('de      ',de,mesh)
      CALL rarray('di      ',di,mesh)
      CALL page
      CALL rarray('e       ',e,mesh)
      CALL rarray('f       ',f,mesh)
      CALL rarray('ge      ',ge,mesh)
      CALL rarray('gi      ',gi,mesh)
      CALL page
      CALL rarray('q2      ',q2,mesh)
      CALL rarray('q4      ',q4,mesh)
      CALL rarray('teite   ',teite,mesh)
      CALL rarray('tiite   ',tiite,mesh)
      CALL page
      CALL rarray('uite    ',uite,mesh)
      IF(kgroup.ne.0)RETURN
c 4. house-keeping
c 4.1 administrative variables
 1200 CALL page
      CALL rarray('piq     ',piq,mesh)
      IF(kgroup.ne.0)RETURN
c 5. input-output
c 5.1 input-output control variables
 1300 CALL page
      CALL rarray('buf1    ',buf1,mesh)
      CALL rarray('buf2    ',buf2,mesh)
      CALL rarray('buf3    ',buf3,mesh)
      CALL rarray('buf4    ',buf4,mesh)
      CALL page
      CALL rarray('buf5    ',buf5,mesh)
      CALL rarray('buf6    ',buf6,mesh)
      RETURN
      END


      SUBROUTINE auxval
      PARAMETER(kk=1001,underf=1.E-30)
c 1.5 set auxiliary values
c c1.1 basic system parameters
      LOGICAL nlend,nlres
      COMMON /combas/ altime,cptime,stime,ndiary,nin,nledge,nonlin,nout,
     &                nprint,npunch,nread,nrec,nresm,nrun,nstep,nlend,
     &                nlres
c end of basic system parameters
c c2.1 hydrodynamics
      DIMENSION dm(kk),p3(kk),r1(kk),r3(kk),r5(kk),rho1(kk),rho3(kk),
     &          rho5(kk),u2(kk),u4(kk),v1(kk),v3(kk),v5(kk),lstat(kk)
      COMMON /comhyd/ dm,p3,pini,r1,r3,r5,rho1,rho3,rho5,rhoini,rhor,
     &                rini,time,u2,u4,uedge,v1,v3,v5
c end of hydrodynamics
c c2.2 thermodynamics
      DIMENSION ddroe1(kk),ddroe3(kk),ddroi1(kk),ddroi3(kk),ddte1(kk),
     &          ddte3(kk),ddti1(kk),ddti3(kk),xappae(kk),xappai(kk),
     &          pe1(kk),pe3(kk),pi1(kk),pi3(kk),te1(kk),te3(kk),ti1(kk),
     &          ti3(kk)
      COMMON /comth / ddroe1,ddroe3,ddroi1,ddroi3,ddte1,ddte3,ddti1,
     &                ddti3,gammae,gammai,xappae,xappai,pe1,pe3,pi1,pi3,
     &                te1,te3,teini,ti1,ti3,tiini
c end of thermodynamics
c c2.3 ions and electrons
c end of ions and electrons
c c2.4 laser variables
      DIMENSION alpha1(kk),xaser1(kk),xl1(kk),xl3(kk)
      COMMON /comlas/ alpha1,elas1,xamda1,xaser1,xecri1,plas1,rabs1,
     &                rocri1,xl1,xl3,nabs1
c end of laser variables
c c2.5 thermonuclear reactions
      DIMENSION f1d(kk),f1h(kk),f1he3(kk),f1he4(kk),f1neu(kk),f1ntrl(kk)
     &          ,f1t(kk),f1x(kk),f1x1(kk),f1x2(kk),f1x3(kk),f1x4(kk),
     &          f3d(kk),f3h(kk),f3he3(kk),f3he4(kk),f3neu(kk),f3ntrl(kk)
     &          ,f3t(kk),f3x(kk),f3x1(kk),f3x2(kk),f3x3(kk),f3x4(kk),
     &          r1dd(kk),r1dhe3(kk),r1dt(kk),r3dd(kk),r3dhe3(kk),
     &          r3dt(kk),ye1(kk),ye3(kk),yi1(kk),yi3(kk)
      COMMON /comfus/ deuter,f1d,f1h,f1he3,f1he4,f1neu,f1ntrl,f1t,f1x,
     &                f1x1,f1x2,f1x3,f1x4,f3d,f3h,f3he3,f3he4,f3neu,
     &                f3ntrl,f3t,f3x,f3x1,f3x2,f3x3,f3x4,heliu3,heliu4,
     &                hydrog,xetral,xtrlms,pneut1,pneut3,r1dd,r1dhe3,
     &                r1dt,r3dd,r3dhe3,r3dt,rneut1,rneut3,tinucl,totneu,
     &                tritiu,xmass,xmass1,xmass2,xmass3,xmass4,xtra,
     &                xtra1,xtra2,xtra3,xtra4,xz,xz1,xz2,xz3,xz4,ye1,
     &                ye3,yi1,yi3,yield
c end of thermonuclear reactions
c c2.7 physics control
      LOGICAL mlbrms,mlecon,mlfuse,mlicon,mlx,nlabs,nlbrms,nlburn,
     &        nlcri1,nldepo,nlecon,nlfuse,nlicon,nlmove,nlpfe,nlpfi,nlx,
     &        nlhf,nltnl
      COMMON /comcon/ mstep,ncase,ngeom,mlbrms,mlecon,mlfuse,mlicon,mlx,
     &                nlabs,nlbrms,nlburn,nlcri1,nldepo,nlecon,nlfuse,
     &                nlicon,nlmove,nlpfe,nlpfi,nlx,nlhf,nltnl,pmin,
     &                rhomin,temin,timin,umin
c end of physics control
c c3.1 numerical control parameters
      LOGICAL break,nlgoon,nlite
      COMMON /comnc / nceldt,ncondt,nit,nitmax,break,nlgoon,nlite,ak0,
     &                ak1,ak2,ak3,ak4,ak5,ak6,ak7,ak8,ak9,bneum,deltat,
     &                dt2,dt3,dt4,dtemax,dtimax,dtmax,dtmin,dumax,rdt2,
     &                rdt3,rdt4
c end of numerical control parameters
c c3.2 mesh and numerical methods
      DIMENSION ae(kk),ai(kk),be(kk),bi(kk),ce(kk),ci(kk),de(kk),di(kk),
     &          e(kk),f(kk),ge(kk),gi(kk),q2(kk),q4(kk),teite(kk),
     &          tiite(kk),uite(kk)
      COMMON /comnum/ mesh,nj,njm1,njp1,nl,nlm1,nlp1,ae,ai,be,bi,ce,ci,
     &                de,di,e,f,ge,gi,q2,q4,teite,tiite,uite
c end of mesh and numerical methods
c c4.1 - administrative variables
      LOGICAL nldump,nlemp
      DIMENSION piq(302)
      COMMON /comadm/ maxdim,maxrun,ndump,nrep,nldump,nlemp,piq,tstop
c c4.1 - end administrative variables
c c5.1 input-output control variables
      LOGICAL nlfilm,nlhcpy,nlprnt
      DIMENSION buf1(kk),buf2(kk),buf3(kk),buf4(kk),buf5(kk),buf6(kk)
      REAL nprnt,nproc
      COMMON /comout/ nfilm,nhdcpy,np1,np2,np3,nlfilm,nlhcpy,nlprnt,
     &        ztnext,nprnt,nproc,buf1,buf2,buf3,buf4,buf5,buf6,scp,
     &                scr,scrho,scte,scti,sctime,nt1,nt2
c end of input-output control variables
c 2. physics
c 2.1 hydrodynamics
      rini=rini/scr
      rhoini=rhoini/scrho
c 2.2 thermodynamics
      tiini=tiini/scti
      teini=teini/scte
c 2.3 ions and electrons
c 2.4 laser variables
      xamda1=xamda1/scr
      zc2psq=(2.0*3.14159265*3.0*1.0E8)**2
c zmi is an average ion mass number
      zmi=deuter*2.0+hydrog*1.0+tritiu*3.0+heliu3*3.0+heliu4*4.0+
     &    xetral*xtrlms+xtra*xmass
      z=deuter+hydrog+tritiu+2.0*(heliu3+heliu4)+xz*xtra
c find the critical density from equation 25
      rocri1=5.2*1.0E-31*zc2psq*zmi/(xamda1*xamda1*z)
      xecri1=1.115E+15/(xamda1*xamda1)
c 2.5 thermonuclear reactions
      tinucl=tinucl/scti
c check that charge number of extra element does not exceed xmass
      IF(xz.gt.xmass)THEN
c 6. faulty input data. terminate.
       CALL mesage('    charge number of extra element too large')
       CALL clist(2,5)
      ELSE
c check that fractions do not exceed unity
       IF(deuter.le.1.0)THEN
        IF(heliu3.le.1.0)THEN
         IF(heliu4.le.1.0)THEN
          IF(hydrog.le.1.0)THEN
           IF(xetral.le.1.0)THEN
            IF(tritiu.le.1.0)THEN
             IF(xtra.le.1.0)THEN
              IF(abs(deuter).lt.underf.and.nlfuse)THEN
               CALL mesage(
     &                    '    no deuterium for thermonuclear reactions'
     &                    )
               CALL clist(2,5)
               CALL clist(2,7)
               GOTO 100
              ELSE
c 2.6 energies
c 2.7 physics control
c threshold values depend on initial values
               rhomin=1.0E-14*rhoini
c ad1
               timin=1.0
               temin=1.0
c use a typical soundspeed to determine umin
               zkm=8255.17
               umin=1.0E-10*sqrt(zkm*(teini/scte+tiini/scti))
c assign inner switches for particle processes
               mlbrms=nlbrms
               mlecon=nlecon
               mlfuse=nlfuse
               mlicon=nlicon
               mlx=nlx
               mlx=nlx
c if fusion is switched off there is no burn-up or deposition
               IF(.not.nlfuse)nlburn=.false.
               IF(.not.nlfuse)nldepo=.false.
               IF(ncase.eq.1)THEN
                IF(((.not.nlabs).and.(.not.nlcri1)).and.
     &             (abs(xaser1(2)).gt.underf))THEN
                 CALL mesage('    no absorption mechanisms')
                 CALL clist(2,4)
                 CALL arrays(2,4)
                 CALL clist(2,7)
                 GOTO 100
                ENDIF
               ENDIF
c 3. numerical scheme
c 3.1 numerical control parameters
               deltat=deltat/sctime
c 3.2 mesh and numerical methods
               nj=mesh
               njp1=nj+1
               njm1=nj-1
               nl=mesh
               nlp1=nl+1
               nlm1=nl-1
c 4. house-keeping
c 4.1 administrative variables
               tstop=tstop/sctime
c check for buffer overflow
               IF(nrun.gt.maxrun.or.mesh.gt.maxdim)THEN
                CALL mesage('    buffer overflow')
                CALL clist(4,1)
                GOTO 100
               ELSE
c 5. input-output
c 5.1 input-output control variables
                np1=1
                np2=mesh
                RETURN
               ENDIF
              ENDIF
             ENDIF
            ENDIF
           ENDIF
          ENDIF
         ENDIF
        ENDIF
       ENDIF
       CALL mesage('    initial fraction(s) larger than unity')
       CALL clist(2,5)
      ENDIF
 100  CALL endrun
      RETURN
      END


      SUBROUTINE basic
c 0.1 initialize basic data
c c1.1 basic system parameters
      CHARACTER label1*80,label2*80,label3*80,label4*80,label5*80,
     &          label6*80,label7*80,label8*80
      LOGICAL nlend,nlres
      COMMON /comlab/ label1,label2,label3,label4,label5,label6,label7,
     &                label8
      COMMON /combas/ altime,cptime,stime,ndiary,nin,nledge,nonlin,nout,
     &                nprint,npunch,nread,nrec,nresm,nrun,nstep,nlend,
     &                nlres
c end of basic system parameters
c c1.9 development and diagnostic parameters
      DIMENSION nadump(20),npdump(20),nvdump(20)
      LOGICAL nlched,nlhead(9),nlomt1(50),nlomt2(50),nlomt3(50),nlrept
      COMMON /comddp/ maxdum,mxdump,nadump,nclass,npdump,npoint,nsub,
     &                nvdump,nlched,nlhead,nlomt1,nlomt2,nlomt3,nlrept
c end of development and diagnostic parameters
      CHARACTER iblank*(80)
      DATA iblank/'
     &                           '/
c 1. general olympus data
c 1.1 basic system parameters
c cpu - time used so far
      cptime=0.0
c clear all 8 label arrays
      CALL resetc(label1,iblank)
c input-output channels
      nledge=15
      nonlin=1
c     npunch    = 48
      npunch=6
      nprint=6
      nread=5
      ndiary=npunch
      nin=nread
      nout=nprint
c timestep control
      nrun=1
      nstep=0
c restart control
      nrec=1
      nresm=nledge
c logical switches
      nlend=.false.
      nlres=.false.
c 1.9 diagnostic and development parameters
c maximum dimensions of dump arrays
      maxdum=20
      mxdump=10
c reset dump arrays
      CALL reseti(nadump,maxdum,0)
      CALL reseti(npdump,maxdum,0)
      CALL reseti(nvdump,maxdum,0)
c tracer variables
      nclass=0
      nsub=1
      npoint=1
c logical switches
      nlched=.false.
      nlrept=.false.
c report heads for classes 1-9
      CALL resetl(nlhead,9,.false.)
c reset class 1, 2, 3 subprogram selector array
      CALL resetl(nlomt1,50,.false.)
      CALL resetl(nlomt2,50,.false.)
      CALL resetl(nlomt3,50,.false.)
c user interface
      CALL modify
      RETURN
      END


      SUBROUTINE blines(k)
c u.3 insert blank lines on output channel
c c1.1 basic system parameters
      LOGICAL nlend,nlres
      COMMON /combas/ altime,cptime,stime,ndiary,nin,nledge,nonlin,nout,
     &                nprint,npunch,nread,nrec,nresm,nrun,nstep,nlend,
     &                nlres
c end of basic system parameters
      DO 100 j=1,k
       WRITE(nout,10100)
 100  CONTINUE
      RETURN
10100 FORMAT(' ')
      END


      SUBROUTINE boundy(tnow)
      PARAMETER(kk=1001,underf=1.E-30)
c 2.25 sets the boundary values at time tnow
c c1.9 development and diagnostic parameters
      DIMENSION nadump(20),npdump(20),nvdump(20)
      LOGICAL nlched,nlhead(9),nlomt1(50),nlomt2(50),nlomt3(50),nlrept
      COMMON /comddp/ maxdum,mxdump,nadump,nclass,npdump,npoint,nsub,
     &                nvdump,nlched,nlhead,nlomt1,nlomt2,nlomt3,nlrept
c end of development and diagnostic parameters
c c1.1 basic system parameters
c end of basic system parameters
c c2.1 hydrodynamics
      DIMENSION dm(kk),p3(kk),r1(kk),r3(kk),r5(kk),rho1(kk),rho3(kk),
     &          rho5(kk),u2(kk),u4(kk),v1(kk),v3(kk),v5(kk),lstat(kk)
      COMMON /comhyd/ dm,p3,pini,r1,r3,r5,rho1,rho3,rho5,rhoini,rhor,
     &                rini,time,u2,u4,uedge,v1,v3,v5
c end of hydrodynamics
c c2.2 thermodynamics
      DIMENSION ddroe1(kk),ddroe3(kk),ddroi1(kk),ddroi3(kk),ddte1(kk),
     &          ddte3(kk),ddti1(kk),ddti3(kk),xappae(kk),xappai(kk),
     &          pe1(kk),pe3(kk),pi1(kk),pi3(kk),te1(kk),te3(kk),ti1(kk),
     &          ti3(kk)
      COMMON /comth / ddroe1,ddroe3,ddroi1,ddroi3,ddte1,ddte3,ddti1,
     &                ddti3,gammae,gammai,xappae,xappai,pe1,pe3,pi1,pi3,
     &                te1,te3,teini,ti1,ti3,tiini
c end of thermodynamics
c c2.3 ions and electrons
c end of ions and electrons
c c2.4 laser variables
c end of laser variables
c c2.7 physics control
      LOGICAL mlbrms,mlecon,mlfuse,mlicon,mlx,nlabs,nlbrms,nlburn,
     &        nlcri1,nldepo,nlecon,nlfuse,nlicon,nlmove,nlpfe,nlpfi,nlx,
     &        nlhf,nltnl
      COMMON /comcon/ mstep,ncase,ngeom,mlbrms,mlecon,mlfuse,mlicon,mlx,
     &                nlabs,nlbrms,nlburn,nlcri1,nldepo,nlecon,nlfuse,
     &                nlicon,nlmove,nlpfe,nlpfi,nlx,nlhf,nltnl,pmin,
     &                rhomin,temin,timin,umin
c end of physics control
c c3.2 mesh and numerical methods
      DIMENSION ae(kk),ai(kk),be(kk),bi(kk),ce(kk),ci(kk),de(kk),di(kk),
     &          e(kk),f(kk),ge(kk),gi(kk),q2(kk),q4(kk),teite(kk),
     &          tiite(kk),uite(kk)
      COMMON /comnum/ mesh,nj,njm1,njp1,nl,nlm1,nlp1,ae,ai,be,bi,ce,ci,
     &                de,di,e,f,ge,gi,q2,q4,teite,tiite,uite
c end of mesh and numerical methods
c c4.1 - administrative variables
      LOGICAL nldump,nlemp
      DIMENSION piq(302)
      COMMON /comadm/ maxdim,maxrun,ndump,nrep,nldump,nlemp,piq,tstop
c c4.1 - end administrative variables
      DATA iclass,isub/2,25/
      IF(nlomt2(isub))RETURN
      IF(ncase.eq.1)THEN
      ELSEIF(ncase.eq.2)THEN
c 2. temperatures at immoveable wall (eq.38)
       ti3(nlp1)=ti3(nl)
c eq.92
       te3(nlp1)=min(piq(5),teini+piq(5)*piq(3)*tnow)
       RETURN
      ELSEIF(ncase.eq.3)THEN
c 3. pressure applied (eq.93)
       zft=1.0-tnow/piq(3)
       IF(zft.lt.0.0)RETURN
       p3(nlp1)=pini*min(piq(6),zft**piq(4))
       RETURN
      ELSEIF(ncase.eq.4)THEN
       GOTO 100
      ELSE
       CALL gotoer(' MEDUSA: boundy: goto error, ncase is',ncase)
      ENDIF
c standard  test  cases
c 1. thermally insulated boundary
      IF(abs(piq(8)).gt.underf)THEN
c the edge pressure is now a function of the pressure just inside
       zpsave=0.
       IF(abs(tnow).lt.underf)zpsave=p3(nl)
       p3(nlp1)=zpsave-(p3(nl)-zpsave)
c when the edge pressure approaches zero we keep it constant
       IF(p3(nlp1).ge.pmin)RETURN
      ENDIF
c the hydrodynamic pressure is zero (or very close to) (eq.38)
      p3(nlp1)=pmin
      RETURN
c 4. velocity determined
 100  u4(njp1)=piq(2)*sqrt(p3(nl)*v3(nl)*gammai)
      RETURN
      END


      SUBROUTINE brems
      PARAMETER(kk=1001)
c 2.10 the loss rate through bremsstrahlung
c c1.9 development and diagnostic parameters
      DIMENSION nadump(20),npdump(20),nvdump(20)
      LOGICAL nlched,nlhead(9),nlomt1(50),nlomt2(50),nlomt3(50),nlrept
      COMMON /comddp/ maxdum,mxdump,nadump,nclass,npdump,npoint,nsub,
     &                nvdump,nlched,nlhead,nlomt1,nlomt2,nlomt3,nlrept
c end of development and diagnostic parameters
c c2.1 hydrodynamics
      DIMENSION dm(kk),p3(kk),r1(kk),r3(kk),r5(kk),rho1(kk),rho3(kk),
     &          rho5(kk),u2(kk),u4(kk),v1(kk),v3(kk),v5(kk),lstat(kk)
      COMMON /comhyd/ dm,p3,pini,r1,r3,r5,rho1,rho3,rho5,rhoini,rhor,
     &                rini,time,u2,u4,uedge,v1,v3,v5
c end of hydrodynamics
c c2.2 thermodynamics
      DIMENSION ddroe1(kk),ddroe3(kk),ddroi1(kk),ddroi3(kk),ddte1(kk),
     &          ddte3(kk),ddti1(kk),ddti3(kk),xappae(kk),xappai(kk),
     &          pe1(kk),pe3(kk),pi1(kk),pi3(kk),te1(kk),te3(kk),ti1(kk),
     &          ti3(kk)
      COMMON /comth / ddroe1,ddroe3,ddroi1,ddroi3,ddte1,ddte3,ddti1,
     &                ddti3,gammae,gammai,xappae,xappai,pe1,pe3,pi1,pi3,
     &                te1,te3,teini,ti1,ti3,tiini
c end of thermodynamics
c c2.3 ions and electrons
      DIMENSION brems1(kk),brems3(kk),ddz1(kk),ddz3(kk),degen(kk),
     &          effz(kk),eixch2(kk),fidash(kk),fz1(kk),fz3(kk),fzsq1(kk)
     &          ,fzsq3(kk),xlc(kk),xmieff(kk),xne(kk),xni(kk),omega1(kk)
      COMMON /comie / brems1,brems3,ddz1,ddz3,degen,degmax,degmin,effz,
     &                eixch2,fidash,fz1,fz3,fzsq1,fzsq3,xlc,xmieff,xne,
     &                xni,omega1,pmass
c end of ions and electrons
c c2.5 thermonuclear reactions
      DIMENSION f1d(kk),f1h(kk),f1he3(kk),f1he4(kk),f1neu(kk),f1ntrl(kk)
     &          ,f1t(kk),f1x(kk),f1x1(kk),f1x2(kk),f1x3(kk),f1x4(kk),
     &          f3d(kk),f3h(kk),f3he3(kk),f3he4(kk),f3neu(kk),f3ntrl(kk)
     &          ,f3t(kk),f3x(kk),f3x1(kk),f3x2(kk),f3x3(kk),f3x4(kk),
     &          r1dd(kk),r1dhe3(kk),r1dt(kk),r3dd(kk),r3dhe3(kk),
     &          r3dt(kk),ye1(kk),ye3(kk),yi1(kk),yi3(kk)
      COMMON /comfus/ deuter,f1d,f1h,f1he3,f1he4,f1neu,f1ntrl,f1t,f1x,
     &                f1x1,f1x2,f1x3,f1x4,f3d,f3h,f3he3,f3he4,f3neu,
     &                f3ntrl,f3t,f3x,f3x1,f3x2,f3x3,f3x4,heliu3,heliu4,
     &                hydrog,xetral,xtrlms,pneut1,pneut3,r1dd,r1dhe3,
     &                r1dt,r3dd,r3dhe3,r3dt,rneut1,rneut3,tinucl,totneu,
     &                tritiu,xmass,xmass1,xmass2,xmass3,xmass4,xtra,
     &                xtra1,xtra2,xtra3,xtra4,xz,xz1,xz2,xz3,xz4,ye1,
     &                ye3,yi1,yi3,yield
c end of thermonuclear reactions
c c2.7 physics control
      LOGICAL mlbrms,mlecon,mlfuse,mlicon,mlx,nlabs,nlbrms,nlburn,
     &        nlcri1,nldepo,nlecon,nlfuse,nlicon,nlmove,nlpfe,nlpfi,nlx,
     &        nlhf,nltnl
      COMMON /comcon/ mstep,ncase,ngeom,mlbrms,mlecon,mlfuse,mlicon,mlx,
     &                nlabs,nlbrms,nlburn,nlcri1,nldepo,nlecon,nlfuse,
     &                nlicon,nlmove,nlpfe,nlpfi,nlx,nlhf,nltnl,pmin,
     &                rhomin,temin,timin,umin
c end of physics control
c c3.1 numerical control parameters
      LOGICAL break,nlgoon,nlite
      COMMON /comnc / nceldt,ncondt,nit,nitmax,break,nlgoon,nlite,ak0,
     &                ak1,ak2,ak3,ak4,ak5,ak6,ak7,ak8,ak9,bneum,deltat,
     &                dt2,dt3,dt4,dtemax,dtimax,dtmax,dtmin,dumax,rdt2,
     &                rdt3,rdt4
c end of numerical control parameters
c c3.2 mesh and numerical methods
      DIMENSION ae(kk),ai(kk),be(kk),bi(kk),ce(kk),ci(kk),de(kk),di(kk),
     &          e(kk),f(kk),ge(kk),gi(kk),q2(kk),q4(kk),teite(kk),
     &          tiite(kk),uite(kk)
      COMMON /comnum/ mesh,nj,njm1,njp1,nl,nlm1,nlp1,ae,ai,be,bi,ce,ci,
     &                de,di,e,f,ge,gi,q2,q4,teite,tiite,uite
c end of mesh and numerical methods
c c4.1 - administrative variables
      LOGICAL nldump,nlemp
      DIMENSION piq(302)
      COMMON /comadm/ maxdim,maxrun,ndump,nrep,nldump,nlemp,piq,tstop
c c4.1 - end administrative variables
      DATA iclass,isub/2,10/
      IF(nlomt2(isub))RETURN
c 1. the bremsstrahlung rate of the electrons
      an=6.0232E23
c set the constant of eq.24 (n.b. brems is negative)
      zg=gammae-1.0
      zj=-9.05*1.0E-14*(1.0+piq(24))
c charge numbers squared
      zzd=1.0
      zzh=1.0
      zzhe3=4.0
      zzhe4=4.0
      zzt=1.0
      zzx=xz*xz
      DO 100 l=1,nl
       IF(.not.nlpfe)THEN
c 1.2 a weakly degenerate electron gas
c for partially degenerate electrons te3 determines the radiation
        IF(degen(l).gt.degmin)THEN
         zt3=te3(l)
c 1.3 a strongly degenerate electron gas
c the departure from the adiabat determines the radiation
         GOTO 50
        ENDIF
       ENDIF
c 1.1 a perfect electron gas
c evaluate the departure from the adiabat (eq.41)
       rho0=rhoini
       tein1=teini
       IF(l.gt.mesh-piq(63).and.piq(63).gt.0) then
            rho0=piq(62)
            tein1=piq(89)
       elseIF(l.gt.mesh-piq(63)-piq(13).and.piq(13).gt.0) then
           rho0=piq(12)
            teini1=piq(86)
       endif
       zt3=te3(l)-tein1*(rho3(l)/rho0)**zg
       IF(zt3.le.0.0)zt3=0.0
c 1.4 the bremsstrahlung rate (eqs.24 or 40)
c the average of z**2
 50    brems3(l)=zj*xne(l)*sqrt(zt3)*fzsq3(l)/xmieff(l)
c 1.5 limitations on the bremsstrahlung rate the energy radiated away
       zerad=-brems3(l)*dt2
c the energy available for radiation
       zavail=ddte3(l)*zt3*0.5
c restrict the bremsstrahlung rate if necessary
       IF(zerad.ge.zavail)brems3(l)=-zavail*rdt2
 100  CONTINUE
      RETURN
      END


      SUBROUTINE burnup
      PARAMETER(kk=1001)
c 2.28 burn-up of deuterium and tritium
c c1.9 development and diagnostic parameters
      DIMENSION nadump(20),npdump(20),nvdump(20)
      LOGICAL nlched,nlhead(9),nlomt1(50),nlomt2(50),nlomt3(50),nlrept
      COMMON /comddp/ maxdum,mxdump,nadump,nclass,npdump,npoint,nsub,
     &                nvdump,nlched,nlhead,nlomt1,nlomt2,nlomt3,nlrept
c end of development and diagnostic parameters
c c2.1 hydrodynamics
      DIMENSION dm(kk),p3(kk),r1(kk),r3(kk),r5(kk),rho1(kk),rho3(kk),
     &          rho5(kk),u2(kk),u4(kk),v1(kk),v3(kk),v5(kk),lstat(kk)
      COMMON /comhyd/ dm,p3,pini,r1,r3,r5,rho1,rho3,rho5,rhoini,rhor,
     &                rini,time,u2,u4,uedge,v1,v3,v5
c end of hydrodynamics
c c2.3 ions and electrons
      DIMENSION brems1(kk),brems3(kk),ddz1(kk),ddz3(kk),degen(kk),
     &          effz(kk),eixch2(kk),fidash(kk),fz1(kk),fz3(kk),fzsq1(kk)
     &          ,fzsq3(kk),xlc(kk),xmieff(kk),xne(kk),xni(kk),omega1(kk)
      COMMON /comie / brems1,brems3,ddz1,ddz3,degen,degmax,degmin,effz,
     &                eixch2,fidash,fz1,fz3,fzsq1,fzsq3,xlc,xmieff,xne,
     &                xni,omega1,pmass
c end of ions and electrons
c c2.5 thermonuclear reactions
      DIMENSION f1d(kk),f1h(kk),f1he3(kk),f1he4(kk),f1neu(kk),f1ntrl(kk)
     &          ,f1t(kk),f1x(kk),f1x1(kk),f1x2(kk),f1x3(kk),f1x4(kk),
     &          f3d(kk),f3h(kk),f3he3(kk),f3he4(kk),f3neu(kk),f3ntrl(kk)
     &          ,f3t(kk),f3x(kk),f3x1(kk),f3x2(kk),f3x3(kk),f3x4(kk),
     &          r1dd(kk),r1dhe3(kk),r1dt(kk),r3dd(kk),r3dhe3(kk),
     &          r3dt(kk),ye1(kk),ye3(kk),yi1(kk),yi3(kk)
      COMMON /comfus/ deuter,f1d,f1h,f1he3,f1he4,f1neu,f1ntrl,f1t,f1x,
     &                f1x1,f1x2,f1x3,f1x4,f3d,f3h,f3he3,f3he4,f3neu,
     &                f3ntrl,f3t,f3x,f3x1,f3x2,f3x3,f3x4,heliu3,heliu4,
     &                hydrog,xetral,xtrlms,pneut1,pneut3,r1dd,r1dhe3,
     &                r1dt,r3dd,r3dhe3,r3dt,rneut1,rneut3,tinucl,totneu,
     &                tritiu,xmass,xmass1,xmass2,xmass3,xmass4,xtra,
     &                xtra1,xtra2,xtra3,xtra4,xz,xz1,xz2,xz3,xz4,ye1,
     &                ye3,yi1,yi3,yield
c end of thermonuclear reactions
c c3.1 numerical control parameters
      LOGICAL break,nlgoon,nlite
      COMMON /comnc / nceldt,ncondt,nit,nitmax,break,nlgoon,nlite,ak0,
     &                ak1,ak2,ak3,ak4,ak5,ak6,ak7,ak8,ak9,bneum,deltat,
     &                dt2,dt3,dt4,dtemax,dtimax,dtmax,dtmin,dumax,rdt2,
     &                rdt3,rdt4
c end of numerical control parameters
c c3.2 mesh and numerical methods
      DIMENSION ae(kk),ai(kk),be(kk),bi(kk),ce(kk),ci(kk),de(kk),di(kk),
     &          e(kk),f(kk),ge(kk),gi(kk),q2(kk),q4(kk),teite(kk),
     &          tiite(kk),uite(kk)
      COMMON /comnum/ mesh,nj,njm1,njp1,nl,nlm1,nlp1,ae,ai,be,bi,ce,ci,
     &                de,di,e,f,ge,gi,q2,q4,teite,tiite,uite
c end of mesh and numerical methods
      DATA iclass,isub/2,28/
      IF(nlomt2(isub))RETURN
c 1. burn-up of deuterium and tritium (eq.30)
c standard ion masses
      zmd=2.0
      zmh=1.0
      zmhe3=3.0
      zmhe4=4.0
      zmntrl=xtrlms
      zmt=3.0
      zmx=xmass
c     the fractional changes at levels 1 and 3 are with respect to the
c     number densities at these levels. (include dt2/2)
c     effective mass number at level 1
      DO 100 l=1,nl
       zm1=zmd*f1d(l)+zmh*f1h(l)+zmhe3*f1he3(l)+zmhe4*f1he4(l)
     &     +zmntrl*f1ntrl(l)+zmt*f1t(l)+zmx*f1x(l)
c     evaluate dt2 / (2*ni1). ni1 is obtained from f1 and rho1
       zni1=0.5*dt2*pmass*zm1/rho1(l)
c     set sums for fractional changes to zero
       zdfd=0.0
       zdfhe3=0.0
       zdfhe4=0.0
       zdfh=0.0
       zdfneu=0.0
       zdft=0.0
c 1.1 d-d process
       zfdd=(r1dd(l)+r3dd(l))*zni1
       zdfd=zdfd-2.0*zfdd
       zdfh=zdfh+0.5*zfdd
       zdfhe3=zdfhe3+0.5*zfdd
       zdfneu=zdfneu+0.5*zfdd
       zdft=zdft+0.5*zfdd
c 1.2 d-t process
       zfdt=(r1dt(l)+r3dt(l))*zni1
       zdfd=zdfd-zfdt
       zdfhe4=zdfhe4+zfdt
       zdfneu=zdfneu+zfdt
       zdft=zdft-zfdt
c 1.3 d-he3 process
       zfdhe3=(r1dhe3(l)+r3dhe3(l))*zni1
       zdfd=zdfd-zfdhe3
       zdfh=zdfh+zfdhe3
       zdfhe3=zdfhe3-zfdhe3
       zdfhe4=zdfhe4+zfdhe3
c 1.4 sum up contributions from all processes
       f3d(l)=f1d(l)+zdfd
       f3h(l)=f1h(l)+zdfh
       f3he3(l)=f1he3(l)+zdfhe3
       f3he4(l)=f1he4(l)+zdfhe4
       f3neu(l)=f1neu(l)+zdfneu
       f3t(l)=f1t(l)+zdft
 100  CONTINUE
c 2. renormalize fractions so that sum(f) = 1
      DO 200 l=1,nl
c apply equation 2, but omit the neutron fraction
       zsumf=f3d(l)+f3h(l)+f3he3(l)+f3he4(l)+f3ntrl(l)+f3t(l)+f3x(l)
     &       +f3x1(l)+f3x2(l)+f3x3(l)+f3x4(l)
c     renormalize fractions
       f3d(l)=f3d(l)/zsumf
       f3h(l)=f3h(l)/zsumf
       f3he3(l)=f3he3(l)/zsumf
       f3he4(l)=f3he4(l)/zsumf
       f3ntrl(l)=f3ntrl(l)/zsumf
       f3t(l)=f3t(l)/zsumf
       f3x(l)=f3x(l)/zsumf
       f3x1(l)=f3x1(l)/zsumf
       f3x2(l)=f3x2(l)/zsumf
       f3x3(l)=f3x3(l)/zsumf
       f3x4(l)=f3x4(l)/zsumf
 200  CONTINUE
      RETURN
      END


      SUBROUTINE carray(kname,ka)
c u.10 print name and value of character string
c c1.1 basic system parameters
      LOGICAL nlend,nlres
      COMMON /combas/ altime,cptime,stime,ndiary,nin,nledge,nonlin,nout,
     &                nprint,npunch,nread,nrec,nresm,nrun,nstep,nlend,
     &                nlres
c end of basic system parameters
      CHARACTER kname*(*),ka*(*)
      CALL blines(1)
      WRITE(nout,10100)kname
      CALL blines(1)
      WRITE(nout,10200)ka
      RETURN
10100 FORMAT(4x,a)
10200 FORMAT(6x,'''',a,'''')
      END


      SUBROUTINE chinic
      PARAMETER(kk=1001)
c 1.10 check that initial conditions are meaningful
c c1.1 basic system parameters
      LOGICAL nlend,nlres
      COMMON /combas/ altime,cptime,stime,ndiary,nin,nledge,nonlin,nout,
     &                nprint,npunch,nread,nrec,nresm,nrun,nstep,nlend,
     &                nlres
c end of basic system parameters
c c2.1 hydrodynamics
      DIMENSION dm(kk),p3(kk),r1(kk),r3(kk),r5(kk),rho1(kk),rho3(kk),
     &          rho5(kk),u2(kk),u4(kk),v1(kk),v3(kk),v5(kk),lstat(kk)
      COMMON /comhyd/ dm,p3,pini,r1,r3,r5,rho1,rho3,rho5,rhoini,rhor,
     &                rini,time,u2,u4,uedge,v1,v3,v5
c end of hydrodynamics
c c2.7 physics control
c end of physics control
c c3.1 numerical control parameters
      LOGICAL break,nlgoon,nlite
      COMMON /comnc / nceldt,ncondt,nit,nitmax,break,nlgoon,nlite,ak0,
     &                ak1,ak2,ak3,ak4,ak5,ak6,ak7,ak8,ak9,bneum,deltat,
     &                dt2,dt3,dt4,dtemax,dtimax,dtmax,dtmin,dumax,rdt2,
     &                rdt3,rdt4
c end of numerical control parameters
c c3.2 mesh and numerical methods
      DIMENSION ae(kk),ai(kk),be(kk),bi(kk),ce(kk),ci(kk),de(kk),di(kk),
     &          e(kk),f(kk),ge(kk),gi(kk),q2(kk),q4(kk),teite(kk),
     &          tiite(kk),uite(kk)
      COMMON /comnum/ mesh,nj,njm1,njp1,nl,nlm1,nlp1,ae,ai,be,bi,ce,ci,
     &                de,di,e,f,ge,gi,q2,q4,teite,tiite,uite
c end of mesh and numerical methods
c 1. first check initial volumes, masses etc.
      l=0
c     the volumes must be positive
 100  l=l+1
      IF(dm(l).gt.0.0)THEN
       i=l
       IF(v1(i).gt.0.0)THEN
        IF(l.lt.nl)GOTO 100
c     it is necessary to copy coordinates to level 5 in order to
c     call subprogram exam
        CALL copyr(r1,1,r5,1,njp1)
c 2. now examine coordinates & temperatures
        CALL exam
c appropriate particle processes have now been switched on,
c were coordinates right?
        IF(.not.(break))GOTO 200
       ENDIF
      ENDIF
c 3. terminate the run
      CALL page
      CALL blines(2)
      CALL mesage('    initial conditions incorrect ***************')
      CALL blines(2)
      IF(.not.break)THEN
       zr1=r1(l-1)
       zr2=r1(l)
       zdm=dm(l)
       zv=v1(l)
       l1=l-1
       WRITE(nout,10100)l1,zr1,l,zr2,l,zdm,i,zv
      ELSE
       CALL mesage('    coordinates incorrectly placed')
      ENDIF
      CALL clist(0,0)
      CALL endrun
c 4. further checks can be made here
 200  RETURN
10100 FORMAT(5x,'r(',i3,')=',e10.3,5x,'r(',i3,')=',e10.3,5x,'v(',i3,
     &       ')=',e10.3,5x,'dm(',i3,')=',e10.3)
      END


      SUBROUTINE clear
      PARAMETER(kk=1001)
c 1.2 clear variables and arrays
c c2.1 hydrodynamics
      DIMENSION dm(kk),p3(kk),r1(kk),r3(kk),r5(kk),rho1(kk),rho3(kk),
     &          rho5(kk),u2(kk),u4(kk),v1(kk),v3(kk),v5(kk),lstat(kk)
      COMMON /comhyd/ dm,p3,pini,r1,r3,r5,rho1,rho3,rho5,rhoini,rhor,
     &                rini,time,u2,u4,uedge,v1,v3,v5
c end of hydrodynamics
      LOGICAL lphi(kk)
      DIMENSION deltae(50),ebdy(51),emid(50),ha(kk),hb(kk),hc(kk),hd(kk)
     &          ,hdiff(kk),he(kk),heta(kk),hf(kk),hg(kk),hgrp(kk),hj(kk)
     &          ,hlas(kk),hloss(50),hncold(kk),hp(kk),hphi(kk),htaug(kk)
     &          ,htherm(kk),htr(kk),xhot1(kk,50),xhot3(kk,50)
      COMMON /comhot/ igrp,ngroup,lphi,deltae,ebdy,ehot,emid,ha,hb,hc,
     &                hd,hdiff,he,heta,hf,hg,hgrp,hj,hlas,hloss,hncold,
     &                hp,hphi,htaug,htherm,htr,xhot1,xhot3,thot
c c2.2 thermodynamics
      DIMENSION ddroe1(kk),ddroe3(kk),ddroi1(kk),ddroi3(kk),ddte1(kk),
     &          ddte3(kk),ddti1(kk),ddti3(kk),xappae(kk),xappai(kk),
     &          pe1(kk),pe3(kk),pi1(kk),pi3(kk),te1(kk),te3(kk),ti1(kk),
     &          ti3(kk)
      COMMON /comth / ddroe1,ddroe3,ddroi1,ddroi3,ddte1,ddte3,ddti1,
     &                ddti3,gammae,gammai,xappae,xappai,pe1,pe3,pi1,pi3,
     &                te1,te3,teini,ti1,ti3,tiini
c end of thermodynamics
c c2.3 ions and electrons
      DIMENSION brems1(kk),brems3(kk),ddz1(kk),ddz3(kk),degen(kk),
     &          effz(kk),eixch2(kk),fidash(kk),fz1(kk),fz3(kk),fzsq1(kk)
     &          ,fzsq3(kk),xlc(kk),xmieff(kk),xne(kk),xni(kk),omega1(kk)
      COMMON /comie / brems1,brems3,ddz1,ddz3,degen,degmax,degmin,effz,
     &                eixch2,fidash,fz1,fz3,fzsq1,fzsq3,xlc,xmieff,xne,
     &                xni,omega1,pmass
c end of ions and electrons
c c2.4 laser variables
      DIMENSION alpha1(kk),xaser1(kk),xl1(kk),xl3(kk)
      COMMON /comlas/ alpha1,elas1,xamda1,xaser1,xecri1,plas1,rabs1,
     &                rocri1,xl1,xl3,nabs1
c end of laser variables
c c2.5 thermonuclear reactions
      DIMENSION f1d(kk),f1h(kk),f1he3(kk),f1he4(kk),f1neu(kk),f1ntrl(kk)
     &          ,f1t(kk),f1x(kk),f1x1(kk),f1x2(kk),f1x3(kk),f1x4(kk),
     &          f3d(kk),f3h(kk),f3he3(kk),f3he4(kk),f3neu(kk),f3ntrl(kk)
     &          ,f3t(kk),f3x(kk),f3x1(kk),f3x2(kk),f3x3(kk),f3x4(kk),
     &          r1dd(kk),r1dhe3(kk),r1dt(kk),r3dd(kk),r3dhe3(kk),
     &          r3dt(kk),ye1(kk),ye3(kk),yi1(kk),yi3(kk)
      COMMON /comfus/ deuter,f1d,f1h,f1he3,f1he4,f1neu,f1ntrl,f1t,f1x,
     &                f1x1,f1x2,f1x3,f1x4,f3d,f3h,f3he3,f3he4,f3neu,
     &                f3ntrl,f3t,f3x,f3x1,f3x2,f3x3,f3x4,heliu3,heliu4,
     &                hydrog,xetral,xtrlms,pneut1,pneut3,r1dd,r1dhe3,
     &                r1dt,r3dd,r3dhe3,r3dt,rneut1,rneut3,tinucl,totneu,
     &                tritiu,xmass,xmass1,xmass2,xmass3,xmass4,xtra,
     &                xtra1,xtra2,xtra3,xtra4,xz,xz1,xz2,xz3,xz4,ye1,
     &                ye3,yi1,yi3,yield
c end of thermonuclear reactions
c c2.6 - energies
      COMMON /comen / eerror,en,einput,eloss,eefuse,eifuse,eindt1,
     &                eindt3,eneutr,pv,usqm,xionen,en0,eions,ecbbs,
     &                ecbfs,edcvs,erbfs,eatis,ezsts,erbbs,erfbs,edfbs
c c2.6 - end energies
c c2.7 physics control
      LOGICAL mlbrms,mlecon,mlfuse,mlicon,mlx,nlabs,nlbrms,nlburn,
     &        nlcri1,nldepo,nlecon,nlfuse,nlicon,nlmove,nlpfe,nlpfi,nlx,
     &        nlhf,nltnl
      COMMON /comcon/ mstep,ncase,ngeom,mlbrms,mlecon,mlfuse,mlicon,mlx,
     &                nlabs,nlbrms,nlburn,nlcri1,nldepo,nlecon,nlfuse,
     &                nlicon,nlmove,nlpfe,nlpfi,nlx,nlhf,nltnl,pmin,
     &                rhomin,temin,timin,umin
c end of physics control
c c3.1 numerical control parameters
      LOGICAL break,nlgoon,nlite
      COMMON /comnc / nceldt,ncondt,nit,nitmax,break,nlgoon,nlite,ak0,
     &                ak1,ak2,ak3,ak4,ak5,ak6,ak7,ak8,ak9,bneum,deltat,
     &                dt2,dt3,dt4,dtemax,dtimax,dtmax,dtmin,dumax,rdt2,
     &                rdt3,rdt4
c end of numerical control parameters
c c3.2 mesh and numerical methods
      DIMENSION ae(kk),ai(kk),be(kk),bi(kk),ce(kk),ci(kk),de(kk),di(kk),
     &          e(kk),f(kk),ge(kk),gi(kk),q2(kk),q4(kk),teite(kk),
     &          tiite(kk),uite(kk)
      COMMON /comnum/ mesh,nj,njm1,njp1,nl,nlm1,nlp1,ae,ai,be,bi,ce,ci,
     &                de,di,e,f,ge,gi,q2,q4,teite,tiite,uite
c end of mesh and numerical methods
c c4.1 - administrative variables
      LOGICAL nldump,nlemp
      DIMENSION piq(302)
      COMMON /comadm/ maxdim,maxrun,ndump,nrep,nldump,nlemp,piq,tstop
c c4.1 - end administrative variables
c c5.1 input-output control variables
      LOGICAL nlfilm,nlhcpy,nlprnt
      DIMENSION buf1(kk),buf2(kk),buf3(kk),buf4(kk),buf5(kk),buf6(kk)
      REAL nprnt,nproc
      COMMON /comout/ nfilm,nhdcpy,np1,np2,np3,nlfilm,nlhcpy,nlprnt,
     &        ztnext,nprnt,nproc,buf1,buf2,buf3,buf4,buf5,buf6,scp,
     &                scr,scrho,scte,scti,sctime,nt1,nt2
c end of input-output control variables
c 2. physics
      maxdim=kk
      maxrun=99999999
c 2.1 hydrodynamics - the length of block comhyd
      i21=13*maxdim+6
      CALL resetr(dm,i21,0.0)
c 2.2 thermodynamics - the length of block comth
      i22=18*maxdim+4
      CALL resetr(ddroe1,i22,0.0)
c 2.3 ions and electrons - the length of block comie
      i23=17*maxdim+3
      CALL resetr(brems1,i23,0.0)
c 2.4 laser variables - the length of block comlas
c     i24       = 3 * maxdim + 15 + 6
      i24=4*maxdim+6
      CALL resetr(alpha1,i24,0.0)
      nabs1=0
c 2.5 thermonuclear reactions - the length of block comfus
      i25=34*maxdim+29
      deuter=0.0
c 2.6 energies - the length of block comen
      i26=12
      eefuse=0.0
c 2.7 physics control - reset different types of variables separately
c integers
      i27=3
      mstep=0
c logicals
      i27=17
      mlbrms=.false.
c reals
      i27=5
      pmin=0.0
c 3. numerical scheme
c 3.1 numerical control parameters - reset different types of
c variables separately
c integers
      i31=4
      nceldt=0
c logicals
      i31=3
      break=.false.
c reals
      i31=19
      ak0=0.0
c 3.2 mesh and numerical methods - the length of block comnum
c integers
      i32=7
      mesh=0
c reals
      i32=17*maxdim
      CALL resetr(ae,i32,0.0)
c 4. house-keeping
c 4.1 administrative variables - reset different types of
c variables separately
c integers
      i41=4
c reset maximum mesh size
      maxdim=kk
      maxrun=99999999
c logicals
      i41=2
      nldump=.false.
c reals
      i41=maxdim+1
      CALL resetr(piq,i41,0.0)
c 5. input-output
c 5.1 input-output control variables - reset different types of
c variables separately
c integers
      i51=6
      nfilm=0
c logicals
      i51=3
      nlfilm=.false.
c reals
      i51=6*maxdim+6
      CALL resetr(buf1,i51,0.0)
c   hot electron common block
c integers
      i62=2
      igrp=0
c logicals
      i63=maxdim
      CALL resetl(lphi,i63,.true.)
c reals
c     maxgrp    = 20
      maxgrp=20
c     i61       = maxdim * (2 * maxgrp + 18) + 4 * maxgrp + 2
      i61=maxdim*(2*maxgrp+18)+3*maxgrp+2+21
      CALL resetr(deltae,i61,0.0)
      RETURN
      END


      SUBROUTINE clist(kgroup,kblock)
      PARAMETER(kk=1001)
c 5.2 print all common variables
c c1.1 basic system parameters
      LOGICAL nlend,nlres
      COMMON /combas/ altime,cptime,stime,ndiary,nin,nledge,nonlin,nout,
     &                nprint,npunch,nread,nrec,nresm,nrun,nstep,nlend,
     &                nlres
c end of basic system parameters
c c1.9 development and diagnostic parameters
      DIMENSION nadump(20),npdump(20),nvdump(20)
      LOGICAL nlched,nlhead(9),nlomt1(50),nlomt2(50),nlomt3(50),nlrept
      COMMON /comddp/ maxdum,mxdump,nadump,nclass,npdump,npoint,nsub,
     &                nvdump,nlched,nlhead,nlomt1,nlomt2,nlomt3,nlrept
c end of development and diagnostic parameters
c c2.1 hydrodynamics
      DIMENSION dm(kk),p3(kk),r1(kk),r3(kk),r5(kk),rho1(kk),rho3(kk),
     &          rho5(kk),u2(kk),u4(kk),v1(kk),v3(kk),v5(kk),lstat(kk)
      COMMON /comhyd/ dm,p3,pini,r1,r3,r5,rho1,rho3,rho5,rhoini,rhor,
     &                rini,time,u2,u4,uedge,v1,v3,v5
c end of hydrodynamics
c c2.2 thermodynamics
      DIMENSION ddroe1(kk),ddroe3(kk),ddroi1(kk),ddroi3(kk),ddte1(kk),
     &          ddte3(kk),ddti1(kk),ddti3(kk),xappae(kk),xappai(kk),
     &          pe1(kk),pe3(kk),pi1(kk),pi3(kk),te1(kk),te3(kk),ti1(kk),
     &          ti3(kk)
      COMMON /comth / ddroe1,ddroe3,ddroi1,ddroi3,ddte1,ddte3,ddti1,
     &                ddti3,gammae,gammai,xappae,xappai,pe1,pe3,pi1,pi3,
     &                te1,te3,teini,ti1,ti3,tiini
c end of thermodynamics
c c2.3 ions and electrons
      DIMENSION brems1(kk),brems3(kk),ddz1(kk),ddz3(kk),degen(kk),
     &          effz(kk),eixch2(kk),fidash(kk),fz1(kk),fz3(kk),fzsq1(kk)
     &          ,fzsq3(kk),xlc(kk),xmieff(kk),xne(kk),xni(kk),omega1(kk)
      COMMON /comie / brems1,brems3,ddz1,ddz3,degen,degmax,degmin,effz,
     &                eixch2,fidash,fz1,fz3,fzsq1,fzsq3,xlc,xmieff,xne,
     &                xni,omega1,pmass
c end of ions and electrons
c c2.4 laser variables
      DIMENSION alpha1(kk),xaser1(kk),xl1(kk),xl3(kk)
      COMMON /comlas/ alpha1,elas1,xamda1,xaser1,xecri1,plas1,rabs1,
     &                rocri1,xl1,xl3,nabs1
c end of laser variables
c c2.5 thermonuclear reactions
      DIMENSION f1d(kk),f1h(kk),f1he3(kk),f1he4(kk),f1neu(kk),f1ntrl(kk)
     &          ,f1t(kk),f1x(kk),f1x1(kk),f1x2(kk),f1x3(kk),f1x4(kk),
     &          f3d(kk),f3h(kk),f3he3(kk),f3he4(kk),f3neu(kk),f3ntrl(kk)
     &          ,f3t(kk),f3x(kk),f3x1(kk),f3x2(kk),f3x3(kk),f3x4(kk),
     &          r1dd(kk),r1dhe3(kk),r1dt(kk),r3dd(kk),r3dhe3(kk),
     &          r3dt(kk),ye1(kk),ye3(kk),yi1(kk),yi3(kk)
      COMMON /comfus/ deuter,f1d,f1h,f1he3,f1he4,f1neu,f1ntrl,f1t,f1x,
     &                f1x1,f1x2,f1x3,f1x4,f3d,f3h,f3he3,f3he4,f3neu,
     &                f3ntrl,f3t,f3x,f3x1,f3x2,f3x3,f3x4,heliu3,heliu4,
     &                hydrog,xetral,xtrlms,pneut1,pneut3,r1dd,r1dhe3,
     &                r1dt,r3dd,r3dhe3,r3dt,rneut1,rneut3,tinucl,totneu,
     &                tritiu,xmass,xmass1,xmass2,xmass3,xmass4,xtra,
     &                xtra1,xtra2,xtra3,xtra4,xz,xz1,xz2,xz3,xz4,ye1,
     &                ye3,yi1,yi3,yield
c end of thermonuclear reactions
c c2.6 - energies
      COMMON /comen / eerror,en,einput,eloss,eefuse,eifuse,eindt1,
     &                eindt3,eneutr,pv,usqm,xionen,en0,eions,ecbbs,
     &                ecbfs,edcvs,erbfs,eatis,ezsts,erbbs,erfbs,edfbs
c c2.6 - end energies
c c2.7 physics control
      LOGICAL mlbrms,mlecon,mlfuse,mlicon,mlx,nlabs,nlbrms,nlburn,
     &        nlcri1,nldepo,nlecon,nlfuse,nlicon,nlmove,nlpfe,nlpfi,nlx,
     &        nlhf,nltnl
      COMMON /comcon/ mstep,ncase,ngeom,mlbrms,mlecon,mlfuse,mlicon,mlx,
     &                nlabs,nlbrms,nlburn,nlcri1,nldepo,nlecon,nlfuse,
     &                nlicon,nlmove,nlpfe,nlpfi,nlx,nlhf,nltnl,pmin,
     &                rhomin,temin,timin,umin
c end of physics control
c c3.1 numerical control parameters
      LOGICAL break,nlgoon,nlite
      COMMON /comnc / nceldt,ncondt,nit,nitmax,break,nlgoon,nlite,ak0,
     &                ak1,ak2,ak3,ak4,ak5,ak6,ak7,ak8,ak9,bneum,deltat,
     &                dt2,dt3,dt4,dtemax,dtimax,dtmax,dtmin,dumax,rdt2,
     &                rdt3,rdt4
c end of numerical control parameters
c c3.2 mesh and numerical methods
      DIMENSION ae(kk),ai(kk),be(kk),bi(kk),ce(kk),ci(kk),de(kk),di(kk),
     &          e(kk),f(kk),ge(kk),gi(kk),q2(kk),q4(kk),teite(kk),
     &          tiite(kk),uite(kk)
      COMMON /comnum/ mesh,nj,njm1,njp1,nl,nlm1,nlp1,ae,ai,be,bi,ce,ci,
     &                de,di,e,f,ge,gi,q2,q4,teite,tiite,uite
c end of mesh and numerical methods
c c4.1 - administrative variables
      LOGICAL nldump,nlemp
      DIMENSION piq(302)
      COMMON /comadm/ maxdim,maxrun,ndump,nrep,nldump,nlemp,piq,tstop
c c4.1 - end administrative variables
c c5.1 input-output control variables
      LOGICAL nlfilm,nlhcpy,nlprnt
      DIMENSION buf1(kk),buf2(kk),buf3(kk),buf4(kk),buf5(kk),buf6(kk)
      REAL nprnt,nproc
      COMMON /comout/ nfilm,nhdcpy,np1,np2,np3,nlfilm,nlhcpy,nlprnt,
     &        ztnext,nprnt,nproc,buf1,buf2,buf3,buf4,buf5,buf6,scp,
     &                scr,scrho,scte,scti,sctime,nt1,nt2
c end of input-output control variables
c 1. general olympus data
c select group to be printed
      IF(kgroup.ne.0)THEN
       IF(kgroup.eq.1)THEN
       ELSEIF(kgroup.eq.2)THEN
        GOTO 300
       ELSEIF(kgroup.eq.3)THEN
        GOTO 1000
       ELSEIF(kgroup.eq.4)THEN
        GOTO 1200
       ELSEIF(kgroup.eq.5)THEN
        GOTO 1300
       ELSE
        CALL gotoer(' MEDUSA: clist: goto error, kgroup is',kgroup)
       ENDIF
      ENDIF
      IF(kblock.ne.0)THEN
       IF(kblock.eq.1)THEN
       ELSEIF(kblock.eq.2.or.kblock.eq.3.or.kblock.eq.4.or.
     &        kblock.eq.5.or.kblock.eq.6.or.kblock.eq.7.or.kblock.eq.8)
     &        THEN
        GOTO 200
       ELSEIF(kblock.eq.9)THEN
        GOTO 100
       ELSE
        CALL gotoer(' MEDUSA: clist: goto error, kblock is',kblock)
       ENDIF
      ENDIF
c 1.1 basic system parameters
      CALL page
      CALL blines(2)
      CALL mesage('      c1.1      combas')
      CALL blines(2)
      CALL rvar('altime  ',altime)
      CALL rvar('cptime  ',cptime)
      CALL ivar('ndiary  ',ndiary)
      CALL ivar('nin     ',nin)
      CALL ivar('nledge  ',nledge)
      CALL lvar('nlend   ',nlend)
      CALL lvar('nlres   ',nlres)
      CALL ivar('nonlin  ',nonlin)
      CALL ivar('nout    ',nout)
      CALL ivar('nprint  ',nprint)
      CALL ivar('npunch  ',npunch)
      CALL ivar('nread   ',nread)
      CALL ivar('nrec    ',nrec)
      CALL ivar('nresm  ',nresm)
      CALL ivar('nrun    ',nrun)
      CALL ivar('nstep   ',nstep)
      CALL rvar('stime   ',stime)
      IF(kblock.ne.0)RETURN
c 1.9 diagnostics and program development
 100  CALL page
      CALL blines(2)
      CALL mesage('      c1.9      comddp')
      CALL blines(2)
      CALL ivar('maxdum  ',maxdum)
      CALL ivar('mxdump  ',mxdump)
      CALL ivar('nclass  ',nclass)
      CALL lvar('nlched  ',nlched)
      CALL lvar('nlrept  ',nlrept)
      CALL ivar('npoint  ',npoint)
      CALL ivar('nsub    ',nsub)
 200  IF(kgroup.ne.0)RETURN
c 2. physics
 300  IF(kblock.ne.0)THEN
       IF(kblock.eq.1)THEN
       ELSEIF(kblock.eq.2)THEN
        GOTO 400
       ELSEIF(kblock.eq.3)THEN
        GOTO 500
       ELSEIF(kblock.eq.4)THEN
        GOTO 600
       ELSEIF(kblock.eq.5)THEN
        GOTO 700
       ELSEIF(kblock.eq.6)THEN
        GOTO 800
       ELSEIF(kblock.eq.7)THEN
        GOTO 900
       ELSE
        CALL gotoer(' MEDUSA: clist: goto error, kblock is',kblock)
       ENDIF
      ENDIF
c 2.1 hydrodynamics
      CALL page
      CALL blines(2)
      CALL mesage('      c2.1      comhyd')
      CALL blines(2)
      CALL rvar('rhoini  ',rhoini)
      CALL rvar('rhor    ',rhor)
      CALL rvar('rini    ',rini)
      CALL rvar('time    ',time)
      CALL rvar('uedge   ',uedge)
      IF(kblock.ne.0)RETURN
c 2.2 thermodynamics
 400  CALL page
      CALL blines(2)
      CALL mesage('      c2.2      comth')
      CALL blines(2)
      CALL rvar('gammae  ',gammae)
      CALL rvar('gammai  ',gammai)
      CALL rvar('teini   ',teini)
      CALL rvar('tiini   ',tiini)
      IF(kblock.ne.0)RETURN
c 2.3 ions and electrons
 500  CALL page
      CALL blines(2)
      CALL mesage('      c2.3      comie')
      CALL blines(2)
      CALL rvar('degmax  ',degmax)
      CALL rvar('degmin  ',degmin)
      CALL rvar('pmass   ',pmass)
      IF(kblock.ne.0)RETURN
c 2.4 laser variables
 600  CALL page
      CALL blines(2)
      CALL mesage('      c2.4      comlas')
      CALL blines(2)
      CALL rvar('elas1   ',elas1)
      CALL rvar('xamda1  ',xamda1)
      CALL rvar('plas1   ',plas1)
      CALL rvar('rabs1   ',rabs1)
      CALL rvar('rocri1  ',rocri1)
      CALL ivar('nabs1   ',nabs1)
      IF(kblock.ne.0)RETURN
c 2.5 thermonuclear reactions
 700  CALL page
      CALL blines(2)
      CALL mesage('      c2.5      comfus')
      CALL blines(2)
      CALL rvar('deuter  ',deuter)
      CALL rvar('heliu3  ',heliu3)
      CALL rvar('heliu4  ',heliu4)
      CALL rvar('hydrog  ',hydrog)
      CALL rvar('xetral  ',xetral)
      CALL rvar('xtrlms  ',xtrlms)
      CALL rvar('pneut1  ',pneut1)
      CALL rvar('pneut3  ',pneut3)
      CALL rvar('rneut1  ',rneut1)
      CALL rvar('rneut3  ',rneut3)
      CALL rvar('tinucl  ',tinucl)
      CALL rvar('totneu  ',totneu)
      CALL rvar('tritiu  ',tritiu)
      CALL rvar('xmass   ',xmass)
      CALL rvar('xtra    ',xtra)
      CALL rvar('xz      ',xz)
      CALL rvar('xmass1  ',xmass1)
      CALL rvar('xmass2  ',xmass2)
      CALL rvar('xmass3  ',xmass3)
      CALL rvar('xmass4  ',xmass4)
      CALL rvar('xz1     ',xz1)
      CALL rvar('xz2     ',xz2)
      CALL rvar('xz3     ',xz3)
      CALL rvar('xz4     ',xz4)
      CALL rvar('xtra1   ',xtra1)
      CALL rvar('xtra2   ',xtra2)
      CALL rvar('xtra3   ',xtra3)
      CALL rvar('xtra4   ',xtra4)
      CALL rvar('yield   ',yield)
      IF(kblock.ne.0)RETURN
c 2.6 energies
 800  CALL page
      CALL blines(2)
      CALL mesage('      c2.6      comen')
      CALL blines(2)
      CALL rvar('eefuse  ',eefuse)
      CALL rvar('eifuse  ',eifuse)
      CALL rvar('eerror  ',eerror)
      CALL rvar('eindt1  ',eindt1)
      CALL rvar('eindt3  ',eindt3)
      CALL rvar('einput  ',einput)
      CALL rvar('eloss   ',eloss)
      CALL rvar('eions     ',eions)
      CALL rvar('ecbbs     ',ecbbs)
      CALL rvar('ecbfs     ',ecbfs)
      CALL rvar('edcvs     ',edcvs)
      CALL rvar('erbfs     ',erbfs)
      CALL rvar('erbbs     ',erbbs)
      CALL rvar('erfbs     ',erfbs)
      CALL rvar('edfbs     ',edfbs)
      CALL rvar('en      ',en)
      CALL rvar('eneutr  ',eneutr)
      CALL rvar('pv      ',pv)
      CALL rvar('usqm    ',usqm)
      IF(kblock.ne.0)RETURN
c 2.7 physics control
 900  CALL page
      CALL blines(2)
      CALL mesage('      c2.7      comcon')
      CALL blines(2)
      CALL rvar('pmin    ',pmin)
      CALL rvar('rhomin  ',rhomin)
      CALL rvar('temin   ',temin)
      CALL rvar('timin   ',timin)
      CALL rvar('umin    ',umin)
      CALL ivar('mstep   ',mstep)
      CALL ivar('ncase   ',ncase)
      CALL ivar('ngeom   ',ngeom)
      CALL lvar('mlbrms  ',mlbrms)
      CALL lvar('mlecon  ',mlecon)
      CALL lvar('mlfuse  ',mlfuse)
      CALL lvar('mlicon  ',mlicon)
      CALL lvar('mlx     ',mlx)
      CALL lvar('nlabs   ',nlabs)
      CALL lvar('nlbrms  ',nlbrms)
      CALL lvar('nlburn  ',nlburn)
      CALL lvar('nlcri1  ',nlcri1)
      CALL lvar('nldepo  ',nldepo)
      CALL lvar('nlecon  ',nlecon)
      CALL lvar('nlfuse  ',nlfuse)
      CALL lvar('nlicon  ',nlicon)
      CALL lvar('nlmove  ',nlmove)
      CALL lvar('nlpfe   ',nlpfe)
      CALL lvar('nlpfi   ',nlpfi)
      CALL lvar('nlx     ',nlx)
      IF(kgroup.ne.0)RETURN
c 3. numerical scheme
 1000 IF(kblock.ne.0)THEN
       IF(kblock.eq.1)THEN
       ELSEIF(kblock.eq.2)THEN
        GOTO 1100
       ELSE
        CALL gotoer(' MEDUSA: clist: goto error, kblock is',kblock)
       ENDIF
      ENDIF
c 3.1 numerical control parameters
      CALL page
      CALL blines(2)
      CALL mesage('      c3.1      comnc')
      CALL blines(2)
      CALL rvar('ak0     ',ak0)
      CALL rvar('ak1     ',ak1)
      CALL rvar('ak2     ',ak2)
      CALL rvar('ak3     ',ak3)
      CALL rvar('ak4     ',ak4)
      CALL rvar('ak5     ',ak5)
      CALL rvar('bneum   ',bneum)
      CALL rvar('delta t ',deltat)
      CALL rvar('dt2     ',dt2)
      CALL rvar('dt3     ',dt3)
      CALL rvar('dt4     ',dt4)
      CALL rvar('dtemax  ',dtemax)
      CALL rvar('dtimax  ',dtimax)
      CALL rvar('dtmax   ',dtmax)
      CALL rvar('dtmin   ',dtmin)
      CALL rvar('dumax   ',dumax)
      CALL rvar('rdt2    ',rdt2)
      CALL rvar('rdt3    ',rdt3)
      CALL rvar('rdt4    ',rdt4)
      CALL ivar('nceldt  ',nceldt)
      CALL ivar('ncondt  ',ncondt)
      CALL ivar('nit     ',nit)
      CALL ivar('nitmax  ',nitmax)
      CALL lvar('break   ',break)
      CALL lvar('nlgoon  ',nlgoon)
      CALL lvar('nlite   ',nlite)
      IF(kblock.ne.0)RETURN
c 3.2 mesh and numerical methods
 1100 CALL page
      CALL blines(2)
      CALL mesage('      c3.2      comnum')
      CALL blines(2)
      CALL ivar('mesh    ',mesh)
      CALL ivar('nj      ',nj)
      CALL ivar('njm1    ',njm1)
      CALL ivar('njp1    ',njp1)
      CALL ivar('nl      ',nl)
      CALL ivar('nlm1    ',nlm1)
      CALL ivar('nlp1    ',nlp1)
      IF(kgroup.ne.0)RETURN
c 4. house-keeping
c 4.1 administrative variables
 1200 CALL page
      CALL blines(2)
      CALL mesage('      c4.1      comadm')
      CALL blines(2)
      CALL rvar('tstop   ',tstop)
      CALL ivar('maxdim  ',maxdim)
      CALL ivar('maxrun  ',maxrun)
      CALL ivar('ndump   ',ndump)
      CALL lvar('nldump  ',nldump)
      CALL lvar('nlemp   ',nlemp)
      CALL ivar('nrep    ',nrep)
      IF(kgroup.ne.0)RETURN
c 5. input-output
c 5.1 input-output control variables
 1300 CALL page
      CALL blines(2)
      CALL mesage('      c5.1      comout')
      CALL blines(2)
      CALL rvar('scp     ',scp)
      CALL rvar('scr     ',scr)
      CALL rvar('scrho   ',scrho)
      CALL rvar('scte    ',scte)
      CALL rvar('scti    ',scti)
      CALL rvar('sctime  ',sctime)
      CALL ivar('nfilm   ',nfilm)
      CALL ivar('nhdcpy  ',nhdcpy)
      CALL ivar('np1     ',np1)
      CALL ivar('np3     ',np3)
      CALL lvar('nlfilm  ',nlfilm)
      CALL lvar('nlhcpy  ',nlhcpy)
      CALL lvar('nlprnt  ',nlprnt)
      RETURN
      END


      SUBROUTINE copyi(ka1,k1,ka2,k2,kdim)
c u.24 copy one integer array into another
      DIMENSION ka1(kdim),ka2(kdim)
      DO 100 j=1,kdim
       i1=k1+j-1
       i2=k2+j-1
       ka2(i2)=ka1(i1)
 100  CONTINUE
      RETURN
      END


      SUBROUTINE copyr(pa1,k1,pa2,k2,kdim)
c u.23 copy one real array into another
      DIMENSION pa1(kdim),pa2(kdim)
      DO 100 j=1,kdim
       i1=k1+j-1
       i2=k2+j-1
       pa2(i2)=pa1(i1)
 100  CONTINUE
      RETURN
      END


      SUBROUTINE cotrol
c 0.3 control the run
c c1.1 basic system parameters
      LOGICAL nlend,nlres
      COMMON /combas/ altime,cptime,stime,ndiary,nin,nledge,nonlin,nout,
     &                nprint,npunch,nread,nrec,nresm,nrun,nstep,nlend,
     &                nlres
c end of basic system parameters
c c1.9 development and diagnostic parameters
c end of development and diagnostic parameters
      DATA iclass,isub/0,3/
c 1. prologue
c 1.7 pick up record, modify required parameters
c label the new run or the continuation run
      CALL labrun
      CALL expert(iclass,isub,8)
c clear variables and arrays
      CALL clear
      CALL expert(iclass,isub,9)
c 1.3 set default values
      CALL preset
      CALL expert(iclass,isub,3)
c read data file
      CALL data
      CALL expert(iclass,isub,11)
c 1.5 set auxiliary values
      CALL auxval
      CALL expert(iclass,isub,5)
c
c        a. new run
c        define physical initial conditions
       CALL inital
       CALL expert(iclass,isub,6)
c        preliminary operations
c        start  the run
       CALL start
       CALL expert(iclass,isub,13)
      IF(nlres)THEN
c        resum from previous run
       CALL resumerun
      ENDIF
c
c initial output
      CALL output(1)
      CALL expert(iclass,isub,14)
c 2. calculation
c 2.1 step on the calculation
 100  CALL stepon
      CALL expert(iclass,isub,15)
c 3. output
c 3.1 periodic production of output
      CALL output(2)
      CALL expert(iclass,isub,16)
c 4. epilogue
c 4.1 test for completion of run
      CALL tesend
      CALL expert(iclass,isub,17)
      IF(.not.nlend)GOTO 100
c final output
c     call output(3)
c     call expert(iclass,isub,18)
c 4.2 terminate the run
      CALL endrun
      RETURN
      END


      SUBROUTINE coulog
      PARAMETER(kk=1001)
c 2.8 evaluates the coulomb logarithm
c c1.9 development and diagnostic parameters
      DIMENSION nadump(20),npdump(20),nvdump(20)
      LOGICAL nlched,nlhead(9),nlomt1(50),nlomt2(50),nlomt3(50),nlrept
      COMMON /comddp/ maxdum,mxdump,nadump,nclass,npdump,npoint,nsub,
     &                nvdump,nlched,nlhead,nlomt1,nlomt2,nlomt3,nlrept
c end of development and diagnostic parameters
c c2.2 thermodynamics
      DIMENSION ddroe1(kk),ddroe3(kk),ddroi1(kk),ddroi3(kk),ddte1(kk),
     &          ddte3(kk),ddti1(kk),ddti3(kk),xappae(kk),xappai(kk),
     &          pe1(kk),pe3(kk),pi1(kk),pi3(kk),te1(kk),te3(kk),ti1(kk),
     &          ti3(kk)
      COMMON /comth / ddroe1,ddroe3,ddroi1,ddroi3,ddte1,ddte3,ddti1,
     &                ddti3,gammae,gammai,xappae,xappai,pe1,pe3,pi1,pi3,
     &                te1,te3,teini,ti1,ti3,tiini
c end of thermodynamics
c c2.3 ions and electrons
      DIMENSION brems1(kk),brems3(kk),ddz1(kk),ddz3(kk),degen(kk),
     &          effz(kk),eixch2(kk),fidash(kk),fz1(kk),fz3(kk),fzsq1(kk)
     &          ,fzsq3(kk),xlc(kk),xmieff(kk),xne(kk),xni(kk),omega1(kk)
      COMMON /comie / brems1,brems3,ddz1,ddz3,degen,degmax,degmin,effz,
     &                eixch2,fidash,fz1,fz3,fzsq1,fzsq3,xlc,xmieff,xne,
     &                xni,omega1,pmass
c end of ions and electrons
c c3.2 mesh and numerical methods
      DIMENSION ae(kk),ai(kk),be(kk),bi(kk),ce(kk),ci(kk),de(kk),di(kk),
     &          e(kk),f(kk),ge(kk),gi(kk),q2(kk),q4(kk),teite(kk),
     &          tiite(kk),uite(kk)
      COMMON /comnum/ mesh,nj,njm1,njp1,nl,nlm1,nlp1,ae,ai,be,bi,ce,ci,
     &                de,di,e,f,ge,gi,q2,q4,teite,tiite,uite
c end of mesh and numerical methods
c c4.1 - administrative variables
      LOGICAL nldump,nlemp
      DIMENSION piq(302)
      COMMON /comadm/ maxdim,maxrun,ndump,nrep,nldump,nlemp,piq,tstop
c c4.1 - end administrative variables

c c2.7 physics control
c end of physics control

      DATA iclass,isub/2,8/
      IF(nlomt2(isub))RETURN
c 1. evaluate the coulomb logarithm (eq.16)
      zc=1.2427*1.0E7*(1.0+piq(20))
      zloge=2.3026
      DO 100 l=1,nl
       zl=te3(l)*te3(l)*te3(l)/xne(l)
       zlam=zc*sqrt(zl)/fz3(l)
       zlam2=6.9617E9*te3(l)/sqrt(xne(l))
       IF(zlam.gt.zlam2)zlam=zlam2
       IF(zlam.le.10.)zlam=10.
       xlc(l)=zloge*log10(zlam)
c ad1 allan: set coulomb log (in heat conduction and e-i exchange) to 5
c        xlc(l)=5.
 100  CONTINUE
      zl=te3(nlp1)*te3(nlp1)*te3(nlp1)/xne(nl)
      zlam=zc*sqrt(zl)/fz3(nl)
      IF(zlam.le.10.)zlam=10.
      xlc(nlp1)=zloge*log10(zlam)
c ad1 allan: set coulomb log at mesh nlp1
c     xlc(nlp1)=5.
c lc has been evaluated in the boundary cell l = nlp1 since
c case 2 requires the conductivity in this cell
      RETURN
      END


      SUBROUTINE cvar(kname,kvalue)
c u.6 print name and value of character string
c c1.1 basic system parameters
      LOGICAL nlend,nlres
      COMMON /combas/ altime,cptime,stime,ndiary,nin,nledge,nonlin,nout,
     &                nprint,npunch,nread,nrec,nresm,nrun,nstep,nlend,
     &                nlres
c end of basic system parameters
      CHARACTER kname*(*),kvalue*(*)
      WRITE(nout,10100)kname,kvalue
      RETURN
10100 FORMAT(4x,a,' = ','''',a,'''')
      END


      SUBROUTINE cverge
      PARAMETER(kk=1001)
c 2.21 examines the convergence of ti, te and u
c c1.9 development and diagnostic parameters
      DIMENSION nadump(20),npdump(20),nvdump(20)
      LOGICAL nlched,nlhead(9),nlomt1(50),nlomt2(50),nlomt3(50),nlrept
      COMMON /comddp/ maxdum,mxdump,nadump,nclass,npdump,npoint,nsub,
     &                nvdump,nlched,nlhead,nlomt1,nlomt2,nlomt3,nlrept
c end of development and diagnostic parameters
c c2.1 hydrodynamics
      DIMENSION dm(kk),p3(kk),r1(kk),r3(kk),r5(kk),rho1(kk),rho3(kk),
     &          rho5(kk),u2(kk),u4(kk),v1(kk),v3(kk),v5(kk),lstat(kk)
      COMMON /comhyd/ dm,p3,pini,r1,r3,r5,rho1,rho3,rho5,rhoini,rhor,
     &                rini,time,u2,u4,uedge,v1,v3,v5
c end of hydrodynamics
c c2.2 thermodynamics
      DIMENSION ddroe1(kk),ddroe3(kk),ddroi1(kk),ddroi3(kk),ddte1(kk),
     &          ddte3(kk),ddti1(kk),ddti3(kk),xappae(kk),xappai(kk),
     &          pe1(kk),pe3(kk),pi1(kk),pi3(kk),te1(kk),te3(kk),ti1(kk),
     &          ti3(kk)
      COMMON /comth / ddroe1,ddroe3,ddroi1,ddroi3,ddte1,ddte3,ddti1,
     &                ddti3,gammae,gammai,xappae,xappai,pe1,pe3,pi1,pi3,
     &                te1,te3,teini,ti1,ti3,tiini
c end of thermodynamics
c c2.3 ions and electrons
      DIMENSION brems1(kk),brems3(kk),ddz1(kk),ddz3(kk),degen(kk),
     &          effz(kk),eixch2(kk),fidash(kk),fz1(kk),fz3(kk),fzsq1(kk)
     &          ,fzsq3(kk),xlc(kk),xmieff(kk),xne(kk),xni(kk),omega1(kk)
      COMMON /comie / brems1,brems3,ddz1,ddz3,degen,degmax,degmin,effz,
     &                eixch2,fidash,fz1,fz3,fzsq1,fzsq3,xlc,xmieff,xne,
     &                xni,omega1,pmass
c end of ions and electrons
c c2.7 physics control
      LOGICAL mlbrms,mlecon,mlfuse,mlicon,mlx,nlabs,nlbrms,nlburn,
     &        nlcri1,nldepo,nlecon,nlfuse,nlicon,nlmove,nlpfe,nlpfi,nlx,
     &        nlhf,nltnl
      COMMON /comcon/ mstep,ncase,ngeom,mlbrms,mlecon,mlfuse,mlicon,mlx,
     &                nlabs,nlbrms,nlburn,nlcri1,nldepo,nlecon,nlfuse,
     &                nlicon,nlmove,nlpfe,nlpfi,nlx,nlhf,nltnl,pmin,
     &                rhomin,temin,timin,umin
c end of physics control
c c3.1 numerical control parameters
      LOGICAL break,nlgoon,nlite
      COMMON /comnc / nceldt,ncondt,nit,nitmax,break,nlgoon,nlite,ak0,
     &                ak1,ak2,ak3,ak4,ak5,ak6,ak7,ak8,ak9,bneum,deltat,
     &                dt2,dt3,dt4,dtemax,dtimax,dtmax,dtmin,dumax,rdt2,
     &                rdt3,rdt4
c end of numerical control parameters
c c3.2 mesh and numerical methods
      DIMENSION ae(kk),ai(kk),be(kk),bi(kk),ce(kk),ci(kk),de(kk),di(kk),
     &          e(kk),f(kk),ge(kk),gi(kk),q2(kk),q4(kk),teite(kk),
     &          tiite(kk),uite(kk)
      COMMON /comnum/ mesh,nj,njm1,njp1,nl,nlm1,nlp1,ae,ai,be,bi,ce,ci,
     &                de,di,e,f,ge,gi,q2,q4,teite,tiite,uite
c end of mesh and numerical methods
      DATA iclass,isub/2,21/
      IF(nlomt2(isub))RETURN
c 1. find max. deviations in ti,te and u (eq.78)
c reset switch before examining convergence
      nlgoon=.false.
      break=.false.
      zcu=0.0
      zcti=0.0
      zcte=0.0
      DO 100 l=1,nl
c notice we add 1 to l to avoid r = 0 and to include the boundary
       j=l+1
       IF(ti3(l).le.0.0)GOTO 200
       IF(te3(l).le.0.0)GOTO 200
       zdti=abs(ti3(l)-tiite(l))/(ti3(l)+tiite(l)+timin)
       zdte=abs(te3(l)-teite(l))/(te3(l)+teite(l)+temin)
c for very strong degeneracy we ignore the convergence check on te
       IF((.not.nlpfe).and.(degen(l).lt.0.01))zdte=0.0
c if the displacement caused by u is below zumin * dt or if the
c deviation in u*dt caused by iterations is too small then we ignore
c the convergence check
       zdr3=r3(j)-r3(j-1)
       zumin=zdr3*rdt3*1.0E-2
       zditu=abs(u4(j)-uite(j))*dt3
       zratio=zditu/zdr3
       IF((zratio.ge.1.0E-4).and.(abs(u4(j)).ge.zumin))THEN
        zdu=abs(u4(j)-uite(j))/(abs(u4(j))+abs(uite(j))+umin)
        IF(zdu.ge.zcu)THEN
         zcu=zdu
         ij=j
        ENDIF
       ENDIF
       IF(zdti.ge.zcti)THEN
        zcti=zdti
        il3=l
       ENDIF
       IF(zdte.ge.zcte)THEN
        zcte=zdte
        il4=l
       ENDIF
 100  CONTINUE
c 2. is convergence achieved ?
c include the factor 2 resulting from averaging
      zcti=2.0*zcti
      zcte=2.0*zcte
      zcu=2.0*zcu
c if no convergence find out where it happened
      IF(zcti.gt.dtimax)THEN
       ncondt=3
       nceldt=il3
      ELSEIF(zcte.gt.dtemax)THEN
       ncondt=4
       nceldt=il4
      ELSEIF(zcu.le.dumax)THEN
c the iteration has converged
       nlgoon=.true.
       RETURN
      ELSE
       ncondt=1
       nceldt=ij
      ENDIF
c 3. no convergence. another iteration ?
c     we check the number of iterations made
      IF(nit.ge.nitmax)GOTO 300
c if the number of iterations exceeds the maximum value permitted,
c then the timestep hasn't converged. other wise  shift the
c last found values into the arrays labelled -ite.
c and do the  next iteration
      CALL copyr(ti3,1,tiite,1,nlp1)
      CALL copyr(te3,1,teite,1,nlp1)
      CALL copyr(u4,1,uite,1,nlp1)
      RETURN
c 3.1 negative temperatures
c    prepare to call subrout revers
 200  break=.true.
      nlgoon=.false.
      nceldt=j-1
      IF(ti3(nceldt).le.0.0)ncondt=3
      IF(te3(nceldt).le.0.0)ncondt=4
      RETURN
 300  CALL mesage(' ****** warning iterations did not converge *****')
      CALL ivar('ncondt  ',ncondt)
      CALL ivar('nceldt  ',nceldt)
c    prepare for the next timestep
      break=.false.
      nlgoon=.true.
c     stop
      RETURN
      END


      SUBROUTINE daytim
c u.13 print date and time
c c1.1 basic system parameters
      LOGICAL nlend,nlres
      COMMON /combas/ altime,cptime,stime,ndiary,nin,nledge,nonlin,nout,
     &                nprint,npunch,nread,nrec,nresm,nrun,nstep,nlend,
     &                nlres
c end of basic system parameters
c for the cray at ral
      WRITE(nout,10300)
      WRITE(nout,10100)
      WRITE(nout,10200)
      WRITE(nout,10300)
c     call clock(time)
c     call date(day)
c     write(nout, 10) time,day
c  10 format(5x,' time ', a8,' date ',a8)
      RETURN
10100 FORMAT(45x,'timing routines are turned off')
10200 FORMAT(49x,'to preserve portability')
10300 FORMAT(45x,'------------------------------')
      END


      SUBROUTINE dumcom(kclass,ksub,kpoint)
c u.27 dump selected common blocks
c c1.1 basic system parameters
      LOGICAL nlend,nlres
      COMMON /combas/ altime,cptime,stime,ndiary,nin,nledge,nonlin,nout,
     &                nprint,npunch,nread,nrec,nresm,nrun,nstep,nlend,
     &                nlres
c end of basic system parameters
c c1.9 development and diagnostic parameters
      DIMENSION nadump(20),npdump(20),nvdump(20)
      LOGICAL nlched,nlhead(9),nlomt1(50),nlomt2(50),nlomt3(50),nlrept
      COMMON /comddp/ maxdum,mxdump,nadump,nclass,npdump,npoint,nsub,
     &                nvdump,nlched,nlhead,nlomt1,nlomt2,nlomt3,nlrept
c end of development and diagnostic parameters
      LOGICAL ilrept
      DIMENSION idigit(8)
      DATA idmax/8/,ilrept/.true./
c 1. initialize and print step number
      icode=10000*kclass+100*ksub+kpoint
      IF(icode.eq.20101)THEN
       CALL blines(1)
       CALL ivar('step    ',nstep)
       CALL mesage('**********************')
      ENDIF
      IF(kclass.ne.0)THEN
       IF(.not.nlhead(kclass))GOTO 100
      ELSEIF(.not.(nlched))THEN
       GOTO 100
      ENDIF
      CALL repthd(kclass,ksub,kpoint)
      ilrept=.false.
c 2. scan over list
 100  DO 400 j=1,mxdump
       IF(npdump(j).ne.icode)GOTO 400
c 3. dumping point found
c print heading only once
       IF(ilrept)CALL repthd(kclass,ksub,kpoint)
       ilrept=.false.
c 4. are variables to be dumped?
       i=1
       id=nvdump(j)
       IF(id.ne.0)THEN
        IF(id.ne.100)GOTO 200
        CALL clist(0,0)
       ENDIF
c 5. are arrays to be dumped
 150   i=2
       id=nadump(j)
       IF(id.eq.0)GOTO 400
       IF(id.eq.100)THEN
        CALL arrays(0,0)
        GOTO 400
       ENDIF
c 6. disentangle code
 200   CALL reseti(idigit,idmax,0)
       DO 250 jd=1,idmax
        idiv=id/10
        idigit(jd)=id-idiv*10
        in=jd
        IF(idiv.eq.0)GOTO 300
        id=idiv
 250   CONTINUE
c 7. issue calls
c make *in* even
 300   in=2*(in/2)
       IF(in.ne.0)THEN
        DO 320 jd=1,in,2
         ij=in-jd+1
         ig=idigit(ij)
         ib=idigit(ij-1)
         IF(ig.ne.0)THEN
          IF(i.eq.1)CALL clist(ig,ib)
          IF(i.eq.2)CALL arrays(ig,ib)
         ENDIF
 320    CONTINUE
        IF(i.eq.1)GOTO 150
       ENDIF
c 8. next entry in list
 400  CONTINUE
      RETURN
      END


      SUBROUTINE endrun
      PARAMETER(kk=1001)
c 4.2 terminate the run
c c1.1 basic system parameters
      LOGICAL nlend,nlres
      COMMON /combas/ altime,cptime,stime,ndiary,nin,nledge,nonlin,nout,
     &                nprint,npunch,nread,nrec,nresm,nrun,nstep,nlend,
     &                nlres
c end of basic system parameters
c c2.4 laser variables
      DIMENSION alpha1(kk),xaser1(kk),xl1(kk),xl3(kk)
      COMMON /comlas/ alpha1,elas1,xamda1,xaser1,xecri1,plas1,rabs1,
     &                rocri1,xl1,xl3,nabs1
c end of laser variables
c c2.5 thermonuclear reactions
c end of thermonuclear reactions
c c3.2 mesh and numerical methods
c end of mesh and numerical methods
      CALL output(3)
      CALL page
      CALL blines(2)
c     call daytim
 100  CALL runtim
      IF(nout.ne.ndiary)CALL blines(2)
      IF(.not.(nlend))CALL mesage(
     &                '...........a b n o r m a l   e x i t   .........'
     &                )
c     call mesage('***********  e n d   o f   j o b   *************')
      IF(nout.ne.ndiary)THEN
c     we repeat these messages and send them to the diary
       nout=ndiary
       zplas1=xaser1(2)
       zelas1=xaser1(7)
       GOTO 100
      ENDIF
      STOP
c     the output channel ndiary is set by basic (olympus)
c     write(nout, 9901) mesh, nstep, yield
10100 FORMAT(4x,'mesh =',i4,4x,'timestep =',i4,4x,'yield =',e10.3)
c     write(nout, 9902) xamda1, zplas1, zelas1
10200 FORMAT(4x,'l =',e10.3,4x,'p =',e10.3,4x,'e =',e10.3)
      END


      SUBROUTINE energy
      PARAMETER(kk=1001,underf=1.E-30)
c 2.22 kinetic and thermal energy. energy input and gain.
c c1.9 development and diagnostic parameters
      DIMENSION nadump(20),npdump(20),nvdump(20)
      LOGICAL nlched,nlhead(9),nlomt1(50),nlomt2(50),nlomt3(50),nlrept
      COMMON /comddp/ maxdum,mxdump,nadump,nclass,npdump,npoint,nsub,
     &                nvdump,nlched,nlhead,nlomt1,nlomt2,nlomt3,nlrept
c end of development and diagnostic parameters
c c2.1 hydrodynamics
      DIMENSION dm(kk),p3(kk),r1(kk),r3(kk),r5(kk),rho1(kk),rho3(kk),
     &          rho5(kk),u2(kk),u4(kk),v1(kk),v3(kk),v5(kk),lstat(kk)
      COMMON /comhyd/ dm,p3,pini,r1,r3,r5,rho1,rho3,rho5,rhoini,rhor,
     &                rini,time,u2,u4,uedge,v1,v3,v5
c end of hydrodynamics
c c2.2 thermodynamics
      DIMENSION ddroe1(kk),ddroe3(kk),ddroi1(kk),ddroi3(kk),ddte1(kk),
     &          ddte3(kk),ddti1(kk),ddti3(kk),xappae(kk),xappai(kk),
     &          pe1(kk),pe3(kk),pi1(kk),pi3(kk),te1(kk),te3(kk),ti1(kk),
     &          ti3(kk)
      COMMON /comth / ddroe1,ddroe3,ddroi1,ddroi3,ddte1,ddte3,ddti1,
     &                ddti3,gammae,gammai,xappae,xappai,pe1,pe3,pi1,pi3,
     &                te1,te3,teini,ti1,ti3,tiini
c end of thermodynamics
c c2.3 ions and electrons
      DIMENSION brems1(kk),brems3(kk),ddz1(kk),ddz3(kk),degen(kk),
     &          effz(kk),eixch2(kk),fidash(kk),fz1(kk),fz3(kk),fzsq1(kk)
     &          ,fzsq3(kk),xlc(kk),xmieff(kk),xne(kk),xni(kk),omega1(kk)
      COMMON /comie / brems1,brems3,ddz1,ddz3,degen,degmax,degmin,effz,
     &                eixch2,fidash,fz1,fz3,fzsq1,fzsq3,xlc,xmieff,xne,
     &                xni,omega1,pmass
c end of ions and electrons
c c2.4 laser variables
      DIMENSION alpha1(kk),xaser1(kk),xl1(kk),xl3(kk)
      COMMON /comlas/ alpha1,elas1,xamda1,xaser1,xecri1,plas1,rabs1,
     &                rocri1,xl1,xl3,nabs1
c end of laser variables
c c2.5 thermonuclear reactions
      DIMENSION f1d(kk),f1h(kk),f1he3(kk),f1he4(kk),f1neu(kk),f1ntrl(kk)
     &          ,f1t(kk),f1x(kk),f1x1(kk),f1x2(kk),f1x3(kk),f1x4(kk),
     &          f3d(kk),f3h(kk),f3he3(kk),f3he4(kk),f3neu(kk),f3ntrl(kk)
     &          ,f3t(kk),f3x(kk),f3x1(kk),f3x2(kk),f3x3(kk),f3x4(kk),
     &          r1dd(kk),r1dhe3(kk),r1dt(kk),r3dd(kk),r3dhe3(kk),
     &          r3dt(kk),ye1(kk),ye3(kk),yi1(kk),yi3(kk)
      COMMON /comfus/ deuter,f1d,f1h,f1he3,f1he4,f1neu,f1ntrl,f1t,f1x,
     &                f1x1,f1x2,f1x3,f1x4,f3d,f3h,f3he3,f3he4,f3neu,
     &                f3ntrl,f3t,f3x,f3x1,f3x2,f3x3,f3x4,heliu3,heliu4,
     &                hydrog,xetral,xtrlms,pneut1,pneut3,r1dd,r1dhe3,
     &                r1dt,r3dd,r3dhe3,r3dt,rneut1,rneut3,tinucl,totneu,
     &                tritiu,xmass,xmass1,xmass2,xmass3,xmass4,xtra,
     &                xtra1,xtra2,xtra3,xtra4,xz,xz1,xz2,xz3,xz4,ye1,
     &                ye3,yi1,yi3,yield
c end of thermonuclear reactions
c c2.6 - energies
      COMMON /comen / eerror,en,einput,eloss,eefuse,eifuse,eindt1,
     &                eindt3,eneutr,pv,usqm,xionen,en0,eions,ecbbs,
     &                ecbfs,edcvs,erbfs,eatis,ezsts,erbbs,erfbs,edfbs
c c2.6 - end energies
c c2.7 physics control
      LOGICAL mlbrms,mlecon,mlfuse,mlicon,mlx,nlabs,nlbrms,nlburn,
     &        nlcri1,nldepo,nlecon,nlfuse,nlicon,nlmove,nlpfe,nlpfi,nlx,
     &        nlhf,nltnl
      COMMON /comcon/ mstep,ncase,ngeom,mlbrms,mlecon,mlfuse,mlicon,mlx,
     &                nlabs,nlbrms,nlburn,nlcri1,nldepo,nlecon,nlfuse,
     &                nlicon,nlmove,nlpfe,nlpfi,nlx,nlhf,nltnl,pmin,
     &                rhomin,temin,timin,umin
c end of physics control
c c3.1 numerical control parameters
      LOGICAL break,nlgoon,nlite
      COMMON /comnc / nceldt,ncondt,nit,nitmax,break,nlgoon,nlite,ak0,
     &                ak1,ak2,ak3,ak4,ak5,ak6,ak7,ak8,ak9,bneum,deltat,
     &                dt2,dt3,dt4,dtemax,dtimax,dtmax,dtmin,dumax,rdt2,
     &                rdt3,rdt4
c end of numerical control parameters
c c3.2 mesh and numerical methods
      DIMENSION ae(kk),ai(kk),be(kk),bi(kk),ce(kk),ci(kk),de(kk),di(kk),
     &          e(kk),f(kk),ge(kk),gi(kk),q2(kk),q4(kk),teite(kk),
     &          tiite(kk),uite(kk)
      COMMON /comnum/ mesh,nj,njm1,njp1,nl,nlm1,nlp1,ae,ai,be,bi,ce,ci,
     &                de,di,e,f,ge,gi,q2,q4,teite,tiite,uite
c end of mesh and numerical methods
c c4.1 - administrative variables
      LOGICAL nldump,nlemp
      DIMENSION piq(302)
      COMMON /comadm/ maxdim,maxrun,ndump,nrep,nldump,nlemp,piq,tstop
c c4.1 - end administrative variables
c c?.?? ??????????
      DIMENSION tfmult(kk)
      COMMON /comtfc/ kstart,tfmult
c end of ??????????
      COMMON /dpi   / ecbb(kk),ecbf(kk),eraddr(kk),erbf(kk),eati(kk),
     &                dcvz(kk),efbx2(kk),efrx2(kk),erbb(kk),erfb(kk),
     &                acfc(10,kk),acfr(10,kk),bcfc(10,kk),bcfr(10,kk),
     &                edfb(kk),ter(kk),dqdtim(kk),nmaxr
c instantaneous rates
      COMMON/wrats/weefus,weifus,wrad,wionz,wcbbz,wcbfz,wrbfz,
     & watiz,wrfbz,wdfbz,wrbbz
      DATA iclass,isub/2,22/
      IF(nlomt2(isub))RETURN
      zk=1.38E-23/(gammae-1.)
c 1. energy contents at level 3
      igeom=ngeom-1
      zrg1m1=1.0/(gammai-1.0)
      zke=0.0
      zthe=0.0
      zionen=0.0
      DO 100 l=1,nl
c ad1
       IF(nlpfe.or.piq(59).le.0.5)THEN
        zthe=zthe+zrg1m1*pe3(l)*v3(l)*dm(l)
       ELSE
c eq.79
        qqtfz=effz(l)
        qqtft=te3(l)
        qqtfv=xmieff(l)*pmass/rho3(l)
        qqtfc=tfmult(l)
        zthe=zthe+dm(l)*v3(l)*tfu(qqtft,qqtfv,qqtfz,qqtfc)
       ENDIF
       zthe=zthe+zrg1m1*pi3(l)*v3(l)*dm(l)
 100  CONTINUE
      DO 200 l=2,nl
c eq.80
       j=l
c the mass is an average of two cell masses
       zdm=0.5*(dm(l-1)+dm(l))
       zu2sq=u2(j)*u2(j)
       zu4sq=u4(j)*u4(j)
       zke=zke+0.25*zdm*(zu2sq+zu4sq)
c the internal energy of ionisation  fidash is in joule / kg, tf case
 200  CONTINUE
c the internal energy of ionisation tf case
c already included in ddroe
c     do 120 l = 1, nl
c        zionen = zionen + fidash(l) * dm(l)
c 120 continue
c 2. find energy at the boundary
c boundary velocity defined here
      uedge=0.5*(u4(njp1)+u2(njp1))
c the inner cell has no kinetic energy, only thermal
c the outer cell has kinetic energy but carries only a :half: mass
      zke=zke+0.25*dm(nl)*0.5*(u2(njp1)*u2(njp1)+u4(njp1)*u4(njp1))
      pv=zthe
      usqm=zke
      xionen=zionen
      en=pv+usqm+xionen
c 3. energy input from level 1 to 3
      IF(ncase.eq.1)THEN
      ELSEIF(ncase.eq.2)THEN
c 3.2 temperatures determined
c in this case heat enters as two fluxes, arising from altering ti
c and te at the boundaries
       z1dr=1.0/((r3(njp1)-r3(nj))*2.0)
c the total heat fluxes are determined by the wall temperature
c gradient and the averaged heat conductivity
       zfi=0.0
       zfe=0.0
       zfi=r3(njp1)**igeom*(ti3(nlp1)-ti3(nl))*z1dr*xappai(njp1)
       zfe=r3(njp1)**igeom*(te3(nlp1)-te3(nl))*z1dr*xappae(njp1)
c eq.82
       eindt3=zfi+zfe
       GOTO 500
      ELSEIF(ncase.eq.3)THEN
       GOTO 400
      ELSEIF(ncase.eq.4)THEN
c 3.4 velocity determined
c the pressure increase is proportional to the acceleration
       zdudt3=(u4(njp1)-u2(njp1))/dt3
       zq3=0.5*(q4(nl)+q2(nl))
       zdp=-0.5*dm(nl)*zdudt3/(r3(njp1)**igeom)
c eq.84
       p3(nlp1)=p3(nl)+zq3+zdp
c the energy input can be worked out as above
       GOTO 400
      ELSE
       CALL gotoer(' MEDUSA: energy: goto error, ncase is',ncase)
      ENDIF
c 3.1 thermally insulated wall
      zex=0.0
c in this case has received an amount xl3(l) per unit mass
      DO 300 l=1,nl
c eq.81
       zex=zex+xl3(l)*dm(l)
 300  CONTINUE
      eindt3=zex
      GOTO 500
c 3.3 pressure determined
c the energy input is the work done by the applied pressure
c work out zpedge at level 2
 400  IF(mstep.eq.0)zpedge=p3(nlp1)
      zpedge=(p3(nlp1)+zpedge)*0.5
c eq.83
      eindt3=-zpedge*u2(njp1)*0.5*(r3(njp1)**igeom+r1(njp1)**igeom)
      zpedge=p3(nlp1)
c 4. the total energy input and gain at level 3
c eq.85
 500  einput=einput+dt2*0.5*(eindt1+eindt3)
c the total energy loss through radiation and the total energy
c deposits from thermonuclear reactions
      zrad=0.0
      zeefus=0.0
      zeifus=0.0
      zionz=0.0
      zcbbz=0.0
      zcbfz=0.0
      zrrdz=0.
      zrbfz=0.0
      zatiz=0.
      edcvs=0.0
      zrbbz=0.0
      zrfbz=0.0
      zdfbz=0.0
      DO 600 l=1,nl
c eq.86
       zeefus=zeefus+(ye3(l)+ye1(l))*dm(l)
       zeifus=zeifus+(yi3(l)+yi1(l))*dm(l)
c eq.87
       zrad=zrad+(brems3(l)+brems1(l))*dm(l)
c ad1
       IF(abs(piq(58)-2.0).lt.underf)THEN
c internal energy of ionisation fe case
c
        zionz=zionz+2.*(efbx2(l))*dm(l)
        zcbbz=zcbbz+2.*(ecbb(l))*dm(l)
        zcbfz=zcbfz+2.*(ecbf(l))*dm(l)
        zrbfz=zrbfz+2.*(erbf(l))*dm(l)
        zatiz=zatiz+2.*(eati(l))*dm(l)
        zrfbz=zrfbz+2.*(erfb(l))*dm(l)
        zdfbz=zdfbz+2.*(edfb(l))*dm(l)
        zrbbz=zrbbz+2.*(erbb(l))*dm(l)
       ENDIF
 600  CONTINUE
c instantaneous rates
        weefus=zeefus/2.
        weifus=zeifus/2.
        wrad=zrad/2.
        wionz=zionz/2.
        wcbbz=zcbbz/2.
        wcbfz=zcbfz/2.
        wrbfz=zrbfz/2.
        watiz=zatiz/2.
        wrfbz=zrfbz/2.
        wdfbz=zdfbz/2.
        wrbbz=zrbbz/2.
c
      eefuse=eefuse+0.5*zeefus*dt2
      eifuse=eifuse+0.5*zeifus*dt2
c the energy carried away by the neutrons
      eneutr=eneutr+dt2*0.5*(pneut1+pneut3)
c the total number of neutrons released
      totneu=totneu+dt2*0.5*(rneut1+rneut3)
      eloss=eloss+0.5*zrad*dt2
      elas1=elas1+plas1*dt2
c the error in the energy calculation (eq.88)
      eerror=(en-eloss)-einput-eefuse-eifuse
c ad1 ionisation energy includes tdcv  as calculation of te3 used zstar
c     from previous timestep
      IF(abs(piq(58)-2.0).lt.underf)THEN
       eions=eions+0.5*dt2*zionz
       ecbbs=ecbbs+0.5*dt2*zcbbz
       ecbfs=ecbfs+0.5*dt2*zcbfz
       erbfs=erbfs+0.5*dt2*zrbfz
       eatis=eatis+0.5*dt2*zatiz
       erbbs=erbbs+0.5*dt2*zrbbz
       erfbs=erfbs+0.5*dt2*zrfbz
       edfbs=edfbs+0.5*dt2*zdfbz
c      write(6,*)eions,.5*dt2*zionz,ecbbs,0.5*dt2*zcbbz,
c    & ecbfs,.5*dt2*zcbfz
      ENDIF
c take account of ionisation energy in eerror (not including this timest
      IF(abs(piq(58)-2.0).lt.underf.and.(nlpfe.or.piq(59).lt.0.5))
     &   eerror=eerror-eions-erbfs-eatis
      zratio=abs(eerror)/en
      IF(zratio.lt.1.0E-7)eerror=0.0
c thermonuclear yield includes neutron energy
      IF(abs(elas1).gt.underf)yield=(eefuse+eifuse+eneutr)/elas1
c 5. the lawson criterion
      zrhor=0.0
      zdm=0.0
      DO 700 l=1,nl
       j=l
       zdm=zdm+dm(l)
       zrhor=zrhor+rho3(l)*(r3(j)+r3(j+1))*dm(l)
 700  CONTINUE
c the factor 0.5 in rhor arises from averaging r
      rhor=0.5*zrhor/zdm
      CALL expert(iclass,isub,1)
      RETURN
      END


      SUBROUTINE exam
      PARAMETER(kk=1001)
c 2.26 examines present state of system
c c1.1 basic system parameters
      LOGICAL nlend,nlres
      COMMON /combas/ altime,cptime,stime,ndiary,nin,nledge,nonlin,nout,
     &                nprint,npunch,nread,nrec,nresm,nrun,nstep,nlend,
     &                nlres
c end of basic system parameters
c c1.9 development and diagnostic parameters
      DIMENSION nadump(20),npdump(20),nvdump(20)
      LOGICAL nlched,nlhead(9),nlomt1(50),nlomt2(50),nlomt3(50),nlrept
      COMMON /comddp/ maxdum,mxdump,nadump,nclass,npdump,npoint,nsub,
     &                nvdump,nlched,nlhead,nlomt1,nlomt2,nlomt3,nlrept
c end of development and diagnostic parameters
c c2.1 hydrodynamics
      DIMENSION dm(kk),p3(kk),r1(kk),r3(kk),r5(kk),rho1(kk),rho3(kk),
     &          rho5(kk),u2(kk),u4(kk),v1(kk),v3(kk),v5(kk),lstat(kk)
      COMMON /comhyd/ dm,p3,pini,r1,r3,r5,rho1,rho3,rho5,rhoini,rhor,
     &                rini,time,u2,u4,uedge,v1,v3,v5
c end of hydrodynamics
c c2.2 thermodynamics
      DIMENSION ddroe1(kk),ddroe3(kk),ddroi1(kk),ddroi3(kk),ddte1(kk),
     &          ddte3(kk),ddti1(kk),ddti3(kk),xappae(kk),xappai(kk),
     &          pe1(kk),pe3(kk),pi1(kk),pi3(kk),te1(kk),te3(kk),ti1(kk),
     &          ti3(kk)
      COMMON /comth / ddroe1,ddroe3,ddroi1,ddroi3,ddte1,ddte3,ddti1,
     &                ddti3,gammae,gammai,xappae,xappai,pe1,pe3,pi1,pi3,
     &                te1,te3,teini,ti1,ti3,tiini
c end of thermodynamics
c c2.3 ions and electrons
c end of ions and electrons
c c2.5 thermonuclear reactions
      DIMENSION f1d(kk),f1h(kk),f1he3(kk),f1he4(kk),f1neu(kk),f1ntrl(kk)
     &          ,f1t(kk),f1x(kk),f1x1(kk),f1x2(kk),f1x3(kk),f1x4(kk),
     &          f3d(kk),f3h(kk),f3he3(kk),f3he4(kk),f3neu(kk),f3ntrl(kk)
     &          ,f3t(kk),f3x(kk),f3x1(kk),f3x2(kk),f3x3(kk),f3x4(kk),
     &          r1dd(kk),r1dhe3(kk),r1dt(kk),r3dd(kk),r3dhe3(kk),
     &          r3dt(kk),ye1(kk),ye3(kk),yi1(kk),yi3(kk)
      COMMON /comfus/ deuter,f1d,f1h,f1he3,f1he4,f1neu,f1ntrl,f1t,f1x,
     &                f1x1,f1x2,f1x3,f1x4,f3d,f3h,f3he3,f3he4,f3neu,
     &                f3ntrl,f3t,f3x,f3x1,f3x2,f3x3,f3x4,heliu3,heliu4,
     &                hydrog,xetral,xtrlms,pneut1,pneut3,r1dd,r1dhe3,
     &                r1dt,r3dd,r3dhe3,r3dt,rneut1,rneut3,tinucl,totneu,
     &                tritiu,xmass,xmass1,xmass2,xmass3,xmass4,xtra,
     &                xtra1,xtra2,xtra3,xtra4,xz,xz1,xz2,xz3,xz4,ye1,
     &                ye3,yi1,yi3,yield
c end of thermonuclear reactions
c c2.7 physics control
      LOGICAL mlbrms,mlecon,mlfuse,mlicon,mlx,nlabs,nlbrms,nlburn,
     &        nlcri1,nldepo,nlecon,nlfuse,nlicon,nlmove,nlpfe,nlpfi,nlx,
     &        nlhf,nltnl
      COMMON /comcon/ mstep,ncase,ngeom,mlbrms,mlecon,mlfuse,mlicon,mlx,
     &                nlabs,nlbrms,nlburn,nlcri1,nldepo,nlecon,nlfuse,
     &                nlicon,nlmove,nlpfe,nlpfi,nlx,nlhf,nltnl,pmin,
     &                rhomin,temin,timin,umin
c end of physics control
c c3.1 numerical control parameters
      LOGICAL break,nlgoon,nlite
      COMMON /comnc / nceldt,ncondt,nit,nitmax,break,nlgoon,nlite,ak0,
     &                ak1,ak2,ak3,ak4,ak5,ak6,ak7,ak8,ak9,bneum,deltat,
     &                dt2,dt3,dt4,dtemax,dtimax,dtmax,dtmin,dumax,rdt2,
     &                rdt3,rdt4
c end of numerical control parameters
c c3.2 mesh and numerical methods
      DIMENSION ae(kk),ai(kk),be(kk),bi(kk),ce(kk),ci(kk),de(kk),di(kk),
     &          e(kk),f(kk),ge(kk),gi(kk),q2(kk),q4(kk),teite(kk),
     &          tiite(kk),uite(kk)
      COMMON /comnum/ mesh,nj,njm1,njp1,nl,nlm1,nlp1,ae,ai,be,bi,ce,ci,
     &                de,di,e,f,ge,gi,q2,q4,teite,tiite,uite
c end of mesh and numerical methods
c c4.1 - administrative variables
      LOGICAL nldump,nlemp
      DIMENSION piq(302)
      COMMON /comadm/ maxdim,maxrun,ndump,nrep,nldump,nlemp,piq,tstop
c c4.1 - end administrative variables
      LOGICAL ls1,ls2,ls3,ls4,ls5
      DATA iclass,isub/2,26/
      break=.false.
      IF(nlomt2(isub))RETURN
c 1. are volumes or temperatures negative ?
      j=0
 100  j=j+1
      l=j
      IF(r5(j+1).le.r5(j))THEN
       icell=j
       i=-2
       GOTO 400
      ELSEIF(ti3(l).gt.-timin)THEN
       IF(ti3(l).lt.timin)ti3(l)=timin
       IF(te3(l).gt.-temin)THEN
        IF(te3(l).lt.temin)te3(l)=temin
        IF(j.lt.nj)GOTO 100
c 2. examine temperatures
        zg=gammae-1.0
        lstop=nlm1
        IF(ncase.eq.2)lstop=nl
c set switches initially
        ls1=.false.
        ls2=.false.
        ls3=.false.
        ls4=.false.
        ls5=.false.
        l=0
       ELSE
        icell=l
        i=-4
        GOTO 400
       ENDIF
      ELSE
       icell=l
       i=-3
       GOTO 400
      ENDIF
 200  l=l+1
c 2.1 are all heat fluxes zero?
      zdti=abs(ti3(l+1)-ti3(l))
      zdte=abs(te3(l+1)-te3(l))
      IF(zdti.ge.timin)ls1=.true.
      IF(zdte.ge.temin)ls2=.true.
c 2.2 is bremsstrahlung negligible?
      ztem=te3(l)
      IF(nlpfe)THEN
       rho0=rhoini
       tein1=teini
       IF(l.gt.mesh-piq(63).and.piq(63).gt.0) then
            rho0=piq(62)
            tein1=piq(89)
       elseIF(l.gt.mesh-piq(13)-piq(63).and.piq(13).gt.0) then
           rho0=piq(12)
            tein1=piq(86)
       endif
       ztadia=tein1*(rho3(l)/rho0)**zg
       zdiff=te3(l)-ztadia
       ztem=0.0
       IF(zdiff.gt.0.0)ztem=zdiff
      ENDIF
      IF(ztem.gt.temin)ls3=.true.
c 2.3 is the exchange rate important?
c ad1
      IF(abs(ti3(l)-te3(l)).gt.timin)ls4=.true.
c 2.4 thermonuclear reactions?
      IF(ti3(l).ge.tinucl)ls5=.true.
      IF(l.lt.lstop)GOTO 200
c 2.5 switch particle processes on as appropriate
      mlicon=ls1
      mlecon=ls2
      mlbrms=ls3
      mlx=ls4
      mlfuse=ls5
c 3. are basic variables below permitted minima?
      DO 300 l=1,nl
       icell=l
       i=2
       IF(rho3(l).lt.rhomin)GOTO 400
       i=3
       IF(ti3(l).lt.timin)GOTO 400
       i=4
       IF(te3(l).lt.temin)GOTO 400
 300  CONTINUE
      RETURN
c 4. terminate the calculation
 400  CALL page
      CALL blines(2)
      CALL mesage('                          calculation terminated')
      CALL mesage('                          ----------------------')
      CALL blines(2)
      IF(i.eq.-2)CALL mesage('      volumes negative')
      IF(i.eq.-3)CALL mesage('      ion temperature negative')
      IF(i.eq.-4)CALL mesage('      electron temperature negative')
      IF(i.eq.2)CALL mesage('      densities too small')
      IF(i.eq.3)CALL mesage('      ion temperature too small')
      IF(i.eq.4)CALL mesage('      electron temperature too small')
      CALL blines(2)
      CALL ivar('cell no ',icell)
      CALL ivar('nstep   ',nstep)
      CALL rvar('time    ',time)
c terminate calculation
c force an abnormal exit.
c     if ( i.gt.0 ) nlend = .true.
      break=.true.
      RETURN
      END


      SUBROUTINE expert(kclass,ksub,kpoint)
c 0.4 modify standard operation of program
c c1.9 development and diagnostic parameters
      DIMENSION nadump(20),npdump(20),nvdump(20)
      LOGICAL nlched,nlhead(9),nlomt1(50),nlomt2(50),nlomt3(50),nlrept
      COMMON /comddp/ maxdum,mxdump,nadump,nclass,npdump,npoint,nsub,
     &                nvdump,nlched,nlhead,nlomt1,nlomt2,nlomt3,nlrept
c end of development and diagnostic parameters
c set tracer variables
      nclass=kclass
      nsub=ksub
      npoint=kpoint
c are diagnostics required?
      IF(nlrept)CALL report(kclass,ksub,kpoint)
      RETURN
      END


      SUBROUTINE findt
      PARAMETER(kk=1001)
c 2.15 the temperature equations (gauss-elimination)
c c1.9 development and diagnostic parameters
      DIMENSION nadump(20),npdump(20),nvdump(20)
      LOGICAL nlched,nlhead(9),nlomt1(50),nlomt2(50),nlomt3(50),nlrept
      COMMON /comddp/ maxdum,mxdump,nadump,nclass,npdump,npoint,nsub,
     &                nvdump,nlched,nlhead,nlomt1,nlomt2,nlomt3,nlrept
c end of development and diagnostic parameters
c c2.2 thermodynamics
      DIMENSION ddroe1(kk),ddroe3(kk),ddroi1(kk),ddroi3(kk),ddte1(kk),
     &          ddte3(kk),ddti1(kk),ddti3(kk),xappae(kk),xappai(kk),
     &          pe1(kk),pe3(kk),pi1(kk),pi3(kk),te1(kk),te3(kk),ti1(kk),
     &          ti3(kk)
      COMMON /comth / ddroe1,ddroe3,ddroi1,ddroi3,ddte1,ddte3,ddti1,
     &                ddti3,gammae,gammai,xappae,xappai,pe1,pe3,pi1,pi3,
     &                te1,te3,teini,ti1,ti3,tiini
c end of thermodynamics
c c3.2 mesh and numerical methods
      DIMENSION ae(kk),ai(kk),be(kk),bi(kk),ce(kk),ci(kk),de(kk),di(kk),
     &          e(kk),f(kk),ge(kk),gi(kk),q2(kk),q4(kk),teite(kk),
     &          tiite(kk),uite(kk)
      COMMON /comnum/ mesh,nj,njm1,njp1,nl,nlm1,nlp1,ae,ai,be,bi,ce,ci,
     &                de,di,e,f,ge,gi,q2,q4,teite,tiite,uite
c end of mesh and numerical methods
      LOGICAL mlbrms,mlecon,mlfuse,mlicon,mlx,nlabs,nlbrms,nlburn,
     &        nlcri1,nldepo,nlecon,nlfuse,nlicon,nlmove,nlpfe,nlpfi,nlx,
     &        nlhf,nltnl
      COMMON /comcon/ mstep,ncase,ngeom,mlbrms,mlecon,mlfuse,mlicon,mlx,
     &                nlabs,nlbrms,nlburn,nlcri1,nldepo,nlecon,nlfuse,
     &                nlicon,nlmove,nlpfe,nlpfi,nlx,nlhf,nltnl,pmin,
     &                rhomin,temin,timin,umin
      LOGICAL nldump,nlemp
      DIMENSION piq(302)
      COMMON /comadm/ maxdim,maxrun,ndump,nrep,nldump,nlemp,piq,tstop
      DATA iclass,isub/2,15/
      IF(nlomt2(isub))RETURN
c 1. solve the ion-temperature equation (62)
c 1.1 loop forwards to find the e and f coeffs
c eq.68
      e(1)=ci(1)/bi(1)
      f(1)=(di(1)+gi(1))/bi(1)
      DO 100 l=2,nl
       znumi=1.0/(bi(l)-ai(l)*e(l-1))
c eq.69
       e(l)=ci(l)*znumi
       f(l)=(di(l)+gi(l)-ai(l)*f(l-1))*znumi
 100  CONTINUE
c 1.2      find ti starting from the boundary
c eq.70
      ti3(nl)=f(nl)
      DO 200 ldash=1,nlm1
       l=nl-ldash
c eq.71
       ti3(l)=f(l)-e(l)*ti3(l+1)
 200  CONTINUE
c 2. solve the electron-temperature equation (62)
c 2.1 loop forwards to find the e and f coeffs
c eq.68
      e(1)=ce(1)/be(1)
      f(1)=(de(1)+ge(1))/be(1)
      DO 300 l=2,nl
       znume=1.0/(be(l)-ae(l)*e(l-1))
c eq.69
       e(l)=ce(l)*znume
       f(l)=(de(l)+ge(l)-ae(l)*f(l-1))*znume
 300  CONTINUE
c 2.2 find te starting from the boundary
c eq.70
      te3(nl)=f(nl)
      DO 400 ldash=1,nlm1
       l=nl-ldash
c eq.71
       te3(l)=f(l)-e(l)*te3(l+1)
 400  CONTINUE
c in the tf case
c acoustic
c     IF(mstep.eq.1)THEN
c     do 26 l=1,10
c      te3(l)=700.
c      ti3(l)=700.
c26   continue
c     ENDIF


      RETURN
      END


      SUBROUTINE formp
      PARAMETER(kk=1001,underf=1.E-30)
c 2.24 forms the hydrodynamic pressure
c c1.9 development and diagnostic parameters
      DIMENSION nadump(20),npdump(20),nvdump(20)
      LOGICAL nlched,nlhead(9),nlomt1(50),nlomt2(50),nlomt3(50),nlrept
      COMMON /comddp/ maxdum,mxdump,nadump,nclass,npdump,npoint,nsub,
     &                nvdump,nlched,nlhead,nlomt1,nlomt2,nlomt3,nlrept
c end of development and diagnostic parameters
      LOGICAL lphi(kk)
      DIMENSION deltae(50),ebdy(51),emid(50),ha(kk),hb(kk),hc(kk),hd(kk)
     &          ,hdiff(kk),he(kk),heta(kk),hf(kk),hg(kk),hgrp(kk),hj(kk)
     &          ,hlas(kk),hloss(50),hncold(kk),hp(kk),hphi(kk),htaug(kk)
     &          ,htherm(kk),htr(kk),xhot1(kk,50),xhot3(kk,50)
      COMMON /comhot/ igrp,ngroup,lphi,deltae,ebdy,ehot,emid,ha,hb,hc,
     &                hd,hdiff,he,heta,hf,hg,hgrp,hj,hlas,hloss,hncold,
     &                hp,hphi,htaug,htherm,htr,xhot1,xhot3,thot
c c4.1 - administrative variables
c c4.1 - end administrative variables
c c2.3 ions and electrons
c end of ions and electrons
c c2.4 laser variables
c end of laser variables
c c2.1 hydrodynamics
      DIMENSION dm(kk),p3(kk),r1(kk),r3(kk),r5(kk),rho1(kk),rho3(kk),
     &          rho5(kk),u2(kk),u4(kk),v1(kk),v3(kk),v5(kk),lstat(kk)
      COMMON /comhyd/ dm,p3,pini,r1,r3,r5,rho1,rho3,rho5,rhoini,rhor,
     &                rini,time,u2,u4,uedge,v1,v3,v5
c end of hydrodynamics
c c2.2 thermodynamics
      DIMENSION ddroe1(kk),ddroe3(kk),ddroi1(kk),ddroi3(kk),ddte1(kk),
     &          ddte3(kk),ddti1(kk),ddti3(kk),xappae(kk),xappai(kk),
     &          pe1(kk),pe3(kk),pi1(kk),pi3(kk),te1(kk),te3(kk),ti1(kk),
     &          ti3(kk)
      COMMON /comth / ddroe1,ddroe3,ddroi1,ddroi3,ddte1,ddte3,ddti1,
     &                ddti3,gammae,gammai,xappae,xappai,pe1,pe3,pi1,pi3,
     &                te1,te3,teini,ti1,ti3,tiini
c end of thermodynamics
c c2.7 physics control
      LOGICAL mlbrms,mlecon,mlfuse,mlicon,mlx,nlabs,nlbrms,nlburn,
     &        nlcri1,nldepo,nlecon,nlfuse,nlicon,nlmove,nlpfe,nlpfi,nlx,
     &        nlhf,nltnl
      COMMON /comcon/ mstep,ncase,ngeom,mlbrms,mlecon,mlfuse,mlicon,mlx,
     &                nlabs,nlbrms,nlburn,nlcri1,nldepo,nlecon,nlfuse,
     &                nlicon,nlmove,nlpfe,nlpfi,nlx,nlhf,nltnl,pmin,
     &                rhomin,temin,timin,umin
c end of physics control
c c3.2 mesh and numerical methods
      DIMENSION ae(kk),ai(kk),be(kk),bi(kk),ce(kk),ci(kk),de(kk),di(kk),
     &          e(kk),f(kk),ge(kk),gi(kk),q2(kk),q4(kk),teite(kk),
     &          tiite(kk),uite(kk)
      COMMON /comnum/ mesh,nj,njm1,njp1,nl,nlm1,nlp1,ae,ai,be,bi,ce,ci,
     &                de,di,e,f,ge,gi,q2,q4,teite,tiite,uite
c end of mesh and numerical methods
c c2.3 ions and electrons
      DIMENSION brems1(kk),brems3(kk),ddz1(kk),ddz3(kk),degen(kk),
     &          effz(kk),eixch2(kk),fidash(kk),fz1(kk),fz3(kk),fzsq1(kk)
     &          ,fzsq3(kk),xlc(kk),xmieff(kk),xne(kk),xni(kk),omega1(kk)
      COMMON /comie / brems1,brems3,ddz1,ddz3,degen,degmax,degmin,effz,
     &                eixch2,fidash,fz1,fz3,fzsq1,fzsq3,xlc,xmieff,xne,
     &                xni,omega1,pmass
c end of ions and electrons
      COMMON /radflg/ yamda,xhnu,xihnu,dihnu,xplas1(kk),xrbf(kk),
     &                atrbf(kk),corhf1(kk),corhf(kk),iradf1,iradf2
c c4.1 - administrative variables
      LOGICAL nldump,nlemp
      DIMENSION piq(302)
      COMMON /comadm/ maxdim,maxrun,ndump,nrep,nldump,nlemp,piq,tstop
c c4.1 - end administrative variables
c c2.4 laser variables
      DIMENSION alpha1(kk),xaser1(kk),xl1(kk),xl3(kk)
      COMMON /comlas/ alpha1,elas1,xamda1,xaser1,xecri1,plas1,rabs1,
     &                rocri1,xl1,xl3,nabs1
c end of laser variables
      DATA iclass,isub/2,24/
      IF(nlomt2(isub))RETURN
c 1. the hydrodynamic pressure (eq.36)
      DO 100 l=1,nl
       p3(l)=pi3(l)+pe3(l)
 100  CONTINUE
c include the momentum of the radiation field for all cells
c farther out than the critical density
c the ponderomotive force requires that we calculate the
c local value of <e*e>   we do this in the w.k.b. approximation
c which is invalid in the case of strong profile modification
c ad1  see speed for alternative
      IF(nabs1.gt.1)THEN
       zprad=plas1/(3.0E+08*r3(nabs1)**(ngeom-1))
       labs1=max(nabs1,1)
       DO 150 l=labs1,nl
        zprad=xplas1(l+1)/(3.0E+08*r3(l)**(ngeom-1))
        zepsln=1.0-xne(l)/xecri1
        IF(zepsln.lt.0.001)zepsln=0.001
        p3(l)=p3(l)+piq(56)*zprad/sqrt(zepsln)
 150   CONTINUE
      ENDIF
c     include the pressure of the thermal radiation
      IF(ngroup.ne.0)THEN
       DO 200 ig=1,ngroup
        DO 160 l=1,nl
         p3(l)=p3(l)+0.333333*xhot1(l,ig)
 160    CONTINUE
 200   CONTINUE
      ENDIF
c     the pressure at the boundary is set here
      IF(ncase.eq.1)THEN
      ELSEIF(ncase.eq.2)THEN
c        case 2 : immoveable wall
       p3(nlp1)=p3(nl)+0.5*(q4(nl)+q2(nl))
       RETURN
      ELSEIF(ncase.eq.3)THEN
c        case 3 : applied pressure
       CALL boundy(time)
       RETURN
      ELSEIF(ncase.eq.4)THEN
       GOTO 300
      ELSE
       CALL gotoer(' MEDUSA: formp: goto error, ncase is',ncase)
      ENDIF
c     case 1 : thermally insulated wall (vacuum)
      CALL boundy(time)
      RETURN
c     case 4 : applied velocity
c     the pressure associated with prescribed boundary velocity at
c     level 3 is proportional to the acceleration from u2 to u4. this
c     pressure is worked out later, as it is not required yet.
 300  RETURN
      END


      SUBROUTINE fusion
      PARAMETER(kk=1001)
c 2.12 calculate the energy released by thermonuclear reactions
c c1.9 development and diagnostic parameters
      DIMENSION nadump(20),npdump(20),nvdump(20)
      LOGICAL nlched,nlhead(9),nlomt1(50),nlomt2(50),nlomt3(50),nlrept
      COMMON /comddp/ maxdum,mxdump,nadump,nclass,npdump,npoint,nsub,
     &                nvdump,nlched,nlhead,nlomt1,nlomt2,nlomt3,nlrept
c end of development and diagnostic parameters
c c2.1 hydrodynamics
      DIMENSION dm(kk),p3(kk),r1(kk),r3(kk),r5(kk),rho1(kk),rho3(kk),
     &          rho5(kk),u2(kk),u4(kk),v1(kk),v3(kk),v5(kk),lstat(kk)
      COMMON /comhyd/ dm,p3,pini,r1,r3,r5,rho1,rho3,rho5,rhoini,rhor,
     &                rini,time,u2,u4,uedge,v1,v3,v5
c end of hydrodynamics
c c2.2 thermodynamics
      DIMENSION ddroe1(kk),ddroe3(kk),ddroi1(kk),ddroi3(kk),ddte1(kk),
     &          ddte3(kk),ddti1(kk),ddti3(kk),xappae(kk),xappai(kk),
     &          pe1(kk),pe3(kk),pi1(kk),pi3(kk),te1(kk),te3(kk),ti1(kk),
     &          ti3(kk)
      COMMON /comth / ddroe1,ddroe3,ddroi1,ddroi3,ddte1,ddte3,ddti1,
     &                ddti3,gammae,gammai,xappae,xappai,pe1,pe3,pi1,pi3,
     &                te1,te3,teini,ti1,ti3,tiini
c end of thermodynamics
c c2.3 ions and electrons
      DIMENSION brems1(kk),brems3(kk),ddz1(kk),ddz3(kk),degen(kk),
     &          effz(kk),eixch2(kk),fidash(kk),fz1(kk),fz3(kk),fzsq1(kk)
     &          ,fzsq3(kk),xlc(kk),xmieff(kk),xne(kk),xni(kk),omega1(kk)
      COMMON /comie / brems1,brems3,ddz1,ddz3,degen,degmax,degmin,effz,
     &                eixch2,fidash,fz1,fz3,fzsq1,fzsq3,xlc,xmieff,xne,
     &                xni,omega1,pmass
c end of ions and electrons
c c2.5 thermonuclear reactions
      DIMENSION f1d(kk),f1h(kk),f1he3(kk),f1he4(kk),f1neu(kk),f1ntrl(kk)
     &          ,f1t(kk),f1x(kk),f1x1(kk),f1x2(kk),f1x3(kk),f1x4(kk),
     &          f3d(kk),f3h(kk),f3he3(kk),f3he4(kk),f3neu(kk),f3ntrl(kk)
     &          ,f3t(kk),f3x(kk),f3x1(kk),f3x2(kk),f3x3(kk),f3x4(kk),
     &          r1dd(kk),r1dhe3(kk),r1dt(kk),r3dd(kk),r3dhe3(kk),
     &          r3dt(kk),ye1(kk),ye3(kk),yi1(kk),yi3(kk)
      COMMON /comfus/ deuter,f1d,f1h,f1he3,f1he4,f1neu,f1ntrl,f1t,f1x,
     &                f1x1,f1x2,f1x3,f1x4,f3d,f3h,f3he3,f3he4,f3neu,
     &                f3ntrl,f3t,f3x,f3x1,f3x2,f3x3,f3x4,heliu3,heliu4,
     &                hydrog,xetral,xtrlms,pneut1,pneut3,r1dd,r1dhe3,
     &                r1dt,r3dd,r3dhe3,r3dt,rneut1,rneut3,tinucl,totneu,
     &                tritiu,xmass,xmass1,xmass2,xmass3,xmass4,xtra,
     &                xtra1,xtra2,xtra3,xtra4,xz,xz1,xz2,xz3,xz4,ye1,
     &                ye3,yi1,yi3,yield
c end of thermonuclear reactions
c c2.7 physics control
      LOGICAL mlbrms,mlecon,mlfuse,mlicon,mlx,nlabs,nlbrms,nlburn,
     &        nlcri1,nldepo,nlecon,nlfuse,nlicon,nlmove,nlpfe,nlpfi,nlx,
     &        nlhf,nltnl
      COMMON /comcon/ mstep,ncase,ngeom,mlbrms,mlecon,mlfuse,mlicon,mlx,
     &                nlabs,nlbrms,nlburn,nlcri1,nldepo,nlecon,nlfuse,
     &                nlicon,nlmove,nlpfe,nlpfi,nlx,nlhf,nltnl,pmin,
     &                rhomin,temin,timin,umin
c end of physics control
c c3.2 mesh and numerical methods
      DIMENSION ae(kk),ai(kk),be(kk),bi(kk),ce(kk),ci(kk),de(kk),di(kk),
     &          e(kk),f(kk),ge(kk),gi(kk),q2(kk),q4(kk),teite(kk),
     &          tiite(kk),uite(kk)
      COMMON /comnum/ mesh,nj,njm1,njp1,nl,nlm1,nlp1,ae,ai,be,bi,ce,ci,
     &                de,di,e,f,ge,gi,q2,q4,teite,tiite,uite
c end of mesh and numerical methods
c c4.1 - administrative variables
      LOGICAL nldump,nlemp
      DIMENSION piq(302)
      COMMON /comadm/ maxdim,maxrun,ndump,nrep,nldump,nlemp,piq,tstop
c c4.1 - end administrative variables
      DATA iclass,isub/2,12/
      IF(nlomt2(isub))RETURN
c 1. constants for energy release & cross sections
c conversion factor for temperatures from degrees kelvin to kev
      zkev=1.0/1.1606E7
      zevjou=1.609*1.0E-19
c energy release from d-d, d-t and d-he3 processes in joule
      IF(nldepo)THEN
       zedd=4.8*1.0E6*zevjou
       zedt=3.6*1.0E6*zevjou
       zedhe3=18.3*1.0E6*zevjou
      ELSE
c set the energies deposited by reaction products to zero
       zedd=0.0
       zedt=0.0
       zedhe3=0.0
      ENDIF
c neutron energies
      zenedd=2.45*1.0E6*zevjou
      zenedt=14.1*1.0E6*zevjou
c apportionment constants
      ztdd=1.2*1.0E9*(1.0+piq(28))
      ztdt=3.71*1.0E8*(1.0+piq(29))
      ztdhe3=2.0*1.0E9*(1.0+piq(30))
c exponentiation factors
      z13=-1.0/3.0
      z23=-2.0/3.0
c cross section factors
      zcdd=2.33*1.0E-20
      zcdt=3.68*1.0E-18
c approximated reaction rates at high temperatures
      zcdt35=7.5*1.0E-22
      zcdhe3=1.0E-23
c set sums to zero
      zneutr=0.0
      zeneut=0.0
c 2. reaction rates
      DO 100 l=1,nl
c if a cell has not reached the ignition temperature then there
c is no burning of material
       IF(ti3(l).ge.tinucl)THEN
c express ti in kev
        ztikev=ti3(l)*zkev
c cross section for the d-d reaction
        zcsdd=zcdd*exp(-19.42*ztikev**z13)*ztikev**z23
c cross section for d-t reaction below 20 kev
        zcsdt=zcdt*exp(-19.94*ztikev**z13)*ztikev**z23
        IF(ztikev.gt.35.0)zcsdt=zcdt35
c we ignore the d-he3 reaction below 50 kev
        zcsdh3=0.0
        IF(ztikev.gt.50.0)zcsdh3=zcdhe3
c 3. number of thermonuclear reactions
c these are per m**3 * sec
        r3dd(l)=zcsdd*f3d(l)*f3d(l)*xni(l)*xni(l)*0.5
        r3dt(l)=zcsdt*f3d(l)*f3t(l)*xni(l)*xni(l)
        r3dhe3(l)=zcsdh3*f3d(l)*f3he3(l)*xni(l)*xni(l)
c the rate of neutron production (total)
        zneutr=zneutr+(0.5*r3dd(l)+r3dt(l))*v3(l)*dm(l)
c the total neutron energy released
        zeneut=zeneut+(0.5*r3dd(l)*zenedd+r3dt(l)*zenedt)*v3(l)*dm(l)
c 4. rate of energy given to ions and electrons
c energy liberated per kg*sec
        zydd=r3dd(l)*zedd*v3(l)
        zydt=r3dt(l)*zedt*v3(l)
        zydhe3=r3dhe3(l)*zedhe3*v3(l)
c apportionment factor for the ions
        zpdd=te3(l)/(te3(l)+ztdd)
        zpdt=te3(l)/(te3(l)+ztdt)
        zpdhe3=te3(l)/(te3(l)+ztdhe3)
c rate of energy given per kg
        yi3(l)=zpdd*zydd+zpdt*zydt+zpdhe3*zydhe3
        ye3(l)=(1.0-zpdd)*zydd+(1.0-zpdt)*zydt+(1.0-zpdhe3)*zydhe3
       ENDIF
 100  CONTINUE
c 5. store neutron data
c store sums of total neutron rates and energies
      rneut3=zneutr
      pneut3=zeneut
      RETURN
      END


      SUBROUTINE gie
      PARAMETER(kk=1001)
c 2.4 the terms of the energy equation at level n-1
c c1.9 development and diagnostic parameters
      DIMENSION nadump(20),npdump(20),nvdump(20)
      LOGICAL nlched,nlhead(9),nlomt1(50),nlomt2(50),nlomt3(50),nlrept
      COMMON /comddp/ maxdum,mxdump,nadump,nclass,npdump,npoint,nsub,
     &                nvdump,nlched,nlhead,nlomt1,nlomt2,nlomt3,nlrept
c end of development and diagnostic parameters
c c2.1 hydrodynamics
      DIMENSION dm(kk),p3(kk),r1(kk),r3(kk),r5(kk),rho1(kk),rho3(kk),
     &          rho5(kk),u2(kk),u4(kk),v1(kk),v3(kk),v5(kk),lstat(kk)
      COMMON /comhyd/ dm,p3,pini,r1,r3,r5,rho1,rho3,rho5,rhoini,rhor,
     &                rini,time,u2,u4,uedge,v1,v3,v5
c end of hydrodynamics
c c2.2 thermodynamics
      DIMENSION ddroe1(kk),ddroe3(kk),ddroi1(kk),ddroi3(kk),ddte1(kk),
     &          ddte3(kk),ddti1(kk),ddti3(kk),xappae(kk),xappai(kk),
     &          pe1(kk),pe3(kk),pi1(kk),pi3(kk),te1(kk),te3(kk),ti1(kk),
     &          ti3(kk)
      COMMON /comth / ddroe1,ddroe3,ddroi1,ddroi3,ddte1,ddte3,ddti1,
     &                ddti3,gammae,gammai,xappae,xappai,pe1,pe3,pi1,pi3,
     &                te1,te3,teini,ti1,ti3,tiini
c end of thermodynamics
c c2.3 ions and electrons
      DIMENSION brems1(kk),brems3(kk),ddz1(kk),ddz3(kk),degen(kk),
     &          effz(kk),eixch2(kk),fidash(kk),fz1(kk),fz3(kk),fzsq1(kk)
     &          ,fzsq3(kk),xlc(kk),xmieff(kk),xne(kk),xni(kk),omega1(kk)
      COMMON /comie / brems1,brems3,ddz1,ddz3,degen,degmax,degmin,effz,
     &                eixch2,fidash,fz1,fz3,fzsq1,fzsq3,xlc,xmieff,xne,
     &                xni,omega1,pmass
c end of ions and electrons
c c2.4 laser variables
      DIMENSION alpha1(kk),xaser1(kk),xl1(kk),xl3(kk)
      COMMON /comlas/ alpha1,elas1,xamda1,xaser1,xecri1,plas1,rabs1,
     &                rocri1,xl1,xl3,nabs1
c end of laser variables
c c2.5 thermonuclear reactions
      DIMENSION f1d(kk),f1h(kk),f1he3(kk),f1he4(kk),f1neu(kk),f1ntrl(kk)
     &          ,f1t(kk),f1x(kk),f1x1(kk),f1x2(kk),f1x3(kk),f1x4(kk),
     &          f3d(kk),f3h(kk),f3he3(kk),f3he4(kk),f3neu(kk),f3ntrl(kk)
     &          ,f3t(kk),f3x(kk),f3x1(kk),f3x2(kk),f3x3(kk),f3x4(kk),
     &          r1dd(kk),r1dhe3(kk),r1dt(kk),r3dd(kk),r3dhe3(kk),
     &          r3dt(kk),ye1(kk),ye3(kk),yi1(kk),yi3(kk)
      COMMON /comfus/ deuter,f1d,f1h,f1he3,f1he4,f1neu,f1ntrl,f1t,f1x,
     &                f1x1,f1x2,f1x3,f1x4,f3d,f3h,f3he3,f3he4,f3neu,
     &                f3ntrl,f3t,f3x,f3x1,f3x2,f3x3,f3x4,heliu3,heliu4,
     &                hydrog,xetral,xtrlms,pneut1,pneut3,r1dd,r1dhe3,
     &                r1dt,r3dd,r3dhe3,r3dt,rneut1,rneut3,tinucl,totneu,
     &                tritiu,xmass,xmass1,xmass2,xmass3,xmass4,xtra,
     &                xtra1,xtra2,xtra3,xtra4,xz,xz1,xz2,xz3,xz4,ye1,
     &                ye3,yi1,yi3,yield
c end of thermonuclear reactions
c c3.1 numerical control parameters
      LOGICAL break,nlgoon,nlite
      COMMON /comnc / nceldt,ncondt,nit,nitmax,break,nlgoon,nlite,ak0,
     &                ak1,ak2,ak3,ak4,ak5,ak6,ak7,ak8,ak9,bneum,deltat,
     &                dt2,dt3,dt4,dtemax,dtimax,dtmax,dtmin,dumax,rdt2,
     &                rdt3,rdt4
c end of numerical control parameters
c c3.2 mesh and numerical methods
      DIMENSION ae(kk),ai(kk),be(kk),bi(kk),ce(kk),ci(kk),de(kk),di(kk),
     &          e(kk),f(kk),ge(kk),gi(kk),q2(kk),q4(kk),teite(kk),
     &          tiite(kk),uite(kk)
      COMMON /comnum/ mesh,nj,njm1,njp1,nl,nlm1,nlp1,ae,ai,be,bi,ce,ci,
     &                de,di,e,f,ge,gi,q2,q4,teite,tiite,uite
c end of mesh and numerical methods
c end of basic system parameters
c end of physics control
c c4.1 - administrative variables
c c4.1 - end administrative variables
      DATA iclass,isub/2,4/
      IF(nlomt2(isub))RETURN
c scan over all nl cells and find g for both ions and electrons, but
c notice a(1)=0.0 . ti(0) and te(0) may be arbitrary (also ti(nlp1))
c the arrays a and c refer to level 1 (previous timestep)
c     do 100 l = 1, nl
c        zqv    = (v3(l) - v1(l)) * q2(l)
c        zhi    = ai(l) * (ti1(l) - ti1(l - 1)) - ci(l) *
c    +                  (ti1(l + 1) - ti1(l))
c        zhe    = ae(l) * (te1(l) - te1(l - 1)) - ce(l) *
c    +                  (te1(l + 1) - te1(l))
c do the start
      zqv=(v3(1)-v1(1))*q2(1)
      zhi=ai(1)*ti1(1)-ci(1)*(ti1(2)-ti1(1))
      zhe=ae(1)*te1(1)-ce(1)*(te1(2)-te1(1))
      gi(1)=0.5*dt2*yi1(1)+zhi-zqv
      ge(1)=0.5*dt2*(ye1(1)+xl1(1)+brems1(1))+zhe
c do the end
      zqv=(v3(nl)-v1(nl))*q2(nl)
      zhi=ai(nl)*(ti1(nl)-ti1(nl-1))-ci(nl)*(-ti1(nl))
      zhe=ae(nl)*(te1(nl)-te1(nl-1))-ce(nl)*(-te1(nl))
      gi(nl)=0.5*dt2*yi1(nl)+zhi-zqv
      ge(nl)=0.5*dt2*(ye1(nl)+xl1(nl)+brems1(nl))+zhe
c do the middle region
      i=nl-1
c zero trip ok
      DO 100 l=2,i
       zqv=(v3(l)-v1(l))*q2(l)
       zhi=ai(l)*(ti1(l)-ti1(l-1))-ci(l)*(ti1(l+1)-ti1(l))
       zhe=ae(l)*(te1(l)-te1(l-1))-ce(l)*(te1(l+1)-te1(l))
c eq.67
c ions : no absorption or bremsstrahlung, but neuman pressure
       gi(l)=0.5*dt2*yi1(l)+zhi-zqv
c electrons : no neuman pressure
       ge(l)=0.5*dt2*(ye1(l)+xl1(l)+brems1(l))+zhe
 100  CONTINUE
      RETURN
      END


      SUBROUTINE gotoer(notes,i)
c a call to this routine should be placed after computed goto's so
c aborting the program (as in fortran 66) when the integer value is
c either less than 1 or greater than the number of statement labels.
c in fortran 77 no abort occurs and execution continues with the next
c executable statement after the computed goto.
c smc 24 feb 1985
c c1.1 basic system parameters
      LOGICAL nlend,nlres
      COMMON /combas/ altime,cptime,stime,ndiary,nin,nledge,nonlin,nout,
     &                nprint,npunch,nread,nrec,nresm,nrun,nstep,nlend,
     &                nlres
c end of basic system parameters
      CHARACTER notes*(*)
      WRITE(nout,*)notes,i
      STOP 1600
      END


      SUBROUTINE hcduct
      PARAMETER(kk=1001)
c 2.13 thermal conductivities
c c1.9 development and diagnostic parameters
      DIMENSION nadump(20),npdump(20),nvdump(20)
      LOGICAL nlched,nlhead(9),nlomt1(50),nlomt2(50),nlomt3(50),nlrept
      COMMON /comddp/ maxdum,mxdump,nadump,nclass,npdump,npoint,nsub,
     &                nvdump,nlched,nlhead,nlomt1,nlomt2,nlomt3,nlrept
c end of development and diagnostic parameters
c c2.1 hydrodynamics
      DIMENSION dm(kk),p3(kk),r1(kk),r3(kk),r5(kk),rho1(kk),rho3(kk),
     &          rho5(kk),u2(kk),u4(kk),v1(kk),v3(kk),v5(kk),lstat(kk)
      COMMON /comhyd/ dm,p3,pini,r1,r3,r5,rho1,rho3,rho5,rhoini,rhor,
     &                rini,time,u2,u4,uedge,v1,v3,v5
c end of hydrodynamics
c c2.2 thermodynamics
      DIMENSION ddroe1(kk),ddroe3(kk),ddroi1(kk),ddroi3(kk),ddte1(kk),
     &          ddte3(kk),ddti1(kk),ddti3(kk),xappae(kk),xappai(kk),
     &          pe1(kk),pe3(kk),pi1(kk),pi3(kk),te1(kk),te3(kk),ti1(kk),
     &          ti3(kk)
      COMMON /comth / ddroe1,ddroe3,ddroi1,ddroi3,ddte1,ddte3,ddti1,
     &                ddti3,gammae,gammai,xappae,xappai,pe1,pe3,pi1,pi3,
     &                te1,te3,teini,ti1,ti3,tiini
c end of thermodynamics
c c2.3 ions and electrons
      DIMENSION brems1(kk),brems3(kk),ddz1(kk),ddz3(kk),degen(kk),
     &          effz(kk),eixch2(kk),fidash(kk),fz1(kk),fz3(kk),fzsq1(kk)
     &          ,fzsq3(kk),xlc(kk),xmieff(kk),xne(kk),xni(kk),omega1(kk)
      COMMON /comie / brems1,brems3,ddz1,ddz3,degen,degmax,degmin,effz,
     &                eixch2,fidash,fz1,fz3,fzsq1,fzsq3,xlc,xmieff,xne,
     &                xni,omega1,pmass
c end of ions and electrons
c c2.5 thermonuclear reactions
c end of thermonuclear reactions
c c2.7 physics control
      LOGICAL mlbrms,mlecon,mlfuse,mlicon,mlx,nlabs,nlbrms,nlburn,
     &        nlcri1,nldepo,nlecon,nlfuse,nlicon,nlmove,nlpfe,nlpfi,nlx,
     &        nlhf,nltnl
      COMMON /comcon/ mstep,ncase,ngeom,mlbrms,mlecon,mlfuse,mlicon,mlx,
     &                nlabs,nlbrms,nlburn,nlcri1,nldepo,nlecon,nlfuse,
     &                nlicon,nlmove,nlpfe,nlpfi,nlx,nlhf,nltnl,pmin,
     &                rhomin,temin,timin,umin
c end of physics control
c c3.2 mesh and numerical methods
      DIMENSION ae(kk),ai(kk),be(kk),bi(kk),ce(kk),ci(kk),de(kk),di(kk),
     &          e(kk),f(kk),ge(kk),gi(kk),q2(kk),q4(kk),teite(kk),
     &          tiite(kk),uite(kk)
      COMMON /comnum/ mesh,nj,njm1,njp1,nl,nlm1,nlp1,ae,ai,be,bi,ce,ci,
     &                de,di,e,f,ge,gi,q2,q4,teite,tiite,uite
c end of mesh and numerical methods
c c4.1 - administrative variables
      LOGICAL nldump,nlemp
      DIMENSION piq(302)
      COMMON /comadm/ maxdim,maxrun,ndump,nrep,nldump,nlemp,piq,tstop
c c4.1 - end administrative variables
      DATA iclass,isub/2,13/
      IF(nlomt2(isub))RETURN
c 1. set constants
c the constant of equation 14
c     zi        = 1.084e-11 * (1.0 + piq(26))
      zi=7.46E-12*(1.0+piq(26))
c the constant of equation 15
c     ze        = 4.892e-9 * (1.0 + piq(27))
      ze=1.95E-9*(1.0+piq(27))
c ad1 allan: constant for heat conduction using ad1 allan's formula
c     ze = 1.782e-10*(1.+piq(27))
c the constant of equation 19
      zc=1.38E-23/9.1096E-31
c switch on the limit set by free streaming
      zfac=piq(10)+1.0E-6
c 2. ion thermal conductivity (eq.14)
c is ion heat conduction switched on?
      IF(nlicon.and.mlicon)THEN
       zk1=0.0
       zk2=0.0
       DO 50 l=1,nl
        j=l
        z2=fzsq3(l)
c average of sqrt(m) / z ** 2
        zsqmz2=1.0/sqrt(xmieff(l))*fzsq3(l)
        zk1=zi*ti3(l)*ti3(l)*sqrt(ti3(l))*zsqmz2/(xlc(l)*z2)
c we store the linear average of xappai at l and l - 1 in position j
        xappai(j)=0.5*(zk1+zk2)
        zk2=zk1
 50    CONTINUE
       IF(ncase.eq.2)THEN
c in case 2 we have to define the conductivity at the wall
        zk1=zi*ti3(nlp1)*ti3(nlp1)*sqrt(ti3(nlp1))*zsqmz2/(xlc(nlp1)*z2)
        xappai(njp1)=0.5*(zk1+zk2)
       ENDIF
      ELSE
       CALL resetr(xappai,nlp1,0.0)
      ENDIF
c 3. electron thermal conductivity (eq.15)
c is electron heat conduction switched on?
      IF(.not.(nlecon.and.mlecon))THEN
       CALL resetr(xappae,nlp1,0.0)
       RETURN
      ENDIF
      zk1=0.0
      zk2=0.0
      DO 100 l=1,nl
c average of z ** 2
       z2=fzsq3(l)
       j=l
       zk1=ze*te3(l)*te3(l)*sqrt(te3(l))*fz3(l)/(xlc(l)*z2)
c ad1 allan: heat conduction using ad1 allan's formula
c        zk1 = ze*te3(l)*te3(l)*sqrt(te3(l))/(xlc(l))
c we introduce the spitzer correction factors epsilon and delta
       zk1=zk1*0.095*(fz3(l)+0.24)/(1.0+0.24*fz3(l))
c we store the average at l and l - 1 in j
       xappae(j)=0.5*(zk1+zk2)
       zk2=zk1
c the heat flow is compared with the saturated
c heat flow and the minimum is taken
       IF(l.ne.1)THEN
        zcheck=-xappae(j)*(te3(l)-te3(l-1))/((r3(j+1)-r3(j-1))/2.0)
c the ne values are averaged over the cell volume
        zvoll1=((r3(j-1)+r3(j))/2.0)**(ngeom-1)*(r3(j)-r3(j-1))
        zvoll=((r3(j)+r3(j+1))/2.0)**(ngeom-1)*(r3(j+1)-r3(j))
        znbar=((zvoll1*xne(l-1))+zvoll*xne(l))/(zvoll1+zvoll)
c the maximum temperature is taken
        zte3j=te3(l)
        IF(te3(l-1).gt.te3(l))zte3j=te3(l-1)
        zsatj=(1.38E-23/zfac)*sqrt(zc*zte3j)*zte3j
        zsatj=zsatj*znbar
        zhtrat=abs(zcheck)/zsatj
        IF(zhtrat.gt.1.0)xappae(j)=xappae(j)/zhtrat
       ENDIF
 100  CONTINUE
      IF(ncase.ne.2)RETURN
c in case 2 we need to define the conductivity at the wall, but
c we do not impose the free streaming limit. further the composition
c coefficients at the wall are those of the plasma inside it.
      zk1=ze*te3(nlp1)*te3(nlp1)*sqrt(te3(nlp1))*fz3(nl)/(xlc(nlp1)*z2)
c ad1 allan: heat conduction using ad1 allan's formula
c     zk1 = ze*te3(nlp1)*te3(nlp1)*sqrt(te3(nlp1))
c    &      /(xlc(nlp1))
      xappae(njp1)=0.5*(zk1+zk2)
      RETURN
      END


      SUBROUTINE hdcopy(k)
c this is a routine to draw smop graphs
c c2.1 hydrodynamics
c end of hydrodynamics
c c2.2 thermodynamics
c end of thermodynamics
c c2.4 laser variables
c end of laser variables
c c3.2 mesh and numerical methods
c end of mesh and numerical methods
c c4.1 - administrative variables
c c4.1 - end administrative variables
c c1.1 basic system parameters
c end of basic system parameters
c zx and zy are the points to be plotted on the graph
c zlbx, zlbd, zlbe, zlbi are used only for labelling
c the routine is called with k = 1 to initialise graphics
c                            k = 2 to plot a frame
c                            k = 3 to close graphics
      RETURN
c     goto (1, 2, 3), k
c     call gotoer(' medusa: hdcopy: goto error, k is', k)
c initialise graphics
c   1 call sp$set(1, 10000000)
c     call frspec
c set maxinum and mininum quantities
c to be plotted on x and y axis
c     zmaxx     = 2.0 * rini
c     zmaxy     = 9.0
c     zminy     = -1.0
c     call limits(-10.0, 10.0, 512.0, 512.0, 1.0, 1.0)
c save the graphic commands for writing labels
c and drawing axises on every frame
c     zlbx(1)   = 0.000015
c     zlbx(2)   = 0.00002
c     zlbd(1)   = 8.0
c     zlbd(2)   = 8.0
c     zlbe(1)   = 7.5
c     zlbe(2)   = 7.5
c     zlbi(1)   = 7.0
c     zlbi(2)   = 7.0
c     call draxop(20, .false.)
c     call draxop(12, 3)
c     call draxop(19, .false.)
c advance one frame
c   2 continue
c copy quantities to be plotted to arrays zx and zy
c     do 10 l = 1, mesh
c        zx(l)  = r3(l)
c        zy(l)  = log10(rho3(l))
c  10 continue
c draw axes and plot  the points with a continuous line for
c log rho vs. radius
c     call advflm
c     call draxes(0.0, zminy, zmaxx, zmaxy, 'radius!',
c    +            'log rho  log te  log  ti!', 'medusa!', '!')
c     call color(1.0, 2.0, 200.0)
c     call drline(zx, zy, mesh)
c copy quantities to be plotted for electron temperature
c     do 15 l = 1, mesh
c        zy(l)  = log10(te3(l))
c  15 continue
c now draw the graph of electron temprature vs radius
c     call color(4.0, 10.0, 240.0)
c     call drdash(zx, zy, mesh, 0.3, 0.3, 0.6)
c copy quantities to be plotted for ion temperature
c     do 20 l = 1, mesh
c        zy(l)  = log10(ti3(l))
c  20 continue
c now draw the graph of ion temperature vs radius
c     call color(2.0, 1.0, 200.0)
c     call drdash(zx, zy, mesh, 2.0, 0.3, 0.6)
c write time and labels above the curves
c     call color(8.0, 2.0, 200.0)
c     call setxy(0.000005, 8.5)
c     call htext(7, 'time = ')
c     call typnmb(time, 0.0, 3.0)
c the time is drawn on all pictures but the key to
c the different lines is drawn only on the first picture
c     if (nstep.gt.300) then
c       return
c     end if
c     call color(1.0, 2.0, 200.0)
c     call setxy(0.000005, 8.0)
c     call htext(8, 'density ')
c     call drline(zlbx, zlbd, 2)
c     call color(4.0, 10.0, 240.0)
c     call setxy(0.000005,7.5)
c     call htext(7, 'temp.e.')
c     call drdash(zlbx, zlbe, 2, 0.3, 0.3, 0.6)
c     call color(2.0, 1.0, 200.0)
c     call setxy(0.000005, 7.0)
c     call htext(7, 'temp.i.')
c     call drdash(zlbx, zlbi, 2, 2.0, 0.3, 0.6)
c     call color(8.0, 2.0, 200.0)
c     if (time .gt. 0.0) then
c       return
c     end if
c  22 continue
c     return
c close the graphics
c   3 call endspr
c     return
      END


      SUBROUTINE hdump(k)
c dummy routine
      RETURN
      END


      SUBROUTINE hotran
c 3.2 mesh and numerical methods
c c3.2 mesh and numerical methods
c end of mesh and numerical methods
c dummy routine
      RETURN
      END


      SUBROUTINE iarray(kname,ka,kdim)
c u.9 print name and values of integer array
c c1.1 basic system parameters
      LOGICAL nlend,nlres
      COMMON /combas/ altime,cptime,stime,ndiary,nin,nledge,nonlin,nout,
     &                nprint,npunch,nread,nrec,nresm,nrun,nstep,nlend,
     &                nlres
c end of basic system parameters
      CHARACTER kname*(*)
      DIMENSION ka(kdim)
      CALL blines(1)
      WRITE(nout,10100)kname
      CALL blines(1)
      WRITE(nout,10200)ka
      RETURN
10100 FORMAT(4x,a)
10200 FORMAT((6x,10(i12)))
      END


      SUBROUTINE inital
      PARAMETER(kk=1001,underf=1.E-30)
c 1.6 define physical conditions
c c1.9 development and diagnostic parameters
      DIMENSION nadump(20),npdump(20),nvdump(20)
      LOGICAL nlched,nlhead(9),nlomt1(50),nlomt2(50),nlomt3(50),nlrept
      COMMON /comddp/ maxdum,mxdump,nadump,nclass,npdump,npoint,nsub,
     &                nvdump,nlched,nlhead,nlomt1,nlomt2,nlomt3,nlrept
c end of development and diagnostic parameters
c c2.1 hydrodynamics
      DIMENSION dm(kk),p3(kk),r1(kk),r3(kk),r5(kk),rho1(kk),rho3(kk),
     &          rho5(kk),u2(kk),u4(kk),v1(kk),v3(kk),v5(kk),lstat(kk)
      COMMON /comhyd/ dm,p3,pini,r1,r3,r5,rho1,rho3,rho5,rhoini,rhor,
     &                rini,time,u2,u4,uedge,v1,v3,v5
c end of hydrodynamics
c c2.2 thermodynamics
      DIMENSION ddroe1(kk),ddroe3(kk),ddroi1(kk),ddroi3(kk),ddte1(kk),
     &          ddte3(kk),ddti1(kk),ddti3(kk),xappae(kk),xappai(kk),
     &          pe1(kk),pe3(kk),pi1(kk),pi3(kk),te1(kk),te3(kk),ti1(kk),
     &          ti3(kk)
      COMMON /comth / ddroe1,ddroe3,ddroi1,ddroi3,ddte1,ddte3,ddti1,
     &                ddti3,gammae,gammai,xappae,xappai,pe1,pe3,pi1,pi3,
     &                te1,te3,teini,ti1,ti3,tiini
c end of thermodynamics
c c2.3 ions and electrons
      DIMENSION brems1(kk),brems3(kk),ddz1(kk),ddz3(kk),degen(kk),
     &          effz(kk),eixch2(kk),fidash(kk),fz1(kk),fz3(kk),fzsq1(kk)
     &          ,fzsq3(kk),xlc(kk),xmieff(kk),xne(kk),xni(kk),omega1(kk)
      COMMON /comie / brems1,brems3,ddz1,ddz3,degen,degmax,degmin,effz,
     &                eixch2,fidash,fz1,fz3,fzsq1,fzsq3,xlc,xmieff,xne,
     &                xni,omega1,pmass
c end of ions and electrons
c c2.5 thermonuclear reactions
      DIMENSION f1d(kk),f1h(kk),f1he3(kk),f1he4(kk),f1neu(kk),f1ntrl(kk)
     &          ,f1t(kk),f1x(kk),f1x1(kk),f1x2(kk),f1x3(kk),f1x4(kk),
     &          f3d(kk),f3h(kk),f3he3(kk),f3he4(kk),f3neu(kk),f3ntrl(kk)
     &          ,f3t(kk),f3x(kk),f3x1(kk),f3x2(kk),f3x3(kk),f3x4(kk),
     &          r1dd(kk),r1dhe3(kk),r1dt(kk),r3dd(kk),r3dhe3(kk),
     &          r3dt(kk),ye1(kk),ye3(kk),yi1(kk),yi3(kk)
      COMMON /comfus/ deuter,f1d,f1h,f1he3,f1he4,f1neu,f1ntrl,f1t,f1x,
     &                f1x1,f1x2,f1x3,f1x4,f3d,f3h,f3he3,f3he4,f3neu,
     &                f3ntrl,f3t,f3x,f3x1,f3x2,f3x3,f3x4,heliu3,heliu4,
     &                hydrog,xetral,xtrlms,pneut1,pneut3,r1dd,r1dhe3,
     &                r1dt,r3dd,r3dhe3,r3dt,rneut1,rneut3,tinucl,totneu,
     &                tritiu,xmass,xmass1,xmass2,xmass3,xmass4,xtra,
     &                xtra1,xtra2,xtra3,xtra4,xz,xz1,xz2,xz3,xz4,ye1,
     &                ye3,yi1,yi3,yield
c end of thermonuclear reactions
c c2.7 physics control
      LOGICAL mlbrms,mlecon,mlfuse,mlicon,mlx,nlabs,nlbrms,nlburn,
     &        nlcri1,nldepo,nlecon,nlfuse,nlicon,nlmove,nlpfe,nlpfi,nlx,
     &        nlhf,nltnl
      COMMON /comcon/ mstep,ncase,ngeom,mlbrms,mlecon,mlfuse,mlicon,mlx,
     &                nlabs,nlbrms,nlburn,nlcri1,nldepo,nlecon,nlfuse,
     &                nlicon,nlmove,nlpfe,nlpfi,nlx,nlhf,nltnl,pmin,
     &                rhomin,temin,timin,umin
c end of physics control
c c3.2 mesh and numerical methods
      DIMENSION ae(kk),ai(kk),be(kk),bi(kk),ce(kk),ci(kk),de(kk),di(kk),
     &          e(kk),f(kk),ge(kk),gi(kk),q2(kk),q4(kk),teite(kk),
     &          tiite(kk),uite(kk)
      COMMON /comnum/ mesh,nj,njm1,njp1,nl,nlm1,nlp1,ae,ai,be,bi,ce,ci,
     &                de,di,e,f,ge,gi,q2,q4,teite,tiite,uite
c end of mesh and numerical methods
c c4.1 - administrative variables
      LOGICAL nldump,nlemp
      DIMENSION piq(302)
      COMMON /comadm/ maxdim,maxrun,ndump,nrep,nldump,nlemp,piq,tstop
c c4.1 - end administrative variables
c x-ray laser common blocks
      COMMON /atom  / pu(10,kk),p(10,kk),pz(10,kk),zs(10,kk),zst(kk),
     &                zstz(kk),zt(kk),a(kk),fn(10),deg(10),qs(10,kk),
     &                zstmi(kk),zstu(kk),engu(10,kk),eng(10,kk),
     &                rdrhe(kk),zstin(kk),nmax
      COMMON /radflg/ yamda,xhnu,xihnu,dihnu,xplas1(kk),xrbf(kk),
     &                atrbf(kk),corhf1(kk),corhf(kk),iradf1,iradf2
c 1. hydrodynamic variables at level 1 & 2
c tracer variables are set here. this enables us to modify the
c initial conditions in subprogram source
      nclass=1
      nsub=6
c 1.1 lagrangian coordinates. all geometries.
      zdr=rini/real(mesh)
      DO 100 j=1,njp1
       r1(j)=real(j-1)*zdr
 100  CONTINUE
      npoint=1
      CALL source
c 1.2 the initial density profile
      DO 200 l=1,nl
       rho1(l)=rhoini
 200  CONTINUE
      npoint=2
      CALL source
c 1.3 the masses of each cell
      npoint=3
      IF(ngeom.eq.1)THEN
      ELSEIF(ngeom.eq.2)THEN
c ngeom = 2 cylinder geometry: a section of unit height & angle
       DO 250 l=1,nl
        j=l
        dm(l)=(r1(j+1)*r1(j+1)-r1(j)*r1(j))*rho1(l)*0.5
 250   CONTINUE
       GOTO 500
      ELSEIF(ngeom.eq.3)THEN
c ngeom = 3 spherical geometry: a steradian
       DO 300 l=1,nl
        j=l
        dm(l)=(r1(j+1)*r1(j+1)*r1(j+1)-r1(j)*r1(j)*r1(j))/3.0*rho1(l)
 300   CONTINUE
       GOTO 500
      ELSE
       CALL gotoer(' MEDUSA: inital: goto error, ngeom is',ngeom)
      ENDIF
c ngeom = 1 cartesian geometry: a slab of unit cross section
      DO 400 l=1,nl
       j=l
       dm(l)=(r1(j+1)-r1(j))*rho1(l)
 400  CONTINUE
c 1.4 initial velocities
c these are all left at zero
 500  npoint=4
      CALL source
c 1.5 specific volumes (eq.6)
      DO 600 l=1,nl
       v1(l)=1.0/rho1(l)
 600  CONTINUE
      npoint=5
      CALL source
c 1.6 component fractions (eq.1)
      IF(abs(deuter).gt.underf)CALL resetr(f1d,mesh,deuter)
      IF(abs(hydrog).gt.underf)CALL resetr(f1h,mesh,hydrog)
      IF(abs(heliu3).gt.underf)CALL resetr(f1he3,mesh,heliu3)
      IF(abs(heliu4).gt.underf)CALL resetr(f1he4,mesh,heliu4)
      IF(abs(xetral).gt.underf)CALL resetr(f1ntrl,mesh,xetral)
      IF(abs(tritiu).gt.underf)CALL resetr(f1t,mesh,tritiu)
      IF(abs(xtra).gt.underf)CALL resetr(f1x,mesh,xtra)
      npoint=6
      CALL source
c 1.7 average mass and charge numbers (eqs.3&4)
c the standard ion masses
      zmd=2.0
      zmh=1.0
      zmhe3=3.0
      zmhe4=4.0
      zmt=3.0
c charge numbers
      zcd=1.0
      zch=1.0
      zche3=2.0
      zche4=2.0
      zct=1.0
c charge numbers squared
      zcd2=1.0
      zch2=1.0
      zche32=4.0
      zche42=4.0
      zct2=1.0
      DO 700 l=1,nl
       xmieff(l)=f1d(l)*zmd+f1h(l)*zmh+f1he3(l)*zmhe3+f1he4(l)
     &           *zmhe4+f1t(l)*zmt+f1x(l)*xmass+f1ntrl(l)*xtrlms+f1x1(l)
     &           *xmass1+f1x2(l)*xmass2+f1x3(l)*xmass3+f1x4(l)*xmass4
       effz(l)=f1d(l)*zcd+f1h(l)*zch+f1he3(l)*zche3+f1he4(l)
     &         *zche4+f1t(l)*zct+f1x(l)*xz+f1x1(l)*xz1+f1x2(l)
     &         *xz2+f1x3(l)*xz3+f1x4(l)*xz4
       fz1(l)=f1d(l)*zcd+f1h(l)*zch+f1he3(l)*zche3+f1he4(l)*zche4+f1t(l)
     &        *zct+f1x(l)*xz+f1x1(l)*xz1+f1x2(l)*xz2+f1x3(l)*xz3+f1x4(l)
     &        *xz4
       fzsq1(l)=f1d(l)*zcd2+f1h(l)*zch2+f1he3(l)*zche32+f1he4(l)
     &          *zche42+f1t(l)*zct2+f1x(l)*xz*xz+f1x1(l)*xz1*xz1+f1x2(l)
     &          *xz2*xz2+f1x3(l)*xz3*xz3+f1x4(l)*xz4*xz4
 700  CONTINUE
c 1.8 particle densities (eqs.5&6)
      DO 800 l=1,nl
       xni(l)=rho1(l)/(xmieff(l)*pmass)
       xne(l)=effz(l)*xni(l)
 800  CONTINUE
c 1.9 copy to level 3 where necessary
      CALL copyr(r1,1,r3,1,njp1)
      CALL copyr(v1,1,v3,1,nl)
      CALL copyr(rho1,1,rho3,1,nl)
      CALL copyr(f1d,1,f3d,1,nl)
      CALL copyr(f1h,1,f3h,1,nl)
      CALL copyr(f1he3,1,f3he3,1,nl)
      CALL copyr(f1he4,1,f3he4,1,nl)
      CALL copyr(f1ntrl,1,f3ntrl,1,nl)
      CALL copyr(f1t,1,f3t,1,nl)
      CALL copyr(f1x,1,f3x,1,nl)
      CALL copyr(f1x1,1,f3x1,1,nl)
      CALL copyr(f1x2,1,f3x2,1,nl)
      CALL copyr(f1x3,1,f3x3,1,nl)
      CALL copyr(f1x4,1,f3x4,1,nl)
      CALL copyr(fz1,1,fz3,1,nl)
      CALL copyr(fzsq1,1,fzsq3,1,nl)
c 2. thermodynamic variables
c 2.1 standard temperatures
c ad1
      njg=nl-nint(piq(13))-nint(piq(63))
      njp=nl-nint(piq(63))
      IF(abs(piq(86)).lt.underf)piq(86)=teini
      IF(abs(piq(89)).lt.underf)piq(89)=teini
      IF(abs(piq(87)).lt.underf)piq(87)=tiini
      IF(abs(piq(90)).lt.underf)piq(90)=tiini
      IF(abs(piq(85)).lt.underf)piq(85)=1.000
      IF(abs(piq(88)).lt.underf)piq(88)=1.000
      IF(abs(piq(91)).lt.underf)piq(91)=1.000
      DO 900 l=1,njg
       ti1(l)=tiini
       ti3(l)=tiini
       te1(l)=teini
       te3(l)=teini
       zstmi(l)=piq(85)
       zstin(l)=piq(92)
 900  CONTINUE
      DO 1000 l=njg+1,njp
       ti1(l)=piq(87)
       ti3(l)=piq(87)
       te1(l)=piq(86)
       te3(l)=piq(86)
       zstmi(l)=piq(88)
       zstin(l)=piq(93)
 1000 CONTINUE
      DO 1100 l=njp+1,nl
       zstin(l)=piq(94)
       zstmi(l)=piq(91)
       ti1(l)=piq(90)
       ti3(l)=piq(90)
       te1(l)=piq(89)
       te3(l)=piq(89)
 1100 CONTINUE
c hf correction
      DO 1200 l=1,nlp1
       corhf(l)=1.
       corhf1(l)=1.
 1200 CONTINUE
      npoint=7
      CALL source
c 2.2 corresponding pressures
c these are found from the equations of state, stored at level 3, so
c we make a copy at level 1
c ad1 following lines added to get proper initial energy input e/ad1t
      CALL radran
c partial derivatives wrt z fom ionbal
      CALL ionbal
c ionbal provides correct values of fz1 and fzsq1 to overide
c values calculated earlier
      CALL copyr(fz3,1,fz1,1,nlp1)
      CALL copyr(fzsq3,1,fzsq1,1,nlp1)
c ad1
      CALL statei(2)
      CALL statee(2)
      CALL copyr(pi3,1,pi1,1,nl)
      CALL copyr(pe3,1,pe1,1,nl)
c form the total pressure
      CALL formp
      npoint=8
      CALL source
c we have now defined the initial conditions except at the boundary
c cell point  l = nl + 1 = nlp1. also to start the calculation it is
c necessary to work out a few more quantities. this is done in
c subrout    start.
c 3. initial boundary values
c which case are we treating?
      IF(ncase.eq.1)THEN
      ELSEIF(ncase.eq.2)THEN
c 3.2 temperature determined
c we set the boundary temperatures at nlp1 equal to those at nl
       ti1(nlp1)=ti1(nl)
       ti3(nlp1)=ti3(nl)
       te1(nlp1)=te1(nl)
       te3(nlp1)=te3(nl)
      ELSEIF(ncase.eq.3)THEN
c 3.3 pressure determined
c the initial boundary pressure is set equal to the pressure just
c inside the boundary and stored in pini
       pini=p3(nl)
      ELSEIF(ncase.ne.4)THEN
c 3.1 thermally insulated wall
c no quantities defined
       CALL gotoer(' MEDUSA: inital: goto error, ncase is',ncase)
      ENDIF
c 3.4 velocity determined
c the boundary velocity is initially zero
c 4. initial conditions can be changed here
      npoint=9
      CALL source
      RETURN
      END


      SUBROUTINE ionbal
      PARAMETER(kk=1001,underf=1.E-30)
c subroutine to calculate saha ionisation balance
c of a mixture of materials
c c2.2 thermodynamics
      DIMENSION ddroe1(kk),ddroe3(kk),ddroi1(kk),ddroi3(kk),ddte1(kk),
     &          ddte3(kk),ddti1(kk),ddti3(kk),xappae(kk),xappai(kk),
     &          pe1(kk),pe3(kk),pi1(kk),pi3(kk),te1(kk),te3(kk),ti1(kk),
     &          ti3(kk)
      COMMON /comth / ddroe1,ddroe3,ddroi1,ddroi3,ddte1,ddte3,ddti1,
     &                ddti3,gammae,gammai,xappae,xappai,pe1,pe3,pi1,pi3,
     &                te1,te3,teini,ti1,ti3,tiini
c end of thermodynamics
c c4.1 - administrative variables
      LOGICAL nldump,nlemp
      DIMENSION piq(302)
      COMMON /comadm/ maxdim,maxrun,ndump,nrep,nldump,nlemp,piq,tstop
c c4.1 - end administrative variables
c c2.3 ions and electrons
      DIMENSION brems1(kk),brems3(kk),ddz1(kk),ddz3(kk),degen(kk),
     &          effz(kk),eixch2(kk),fidash(kk),fz1(kk),fz3(kk),fzsq1(kk)
     &          ,fzsq3(kk),xlc(kk),xmieff(kk),xne(kk),xni(kk),omega1(kk)
      COMMON /comie / brems1,brems3,ddz1,ddz3,degen,degmax,degmin,effz,
     &                eixch2,fidash,fz1,fz3,fzsq1,fzsq3,xlc,xmieff,xne,
     &                xni,omega1,pmass
c end of ions and electrons
c c3.1 numerical control parameters
c end of numerical control parameters
c c3.2 mesh and numerical methods
      DIMENSION ae(kk),ai(kk),be(kk),bi(kk),ce(kk),ci(kk),de(kk),di(kk),
     &          e(kk),f(kk),ge(kk),gi(kk),q2(kk),q4(kk),teite(kk),
     &          tiite(kk),uite(kk)
      COMMON /comnum/ mesh,nj,njm1,njp1,nl,nlm1,nlp1,ae,ai,be,bi,ce,ci,
     &                de,di,e,f,ge,gi,q2,q4,teite,tiite,uite
c end of mesh and numerical methods
c c2.1 hydrodynamics
      DIMENSION dm(kk),p3(kk),r1(kk),r3(kk),r5(kk),rho1(kk),rho3(kk),
     &          rho5(kk),u2(kk),u4(kk),v1(kk),v3(kk),v5(kk),lstat(kk)
      COMMON /comhyd/ dm,p3,pini,r1,r3,r5,rho1,rho3,rho5,rhoini,rhor,
     &                rini,time,u2,u4,uedge,v1,v3,v5
c end of hydrodynamics
c c2.5 thermonuclear reactions
      DIMENSION f1d(kk),f1h(kk),f1he3(kk),f1he4(kk),f1neu(kk),f1ntrl(kk)
     &          ,f1t(kk),f1x(kk),f1x1(kk),f1x2(kk),f1x3(kk),f1x4(kk),
     &          f3d(kk),f3h(kk),f3he3(kk),f3he4(kk),f3neu(kk),f3ntrl(kk)
     &          ,f3t(kk),f3x(kk),f3x1(kk),f3x2(kk),f3x3(kk),f3x4(kk),
     &          r1dd(kk),r1dhe3(kk),r1dt(kk),r3dd(kk),r3dhe3(kk),
     &          r3dt(kk),ye1(kk),ye3(kk),yi1(kk),yi3(kk)
      COMMON /comfus/ deuter,f1d,f1h,f1he3,f1he4,f1neu,f1ntrl,f1t,f1x,
     &                f1x1,f1x2,f1x3,f1x4,f3d,f3h,f3he3,f3he4,f3neu,
     &                f3ntrl,f3t,f3x,f3x1,f3x2,f3x3,f3x4,heliu3,heliu4,
     &                hydrog,xetral,xtrlms,pneut1,pneut3,r1dd,r1dhe3,
     &                r1dt,r3dd,r3dhe3,r3dt,rneut1,rneut3,tinucl,totneu,
     &                tritiu,xmass,xmass1,xmass2,xmass3,xmass4,xtra,
     &                xtra1,xtra2,xtra3,xtra4,xz,xz1,xz2,xz3,xz4,ye1,
     &                ye3,yi1,yi3,yield
c end of thermonuclear reactions
c c2.7 physics control
c end of physics control
      COMMON /atom  / pu(10,kk),p(10,kk),pz(10,kk),zs(10,kk),zst(kk),
     &                zstz(kk),zt(kk),a(kk),fn(10),deg(10),qs(10,kk),
     &                zstmi(kk),zstu(kk),engu(10,kk),eng(10,kk),
     &                rdrhe(kk),zstin(kk),nmax
c x-ray laser common blocks
      zgkm=8255.17/(gammae-1.)
c saha=0 take the plasma to be always fully ionised
      IF(abs(piq(58)).lt.underf)THEN
       DO 50 l=1,nl
        fz3(l)=fz1(l)
        fzsq3(l)=fzsq1(l)
        ddz3(l)=0.0
        fidash(l)=0.0
 50    CONTINUE
       RETURN
      ENDIF
c take the plasma to have the degree of ionisation from the nlte
c package
c saha=2  then use ionbal routine (sjr) for xrl calculations
      IF(abs(piq(58)-2.0).lt.underf)THEN
       DO 100 l=1,nl
c ad1
        IF(zst(l).lt.zstmi(l).or.zst(l).gt.zt(l))zst(l)=zstmi(l)
        fz3(l)=zst(l)
        fzsq3(l)=fz3(l)*fz3(l)
        ddz3(l)=0.0
        fidash(l)=0.0
 100   CONTINUE
       RETURN
      ENDIF
c saha=1   old ionbal routines invoked
      zpmass=1.6726E-27
      DO 300 l=1,nl
       fz3(l)=0.0
       fzsq3(l)=0.0
       ddz3(l)=0.0
       fidash(l)=0.0
       zte=te3(l)/11605.0
       zne=xne(l)
       IF(zte.lt.5.0)zte=5.0
       zfh=f3h(l)+f3d(l)+f3t(l)
       IF(abs(zfh).gt.underf)THEN
        z=1.0
        CALL saha(z,zne,zte,zbar,nite,zphi,zchi,xmieff(l))
        fz3(l)=fz3(l)+zfh*zbar
        fzsq3(l)=fzsq3(l)+zfh*zbar*zbar
        ddz3(l)=ddz3(l)+zfh*zphi
        fidash(l)=fidash(l)+zfh*zchi
       ENDIF
       zfhe=f3he3(l)+f3he4(l)
       IF(abs(zfhe).gt.underf)THEN
        z=2.0
        CALL saha(z,zne,zte,zbar,nite,zphi,zchi,xmieff(l))
        fz3(l)=fz3(l)+zfhe*zbar
        fzsq3(l)=fzsq3(l)+zfhe*zbar*zbar
        ddz3(l)=ddz3(l)+zfhe*zphi
        fidash(l)=fidash(l)+zfhe*zchi
       ENDIF
       IF(abs(f3x(l)).gt.underf)THEN
        IF(xz.gt.18.0)GOTO 150
        z=xz
        CALL saha(z,zne,zte,zbar,nite,zphi,zchi,xmieff(l))
        fz3(l)=fz3(l)+f3x(l)*zbar
        fzsq3(l)=fzsq3(l)+f3x(l)*zbar*zbar
        ddz3(l)=ddz3(l)+f3x(l)*zphi
        fidash(l)=fidash(l)+f3x(l)*zchi
       ENDIF
       IF(abs(f3x1(l)).gt.underf)THEN
        IF(xz1.gt.18.0)GOTO 150
        CALL saha(xz1,zne,zte,zbar,nite,zphi,zchi,xmieff(l))
        fz3(l)=fz3(l)+f3x1(l)*zbar
        fzsq3(l)=fzsq3(l)+f3x1(l)*zbar*zbar
        ddz3(l)=ddz3(l)+f3x1(l)*zphi
        fidash(l)=fidash(l)+f3x1(l)*zchi
       ENDIF
       IF(abs(f3x2(l)).gt.underf)THEN
        IF(xz2.gt.18.0)GOTO 150
        CALL saha(xz2,zne,zte,zbar,nite,zphi,zchi,xmieff(l))
        fz3(l)=fz3(l)+f3x2(l)*zbar
        fzsq3(l)=fzsq3(l)+f3x2(l)*zbar*zbar
        ddz3(l)=ddz3(l)+f3x2(l)*zphi
        fidash(l)=fidash(l)+f3x2(l)*zchi
       ENDIF
       IF(abs(f3x3(l)).gt.underf)THEN
        IF(xz3.gt.18.0)GOTO 150
        CALL saha(xz3,zne,zte,zbar,nite,zphi,zchi,xmieff(l))
        fz3(l)=fz3(l)+f3x3(l)*zbar
        fzsq3(l)=fzsq3(l)+f3x3(l)*zbar*zbar
        ddz3(l)=ddz3(l)+f3x3(l)*zphi
        fidash(l)=fidash(l)+f3x3(l)*zchi
       ENDIF
       IF(abs(f3x4(l)).lt.underf)GOTO 200
       IF(xz4.le.18.0)THEN
        CALL saha(xz4,zne,zte,zbar,nite,zphi,zchi,xmieff(l))
        fz3(l)=fz3(l)+f3x4(l)*zbar
        fzsq3(l)=fzsq3(l)+f3x4(l)*zbar*zbar
        ddz3(l)=ddz3(l)+f3x4(l)*zphi
        fidash(l)=fidash(l)+f3x4(l)*zchi
        GOTO 200
       ENDIF
c for z gt 18 we calculate fz3 from the thomas fermi pressure
c150     fz3(l) = pe3(l)*xmieff(l)/(8255.17*rho3(l)*te3(l))
c  ad1    use3 dick more's tf zstar
 150   rho=rho3(l)/1000.
       theta=te3(l)/11605.
       CALL zbar1(rho,theta,effz(l),xmieff(l),zstar)
       fz3(l)=zstar
       IF(fz3(l).lt.1.0)fz3(l)=1.0
       fzsq3(l)=fz3(l)*fz3(l)
c because of convergence problems in this version of ionbal
c we put ddz = 0.0
 200   ddz3(l)=0.0
       fidash(l)=0.0
c ad1 this is actually because the ionisation energy is already included
c in ddroe3 and ddte3 in  statee
 300  CONTINUE
      RETURN
      END


      SUBROUTINE zbar1(rho,theta,z,a,zstar)
c
c this function calculates an approximation to the thomas-fermi
c lte degree of ionisation
c
c form tzero,r,tf
c
      tev=theta
      fthrds=1.3333333333333333
      tzero=tev/(z**fthrds)
      r=rho/(z*a)
      tf=tzero/(1.0+tzero)
c
c setup constants
c
      a1=0.003323467
      a2=0.97183224
      a3=9.26148E-05
      a4=3.1016524
      b0=-1.762999
      b1=1.4317567
      b2=0.31546338
      c1=-0.36666667
      c2=0.98333333
      alpha=14.3139316
      beta=0.66240046
c
c calculate a,b and c
c
      aa=(a1*(tzero**a2))+(a3*(tzero**a4))
      b=-exp(b0+(b1*tf)+(b2*(tf**7)))
      c=c2+(c1*tf)
c
c calculate q1 and thereby q
c
      q1=aa*(r**b)
      q=(r**c)+(q1**c)
      cm=1.0/c
      q=q**cm
c
c calculate x
c
      x=alpha*(q**beta)
c
c calculate zstar
c
      f=x/(1.0+x+(sqrt(1.0+(2.0*x))))
      zstar=f*z
      RETURN
      END


      SUBROUTINE ivar(kname,kvalue)
c u.5 print name and value of integer variable
c c1.1 basic system parameters
      LOGICAL nlend,nlres
      COMMON /combas/ altime,cptime,stime,ndiary,nin,nledge,nonlin,nout,
     &                nprint,npunch,nread,nrec,nresm,nrun,nstep,nlend,
     &                nlres
c end of basic system parameters
      CHARACTER kname*(*)
      WRITE(nout,10100)kname,kvalue
      RETURN
10100 FORMAT(4x,a,' = ',i12)
      END


      SUBROUTINE jobtim(ptime)
c u.17 fetch allocated jobtime (secs)
c cpulft(0)     returns the amount of cpu time (in minutes) remaining
c               before reaching time limit.
c cpulft(1)     returns the amount of cpu time (in minutes) used so far.
c     second(real) is cray specific - used for timing
c     ptime     = second(0) + second(1)
      RETURN
      END


      SUBROUTINE labrun
c 1.1 label the run
c c1.1 basic system parameters
      CHARACTER label1*80,label2*80,label3*80,label4*80,label5*80,
     &          label6*80,label7*80,label8*80
      LOGICAL nlend,nlres
      COMMON /comlab/ label1,label2,label3,label4,label5,label6,label7,
     &                label8
      COMMON /combas/ altime,cptime,stime,ndiary,nin,nledge,nonlin,nout,
     &                nprint,npunch,nread,nrec,nresm,nrun,nstep,nlend,
     &                nlres
c end of basic system parameters
c read in the date and descriptor sentences
      READ(nread,10300,end=100)label1
 100  READ(nread,10300,end=200)label2
 200  READ(nread,10300,end=300)label3
 300  READ(nread,10300,end=400)label4
c write headings
 400  CALL blines(2)
      WRITE(nout,10100)
      WRITE(nout,10200)
      CALL blines(2)
      WRITE(nout,10400)label1
      CALL blines(2)
      WRITE(nout,10400)label2
      CALL blines(2)
      WRITE(nout,10400)label3
      CALL blines(2)
      WRITE(nout,10400)label4
c write the run descriptor sentences
c     write(ndiary, 9905) label1
      RETURN
10100 FORMAT(50x,'p r o g r a m   M E D U S A vs. 103')
10200 FORMAT(48x,38('+'))
10300 FORMAT(a)
10400 FORMAT(26x,a)
10500 FORMAT(4x,a)
      END


      SUBROUTINE larray(kname,kla,kdim)
c u.18 print name and values of logical array
c c1.1 basic system parameters
      LOGICAL nlend,nlres
      COMMON /combas/ altime,cptime,stime,ndiary,nin,nledge,nonlin,nout,
     &                nprint,npunch,nread,nrec,nresm,nrun,nstep,nlend,
     &                nlres
c end of basic system parameters
      CHARACTER kname*(*)
      LOGICAL kla(kdim)
      CALL blines(1)
      WRITE(nout,10100)kname
      CALL blines(1)
      WRITE(nout,10200)kla
      RETURN
10100 FORMAT(4x,a)
10200 FORMAT((6x,10(l12)))
      END


      SUBROUTINE laser(tnow)
      PARAMETER(kk=1001,underf=1.E-30)
c 2.3 produces the laser power
c c1.1 basic system parameters
      LOGICAL nlend,nlres
      COMMON /combas/ altime,cptime,stime,ndiary,nin,nledge,nonlin,nout,
     &                nprint,npunch,nread,nrec,nresm,nrun,nstep,nlend,
     &                nlres
c end of basic system parameters
c c3.1 numerical control parameters
      LOGICAL break,nlgoon,nlite
      COMMON /comnc / nceldt,ncondt,nit,nitmax,break,nlgoon,nlite,ak0,
     &                ak1,ak2,ak3,ak4,ak5,ak6,ak7,ak8,ak9,bneum,deltat,
     &                dt2,dt3,dt4,dtemax,dtimax,dtmax,dtmin,dumax,rdt2,
     &                rdt3,rdt4
c end of numerical control parameters
c c1.9 development and diagnostic parameters
      DIMENSION nadump(20),npdump(20),nvdump(20)
      LOGICAL nlched,nlhead(9),nlomt1(50),nlomt2(50),nlomt3(50),nlrept
      COMMON /comddp/ maxdum,mxdump,nadump,nclass,npdump,npoint,nsub,
     &                nvdump,nlched,nlhead,nlomt1,nlomt2,nlomt3,nlrept
c end of development and diagnostic parameters
c c2.1 hydrodynamics
c end of hydrodynamics
c c2.4 laser variables
      DIMENSION alpha1(kk),xaser1(kk),xl1(kk),xl3(kk),
     &	 	xlres3(kk),xlhe3(kk)
      COMMON /comlas/ alpha1,elas1,xamda1,xaser1,xecri1,plas1,rabs1,
     &                rocri1,xl1,xl3,nabs1
c end of laser variables
c c4.1 - administrative variables
      LOGICAL nldump,nlemp
      DIMENSION piq(302)
      COMMON /comadm/ maxdim,maxrun,ndump,nrep,nldump,nlemp,piq,tstop
c c4.1 - end administrative variables
c x-ray laser common blocks
      DIMENSION fl(2)
      COMMON /xrl   / ropmul,tpon,tpoff,rmpd,ztrtop,dlamda,fl,fwide,
     &                efbmul,drmul,bbtrap,fbtrap,icxrl,istage,idflg,
     &                itbflg,istflg,idrflg,ipuflg,nlp,nup,llp,lup,nst,
     &                nfl,nlmax,ifrsta,ilosta,ihista,igstat
      DATA iclass,isub/2,3/
      data icall/0/
      save icall
      IF(nlomt2(isub))RETURN
c this version of laser allows a choice of laser pulse shapes.
c singular pulse: store pulse parameters in array xaser1
c as shown below.
c series of 1-5 gaussian pulses: store parameters in xaser1 and use
c xaser1(38) - xaser1(41) to store  maximum power at which each
c pulse peaks (in sequence) ensure piq(37) is non zero (usually 1.0).
c current size of :xaser1: limits number of separate pulses to 5.
c xaser1(42) - xaser1(45) to store the pulse width of each of the npuls-1
c gaussian pulses (notice first pulse uses plenth=xaser1(4))
c xaser1(46) - xaser1((49) to store number of plenth(l=2,5) before the peak
c of the corresponding pulse (notice first pulse uses pmult=xaser1(2)
c element       quantity (where different, 'gaussian' in brackets)
c  1            time to switch on laser
c  2            initial power level
c  3            pulse duration (no. of pulses)
c  4            pulse shape factor (pulse width)
c  5            time to switch off laser
c  6            maximum permissible power level
c  7            maximum energy available (may be omitted)
c  8            spare variable (pulse separation)
c  9            signal indicating laser is switched on
c 10            signal indicating laser is switched off
c 11-15         messages
c ............  notice that these quantities are in mksa units ....
c 1. laser 1
c 1.1 evaluate laser power as a function of time
c should the laser be switched on ?
      IF(tnow.ge.xaser1(1))THEN
c      piq(37)=0.0 isentropic pulse
       IF(abs(piq(37)).lt.underf)THEN
        IF(xaser1(2).le.0.0)GOTO 100
        IF(abs(xaser1(3)).gt.underf)zt1=1.0-tnow/xaser1(3)
        IF(zt1.gt.0.0)THEN
c        has the energy available been delivered ?
         IF(elas1.lt.xaser1(7))THEN
c         is the maximum power exceeded ?
          IF(plas1.ge.xaser1(6))THEN
c          1.2 keep laser power constant
           plas1=xaser1(6)
c          signal constant power level
           xaser1(10)=1.0
          ELSE
c          signals xaser1 switched on
           xaser1(9)=1.0
           xaser1(10)=0.0
           z4=xaser1(4)
c          singular laser pulse
           plas1=xaser1(2)*zt1**z4
          ENDIF
          GOTO 100
         ENDIF
        ENDIF
c       2. modifications
c       should laser be switched off ?
       ELSEIF(tnow.lt.xaser1(5))THEN
        xaser1(9)=1.0
        xaser1(10)=0.0
c       piq(37) < 0 for triangular pulse
c       ad1 gauss = -1 (triangular) -2 (triangular with flat top)
        IF(abs(piq(37)+1.0).lt.underf)THEN
         IF(tnow.gt.xaser1(4))THEN
          plas1=xaser1(6)*(xaser1(5)-tnow)/(xaser1(5)-xaser1(4))
          RETURN
         ELSE
          plas1=xaser1(6)*tnow/xaser1(4)
          RETURN
         ENDIF
        ELSEIF(abs(piq(37)+2.0).lt.underf)THEN
c        ad1
c        triangular pulse with prepulse and flat top can be set as
c prepulse power from 0 at time 0 to xaser1(21) w/m**2 at time xaser1(20
c linear incresase of power to pmax=xaser1(6) for plenth=xaser1(4)
c flat top at xaser1(6) for xaser1(22)
c then linear decrease to zero at xaser1(38)
c        gauss must be = -2.
          plas1=0.
         IF(tnow.lt.xaser1(20))THEN
          plas1=xaser1(21)*tnow/xaser1(20)
          goto 172
         ELSEIF(tnow.lt.xaser1(4)+xaser1(20))THEN
          plas1=xaser1(21)+
     &    (xaser1(6)-xaser1(21))*(tnow-xaser1(20))/xaser1(4)
          goto 172
        ELSEif(tnow.lt.xaser1(4)+xaser1(20)+xaser1(22))THEN
          plas1=xaser1(6)
          goto 172
         ELSEIF(tnow.lt.xaser1(38))THEN
          plas1=xaser1(6)*(xaser1(38)-tnow)
     &          /(xaser1(38)-(xaser1(4)+xaser1(20)+xaser1(22)))
         ENDIF
c add power from the second flat tringular pulse
172      if(nint(xaser1(3)).eq.2) then
c   2nd     triangular pulse with prepulse and flat top can be set as
c prepulse power from 0 at time 0 to xaser1(45) w/m**2 at time xaser1(40
c linear incresase of power to p=xaser1(46) w/m2 for xaser1(41)
c flat top at xaser1(46) for xaser1(42)
c then linear decrease to zero at xaser1(43)
c        anpuls  must be = 2
         IF(tnow.lt.xaser1(40))THEN
          plas1=plas1+xaser1(45)*tnow/xaser1(40)
          return
         ELSEIF(tnow.lt.xaser1(41)+xaser1(40))THEN
          plas1=plas1+xaser1(45)+
     &    (xaser1(46)-xaser1(45))*(tnow-xaser1(40))/xaser1(41)
          return
        ELSEif(tnow.lt.xaser1(41)+xaser1(40)+xaser1(42))THEN
          plas1=plas1+xaser1(46)
          return
         ELSEIF(tnow.lt.xaser1(43))THEN
          plas1=plas1+xaser1(46)*(xaser1(43)-tnow)
     &          /(xaser1(43)-(xaser1(41)+xaser1(40)+xaser1(42)))
         ENDIF
         endif
         return
        ELSE
c        gaussian : set first pulse shape parameters before loop
         zmult=xaser1(2)
         zwth=xaser1(4)
         zt=tnow-zmult*zwth
         ztrtop=zt
         npuls=nint(xaser1(3))
         IF(npuls.lt.1)npuls=1
c        max. power depends on geometry of target
         zpowmx=xaser1(6)
c        zsep = xaser1(8)
c        gaussian pulses
         plas1=0.0
         DO 10 l=1,npuls
c         zetst = zt - (l-1)*zsep
          zetst=zt
          IF(abs(zetst).le.zmult*zwth)THEN
           zetst=(zetst*zetst)/(zwth*zwth)
           IF(zetst.le.170.0)plas1=plas1+zpowmx*exp(-zetst)
          ENDIF
c         set parameters for next pulse if any (npuls?)
          zmult=xaser1(45+l)
          zwth=xaser1(41+l)
          zt=tnow-zmult*zwth
          j=37+l
          zpowmx=xaser1(j)
 10      CONTINUE

c  write out instantaneous intensity and timestep (5/97)
        if (mod(icall,5).eq.0) then
         write (90,*) zt*1.e12,max(plas1,1.e-1),dt2*1.e15
        endif
        icall = icall + 1
         RETURN
        ENDIF
       ENDIF
c      1.3 switch off the laser
       plas1=0.0
       IF(abs(xaser1(9)).gt.underf)THEN
        xaser1(9)=0.0
        xaser1(10)=1.0
c       record laser variables at switch-off time
        xaser1(11)=real(nstep)
        xaser1(12)=tnow
        xaser1(13)=plas1
        xaser1(14)=elas1
        xaser1(15)=rabs1
       ENDIF
      ENDIF

 100  RETURN
      END


      SUBROUTINE lvar(kname,klval)
c u.7 print name and value of logical variable
c c1.1 basic system parameters
      LOGICAL nlend,nlres
      COMMON /combas/ altime,cptime,stime,ndiary,nin,nledge,nonlin,nout,
     &                nprint,npunch,nread,nrec,nresm,nrun,nstep,nlend,
     &                nlres
c end of basic system parameters
      CHARACTER kname*(*)
      LOGICAL klval
      IF(klval)THEN
       WRITE(nout,10100)kname
      ELSE
       WRITE(nout,10200)kname
      ENDIF
      RETURN
10100 FORMAT(4x,a,' =  .true.')
10200 FORMAT(4x,a,' = .false.')
      END


      SUBROUTINE mesage(kmess)
c u.1 print message on output channel
c c1.1 basic system parameters
      LOGICAL nlend,nlres
      COMMON /combas/ altime,cptime,stime,ndiary,nin,nledge,nonlin,nout,
     &                nprint,npunch,nread,nrec,nresm,nrun,nstep,nlend,
     &                nlres
c end of basic system parameters
      CHARACTER kmess*(*)
      WRITE(nout,10100)kmess
      RETURN
10100 FORMAT(4x,a)
      END


      SUBROUTINE modify
c 0.2 modify basic data if required
      RETURN
      END


      SUBROUTINE motion
      PARAMETER(kk=1001,underf=1.E-30)
c 2.2 controls the motion
c c1.9 development and diagnostic parameters
      DIMENSION nadump(20),npdump(20),nvdump(20)
      LOGICAL nlched,nlhead(9),nlomt1(50),nlomt2(50),nlomt3(50),nlrept
      COMMON /comddp/ maxdum,mxdump,nadump,nclass,npdump,npoint,nsub,
     &                nvdump,nlched,nlhead,nlomt1,nlomt2,nlomt3,nlrept
c end of development and diagnostic parameters
c c2.1 hydrodynamics
      DIMENSION dm(kk),p3(kk),r1(kk),r3(kk),r5(kk),rho1(kk),rho3(kk),
     &          rho5(kk),u2(kk),u4(kk),v1(kk),v3(kk),v5(kk),lstat(kk)
      COMMON /comhyd/ dm,p3,pini,r1,r3,r5,rho1,rho3,rho5,rhoini,rhor,
     &                rini,time,u2,u4,uedge,v1,v3,v5
c end of hydrodynamics
c c2.7 physics control
      LOGICAL mlbrms,mlecon,mlfuse,mlicon,mlx,nlabs,nlbrms,nlburn,
     &        nlcri1,nldepo,nlecon,nlfuse,nlicon,nlmove,nlpfe,nlpfi,nlx,
     &        nlhf,nltnl
      COMMON /comcon/ mstep,ncase,ngeom,mlbrms,mlecon,mlfuse,mlicon,mlx,
     &                nlabs,nlbrms,nlburn,nlcri1,nldepo,nlecon,nlfuse,
     &                nlicon,nlmove,nlpfe,nlpfi,nlx,nlhf,nltnl,pmin,
     &                rhomin,temin,timin,umin
c end of physics control
c c3.1 numerical control parameters
      LOGICAL break,nlgoon,nlite
      COMMON /comnc / nceldt,ncondt,nit,nitmax,break,nlgoon,nlite,ak0,
     &                ak1,ak2,ak3,ak4,ak5,ak6,ak7,ak8,ak9,bneum,deltat,
     &                dt2,dt3,dt4,dtemax,dtimax,dtmax,dtmin,dumax,rdt2,
     &                rdt3,rdt4
c end of numerical control parameters
c c4.1 - administrative variables
      LOGICAL nldump,nlemp
      DIMENSION piq(302)
      COMMON /comadm/ maxdim,maxrun,ndump,nrep,nldump,nlemp,piq,tstop
c c4.1 - end administrative variables
      DATA iclass,isub/2,2/
      IF(nlomt2(isub))RETURN
c this routine assumes the following quantities known
c r1, r3, u2, q2, ti1, te1, xl1, yi1, ye1, brems1, omega1
c by calling appropriate routines it advances all quantities one
c level such that at the end of the routine we now know
c r1, r3, r5, u2, u4, q2, q4, ti1, ti3, te1, te3, xl1, xl3,
c brems1, brems3, yi1, yi3, ye1, ye3, eixch2
c at the end we shift all levels back and form derived variables
c where necessary.
c the thermodynamics have been advanced from level 1 to 3 using dt2.
c the hydrodynamics  have been advanced from level 3 to 5 using dt4.
c the time is however advanced from 1 to 3, i.e. time   = time+dt2.
c 1. calculate non-iterative quantities
c initialise for iteration
      CALL stit
      CALL expert(iclass,isub,3)
c evaluate the laser power at level 3
      IF(ncase.eq.1)CALL laser(time)
      CALL expert(iclass,isub,1)
c we precalculate the terms that are not affected by iteration
      CALL gie
c precalculate the time constant for the ion-electron energy
c exchange at level 1
      IF(nlx.and.mlx)CALL xchang(1)
      CALL expert(iclass,isub,2)
c 2. advance ti, te from 1 to 3 & u from 2 to 4
 100  nit=nit+1
c 2.1 boundary conditions at level 3
      IF(ncase.eq.2)CALL boundy(time)
c 2.2 derivatives of energy
c ions:
      CALL statei(1)
c electrons:
      CALL statee(1)
      CALL expert(iclass,isub,4)
c 2.3 particle processes
c the energy absorbed by the electrons (only case 1)
      IF(ncase.eq.1)CALL absorb
c       find the coulomb logarithm
      CALL coulog
      CALL expert(iclass,isub,5)
c in the laser case calculate the hot electron transport
c     if ( ncase.eq.1 ) call hotran
      CALL expert(iclass,isub,6)
c the energy loss by bremsstrahlung
      IF(nlbrms.and.mlbrms)CALL brems
      CALL expert(iclass,isub,7)
c the energy released by thermonuclear reactions
      IF(nlfuse.and.mlfuse)CALL fusion
      CALL expert(iclass,isub,8)
c the exchange of energy between ions and electrons
      IF(nlx.and.mlx)CALL xchang(2)
      CALL expert(iclass,isub,9)
c 2.4 thermal heat conduction
      CALL hcduct
      CALL expert(iclass,isub,10)
c ad1* nlte ionisation?
      IF(abs(piq(58)-2.0).lt.underf)CALL radran
c 2.5 solve energy equation
c first get coefficients
      CALL abcd
      CALL expert(iclass,isub,11)
c find ion and electron temperatures
      CALL findt
c update the ionisation balance
      CALL ionbal
      CALL expert(iclass,isub,12)
c from equations of state we find the corresponding pressures
c ions:
      CALL statei(2)
c electrons:
      CALL statee(2)
      CALL expert(iclass,isub,13)
c then form the hydrodynamic pressure
      CALL formp
      CALL expert(iclass,isub,14)
c 2.6 calculate the next timestep
      CALL timstp
c ad1 fix dt and compare runs
c     dt4=1.156139661*dt2

      dt3=0.5*(dt2+dt4)
      rdt3=1.0/dt3
      rdt4=1.0/dt4
c is the time centering damaged ?
      IF(.not.(break))THEN
       CALL expert(iclass,isub,15)
c 2.7 the hydrodynamic velocities
       IF(nlmove)CALL speed
c if boundary velocity is determined find it here
       tdash=time+0.5*dt4
       IF(ncase.eq.4)CALL boundy(tdash)
       CALL expert(iclass,isub,16)
c the artificial viscous pressure is found from the velocities
       IF(nlmove)CALL neuman
       CALL expert(iclass,isub,17)
c 2.8 move the points to level 5
       IF(nlmove)CALL moveon
c calculate the burn-up of deuterium and tritium
       IF(nlburn)CALL burnup
       CALL expert(iclass,isub,18)
c find volumes, densities etc. at level 5
       CALL volume
       CALL expert(iclass,isub,19)
c 3. examine convergence of iterations
c are iterations required at all?
       IF((.not.nlite))GOTO 200
c is convergence achieved?
       CALL cverge
       CALL expert(iclass,isub,20)
c can we go on?
       IF(nlgoon)GOTO 200
c or is it necessary to break the calculation and go back one level
c using a new value of dt2.
c if we do not break we just perform another iteration
       IF((.not.break))GOTO 100
      ENDIF
c 4. break the calculation if timecentering damaged
c otherwise reset all quantities to level 1 and iterate
      CALL revers
      CALL expert(iclass,isub,21)
c and carry on !!!!!!!!!!!!!!!!!!!????????????????????????
c 5. prepare for the next timestep
c find all energies: kinetic, thermal, energy input and output
 200  CALL energy
c then check that it is reasonable to go on
      CALL exam
c appropriate particle processes are now switched on or off
c if something has gone wrong then end the calculation
      IF(break)CALL endrun
c all levels are rearranged in subprogram shift called by stepon
      RETURN
      END


      SUBROUTINE moveon
      PARAMETER(kk=1001)
c 2.19 advances the lagrangian coordinates
c c1.9 development and diagnostic parameters
      DIMENSION nadump(20),npdump(20),nvdump(20)
      LOGICAL nlched,nlhead(9),nlomt1(50),nlomt2(50),nlomt3(50),nlrept
      COMMON /comddp/ maxdum,mxdump,nadump,nclass,npdump,npoint,nsub,
     &                nvdump,nlched,nlhead,nlomt1,nlomt2,nlomt3,nlrept
c end of development and diagnostic parameters
c c2.1 hydrodynamics
      DIMENSION dm(kk),p3(kk),r1(kk),r3(kk),r5(kk),rho1(kk),rho3(kk),
     &          rho5(kk),u2(kk),u4(kk),v1(kk),v3(kk),v5(kk),lstat(kk)
      COMMON /comhyd/ dm,p3,pini,r1,r3,r5,rho1,rho3,rho5,rhoini,rhor,
     &                rini,time,u2,u4,uedge,v1,v3,v5
c end of hydrodynamics
c c3.1 numerical control parameters
      LOGICAL break,nlgoon,nlite
      COMMON /comnc / nceldt,ncondt,nit,nitmax,break,nlgoon,nlite,ak0,
     &                ak1,ak2,ak3,ak4,ak5,ak6,ak7,ak8,ak9,bneum,deltat,
     &                dt2,dt3,dt4,dtemax,dtimax,dtmax,dtmin,dumax,rdt2,
     &                rdt3,rdt4
c end of numerical control parameters
c c3.2 mesh and numerical methods
      DIMENSION ae(kk),ai(kk),be(kk),bi(kk),ce(kk),ci(kk),de(kk),di(kk),
     &          e(kk),f(kk),ge(kk),gi(kk),q2(kk),q4(kk),teite(kk),
     &          tiite(kk),uite(kk)
      COMMON /comnum/ mesh,nj,njm1,njp1,nl,nlm1,nlp1,ae,ai,be,bi,ce,ci,
     &                de,di,e,f,ge,gi,q2,q4,teite,tiite,uite
c end of mesh and numerical methods
      DATA iclass,isub/2,19/
      IF(nlomt2(isub))RETURN
c 1. move r from level 3 to 5 (eq.74)
      DO 100 j=1,njp1
       r5(j)=r3(j)+u4(j)*dt4
 100  CONTINUE
c the point at r = 0 is at rest
      RETURN
      END


      SUBROUTINE neuman
      PARAMETER(kk=1001)
c 2.18 the artificial viscous pressure
c c1.9 development and diagnostic parameters
      DIMENSION nadump(20),npdump(20),nvdump(20)
      LOGICAL nlched,nlhead(9),nlomt1(50),nlomt2(50),nlomt3(50),nlrept
      COMMON /comddp/ maxdum,mxdump,nadump,nclass,npdump,npoint,nsub,
     &                nvdump,nlched,nlhead,nlomt1,nlomt2,nlomt3,nlrept
c end of development and diagnostic parameters
c c2.1 hydrodynamics
      DIMENSION dm(kk),p3(kk),r1(kk),r3(kk),r5(kk),rho1(kk),rho3(kk),
     &          rho5(kk),u2(kk),u4(kk),v1(kk),v3(kk),v5(kk),lstat(kk)
      COMMON /comhyd/ dm,p3,pini,r1,r3,r5,rho1,rho3,rho5,rhoini,rhor,
     &                rini,time,u2,u4,uedge,v1,v3,v5
c end of hydrodynamics
c c2.2 thermodynamics
      DIMENSION ddroe1(kk),ddroe3(kk),ddroi1(kk),ddroi3(kk),ddte1(kk),
     &          ddte3(kk),ddti1(kk),ddti3(kk),xappae(kk),xappai(kk),
     &          pe1(kk),pe3(kk),pi1(kk),pi3(kk),te1(kk),te3(kk),ti1(kk),
     &          ti3(kk)
      COMMON /comth / ddroe1,ddroe3,ddroi1,ddroi3,ddte1,ddte3,ddti1,
     &                ddti3,gammae,gammai,xappae,xappai,pe1,pe3,pi1,pi3,
     &                te1,te3,teini,ti1,ti3,tiini
c end of thermodynamics
c c2.7 physics control
      LOGICAL mlbrms,mlecon,mlfuse,mlicon,mlx,nlabs,nlbrms,nlburn,
     &        nlcri1,nldepo,nlecon,nlfuse,nlicon,nlmove,nlpfe,nlpfi,nlx,
     &        nlhf,nltnl
      COMMON /comcon/ mstep,ncase,ngeom,mlbrms,mlecon,mlfuse,mlicon,mlx,
     &                nlabs,nlbrms,nlburn,nlcri1,nldepo,nlecon,nlfuse,
     &                nlicon,nlmove,nlpfe,nlpfi,nlx,nlhf,nltnl,pmin,
     &                rhomin,temin,timin,umin
c end of physics control
c c3.1 numerical control parameters
      LOGICAL break,nlgoon,nlite
      COMMON /comnc / nceldt,ncondt,nit,nitmax,break,nlgoon,nlite,ak0,
     &                ak1,ak2,ak3,ak4,ak5,ak6,ak7,ak8,ak9,bneum,deltat,
     &                dt2,dt3,dt4,dtemax,dtimax,dtmax,dtmin,dumax,rdt2,
     &                rdt3,rdt4
c end of numerical control parameters
c c3.2 mesh and numerical methods
      DIMENSION ae(kk),ai(kk),be(kk),bi(kk),ce(kk),ci(kk),de(kk),di(kk),
     &          e(kk),f(kk),ge(kk),gi(kk),q2(kk),q4(kk),teite(kk),
     &          tiite(kk),uite(kk)
      COMMON /comnum/ mesh,nj,njm1,njp1,nl,nlm1,nlp1,ae,ai,be,bi,ce,ci,
     &                de,di,e,f,ge,gi,q2,q4,teite,tiite,uite
c end of mesh and numerical methods
      DATA iclass,isub/2,18/
      IF(nlomt2(isub))RETURN
c 1. calculate grad(u) and div(u)
c the viscous pressure is only applied if du/dr .lt. 0 and divu .lt. 0
      igeom=ngeom-1
      zbsq=bneum*bneum
      DO 100 l=1,nl
       j=l
       q4(l)=0.0
c is the velocity gradient negative?
       zdu4=u4(j+1)-u4(j)
       IF(zdu4.lt.0.0)THEN
        IF(ngeom.eq.1)THEN
        ELSEIF(ngeom.eq.2)THEN
         zdivu4=2.0*(u4(j+1)*(r5(j+1)+r3(j+1))-u4(j)*(r5(j)+r3(j)))
     &          /(r5(j+1)+r5(j)+r3(j+1)+r3(j))
         GOTO 20
        ELSEIF(ngeom.eq.3)THEN
         zdivu4=8.0*(u4(j+1)*(r5(j+1)*r5(j+1)+r3(j+1)*r3(j+1))-u4(j)
     &          *(r5(j)*r5(j)+r3(j)*r3(j)))
     &          /((r5(j+1)+r5(j)+r3(j+1)+r3(j))
     &          *(r5(j+1)+r5(j)+r3(j+1)+r3(j)))
         GOTO 20
        ELSE
         CALL gotoer(' MEDUSA: neuman: goto error, ngeom is',ngeom)
        ENDIF
        zdivu4=u4(j+1)-u4(j)
c 2. is grad(u) and div(u) negative?
 20     IF(zdivu4.lt.0.0)THEN
c only radial compression causes heating (q * grad(u) / div(u) is the
c shock heating)
         zradia=zdu4/zdivu4
c for outward motions include total compression
         IF(zradia.gt.1.0)zradia=1.0
         zrho4=0.5*(rho3(l)+rho5(l))
c the applied viscous pressure
         q4(l)=zbsq*zrho4*zdu4*zdu4*zradia
c acoustic(planar)
c        zsound=sqrt(gammai*p3(l)/zrho4)
c        q4(l)=q4(l)+bneum*zrho4*zsound*abs(zdu4)
        ENDIF
       ENDIF
 100  CONTINUE
      RETURN
      END


c  
c   Output routine modified to give excited state pops
c
      SUBROUTINE output(k)
      PARAMETER(kk=1001,underf=1.E-30)
c 3.1 control the output
c c2.2 thermodynamics
      DIMENSION ddroe1(kk),ddroe3(kk),ddroi1(kk),ddroi3(kk),ddte1(kk),
     &          ddte3(kk),ddti1(kk),ddti3(kk),xappae(kk),xappai(kk),
     &          pe1(kk),pe3(kk),pi1(kk),pi3(kk),te1(kk),te3(kk),ti1(kk),
     &          ti3(kk)
      COMMON /comth / ddroe1,ddroe3,ddroi1,ddroi3,ddte1,ddte3,ddti1,
     &                ddti3,gammae,gammai,xappae,xappai,pe1,pe3,pi1,pi3,
     &                te1,te3,teini,ti1,ti3,tiini
c end of thermodynamics
c c2.1 hydrodynamics
      DIMENSION dm(kk),p3(kk),r1(kk),r3(kk),r5(kk),rho1(kk),rho3(kk),
     &          rho5(kk),u2(kk),u4(kk),v1(kk),v3(kk),v5(kk),lstat(kk)
      COMMON /comhyd/ dm,p3,pini,r1,r3,r5,rho1,rho3,rho5,rhoini,rhor,
     &                rini,time,u2,u4,uedge,v1,v3,v5
c end of hydrodynamics
c c1.1 basic system parameters
      LOGICAL nlend,nlres
      COMMON /combas/ altime,cptime,stime,ndiary,nin,nledge,nonlin,nout,
     &                nprint,npunch,nread,nrec,nresm,nrun,nstep,nlend,
     &                nlres
c end of basic system parameters
c c3.1 numerical control parameters
      LOGICAL break,nlgoon,nlite
      COMMON /comnc / nceldt,ncondt,nit,nitmax,break,nlgoon,nlite,ak0,
     &                ak1,ak2,ak3,ak4,ak5,ak6,ak7,ak8,ak9,bneum,deltat,
     &                dt2,dt3,dt4,dtemax,dtimax,dtmax,dtmin,dumax,rdt2,
     &                rdt3,rdt4
c end of numerical control parameters
c c4.1 - administrative variables
      LOGICAL nldump,nlemp
      DIMENSION piq(302)
      COMMON /comadm/ maxdim,maxrun,ndump,nrep,nldump,nlemp,piq,tstop
c c4.1 - end administrative variables
c c5.1 input-output control variables
      LOGICAL nlfilm,nlhcpy,nlprnt
      DIMENSION buf1(kk),buf2(kk),buf3(kk),buf4(kk),buf5(kk),buf6(kk)
      REAL nprnt,nproc
      COMMON /comout/ nfilm,nhdcpy,np1,np2,np3,nlfilm,nlhcpy,nlprnt,
     &        ztnext,nprnt,nproc,buf1,buf2,buf3,buf4,buf5,buf6,scp,
     &                scr,scrho,scte,scti,sctime,nt1,nt2
c end of input-output control variables
c c3.2 mesh and numerical methods
      DIMENSION ae(kk),ai(kk),be(kk),bi(kk),ce(kk),ci(kk),de(kk),di(kk),
     &          e(kk),f(kk),ge(kk),gi(kk),q2(kk),q4(kk),teite(kk),
     &          tiite(kk),uite(kk)
      COMMON /comnum/ mesh,nj,njm1,njp1,nl,nlm1,nlp1,ae,ai,be,bi,ce,ci,
     &                de,di,e,f,ge,gi,q2,q4,teite,tiite,uite
c end of mesh and numerical methods
      COMMON /atom  / pu(10,kk),p(10,kk),pz(10,kk),zs(10,kk),zst(kk),
     &                zstz(kk),z(kk),a(kk),fn(10),deg(10),qs(10,kk),
     &                zstmi(kk),zstu(kk),engu(10,kk),eng(10,kk),
     &                rdrhe(kk),zstin(kk),nmax
c c2.3 ions and electrons
      DIMENSION brems1(kk),brems3(kk),ddz1(kk),ddz3(kk),degen(kk),
     &          effz(kk),eixch2(kk),fidash(kk),fz1(kk),fz3(kk),fzsq1(kk)
     &          ,fzsq3(kk),xlc(kk),xmieff(kk),xne(kk),xni(kk),omega1(kk)
      COMMON /comie / brems1,brems3,ddz1,ddz3,degen,degmax,degmin,effz,
     &                eixch2,fidash,fz1,fz3,fzsq1,fzsq3,xlc,xmieff,xne,
     &                xni,omega1,pmass
c end of ions and electrons
c c2.7 physics control
      LOGICAL mlbrms,mlecon,mlfuse,mlicon,mlx,nlabs,nlbrms,nlburn,
     &        nlcri1,nldepo,nlecon,nlfuse,nlicon,nlmove,nlpfe,nlpfi,nlx,
     &        nlhf,nltnl
      COMMON /comcon/ mstep,ncase,ngeom,mlbrms,mlecon,mlfuse,mlicon,mlx,
     &                nlabs,nlbrms,nlburn,nlcri1,nldepo,nlecon,nlfuse,
     &                nlicon,nlmove,nlpfe,nlpfi,nlx,nlhf,nltnl,pmin,
     &                rhomin,temin,timin,umin
c end of physics control
c c2.4 laser variables
      DIMENSION alpha1(kk),xaser1(kk),xl1(kk),xl3(kk)
      COMMON /comlas/ alpha1,elas1,xamda1,xaser1,xecri1,plas1,rabs1,
     &                rocri1,xl1,xl3,nabs1
c end of laser variables
c x-ray laser common blocks
      DIMENSION fl(2)
      COMMON /xrl   / ropmul,tpon,tpoff,rmpd,ztrtop,dlamda,fl,fwide,
     &                efbmul,drmul,bbtrap,fbtrap,icxrl,istage,idflg,
     &                itbflg,istflg,idrflg,ipuflg,nlp,nup,llp,lup,nst,
     &                nfl,nlmax,ifrsta,ilosta,ihista,igstat
      DIMENSION rh(kk),te(kk),ti(kk),rb(kk),uu(kk)
      DIMENSION xdta(4000),ydta(6,4000)
      LOGICAL ldump,lemp,lfilm,lhcpy,lprnt

      dimension dnigib(50),teev(kk),tiev(kk),deni(kk),nph(10)
      data iout18/1/

 
      avogadro=6.022045e23
c 1. test if output is required
c preset inner switches
      ldump=.false.
      lemp=.false.
      lfilm=.false.
      lhcpy=.false.
      lprnt=.false.
c frequency testing. couple inner to outer switches
      IF(nstep.eq.(nstep/ndump)*ndump)ldump=nldump
      IF(nstep.eq.(nstep/nfilm)*nfilm)lfilm=nlfilm
      IF(nstep.eq.(nstep/nhdcpy)*nhdcpy)lhcpy=nlhcpy
      IF(nprnt.lt.0.)THEN
       IF(k.ne.2)THEN
        lprnt=nlprnt
        ztlast=time+nprnt*1.0E-12
        ztnext=ztlast-nprnt*1.0E-12
       ENDIF
       IF(time.ge.ztnext)THEN
        lprnt=nlprnt
        lhcpy=nlhcpy
        ztnext=ztnext-nprnt*1.0E-12
       ENDIF
      ELSEIF(nprnt.ge.1..and.nstep.eq.(nstep/int(nprnt))*int(nprnt))THEN
       lprnt=nlprnt
      ENDIF
c
c return here if all switches are off
      IF((.not.lhcpy.and.(.not.lfilm)).and.(.not.lprnt.and.(.not.ldump))
     &   )THEN
c prepare for emergency printing
c we can only get emergency printing if break has occurred
       IF(break.and.nlemp)CALL arrays(0,0)
       RETURN
      ELSE
c 2. selection of output
c construct output and transfer to output buffers
       CALL select
c 3. output
c is printing required?
       IF(lprnt)THEN
        CALL mprint(k)
        CALL dumpstat
cttttttttttttttttttttttttttttttttttttttttttttttttttttttttt
        if(k.ne.3)then
c
c time print
         DO 101 l=1,nl
          rb(l)=r3(l)*100.
          deni(l)=1.e-3*rho3(l)*avogadro/xmieff(l)
          teev(l)=te3(l)/11604.5
          tiev(l)=ti3(l)/11604.5
 101     CONTINUE
         rb(nl+1)=100.0*r3(nl+1)
c which cells to include in printout 
c -  included in input deck + common block /comout/ 28/11/95
c         nt1=1
c         nt2=nl
c
      IF(iout18.EQ.1) THEN
c for DEC alpha
c        OPEN(unit=18,file='week_disk:[public.week.ad1]tpgib.out',
c     &       status='new',recl=200)
c for unix
         OPEN(unit=18,file='med.expops')
c write title once
        WRITE(18,*)' time(s), time wrt peak, Power ',
     & ' then following data for ',nt2-nt1+1,' cells + 1'
        write(18,*)   ' coord(cm)  ',' vel(cm/s)  ',' rho(g/cc)  ',
     & '   zstar    ','   Te(eV)   ','    Ti(eV)  ','  H-(n=1/cc)',
     & '   H-(n=2)  ','   H-(n=3)  ','   H-(n=4)  ',
     & '   H-(n=5)  ','   H-(n=6)  ','   H-(n=7)  ','  H-(n=8)   ',
     & '   H-(n=9)  ','   H-(n=10) '
      endif
      iout18=2
c
c nt1,nt2 first and last cell for which data is printed
      write(18,8444)time,ztrtop,plas1
      do 844 l=nt1,nt2
          DO 54 i=1,10

c       ground and excited confs (see config data)
c for this particular case configuration 1,7 are set to h-like ground,
c 1st, 2nd ... excited states and configs 8,11 for He-like ground,
c 1st, 2nd, etc.. excited states

             icf=i
              CALL loadcg(nph,icf)
              dnigib(i)=frac(nph,l)*deni(l)
   54     CONTINUE
c
      WRITE(18,8444)rb(l),u2(l)*100.,rho3(l)*1.e-3,zst(l),teev(l),
     & tiev(l),(dnigib(i),i=1,10)

c i=1,7 therefore means H-like ground to 6th excited state
c and i=8,11 means He-like ground to 3rd excited state

  844 continue
c finally write boundary velocity
      write(18,8444)rb(nt2+1),u2(nt2+1)*100.
8444  format(16(1pe12.4))
c time print
      endif
ctttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttt
        IF(abs(piq(58)-2.0).lt.underf)THEN
c gain printing section
         DO 10 l=1,nl
          rh(l)=0.001*rho3(l)
          ti(l)=ti3(l)/11605000.0
          te(l)=te3(l)/11605000.0
          rb(l)=100.0*r3(l)
          uu(l)=100.0*u2(l)
 10      CONTINUE
         rb(nl+1)=100.0*r3(nl+1)
         uu(nl+1)=100.0*u2(nl+1)
c       write (6,99006) nstep , time
c x-ray laser common blocks
c99006 format (//15x,':::::: timestep number:',i6,5x,'time : ',1pe14.7)
         IF(icxrl.ne.0)THEN
          IF(istage.eq.5)THEN
c          ni-like gain calc  no implemented yet
c          CALL nigain(rb,rh,ti,te,nst,nfl,nstep,uu,ngeom,idflg,ropmul,
c    &                 ipuflg,llp,lup,rmpd,time,tpon,tpoff)
c          calculate refraction effects
c          CALL nirefr(rh,rb,nst,nfl,fl,nlmax,fwide)
c          calculate axial line intensities
c          CALL niline(rh,ti,rb,time,npt,ny,xdta,ydta,ngeom)
c          tstep=dt2
c          CALL page
c          CALL print(rh,te,rb,nstep,time,tstep,1)
c          CALL niprnt(rh,te,rb,nstep,time,tstep,nprnt)
          ELSEIF(istage.eq.4)THEN
c          ne-like gain calc
           CALL negain(rb,rh,ti,te,nst,nfl,nstep,uu,ngeom,idflg,ropmul,
     &                 ipuflg,llp,lup,rmpd,time,tpon,tpoff)
c          calculate refraction effects
           CALL nerefr(rh,rb,nst,nfl,fl,nlmax,fwide)
c          calculate axial line intensities
           CALL neline(rh,ti,rb,time,npt,ny,xdta,ydta,ngeom)
           tstep=dt2
           CALL page
           CALL print(rh,te,rb,nstep,time,tstep,1)
           CALL neprnt(rh,te,rb,nstep,time,tstep,nprnt)
          ELSEIF(istage.le.3)THEN
c ad1 neon
c          il1=il
c          iu1=iu
c          if (istage.eq.1) then
c           il=1
c           iu=2
c           ed=0.0136*z(nst)**2*(1.0/real(il*il)-1.0/real(iu*iu))
c          elseif (istage.eq.2) then
c           il=2
c           iu=3
c           ed=0.0136*(z(nst)-2.)**2*(1.0/real(il*il)-1.0/real(iu*iu))
c          elseif (istage.eq.3) then
c           il=3
c           iu=4
c           ed=0.0136*(z(nst)-10.)
c    &         **2*(1.0/real(il*il)-1.0/real(iu*iu))
c          endif
c     dlamda in cm, de in kev
c          dlamda=1.24e-07/ed
c ad1 first rec calc
           call gain(istage,dlamda,rh,ti,nst,nfl,istflg)
c          calculate axial line intensities
           call line(nst,nfl,rh,ti,rb,fl,nlmax,time,npt,ny,xdta,ydta,
     &     fwide,istflg,istage,ngeom)
c          calculate refraction effects
           call refr(dlamda,rh,rb,nst,nfl,fl,nlmax)
           tstep=dt2
           call page
           call print(rh,te,rb,nstep,time,tstep,1)
c          gain1=alpha(1,1)
c          if (gain1.le.0.01) gain1=0.01
c          il=il1
c          iu=iu1
c          if (istage.eq.1) then
c           ed=0.0136*z(nst)**2*(1.0/real(il*il)-1.0/real(iu*iu))
c          elseif (istage.eq.2) then
c           ed=0.0136*(z(nst)-2.)**2*(1.0/real(il*il)-1.0/real(iu*iu))
c          elseif (istage.eq.3) then
c           ed=0.0136*(z(nst)-10.)
c    &         **2*(1.0/real(il*il)-1.0/real(iu*iu))
c          endif
c     dlamda in cm, de in kev
c          dlamda=1.24e-07/ed
 
c   ad1 neon       2nd recombination laser gain calculation
c          CALL gain(istage,dlamda,rh,ti,nst,nfl,istflg)
c          calculate axial line intensities
c          CALL line(nst,nfl,rh,ti,rb,fl,nlmax,time,npt,ny,xdta,ydta,
c    &               fwide,istflg,istage,ngeom)
c          calculate refraction effects
c          CALL refr(dlamda,rh,rb,nst,nfl,fl,nlmax)
c          tstep=dt2
c          CALL page
c          CALL print(rh,te,rb,nstep,time,tstep,2)
c ad1 neon
c          gain2=alpha(1,1)
c          if (gain2.le.0.01) gain2=0.01
c          an=6.023e23
c          dnion=(rh(1)*an)/a(1)
c          dne=zst(1)*dnion
c          fff=rh(1)/dne/(1.6e-04)
c          write (22,10100) time,gain1,gain2,dne,te(1)*1000.,xplas1(2)
c    &                      *1.e-4,xl1(1)*fff,eati(1)*fff,zst(1)
c
          ENDIF
         ELSEIF(icxrl.eq.0)THEN
          tstep=dt2
          CALL page
          CALL print(rh,te,rb,nstep,time,tstep,1)
         ENDIF
c end of gain printing section
        ENDIF
       ENDIF
c is graphical output required?
       IF(lhcpy)CALL hdcopy(k)
       RETURN
      ENDIF
10100 FORMAT(1x,10(1pe10.3))
      END


      SUBROUTINE output1(k)
      PARAMETER(kk=1001,underf=1.E-30)
c 3.1 control the output
c c2.2 thermodynamics
      DIMENSION ddroe1(kk),ddroe3(kk),ddroi1(kk),ddroi3(kk),ddte1(kk),
     &          ddte3(kk),ddti1(kk),ddti3(kk),xappae(kk),xappai(kk),
     &          pe1(kk),pe3(kk),pi1(kk),pi3(kk),te1(kk),te3(kk),ti1(kk),
     &          ti3(kk)
      COMMON /comth / ddroe1,ddroe3,ddroi1,ddroi3,ddte1,ddte3,ddti1,
     &                ddti3,gammae,gammai,xappae,xappai,pe1,pe3,pi1,pi3,
     &                te1,te3,teini,ti1,ti3,tiini
c end of thermodynamics
c c2.1 hydrodynamics
      DIMENSION dm(kk),p3(kk),r1(kk),r3(kk),r5(kk),rho1(kk),rho3(kk),
     &          rho5(kk),u2(kk),u4(kk),v1(kk),v3(kk),v5(kk),lstat(kk)
      COMMON /comhyd/ dm,p3,pini,r1,r3,r5,rho1,rho3,rho5,rhoini,rhor,
     &                rini,time,u2,u4,uedge,v1,v3,v5
c end of hydrodynamics
c c1.1 basic system parameters
      LOGICAL nlend,nlres
      COMMON /combas/ altime,cptime,stime,ndiary,nin,nledge,nonlin,nout,
     &                nprint,npunch,nread,nrec,nresm,nrun,nstep,nlend,
     &                nlres
c end of basic system parameters
c c3.1 numerical control parameters
      LOGICAL break,nlgoon,nlite
      COMMON /comnc / nceldt,ncondt,nit,nitmax,break,nlgoon,nlite,ak0,
     &                ak1,ak2,ak3,ak4,ak5,ak6,ak7,ak8,ak9,bneum,deltat,
     &                dt2,dt3,dt4,dtemax,dtimax,dtmax,dtmin,dumax,rdt2,
     &                rdt3,rdt4
c end of numerical control parameters
c c4.1 - administrative variables
      LOGICAL nldump,nlemp
      DIMENSION piq(302)
      COMMON /comadm/ maxdim,maxrun,ndump,nrep,nldump,nlemp,piq,tstop
c c4.1 - end administrative variables
c c5.1 input-output control variables
      LOGICAL nlfilm,nlhcpy,nlprnt
      DIMENSION buf1(kk),buf2(kk),buf3(kk),buf4(kk),buf5(kk),buf6(kk)
      REAL nprnt,nproc
      COMMON /comout/ nfilm,nhdcpy,np1,np2,np3,nlfilm,nlhcpy,nlprnt,
     &        ztnext,nprnt,nproc,buf1,buf2,buf3,buf4,buf5,buf6,scp,
     &                scr,scrho,scte,scti,sctime,nt1,nt2
c end of input-output control variables
c c3.2 mesh and numerical methods
      DIMENSION ae(kk),ai(kk),be(kk),bi(kk),ce(kk),ci(kk),de(kk),di(kk),
     &          e(kk),f(kk),ge(kk),gi(kk),q2(kk),q4(kk),teite(kk),
     &          tiite(kk),uite(kk)
      COMMON /comnum/ mesh,nj,njm1,njp1,nl,nlm1,nlp1,ae,ai,be,bi,ce,ci,
     &                de,di,e,f,ge,gi,q2,q4,teite,tiite,uite
c end of mesh and numerical methods
      COMMON /atom  / pu(10,kk),p(10,kk),pz(10,kk),zs(10,kk),zst(kk),
     &                zstz(kk),z(kk),a(kk),fn(10),deg(10),qs(10,kk),
     &                zstmi(kk),zstu(kk),engu(10,kk),eng(10,kk),
     &                rdrhe(kk),zstin(kk),nmax
c c2.3 ions and electrons
      DIMENSION brems1(kk),brems3(kk),ddz1(kk),ddz3(kk),degen(kk),
     &          effz(kk),eixch2(kk),fidash(kk),fz1(kk),fz3(kk),fzsq1(kk)
     &          ,fzsq3(kk),xlc(kk),xmieff(kk),xne(kk),xni(kk),omega1(kk)
      COMMON /comie / brems1,brems3,ddz1,ddz3,degen,degmax,degmin,effz,
     &                eixch2,fidash,fz1,fz3,fzsq1,fzsq3,xlc,xmieff,xne,
     &                xni,omega1,pmass
c end of ions and electrons
c c2.7 physics control
      LOGICAL mlbrms,mlecon,mlfuse,mlicon,mlx,nlabs,nlbrms,nlburn,
     &        nlcri1,nldepo,nlecon,nlfuse,nlicon,nlmove,nlpfe,nlpfi,nlx,
     &        nlhf,nltnl
      COMMON /comcon/ mstep,ncase,ngeom,mlbrms,mlecon,mlfuse,mlicon,mlx,
     &                nlabs,nlbrms,nlburn,nlcri1,nldepo,nlecon,nlfuse,
     &                nlicon,nlmove,nlpfe,nlpfi,nlx,nlhf,nltnl,pmin,
     &                rhomin,temin,timin,umin
c end of physics control
c c2.4 laser variables
c end of laser variables
c x-ray laser common blocks
      DIMENSION fl(2)
      COMMON /xrl   / ropmul,tpon,tpoff,rmpd,ztrtop,dlamda,fl,fwide,
     &                efbmul,drmul,bbtrap,fbtrap,icxrl,istage,idflg,
     &                itbflg,istflg,idrflg,ipuflg,nlp,nup,llp,lup,nst,
     &                nfl,nlmax,ifrsta,ilosta,ihista,igstat
      DIMENSION rh(kk),te(kk),ti(kk),rb(kk),uu(kk)
      DIMENSION xdta(4000),ydta(6,4000)
      LOGICAL ldump,lemp,lfilm,lhcpy,lprnt


      avogadro=6.022045e23
c 1. test if output is required
c preset inner switches
      ldump=.false.
      lemp=.false.
      lfilm=.false.
      lhcpy=.false.
      lprnt=.false.
c frequency testing. couple inner to outer switches
      IF(nstep.eq.(nstep/ndump)*ndump)ldump=nldump
      IF(nstep.eq.(nstep/nfilm)*nfilm)lfilm=nlfilm
      IF(nstep.eq.(nstep/nhdcpy)*nhdcpy)lhcpy=nlhcpy
      IF(nprnt.lt.0.)THEN
       IF(k.ne.2)THEN
        lprnt=nlprnt
        ztlast=time+nprnt*1.0E-12
        ztnext=ztlast-nprnt*1.0E-12
       ENDIF
       IF(time.ge.ztnext)THEN
        lprnt=nlprnt
        lhcpy=nlhcpy
        ztnext=ztnext-nprnt*1.0E-12
       ENDIF
      ELSEIF(nprnt.ge.1..and.nstep.eq.(nstep/int(nprnt))*int(nprnt))THEN
       lprnt=nlprnt
      ENDIF
c

c return here if all switches are off
      IF((.not.lhcpy.and.(.not.lfilm)).and.(.not.lprnt.and.(.not.ldump))
     &   )THEN
c prepare for emergency printing
c we can only get emergency printing if break has occurred
       IF(break.and.nlemp)CALL arrays(0,0)
       RETURN
      ELSE
c 2. selection of output
c construct output and transfer to output buffers
       CALL select
c 3. output
c is printing required?
       IF(lprnt)THEN
        CALL mprint(k)
        CALL dumpstat
        IF(abs(piq(58)-2.0).lt.underf)THEN
c gain printing section
         DO 10 l=1,nl
          rh(l)=0.001*rho3(l)
          ti(l)=ti3(l)/11605000.0
          te(l)=te3(l)/11605000.0
          rb(l)=100.0*r3(l)
          uu(l)=100.0*u2(l)
 10      CONTINUE
         rb(nl+1)=100.0*r3(nl+1)
         uu(nl+1)=100.0*u2(nl+1)
c       write (6,99006) nstep , time
c x-ray laser common blocks
c99006 format (//15x,':::::: timestep number:',i6,5x,'time : ',1pe14.7)
         IF(icxrl.ne.0)THEN
          IF(istage.eq.5)THEN
c          ni-like gain calc  no implemented yet
c          CALL nigain(rb,rh,ti,te,nst,nfl,nstep,uu,ngeom,idflg,ropmul,
c    &                 ipuflg,llp,lup,rmpd,time,tpon,tpoff)
c          calculate refraction effects
c          CALL nirefr(rh,rb,nst,nfl,fl,nlmax,fwide)
c          calculate axial line intensities
c          CALL niline(rh,ti,rb,time,npt,ny,xdta,ydta,ngeom)
c          tstep=dt2
c          CALL page
c          CALL print(rh,te,rb,nstep,time,tstep,1)
c          CALL niprnt(rh,te,rb,nstep,time,tstep,nprnt)
          ELSEIF(istage.eq.4)THEN
c          ne-like gain calc
           CALL negain(rb,rh,ti,te,nst,nfl,nstep,uu,ngeom,idflg,ropmul,
     &                 ipuflg,llp,lup,rmpd,time,tpon,tpoff)
c          calculate refraction effects
           CALL nerefr(rh,rb,nst,nfl,fl,nlmax,fwide)
c          calculate axial line intensities
           CALL neline(rh,ti,rb,time,npt,ny,xdta,ydta,ngeom)
           tstep=dt2
           CALL page
           CALL print(rh,te,rb,nstep,time,tstep,1)
           CALL neprnt(rh,te,rb,nstep,time,tstep,nprnt)
          ELSEIF(istage.le.3)THEN
c ad1 neon
c          il1=il
c          iu1=iu
c          if (istage.eq.1) then
c           il=1
c           iu=2
c           ed=0.0136*z(nst)**2*(1.0/real(il*il)-1.0/real(iu*iu))
c          elseif (istage.eq.2) then
c           il=2
c           iu=3
c           ed=0.0136*(z(nst)-2.)**2*(1.0/real(il*il)-1.0/real(iu*iu))
c          elseif (istage.eq.3) then
c           il=3
c           iu=4
c           ed=0.0136*(z(nst)-10.)
c    &         **2*(1.0/real(il*il)-1.0/real(iu*iu))
c          endif
c     dlamda in cm, de in kev
c          dlamda=1.24e-07/ed
c ad1 first rec calc
           call gain(istage,dlamda,rh,ti,nst,nfl,istflg)
c          calculate axial line intensities
           call line(nst,nfl,rh,ti,rb,fl,nlmax,time,npt,ny,xdta,ydta,
     &     fwide,istflg,istage,ngeom)
c          calculate refraction effects
           call refr(dlamda,rh,rb,nst,nfl,fl,nlmax)
           tstep=dt2
           call page
           call print(rh,te,rb,nstep,time,tstep,1)
c          gain1=alpha(1,1)
c          if (gain1.le.0.01) gain1=0.01
c          il=il1
c          iu=iu1
c          if (istage.eq.1) then
c           ed=0.0136*z(nst)**2*(1.0/real(il*il)-1.0/real(iu*iu))
c          elseif (istage.eq.2) then
c           ed=0.0136*(z(nst)-2.)**2*(1.0/real(il*il)-1.0/real(iu*iu))
c          elseif (istage.eq.3) then
c           ed=0.0136*(z(nst)-10.)
c    &         **2*(1.0/real(il*il)-1.0/real(iu*iu))
c          endif
c     dlamda in cm, de in kev
c          dlamda=1.24e-07/ed

c   ad1 neon       2nd recombination laser gain calculation
c          CALL gain(istage,dlamda,rh,ti,nst,nfl,istflg)
c          calculate axial line intensities
c          CALL line(nst,nfl,rh,ti,rb,fl,nlmax,time,npt,ny,xdta,ydta,
c    &               fwide,istflg,istage,ngeom)
c          calculate refraction effects
c          CALL refr(dlamda,rh,rb,nst,nfl,fl,nlmax)
c          tstep=dt2
c          CALL page
c          CALL print(rh,te,rb,nstep,time,tstep,2)
c ad1 neon
c          gain2=alpha(1,1)
c          if (gain2.le.0.01) gain2=0.01
c          an=6.023e23
c          dnion=(rh(1)*an)/a(1)
c          dne=zst(1)*dnion
c          fff=rh(1)/dne/(1.6e-04)
c          write (22,10100) time,gain1,gain2,dne,te(1)*1000.,xplas1(2)
c    &                      *1.e-4,xl1(1)*fff,eati(1)*fff,zst(1)
c
          ENDIF
         ELSEIF(icxrl.eq.0)THEN
          tstep=dt2
          CALL page
          CALL print(rh,te,rb,nstep,time,tstep,1)
         ENDIF
c end of gain printing section
        ENDIF
       ENDIF
c is graphical output required?
       IF(lhcpy)CALL hdcopy(k)
       RETURN
      ENDIF
10100 FORMAT(1x,10(1pe10.3))
      END


      SUBROUTINE page
c u.2 fetch new page on output channel
c c1.1 basic system parameters
      LOGICAL nlend,nlres
      COMMON /combas/ altime,cptime,stime,ndiary,nin,nledge,nonlin,nout,
     &                nprint,npunch,nread,nrec,nresm,nrun,nstep,nlend,
     &                nlres
c end of basic system parameters
      WRITE(nout,10100)
      RETURN
10100 FORMAT('1')
      END


      SUBROUTINE polint(xa,ya,n,x,y,dy)
      PARAMETER(nmx=10,underf=1.E-30)
      DIMENSION xa(n),ya(n),c(nmx),d(nmx)
      ns=1
      dif=abs(x-xa(1))
      DO 100 i=1,n
       dift=abs(x-xa(i))
       IF(dift.lt.dif)THEN
        ns=i
        dif=dift
       ENDIF
       c(i)=ya(i)
       d(i)=ya(i)
 100  CONTINUE
      y=ya(ns)
      ns=ns-1
      DO 200 m=1,n-1
       DO 150 i=1,n-m
        ho=xa(i)-x
        hp=xa(i+m)-x
        w=c(i+1)-d(i)
        den=ho-hp
        IF(abs(den).lt.underf)RETURN
        den=w/den
        d(i)=hp*den
        c(i)=ho*den
 150   CONTINUE
       IF(2*ns.lt.n-m)THEN
        dy=c(ns+1)
       ELSE
        dy=d(ns)
        ns=ns-1
       ENDIF
       y=y+dy
 200  CONTINUE
      RETURN
      END


      SUBROUTINE preset
      PARAMETER(kk=1001,underf=1.E-30)
c 1.3 set default values
c c1.1 basic system parameters
      LOGICAL nlend,nlres
      COMMON /combas/ altime,cptime,stime,ndiary,nin,nledge,nonlin,nout,
     &                nprint,npunch,nread,nrec,nresm,nrun,nstep,nlend,
     &                nlres
c end of basic system parameters
c c2.1 hydrodynamics
      DIMENSION dm(kk),p3(kk),r1(kk),r3(kk),r5(kk),rho1(kk),rho3(kk),
     &          rho5(kk),u2(kk),u4(kk),v1(kk),v3(kk),v5(kk),lstat(kk)
      COMMON /comhyd/ dm,p3,pini,r1,r3,r5,rho1,rho3,rho5,rhoini,rhor,
     &                rini,time,u2,u4,uedge,v1,v3,v5
c end of hydrodynamics
c c2.2 thermodynamics
      DIMENSION ddroe1(kk),ddroe3(kk),ddroi1(kk),ddroi3(kk),ddte1(kk),
     &          ddte3(kk),ddti1(kk),ddti3(kk),xappae(kk),xappai(kk),
     &          pe1(kk),pe3(kk),pi1(kk),pi3(kk),te1(kk),te3(kk),ti1(kk),
     &          ti3(kk)
      COMMON /comth / ddroe1,ddroe3,ddroi1,ddroi3,ddte1,ddte3,ddti1,
     &                ddti3,gammae,gammai,xappae,xappai,pe1,pe3,pi1,pi3,
     &                te1,te3,teini,ti1,ti3,tiini
c end of thermodynamics
c c2.3 ions and electrons
      DIMENSION brems1(kk),brems3(kk),ddz1(kk),ddz3(kk),degen(kk),
     &          effz(kk),eixch2(kk),fidash(kk),fz1(kk),fz3(kk),fzsq1(kk)
     &          ,fzsq3(kk),xlc(kk),xmieff(kk),xne(kk),xni(kk),omega1(kk)
      COMMON /comie / brems1,brems3,ddz1,ddz3,degen,degmax,degmin,effz,
     &                eixch2,fidash,fz1,fz3,fzsq1,fzsq3,xlc,xmieff,xne,
     &                xni,omega1,pmass
c end of ions and electrons
c c2.4 laser variables
      DIMENSION alpha1(kk),xaser1(kk),xl1(kk),xl3(kk)
      COMMON /comlas/ alpha1,elas1,xamda1,xaser1,xecri1,plas1,rabs1,
     &                rocri1,xl1,xl3,nabs1
c end of laser variables
c c2.5 thermonuclear reactions
      DIMENSION f1d(kk),f1h(kk),f1he3(kk),f1he4(kk),f1neu(kk),f1ntrl(kk)
     &          ,f1t(kk),f1x(kk),f1x1(kk),f1x2(kk),f1x3(kk),f1x4(kk),
     &          f3d(kk),f3h(kk),f3he3(kk),f3he4(kk),f3neu(kk),f3ntrl(kk)
     &          ,f3t(kk),f3x(kk),f3x1(kk),f3x2(kk),f3x3(kk),f3x4(kk),
     &          r1dd(kk),r1dhe3(kk),r1dt(kk),r3dd(kk),r3dhe3(kk),
     &          r3dt(kk),ye1(kk),ye3(kk),yi1(kk),yi3(kk)
      COMMON /comfus/ deuter,f1d,f1h,f1he3,f1he4,f1neu,f1ntrl,f1t,f1x,
     &                f1x1,f1x2,f1x3,f1x4,f3d,f3h,f3he3,f3he4,f3neu,
     &                f3ntrl,f3t,f3x,f3x1,f3x2,f3x3,f3x4,heliu3,heliu4,
     &                hydrog,xetral,xtrlms,pneut1,pneut3,r1dd,r1dhe3,
     &                r1dt,r3dd,r3dhe3,r3dt,rneut1,rneut3,tinucl,totneu,
     &                tritiu,xmass,xmass1,xmass2,xmass3,xmass4,xtra,
     &                xtra1,xtra2,xtra3,xtra4,xz,xz1,xz2,xz3,xz4,ye1,
     &                ye3,yi1,yi3,yield
c end of thermonuclear reactions
      LOGICAL lphi(kk)
      DIMENSION deltae(50),ebdy(51),emid(50),ha(kk),hb(kk),hc(kk),hd(kk)
     &          ,hdiff(kk),he(kk),heta(kk),hf(kk),hg(kk),hgrp(kk),hj(kk)
     &          ,hlas(kk),hloss(50),hncold(kk),hp(kk),hphi(kk),htaug(kk)
     &          ,htherm(kk),htr(kk),xhot1(kk,50),xhot3(kk,50)
      COMMON /comhot/ igrp,ngroup,lphi,deltae,ebdy,ehot,emid,ha,hb,hc,
     &                hd,hdiff,he,heta,hf,hg,hgrp,hj,hlas,hloss,hncold,
     &                hp,hphi,htaug,htherm,htr,xhot1,xhot3,thot
c c2.6 - energies
      COMMON /comen / eerror,en,einput,eloss,eefuse,eifuse,eindt1,
     &                eindt3,eneutr,pv,usqm,xionen,en0,eions,ecbbs,
     &                ecbfs,edcvs,erbfs,eatis,ezsts,erbbs,erfbs,edfbs
c c2.6 - end energies
c c2.7 physics control
      LOGICAL mlbrms,mlecon,mlfuse,mlicon,mlx,nlabs,nlbrms,nlburn,
     &        nlcri1,nldepo,nlecon,nlfuse,nlicon,nlmove,nlpfe,nlpfi,nlx,
     &        nlhf,nltnl
      COMMON /comcon/ mstep,ncase,ngeom,mlbrms,mlecon,mlfuse,mlicon,mlx,
     &                nlabs,nlbrms,nlburn,nlcri1,nldepo,nlecon,nlfuse,
     &                nlicon,nlmove,nlpfe,nlpfi,nlx,nlhf,nltnl,pmin,
     &                rhomin,temin,timin,umin
c end of physics control
c c3.1 numerical control parameters
      LOGICAL break,nlgoon,nlite
      COMMON /comnc / nceldt,ncondt,nit,nitmax,break,nlgoon,nlite,ak0,
     &                ak1,ak2,ak3,ak4,ak5,ak6,ak7,ak8,ak9,bneum,deltat,
     &                dt2,dt3,dt4,dtemax,dtimax,dtmax,dtmin,dumax,rdt2,
     &                rdt3,rdt4
c end of numerical control parameters
c c3.2 mesh and numerical methods
      DIMENSION ae(kk),ai(kk),be(kk),bi(kk),ce(kk),ci(kk),de(kk),di(kk),
     &          e(kk),f(kk),ge(kk),gi(kk),q2(kk),q4(kk),teite(kk),
     &          tiite(kk),uite(kk)
      COMMON /comnum/ mesh,nj,njm1,njp1,nl,nlm1,nlp1,ae,ai,be,bi,ce,ci,
     &                de,di,e,f,ge,gi,q2,q4,teite,tiite,uite
c end of mesh and numerical methods
c c4.1 - administrative variables
      LOGICAL nldump,nlemp
      DIMENSION piq(302)
      COMMON /comadm/ maxdim,maxrun,ndump,nrep,nldump,nlemp,piq,tstop
c c4.1 - end administrative variables
c c5.1 input-output control variables
      LOGICAL nlfilm,nlhcpy,nlprnt
      DIMENSION buf1(kk),buf2(kk),buf3(kk),buf4(kk),buf5(kk),buf6(kk)
      REAL nprnt,nproc
      COMMON /comout/ nfilm,nhdcpy,np1,np2,np3,nlfilm,nlhcpy,nlprnt,
     &        ztnext,nprnt,nproc,buf1,buf2,buf3,buf4,buf5,buf6,scp,
     &                scr,scrho,scte,scti,sctime,nt1,nt2
c x-ray laser common blocks
      DIMENSION fl(2)
      COMMON /xrl   / ropmul,tpon,tpoff,rmpd,ztrtop,dlamda,fl,fwide,
     &                efbmul,drmul,bbtrap,fbtrap,icxrl,istage,idflg,
     &                itbflg,istflg,idrflg,ipuflg,nlp,nup,llp,lup,nst,
     &                nfl,nlmax,ifrsta,ilosta,ihista,igstat
      COMMON /sornew/ rf1,rf2,rf3,fne2
c end of input-output control variables
c 1. general olympus data
c 1.1 basic system parameters
      nrun=100
c 1.9 diagnostics and program development
c 2. physics
c 2.1 hydrodynamics
      rhoini=124.0
      rhor=0.0
      rini=0.00048
      time=0.0
      uedge=0.0
      IF(abs(piq(11)).lt.underf)zglas=0.0
      IF(abs(piq(61)).lt.underf)zplas=0.0
      fne2=1.00
c parameters for shell targets
      dr1=0.0
      dr2=0.0
      dr3=0.0
      dr4=0.0
c 2.2 thermodynamics
      gammae=5.0/3.0
      gammai=5.0/3.0
      teini=1000.
      tiini=1000.
c 2.3 ions and electrons
      degmax=100.0
      degmin=0.38
      pmass=1.67216*1.0E-27
c 2.4 laser variables
      elas1=0.0
c     laser light of 10 000 angstroem
      xamda1=1.0E-5
      plas1=0.0
      rabs1=0.0
      nabs1=0
c 2.5 thermonuclear reactions
      deuter=0.5
      heliu3=0.0
      heliu4=0.0
      hydrog=0.0
      xetral=0.0
      xtrlms=0.0
      pneut1=0.0
      pneut3=0.0
      rneut1=0.0
      rneut3=0.0
      tinucl=1.0E7
      totneu=0.0
      tritiu=0.5
      xmass=1.0
      xtra=0.0
      xz=1.0
      yield=0.0
      xtra1=0.0
      xtra2=0.0
      xtra3=0.0
      xtra4=0.0
      xmass1=1.0
      xmass2=1.0
      xmass3=1.0
      xmass4=1.0
      xz1=1.0
      xz2=1.0
      xz3=1.0
      xz4=1.0
c 2.6 energies
      eefuse=0.0
      eerror=0.0
      eifuse=0.0
      eindt1=0.0
      eindt3=0.0
      einput=0.0
      eloss=0.0
      eions=0.0
      ecbbs=0.0
      ecbfs=0.0
      edcvs=0.
      erbfs=0.0
      eatis=0.
      erbbs=0.0
      erfbs=0.0
      edfbs=0.0
      en=0.0
      eneutr=0.0
      pv=0.0
      usqm=0.0
c 2.7 physics control
      pmin=1.0E-35
      rhomin=1.0E-35
      temin=1.0
      timin=1.0
      umin=1.0E-6
      mstep=0
      ncase=1
      ngeom=3
      mlbrms=.true.
      mlecon=.true.
      mlfuse=.true.
      mlicon=.true.
      mlx=.true.
      nlabs=.true.
      nlbrms=.true.
      nlburn=.true.
      nlcri1=.true.
      nldepo=.true.
      nlecon=.true.
      nlfuse=.true.
      nlicon=.true.
      nlmove=.true.
      nlpfe=.true.
      nlpfi=.true.
      nlx=.true.
      nlhf=.false.
      nltnl=.false.
c 3. numerical scheme
c 3.1 numerical control parameters
      ak0=1.E08
      ak1=0.25
      ak2=0.25
      ak3=0.25
      ak4=0.25
      ak7=1.0
      ak8=1.00
      ak9=1.00
      bneum=1.0
      deltat=1.0E-18
      dt2=0.0
      dt3=0.0
      dt4=0.0
      dtemax=0.10
      dtimax=0.10
      dtmax=1.E-18
      dtmin=1.E-18
      dumax=0.10
      nceldt=0
      ncondt=0
      nit=0
      nitmax=06
      break=.false.
      nlgoon=.true.
      nlite=.true.
c 3.2 mesh and numerical methods
      mesh=40
c 4. house-keeping
c 4.1 administrative variables
      tstop=1.0E-6
      nrep=0
      ndump=10
      nldump=.false.
      nlemp=.true.
c 5. input-output
c 5.1 input-output control variables
      scp=1.0
      scr=1.0
      scrho=1.0
      scte=1.0
      scti=1.0
      sctime=1.0
      nfilm=1
      nhdcpy=100
      np1=1
      np2=mesh
      np3=mesh/20
      nprnt=100.0
      nlfilm=.true.
      nlhcpy=.false.
      nlprnt=.true.
      nt1=1
      nt2=mesh
c hot electron energy groups
      ngroup=5
      ebdy(1)=0.0
      ebdy(2)=3000.0
      ebdy(3)=6000.0
      ebdy(4)=12000.0
      ebdy(5)=20000.0
      ebdy(6)=3.0E+6
      emid(1)=1500.0
      emid(2)=4500.0
      emid(3)=6000.0
      emid(4)=16000.0
      emid(5)=25000.0
      deltae(1)=1500.0
      deltae(2)=3000.0
      deltae(3)=4500.0
      deltae(4)=7000.0
      deltae(5)=9000.0
c     xrl constants
      icxrl=0
      ifrsta=0
      ilosta=1
      ihista=2
      igstat=0
      rf1=0.99999
      rf2=0.99999
      rf3=0.99999
      istage=0
      idflg=0
      idrflg=0
      drmul=1.
      efbmul=1.0
      bbtrap=0.0
      fbtrap=0.0
      IF(idflg.eq.0)ropmul=0.0
      itbflg=0
      istflg=0
      ipuflg=0
c in cm   dlamda = 81.0e-8
      dlamda=0.
      nlmax=2
      fl(1)=0.001
      fl(nlmax)=0.1
      fwide=0.
      RETURN
      END


c**********************************************************************
      SUBROUTINE depth(rho,l,dr,vel,istage,ropmul)
      PARAMETER(kk=1001,underf=1.E-30)
      COMMON /atom  / pu(10,kk),p(10,kk),pz(10,kk),zs(10,kk),zst(kk),
     &                zstz(kk),z(kk),a(kk),fn(10),deg(10),qs(10,kk),
     &                zstmi(kk),zstu(kk),engu(10,kk),eng(10,kk),
     &                rdrhe(kk),zstin(kk),nmax
      COMMON /atomdt/ foss(10,10),sig(10,10)
      COMMON /opdep / tau(10,10)
      COMMON /npp   / npp(10,10),il,iu
      DIMENSION np(10)
      if(istage.gt.3) return
      pi=3.1415927
      an=6.0232E+23
      const=2.654E-02
      ione=1
      ig=il-1
c calculate optical depth on the resonance lines of the hydrogen-like
c system
c setup population array for the ground state
      DO 100 i=1,nmax
       np(i)=npp(i,ione)
 100  CONTINUE
c if velocity of outside cell is zero (as it is for the first tstep)
c then set the optical depth to zero
      DO 200 j=2,nmax
       DO 150 i=1,j-1
        tau(j,i)=0.0
 150   CONTINUE
 200  CONTINUE
      IF(abs(vel).gt.underf)THEN
       DO 250 j=il,nmax
        IF(istage.eq.1)de=0.0136*(z(l)**2)
     &                    *((1.0/real(ig*ig))-(1.0/real(j*j)))
        IF(istage.eq.2)de=0.0136*(z(l)-2.)
     &                    **2*((1.0/real(ig*ig))-(1.0/real(j*j)))
        IF(istage.eq.3)de=0.0136*(z(l)-10.)
     &                    **2*(1.0/real(ig*ig)-1.0/real(j*j))
        xl=1.24E-07/de
        phi0=xl/vel
        phi0=phi0/sqrt(pi)
        dnion=(rho*an)/a(l)
        tau(j,ig)=const*dnion*frac(np,l)*foss(ig,j)*phi0*dr
c account for fine structure on lyman alpha in h-like systems
      if(istage.eq.1.and.j.eq.2) tau(j,ig)=0.6666666*tau(j,ig)
c insert opacity multiplier
        tau(j,ig)=ropmul*tau(j,ig)
 250   CONTINUE
      ENDIF
      RETURN
      END
c**********************************************************************
      FUNCTION eintm(x)
c this function evaluates the exponential integral * exp(x)
c from abrmowitz and stegun
c function set to zero for x less than or equal to zero
      eintm=0.0
c no evaluation for x less than or equal to zero
      IF(x.le.0.0)THEN
       WRITE(6,10100)x
c evaluation for x less than 1
      ELSEIF(x.gt.1.0)THEN
c evaluation for x greater than 1
       top=0.2677737343+(8.6347608925*x)+(18.0590169730*x*x)
     &     +(8.573328740*x*x*x)+(x*x*x*x)
       bot=3.9584969228+(21.0996530827*x)+(25.6329561486*x*x)
     &     +(9.57332234554*x*x*x)+(x*x*x*x)
       eintm=(top/bot)/x
      ELSE
       eint=-0.57721566-log(x)+(0.99999193*x)-(0.24991055*x*x)
     &      +(0.05519968*x*x*x)-(0.00976004*x*x*x*x)
     &      +(0.00107857*x*x*x*x*x)
       eintm=eint*exp(x)
      ENDIF
      RETURN
10100 FORMAT(6x,'arg of eintm is le 0.0, x = ',e15.8)
      END


      SUBROUTINE ltepop(rho,theta,l,istage,icxrl)
      PARAMETER(kk=1001)
      COMMON /atom  / pu(10,kk),p(10,kk),pz(10,kk),zs(10,kk),zst(kk),
     &                zstz(kk),z(kk),a(kk),fn(10),deg(10),qs(10,kk),
     &                zstmi(kk),zstu(kk),engu(10,kk),eng(10,kk),
     &                rdrhe(kk),zstin(kk),nmax
      DIMENSION px(10)
c this subroutine calculates the lte populations corresponding
c to the present energy levels
      ic=0
      icmax=100
      xtest=5.0E-06
      alp=0.5
c zero px array
      DO 100 i=1,nmax
       px(i)=0.0
 100  CONTINUE
c guess initial estimate of zstar
      CALL zbar(rho,theta,l,zstar)
c iterate to obtain solution
 200  CALL ztp(rho,theta,l,zstar,px)
      zstarn=z(l)
      DO 300 i=1,nmax
       zstarn=zstarn-px(i)
 300  CONTINUE
      IF(zstarn.gt.0.0)THEN
       diff=abs(zstarn-zstar)
       zstar=((1.0-alp)*zstar)+(alp*zstarn)
       ic=ic+1
       IF(ic.gt.icmax)THEN
        WRITE(6,10100)l,diff
       ELSEIF(diff.gt.xtest)THEN
        GOTO 200
       ENDIF
      ELSE
c if iteration produces more than z-.01 electrons reset populations
       ifg=0
       tot=0.0
       DO 350 i=1,nmax
        IF(ifg.ne.0)THEN
         px(i)=0.0
        ELSE
         tot=tot+px(i)
         IF(tot.ge.(z(l)-zstmi(l)))THEN
          px(i)=(z(l)-zstmi(l))-(tot-px(i))
          ifg=ifg+1
         ENDIF
        ENDIF
 350   CONTINUE
      ENDIF
      DO 400 i=1,nmax
       pz(i,l)=px(i)
 400  CONTINUE
      zstz(l)=z(l)
      DO 500 i=1,nmax
       zstz(l)=zstz(l)-pz(i,l)
 500  CONTINUE
      RETURN
10100 FORMAT(6x,'lte pop failure to converge in cell',i4,'diff=',f15.9)
      END


c**********************************************************************
      SUBROUTINE ztp(rho,theta,l,zstar,px)
      PARAMETER(kk=1001)
      COMMON /atom  / pu(10,kk),p(10,kk),pz(10,kk),zs(10,kk),zst(kk),
     &                zstz(kk),z(kk),a(kk),fn(10),deg(10),qs(10,kk),
     &                zstmi(kk),zstu(kk),engu(10,kk),eng(10,kk),
     &                rdrhe(kk),zstin(kk),nmax
      DIMENSION px(10)
c this subroutine calculates the lte populations given zstar and
c the energy levels
      DO 100 i=1,nmax
       px(i)=deg(i)
       eth=eng(i,l)/theta
       IF(eth.le.170.0)THEN
        wa=exp(-eth)*theta*sqrt(theta)
        wb=(317.0*a(l))/(rho*zstar)
        px(i)=px(i)/(1.0+(wa*wb))
       ENDIF
 100  CONTINUE
      RETURN
      END


c**********************************************************************
      SUBROUTINE zbar(rho,theta,l,zstar)
      PARAMETER(kk=1001)
      COMMON /atom  / pu(10,kk),p(10,kk),pz(10,kk),zs(10,kk),zst(kk),
     &                zstz(kk),z(kk),a(kk),fn(10),deg(10),qs(10,kk),
     &                zstmi(kk),zstu(kk),engu(10,kk),eng(10,kk),
     &                rdrhe(kk),zstin(kk),nmax
c this function calculates an approximation to the thomas-fermi
c lte degree of ionisation
      an=6.0232E+23
c form tzero,r,tf
      tev=1000.0*theta
      fthrds=1.3333333333333333
      tzero=tev/(z(l)**fthrds)
      r=rho/(z(l)*a(l))
      tf=tzero/(1.0+tzero)
c setup constants
      a1=0.003323467
      a2=0.97183224
      a3=9.26148E-05
      a4=3.1016524
      b0=-1.762999
      b1=1.4317567
      b2=0.31546338
      c1=-0.36666667
      c2=0.98333333
      alpha=14.3139316
      beta=0.66240046
c calculate a,b and c
      aa=(a1*(tzero**a2))+(a3*(tzero**a4))
      b=-exp(b0+(b1*tf)+(b2*(tf**7)))
      c=c2+(c1*tf)
c calculate q1 and thereby q
      q1=aa*(r**b)
      q=(r**c)+(q1**c)
      cm=1.0/c
      q=q**cm
c calculate x
      x=alpha*(q**beta)
c calculate zstar
      f=x/(1.0+x+(sqrt(1.0+(2.0*x))))
      zstar=f*z(l)
      RETURN
      END


c**********************************************************************
      SUBROUTINE print(rh,te,rb,nstep,time,tstep,kprt)
      PARAMETER(kk=1001)
c this subroutine prints the lte and non-lte populations
c and the gain coefficients for the chosen line
      COMMON /refrac/ refrac(2,kk)
      COMMON /atom  / pu(10,kk),p(10,kk),pz(10,kk),zs(10,kk),zst(kk),
     &                zstz(kk),z(kk),a(kk),fn(10),deg(10),qs(10,kk),
     &                zstmi(kk),zstu(kk),engu(10,kk),eng(10,kk),
     &                rdrhe(kk),zstin(kk),nmax
      COMMON /alph  / rwi(kk,3),frinv(kk,3),alpha(kk,3),ideg(10),
     &                idegm(10)
      COMMON /npp   / npp(10,10),il,iu
      DIMENSION rh(kk),te(kk),rb(kk),ealph(10)
      COMMON /axint / axint(10,6,2),xt(10,6)
c x-ray laser common blocks
      DIMENSION fl(2)
      COMMON /xrl   / ropmul,tpon,tpoff,rmpd,ztrtop,dlamda,fl,fwide,
     &                efbmul,drmul,bbtrap,fbtrap,icxrl,istage,idflg,
     &                itbflg,istflg,idrflg,ipuflg,nlp,nup,llp,lup,nst,
     &                nfl,nlmax,ifrsta,ilosta,ihista,igstat
      INTEGER nabs1
      REAL alpha1(kk),elas1,xamda1,xaser1(kk),xecri1,plas1,rabs1,rocri1,
     &     xl1(kk),xl3(kk)
      COMMON /comlas/ alpha1,elas1,xamda1,xaser1,xecri1,plas1,rabs1,
     &                rocri1,xl1,xl3,nabs1
c end of laser variables
      COMMON /dpi   / ecbb(kk),ecbf(kk),eraddr(kk),erbf(kk),eati(kk),
     &                dcvz(kk),efbx2(kk),efrx2(kk),erbb(kk),erfb(kk),
     &                acfc(10,kk),acfr(10,kk),bcfc(10,kk),bcfr(10,kk),
     &                edfb(kk),ter(kk),dqdtim(kk),nmaxr
      COMMON /radflg/ yamda,xhnu,xihnu,dihnu,xplas1(kk),xrbf(kk),
     &                atrbf(kk),corhf1(kk),corhf(kk),iradf1,iradf2
      COMMON /ratt  / tcbf(kk),tcfb(kk),trfb(kk),trfbst(kk),trbf(kk),
     &                tatir(kk),trdr(kk)
ccccccccccccccc
      CHARACTER*10 shname(12),elect,rtname(12),lasp,ainvb,aati
      COMMON /sname / shname
      DATA elect/'   elec T:'/
      DATA lasp,ainvb,aati/'  P(W/m2):','I-B(W/kg):','ATI(W/kg):'/
      DATA rtname/'  Col.ion/','  Col.rec/','  Rad.ion/','Sp.rd.rec/',
     &     'St.rd.rec/','  Die.rec/','  Tun.ion/','Col bb.in/',
     &     ' Rd bb.in/','          ','   =======','          '/
      DATA iheader/0/

c       population printing
      an=6.0232E+23
      ztrtop=time-xaser1(2)*xaser1(4)
      IF(kprt.ne.2)THEN
       WRITE(6,10900)nstep,time,ztrtop,tstep
       WRITE(6,11000)
c       print non-lte populations
       WRITE(6,10100)
       WRITE(6,10200)(shname(i),i=1,nmax),shname(12)
       WRITE(6,10300)(shname(11),i=1,nmax+1)
       DO 50 l=nst,nfl
        WRITE(6,10400)l,(p(i,l),i=1,nmax),zst(l)
 50    CONTINUE
c       print lte populations
c      CALL page
c      WRITE(6,10600)
c      WRITE(6,10200)(shname(i),i=1,nmax),shname(12)
c      WRITE(6,10300)(shname(11),i=1,nmax+1)
c      DO 100 l=nst,nfl
c       rho=rh(l)
c       theta=te(l)
c       CALL ltepop(rho,theta,l,istage,icxrl)
c       WRITE(6,10400)l,(pz(i,l),i=1,nmax),zstz(l)
c100   CONTINUE
c       print level energies
c      CALL page
c      WRITE(6,10700)
c      WRITE(6,10200)(shname(i),i=1,nmax),elect
c      WRITE(6,10300)(shname(11),i=1,nmax+1)
c      DO 150 l=nst,nfl
c       WRITE(6,10500)l,(1.E3*eng(i,l),i=1,nmax),te(l)*1.E3
c150   CONTINUE
c       total rates per atom
c      CALL page
c      WRITE(6,10800)
c      WRITE(6,10200)(rtname(i),i=1,7),lasp,ainvb,aati
c      WRITE(6,10300)(shname(11),i=1,10)
c      DO 200 l=nst,nfl
c       WRITE(6,10500)l,tcbf(l),tcfb(l),trbf(l),trfb(l),trfbst(l),
c    &                trdr(l),tatir(l),xplas1(l+1),xl1(l),eati(l)
c200   CONTINUE
      ENDIF
c ad1
      IF(nstep.eq.0)RETURN
      IF(icxrl.ne.0)THEN
       IF(istage.gt.3)RETURN
       IF(iheader.eq.0)THEN
        idum=3
        CALL header(11)
        WRITE(11,11100)nfl-nst+1,idum
        iheader=iheader+1
       ENDIF
       WRITE(11,11200)time,ztrtop
c print gain and summaryrun of plasma conditions
       CALL page

       WRITE(6,*)' Gain for lamda(Angs) = ',dlamda*1.E08,' iu=',iu,
     &           ' il=',il

       WRITE(6,11300)
       WRITE(6,11400)
       alphmax=-1.0E30
       lmax=nst
       IF(kprt.eq.2)WRITE(24,12900)time
       DO 250 l=nst,nfl
        dne=(zst(l)*rh(l)*an)/a(l)
        rc=0.5*(rb(l)+rb(l+1))
        IF(alpha(l,1).gt.alphmax)THEN
         lmax=l
         alphmax=alpha(l,1)
        ENDIF
        WRITE(6,11800)l,rc,dne,te(l),alpha(l,1),frinv(l,1),
     &                     (refrac(k,l),k=1,nlmax)
        WRITE(11,11500)rc,alpha(l,1)
        IF(kprt.eq.2)THEN
         fff=rh(l)/dne/(1.6E-04)
         WRITE(24,12900)rc,dne,te(l),alpha(l,1),xplas1(l+1)*1.E-4,xl1(l)
     &                  *fff,eati(l)*fff
        ENDIF
 250   CONTINUE
       dnemax=(zst(lmax)*rh(lmax)*an)/a(lmax)
       rcmax=0.5*(rb(lmax)+rb(lmax+1))
       WRITE(6,12400)
       WRITE(6,11600)
       WRITE(6,11700)
       WRITE(6,11800)lmax,rcmax,dnemax,te(lmax),alpha(lmax,1),
     &  frinv(lmax,1),(refrac(k,lmax),k=1,nlmax)
       WRITE(6,12400)
c ad1 neon
c       print axial intensities and effective temperatures
        iupt = iu + 2
        do 300 k = 1 , nlmax
           write (6,11900) fl(k)
           do 260 m = 1 , 2
              if ( m.eq.1 ) write (6,12000)
              if ( m.eq.2 ) write (6,12100)
              write (6,12200) (axint(j,k,m),j=iu,iupt)
 260       continue
           write (6,12300)
           write (6,12200) (xt(j,k),j=iu,iupt)
           write (6,12400)
 300    continue
c       calculate effective alpha
        if ( fl(1).gt.fl(nlmax) ) then
            ilo = 1
            ish = nlmax
            rlong = fl(1)
            rshort = fl(nlmax)
        else
            ilo = nlmax
            ish = 1
            rlong = fl(nlmax)
            rshort = fl(1)
        endif
        do 350 i1 = iu , iupt
         call effg(axint(i1,ilo,2),axint(i1,ish,2),rlong,rshort,
     &                    ealph(i1))
 350    continue
        write (6,12400)
        write (6,12600)
        write (6,12700) rlong , rshort
        write(6,12800) (ealph(i),i=iu,iupt)
        write(11,12500) (ealph(i),i=iu,iupt)
        write (6,12400)
        write (6,*) '  '
      ENDIF
      IF(kprt.eq.0)RETURN
c print the ground state h-like and he-like number densities
      IF(ifrsta.eq.1)THEN
       CALL page
       CALL gprint(rh,rb,time)
      ENDIF
      RETURN
10100 FORMAT(/7x,'nlte populations:',/)
10200 FORMAT(5x,' cn:',11(a10:))
10300 FORMAT(5x,'====',11(a10:))
10400 FORMAT(5x,i3,':',11F10.5)
10500 FORMAT(5x,i3,':',11(1PE10.3))
10600 FORMAT(/7x,'lte populations:',/)
10700 FORMAT(/7x,'level energies (eV):',/)
10800 FORMAT(/7x,'total rates  (per atom per sec):',/)
10900 FORMAT(15x,'timestep number ',i6,5x,'time = ',1pe14.7,6x,
     &       'time r-to-p = ',1pe14.7,6x,'delta t = ',1pe14.7)
11000 FORMAT(15x,21('-'),6x,21('-'),6x,28('-'),6x,24('-'))
11100 FORMAT(i3,i3)
11200 FORMAT(e15.5,e15.5)
11300 FORMAT(///4x,'    l     radius(cm)     ne(cm**-3)      te(kev)',
     &       '      gain(cm**-1)      inv frac         refraction tol')
11400 FORMAT(4x,'  ===     ==========     ==========     ========',
     &  '==    ============    ==========     ========================='
     &  )
11500 FORMAT(e15.5,e15.5)
11600 FORMAT(5x,'maximum gain point at:')
11700 FORMAT(5x,
     &'  l      radius(cm)     ne(cm**-3)     te(kev)         gain(cm**-
     &1)   inv frac             refraction tol')
11800 FORMAT(5x,i4,7(1pe15.4))
11900 FORMAT(/2x,'fibre length (cm) = ',1pe11.2)
12000 FORMAT(2x,'wide ph s**-1 st**-1   alpha      beta       gamma')
12100 FORMAT(2x,'narr ph s**-1 st**-1   alpha      beta       gamma')
12200 FORMAT(20x,3(1pe11.2))
12300 FORMAT(2x,'saturation factor      alpha      beta       gamma')
12400 FORMAT(///1x,'  ')
12500 FORMAT(e15.5)
12600 FORMAT(1x,'effective gain calculated from narrow line')
12700 FORMAT(1x,'intensities using lengths :',1pe11.5,' / ',1pe11.5,
     &       ' cm')
12800 FORMAT(1x,' alpha : ',1pe15.5,'  beta : ',1pe15.5,'  gamma : ',
     &       1pe15.5)
12900 FORMAT(10(1pe10.3))
      END


c**********************************************************************
      SUBROUTINE gprint(rh,rb,time)
      PARAMETER(kk=1001)
c this subroutine prints the ion number densities of the ground states
c of specified ion species
c ilosta gives the lowest ionisation stage considered
c ihista gives the highest ionisation stage considered
c e.g. ilosta=1 is h-like, ilosta=2 is he-like, etc
c igstat=1 gives ground state population occupation statistics
      COMMON /atom  / pu(10,kk),p(10,kk),pz(10,kk),zs(10,kk),zst(kk),
     &                zstz(kk),z(kk),a(kk),fn(10),deg(10),qs(10,kk),
     &                zstmi(kk),zstu(kk),engu(10,kk),eng(10,kk),
     &                rdrhe(kk),zstin(kk),nmax
c x-ray laser common blocks
      DIMENSION fl(2)
      COMMON /xrl   / ropmul,tpon,tpoff,rmpd,ztrtop,dlamda,fl,fwide,
     &                efbmul,drmul,bbtrap,fbtrap,icxrl,istage,idflg,
     &                itbflg,istflg,idrflg,ipuflg,nlp,nup,llp,lup,nst,
     &                nfl,nlmax,ifrsta,ilosta,ihista,igstat
      DIMENSION rh(kk),rb(kk)
      DIMENSION nph(10),frach(kk,0:50),dnh(kk,0:50)
      CHARACTER*2 aname(50)
      COMMON /atname/ aname
c c2.7 physics control
      LOGICAL mlbrms,mlecon,mlfuse,mlicon,mlx,nlabs,nlbrms,nlburn,
     &        nlcri1,nldepo,nlecon,nlfuse,nlicon,nlmove,nlpfe,nlpfi,nlx,
     &        nlhf,nltnl
      COMMON /comcon/ mstep,ncase,ngeom,mlbrms,mlecon,mlfuse,mlicon,mlx,
     &                nlabs,nlbrms,nlburn,nlcri1,nldepo,nlecon,nlfuse,
     &                nlicon,nlmove,nlpfe,nlpfi,nlx,nlhf,nltnl,pmin,
     &                rhomin,temin,timin,umin
c end of physics control
      dimension sumdnh(0:50)

      INTEGER icount,ifirst
      DATA ifirst/0/
      an=6.0232E+23
c calculate ground state ion number densities
      WRITE(6,10700)
      IF(ilosta.lt.0)ilosta=0
      IF(ihista.gt.50)ihista=50
      DO 100 ionstb=ilosta,ihista,2
       icount=0
       sumdnh(ionstb)=0.
       sumdnh(ionstb+1)=0.
       DO 50 l=nst,nfl
        dnion=(rh(l)*an)/a(l)
        DO 20 ionsta=ionstb,ionstb+1
         IF(ionsta.le.z(l))THEN
          CALL loadcg(nph,ionsta)
          frach(l,ionsta)=frac(nph,l)
          dnh(l,ionsta)=dnion*frach(l,ionsta)
          IF(ionsta.gt.icount)icount=ionsta-ionstb+1
         ELSE
          frach(l,ionsta)=0.0
          dnh(l,ionsta)=0.0
         ENDIF
 20     CONTINUE
c
       IF(ngeom.eq.1)THEN
       ELSEIF(ngeom.eq.2)THEN
        voll=(rb(l+1)*rb(l+1)-rb(l)*rb(l))/2.
        GOTO 56
       ELSEIF(ngeom.eq.3)THEN
        voll=(rb(l+1)*rb(l+1)*rb(l+1)-rb(l)*rb(l)*rb(l))/3.0
        GOTO 56
       ELSE
        CALL gotoer(' MEDUSA: volume: goto error, ngeom is',ngeom)
       ENDIF
       voll=rb(l+1)-rb(l)
 56     continue
        if(abs(z(l)-z(nfl)).lt.0.1) then
         sumdnh(ionstb)=sumdnh(ionstb)+dnh(l,ionstb)*voll
         sumdnh(ionstb+1)=sumdnh(ionstb+1)+dnh(l,ionstb+1)*voll
        endif
c
 50    CONTINUE
       IF(icount.eq.2)THEN
        IF(ionstb.eq.0)THEN
         WRITE(6,10900)
        ELSE
         WRITE(6,10800)aname(ionstb),aname(ionstb),aname(ionstb+1),
     &                 aname(ionstb+1)
        ENDIF
        WRITE(6,11100)
        WRITE(6,11200)
        DO 60 l=nst,nfl
         rc=0.5*(rb(l)+rb(l+1))
         WRITE(6,11300)l,rc,frach(l,ionstb),dnh(l,ionstb),
     &                 frach(l,ionstb+1),dnh(l,ionstb+1)
 60     CONTINUE
        write(6,11302)sumdnh(ionstb),sumdnh(ionstb+1)
       ELSEIF(icount.eq.1)THEN
        IF(ionstb.eq.0)THEN
         WRITE(6,11000)
        ELSE
         WRITE(6,11400)aname(ionstb),aname(ionstb)
        ENDIF
        WRITE(6,11500)
        WRITE(6,11600)
        DO 80 l=nst,nfl
         rc=0.5*(rb(l)+rb(l+1))
         WRITE(6,11700)l,rc,frach(l,ionstb),dnh(l,ionstb)
 80     CONTINUE
         write(6,11702)sumdnh(ionstb)
       ENDIF
c        call page
 100  CONTINUE
      write(17,*)time,(sumdnh(i),i=ilosta,ihista)

c     write(6,*)' i,rc, rdrtot(n=1), idenhe, rdrtot(n=1) * idenhe'
c     write number densities to file
      IF(ifirst.eq.0)THEN
       CALL header(12)
       WRITE(12,*)(nfl-nst+1),ilosta,ihista
       ifirst=1
      ENDIF
      WRITE(12,*)time,ztrtop
      DO 200 l=nst,nfl
       rc=0.5*(rb(l)+rb(l+1))
       WRITE(12,10200)rc,(dnh(l,i),i=ilosta,ihista)
c        write (6,99013) l, rc , rdrhe(l),dnh(l,2), rdrhe(l)*dnh(l,2)
 200  CONTINUE
      IF(igstat.eq.1)THEN
       totmin=1.0
       imin=100
       totmax=0.0
       imax=0
       runtot=0.0
       DO 250 l=nst,nfl
        totfx=0.0
        DO 220 ionsta=0,int(z(l))
         CALL loadcg(nph,ionsta)
         fx=frac(nph,l)
         totfx=totfx+fx
 220    CONTINUE
        IF(totfx.gt.totmax)THEN
         totmax=totfx
         imax=l
        ENDIF
        IF(totfx.lt.totmin)THEN
         totmin=totfx
         imin=l
        ENDIF
        runtot=runtot+totfx
 250   CONTINUE
       WRITE(6,10300)
       WRITE(6,10400)totmax*100.0,imax
       WRITE(6,10500)totmin*100.0,imin
       WRITE(6,10600)(runtot/real(nfl-nst+1))*100.0
      ENDIF
      RETURN
10100 FORMAT(3I5)
10200 FORMAT(1pe12.5,10(1pe10.3))
10300 FORMAT(6x,'ground state occupation statistics:')
10400 FORMAT(6x,'maximum ground state occupation of all stages:',f6.2,
     &       '% at cell:',i2)
10500 FORMAT(6x,'minimum ground state occupation of all stages:',f6.2,
     &       '% at cell:',i2)
10600 FORMAT(6x,'average occupation over all cells and stages:',f6.2,
     &       '%')
10700 FORMAT(/2x,'ground-state number fractions and densities')
10800 FORMAT(/4x,'cell','            c-of-c       ','     ',a2,
     &       '-like        ','     ',a2,'-like        ','     ',a2,
     &       '-like        ','     ',a2,'-like        ')
10900 FORMAT(/4x,'cell            c-of-c       ',
     &       '     stripped            stripped       ',
     &       '      h-like        ','      h-like        ')
11000 FORMAT(/4x,'cell            c-of-c       ',
     &       '     stripped            stripped       ')
11100 FORMAT(4x,' no ','          coordinate     ',
     &       '    no fraction         no.density      ',
     &       '    no fraction         no.density      ')
11200 FORMAT(4x,'====','         ===========     ',
     &       '    ===========         ===========     ',
     &       '    ===========         ===========     ')
11300 FORMAT(4x,i4,5(4x,1pe16.5:))
11302 FORMAT(4x,'total number of ions:',2(4x,1pe16.5:))
11400 FORMAT(/4x,'cell','            c-of-c       ','     ',a2,
     &       '-like        ','     ',a2,'-like        ')
11500 FORMAT(4x,' no ','          coordinate     ',
     &       '    no fraction         no.density      ')
11600 FORMAT(4x,'====','         ===========     ',
     &       '    ===========         ===========     ')
11700 FORMAT(4x,i4,3(4x,1pe16.5:))
11702 FORMAT(4x,'total number of ions:',2(4x,1pe16.5:))
      END


c**********************************************************************
      SUBROUTINE negain(rb,rh,ti,te,nst,nfl,nstep,uu,ngeom,idflg,ropmul,
     &                  ipuflg,llp,lup,rmpd,time,tpon,tpoff)
      PARAMETER(kk=1001,underf=1.E-30)
c given density rh(g/cc), rb in cm , ti,te (kev) for nst to nfl zones,
c calculate gain for most important (delta n(=3) = 0) transitions
c in ne-like germanium.
      PARAMETER(nlvl=27,ntrs=35,n1max=3)
      DIMENSION rates(nlvl,nlvl),fx(nlvl-1),fy(nlvl-1)
      DIMENSION sw(nlvl),icol(nlvl),energ(nlvl),de(ntrs)
      DIMENSION colst(nlvl,nlvl)
      DIMENSION u(nlvl-1,nlvl-1),w(nlvl-1),v(nlvl-1,nlvl-1)
      DIMENSION uu(kk),tiev(kk),teev(kk),fne(10,kk),ff1(kk)
      DIMENSION besc(nlvl,nlvl)
      DIMENSION rb(kk),rh(kk),ti(kk),te(kk),rdc(kk),npne(10),npf(10)
      DIMENSION dne(kk),dni(kk)
      CHARACTER*9 level,slevel,alevel,tlevel,flevel,glevel,clevel
c                                aa quatities
      COMMON /atom  / pu(10,kk),p(10,kk),pz(10,kk),zs(10,kk),zst(kk),
     &                zstz(kk),z(kk),a(kk),fn(10),deg(10),qs(10,kk),
     &                zstmi(kk),zstu(kk),engu(10,kk),eng(10,kk),
     &                rdrhe(kk),zstin(kk),nmax
      COMMON /npp   / npp(10,10),il,iu
      COMMON /rats  / rrc(10,10),cbb(10,10),rbb(10,10),rd(10,10),cbf(10)
     &                ,cfb(10),rfb(10),rfbst(10),rbf(10),atir(10),
     &                atie(10),rdr(10),edrbb
c
      COMMON /ne01  / gain(kk,ntrs),df(kk,ntrs),xl(ntrs),flj(kk,nlvl)
      COMMON /nevoig/ gamt(nlvl,kk),geva(ntrs,kk)
      COMMON /neatom/ raddec(nlvl,nlvl),colstr(nlvl,nlvl),energy(nlvl),
     &                key(nlvl),jtot(nlvl),klow(ntrs),kup(ntrs)
      COMMON /nelevel/ level(nlvl)
      COMMON /siatom/ siraddec(nlvl,nlvl),sicolstr(nlvl,nlvl),
     &                sienergy(nlvl),keysi(nlvl),jtotsi(nlvl),
     &                klowsi(ntrs),kupsi(ntrs)
      COMMON /silevel/ slevel(nlvl)
      COMMON /aratom/ arraddec(nlvl,nlvl),arcolstr(nlvl,nlvl),
     &                arenergy(nlvl),keyar(nlvl),jtotar(nlvl),
     &                klowar(ntrs),kupar(ntrs)
      COMMON /arlevel/ alevel(nlvl)
      COMMON /tiatom/ tiraddec(nlvl,nlvl),ticolstr(nlvl,nlvl),
     &                tienergy(nlvl),keyti(nlvl),jtotti(nlvl),
     &                klowti(ntrs),kupti(ntrs)
      COMMON /tilevel/ tlevel(nlvl)
      COMMON /featom/ feraddec(nlvl,nlvl),fecolstr(nlvl,nlvl),
     &                feenergy(nlvl),keyfe(nlvl),jtotfe(nlvl),
     &                klowfe(ntrs),kupfe(ntrs)
      COMMON /felevel/ flevel(nlvl)
      COMMON /geatom/ geraddec(nlvl,nlvl),gecolstr(nlvl,nlvl),
     &                geenergy(nlvl),keyge(nlvl),jtotge(nlvl),
     &                klowge(ntrs),kupge(ntrs)
      COMMON /gelevel/ glevel(nlvl)
      COMMON /cratom/ crraddec(nlvl,nlvl),crcolstr(nlvl,nlvl),
     &                crenergy(nlvl),keykr(nlvl),jtotkr(nlvl),
     &                klowkr(ntrs),kupkr(ntrs)
      COMMON /crlevel/ clevel(nlvl)
      DATA pi,spl,an/3.1415927,2.9979E+10,6.022045E+23/

c given z of element get correct atomic data
      IF(abs(z(nst)-32.).lt.underf)THEN
       DO 50 i=1,nlvl
        energy(i)=geenergy(i)
        key(i)=keyge(i)
        jtot(i)=jtotge(i)
        level(i)=glevel(i)
        DO 20 j=i+1,nlvl
         raddec(j,i)=geraddec(j,i)
         colstr(i,j)=gecolstr(i,j)
 20     CONTINUE
 50    CONTINUE
       DO 100 i=1,ntrs
        klow(i)=klowge(i)
        kup(i)=kupge(i)
 100   CONTINUE
      ELSEIF(abs(z(nst)-14.).lt.underf)THEN
       DO 150 i=1,nlvl
        energy(i)=sienergy(i)
        key(i)=keysi(i)
        jtot(i)=jtotsi(i)
        level(i)=slevel(i)
        DO 120 j=i+1,nlvl
         raddec(j,i)=siraddec(j,i)
         colstr(i,j)=sicolstr(i,j)
 120    CONTINUE
 150   CONTINUE
       DO 200 i=1,ntrs
        klow(i)=klowsi(i)
        kup(i)=kupsi(i)
 200   CONTINUE
      ELSEIF(abs(z(nst)-18.).lt.underf)THEN
       DO 250 i=1,nlvl
        energy(i)=arenergy(i)
        key(i)=keyar(i)
        jtot(i)=jtotar(i)
        level(i)=alevel(i)
        DO 220 j=i+1,nlvl
         raddec(j,i)=arraddec(j,i)
         colstr(i,j)=arcolstr(i,j)
 220    CONTINUE
 250   CONTINUE
       DO 300 i=1,ntrs
        klow(i)=klowar(i)
        kup(i)=kupar(i)
 300   CONTINUE
      ELSEIF(abs(z(nst)-22.).lt.underf)THEN
       DO 350 i=1,nlvl
        energy(i)=tienergy(i)
        key(i)=keyti(i)
        jtot(i)=jtotti(i)
        level(i)=tlevel(i)
        DO 320 j=i+1,nlvl
         raddec(j,i)=tiraddec(j,i)
         colstr(i,j)=ticolstr(i,j)
 320    CONTINUE
 350   CONTINUE
       DO 400 i=1,ntrs
        klow(i)=klowti(i)
        kup(i)=kupti(i)
 400   CONTINUE
      ELSEIF(abs(z(nst)-26.).lt.underf)THEN
       DO 450 i=1,nlvl
        energy(i)=feenergy(i)
        key(i)=keyfe(i)
        jtot(i)=jtotfe(i)
        level(i)=flevel(i)
        DO 420 j=i+1,nlvl
         raddec(j,i)=feraddec(j,i)
         colstr(i,j)=fecolstr(i,j)
 420    CONTINUE
 450   CONTINUE
       DO 500 i=1,ntrs
        klow(i)=klowfe(i)
        kup(i)=kupfe(i)
 500   CONTINUE
      ELSEIF(abs(z(nst)-36.).lt.underf)THEN
       DO 550 i=1,nlvl
        energy(i)=crenergy(i)
        key(i)=keykr(i)
        jtot(i)=jtotkr(i)
        level(i)=clevel(i)
        DO 520 j=i+1,nlvl
         raddec(j,i)=crraddec(j,i)
         colstr(i,j)=crcolstr(i,j)
 520    CONTINUE
 550   CONTINUE
       DO 600 i=1,ntrs
        klow(i)=klowkr(i)
        kup(i)=kupkr(i)
 600   CONTINUE
      ENDIF
      swtot=0.0
      DO 700 i=1,nlvl
c       convert energies from cm**-1 to kev
       energ(i)=energy(i)*1.24E-07
c       statistical weight swt of levels  =  (2j+1)
       sw(i)=2.0*jtot(i)+1.
       swtot=swtot+sw(i)
       DO 650 l=1,kk
        gamt(i,l)=0.0
        flj(l,i)=0.0
 650   CONTINUE
 700  CONTINUE
c     f like g-state ions
      ionsta=9
      CALL loadcg(npf,ionsta)
      DO 800 j=1,ntrs
c       de in kev
       de(j)=energ(kup(j))-energ(klow(j))
c       lambda (xl) in cm
       xl(j)=1.2400E-07/de(j)
       DO 750 l=nst,nfl
        gain(l,j)=0.0
 750   CONTINUE
 800  CONTINUE
c     write(6,*)npp
      DO 900 l=nst,nfl
       rdc(l)=0.5*(rb(l)+rb(l+1))
       dni(l)=rh(l)*an/a(l)
       dne(l)=zst(l)*dni(l)
       tiev(l)=ti(l)*1000.
       teev(l)=te(l)*1000.
c       ne-like ions
       DO 850 iex=1,nmax
        DO 820 il=1,nmax
         npne(il)=npp(il,iex)
 820    CONTINUE
        fne(iex,l)=frac(npne,l)
 850   CONTINUE
c fluorine-like ions
       ff1(l)=frac(npf,l)
 900  CONTINUE
c     calculate gain in every mesh and every transition
c     average flow velocity
c     vflow=abs(uu(nfl+1)-uu(nst)
c     dr2=rb(nfl+1)-rb(nst)
      gmax=-1.E20
c     write(6,*)'fne',(fne(1,l),l=nst,nfl)
      DO 1100 l=nst,nfl
       IF(fne(1,l).gt.1.E-6)THEN
c       line trapping (delat n=0) :  escape factors besc(j,k)
c       from ground state only
        DO 920 j=2,nlvl
         DO 910 k=1,j-1
          besc(j,k)=1.
          raddec(k,j)=0.0
c        include pumping on the llp -> lup transition for specified
c        period if ipuflg =/ 0
          IF(ipuflg.ne.0)THEN
           IF(k.eq.llp.and.j.eq.lup)THEN
            IF(time.gt.tpon.and.time.lt.tpoff)THEN
c        include pumping on the llp   /k/ -> lup  /j/  transition
c        stimulated emission
             besc(j,k)=(1.0+rmpd)
             raddec(k,j)=raddec(j,k)*sw(j)/sw(k)*rmpd
             IF(abs(time-timeol).gt.underf)WRITE(6,10100)time,k,j,rmpd
             timeol=time
             GOTO 910
            ENDIF
           ENDIF
          ENDIF
c        optical depth factor for resonance lines
c        skip optical depth correction if idflg = 0
          IF(idflg.ne.0)THEN
c        trapping on resonance line into lower levels of gain only
           IF(k.eq.1)THEN
            IF(j.eq.3.or.j.eq.5)THEN
             dele=energ(j)-energ(k)
             xlj=1.24E-07/dele
             veli=1.38E06*sqrt(tiev(l)/a(l))
c        dopw=veli/xlj
             tauj=6.032E-16*raddec(j,k)*sw(j)/sw(k)/dele**2*xlj/sqrt(pi)
c        sobolev escape factor
             vflow=abs(uu(nfl+1)-uu(nst))
c                          average cord length drx
             drx=2.*(rb(nfl+1)-rb(nst))
             IF((vflow).gt.underf)THEN
              drmax=2.*(rb(nfl+1)-rb(nst))
              IF((veli*drx/vflow).gt.drmax)drx=drmax*vflow/veli
              taux=tauj*dni(l)*fne(k,l)*drx/vflow*ropmul
              besc(j,k)=hg(taux)
c        static y  direction
c                              dr = fwide
c                              veli = 1.38e06*sqrt(tiev(l)/a(l))
c                              taus = tauj*dni(l)*fne(k,l)/veli*dr
c        assume an optical thickness in y dir = taus
c                              tauy = taus
c                              besc(j,k) = 0.5*(hg(taux)+hg(tauy))
             ENDIF
            ENDIF
           ENDIF
          ENDIF
 910     CONTINUE
 920    CONTINUE
c       electron collisional excitation rates
        DO 940 j=1,nlvl-1
         DO 930 k=j+1,nlvl
c       collisional excitation ratesin cm**3/s
          exptrm=(energ(k)-energ(j))/te(l)
          IF(exptrm.gt.170.)THEN
           colst(j,k)=0.
           colst(k,j)=0.
          ELSE
           exptrm=exp(-exptrm)
           colst(j,k)=dne(l)*8.0112E-08*colstr(j,k)
     &                *exptrm/(sw(j)*sqrt(teev(l)))
c       de-excitation rate is given by detailed balance
           colst(k,j)=colst(j,k)*sw(j)/(sw(k)*exptrm)
          ENDIF
 930     CONTINUE
 940    CONTINUE
c       elements of rates matrix
        DO 980 j=1,nlvl
         skltj=0.
         skgtj=0.
         DO 950 k=1,j-1
          skltj=skltj+raddec(j,k)*besc(j,k)+colst(j,k)
          rates(j,k)=colst(k,j)+raddec(k,j)
 950     CONTINUE
         DO 960 k=j+1,nlvl
          skgtj=skgtj+colst(j,k)+raddec(j,k)
          rates(j,k)=raddec(k,j)*besc(k,j)+colst(k,j)
 960     CONTINUE
         rates(j,j)=-1.*(skltj+skgtj)
         gamt(j,l)=-rates(j,j)
 980    CONTINUE
c       reduced matrix (2,nlvl) for steady state
c       and solution fx(i)=ni/n1
c       solve the steady state rate equations for
c       fx(i)=ni/n1 in ne-like ions
        DO 1000 j=1,nlvl-1
         fx(j)=0.0
         DO 990 k=1,nlvl-1
          u(j,k)=rates(j+1,k+1)
 990     CONTINUE
         fy(j)=-rates(j+1,1)
 1000   CONTINUE
        mlev=nlvl-1
        nlev=nlvl-1
        mpev=nlvl-1
        npev=nlvl-1
        CALL svdcmp(u,mlev,nlev,mpev,npev,w,v)
        wmax=0.E0
        DO 1020 j=1,nlev
         IF(w(j).gt.wmax)wmax=w(j)
 1020   CONTINUE
        wmin=wmax*1.E-12
        DO 1040 j=1,nlev
         IF(w(j).lt.wmin)w(j)=0.E0
 1040   CONTINUE
        CALL svbksb(u,w,v,mlev,nlev,mpev,npev,fy,fx)
c           write(6,*)l,'fx',fx
c ad1
        DO 1060 j=nlev,1,-1
c              if ( fx(j).lt.0.0 ) then
c                 write (6,*) ' *** error in sv *** fx(' , j , ') = ' ,
c    &                        fx(j)
c                 return
c              endif
         flj(l,j+1)=fx(j)*fne(1,l)
 1060   CONTINUE
        flj(l,1)=fne(1,l)
c***********::::::::::::::::::::
        DO 1080 j=1,ntrs
c         calculate gain for every transition j
         iu=kup(j)
         il=klow(j)
cad1
c           write(6,6611) xl(j)*1e8,flj(l,iu)*raddec(iu,il)
c6611  format(1x,f6.2,1pe9.2)
         dnu=flj(l,iu)*dni(l)
         dnl=flj(l,il)*dni(l)
         IF(abs(dnu).gt.underf)THEN
          df(l,j)=1.-(dnl/dnu)*(sw(iu)/sw(il))
         ELSE
          df(l,j)=-1.E+20
         ENDIF
         wa=xl(j)**2/(8.*pi)
         wb=raddec(iu,il)
         dopw=1.116E+13*de(j)*sqrt(tiev(l)/a(l))
         geva(j,l)=(gamt(il,l)+gamt(iu,l))/(4.0*pi*dopw)
         CALL voigt(geva(j,l),0.0,vfac)
cad1
c     if((il.eq.5.and.iu.eq.15).or.(il.eq.3.and.iu.eq.10).or.(il.eq.5.
c    $ .and.iu.eq.14).or.(il.eq.3.and.iu.eq.9).or.(il.eq.3.and.iu.eq.7)
c    $ ) then
c     write(6,2928)l,xl(j)*1.e8,
c    $ gamt(il,l)+gamt(iu,l),dopw,vfac
c2928 format(1x,'cell:',i2,' ilow:',i2,' iup:',i2,' lamda:',f6.2,
c    & ' gamtot:',1pe12.4,' dopw:',1pe12.4,' vfac:',1pe12.4)
c     endif
c
         wc=1./(sqrt(pi)*dopw)
         IF(abs(dnu).gt.underf)THEN
          wd=df(l,j)*dnu
         ELSE
          wd=-dnl*sw(iu)/sw(il)
         ENDIF
         gain(l,j)=wa*wb*wc*wd*vfac
         IF(gain(l,j).gt.gmax)THEN
          gmax=gain(l,j)
          lmax=l
         ENDIF
 1080   CONTINUE
       ENDIF
 1100 CONTINUE
c output kinetics where max gain
      l=lmax
      IF(gmax.lt.1.)RETURN
      CALL page
c     write(6,*)' l=',l,' step:', nstep,'  kinetics of level population'
c     do 1200 i = 1 , nlvl - 1
c        do 1150 j = i + 1 , nlvl
c           trup = rates(j,i)*flj(l,i)
c           trdo = rates(i,j)*flj(l,j)
c           diff = trdo - trup
c           rat = diff/(trup+trdo)*2.
c1150    continue
c1200 continue
      RETURN
10100 FORMAT(/2x,'pump on time',1pe15.4,' levels',2I4,' rmpd',1pe20.9)
c     write(6,*)' level',i,' lifetime',1./(-rates(i,i)),
c    &           'frac ionic pop ',flj(l,i)
c     write(6,1004)
10200 FORMAT(1x,'  l  u/ rt(l2u) / rt(u2l) / rr(u2l) / fu*rtu2l/',
     &       ' fl*rl2u /  diffe')
c     write(6,1003) i,j,rates(j,i),rates(i,j),raddec(j,i),trdo,trup,diff
10300 FORMAT(1x,i3,i3,7(1pe10.2))
      END


c **********************************************************************
      SUBROUTINE neline(rh,ti,rb,time,npt,ny,xdta,ydta,ngeom)
      PARAMETER(kk=1001)
c units are g/cc , cm, kev
      PARAMETER(nlvl=27,ntrs=35)
      COMMON /atom  / pu(10,kk),p(10,kk),pz(10,kk),zs(10,kk),zst(kk),
     &                zstz(kk),z(kk),a(kk),fn(10),deg(10),qs(10,kk),
     &                zstmi(kk),zstu(kk),engu(10,kk),eng(10,kk),
     &                rdrhe(kk),zstin(kk),nmax
      COMMON /neaxin/ axint(ntrs,6,2),xt(ntrs,6)
      COMMON /grap  / ipt
      COMMON /ne01  / gain(kk,ntrs),df(kk,ntrs),xl(ntrs),flj(kk,nlvl)
      COMMON /nevoig/ gamt(nlvl,kk),geva(ntrs,kk)
      DIMENSION rh(kk),ti(kk),rb(kk)
      DIMENSION sum(ntrs,2)
      DIMENSION tx(501)
      DIMENSION xdta(4000),ydta(6,4000)
c x-ray laser common blocks
      DIMENSION fl(2)
      COMMON /xrl   / ropmul,tpon,tpoff,rmpd,ztrtop,dlamda,fl,fwide,
     &                efbmul,drmul,bbtrap,fbtrap,icxrl,istage,idflg,
     &                itbflg,istflg,idrflg,ipuflg,nlp,nup,llp,lup,nst,
     &                nfl,nlmax,ifrsta,ilosta,ihista,igstat
c this subroutine calculates the axial photon flux in ntrs lines
      pi=3.1415927
      an=6.0232E+23
      const=2.654E-02
      spl=2.9979E10
c ad1  convert fwide to radians (cylindrical geometry only)
      IF(ngeom.gt.1)THEN
c     for cyl geometry fwide in radians
       awide=fwide/rb(nfl+1)
       IF(awide.le.0.0.or.awide.gt.3.141592655)awide=6.283185307
      ELSE
       awide=fwide
       IF(awide.le.0.0)awide=1.
      ENDIF
      glmax=170.

c zero arrays
      DO 100 k=1,nlmax
       DO 50 j=1,ntrs
        DO 20 m=1,2
         axint(j,k,m)=0.0
 20     CONTINUE
        xt(j,k)=0.0
        sum(j,k)=0.0
 50    CONTINUE
 100  CONTINUE
c     loop over the different fibre lengths
      DO 300 k=1,nlmax
c     perform normal and line-narrowed calculations
       DO 150 m=1,2
c     loop over the different cells
        DO 120 l=nst,nfl
c     loop over the different transitions
         DO 110 j=1,ntrs
          IF(gain(l,j).gt.0.0)THEN
           effl=fl(k)
           de=1.24E-07/xl(j)
           veli=1.38E+06*sqrt((1000.0*ti(l))/a(l))
           wa=3.147E+31*(de**2)
           wb=df(l,j)
           wc=gain(l,j)*effl
c           is saturation limit reached
           IF(wc.gt.glmax)THEN
            dxp=exp(glmax)+wc-glmax
c     write(6,*)' ++++++ saturation limit reached in mesh ',l,'+++++++'
           ELSE
            dxp=exp(wc)
           ENDIF
           wd=(wa/wb)*(dxp-1.0)
           IF(ngeom.eq.2)THEN
            sum(j,k)=sum(j,k)+0.5*awide*((rb(l+1)**2)-(rb(l)**2))*wd
           ELSEIF(ngeom.eq.1)THEN
            sum(j,k)=sum(j,k)+awide*(rb(l+1)-rb(l))*wd
           ENDIF
           we=(veli*de)/spl
c          perform integration over the line profile if m=2
c          number of intervals must be exactly divisable by the number
c          of half-widths
           IF(m.eq.1)THEN
c          perform normal line width calculation if m=1
            wd=wd*we*sqrt(pi)
           ELSE
            ip=500
            nhw=10
            IF(mod(nhw,2).ne.0)GOTO 500
            IF(mod(ip,2).ne.0)GOTO 500
            IF(mod(ip,nhw).ne.0)GOTO 500
            itp=ip+1
            DO 102 i=1,itp
             tx(i)=0.0
 102        CONTINUE
            step=we/real(ip/nhw)
            ihp=ip/2
            icp=ihp+1
            tx(icp)=wd
            denp=de
            denm=de
            DO 104 i=1,ihp
             denp=denp+step
             denm=denm-step
c            fill array for integration
             wap=3.147E+31*(denp**2)
             wam=3.147E+31*(denm**2)
             vv=(real(i)*step)/we
             aa=geva(j,l)
             CALL voigt(aa,vv,fexp)
             wc=gain(l,j)*effl*fexp
c            is saturation limit reached
             IF(wc.gt.glmax)THEN
              dxp=exp(glmax)+wc-glmax
c     write(6,*)' ++++++ saturation limit reached in mesh ',l,'+++++++'
             ELSE
              dxp=exp(wc)
             ENDIF
             tx(icp+i)=(wap/wb)*(dxp-1.0)
             tx(icp-i)=(wam/wb)*(dxp-1.0)
 104        CONTINUE
            wd=sint(itp,step,tx)
           ENDIF
c     accumulate area factor
           IF(ngeom.eq.2)THEN
            axint(j,k,m)=axint(j,k,m)
     &                   +0.5*awide*((rb(l+1)**2)-(rb(l)**2))*wd
           ELSEIF(ngeom.eq.1)THEN
            axint(j,k,m)=axint(j,k,m)+awide*(rb(l+1)-rb(l))*wd
           ENDIF
          ENDIF
c         end if gain
 110     CONTINUE
 120    CONTINUE
 150   CONTINUE
c     calculate the effective temperature in each line
       DO 200 n=1,ntrs
        IF(ngeom.eq.2)THEN
         xa=0.5*awide*((rb(nfl+1)**2)-(rb(nst)**2))
        ELSEIF(ngeom.eq.1)THEN
         xa=awide*(rb(nfl+1)-rb(nst))
        ENDIF
        xb=1.24E-07/xl(n)
        xc=3.147E+31*(xb**2)
        xt(n,k)=0.0
        IF(sum(n,k).gt.0.0)THEN
         xd=((xa*xc)/sum(n,k))+1.0
         IF(xd.gt.1.0)THEN
          xt(n,k)=xb/log(xd)
c     assuming an angular spread of 10**-4 steradians calculate
c     saturation factor for each line
          domega=1.0E-04
          xe=domega/(4.0*pi)
          xf=exp(xb/xt(n,k))-1.0
          xt(n,k)=xe/xf
         ENDIF
        ENDIF
 200   CONTINUE
 300  CONTINUE
c fill next element in array for graphics
c draw only the line-narrowed case
      m=2
      iy=0
      ipt=ipt+1
      xdta(ipt)=time
      DO 400 k=1,nlmax
       DO 350 j=1,ntrs
        iy=iy+1
        ydta(iy,ipt)=axint(j,k,m)
 350   CONTINUE
 400  CONTINUE
      npt=ipt
      ny=iy
 500  RETURN
      END


C**********************************************************************
      SUBROUTINE voigt(a,b,h)
C**********************************************************************
C   computes a voigt function  h=h(a,v)
      REAL ak(19),a1(5)
      PARAMETER(underf=1.E-30)
      DATA ak/-1.12470432,-0.15516677,3.28867591,-2.34367915,0.42139162,
     &     -4.48480194,9.39456063,-6.61487486,1.98919585,-0.22041650,
     &     0.554153432,0.278711796,-0.188325687,0.042991293,
     &     -0.003278278,0.979895023,-0.962846325,0.532770573,
     &     -0.122727278/
      DATA sqp/1.772453851/,sq2/1.414213562/
      v=abs(b)
      u=a+v
      v2=v*v
      IF(abs(a).gt.underf)THEN
       IF(a.gt.0.2)THEN
        IF(a.gt.1.4.or.u.gt.3.2)THEN
         a2=a*a
         u=sq2*(a2+v2)
         u2=1.0/(u*u)
C     a gt 1.4  or  a+v gt 3.2
         h=sq2/sqp*a/u*(1.0+u2*(3.0*v2-a2)
     &     +u2*u2*(15.0*v2*v2-30.0*v2*a2+3.0*a2*a2))
         RETURN
        ELSE
         ex=0.0
         IF(v.lt.100.0)ex=exp(-v2)
         k=2
        ENDIF
       ELSEIF(v.ge.5.0)THEN
C     a lt 0.2  and  v ge 5.0
        h=a*(15.0+6.0*v2+4.0*v2*v2)/(4.0*v2*v2*v2*sqp)
        RETURN
       ELSE
        ex=0
        IF(v2.lt.100.0)ex=exp(-v2)
        k=1
       ENDIF
       quo=1.
       IF(v.lt.2.4)THEN
        m=6
        IF(v.lt.1.3)m=1
       ELSE
        quo=1.0/(v2-1.5)
        m=11
       ENDIF
       DO 50 i=1,5
        a1(i)=ak(m)
        m=m+1
 50    CONTINUE
       h1=quo*(a1(1)+v*(a1(2)+v*(a1(3)+v*(a1(4)+v*a1(5)))))
       IF(k.gt.1)THEN
        pqs=2.0/sqp
        h1p=h1+pqs*ex
        h2p=pqs*h1p-2.0*v2*ex
        h3p=(pqs*(1.0-ex*(1.0-2.0*v2))-2.0*v2*h1p)/3.0+pqs*h2p
        h4p=(2.0*v2*v2*ex-pqs*h1p)/3.0+pqs*h3p
        psi=ak(16)+a*(ak(17)+a*(ak(18)+a*ak(19)))
C     0.2 lt a le 1.4  and  a+v le 3.2
        h=psi*(ex+a*(h1p+a*(h2p+a*(h3p+a*h4p))))
        RETURN
       ELSE
C     a le 0.2 and v lt 5
        h=h1*a+ex*(1.0+a*a*(1.0-2.0*v2))
        RETURN
       ENDIF
      ENDIF
C     a eq 0.0
      h=0.0
      IF(v2.lt.100.0)h=exp(-v2)
      RETURN
      END


c**********************************************************************
      SUBROUTINE neeffg(rilong,rishort,rlong,rshort,alpha)
c     this calculates the effective gain as seen from using
c     different target lengths  i.e. the experimental observable.
      alpha=0.
      IF((rilong.lt.1.E-06).or.(rishort.lt.1.E-06))THEN
       WRITE(6,*)' itensities are zero'
       RETURN
      ENDIF
      r=(rilong/rishort)**(2./3.)
      alphmax=170.0/rlong
      alph0=0.0001
      DO 100 i1=1,999
       IF(alph0.gt.alphmax)alph0=alphmax
       IF(alph0*rlong.lt.-170.)THEN
        el=0.0
       ELSE
        el=exp(alph0*rlong)
       ENDIF
       IF(alph0*rshort.lt.-170.)THEN
        es=0.0
       ELSE
        es=exp(alph0*rshort)
       ENDIF
       rdifl3=(rshort-rlong)/3.
       IF(alph0*rdifl3.lt.-170.)THEN
        ed3=0.0
       ELSE
        ed3=exp(alph0*rdifl3)
       ENDIF
       rls3=(rshort/rlong)**(1./3.)
       x1=rls3*ed3*(el-1.)/(es-1.)
       alph1=alph0-(x1-r)/rls3*(es-1.)
     &       /((rshort-rlong)*ed3*(el-1.)/3.+ed3*rlong*el-
     &       rshort*es*x1/rls3)
       diff=abs(alph1/alph0)
       IF((diff.le.1.001).and.(diff.gt.0.999))GOTO 200
       alph0=alph1
 100  CONTINUE
      WRITE(6,*)'not converging in effective alpha'
      WRITE(6,*)'error (%): ',abs(1.0-diff)*100.0
 200  alpha=alph1
      RETURN
      END


      SUBROUTINE nerefr(rh,rb,nst,nfl,fl,nlmax,fwide)
      PARAMETER(kk=1001)
c this subroutine calculates the radius of curvature due to refraction
c from the free electrons in the plasma and calculates the maximum
c transverse dimension which can thereby be tolerated.
      PARAMETER(nlvl=27,ntrs=35)
      DIMENSION rh(kk),rb(kk),rc(kk),fl(2)
      COMMON /atom  / pu(10,kk),p(10,kk),pz(10,kk),zs(10,kk),zst(kk),
     &                zstz(kk),z(kk),a(kk),fn(10),deg(10),qs(10,kk),
     &                zstmi(kk),zstu(kk),engu(10,kk),eng(10,kk),
     &                rdrhe(kk),zstin(kk),nmax
      COMMON /nefrac/ refrac(ntrs,2,kk)
      COMMON /ne01  / gain(kk,ntrs),df(kk,ntrs),xl(ntrs),flj(kk,nlvl)

      DIMENSION rind(kk),drind(kk)
      pi=3.141592654
      spl=2.9979E10
      an=6.0232E+23
      nstp=nst+1
      nflm=nfl-1
      nflm2=nfl-2
      DO 100 l=nst,nfl
       rc(l)=0.5*(rb(l)+rb(l+1))
 100  CONTINUE
      do 600 j=1,ntrs
      DO 200 l=nst,nfl
       dne=(zst(l)*rh(l)*an)/a(l)
       w2=((2.0*pi*spl)/xl(j))**2
       wp2=3.1826E+09*dne
       rind(l)=1.0-(wp2/w2)
       IF(rind(l).lt.0.0)rind(l)=0.0
       rind(l)=sqrt(rind(l))
 200  CONTINUE
      DO 300 l=nst,nflm2
       drind(l+1)=((rind(l+1)-rind(l))/(rc(l+1)-rc(l)))
     &            +((rind(l+2)-rind(l+1))/(rc(l+2)-rc(l+1)))
       drind(l+1)=0.5*drind(l+1)
 300  CONTINUE
      DO 400 l=nstp,nflm
       if(rind(l).gt.0.)radcur=drind(l)/rind(l)
       IF(abs(radcur).gt.0.)radcur=1.0/radcur
       DO 350 k=1,nlmax
        refrac(j,k,l)=radcur
c       refrac(j,k,l)=0.0
c       IF(abs(radcur).gt.0.)refrac(j,k,l)=fl(k)**2/(8.0*radcur)
 350   CONTINUE
 400  CONTINUE
 600  continue

      RETURN
      END


c**********************************************************************
      SUBROUTINE neprnt(rh,te,rb,nstep,time,tstep,nprnt)
      PARAMETER(kk=1001,underf=1.E-30)
c this subroutine prints
c the gain coefficients when positive for the ne-like ge lines
      PARAMETER(nlvl=27,ntrs=35,n1max=3)
      REAL nprnt,nproc
      CHARACTER*4 ss(nlvl)
      CHARACTER*9 level
      DIMENSION rh(kk),te(kk),rb(kk),rc(kk)
      DIMENSION ealph(ntrs),dne(kk),dni(kk),jprint(ntrs)
      COMMON /atom  / pu(10,kk),p(10,kk),pz(10,kk),zs(10,kk),zst(kk),
     &                zstz(kk),z(kk),a(kk),fn(10),deg(10),qs(10,kk),
     &                zstmi(kk),zstu(kk),engu(10,kk),eng(10,kk),
     &                rdrhe(kk),zstin(kk),nmax
      COMMON /ne01  / gain(kk,ntrs),df(kk,ntrs),xl(ntrs),flj(kk,nlvl)
      COMMON /neatom/ raddec(nlvl,nlvl),colstr(nlvl,nlvl),energy(nlvl),
     &                key(nlvl),jtot(nlvl),klow(ntrs),kup(ntrs)
      COMMON /nelevel/ level(nlvl)
      COMMON /neaxin/ axint(ntrs,6,2),xt(ntrs,6)
      COMMON /nevoig/ gamt(nlvl,kk),geva(ntrs,kk)
      COMMON /cvoigt/ bcfl(10,kk),va(kk)
      COMMON /nefrac/ refrac(ntrs,2,kk)
      COMMON /netint/ taxint(ntrs,6),tgain(ntrs)
c x-ray laser common blocks
      DIMENSION fl(2)
      COMMON /xrl   / ropmul,tpon,tpoff,rmpd,ztrtop,dlamda,fl,fwide,
     &                efbmul,drmul,bbtrap,fbtrap,icxrl,istage,idflg,
     &                itbflg,istflg,idrflg,ipuflg,nlp,nup,llp,lup,nst,
     &                nfl,nlmax,ifrsta,ilosta,ihista,igstat
      INTEGER nabs1
      REAL alpha1(kk),elas1,xamda1,xaser1(kk),xecri1,plas1,rabs1,rocri1,
     &     xl1(kk),xl3(kk)
      COMMON /comlas/ alpha1,elas1,xamda1,xaser1,xecri1,plas1,rabs1,
     &                rocri1,xl1,xl3,nabs1
c end of laser variables
      DATA an/6.0232E23/
      DATA iheader/0/
      IF(iheader.eq.0)THEN
       dlamda=0.0
       CALL header(11)
       WRITE(11,10100)nfl-nst+1,ntrs
       WRITE(11,10200)(xl(j)*1.E08,j=1,ntrs)
       iheader=iheader+1
      ENDIF
      WRITE(11,10300)time,ztrtop
c     print gain and summaryrun of plasma conditions
      DO 100 j=1,ntrs
       jprint(j)=0
c      ad1
       DO 50 l=nst,nfl
        IF(gain(l,j).gt.0.0)jprint(j)=1
 50    CONTINUE
 100  CONTINUE
      DO 200 l=nst,nfl
       dni(l)=rh(l)*an/a(l)
       dne(l)=zst(l)*dni(l)
       rc(l)=0.5*(rb(l)+rb(l+1))
 200  CONTINUE
      DO 300 i=1,nlvl
       IF(level(i)(8:8).eq.'1')ss(i)='  0 '
       IF(level(i)(8:8).eq.'2')ss(i)=' 1/2'
       IF(level(i)(8:8).eq.'3')ss(i)='  1 '
       IF(level(i)(8:8).eq.'4')ss(i)=' 3/2'
       IF(level(i)(8:8).eq.'5')ss(i)='  2 '
 300  CONTINUE
      kprint=0
c
      DO 400 j=1,ntrs
       ealph(j)=0.
       tgain(j)=0.
       IF(jprint(j).eq.1)THEN
        kprint=1
        ju=kup(j)
        jl=klow(j)
        jtu=jtot(ju)
        jtl=jtot(jl)
        CALL page
        WRITE(6,10600)nstep,time,ztrtop,tstep
        WRITE(6,10700)
        WRITE(6,10800)j,level(ju),jtu,level(jl),jtl,jtu,jtl,level(ju)
     &                (6:6),level(jl)(6:6),level(ju)(9:9),level(jl)(9:9)
     &                ,ss(ju),ss(jl),ju,jl,xl(j)*1.E+08
        WRITE(6,10900)
        WRITE(6,11000)
        alphmax=-1.0E+06
        lmax=nst
        DO 320 l=nst,nfl
c        ad1
         IF(gain(l,j).gt.0.0)THEN
          IF(gain(l,j).gt.alphmax)THEN
           lmax=l
           alphmax=gain(l,j)
          ENDIF
          WRITE(6,11200)l,rc(l),dne(l),te(l),gain(l,j),geva(j,l),df(l,j)
     &                  ,refrac(j,nlmax,l)
c        ad1
         ENDIF
 320    CONTINUE
c       ad1
c       if(kprint.ne.12345) goto 400
        dnemax=dne(lmax)
        dnimax=dni(lmax)
        rcmax=0.5*(rb(lmax)+rb(lmax+1))
        WRITE(6,11800)
        WRITE(6,11100)
        WRITE(6,10900)
        WRITE(6,11000)
        WRITE(6,11200)lmax,rcmax,dnemax,te(lmax),gain(lmax,j),
     &                geva(j,lmax),df(lmax,j),refrac(j,nlmax,lmax)
c       print axial intensities and effective temperatures
c       photon energy in  joules
        dtint=-(nprnt)*1.E-12
        DO 340 k=1,nlmax
         taxint(j,k)=taxint(j,k)+axint(j,k,2)*dtint
 340    CONTINUE
        demj=1.24E-04/xl(j)*1.6E-19
        DO 360 k=1,nlmax
         WRITE(6,11800)
         WRITE(6,11300)fl(k)
         DO 350 m=1,2
          IF(m.eq.1)WRITE(6,11400)axint(j,k,m),axint(j,k,m)*demj
          IF(m.eq.2)WRITE(6,11500)axint(j,k,m),axint(j,k,m)*demj
 350     CONTINUE
         WRITE(6,11700)xt(j,k)
         IF(k.ne.1)THEN
c        calculate effective alpha
          IF(fl(1).gt.fl(k))THEN
           il=1
           is=k
           rlong=fl(1)
           rshort=fl(k)
          ELSE
           il=k
           is=1
           rlong=fl(k)
           rshort=fl(1)
          ENDIF
          CALL neeffg(axint(j,il,2),axint(j,is,2),rlong,rshort,ealph(j))
          WRITE(6,12200)
          WRITE(6,12300)rlong,rshort,ealph(j)
c         write time integrated x-ray energy flux and gain
          WRITE(6,11900)
          WRITE(6,12000)fl(k),taxint(j,k)*demj
          CALL neeffg(taxint(j,il),taxint(j,is),rlong,rshort,tgain(j))
          WRITE(6,12100)tgain(j)
         ENDIF
 360    CONTINUE
       ENDIF
 400  CONTINUE
c ad1
c     if(kprint.ne.12345) goto 505
      DO 500 l=nst,nfl
       WRITE(11,10400)rc(l),plas1
       WRITE(11,10400)(gain(l,j),j=1,ntrs)
       WRITE(11,10400)(refrac(j,1,l),j=1,ntrs)
 500  CONTINUE
      WRITE(11,10400)(ealph(j),j=1,ntrs)
      WRITE(11,10400)(tgain(j),j=1,ntrs)
      WRITE(6,11800)
      WRITE(6,*)'  '
c print the ground state + 26 excited state densities of ne-like ge ions
      IF(kprint.eq.1)THEN
       CALL page
       WRITE(6,10600)nstep,time,ztrtop,tstep
       WRITE(6,12400)
       DO 550 l=nst,nfl
        IF(abs(flj(l,1)).gt.underf)WRITE(6,10500)l,
     &     (flj(l,k)*dni(l),k=1,nlvl)
 550   CONTINUE
       CALL page
       WRITE(6,10600)nstep,time,ztrtop,tstep
       WRITE(6,12500)
       ig=n1max-1
       DO 600 l=nst,nfl
        IF(abs(flj(l,1)).gt.underf)WRITE(6,10500)l,bcfl(ig,l)*p(ig,l),
     &     (gamt(k,l),k=2,nlvl)
 600   CONTINUE
      ENDIF
      IF(ifrsta.eq.1)THEN
       CALL page
       WRITE(6,10600)nstep,time,ztrtop,tstep
       WRITE(6,10700)
       CALL gprint(rh,rb,time)
      ENDIF
      RETURN
10100 FORMAT(i3,i3)
10200 FORMAT(6(1x,'''',f6.2,''''))
10300 FORMAT(e15.5,e15.5)
10400 FORMAT(6(1x,e12.5))
10500 FORMAT(/6x,i2,6x,9(1pe12.5),/14x,9(1pe12.5),/14x,9(1pe12.5))
10600 FORMAT(15x,'timestep number ',i6,5x,'time = ',1pe14.7,6x,
     &       'time r-to-p = ',1pe14.7,6x,'delta t = ',1pe14.7)
10700 FORMAT(15x,21('-'),6x,21('-'),6x,28('-'),6x,24('-'))
10800 FORMAT(/1x,'Transition',i3,':  ',a9,i1,' --> ',a9,i1,',    J=',i1,
     &       '-->',i1,',    l=',a1,'-->',a1,',   L=',a1,'-->',a1,
     &       ',    S=',a4,'-->',a4,/15x,'  level ',i2,' --> ',i2,11x,
     &       'Wavelength in Angstrom : ',f6.1)
10900 FORMAT(/4x,'    l     radius(cm)     ne(cm**-3)      te(kev)',
     &  '    gain(cm**-1)      Voigt  a     inversion f    ref tol'
     &  )
11000 FORMAT(4x,'  ===     ==========     ==========     ========',
     &  '==  ============    ==========     ============   ============'
     &  )
11100 FORMAT(5x,'maximum gain point at:')
11200 FORMAT(5x,i4,7(1pe15.7))
11300 FORMAT(/2x,'Ge stripe length (cm) = ',1pe11.2)
11400 FORMAT(2x,'wide ph                ',1pe11.2,' s**-1 sr**-1 ','or',
     &       1x,1pe11.2,' W sr**-1')
11500 FORMAT(2x,'narr ph (voigt profile)',1pe11.2,' s**-1 sr**-1 ','or',
     &       1x,1pe11.2,' W sr**-1')
11600 FORMAT(2x,'narr ph(simple formula)',1pe11.2,' s**-1 sr**-1 ','or',
     &       1x,1pe11.2,' W sr**-1')
11700 FORMAT(2x,'saturation factor      ',1pe11.2)
11800 FORMAT(1x,'  ')
11900 FORMAT(2x,'time integrated x-ray energy flux ',
     &       'calculated from narrow lines :')
12000 FORMAT(2x,'Ge stripe length (cm) = ',1pe11.2,3x,' flux = ',
     &       1pe11.2,' J sr**-1')
12100 FORMAT(2x,'time integrated gain   from  ',
     &       'time integrated intensities  = ',1pe11.2)
12200 FORMAT(/2x,'effective gain calculated from narrow line')
12300 FORMAT(2x,'intensity   using lengths :',1pe11.5,' / ',1pe11.5,
     &       ' cm',2x,1pe15.5)
12400 FORMAT(2x,'Number density of Ne-like Ge ions (in ground + first',
     &       ' 26 excited configurations)',//6x,'l=',/6x,'--')
12500 FORMAT(2x,'relaxation rate of ground level (from AA) and  first',
     &       ' 26 excited configurations (CR)',//6x,'l=',/6x,'--')
      END


      SUBROUTINE svdcmp(a,m,n,mp,np,w,v)
c given a matrix a, with logical dimension m by n and physical
c dimensions mp by np, this routine computes its singular value
c decomposition a=u.w.v(t). the matrix u replaces a on output
c the diagonal matrix of singular values w is output as a vector w
c the matrix v (not the transpose v(t)) is output as v. m must be grater
c or equal to n. if it is saller then a must should be filled up
c to square with zero rows.
      PARAMETER(nmx=100,underf=1.E-30)
      DIMENSION a(mp,np),w(np),v(np,np),rv1(nmx)
      g=0.E0
      scale=0.E0
      anorm=0.E0
      DO 100 i=1,n
       l=i+1
       rv1(i)=scale*g
       g=0.E0
       s=0.E0
       scale=0.E0
       IF(i.le.m)THEN
        DO 20 k=i,m
         scale=scale+abs(a(k,i))
 20     CONTINUE
        IF(abs(scale).gt.underf)THEN
         DO 30 k=i,m
          a(k,i)=a(k,i)/scale
          s=s+a(k,i)*a(k,i)
 30      CONTINUE
         f=a(i,i)
         g=-sign(sqrt(s),f)
         h=f*g-s
         a(i,i)=f-g
         IF(i.ne.n)THEN
          DO 35 j=l,n
           s=0.E0
           DO 32 k=i,m
            s=s+a(k,i)*a(k,j)
 32        CONTINUE
           f=s/h
           DO 34 k=i,m
            a(k,j)=a(k,j)+f*a(k,i)
 34        CONTINUE
 35       CONTINUE
         ENDIF
         DO 40 k=i,m
          a(k,i)=scale*a(k,i)
 40      CONTINUE
        ENDIF
       ENDIF
       w(i)=scale*g
       g=0.E0
       s=0.E0
       scale=0.E0
       IF((i.le.m).and.(i.ne.n))THEN
        DO 60 k=l,n
         scale=scale+abs(a(i,k))
 60     CONTINUE
        IF(abs(scale).gt.underf)THEN
         DO 70 k=l,n
          a(i,k)=a(i,k)/scale
          s=s+a(i,k)*a(i,k)
 70      CONTINUE
         f=a(i,l)
         g=-sign(sqrt(s),f)
         h=f*g-s
         a(i,l)=f-g
         DO 80 k=l,n
          rv1(k)=a(i,k)/h
 80      CONTINUE
         IF(i.ne.m)THEN
          DO 85 j=l,m
           s=0.E0
           DO 82 k=l,n
            s=s+a(j,k)*a(i,k)
 82        CONTINUE
           DO 84 k=l,n
            a(j,k)=a(j,k)+s*rv1(k)
 84        CONTINUE
 85       CONTINUE
         ENDIF
         DO 90 k=l,n
          a(i,k)=scale*a(i,k)
 90      CONTINUE
        ENDIF
       ENDIF
       anorm=max(anorm,(abs(w(i))+abs(rv1(i))))
 100  CONTINUE
      DO 200 i=n,1,-1
       IF(i.lt.n)THEN
        IF(abs(g).gt.underf)THEN
         DO 110 j=l,n
          v(j,i)=(a(i,j)/a(i,l))/g
 110     CONTINUE
         DO 130 j=l,n
          s=0.E0
          DO 115 k=l,n
           s=s+a(i,k)*v(k,j)
 115      CONTINUE
          DO 120 k=l,n
           v(k,j)=v(k,j)+s*v(k,i)
 120      CONTINUE
 130     CONTINUE
        ENDIF
        DO 140 j=l,n
         v(i,j)=0.E0
         v(j,i)=0.E0
 140    CONTINUE
       ENDIF
       v(i,i)=1.E0
       g=rv1(i)
       l=i
 200  CONTINUE
      DO 300 i=n,1,-1
       l=i+1
       g=w(i)
       IF(i.lt.n)THEN
        DO 220 j=l,n
         a(i,j)=0.E0
 220    CONTINUE
       ENDIF
       IF(abs(g).gt.underf)THEN
        g=1.E0/g
        IF(i.ne.n)THEN
         DO 240 j=l,n
          s=0.E0
          DO 225 k=l,m
           s=s+a(k,i)*a(k,j)
 225      CONTINUE
          f=(s/a(i,i))*g
          DO 230 k=i,m
           a(k,j)=a(k,j)+f*a(k,i)
 230      CONTINUE
 240     CONTINUE
        ENDIF
        DO 260 j=i,m
         a(j,i)=a(j,i)*g
 260    CONTINUE
       ELSE
        DO 280 j=i,m
         a(j,i)=0.E0
 280    CONTINUE
       ENDIF
       a(i,i)=a(i,i)+1.E0
 300  CONTINUE
      DO 500 k=n,1,-1
       DO 450 its=1,30
        DO 320 l=k,1,-1
         nm=l-1
         IF(abs(rv1(l)).lt.underf)GOTO 380
         IF(abs(w(nm)).lt.underf)GOTO 340
 320    CONTINUE
 340    c=0.E0
        s=1.E0
        DO 360 i=l,k
         f=s*rv1(i)
         IF(abs(f).gt.underf)THEN
          g=w(i)
          h=sqrt(f*f+g*g)
          w(i)=h
          h=1.E0/h
          c=(g*h)
          s=-(f*h)
          DO 345 j=1,m
           y=a(j,nm)
           z=a(j,i)
           a(j,nm)=(y*c)+(z*s)
           a(j,i)=-(y*s)+(z*c)
 345      CONTINUE
         ENDIF
 360    CONTINUE
 380    z=w(k)
        IF(l.eq.k)THEN
         IF(z.lt.0.E0)THEN
          w(k)=-z
          DO 385 j=1,n
           v(j,k)=-v(j,k)
 385      CONTINUE
         ENDIF
         GOTO 500
        ENDIF
        IF(its.eq.30)WRITE(6,*)'svdcmp:no convergence in 30 iterations'
        x=w(l)
        nm=k-1
        y=w(nm)
        g=rv1(nm)
        h=rv1(k)
        f=((y-z)*(y+z)+(g-h)*(g+h))/(2.E0*h*y)
        g=sqrt(f*f+1.E0)
        f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x
        c=1.E0
        s=1.E0
        DO 420 j=l,nm
         i=j+1
         g=rv1(i)
         y=w(i)
         h=s*g
         g=c*g
         z=sqrt(f*f+h*h)
         rv1(j)=z
         c=f/z
         s=h/z
         f=(x*c)+(g*s)
         g=-(x*s)+(g*c)
         h=y*s
         y=y*c
         DO 390 nm=1,n
          x=v(nm,j)
          z=v(nm,i)
          v(nm,j)=(x*c)+(z*s)
          v(nm,i)=-(x*s)+(z*c)
 390     CONTINUE
         z=sqrt(f*f+h*h)
         w(j)=z
         IF(abs(z).gt.underf)THEN
          z=1.E0/z
          c=f*z
          s=h*z
         ENDIF
         f=(c*g)+(s*y)
         x=-(s*g)+(c*y)
         DO 400 nm=1,m
          y=a(nm,j)
          z=a(nm,i)
          a(nm,j)=(y*c)+(z*s)
          a(nm,i)=-(y*s)+(z*c)
 400     CONTINUE
 420    CONTINUE
        rv1(l)=0.E0
        rv1(k)=f
        w(k)=x
 450   CONTINUE
 500  CONTINUE
      RETURN
      END


      SUBROUTINE svbksb(u,w,v,m,n,mp,np,b,x)
c solves a.x=b for a vector x where a is specified by the arrays u,w,v,
c as returned by svdcmp. m and n are the logical dimensions of a, and
c will be equal for square matrices. mp and np are the physical
c dimensions of a. b is the inpu right hand side. x is the output
c solution vector. no input quatities are destroyed, so the routine
c may be called sequentially with different b's
      PARAMETER(nmx=100,underf=1.E-30)
      DIMENSION u(mp,np),w(np),v(np,np),b(mp),x(np),tmp(nmx)
      DO 100 j=1,n
       s=0.E0
       IF(abs(w(j)).gt.underf)THEN
        DO 20 i=1,m
         s=s+u(i,j)*b(i)
 20     CONTINUE
        s=s/w(j)
       ENDIF
       tmp(j)=s
 100  CONTINUE
      DO 200 j=1,n
       s=0.E0
       DO 150 jj=1,n
        s=s+v(j,jj)*tmp(jj)
 150   CONTINUE
       x(j)=s
 200  CONTINUE
      RETURN
      END
c     ad1
c     gexxiii atomic data from a.k.bhatia et al atomic data & nuclear
c                                              tables vol. 32, no. 3.
      BLOCK DATA si
      PARAMETER(nlvl=27,ntrs=35)
      CHARACTER*9 slevel
      COMMON /siatom/ siraddec(nlvl,nlvl),sicolstr(nlvl,nlvl),
     &                sienergy(nlvl),keysi(nlvl),jtotsi(nlvl),
     &                klowsi(ntrs),kupsi(ntrs)
      COMMON /silevel/ slevel(nlvl)

c level configurations and energies (in cm**-1)

      DATA(keysi(i),slevel(i),jtotsi(i),sienergy(i),i=1,nlvl)/1,
     &     '2p6    1S',0,0,2,'2p5 3s 3P',2,838017,3,'2p5 3s 3P',1,
     &     840590,4,'2p5 3s 3P',0,843071,5,'2p5 3s 1P',1,848511,6,
     &     '2p5 3p 3S',1,906252,7,'2p5 3p 3D',3,917929,8,'2p5 3p 3D',2,
     &     918959,9,'2p5 3p 3D',1,920864,10,'2p5 3p 3P',2,924292,11,
     &     '2p5 3p 1P',1,925947,12,'2p5 3p 1D',2,927398,13,'2p5 3p 3P',
     &     0,927806,14,'2p5 3p 3P',1,928405,15,'2p5 3p 1S',0,962950,16,
     &     '2p5 3d 3P',0,1017629,17,'2p5 3d 3P',1,1018236,18,
     &     '2p5 3d 3P',2,1019537,19,'2p5 3d 3F',4,1021384,20,
     &     '2p5 3d 3F',3,1022351,21,'2p5 3d 3F',2,1024240,22,
     &     '2p5 3d 1F',3,1025526,23,'2p5 3d 1D',2,1029340,24,
     &     '2p5 3d 3D',1,1029407,25,'2p5 3d 3D',3,1029875,26,
     &     '2p5 3d 3D',2,1030414,27,'2p5 3d 1P',1,1036915/
c most important transitions
c     data klow/4*1,3,5,6,6,6,7,11,8,9,12,14,10,10,2,2,3,3,2,4,3,4,2,3,
      DATA klowsi/4*1,3,5,6,6,7,11,8,9,12,14,10,10,2,2,3,2,4,3,4,2,3,2,
     &     5,5,3,4,5,5,2,3,4/
c     data kup/27,24,5,3,15,15,18,17,16,19,23,20,21,25,26,22,18,14,12,
c    +         13,12,10,14,10,11,8,9,7,14,12,8,9,11,10,6,6,6/
      DATA kupsi/27,24,5,3,15,15,18,17,19,23,20,21,25,26,22,18,14,12,12,
     &     10,14,10,11,8,9,7,14,12,8,9,11,10,6,6,6/
c radiative decay /raddec/ rates aji (/s) and dimensionless electron
c collision strength /colstr/ at an energy of 110 rydbergs.
      DATA(siraddec(j,1),sicolstr(1,j),j=2,nlvl)/0.000E+00,7.080E-03,
     &     3.700E+09,1.680E-02,0.000E+00,1.420E-03,4.510E+10,6.260E-02,
     &     2.740E-01,1.300E-02,0.000E+00,1.560E-02,4.480E+05,1.200E-02,
     &     0.000E+00,1.560E-02,2.280E+06,1.400E-02,6.540E+00,5.360E-03,
     &     2.260E+06,1.420E-02,0.000E+00,1.860E-03,1.540E+01,5.070E-03,
     &     0.000E+00,2.810E-01,0.000E+00,3.760E-03,6.040E+08,1.160E-02,
     &     0.000E+00,1.790E-02,0.000E+00,1.370E-02,0.000E+00,9.980E-03,
     &     0.000E+00,6.910E-03,0.000E+00,7.200E-03,0.000E+00,4.930E-03,
     &     1.540E+10,1.450E-02,0.000E+00,6.910E-03,0.000E+00,5.190E-03,
     &     3.480E+11,2.610E-01/
      DATA(siraddec(j,2),sicolstr(2,j),j=3,nlvl)/3.970E-01,1.900E-01,
     &     0.000E+00,4.910E-02,2.740E+00,1.630E-01,4.400E+08,1.260E+01,
     &     9.340E+08,4.010E+01,3.740E+08,1.100E+01,8.840E+07,1.450E+00,
     &     5.210E+08,1.210E+01,8.760E+07,1.150E+00,1.980E+08,4.060E+00,
     &     0.000E+00,4.090E-03,2.310E+08,2.750E+00,0.000E+00,6.900E-03,
     &     0.000E+00,1.100E-01,4.350E+04,2.810E-01,2.760E+04,2.920E-01,
     &     5.330E+04,9.490E-01,2.170E+04,3.180E-01,1.130E+04,1.200E-01,
     &     2.300E+04,3.060E-01,5.920E+03,5.530E-02,9.190E+03,5.410E-02,
     &     1.260E+04,1.400E-01,1.110E+04,8.940E-02,1.350E+02,1.190E-02/
      DATA(siraddec(j,3),sicolstr(3,j),j=4,nlvl)/6.880E-01,6.280E-02,
     &     6.800E-01,9.530E-02,1.560E+08,5.130E+00,0.000E+00,2.690E-02,
     &     5.340E+08,1.750E+01,5.620E+08,1.020E+01,2.350E+08,6.080E+00,
     &     6.910E+07,1.020E+00,2.040E+08,4.640E+00,1.100E+09,4.960E+00,
     &     9.150E+07,1.220E+00,4.080E+08,2.740E-01,0.000E+00,2.620E-03,
     &     7.200E+03,5.520E-02,1.850E+04,2.080E-01,0.000E+00,2.350E-02,
     &     2.820E+04,4.010E-01,3.100E+04,3.000E-01,5.540E+03,8.820E-02,
     &     1.190E+02,1.610E-02,5.030E+04,2.490E-01,2.080E+04,2.490E-01,
     &     2.490E+03,3.360E-02,5.740E+02,1.190E-02/
      DATA(siraddec(j,4),sicolstr(4,j),j=5,nlvl)/3.790E-01,3.480E-02,
     &     4.000E+07,1.490E+00,0.000E+00,1.170E-03,0.000E+00,5.440E-03,
     &     2.680E+08,5.410E+00,0.000E+00,1.290E-03,3.470E+08,5.560E+00,
     &     0.000E+00,1.440E-02,0.000E+00,2.310E-03,3.150E+08,4.570E+00,
     &     0.000E+00,1.540E-03,0.000E+00,1.850E-04,8.350E-03,7.550E-04,
     &     4.880E+03,5.820E-02,0.000E+00,1.480E-04,0.000E+00,1.840E-03,
     &     9.840E+03,1.480E-04,0.000E+00,2.610E-03,2.280E+04,2.030E-01,
     &     1.410E-03,2.980E-03,0.000E+00,1.520E-02,1.730E+04,1.540E-01,
     &     3.920E-02,4.940E-03/
      DATA(siraddec(j,5),sicolstr(5,j),j=6,nlvl)/3.600E+06,2.020E-01,
     &     0.000E+00,1.700E-02,5.190E+06,2.740E-01,1.960E+06,6.690E-02,
     &     2.530E+08,9.610E+00,4.140E+08,8.770E+00,5.840E+08,1.930E+01,
     &     5.200E+07,3.450E-01,4.010E+08,7.640E+00,5.060E+09,4.160E+00,
     &     4.420E-02,1.980E-03,1.790E+02,8.400E-03,2.230E+01,1.450E-02,
     &     0.000E+00,1.450E-02,2.020E+03,5.010E-02,7.170E+02,2.310E-02,
     &     2.200E+04,3.700E-01,2.280E+04,2.600E-01,8.010E+01,9.190E-03,
     &     2.200E+04,3.440E-01,2.290E+04,2.550E-01,6.700E+04,3.100E-01/
      DATA(siraddec(j,6),sicolstr(6,j),j=7,nlvl)/9.930E-03,5.120E-01,
     &     1.170E-02,2.670E-01,1.960E-02,1.220E-01,4.000E-01,6.020E-02,
     &     8.390E-02,3.450E-02,5.110E-01,6.080E-02,5.140E+00,1.460E-02,
     &     2.740E+00,3.220E-02,4.590E-01,2.640E-02,2.290E+09,4.320E+00,
     &     2.130E+09,1.180E+01,1.790E+09,1.600E+00,0.000E+00,1.080E-01,
     &     0.000E+00,6.390E-02,2.640E+05,3.730E-02,0.000E+00,2.230E-02,
     &     4.590E+05,1.250E-02,5.090E+06,2.630E-02,0.000E+00,1.260E-02,
     &     3.430E+07,2.280E-01,4.220E+04,2.220E-02/
      DATA(siraddec(j,7),sicolstr(7,j),j=8,nlvl)/3.000E-02,3.400E-01,
     &     9.910E-06,1.250E-01,5.530E-02,7.840E-01,1.890E-03,1.520E-01,
     &     2.920E+00,4.350E-01,0.000E+00,1.140E-02,1.640E-02,3.310E-01,
     &     0.000E+00,1.650E-02,0.000E+00,4.110E-02,0.000E+00,9.770E-02,
     &     7.710E+07,1.020E+00,2.960E+09,6.250E+01,0.000E+00,1.280E+00,
     &     4.320E+07,5.130E-01,4.450E+08,6.390E+00,2.380E+07,2.410E-01,
     &     0.000E+00,3.170E-02,2.540E+08,3.180E+00,4.990E+07,4.720E-01,
     &     0.000E+00,1.510E-02/
      DATA(siraddec(j,8),sicolstr(8,j),j=9,nlvl)/2.110E-01,3.400E-01,
     &     1.820E-01,5.440E-01,1.090E+00,1.860E-01,3.080E-03,2.610E-01,
     &     2.210E-02,2.950E-01,2.240E-01,4.260E-02,6.770E+01,2.650E-02,
     &     0.000E+00,6.460E-03,1.290E+08,1.060E+00,4.000E+07,5.880E-01,
     &     0.000E+00,1.180E-01,2.490E+09,4.090E+01,8.330E+08,9.240E+00,
     &     7.180E+04,8.910E-02,1.450E+08,1.370E+00,2.050E+08,1.200E+00,
     &     1.010E+08,1.340E+00,2.150E+08,1.950E+00,7.860E+03,1.620E-02/
      DATA(siraddec(j,9),sicolstr(9,j),j=10,nlvl)/5.600E-02,2.820E-01,
     &     4.440E-01,2.750E-01,3.480E-01,1.100E-01,1.040E+00,2.550E-02,
     &     3.130E-01,2.470E-01,1.210E+00,8.580E-03,1.860E+08,5.360E-01,
     &     9.150E+06,8.540E-02,3.180E+07,4.840E-01,0.000E+00,2.530E-02,
     &     0.000E+00,5.070E-02,2.100E+09,2.470E+01,0.000E+00,5.010E-02,
     &     1.750E+08,1.760E+00,1.060E+09,5.260E+01,0.000E+00,5.070E-02,
     &     6.450E+07,6.540E-01,2.370E+07,1.130E-01/
      DATA(siraddec(j,10),sicolstr(10,j),j=11,nlvl)/1.040E-02,3.630E-01,
     &     2.650E-01,2.370E-01,1.700E-05,2.470E-02,4.880E-01,1.450E-01,
     &     2.600E+02,9.370E-02,0.000E+00,4.450E-03,2.580E+08,2.530E+00,
     &     6.420E+08,9.980E+00,0.000E+00,1.180E-01,4.280E+07,9.680E-01,
     &     2.650E+05,9.030E-02,2.240E+09,3.980E+01,8.780E+06,1.480E-01,
     &     1.360E+07,1.310E-01,6.450E+07,1.030E+00,4.210E+08,4.580E+00,
     &     3.600E+07,2.150E-01/
      DATA(siraddec(j,11),sicolstr(11,j),j=12,nlvl)/1.260E-02,3.170E-01,
     &     7.710E-02,3.370E-02,3.340E-02,3.970E-02,4.920E+00,1.470E-02,
     &     1.770E+08,6.270E-01,7.220E+07,7.520E-01,8.690E+07,1.450E+00,
     &     0.000E+00,2.480E-02,0.000E+00,2.550E-02,1.420E+07,2.440E-01,
     &     0.000E+00,7.440E-02,2.240E+09,2.640E+01,1.450E+08,1.040E+00,
     &     0.000E+00,5.500E-02,1.180E+05,7.640E-02,9.310E+08,4.700E+00/
      DATA(siraddec(j,12),sicolstr(12,j),j=13,nlvl)/0.000E+00,1.570E-01,
     &     4.930E-03,5.820E-01,2.180E+02,1.010E-01,0.000E+00,1.930E-03,
     &     1.130E+08,1.250E+00,2.750E+08,4.790E+00,0.000E+00,4.080E-02,
     &     1.500E+07,3.670E-01,2.850E+06,7.700E-02,1.420E+08,2.870E+00,
     &     3.870E+08,4.900E+00,1.060E+03,3.680E-02,2.520E+09,4.330E+01,
     &     1.140E+08,1.480E+00,3.980E+07,2.880E-01/
      DATA(siraddec(j,13),sicolstr(13,j),j=14,nlvl)/3.880E-03,3.970E-02,
     &     0.000E+00,1.680E-03,0.000E+00,1.910E-03,1.810E+08,2.000E+00,
     &     0.000E+00,2.170E-03,0.000E+00,4.760E-03,0.000E+00,3.070E-02,
     &     0.000E+00,4.050E-03,0.000E+00,2.790E-03,0.000E+00,7.700E-03,
     &     1.280E+09,9.570E+00,0.000E+00,4.080E-02,0.000E+00,1.810E-02,
     &     6.740E+07,3.620E-01/
      DATA(siraddec(j,14),sicolstr(14,j),j=15,nlvl)/1.020E+01,1.040E-02,
     &     3.770E+08,1.460E+00,1.240E+08,1.410E+00,4.440E+07,8.050E-01,
     &     0.000E+00,3.650E-02,0.000E+00,1.030E-02,1.220E+06,4.080E-02,
     &     0.000E+00,1.800E-02,4.830E+05,9.460E-02,2.850E+08,2.220E+00,
     &     0.000E+00,9.680E-02,2.060E+09,2.260E+01,7.590E+08,4.140E+00/
      DATA(siraddec(j,15),sicolstr(15,j),j=16,nlvl)/0.000E+00,2.400E-03,
     &     2.360E+04,2.060E-02,0.000E+00,1.190E-02,0.000E+00,9.490E-03,
     &     0.000E+00,1.440E-02,0.000E+00,5.060E-03,0.000E+00,3.140E-02,
     &     0.000E+00,3.800E-03,2.470E+06,4.850E-01,0.000E+00,1.590E-02,
     &     0.000E+00,3.620E-03,1.550E+08,1.360E+01/
      DATA(siraddec(j,16),sicolstr(16,j),j=17,nlvl)/4.350E-03,6.610E-02,
     &     1.600E-06,1.330E-01,0.000E+00,2.490E-02,0.000E+00,5.070E-02,
     &     1.280E-04,3.920E-02,0.000E+00,7.770E-03,2.810E-03,6.090E-02,
     &     3.360E-01,1.950E-02,0.000E+00,5.750E-03,6.390E-03,4.760E-02,
     &     8.890E-01,1.630E-02/
      DATA(siraddec(j,17),sicolstr(17,j),j=18,nlvl)/2.920E-02,2.920E-01,
     &     0.000E+00,9.540E-02,8.920E-06,1.290E-01,4.340E-03,7.950E-02,
     &     5.840E-04,1.260E-01,1.410E-01,6.070E-02,9.030E-01,1.610E-01,
     &     5.920E-03,1.520E-01,3.360E-01,3.980E-02,1.140E+00,4.790E-02/
      DATA(siraddec(j,18),sicolstr(18,j),j=19,nlvl)/0.000E+00,2.260E-01,
     &     2.060E-04,1.060E-01,5.730E-04,1.200E-01,4.410E-02,4.670E-01,
     &     4.370E-02,1.200E-01,4.070E-01,8.250E-02,4.280E-01,3.400E-01,
     &     1.730E+00,1.770E-01,2.490E-02,6.940E-02/
      DATA(siraddec(j,19),sicolstr(19,j),j=20,nlvl)/2.190E-02,5.840E-01,
     &     3.530E-06,1.310E-01,1.670E-01,8.040E-01,4.060E-04,7.930E-02,
     &     0.000E+00,5.480E-02,2.040E+00,5.890E-01,1.680E-03,2.550E-01,
     &     0.000E+00,6.310E-02/
      DATA(siraddec(j,20),sicolstr(20,j),j=21,nlvl)/1.430E-01,7.600E-01,
     &     2.860E-02,2.530E-01,1.420E+00,3.970E-01,1.460E-03,2.330E-01,
     &     1.870E-02,1.420E-01,1.400E-01,9.380E-02,2.820E-06,8.000E-02/
      DATA(siraddec(j,21),sicolstr(21,j),j=22,nlvl)/1.500E-02,1.330E-01,
     &     6.000E-01,2.170E-01,4.650E-01,5.170E-01,7.140E-02,6.940E-02,
     &     5.600E-02,1.890E-01,5.700E-01,6.550E-02/
      DATA(siraddec(j,22),sicolstr(22,j),j=23,nlvl)/3.640E-02,1.190E-01,
     &     6.760E-06,6.910E-02,4.870E-01,2.660E-01,4.790E-01,6.170E-01,
     &     6.450E-03,5.600E-02/
      DATA(siraddec(j,23),sicolstr(23,j),j=24,nlvl)/1.830E-06,1.280E-01,
     &     9.370E-04,4.360E-01,9.120E-04,9.910E-02,6.440E-01,2.600E-01/
      DATA(siraddec(j,24),sicolstr(24,j),j=25,nlvl)/0.000E+00,7.040E-02,
     &     1.350E-02,1.860E-01,5.840E-01,2.540E-02/
      DATA(siraddec(j,25),sicolstr(25,j),j=26,nlvl)/1.520E-03,3.500E-01,
     &     7.080E-04,6.610E-02/
      DATA siraddec(nlvl,26),sicolstr(26,nlvl)/9.900E-01,2.590E-01/

      END

      BLOCK DATA ar
      PARAMETER(nlvl=27,ntrs=35)
      CHARACTER*9 alevel
      COMMON /aratom/ arraddec(nlvl,nlvl),arcolstr(nlvl,nlvl),
     &                arenergy(nlvl),keyar(nlvl),jtotar(nlvl),
     &                klowar(ntrs),kupar(ntrs)
      COMMON /arlevel/ alevel(nlvl)
c level configurations and energies (in cm**-1)

      DATA(keyar(i),alevel(i),jtotar(i),arenergy(i),i=1,nlvl)/1,
     &     '2p6    1S',0,0,2,'2p5 3s 3P',2,2026503,3,'2p5 3s 3P',1,
     &     2033140,4,'2p5 3s 3P',0,2044512,5,'2p5 3s 1P',1,2051750,6,
     &     '2p5 3p 3S',1,2149340,7,'2p5 3p 3D',3,2169872,8,'2p5 3p 3D',
     &     2,2170887,9,'2p5 3p 3D',1,2176715,10,'2p5 3p 3P',2,2182194,
     &     11,'2p5 3p 1P',1,2189188,12,'2p5 3p 3P',0,2192400,13,
     &     '2p5 3p 1D',2,2195119,14,'2p5 3p 3P',1,2196143,15,
     &     '2p5 3p 1S',0,2263884,16,'2p5 3d 3P',0,2349300,17,
     &     '2p5 3d 3P',1,2351360,18,'2p5 3d 3P',2,2355545,19,
     &     '2p5 3d 3F',4,2358765,20,'2p5 3d 3F',3,2361737,21,
     &     '2p5 3d 3F',2,2366927,22,'2p5 3d 3D',3,2370624,23,
     &     '2p5 3d 3D',1,2380781,24,'2p5 3d 3D',2,2382313,25,
     &     '2p5 3d 1F',3,2384262,26,'2p5 3d 1D',2,2384828,27,
     &     '2p5 3d 1P',1,2411016/
c most important transitions
c     data klow/5*1,3,5,6,6,6,8,11,8,9,13,7,14,10,10,2,3,3,2,4,2,3,4,5,
      DATA klowar/5*1,3,5,6,6,6,8,11,8,9,13,7,14,10,10,2,2,4,2,3,4,5,2,
     &     3,5,2,3,5,5,2,3/
c    +         14,13,12,10,14,9,10,11,14,8,9,13,7,8,11,10,6,6/
      DATA kupar/27,23,17,5,3,15,15,18,17,16,21,24,20,21,25,19,26,22,18,
     &     14,10,14,9,10,11,14,8,9,13,7,8,11,10,6,6/
c radiative decay /raddec/ rates aji (/s) and dimensionless electron
c collision strength /colstr/ at an energy of 110 rydbergs.
      DATA(arraddec(j,1),arcolstr(1,j),j=2,nlvl)/2.800E+03,4.700E-03,
     &     6.410E+10,5.720E-03,0.000E+00,9.140E-04,1.830E+11,1.070E-02,
     &     4.500E+01,9.650E-03,0.000E+00,1.180E-02,1.260E+07,8.840E-03,
     &     1.800E+02,4.480E-03,2.480E+07,7.770E-03,1.440E+02,4.270E-03,
     &     0.000E+00,1.680E-03,3.320E+07,8.930E-03,7.320E+02,3.860E-03,
     &     0.000E+00,1.250E-01,0.000E+00,4.370E-03,4.800E+09,1.310E-02,
     &     0.000E+00,2.000E-02,0.000E+00,1.610E-02,0.000E+00,1.110E-02,
     &     0.000E+00,7.700E-03,0.000E+00,7.000E-03,1.280E+11,1.080E-02,
     &     0.000E+00,6.680E-03,0.000E+00,7.980E-03,0.000E+00,7.240E-03,
     &     3.000E+12,1.690E-01/
      DATA(arraddec(j,2),arcolstr(2,j),j=3,nlvl)/5.120E+00,1.200E-01,
     &     1.340E-03,2.850E-02,1.060E+02,9.840E-02,1.250E+09,2.320E+00,
     &     2.350E+09,1.740E+01,1.030E+09,6.350E+00,3.100E+08,9.060E-01,
     &     1.520E+09,6.350E+00,4.610E+07,1.050E-01,0.000E+00,3.600E-03,
     &     1.660E+08,5.300E-01,4.430E+08,8.410E-01,0.000E+00,5.720E-03,
     &     1.440E+05,5.910E-02,1.300E+05,1.530E-01,9.080E+04,1.730E-01,
     &     1.590E+05,4.550E-01,7.260E+04,1.880E-01,4.180E+04,8.640E-02,
     &     8.700E+04,1.970E-01,2.540E+04,3.040E-02,1.080E+04,1.810E-01,
     &     2.040E+04,3.260E-02,2.880E+04,3.580E-02,3.850E+02,1.260E-02/
      DATA(arraddec(j,3),arcolstr(3,j),j=4,nlvl)/5.500E+01,3.620E-02,
     &     1.840E+01,5.650E-02,2.590E+08,1.570E+00,0.000E+00,2.930E-02,
     &     1.180E+09,7.330E+00,1.800E+09,6.020E+00,9.670E+08,4.620E+00,
     &     5.770E+07,1.520E-01,2.460E+09,1.930E+00,2.080E+08,7.490E-01,
     &     5.460E+07,1.240E-01,3.260E+09,4.710E-01,1.660E+00,3.980E-03,
     &     1.400E+04,2.760E-02,4.670E+04,1.000E-01,0.000E+00,3.400E-02,
     &     8.410E+04,2.050E-01,1.040E+05,1.720E-01,5.610E+04,1.300E-01,
     &     1.560E+05,1.210E-01,1.790E+03,1.160E-02,2.310E+04,4.640E-02,
     &     8.420E+02,1.000E-02,4.010E+04,2.510E-02/
      DATA(arraddec(j,4),arcolstr(4,j),j=5,nlvl)/2.350E+00,2.890E-02,
     &     5.410E+07,4.560E-01,0.000E+00,4.050E-04,0.000E+00,1.960E-03,
     &     1.760E+08,7.740E-01,0.000E+00,4.080E-04,1.160E+09,3.750E+00,
     &     0.000E+00,1.960E-03,0.000E+00,1.520E-02,1.020E+09,2.810E+00,
     &     0.000E+00,1.480E-03,0.000E+00,2.600E-05,1.010E-01,6.180E-04,
     &     9.310E+03,2.330E-02,0.000E+00,5.100E-05,0.000E+00,1.390E-03,
     &     1.460E+04,2.930E-02,0.000E+00,8.500E-04,8.920E-03,3.420E-03,
     &     8.190E+04,1.280E-03,0.000E+00,1.930E-02,6.080E+04,9.430E-02,
     &     1.140E+00,5.500E-03/
      DATA(arraddec(j,5),arcolstr(5,j),j=6,nlvl)/2.490E+07,2.750E-01,
     &     0.000E+00,6.690E-03,3.290E+06,4.090E-02,3.190E+06,1.150E-02,
     &     1.920E+08,2.060E+00,9.390E+08,3.670E+00,4.390E+08,5.380E-01,
     &     2.070E+09,1.130E+01,1.120E+09,3.700E+00,9.520E+09,1.770E+00,
     &     4.820E-01,9.070E-04,2.140E+03,7.650E-03,1.640E+03,1.510E-02,
     &     0.000E+00,6.460E-03,3.300E+02,8.880E-03,1.360E+02,1.050E-02,
     &     2.660E+04,7.520E-02,1.240E+04,2.220E-02,6.780E+04,1.400E-01,
     &     1.270E+05,3.070E-01,8.150E+04,1.550E-01,2.430E+05,1.450E-01/
      DATA(arraddec(j,6),arcolstr(6,j),j=7,nlvl)/2.540E-02,3.320E-01,
     &     4.530E-01,1.430E-01,3.820E-01,4.500E-02,9.110E+00,4.130E-02,
     &     3.320E+00,2.720E-02,1.350E+02,1.030E-02,9.680E+00,3.290E-02,
     &     6.600E+01,1.850E-02,3.270E+01,1.730E-02,5.200E+09,1.930E+00,
     &     4.630E+09,4.970E+00,3.450E+09,5.760E+00,0.000E+00,6.160E-02,
     &     0.000E+00,3.800E-02,8.190E+06,2.880E-02,0.000E+00,2.060E-02,
     &     1.710E+06,8.400E-03,4.600E+06,1.120E-02,0.000E+00,8.770E-03,
     &     4.660E+07,5.280E-02,6.090E+06,2.330E-02/
      DATA(arraddec(j,7),arcolstr(7,j),j=8,nlvl)/2.670E-02,2.980E-01,
     &     1.460E-04,1.300E-01,1.590E-01,4.890E-01,5.330E-03,3.400E-02,
     &     0.000E+00,7.530E-03,1.060E+02,1.720E-01,1.740E-01,1.600E-01,
     &     0.000E+00,1.060E-02,0.000E+00,1.970E-02,0.000E+00,5.000E-02,
     &     1.910E+08,4.760E-01,6.260E+09,2.290E+01,9.600E+08,2.780E+00,
     &     1.130E+08,2.590E-01,1.180E+09,2.900E+00,0.000E+00,2.240E-02,
     &     2.980E+07,5.440E-02,2.840E+08,5.490E-01,9.190E+07,1.410E-01,
     &     0.000E+00,1.640E-02/
      DATA(arraddec(j,8),arcolstr(8,j),j=9,nlvl)/3.890E+00,2.890E-01,
     &     2.940E+00,3.590E-01,5.960E+01,1.200E-01,2.080E-01,1.380E-01,
     &     2.690E-02,5.440E-02,4.020E+00,1.770E-02,2.770E+02,2.110E-02,
     &     0.000E+00,7.210E-03,3.660E+08,5.600E-01,1.960E+08,2.350E-01,
     &     0.000E+00,9.040E-02,5.440E+09,1.540E+01,2.110E+09,4.020E+00,
     &     5.510E+07,2.020E-01,4.200E+08,3.980E-01,2.500E+08,3.750E-01,
     &     1.730E+07,4.850E-02,2.980E+08,4.220E-01,4.840E+06,1.790E-02/
      DATA(arraddec(j,9),arcolstr(9,j),j=10,nlvl)/3.420E-01,2.770E-01,
     &     1.120E+01,7.730E-02,3.320E+01,1.830E-02,7.440E+00,2.800E-02,
     &     1.830E+01,7.070E-03,7.340E+01,6.260E-03,5.440E+08,3.100E-01,
     &     6.600E+07,1.180E-01,2.160E+08,5.750E-01,0.000E+00,3.110E-02,
     &     0.000E+00,3.790E-02,4.300E+09,9.050E+00,0.000E+00,5.770E+00,
     &     2.880E+09,2.830E+00,2.260E+07,4.870E-02,0.000E+00,1.540E-02,
     &     8.050E+07,1.350E-01,7.010E+08,3.950E-01/
      DATA(arraddec(j,10),arcolstr(10,j),j=11,nlvl)/4.640E-01,6.850E-02,
     &     9.760E-04,1.890E-02,1.630E+01,9.170E-02,2.520E+01,7.020E-02,
     &     3.970E+02,3.700E-02,0.000E+00,5.490E-03,6.230E+08,7.390E-01,
     &     1.710E+09,4.830E+00,0.000E+00,8.110E-02,3.730E+07,1.970E-01,
     &     1.080E+06,7.040E-02,5.070E+09,1.500E+01,2.340E+07,5.040E-02,
     &     2.700E+06,2.020E-02,6.570E+07,1.623E-01,6.850E+08,1.135E+00,
     &     6.200E+07,5.480E-02/
      DATA(arraddec(j,11),arcolstr(11,j),j=12,nlvl)/1.550E-01,3.340E-02,
     &     1.180E+00,2.700E-01,2.400E-01,2.660E-02,4.110E-01,8.380E-03,
     &     7.810E+07,5.860E-02,6.340E+07,1.370E-01,8.970E+07,3.030E-01,
     &     0.000E+00,6.810E-03,0.000E+00,9.780E-03,1.810E+07,6.240E-02,
     &     0.000E+00,1.670E-02,1.850E+05,1.160E-02,5.610E+09,1.120E+01,
     &     0.000E+00,7.640E-02,2.040E+05,5.310E-02,2.490E+09,1.660E+00/
      DATA(arraddec(j,12),arcolstr(12,j),j=13,nlvl)/2.750E-06,1.140E-01,
     &     7.680E-01,3.340E-02,0.000E+00,1.050E-03,0.000E+00,1.810E-03,
     &     2.810E+08,6.330E-01,0.000E+00,2.030E-03,0.000E+00,5.950E-03,
     &     0.000E+00,1.160E-02,0.000E+00,2.580E-03,0.000E+00,1.120E-03,
     &     2.990E+09,2.850E+00,0.000E+00,8.700E-03,0.000E+00,2.020E-02,
     &     0.000E+00,1.770E-02,1.260E+08,9.140E-02/
      DATA(arraddec(j,13),arcolstr(13,j),j=14,nlvl)/3.310E-03,3.350E-01,
     &     3.470E+02,6.060E-03,0.000E+00,6.740E-04,8.720E+07,2.150E-01,
     &     2.550E+08,9.490E-01,0.000E+00,9.310E-03,1.670E+06,2.030E-02,
     &     5.890E+03,1.870E-02,1.280E+08,4.930E-01,1.890E+07,5.060E-02,
     &     7.160E+08,1.670E+00,6.070E+09,1.760E+01,5.390E+08,1.190E+00,
     &     1.070E+08,1.200E+00/
      DATA(arraddec(j,14),arcolstr(14,j),j=15,nlvl)/1.690E+02,6.350E-03,
     &     4.620E+08,3.950E-01,1.440E+08,3.570E-01,8.110E+07,3.060E-01,
     &     0.000E+00,1.260E-02,0.000E+00,7.190E-03,6.820E+06,2.890E-02,
     &     0.000E+00,5.220E-03,5.430E+08,7.750E-01,2.830E+07,1.190E-01,
     &     0.000E+00,6.670E-02,4.850E+09,1.040E+01,1.880E+09,1.390E+00/
      DATA(arraddec(j,15),arcolstr(15,j),j=16,nlvl)/0.000E+00,2.400E-03,
     &     1.940E+04,7.520E-03,0.000E+00,1.300E-02,0.000E+00,9.970E-03,
     &     0.000E+00,9.400E-03,0.000E+00,5.000E-03,0.000E+00,1.040E-02,
     &     5.980E+06,1.060E-01,0.000E+00,4.980E-03,0.000E+00,1.190E-02,
     &     0.000E+00,5.120E-03,8.750E+08,4.580E+00/
      DATA(arraddec(j,16),arcolstr(16,j),j=17,nlvl)/1.580E-01,5.990E-02,
     &     5.860E-05,7.400E-02,0.000E+00,2.020E-02,0.000E+00,4.540E-02,
     &     1.770E-03,2.470E-02,0.000E+00,1.580E-02,1.240E+01,1.920E-02,
     &     3.570E-02,2.670E-02,0.000E+00,3.850E-03,7.890E-02,1.650E-01,
     &     3.350E+01,1.180E-02/
      DATA(arraddec(j,17),arcolstr(17,j),j=18,nlvl)/9.660E-01,1.870E-01,
     &     0.000E+00,8.040E-02,9.620E-05,9.970E-02,1.940E-01,7.480E-02,
     &     8.860E-03,1.010E-01,3.100E+01,7.830E-02,4.520E+00,3.600E-02,
     &     7.510E-02,6.010E-02,1.280E+01,2.830E-02,4.660E+01,3.640E-02/
      DATA(arraddec(j,18),arcolstr(18,j),j=19,nlvl)/0.000E+00,1.800E-01,
     &     8.240E-03,1.020E-01,9.320E-03,1.010E-01,1.910E+00,3.490E-01,
     &     1.280E+01,5.130E-02,1.590E+00,6.400E-02,1.290E+01,1.300E-01,
     &     5.950E+01,9.100E-02,3.570E-01,6.250E-02/
      DATA(arraddec(j,19),arcolstr(19,j),j=20,nlvl)/4.920E-01,4.310E-01,
     &     7.720E-05,1.570E-01,1.660E+00,5.300E-01,0.000E+00,5.780E-02,
     &     5.690E-03,4.130E-02,9.030E+01,2.360E-01,3.770E-02,1.420E-01,
     &     0.000E+00,7.770E-02/
      DATA(arraddec(j,20),arcolstr(20,j),j=21,nlvl)/2.730E+00,4.600E-01,
     &     8.970E-01,2.180E-01,2.220E-02,1.250E-01,6.270E+01,2.000E-01,
     &     5.190E-01,5.790E-02,5.230E+00,5.580E-02,3.450E-02,6.260E-02/
      DATA(arraddec(j,21),arcolstr(21,j),j=22,nlvl)/2.760E-01,1.700E-01,
     &     1.590E+01,2.570E-01,2.220E+01,1.150E-01,2.890E+00,5.160E-02,
     &     2.910E+00,9.820E-02,3.380E+01,5.480E-02/
      DATA(arraddec(j,22),arcolstr(22,j),j=23,nlvl)/7.460E-05,5.290E-02,
     &     1.700E+00,4.120E-02,2.100E+01,1.460E-01,1.990E+01,2.520E-01,
     &     6.620E-02,2.900E-02/
      DATA(arraddec(j,23),arcolstr(23,j),j=24,nlvl)/1.180E-02,8.970E-02,
     &     3.600E-07,6.840E-02,1.110E+00,1.670E-01,2.420E+01,2.540E-02/
      DATA(arraddec(j,24),arcolstr(24,j),j=25,nlvl)/1.490E-01,3.930E-01,
     &     1.530E-02,1.030E-01,1.950E+01,1.310E-01/
      DATA(arraddec(j,25),arcolstr(25,j),j=26,nlvl)/2.910E-04,3.170E-01,
     &     2.040E-02,6.040E-02/
      DATA arraddec(nlvl,26),arcolstr(26,nlvl)/4.230E+01,1.370E-01/

      END

      BLOCK DATA ti
      PARAMETER(nlvl=27,ntrs=35)
      CHARACTER*9 tlevel
      COMMON /tiatom/ tiraddec(nlvl,nlvl),ticolstr(nlvl,nlvl),
     &                tienergy(nlvl),keyti(nlvl),jtotti(nlvl),
     &                klowti(ntrs),kupti(ntrs)
      COMMON /tilevel/ tlevel(nlvl)
c level configurations and energies (in cm**-1)

      DATA(keyti(i),tlevel(i),jtotti(i),tienergy(i),i=1,nlvl)/1,
     &     '2p6    1S',0,0,2,'2p5 3s 3P',2,3698182,3,'2p5 3s 3P',1,
     &     3709200,4,'2p5 3s 3P',0,3745238,5,'2p5 3s 1P',1,3753600,6,
     &     '2p5 3p 3P',1,3879304,7,'2p5 3p 3D',2,3906182,8,'2p5 3p 3D',
     &     3,3908881,9,'2p5 3p 3D',1,3918090,10,'2p5 3p 3P',2,3926905,
     &     11,'2p5 3p 3P',0,3949860,12,'2p5 3p 1P',1,3951159,13,
     &     '2p5 3p 3P',1,3964847,14,'2p5 3p 1D',2,3965425,15,
     &     '2p5 3p 1S',0,4063885,16,'2p5 3d 3P',0,4163645,17,
     &     '2p5 3d 3P',1,4168200,18,'2p5 3d 3P',2,4176898,19,
     &     '2p5 3d 3F',4,4179494,20,'2p5 3d 3F',3,4184493,21,
     &     '2p5 3d 3F',2,4193938,22,'2p5 3d 3D',3,4199705,23,
     &     '2p5 3d 3D',1,4219800,24,'2p5 3d 3D',2,4233017,25,
     &     '2p5 3d 1D',2,4237394,26,'2p5 3d 1F',3,4238836,27,
     &     '2p5 3d 1P',1,4281600/
c most important transitions
c     data klow/7*1,3,5,6,6,7,6,12,7,9,14,10,13,8,2,10,3,2,4,3,5,5,2,
      DATA klowti/6*1,3,5,6,6,7,6,12,7,9,14,10,13,8,2,10,3,2,4,3,5,5,2,
     &     3,2,4,5,3,2,3/
c     data kup/27,23,17,10,7,5,3,15,15,18,17,21,16,24,20,21,26,22,25,
      DATA kupti/27,23,17,7,5,3,15,15,18,17,21,16,24,20,21,26,22,25,19,
     &     13,18,11,10,13,10,14,13,8,9,7,12,12,7,6,6/
c radiative decay /raddec/ rates aji (/s) and dimensionless electron
c collision strength /colstr/ at an energy of 110 rydbergs.
      DATA(tiraddec(j,1),ticolstr(1,j),j=2,nlvl)/2.800E+04,2.350E-03,
     &     3.130E+11,3.520E-03,0.000E+00,4.710E-04,4.250E+11,4.060E-03,
     &     1.650E+03,5.640E-03,1.020E+08,5.030E-03,0.000E+00,6.640E-03,
     &     3.310E+03,2.370E-03,1.380E+08,4.230E-03,0.000E+00,1.740E-03,
     &     5.630E+02,2.550E-03,1.060E+04,2.360E-03,1.890E+08,5.020E-03,
     &     0.000E+00,7.500E-02,0.000E+00,2.730E-03,2.330E+10,8.170E-03,
     &     0.000E+00,1.160E-02,0.000E+00,9.780E-03,0.000E+00,6.580E-03,
     &     0.000E+00,4.230E-03,0.000E+00,4.350E-03,9.750E+11,1.470E-03,
     &     0.000E+00,4.370E-03,0.000E+00,5.210E-03,0.000E+00,5.140E-03,
     &     1.080E+13,1.330E-01/
      DATA(tiraddec(j,2),ticolstr(2,j),j=3,nlvl)/1.740E+01,6.220E-02,
     &     6.600E-02,1.660E-02,1.730E+03,5.580E-02,2.490E+09,4.140E+00,
     &     1.810E+09,3.460E+00,4.130E+09,1.030E+01,3.850E+08,3.840E-01,
     &     2.720E+09,3.830E+00,0.000E+00,1.870E-03,1.250E+06,1.590E-03,
     &     5.190E+08,2.640E-01,1.110E+08,9.310E-02,0.000E+00,2.390E-03,
     &     2.450E+05,3.360E-02,2.340E+05,9.870E-02,1.840E+05,1.090E-01,
     &     2.790E+05,2.710E-01,1.310E+05,1.080E-01,7.850E+04,5.020E-02,
     &     1.730E+05,1.190E-01,3.370E+04,1.390E-02,6.630E+03,3.110E-03,
     &     3.410E+04,1.100E-02,1.640E+04,7.000E-03,1.320E+03,4.870E-03/
      DATA(tiraddec(j,3),ticolstr(3,j),j=4,nlvl)/1.330E+03,1.720E-02,
     &     2.900E+02,3.100E-02,2.390E+08,4.930E-01,1.850E+09,4.190E+00,
     &     0.000E+00,1.510E-02,3.500E+09,3.120E+00,1.980E+09,3.250E+00,
     &     4.140E+09,9.840E-01,2.520E+07,1.910E-02,9.960E+06,7.120E-03,
     &     1.200E+08,1.140E-01,8.540E+09,4.420E-01,1.510E+01,2.150E-03,
     &     1.430E+04,1.180E-02,6.800E+04,4.920E-02,0.000E+00,1.820E-02,
     &     1.440E+05,1.200E-01,2.000E+05,1.070E-01,1.300E+05,9.490E-02,
     &     2.920E+05,6.860E-02,2.080E+03,2.430E-03,1.000E+02,1.940E-03,
     &     1.080E+04,6.520E-03,1.390E+05,1.820E-02/
      DATA(tiraddec(j,4),ticolstr(4,j),j=5,nlvl)/5.370E+00,1.750E-02,
     &     3.660E+07,1.580E-01,0.000E+00,4.310E-04,0.000E+00,1.150E-04,
     &     4.010E+07,8.630E-02,0.000E+00,1.140E-03,0.000E+00,8.140E-04,
     &     1.800E+09,2.150E+00,2.280E+09,2.200E+00,0.000E+00,7.510E-03,
     &     0.000E+00,7.800E-04,0.000E+00,1.000E-05,4.790E-01,1.900E-04,
     &     7.570E+03,7.710E-03,0.000E+00,1.700E-05,0.000E+00,3.300E-04,
     &     7.010E+03,5.750E-03,0.000E+00,1.250E-04,2.800E-02,1.340E-03,
     &     1.510E+05,7.930E-02,1.270E+05,6.440E-02,0.000E+00,9.820E-03,
     &     1.170E+01,2.990E-03/
      DATA(tiraddec(j,5),ticolstr(5,j),j=6,nlvl)/3.220E+07,1.740E-01,
     &     3.560E+05,1.760E-02,0.000E+00,1.380E-03,3.160E+05,1.260E-02,
     &     7.440E+07,2.620E-01,1.020E+09,4.840E-01,1.740E+09,2.420E+00,
     &     1.780E+09,1.970E+00,4.060E+09,7.250E+00,1.210E+10,1.540E+00,
     &     0.000E+00,1.470E-04,3.530E+03,3.400E-03,3.140E+03,5.980E-03,
     &     0.000E+00,1.010E-03,1.310E+01,1.620E-03,1.000E+02,1.920E-03,
     &     1.180E+04,1.330E-02,3.410E+04,1.620E-03,1.220E+05,8.310E-02,
     &     1.440E+05,8.940E-02,2.720E+05,2.090E-01,4.240E+05,7.880E-02/
      DATA(tiraddec(j,6),ticolstr(6,j),j=7,nlvl)/3.580E+00,6.960E-02,
     &     5.320E-02,2.030E-01,9.660E+00,1.440E-02,7.120E+01,1.910E-02,
     &     1.550E+03,6.160E-03,1.090E+02,1.490E-02,6.710E+02,7.780E-03,
     &     9.040E+01,1.160E-02,8.380E+02,4.740E-03,8.500E+09,1.200E+00,
     &     7.220E+09,2.890E+00,4.610E+09,2.790E+00,0.000E+00,3.260E-02,
     &     0.000E+00,1.840E-02,7.930E+07,6.830E-02,0.000E+00,9.920E-03,
     &     3.080E+06,5.680E-03,1.310E+07,6.320E-03,6.220E+07,2.110E-02,
     &     0.000E+00,2.430E-03,1.380E+07,9.200E-03/
      DATA(tiraddec(j,7),ticolstr(7,j),j=8,nlvl)/1.620E-01,3.020E-01,
     &     2.100E+01,1.330E-01,2.420E+01,2.000E-01,1.930E+00,6.980E-02,
     &     1.250E+03,6.800E-02,1.930E+01,4.950E-03,3.200E-01,1.010E-02,
     &     6.270E+02,1.200E-02,0.000E+00,3.450E-03,6.830E+08,3.660E-01,
     &     5.220E+08,4.290E-01,0.000E+00,4.280E-02,8.310E+09,8.260E+00,
     &     3.530E+09,2.310E+00,9.280E+07,1.110E-01,5.610E+08,1.710E-01,
     &     1.530E+08,6.730E-02,2.640E+08,1.080E-01,4.730E+06,3.140E-03,
     &     5.500E+06,7.100E-03/
      DATA(tiraddec(j,8),ticolstr(8,j),j=9,nlvl)/1.490E-04,7.460E-02,
     &     2.250E-01,2.720E-01,0.000E+00,4.160E-03,5.450E-03,4.740E-03,
     &     1.530E+00,7.420E-02,1.620E+03,8.400E-02,0.000E+00,4.680E-03,
     &     0.000E+00,9.330E-03,0.000E+00,2.440E-02,3.430E+08,3.040E-01,
     &     9.170E+09,1.250E+01,1.450E+09,1.540E+00,1.780E+08,9.470E-01,
     &     1.950E+09,1.720E+00,0.000E+00,1.170E-02,1.510E+07,8.650E-03,
     &     9.660E+07,4.440E-02,1.920E+08,1.100E-01,0.000E+00,5.970E-03/
      DATA(tiraddec(j,9),ticolstr(9,j),j=10,nlvl)/1.430E+00,1.750E-01,
     &     3.600E+02,9.050E-03,1.450E+02,1.660E-02,4.350E+02,2.170E-02,
     &     9.130E+01,1.140E-02,9.630E+02,2.610E-03,5.770E+08,1.240E-01,
     &     7.650E+07,5.110E-02,5.720E+08,5.280E-01,0.000E+00,1.630E-02,
     &     0.000E+00,1.930E-02,6.260E+09,4.680E+00,0.000E+00,3.090E-02,
     &     4.790E+09,1.550E+00,3.340E+05,2.440E-03,4.960E+07,2.510E-02,
     &     0.000E+00,2.460E-03,1.890E+09,3.020E-01/
      DATA(tiraddec(j,10),ticolstr(10,j),j=11,nlvl)/2.550E-02,1.250E-02,
     &     4.070E+00,9.860E-03,5.290E+02,3.650E-02,3.520E+02,4.480E-02,
     &     5.820E+02,2.200E-04,0.000E+00,2.500E-03,9.360E+08,6.440E-01,
     &     2.880E+09,2.910E+00,0.000E+00,3.740E-02,6.570E+07,1.160E-01,
     &     2.570E+05,3.580E-02,7.680E+09,8.110E+00,9.410E+06,1.650E-02,
     &     4.710E+06,5.090E-03,6.350E+08,3.180E-01,2.820E+07,2.110E-02,
     &     8.100E+07,2.040E-02/
      DATA(tiraddec(j,11),ticolstr(11,j),j=12,nlvl)/1.580E-03,5.430E-03,
     &     3.880E+02,1.650E-02,2.400E-03,6.760E-02,0.000E+00,5.120E-04,
     &     0.000E+00,9.000E-04,2.520E+08,2.400E-01,0.000E+00,1.160E-03,
     &     0.000E+00,3.270E-03,0.000E+00,5.230E-03,0.000E+00,9.330E-04,
     &     0.000E+00,2.390E-03,4.560E+09,2.140E+00,0.000E+00,3.900E-03,
     &     0.000E+00,7.340E-03,0.000E+00,8.740E-03,3.230E+08,7.160E-02/
      DATA(tiraddec(j,12),ticolstr(12,j),j=13,nlvl)/6.480E+01,1.720E-01,
     &     1.990E+01,1.680E-01,8.310E+01,3.630E-03,5.240E+06,2.000E-03,
     &     2.050E+07,2.070E-02,3.860E+07,5.830E-02,0.000E+00,9.130E-04,
     &     0.000E+00,2.100E-03,2.660E+07,3.290E-02,0.000E+00,2.350E-03,
     &     1.470E+08,7.380E-02,8.760E+09,6.120E+00,2.320E+06,3.130E-02,
     &     0.000E+00,4.240E-02,4.020E+09,8.815E-01/
      DATA(tiraddec(j,13),ticolstr(13,j),j=14,nlvl)/6.010E-04,1.630E-01,
     &     1.130E+03,2.670E-03,2.880E+08,1.240E-01,8.980E+07,1.090E-01,
     &     5.790E+07,1.010E-02,0.000E+00,3.260E-03,0.000E+00,1.870E-03,
     &     7.320E+06,1.180E-02,0.000E+00,1.080E-03,4.490E+08,2.640E-01,
     &     1.480E+08,1.520E-01,7.660E+09,5.930E+00,0.000E+00,3.300E-02,
     &     2.720E+09,6.890E-01/
      DATA(tiraddec(j,14),ticolstr(14,j),j=15,nlvl)/3.090E+02,3.610E-02,
     &     0.000E+00,1.130E-04,4.470E+07,5.450E-02,1.380E+08,2.470E-01,
     &     0.000E+00,1.890E-03,7.570E+06,1.940E-02,2.440E+05,3.890E-03,
     &     4.360E+07,7.890E-02,5.030E+07,3.840E-02,1.030E+09,9.000E-01,
     &     1.010E+09,8.160E-01,9.320E+09,9.850E+00,1.810E+08,6.820E-02/
      DATA(tiraddec(j,15),ticolstr(15,j),j=16,nlvl)/0.000E+00,1.010E-03,
     &     3.260E+03,3.720E-03,0.000E+00,4.850E-03,0.000E+00,4.150E-03,
     &     0.000E+00,4.140E-03,0.000E+00,1.800E-03,0.000E+00,3.960E-03,
     &     1.580E+07,7.860E-02,0.000E+00,2.910E-03,0.000E+00,3.300E-03,
     &     0.000E+00,7.170E-03,1.930E+09,2.310E+00/
      DATA(tiraddec(j,16),ticolstr(16,j),j=17,nlvl)/1.670E+00,2.930E-02,
     &     7.120E-04,3.990E-02,0.000E+00,9.710E-03,0.000E+00,2.370E-02,
     &     6.770E-03,1.220E-02,0.000E+00,7.630E-03,1.360E+02,9.910E-03,
     &     2.120E-01,1.170E-02,4.010E-01,4.310E-03,0.000E+00,8.760E-04,
     &     4.880E+02,7.630E-03/
      DATA(tiraddec(j,17),ticolstr(17,j),j=18,nlvl)/9.070E+00,9.480E-02,
     &     0.000E+00,3.950E-02,2.680E-04,5.400E-02,1.850E+00,4.180E-02,
     &     2.810E-02,5.300E-02,3.210E+02,3.290E-02,6.410E+01,1.690E-02,
     &     1.910E+02,1.020E-02,4.670E-01,2.100E-02,7.770E+02,1.620E-02/
      DATA(tiraddec(j,18),ticolstr(18,j),j=19,nlvl)/0.000E+00,9.390E-02,
     &     4.130E-02,6.080E-02,8.930E-03,6.500E-02,1.490E+01,2.020E-01,
     &     1.340E+02,2.500E-02,3.960E+01,1.870E-02,8.030E+02,3.650E-02,
     &     1.630E+02,4.190E-02,4.390E+01,2.550E-02/
      DATA(tiraddec(j,19),ticolstr(19,j),j=20,nlvl)/1.820E+00,2.080E-01,
     &     3.380E-04,8.950E-02,5.170E+00,2.830E-01,0.000E+00,3.380E-02,
     &     3.140E-02,8.610E-03,4.270E-01,7.000E-02,1.470E+03,1.090E-01,
     &     0.000E+00,2.410E-02/
      DATA(tiraddec(j,20),ticolstr(20,j),j=21,nlvl)/1.230E+01,2.480E-01,
     &     6.160E+00,1.080E-01,1.230E-01,5.890E-02,1.160E+03,9.370E-02,
     &     5.360E+01,1.540E-02,4.230E+00,1.340E-02,1.090E+00,3.260E-02/
      DATA(tiraddec(j,21),ticolstr(21,j),j=22,nlvl)/1.160E+00,9.100E-02,
     &     1.630E+02,1.330E-01,3.590E+02,4.700E-02,7.390E+01,3.250E-02,
     &     4.000E+01,1.550E-02,6.840E+02,2.430E-02/
      DATA(tiraddec(j,22),ticolstr(22,j),j=23,nlvl)/4.850E-04,2.540E-02,
     &     2.380E+01,8.330E-03,3.910E+02,1.020E-01,4.080E+02,7.070E-02,
     &     1.400E-01,1.140E-02/
      DATA(tiraddec(j,23),ticolstr(23,j),j=24,nlvl)/3.970E+00,4.000E-02,
     &     5.490E+01,7.570E-02,2.940E-04,3.680E-02,2.820E+02,1.220E-02/
      DATA(tiraddec(j,24),ticolstr(24,j),j=25,nlvl)/2.510E-02,6.350E-02,
     &     2.580E+00,2.030E-01,1.070E+02,6.970E-02/
      DATA(tiraddec(j,25),ticolstr(25,j),j=26,nlvl)/8.850E-03,1.810E-01,
     &     3.900E+02,7.300E-02/
      DATA tiraddec(nlvl,26),ticolstr(26,nlvl)/5.700E-02,3.250E-02/

      END

      BLOCK DATA fe
      PARAMETER(nlvl=27,ntrs=35)
      CHARACTER*9 flevel
      COMMON /featom/ feraddec(nlvl,nlvl),fecolstr(nlvl,nlvl),
     &                feenergy(nlvl),keyfe(nlvl),jtotfe(nlvl),
     &                klowfe(ntrs),kupfe(ntrs)
      COMMON /felevel/ flevel(nlvl)
c level configurations and energies (in cm**-1)

      DATA(keyfe(i),flevel(i),jtotfe(i),feenergy(i),i=1,nlvl)/1,
     &     '2p6    1S',0,0,2,'2p5 3s 3P',2,5849320,3,'2p5 3s 1P',1,
     &     5864590,4,'2p5 3s 3P',0,5951212,5,'2p5 3s 3P',1,5960870,6,
     &     '2p5 3p 3S',1,6093407,7,'2p5 3p 3D',2,6121606,8,'2p5 3p 3D',
     &     3,6134630,9,'2p5 3p 3P',1,6143730,10,'2p5 3p 3P',2,6158360,
     &     11,'2p5 3p 3P',0,6202450,12,'2p5 3p 1P',1,6219114,13,
     &     '2p5 3p 3D',1,6245225,14,'2p5 3p 1D',2,6248350,15,
     &     '2p5 3p 1S',0,6353230,16,'2p5 3d 3P',0,6463942,17,
     &     '2p5 3d 3P',1,6472500,18,'2p5 3d 3F',4,6486530,19,
     &     '2p5 3d 3P',2,6486288,20,'2p5 3d 3F',3,6492788,21,
     &     '2p5 3d 3D',2,6506650,22,'2p5 3d 3D',3,6515320,23,
     &     '2p5 3d 3D',1,6552200,24,'2p5 3d 3P',2,6594461,25,
     &     '2p5 3d 1D',2,6606500,26,'2p5 3d 1F',3,6606500,27,
     &     '2p5 3d 1P',1,6660000/
c most important transitions
      DATA klowfe/8*1,3,6,5,7,12,7,6,9,13,14,10,8,3,10,2,4,3,5,2,5,3,2,
     &     4,5,3,2,3/
      DATA kupfe/27,23,17,10,7,5,3,2,15,19,15,21,24,20,16,21,25,26,22,
     &     18,11,19,10,13,10,14,8,13,9,7,12,12,7,6,6/
c radiative decay /raddec/ rates aji (/s) and dimensionless electron
c collision strength /colstr/ at an energy of 110 rydbergs.
      DATA(feraddec(j,1),fecolstr(1,j),j=2,nlvl)/2.480E+05,1.460E-03,
     &     8.490E+11,2.190E-03,0.000E+00,2.910E-04,8.400E+11,2.040E-03,
     &     2.490E+04,3.460E-03,4.480E+08,3.230E-03,0.000E+00,4.240E-03,
     &     2.240E+04,1.540E-03,5.280E+08,2.670E-03,0.000E+00,2.520E-03,
     &     9.860E+02,1.670E-03,8.230E+04,1.640E-03,7.150E+08,3.200E-03,
     &     0.000E+00,4.810E-02,0.000E+00,1.840E-03,6.860E+10,5.420E-03,
     &     0.000E+00,6.500E-03,0.000E+00,7.130E-03,0.000E+00,4.310E-03,
     &     0.000E+00,2.670E-03,0.000E+00,2.920E-03,5.000E+12,2.060E-02,
     &     0.000E+00,3.050E-03,0.000E+00,4.050E-03,0.000E+00,3.590E-03,
     &     2.590E+13,9.760E-02/
      DATA(feraddec(j,2),fecolstr(2,j),j=3,nlvl)/3.690E+01,3.780E-02,
     &     1.480E+00,1.100E-02,1.670E+04,3.640E-02,4.160E+09,3.000E+00,
     &     2.650E+09,2.400E+00,6.490E+09,7.080E+00,2.050E+08,9.360E-02,
     &     4.220E+09,2.550E+00,0.000E+00,1.280E-03,5.610E+06,1.340E-03,
     &     5.160E+08,8.180E-02,8.150E+07,2.140E-02,0.000E+00,1.150E-03,
     &     3.650E+05,2.180E-02,3.660E+05,6.140E-02,4.280E+05,1.900E-01,
     &     3.000E+05,8.130E-02,2.000E+05,7.400E-02,1.030E+05,2.980E-02,
     &     2.810E+05,8.340E-02,3.480E+04,7.770E-03,3.200E+03,5.930E-04,
     &     3.340E+04,3.580E-03,1.310E+04,1.920E-03,3.520E+03,2.050E-03/
      DATA(feraddec(j,3),fecolstr(3,j),j=4,nlvl)/1.550E+04,9.860E-03,
     &     2.920E+03,2.000E-02,1.200E+08,1.090E-01,2.540E+09,2.790E+00,
     &     0.000E+00,9.260E-03,5.700E+09,2.860E+00,3.260E+09,2.330E+00,
     &     7.380E+09,6.510E-01,1.680E+07,4.160E-03,7.410E+05,5.040E-04,
     &     7.410E+07,2.170E-02,1.380E+10,3.020E-01,8.230E+01,1.310E-03,
     &     1.320E+04,6.100E-03,0.000E+00,1.140E-02,8.270E+04,2.700E-02,
     &     2.170E+05,8.170E-02,3.310E+05,7.750E-02,2.100E+05,6.670E-02,
     &     5.220E+05,4.880E-02,1.730E+03,5.850E-04,2.610E+02,5.170E-04,
     &     5.230E+03,1.200E-03,2.310E+05,1.000E-02/
      DATA(feraddec(j,4),fecolstr(4,j),j=5,nlvl)/9.320E+00,1.180E-02,
     &     1.490E+07,5.900E-02,0.000E+00,1.360E-04,0.000E+00,4.500E-05,
     &     5.710E+06,9.770E-03,0.000E+00,4.600E-05,0.000E+00,3.930E-04,
     &     2.190E+09,1.250E+00,4.110E+09,1.720E+00,0.000E+00,4.430E-03,
     &     0.000E+00,5.300E-04,0.000E+00,6.000E-06,1.440E+00,6.500E-05,
     &     0.000E+00,7.000E-06,4.410E+03,2.690E-03,0.000E+00,9.200E-05,
     &     2.320E+03,1.170E-03,0.000E+00,2.900E-05,2.720E+00,5.430E-04,
     &     2.130E+05,4.950E-02,2.270E+05,4.930E-02,0.000E+00,5.880E-03,
     &     6.740E+01,2.030E-03/
      DATA(feraddec(j,5),fecolstr(5,j),j=6,nlvl)/1.560E+07,7.790E-02,
     &     8.260E+03,3.920E-04,0.000E+00,4.270E-04,1.530E+06,3.340E-03,
     &     2.370E+07,6.280E-02,1.210E+09,3.260E-01,2.800E+09,1.810E+00,
     &     2.440E+09,1.150E+00,6.570E+09,4.410E+00,1.480E+10,6.480E-01,
     &     7.290E+00,3.600E-05,2.730E+03,1.480E-03,0.000E+00,2.380E-04,
     &     2.380E+03,2.290E-03,0.000E+00,4.810E-04,9.040E+01,4.740E-04,
     &     4.440E+03,3.040E-03,3.720E+04,8.830E-03,2.030E+05,5.930E-02,
     &     2.060E+05,5.660E-02,4.410E+05,1.460E-01,6.500E+05,5.410E-02/
      DATA(feraddec(j,6),fecolstr(6,j),j=7,nlvl)/1.090E+01,3.970E-02,
     &     1.570E-01,1.420E-01,1.010E+02,7.070E-03,3.610E+02,1.370E-02,
     &     1.120E+04,4.830E-03,1.720E+03,1.070E-02,3.860E+03,3.480E-03,
     &     5.130E+02,4.580E-03,1.130E+04,4.220E-03,1.240E+10,8.390E-01,
     &     1.020E+10,1.920E+00,0.000E+00,2.040E-02,5.560E+09,1.540E+00,
     &     0.000E+00,1.060E-02,3.630E+08,8.980E-02,0.000E+00,6.640E-03,
     &     8.240E+07,1.220E-02,2.210E+07,3.340E-03,6.970E+07,8.390E-03,
     &     0.000E+00,7.140E-04,2.840E+06,3.010E-03/
      DATA(feraddec(j,7),fecolstr(7,j),j=8,nlvl)/2.030E+01,1.080E-01,
     &     1.000E+02,9.000E-02,1.560E+02,1.290E-01,1.640E+01,4.300E-02,
     &     1.360E+04,4.320E-02,3.200E+01,1.420E-03,2.620E+00,2.430E-03,
     &     1.140E+03,6.640E-03,0.000E+00,2.060E-03,1.100E+09,2.600E-01,
     &     0.000E+00,2.850E-02,1.140E+09,3.940E-01,1.170E+10,5.260E+00,
     &     4.890E+09,1.400E+00,9.930E+07,5.830E-02,6.810E+08,8.580E-02,
     &     8.100E+07,1.250E-02,2.220E+08,3.100E-02,4.580E+06,1.690E-03,
     &     7.470E+05,2.780E-03/
      DATA(feraddec(j,8),fecolstr(8,j),j=9,nlvl)/2.850E-05,4.220E-02,
     &     2.330E-01,1.770E-01,0.000E+00,2.910E-03,2.220E-04,1.020E-03,
     &     1.190E+01,3.970E-02,1.500E+04,5.080E-02,0.000E+00,2.300E-03,
     &     0.000E+00,5.450E-03,0.000E+00,1.460E-02,1.210E+10,8.240E+00,
     &     5.390E+08,2.210E-01,1.900E+09,9.850E-01,2.030E+08,8.030E-02,
     &     2.710E+09,1.130E+00,0.000E+00,7.270E-03,5.890E+06,1.420E-03,
     &     8.770E+07,1.460E-02,1.320E+08,2.780E-02,0.000E+00,2.500E-03/
      DATA(feraddec(j,9),fecolstr(9,j),j=10,nlvl)/4.720E+00,1.130E-01,
     &     1.850E+03,4.720E-03,9.650E+02,6.810E-03,5.660E+03,1.170E-02,
     &     8.910E+02,9.120E-03,5.950E+03,1.350E-03,2.730E+08,2.980E-02,
     &     9.490E+06,6.000E-03,0.000E+00,9.530E-03,1.030E+09,4.310E-01,
     &     0.000E+00,1.200E-02,8.440E+09,2.930E+00,0.000E+00,1.760E-02,
     &     7.340E+09,1.010E+00,1.940E+06,8.740E-04,2.340E+07,4.500E-03,
     &     0.000E+00,5.950E-04,2.650E+09,1.600E-01/
      DATA(feraddec(j,10),fecolstr(10,j),j=11,nlvl)/3.470E-01,9.820E-03,
     &     8.540E+00,2.110E-03,6.050E+03,2.230E-02,3.940E+03,2.880E-02,
     &     8.750E+02,8.490E-03,0.000E+00,1.470E-03,1.290E+09,4.070E-01,
     &     0.000E+00,2.250E-02,4.060E+09,1.920E+00,1.150E+08,9.280E-02,
     &     2.570E+07,3.130E-02,1.030E+10,5.190E+00,7.200E+06,8.570E-03,
     &     9.780E+06,2.480E-03,5.080E+08,9.090E-02,1.360E+07,3.840E-03,
     &     9.140E+07,8.460E-03/
      DATA(feraddec(j,11),fecolstr(11,j),j=12,nlvl)/5.510E-01,2.330E-03,
     &     7.990E+02,9.620E-03,1.760E-01,4.130E-02,7.920E+01,2.400E-03,
     &     0.000E+00,6.120E-04,1.790E+08,9.870E-02,0.000E+00,2.330E-03,
     &     0.000E+00,1.090E-03,0.000E+00,3.260E-03,0.000E+00,6.290E-04,
     &     0.000E+00,1.820E-03,6.900E+09,1.370E+00,0.000E+00,1.940E-03,
     &     0.000E+00,3.490E-03,0.000E+00,4.280E-03,7.090E+08,6.380E-02/
      DATA(feraddec(j,12),fecolstr(12,j),j=13,nlvl)/8.870E+00,1.100E-01,
     &     1.780E+02,1.130E-01,0.000E+00,0.000E+00,3.710E+04,6.400E-05,
     &     5.140E+06,3.790E-03,0.000E+00,2.280E-02,1.450E+07,1.470E-02,
     &     0.000E+00,6.500E-04,1.010E+07,8.150E-03,0.000E+00,6.380E-04,
     &     2.180E+08,5.950E-02,1.210E+10,3.790E+00,8.490E+06,2.100E-02,
     &     0.000E+00,2.590E-02,6.140E+09,6.120E-01/
      DATA(feraddec(j,13),fecolstr(13,j),j=14,nlvl)/7.860E-02,9.320E-02,
     &     3.820E+03,1.760E-03,1.190E+08,4.240E-02,3.710E+07,3.570E-02,
     &     0.000E+00,9.500E-04,2.600E+07,3.380E-02,0.000E+00,5.570E-04,
     &     2.850E+06,3.310E-03,0.000E+00,3.210E-04,2.420E+08,8.840E-02,
     &     2.630E+08,1.240E-01,1.060E+10,3.890E+00,0.000E+00,1.980E-02,
     &     3.320E+09,4.090E-01/
      DATA(feraddec(j,14),fecolstr(14,j),j=15,nlvl)/1.890E+02,2.710E-02,
     &     0.000E+00,2.800E-05,1.860E+07,1.890E-02,0.000E+00,5.930E-04,
     &     5.380E+07,7.360E-02,5.050E+06,9.540E-03,5.430E+05,1.670E-03,
     &     1.290E+07,5.620E-03,5.770E+07,2.460E-02,1.430E+09,6.020E-01,
     &     1.350E+09,5.280E-01,1.250E+10,6.360E+00,2.760E+08,4.920E-02/
      DATA(feraddec(j,15),fecolstr(15,j),j=16,nlvl)/0.000E+00,4.890E-04,
     &     6.460E+03,1.560E-03,0.000E+00,2.030E-03,0.000E+00,2.230E-03,
     &     0.000E+00,2.090E-03,0.000E+00,8.090E-04,0.000E+00,1.810E-03,
     &     3.040E+07,6.770E-02,0.000E+00,2.180E-03,0.000E+00,2.740E-03,
     &     0.000E+00,2.400E-04,3.370E+09,1.400E+00/
      DATA(feraddec(j,16),fecolstr(16,j),j=17,nlvl)/8.770E+00,1.730E-02,
     &     0.000E+00,5.780E-03,4.070E-03,2.560E-02,0.000E+00,1.460E-02,
     &     1.200E-02,6.660E-03,0.000E+00,4.520E-03,7.630E+02,6.240E-03,
     &     1.310E+00,6.680E-03,1.730E+00,1.320E-03,0.000E+00,2.400E-04,
     &     4.370E+03,3.090E-03/
      DATA(feraddec(j,17),fecolstr(17,j),j=18,nlvl)/0.000E+00,2.350E-02,
     &     4.390E+01,6.020E-02,4.400E-04,3.470E-02,4.740E+00,2.570E-02,
     &     5.820E-02,3.310E-02,1.720E+03,1.790E-02,6.610E+02,1.110E-02,
     &     1.600E+03,4.520E-03,2.680E+00,9.140E-03,7.740E+03,8.760E-03/
      DATA(feraddec(j,18),fecolstr(18,j),j=19,nlvl)/0.000E+00,5.800E-02,
     &     2.050E+00,1.430E-01,5.040E-04,5.860E-02,1.160E+01,1.850E-01,
     &     0.000E+00,2.470E-02,1.230E-01,2.100E-03,4.100E+00,4.200E-02,
     &     1.380E+04,6.520E-02,0.000E+00,1.000E-02/
      DATA(feraddec(j,19),fecolstr(19,j),j=20,nlvl)/4.650E+00,4.490E-02,
     &     1.280E-01,4.470E-02,5.690E+01,1.350E-01,7.510E+02,1.620E-02,
     &     6.720E+02,8.850E-03,6.640E+03,1.850E-02,1.340E+03,1.980E-02,
     &     1.380E+03,1.210E-02/
      DATA(feraddec(j,20),fecolstr(20,j),j=21,nlvl)/3.080E+01,1.570E-01,
     &     2.480E+01,6.760E-02,4.580E-01,3.200E-02,1.210E+04,5.500E-02,
     &     2.520E+02,4.650E-03,2.530E+01,3.740E-03,1.190E+01,2.150E-02/
      DATA(feraddec(j,21),fecolstr(21,j),j=22,nlvl)/4.040E+00,5.910E-02,
     &     8.270E+02,8.260E-02,3.460E+03,2.600E-02,1.230E+03,1.410E-02,
     &     4.120E+02,8.260E-03,6.720E+03,1.230E-02/
      DATA(feraddec(j,22),fecolstr(22,j),j=23,nlvl)/1.780E-03,1.600E-02,
     &     1.470E+02,2.260E-03,4.460E+03,5.180E-02,4.430E+03,4.340E-02,
     &     2.480E-01,4.700E-03/
      DATA(feraddec(j,23),fecolstr(23,j),j=24,nlvl)/9.120E+01,1.960E-02,
     &     1.210E+03,3.820E-02,2.910E-02,2.280E-02,1.940E+03,6.650E-03/
      DATA(feraddec(j,24),fecolstr(24,j),j=25,nlvl)/2.800E-04,4.270E-02,
     &     1.770E+01,1.260E-01,2.540E+02,5.010E-02/
      DATA(feraddec(j,25),fecolstr(25,j),j=26,nlvl)/2.140E-01,1.150E-01,
     &     1.570E+03,5.010E-01/
      DATA feraddec(nlvl,26),fecolstr(26,nlvl)/8.760E-02,2.240E-02/

      END

      BLOCK DATA ge
      PARAMETER(nlvl=27,ntrs=35)
      CHARACTER*9 glevel
      COMMON /geatom/ geraddec(nlvl,nlvl),gecolstr(nlvl,nlvl),
     &                geenergy(nlvl),keyge(nlvl),jtotge(nlvl),
     &                klowge(ntrs),kupge(ntrs)
      COMMON /gelevel/ glevel(nlvl)
c level configurations and energies (in cm**-1)

      DATA(keyge(i),glevel(i),jtotge(i),geenergy(i),i=1,nlvl)/1,
     &     '2p6    1S',0,0,2,'2p5 3s 3P',2,10005993,3,'2p5 3s 1P',1,
     &     10027116,4,'2p5 3s 3P',0,10262398,5,'2p5 3s 3P',1,10274238,6,
     &     '2p5 3p 3P',1,10354527,7,'2p5 3p 1D',2,10376205,8,
     &     '2p5 3p 3D',3,10422831,9,'2p5 3p 3S',1,10431450,10,
     &     '2p5 3p 3P',2,10457705,11,'2p5 3p 3P',0,10546284,12,
     &     '2p5 3p 1P',1,10624977,13,'2p5 3p 3D',1,10688815,14,
     &     '2p5 3p 3D',2,10697501,15,'2p5 3p 1S',0,10784286,16,
     &     '2p5 3d 3P',0,10867746,17,'2p5 3d 3P',1,10882567,18,
     &     '2p5 3d 3F',4,10907722,19,'2p5 3d 3F',3,10907789,20,
     &     '2p5 3d 3P',2,10909057,21,'2p5 3d 3D',2,10929675,22,
     &     '2p5 3d 3D',3,10948407,23,'2p5 3d 3D',1,11018449,24,
     &     '2p5 3d 3F',2,11160400,25,'2p5 3d 1D',2,11178566,26,
     &     '2p5 3d 1F',3,11186849,27,'2p5 3d 1P',1,11259601/
c most important transitions
c     data klow/9*1,3,6,7,12,7,7,3,6,5,9,10,14,13,8,10,2,3,4,5,
      DATA klowge/9*1,3,6,12,7,7,3,6,5,9,10,14,13,8,10,2,3,4,5,2,5,3,2,
     &     4,5,3,2/
c     data kup/27,23,17,14,10,7,5,3,2,15,20,21,24,20,19,11,16,15,
      DATA kupge/27,23,17,14,10,7,5,3,2,15,20,24,20,19,11,16,15,21,22,
     &     26,25,18,20,10,10,13,14,8,13,9,7,12,12,7,6/
c radiative decay /raddec/ rates aji (/s) and dimensionless electron
c collision strength /colstr/ at an energy of 110 rydbergs.
      DATA(geraddec(j,1),gecolstr(1,j),j=2,nlvl)/4.000E+07,9.030E-04,
     &     2.360E+12,1.130E-03,0.000E+00,1.800E-04,2.060E+12,9.760E-04,
     &     4.980E+05,1.930E-03,2.290E+09,2.040E-03,0.000E+00,2.680E-03,
     &     1.350E+05,1.110E-03,2.550E+09,1.650E-03,0.000E+00,4.370E-03,
     &     1.210E+03,1.090E-03,2.290E+05,1.130E-03,3.550E+09,1.980E-03,
     &     0.000E+00,2.560E-02,0.000E+00,1.230E-03,1.630E+11,3.450E-03,
     &     0.000E+00,4.340E-03,0.000E+00,2.790E-03,0.000E+00,4.160E-03,
     &     0.000E+00,1.890E-03,0.000E+00,1.920E-03,2.880E+13,2.550E-02,
     &     0.000E+00,2.120E-03,0.000E+00,3.100E-03,0.000E+00,2.440E-03,
     &     6.450E+13,5.170E-02/
      DATA(geraddec(j,2),gecolstr(2,j),j=3,nlvl)/8.340E+01,2.320E-02,
     &     6.190E+01,6.900E-03,2.630E+05,2.300E-02,2.980E+09,1.930E+00,
     &     4.020E+09,1.520E+00,1.200E+10,4.300E+00,4.110E+08,7.040E-03,
     &     7.580E+09,1.530E+00,0.000E+00,1.080E-03,2.090E+07,9.440E-04,
     &     5.190E+08,1.600E-02,6.600E+07,3.420E-03,0.000E+00,4.580E-04,
     &     6.230E+05,1.350E-02,6.570E+05,3.900E-02,7.680E+05,1.170E-01,
     &     3.400E+05,4.560E-02,6.970E+05,5.930E-02,7.320E+04,1.160E-02,
     &     5.220E+05,5.340E-02,3.000E+04,4.520E-03,1.360E+03,8.500E-05,
     &     3.340E+04,7.380E-04,1.170E+04,3.600E-04,9.550E+03,6.540E-04/
      DATA(geraddec(j,3),gecolstr(3,j),j=4,nlvl)/2.770E+05,5.690E-03,
     &     4.860E+04,1.240E-02,5.250E+08,4.600E-02,3.660E+09,1.690E+00,
     &     0.000E+00,6.060E-03,1.050E+10,1.800E+00,6.190E+09,1.470E+00,
     &     1.870E+10,4.640E-01,1.280E+07,6.870E-04,1.520E+05,9.900E-05,
     &     5.100E+07,2.770E-03,1.920E+10,1.270E-01,6.290E+02,8.460E-04,
     &     1.090E+04,3.170E-03,0.000E+00,7.370E-03,3.800E+05,5.130E-02,
     &     5.920E+04,9.190E-03,6.829E+05,5.570E-02,3.780E+05,4.130E-02,
     &     1.150E+06,3.430E-02,1.460E+03,1.110E-04,4.860E+02,1.100E-04,
     &     2.820E+03,1.760E-04,3.140E+05,3.120E-03/
      DATA(geraddec(j,4),gecolstr(4,j),j=5,nlvl)/1.810E+01,7.920E-03,
     &     7.840E+05,1.650E-02,0.000E+00,4.600E-05,0.000E+00,1.700E-05,
     &     4.630E+04,1.160E-04,0.000E+00,1.900E-05,0.000E+00,1.530E-04,
     &     2.620E+09,6.880E-01,8.150E+09,1.160E+00,0.000E+00,2.880E-03,
     &     0.000E+00,4.410E-04,0.000E+00,4.000E-06,5.180E+00,1.700E-05,
     &     0.000E+00,3.000E-06,0.000E+00,2.200E-05,1.460E+03,7.030E-04,
     &     2.190E+02,1.030E-04,0.000E+00,7.000E-06,7.960E+01,1.540E-04,
     &     3.270E+05,2.790E-02,4.600E+05,3.490E-02,0.000E+00,3.810E-03,
     &     4.760E+02,1.520E-03/
      DATA(geraddec(j,5),gecolstr(5,j),j=6,nlvl)/6.500E+05,2.140E-02,
     &     3.530E+02,9.700E-05,0.000E+00,1.440E-04,1.580E+06,6.120E-03,
     &     2.920E+06,1.190E-02,6.490E+08,1.360E-01,4.490E+09,1.220E+00,
     &     4.060E+09,6.380E-01,1.230E+10,3.050E+00,2.030E+10,4.790E-01,
     &     0.000E+00,1.000E-05,9.490E+02,4.380E-04,0.000E+00,6.500E-05,
     &     4.700E+01,1.270E-04,7.550E+02,5.640E-04,2.080E+02,1.720E-04,
     &     9.770E+02,5.120E-04,1.770E+04,2.680E-03,3.880E+05,4.040E-02,
     &     3.320E+05,3.300E-02,8.070E+05,9.160E-02,1.140E+06,3.680E-02/
      DATA(geraddec(j,6),gecolstr(6,j),j=7,nlvl)/2.180E+01,2.300E-02,
     &     1.690E+00,9.470E-02,1.650E+03,6.170E-03,3.050E+03,1.500E-02,
     &     1.080E+05,3.740E-03,3.970E+04,9.200E-03,2.370E+04,9.600E-04,
     &     3.170E+03,1.210E-03,2.550E+05,2.830E-03,1.980E+10,5.030E-01,
     &     1.580E+10,1.080E+00,0.000E+00,1.290E-02,0.000E+00,6.300E-03,
     &     7.780E+09,7.960E-01,1.190E+09,1.100E-01,0.000E+00,5.470E-03,
     &     8.580E+08,3.220E-02,3.130E+07,1.060E-03,6.960E+07,1.960E-03,
     &     0.000E+00,1.450E-04,2.680E+07,1.120E-03/
      DATA(geraddec(j,7),gecolstr(7,j),j=8,nlvl)/1.020E+03,6.960E-02,
     &     9.920E+02,5.930E-02,1.870E+03,7.990E-02,2.510E+02,2.640E-02,
     &     2.350E+05,2.620E-02,9.110E+00,2.900E-04,1.980E+01,4.440E-04,
     &     2.230E+03,3.160E-03,0.000E+00,1.350E-03,2.020E+09,1.610E-01,
     &     0.000E+00,1.780E-02,1.870E+10,3.060E+00,3.660E+09,4.310E-01,
     &     6.260E+09,6.460E-01,9.410E+07,2.570E-02,8.640E+08,3.750E-02,
     &     4.200E+07,1.580E-03,1.860E+08,5.760E-03,7.680E+06,4.280E-04,
     &     7.680E+06,8.960E-04/
      DATA(geraddec(j,8),gecolstr(8,j),j=9,nlvl)/1.000E-05,2.470E-02,
     &     1.920E-01,1.090E-01,0.000E+00,2.300E-03,2.300E-01,2.350E-04,
     &     1.990E+02,2.030E-02,2.300E+05,3.030E-02,0.000E+00,8.520E-04,
     &     0.000E+00,3.170E-03,0.000E+00,8.730E-03,1.690E+10,4.860E+00,
     &     2.510E+09,5.770E-01,9.510E+08,1.620E-01,1.190E+08,2.160E-02,
     &     3.970E+09,6.890E-01,0.000E+00,5.500E-03,1.530E+06,1.690E-04,
     &     7.700E+07,4.620E-03,9.270E+07,4.980E-03,0.000E+00,7.770E-04/
      DATA(geraddec(j,9),gecolstr(9,j),j=10,nlvl)/1.200E+01,6.580E-02,
     &     6.410E+03,2.140E-03,6.570E+03,2.150E-03,1.160E+05,7.510E-03,
     &     1.580E+04,7.760E-03,3.880E+04,6.350E-04,1.500E+08,1.500E-03,
     &     1.690E+08,2.180E-02,0.000E+00,6.470E-03,0.000E+00,7.810E-03,
     &     1.230E+09,2.850E-01,1.250E+10,1.800E+00,0.000E+00,9.710E-03,
     &     1.220E+10,5.810E-01,4.020E+05,1.250E-04,5.300E+06,4.000E-04,
     &     0.000E+00,1.610E-04,2.860E+09,4.760E-02/
      DATA(geraddec(j,10),gecolstr(10,j),j=11,nlvl)/5.870E+00,8.260E-03,
     &     2.330E+00,3.800E-04,1.100E+05,1.320E-02,6.950E+04,1.850E-02,
     &     1.770E+03,3.160E-03,0.000E+00,9.760E-04,1.620E+09,2.350E-01,
     &     0.000E+00,1.390E-02,2.170E+08,7.250E-02,5.500E+09,1.080E+00,
     &     4.700E+08,9.210E-02,1.460E+10,3.060E+00,3.200E+07,7.200E-03,
     &     1.530E+07,8.460E-04,3.840E+08,1.660E-02,6.690E+06,5.210E-04,
     &     8.710E+07,2.200E-03/
      DATA(geraddec(j,11),gecolstr(11,j),j=12,nlvl)/5.600E+00,7.900E-04,
     &     2.520E+04,5.400E-03,1.650E+01,2.020E-02,0.000E+00,1.390E-04,
     &     0.000E+00,5.050E-04,8.560E+07,1.190E-02,0.000E+00,2.030E-03,
     &     0.000E+00,2.310E-03,0.000E+00,1.170E-03,0.000E+00,6.320E-04,
     &     0.000E+00,1.460E-03,8.070E+09,6.660E-01,0.000E+00,1.170E-03,
     &     0.000E+00,1.340E-03,0.000E+00,1.560E-03,1.360E+09,3.700E-02/
      DATA(geraddec(j,12),gecolstr(12,j),j=13,nlvl)/3.970E+02,6.650E-02,
     &     2.730E+03,7.350E-02,3.560E+01,1.900E-03,1.390E+05,6.480E-05,
     &     4.310E+05,4.450E-04,0.000E+00,5.400E-05,0.000E+00,1.860E-04,
     &     3.220E+06,3.330E-03,9.800E+05,8.560E-04,0.000E+00,1.650E-04,
     &     9.210E+07,1.790E-02,1.910E+10,2.190E+00,4.350E+06,1.200E-02,
     &     0.000E+00,1.650E-02,1.030E+10,4.070E-01/
      DATA(geraddec(j,13),gecolstr(13,j),j=14,nlvl)/9.910E-01,5.410E-02,
     &     7.070E+03,1.680E-03,1.340E+07,1.120E-02,4.400E+06,8.620E-03,
     &     0.000E+00,2.110E-04,0.000E+00,1.310E-04,3.720E+06,7.930E-03,
     &     7.590E+04,2.110E-04,0.000E+00,9.300E-05,4.370E+07,1.720E-02,
     &     3.190E+08,6.720E-02,1.500E+10,5.130E-01,0.000E+00,1.270E-02,
     &     4.150E+09,2.380E-01/
      DATA(geraddec(j,14),gecolstr(14,j),j=15,nlvl)/3.230E+01,2.370E-02,
     &     0.000E+00,7.000E-06,2.560E+06,5.620E-03,0.000E+00,1.760E-04,
     &     8.160E+05,2.840E-03,6.210E+06,1.480E-02,6.350E+05,1.340E-03,
     &     1.520E+06,3.170E-03,3.130E+07,1.230E-02,2.110E+09,4.020E-01,
     &     1.780E+09,2.990E-01,1.770E+10,3.770E+00,4.460E+08,3.620E-02/
      DATA(geraddec(j,15),gecolstr(15,j),j=16,nlvl)/0.000E+00,1.900E-04,
     &     2.830E+04,1.500E-03,0.000E+00,8.120E-04,0.000E+00,7.920E-04,
     &     0.000E+00,8.000E-04,0.000E+00,3.470E-04,0.000E+00,6.110E-04,
     &     2.730E+07,3.750E-02,0.000E+00,1.970E-03,0.000E+00,2.680E-03,
     &     0.000E+00,4.020E-03,7.020E+09,8.120E-01/
      DATA(geraddec(j,16),gecolstr(16,j),j=17,nlvl)/5.070E+01,1.070E-02,
     &     0.000E+00,3.790E-03,0.000E+00,9.400E-03,2.630E-02,1.650E-02,
     &     7.830E-03,3.270E-03,0.000E+00,3.100E-03,4.550E+03,4.290E-03,
     &     2.160E+01,3.910E-03,1.370E+01,3.250E-04,0.000E+00,5.400E-05,
     &     6.450E+04,1.930E-03/
      DATA(geraddec(j,17),gecolstr(17,j),j=18,nlvl)/0.000E+00,1.490E-02,
     &     5.710E-04,2.350E-02,2.520E+02,4.160E-02,6.960E-01,1.410E-02,
     &     1.530E-01,2.140E-02,9.590E+03,1.100E-02,1.210E+04,7.940E-03,
     &     1.990E+04,1.920E-03,3.200E+01,3.640E-03,1.260E+05,5.100E-03/
      DATA(geraddec(j,18),gecolstr(18,j),j=19,nlvl)/4.730E-06,9.030E-02,
     &     0.000E+00,3.520E-02,2.380E-04,4.200E-02,3.050E+01,1.160E-01,
     &     0.000E+00,1.980E-02,5.820E-01,3.650E-04,9.000E+01,2.450E-02,
     &     2.130E+05,3.810E-02,0.000E+00,2.970E-03/
      DATA(geraddec(j,19),gecolstr(19,j),j=20,nlvl)/4.350E-03,4.290E-02,
     &     9.040E+01,8.460E-02,1.420E+02,4.280E-02,1.770E+00,1.630E-02,
     &     2.040E+05,3.220E-02,1.090E+03,9.970E-04,2.300E+02,7.880E-04,
     &     2.200E+02,1.480E-02/
      DATA(geraddec(j,20),gecolstr(20,j),j=21,nlvl)/7.410E+00,2.880E-02,
     &     2.910E+02,8.020E-02,5.100E+03,1.190E-02,2.630E+04,7.930E-03,
     &     7.680E+04,8.080E-03,1.560E+04,9.070E-03,5.470E+04,6.390E-03/
      DATA(geraddec(j,21),gecolstr(21,j),j=22,nlvl)/1.830E+01,4.680E-02,
     &     3.260E+03,4.780E-02,4.360E+04,1.250E-02,4.420E+04,8.220E-03,
     &     1.050E+04,7.180E-03,8.420E+04,5.590E-03/
      DATA(geraddec(j,22),gecolstr(22,j),j=23,nlvl)/4.490E-03,1.120E-02,
     &     9.260E+03,5.090E-04,8.190E+04,2.550E-02,7.740E+04,2.660E-02,
     &     5.360E-01,1.380E-03/
      DATA(geraddec(j,23),gecolstr(23,j),j=24,nlvl)/1.840E+03,7.180E-03,
     &     3.960E+04,1.710E-02,4.260E+00,1.340E-02,2.160E+04,3.180E-03/
      DATA(geraddec(j,24),gecolstr(24,j),j=25,nlvl)/1.640E+00,2.940E-02,
     &     1.770E+02,8.200E-02,4.470E+02,3.760E-02/
      DATA(geraddec(j,25),gecolstr(25,j),j=26,nlvl)/1.980E+00,7.310E-02,
     &     4.460E+03,3.690E-02/
      DATA geraddec(nlvl,26),gecolstr(26,nlvl)/7.630E-02,1.770E-02/

      END

      BLOCK DATA kr
      PARAMETER(nlvl=27,ntrs=35)
      CHARACTER*9 clevel
      COMMON /cratom/ crraddec(nlvl,nlvl),crcolstr(nlvl,nlvl),
     &                crenergy(nlvl),keykr(nlvl),jtotkr(nlvl),
     &                klowkr(ntrs),kupkr(ntrs)
      COMMON /crlevel/ clevel(nlvl)
c level configurations and energies (in cm**-1)

      DATA(keykr(i),clevel(i),jtotkr(i),crenergy(i),i=1,nlvl)/1,
     &     '2p6    1S',0,0,2,'2p5 3s 3P',2,13367620,3,'2p5 3s 1P',1,
     &     13392937,4,'2p5 3p 3P',1,13783244,5,'2p5 3s 3P',0,13798566,6,
     &     '2p5 3p 1D',2,13803571,7,'2p5 3s 3P',1,13812078,8,
     &     '2p5 3p 3D',3,13895333,9,'2p5 3p 3S',1,13899903,10,
     &     '2p5 3p 3P',2,13931512,11,'2p5 3p 3P',0,14059763,12,
     &     '2p5 3p 1P',1,14223809,13,'2p5 3p 3D',1,14335770,14,
     &     '2p5 3p 3D',2,14343874,15,'2p5 3p 1S',0,14411330,16,
     &     '2p5 3d 3P',0,14413897,17,'2p5 3d 3P',1,14434401,18,
     &     '2p5 3d 1F',3,14463043,19,'2p5 3d 3P',2,14470024,20,
     &     '2p5 3d 3F',4,14472024,21,'2p5 3d 3D',2,14494151,22,
     &     '2p5 3d 3D',3,14520915,23,'2p5 3d 3D',1,14611481,24,
     &     '2p5 3d 3F',2,14887987,25,'2p5 3d 1D',2,14918502,26,
     &     '2p5 3d 3F',3,14929785,27,'2p5 3d 1P',1,14997775/

c most important transitions
c     data klow/9*1,3,6,4,3,6,12,6,4,8,7,9,10,14,13,8,2,3,10,5,7,2,
      DATA klowkr/9*1,3,4,3,6,12,6,4,7,9,10,14,13,8,2,3,5,7,2,7,3,2,5,2,
     &     7,3,3/
c     data kup/27,23,17,14,10,7,6,3,2,15,21,19,11,19,24,18,16,22,15,21,
c    +         22,26,25,20,10,10,19,13,14,8,13,9,6,12,4,12,6,5/
      DATA kupkr/27,23,17,14,10,7,6,3,2,15,19,11,19,24,18,16,15,21,22,
     &     26,25,20,10,10,13,14,8,13,9,6,12,4,12,6,5/
c radiative decay /raddec/ rates aji (/s) and dimensionless electron
c collision strength /colstr/ at an energy of 110 rydbergs.

      DATA(crraddec(j,1),crcolstr(1,j),j=2,nlvl)/8.000E+08,7.010E-04,
     &     3.910E+12,8.050E-04,2.280E+06,1.390E-03,0.000E+00,1.400E-04,
     &     5.380E+09,1.590E-03,3.520E+12,7.050E-04,0.000E+00,2.100E-03,
     &     2.890E+05,9.420E-04,5.950E+09,1.280E-03,0.000E+00,5.280E-03,
     &     1.180E+03,8.660E-04,3.560E+06,9.100E-04,8.530E+09,1.540E-03,
     &     0.000E+00,1.700E-02,0.000E+00,1.190E-03,1.790E+11,3.210E-03,
     &     0.000E+00,2.670E-03,0.000E+00,3.550E-03,0.000E+00,4.280E-03,
     &     0.000E+00,2.250E-03,0.000E+00,1.900E-03,6.550E+13,2.860E-02,
     &     0.000E+00,2.130E-03,0.000E+00,2.130E-03,0.000E+00,2.400E-03,
     &     1.030E+14,4.040E-02/
      DATA(crraddec(j,2),crcolstr(2,j),j=3,nlvl)/1.290E+02,1.810E-03,
     &     8.740E+09,1.460E+00,4.890E+02,5.420E-03,5.010E+09,1.200E+00,
     &     1.230E+06,1.800E-02,1.800E+10,3.200E+00,4.100E+08,3.360E-02,
     &     1.110E+10,1.160E+00,0.000E+00,9.950E-04,2.820E+07,4.830E-04,
     &     5.580E+08,6.210E-03,6.580E+07,1.270E-03,0.000E+00,2.320E-04,
     &     8.900E+05,9.070E-03,9.590E+05,2.660E-02,4.790E+05,3.130E-02,
     &     1.120E+06,4.350E-02,1.140E+06,7.970E-02,9.350E+03,5.750E-03,
     &     7.850E+05,3.740E-02,2.240E+04,3.510E-03,1.050E+03,3.300E-05,
     &     3.730E+04,2.680E-04,1.270E+04,1.300E-04,1.540E+04,3.170E-04/
      DATA(crraddec(j,3),crcolstr(3,j),j=4,nlvl)/1.440E+04,3.210E-02,
     &     1.350E+06,4.290E-03,4.470E+09,1.290E+00,2.320E+05,9.580E-03,
     &     0.000E+00,4.880E-03,1.530E+10,1.340E+00,9.360E+09,1.130E+00,
     &     3.310E+10,3.100E-01,1.140E+07,9.600E-05,4.140E+05,5.200E-05,
     &     4.740E+07,8.940E-04,2.080E+10,6.180E-02,2.050E+03,6.890E-04,
     &     8.860E+03,2.290E-03,5.510E+05,3.510E-02,8.160E+03,3.640E-03,
     &     0.000E+00,6.040E-03,1.090E+06,4.060E-02,5.550E+05,2.800E-02,
     &     1.850E+06,2.420E-02,1.470E+03,4.000E-05,7.120E+02,4.500E-05,
     &     2.470E+03,6.400E-05,3.700E+05,1.200E-03/
      DATA(crraddec(j,4),crcolstr(4,j),j=5,nlvl)/3.760E+03,7.270E-04,
     &     2.580E+01,1.810E-03,9.810E+03,8.800E-04,1.000E+01,7.410E-02,
     &     8.210E+03,6.220E-03,1.130E+04,1.580E-02,3.330E+05,2.910E-03,
     &     2.080E+05,8.190E-03,5.530E+04,3.850E-04,7.440E+03,1.450E-02,
     &     1.360E+06,2.530E-03,2.620E+10,3.750E-01,2.100E+10,8.140E-01,
     &     0.000E+00,1.740E-04,1.130E+10,6.140E-01,0.000E+00,8.630E-03,
     &     1.240E+09,6.290E-02,0.000E+00,4.930E-03,2.040E+09,3.830E-02,
     &     3.510E+07,4.570E-04,6.900E+07,7.420E-04,0.000E+00,5.100E-05,
     &     8.380E+07,7.440E-04/
      DATA(crraddec(j,5),crcolstr(5,j),j=6,nlvl)/0.000E+00,2.700E-05,
     &     2.620E+01,6.370E-03,0.000E+00,1.100E-05,6.640E+02,2.200E-05,
     &     0.000E+00,1.200E-05,0.000E+00,7.500E-05,3.350E+09,5.140E-01,
     &     1.240E+10,8.990E-01,0.000E+00,2.320E-03,0.000E+00,4.000E-04,
     &     0.000E+00,3.000E-06,1.140E+01,7.000E-06,0.000E+00,1.100E-05,
     &     5.590E+02,2.740E-04,0.000E+00,3.000E-06,8.400E+00,1.600E-05,
     &     0.000E+00,4.000E-06,3.340E+02,6.900E-05,4.460E+05,1.830E-02,
     &     7.050E+05,2.440E-02,0.000E+00,3.080E-03,1.360E+03,1.290E-03/
      DATA(crraddec(j,6),crcolstr(6,j),j=7,nlvl)/2.180E-01,5.700E-05,
     &     6.740E+03,5.520E-02,3.980E+03,4.650E-02,8.060E+03,6.190E-02,
     &     1.090E+03,2.010E-02,1.140E+06,1.990E-02,3.950E+01,1.180E-04,
     &     5.340E+01,1.700E-04,2.740E+03,2.450E-03,0.000E+00,1.150E-03,
     &     2.970E+09,1.310E-01,2.570E+10,2.230E+00,8.000E+09,4.830E-01,
     &     0.000E+00,1.410E-02,5.700E+09,3.090E-01,8.980E+07,1.640E-02,
     &     1.010E+09,2.300E-02,3.470E+07,5.740E-04,1.810E+08,2.120E-03,
     &     1.450E+07,3.490E-04,2.810E+07,5.200E-04/
      DATA(crraddec(j,7),crcolstr(7,j),j=8,nlvl)/0.000E+00,8.700E-05,
     &     2.240E+05,5.900E-04,3.160E+05,5.310E-03,2.150E+08,6.600E-02,
     &     5.650E+09,9.630E-01,5.930E+09,4.720E-01,1.850E+10,2.320E+00,
     &     2.430E+10,4.110E-01,0.000E+00,7.000E-06,3.330E+02,1.840E-04,
     &     1.800E+01,6.700E-05,2.160E+02,1.940E-04,0.000E+00,5.000E-05,
     &     2.070E+02,1.380E-04,3.300E+02,1.700E-04,7.300E+03,1.080E-03,
     &     5.730E+05,2.870E-02,4.750E+05,2.250E-02,1.200E+06,6.250E-02,
     &     1.610E+06,2.590E-02/
      DATA(crraddec(j,8),crcolstr(8,j),j=9,nlvl)/1.000E-05,2.060E-02,
     &     1.580E-01,8.450E-02,0.000E+00,2.020E-03,1.320E+00,9.700E-05,
     &     1.090E+03,1.450E-02,1.000E+06,2.310E-02,0.000E+00,4.100E-04,
     &     0.000E+00,2.300E-03,0.000E+00,6.400E-03,2.910E+09,4.240E-01,
     &     1.280E+09,1.330E-01,2.080E+10,3.600E+00,1.660E+07,8.300E-03,
     &     4.950E+09,5.180E-01,0.000E+00,4.830E-03,7.810E+05,5.900E-05,
     &     7.600E+07,1.380E-03,8.370E+07,1.900E-03,0.000E+00,4.040E-04/
      DATA(crraddec(j,9),crcolstr(9,j),j=10,nlvl)/1.390E+01,4.770E-02,
     &     8.140E+03,1.490E-03,1.580E+04,9.370E-04,5.870E+05,5.990E-03,
     &     7.660E+04,6.770E-03,9.150E+04,3.550E-04,3.080E+08,9.460E-03,
     &     5.570E+08,4.370E-02,0.000E+00,6.050E-03,5.170E+08,5.590E-02,
     &     0.000E+00,5.770E-03,1.620E+10,1.410E+00,0.000E+00,6.960E-03,
     &     1.570E+10,4.530E-01,3.300E+03,5.500E-05,1.700E+06,1.060E-04,
     &     0.000E+00,7.200E-05,2.790E+09,2.060E-02/
      DATA(crraddec(j,10),crcolstr(10,j),j=11,nlvl)/2.050E+01,7.390E-03,
     &     4.720E+00,1.500E-04,5.430E+05,1.000E-02,3.380E+05,1.450E-02,
     &     3.040E+03,1.510E-03,0.000E+00,8.580E-04,1.890E+09,1.720E-01,
     &     2.990E+08,6.210E-02,5.640E+09,6.820E-01,0.000E+00,1.120E-02,
     &     1.750E+09,1.500E+00,1.800E+10,2.260E+00,9.020E+07,7.340E-03,
     &     1.690E+07,4.040E-04,3.540E+08,6.060E-03,5.160E+06,9.700E-05,
     &     7.770E+07,8.520E-04/
      DATA(crraddec(j,11),crcolstr(11,j),j=12,nlvl)/1.160E+01,3.640E-04,
     &     1.530E+05,3.790E-03,1.600E+02,1.230E-02,0.000E+00,7.700E-05,
     &     0.000E+00,4.800E-04,4.640E+07,1.210E-02,0.000E+00,2.000E-03,
     &     0.000E+00,1.200E-03,0.000E+00,2.180E-03,0.000E+00,8.430E-04,
     &     0.000E+00,1.330E-03,9.620E+09,6.390E-01,0.000E+00,3.730E-04,
     &     0.000E+00,6.170E-04,0.000E+00,6.860E-04,1.560E+09,1.820E-02/
      DATA(crraddec(j,12),crcolstr(12,j),j=13,nlvl)/2.490E+03,5.120E-02,
     &     1.290E+04,5.770E-02,1.640E+01,1.670E-03,8.220E+03,1.500E-05,
     &     2.970E+04,7.800E-05,0.000E+00,9.700E-05,6.380E+05,1.610E-03,
     &     0.000E+00,2.900E-05,6.870E+04,8.000E-05,0.000E+00,1.010E-04,
     &     2.810E+07,6.250E-03,2.610E+10,1.600E+00,1.880E+06,8.750E-03,
     &     0.000E+00,1.290E-02,1.400E+10,3.110E-01/
      DATA(crraddec(j,13),crcolstr(13,j),j=14,nlvl)/2.510E+00,4.140E-02,
     &     3.460E+03,1.760E-03,5.130E+05,6.150E-03,2.450E+05,3.850E-02,
     &     0.000E+00,6.700E-05,3.170E+05,3.550E-03,0.000E+00,9.600E-05,
     &     3.440E+03,5.100E-05,0.000E+00,5.500E-05,8.120E+06,5.650E-03,
     &     3.440E+08,4.660E-02,1.910E+10,1.770E+00,0.000E+00,1.040E-02,
     &     4.710E+09,1.770E-01/
      DATA(crraddec(j,14),crcolstr(14,j),j=15,nlvl)/2.190E+00,2.220E-02,
     &     0.000E+00,4.000E-06,1.430E+05,4.140E-04,6.030E+04,2.360E-04,
     &     3.990E+05,5.630E-03,0.000E+00,9.600E-05,1.850E-05,4.800E-04,
     &     1.770E+05,3.920E-04,1.200E+07,5.630E-03,2.540E+09,3.080E-01,
     &     2.120E+09,2.170E-01,2.180E+10,2.790E+00,5.590E+08,2.920E-02/
      DATA(crraddec(j,15),crcolstr(15,j),j=16,nlvl)/0.000E+00,1.180E-04,
     &     8.430E+02,4.670E-04,0.000E+00,4.980E-04,0.000E+00,4.660E-04,
     &     0.000E+00,6.160E-04,0.000E+00,3.150E-04,0.000E+00,3.660E-04,
     &     9.730E+06,1.820E-02,0.000E+00,2.210E-03,0.000E+00,2.860E-03,
     &     0.000E+00,3.680E-03,1.120E+10,6.030E-01/
      DATA(crraddec(j,16),crcolstr(16,j),j=17,nlvl)/1.240E+02,9.020E-03,
     &     0.000E+00,8.840E-03,6.240E-02,1.150E-02,0.000E+00,4.090E-03,
     &     1.060E-04,3.240E-03,0.000E+00,3.450E-03,1.110E+04,4.200E-03,
     &     1.290E+02,2.940E-03,5.080E+01,1.250E-04,0.000E+00,3.000E-05,
     &     2.900E+05,1.310E-03/
      DATA(crraddec(j,17),crcolstr(17,j),j=18,nlvl)/5.620E-04,2.160E-02,
     &     5.710E+02,3.570E-02,0.000E+00,1.520E-02,1.160E-01,1.150E-02,
     &     3.070E-01,1.850E-02,2.220E+04,1.040E-02,5.990E+04,6.630E-03,
     &     7.970E+04,1.190E-03,1.480E+02,2.160E-03,5.900E+05,3.570E-03/
      DATA(crraddec(j,18),crcolstr(18,j),j=19,nlvl)/1.300E+00,5.000E-02,
     &     0.000E+00,7.750E-02,2.000E+02,5.290E-02,4.140E+02,3.930E-02,
     &     3.350E+00,1.440E-02,9.780E+05,2.450E-02,2.310E+03,4.720E-04,
     &     7.550E+02,3.680E-04,1.210E+03,1.230E-02/
      DATA(crraddec(j,19),crcolstr(19,j),j=20,nlvl)/0.000E+00,3.050E-02,
     &     4.030E+01,2.580E-02,8.580E+02,5.680E-02,1.400E+04,1.300E-02,
     &     2.040E+05,9.460E-03,2.310E+05,3.850E-03,4.660E+04,4.720E-03,
     &     4.040E+05,5.590E-03/
      DATA(crraddec(j,20),crcolstr(20,j),j=21,nlvl)/5.960E-05,4.440E-02,
     &     5.230E+01,9.130E-02,0.000E+00,2.190E-02,1.390E+00,1.870E-04,
     &     5.700E+02,1.830E-02,9.850E+05,2.880E-02,0.000E+00,2.000E-03/
      DATA(crraddec(j,21),crcolstr(21,j),j=22,nlvl)/2.480E+01,5.020E-02,
     &     4.180E+03,3.110E-02,1.350E+05,6.350E-03,3.350E+05,8.220E-03,
     &     7.350E+04,8.180E-03,2.520E+05,3.170E-03/
      DATA(crraddec(j,22),crcolstr(22,j),j=23,nlvl)/4.780E-03,1.120E-02,
     &     2.430E+03,2.520E-04,4.050E+05,8.220E-03,3.760E+05,2.090E-02,
     &     6.920E-01,9.430E-04/
      DATA(crraddec(j,23),crcolstr(23,j),j=24,nlvl)/8.330E+03,3.870E-03,
     &     2.400E+05,1.150E-02,5.610E+01,1.000E-02,8.570E+04,3.400E-02/
      DATA(crraddec(j,24),crcolstr(24,j),j=25,nlvl)/1.380E+01,2.950E-02,
     &     6.880E+02,7.290E-02,5.650E+02,3.320E-02/
      DATA(crraddec(j,25),crcolstr(25,j),j=26,nlvl)/4.780E+00,6.680E-02,
     &     5.520E+03,3.400E-02/
      DATA crraddec(nlvl,26),crcolstr(26,nlvl)/4.260E-02,2.160E-02/
      END


c *********************************************************************
      SUBROUTINE gain(istage,dlamda,rh,ti,nst,nfl,istflg)
      PARAMETER(kk=1001)
c this subroutine calculates the gain between specified levels of
c hydrogen-like ions.
      COMMON /atom  / pu(10,kk),p(10,kk),pz(10,kk),zs(10,kk),zst(kk),
     &                zstz(kk),z(kk),a(kk),fn(10),deg(10),qs(10,kk),
     &                zstmi(kk),zstu(kk),engu(10,kk),eng(10,kk),
     &                rdrhe(kk),zstin(kk),nmax
      COMMON /atomdt/ foss(10,10),sig(10,10)
      COMMON /alph  / rwi(kk,3),frinv(kk,3),alpha(kk,3),ideg(10),
     &                idegm(10)
      COMMON /cvoigt/ bcfl(10,kk),va(kk)
      COMMON /npp   / npp(10,10),il,iu
      DIMENSION rh(kk),ti(kk)
      DIMENSION npl(10),npu(10)
      DIMENSION denne(4),swidth(4)
      DATA denne/5.E19,1.E20,2.5E20,5.E20/
      DATA swidth/1.975E-4,2.37E-4,3.95E-4,7.11E-4/

c setup populations
      nl=il-(istage-1)
      nu=iu-(istage-1)
      DO 100 i=1,nmax
       npl(i)=npp(i,nl)
       npu(i)=npp(i,nu)
 100  CONTINUE
c calculate the first term in gain formula
      pi=3.1415927
      spl=2.9979E10
      an=6.022045E23

      wa=(dlamda**2)/(8.0*pi)
      za=dlamda**3/(8.*pi*spl)
c calculate the a value
      engd=1.24E-07/dlamda
      wb=4.4E+13*foss(il,iu)*engd*engd
      degcor=deg(il)/deg(iu)
      wb=wb*degcor
c ad1 neon
c     if (istage.eq.2) then
c      if (iu.eq.3.and.il.eq.2) then
c        correction for term splittingin case of 3d5/2 - 2p3/2.
c        in li-like  al,ne  (2.419/8/.637) see grasp
c       wb=.4706*wb
c      elseif (iu.eq.4.and.il.eq.3) then
c        correction for term splittingin case of 4f7/2 - 3d5/2.
c        in li-like al and ne (5.8/10/.8408) see grasp
c       wb=.3836*wb
c      endif
c     endif
c     zb=wb
c loop over the different cells in the problem
c calculate the doppler width
      DO 200 l=nst,nfl
       tionev=1000.0*ti(l)
       dnion=(rh(l)*an)/a(l)
       dne=zst(l)*dnion

       dopw=1.11E+13*engd*sqrt(tionev/a(l))
c calculate the voigt a parameter
       va(l)=(bcfl(il,l)+bcfl(iu,l))/(4.0*pi*dopw)
       wc=1.0/(sqrt(pi)*dopw)
c ad1 neon 1  / (delta lamda/lamda)
c      zc=4.e-11*(iu*iu-il*il)*dne**(2./3.)*zst(l)**(.3333333)
c      zc=z(l)/(zc*dlamda)

c ad1 neon case
c      if (il.eq.3.and.iu.eq.4) then
c       zc=1./3.e-4
c      elseif (il.eq.2.and.iu.eq.3) then
c ad1 neon case
c       call linint(denne,swidth,dne,rstark,4)
c       zc=1./rstark
c      endif
c      rwi(l,1)=1./zc

c calculate the ion number density and inversion
       dnu=frac(npu,l)*dnion
       dnl=frac(npl,l)*dnion
       wd=dnu-((deg(iu)/deg(il))*dnl)
       frinv(l,1)=0.
       if(dnu.gt.0.)frinv(l,1)=wd/dnu
       zd=wd
c calculate gain
       alpha(l,1)=wa*wb*wc*wd
c        calculate stark broadening and corresponding small signal gain
c      alpha(l,1)=za*zb*zc*zd
c simple correction for fine-structure on balmer alpha in h-like scheme
       IF(istage.eq.1)THEN
        IF(il.eq.2.and.iu.eq.3)alpha(l,1)=0.5*alpha(l,1)
       ENDIF
c insert stark broadening correction factor for h-like balmer alpha
       IF(istflg.ne.0.and.istage.eq.1)THEN
        IF(il.eq.2.and.iu.eq.3)THEN
         tiev=1000.0*ti(l)
         stw=1.2019*spl*stark(dne,tiev)
         corr=1.0+((stw/dopw)**2)
         corr=1.0/sqrt(corr)
         alpha(l,1)=corr*alpha(l,1)
        ENDIF
       ENDIF
 200  CONTINUE
      RETURN
      END

c**********************************************************************
c**********************************************************************
      SUBROUTINE line(nst,nfl,rh,ti,rb,fl,nlmax,time,npt,ny,xdta,ydta,
     &                fwide,istflg,istage,ngeom)
      PARAMETER(kk=1001,underf=1.E-30)
      COMMON /atom  / pu(10,kk),p(10,kk),pz(10,kk),zs(10,kk),zst(kk),
     &                zstz(kk),z(kk),a(kk),fn(10),deg(10),qs(10,kk),
     &                zstmi(kk),zstu(kk),engu(10,kk),eng(10,kk),
     &                rdrhe(kk),zstin(kk),nmax
      COMMON /atomdt/ foss(10,10),sig(10,10)
      COMMON /axint / axint(10,6,2),xt(10,6)
      COMMON /grap  / ipt
      COMMON /npp   / npp(10,10),il,iu
      DIMENSION rh(kk),ti(kk),rb(kk)
      DIMENSION fl(2),sum(10,2)
      DIMENSION tx(501)
      DIMENSION xdta(4000),ydta(6,4000)
      DIMENSION npl(10),npu(10)
c this subroutine calculates the axial photon flux in alpha, beta
c and gamma lines
      pi=3.1415927
      an=6.0232E+23
      const=2.654E-02
      spl=2.9979E10
      iupt=iu+2
c ad1  convert fwide to radians (cylindrical geometry only)
      IF(ngeom.gt.1)THEN
c     for cyl geometry fwide in radians
       awide=fwide/rb(nfl+1)
       IF(awide.le.0.0.or.awide.gt.3.141592655)awide=6.283185307
      ELSE
       awide=fwide
       IF(awide.le.0.0)awide=1.
      ENDIF
c zero arrays
      DO 100 k=1,nlmax
       DO 50 j=1,nmax
        DO 20 m=1,2
         axint(j,k,m)=0.0
 20     CONTINUE
        xt(j,k)=0.0
        sum(j,k)=0.0
 50    CONTINUE
 100  CONTINUE
c loop over the different fibre lengths
      DO 300 k=1,nlmax
c perform normal and line-narrowed calculations
       DO 150 m=1,2
c loop over the different cells
        DO 140 l=nst,nfl
c setup population of the lower level
         nl=il-(istage-1)
         DO 110 ii=1,nmax
          npl(ii)=npp(ii,nl)
 110     CONTINUE
         frl=frac(npl,l)
         gl=deg(il)
         dnion=(rh(l)*an)/a(l)
c loop over the different transitions
         DO 120 j=iu,iupt
          IF(istage.eq.1)ed=0.0136*z(l)
     &                      **2*(1.0/real(il*il)-1.0/real(j*j))
          IF(istage.eq.2)ed=0.0136*(z(l)-2.)
     &                      **2*(1.0/real(il*il)-1.0/real(j*j))
          IF(istage.eq.3)ed=0.0136*(z(l)-10.)
     &                      **2*(1.0/real(il*il)-1.0/real(j*j))
          xl=1.24E-07/ed
          veli=1.38E+06*sqrt((1000.0*ti(l))/a(l))
          phi0=xl/veli
          phi0=phi0/sqrt(pi)
c make the stark broadening correction only for h-like balmer alpha
          IF(istflg.ne.0.and.istage.eq.1)THEN
           IF(j.eq.3)THEN
            dne=zst(l)*dnion
            tiev=1000.0*ti(l)
            dopw=1.11E+13*ed*sqrt(tiev/a(l))
            stw=1.2019*spl*stark(dne,tiev)
            corr=1.0+((stw/dopw)**2)
            corr=1.0/sqrt(corr)
            phi0=phi0*corr
           ENDIF
          ENDIF
c set up populations for the upper level
          nj=j-(istage-1)
          DO 115 ii=1,nmax
           npu(ii)=npp(ii,nj)
 115      CONTINUE
          fru=frac(npu,l)
          IF(abs(fru).gt.underf)THEN
           gu=deg(j)
           wa=3.147E+31*(ed**2)
           wb=((frl*gu)/(fru*gl))-1.0
           wc=const*foss(il,j)*(frl-((gl/gu)*fru))*dnion*phi0*fl(k)
c simple correction for fine-structure on h-like balmer alpha
           IF(istage.eq.1.and.j.eq.3)wc=0.5*wc
           dxp=0.0
           IF(wc.le.170.0)dxp=exp(-wc)
           wd=(wa/wb)*(1.0-dxp)
           IF(ngeom.eq.2)THEN
            sum(j,k)=sum(j,k)+0.5*awide*((rb(l+1)**2)-(rb(l)**2))*wd
           ELSEIF(ngeom.eq.1)THEN
            sum(j,k)=sum(j,k)+awide*(rb(l+1)-rb(l))*wd
           ENDIF
           we=(veli*ed)/spl
c make stark broadening correction only for h-like balmer alpha
           IF(istflg.ne.0.and.istage.eq.1)THEN
            IF(j.eq.3)we=we/corr
           ENDIF
c perform integration over the line profile if m=2
c zero array
c number of intervals must be exactly divisable by the number of
c half-widths
           IF(m.eq.1)THEN
c perform normal line width calculation if m=1
            wd=wd*we*sqrt(pi)
           ELSE
            ip=500
            nhw=10
            IF(mod(nhw,2).ne.0)GOTO 500
            IF(mod(ip,2).ne.0)GOTO 500
            IF(mod(ip,nhw).ne.0)GOTO 500
            itp=ip+1
            DO 116 i=1,itp
             tx(i)=0.0
 116        CONTINUE
            step=we/real(ip/nhw)
            ihp=ip/2
            icp=ihp+1
            tx(icp)=wd
            denp=ed
            denm=ed
            DO 118 i=1,ihp
             denp=denp+step
             denm=denm-step
c fill array for integration
             wap=3.147E+31*(denp**2)
             wam=3.147E+31*(denm**2)
             fexp=exp(-(((real(i)*step)/we)**2))
             wc=const*foss(il,j)*(frl-((gl/gu)*fru))
     &          *dnion*phi0*fexp*fl(k)
c simple correction for fine-structure on h-like balmer alpha
             IF(istage.eq.1.and.j.eq.3)wc=0.5*wc
             dxp=0.0
             IF(wc.le.170.0)dxp=exp(-wc)
             tx(icp+i)=(wap/wb)*(1.0-dxp)
             tx(icp-i)=(wam/wb)*(1.0-dxp)
 118        CONTINUE
            wd=sint(itp,step,tx)
           ENDIF
c accumulate area factor
           IF(ngeom.eq.2)THEN
            axint(j,k,m)=axint(j,k,m)
     &                   +0.5*awide*((rb(l+1)**2)-(rb(l)**2))*wd
           ELSEIF(ngeom.eq.1)THEN
            axint(j,k,m)=axint(j,k,m)+awide*(rb(l+1)-rb(l))*wd
           ENDIF
          ENDIF
 120     CONTINUE
 140    CONTINUE
 150   CONTINUE
c calculate the effective temperature in each line
       DO 200 j=iu,iupt
        IF(ngeom.eq.2)THEN
         xa=0.5*awide*((rb(nfl+1)**2)-(rb(nst)**2))
        ELSEIF(ngeom.eq.1)THEN
         xa=awide*(rb(nfl+1)-rb(nst))
        ENDIF
        IF(istage.eq.1)xb=0.0136*z(l)
     &                    **2*(1.0/real(il*il)-1.0/real(j*j))
        IF(istage.eq.2)xb=0.0136*(z(l)-2.)
     &                    **2*(1.0/real(il*il)-1.0/real(j*j))
        IF(istage.eq.3)xb=0.0136*(z(l)-10.)
     &                    **2*(1.0/real(il*il)-1.0/real(j*j))
        xc=3.147E+31*(xb**2)
        xt(j,k)=0.0
        IF(sum(j,k).gt.0.0)THEN
         xd=((xa*xc)/sum(j,k))+1.0
         IF(xd.gt.1.0)THEN
          xt(j,k)=xb/log(xd)
c assuming an angular spread of 10**-4 steradians calculate
c saturation factor for each line
          domega=1.0E-04
          xe=domega/(4.0*pi)
          xf=exp(xb/xt(j,k))-1.0
          xt(j,k)=xe/xf
         ENDIF
        ENDIF
 200   CONTINUE
 300  CONTINUE
c fill next element in array for graphics
c draw only the line-narrowed case
      m=2
      iy=0
      ipt=ipt+1
      xdta(ipt)=time
      DO 400 k=1,nlmax
       DO 350 j=iu,iupt
        iy=iy+1
        ydta(iy,ipt)=axint(j,k,m)
 350   CONTINUE
 400  CONTINUE
      npt=ipt
      ny=iy
 500  RETURN
      END
      FUNCTION sint(n,h,tx)
c this function performs a simpson's rule integration
c n is the number of grid points in the integration
c h is the step size
c tx is the value of the integrand at grid point i
      DIMENSION tx(501)
      nn=n-2
      wa=0.0
      hdt=h/3.0
      DO 100 i=1,nn,2
       wa=wa+(hdt*(tx(i)+(4.0*tx(i+1))+tx(i+2)))
 100  CONTINUE
      sint=wa
      RETURN
      END

c**********************************************************************
      SUBROUTINE refr(dlamda,rh,rb,nst,nfl,fl,nlmax)
      PARAMETER(kk=1001,underf=1.E-30)
c this subroutine calculates the radius of curvature due to refraction
c from the free electrons in the plasma and calculates the maximum
c transverse dimension which can thereby be tolerated.
      COMMON /atom  / pu(10,kk),p(10,kk),pz(10,kk),zs(10,kk),zst(kk),
     &                zstz(kk),z(kk),a(kk),fn(10),deg(10),qs(10,kk),
     &                zstmi(kk),zstu(kk),engu(10,kk),eng(10,kk),
     &                rdrhe(kk),zstin(kk),nmax
      COMMON /refrac/ refrac(2,kk)
      DIMENSION rh(kk),rb(kk),rc(kk),fl(2)
      DIMENSION rind(kk),drind(kk)
      pi=3.141592654
      spl=2.9979E10
      an=6.0232E+23
      nstp=nst+1
      nflm=nfl-1
      nflm2=nfl-2
      DO 100 l=nst,nfl
       rc(l)=0.5*(rb(l)+rb(l+1))
 100  CONTINUE
      DO 200 l=nst,nfl
       dne=(zst(l)*rh(l)*an)/a(l)
       w2=((2.0*pi*spl)/dlamda)**2
       wp2=3.1826E+09*dne
       rind(l)=1.0-(wp2/w2)
       IF(rind(l).lt.0.0)rind(l)=0.0
       rind(l)=sqrt(rind(l))
 200  CONTINUE
      DO 300 l=nst,nflm2
       drind(l+1)=((rind(l+1)-rind(l))/(rc(l+1)-rc(l)))
     &            +((rind(l+2)-rind(l+1))/(rc(l+2)-rc(l+1)))
       drind(l+1)=0.5*drind(l+1)
 300  CONTINUE
      DO 400 l=nstp,nflm
       radcur=drind(l)/rind(l)
       IF(abs(radcur).gt.underf)radcur=1.0/radcur
       DO 350 k=1,nlmax
        refrac(k,l)=0.0
        IF(abs(radcur).gt.underf)refrac(k,l)=fl(k)**2/(8.0*radcur)
 350   CONTINUE
 400  CONTINUE
      RETURN
      END
c**********************************************************************
      FUNCTION hg(x)
c  calculates the holstein escape   factor
      IF(x.le.0.0125)THEN
       hg=1.0-x/1.414
       RETURN
      ELSEIF(x.lt.3.0)THEN
       g1=1.0
       y=1.0
       a=0.0
 50    prevg=g1
       a=a+1.0
       y=-y*x/a
       g1=g1+y/sqrt(a+1.0)
       prevg=abs((prevg-g1)/g1)
       IF(prevg.gt.1.0E-2)GOTO 50
       hg=g1
       RETURN
      ENDIF
      hg=1.0/(x*sqrt(3.142*log(x)))
      RETURN
      END
C********************************************************************
C     FUNCTION HG(X)
C THIS FUNCTION CALCULATES THE escape   FACTOR ACCORDING TO HUMMER
C AND RYBICKI
C     HG=1.0
c     pi=3.141592654
C     IF(X.LT.1.0E-5) GOTO 1
C     WA=EXP(-X)
C     WA=WA*((2.0*X)-1.0)
C     WB=ERFN(sqrt(X))
C     WC=2.0*sqrt(PI)*(X**1.5)
C     WD=1.0/(3.0*X)
C     WE=1.0+WA+WC*(WB-1.0)
C     HG=WD*WE
C   1 RETURN
C     END
C**********************************************************************
c     function erfn(xx)
c     dimension a(6)
c this subroutine calculates the error function erf(x)
c it employs the rational approximation in 'hastings,
c approximations for digital computers',1955.
c set up coefficients
c     data a/0.0705230784e 00,0.0422820123e 00,0.0092705272e 00,
c    @0.0001520143e 00,0.0002765672e 00,0.0000430638e 00/
c     x=abs(xx)
c     erfn=1.0e 00
c     if(x.ge.4.0e 00) goto 1
c     cor=1.0e 00
c     do 2 i=1,6
c     cor=cor+(a(i)*(x**i))
c  2  continue
c     cor=1.0e 00/(cor**16)
c     cor=1.0e 00-cor
c     erfn=cor
c  1  if(xx.lt.0.0e 00)erfn=-erfn
c     return
c     end
c**********************************************************************
      FUNCTION frac(np,l)
      PARAMETER(kk=1001)
c this function calculates the fraction of ions in a given
c configuration
      COMMON /atom  / pu(10,kk),p(10,kk),pz(10,kk),zs(10,kk),zst(kk),
     &                zstz(kk),z(kk),a(kk),fn(10),deg(10),qs(10,kk),
     &                zstmi(kk),zstu(kk),engu(10,kk),eng(10,kk),
     &                rdrhe(kk),zstin(kk),nmax
      COMMON /alph  / rwi(kk,3),frinv(kk,3),alpha(kk,3),ideg(10),
     &                idegm(10)
      DIMENSION np(10)
      prob=1.0
      DO 100 i=1,nmax
       ni=ideg(i)-np(i)
C calculate combinatorial factor
       comfac=factr(ideg(i),ni)
       comfac=comfac/factorial(np(i))
C calculate probabilities factor
       pi=p(i,l)/deg(i)
       pie=1.0
       IF(np(i).ne.0)pie=pi**np(i)
       qi=1.0-pi
       qie=1.0
       IF(ni.ne.0)qie=qi**ni
       prob=prob*comfac*pie*qie
 100  CONTINUE
      frac=prob
      RETURN
      END
c**********************************************************************
      FUNCTION factorial(n)
c this subroutine caculates n! and returns the value as a real
      factorial=1.0
      IF(n.gt.0)THEN
       DO 50 i=1,n
        factorial=factorial*real(i)
 50    CONTINUE
      ENDIF
      RETURN
      END
c**********************************************************************
      FUNCTION factr(nu,nl)
c this subroutine caculates nu!/nl! (where nu ge nl)
c and returns the value as a real
      nlp=nl+1
      factr=1.0
      IF(nu.gt.nl)THEN
       factr=nu
       IF(nu.ne.nlp)THEN
        factr=1.0
        DO 20 i=nlp,nu
         factr=factr*real(i)
 20     CONTINUE
       ENDIF
      ENDIF
      RETURN
      END
c**********************************************************************
      FUNCTION stark(xa,ya)
c this function calculates the value of the motional stark broadening
c half width in cm**-1 by interpolation in values calculated by oza,
c greene and kelleher (1986)
c parameters      : xa     : the electron density (in cm**-3)
c                 : ya     : the temperature (in e.v.)
      REAL x,y,z,hwhh1(31),hwhh2(31),hwhh3(31),eden(31),t(3),yfit(3),xa,
     &     ya
      DATA t/1.30103,2.0,2.4771213/
      DATA eden/17.00,17.10,17.20,17.30,17.40,17.50,17.60,17.70,17.80,
     &     17.90,18.00,18.10,18.20,18.30,18.40,18.50,18.60,18.70,18.80,
     &     18.90,19.00,19.10,19.20,19.30,19.40,19.50,19.60,19.70,19.80,
     &     19.90,20.00/
      DATA hwhh1/0.83,0.87,0.90,0.93,0.96,0.99,1.01,1.04,1.08,1.10,1.14,
     &     1.17,1.20,1.24,1.28,1.32,1.35,1.39,1.44,1.48,1.53,1.58,1.63,
     &     1.69,1.74,1.81,1.88,1.95,2.02,2.09,2.19/
      DATA hwhh2/0.79,0.87,0.94,0.99,1.06,1.11,1.18,1.24,1.29,1.34,1.39,
     &     1.44,1.49,1.54,1.59,1.63,1.68,1.71,1.75,1.79,1.84,1.89,1.93,
     &     1.97,2.00,2.04,2.08,2.12,2.16,2.20,2.24/
      DATA hwhh3/0.60,0.70,0.80,0.90,0.99,1.08,1.18,1.26,1.34,1.43,1.51,
     &     1.59,1.68,1.74,1.81,1.88,1.94,1.99,2.05,2.10,2.15,2.19,2.23,
     &     2.26,2.29,2.33,2.35,2.38,2.39,2.40,2.41/
c the data is held in log form in the look up tables to decrease
c errors, so convert input quantities to logs
      x=log10(xa)
      y=log10(ya)
c interpolate value on each of the three known temperature lines
c for the given electron density
      CALL interp(eden,hwhh1,x,yinterp,31)
      yfit(1)=yinterp
      CALL interp(eden,hwhh2,x,yinterp,31)
      yfit(2)=yinterp
      CALL interp(eden,hwhh3,x,yinterp,31)
      yfit(3)=yinterp
c now interpolate linearly across these three temperatures for the
c specified temperature y
      CALL interp(t,yfit,y,z,3)
c convert interpolated value from log to actual
      stark=10.0**z
c introduce minimum in the interpolation set by static stark broadening
      starkm=4.8E-17*(xa**0.9232)
      IF(stark.lt.starkm)stark=starkm
      RETURN
      END

c**********************************************************************
      SUBROUTINE effg(rilong,rishort,rlong,rshort,alpha)
c     this calculates the effective gain as seen from axial
c     and transverse viewpoints i.e. the experimental observable.
      alpha=0.
      IF((rilong.lt.1.E-06).or.(rishort.lt.1.E-06))THEN
       WRITE(6,*)' itensities are zero'
       RETURN
      ENDIF
      r=rilong/rishort
      alphmax=170.0/rlong
      alph0=0.0001
      DO 100 i1=1,999
       IF(alph0.gt.alphmax)alph0=alphmax
       IF(alph0*rlong.lt.-170.)THEN
        el=0.0
       ELSE
        el=exp(alph0*rlong)
       ENDIF
       IF(alph0*rshort.lt.-170.)THEN
        es=0.0
       ELSE
        es=exp(alph0*rshort)
       ENDIF
       x1=(el-1.0)/(es-1.0)
       alph1=alph0-(((x1-r)*(es-1.0))/((rlong*el)-x1*rshort*es))
       diff=abs(alph1/alph0)
       IF((diff.le.1.001).and.(diff.gt.0.999))GOTO 200
       alph0=alph1
 100  CONTINUE
      WRITE(6,*)'not converging in effective alpha'
      WRITE(6,*)'error (%): ',abs(1.0-diff)*100.0
 200  alpha=alph1
      RETURN
      END
c**********************************************************************

      BLOCK DATA bdat
      COMMON /atomdt/ foss(10,10),sig(10,10)
      DATA sig/0.6250,0.9380,0.9810,0.9870,0.9940,0.9970,0.9990,1.0000,
     &     1.0000,1.0000,0.2350,0.6900,0.8930,0.9400,0.9700,0.9840,
     &     0.9900,0.9930,0.9950,1.0000,0.1090,0.3970,0.7020,0.8500,
     &     0.9200,0.9550,0.9700,0.9800,0.9900,1.0000,0.0617,0.2350,
     &     0.4780,0.7050,0.8300,0.9000,0.9500,0.9700,0.9800,0.9900,
     &     0.0398,0.1550,0.3310,0.5310,0.7200,0.8300,0.9000,0.9500,
     &     0.9700,0.9800,0.0277,0.1090,0.2390,0.4000,0.4800,0.7350,
     &     0.8300,0.9000,0.9500,0.9700,0.0204,0.0808,0.1780,0.3100,
     &     0.4590,0.6100,0.7450,0.8300,0.9000,0.9500,0.0156,0.0625,
     &     0.1380,0.2430,0.3710,0.5060,0.6350,0.7500,0.8300,0.9000,
     &     0.0123,0.0494,0.1110,0.1940,0.2990,0.4310,0.5440,0.6560,
     &     0.7000,0.8300,0.0100,0.0400,0.0900,0.1580,0.2450,0.3530,
     &     0.4600,0.5760,0.6700,0.7650/
c screening constants from jqsrt vol 43, n 2, pp. 149-154
c     data sig/   0.5966,0.8597,0.9923,0.9800,0.9725,
c    1    0.9970,0.9990,0.9999,0.9999,0.9999,
c    1   0.2345,0.6888,0.8877,0.9640,1.0000,
c    1    0.9880,0.9900,0.9990,0.9999,0.9999,
c    1   0.1093,0.4018,0.7322,0.9415,0.9897,
c    1    0.9820,0.9860,0.9900,0.9920,0.9999,
c    1   0.0622,0.2430,0.5150,0.6986,0.8590,
c    1    0.9600,0.9750,0.9830,0.9860,0.9900,
c    1   0.0399,0.1597,0.3527,0.5888,0.8502,
c    1    0.8300,0.9000,0.9500,0.9700,0.9800,
c    1   0.0277,0.1098,0.2455,0.4267,0.5764,
c    1    0.7248,0.8300,0.9000,0.9500,0.9700,
c    1   0.0204,0.0808,0.1811,0.3184,0.4592,
c    1    0.6098,0.7374,0.8300,0.9000,0.9500,
c    1   0.0156,0.0624,0.1392,0.2457,0.3711,
c    1    0.5062,0.6355,0.7441,0.8300,0.9000,
c    1   0.0123,0.0493,0.1102,0.1948,0.2994,
c    1    0.4222,0.5444,0.6558,0.7553,0.8300,
c    1   0.0100,0.0400,0.0900,0.1584,0.2450,
c    1    0.3492,0.4655,0.5760,0.6723,0.7612/
      END
c**********************************************************************
      BLOCK DATA bdnp
c this block data sets up populations of ground and excited
c configurations of h-like, he-like and li-like states
      COMMON /npist / npp1(10,10),npp2(10,10),npp3(10,10),npp4(10,10),
     &                npp5(10,10)
c h-like
      DATA npp1/1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,
     &     0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,
     &     0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,
     &     1,0,0,0,0,0,0,0,0,0,0,1/
c li-like
      DATA npp2/2,1,0,0,0,0,0,0,0,0,2,0,1,0,0,0,0,0,0,0,2,0,0,1,0,0,0,0,
     &     0,0,2,0,0,0,1,0,0,0,0,0,2,0,0,0,0,1,0,0,0,0,2,0,0,0,0,0,1,0,
     &     0,0,2,0,0,0,0,0,0,1,0,0,2,0,0,0,0,0,0,0,1,0,2,0,0,0,0,0,0,0,
     &     0,1,0,0,0,0,0,0,0,0,0,0/
c na-like
      DATA npp3/2,8,1,0,0,0,0,0,0,0,2,8,0,1,0,0,0,0,0,0,2,8,0,0,1,0,0,0,
     &     0,0,2,8,0,0,0,1,0,0,0,0,2,8,0,0,0,0,1,0,0,0,2,8,0,0,0,0,0,1,
     &     0,0,2,8,0,0,0,0,0,0,1,0,2,8,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,
     &     0,0,0,0,0,0,0,0,0,0,0,0/
c ne-like
      DATA npp4/2,8,0,0,0,0,0,0,0,0,2,7,1,0,0,0,0,0,0,0,2,7,0,1,0,0,0,0,
     &     0,0,2,7,0,0,1,0,0,0,0,0,2,7,0,0,0,1,0,0,0,0,2,7,0,0,0,0,1,0,
     &     0,0,2,7,0,0,0,0,0,1,0,0,2,7,0,0,0,0,0,0,1,0,2,7,0,0,0,0,0,0,
     &     0,1,0,0,0,0,0,0,0,0,0,0/
c ni-like
      DATA npp5/2,8,18,0,0,0,0,0,0,0,2,8,17,1,0,0,0,0,0,0,2,8,17,0,1,0,
     &     0,0,0,0,2,8,17,0,0,1,0,0,0,0,2,8,17,0,0,0,1,0,0,0,2,8,17,0,0,
     &     0,0,1,0,0,2,8,17,0,0,0,0,0,1,0,2,8,17,0,0,0,0,0,0,1,0,0,0,0,
     &     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/
      END
c********************************************************************
      BLOCK DATA bdcf
      COMMON /ioncfg/ iconfig(10,0:50)
      DATA(iconfig(i,0),i=1,10)/0,0,0,0,0,0,0,0,0,0/
      DATA(iconfig(i,1),i=1,10)/1,0,0,0,0,0,0,0,0,0/
      DATA(iconfig(i,2),i=1,10)/0,1,0,0,0,0,0,0,0,0/
      DATA(iconfig(i,3),i=1,10)/0,0,1,0,0,0,0,0,0,0/
      DATA(iconfig(i,4),i=1,10)/0,0,0,1,0,0,0,0,0,0/
      DATA(iconfig(i,5),i=1,10)/0,0,0,0,1,0,0,0,0,0/
      DATA(iconfig(i,6),i=1,10)/0,0,0,0,0,1,0,0,0,0/
      DATA(iconfig(i,7),i=1,10)/0,0,0,0,0,0,1,0,0,0/
      DATA(iconfig(i,8),i=1,10)/0,0,0,0,0,0,0,1,0,0/
      DATA(iconfig(i,9),i=1,10)/0,0,0,0,0,0,0,0,1,0/
      DATA(iconfig(i,10),i=1,10)/0,0,0,0,0,0,0,0,0,1/
      DATA(iconfig(i,11),i=1,10)/2,0,0,0,1,0,0,0,0,0/
c      DATA(iconfig(i,2),i=1,10)/2,0,0,0,0,0,0,0,0,0/
c      DATA(iconfig(i,3),i=1,10)/2,1,0,0,0,0,0,0,0,0/
c      DATA(iconfig(i,4),i=1,10)/2,2,0,0,0,0,0,0,0,0/
c      DATA(iconfig(i,5),i=1,10)/2,3,0,0,0,0,0,0,0,0/
c      DATA(iconfig(i,6),i=1,10)/2,4,0,0,0,0,0,0,0,0/
c      DATA(iconfig(i,7),i=1,10)/2,5,0,0,0,0,0,0,0,0/
c      DATA(iconfig(i,8),i=1,10)/2,6,0,0,0,0,0,0,0,0/
c      DATA(iconfig(i,9),i=1,10)/2,7,0,0,0,0,0,0,0,0/
c      DATA(iconfig(i,10),i=1,10)/2,8,0,0,0,0,0,0,0,0/
c      DATA(iconfig(i,11),i=1,10)/2,8,1,0,0,0,0,0,0,0/
      DATA(iconfig(i,12),i=1,10)/2,8,2,0,0,0,0,0,0,0/
      DATA(iconfig(i,13),i=1,10)/2,8,3,0,0,0,0,0,0,0/
      DATA(iconfig(i,14),i=1,10)/2,8,4,0,0,0,0,0,0,0/
      DATA(iconfig(i,15),i=1,10)/2,8,5,0,0,0,0,0,0,0/
      DATA(iconfig(i,16),i=1,10)/2,8,6,0,0,0,0,0,0,0/
      DATA(iconfig(i,17),i=1,10)/2,8,7,0,0,0,0,0,0,0/
      DATA(iconfig(i,18),i=1,10)/2,8,8,0,0,0,0,0,0,0/
      DATA(iconfig(i,19),i=1,10)/2,8,9,0,0,0,0,0,0,0/
      DATA(iconfig(i,20),i=1,10)/2,8,10,0,0,0,0,0,0,0/
      DATA(iconfig(i,21),i=1,10)/2,8,11,0,0,0,0,0,0,0/
      DATA(iconfig(i,22),i=1,10)/2,8,12,0,0,0,0,0,0,0/
      DATA(iconfig(i,23),i=1,10)/2,8,13,0,0,0,0,0,0,0/
      DATA(iconfig(i,24),i=1,10)/2,8,14,0,0,0,0,0,0,0/
      DATA(iconfig(i,25),i=1,10)/2,8,15,0,0,0,0,0,0,0/
      DATA(iconfig(i,26),i=1,10)/2,8,16,0,0,0,0,0,0,0/
      DATA(iconfig(i,27),i=1,10)/2,8,17,0,0,0,0,0,0,0/
      DATA(iconfig(i,28),i=1,10)/2,8,18,0,0,0,0,0,0,0/
      DATA(iconfig(i,29),i=1,10)/2,8,18,1,0,0,0,0,0,0/
      DATA(iconfig(i,30),i=1,10)/2,8,18,2,0,0,0,0,0,0/
c     DATA(iconfig(i,29),i=1,10)/2,8,16,1,0,0,0,0,0,0/
c     DATA(iconfig(i,30),i=1,10)/2,8,17,1,0,0,0,0,0,0/

      DATA(iconfig(i,31),i=1,10)/2,8,18,3,0,0,0,0,0,0/
      DATA(iconfig(i,32),i=1,10)/2,8,18,4,0,0,0,0,0,0/
      DATA(iconfig(i,33),i=1,10)/2,8,18,5,0,0,0,0,0,0/
      DATA(iconfig(i,34),i=1,10)/2,8,18,6,0,0,0,0,0,0/
      DATA(iconfig(i,35),i=1,10)/2,8,18,7,0,0,0,0,0,0/
      DATA(iconfig(i,36),i=1,10)/2,8,18,8,0,0,0,0,0,0/
      DATA(iconfig(i,37),i=1,10)/2,8,18,9,0,0,0,0,0,0/
      DATA(iconfig(i,38),i=1,10)/2,8,18,10,0,0,0,0,0,0/
      DATA(iconfig(i,39),i=1,10)/2,8,18,11,0,0,0,0,0,0/
      DATA(iconfig(i,40),i=1,10)/2,8,18,12,0,0,0,0,0,0/
      DATA(iconfig(i,41),i=1,10)/2,8,18,13,0,0,0,0,0,0/
      DATA(iconfig(i,42),i=1,10)/2,8,18,14,0,0,0,0,0,0/
      DATA(iconfig(i,43),i=1,10)/2,8,18,15,0,0,0,0,0,0/
      DATA(iconfig(i,44),i=1,10)/2,8,18,16,0,0,0,0,0,0/
      DATA(iconfig(i,45),i=1,10)/2,8,18,17,0,0,0,0,0,0/
      DATA(iconfig(i,46),i=1,10)/2,8,18,18,0,0,0,0,0,0/
      DATA(iconfig(i,47),i=1,10)/2,8,18,19,0,0,0,0,0,0/
      DATA(iconfig(i,48),i=1,10)/2,8,18,20,0,0,0,0,0,0/
      DATA(iconfig(i,49),i=1,10)/2,8,18,21,0,0,0,0,0,0/
      DATA(iconfig(i,50),i=1,10)/2,8,18,22,0,0,0,0,0,0/
      END
c**********************************************************************
      BLOCK DATA bdatn
      CHARACTER*2 aname(50)
      CHARACTER*10 shname(12)
      COMMON /atname/ aname
      COMMON /sname / shname
      DATA aname/' H','He','Li','Be',' B',' C',' N',' O',' F','Ne','Na',
     &     'Mg','Al','Si',' P',' S','Cl',' A',' K','Ca','Sc','Ti','Va',
     &     'Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br',
     &     'Kr','Rb','Sr',' Y','Zr','Cb','Mo','Ma','Ru','Rh','Pd','Ag',
     &     'Cd','In','Sn'/
      DATA shname/'   K-shell','   L-shell','   M-shell','   N-shell',
     &     '   O-shell','   P-shell','   Q-shell','   R-shell',
     &     '   S-shell','   T-shell','   =======','     Z*   '/
      END

c**********************************************************************
      SUBROUTINE loadcg(np,ionstage)
      INTEGER ionstage,np(10)
      COMMON /ioncfg/ iconfig(10,0:50)
      DO 100 i=1,10
       np(i)=iconfig(i,ionstage)
 100  CONTINUE
      RETURN
      END


      SUBROUTINE raray2(kname,pa,kdimx,kx,ky)
c u.20 print doubly-subscripted array
c c1.1 basic system parameters
      LOGICAL nlend,nlres
      COMMON /combas/ altime,cptime,stime,ndiary,nin,nledge,nonlin,nout,
     &                nprint,npunch,nread,nrec,nresm,nrun,nstep,nlend,
     &                nlres
c end of basic system parameters
      CHARACTER kname*(*)
      DIMENSION pa(*)
c output array name
      CALL blines(1)
      WRITE(nout,10100)kname
      CALL blines(1)
c output array
      DO 100 jy=1,ky
c index of first element of line
       i1=1+(jy-1)+kdimx
c index of last element
       i2=i1+kx-1
c output line
       WRITE(nout,10200)(pa(j),j=i1,i2)
 100  CONTINUE
      RETURN
10100 FORMAT(4x,a)
10200 FORMAT((6x,10(1pe12.3)/))
      END


      SUBROUTINE rarray(kname,pa,kdim)
c u.8 print name and values of real array
c c1.1 basic system parameters
      LOGICAL nlend,nlres
      COMMON /combas/ altime,cptime,stime,ndiary,nin,nledge,nonlin,nout,
     &                nprint,npunch,nread,nrec,nresm,nrun,nstep,nlend,
     &                nlres
c end of basic system parameters
      CHARACTER kname*(*)
      DIMENSION pa(kdim)
      CALL blines(1)
      WRITE(nout,10100)kname
      CALL blines(1)
      WRITE(nout,10200)pa
      RETURN
10100 FORMAT(4x,a)
10200 FORMAT((6x,10(1pe12.3)))
      END


      SUBROUTINE report(kclass,ksub,kpoint)
c 5.1 control the diagnostics
      CALL dumcom(kclass,ksub,kpoint)
      RETURN
      END


      SUBROUTINE repthd(kclass,ksub,kpoint)
c u.11 print heading for diagnostics report
c c1.1 basic system parameters
      LOGICAL nlend,nlres
      COMMON /combas/ altime,cptime,stime,ndiary,nin,nledge,nonlin,nout,
     &                nprint,npunch,nread,nrec,nresm,nrun,nstep,nlend,
     &                nlres
c end of basic system parameters
      WRITE(nout,10100)kclass,ksub,kpoint
      RETURN
10100 FORMAT(4x,'class = ',i6,', subprogram = ',i6,', point = ',i6)
      END


      SUBROUTINE resetc(ka,kvalue)
c u.16 reset character string to specified value
      CHARACTER ka*(*),kvalue*(*)
      i=len(ka)
      ka=kvalue(1:i)
      RETURN
      END


      SUBROUTINE reseti(ka,kdim,kvalue)
c u.15 reset integer array to specified value
      DIMENSION ka(kdim)
      DO 100 j=1,kdim
       ka(j)=kvalue
 100  CONTINUE
      RETURN
      END


      SUBROUTINE resetl(kla,kdim,klval)
c u.19 reset logical array to specified value
      LOGICAL kla(kdim),klval
      DO 100 j=1,kdim
       kla(j)=klval
 100  CONTINUE
      RETURN
      END


      SUBROUTINE resetr(pa,kdim,pvalue)
c u.14 reset real array to specified value
      DIMENSION pa(kdim)
      DO 100 j=1,kdim
       pa(j)=pvalue
 100  CONTINUE
      RETURN
      END


      SUBROUTINE resumerun
      PARAMETER(kk=1001)
      PARAMETER(ntrs=35)
c  reads  restart data file for continuation
c c1.1 basic system parameters
      LOGICAL nlend,nlres
      COMMON /combas/ altime,cptime,stime,ndiary,nin,nledge,nonlin,nout,
     &                nprint,npunch,nread,nrec,nresm,nrun,nstep,nlend,
     &                nlres
c end of basic system parameters
c c5.1 input-output control variables
      LOGICAL nlfilm,nlhcpy,nlprnt
      DIMENSION buf1(kk),buf2(kk),buf3(kk),buf4(kk),buf5(kk),buf6(kk)
      REAL nprnt,nproc
      COMMON /comout/ nfilm,nhdcpy,np1,np2,np3,nlfilm,nlhcpy,nlprnt,
     &        ztnext,nprnt,nproc,buf1,buf2,buf3,buf4,buf5,buf6,scp,
     &                scr,scrho,scte,scti,sctime,nt1,nt2
c end of input-output control variables
c c1.9 development and diagnostic parameters
c end of development and diagnostic parameters
c c2.1 hydrodynamics
      DIMENSION dm(kk),p3(kk),r1(kk),r3(kk),r5(kk),rho1(kk),rho3(kk),
     &          rho5(kk),u2(kk),u4(kk),v1(kk),v3(kk),v5(kk),lstat(kk)
      COMMON /comhyd/ dm,p3,pini,r1,r3,r5,rho1,rho3,rho5,rhoini,rhor,
     &                rini,time,u2,u4,uedge,v1,v3,v5
c end of hydrodynamics
c c2.2 thermodynamics
      DIMENSION ddroe1(kk),ddroe3(kk),ddroi1(kk),ddroi3(kk),ddte1(kk),
     &          ddte3(kk),ddti1(kk),ddti3(kk),xappae(kk),xappai(kk),
     &          pe1(kk),pe3(kk),pi1(kk),pi3(kk),te1(kk),te3(kk),ti1(kk),
     &          ti3(kk)
      COMMON /comth / ddroe1,ddroe3,ddroi1,ddroi3,ddte1,ddte3,ddti1,
     &                ddti3,gammae,gammai,xappae,xappai,pe1,pe3,pi1,pi3,
     &                te1,te3,teini,ti1,ti3,tiini
c end of thermodynamics
c c2.3 ions and electrons
      DIMENSION brems1(kk),brems3(kk),ddz1(kk),ddz3(kk),degen(kk),
     &          effz(kk),eixch2(kk),fidash(kk),fz1(kk),fz3(kk),fzsq1(kk)
     &          ,fzsq3(kk),xlc(kk),xmieff(kk),xne(kk),xni(kk),omega1(kk)
      COMMON /comie / brems1,brems3,ddz1,ddz3,degen,degmax,degmin,effz,
     &                eixch2,fidash,fz1,fz3,fzsq1,fzsq3,xlc,xmieff,xne,
     &                xni,omega1,pmass
c end of ions and electrons
c c2.4 laser variables
      DIMENSION alpha1(kk),xaser1(kk),xl1(kk),xl3(kk)
      COMMON /comlas/ alpha1,elas1,xamda1,xaser1,xecri1,plas1,rabs1,
     &                rocri1,xl1,xl3,nabs1
c end of laser variables
c c2.5 thermonuclear reactions
      DIMENSION f1d(kk),f1h(kk),f1he3(kk),f1he4(kk),f1neu(kk),f1ntrl(kk)
     &          ,f1t(kk),f1x(kk),f1x1(kk),f1x2(kk),f1x3(kk),f1x4(kk),
     &          f3d(kk),f3h(kk),f3he3(kk),f3he4(kk),f3neu(kk),f3ntrl(kk)
     &          ,f3t(kk),f3x(kk),f3x1(kk),f3x2(kk),f3x3(kk),f3x4(kk),
     &          r1dd(kk),r1dhe3(kk),r1dt(kk),r3dd(kk),r3dhe3(kk),
     &          r3dt(kk),ye1(kk),ye3(kk),yi1(kk),yi3(kk)
      COMMON /comfus/ deuter,f1d,f1h,f1he3,f1he4,f1neu,f1ntrl,f1t,f1x,
     &                f1x1,f1x2,f1x3,f1x4,f3d,f3h,f3he3,f3he4,f3neu,
     &                f3ntrl,f3t,f3x,f3x1,f3x2,f3x3,f3x4,heliu3,heliu4,
     &                hydrog,xetral,xtrlms,pneut1,pneut3,r1dd,r1dhe3,
     &                r1dt,r3dd,r3dhe3,r3dt,rneut1,rneut3,tinucl,totneu,
     &                tritiu,xmass,xmass1,xmass2,xmass3,xmass4,xtra,
     &                xtra1,xtra2,xtra3,xtra4,xz,xz1,xz2,xz3,xz4,ye1,
     &                ye3,yi1,yi3,yield
c end of thermonuclear reactions
c c2.6 - energies
      COMMON /comen / eerror,en,einput,eloss,eefuse,eifuse,eindt1,
     &                eindt3,eneutr,pv,usqm,xionen,en0,eions,ecbbs,
     &                ecbfs,edcvs,erbfs,eatis,ezsts,erbbs,erfbs,edfbs
c c2.6 - end energies
c c2.7 physics control
      LOGICAL mlbrms,mlecon,mlfuse,mlicon,mlx,nlabs,nlbrms,nlburn,
     &        nlcri1,nldepo,nlecon,nlfuse,nlicon,nlmove,nlpfe,nlpfi,nlx,
     &        nlhf,nltnl
      COMMON /comcon/ mstep,ncase,ngeom,mlbrms,mlecon,mlfuse,mlicon,mlx,
     &                nlabs,nlbrms,nlburn,nlcri1,nldepo,nlecon,nlfuse,
     &                nlicon,nlmove,nlpfe,nlpfi,nlx,nlhf,nltnl,pmin,
     &                rhomin,temin,timin,umin
c end of physics control
c c3.1 numerical control parameters
      LOGICAL break,nlgoon,nlite
      COMMON /comnc / nceldt,ncondt,nit,nitmax,break,nlgoon,nlite,ak0,
     &                ak1,ak2,ak3,ak4,ak5,ak6,ak7,ak8,ak9,bneum,deltat,
     &                dt2,dt3,dt4,dtemax,dtimax,dtmax,dtmin,dumax,rdt2,
     &                rdt3,rdt4
c end of numerical control parameters
c c3.2 mesh and numerical methods
      DIMENSION ae(kk),ai(kk),be(kk),bi(kk),ce(kk),ci(kk),de(kk),di(kk),
     &          e(kk),f(kk),ge(kk),gi(kk),q2(kk),q4(kk),teite(kk),
     &          tiite(kk),uite(kk)
      COMMON /comnum/ mesh,nj,njm1,njp1,nl,nlm1,nlp1,ae,ai,be,bi,ce,ci,
     &                de,di,e,f,ge,gi,q2,q4,teite,tiite,uite
c end of mesh and numerical methods
c c4.1 - administrative variables
c c4.1 - end administrative variables
      COMMON /atom  / pu(10,kk),p(10,kk),pz(10,kk),zs(10,kk),zst(kk),
     &                zstz(kk),z(kk),a(kk),fn(10),deg(10),qs(10,kk),
     &                zstmi(kk),zstu(kk),engu(10,kk),eng(10,kk),
     &                rdrhe(kk),zstin(kk),nmax
      COMMON /dpi   / ecbb(kk),ecbf(kk),eraddr(kk),erbf(kk),eati(kk),
     &                dcvz(kk),efbx2(kk),efrx2(kk),erbb(kk),erfb(kk),
     &                acfc(10,kk),acfr(10,kk),bcfc(10,kk),bcfr(10,kk),
     &                edfb(kk),ter(kk),dqdtim(kk),nmaxr
      COMMON /alph  / rwi(kk,3),frinv(kk,3),alpha(kk,3),ideg(10),
     &                idegm(10)
      COMMON /atomdt/ foss(10,10),sig(10,10)
      COMMON /npp   / npp(10,10),il,iu
      COMMON /netint/ taxint(ntrs,6),tgain(ntrs)
c x-ray laser common blocks
      DIMENSION fl(2)
      COMMON /xrl   / ropmul,tpon,tpoff,rmpd,ztrtop,dlamda,fl,fwide,
     &                efbmul,drmul,bbtrap,fbtrap,icxrl,istage,idflg,
     &                itbflg,istflg,idrflg,ipuflg,nlp,nup,llp,lup,nst,
     &                nfl,nlmax,ifrsta,ilosta,ihista,igstat
      COMMON /radflg/ yamda,xhnu,xihnu,dihnu,xplas1(kk),xrbf(kk),
     &                atrbf(kk),corhf1(kk),corhf(kk),iradf1,iradf2
c ad1 line transfer
      PARAMETER(nangs=6,nfreqs=200,nlvls=3,nzons=301)
      logical lproces
      COMMON /detect/ photons(nangs,nfreqs),solid,area
     & ,prnext,prlast,proces,lproces
      LOGICAL inuse
c
c     INQUIRE(unit=nresm,opened=inuse)
c     IF(inuse)THEN
c rewind restart unit
c      REWIND nresm
c     ELSE
c      OPEN(unit=nresm,form='unformatted',status='unknown')
c     ENDIF
      READ(nresm)nstep,mstep,time,uedge,dt2,rdt2,dtmin,dtmax1,dt3
      READ(nresm)rdt3,dt4,rdt4,ztnext,nprnt1
     & ,ncondt,nceldt,plas0,elas1
      READ(nresm)xamda0,nabs0,rabs0,rocri0,xecri0,rhor,pneut1,rneut1
      READ(nresm)eerror,elas1,einput,eloss,eefuse,eifuse,eindt1,eindt3
      READ(nresm)eneutr,pv,usqm,xionen,en,en0,eions,ecbbs,ecbfs,edcvs
      READ(nresm)erbfs,eatis,ezsts,erbbs,erfbs,edfbs,ncase,mlbrms
      READ(nresm)mlecon,mlfuse,mlicon,mlx,pmin,rhomin,temin,timin
      READ(nresm)umin,mesh,nj,njm1,njp1,nl,nlm1,nlp1,np1,np2,np3,nst
      READ(nresm)nfl,nmax,nmaxr,il,iu
      READ(nresm)nang,nfreq
      if(dt2.gt.dtmax) then
c     changes input data
       dt2=dtmax
       rdt2=1./dt2
       dt3=dt2
       rdt3=rdt2
      endif
      if(nprnt.ne.nprnt1) then
      IF(nprnt1.lt.0.)THEN
        ztlast=time+nprnt1*1.0E-12
        if(nprnt.lt.0.)ztnext=ztlast-nprnt*1.0E-12
      ENDIF
      endif


c 2. physics
c 2.1 hydrodynamics
      READ(nresm)(r1(j),j=1,njp1)
      READ(nresm)(r3(j),j=1,njp1)
      READ(nresm)(v1(l),l=1,nlp1)
      READ(nresm)(v3(l),l=1,nlp1)
      READ(nresm)(dm(l),l=1,nlp1)
      READ(nresm)(rho1(l),l=1,nlp1)
      READ(nresm)(rho3(l),l=1,nlp1)
      READ(nresm)(u2(j),j=1,njp1)
      READ(nresm)(lstat(j),j=1,njp1)
c 2.2 thermodynamics
      READ(nresm)(ddroe1(l),l=1,nlp1)
      READ(nresm)(ddroi1(l),l=1,nlp1)
      READ(nresm)(ddte1(l),l=1,nlp1)
      READ(nresm)(ddti1(l),l=1,nlp1)
      READ(nresm)(ti1(l),l=1,nlp1)
      READ(nresm)(te1(l),l=1,nlp1)
      READ(nresm)(pi1(l),l=1,nlp1)
      READ(nresm)(pe1(l),l=1,nlp1)
      READ(nresm)(p3(l),l=1,nlp1)
      READ(nresm)(q2(l),l=1,nlp1)
      READ(nresm)(fz1(l),l=1,nlp1)
      READ(nresm)(fzsq1(l),l=1,nlp1)
      READ(nresm)(ddz1(l),l=1,nlp1)
      READ(nresm)(xappai(l),l=1,nlp1)
      READ(nresm)(xappae(l),l=1,nlp1)
c 2.3 ions and electrons
      READ(nresm)(ae(l),l=1,nlp1)
      READ(nresm)(ai(l),l=1,nlp1)
      READ(nresm)(be(l),l=1,nlp1)
      READ(nresm)(bi(l),l=1,nlp1)
      READ(nresm)(ce(l),l=1,nlp1)
      READ(nresm)(ci(l),l=1,nlp1)
      READ(nresm)(de(l),l=1,nlp1)
      READ(nresm)(di(l),l=1,nlp1)
      READ(nresm)(ge(l),l=1,nlp1)
      READ(nresm)(gi(l),l=1,nlp1)
      READ(nresm)(brems1(l),l=1,nlp1)
      READ(nresm)(xne(l),l=1,nlp1)
      READ(nresm)(xni(l),l=1,nlp1)
      READ(nresm)(xmieff(l),l=1,nlp1)
      READ(nresm)(effz(l),l=1,nlp1)
      READ(nresm)(xlc(l),l=1,nlp1)
c 2.4 laser variables
      READ(nresm)(xl1(l),l=1,nlp1)
c 2.5      thermonuclear reactions
      READ(nresm)(f1d(l),l=1,nlp1)
      READ(nresm)(f1h(l),l=1,nlp1)
      READ(nresm)(f1he3(l),l=1,nlp1)
      READ(nresm)(f1he4(l),l=1,nlp1)
      READ(nresm)(f1neu(l),l=1,nlp1)
      READ(nresm)(f1ntrl(l),l=1,nlp1)
      READ(nresm)(f1t(l),l=1,nlp1)
      READ(nresm)(f1x(l),l=1,nlp1)
      READ(nresm)(f1x1(l),l=1,nlp1)
      READ(nresm)(f1x2(l),l=1,nlp1)
      READ(nresm)(f1x3(l),l=1,nlp1)
      READ(nresm)(f1x4(l),l=1,nlp1)
      READ(nresm)(r1dd(l),l=1,nlp1)
      READ(nresm)(r1dhe3(l),l=1,nlp1)
      READ(nresm)(r1dt(l),l=1,nlp1)
      READ(nresm)(yi1(l),l=1,nlp1)
      READ(nresm)(ye1(l),l=1,nlp1)
c hot electrons
c
c ionisation
      READ(nresm)(zst(l),l=1,nlp1)
      READ(nresm)(zstz(l),l=1,nlp1)
      READ(nresm)(z(l),l=1,nlp1)
      READ(nresm)(a(l),l=1,nlp1)
      READ(nresm)(zstmi(l),l=1,nlp1)
      READ(nresm)(ecbb(l),l=1,nlp1)
      READ(nresm)(ecbf(l),l=1,nlp1)
      READ(nresm)(eraddr(l),l=1,nlp1)
      READ(nresm)(erbf(l),l=1,nlp1)
      READ(nresm)(eati(l),l=1,nlp1)
      READ(nresm)(corhf(l),l=1,nlp1)
      READ(nresm)(corhf1(l),l=1,nlp1)
      READ(nresm)(dcvz(l),l=1,nlp1)
      READ(nresm)(efbx2(l),l=1,nlp1)
      READ(nresm)(efrx2(l),l=1,nlp1)
      READ(nresm)(erbb(l),l=1,nlp1)
      READ(nresm)(erfb(l),l=1,nlp1)
      READ(nresm)(edfb(l),l=1,nlp1)
      READ(nresm)fn,deg,ideg,idegm,npp,foss,sig
      DO 100 i=1,nmax
       READ(nresm)(p(i,l),l=1,nlp1)
       READ(nresm)(eng(i,l),l=1,nlp1)
       READ(nresm)(zs(i,l),l=1,nlp1)
       READ(nresm)(qs(i,l),l=1,nlp1)
 100  CONTINUE
      READ(nresm)taxint,tgain
      READ(nresm)prnext,prlast,photons
c
c prepare for next iteration, copy to next level where necessary.
      CALL copyr(ti1,1,ti3,1,nlp1)
      CALL copyr(te1,1,te3,1,nlp1)
      CALL copyr(pe1,1,pe3,1,nlp1)
      CALL copyr(pi1,1,pi3,1,nlp1)
      CALL copyr(u2,1,u4,1,njp1)
      CALL copyr(q2,1,q4,1,nlp1)
      CALL copyr(r3,1,r5,1,njp1)
      CALL copyr(v3,1,v5,1,nlp1)
      CALL copyr(rho3,1,rho5,1,nlp1)
c
      CALL copyr(f1d,1,f3d,1,nlp1)
      CALL copyr(f1h,1,f3h,1,nlp1)
      CALL copyr(f1he3,1,f3he3,1,nlp1)
      CALL copyr(f1he4,1,f3he4,1,nlp1)
      CALL copyr(f1ntrl,1,f3ntrl,1,nlp1)
      CALL copyr(f1t,1,f3t,1,nlp1)
      CALL copyr(f1x,1,f3x,1,nlp1)
      CALL copyr(f1x1,1,f3x1,1,nlp1)
      CALL copyr(f1x2,1,f3x2,1,nlp1)
      CALL copyr(f1x3,1,f3x3,1,nlp1)
      CALL copyr(f1x4,1,f3x4,1,nlp1)
c
      CALL copyr(fz1,1,fz3,1,nlp1)
      CALL copyr(fzsq1,1,fzsq3,1,nlp1)
c ionisation, keep old values
      DO 200 l=1,mesh
       zstu(l)=zst(l)
       DO 150 i=1,nmax
        pu(i,l)=p(i,l)
        engu(i,l)=eng(i,l)
 150   CONTINUE
 200  CONTINUE
c hot electrons if implemented
c
      RETURN
      END


      SUBROUTINE dumpstat
      PARAMETER(kk=1001)
      PARAMETER(ntrs=35)
c  dump state of simulation for use in restart case if needed
c c1.1 basic system parameters
      LOGICAL nlend,nlres
      COMMON /combas/ altime,cptime,stime,ndiary,nin,nledge,nonlin,nout,
     &                nprint,npunch,nread,nrec,nresm,nrun,nstep,nlend,
     &                nlres
c end of basic system parameters
c c5.1 input-output control variables
      LOGICAL nlfilm,nlhcpy,nlprnt
      DIMENSION buf1(kk),buf2(kk),buf3(kk),buf4(kk),buf5(kk),buf6(kk)
      REAL nprnt,nproc
      COMMON /comout/ nfilm,nhdcpy,np1,np2,np3,nlfilm,nlhcpy,nlprnt,
     &        ztnext,nprnt,nproc,buf1,buf2,buf3,buf4,buf5,buf6,scp,
     &                scr,scrho,scte,scti,sctime,nt1,nt2
c end of input-output control variables
c c1.9 development and diagnostic parameters
c end of development and diagnostic parameters
c c2.1 hydrodynamics
      DIMENSION dm(kk),p3(kk),r1(kk),r3(kk),r5(kk),rho1(kk),rho3(kk),
     &          rho5(kk),u2(kk),u4(kk),v1(kk),v3(kk),v5(kk),lstat(kk)
      COMMON /comhyd/ dm,p3,pini,r1,r3,r5,rho1,rho3,rho5,rhoini,rhor,
     &                rini,time,u2,u4,uedge,v1,v3,v5
c end of hydrodynamics
c c2.2 thermodynamics
      DIMENSION ddroe1(kk),ddroe3(kk),ddroi1(kk),ddroi3(kk),ddte1(kk),
     &          ddte3(kk),ddti1(kk),ddti3(kk),xappae(kk),xappai(kk),
     &          pe1(kk),pe3(kk),pi1(kk),pi3(kk),te1(kk),te3(kk),ti1(kk),
     &          ti3(kk)
      COMMON /comth / ddroe1,ddroe3,ddroi1,ddroi3,ddte1,ddte3,ddti1,
     &                ddti3,gammae,gammai,xappae,xappai,pe1,pe3,pi1,pi3,
     &                te1,te3,teini,ti1,ti3,tiini
c end of thermodynamics
c c2.3 ions and electrons
      DIMENSION brems1(kk),brems3(kk),ddz1(kk),ddz3(kk),degen(kk),
     &          effz(kk),eixch2(kk),fidash(kk),fz1(kk),fz3(kk),fzsq1(kk)
     &          ,fzsq3(kk),xlc(kk),xmieff(kk),xne(kk),xni(kk),omega1(kk)
      COMMON /comie / brems1,brems3,ddz1,ddz3,degen,degmax,degmin,effz,
     &                eixch2,fidash,fz1,fz3,fzsq1,fzsq3,xlc,xmieff,xne,
     &                xni,omega1,pmass
c end of ions and electrons
c c2.4 laser variables
      DIMENSION alpha1(kk),xaser1(kk),xl1(kk),xl3(kk)
      COMMON /comlas/ alpha1,elas1,xamda1,xaser1,xecri1,plas1,rabs1,
     &                rocri1,xl1,xl3,nabs1
c end of laser variables
c c2.5 thermonuclear reactions
      DIMENSION f1d(kk),f1h(kk),f1he3(kk),f1he4(kk),f1neu(kk),f1ntrl(kk)
     &          ,f1t(kk),f1x(kk),f1x1(kk),f1x2(kk),f1x3(kk),f1x4(kk),
     &          f3d(kk),f3h(kk),f3he3(kk),f3he4(kk),f3neu(kk),f3ntrl(kk)
     &          ,f3t(kk),f3x(kk),f3x1(kk),f3x2(kk),f3x3(kk),f3x4(kk),
     &          r1dd(kk),r1dhe3(kk),r1dt(kk),r3dd(kk),r3dhe3(kk),
     &          r3dt(kk),ye1(kk),ye3(kk),yi1(kk),yi3(kk)
      COMMON /comfus/ deuter,f1d,f1h,f1he3,f1he4,f1neu,f1ntrl,f1t,f1x,
     &                f1x1,f1x2,f1x3,f1x4,f3d,f3h,f3he3,f3he4,f3neu,
     &                f3ntrl,f3t,f3x,f3x1,f3x2,f3x3,f3x4,heliu3,heliu4,
     &                hydrog,xetral,xtrlms,pneut1,pneut3,r1dd,r1dhe3,
     &                r1dt,r3dd,r3dhe3,r3dt,rneut1,rneut3,tinucl,totneu,
     &                tritiu,xmass,xmass1,xmass2,xmass3,xmass4,xtra,
     &                xtra1,xtra2,xtra3,xtra4,xz,xz1,xz2,xz3,xz4,ye1,
     &                ye3,yi1,yi3,yield
c end of thermonuclear reactions
c c2.6 - energies
      COMMON /comen / eerror,en,einput,eloss,eefuse,eifuse,eindt1,
     &                eindt3,eneutr,pv,usqm,xionen,en0,eions,ecbbs,
     &                ecbfs,edcvs,erbfs,eatis,ezsts,erbbs,erfbs,edfbs
c c2.6 - end energies
c c2.7 physics control
      LOGICAL mlbrms,mlecon,mlfuse,mlicon,mlx,nlabs,nlbrms,nlburn,
     &        nlcri1,nldepo,nlecon,nlfuse,nlicon,nlmove,nlpfe,nlpfi,nlx,
     &        nlhf,nltnl
      COMMON /comcon/ mstep,ncase,ngeom,mlbrms,mlecon,mlfuse,mlicon,mlx,
     &                nlabs,nlbrms,nlburn,nlcri1,nldepo,nlecon,nlfuse,
     &                nlicon,nlmove,nlpfe,nlpfi,nlx,nlhf,nltnl,pmin,
     &                rhomin,temin,timin,umin
c end of physics control
c c3.1 numerical control parameters
      LOGICAL break,nlgoon,nlite
      COMMON /comnc / nceldt,ncondt,nit,nitmax,break,nlgoon,nlite,ak0,
     &                ak1,ak2,ak3,ak4,ak5,ak6,ak7,ak8,ak9,bneum,deltat,
     &                dt2,dt3,dt4,dtemax,dtimax,dtmax,dtmin,dumax,rdt2,
     &                rdt3,rdt4
c end of numerical control parameters
c c3.2 mesh and numerical methods
      DIMENSION ae(kk),ai(kk),be(kk),bi(kk),ce(kk),ci(kk),de(kk),di(kk),
     &          e(kk),f(kk),ge(kk),gi(kk),q2(kk),q4(kk),teite(kk),
     &          tiite(kk),uite(kk)
      COMMON /comnum/ mesh,nj,njm1,njp1,nl,nlm1,nlp1,ae,ai,be,bi,ce,ci,
     &                de,di,e,f,ge,gi,q2,q4,teite,tiite,uite
c end of mesh and numerical methods
c c4.1 - administrative variables
c c4.1 - end administrative variables
      COMMON /atom  / pu(10,kk),p(10,kk),pz(10,kk),zs(10,kk),zst(kk),
     &                zstz(kk),z(kk),a(kk),fn(10),deg(10),qs(10,kk),
     &                zstmi(kk),zstu(kk),engu(10,kk),eng(10,kk),
     &                rdrhe(kk),zstin(kk),nmax
      COMMON /dpi   / ecbb(kk),ecbf(kk),eraddr(kk),erbf(kk),eati(kk),
     &                dcvz(kk),efbx2(kk),efrx2(kk),erbb(kk),erfb(kk),
     &                acfc(10,kk),acfr(10,kk),bcfc(10,kk),bcfr(10,kk),
     &                edfb(kk),ter(kk),dqdtim(kk),nmaxr
      COMMON /alph  / rwi(kk,3),frinv(kk,3),alpha(kk,3),ideg(10),
     &                idegm(10)
      COMMON /atomdt/ foss(10,10),sig(10,10)
      COMMON /npp   / npp(10,10),il,iu
      COMMON /netint/ taxint(ntrs,6),tgain(ntrs)
c x-ray laser common blocks
      DIMENSION fl(2)
      COMMON /xrl   / ropmul,tpon,tpoff,rmpd,ztrtop,dlamda,fl,fwide,
     &                efbmul,drmul,bbtrap,fbtrap,icxrl,istage,idflg,
     &                itbflg,istflg,idrflg,ipuflg,nlp,nup,llp,lup,nst,
     &                nfl,nlmax,ifrsta,ilosta,ihista,igstat
      COMMON /radflg/ yamda,xhnu,xihnu,dihnu,xplas1(kk),xrbf(kk),
     &                atrbf(kk),corhf1(kk),corhf(kk),iradf1,iradf2
c ad1 line transfer
      PARAMETER(nangs=6,nfreqs=200,nlvls=3,nzons=301)
      logical lproces
      COMMON /detect/ photons(nangs,nfreqs),solid,area
     & ,prnext,prlast,proces,lproces
      LOGICAL inuse
      nrst=nresm+1
c      INQUIRE(unit=nrst,opened=inuse)
c      IF(inuse)THEN
c rewind restart unit
c       REWIND nrst
c      ELSE
c       OPEN(unit=nrst,form='unformatted',status='unknown')
c      ENDIF
      WRITE(nrst)nstep,mstep,time,uedge,dt2,rdt2,dtmin,dtmax,dt3
      WRITE(nrst)rdt3,dt4,rdt4,ztnext,nprnt
     & ,ncondt,nceldt,plas1,elas1
      WRITE(nrst)xamda1,nabs1,rabs1,rocri1,xecri1,rhor,pneut1,rneut1
      WRITE(nrst)eerror,elas1,einput,eloss,eefuse,eifuse,eindt1,eindt3
      WRITE(nrst)eneutr,pv,usqm,xionen,en,en0,eions,ecbbs,ecbfs,edcvs
      WRITE(nrst)erbfs,eatis,ezsts,erbbs,erfbs,edfbs,ncase,mlbrms
      WRITE(nrst)mlecon,mlfuse,mlicon,mlx,pmin,rhomin,temin,timin
      WRITE(nrst)umin,mesh,nj,njm1,njp1,nl,nlm1,nlp1,np1,np2,np3,nst
      WRITE(nrst)nfl,nmax,nmaxr,il,iu
      WRITE(nrst)nangs,nfreqs
c 2. physics
c 2.1 hydrodynamics
      WRITE(nrst)(r1(j),j=1,njp1)
      WRITE(nrst)(r3(j),j=1,njp1)
      WRITE(nrst)(v1(l),l=1,nlp1)
      WRITE(nrst)(v3(l),l=1,nlp1)
      WRITE(nrst)(dm(l),l=1,nlp1)
      WRITE(nrst)(rho1(l),l=1,nlp1)
      WRITE(nrst)(rho3(l),l=1,nlp1)
      WRITE(nrst)(u2(j),j=1,njp1)
      WRITE(nrst)(lstat(j),j=1,njp1)
c 2.2 thermodynamics
      WRITE(nrst)(ddroe1(l),l=1,nlp1)
      WRITE(nrst)(ddroi1(l),l=1,nlp1)
      WRITE(nrst)(ddte1(l),l=1,nlp1)
      WRITE(nrst)(ddti1(l),l=1,nlp1)
      WRITE(nrst)(ti1(l),l=1,nlp1)
      WRITE(nrst)(te1(l),l=1,nlp1)
      WRITE(nrst)(pi1(l),l=1,nlp1)
      WRITE(nrst)(pe1(l),l=1,nlp1)
      WRITE(nrst)(p3(l),l=1,nlp1)
      WRITE(nrst)(q2(l),l=1,nlp1)
      WRITE(nrst)(fz1(l),l=1,nlp1)
      WRITE(nrst)(fzsq1(l),l=1,nlp1)
      WRITE(nrst)(ddz1(l),l=1,nlp1)
      WRITE(nrst)(xappai(l),l=1,nlp1)
      WRITE(nrst)(xappae(l),l=1,nlp1)
c 2.3 ions and electrons
      WRITE(nrst)(ae(l),l=1,nlp1)
      WRITE(nrst)(ai(l),l=1,nlp1)
      WRITE(nrst)(be(l),l=1,nlp1)
      WRITE(nrst)(bi(l),l=1,nlp1)
      WRITE(nrst)(ce(l),l=1,nlp1)
      WRITE(nrst)(ci(l),l=1,nlp1)
      WRITE(nrst)(de(l),l=1,nlp1)
      WRITE(nrst)(di(l),l=1,nlp1)
      WRITE(nrst)(ge(l),l=1,nlp1)
      WRITE(nrst)(gi(l),l=1,nlp1)
      WRITE(nrst)(brems1(l),l=1,nlp1)
      WRITE(nrst)(xne(l),l=1,nlp1)
      WRITE(nrst)(xni(l),l=1,nlp1)
      WRITE(nrst)(xmieff(l),l=1,nlp1)
      WRITE(nrst)(effz(l),l=1,nlp1)
      WRITE(nrst)(xlc(l),l=1,nlp1)
c 2.4 laser variables
      WRITE(nrst)(xl1(l),l=1,nlp1)
c 2.5      thermonuclear reactions
      WRITE(nrst)(f1d(l),l=1,nlp1)
      WRITE(nrst)(f1h(l),l=1,nlp1)
      WRITE(nrst)(f1he3(l),l=1,nlp1)
      WRITE(nrst)(f1he4(l),l=1,nlp1)
      WRITE(nrst)(f1neu(l),l=1,nlp1)
      WRITE(nrst)(f1ntrl(l),l=1,nlp1)
      WRITE(nrst)(f1t(l),l=1,nlp1)
      WRITE(nrst)(f1x(l),l=1,nlp1)
      WRITE(nrst)(f1x1(l),l=1,nlp1)
      WRITE(nrst)(f1x2(l),l=1,nlp1)
      WRITE(nrst)(f1x3(l),l=1,nlp1)
      WRITE(nrst)(f1x4(l),l=1,nlp1)
      WRITE(nrst)(r1dd(l),l=1,nlp1)
      WRITE(nrst)(r1dhe3(l),l=1,nlp1)
      WRITE(nrst)(r1dt(l),l=1,nlp1)
      WRITE(nrst)(yi1(l),l=1,nlp1)
      WRITE(nrst)(ye1(l),l=1,nlp1)
c hot electrons
c
c ionisation
      WRITE(nrst)(zst(l),l=1,nlp1)
      WRITE(nrst)(zstz(l),l=1,nlp1)
      WRITE(nrst)(z(l),l=1,nlp1)
      WRITE(nrst)(a(l),l=1,nlp1)
      WRITE(nrst)(zstmi(l),l=1,nlp1)
      WRITE(nrst)(ecbb(l),l=1,nlp1)
      WRITE(nrst)(ecbf(l),l=1,nlp1)
      WRITE(nrst)(eraddr(l),l=1,nlp1)
      WRITE(nrst)(erbf(l),l=1,nlp1)
      WRITE(nrst)(eati(l),l=1,nlp1)
      WRITE(nrst)(corhf(l),l=1,nlp1)
      WRITE(nrst)(corhf1(l),l=1,nlp1)
      WRITE(nrst)(dcvz(l),l=1,nlp1)
      WRITE(nrst)(efbx2(l),l=1,nlp1)
      WRITE(nrst)(efrx2(l),l=1,nlp1)
      WRITE(nrst)(erbb(l),l=1,nlp1)
      WRITE(nrst)(erfb(l),l=1,nlp1)
      WRITE(nrst)(edfb(l),l=1,nlp1)
      WRITE(nrst)fn,deg,ideg,idegm,npp,foss,sig
      DO 100 i=1,nmax
       WRITE(nrst)(p(i,l),l=1,nlp1)
       WRITE(nrst)(eng(i,l),l=1,nlp1)
       WRITE(nrst)(zs(i,l),l=1,nlp1)
       WRITE(nrst)(qs(i,l),l=1,nlp1)
 100  CONTINUE
      WRITE(nrst)taxint,tgain
      WRITE(nrst)prnext,prlast,photons
      RETURN
      END


      SUBROUTINE revers
      PARAMETER(kk=1001)
c 2.23 shifts the levels 1-5 back
c @(#)comhot    1.2 20/5/85 19:47:07 (20/5/85 19:48:45) @(#)
c c1.9 development and diagnostic parameters
      DIMENSION nadump(20),npdump(20),nvdump(20)
      LOGICAL nlched,nlhead(9),nlomt1(50),nlomt2(50),nlomt3(50),nlrept
      COMMON /comddp/ maxdum,mxdump,nadump,nclass,npdump,npoint,nsub,
     &                nvdump,nlched,nlhead,nlomt1,nlomt2,nlomt3,nlrept
c end of development and diagnostic parameters
c c2.1 hydrodynamics
      DIMENSION dm(kk),p3(kk),r1(kk),r3(kk),r5(kk),rho1(kk),rho3(kk),
     &          rho5(kk),u2(kk),u4(kk),v1(kk),v3(kk),v5(kk),lstat(kk)
      COMMON /comhyd/ dm,p3,pini,r1,r3,r5,rho1,rho3,rho5,rhoini,rhor,
     &                rini,time,u2,u4,uedge,v1,v3,v5
c end of hydrodynamics
c c2.2 thermodynamics
      DIMENSION ddroe1(kk),ddroe3(kk),ddroi1(kk),ddroi3(kk),ddte1(kk),
     &          ddte3(kk),ddti1(kk),ddti3(kk),xappae(kk),xappai(kk),
     &          pe1(kk),pe3(kk),pi1(kk),pi3(kk),te1(kk),te3(kk),ti1(kk),
     &          ti3(kk)
      COMMON /comth / ddroe1,ddroe3,ddroi1,ddroi3,ddte1,ddte3,ddti1,
     &                ddti3,gammae,gammai,xappae,xappai,pe1,pe3,pi1,pi3,
     &                te1,te3,teini,ti1,ti3,tiini
c end of thermodynamics
c c2.3 ions and electrons
      DIMENSION brems1(kk),brems3(kk),ddz1(kk),ddz3(kk),degen(kk),
     &          effz(kk),eixch2(kk),fidash(kk),fz1(kk),fz3(kk),fzsq1(kk)
     &          ,fzsq3(kk),xlc(kk),xmieff(kk),xne(kk),xni(kk),omega1(kk)
      COMMON /comie / brems1,brems3,ddz1,ddz3,degen,degmax,degmin,effz,
     &                eixch2,fidash,fz1,fz3,fzsq1,fzsq3,xlc,xmieff,xne,
     &                xni,omega1,pmass
c end of ions and electrons
c c2.4 laser variables
      DIMENSION alpha1(kk),xaser1(kk),xl1(kk),xl3(kk)
      COMMON /comlas/ alpha1,elas1,xamda1,xaser1,xecri1,plas1,rabs1,
     &                rocri1,xl1,xl3,nabs1
c end of laser variables
c c2.5 thermonuclear reactions
      DIMENSION f1d(kk),f1h(kk),f1he3(kk),f1he4(kk),f1neu(kk),f1ntrl(kk)
     &          ,f1t(kk),f1x(kk),f1x1(kk),f1x2(kk),f1x3(kk),f1x4(kk),
     &          f3d(kk),f3h(kk),f3he3(kk),f3he4(kk),f3neu(kk),f3ntrl(kk)
     &          ,f3t(kk),f3x(kk),f3x1(kk),f3x2(kk),f3x3(kk),f3x4(kk),
     &          r1dd(kk),r1dhe3(kk),r1dt(kk),r3dd(kk),r3dhe3(kk),
     &          r3dt(kk),ye1(kk),ye3(kk),yi1(kk),yi3(kk)
      COMMON /comfus/ deuter,f1d,f1h,f1he3,f1he4,f1neu,f1ntrl,f1t,f1x,
     &                f1x1,f1x2,f1x3,f1x4,f3d,f3h,f3he3,f3he4,f3neu,
     &                f3ntrl,f3t,f3x,f3x1,f3x2,f3x3,f3x4,heliu3,heliu4,
     &                hydrog,xetral,xtrlms,pneut1,pneut3,r1dd,r1dhe3,
     &                r1dt,r3dd,r3dhe3,r3dt,rneut1,rneut3,tinucl,totneu,
     &                tritiu,xmass,xmass1,xmass2,xmass3,xmass4,xtra,
     &                xtra1,xtra2,xtra3,xtra4,xz,xz1,xz2,xz3,xz4,ye1,
     &                ye3,yi1,yi3,yield
c end of thermonuclear reactions
c c2.6 - energies
      COMMON /comen / eerror,en,einput,eloss,eefuse,eifuse,eindt1,
     &                eindt3,eneutr,pv,usqm,xionen,en0,eions,ecbbs,
     &                ecbfs,edcvs,erbfs,eatis,ezsts,erbbs,erfbs,edfbs
c c2.6 - end energies
c c2.7 physics control
      LOGICAL mlbrms,mlecon,mlfuse,mlicon,mlx,nlabs,nlbrms,nlburn,
     &        nlcri1,nldepo,nlecon,nlfuse,nlicon,nlmove,nlpfe,nlpfi,nlx,
     &        nlhf,nltnl
      COMMON /comcon/ mstep,ncase,ngeom,mlbrms,mlecon,mlfuse,mlicon,mlx,
     &                nlabs,nlbrms,nlburn,nlcri1,nldepo,nlecon,nlfuse,
     &                nlicon,nlmove,nlpfe,nlpfi,nlx,nlhf,nltnl,pmin,
     &                rhomin,temin,timin,umin
c end of physics control
c c3.1 numerical control parameters
      LOGICAL break,nlgoon,nlite
      COMMON /comnc / nceldt,ncondt,nit,nitmax,break,nlgoon,nlite,ak0,
     &                ak1,ak2,ak3,ak4,ak5,ak6,ak7,ak8,ak9,bneum,deltat,
     &                dt2,dt3,dt4,dtemax,dtimax,dtmax,dtmin,dumax,rdt2,
     &                rdt3,rdt4
c end of numerical control parameters
c c3.2 mesh and numerical methods
      DIMENSION ae(kk),ai(kk),be(kk),bi(kk),ce(kk),ci(kk),de(kk),di(kk),
     &          e(kk),f(kk),ge(kk),gi(kk),q2(kk),q4(kk),teite(kk),
     &          tiite(kk),uite(kk)
      COMMON /comnum/ mesh,nj,njm1,njp1,nl,nlm1,nlp1,ae,ai,be,bi,ce,ci,
     &                de,di,e,f,ge,gi,q2,q4,teite,tiite,uite
c end of mesh and numerical methods
c c4.1 - administrative variables
c c4.1 - end administrative variables
      COMMON /atom  / pu(10,kk),p(10,kk),pz(10,kk),zs(10,kk),zst(kk),
     &                zstz(kk),z(kk),a(kk),fn(10),deg(10),qs(10,kk),
     &                zstmi(kk),zstu(kk),engu(10,kk),eng(10,kk),
     &                rdrhe(kk),zstin(kk),nmax
      COMMON /dpi   / ecbb(kk),ecbf(kk),eraddr(kk),erbf(kk),eati(kk),
     &                dcvz(kk),efbx2(kk),efrx2(kk),erbb(kk),erfb(kk),
     &                acfc(10,kk),acfr(10,kk),bcfc(10,kk),bcfr(10,kk),
     &                edfb(kk),ter(kk),dqdtim(kk),nmaxr
      SAVE nrever
      DATA iclass,isub/2,27/
      DATA nrever/0/

      IF(nlomt2(isub))RETURN
c in this version of revers we terminate the calculation
c 1. record reason for termination
c     call page
      CALL blines(2)
c are basic variables negative or too small
      IF(ncondt.eq.7)THEN
c or is  time centering damaged?
       CALL mesage('      time-centering damaged  ******************')
       CALL blines(2)
       CALL rvar('dt2     ',dt2)
       CALL rvar('dt4     ',dt4)
       CALL endrun
      ELSE
       CALL mesage(
     &  '+++ warning : variable too small , timestep reset to dtmin ???'
     &  )
c    &'+++ variable too small , calculation stoped      ???')
       CALL ivar('variable',ncondt)
       CALL ivar('cell no ',nceldt)
       CALL rvar('time    ',time)
       CALL ivar('iteration',nit)
c         call endrun
      ENDIF
c ad1 reset level 3 to 1 and continue               !!!!!!!!!
c 2.1 hydrodynamics
      CALL copyr(r3,1,r5,1,njp1)
      CALL copyr(v3,1,v5,1,nlp1)
      CALL copyr(rho3,1,rho5,1,nlp1)
      CALL copyr(u2,1,u4,1,njp1)
      CALL copyr(q2,1,q4,1,njp1)
c 2.2 thermodynamics
      CALL copyr(ddroe1,1,ddroe3,1,nlp1)
      CALL copyr(ddroi1,1,ddroi3,1,nlp1)
      CALL copyr(ddte1,1,ddte3,1,nlp1)
      CALL copyr(ddti1,1,ddti3,1,nlp1)
      CALL copyr(ti1,1,ti3,1,nlp1)
      CALL copyr(te1,1,te3,1,nlp1)
      CALL copyr(pi1,1,pi3,1,nlp1)
      CALL copyr(pe1,1,pe3,1,nlp1)
      CALL copyr(fz1,1,fz3,1,nlp1)
      CALL copyr(fzsq1,1,fzsq3,1,nlp1)
c 2.3 ions and electrons
      CALL copyr(brems1,1,brems3,1,nlp1)
c 2.4 laser variables
      IF(ncase.eq.1)CALL copyr(xl1,1,xl3,1,nlp1)
c 2.5      thermonuclear reactions
c is fusion switched on ?
      IF(nlfuse)THEN
c is burning of deuterium and tritium switched on ?
       IF(nlburn)THEN
        CALL copyr(f1d,1,f3d,1,nlp1)
        CALL copyr(f1h,1,f3h,1,nlp1)
        CALL copyr(f1he3,1,f3he3,1,nlp1)
        CALL copyr(f1he4,1,f3he4,1,nlp1)
        CALL copyr(f1neu,1,f3neu,1,nlp1)
        CALL copyr(f1ntrl,1,f3ntrl,1,nlp1)
        CALL copyr(f1t,1,f3t,1,nlp1)
        CALL copyr(f1x,1,f3x,1,nlp1)
        CALL copyr(f1x1,1,f3x1,1,nlp1)
        CALL copyr(f1x2,1,f3x2,1,nlp1)
        CALL copyr(f1x3,1,f3x3,1,nlp1)
        CALL copyr(f1x4,1,f3x4,1,nlp1)
       ENDIF
       pneut3=pneut1
       rneut3=rneut1
       CALL copyr(r1dd,1,r3dd,1,nlp1)
       CALL copyr(r1dhe3,1,r3dhe3,1,nlp1)
       CALL copyr(r1dt,1,r3dt,1,nlp1)
       CALL copyr(yi1,1,yi3,1,nlp1)
       CALL copyr(ye1,1,ye3,1,nlp1)
      ENDIF
c 2.6 energies
      eindt3=eindt1
c
      CALL copyr(u2,1,uite,1,njp1)
      CALL copyr(ti1,1,tiite,1,nlp1)
      CALL copyr(te1,1,teite,1,nlp1)
      CALL resetr(efbx2,nl,0.0)
      CALL resetr(efrx2,nl,0.0)
      CALL resetr(dcvz,nl,0.0)
c ionisation
      DO 100 l=1,mesh
c     shift
       zst(l)=zstu(l)
       DO 50 i=1,nmax
        p(i,l)=pu(i,l)
        eng(i,l)=engu(i,l)
 50    CONTINUE
 100  CONTINUE
      dt4=dtmin
      time=time-dt2+dt4
      dt2=dt4
      dt3=0.5*(dt2+dt4)
      rdt3=1.0/dt3
      rdt4=1.0/dt4
      rdt2=1./dt2
c maximum number of iterations exceeded, end the run
      nrever=nrever+1
      IF(nrever.gt.10)CALL endrun
      RETURN
      END


      SUBROUTINE runtim
c u.12 update cpu time (secs) and print it
c c1.1 basic system parameters
c end of basic system parameters
c cpulft(1)     returns the amount of cpu time (in minutes) used so far.
c     cpu       = second(1)
c     write(nout, 9900) cpu
c9900 format(5x, g10.3, ' ibm cp seconds used so far')
      RETURN
      END

      SUBROUTINE rvar(kname,pvalue)
c u.4 print name and value of real variable
c c1.1 basic system parameters
      LOGICAL nlend,nlres
      COMMON /combas/ altime,cptime,stime,ndiary,nin,nledge,nonlin,nout,
     &                nprint,npunch,nread,nrec,nresm,nrun,nstep,nlend,
     &                nlres
c end of basic system parameters
      CHARACTER kname*(*)
      WRITE(nout,10100)kname,pvalue
      RETURN
10100 FORMAT(4x,a,' = ',1pe12.4)
      END


      SUBROUTINE saha(z,ane,te,effz,nite,zphi,zchi,zm)
      DIMENSION phi(189)
      DATA(phi(i),i=1,80)/0.0,13.6,0.0,24.6,54.4,0.0,5.4,75.6,122.4,0.0,
     &     9.3,18.2,153.8,217.6,0.0,8.3,25.1,37.9,259.3,340.1,0.0,11.2,
     &     24.4,47.9,64.5,392.0,490.0,0.0,14.5,29.6,47.4,77.5,97.9,
     &     552.0,667.0,0.0,13.6,35.1,54.9,77.4,113.9,138.1,739.0,871.0,
     &     0.0,17.4,35.0,62.6,87.1,114.2,157.1,185.0,954.0,1102.0,0.0,
     &     21.6,41.1,63.5,97.0,126.3,157.9,207.0,239.0,1195.0,1360.0,
     &     0.0,5.1,47.3,71.6,98.9,138.4,172.1,208.4,264.0,300.0,1465.0,
     &     1646.0,0.0,7.64,15.0/
      DATA(phi(i),i=81,160)/80.1,109.3,141.2,186.5,224.9,266.0,327.9,
     &     367.0,1761.0,1959.0,0.0,6.0,18.8,28.4,120.0,153.8,190.0,
     &     241.0,284.0,330.0,398.0,442.0,2085.0,2299.0,0.0,8.1,16.3,
     &     33.5,45.1,166.7,205.0,246.0,303.0,351.0,401.0,476.0,523.0,
     &     2436.0,2666.0,0.0,10.5,19.7,30.1,51.3,65.0,220.0,263.0,309.0,
     &     372.0,424.0,479.0,560.0,611.0,2815.0,3061.0,0.0,10.4,23.4,
     &     35.0,47.3,72.5,88.0,281.0,329.0,379.0,447.0,506.0,566.0,
     &     651.0,706.0,3220.0,3482.0,0.0,13.0,23.8,39.9,53.5,67.8,96.7,
     &     114.0/
      DATA(phi(i),i=161,189)/348.0,400.0,455.0,531.0,593.0,663.0,749.0,
     &     807.0,3654.0,3931.0,0.0,15.8,27.6,40.9,59.8,75.0,91.3,124.0,
     &     143.0,422.0,479.0,539.0,621.0,687.0,755.0,854.0,956.0,4115.0,
     &     4407.0/
      zpmass=1.6726E-27
      nz=nint(z)
      IF(nz.gt.18)GOTO 400
      nzp1=nz+1
      ntable=(nz*nzp1/2)-1
      nite=1
      test=log(6.0E27*te*sqrt(te)/ane)
      IF(test.gt.4.0)test=4.0
      IF(test.lt.0.5)test=0.5
      zsum=0.0
      test=test*te
      IF(test.lt.20.0)test=20.0
      DO 100 i=1,nzp1
       itable=ntable+i
       IF(test.lt.phi(itable))GOTO 200
       zsum=zsum+phi(itable)
 100  CONTINUE
      aphi=nzp1
      GOTO 300
 200  p=(phi(itable)-test)/(phi(itable)-phi(itable-1))
      q=1.0-p
      aphi=p*(i-1)+q*i
      zsum=zsum+q*phi(itable)
 300  effz=aphi-1.
      zchi=zsum*1.6E-19/(zm*zpmass)
      zphi=phi(itable)*1.6E-19/(zm*zpmass)
      RETURN
c   for high z materials we assume a z of about 20
 400  effz=20.0
      RETURN
      END


      SUBROUTINE scalei(ka,kdim,kc)
c u.22 scale an integer array by an integer value
      DIMENSION ka(kdim)
      DO 100 j=1,kdim
       ka(j)=ka(j)*kc
 100  CONTINUE
      RETURN
      END


      SUBROUTINE scaler(pa,kdim,pc)
c u.21 scale a real array by a real value
      DIMENSION pa(kdim)
      DO 100 j=1,kdim
       pa(j)=pa(j)*pc
 100  CONTINUE
      RETURN
      END


      SUBROUTINE select
      PARAMETER(kk=1001)
c 3.3 select output and transfer to buffers
c c2.1 hydrodynamics
      DIMENSION dm(kk),p3(kk),r1(kk),r3(kk),r5(kk),rho1(kk),rho3(kk),
     &          rho5(kk),u2(kk),u4(kk),v1(kk),v3(kk),v5(kk),lstat(kk)
      COMMON /comhyd/ dm,p3,pini,r1,r3,r5,rho1,rho3,rho5,rhoini,rhor,
     &                rini,time,u2,u4,uedge,v1,v3,v5
c end of hydrodynamics
c c2.2 thermodynamics
      DIMENSION ddroe1(kk),ddroe3(kk),ddroi1(kk),ddroi3(kk),ddte1(kk),
     &          ddte3(kk),ddti1(kk),ddti3(kk),xappae(kk),xappai(kk),
     &          pe1(kk),pe3(kk),pi1(kk),pi3(kk),te1(kk),te3(kk),ti1(kk),
     &          ti3(kk)
      COMMON /comth / ddroe1,ddroe3,ddroi1,ddroi3,ddte1,ddte3,ddti1,
     &                ddti3,gammae,gammai,xappae,xappai,pe1,pe3,pi1,pi3,
     &                te1,te3,teini,ti1,ti3,tiini
c end of thermodynamics
c c3.2 mesh and numerical methods
      DIMENSION ae(kk),ai(kk),be(kk),bi(kk),ce(kk),ci(kk),de(kk),di(kk),
     &          e(kk),f(kk),ge(kk),gi(kk),q2(kk),q4(kk),teite(kk),
     &          tiite(kk),uite(kk)
      COMMON /comnum/ mesh,nj,njm1,njp1,nl,nlm1,nlp1,ae,ai,be,bi,ce,ci,
     &                de,di,e,f,ge,gi,q2,q4,teite,tiite,uite
c end of mesh and numerical methods
c c5.1 input-output control variables
      LOGICAL nlfilm,nlhcpy,nlprnt
      DIMENSION buf1(kk),buf2(kk),buf3(kk),buf4(kk),buf5(kk),buf6(kk)
      REAL nprnt,nproc
      COMMON /comout/ nfilm,nhdcpy,np1,np2,np3,nlfilm,nlhcpy,nlprnt,
     &        ztnext,nprnt,nproc,buf1,buf2,buf3,buf4,buf5,buf6,scp,
     &                scr,scrho,scte,scti,sctime,nt1,nt2
c end of input-output control variables
c 1. scale variables. construct logarithms.
      DO 100 l=1,nl
       j=l
       buf1(l)=r1(j)*scr
       buf2(l)=u2(j)*scr/sctime
       z3=rho1(l)*scrho
       buf3(l)=log10(z3)
       z4=p3(l)*scp
c ad1
       buf4(l)=z4
       z5=ti1(l)*scti
       buf5(l)=log10(z5)
       z6=te1(l)*scte
       buf6(l)=log10(z6)
 100  CONTINUE
      RETURN
      END


      SUBROUTINE shifts
      PARAMETER(kk=1001)
c 2.23 shifts the levels 1-5 back
c @(#)comhot    1.2 20/5/85 19:47:07 (20/5/85 19:48:45) @(#)
c c1.9 development and diagnostic parameters
      DIMENSION nadump(20),npdump(20),nvdump(20)
      LOGICAL nlched,nlhead(9),nlomt1(50),nlomt2(50),nlomt3(50),nlrept
      COMMON /comddp/ maxdum,mxdump,nadump,nclass,npdump,npoint,nsub,
     &                nvdump,nlched,nlhead,nlomt1,nlomt2,nlomt3,nlrept
c end of development and diagnostic parameters
c c2.1 hydrodynamics
      DIMENSION dm(kk),p3(kk),r1(kk),r3(kk),r5(kk),rho1(kk),rho3(kk),
     &          rho5(kk),u2(kk),u4(kk),v1(kk),v3(kk),v5(kk),lstat(kk)
      COMMON /comhyd/ dm,p3,pini,r1,r3,r5,rho1,rho3,rho5,rhoini,rhor,
     &                rini,time,u2,u4,uedge,v1,v3,v5
c end of hydrodynamics
c c2.2 thermodynamics
      DIMENSION ddroe1(kk),ddroe3(kk),ddroi1(kk),ddroi3(kk),ddte1(kk),
     &          ddte3(kk),ddti1(kk),ddti3(kk),xappae(kk),xappai(kk),
     &          pe1(kk),pe3(kk),pi1(kk),pi3(kk),te1(kk),te3(kk),ti1(kk),
     &          ti3(kk)
      COMMON /comth / ddroe1,ddroe3,ddroi1,ddroi3,ddte1,ddte3,ddti1,
     &                ddti3,gammae,gammai,xappae,xappai,pe1,pe3,pi1,pi3,
     &                te1,te3,teini,ti1,ti3,tiini
c end of thermodynamics
c c2.3 ions and electrons
      DIMENSION brems1(kk),brems3(kk),ddz1(kk),ddz3(kk),degen(kk),
     &          effz(kk),eixch2(kk),fidash(kk),fz1(kk),fz3(kk),fzsq1(kk)
     &          ,fzsq3(kk),xlc(kk),xmieff(kk),xne(kk),xni(kk),omega1(kk)
      COMMON /comie / brems1,brems3,ddz1,ddz3,degen,degmax,degmin,effz,
     &                eixch2,fidash,fz1,fz3,fzsq1,fzsq3,xlc,xmieff,xne,
     &                xni,omega1,pmass
c end of ions and electrons
c c2.4 laser variables
      DIMENSION alpha1(kk),xaser1(kk),xl1(kk),xl3(kk)
      COMMON /comlas/ alpha1,elas1,xamda1,xaser1,xecri1,plas1,rabs1,
     &                rocri1,xl1,xl3,nabs1
c end of laser variables
c c2.5 thermonuclear reactions
      DIMENSION f1d(kk),f1h(kk),f1he3(kk),f1he4(kk),f1neu(kk),f1ntrl(kk)
     &          ,f1t(kk),f1x(kk),f1x1(kk),f1x2(kk),f1x3(kk),f1x4(kk),
     &          f3d(kk),f3h(kk),f3he3(kk),f3he4(kk),f3neu(kk),f3ntrl(kk)
     &          ,f3t(kk),f3x(kk),f3x1(kk),f3x2(kk),f3x3(kk),f3x4(kk),
     &          r1dd(kk),r1dhe3(kk),r1dt(kk),r3dd(kk),r3dhe3(kk),
     &          r3dt(kk),ye1(kk),ye3(kk),yi1(kk),yi3(kk)
      COMMON /comfus/ deuter,f1d,f1h,f1he3,f1he4,f1neu,f1ntrl,f1t,f1x,
     &                f1x1,f1x2,f1x3,f1x4,f3d,f3h,f3he3,f3he4,f3neu,
     &                f3ntrl,f3t,f3x,f3x1,f3x2,f3x3,f3x4,heliu3,heliu4,
     &                hydrog,xetral,xtrlms,pneut1,pneut3,r1dd,r1dhe3,
     &                r1dt,r3dd,r3dhe3,r3dt,rneut1,rneut3,tinucl,totneu,
     &                tritiu,xmass,xmass1,xmass2,xmass3,xmass4,xtra,
     &                xtra1,xtra2,xtra3,xtra4,xz,xz1,xz2,xz3,xz4,ye1,
     &                ye3,yi1,yi3,yield
c end of thermonuclear reactions
c c2.6 - energies
      COMMON /comen / eerror,en,einput,eloss,eefuse,eifuse,eindt1,
     &                eindt3,eneutr,pv,usqm,xionen,en0,eions,ecbbs,
     &                ecbfs,edcvs,erbfs,eatis,ezsts,erbbs,erfbs,edfbs
c c2.6 - end energies
c c2.7 physics control
      LOGICAL mlbrms,mlecon,mlfuse,mlicon,mlx,nlabs,nlbrms,nlburn,
     &        nlcri1,nldepo,nlecon,nlfuse,nlicon,nlmove,nlpfe,nlpfi,nlx,
     &        nlhf,nltnl
      COMMON /comcon/ mstep,ncase,ngeom,mlbrms,mlecon,mlfuse,mlicon,mlx,
     &                nlabs,nlbrms,nlburn,nlcri1,nldepo,nlecon,nlfuse,
     &                nlicon,nlmove,nlpfe,nlpfi,nlx,nlhf,nltnl,pmin,
     &                rhomin,temin,timin,umin
c end of physics control
c c3.1 numerical control parameters
      LOGICAL break,nlgoon,nlite
      COMMON /comnc / nceldt,ncondt,nit,nitmax,break,nlgoon,nlite,ak0,
     &                ak1,ak2,ak3,ak4,ak5,ak6,ak7,ak8,ak9,bneum,deltat,
     &                dt2,dt3,dt4,dtemax,dtimax,dtmax,dtmin,dumax,rdt2,
     &                rdt3,rdt4
c end of numerical control parameters
c c3.2 mesh and numerical methods
      DIMENSION ae(kk),ai(kk),be(kk),bi(kk),ce(kk),ci(kk),de(kk),di(kk),
     &          e(kk),f(kk),ge(kk),gi(kk),q2(kk),q4(kk),teite(kk),
     &          tiite(kk),uite(kk)
      COMMON /comnum/ mesh,nj,njm1,njp1,nl,nlm1,nlp1,ae,ai,be,bi,ce,ci,
     &                de,di,e,f,ge,gi,q2,q4,teite,tiite,uite
c end of mesh and numerical methods
c c4.1 - administrative variables
c c4.1 - end administrative variables
      COMMON /atom  / pu(10,kk),p(10,kk),pz(10,kk),zs(10,kk),zst(kk),
     &                zstz(kk),z(kk),a(kk),fn(10),deg(10),qs(10,kk),
     &                zstmi(kk),zstu(kk),engu(10,kk),eng(10,kk),
     &                rdrhe(kk),zstin(kk),nmax
      COMMON /dpi   / ecbb(kk),ecbf(kk),eraddr(kk),erbf(kk),eati(kk),
     &                dcvz(kk),efbx2(kk),efrx2(kk),erbb(kk),erfb(kk),
     &                acfc(10,kk),acfr(10,kk),bcfc(10,kk),bcfr(10,kk),
     &                edfb(kk),ter(kk),dqdtim(kk),nmaxr
      COMMON /radflg/ yamda,xhnu,xihnu,dihnu,xplas1(kk),xrbf(kk),
     &                atrbf(kk),corhf1(kk),corhf(kk),iradf1,iradf2
      DATA iclass,isub/2,23/
      IF(nlomt2(isub))RETURN
c ionisation
      CALL copyr(zst,1,zstu,1,nlp1)
      DO 100 l=1,mesh
       DO 50 i=1,nmax
        pu(i,l)=p(i,l)
        engu(i,l)=eng(i,l)
 50    CONTINUE
c hf correction
       CALL copyr(corhf,1,corhf1,1,nlp1)
 100  CONTINUE
      CALL resetr(dqdtim,nl,0.0)
      CALL resetr(efbx2,nl,0.0)
      CALL resetr(efrx2,nl,0.0)
      CALL resetr(dcvz,nl,0.0)
      CALL resetr(ecbb,nl,0.0)
      CALL resetr(ecbf,nl,0.0)
      CALL resetr(eraddr,nl,0.0)
      CALL resetr(erbf,nl,0.0)
      CALL resetr(erbb,nl,0.0)
      CALL resetr(erfb,nl,0.0)
      CALL resetr(edfb,nl,0.0)
c 2. physics
c 2.1 hydrodynamics
      CALL copyr(r3,1,r1,1,njp1)
      CALL copyr(r5,1,r3,1,njp1)
      CALL copyr(v3,1,v1,1,nlp1)
      CALL copyr(v5,1,v3,1,nlp1)
      CALL copyr(rho3,1,rho1,1,nlp1)
      CALL copyr(rho5,1,rho3,1,nlp1)
      CALL copyr(u4,1,u2,1,njp1)
c 2.2 thermodynamics
      CALL copyr(ddroe3,1,ddroe1,1,nlp1)
      CALL copyr(ddroi3,1,ddroi1,1,nlp1)
      CALL copyr(ddte3,1,ddte1,1,nlp1)
      CALL copyr(ddti3,1,ddti1,1,nlp1)
      CALL copyr(ti3,1,ti1,1,nlp1)
      CALL copyr(te3,1,te1,1,nlp1)
      CALL copyr(pi3,1,pi1,1,nlp1)
      CALL copyr(pe3,1,pe1,1,nlp1)
      CALL copyr(fz3,1,fz1,1,nlp1)
      CALL copyr(fzsq3,1,fzsq1,1,nlp1)
      CALL copyr(ddz3,1,ddz1,1,nlp1)
c if heat conduction is switched off clear xappa
      IF(.not.(nlicon.and.mlicon))CALL resetr(xappai,nl,0.0)
      IF(.not.(nlecon.and.mlecon))CALL resetr(xappae,nl,0.0)
c 2.3 ions and electrons
      CALL copyr(brems3,1,brems1,1,nlp1)
c clear the array because radiation may change
      CALL resetr(brems3,nl,0.0)
      IF(.not.(nlx.and.mlx))THEN
       CALL resetr(eixch2,nl,0.0)
       CALL resetr(omega1,nl,0.0)
      ENDIF
c multi group hot electrons
c 2.4 laser variables
      IF(ncase.eq.1)THEN
       CALL copyr(xl3,1,xl1,1,nlp1)
c clear the array because the absorption takes place at different
c positions according to the density profile
       CALL resetr(xl3,nl,0.0)
       CALL resetr(alpha1,nlp1,0.0)
      ENDIF
c 2.5      thermonuclear reactions
c is fusion switched on ?
      IF(nlfuse)THEN
c is burning of deuterium and tritium switched on ?
       IF(nlburn)THEN
        CALL copyr(f3d,1,f1d,1,nlp1)
        CALL copyr(f3h,1,f1h,1,nlp1)
        CALL copyr(f3he3,1,f1he3,1,nlp1)
        CALL copyr(f3he4,1,f1he4,1,nlp1)
        CALL copyr(f3neu,1,f1neu,1,nlp1)
        CALL copyr(f3ntrl,1,f1ntrl,1,nlp1)
        CALL copyr(f3t,1,f1t,1,nlp1)
        CALL copyr(f3x,1,f1x,1,nlp1)
        CALL copyr(f3x1,1,f1x1,1,nlp1)
        CALL copyr(f3x2,1,f1x2,1,nlp1)
        CALL copyr(f3x3,1,f1x3,1,nlp1)
        CALL copyr(f3x4,1,f1x4,1,nlp1)
       ENDIF
       pneut1=pneut3
       pneut3=0.0
       rneut1=rneut3
       rneut3=0.0
       CALL copyr(r3dd,1,r1dd,1,nlp1)
       CALL copyr(r3dhe3,1,r1dhe3,1,nlp1)
       CALL copyr(r3dt,1,r1dt,1,nlp1)
       CALL copyr(yi3,1,yi1,1,nlp1)
       CALL copyr(ye3,1,ye1,1,nlp1)
c clear the arrays because some cells may increase/decrease their
c ion temperature above/below the ignition temperature (tinucl)
       CALL resetr(r3dd,nlp1,0.0)
       CALL resetr(r3dhe3,nlp1,0.0)
       CALL resetr(r3dt,nlp1,0.0)
       CALL resetr(ye3,nlp1,0.0)
       CALL resetr(yi3,nlp1,0.0)
      ENDIF
c 2.6 energies
      eindt1=eindt3
      eindt3=0.0
c 3. numerical scheme
c 3.1 numerical control parameters
      dt2=dt4
      rdt2=rdt4
c 3.2 mesh and numerical methods
      CALL copyr(q4,1,q2,1,nlp1)
      RETURN
      END


      SUBROUTINE signi(ka,kdim)
c u.26 change the sign of an integer array
      DIMENSION ka(kdim)
      DO 100 j=1,kdim
       ka(j)=-ka(j)
 100  CONTINUE
      RETURN
      END


      SUBROUTINE signr(pa,kdim)
c u.25 change the sign of a real array
      DIMENSION pa(kdim)
      DO 100 j=1,kdim
       pa(j)=-pa(j)
 100  CONTINUE
      RETURN
      END


      SUBROUTINE source
      PARAMETER(kk=1001,underf=1.E-30)
c 1.9 change the standard set of initial conditions
c c1.1 basic system parameters
c end of basic system parameters
c c1.9 development and diagnostic parameters
      DIMENSION nadump(20),npdump(20),nvdump(20)
      LOGICAL nlched,nlhead(9),nlomt1(50),nlomt2(50),nlomt3(50),nlrept
      COMMON /comddp/ maxdum,mxdump,nadump,nclass,npdump,npoint,nsub,
     &                nvdump,nlched,nlhead,nlomt1,nlomt2,nlomt3,nlrept
c end of development and diagnostic parameters
c c2.1 hydrodynamics
      DIMENSION dm(kk),p3(kk),r1(kk),r3(kk),r5(kk),rho1(kk),rho3(kk),
     &          rho5(kk),u2(kk),u4(kk),v1(kk),v3(kk),v5(kk),lstat(kk)
      COMMON /comhyd/ dm,p3,pini,r1,r3,r5,rho1,rho3,rho5,rhoini,rhor,
     &                rini,time,u2,u4,uedge,v1,v3,v5
c end of hydrodynamics
c c2.2 thermodynamics
c end of thermodynamics
c c2.5 thermonuclear reactions
      DIMENSION f1d(kk),f1h(kk),f1he3(kk),f1he4(kk),f1neu(kk),f1ntrl(kk)
     &          ,f1t(kk),f1x(kk),f1x1(kk),f1x2(kk),f1x3(kk),f1x4(kk),
     &          f3d(kk),f3h(kk),f3he3(kk),f3he4(kk),f3neu(kk),f3ntrl(kk)
     &          ,f3t(kk),f3x(kk),f3x1(kk),f3x2(kk),f3x3(kk),f3x4(kk),
     &          r1dd(kk),r1dhe3(kk),r1dt(kk),r3dd(kk),r3dhe3(kk),
     &          r3dt(kk),ye1(kk),ye3(kk),yi1(kk),yi3(kk)
      COMMON /comfus/ deuter,f1d,f1h,f1he3,f1he4,f1neu,f1ntrl,f1t,f1x,
     &                f1x1,f1x2,f1x3,f1x4,f3d,f3h,f3he3,f3he4,f3neu,
     &                f3ntrl,f3t,f3x,f3x1,f3x2,f3x3,f3x4,heliu3,heliu4,
     &                hydrog,xetral,xtrlms,pneut1,pneut3,r1dd,r1dhe3,
     &                r1dt,r3dd,r3dhe3,r3dt,rneut1,rneut3,tinucl,totneu,
     &                tritiu,xmass,xmass1,xmass2,xmass3,xmass4,xtra,
     &                xtra1,xtra2,xtra3,xtra4,xz,xz1,xz2,xz3,xz4,ye1,
     &                ye3,yi1,yi3,yield
c end of thermonuclear reactions
c c3.2 mesh and numerical methods
      DIMENSION ae(kk),ai(kk),be(kk),bi(kk),ce(kk),ci(kk),de(kk),di(kk),
     &          e(kk),f(kk),ge(kk),gi(kk),q2(kk),q4(kk),teite(kk),
     &          tiite(kk),uite(kk)
      COMMON /comnum/ mesh,nj,njm1,njp1,nl,nlm1,nlp1,ae,ai,be,bi,ce,ci,
     &                de,di,e,f,ge,gi,q2,q4,teite,tiite,uite
c end of mesh and numerical methods
c c4.1 - administrative variables
      LOGICAL nldump,nlemp
      DIMENSION piq(302)
      COMMON /comadm/ maxdim,maxrun,ndump,nrep,nldump,nlemp,piq,tstop
c c4.1 - end administrative variables
c c5.1 input-output control variables
c end of input-output control variables
c   new subroutine source control factors and layer materials
      COMMON /sornew/ rf1,rf2,rf3,fne2
      njg=njp1-nint(piq(13))-nint(piq(63))
      njp=njp1-nint(piq(63))
      njpp1=njp+1
      njpm1=njp-1
      njgp1=njg+1
      njgm1=njg-1
      IF(abs(piq(13)).lt.underf.and.abs(piq(63)).lt.underf.and.
     &   abs(piq(14)).lt.underf)RETURN
      IF(npoint.eq.1)THEN
      ELSEIF(npoint.eq.2)THEN
c-------------------------------------------------------------------
c set up the densities
       DO 50 l=1,njgm1
        rho1(l)=rhoini
 50    CONTINUE
       IF(abs(piq(11)).gt.underf)THEN
        DO 60 l=njg,njpm1
         rho1(l)=piq(12)
 60     CONTINUE
       ENDIF
       IF(abs(piq(61)).gt.underf)THEN
        DO 80 l=njp,nl
         rho1(l)=piq(62)
 80     CONTINUE
       ENDIF
       GOTO 400
      ELSEIF(npoint.eq.3.or.npoint.eq.4.or.npoint.eq.5.or.
     &       npoint.eq.7.or.npoint.eq.8.or.npoint.eq.9)THEN
       GOTO 400
      ELSEIF(npoint.eq.6)THEN
c set up the chemical composition
c since we are changing the initial conditions from inital
c we first reset the standard medusa variables to zero
       CALL resetr(f1d,nlp1,0.0)
       CALL resetr(f1h,nlp1,0.0)
       CALL resetr(f1t,nlp1,0.0)
       CALL resetr(f1x,nlp1,0.0)
c set up the atomic and mass numbers for silicon, oxygen and carbon
c     if (xz .lt. 1.1) xz = 10.0
c     if (xz1 .lt. 1.1) xz1 = 6.0
c     if (xz2 .lt. 1.1) xz2 = 8.0
       xz3=14.0
       xz4=6.0
c     if (xmass .lt. 1.1) xmass = 20.179
c     if (xmass1 .lt. 1.1) xmass1 = 12.011
c     if (xmass2 .lt. 1.1) xmass2 = 16.0
       xmass3=28.0855
       xmass4=12.011
       zf1d=deuter-piq(14)/2.0
       zf1t=tritiu-piq(14)/2.0
       zf1x=piq(14)
       DO 100 l=1,njgm1
        f1d(l)=zf1d
        f1t(l)=zf1t
        f1x(l)=zf1x
 100   CONTINUE
       IF(abs(piq(11)).gt.underf)THEN
        DO 120 l=njg,njpm1
         f1x2(l)=fne2
         f1x3(l)=1.0-fne2
 120    CONTINUE
       ENDIF
       IF(abs(piq(61)).gt.underf)THEN
        DO 140 l=njp,nl
         f1h(l)=hydrog
         f1x4(l)=1.0-hydrog
         IF(hydrog.le.0.0)THEN
c in this case the material is xz1
          f1x4(l)=0.0
          f1h(l)=0.0
          f1x1(l)=1.0
         ENDIF
 140    CONTINUE
       ENDIF
       GOTO 500
      ELSE
       CALL gotoer(' MEDUSA: source: goto error, npoint is',npoint)
      ENDIF
c--------------------------------------------------------------------
c                      ************
c                      * gridding *
c                      ************
c--------------------------------------------------------------------
c     layer 1 (the innermost layer or gass fill)
      IF(1.gt.njg)THEN
       idodoa=1
      ELSE
       idodoa=njgm1
      ENDIF
      IF((rf1.le.0.0).or.((rf1.gt.0.99999).and.(rf1.lt.1.00001)))
     &   rf1=0.99999
      ncells=idodoa
      r1(1)=0.0
      DO 200 j=1,idodoa
       fact=((rf1**real(j))-1.0)/((rf1**real(ncells))-1.0)
       r1(j+1)=rini*fact
 200  CONTINUE
c     layer 2 (the middle layer or glass layer)
      IF(abs(piq(11)).gt.underf)THEN
       IF((rf2.le.0.0).or.((rf2.gt.0.99999).and.(rf2.lt.1.00001)))
     &    rf2=0.99999
       IF(njgp1.gt.njp)THEN
        idodoa=njg
       ELSE
        idodoa=njp-njg
       ENDIF
       ioffst=njg
       roffst=r1(ioffst)
       ncells=idodoa
       DO 250 j=1,idodoa
        fact=((rf2**real(j))-1.0)/((rf2**real(ncells))-1.0)
        r1(j+ioffst)=(piq(11)*fact)+roffst
 250   CONTINUE
      ENDIF
c     layer 3 (the outer layer or plastic layer)
      IF(abs(piq(61)).gt.underf)THEN
       IF(njpp1.gt.njp1)THEN
        idodoa=njp
       ELSE
        idodoa=njp1-njp
       ENDIF
       IF((rf3.le.0.0).or.((rf3.gt.0.99999).and.(rf3.lt.1.00001)))
     &    rf3=0.99999
       ioffst=njp
       roffst=r1(ioffst)
       ncells=idodoa
       DO 300 j=1,idodoa
        fact=((rf3**real(j))-1.0)/((rf3**real(ncells))-1.0)
        r1(j+ioffst)=(piq(61)*fact)+roffst
 300   CONTINUE
      ENDIF
      RETURN
 400  RETURN
c copy all quantities to to time level 3
 500  CALL copyr(f1d,1,f3d,1,nlp1)
      CALL copyr(f1t,1,f3t,1,nlp1)
      CALL copyr(f1h,1,f3h,1,nlp1)
      CALL copyr(f1x,1,f3x,1,nlp1)
      CALL copyr(f1x1,1,f3x1,1,nlp1)
      CALL copyr(f1x2,1,f3x2,1,nlp1)
      CALL copyr(f1x3,1,f3x3,1,nlp1)
      CALL copyr(f1x4,1,f3x4,1,nlp1)
      RETURN
      END


      SUBROUTINE speed
      PARAMETER(kk=1001)
c 2.17 the hydrodynamic velocities
c c1.9 development and diagnostic parameters
      DIMENSION nadump(20),npdump(20),nvdump(20)
      LOGICAL nlched,nlhead(9),nlomt1(50),nlomt2(50),nlomt3(50),nlrept
      COMMON /comddp/ maxdum,mxdump,nadump,nclass,npdump,npoint,nsub,
     &                nvdump,nlched,nlhead,nlomt1,nlomt2,nlomt3,nlrept
c end of development and diagnostic parameters
c c4.1 - administrative variables
      LOGICAL nldump,nlemp
      DIMENSION piq(302)
      COMMON /comadm/ maxdim,maxrun,ndump,nrep,nldump,nlemp,piq,tstop
c c4.1 - end administrative variables
c c2.1 hydrodynamics
      DIMENSION dm(kk),p3(kk),r1(kk),r3(kk),r5(kk),rho1(kk),rho3(kk),
     &          rho5(kk),u2(kk),u4(kk),v1(kk),v3(kk),v5(kk),lstat(kk)
      COMMON /comhyd/ dm,p3,pini,r1,r3,r5,rho1,rho3,rho5,rhoini,rhor,
     &                rini,time,u2,u4,uedge,v1,v3,v5
c end of hydrodynamics
c c2.7 physics control
      LOGICAL mlbrms,mlecon,mlfuse,mlicon,mlx,nlabs,nlbrms,nlburn,
     &        nlcri1,nldepo,nlecon,nlfuse,nlicon,nlmove,nlpfe,nlpfi,nlx,
     &        nlhf,nltnl
      COMMON /comcon/ mstep,ncase,ngeom,mlbrms,mlecon,mlfuse,mlicon,mlx,
     &                nlabs,nlbrms,nlburn,nlcri1,nldepo,nlecon,nlfuse,
     &                nlicon,nlmove,nlpfe,nlpfi,nlx,nlhf,nltnl,pmin,
     &                rhomin,temin,timin,umin
c end of physics control
c c3.1 numerical control parameters
      LOGICAL break,nlgoon,nlite
      COMMON /comnc / nceldt,ncondt,nit,nitmax,break,nlgoon,nlite,ak0,
     &                ak1,ak2,ak3,ak4,ak5,ak6,ak7,ak8,ak9,bneum,deltat,
     &                dt2,dt3,dt4,dtemax,dtimax,dtmax,dtmin,dumax,rdt2,
     &                rdt3,rdt4
c end of numerical control parameters
c c2.2 thermodynamics
      DIMENSION ddroe1(kk),ddroe3(kk),ddroi1(kk),ddroi3(kk),ddte1(kk),
     &          ddte3(kk),ddti1(kk),ddti3(kk),xappae(kk),xappai(kk),
     &          pe1(kk),pe3(kk),pi1(kk),pi3(kk),te1(kk),te3(kk),ti1(kk),
     &          ti3(kk)
      COMMON /comth / ddroe1,ddroe3,ddroi1,ddroi3,ddte1,ddte3,ddti1,
     &                ddti3,gammae,gammai,xappae,xappai,pe1,pe3,pi1,pi3,
     &                te1,te3,teini,ti1,ti3,tiini
c end of thermodynamics
c c3.2 mesh and numerical methods
      DIMENSION ae(kk),ai(kk),be(kk),bi(kk),ce(kk),ci(kk),de(kk),di(kk),
     &          e(kk),f(kk),ge(kk),gi(kk),q2(kk),q4(kk),teite(kk),
     &          tiite(kk),uite(kk)
      COMMON /comnum/ mesh,nj,njm1,njp1,nl,nlm1,nlp1,ae,ai,be,bi,ce,ci,
     &                de,di,e,f,ge,gi,q2,q4,teite,tiite,uite
c end of mesh and numerical methods
c c1.1 basic system parameters
c end of basic system parameters
c c2.3 ions and electrons
c end of ions and electrons
c c2.4 laser variables
c end of laser variables

      LOGICAL lstat
      DATA iclass,isub/2,17/

      if(mstep.le.1) then
      do 72 l=1,kk
       lstat(l)=.false.
c acoustic
c      lstat(l)=.true.
 72   continue
      endif
c the array lstat is used to prevent mesh points from
c moving until a shock wave or thermal front arrives
c this condition does not apply to the last cell
c acoustic
c     lstat(njp1)=.false.
c     lstat(1)=.false.
      IF(nlomt2(isub))RETURN
c 1. find u from equation 73
      igeom=ngeom-1
      njg=nl-nint(piq(13))-nint(piq(63))
      njp=nl-nint(piq(63))
      DO 100 j=2,njp1
c acoustic
c     DO 100 j=2,nl
       l=j
       zrdt=dt3*r3(j)**igeom*2.0/(dm(l)+dm(l-1))
c      the average values of q at level 3 in cells l - 1 and l
       zqlm1=0.5*(q4(l-1)+q2(l-1))
       zql=0.5*(q4(l)+q2(l))
       u4(j)=u2(j)-zrdt*(p3(l)-p3(l-1)+zql-zqlm1)
c ad1  add ponderomotive force see formp for alternative rad pres
c      gprad1=xplas1(l+1)
c      gprad2=xplas1(l)
c      gprad=(gprad1-gprad2)/(6.e08)
c      zepsln = 0.5*(xne(l)+xne(l-1))/xecri1
c      u4(j)=u4(j)- zrdt  * piq(56)*zepsln*gprad

c if the displacement caused by u4 is too small (below rounding
c off level) we set u4 to zero in order to avoid a meaningless
c check on convergence of iterations
       zumin=(r3(j)-r3(j-1))*rdt3*1.0E-09
       IF(abs(u4(j)).lt.zumin)u4(j)=0.0
c for use with layered or shell targets we prevent the
c mesh from moving until the electron or ion temperatures
c have changed from the initial values
       tein1=teini
       tiin1=tiini
       if(j.ge.njp)then
        tein1=piq(89)
        tiin1=piq(90)
       elseif(l.ge.njg)then
        tein1=piq(86)
        tiin1=piq(87)
       endif
c acoustic
       IF(.not.(lstat(l))) lstat(l)=(abs(te3(l-1)-tein1)/(tein1+temin)
c     &    +abs(ti3(l-1)-tiin1)/(tiin1+timin)).gt.10.01
     &    +abs(ti3(l-1)-tiin1)/(tiin1+timin)).gt.1.01
       IF(.not.lstat(l))u4(j)=0.0
 100  CONTINUE
      IF(ncase.eq.2)u4(njp1)=0.0
      itest=nint(piq(55))
      IF(itest.eq.1.or.itest.eq.3)u4(njp1)=0.0
      IF(itest.lt.2.or.ngeom.ne.1)RETURN
      zrdt=dt3*2.0/dm(1)
      u4(1)=u2(1)-zrdt*(p3(1)+0.5*(q2(1)+q4(1)))
      zumin=(r3(2)-r3(1))*rdt3*1.0E-09
      IF(abs(u4(1)).lt.zumin)u4(1)=0.0
       tein1=teini
       tiin1=tiini
       if(l.gt.njp)then
        tein1=piq(89)
        tiin1=piq(90)
       elseif(l.gt.njg)then
        tein1=piq(86)
        tiin1=piq(87)
       endif
      IF(.not.(lstat(1)))lstat(1)=(abs(te3(1)-tein1)/(tein1+temin)
c     &  +abs(ti3(1)- tiin1)/(tiin1+timin)).gt.10.01
     &  +abs(ti3(1)- tiin1)/(tiin1+timin)).gt.1.01
      IF(.not.lstat(1))u4(1)=0.0
      RETURN
      END


      SUBROUTINE start
      PARAMETER(kk=1001,underf=1.E-30)
c 1.8 start the calculation
c c1.1 basic system parameters
c end of basic system parameters
c c2.1 hydrodynamics
      DIMENSION dm(kk),p3(kk),r1(kk),r3(kk),r5(kk),rho1(kk),rho3(kk),
     &          rho5(kk),u2(kk),u4(kk),v1(kk),v3(kk),v5(kk),lstat(kk)
      COMMON /comhyd/ dm,p3,pini,r1,r3,r5,rho1,rho3,rho5,rhoini,rhor,
     &                rini,time,u2,u4,uedge,v1,v3,v5
c end of hydrodynamics
c c2.2 thermodynamics
      DIMENSION ddroe1(kk),ddroe3(kk),ddroi1(kk),ddroi3(kk),ddte1(kk),
     &          ddte3(kk),ddti1(kk),ddti3(kk),xappae(kk),xappai(kk),
     &          pe1(kk),pe3(kk),pi1(kk),pi3(kk),te1(kk),te3(kk),ti1(kk),
     &          ti3(kk)
      COMMON /comth / ddroe1,ddroe3,ddroi1,ddroi3,ddte1,ddte3,ddti1,
     &                ddti3,gammae,gammai,xappae,xappai,pe1,pe3,pi1,pi3,
     &                te1,te3,teini,ti1,ti3,tiini
c end of thermodynamics
c c2.3 ions and electrons
      DIMENSION brems1(kk),brems3(kk),ddz1(kk),ddz3(kk),degen(kk),
     &          effz(kk),eixch2(kk),fidash(kk),fz1(kk),fz3(kk),fzsq1(kk)
     &          ,fzsq3(kk),xlc(kk),xmieff(kk),xne(kk),xni(kk),omega1(kk)
      COMMON /comie / brems1,brems3,ddz1,ddz3,degen,degmax,degmin,effz,
     &                eixch2,fidash,fz1,fz3,fzsq1,fzsq3,xlc,xmieff,xne,
     &                xni,omega1,pmass
c end of ions and electrons
c c2.4 laser variables
      DIMENSION alpha1(kk),xaser1(kk),xl1(kk),xl3(kk)
      COMMON /comlas/ alpha1,elas1,xamda1,xaser1,xecri1,plas1,rabs1,
     &                rocri1,xl1,xl3,nabs1
c end of laser variables
c c2.5 thermonuclear reactions
      DIMENSION f1d(kk),f1h(kk),f1he3(kk),f1he4(kk),f1neu(kk),f1ntrl(kk)
     &          ,f1t(kk),f1x(kk),f1x1(kk),f1x2(kk),f1x3(kk),f1x4(kk),
     &          f3d(kk),f3h(kk),f3he3(kk),f3he4(kk),f3neu(kk),f3ntrl(kk)
     &          ,f3t(kk),f3x(kk),f3x1(kk),f3x2(kk),f3x3(kk),f3x4(kk),
     &          r1dd(kk),r1dhe3(kk),r1dt(kk),r3dd(kk),r3dhe3(kk),
     &          r3dt(kk),ye1(kk),ye3(kk),yi1(kk),yi3(kk)
      COMMON /comfus/ deuter,f1d,f1h,f1he3,f1he4,f1neu,f1ntrl,f1t,f1x,
     &                f1x1,f1x2,f1x3,f1x4,f3d,f3h,f3he3,f3he4,f3neu,
     &                f3ntrl,f3t,f3x,f3x1,f3x2,f3x3,f3x4,heliu3,heliu4,
     &                hydrog,xetral,xtrlms,pneut1,pneut3,r1dd,r1dhe3,
     &                r1dt,r3dd,r3dhe3,r3dt,rneut1,rneut3,tinucl,totneu,
     &                tritiu,xmass,xmass1,xmass2,xmass3,xmass4,xtra,
     &                xtra1,xtra2,xtra3,xtra4,xz,xz1,xz2,xz3,xz4,ye1,
     &                ye3,yi1,yi3,yield
c end of thermonuclear reactions
c c2.6 - energies
      COMMON /comen / eerror,en,einput,eloss,eefuse,eifuse,eindt1,
     &                eindt3,eneutr,pv,usqm,xionen,en0,eions,ecbbs,
     &                ecbfs,edcvs,erbfs,eatis,ezsts,erbbs,erfbs,edfbs
c c2.6 - end energies
c c2.7 physics control
      LOGICAL mlbrms,mlecon,mlfuse,mlicon,mlx,nlabs,nlbrms,nlburn,
     &        nlcri1,nldepo,nlecon,nlfuse,nlicon,nlmove,nlpfe,nlpfi,nlx,
     &        nlhf,nltnl
      COMMON /comcon/ mstep,ncase,ngeom,mlbrms,mlecon,mlfuse,mlicon,mlx,
     &                nlabs,nlbrms,nlburn,nlcri1,nldepo,nlecon,nlfuse,
     &                nlicon,nlmove,nlpfe,nlpfi,nlx,nlhf,nltnl,pmin,
     &                rhomin,temin,timin,umin
c end of physics control
c c3.1 numerical control parameters
      LOGICAL break,nlgoon,nlite
      COMMON /comnc / nceldt,ncondt,nit,nitmax,break,nlgoon,nlite,ak0,
     &                ak1,ak2,ak3,ak4,ak5,ak6,ak7,ak8,ak9,bneum,deltat,
     &                dt2,dt3,dt4,dtemax,dtimax,dtmax,dtmin,dumax,rdt2,
     &                rdt3,rdt4
c end of numerical control parameters
c c3.2 mesh and numerical methods
      DIMENSION ae(kk),ai(kk),be(kk),bi(kk),ce(kk),ci(kk),de(kk),di(kk),
     &          e(kk),f(kk),ge(kk),gi(kk),q2(kk),q4(kk),teite(kk),
     &          tiite(kk),uite(kk)
      COMMON /comnum/ mesh,nj,njm1,njp1,nl,nlm1,nlp1,ae,ai,be,bi,ce,ci,
     &                de,di,e,f,ge,gi,q2,q4,teite,tiite,uite
c end of mesh and numerical methods
c c?.?? ??????????

c x-ray laser common blocks
c end of ??????????
      iclass=1
      isub=8
c in order to start this calculation it is necessary to work out
c several derived variables at time level 1. to facilitate this we
c use the routines of class 2 which store the derived variables at
c level 3. hence we copy to level 1. however, it is first necessary
c to calculate the initial timestep dt2.
c 1. find initial timestep
c before we do this the initial conditions are checked
      CALL chinic
c the only restriction on dt2 is now set by the c-f-l condition
c (see eq.75)
      ncondt=1
      zd1=0.0
      DO 100 l=1,nl
       j=l
       zz1=sqrt(abs(p3(l))*v1(l))/(r1(j+1)-r1(j))
       IF(zz1.ge.zd1)THEN
        zd1=zz1
        nceldt=l
       ENDIF
 100  CONTINUE
      dt2=ak1/(sqrt(gammai)*zd1)
      zdt2=dt2
      CALL expert(iclass,isub,1)
c 2. take account of changes at the boundary
      IF(ncase.ne.1)THEN
       IF(ncase.eq.2)THEN
c 2.2 temperature determined
c find the boundary temperatures at the first step
        CALL boundy(dt2)
c include threshold values for safety
        zd3=(ti3(nlp1)+timin)/(ti1(nlp1)+timin)
        zd3=max(zd3,1.0/zd3)
        zd4=(te3(nlp1)+temin)/(te1(nlp1)+temin)
        zd4=max(zd4,1.0/zd4)
        dt2=min(dt2,1.0/zd3*dt2,1.0/zd4*dt2)
c put back the old values
        ti3(nlp1)=ti1(nlp1)
        te3(nlp1)=te1(nlp1)
        IF(abs(dt2-zdt2).gt.underf)THEN
         nceldt=nlp1
         ncondt=3
         IF(zd4.gt.zd3)ncondt=4
        ENDIF
       ELSEIF(ncase.eq.3)THEN
c 2.3 pressure determined
        zpsave=p3(nlp1)
c find the boundary pressure at the first step
        CALL boundy(dt2)
        zd1=sqrt(p3(nlp1)*v3(nl)*gammai)/(r3(njp1)-r3(nj))
        dt2=min(dt2,ak1/zd1)
c put back the old value
        p3(nlp1)=zpsave
        IF(abs(dt2-zdt2).gt.underf)THEN
         nceldt=nlp1
         ncondt=5
        ENDIF
       ELSEIF(ncase.eq.4)THEN
c 2.4 velocity determined
c find the boundary velocity at the first step
        tdash=0.5*dt2
        CALL boundy(tdash)
        zd1=abs(u4(njp1))/(r3(njp1)-r3(nj))
        dt2=min(dt2,1.0/zd1)
        IF(abs(dt2-zdt2).gt.underf)THEN
         nceldt=nlp1
         ncondt=5
        ENDIF
       ELSE
        CALL gotoer(' MEDUSA: start: goto error, ncase is',ncase)
        GOTO 200
       ENDIF
       GOTO 300
      ENDIF
c 2.1 thermally insulated wall
 200  CALL laser(dt2)
      IF(abs(plas1).gt.underf)THEN
c find out where the laser light is absorbed
       CALL absorb
       IF(nabs1.gt.1)THEN
c the pressure increase caused by the absorption of plaser is
c thought to cause the following displacement of a point near the
c critical density. (notice spacings dr are assumed constant)
        zdr=r3(nabs1)**(ngeom-1)
     &      *plas1*dt2*dt2*dt2/(dm(nabs1)*(r3(nabs1)-r3(nabs1-1)))
c make sure that the boundaries of the cell nabs can not cross
c the neighbouring cell-boundaries
        zdrabs=min((r3(nabs1+1)-r3(nabs1)),(r3(nabs1)-r3(nabs1-1)))
        dt2=min(dt2,zdrabs/zdr*dt2)
        IF(abs(dt2-zdt2).gt.underf)THEN
         nceldt=nabs1
         ncondt=6
        ENDIF
       ENDIF
c reset the laser power and nabs1
       CALL resetr(alpha1,nlp1,0.0)
       CALL resetr(xl3,nlp1,0.0)
       plas1=xaser1(2)
       nabs1=0
      ENDIF
c 2.5 select the proper value of dt2
c if the found value of dt2 is larger than the value specified in
c subroutine    data we choose the latter
 300  IF(dt2.gt.deltat.and.deltat.gt.0.0)dt2=deltat
      IF(abs(dt2-deltat).lt.underf)ncondt=-1
c 3. precalculate the artificial viscous pressure
c first copy level 3 to 5 and 2 to 4
      CALL copyr(rho3,1,rho5,1,nl)
      IF(nlmove)THEN
       CALL copyr(u2,1,u4,1,njp1)
       CALL neuman
c then bring q4 back to q2
       CALL copyr(q4,1,q2,1,njp1)
c the arrays u4 and rho5 can be left as they are
       CALL expert(iclass,isub,2)
      ENDIF
c absorption of energy
      IF(ncase.eq.1)CALL absorb
c 4. precalculate particle terms and shift to 1
      CALL coulog
c bremsstrahlung
      IF(nlbrms.and.mlbrms)CALL brems
c thermonuclear reactions
      IF(nlfuse.and.mlfuse)CALL fusion
c now shift all terms to level 1
      CALL copyr(xl3,1,xl1,1,nlp1)
      CALL copyr(brems3,1,brems1,1,nlp1)
      CALL copyr(yi3,1,yi1,1,nlp1)
      CALL copyr(ye3,1,ye1,1,nlp1)
c we leave the level 3 arrays as they are
      CALL expert(iclass,isub,3)
c 5. precalculate the thermodynamics
c heat conductivities
      CALL hcduct
c partial derivatives of energy
      CALL statei(1)
      CALL statee(1)
c ion-electron exchange
c the time constant
      IF(nlx.and.mlx)CALL xchang(1)
c the actual rate of energy transferred
      IF(nlx.and.mlx)CALL xchang(2)
c partial derivatives wrt z
      CALL radran
      CALL ionbal
      CALL copyr(fz3,1,fz1,1,nlp1)
      CALL copyr(fzsq3,1,fzsq1,1,nlp1)
c ad1
      CALL statei(2)
      CALL statee(2)
      CALL formp
c now copy derivatives from level 3 to 1
      CALL copyr(ddroe3,1,ddroe1,1,nlp1)
      CALL copyr(ddroi3,1,ddroi1,1,nlp1)
      CALL copyr(ddti3,1,ddti1,1,nlp1)
      CALL copyr(ddte3,1,ddte1,1,nlp1)
      CALL copyr(ddz3,1,ddz1,1,nlp1)
c finally evaluate the coefficients a,b,c and d at level 1. this is
c done so that the coefficients gi and ge can be worked out.
      rdt2=1.0/dt2
      dt3=dt2
      rdt3=rdt2
c ad1 set default dtmax
      IF(abs(dtmax-1.E-18).lt.underf.or.dtmax.le.dt2)dtmax=dt2*1.0E6
c     if ( abs(dtmin).lt.underf .or. dtmin.ge.dt2 ) dtmin = dt2
      CALL abcd
c we are now at time = 0.0, and we can start the calculation from
c level 1, since all quantities are known. we have assumed that no
c motion takes place in advancing energies from level 1 to 3, i.e.
c u2 = 0.0  and  r3 = r1.
      IF(mstep.le.0)THEN
c initially the energy input is the energy contents. gain is zero.
       CALL energy
       einput=en
       en0=en
       eerror=0.
      ENDIF
c the lawson criterion
      zrhor=0.0
      zdm=0.0
      DO 400 l=1,nl
       j=l
       zdm=zdm+dm(l)
       zrhor=zrhor+rho1(l)*(r1(j)+r1(j+1))*dm(l)
 400  CONTINUE
c the factor 0.5 arises from averaging r
      rhor=0.5*zrhor/zdm
      CALL expert(iclass,isub,4)
      RETURN
      END


      SUBROUTINE statee(k)
      PARAMETER(kk=1001,underf=1.E-30)
c 2.7 describes the state of electrons including degeneracy
c c1.9 development and diagnostic parameters
      DIMENSION nadump(20),npdump(20),nvdump(20)
      LOGICAL nlched,nlhead(9),nlomt1(50),nlomt2(50),nlomt3(50),nlrept
      COMMON /comddp/ maxdum,mxdump,nadump,nclass,npdump,npoint,nsub,
     &                nvdump,nlched,nlhead,nlomt1,nlomt2,nlomt3,nlrept
c end of development and diagnostic parameters
c c2.1 hydrodynamics
      DIMENSION dm(kk),p3(kk),r1(kk),r3(kk),r5(kk),rho1(kk),rho3(kk),
     &          rho5(kk),u2(kk),u4(kk),v1(kk),v3(kk),v5(kk),lstat(kk)
      COMMON /comhyd/ dm,p3,pini,r1,r3,r5,rho1,rho3,rho5,rhoini,rhor,
     &                rini,time,u2,u4,uedge,v1,v3,v5
c end of hydrodynamics
c c2.2 thermodynamics
      DIMENSION ddroe1(kk),ddroe3(kk),ddroi1(kk),ddroi3(kk),ddte1(kk),
     &          ddte3(kk),ddti1(kk),ddti3(kk),xappae(kk),xappai(kk),
     &          pe1(kk),pe3(kk),pi1(kk),pi3(kk),te1(kk),te3(kk),ti1(kk),
     &          ti3(kk)
      COMMON /comth / ddroe1,ddroe3,ddroi1,ddroi3,ddte1,ddte3,ddti1,
     &                ddti3,gammae,gammai,xappae,xappai,pe1,pe3,pi1,pi3,
     &                te1,te3,teini,ti1,ti3,tiini
c end of thermodynamics
c c2.3 ions and electrons
      DIMENSION brems1(kk),brems3(kk),ddz1(kk),ddz3(kk),degen(kk),
     &          effz(kk),eixch2(kk),fidash(kk),fz1(kk),fz3(kk),fzsq1(kk)
     &          ,fzsq3(kk),xlc(kk),xmieff(kk),xne(kk),xni(kk),omega1(kk)
      COMMON /comie / brems1,brems3,ddz1,ddz3,degen,degmax,degmin,effz,
     &                eixch2,fidash,fz1,fz3,fzsq1,fzsq3,xlc,xmieff,xne,
     &                xni,omega1,pmass
c end of ions and electrons
c c2.7 physics control
      LOGICAL mlbrms,mlecon,mlfuse,mlicon,mlx,nlabs,nlbrms,nlburn,
     &        nlcri1,nldepo,nlecon,nlfuse,nlicon,nlmove,nlpfe,nlpfi,nlx,
     &        nlhf,nltnl
      COMMON /comcon/ mstep,ncase,ngeom,mlbrms,mlecon,mlfuse,mlicon,mlx,
     &                nlabs,nlbrms,nlburn,nlcri1,nldepo,nlecon,nlfuse,
     &                nlicon,nlmove,nlpfe,nlpfi,nlx,nlhf,nltnl,pmin,
     &                rhomin,temin,timin,umin
c end of physics control
c c3.2 mesh and numerical methods
      DIMENSION ae(kk),ai(kk),be(kk),bi(kk),ce(kk),ci(kk),de(kk),di(kk),
     &          e(kk),f(kk),ge(kk),gi(kk),q2(kk),q4(kk),teite(kk),
     &          tiite(kk),uite(kk)
      COMMON /comnum/ mesh,nj,njm1,njp1,nl,nlm1,nlp1,ae,ai,be,bi,ce,ci,
     &                de,di,e,f,ge,gi,q2,q4,teite,tiite,uite
c end of mesh and numerical methods
c c4.1 - administrative variables
      LOGICAL nldump,nlemp
      DIMENSION piq(302)
      COMMON /comadm/ maxdim,maxrun,ndump,nrep,nldump,nlemp,piq,tstop
c c4.1 - end administrative variables
c c?.?? ??????????
      DIMENSION tfmult(kk)
      COMMON /comtfc/ kstart,tfmult
c     end of ??????????
      DATA iclass,isub/2,7/
      DATA kdum/0/
      kstart=kdum
      IF(nlomt2(isub))RETURN
c     1. set constants
c     the factor  8255.17  is  k / m(proton)
      zkm=8255.17
      z53=5.0/3.0
      z23=gammae-1.0
      zg=1.0/z23
      zdeg=1.05478359*1.0E7
      zgkm=zg*zkm
c     coefficients for expansions in weakly or strongly degenerate cases
      zp1=0.4
      zp2=1.6449341
      zp3=-2.4352273
      zpw=0.1329808
      zc1=2.0*zp2
      zc2=4.0*zp3
      zcw=-zpw/2.0
c     the coeff zbmp1 - zbmp3 arise from taking a difference
c     as explained below
      zbmp1=0.0
      IF(abs(gammae-z53).gt.underf)zbmp1=1.0-2.0*zg/3.0
      zbmp2=(2.0*zg+3.0)*1.3707784
      zbmp3=-(2.0*zg+1.0)*6.0880682
      zwbmp1=(1.0-zg)*zpw

c     if 1
      IF(.not.nlpfe.and.kstart.ne.1143.and.piq(59).gt.0.5)THEN
c     calculate TF correction factors
c     write(6, 81) piq(59)
       IF(piq(59).le.1.5)THEN
        DO 20 l=1,nl
         tfmult(l)=-1.0
 20     CONTINUE
       ELSEIF(piq(59).le.2.5)THEN
        tmin=1.E10
        DO 40 l=1,nl
         IF(te3(l).lt.tmin)tmin=te3(l)
         tfmult(l)=1.0
c     if (tmin .lt. 9.9e3) write(6, 1031)
c     if (tmin .lt. 9.9e3) write(6, 1032)
 40     CONTINUE
       ELSE
        tmin=1.E10
        tpini=1.0E09
        DO 60 l=1,nl
         tfrho=rho3(l)
         IF(tfrho.lt.101.0)THEN
          tfmult(l)=-1.0
         ELSE
          tft=te3(l)
          IF(tft.lt.tmin)tmin=tft
          tfv=xmieff(l)*pmass/tfrho
          tfz=effz(l)
          tfp1=tfp(tft,tfv,tfz,-1.0)
          tfp2=tfp(tft,tfv,tfz,1.0)
          tfpion=zkm*tfrho*ti3(l)/xmieff(l)
          tfp1=tfp1+tfpion
          tfp2=tfp2+tfpion
          tfmult(l)=4.0
          IF(abs(tfp1-tfp2).gt.1.0)tfmult(l)=(tfp1-tpini)/(tfp1-tfp2)
          IF(tfmult(l).gt.3.0)tfmult(l)=3.0
         ENDIF
 60     CONTINUE
c     if (tmin .lt. 9.9e3) write(6, 1031)
c     if (tmin .lt. 9.9e3) write(6, 1032)
        WRITE(6,10600)
        WRITE(6,10300)(tfmult(l),l=1,nl)
       ENDIF
      ENDIF
c     end if 1

      kstart=1143
      kdum=1143
c     evaluate fermi temp (eq. 43)
      DO 100 l=1,nl
c      zz             = (effz(l) * rho3(l) / xmieff(l)) ** z23
c      ad1 fermi energy
       zz=(fz3(l)*rho3(l)/xmieff(l))**z23
       zfermi=2.5*zdeg/zkm*zz
c      the parameter measuring degeneracy (eq.42)
       degen(l)=te3(l)/zfermi
 100  CONTINUE

c     perfect gas eos?
c     if 2
      IF(nlpfe)THEN
       IF(k.eq.1)THEN
c      5. ideal gas laws. derivatives
        DO 120 l=1,nl
         ddte3(l)=zgkm*fz3(l)/xmieff(l)
         ddroe3(l)=0.0
 120    CONTINUE
        RETURN
       ELSEIF(k.eq.2)THEN
c      6. ideal gas laws. pressure.
        DO 140 l=1,nl
         pe3(l)=zkm*rho3(l)*te3(l)*fz3(l)/xmieff(l)
 140    CONTINUE
        GOTO 300
       ELSE
        CALL gotoer(' MEDUSA: statee: goto error, k is',k)
       ENDIF
       RETURN
      ENDIF
c     end if 2

c     NON-perfect gas eos
c     if 3
      IF(piq(59).lt.0.5)THEN
c     are electrons degenerate?
       IF(k.eq.1)THEN
       ELSEIF(k.eq.2)THEN
c      4. degenerate case. pressure from eq.44
        DO 160 l=1,nl
         zd=degen(l)
         zpe=zkm*rho3(l)*te3(l)*fz3(l)/xmieff(l)
         IF(zd.ge.degmax)THEN
c         4.3 the non-degenerate case
          pe3(l)=zpe
         ELSEIF(zd.ge.degmin)THEN
c         4.2 the weakly degenerate case
          pe3(l)=zpe*(1.0+zpw/(zd*sqrt(zd)))
         ELSE
c         4.1 the strongly degenerate case
          pe3(l)=zpe*(zp1/zd+zp2*zd+zp3*zd*zd*zd)
         ENDIF
 160    CONTINUE
        GOTO 300
       ELSE
        CALL gotoer(' MEDUSA: statee: goto error, k is',k)
       ENDIF
c      3. degenerate case. derivatives from eqs. 45 & 46
c      in the degenerate case we combine the terms ddroe * drho / dt and
c      pe * dv / dt. the former is written as -ddroe / (v * v) * dv / dt
c      and the difference pe - ddroe / (v * v) is stored in ddroe.
       DO 200 l=1,nl
        zd=degen(l)
        zzr=fz3(l)*rho3(l)/xmieff(l)
c       the pressure of a completely degenerate gas (eq.47)
        zpdeg=zdeg*zzr**gammae
        zddte=zg*zkm*fz3(l)/xmieff(l)
        zpe=zkm*rho3(l)*te3(l)*fz3(l)/xmieff(l)
        IF(zd.ge.degmax)THEN
c        3.3 the non-degenerate case
         ddte3(l)=zddte
c        in the non-degenerate case only pe is non-zero
         ddroe3(l)=zpe
        ELSEIF(zd.ge.degmin)THEN
c        3.2 the weakly degenerate case
         ddte3(l)=zddte*(1.0+zcw/(zd*sqrt(zd)))
         ddroe3(l)=zpe*(1.0+zwbmp1/(zd*sqrt(zd)))
        ELSE
c        3.1 the strongly degenerate case
         ddte3(l)=zddte*(zc1*zd+zc2*zd*zd*zd)
         ddroe3(l)=zpdeg*(zbmp1+zbmp2*zd*zd+zbmp3*zd*zd*zd*zd)
        ENDIF
 200   CONTINUE
       RETURN
      ELSE
c     Thomas Fermi
       DO 250 l=1,nl
        tft=te3(l)
        tfv=xmieff(l)*pmass/rho3(l)
        tfz=effz(l)
        tfc=tfmult(l)
        IF(k.ne.2)THEN
         CALL tfdu(tft,tfv,tfz,tfc,dudt,dudv)
c        ddroe3 contains the pressure term too. just as in degen case
         ddroe3(l)=dudv
         ddte3(l)=dudt/(pmass*xmieff(l))
        ELSE
         pe3(l)=tfp(tft,tfv,tfz,tfc)
        ENDIF
 250   CONTINUE
       RETURN
      ENDIF
c     end if 3

 300  RETURN
10100 FORMAT(' piq(59) =',f5.2/)
10200 FORMAT('  kstart, k, piq(59)',2I10,f10.2)
10300 FORMAT(' ',10F10.4)
10400 FORMAT('   *************************************************'/
     &       '   *************************************************'/
     &       '   *************************************************'/
     &       '   ***                                           ***'/
     &       '   ***     if this run fails to converge try     ***'/
     &       '   ***  raising your starting temp. to a value   ***'/
     &       '   ***   nearer to that required for marginal    ***'/
     &       '   ***     degeneracy (1.e4 is a safe value)     ***'/
     &       '   ***                                           ***')
10500 FORMAT('   *************************************************'/
     &       '   *************************************************'/
     &       '   *************************************************'////)
10600 FORMAT(' the tfc quantum',
     &       ' corrections have been added to the basic'/
     &       ' tf quantities in the following amounts in each cell.'/
     &       ' in the exact tfc theory the following numbers would all'/
     &       ' be unity, but the following values are chosen to obtain'/
     &       ' the required solid density.  any values differing from'/
     &     ' unity by more than a factor of two should be investigated.'
     &     /' a value of -1 means that the corrections have been'/
     &     ' left out completely because the material is gaseous.'///)
      END


      SUBROUTINE statei(k)
      PARAMETER(kk=1001)
c 2.6 describes the state of ions (hydrogen)
c c1.9 development and diagnostic parameters
      DIMENSION nadump(20),npdump(20),nvdump(20)
      LOGICAL nlched,nlhead(9),nlomt1(50),nlomt2(50),nlomt3(50),nlrept
      COMMON /comddp/ maxdum,mxdump,nadump,nclass,npdump,npoint,nsub,
     &                nvdump,nlched,nlhead,nlomt1,nlomt2,nlomt3,nlrept
c end of development and diagnostic parameters
c c2.1 hydrodynamics
      DIMENSION dm(kk),p3(kk),r1(kk),r3(kk),r5(kk),rho1(kk),rho3(kk),
     &          rho5(kk),u2(kk),u4(kk),v1(kk),v3(kk),v5(kk),lstat(kk)
      COMMON /comhyd/ dm,p3,pini,r1,r3,r5,rho1,rho3,rho5,rhoini,rhor,
     &                rini,time,u2,u4,uedge,v1,v3,v5
c end of hydrodynamics
c c2.2 thermodynamics
      DIMENSION ddroe1(kk),ddroe3(kk),ddroi1(kk),ddroi3(kk),ddte1(kk),
     &          ddte3(kk),ddti1(kk),ddti3(kk),xappae(kk),xappai(kk),
     &          pe1(kk),pe3(kk),pi1(kk),pi3(kk),te1(kk),te3(kk),ti1(kk),
     &          ti3(kk)
      COMMON /comth / ddroe1,ddroe3,ddroi1,ddroi3,ddte1,ddte3,ddti1,
     &                ddti3,gammae,gammai,xappae,xappai,pe1,pe3,pi1,pi3,
     &                te1,te3,teini,ti1,ti3,tiini
c end of thermodynamics
c c2.3 ions and electrons
      DIMENSION brems1(kk),brems3(kk),ddz1(kk),ddz3(kk),degen(kk),
     &          effz(kk),eixch2(kk),fidash(kk),fz1(kk),fz3(kk),fzsq1(kk)
     &          ,fzsq3(kk),xlc(kk),xmieff(kk),xne(kk),xni(kk),omega1(kk)
      COMMON /comie / brems1,brems3,ddz1,ddz3,degen,degmax,degmin,effz,
     &                eixch2,fidash,fz1,fz3,fzsq1,fzsq3,xlc,xmieff,xne,
     &                xni,omega1,pmass
c end of ions and electrons
c c3.2 mesh and numerical methods
      DIMENSION ae(kk),ai(kk),be(kk),bi(kk),ce(kk),ci(kk),de(kk),di(kk),
     &          e(kk),f(kk),ge(kk),gi(kk),q2(kk),q4(kk),teite(kk),
     &          tiite(kk),uite(kk)
      COMMON /comnum/ mesh,nj,njm1,njp1,nl,nlm1,nlp1,ae,ai,be,bi,ce,ci,
     &                de,di,e,f,ge,gi,q2,q4,teite,tiite,uite
c end of mesh and numerical methods
c end of physics control
      DATA iclass,isub/2,6/
      IF(nlomt2(isub))RETURN
c the factor  8255.17  is  k / m(proton)
      zg=8255.17/(gammai-1.0)
      IF(k.eq.1)THEN
      ELSEIF(k.eq.2)THEN
       GOTO 200
      ELSE
       CALL gotoer(' MEDUSA: statei: goto error, k is',k)
      ENDIF
c 1. derivatives of internal energy (eq.39)
      DO 100 l=1,nl
       ddti3(l)=zg/xmieff(l)
 100  CONTINUE
c the coefficient ddroi3 is zero
      RETURN
c 2. the ion pressure (eq.39)
 200  DO 300 l=1,nl
       pi3(l)=8255.17*ti3(l)*rho3(l)/xmieff(l)
 300  CONTINUE
      RETURN
      END


      SUBROUTINE stepon
      PARAMETER(kk=1001)
c 2.1 step on the calculation
c c1.1 basic system parameters
      LOGICAL nlend,nlres
      COMMON /combas/ altime,cptime,stime,ndiary,nin,nledge,nonlin,nout,
     &                nprint,npunch,nread,nrec,nresm,nrun,nstep,nlend,
     &                nlres
c end of basic system parameters
c c1.9 development and diagnostic parameters
      DIMENSION nadump(20),npdump(20),nvdump(20)
      LOGICAL nlched,nlhead(9),nlomt1(50),nlomt2(50),nlomt3(50),nlrept
      COMMON /comddp/ maxdum,mxdump,nadump,nclass,npdump,npoint,nsub,
     &                nvdump,nlched,nlhead,nlomt1,nlomt2,nlomt3,nlrept
c end of development and diagnostic parameters
c c2.1 hydrodynamics
      DIMENSION dm(kk),p3(kk),r1(kk),r3(kk),r5(kk),rho1(kk),rho3(kk),
     &          rho5(kk),u2(kk),u4(kk),v1(kk),v3(kk),v5(kk),lstat(kk)
      COMMON /comhyd/ dm,p3,pini,r1,r3,r5,rho1,rho3,rho5,rhoini,rhor,
     &                rini,time,u2,u4,uedge,v1,v3,v5
c end of hydrodynamics
c c2.7 physics control
      LOGICAL mlbrms,mlecon,mlfuse,mlicon,mlx,nlabs,nlbrms,nlburn,
     &        nlcri1,nldepo,nlecon,nlfuse,nlicon,nlmove,nlpfe,nlpfi,nlx,
     &        nlhf,nltnl
      COMMON /comcon/ mstep,ncase,ngeom,mlbrms,mlecon,mlfuse,mlicon,mlx,
     &                nlabs,nlbrms,nlburn,nlcri1,nldepo,nlecon,nlfuse,
     &                nlicon,nlmove,nlpfe,nlpfi,nlx,nlhf,nltnl,pmin,
     &                rhomin,temin,timin,umin
c end of physics control
c c3.1 numerical control parameters
      LOGICAL break,nlgoon,nlite
      COMMON /comnc / nceldt,ncondt,nit,nitmax,break,nlgoon,nlite,ak0,
     &                ak1,ak2,ak3,ak4,ak5,ak6,ak7,ak8,ak9,bneum,deltat,
     &                dt2,dt3,dt4,dtemax,dtimax,dtmax,dtmin,dumax,rdt2,
     &                rdt3,rdt4
c end of numerical control parameters
      DATA iclass,isub/2,1/
c time is now at level 1
c ====-==-===-==-=====-=
c update time and counters to level 3
      mstep=mstep+1
      nstep=nstep+1
      IF(nlomt2(isub))RETURN
      time=time+dt2
c advance one timestep
      CALL motion
c time is now at level 3
c ====-==-===-==-=====-=
c shift levels 5 to 3, 4 to 2 and 3 to 1
      CALL shifts
c the next step can now be implemented
      RETURN
      END


      SUBROUTINE stit
      PARAMETER(kk=1001)
c 2.5 starts the iteration by assuming level 3  = level1
c c1.9 development and diagnostic parameters
      DIMENSION nadump(20),npdump(20),nvdump(20)
      LOGICAL nlched,nlhead(9),nlomt1(50),nlomt2(50),nlomt3(50),nlrept
      COMMON /comddp/ maxdum,mxdump,nadump,nclass,npdump,npoint,nsub,
     &                nvdump,nlched,nlhead,nlomt1,nlomt2,nlomt3,nlrept
c end of development and diagnostic parameters
c c2.1 hydrodynamics
      DIMENSION dm(kk),p3(kk),r1(kk),r3(kk),r5(kk),rho1(kk),rho3(kk),
     &          rho5(kk),u2(kk),u4(kk),v1(kk),v3(kk),v5(kk),lstat(kk)
      COMMON /comhyd/ dm,p3,pini,r1,r3,r5,rho1,rho3,rho5,rhoini,rhor,
     &                rini,time,u2,u4,uedge,v1,v3,v5
c end of hydrodynamics
c c2.2 thermodynamics
      DIMENSION ddroe1(kk),ddroe3(kk),ddroi1(kk),ddroi3(kk),ddte1(kk),
     &          ddte3(kk),ddti1(kk),ddti3(kk),xappae(kk),xappai(kk),
     &          pe1(kk),pe3(kk),pi1(kk),pi3(kk),te1(kk),te3(kk),ti1(kk),
     &          ti3(kk)
      COMMON /comth / ddroe1,ddroe3,ddroi1,ddroi3,ddte1,ddte3,ddti1,
     &                ddti3,gammae,gammai,xappae,xappai,pe1,pe3,pi1,pi3,
     &                te1,te3,teini,ti1,ti3,tiini
c end of thermodynamics
c c3.1 numerical control parameters
      LOGICAL break,nlgoon,nlite
      COMMON /comnc / nceldt,ncondt,nit,nitmax,break,nlgoon,nlite,ak0,
     &                ak1,ak2,ak3,ak4,ak5,ak6,ak7,ak8,ak9,bneum,deltat,
     &                dt2,dt3,dt4,dtemax,dtimax,dtmax,dtmin,dumax,rdt2,
     &                rdt3,rdt4
c end of numerical control parameters
c c3.2 mesh and numerical methods
      DIMENSION ae(kk),ai(kk),be(kk),bi(kk),ce(kk),ci(kk),de(kk),di(kk),
     &          e(kk),f(kk),ge(kk),gi(kk),q2(kk),q4(kk),teite(kk),
     &          tiite(kk),uite(kk)
      COMMON /comnum/ mesh,nj,njm1,njp1,nl,nlm1,nlp1,ae,ai,be,bi,ce,ci,
     &                de,di,e,f,ge,gi,q2,q4,teite,tiite,uite
c end of mesh and numerical methods
      DATA iclass,isub/2,5/
      IF(nlomt2(isub))RETURN
c 1. reset quantities
      nit=0
      IF(.not.(nlite))RETURN
c 2. copy to -ite arrays for convergence check
      nlgoon=.false.
      break=.false.
      CALL copyr(u2,1,uite,1,njp1)
      CALL copyr(ti1,1,tiite,1,nlp1)
      CALL copyr(te1,1,teite,1,nlp1)
      RETURN
      END


      SUBROUTINE tesend
      PARAMETER(kk=1001)
c 4.1 test for completion of run
c c1.1 basic system parameters
      LOGICAL nlend,nlres
      COMMON /combas/ altime,cptime,stime,ndiary,nin,nledge,nonlin,nout,
     &                nprint,npunch,nread,nrec,nresm,nrun,nstep,nlend,
     &                nlres
c end of basic system parameters
c c2.1 hydrodynamics
      DIMENSION dm(kk),p3(kk),r1(kk),r3(kk),r5(kk),rho1(kk),rho3(kk),
     &          rho5(kk),u2(kk),u4(kk),v1(kk),v3(kk),v5(kk),lstat(kk)
      COMMON /comhyd/ dm,p3,pini,r1,r3,r5,rho1,rho3,rho5,rhoini,rhor,
     &                rini,time,u2,u4,uedge,v1,v3,v5
c end of hydrodynamics
c c4.1 - administrative variables
      LOGICAL nldump,nlemp
      DIMENSION piq(302)
      COMMON /comadm/ maxdim,maxrun,ndump,nrep,nldump,nlemp,piq,tstop
c c4.1 - end administrative variables
c the timestep limit
      IF(nstep.ge.nrun)nlend=.true.
c the time limit
      IF(time.ge.tstop)nlend=.true.
      RETURN
      END


      SUBROUTINE tfdu(t1,v1,z,cmult,dudt,dudv)
      DIMENSION a(6),b(7)
      DATA a/0.000000E+00,0.480750E+00,0.000000E+00,6.934000E-02,
     &     9.700000E-03,3.370400E-03/
      DATA b/0.000000E+00,0.480750E+00,0.434620E+00,6.920300E-02,
     &     5.947200E-02,-4.968800E-03,4.338600E-04/
      DATA cb1/4.486000E+10/
      DATA cb3/8.538300E+17/
      DATA cb4/1.381100E+21/
      DATA cb5/5.570700E+24/
      DATA cc1/3.980000E+23/
      DATA cc2/2.188000E+24/
      DATA ck1/2.110157E-11/
      DATA ck2/2.813543E-12/
      DATA ck3/3.283000E+07/
      DATA ck4/1.805000E+08/
      DATA cm1/2.550740E+10/
      DATA cm2/3.511074E+11/
      DATA cm3/4.365160E+10/
      DATA cm4/1.237000E-20/
      DATA r/1.324000E+08/
      pfgas=0.2
      zlog=log(z)
      z1=exp(0.333333333*zlog)
      z2=z*z1
      tdum=t1
      IF(t1.lt.100.0)tdum=100.0
      t=tdum/(1.161E4*z2)
      v=1.0E6*z*v1
      tlog=log(t)
      v1log=log(v1)
      v2=sqrt(v)
      v3=1.0E2*exp(0.333333333*(v1log+zlog))
      v6=sqrt(v3)
      dudt=pfgas*2.07E-23*z
c energy derived by integration w.r.t. v of t * dp / dt - p
c at non-zero temp
      c1=ck3*exp(tlog*0.1466666666)
      c2=ck4/exp(tlog*1.0733333333)
      c8=c1+c2
      c9=c8*sqrt(c8)
      v24=v*1.0E24
      du1dt=(1.61*c2-0.22*c1)*(1.0/sqrt(v24*sqrt(v24))+log(v*1.0E22))
     &      /(t*c8*c9)
      du1dv=(1.0-0.75/sqrt(v24*sqrt(v24)))/(c9*v)
c derivative of energy as a function of t along line v = 1.0e-20 cc
      du0dt=2.403E-12*t*(2.2+t)/((1.1+t)*(1.1+t))
c sum components of derivatives
      dudt=dudt+(1.0-pfgas)*(du0dt+du1dt)*8.6133E-12*z
      dudv=(1.0-pfgas)*du1dv*0.1*z*z*z2
c add pressure to dudv
      p2=cc1/(t*t)+cc2/exp(tlog*3.22)
      v24=v*1.0E24
      p2=1.0/(v*sqrt(p2))-(1.0/sqrt(p2)-t/sqrt(cc1))
     &   *0.75/(v*sqrt(v24*sqrt(v24)))
      tfp=(1.0-pfgas)*0.1*p2*z*z*z2
      tfp=tfp+pfgas*1.38E-23*t1*z/v1
      dudv=tfp+dudv
      RETURN
      END
      FUNCTION tfp(t1,v1,z,cmult)
      DATA cb1/4.4860E10/
      DATA cb3/8.5383E17/
      DATA cb4/1.3811E21/
      DATA cb5/5.5707E24/
      DATA cc1/3.9800E23/
      DATA cc2/2.1880E24/
      pfgas=0.2
      zlog=log(z)
      z1=exp(0.333333333*zlog)
      z2=z*z1
      tdum=t1
      IF(t1.lt.100.0)tdum=100.0
      t=tdum/(1.161E4*z2)
      v=1.0E6*z*v1
      tlog=log(t)
      v1log=log(v1)
      v2=sqrt(v)
      v3=1.0E2*exp(0.333333333*(v1log+zlog))
      v6=sqrt(v3)
      tfp=pfgas*1.38E-23*t1*z/v1
c degen curve from 'march, adv phys, 6, 1'
      p1=cb1*v2*v6+v*(cb3+v6*(cb4+v6*cb5))
      p1=1.0/(p1*p1*sqrt(p1))
c perfect gas pressure with ionisation
      p2=cc1/(t*t)+cc2/exp(tlog*3.22)
      v24=v*1.0E24
      p2=1.0/(v*sqrt(p2))-(1.0/sqrt(p2)-t/sqrt(cc1))
     &   *0.75/(v*sqrt(v24*sqrt(v24)))
c sum pressures
      tfp=tfp+0.1*(p1+(1.0-pfgas)*p2)*z*z*z2
      IF(cmult.lt.-0.5)RETURN
c tfc correction term giving binding
c zero temp part of correection
      p1=1.0229279E14*exp(0.6645*v1log-0.66883335*zlog)
      p2=3.1077321E28*exp(1.185*v1log-0.14833333*zlog)
      dp=1.0/(p1+p2)
      dp=dp*dp
      tfp=tfp-dp*cmult
      RETURN
      END
      FUNCTION tfu(t1,v1,z,cmult)
      DIMENSION b(7)
      DIMENSION a(6)
      DATA a/0.00000,0.48075,0.00000,0.06934,9.70000E-3,3.37040E-3/
      DATA b/0.00000,0.48075,0.43462,6.92030E-2,5.94720E-2,-4.96880E-3,
     &     4.33860E-4/
      pfgas=0.2
      zlog=log(z)
      v=1.0E6*v1*z
      vlog=log(v)
      tdum=t1
      t1log=log(tdum)
      tlog=t1log-1.33333333*zlog-9.3596221
      t=exp(tlog)
      tfu=pfgas*2.07E-23*t1*z/v1
      upf=tfu*v1/1.602E-19
c degeneracy energy from 'march, adv phys, 6, 1'
      x0log=0.333333333*(vlog+56.104177)
      x0=exp(x0log)
      sum=1.5311E-6*exp(7.772*x0log)
      x1=x0
      DO 100 i=1,7
       sum=sum+b(i)*x1
       x1=x1*x0
 100  CONTINUE
      ad=1.0/sum
      sum=0.0
      x1=sqrt(x0)
      x2=1.0
      DO 200 i=1,6
       x2=x2*x1
       sum=sum+a(i)*x2
 200  CONTINUE
      phi=1.0/sum
      ud=2.1101E-11*ad+2.8135E-12*sqrt(phi*x0)*phi*phi
      ud=ud/1.602E-12
      t0=1.0
c energy derived by integration w.r.t. v of
c t * dp / dt - p at non-zero temp
      c1=3.283E7*exp(0.14666666*tlog)
      c2=1.805E8/exp(1.07333333*tlog)
      c3=c1+c2
      v24=v*1.0E24
      u1=(1.0/sqrt(v24*sqrt(v24))+log(v*1.0E22))/(c3*sqrt(c3))
      u1=u1/1.602E-12
c energy along line v = 1.0e-20 cc which arose as constant of
c integration w.r.t. v of t * dp / dp - p
      u0=1.5*t*t/(1.1+t)
c sum the component energies
      u=ud+(u1+u0)*(1.0-pfgas)
      tfu=tfu+(u*1.602E-19*exp(2.3333333*zlog))/v1
      IF(cmult.lt.-0.5)RETURN
c tfc energy corections
c degen curve correction
      v2=v1*z
      uct0=-uc0(v2)*exp(1.666666*zlog)
      uct=uct0
      uctt=uct/1.602E-19
      tfu=tfu+cmult*(uct)/v1
      RETURN
      END


      FUNCTION uc0(v)
      DIMENSION udc(101)
      DATA(udc(i),i=1,40)/24.2417,24.1071,23.9723,23.8373,23.7021,
     &     23.5669,23.4317,23.2966,23.1617,23.0269,22.8923,22.7580,
     &     22.6240,22.4904,22.3571,22.2243,22.0920,21.9603,21.8291,
     &     21.6986,21.5689,21.4401,21.3121,21.1852,21.0594,20.9348,
     &     20.8115,20.6896,20.5694,20.4507,20.3339,20.2190,20.1060,
     &     19.9953,19.8867,19.7806,19.6768,19.5755,19.4768,19.3808/
      DATA(udc(i),i=41,80)/19.2874,19.1966,19.1085,19.0231,18.9403,
     &     18.8601,18.7824,18.7073,18.6345,18.5640,18.4957,18.4296,
     &     18.3656,18.3034,18.2432,18.1846,18.1278,18.0725,18.0186,
     &     17.9662,17.9150,17.8651,17.8163,17.7686,17.7219,17.6761,
     &     17.6312,17.5871,17.5438,17.5012,17.4592,17.4179,17.3771,
     &     17.3369,17.2972,17.2580,17.2192,17.1808,17.1428,17.1051/
      DATA(udc(i),i=81,101)/17.0678,17.0307,16.9940,16.9575,16.9213,
     &     16.8853,16.8495,16.8139,16.7785,16.7433,16.7083,16.6734,
     &     16.6386,16.6040,16.5694,16.5350,16.5007,16.4666,16.4324,
     &     16.3984,16.3645/
      vlog=log10(v)
      i1=nint((vlog+23.9)*10.0)
      i1=-i1
      IF(i1.le.0)i1=1
      IF(i1.gt.100)i1=100
      i2=i1+1
      dvlog1=(vlog+23.9+i1*0.1)*10.0
      dvlog1=-dvlog1
      dvlog2=1.0-dvlog1
      udclog=dvlog1*udc(i2)+dvlog2*udc(i1)
      uc0=0.1**udclog
      RETURN
      END


      SUBROUTINE user(k)
      RETURN
      END


      SUBROUTINE volume
      PARAMETER(kk=1001,underf=1.E-30)
c 2.20 volumes and densities
c c1.9 development and diagnostic parameters
      DIMENSION nadump(20),npdump(20),nvdump(20)
      LOGICAL nlched,nlhead(9),nlomt1(50),nlomt2(50),nlomt3(50),nlrept
      COMMON /comddp/ maxdum,mxdump,nadump,nclass,npdump,npoint,nsub,
     &                nvdump,nlched,nlhead,nlomt1,nlomt2,nlomt3,nlrept
c end of development and diagnostic parameters
c c2.1 hydrodynamics
      DIMENSION dm(kk),p3(kk),r1(kk),r3(kk),r5(kk),rho1(kk),rho3(kk),
     &          rho5(kk),u2(kk),u4(kk),v1(kk),v3(kk),v5(kk),lstat(kk)
      COMMON /comhyd/ dm,p3,pini,r1,r3,r5,rho1,rho3,rho5,rhoini,rhor,
     &                rini,time,u2,u4,uedge,v1,v3,v5
c end of hydrodynamics
c c2.3 ions and electrons
      DIMENSION brems1(kk),brems3(kk),ddz1(kk),ddz3(kk),degen(kk),
     &          effz(kk),eixch2(kk),fidash(kk),fz1(kk),fz3(kk),fzsq1(kk)
     &          ,fzsq3(kk),xlc(kk),xmieff(kk),xne(kk),xni(kk),omega1(kk)
      COMMON /comie / brems1,brems3,ddz1,ddz3,degen,degmax,degmin,effz,
     &                eixch2,fidash,fz1,fz3,fzsq1,fzsq3,xlc,xmieff,xne,
     &                xni,omega1,pmass
c end of ions and electrons
c c2.5 thermonuclear reactions
      DIMENSION f1d(kk),f1h(kk),f1he3(kk),f1he4(kk),f1neu(kk),f1ntrl(kk)
     &          ,f1t(kk),f1x(kk),f1x1(kk),f1x2(kk),f1x3(kk),f1x4(kk),
     &          f3d(kk),f3h(kk),f3he3(kk),f3he4(kk),f3neu(kk),f3ntrl(kk)
     &          ,f3t(kk),f3x(kk),f3x1(kk),f3x2(kk),f3x3(kk),f3x4(kk),
     &          r1dd(kk),r1dhe3(kk),r1dt(kk),r3dd(kk),r3dhe3(kk),
     &          r3dt(kk),ye1(kk),ye3(kk),yi1(kk),yi3(kk)
      COMMON /comfus/ deuter,f1d,f1h,f1he3,f1he4,f1neu,f1ntrl,f1t,f1x,
     &                f1x1,f1x2,f1x3,f1x4,f3d,f3h,f3he3,f3he4,f3neu,
     &                f3ntrl,f3t,f3x,f3x1,f3x2,f3x3,f3x4,heliu3,heliu4,
     &                hydrog,xetral,xtrlms,pneut1,pneut3,r1dd,r1dhe3,
     &                r1dt,r3dd,r3dhe3,r3dt,rneut1,rneut3,tinucl,totneu,
     &                tritiu,xmass,xmass1,xmass2,xmass3,xmass4,xtra,
     &                xtra1,xtra2,xtra3,xtra4,xz,xz1,xz2,xz3,xz4,ye1,
     &                ye3,yi1,yi3,yield
c end of thermonuclear reactions
c c2.7 physics control
      LOGICAL mlbrms,mlecon,mlfuse,mlicon,mlx,nlabs,nlbrms,nlburn,
     &        nlcri1,nldepo,nlecon,nlfuse,nlicon,nlmove,nlpfe,nlpfi,nlx,
     &        nlhf,nltnl
      COMMON /comcon/ mstep,ncase,ngeom,mlbrms,mlecon,mlfuse,mlicon,mlx,
     &                nlabs,nlbrms,nlburn,nlcri1,nldepo,nlecon,nlfuse,
     &                nlicon,nlmove,nlpfe,nlpfi,nlx,nlhf,nltnl,pmin,
     &                rhomin,temin,timin,umin
c end of physics control
c c3.1 numerical control parameters
      LOGICAL break,nlgoon,nlite
      COMMON /comnc / nceldt,ncondt,nit,nitmax,break,nlgoon,nlite,ak0,
     &                ak1,ak2,ak3,ak4,ak5,ak6,ak7,ak8,ak9,bneum,deltat,
     &                dt2,dt3,dt4,dtemax,dtimax,dtmax,dtmin,dumax,rdt2,
     &                rdt3,rdt4
c end of numerical control parameters
c c3.2 mesh and numerical methods
      DIMENSION ae(kk),ai(kk),be(kk),bi(kk),ce(kk),ci(kk),de(kk),di(kk),
     &          e(kk),f(kk),ge(kk),gi(kk),q2(kk),q4(kk),teite(kk),
     &          tiite(kk),uite(kk)
      COMMON /comnum/ mesh,nj,njm1,njp1,nl,nlm1,nlp1,ae,ai,be,bi,ce,ci,
     &                de,di,e,f,ge,gi,q2,q4,teite,tiite,uite
c end of mesh and numerical methods
      DATA iclass,isub/2,20/
      IF(nlomt2(isub))RETURN
c 1. mass and charge numbers. particle densities.
      zdm2=0.0
      zdni2=0.0
c standard ion masses
      zmd=2.0
      zmh=1.0
      zmhe3=3.0
      zmhe4=4.0
      zmt=3.0
c standard charge numbers
      zcd=1.0
      zch=1.0
      zche3=2.0
      zche4=2.0
      zct=1.0
c these quantities are at level 3
      DO 100 l=1,nl
       j=l
c eq.3
       xmieff(l)=zmd*f3d(l)+zmh*f3h(l)+zmhe3*f3he3(l)+zmhe4*f3he4(l)
     &           +xtrlms*f3ntrl(l)+zmt*f3t(l)+xmass*f3x(l)
     &           +xmass1*f3x1(l)+xmass2*f3x2(l)+xmass3*f3x3(l)
     &           +xmass4*f3x4(l)
c eq.4
c we only wish to recompute fz3 if there is tritium in the cell
       IF(abs(f3t(l)).gt.underf)fz3(l)=zcd*f3d(l)+zch*f3h(l)
     &                                 +zche3*f3he3(l)+zche4*f3he4(l)
     &                                 +zct*f3t(l)+xz*f3x(l)+xz1*f3x1(l)
     &                                 +xz2*f3x2(l)+xz3*f3x3(l)
     &                                 +xz4*f3x4(l)
c find total volume of a cell at level 3
       IF(ngeom.eq.1)THEN
       ELSEIF(ngeom.eq.2)THEN
        zvol=0.5*(r3(j+1)*r3(j+1)-r3(j)*r3(j))
        GOTO 50
       ELSEIF(ngeom.eq.3)THEN
        zvol=(r3(j+1)*r3(j+1)*r3(j+1)-r3(j)*r3(j)*r3(j))/3.0
        GOTO 50
       ELSE
        CALL gotoer(' MEDUSA: volume: goto error, ngeom is',ngeom)
       ENDIF
       zvol=r3(j+1)-r3(j)
c define mass number at level 1
 50    zmi1=zmd*f1d(l)+zmh*f1h(l)+zmhe3*f1he3(l)+zmhe4*f1he4(l)
     &      +xtrlms*f1ntrl(l)+zmt*f1t(l)+xmass*f1x(l)+xmass1*f1x1(l)
     &      +xmass2*f1x2(l)+xmass3*f1x3(l)+xmass4*f1x4(l)
c if no burning takes place there are no changes in ni and dm due
c to the neutron loss
       IF(nlburn)THEN
c the change in ni due to d-t and 1 of every 2 d-d reactions
        zdni2=-0.5*dt2*(r3dt(l)+r1dt(l)+0.5*(r3dd(l)+r1dd(l)))
c the change in mass equals the number of neutrons lost times pmass
        zdm2=zdni2*zvol*pmass
       ENDIF
c eq.6
       xni(l)=rho3(l)/(pmass*zmi1)+zdni2
c eq.5
       xne(l)=fz3(l)*xni(l)
c the new mass at level 3
       dm(l)=zvol*rho3(l)+zdm2
 100  CONTINUE
c 2. specific volumes at level 5 (eq.48)
      IF(ngeom.eq.1)THEN
      ELSEIF(ngeom.eq.2)THEN
c 2.2 cylinder geometry
c we are treating a cylindrical section: angle 1, height 1
       DO 150 l=1,nl
        j=l
        v5(l)=(r5(j+1)*r5(j+1)-r5(j)*r5(j))/(dm(l)*2.0)
        rho5(l)=1.0/v5(l)
 150   CONTINUE
       RETURN
      ELSEIF(ngeom.eq.3)THEN
       GOTO 300
      ELSE
       CALL gotoer(' MEDUSA: volume: goto error, ngeom is',ngeom)
      ENDIF
c the specific volumes are at level 5, but the mass is at level 3
c 2.1 slab geometry
c we are treating a slab of unit cross section
      DO 200 l=1,nl
       j=l
       v5(l)=(r5(j+1)-r5(j))/dm(l)
       rho5(l)=1.0/v5(l)
 200  CONTINUE
      RETURN
c 2.3 spherical geometry
c we are considering one steradian
 300  DO 400 l=1,nl
       j=l
       v5(l)=(r5(j+1)*r5(j+1)*r5(j+1)-r5(j)*r5(j)*r5(j))/(3.0*dm(l))
       rho5(l)=1.0/v5(l)
 400  CONTINUE
      RETURN
      END


      SUBROUTINE xchang(k)
      PARAMETER(kk=1001)
c 2.11 the rate of exchange of energy
c c1.9 development and diagnostic parameters
      DIMENSION nadump(20),npdump(20),nvdump(20)
      LOGICAL nlched,nlhead(9),nlomt1(50),nlomt2(50),nlomt3(50),nlrept
      COMMON /comddp/ maxdum,mxdump,nadump,nclass,npdump,npoint,nsub,
     &                nvdump,nlched,nlhead,nlomt1,nlomt2,nlomt3,nlrept
c end of development and diagnostic parameters
c c2.1 hydrodynamics
c end of hydrodynamics
c c2.2 thermodynamics
      DIMENSION ddroe1(kk),ddroe3(kk),ddroi1(kk),ddroi3(kk),ddte1(kk),
     &          ddte3(kk),ddti1(kk),ddti3(kk),xappae(kk),xappai(kk),
     &          pe1(kk),pe3(kk),pi1(kk),pi3(kk),te1(kk),te3(kk),ti1(kk),
     &          ti3(kk)
      COMMON /comth / ddroe1,ddroe3,ddroi1,ddroi3,ddte1,ddte3,ddti1,
     &                ddti3,gammae,gammai,xappae,xappai,pe1,pe3,pi1,pi3,
     &                te1,te3,teini,ti1,ti3,tiini
c end of thermodynamics
c c2.3 ions and electrons
      DIMENSION brems1(kk),brems3(kk),ddz1(kk),ddz3(kk),degen(kk),
     &          effz(kk),eixch2(kk),fidash(kk),fz1(kk),fz3(kk),fzsq1(kk)
     &          ,fzsq3(kk),xlc(kk),xmieff(kk),xne(kk),xni(kk),omega1(kk)
      COMMON /comie / brems1,brems3,ddz1,ddz3,degen,degmax,degmin,effz,
     &                eixch2,fidash,fz1,fz3,fzsq1,fzsq3,xlc,xmieff,xne,
     &                xni,omega1,pmass
c end of ions and electrons
c c2.5 thermonuclear reactions
      DIMENSION f1d(kk),f1h(kk),f1he3(kk),f1he4(kk),f1neu(kk),f1ntrl(kk)
     &          ,f1t(kk),f1x(kk),f1x1(kk),f1x2(kk),f1x3(kk),f1x4(kk),
     &          f3d(kk),f3h(kk),f3he3(kk),f3he4(kk),f3neu(kk),f3ntrl(kk)
     &          ,f3t(kk),f3x(kk),f3x1(kk),f3x2(kk),f3x3(kk),f3x4(kk),
     &          r1dd(kk),r1dhe3(kk),r1dt(kk),r3dd(kk),r3dhe3(kk),
     &          r3dt(kk),ye1(kk),ye3(kk),yi1(kk),yi3(kk)
      COMMON /comfus/ deuter,f1d,f1h,f1he3,f1he4,f1neu,f1ntrl,f1t,f1x,
     &                f1x1,f1x2,f1x3,f1x4,f3d,f3h,f3he3,f3he4,f3neu,
     &                f3ntrl,f3t,f3x,f3x1,f3x2,f3x3,f3x4,heliu3,heliu4,
     &                hydrog,xetral,xtrlms,pneut1,pneut3,r1dd,r1dhe3,
     &                r1dt,r3dd,r3dhe3,r3dt,rneut1,rneut3,tinucl,totneu,
     &                tritiu,xmass,xmass1,xmass2,xmass3,xmass4,xtra,
     &                xtra1,xtra2,xtra3,xtra4,xz,xz1,xz2,xz3,xz4,ye1,
     &                ye3,yi1,yi3,yield
c end of thermonuclear reactions
c c2.7 physics control
      LOGICAL mlbrms,mlecon,mlfuse,mlicon,mlx,nlabs,nlbrms,nlburn,
     &        nlcri1,nldepo,nlecon,nlfuse,nlicon,nlmove,nlpfe,nlpfi,nlx,
     &        nlhf,nltnl
      COMMON /comcon/ mstep,ncase,ngeom,mlbrms,mlecon,mlfuse,mlicon,mlx,
     &                nlabs,nlbrms,nlburn,nlcri1,nldepo,nlecon,nlfuse,
     &                nlicon,nlmove,nlpfe,nlpfi,nlx,nlhf,nltnl,pmin,
     &                rhomin,temin,timin,umin
c end of physics control
c c3.1 numerical control parameters
      LOGICAL break,nlgoon,nlite
      COMMON /comnc / nceldt,ncondt,nit,nitmax,break,nlgoon,nlite,ak0,
     &                ak1,ak2,ak3,ak4,ak5,ak6,ak7,ak8,ak9,bneum,deltat,
     &                dt2,dt3,dt4,dtemax,dtimax,dtmax,dtmin,dumax,rdt2,
     &                rdt3,rdt4
c end of numerical control parameters
c c3.2 mesh and numerical methods
      DIMENSION ae(kk),ai(kk),be(kk),bi(kk),ce(kk),ci(kk),de(kk),di(kk),
     &          e(kk),f(kk),ge(kk),gi(kk),q2(kk),q4(kk),teite(kk),
     &          tiite(kk),uite(kk)
      COMMON /comnum/ mesh,nj,njm1,njp1,nl,nlm1,nlp1,ae,ai,be,bi,ce,ci,
     &                de,di,e,f,ge,gi,q2,q4,teite,tiite,uite
c end of mesh and numerical methods
c c4.1 - administrative variables
      LOGICAL nldump,nlemp
      DIMENSION piq(302)
      COMMON /comadm/ maxdim,maxrun,ndump,nrep,nldump,nlemp,piq,tstop
c c4.1 - end administrative variables

      COMMON /radflg/ yamda,xhnu,xihnu,dihnu,xplas1(kk),xrbf(kk),
     &                atrbf(kk),corhf1(kk),corhf(kk),iradf1,iradf2

      DATA iclass,isub/2,11/
      IF(nlomt2(isub))RETURN
c 1. set constants
c set the constant of equation 22
c     zc        = 0.5936e-8 * (1.0 + piq(25))
      zc=0.4E-8*(1.0+piq(25))
      zg=gammai-1.0
c for perfect ions-electrons zcvie = 1 + ddte / ddti=2
      zcvie=2.0
c define products of z ** 2 and m ** (-1)
      zzmd=1.0/2.0
      zzmh=1.0/1.0
      zzmhe3=4.0/3.0
      zzmhe4=4.0/4.0
      zzmt=1.0/3.0
      zzmx=xz*xz/xmass
c are we calculating omega or the rate of energy?
      IF(k.eq.1)THEN
      ELSEIF(k.eq.2)THEN
       GOTO 200
      ELSE
       CALL gotoer(' MEDUSA: xchang: goto error, k is',k)
      ENDIF
c 2.  find omega at level 1 from eq.20
      DO 100 l=1,nl
c calculate the average of z ** 2 and m ** (-1) at level 1
c        zzm            = fzsq1(l) / (xmieff(l) * xmieff(l))
       zzm=fzsq1(l)/xmieff(l)
       omega1(l)=zc*xlc(l)*xni(l)*zzm/(te1(l)*sqrt(te1(l)))
c ad1 allan
       IF(nlhf)omega1(l)=omega1(l)*corhf1(l)
 100  CONTINUE
      RETURN
c 3. the rate of exchanged energy at level 2
 200  DO 300 l=1,nl
c       evaluate zcvie for non-perfect electrons
c       if ( .not.nlpfe )
       zcvie=1.0+(ddti1(l)+ddti3(l))/(ddte1(l)+ddte3(l))
c       calculate the average of z ** 2 and m ** (-1) at level 3
       zzm=fzsq3(l)/xmieff(l)
c       the time constant at level 3 (eq.20)
       zomeg3=zc*xlc(l)*xni(l)*zzm/(te3(l)*sqrt(te3(l)))
c ad1 allan
       IF(nlhf)zomeg3=zomeg3*corhf(l)
c       the average time constant at level 2 (eq.55)
       zomeg2=0.5*(omega1(l)+zomeg3)
c       zcg is a factor converting to rate of energy per unit mass
c       zcg    = 8255.17 / (zg * xmieff(l))
       zcg=(fz3(l)*8255.17)/(zg*xmieff(l))
c       the difference between the ion and electron source terms
c       (specific heats are not necessarily equal)
c       take into account the factor already contained in eixch2 (eq.58)
       zsdif2=rdt2*(ti3(l)-ti1(l)-te3(l)+te1(l))+eixch2(l)
     &        /zcg*(zcvie/2.0)
c       the relaxation factor (numerator of eq.57)
       zomegt=-zcvie*zomeg2*dt2
       IF(abs(zomegt).gt.10.0)zomegt=-10.0
       IF(abs(zomegt).lt.1.0E-03)THEN
c         3.2 omega * dt small
c         we expand the exponentiel term and approximate for ztie2 eq.56
        ztie2=(ti1(l)-te1(l))*zomeg2
       ELSE
        zrelax=1.0-exp(zomegt)
c         3.1 omega * dt large
c         the correct averaged product of ti2 - te2 and omega (eq.56)
        ztie2=(ti1(l)-te1(l)-zsdif2/(zcvie*zomeg2))*zrelax/(zcvie*dt2)
     &        +zsdif2/zcvie
       ENDIF
c       we include a  2.0  in the exchange term because of the averaging
c       made in evaluating the coefficients a, b, c and d
       eixch2(l)=2.0*ztie2*zcg

c      ad1
c      check rate of energy exch does not exceed maximum
       tteq=(ddte1(l)+ddte3(l))*te1(l)+(ddti1(l)+ddti3(l))
     &      *ti1(l)
       tteq=tteq/(ddte1(l)+ddte3(l)+ddti1(l)+ddti3(l))
       eimax=0.5*(ddte1(l)+ddte3(l))*abs(te1(l)-tteq)
       IF(abs(eixch2(l)*dt2).gt.eimax)THEN
        efac=eimax/abs(eixch2(l)*dt2)
        eixch2(l)=eixch2(l)*efac
       ENDIF

 300  CONTINUE
      RETURN
      END


      SUBROUTINE data
      PARAMETER(kk=1001)
c 1.4 define data specific to run
c c1.1 basic system parameters
      LOGICAL nlend,nlres
      COMMON /combas/ altime,cptime,stime,ndiary,nin,nledge,nonlin,nout,
     &                nprint,npunch,nread,nrec,nresm,nrun,nstep,nlend,
     &                nlres
c end of basic system parameters
c c1.9 development and diagnostic parameters
      DIMENSION nadump(20),npdump(20),nvdump(20)
      LOGICAL nlched,nlhead(9),nlomt1(50),nlomt2(50),nlomt3(50),nlrept
      COMMON /comddp/ maxdum,mxdump,nadump,nclass,npdump,npoint,nsub,
     &                nvdump,nlched,nlhead,nlomt1,nlomt2,nlomt3,nlrept
c end of development and diagnostic parameters
c c2.1 hydrodynamics
      DIMENSION dm(kk),p3(kk),r1(kk),r3(kk),r5(kk),rho1(kk),rho3(kk),
     &          rho5(kk),u2(kk),u4(kk),v1(kk),v3(kk),v5(kk),lstat(kk)
      COMMON /comhyd/ dm,p3,pini,r1,r3,r5,rho1,rho3,rho5,rhoini,rhor,
     &                rini,time,u2,u4,uedge,v1,v3,v5
c end of hydrodynamics
c c2.2 thermodynamics
      REAL ddroe1(kk),ddroe3(kk),ddroi1(kk),ddroi3(kk),ddte1(kk),
     &     ddte3(kk),ddti1(kk),ddti3(kk),gammae,gammai,xappae(kk),
     &     xappai(kk),pe1(kk),pe3(kk),pi1(kk),pi3(kk),te1(kk),te3(kk),
     &     teini,ti1(kk),ti3(kk),tiini
      COMMON /comth / ddroe1,ddroe3,ddroi1,ddroi3,ddte1,ddte3,ddti1,
     &                ddti3,gammae,gammai,xappae,xappai,pe1,pe3,pi1,pi3,
     &                te1,te3,teini,ti1,ti3,tiini
c end of thermodynamics
c c2.4 laser variables
      DIMENSION alpha1(kk),xaser1(kk),xl1(kk),xl3(kk)
      COMMON /comlas/ alpha1,elas1,xamda1,xaser1,xecri1,plas1,rabs1,
     &                rocri1,xl1,xl3,nabs1
c end of laser variables
c c2.5 thermonuclear reactions
      REAL deuter,f1d(kk),f1h(kk),f1he3(kk),f1he4(kk),f1neu(kk),
     &     f1ntrl(kk),f1t(kk),f1x(kk),f1x1(kk),f1x2(kk),f1x3(kk),
     &     f1x4(kk),f3d(kk),f3h(kk),f3he3(kk),f3he4(kk),f3neu(kk),
     &     f3ntrl(kk),f3t(kk),f3x(kk),f3x1(kk),f3x2(kk),f3x3(kk),
     &     f3x4(kk),heliu3,heliu4,hydrog,xetral,xtrlms,pneut1,pneut3,
     &     r1dd(kk),r1dhe3(kk),r1dt(kk),r3dd(kk),r3dhe3(kk),r3dt(kk),
     &     rneut1,rneut3,tinucl,totneu,tritiu,xmass,xmass1,xmass2,
     &     xmass3,xmass4,xtra,xtra1,xtra2,xtra3,xtra4,xz,xz1,xz2,xz3,
     &     xz4,ye1(kk),ye3(kk),yi1(kk),yi3(kk),yield
      COMMON /comfus/ deuter,f1d,f1h,f1he3,f1he4,f1neu,f1ntrl,f1t,f1x,
     &                f1x1,f1x2,f1x3,f1x4,f3d,f3h,f3he3,f3he4,f3neu,
     &                f3ntrl,f3t,f3x,f3x1,f3x2,f3x3,f3x4,heliu3,heliu4,
     &                hydrog,xetral,xtrlms,pneut1,pneut3,r1dd,r1dhe3,
     &                r1dt,r3dd,r3dhe3,r3dt,rneut1,rneut3,tinucl,totneu,
     &                tritiu,xmass,xmass1,xmass2,xmass3,xmass4,xtra,
     &                xtra1,xtra2,xtra3,xtra4,xz,xz1,xz2,xz3,xz4,ye1,
     &                ye3,yi1,yi3,yield
c end of thermonuclear reactions
c c2.7 physics control
      INTEGER mstep,ncase,ngeom
      LOGICAL mlbrms,mlecon,mlfuse,mlicon,mlx,nlabs,nlbrms,nlburn,
     &        nlcri1,nldepo,nlecon,nlfuse,nlicon,nlmove,nlpfe,nlpfi,nlx,
     &        nlhf,nltnl
      REAL pmin,rhomin,temin,timin,umin
      COMMON /comcon/ mstep,ncase,ngeom,mlbrms,mlecon,mlfuse,mlicon,mlx,
     &                nlabs,nlbrms,nlburn,nlcri1,nldepo,nlecon,nlfuse,
     &                nlicon,nlmove,nlpfe,nlpfi,nlx,nlhf,nltnl,pmin,
     &                rhomin,temin,timin,umin
c end of physics control
c c3.1 numerical control parameters
      INTEGER nceldt,ncondt,nit,nitmax
      LOGICAL break,nlgoon,nlite
      REAL ak0,ak1,ak2,ak3,ak4,ak5,ak6,ak7,ak8,ak9,bneum,deltat,dt2,dt3,
     &     dt4,dtemax,dtimax,dtmax,dtmin,dumax,rdt2,rdt3,rdt4
      COMMON /comnc / nceldt,ncondt,nit,nitmax,break,nlgoon,nlite,ak0,
     &                ak1,ak2,ak3,ak4,ak5,ak6,ak7,ak8,ak9,bneum,deltat,
     &                dt2,dt3,dt4,dtemax,dtimax,dtmax,dtmin,dumax,rdt2,
     &                rdt3,rdt4
c end of numerical control parameters
      LOGICAL lphi(kk)
      DIMENSION deltae(50),ebdy(51),emid(50),ha(kk),hb(kk),hc(kk),hd(kk)
     &          ,hdiff(kk),he(kk),heta(kk),hf(kk),hg(kk),hgrp(kk),hj(kk)
     &          ,hlas(kk),hloss(50),hncold(kk),hp(kk),hphi(kk),htaug(kk)
     &          ,htherm(kk),htr(kk),xhot1(kk,50),xhot3(kk,50)
      COMMON /comhot/ igrp,ngroup,lphi,deltae,ebdy,ehot,emid,ha,hb,hc,
     &                hd,hdiff,he,heta,hf,hg,hgrp,hj,hlas,hloss,hncold,
     &                hp,hphi,htaug,htherm,htr,xhot1,xhot3,thot
c c3.2 mesh and numerical methods
      DIMENSION ae(kk),ai(kk),be(kk),bi(kk),ce(kk),ci(kk),de(kk),di(kk),
     &          e(kk),f(kk),ge(kk),gi(kk),q2(kk),q4(kk),teite(kk),
     &          tiite(kk),uite(kk)
      COMMON /comnum/ mesh,nj,njm1,njp1,nl,nlm1,nlp1,ae,ai,be,bi,ce,ci,
     &                de,di,e,f,ge,gi,q2,q4,teite,tiite,uite
c end of mesh and numerical methods
c c4.1 - administrative variables
      LOGICAL nldump,nlemp
      DIMENSION piq(302)
      COMMON /comadm/ maxdim,maxrun,ndump,nrep,nldump,nlemp,piq,tstop
c c4.1 - end administrative variables
c c5.1 input-output control variables
      INTEGER nfilm,nhdcpy,np1,np2,np3
      LOGICAL nlfilm,nlhcpy,nlprnt
      REAL buf1(kk),buf2(kk),buf3(kk),buf4(kk),buf5(kk),buf6(kk),scp,
     &     scr,scrho,scte,scti,sctime
      REAL nprnt,nproc
      COMMON /comout/ nfilm,nhdcpy,np1,np2,np3,nlfilm,nlhcpy,nlprnt,
     &        ztnext,nprnt,nproc,buf1,buf2,buf3,buf4,buf5,buf6,scp,
     &                scr,scrho,scte,scti,sctime,nt1,nt2
c end of input-output control variables
c x-ray laser common blocks
      DIMENSION fl(2)
      COMMON /xrl   / ropmul,tpon,tpoff,rmpd,ztrtop,dlamda,fl,fwide,
     &                efbmul,drmul,bbtrap,fbtrap,icxrl,istage,idflg,
     &                itbflg,istflg,idrflg,ipuflg,nlp,nup,llp,lup,nst,
     &                nfl,nlmax,ifrsta,ilosta,ihista,igstat
c end of input-output control variables
c   new subroutine source control factors and layer materials
      COMMON /sornew/ rf1,rf2,rf3,fne2
      EQUIVALENCE(piq(63),zplas),(piq(62),roplas),(piq(61),drplas),
     &            (piq(59),state),(piq(58),saha),(piq(56),pondf),
     &            (piq(57),rhot),(piq(37),gauss),(piq(52),anabs),
     &            (piq(22),fthot),(piq(21),fhot),(piq(14),fne),
     &            (piq(13),zglas),(piq(12),roglas),(piq(11),drglas),
     &            (piq(10),flimit),(rhoini,rhogas),(xaser1(1),ton),
     &            (xaser1(2),pmult),(xaser1(3),anpuls),
     &            (xaser1(4),plenth),(xaser1(5),toff),(xaser1(6),pmax),
     &            (fl(1),flshort),(fl(2),fllong)
      CHARACTER*80 line
      NAMELIST /newrun/ ak0,ak1,ak2,ak3,ak4,ak5,ak6,ak7,ak8,ak9,icxrl,
     & xz2,xz1,xmass2,xmass1,fne2,bneum,deltat,deuter,dtemax,dtimax,
     & dumax,dtmax,istage,gammae,gammai,heliu3,heliu4,hydrog,xamda1,
     & idflg,xetral,xtrlms,rhoini,rini,scp,scr,itbflg,drmul,scrho,scte,
     & scti,sctime,teini,tiini,istflg,idrflg,tinucl,tritiu,tstop,xmass,
     & xtra,xz,ipuflg,mesh,ncase,ndump,nfilm,ngeom,nhdcpy,nlp,llp,
     & nitmax,np1,np2,np3,nprnt,nproc,nrep,nup,lup
     & ,nlabs,nlbrms,nlburn,nlcri1,
     & nldepo,nldump,nlmax,nlecon,nlemp,nlfilm,nlfuse,nlhcpy,nlicon,
     & ropmul,bbtrap,nlite,nlmove,nlpfe,nlpfi,nlprnt,nlx,nlhf,nltnl,
     & tpon,fbtrap,xaser1,piq,emid,ngroup,ebdy,tpoff,rmpd,nlres,nresm,
     & nonlin,nout,nprint,nin,npunch,dlamda,nrun,mxdump,nadump,npdump,
     & nvdump,flshort,fllong,nlched,nlhead,nlomt1,nlomt2,nlomt3,nlrept,
     & ifrsta,flimit,drglas,roglas,zglas,fne,fhot,fl,ilosta,fthot,anabs,
     & saha,state,pondf,ton,fwide,efbmul,ihista,pmult,anpuls,plenth,
     & toff,pmax,rhogas,igstat,gauss,rhot,drplas,roplas,zplas,rf1,rf2,
     & rf3,nst,nfl,nt1,nt2
      READ(nread,newrun,end=300)
c 100  continue
c      REWIND(unit=nread)
c      WRITE(nprint,*)'  III =============================== III'
c      WRITE(nprint,*)'  II  Input data file was as follows:  II'
c      WRITE(nprint,*)'  III =============================== III'
c 200  READ(nread,10100,end=300)line
c      WRITE(nprint,*)line
c      GOTO 200
 300  IF(mesh.ge.kk-1)THEN
       WRITE(nprint,*)' warning:- mesh is too large'
       WRITE(nprint,*)'           should be less than ',kk-1
       WRITE(nprint,*)'           increase kk in parameter statmnt'
       WRITE(nprint,*)'           ================================'
      ENDIF
      close(unit=nread)
      RETURN
10100 FORMAT(a)
10200 FORMAT(1x,a)
      END


      SUBROUTINE mprint(k)
      PARAMETER(kk=1001,underf=1.E-30)
      LOGICAL nlend,nlres
      REAL altime,cptime,stime
      COMMON /combas/ altime,cptime,stime,ndiary,nin,nledge,nonlin,nout,
     &                nprint,npunch,nread,nrec,nresm,nrun,nstep,nlend,
     &                nlres
c end of basic system parameters
      REAL brems1(kk),brems3(kk),ddz1(kk),ddz3(kk),degen(kk),degmax,
     &     degmin,effz(kk),eixch2(kk),fidash(kk),fz1(kk),fz3(kk),
     &     fzsq1(kk),fzsq3(kk),xlc(kk),xmieff(kk),xne(kk),xni(kk),
     &     omega1(kk),pmass
      COMMON /comie / brems1,brems3,ddz1,ddz3,degen,degmax,degmin,effz,
     &                eixch2,fidash,fz1,fz3,fzsq1,fzsq3,xlc,xmieff,xne,
     &                xni,omega1,pmass
c end of ions and electrons
      REAL dm(kk),p3(kk),pini,r1(kk),r3(kk),r5(kk),rho1(kk),rho3(kk),
     &     rho5(kk),rhoini,rhor,rini,time,u2(kk),u4(kk),uedge,v1(kk),
     &     v3(kk),v5(kk)
      COMMON /comhyd/ dm,p3,pini,r1,r3,r5,rho1,rho3,rho5,rhoini,rhor,
     &                rini,time,u2,u4,uedge,v1,v3,v5
c end of hydrodynamics
      REAL ddroe1(kk),ddroe3(kk),ddroi1(kk),ddroi3(kk),ddte1(kk),
     &     ddte3(kk),ddti1(kk),ddti3(kk),gammae,gammai,xappae(kk),
     &     xappai(kk),pe1(kk),pe3(kk),pi1(kk),pi3(kk),te1(kk),te3(kk),
     &     teini,ti1(kk),ti3(kk),tiini
      COMMON /comth / ddroe1,ddroe3,ddroi1,ddroi3,ddte1,ddte3,ddti1,
     &                ddti3,gammae,gammai,xappae,xappai,pe1,pe3,pi1,pi3,
     &                te1,te3,teini,ti1,ti3,tiini
c end of thermodynamics
      INTEGER nabs1
      REAL alpha1(kk),elas1,xamda1,xaser1(kk),xecri1,plas1,rabs1,rocri1,
     &     xl1(kk),xl3(kk)
      COMMON /comlas/ alpha1,elas1,xamda1,xaser1,xecri1,plas1,rabs1,
     &                rocri1,xl1,xl3,nabs1
c end of laser variables
      REAL deuter,f1d(kk),f1h(kk),f1he3(kk),f1he4(kk),f1neu(kk),
     &     f1ntrl(kk),f1t(kk),f1x(kk),f1x1(kk),f1x2(kk),f1x3(kk),
     &     f1x4(kk),f3d(kk),f3h(kk),f3he3(kk),f3he4(kk),f3neu(kk),
     &     f3ntrl(kk),f3t(kk),f3x(kk),f3x1(kk),f3x2(kk),f3x3(kk),
     &     f3x4(kk),heliu3,heliu4,hydrog,xetral,xtrlms,pneut1,pneut3,
     &     r1dd(kk),r1dhe3(kk),r1dt(kk),r3dd(kk),r3dhe3(kk),r3dt(kk),
     &     rneut1,rneut3,tinucl,totneu,tritiu,xmass,xmass1,xmass2,
     &     xmass3,xmass4,xtra,xtra1,xtra2,xtra3,xtra4,xz,xz1,xz2,xz3,
     &     xz4,ye1(kk),ye3(kk),yi1(kk),yi3(kk),yield
      COMMON /comfus/ deuter,f1d,f1h,f1he3,f1he4,f1neu,f1ntrl,f1t,f1x,
     &                f1x1,f1x2,f1x3,f1x4,f3d,f3h,f3he3,f3he4,f3neu,
     &                f3ntrl,f3t,f3x,f3x1,f3x2,f3x3,f3x4,heliu3,heliu4,
     &                hydrog,xetral,xtrlms,pneut1,pneut3,r1dd,r1dhe3,
     &                r1dt,r3dd,r3dhe3,r3dt,rneut1,rneut3,tinucl,totneu,
     &                tritiu,xmass,xmass1,xmass2,xmass3,xmass4,xtra,
     &                xtra1,xtra2,xtra3,xtra4,xz,xz1,xz2,xz3,xz4,ye1,
     &                ye3,yi1,yi3,yield
c end of thermonuclear reactions
c c2.6 - energies
      COMMON /comen / eerror,en,einput,eloss,eefuse,eifuse,eindt1,
     &                eindt3,eneutr,pv,usqm,xionen,en0,eions,ecbbs,
     &                ecbfs,edcvs,erbfs,eatis,ezsts,erbbs,erfbs,edfbs
c c2.6 - end energies
      INTEGER mstep,ncase,ngeom
      LOGICAL mlbrms,mlecon,mlfuse,mlicon,mlx,nlabs,nlbrms,nlburn,
     &        nlcri1,nldepo,nlecon,nlfuse,nlicon,nlmove,nlpfe,nlpfi,nlx,
     &        nlhf,nltnl
      REAL pmin,rhomin,temin,timin,umin
      COMMON /comcon/ mstep,ncase,ngeom,mlbrms,mlecon,mlfuse,mlicon,mlx,
     &                nlabs,nlbrms,nlburn,nlcri1,nldepo,nlecon,nlfuse,
     &                nlicon,nlmove,nlpfe,nlpfi,nlx,nlhf,nltnl,pmin,
     &                rhomin,temin,timin,umin
c end of physics control
      INTEGER nceldt,ncondt,nit,nitmax
      LOGICAL break,nlgoon,nlite
      REAL ak0,ak1,ak2,ak3,ak4,ak5,ak6,ak7,ak8,ak9,bneum,deltat,dt2,dt3,
     &     dt4,dtemax,dtimax,dtmax,dtmin,dumax,rdt2,rdt3,rdt4
      COMMON /comnc / nceldt,ncondt,nit,nitmax,break,nlgoon,nlite,ak0,
     &                ak1,ak2,ak3,ak4,ak5,ak6,ak7,ak8,ak9,bneum,deltat,
     &                dt2,dt3,dt4,dtemax,dtimax,dtmax,dtmin,dumax,rdt2,
     &                rdt3,rdt4
c end of numerical control parameters
      INTEGER mesh,nj,njm1,njp1,nl,nlm1,nlp1
      REAL ae(kk),ai(kk),be(kk),bi(kk),ce(kk),ci(kk),de(kk),di(kk),
     &     e(kk),f(kk),ge(kk),gi(kk),q2(kk),q4(kk),teite(kk),tiite(kk),
     &     uite(kk)
      COMMON /comnum/ mesh,nj,njm1,njp1,nl,nlm1,nlp1,ae,ai,be,bi,ce,ci,
     &                de,di,e,f,ge,gi,q2,q4,teite,tiite,uite
c end of mesh and numerical methods
      INTEGER maxdim,maxrun,ndump,nrep
      LOGICAL nldump,nlemp
      REAL piq(302),tstop
      COMMON /comadm/ maxdim,maxrun,ndump,nrep,nldump,nlemp,piq,tstop
c c4.1 - end administrative variables
      INTEGER nfilm,nhdcpy,np1,np2,np3
      LOGICAL nlfilm,nlhcpy,nlprnt
      REAL buf1(kk),buf2(kk),buf3(kk),buf4(kk),buf5(kk),buf6(kk),scp,
     &     scr,scrho,scte,scti,sctime
      REAL nprnt,nproc
      COMMON /comout/ nfilm,nhdcpy,np1,np2,np3,nlfilm,nlhcpy,nlprnt,
     &        ztnext,nprnt,nproc,buf1,buf2,buf3,buf4,buf5,buf6,scp,
     &                scr,scrho,scte,scti,sctime,nt1,nt2
c end of input-output control variables
c end of input-output control variables
c x-ray laser common blocks
      DIMENSION fl(2)
      COMMON /xrl   / ropmul,tpon,tpoff,rmpd,ztrtop,dlamda,fl,fwide,
     &                efbmul,drmul,bbtrap,fbtrap,icxrl,istage,idflg,
     &                itbflg,istflg,idrflg,ipuflg,nlp,nup,llp,lup,nst,
     &                nfl,nlmax,ifrsta,ilosta,ihista,igstat
c instantaneous rates
      COMMON/wrats/weefus,weifus,wrad,wionz,wcbbz,wcbfz,wrbfz,
     & watiz,wrfbz,wdfbz,wrbbz
      COMMON /atom  / pu(10,kk),p(10,kk),pz(10,kk),zs(10,kk),zst(kk),
     &                zstz(kk),z(kk),a(kk),fn(10),deg(10),qs(10,kk),
     &                zstmi(kk),zstu(kk),engu(10,kk),eng(10,kk),
     &                rdrhe(kk),zstin(kk),nmax
c
      DIMENSION buff1(kk),buff2(kk),buff3(kk),buff4(kk),buff5(kk),
     &          buff6(kk),buff22(kk)
      DATA iburn,ifirst/0,0/
c
      avog=6.0232E23
c 0.5 print summaryrun sheet.
      IF(nstep.eq.0)CALL summaryrun
c 1. print headings, times etc.
c avoid duplication of final printing
c     iavoid = (nstep/nprnt)*nprnt - nstep
c     if ( (k.ne.3) .or. (iavoid.ne.0) ) then
      IF(.not.(k.eq.3.and.break))THEN
       CALL page
       zdt2=dt2*sctime
       ztime=time*sctime
       WRITE(nout,10400)nstep,ztime,ztrtop*sctime,zdt2
       WRITE(80,10401)nstep,ztime*1e12
       WRITE(nout,10500)
       IF(ncondt.gt.0)WRITE(nout,11700)ncondt,nceldt
       CALL blines(1)
c      1.1 boundary values
       zredge=r1(njp1)*scr
       zuedge=uedge*scr/sctime
c**    change from m to cm
       ab1=zredge*1.0E+02
       ab2=zuedge*1.0E+02
c      notice that at this point p3 really refers to level 1
       zpedge=p3(nlp1)*scp
       zsound=sqrt(abs(p3(nl))*v1(nl)*gammai)*scr/sctime
       zteb=te1(nlp1)*scte
       ztib=ti1(nlp1)*scti
       WRITE(nout,10600)ab1,ab2
       CALL blines(1)
c 0.   change units :
c**    for table elements
c**    m to cm (*1.0e+02)
c**    kg to g (*1.0e+03) so kg/m**3 to g/cm**3 (*1.0e-03)
c**    j/m**3 to megabar (*1.0e-11)
c**    degree k to ev (*8.6173e-05)
c**    no longer take logs(base 10) of buf3,buf4,buf5,buf6
       IF(nstep.eq.0)THEN
        CALL header(13)
        WRITE(13,*)mesh
       ENDIF
       WRITE(13,*)time*sctime,ztrtop*sctime
       ldum=np1
       IF(abs(piq(31)).lt.underf)THEN
        buf1(np2+np3)=ab1*1.0E-2
        buf2(np2+np3)=ab2*1.0E-2
        DO 20 j=np1,np2,np3
         buff1(j)=1.0E+02*buf1(j)
         buff22(j)=(buf1(j+np3)+buf1(j))*0.5*1.0E2
         buff2(j)=1.0E+02*buf2(j)
         buff3(j)=1.0E-03*(10**buf3(j))
         buff4(j)=1.0E-11*buf4(j)
c       buff4(j)=1.0e-11*(10**buf4(j))
         buff5(j)=8.6173E-05*(10**buf5(j))
         buff6(j)=8.6173E-05*(10**buf6(j))
         IF(abs(xmieff(j)).gt.underf)THEN
          dene=buff3(j)*fz3(j)*avog/xmieff(j)
         ELSE
          dene=0.0
         ENDIF
         WRITE(13,10300)buff1(j),buff2(j),buff22(j),buff3(j),buff4(j),
     &                  buff6(j),buff5(j),fz3(j),dene
 20     CONTINUE
       ELSEIF(piq(31).gt.0.0)THEN
        buf1(mesh+1)=ab1*1.0E-2
        buf2(mesh+1)=ab2*1.0E-2
        DO 40 j=1,mesh
         buff1(j)=1.0E+02*buf1(j)
         buff22(j)=(buf1(j)+buf1(j+1))*0.5*1.0E2
         buff2(j)=1.0E+02*buf2(j)
         buff3(j)=1.0E-03*(10**buf3(j))
c         buff4(j)=1.0e-11*(10**buf4(j))
         buff4(j)=1.0E-11*buf4(j)
         buff5(j)=8.6173E-05*(10**buf5(j))
         buff6(j)=8.6173E-05*(10**buf6(j))
         IF(abs(xmieff(j)).gt.underf)THEN
          dene=buff3(j)*fz3(j)*avog/xmieff(j)
         ELSE
          dene=0.0
         ENDIF
         WRITE(13,10300)buff1(j),buff2(j),buff22(j),buff3(j),buff4(j),
     &                  buff6(j),buff5(j),fz3(ldum),dene
         ldum=ldum+np3
 40     CONTINUE
       ENDIF
       WRITE(13,10200)ab1,ab2
c**   input value of pi for conversion purposes
       pi=3.1415927E00
c**   conversion factors (depend on value of ngeom) for output
       zng1=1.0E-04
       zng2=2.0*pi*1.0E-02
       zng3=4.0*pi
c 1.2 energies
       zef=eefuse+eifuse+eneutr
c energy scaling factor
       zsce=scrho*scr**5*sctime**2
       zpv=pv*zsce
       zusqm=usqm*zsce
       zef=(eefuse+eifuse+eneutr)*zsce
       zerror=eerror*zsce
       zrhor=rhor*scrho*scr
c**   change units :
c**   from j/m**2 to j/cm**2
       IF(ngeom.eq.1)THEN
        cd0=zsce*zng1
        cd1=zpv*zng1
        cd2=zusqm*zng1
        cd3=zef*zng1
        cd4=zerror*zng1
        WRITE(nout,12700)cd1,cd2,cd3,cd4
c**   from j/m/rad to j/cm
       ELSEIF(ngeom.eq.2)THEN
        cd0=zsce*zng2
        cd1=zpv*zng2
        cd2=zusqm*zng2
        cd3=zef*zng2
        cd4=zerror*zng2
        WRITE(nout,12800)cd1,cd2,cd3,cd4
c**   from j/st to j
       ELSEIF(ngeom.eq.3)THEN
        cd0=zsce*zng3
        cd1=zpv*zng3
        cd2=zusqm*zng3
        cd3=zef*zng3
        cd4=zerror*zng3
        WRITE(nout,12900)cd1,cd2,cd3,cd4
       ENDIF
       CALL blines(1)
c 1.3 neutrons
c have thermonuclear reactions started?
       IF(yield.gt.0.0)THEN
c have they finished?
        IF(.not.((.not.mlfuse).and.nlfuse))THEN
         zeneut=eneutr*zsce
         WRITE(nout,11600)yield,totneu,rneut1,zeneut
         CALL blines(1)
        ENDIF
       ENDIF
c 1.4 the laser
       IF(ncase.eq.1)THEN
c is the laser switched on?
        IF(abs(xaser1(9)).gt.underf)THEN
c critical density at this distance
         zrcrit=(r1(nabs1)+r1(nabs1+1))*0.5*scr
c**   change from m to cm
         ef3=zrcrit*1.0E+02
c notice that laser power and energy are not scaled
c**   change units :
c**   from /m**2 to /cm**2
         IF(ngeom.eq.1)THEN
          ef1=plas1*zng1*zsce
          ef2=elas1*zng1*zsce
          WRITE(nout,13000)ef1,ef2,(einput-en0)*zng1*zsce
c**   from /m/rad to /cm
         ELSEIF(ngeom.eq.2)THEN
          ef1=plas1*zng2*zsce
          ef2=elas1*zng2*zsce
          WRITE(nout,13100)ef1,ef2,(einput-en0)*zng2*zsce
c**   from /st to /1
         ELSEIF(ngeom.eq.3)THEN
          ef1=plas1*zng3*zsce
          ef2=elas1*zng3*zsce
          WRITE(nout,13200)ef1,ef2,(einput-en0)*zng3*zsce
         ENDIF
         WRITE(nout,13300)ef3,nabs1
         CALL blines(1)
c 1.5 record details of switch off of laser
c has the laser been switched off since last output?
        ELSEIF(abs(xaser1(10)).gt.underf)THEN
         il1=int(xaser1(11))
         zl2=xaser1(12)*sctime
         zl3=xaser1(13)
         zl4=xaser1(14)
         zl5=xaser1(15)*scr
         WRITE(nout,11800)il1,zl2,zl3,zl4,zl5
c reset signal after message has been recorded
         xaser1(10)=0.0
        ENDIF
       ENDIF
c 1.6 radiation (bremsstrahlung)
c in the current version we only print the radiation loss
c     if (.not. nlbrms) goto 200
c     if (.not. mlbrms) goto 200
       zeloss=eloss*cd0
       zeions=eions*cd0
       zecbbs=ecbbs*cd0
       zecbfs=ecbfs*cd0
       zedcvs=edcvs*cd0
       zerbfs=erbfs*cd0
       zeatis=eatis*cd0
       zerbbs=erbbs*cd0
       zerfbs=erfbs*cd0
       zedfbs=edfbs*cd0
       IF(.not.nlpfe.and.piq(58).lt.1.5)THEN
        WRITE(nout,11900)zeloss
       ELSE
        WRITE(nout,12100)zeloss,zerbbs,zerfbs,zeions
        WRITE(nout,12200)zecbbs,zecbfs,zerbfs,zeatis
        wrtot=wrad+wrfbz+wrbbz
        WRITE(22,3333)time,-wrad,-wrfbz,-wrad-wrfbz,-wrbbz,-wrtot
c write input files for fly
        if(ifirst.eq.0) then
         ifirst=1
         do 62 j=1,mesh
          nfly=30+j
c         write(nfly,*)'time  te  ne  size'
62       continue
        endif
        do 63 j=1,mesh
          nfly=30+j
          teevs=te3(j)/11604.5
          dnecm3=zst(j)*avog/xmieff(j)*rho3(j)*1.e-3
          size=(r3(j+1)-r3(j))*100.
c         write(nfly,3333)time,teevs,dnecm3,size
63      continue
       ENDIF
c hot electron parameters
       CALL blines(1)
       hetev=8.6173E-05*xaser1(kk)
       IF(abs(piq(21)).gt.underf)WRITE(nout,13600)hetev,
     &                                 (xaser1(j),j=52,61)
c 2. unload output buffers
c the reciprocal scaling factors are being printed
       zscp=1.0/scp
       zscr=1.0/scr
       zscrho=1.0/scrho
       zscte=1.0/scte
       zscti=1.0/scti
       zscu=sctime/scr
c explanation of what the columns are
c**   table consists of 8 columns
c**   describe what each variable is and give units used
c**   this is printed only at timestep 0
       CALL blines(1)
       IF(nstep.eq.0)THEN
        WRITE(nout,10700)
        IF(piq(31).ge.0.0)WRITE(nout,10800)
        IF(piq(32).ge.0.0)WRITE(nout,10900)
        IF(piq(31).ge.0.0)WRITE(nout,11000)
        IF(piq(33).ge.0.0)WRITE(nout,11100)
        IF(piq(34).ge.0.0)WRITE(nout,11200)
        IF(piq(36).ge.0.0)WRITE(nout,11300)
        IF(piq(35).ge.0.0)WRITE(nout,11400)
        WRITE(nout,11500)
        CALL blines(2)
       ENDIF
c**   the following assumes that the values of piq(x) , x from 31 to 36,
c**   are the same at any given instant in time
c print table headings
       IF(piq(31).ge.0.0)THEN
        WRITE(nout,12500)
        WRITE(nout,12600)
        CALL blines(1)
       ENDIF
c allow for full-partial-no printing of output buffers
c**   remember that buffer fz3 (i.e.average z) is printed as follows
c**   fz3(l),l=np1,np2,np3
       l=np1
c**   the following assumes that the values of piq(x) , x from 31 to 36,
c**   are the same at any given instant in time
c**   this partially prints the results
       IF(abs(piq(31)).lt.underf)THEN
        buff1(np2+np3)=ab1
        buff2(np2+np3)=ab2
        DO 60 i=np1,np2,np3
         buff22(i)=(buff1(i)+buff1(i+np3))*0.5
         WRITE(nout,13400)i,buff1(i),buff2(i),buff22(i),buff3(i),
     &                    buff4(i),buff6(i),buff5(i),fz3(i),effz(i)
 60     CONTINUE
        WRITE(nout,13500)np2+np3,buff1(np2+np3),buff2(np2+np3)
c**   this fully prints the results
       ELSEIF(piq(31).gt.0.0)THEN
        buff1(mesh+1)=ab1
        buff2(mesh+1)=ab2
        DO 80 i=1,mesh
         buff22(i)=(buff1(i+1)+buff1(i))*0.5
         WRITE(nout,13400)i,buff1(i),buff2(i),buff22(i),buff3(i),
     &                    buff4(i),buff6(i),buff5(i),fz3(l),effz(i)
c**   remember that buffer fz3 is not printed from 1 to mesh
         l=l+np3
 80     CONTINUE
        WRITE(nout,13500)mesh+1,buff1(mesh+1),buff2(mesh+1)
       ENDIF
c hot electron printout
       CALL user(k)
       CALL hdump(1)
       IF(nlemp)CALL hdump(2)
       IF(.not.((k.eq.3).and.nlburn))THEN
        IF(.not.((nlfuse.and.(.not.mlfuse)).and.((yield.gt.0.0).and.
     &     nlburn)))RETURN
       ENDIF
      ENDIF
c     endif
c 3. special print-out of burn-up of material
      IF((.not.nlburn))RETURN
      IF(yield.le.0.0)RETURN
c has the burn-up been printed previously?
      IF(iburn.ne.0)RETURN
      iburn=1
      j1=1
      j2=40
c     call page
 100  WRITE(nout,12300)
      CALL blines(1)
      DO 200 l=j1,j2
       z1=f1d(l)
       z2=f1h(l)
       z3=f1he3(l)
       z4=f1he4(l)
       z5=f1neu(l)
       z6=f1ntrl(l)
       z7=f1t(l)
       z8=f1x(l)
       WRITE(nout,12400)z1,z2,z3,z4,z5,z6,z7,z8
 200  CONTINUE
      IF(j2.ge.mesh)RETURN
      j1=j1+40
      j2=j2+40
      IF(j2.gt.maxdim)j2=maxdim-1
      GOTO 100
10100 FORMAT(i5)
10200 FORMAT(1p,2E11.4)
 3333 format(1x,1pe14.7,7(1x,1pe11.4))
10300 FORMAT(1p,3E12.4,2E11.3,3E10.2,3e12.4)
c 99. format statements
c**   9905-9912 are for explanation of columns in table
c**   9920,9921 are for table headings
c**   9930-9932 are for energies output (depends on value of ngeom)
c**   9940-9943 are for laser power output (depends on value of ngeom)
c**   9950 is for printing results in table
10400 FORMAT(15x,'timestep number ',i6,5x,'time = ',1pe14.7,6x,
     &       'time r-to-p = ',1pe14.7,6x,'delta t = ',1pe14.7)
10401 FORMAT(15x,'timestep number ',i6,5x,' = ',f14.2,' ps')
10500 FORMAT(15x,21('-'),6x,21('-'),6x,28('-'),6x,24('-'))
10600 FORMAT(10x,'boundary :  r (cm) = ',1pe10.4,'  u (cm/sec) = ',
     &       1pe11.4)
10700 FORMAT(10x,'column 1 : cell number')
10800 FORMAT(10x,'column 2 : cell edge (cm)')
10900 FORMAT(10x,'column 3 : hydrodynamic velocities (cm/sec)')
11000 FORMAT(10x,'column 4 : cell centre(cm)')
11100 FORMAT(10x,'column 5 : density (g/cm**3)')
11200 FORMAT(10x,'column 6 : hydrodynamic pressure (Mb)')
11300 FORMAT(10x,'column 7 : electron temperature (ev)')
11400 FORMAT(10x,'column 8 : ion temperature (ev)')
11500 FORMAT(10x,'column 9 : average z')
11600 FORMAT(10x,'fusion   :  yield ',1pe12.5,'   neutrons ',1pe12.5,
     &       '   rate ',1pe12.5,'   energy ',1pe12.5)
11700 FORMAT(41x,'delta t determined by condition ',i3,
     &       '  at meshpoint ',i3)
11800 FORMAT(3x,'***   laser switched off at step ',i5,'   time ',
     &       1pe14.7,'   power ',1pe12.5,'   energy ',1pe12.5,
     &       '   r(abs) ',1pe12.5)
11900 FORMAT(10x,'radiation loss      ff:',9x,1pe12.5,9x)
12000 FORMAT(10x,'radiation loss      ff:',9x,1pe12.5,9x,'bb:',1pe12.5,
     &       9x,'fb:',1pe12.5)
12100 FORMAT(10x,'radiation loss      ff:',9x,1pe12.5,9x,'bb:',1pe12.5,
     &       9x,'fb:',1pe12.5,/10x,'Internal energy of ionisation : ',
     &       1pe12.5)
12200 FORMAT(10x,'coll deexcit-excit  :',3x,1pe12.5,3x,'coll rec-ion :',
     &       1pe12.5,3x,'Photoabs  :',1pe12.5,' eatis : ',1pe12.5)
12300 FORMAT(8x,'deuterium',6x,'hydrogen',7x,'helium3',8x,'helium4',8x,
     &       'neutrons',7x,'neutrals',7x,'tritium',8x,'extra')
12400 FORMAT(7x,8(e12.5,3x))
c9920 format(10x,'cell',6x,'coords',9x,'velocity',9x,'density',8x,
c    +'pressure',8x,'elec temp',7x,'ion temp',10x,'aver z')
12500 FORMAT(10x,'cell',3x,'cell edge',4x,'velocity',5x,'cell centre',
     &       2x,'density',6x,'pressure',5x,'elec temp',4x,'ion temp',5x,
     &       '   z*         z')
12600 FORMAT(10x,4('='),3x,9('='),4x,8('='),5x,11('='),2x,7('='),6x,
     &       8('='),5x,9('='),4x,8('='),5x,9('='),5x,5('='))
12700 FORMAT(10x,'energies (j/cm**2) :  thermal ',1pe12.5,'   kinetic ',
     &       1pe12.5,'   nuclear ',1pe12.5,'   error ',1pe12.5)
12800 FORMAT(10x,'energies (j/cm) :  thermal ',1pe12.5,'   kinetic ',
     &       1pe12.5,'   nuclear ',1pe12.5,'   error ',1pe12.5)
12900 FORMAT(10x,'energies (j) :  thermal ',1pe12.5,'   kinetic ',
     &       1pe12.5,'   nuclear ',1pe12.5,'   error ',1pe12.5)
13000 FORMAT(10x,'laser power (w/cm**2) ',1pe12.5,
     &       '   total energy input from laser (j/cm**2) ',1pe12.5,1x,
     &       ' absorbed (j/cm**2) ',1pe12.5)
13100 FORMAT(10x,'laser power (w/cm) ',1pe12.5,
     &       '   total energy input from laser (j/cm) ',1pe12.5,1x,
     &       ' absorbed (j/cm) ',1pe12.5)
13200 FORMAT(10x,'laser power (w) ',1pe12.5,
     &       '   total energy input from laser (j) ',1pe12.5,1x,
     &       ' absorbed (j) ',1pe12.5)
13300 FORMAT(/10x,'absorption at r (cm) = ',1pe11.4,5x,'at meshpoint ',
     &       i3)
c9950 format(10x,i3,4x,7(1pe13.6,3x))
13400 FORMAT(10x,i3,3x,8(1pe11.4,2x),1pg11.4)
13500 FORMAT(10x,i3,3x,2(1pe11.4,2x))
13600 FORMAT(10x,'hot electron temperature',1pe12.5,' ev',/,10x,
     &       'fractional energy loss in each group',/,
     &       (10x,5(1pe10.3,1x)))
      END


      SUBROUTINE header(istream)
      PARAMETER(kk=1001)
c 1.4 define data specific to run
c c1.1 basic system parameters
      CHARACTER label1*80,label2*80,label3*80,label4*80,label5*80,
     &          label6*80,label7*80,label8*80
      COMMON /comlab/ label1,label2,label3,label4,label5,label6,label7,
     &                label8
c end of basic system parameters
c c1.9 development and diagnostic parameters
c end of development and diagnostic parameters
c c2.1 hydrodynamics
      DIMENSION dm(kk),p3(kk),r1(kk),r3(kk),r5(kk),rho1(kk),rho3(kk),
     &          rho5(kk),u2(kk),u4(kk),v1(kk),v3(kk),v5(kk),lstat(kk)
      COMMON /comhyd/ dm,p3,pini,r1,r3,r5,rho1,rho3,rho5,rhoini,rhor,
     &                rini,time,u2,u4,uedge,v1,v3,v5
c end of hydrodynamics
c c2.2 thermodynamics
c end of thermodynamics
c c2.4 laser variables
      DIMENSION alpha1(kk),xaser1(kk),xl1(kk),xl3(kk)
      COMMON /comlas/ alpha1,elas1,xamda1,xaser1,xecri1,plas1,rabs1,
     &                rocri1,xl1,xl3,nabs1
c end of laser variables
c c2.5 thermonuclear reactions
      DIMENSION f1d(kk),f1h(kk),f1he3(kk),f1he4(kk),f1neu(kk),f1ntrl(kk)
     &          ,f1t(kk),f1x(kk),f1x1(kk),f1x2(kk),f1x3(kk),f1x4(kk),
     &          f3d(kk),f3h(kk),f3he3(kk),f3he4(kk),f3neu(kk),f3ntrl(kk)
     &          ,f3t(kk),f3x(kk),f3x1(kk),f3x2(kk),f3x3(kk),f3x4(kk),
     &          r1dd(kk),r1dhe3(kk),r1dt(kk),r3dd(kk),r3dhe3(kk),
     &          r3dt(kk),ye1(kk),ye3(kk),yi1(kk),yi3(kk)
      COMMON /comfus/ deuter,f1d,f1h,f1he3,f1he4,f1neu,f1ntrl,f1t,f1x,
     &                f1x1,f1x2,f1x3,f1x4,f3d,f3h,f3he3,f3he4,f3neu,
     &                f3ntrl,f3t,f3x,f3x1,f3x2,f3x3,f3x4,heliu3,heliu4,
     &                hydrog,xetral,xtrlms,pneut1,pneut3,r1dd,r1dhe3,
     &                r1dt,r3dd,r3dhe3,r3dt,rneut1,rneut3,tinucl,totneu,
     &                tritiu,xmass,xmass1,xmass2,xmass3,xmass4,xtra,
     &                xtra1,xtra2,xtra3,xtra4,xz,xz1,xz2,xz3,xz4,ye1,
     &                ye3,yi1,yi3,yield
c end of thermonuclear reactions
c c2.7 physics control
      LOGICAL mlbrms,mlecon,mlfuse,mlicon,mlx,nlabs,nlbrms,nlburn,
     &        nlcri1,nldepo,nlecon,nlfuse,nlicon,nlmove,nlpfe,nlpfi,nlx,
     &        nlhf,nltnl
      COMMON /comcon/ mstep,ncase,ngeom,mlbrms,mlecon,mlfuse,mlicon,mlx,
     &                nlabs,nlbrms,nlburn,nlcri1,nldepo,nlecon,nlfuse,
     &                nlicon,nlmove,nlpfe,nlpfi,nlx,nlhf,nltnl,pmin,
     &                rhomin,temin,timin,umin
c end of physics control
c c3.1 numerical control parameters
c end of numerical control parameters
c c3.2 mesh and numerical methods
      DIMENSION ae(kk),ai(kk),be(kk),bi(kk),ce(kk),ci(kk),de(kk),di(kk),
     &          e(kk),f(kk),ge(kk),gi(kk),q2(kk),q4(kk),teite(kk),
     &          tiite(kk),uite(kk)
      COMMON /comnum/ mesh,nj,njm1,njp1,nl,nlm1,nlp1,ae,ai,be,bi,ce,ci,
     &                de,di,e,f,ge,gi,q2,q4,teite,tiite,uite
c end of mesh and numerical methods
c c4.1 - administrative variables
      LOGICAL nldump,nlemp
      DIMENSION piq(302)
      COMMON /comadm/ maxdim,maxrun,ndump,nrep,nldump,nlemp,piq,tstop
c c4.1 - end administrative variables
c c5.1 input-output control variables
c x-ray laser common blocks
      DIMENSION fl(2)
      COMMON /xrl   / ropmul,tpon,tpoff,rmpd,ztrtop,dlamda,fl,fwide,
     &                efbmul,drmul,bbtrap,fbtrap,icxrl,istage,idflg,
     &                itbflg,istflg,idrflg,ipuflg,nlp,nup,llp,lup,nst,
     &                nfl,nlmax,ifrsta,ilosta,ihista,igstat
c end of input-output control variables
c   new subroutine source control factors and layer materials
      COMMON /sornew/ rf1,rf2,rf3,fne2
      EQUIVALENCE(piq(63),zplas),(piq(62),roplas),(piq(61),drplas),
     &            (piq(59),state),(piq(58),saha),(piq(56),pondf),
     &            (piq(57),rhot),(piq(37),gauss),(piq(52),anabs),
     &            (piq(22),fthot),(piq(21),fhot),(piq(14),fne),
     &            (piq(13),zglas),(piq(12),roglas),(piq(11),drglas),
     &            (piq(10),flimit),(rhoini,rhogas),(xaser1(1),ton),
     &            (xaser1(2),pmult),(xaser1(3),anpuls),
     &            (xaser1(4),plenth),(xaser1(5),toff),(xaser1(6),pmax),
     &            (fl(1),flshort),(fl(2),fllong)
      CHARACTER*80 buff
      CHARACTER c4*4,c10*10,c5*5
      buff(1:48)=label1(1:48)
      buff(49:54)='ngeom:'
      CALL itoc(ngeom,c4)
      buff(55:64)=c4
      buff(65:70)=' pmax:'
      CALL rtoc(pmax,c10)
      buff(71:80)=c10
      WRITE(istream,10100)buff
      buff(1:7)='xamda1:'
      CALL rtoc(xamda1,c10)
      buff(8:17)=c10
      buff(18:25)=' plenth:'
      CALL rtoc(plenth,c10)
      buff(26:35)=c10
      buff(36:42)=' pmult:'
      CALL rtoc(pmult,c10)
      buff(43:52)=c10
      buff(53:58)=' fhot:'
      CALL rtoc(fhot,c10)
      buff(59:68)=c10
      buff(69:72)=' xz:'
      CALL rtoc5(xz,c5)
      buff(73:77)=c5
      buff(78:80)='   '
      WRITE(istream,10100)buff
      buff(1:6)='xmass:'
      CALL rtoc5(xmass,c5)
      buff(7:11)=c5
      buff(12:17)=' rini:'
      CALL rtoc(rini,c10)
      buff(18:27)=c10
      buff(28:32)=' xz2:'
      CALL rtoc5(xz2,c5)
      buff(33:37)=c5
      buff(38:45)=' xmass2:'
      CALL rtoc5(xmass2,c5)
      buff(46:50)=c5
      buff(51:56)=' fne2:'
      CALL rtoc5(fne2,c5)
      buff(57:61)=c5
      buff(62:69)=' drglas:'
      CALL rtoc(piq(11),c10)
      buff(70:79)=c10
      buff(80:80)=' '
      WRITE(istream,10100)buff

      buff(1:4)='xz1:'
      CALL rtoc5(xz1,c5)
      buff(5:9)=c5
      buff(10:17)=' xmass1:'
      CALL rtoc5(xmass1,c5)
      buff(18:22)=c5
      buff(23:30)=' drplas:'
      CALL rtoc(piq(61),c10)
      buff(31:40)=c10
      buff(41:46)=' mesh:'
      CALL itoc(mesh,c4)
      buff(47:50)=c4
      buff(51:57)=' zglas:'
      CALL rtoc5(piq(13),c5)
      buff(58:62)=c5
      buff(63:69)=' zplas:'
      CALL rtoc5(piq(63),c5)
      buff(70:74)=c5
      buff(75:80)='      '
      WRITE(istream,10100)buff
      buff(1:4)='rf1:'
      CALL rtoc(rf1,c10)
      buff(5:14)=c10
      buff(15:19)=' rf2:'
      CALL rtoc(rf2,c10)
      buff(20:29)=c10
      buff(30:34)=' rf3:'
      CALL rtoc(rf3,c10)
      buff(35:44)=c10
      IF(icxrl.eq.0)buff(45:52)=' xrl off'
      IF(icxrl.eq.1)buff(45:52)=' xrl on '
      IF(istage.le.3)THEN
       buff(53:60)=' dlamda:'
       CALL rtoc(1.0E8*dlamda,c10)
       buff(61:70)=c10
       buff(71:80)='Angstroms '
      ENDIF
      WRITE(istream,10100)buff
      rtot=rini+piq(11)+piq(61)
      IF(ngeom.eq.2)THEN
       pequiv=pmax/rtot
      ELSEIF(ngeom.eq.3)THEN
       pequiv=pmax/(rtot*rtot)
      ELSEIF(ngeom.eq.1)THEN
       pequiv=pmax
      ENDIF
      buff(1:7)='pequiv:'
      CALL rtoc(pequiv,c10)
      buff(8:17)=c10
      buff(18:24)=' w/m**2'
      IF(idflg.eq.0)THEN
       buff(25:48)=' opacity corrections off'
      ELSE
       buff(25:39)='        ropmul:'
       CALL rtoc(ropmul,c10)
       buff(39:48)=c10
      ENDIF
      buff(49:65)='                 '
      buff(66:80)='               '
      WRITE(istream,10100)buff
      RETURN
10100 FORMAT(a80)
      END

c*******************************************************************
      SUBROUTINE itoc(i,c)
      CHARACTER*(*) c
      INTEGER i
      WRITE(c,10100)i
      RETURN
10100 FORMAT(i4)
      END


      SUBROUTINE rtoc(r,c)
      CHARACTER*(*) c
      REAL r
      WRITE(c,10100)r
      RETURN
10100 FORMAT(1pe10.3)
      END


      SUBROUTINE rtoc5(r,c)
      CHARACTER*(*) c
      REAL r
      WRITE(c,10100)r
      RETURN
10100 FORMAT(f5.0)
      END


      SUBROUTINE summaryrun
      PARAMETER(kk=1001,underf=1.E-30)
c 1.4 define data specific to run
c c1.1 basic system parameters
      LOGICAL nlend,nlres
      COMMON /combas/ altime,cptime,stime,ndiary,nin,nledge,nonlin,nout,
     &                nprint,npunch,nread,nrec,nresm,nrun,nstep,nlend,
     &                nlres
c end of basic system parameters
c c1.9 development and diagnostic parameters
c end of development and diagnostic parameters
c c2.1 hydrodynamics
      DIMENSION dm(kk),p3(kk),r1(kk),r3(kk),r5(kk),rho1(kk),rho3(kk),
     &          rho5(kk),u2(kk),u4(kk),v1(kk),v3(kk),v5(kk),lstat(kk)
      COMMON /comhyd/ dm,p3,pini,r1,r3,r5,rho1,rho3,rho5,rhoini,rhor,
     &                rini,time,u2,u4,uedge,v1,v3,v5
c end of hydrodynamics
c c2.2 thermodynamics
      DIMENSION ddroe1(kk),ddroe3(kk),ddroi1(kk),ddroi3(kk),ddte1(kk),
     &          ddte3(kk),ddti1(kk),ddti3(kk),xappae(kk),xappai(kk),
     &          pe1(kk),pe3(kk),pi1(kk),pi3(kk),te1(kk),te3(kk),ti1(kk),
     &          ti3(kk)
      COMMON /comth / ddroe1,ddroe3,ddroi1,ddroi3,ddte1,ddte3,ddti1,
     &                ddti3,gammae,gammai,xappae,xappai,pe1,pe3,pi1,pi3,
     &                te1,te3,teini,ti1,ti3,tiini
c end of thermodynamics
c c2.4 laser variables
      DIMENSION alpha1(kk),xaser1(kk),xl1(kk),xl3(kk)
      COMMON /comlas/ alpha1,elas1,xamda1,xaser1,xecri1,plas1,rabs1,
     &                rocri1,xl1,xl3,nabs1
c end of laser variables
c c2.5 thermonuclear reactions
      DIMENSION f1d(kk),f1h(kk),f1he3(kk),f1he4(kk),f1neu(kk),f1ntrl(kk)
     &          ,f1t(kk),f1x(kk),f1x1(kk),f1x2(kk),f1x3(kk),f1x4(kk),
     &          f3d(kk),f3h(kk),f3he3(kk),f3he4(kk),f3neu(kk),f3ntrl(kk)
     &          ,f3t(kk),f3x(kk),f3x1(kk),f3x2(kk),f3x3(kk),f3x4(kk),
     &          r1dd(kk),r1dhe3(kk),r1dt(kk),r3dd(kk),r3dhe3(kk),
     &          r3dt(kk),ye1(kk),ye3(kk),yi1(kk),yi3(kk)
      COMMON /comfus/ deuter,f1d,f1h,f1he3,f1he4,f1neu,f1ntrl,f1t,f1x,
     &                f1x1,f1x2,f1x3,f1x4,f3d,f3h,f3he3,f3he4,f3neu,
     &                f3ntrl,f3t,f3x,f3x1,f3x2,f3x3,f3x4,heliu3,heliu4,
     &                hydrog,xetral,xtrlms,pneut1,pneut3,r1dd,r1dhe3,
     &                r1dt,r3dd,r3dhe3,r3dt,rneut1,rneut3,tinucl,totneu,
     &                tritiu,xmass,xmass1,xmass2,xmass3,xmass4,xtra,
     &                xtra1,xtra2,xtra3,xtra4,xz,xz1,xz2,xz3,xz4,ye1,
     &                ye3,yi1,yi3,yield
c end of thermonuclear reactions
c c2.7 physics control
      LOGICAL mlbrms,mlecon,mlfuse,mlicon,mlx,nlabs,nlbrms,nlburn,
     &        nlcri1,nldepo,nlecon,nlfuse,nlicon,nlmove,nlpfe,nlpfi,nlx,
     &        nlhf,nltnl
      COMMON /comcon/ mstep,ncase,ngeom,mlbrms,mlecon,mlfuse,mlicon,mlx,
     &                nlabs,nlbrms,nlburn,nlcri1,nldepo,nlecon,nlfuse,
     &                nlicon,nlmove,nlpfe,nlpfi,nlx,nlhf,nltnl,pmin,
     &                rhomin,temin,timin,umin
c end of physics control
c c3.1 numerical control parameters
c end of numerical control parameters
c c3.2 mesh and numerical methods
      DIMENSION ae(kk),ai(kk),be(kk),bi(kk),ce(kk),ci(kk),de(kk),di(kk),
     &          e(kk),f(kk),ge(kk),gi(kk),q2(kk),q4(kk),teite(kk),
     &          tiite(kk),uite(kk)
      COMMON /comnum/ mesh,nj,njm1,njp1,nl,nlm1,nlp1,ae,ai,be,bi,ce,ci,
     &                de,di,e,f,ge,gi,q2,q4,teite,tiite,uite
c end of mesh and numerical methods
c c4.1 - administrative variables
      LOGICAL nldump,nlemp
      DIMENSION piq(302)
      COMMON /comadm/ maxdim,maxrun,ndump,nrep,nldump,nlemp,piq,tstop
c c4.1 - end administrative variables
c c5.1 input-output control variables
c x-ray laser common blocks
      DIMENSION fl(2)
      COMMON /xrl   / ropmul,tpon,tpoff,rmpd,ztrtop,dlamda,fl,fwide,
     &                efbmul,drmul,bbtrap,fbtrap,icxrl,istage,idflg,
     &                itbflg,istflg,idrflg,ipuflg,nlp,nup,llp,lup,nst,
     &                nfl,nlmax,ifrsta,ilosta,ihista,igstat
      COMMON /atom  / pu(10,kk),p(10,kk),pz(10,kk),zs(10,kk),zst(kk),
     &                zstz(kk),z(kk),a(kk),fn(10),deg(10),qs(10,kk),
     &                zstmi(kk),zstu(kk),engu(10,kk),eng(10,kk),
     &                rdrhe(kk),zstin(kk),nmax
      CHARACTER*10 shname(12)
      COMMON /sname / shname
c end of input-output control variables
c   new subroutine source control factors and layer materials
      COMMON /sornew/ rf1,rf2,rf3,fne2
      EQUIVALENCE(piq(63),zplas),(piq(62),roplas),(piq(61),drplas),
     &            (piq(59),state),(piq(58),saha),(piq(56),pondf),
     &            (piq(57),rhot),(piq(37),gauss),(piq(52),anabs),
     &            (piq(22),fthot),(piq(21),fhot),(piq(14),fne),
     &            (piq(13),zglas),(piq(12),roglas),(piq(11),drglas),
     &            (piq(10),flimit),(rhoini,rhogas),(xaser1(1),ton),
     &            (xaser1(2),pmult),(xaser1(3),anpuls),
     &            (xaser1(4),plenth),(xaser1(5),toff),(xaser1(6),pmax),
     &            (fl(1),flshort),(fl(2),fllong)
c        call page
      WRITE(nprint,10500)
      CALL blines(1)
c        call daytim
      CALL blines(1)
      WRITE(nprint,10600)
      WRITE(nprint,10700)
      CALL blines(1)
      WRITE(nprint,10800)
      IF(ngeom.eq.1)WRITE(nprint,10900)
      IF(ngeom.eq.2)WRITE(nprint,11000)
      IF(ngeom.eq.3)WRITE(nprint,11100)
      CALL blines(1)
      WRITE(nprint,11300)
      CALL blines(1)
      IF(abs(gauss+1.0).lt.underf)WRITE(nprint,11500)
      IF(abs(gauss).lt.underf)WRITE(nprint,11600)
      IF(abs(gauss-1.0).lt.underf)WRITE(nprint,11700)
      WRITE(nprint,11400)xamda1
      WRITE(nprint,11800)toff
      IF(abs(gauss-1.0).lt.underf)WRITE(nprint,15400)anpuls
      IF(ngeom.eq.1)WRITE(nprint,11900)pmax
      IF(ngeom.eq.2)WRITE(nprint,15600)pmax
      IF(ngeom.eq.3)WRITE(nprint,15700)pmax
      rtot=rini+piq(11)+piq(61)
      IF(ngeom.eq.2)THEN
       pequiv=pmax/rtot
       WRITE(nprint,18200)pequiv
      ELSEIF(ngeom.eq.3)THEN
       pequiv=pmax/(rtot*rtot)
       WRITE(nprint,18200)pequiv
      ENDIF
      WRITE(nprint,12000)plenth
      WRITE(nprint,15300)pmult
      CALL blines(1)
      WRITE(nprint,12100)
      WRITE(nprint,12200)teini
      WRITE(nprint,12300)tiini
      WRITE(nprint,12400)fhot
      IF(nlpfe)THEN
       WRITE(nprint,12500)
      ELSE
       IF(abs(state).lt.underf)WRITE(nprint,12600)
       IF(abs(state-1.0).lt.underf)WRITE(nprint,12700)
       IF(abs(state-2.0).lt.underf)WRITE(nprint,12800)
       IF(abs(state-3.0).lt.underf)WRITE(nprint,12900)
      ENDIF
      CALL blines(1)
      WRITE(nprint,13000)
      WRITE(nprint,13100)xz
      WRITE(nprint,13200)xmass
      WRITE(nprint,11200)rini
      IF(abs(fne-1.0).lt.underf)THEN
       WRITE(nprint,13900)
      ELSEIF(abs(fne).lt.underf)THEN
       WRITE(nprint,14000)
      ELSE
       WRITE(nprint,14100)(fne*100.0)
       WRITE(nprint,14200)((1.0-fne)/2.0)*100.0,((1.0-fne)/2.0)*100.0
      ENDIF
      WRITE(nprint,13800)rhogas
      WRITE(nprint,18300)rf1
      CALL blines(1)
      WRITE(nprint,18800)xz2
      WRITE(nprint,18900)xmass2
      WRITE(nprint,13500)drglas
      WRITE(nprint,13600)roglas
      WRITE(nprint,18400)rf2
      WRITE(nprint,18500)fne2*100.0
      IF(abs(fne2-1.0).gt.underf)WRITE(nprint,18600)(1.0-fne2)*100.0
      CALL blines(1)
      WRITE(nprint,19000)xz1
      WRITE(nprint,19100)xmass1
      WRITE(nprint,13300)drplas
      WRITE(nprint,13400)roplas
      WRITE(nprint,18700)rf3
      WRITE(nprint,13700)hydrog
      CALL blines(1)
      WRITE(nprint,14300)
      IF(nlbrms)WRITE(nprint,14400)
      WRITE(nprint,14500)anabs
      IF(abs(saha-1.0).lt.underf)WRITE(nprint,14600)
      IF(abs(saha-2.0).lt.underf)WRITE(nprint,14700)
      CALL blines(1)
      WRITE(nprint,14800)
      WRITE(nprint,14900)mesh
      WRITE(nprint,15000)zplas
      WRITE(nprint,15100)zglas
      WRITE(nprint,15200)tstop
      WRITE(nprint,18100)
      IF(icxrl.eq.1)THEN
       WRITE(nprint,15800)
       IF(istage.eq.4)WRITE(nprint,16000)
       WRITE(nprint,18100)
       IF(istage.eq.1)WRITE(nprint,16100)
       IF(istage.eq.2)WRITE(nprint,16200)
       IF(istage.eq.3)WRITE(nprint,16300)
       IF(idflg.eq.0)WRITE(nprint,16400)
       IF(idflg.eq.1)WRITE(nprint,16500)ropmul
       IF(itbflg.eq.0)WRITE(nprint,16600)
       IF(itbflg.eq.1)WRITE(nprint,16700)
       IF(istflg.eq.0)WRITE(nprint,16800)
       IF(istflg.eq.1)WRITE(nprint,16900)
       IF(ipuflg.eq.0)WRITE(nprint,17000)
       IF(ipuflg.eq.1)THEN
        WRITE(nprint,17100)
        WRITE(nprint,17200)tpon
        WRITE(nprint,17300)tpoff
        WRITE(nprint,17400)nlp
        WRITE(nprint,17500)nup
        WRITE(nprint,17600)rmpd
       ENDIF
       WRITE(nprint,17700)dlamda*1.0E8
       WRITE(nprint,17800)nlmax
       DO 50 i=1,nlmax
        WRITE(nprint,17900)i,fl(i)/100.0
 50    CONTINUE
       WRITE(6,18000)fwide
      ELSE
       WRITE(nprint,15900)
      ENDIF
      WRITE(nprint,18100)
c ad1
      IF(abs(piq(58)-2.0).lt.underf)THEN
c write out initial nlte populations
       CALL page
       WRITE(6,10300)
       WRITE(6,10100)(shname(i),i=1,nmax),shname(12)
       WRITE(6,10200)(shname(11),i=1,nmax+1)
       DO 100 l=nst,nfl
        WRITE(6,10400)l,(p(i,l),i=1,nmax),zst(l)
 100   CONTINUE
      ENDIF
      RETURN
10100 FORMAT(6x,' cn',11(a10:))
10200 FORMAT(6x,'===',11(a10:))
10300 FORMAT(//7x,'initial nlte populations'/)
10400 FORMAT(6x,i3,11F10.5)
c---------------------------------------------------------------------
cl              3.         format statements
10500 FORMAT(48x,'summary of the input data'/48x,
     &       '-------------------------')
10600 FORMAT(31x,'version 103 supporting xrl H-, Li-, Na-like ',
     &       'recombination',//39x,
     &       'and Ne-like collisional Si(14), Ar(18), ',
     &       'Ti(22), Fe(26), Ge(32), Kr(36)')
10700 FORMAT(/38x,' :: Djaoui-Rose implementation, Aug 1991::'/38x,
     &       ' -----------------------------------------')
10800 FORMAT(55x,'geometry')
10900 FORMAT(10x,'ngeom  :',2x,'plane geometry problem')
11000 FORMAT(10x,'ngeom  :',2x,'cylindrical geometry problem')
11100 FORMAT(10x,'ngeom  :',2x,'spherical geometry problem')
11200 FORMAT(10x,'rini   :',2x,'radius of the shell filled with',
     &       ' the gas',17x,'=',e11.5,10x,'meters')
11300 FORMAT(55x,'laser')
11400 FORMAT(10x,'xamda1 :',2x,'the laser wavelength',36x,'=',e11.5,10x,
     &       'meters')
11500 FORMAT(10x,'gauss  :',2x,'the laser pulse is triangular')
11600 FORMAT(10x,'gauss  :',2x,'the laser pulse is isentropic')
11700 FORMAT(10x,'gauss  :',2x,'the laser pulse is gaussian')
11800 FORMAT(10x,'toff   :',2x,'time to switch off laser',32x,'=',e11.5,
     &       10x,'seconds')
11900 FORMAT(10x,'pmax   :',2x,'the maximum laser power',33x,'=',e11.5,
     &       10x,'watts/square meters')
12000 FORMAT(10x,'plenth :',2x,'the half width of the pulse in time',
     &       21x,'=',e11.5,10x,'seconds')
12100 FORMAT(55x,'ions and electrons')
12200 FORMAT(10x,'teini  :',2x,'initial electrons temperature',27x,'=',
     &       e11.5,10x,'degree kelvin')
12300 FORMAT(10x,'tiini  :',2x,'initial ions temperature',32x,'=',e11.5,
     &       10x,'degree kelvin')
12400 FORMAT(10x,'fhot   :',2x,'fraction of anomalous absorption',
     &       ' into hot electrons',5x,'=',f7.2)
12500 FORMAT(10x,'nlpfe  :',2x,'electron equation of state is',
     &       ' perfect gas ')
12600 FORMAT(10x,'state  :',2x,'electron equation of state is fermi',
     &       ' dirac')
12700 FORMAT(10x,'state  :',2x,'electron equation of state is',
     &       ' thomas fermi')
12800 FORMAT(10x,'state  :',2x,'electron equation of state is',
     &       ' thomas fermi with quantum corrections')
12900 FORMAT(10x,'state  :',2x,'electron equation of state is',
     &       ' thomas fermi with modified corrections')
13000 FORMAT(55x,'target')
13100 FORMAT(10x,'xz     :',2x,'charge number of extra element',26x,'=',
     &       f7.2)
13200 FORMAT(10x,'xmass  :',2x,'mass number of extra element',28x,'=',
     &       f7.2)
13300 FORMAT(10x,'drplas :',2x,'thickness of plastic coating',28x,'=',
     &       e11.5,10x,'meters')
13400 FORMAT(10x,'roplas :',2x,'density of plastic',38x,'=',e11.5,10x,
     &       'kg/cubic meters')
13500 FORMAT(10x,'drglas :',2x,'thickness of glass',38x,'=',e11.5,10x,
     &       'meters')
13600 FORMAT(10x,'roglas :',2x,'density of glass',40x,'=',e11.5,10x,
     &       'kg/cubic meters')
13700 FORMAT(10x,'hydrog :',2x,'fractional number density of ',
     &       'hydrogen in the target',5x,'=',f7.2)
13800 FORMAT(10x,'rhogas :',2x,'the gas density',41x,'=',e11.5,10x,
     &       'kg/cubic meters')
13900 FORMAT(10x,'fne    :',2x,'the composition of gas is pure xz  ')
14000 FORMAT(10x,'fne    :',2x,'the composition of gas is 50%',
     &       ' deuterium and 50% tritium')
14100 FORMAT(10x,'fne    :',2x,'the composition of gas is ',f10.5,
     &       '% xz ')
14200 FORMAT(10x,'        ',2x,f10.5,'% deuterium and ',f10.5,
     &       '% tritium')
14300 FORMAT(55x,'radiation')
14400 FORMAT(10x,'nlbrms :',2x,'bremsstrahlung radiation is',
     &       ' included')
14500 FORMAT(10x,'anabs  :',2x,'fractional anomalous absorption at',
     &       ' critical',13x,'=',f7.2)
14600 FORMAT(10x,'saha   :',2x,'the ionisation equilibrium is',
     &       ' calculated with TF approximation')
14700 FORMAT(10x,'saha   :',2x,'the ionisation equilibrium is',
     &       ' calculated with nlte package')
14800 FORMAT(55x,'numerical')
14900 FORMAT(10x,'mesh   :',2x,'the total number of mesh points',25x,
     &       '=',i6)
15000 FORMAT(10x,'zplas  :',2x,'number of mesh points in plastic',24x,
     &       '=',f7.0)
15100 FORMAT(10x,'zglas  :',2x,'number of mesh points in glass',26x,'=',
     &       f7.0)
15200 FORMAT(10x,'tstop  :',2x,'maximum value of time permited',26x,'=',
     &       e11.5,10x,'seconds')
15300 FORMAT(10x,'pmult  :',2x,'initial power level',37x,'=',e11.5,10x,
     &       'watts')
15400 FORMAT(10x,'anpuls :',2x,'number of gaussian pulses',31x,'=',
     &       e11.5)
15500 FORMAT(10x,'nozon  :',2x,'total number of zones',35x,'=',i6)
15600 FORMAT(10x,'pmax   :',2x,'the maximum laser power',33x,'=',e11.5,
     &       10x,'watts/meters/rad')
15700 FORMAT(10x,'pmax   :',2x,'the maximum laser power',33x,'=',e11.5,
     &       10x,'watts/sterad')
15800 FORMAT(10x,'x-ray laser calculations are on')
15900 FORMAT(10x,'x-ray laser calculations are off')
16000 FORMAT(10x,'* Ne-like scheme option *  ')
16100 FORMAT(10x,'istage :',2x,'AA model ionisation stage : H-like')
16200 FORMAT(10x,'istage :',2x,'AA model ionisation stage : Li-like')
16300 FORMAT(10x,'istage :',2x,'AA model ionisation stage : Na-like')
16400 FORMAT(10x,'idflg  :',2x,'no optical depth correction      ')
16500 FORMAT(10x,'idflg  :',2x,'optical depth correction  of : ',
     &       1pe15.5)
16600 FORMAT(10x,'itbflg :',2x,'no thermal band in AA')
16700 FORMAT(10x,'itbflg :',2x,
     &       'population of highest level forced into lte in AA')
16800 FORMAT(10x,'istflg :',2x,
     &       'motional stark broadening not included in AA')
16900 FORMAT(10x,'istflg :',2x,
     &       'motional stark broadening on balmer alpha in AA')
17000 FORMAT(10x,'ipuflg :',2x,'photopumping off in AA')
17100 FORMAT(10x,'ipuflg :',2x,'photopumping on in AA')
17200 FORMAT(10x,'tpon   :',2x,'time pump on   ',41x,'=',1pe15.5,6x,
     &       'seconds')
17300 FORMAT(10x,'tpoff  :',2x,'time pump off  ',41x,'=',1pe15.5,6x,
     &       'seconds')
17400 FORMAT(10x,'nlp    :',2x,'lower pump level : ',37x,'=',i2)
17500 FORMAT(10x,'nup    :',2x,'upper pump level : ',37x,'=',i2)
17600 FORMAT(10x,'rmpd   :',2x,'modal photon density : ',33x,'=',
     &       1pe15.5)
17700 FORMAT(10x,'dlamda :',2x,'x-ray laser wavelength   ',31x,'=',
     &       1pe15.5,6x,'Angstroms')
17800 FORMAT(10x,'nlmax  :',2x,
     &       'no of lengths used in intensity calculation   ',10x,'=',
     &       i2)
17900 FORMAT(10x,'fl(',i1,')  :',2x,'length',50x,'=',1pe15.5,6x,
     &       'meters')
18000 FORMAT(10x,'width  :',2x,'width of line focus  (cm)   ',28x,'=',
     &       1pe15.5)
18100 FORMAT(10x,'   ')
18200 FORMAT(10x,'pequiv :',2x,'equivalent surface irradiance',27x,'=',
     &       e11.5,10x,'watts/square meter')
18300 FORMAT(10x,'rf1    :',2x,'gas fill grid factor         ',27x,'=',
     &       f10.8)
18400 FORMAT(10x,'rf2    :',2x,'glass shell grid factor      ',27x,'=',
     &       f10.8)
18500 FORMAT(10x,'fne2   :',2x,'fraction of xz2              ',27x,'=',
     &       f10.6,11x,'%')
18600 FORMAT(10x,'        ',2x,'fraction of silicon          ',27x,'=',
     &       f10.6,11x,'%')
18700 FORMAT(10x,'rf3    :',2x,'plastic shell grid factor    ',27x,'=',
     &       f10.8)
18800 FORMAT(10x,'xz2    :',2x,'charge number of glass shell element',
     &       20x,'=',f7.2)
18900 FORMAT(10x,'xmass2 :',2x,'mass number of glass shell element',22x,
     &       '=',f7.2)
19000 FORMAT(10x,'xz1    :',2x,'charge number of plastic shell element',
     &       18x,'=',f7.2)
19100 FORMAT(10x,'xmass1 :',2x,'mass number of plastic shell element',
     &       20x,'=',f7.2)
      END


      SUBROUTINE radran
      PARAMETER(kk=1001,underf=1.E-30)
c this subroutine is an interface between medusa
c and the non-lte routines
c medusa common blocks
      LOGICAL nldump,nlemp
      DIMENSION piq(302)
      COMMON /comadm/ maxdim,maxrun,ndump,nrep,nldump,nlemp,piq,tstop
c c4.1 - end administrative variables
      INTEGER nabs1
      REAL alpha1(kk),elas1,xamda1,xaser1(kk),xecri1,plas1,rabs1,rocri1,
     &     xl1(kk),xl3(kk)
      COMMON /comlas/ alpha1,elas1,xamda1,xaser1,xecri1,plas1,rabs1,
     &                rocri1,xl1,xl3,nabs1
c end of laser variables
      DIMENSION ae(kk),ai(kk),be(kk),bi(kk),ce(kk),ci(kk),de(kk),di(kk),
     &          e(kk),f(kk),ge(kk),gi(kk),q2(kk),q4(kk),teite(kk),
     &          tiite(kk),uite(kk)
      COMMON /comnum/ mesh,nj,njm1,njp1,nl,nlm1,nlp1,ae,ai,be,bi,ce,ci,
     &                de,di,e,f,ge,gi,q2,q4,teite,tiite,uite
c end of mesh and numerical methods
      LOGICAL nlfilm,nlhcpy,nlprnt
      DIMENSION buf1(kk),buf2(kk),buf3(kk),buf4(kk),buf5(kk),buf6(kk)
      REAL nprnt,nproc,nprnti
      COMMON /comout/ nfilm,nhdcpy,np1,np2,np3,nlfilm,nlhcpy,nlprnt,
     &        ztnext,nprnt,nproc,buf1,buf2,buf3,buf4,buf5,buf6,scp,
     &                scr,scrho,scte,scti,sctime,nt1,nt2
c end of input-output control variables
      DIMENSION dm(kk),p3(kk),r1(kk),r3(kk),r5(kk),rho1(kk),rho3(kk),
     &          rho5(kk),u2(kk),u4(kk),v1(kk),v3(kk),v5(kk),lstat(kk)
      COMMON /comhyd/ dm,p3,pini,r1,r3,r5,rho1,rho3,rho5,rhoini,rhor,
     &                rini,time,u2,u4,uedge,v1,v3,v5
c end of hydrodynamics
      LOGICAL nlend,nlres
      COMMON /combas/ altime,cptime,stime,ndiary,nin,nledge,nonlin,nout,
     &                nprint,npunch,nread,nrec,nresm,nrun,nstep,nlend,
     &                nlres
c end of basic system parameters
      LOGICAL break,nlgoon,nlite
      COMMON /comnc / nceldt,ncondt,nit,nitmax,break,nlgoon,nlite,ak0,
     &                ak1,ak2,ak3,ak4,ak5,ak6,ak7,ak8,ak9,bneum,deltat,
     &                dt2,dt3,dt4,dtemax,dtimax,dtmax,dtmin,dumax,rdt2,
     &                rdt3,rdt4
c end of numerical control parameters
      DIMENSION ddroe1(kk),ddroe3(kk),ddroi1(kk),ddroi3(kk),ddte1(kk),
     &          ddte3(kk),ddti1(kk),ddti3(kk),xappae(kk),xappai(kk),
     &          pe1(kk),pe3(kk),pi1(kk),pi3(kk),te1(kk),te3(kk),ti1(kk),
     &          ti3(kk)
      COMMON /comth / ddroe1,ddroe3,ddroi1,ddroi3,ddte1,ddte3,ddti1,
     &                ddti3,gammae,gammai,xappae,xappai,pe1,pe3,pi1,pi3,
     &                te1,te3,teini,ti1,ti3,tiini
c end of thermodynamics
      DIMENSION brems1(kk),brems3(kk),ddz1(kk),ddz3(kk),degen(kk),
     &          effz(kk),eixch2(kk),fidash(kk),fz1(kk),fz3(kk),fzsq1(kk)
     &          ,fzsq3(kk),xlc(kk),xmieff(kk),xne(kk),xni(kk),omega1(kk)
      COMMON /comie / brems1,brems3,ddz1,ddz3,degen,degmax,degmin,effz,
     &                eixch2,fidash,fz1,fz3,fzsq1,fzsq3,xlc,xmieff,xne,
     &                xni,omega1,pmass
c end of ions and electrons
c end of thermonuclear reactions
c c2.6 - energies
c c2.6 - end energies
c c2.7 physics control
      LOGICAL mlbrms,mlecon,mlfuse,mlicon,mlx,nlabs,nlbrms,nlburn,
     &        nlcri1,nldepo,nlecon,nlfuse,nlicon,nlmove,nlpfe,nlpfi,nlx,
     &        nlhf,nltnl
      COMMON /comcon/ mstep,ncase,ngeom,mlbrms,mlecon,mlfuse,mlicon,mlx,
     &                nlabs,nlbrms,nlburn,nlcri1,nldepo,nlecon,nlfuse,
     &                nlicon,nlmove,nlpfe,nlpfi,nlx,nlhf,nltnl,pmin,
     &                rhomin,temin,timin,umin
c end of physics control
c non-lte dimension statements
      DIMENSION zm(kk),am(kk),rh(kk),te(kk),ti(kk),rb(kk),uu(kk)
c x-ray laser common blocks
      DIMENSION fl(2)
      COMMON /xrl   / ropmul,tpon,tpoff,rmpd,ztrtop,dlamda,fl,fwide,
     &                efbmul,drmul,bbtrap,fbtrap,icxrl,istage,idflg,
     &                itbflg,istflg,idrflg,ipuflg,nlp,nup,llp,lup,nst,
     &                nfl,nlmax,ifrsta,ilosta,ihista,igstat
      COMMON /dpi   / ecbb(kk),ecbf(kk),eraddr(kk),erbf(kk),eati(kk),
     &                dcvz(kk),efbx2(kk),efrx2(kk),erbb(kk),erfb(kk),
     &                acfc(10,kk),acfr(10,kk),bcfc(10,kk),bcfr(10,kk),
     &                edfb(kk),ter(kk),dqdtim(kk),nmaxr
      COMMON /atom  / pu(10,kk),p(10,kk),pz(10,kk),zs(10,kk),zst(kk),
     &                zstz(kk),z(kk),a(kk),fn(10),deg(10),qs(10,kk),
     &                zstmi(kk),zstu(kk),engu(10,kk),eng(10,kk),
     &                rdrhe(kk),zstin(kk),nmax


      COMMON /radflg/ yamda,xhnu,xihnu,dihnu,xplas1(kk),xrbf(kk),
     &                atrbf(kk),corhf1(kk),corhf(kk),iradf1,iradf2
      DIMENSION dfac(kk)
      SAVE ist,itbfli
      DATA ist/0/

      IF(abs(piq(58)-2.).gt.underf)RETURN
      an=6.0232E23
      zk=1.38E-23/(gammae-1.)
c setup values required at the start of the calculation
      IF(ist.eq.0)THEN
c this section is only gone through once at beginning of calculation
c setup numbers of cells to be included in the calculation
c ad1  if not specified in input use 1 to mesh
       IF((nfl-nst).le.0.or.nst.lt.1.or.nfl.gt.mesh)THEN
        nst=1
        nfl=mesh
       ENDIF
c copy initial z, a and density and temperature
       DO 50 l=1,mesh
        zm(l)=effz(l)
        am(l)=xmieff(l)
        rh(l)=0.001*rho3(l)
        ter(l)=te3(l)
        te(l)=ter(l)/11605000.0
        rb(l)=100.0*r3(l)
        uu(l)=100.0*u2(l)
 50    CONTINUE
c       setup initial conditions

       iradf1=0
       iradf2=0
       yamda=xamda1
c       xhnu is xray photon energy in kev
       xhnu=1.24E-7/(xamda1*100.)

       CALL setup(zm,am,rh,te,dt2,mesh,icxrl,istage,rb,uu,nlres)
       nprnti=nprnt
       itbfli=itbflg
       ist=ist+1
       RETURN
      ENDIF
c=======================================================================
c this section is gone through every timestep
c=======================================================================
c     advance the populations by one time-step
      IF(nstep.le.1)RETURN
c select max timestep and frequency of printout for 0 fs pulse
c     if(time .lt. 1.e-14) then
c        nprnt=-1.e-3
c        dtmax=1.e-15/2.
c     elseif(time.lt.1.e-11) then
c        nprnt=-2.e-2
c        dtmax=2.e-14/2.
c     elseif(time.lt.1.e-10) then
c        nprnt=-5.e-1
c        dtmax=5.e-13/2.
c     elseif(time.lt.1.e-09) then
c        nprnt=-5.
c        dtmax=5.e-12/2.
c     else
c         nprnt=nprnti
c     endif
c select max timestep and frequency of printout for 100 fs pulse
c     if(time .lt. 0.7e-12) then
c        nprnt=-1.e-1
c        dtmax=1.e-13/2.
c     elseif(time.lt.1.3e-12) then
c        nprnt=-1.e-2
c        dtmax=1.e-14/2.
c     elseif(time.lt.1.e-11) then
c        nprnt=-1.e-1
c        dtmax=1.e-13/2.
c     elseif(time.lt.1.e-09) then
c        nprnt=-5.
c        dtmax=5.e-12/2.
c     else
c         nprnt=nprnti
c     endif
c select max timestep and frequency of printout for 1 ps   pulse
c     if(time.lt.1.e-10) then
c        nprnt=-1.e00
c        dtmax=5.e-13/2.
c     elseif(time.lt.1.e-09) then
c        nprnt=-5.
c        dtmax=5.e-12/2.
c     else
c         nprnt=nprnti
c     endif
c select max timestep and frequency of printout for 10 ps   pulse
c     if(time .lt. 7e-12) then
c        nprnt=-1.
c        dtmax=1.e-12/2.
c     elseif(time.lt.1.5e-12) then
c        nprnt=-5.e-2
c        dtmax=5.e-14/2.
c     elseif(time.lt.1.e-10) then
c        nprnt=-5.e-1
c        dtmax=5.e-13/2.
c     elseif(time.lt.1.e-09) then
c        nprnt=-5.
c        dtmax=5.e-12/2.
c     else
c         nprnt=nprnti
c     endif
c
      tstep=dt2
      IF(piq(95).le.0.0)piq(95)=0.150
      fcvt=piq(95)
      DO 100 l=1,mesh
       rh(l)=0.001*rho3(l)
       ti(l)=ti3(l)/11605000.0
       rb(l)=100.0*r3(l)
       uu(l)=100.0*u2(l)
       dfac(l)=1.
 100  CONTINUE
      rb(mesh+1)=100.0*r3(mesh+1)
      uu(mesh+1)=100.0*u2(mesh+1)

c11111111111111111111111111111111111111
      IF(nit.eq.1)THEN
c     for first iteration  solve for temp without ionization
       CALL abcd
       e(1)=ce(1)/be(1)
       f(1)=(de(1)+ge(1))/be(1)
       DO 150 l=2,nl
        znume=1.0/(be(l)-ae(l)*e(l-1))
        e(l)=ce(l)*znume
        f(l)=(de(l)+ge(l)-ae(l)*f(l-1))*znume
 150   CONTINUE
c         find te starting from the boundary
       te3(nl)=f(nl)
       DO 200 ldash=1,nlm1
        l=nl-ldash
        te3(l)=f(l)-e(l)*te3(l+1)
 200   CONTINUE
      ENDIF
c     111111111111111111111111111111111111111

c     ionisation energy in each cell. every iteration
      DO 300 l=mesh,1,-1
c       energy input
        einc=ddte3(l)*te3(l)-ddte1(l)*te1(l)-efbx2(l)*dt2
        einc2=einc/dt2
c       thermal energy
        cvt=ddte1(l)*te1(l)
c       22222222222
c       temperature at which rate of ionisation is calculated
        te3t=te1(l)+0.5*einc/(ddte1(l)+ddte3(l))
        ter(l)=0.5*te1(l)+0.5*max(te1(l),te3t)
        te(l)=ter(l)/11605000.0
        rho=rh(l)
        theta=te(l)

        thetar=theta
c       average x-ray irradiance to kev/cm**2/sec
        xihnu=xplas1(l+1)*1.E-04/(1.60219E-16)/12.56637
c       ?number of photons with specific direction of polarization 1./2
c       xihnu=xihnu/2.
c       setup radiation field
        IF(iradf1.ne.0)CALL setr(thetar)
        IF(iradf2.ne.0)CALL setzi(thetar)
c       set itbflg = 0 if laser is on
        itbflg=itbfli

c       ad1 allan: if tunnel ionisation is included and laser is on -
c       make sure thermal band flag is zero
        IF(nltnl.and.(xplas1(l+1).gt.1.0))itbflg=0

c       calculate new populations, energy rate of relaxation towards eq
        xtest=5.e-6
        if(icxrl.eq.1.and.istage.le.3)xtest=5.E-7
        CALL rate(rho,theta,nstep,time,tstep,l,rb,uu,xtest)
c if TF do not double count ionisation energy
        if(.not.nlpfe.and.piq(59).gt.0.5) goto 300
c     if(l.eq.34)then
c     write(6,*)mstep,dqdtim(l),ecbb(l),ecbf(l)
c    &,eraddr(l),erfb(l),edfb(l),erbb(l),efrx2(l),erbf(l),
c    & zst(l)
c     endif

        efbx2(l)=dqdtim(l)

c       222222222222
c       if tf eos used do not double count ionisation energy
         xfbmax=0.0
c        if system recombines no limit on e ->dump into free elec
         IF (dqdtim(l).gt.0.) goto 300
c        system is ionising how much energy used in ionisation
         IF (abs(dqdtim(l)).gt.underf) THEN
c          ionisation limit to ;
           if(einc2.gt.underf) then
             fionz=-dqdtim(l)/(einc2-dqdtim(l))
             xfbmax=fionz*einc+fcvt*cvt
           else
            xfbmax=fcvt*cvt
           endif
c          limit energy loss through ionisation + excitation
           dfac(l)=1.0
           IF ((-dqdtim(l)*tstep).gt.xfbmax) THEN
             IF ((-dqdtim(l)*tstep).gt.underf)
     &       dfac(l)=-xfbmax /(dqdtim(l)*tstep)
           ENDIF
         ENDIF
c        end if abs(dqdtim(l))
c        xxxxxxxxxxxxxxxxxxx  rescaling col rates
         IF (dfac(l).lt.1.0.and.dfac(l).ge.0.0) THEN
           efbx2(l)=efbx2(l)*dfac(l)
           ecbb(l)=ecbb(l)*dfac(l)
           ecbf(l)=ecbf(l)*dfac(l)
           eraddr(l)=eraddr(l)*dfac(l)
           erbb(l)=erbb(l)*dfac(l)
           erfb(l)=erfb(l)*dfac(l)
           edfb(l)=edfb(l)*dfac(l)
           DO 310 i=1,nmaxr
c            calculate dp with rescaled collisional rates
             acf=acfc(i,l)*dfac(l)+acfr(i,l)
             bcf=bcfc(i,l)*dfac(l)+bcfr(i,l)
           p(i,l)=(pu(i,l)+(acf*tstep))/((((acf/deg(i))+bcf)*tstep)+1.0)
 310       CONTINUE
c          get corresponding eng and z*
           CALL pe(rho,theta,l)
         ENDIF
c        if(l.eq.34) then
c        write(6,*)nit,te1(l),te3(l),cvt,einc2,efbx2(l),
c    &    zst(l),dfac(l)
c        endif

 300  CONTINUE
      RETURN
      END


      SUBROUTINE timstp
      PARAMETER(kk=1001,underf=1.E-30)
c 2.16 most restrictive value of timestep
c c1.9 development and diagnostic parameters
      DIMENSION nadump(20),npdump(20),nvdump(20)
      LOGICAL nlched,nlhead(9),nlomt1(50),nlomt2(50),nlomt3(50),nlrept
      COMMON /comddp/ maxdum,mxdump,nadump,nclass,npdump,npoint,nsub,
     &                nvdump,nlched,nlhead,nlomt1,nlomt2,nlomt3,nlrept
c end of development and diagnostic parameters
c c2.1 hydrodynamics
      DIMENSION dm(kk),p3(kk),r1(kk),r3(kk),r5(kk),rho1(kk),rho3(kk),
     &          rho5(kk),u2(kk),u4(kk),v1(kk),v3(kk),v5(kk),lstat(kk)
      COMMON /comhyd/ dm,p3,pini,r1,r3,r5,rho1,rho3,rho5,rhoini,rhor,
     &                rini,time,u2,u4,uedge,v1,v3,v5
c end of hydrodynamics
c c2.2 thermodynamics
      DIMENSION ddroe1(kk),ddroe3(kk),ddroi1(kk),ddroi3(kk),ddte1(kk),
     &          ddte3(kk),ddti1(kk),ddti3(kk),xappae(kk),xappai(kk),
     &          pe1(kk),pe3(kk),pi1(kk),pi3(kk),te1(kk),te3(kk),ti1(kk),
     &          ti3(kk)
      COMMON /comth / ddroe1,ddroe3,ddroi1,ddroi3,ddte1,ddte3,ddti1,
     &                ddti3,gammae,gammai,xappae,xappai,pe1,pe3,pi1,pi3,
     &                te1,te3,teini,ti1,ti3,tiini
c end of thermodynamics
c c2.3 ions and electrons
      DIMENSION brems1(kk),brems3(kk),ddz1(kk),ddz3(kk),degen(kk),
     &          effz(kk),eixch2(kk),fidash(kk),fz1(kk),fz3(kk),fzsq1(kk)
     &          ,fzsq3(kk),xlc(kk),xmieff(kk),xne(kk),xni(kk),omega1(kk)
      COMMON /comie / brems1,brems3,ddz1,ddz3,degen,degmax,degmin,effz,
     &                eixch2,fidash,fz1,fz3,fzsq1,fzsq3,xlc,xmieff,xne,
     &                xni,omega1,pmass
c end of ions and electrons
c c2.4 laser variables
      DIMENSION alpha1(kk),xaser1(kk),xl1(kk),xl3(kk)
      COMMON /comlas/ alpha1,elas1,xamda1,xaser1,xecri1,plas1,rabs1,
     &                rocri1,xl1,xl3,nabs1
c end of laser variables
c c2.7 physics control
      LOGICAL mlbrms,mlecon,mlfuse,mlicon,mlx,nlabs,nlbrms,nlburn,
     &        nlcri1,nldepo,nlecon,nlfuse,nlicon,nlmove,nlpfe,nlpfi,nlx,
     &        nlhf,nltnl
      COMMON /comcon/ mstep,ncase,ngeom,mlbrms,mlecon,mlfuse,mlicon,mlx,
     &                nlabs,nlbrms,nlburn,nlcri1,nldepo,nlecon,nlfuse,
     &                nlicon,nlmove,nlpfe,nlpfi,nlx,nlhf,nltnl,pmin,
     &                rhomin,temin,timin,umin
c end of physics control
c c3.1 numerical control parameters
      LOGICAL break,nlgoon,nlite
      COMMON /comnc / nceldt,ncondt,nit,nitmax,break,nlgoon,nlite,ak0,
     &                ak1,ak2,ak3,ak4,ak5,ak6,ak7,ak8,ak9,bneum,deltat,
     &                dt2,dt3,dt4,dtemax,dtimax,dtmax,dtmin,dumax,rdt2,
     &                rdt3,rdt4
c end of numerical control parameters
c c3.2 mesh and numerical methods
      DIMENSION ae(kk),ai(kk),be(kk),bi(kk),ce(kk),ci(kk),de(kk),di(kk),
     &          e(kk),f(kk),ge(kk),gi(kk),q2(kk),q4(kk),teite(kk),
     &          tiite(kk),uite(kk)
      COMMON /comnum/ mesh,nj,njm1,njp1,nl,nlm1,nlp1,ae,ai,be,bi,ce,ci,
     &                de,di,e,f,ge,gi,q2,q4,teite,tiite,uite
c end of mesh and numerical methods
c c5.1 input-output control variables
      LOGICAL nlfilm,nlhcpy,nlprnt
      DIMENSION buf1(kk),buf2(kk),buf3(kk),buf4(kk),buf5(kk),buf6(kk)
      REAL nprnt,nproc
      COMMON /comout/ nfilm,nhdcpy,np1,np2,np3,nlfilm,nlhcpy,nlprnt,
     &        ztnext,nprnt,nproc,buf1,buf2,buf3,buf4,buf5,buf6,scp,
     &                scr,scrho,scte,scti,sctime,nt1,nt2
c end of input-output control variables
c c2.6 - energies
      COMMON /comen / eerror,en,einput,eloss,eefuse,eifuse,eindt1,
     &                eindt3,eneutr,pv,usqm,xionen,en0,eions,ecbbs,
     &                ecbfs,edcvs,erbfs,eatis,ezsts,erbbs,erfbs,edfbs
c c2.6 - end energies
c c4.1 - administrative variables
      LOGICAL nldump,nlemp
      DIMENSION piq(302)
      COMMON /comadm/ maxdim,maxrun,ndump,nrep,nldump,nlemp,piq,tstop
c c4.1 - end administrative variables

      COMMON /atom  / pu(10,kk),p(10,kk),pz(10,kk),zs(10,kk),zst(kk),
     &                zstz(kk),z(kk),a(kk),fn(10),deg(10),qs(10,kk),
     &                zstmi(kk),zstu(kk),engu(10,kk),eng(10,kk),
     &                rdrhe(kk),zstin(kk),nmax
      COMMON /dpi   / ecbb(kk),ecbf(kk),eraddr(kk),erbf(kk),eati(kk),
     &                dcvz(kk),efbx2(kk),efrx2(kk),erbb(kk),erfb(kk),
     &                acfc(10,kk),acfr(10,kk),bcfc(10,kk),bcfr(10,kk),
     &                edfb(kk),ter(kk),dqdtim(kk),nmaxr

      COMMON /radflg/ yamda,xhnu,xihnu,dihnu,xplas1(kk),xrbf(kk),
     &                atrbf(kk),corhf1(kk),corhf(kk),iradf1,iradf2
      DATA iclass,isub/2,16/
      break=.false.
      IF(nlomt2(isub))RETURN
c 1. maximum sound speed. changes in v and te, ti.
      zd1=0.0
      zd2=1.0E-7
      zd3=1.0E-7
      zd4=1.0E-7
      DO 100 l=1,nl
c        if(te3(l) .le. 0.0) goto 400
c        if(ti3(l) .le. 0.0) goto 400
       j=l
c condition 75
       zz1=sqrt(abs(p3(l))*v3(l))/(r3(j+1)-r3(j))
       IF(zz1.ge.zd1)THEN
        zd1=zz1
        il1=l
       ENDIF
       IF(l.ne.1)THEN
        zz1=0.5*sqrt(abs(p3(l-1))*v3(l-1))/(r3(j+1)-r3(j))
        IF(zz1.ge.zd1)THEN
         zd1=zz1
         il1=l
        ENDIF
       ENDIF
       IF(l.ne.nl)THEN
        zz1=0.5*sqrt(abs(p3(l+1))*v3(l+1))/(r3(j+1)-r3(j))
        IF(zz1.ge.zd1)THEN
         zd1=zz1
         il1=l
        ENDIF
       ENDIF
c condition 76
       zz2=abs(v3(l)-v1(l))/(v3(l)+v1(l))
       IF(zz2.ge.zd2)THEN
        zd2=zz2
        il2=l
       ENDIF
       zz3=abs(ti3(l)-ti1(l))/(ti3(l)+ti1(l)+timin)
       IF(zz3.ge.zd3)THEN
        zd3=zz3
        il3=l
       ENDIF
       zz4=abs(te3(l)-te1(l))/(te3(l)+te1(l)+temin)
c in the strongly degenerate case it is the change in pe which
c determines dt4
       IF((.not.nlpfe).and.(degen(l).lt.0.01))zz4=abs(pe3(l)-pe1(l))
     &    /(pe3(l)+pe1(l)+pmin)
       IF(zz4.ge.zd4)THEN
        zd4=zz4
        il4=l
       ENDIF
 100  CONTINUE
c remember the gamma value in the formulae for the soundspeed
      z1=ak1/(sqrt(gammai)*zd1)
      z2=0.5*ak2/zd2*dt2
      z3=0.5*ak3/zd3*dt2
      z4=0.5*ak4/zd4*dt2
c if no motion takes place we ignore the c-f-l condition
      IF(.not.nlmove)z1=max(2.0*z2,2.0*z3,2.0*z4)
      zmin=min(z1,z2,z3,z4)
c now choose the most restrictive value
      dt4=max(dtmin,min(dtmax,zmin))
c which cell puts the restriction on dt?
      IF(abs(dt4-z1).lt.underf)THEN
       ncondt=1
       nceldt=il1
      ELSEIF(abs(dt4-z2).lt.underf)THEN
       ncondt=2
       nceldt=il2
      ELSEIF(abs(dt4-z3).lt.underf)THEN
       ncondt=3
       nceldt=il3
      ELSEIF(abs(dt4-z4).lt.underf)THEN
       ncondt=4
       nceldt=il4
      ELSE
       ncondt=0
       nceldt=0
      ENDIF
      zdt4=dt4
c 2. take account of changes at the boundary
      IF(ncase.ne.1)THEN
       IF(ncase.eq.2)THEN
c 2.2 temperatures determined
c find the boundary temperatures one step ahead, but save old values
        zesave=te3(nlp1)
        zisave=ti3(nlp1)
        tdash=time+dt4
        CALL boundy(tdash)
        zd3=abs(ti3(nlp1)-ti1(nlp1))/(ti3(nlp1)+ti1(nlp1)+timin)
        zd4=abs(te3(nlp1)-te1(nlp1))/(te3(nlp1)+te1(nlp1)+temin)
c if there are no changes we avoid dividing by zero
        IF(abs(zd3).gt.underf.or.abs(zd4).gt.underf)THEN
         IF(abs(zd3).lt.underf)zd3=1.0E-10
         IF(abs(zd4).lt.underf)zd4=1.0E-10
         z3=ak3/zd3*dt4
         z4=ak4/zd4*dt4
         dt4=min(dt4,z3,z4)
c put back old values
         ti3(nlp1)=zisave
         te3(nlp1)=zesave
         IF(abs(dt4-zdt4).gt.underf)THEN
          ncondt=4
          IF(z4.gt.z3)ncondt=3
          nceldt=nlp1
         ENDIF
        ENDIF
       ELSEIF(ncase.eq.3)THEN
c 2.3 pressure determined
c find the boundary pressure one step ahead, but save old value
        zpsave=p3(nlp1)
        tdash=time+dt4
        CALL boundy(tdash)
        zd1=sqrt(abs(p3(nlp1))*v3(nl))/(r3(njp1)-r3(nj))
        zd1=zd1*sqrt(gammai)
        dt4=min(dt4,ak1/zd1)
c put back old value
        p3(nlp1)=zpsave
        IF(abs(dt4-zdt4).gt.underf)THEN
         nceldt=nlp1
         ncondt=5
        ENDIF
       ELSEIF(ncase.eq.4)THEN
c 2.4 velocity determined
c find the boundary velocity one step ahead, but save old value
        zusave=u4(njp1)
        tdash=time+0.5*dt4
        CALL boundy(tdash)
        zd1=abs(u4(njp1))/(r3(njp1)-r3(nj))
c avoid dividing by zero
        IF(abs(zd1).lt.underf)zd1=ak1/dt4
        dt4=min(dt4,ak1/zd1)
c put back old value
        u4(njp1)=zusave
        IF(abs(dt4-zdt4).gt.underf)THEN
         nceldt=nlp1
         ncondt=5
        ENDIF
       ELSE
        CALL gotoer(' MEDUSA: timstp: goto error, ncase is',ncase)
        GOTO 200
       ENDIF
       GOTO 400
      ENDIF
c 2.1 thermally insulated wall
c save current values of nabs1 and plas1
 200  iabs1=nabs1
      zlas1=plas1
      tdash=time+dt4
      CALL laser(tdash)

      IF(abs(plas1).gt.underf)THEN
c we save the previous values of xl3 and alpha1 as they are required
c at the next timestep
       CALL copyr(xl3,1,buf1,1,nlp1)
       CALL copyr(alpha1,1,buf2,1,nlp1)
       CALL copyr(xplas1,1,buf3,1,nlp1)
c find out where the laser light is absorbed
       CALL absorb
       IF(nabs1.gt.1)THEN
c the pressure increase caused by the absorption of laser light
c will approximately produce a displacement of the nearby point as
c (nabs is the cell where the light is absorbed)
        zdr=r3(nabs1)**(ngeom-1)
     &      *plas1*dt4*dt4*dt4/(dm(nabs1)*(r3(nabs1)-r3(nabs1-1)))
c now make sure that the boundaries of the cell nabs can not cross
c the neighbouring cell-boundaries
        zdrabs=min(r3(nabs1+1)-r3(nabs1),r3(nabs1)-r3(nabs1-1))
        IF(abs(zdr).gt.underf)dt4=min(dt4,zdrabs/zdr*dt4)
        IF(abs(dt4-zdt4).gt.underf)THEN
         nceldt=nabs1
         ncondt=6
        ENDIF
       ENDIF
c bring back the old values of xl3 and alpha1
       CALL copyr(buf1,1,xl3,1,nlp1)
       CALL copyr(buf2,1,alpha1,1,nlp1)
       CALL copyr(buf3,1,xplas1,1,nlp1)
      ENDIF
c reset plas1 and nabs1
      plas1=zlas1
      nabs1=iabs1

c :::: ad1 restrict timestep by  applied laser power. This is not needed
c for the the very high powers where the absorption is small
c     zdt4=dt4
c     zlas1=plas1
c     tdash=time+dt4
c     CALL laser(tdash)
c     IF(plas1.gt.underf)THEN
c      zd7=0.
c      IF(pv.gt.underf)THEN
c       zd7=(zlas1+plas1)*dt4/pv
c       dt4=min(dt4,ak7/zd7*dt4)
c      ENDIF
c      IF(abs(dt4-zdt4).gt.underf)THEN
c       nceldt=nlp1
c       ncondt=7
c      ENDIF
c     ENDIF
c     dt4=max(dtmin,min(dtmax,dt4))
c     put back old value
c     plas1=zlas1

      IF(abs(piq(58)-2.).lt.underf.and.dt4.gt.10*dtmin)THEN
c      :::: ad1 restrict timestep by variation in  energy in
c      each zone to a fraction of thermal energy
       zd8=1.E-7
       DO 250 l=1,nl
        zther=0.5*(ddte1(l)+ddte3(l))*(te1(l)+te3(l))
        zz8=abs(xl3(l)/2.+eixch2(l)/2.+brems3(l)/2.+efbx2(l)+efrx2(l))
        zz8=zz8*dt4/zther
        IF(zz8.gt.zd8)THEN
         zd8=zz8
         il8=l
        ENDIF
 250   CONTINUE
c       write(6,*)mstep,il8,xplas1(il8+1),zd8
       z8=ak8*dt4/zd8
       zdt4=dt4
       dt4=min(dt4,z8)
       IF(abs(dt4-zdt4).gt.underf)THEN
        nceldt=il8
        ncondt=8
       ENDIF
       dt4=max(dtmin,min(dtmax,dt4))
c:::: restrict timestep by rate of change of ionisation state
       zd9=1.E-7
       DO 300 l=1,nl
c allow 1.0*ak9 electrons approximately for dz*
        dzal=1.0
        zz9=abs(zst(l)-zstu(l))*dt4/dt2
        zz9=zz9/dzal
        IF(zz9.gt.zd9)THEN
         zd9=zz9
         il9=l
        ENDIF
 300   CONTINUE
       z9=ak9*dt4/zd9
c       write(6,*)il9,xplas1(il9+1),zd9

       zdt4=dt4
       dt4=min(dt4,z9)
       IF(abs(dt4-zdt4).gt.underf)THEN
        nceldt=il9
        ncondt=9
       ENDIF
       dt4=max(dtmin,min(dtmax,dt4))
      ENDIF
c
c 3. compare dt4 and dt2 (eq.77)
c if the new value dt4 deviates too much from the old value dt2 then
c the time centering is damaged as regards the calculation of u
 400  zratio=dt4/dt2
      zak0=1.0/ak0
      IF(zratio.ge.zak0)THEN
       IF(zratio.gt.2.0)dt4=2.0*dt2
       IF(dt4.ge.dtmin)RETURN
      ENDIF
c the required value of dt4 is now too small
c if we iterate and if we are at the last iteration then break and
c revers the calculation. if we do not iterate then break and revers
      IF(.not.(nlite.and.(nit.ge.nitmax)))THEN
       IF(nlite)RETURN
      ENDIF
c break
      break=.true.
      ncondt=7
      RETURN
      END


      SUBROUTINE hfcor(zsr,p,w,beta,t,cor)
c     ad1
c     subroutine to calculate high field correction for inverse
c     bremsstrahlung absorption when vosc >=< vthermal
c     see l. schlessinger and j. wright phys. rev. a v20 n5 1934-1945
c     p : laser power in w/m**2
c     w : laser wavelength in meters
c     beta: omegap**2/omega**2
c     t : temperature in  k
      cor=1.+(0.4334*p*w*w/t/6./sqrt(1.-beta))
      cor=1./cor**1.5
c langdon correction
c     flang=1.
c     alf=zsr*0.4334*p*w*w/t
c     if(alf.gt.0.001)flang= 1.-0.553/(1.+(0.25/alf)**0.75)
c     cor=min(cor,flang)
      RETURN
      END


      SUBROUTINE absorb
      PARAMETER(kk=1001,underf=1.E-30)
c 2.9 evaluates the absorption of laser light
c c1.9 development and diagnostic parameters
      DIMENSION nadump(20),npdump(20),nvdump(20)
      LOGICAL nlched,nlhead(9),nlomt1(50),nlomt2(50),nlomt3(50),nlrept
      COMMON /comddp/ maxdum,mxdump,nadump,nclass,npdump,npoint,nsub,
     &                nvdump,nlched,nlhead,nlomt1,nlomt2,nlomt3,nlrept
c end of development and diagnostic parameters
c c2.1 hydrodynamics
      DIMENSION dm(kk),p3(kk),r1(kk),r3(kk),r5(kk),rho1(kk),rho3(kk),
     &          rho5(kk),u2(kk),u4(kk),v1(kk),v3(kk),v5(kk),lstat(kk)
      COMMON /comhyd/ dm,p3,pini,r1,r3,r5,rho1,rho3,rho5,rhoini,rhor,
     &                rini,time,u2,u4,uedge,v1,v3,v5
c end of hydrodynamics
c c2.2 thermodynamics
      DIMENSION ddroe1(kk),ddroe3(kk),ddroi1(kk),ddroi3(kk),ddte1(kk),
     &          ddte3(kk),ddti1(kk),ddti3(kk),xappae(kk),xappai(kk),
     &          pe1(kk),pe3(kk),pi1(kk),pi3(kk),te1(kk),te3(kk),ti1(kk),
     &          ti3(kk)
      COMMON /comth / ddroe1,ddroe3,ddroi1,ddroi3,ddte1,ddte3,ddti1,
     &                ddti3,gammae,gammai,xappae,xappai,pe1,pe3,pi1,pi3,
     &                te1,te3,teini,ti1,ti3,tiini
c end of thermodynamics
c c2.4 laser variables
      DIMENSION alpha1(kk),xaser1(kk),xl1(kk),xl3(kk)
      COMMON /comlas/ alpha1,elas1,xamda1,xaser1,xecri1,plas1,rabs1,
     &                rocri1,xl1,xl3,nabs1
c end of laser variables
c c2.5 thermonuclear reactions
      DIMENSION f1d(kk),f1h(kk),f1he3(kk),f1he4(kk),f1neu(kk),f1ntrl(kk)
     &          ,f1t(kk),f1x(kk),f1x1(kk),f1x2(kk),f1x3(kk),f1x4(kk),
     &          f3d(kk),f3h(kk),f3he3(kk),f3he4(kk),f3neu(kk),f3ntrl(kk)
     &          ,f3t(kk),f3x(kk),f3x1(kk),f3x2(kk),f3x3(kk),f3x4(kk),
     &          r1dd(kk),r1dhe3(kk),r1dt(kk),r3dd(kk),r3dhe3(kk),
     &          r3dt(kk),ye1(kk),ye3(kk),yi1(kk),yi3(kk)
      COMMON /comfus/ deuter,f1d,f1h,f1he3,f1he4,f1neu,f1ntrl,f1t,f1x,
     &                f1x1,f1x2,f1x3,f1x4,f3d,f3h,f3he3,f3he4,f3neu,
     &                f3ntrl,f3t,f3x,f3x1,f3x2,f3x3,f3x4,heliu3,heliu4,
     &                hydrog,xetral,xtrlms,pneut1,pneut3,r1dd,r1dhe3,
     &                r1dt,r3dd,r3dhe3,r3dt,rneut1,rneut3,tinucl,totneu,
     &                tritiu,xmass,xmass1,xmass2,xmass3,xmass4,xtra,
     &                xtra1,xtra2,xtra3,xtra4,xz,xz1,xz2,xz3,xz4,ye1,
     &                ye3,yi1,yi3,yield
c end of thermonuclear reactions
c c2.7 physics control
      LOGICAL mlbrms,mlecon,mlfuse,mlicon,mlx,nlabs,nlbrms,nlburn,
     &        nlcri1,nldepo,nlecon,nlfuse,nlicon,nlmove,nlpfe,nlpfi,nlx,
     &        nlhf,nltnl
      COMMON /comcon/ mstep,ncase,ngeom,mlbrms,mlecon,mlfuse,mlicon,mlx,
     &                nlabs,nlbrms,nlburn,nlcri1,nldepo,nlecon,nlfuse,
     &                nlicon,nlmove,nlpfe,nlpfi,nlx,nlhf,nltnl,pmin,
     &                rhomin,temin,timin,umin
c end of physics control
c c3.1 numerical control parameters
      LOGICAL break,nlgoon,nlite
      COMMON /comnc / nceldt,ncondt,nit,nitmax,break,nlgoon,nlite,ak0,
     &                ak1,ak2,ak3,ak4,ak5,ak6,ak7,ak8,ak9,bneum,deltat,
     &                dt2,dt3,dt4,dtemax,dtimax,dtmax,dtmin,dumax,rdt2,
     &                rdt3,rdt4
c end of numerical control parameters
c c3.2 mesh and numerical methods
      DIMENSION ae(kk),ai(kk),be(kk),bi(kk),ce(kk),ci(kk),de(kk),di(kk),
     &          e(kk),f(kk),ge(kk),gi(kk),q2(kk),q4(kk),teite(kk),
     &          tiite(kk),uite(kk)
      COMMON /comnum/ mesh,nj,njm1,njp1,nl,nlm1,nlp1,ae,ai,be,bi,ce,ci,
     &                de,di,e,f,ge,gi,q2,q4,teite,tiite,uite
c end of mesh and numerical methods
c c4.1 - administrative variables
      LOGICAL nldump,nlemp
      DIMENSION piq(302)
      COMMON /comadm/ maxdim,maxrun,ndump,nrep,nldump,nlemp,piq,tstop
c c4.1 - end administrative variables
c c2.3 ions and electrons
      DIMENSION brems1(kk),brems3(kk),ddz1(kk),ddz3(kk),degen(kk),
     &          effz(kk),eixch2(kk),fidash(kk),fz1(kk),fz3(kk),fzsq1(kk)
     &          ,fzsq3(kk),xlc(kk),xmieff(kk),xne(kk),xni(kk),omega1(kk)
      COMMON /comie / brems1,brems3,ddz1,ddz3,degen,degmax,degmin,effz,
     &                eixch2,fidash,fz1,fz3,fzsq1,fzsq3,xlc,xmieff,xne,
     &                xni,omega1,pmass
c end of ions and electrons
      COMMON /radflg/ yamda,xhnu,xihnu,dihnu,xplas1(kk),xrbf(kk),
     &                atrbf(kk),corhf1(kk),corhf(kk),iradf1,iradf2
      DIMENSION zneg(kk),zteg(kk),xa(kk),ya(kk)
      DIMENSION ze(10),zn(10)
      DATA ze,zn/0.2000,0.6000,1.0000,1.4000,1.8000,2.2500,2.7500,
     &     3.5000,4.5000,6.0000,0.0615,0.1297,0.1461,0.1377,0.1190,
     &     0.1187,0.0882,0.1075,0.0512,0.0404/
      DATA iclass,isub/2,9/
      IF(nlomt2(isub))RETURN
c     1. simplified absorption model
      CALL resetr(xl3,nl,0.0)
      CALL resetr(xplas1,nlp1,0.0)
c     if(plas1.lt.underf) return

      IF(nlabs)THEN
c      2. absorption by inverse bremsstrahlung in subcritical zones
c      clear the array xl3 if iterations are on
       IF(nlite)CALL resetr(xl3,nl,0.0)
       IF(nlite)CALL resetr(alpha1,nl,0.0)
c      2.1 the absorption coefficient (eq.26)
c      charge numbers
       zc2d=1.0
       zc2h=1.0
       zc2he3=2.0
       zc2he4=2.0
       zc2t=1.0
       zc2x=xz
       zloge=2.3026
       zc1=13.49/(xamda1*xamda1)*(1.0+piq(19))
       za1=5.05+zloge*log10(xamda1)
       l=nlp1
 50    l=l-1
       IF(l.ne.0)THEN
        zbeta1=xne(l)/xecri1
        IF(zbeta1.lt.1.00)THEN
         z2=fzsq3(l)/fz3(l)
         zb1=z2*zbeta1*zbeta1/(sqrt(1.0-zbeta1))
         zpmin1=1.67E-3*fz3(l)/te3(l)
         zpmin=2.97E-6/sqrt(te3(l))
         IF(zpmin1.gt.zpmin)zpmin=zpmin1
         zfac1=log(2.067E-4*xamda1*sqrt(te3(l))/zpmin)
         IF(zfac1.lt.1.0)zfac1=1.0
         alpha1(l)=zc1*zb1*zfac1/(te3(l)*sqrt(te3(l)))
c       ff absorption coef
c       theta=te3(l)/11605000.
c       srth = sqrt(theta)
c       rho=0.001*rho3(l)
c       dnion = (6.022e23*rho)/xmieff(l)
c       dnele= fz3(l)*dnion
c       zc=7.7e-46*dnele*dnion*fz3(l)*fz3(l)
c       un=xhnu/theta
c       eun = 0.0
c       if ( un.le.170.0 ) eun = exp(-un)
c       zrff=zc*(1.-eun)/(un*un*un*theta*theta*theta*srth)
c       zrff=zrff/sqrt(1.-zbeta1)
c       alpha1(l)=zrff
         IF(l.gt.0)GOTO 50
        ENDIF
       ENDIF
       nabs1=l
       rabs1=r3(nabs1+1)
       IF(nabs1.ne.0.and.nabs1.ne.nl)THEN
c      interpolation for radius of critical density
        nri=nl-nabs1
        nle=nabs1
        IF(nri.ge.2.and.nle.ge.2)THEN
         DO 60 k=nabs1-1,nabs1+2
          ki=k+1-(nabs1-1)
          xa(ki)=xne(k)
          ya(ki)=0.5*(r3(k+1)+r3(k))
 60      CONTINUE
        ELSEIF(nri.lt.2.and.nl.ge.4)THEN
         DO 70 k=nl-3,nl
          ki=k+1-(nl-3)
          xa(ki)=xne(k)
          ya(ki)=0.5*(r3(k+1)+r3(k))
 70      CONTINUE
        ELSEIF(nle.lt.2.and.nl.ge.4)THEN
         DO 80 k=1,4
          xa(k)=xne(k)
          ya(k)=0.5*(r3(k+1)+r3(k))
 80      CONTINUE
        ELSE
         WRITE(6,*)' less than 4 points for interpolation in absorb'
         STOP 1200
        ENDIF
        CALL polint(xa,ya,4,xecri1,rabs1,drabs1)
        IF(nabs1.lt.nl.and.rabs1.gt.r3(nabs1+2))THEN
         rabs1=r3(nabs1+1)
        ELSEIF(rabs1.gt.r3(nl+1))THEN
         rabs1=r3(nl+1)
        ELSEIF(rabs1.lt.r3(nabs1))THEN
         rabs1=r3(nabs1+1)
        ENDIF
       ENDIF
c      2.2 the rate of absorption
       zplas1=plas1
       xplas1(nlp1)=plas1
       l=nlp1
 100   l=l-1
       j=l
c      first for subcritical zones
       IF(l.gt.nabs1)THEN
c       the attenuation factor of the laser power through absorption
c       (if the rate is below truncation level we omit it)
c       high field  correction to absorption nlhf=.true.
        corhf(l)=1.0
        IF(nlhf)THEN
c        recalculate absorption coef
         alpha1(l)=0.
         zbeta1=xne(l)/xecri1
         IF(zbeta1.lt.1.00)THEN
          CALL hfcor(fz3(l),zplas1,xamda1,zbeta1,te3(l),corhf(l))
          z2=fzsq3(l)/fz3(l)
          zb1=z2*zbeta1*zbeta1/(sqrt(1.0-zbeta1))
          zpmin1=1.67E-3*fz3(l)/te3(l)
          zpmin=2.97E-6/sqrt(te3(l))
          IF(zpmin1.gt.zpmin)zpmin=zpmin1
          zfac1=log(2.067E-4*xamda1*sqrt(te3(l))/zpmin)
          IF(zfac1.lt.1.0)zfac1=1.0
          alpha1(l)=zc1*zb1*zfac1/(te3(l)*sqrt(te3(l)))*corhf(l)
         ENDIF
        ENDIF
        alphtt=alpha1(l)+xrbf(l)+atrbf(l)
        zal1dr=-alphtt*(r3(j+1)-r3(j))
c       take into account exact position for zone next to critical
        IF(l.eq.(nabs1+1).and.rabs1.gt.r3(j).and.rabs1.lt.r3(j+1))THEN
         alphtt=alpha1(l)+xrbf(l)+atrbf(l)
         zal1dr=-alphtt*(r3(j+1)-rabs1)
        ENDIF
        IF(abs(zal1dr).ge.1.0E-7)THEN
         IF(abs(zal1dr).gt.20.0)zal1dr=-20.0
         zabs1=1.0-exp(zal1dr)
c        to improve the accuracy we taylor expand the
c        absorp factor if zal1dr is small
         IF(abs(zal1dr).lt.1.0E-4)zabs1=-zal1dr
         zdpl1=zplas1*zabs1
c        rate of absorption in cell l (eq.28)
         ratio=alpha1(l)/(alpha1(l)+xrbf(l)+atrbf(l))
         xl3(l)=ratio*zdpl1/dm(l)
c        the attenuation of the laser beam (eq.27)
c        subtract amount of power that is absorbed in cell  l
         zplas1=zplas1-zdpl1
         xplas1(j)=zplas1
        ENDIF
        GOTO 100
       ENDIF
c      if critical dens in cell nabs1 then dump some energy there
       IF(nabs1.ne.0)THEN
c      if critical surface inside nabs1 zone add inv brems in nabs1
        IF(rabs1.lt.r3(nabs1+1))THEN
         l=nabs1
         j=l
         IF(nabs1.eq.nl)THEN
          zbeta1=0.98
         ELSE
          zbeta1=(xne(nabs1+1)+xecri1)/2./xecri1
         ENDIF
         z2=fzsq3(l)/fz3(l)
         zb1=z2*zbeta1*zbeta1/(sqrt(1.0-zbeta1))
         zpmin1=1.67E-3*fz3(l)/te3(l)
         zpmin=2.97E-6/sqrt(te3(l))
         IF(zpmin1.gt.zpmin)zpmin=zpmin1
         zfac1=log(2.067E-4*xamda1*sqrt(te3(l))/zpmin)
         IF(zfac1.lt.1.0)zfac1=1.0
         alpha1(l)=zc1*zb1*zfac1/(te3(l)*sqrt(te3(l)))
         alphtt=alpha1(l)+xrbf(l)+atrbf(l)
         zal1dr=-alphtt*(r3(j+1)-rabs1)
         IF(abs(zal1dr).ge.1.0E-7)THEN
          IF(abs(zal1dr).gt.20.0)zal1dr=-20.0
          zabs1=1.0-exp(zal1dr)
          IF(abs(zal1dr).lt.1.0E-4)zabs1=-zal1dr
          zdpl1=zplas1*zabs1
          ratio=alpha1(l)/(alpha1(l)+xrbf(l)+atrbf(l))
          xl3(nabs1)=ratio*zdpl1/dm(nabs1)
          zplas1=zplas1-zdpl1
          xplas1(j)=zplas1
         ENDIF
        ENDIF
       ENDIF
c     else if .not. nlabs
      ELSE
       IF(.not.nlcri1)RETURN
c      find critical density
       l=nlp1
 150   l=l-1
       IF(xne(l).lt.xecri1)THEN
        IF(l.gt.1)GOTO 150
       ENDIF
       nabs1=l
       zplas1=plas1
       xplas1(nabs1+1)=plas1
      ENDIF
c     endif nlabs
c     *************************************************************
c     this part of the code attempts to simulate heat transport
c     by fast electrons
c     the hot electron temperature may be calculated as a
c     function of (i lamda squared) or as a multiple of tcold
c     if nlcri1 and hot electron
      IF(abs(piq(21)).gt.underf)THEN
       IF(piq(22).lt.0.0)THEN
        zr=r3(nabs1)
        IF(zr.lt.1.0E-06)zr=1.0E-06
        igeom=ngeom-1
        zil=zplas1*xamda1*xamda1/(zr**igeom)
        IF(zil.gt.1.0E7)THEN
         zthot=-piq(22)*1.89E6*zil**0.25
        ELSE
         zthot=-piq(22)*2.28E3*zil**0.67
        ENDIF
       ELSE
        zthot=piq(22)*te3(nabs1)
       ENDIF
       IF(zthot.lt.10000.0)zthot=10000.0
c      we put zthot into common for printing out later
       xaser1(kk)=zthot
c      the absorbed laser power is divided into a fraction piq(21)
c      which goes into fast electrons  and  ( 1.0 - piq(21) ) which
c      goes into thermal electrons
c
c      a fraction  piq(52)  of the incident laser energy
c      is absorbed    this fraction is divided between thermal
c      and supra-thermal (i.e. fast ) electrons
       zplas1=piq(52)*zplas1
c      the fraction of laser energy into thermal electrons
       xl3(nabs1)=(1.0-piq(21))*zplas1/dm(nabs1)
c      the fast electrons are described by a
c      multi-group treatment
c      we use ten groups of fast electrons
c      the number and energy of each group being specified
c      by the arrays  zn  and  ze
c      the term 'number density' is used loosely , in fact
c      it is defined such that the energy content of each group
c      is given by 'ne *te'
c      we define number density and temperature for each group
       DO 250 igroup=1,10
        ztef=ze(igroup)*zthot
        znef=zn(igroup)*zplas1*piq(21)/ztef
c       the approximate coulomb logarithm for fast electrons
        zcoulo=5.0*(1.0+piq(53))
c       we follow the fast electrons inwards , diminishing znef
c       at each cell by a geometric factor to prevent excessive preheat
c       at the centre . the amount of energy ' lost ' in this way is
c       stored in arrays zneg and zteg for use later
        CALL resetr(zneg,nlp1,0.0)
        CALL resetr(zteg,nlp1,0.0)
c       although only half the fast electrons go inwards
c       this has the effect of making the absorbed power a function
c       of piq(21)   so for this version of absorb we assume that
c       all the fastelectrons are directed inwards
c       alternatively one could say that the outward going fast electrs
c       are turned back by electric fields in the plasma sheath
        DO 160 ll=1,nabs1
         l=nabs1-ll+1
c        the fractional energy deposited in cell  l  is calculated
c        from the energy loss equation for fast test particles
c        ' nrl plasma formulary p5 '
         zfe=3.5E-09*xni(l)*effz(l)*zcoulo*(r3(l+1)-r3(l))/(ztef*ztef)
         IF(zfe.gt.1.0)zfe=1.0
         zeabs=zfe*znef*ztef
         xl3(l)=xl3(l)+zeabs/dm(l)
c        the energy of the fast electrons is reduced
         ztef=ztef*(1.0-zfe)
         IF(ztef.lt.100.0)GOTO 180
c        the number of fast electrons is reduced by a geom factor
         IF(l.ne.1)THEN
          igeom=ngeom-1
          zfgeom=(r3(l)/r3(l+1))**igeom
          zneg(l)=znef*(1.0-zfgeom)
          zteg(l)=ztef
          znef=znef*zfgeom
         ENDIF
 160    CONTINUE
c       having followed the fast electrons into the center
c       we now go back outwards and pick up the energy 'lost' in the
c       geometrical factor
 180    DO 200 l=2,nl
         IF(ztef.ge.100.0)THEN
          zfe=3.5E-09*xni(l)*effz(l)*zcoulo*(r3(l+1)-r3(l))/(ztef*ztef)
          IF(zfe.gt.1.0)zfe=1.0
          zeabs=zfe*znef*ztef
          xl3(l)=xl3(l)+zeabs/dm(l)
          ztef=ztef*(1.0-zfe)
         ENDIF
         zneft=znef+zneg(l)
         IF(abs(zneft).gt.underf)THEN
          ztef=(znef*ztef+zneg(l)*zteg(l))/zneft
          znef=zneft
         ENDIF
 200    CONTINUE
c       we wish to print out the energy not absorbed so we put
c       it into a common array
        igrpd=51+igroup
        IF(abs(zplas1).lt.underf)THEN
         xaser1(igrpd)=0.0
        ELSE
         xaser1(igrpd)=znef*ztef/(zn(igroup)*zplas1)
        ENDIF
 250   CONTINUE
c      any fast electron energy that has not been absorbed at
c      this stage is assumed to be distributed in proportion
c      to the cell masses
c      the energy that has not been absorbed
       zenabs=0.0
       DO 300 igroup=1,10
        igrpd=51+igroup
        zenabs=zenabs+xaser1(igrpd)*zn(igroup)*zplas1
 300   CONTINUE
c      the total target mass
       zmass=0.0
       DO 350 l=1,nl
        zmass=zmass+dm(l)
 350   CONTINUE
       DO 400 l=1,nl
        xl3(l)=xl3(l)+piq(57)*zenabs/zmass
 400   CONTINUE
c in case of fast electrons return
       RETURN
      ENDIF
c     ***************************************************************
      IF(.not.nlcri1)RETURN
c     else if nlcri and no hot electrons
c     dump laser power in cell nabs1
      IF(nabs1.le.1)return
      xl3(nabs1)=xl3(nabs1)+piq(52)*zplas1/dm(nabs1)
      RETURN
      END


      SUBROUTINE setup(zm,am,rh,te,dt2,mesh,icxrl,istage,rb,uu,nlres)
      PARAMETER(kk=1001,underf=1.E-30)
      PARAMETER(ntrs=35)
      LOGICAL nlres
      COMMON /atom  / pu(10,kk),p(10,kk),pz(10,kk),zs(10,kk),zst(kk),
     &                zstz(kk),z(kk),a(kk),fn(10),deg(10),qs(10,kk),
     &                zstmi(kk),zstu(kk),engu(10,kk),eng(10,kk),
     &                rdrhe(kk),zstin(kk),nmax
      COMMON /atomdt/ foss(10,10),sig(10,10)
      COMMON /alph  / rwi(kk,3),frinv(kk,3),alpha(kk,3),ideg(10),
     &                idegm(10)
      COMMON /opdep / tau(10,10)
      COMMON /grap  / ipt
      COMMON /npist / npp1(10,10),npp2(10,10),npp3(10,10),npp4(10,10),
     &                npp5(10,10)
      COMMON /npp   / npp(10,10),il,iu
      COMMON /netint/ taxint(ntrs,6),tgain(ntrs)

      DIMENSION zm(kk),am(kk),rh(kk),te(kk),rb(kk),uu(kk)
      COMMON /dpi   / ecbb(kk),ecbf(kk),eraddr(kk),erbf(kk),eati(kk),
     &                dcvz(kk),efbx2(kk),efrx2(kk),erbb(kk),erfb(kk),
     &                acfc(10,kk),acfr(10,kk),bcfc(10,kk),bcfr(10,kk),
     &                edfb(kk),ter(kk),dqdtim(kk),nmaxr
      ifive=5
      iten=10
c setup default maximum principal quantum number
      nmax=ifive+1
c setup different populations for h-like, li-like and na-like schemes
      IF(istage.eq.1)THEN
c h-like scheme istage=1
       DO 50 i=1,iten
        DO 20 j=1,iten
         npp(i,j)=npp1(i,j)
 20     CONTINUE
 50    CONTINUE
       il=2
       iu=3
       nmax=6
      ELSEIF(istage.eq.2)THEN
c li-like scheme istage=2
       DO 100 i=1,iten
        DO 60 j=1,iten
         npp(i,j)=npp2(i,j)
 60     CONTINUE
 100   CONTINUE
       il=3
       iu=4
       nmax=6
      ELSEIF(istage.eq.3)THEN
c na-like scheme istage=3
       DO 150 i=1,iten
        DO 120 j=1,iten
         npp(i,j)=npp3(i,j)
 120    CONTINUE
 150   CONTINUE
       il=4
       iu=5
       nmax=7
      ELSEIF(istage.eq.4)THEN
c ne-like scheme istage=4
       DO 200 i=1,iten
        DO 160 j=1,iten
         npp(i,j)=npp4(i,j)
 160    CONTINUE
 200   CONTINUE
c        write(6,*)'setup',npp
       il=4
       iu=5
       nmax=8
      ELSEIF(istage.eq.5)THEN
c ni-like scheme   istage=5
       DO 250 i=1,iten
        DO 220 j=1,iten
         npp(i,j)=npp5(i,j)
 220    CONTINUE
 250   CONTINUE
       il=4
       iu=5
       nmax=9
      ELSE
       WRITE(6,10300)istage,nmax
c        stop 1800
      ENDIF
c     nmax=ifive+(istage+1)/2
c ad1 force nmax
      nmax=9
c     setup average z and a in each cell
      DO 300 l=1,mesh
       z(l)=zm(l)
       a(l)=am(l)
 300  CONTINUE
c     force   nmax to
c     nmax=10
c     setup oscillator strengths
      DO 400 i=1,nmax
       DO 350 j=1,nmax
        foss(i,j)=0.0
        IF(j.gt.i)THEN
         fi=real(i)
         fj=real(j)
         wa=(fi**5)*(fj**3)
         wa=1.0/wa
         wb=1.0/(fi*fi)
         wc=1.0/(fj*fj)
         wd=wb-wc
         wd=1.0/(wd*wd*wd)
         foss(i,j)=1.96*wa*wd
c set optical depth to zero
         tau(i,j)=0.0
        ENDIF
 350   CONTINUE
 400  CONTINUE
      foss(1,2)=0.4161
      foss(1,3)=0.0792
      foss(1,4)=0.0290
      foss(1,5)=0.0139
      foss(2,3)=0.6370
      foss(2,4)=0.1190
      foss(2,5)=0.0443
      foss(3,4)=0.8408
      foss(3,5)=0.1499
      foss(4,5)=1.037
c     set graph counter to zero
      ipt=0
c setup integer values of principal quantum numbers and degeneracies
      DO 500 i=1,nmax
       fn(i)=real(i)
       deg(i)=2.0*fn(i)*fn(i)
       ideg(i)=2*i*i
       idegm(i)=ideg(i)-1
 500  CONTINUE
      IF(nlres)RETURN
c
c     setup initial populations
      DO 700 l=1,mesh
       rhoi=rh(l)
       thetai=te(l)
c      calculate approximation to the lte degree of ionisation
c      if(initial ionization is 0 calculate init ioniz
       if(zstin(l).lt.underf) then
         CALL zbar(rhoi,thetai,l,zstini)
         zstt=z(l)-zstini
c        calculate initial nlte populations by aufbau principle
c        using lte degree of ionisation
         ifg=0
         bel=0.0
         DO 520 j=1,nmax
           IF(ifg.gt.0)THEN
             p(j,l)=0.0
           ELSE
             belt=bel+deg(j)
             IF(belt.ge.zstt)THEN
               ifg=ifg+1
               p(j,l)=zstt-bel
             ELSE
               bel=belt
               p(j,l)=deg(j)
             ENDIF
           ENDIF
 520     CONTINUE
c ad1
         CALL pe(rhoi,thetai,l)
         DO 540 i=1,nmax
           pu(i,l)=p(i,l)
           engu(i,l)=eng(i,l)
 540     CONTINUE

c        does rate equations give reasonable initial zstar
         xtesti=1.E-6
         DO 560 it=1,100
           theta1=thetai+(81-it)*(0.1*thetai)
           IF(theta1.lt.thetai)theta1=thetai
           CALL rate(rhoi,theta1,nstep,time,1.00,l,rb,uu,xtesti)
           dif=0.0
           DO 550 i=1,nmax
             pu(i,l)=p(i,l)
             engu(i,l)=eng(i,l)
             diff=diff+abs(p(i,l)-pu(i,l))
 550       CONTINUE
           zstu(l)=zst(l)
 560     CONTINUE

c        if no convergence zstart = zstmi reset  to zstmi
         IF(abs(zst(l)-zstmi(l)).lt.underf)THEN
           IF(zstini.lt.zstmi(l))THEN
             zstini=zstmi(l)
           ELSEIF(zstini.gt.zstmi(l))THEN
             zstmi(l)=zstini
           ENDIF
         ENDIF
         zstt=z(l)-zstini
c recalculate initial nlte populations by aufbau principle using lte
c degree of ionisation
         ifg=0
         bel=0.0
         DO 601 j=1,nmax
           IF(ifg.gt.0)THEN
             p(j,l)=0.0
           ELSE
             belt=bel+deg(j)
             IF(belt.ge.zstt)THEN
               ifg=ifg+1
               p(j,l)=zstt-bel
             ELSE
               bel=belt
               p(j,l)=deg(j)
             ENDIF
           ENDIF
 601     CONTINUE
         CALL pe(rhoi,thetai,l)
         DO 651 i=1,nmax
           pu(i,l)=p(i,l)
           engu(i,l)=eng(i,l)
 651     CONTINUE
         zstu(l)=zst(l)
       else
c        else if zstin(l) set in input data (use it)
         zstini=zstin(l)
         zstt=z(l)-zstini
c recalculate initial nlte populations by aufbau principle using lte
c degree of ionisation
         ifg=0
         bel=0.0
         DO 600 j=1,nmax
           IF(ifg.gt.0)THEN
             p(j,l)=0.0
           ELSE
             belt=bel+deg(j)
             IF(belt.ge.zstt)THEN
               ifg=ifg+1
               p(j,l)=zstt-bel
             ELSE
               bel=belt
               p(j,l)=deg(j)
             ENDIF
           ENDIF
 600     CONTINUE
         CALL pe(rhoi,thetai,l)
         DO 650 i=1,nmax
           pu(i,l)=p(i,l)
           engu(i,l)=eng(i,l)
 650     CONTINUE
         zstu(l)=zst(l)
       endif
 700  CONTINUE
      ecbb(l)=0.
      ecbf(l)=0.0
      eraddr(l)=0.
      erbf(l)=0.
      eati(l)=0.
      dcvz(l)=0.0
      efbx2(l)=0.0
      efrx2(l)=0.0
      erbb(l)=0.0
      erfb(l)=0.0
      edfb(l)=0.0
      zstu(l)=zstmi(l)
      DO 800 j=1,ntrs
       DO 750 k=1,2
        taxint(j,k)=0.0
 750   CONTINUE
 800  CONTINUE
      RETURN
10100 FORMAT(1x,'###: can not find initial solution from rate ',
     &       'equations in cell ',i4,/10x,'using minimum  Z*=',f10.3,
     &       ' for Zstart')
10200 FORMAT(1x,'###: can not find initial solution from rate ',
     &       'equations in cell ',i4,/10x,'using equilibrium  Z*=',
     &       f10.3,' for Zstart and Zstmin ')
10300 FORMAT(/2x,'error in istage = ',i5,' using nmax = ',i5)
      END


      SUBROUTINE rate(rho,theta,nstep,time,tstep,l,rb,uu,xtest)
      PARAMETER(kk=1001,underf=1.E-30)
c     x-ray laser common blocks
      DIMENSION fl(2)
      COMMON /xrl   / ropmul,tpon,tpoff,rmpd,ztrtop,dlamda,fl,fwide,
     &                efbmul,drmul,bbtrap,fbtrap,icxrl,istage,idflg,
     &                itbflg,istflg,idrflg,ipuflg,nlp,nup,llp,lup,nst,
     &                nfl,nlmax,ifrsta,ilosta,ihista,igstat
      COMMON /atom  / pu(10,kk),p(10,kk),pz(10,kk),zs(10,kk),zst(kk),
     &                zstz(kk),z(kk),a(kk),fn(10),deg(10),qs(10,kk),
     &                zstmi(kk),zstu(kk),engu(10,kk),eng(10,kk),
     &                rdrhe(kk),zstin(kk),nmax
      COMMON /rats  / rrc(10,10),cbb(10,10),rbb(10,10),rd(10,10),cbf(10)
     &                ,cfb(10),rfb(10),rfbst(10),rbf(10),atir(10),
     &                atie(10),rdr(10),edrbb
      COMMON /ratt  / tcbf(kk),tcfb(kk),trfb(kk),trfbst(kk),trbf(kk),
     &                tatir(kk),trdr(kk)
      COMMON /cvoigt/ bcfl(10,kk),va(kk)
      COMMON /dpi   / ecbb(kk),ecbf(kk),eraddr(kk),erbf(kk),eati(kk),
     &                dcvz(kk),efbx2(kk),efrx2(kk),erbb(kk),erfb(kk),
     &                acfc(10,kk),acfr(10,kk),bcfc(10,kk),bcfr(10,kk),
     &                edfb(kk),ter(kk),dqdtim(kk),nmaxr
c     c3.2 mesh and numerical methods
c end of mesh and numerical methods
c c2.2 thermodynamics
      DIMENSION ddroe1(kk),ddroe3(kk),ddroi1(kk),ddroi3(kk),ddte1(kk),
     &          ddte3(kk),ddti1(kk),ddti3(kk),xappae(kk),xappai(kk),
     &          pe1(kk),pe3(kk),pi1(kk),pi3(kk),te1(kk),te3(kk),ti1(kk),
     &          ti3(kk)
      COMMON /comth / ddroe1,ddroe3,ddroi1,ddroi3,ddte1,ddte3,ddti1,
     &                ddti3,gammae,gammai,xappae,xappai,pe1,pe3,pi1,pi3,
     &                te1,te3,teini,ti1,ti3,tiini
c end of thermodynamics
c c2.3 ions and electrons
      DIMENSION brems1(kk),brems3(kk),ddz1(kk),ddz3(kk),degen(kk),
     &          effz(kk),eixch2(kk),fidash(kk),fz1(kk),fz3(kk),fzsq1(kk)
     &          ,fzsq3(kk),xlc(kk),xmieff(kk),xne(kk),xni(kk),omega1(kk)
      COMMON /comie / brems1,brems3,ddz1,ddz3,degen,degmax,degmin,effz,
     &                eixch2,fidash,fz1,fz3,fzsq1,fzsq3,xlc,xmieff,xne,
     &                xni,omega1,pmass
c end of ions and electrons
c c4.1 - administrative variables
c c4.1 - end administrative variables
c c3.1 numerical control parameters
c end of numerical control parameters
c c2.7 physics control
      LOGICAL mlbrms,mlecon,mlfuse,mlicon,mlx,nlabs,nlbrms,nlburn,
     &        nlcri1,nldepo,nlecon,nlfuse,nlicon,nlmove,nlpfe,nlpfi,nlx,
     &        nlhf,nltnl
      COMMON /comcon/ mstep,ncase,ngeom,mlbrms,mlecon,mlfuse,mlicon,mlx,
     &                nlabs,nlbrms,nlburn,nlcri1,nldepo,nlecon,nlfuse,
     &                nlicon,nlmove,nlpfe,nlpfi,nlx,nlhf,nltnl,pmin,
     &                rhomin,temin,timin,umin
c end of physics control

c c4.1 - administrative variables
      LOGICAL nldump,nlemp
      DIMENSION piq(302)
      COMMON /comadm/ maxdim,maxrun,ndump,nrep,nldump,nlemp,piq,tstop
c c4.1 - end administrative variables
      COMMON /npp   / npp(10,10),il,iu
      COMMON /radflg/ yamda,xhnu,xihnu,dihnu,xplas1(kk),xrbf(kk),
     &                atrbf(kk),corhf1(kk),corhf(kk),iradf1,iradf2

      DIMENSION pn(10),po(10),eo(10),ppr(10),rb(kk),uu(kk)
c updates pupulations by solving rate
c equations
      zk=1.38E-23
      an=6.0232E+23
      ic=0
      icmax=600
      alp=0.5
c     IF(icxrl.eq.0.or.istage.gt.3)xtest=1.E-04
      IF(zst(l).le.zstmi(l).or.zstu(l).le.zstmi(l))xtest=1.E-03
      dpito=1.E-06*xtest
cc
      dr=(rb(nfl+1)-rb(l))
      vel=abs(uu(nfl+1)-uu(l))
      veli=1.38E06*sqrt(8.61697E-05*ti3(l)/a(l))
      drmax=2.0*(rb(nfl+1)-rb(l))
      IF(abs(vel*veli).gt.underf)then
      IF((veli*dr/vel).gt.drmax)
     & dr=drmax*vel/veli
      endif
c     set old population to that obtained at previous timestep
c     and new populations to zero
c     ad1 po(i)
      zstu(l)=z(l)
      DO 100 i=1,nmax
       po(i)=pu(i,l)
       zstu(l)=zstu(l)-po(i)
       eo(i)=engu(i,l)
       pn(i)=0.0
 100  CONTINUE
c     if thermal band is required, alter the number of levels
c     for which the rate equations are solved
      nmaxr=nmax
      IF(itbflg.ne.0)nmaxr=nmax-1
c     iterate to solve rate equation
 200  CALL pe(rho,theta,l)
      IF(idrflg.ne.0)CALL dier(drmul,theta,rho,l)
c     skip calculation of optical depth if idflg = 0
      IF(idflg.ne.0.and.abs(ropmul).gt.underf)
     &   CALL depth(rho,l,dr,vel,istage,ropmul)
      CALL ratset(rho,theta,l,idflg,ipuflg,time,bbtrap,tpon,tpoff,nup,
     &            nlp,rmpd,istage,il,iu)
      diff=0.0
      xrbf(l)=0.
      atrbf(l)=0.
c     if(l.eq.50)write(6,*)ic,(p(i,l),i=1,nmax),zst(l)
      tcbf(l)=0.
      tcfb(l)=0.
      trbf(l)=0.
      trfb(l)=0.
      trfbst(l)=0.
      tatir(l)=0.
      trdr(l)=0.
      DO 300 i=1,nmaxr
       CALL coff(rho,theta,i,l,acf,bcf,acf2,bcf2,fbtrap,idrflg)
       acf=acf+acf2
       bcf=bcf+bcf2
       bcfl(i,l)=bcf
       pn(i)=(po(i)+(acf*tstep))/((((acf/deg(i))+bcf)*tstep)+1.0)

       diff=abs(pn(i)-p(i,l))+diff
c       if(l.eq.50 .and. nstep.gt.115)write(6,*)i,acf,bcf,pn(i)

       qi=1.-p(i,l)/deg(i)
       tcbf(l)=tcbf(l)+p(i,l)*cbf(i)
       tcfb(l)=tcfb(l)+zst(l)*cfb(i)*qi
       trbf(l)=trbf(l)+p(i,l)*rbf(i)
       trfb(l)=trfb(l)+zst(l)*rfb(i)*qi
       trfbst(l)=trfbst(l)+zst(l)*rfbst(i)*qi
       trdr(l)=trdr(l)+zst(l)*rdr(i)*qi
       tatir(l)=tatir(l)+p(i,l)*atir(i)

c      if(l.eq.nfl) then
c       write(6,*)'i,pi,acf,bcf',i,pn(i),acf,bcf
c      endif

 300  CONTINUE
c     if(l.eq.nfl)write(6,*)'itbflg',itbflg
c add change in last level if the thermal band is used
      IF(itbflg.ne.0)THEN
       pn(nmax)=deg(nmax)
       eth=eng(nmax,l)/theta
       IF(eth.le.170.0)THEN
        wa=exp(-eth)*theta*sqrt(theta)
        wb=(317.0*a(l))/(rho*zst(l))
        pn(nmax)=pn(nmax)/(1.0+(wa*wb))
       ENDIF
       diff=abs(pn(nmax)-p(nmax,l))+diff
      ENDIF
      ic=ic+1
c damp changes in populations and reset old population for next
c iteration
      zstp=zst(l)
      tot=0.0
      DO 400 i=1,nmaxr
       ppr(i)=p(i,l)
       p(i,l)=((1.0-alp)*p(i,l))+(alp*pn(i))
       tot=tot+p(i,l)
 400  CONTINUE
c no damping for last level if thermal band is being used
      IF(itbflg.ne.0)ppr(nmax)=p(nmax,l)
      IF(itbflg.ne.0)p(nmax,l)=pn(nmax)
      IF(itbflg.ne.0)tot=tot+p(nmax,l)
      IF(tot.gt.(z(l)-zstmi(l)))THEN
c ensure number of bound electrons is less than z-zstmi
       tot=0.
       ifg=0
       DO 450 i=1,nmax
        IF(ifg.ne.0)THEN
         p(i,l)=0.0
        ELSE
         tot=tot+p(i,l)
         IF(tot.ge.(z(l)-zstmi(l)))THEN
          p(i,l)=(z(l)-zstmi(l))-(tot-p(i,l))
          ifg=ifg+1
         ENDIF
        ENDIF
 450   CONTINUE
       zst(l)=zstmi(l)
c test for convergence in population changes
      ELSEIF(ic.gt.icmax)THEN
       WRITE(6,10100)nstep,zstu(l),l,diff
      ELSEIF(diff.gt.xtest)THEN
       GOTO 200
      ENDIF
c
c
c     ionisation energy:
c     ==================
      ecbb(l)=0.
      ecbf(l)=0.0
      eraddr(l)=0.
      erbf(l)=0.
      eati(l)=0.
      dcvz(l)=0.0
      dqdtim(l)=0.0
      efrx2(l)=0.0
      erbb(l)=0.0
      erfb(l)=0.0
      edfb(l)=0.0
c if no ionization energy
c      return
      if(.not.nlpfe.and.piq(59).gt.0.5) return
      if(ic.gt.icmax) return

      zk=1.38E-23
      an=6.0232E+23
      te3k=0.5*(te1(l)+te3(l))/11604000.
      edrbbk=edrbb
      xhnuk=xhnu
      IF(nstep.le.1)RETURN
      call pe(rho,theta,l)
      IF((zst(l)-zstmi(l)).gt.xtest.and.(zstu(l)-zstmi(l)).gt.xtest)
     & THEN
      DO 500 i=1,nmaxr
       qi=(1.-p(i,l)/deg(i))
       engik=eng(i,l)
c      tunnel ionisation and ati energy
       atiek=0.
       IF(nltnl)THEN
        IF(p(i,l).gt.xtest)atiek=atie(i)*1.e-3
        eati(l)=eati(l)+p(i,l)*atir(i)*atiek
       ENDIF
c      photoionisation corrected for stimulated emission.
       IF(xhnuk.ge.engik)erbf(l)=erbf(l)
     &                           +(p(i,l)*rbf(i)-zst(l)*rfbst(i)*qi)
     &                           *(xhnuk-engik)
 500  CONTINUE
c     additional energy input into f e from rads
      efrx2(l)=erbf(l)+eati(l)

c loss of free e- energy through coll processes
       ecbfp=0.0
       ecbfn=0.0
       ecbbp=0.0
       ecbbn=0.0
c xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
c     if(l.eq.50)write(6,*)'in ioni',(p(i,l),i=1,nmax),zst(l)
      DO 550 i=1,nmax
       qi=(1.-p(i,l)/deg(i))
c      keep rates  for use in rescaling
       acfc(i,l)=0.
       bcfc(i,l)=0.
       acfr(i,l)=0.
       bcfr(i,l)=0.
       DO 503 j=1,nmax
        IF (j.ne.i) THEN
         acfc(i,l)=acfc(i,l)+p(j,l)*(cbb(j,i)+rbb(j,i))
         qj=1.-p(j,l)/deg(j)
         bcfc(i,l)=bcfc(i,l)+(cbb(i,j)+rbb(i,j))*qj
        ENDIF
 503   CONTINUE
       acfc(i,l)=acfc(i,l)+zst(l)*(cfb(i)+rfb(i)+rdr(i))
       acfr(i,l)=acfr(i,l)+zst(l)*rfbst(i)
       bcfc(i,l)=bcfc(i,l)+cbf(i)
       bcfr(i,l)=bcfr(i,l)+rbf(i)+atir(i)
c      if(l.eq.nfl)then
c      write(6,*)'i,acfcr,bcfcr',i,acfc(i,l)+acfr(i,l),bcfc(i,l)+
c    & bcfr(i,l)
c      endif
        engik=eng(i,l)
c        now  calculate energies
c        arising from transitions from level i
        DO 520 j=i+1,nmax
c     collision de-excitation  /   excitation energies
         engjk=eng(j,l)
         qj=1.-p(j,l)/deg(j)
         delte=engik-engjk
         ecbbp=ecbbp+p(j,l)*cbb(j,i)*qi*delte
         ecbbn=ecbbn+p(i,l)*cbb(i,j)*qj*delte
 520    CONTINUE
c      3b recomb - collis ioniz
        ecbfp=ecbfp+zst(l)*cfb(i)*qi*engik
        ecbfn=ecbfn+p(i,l)*cbf(i)*engik
c     free-bound radiation
c        srth   = sqrt(theta)
c        cf     = (2.0*fn(i)*fn(i)*rho)/(317.0*a(l)*theta*srth)
c        wa     = 7.86864e+09
c        wb     = (zs(i,l)**4)/(fn(i)**5)
c        wc     = cf
c     erfb(l)   = erfb(l) - zst(l)*wa*wb*wc*qi*te3k*efbmul
c     assume energy of recomb electrons 3/2 kt
        rarfbi=zst(l)*rfb(i)*qi
        erfbn=rarfbi*(1.5*te3k*efbmul+engik)
        erfb(l)=erfb(l)-erfbn
c       bound- bound radiation (does not enter energy eq.)
        DO 540 j=i+1,nmax
         engjk=eng(j,l)
         erbbn=p(j,l)*rbb(j,i)*qi*(engik-engjk)
         erbb(l)=erbb(l)-erbbn
 540    CONTINUE
 550   CONTINUE

       ecbb(l)=(ecbbp-ecbbn)
c      add collis ion / 3b recomb to collis exc/dexcit
       ecbf(l)=(ecbfp-ecbfn)
c      rad and dr part contribution to q
c      add bb dier to bb radiation
       edfb(l)=-edrbbk
       erbb(l)=erbb(l)-edrbbk
       eraddr(l)=erfb(l)+edfb(l)
c      put all collisiona energy consumption in one term
       dqdtim(l)=ecbb(l)+ecbf(l)+eraddr(l)
c      change in cv  corresponding to dz*
       dzo=zst(l)-zstu(l)
       dcvz(l)=(ddte3(l)-ddte1(l))

c      convert energies to joules / unit mass
       denom=pmass*xmieff(l)/1.6e-16
       dqdtim(l)=dqdtim(l)/denom
       efrx2(l)=efrx2(l)/denom
       ecbb(l)=ecbb(l)/denom
       ecbf(l)=ecbf(l)/denom
       eraddr(l)=eraddr(l)/denom
       erbf(l)=erbf(l)/denom
       eati(l)=eati(l)/denom
c      dcvz(l)=dcvz(l)
       erbb(l)=erbb(l)/denom
       erfb(l)=erfb(l)/denom
       edfb(l)=edfb(l)/denom
      ENDIF
      RETURN
10100 FORMAT(1x,i5,1pe12.3,'populations failed to converge in cell',i4,
     &       'diff=',f15.9)
10200 FORMAT(1x,'populations failed to converge in cell',i4,'diff=',
     &       f15.9,' time=',1pe14.7)
      END


      SUBROUTINE dier(drmul,theta,rho,l)
      PARAMETER(kk=1001,underf=1.E-30)
      COMMON /atom  / pu(10,kk),p(10,kk),pz(10,kk),zs(10,kk),zst(kk),
     &                zstz(kk),z(kk),a(kk),fn(10),deg(10),qs(10,kk),
     &                zstmi(kk),zstu(kk),engu(10,kk),eng(10,kk),
     &                rdrhe(kk),zstin(kk),nmax
      COMMON /atomdt/ foss(10,10),sig(10,10)
      COMMON /rats  / rrc(10,10),cbb(10,10),rbb(10,10),rd(10,10),cbf(10)
     &                ,cfb(10),rfb(10),rfbst(10),rbf(10),atir(10),
     &                atie(10),rdr(10),edrbb
      DIMENSION drfac(10)
c calculate dielectronic rates for he-like and lower ions
      accy=1.0E-06
      srthet=sqrt(theta)
      dnele=(6.022E23*zst(l)*rho)/a(l)
c zero dielectronic excitation array
      nmaxm=nmax-1
      DO 100 i=1,nmax
       drfac(i)=0.0
       rdr(i)=0.0
       DO 50 j=1,nmax
        rd(i,j)=0.0
 50    CONTINUE
 100  CONTINUE
      edrbb=0.0
      IF(zst(l).gt.(z(l)-2.0))RETURN
c dielectronic upwards rate
      DO 200 i=1,nmaxm
       ip=i+1
       tn=4.77E18/dnele*zs(i,l)**6*srthet
       tn=tn**(1./7.)
       nt=nint(tn)
       nt=min0(nmax,ip+nt)
       DO 150 j=ip,nt
        une=0.0136*(zs(i,l)**2/real(i)**2-zs(j,l)**2/real(j)**2)
        az=1.+0.015*zs(i,l)**3/(zs(i,l)+1.)**2
        un=une/(az*theta)
        eun=0.0
        IF(un.le.170.0)THEN
         eun=exp(-un)
c dielectronic upwards rate
         bz=sqrt(zs(i,l)/(zs(i,l)**2+13.4))
         bz=bz*(zs(i,l)+1.0)**2.5
         y=(zs(i,l)+1.0)/(1./real(i)**2-1./real(j)**2)
         ay=0.5*sqrt(y)/(1.0+(0.210*y)+(0.030*y*y))
         dzt=0.0015*(tn*(zs(i,l)+1.))**2
         dzt=dzt/(1.+dzt)
         xa=7.59E-14*dnele*foss(i,j)*bz*dzt*ay*eun
         xb=theta*srthet
         rd(i,j)=xa/xb*drmul
        ENDIF
 150   CONTINUE
 200  CONTINUE
c sum excitation rates to provide dielectronic recombination rates
      rdrtot=0.0
      DO 300 i=1,nmaxm
       ip=i+1
       DO 250 j=ip,nmax
        za=p(i,l)*(1.0-p(j,l)/deg(j))*rd(i,j)
        rdr(i)=rdr(i)+za
c energy emitted as bb radiation following down transition
c of excited electron
        une=0.0136*zs(i,l)**2*(1./real(i)**2-1./real(j)**2)
        edrbb=edrbb+2.*une*za
 250   CONTINUE
       rdrtot=rdrtot+rdr(i)
 300  CONTINUE
      IF(abs(rdrtot).lt.underf)RETURN
c calculate fraction of dielectronic recombination in each shell
      idflg=0
      qtot=0.0
      DO 400 j=1,nmax
       za=1.0-p(j,l)/deg(j)
       qtot=qtot+za
       IF(qtot.le.1.0)THEN
        drfac(j)=1.0
       ELSEIF(idflg.eq.0)THEN
        idflg=idflg+1
        drfac(j)=(1.0-(qtot-za))/za
       ENDIF
 400  CONTINUE
c adjust dielectronic recombination rates into each shell
      rdrt=0.0
      DO 500 j=1,nmax
       IF(drfac(j).le.1.0)THEN
        rdr(j)=drfac(j)*rdrtot
        za=1.0-p(j,l)/deg(j)
        rdrt=rdrt+(rdr(j)*za)
       ELSE
        WRITE(6,10100)j,drfac(j)
        STOP 2000
       ENDIF
       rdr(j)=rdr(j)/zst(l)
 500  CONTINUE
      rdtest=(rdrt-rdrtot)/rdrtot
      IF(abs(rdtest).ge.accy)THEN
       WRITE(6,10200)rdrt,rdrtot
       STOP 2200
      ENDIF
      RETURN
10100 FORMAT(4x,'error 1 in dr ',i4,1pe15.4)
10200 FORMAT(2x,'error 2 in dr ',2(1pe15.4))
      END


      SUBROUTINE pe(rho,theta,l)
      PARAMETER(kk=1001)
      COMMON /atom  / pu(10,kk),p(10,kk),pz(10,kk),zs(10,kk),zst(kk),
     &                zstz(kk),z(kk),a(kk),fn(10),deg(10),qs(10,kk),
     &                zstmi(kk),zstu(kk),engu(10,kk),eng(10,kk),
     &                rdrhe(kk),zstin(kk),nmax
      COMMON /atomdt/ foss(10,10),sig(10,10)
c this subroutine calculates the non-lte screened nuclear charges
c and energy levels given the populations
c calculate the inner screened nuclear charges
      DO 100 i=1,nmax
       im=i-1
       ip=i+1
       qs(i,l)=z(l)
       IF(i.ne.1)THEN
        DO 20 j=1,im
         qs(i,l)=qs(i,l)-(sig(i,j)*p(j,l))
 20     CONTINUE
       ENDIF
       zs(i,l)=qs(i,l)
c setup total screened nuclear charges
       DO 50 j=ip,nmax
        zs(i,l)=zs(i,l)-(sig(i,j)*p(j,l))
 50    CONTINUE
c
c add contributions to both screened charges from the
c same shell
c
       wa=sig(i,i)*p(i,l)
       wb=0.5*(1.0-(1.0/deg(i)))*wa
       qs(i,l)=qs(i,l)-wb
       zs(i,l)=zs(i,l)-2.*wb
 100  CONTINUE
c calculate the ionisation energy from each shell
      al=0.0072971
      DO 200 i=1,nmax
       wa=0.0068*((zs(i,l)/fn(i))**2)
       wb=((al*zs(i,l))**2)/fn(i)
       wc=1.0-(0.75/fn(i))
       wd=0.25/fn(i)
       eng(i,l)=wa*((1.0+(wb*wc))+(1.0+(wb*wd)))
 200  CONTINUE
c recalculate the total number of free electrons
      zst(l)=z(l)
      DO 300 i=1,nmax
       zst(l)=zst(l)-p(i,l)
 300  CONTINUE
      RETURN
      END


      SUBROUTINE ratset(rho,theta,l,idflg,ipuflg,time,bbtrap,tpon,tpoff,
     &                  nup,nlp,rmpd,istage,il,iu)
      PARAMETER(kk=1001,underf=1.E-30)
      COMMON /atom  / pu(10,kk),p(10,kk),pz(10,kk),zs(10,kk),zst(kk),
     &                zstz(kk),z(kk),a(kk),fn(10),deg(10),qs(10,kk),
     &                zstmi(kk),zstu(kk),engu(10,kk),eng(10,kk),
     &                rdrhe(kk),zstin(kk),nmax
      COMMON /atomdt/ foss(10,10),sig(10,10)
      COMMON /rats  / rrc(10,10),cbb(10,10),rbb(10,10),rd(10,10),cbf(10)
     &                ,cfb(10),rfb(10),rfbst(10),rbf(10),atir(10),
     &                atie(10),rdr(10),edrbb
      COMMON /opdep / tau(10,10)

      COMMON /radflg/ yamda,xhnu,xihnu,dihnu,xplas1(kk),xrbf(kk),
     &                atrbf(kk),corhf1(kk),corhf(kk),iradf1,iradf2

      srth=sqrt(theta)
      feln=(zst(l)*rho)/a(l)
c this subroutine calculates the bound-bound collisional and
c radiative rates
      nmaxm=nmax-1
      DO 100 i=1,nmax
       cbf(i)=0.0
       cfb(i)=0.0
       rfb(i)=0.0
       rfbst(i)=0.0
       rbf(i)=0.0
       atir(i)=0.
       atie(i)=0.
       DO 50 j=1,nmax
        rrc(i,j)=0.0
        cbb(i,j)=0.0
        rbb(i,j)=0.0
 50    CONTINUE
 100  CONTINUE
c ad1 allan
c        b-b transitions (should be 0) laser on and nltnl
c     if( (xplas1(l+1).gt.-1.e00)) return
      DO 200 i=1,nmaxm
       ip=i+1
       DO 150 j=ip,nmax
        une=eng(i,l)-eng(j,l)
        un=une/theta
        eun=0.0
c           write(6,*)i,eng(i,l),j,eng(j,l),une,theta,un
        IF(un.le.170.0)eun=exp(-un)
c calculate the correction to the collisional excitation rate gfac
        dn=real(i)
        dm=real(j)
        za=(1.0-(2.0/z(l)))*un
        za=1.0-za
        zb=(dm*(dm-dn))/20.0
        zc=1.0+(za*zb)
        ze=eintm(un)
        zc=0.9*zc*ze
        zd=1.0+zc
        gfac=0.19*zd
        IF(gfac.le.0.0)gfac=0.0
        wa=3.01E+14*foss(i,j)*feln*eun*gfac
        wb=srth*une
c[
        rrc(i,j)=wa/wb
        cbb(i,j)=rrc(i,j)
        wc=3.01E+14*foss(i,j)*feln*(deg(i)/deg(j))*gfac
        IF(abs(eun).lt.underf)wc=0.0
        wd=srth*une
c[
        rrc(j,i)=wc/wd
        cbb(j,i)=rrc(j,i)
 150   CONTINUE
 200  CONTINUE
      DO 300 i=1,nmaxm
       ip=i+1
       DO 250 j=ip,nmax
c add radiative (stimulated and spontaneous) rates
        une=eng(i,l)-eng(j,l)
c
       if(iradf1.eq.0) then
c
        rr=4.34327E+13*foss(i,j)*une*une*((fn(i)/fn(j))**2)
        rr=rr*(1.-bbtrap)
c
c         include pumping on the nlp -> nup
c         transition for specified period if
c         ipuflg =/ 0
        IF (istage.le.3.and.ipuflg.ne.0) THEN
         IF (i.eq.nlp.and.j.eq.nup) THEN
          IF (time.gt.tpon.and.time.lt.tpoff) THEN
           rr=4.34327E+13*foss(i,j)*une*une*((fn(i)/fn(j))**2)
           rr=rr*(1.0+rmpd)
           rrc(j,i)=rrc(j,i)+rr
           rbb(j,i)=rr
           rr1=4.34327E+13*foss(i,j)*une*une*rmpd
           rrc(i,j)=rrc(i,j)+rr1
           rbb(i,j)=rr1
           IF (abs(time-timeol).gt.underf) WRITE (6,11100) time,i,j,rmpd
           timeol=time
           GOTO 250
          ENDIF
         ENDIF
        ENDIF

c       insert optical depth factor for resonance lines
c       skip optical depth correction if idflg = 0
        IF (istage.le.3.and.idflg.ne.0) THEN
         ig=il-1
         IF (i.eq.ig) THEN
          rrmult=hg(tau(j,i))
          rr=4.34327E+13*foss(i,j)*une*une*((fn(i)/fn(j))**2)
          rr=rr*rrmult
c        if(l.eq.210)write(6,*)' ratset',j,i,tau(j,i),rrmult
         ENDIF
        ENDIF
c
        rrc(j,i)=rrc(j,i)+rr
        rbb(j,i)=rr
11100 FORMAT (/2x,'pump on time',1pe15.4,' levels',2I4,' rmpd',1pe20.9)
c
       else
c
        un=une/theta
        radtb=0.0
        IF(iradf1.ne.0)CALL radb(un,theta,radtb)
        rr=4.34327E+13*foss(i,j)*une*une*((fn(i)/fn(j))**2)
        rr=rr*(1.+radtb)
c[
        rrc(j,i)=rrc(j,i)+rr
        rbb(j,i)=rr
        rr1=4.34327E13*foss(i,j)*une*une*radtb
        rbb(i,j)=rr1
c[
        rrc(i,j)=rrc(i,j)+rr1
       endif
 250   CONTINUE
 300  CONTINUE
      RETURN
10100 FORMAT(/2x,'pump on time',1pe15.4,' levels',2I4,' rmpd',1pe20.9)
      END


      SUBROUTINE coff(rho,theta,i,l,acf,bcf,acf2,bcf2,fbtrap,idrflg)
      PARAMETER(kk=1001,underf=1.E-30)
      COMMON /atom  / pu(10,kk),p(10,kk),pz(10,kk),zs(10,kk),zst(kk),
     &                zstz(kk),z(kk),a(kk),fn(10),deg(10),qs(10,kk),
     &                zstmi(kk),zstu(kk),engu(10,kk),eng(10,kk),
     &                rdrhe(kk),zstin(kk),nmax
      COMMON /rats  / rrc(10,10),cbb(10,10),rbb(10,10),rd(10,10),cbf(10)
     &                ,cfb(10),rfb(10),rfbst(10),rbf(10),atir(10),
     &                atie(10),rdr(10),edrbb
      COMMON /radflg/ yamda,xhnu,xihnu,dihnu,xplas1(kk),xrbf(kk),
     &                atrbf(kk),corhf1(kk),corhf(kk),iradf1,iradf2

c     physics control
      LOGICAL mlbrms,mlecon,mlfuse,mlicon,mlx,nlabs,nlbrms,nlburn,
     &        nlcri1,nldepo,nlecon,nlfuse,nlicon,nlmove,nlpfe,nlpfi,nlx,
     &        nlhf,nltnl
      COMMON /comcon/ mstep,ncase,ngeom,mlbrms,mlecon,mlfuse,mlicon,mlx,
     &                nlabs,nlbrms,nlburn,nlcri1,nldepo,nlecon,nlfuse,
     &                nlicon,nlmove,nlpfe,nlpfi,nlx,nlhf,nltnl,pmin,
     &                rhomin,temin,timin,umin
c end of physics control

c this subroutine sets up the coefficients in the rate equation
c acf and bcf
      acf=0.0
      bcf=0.
c calculate coll+  rad  excitation and de-excitation contribution
      xc=0.0
      xd=0.0
      IF(i.ne.1)THEN
       im=i-1
       DO 50 k=1,im
        xc=xc+(p(k,l)*rrc(k,i))
        facq1=(1.0-(p(k,l)/deg(k)))*rrc(i,k)
        IF(facq1.lt.0.0)facq1=0.0
        xd=xd+facq1
 50    CONTINUE
      ENDIF
      IF(i.ne.nmax)THEN
       ip=i+1
       DO 100 k=ip,nmax
        xc=xc+(p(k,l)*rrc(k,i))
        facq2=((1.0-(p(k,l)/deg(k)))*rrc(i,k))
        IF(facq2.lt.0.0)facq2=0.0
        xd=xd+facq2
 100   CONTINUE
      ENDIF
c ad1 allan: if tunnel ionisation is included and laser is on - set
c        coll excitation/deexcitation rates to zero.
c     if(nltnl.and. (xplas1(l+1).gt.1.e00))xc=0.
c     if(nltnl.and.(xplas1(l+1).gt.1.e00))xd=0.

      acf=xc
      bcf=xd
c
      acf2=0.
      bcf2=0.
      srth=sqrt(theta)
      feln=(zst(l)*rho)/a(l)
      un=eng(i,l)/theta
      eun=0.0
      IF(un.le.170.0)eun=exp(-un)
c calcualte collisional ionisation and three-body recombination
c contributions
c calculate the function gmnfac
      zn=real(i)
      y=log10(0.25*un)
      y2=y*y
      y3=y2*y
      y4=y3*y
      fy=0.23+(0.046*y)+(0.1074*y2)-(0.0459*y3)-(0.01505*y4)
      za=12.18*exp(-(zn/(zn+5.0)))*(1.0+(0.0335*zn))
      zb=1.0-(0.622/qs(i,l))-(0.0745/(qs(i,l)*qs(i,l)))
      zc=1.0
      IF(un.le.170.0)zc=zc-exp(-un)
      gmnfac=za*zb*zc*abs(fy)
      IF(gmnfac.le.0.0)gmnfac=0.0
      wa=2.08E+13*feln*srth*eun*gmnfac
      wb=1.0/(eng(i,l)**2)
      xa=wa*wb
c[
      waa=2.08E+13*feln*srth*gmnfac
      cf=(deg(i)*rho)/(317.0*a(l)*theta*srth)
      waa=waa*cf
      xb=waa*wb
c
c ad1 allan: if tunnel ionisation is included and laser is on - set
c        coll ionisatio/recombination rates to zero.
c     if(nltnl.and.(xplas1(l+1).gt.1.0e00))xa=0.
c     if(nltnl.and.(xplas1(l+1).gt.1.0e00))xb=0.

      cbf(i)=xa
      cfb(i)=xb
c[c[
      acf2=acf2+(zst(l)*xb)
      bcf2=bcf2+xa

c add dielectronic recombination if requested
c do not include dielectronic recombination if radiation field
c is set
      IF(idrflg.ne.0)acf2=acf2+zst(l)*rdr(i)
c
      atie(i)=0.
      atir(i)=0.
      IF(nltnl)THEN
c        calculate tunnel ionization rate and ati energy from level i
c        tunneling ionization rate according to ammosov et al from paper
c        by penetrante
c        tunnel ionization rate(per s) and ati energy in eV

c        find position of ionisation pot in array ui
         mz=int(z(l)+0.5)
         mzm1=mz-1
         ibot=(mz*mzm1)/2
c          zleft = sum pops in and inside level i
           zleft=0.
           do 61 i1=1,i
            zleft=zleft+pu(i1,l)
61         continue
           if(zleft.gt.1.e-8) then
            igrd=0
            IF(zleft.le.2.)THEN
             igrd=1
            ELSEIF(zleft.le.10.)THEN
             igrd=2
            ELSEIF(zleft.le.28.)Then
             igrd=3
            ELSEIF(zleft.le.60.)Then
             igrd=4
            ELSEIF(zleft.le.103.)Then
             igrd=5
            ENDIF
c           ionization of electron nele
            nele=int(z(l)-zleft) + 1
            if(nele.lt.1)nele=1
            if(nele.gt.mz)nele=mz
      CALL tunnel(z(l),ibot,nele,igrd,i,xplas1(l+1),yamda,
     & atir1,atie1)
            atir(i)=atir1
            atie(i)=atie1
c ad1 allan: option to set ati energy to 0
c           atie(i)=0.
           ENDIF
           bcf2=bcf2+atir(i)
      ENDIF
c radiative ionisation/ recombination rates
      wa=7.86864E+09
      wb=(zs(i,l)**4)/(fn(i)**5)
      za=0.
      zb=wa*wb*cf*eintm(un)
      zbs=0.
      IF(iradf1.ne.0)THEN
       CALL radint(un,theta,radti,radtr)
       IF(radtr.gt.0.0.and.radti.gt.0.0)THEN
        zbs=(wa*wb*cf*radtr)-zb
        za=wa*wb*radti
       ELSE
        zbs=0.
        za=0.0
       ENDIF
      ELSEIF(iradf2.eq.0)THEN
       za=0.0
      ELSE
       radtsr=0.
       radti=0.
       IF(xhnu.ge.eng(i,l))THEN
        radti=dihnu/xhnu
        unn=(xhnu-eng(i,l))/theta
        eun=0.0
        IF(unn.le.170.0)THEN
         eun=exp(-unn)
         radtsr=dihnu*eun/xhnu
        ENDIF
       ENDIF
       zbs=wa*wb*cf*radtsr
       za=wa*wb*radti
      ENDIF
c
c ad1 allan: if tunnel ionisation is included and laser is on - set
c        radiative rates to zero
c     if(nltnl.and.(xplas1(l+1).gt.1.0e00))za=0.
c     if(nltnl.and.(xplas1(l+1).gt.1.0e00))zb=0.
c     if(nltnl.and.(xplas1(l+1).gt.1.0e00))zbs=0.

c       spontaneous radiative recombination
      rfb(i)=zb*(1.-fbtrap)
c       stimulated
      rfbst(i)=zbs
      IF(rfbst(i).lt.0.)rfbst(i)=0.
      acf2=acf2+zst(l)*(rfb(i)+rfbst(i))
c       bf ionization for use in rate equation (per electron
c       per atom
      rbf(i)=za
      bcf2=bcf2+za

c     x ray bf ionization for une in absorb (should total and
c     per metre
      xrbf(l)=0.
      IF(abs(rbf(i)).gt.underf)THEN
       an=6.022045E23
       dnion=rho*an/a(l)
c       given laser intensity in (w/m**2), xhnu in kev calculate bf
c       absorption in units of   (per m) for use in absorb
c       first per cm
       wa=1.99E-23
       rncru=wa*wb/(xhnu**3)
       unn=(xhnu-eng(i,l))/theta
       eun=0.0
       IF(unn.le.170.0)eun=exp(-unn)
       rcnrd=wa*wb*cf*eun/(xhnu**3)
       IF(xhnu.ge.eng(i,l))THEN
        qi=1.-p(i,l)/deg(i)
c       multiply by 100. to convert from /cm to /m
        xrbf(l)=xrbf(l)+dnion*max(0.,(p(i,l)*rncru-zst(l)*rcnrd*qi))
       ENDIF
      ENDIF
      RETURN
      END


      SUBROUTINE tunnel(zz,ibot,nele,igrd,ilev,ai,amda,rate,eati)
      PARAMETER(underf=1.E-30)
c   given element zz
c   ground state ionization pot (eV) is ui(ibot+1)
c   laser intensity ai(w/m**2), wavelength amda (m)
c   calculate ionisation rate (per s) and ati energy in eV
c   of nele'th electron from excited level ilev (ground igrd)
c   in this version zz= 36 maximum
c   ionization potential and ATI are now in eV
c
      DIMENSION ui(666)
      DATA(ui(i),i=1,105)/13.598,24.587,54.416,5.392,75.638,122.451,
     &     9.322,18.211,153.893,217.713,8.298,25.154,37.930,259.368,
     &     340.217,11.260,24.383,47.887,64.492,392.077,489.983,14.534,
     &     29.601,47.448,77.742,97.888,552.057,667.029,13.618,35.116,
     &     54.934,77.412,113.896,138.116,739.315,871.39,17.422,34.970,
     &     62.707,87.138,114.24,157.161,185.182,953.886,1103.09,21.564,
     &     40.962,63.45,97.11,126.21,157.93,207.27,239.09,1195.797,
     &     1362.16,5.139,47.286,71.64,98.91,138.39,172.15,208.47,264.18,
     &     299.87,1465.091,1648.66,7.646,15.035,80.143,109.24,141.26,
     &     186.50,224.94,265.90,327.95,367.53,1761.802,1962.61,5.986,
     &     18.828,28.447,119.99,153.71,190.47,241.43,284.59,330.21,
     &     398.57,442.07,2085.983,2304.08,8.151,16.345,33.492,45.141,
     &     166.77,205.05,246.52,303.17,351.10,410.43,476.06,523.5,
     &     2437.676,2673.11/
      DATA(ui(i),i=106,210)/10.486,19.725,30.18,51.37,65.023,220.43,
     &     263.22,309.41,371.73,424.50,479.57,560.41,611.85,2816.943,
     &     3069.76,10.36,23.33,34.83,47.304,72.74,88.049,280.93,328.23,
     &     379.1,447.09,504.78,546.65,651.63,707.14,3223.836,3494.1,
     &     12.967,23.81,39.61,53.46,67.8,96.98,114.193,348.28,400.05,
     &     455.62,529.26,591.97,656.69,769.74,809.39,3658.425,3946.19,
     &     15.759,27.624,40.79,59.81,75.02,91.007,124.319,143.406,
     &     422.44,478.68,538.95,618.24,686.09,755.73,854.75,918.00,
     &     4120.788,4426.11,4.341,31.625,45.72,60.91,82.66,99.89,117.56,
     &     154.86,175.814,503.44,564.13,629.09,714.02,787.13,861.77,
     &     966.0,1034.0,4610.995,4933.93,6.113,11.871,50.908,67.10,
     &     84.41,108.78,127.7,147.24,188.54,211.270,591.25,656.39,
     &     726.03,816.61,895.12,974.0,1085.0,1157.0,5129.045,5469.74/
      DATA(ui(i),i=211,231)/6.54,12.8,24.76,73.47,91.66,111.1,137.,
     &     158.7,180.02,225.32,249.832,685.89,755.47,829.79,926.,1009.,
     &     1094.,1210.,1288.,5675.,6033.6/
      DATA(ui(i),i=232,253)/6.84,13.58,27.491,43.266,99.2,119.36,140.8,
     &     168.5,193.,215.91,265.31,291.497,787.33,861.33,940.36,1042.,
     &     1131.,1220.,1342.,1425.,6249.,6625.6/
      DATA(ui(i),i=254,276)/6.74,14.65,29.31,46.707,56.23,128.12,150.17,
     &     173.7,204.,230.,255.04,308.25,336.267,895.58,974.02,1058.,
     &     1165.,1259.,1353.,1482.,1569.,6851.,7245.9/
      DATA(ui(i),i=277,300)/6.766,16.5,30.96,49.1,71.,90.56,161.1,184.2,
     &     209.8,242.,270.,297.,354.,384.3,1010.64,1093.,1182.,1295.,
     &     1394.,1493.,1627.,1720.,7482.,7984.5/
      DATA(ui(i),i=301,325)/7.435,15.64,33.667,53.,73.,97.,119.27,196.2,
     &     221.4,248.5,283.,313.,342.,404.,435.3,1136.2,1220.,1313.,
     &     1431.,1563.,1640.,1780.,1879.,8141.,8571.5/
c Fe
      DATA(ui(i),i=326,351)/7.87,16.18,30.651,56.,77.,101.,126.,151.06,
     &     235.04,262.1,290.,328.,360.,391.,456.,489.5,1266.1,1353.,
     &     1451.,1575.,1685.,1794.,1940.,2045.,8828.,9277.2/
c Co
      DATA(ui(i),i=352,378)/7.86,17.06,33.5,57.52,81.,105.,131.,159.,
     & 186.13,276.6,305.,335.,375.,409.,441.,512.,546.8,1403.,1493.,
     & 1595.,1724.,1841.,1955.,2106.,2218.,9544.,10011./
c Ni
      DATA(ui(i),i=379,406)/7.635,18.168,35.17,59.,82.8,110.,136.,
     & 165.,195.,224.5,321.2,352.5,384.,426.,461.,495.,570.,607.2,
     & 1546.9,1639.,1747.,1881.,2003.,2123.,2279.,2398.,10288.,10775./
c Cu
      DATA(ui(i),i=407,435)/7.726,20.292,36.83,55.2,85.,111.,141.,
     & 170.,201.,234.,266.,370.,402.,435.,479.,517.,553.,631.,670.9,
     & 1697.9,1793.,1905.,2045.,2173.,2298.,2459.,2585.,11062.,11566./
c Zn
      DATA(ui(i),i=436,465)/9.394,17.964,39.722,63.,87.,114.,144.,176.,
     & 207.,241.,276.,310.,421.,454.,489.,536.,575.,613.,696.,737.8,
     & 1855.9,1953.,2070.,2215.,2349.,2479.,2646.,2780.,11864.,12387./
c Ga
      DATA(ui(i),i=466,496)/5.999,20.51,30.7,64.,90.,117.,147.,179.,
     & 214.,248.,284.,321.,358.,457.,510.,546.,596.,637.,676.,764.,
     & 808.,2020.,2120.,2242.,2392.,2533.,2668.,2840.,2982.,
     & 12695.,13266./
c Ge
      DATA(ui(i),i=497,528)/7.899,15.934,34.22,45.71,93.5,120.,151.,
     & 183.,217.,255.,291.,330.,369.,409.,533.,568.,607.,658.,701.,
     & 743.,834.,881.,2192.,2294.,2420.,2576.,2723.,2863.,3041.,3192.,
     & 13556.,14117./
c As
      DATA(ui(i),i=529,561)/9.81,18.633,28.351,50.13,62.63,127.6,154.,
     & 187.,222.,259.,300.,338.,379.,421.,463.,594.,630.,670.,724.,769.,
     & 812.,908.,953.,2372.,2474.,2606.,2767.,2920.,3065.,3248.,3409.,
     & 14444.,15026./
c Se
      DATA(ui(i),i=562,595)/9.752,21.19,30.820,42.944,68.4,81.7,155.4,
     & 191.,227.,264.,304.,347.,388.,431.,475.,520.,657.,695.,736.,
     & 793.,839.,884.,985.,1037.,2558.,2661.,2798.,2964.,3123.,3274.,
     & 3463.,3633.,15365.,15964./
c Br
      DATA(ui(i),i=596,630)/11.814,21.8,35.5,48.,60.,89.,101.,
     & 192.8,232.,270.,310.,352.,398.,441.,485.,533.,580.,724.,
     & 762.,806.,864.,913.,960.,1064.,1120.,2752.,2855.,2997.,
     & 3169.,3334.,3490.,3664.,3865.,16313.,16933./
c Kr
      DATA(ui(i),i=631,666)/13.999,24.359,36.95,53.,65.,78.,111.,
     & 123.,230.9,275.,316.,358.,403.,451.,497.,545.,593.,643.,
     & 794.,833.,878.,939.,989.,1038.,1147.,1206.,2953.,3056.,3203.,
     & 3380.,3551.,3712.,3912.,4105.,17292.,17932./
c
c
      rate=0.
      eati=0.

      IF(zz.gt.36.0) Then
       WRITE(6,*)' *** TUNNEL: element ',zz,' not implemented'
       rate=0.
       eati=0.
       stop 36
      ENDIF
      if(ilev.lt.igrd.or.igrd.lt.1) return
      uig=ui(ibot+nele)
      uion=uig*igrd**2/ilev**2
c ionization of nele'th  electron
      zn=real(nele)
c
      IF(ai.le.1.e10)THEN
       rate=0.
       eati=0.
      ELSE
c       atomic field in v/m
       eau=5.142E+11
       oau=4.134E+16
       e0=sqrt(753.46*ai)
c effective n
       an=sqrt(zn*zn*13.598/uion)
       equo=eau/e0
       bb=2*an-1.5
       t1=(10.873*(zn**3)*equo/(an**4))**bb
       cc=(2*zn**3*equo/3./an**3)
       t2=0.
       IF(cc.lt.175.)t2=exp(-cc)
       rate=1.611*oau*zn*zn*t1*t2/an**4.5
c       assuming plane polarized laser calculate ati energy
       pi=acos(-1.)
c       equiv in ev
       equiv=9.3372E-6*ai*amda**2
       t1=0.
       t2=0.
       dfi=pi/2./100.
c       integrate
       fi=-0.5*dfi
       DO 50 in=1,100
        fi=fi+dfi
        xx=sin(fi)
        yy=1.-xx*xx
        tt=0.
        cox=cc/xx
        IF(cox.lt.175.)tt=xx**(-bb)*exp(-cc/xx)*dfi
        t1=t1+tt*yy
        t2=t2+tt
 50    CONTINUE
       rat2=0.
       IF(abs(t2).gt.underf)rat2=2.*t1/t2
       eati=equiv*rat2
      ENDIF
      RETURN
      END


      SUBROUTINE interp(xdat,ydat,x,y,n)
c parameters      : xdat   : x points with a known y value
c                 : ydat   : corresponding y values
c                 : x      : point at which interpolation is required
c                 : y      : interpolated value
c                 : n      : number of points at which values are known
      REAL xdat(n),ydat(n),x,y
      IF(x.lt.xdat(1))THEN
       i=1
       j=2
      ELSEIF(x.gt.xdat(n))THEN
       i=n-1
       j=n
      ELSE
       DO 50 k=1,n-1
        IF((x.ge.xdat(k)).and.(x.le.xdat(k+1)))THEN
         i=k
         j=k+1
        ENDIF
 50    CONTINUE
      ENDIF
      y=ydat(i)+(ydat(j)-ydat(i))*(x-xdat(i))/(xdat(j)-xdat(i))
      RETURN
      END


      SUBROUTINE radint(un,theta,radti,radtr)
c
      COMMON /rad   / radf(1501),hrad,nrad
      DIMENSION ty(1501),tz(1501)
c
c check that the ionisation energy lies within the specified radiation
c field
c
      radhi=(real(nrad)-10.0)*hrad
      radti=0.0
      radtr=0.0
      eionz=theta*un
      IF(eionz.le.radhi)THEN
c
c find the grid point corresponding to the ionisation energy
c
       rnie=((theta*un)/hrad)+1.0
       nieint=int(rnie)
       deltan=rnie-real(nieint)
       nieint=nieint+1
       IF(nieint.eq.1)nieint=2
c
c calculate the initial and final points of the integration for
c simpson's rule
c
       nst=nieint
       nfin=nrad
       ndiff=nrad-nieint
       IF(mod(ndiff,2).ne.0)nfin=nfin-1
c
c fill the integrand array
c
       ncount=nst-1
 50    ekev=real(ncount)*hrad
       ncount=ncount+1
       wa=1.0/ekev
       unn=ekev/theta
       wb=0.0
       IF(unn-un.le.170.0)wb=exp(un-unn)
       wc=radf(ncount)
       ty(ncount)=wa*wc
       tz(ncount)=wa*wb*(1.0+wc)
       IF(ncount.lt.nfin)GOTO 50
c
c perform integration from ionisation energy to infinity
c
       radti=sint1(nst,nfin,hrad,ty,deltan)
       radtr=sint1(nst,nfin,hrad,tz,deltan)
      ENDIF
      RETURN
      END
c
      FUNCTION sint1(nst,nfin,h,ty,deltan)
c
c this function performs a euler-maclaurin summation  between two
c specified points
c nst is the first point
c nfin is the last point
c h is the step size
c ty is the value of the integrand at grid point i
c
      DIMENSION ty(1501)
      nfinm=nfin-2
      wa=0.0
      hdt=h/3.0
      DO 100 i=nst,nfinm,2
       wa=wa+h*(0.5*ty(i)+ty(i+1)+0.5*ty(i+2))
 100  CONTINUE
      wb=ty(nst)*(1.0-deltan)*h
      grad=(ty(nst+1)-ty(nst))/h
      wc=-0.5*grad*(((1.0-deltan)*h)**2)
      sint1=wa+wb+wc
      RETURN
      END


      SUBROUTINE radb(un,theta,radtb)
c
      COMMON /rad   / radf(1501),hrad,nrad
c
c check that the excitation energy lies within the specified radiation
c field
      radhi=(real(nrad)-10.0)*hrad
      radtb=0.0
      eez=theta*un
      IF(eez.le.radhi)THEN
c
c find the grid points corresponding to the excitation energy
c
       rnee=(eez/hrad)+1.0
       neeint=int(rnee)
       deltan=rnee-real(neeint)
c
c interpolate to find the radiation field
c
       grad=radf(neeint+1)-radf(neeint)
       radtb=radf(neeint)+(grad*deltan)
      ENDIF
      RETURN
      END

c
      SUBROUTINE setr(thetar)
c
      PARAMETER(kk=1001)
      COMMON /rad   / radf(1501),hrad,nrad
c

      COMMON /radflg/ yamda,xhnu,xihnu,dihnu,xplas1(kk),xrbf(kk),
     &                atrbf(kk),corhf1(kk),corhf(kk),iradf1,iradf2

c specify the modal photon density of the radiation field
c set up for a planckian at the specified radiation temperature
c
      nrad=1501
      nradm=nrad-1
      hrad=max(30.0*thetar,1.1*xhnu)/real(nradm)
      ekev=hrad
      radf(1)=0.0
      DO 100 i=1,nradm
       unr=ekev/thetar
       IF(unr.lt.170.)THEN
        radf(i+1)=1.0/((exp(unr))-1.0)
       ELSE
        radf(i+1)=0.
       ENDIF
       ekev=ekev+hrad
 100  CONTINUE
      RETURN
      END


      SUBROUTINE setzi(thetar)
      PARAMETER(kk=1001)
      COMMON /rad   / radf(1501),hrad,nrad
c

      COMMON /radflg/ yamda,xhnu,xihnu,dihnu,xplas1(kk),xrbf(kk),
     &                atrbf(kk),corhf1(kk),corhf(kk),iradf1,iradf2

c specify the modal photon density of the x-ray radiation field
c
      dihnu=0.
      IF(iradf1.eq.0)THEN
       nrad=1501
       nradm=nrad-1
       hrad=max(30.0*thetar,1.1*xhnu)/real(nradm)
       radf(1)=0.0
       DO 50 i=1,nradm
        radf(i+1)=0.
 50    CONTINUE
      ENDIF
c check that the x-ray photon energy lies within the energy range
c
      radhi=(real(nrad)-10.0)*hrad
      IF(xhnu.gt.radhi)THEN
       WRITE(6,*)' xhnu is .gt. radhi program stoped'
       STOP
      ENDIF
c
c find the grid points corresponding to the x-ray photon energy
c
      rnee=(xhnu/hrad)+1.0
      neeint=int(rnee)
c
c add x-ray photon modal density  at appropriate grid point
c
      dihnu=3.17868E-32*xihnu/(xhnu**3)
      radf(neeint)=radf(neeint)+dihnu/hrad
c     write(6,*)yamda,xhnu,xihnu,hrad,neeint,radf(neeint-1),radf(neeint
c    & radf(neeint+1)
      RETURN
      END


      SUBROUTINE linint(xdat,ydat,x,y,n)
      INTEGER n,i,j,k
      REAL xdat(n),ydat(n),x,y
      IF(x.lt.xdat(1))THEN
       i=1
       j=2
      ELSEIF(x.gt.xdat(n))THEN
       i=n-1
       j=n
      ELSE
       DO 50 k=1,n-1
        IF((x.ge.xdat(k)).and.(x.le.xdat(k+1)))THEN
         i=k
         j=k+1
        ENDIF
 50    CONTINUE
      ENDIF
      y=ydat(i)+((ydat(j)-ydat(i))*((x-xdat(i))/(xdat(j)-xdat(i))))
      END

