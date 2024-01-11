      program drelp
      
c     ***************************************************************
c     
c     Dynamic RElativistic Laser Propagation
c     
c     Version 3.0
c     
c     ***************************************************************
c     
c     Paul Gibbon
c     
c     Saclay,  June 1993
c     Revision and public release, Juelich, March 2004
c     
c     ***************************************************************
c     
c     Features:
c     
c     -  2D cylindrical em envelope propagation
c     -           -> relativistic self-focussing
c     -  adiabatic electrons -> ponderomotive channel formation
c     -  mobile ions treated as lagrangian fluid slices
c     -  field ionisation using Keldysh model for H and He
c     
c     ***************************************************************
c     
c     Fortran file units:
c     10 - namelist inputs
c     15 - printed o/p
c     50 - 1D graphics header  .oddata
c     70 - 2D snapshots        .snap
c
c  **  See README for detailed description of input parameters **
c  **  and ouput files 
c
      include 'dynafoc.h'

      call stamp(15,1)
      call stamp(6,1)
      call init
      call setup
      call arinit
      call grid

c     initialise pulse
      if (ilas.ge.4) call pulinit

c     initialise density
      call ioniz(0.,t0)
      call density(t0)

c     output initial parameters
      call outini

c     absorber at RHB
      call absorb

c     call grinp
      call diagnostics(t0)

c     timing stuff
      call cputime(time,iv)
      tim1=time
      tdiag=0.
      tem=0.
      tspline=0.
      tplas=0.
      tioniz=0.
 
c     Outer loop over em timesteps

      do itime=1,nt
         if (ilas.le.3) then
c     antenna at LHB - fixed pulse shape
            call launch
         else 
c     NLSE solution
            call emprop(dt,tem)
         endif
         

c     ionization/gas density profile
         call ioniz(dts,tioniz)


c     static/adiabatic electrons 
         
         call density(tplas)

c     rezoning: push plasma back through EM grid with dz=v_g.dt            
         if (nz.gt.1) then 
            call rezone(dt)
         else
c  shift ion profile back even if steady state
            zedge=zedge-dt
            zend = zend-dt
         endif

         call diagnostics(tdiag)


c     update timestep - not needed with implicit rate equations
         if (mioniz.eq.1) then
            dts=min(dt,dtioniz)
c     allow up to 5 subcycles for ionization
            ids = dt/dts
            icycem = min(ids,isubm)
c            icycem=1
            dts = dt/icycem
         endif
      end do

c     Timing

      call cputime(time,iflag)
      tim2=time
      write (15,111)
      write (6,111)
 111  format(//'Loop timing (in seconds):'
     :     /'------------------------')
      ttot=tim2-tim1
      tloop=ttot/nt
      call r0form('total ',ttot,'f12.2')
      call r0form('    em',tem,'f12.2')
      call r0form('plasma',tplas,'f12.2')
      call r0form('ioniz ',tioniz,'f12.2')
      call r0form('  diag',tdiag,'f12.2')
      call r0form('loop  ',tloop,'f13.3')
      call blank
      
      call ghist
      call powcon
      call stamp(6,2)
      call stamp(15,2)

      close(11)
      end
