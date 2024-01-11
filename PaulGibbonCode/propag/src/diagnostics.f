
      
c     =============================================================
c     
c     Diagnostic output
c     
c     =============================================================


      subroutine diagnostics(tim)
      include 'dynafoc.h'
      character ct*40,cbl*16

      call cputime(tim1,v)
      

c   offset for focussing geometry
      if (mfoc.eq.1) then
         zoff = Rfoc
      else
         zoff = 0.
      endif


c   convert cpu time to char variable
      cbl='                '
      tlb = (timec-zoff)*rfac

      call chr(tlb,1,ct,lct)
      lctime=12-lct
      ctime='z='//ct(1:lct)//'mu'


c     graphical snapshots every igr timesteps
      if (mod(itime,igr).eq.0) then
         call grafout
      endif


c     surface movie plots: 2D snapshots
      if (mod(itime,imovie).eq.0) then
         call gmovie
      endif


c     r,z plots - time history of radial slice
      if (mod(itime,izsurf).eq.0) then
         call grzsurf
      endif


c     written o/p every iout timesteps
      if (mod(itime,iout).eq.0) then
         if (itime.eq.0)
     :        write (6,'(4a15)') ' # dt','ct/mu','edge/mu','% runtime'
         write (6,'(i15,3f15.2)') itime, tlb, zedge*rfac, timec/trun*100
         write (15,*) 'plasma front edge (mu)',zedge*rfac
     :,' end ',zend*rfac
      endif


c     store time-histories every ihist
      if (mod(itime,ihist).eq.0) then       
         call flhist
         call emhist
      endif

c     em spectral probe - sample every icycem timesteps 
c     start in far-field (2*zfoc)

      if (timec.ge.tprobe .and. nzp.gt.1) then
         call spec(dt)
         if (mod(itime,iout).eq.0) 
     :        write (15,*) 'Probe at ',zdens(kprobe)
      endif
      call cputime(tim2,v)
      tim = tim + tim2-tim1
      timec = timec+dt


      end
