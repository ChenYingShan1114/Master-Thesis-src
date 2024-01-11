      program medgraf
      
      parameter(nshot=100,Avdro=6.022e23)
      real edge(1000), vel(1000), centre(1000)
     :, density(1000), press(1000), Ti(1000)
     :, Te(1000), Zave(1000), tim(nshot),piq(1000)
     :, xni(1000), xne(1000),xaser1(1000),nprnt
     
      character*10 char(6)*80,crun*80,ct*40,ct2*40,color(6)
      data color/'black','blue','green','red','cyan','magenta'/
      logical nlcri1,nlbrms,nlpfe,nlemp,nlabs,nlburn,nldepo
     :,nlecon,nlicon,nlfuse,nlpfi,nlx,nlmove,nlhf,nltnl,nlomt3(50)
      integer nmesh

      namelist /newrun/ 
     :  xamda1,   gauss,     anpuls,     toff,
     :  plenth,  pmax,    pmult,
     :  xaser1,
     :  ngeom,         piq,    teini,     tiini,
     :  mesh,         rini,   rhogas,   rf1,
     :  xz,          xmass,     fne,
     :  zglas,       drglas,  roglas,   rf2,
     :  xz2,         xmass2,    fne2,
     :  zplas,        drplas,  roplas,      rf3,
     :  xz1,          xmass1,     hydrog,
     :  nprnt,	tstop,	nrun,
     :  np3,	nlemp,  pondf, nlcri1,
     :  nlbrms,	flimit,	saha,
     :  anabs,	fhot,	fthot,	rhot,
     :  state,	nlpfe,  dtmax,
     :  ak4, nlabs, nlburn, nldepo, nlecon, nlicon, nlfuse,
     :  nlpfi, nlx, nlmove, nlhf, nltnl,
     :  ilosta, ihista, ifrsta, nt1, nt2, nlomt3
 

c  fetch namelist
      open (10,file='med.ind')
      read (10,newrun)

      open (20,file='med_graph.asc')
      open (25,file='graph.id')

c  skip header
      do i=1,6
        read (20,'(a)') char(i)
      end do

c  # mesh points
      read (20,*) nmesh

      i = 0
  1   continue
      i = i+1
      read (20,*,end=2) tim(i),tim2

      if (i.ge.2 .and. tim(i).eq.tim(i-1)) goto 2

      do j=1,nmesh
        read (20,*) edge(j), vel(j), centre(j), density(j)
     :               ,press(j), Te(j), Ti(j), Zave(j)
      end do
 101  format (1p,3e11.4,2e10.3,3e9.2)

      read (20,*) edge(nmesh+1), vel(nmesh+1)

c  compute electron and ion number densities

c  1st layer
      do j=1,nmesh-zglas
         xni(j) = density(j)/xmass*Avdro
      end do

c  2nd layer
      do j=nmesh-zglas+1,nmesh
         xni(j) = density(j)/xmass2*Avdro
      end do

      do j=1,nmesh
         xne(j) = xni(j)*Zave(j)
      end do

      
c  compute correct time relative to pulse max
      timtomax = tim(i)*1.e12 - pmult*1.e12*plenth
      xpr = abs(nprnt)/10.
      tim1 = xpr*nint(tim(i)*1.e12/xpr)
      tps = xpr*nint(timtomax/xpr)
c      write(*,*) tim(i)*1.e12,timtomax,npr,tim1,tps

      call chr(tps,3,ct2,lct2)
      call chr(tim1,3,ct,lct)
      write (25,'(2(a,2x),a)') ct(1:lct),ct2(1:lct2),color(mod(i-1,6)+1)
     
      open (30,file='rho'//ct(1:lct))
      write (30,'(2e12.4)') (centre(j)*1.e4,density(j),j=1,nmesh)
      close (30)

      open (30,file='ni'//ct(1:lct))
      write (30,'(2e12.4)') (centre(j)*1.e4,xni(j),j=1,nmesh)
      close (30)

      open (30,file='ne'//ct(1:lct))
      write (30,'(2e12.4)') (centre(j)*1.e4,xne(j),j=1,nmesh)
      close (30)

      open (30,file='ti'//ct(1:lct))
      write (30,'(2e12.4)') (centre(j)*1.e4,Ti(j),j=1,nmesh)
      close (30)

      open (30,file='te'//ct(1:lct))
      write (30,'(2e12.4)') (centre(j)*1.e4,Te(j),j=1,nmesh)
      close (30)

      open (30,file='p'//ct(1:lct))
      write (30,'(2e12.4)') (centre(j)*1.e4,press(j),j=1,nmesh)
      close (30)

      open (30,file='vel'//ct(1:lct))
      write (30,'(2e13.4)') (centre(j)*1.e4,vel(j),j=1,nmesh)
      close (30)

      open (30,file='z'//ct(1:lct))
      write (30,'(2e12.4)') (centre(j)*1.e4,Zave(j),j=1,nmesh)
      close (30)



      goto 1
   2  continue
      nstep = i-1
      end

      
