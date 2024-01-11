      program timehist
      
      parameter(nshot=200,Avdro=6.022e23,nposm=5)
      real edge(300), vel(300), centre(300)
     :, density(300), press(300), Ti(300)
     :, Te(300), Zave(300),piq(300)
     :, xni(300), xne(300),xaser1(300),nprnt,xpos(nposm)
     :,  tim(nshot),Tehist(nshot), nehist(nshot), zhist(nshot)
     
      character*10 char(6)*80,crun*80,ct*40,ct2*40,color(6)
      data color/'black','blue','green','red','cyan','magenta'/
      logical nlcri1,nlbrms,nlpfe,nlemp,nlabs,nlburn,nldepo
     :,nlecon,nlicon,nlfuse,nlpfi,nlx,nlmove,nlhf,nltnl,found
      integer nmesh

c  Medusa namelist: variable names must match those in med.ind (input deck)

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
     :  ilosta, ihista, ifrsta, nt1, nt2
 

c  fetch namelist
      open (10,file='med.ind')
      read (10,newrun)

      open (20,file='med_graph.asc')
      open (25,file='graph.id')

c  get probe positions for time histories
      write (6,*) 'Give number of probe positions:'
      read (5,*) npos
      write (6,*) 'Give probe positions (mu-m):'
      
      do i=1,npos
	read(5,*) xpos(i)
c convert to cm
	xpos(i)=xpos(i)*1.e-4
      end do

c  skip header
      do i=1,6
        read (20,'(a)') char(i)
      end do

      open (30,file='hist.dat')

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
      read (20,*) edge(nmesh+1), vel(nmesh+1)

      
c  compute correct time relative to pulse max
      timtomax = tim(i)*1.e12 - pmult*1.e12*plenth
      xpr = abs(nprnt)/10.
      tim1 = xpr*nint(tim(i)*1.e12/xpr)
      tps = xpr*nint(timtomax/xpr)
      tns = tps/1000.
      write(*,*) timtomax

c      call chr(tps,3,ct2,lct2)
c      call chr(tim1,3,ct,lct)
c      write (25,'(2(a,2x),a)') ct(1:lct),ct2(1:lct2),color(mod(i-1,6)+1)
     
c 	find temp, density at probe positions

       do iprobe=1,npos
c  start indices
	icm = 1
	icp = 2
	found = .false.

c  sweep down
	is = icm
	do while (.not.found .and. is.le.nmesh)
	   if (centre(is).gt.xpos(iprobe)) then
             found = .true.
             frac=(xpos(iprobe)-centre(is-1))/(centre(is)-centre(is-1))
             tehist(iprobe) = te(is-1)+frac*(te(is)-te(is-1))
             nehist(iprobe) = xne(is-1)+frac*(xne(is)-xne(is-1))
             zhist(iprobe) = zave(is-1)+frac*(zave(is)-zave(is-1))
	   endif	  
   	   is = is+1
	end do
	if (.not.found.or.is.gt.nmesh) then
	  tehist(iprobe) = te(nmesh)
	  nehist(iprobe) = xne(nmesh)
	  zhist(iprobe) = zave(nmesh)
	endif
       end do

      call chr(1+npos*3.,0,ct2,lct2)
      ct = '('//ct2(1:lct2)//'(1pe12.4))'
      lct = lench(ct)
       write (30,ct(1:lct)) 
c     : (timtomax,tehist(iprobe),nehist(iprobe),iprobe=1,npos)
     : tns,(tehist(iprobe),nehist(iprobe),zhist(iprobe),iprobe=1,npos)

      goto 1
   2  continue

      close (30)


      nstep = i-1
      end

      
