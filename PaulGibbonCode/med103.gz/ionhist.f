      program timehist
      
      parameter(nshot=200,Avdro=6.022e23,nposm=5)
      real edge(300), vel(300), centre(300)
     :, density(300), press(300), Ti(300)
     :, Te(300), Zave(300),piq(300)
     :, xni(300), xne(300),xaser1(300),nprnt,xpos(nposm)
     :,  tim(nshot),Tehist(nshot), nehist(nshot)
     :,  hion1(nshot),hion2(nshot)
     :, hion(10,300)
      character*10 char(6)*80,crun*80,ct*40,ct2*40,color(6)
      data color/'black','blue','green','red','cyan','magenta'/
      logical nlcri1,nlbrms,nlpfe,nlemp,nlabs,nlburn,nldepo
     :,nlecon,nlicon,nlfuse,nlpfi,nlx,nlmove,nlhf,nltnl,found
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
     :  ilosta, ihista, ifrsta, nt1, nt2
 

c  fetch namelist
      open (10,file='med.ind')
      read (10,newrun)

      open (20,file='med_graph.asc')
      open (25,file='med.expops')

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
c  skip header
      do i=1,4
        read (25,'(a)') char(i)
      end do

      open (30,file='ionhist.dat')

c  # mesh points
      read (20,*) nmesh

      i = 0
  1   continue
      i = i+1
      read (25,*,end=2) tim(i),tim2,power

      if (i.ge.2 .and. tim(i).eq.tim(i-1)) goto 2

      do j=1,nmesh
        read (25,*) edge(j), vel(j), density(j), Zave(j)
     : , Te(j), Ti(j), (hion(iion,j),iion=1,10)
      end do

      do j=1,nmesh
c        write (6,'(i4,4(1pe12.3))') j,edge(j), density(j),
c     :(hion(iion,j),iion=1,2)
      end do

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
      read (25,*) edge(nmesh+1), vel(nmesh+1)

      
c  compute correct time relative to pulse max
      timtomax = tim(i)*1.e12 - pmult*1.e12*plenth
      xpr = abs(nprnt)/10.
      tim1 = xpr*nint(tim(i)*1.e12/xpr)
      tps = xpr*nint(timtomax/xpr)
      tns = tps/1000.
      write (6,'(f12.4)') tns

c 	find temp, density at probe positions

       do iprobe=1,npos
c  start indices
	icm = 1
	icp = 2
	found = .false.

c  sweep down
	is = icm
	centre(is)=(edge(is)+edge(is+1))/2.
	do while (.not.found .and. is.le.nmesh)
	   centre(is)=(edge(is)+edge(is+1))/2.
	   if (centre(is).gt.xpos(iprobe)) then
             found = .true.
             frac=(xpos(iprobe)-centre(is-1))/(centre(is)-centre(is-1))
             hion1(iprobe) = hion(1,is-1)+frac*(hion(1,is)-hion(1,is-1))
             hion2(iprobe) = hion(2,is-1)+frac*(hion(2,is)-hion(2,is-1))
	   endif	  
   	   is = is+1
	end do
	if (.not.found.or.is.gt.nmesh) then
	  hion1(iprobe) = hion(1,nmesh)
	  hion2(iprobe) = hion(2,nmesh)
	endif
       end do

      call chr(1+npos*2.,0,ct2,lct2)
      ct = '('//ct2(1:lct2)//'(1pe12.4))'
      lct = lench(ct)
       write (30,ct(1:lct)) 
c     : (timtomax,tehist(iprobe),nehist(iprobe),iprobe=1,npos)
     : tns,(hion1(iprobe),hion2(iprobe),iprobe=1,npos)

      goto 1
   2  continue

      close (30)


      nstep = i-1
      end

      
