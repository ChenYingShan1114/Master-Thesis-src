
           
c  ================================================================
c
c     em spectral probe
c
c  ================================================================

      subroutine spec(dtsamp)
      include 'dynafoc.h'

      integer iwk(6*nfmax+150),ipow2(15)
      real wk(6*nfmax+150),work1(0:nfmax),xw(0:nfmax)
      complex apr(nfmax)
      data ipow2/1,2,4,8,16,32,64,128,256,512,1024,2048,4196,8192
     :,16384/
      data nsamp,nplot/1,1/
      save nsamp,nplot


c  probe position  kprobe  follows plasma grid (zdens)
c  - gets reset in REZONE every time zdens(1) reaches LHB
      if (nsamp.eq.1) then
        nsmax = amin1(nfmax*1.,zl/dtsamp)
        write (15,*) 'Started spectral sampling at t = ',ctime
        write (15,*) '# samples: ',nsmax
      endif

      apr(nsamp) = a(0,kprobe)
      nsamp = nsamp+1

c  do FT when probe gone right across grid      
      if (nsamp.gt.nsmax) then

        call fftcc(apr,nsamp,iwk,wk)

        dw = 2*pi/nsamp/dtsamp
        uwmax = 0.
        n2 = nsamp**2
        nso2 = nsamp/2

        do i=1,nso2
c  +ve frequencies
          xw(nso2+i) = dw*i
          work1(nso2+i) = abs(apr(i+1))**2/nsamp**2
c  -ve frequencies
          xw(nso2-i) = -dw*i
          work1(nso2-i) = abs(apr(nsamp-i))**2/nsamp**2
          uwmax = amax1(uwmax,work1(nso2+i),work1(nso2-i))
        end do

        xw(nso2) = 0.
        work1(nso2) = abs(apr(1))**2/nsamp**2
        uwmax = amax1(uwmax,work1(nso2))

        work1(0)=uwmin

        ngf=1+nsamp/5000

c  spectrum limits
        if (omemin.eq.omemax) then
          wmax = 2*pi/tau
          wmin = -omemax
        else
          wmin = omemin
          wmax = omemax
        endif
        
        call grxy(xw,work1,nsamp-1,7000+nplot,ngf,-1
     :,' w             ',' F(w)**2       ','spec'//ctime(1:12)
     :,wmin,wmax,wmax/2.,0.,uwmax,10.)

        nsamp=1
        nplot=nplot+1
        kprobe = nzp
      endif
        

      end
