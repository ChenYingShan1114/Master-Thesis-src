c  ====================
c  **  i/o routines  **
c ====================

ci0prnt
      subroutine i0prnt(a,i0)
      parameter(ichan=15)
      character a*6
      write(15,10) a,i0
      write(6,10) a,i0
   10 format(1x,a6,i8)
      return
      end


c     ==========
ci0prn6
      subroutine i0prn6(a,i0)
      character a*6
      write(6,10) a,i0
   10 format(1x,a6,i8)
      return
      end


c     ==========

      subroutine blank6
      write (6,10)
  10  format(/)
      end

c     ==========


      subroutine blank
      parameter(ichan=15)
      write (ichan,10)
  10  format(/)
      end


cr0form
c
c  formatted real number

      subroutine r0form(a,r0,cfm)
      character a*6,cfm*10,cform*20
      if (cfm(1:1).eq.'f') then
        cform='(1x,a6,'//cfm(1:5)//')'
        lfm=13
      else
	cform='(1x,a6,1pe18.8)'
	lfm=15
      endif
      write(15,cform(1:lfm)) a,r0
      write(6,cform(1:lfm)) a,r0
      return
      end

      subroutine ruform(a,r0,cfm,u)
      character a*6,cfm*10,cform*25,u
      if (cfm(1:1).eq.'f') then
        cform='(1x,a6,'//cfm(1:5)//',1x,a6)'
        lfm=19
      else
	cform='(1x,a6,1pe18.8,1x,a6)'
	lfm=21
      endif
      write(15,cform(1:lfm)) a,r0,u
      return
      end

cr0prn6
      subroutine r0prn6(a,r0)
      character a*6
      write(6,10) a,r0
   10 format(1x,a6,10(1pe12.4))
      return
      end

cr0prnt
      subroutine r0prnt(a,r0)
      character a*6
      write(15,10) a,r0
   10 format(1x,a6,10(1pe12.4))
      return
      end


cr1prnt
      subroutine r1prnt(a,r1,m)
      character a*6
      dimension r1(m)
      write(15,10) a,(r1(i),i=1,m)
   10 format(1x/1x,a6/(1x,6(1pe12.3)))
      return
      end


ci1prnt
      subroutine i1prnt(a,i1,m)
      parameter(ichan=15)
      character a*6
      dimension i1(m)
      write(ichan,10) a,(i1(i),i=1,m)
   10 format(1x,a6/(1x,8i8))
      return
      end


ci3prnt
      subroutine i1prnt3(a,i1,m)
      parameter(ichan=15)
      character a*6
      dimension i1(m)
      write(ichan,10) a,(i1(i),i=1,m)
   10 format(1x/1x,a6/(1x,3i8))
      return
      end


      subroutine i2list(a,i1,m1,m2,is,ie)
      parameter(ichan=15)
      character a*6
      dimension i1(m1:m2)
      write(ichan,10) a,(i,i1(i),i=is,ie)
   10 format(1x/1x,a6/(1x,2i8))
      return
      end


      subroutine i3list(a1,i1,a2,i2,m1,m2,is,ie)
      parameter(ichan=15)
      character a1*6,a2*6
      dimension i1(m1:m2),i2(m1:m2)
      write(ichan,10) a1,a2,(i,i1(i),i2(i),i=is,ie)
   10 format(1x/12x,2(2x,a6)/(1x,3i8))
      return
      end

      subroutine r2list(a,r1,r2,m1,m2,is,ie)
      parameter(ichan=15)
      character a*6
      real r1(m1:m2),r2(m1:m2)
      write(ichan,10) a,(i,r1(i),r2(i),i=is,ie)
   10 format(1x/1x,a6/(1x,(i8,2(1pe12.3))))
      return
      end


      subroutine r3list(a,r1,r2,r3,m1,m2,is,ie)
      parameter(ichan=15)
      character a*6
      real r1(m1:m2),r2(m1:m2),r3(m1:m2)
      write(ichan,10) a,(i,r1(i),r2(i),r3(i),i=is,ie)
   10 format(1x/1x,a6/(1x,(i8,3(1pe12.3))))
      return
      end


      subroutine d2list(a,r1,r2,m1,m2,is,ie)
      parameter(ichan=15)
      character a*6
      double precision r1(m1:m2),r2(m1:m2)
      write(ichan,10) a,(i,r1(i),r2(i),i=is,ie)
   10 format(1x/1x,a6/(1x,(i8,2(f12.3))))
      return
      end

      subroutine d3list(a,r1,r2,r3,m1,m2,is,ie)
      parameter(ichan=15)
      character a*6
      double precision r1(m1:m2),r2(m1:m2),r3(m1:m2)
      write(ichan,10) a,(i,r1(i),r2(i),r3(i),i=is,ie)
   10 format(1x/1x,a6/(1x,(i8,3(f12.3))))
      return
      end


      subroutine r4list(a,r1,r2,r3,r4,m1,m2,is,ie)
      parameter(ichan=15)
      character a*6
      integer r1
      dimension r1(m1:m2),r2(m1:m2),r3(m1:m2),r4(m1:m2)
      write(ichan,10) a,(i,r1(i),r2(i),r3(i),r4(i),i=is,ie)
   10 format(1x/1x,a6/(1x,(2i8,3(f12.3))))
      return
      end


cc0prnt
      subroutine c0prnt(a,c0)
      complex c0
      real a
      c0r=real(c0)
      c0i=aimag(c0)
      write(6,10) a,c0r,c0i
   10 format(1x,a6,'   (',1pe12.3,',',1pe12.3,')')
      return
      end


c     c1prnt
      subroutine c1prnt(a,c1,m)
      complex c1
      real a
      dimension c1(m)
      dimension c1r(300),c1i(300)
      do  i=1,m
         c1r(i)=real(c1(i))
         c1i(i)=aimag(c1(i))
      end do
      write(15,10) a,(c1r(i),i=1,m)
      write(15,11) a,(c1i(i),i=1,m)
 10   format(1x/1x,a6,' (real part)'/(1x,6(1pe12.3)))
 11   format(1x/1x,a6,' (imag part)'/(1x,6(1pe12.3)))
      return
      end

cd0prnt
      subroutine d0prnt(a,r0)
      real a,r0
      write(15,10) a,r0
   10 format(1x,a6,10(1pe12.3))
      return
      end


cd1prnt
      subroutine d1prnt(a,r1,m)
      real a,r1
      dimension r1(m)
      write(15,10) a,(r1(i),i=1,m)
   10 format(1x/1x,a6/(1x,6(1pe12.3)))
      return
      end

c  ========================================
c
c      function   PHA
c
c  Returns phase of complex number -pi -> pi
c
c  ========================================

      real function pha(w)
      parameter (pi=3.1415926)
      complex w
      x=real(w)
      y=aimag(w)
      sx=sign(1.,x)
      sy=sign(1.,y)
      itx=1
      if (sx.eq.1.) itx=0
      if ((x.ne.0.).and.(y.ne.0)) then
         pha=atan(y/x)+itx*sy*pi
      else
	 pha=0.
      endif
      return
      end

c  ========================================
c
c         function ZE
c
c   returns exponent of v
c
c  ========================================


      real function ze(v)
      call manex(v,q,ia)
      ze=ia
      end

c     ========================================
c     
c     MANEX
c     
c     returns   mantissa  q
c     exponent  a
c     
c     ========================================

      subroutine manex(z,q,ia)
      double precision xlog,x

      s=sign(1.,z)
      x=abs(z)
      if (x.eq.0) then
         q=0.
         a=0.
	 return
      endif
      xlog=log10(x*1.000001d0)
      ia=xlog
      if (xlog.lt.0) ia=ia-1
      if (ia.ne.xlog) then
         b=xlog-ia
      else
         b=xlog
      endif
      q=s*10**b
      end


c  ========================================
c
c         ICSIZ
c
c  Returns length of character string
c
c  ========================================

      integer function icsiz(cm)
      character cm*40,ch*2
      l=0
20    l=l+1
      ch=cm(l:l)
      if ((ch.ne.'$').and.(l.lt.40)) goto 20
      if (ch.eq.' ') l=0
      icsiz=max(l-1,1)
      end


c     ==========

      subroutine warn(a)
      parameter(ichan=15)
      character a*30
      write (ichan,10) a
      write (6,10) a
  10  format(/1x,a30/1x,'------------------------------'/)
      end


c  ======================================
c
c            CHR
c
c    converts real number to character
c
c  ======================================

      subroutine chr(z,ndp,ch,l)
      character ch*40,chnum*10
      data chnum/'0123456789'/
      real z
c
c   z     = real number
c   ndp   = # decimal points (0-7)
c   ch*40 = returned character string
c   l     = length of string
c
c  multiply up to remove decimal points
      za = max(1.,abs(z))
      ia=log10(za)
      e=float(ia)
      ie=10**ndp
      iz=int(abs(z)*ie+0.5)
c  add - sign for negatives
      is=sign(1.,z)
      l=0
      if (is.lt.0) then
        l=l+1
        ch(l:l)='-'
      endif
c  left of decimal point
      do idp=ia,0,-1
	l=l+1
	iunit=10**(idp+ndp)
	idigit=iz/iunit
	ic=max(min(idigit+1,10),1)
	ch(l:l)=chnum(ic:ic)
        iz=iz-idigit*iunit
      end do
c  decimal point
      if (ndp.gt.0) then
	l=l+1
	ch(l:l)='.'
      endif
c  right of decimal point
      do idp=-1,-ndp,-1
	l=l+1
	iunit=10**(idp+ndp)
	idigit=iz/iunit
	ic=max(min(idigit+1,10),1)
	ch(l:l)=chnum(ic:ic)
        iz=iz-idigit*iunit
      end do
      end


      function lench(c1)
      character c1*80,c0*1
      l=80
  10  l=l-1
      c0=c1(l:l)
      if ((c0.eq.' ').and.(l.gt.1)) goto 10
      if (c0.eq.' ') l=0
      lench=l
      end
c     Random number scrambler
c     =======================
c     
c     called with:
c     x=rano(iseed)
c     
c     returns real number in the interval (0.0 - 1.0)
c     set iseed = -ve integer before using call
c     
c     Example of calling routine:
c     
c     subroutine coords
c     include 'common.h'
c     
c     iseed1 = -11
c     iseed2 = -7
c     
c     
c     do i=1,n
c     x(i)=xlen*rano(iseed1)
c     y(i)=ylen*rano(iseed2)
c     end do
c     
c     
c     end




c
      function rano(idum)
      
      double precision dseed
      real v(97), y
      integer iff, icall
      data iff,icall/0,0/
      save v, y, dseed

      if (idum.lt.0.or.iff.eq.0) then
         iff = 1
         dseed=abs(idum)*1.d0
	 idum=1
	 do j=1,97
            dum=ggubfs(dseed)
         end do

	 do j=1,97
	   v(j)=ggubfs(dseed)
         end do
	 y=ggubfs(dseed)
      endif

c  next index - make sure we don't overstep array bounds if
c  generator returns a 0.0 or 1.0

      j=max(mod(1+int(97.*y),98),1)
      if(j.gt.97.or.j.lt.1) then
	write (6,*) 'Call: ',icall
        write (6,*) 'idum = ',idum,'j = ',j,' y= ',y
	write (6,*) 'Random No. generator not initialised properly'
	write (6,*) 'dummy =',dum,' dseed=',dseed
	write (6,100) (i,v(i),i=1,97)
 100    format (i4,f10.6)
	stop
      endif
c  get next variate and generate new one to fill gap

      y=v(j)
      rano=y
      v(j)=ggubfs(dseed)
      icall = icall + 1

      return
      end

c   imsl routine name   - ggubfs
c
c-----------------------------------------------------------------------
c

c   purpose             - basic uniform (0,1) random number generator -
c                           function form of ggubs
c
c   usage               - function ggubfs (dseed)
c
c   arguments    ggubfs - resultant deviate.
c                dseed  - input/output double precision variable
c                           assigned an integer value in the
c                           exclusive range (1.d0, 2147483647.d0).
c                           dseed is replaced by a new value to be
c                           used in a subsequent call.
c
c   precision/hardware  - single/all
c
c   reqd. imsl routines - none required
c
c   notation            - information on special notation and
c                           conventions is available in the manual
c                           introduction or through imsl routine uhelp
c
c
c-----------------------------------------------------------------------
c
      real function ggubfs (dseed)
c                                  specifications for arguments
      double precision   dseed
c                                  specifications for local variables
      double precision   d2p31m,d2p31
c                                  d2p31m=(2**31) - 1
c                                  d2p31 =(2**31)(or an adjusted value)
      data               d2p31m/2147483647.d0/
      data               d2p31 /2147483648.d0/
c                                  first executable statement
      dseed = dmod(16807.d0*dseed,d2p31m)
      ggubfs = dseed / d2p31
      return
      end




