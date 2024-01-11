

c ======================

      subroutine powcon

c ======================

c  writes out 2D power contours for postprocessing by CPLOT



      include 'dynafoc.h'
      character chead*80,cp*80,csym*40,cxmin*80,cxmax*80,cymin*80
     :,cmima*60,cfile*80,csnap*80,cymax*80
      real xt(0:ntmax)

      nhist = nt/ihist

      if (mfoc.eq.1) then
        zoff = Rfoc
      else
        zoff = 0.
      endif

      do i=1,nhist
        xt(i) = (i*dt*ihist - zoff)*rfac
      end do


      do icon = 1,ncon
c  fetch data
        call chr(1.*icon,0,csym,lc)
        cfile='pow'//csym(1:lc)
        open(51,file=cfile)
        do i=1,nhist
          write (51,'(2f15.4)') xt(i),rcon(icon,i)*rfac
        end do

        close(51)
      end do


      end
