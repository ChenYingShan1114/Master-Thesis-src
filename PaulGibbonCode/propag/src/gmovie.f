

c     ==================================================================
c     
c     2D surface/contour snapshots for (r,z,t) mode
c     
c     ================================================================== 

      subroutine gmovie
      include 'dynafoc.h'
      real surf(0:nrmax,0:nzmax),reno(0:nrmax),work2(0:nrmax)
      character chars(1)*8,chx*15,chy*15,ctitle*16,cfile*12,cid*2,cid2*2
      data itc/0/
      save itc
      if (nz.eq.1) return


      rlno = rl*rfac
      zlno = zl*rfac
      nrplot = rplmax/dr
      tnorm = dt*itime*tfac

      if (mfoc.eq.1) then
         zoff = Rfoc
      else
         zoff = 0.
      endif


c     compute parameters for graphics grid interpolation                         
      if(mod(nzp,ksurf).eq.0) then
         npz=nz/ksurf
      else
         npz=nz/ksurf+1
      endif

      if(mod(nr,isurf).eq.0) then
         npr=nrplot/isurf
      else
         npr=nrplot/isurf+1
      endif
      
c     double up in radial dirn                                                     
      npx=2*npr     
      
c     laser intensity
c     ---------------

      call chr(1.*itc,1,cid,lct)  
      if (itc.lt.10) then
         cid2 = '0'//cid(1:1)
      else
         cid2 = cid
      endif

      cfile='pulse'//cid2//'.snap'
      open(60,file=cfile)
      write (15,'(a,2i6)') cfile,npz,npx
c      write (60,'(a,2i6)') '! '//cfile,npz,npx


      do k=1,nz,ksurf
         kg = (k+ksurf-1)/ksurf 
         do i=nrplot,1,-isurf
            ig=(i+isurf-1)/isurf
            surf(ig,kg) = (abs(a(i,k)))**2*sfac
c             surf(ig,kg) = a(i,k)*conjg(a(i,k))/a0**2
c     write (60,103) zdens(k)*rfac,-i*dr*rfac,surf(ig,kg)
            write (60,105) surf(ig,kg)*xi0
         end do
         do i=0,nrplot,isurf
            ig=(i+isurf-1)/isurf+npr
            surf(ig,kg) = (abs(a(i,k))**2)*sfac
c             surf(ig,kg) = a(i,k)*conjg(a(i,k))/a0**2
c     write (60,103) zdens(k)*rfac,i*dr*rfac,surf(ig,kg)
            write (60,105) surf(ig,kg)*xi0
         end do
c  gnuplot EOL
      write (60,'(a)')
      end do

      close(60)

c     gamma
c     ---------------

      cfile='gamma'//cid2//'.sdat'
c     open(60,file=cfile)
c     write (15,'(a,2i6)') cfile,npz,npx


      do k=1,nz,ksurf
         kg = (k+ksurf-1)/ksurf 
         do i=nrplot,1,-isurf
            ig=(i+isurf-1)/isurf
            surf(ig,kg) = gamma(i,k)
c     write (60,103) zdens(k)*rfac,-i*dr*rfac,surf(ig,kg)
         end do
         do i=0,nrplot,isurf
            ig=(i+isurf-1)/isurf+npr
            surf(ig,kg) = gamma(i,k)
c     write (60,103) zdens(k)*rfac,i*dr*rfac,surf(ig,kg)
         end do
c     write (60,'(a)')
      end do

c     close(60)


c     ion density 
c     -------------------------------

      cfile='rhoi'//cid2//'.snap'
c     open(60,file=cfile)
      write (15,'(a,2i6)') cfile,npz,npr
      write (60,'(a,2i6)') '! '//cfile,npz,npx


      do k=1,nz,ksurf
         kg = (k+ksurf-1)/ksurf 
         do i=nrplot,0,-isurf
            ig=(i+isurf-1)/isurf
            surf(ig,kg) = rhoi(i,k,1)*frion(i)/rho0
c     write (60,103) zdens(k)*rfac,-i*dr*rfac,surf(ig,kg)
         end do
         do i=0,nrplot,isurf
            ig=(i+isurf-1)/isurf+npr
            surf(ig,kg) = rhoi(i,k,1)*frion(i)/rho0
c     write (60,103) zdens(k)*rfac,-i*dr*rfac,surf(ig,kg)
            write (60,105) surf(ig,kg)
         end do
      write (60,'(a)')
      end do


c     close(60)


c     electron density 
c     -------------------------------

      cfile='rhoe'//cid2//'.snap'
      open(60,file=cfile)


      do k=1,nz,ksurf
         kg = (k+ksurf-1)/ksurf 
         do i=nrplot,1,-isurf
            ig=(i+isurf-1)/isurf
            surf(ig,kg) = rhoe(i,k)/rhoe0*rhoe0*den0/rho0
            write (60,103) (itime*dt+k*dz-zoff)*rfac,
     :        -i*dr*rfac,surf(ig,kg)
         end do
         do i=0,nrplot,isurf
            ig=(i+isurf-1)/isurf+npr
            surf(ig,kg) = rhoe(i,k)/rhoe0*rhoe0*den0/rho0
            write (60,103) (itime*dt+k*dz-zoff)*rfac,
     :        i*dr*rfac,surf(ig,kg)
         end do
         write (60,'(a)')
      end do

 103  format (2f12.2,(1pe12.4))
 105  format (1pe12.4)

      close(60)




      itc=itc+1
      end










