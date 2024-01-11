
      

c =============================================================
 

      subroutine pulinit
      include 'dynafoc.h'
      complex w
      
      write (15,*) 'Initial phase',xsi

      if (mfoc.eq.1) then

c  Focussing lens - replace Rfoc by z 
c   radius of curvature = z in far field
        do k=1,nz
         z=k*dz
c  focal length
         zf = Rfoc - (z - t0)  

c  local spot size and intensity at zf
	 sig = sigma0*sqrt(1.+Rfoc**2/Zray**2)
         af = a0*sigma0/sig
c  radius of curvature of phase front
         Rc = zf + Zray**2/zf
c  axial phase
         w = cmplx(Zray,-zf)
         xsi = pha(w)

         do i=0,nr
          argr = min( r(i)**2/2/sig**2,20.)
          argz = min( (z-t0)**2/2/tau**2,20.)

          argf = r(i)**2/2./Rc + xsi
          a(i,k) = af*cexp(-argr - argz + yi*argf)
         end do
        end do

      else
c  Flat-phase Gaussian - initialised in centre of grid
        do k=1,nz
         z=k*dz
         do i=0,nr
          argr = min( r(i)**2/2/sigma0**2,20.)
          argz = min( (z-t0)**2/2/tau**2,20.)

          a(i,k) = a0*exp(-argr - argz)
         end do
        end do
      endif

      end
