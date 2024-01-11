

c     =============================================================
      

      subroutine arinit
      include 'dynafoc.h'


      if (nz.gt.nzmax .or. nr.gt.nrmax) then
         call warn('!!!! Not enough array space for fields  ')
         stop
      endif
      
c     plasma grid
      do k=0,nz+1
c     initial z-positions
         zdens(k) = k*dz
         z=k*dz
         do j=0,nr+1
c     fixed grid quantities
            rhoe(j,k) = 0.
            if (mioniz.eq.1) then
               do is = 1,nion
                  rhoi(j,k,is) = 0.
               end do
              if (z.ge.zedge .and. z.lt.zend) then
                 rhoatom(j,k) = rho0
                 rhoi(j,k,0) = rho0
              else
                 rhoatom(j,k) = 0.
                 rhoi(j,k,0) = 0.
              endif
            else
               rhoi(j,k,1) = 0.
            endif
            Er(j,k)=0.
            epond(j,k) = 0.

         end do
      end do

c     em grid
      do k=0,nz+1
         do j=0,nr+1
            a(j,k) = (0.,0.)
         end do
      end do

      end
