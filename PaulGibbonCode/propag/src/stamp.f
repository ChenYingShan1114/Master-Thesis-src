
c     =========================
c     
c     Time stamp 
c     
c     =========================
      
      subroutine stamp(istream,ibegin)
      character fdate*24
      character cdate*8, ctime*10, czone*5
      integer vals(4)
      integer ibegin

c      call DATE_AND_TIME(cdate,ctime,czone,vals)
       call DATE_AND_TIME(cdate,ctime,czone)

      if (ibegin.eq.1) then

         write(istream,'(///a15,a12/a12,a6,a8,a5)') 
     : 'DRELP  run on'
     :        ,cdate(7:8)//'/'//cdate(5:6)//'/'//cdate(1:4)
     :        ,'Start time: ',ctime(1:2)//':'//ctime(3:4),'GMT',czone

      else 
         write(istream,'(a,a6)') 'Finished run at time: '
     :        ,ctime(1:2)//':'//ctime(3:4)
      endif
c     =========================
c     
c     SUN date and time
c     
c     =========================
      
c     write(istream,'(/a24/)') fdate()

      end


