      
c     =================================
c     
c     Write xy graph in odplot format
c     
c     =================================

      subroutine grxy(xx,yy,ndp,id,ig,iman,chx,chy,ctitle)
      real x(50002),y(50002),store(15),xx(0:ndp),yy(0:ndp)
      real xmin,xmax,ymin,ymax,x1,x2,y1,y2,xl,yl,dax,day
      integer iman
      character chars(1)*8,chx*15,chy*15,ctitle*16,cfile*12,cid*1

c     open data file
      call chr(mod(abs(id),10)*1.0,0,cid,lc)
      cfile=ctitle(1:4)//cid//'.xy'
      if (ndp.le.1) then
         write (15,*) 'No data for ',cfile,' id=',id
         return
      endif
      write (15,'(a)') cfile
      open (51,file=cfile)

      if (ig.eq.0) then 
         igx=1
         ifirst=0
         nd = ndp+1
      else
         ifirst=1
         nd = ndp
         igx = ig
      endif

      n=nd/igx
      if (mod(nd,igx).ne.0) n=n+1

      do i=ifirst,nd,igx
         ip=(i+igx-1)/igx + 1-ifirst
         x(ip)=xx(i)
         y(ip)=yy(i)
      end do

      idp=id
      ispage=1
      idash=1

      if (iman.lt.0) then
c     manual axes: ityp=-iman; id=-id
         ityp=-iman
         idp=-id
      else if (iman.gt.0) then
c     auto axes: ityp=iman
         ityp=iman
      endif

c     write out to 'ODPLOT' datafile
c     header
      write (50,101)idp,n,ityp,idash,0,ispage,chx,chy,ctitle
c     manual axes
      if (iman.lt.0) then
         write(50,102) xmin,xmax,dax,ymin,ymax,day
      endif
c     data
      write (51,103) (x(i),y(i),i=1,n)
      close(51)

 101  format (2i6,4i4,2a15,1a16)
 102  format (6(1pe12.4))
 103  format (2(1pe11.4))
      end


      subroutine grxy2(xx,yy,ndp,id,ig,iman,chx,chy,ctitle
     :     ,xmin,xmax,dax,ymin,ymax,day)
      real x(50002),y(50002),store(15),xx(0:ndp),yy(0:ndp)
      real xmin,xmax,ymin,ymax,x1,x2,y1,y2,xl,yl,dax,day
      integer iman
      character chars(1)*8,chx*15,chy*15,ctitle*16,cfile*12,cid*1

c     open data file
      call chr(mod(abs(id),10)*1.0,0,cid,lc)
      cfile=ctitle(1:4)//cid//'.xy'
      if (ndp.le.1) then
         write (15,*) 'No data for ',cfile,' id=',id
         return
      endif
      write (15,'(a)') cfile
      open (51,file=cfile)

      if (ig.eq.0) then 
         igx=1
         ifirst=0
         nd = ndp+1
      else
         ifirst=1
         nd = ndp
         igx = ig
      endif

      n=nd/igx
      if (mod(nd,igx).ne.0) n=n+1

      do i=ifirst,nd,igx
         ip=(i+igx-1)/igx + 1-ifirst
         x(ip)=xx(i)
         y(ip)=yy(i)
      end do

      idp=id
      ispage=1
      idash=1

      if (iman.lt.0) then
c     manual axes: ityp=-iman; id=-id
         ityp=-iman
         idp=-id
      else if (iman.gt.0) then
c     auto axes: ityp=iman
         ityp=iman
      endif

c     write out to 'ODPLOT' datafile
c     header
      write (50,101)idp,n,ityp,idash,0,ispage,chx,chy,ctitle
c     manual axes
      if (iman.lt.0) then
         write(50,102) xmin,xmax,dax,ymin,ymax,day
      endif
c     data
      write (51,103) (x(i),y(i),i=1,n)
      close(51)

 101  format (2i6,4i4,2a15,1a16)
 102  format (6(1pe12.4))
 103  format (2(1pe11.4))
      end
