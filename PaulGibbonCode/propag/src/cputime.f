
      subroutine cputime(t,i1)
      real tar(2)

c  IBM AIX
       call cpu_time(t) 

c  Linux i386, GNU f77 
c      tuser=etime(tar)
c      t=tar(1)

c  Cray t3e 
c       call second(t1)
c       t=t1
      end
