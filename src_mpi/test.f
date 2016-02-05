       PROGRAM test
       
C test ANINT
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       REAL*8  r,cube

       WRITE(*,*) 'input r='
       READ(*,*) r
       WRITE(*,*) 'input cube='
       READ(*,*) cube

       r=r-cube*ANINT(r/cube)

       WRITE(*,*)'final r=',r

       CALL write_error
       END

