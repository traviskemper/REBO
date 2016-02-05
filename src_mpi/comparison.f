       PROGRAM compare
       
C This program compare two files
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Local variable
       PARAMETER (NPMAX=20000)
       INTEGER K1(NPMAX),K2(NPMAX),NP1,NP2
       REAL*8 V1(3,NPMAX),V2(3,NPMAX),D1,D2,D3,Q
       CHARACTER*40 NAME1,NAME2

       WRITE(*,*) 'input the standard file name'
       READ(*,*) name1
       WRITE(*,*) 'input total # of lines in standard file'
       READ(*,*) NP1
       WRITE(*,*) 'input the compared file name'
       READ(*,*) name2
       WRITE(*,*) 'input total # of lines in compared file'
       READ(*,*) NP2
       WRITE(*,*) 'input the critiria Q='
       READ(*,*) Q

       OPEN(11,file=name1,status='old')
       OPEN(12,file=name2,status='old')

       DO i=1,NP1
         READ(11,999) K1(i),V1(1,i),V1(2,i),V1(3,i)
       ENDDO
       DO i=1,NP2
         READ(12,999) K2(i),V2(1,i),V2(2,i),V2(3,i)
       ENDDO
       DO j=1,NP2
         DO i=1,NP1
           IF (K2(j).EQ.K1(i)) THEN
             D1=ABS(V2(1,j)-V1(1,i))
             D2=ABS(V2(2,j)-V1(2,i))
             D3=ABS(V2(3,j)-V1(3,i))
             WRITE(99,999) K2(j),D1,D2,D3
             IF ((D1.GT.Q).OR.(D2.GT.Q).OR.(D3.GT.Q)) THEN
               WRITE(*,999) K2(j),D1,D2,D3
             ENDIF
           ENDIF
         ENDDO
       ENDDO
999    FORMAT(I5,3E20.11)
       CALL write_error
       END

