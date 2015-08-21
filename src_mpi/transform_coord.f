       PROGRAM trans_coord

C This program transform coord.d format for MPI_FCH_version3_v04
C In this version of parallel REBO,
C local node only know local information
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C Local variable
       PARAMETER (NPMAX=100000)
       INTEGER NP,K,KTYPE(NPMAX),itr(NPMAX)
       REAL*8 R0(3,NPMAX),R1(3,NPMAX),R2(3,NPMAX),R3(3,NPMAX),
     & TTIME,DELTA,CUBE(3)
       CHARACTER*80 HEAD
C Read original coord.d file
       READ(*,*) HEAD
       READ(*,*) NP
       READ(*,*) TTIME,DELTA
       READ(*,*) CUBE(1),CUBE(2),CUBE(3)
       DO i=1,NP
         READ(*,*) K,KTYPE(K),(R0(N,K),N=1,3),itr(K)
       ENDDO
       DO i=1,NP
         READ(*,*) K,(R1(N,K),N=1,3)
       ENDDO
       DO i=1,NP
         READ(*,*) K,(R2(N,K),N=1,3)
       ENDDO
       DO i=1,NP
         READ(*,*) K,(R3(N,K),N=1,3)
       ENDDO
C Write in ner format coord.d
       WRITE(*,100) HEAD
       WRITE(*,200) NP
       WRITE(*,300) TTIME,DELTA
       WRITE(*,400) CUBE(1),CUBE(2),CUBE(3)
       DO i=1,NP
         K=i
         WRITE(*,500) K,KTYPE(K),(R0(N,K),N=1,3),itr(K)
         WRITE(*,600) K,(R1(N,K),N=1,3)
         WRITE(*,600) K,(R2(N,K),N=1,3)
         WRITE(*,600) K,(R3(N,K),N=1,3)
       ENDDO
100    FORMAT(A80)
200    FORMAT(I10)
300    FORMAT(2(E20.11))
400    FORMAT(3(E20.11))
500    FORMAT(I10,2x,I2,3E20.11,I3)
600    FORMAT(I10,4x,3E20.11)
C       WRITE(*,*) 'finish transform coord.d'
       CALL write_error
       END

