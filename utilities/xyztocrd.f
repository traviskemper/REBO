      
C     SYNTAX > xyz_to_coord.x < [].xyz > coord.d.[]

      INTEGER NA,NLAYER,PARLAN,PARAMNOR,iii,PARARIGID
     & IDUM,NRA,NLA
      CHARACTER*80 TITLE 
      Dimension R(200000,3),R1(200000,3),R2(200000,3)
     & ,R3(200000,3),R4(200000,3),kt(200000)


       ttime=0.0
       delta=0.2
       cube1=1.0E+14
       cube2=1.0E+14
       cube3=1.0E+14
	IDUM=0
	NRA=0
	NLA=0 

       READ (5,*) NA
       READ (5,3) title 
  3   FORMAT (A80)  


      Do 11 J=1,NA
      READ (5,*)  KT(J), R(J,1),R(J,2),R(J,3)
  11  Continue

      write (6,*) 'title'
      WRITE (6,*) na
      write (6,*) ttime, delta
      write (6,*) cube1, cube2, cube3
      write (6,*) cube1, cube2, cube3

	iii=1
        DO 69 JJ=1,na
       WRITE(6,68) JJ,KT(JJ),R(JJ,1),R(JJ,2),R(JJ,3),iii  
       R1(JJ,1)=0.0
        R1(JJ,2)=0.0
        R1(JJ,3)=0.0   
       WRITE(6,831) JJ,  R1(JJ,1),R1(JJ,2),R1(JJ,3)
       R2(JJ,1)=0.0
        R2(JJ,2)=0.0
        R2(JJ,3)=0.0   
       WRITE(6,831) JJ,  R2(JJ,1),R2(JJ,2),R2(JJ,3)
       R3(JJ,1)=0.0
        R3(JJ,2)=0.0
        R3(JJ,3)=0.0   
       WRITE(6,831) JJ,  R3(JJ,1),R3(JJ,2),R3(JJ,3)

  68    FORMAT(2x,I5,3x,i4,3(3x,F12.6),3x,I5) 
  69     CONTINUE

  
 831   FORMAT(3x,I5, 3(3x,3F15.7))  
				
        stop
         end
