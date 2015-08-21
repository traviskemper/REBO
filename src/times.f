C Subroutines for calculating run time

      SUBROUTINE TIMEI

      USE TM

      CALL CPU_TIME(TTI)

      END SUBROUTINE TIMEI

C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE TIMEF

      USE TM

      IMPLICIT none
 
      INTEGER :: jtime,jday,jhr,jmin,rssec

      CALL CPU_TIME(TTF)
      jtime=int(TTF-TTI)
      jday=int(jtime/86400)
      jhr=int((jtime-jday*86400)/3600)
      jmin=int((jtime-jhr*3600)/60)
      rssec=((jtime-jmin*60))
      WRITE(*,*)''
      WRITE(*,*)'Normal Termination'
      WRITE(*,102) 'Total time:',jday,'days',jhr,'hours',jmin,
     &                  'minutes',rssec,'seconds'

 102  FORMAT(a16,4(i3,a8,1x))


      END SUBROUTINE TIMEF

      
