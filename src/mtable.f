      subroutine mtable


      USE POTS
      USE TABS

      IMPLICIT none

      INTEGER :: KI,KJ,I
      REAL*8 :: VA,DVA,VB,DVB,VC,DVC,VV,DVV
     & ,RC,RSQ,FC,DFC,DTEMP,FF1,FF2,DF1,DF2,DVM


C
C generate lookup tables for bond-order potentials
C

      do ki=1,RTYPES
            do kj = ki,RTYPES
C set table divisions based on outer cutoff
                 ddtab(ki,kj) = rb2(ki,kj)/float(ntab - 2)
                 ddtab(kj,ki) = ddtab(ki,kj)
                 rc = 0.0d0
C Loop over all vales of r from 0 - r_{cut}^{outer}
                 do i=2,ntab-1
                   if(ddtab(ki,kj).ne.0.0d0) then
                      rc = rc + ddtab(ki,kj)
                      rsq = rc*rc
c cut-off function
                      FC=0.0d0
                      DFC=0.0d0
C 
                      IF(RC.LT.RB2(KI,KJ)) THEN
                           DTEMP=PID(KI,KJ)*(RC-RB1(KI,KJ))
                           FC=(1.0d0+COS(DTEMP))/2.0d0
                           DFC=-PID(KI,KJ)/2.0d0*SIN(DTEMP)
                      ENDIF
C
                      IF(RC.LE.RB1(KI,KJ)) THEN
                           FC=1.0d0
                           DFC=0.0d0
                      ENDIF
C Cut off function
                      tabfc(ki,kj,i) = fc
                      tabfc(kj,ki,i) = tabfc(ki,kj,i)
                      tabdfc(ki,kj,i) = dfc
                      tabdfc(kj,ki,i) = tabdfc(ki,kj,i)

c attractive pair terms
                      VA=AD(KI,KJ)*EXP(-AXL(KI,KJ)*RC)
                      DVA=-AXL(KI,KJ)*VA
C
                      VB=BD(KI,KJ)*EXP(-BXL(KI,KJ)*RC)
                      DVB=-BXL(KI,KJ)*VB
C
                      VC=CD(KI,KJ)*EXP(-CXL(KI,KJ)*RC)
                      DVC=-CXL(KI,KJ)*VC
C
                      VV=(VA+VB+VC)/2.0D0        !1/2 for bij ave
                      DVV=(DVA+DVB+DVC)/2.0D0    !1/2 for bij ave
                      atable(ki,kj,i) = FC*VV
                      atable(kj,ki,i) = atable(ki,kj,i)
                      datable(ki,kj,i) = (FC*DVV+DFC*VV)/RC
                      datable(kj,ki,i) = datable(ki,kj,i)
c repulsive pair terms
                      FF1=DD(KI,KJ)*EXP(-DXL(KI,KJ)*RC)
                      DF1=-DXL(KI,KJ)*FF1
C
                      FF2=(1.0D0+ED(KI,KJ)/RC)
                      DF2=-ED(KI,KJ)/RSQ
C
                      VV=FF1*FF2
                      DVM=(DF1*FF2 + FF1*DF2)
                      rtable(ki,kj,i) = vv*fc
                      rtable(kj,ki,i) = rtable(ki,kj,i)
                      drtable(ki,kj,i) = -(FC*DVM+DFC*VV)/RC
                      drtable(kj,ki,i) = drtable(ki,kj,i)
                   else
                      tabfc(ki,kj,i)=0.0d0
                      tabfc(kj,ki,i) = tabfc(ki,kj,i)
                      tabdfc(ki,kj,i) = 0.0d0
                      tabdfc(kj,ki,i) = tabdfc(ki,kj,i)
                      atable(ki,kj,i) = 0.0d0
                      atable(kj,ki,i) = atable(ki,kj,i)
                      datable(ki,kj,i) = 0.0d0
                      datable(kj,ki,i) = datable(ki,kj,i)
                      rtable(ki,kj,i) = 0.0d0
                      rtable(kj,ki,i) = rtable(ki,kj,i)
                      drtable(ki,kj,i) = 0.0d0
                      drtable(kj,ki,i) = drtable(ki,kj,i)
                   endif 
                 enddo

                 atable(ki,kj,1) = atable(ki,kj,2)
                 atable(kj,ki,1) = atable(ki,kj,1)
                 datable(ki,kj,1) = datable(ki,kj,2)
                 datable(kj,ki,1) = datable(ki,kj,1)
                 rtable(ki,kj,1) = rtable(ki,kj,2)
                 rtable(kj,ki,1) = rtable(ki,kj,1)
                 drtable(ki,kj,1) = drtable(ki,kj,2)
                 drtable(kj,ki,1) = drtable(ki,kj,1)
                 tabfc(ki,kj,1) = 0.0d0
                 tabfc(kj,ki,1) = 0.0d0
                 tabdfc(ki,kj,1) = 0.0d0
                 tabdfc(kj,ki,1) = 0.0d0

                 atable(ki,kj,ntab) = 0.0d0
                 atable(kj,ki,ntab) = atable(ki,kj,ntab)
                 datable(ki,kj,ntab) = 0.0d0
                 datable(kj,ki,ntab) = datable(ki,kj,ntab)
                 rtable(ki,kj,ntab) = 0.0d0
                 rtable(kj,ki,ntab) = rtable(ki,kj,ntab)
                 drtable(ki,kj,ntab) = 0.0d0
                 drtable(kj,ki,ntab) = drtable(ki,kj,ntab)
                 tabfc(ki,kj,ntab) = 0.0d0
                 tabfc(kj,ki,ntab) = tabfc(ki,kj,ntab)
                 tabdfc(ki,kj,ntab) = 0.0d0
                 tabdfc(kj,ki,ntab) = tabdfc(ki,kj,ntab)

            enddo
       enddo
C
       return
       end

