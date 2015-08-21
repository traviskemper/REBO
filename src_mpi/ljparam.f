      SUBROUTINE LJPARAM

      USE dyn_array
      USE STRUCTR
      USE MDSTP
      USE PARAMS
      USE POTS
      USE TABS

      IMPLICIT none

      INTEGER :: I,J,K,L,M,N,KI,KJ,IT,kli,llj,II
      REAL*8 :: RR(3),RI(3)
      REAL*8 :: RLIS,RSQ,RC,RT,VV,RP
     &,R,vdw,dvdw,RSQS,dr,r6,vlj,dvlj,kmax


c Evaluate the LJ constants: units are in Ang. and eV

c convert K to eV
      do i=1,ktmax
           eps(i,i)=4.0d0*eps(i,i)/11605.0d0
      enddo
c
      do i=1,ktmax
           do j=1,ktmax
                eps(i,j) = sqrt(eps(i,i)*eps(j,j))
                sig(i,j) = (sig(i,i)+sig(j,j))/2.0d0
                XMM(i,j)=0.0d0
                RSPL(i,j)=0.0d0
           enddo
      enddo
c
c inner and outer points for spline for C-C,H-H,C-H
c           C-F, F-F, F-H
c
      do i=1,RTYPES
           do j=1,RTYPES
                if(eps(i,j).ne.0.0d0) then
                     XMM(i,j)=rb2(i,j)
                     RSPL(i,j) = 0.95d0*sig(i,j)
                endif
           enddo
      enddo
c
c inner and outer points for spline for Si-Si,Ge-Ge,Si-Ge
c
c$$$      do i=4,5
c$$$           do j=4,5
c$$$                if(eps(i,j).ne.0.0d0) then
c$$$                     XMM(i,j)=rb2(i,j)
c$$$c
c$$$                     RSPL(i,j) = 0.95d0*sig(i,j)
c$$$                endif
c$$$           enddo
c$$$      enddo

c
      do i=1,ktmax
           do j=1,ktmax
                if(eps(i,j).ne.0.0d0) then
                     XMMS(i,j) = XMM(i,j) - RLL
                     XM(i,j) = XMM(i,j)
                     RMAXLJ(i,j) = 2.50d0*sig(i,j)
                     RSLJ(i,j) = RMAXLJ(i,j)+RLL
!
                     RMAXLJ(i,j) = RMAXLJ(i,j)*RMAXLJ(i,j)
                     RSLJ(i,j) = RSLJ(i,j)*RSLJ(i,j)
                     XMM(i,j) = XMM(i,j)*XMM(i,j)
                     XMMS(i,j) = XMMS(i,j)*XMMS(i,j)
                     RSPLS(i,j) = RSPL(i,j)*RSPL(i,j)
                     sig(i,j) = sig(i,j)*sig(i,j)
                else
                     RMAXLJ(i,j) = 0.0d0
                     RSLJ(i,j) = 0.0d0
                endif
           enddo
      enddo
c
c find spline parameters for C-C, H-H, and C-H
c            C-F, F-F, F-H

      do i= 1,RTYPES
           do j=1,RTYPES
               if(eps(i,j).ne.0.0d0) then
                    dr = rspl(i,j) - rb2(i,j)
                    r6 = (sig(i,j)/rspls(i,j))**3
                    vlj = eps(i,j)*r6*(r6-1.0d0)
                    dvlj = -eps(i,j)/rspl(i,j)*r6*(12.0d0*r6 - 6.0d0)
                    c2(i,j) = (3.0d0/dr*vlj - dvlj)/dr
                    c3(i,j) = (vlj/(dr**2) -c2(i,j))/dr
               endif
           enddo
      enddo
c
c
c find spline parameters for Si-Si,Ge-Ge,Si-Ge
c
c$$$      do i= 4,5
c$$$           do j=4,5
c$$$
c$$$               if(eps(i,j).ne.0.0d0) then
c$$$                    dr = rspl(i,j) - rb2(i,j)
c$$$                    r6 = (sig(i,j)/rspls(i,j))**3
c$$$                    vlj = eps(i,j)*r6*(r6-1.0d0)
c$$$                    dvlj = -eps(i,j)/rspl(i,j)*r6*(12.0d0*r6 - 6.0d0)
c$$$                    c2(i,j) = (3.0d0/dr*vlj - dvlj)/dr
c$$$                    c3(i,j) = (vlj/(dr**2) -c2(i,j))/dr
c$$$               endif
c$$$           enddo
c$$$      enddo
c

c generate table look up for pair interactions
C
      dellj = 0.001d0
      do i=1,ktmax
           do j=1,ktmax
                 kmax = (sqrt(rmaxlj(i,j)) - xm(i,j))/dellj
                 if(kmax+1.gt.10000) then
                      write(6,*) 'kmax,i,j: ',kmax,i,j
                      include 'close.inc'
                      CALL write_error
                 endif
C
                 do k=kmax+1,10000
                      vlook(k,i,j) = 0.0d0
                      dlook(k,i,j) = 0.0d0
                 enddo
C
                 do k=kmax,1,-1
                           r = (k-1)*dellj + xm(i,j)
                           if(((i.gt.2).or.(j.gt.2)).and.(r.eq.0.0d0))
     &                           r=dellj
                           rsqs = r*r
                           if(rsqs.lt.rmaxlj(i,j)) then
                                if(rsqs.ge.rspls(i,j)) then
                                     r6 = (sig(i,j)/rsqs)**3
                                     vdw = eps(i,j) * r6 * (r6 - 1.0d0)
                                     dvdw = eps(i,j)/rsqs *
     &                                     r6*(12.0d0*r6 - 6.0d0)
                                else
                                     dr = r  - rb2(i,j)
                                     vdw = dr*dr*(dr*c3(i,j) +
     &                                                  c2(i,j))
                                     dvdw = -dr*(3.0d0*dr*c3(i,j)
     &                                      + 2.0d0*c2(i,j))/r
                                endif
                           else
                                vdw = 0.0d0
                                dvdw = 0.0d0
                           endif
                      vlook(k,i,j) =  vdw
                      dlook(k,i,j) = dvdw

                 enddo
            enddo
      enddo
c
      return
      end

