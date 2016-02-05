       SUBROUTINE allocate_array
      
       USE MPIvars 
       USE dyn_array
       USE PARAMS
       USE STRUCTR 
       USE specif 
       USE POTS
       USE ANALYSIS
       USE MDSTP
!
       IMPLICIT none
   
       INTEGER :: i

       REAL*8 PBE(3),density,PBE2(3)
       INTEGER NPMAX2,NbufMAX2,NCMAXlj

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C assign value to variables
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C NTMAX : total # of atoms in the system
C       NTMAX=np
       PBE=0.0d0
       density=2.0d0*float(NP)/(CUBE_r(1)*CUBE_r(2)*CUBE_r(3))
C       WRITE(*,*) 'density=',density
C Estimate # of atoms in the buffer region
       IF (NNDP.EQ.1) THEN
         DO i=1,3
           IF (i.EQ.NPD_i) THEN
             PBE(i)=RLJMAX
             GOTO 10
           ENDIF
           PBE(i)=CUBE_r(i)
10         CONTINUE
         ENDDO
         NbufMAX=INT(density*2.0d0*PBE(1)*PBE(2)*PBE(3))
C         WRITE(*,*)'NbufMAX=',NbufMAX,density,PBE(1),PBE(2),PBE(3)
C NPMAX : # of atoms per node
         NPMAX=INT(1.10d0*(NP/nprocs))+NbufMAX
         WRITE(*,*) 'NPMAX',NPMAX
         NPMAX2=INT((NP/nprocs))+
     &          INT(NbufMAX*(2.0d0*RLIMAX/RLJMAX))
C NBMAX1: # of atoms per buffer region per direction (LL,RR)
         NBMAX1=INT(0.5d0*NbufMAX)
C Increase NBMAX1 but real problem may be else where
C -travisk
C         NBMAX1=INT(1.d0*NbufMAX)
       ELSE
C Estimate # of atoms in the buffer region
         DO i=1,3
           IF (i.EQ.NPD_i) THEN
             PBE(i)=(CUBE_r(i)/NDP(i))+2.0d0*RLJMAX
             PBE2(i)=(CUBE_r(i)/NDP(i))+4.0d0*RLIMAX
             GOTO 20
           ENDIF
           IF (i.EQ.NPD_j) THEN
             PBE(i)=(CUBE_r(i)/NDP(i))+2.0d0*RLJMAX
             PBE2(i)=(CUBE_r(i)/NDP(i))+4.0d0*RLIMAX
             GOTO 20
           ENDIF
           PBE(i)=CUBE_r(i)
           PBE2(i)=CUBE_r(i)
20         CONTINUE
         ENDDO
         NbufMAX=INT(density*((PBE(NPD_i)*PBE(NPD_j))-
     &              ((CUBE_r(NPD_i)/NDP(NPD_i))*
     &              (CUBE_r(NPD_j)/NDP(NPD_j))))*PBE(NPD_k))
         NbufMAX2=INT(density*((PBE2(NPD_i)*PBE2(NPD_j))-
     &              ((CUBE_r(NPD_i)/NDP(NPD_i))*
     &              (CUBE_r(NPD_j)/NDP(NPD_j))))*PBE2(NPD_k))
C       WRITE(*,*)'NbufMAX=',NbufMAX,NP,nprocs
C NPMAX : # of atoms per node
         NPMAX=INT(1.10d0*(NP/nprocs))+NbufMAX
!         WRITE(*,*) 'NPMAX =',NPMAX,NbufMAX,NP 
         NPMAX2=INT((NP/nprocs))+NbufMAX2
C NBMAX1: # of atoms per buffer region per direction (LL,RR)
         NBMAX1=INT(density*(CUBE_r(NPD_j)/NDP(NPD_j))*RLJMAX*
     &              CUBE_r(NPD_k))
C NBMAX2: # of atoms per buffer region per direction (DD,UU)
         NBMAX2=INT(density*(CUBE_r(NPD_i)/NDP(NPD_i))*RLJMAX*
     &              CUBE_r(NPD_k))
C NBMAX3: # of atoms per buffer region per direction (LD,RD,LU,RU)
         NBMAX3=INT(density*RLJMAX*RLJMAX*CUBE_r(NPD_k))
C         WRITE(*,*) NPMAX,NBMAX1,NBMAX2,NBMAX3
       ENDIF
C NCMAX: total # of REBO cells in the system
!       WRITE(*,*) nncell(:)
       NCMAX=nncell(1)*nncell(2)*nncell(3)
C NCMAXlj: total # of LJ cells in the system
       NCMAXlj=nncell_lj(1)*nncell_lj(2)*nncell_lj(3)
C NPM1 : # of atoms per node +1
       NPM1=NPMAX+1
C NLMAX : # of atoms in the REBO neighbor list
       NLMAX=NPMAX2*INT((1.5d0*density*3.14*(RLIMAX**3)*4)/3-1)
C NMABIG : # of atoms in the LJ neighbor list
C -travisk increase the size to prevent error 
       NMABIG=INT((NP/nprocs))*
     &        INT((1.0*density*3.14*(RLJMAX**3)*4)/3-1)
C     &        INT((0.5*density*3.14*(RLJMAX**3)*4)/3-1)
C NSMAX: total # of MD steps
       NSMAX=KVC
C Nsendtime: # of times to send read coord.d to sub-nodes
       Nsendtime=INT(NP/NSendMAX)
       IF (NP.LT.NsendMAX) THEN
         Nsendtime=1
       ENDIF       
!       WRITE(*,*) 'allocate',NPMAX,NLMAX,NMABIG,NCMAX
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C ALLOCATE arrays
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C from common_ch.inc
       ALLOCATE (NABORS(NPM1),nclist_node(NPMAX),nchead_node(NCMAX))
       ALLOCATE (icneigh(27*NCMAX),n_icell_node(NPMAX))
       ALLOCATE (XHC(NPMAX,RTYPES))
       ALLOCATE (WW(NLMAX),DWW(NLMAX),RCOR(NLMAX))
       ALLOCATE (COR(NLMAX,3))
       ALLOCATE (LIST(NLMAX),LCHECK(NLMAX))
       ALLOCATE (EXX1(NLMAX),DEXX1(NLMAX),EXX2(NLMAX))

C from common_md.inc
       ALLOCATE (R0_node(3,NPMAX),R1_node(3,NPMAX),R2_node(3,NPMAX),
     &           R3_node(3,NPMAX))
       ALLOCATE (KTYPE_node(NPMAX),itr_node(NPMAX),NA(NPMAX))
       ALLOCATE (RNP(3,NP),RNP_node(3,NPMAX),
     &           RNP_lj(3,NPMAX),RPP(3,NLMAX),RNP_buf(3,NP))
       ALLOCATE (R0L(3,NPMAX),nlist(NP))
       ALLOCATE (GL(3*NP))
C       ALLOCATE (eatom(NP))
       ALLOCATE (HDEL(NSMAX),EN(NSMAX),ETOT(NSMAX))
C       ALLOCATE (XKE(NLMAX))
      ALLOCATE ( noa_node(NTYPES)
     & )
C from commmon_lj_new.inc
       ALLOCATE (ljlist_node(NPMAX),ljhead_node(NCMAXlj))
       ALLOCATE (ljnbc(27*NCMAXlj),n_icell_lj(NPMAX),ljnabr(NPM1))
       ALLOCATE (ljvlit(NMABIG))

C from passing and receiving
       ALLOCATE (np_send(Nsendtime),np_remain(Nsendtime))
       ALLOCATE (buf_xmol(5*NPMAX),buf_xmol_recv(5*NPMAX))
       ALLOCATE (buf_cor(15*NPMAX),buf_cor_recv(15*NPMAX))
       IF (NNDP.EQ.1) THEN
         ALLOCATE (buf_atom_passLL(15*NBMAX1),
     &             buf_atom_passRR(15*NBMAX1))
         ALLOCATE (buf_force_LL(4*NBMAX1),buf_force_RR(4*NBMAX1))
       ELSE
         ALLOCATE (buf_atom_passLL(15*NBMAX1),
     &             buf_atom_passRR(15*NBMAX1),
     &             buf_atom_passDD(15*NBMAX2),
     &             buf_atom_passUU(15*NBMAX2))
         ALLOCATE (buf_atom_passLD(15*NBMAX3),
     &             buf_atom_passLU(15*NBMAX3),
     &             buf_atom_passRD(15*NBMAX3),
     &             buf_atom_passRU(15*NBMAX3))
         ALLOCATE (buf_force_LL(4*NBMAX1),buf_force_RR(4*NBMAX1),
     &             buf_force_DD(4*NBMAX2),buf_force_UU(4*NBMAX2))
         ALLOCATE (buf_force_LD(4*NBMAX3),buf_force_LU(4*NBMAX3),
     &             buf_force_RD(4*NBMAX3),buf_force_RU(4*NBMAX3))
       ENDIF
       ALLOCATE (buf_atom_recv(15*NbufMAX))
       ALLOCATE (buf_force_recv(4*NbufMAX))
!    
! Allocate structural information which was in setin
!  -travisk
      ALLOCATE(
     & KA(2*NSendMAX)
     &,KTYPE(2*NSendMAX)
     &,R0(3,2*NSendMAX) 
     &,R1(3,2*NSendMAX)
     &,R2(3,2*NSendMAX)
     &,R3(3,2*NSendMAX)
     &,R4(3,2*NSendMAX)
     &,ITR(2*NSendMAX)
C     &,MLIST(2*NSendMAX)
C     &,R0L(3,2*NSendMAX) !in allocate_arrrays
     &)
! Allocate analysis arrays
      IF ( mynode.EQ.0 ) THEN
      ALLOCATE(
     &  MOLSTi(NP),MPNTi(NP),ASUBi(NP)
     & ,MOLSTo(NP),MPNTo(NP),ASUBo(NP)
     & ,NABORSi(NLMAX),LISTi(NLMAX)
     & ,NABORSo(NLMAX),LISTo(NLMAX),RECTI(NP)
     & ,KTYPEi(NP),ITRi(NP)
     & ,ELISTi(NP),APLIST(NP),ELISTo(NP),ARLIST(NP)
     & ,R0o(3,NP),R0i(3,NP),R1i(3,NP)
     & ,R2i(3,NP),R3i(3,NP)
     & ,RCLIST(NP),IMOLo(NP),IMOLi(NP),MIMOLi(NP),SIMOLi(NP)
     & ,SMPNTi(NP),SMOLSTi(NP),SMIMOLi(NP)
     & )
      ENDIF

C       WRITE(*,*) 'finish allocate_array'
       RETURN
       END

