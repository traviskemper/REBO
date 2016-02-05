       SUBROUTINE deallocate_array
       
       USE dyn_array
       USE MPIvars
       USE SPECIF
       USE STRUCTR
       USE POTS
       USE TABS
       USE PARAMS
       USE SPECIF 
       USE MDSTP
       USE CONTM
       USE ANALYSIS
!
       IMPLICIT none
!
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C DEALLOCATE arrays
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C from common_ch.inc
       DEALLOCATE (NABORS,nclist_node,nchead_node)
       DEALLOCATE (icneigh,n_icell_node)
       DEALLOCATE (XHC)
       DEALLOCATE (WW,DWW,RCOR,COR)
       DEALLOCATE (LIST,LCHECK)
       DEALLOCATE (EXX1,DEXX1,EXX2)

C from commmon_lj_new.inc
       DEALLOCATE (ljlist_node,ljhead_node)
       DEALLOCATE (ljnbc,n_icell_lj,ljnabr)
       DEALLOCATE (ljvlit)

C from passing and receiving
       DEALLOCATE (np_send,np_remain)
       DEALLOCATE (buf_xmol,buf_xmol_recv)
       IF (NNDP.EQ.1) THEN
         DEALLOCATE (buf_atom_passLL,buf_atom_passRR)
         DEALLOCATE (buf_force_LL,buf_force_RR)
       ELSE
         DEALLOCATE (buf_atom_passLL,buf_atom_passRR,buf_atom_passDD,
     &               buf_atom_passUU)
         DEALLOCATE (buf_atom_passLD,buf_atom_passLU,buf_atom_passRD,
     &               buf_atom_passRU)
         DEALLOCATE (buf_force_LL,buf_force_RR,buf_force_DD,
     &               buf_force_UU)
         DEALLOCATE (buf_force_LD,buf_force_LU,buf_force_RD,
     &               buf_force_RU)
       ENDIF
       DEALLOCATE (buf_atom_recv)
       DEALLOCATE (buf_force_recv)

C from common_md.inc
       DEALLOCATE (R0_node,R1_node,R2_node,R3_node)
       DEALLOCATE (KTYPE_node,itr_node,NA)
       DEALLOCATE (RNP,RNP_node,RNP_buf,RNP_lj)
       DEALLOCATE (R0L,NLIST)
       DEALLOCATE (GL)
       DEALLOCATE (HDEL,EN,ETOT)
C       DEALLOCATE (XKE)      
      DEALLOCATE(AD,AXL,BD,BXL,CD,CXL,DD,DXL,ED,RB1,RB2,PID
     &,RMAX,RMIN,CHI,RLIST)
      DEALLOCATE(DDTAB,tabfc,tabdfc,atable,datable,rtable,drtable )

C Allocate noa since it is initialized here
C -travisk
      DEALLOCATE( NOA )
!
! LJ arrays
      DEALLOCATE(
     & XM,XMM,XMMS,RSPL,RMAXLJ,RSLJ,RSPLS,C2,vlook
     &,dlook,SIG,EPS,iv,jv,TAU,EPSS,SIGS )
!       
      IF( mynode.EQ.0 )THEN   
      DEALLOCATE( 
     &   BRB,FRB,EXB,BRBDP,FRBDP,EXBDP,BRBGS
     &   ,FRBGS,EXBGS,BID
     &   )
      DEALLOCATE( GSA,SBA,DSBA,MXDEP )
      ENDIF 
!       WRITE(*,*) 'finish DEALLOCATE_array'
       RETURN
       END
