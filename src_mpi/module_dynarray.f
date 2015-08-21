       MODULE dyn_array

C Constants
      INTEGER :: NTMAX,NBMAX1,NBMAX2,NBMAX3,NCMAX,NPM1,
     &NSMAX,Nsendtime,NbufMAX,NLMAX,NMABIG,NPMAX

C from common_ch.inc
         INTEGER,ALLOCATABLE ::NABORS(:),nclist_node(:),nchead_node(:),
     &                         icneigh(:),n_icell_node(:)
         INTEGER,ALLOCATABLE ::LIST(:),LCHECK(:)
         REAL*8,ALLOCATABLE  ::XHC(:,:)
         REAL*8,ALLOCATABLE  ::WW(:),DWW(:),RCOR(:),COR(:,:),
     &                         EXX1(:),DEXX1(:),EXX2(:)

C from common_md.inc
         INTEGER,ALLOCATABLE ::KTYPE_node(:),itr_node(:),NA(:),nlist(:)
     & ,noa_node(:)
         REAL*8,ALLOCATABLE  ::R0_node(:,:),R1_node(:,:),R2_node(:,:),
     &                         R3_node(:,:)
         REAL*8,ALLOCATABLE  ::RNP(:,:),RNP_node(:,:),
     &                         RNP_lj(:,:),RPP(:,:),RNP_buf(:,:)
         REAL*8,ALLOCATABLE  ::HDEL(:),EN(:),ETOT(:)
C         REAL*8,ALLOCATABLE  ::XKE(:)

C from common_lj_new
         INTEGER,ALLOCATABLE ::ljlist_node(:),ljhead_node(:),ljnbc(:),
     &                         n_icell_lj(:),ljnabr(:),ljvlit(:)

C from  passing and receiving
         INTEGER,ALLOCATABLE ::np_send(:),np_remain(:)
         REAL*8,ALLOCATABLE  ::buf_xmol(:),buf_xmol_recv(:)
         REAL*8,ALLOCATABLE  ::buf_cor(:),buf_cor_recv(:)
         REAL*8,ALLOCATABLE  ::buf_atom_passLL(:),buf_atom_passRR(:),
     &                         buf_atom_passDD(:),buf_atom_passUU(:),
     &                         buf_atom_passLD(:),buf_atom_passRD(:),
     &                         buf_atom_passLU(:),buf_atom_passRU(:),
     &                         buf_atom_recv(:)
         REAL*8,ALLOCATABLE  ::buf_force_LL(:),buf_force_RR(:),
     &                         buf_force_DD(:),buf_force_UU(:),
     &                         buf_force_LD(:),buf_force_RD(:),
     &                         buf_force_LU(:),buf_force_RU(:),
     &                         buf_force_recv(:)

       END MODULE dyn_array
