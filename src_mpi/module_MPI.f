
      MODULE MPIvars

      INTEGER ::  mynode,nprocs,NDP(3)
     &,nabor_node(8),NNDP,NPD_i,NPD_j,NPD_k
     & ,node_cell_min(3)
     & ,node_cell_max(3),nncell(3),nncell_lj(3)
     & ,NP_node,NPtot_node,NPtot3_node
     & ,NPmax_node,NPtot4_node
     &,N_atom_passLL(3),N_atom_passRR(3)
     &,N_atom_passDD(3),N_atom_passUU(3),N_atom_passLD(3)
     &,N_atom_passLU(3),N_atom_passRD(3),N_atom_passRU(3)
     &,N_atom_recvLL(3),N_atom_recvRR(3)
     &,N_atom_recvDD(3),N_atom_recvUU(3),N_atom_recvLD(3)
     &,N_atom_recvLU(3),N_atom_recvRD(3),N_atom_recvRU(3)
     &,N_force_passLL(3),N_force_passRR(3)
     &,N_force_recvLL(3),N_force_recvRR(3)
     &,N_force,N_force_corner,N_force_UD,N_force_RL
     &,N_force_passDD(3),N_force_passUU(3),N_force_passLD(3)
     &,N_force_passLU(3),N_force_passRD(3),N_force_passRU(3)
     &,N_force_recvDD(3),N_force_recvUU(3),N_force_recvLD(3)
     &,N_force_recvLU(3),N_force_recvRD(3),N_force_recvRU(3)
     &,NDD,NUU,NRU,NLU,NRR,NRD,NLD,NLL

      REAL*8 :: RLIMAX,RLJMAX,fload_node(3)
     &,pmax_node(3),pmin_node(3)


      END MODULE MPIvars
