  h  )   k820309    9          19.0        ���]                                                                                                          
       /home/albd/research/FHIaims/running_copy/FHIaims/src/ipc/ipc_hamiltonian_and_ovl_transfer.f90 IPC_HAMILTONIAN_AND_OVL_TRANSFER                      @                              
                            @                              
                                                           
                                                           
                                                           
       IPC_WRITE_ARRAY IPC_START_TRANSACTION IPC_READ_ARRAY                                                     
       CHECK_ALLOCATION #         @                                                      #A    #B 	             
                                                    1        @                                  	                    
     p          1     1                   %         @                                
                           #TRANSACTION_NAME                                                                  1 #         @                                                      #A    #B              
                                                    1        @                                                      
     p          1     1                   #         @                                                      #INFO    #NAME    #CALLER    #DIMENSION1    #DIMENSION2    #DIMENSION3    #DIMENSION4              
                                                       
                                                    1           
                                                   1           
                                                      
                                                      
                                                      
                                                       @                                                      @                                                                   &                   &                                                                                                                                                                                  &                   &                                                                                                                                                                                                                                                                                                                        0                                                                  @@                                                  
                &                   &                                           #         @                                   !                     #         @                                  "                                                               #         @                                  #                    #I_K_POINT $             D @                               $            #         @                                  %                    #HAMILTONIAN_W_COMPLEX &            
                                 &                     	     p          5 r    '  5 r    n                                           1n                                      2p            5 r    '  5 r    n                                      1n                                      2  5 r          5 r    '  5 r    n                                      1n                                      2  5 r                                               #         @                                  '                    #OVERLAP_MATRIX_W_COMPLEX (            
                                 (                        p            5 r    '  5 r    n                                       1n                                      2      5 r    '  5 r    n                                      1n                                      2                                        �   �      fn#fn    '  @   j   DIMENSIONS    g  @   J   PBC_LISTS     �  @   J   RUNTIME_CHOICES    �  @   J   PHYSICS    '  u   J  IPC    �  Q   J  MPI_TASKS $   �  V       IPC_WRITE_ARRAY+IPC &   C  L   a   IPC_WRITE_ARRAY%A+IPC &   �  �   a   IPC_WRITE_ARRAY%B+IPC *     f       IPC_START_TRANSACTION+IPC ;   y  L   a   IPC_START_TRANSACTION%TRANSACTION_NAME+IPC #   �  V       IPC_READ_ARRAY+IPC %     L   a   IPC_READ_ARRAY%A+IPC %   g  �   a   IPC_READ_ARRAY%B+IPC +   �  �       CHECK_ALLOCATION+MPI_TASKS 0   �  @   a   CHECK_ALLOCATION%INFO+MPI_TASKS 0   �  L   a   CHECK_ALLOCATION%NAME+MPI_TASKS 2     L   a   CHECK_ALLOCATION%CALLER+MPI_TASKS 6   k  @   a   CHECK_ALLOCATION%DIMENSION1+MPI_TASKS 6   �  @   a   CHECK_ALLOCATION%DIMENSION2+MPI_TASKS 6   �  @   a   CHECK_ALLOCATION%DIMENSION3+MPI_TASKS 6   +  @   a   CHECK_ALLOCATION%DIMENSION4+MPI_TASKS &   k  @       N_K_POINTS+DIMENSIONS "   �  �       K_PHASE+PBC_LISTS "   O	  @       N_CELLS+PBC_LISTS %   �	  �       CELL_INDEX+PBC_LISTS #   3
  @       N_BASIS+DIMENSIONS "   s
  @       N_SPIN+DIMENSIONS 5   �
  @       PACKED_MATRIX_FORMAT+RUNTIME_CHOICES (   �
  q       PM_NONE+RUNTIME_CHOICES -   d  @       N_CENTERS_BASIS_I+DIMENSIONS !   �  �       ATK_K_POINT_LIST :   H  H       INITIATE_HAMILTONIAN_AND_OVL_IPC_TRANSFER 2   �  r       SET_K_PHASE_FROM_ATK_K_POINT_LIST ;     W       CALCULATE_AND_TRANSFER_HAMILTONIAN_AND_OVL E   Y  @   a   CALCULATE_AND_TRANSFER_HAMILTONIAN_AND_OVL%I_K_POINT /   �  c       TRANSFER_HAMILTONIAN_TO_MEMORY E   �  f  a   TRANSFER_HAMILTONIAN_TO_MEMORY%HAMILTONIAN_W_COMPLEX +   b  f       TRANSFER_OVERLAP_TO_MEMORY D   �  �  a   TRANSFER_OVERLAP_TO_MEMORY%OVERLAP_MATRIX_W_COMPLEX 