  �%  K   k820309    9          19.0        h��]                                                                                                          
       /home/albd/research/FHIaims/running_copy/FHIaims/src/friction.f90 FRICTION               FRICTION_CALCULATION ALLOCATE_FRICTION DEALLOCATE_FRICTION FRICTION_TENSOR FRICTION_EIGVECS FRICTION_EIGENVALUES FRICTION_N_ACTIVE_ATOMS FRICTION_ITER_LIMIT FRICTION_ACCURACY_ETOT FRICTION_ACCURACY_EEV FRICTION_ACCURACY_RHO FRICTION_ACCURACY_POTJUMP FRICTION_TEMPERATURE NUMERICAL_FRICTION FRICTION_ACTIVE_ATOMS FRICTION_NUMERIC_DISP OUTPUT_FIRST_ORDER_MATRICES OUTPUT_FRICTION_EIGENVECTORS FRICTION_ATOMS_LIST FRICTION_INDEX_LIST FRICTION_KNOTK FRICTION_READ_MATRICES FRICTION_OUTPUT_SPECTRUM FRICTION_DELTA_TYPE FRICTION_WINDOW_SIZE FRICTION_BROADENING_WIDTH FRICTION_PERTURBATION FRICTION_MAX_ENERGY FRICTION_DISCRETIZATION_LENGTH FRICTION_USE_COMPLEX_MATRICES FRICTION_COUPLING_MATRIX_MODE FRICTION_N_Q_POINTS                      @                              
                                                           
                            @                              
                            @                              
                                                           
                                                           
                            @                              
                                                           
                @       �                             	     
                @       �                             
     
                                                           
                                                           
                            @                              
       MY_K_POINT EIGENVEC EIGENVEC_COMPLEX L_COL L_ROW MXLD MXCOL OVLP OVLP_COMPLEX HAM HAM_COMPLEX                      @                              
       MODULE_IS_DEBUGGED DEBUGPRINT                @  @               A                '(                   #METHOD    #COMM    #DEBUG    #CALCULATE_FORCES    #CALCULATE_SPECTRUM    #DO_RPA    #RPA_ORDERS    #RPA_RESCALE_EIGS    #TS_ENE_ACC    #TS_F_ACC    #N_OMEGA_GRID    #K_GRID_SHIFT    #ZERO_NEGATIVE_EIGVALS    #XC    #TS_D    #TS_SR    #MBD_A     #MBD_BETA !   #VDW_PARAMS_KIND "   #ATOM_TYPES #   #FREE_VALUES $   #COORDS %   #LATTICE_VECTORS &   #K_GRID '   #PARALLEL_MODE (               � $                                                                               �                                       �              Cmbd-rsscs                                                     � $                                                                        �                                         ��������                        � $                                   $                                    �                                                                         � $                                   (                                    �                                         ��������                        � $                                   ,                                    �                                                                         � $                                   0                                    �                                                                         � $                                   4                                    �                                                                         � $                                   8                                    �                                                                         � $                                  @       	  
                           �                     
                 �����ư>        1D-6                � $                                  H       
  
                           �                     
                 H�����z>        1D-7                � $                                   P                                    �                                                     15                � $                                  X         
                           �                     
                       �?        0.5D0                � $                                   `                                    �                                                                         � $                                         d                                     �                                       �              C                                                    � $                                  x         
                           �                     
                       4@        20D0                � $                                  �         
                           �                       
                        �                        � $                                   �         
                           �                     
                       @        6D0                � $                             !     �         
                           �                       
                        �                        � $                             "     
       �                                     �                                       �              Cts                            .           � $                             #            �                             &                                                              � $                             $            �                 
            &                   &                                                      � $                             %            P                
            &                   &                                                      � $                             &            �                
            &                   &                                                       � $                              '                             p          p            p                                                  �                                                          ��������                                             ��������                                             ��������                                � $                             (     
                                           �                                       �              Cauto                                   @ @                              )                                   &                   &                   &                                                    @ @                              *                                   &                   &                   &                                                    @ @                              +                   
                &                   &                                                      @                               ,                                                        -                                                       .     
                                                  /     
                                                  0     
                                                  1     
                   @                              2     
                 @                                 3                                                        4                                                       5     
                                                   6                                                        7                     @ @                               8                                   &                                                    @ @                               9                                   &                                                     @                                 :                                                        ;                                                        <                                                       =                        @                              >     
                   @                              ?     
                   @                              @     
                                                  A     
                                                  B     
                 @                                 C                                                        D                      @                                 E            #         @                                   F                    #CONVERGED_SCF G             
                                  G           #         @                                   H                     #         @                                   I                        �   S      fn#fn    �   �  b   uapp(FRICTION    �  @   j  DIMENSIONS       @   J  RUNTIME_CHOICES    M  @   J  PBC_LISTS    �  @   J  TIMING    �  @   J  PHYSICS      @   J  SPECIES_DATA    M  @   J  LOCALORB_IO    �  @   J  GEOMETRY    �  @   J  CONSTANTS      @   J  MPI_TASKS &   M  @   J  SYNCHRONIZE_MPI_BASIC    �  @   J  BASIS "   �  �   J  SCALAPACK_WRAPPER    k  ^   J  DEBUGMANAGER     �  �     MBD_INPUT_T+MBD '   �	  �   a   MBD_INPUT_T%METHOD+MBD %   ~
  �   a   MBD_INPUT_T%COMM+MBD &   "  �   a   MBD_INPUT_T%DEBUG+MBD 1   �  �   a   MBD_INPUT_T%CALCULATE_FORCES+MBD 3   j  �   a   MBD_INPUT_T%CALCULATE_SPECTRUM+MBD '     �   a   MBD_INPUT_T%DO_RPA+MBD +   �  �   a   MBD_INPUT_T%RPA_ORDERS+MBD 1   V  �   a   MBD_INPUT_T%RPA_RESCALE_EIGS+MBD +   �  �   a   MBD_INPUT_T%TS_ENE_ACC+MBD )   �  �   a   MBD_INPUT_T%TS_F_ACC+MBD -   J  �   a   MBD_INPUT_T%N_OMEGA_GRID+MBD -   �  �   a   MBD_INPUT_T%K_GRID_SHIFT+MBD 6   �  �   a   MBD_INPUT_T%ZERO_NEGATIVE_EIGVALS+MBD #   =  �   a   MBD_INPUT_T%XC+MBD %     �   a   MBD_INPUT_T%TS_D+MBD &   �  �   a   MBD_INPUT_T%TS_SR+MBD &   Z  �   a   MBD_INPUT_T%MBD_A+MBD )     �   a   MBD_INPUT_T%MBD_BETA+MBD 0   �  �   a   MBD_INPUT_T%VDW_PARAMS_KIND+MBD +   l  �   a   MBD_INPUT_T%ATOM_TYPES+MBD ,     �   a   MBD_INPUT_T%FREE_VALUES+MBD '   �  �   a   MBD_INPUT_T%COORDS+MBD 0   `  �   a   MBD_INPUT_T%LATTICE_VECTORS+MBD '     �  a   MBD_INPUT_T%K_GRID+MBD .   �  �   a   MBD_INPUT_T%PARALLEL_MODE+MBD     W  �       FRICTION_TENSOR !     �       FRICTION_EIGVECS %   �  �       FRICTION_EIGENVALUES (   s  @       FRICTION_N_ACTIVE_ATOMS $   �  @       FRICTION_ITER_LIMIT '   �  @       FRICTION_ACCURACY_ETOT &   3  @       FRICTION_ACCURACY_EEV &   s  @       FRICTION_ACCURACY_RHO *   �  @       FRICTION_ACCURACY_POTJUMP %   �  @       FRICTION_TEMPERATURE #   3  @       NUMERICAL_FRICTION &   s  @       FRICTION_ACTIVE_ATOMS &   �  @       FRICTION_NUMERIC_DISP ,   �  @       OUTPUT_FIRST_ORDER_MATRICES -   3   @       OUTPUT_FRICTION_EIGENVECTORS $   s   �       FRICTION_ATOMS_LIST $   �   �       FRICTION_INDEX_LIST    �!  @       FRICTION_KNOTK '   �!  @       FRICTION_READ_MATRICES )   "  @       FRICTION_OUTPUT_SPECTRUM $   K"  @       FRICTION_DELTA_TYPE %   �"  @       FRICTION_WINDOW_SIZE *   �"  @       FRICTION_BROADENING_WIDTH &   #  @       FRICTION_PERTURBATION $   K#  @       FRICTION_MAX_ENERGY /   �#  @       FRICTION_DISCRETIZATION_LENGTH .   �#  @       FRICTION_USE_COMPLEX_MATRICES .   $  @       FRICTION_COUPLING_MATRIX_MODE $   K$  @       FRICTION_N_Q_POINTS %   �$  [       FRICTION_CALCULATION 3   �$  @   a   FRICTION_CALCULATION%CONVERGED_SCF "   &%  H       ALLOCATE_FRICTION $   n%  H       DEALLOCATE_FRICTION 