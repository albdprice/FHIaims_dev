  �     k820309    9          19.0        ���]                                                                                                          
       /home/albd/research/FHIaims/running_copy/FHIaims/src/v_hartree_multipole_evaluation.f90 V_HARTREE_MULTIPOLE_EVALUATION              GET_V_HARTREE_MULTIPOLE_AND_GRADIENT GET_V_NUC_AND_GRADIENT                      @                              
       DP                                                     
       AIMS_STOP                                                     
       COMMUNICATION_TYPE SHMEM_COMM FORCE_NEW_FUNCTIONAL USE_HARTREE_NON_PERIODIC_EWALD                                                     
       MULTIPOLE_RADIUS_SQ OUTER_POTENTIAL_RADIUS L_HARTREE_MAX_FAR_DISTANCE MULTIPOLE_MOMENTS          @       �   @                              
       INITIALIZE_ANALYTIC_MULTIPOLE_COEFFICIENTS CLEANUP_ANALYTIC_MULTIPOLE_COEFFICIENTS                      @                              
       GET_RHO_MULTIPOLE_SPL                      @                              
  	     HARTREE_POTENTIAL_REAL_COEFF UPDATE_OUTER_RADIUS_L HARTREE_FORCE_L_ADD FAR_DISTANCE_HARTREE_FP_PERIODIC_SINGLE_ATOM FAR_DISTANCE_HARTREE_FP_CLUSTER_SINGLE_ATOM_P2 FAR_DISTANCE_REAL_HARTREE_POTENTIAL_SINGLE_ATOM FAR_DISTANCE_REAL_HARTREE_POTENTIAL_SINGLE_ATOM_P2 FAR_DISTANCE_REAL_GRADIENT_HARTREE_POTENTIAL_SINGLE_ATOM FAR_DISTANCE_REAL_GRADIENT_HARTREE_POTENTIAL_SINGLE_ATOM_P2          @       �   @                              
       L_POT_MAX N_ATOMS N_CENTERS_HARTREE_POTENTIAL N_MAX_SPLINE N_MAX_RADIAL N_HARTREE_GRID N_PERIODIC USE_FORCES                                                	     
       N_GRID N_RADIAL                                                
     
       SPECIES EMPTY                                                     
       FREE_POT_ES_SPL                      @                              
       SPLINE_VECTOR_V2 SPLINE_DERIV_VECTOR_V2 VAL_SPLINE VAL_SPLINE_DERIV                      @                              
       CENTERS_HARTREE_POTENTIAL SPECIES_CENTER CENTER_TO_ATOM                                                     
       SPECIES_Z L_HARTREE #         @                                                       #N_POINTS    #POINTS    #V_HARTREE_MP    #V_HARTREE_MP_GRADIENT              
                                                      
                                                     
    p          p          5 � p        r        p          5 � p        r                               D                                                    
     p          5 � p        r        5 � p        r                               F @                                                  
     p          p          5 � p        r        p          5 � p        r                      #         @                                                       #N_POINTS    #POINTS    #V_NUC    #V_NUC_GRADIENT              
                                                      
@ @                                                  
    p          p          5 � p        r        p          5 � p        r                               D                                                    
     p          5 � p        r        5 � p        r                               D                                                    
     p          p          5 � p        r        p          5 � p        r                         �         fn#fn 4     L   b   uapp(V_HARTREE_MULTIPOLE_EVALUATION    k  C   J  TYPES    �  J   J  MPI_TASKS     �  �   J  RUNTIME_CHOICES    �  �   J  PHYSICS 0   "  �   J  ANALYTIC_MULTIPOLE_COEFFICIENTS *   �  V   J  HARTREE_POTENTIAL_STORAGE *     �  J  HARTREE_POTENTIAL_REAL_P0    �  �   j  DIMENSIONS    s  P   J  GRIDS    �  N   J  GEOMETRY      P   J  FREE_ATOMS    a  �   J  SPLINE    �  x   J  PBC_LISTS    ]  T   J  SPECIES_DATA 5   �  �       GET_V_HARTREE_MULTIPOLE_AND_GRADIENT >   @	  @   a   GET_V_HARTREE_MULTIPOLE_AND_GRADIENT%N_POINTS <   �	  �   a   GET_V_HARTREE_MULTIPOLE_AND_GRADIENT%POINTS B   T
  �   a   GET_V_HARTREE_MULTIPOLE_AND_GRADIENT%V_HARTREE_MP K     �   a   GET_V_HARTREE_MULTIPOLE_AND_GRADIENT%V_HARTREE_MP_GRADIENT '   �  �       GET_V_NUC_AND_GRADIENT 0   ]  @   a   GET_V_NUC_AND_GRADIENT%N_POINTS .   �  �   a   GET_V_NUC_AND_GRADIENT%POINTS -   q  �   a   GET_V_NUC_AND_GRADIENT%V_NUC 6   %  �   a   GET_V_NUC_AND_GRADIENT%V_NUC_GRADIENT 