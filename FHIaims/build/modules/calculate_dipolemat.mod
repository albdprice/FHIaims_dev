  µh  ª   k820309    9          19.0        ×ÛÁ]                                                                                                          
       /home/albd/research/FHIaims/running_copy/FHIaims/src/calculate_dipolemat.f90 CALCULATE_DIPOLEMAT          @@                                                 
                &                                                    @@                                                 
                &                                                    @@                                                 
                &                                                    @@                                                 
                &                                                    @@                                                 
                &                                                    @@                                                 
                &                                                    @@                                                 
                &                                                    @@                                                 
                &                   &                                                    @@                              	                                   &                   &                                                    @@                              
                                   &                                                    @@                                                                 &                                                    @@                                                                 &                                                                                          
                                                      
                                                      
                  @                                                       @                                            #         @                                                                                                                                                                              #         @                                                                                                                                                                                          #         @                                                                                                                                                                                                          #         @                                                                                                                                                                                                  #         @                                                                                                                                                                              #         @                                                                                                                                                                                                      #         @                                                      #N_POINTS    #PARTITION    #N_COMPUTE    #GRADIENT_BASIS_WAVE    #N_BASIS_LIST    #WAVE    #DIPOLEMAT_SHELL    #I_COORD                                                                                                                                                 D @                                                                                                        
     p          5  p        r        5  p        r                                D @                                                                                                        
         p        p        p        5  p        r    p          5  p        r      p          5  p        r        5  p        r      p          5  p        r                                                                                                                                          
       p        5  p        r    p          5  p        r      5  p        r        5  p        r      5  p        r                               D @                                                  
       p        5  p        r    p          5  p        r      5  p        r        5  p        r      5  p        r                                                                              %         @                                !                   
       #ERGEBNIS%N_GREENENERGY "   #ABTEW #   #ABTEWSEEBECK $   #ZELLGR %   #SWITCH &                                                                                                                                     "                                                     #                    
     p          5 r "       5 r "                                                              $                    
     p          5 r "       5 r "                                                               %     
                                                  &            #         @                                   '                   #CALC_FERMIDERIV%N_GREENENERGY (   #FERMIDERIV )   #CHEMICAL_POTENTIAL *                                                                                                                                                                 (                     D                                )                    
     p          5 r (       5 r (                                                               *     
       #         @                                   +                   #CALC_DIEL%N_K_POINTS ,   #CALC_DIEL%N_SPIN -   #CALC_DIEL%N_STATES .   #CALC_DIEL%N_GREENENERGY /   #CALC_DIEL%N_OMEGA 0   #DIPELEMENTXI 1   #DIPELEMENTXJ 2   #DIE_EL 3   #SEEBECK 4   #ABTEW 5   #ABTEWSEEBECK 6   #FERMIDERIV 7   #OMEGAPL 8   #CHEMICAL_POTENTIAL 9   #KS_EIGEN :   #KS_EIGEN_FULL ;   #K_WEIGHT <                                                                                                                                                           ,                                                        -                                                        .                                                        /                                                        0                                                     1                         p            '  5 r .   n                                           15 r .   n                                          2      '  5 r .   n                                          15 r .   n                                          2                                                              2                         p            '  5 r .   n                                           15 r .   n                                          2      '  5 r .   n                                          15 r .   n                                          2                              D                                3                    
     p          5 r 0       5 r 0                              D                                4                    
     p          5 r 0       5 r 0                              D                                5                    
     p          5 r /       5 r /                              D                                6                    
     p          5 r /       5 r /                                                              7                    
     p          5 r /       5 r /                               D                                8     
                                                 9     
                
                                 :                    
      p        5 r .   p          5 r .     5 r -       5 r .     5 r -                              
                                 ;                    
        p        5 r -   p        5 r .   p          5 r .     5 r -     5 r ,       5 r .     5 r -     5 r ,                                                               <     
       #         @                                   =                	   #OUT_DIE_EL%N_GREENENERGY >   #OUT_DIE_EL%N_OMEGA ?   #DIE_EL @   #SEEBECK A   #ABTEW B   #ABTEWSEEBECK C   #FERMIDERIV D   #EP1_IN E   #EP2_IN F   #ZELLGR G   #CHEMPOT H                                                                                                                                                                >                                                        ?                     
                                 @                    
 !   p          5 r ?       5 r ?                              
                                 A                    
 "   p          5 r ?       5 r ?                              
                                 B                    
 #   p          5 r >       5 r >                              
                                 C                    
 $   p          5 r >       5 r >                              
                                 D                    
 %   p          5 r >       5 r >                               
  @                              E                                     
  @                              F                                                                     G     
                                                 H     
       #         @                                   I                   #CONSTRUCT_DIPOLEMAT_P1%N_CENTERS_BASIS_I J   #CONSTRUCT_DIPOLEMAT_P1%N_BASIS K   #CONSTRUCT_DIPOLEMAT_P1%N_CELLS_IN_HAMILTONIAN L   #DIPOLE_MAT_FULL_ONED M   #DIPOLEMAT_FULL_W N   #DIPOLEMAT_FULL_W_COMPLEX O   #K_PHASE_AC P   #WORK_HAM Q                                                                                                                                                                                             J                                                        K                                                        L                                                     M                    
 &    p            5 r L   5 r K   5 r K         5 r L   5 r K   5 r K                              D                                N                    
 '      p        5 r K   p          5 r K     5 r K       5 r K     5 r K                              D                                O                     (      p        5 r K   p          5 r K     5 r K       5 r K     5 r K                                                              P                     *    p          5 r L       5 r L                                                              Q                    
 )      p        5 r J   p          5 r J     5 r J       5 r J     5 r J                     #         @                                   R                   #CALC_DIPELEMENT_P0%N_SPIN S   #CALC_DIPELEMENT_P0%N_STATES T   #CALC_DIPELEMENT_P0%N_BASIS U   #DIPELEMENT V   #DIPOLEMAT_FULL_W W   #DIPOLEMAT_FULL_W_COMPLEX X   #KS_VEC Y   #KS_VEC_COMPLEX Z   #K_POINT [                                                                                                                                                                             S                                                        T                       @ @                              U                     D @                              V                     0    p            '  5 r T   n                                           15 r T   n                                          2      '  5 r T   n                                          15 r T   n                                          2                              
@ @                              W                    
 .     p        5 r U   p          5 r U     5 r U       5 r U     5 r U                              
@ @                              X                     /     p        5 r U   p          5 r U     5 r U       5 r U     5 r U                              
@ @                              Y                    
 -       p        5 r T   p        5 r U   p          5 r U     5 r T     5 r S       5 r U     5 r T     5 r S                              
                                 Z                     ,       p        5 r T   p        5 r U   p          5 r U     5 r T     5 r S       5 r U     5 r T     5 r S                               
                                  [           #         @                                   \                   #GET_STATE_MINMAX%N_K_POINTS ]   #GET_STATE_MINMAX%N_SPIN ^   #GET_STATE_MINMAX%N_STATES _   #KS_EIGEN `                                                                                                                                                                     ]                                                        ^                                                        _                     
                                 `                    
 3       p        5 r ^   p        5 r _   p          5 r _     5 r ^     5 r ]       5 r _     5 r ^     5 r ]                     #         @                                   a                   #OUT_DIPELEMENT%N_K_POINTS b   #OUT_DIPELEMENT%N_SPIN c   #OUT_DIPELEMENT%N_STATES d   #DIPELEMENT_ONE e   #DIPELEMENT_TWO f   #DIPELEMENT_THREE g   #KS_EIGENVALUE h   #K_POINT i                                                                                                                                                             b                                                        c                                                        d                     
                                 e                     5   p            '  5 r d   n                                           15 r d   n                                          2      '  5 r d   n                                          15 r d   n                                          2                              
                                 f                     6   p            '  5 r d   n                                           15 r d   n                                          2      '  5 r d   n                                          15 r d   n                                          2                              
                                 g                     7   p            '  5 r d   n                                           15 r d   n                                          2      '  5 r d   n                                          15 r d   n                                          2                              
                                 h                    
 8       p        5 r c   p        5 r d   p          5 r d     5 r c     5 r b       5 r d     5 r c     5 r b                               
                                  i           #         @                                   j     
            #CALCULATE_DIPOLEMAT_P0%N_MAX_COMPUTE_HAM k   #CALCULATE_DIPOLEMAT_P0%N_MAX_COMPUTE_FNS_HAM l   #CALCULATE_DIPOLEMAT_P0%N_MAX_COMPUTE_ATOMS m   #CALCULATE_DIPOLEMAT_P0%N_CENTERS n   #CALCULATE_DIPOLEMAT_P0%N_BASIS_FNS o   #CALCULATE_DIPOLEMAT_P0%N_MAX_BATCH_SIZE p   #CALCULATE_DIPOLEMAT_P0%N_CENTERS_INTEGRALS q   #CALCULATE_DIPOLEMAT_P0%N_SPECIES r   #CALCULATE_DIPOLEMAT_P0%N_FULL_POINTS s   #CALCULATE_DIPOLEMAT_P0%N_CENTERS_BASIS_I t   #CALCULATE_DIPOLEMAT_P0%GRID_POINT u   #CALCULATE_DIPOLEMAT_P0%BATCH_OF_POINTS z   #CALCULATE_DIPOLEMAT_P0%BATCH_PERMUTATION    #MPI_TASKS!CALCULATE_DIPOLEMAT_P0%MPIFCMB5    #MPI_TASKS!CALCULATE_DIPOLEMAT_P0%MPIFCMB9    #MPI_TASKS!CALCULATE_DIPOLEMAT_P0%MPIPRIV1    #MPI_TASKS!CALCULATE_DIPOLEMAT_P0%MPIPRIV2    #MPI_TASKS!CALCULATE_DIPOLEMAT_P0%MPIPRIVC    #PARTITION_TAB_STD    #BASIS_L_MAX    #DIPOLE_MAT_FULL_ONED    #DIRECTION    #COUNT1                                                                                                                                                                                                                                        @                          u     '(                    #COORDS v   #INDEX_ATOM w   #INDEX_RADIAL x   #INDEX_ANGULAR y                                              v                              
  p          p            p                                                                      w                                                              x                                                              y                                     @                          z     '                     #SIZE {   #POINTS |   #BATCH_N_COMPUTE }   #BATCH_I_BASIS ~                                               {                                                             |                   (             #CALCULATE_DIPOLEMAT_P0%GRID_POINT u             &                                                                                       }     P                                                       ~            X                             &                                                            @                                '                    #INITIALIZED    #N_MY_BATCHES    #N_FULL_POINTS    #BATCHES    #PERM_BATCH_OWNER    #POINT_SEND_CNT    #POINT_SEND_OFF    #POINT_RECV_CNT    #POINT_RECV_OFF    #N_BASIS_LOCAL    #N_LOCAL_MATRIX_SIZE    #I_BASIS_LOCAL    #I_BASIS_GLB_TO_LOC    #PARTITION_TAB                  $                                                               $                                                              $                                                            $                                                               #CALCULATE_DIPOLEMAT_P0%BATCH_OF_POINTS z             &                                                       $                                          X                             &                                                       $                                                                        &                                                       $                                          è                             &                                                       $                                          0                            &                                                       $                                          x             	               &                                                         $                                   À      
                    $                                   Ä                        $                                          È                            &                                                       $                                                                      &                                                       $                                         X                
            &                                                      @ @                              k                       @ @                              l                       @ @                              m                                                        n                                                        o                                                        p                       @ @                              q                                                        r                                                        s                       @ @                              t                                                                              #CALCULATE_DIPOLEMAT_P0%MPIFCMB5%MPI_UNWEIGHTED                                                                                                                                  #CALCULATE_DIPOLEMAT_P0%MPIFCMB9%MPI_WEIGHTS_EMPTY                                                                                                                                  #CALCULATE_DIPOLEMAT_P0%MPIPRIV1%MPI_BOTTOM    #CALCULATE_DIPOLEMAT_P0%MPIPRIV1%MPI_IN_PLACE    #CALCULATE_DIPOLEMAT_P0%MPIPRIV1%MPI_STATUS_IGNORE                                                                                                                                                                                                   p          p            p                                                                                                    #CALCULATE_DIPOLEMAT_P0%MPIPRIV2%MPI_STATUSES_IGNORE    #CALCULATE_DIPOLEMAT_P0%MPIPRIV2%MPI_ERRCODES_IGNORE                                                                             p          p          p            p          p                                                                                                          p          p            p                                                                                                    #CALCULATE_DIPOLEMAT_P0%MPIPRIVC%MPI_ARGVS_NULL    #CALCULATE_DIPOLEMAT_P0%MPIPRIVC%MPI_ARGV_NULL    -                                                                         p          p          p            p          p                                  -                                                                        p          p            p                                                                                               
 9    p          5 r s       5 r s                              D @                                                    :    p          5 r r       5 r r                            @  D @                                                  
 ;    p          1     1                                                                                                                                                  #         @                                  ¡                   #UPDATE_FULL_MATRIX_P0XXX%N_CELLS ¢   #UPDATE_FULL_MATRIX_P0XXX%N_BASIS £   #UPDATE_FULL_MATRIX_P0XXX%N_CELLS_IN_HAMILTONIAN ¤   #N_COMPUTE_C ¥   #N_COMPUTE_A ¦   #I_BASIS §   #MATRIX_SHELL ¨   #MATRIX ©                                                                                                                                                                                                     ¢                                                        £                                                        ¤                                                       ¥                                                       ¦                                                      §                     e    p          5  p        r ¥       5  p        r ¥                                                              ¨                    
 f      p        5  p        r ¥   p          5  p        r ¥     5  p        r ¦       5  p        r ¥     5  p        r ¦                              D                                ©                    
 g    p            5 r ¤   5 r £   5 r £         5 r ¤   5 r £   5 r £                            i      fn#fn %   	         DIPOLE_MAT_FULL_ONED )            DIPOLE_MAT_FULL_ONED_TWO    !         DIE_EL    ­         SEEBECK    9         ABTEW    Å         ABTEWSEEBECK    Q         FERMIDERIV !   Ý  ¤       DIPOLEMAT_FULL_W )     ¤       DIPOLEMAT_FULL_W_COMPLEX    %         DIPELEMENT_ONE    ±         DIPELEMENT_TWO !   =         DIPELEMENT_THREE    É  @       OMEGAPL #   	  @       ABTEWINTEGRAL_COND #   I  @       ABTEWINTEGRAL_SEEB      @       N_STATE_MIN    É  @       N_STATE_MAX !   		  ¾       ALLOCATE_SPECTRA $   Ç	  Ê       ALLOCATE_DIPOLE_MAT (   
  Ú       ALLOCATE_DIPOLE_MAT_TWO &   k  Ò       ALLOCATE_DIPOLE_MAT_K !   =  ¾       CLEAN_DIPOLE_MAT '   û  Ö       CLEAN_DIPOLE_MAT_FINAL $   Ñ  M      EVALUATE_DIPOLE_MAT -     @   a   EVALUATE_DIPOLE_MAT%N_POINTS .   ^  ´   a   EVALUATE_DIPOLE_MAT%PARTITION .     @   a   EVALUATE_DIPOLE_MAT%N_COMPUTE 8   R  d  a   EVALUATE_DIPOLE_MAT%GRADIENT_BASIS_WAVE 1   ¶  @   a   EVALUATE_DIPOLE_MAT%N_BASIS_LIST )   ö  $  a   EVALUATE_DIPOLE_MAT%WAVE 4     $  a   EVALUATE_DIPOLE_MAT%DIPOLEMAT_SHELL ,   >  @   a   EVALUATE_DIPOLE_MAT%I_COORD    ~  ÷       ERGEBNIS 7   u  @     ERGEBNIS%N_GREENENERGY+RUNTIME_CHOICES    µ     a   ERGEBNIS%ABTEW &   I     a   ERGEBNIS%ABTEWSEEBECK     Ý  @   a   ERGEBNIS%ZELLGR       @   a   ERGEBNIS%SWITCH     ]        CALC_FERMIDERIV >   b  @     CALC_FERMIDERIV%N_GREENENERGY+RUNTIME_CHOICES +   ¢     a   CALC_FERMIDERIV%FERMIDERIV 3   6  @   a   CALC_FERMIDERIV%CHEMICAL_POTENTIAL    v  î      CALC_DIEL 0   d  @     CALC_DIEL%N_K_POINTS+DIMENSIONS ,   ¤  @     CALC_DIEL%N_SPIN+DIMENSIONS .   ä  @     CALC_DIEL%N_STATES+DIMENSIONS 8   $  @     CALC_DIEL%N_GREENENERGY+RUNTIME_CHOICES 2   d  @     CALC_DIEL%N_OMEGA+RUNTIME_CHOICES '   ¤     a   CALC_DIEL%DIPELEMENTXI '   D     a   CALC_DIEL%DIPELEMENTXJ !   ä     a   CALC_DIEL%DIE_EL "   x      a   CALC_DIEL%SEEBECK     !     a   CALC_DIEL%ABTEW '    !     a   CALC_DIEL%ABTEWSEEBECK %   4"     a   CALC_DIEL%FERMIDERIV "   È"  @   a   CALC_DIEL%OMEGAPL -   #  @   a   CALC_DIEL%CHEMICAL_POTENTIAL #   H#  Ô   a   CALC_DIEL%KS_EIGEN (   $    a   CALC_DIEL%KS_EIGEN_FULL #   0%  @   a   CALC_DIEL%K_WEIGHT    p%  f      OUT_DIE_EL 9   Ö&  @     OUT_DIE_EL%N_GREENENERGY+RUNTIME_CHOICES 3   '  @     OUT_DIE_EL%N_OMEGA+RUNTIME_CHOICES "   V'     a   OUT_DIE_EL%DIE_EL #   ê'     a   OUT_DIE_EL%SEEBECK !   ~(     a   OUT_DIE_EL%ABTEW (   )     a   OUT_DIE_EL%ABTEWSEEBECK &   ¦)     a   OUT_DIE_EL%FERMIDERIV "   :*  P   a   OUT_DIE_EL%EP1_IN "   *  P   a   OUT_DIE_EL%EP2_IN "   Ú*  @   a   OUT_DIE_EL%ZELLGR #   +  @   a   OUT_DIE_EL%CHEMPOT '   Z+  Ç      CONSTRUCT_DIPOLEMAT_P1 D   !-  @     CONSTRUCT_DIPOLEMAT_P1%N_CENTERS_BASIS_I+DIMENSIONS :   a-  @     CONSTRUCT_DIPOLEMAT_P1%N_BASIS+DIMENSIONS H   ¡-  @     CONSTRUCT_DIPOLEMAT_P1%N_CELLS_IN_HAMILTONIAN+PBC_LISTS <   á-  Ô   a   CONSTRUCT_DIPOLEMAT_P1%DIPOLE_MAT_FULL_ONED 8   µ.  Ô   a   CONSTRUCT_DIPOLEMAT_P1%DIPOLEMAT_FULL_W @   /  Ô   a   CONSTRUCT_DIPOLEMAT_P1%DIPOLEMAT_FULL_W_COMPLEX 2   ]0     a   CONSTRUCT_DIPOLEMAT_P1%K_PHASE_AC 0   ñ0  Ô   a   CONSTRUCT_DIPOLEMAT_P1%WORK_HAM #   Å1        CALC_DIPELEMENT_P0 5   \3  @     CALC_DIPELEMENT_P0%N_SPIN+DIMENSIONS 7   3  @     CALC_DIPELEMENT_P0%N_STATES+DIMENSIONS 6   Ü3  @     CALC_DIPELEMENT_P0%N_BASIS+DIMENSIONS .   4     a   CALC_DIPELEMENT_P0%DIPELEMENT 4   ¼5  Ô   a   CALC_DIPELEMENT_P0%DIPOLEMAT_FULL_W <   6  Ô   a   CALC_DIPELEMENT_P0%DIPOLEMAT_FULL_W_COMPLEX *   d7    a   CALC_DIPELEMENT_P0%KS_VEC 2   x8    a   CALC_DIPELEMENT_P0%KS_VEC_COMPLEX +   9  @   a   CALC_DIPELEMENT_P0%K_POINT !   Ì9  )      GET_STATE_MINMAX 7   õ:  @     GET_STATE_MINMAX%N_K_POINTS+DIMENSIONS 3   5;  @     GET_STATE_MINMAX%N_SPIN+DIMENSIONS 5   u;  @     GET_STATE_MINMAX%N_STATES+DIMENSIONS *   µ;    a   GET_STATE_MINMAX%KS_EIGEN    É<  k      OUT_DIPELEMENT 5   4>  @     OUT_DIPELEMENT%N_K_POINTS+DIMENSIONS 1   t>  @     OUT_DIPELEMENT%N_SPIN+DIMENSIONS 3   ´>  @     OUT_DIPELEMENT%N_STATES+DIMENSIONS .   ô>     a   OUT_DIPELEMENT%DIPELEMENT_ONE .   @     a   OUT_DIPELEMENT%DIPELEMENT_TWO 0   4B     a   OUT_DIPELEMENT%DIPELEMENT_THREE -   ÔC    a   OUT_DIPELEMENT%KS_EIGENVALUE '   èD  @   a   OUT_DIPELEMENT%K_POINT '   (E        CALCULATE_DIPOLEMAT_P0 8   ÃI        CALCULATE_DIPOLEMAT_P0%GRID_POINT+GRIDS ?   TJ     a   CALCULATE_DIPOLEMAT_P0%GRID_POINT%COORDS+GRIDS C   ðJ  H   a   CALCULATE_DIPOLEMAT_P0%GRID_POINT%INDEX_ATOM+GRIDS E   8K  H   a   CALCULATE_DIPOLEMAT_P0%GRID_POINT%INDEX_RADIAL+GRIDS F   K  H   a   CALCULATE_DIPOLEMAT_P0%GRID_POINT%INDEX_ANGULAR+GRIDS =   ÈK        CALCULATE_DIPOLEMAT_P0%BATCH_OF_POINTS+GRIDS B   VL  H   a   CALCULATE_DIPOLEMAT_P0%BATCH_OF_POINTS%SIZE+GRIDS D   L  »   a   CALCULATE_DIPOLEMAT_P0%BATCH_OF_POINTS%POINTS+GRIDS M   YM  H   a   CALCULATE_DIPOLEMAT_P0%BATCH_OF_POINTS%BATCH_N_COMPUTE+GRIDS K   ¡M     a   CALCULATE_DIPOLEMAT_P0%BATCH_OF_POINTS%BATCH_I_BASIS+GRIDS H   5N  c     CALCULATE_DIPOLEMAT_P0%BATCH_PERMUTATION+LOAD_BALANCING T   O  H   a   CALCULATE_DIPOLEMAT_P0%BATCH_PERMUTATION%INITIALIZED+LOAD_BALANCING U   àO  H   a   CALCULATE_DIPOLEMAT_P0%BATCH_PERMUTATION%N_MY_BATCHES+LOAD_BALANCING V   (P  H   a   CALCULATE_DIPOLEMAT_P0%BATCH_PERMUTATION%N_FULL_POINTS+LOAD_BALANCING P   pP  À   a   CALCULATE_DIPOLEMAT_P0%BATCH_PERMUTATION%BATCHES+LOAD_BALANCING Y   0Q     a   CALCULATE_DIPOLEMAT_P0%BATCH_PERMUTATION%PERM_BATCH_OWNER+LOAD_BALANCING W   ÄQ     a   CALCULATE_DIPOLEMAT_P0%BATCH_PERMUTATION%POINT_SEND_CNT+LOAD_BALANCING W   XR     a   CALCULATE_DIPOLEMAT_P0%BATCH_PERMUTATION%POINT_SEND_OFF+LOAD_BALANCING W   ìR     a   CALCULATE_DIPOLEMAT_P0%BATCH_PERMUTATION%POINT_RECV_CNT+LOAD_BALANCING W   S     a   CALCULATE_DIPOLEMAT_P0%BATCH_PERMUTATION%POINT_RECV_OFF+LOAD_BALANCING V   T  H   a   CALCULATE_DIPOLEMAT_P0%BATCH_PERMUTATION%N_BASIS_LOCAL+LOAD_BALANCING \   \T  H   a   CALCULATE_DIPOLEMAT_P0%BATCH_PERMUTATION%N_LOCAL_MATRIX_SIZE+LOAD_BALANCING V   ¤T     a   CALCULATE_DIPOLEMAT_P0%BATCH_PERMUTATION%I_BASIS_LOCAL+LOAD_BALANCING [   8U     a   CALCULATE_DIPOLEMAT_P0%BATCH_PERMUTATION%I_BASIS_GLB_TO_LOC+LOAD_BALANCING V   ÌU     a   CALCULATE_DIPOLEMAT_P0%BATCH_PERMUTATION%PARTITION_TAB+LOAD_BALANCING D   `V  @     CALCULATE_DIPOLEMAT_P0%N_MAX_COMPUTE_HAM+DIMENSIONS H    V  @     CALCULATE_DIPOLEMAT_P0%N_MAX_COMPUTE_FNS_HAM+DIMENSIONS F   àV  @     CALCULATE_DIPOLEMAT_P0%N_MAX_COMPUTE_ATOMS+DIMENSIONS <    W  @     CALCULATE_DIPOLEMAT_P0%N_CENTERS+DIMENSIONS >   `W  @     CALCULATE_DIPOLEMAT_P0%N_BASIS_FNS+DIMENSIONS C    W  @     CALCULATE_DIPOLEMAT_P0%N_MAX_BATCH_SIZE+DIMENSIONS F   àW  @     CALCULATE_DIPOLEMAT_P0%N_CENTERS_INTEGRALS+DIMENSIONS <    X  @     CALCULATE_DIPOLEMAT_P0%N_SPECIES+DIMENSIONS @   `X  @     CALCULATE_DIPOLEMAT_P0%N_FULL_POINTS+DIMENSIONS D    X  @     CALCULATE_DIPOLEMAT_P0%N_CENTERS_BASIS_I+DIMENSIONS M   àX       MPI_TASKS!CALCULATE_DIPOLEMAT_P0%MPIFCMB5+MPI_TASKS=MPIFCMB5 I   dY  H     CALCULATE_DIPOLEMAT_P0%MPIFCMB5%MPI_UNWEIGHTED+MPI_TASKS M   ¬Y       MPI_TASKS!CALCULATE_DIPOLEMAT_P0%MPIFCMB9+MPI_TASKS=MPIFCMB9 L   3Z  H     CALCULATE_DIPOLEMAT_P0%MPIFCMB9%MPI_WEIGHTS_EMPTY+MPI_TASKS M   {Z  é     MPI_TASKS!CALCULATE_DIPOLEMAT_P0%MPIPRIV1+MPI_TASKS=MPIPRIV1 E   d[  H     CALCULATE_DIPOLEMAT_P0%MPIPRIV1%MPI_BOTTOM+MPI_TASKS G   ¬[  H     CALCULATE_DIPOLEMAT_P0%MPIPRIV1%MPI_IN_PLACE+MPI_TASKS L   ô[  ¤     CALCULATE_DIPOLEMAT_P0%MPIPRIV1%MPI_STATUS_IGNORE+MPI_TASKS M   \  Â     MPI_TASKS!CALCULATE_DIPOLEMAT_P0%MPIPRIV2+MPI_TASKS=MPIPRIV2 N   Z]  Ä     CALCULATE_DIPOLEMAT_P0%MPIPRIV2%MPI_STATUSES_IGNORE+MPI_TASKS N   ^  ¤     CALCULATE_DIPOLEMAT_P0%MPIPRIV2%MPI_ERRCODES_IGNORE+MPI_TASKS M   Â^  ·     MPI_TASKS!CALCULATE_DIPOLEMAT_P0%MPIPRIVC+MPI_TASKS=MPIPRIVC I   y_  Ä     CALCULATE_DIPOLEMAT_P0%MPIPRIVC%MPI_ARGVS_NULL+MPI_TASKS H   =`  ¤     CALCULATE_DIPOLEMAT_P0%MPIPRIVC%MPI_ARGV_NULL+MPI_TASKS 9   á`     a   CALCULATE_DIPOLEMAT_P0%PARTITION_TAB_STD 3   ua     a   CALCULATE_DIPOLEMAT_P0%BASIS_L_MAX <   	b     a   CALCULATE_DIPOLEMAT_P0%DIPOLE_MAT_FULL_ONED 1   b  P   a   CALCULATE_DIPOLEMAT_P0%DIRECTION .   Ýb  @   a   CALCULATE_DIPOLEMAT_P0%COUNT1 )   c  ¬      UPDATE_FULL_MATRIX_P0XXX ;   Éd  @     UPDATE_FULL_MATRIX_P0XXX%N_CELLS+PBC_LISTS <   	e  @     UPDATE_FULL_MATRIX_P0XXX%N_BASIS+DIMENSIONS J   Ie  @     UPDATE_FULL_MATRIX_P0XXX%N_CELLS_IN_HAMILTONIAN+PBC_LISTS 5   e  @   a   UPDATE_FULL_MATRIX_P0XXX%N_COMPUTE_C 5   Ée  @   a   UPDATE_FULL_MATRIX_P0XXX%N_COMPUTE_A 1   	f  ´   a   UPDATE_FULL_MATRIX_P0XXX%I_BASIS 6   ½f  $  a   UPDATE_FULL_MATRIX_P0XXX%MATRIX_SHELL 0   ág  Ô   a   UPDATE_FULL_MATRIX_P0XXX%MATRIX 