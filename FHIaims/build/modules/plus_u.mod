  �@  p   k820309    9          19.0        h��]                                                                                                          
       /home/albd/research/FHIaims/running_copy/FHIaims/src/plus_u.f90 PLUS_U              ALLOCATE_PLUS_U PLUS_U_INIT_IDX DEALLOCATE_PLUS_U OCC_NUMBERS_PLUS_U ADD_PLUS_U_TO_HAMILTONIAN PLUS_U_ENERGY_CORRECTION_TERM PLUS_U_WRITE_OCC_MAT PLUS_U_READ_OCC_MAT PLUS_U_MATRIX_ERROR PLUS_U_CHECK_OCCUPATION PLUS_U_EIGENVALUES OCC_NUMBERS_PLUS_U_MULLIKEN ADD_PLUS_U_MULLIKEN_TO_HAMILTONIAN OCC_NUMBERS_PLUS_U_FULL ADD_PLUS_U_FULL_TO_HAMILTONIAN PLUS_U_RAMPING_WORK PLUS_U_PETUKHOV_MIXING_DEFINED PLUS_U_RAMPING_DEFINED PLUS_U_MATRIX_RELEASE_DEFINED PLUS_U_EIGENVALUES_DEFINED PLUS_U_OCCUPATION_MATRIX_CONTROL_READ PLUS_U_OCCUPATION_MATRIX_CONTROL_WRITE PLUS_U_HYDROS_DEFINED PLUS_U_MATRIX_ERROR_DEFINED OCC_MAT_FILE_EXISTS PLUS_U_MATRIX_RELEASE PLUS_U_RAMPING_ACCURACY PLUS_U_PETUKHOV_MIXING PLUS_U_ENERGY_CORRECTION                                                     
                            @                              
                                                                           @                                                                                                                                                                                                                                                                                                                              	                                                         
                                                                                                                      
                                                        
                                                        
                  @                                     
       #         @                                                                                     #         @                                                                                                                                                #         @                                                        #         @                                                     #OCC_NUMBERS_PLUS_U%N_HAMILTONIAN_MATRIX_SIZE    #MPI_TASKS!OCC_NUMBERS_PLUS_U%MPIFCMB5    #MPI_TASKS!OCC_NUMBERS_PLUS_U%MPIFCMB9    #MPI_TASKS!OCC_NUMBERS_PLUS_U%MPIPRIV1    #MPI_TASKS!OCC_NUMBERS_PLUS_U%MPIPRIV2    #MPI_TASKS!OCC_NUMBERS_PLUS_U%MPIPRIVC     #DENSITY_MATRIX_SPARSE #   #DENSITY_MATRIX $   #I_SPIN %                                                                                                                                                                                                                                                              #OCC_NUMBERS_PLUS_U%MPIFCMB5%MPI_UNWEIGHTED              �            �                                                                                                        #OCC_NUMBERS_PLUS_U%MPIFCMB9%MPI_WEIGHTS_EMPTY              �            �                                                                                                        #OCC_NUMBERS_PLUS_U%MPIPRIV1%MPI_BOTTOM    #OCC_NUMBERS_PLUS_U%MPIPRIV1%MPI_IN_PLACE    #OCC_NUMBERS_PLUS_U%MPIPRIV1%MPI_STATUS_IGNORE              �            �                                                �            �                                               �            �                                                  p          p            p                                                                                                    #OCC_NUMBERS_PLUS_U%MPIPRIV2%MPI_STATUSES_IGNORE    #OCC_NUMBERS_PLUS_U%MPIPRIV2%MPI_ERRCODES_IGNORE              �            �                                                   p          p          p            p          p                                            �            �                                                  p          p            p                                                                                                     #OCC_NUMBERS_PLUS_U%MPIPRIVC%MPI_ARGVS_NULL !   #OCC_NUMBERS_PLUS_U%MPIPRIVC%MPI_ARGV_NULL "   -          �            �                  !                                 p          p          p            p          p                                  -          �            �                  "                                p          p            p                                           
                                 #                    
    p          5 r        5 r                                
                                $                   
               &                   &                                                     
                                  %           #         @                                   &                   #ADD_PLUS_U_TO_HAMILTONIAN%N_HAMILTONIAN_MATRIX_SIZE '   #ADD_PLUS_U_TO_HAMILTONIAN%N_SPIN (   #HAMILTONIAN )                                                                                                                                                                                                            '                                                        (                     
D                                )                    
       p        5 r '   p          5 r '     5 r (       5 r '     5 r (                     #         @                                   *                     #         @                                   +                     #         @                                   ,                    #MPI_TASKS!PLUS_U_READ_OCC_MAT%MPIFCMB5 -   #MPI_TASKS!PLUS_U_READ_OCC_MAT%MPIFCMB9 /   #MPI_TASKS!PLUS_U_READ_OCC_MAT%MPIPRIV1 1   #MPI_TASKS!PLUS_U_READ_OCC_MAT%MPIPRIV2 5   #MPI_TASKS!PLUS_U_READ_OCC_MAT%MPIPRIVC 8                                                                                                                                               -                          #PLUS_U_READ_OCC_MAT%MPIFCMB5%MPI_UNWEIGHTED .             �            �                  .                                                            /                          #PLUS_U_READ_OCC_MAT%MPIFCMB9%MPI_WEIGHTS_EMPTY 0             �            �                  0                                                            1                          #PLUS_U_READ_OCC_MAT%MPIPRIV1%MPI_BOTTOM 2   #PLUS_U_READ_OCC_MAT%MPIPRIV1%MPI_IN_PLACE 3   #PLUS_U_READ_OCC_MAT%MPIPRIV1%MPI_STATUS_IGNORE 4             �            �                  2                              �            �                  3                             �            �                  4                                p          p            p                                                                          5                          #PLUS_U_READ_OCC_MAT%MPIPRIV2%MPI_STATUSES_IGNORE 6   #PLUS_U_READ_OCC_MAT%MPIPRIV2%MPI_ERRCODES_IGNORE 7             �            �                  6                                 p          p          p            p          p                                            �            �                  7                                p          p            p                                                                          8                          #PLUS_U_READ_OCC_MAT%MPIPRIVC%MPI_ARGVS_NULL 9   #PLUS_U_READ_OCC_MAT%MPIPRIVC%MPI_ARGV_NULL :   -          �            �                  9                                 p          p          p            p          p                                  -          �            �                  :                                p          p            p                                  #         @                                   ;                     #         @                                   <                     #         @                                   =                     #         @                                   >                  #OCC_NUMBERS_PLUS_U_MULLIKEN%N_HAMILTONIAN_MATRIX_SIZE ?   #MPI_TASKS!OCC_NUMBERS_PLUS_U_MULLIKEN%MPIFCMB5 @   #MPI_TASKS!OCC_NUMBERS_PLUS_U_MULLIKEN%MPIFCMB9 B   #MPI_TASKS!OCC_NUMBERS_PLUS_U_MULLIKEN%MPIPRIV1 D   #MPI_TASKS!OCC_NUMBERS_PLUS_U_MULLIKEN%MPIPRIV2 H   #MPI_TASKS!OCC_NUMBERS_PLUS_U_MULLIKEN%MPIPRIVC K   #DENSITY_MATRIX_SPARSE N   #DENSITY_MATRIX O   #I_SPIN P                                                                                                                                                                                                                    ?                                                    @                          #OCC_NUMBERS_PLUS_U_MULLIKEN%MPIFCMB5%MPI_UNWEIGHTED A             �            �                  A                                                            B                          #OCC_NUMBERS_PLUS_U_MULLIKEN%MPIFCMB9%MPI_WEIGHTS_EMPTY C             �            �                  C                                                            D                          #OCC_NUMBERS_PLUS_U_MULLIKEN%MPIPRIV1%MPI_BOTTOM E   #OCC_NUMBERS_PLUS_U_MULLIKEN%MPIPRIV1%MPI_IN_PLACE F   #OCC_NUMBERS_PLUS_U_MULLIKEN%MPIPRIV1%MPI_STATUS_IGNORE G             �            �                  E                              �            �                  F                             �            �                  G                                p          p            p                                                                          H                          #OCC_NUMBERS_PLUS_U_MULLIKEN%MPIPRIV2%MPI_STATUSES_IGNORE I   #OCC_NUMBERS_PLUS_U_MULLIKEN%MPIPRIV2%MPI_ERRCODES_IGNORE J             �            �                  I                                 p          p          p            p          p                                            �            �                  J                                p          p            p                                                                          K                          #OCC_NUMBERS_PLUS_U_MULLIKEN%MPIPRIVC%MPI_ARGVS_NULL L   #OCC_NUMBERS_PLUS_U_MULLIKEN%MPIPRIVC%MPI_ARGV_NULL M   -          �            �                  L                                 p          p          p            p          p                                  -          �            �                  M                                p          p            p                                           
                                 N                    
    p          5 r ?       5 r ?                               
                                O                   
               &                   &                                                     
                                  P           #         @                                   Q                   #ADD_PLUS_U_MULLIKEN_TO_HAMILTONIAN%N_HAMILTONIAN_MATRIX_SIZE R   #ADD_PLUS_U_MULLIKEN_TO_HAMILTONIAN%N_SPIN S   #HAMILTONIAN T                                                                                                                                                                                                                                                R                                                        S                     
D                                T                    
       p        5 r R   p          5 r R     5 r S       5 r R     5 r S                     #         @                                   U                  #OCC_NUMBERS_PLUS_U_FULL%N_HAMILTONIAN_MATRIX_SIZE V   #MPI_TASKS!OCC_NUMBERS_PLUS_U_FULL%MPIFCMB5 W   #MPI_TASKS!OCC_NUMBERS_PLUS_U_FULL%MPIFCMB9 Y   #MPI_TASKS!OCC_NUMBERS_PLUS_U_FULL%MPIPRIV1 [   #MPI_TASKS!OCC_NUMBERS_PLUS_U_FULL%MPIPRIV2 _   #MPI_TASKS!OCC_NUMBERS_PLUS_U_FULL%MPIPRIVC b   #DENSITY_MATRIX_SPARSE e   #DENSITY_MATRIX f   #I_SPIN g                                                                                                                                                                                                    V                                                    W                          #OCC_NUMBERS_PLUS_U_FULL%MPIFCMB5%MPI_UNWEIGHTED X             �            �                  X                                                            Y                          #OCC_NUMBERS_PLUS_U_FULL%MPIFCMB9%MPI_WEIGHTS_EMPTY Z             �            �                  Z                                                            [                          #OCC_NUMBERS_PLUS_U_FULL%MPIPRIV1%MPI_BOTTOM \   #OCC_NUMBERS_PLUS_U_FULL%MPIPRIV1%MPI_IN_PLACE ]   #OCC_NUMBERS_PLUS_U_FULL%MPIPRIV1%MPI_STATUS_IGNORE ^             �            �                  \                              �            �                  ]                             �            �                  ^                                p          p            p                                                                          _                          #OCC_NUMBERS_PLUS_U_FULL%MPIPRIV2%MPI_STATUSES_IGNORE `   #OCC_NUMBERS_PLUS_U_FULL%MPIPRIV2%MPI_ERRCODES_IGNORE a             �            �                  `                                 p          p          p            p          p                                            �            �                  a                                p          p            p                                                                          b                          #OCC_NUMBERS_PLUS_U_FULL%MPIPRIVC%MPI_ARGVS_NULL c   #OCC_NUMBERS_PLUS_U_FULL%MPIPRIVC%MPI_ARGV_NULL d   -          �            �                  c                                 p          p          p            p          p                                  -          �            �                  d                                p          p            p                                           
                                 e                    
    p          5 r V       5 r V                               
                                f                   
               &                   &                                                     
                                  g           #         @                                   h                   #ADD_PLUS_U_FULL_TO_HAMILTONIAN%N_HAMILTONIAN_MATRIX_SIZE i   #ADD_PLUS_U_FULL_TO_HAMILTONIAN%N_SPIN j   #HAMILTONIAN k                                                                                                                                                                                                                                i                                                        j                     
D                                k                    
       p        5 r i   p          5 r i     5 r j       5 r i     5 r j                     #         @                                   l                    #ETOT m   #ETOT_PREV n                                              
                                 m     
                
                                 n     
         �   O      fn#fn    �   �  b   uapp(PLUS_U    �  @   J  PHYSICS      @   J  PBC_LISTS /   N  @       PLUS_U_PETUKHOV_MIXING_DEFINED '   �  @       PLUS_U_RAMPING_DEFINED .   �  @       PLUS_U_MATRIX_RELEASE_DEFINED +     @       PLUS_U_EIGENVALUES_DEFINED 6   N  @       PLUS_U_OCCUPATION_MATRIX_CONTROL_READ 7   �  @       PLUS_U_OCCUPATION_MATRIX_CONTROL_WRITE &   �  @       PLUS_U_HYDROS_DEFINED ,     @       PLUS_U_MATRIX_ERROR_DEFINED $   N  @       OCC_MAT_FILE_EXISTS &   �  @       PLUS_U_MATRIX_RELEASE (   �  @       PLUS_U_RAMPING_ACCURACY '     @       PLUS_U_PETUKHOV_MIXING )   N  @       PLUS_U_ENERGY_CORRECTION     �  e       ALLOCATE_PLUS_U     �  �       PLUS_U_INIT_IDX "   �  H       DEALLOCATE_PLUS_U #   �        OCC_NUMBERS_PLUS_U H   �
  @     OCC_NUMBERS_PLUS_U%N_HAMILTONIAN_MATRIX_SIZE+DIMENSIONS I   (  �   �  MPI_TASKS!OCC_NUMBERS_PLUS_U%MPIFCMB5+MPI_TASKS=MPIFCMB5 E   �  H     OCC_NUMBERS_PLUS_U%MPIFCMB5%MPI_UNWEIGHTED+MPI_TASKS I   �  �   �  MPI_TASKS!OCC_NUMBERS_PLUS_U%MPIFCMB9+MPI_TASKS=MPIFCMB9 H   s  H     OCC_NUMBERS_PLUS_U%MPIFCMB9%MPI_WEIGHTS_EMPTY+MPI_TASKS I   �  �   �  MPI_TASKS!OCC_NUMBERS_PLUS_U%MPIPRIV1+MPI_TASKS=MPIPRIV1 A   �  H     OCC_NUMBERS_PLUS_U%MPIPRIV1%MPI_BOTTOM+MPI_TASKS C   �  H     OCC_NUMBERS_PLUS_U%MPIPRIV1%MPI_IN_PLACE+MPI_TASKS H   (  �     OCC_NUMBERS_PLUS_U%MPIPRIV1%MPI_STATUS_IGNORE+MPI_TASKS I   �  �   �  MPI_TASKS!OCC_NUMBERS_PLUS_U%MPIPRIV2+MPI_TASKS=MPIPRIV2 J   �  �     OCC_NUMBERS_PLUS_U%MPIPRIV2%MPI_STATUSES_IGNORE+MPI_TASKS J   J  �     OCC_NUMBERS_PLUS_U%MPIPRIV2%MPI_ERRCODES_IGNORE+MPI_TASKS I   �  �   �  MPI_TASKS!OCC_NUMBERS_PLUS_U%MPIPRIVC+MPI_TASKS=MPIPRIVC E   �  �     OCC_NUMBERS_PLUS_U%MPIPRIVC%MPI_ARGVS_NULL+MPI_TASKS D   a  �     OCC_NUMBERS_PLUS_U%MPIPRIVC%MPI_ARGV_NULL+MPI_TASKS 9     �   a   OCC_NUMBERS_PLUS_U%DENSITY_MATRIX_SPARSE 2   �  �   a   OCC_NUMBERS_PLUS_U%DENSITY_MATRIX *   =  @   a   OCC_NUMBERS_PLUS_U%I_SPIN *   }  U      ADD_PLUS_U_TO_HAMILTONIAN O   �  @     ADD_PLUS_U_TO_HAMILTONIAN%N_HAMILTONIAN_MATRIX_SIZE+DIMENSIONS <     @     ADD_PLUS_U_TO_HAMILTONIAN%N_SPIN+DIMENSIONS 6   R  �   a   ADD_PLUS_U_TO_HAMILTONIAN%HAMILTONIAN .   &  H       PLUS_U_ENERGY_CORRECTION_TERM %   n  H       PLUS_U_WRITE_OCC_MAT $   �  �      PLUS_U_READ_OCC_MAT J   >  �   �  MPI_TASKS!PLUS_U_READ_OCC_MAT%MPIFCMB5+MPI_TASKS=MPIFCMB5 F   �  H     PLUS_U_READ_OCC_MAT%MPIFCMB5%MPI_UNWEIGHTED+MPI_TASKS J     �   �  MPI_TASKS!PLUS_U_READ_OCC_MAT%MPIFCMB9+MPI_TASKS=MPIFCMB9 I   �  H     PLUS_U_READ_OCC_MAT%MPIFCMB9%MPI_WEIGHTS_EMPTY+MPI_TASKS J   �  �   �  MPI_TASKS!PLUS_U_READ_OCC_MAT%MPIPRIV1+MPI_TASKS=MPIPRIV1 B   �  H     PLUS_U_READ_OCC_MAT%MPIPRIV1%MPI_BOTTOM+MPI_TASKS D   �  H     PLUS_U_READ_OCC_MAT%MPIPRIV1%MPI_IN_PLACE+MPI_TASKS I   C  �     PLUS_U_READ_OCC_MAT%MPIPRIV1%MPI_STATUS_IGNORE+MPI_TASKS J   �  �   �  MPI_TASKS!PLUS_U_READ_OCC_MAT%MPIPRIV2+MPI_TASKS=MPIPRIV2 K   �  �     PLUS_U_READ_OCC_MAT%MPIPRIV2%MPI_STATUSES_IGNORE+MPI_TASKS K   g  �     PLUS_U_READ_OCC_MAT%MPIPRIV2%MPI_ERRCODES_IGNORE+MPI_TASKS J     �   �  MPI_TASKS!PLUS_U_READ_OCC_MAT%MPIPRIVC+MPI_TASKS=MPIPRIVC F   �  �     PLUS_U_READ_OCC_MAT%MPIPRIVC%MPI_ARGVS_NULL+MPI_TASKS E   �   �     PLUS_U_READ_OCC_MAT%MPIPRIVC%MPI_ARGV_NULL+MPI_TASKS $   $!  H       PLUS_U_MATRIX_ERROR (   l!  H       PLUS_U_CHECK_OCCUPATION #   �!  H       PLUS_U_EIGENVALUES ,   �!  g      OCC_NUMBERS_PLUS_U_MULLIKEN Q   c$  @     OCC_NUMBERS_PLUS_U_MULLIKEN%N_HAMILTONIAN_MATRIX_SIZE+DIMENSIONS R   �$  �   �  MPI_TASKS!OCC_NUMBERS_PLUS_U_MULLIKEN%MPIFCMB5+MPI_TASKS=MPIFCMB5 N   ,%  H     OCC_NUMBERS_PLUS_U_MULLIKEN%MPIFCMB5%MPI_UNWEIGHTED+MPI_TASKS R   t%  �   �  MPI_TASKS!OCC_NUMBERS_PLUS_U_MULLIKEN%MPIFCMB9+MPI_TASKS=MPIFCMB9 Q    &  H     OCC_NUMBERS_PLUS_U_MULLIKEN%MPIFCMB9%MPI_WEIGHTS_EMPTY+MPI_TASKS R   H&  �   �  MPI_TASKS!OCC_NUMBERS_PLUS_U_MULLIKEN%MPIPRIV1+MPI_TASKS=MPIPRIV1 J   @'  H     OCC_NUMBERS_PLUS_U_MULLIKEN%MPIPRIV1%MPI_BOTTOM+MPI_TASKS L   �'  H     OCC_NUMBERS_PLUS_U_MULLIKEN%MPIPRIV1%MPI_IN_PLACE+MPI_TASKS Q   �'  �     OCC_NUMBERS_PLUS_U_MULLIKEN%MPIPRIV1%MPI_STATUS_IGNORE+MPI_TASKS R   t(  �   �  MPI_TASKS!OCC_NUMBERS_PLUS_U_MULLIKEN%MPIPRIV2+MPI_TASKS=MPIPRIV2 S   @)  �     OCC_NUMBERS_PLUS_U_MULLIKEN%MPIPRIV2%MPI_STATUSES_IGNORE+MPI_TASKS S   *  �     OCC_NUMBERS_PLUS_U_MULLIKEN%MPIPRIV2%MPI_ERRCODES_IGNORE+MPI_TASKS R   �*  �   �  MPI_TASKS!OCC_NUMBERS_PLUS_U_MULLIKEN%MPIPRIVC+MPI_TASKS=MPIPRIVC N   i+  �     OCC_NUMBERS_PLUS_U_MULLIKEN%MPIPRIVC%MPI_ARGVS_NULL+MPI_TASKS M   -,  �     OCC_NUMBERS_PLUS_U_MULLIKEN%MPIPRIVC%MPI_ARGV_NULL+MPI_TASKS B   �,  �   a   OCC_NUMBERS_PLUS_U_MULLIKEN%DENSITY_MATRIX_SPARSE ;   e-  �   a   OCC_NUMBERS_PLUS_U_MULLIKEN%DENSITY_MATRIX 3   	.  @   a   OCC_NUMBERS_PLUS_U_MULLIKEN%I_SPIN 3   I.  �      ADD_PLUS_U_MULLIKEN_TO_HAMILTONIAN X   �/  @     ADD_PLUS_U_MULLIKEN_TO_HAMILTONIAN%N_HAMILTONIAN_MATRIX_SIZE+DIMENSIONS E   0  @     ADD_PLUS_U_MULLIKEN_TO_HAMILTONIAN%N_SPIN+DIMENSIONS ?   T0  �   a   ADD_PLUS_U_MULLIKEN_TO_HAMILTONIAN%HAMILTONIAN (   (1  ?      OCC_NUMBERS_PLUS_U_FULL M   g3  @     OCC_NUMBERS_PLUS_U_FULL%N_HAMILTONIAN_MATRIX_SIZE+DIMENSIONS N   �3  �   �  MPI_TASKS!OCC_NUMBERS_PLUS_U_FULL%MPIFCMB5+MPI_TASKS=MPIFCMB5 J   ,4  H     OCC_NUMBERS_PLUS_U_FULL%MPIFCMB5%MPI_UNWEIGHTED+MPI_TASKS N   t4  �   �  MPI_TASKS!OCC_NUMBERS_PLUS_U_FULL%MPIFCMB9+MPI_TASKS=MPIFCMB9 M   �4  H     OCC_NUMBERS_PLUS_U_FULL%MPIFCMB9%MPI_WEIGHTS_EMPTY+MPI_TASKS N   D5  �   �  MPI_TASKS!OCC_NUMBERS_PLUS_U_FULL%MPIPRIV1+MPI_TASKS=MPIPRIV1 F   06  H     OCC_NUMBERS_PLUS_U_FULL%MPIPRIV1%MPI_BOTTOM+MPI_TASKS H   x6  H     OCC_NUMBERS_PLUS_U_FULL%MPIPRIV1%MPI_IN_PLACE+MPI_TASKS M   �6  �     OCC_NUMBERS_PLUS_U_FULL%MPIPRIV1%MPI_STATUS_IGNORE+MPI_TASKS N   d7  �   �  MPI_TASKS!OCC_NUMBERS_PLUS_U_FULL%MPIPRIV2+MPI_TASKS=MPIPRIV2 O   (8  �     OCC_NUMBERS_PLUS_U_FULL%MPIPRIV2%MPI_STATUSES_IGNORE+MPI_TASKS O   �8  �     OCC_NUMBERS_PLUS_U_FULL%MPIPRIV2%MPI_ERRCODES_IGNORE+MPI_TASKS N   �9  �   �  MPI_TASKS!OCC_NUMBERS_PLUS_U_FULL%MPIPRIVC+MPI_TASKS=MPIPRIVC J   I:  �     OCC_NUMBERS_PLUS_U_FULL%MPIPRIVC%MPI_ARGVS_NULL+MPI_TASKS I   ;  �     OCC_NUMBERS_PLUS_U_FULL%MPIPRIVC%MPI_ARGV_NULL+MPI_TASKS >   �;  �   a   OCC_NUMBERS_PLUS_U_FULL%DENSITY_MATRIX_SPARSE 7   E<  �   a   OCC_NUMBERS_PLUS_U_FULL%DENSITY_MATRIX /   �<  @   a   OCC_NUMBERS_PLUS_U_FULL%I_SPIN /   )=  s      ADD_PLUS_U_FULL_TO_HAMILTONIAN T   �>  @     ADD_PLUS_U_FULL_TO_HAMILTONIAN%N_HAMILTONIAN_MATRIX_SIZE+DIMENSIONS A   �>  @     ADD_PLUS_U_FULL_TO_HAMILTONIAN%N_SPIN+DIMENSIONS ;   ?  �   a   ADD_PLUS_U_FULL_TO_HAMILTONIAN%HAMILTONIAN $   �?  �       PLUS_U_RAMPING_WORK )   r@  @   a   PLUS_U_RAMPING_WORK%ETOT .   �@  @   a   PLUS_U_RAMPING_WORK%ETOT_PREV 