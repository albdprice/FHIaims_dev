  �=  �   k820309    9          19.0        <��]                                                                                                          
       /home/albd/research/FHIaims/running_copy/FHIaims/src/esp_grids.f90 ESP_GRIDS                                                     
       BOHR                      @                              
       USE_UNIT                                                  
                 
                 ��W��?        0.52917721D0           @ @                                                    @ @                                                 
                &                                                    @ @                                                 
                &                                                    @ @                                                 
                &                                                    @ @                                                 
                &                                                    @ @                               	                                   &                                                    @ @                              
                   
                &                   &                                                    @ @                                                                  &                                                    @ @                                                 
                &                                                    @ @                                                                  &                                                    @ @                                                                  &                                                    @ @                                                 
                &                                                    @ @                                                 
                &                   &                                                    @ @                                                 
                &                   &                                                    @ @                                                                  &                   &                                                    @ @                                                 
                &                   &                   &                   &                                                    @ @                                                 
                &                   &                   &                                                    @ @                                                                  &                                                    @ @                                                 
                &                   &                   &                                                    @ @                                                 
                &                   &                                                    @ @                                                                  &                                                    @ @                                                 
                &                   &                   &                                                    @ @                                                 
                &                   &                                                    @ @                                                                  &                                                    @ @                                                                  &                   &                                                    @ @                                                                  &                   &                                                    @ @                                                                  &                   &                                                    @ @                                                                  &                   &                                                    @ @                                                                   &                   &                   &                                                    @ @                              !                   
                &                   &                   &                                                      @                                 "                       @                                 #                                                        $                              @                          %     '(                    #COORDS_ESP &   #INDEX_ATOM_ESP '   #INDEX_RADIAL_ESP (   #INDEX_ANGULAR_ESP )                �                              &                              
  p          p            p                                       �                               '                               �                               (                               �                               )                                     @                           *     '�                    #SIZE_ESP +   #POINTS_ESP ,   #BATCH_N_COMPUTE_ESP -   #BATCH_I_BASIS_ESP .                �                               +                               �                              ,                   (             #GRID_POINT_ESP %             &                                                        �                               -     P                         �                              .            X                             &                                                     @                                /            �                        &                                           #BATCH_OF_POINTS_ESP *                                               0                       @                                 1                                                         2                                                        3                                                        4                       @ @                               5                       @ @                               6                       @ @                               7                                                        8            #         @                                   9                     #         @                                   :                   #GET_GRIDS_ESP%N_MAX_ANGULAR ;   #GET_GRIDS_ESP%N_SPECIES <   #RADIUS_ESP_MIN_NEW =   #RADIUS_ESP_MAX_NEW >   #OUT_GRIDS ?   #ESP_LOG_GRID_IN @                                               ;                                                        <                      @                              =                    
 "    p          5 r <       5 r <                                                              >                    
 #    p          5 r <       5 r <                                                                ?                                                       @            #         @                                  A                    #I_SPECIES B                                              B            %         @                               C                           #N_ANG D                                              D            #         @                                   E                   #GET_LEBEDEV_ESP%N_MAX_ANGULAR F   #N_ANG G   #I_SPECIES H   #I_RADIAL I                                               F                      D                                 G                                                       H                                                       I            %         @                                J                           #N_ANG K                                              K            %         @                                L                           #N_ANG M                                              M            %         @                                N                           #N_ANG O                                              O            #         @                                   P                    #TEST_RADIAL_GRID_ESP%N_MAX_RADIAL Q                                               Q            #         @                                   R                     %         @                               S                    
       #R_CURRENT T   #N_SCALE U   #R_SCALE V                                             T     
                                                  U                                                      V     
       %         @                               W                    
       #R_CURRENT X   #R_MIN Y   #SCALE Z                                             X     
                                                 Y     
                  @                              Z     
       %         @                               [                    
       #R_CURRENT \   #I_SPECIES ]                                             \     
                                                  ]            %         @                               ^                    
       #INDEX _   #SCALE `   #N_MAX a                                             _     
                                                 `     
                                                  a            #         @                                   b                   #GET_LOCAL_YLM_TAB_ESP%N_SPECIES c   #L_HARTREE d                                               c                     D @                               d                     -    p          5 r c       5 r c                     #         @                                   e                    #R_CURRENT f   #R_MIN h   #SCALE i   #N_POINTS g   #OUT_VALUE j                                            f                    
 0    p          5 � p        r g       5 � p        r g                                                               h     
                  @                              i     
                                                  g                     D                                j                    
 1    p          5 � p        r g       5 � p        r g                     #         @                                   k                   #ESP_GET_N_COMPUTE_DENS%BATCH_PERM l   #ESP_GET_N_COMPUTE_DENS%N_CENTERS_INTEGRALS m   #ESP_GET_N_COMPUTE_DENS%N_CENTERS n   #ESP_GET_N_COMPUTE_DENS%N_BASIS_FNS o   #PARTITION_TAB q                                                                                                                    l            �             p          p            p                          #ESP_GET_N_COMPUTE_DENS%BATCH_PERMUTATION p              @ @                              m                                                        n                                                        o                                                     q                    
 2    p          5 r 4       5 r 4                                      @                          p     '�                   #INITIALIZED r   #N_MY_BATCHES s   #N_FULL_POINTS t   #BATCHES u   #PERM_BATCH_OWNER �   #POINT_SEND_CNT �   #POINT_SEND_OFF �   #POINT_RECV_CNT �   #POINT_RECV_OFF �   #N_BASIS_LOCAL �   #N_LOCAL_MATRIX_SIZE �   #I_BASIS_LOCAL �   #I_BASIS_GLB_TO_LOC �   #PARTITION_TAB �                � $                              r                                � $                              s                               � $                              t                              �$                              u                   �             #ESP_GET_N_COMPUTE_DENS%BATCH_OF_POINTS v             &                                                          @  @                         v     '�                    #SIZE w   #POINTS x   #BATCH_N_COMPUTE ~   #BATCH_I_BASIS                 �                               w                               �                              x                   (             #ESP_GET_N_COMPUTE_DENS%GRID_POINT y             &                                                          @  @                         y     '(                    #COORDS z   #INDEX_ATOM {   #INDEX_RADIAL |   #INDEX_ANGULAR }                �                              z                              
  p          p            p                                       �                               {                               �                               |                               �                               }                                �                               ~     P                         �                                          X                             &                                                       �$                              �            X                             &                                                       �$                              �            �                             &                                                       �$                              �            �                             &                                                       �$                              �            0                            &                                                       �$                              �            x             	               &                                                        � $                              �     �      
                   � $                              �     �                        �$                              �            �                            &                                                       �$                              �                                        &                                                       �$                             �            X                
            &                                              �   U      fn#fn    �   E   J  CONSTANTS    :  I   J  LOCALORB_IO    �  |       BOHR+CONSTANTS %   �  @       USE_UNIT+LOCALORB_IO    ?  �       R_GRID_MIN_ESP    �  �       R_GRID_MAX_ESP    W  �       R_GRID_INC_ESP #   �  �       LOG_R_GRID_INC_ESP    o  �       N_GRID_ESP    �  �       R_GRID_ESP    �  �       N_RADIAL_ESP !   +  �       SCALE_RADIAL_ESP "   �  �       ANGULAR_LIMIT_ESP     C  �       ANGULAR_MIN_ESP     �  �       ANGULAR_ACC_ESP    [  �       R_RADIAL_ESP    �  �       W_RADIAL_ESP    �	  �       N_ANGULAR_ESP    G
  �       R_ANGULAR_ESP      �       W_ANGULAR_ESP "   �  �       N_ANG_LEBEDEV_ESP "   c  �       R_ANG_LEBEDEV_ESP "     �       W_ANG_LEBEDEV_ESP &   �  �       N_ANGULAR_LEBEDEV_ESP &   O  �       R_ANGULAR_LEBEDEV_ESP &     �       W_ANGULAR_LEBEDEV_ESP '   �  �       N_DIVISION_LEBEDEV_ESP 0   ;  �       DIVISION_BOUNDARIES_LEBEDEV_ESP %   �  �       FIXED_GRID_INDEX_ESP '   �  �       LEBEDEV_GRID_INDEX_ESP    '  �       N_DIVISION_ESP (   �  �       DIVISION_BOUNDARIES_ESP "   �  �       LOCAL_YLM_TAB_ESP "   C  @       N_MAX_LEBEDEV_ESP    �  @       N_GRIDS_ESP %   �  @       GRID_PARTITIONED_ESP      �       GRID_POINT_ESP *   �  �   a   GRID_POINT_ESP%COORDS_ESP .   @  H   a   GRID_POINT_ESP%INDEX_ATOM_ESP 0   �  H   a   GRID_POINT_ESP%INDEX_RADIAL_ESP 1   �  H   a   GRID_POINT_ESP%INDEX_ANGULAR_ESP $     �       BATCH_OF_POINTS_ESP -   �  H   a   BATCH_OF_POINTS_ESP%SIZE_ESP /   �  �   a   BATCH_OF_POINTS_ESP%POINTS_ESP 8   �  H   a   BATCH_OF_POINTS_ESP%BATCH_N_COMPUTE_ESP 6   �  �   a   BATCH_OF_POINTS_ESP%BATCH_I_BASIS_ESP    �  �       BATCHES_ESP #   '  @       N_GRID_BATCHES_ESP &   g  @       N_POINTS_IN_BATCH_ESP !   �  @       N_MY_BATCHES_ESP %   �  @       N_MAX_BATCH_SIZE_ESP "   '  @       N_FULL_POINTS_ESP (   g  @       N_MAX_COMPUTE_ATOMS_ESP '   �  @       N_MAX_COMPUTE_DENS_ESP +   �  @       N_MAX_COMPUTE_FNS_DENS_ESP (   '  @       GOT_N_COMPUTE_MAXES_ESP #   g  H       ALLOCATE_GRIDS_ESP    �  �       GET_GRIDS_ESP 7   �  @     GET_GRIDS_ESP%N_MAX_ANGULAR+DIMENSIONS 3   �  @     GET_GRIDS_ESP%N_SPECIES+DIMENSIONS 1   	  �   a   GET_GRIDS_ESP%RADIUS_ESP_MIN_NEW 1   �  �   a   GET_GRIDS_ESP%RADIUS_ESP_MAX_NEW (   1  @   a   GET_GRIDS_ESP%OUT_GRIDS .   q  @   a   GET_GRIDS_ESP%ESP_LOG_GRID_IN )   �  W       GET_LOGARITHMIC_GRID_ESP 3      @   a   GET_LOGARITHMIC_GRID_ESP%I_SPECIES '   H   [       LEBEDEV_GRID_FLOOR_ESP -   �   @   a   LEBEDEV_GRID_FLOOR_ESP%N_ANG     �   �       GET_LEBEDEV_ESP 9   v!  @     GET_LEBEDEV_ESP%N_MAX_ANGULAR+DIMENSIONS &   �!  @   a   GET_LEBEDEV_ESP%N_ANG *   �!  @   a   GET_LEBEDEV_ESP%I_SPECIES )   6"  @   a   GET_LEBEDEV_ESP%I_RADIAL &   v"  [       LEBEDEV_GRID_CEIL_ESP ,   �"  @   a   LEBEDEV_GRID_CEIL_ESP%N_ANG    #  [       GRID_CEIL_ESP $   l#  @   a   GRID_CEIL_ESP%N_ANG    �#  [       GRID_FLOOR_ESP %   $  @   a   GRID_FLOOR_ESP%N_ANG %   G$  o       TEST_RADIAL_GRID_ESP =   �$  @     TEST_RADIAL_GRID_ESP%N_MAX_RADIAL+DIMENSIONS "   �$  H       CLEANUP_GRIDS_ESP '   >%  y       INVERT_RADIAL_GRID_ESP 1   �%  @   a   INVERT_RADIAL_GRID_ESP%R_CURRENT /   �%  @   a   INVERT_RADIAL_GRID_ESP%N_SCALE /   7&  @   a   INVERT_RADIAL_GRID_ESP%R_SCALE $   w&  u       INVERT_LOG_GRID_ESP .   �&  @   a   INVERT_LOG_GRID_ESP%R_CURRENT *   ,'  @   a   INVERT_LOG_GRID_ESP%R_MIN *   l'  @   a   INVERT_LOG_GRID_ESP%SCALE '   �'  n       INVERT_LOG_GRID_P2_ESP 1   (  @   a   INVERT_LOG_GRID_P2_ESP%R_CURRENT 1   Z(  @   a   INVERT_LOG_GRID_P2_ESP%I_SPECIES &   �(  q       GET_RADIAL_WEIGHT_ESP ,   )  @   a   GET_RADIAL_WEIGHT_ESP%INDEX ,   K)  @   a   GET_RADIAL_WEIGHT_ESP%SCALE ,   �)  @   a   GET_RADIAL_WEIGHT_ESP%N_MAX &   �)  |       GET_LOCAL_YLM_TAB_ESP ;   G*  @     GET_LOCAL_YLM_TAB_ESP%N_SPECIES+DIMENSIONS 0   �*  �   a   GET_LOCAL_YLM_TAB_ESP%L_HARTREE +   +  �       INVERT_LOG_GRID_VECTOR_ESP 5   �+  �   a   INVERT_LOG_GRID_VECTOR_ESP%R_CURRENT 1   Y,  @   a   INVERT_LOG_GRID_VECTOR_ESP%R_MIN 1   �,  @   a   INVERT_LOG_GRID_VECTOR_ESP%SCALE 4   �,  @   a   INVERT_LOG_GRID_VECTOR_ESP%N_POINTS 5   -  �   a   INVERT_LOG_GRID_VECTOR_ESP%OUT_VALUE '   �-  G      ESP_GET_N_COMPUTE_DENS L   /  �     ESP_GET_N_COMPUTE_DENS%BATCH_PERM+LOAD_BALANCING=BATCH_PERM F   �/  @     ESP_GET_N_COMPUTE_DENS%N_CENTERS_INTEGRALS+DIMENSIONS <   0  @     ESP_GET_N_COMPUTE_DENS%N_CENTERS+DIMENSIONS >   V0  @     ESP_GET_N_COMPUTE_DENS%N_BASIS_FNS+DIMENSIONS 5   �0  �   a   ESP_GET_N_COMPUTE_DENS%PARTITION_TAB Z   *1  c     ESP_GET_N_COMPUTE_DENS%BATCH_PERMUTATION+LOAD_BALANCING=BATCH_PERMUTATION T   �2  H   a   ESP_GET_N_COMPUTE_DENS%BATCH_PERMUTATION%INITIALIZED+LOAD_BALANCING U   �2  H   a   ESP_GET_N_COMPUTE_DENS%BATCH_PERMUTATION%N_MY_BATCHES+LOAD_BALANCING V   3  H   a   ESP_GET_N_COMPUTE_DENS%BATCH_PERMUTATION%N_FULL_POINTS+LOAD_BALANCING P   e3  �   a   ESP_GET_N_COMPUTE_DENS%BATCH_PERMUTATION%BATCHES+LOAD_BALANCING M   %4  �      ESP_GET_N_COMPUTE_DENS%BATCH_OF_POINTS+GRIDS=BATCH_OF_POINTS B   �4  H   a   ESP_GET_N_COMPUTE_DENS%BATCH_OF_POINTS%SIZE+GRIDS D   �4  �   a   ESP_GET_N_COMPUTE_DENS%BATCH_OF_POINTS%POINTS+GRIDS C   �5  �      ESP_GET_N_COMPUTE_DENS%GRID_POINT+GRIDS=GRID_POINT ?   G6  �   a   ESP_GET_N_COMPUTE_DENS%GRID_POINT%COORDS+GRIDS C   �6  H   a   ESP_GET_N_COMPUTE_DENS%GRID_POINT%INDEX_ATOM+GRIDS E   +7  H   a   ESP_GET_N_COMPUTE_DENS%GRID_POINT%INDEX_RADIAL+GRIDS F   s7  H   a   ESP_GET_N_COMPUTE_DENS%GRID_POINT%INDEX_ANGULAR+GRIDS M   �7  H   a   ESP_GET_N_COMPUTE_DENS%BATCH_OF_POINTS%BATCH_N_COMPUTE+GRIDS K   8  �   a   ESP_GET_N_COMPUTE_DENS%BATCH_OF_POINTS%BATCH_I_BASIS+GRIDS Y   �8  �   a   ESP_GET_N_COMPUTE_DENS%BATCH_PERMUTATION%PERM_BATCH_OWNER+LOAD_BALANCING W   +9  �   a   ESP_GET_N_COMPUTE_DENS%BATCH_PERMUTATION%POINT_SEND_CNT+LOAD_BALANCING W   �9  �   a   ESP_GET_N_COMPUTE_DENS%BATCH_PERMUTATION%POINT_SEND_OFF+LOAD_BALANCING W   S:  �   a   ESP_GET_N_COMPUTE_DENS%BATCH_PERMUTATION%POINT_RECV_CNT+LOAD_BALANCING W   �:  �   a   ESP_GET_N_COMPUTE_DENS%BATCH_PERMUTATION%POINT_RECV_OFF+LOAD_BALANCING V   {;  H   a   ESP_GET_N_COMPUTE_DENS%BATCH_PERMUTATION%N_BASIS_LOCAL+LOAD_BALANCING \   �;  H   a   ESP_GET_N_COMPUTE_DENS%BATCH_PERMUTATION%N_LOCAL_MATRIX_SIZE+LOAD_BALANCING V   <  �   a   ESP_GET_N_COMPUTE_DENS%BATCH_PERMUTATION%I_BASIS_LOCAL+LOAD_BALANCING [   �<  �   a   ESP_GET_N_COMPUTE_DENS%BATCH_PERMUTATION%I_BASIS_GLB_TO_LOC+LOAD_BALANCING V   3=  �   a   ESP_GET_N_COMPUTE_DENS%BATCH_PERMUTATION%PARTITION_TAB+LOAD_BALANCING 