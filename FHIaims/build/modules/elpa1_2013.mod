  ŚL  °   k820309    9          19.0        ĂŘÁ]                                                                                                          
       /home/albd/research/FHIaims/running_copy/FHIaims/src/elpa/elpa1_st.f90 ELPA1_2013              GET_ELPA_ROW_COL_COMMS_2013 SOLVE_EVP_REAL_2013 SOLVE_EVP_COMPLEX_2013 TRIDIAG_REAL_2013 TRANS_EV_REAL_2013 MULT_AT_B_REAL_2013 TRIDIAG_COMPLEX_2013 TRANS_EV_COMPLEX_2013 MULT_AH_B_COMPLEX_2013 SOLVE_TRIDI_2013 CHOLESKY_REAL_2013 INVERT_TRM_REAL_2013 CHOLESKY_COMPLEX_2013 INVERT_TRM_COMPLEX_2013 LOCAL_INDEX LEAST_COMMON_MULTIPLE HH_TRANSFORM_REAL HH_TRANSFORM_COMPLEX TIME_EVP_FWD TIME_EVP_SOLVE TIME_EVP_BACK ELPA_PRINT_TIMES #         @                                                       #MPI_COMM_WORLD    #MY_PROW    #MY_PCOL    #MPI_COMM_ROWS    #MPI_COMM_COLS              
@ @                                                    
@ @                                                    
@ @                                                    D @                                                     D @                                           #         @                                                       #NA    #NEV 	   #A 
   #LDA    #EV    #Q    #LDQ    #NBLK    #MPI_COMM_ROWS    #MPI_COMM_COLS    #MPI_COMM_GLOBAL              
@ @                                                    
@ @                               	                  B  D @                              
                    
       p        5  p        r    p          5  p        r      1     5  p        r      1                             
@ @                                                   D @                                                  
     p          5  p        r        5  p        r                             B  D @                                                  
       p        5  p        r    p          5  p        r      1     5  p        r      1                             
@ @                                                    
@ @                                                    
@ @                                                    
@ @                                                                                                 #         @                                                    
   #NA    #NEV    #A    #LDA    #EV    #Q    #LDQ    #NBLK    #MPI_COMM_ROWS    #MPI_COMM_COLS              
@ @                                                    
@ @                                                 B  D @                                                         p        5  p        r    p          5  p        r      1     5  p        r      1                             
@ @                                                   D @                                                  
     p          5  p        r        5  p        r                             B  D @                                                         p        5  p        r    p          5  p        r      1     5  p        r      1                             
@ @                                                    
@ @                                                    
@ @                                                    
@ @                                          #         @                                                   	   #NA    #A     #LDA !   #NBLK "   #MPI_COMM_ROWS #   #MPI_COMM_COLS $   #D %   #E &   #TAU '             D @                                                  B  D @                                                   
       p        5  p        r !   p          5  p        r !     1     5  p        r !     1                             D @                               !                      D @                               "                      D @                               #                      D @                               $                     D @                              %                    
     p          5  p        r        5  p        r                               D @                              &                    
     p          5  p        r        5  p        r                               D @                              '                    
     p          5  p        r        5  p        r                      #         @                                  (                 
   #NA )   #NQC *   #A +   #LDA ,   #TAU -   #Q .   #LDQ /   #NBLK 0   #MPI_COMM_ROWS 1   #MPI_COMM_COLS 2              @                               )                      D @                               *                   B                                  +                    
       p        5  p        r ,   p          5  p        r ,     1     5  p        r ,     1                                                              ,                                                     -                    
 !    p          5  p        r )       5  p        r )                            B  D @                              .                    
        p        5  p        r /   p          5  p        r /     1     5  p        r /     1                             D @                               /                      D @                               0                      D @                               1                      D @                               2            #         @                                   3                    #UPLO_A 4   #UPLO_C 5   #NA 6   #NCB 7   #A 8   #LDA 9   #B :   #LDB ;   #NBLK <   #MPI_COMM_ROWS =   #MPI_COMM_COLS >   #C ?   #LDC @                                              4                                                                       5                                      D @                               6                      D @                               7                   B                                  8                    
 )      p        5  p        r 9   p          5  p        r 9     1     5  p        r 9     1                                                              9                   B  D @                              :                    
 *      p        5  p        r ;   p          5  p        r ;     1     5  p        r ;     1                             D @                               ;                      D @                               <                      D @                               =                      D @                               >                   B  D                                ?                    
 +      p        5  p        r @   p          5  p        r @     1     5  p        r @     1                                                              @            #         @                                  A                 	   #NA B   #A C   #LDA D   #NBLK E   #MPI_COMM_ROWS F   #MPI_COMM_COLS G   #D H   #E I   #TAU J             D @                               B                   B  D @                              C                     2      p        5  p        r D   p          5  p        r D     1     5  p        r D     1                             D @                               D                      D @                               E                      D @                               F                      D @                               G                     D @                              H                    
 4    p          5  p        r B       5  p        r B                              D @                              I                    
 5    p          5  p        r B       5  p        r B                              D @                              J                     3    p          5  p        r B       5  p        r B                     #         @                                  K                 
   #NA L   #NQC M   #A N   #LDA O   #TAU P   #Q Q   #LDQ R   #NBLK S   #MPI_COMM_ROWS T   #MPI_COMM_COLS U              @                               L                      D @                               M                   B                                  N                     A      p        5  p        r O   p          5  p        r O     1     5  p        r O     1                                                              O                                                     P                     C    p          5  p        r L       5  p        r L                            B  D @                              Q                     B      p        5  p        r R   p          5  p        r R     1     5  p        r R     1                             D @                               R                      D @                               S                      D @                               T                      D @                               U            #         @                                   V                    #UPLO_A W   #UPLO_C X   #NA Y   #NCB Z   #A [   #LDA \   #B ]   #LDB ^   #NBLK _   #MPI_COMM_ROWS `   #MPI_COMM_COLS a   #C b   #LDC c                                              W                                                                       X                                      D @                               Y                      D @                               Z                   B                                  [                     K      p        5  p        r \   p          5  p        r \     1     5  p        r \     1                                                              \                   B  D @                              ]                     L      p        5  p        r ^   p          5  p        r ^     1     5  p        r ^     1                             D @                               ^                      D @                               _                      D @                               `                      D @                               a                   B  D                                b                     M      p        5  p        r c   p          5  p        r c     1     5  p        r c     1                                                              c            #         @                                  d                 	   #NA e   #NEV f   #D g   #E h   #Q i   #LDQ j   #NBLK k   #MPI_COMM_ROWS l   #MPI_COMM_COLS m             D @                               e                       @                               f                     D @                              g                    
 T    p          5  p        r e       5  p        r e                              D @                              h                    
 U    p          5  p        r e       5  p        r e                            B  D @                              i                    
 V      p        5  p        r j   p          5  p        r j     1     5  p        r j     1                             D @                               j                      D @                               k                      D @                               l                      D @                               m            #         @                                   n                    #NA o   #A p   #LDA q   #NBLK r   #MPI_COMM_ROWS s   #MPI_COMM_COLS t             D @                               o                   B  D @                              p                    
       p        5  p        r q   p          5  p        r q     1     5  p        r q     1                             D @                               q                      D @                               r                      D @                               s                      D @                               t            #         @                                   u                    #NA v   #A w   #LDA x   #NBLK y   #MPI_COMM_ROWS z   #MPI_COMM_COLS {             D @                               v                   B  D @                              w                    
       p        5  p        r x   p          5  p        r x     1     5  p        r x     1                             D @                               x                      D @                               y                      D @                               z                      D @                               {            #         @                                   |                    #NA }   #A ~   #LDA    #NBLK    #MPI_COMM_ROWS    #MPI_COMM_COLS              D @                               }                   B  D @                              ~                     Ł      p        5  p        r    p          5  p        r      1     5  p        r      1                             D @                                                     D @                                                     D @                                                     D @                                           #         @                                                       #NA    #A    #LDA    #NBLK    #MPI_COMM_ROWS    #MPI_COMM_COLS              D @                                                  B  D @                                                   ¨      p        5  p        r    p          5  p        r      1     5  p        r      1                             D @                                                     D @                                                     D @                                                     D @                                           %         @                                                          #IDX    #MY_PROC    #NUM_PROCS    #NBLK    #IFLAG                                                                                                                             @                                                      @                                                                                                  %         @                                                          #A    #B              
                                                       
  @                                          #         @                                                      #ALPHA    #XNORM_SQ    #XF    #TAU              
D @                                   
                 
                                      
                D                                     
                 D                                     
       #         @                                                      #ALPHA    #XNORM_SQ    #XF    #TAU              
D @                                                    
                                      
                D                                                      D                                                       @                                     
                  @                                     
                  @                                     
                                                                              @                           Ą                          #MPI_UNWEIGHTED ˘                @                           ˘                                  @                           Ł                          #MPI_WEIGHTS_EMPTY ¤                @                           ¤                                  @                           Ľ                          #MPI_BOTTOM Ś   #MPI_IN_PLACE §   #MPI_STATUS_IGNORE ¨                @                           Ś                                 @                           §                                @                           ¨                                p          p            p                                                @                           Š                          #MPI_STATUSES_IGNORE Ş   #MPI_ERRCODES_IGNORE Ť                @                           Ş                                 p          p          p            p          p                                               @                           Ť                                p          p            p                                                @                           Ź                          #MPI_ARGVS_NULL ­   #MPI_ARGV_NULL Ž   -             @                           ­                                 p          p          p            p          p                                  -             @                           Ž                                p          p            p                                         Z      fn#fn     ú   ˝  b   uapp(ELPA1_2013 ,   ˇ         GET_ELPA_ROW_COL_COMMS_2013 ;   S  @   a   GET_ELPA_ROW_COL_COMMS_2013%MPI_COMM_WORLD 4     @   a   GET_ELPA_ROW_COL_COMMS_2013%MY_PROW 4   Ó  @   a   GET_ELPA_ROW_COL_COMMS_2013%MY_PCOL :     @   a   GET_ELPA_ROW_COL_COMMS_2013%MPI_COMM_ROWS :   S  @   a   GET_ELPA_ROW_COL_COMMS_2013%MPI_COMM_COLS $     Ć       SOLVE_EVP_REAL_2013 '   Y  @   a   SOLVE_EVP_REAL_2013%NA (     @   a   SOLVE_EVP_REAL_2013%NEV &   Ů  ô   a   SOLVE_EVP_REAL_2013%A (   Í  @   a   SOLVE_EVP_REAL_2013%LDA '     ´   a   SOLVE_EVP_REAL_2013%EV &   Á  ô   a   SOLVE_EVP_REAL_2013%Q (   ľ  @   a   SOLVE_EVP_REAL_2013%LDQ )   ő  @   a   SOLVE_EVP_REAL_2013%NBLK 2   5	  @   a   SOLVE_EVP_REAL_2013%MPI_COMM_ROWS 2   u	  @   a   SOLVE_EVP_REAL_2013%MPI_COMM_COLS 4   ľ	  @   a   SOLVE_EVP_REAL_2013%MPI_COMM_GLOBAL '   ő	  ą       SOLVE_EVP_COMPLEX_2013 *   Ś
  @   a   SOLVE_EVP_COMPLEX_2013%NA +   ć
  @   a   SOLVE_EVP_COMPLEX_2013%NEV )   &  ô   a   SOLVE_EVP_COMPLEX_2013%A +     @   a   SOLVE_EVP_COMPLEX_2013%LDA *   Z  ´   a   SOLVE_EVP_COMPLEX_2013%EV )     ô   a   SOLVE_EVP_COMPLEX_2013%Q +     @   a   SOLVE_EVP_COMPLEX_2013%LDQ ,   B  @   a   SOLVE_EVP_COMPLEX_2013%NBLK 5     @   a   SOLVE_EVP_COMPLEX_2013%MPI_COMM_ROWS 5   Â  @   a   SOLVE_EVP_COMPLEX_2013%MPI_COMM_COLS "     §       TRIDIAG_REAL_2013 %   Š  @   a   TRIDIAG_REAL_2013%NA $   é  ô   a   TRIDIAG_REAL_2013%A &   Ý  @   a   TRIDIAG_REAL_2013%LDA '     @   a   TRIDIAG_REAL_2013%NBLK 0   ]  @   a   TRIDIAG_REAL_2013%MPI_COMM_ROWS 0     @   a   TRIDIAG_REAL_2013%MPI_COMM_COLS $   Ý  ´   a   TRIDIAG_REAL_2013%D $     ´   a   TRIDIAG_REAL_2013%E &   E  ´   a   TRIDIAG_REAL_2013%TAU #   ů  ˛       TRANS_EV_REAL_2013 &   Ť  @   a   TRANS_EV_REAL_2013%NA '   ë  @   a   TRANS_EV_REAL_2013%NQC %   +  ô   a   TRANS_EV_REAL_2013%A '     @   a   TRANS_EV_REAL_2013%LDA '   _  ´   a   TRANS_EV_REAL_2013%TAU %     ô   a   TRANS_EV_REAL_2013%Q '     @   a   TRANS_EV_REAL_2013%LDQ (   G  @   a   TRANS_EV_REAL_2013%NBLK 1     @   a   TRANS_EV_REAL_2013%MPI_COMM_ROWS 1   Ç  @   a   TRANS_EV_REAL_2013%MPI_COMM_COLS $     Ń       MULT_AT_B_REAL_2013 +   Ř  P   a   MULT_AT_B_REAL_2013%UPLO_A +   (  P   a   MULT_AT_B_REAL_2013%UPLO_C '   x  @   a   MULT_AT_B_REAL_2013%NA (   ¸  @   a   MULT_AT_B_REAL_2013%NCB &   ř  ô   a   MULT_AT_B_REAL_2013%A (   ě  @   a   MULT_AT_B_REAL_2013%LDA &   ,  ô   a   MULT_AT_B_REAL_2013%B (      @   a   MULT_AT_B_REAL_2013%LDB )   `  @   a   MULT_AT_B_REAL_2013%NBLK 2      @   a   MULT_AT_B_REAL_2013%MPI_COMM_ROWS 2   ŕ  @   a   MULT_AT_B_REAL_2013%MPI_COMM_COLS &      ô   a   MULT_AT_B_REAL_2013%C (     @   a   MULT_AT_B_REAL_2013%LDC %   T  §       TRIDIAG_COMPLEX_2013 (   ű  @   a   TRIDIAG_COMPLEX_2013%NA '   ;   ô   a   TRIDIAG_COMPLEX_2013%A )   /!  @   a   TRIDIAG_COMPLEX_2013%LDA *   o!  @   a   TRIDIAG_COMPLEX_2013%NBLK 3   Ż!  @   a   TRIDIAG_COMPLEX_2013%MPI_COMM_ROWS 3   ď!  @   a   TRIDIAG_COMPLEX_2013%MPI_COMM_COLS '   /"  ´   a   TRIDIAG_COMPLEX_2013%D '   ă"  ´   a   TRIDIAG_COMPLEX_2013%E )   #  ´   a   TRIDIAG_COMPLEX_2013%TAU &   K$  ˛       TRANS_EV_COMPLEX_2013 )   ý$  @   a   TRANS_EV_COMPLEX_2013%NA *   =%  @   a   TRANS_EV_COMPLEX_2013%NQC (   }%  ô   a   TRANS_EV_COMPLEX_2013%A *   q&  @   a   TRANS_EV_COMPLEX_2013%LDA *   ą&  ´   a   TRANS_EV_COMPLEX_2013%TAU (   e'  ô   a   TRANS_EV_COMPLEX_2013%Q *   Y(  @   a   TRANS_EV_COMPLEX_2013%LDQ +   (  @   a   TRANS_EV_COMPLEX_2013%NBLK 4   Ů(  @   a   TRANS_EV_COMPLEX_2013%MPI_COMM_ROWS 4   )  @   a   TRANS_EV_COMPLEX_2013%MPI_COMM_COLS '   Y)  Ń       MULT_AH_B_COMPLEX_2013 .   **  P   a   MULT_AH_B_COMPLEX_2013%UPLO_A .   z*  P   a   MULT_AH_B_COMPLEX_2013%UPLO_C *   Ę*  @   a   MULT_AH_B_COMPLEX_2013%NA +   
+  @   a   MULT_AH_B_COMPLEX_2013%NCB )   J+  ô   a   MULT_AH_B_COMPLEX_2013%A +   >,  @   a   MULT_AH_B_COMPLEX_2013%LDA )   ~,  ô   a   MULT_AH_B_COMPLEX_2013%B +   r-  @   a   MULT_AH_B_COMPLEX_2013%LDB ,   ˛-  @   a   MULT_AH_B_COMPLEX_2013%NBLK 5   ň-  @   a   MULT_AH_B_COMPLEX_2013%MPI_COMM_ROWS 5   2.  @   a   MULT_AH_B_COMPLEX_2013%MPI_COMM_COLS )   r.  ô   a   MULT_AH_B_COMPLEX_2013%C +   f/  @   a   MULT_AH_B_COMPLEX_2013%LDC !   Ś/  §       SOLVE_TRIDI_2013 $   M0  @   a   SOLVE_TRIDI_2013%NA %   0  @   a   SOLVE_TRIDI_2013%NEV #   Í0  ´   a   SOLVE_TRIDI_2013%D #   1  ´   a   SOLVE_TRIDI_2013%E #   52  ô   a   SOLVE_TRIDI_2013%Q %   )3  @   a   SOLVE_TRIDI_2013%LDQ &   i3  @   a   SOLVE_TRIDI_2013%NBLK /   Š3  @   a   SOLVE_TRIDI_2013%MPI_COMM_ROWS /   é3  @   a   SOLVE_TRIDI_2013%MPI_COMM_COLS #   )4         CHOLESKY_REAL_2013 &   š4  @   a   CHOLESKY_REAL_2013%NA %   ů4  ô   a   CHOLESKY_REAL_2013%A '   í5  @   a   CHOLESKY_REAL_2013%LDA (   -6  @   a   CHOLESKY_REAL_2013%NBLK 1   m6  @   a   CHOLESKY_REAL_2013%MPI_COMM_ROWS 1   ­6  @   a   CHOLESKY_REAL_2013%MPI_COMM_COLS %   í6         INVERT_TRM_REAL_2013 (   }7  @   a   INVERT_TRM_REAL_2013%NA '   ˝7  ô   a   INVERT_TRM_REAL_2013%A )   ą8  @   a   INVERT_TRM_REAL_2013%LDA *   ń8  @   a   INVERT_TRM_REAL_2013%NBLK 3   19  @   a   INVERT_TRM_REAL_2013%MPI_COMM_ROWS 3   q9  @   a   INVERT_TRM_REAL_2013%MPI_COMM_COLS &   ą9         CHOLESKY_COMPLEX_2013 )   A:  @   a   CHOLESKY_COMPLEX_2013%NA (   :  ô   a   CHOLESKY_COMPLEX_2013%A *   u;  @   a   CHOLESKY_COMPLEX_2013%LDA +   ľ;  @   a   CHOLESKY_COMPLEX_2013%NBLK 4   ő;  @   a   CHOLESKY_COMPLEX_2013%MPI_COMM_ROWS 4   5<  @   a   CHOLESKY_COMPLEX_2013%MPI_COMM_COLS (   u<         INVERT_TRM_COMPLEX_2013 +   =  @   a   INVERT_TRM_COMPLEX_2013%NA *   E=  ô   a   INVERT_TRM_COMPLEX_2013%A ,   9>  @   a   INVERT_TRM_COMPLEX_2013%LDA -   y>  @   a   INVERT_TRM_COMPLEX_2013%NBLK 6   š>  @   a   INVERT_TRM_COMPLEX_2013%MPI_COMM_ROWS 6   ů>  @   a   INVERT_TRM_COMPLEX_2013%MPI_COMM_COLS    9?         LOCAL_INDEX     Ă?  @   a   LOCAL_INDEX%IDX $   @  @   a   LOCAL_INDEX%MY_PROC &   C@  @   a   LOCAL_INDEX%NUM_PROCS !   @  @   a   LOCAL_INDEX%NBLK "   Ă@  @   a   LOCAL_INDEX%IFLAG &   A  ^       LEAST_COMMON_MULTIPLE (   aA  @   a   LEAST_COMMON_MULTIPLE%A (   ĄA  @   a   LEAST_COMMON_MULTIPLE%B "   áA  r       HH_TRANSFORM_REAL (   SB  @   a   HH_TRANSFORM_REAL%ALPHA +   B  @   a   HH_TRANSFORM_REAL%XNORM_SQ %   ÓB  @   a   HH_TRANSFORM_REAL%XF &   C  @   a   HH_TRANSFORM_REAL%TAU %   SC  r       HH_TRANSFORM_COMPLEX +   ĹC  @   a   HH_TRANSFORM_COMPLEX%ALPHA .   D  @   a   HH_TRANSFORM_COMPLEX%XNORM_SQ (   ED  @   a   HH_TRANSFORM_COMPLEX%XF )   D  @   a   HH_TRANSFORM_COMPLEX%TAU    ĹD  @       TIME_EVP_FWD    E  @       TIME_EVP_SOLVE    EE  @       TIME_EVP_BACK !   E  @       ELPA_PRINT_TIMES $   ĹE  d      ELPA1_2013!MPIFCMB5    )F  H      MPI_UNWEIGHTED $   qF  g      ELPA1_2013!MPIFCMB9 "   ŘF  H      MPI_WEIGHTS_EMPTY $    G        ELPA1_2013!MPIPRIV1    ŠG  H      MPI_BOTTOM    ńG  H      MPI_IN_PLACE "   9H  ¤      MPI_STATUS_IGNORE $   ÝH        ELPA1_2013!MPIPRIV2 $   _I  Ä      MPI_STATUSES_IGNORE $   #J  ¤      MPI_ERRCODES_IGNORE $   ÇJ  w      ELPA1_2013!MPIPRIVC    >K  Ä      MPI_ARGVS_NULL    L  ¤      MPI_ARGV_NULL 