  C&  Z   k820309    l          18.0        7�]                                                                                                          
       /home/lleblanc/src/psiesta-4.1-b4/Src/class_Distribution.F90 CLASS_DISTRIBUTION              DISTRIBUTION TYPE_NULL TYPE_BLOCK_CYCLIC TYPE_PEXSI gen@NEWDISTRIBUTION gen@DISTTYPE gen@REF_COMM gen@GET_RANKS_IN_REF_COMM gen@NUM_LOCAL_ELEMENTS gen@INDEX_LOCAL_TO_GLOBAL gen@INDEX_GLOBAL_TO_LOCAL gen@NODE_HANDLING_ELEMENT i@| gen@INIT gen@DELETE gen@REFCOUNT gen@ID gen@NAME gen@INITIALIZED gen@SAME                                                     
                                                              u #NEWDISTRIBUTION_    #         @     @X                                                 #THIS    #REF_COMM    #RANKS_IN_REF_COMM    #DIST_TYPE    #BLOCKSIZE    #NAME 	             
D @                                                   #DISTRIBUTION              
  @                                                    
 @                                                                &                                                     
                                                       
                                                       
 @                             	                    1                                                        u #DIST_TYPE 
   %         @    @X                            
                           #THIS                                                                   #DISTRIBUTION                                                           u #REF_COMM_    %         @    @X                                                       #THIS              
                                                     #DISTRIBUTION                                                           u #GET_RANKS_IN_REF_COMM_    #         @     @X                                                 #THIS    #RANKS              
                                                     #DISTRIBUTION            D                                                                   &                                                                                                  u #NUM_LOCAL_ELEMENTS_    %         @    @X                                                       #THIS    #NELS    #NODE              
                                                     #DISTRIBUTION              
@ @                                                    
@ @                                                                                                 u #INDEX_LOCAL_TO_GLOBAL_    %         @    @X                                                       #THIS    #IL    #NODE              
                                                     #DISTRIBUTION              
@ @                                                    
@ @                                                                                                 u #INDEX_GLOBAL_TO_LOCAL_    %         @    @X                                                       #THIS    #IG    #NODE              
  @                                                  #DISTRIBUTION              
@ @                                                    
B @                                                                                                 u #NODE_HANDLING_ELEMENT_    %         @    @X                                                      #THIS    #IG              
                                                     #DISTRIBUTION              
@ @                                                                                                    |  #ASSIGN_     #         @     @X                                                  #THIS !   #OTHER "             
D @                               !                    #DISTRIBUTION              
  @                               "                   #DISTRIBUTION                                                          u #INIT_ #   #         @     @X                            #                    #THIS $             
D @                               $                    #DISTRIBUTION                                                          u #DELETE_ %   #         @     @X                            %                    #THIS &             
D @                               &                    #DISTRIBUTION                                                           u #REFCOUNT_ '   %         @    @X                            '                           #THIS (             
                                  (                   #DISTRIBUTION                                                           u #ID_ )   $         @    @X                            )     $                      #THIS *                     
                                  *                   #DISTRIBUTION                                                           u #NAME_ +   $        @   @X                            +                            #THIS ,   I r -     5 8 8#DISTRIBUTION_ .    O#DISTRIBUTION     p        U #DISTRIBUTION_ .       /   U  .   0              
                                  ,                   #DISTRIBUTION                                                           u #INITIALIZED_ 1   %         @    @X                          1                           #THIS 2             
                                  2                   #DISTRIBUTION                                                           u #SAME_ 3   %         @    @X                           3                           #THIS1 4   #THIS2 5             
  @                               4                   #DISTRIBUTION              
  @                               5                   #DISTRIBUTION                                                 6                                                      2                                             7                                                       0                                             8                                          �������                                                    9                                                         :                                                       0                                             ;                                          ��������                                                     <                                                      1                                             =                                                      2              @  @              A           .     '�                   #REFCOUNT >   #ID ?   #NAME 0   #DIST_TYPE @   #REF_COMM A   #GROUP B   #NODE C   #RANKS_IN_REF_COMM D   #NODES E   #NODE_IO F   #BLOCKSIZE G   #ISRCPROC H               �                               >                                          �                                                      0                 �                              ?     $                                �                              0            (                                     �                                                      Cnull Dist                                                                                                                                                                                                                                                                                        �                               @     (                        �                               A     ,                                   �                                                      2                �                               B     0                                   �                                                       0                �                               C     4                                   �                                         �������                      �                               D            8                            &                                                       �                               E     �      	                             �                                                      0                �                               F     �      
                             �                                         ��������                        �                               G     �                                   �                                                      0                �                               H     �                                   �                                                      0                      @                               '                    #DATA /              �$                              /     �                    #DISTRIBUTION_ .                                         y#DISTRIBUTION_ .                                                                @                           -     LEN_TRIM    �   X      fn#fn (   �   ?  b   uapp(CLASS_DISTRIBUTION    7  @   J  MPI_SIESTA $   w  V       gen@NEWDISTRIBUTION !   �  �      NEWDISTRIBUTION_ &   l  Z   a   NEWDISTRIBUTION_%THIS *   �  @   a   NEWDISTRIBUTION_%REF_COMM 3     �   a   NEWDISTRIBUTION_%RANKS_IN_REF_COMM +   �  @   a   NEWDISTRIBUTION_%DIST_TYPE +   �  @   a   NEWDISTRIBUTION_%BLOCKSIZE &     L   a   NEWDISTRIBUTION_%NAME    ^  O       gen@DISTTYPE    �  Z      DIST_TYPE      Z   a   DIST_TYPE%THIS    a  O       gen@REF_COMM    �  Z      REF_COMM_    
  Z   a   REF_COMM_%THIS *   d  \       gen@GET_RANKS_IN_REF_COMM '   �  ]      GET_RANKS_IN_REF_COMM_ ,     Z   a   GET_RANKS_IN_REF_COMM_%THIS -   w  �   a   GET_RANKS_IN_REF_COMM_%RANKS '   	  Y       gen@NUM_LOCAL_ELEMENTS $   \	  n      NUM_LOCAL_ELEMENTS_ )   �	  Z   a   NUM_LOCAL_ELEMENTS_%THIS )   $
  @   a   NUM_LOCAL_ELEMENTS_%NELS )   d
  @   a   NUM_LOCAL_ELEMENTS_%NODE *   �
  \       gen@INDEX_LOCAL_TO_GLOBAL '      l      INDEX_LOCAL_TO_GLOBAL_ ,   l  Z   a   INDEX_LOCAL_TO_GLOBAL_%THIS *   �  @   a   INDEX_LOCAL_TO_GLOBAL_%IL ,     @   a   INDEX_LOCAL_TO_GLOBAL_%NODE *   F  \       gen@INDEX_GLOBAL_TO_LOCAL '   �  l      INDEX_GLOBAL_TO_LOCAL_ ,     Z   a   INDEX_GLOBAL_TO_LOCAL_%THIS *   h  @   a   INDEX_GLOBAL_TO_LOCAL_%IG ,   �  @   a   INDEX_GLOBAL_TO_LOCAL_%NODE *   �  \       gen@NODE_HANDLING_ELEMENT '   D  b      NODE_HANDLING_ELEMENT_ ,   �  Z   a   NODE_HANDLING_ELEMENT_%THIS *      @   a   NODE_HANDLING_ELEMENT_%IG    @  M      i@|    �  ]      ASSIGN_    �  Z   a   ASSIGN_%THIS    D  Z   a   ASSIGN_%OTHER    �  K       gen@INIT    �  R      INIT_    ;  Z   a   INIT_%THIS    �  M       gen@DELETE    �  R      DELETE_    4  Z   a   DELETE_%THIS    �  O       gen@REFCOUNT    �  Z      REFCOUNT_    7  Z   a   REFCOUNT_%THIS    �  I       gen@ID    �  b      ID_    <  Z   a   ID_%THIS    �  K       gen@NAME    �  �      NAME_    �  Z   a   NAME_%THIS       R       gen@INITIALIZED    g  Z      INITIALIZED_ "   �  Z   a   INITIALIZED_%THIS      K       gen@SAME    f  f      SAME_    �  Z   a   SAME_%THIS1    &  Z   a   SAME_%THIS2 +   �  q       MPI_COMM_NULL+MPI__INCLUDE ,   �  q       MPI_GROUP_NULL+MPI__INCLUDE +   b  p       MPI_UNDEFINED+MPI__INCLUDE *   �  @       MPI_COMM_WORLD+MPI_SIESTA @     q       TRUE_MPI_COMM_WORLD+MPI__INCLUDE=MPI_COMM_WORLD    �  p       TYPE_NULL "   �  q       TYPE_BLOCK_CYCLIC    d  q       TYPE_PEXSI    �  �       DISTRIBUTION_ '   �  �   a   DISTRIBUTION_%REFCOUNT !   h  P   a   DISTRIBUTION_%ID #   �  �  a   DISTRIBUTION_%NAME (   u  H   a   DISTRIBUTION_%DIST_TYPE '   �  �   a   DISTRIBUTION_%REF_COMM $   b   �   a   DISTRIBUTION_%GROUP #   !  �   a   DISTRIBUTION_%NODE 0   �!  �   a   DISTRIBUTION_%RANKS_IN_REF_COMM $   ?"  �   a   DISTRIBUTION_%NODES &   �"  �   a   DISTRIBUTION_%NODE_IO (   �#  �   a   DISTRIBUTION_%BLOCKSIZE '   -$  �   a   DISTRIBUTION_%ISRCPROC    �$  Z       DISTRIBUTION "   ,%  �   a   DISTRIBUTION%DATA    &  A      NAME_%LEN_TRIM 