  l  ì   k820309    l          18.0        ©¤]                                                                                                          
       /home/lleblanc/src/psiesta-4.1-b4/Src/class_SpData2D.F90 CLASS_SSPDATA2D              SSPDATA2D gen@VAL gen@INIT_VAL gen@NROWS gen@NROWS_G gen@NNZS gen@N_COL gen@LIST_PTR gen@LIST_COL gen@PRINT_TYPE gen@SIZE gen@NEWSSPDATA2D gen@SPAR gen@DIST gen@SPAR_DIM gen@INIT gen@DELETE i@| gen@REFCOUNT gen@ID gen@NAME gen@SAME gen@INITIALIZED                      @                              
                            @                              
                            @                              
                                                             u #VALDATA2D    #VALDATA2D_IDX    #VALSPDATA    #VALSPDATA_IDX    (         D   @                                                             	    #THIS              &                   &                                                     
                                                     #SDATA2D    %         @    @                                              	       #THIS    #IDX1 	   #IDX2 
             
                                                     #SDATA2D              
                                  	                     
                                  
           (         D   @X                                                               	    #THIS              &                   &                                                     
  @                                                  #SSPDATA2D    %         @    @X                                                	       #THIS    #IDX1    #IDX2              
  @                                                  #SSPDATA2D              
  @                                                    
  @                                                                                               u #INITIALIZEDATA2D    #INITIALIZESPDATA    #         @     @                                               #THIS              
                                                     #SDATA2D    #         @     @X                                                 #THIS              
 @                                                   #SSPDATA2D                                                          u #NROWSSPARSITY    #NROWSSPDATA    %         @    @                                                    #THIS              
                                                     #SPARSITY    %         @    @X                                                       #THIS              
  @                                                  #SSPDATA2D                                                          u #NROWS_GSPARSITY    #NROWS_GSPDATA    %         @    @                                                    #THIS              
                                                     #SPARSITY    %         @    @X                                                       #THIS              
  @                                                  #SSPDATA2D                                                          u #NNZSSPARSITY    #NNZSSPDATA !   %         @    @                                                    #THIS               
                                                      #SPARSITY    %         @    @X                           !                           #THIS "             
  @                               "                   #SSPDATA2D                                                          u #N_COLSPARSITY #   #N_COLSPDATA %   #N_COLSPARSITYI '   (         D   @                          #                                       #THIS $             &                                                     
                                  $                   #SPARSITY    (         D   @X                            %                                       #THIS &             &                                                     
  @                               &                   #SSPDATA2D    %         H    @                          '                           #THIS (   #ROW )             
                                  (                   #SPARSITY              
                                  )                                                                 u #LIST_PTRSPARSITY *   #LIST_PTRSPDATA ,   #LIST_PTRSPARSITYI .   (         D   @                          *                                       #THIS +             &                                                     
                                  +                   #SPARSITY    (         D   @X                            ,                                       #THIS -             &                                                     
  @                               -                   #SSPDATA2D    %         H    @                          .                           #THIS /   #I 0             
                                  /                   #SPARSITY              
                                  0                                                                 u #LIST_COLSPARSITY 1   #LIST_COLSPDATA 3   #LIST_COLSPARSITYI 5   (         D   @                          1                   	                    #THIS 2             &                                                     
                                  2                   #SPARSITY    (         D   @X                            3                                       #THIS 4             &                                                     
  @                               4                   #SSPDATA2D    %         H    @                          5                           #THIS 6   #I 7             
                                  6                   #SPARSITY              
                                  7                                                                u #PRINTDATA2D 8   #PRINTSPARSITY :   #PRINTORBITALDISTRIBUTION <   #PRINTSPDATA ?   #         @     @                           8                    #THIS 9             
                                  9                   #SDATA2D    #         @     @                           :                    #SP ;             
                                  ;                   #SPARSITY    #         @     @                            <                    #THIS =             
                                  =                   #ORBITALDISTRIBUTION >   #         @     @X                             ?                    #THIS @             
  @                               @                   #SSPDATA2D                                                        u #SIZEDATA2D A   #SIZEDATA2D_DIM C   #SIZESPDATA F   %         @    @                           A                           #THIS B             
                                  B                   #SDATA2D    %         @    @                           C                           #THIS D   #DIM E             
                                  D                   #SDATA2D              
                                  E           %         @    @X                            F                           #THIS G   #DIM H             
  @                               G                   #SSPDATA2D              
 @                               H                                                                  u #NEWSPDATAFROMDATA I   #NEWSPDATAFROMDIMS O   #         @     @X                             I                    #SP J   #A K   #DIST L   #THIS M   #NAME N             
                                  J                   #SPARSITY              
                                  K                   #SDATA2D              
                                  L                   #ORBITALDISTRIBUTION >             
D @                               M                    #SSPDATA2D              
 @                             N                    1 #         @     @X                             O                    #SP P   #DIM Q   #DIST R   #THIS S   #NAME T   #SPARSITY_DIM U             
  @                               P                   #SPARSITY              
  @                               Q                     
                                  R                   #ORBITALDISTRIBUTION >             
D @                               S                    #SSPDATA2D              
 @                             T                    1           
 @                               U                                                                  u #SPARSPDATA V   &         @   @X                            V                           #THIS W   #SPARSITY              
                                  W                   #SSPDATA2D                                                           u #DISTSPDATA X   &         @   @X                            X                           #THIS Y   #ORBITALDISTRIBUTION >             
                                  Y                   #SSPDATA2D                                                           u #SPAR_DIMSPDATA Z   %         @    @X                            Z                           #THIS [             
  @                               [                   #SSPDATA2D                                                         u #INIT_ \   #INIT_ ^   #INIT_ `   #INIT_ b   #         @     @                            \                    #THIS ]             
                                 ]                    #SDATA2D    #         @     @                            ^                    #THIS _             
                                 _                    #SPARSITY    #         @     @                            `                    #THIS a             
                                 a                    #ORBITALDISTRIBUTION >   #         @     @X                            b                    #THIS c             
D @                               c                    #SSPDATA2D                                                         u #DELETE_ d   #DELETE_ f   #DELETE_ h   #DELETE_ j   #         @     @                           d                    #THIS e             
                                 e                    #SDATA2D    #         @     @                           f                    #THIS g             
                                 g                    #SPARSITY    #         @     @                           h                    #THIS i             
                                 i                    #ORBITALDISTRIBUTION >   #         @     @X                            j                    #THIS k             
D @                               k                    #SSPDATA2D                                                              |  #ASSIGN_ l   #ASSIGN_ o   #         @     @                           l                    #THIS m   #OTHER n             
                                 m                    #SDATA2D              
                                  n                   #SDATA2D    #         @     @X                             o                    #THIS p   #OTHER q             
D @                               p                    #SSPDATA2D              
  @                               q                   #SSPDATA2D                                                          u #REFCOUNT_ r   #REFCOUNT_ t   #REFCOUNT_ v   #REFCOUNT_ x   %         @    @                           r                           #THIS s             
                                  s                   #SDATA2D    %         @    @                           t                           #THIS u             
                                  u                   #SPARSITY    %         @    @                           v                           #THIS w             
                                  w                   #ORBITALDISTRIBUTION >   %         @    @X                           x                           #THIS y             
                                  y                   #SSPDATA2D                                                          u #ID_ z   #ID_ |   #ID_ ~   #ID_    $         @    @                           z     $                      #THIS {                     
                                  {                   #SDATA2D    $         @    @                           |     $                      #THIS }                     
                                  }                   #SPARSITY    $         @    @                           ~     $                      #THIS                      
                                                     #ORBITALDISTRIBUTION >   $         @    @X                                 $                      #THIS                      
                                                     #SSPDATA2D                                                          u #NAME_    #NAME_    #NAME_    #NAME_    $        @   @                                                       #THIS    I r      5 8 8#SDATA2D_     O#SDATA2D     p        U #SDATA2D_           U                    
                                                     #SDATA2D    $        @   @                                                       #THIS    I r      5 8 8#SPARSITY_     O#SPARSITY     p        U #SPARSITY_           U                    
                                                     #SPARSITY    $        @   @                                                       #THIS    I r      5 8 8#ORBITALDISTRIBUTION_     O#ORBITALDISTRIBUTION >    p        U #ORBITALDISTRIBUTION_     >      U                    
                                                     #ORBITALDISTRIBUTION >   $        @   @X                                                        #THIS    I r      5 8 8#SSPDATA2D_     O#SSPDATA2D     p        U #SSPDATA2D_           U                   
                                                     #SSPDATA2D                                                          u #SAME_    #SAME_    #SAME_     #SAME_ £   %         @    @                                                     #THIS1    #THIS2              
                                                     #SDATA2D              
                                                     #SDATA2D    %         @    @                                                     #THIS1    #THIS2              
                                                     #SPARSITY              
                                                     #SPARSITY    %         @    @                                                      #THIS1 ¡   #THIS2 ¢             
                                  ¡                   #ORBITALDISTRIBUTION >             
                                  ¢                   #ORBITALDISTRIBUTION >   %         @    @X                           £                           #THIS1 ¤   #THIS2 ¥             
  @                               ¤                   #SSPDATA2D              
  @                               ¥                   #SSPDATA2D                                                          u #INITIALIZED_ ¦   #INITIALIZED_ ¨   #INITIALIZED_ ª   #INITIALIZED_ ¬   %         @    @                          ¦                           #THIS §             
                                  §                   #SDATA2D    %         @    @                          ¨                           #THIS ©             
                                  ©                   #SPARSITY    %         @    @                          ª                           #THIS «             
                                  «                   #ORBITALDISTRIBUTION >   %         @    @X                          ¬                           #THIS ­             
                                  ­                   #SSPDATA2D                                                              |  #ASSIGN_ ®   #ASSIGN_ ±   #         @     @                           ®                    #THIS ¯   #OTHER °             
                                 ¯                    #SPARSITY              
                                  °                   #SPARSITY    #         @     @                           ±                    #THIS ²   #OTHER ³             
                                 ²                    #ORBITALDISTRIBUTION >             
                                  ³                   #ORBITALDISTRIBUTION >                     @                              '                    #DATA                $                                                       #SPARSITY_                                          y#SPARSITY_                                                                 @  @                              '                   #REFCOUNT ´   #ID µ   #NAME    #NROWS ¶   #NROWS_G ·   #NCOLS ¸   #NCOLS_G ¹   #NNZS º   #N_COL »   #LIST_COL ¼   #LIST_PTR ½                                              ´                                                                                                0                                              µ     $                                                                            %                     Cnull_id                                                                                                       (                                                                           ¢              Cnull_sparsity                                                                                                                                                                                                                                                                                                                  ¶     (                                                                                         0                                               ·     ,                                                                                         0                                               ¸     0                                                                                         0                                               ¹     4                                                                                         0                                               º     8                                                                                         0                                             »            @             	              &                                                                                 y                                                                                        ¼                         
              &                                                                                 y                                                                                        ½            Ð                           &                                                                                 y                                                                 @                              '                    #DATA                $                                                       #SDATA2D_                                          y#SDATA2D_                                                                 @  @                              '                   #REFCOUNT ¾   #ID ¿   #NAME    #VAL À                                              ¾                                                                                                0                                              ¿     $                                                                            %                     Cnull_id                                                                                                       (                                                                                         Cnull sData2D                                                                                                                                                                                                                                                                                                               À            (               	            &                   &                                                                                 y	                                                                 @                         >     '                    #DATA                $                                   h                    #ORBITALDISTRIBUTION_                                          y#ORBITALDISTRIBUTION_                                                                 @  @                              'h                   #REFCOUNT Á   #ID Â   #NAME    #COMM Ã   #NODE Ä   #NODES Å   #NODE_IO Æ   #BLOCKSIZE Ç   #ISRCPROC È   #N É   #NROC_PROC Ê   #NL2G Ë   #NG2L Ì   #NG2P Í                                              Á                                                                                                0                                               Â     $                                                                           (                                                                           Ê              Cnull OrbitalDistribution                                                                                                                                                                                                                                                                                                       Ã     (                                                                            ÿÿÿÿÿÿÿÿ                                                       Ä     ,                                                                            ÿÿÿÿÿÿÿÿ                                                       Å     0                                                                                         0                                               Æ     4                                                                            ÿÿÿÿÿÿÿÿ                                                       Ç     8                                                                                         0                                               È     <      	                                                                                   0                                               É     @      
                                                                      ÿÿÿÿÿÿÿÿ                                                     Ê            H                           &                                                                                 y                                                                                        Ë                                       &                                                                                 y                                                                                        Ì            Ø                           &                                                                                 y                                                                                        Í                                        &                                                                                 y                                                                                          Î     SIZE               @  @                              '@                   #REFCOUNT Ï   #ID Ð   #NAME    #SP Ñ   #A Ò   #DIST Ó                                              Ï                                                                                                0                                              Ð     $                                                                            %                       Cnull_id                                                                                                      (                                                                                           Cnull sSpData2D                                                                                                                                                                                                                                                                                                                  Ñ            (             #SPARSITY                                                Ò            0             #SDATA2D                                                Ó            8             #ORBITALDISTRIBUTION >                     @                               '                    #DATA               $                                   @                    #SSPDATA2D_                                          y#SSPDATA2D_                                                                 @                                LEN_TRIM               @                                LEN_TRIM               @                                LEN_TRIM               @                                LEN_TRIM        Q      fn#fn %   ñ     b   uapp(CLASS_SSPDATA2D    ù  @   J  CLASS_SDATA2D    9  @   J  CLASS_SPARSITY *   y  @   J  CLASS_ORBITALDISTRIBUTION    ¹         gen@VAL (   =  ¾      VALDATA2D+CLASS_SDATA2D -   û  U   a   VALDATA2D%THIS+CLASS_SDATA2D ,   P  n      VALDATA2D_IDX+CLASS_SDATA2D 1   ¾  U   a   VALDATA2D_IDX%THIS+CLASS_SDATA2D 1     @   a   VALDATA2D_IDX%IDX1+CLASS_SDATA2D 1   S  @   a   VALDATA2D_IDX%IDX2+CLASS_SDATA2D      ¾      VALSPDATA    Q  W   a   VALSPDATA%THIS    ¨  n      VALSPDATA_IDX #     W   a   VALSPDATA_IDX%THIS #   m  @   a   VALSPDATA_IDX%IDX1 #   ­  @   a   VALSPDATA_IDX%IDX2    í  l       gen@INIT_VAL /   Y  R      INITIALIZEDATA2D+CLASS_SDATA2D 4   «  U   a   INITIALIZEDATA2D%THIS+CLASS_SDATA2D !    	  R      INITIALIZESPDATA &   R	  W   a   INITIALIZESPDATA%THIS    ©	  d       gen@NROWS -   
  Z      NROWSSPARSITY+CLASS_SPARSITY 2   g
  V   a   NROWSSPARSITY%THIS+CLASS_SPARSITY    ½
  Z      NROWSSPDATA !     W   a   NROWSSPDATA%THIS    n  h       gen@NROWS_G /   Ö  Z      NROWS_GSPARSITY+CLASS_SPARSITY 4   0  V   a   NROWS_GSPARSITY%THIS+CLASS_SPARSITY      Z      NROWS_GSPDATA #   à  W   a   NROWS_GSPDATA%THIS    7  b       gen@NNZS ,     Z      NNZSSPARSITY+CLASS_SPARSITY 1   ó  V   a   NNZSSPARSITY%THIS+CLASS_SPARSITY    I  Z      NNZSSPDATA     £  W   a   NNZSSPDATA%THIS    ú  x       gen@N_COL -   r  ¦      N_COLSPARSITY+CLASS_SPARSITY 2     V   a   N_COLSPARSITY%THIS+CLASS_SPARSITY    n  ¦      N_COLSPDATA !     W   a   N_COLSPDATA%THIS .   k  c      N_COLSPARSITYI+CLASS_SPARSITY 3   Î  V   a   N_COLSPARSITYI%THIS+CLASS_SPARSITY 2   $  @   a   N_COLSPARSITYI%ROW+CLASS_SPARSITY    d         gen@LIST_PTR 0   å  ¦      LIST_PTRSPARSITY+CLASS_SPARSITY 5     V   a   LIST_PTRSPARSITY%THIS+CLASS_SPARSITY    á  ¦      LIST_PTRSPDATA $     W   a   LIST_PTRSPDATA%THIS 1   Þ  a      LIST_PTRSPARSITYI+CLASS_SPARSITY 6   ?  V   a   LIST_PTRSPARSITYI%THIS+CLASS_SPARSITY 3     @   a   LIST_PTRSPARSITYI%I+CLASS_SPARSITY    Õ         gen@LIST_COL 0   V  ¦      LIST_COLSPARSITY+CLASS_SPARSITY 5   ü  V   a   LIST_COLSPARSITY%THIS+CLASS_SPARSITY    R  ¦      LIST_COLSPDATA $   ø  W   a   LIST_COLSPDATA%THIS 1   O  a      LIST_COLSPARSITYI+CLASS_SPARSITY 6   °  V   a   LIST_COLSPARSITYI%THIS+CLASS_SPARSITY 3     @   a   LIST_COLSPARSITYI%I+CLASS_SPARSITY    F         gen@PRINT_TYPE *   Ù  R      PRINTDATA2D+CLASS_SDATA2D /   +  U   a   PRINTDATA2D%THIS+CLASS_SDATA2D -     P      PRINTSPARSITY+CLASS_SPARSITY 0   Ð  V   a   PRINTSPARSITY%SP+CLASS_SPARSITY C   &  R      PRINTORBITALDISTRIBUTION+CLASS_ORBITALDISTRIBUTION H   x  a   a   PRINTORBITALDISTRIBUTION%THIS+CLASS_ORBITALDISTRIBUTION    Ù  R      PRINTSPDATA !   +  W   a   PRINTSPDATA%THIS      t       gen@SIZE )   ö  Z      SIZEDATA2D+CLASS_SDATA2D .   P  U   a   SIZEDATA2D%THIS+CLASS_SDATA2D -   ¥  c      SIZEDATA2D_DIM+CLASS_SDATA2D 2     U   a   SIZEDATA2D_DIM%THIS+CLASS_SDATA2D 1   ]  @   a   SIZEDATA2D_DIM%DIM+CLASS_SDATA2D      c      SIZESPDATA        W   a   SIZESPDATA%THIS    W  @   a   SIZESPDATA%DIM !     n       gen@NEWSSPDATA2D "      u      NEWSPDATAFROMDATA %   z   V   a   NEWSPDATAFROMDATA%SP $   Ð   U   a   NEWSPDATAFROMDATA%A '   %!  a   a   NEWSPDATAFROMDATA%DIST '   !  W   a   NEWSPDATAFROMDATA%THIS '   Ý!  L   a   NEWSPDATAFROMDATA%NAME "   )"        NEWSPDATAFROMDIMS %   ²"  V   a   NEWSPDATAFROMDIMS%SP &   #  @   a   NEWSPDATAFROMDIMS%DIM '   H#  a   a   NEWSPDATAFROMDIMS%DIST '   ©#  W   a   NEWSPDATAFROMDIMS%THIS '    $  L   a   NEWSPDATAFROMDIMS%NAME /   L$  @   a   NEWSPDATAFROMDIMS%SPARSITY_DIM    $  P       gen@SPAR    Ü$  h      SPARSPDATA     D%  W   a   SPARSPDATA%THIS    %  P       gen@DIST    ë%  s      DISTSPDATA     ^&  W   a   DISTSPDATA%THIS    µ&  T       gen@SPAR_DIM    	'  Z      SPAR_DIMSPDATA $   c'  W   a   SPAR_DIMSPDATA%THIS    º'  l       gen@INIT $   &(  R      INIT_+CLASS_SDATA2D )   x(  U   a   INIT_%THIS+CLASS_SDATA2D %   Í(  R      INIT_+CLASS_SPARSITY *   )  V   a   INIT_%THIS+CLASS_SPARSITY 0   u)  R      INIT_+CLASS_ORBITALDISTRIBUTION 5   Ç)  a   a   INIT_%THIS+CLASS_ORBITALDISTRIBUTION    (*  R      INIT_    z*  W   a   INIT_%THIS    Ñ*  t       gen@DELETE &   E+  R      DELETE_+CLASS_SDATA2D +   +  U   a   DELETE_%THIS+CLASS_SDATA2D '   ì+  R      DELETE_+CLASS_SPARSITY ,   >,  V   a   DELETE_%THIS+CLASS_SPARSITY 2   ,  R      DELETE_+CLASS_ORBITALDISTRIBUTION 7   æ,  a   a   DELETE_%THIS+CLASS_ORBITALDISTRIBUTION    G-  R      DELETE_    -  W   a   DELETE_%THIS    ð-  Z      i@| &   J.  ]      ASSIGN_+CLASS_SDATA2D +   §.  U   a   ASSIGN_%THIS+CLASS_SDATA2D ,   ü.  U   a   ASSIGN_%OTHER+CLASS_SDATA2D    Q/  ]      ASSIGN_    ®/  W   a   ASSIGN_%THIS    0  W   a   ASSIGN_%OTHER    \0  |       gen@REFCOUNT (   Ø0  Z      REFCOUNT_+CLASS_SDATA2D -   21  U   a   REFCOUNT_%THIS+CLASS_SDATA2D )   1  Z      REFCOUNT_+CLASS_SPARSITY .   á1  V   a   REFCOUNT_%THIS+CLASS_SPARSITY 4   72  Z      REFCOUNT_+CLASS_ORBITALDISTRIBUTION 9   2  a   a   REFCOUNT_%THIS+CLASS_ORBITALDISTRIBUTION    ò2  Z      REFCOUNT_    L3  W   a   REFCOUNT_%THIS    £3  d       gen@ID "   4  b      ID_+CLASS_SDATA2D '   i4  U   a   ID_%THIS+CLASS_SDATA2D #   ¾4  b      ID_+CLASS_SPARSITY (    5  V   a   ID_%THIS+CLASS_SPARSITY .   v5  b      ID_+CLASS_ORBITALDISTRIBUTION 3   Ø5  a   a   ID_%THIS+CLASS_ORBITALDISTRIBUTION    96  b      ID_    6  W   a   ID_%THIS    ò6  l       gen@NAME $   ^7  Ë      NAME_+CLASS_SDATA2D )   )8  U   a   NAME_%THIS+CLASS_SDATA2D %   ~8  Î      NAME_+CLASS_SPARSITY *   L9  V   a   NAME_%THIS+CLASS_SPARSITY 0   ¢9  ï      NAME_+CLASS_ORBITALDISTRIBUTION 5   :  a   a   NAME_%THIS+CLASS_ORBITALDISTRIBUTION    ò:  Ñ      NAME_    Ã;  W   a   NAME_%THIS    <  l       gen@SAME $   <  f      SAME_+CLASS_SDATA2D *   ì<  U   a   SAME_%THIS1+CLASS_SDATA2D *   A=  U   a   SAME_%THIS2+CLASS_SDATA2D %   =  f      SAME_+CLASS_SPARSITY +   ü=  V   a   SAME_%THIS1+CLASS_SPARSITY +   R>  V   a   SAME_%THIS2+CLASS_SPARSITY 0   ¨>  f      SAME_+CLASS_ORBITALDISTRIBUTION 6   ?  a   a   SAME_%THIS1+CLASS_ORBITALDISTRIBUTION 6   o?  a   a   SAME_%THIS2+CLASS_ORBITALDISTRIBUTION    Ð?  f      SAME_    6@  W   a   SAME_%THIS1    @  W   a   SAME_%THIS2     ä@         gen@INITIALIZED +   lA  Z      INITIALIZED_+CLASS_SDATA2D 0   ÆA  U   a   INITIALIZED_%THIS+CLASS_SDATA2D ,   B  Z      INITIALIZED_+CLASS_SPARSITY 1   uB  V   a   INITIALIZED_%THIS+CLASS_SPARSITY 7   ËB  Z      INITIALIZED_+CLASS_ORBITALDISTRIBUTION <   %C  a   a   INITIALIZED_%THIS+CLASS_ORBITALDISTRIBUTION    C  Z      INITIALIZED_ "   àC  W   a   INITIALIZED_%THIS    7D  Z      i@| '   D  ]      ASSIGN_+CLASS_SPARSITY ,   îD  V   a   ASSIGN_%THIS+CLASS_SPARSITY -   DE  V   a   ASSIGN_%OTHER+CLASS_SPARSITY 2   E  ]      ASSIGN_+CLASS_ORBITALDISTRIBUTION 7   ÷E  a   a   ASSIGN_%THIS+CLASS_ORBITALDISTRIBUTION 8   XF  a   a   ASSIGN_%OTHER+CLASS_ORBITALDISTRIBUTION (   ¹F  Z       SPARSITY+CLASS_SPARSITY -   G  Î   a   SPARSITY%DATA+CLASS_SPARSITY )   áG  Ñ      SPARSITY_+CLASS_SPARSITY 2   ²H  ¥   a   SPARSITY_%REFCOUNT+CLASS_SPARSITY ,   WI  á   a   SPARSITY_%ID+CLASS_SPARSITY .   8J  ½  a   SPARSITY_%NAME+CLASS_SPARSITY /   õK  ¥   a   SPARSITY_%NROWS+CLASS_SPARSITY 1   L  ¥   a   SPARSITY_%NROWS_G+CLASS_SPARSITY /   ?M  ¥   a   SPARSITY_%NCOLS+CLASS_SPARSITY 1   äM  ¥   a   SPARSITY_%NCOLS_G+CLASS_SPARSITY .   N  ¥   a   SPARSITY_%NNZS+CLASS_SPARSITY /   .O  ô   a   SPARSITY_%N_COL+CLASS_SPARSITY 2   "P  ô   a   SPARSITY_%LIST_COL+CLASS_SPARSITY 2   Q  ô   a   SPARSITY_%LIST_PTR+CLASS_SPARSITY &   
R  Z       SDATA2D+CLASS_SDATA2D +   dR  Ì   a   SDATA2D%DATA+CLASS_SDATA2D '   0S  y      SDATA2D_+CLASS_SDATA2D 0   ©S  ¥   a   SDATA2D_%REFCOUNT+CLASS_SDATA2D *   NT  á   a   SDATA2D_%ID+CLASS_SDATA2D ,   /U  ½  a   SDATA2D_%NAME+CLASS_SDATA2D +   ìV    a   SDATA2D_%VAL+CLASS_SDATA2D >   øW  Z       ORBITALDISTRIBUTION+CLASS_ORBITALDISTRIBUTION C   RX  ä   a   ORBITALDISTRIBUTION%DATA+CLASS_ORBITALDISTRIBUTION ?   6Y  í      ORBITALDISTRIBUTION_+CLASS_ORBITALDISTRIBUTION H   #Z  ¥   a   ORBITALDISTRIBUTION_%REFCOUNT+CLASS_ORBITALDISTRIBUTION B   ÈZ  P   a   ORBITALDISTRIBUTION_%ID+CLASS_ORBITALDISTRIBUTION D   [  ½  a   ORBITALDISTRIBUTION_%NAME+CLASS_ORBITALDISTRIBUTION D   Õ\  ¤   a   ORBITALDISTRIBUTION_%COMM+CLASS_ORBITALDISTRIBUTION D   y]  ¤   a   ORBITALDISTRIBUTION_%NODE+CLASS_ORBITALDISTRIBUTION E   ^  ¥   a   ORBITALDISTRIBUTION_%NODES+CLASS_ORBITALDISTRIBUTION G   Â^  ¤   a   ORBITALDISTRIBUTION_%NODE_IO+CLASS_ORBITALDISTRIBUTION I   f_  ¥   a   ORBITALDISTRIBUTION_%BLOCKSIZE+CLASS_ORBITALDISTRIBUTION H   `  ¥   a   ORBITALDISTRIBUTION_%ISRCPROC+CLASS_ORBITALDISTRIBUTION A   °`  ¤   a   ORBITALDISTRIBUTION_%N+CLASS_ORBITALDISTRIBUTION I   Ta  ô   a   ORBITALDISTRIBUTION_%NROC_PROC+CLASS_ORBITALDISTRIBUTION D   Hb  ô   a   ORBITALDISTRIBUTION_%NL2G+CLASS_ORBITALDISTRIBUTION D   <c  ô   a   ORBITALDISTRIBUTION_%NG2L+CLASS_ORBITALDISTRIBUTION D   0d  ô   a   ORBITALDISTRIBUTION_%NG2P+CLASS_ORBITALDISTRIBUTION #   $e  =       SIZE+CLASS_SDATA2D    ae         SSPDATA2D_ $   êe  ¥   a   SSPDATA2D_%REFCOUNT    f  á   a   SSPDATA2D_%ID     pg  ½  a   SSPDATA2D_%NAME    -i  ^   a   SSPDATA2D_%SP    i  ]   a   SSPDATA2D_%A     èi  i   a   SSPDATA2D_%DIST    Qj  Z       SSPDATA2D    «j  Ð   a   SSPDATA2D%DATA    {k  A      NAME_%LEN_TRIM B   ¼k  A      NAME_%LEN_TRIM+CLASS_ORBITALDISTRIBUTION=LEN_TRIM 7   ýk  A      NAME_%LEN_TRIM+CLASS_SPARSITY=LEN_TRIM 6   >l  A      NAME_%LEN_TRIM+CLASS_SDATA2D=LEN_TRIM 