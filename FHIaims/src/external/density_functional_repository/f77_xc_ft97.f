      subroutine uks_xc_ft97
     & (ideriv,npt,rhoa1,rhob1,sigmaaa1,sigmabb1,sigmaab1,
     &  zk,vrhoa,vrhob,vsigmaaa,vsigmabb,vsigmaab,
     &  v2rhoa2,v2rhob2,v2rhoab,
     &  v2rhoasigmaaa,v2rhoasigmaab,v2rhoasigmabb,
     &  v2rhobsigmabb,v2rhobsigmaab,v2rhobsigmaaa,
     &  v2sigmaaa2,v2sigmaaaab,v2sigmaaabb,
     &  v2sigmaab2,v2sigmaabbb,v2sigmabb2)
c
c     M. Filatov, and W. Thiel
c     A nonlocal correlation energy density functional from a Coulomb
c     hole model
c     Int. J. Quant. Chem. 62 (1997) 603-616
c
c     M. Filatov, and W. Thiel
c     A new gradient-corrected exchange-correlation density functional
c     Mol. Phys. 91 (1997) 847-859
c
c
c     CITATION:
c
c     Functionals were obtained from the Density Functional Repository 
c     as developed and distributed by the Quantum Chemistry Group, 
c     CCLRC Daresbury Laboratory, Daresbury, Cheshire, WA4 4AD 
c     United Kingdom. Contact Huub van Dam (h.j.j.vandam@dl.ac.uk) or 
c     Paul Sherwood for further information.
c
c     COPYRIGHT:
c
c     Users may incorporate the source code into software packages and
c     redistribute the source code provided the source code is not
c     changed in anyway and is properly cited in any documentation or
c     publication related to its use.
c
c     ACKNOWLEDGEMENT:
c
c     The source code was generated using Maple 8 through a modified
c     version of the dfauto script published in:
c
c        R. Strange, F.R. Manby, P.J. Knowles
c        Automatic code generation in density functional theory
c        Comp. Phys. Comm. 136 (2001) 310-318.
c
      implicit none
      integer ideriv,npt
      real*8 rhoa1(npt),rhob1(npt)
      real*8 sigmaaa1(npt),sigmabb1(npt),sigmaab1(npt)
      real*8 zk(npt),vrhoa(npt),vrhob(npt)
      real*8 vsigmaaa(npt),vsigmabb(npt),vsigmaab(npt)
      real*8 v2rhoa2(npt),v2rhob2(npt),v2rhoab(npt)
      real*8 v2rhoasigmaaa(npt),v2rhoasigmaab(npt)
      real*8 v2rhoasigmabb(npt),v2rhobsigmabb(npt)
      real*8 v2rhobsigmaab(npt),v2rhobsigmaaa(npt)
      real*8 v2sigmaaa2(npt),v2sigmaaaab(npt),v2sigmaaabb(npt)
      real*8 v2sigmaab2(npt),v2sigmaabbb(npt),v2sigmabb2(npt)
c
      integer mpt,i
      parameter(mpt=200)
      real*8 tk(mpt),trhoa(mpt),trhob(mpt)
      real*8 tsigmaaa(mpt),tsigmabb(mpt),tsigmaab(mpt)
      real*8 t2rhoa2(mpt),t2rhob2(mpt),t2rhoab(mpt)
      real*8 t2rhoasigmaaa(mpt),t2rhoasigmaab(mpt)
      real*8 t2rhoasigmabb(mpt),t2rhobsigmabb(mpt)
      real*8 t2rhobsigmaab(mpt),t2rhobsigmaaa(mpt)
      real*8 t2sigmaaa2(mpt),t2sigmaaaab(mpt),t2sigmaaabb(mpt)
      real*8 t2sigmaab2(mpt),t2sigmaabbb(mpt),t2sigmabb2(mpt)
c
      if (npt.gt.mpt) then
         write(*,*)'*** ERROR: npt.gt.mpt re-parametrize'
      endif
c
c     Clearly there is a problem here. The above mucking about
c     with temporary arrays is a terrible hack. Please adjust this
c     in a way suitable to your own code.
c
      call uks_x_ft97b
     & (ideriv,npt,rhoa1,rhob1,sigmaaa1,sigmabb1,sigmaab1,
     &  zk,vrhoa,vrhob,vsigmaaa,vsigmabb,vsigmaab,
     &  v2rhoa2,v2rhob2,v2rhoab,
     &  v2rhoasigmaaa,v2rhoasigmaab,v2rhoasigmabb,
     &  v2rhobsigmabb,v2rhobsigmaab,v2rhobsigmaaa,
     &  v2sigmaaa2,v2sigmaaaab,v2sigmaaabb,
     &  v2sigmaab2,v2sigmaabbb,v2sigmabb2)
      call uks_c_ft97
     & (ideriv,npt,rhoa1,rhob1,sigmaaa1,sigmabb1,sigmaab1,
     &  tk,trhoa,trhob,tsigmaaa,tsigmabb,tsigmaab,
     &  t2rhoa2,t2rhob2,t2rhoab,
     &  t2rhoasigmaaa,t2rhoasigmaab,t2rhoasigmabb,
     &  t2rhobsigmabb,t2rhobsigmaab,t2rhobsigmaaa,
     &  t2sigmaaa2,t2sigmaaaab,t2sigmaaabb,
     &  t2sigmaab2,t2sigmaabbb,t2sigmabb2)
c
      do i=1,npt
         zk(i)=zk(i)+tk(i)
         if (ideriv.ge.1) then
            vrhoa(i)=vrhoa(i)+trhoa(i)
            vrhob(i)=vrhob(i)+trhob(i)
            vsigmaaa(i)=vsigmaaa(i)+tsigmaaa(i)
            vsigmaab(i)=vsigmaab(i)+tsigmaab(i)
            vsigmabb(i)=vsigmabb(i)+tsigmabb(i)
         endif
         if (ideriv.ge.2) then
            v2rhoa2(i)=v2rhoa2(i)+t2rhoa2(i)
            v2rhob2(i)=v2rhob2(i)+t2rhob2(i)
            v2rhoab(i)=v2rhoab(i)+t2rhoab(i)
            v2rhoasigmaaa(i)=v2rhoasigmaaa(i)+t2rhoasigmaaa(i)
            v2rhoasigmaab(i)=v2rhoasigmaab(i)+t2rhoasigmaab(i)
            v2rhoasigmabb(i)=v2rhoasigmabb(i)+t2rhoasigmabb(i)
            v2rhobsigmaaa(i)=v2rhobsigmaaa(i)+t2rhobsigmaaa(i)
            v2rhobsigmaab(i)=v2rhobsigmaab(i)+t2rhobsigmaab(i)
            v2rhobsigmabb(i)=v2rhobsigmabb(i)+t2rhobsigmabb(i)
            v2sigmaaa2(i)=v2sigmaaa2(i)+t2sigmaaa2(i)
            v2sigmaab2(i)=v2sigmaab2(i)+t2sigmaab2(i)
            v2sigmabb2(i)=v2sigmabb2(i)+t2sigmabb2(i)
            v2sigmaaaab(i)=v2sigmaaaab(i)+t2sigmaaaab(i)
            v2sigmaaabb(i)=v2sigmaaabb(i)+t2sigmaaabb(i)
            v2sigmaabbb(i)=v2sigmaabbb(i)+t2sigmaabbb(i)
         endif
      enddo
c
      return
      end
c
      subroutine rks_xc_ft97
     & (ideriv,npt,rhoa1,sigmaaa1,
     &  zk,vrhoa,vsigmaaa,
     &  v2rhoa2,v2rhoasigmaaa,v2sigmaaa2)
c
c     M. Filatov, and W. Thiel
c     A nonlocal correlation energy density functional from a
c     Coulomb hole model
c     Int. J. Quant. Chem. 62 (1997) 603-616
c
c     M. Filatov, and W. Thiel
c     A new gradient-corrected exchange-correlation density functional
c     Mol. Phys. 91 (1997) 847-859
c
c
c     CITATION:
c
c     Functionals were obtained from the Density Functional Repository 
c     as developed and distributed by the Quantum Chemistry Group, 
c     CCLRC Daresbury Laboratory, Daresbury, Cheshire, WA4 4AD 
c     United Kingdom. Contact Huub van Dam (h.j.j.vandam@dl.ac.uk) or 
c     Paul Sherwood for further information.
c
c     COPYRIGHT:
c
c     Users may incorporate the source code into software packages and
c     redistribute the source code provided the source code is not
c     changed in anyway and is properly cited in any documentation or
c     publication related to its use.
c
c     ACKNOWLEDGEMENT:
c
c     The source code was generated using Maple 8 through a modified
c     version of the dfauto script published in:
c
c        R. Strange, F.R. Manby, P.J. Knowles
c        Automatic code generation in density functional theory
c        Comp. Phys. Comm. 136 (2001) 310-318.
c
      implicit none
      integer ideriv,npt
      real*8 rhoa1(npt)
      real*8 sigmaaa1(npt)
      real*8 zk(npt),vrhoa(npt),vsigmaaa(npt)
      real*8 v2rhoa2(npt),v2rhoasigmaaa(npt),v2sigmaaa2(npt)
c
      integer mpt,i
      parameter(mpt=200)
      real*8 tk(mpt),trhoa(mpt),tsigmaaa(mpt)
      real*8 t2rhoa2(mpt),t2rhoasigmaaa(mpt),t2sigmaaa2(mpt)
c
      if (npt.gt.mpt) then
         write(*,*)'*** ERROR: npt.gt.mpt re-parametrize'
      endif
c
c     Clearly there is a problem here. The above mucking about
c     with temporary arrays is a terrible hack. Please adjust this
c     in a way suitable to your own code.
c
      call rks_x_ft97b
     & (ideriv,npt,rhoa1,sigmaaa1,
     &  zk,vrhoa,vsigmaaa,
     &  v2rhoa2,v2rhoasigmaaa,v2sigmaaa2)
      call rks_c_ft97
     & (ideriv,npt,rhoa1,sigmaaa1,
     &  tk,trhoa,tsigmaaa,
     &  t2rhoa2,t2rhoasigmaaa,t2sigmaaa2)
c
      do i=1,npt
         zk(i)=zk(i)+tk(i)
         if (ideriv.ge.1) then
            vrhoa(i)=vrhoa(i)+trhoa(i)
            vsigmaaa(i)=vsigmaaa(i)+tsigmaaa(i)
         endif
         if (ideriv.ge.2) then
            v2rhoa2(i)=v2rhoa2(i)+t2rhoa2(i)
            v2rhoasigmaaa(i)=v2rhoasigmaaa(i)+t2rhoasigmaaa(i)
            v2sigmaaa2(i)=v2sigmaaa2(i)+t2sigmaaa2(i)
         endif
      enddo
c
      return
      end
c
c----------------------------------------------------------------------
c
      SUBROUTINE CALCEI(ARG,RESULT,INT)
C----------------------------------------------------------------------
C
C This Fortran 77 packet computes the exponential integrals Ei(x),
C  E1(x), and  exp(-x)*Ei(x)  for real arguments  x  where
C
C           integral (from t=-infinity to t=x) (exp(t)/t),  x > 0,
C  Ei(x) =
C          -integral (from t=-x to t=infinity) (exp(t)/t),  x < 0,
C
C  and where the first integral is a principal value integral.
C  The packet contains three function type subprograms: EI, EONE,
C  and EXPEI;  and one subroutine type subprogram: CALCEI.  The
C  calling statements for the primary entries are
C
C                 Y = EI(X),            where  X .NE. 0,
C
C                 Y = EONE(X),          where  X .GT. 0,
C  and
C                 Y = EXPEI(X),         where  X .NE. 0,
C
C  and where the entry points correspond to the functions Ei(x),
C  E1(x), and exp(-x)*Ei(x), respectively.  The routine CALCEI
C  is intended for internal packet use only, all computations within
C  the packet being concentrated in this routine.  The function
C  subprograms invoke CALCEI with the Fortran statement
C         CALL CALCEI(ARG,RESULT,INT)
C  where the parameter usage is as follows
C
C     Function                  Parameters for CALCEI
C       Call                 ARG             RESULT         INT
C
C      EI(X)              X .NE. 0          Ei(X)            1
C      EONE(X)            X .GT. 0         -Ei(-X)           2
C      EXPEI(X)           X .NE. 0          exp(-X)*Ei(X)    3
C
C  The main computation involves evaluation of rational Chebyshev
C  approximations published in Math. Comp. 22, 641-649 (1968), and
C  Math. Comp. 23, 289-303 (1969) by Cody and Thacher.  This
C  transportable program is patterned after the machine-dependent
C  FUNPACK packet  NATSEI,  but cannot match that version for
C  efficiency or accuracy.  This version uses rational functions
C  that theoretically approximate the exponential integrals to
C  at least 18 significant decimal digits.  The accuracy achieved
C  depends on the arithmetic system, the compiler, the intrinsic
C  functions, and proper selection of the machine-dependent
C  constants.
C
C
C*******************************************************************
C*******************************************************************
C
C Explanation of machine-dependent constants
C
C   beta = radix for the floating-point system.
C   minexp = smallest representable power of beta.
C   maxexp = smallest power of beta that overflows.
C   XBIG = largest argument acceptable to EONE; solution to
C          equation:
C                     exp(-x)/x * (1 + 1/x) = beta ** minexp.
C   XINF = largest positive machine number; approximately
C                     beta ** maxexp
C   XMAX = largest argument acceptable to EI; solution to
C          equation:  exp(x)/x * (1 + 1/x) = beta ** maxexp.
C
C     Approximate values for some important machines are:
C
C                           beta      minexp      maxexp
C
C  CRAY-1        (S.P.)       2       -8193        8191
C  Cyber 180/185 
C    under NOS   (S.P.)       2        -975        1070
C  IEEE (IBM/XT,
C    SUN, etc.)  (S.P.)       2        -126         128
C  IEEE (IBM/XT,
C    SUN, etc.)  (D.P.)       2       -1022        1024
C  IBM 3033      (D.P.)      16         -65          63
C  VAX D-Format  (D.P.)       2        -128         127
C  VAX G-Format  (D.P.)       2       -1024        1023
C
C                           XBIG       XINF       XMAX
C
C  CRAY-1        (S.P.)    5670.31  5.45E+2465   5686.21
C  Cyber 180/185 
C    under NOS   (S.P.)     669.31  1.26E+322     748.28
C  IEEE (IBM/XT,
C    SUN, etc.)  (S.P.)      82.93  3.40E+38       93.24
C  IEEE (IBM/XT,
C    SUN, etc.)  (D.P.)     701.84  1.79D+308     716.35
C  IBM 3033      (D.P.)     175.05  7.23D+75      179.85
C  VAX D-Format  (D.P.)      84.30  1.70D+38       92.54
C  VAX G-Format  (D.P.)     703.22  8.98D+307     715.66
C
C*******************************************************************
C*******************************************************************
C
C error returns
C
C  The following table shows the types of error that may be
C  encountered in this routine and the function value supplied
C  in each case.
C
C       error       Argument         Function values for
C                    Range         EI      EXPEI     EONE
C
C     UNDERFLOW  (-)X .GT. XBIG     0        -         0
C     OVERFLOW      X .GE. XMAX    XINF      -         -
C     ILLEGAL X       X = 0       -XINF    -XINF     XINF
C     ILLEGAL X      X .LT. 0       -        -     USE ABS(X)
C
C Intrinsic functions required are:
C
C     ABS, SQRT, EXP
C
C
C  Author: W. J. Cody
C          Mathematics abd Computer Science Division
C          Argonne National Laboratory
C          Argonne, IL 60439
C
C  Latest modification: September 9, 1988
C
C----------------------------------------------------------------------
      INTEGER I,INT
CS    real*8
CD    DOUBLE PRECISION 
      real*8
     1       A,ARG,B,C,D,EXP40,E,EI,F,FOUR,FOURTY,FRAC,HALF,ONE,P,
     2       PLG,PX,P037,P1,P2,q,QLG,QX,Q1,Q2,R,RESULT,S,SIX,SUMP,
     3       SUMQ,T,THREE,TWELVE,TWO,TWO4,W,X,XBIG,XINF,XMAX,XMX0,
     4       X0,X01,X02,X11,Y,YSQ,ZERO
      DIMENSION  A(7),B(6),C(9),D(9),E(10),F(10),P(10),q(10),R(10),
     1   S(9),P1(10),Q1(9),P2(10),Q2(9),PLG(4),QLG(4),PX(10),QX(10)
C----------------------------------------------------------------------
C  Mathematical constants
C   EXP40 = exp(40)
C   X0 = zero of Ei
C   X01/X11 + X02 = zero of Ei to extra precision
C----------------------------------------------------------------------
CS    DATA ZERO,P037,HALF,ONE,TWO/0.0E0,0.037E0,0.5E0,1.0E0,2.0E0/,
CS   1     THREE,FOUR,SIX,TWELVE,TWO4/3.0E0,4.0E0,6.0E0,12.E0,24.0E0/,
CS   2     FOURTY,EXP40/40.0E0,2.3538526683701998541E17/,
CS   3     X01,X11,X02/381.5E0,1024.0E0,-5.1182968633365538008E-5/,
CS   4     X0/3.7250741078136663466E-1/
      DATA ZERO,P037,HALF,ONE,TWO/0.0D0,0.037D0,0.5D0,1.0D0,2.0D0/,
     1     THREE,FOUR,SIX,TWELVE,TWO4/3.0D0,4.0D0,6.0D0,12.D0,24.0D0/,
     2     FOURTY,EXP40/40.0D0,2.3538526683701998541D17/,
     3     X01,X11,X02/381.5D0,1024.0D0,-5.1182968633365538008D-5/,
     4     X0/3.7250741078136663466D-1/
C----------------------------------------------------------------------
C Machine-dependent constants
C----------------------------------------------------------------------
CS    DATA XINF/3.40E+38/,XMAX/93.246E0/,XBIG/82.93E0/
      DATA XINF/1.79D+308/,XMAX/716.351D0/,XBIG/701.84D0/
C----------------------------------------------------------------------
C Coefficients  for -1.0 <= X < 0.0
C----------------------------------------------------------------------
CS    DATA A/1.1669552669734461083368E2, 2.1500672908092918123209E3,
CS   1       1.5924175980637303639884E4, 8.9904972007457256553251E4,
CS   2       1.5026059476436982420737E5,-1.4815102102575750838086E5,
CS   3       5.0196785185439843791020E0/
CS    DATA B/4.0205465640027706061433E1, 7.5043163907103936624165E2,
CS   1       8.1258035174768735759855E3, 5.2440529172056355429883E4,
CS   2       1.8434070063353677359298E5, 2.5666493484897117319268E5/
      DATA A/1.1669552669734461083368D2, 2.1500672908092918123209D3,
     1       1.5924175980637303639884D4, 8.9904972007457256553251D4,
     2       1.5026059476436982420737D5,-1.4815102102575750838086D5,
     3       5.0196785185439843791020D0/
      DATA B/4.0205465640027706061433D1, 7.5043163907103936624165D2,
     1       8.1258035174768735759855D3, 5.2440529172056355429883D4,
     2       1.8434070063353677359298D5, 2.5666493484897117319268D5/
C----------------------------------------------------------------------
C Coefficients for -4.0 <= X < -1.0
C----------------------------------------------------------------------
CS    DATA C/3.828573121022477169108E-1, 1.107326627786831743809E+1,
CS   1       7.246689782858597021199E+1, 1.700632978311516129328E+2,
CS   2       1.698106763764238382705E+2, 7.633628843705946890896E+1,
CS   3       1.487967702840464066613E+1, 9.999989642347613068437E-1,
CS   4       1.737331760720576030932E-8/
CS    DATA D/8.258160008564488034698E-2, 4.344836335509282083360E+0,
CS   1       4.662179610356861756812E+1, 1.775728186717289799677E+2,
CS   2       2.953136335677908517423E+2, 2.342573504717625153053E+2,
CS   3       9.021658450529372642314E+1, 1.587964570758947927903E+1,
CS   4       1.000000000000000000000E+0/
      DATA C/3.828573121022477169108D-1, 1.107326627786831743809D+1,
     1       7.246689782858597021199D+1, 1.700632978311516129328D+2,
     2       1.698106763764238382705D+2, 7.633628843705946890896D+1,
     3       1.487967702840464066613D+1, 9.999989642347613068437D-1,
     4       1.737331760720576030932D-8/
      DATA D/8.258160008564488034698D-2, 4.344836335509282083360D+0,
     1       4.662179610356861756812D+1, 1.775728186717289799677D+2,
     2       2.953136335677908517423D+2, 2.342573504717625153053D+2,
     3       9.021658450529372642314D+1, 1.587964570758947927903D+1,
     4       1.000000000000000000000D+0/
C----------------------------------------------------------------------
C Coefficients for X < -4.0
C----------------------------------------------------------------------
CS    DATA E/1.3276881505637444622987E+2,3.5846198743996904308695E+4,
CS   1       1.7283375773777593926828E+5,2.6181454937205639647381E+5,
CS   2       1.7503273087497081314708E+5,5.9346841538837119172356E+4,
CS   3       1.0816852399095915622498E+4,1.0611777263550331766871E03,
CS   4       5.2199632588522572481039E+1,9.9999999999999999087819E-1/
CS    DATA F/3.9147856245556345627078E+4,2.5989762083608489777411E+5,
CS   1       5.5903756210022864003380E+5,5.4616842050691155735758E+5,
CS   2       2.7858134710520842139357E+5,7.9231787945279043698718E+4,
CS   3       1.2842808586627297365998E+4,1.1635769915320848035459E+3,
CS   4       5.4199632588522559414924E+1,1.0E0/
      DATA E/1.3276881505637444622987D+2,3.5846198743996904308695D+4,
     1       1.7283375773777593926828D+5,2.6181454937205639647381D+5,
     2       1.7503273087497081314708D+5,5.9346841538837119172356D+4,
     3       1.0816852399095915622498D+4,1.0611777263550331766871D03,
     4       5.2199632588522572481039D+1,9.9999999999999999087819D-1/
      DATA F/3.9147856245556345627078D+4,2.5989762083608489777411D+5,
     1       5.5903756210022864003380D+5,5.4616842050691155735758D+5,
     2       2.7858134710520842139357D+5,7.9231787945279043698718D+4,
     3       1.2842808586627297365998D+4,1.1635769915320848035459D+3,
     4       5.4199632588522559414924D+1,1.0D0/
C----------------------------------------------------------------------
C  Coefficients for rational approximation to ln(x/a), |1-x/a| < .1
C----------------------------------------------------------------------
CS    DATA PLG/-2.4562334077563243311E+01,2.3642701335621505212E+02,
CS   1         -5.4989956895857911039E+02,3.5687548468071500413E+02/
CS    DATA QLG/-3.5553900764052419184E+01,1.9400230218539473193E+02,
CS   1         -3.3442903192607538956E+02,1.7843774234035750207E+02/
      DATA PLG/-2.4562334077563243311D+01,2.3642701335621505212D+02,
     1         -5.4989956895857911039D+02,3.5687548468071500413D+02/
      DATA QLG/-3.5553900764052419184D+01,1.9400230218539473193D+02,
     1         -3.3442903192607538956D+02,1.7843774234035750207D+02/
C----------------------------------------------------------------------
C Coefficients for  0.0 < X < 6.0,
C  ratio of Chebyshev polynomials
C----------------------------------------------------------------------
CS    DATA P/-1.2963702602474830028590E01,-1.2831220659262000678155E03,
CS   1       -1.4287072500197005777376E04,-1.4299841572091610380064E06,
CS   2       -3.1398660864247265862050E05,-3.5377809694431133484800E08,
CS   3        3.1984354235237738511048E08,-2.5301823984599019348858E10,
CS   4        1.2177698136199594677580E10,-2.0829040666802497120940E11/
CS    DATA Q/ 7.6886718750000000000000E01,-5.5648470543369082846819E03,
CS   1        1.9418469440759880361415E05,-4.2648434812177161405483E06,
CS   2        6.4698830956576428587653E07,-7.0108568774215954065376E08,
CS   3        5.4229617984472955011862E09,-2.8986272696554495342658E10,
CS   4        9.8900934262481749439886E10,-8.9673749185755048616855E10/
      DATA P/-1.2963702602474830028590D01,-1.2831220659262000678155D03,
     1       -1.4287072500197005777376D04,-1.4299841572091610380064D06,
     2       -3.1398660864247265862050D05,-3.5377809694431133484800D08,
     3        3.1984354235237738511048D08,-2.5301823984599019348858D10,
     4        1.2177698136199594677580D10,-2.0829040666802497120940D11/
      DATA q/ 7.6886718750000000000000D01,-5.5648470543369082846819D03,
     1        1.9418469440759880361415D05,-4.2648434812177161405483D06,
     2        6.4698830956576428587653D07,-7.0108568774215954065376D08,
     3        5.4229617984472955011862D09,-2.8986272696554495342658D10,
     4        9.8900934262481749439886D10,-8.9673749185755048616855D10/
C----------------------------------------------------------------------
C J-fraction coefficients for 6.0 <= X < 12.0
C----------------------------------------------------------------------
CS    DATA R/-2.645677793077147237806E00,-2.378372882815725244124E00,
CS   1       -2.421106956980653511550E01, 1.052976392459015155422E01,
CS   2        1.945603779539281810439E01,-3.015761863840593359165E01,
CS   3        1.120011024227297451523E01,-3.988850730390541057912E00,
CS   4        9.565134591978630774217E00, 9.981193787537396413219E-1/
CS    DATA S/ 1.598517957704779356479E-4, 4.644185932583286942650E00,
CS   1        3.697412299772985940785E02,-8.791401054875438925029E00,
CS   2        7.608194509086645763123E02, 2.852397548119248700147E01,
CS   3        4.731097187816050252967E02,-2.369210235636181001661E02,
CS   4        1.249884822712447891440E00/
      DATA R/-2.645677793077147237806D00,-2.378372882815725244124D00,
     1       -2.421106956980653511550D01, 1.052976392459015155422D01,
     2        1.945603779539281810439D01,-3.015761863840593359165D01,
     3        1.120011024227297451523D01,-3.988850730390541057912D00,
     4        9.565134591978630774217D00, 9.981193787537396413219D-1/
      DATA S/ 1.598517957704779356479D-4, 4.644185932583286942650D00,
     1        3.697412299772985940785D02,-8.791401054875438925029D00,
     2        7.608194509086645763123D02, 2.852397548119248700147D01,
     3        4.731097187816050252967D02,-2.369210235636181001661D02,
     4        1.249884822712447891440D00/
C----------------------------------------------------------------------
C J-fraction coefficients for 12.0 <= X < 24.0
C----------------------------------------------------------------------
CS    DATA P1/-1.647721172463463140042E00,-1.860092121726437582253E01,
CS   1        -1.000641913989284829961E01,-2.105740799548040450394E01,
CS   2        -9.134835699998742552432E-1,-3.323612579343962284333E01,
CS   3         2.495487730402059440626E01, 2.652575818452799819855E01,
CS   4        -1.845086232391278674524E00, 9.999933106160568739091E-1/
CS    DATA Q1/ 9.792403599217290296840E01, 6.403800405352415551324E01,
CS   1         5.994932325667407355255E01, 2.538819315630708031713E02,
CS   2         4.429413178337928401161E01, 1.192832423968601006985E03,
CS   3         1.991004470817742470726E02,-1.093556195391091143924E01,
CS   4         1.001533852045342697818E00/
      DATA P1/-1.647721172463463140042D00,-1.860092121726437582253D01,
     1        -1.000641913989284829961D01,-2.105740799548040450394D01,
     2        -9.134835699998742552432D-1,-3.323612579343962284333D01,
     3         2.495487730402059440626D01, 2.652575818452799819855D01,
     4        -1.845086232391278674524D00, 9.999933106160568739091D-1/
      DATA Q1/ 9.792403599217290296840D01, 6.403800405352415551324D01,
     1         5.994932325667407355255D01, 2.538819315630708031713D02,
     2         4.429413178337928401161D01, 1.192832423968601006985D03,
     3         1.991004470817742470726D02,-1.093556195391091143924D01,
     4         1.001533852045342697818D00/
C----------------------------------------------------------------------
C J-fraction coefficients for  X .GE. 24.0
C----------------------------------------------------------------------
CS    DATA P2/ 1.75338801265465972390E02,-2.23127670777632409550E02,
CS   1        -1.81949664929868906455E01,-2.79798528624305389340E01,
CS   2        -7.63147701620253630855E00,-1.52856623636929636839E01,
CS   3        -7.06810977895029358836E00,-5.00006640413131002475E00,
CS   4        -3.00000000320981265753E00, 1.00000000000000485503E00/
CS    DATA Q2/ 3.97845977167414720840E04, 3.97277109100414518365E00,
CS   1         1.37790390235747998793E02, 1.17179220502086455287E02,
CS   2         7.04831847180424675988E01,-1.20187763547154743238E01,
CS   3        -7.99243595776339741065E00,-2.99999894040324959612E00,
CS   4         1.99999999999048104167E00/
      DATA P2/ 1.75338801265465972390D02,-2.23127670777632409550D02,
     1        -1.81949664929868906455D01,-2.79798528624305389340D01,
     2        -7.63147701620253630855D00,-1.52856623636929636839D01,
     3        -7.06810977895029358836D00,-5.00006640413131002475D00,
     4        -3.00000000320981265753D00, 1.00000000000000485503D00/
      DATA Q2/ 3.97845977167414720840D04, 3.97277109100414518365D00,
     1         1.37790390235747998793D02, 1.17179220502086455287D02,
     2         7.04831847180424675988D01,-1.20187763547154743238D01,
     3        -7.99243595776339741065D00,-2.99999894040324959612D00,
     4         1.99999999999048104167D00/
C----------------------------------------------------------------------
      X = ARG
      IF (X .EQ. ZERO) THEN
            EI = -XINF
            IF (INT .EQ. 2) EI = -EI
         ELSE IF ((X .LT. ZERO) .OR. (INT .EQ. 2)) THEN 
C----------------------------------------------------------------------
C Calculate EI for negative argument or for E1.
C----------------------------------------------------------------------
            Y = ABS(X)
            IF (Y .LE. ONE) THEN
                  SUMP = A(7) * Y + A(1)
                  SUMQ = Y + B(1)
                  DO 110 I = 2, 6
                     SUMP = SUMP * Y + A(I)
                     SUMQ = SUMQ * Y + B(I)
  110             CONTINUE
                  EI = LOG(Y) - SUMP / SUMQ
                  IF (INT .EQ. 3) EI = EI * EXP(Y)
               ELSE IF (Y .LE. FOUR) THEN
                  W = ONE / Y
                  SUMP = C(1)
                  SUMQ = D(1)
                  DO 130 I = 2, 9
                     SUMP = SUMP * W + C(I)
                     SUMQ = SUMQ * W + D(I)
  130             CONTINUE
                  EI = - SUMP / SUMQ
                  IF (INT .NE. 3) EI = EI * EXP(-Y)
               ELSE
                  IF ((Y .GT. XBIG) .AND. (INT .LT. 3)) THEN
                        EI = ZERO
                     ELSE
                        W = ONE / Y
                        SUMP = E(1) 
                        SUMQ = F(1)
                        DO 150 I = 2, 10
                           SUMP = SUMP * W + E(I)
                           SUMQ = SUMQ * W + F(I)
  150                   CONTINUE
                        EI = -W * (ONE - W * SUMP / SUMQ )
                        IF (INT .NE. 3) EI = EI * EXP(-Y)
                  END IF
            END IF
            IF (INT .EQ. 2) EI = -EI
         ELSE IF (X .LT. SIX) THEN
C----------------------------------------------------------------------
C  To improve conditioning, rational approximations are expressed
C    in terms of Chebyshev polynomials for 0 <= X < 6, and in
C    continued fraction form for larger X.
C----------------------------------------------------------------------
            T = X + X
            T = T / THREE - TWO
            PX(1) = ZERO
            QX(1) = ZERO
            PX(2) = P(1)
            QX(2) = q(1)
            DO 210 I = 2, 9
               PX(I+1) = T * PX(I) - PX(I-1) + P(I)
               QX(I+1) = T * QX(I) - QX(I-1) + q(I)
  210       CONTINUE
            SUMP = HALF * T * PX(10) - PX(9) + P(10)
            SUMQ = HALF * T * QX(10) - QX(9) + q(10)
            FRAC = SUMP / SUMQ
            XMX0 = (X - X01/X11) - X02
            IF (ABS(XMX0) .GE. P037) THEN
                  EI = LOG(X/X0) + XMX0 * FRAC
                  IF (INT .EQ. 3) EI = EXP(-X) * EI
               ELSE
C----------------------------------------------------------------------
C Special approximation to  ln(X/X0)  for X close to X0
C----------------------------------------------------------------------
                  Y = XMX0 / (X + X0)
                  YSQ = Y*Y
                  SUMP = PLG(1)
                  SUMQ = YSQ + QLG(1)
                  DO 220 I = 2, 4
                     SUMP = SUMP*YSQ + PLG(I)
                     SUMQ = SUMQ*YSQ + QLG(I)
  220             CONTINUE
                  EI = (SUMP / (SUMQ*(X+X0)) + FRAC) * XMX0
                  IF (INT .EQ. 3) EI = EXP(-X) * EI
            END IF
         ELSE IF (X .LT. TWELVE) THEN
            FRAC = ZERO
            DO 230 I = 1, 9
               FRAC = S(I) / (R(I) + X + FRAC)
  230       CONTINUE
            EI = (R(10) + FRAC) / X
            IF (INT .NE. 3) EI = EI * EXP(X)
         ELSE IF (X .LE. TWO4) THEN
            FRAC = ZERO
            DO 240 I = 1, 9
               FRAC = Q1(I) / (P1(I) + X + FRAC)
  240       CONTINUE
            EI = (P1(10) + FRAC) / X
            IF (INT .NE. 3) EI = EI * EXP(X)
         ELSE
            IF ((X .GE. XMAX) .AND. (INT .LT. 3)) THEN
                  EI = XINF
               ELSE
                  Y = ONE / X
                  FRAC = ZERO
                  DO 250 I = 1, 9
                     FRAC = Q2(I) / (P2(I) + X + FRAC)
  250             CONTINUE
                  FRAC = P2(10) + FRAC
                  EI = Y + Y * Y * FRAC
                  IF (INT .NE. 3) THEN
                        IF (X .LE. XMAX-TWO4) THEN
                              EI = EI * EXP(X)
                           ELSE
C----------------------------------------------------------------------
C Calculation reformulated to avoid premature overflow
C----------------------------------------------------------------------
                              EI = (EI * EXP(X-FOURTY)) * EXP40
                        END IF
                  END IF
            END IF
      END IF
      RESULT = EI
      RETURN
C---------- Last line of CALCEI ----------
      END
!      FUNCTION EI(X)
!C--------------------------------------------------------------------
!C
!C This function program computes approximate values for the
!C   exponential integral  Ei(x), where  x  is real.
!C
!C  Author: W. J. Cody
!C
!C  Latest modification: January 12, 1988
!C
!C--------------------------------------------------------------------
!      INTEGER INT
!CS    real*8  EI, X, RESULT
!CD    DOUBLE PRECISION  EI, X, RESULT
!      real*8  EI, X, RESULT
!C--------------------------------------------------------------------
!      INT = 1
!      CALL CALCEI(X,RESULT,INT)
!      EI = RESULT
!      RETURN
!C---------- Last line of EI ----------
!      END
      FUNCTION EXPEI(X)
C--------------------------------------------------------------------
C
C This function program computes approximate values for the
C   function  exp(-x) * Ei(x), where  Ei(x)  is the exponential
C   integral, and  x  is real.
C
C  Author: W. J. Cody
C
C  Latest modification: January 12, 1988
C
C--------------------------------------------------------------------
      INTEGER INT
CS    real*8  EXPEI, X, RESULT
      real*8  EXPEI, X, RESULT
CD    DOUBLE PRECISION  EXPEI, X, RESULT
C--------------------------------------------------------------------
      INT = 3
      CALL CALCEI(X,RESULT,INT)
      EXPEI = RESULT
      RETURN
C---------- Last line of EXPEI ----------
      END
!      FUNCTION E1(X)
!C--------------------------------------------------------------------
!C
!C This function program computes approximate values for the
!C   exponential integral E1(x), where  x  is real.
!C
!C  Author: W. J. Cody
!C
!C  Latest modification: January 12, 1988
!C
!C--------------------------------------------------------------------
!      INTEGER INT
!CS    real*8  EONE, X, RESULT
!CD    DOUBLE PRECISION  EONE, X, RESULT
!      real*8  E1, X, RESULT
!C--------------------------------------------------------------------
!      INT = 2
!      CALL CALCEI(X,RESULT,INT)
!      E1 = RESULT
!      RETURN
!C---------- Last line of E1 ----------
!      END
c:X_FT97Bsubrstart

c    Generated: Thu Jan 30 10:36:52 GMT 2003

      subroutine uks_x_ft97b
     & (ideriv,npt,rhoa1,rhob1,sigmaaa1,sigmabb1,sigmaab1,
     &  zk,vrhoa,vrhob,vsigmaaa,vsigmabb,vsigmaab,
     &  v2rhoa2,v2rhob2,v2rhoab,
     &  v2rhoasigmaaa,v2rhoasigmaab,v2rhoasigmabb,
     &  v2rhobsigmabb,v2rhobsigmaab,v2rhobsigmaaa,
     &  v2sigmaaa2,v2sigmaaaab,v2sigmaaabb,
     &  v2sigmaab2,v2sigmaabbb,v2sigmabb2)
c
c     M. Filatov, and W. Thiel
c     A new gradient-corrected exchange-correlation density functional
c     Mol. Phys. 91 (1997) 847-859
c
c
c     CITATION:
c
c     Functionals were obtained from the Density Functional Repository 
c     as developed and distributed by the Quantum Chemistry Group, 
c     CCLRC Daresbury Laboratory, Daresbury, Cheshire, WA4 4AD 
c     United Kingdom. Contact Huub van Dam (h.j.j.vandam@dl.ac.uk) or 
c     Paul Sherwood for further information.
c
c     COPYRIGHT:
c
c     Users may incorporate the source code into software packages and
c     redistribute the source code provided the source code is not
c     changed in anyway and is properly cited in any documentation or
c     publication related to its use.
c
c     ACKNOWLEDGEMENT:
c
c     The source code was generated using Maple 8 through a modified
c     version of the dfauto script published in:
c
c        R. Strange, F.R. Manby, P.J. Knowles
c        Automatic code generation in density functional theory
c        Comp. Phys. Comm. 136 (2001) 310-318.
c
      implicit real*8 (a-h,o-z)
      integer ideriv,npt
      real*8 rhoa1(npt),rhob1(npt)
      real*8 sigmaaa1(npt),sigmabb1(npt),sigmaab1(npt)
      real*8 zk(npt),vrhoa(npt),vrhob(npt)
      real*8 vsigmaaa(npt),vsigmabb(npt),vsigmaab(npt)
      real*8 v2rhoa2(npt),v2rhob2(npt),v2rhoab(npt)
      real*8 v2rhoasigmaaa(npt),v2rhoasigmaab(npt)
      real*8 v2rhoasigmabb(npt),v2rhobsigmabb(npt)
      real*8 v2rhobsigmaab(npt),v2rhobsigmaaa(npt)
      real*8 v2sigmaaa2(npt),v2sigmaaaab(npt),v2sigmaaabb(npt)
      real*8 v2sigmaab2(npt),v2sigmaabbb(npt),v2sigmabb2(npt)
      parameter(tol=1.0d-20)
      
      if (ideriv.eq.0) then
      
      do i=1,npt
      rhoa = dmax1(0.D0,rhoa1(i))
      rhob = dmax1(0.D0,rhob1(i))
      rho = rhoa+rhob
      if(rho.gt.tol) then
      if(rhoa.lt.tol) then
      rho = rhob
      sigmabb = dmax1(0.D0,sigmabb1(i))
      sigma = sigmabb
      t2 = rhob**(1.D0/3.D0)
      t8 = 0.2913644D-2+0.9474169D-3*sigmabb/(0.6255746320201D7+sigmabb)
      t10 = rhob**2
      t11 = t2**2
      t13 = 1/t11/t10
      t14 = t8**2
      t17 = dlog(sigmabb*t13+dsqrt(1+sigmabb**2*t13**2))
      t18 = t17**2
      t23 = dsqrt(1.D0+9.D0*t14*sigmabb*t13*t18)
      zk(i) = -0.9305257363491D0*t2*rhob*(1.D0+0.1074661302677646D1*t8
     &*sigmabb*t13/t23)
      elseif(rhob.lt.tol) then
      rho = rhoa
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigma = sigmaaa
      t2 = rhoa**(1.D0/3.D0)
      t8 = 0.2913644D-2+0.9474169D-3*sigmaaa/(0.6255746320201D7+sigmaaa)
      t10 = rhoa**2
      t11 = t2**2
      t13 = 1/t11/t10
      t14 = t8**2
      t17 = dlog(sigmaaa*t13+dsqrt(1+sigmaaa**2*t13**2))
      t18 = t17**2
      t23 = dsqrt(1.D0+9.D0*t14*sigmaaa*t13*t18)
      zk(i) = -0.9305257363491D0*t2*rhoa*(1.D0+0.1074661302677646D1*t8
     &*sigmaaa*t13/t23)
      else ! (.not.(rhoa.lt.tol).and..not.(rhob.lt.tol))
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigmaab = sigmaab1(i)
      sigmabb = dmax1(0.D0,sigmabb1(i))
      sigma = sigmaaa+sigmabb+2.D0*sigmaab
      t4 = rhoa**(1.D0/3.D0)
      t10 = 0.2913644D-2+0.9474169D-3*sigmaaa/(0.6255746320201D7
     &+sigmaaa)
      t12 = rhoa**2
      t13 = t4**2
      t15 = 1/t13/t12
      t16 = t10**2
      t19 = dlog(sigmaaa*t15+dsqrt(1+sigmaaa**2*t15**2))
      t20 = t19**2
      t25 = dsqrt(1.D0+9.D0*t16*sigmaaa*t15*t20)
      t33 = rhob**(1.D0/3.D0)
      t39 = 0.2913644D-2+0.9474169D-3*sigmabb/(0.6255746320201D7
     &+sigmabb)
      t41 = rhob**2
      t42 = t33**2
      t44 = 1/t42/t41
      t45 = t39**2
      t48 = dlog(sigmabb*t44+dsqrt(1+sigmabb**2*t44**2))
      t49 = t48**2
      t54 = dsqrt(1.D0+9.D0*t45*sigmabb*t44*t49)
      zk(i) = -0.9305257363491D0*t4*rhoa*(1.D0+0.1074661302677646D1
     &*t10*sigmaaa*t15/t25)-0.9305257363491D0*t33*rhob*(1.D0
     &+0.1074661302677646D1*t39*sigmabb*t44/t54)
      endif ! rhoa,rhob
      else ! rho
      zk(i) = 0.0d0
      endif ! rho
      enddo
      
      else if(ideriv.eq.1) then
      
      do i=1,npt
      rhoa = dmax1(0.D0,rhoa1(i))
      rhob = dmax1(0.D0,rhob1(i))
      rho = rhoa+rhob
      if(rho.gt.tol) then
      if(rhoa.lt.tol) then
      rho = rhob
      sigmabb = dmax1(0.D0,sigmabb1(i))
      sigma = sigmabb
      t2 = rhob**(1.D0/3.D0)
      t3 = t2*rhob
      t4 = 0.6255746320201D7+sigmabb
      t5 = 1/t4
      t8 = 0.2913644D-2+0.9474169D-3*sigmabb*t5
      t9 = t8*sigmabb
      t10 = rhob**2
      t11 = t2**2
      t13 = 1/t11/t10
      t14 = t8**2
      t15 = t14*sigmabb
      t17 = dlog(sigmabb*t13+dsqrt(1+sigmabb**2*t13**2))
      t18 = t17**2
      t19 = t13*t18
      t22 = 1.D0+9.D0*t15*t19
      t23 = dsqrt(t22)
      t24 = 1/t23
      t25 = t13*t24
      t28 = 1.D0+0.1074661302677646D1*t9*t25
      zk(i) = -0.9305257363491D0*t3*t28
      vrhoa(i) = 0.D0
      t35 = 1/t11/t10/rhob
      t41 = t13/t23/t22
      t45 = sigmabb**2
      t47 = t10**2
      t54 = 1/t2/t47/rhob
      t57 = dsqrt(1.D0+t45*t54)
      t58 = 1/t57
      vrhob(i) = -0.12407009817988D1*t2*t28-0.9305257363491D0*t3*(
     &-0.2865763473807057D1*t9*t35*t24-0.5373306513388232D0*t9*t41*(
     &-24.D0*t15*t35*t18-48.D0*t14*t45/t2/t47/t10*t17*t58))
      vsigmaaa(i) = 0.D0
      vsigmaab(i) = 0.D0
      t70 = t4**2
      t74 = 0.9474169D-3*t5-0.9474169D-3*sigmabb/t70
      vsigmabb(i) = -0.9305257363491D0*t3*(0.1074661302677646D1*t74
     &*sigmabb*t25+0.1074661302677646D1*t8*t13*t24
     &-0.5373306513388232D0*t9*t41*(18.D0*t9*t19*t74+9.D0*t14*t13*t18
     &+18.D0*t15*t54*t17*t58))
      elseif(rhob.lt.tol) then
      rho = rhoa
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigma = sigmaaa
      t2 = rhoa**(1.D0/3.D0)
      t3 = t2*rhoa
      t4 = 0.6255746320201D7+sigmaaa
      t5 = 1/t4
      t8 = 0.2913644D-2+0.9474169D-3*sigmaaa*t5
      t9 = t8*sigmaaa
      t10 = rhoa**2
      t11 = t2**2
      t13 = 1/t11/t10
      t14 = t8**2
      t15 = t14*sigmaaa
      t17 = dlog(sigmaaa*t13+dsqrt(1+sigmaaa**2*t13**2))
      t18 = t17**2
      t19 = t13*t18
      t22 = 1.D0+9.D0*t15*t19
      t23 = dsqrt(t22)
      t24 = 1/t23
      t25 = t13*t24
      t28 = 1.D0+0.1074661302677646D1*t9*t25
      zk(i) = -0.9305257363491D0*t3*t28
      t35 = 1/t11/t10/rhoa
      t41 = t13/t23/t22
      t45 = sigmaaa**2
      t47 = t10**2
      t54 = 1/t2/t47/rhoa
      t57 = dsqrt(1.D0+t45*t54)
      t58 = 1/t57
      vrhoa(i) = -0.12407009817988D1*t2*t28-0.9305257363491D0*t3*(
     &-0.2865763473807057D1*t9*t35*t24-0.5373306513388232D0*t9*t41*(
     &-24.D0*t15*t35*t18-48.D0*t14*t45/t2/t47/t10*t17*t58))
      vrhob(i) = 0.D0
      t70 = t4**2
      t74 = 0.9474169D-3*t5-0.9474169D-3*sigmaaa/t70
      vsigmaaa(i) = -0.9305257363491D0*t3*(0.1074661302677646D1*t74
     &*sigmaaa*t25+0.1074661302677646D1*t8*t13*t24
     &-0.5373306513388232D0*t9*t41*(18.D0*t9*t19*t74+9.D0*t14*t13*t18
     &+18.D0*t15*t54*t17*t58))
      vsigmaab(i) = 0.D0
      vsigmabb(i) = 0.D0
      else ! (.not.(rhoa.lt.tol).and..not.(rhob.lt.tol))
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigmaab = sigmaab1(i)
      sigmabb = dmax1(0.D0,sigmabb1(i))
      sigma = sigmaaa+sigmabb+2.D0*sigmaab
      t4 = rhoa**(1.D0/3.D0)
      t5 = t4*rhoa
      t6 = 0.6255746320201D7+sigmaaa
      t7 = 1/t6
      t10 = 0.2913644D-2+0.9474169D-3*sigmaaa*t7
      t11 = t10*sigmaaa
      t12 = rhoa**2
      t13 = t4**2
      t15 = 1/t13/t12
      t16 = t10**2
      t17 = t16*sigmaaa
      t19 = dlog(sigmaaa*t15+dsqrt(1+sigmaaa**2*t15**2))
      t20 = t19**2
      t21 = t15*t20
      t24 = 1.D0+9.D0*t17*t21
      t25 = dsqrt(t24)
      t26 = 1/t25
      t27 = t15*t26
      t30 = 1.D0+0.1074661302677646D1*t11*t27
      t33 = rhob**(1.D0/3.D0)
      t34 = t33*rhob
      t35 = 0.6255746320201D7+sigmabb
      t36 = 1/t35
      t39 = 0.2913644D-2+0.9474169D-3*sigmabb*t36
      t40 = t39*sigmabb
      t41 = rhob**2
      t42 = t33**2
      t44 = 1/t42/t41
      t45 = t39**2
      t46 = t45*sigmabb
      t48 = dlog(sigmabb*t44+dsqrt(1+sigmabb**2*t44**2))
      t49 = t48**2
      t50 = t44*t49
      t53 = 1.D0+9.D0*t46*t50
      t54 = dsqrt(t53)
      t55 = 1/t54
      t56 = t44*t55
      t59 = 1.D0+0.1074661302677646D1*t40*t56
      zk(i) = -0.9305257363491D0*t5*t30-0.9305257363491D0*t34*t59
      t66 = 1/t13/t12/rhoa
      t72 = t15/t25/t24
      t76 = sigmaaa**2
      t78 = t12**2
      t85 = 1/t4/t78/rhoa
      t88 = dsqrt(1.D0+t76*t85)
      t89 = 1/t88
      vrhoa(i) = -0.12407009817988D1*t4*t30-0.9305257363491D0*t5*(
     &-0.2865763473807057D1*t11*t66*t26-0.5373306513388232D0*t11*t72*(
     &-24.D0*t17*t66*t20-48.D0*t16*t76/t4/t78/t12*t19*t89))
      t104 = 1/t42/t41/rhob
      t110 = t44/t54/t53
      t114 = sigmabb**2
      t116 = t41**2
      t123 = 1/t33/t116/rhob
      t126 = dsqrt(1.D0+t114*t123)
      t127 = 1/t126
      vrhob(i) = -0.12407009817988D1*t33*t59-0.9305257363491D0*t34*(
     &-0.2865763473807057D1*t40*t104*t55-0.5373306513388232D0*t40*t110
     &*(-24.D0*t46*t104*t49-48.D0*t45*t114/t33/t116/t41*t48*t127))
      t139 = t6**2
      t143 = 0.9474169D-3*t7-0.9474169D-3*sigmaaa/t139
      vsigmaaa(i) = -0.9305257363491D0*t5*(0.1074661302677646D1*t143
     &*sigmaaa*t27+0.1074661302677646D1*t10*t15*t26
     &-0.5373306513388232D0*t11*t72*(18.D0*t11*t21*t143+9.D0*t16*t15
     &*t20+18.D0*t17*t85*t19*t89))
      vsigmaab(i) = 0.D0
      t168 = t35**2
      t172 = 0.9474169D-3*t36-0.9474169D-3*sigmabb/t168
      vsigmabb(i) = -0.9305257363491D0*t34*(0.1074661302677646D1*t172
     &*sigmabb*t56+0.1074661302677646D1*t39*t44*t55
     &-0.5373306513388232D0*t40*t110*(18.D0*t40*t50*t172+9.D0*t45*t44
     &*t49+18.D0*t46*t123*t48*t127))
      endif ! rhoa,rhob
      else ! rho
      zk(i) = 0.0d0
      vrhoa(i) = 0.0d0
      vrhob(i) = 0.0d0
      vsigmaaa(i) = 0.0d0
      vsigmaab(i) = 0.0d0
      vsigmabb(i) = 0.0d0
      endif ! rho
      enddo
      
      else if(ideriv.eq.2) then
      
      do i=1,npt
      rhoa = dmax1(0.D0,rhoa1(i))
      rhob = dmax1(0.D0,rhob1(i))
      rho = rhoa+rhob
      if(rho.gt.tol) then
      if(rhoa.lt.tol) then
      rho = rhob
      sigmabb = dmax1(0.D0,sigmabb1(i))
      sigma = sigmabb
      t2 = rhob**(1.D0/3.D0)
      t3 = t2*rhob
      t4 = 0.6255746320201D7+sigmabb
      t5 = 1/t4
      t8 = 0.2913644D-2+0.9474169D-3*sigmabb*t5
      t9 = t8*sigmabb
      t10 = rhob**2
      t11 = t2**2
      t13 = 1/t11/t10
      t14 = t8**2
      t15 = t14*sigmabb
      t17 = dlog(sigmabb*t13+dsqrt(1+sigmabb**2*t13**2))
      t18 = t17**2
      t19 = t13*t18
      t22 = 1.D0+9.D0*t15*t19
      t23 = dsqrt(t22)
      t24 = 1/t23
      t25 = t13*t24
      t28 = 1.D0+0.1074661302677646D1*t9*t25
      zk(i) = -0.9305257363491D0*t3*t28
      vrhoa(i) = 0.D0
      t33 = t10*rhob
      t35 = 1/t11/t33
      t40 = 1/t23/t22
      t41 = t13*t40
      t45 = sigmabb**2
      t46 = t14*t45
      t47 = t10**2
      t54 = 1/t2/t47/rhob
      t56 = 1.D0+t45*t54
      t57 = dsqrt(t56)
      t58 = 1/t57
      t62 = -24.D0*t15*t35*t18-48.D0*t46/t2/t47/t10*t17*t58
      t66 = -0.2865763473807057D1*t9*t35*t24-0.5373306513388232D0*t9
     &*t41*t62
      vrhob(i) = -0.12407009817988D1*t2*t28-0.9305257363491D0*t3*t66
      vsigmaaa(i) = 0.D0
      vsigmaab(i) = 0.D0
      t70 = t4**2
      t71 = 1/t70
      t74 = 0.9474169D-3*t5-0.9474169D-3*sigmabb*t71
      t75 = t74*sigmabb
      t78 = t8*t13
      t91 = 18.D0*t9*t19*t74+9.D0*t14*t13*t18+18.D0*t15*t54*t17*t58
      t92 = t41*t91
      vsigmabb(i) = -0.9305257363491D0*t3*(0.1074661302677646D1*t75
     &*t25+0.1074661302677646D1*t78*t24-0.5373306513388232D0*t9*t92)
      v2rhoa2(i) = 0.D0
      v2rhoab(i) = 0.D0
      t104 = 1/t11/t47
      t112 = t22**2
      t115 = t13/t23/t112
      t116 = t62**2
      t132 = t47**2
      t133 = t132*t10
      t135 = 1/t56
      t139 = t45**2
      t146 = 1/t57/t56
      v2rhob2(i) = -0.4135669939329333D0/t11*t28-0.24814019635976D1*t2
     &*t66-0.9305257363491D0*t3*(0.1050779940395921D2*t9*t104*t24
     &+0.2865763473807057D1*t9*t35*t40*t62+0.8059959770082348D0*t9
     &*t115*t116-0.5373306513388232D0*t9*t41*(88.D0*t15*t104*t18
     &+432.D0*t46/t2/t47/t33*t17*t58+128.D0*t14*t45*sigmabb/t133*t135
     &-128.D0*t14*t139/t11/t132/t47*t17*t146))
      v2sigmaaa2(i) = 0.D0
      v2sigmaaaab(i) = 0.D0
      v2sigmaaabb(i) = 0.D0
      v2sigmaab2(i) = 0.D0
      v2sigmaabbb(i) = 0.D0
      t162 = -0.18948338D-2*t71+0.18948338D-2*sigmabb/t70/t4
      t174 = t91**2
      t178 = t74**2
      v2sigmabb2(i) = -0.9305257363491D0*t3*(0.1074661302677646D1*t162
     &*sigmabb*t25+0.2149322605355293D1*t74*t13*t24
     &-0.1074661302677646D1*t75*t92-0.1074661302677646D1*t78*t40*t91
     &+0.8059959770082348D0*t9*t115*t174-0.5373306513388232D0*t9*t41*
     &(18.D0*t178*sigmabb*t19+36.D0*t78*t18*t74+72.D0*t9*t54*t17*t74
     &*t58+18.D0*t9*t19*t162+36.D0*t14*t54*t17*t58+18.D0*t15/t132*t135
     &-18.D0*t46/t11/t133*t17*t146))
      elseif(rhob.lt.tol) then
      rho = rhoa
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigma = sigmaaa
      t2 = rhoa**(1.D0/3.D0)
      t3 = t2*rhoa
      t4 = 0.6255746320201D7+sigmaaa
      t5 = 1/t4
      t8 = 0.2913644D-2+0.9474169D-3*sigmaaa*t5
      t9 = t8*sigmaaa
      t10 = rhoa**2
      t11 = t2**2
      t13 = 1/t11/t10
      t14 = t8**2
      t15 = t14*sigmaaa
      t17 = dlog(sigmaaa*t13+dsqrt(1+sigmaaa**2*t13**2))
      t18 = t17**2
      t19 = t13*t18
      t22 = 1.D0+9.D0*t15*t19
      t23 = dsqrt(t22)
      t24 = 1/t23
      t25 = t13*t24
      t28 = 1.D0+0.1074661302677646D1*t9*t25
      zk(i) = -0.9305257363491D0*t3*t28
      t33 = t10*rhoa
      t35 = 1/t11/t33
      t40 = 1/t23/t22
      t41 = t13*t40
      t45 = sigmaaa**2
      t46 = t14*t45
      t47 = t10**2
      t54 = 1/t2/t47/rhoa
      t56 = 1.D0+t45*t54
      t57 = dsqrt(t56)
      t58 = 1/t57
      t62 = -24.D0*t15*t35*t18-48.D0*t46/t2/t47/t10*t17*t58
      t66 = -0.2865763473807057D1*t9*t35*t24-0.5373306513388232D0*t9
     &*t41*t62
      vrhoa(i) = -0.12407009817988D1*t2*t28-0.9305257363491D0*t3*t66
      vrhob(i) = 0.D0
      t70 = t4**2
      t71 = 1/t70
      t74 = 0.9474169D-3*t5-0.9474169D-3*sigmaaa*t71
      t75 = t74*sigmaaa
      t78 = t8*t13
      t91 = 18.D0*t9*t19*t74+9.D0*t14*t13*t18+18.D0*t15*t54*t17*t58
      t92 = t41*t91
      vsigmaaa(i) = -0.9305257363491D0*t3*(0.1074661302677646D1*t75
     &*t25+0.1074661302677646D1*t78*t24-0.5373306513388232D0*t9*t92)
      vsigmaab(i) = 0.D0
      vsigmabb(i) = 0.D0
      t104 = 1/t11/t47
      t112 = t22**2
      t115 = t13/t23/t112
      t116 = t62**2
      t132 = t47**2
      t133 = t132*t10
      t135 = 1/t56
      t139 = t45**2
      t146 = 1/t57/t56
      v2rhoa2(i) = -0.4135669939329333D0/t11*t28-0.24814019635976D1*t2
     &*t66-0.9305257363491D0*t3*(0.1050779940395921D2*t9*t104*t24
     &+0.2865763473807057D1*t9*t35*t40*t62+0.8059959770082348D0*t9
     &*t115*t116-0.5373306513388232D0*t9*t41*(88.D0*t15*t104*t18
     &+432.D0*t46/t2/t47/t33*t17*t58+128.D0*t14*t45*sigmaaa/t133*t135
     &-128.D0*t14*t139/t11/t132/t47*t17*t146))
      v2rhob2(i) = 0.D0
      v2rhoab(i) = 0.D0
      t162 = -0.18948338D-2*t71+0.18948338D-2*sigmaaa/t70/t4
      t174 = t91**2
      t178 = t74**2
      v2sigmaaa2(i) = -0.9305257363491D0*t3*(0.1074661302677646D1*t162
     &*sigmaaa*t25+0.2149322605355293D1*t74*t13*t24
     &-0.1074661302677646D1*t75*t92-0.1074661302677646D1*t78*t40*t91
     &+0.8059959770082348D0*t9*t115*t174-0.5373306513388232D0*t9*t41*
     &(18.D0*t178*sigmaaa*t19+36.D0*t78*t18*t74+72.D0*t9*t54*t17*t74
     &*t58+18.D0*t9*t19*t162+36.D0*t14*t54*t17*t58+18.D0*t15/t132*t135
     &-18.D0*t46/t11/t133*t17*t146))
      v2sigmaaaab(i) = 0.D0
      v2sigmaaabb(i) = 0.D0
      v2sigmaab2(i) = 0.D0
      v2sigmaabbb(i) = 0.D0
      v2sigmabb2(i) = 0.D0
      else ! (.not.(rhoa.lt.tol).and..not.(rhob.lt.tol))
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigmaab = sigmaab1(i)
      sigmabb = dmax1(0.D0,sigmabb1(i))
      sigma = sigmaaa+sigmabb+2.D0*sigmaab
      t4 = rhoa**(1.D0/3.D0)
      t5 = t4*rhoa
      t6 = 0.6255746320201D7+sigmaaa
      t7 = 1/t6
      t10 = 0.2913644D-2+0.9474169D-3*sigmaaa*t7
      t11 = t10*sigmaaa
      t12 = rhoa**2
      t13 = t4**2
      t15 = 1/t13/t12
      t16 = t10**2
      t17 = t16*sigmaaa
      t19 = dlog(sigmaaa*t15+dsqrt(1+sigmaaa**2*t15**2))
      t20 = t19**2
      t21 = t15*t20
      t24 = 1.D0+9.D0*t17*t21
      t25 = dsqrt(t24)
      t26 = 1/t25
      t27 = t15*t26
      t30 = 1.D0+0.1074661302677646D1*t11*t27
      t33 = rhob**(1.D0/3.D0)
      t34 = t33*rhob
      t35 = 0.6255746320201D7+sigmabb
      t36 = 1/t35
      t39 = 0.2913644D-2+0.9474169D-3*sigmabb*t36
      t40 = t39*sigmabb
      t41 = rhob**2
      t42 = t33**2
      t44 = 1/t42/t41
      t45 = t39**2
      t46 = t45*sigmabb
      t48 = dlog(sigmabb*t44+dsqrt(1+sigmabb**2*t44**2))
      t49 = t48**2
      t50 = t44*t49
      t53 = 1.D0+9.D0*t46*t50
      t54 = dsqrt(t53)
      t55 = 1/t54
      t56 = t44*t55
      t59 = 1.D0+0.1074661302677646D1*t40*t56
      zk(i) = -0.9305257363491D0*t5*t30-0.9305257363491D0*t34*t59
      t64 = t12*rhoa
      t66 = 1/t13/t64
      t67 = t66*t26
      t71 = 1/t25/t24
      t72 = t15*t71
      t73 = t66*t20
      t76 = sigmaaa**2
      t77 = t16*t76
      t78 = t12**2
      t81 = 1/t4/t78/t12
      t85 = 1/t4/t78/rhoa
      t87 = 1.D0+t76*t85
      t88 = dsqrt(t87)
      t89 = 1/t88
      t90 = t81*t19*t89
      t93 = -24.D0*t17*t73-48.D0*t77*t90
      t94 = t72*t93
      t97 = -0.2865763473807057D1*t11*t67-0.5373306513388232D0*t11*t94
      vrhoa(i) = -0.12407009817988D1*t4*t30-0.9305257363491D0*t5*t97
      t102 = t41*rhob
      t104 = 1/t42/t102
      t105 = t104*t55
      t109 = 1/t54/t53
      t110 = t44*t109
      t111 = t104*t49
      t114 = sigmabb**2
      t115 = t45*t114
      t116 = t41**2
      t119 = 1/t33/t116/t41
      t123 = 1/t33/t116/rhob
      t125 = 1.D0+t114*t123
      t126 = dsqrt(t125)
      t127 = 1/t126
      t128 = t119*t48*t127
      t131 = -24.D0*t46*t111-48.D0*t115*t128
      t132 = t110*t131
      t135 = -0.2865763473807057D1*t40*t105-0.5373306513388232D0*t40
     &*t132
      vrhob(i) = -0.12407009817988D1*t33*t59-0.9305257363491D0*t34*t135
      t139 = t6**2
      t140 = 1/t139
      t143 = 0.9474169D-3*t7-0.9474169D-3*sigmaaa*t140
      t144 = t143*sigmaaa
      t147 = t10*t15
      t160 = 18.D0*t11*t21*t143+9.D0*t16*t15*t20+18.D0*t17*t85*t19*t89
      t161 = t72*t160
      t164 = 0.1074661302677646D1*t144*t27+0.1074661302677646D1*t147
     &*t26-0.5373306513388232D0*t11*t161
      vsigmaaa(i) = -0.9305257363491D0*t5*t164
      vsigmaab(i) = 0.D0
      t168 = t35**2
      t169 = 1/t168
      t172 = 0.9474169D-3*t36-0.9474169D-3*sigmabb*t169
      t173 = t172*sigmabb
      t176 = t39*t44
      t189 = 18.D0*t40*t50*t172+9.D0*t45*t44*t49+18.D0*t46*t123*t48*t127
      t190 = t110*t189
      t193 = 0.1074661302677646D1*t173*t56+0.1074661302677646D1*t176
     &*t55-0.5373306513388232D0*t40*t190
      vsigmabb(i) = -0.9305257363491D0*t34*t193
      t202 = 1/t13/t78
      t206 = t66*t71
      t210 = t24**2
      t212 = 1/t25/t210
      t213 = t15*t212
      t214 = t93**2
      t229 = t16*t76*sigmaaa
      t230 = t78**2
      t231 = t230*t12
      t233 = 1/t87
      t237 = t76**2
      t244 = 1/t88/t87
      v2rhoa2(i) = -0.4135669939329333D0/t13*t30-0.24814019635976D1*t4
     &*t97-0.9305257363491D0*t5*(0.1050779940395921D2*t11*t202*t26
     &+0.2865763473807057D1*t11*t206*t93+0.8059959770082348D0*t11*t213
     &*t214-0.5373306513388232D0*t11*t72*(88.D0*t17*t202*t20+432.D0
     &*t77/t4/t78/t64*t19*t89+128.D0*t229/t231*t233-128.D0*t16*t237
     &/t13/t230/t78*t19*t244))
      t261 = 1/t42/t116
      t265 = t104*t109
      t269 = t53**2
      t271 = 1/t54/t269
      t272 = t44*t271
      t273 = t131**2
      t288 = t45*t114*sigmabb
      t289 = t116**2
      t290 = t289*t41
      t292 = 1/t125
      t296 = t114**2
      t303 = 1/t126/t125
      v2rhob2(i) = -0.4135669939329333D0/t42*t59-0.24814019635976D1
     &*t33*t135-0.9305257363491D0*t34*(0.1050779940395921D2*t40*t261
     &*t55+0.2865763473807057D1*t40*t265*t131+0.8059959770082348D0*t40
     &*t272*t273-0.5373306513388232D0*t40*t110*(88.D0*t46*t261*t49
     &+432.D0*t115/t33/t116/t102*t48*t127+128.D0*t288/t290*t292-128.D0
     &*t45*t296/t42/t289/t116*t48*t303))
      v2rhoab(i) = 0.D0
      t344 = t19*t89
      t345 = t344*t143
      v2rhoasigmaaa(i) = -0.12407009817988D1*t4*t164-0.9305257363491D0
     &*t5*(-0.2865763473807057D1*t144*t67-0.2865763473807057D1*t10*t66
     &*t26+0.1432881736903529D1*t11*t206*t160-0.5373306513388232D0
     &*t144*t94-0.5373306513388232D0*t147*t71*t93+0.8059959770082348D0
     &*t11*t15*t212*t93*t160-0.5373306513388232D0*t11*t72*(-48.D0*t11
     &*t73*t143-24.D0*t16*t66*t20-144.D0*t17*t90-96.D0*t10*t76*t81
     &*t345-48.D0*t77/t230/rhoa*t233+48.D0*t229/t13/t230/t64*t19*t244))
      v2rhoasigmaab(i) = 0.D0
      v2rhoasigmabb(i) = 0.D0
      v2rhobsigmaaa(i) = 0.D0
      v2rhobsigmaab(i) = 0.D0
      t397 = t48*t127
      t398 = t397*t172
      v2rhobsigmabb(i) = -0.12407009817988D1*t33*t193
     &-0.9305257363491D0*t34*(-0.2865763473807057D1*t173*t105
     &-0.2865763473807057D1*t39*t104*t55+0.1432881736903529D1*t40*t265
     &*t189-0.5373306513388232D0*t173*t132-0.5373306513388232D0*t176
     &*t109*t131+0.8059959770082348D0*t40*t44*t271*t131*t189
     &-0.5373306513388232D0*t40*t110*(-48.D0*t40*t111*t172-24.D0*t45
     &*t104*t49-144.D0*t46*t128-96.D0*t39*t114*t119*t398-48.D0*t115
     &/t289/rhob*t292+48.D0*t288/t42/t289/t102*t48*t303))
      t425 = -0.18948338D-2*t140+0.18948338D-2*sigmaaa/t139/t6
      t437 = t160**2
      t441 = t143**2
      v2sigmaaa2(i) = -0.9305257363491D0*t5*(0.1074661302677646D1*t425
     &*sigmaaa*t27+0.2149322605355293D1*t143*t15*t26
     &-0.1074661302677646D1*t144*t161-0.1074661302677646D1*t147*t71
     &*t160+0.8059959770082348D0*t11*t213*t437-0.5373306513388232D0
     &*t11*t72*(18.D0*t441*sigmaaa*t21+36.D0*t147*t20*t143+72.D0*t11
     &*t85*t345+18.D0*t11*t21*t425+36.D0*t16*t85*t344+18.D0*t17/t230
     &*t233-18.D0*t77/t13/t231*t19*t244))
      v2sigmaaaab(i) = 0.D0
      v2sigmaaabb(i) = 0.D0
      v2sigmaab2(i) = 0.D0
      v2sigmaabbb(i) = 0.D0
      t479 = -0.18948338D-2*t169+0.18948338D-2*sigmabb/t168/t35
      t491 = t189**2
      t495 = t172**2
      v2sigmabb2(i) = -0.9305257363491D0*t34*(0.1074661302677646D1
     &*t479*sigmabb*t56+0.2149322605355293D1*t172*t44*t55
     &-0.1074661302677646D1*t173*t190-0.1074661302677646D1*t176*t109
     &*t189+0.8059959770082348D0*t40*t272*t491-0.5373306513388232D0
     &*t40*t110*(18.D0*t495*sigmabb*t50+36.D0*t176*t49*t172+72.D0*t40
     &*t123*t398+18.D0*t40*t50*t479+36.D0*t45*t123*t397+18.D0*t46/t289
     &*t292-18.D0*t115/t42/t290*t48*t303))
      endif ! rhoa,rhob
      else ! rho
      zk(i) = 0.0d0
      vrhoa(i) = 0.0d0
      vrhob(i) = 0.0d0
      v2rhoa2(i) = 0.0d0
      v2rhob2(i) = 0.0d0
      v2rhoab(i) = 0.0d0
      vsigmaaa(i) = 0.0d0
      vsigmaab(i) = 0.0d0
      vsigmabb(i) = 0.0d0
      v2rhoasigmaaa(i) = 0.0d0
      v2rhoasigmaab(i) = 0.0d0
      v2rhoasigmabb(i) = 0.0d0
      v2rhobsigmaaa(i) = 0.0d0
      v2rhobsigmaab(i) = 0.0d0
      v2rhobsigmabb(i) = 0.0d0
      v2sigmaaa2(i) = 0.0d0
      v2sigmaab2(i) = 0.0d0
      v2sigmabb2(i) = 0.0d0
      v2sigmaaaab(i) = 0.0d0
      v2sigmaaabb(i) = 0.0d0
      v2sigmaabbb(i) = 0.0d0
      endif ! rho
      enddo
      
      endif ! ideriv
      return
      end
      
      
      subroutine rks_x_ft97b
     & (ideriv,npt,rhoa1,sigmaaa1,
     &  zk,vrhoa,vsigmaaa,
     &  v2rhoa2,v2rhoasigmaaa,v2sigmaaa2)
c
c     M. Filatov, and W. Thiel
c     A new gradient-corrected exchange-correlation density functional
c     Mol. Phys. 91 (1997) 847-859
c
c
c     CITATION:
c
c     Functionals were obtained from the Density Functional Repository 
c     as developed and distributed by the Quantum Chemistry Group, 
c     CCLRC Daresbury Laboratory, Daresbury, Cheshire, WA4 4AD 
c     United Kingdom. Contact Huub van Dam (h.j.j.vandam@dl.ac.uk) or 
c     Paul Sherwood for further information.
c
c     COPYRIGHT:
c
c     Users may incorporate the source code into software packages and
c     redistribute the source code provided the source code is not
c     changed in anyway and is properly cited in any documentation or
c     publication related to its use.
c
c     ACKNOWLEDGEMENT:
c
c     The source code was generated using Maple 8 through a modified
c     version of the dfauto script published in:
c
c        R. Strange, F.R. Manby, P.J. Knowles
c        Automatic code generation in density functional theory
c        Comp. Phys. Comm. 136 (2001) 310-318.
c
      implicit real*8 (a-h,o-z)
      integer ideriv,npt
      real*8 rhoa1(npt)
      real*8 sigmaaa1(npt)
      real*8 zk(npt),vrhoa(npt),vsigmaaa(npt)
      real*8 v2rhoa2(npt),v2rhoasigmaaa(npt),v2sigmaaa2(npt)
      parameter(tol=1.0d-20)
      
      if(ideriv.eq.0) then
      
      do i=1,npt
      rho = dmax1(0.D0,rhoa1(i))
      if(rho.gt.tol) then
      sigma = dmax1(0.D0,sigmaaa1(i))
      t2 = rho**(1.D0/3.D0)
      t9 = 0.2913644D-2+0.236854225D-3*sigma/(0.6255746320201D7+0.25D0
     &*sigma)
      t11 = rho**2
      t12 = t2**2
      t14 = 1/t12/t11
      t15 = t9**2
      t19 = dlog(0.1587401051968199D1*sigma*t14+dsqrt(1
     &+0.2519842099789745D1*sigma**2*t14**2))
      t20 = t19**2
      t25 = dsqrt(1.D0+0.142866094677138D2*t15*sigma*t14*t20)
      zk(i) = -0.7385587663820224D0*t2*rho*(1.D0+0.1705918482380012D1
     &*t9*sigma*t14/t25)
      else ! rho
      zk(i) = 0.0d0
      endif ! rho
      enddo
      
      else if(ideriv.eq.1) then
      
      do i=1,npt
      rho = dmax1(0.D0,rhoa1(i))
      if(rho.gt.tol) then
      sigma = dmax1(0.D0,sigmaaa1(i))
      t2 = rho**(1.D0/3.D0)
      t3 = t2*rho
      t5 = 0.6255746320201D7+0.25D0*sigma
      t6 = 1/t5
      t9 = 0.2913644D-2+0.236854225D-3*sigma*t6
      t10 = t9*sigma
      t11 = rho**2
      t12 = t2**2
      t14 = 1/t12/t11
      t15 = t9**2
      t16 = t15*sigma
      t19 = dlog(0.1587401051968199D1*sigma*t14+dsqrt(1
     &+0.2519842099789745D1*sigma**2*t14**2))
      t20 = t19**2
      t21 = t14*t20
      t24 = 1.D0+0.142866094677138D2*t16*t21
      t25 = dsqrt(t24)
      t26 = 1/t25
      t27 = t14*t26
      t30 = 1.D0+0.1705918482380012D1*t10*t27
      zk(i) = -0.7385587663820224D0*t3*t30
      t37 = 1/t12/t11/rho
      t43 = t14/t25/t24
      t47 = sigma**2
      t49 = t11**2
      t56 = 1/t2/t49/rho
      t60 = dsqrt(1.D0+0.2519842099789746D1*t47*t56)
      t61 = 1/t60
      vrhoa(i) = -0.9847450218426965D0*t2*t30-0.3692793831910112D0*t3*
     &(-0.9098231906026728D1*t10*t37*t26-0.8529592411900058D0*t10*t43*
     &(-0.7619525049447357D2*t16*t37*t20-0.2419048415798156D3*t15*t47
     &/t2/t49/t11*t19*t61))
      t73 = t5**2
      t77 = 0.9474169D-3*t6-0.236854225D-3*sigma/t73
      vsigmaaa(i) = -0.7385587663820224D0*t3*(0.1705918482380012D1*t77
     &*sigma*t27+0.6823673929520046D1*t9*t14*t26-0.8529592411900058D0
     &*t10*t43*(0.2857321893542759D2*t10*t21*t77+0.5714643787085518D2
     &*t15*t14*t20+0.1814286311848617D3*t16*t56*t19*t61))
      else ! rho
      zk(i) = 0.0d0
      vrhoa(i) = 0.0d0
      vsigmaaa(i) = 0.0d0
      endif ! rho
      enddo
      
      else if(ideriv.eq.2) then
      
      do i=1,npt
      rho = dmax1(0.D0,rhoa1(i))
      if(rho.gt.tol) then
      sigma = dmax1(0.D0,sigmaaa1(i))
      t2 = rho**(1.D0/3.D0)
      t3 = t2*rho
      t5 = 0.6255746320201D7+0.25D0*sigma
      t6 = 1/t5
      t9 = 0.2913644D-2+0.236854225D-3*sigma*t6
      t10 = t9*sigma
      t11 = rho**2
      t12 = t2**2
      t14 = 1/t12/t11
      t15 = t9**2
      t16 = t15*sigma
      t19 = dlog(0.1587401051968199D1*sigma*t14+dsqrt(1
     &+0.2519842099789745D1*sigma**2*t14**2))
      t20 = t19**2
      t21 = t14*t20
      t24 = 1.D0+0.142866094677138D2*t16*t21
      t25 = dsqrt(t24)
      t26 = 1/t25
      t27 = t14*t26
      t30 = 1.D0+0.1705918482380012D1*t10*t27
      zk(i) = -0.7385587663820224D0*t3*t30
      t35 = t11*rho
      t37 = 1/t12/t35
      t38 = t37*t26
      t42 = 1/t25/t24
      t43 = t14*t42
      t44 = t37*t20
      t47 = sigma**2
      t48 = t15*t47
      t49 = t11**2
      t52 = 1/t2/t49/t11
      t56 = 1/t2/t49/rho
      t59 = 1.D0+0.2519842099789746D1*t47*t56
      t60 = dsqrt(t59)
      t61 = 1/t60
      t62 = t52*t19*t61
      t65 = -0.7619525049447357D2*t16*t44-0.2419048415798156D3*t48*t62
      t66 = t43*t65
      t69 = -0.9098231906026728D1*t10*t38-0.8529592411900058D0*t10*t66
      vrhoa(i) = -0.9847450218426965D0*t2*t30-0.3692793831910112D0*t3
     &*t69
      t73 = t5**2
      t74 = 1/t73
      t77 = 0.9474169D-3*t6-0.236854225D-3*sigma*t74
      t78 = t77*sigma
      t81 = t9*t14
      t94 = 0.2857321893542759D2*t10*t21*t77+0.5714643787085518D2*t15
     &*t14*t20+0.1814286311848617D3*t16*t56*t19*t61
      t95 = t43*t94
      t98 = 0.1705918482380012D1*t78*t27+0.6823673929520046D1*t81*t26
     &-0.8529592411900058D0*t10*t95
      vsigmaaa(i) = -0.7385587663820224D0*t3*t98
      t107 = 1/t12/t49
      t111 = t37*t42
      t115 = t24**2
      t117 = 1/t25/t115
      t118 = t14*t117
      t119 = t65**2
      t134 = t15*t47*sigma
      t135 = t49**2
      t136 = t135*t11
      t138 = 1/t59
      t142 = t47**2
      t149 = 1/t60/t59
      v2rhoa2(i) = -0.6564966812284644D0/t12*t30-0.1969490043685393D1
     &*t2*t69-0.3692793831910112D0*t3*(0.6672036731086267D2*t10*t107
     &*t26+0.9098231906026728D1*t10*t111*t65+0.1279438861785009D1*t10
     &*t118*t119-0.8529592411900058D0*t10*t43*(0.5587651702928062D3
     &*t16*t107*t20+0.4354287148436682D4*t48/t2/t49/t35*t19*t61
     &+2048.D0*t134/t136*t138-0.3250997354430873D4*t15*t142/t12/t135
     &/t49*t19*t149))
      t190 = t19*t61
      t191 = t190*t77
      v2rhoasigmaaa(i) = -0.9847450218426965D0*t2*t98
     &-0.3692793831910112D0*t3*(-0.9098231906026728D1*t78*t38
     &-0.3639292762410691D2*t9*t37*t26+0.4549115953013364D1*t10*t111
     &*t94-0.8529592411900058D0*t78*t66-0.3411836964760023D1*t81*t42
     &*t65+0.1279438861785009D1*t10*t14*t117*t65*t94
     &-0.8529592411900058D0*t10*t43*(-0.1523905009889471D3*t10*t44*t77
     &-0.3047810019778943D3*t15*t37*t20-0.2902858098957788D4*t16*t62
     &-0.4838096831596313D3*t9*t47*t52*t191-1536.D0*t48/t135/rho*t138
     &+0.2438248015823154D4*t134/t12/t135/t35*t19*t149))
      t218 = -0.18948338D-2*t74+0.47370845D-3*sigma/t73/t5
      t230 = t94**2
      t234 = t77**2
      v2sigmaaa2(i) = -0.7385587663820224D0*t3*(0.1705918482380012D1
     &*t218*sigma*t27+0.1364734785904009D2*t77*t14*t26
     &-0.1705918482380012D1*t78*t95-0.6823673929520046D1*t81*t42*t94
     &+0.1279438861785009D1*t10*t118*t230-0.8529592411900058D0*t10*t43
     &*(0.2857321893542759D2*t234*sigma*t21+0.2285857514834207D3*t81
     &*t20*t77+0.7257145247394469D3*t10*t56*t191+0.2857321893542759D2
     &*t10*t21*t218+0.1451429049478894D4*t15*t56*t190+1152.D0*t16/t135
     &*t138-0.1828686011867366D4*t48/t12/t136*t19*t149))
      else ! rho
      zk(i) = 0.0d0
      vrhoa(i) = 0.0d0
      v2rhoa2(i) = 0.0d0
      vsigmaaa(i) = 0.0d0
      v2rhoasigmaaa(i) = 0.0d0
      v2sigmaaa2(i) = 0.0d0
      endif ! rho
      enddo
      
      endif ! ideriv
      return
      end

c:X_FT97Bsubrend
c:C_FT97subrstart

c    Generated: Tue Nov  4 11:48:09 GMT 2003

      subroutine uks_c_ft97
     & (ideriv,npt,rhoa1,rhob1,sigmaaa1,sigmabb1,sigmaab1,
     &  zk,vrhoa,vrhob,vsigmaaa,vsigmabb,vsigmaab,
     &  v2rhoa2,v2rhob2,v2rhoab,
     &  v2rhoasigmaaa,v2rhoasigmaab,v2rhoasigmabb,
     &  v2rhobsigmabb,v2rhobsigmaab,v2rhobsigmaaa,
     &  v2sigmaaa2,v2sigmaaaab,v2sigmaaabb,
     &  v2sigmaab2,v2sigmaabbb,v2sigmabb2)
c
c     M. Filatov, and W. Thiel
c     A nonlocal correlation energy density functional from a 
c     Coulomb hole model
c     Int. J. Quant. Chem. 62 (1997) 603-616
c
c     M. Filatov, and W. Thiel
c     A new gradient-corrected exchange-correlation density functional
c     Mol. Phys. 91 (1997) 847-859
c
c
c     CITATION:
c
c     Functionals were obtained from the Density Functional Repository 
c     as developed and distributed by the Quantum Chemistry Group, 
c     CCLRC Daresbury Laboratory, Daresbury, Cheshire, WA4 4AD 
c     United Kingdom. Contact Huub van Dam (h.j.j.vandam@dl.ac.uk) or 
c     Paul Sherwood for further information.
c
c     COPYRIGHT:
c
c     Users may incorporate the source code into software packages and
c     redistribute the source code provided the source code is not
c     changed in anyway and is properly cited in any documentation or
c     publication related to its use.
c
c     ACKNOWLEDGEMENT:
c
c     The source code was generated using Maple 8 through a modified
c     version of the dfauto script published in:
c
c        R. Strange, F.R. Manby, P.J. Knowles
c        Automatic code generation in density functional theory
c        Comp. Phys. Comm. 136 (2001) 310-318.
c
      implicit real*8 (a-h,o-z)
      integer ideriv,npt
      real*8 rhoa1(npt),rhob1(npt)
      real*8 sigmaaa1(npt),sigmabb1(npt),sigmaab1(npt)
      real*8 zk(npt),vrhoa(npt),vrhob(npt)
      real*8 vsigmaaa(npt),vsigmabb(npt),vsigmaab(npt)
      real*8 v2rhoa2(npt),v2rhob2(npt),v2rhoab(npt)
      real*8 v2rhoasigmaaa(npt),v2rhoasigmaab(npt)
      real*8 v2rhoasigmabb(npt),v2rhobsigmabb(npt)
      real*8 v2rhobsigmaab(npt),v2rhobsigmaaa(npt)
      real*8 v2sigmaaa2(npt),v2sigmaaaab(npt),v2sigmaaabb(npt)
      real*8 v2sigmaab2(npt),v2sigmaabbb(npt),v2sigmabb2(npt)
      parameter(tolmin=1.0d-20)
      parameter(tolmax=1.0d+5)
      
      if (ideriv.eq.0) then
      
      do i=1,npt
      rhoa = dmax1(0.D0,rhoa1(i))
      rhob = dmax1(0.D0,rhob1(i))
      rho = rhoa+rhob
      zk(i) = 0.D0
      if(rho.gt.tolmin) then
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigmaab = sigmaab1(i)
      sigmabb = dmax1(0.D0,sigmabb1(i))
      sigma = sigmaaa+sigmabb+2.D0*sigmaab
      if(rhoa.gt.tolmin) then
      t1 = 1/rhoa
      t2 = t1**(1.D0/3.D0)
      trsa = 0.1285808356226127D-1*t2
      t3 = t1**(1.D0/15.D0)
      t4 = t3**2
      t6 = dexp(-0.5416509560827845D0*t4)
      t8 = t1**(1.D0/6.D0)
      t10 = dexp(-0.8579886777788441D0*t8)
      t13 = (0.1247511874D1+0.812904345D0*t6-0.859614445D0*t10)**2
      t14 = sigmaaa**2
      t15 = rhoa**2
      t16 = t15**2
      t18 = rhoa**(1.D0/3.D0)
      t21 = t14/t18/t16/rhoa
      t24 = (0.1D1+0.4473229026479948D-1*t21)**2
      t26 = t18**2
      t36 = dexp(-0.4473229026479948D-1*t21)
      t37 = t36**2
      trcaa = t13*t24/(0.1D1+0.2481823842598833D0*sigmaaa/t26/t15/t2)
     &*t37
      if(trsa.lt.tolmax*trcaa) then
      t2 = 1/rhoa
      t3 = t2**(1.D0/3.D0)
      t4 = t3**2
      t5 = t2**(1.D0/6.D0)
      t9 = (0.7395908974809453D0*t5+0.1075172860312113D1*t3)**2
      t13 = dexp(-0.3848347315591266D0*t4/t9)
      t15 = 0.3141592653589793D1**2
      t18 = 1/0.3141592653589793D1*t2
      t19 = t18**(1.D0/3.D0)
      t21 = t18**(1.D0/15.D0)
      t22 = t21**2
      t24 = dexp(-0.630966299458536D0*t22)
      t26 = t18**(1.D0/6.D0)
      t28 = dexp(-0.1038340679664977D1*t26)
      t31 = (0.1247511874D1+0.812904345D0*t24-0.859614445D0*t28)**2
      t34 = 0.3141592653589793D1**(1.D0/3.D0)
      t37 = sigmaaa**2
      t39 = rhoa**2
      t40 = t39**2
      t42 = rhoa**(1.D0/3.D0)
      t44 = 1/t42/t40/rhoa
      t45 = 1/t34/0.3141592653589793D1*t37*t44
      t48 = (0.1D1+0.2058200272046996D0*t45)**2
      t50 = t34**2
      t53 = t42**2
      t55 = 1/t53/t39
      t63 = dexp(0.4116400544093991D0*t45)
      t67 = expei(-0.1858628590577086D0/t15*t19/t31/t48*(0.1D1
     &+0.3634859066227017D0/t50*sigmaaa*t55/t19)*t63)
      t69 = t2**(1.D0/15.D0)
      t70 = t69**2
      t72 = dexp(-0.5416509560827845D0*t70)
      t75 = dexp(-0.8579886777788441D0*t5)
      t78 = (0.1247511874D1+0.812904345D0*t72-0.859614445D0*t75)**2
      t80 = t3/t78
      t81 = t37*t44
      t84 = (0.1D1+0.4473229026479948D-1*t81)**2
      t85 = 1/t84
      t90 = 0.1D1+0.2481823842598833D0*sigmaaa*t55/t3
      t93 = dexp(0.8946458052959895D-1*t81)
      t95 = t80*t85*t90*t93
      t96 = dsqrt(t95)
      zk(i) = zk(i)+0.1554534543482745D-1*rhoa*t13*(1.D0*t67+(0.6D1
     &+0.4535739597862518D0*t96+0.5143233424904509D-1*t95)/(0.3D1
     &+0.6803609396793777D0*t96+0.7714850137356763D-1*t95)*
     &(0.1285808356226127D-1*t80*t85*t90*t67*t93+0.1D1))
      endif
      endif ! rhoa
      if(rhob.gt.tolmin) then
      t1 = 1/rhob
      t2 = t1**(1.D0/3.D0)
      trsb = 0.1285808356226127D-1*t2
      t3 = t1**(1.D0/15.D0)
      t4 = t3**2
      t6 = dexp(-0.5416509560827845D0*t4)
      t8 = t1**(1.D0/6.D0)
      t10 = dexp(-0.8579886777788441D0*t8)
      t13 = (0.1247511874D1+0.812904345D0*t6-0.859614445D0*t10)**2
      t14 = sigmabb**2
      t15 = rhob**2
      t16 = t15**2
      t18 = rhob**(1.D0/3.D0)
      t21 = t14/t18/t16/rhob
      t24 = (0.1D1+0.4473229026479948D-1*t21)**2
      t26 = t18**2
      t36 = dexp(-0.4473229026479948D-1*t21)
      t37 = t36**2
      trcbb = t13*t24/(0.1D1+0.2481823842598833D0*sigmabb/t26/t15/t2)
     &*t37
      if(trsb.lt.tolmax*trcbb) then
      t2 = 1/rhob
      t3 = t2**(1.D0/3.D0)
      t4 = t3**2
      t5 = t2**(1.D0/6.D0)
      t9 = (0.7395908974809453D0*t5+0.1075172860312113D1*t3)**2
      t13 = dexp(-0.3848347315591266D0*t4/t9)
      t15 = 0.3141592653589793D1**2
      t18 = 1/0.3141592653589793D1*t2
      t19 = t18**(1.D0/3.D0)
      t21 = t18**(1.D0/15.D0)
      t22 = t21**2
      t24 = dexp(-0.630966299458536D0*t22)
      t26 = t18**(1.D0/6.D0)
      t28 = dexp(-0.1038340679664977D1*t26)
      t31 = (0.1247511874D1+0.812904345D0*t24-0.859614445D0*t28)**2
      t34 = 0.3141592653589793D1**(1.D0/3.D0)
      t37 = sigmabb**2
      t39 = rhob**2
      t40 = t39**2
      t42 = rhob**(1.D0/3.D0)
      t44 = 1/t42/t40/rhob
      t45 = 1/t34/0.3141592653589793D1*t37*t44
      t48 = (0.1D1+0.2058200272046996D0*t45)**2
      t50 = t34**2
      t53 = t42**2
      t55 = 1/t53/t39
      t63 = dexp(0.4116400544093991D0*t45)
      t67 = expei(-0.1858628590577086D0/t15*t19/t31/t48*(0.1D1
     &+0.3634859066227017D0/t50*sigmabb*t55/t19)*t63)
      t69 = t2**(1.D0/15.D0)
      t70 = t69**2
      t72 = dexp(-0.5416509560827845D0*t70)
      t75 = dexp(-0.8579886777788441D0*t5)
      t78 = (0.1247511874D1+0.812904345D0*t72-0.859614445D0*t75)**2
      t80 = t3/t78
      t81 = t37*t44
      t84 = (0.1D1+0.4473229026479948D-1*t81)**2
      t85 = 1/t84
      t90 = 0.1D1+0.2481823842598833D0*sigmabb*t55/t3
      t93 = dexp(0.8946458052959895D-1*t81)
      t95 = t80*t85*t90*t93
      t96 = dsqrt(t95)
      zk(i) = zk(i)+0.1554534543482745D-1*rhob*t13*(1.D0*t67+(0.6D1
     &+0.4535739597862518D0*t96+0.5143233424904509D-1*t95)/(0.3D1
     &+0.6803609396793777D0*t96+0.7714850137356763D-1*t95)*
     &(0.1285808356226127D-1*t80*t85*t90*t67*t93+0.1D1))
      endif
      endif ! rhob
      if ((rhoa.gt.tolmin).and.(rhob.gt.tolmin)) then
      t1 = 1/rhob
      t2 = t1**(1.D0/3.D0)
      trsb = 0.1285808356226127D-1*t2
      t3 = t1**(1.D0/15.D0)
      t4 = t3**2
      t5 = t4**2
      t7 = dexp(-0.5683671053580832D-1*t5)
      t10 = (0.942486901D0+0.349064173D0*t7)**2
      t11 = rhob**2
      t12 = rhob**(1.D0/3.D0)
      t13 = t12**2
      t16 = sigmabb/t13/t11
      t18 = sigmabb**2
      t19 = t11**2
      t23 = t18/t12/t19/rhob
      t26 = (0.1D1+0.6936084891727404D-1*t16+0.4389159297635699D-3*t23
     &)**2
      t34 = dexp(-0.4389159297635699D-3*t23)
      t35 = t34**2
      trcab = t10*t26/(0.1D1+0.9505299311782149D-1*t16/t2)*t35
      if(trsb.lt.tolmax*trcab) then
      t2 = 0.3141592653589793D1**2
      t5 = 1/rhob
      t6 = 1/0.3141592653589793D1*t5
      t7 = t6**(1.D0/3.D0)
      t9 = t6**(1.D0/15.D0)
      t10 = t9**2
      t11 = t10**2
      t13 = dexp(-0.7712625328179681D-1*t11)
      t16 = (0.942486901D0+0.349064173D0*t13)**2
      t19 = 0.3141592653589793D1**(1.D0/3.D0)
      t20 = t19**2
      t22 = 1/t20*sigmabb
      t23 = rhob**2
      t24 = rhob**(1.D0/3.D0)
      t25 = t24**2
      t27 = 1/t25/t23
      t32 = sigmabb**2
      t34 = t23**2
      t37 = 1/t24/t34/rhob
      t38 = 1/t19/0.3141592653589793D1*t32*t37
      t41 = (0.1D1+0.1487810599361293D0*t22*t27+0.2019518519390501D-2
     &*t38)**2
      t50 = dexp(0.4039037038781002D-2*t38)
      t54 = expei(-0.1858628590577086D0/t2*t7/t16/t41*(0.1D1
     &+0.1392138426088027D0*t22*t27/t7)*t50)
      t56 = t5**(1.D0/3.D0)
      t57 = t5**(1.D0/15.D0)
      t58 = t57**2
      t59 = t58**2
      t61 = dexp(-0.5683671053580832D-1*t59)
      t64 = (0.942486901D0+0.349064173D0*t61)**2
      t66 = t56/t64
      t67 = sigmabb*t27
      t69 = t32*t37
      t72 = (0.1D1+0.6936084891727404D-1*t67+0.4389159297635699D-3*t69
     &)**2
      t73 = 1/t72
      t77 = 0.1D1+0.9505299311782149D-1*t67/t56
      t80 = dexp(0.8778318595271399D-3*t69)
      t82 = t66*t73*t77*t80
      t83 = dsqrt(t82)
      zk(i) = zk(i)+0.1554534543482745D-1*rhoa*(1.D0*t54+(0.6D1
     &+0.4535739597862518D0*t83+0.5143233424904509D-1*t82)/(0.3D1
     &+0.6803609396793777D0*t83+0.7714850137356763D-1*t82)*
     &(0.1285808356226127D-1*t66*t73*t77*t54*t80+0.1D1))
      endif
      t1 = 1/rhoa
      t2 = t1**(1.D0/3.D0)
      trsa = 0.1285808356226127D-1*t2
      t3 = t1**(1.D0/15.D0)
      t4 = t3**2
      t5 = t4**2
      t7 = dexp(-0.5683671053580832D-1*t5)
      t10 = (0.942486901D0+0.349064173D0*t7)**2
      t11 = rhoa**2
      t12 = rhoa**(1.D0/3.D0)
      t13 = t12**2
      t16 = sigmaaa/t13/t11
      t18 = sigmaaa**2
      t19 = t11**2
      t23 = t18/t12/t19/rhoa
      t26 = (0.1D1+0.6936084891727404D-1*t16+0.4389159297635699D-3*t23
     &)**2
      t34 = dexp(-0.4389159297635699D-3*t23)
      t35 = t34**2
      trcba = t10*t26/(0.1D1+0.9505299311782149D-1*t16/t2)*t35
      if(trsa.lt.tolmax*trcba) then
      t2 = 0.3141592653589793D1**2
      t5 = 1/rhoa
      t6 = 1/0.3141592653589793D1*t5
      t7 = t6**(1.D0/3.D0)
      t9 = t6**(1.D0/15.D0)
      t10 = t9**2
      t11 = t10**2
      t13 = dexp(-0.7712625328179681D-1*t11)
      t16 = (0.942486901D0+0.349064173D0*t13)**2
      t19 = 0.3141592653589793D1**(1.D0/3.D0)
      t20 = t19**2
      t22 = 1/t20*sigmaaa
      t23 = rhoa**2
      t24 = rhoa**(1.D0/3.D0)
      t25 = t24**2
      t27 = 1/t25/t23
      t32 = sigmaaa**2
      t34 = t23**2
      t37 = 1/t24/t34/rhoa
      t38 = 1/t19/0.3141592653589793D1*t32*t37
      t41 = (0.1D1+0.1487810599361293D0*t22*t27+0.2019518519390501D-2
     &*t38)**2
      t50 = dexp(0.4039037038781002D-2*t38)
      t54 = expei(-0.1858628590577086D0/t2*t7/t16/t41*(0.1D1
     &+0.1392138426088027D0*t22*t27/t7)*t50)
      t56 = t5**(1.D0/3.D0)
      t57 = t5**(1.D0/15.D0)
      t58 = t57**2
      t59 = t58**2
      t61 = dexp(-0.5683671053580832D-1*t59)
      t64 = (0.942486901D0+0.349064173D0*t61)**2
      t66 = t56/t64
      t67 = sigmaaa*t27
      t69 = t32*t37
      t72 = (0.1D1+0.6936084891727404D-1*t67+0.4389159297635699D-3*t69
     &)**2
      t73 = 1/t72
      t77 = 0.1D1+0.9505299311782149D-1*t67/t56
      t80 = dexp(0.8778318595271399D-3*t69)
      t82 = t66*t73*t77*t80
      t83 = dsqrt(t82)
      zk(i) = zk(i)+0.1554534543482745D-1*rhob*(1.D0*t54+(0.6D1
     &+0.4535739597862518D0*t83+0.5143233424904509D-1*t82)/(0.3D1
     &+0.6803609396793777D0*t83+0.7714850137356763D-1*t82)*
     &(0.1285808356226127D-1*t66*t73*t77*t54*t80+0.1D1))
      endif
      endif ! rhoa,rhob
      endif ! rho
      enddo
      
      else if(ideriv.eq.1) then
      
      do i=1,npt
      rhoa = dmax1(0.D0,rhoa1(i))
      rhob = dmax1(0.D0,rhob1(i))
      rho = rhoa+rhob
      zk(i) = 0.D0
      vrhoa(i) = 0.D0
      vrhob(i) = 0.D0
      vsigmaaa(i) = 0.D0
      vsigmaab(i) = 0.D0
      vsigmabb(i) = 0.D0
      if(rho.gt.tolmin) then
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigmaab = sigmaab1(i)
      sigmabb = dmax1(0.D0,sigmabb1(i))
      sigma = sigmaaa+sigmabb+2.D0*sigmaab
      if(rhoa.gt.tolmin) then
      t1 = 1/rhoa
      t2 = t1**(1.D0/3.D0)
      trsa = 0.1285808356226127D-1*t2
      t3 = t1**(1.D0/15.D0)
      t4 = t3**2
      t6 = dexp(-0.5416509560827845D0*t4)
      t8 = t1**(1.D0/6.D0)
      t10 = dexp(-0.8579886777788441D0*t8)
      t13 = (0.1247511874D1+0.812904345D0*t6-0.859614445D0*t10)**2
      t14 = sigmaaa**2
      t15 = rhoa**2
      t16 = t15**2
      t18 = rhoa**(1.D0/3.D0)
      t21 = t14/t18/t16/rhoa
      t24 = (0.1D1+0.4473229026479948D-1*t21)**2
      t26 = t18**2
      t36 = dexp(-0.4473229026479948D-1*t21)
      t37 = t36**2
      trcaa = t13*t24/(0.1D1+0.2481823842598833D0*sigmaaa/t26/t15/t2)
     &*t37
      if(trsa.lt.tolmax*trcaa) then
      t2 = 1/rhoa
      t3 = t2**(1.D0/3.D0)
      t4 = t3**2
      t5 = t2**(1.D0/6.D0)
      t8 = 0.7395908974809453D0*t5+0.1075172860312113D1*t3
      t9 = t8**2
      t10 = 1/t9
      t13 = dexp(-0.3848347315591266D0*t4*t10)
      t14 = rhoa*t13
      t15 = 0.3141592653589793D1**2
      t18 = 1/0.3141592653589793D1*t2
      t19 = t18**(1.D0/3.D0)
      t21 = t18**(1.D0/15.D0)
      t22 = t21**2
      t24 = dexp(-0.630966299458536D0*t22)
      t26 = t18**(1.D0/6.D0)
      t28 = dexp(-0.1038340679664977D1*t26)
      t31 = (0.1247511874D1+0.812904345D0*t24-0.859614445D0*t28)**2
      t34 = 0.3141592653589793D1**(1.D0/3.D0)
      t37 = sigmaaa**2
      t39 = rhoa**2
      t40 = t39**2
      t42 = rhoa**(1.D0/3.D0)
      t44 = 1/t42/t40/rhoa
      t45 = 1/t34/0.3141592653589793D1*t37*t44
      t48 = (0.1D1+0.2058200272046996D0*t45)**2
      t50 = t34**2
      t53 = t42**2
      t55 = 1/t53/t39
      t63 = dexp(0.4116400544093991D0*t45)
      t67 = expei(-0.1858628590577086D0/t15*t19/t31/t48*(0.1D1
     &+0.3634859066227017D0/t50*sigmaaa*t55/t19)*t63)
      t69 = t2**(1.D0/15.D0)
      t70 = t69**2
      t72 = dexp(-0.5416509560827845D0*t70)
      t75 = dexp(-0.8579886777788441D0*t5)
      t77 = 0.1247511874D1+0.812904345D0*t72-0.859614445D0*t75
      t78 = t77**2
      t79 = 1/t78
      t80 = t3*t79
      t81 = t37*t44
      t83 = 0.1D1+0.4473229026479948D-1*t81
      t84 = t83**2
      t85 = 1/t84
      t87 = 1/t3
      t90 = 0.1D1+0.2481823842598833D0*sigmaaa*t55*t87
      t91 = t85*t90
      t93 = dexp(0.8946458052959895D-1*t81)
      t95 = t80*t91*t93
      t96 = dsqrt(t95)
      t99 = 0.6D1+0.4535739597862518D0*t96+0.5143233424904509D-1*t95
      t102 = 0.3D1+0.6803609396793777D0*t96+0.7714850137356763D-1*t95
      t103 = 1/t102
      t104 = t99*t103
      t105 = t80*t85
      t106 = t90*t67
      t110 = 0.1285808356226127D-1*t105*t106*t93+0.1D1
      t112 = 1.D0*t67+t104*t110
      zk(i) = zk(i)+0.1554534543482745D-1*t14*t112
      t116 = t13*t112
      t119 = 1/t39
      t125 = t5**2
      t126 = t125**2
      t129 = 1/t126/t5*t119
      t131 = 1/t4
      t142 = t131*t79*t85
      t143 = t90*t93
      t145 = t142*t143*t119
      t146 = 0.4286027854087091D-2*t145
      t150 = t3/t78/t77*t85
      t151 = t70**2
      t153 = t151**2
      t161 = 0.5870805542307996D-1/t153/t151/t69*t119*t72
     &-0.1229232435108575D0*t129*t75
      t163 = t150*t143*t161
      t164 = 0.2571616712452254D-1*t163
      t166 = 1/t84/t83
      t167 = t80*t166
      t170 = 1/t42/t40/t39
      t172 = t143*t37*t170
      t173 = t167*t172
      t174 = 0.613516294566526D-2*t173
      t188 = -0.6618196913596887D0*sigmaaa/t53/t39/rhoa*t87
     &+0.8272746141996109D-1*sigmaaa/t53/t40/t3/t2
      t191 = t80*t85*t188*t93
      t192 = 0.1285808356226127D-1*t191
      t193 = t105*t172
      t194 = 0.613516294566526D-2*t193
      t195 = -t146-t164+t174+t192-t194
      t204 = dexp(-0.8946458052959895D-1*t81)
      t205 = t84/t90*t204
      t208 = 1/t96
      t215 = t208*(-0.2306016315371713D-1*t145-0.1383609789223028D0
     &*t163+0.3300908517586267D-1*t173+0.6918048946115138D-1*t191
     &-0.3300908517586267D-1*t193)
      t225 = t102**2
      t227 = t99/t225
      t247 = t80*t166*t90
      t250 = t67*t37*t170*t93
      t257 = t80*t91
      t261 = t67*t93
      t265 = -0.4286027854087091D-2*t142*t106*t119*t93
     &-0.2571616712452254D-1*t150*t106*t161*t93+0.613516294566526D-2
     &*t247*t250+0.1285808356226127D-1*t105*t188*t67*t93
     &-0.613516294566526D-2*t257*t250+0.1285808356226127D-1*t105*t90
     &*t195*t261-t146-t164+t174+t192-t194
      vrhoa(i) = vrhoa(i)+0.1554534543482745D-1*t116
     &+0.1554534543482745D-1*rhoa*(0.2565564877060844D0*t87*t10*t119
     &+0.7696694631182532D0*t4/t9/t8*(-0.1232651495801576D0*t129
     &-0.358390953437371D0*t131*t119))*t116+0.1554534543482745D-1*t14*
     &(1.D0*t195*t67+0.7777208750882749D2*t195*t87*t78*t205+
     &(0.3278192763011299D1*t215-0.1714411141634836D-1*t145
     &-0.1028646684980902D0*t163+0.2454065178266104D-1*t173
     &+0.5143233424904509D-1*t191-0.2454065178266104D-1*t193)*t103
     &*t110-1.D0*t227*t110*(0.4917289144516948D1*t215
     &-0.2571616712452254D-1*t145-0.1542970027471353D0*t163
     &+0.3681097767399156D-1*t173+0.7714850137356763D-1*t191
     &-0.3681097767399156D-1*t193)+t104*t265)
      t272 = t143*sigmaaa*t44
      t273 = t167*t272
      t274 = 0.2300686104624472D-2*t273
      t275 = t79*t85
      t277 = t275*t55*t93
      t278 = 0.3191149835494816D-2*t277
      t279 = t105*t272
      t280 = 0.2300686104624472D-2*t279
      t281 = -t274+t278+t280
      t293 = t208*(-0.123784069409485D-1*t273+0.1716937881873428D-1
     &*t277+0.123784069409485D-1*t279)
      t311 = t67*sigmaaa*t44*t93
      vsigmaaa(i) = vsigmaaa(i)+0.1554534543482745D-1*t14*(1.D0*t281
     &*t67+0.7777208750882749D2*t281*t87*t78*t205+
     &(0.3278192763011299D1*t293-0.9202744418497889D-2*t273
     &+0.1276459934197926D-1*t277+0.9202744418497889D-2*t279)*t103
     &*t110-1.D0*t227*t110*(0.4917289144516948D1*t293
     &-0.1380411662774683D-1*t273+0.1914689901296889D-1*t277
     &+0.1380411662774683D-1*t279)+t104*(-0.2300686104624472D-2*t247
     &*t311+0.3191149835494816D-2*t275*t55*t67*t93
     &+0.2300686104624472D-2*t257*t311+0.1285808356226127D-1*t105*t90
     &*t281*t261-t274+t278+t280))
      endif
      endif ! rhoa
      if(rhob.gt.tolmin) then
      t1 = 1/rhob
      t2 = t1**(1.D0/3.D0)
      trsb = 0.1285808356226127D-1*t2
      t3 = t1**(1.D0/15.D0)
      t4 = t3**2
      t6 = dexp(-0.5416509560827845D0*t4)
      t8 = t1**(1.D0/6.D0)
      t10 = dexp(-0.8579886777788441D0*t8)
      t13 = (0.1247511874D1+0.812904345D0*t6-0.859614445D0*t10)**2
      t14 = sigmabb**2
      t15 = rhob**2
      t16 = t15**2
      t18 = rhob**(1.D0/3.D0)
      t21 = t14/t18/t16/rhob
      t24 = (0.1D1+0.4473229026479948D-1*t21)**2
      t26 = t18**2
      t36 = dexp(-0.4473229026479948D-1*t21)
      t37 = t36**2
      trcbb = t13*t24/(0.1D1+0.2481823842598833D0*sigmabb/t26/t15/t2)
     &*t37
      if(trsb.lt.tolmax*trcbb) then
      t2 = 1/rhob
      t3 = t2**(1.D0/3.D0)
      t4 = t3**2
      t5 = t2**(1.D0/6.D0)
      t8 = 0.7395908974809453D0*t5+0.1075172860312113D1*t3
      t9 = t8**2
      t10 = 1/t9
      t13 = dexp(-0.3848347315591266D0*t4*t10)
      t14 = rhob*t13
      t15 = 0.3141592653589793D1**2
      t18 = 1/0.3141592653589793D1*t2
      t19 = t18**(1.D0/3.D0)
      t21 = t18**(1.D0/15.D0)
      t22 = t21**2
      t24 = dexp(-0.630966299458536D0*t22)
      t26 = t18**(1.D0/6.D0)
      t28 = dexp(-0.1038340679664977D1*t26)
      t31 = (0.1247511874D1+0.812904345D0*t24-0.859614445D0*t28)**2
      t34 = 0.3141592653589793D1**(1.D0/3.D0)
      t37 = sigmabb**2
      t39 = rhob**2
      t40 = t39**2
      t42 = rhob**(1.D0/3.D0)
      t44 = 1/t42/t40/rhob
      t45 = 1/t34/0.3141592653589793D1*t37*t44
      t48 = (0.1D1+0.2058200272046996D0*t45)**2
      t50 = t34**2
      t53 = t42**2
      t55 = 1/t53/t39
      t63 = dexp(0.4116400544093991D0*t45)
      t67 = expei(-0.1858628590577086D0/t15*t19/t31/t48*(0.1D1
     &+0.3634859066227017D0/t50*sigmabb*t55/t19)*t63)
      t69 = t2**(1.D0/15.D0)
      t70 = t69**2
      t72 = dexp(-0.5416509560827845D0*t70)
      t75 = dexp(-0.8579886777788441D0*t5)
      t77 = 0.1247511874D1+0.812904345D0*t72-0.859614445D0*t75
      t78 = t77**2
      t79 = 1/t78
      t80 = t3*t79
      t81 = t37*t44
      t83 = 0.1D1+0.4473229026479948D-1*t81
      t84 = t83**2
      t85 = 1/t84
      t87 = 1/t3
      t90 = 0.1D1+0.2481823842598833D0*sigmabb*t55*t87
      t91 = t85*t90
      t93 = dexp(0.8946458052959895D-1*t81)
      t95 = t80*t91*t93
      t96 = dsqrt(t95)
      t99 = 0.6D1+0.4535739597862518D0*t96+0.5143233424904509D-1*t95
      t102 = 0.3D1+0.6803609396793777D0*t96+0.7714850137356763D-1*t95
      t103 = 1/t102
      t104 = t99*t103
      t105 = t80*t85
      t106 = t90*t67
      t110 = 0.1285808356226127D-1*t105*t106*t93+0.1D1
      t112 = 1.D0*t67+t104*t110
      zk(i) = zk(i)+0.1554534543482745D-1*t14*t112
      t116 = t13*t112
      t119 = 1/t39
      t125 = t5**2
      t126 = t125**2
      t129 = 1/t126/t5*t119
      t131 = 1/t4
      t142 = t131*t79*t85
      t143 = t90*t93
      t145 = t142*t143*t119
      t146 = 0.4286027854087091D-2*t145
      t150 = t3/t78/t77*t85
      t151 = t70**2
      t153 = t151**2
      t161 = 0.5870805542307996D-1/t153/t151/t69*t119*t72
     &-0.1229232435108575D0*t129*t75
      t163 = t150*t143*t161
      t164 = 0.2571616712452254D-1*t163
      t166 = 1/t84/t83
      t167 = t80*t166
      t170 = 1/t42/t40/t39
      t172 = t143*t37*t170
      t173 = t167*t172
      t174 = 0.613516294566526D-2*t173
      t188 = -0.6618196913596887D0*sigmabb/t53/t39/rhob*t87
     &+0.8272746141996109D-1*sigmabb/t53/t40/t3/t2
      t191 = t80*t85*t188*t93
      t192 = 0.1285808356226127D-1*t191
      t193 = t105*t172
      t194 = 0.613516294566526D-2*t193
      t195 = -t146-t164+t174+t192-t194
      t204 = dexp(-0.8946458052959895D-1*t81)
      t205 = t84/t90*t204
      t208 = 1/t96
      t215 = t208*(-0.2306016315371713D-1*t145-0.1383609789223028D0
     &*t163+0.3300908517586267D-1*t173+0.6918048946115138D-1*t191
     &-0.3300908517586267D-1*t193)
      t225 = t102**2
      t227 = t99/t225
      t247 = t80*t166*t90
      t250 = t67*t37*t170*t93
      t257 = t80*t91
      t261 = t67*t93
      t265 = -0.4286027854087091D-2*t142*t106*t119*t93
     &-0.2571616712452254D-1*t150*t106*t161*t93+0.613516294566526D-2
     &*t247*t250+0.1285808356226127D-1*t105*t188*t67*t93
     &-0.613516294566526D-2*t257*t250+0.1285808356226127D-1*t105*t90
     &*t195*t261-t146-t164+t174+t192-t194
      vrhob(i) = vrhob(i)+0.1554534543482745D-1*t116
     &+0.1554534543482745D-1*rhob*(0.2565564877060844D0*t87*t10*t119
     &+0.7696694631182532D0*t4/t9/t8*(-0.1232651495801576D0*t129
     &-0.358390953437371D0*t131*t119))*t116+0.1554534543482745D-1*t14*
     &(1.D0*t195*t67+0.7777208750882749D2*t195*t87*t78*t205+
     &(0.3278192763011299D1*t215-0.1714411141634836D-1*t145
     &-0.1028646684980902D0*t163+0.2454065178266104D-1*t173
     &+0.5143233424904509D-1*t191-0.2454065178266104D-1*t193)*t103
     &*t110-1.D0*t227*t110*(0.4917289144516948D1*t215
     &-0.2571616712452254D-1*t145-0.1542970027471353D0*t163
     &+0.3681097767399156D-1*t173+0.7714850137356763D-1*t191
     &-0.3681097767399156D-1*t193)+t104*t265)
      t272 = t143*sigmabb*t44
      t273 = t167*t272
      t274 = 0.2300686104624472D-2*t273
      t275 = t79*t85
      t277 = t275*t55*t93
      t278 = 0.3191149835494816D-2*t277
      t279 = t105*t272
      t280 = 0.2300686104624472D-2*t279
      t281 = -t274+t278+t280
      t293 = t208*(-0.123784069409485D-1*t273+0.1716937881873428D-1
     &*t277+0.123784069409485D-1*t279)
      t311 = t67*sigmabb*t44*t93
      vsigmabb(i) = vsigmabb(i)+0.1554534543482745D-1*t14*(1.D0*t281
     &*t67+0.7777208750882749D2*t281*t87*t78*t205+
     &(0.3278192763011299D1*t293-0.9202744418497889D-2*t273
     &+0.1276459934197926D-1*t277+0.9202744418497889D-2*t279)*t103
     &*t110-1.D0*t227*t110*(0.4917289144516948D1*t293
     &-0.1380411662774683D-1*t273+0.1914689901296889D-1*t277
     &+0.1380411662774683D-1*t279)+t104*(-0.2300686104624472D-2*t247
     &*t311+0.3191149835494816D-2*t275*t55*t67*t93
     &+0.2300686104624472D-2*t257*t311+0.1285808356226127D-1*t105*t90
     &*t281*t261-t274+t278+t280))
      endif
      endif ! rhob
      if ((rhoa.gt.tolmin).and.(rhob.gt.tolmin)) then
      t1 = 1/rhob
      t2 = t1**(1.D0/3.D0)
      trsb = 0.1285808356226127D-1*t2
      t3 = t1**(1.D0/15.D0)
      t4 = t3**2
      t5 = t4**2
      t7 = dexp(-0.5683671053580832D-1*t5)
      t10 = (0.942486901D0+0.349064173D0*t7)**2
      t11 = rhob**2
      t12 = rhob**(1.D0/3.D0)
      t13 = t12**2
      t16 = sigmabb/t13/t11
      t18 = sigmabb**2
      t19 = t11**2
      t23 = t18/t12/t19/rhob
      t26 = (0.1D1+0.6936084891727404D-1*t16+0.4389159297635699D-3*t23
     &)**2
      t34 = dexp(-0.4389159297635699D-3*t23)
      t35 = t34**2
      trcab = t10*t26/(0.1D1+0.9505299311782149D-1*t16/t2)*t35
      if(trsb.lt.tolmax*trcab) then
      t2 = 0.3141592653589793D1**2
      t5 = 1/rhob
      t6 = 1/0.3141592653589793D1*t5
      t7 = t6**(1.D0/3.D0)
      t9 = t6**(1.D0/15.D0)
      t10 = t9**2
      t11 = t10**2
      t13 = dexp(-0.7712625328179681D-1*t11)
      t16 = (0.942486901D0+0.349064173D0*t13)**2
      t19 = 0.3141592653589793D1**(1.D0/3.D0)
      t20 = t19**2
      t22 = 1/t20*sigmabb
      t23 = rhob**2
      t24 = rhob**(1.D0/3.D0)
      t25 = t24**2
      t27 = 1/t25/t23
      t32 = sigmabb**2
      t34 = t23**2
      t37 = 1/t24/t34/rhob
      t38 = 1/t19/0.3141592653589793D1*t32*t37
      t41 = (0.1D1+0.1487810599361293D0*t22*t27+0.2019518519390501D-2
     &*t38)**2
      t50 = dexp(0.4039037038781002D-2*t38)
      t54 = expei(-0.1858628590577086D0/t2*t7/t16/t41*(0.1D1
     &+0.1392138426088027D0*t22*t27/t7)*t50)
      t56 = t5**(1.D0/3.D0)
      t57 = t5**(1.D0/15.D0)
      t58 = t57**2
      t59 = t58**2
      t60 = 0.5683671053580832D-1*t59
      t61 = dexp(-t60)
      t63 = 0.942486901D0+0.349064173D0*t61
      t64 = t63**2
      t65 = 1/t64
      t66 = t56*t65
      t67 = sigmabb*t27
      t69 = t32*t37
      t71 = 0.1D1+0.6936084891727404D-1*t67+0.4389159297635699D-3*t69
      t72 = t71**2
      t73 = 1/t72
      t74 = 1/t56
      t77 = 0.1D1+0.9505299311782149D-1*t67*t74
      t78 = t73*t77
      t79 = 0.8778318595271399D-3*t69
      t80 = dexp(t79)
      t82 = t66*t78*t80
      t83 = dsqrt(t82)
      t86 = 0.6D1+0.4535739597862518D0*t83+0.5143233424904509D-1*t82
      t89 = 0.3D1+0.6803609396793777D0*t83+0.7714850137356763D-1*t82
      t90 = 1/t89
      t91 = t86*t90
      t92 = t66*t73
      t93 = t77*t54
      t97 = 0.1285808356226127D-1*t92*t93*t80+0.1D1
      t98 = t91*t97
      zk(i) = zk(i)+0.1554534543482745D-1*rhoa*(1.D0*t54+t98)
      vrhoa(i) = vrhoa(i)+0.1554534543482745D-1*t54
     &+0.1554534543482745D-1*t98
      t106 = t56**2
      t109 = 1/t106*t65*t73
      t110 = t77*t80
      t111 = 1/t23
      t113 = t109*t110*t111
      t114 = 0.4286027854087091D-2*t113
      t115 = t5**(1.D0/5.D0)
      t116 = t115**2
      t121 = 1/t116/t64/t63*t73
      t124 = dexp(t79-t60)
      t126 = t121*t77*t111*t124
      t127 = 0.1360533322067624D-3*t126
      t130 = t66/t72/t71
      t134 = sigmabb/t25/t23/rhob
      t138 = 1/t24/t34/t23
      t139 = t32*t138
      t141 = -0.1849622637793975D0*t134-0.234088495873904D-2*t139
      t143 = t130*t110*t141
      t144 = 0.2571616712452254D-1*t143
      t154 = -0.2534746483141907D0*t134*t74+0.3168433103927383D-1
     &*sigmabb/t25/t34/t56/t5
      t157 = t66*t73*t154*t80
      t158 = 0.1285808356226127D-1*t157
      t160 = t92*t110*t139
      t161 = 0.601985888182142D-4*t160
      t162 = -t114-t127-t144+t158-t161
      t171 = dexp(-0.8778318595271399D-3*t69)
      t172 = t72/t77*t171
      t175 = 1/t83
      t182 = t175*(-0.2306016315371713D-1*t113-0.7320092507805402D-3
     &*t126-0.1383609789223028D0*t143+0.6918048946115138D-1*t157
     &-0.3238871344356278D-3*t160)
      t192 = t89**2
      t194 = t86/t192
      t221 = t66*t78
      t228 = t54*t80
      t232 = -0.4286027854087091D-2*t109*t93*t111*t80
     &-0.1360533322067624D-3*t121*t93*t111*t124-0.2571616712452254D-1
     &*t130*t93*t141*t80+0.1285808356226127D-1*t92*t154*t54*t80
     &-0.601985888182142D-4*t221*t54*t32*t138*t80
     &+0.1285808356226127D-1*t92*t77*t162*t228-t114-t127-t144+t158-t161
      vrhob(i) = vrhob(i)+0.1554534543482745D-1*rhoa*(1.D0*t162*t54
     &+0.7777208750882749D2*t162*t74*t64*t172+(0.3278192763011299D1
     &*t182-0.1714411141634836D-1*t113-0.5442133288270495D-3*t126
     &-0.1028646684980902D0*t143+0.5143233424904509D-1*t157
     &-0.2407943552728568D-3*t160)*t90*t97-1.D0*t194*t97*
     &(0.4917289144516948D1*t182-0.2571616712452254D-1*t113
     &-0.8163199932405743D-3*t126-0.1542970027471353D0*t143
     &+0.7714850137356763D-1*t157-0.3611915329092852D-3*t160)+t91*t232)
      t239 = sigmabb*t37
      t241 = 0.6936084891727404D-1*t27+0.8778318595271399D-3*t239
      t243 = t130*t110*t241
      t244 = 0.2571616712452254D-1*t243
      t245 = t65*t73
      t247 = t245*t27*t80
      t248 = 0.1222199328351994D-2*t247
      t250 = t92*t110*t239
      t251 = 0.2257447080683033D-4*t250
      t252 = -t244+t248+t251
      t264 = t175*(-0.1383609789223028D0*t243+0.6575812588638345D-2
     &*t247+0.1214576754133604D-3*t250)
      vsigmabb(i) = vsigmabb(i)+0.1554534543482745D-1*rhoa*(1.D0*t252
     &*t54+0.7777208750882749D2*t252*t74*t64*t172+
     &(0.3278192763011299D1*t264-0.1028646684980902D0*t243
     &+0.4888797313407978D-2*t247+0.9029788322732131D-4*t250)*t90*t97
     &-1.D0*t194*t97*(0.4917289144516948D1*t264-0.1542970027471353D0
     &*t243+0.7333195970111966D-2*t247+0.135446824840982D-3*t250)+t91*
     &(-0.2571616712452254D-1*t130*t93*t241*t80+0.1222199328351994D-2
     &*t245*t27*t54*t80+0.2257447080683033D-4*t221*t54*sigmabb*t37*t80
     &+0.1285808356226127D-1*t92*t77*t252*t228-t244+t248+t251))
      endif
      t1 = 1/rhoa
      t2 = t1**(1.D0/3.D0)
      trsa = 0.1285808356226127D-1*t2
      t3 = t1**(1.D0/15.D0)
      t4 = t3**2
      t5 = t4**2
      t7 = dexp(-0.5683671053580832D-1*t5)
      t10 = (0.942486901D0+0.349064173D0*t7)**2
      t11 = rhoa**2
      t12 = rhoa**(1.D0/3.D0)
      t13 = t12**2
      t16 = sigmaaa/t13/t11
      t18 = sigmaaa**2
      t19 = t11**2
      t23 = t18/t12/t19/rhoa
      t26 = (0.1D1+0.6936084891727404D-1*t16+0.4389159297635699D-3*t23
     &)**2
      t34 = dexp(-0.4389159297635699D-3*t23)
      t35 = t34**2
      trcba = t10*t26/(0.1D1+0.9505299311782149D-1*t16/t2)*t35
      if(trsa.lt.tolmax*trcba) then
      t2 = 0.3141592653589793D1**2
      t5 = 1/rhoa
      t6 = 1/0.3141592653589793D1*t5
      t7 = t6**(1.D0/3.D0)
      t9 = t6**(1.D0/15.D0)
      t10 = t9**2
      t11 = t10**2
      t13 = dexp(-0.7712625328179681D-1*t11)
      t16 = (0.942486901D0+0.349064173D0*t13)**2
      t19 = 0.3141592653589793D1**(1.D0/3.D0)
      t20 = t19**2
      t22 = 1/t20*sigmaaa
      t23 = rhoa**2
      t24 = rhoa**(1.D0/3.D0)
      t25 = t24**2
      t27 = 1/t25/t23
      t32 = sigmaaa**2
      t34 = t23**2
      t37 = 1/t24/t34/rhoa
      t38 = 1/t19/0.3141592653589793D1*t32*t37
      t41 = (0.1D1+0.1487810599361293D0*t22*t27+0.2019518519390501D-2
     &*t38)**2
      t50 = dexp(0.4039037038781002D-2*t38)
      t54 = expei(-0.1858628590577086D0/t2*t7/t16/t41*(0.1D1
     &+0.1392138426088027D0*t22*t27/t7)*t50)
      t56 = t5**(1.D0/3.D0)
      t57 = t5**(1.D0/15.D0)
      t58 = t57**2
      t59 = t58**2
      t60 = 0.5683671053580832D-1*t59
      t61 = dexp(-t60)
      t63 = 0.942486901D0+0.349064173D0*t61
      t64 = t63**2
      t65 = 1/t64
      t66 = t56*t65
      t67 = sigmaaa*t27
      t69 = t32*t37
      t71 = 0.1D1+0.6936084891727404D-1*t67+0.4389159297635699D-3*t69
      t72 = t71**2
      t73 = 1/t72
      t74 = 1/t56
      t77 = 0.1D1+0.9505299311782149D-1*t67*t74
      t78 = t73*t77
      t79 = 0.8778318595271399D-3*t69
      t80 = dexp(t79)
      t82 = t66*t78*t80
      t83 = dsqrt(t82)
      t86 = 0.6D1+0.4535739597862518D0*t83+0.5143233424904509D-1*t82
      t89 = 0.3D1+0.6803609396793777D0*t83+0.7714850137356763D-1*t82
      t90 = 1/t89
      t91 = t86*t90
      t92 = t66*t73
      t93 = t77*t54
      t97 = 0.1285808356226127D-1*t92*t93*t80+0.1D1
      t98 = t91*t97
      zk(i) = zk(i)+0.1554534543482745D-1*rhob*(1.D0*t54+t98)
      t103 = t56**2
      t106 = 1/t103*t65*t73
      t107 = t77*t80
      t108 = 1/t23
      t110 = t106*t107*t108
      t111 = 0.4286027854087091D-2*t110
      t112 = t5**(1.D0/5.D0)
      t113 = t112**2
      t118 = 1/t113/t64/t63*t73
      t121 = dexp(t79-t60)
      t123 = t118*t77*t108*t121
      t124 = 0.1360533322067624D-3*t123
      t127 = t66/t72/t71
      t131 = sigmaaa/t25/t23/rhoa
      t135 = 1/t24/t34/t23
      t136 = t32*t135
      t138 = -0.1849622637793975D0*t131-0.234088495873904D-2*t136
      t140 = t127*t107*t138
      t141 = 0.2571616712452254D-1*t140
      t151 = -0.2534746483141907D0*t131*t74+0.3168433103927383D-1
     &*sigmaaa/t25/t34/t56/t5
      t154 = t66*t73*t151*t80
      t155 = 0.1285808356226127D-1*t154
      t157 = t92*t107*t136
      t158 = 0.601985888182142D-4*t157
      t159 = -t111-t124-t141+t155-t158
      t168 = dexp(-0.8778318595271399D-3*t69)
      t169 = t72/t77*t168
      t172 = 1/t83
      t179 = t172*(-0.2306016315371713D-1*t110-0.7320092507805402D-3
     &*t123-0.1383609789223028D0*t140+0.6918048946115138D-1*t154
     &-0.3238871344356278D-3*t157)
      t189 = t89**2
      t191 = t86/t189
      t218 = t66*t78
      t225 = t54*t80
      t229 = -0.4286027854087091D-2*t106*t93*t108*t80
     &-0.1360533322067624D-3*t118*t93*t108*t121-0.2571616712452254D-1
     &*t127*t93*t138*t80+0.1285808356226127D-1*t92*t151*t54*t80
     &-0.601985888182142D-4*t218*t54*t32*t135*t80
     &+0.1285808356226127D-1*t92*t77*t159*t225-t111-t124-t141+t155-t158
      vrhoa(i) = vrhoa(i)+0.1554534543482745D-1*rhob*(1.D0*t159*t54
     &+0.7777208750882749D2*t159*t74*t64*t169+(0.3278192763011299D1
     &*t179-0.1714411141634836D-1*t110-0.5442133288270495D-3*t123
     &-0.1028646684980902D0*t140+0.5143233424904509D-1*t154
     &-0.2407943552728568D-3*t157)*t90*t97-1.D0*t191*t97*
     &(0.4917289144516948D1*t179-0.2571616712452254D-1*t110
     &-0.8163199932405743D-3*t123-0.1542970027471353D0*t140
     &+0.7714850137356763D-1*t154-0.3611915329092852D-3*t157)+t91*t229)
      vrhob(i) = vrhob(i)+0.1554534543482745D-1*t54
     &+0.1554534543482745D-1*t98
      t239 = sigmaaa*t37
      t241 = 0.6936084891727404D-1*t27+0.8778318595271399D-3*t239
      t243 = t127*t107*t241
      t244 = 0.2571616712452254D-1*t243
      t245 = t65*t73
      t247 = t245*t27*t80
      t248 = 0.1222199328351994D-2*t247
      t250 = t92*t107*t239
      t251 = 0.2257447080683033D-4*t250
      t252 = -t244+t248+t251
      t264 = t172*(-0.1383609789223028D0*t243+0.6575812588638345D-2
     &*t247+0.1214576754133604D-3*t250)
      vsigmaaa(i) = vsigmaaa(i)+0.1554534543482745D-1*rhob*(1.D0*t252
     &*t54+0.7777208750882749D2*t252*t74*t64*t169+
     &(0.3278192763011299D1*t264-0.1028646684980902D0*t243
     &+0.4888797313407978D-2*t247+0.9029788322732131D-4*t250)*t90*t97
     &-1.D0*t191*t97*(0.4917289144516948D1*t264-0.1542970027471353D0
     &*t243+0.7333195970111966D-2*t247+0.135446824840982D-3*t250)+t91*
     &(-0.2571616712452254D-1*t127*t93*t241*t80+0.1222199328351994D-2
     &*t245*t27*t54*t80+0.2257447080683033D-4*t218*t54*sigmaaa*t37*t80
     &+0.1285808356226127D-1*t92*t77*t252*t225-t244+t248+t251))
      endif
      endif ! rhoa,rhob
      endif ! rho
      enddo
      
      else if(ideriv.eq.2) then
      
      do i=1,npt
      rhoa = dmax1(0.D0,rhoa1(i))
      rhob = dmax1(0.D0,rhob1(i))
      rho = rhoa+rhob
      zk(i) = 0.D0
      vrhoa(i) = 0.D0
      vrhob(i) = 0.D0
      vsigmaaa(i) = 0.D0
      vsigmaab(i) = 0.D0
      vsigmabb(i) = 0.D0
      v2rhoa2(i) = 0.D0
      v2rhob2(i) = 0.D0
      v2rhoab(i) = 0.D0
      v2rhoasigmaaa(i) = 0.D0
      v2rhoasigmaab(i) = 0.D0
      v2rhoasigmabb(i) = 0.D0
      v2rhobsigmaaa(i) = 0.D0
      v2rhobsigmaab(i) = 0.D0
      v2rhobsigmabb(i) = 0.D0
      v2sigmaaa2(i) = 0.D0
      v2sigmaab2(i) = 0.D0
      v2sigmabb2(i) = 0.D0
      v2sigmaaaab(i) = 0.D0
      v2sigmaaabb(i) = 0.D0
      v2sigmaabbb(i) = 0.D0
      if(rho.gt.tolmin) then
      sigmaaa = dmax1(0.D0,sigmaaa1(i))
      sigmaab = sigmaab1(i)
      sigmabb = dmax1(0.D0,sigmabb1(i))
      sigma = sigmaaa+sigmabb+2.D0*sigmaab
      if(rhoa.gt.tolmin) then
      t1 = 1/rhoa
      t2 = t1**(1.D0/3.D0)
      trsa = 0.1285808356226127D-1*t2
      t3 = t1**(1.D0/15.D0)
      t4 = t3**2
      t6 = dexp(-0.5416509560827845D0*t4)
      t8 = t1**(1.D0/6.D0)
      t10 = dexp(-0.8579886777788441D0*t8)
      t13 = (0.1247511874D1+0.812904345D0*t6-0.859614445D0*t10)**2
      t14 = sigmaaa**2
      t15 = rhoa**2
      t16 = t15**2
      t18 = rhoa**(1.D0/3.D0)
      t21 = t14/t18/t16/rhoa
      t24 = (0.1D1+0.4473229026479948D-1*t21)**2
      t26 = t18**2
      t36 = dexp(-0.4473229026479948D-1*t21)
      t37 = t36**2
      trcaa = t13*t24/(0.1D1+0.2481823842598833D0*sigmaaa/t26/t15/t2)
     &*t37
      if(trsa.lt.tolmax*trcaa) then
      t2 = 1/rhoa
      t3 = t2**(1.D0/3.D0)
      t4 = t3**2
      t5 = t2**(1.D0/6.D0)
      t8 = 0.7395908974809453D0*t5+0.1075172860312113D1*t3
      t9 = t8**2
      t10 = 1/t9
      t13 = dexp(-0.3848347315591266D0*t4*t10)
      t14 = rhoa*t13
      t15 = 0.3141592653589793D1**2
      t18 = 1/0.3141592653589793D1*t2
      t19 = t18**(1.D0/3.D0)
      t21 = t18**(1.D0/15.D0)
      t22 = t21**2
      t24 = dexp(-0.630966299458536D0*t22)
      t26 = t18**(1.D0/6.D0)
      t28 = dexp(-0.1038340679664977D1*t26)
      t31 = (0.1247511874D1+0.812904345D0*t24-0.859614445D0*t28)**2
      t34 = 0.3141592653589793D1**(1.D0/3.D0)
      t37 = sigmaaa**2
      t39 = rhoa**2
      t40 = t39**2
      t41 = t40*rhoa
      t42 = rhoa**(1.D0/3.D0)
      t44 = 1/t42/t41
      t45 = 1/t34/0.3141592653589793D1*t37*t44
      t48 = (0.1D1+0.2058200272046996D0*t45)**2
      t50 = t34**2
      t53 = t42**2
      t55 = 1/t53/t39
      t63 = dexp(0.4116400544093991D0*t45)
      t67 = expei(-0.1858628590577086D0/t15*t19/t31/t48*(0.1D1
     &+0.3634859066227017D0/t50*sigmaaa*t55/t19)*t63)
      t69 = t2**(1.D0/15.D0)
      t70 = t69**2
      t72 = dexp(-0.5416509560827845D0*t70)
      t75 = dexp(-0.8579886777788441D0*t5)
      t77 = 0.1247511874D1+0.812904345D0*t72-0.859614445D0*t75
      t78 = t77**2
      t79 = 1/t78
      t80 = t3*t79
      t81 = t37*t44
      t83 = 0.1D1+0.4473229026479948D-1*t81
      t84 = t83**2
      t85 = 1/t84
      t87 = 1/t3
      t90 = 0.1D1+0.2481823842598833D0*sigmaaa*t55*t87
      t91 = t85*t90
      t93 = dexp(0.8946458052959895D-1*t81)
      t95 = t80*t91*t93
      t96 = dsqrt(t95)
      t99 = 0.6D1+0.4535739597862518D0*t96+0.5143233424904509D-1*t95
      t102 = 0.3D1+0.6803609396793777D0*t96+0.7714850137356763D-1*t95
      t103 = 1/t102
      t104 = t99*t103
      t105 = t80*t85
      t106 = t90*t67
      t110 = 0.1285808356226127D-1*t105*t106*t93+0.1D1
      t112 = 1.D0*t67+t104*t110
      zk(i) = zk(i)+0.1554534543482745D-1*t14*t112
      t116 = t13*t112
      t118 = t87*t10
      t119 = 1/t39
      t123 = 1/t9/t8
      t124 = t4*t123
      t125 = t5**2
      t126 = t125**2
      t127 = t126*t5
      t128 = 1/t127
      t129 = t128*t119
      t131 = 1/t4
      t134 = -0.1232651495801576D0*t129-0.358390953437371D0*t131*t119
      t137 = 0.2565564877060844D0*t118*t119+0.7696694631182532D0*t124
     &*t134
      t138 = rhoa*t137
      t141 = t131*t79
      t142 = t141*t85
      t143 = t90*t93
      t145 = t142*t143*t119
      t146 = 0.4286027854087091D-2*t145
      t148 = 1/t78/t77
      t149 = t3*t148
      t150 = t149*t85
      t151 = t70**2
      t153 = t151**2
      t154 = t153*t151*t69
      t155 = 1/t154
      t161 = 0.5870805542307996D-1*t155*t119*t72-0.1229232435108575D0
     &*t129*t75
      t163 = t150*t143*t161
      t164 = 0.2571616712452254D-1*t163
      t166 = 1/t84/t83
      t167 = t80*t166
      t168 = t40*t39
      t170 = 1/t42/t168
      t171 = t37*t170
      t172 = t143*t171
      t173 = t167*t172
      t174 = 0.613516294566526D-2*t173
      t175 = t39*rhoa
      t177 = 1/t53/t175
      t182 = 1/t53/t40
      t183 = sigmaaa*t182
      t185 = 1/t3/t2
      t188 = -0.6618196913596887D0*sigmaaa*t177*t87
     &+0.8272746141996109D-1*t183*t185
      t189 = t85*t188
      t191 = t80*t189*t93
      t192 = 0.1285808356226127D-1*t191
      t193 = t105*t172
      t194 = 0.613516294566526D-2*t193
      t195 = -t146-t164+t174+t192-t194
      t196 = t195*t67
      t198 = -t195
      t199 = t198*t87
      t200 = t199*t78
      t201 = 1/t90
      t202 = t84*t201
      t204 = dexp(-0.8946458052959895D-1*t81)
      t205 = t202*t204
      t208 = 1/t96
      t214 = -0.2306016315371713D-1*t145-0.1383609789223028D0*t163
     &+0.3300908517586267D-1*t173+0.6918048946115138D-1*t191
     &-0.3300908517586267D-1*t193
      t215 = t208*t214
      t222 = 0.3278192763011299D1*t215-0.1714411141634836D-1*t145
     &-0.1028646684980902D0*t163+0.2454065178266104D-1*t173
     &+0.5143233424904509D-1*t191-0.2454065178266104D-1*t193
      t223 = t222*t103
      t225 = t102**2
      t226 = 1/t225
      t227 = t99*t226
      t234 = 0.4917289144516948D1*t215-0.2571616712452254D-1*t145
     &-0.1542970027471353D0*t163+0.3681097767399156D-1*t173
     &+0.7714850137356763D-1*t191-0.3681097767399156D-1*t193
      t235 = t110*t234
      t238 = t119*t93
      t242 = t161*t93
      t246 = t166*t90
      t247 = t80*t246
      t248 = t67*t37
      t249 = t170*t93
      t250 = t248*t249
      t253 = t188*t67
      t257 = t80*t91
      t261 = t67*t93
      t265 = -0.4286027854087091D-2*t142*t106*t238
     &-0.2571616712452254D-1*t150*t106*t242+0.613516294566526D-2*t247
     &*t250+0.1285808356226127D-1*t105*t253*t93-0.613516294566526D-2
     &*t257*t250+0.1285808356226127D-1*t105*t90*t195*t261-t146-t164
     &+t174+t192-t194
      t267 = 1.D0*t196-0.7777208750882749D2*t200*t205+t223*t110-1.D0
     &*t227*t235+t104*t265
      vrhoa(i) = vrhoa(i)+0.1554534543482745D-1*t116
     &+0.1554534543482745D-1*t138*t116+0.1554534543482745D-1*t14*t267
      t271 = sigmaaa*t44
      t272 = t143*t271
      t273 = t167*t272
      t274 = 0.2300686104624472D-2*t273
      t275 = t79*t85
      t276 = t55*t93
      t277 = t275*t276
      t278 = 0.3191149835494816D-2*t277
      t279 = t105*t272
      t280 = 0.2300686104624472D-2*t279
      t281 = -t274+t278+t280
      t282 = t281*t67
      t284 = -t281
      t285 = t284*t87
      t292 = -0.123784069409485D-1*t273+0.1716937881873428D-1*t277
     &+0.123784069409485D-1*t279
      t293 = t208*t292
      t298 = 0.3278192763011299D1*t293-0.9202744418497889D-2*t273
     &+0.1276459934197926D-1*t277+0.9202744418497889D-2*t279
      t299 = t298*t103
      t305 = 0.4917289144516948D1*t293-0.1380411662774683D-1*t273
     &+0.1914689901296889D-1*t277+0.1380411662774683D-1*t279
      t306 = t110*t305
      t309 = t67*sigmaaa
      t310 = t44*t93
      t311 = t309*t310
      t324 = -0.2300686104624472D-2*t247*t311+0.3191149835494816D-2
     &*t275*t55*t67*t93+0.2300686104624472D-2*t257*t311
     &+0.1285808356226127D-1*t105*t90*t281*t261-t274+t278+t280
      t326 = 1.D0*t282-0.7777208750882749D2*t285*t78*t205+t299*t110
     &-1.D0*t227*t306+t104*t324
      vsigmaaa(i) = vsigmaaa(i)+0.1554534543482745D-1*t14*t326
      t333 = t13*t267
      t336 = 1/t40
      t343 = 1/t175
      t346 = t9**2
      t349 = t134**2
      t354 = 1/t127/t2*t336
      t356 = t128*t343
      t359 = 1/t4/t2
      t360 = t359*t336
      t371 = t137**2
      t377 = t141*t166
      t378 = t40**2
      t380 = 1/t42/t378
      t382 = t143*t380*t37
      t383 = t377*t382
      t384 = 0.4090108630443506D-2*t383
      t385 = t188*t93
      t387 = t150*t385*t161
      t388 = 0.5143233424904509D-1*t387
      t410 = 0.5088031470000263D-1/t154/t2*t336*t72
     &-0.1174161108461599D0*t155*t343*t72+0.4239903246622982D-2/t153
     &/t70/t69/t2*t336*t72-0.1024360362590479D0*t354*t75
     &+0.245846487021715D0*t356*t75-0.1757779186136125D-1*t360*t75
      t412 = t150*t143*t410
      t413 = 0.2571616712452254D-1*t412
      t415 = t142*t385*t119
      t416 = 0.8572055708174181D-2*t415
      t418 = t142*t143*t343
      t419 = 0.8572055708174181D-2*t418
      t421 = t359*t79*t85
      t423 = t421*t143*t336
      t424 = 0.2857351902724727D-2*t423
      t425 = t142*t382
      t426 = 0.4090108630443506D-2*t425
      t427 = t78**2
      t430 = t3/t427*t85
      t431 = t161**2
      t433 = t430*t143*t431
      t434 = 0.7714850137356763D-1*t433
      t435 = t131*t148
      t439 = t435*t85*t143*t119*t161
      t440 = 0.1714411141634836D-1*t439
      t443 = 1/t42/t40/t175
      t445 = t143*t37*t443
      t446 = t105*t445
      t447 = 0.3885603198921331D-1*t446
      t448 = t37**2
      t451 = 1/t53/t378/t40
      t453 = t143*t448*t451
      t454 = t105*t453
      t455 = 0.2927358823544966D-2*t454
      t470 = 0.2426672201652192D1*t183*t87-0.606668050413048D0*sigmaaa
     &/t53/t41*t185+0.1103032818932815D0*sigmaaa/t53/t168/t3/t119
      t473 = t80*t85*t470*t93
      t474 = 0.1285808356226127D-1*t473
      t475 = t385*t171
      t476 = t105*t475
      t477 = 0.1227032589133052D-1*t476
      t478 = t167*t453
      t479 = 0.5854717647089932D-2*t478
      t480 = t167*t445
      t481 = 0.3885603198921331D-1*t480
      t482 = t84**2
      t483 = 1/t482
      t484 = t80*t483
      t485 = t484*t453
      t486 = 0.4391038235317449D-2*t485
      t487 = t167*t475
      t488 = 0.1227032589133052D-1*t487
      t489 = t149*t91
      t490 = t242*t171
      t491 = t489*t490
      t492 = 0.2454065178266104D-1*t491
      t493 = t149*t246
      t494 = t493*t490
      t495 = 0.2454065178266104D-1*t494
      t496 = -t384-t388-t413-t416+t419-t424+t426+t434+t440+t447+t455
     &+t474-t477-t479-t481+t486+t488+t492-t495
      t499 = t195**2
      t502 = t195*t198
      t504 = t78*t84
      t506 = t504*t201*t204
      t514 = t198**2
      t530 = t78*t83
      t531 = t199*t530
      t534 = t201*t37*t170*t204
      t537 = t90**2
      t539 = t84/t537
      t544 = t199*t504
      t586 = -0.2200605678390844D-1*t383-0.2767219578446055D0*t387
     &-0.1383609789223028D0*t412-0.4612032630743426D-1*t415
     &+0.4612032630743426D-1*t418-0.1537344210247809D-1*t423
     &+0.2200605678390844D-1*t425+0.4150829367669083D0*t433
     &+0.9224065261486851D-1*t439+0.2090575394471302D0*t446
     &+0.1575010111426323D-1*t454+0.6918048946115138D-1*t473
     &-0.6601817035172533D-1*t476-0.3150020222852647D-1*t478
     &-0.2090575394471302D0*t480+0.2362515167139485D-1*t485
     &+0.6601817035172533D-1*t487+0.1320363407034507D0*t491
     &-0.1320363407034507D0*t494
      t587 = t208*t586
      t590 = 1/t96/t95
      t591 = t214**2
      t592 = t590*t591
      t594 = 0.1170943529417986D-1*t454+0.5143233424904509D-1*t473
     &-0.4908130356532208D-1*t476-0.2341887058835973D-1*t478
     &-0.1554241279568532D0*t480+0.175641529412698D-1*t485
     &+0.4908130356532208D-1*t487+0.9816260713064415D-1*t491
     &-0.9816260713064415D-1*t494+0.3278192763011299D1*t587
     &-0.2369304401100098D2*t592
      t598 = t222*t226
      t605 = t99/t225/t102
      t606 = t234**2
      t635 = 0.175641529412698D-1*t454+0.7714850137356763D-1*t473
     &-0.7362195534798311D-1*t476-0.3512830588253959D-1*t478
     &-0.2331361919352799D0*t480+0.263462294119047D-1*t485
     &+0.7362195534798311D-1*t487+0.1472439106959662D0*t491
     &-0.1472439106959662D0*t494+0.4917289144516948D1*t587
     &-0.3553956601650148D2*t592
      t640 = t67*t161
      t641 = t171*t93
      t642 = t640*t641
      t656 = t141*t91
      t662 = t67*t448*t451*t93
      t668 = t80*t166*t188
      t677 = t248*t443*t93
      t681 = t80*t483*t90
      t686 = t474+0.2454065178266104D-1*t489*t642
     &+0.8572055708174181D-2*t142*t106*t343*t93-0.8572055708174181D-2
     &*t142*t253*t238+0.1285808356226127D-1*t105*t470*t67*t93
     &-0.8572055708174181D-2*t656*t196*t238-0.5854717647089932D-2*t247
     &*t662+0.2927358823544966D-2*t257*t662+0.1227032589133052D-1*t668
     &*t250+0.1714411141634836D-1*t435*t91*t67*t119*t242
     &+0.3885603198921331D-1*t257*t677+0.4391038235317449D-2*t681*t662
     &-0.3885603198921331D-1*t247*t677
      t687 = t196*t641
      t690 = t80*t189
      t701 = 0.1227032589133052D-1*t247*t687-0.1227032589133052D-1
     &*t690*t250-0.5143233424904509D-1*t489*t196*t242-t384-t388-t413
     &-t416+t419-0.2454065178266104D-1*t493*t642-t424+t426
     &-0.5143233424904509D-1*t150*t253*t242+t434
      t713 = 0.1285808356226127D-1*t105*t90*t499*t261+t440+t447+t455
     &-t477-t479-0.1227032589133052D-1*t257*t687-t481+t486+t488+t492
     &-t495+0.7714850137356763D-1*t430*t106*t431*t93
      t729 = 1/t77
      t741 = 1/t83
      t747 = t67*t380*t37*t93
      t753 = t141*t246
      t756 = -0.2857351902724727D-2*t421*t106*t336*t93
     &+0.2571616712452254D-1*t105*t188*t195*t261+0.1285808356226127D-1
     &*t105*t90*t496*t261+0.4771444294911944D0*t198*t37*t170+0.2D1
     &*t729*t198*t161-0.1D1*t188*t198*t201-0.2571616712452254D-1*t150
     &*t106*t410*t93-0.1D1*t514-0.4771444294911944D0*t741*t198*t171
     &+0.4090108630443506D-2*t656*t747+0.3333333333333333D0*t2*t198
     &-0.2D1*t502-0.4090108630443506D-2*t753*t747
      s1 = 1.D0*t496*t67+1.D0*t499*t67-0.155544175017655D3*t502*t87
     &*t506+0.7777208750882749D2*t496*t87*t78*t205
     &-0.7777208750882749D2*t514*t87*t78*t205-0.2592402916960916D2
     &*t198*t185*t78*t202*t119*t204-0.155544175017655D3*t199*t77*t202
     &*t161*t204+0.3710851832473874D2*t531*t534
      t760 = s1+0.7777208750882749D2*t200*t539*t188*t204
     &-0.3710851832473874D2*t544*t534+(-0.1636043452177403D-1*t383
     &-0.2057293369961804D0*t387-0.1028646684980902D0*t412
     &-0.3428822283269673D-1*t415+0.3428822283269673D-1*t418
     &-0.1142940761089891D-1*t423+0.1636043452177403D-1*t425
     &+0.3085940054942705D0*t433+0.6857644566539345D-1*t439
     &+0.1554241279568532D0*t446+t594)*t103*t110-2.D0*t598*t235+2.D0
     &*t223*t265+2.D0*t605*t110*t606-2.D0*t227*t265*t234-1.D0*t227
     &*t110*(-0.2454065178266104D-1*t383-0.3085940054942705D0*t387
     &-0.1542970027471353D0*t412-0.5143233424904509D-1*t415
     &+0.5143233424904509D-1*t418-0.1714411141634836D-1*t423
     &+0.2454065178266104D-1*t425+0.4628910082414058D0*t433
     &+0.1028646684980902D0*t439+0.2331361919352799D0*t446+t635)+t104*
     &(t686+t701+t713+t756)
      v2rhoa2(i) = v2rhoa2(i)+0.310906908696549D-1*t137*t13*t112
     &+0.310906908696549D-1*t333+0.1554534543482745D-1*rhoa*
     &(0.8551882923536146D-1*t185*t10*t336-0.1026225950824338D1*t87
     &*t123*t119*t134-0.5131129754121688D0*t118*t343
     &-0.2309008389354759D1*t4/t346*t349+0.7696694631182532D0*t124*(
     &-0.1027209579834646D0*t354+0.2465302991603151D0*t356
     &-0.2389273022915807D0*t360+0.7167819068747421D0*t131*t343))*t116
     &+0.1554534543482745D-1*rhoa*t371*t116+0.310906908696549D-1*t138
     &*t333+0.1554534543482745D-1*t14*t760
      t764 = t13*t326
      t769 = t143*t443*sigmaaa
      t770 = t377*t769
      t771 = 0.7668953682081574D-3*t770
      t772 = t177*t79
      t774 = t772*t85*t93
      t775 = 0.1063716611831605D-2*t774
      t776 = t142*t769
      t777 = 0.7668953682081574D-3*t776
      t778 = t242*t271
      t779 = t493*t778
      t780 = 0.4601372209248945D-2*t779
      t781 = t148*t85
      t783 = t781*t276*t161
      t784 = 0.6382299670989632D-2*t783
      t785 = t489*t778
      t786 = 0.4601372209248945D-2*t785
      t787 = t37*sigmaaa
      t790 = 1/t53/t378/t175
      t792 = t143*t787*t790
      t793 = t484*t792
      t794 = 0.1646639338244043D-2*t793
      t795 = t79*t166
      t797 = 1/t378/rhoa
      t799 = t797*t93*t37
      t800 = t795*t799
      t801 = 0.1522639367678093D-2*t800
      t802 = t167*t792
      t803 = 0.2195519117658725D-2*t802
      t805 = t143*sigmaaa*t170
      t806 = t167*t805
      t807 = 0.1227032589133052D-1*t806
      t808 = t385*t271
      t809 = t167*t808
      t810 = 0.2300686104624472D-2*t809
      t815 = -0.6618196913596887D0*t177*t87+0.8272746141996109D-1*t182
     &*t185
      t818 = t80*t85*t815*t93
      t819 = 0.1285808356226127D-1*t818
      t820 = t105*t808
      t821 = 0.2300686104624472D-2*t820
      t822 = t275*t799
      t823 = 0.1522639367678093D-2*t822
      t824 = t105*t792
      t825 = 0.1097759558829362D-2*t824
      t826 = t105*t805
      t827 = 0.1227032589133052D-1*t826
      t828 = t771-t775-t777+t780-t784-t786-t794+t801+t803+t807-t810
     &+t819+t821-t823-t825-t827
      t831 = t195*t281
      t834 = t195*t284
      t838 = t281*t198
      t847 = t198*t284
      t853 = t201*sigmaaa*t44*t204
      t859 = t539*t55*t204
      t865 = t590*t214*t292
      t883 = 0.4126135646982833D-2*t770-0.5723126272911426D-2*t774
     &-0.4126135646982833D-2*t776+0.24756813881897D-1*t779
     &-0.3433875763746856D-1*t783-0.24756813881897D-1*t785
     &-0.8859431876773069D-2*t793+0.8192273461183165D-2*t800
     &+0.1181257583569743D-1*t802+0.6601817035172533D-1*t806
     &-0.123784069409485D-1*t809+0.6918048946115138D-1*t818
     &+0.123784069409485D-1*t820-0.8192273461183165D-2*t822
     &-0.5906287917848713D-2*t824-0.6601817035172533D-1*t826
      t884 = t208*t883
      t902 = -0.2369304401100098D2*t865+0.3278192763011299D1*t884
     &+0.306758147283263D-2*t770-0.4254866447326421D-2*t774
     &-0.306758147283263D-2*t776+0.1840548883699578D-1*t779
     &-0.2552919868395853D-1*t783-0.1840548883699578D-1*t785
     &-0.6586557352976174D-2*t793+0.6090557470712371D-2*t800
     &+0.8782076470634898D-2*t802+0.4908130356532208D-1*t806
     &-0.9202744418497889D-2*t809+0.5143233424904509D-1*t818
     &+0.9202744418497889D-2*t820-0.6090557470712371D-2*t822
     &-0.4391038235317449D-2*t824-0.4908130356532208D-1*t826
      t908 = t298*t226
      t935 = -0.3553956601650148D2*t865+0.4917289144516948D1*t884
     &+0.4601372209248945D-2*t770-0.6382299670989632D-2*t774
     &-0.4601372209248945D-2*t776+0.2760823325549367D-1*t779
     &-0.3829379802593779D-1*t783-0.2760823325549367D-1*t785
     &-0.9879836029464261D-2*t793+0.9135836206068557D-2*t800
     &+0.1317311470595235D-1*t802+0.7362195534798311D-1*t806
     &-0.1380411662774683D-1*t809+0.7714850137356763D-1*t818
     &+0.1380411662774683D-1*t820-0.9135836206068557D-2*t822
     &-0.6586557352976174D-2*t824-0.7362195534798311D-1*t826
      t943 = t271*t93
      t944 = t196*t943
      t947 = t275*t55
      t951 = -0.2300686104624472D-2*t247*t944+t780+t771
     &+0.3191149835494816D-2*t947*t196*t93-t777-t775-t794+t801+t803
     &+t821-t823-t825
      t952 = t640*t943
      t963 = t67*t787*t790*t93
      t966 = t309*t249
      t974 = -t786-t784-t827-0.4601372209248945D-2*t489*t952
     &+0.4601372209248945D-2*t493*t952-0.1063716611831605D-2*t772*t85
     &*t67*t93+t819-0.1097759558829362D-2*t257*t963-t810+t807
     &-0.1227032589133052D-1*t257*t966+0.1285808356226127D-1*t257*t831
     &*t261+0.2195519117658725D-2*t247*t963
      t978 = t67*t443*sigmaaa*t93
      t989 = t282*t641
      t1004 = t741*t284
      t1010 = 0.7668953682081574D-3*t753*t978-0.6382299670989632D-2
     &*t781*t55*t640*t93+0.1285808356226127D-1*t105*t815*t67*t93
     &+0.613516294566526D-2*t247*t989-0.1D1*t834-0.613516294566526D-2
     &*t257*t989+0.2D1*t729*t284*t161-0.1D1*t188*t284*t201
     &+0.4771444294911944D0*t284*t37*t170-0.4771444294911944D0*t1004
     &*t171+0.3333333333333333D0*t2*t284-0.1D1*t838
      t1020 = t248*t93
      t1045 = -0.1D1*t847-0.2571616712452254D-1*t489*t282*t242
     &+0.1227032589133052D-1*t247*t966-0.1646639338244043D-2*t681*t963
     &-0.1522639367678093D-2*t275*t797*t1020+0.1285808356226127D-1
     &*t105*t188*t281*t261-0.7668953682081574D-3*t656*t978
     &-0.4286027854087091D-2*t656*t282*t238+0.2300686104624472D-2*t690
     &*t311+0.1285808356226127D-1*t105*t90*t828*t261
     &+0.1522639367678093D-2*t795*t797*t1020+0.2300686104624472D-2
     &*t257*t944-0.2300686104624472D-2*t668*t311
      t1049 = 1.D0*t828*t67+1.D0*t831*t67-0.7777208750882749D2*t834
     &*t87*t506-0.7777208750882749D2*t838*t87*t506
     &+0.7777208750882749D2*t828*t87*t78*t205-0.7777208750882749D2
     &*t847*t87*t506-0.1391569437177703D2*t531*t853
     &+0.1930166210680909D2*t198*t131*t78*t859+0.1391569437177703D2
     &*t544*t853+t902*t103*t110-1.D0*t598*t306+t223*t324-1.D0*t908
     &*t235+2.D0*t605*t235*t305-1.D0*t227*t324*t234-1.D0*t227*t110
     &*t935+t299*t265-1.D0*t227*t265*t305+t104*(t951+t974+t1010+t1045)
      v2rhoasigmaaa(i) = v2rhoasigmaaa(i)+0.1554534543482745D-1*t764
     &+0.1554534543482745D-1*t138*t764+0.1554534543482745D-1*t14*t1049
      t1055 = 1/t53/t378/t39
      t1057 = t143*t37*t1055
      t1058 = t484*t1057
      t1059 = 0.6174897518415163D-3*t1058
      t1060 = 1/t378
      t1062 = t1060*t93*sigmaaa
      t1063 = t795*t1062
      t1064 = 0.114197952575857D-2*t1063
      t1065 = t167*t1057
      t1066 = 0.8233196691220217D-3*t1065
      t1067 = t143*t44
      t1068 = t167*t1067
      t1069 = 0.2300686104624472D-2*t1068
      t1070 = t275*t1062
      t1071 = 0.114197952575857D-2*t1070
      t1072 = t105*t1057
      t1073 = 0.4116598345610109D-3*t1072
      t1074 = t105*t1067
      t1075 = 0.2300686104624472D-2*t1074
      t1076 = t1059-t1064-t1066-t1069+t1071+t1073+t1075
      t1079 = t281**2
      t1082 = t281*t284
      t1091 = t284**2
      t1106 = t292**2
      t1107 = t590*t1106
      t1117 = t208*(0.3322286953789901D-2*t1058-0.6144205095887373D-2
     &*t1063-0.4429715938386535D-2*t1065-0.123784069409485D-1*t1068
     &+0.6144205095887373D-2*t1070+0.2214857969193267D-2*t1072
     &+0.123784069409485D-1*t1074)
      t1133 = t305**2
      t1154 = t309*t93
      t1168 = t248*t1055*t93
      t1171 = t282*t943
      t1180 = t106*t310
      t1185 = 0.114197952575857D-2*t275*t1060*t1154
     &+0.6382299670989632D-2*t947*t282*t93-0.2481823842598833D0*t55
     &*t284*t87*t201-0.114197952575857D-2*t795*t1060*t1154
     &-0.8233196691220217D-3*t247*t1168-0.4601372209248945D-2*t247
     &*t1171+0.1789291610591979D0*t1004*t271+0.6174897518415163D-3
     &*t681*t1168+0.4116598345610109D-3*t257*t1168
     &-0.2300686104624472D-2*t167*t1180-0.2D1*t1082-0.1D1*t1091
      t1201 = 0.1285808356226127D-1*t105*t90*t1079*t261+t1059-t1064
     &-t1066-t1069+t1071+t1073+t1075+0.4601372209248945D-2*t257*t1171
     &-0.1789291610591979D0*t284*sigmaaa*t44+0.1285808356226127D-1
     &*t105*t90*t1076*t261+0.2300686104624472D-2*t105*t1180
      s1 = 1.D0*t1076*t67+1.D0*t1079*t67-0.155544175017655D3*t1082*t87
     &*t506+0.7777208750882749D2*t1076*t87*t78*t205
     &-0.7777208750882749D2*t1091*t87*t78*t205-0.1391569437177703D2
     &*t285*t530*t853+0.1930166210680909D2*t284*t131*t78*t859
      t1204 = s1+0.1391569437177703D2*t285*t504*t853+(
     &-0.2369304401100098D2*t1107+0.3278192763011299D1*t1117
     &+0.2469959007366065D-2*t1058-0.4567918103034278D-2*t1063
     &-0.3293278676488087D-2*t1065-0.9202744418497889D-2*t1068
     &+0.4567918103034278D-2*t1070+0.1646639338244043D-2*t1072
     &+0.9202744418497889D-2*t1074)*t103*t110-2.D0*t908*t306+2.D0*t299
     &*t324+2.D0*t605*t110*t1133-2.D0*t227*t324*t305-1.D0*t227*t110*(
     &-0.3553956601650148D2*t1107+0.4917289144516948D1*t1117
     &+0.3704938511049098D-2*t1058-0.6851877154551418D-2*t1063
     &-0.493991801473213D-2*t1065-0.1380411662774683D-1*t1068
     &+0.6851877154551418D-2*t1070+0.2469959007366065D-2*t1072
     &+0.1380411662774683D-1*t1074)+t104*(t1185+t1201)
      v2sigmaaa2(i) = v2sigmaaa2(i)+0.1554534543482745D-1*t14*t1204
      endif
      endif ! rhoa
      if(rhob.gt.tolmin) then
      t1 = 1/rhob
      t2 = t1**(1.D0/3.D0)
      trsb = 0.1285808356226127D-1*t2
      t3 = t1**(1.D0/15.D0)
      t4 = t3**2
      t6 = dexp(-0.5416509560827845D0*t4)
      t8 = t1**(1.D0/6.D0)
      t10 = dexp(-0.8579886777788441D0*t8)
      t13 = (0.1247511874D1+0.812904345D0*t6-0.859614445D0*t10)**2
      t14 = sigmabb**2
      t15 = rhob**2
      t16 = t15**2
      t18 = rhob**(1.D0/3.D0)
      t21 = t14/t18/t16/rhob
      t24 = (0.1D1+0.4473229026479948D-1*t21)**2
      t26 = t18**2
      t36 = dexp(-0.4473229026479948D-1*t21)
      t37 = t36**2
      trcbb = t13*t24/(0.1D1+0.2481823842598833D0*sigmabb/t26/t15/t2)
     &*t37
      if(trsb.lt.tolmax*trcbb) then
      t2 = 1/rhob
      t3 = t2**(1.D0/3.D0)
      t4 = t3**2
      t5 = t2**(1.D0/6.D0)
      t8 = 0.7395908974809453D0*t5+0.1075172860312113D1*t3
      t9 = t8**2
      t10 = 1/t9
      t13 = dexp(-0.3848347315591266D0*t4*t10)
      t14 = rhob*t13
      t15 = 0.3141592653589793D1**2
      t18 = 1/0.3141592653589793D1*t2
      t19 = t18**(1.D0/3.D0)
      t21 = t18**(1.D0/15.D0)
      t22 = t21**2
      t24 = dexp(-0.630966299458536D0*t22)
      t26 = t18**(1.D0/6.D0)
      t28 = dexp(-0.1038340679664977D1*t26)
      t31 = (0.1247511874D1+0.812904345D0*t24-0.859614445D0*t28)**2
      t34 = 0.3141592653589793D1**(1.D0/3.D0)
      t37 = sigmabb**2
      t39 = rhob**2
      t40 = t39**2
      t41 = t40*rhob
      t42 = rhob**(1.D0/3.D0)
      t44 = 1/t42/t41
      t45 = 1/t34/0.3141592653589793D1*t37*t44
      t48 = (0.1D1+0.2058200272046996D0*t45)**2
      t50 = t34**2
      t53 = t42**2
      t55 = 1/t53/t39
      t63 = dexp(0.4116400544093991D0*t45)
      t67 = expei(-0.1858628590577086D0/t15*t19/t31/t48*(0.1D1
     &+0.3634859066227017D0/t50*sigmabb*t55/t19)*t63)
      t69 = t2**(1.D0/15.D0)
      t70 = t69**2
      t72 = dexp(-0.5416509560827845D0*t70)
      t75 = dexp(-0.8579886777788441D0*t5)
      t77 = 0.1247511874D1+0.812904345D0*t72-0.859614445D0*t75
      t78 = t77**2
      t79 = 1/t78
      t80 = t3*t79
      t81 = t37*t44
      t83 = 0.1D1+0.4473229026479948D-1*t81
      t84 = t83**2
      t85 = 1/t84
      t87 = 1/t3
      t90 = 0.1D1+0.2481823842598833D0*sigmabb*t55*t87
      t91 = t85*t90
      t93 = dexp(0.8946458052959895D-1*t81)
      t95 = t80*t91*t93
      t96 = dsqrt(t95)
      t99 = 0.6D1+0.4535739597862518D0*t96+0.5143233424904509D-1*t95
      t102 = 0.3D1+0.6803609396793777D0*t96+0.7714850137356763D-1*t95
      t103 = 1/t102
      t104 = t99*t103
      t105 = t80*t85
      t106 = t90*t67
      t110 = 0.1285808356226127D-1*t105*t106*t93+0.1D1
      t112 = 1.D0*t67+t104*t110
      zk(i) = zk(i)+0.1554534543482745D-1*t14*t112
      t116 = t13*t112
      t118 = t87*t10
      t119 = 1/t39
      t123 = 1/t9/t8
      t124 = t4*t123
      t125 = t5**2
      t126 = t125**2
      t127 = t126*t5
      t128 = 1/t127
      t129 = t128*t119
      t131 = 1/t4
      t134 = -0.1232651495801576D0*t129-0.358390953437371D0*t131*t119
      t137 = 0.2565564877060844D0*t118*t119+0.7696694631182532D0*t124
     &*t134
      t138 = rhob*t137
      t141 = t131*t79
      t142 = t141*t85
      t143 = t90*t93
      t145 = t142*t143*t119
      t146 = 0.4286027854087091D-2*t145
      t148 = 1/t78/t77
      t149 = t3*t148
      t150 = t149*t85
      t151 = t70**2
      t153 = t151**2
      t154 = t153*t151*t69
      t155 = 1/t154
      t161 = 0.5870805542307996D-1*t155*t119*t72-0.1229232435108575D0
     &*t129*t75
      t163 = t150*t143*t161
      t164 = 0.2571616712452254D-1*t163
      t166 = 1/t84/t83
      t167 = t80*t166
      t168 = t40*t39
      t170 = 1/t42/t168
      t171 = t37*t170
      t172 = t143*t171
      t173 = t167*t172
      t174 = 0.613516294566526D-2*t173
      t175 = t39*rhob
      t177 = 1/t53/t175
      t182 = 1/t53/t40
      t183 = sigmabb*t182
      t185 = 1/t3/t2
      t188 = -0.6618196913596887D0*sigmabb*t177*t87
     &+0.8272746141996109D-1*t183*t185
      t189 = t85*t188
      t191 = t80*t189*t93
      t192 = 0.1285808356226127D-1*t191
      t193 = t105*t172
      t194 = 0.613516294566526D-2*t193
      t195 = -t146-t164+t174+t192-t194
      t196 = t195*t67
      t198 = -t195
      t199 = t198*t87
      t200 = t199*t78
      t201 = 1/t90
      t202 = t84*t201
      t204 = dexp(-0.8946458052959895D-1*t81)
      t205 = t202*t204
      t208 = 1/t96
      t214 = -0.2306016315371713D-1*t145-0.1383609789223028D0*t163
     &+0.3300908517586267D-1*t173+0.6918048946115138D-1*t191
     &-0.3300908517586267D-1*t193
      t215 = t208*t214
      t222 = 0.3278192763011299D1*t215-0.1714411141634836D-1*t145
     &-0.1028646684980902D0*t163+0.2454065178266104D-1*t173
     &+0.5143233424904509D-1*t191-0.2454065178266104D-1*t193
      t223 = t222*t103
      t225 = t102**2
      t226 = 1/t225
      t227 = t99*t226
      t234 = 0.4917289144516948D1*t215-0.2571616712452254D-1*t145
     &-0.1542970027471353D0*t163+0.3681097767399156D-1*t173
     &+0.7714850137356763D-1*t191-0.3681097767399156D-1*t193
      t235 = t110*t234
      t238 = t119*t93
      t242 = t161*t93
      t246 = t166*t90
      t247 = t80*t246
      t248 = t67*t37
      t249 = t170*t93
      t250 = t248*t249
      t253 = t188*t67
      t257 = t80*t91
      t261 = t67*t93
      t265 = -0.4286027854087091D-2*t142*t106*t238
     &-0.2571616712452254D-1*t150*t106*t242+0.613516294566526D-2*t247
     &*t250+0.1285808356226127D-1*t105*t253*t93-0.613516294566526D-2
     &*t257*t250+0.1285808356226127D-1*t105*t90*t195*t261-t146-t164
     &+t174+t192-t194
      t267 = 1.D0*t196-0.7777208750882749D2*t200*t205+t223*t110-1.D0
     &*t227*t235+t104*t265
      vrhob(i) = vrhob(i)+0.1554534543482745D-1*t116
     &+0.1554534543482745D-1*t138*t116+0.1554534543482745D-1*t14*t267
      t271 = sigmabb*t44
      t272 = t143*t271
      t273 = t167*t272
      t274 = 0.2300686104624472D-2*t273
      t275 = t79*t85
      t276 = t55*t93
      t277 = t275*t276
      t278 = 0.3191149835494816D-2*t277
      t279 = t105*t272
      t280 = 0.2300686104624472D-2*t279
      t281 = -t274+t278+t280
      t282 = t281*t67
      t284 = -t281
      t292 = -0.123784069409485D-1*t273+0.1716937881873428D-1*t277
     &+0.123784069409485D-1*t279
      t293 = t208*t292
      t298 = 0.3278192763011299D1*t293-0.9202744418497889D-2*t273
     &+0.1276459934197926D-1*t277+0.9202744418497889D-2*t279
      t299 = t298*t103
      t305 = 0.4917289144516948D1*t293-0.1380411662774683D-1*t273
     &+0.1914689901296889D-1*t277+0.1380411662774683D-1*t279
      t306 = t110*t305
      t309 = t67*sigmabb
      t311 = t309*t44*t93
      t324 = -0.2300686104624472D-2*t247*t311+0.3191149835494816D-2
     &*t275*t55*t67*t93+0.2300686104624472D-2*t257*t311
     &+0.1285808356226127D-1*t105*t90*t281*t261-t274+t278+t280
      t326 = 1.D0*t282-0.7777208750882749D2*t284*t87*t78*t205+t299
     &*t110-1.D0*t227*t306+t104*t324
      vsigmabb(i) = vsigmabb(i)+0.1554534543482745D-1*t14*t326
      t333 = t13*t267
      t336 = 1/t40
      t343 = 1/t175
      t346 = t9**2
      t349 = t134**2
      t354 = 1/t127/t2*t336
      t356 = t128*t343
      t359 = 1/t4/t2
      t360 = t359*t336
      t371 = t137**2
      t377 = t78**2
      t380 = t3/t377*t85
      t381 = t161**2
      t383 = t380*t143*t381
      t384 = 0.7714850137356763D-1*t383
      t385 = t40**2
      t387 = 1/t42/t385
      t389 = t143*t387*t37
      t390 = t142*t389
      t391 = 0.4090108630443506D-2*t390
      t392 = t141*t166
      t393 = t392*t389
      t394 = 0.4090108630443506D-2*t393
      t396 = t359*t79*t85
      t398 = t396*t143*t336
      t399 = 0.2857351902724727D-2*t398
      t401 = t142*t143*t343
      t402 = 0.8572055708174181D-2*t401
      t403 = t188*t93
      t405 = t142*t403*t119
      t406 = 0.8572055708174181D-2*t405
      t407 = t131*t148
      t411 = t407*t85*t143*t119*t161
      t412 = 0.1714411141634836D-1*t411
      t413 = t149*t246
      t414 = t242*t171
      t415 = t413*t414
      t416 = 0.2454065178266104D-1*t415
      t438 = 0.5088031470000263D-1/t154/t2*t336*t72
     &-0.1174161108461599D0*t155*t343*t72+0.4239903246622982D-2/t153
     &/t70/t69/t2*t336*t72-0.1024360362590479D0*t354*t75
     &+0.245846487021715D0*t356*t75-0.1757779186136125D-1*t360*t75
      t440 = t150*t143*t438
      t441 = 0.2571616712452254D-1*t440
      t443 = t150*t403*t161
      t444 = 0.5143233424904509D-1*t443
      t445 = t37**2
      t448 = 1/t53/t385/t40
      t450 = t143*t445*t448
      t451 = t105*t450
      t452 = 0.2927358823544966D-2*t451
      t455 = 1/t42/t40/t175
      t457 = t143*t37*t455
      t458 = t105*t457
      t459 = 0.3885603198921331D-1*t458
      t460 = t403*t171
      t461 = t105*t460
      t462 = 0.1227032589133052D-1*t461
      t477 = 0.2426672201652192D1*t183*t87-0.606668050413048D0*sigmabb
     &/t53/t41*t185+0.1103032818932815D0*sigmabb/t53/t168/t3/t119
      t480 = t80*t85*t477*t93
      t481 = 0.1285808356226127D-1*t480
      t482 = t84**2
      t483 = 1/t482
      t484 = t80*t483
      t485 = t484*t450
      t486 = 0.4391038235317449D-2*t485
      t487 = t167*t450
      t488 = 0.5854717647089932D-2*t487
      t489 = t167*t457
      t490 = 0.3885603198921331D-1*t489
      t491 = t167*t460
      t492 = 0.1227032589133052D-1*t491
      t493 = t149*t91
      t494 = t493*t414
      t495 = 0.2454065178266104D-1*t494
      t496 = t384+t391-t394-t399+t402-t406+t412-t416-t441-t444+t452
     &+t459-t462+t481+t486-t488-t490+t492+t495
      t499 = t195**2
      t502 = t195*t198
      t504 = t78*t84
      t506 = t504*t201*t204
      t514 = t198**2
      t531 = t199*t78*t83
      t534 = t201*t37*t170*t204
      t537 = t90**2
      t539 = t84/t537
      t544 = t199*t504
      t568 = 1/t96/t95
      t569 = t214**2
      t570 = t568*t569
      t591 = 0.4150829367669083D0*t383+0.2200605678390844D-1*t390
     &-0.2200605678390844D-1*t393-0.1537344210247809D-1*t398
     &+0.4612032630743426D-1*t401-0.4612032630743426D-1*t405
     &+0.9224065261486851D-1*t411-0.1320363407034507D0*t415
     &-0.1383609789223028D0*t440-0.2767219578446055D0*t443
     &+0.1575010111426323D-1*t451+0.2090575394471302D0*t458
     &-0.6601817035172533D-1*t461+0.6918048946115138D-1*t480
     &+0.2362515167139485D-1*t485-0.3150020222852647D-1*t487
     &-0.2090575394471302D0*t489+0.6601817035172533D-1*t491
     &+0.1320363407034507D0*t494
      t592 = t208*t591
      t594 = 0.1170943529417986D-1*t451+0.1554241279568532D0*t458
     &-0.4908130356532208D-1*t461+0.5143233424904509D-1*t480
     &+0.175641529412698D-1*t485-0.2341887058835973D-1*t487
     &-0.1554241279568532D0*t489+0.4908130356532208D-1*t491
     &+0.9816260713064415D-1*t494-0.2369304401100098D2*t570
     &+0.3278192763011299D1*t592
      t598 = t222*t226
      t605 = t99/t225/t102
      t606 = t234**2
      t635 = 0.175641529412698D-1*t451+0.2331361919352799D0*t458
     &-0.7362195534798311D-1*t461+0.7714850137356763D-1*t480
     &+0.263462294119047D-1*t485-0.3512830588253959D-1*t487
     &-0.2331361919352799D0*t489+0.7362195534798311D-1*t491
     &+0.1472439106959662D0*t494-0.3553956601650148D2*t570
     &+0.4917289144516948D1*t592
      t648 = t141*t91
      t651 = t67*t387*t37*t93
      t661 = t141*t246
      t673 = t248*t455*t93
      t677 = t80*t166*t188
      t683 = t80*t189
      t686 = 0.8572055708174181D-2*t142*t106*t343*t93
     &+0.7714850137356763D-1*t380*t106*t381*t93+0.4090108630443506D-2
     &*t648*t651-0.8572055708174181D-2*t142*t253*t238
     &-0.5143233424904509D-1*t150*t253*t242-0.2D1*t502
     &-0.4090108630443506D-2*t661*t651-0.2857351902724727D-2*t396*t106
     &*t336*t93-0.2571616712452254D-1*t150*t106*t438*t93
     &+0.3885603198921331D-1*t257*t673+0.1227032589133052D-1*t677*t250
     &-0.8572055708174181D-2*t648*t196*t238-0.1227032589133052D-1*t683
     &*t250
      t688 = t80*t483*t90
      t691 = t67*t445*t448*t93
      t694 = t171*t93
      t695 = t196*t694
      t709 = t67*t161
      t710 = t709*t694
      t727 = 0.4391038235317449D-2*t688*t691-0.1227032589133052D-1
     &*t257*t695+0.2571616712452254D-1*t105*t188*t195*t261
     &+0.2927358823544966D-2*t257*t691+0.1714411141634836D-1*t407*t91
     &*t67*t119*t242-0.2454065178266104D-1*t413*t710+t391
     &+0.2454065178266104D-1*t493*t710-0.3885603198921331D-1*t247*t673
     &+0.3333333333333333D0*t2*t198+0.1285808356226127D-1*t105*t90
     &*t499*t261+0.1227032589133052D-1*t247*t695-0.5854717647089932D-2
     &*t247*t691
      t736 = 0.1285808356226127D-1*t105*t477*t67*t93-t394-t399+t402
     &-0.5143233424904509D-1*t493*t196*t242+t384-t441-t444+t452-t416
     &+t412-t406-t488
      t745 = 1/t77
      t752 = 1/t83
      t756 = -t490+t492+t459-t462+t481+t486+t495+0.1285808356226127D-1
     &*t105*t90*t496*t261-0.1D1*t514-0.1D1*t188*t198*t201+0.2D1*t745
     &*t198*t161+0.4771444294911944D0*t198*t37*t170
     &-0.4771444294911944D0*t752*t198*t171
      s1 = 1.D0*t496*t67+1.D0*t499*t67-0.155544175017655D3*t502*t87
     &*t506+0.7777208750882749D2*t496*t87*t78*t205
     &-0.7777208750882749D2*t514*t87*t78*t205-0.2592402916960916D2
     &*t198*t185*t78*t202*t119*t204-0.155544175017655D3*t199*t77*t202
     &*t161*t204+0.3710851832473874D2*t531*t534
      t760 = s1+0.7777208750882749D2*t200*t539*t188*t204
     &-0.3710851832473874D2*t544*t534+(0.3085940054942705D0*t383
     &+0.1636043452177403D-1*t390-0.1636043452177403D-1*t393
     &-0.1142940761089891D-1*t398+0.3428822283269673D-1*t401
     &-0.3428822283269673D-1*t405+0.6857644566539345D-1*t411
     &-0.9816260713064415D-1*t415-0.1028646684980902D0*t440
     &-0.2057293369961804D0*t443+t594)*t103*t110-2.D0*t598*t235+2.D0
     &*t223*t265+2.D0*t605*t110*t606-2.D0*t227*t265*t234-1.D0*t227
     &*t110*(0.4628910082414058D0*t383+0.2454065178266104D-1*t390
     &-0.2454065178266104D-1*t393-0.1714411141634836D-1*t398
     &+0.5143233424904509D-1*t401-0.5143233424904509D-1*t405
     &+0.1028646684980902D0*t411-0.1472439106959662D0*t415
     &-0.1542970027471353D0*t440-0.3085940054942705D0*t443+t635)+t104*
     &(t686+t727+t736+t756)
      v2rhob2(i) = v2rhob2(i)+0.310906908696549D-1*t137*t13*t112
     &+0.310906908696549D-1*t333+0.1554534543482745D-1*rhob*
     &(0.8551882923536146D-1*t185*t10*t336-0.1026225950824338D1*t87
     &*t123*t119*t134-0.5131129754121688D0*t118*t343
     &-0.2309008389354759D1*t4/t346*t349+0.7696694631182532D0*t124*(
     &-0.1027209579834646D0*t354+0.2465302991603151D0*t356
     &-0.2389273022915807D0*t360+0.7167819068747421D0*t131*t343))*t116
     &+0.1554534543482745D-1*rhob*t371*t116+0.310906908696549D-1*t138
     &*t333+0.1554534543482745D-1*t14*t760
      t764 = t13*t326
      t769 = t143*t455*sigmabb
      t770 = t392*t769
      t771 = 0.7668953682081574D-3*t770
      t772 = t177*t79
      t774 = t772*t85*t93
      t775 = 0.1063716611831605D-2*t774
      t776 = t142*t769
      t777 = 0.7668953682081574D-3*t776
      t778 = t242*t271
      t779 = t413*t778
      t780 = 0.4601372209248945D-2*t779
      t781 = t148*t85
      t783 = t781*t276*t161
      t784 = 0.6382299670989632D-2*t783
      t785 = t493*t778
      t786 = 0.4601372209248945D-2*t785
      t787 = t37*sigmabb
      t790 = 1/t53/t385/t175
      t792 = t143*t787*t790
      t793 = t484*t792
      t794 = 0.1646639338244043D-2*t793
      t795 = t79*t166
      t797 = 1/t385/rhob
      t799 = t797*t93*t37
      t800 = t795*t799
      t801 = 0.1522639367678093D-2*t800
      t802 = t167*t792
      t803 = 0.2195519117658725D-2*t802
      t805 = t143*sigmabb*t170
      t806 = t167*t805
      t807 = 0.1227032589133052D-1*t806
      t808 = t403*t271
      t809 = t167*t808
      t810 = 0.2300686104624472D-2*t809
      t815 = -0.6618196913596887D0*t177*t87+0.8272746141996109D-1*t182
     &*t185
      t818 = t80*t85*t815*t93
      t819 = 0.1285808356226127D-1*t818
      t820 = t105*t808
      t821 = 0.2300686104624472D-2*t820
      t822 = t275*t799
      t823 = 0.1522639367678093D-2*t822
      t824 = t105*t792
      t825 = 0.1097759558829362D-2*t824
      t826 = t105*t805
      t827 = 0.1227032589133052D-1*t826
      t828 = t771-t775-t777+t780-t784-t786-t794+t801+t803+t807-t810
     &+t819+t821-t823-t825-t827
      t831 = t195*t281
      t834 = t195*t284
      t838 = t281*t198
      t847 = t198*t284
      t853 = t201*sigmabb*t44*t204
      t865 = t568*t214*t292
      t883 = 0.4126135646982833D-2*t770-0.5723126272911426D-2*t774
     &-0.4126135646982833D-2*t776+0.24756813881897D-1*t779
     &-0.3433875763746856D-1*t783-0.24756813881897D-1*t785
     &-0.8859431876773069D-2*t793+0.8192273461183165D-2*t800
     &+0.1181257583569743D-1*t802+0.6601817035172533D-1*t806
     &-0.123784069409485D-1*t809+0.6918048946115138D-1*t818
     &+0.123784069409485D-1*t820-0.8192273461183165D-2*t822
     &-0.5906287917848713D-2*t824-0.6601817035172533D-1*t826
      t884 = t208*t883
      t902 = -0.2369304401100098D2*t865+0.3278192763011299D1*t884
     &+0.306758147283263D-2*t770-0.4254866447326421D-2*t774
     &-0.306758147283263D-2*t776+0.1840548883699578D-1*t779
     &-0.2552919868395853D-1*t783-0.1840548883699578D-1*t785
     &-0.6586557352976174D-2*t793+0.6090557470712371D-2*t800
     &+0.8782076470634898D-2*t802+0.4908130356532208D-1*t806
     &-0.9202744418497889D-2*t809+0.5143233424904509D-1*t818
     &+0.9202744418497889D-2*t820-0.6090557470712371D-2*t822
     &-0.4391038235317449D-2*t824-0.4908130356532208D-1*t826
      t935 = -0.3553956601650148D2*t865+0.4917289144516948D1*t884
     &+0.4601372209248945D-2*t770-0.6382299670989632D-2*t774
     &-0.4601372209248945D-2*t776+0.2760823325549367D-1*t779
     &-0.3829379802593779D-1*t783-0.2760823325549367D-1*t785
     &-0.9879836029464261D-2*t793+0.9135836206068557D-2*t800
     &+0.1317311470595235D-1*t802+0.7362195534798311D-1*t806
     &-0.1380411662774683D-1*t809+0.7714850137356763D-1*t818
     &+0.1380411662774683D-1*t820-0.9135836206068557D-2*t822
     &-0.6586557352976174D-2*t824-0.7362195534798311D-1*t826
      t950 = t282*t694
      t959 = t67*t787*t790*t93
      t972 = -0.4771444294911944D0*t752*t284*t171+0.3333333333333333D0
     &*t2*t284-0.1D1*t847-0.1D1*t838-0.613516294566526D-2*t257*t950
     &-0.1D1*t834+0.4771444294911944D0*t284*t37*t170
     &-0.1646639338244043D-2*t688*t959-0.1063716611831605D-2*t772*t85
     &*t67*t93+0.2D1*t745*t284*t161-0.1D1*t188*t284*t201-t775
      t977 = t271*t93
      t978 = t196*t977
      t992 = t248*t93
      t995 = t309*t249
      t998 = 0.2300686104624472D-2*t683*t311-0.1097759558829362D-2
     &*t257*t959-0.2300686104624472D-2*t247*t978+0.613516294566526D-2
     &*t247*t950+0.1285808356226127D-1*t105*t815*t67*t93
     &+0.1285808356226127D-1*t105*t188*t281*t261-0.1522639367678093D-2
     &*t275*t797*t992-0.1227032589133052D-1*t257*t995-t786-t794+t801
     &+t803+t807
      t1009 = t709*t977
      t1022 = -t810+t819+0.1285808356226127D-1*t105*t90*t828*t261
     &+0.1285808356226127D-1*t257*t831*t261+0.2300686104624472D-2*t257
     &*t978+0.4601372209248945D-2*t413*t1009-t784
     &-0.6382299670989632D-2*t781*t55*t709*t93+t780
     &-0.2571616712452254D-1*t493*t282*t242+t771+0.1522639367678093D-2
     &*t795*t797*t992
      t1032 = t67*t455*sigmabb*t93
      t1045 = -t777-0.4286027854087091D-2*t648*t282*t238+t821-t823
     &-t825-t827-0.4601372209248945D-2*t493*t1009
     &+0.1227032589133052D-1*t247*t995-0.7668953682081574D-3*t648
     &*t1032+0.7668953682081574D-3*t661*t1032+0.2195519117658725D-2
     &*t247*t959+0.3191149835494816D-2*t275*t55*t196*t93
     &-0.2300686104624472D-2*t677*t311
      t1049 = 1.D0*t828*t67+1.D0*t831*t67-0.7777208750882749D2*t834
     &*t87*t506-0.7777208750882749D2*t838*t87*t506
     &+0.7777208750882749D2*t828*t87*t78*t205-0.7777208750882749D2
     &*t847*t87*t506-0.1391569437177703D2*t531*t853
     &+0.1930166210680909D2*t198*t131*t78*t539*t55*t204
     &+0.1391569437177703D2*t544*t853+t902*t103*t110-1.D0*t598*t306
     &+t223*t324-1.D0*t298*t226*t235+2.D0*t605*t235*t305-1.D0*t227
     &*t324*t234-1.D0*t227*t110*t935+t299*t265-1.D0*t227*t265*t305
     &+t104*(t972+t998+t1022+t1045)
      v2rhobsigmabb(i) = v2rhobsigmabb(i)+0.1554534543482745D-1*t764
     &+0.1554534543482745D-1*t138*t764+0.1554534543482745D-1*t14*t1049
      v2sigmabb2(i) = v2sigmabb2(i)
      endif
      endif ! rhob
      if ((rhoa.gt.tolmin).and.(rhob.gt.tolmin)) then
      t1 = 1/rhob
      t2 = t1**(1.D0/3.D0)
      trsb = 0.1285808356226127D-1*t2
      t3 = t1**(1.D0/15.D0)
      t4 = t3**2
      t5 = t4**2
      t7 = dexp(-0.5683671053580832D-1*t5)
      t10 = (0.942486901D0+0.349064173D0*t7)**2
      t11 = rhob**2
      t12 = rhob**(1.D0/3.D0)
      t13 = t12**2
      t16 = sigmabb/t13/t11
      t18 = sigmabb**2
      t19 = t11**2
      t23 = t18/t12/t19/rhob
      t26 = (0.1D1+0.6936084891727404D-1*t16+0.4389159297635699D-3*t23
     &)**2
      t34 = dexp(-0.4389159297635699D-3*t23)
      t35 = t34**2
      trcab = t10*t26/(0.1D1+0.9505299311782149D-1*t16/t2)*t35
      if(trsb.lt.tolmax*trcab) then
      t2 = 0.3141592653589793D1**2
      t5 = 1/rhob
      t6 = 1/0.3141592653589793D1*t5
      t7 = t6**(1.D0/3.D0)
      t9 = t6**(1.D0/15.D0)
      t10 = t9**2
      t11 = t10**2
      t13 = dexp(-0.7712625328179681D-1*t11)
      t16 = (0.942486901D0+0.349064173D0*t13)**2
      t19 = 0.3141592653589793D1**(1.D0/3.D0)
      t20 = t19**2
      t22 = 1/t20*sigmabb
      t23 = rhob**2
      t24 = rhob**(1.D0/3.D0)
      t25 = t24**2
      t27 = 1/t25/t23
      t32 = sigmabb**2
      t34 = t23**2
      t35 = t34*rhob
      t37 = 1/t24/t35
      t38 = 1/t19/0.3141592653589793D1*t32*t37
      t41 = (0.1D1+0.1487810599361293D0*t22*t27+0.2019518519390501D-2
     &*t38)**2
      t50 = dexp(0.4039037038781002D-2*t38)
      t54 = expei(-0.1858628590577086D0/t2*t7/t16/t41*(0.1D1
     &+0.1392138426088027D0*t22*t27/t7)*t50)
      t56 = t5**(1.D0/3.D0)
      t57 = t5**(1.D0/15.D0)
      t58 = t57**2
      t59 = t58**2
      t60 = 0.5683671053580832D-1*t59
      t61 = dexp(-t60)
      t63 = 0.942486901D0+0.349064173D0*t61
      t64 = t63**2
      t65 = 1/t64
      t66 = t56*t65
      t67 = sigmabb*t27
      t69 = t32*t37
      t71 = 0.1D1+0.6936084891727404D-1*t67+0.4389159297635699D-3*t69
      t72 = t71**2
      t73 = 1/t72
      t74 = 1/t56
      t77 = 0.1D1+0.9505299311782149D-1*t67*t74
      t78 = t73*t77
      t79 = 0.8778318595271399D-3*t69
      t80 = dexp(t79)
      t82 = t66*t78*t80
      t83 = dsqrt(t82)
      t86 = 0.6D1+0.4535739597862518D0*t83+0.5143233424904509D-1*t82
      t89 = 0.3D1+0.6803609396793777D0*t83+0.7714850137356763D-1*t82
      t90 = 1/t89
      t91 = t86*t90
      t92 = t66*t73
      t93 = t77*t54
      t97 = 0.1285808356226127D-1*t92*t93*t80+0.1D1
      t98 = t91*t97
      zk(i) = zk(i)+0.1554534543482745D-1*rhoa*(1.D0*t54+t98)
      vrhoa(i) = vrhoa(i)+0.1554534543482745D-1*t54
     &+0.1554534543482745D-1*t98
      t106 = t56**2
      t107 = 1/t106
      t108 = t107*t65
      t109 = t108*t73
      t110 = t77*t80
      t111 = 1/t23
      t113 = t109*t110*t111
      t114 = 0.4286027854087091D-2*t113
      t115 = t5**(1.D0/5.D0)
      t116 = t115**2
      t119 = 1/t64/t63
      t120 = 1/t116*t119
      t121 = t120*t73
      t122 = t77*t111
      t124 = dexp(t79-t60)
      t126 = t121*t122*t124
      t127 = 0.1360533322067624D-3*t126
      t129 = 1/t72/t71
      t130 = t66*t129
      t131 = t23*rhob
      t133 = 1/t25/t131
      t134 = sigmabb*t133
      t136 = t34*t23
      t138 = 1/t24/t136
      t139 = t32*t138
      t141 = -0.1849622637793975D0*t134-0.234088495873904D-2*t139
      t143 = t130*t110*t141
      t144 = 0.2571616712452254D-1*t143
      t148 = 1/t25/t34
      t149 = sigmabb*t148
      t151 = 1/t56/t5
      t154 = -0.2534746483141907D0*t134*t74+0.3168433103927383D-1*t149
     &*t151
      t155 = t73*t154
      t157 = t66*t155*t80
      t158 = 0.1285808356226127D-1*t157
      t160 = t92*t110*t139
      t161 = 0.601985888182142D-4*t160
      t162 = -t114-t127-t144+t158-t161
      t163 = t162*t54
      t165 = -t162
      t166 = t165*t74
      t167 = t166*t64
      t168 = 1/t77
      t169 = t72*t168
      t170 = 0.8778318595271399D-3*t69
      t171 = dexp(-t170)
      t172 = t169*t171
      t173 = t167*t172
      t175 = 1/t83
      t181 = -0.2306016315371713D-1*t113-0.7320092507805402D-3*t126
     &-0.1383609789223028D0*t143+0.6918048946115138D-1*t157
     &-0.3238871344356278D-3*t160
      t182 = t175*t181
      t189 = 0.3278192763011299D1*t182-0.1714411141634836D-1*t113
     &-0.5442133288270495D-3*t126-0.1028646684980902D0*t143
     &+0.5143233424904509D-1*t157-0.2407943552728568D-3*t160
      t190 = t189*t90
      t191 = t190*t97
      t192 = t89**2
      t193 = 1/t192
      t194 = t86*t193
      t201 = 0.4917289144516948D1*t182-0.2571616712452254D-1*t113
     &-0.8163199932405743D-3*t126-0.1542970027471353D0*t143
     &+0.7714850137356763D-1*t157-0.3611915329092852D-3*t160
      t202 = t97*t201
      t203 = t194*t202
      t205 = t111*t80
      t209 = t111*t124
      t213 = t141*t80
      t217 = t154*t54
      t221 = t66*t78
      t222 = t54*t32
      t223 = t138*t80
      t224 = t222*t223
      t228 = t54*t80
      t232 = -0.4286027854087091D-2*t109*t93*t205
     &-0.1360533322067624D-3*t121*t93*t209-0.2571616712452254D-1*t130
     &*t93*t213+0.1285808356226127D-1*t92*t217*t80
     &-0.601985888182142D-4*t221*t224+0.1285808356226127D-1*t92*t77
     &*t162*t228-t114-t127-t144+t158-t161
      t233 = t91*t232
      vrhob(i) = vrhob(i)+0.1554534543482745D-1*rhoa*(1.D0*t163
     &-0.7777208750882749D2*t173+t191-1.D0*t203+t233)
      t239 = sigmabb*t37
      t241 = 0.6936084891727404D-1*t27+0.8778318595271399D-3*t239
      t243 = t130*t110*t241
      t244 = 0.2571616712452254D-1*t243
      t245 = t65*t73
      t246 = t27*t80
      t247 = t245*t246
      t248 = 0.1222199328351994D-2*t247
      t250 = t92*t110*t239
      t251 = 0.2257447080683033D-4*t250
      t252 = -t244+t248+t251
      t253 = t252*t54
      t255 = -t252
      t256 = t255*t74
      t257 = t256*t64
      t258 = t257*t172
      t263 = -0.1383609789223028D0*t243+0.6575812588638345D-2*t247
     &+0.1214576754133604D-3*t250
      t264 = t175*t263
      t269 = 0.3278192763011299D1*t264-0.1028646684980902D0*t243
     &+0.4888797313407978D-2*t247+0.9029788322732131D-4*t250
      t270 = t269*t90
      t271 = t270*t97
      t276 = 0.4917289144516948D1*t264-0.1542970027471353D0*t243
     &+0.7333195970111966D-2*t247+0.135446824840982D-3*t250
      t277 = t97*t276
      t278 = t194*t277
      t280 = t241*t80
      t288 = t54*sigmabb
      t289 = t37*t80
      t290 = t288*t289
      t297 = -0.2571616712452254D-1*t130*t93*t280
     &+0.1222199328351994D-2*t245*t27*t54*t80+0.2257447080683033D-4
     &*t221*t290+0.1285808356226127D-1*t92*t77*t252*t228-t244+t248+t251
      t298 = t91*t297
      vsigmabb(i) = vsigmabb(i)+0.1554534543482745D-1*rhoa*(1.D0*t253
     &-0.7777208750882749D2*t258+t271-1.D0*t278+t298)
      v2rhoa2(i) = v2rhoa2(i)
      v2rhoab(i) = v2rhoab(i)+0.1554534543482745D-1*t163
     &-0.1208993965512352D1*t173+0.1554534543482745D-1*t191
     &-0.1554534543482745D-1*t203+0.1554534543482745D-1*t233
      t311 = t121*t154*t111*t124
      t312 = 0.2721066644135248D-3*t311
      t316 = 1/t24/t34/t131
      t317 = t32*t316
      t319 = 0.678194967191124D0*t149+0.1482560473868058D-1*t317
      t321 = t130*t110*t319
      t322 = 0.2571616712452254D-1*t321
      t323 = t129*t77
      t324 = t66*t323
      t326 = t324*t213*t139
      t327 = 0.2407943552728568D-3*t326
      t342 = 0.9294070438186991D0*t149*t74-0.2323517609546748D0
     &*sigmabb/t25/t35*t151+0.4224577471903178D-1*sigmabb/t25/t136/t56
     &/t111
      t345 = t66*t73*t342*t80
      t346 = 0.1285808356226127D-1*t345
      t347 = t32**2
      t348 = t34**2
      t351 = 1/t25/t348/t34
      t354 = t92*t110*t347*t351
      t355 = 0.2818359422037475D-6*t354
      t357 = t92*t110*t317
      t358 = 0.3812577291820233D-3*t357
      t359 = t154*t80
      t361 = t92*t359*t139
      t362 = 0.1203971776364284D-3*t361
      t364 = t130*t359*t141
      t365 = 0.5143233424904509D-1*t364
      t367 = 1/t58/t5
      t368 = t64**2
      t371 = t367/t368*t73
      t372 = 1/t34
      t373 = t77*t372
      t376 = dexp(t79-0.1136734210716166D0*t59)
      t378 = t371*t373*t376
      t379 = 0.2159401412535421D-5*t378
      t381 = 1/t24/t348
      t383 = t32*t124
      t385 = t121*t77*t381*t383
      t386 = 0.1273940795796543D-5*t385
      t388 = t367*t119*t73
      t389 = t373*t124
      t390 = t388*t389
      t391 = 0.2062086362684779D-5*t390
      t392 = t72**2
      t393 = 1/t392
      t394 = t66*t393
      t395 = t141**2
      t397 = t394*t110*t395
      t398 = 0.7714850137356763D-1*t397
      t399 = 1/t131
      t402 = t121*t77*t399*t124
      t403 = 0.2721066644135248D-3*t402
      t404 = t120*t129
      t405 = t141*t124
      t407 = t404*t122*t405
      t408 = 0.5442133288270495D-3*t407
      t410 = t109*t110*t399
      t411 = 0.8572055708174181D-2*t410
      t415 = 1/t106/t5*t65*t73
      t417 = t415*t110*t372
      t418 = 0.2857351902724727D-2*t417
      t422 = 1/t116/t5*t119*t73
      t423 = t422*t389
      t424 = 0.9070222147117492D-5*t423
      t427 = t109*t110*t381*t32
      t428 = 0.4013239254547614D-4*t427
      t429 = t108*t129
      t432 = t429*t110*t111*t141
      t433 = 0.1714411141634836D-1*t432
      t435 = t109*t359*t111
      t436 = 0.8572055708174181D-2*t435
      t437 = -t312-t322+t327+t346+t355+t358-t362-t365+t379+t386-t391
     &+t398+t403+t408+t411-t418-t424+t428+t433-t436
      t440 = t162**2
      t443 = t162*t165
      t445 = t64*t72
      t447 = t445*t168*t171
      t455 = t165**2
      t471 = dexp(-t170-t60)
      t476 = t71*t168
      t481 = t77**2
      t483 = t72/t481
      t488 = t166*t445
      t505 = -0.1088426657654099D-2*t311-0.1028646684980902D0*t321
     &+0.9631774210914273D-3*t326+0.5143233424904509D-1*t345
     &+0.112734376881499D-5*t354+0.1525030916728093D-2*t357
     &-0.4815887105457136D-3*t361-0.2057293369961804D0*t364
     &+0.8637605650141684D-5*t378+0.5095763183186173D-5*t385
     &-0.8248345450739115D-5*t390
      t516 = 1/t83/t82
      t517 = t181**2
      t518 = t516*t517
      t540 = -0.146401850156108D-2*t311-0.1383609789223028D0*t321
     &+0.1295548537742511D-2*t326+0.6918048946115138D-1*t345
     &+0.1516365042658901D-5*t354+0.2051285184758976D-2*t357
     &-0.6477742688712557D-3*t361-0.2767219578446055D0*t364
     &+0.1161825134662837D-4*t378+0.6854197779240001D-5*t385
     &-0.1109466610563933D-4*t390+0.4150829367669083D0*t397
     &+0.146401850156108D-2*t402+0.2928037003122161D-2*t407
     &+0.4612032630743426D-1*t410-0.1537344210247809D-1*t417
     &-0.4880061671870268D-4*t423+0.2159247562904186D-3*t427
     &+0.9224065261486851D-1*t432-0.4612032630743426D-1*t435
      t541 = t175*t540
      t543 = 0.3085940054942705D0*t397+0.1088426657654099D-2*t402
     &+0.2176853315308198D-2*t407+0.3428822283269673D-1*t410
     &-0.1142940761089891D-1*t417-0.3628088858846997D-4*t423
     &+0.1605295701819045D-3*t427+0.6857644566539345D-1*t432
     &-0.3428822283269673D-1*t435-0.2369304401100098D2*t518
     &+0.3278192763011299D1*t541
      t547 = t189*t193
      t554 = t86/t192/t89
      t555 = t201**2
      t573 = -0.1632639986481149D-2*t311-0.1542970027471353D0*t321
     &+0.1444766131637141D-2*t326+0.7714850137356763D-1*t345
     &+0.1691015653222485D-5*t354+0.228754637509214D-2*t357
     &-0.7223830658185704D-3*t361-0.3085940054942705D0*t364
     &+0.1295640847521253D-4*t378+0.764364477477926D-5*t385
     &-0.1237251817610867D-4*t390
      t585 = 0.4628910082414058D0*t397+0.1632639986481149D-2*t402
     &+0.3265279972962297D-2*t407+0.5143233424904509D-1*t410
     &-0.1714411141634836D-1*t417-0.5442133288270495D-4*t423
     &+0.2407943552728568D-3*t427+0.1028646684980902D0*t432
     &-0.5143233424904509D-1*t435-0.3553956601650148D2*t518
     &+0.4917289144516948D1*t541
      t590 = 1/t71
      t598 = t59**2
      t600 = 1/t598/t58/t57
      t602 = t600/t63
      t613 = t93*t372*t124
      t628 = t120*t78
      t632 = t54*t381
      t636 = t54*t141
      t637 = t139*t80
      t641 = 0.2D1*t590*t165*t141+0.4681769917478079D-2*t165*t32*t138
     &+0.105811516582519D-1*t602*t165*t111*t61+0.2721066644135248D-3
     &*t121*t93*t399*t124-0.2D1*t443-0.2062086362684779D-5*t388*t613
     &+0.1285808356226127D-1*t92*t77*t437*t228-0.2571616712452254D-1
     &*t130*t93*t319*t80+0.7714850137356763D-1*t394*t93*t395*t80
     &-0.2721066644135248D-3*t628*t163*t209+0.1273940795796543D-5*t628
     &*t632*t383+0.2407943552728568D-3*t324*t636*t637+t327
      t646 = t108*t323
      t647 = t54*t111
      t667 = 0.1285808356226127D-1*t92*t342*t54*t80
     &+0.1714411141634836D-1*t646*t647*t213-t322-0.2721066644135248D-3
     &*t121*t217*t209-t312-0.5143233424904509D-1*t130*t217*t213
     &+0.2159401412535421D-5*t371*t93*t372*t376-0.5143233424904509D-1
     &*t324*t163*t213+t379-0.1D1*t154*t165*t168-t365+t358-t362+t346
      t678 = t108*t78
      t679 = t32*t80
      t683 = t120*t323
      t691 = t355+t408+t398-0.8572055708174181D-2*t109*t217*t205
     &-0.2857351902724727D-2*t415*t93*t372*t80+t403-t391
     &-0.9070222147117492D-5*t422*t613+0.4013239254547614D-4*t678*t632
     &*t679+t386+0.5442133288270495D-3*t683*t647*t405
     &+0.1285808356226127D-1*t92*t77*t440*t228+t428
      t718 = t66*t155
      t721 = t433-t424-t418+t411-t436+0.8572055708174181D-2*t109*t93
     &*t399*t80-0.8572055708174181D-2*t678*t163*t205
     &+0.2818359422037475D-6*t221*t54*t347*t351*t80
     &+0.2571616712452254D-1*t92*t154*t162*t228-0.1203971776364284D-3
     &*t221*t163*t637-0.1D1*t455+0.3333333333333333D0*t5*t165
     &+0.3812577291820233D-3*t221*t222*t316*t80-0.1203971776364284D-3
     &*t718*t224
      t725 = 1.D0*t437*t54+1.D0*t440*t54-0.155544175017655D3*t443*t74
     &*t447+0.7777208750882749D2*t437*t74*t64*t172
     &-0.7777208750882749D2*t455*t74*t64*t172-0.2592402916960916D2
     &*t165*t151*t64*t169*t111*t171-0.8229182527097421D0*t165/t57/t5
     &*t63*t169*t111*t471-0.155544175017655D3*t167*t476*t141*t171
     &+0.7777208750882749D2*t167*t483*t154*t171-0.3641110197183013D0
     &*t488*t168*t32*t138*t171+(t505+t543)*t90*t97-2.D0*t547*t202+2.D0
     &*t190*t232+2.D0*t554*t97*t555-2.D0*t194*t232*t201-1.D0*t194*t97*
     &(t573+t585)+t91*(t641+t667+t691+t721)
      v2rhob2(i) = v2rhob2(i)+0.1554534543482745D-1*rhoa*t725
      v2rhoasigmabb(i) = v2rhoasigmabb(i)+0.1554534543482745D-1*t253
     &-0.1208993965512352D1*t258+0.1554534543482745D-1*t271
     &-0.1554534543482745D-1*t278+0.1554534543482745D-1*t298
      t737 = t429*t110*t111*t241
      t738 = 0.8572055708174181D-2*t737
      t739 = t133*t65
      t741 = t739*t73*t80
      t742 = 0.4073997761173315D-3*t741
      t745 = t109*t110*t316*sigmabb
      t746 = 0.7524823602276776D-5*t745
      t747 = t241*t124
      t749 = t404*t122*t747
      t750 = 0.2721066644135248D-3*t749
      t751 = t600*t119
      t754 = t751*t73*t148*t124
      t755 = 0.1293227644990607D-4*t754
      t757 = sigmabb*t124
      t759 = t121*t77*t316*t757
      t760 = 0.2388638992118519D-6*t759
      t763 = t394*t110*t141*t241
      t764 = 0.7714850137356763D-1*t763
      t765 = t65*t129
      t767 = t765*t246*t141
      t768 = 0.2444398656703989D-2*t767
      t770 = t324*t213*t239
      t771 = 0.4514894161366065D-4*t770
      t773 = sigmabb*t138
      t775 = -0.1849622637793975D0*t133-0.4681769917478079D-2*t773
      t777 = t130*t110*t775
      t778 = 0.2571616712452254D-1*t777
      t780 = t130*t359*t241
      t781 = 0.2571616712452254D-1*t780
      t786 = -0.2534746483141907D0*t133*t74+0.3168433103927383D-1*t148
     &*t151
      t789 = t66*t73*t786*t80
      t790 = 0.1285808356226127D-1*t789
      t792 = t92*t359*t239
      t793 = 0.2257447080683033D-4*t792
      t794 = t138*t241
      t796 = t324*t679*t794
      t797 = 0.1203971776364284D-3*t796
      t799 = 1/t348/rhob
      t802 = t245*t799*t80*t32
      t803 = 0.5722056048640281D-5*t802
      t804 = t32*sigmabb
      t807 = 1/t25/t348/t131
      t810 = t92*t110*t804*t807
      t811 = 0.1056884783264053D-6*t810
      t813 = t92*t110*t773
      t814 = 0.1203971776364284D-3*t813
      t815 = t738-t742-t746+t750-t755-t760+t764-t768-t771-t778-t781
     &+t790+t793+t797-t803-t811-t814
      t818 = t162*t252
      t821 = t162*t255
      t825 = t252*t165
      t834 = t165*t255
      t839 = t476*t241*t171
      t845 = t483*t27*t171
      t850 = t168*sigmabb*t37*t171
      t854 = t516*t181*t263
      t873 = 0.4612032630743426D-1*t737-0.2191937529546115D-2*t741
     &-0.4048589180445348D-4*t745+0.146401850156108D-2*t749
     &-0.6957967027662436D-4*t754-0.12851620836075D-5*t759
     &+0.4150829367669083D0*t763-0.1315162517727669D-1*t767
     &-0.2429153508267209D-3*t770-0.1383609789223028D0*t777
     &-0.1383609789223028D0*t780+0.6918048946115138D-1*t789
     &+0.1214576754133604D-3*t792+0.6477742688712557D-3*t796
     &-0.3078644156046066D-4*t802-0.5686368909970879D-6*t810
     &-0.6477742688712557D-3*t813
      t874 = t175*t873
      t893 = -0.2369304401100098D2*t854+0.3278192763011299D1*t874
     &+0.3428822283269673D-1*t737-0.1629599104469326D-2*t741
     &-0.300992944091071D-4*t745+0.1088426657654099D-2*t749
     &-0.5172910579962427D-4*t754-0.9554555968474075D-6*t759
     &+0.3085940054942705D0*t763-0.9777594626815955D-2*t767
     &-0.1805957664546426D-3*t770-0.1028646684980902D0*t777
     &-0.1028646684980902D0*t780+0.5143233424904509D-1*t789
     &+0.9029788322732131D-4*t792+0.4815887105457136D-3*t796
     &-0.2288822419456112D-4*t802-0.4227539133056213D-6*t810
     &-0.4815887105457136D-3*t813
      t899 = t269*t193
      t927 = -0.3553956601650148D2*t854+0.4917289144516948D1*t874
     &+0.5143233424904509D-1*t737-0.2444398656703989D-2*t741
     &-0.4514894161366065D-4*t745+0.1632639986481149D-2*t749
     &-0.775936586994364D-4*t754-0.1433183395271111D-5*t759
     &+0.4628910082414058D0*t763-0.1466639194022393D-1*t767
     &-0.2708936496819639D-3*t770-0.1542970027471353D0*t777
     &-0.1542970027471353D0*t780+0.7714850137356763D-1*t789
     &+0.135446824840982D-3*t792+0.7223830658185704D-3*t796
     &-0.3433233629184168D-4*t802-0.634130869958432D-6*t810
     &-0.7223830658185704D-3*t813
      t941 = t590*t255
      t951 = t54*t316
      t974 = t245*t27
      t978 = 0.4681769917478079D-2*t255*t32*t138-0.1D1*t154*t255*t168
     &+0.2D1*t941*t141+0.3333333333333333D0*t5*t255
     &-0.1293227644990607D-4*t751*t73*t148*t54*t124
     &-0.7524823602276776D-5*t678*t951*sigmabb*t80-t742
     &+0.8572055708174181D-2*t646*t647*t280-0.4073997761173315D-3*t739
     &*t73*t54*t80-0.1360533322067624D-3*t628*t253*t209
     &-0.2388638992118519D-6*t628*t951*t757+0.7714850137356763D-1*t66
     &*t393*t77*t636*t280+0.1222199328351994D-2*t974*t163*t80
      t990 = -0.1D1*t834-0.1D1*t821-0.1D1*t825+t738
     &-0.5722056048640281D-5*t245*t799*t222*t80+t790+t793-t746
     &+0.105811516582519D-1*t602*t255*t111*t61-t768+t764-t760-t755
      t1009 = -0.1203971776364284D-3*t221*t288*t223+t750+t797-t811
     &-t814-t781-t771-t778-t803+0.2257447080683033D-4*t718*t290
     &-0.2571616712452254D-1*t130*t217*t280+0.1285808356226127D-1*t92
     &*t786*t54*t80-0.1056884783264053D-6*t221*t54*t804*t807*t80
      t1017 = t765*t27
      t1027 = t239*t80
      t1055 = 0.1285808356226127D-1*t92*t77*t815*t228
     &-0.2571616712452254D-1*t324*t253*t213-0.2444398656703989D-2
     &*t1017*t636*t80-0.601985888182142D-4*t221*t253*t637
     &-0.2571616712452254D-1*t324*t163*t280-0.4514894161366065D-4*t324
     &*t636*t1027+0.1203971776364284D-3*t324*t222*t794*t80
     &+0.1285808356226127D-1*t92*t154*t252*t228+0.2257447080683033D-4
     &*t221*t163*t1027+0.2721066644135248D-3*t683*t647*t747
     &+0.1285808356226127D-1*t221*t818*t228-0.2571616712452254D-1*t130
     &*t93*t775*t80-0.4286027854087091D-2*t678*t253*t205
      t1059 = 1.D0*t815*t54+1.D0*t818*t54-0.7777208750882749D2*t821
     &*t74*t447-0.7777208750882749D2*t825*t74*t447
     &+0.7777208750882749D2*t815*t74*t64*t172-0.7777208750882749D2
     &*t834*t74*t447-0.155544175017655D3*t167*t839
     &+0.7392469698735191D1*t165*t107*t64*t845+0.136541632394363D0
     &*t488*t850+t893*t90*t97-1.D0*t547*t277+t190*t297-1.D0*t899*t202
     &+2.D0*t554*t202*t276-1.D0*t194*t297*t201-1.D0*t194*t97*t927+t270
     &*t232-1.D0*t194*t232*t276+t91*(t978+t990+t1009+t1055)
      v2rhobsigmabb(i) = v2rhobsigmabb(i)+0.1554534543482745D-1*rhoa
     &*t1059
      t1063 = t241**2
      t1065 = t394*t110*t1063
      t1066 = 0.7714850137356763D-1*t1065
      t1068 = t765*t246*t241
      t1069 = 0.4888797313407978D-2*t1068
      t1071 = t324*t280*t239
      t1072 = 0.9029788322732131D-4*t1071
      t1073 = t110*t37
      t1074 = t130*t1073
      t1075 = 0.2257447080683033D-4*t1074
      t1076 = 1/t348
      t1079 = t245*t1076*t80*sigmabb
      t1080 = 0.4291542036480211D-5*t1079
      t1083 = 1/t25/t348/t23
      t1086 = t92*t110*t32*t1083
      t1087 = 0.39633179372402D-7*t1086
      t1088 = t92*t1073
      t1089 = 0.2257447080683033D-4*t1088
      t1090 = t1066-t1069-t1072-t1075+t1080+t1087+t1089
      t1093 = t252**2
      t1096 = t252*t255
      t1105 = t255**2
      t1119 = t263**2
      t1120 = t516*t1119
      t1130 = t175*(0.4150829367669083D0*t1065-0.2630325035455338D-1
     &*t1068-0.4858307016534418D-3*t1071-0.1214576754133604D-3*t1074
     &+0.2308983117034549D-4*t1079+0.2132388341239079D-6*t1086
     &+0.1214576754133604D-3*t1088)
      t1146 = t276**2
      t1170 = t93*t289
      t1173 = t54*t241
      t1207 = 0.7714850137356763D-1*t394*t93*t1063*t80
     &-0.2257447080683033D-4*t130*t1170-0.9029788322732131D-4*t324
     &*t1173*t1027-0.5143233424904509D-1*t324*t253*t280
     &-0.4888797313407978D-2*t1017*t1173*t80+0.2D1*t941*t241
     &-0.9505299311782149D-1*t27*t255*t74*t168+0.4291542036480211D-5
     &*t245*t1076*t288*t80+0.2444398656703989D-2*t974*t253*t80
     &+0.1285808356226127D-1*t92*t77*t1090*t228+0.39633179372402D-7
     &*t221*t222*t1083*t80-0.175566371905428D-2*t255*sigmabb*t37
      t1219 = 0.4514894161366065D-4*t221*t253*t1027
     &+0.2257447080683033D-4*t92*t1170+t1066-t1069-t1072-t1075+t1080
     &+t1087+t1089+0.1285808356226127D-1*t92*t77*t1093*t228-0.2D1
     &*t1096-0.1D1*t1105
      s1 = 1.D0*t1090*t54+1.D0*t1093*t54-0.155544175017655D3*t1096*t74
     &*t447+0.7777208750882749D2*t1090*t74*t64*t172
     &-0.7777208750882749D2*t1105*t74*t64*t172-0.155544175017655D3
     &*t257*t839+0.7392469698735191D1*t255*t107*t64*t845
      t1222 = s1+0.136541632394363D0*t256*t445*t850+(
     &-0.2369304401100098D2*t1120+0.3278192763011299D1*t1130
     &+0.3085940054942705D0*t1065-0.1955518925363191D-1*t1068
     &-0.3611915329092852D-3*t1071-0.9029788322732131D-4*t1074
     &+0.1716616814592084D-4*t1079+0.158532717489608D-6*t1086
     &+0.9029788322732131D-4*t1088)*t90*t97-2.D0*t899*t277+2.D0*t270
     &*t297+2.D0*t554*t97*t1146-2.D0*t194*t297*t276-1.D0*t194*t97*(
     &-0.3553956601650148D2*t1120+0.4917289144516948D1*t1130
     &+0.4628910082414058D0*t1065-0.2933278388044787D-1*t1068
     &-0.5417872993639278D-3*t1071-0.135446824840982D-3*t1074
     &+0.2574925221888126D-4*t1079+0.237799076234412D-6*t1086
     &+0.135446824840982D-3*t1088)+t91*(t1207+t1219)
      v2sigmabb2(i) = v2sigmabb2(i)+0.1554534543482745D-1*rhoa*t1222
      endif
      t1 = 1/rhoa
      t2 = t1**(1.D0/3.D0)
      trsa = 0.1285808356226127D-1*t2
      t3 = t1**(1.D0/15.D0)
      t4 = t3**2
      t5 = t4**2
      t7 = dexp(-0.5683671053580832D-1*t5)
      t10 = (0.942486901D0+0.349064173D0*t7)**2
      t11 = rhoa**2
      t12 = rhoa**(1.D0/3.D0)
      t13 = t12**2
      t16 = sigmaaa/t13/t11
      t18 = sigmaaa**2
      t19 = t11**2
      t23 = t18/t12/t19/rhoa
      t26 = (0.1D1+0.6936084891727404D-1*t16+0.4389159297635699D-3*t23
     &)**2
      t34 = dexp(-0.4389159297635699D-3*t23)
      t35 = t34**2
      trcba = t10*t26/(0.1D1+0.9505299311782149D-1*t16/t2)*t35
      if(trsa.lt.tolmax*trcba) then
      t2 = 0.3141592653589793D1**2
      t5 = 1/rhoa
      t6 = 1/0.3141592653589793D1*t5
      t7 = t6**(1.D0/3.D0)
      t9 = t6**(1.D0/15.D0)
      t10 = t9**2
      t11 = t10**2
      t13 = dexp(-0.7712625328179681D-1*t11)
      t16 = (0.942486901D0+0.349064173D0*t13)**2
      t19 = 0.3141592653589793D1**(1.D0/3.D0)
      t20 = t19**2
      t22 = 1/t20*sigmaaa
      t23 = rhoa**2
      t24 = rhoa**(1.D0/3.D0)
      t25 = t24**2
      t27 = 1/t25/t23
      t32 = sigmaaa**2
      t34 = t23**2
      t35 = t34*rhoa
      t37 = 1/t24/t35
      t38 = 1/t19/0.3141592653589793D1*t32*t37
      t41 = (0.1D1+0.1487810599361293D0*t22*t27+0.2019518519390501D-2
     &*t38)**2
      t50 = dexp(0.4039037038781002D-2*t38)
      t54 = expei(-0.1858628590577086D0/t2*t7/t16/t41*(0.1D1
     &+0.1392138426088027D0*t22*t27/t7)*t50)
      t56 = t5**(1.D0/3.D0)
      t57 = t5**(1.D0/15.D0)
      t58 = t57**2
      t59 = t58**2
      t60 = 0.5683671053580832D-1*t59
      t61 = dexp(-t60)
      t63 = 0.942486901D0+0.349064173D0*t61
      t64 = t63**2
      t65 = 1/t64
      t66 = t56*t65
      t67 = sigmaaa*t27
      t69 = t32*t37
      t71 = 0.1D1+0.6936084891727404D-1*t67+0.4389159297635699D-3*t69
      t72 = t71**2
      t73 = 1/t72
      t74 = 1/t56
      t77 = 0.1D1+0.9505299311782149D-1*t67*t74
      t78 = t73*t77
      t79 = 0.8778318595271399D-3*t69
      t80 = dexp(t79)
      t82 = t66*t78*t80
      t83 = dsqrt(t82)
      t86 = 0.6D1+0.4535739597862518D0*t83+0.5143233424904509D-1*t82
      t89 = 0.3D1+0.6803609396793777D0*t83+0.7714850137356763D-1*t82
      t90 = 1/t89
      t91 = t86*t90
      t92 = t66*t73
      t93 = t77*t54
      t97 = 0.1285808356226127D-1*t92*t93*t80+0.1D1
      t98 = t91*t97
      zk(i) = zk(i)+0.1554534543482745D-1*rhob*(1.D0*t54+t98)
      t103 = t56**2
      t104 = 1/t103
      t105 = t104*t65
      t106 = t105*t73
      t107 = t77*t80
      t108 = 1/t23
      t110 = t106*t107*t108
      t111 = 0.4286027854087091D-2*t110
      t112 = t5**(1.D0/5.D0)
      t113 = t112**2
      t116 = 1/t64/t63
      t117 = 1/t113*t116
      t118 = t117*t73
      t119 = t77*t108
      t121 = dexp(t79-t60)
      t123 = t118*t119*t121
      t124 = 0.1360533322067624D-3*t123
      t126 = 1/t72/t71
      t127 = t66*t126
      t128 = t23*rhoa
      t130 = 1/t25/t128
      t131 = sigmaaa*t130
      t133 = t34*t23
      t135 = 1/t24/t133
      t136 = t32*t135
      t138 = -0.1849622637793975D0*t131-0.234088495873904D-2*t136
      t140 = t127*t107*t138
      t141 = 0.2571616712452254D-1*t140
      t145 = 1/t25/t34
      t146 = sigmaaa*t145
      t148 = 1/t56/t5
      t151 = -0.2534746483141907D0*t131*t74+0.3168433103927383D-1*t146
     &*t148
      t152 = t73*t151
      t154 = t66*t152*t80
      t155 = 0.1285808356226127D-1*t154
      t157 = t92*t107*t136
      t158 = 0.601985888182142D-4*t157
      t159 = -t111-t124-t141+t155-t158
      t160 = t159*t54
      t162 = -t159
      t163 = t162*t74
      t164 = t163*t64
      t165 = 1/t77
      t166 = t72*t165
      t167 = 0.8778318595271399D-3*t69
      t168 = dexp(-t167)
      t169 = t166*t168
      t170 = t164*t169
      t172 = 1/t83
      t178 = -0.2306016315371713D-1*t110-0.7320092507805402D-3*t123
     &-0.1383609789223028D0*t140+0.6918048946115138D-1*t154
     &-0.3238871344356278D-3*t157
      t179 = t172*t178
      t186 = 0.3278192763011299D1*t179-0.1714411141634836D-1*t110
     &-0.5442133288270495D-3*t123-0.1028646684980902D0*t140
     &+0.5143233424904509D-1*t154-0.2407943552728568D-3*t157
      t187 = t186*t90
      t188 = t187*t97
      t189 = t89**2
      t190 = 1/t189
      t191 = t86*t190
      t198 = 0.4917289144516948D1*t179-0.2571616712452254D-1*t110
     &-0.8163199932405743D-3*t123-0.1542970027471353D0*t140
     &+0.7714850137356763D-1*t154-0.3611915329092852D-3*t157
      t199 = t97*t198
      t200 = t191*t199
      t202 = t108*t80
      t206 = t108*t121
      t210 = t138*t80
      t214 = t151*t54
      t218 = t66*t78
      t219 = t54*t32
      t220 = t135*t80
      t221 = t219*t220
      t225 = t54*t80
      t229 = -0.4286027854087091D-2*t106*t93*t202
     &-0.1360533322067624D-3*t118*t93*t206-0.2571616712452254D-1*t127
     &*t93*t210+0.1285808356226127D-1*t92*t214*t80
     &-0.601985888182142D-4*t218*t221+0.1285808356226127D-1*t92*t77
     &*t159*t225-t111-t124-t141+t155-t158
      t230 = t91*t229
      vrhoa(i) = vrhoa(i)+0.1554534543482745D-1*rhob*(1.D0*t160
     &-0.7777208750882749D2*t170+t188-1.D0*t200+t230)
      vrhob(i) = vrhob(i)+0.1554534543482745D-1*t54
     &+0.1554534543482745D-1*t98
      t239 = sigmaaa*t37
      t241 = 0.6936084891727404D-1*t27+0.8778318595271399D-3*t239
      t243 = t127*t107*t241
      t244 = 0.2571616712452254D-1*t243
      t245 = t65*t73
      t246 = t27*t80
      t247 = t245*t246
      t248 = 0.1222199328351994D-2*t247
      t250 = t92*t107*t239
      t251 = 0.2257447080683033D-4*t250
      t252 = -t244+t248+t251
      t253 = t252*t54
      t255 = -t252
      t256 = t255*t74
      t257 = t256*t64
      t258 = t257*t169
      t263 = -0.1383609789223028D0*t243+0.6575812588638345D-2*t247
     &+0.1214576754133604D-3*t250
      t264 = t172*t263
      t269 = 0.3278192763011299D1*t264-0.1028646684980902D0*t243
     &+0.4888797313407978D-2*t247+0.9029788322732131D-4*t250
      t270 = t269*t90
      t271 = t270*t97
      t276 = 0.4917289144516948D1*t264-0.1542970027471353D0*t243
     &+0.7333195970111966D-2*t247+0.135446824840982D-3*t250
      t277 = t97*t276
      t278 = t191*t277
      t280 = t241*t80
      t288 = t54*sigmaaa
      t289 = t37*t80
      t290 = t288*t289
      t297 = -0.2571616712452254D-1*t127*t93*t280
     &+0.1222199328351994D-2*t245*t27*t54*t80+0.2257447080683033D-4
     &*t218*t290+0.1285808356226127D-1*t92*t77*t252*t225-t244+t248+t251
      t298 = t91*t297
      vsigmaaa(i) = vsigmaaa(i)+0.1554534543482745D-1*rhob*(1.D0*t253
     &-0.7777208750882749D2*t258+t271-1.D0*t278+t298)
      t305 = 1/t24/t34/t128
      t306 = t32*t305
      t308 = t92*t107*t306
      t309 = 0.3812577291820233D-3*t308
      t310 = t32**2
      t311 = t34**2
      t314 = 1/t25/t311/t34
      t317 = t92*t107*t310*t314
      t318 = 0.2818359422037475D-6*t317
      t319 = t151*t80
      t321 = t92*t319*t136
      t322 = 0.1203971776364284D-3*t321
      t325 = 0.678194967191124D0*t146+0.1482560473868058D-1*t306
      t327 = t127*t107*t325
      t328 = 0.2571616712452254D-1*t327
      t329 = t126*t77
      t330 = t66*t329
      t332 = t330*t210*t136
      t333 = 0.2407943552728568D-3*t332
      t348 = 0.9294070438186991D0*t146*t74-0.2323517609546748D0
     &*sigmaaa/t25/t35*t148+0.4224577471903178D-1*sigmaaa/t25/t133/t56
     &/t108
      t351 = t66*t73*t348*t80
      t352 = 0.1285808356226127D-1*t351
      t353 = 1/t128
      t356 = t118*t77*t353*t121
      t357 = 0.2721066644135248D-3*t356
      t359 = 1/t58/t5
      t360 = t64**2
      t363 = t359/t360*t73
      t364 = 1/t34
      t365 = t77*t364
      t368 = dexp(t79-0.1136734210716166D0*t59)
      t370 = t363*t365*t368
      t371 = 0.2159401412535421D-5*t370
      t373 = 1/t24/t311
      t375 = t32*t121
      t377 = t118*t77*t373*t375
      t378 = 0.1273940795796543D-5*t377
      t380 = t359*t116*t73
      t381 = t365*t121
      t382 = t380*t381
      t383 = 0.2062086362684779D-5*t382
      t384 = t72**2
      t385 = 1/t384
      t386 = t66*t385
      t387 = t138**2
      t389 = t386*t107*t387
      t390 = 0.7714850137356763D-1*t389
      t392 = t127*t319*t138
      t393 = 0.5143233424904509D-1*t392
      t396 = t118*t151*t108*t121
      t397 = 0.2721066644135248D-3*t396
      t399 = t106*t107*t353
      t400 = 0.8572055708174181D-2*t399
      t404 = 1/t103/t5*t65*t73
      t406 = t404*t107*t364
      t407 = 0.2857351902724727D-2*t406
      t411 = 1/t113/t5*t116*t73
      t412 = t411*t381
      t413 = 0.9070222147117492D-5*t412
      t416 = t106*t107*t373*t32
      t417 = 0.4013239254547614D-4*t416
      t418 = t117*t126
      t419 = t138*t121
      t421 = t418*t119*t419
      t422 = 0.5442133288270495D-3*t421
      t423 = t105*t126
      t426 = t423*t107*t108*t138
      t427 = 0.1714411141634836D-1*t426
      t429 = t106*t319*t108
      t430 = 0.8572055708174181D-2*t429
      t431 = t309+t318-t322-t328+t333+t352+t357+t371+t378-t383+t390
     &-t393-t397+t400-t407-t413+t417+t422+t427-t430
      t434 = t159**2
      t437 = t159*t162
      t439 = t64*t72
      t441 = t439*t165*t168
      t449 = t162**2
      t465 = dexp(-t167-t60)
      t470 = t71*t165
      t475 = t77**2
      t477 = t72/t475
      t482 = t163*t439
      t499 = 0.1525030916728093D-2*t308+0.112734376881499D-5*t317
     &-0.4815887105457136D-3*t321-0.1028646684980902D0*t327
     &+0.9631774210914273D-3*t332+0.5143233424904509D-1*t351
     &+0.1088426657654099D-2*t356+0.8637605650141684D-5*t370
     &+0.5095763183186173D-5*t377-0.8248345450739115D-5*t382
     &+0.3085940054942705D0*t389
      t529 = 0.2051285184758976D-2*t308+0.1516365042658901D-5*t317
     &-0.6477742688712557D-3*t321-0.1383609789223028D0*t327
     &+0.1295548537742511D-2*t332+0.6918048946115138D-1*t351
     &+0.146401850156108D-2*t356+0.1161825134662837D-4*t370
     &+0.6854197779240001D-5*t377-0.1109466610563933D-4*t382
     &+0.4150829367669083D0*t389-0.2767219578446055D0*t392
     &-0.146401850156108D-2*t396+0.4612032630743426D-1*t399
     &-0.1537344210247809D-1*t406-0.4880061671870268D-4*t412
     &+0.2159247562904186D-3*t416+0.2928037003122161D-2*t421
     &+0.9224065261486851D-1*t426-0.4612032630743426D-1*t429
      t530 = t172*t529
      t533 = 1/t83/t82
      t534 = t178**2
      t535 = t533*t534
      t537 = -0.2057293369961804D0*t392-0.1088426657654099D-2*t396
     &+0.3428822283269673D-1*t399-0.1142940761089891D-1*t406
     &-0.3628088858846997D-4*t412+0.1605295701819045D-3*t416
     &+0.2176853315308198D-2*t421+0.6857644566539345D-1*t426
     &-0.3428822283269673D-1*t429+0.3278192763011299D1*t530
     &-0.2369304401100098D2*t535
      t541 = t186*t190
      t548 = t86/t189/t89
      t549 = t198**2
      t567 = 0.228754637509214D-2*t308+0.1691015653222485D-5*t317
     &-0.7223830658185704D-3*t321-0.1542970027471353D0*t327
     &+0.1444766131637141D-2*t332+0.7714850137356763D-1*t351
     &+0.1632639986481149D-2*t356+0.1295640847521253D-4*t370
     &+0.764364477477926D-5*t377-0.1237251817610867D-4*t382
     &+0.4628910082414058D0*t389
      t579 = -0.3085940054942705D0*t392-0.1632639986481149D-2*t396
     &+0.5143233424904509D-1*t399-0.1714411141634836D-1*t406
     &-0.5442133288270495D-4*t412+0.2407943552728568D-3*t416
     &+0.3265279972962297D-2*t421+0.1028646684980902D0*t426
     &-0.5143233424904509D-1*t429+0.4917289144516948D1*t530
     &-0.3553956601650148D2*t535
      t584 = t309+t318-t328+t371+t357-t383+t333+t352-t322+t400-t397
     &-t393+t390
      t586 = t117*t78
      t606 = t54*t373
      t610 = t378-t430+t427+t422+t417-t413-t407-0.1D1*t449
     &-0.2721066644135248D-3*t586*t160*t206+0.3812577291820233D-3*t218
     &*t219*t305*t80+0.2818359422037475D-6*t218*t54*t310*t314*t80
     &-0.8572055708174181D-2*t106*t214*t202+0.2721066644135248D-3*t118
     &*t93*t353*t121+0.1273940795796543D-5*t586*t606*t375
      t612 = t105*t329
      t613 = t54*t108
      t634 = 1/t71
      t638 = t66*t152
      t650 = t59**2
      t652 = 1/t650/t58/t57
      t654 = t652/t63
      t659 = t136*t80
      t663 = 0.1714411141634836D-1*t612*t613*t210
     &+0.1285808356226127D-1*t92*t77*t431*t225+0.4681769917478079D-2
     &*t162*t32*t135-0.1D1*t151*t162*t165+0.1285808356226127D-1*t92
     &*t77*t434*t225+0.3333333333333333D0*t5*t162-0.2D1*t437+0.2D1
     &*t634*t162*t138-0.1203971776364284D-3*t638*t221
     &+0.2159401412535421D-5*t363*t93*t364*t368-0.2857351902724727D-2
     &*t404*t93*t364*t80+0.105811516582519D-1*t654*t162*t108*t61
     &-0.1203971776364284D-3*t218*t160*t659
      t667 = t117*t329
      t671 = t105*t78
      t675 = t54*t138
      t680 = t93*t364*t121
      t683 = t32*t80
      t715 = -0.2721066644135248D-3*t118*t214*t206
     &+0.5442133288270495D-3*t667*t613*t419-0.8572055708174181D-2*t671
     &*t160*t202+0.2407943552728568D-3*t330*t675*t659
     &-0.9070222147117492D-5*t411*t680+0.4013239254547614D-4*t671*t606
     &*t683-0.2571616712452254D-1*t127*t93*t325*t80
     &+0.7714850137356763D-1*t386*t93*t387*t80-0.5143233424904509D-1
     &*t330*t160*t210+0.1285808356226127D-1*t92*t348*t54*t80
     &+0.2571616712452254D-1*t92*t151*t159*t225-0.5143233424904509D-1
     &*t127*t214*t210-0.2062086362684779D-5*t380*t680
     &+0.8572055708174181D-2*t106*t93*t353*t80
      t719 = 1.D0*t431*t54+1.D0*t434*t54-0.155544175017655D3*t437*t74
     &*t441+0.7777208750882749D2*t431*t74*t64*t169
     &-0.7777208750882749D2*t449*t74*t64*t169-0.2592402916960916D2
     &*t162*t148*t64*t166*t108*t168-0.8229182527097421D0*t162/t57/t5
     &*t63*t166*t108*t465-0.155544175017655D3*t164*t470*t138*t168
     &+0.7777208750882749D2*t164*t477*t151*t168-0.3641110197183013D0
     &*t482*t165*t32*t135*t168+(t499+t537)*t90*t97-2.D0*t541*t199+2.D0
     &*t187*t229+2.D0*t548*t97*t549-2.D0*t191*t229*t198-1.D0*t191*t97*
     &(t567+t579)+t91*(t584+t610+t663+t715)
      v2rhoa2(i) = v2rhoa2(i)+0.1554534543482745D-1*rhob*t719
      v2rhoab(i) = v2rhoab(i)+0.1554534543482745D-1*t160
     &-0.1208993965512352D1*t170+0.1554534543482745D-1*t188
     &-0.1554534543482745D-1*t200+0.1554534543482745D-1*t230
      v2rhob2(i) = v2rhob2(i)
      t731 = t423*t107*t108*t241
      t732 = 0.8572055708174181D-2*t731
      t733 = t130*t65
      t735 = t733*t73*t80
      t736 = 0.4073997761173315D-3*t735
      t739 = t106*t107*t305*sigmaaa
      t740 = 0.7524823602276776D-5*t739
      t741 = t241*t121
      t743 = t418*t119*t741
      t744 = 0.2721066644135248D-3*t743
      t745 = t652*t116
      t748 = t745*t73*t145*t121
      t749 = 0.1293227644990607D-4*t748
      t751 = sigmaaa*t121
      t753 = t118*t77*t305*t751
      t754 = 0.2388638992118519D-6*t753
      t757 = t386*t107*t138*t241
      t758 = 0.7714850137356763D-1*t757
      t759 = t65*t126
      t761 = t759*t246*t138
      t762 = 0.2444398656703989D-2*t761
      t764 = t330*t210*t239
      t765 = 0.4514894161366065D-4*t764
      t767 = sigmaaa*t135
      t769 = -0.1849622637793975D0*t130-0.4681769917478079D-2*t767
      t771 = t127*t107*t769
      t772 = 0.2571616712452254D-1*t771
      t774 = t127*t319*t241
      t775 = 0.2571616712452254D-1*t774
      t780 = -0.2534746483141907D0*t130*t74+0.3168433103927383D-1*t145
     &*t148
      t783 = t66*t73*t780*t80
      t784 = 0.1285808356226127D-1*t783
      t786 = t92*t319*t239
      t787 = 0.2257447080683033D-4*t786
      t788 = t135*t241
      t790 = t330*t683*t788
      t791 = 0.1203971776364284D-3*t790
      t793 = 1/t311/rhoa
      t796 = t245*t793*t80*t32
      t797 = 0.5722056048640281D-5*t796
      t798 = t32*sigmaaa
      t801 = 1/t25/t311/t128
      t804 = t92*t107*t798*t801
      t805 = 0.1056884783264053D-6*t804
      t807 = t92*t107*t767
      t808 = 0.1203971776364284D-3*t807
      t809 = t732-t736-t740+t744-t749-t754+t758-t762-t765-t772-t775
     &+t784+t787+t791-t797-t805-t808
      t812 = t159*t252
      t815 = t159*t255
      t819 = t252*t162
      t828 = t162*t255
      t833 = t470*t241*t168
      t839 = t477*t27*t168
      t844 = t165*sigmaaa*t37*t168
      t848 = t533*t178*t263
      t867 = 0.4612032630743426D-1*t731-0.2191937529546115D-2*t735
     &-0.4048589180445348D-4*t739+0.146401850156108D-2*t743
     &-0.6957967027662436D-4*t748-0.12851620836075D-5*t753
     &+0.4150829367669083D0*t757-0.1315162517727669D-1*t761
     &-0.2429153508267209D-3*t764-0.1383609789223028D0*t771
     &-0.1383609789223028D0*t774+0.6918048946115138D-1*t783
     &+0.1214576754133604D-3*t786+0.6477742688712557D-3*t790
     &-0.3078644156046066D-4*t796-0.5686368909970879D-6*t804
     &-0.6477742688712557D-3*t807
      t868 = t172*t867
      t887 = -0.2369304401100098D2*t848+0.3278192763011299D1*t868
     &+0.3428822283269673D-1*t731-0.1629599104469326D-2*t735
     &-0.300992944091071D-4*t739+0.1088426657654099D-2*t743
     &-0.5172910579962427D-4*t748-0.9554555968474075D-6*t753
     &+0.3085940054942705D0*t757-0.9777594626815955D-2*t761
     &-0.1805957664546426D-3*t764-0.1028646684980902D0*t771
     &-0.1028646684980902D0*t774+0.5143233424904509D-1*t783
     &+0.9029788322732131D-4*t786+0.4815887105457136D-3*t790
     &-0.2288822419456112D-4*t796-0.4227539133056213D-6*t804
     &-0.4815887105457136D-3*t807
      t893 = t269*t190
      t921 = -0.3553956601650148D2*t848+0.4917289144516948D1*t868
     &+0.5143233424904509D-1*t731-0.2444398656703989D-2*t735
     &-0.4514894161366065D-4*t739+0.1632639986481149D-2*t743
     &-0.775936586994364D-4*t748-0.1433183395271111D-5*t753
     &+0.4628910082414058D0*t757-0.1466639194022393D-1*t761
     &-0.2708936496819639D-3*t764-0.1542970027471353D0*t771
     &-0.1542970027471353D0*t774+0.7714850137356763D-1*t783
     &+0.135446824840982D-3*t786+0.7223830658185704D-3*t790
     &-0.3433233629184168D-4*t796-0.634130869958432D-6*t804
     &-0.7223830658185704D-3*t807
      t941 = t54*t305
      t959 = 0.7714850137356763D-1*t66*t385*t77*t675*t280-t736
     &-0.1360533322067624D-3*t586*t253*t206+0.1203971776364284D-3*t330
     &*t219*t788*t80-0.2388638992118519D-6*t586*t941*t751
     &-0.1293227644990607D-4*t745*t73*t145*t54*t121-t749
     &-0.2571616712452254D-1*t127*t214*t280-0.2571616712452254D-1*t330
     &*t253*t210-t762-t775+t784+0.4681769917478079D-2*t255*t32*t135
      t967 = t759*t27
      t974 = t239*t80
      t978 = 0.8572055708174181D-2*t612*t613*t280-t772
     &-0.7524823602276776D-5*t671*t941*sigmaaa*t80-t765
     &-0.2444398656703989D-2*t967*t675*t80-0.4286027854087091D-2*t671
     &*t253*t202+t758-t754+t744-t740-0.4514894161366065D-4*t330*t675
     &*t974+t732-t797
      t1011 = -0.1203971776364284D-3*t218*t288*t220
     &+0.2721066644135248D-3*t667*t613*t741-t805+0.1285808356226127D-1
     &*t92*t151*t252*t225+t791+t787-0.2571616712452254D-1*t330*t160
     &*t280-t808+0.2257447080683033D-4*t638*t290-0.4073997761173315D-3
     &*t733*t73*t54*t80-0.5722056048640281D-5*t245*t793*t219*t80
     &-0.601985888182142D-4*t218*t253*t659-0.1056884783264053D-6*t218
     &*t54*t798*t801*t80
      t1026 = t245*t27
      t1035 = t634*t255
      t1049 = 0.1285808356226127D-1*t92*t77*t809*t225
     &+0.1285808356226127D-1*t218*t812*t225-0.2571616712452254D-1*t127
     &*t93*t769*t80+0.2257447080683033D-4*t218*t160*t974
     &+0.1222199328351994D-2*t1026*t160*t80-0.1D1*t828
     &+0.1285808356226127D-1*t92*t780*t54*t80+0.2D1*t1035*t138-0.1D1
     &*t151*t255*t165+0.3333333333333333D0*t5*t255-0.1D1*t815-0.1D1
     &*t819+0.105811516582519D-1*t654*t255*t108*t61
      t1053 = 1.D0*t809*t54+1.D0*t812*t54-0.7777208750882749D2*t815
     &*t74*t441-0.7777208750882749D2*t819*t74*t441
     &+0.7777208750882749D2*t809*t74*t64*t169-0.7777208750882749D2
     &*t828*t74*t441-0.155544175017655D3*t164*t833
     &+0.7392469698735191D1*t162*t104*t64*t839+0.136541632394363D0
     &*t482*t844+t887*t90*t97-1.D0*t541*t277+t187*t297-1.D0*t893*t199
     &+2.D0*t548*t199*t276-1.D0*t191*t297*t198-1.D0*t191*t97*t921+t270
     &*t229-1.D0*t191*t229*t276+t91*(t959+t978+t1011+t1049)
      v2rhoasigmaaa(i) = v2rhoasigmaaa(i)+0.1554534543482745D-1*rhob
     &*t1053
      v2rhobsigmaaa(i) = v2rhobsigmaaa(i)+0.1554534543482745D-1*t253
     &-0.1208993965512352D1*t258+0.1554534543482745D-1*t271
     &-0.1554534543482745D-1*t278+0.1554534543482745D-1*t298
      t1063 = t241**2
      t1065 = t386*t107*t1063
      t1066 = 0.7714850137356763D-1*t1065
      t1068 = t759*t246*t241
      t1069 = 0.4888797313407978D-2*t1068
      t1071 = t330*t280*t239
      t1072 = 0.9029788322732131D-4*t1071
      t1073 = t107*t37
      t1074 = t127*t1073
      t1075 = 0.2257447080683033D-4*t1074
      t1076 = 1/t311
      t1079 = t245*t1076*t80*sigmaaa
      t1080 = 0.4291542036480211D-5*t1079
      t1083 = 1/t25/t311/t23
      t1086 = t92*t107*t32*t1083
      t1087 = 0.39633179372402D-7*t1086
      t1088 = t92*t1073
      t1089 = 0.2257447080683033D-4*t1088
      t1090 = t1066-t1069-t1072-t1075+t1080+t1087+t1089
      t1093 = t252**2
      t1096 = t252*t255
      t1105 = t255**2
      t1119 = t263**2
      t1120 = t533*t1119
      t1130 = t172*(0.4150829367669083D0*t1065-0.2630325035455338D-1
     &*t1068-0.4858307016534418D-3*t1071-0.1214576754133604D-3*t1074
     &+0.2308983117034549D-4*t1079+0.2132388341239079D-6*t1086
     &+0.1214576754133604D-3*t1088)
      t1146 = t276**2
      t1176 = t93*t289
      t1196 = t54*t241
      t1203 = -0.175566371905428D-2*t255*sigmaaa*t37
     &+0.4514894161366065D-4*t218*t253*t974+0.1285808356226127D-1*t92
     &*t77*t1090*t225+0.2257447080683033D-4*t92*t1176
     &+0.39633179372402D-7*t218*t219*t1083*t80-0.2D1*t1096-0.1D1*t1105
     &+0.4291542036480211D-5*t245*t1076*t288*t80+0.2444398656703989D-2
     &*t1026*t253*t80-0.9505299311782149D-1*t27*t255*t74*t165
     &-0.9029788322732131D-4*t330*t1196*t974-0.5143233424904509D-1
     &*t330*t253*t280
      t1219 = -0.2257447080683033D-4*t127*t1176-0.4888797313407978D-2
     &*t967*t1196*t80+0.2D1*t1035*t241+0.7714850137356763D-1*t386*t93
     &*t1063*t80+t1066-t1069-t1072-t1075+t1080+t1087+t1089
     &+0.1285808356226127D-1*t92*t77*t1093*t225
      s1 = 1.D0*t1090*t54+1.D0*t1093*t54-0.155544175017655D3*t1096*t74
     &*t441+0.7777208750882749D2*t1090*t74*t64*t169
     &-0.7777208750882749D2*t1105*t74*t64*t169-0.155544175017655D3
     &*t257*t833+0.7392469698735191D1*t255*t104*t64*t839
      t1222 = s1+0.136541632394363D0*t256*t439*t844+(
     &-0.2369304401100098D2*t1120+0.3278192763011299D1*t1130
     &+0.3085940054942705D0*t1065-0.1955518925363191D-1*t1068
     &-0.3611915329092852D-3*t1071-0.9029788322732131D-4*t1074
     &+0.1716616814592084D-4*t1079+0.158532717489608D-6*t1086
     &+0.9029788322732131D-4*t1088)*t90*t97-2.D0*t893*t277+2.D0*t270
     &*t297+2.D0*t548*t97*t1146-2.D0*t191*t297*t276-1.D0*t191*t97*(
     &-0.3553956601650148D2*t1120+0.4917289144516948D1*t1130
     &+0.4628910082414058D0*t1065-0.2933278388044787D-1*t1068
     &-0.5417872993639278D-3*t1071-0.135446824840982D-3*t1074
     &+0.2574925221888126D-4*t1079+0.237799076234412D-6*t1086
     &+0.135446824840982D-3*t1088)+t91*(t1203+t1219)
      v2sigmaaa2(i) = v2sigmaaa2(i)+0.1554534543482745D-1*rhob*t1222
      endif
      endif ! rhoa,rhob
      endif ! rho
      enddo
      
      endif ! ideriv
      return
      end
      
      
      subroutine rks_c_ft97
     & (ideriv,npt,rhoa1,sigmaaa1,
     &  zk,vrhoa,vsigmaaa,
     &  v2rhoa2,v2rhoasigmaaa,v2sigmaaa2)
c
c     M. Filatov, and W. Thiel
c     A nonlocal correlation energy density functional from a 
c     Coulomb hole model
c     Int. J. Quant. Chem. 62 (1997) 603-616
c
c     M. Filatov, and W. Thiel
c     A new gradient-corrected exchange-correlation density functional
c     Mol. Phys. 91 (1997) 847-859
c
c
c     CITATION:
c
c     Functionals were obtained from the Density Functional Repository 
c     as developed and distributed by the Quantum Chemistry Group, 
c     CCLRC Daresbury Laboratory, Daresbury, Cheshire, WA4 4AD 
c     United Kingdom. Contact Huub van Dam (h.j.j.vandam@dl.ac.uk) or 
c     Paul Sherwood for further information.
c
c     COPYRIGHT:
c
c     Users may incorporate the source code into software packages and
c     redistribute the source code provided the source code is not
c     changed in anyway and is properly cited in any documentation or
c     publication related to its use.
c
c     ACKNOWLEDGEMENT:
c
c     The source code was generated using Maple 8 through a modified
c     version of the dfauto script published in:
c
c        R. Strange, F.R. Manby, P.J. Knowles
c        Automatic code generation in density functional theory
c        Comp. Phys. Comm. 136 (2001) 310-318.
c
      implicit real*8 (a-h,o-z)
      integer ideriv,npt
      real*8 rhoa1(npt)
      real*8 sigmaaa1(npt)
      real*8 zk(npt),vrhoa(npt),vsigmaaa(npt)
      real*8 v2rhoa2(npt),v2rhoasigmaaa(npt),v2sigmaaa2(npt)
      parameter(tolmin=1.0d-20)
      parameter(tolmax=1.0d+5)
      
      if(ideriv.eq.0) then
      
      do i=1,npt
      rho = dmax1(0.D0,rhoa1(i))
      zk(i) = 0.D0
      if(rho.gt.tolmin) then
      sigma = dmax1(0.D0,sigmaaa1(i))
      t1 = 1/rho
      t2 = t1**(1.D0/3.D0)
      trsa = 0.1620017014140023D-1*t2
      t3 = t1**(1.D0/15.D0)
      t4 = t3**2
      t6 = dexp(-0.5940962989070748D0*t4)
      t8 = t1**(1.D0/6.D0)
      t10 = dexp(-0.9630597286858919D0*t8)
      t13 = (0.1247511874D1+0.812904345D0*t6-0.859614445D0*t10)**2
      t14 = sigma**2
      t15 = rho**2
      t16 = t15**2
      t18 = rho**(1.D0/3.D0)
      t21 = t14/t18/t16/rho
      t24 = (0.1D1+0.1127183082292567D0*t21)**2
      t26 = t18**2
      t36 = dexp(-0.1127183082292567D0*t21)
      t37 = t36**2
      trcaa = t13*t24/(0.1D1+0.312690210142125D0*sigma/t26/t15/t2)*t37
      if(trsa.lt.tolmax*trcaa) then
      t2 = 1/rho
      t3 = t2**(1.D0/3.D0)
      t4 = t3**2
      t5 = t2**(1.D0/6.D0)
      t9 = (0.8301627136974294D0*t5+0.1354632918982911D1*t3)**2
      t13 = dexp(-0.6108870577108572D0*t4/t9)
      t15 = 0.3141592653589793D1**2
      t17 = 2**(1.D0/3.D0)
      t20 = 1/0.3141592653589793D1*t2
      t21 = t20**(1.D0/3.D0)
      t23 = 2**(1.D0/15.D0)
      t24 = t23**2
      t25 = t20**(1.D0/15.D0)
      t26 = t25**2
      t29 = dexp(-0.630966299458536D0*t24*t26)
      t31 = 2**(1.D0/6.D0)
      t32 = t20**(1.D0/6.D0)
      t35 = dexp(-0.1038340679664977D1*t31*t32)
      t38 = (0.1247511874D1+0.812904345D0*t29-0.859614445D0*t35)**2
      t40 = 0.3141592653589793D1**(1.D0/3.D0)
      t43 = sigma**2
      t45 = rho**2
      t46 = t45**2
      t48 = rho**(1.D0/3.D0)
      t50 = 1/t48/t46/rho
      t52 = 1/t40/0.3141592653589793D1*t43*t17*t50
      t55 = (0.1D1+0.4116400544093991D0*t52)**2
      t58 = t40**2
      t61 = t48**2
      t63 = 1/t61/t45
      t71 = dexp(0.8232801088187983D0*t52)
      t76 = expei(-0.1858628590577086D0/t15*t17*t21/t38/t55*(0.1D1
     &+0.3634859066227017D0/t58*sigma*t17*t63/t21)*t71)
      t78 = t2**(1.D0/15.D0)
      t79 = t78**2
      t81 = dexp(-0.5940962989070748D0*t79)
      t84 = dexp(-0.9630597286858919D0*t5)
      t87 = (0.1247511874D1+0.812904345D0*t81-0.859614445D0*t84)**2
      t89 = t3/t87
      t90 = t43*t50
      t93 = (0.1D1+0.1127183082292567D0*t90)**2
      t94 = 1/t93
      t99 = 0.1D1+0.312690210142125D0*sigma*t63/t3
      t102 = dexp(0.2254366164585135D0*t90)
      t104 = t89*t94*t99*t102
      t105 = dsqrt(t104)
      t123 = 0.7772672717413724D-2*rho*t13*(1.D0*t76+(0.6D1
     &+0.5091195559614694D0*t105+0.6480068056560093D-1*t104)/(0.3D1
     &+0.763679333942204D0*t105+0.972010208484014D-1*t104)*
     &(0.1620017014140023D-1*t89*t94*t99*t76*t102+0.1D1))
      zk(i) = zk(i)+2*t123
      endif
      t1 = 1/rho
      t2 = t1**(1.D0/3.D0)
      trsa = 0.1620017014140023D-1*t2
      t3 = t1**(1.D0/15.D0)
      t4 = t3**2
      t5 = t4**2
      t7 = dexp(-0.6837598574312962D-1*t5)
      t10 = (0.942486901D0+0.349064173D0*t7)**2
      t11 = rho**2
      t12 = rho**(1.D0/3.D0)
      t13 = t12**2
      t16 = sigma/t13/t11
      t18 = sigma**2
      t19 = t11**2
      t23 = t18/t12/t19/rho
      t26 = (0.1D1+0.1101034845366882D0*t16+0.1105998838086603D-2*t23)
     &**2
      t34 = dexp(-0.1105998838086603D-2*t23)
      t35 = t34**2
      trcba = t10*t26/(0.1D1+0.1197592668846558D0*t16/t2)*t35
      if(trsa.lt.tolmax*trcba) then
      t2 = 0.3141592653589793D1**2
      t4 = 2**(1.D0/3.D0)
      t7 = 1/rho
      t8 = 1/0.3141592653589793D1*t7
      t9 = t8**(1.D0/3.D0)
      t11 = 2**(1.D0/15.D0)
      t12 = t11**2
      t13 = t12**2
      t14 = t8**(1.D0/15.D0)
      t15 = t14**2
      t16 = t15**2
      t19 = dexp(-0.7712625328179681D-1*t13*t16)
      t22 = (0.942486901D0+0.349064173D0*t19)**2
      t24 = 0.3141592653589793D1**(1.D0/3.D0)
      t25 = t24**2
      t27 = 1/t25*sigma
      t28 = t4**2
      t29 = rho**2
      t30 = rho**(1.D0/3.D0)
      t31 = t30**2
      t33 = 1/t31/t29
      t39 = sigma**2
      t41 = t29**2
      t44 = 1/t30/t41/rho
      t46 = 1/t24/0.3141592653589793D1*t39*t4*t44
      t49 = (0.1D1+0.1487810599361293D0*t27*t28*t33
     &+0.4039037038781002D-2*t46)**2
      t59 = dexp(0.8078074077562004D-2*t46)
      t64 = expei(-0.1858628590577086D0/t2*t4*t9/t22/t49*(0.1D1
     &+0.1392138426088027D0*t27*t4*t33/t9)*t59)
      t66 = t7**(1.D0/3.D0)
      t67 = t7**(1.D0/15.D0)
      t68 = t67**2
      t69 = t68**2
      t71 = dexp(-0.6837598574312962D-1*t69)
      t74 = (0.942486901D0+0.349064173D0*t71)**2
      t76 = t66/t74
      t77 = sigma*t33
      t79 = t39*t44
      t82 = (0.1D1+0.1101034845366882D0*t77+0.1105998838086603D-2*t79)
     &**2
      t83 = 1/t82
      t87 = 0.1D1+0.1197592668846558D0*t77/t66
      t90 = dexp(0.2211997676173206D-2*t79)
      t92 = t76*t83*t87*t90
      t93 = dsqrt(t92)
      t111 = 0.7772672717413724D-2*rho*(1.D0*t64+(0.6D1
     &+0.5091195559614694D0*t93+0.6480068056560093D-1*t92)/(0.3D1
     &+0.763679333942204D0*t93+0.972010208484014D-1*t92)*
     &(0.1620017014140023D-1*t76*t83*t87*t64*t90+0.1D1))
      zk(i) = zk(i)+2*t111
      endif
      endif ! rho
      enddo
      
      else if(ideriv.eq.1) then
      
      do i=1,npt
      rho = dmax1(0.D0,rhoa1(i))
      zk(i) = 0.D0
      vrhoa(i) = 0.D0
      vsigmaaa(i) = 0.D0
      if(rho.gt.tolmin) then
      sigma = dmax1(0.D0,sigmaaa1(i))
      t1 = 1/rho
      t2 = t1**(1.D0/3.D0)
      trsa = 0.1620017014140023D-1*t2
      t3 = t1**(1.D0/15.D0)
      t4 = t3**2
      t6 = dexp(-0.5940962989070748D0*t4)
      t8 = t1**(1.D0/6.D0)
      t10 = dexp(-0.9630597286858919D0*t8)
      t13 = (0.1247511874D1+0.812904345D0*t6-0.859614445D0*t10)**2
      t14 = sigma**2
      t15 = rho**2
      t16 = t15**2
      t18 = rho**(1.D0/3.D0)
      t21 = t14/t18/t16/rho
      t24 = (0.1D1+0.1127183082292567D0*t21)**2
      t26 = t18**2
      t36 = dexp(-0.1127183082292567D0*t21)
      t37 = t36**2
      trcaa = t13*t24/(0.1D1+0.312690210142125D0*sigma/t26/t15/t2)*t37
      if(trsa.lt.tolmax*trcaa) then
      t2 = 1/rho
      t3 = t2**(1.D0/3.D0)
      t4 = t3**2
      t5 = t2**(1.D0/6.D0)
      t8 = 0.8301627136974294D0*t5+0.1354632918982911D1*t3
      t9 = t8**2
      t10 = 1/t9
      t13 = dexp(-0.6108870577108572D0*t4*t10)
      t14 = rho*t13
      t15 = 0.3141592653589793D1**2
      t17 = 2**(1.D0/3.D0)
      t20 = 1/0.3141592653589793D1*t2
      t21 = t20**(1.D0/3.D0)
      t23 = 2**(1.D0/15.D0)
      t24 = t23**2
      t25 = t20**(1.D0/15.D0)
      t26 = t25**2
      t29 = dexp(-0.630966299458536D0*t24*t26)
      t31 = 2**(1.D0/6.D0)
      t32 = t20**(1.D0/6.D0)
      t35 = dexp(-0.1038340679664977D1*t31*t32)
      t38 = (0.1247511874D1+0.812904345D0*t29-0.859614445D0*t35)**2
      t40 = 0.3141592653589793D1**(1.D0/3.D0)
      t43 = sigma**2
      t45 = rho**2
      t46 = t45**2
      t48 = rho**(1.D0/3.D0)
      t50 = 1/t48/t46/rho
      t52 = 1/t40/0.3141592653589793D1*t43*t17*t50
      t55 = (0.1D1+0.4116400544093991D0*t52)**2
      t58 = t40**2
      t61 = t48**2
      t63 = 1/t61/t45
      t71 = dexp(0.8232801088187983D0*t52)
      t76 = expei(-0.1858628590577086D0/t15*t17*t21/t38/t55*(0.1D1
     &+0.3634859066227017D0/t58*sigma*t17*t63/t21)*t71)
      t78 = t2**(1.D0/15.D0)
      t79 = t78**2
      t81 = dexp(-0.5940962989070748D0*t79)
      t84 = dexp(-0.9630597286858919D0*t5)
      t86 = 0.1247511874D1+0.812904345D0*t81-0.859614445D0*t84
      t87 = t86**2
      t88 = 1/t87
      t89 = t3*t88
      t90 = t43*t50
      t92 = 0.1D1+0.1127183082292567D0*t90
      t93 = t92**2
      t94 = 1/t93
      t96 = 1/t3
      t99 = 0.1D1+0.312690210142125D0*sigma*t63*t96
      t100 = t94*t99
      t102 = dexp(0.2254366164585135D0*t90)
      t104 = t89*t100*t102
      t105 = dsqrt(t104)
      t108 = 0.6D1+0.5091195559614694D0*t105+0.6480068056560093D-1*t104
      t111 = 0.3D1+0.763679333942204D0*t105+0.972010208484014D-1*t104
      t112 = 1/t111
      t113 = t108*t112
      t114 = t89*t94
      t115 = t99*t76
      t119 = 0.1620017014140023D-1*t114*t115*t102+0.1D1
      t121 = 1.D0*t76+t113*t119
      t123 = 0.7772672717413724D-2*t14*t121
      zk(i) = zk(i)+2*t123
      t125 = t13*t121
      t128 = 1/t45
      t134 = t5**2
      t135 = t134**2
      t138 = 1/t135/t5*t128
      t140 = 1/t4
      t151 = t140*t88*t94
      t152 = t99*t102
      t154 = t151*t152*t128
      t155 = 0.1080011342760016D-1*t154
      t159 = t3/t87/t86*t94
      t160 = t79**2
      t162 = t160**2
      t170 = 0.1287849233946613D0/t162/t160/t78*t128*t81
     &-0.2759533513920579D0*t138*t84
      t172 = t159*t152*t170
      t173 = 0.3240034028280047D-1*t172
      t175 = 1/t93/t92
      t176 = t89*t175
      t179 = 1/t48/t46/t45
      t181 = t152*t43*t179
      t182 = t176*t181
      t183 = 0.389558564557814D-1*t182
      t197 = -0.1667681120758D1*sigma/t61/t45/rho*t96
     &+0.20846014009475D0*sigma/t61/t46/t3/t2
      t200 = t89*t94*t197*t102
      t201 = 0.1620017014140023D-1*t200
      t202 = t114*t181
      t203 = 0.389558564557814D-1*t202
      t204 = -t155-t173+t183+t201-t203
      t213 = dexp(-0.2254366164585135D0*t90)
      t214 = t93/t99*t213
      t217 = 1/t105
      t224 = t217*(-0.5810796994275671D-1*t154-0.1743239098282701D0
     &*t172+0.2095946261306892D0*t182+0.8716195491413506D-1*t200
     &-0.2095946261306892D0*t202)
      t234 = t111**2
      t236 = t108/t234
      t256 = t89*t175*t99
      t259 = t76*t43*t179*t102
      t266 = t89*t100
      t270 = t76*t102
      t274 = -0.1080011342760016D-1*t151*t115*t128*t102
     &-0.3240034028280047D-1*t159*t115*t170*t102+0.389558564557814D-1
     &*t256*t259+0.1620017014140023D-1*t114*t197*t76*t102
     &-0.389558564557814D-1*t266*t259+0.1620017014140023D-1*t114*t99
     &*t204*t270-t155-t173+t183+t201-t203
      vrhoa(i) = vrhoa(i)+0.1554534543482745D-1*t125
     &+0.7772672717413724D-2*rho*(0.8145160769478096D0*t96*t10*t128
     &+0.1221774115421714D1*t4/t9/t8*(-0.2767209045658098D0*t138
     &-0.9030886126552743D0*t140*t128))*t125+0.7772672717413724D-2*t14
     &*(1.D0*t204*t76+0.6172774676263781D2*t204*t96*t87*t214+
     &(0.2920537730383703D1*t224-0.4320045371040062D-1*t154
     &-0.1296013611312019D0*t172+0.1558234258231256D0*t182
     &+0.6480068056560093D-1*t200-0.1558234258231256D0*t202)*t112*t119
     &-1.D0*t236*t119*(0.4380806595575555D1*t224-0.6480068056560093D-1
     &*t154-0.1944020416968028D0*t172+0.2337351387346884D0*t182
     &+0.972010208484014D-1*t200-0.2337351387346884D0*t202)+t113*t274)
      t281 = t152*sigma*t50
      t282 = t176*t281
      t283 = 0.2921689234183605D-1*t282
      t284 = t88*t94
      t286 = t284*t63*t102
      t287 = 0.2026253842341047D-1*t286
      t288 = t114*t281
      t289 = 0.2921689234183605D-1*t288
      t290 = -t283+t287+t289
      t302 = t217*(-0.1571959695980169D0*t282+0.1090187599939973D0
     &*t286+0.1571959695980169D0*t288)
      t320 = t76*sigma*t50*t102
      vsigmaaa(i) = vsigmaaa(i)+0.1554534543482745D-1*t14*(1.D0*t290
     &*t76+0.6172774676263781D2*t290*t96*t87*t214+
     &(0.2920537730383703D1*t302-0.1168675693673442D0*t282
     &+0.8105015369364188D-1*t286+0.1168675693673442D0*t288)*t112*t119
     &-1.D0*t236*t119*(0.4380806595575555D1*t302-0.1753013540510163D0
     &*t282+0.1215752305404628D0*t286+0.1753013540510163D0*t288)+t113*
     &(-0.2921689234183605D-1*t256*t320+0.2026253842341047D-1*t284*t63
     &*t76*t102+0.2921689234183605D-1*t266*t320+0.1620017014140023D-1
     &*t114*t99*t290*t270-t283+t287+t289))
      endif
      t1 = 1/rho
      t2 = t1**(1.D0/3.D0)
      trsa = 0.1620017014140023D-1*t2
      t3 = t1**(1.D0/15.D0)
      t4 = t3**2
      t5 = t4**2
      t7 = dexp(-0.6837598574312962D-1*t5)
      t10 = (0.942486901D0+0.349064173D0*t7)**2
      t11 = rho**2
      t12 = rho**(1.D0/3.D0)
      t13 = t12**2
      t16 = sigma/t13/t11
      t18 = sigma**2
      t19 = t11**2
      t23 = t18/t12/t19/rho
      t26 = (0.1D1+0.1101034845366882D0*t16+0.1105998838086603D-2*t23)
     &**2
      t34 = dexp(-0.1105998838086603D-2*t23)
      t35 = t34**2
      trcba = t10*t26/(0.1D1+0.1197592668846558D0*t16/t2)*t35
      if(trsa.lt.tolmax*trcba) then
      t2 = 0.3141592653589793D1**2
      t4 = 2**(1.D0/3.D0)
      t7 = 1/rho
      t8 = 1/0.3141592653589793D1*t7
      t9 = t8**(1.D0/3.D0)
      t11 = 2**(1.D0/15.D0)
      t12 = t11**2
      t13 = t12**2
      t14 = t8**(1.D0/15.D0)
      t15 = t14**2
      t16 = t15**2
      t19 = dexp(-0.7712625328179681D-1*t13*t16)
      t22 = (0.942486901D0+0.349064173D0*t19)**2
      t24 = 0.3141592653589793D1**(1.D0/3.D0)
      t25 = t24**2
      t27 = 1/t25*sigma
      t28 = t4**2
      t29 = rho**2
      t30 = rho**(1.D0/3.D0)
      t31 = t30**2
      t33 = 1/t31/t29
      t39 = sigma**2
      t41 = t29**2
      t44 = 1/t30/t41/rho
      t46 = 1/t24/0.3141592653589793D1*t39*t4*t44
      t49 = (0.1D1+0.1487810599361293D0*t27*t28*t33
     &+0.4039037038781002D-2*t46)**2
      t59 = dexp(0.8078074077562004D-2*t46)
      t64 = expei(-0.1858628590577086D0/t2*t4*t9/t22/t49*(0.1D1
     &+0.1392138426088027D0*t27*t4*t33/t9)*t59)
      t66 = t7**(1.D0/3.D0)
      t67 = t7**(1.D0/15.D0)
      t68 = t67**2
      t69 = t68**2
      t70 = 0.6837598574312962D-1*t69
      t71 = dexp(-t70)
      t73 = 0.942486901D0+0.349064173D0*t71
      t74 = t73**2
      t75 = 1/t74
      t76 = t66*t75
      t77 = sigma*t33
      t79 = t39*t44
      t81 = 0.1D1+0.1101034845366882D0*t77+0.1105998838086603D-2*t79
      t82 = t81**2
      t83 = 1/t82
      t84 = 1/t66
      t87 = 0.1D1+0.1197592668846558D0*t77*t84
      t88 = t83*t87
      t89 = 0.2211997676173206D-2*t79
      t90 = dexp(t89)
      t92 = t76*t88*t90
      t93 = dsqrt(t92)
      t96 = 0.6D1+0.5091195559614694D0*t93+0.6480068056560093D-1*t92
      t99 = 0.3D1+0.763679333942204D0*t93+0.972010208484014D-1*t92
      t100 = 1/t99
      t101 = t96*t100
      t102 = t76*t83
      t103 = t87*t64
      t107 = 0.1620017014140023D-1*t102*t103*t90+0.1D1
      t108 = t101*t107
      t111 = 0.7772672717413724D-2*rho*(1.D0*t64+t108)
      zk(i) = zk(i)+2*t111
      t113 = t66**2
      t116 = 1/t113*t75*t83
      t117 = t87*t90
      t118 = 1/t29
      t120 = t116*t117*t118
      t121 = 0.1080011342760016D-1*t120
      t122 = t7**(1.D0/5.D0)
      t123 = t122**2
      t128 = 1/t123/t74/t73*t83
      t131 = dexp(t89-t70)
      t133 = t128*t87*t118*t131
      t134 = 0.4124365791094649D-3*t133
      t137 = t76/t82/t81
      t141 = sigma/t31/t29/rho
      t145 = 1/t30/t41/t29
      t146 = t39*t145
      t148 = -0.5872185841956702D0*t141-0.1179732093959043D-1*t146
      t150 = t137*t117*t148
      t151 = 0.3240034028280047D-1*t150
      t161 = -0.6387160900514977D0*t141*t84+0.7983951125643721D-1
     &*sigma/t31/t41/t66/t7
      t164 = t76*t83*t161*t90
      t165 = 0.1620017014140023D-1*t164
      t167 = t102*t117*t146
      t168 = 0.3822372128681373D-3*t167
      t169 = -t121-t134-t151+t165-t168
      t178 = dexp(-0.2211997676173206D-2*t79)
      t179 = t82/t87*t178
      t182 = 1/t93
      t189 = t182*(-0.5810796994275671D-1*t120-0.2219037096493859D-2
     &*t133-0.1743239098282701D0*t150+0.8716195491413506D-1*t164
     &-0.2056555111688325D-2*t167)
      t199 = t99**2
      t201 = t96/t199
      t228 = t76*t88
      t235 = t64*t90
      t239 = -0.1080011342760016D-1*t116*t103*t118*t90
     &-0.4124365791094649D-3*t128*t103*t118*t131-0.3240034028280047D-1
     &*t137*t103*t148*t90+0.1620017014140023D-1*t102*t161*t64*t90
     &-0.3822372128681373D-3*t228*t64*t39*t145*t90
     &+0.1620017014140023D-1*t102*t87*t169*t235-t121-t134-t151+t165-t168
      vrhoa(i) = vrhoa(i)+0.7772672717413724D-2*rho*(1.D0*t169*t64
     &+0.6172774676263781D2*t169*t84*t74*t179+(0.2920537730383703D1
     &*t189-0.4320045371040062D-1*t120-0.164974631643786D-2*t133
     &-0.1296013611312019D0*t150+0.6480068056560093D-1*t164
     &-0.1528948851472549D-2*t167)*t100*t107-1.D0*t201*t107*
     &(0.4380806595575555D1*t189-0.6480068056560093D-1*t120
     &-0.2474619474656789D-2*t133-0.1944020416968028D0*t150
     &+0.972010208484014D-1*t164-0.2293423277208824D-2*t167)+t101*t239
     &)+0.1554534543482745D-1*t64+0.1554534543482745D-1*t108
      t248 = sigma*t44
      t250 = 0.4404139381467527D0*t33+0.8847990704692823D-2*t248
      t252 = t137*t117*t250
      t253 = 0.3240034028280047D-1*t252
      t254 = t75*t83
      t256 = t254*t33*t90
      t257 = 0.7760481998163131D-2*t256
      t259 = t102*t117*t248
      t260 = 0.286677909651103D-3*t259
      t261 = -t253+t257+t260
      t273 = t182*(-0.1743239098282701D0*t252+0.4175380728300095D-1
     &*t256+0.1542416333766244D-2*t259)
      vsigmaaa(i) = vsigmaaa(i)+0.1554534543482745D-1*rho*(1.D0*t261
     &*t64+0.6172774676263781D2*t261*t84*t74*t179+
     &(0.2920537730383703D1*t273-0.1296013611312019D0*t252
     &+0.3104192799265252D-1*t256+0.1146711638604412D-2*t259)*t100
     &*t107-1.D0*t201*t107*(0.4380806595575555D1*t273
     &-0.1944020416968028D0*t252+0.4656289198897879D-1*t256
     &+0.1720067457906618D-2*t259)+t101*(-0.3240034028280047D-1*t137
     &*t103*t250*t90+0.7760481998163131D-2*t254*t33*t64*t90
     &+0.286677909651103D-3*t228*t64*sigma*t44*t90
     &+0.1620017014140023D-1*t102*t87*t261*t235-t253+t257+t260))
      endif
      endif ! rho
      enddo
      
      else if(ideriv.eq.2) then
      
      do i=1,npt
      rho = dmax1(0.D0,rhoa1(i))
      zk(i) = 0.D0
      vrhoa(i) = 0.D0
      vsigmaaa(i) = 0.D0
      v2rhoa2(i) = 0.D0
      v2rhoasigmaaa(i) = 0.D0
      v2sigmaaa2(i) = 0.D0
      if(rho.gt.tolmin) then
      sigma = dmax1(0.D0,sigmaaa1(i))
      t1 = 1/rho
      t2 = t1**(1.D0/3.D0)
      trsa = 0.1620017014140023D-1*t2
      t3 = t1**(1.D0/15.D0)
      t4 = t3**2
      t6 = dexp(-0.5940962989070748D0*t4)
      t8 = t1**(1.D0/6.D0)
      t10 = dexp(-0.9630597286858919D0*t8)
      t13 = (0.1247511874D1+0.812904345D0*t6-0.859614445D0*t10)**2
      t14 = sigma**2
      t15 = rho**2
      t16 = t15**2
      t18 = rho**(1.D0/3.D0)
      t21 = t14/t18/t16/rho
      t24 = (0.1D1+0.1127183082292567D0*t21)**2
      t26 = t18**2
      t36 = dexp(-0.1127183082292567D0*t21)
      t37 = t36**2
      trcaa = t13*t24/(0.1D1+0.312690210142125D0*sigma/t26/t15/t2)*t37
      if(trsa.lt.tolmax*trcaa) then
      t2 = 1/rho
      t3 = t2**(1.D0/3.D0)
      t4 = t3**2
      t5 = t2**(1.D0/6.D0)
      t8 = 0.8301627136974294D0*t5+0.1354632918982911D1*t3
      t9 = t8**2
      t10 = 1/t9
      t13 = dexp(-0.6108870577108572D0*t4*t10)
      t14 = rho*t13
      t15 = 0.3141592653589793D1**2
      t17 = 2**(1.D0/3.D0)
      t20 = 1/0.3141592653589793D1*t2
      t21 = t20**(1.D0/3.D0)
      t23 = 2**(1.D0/15.D0)
      t24 = t23**2
      t25 = t20**(1.D0/15.D0)
      t26 = t25**2
      t29 = dexp(-0.630966299458536D0*t24*t26)
      t31 = 2**(1.D0/6.D0)
      t32 = t20**(1.D0/6.D0)
      t35 = dexp(-0.1038340679664977D1*t31*t32)
      t38 = (0.1247511874D1+0.812904345D0*t29-0.859614445D0*t35)**2
      t40 = 0.3141592653589793D1**(1.D0/3.D0)
      t43 = sigma**2
      t45 = rho**2
      t46 = t45**2
      t47 = t46*rho
      t48 = rho**(1.D0/3.D0)
      t50 = 1/t48/t47
      t52 = 1/t40/0.3141592653589793D1*t43*t17*t50
      t55 = (0.1D1+0.4116400544093991D0*t52)**2
      t58 = t40**2
      t61 = t48**2
      t63 = 1/t61/t45
      t71 = dexp(0.8232801088187983D0*t52)
      t76 = expei(-0.1858628590577086D0/t15*t17*t21/t38/t55*(0.1D1
     &+0.3634859066227017D0/t58*sigma*t17*t63/t21)*t71)
      t78 = t2**(1.D0/15.D0)
      t79 = t78**2
      t81 = dexp(-0.5940962989070748D0*t79)
      t84 = dexp(-0.9630597286858919D0*t5)
      t86 = 0.1247511874D1+0.812904345D0*t81-0.859614445D0*t84
      t87 = t86**2
      t88 = 1/t87
      t89 = t3*t88
      t90 = t43*t50
      t92 = 0.1D1+0.1127183082292567D0*t90
      t93 = t92**2
      t94 = 1/t93
      t96 = 1/t3
      t99 = 0.1D1+0.312690210142125D0*sigma*t63*t96
      t100 = t94*t99
      t102 = dexp(0.2254366164585135D0*t90)
      t104 = t89*t100*t102
      t105 = dsqrt(t104)
      t108 = 0.6D1+0.5091195559614694D0*t105+0.6480068056560093D-1*t104
      t111 = 0.3D1+0.763679333942204D0*t105+0.972010208484014D-1*t104
      t112 = 1/t111
      t113 = t108*t112
      t114 = t89*t94
      t115 = t99*t76
      t119 = 0.1620017014140023D-1*t114*t115*t102+0.1D1
      t121 = 1.D0*t76+t113*t119
      t123 = 0.7772672717413724D-2*t14*t121
      zk(i) = zk(i)+2*t123
      t125 = t13*t121
      t127 = t96*t10
      t128 = 1/t45
      t132 = 1/t9/t8
      t133 = t4*t132
      t134 = t5**2
      t135 = t134**2
      t136 = t135*t5
      t137 = 1/t136
      t138 = t137*t128
      t140 = 1/t4
      t143 = -0.2767209045658098D0*t138-0.9030886126552743D0*t140*t128
      t146 = 0.8145160769478096D0*t127*t128+0.1221774115421714D1*t133
     &*t143
      t147 = rho*t146
      t150 = t140*t88
      t151 = t150*t94
      t152 = t99*t102
      t154 = t151*t152*t128
      t155 = 0.1080011342760016D-1*t154
      t157 = 1/t87/t86
      t158 = t3*t157
      t159 = t158*t94
      t160 = t79**2
      t162 = t160**2
      t163 = t162*t160*t78
      t164 = 1/t163
      t170 = 0.1287849233946613D0*t164*t128*t81-0.2759533513920579D0
     &*t138*t84
      t172 = t159*t152*t170
      t173 = 0.3240034028280047D-1*t172
      t175 = 1/t93/t92
      t176 = t89*t175
      t177 = t46*t45
      t179 = 1/t48/t177
      t180 = t43*t179
      t181 = t152*t180
      t182 = t176*t181
      t183 = 0.389558564557814D-1*t182
      t184 = t45*rho
      t186 = 1/t61/t184
      t191 = 1/t61/t46
      t192 = sigma*t191
      t194 = 1/t3/t2
      t197 = -0.1667681120758D1*sigma*t186*t96+0.20846014009475D0*t192
     &*t194
      t198 = t94*t197
      t200 = t89*t198*t102
      t201 = 0.1620017014140023D-1*t200
      t202 = t114*t181
      t203 = 0.389558564557814D-1*t202
      t204 = -t155-t173+t183+t201-t203
      t205 = t204*t76
      t207 = -t204
      t208 = t207*t96
      t209 = t208*t87
      t210 = 1/t99
      t211 = t93*t210
      t213 = dexp(-0.2254366164585135D0*t90)
      t214 = t211*t213
      t217 = 1/t105
      t223 = -0.5810796994275671D-1*t154-0.1743239098282701D0*t172
     &+0.2095946261306892D0*t182+0.8716195491413506D-1*t200
     &-0.2095946261306892D0*t202
      t224 = t217*t223
      t231 = 0.2920537730383703D1*t224-0.4320045371040062D-1*t154
     &-0.1296013611312019D0*t172+0.1558234258231256D0*t182
     &+0.6480068056560093D-1*t200-0.1558234258231256D0*t202
      t232 = t231*t112
      t234 = t111**2
      t235 = 1/t234
      t236 = t108*t235
      t243 = 0.4380806595575555D1*t224-0.6480068056560093D-1*t154
     &-0.1944020416968028D0*t172+0.2337351387346884D0*t182
     &+0.972010208484014D-1*t200-0.2337351387346884D0*t202
      t244 = t119*t243
      t247 = t128*t102
      t251 = t170*t102
      t255 = t175*t99
      t256 = t89*t255
      t257 = t76*t43
      t258 = t179*t102
      t259 = t257*t258
      t262 = t197*t76
      t266 = t89*t100
      t270 = t76*t102
      t274 = -0.1080011342760016D-1*t151*t115*t247
     &-0.3240034028280047D-1*t159*t115*t251+0.389558564557814D-1*t256
     &*t259+0.1620017014140023D-1*t114*t262*t102-0.389558564557814D-1
     &*t266*t259+0.1620017014140023D-1*t114*t99*t204*t270-t155-t173
     &+t183+t201-t203
      t276 = 1.D0*t205-0.6172774676263781D2*t209*t214+t232*t119-1.D0
     &*t236*t244+t113*t274
      vrhoa(i) = vrhoa(i)+0.1554534543482745D-1*t125
     &+0.7772672717413724D-2*t147*t125+0.7772672717413724D-2*t14*t276
      t280 = sigma*t50
      t281 = t152*t280
      t282 = t176*t281
      t283 = 0.2921689234183605D-1*t282
      t284 = t88*t94
      t285 = t63*t102
      t286 = t284*t285
      t287 = 0.2026253842341047D-1*t286
      t288 = t114*t281
      t289 = 0.2921689234183605D-1*t288
      t290 = -t283+t287+t289
      t291 = t290*t76
      t293 = -t290
      t294 = t293*t96
      t301 = -0.1571959695980169D0*t282+0.1090187599939973D0*t286
     &+0.1571959695980169D0*t288
      t302 = t217*t301
      t307 = 0.2920537730383703D1*t302-0.1168675693673442D0*t282
     &+0.8105015369364188D-1*t286+0.1168675693673442D0*t288
      t308 = t307*t112
      t314 = 0.4380806595575555D1*t302-0.1753013540510163D0*t282
     &+0.1215752305404628D0*t286+0.1753013540510163D0*t288
      t315 = t119*t314
      t318 = t76*sigma
      t319 = t50*t102
      t320 = t318*t319
      t333 = -0.2921689234183605D-1*t256*t320+0.2026253842341047D-1
     &*t284*t63*t76*t102+0.2921689234183605D-1*t266*t320
     &+0.1620017014140023D-1*t114*t99*t290*t270-t283+t287+t289
      t335 = 1.D0*t291-0.6172774676263781D2*t294*t87*t214+t308*t119
     &-1.D0*t236*t315+t113*t333
      vsigmaaa(i) = vsigmaaa(i)+0.1554534543482745D-1*t14*t335
      t342 = t13*t276
      t345 = 1/t46
      t352 = 1/t184
      t355 = t9**2
      t358 = t143**2
      t363 = 1/t136/t2*t345
      t365 = t137*t352
      t368 = 1/t4/t2
      t369 = t368*t345
      t380 = t146**2
      t386 = t140*t157
      t390 = t386*t94*t152*t128*t170
      t391 = 0.4320045371040062D-1*t390
      t392 = t197*t102
      t394 = t151*t392*t128
      t395 = 0.2160022685520031D-1*t394
      t397 = t151*t152*t352
      t398 = 0.4320045371040062D-1*t397
      t400 = t368*t88*t94
      t402 = t400*t152*t345
      t403 = 0.1440015123680021D-1*t402
      t404 = t150*t175
      t405 = t46**2
      t407 = 1/t48/t405
      t409 = t152*t407*t43
      t410 = t404*t409
      t411 = 0.5194114194104187D-1*t410
      t412 = t151*t409
      t413 = 0.5194114194104187D-1*t412
      t414 = t87**2
      t417 = t3/t414*t94
      t418 = t170**2
      t420 = t417*t152*t418
      t421 = 0.972010208484014D-1*t420
      t423 = t159*t392*t170
      t424 = 0.6480068056560093D-1*t423
      t446 = 0.2232272005507463D0/t163/t2*t345*t81
     &-0.5151396935786452D0*t164*t352*t81+0.2040283902501318D-1/t162
     &/t79/t78/t2*t345*t81-0.4599222523200964D0*t363*t84
     &+0.1103813405568231D1*t365*t84-0.8858651990719928D-1*t369*t84
      t448 = t159*t152*t446
      t449 = 0.3240034028280047D-1*t448
      t450 = t158*t255
      t451 = t251*t180
      t452 = t450*t451
      t453 = 0.1558234258231256D0*t452
      t454 = t158*t100
      t455 = t454*t451
      t456 = 0.1558234258231256D0*t455
      t457 = t392*t180
      t458 = t176*t457
      t459 = 0.779117129115628D-1*t458
      t462 = 1/t48/t46/t184
      t464 = t152*t43*t462
      t465 = t176*t464
      t466 = 0.4934408484398978D0*t465
      t467 = t93**2
      t468 = 1/t467
      t469 = t89*t468
      t470 = t43**2
      t473 = 1/t61/t405/t46
      t475 = t152*t470*t473
      t476 = t469*t475
      t477 = 0.1405132235301584D0*t476
      t478 = t176*t475
      t479 = 0.1873509647068778D0*t478
      t494 = 0.1222966155222533D2*t192*t96-0.3057415388056333D1*sigma
     &/t61/t47*t194+0.5558937069193333D0*sigma/t61/t177/t3/t128
      t497 = t89*t94*t494*t102
      t498 = 0.1620017014140023D-1*t497
      t499 = t114*t457
      t500 = 0.779117129115628D-1*t499
      t501 = t114*t464
      t502 = 0.4934408484398978D0*t501
      t503 = t114*t475
      t504 = 0.9367548235343892D-1*t503
      t505 = t391-t395+t398-t403-t411+t413+t421-t424-t449-t453+t456
     &+t459-t466+t477-t479+t498-t500+t502+t504
      t508 = t204**2
      t511 = t204*t207
      t513 = t87*t93
      t515 = t513*t210*t213
      t523 = t207**2
      t539 = t87*t92
      t540 = t208*t539
      t543 = t210*t43*t179*t213
      t546 = t99**2
      t548 = t93/t546
      t553 = t208*t513
      t577 = 1/t105/t104
      t578 = t223**2
      t579 = t577*t578
      t600 = 0.2324318797710268D0*t390-0.1162159398855134D0*t394
     &+0.2324318797710268D0*t397-0.7747729325700894D-1*t402
     &-0.2794595015075856D0*t410+0.2794595015075856D0*t412
     &+0.5229717294848104D0*t420-0.3486478196565402D0*t423
     &-0.1743239098282701D0*t448-0.8383785045227567D0*t452
     &+0.8383785045227567D0*t455+0.4191892522613783D0*t458
     &-0.2654865264322063D1*t465+0.7560048534846353D0*t476
     &-0.1008006471312847D1*t478+0.8716195491413506D-1*t497
     &-0.4191892522613783D0*t499+0.2654865264322063D1*t501
     &+0.5040032356564235D0*t503
      t601 = t217*t600
      t603 = 0.6232937032925024D0*t455+0.3116468516462512D0*t458
     &-0.1973763393759591D1*t465+0.5620528941206335D0*t476
     &-0.7494038588275113D0*t478+0.6480068056560093D-1*t497
     &-0.3116468516462512D0*t499+0.1973763393759591D1*t501
     &+0.3747019294137557D0*t503-0.1675351208713011D2*t579
     &+0.2920537730383703D1*t601
      t607 = t231*t235
      t614 = t108/t234/t111
      t615 = t243**2
      t644 = 0.9349405549387536D0*t455+0.4674702774693768D0*t458
     &-0.2960645090639387D1*t465+0.8430793411809502D0*t476
     &-0.1124105788241267D1*t478+0.972010208484014D-1*t497
     &-0.4674702774693768D0*t499+0.2960645090639387D1*t501
     &+0.5620528941206335D0*t503-0.2513026813069517D2*t579
     &+0.4380806595575555D1*t601
      t659 = 1/t86
      t666 = t150*t100
      t670 = 1/t92
      t675 = t89*t175*t197
      t679 = t257*t462*t102
      t682 = -0.6480068056560093D-1*t159*t262*t251
     &+0.972010208484014D-1*t417*t115*t418*t102-0.6480068056560093D-1
     &*t454*t205*t251+0.2D1*t659*t207*t170-0.1D1*t197*t207*t210
     &-0.2160022685520031D-1*t666*t205*t247-t403-t411-t395
     &-0.2404657242224144D1*t670*t207*t180+t391+0.779117129115628D-1
     &*t675*t259+0.4934408484398978D0*t266*t679
      t683 = t76*t170
      t684 = t180*t102
      t685 = t683*t684
      t696 = t76*t470*t473*t102
      t702 = t150*t255
      t705 = t76*t407*t43*t102
      t708 = t89*t198
      t711 = t398+t413+t421-0.1558234258231256D0*t450*t685
     &+0.1620017014140023D-1*t114*t99*t505*t270+0.1558234258231256D0
     &*t454*t685+t477+0.9367548235343892D-1*t266*t696
     &+0.2404657242224144D1*t207*t43*t179-0.5194114194104187D-1*t702
     &*t705+t459-t466-0.779117129115628D-1*t708*t259
      t716 = t89*t468*t99
      t723 = -t449-t453+t456-t424-t500+t502+t504-0.4934408484398978D0
     &*t256*t679+0.1405132235301584D0*t716*t696-t479
     &-0.1873509647068778D0*t256*t696+t498+0.5194114194104187D-1*t666
     &*t705
      t724 = t205*t684
      t765 = 0.779117129115628D-1*t256*t724-0.779117129115628D-1*t266
     &*t724-0.1D1*t523-0.2160022685520031D-1*t151*t262*t247
     &+0.4320045371040062D-1*t386*t100*t76*t128*t251
     &+0.1620017014140023D-1*t114*t99*t508*t270-0.2D1*t511
     &+0.6666666666666667D0*t2*t207+0.4320045371040062D-1*t151*t115
     &*t352*t102-0.1440015123680021D-1*t400*t115*t345*t102
     &+0.3240034028280047D-1*t114*t197*t204*t270+0.1620017014140023D-1
     &*t114*t494*t76*t102-0.3240034028280047D-1*t159*t115*t446*t102
      s1 = 1.D0*t505*t76+1.D0*t508*t76-0.1234554935252756D3*t511*t96
     &*t515+0.6172774676263781D2*t505*t96*t87*t214
     &-0.6172774676263781D2*t523*t96*t87*t214-0.4115183117509188D2
     &*t207*t194*t87*t211*t128*t213-0.1234554935252756D3*t208*t86*t211
     &*t170*t213+0.148434073298955D3*t540*t543
      t769 = s1+0.6172774676263781D2*t209*t548*t197*t213
     &-0.148434073298955D3*t553*t543+(0.1728018148416025D0*t390
     &-0.8640090742080124D-1*t394+0.1728018148416025D0*t397
     &-0.5760060494720083D-1*t402-0.2077645677641675D0*t410
     &+0.2077645677641675D0*t412+0.3888040833936056D0*t420
     &-0.2592027222624037D0*t423-0.1296013611312019D0*t448
     &-0.6232937032925024D0*t452+t603)*t112*t119-2.D0*t607*t244+2.D0
     &*t232*t274+2.D0*t614*t119*t615-2.D0*t236*t274*t243-1.D0*t236
     &*t119*(0.2592027222624037D0*t390-0.1296013611312019D0*t394
     &+0.2592027222624037D0*t397-0.8640090742080124D-1*t402
     &-0.3116468516462512D0*t410+0.3116468516462512D0*t412
     &+0.5832061250904084D0*t420-0.3888040833936056D0*t423
     &-0.1944020416968028D0*t448-0.9349405549387536D0*t452+t644)+t113*
     &(t682+t711+t723+t765)
      v2rhoa2(i) = v2rhoa2(i)+0.310906908696549D-1*t146*t13*t121
     &+0.310906908696549D-1*t342+0.7772672717413724D-2*rho*
     &(0.5430107179652064D0*t194*t10*t345-0.3258064307791238D1*t96
     &*t132*t128*t143-0.3258064307791238D1*t127*t352
     &-0.3665322346265143D1*t4/t355*t358+0.1221774115421714D1*t133*(
     &-0.461201507609683D0*t363+0.1106883618263239D1*t365
     &-0.1204118150207032D1*t369+0.3612354450621097D1*t140*t352))*t125
     &+0.7772672717413724D-2*rho*t380*t125+0.1554534543482745D-1*t147
     &*t342+0.7772672717413724D-2*t14*t769
      t773 = t13*t335
      t778 = t152*t462*sigma
      t779 = t404*t778
      t780 = 0.194779282278907D-1*t779
      t781 = t186*t88
      t783 = t781*t94*t102
      t784 = 0.1350835894894031D-1*t783
      t785 = t151*t778
      t786 = 0.194779282278907D-1*t785
      t787 = t251*t280
      t788 = t450*t787
      t789 = 0.584337846836721D-1*t788
      t790 = t157*t94
      t792 = t790*t285*t170
      t793 = 0.4052507684682094D-1*t792
      t794 = t454*t787
      t795 = 0.584337846836721D-1*t794
      t796 = t43*sigma
      t799 = 1/t61/t405/t184
      t801 = t152*t796*t799
      t802 = t469*t801
      t803 = 0.1053849176476188D0*t802
      t804 = t88*t175
      t806 = 1/t405/rho
      t808 = t806*t102*t43
      t809 = t804*t808
      t810 = 0.4872445976569897D-1*t809
      t811 = t176*t801
      t812 = 0.1405132235301584D0*t811
      t814 = t152*sigma*t179
      t815 = t176*t814
      t816 = 0.3116468516462512D0*t815
      t817 = t392*t280
      t818 = t176*t817
      t819 = 0.2921689234183605D-1*t818
      t824 = -0.6670724483032D1*t186*t96+0.8338405603789999D0*t191*t194
      t827 = t89*t94*t824*t102
      t828 = 0.1620017014140023D-1*t827
      t829 = t114*t817
      t830 = 0.2921689234183605D-1*t829
      t831 = t284*t808
      t832 = 0.4872445976569897D-1*t831
      t833 = t114*t801
      t834 = 0.7025661176507919D-1*t833
      t835 = t114*t814
      t836 = 0.3116468516462512D0*t835
      t837 = t780-t784-t786+t789-t793-t795-t803+t810+t812+t816-t819
     &+t828+t830-t832-t834-t836
      t840 = t204*t290
      t843 = t204*t293
      t847 = t290*t207
      t856 = t207*t293
      t862 = t210*sigma*t50*t213
      t868 = t548*t63*t213
      t874 = t577*t223*t301
      t892 = 0.1047973130653446D0*t779-0.7267917332933151D-1*t783
     &-0.1047973130653446D0*t785+0.3143919391960338D0*t788
     &-0.2180375199879945D0*t792-0.3143919391960338D0*t794
     &-0.5670036401134764D0*t802+0.2621527507578613D0*t809
     &+0.7560048534846353D0*t811+0.1676757009045513D1*t815
     &-0.1571959695980169D0*t818+0.8716195491413506D-1*t827
     &+0.1571959695980169D0*t829-0.2621527507578613D0*t831
     &-0.3780024267423176D0*t833-0.1676757009045513D1*t835
      t893 = t217*t892
      t911 = -0.1675351208713011D2*t874+0.2920537730383703D1*t893
     &+0.779117129115628D-1*t779-0.5403343579576125D-1*t783
     &-0.779117129115628D-1*t785+0.2337351387346884D0*t788
     &-0.1621003073872838D0*t792-0.2337351387346884D0*t794
     &-0.4215396705904751D0*t802+0.1948978390627959D0*t809
     &+0.5620528941206335D0*t811+0.1246587406585005D1*t815
     &-0.1168675693673442D0*t818+0.6480068056560093D-1*t827
     &+0.1168675693673442D0*t829-0.1948978390627959D0*t831
     &-0.2810264470603167D0*t833-0.1246587406585005D1*t835
      t917 = t307*t235
      t944 = -0.2513026813069517D2*t874+0.4380806595575555D1*t893
     &+0.1168675693673442D0*t779-0.8105015369364188D-1*t783
     &-0.1168675693673442D0*t785+0.3506027081020326D0*t788
     &-0.2431504610809256D0*t792-0.3506027081020326D0*t794
     &-0.6323095058857127D0*t802+0.2923467585941938D0*t809
     &+0.8430793411809502D0*t811+0.1869881109877507D1*t815
     &-0.1753013540510163D0*t818+0.972010208484014D-1*t827
     &+0.1753013540510163D0*t829-0.2923467585941938D0*t831
     &-0.4215396705904751D0*t833-0.1869881109877507D1*t835
      t958 = t670*t293
      t990 = 0.2D1*t659*t293*t170-0.1D1*t197*t293*t210-t784
     &-0.2404657242224144D1*t958*t180+0.1620017014140023D-1*t114*t197
     &*t290*t270+0.1620017014140023D-1*t114*t824*t76*t102
     &+0.1620017014140023D-1*t114*t99*t837*t270-0.1350835894894031D-1
     &*t781*t94*t76*t102-0.3240034028280047D-1*t454*t291*t251
     &-0.1080011342760016D-1*t666*t291*t247-0.4052507684682094D-1*t790
     &*t63*t683*t102+0.1620017014140023D-1*t266*t840*t270
      t993 = t284*t63
      t999 = t76*t796*t799*t102
      t1005 = t280*t102
      t1006 = t683*t1005
      t1011 = t76*t462*sigma*t102
      t1015 = t257*t102
      t1018 = t205*t1005
      t1021 = -0.1D1*t843-0.1D1*t847+0.2026253842341047D-1*t993*t205
     &*t102+0.1405132235301584D0*t256*t999-0.1D1*t856
     &+0.6666666666666667D0*t2*t293+0.584337846836721D-1*t450*t1006
     &-0.194779282278907D-1*t666*t1011-t793+0.4872445976569897D-1*t804
     &*t806*t1015-t795+0.2921689234183605D-1*t266*t1018-t803
      t1025 = t318*t258
      t1028 = t291*t684
      t1031 = t810+t812+t789+t780-0.2921689234183605D-1*t256*t1018
     &-t786-t836+0.3116468516462512D0*t256*t1025+0.389558564557814D-1
     &*t256*t1028+t830-t832-t834
      t1054 = t816+0.2921689234183605D-1*t708*t320-t819
     &-0.584337846836721D-1*t454*t1006+t828-0.7025661176507919D-1*t266
     &*t999-0.3116468516462512D0*t266*t1025-0.4872445976569897D-1*t284
     &*t806*t1015-0.1053849176476188D0*t716*t999+0.2404657242224144D1
     &*t293*t43*t179-0.2921689234183605D-1*t675*t320
     &+0.194779282278907D-1*t702*t1011-0.389558564557814D-1*t266*t1028
      t1058 = 1.D0*t837*t76+1.D0*t840*t76-0.6172774676263781D2*t843
     &*t96*t515-0.6172774676263781D2*t847*t96*t515
     &+0.6172774676263781D2*t837*t96*t87*t214-0.6172774676263781D2
     &*t856*t96*t515-0.1113255549742162D3*t540*t862
     &+0.7720664842723637D2*t207*t140*t87*t868+0.1113255549742162D3
     &*t553*t862+t911*t112*t119-1.D0*t607*t315+t232*t333-1.D0*t917
     &*t244+2.D0*t614*t244*t314-1.D0*t236*t333*t243-1.D0*t236*t119
     &*t944+t308*t274-1.D0*t236*t274*t314+t113*(t990+t1021+t1031+t1054)
      v2rhoasigmaaa(i) = v2rhoasigmaaa(i)+0.1554534543482745D-1*t773
     &+0.7772672717413724D-2*t147*t773+0.7772672717413724D-2*t14*t1058
      t1064 = 1/t61/t405/t45
      t1066 = t152*t43*t1064
      t1067 = t469*t1066
      t1068 = 0.7903868823571409D-1*t1067
      t1069 = 1/t405
      t1071 = t1069*t102*sigma
      t1072 = t804*t1071
      t1073 = 0.7308668964854846D-1*t1072
      t1074 = t176*t1066
      t1075 = 0.1053849176476188D0*t1074
      t1076 = t152*t50
      t1077 = t176*t1076
      t1078 = 0.1168675693673442D0*t1077
      t1079 = t284*t1071
      t1080 = 0.7308668964854846D-1*t1079
      t1081 = t114*t1066
      t1082 = 0.5269245882380939D-1*t1081
      t1083 = t114*t1076
      t1084 = 0.1168675693673442D0*t1083
      t1085 = t1068-t1073-t1075-t1078+t1080+t1082+t1084
      t1088 = t290**2
      t1091 = t290*t293
      t1100 = t293**2
      t1115 = t301**2
      t1116 = t577*t1115
      t1126 = t217*(0.4252527300851073D0*t1067-0.3932291261367919D0
     &*t1072-0.5670036401134764D0*t1074-0.6287838783920675D0*t1077
     &+0.3932291261367919D0*t1079+0.2835018200567382D0*t1081
     &+0.6287838783920675D0*t1083)
      t1142 = t314**2
      t1163 = t318*t102
      t1167 = t257*t1064*t102
      t1170 = t291*t1005
      t1175 = t115*t319
      t1190 = -0.7308668964854846D-1*t804*t1069*t1163
     &-0.1053849176476188D0*t256*t1167-0.584337846836721D-1*t256*t1170
     &+0.1803492931668108D1*t958*t280-0.1168675693673442D0*t176*t1175
     &+0.7903868823571409D-1*t716*t1167-0.12507608405685D1*t63*t293
     &*t96*t210+0.1620017014140023D-1*t114*t99*t1088*t270-0.2D1*t1091
     &-0.1D1*t1100+t1068-t1073
      t1210 = -t1075-t1078+t1080+t1082+t1084+0.1620017014140023D-1
     &*t114*t99*t1085*t270+0.5269245882380939D-1*t266*t1167
     &+0.584337846836721D-1*t266*t1170-0.1803492931668108D1*t293*sigma
     &*t50+0.1168675693673442D0*t114*t1175+0.7308668964854846D-1*t284
     &*t1069*t1163+0.4052507684682094D-1*t993*t291*t102
      s1 = 1.D0*t1085*t76+1.D0*t1088*t76-0.1234554935252756D3*t1091
     &*t96*t515+0.6172774676263781D2*t1085*t96*t87*t214
     &-0.6172774676263781D2*t1100*t96*t87*t214-0.1113255549742162D3
     &*t294*t539*t862+0.7720664842723637D2*t293*t140*t87*t868
      t1213 = s1+0.1113255549742162D3*t294*t513*t862+(
     &-0.1675351208713011D2*t1116+0.2920537730383703D1*t1126
     &+0.3161547529428563D0*t1067-0.2923467585941938D0*t1072
     &-0.4215396705904751D0*t1074-0.4674702774693768D0*t1077
     &+0.2923467585941938D0*t1079+0.2107698352952376D0*t1081
     &+0.4674702774693768D0*t1083)*t112*t119-2.D0*t917*t315+2.D0*t308
     &*t333+2.D0*t614*t119*t1142-2.D0*t236*t333*t314-1.D0*t236*t119*(
     &-0.2513026813069517D2*t1116+0.4380806595575555D1*t1126
     &+0.4742321294142845D0*t1067-0.4385201378912907D0*t1072
     &-0.6323095058857127D0*t1074-0.7012054162040652D0*t1077
     &+0.4385201378912907D0*t1079+0.3161547529428563D0*t1081
     &+0.7012054162040652D0*t1083)+t113*(t1190+t1210)
      v2sigmaaa2(i) = v2sigmaaa2(i)+0.1554534543482745D-1*t14*t1213
      endif
      t1 = 1/rho
      t2 = t1**(1.D0/3.D0)
      trsa = 0.1620017014140023D-1*t2
      t3 = t1**(1.D0/15.D0)
      t4 = t3**2
      t5 = t4**2
      t7 = dexp(-0.6837598574312962D-1*t5)
      t10 = (0.942486901D0+0.349064173D0*t7)**2
      t11 = rho**2
      t12 = rho**(1.D0/3.D0)
      t13 = t12**2
      t16 = sigma/t13/t11
      t18 = sigma**2
      t19 = t11**2
      t23 = t18/t12/t19/rho
      t26 = (0.1D1+0.1101034845366882D0*t16+0.1105998838086603D-2*t23)
     &**2
      t34 = dexp(-0.1105998838086603D-2*t23)
      t35 = t34**2
      trcba = t10*t26/(0.1D1+0.1197592668846558D0*t16/t2)*t35
      if(trsa.lt.tolmax*trcba) then
      t2 = 0.3141592653589793D1**2
      t4 = 2**(1.D0/3.D0)
      t7 = 1/rho
      t8 = 1/0.3141592653589793D1*t7
      t9 = t8**(1.D0/3.D0)
      t11 = 2**(1.D0/15.D0)
      t12 = t11**2
      t13 = t12**2
      t14 = t8**(1.D0/15.D0)
      t15 = t14**2
      t16 = t15**2
      t19 = dexp(-0.7712625328179681D-1*t13*t16)
      t22 = (0.942486901D0+0.349064173D0*t19)**2
      t24 = 0.3141592653589793D1**(1.D0/3.D0)
      t25 = t24**2
      t27 = 1/t25*sigma
      t28 = t4**2
      t29 = rho**2
      t30 = rho**(1.D0/3.D0)
      t31 = t30**2
      t33 = 1/t31/t29
      t39 = sigma**2
      t41 = t29**2
      t42 = t41*rho
      t44 = 1/t30/t42
      t46 = 1/t24/0.3141592653589793D1*t39*t4*t44
      t49 = (0.1D1+0.1487810599361293D0*t27*t28*t33
     &+0.4039037038781002D-2*t46)**2
      t59 = dexp(0.8078074077562004D-2*t46)
      t64 = expei(-0.1858628590577086D0/t2*t4*t9/t22/t49*(0.1D1
     &+0.1392138426088027D0*t27*t4*t33/t9)*t59)
      t66 = t7**(1.D0/3.D0)
      t67 = t7**(1.D0/15.D0)
      t68 = t67**2
      t69 = t68**2
      t70 = 0.6837598574312962D-1*t69
      t71 = dexp(-t70)
      t73 = 0.942486901D0+0.349064173D0*t71
      t74 = t73**2
      t75 = 1/t74
      t76 = t66*t75
      t77 = sigma*t33
      t79 = t39*t44
      t81 = 0.1D1+0.1101034845366882D0*t77+0.1105998838086603D-2*t79
      t82 = t81**2
      t83 = 1/t82
      t84 = 1/t66
      t87 = 0.1D1+0.1197592668846558D0*t77*t84
      t88 = t83*t87
      t89 = 0.2211997676173206D-2*t79
      t90 = dexp(t89)
      t92 = t76*t88*t90
      t93 = dsqrt(t92)
      t96 = 0.6D1+0.5091195559614694D0*t93+0.6480068056560093D-1*t92
      t99 = 0.3D1+0.763679333942204D0*t93+0.972010208484014D-1*t92
      t100 = 1/t99
      t101 = t96*t100
      t102 = t76*t83
      t103 = t87*t64
      t107 = 0.1620017014140023D-1*t102*t103*t90+0.1D1
      t108 = t101*t107
      t111 = 0.7772672717413724D-2*rho*(1.D0*t64+t108)
      zk(i) = zk(i)+2*t111
      t113 = t66**2
      t114 = 1/t113
      t115 = t114*t75
      t116 = t115*t83
      t117 = t87*t90
      t118 = 1/t29
      t120 = t116*t117*t118
      t121 = 0.1080011342760016D-1*t120
      t122 = t7**(1.D0/5.D0)
      t123 = t122**2
      t126 = 1/t74/t73
      t127 = 1/t123*t126
      t128 = t127*t83
      t129 = t87*t118
      t131 = dexp(t89-t70)
      t133 = t128*t129*t131
      t134 = 0.4124365791094649D-3*t133
      t136 = 1/t82/t81
      t137 = t76*t136
      t138 = t29*rho
      t140 = 1/t31/t138
      t141 = sigma*t140
      t143 = t41*t29
      t145 = 1/t30/t143
      t146 = t39*t145
      t148 = -0.5872185841956702D0*t141-0.1179732093959043D-1*t146
      t150 = t137*t117*t148
      t151 = 0.3240034028280047D-1*t150
      t155 = 1/t31/t41
      t156 = sigma*t155
      t158 = 1/t66/t7
      t161 = -0.6387160900514977D0*t141*t84+0.7983951125643721D-1*t156
     &*t158
      t162 = t83*t161
      t164 = t76*t162*t90
      t165 = 0.1620017014140023D-1*t164
      t167 = t102*t117*t146
      t168 = 0.3822372128681373D-3*t167
      t169 = -t121-t134-t151+t165-t168
      t170 = t169*t64
      t172 = -t169
      t173 = t172*t84
      t174 = t173*t74
      t175 = 1/t87
      t176 = t82*t175
      t177 = 0.2211997676173206D-2*t79
      t178 = dexp(-t177)
      t179 = t176*t178
      t180 = t174*t179
      t182 = 1/t93
      t188 = -0.5810796994275671D-1*t120-0.2219037096493859D-2*t133
     &-0.1743239098282701D0*t150+0.8716195491413506D-1*t164
     &-0.2056555111688325D-2*t167
      t189 = t182*t188
      t196 = 0.2920537730383703D1*t189-0.4320045371040062D-1*t120
     &-0.164974631643786D-2*t133-0.1296013611312019D0*t150
     &+0.6480068056560093D-1*t164-0.1528948851472549D-2*t167
      t197 = t196*t100
      t198 = t197*t107
      t199 = t99**2
      t200 = 1/t199
      t201 = t96*t200
      t208 = 0.4380806595575555D1*t189-0.6480068056560093D-1*t120
     &-0.2474619474656789D-2*t133-0.1944020416968028D0*t150
     &+0.972010208484014D-1*t164-0.2293423277208824D-2*t167
      t209 = t107*t208
      t210 = t201*t209
      t212 = t118*t90
      t216 = t118*t131
      t220 = t148*t90
      t224 = t161*t64
      t228 = t76*t88
      t229 = t64*t39
      t230 = t145*t90
      t231 = t229*t230
      t235 = t64*t90
      t239 = -0.1080011342760016D-1*t116*t103*t212
     &-0.4124365791094649D-3*t128*t103*t216-0.3240034028280047D-1*t137
     &*t103*t220+0.1620017014140023D-1*t102*t224*t90
     &-0.3822372128681373D-3*t228*t231+0.1620017014140023D-1*t102*t87
     &*t169*t235-t121-t134-t151+t165-t168
      t240 = t101*t239
      vrhoa(i) = vrhoa(i)+0.7772672717413724D-2*rho*(1.D0*t170
     &-0.6172774676263781D2*t180+t198-1.D0*t210+t240)
     &+0.1554534543482745D-1*t64+0.1554534543482745D-1*t108
      t248 = sigma*t44
      t250 = 0.4404139381467527D0*t33+0.8847990704692823D-2*t248
      t252 = t137*t117*t250
      t253 = 0.3240034028280047D-1*t252
      t254 = t75*t83
      t255 = t33*t90
      t256 = t254*t255
      t257 = 0.7760481998163131D-2*t256
      t259 = t102*t117*t248
      t260 = 0.286677909651103D-3*t259
      t261 = -t253+t257+t260
      t262 = t261*t64
      t264 = -t261
      t265 = t264*t84
      t266 = t265*t74
      t267 = t266*t179
      t272 = -0.1743239098282701D0*t252+0.4175380728300095D-1*t256
     &+0.1542416333766244D-2*t259
      t273 = t182*t272
      t278 = 0.2920537730383703D1*t273-0.1296013611312019D0*t252
     &+0.3104192799265252D-1*t256+0.1146711638604412D-2*t259
      t279 = t278*t100
      t280 = t279*t107
      t285 = 0.4380806595575555D1*t273-0.1944020416968028D0*t252
     &+0.4656289198897879D-1*t256+0.1720067457906618D-2*t259
      t286 = t107*t285
      t287 = t201*t286
      t289 = t250*t90
      t297 = t64*sigma
      t298 = t44*t90
      t299 = t297*t298
      t306 = -0.3240034028280047D-1*t137*t103*t289
     &+0.7760481998163131D-2*t254*t33*t64*t90+0.286677909651103D-3
     &*t228*t299+0.1620017014140023D-1*t102*t87*t261*t235-t253+t257+t260
      t307 = t101*t306
      vsigmaaa(i) = vsigmaaa(i)+0.1554534543482745D-1*rho*(1.D0*t262
     &-0.6172774676263781D2*t267+t280-1.D0*t287+t307)
      t314 = 1/t30/t41/t138
      t315 = t39*t314
      t317 = t102*t117*t315
      t318 = 0.4841671362996405D-2*t317
      t319 = t39**2
      t320 = t41**2
      t323 = 1/t31/t320/t41
      t326 = t102*t117*t319*t323
      t327 = 0.9018750150519921D-5*t326
      t342 = 0.4683917993710983D1*t156*t84-0.1170979498427746D1*sigma
     &/t31/t42*t158+0.2129053633504992D0*sigma/t31/t143/t66/t118
      t345 = t76*t83*t342*t90
      t346 = 0.1620017014140023D-1*t345
      t347 = t161*t90
      t349 = t102*t347*t146
      t350 = 0.7644744257362745D-3*t349
      t351 = t136*t87
      t352 = t76*t351
      t354 = t352*t220*t146
      t355 = 0.1528948851472549D-2*t354
      t357 = 1/t68/t7
      t359 = t357*t126*t83
      t360 = 1/t41
      t361 = t87*t360
      t362 = t361*t131
      t363 = t359*t362
      t364 = 0.1504040408167142D-4*t363
      t365 = t82**2
      t366 = 1/t365
      t367 = t76*t366
      t368 = t148**2
      t370 = t367*t117*t368
      t371 = 0.972010208484014D-1*t370
      t373 = t137*t347*t148
      t374 = 0.6480068056560093D-1*t373
      t377 = 0.4306269617434915D1*t156+0.1494327319014788D0*t315
      t379 = t137*t117*t377
      t380 = 0.3240034028280047D-1*t379
      t381 = 1/t138
      t384 = t128*t87*t381*t131
      t385 = 0.164974631643786D-2*t384
      t386 = t74**2
      t389 = t357/t386*t83
      t392 = dexp(t89-0.1367519714862592D0*t69)
      t394 = t389*t361*t392
      t395 = 0.1575019863706338D-4*t394
      t397 = 1/t30/t320
      t399 = t39*t131
      t401 = t128*t87*t397*t399
      t402 = 0.1946258676392454D-4*t401
      t403 = t127*t136
      t404 = t148*t131
      t406 = t403*t129*t404
      t407 = 0.164974631643786D-2*t406
      t410 = t128*t161*t118*t131
      t411 = 0.8248731582189298D-3*t410
      t412 = t115*t136
      t415 = t412*t117*t118*t148
      t416 = 0.4320045371040062D-1*t415
      t418 = t116*t347*t118
      t419 = 0.2160022685520031D-1*t418
      t421 = t116*t117*t381
      t422 = 0.4320045371040062D-1*t421
      t426 = 1/t113/t7*t75*t83
      t428 = t426*t117*t360
      t429 = 0.1440015123680021D-1*t428
      t433 = 1/t123/t7*t126*t83
      t434 = t433*t362
      t435 = 0.5499154388126198D-4*t434
      t438 = t116*t117*t397*t39
      t439 = 0.5096496171575164D-3*t438
      t440 = t318+t327+t346-t350+t355-t364+t371-t374-t380+t385+t395
     &+t402+t407-t411+t416-t419+t422-t429-t435+t439
      t443 = t169**2
      t446 = t169*t172
      t448 = t74*t82
      t450 = t448*t175*t178
      t458 = t172**2
      t474 = dexp(-t177-t70)
      t479 = t81*t175
      t484 = t87**2
      t486 = t82/t484
      t491 = t173*t448
      t517 = 0.2604969808138545D-1*t317+0.4852368136508483D-4*t326
     &+0.8716195491413506D-1*t345-0.411311022337665D-2*t349
     &+0.8226220446753301D-2*t354-0.809220527324479D-4*t363
     &+0.5229717294848104D0*t370-0.3486478196565402D0*t373
     &-0.1743239098282701D0*t379+0.8876148385975435D-2*t384
     &+0.8474096824354295D-4*t394+0.1047147712167798D-3*t401
     &+0.8876148385975435D-2*t406-0.4438074192987717D-2*t410
     &+0.2324318797710268D0*t415-0.1162159398855134D0*t418
     &+0.2324318797710268D0*t421-0.7747729325700894D-1*t428
     &-0.2958716128658478D-3*t434+0.27420734822511D-2*t438
      t518 = t182*t517
      t530 = 0.2920537730383703D1*t518+0.1936668545198562D-1*t317
     &+0.3607500060207969D-4*t326+0.6480068056560093D-1*t345
     &-0.3057897702945098D-2*t349+0.6115795405890196D-2*t354
     &-0.601616163266857D-4*t363+0.3888040833936056D0*t370
     &-0.2592027222624037D0*t373-0.1296013611312019D0*t379
     &+0.6598985265751438D-2*t384
      t542 = 1/t93/t92
      t543 = t188**2
      t544 = t542*t543
      t546 = 0.6300079454825353D-4*t394+0.7785034705569816D-4*t401
     &+0.6598985265751438D-2*t406-0.3299492632875719D-2*t410
     &+0.1728018148416025D0*t415-0.8640090742080124D-1*t418
     &+0.1728018148416025D0*t421-0.5760060494720083D-1*t428
     &-0.2199661755250479D-3*t434+0.2038598468630065D-2*t438
     &-0.1675351208713011D2*t544
      t550 = t196*t200
      t557 = t96/t199/t99
      t558 = t208**2
      t576 = 0.4380806595575555D1*t518+0.2905002817797843D-1*t317
     &+0.5411250090311953D-4*t326+0.972010208484014D-1*t345
     &-0.4586846554417647D-2*t349+0.9173693108835294D-2*t354
     &-0.9024242449002855D-4*t363+0.5832061250904084D0*t370
     &-0.3888040833936056D0*t373-0.1944020416968028D0*t379
     &+0.9898477898627157D-2*t384
      t588 = 0.9450119182238029D-4*t394+0.1167755205835472D-3*t401
     &+0.9898477898627157D-2*t406-0.4949238949313579D-2*t410
     &+0.2592027222624037D0*t415-0.1296013611312019D0*t418
     &+0.2592027222624037D0*t421-0.8640090742080124D-1*t428
     &-0.3299492632875719D-3*t434+0.3057897702945098D-2*t438
     &-0.2513026813069517D2*t544
      t596 = 1/t81
      t608 = t103*t360*t131
      t611 = t115*t351
      t612 = t64*t118
      t619 = -0.1D1*t161*t172*t175+0.2D1*t596*t172*t148-t364
     &-0.8248731582189298D-3*t128*t224*t216-0.1440015123680021D-1*t426
     &*t103*t360*t90+t355-0.5499154388126198D-4*t433*t608-t350
     &+0.4320045371040062D-1*t611*t612*t220+t346-0.2160022685520031D-1
     &*t116*t224*t212+t327+t318
      t633 = t402+t395+t385-t380+0.1620017014140023D-1*t102*t342*t64
     &*t90-t374+0.972010208484014D-1*t367*t103*t368*t90+t371+t422-t419
     &-0.1504040408167142D-4*t359*t608+t416-t411-0.6480068056560093D-1
     &*t137*t224*t220
      t636 = t69**2
      t638 = 1/t636/t68/t67
      t640 = t638/t73
      t670 = t127*t351
      t674 = t407+0.2545878071091769D-1*t640*t172*t118*t71
     &+0.4320045371040062D-1*t116*t103*t381*t90+0.1620017014140023D-1
     &*t102*t87*t443*t235-t435+t439-t429+0.1575019863706338D-4*t389
     &*t103*t360*t392-0.3240034028280047D-1*t137*t103*t377*t90-0.1D1
     &*t458+0.1620017014140023D-1*t102*t87*t440*t235
     &+0.164974631643786D-2*t128*t103*t381*t131+0.164974631643786D-2
     &*t670*t612*t404
      t679 = t115*t88
      t686 = t127*t88
      t693 = t64*t397
      t694 = t39*t90
      t713 = t76*t162
      t716 = t64*t148
      t717 = t146*t90
      t724 = 0.3240034028280047D-1*t102*t161*t169*t235
     &-0.2160022685520031D-1*t679*t170*t212-0.6480068056560093D-1*t352
     &*t170*t220-0.8248731582189298D-3*t686*t170*t216
     &+0.2359464187918086D-1*t172*t39*t145+0.5096496171575164D-3*t679
     &*t693*t694+0.6666666666666667D0*t7*t172-0.2D1*t446
     &+0.1946258676392454D-4*t686*t693*t399+0.9018750150519921D-5*t228
     &*t64*t319*t323*t90+0.4841671362996405D-2*t228*t229*t314*t90
     &-0.7644744257362745D-3*t713*t231+0.1528948851472549D-2*t352*t716
     &*t717-0.7644744257362745D-3*t228*t170*t717
      t728 = 1.D0*t440*t64+1.D0*t443*t64-0.1234554935252756D3*t446*t84
     &*t450+0.6172774676263781D2*t440*t84*t74*t179
     &-0.6172774676263781D2*t458*t84*t74*t179-0.4115183117509188D2
     &*t172*t158*t74*t176*t118*t178-0.1571513168609055D1*t172/t67/t7
     &*t73*t176*t118*t474-0.1234554935252756D3*t174*t479*t148*t178
     &+0.6172774676263781D2*t174*t486*t161*t178-0.1456444078873205D1
     &*t491*t175*t39*t145*t178+(t530+t546)*t100*t107-2.D0*t550*t209
     &+2.D0*t197*t239+2.D0*t557*t107*t558-2.D0*t201*t239*t208-1.D0
     &*t201*t107*(t576+t588)+t101*(t619+t633+t674+t724)
      v2rhoa2(i) = v2rhoa2(i)+0.7772672717413724D-2*rho*t728
     &+0.310906908696549D-1*t170-0.1919158292677513D1*t180
     &+0.310906908696549D-1*t198-0.310906908696549D-1*t210
     &+0.310906908696549D-1*t240
      t739 = t412*t117*t118*t250
      t740 = 0.2160022685520031D-1*t739
      t741 = t140*t75
      t743 = t741*t83*t90
      t744 = 0.5173654665442087D-2*t743
      t747 = t116*t117*t314*sigma
      t748 = 0.1911186064340686D-3*t747
      t749 = t250*t131
      t751 = t403*t129*t749
      t752 = 0.8248731582189298D-3*t751
      t753 = t638*t126
      t756 = t753*t83*t155*t131
      t757 = 0.1975724094022595D-3*t756
      t759 = sigma*t131
      t761 = t128*t87*t314*t759
      t762 = 0.7298470036471703D-5*t761
      t765 = t367*t117*t148*t250
      t766 = 0.972010208484014D-1*t765
      t767 = t75*t136
      t769 = t767*t255*t148
      t770 = 0.1552096399632626D-1*t769
      t772 = t352*t220*t248
      t773 = 0.5733558193022059D-3*t772
      t775 = sigma*t145
      t777 = -0.2348874336782681D1*t140-0.9437856751672345D-1*t775
      t779 = t137*t117*t777
      t780 = 0.3240034028280047D-1*t779
      t782 = t137*t347*t250
      t783 = 0.3240034028280047D-1*t782
      t788 = -0.2554864360205991D1*t140*t84+0.3193580450257488D0*t155
     &*t158
      t791 = t76*t83*t788*t90
      t792 = 0.1620017014140023D-1*t791
      t794 = t102*t347*t248
      t795 = 0.286677909651103D-3*t794
      t796 = t145*t250
      t798 = t352*t694*t796
      t799 = 0.7644744257362745D-3*t798
      t801 = 1/t320/rho
      t804 = t254*t801*t90*t39
      t805 = 0.183105793556489D-3*t804
      t806 = t39*sigma
      t809 = 1/t31/t320/t138
      t812 = t102*t117*t806*t809
      t813 = 0.6764062612889941D-5*t812
      t815 = t102*t117*t775
      t816 = 0.3057897702945098D-2*t815
      t817 = t740-t744-t748+t752-t757-t762+t766-t770-t773-t780-t783
     &+t792+t795+t799-t805-t813-t816
      t820 = t169*t261
      t823 = t169*t264
      t827 = t261*t172
      t836 = t172*t264
      t841 = t479*t250*t178
      t847 = t486*t33*t178
      t852 = t175*sigma*t44*t178
      t856 = t542*t188*t272
      t875 = 0.1162159398855134D0*t739-0.2783587152200063D-1*t743
     &-0.1028277555844163D-2*t747+0.4438074192987717D-2*t751
     &-0.1063001023463839D-2*t756-0.3926803920629242D-4*t761
     &+0.5229717294848104D0*t765-0.835076145660019D-1*t769
     &-0.3084832667532488D-2*t772-0.1743239098282701D0*t779
     &-0.1743239098282701D0*t782+0.8716195491413506D-1*t791
     &+0.1542416333766244D-2*t794+0.411311022337665D-2*t798
     &-0.9851661299347411D-3*t804-0.3639276102381362D-4*t812
     &-0.164524408935066D-1*t815
      t876 = t182*t875
      t895 = -0.1675351208713011D2*t856+0.2920537730383703D1*t876
     &+0.8640090742080124D-1*t739-0.2069461866176835D-1*t743
     &-0.7644744257362745D-3*t747+0.3299492632875719D-2*t751
     &-0.7902896376090378D-3*t756-0.2919388014588681D-4*t761
     &+0.3888040833936056D0*t765-0.6208385598530505D-1*t769
     &-0.2293423277208824D-2*t772-0.1296013611312019D0*t779
     &-0.1296013611312019D0*t782+0.6480068056560093D-1*t791
     &+0.1146711638604412D-2*t794+0.3057897702945098D-2*t798
     &-0.7324231742259559D-3*t804-0.2705625045155976D-4*t812
     &-0.1223159081178039D-1*t815
      t901 = t278*t200
      t929 = -0.2513026813069517D2*t856+0.4380806595575555D1*t876
     &+0.1296013611312019D0*t739-0.3104192799265252D-1*t743
     &-0.1146711638604412D-2*t747+0.4949238949313579D-2*t751
     &-0.1185434456413557D-2*t756-0.4379082021883022D-4*t761
     &+0.5832061250904084D0*t765-0.9312578397795757D-1*t769
     &-0.3440134915813235D-2*t772-0.1944020416968028D0*t779
     &-0.1944020416968028D0*t782+0.972010208484014D-1*t791
     &+0.1720067457906618D-2*t794+0.4586846554417647D-2*t798
     &-0.1098634761338934D-2*t804-0.4058437567733965D-4*t812
     &-0.1834738621767059D-1*t815
      t942 = t64*t314
      t970 = t248*t90
      t974 = 0.1620017014140023D-1*t228*t820*t235-t744-0.1D1*t836
     &-0.1D1*t827-0.7298470036471703D-5*t686*t942*t759
     &-0.1911186064340686D-3*t679*t942*sigma*t90-0.6764062612889941D-5
     &*t228*t64*t806*t809*t90-0.183105793556489D-3*t254*t801*t229*t90
     &-0.3057897702945098D-2*t228*t297*t230-0.3822372128681373D-3*t228
     &*t262*t717+0.2359464187918086D-1*t264*t39*t145
     &+0.286677909651103D-3*t713*t299+0.286677909651103D-3*t228*t170
     &*t970
      t997 = t254*t33
      t1019 = 0.7644744257362745D-3*t352*t229*t796*t90
     &-0.5733558193022059D-3*t352*t716*t970-0.1080011342760016D-1*t679
     &*t262*t212-0.3240034028280047D-1*t137*t103*t777*t90
     &-0.4124365791094649D-3*t686*t262*t216+0.972010208484014D-1*t76
     &*t366*t87*t716*t289+0.7760481998163131D-2*t997*t170*t90
     &+0.8248731582189298D-3*t670*t612*t749-0.3240034028280047D-1*t137
     &*t224*t289-0.3240034028280047D-1*t352*t262*t220
     &-0.1975724094022595D-3*t753*t83*t155*t64*t131+t740
     &+0.1620017014140023D-1*t102*t161*t261*t235
      t1036 = -t748+t752-t757-0.5173654665442087D-2*t741*t83*t64*t90
     &-t762+0.1620017014140023D-1*t102*t788*t64*t90+t766-t770
     &+0.2160022685520031D-1*t611*t612*t289-t773+0.2545878071091769D-1
     &*t640*t264*t118*t71-t780-t783
      t1044 = t767*t33
      t1048 = t596*t264
      t1057 = t792+t795+t799-t805-t813-t816+0.6666666666666667D0*t7
     &*t264-0.1D1*t823+0.1620017014140023D-1*t102*t87*t817*t235
     &-0.1552096399632626D-1*t1044*t716*t90+0.2D1*t1048*t148-0.1D1
     &*t161*t264*t175-0.3240034028280047D-1*t352*t170*t289
      t1061 = 1.D0*t817*t64+1.D0*t820*t64-0.6172774676263781D2*t823
     &*t84*t450-0.6172774676263781D2*t827*t84*t450
     &+0.6172774676263781D2*t817*t84*t74*t179-0.6172774676263781D2
     &*t836*t84*t450-0.1234554935252756D3*t174*t841
     &+0.2956987879494076D2*t172*t114*t74*t847+0.1092333059154904D1
     &*t491*t852+t895*t100*t107-1.D0*t550*t286+t197*t306-1.D0*t901
     &*t209+2.D0*t557*t209*t285-1.D0*t201*t306*t208-1.D0*t201*t107
     &*t929+t279*t239-1.D0*t201*t239*t285+t101*(t974+t1019+t1036+t1057)
      v2rhoasigmaaa(i) = v2rhoasigmaaa(i)+0.7772672717413724D-2*rho
     &*t1061+0.1554534543482745D-1*t262-0.9595791463387565D0*t267
     &+0.1554534543482745D-1*t280-0.1554534543482745D-1*t287
     &+0.1554534543482745D-1*t307
      t1070 = t250**2
      t1072 = t367*t117*t1070
      t1073 = 0.972010208484014D-1*t1072
      t1075 = t767*t255*t250
      t1076 = 0.3104192799265252D-1*t1075
      t1078 = t352*t289*t248
      t1079 = 0.1146711638604412D-2*t1078
      t1080 = t117*t44
      t1081 = t137*t1080
      t1082 = 0.1146711638604412D-2*t1081
      t1083 = 1/t320
      t1086 = t254*t1083*t90*sigma
      t1087 = 0.2746586903347335D-3*t1086
      t1090 = 1/t31/t320/t29
      t1093 = t102*t117*t39*t1090
      t1094 = 0.5073046959667456D-5*t1093
      t1095 = t102*t1080
      t1096 = 0.1146711638604412D-2*t1095
      t1097 = t1073-t1076-t1079-t1082+t1087+t1094+t1096
      t1100 = t261**2
      t1103 = t261*t264
      t1112 = t264**2
      t1126 = t272**2
      t1127 = t542*t1126
      t1137 = t182*(0.5229717294848104D0*t1072-0.1670152291320038D0
     &*t1075-0.6169665335064975D-2*t1078-0.6169665335064975D-2*t1081
     &+0.1477749194902112D-2*t1086+0.2729457076786022D-4*t1093
     &+0.6169665335064975D-2*t1095)
      t1153 = t285**2
      t1181 = t103*t298
      t1204 = t64*t250
      t1211 = 0.972010208484014D-1*t367*t103*t1070*t90
     &+0.1620017014140023D-1*t102*t87*t1097*t235+0.1146711638604412D-2
     &*t102*t1181+0.5073046959667456D-5*t228*t229*t1090*t90
     &-0.1769598140938565D-1*t264*sigma*t44+0.2746586903347335D-3*t254
     &*t1083*t297*t90+0.1552096399632626D-1*t997*t262*t90
     &-0.4790370675386232D0*t33*t264*t84*t175-0.2D1*t1103-0.1D1*t1112
     &-0.3104192799265252D-1*t1044*t1204*t90-0.1146711638604412D-2
     &*t352*t1204*t970
      t1226 = -0.6480068056560093D-1*t352*t262*t289+0.2D1*t1048*t250
     &-0.1146711638604412D-2*t137*t1181+0.5733558193022059D-3*t228
     &*t262*t970+0.1620017014140023D-1*t102*t87*t1100*t235+t1073-t1076
     &-t1079-t1082+t1087+t1094+t1096
      s1 = 1.D0*t1097*t64+1.D0*t1100*t64-0.1234554935252756D3*t1103
     &*t84*t450+0.6172774676263781D2*t1097*t84*t74*t179
     &-0.6172774676263781D2*t1112*t84*t74*t179-0.1234554935252756D3
     &*t266*t841+0.2956987879494076D2*t264*t114*t74*t847
      t1229 = s1+0.1092333059154904D1*t265*t448*t852+(
     &-0.1675351208713011D2*t1127+0.2920537730383703D1*t1137
     &+0.3888040833936056D0*t1072-0.1241677119706101D0*t1075
     &-0.4586846554417647D-2*t1078-0.4586846554417647D-2*t1081
     &+0.1098634761338934D-2*t1086+0.2029218783866982D-4*t1093
     &+0.4586846554417647D-2*t1095)*t100*t107-2.D0*t901*t286+2.D0*t279
     &*t306+2.D0*t557*t107*t1153-2.D0*t201*t306*t285-1.D0*t201*t107*(
     &-0.2513026813069517D2*t1127+0.4380806595575555D1*t1137
     &+0.5832061250904084D0*t1072-0.1862515679559151D0*t1075
     &-0.6880269831626471D-2*t1078-0.6880269831626471D-2*t1081
     &+0.1647952142008401D-2*t1086+0.3043828175800473D-4*t1093
     &+0.6880269831626471D-2*t1095)+t101*(t1211+t1226)
      v2sigmaaa2(i) = v2sigmaaa2(i)+0.1554534543482745D-1*rho*t1229
      endif
      endif ! rho
      enddo
      
      endif ! ideriv
      return
      end

c:C_FT97subrend
