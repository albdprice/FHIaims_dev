!!@LICENSE
!
! *********************************************************************
! MODULE m_radfft
!
!   Public procedures provided:
! subroutine radfft              ! Radial fast Fourier transform
! 
!   Public parameters, variables, and arrays:
! none
!
!   Used module procedures:
! use m_bessph,  only: bessph    ! Spherical Bessel functions
! use m_fft_gpfa,only: fft_gpfa_ez ! Fourier transform
! use alloc,     only: de_alloc  ! Deallocation routines
! use alloc,     only: re_alloc  ! (Re)allocation routines
!
!   Used module parameters:
! use precision, only: dp        ! Double precision real kind
!
! *********************************************************************
! SUBROUTINE RADFFT( L, NR, RMAX, F, G )
! *********************************************************************
! Makes a fast Fourier transform of a radial function.
! If function f is of the form
!   f(r_vec) = F(r_mod) * Ylm(theta,phi)
! where Ylm is a spherical harmonic with l = argument L, and
! argument F contains on input the real function F(r_mod), in a uniform
! radial grid:
!   r_mod = ir * RMAX / NR,  ir = 0,1,...,NR,
! and if g is the 3-dimensional Fourier transform of f:
!   g(k_vec) = 1/(2*pi)**(3/2) *
!              Integral( d3_r_vec * exp(-i * k_vec * r_vec) * f(r_vec) )
! then g has the form
!   g(k_vec) = (-i)**L * G(k_mod) * Ylm(theta,phi)
! where argument G contains on output the real function G(k_mod) in
! a uniform radial grid:
!   k_mod = ik * k_max / NR, ik = 0,1,...,NR,  k_max = NR*pi/RMAX
! Ref: J.M.Soler notes of 16/08/95.
! *************** INPUT ***********************************************
! INTEGER L       : Angular momentum quantum number
! INTEGER NR      : Number of radial intervals.
!                   2*NR must be an acceptable number of points for the
!                   FFT routine used.
! REAL*8  RMAX    : Maximum radius
! REAL*8  F(0:NR) : Function to be tranformed, in a radial mesh
! *************** OUTPUT **********************************************
! REAL*8  G(0:NR) : Fourier transform of F (but see point 5 below)
! *************** UNITS ***********************************************
! Units of RMAX and F are arbitrary.
! Units of k_max and G are related with those of RMAX and F in the
!   obvious way (see above).
! *************** BEHAVIOUR *******************************************
! 1) F and G may be the same physical array, i.e. it is allowed:
!      CALL RADFFT( L, NR, RMAX, F, F )
! 2) It also works in the opposite direction, but then the factor
!    multiplying the output is (+i)**L. Thus, the following two calls
!      CALL RADFFT( L, NR, RMAX, F, G )
!      CALL RADFFT( L, NR, NR*PI/RMAX, G, H )
!    make H = F
! 3) If you will divide the output by q**l, truncation errors may be
!    quite large for small k's if L and NR are large. Therefore, these
!    components are calculated by direct integration rather than FFT.
!    Parameter ERRFFT is the typical truncation error in the FFT, and
!    controls which k's are integrated directly. A good value is 1e-8.
!    If you will not divide by k**l, make ERRFFT=1.e-30.
! 4) The function F is assumed to be zero at and beyond RMAX. The last
!    point F(NR) is not used to find G, except G(0) for L=0 (see 5)
! 5) Because of the 'uncertainty principle', if f(r) is strictly zero
!    for r>RMAX, then g(k) cannot be strictly zero for k>kmax.
!    Therefore G(NR), which should be exactly zero, is used (for L=0)
!    as a 'reminder' term for the integral of G beyond kmax, to ensure 
!    that F(0)=Sum[4*pi*r**2*dr*G(IR)] (this allows to recover F(0)
!    when called again in the inverse direction). Thus, the last value
!    G(NR) should be replaced by zero for any other use. NOTICE: this
!    is commented out in this version!
! *********************************************************************
! Written by J.M.Soler. August 1996.
! Work arrays handling by Rogeli Grima, ca 2009
! *********************************************************************

      MODULE m_radfft

      USE precision, only: dp        ! Double precision real kind
      USE m_bessph,  only: bessph    ! Spherical Bessel functions
      use m_fft_gpfa,only: fft_gpfa_ez     ! 1D fast Fourier transform
      USE alloc,     only: re_alloc, de_alloc

      implicit none

      PUBLIC :: radfft               ! Radial fast Fourier transform
      PUBLIC :: reset_radfft         ! Deallocates work arrays

      PRIVATE

      ! Work arrays held in module to minimize reallocations
      ! Note that we avoid "automatic" arrays, which may cause stack problems
      real(dp), pointer :: GG(:)
      real(dp), pointer :: FN(:,:)
      real(dp), pointer :: P(:,:,:)
      integer           :: MAXL = -1
      integer           :: MAXNR = -1

      CONTAINS

      SUBROUTINE RADFFT( L, NR, RMAX, F, G )

! Declare argument types and dimensions -----------------------------
      INTEGER, intent(in) :: L       ! Angular momentum of function
      INTEGER, intent(in) :: NR      ! Number of radial points
      real(dp),intent(in) :: RMAX    ! Radius of last point
      real(dp),intent(in) :: F(0:NR) ! Function to Fourier-transform
      real(dp),intent(out):: G(0:NR) ! Fourier transform of F(r)
! -------------------------------------------------------------------

! ERRFFT is the typical truncation error in the FFT routine ---------
      real(dp),   PARAMETER ::    ERRFFT = 1.0E-8_dp
! -------------------------------------------------------------------

! Internal variable types and dimensions ----------------------------
      INTEGER  ::  I, IQ, IR, JR, M, MQ, N, NQ
      real(dp) ::  C, DQ, DR, FR, PI, R, RN, Q, QMAX
!!      real(dp) ::  GG(0:2*NR), FN(2,0:2*NR), P(2,0:L,0:L)
! -------------------------------------------------------------------

!
!     Allocate local memory 
      if (MAXL.eq.-1) nullify(P)
      if (L.GT.MAXL) then
        call re_alloc( P, 1, 2, 0, L, 0, L, 'P', 'RADFFT' )
        MAXL=L
      endif
      if (MAXNR.eq.-1) nullify(FN,GG)
      if (NR.GT.MAXNR) then
        call re_alloc( FN, 1, 2, 0, 2*NR, 'FN', 'RADFFT' )
        call re_alloc( GG, 0, 2*NR, 'GG', 'RADFFT' )
        MAXNR=NR
      endif

! Find some constants -----------------------------------------------
      PI = 4.D0 * ATAN( 1.D0 )
      NQ = NR
      DR = RMAX / NR
      DQ = PI / RMAX
      QMAX = NQ * DQ
      C = DR / SQRT( 2.D0*PI )
! -------------------------------------------------------------------

! Set up a complex polynomial such that the spherical Bessel function:
!   j_l(x) = Real( Sum_n( P(n,l) * x**n ) * exp(i*x) ) / x**(l+1)
      P(1,0,0) =  0.D0
      P(2,0,0) = -1.D0
      if (l.gt.0) then
        P(1,0,1) =  0.D0
        P(2,0,1) = -1.D0
        P(1,1,1) = -1.D0
        P(2,1,1) =  0.D0
        if (l.gt.1) then
          DO M = 2,L
          DO N = 0,M
          DO I = 1,2
            P(I,N,M) = 0.D0
            IF (N .LT. M) P(I,N,M) = P(I,N,M) + (2*M-1) * P(I,N,M-1)
            IF (N .GE. 2) P(I,N,M) = P(I,N,M) - P(I,N-2,M-2)
          ENDDO
          ENDDO
          ENDDO
        endif
      endif
! -------------------------------------------------------------------

! Initialize accumulation array -------------------------------------
      DO IQ = 0,NQ
        GG(IQ) = 0.D0
      ENDDO
! -------------------------------------------------------------------

! Iterate on terms of the j_l(q*r) polynomial -----------------------
      DO N = 0,L

!       Set up function to be fast fourier transformed
        FN(1,0) = 0.D0
        FN(2,0) = 0.D0
        DO JR = 1, 2*NR-1

          IF (JR .LT. NR) THEN
            IR = JR
            R = IR * DR
            FR = F(IR)
          ELSEIF (JR .EQ. NR) THEN
            IR = JR
            R = IR * DR
            FR = 0.D0
          ELSE
            IR = 2*NR - JR
            R = - (IR * DR)
            FR = F(IR) * (-1.D0)**L
          ENDIF

!         Find  r**2 * r**n / r**(l+1)
          RN = R**(N-L+1)

          FN(1,JR) = C * FR * RN * P(1,N,L)
          FN(2,JR) = C * FR * RN * P(2,N,L)
        ENDDO

!       Perform one-dimensional complex FFT
!
!       Only the elements from 0 to 2*NR-1 of FN are used.
!       (a total of 2*NR). The fft routine will receive a one-dimensional
!       array of size 2*NR.
!
        CALL fft_gpfa_ez( FN, 2*NR, +1 )

!       Accumulate contribution
        DO IQ = 1,NQ
          Q = IQ * DQ
          GG(IQ) = ( GG(IQ) + FN(1,IQ) ) / Q
        ENDDO

      ENDDO
! -------------------------------------------------------------------

! Special case for Q=0 ---------------------------------------------
      GG(0) = 0.D0
      IF ( L .EQ. 0 ) THEN
        DO IR = 1,NR
          R = IR * DR
          GG(0) = GG(0) + R*R * F(IR)
        ENDDO
        GG(0) = GG(0) * 2.D0 * C
      ENDIF
! -------------------------------------------------------------------

! Direct integration for the smallest Q's ---------------------------
      IF (L.EQ.0) THEN
        MQ = 0
      ELSE
        MQ = NQ * ERRFFT**(1.D0/L)
      ENDIF
      DO IQ = 1,MQ
        Q = IQ * DQ
        GG(IQ) = 0.D0
        DO IR = 1,NR
          R = IR * DR
          GG(IQ) = GG(IQ) + R*R * F(IR) * BESSPH(L,Q*R)
        ENDDO
        GG(IQ) = GG(IQ) * 2.D0 * C
      ENDDO
! -------------------------------------------------------------------

! Special case for Q=QMAX -------------------------------------------
!     IF (L.EQ.0) THEN
!       GSUM = 0.D0
!       DO IQ = 1,NQ-1
!         Q = IQ * DQ
!         GSUM = GSUM + Q*Q * GG(IQ)
!       ENDDO
!       GSUM = GSUM * 4.D0 * PI * DQ
!       GG(NQ) = (2.D0*PI)**1.5D0 * F(0) - GSUM
!       GG(NQ) = GG(NQ) / (4.D0 * PI * DQ * QMAX**2)
!     ENDIF
! -------------------------------------------------------------------

! Copy from local to output array -----------------------------------
      DO IQ = 0,NQ
        G(IQ) = GG(IQ)
      ENDDO
! -------------------------------------------------------------------

      END SUBROUTINE radfft

      SUBROUTINE RESET_RADFFT( )
      call de_alloc( P, 'P', 'RADFFT' )
      call de_alloc( FN, 'FN', 'RADFFT' )
      call de_alloc( GG, 'GG', 'RADFFT' )
      MAXL  = -1
      MAXNR = -1
      END SUBROUTINE RESET_RADFFT

      END MODULE m_radfft

