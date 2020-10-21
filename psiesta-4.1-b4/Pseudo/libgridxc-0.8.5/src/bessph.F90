!!@LICENSE

module m_bessph

  implicit none

  private

  public :: bessph

contains

  function bessph(L,X) result(B)

    ! RETURNS THE SPHERICAL BESSEL FUNCTION JL(X).
    ! REF: ABRAMOWITZ AND STEGUN, FORMULAS 10.1.2 AND 10.1.19
    ! WRITTEN BY J.SOLER. NOV/89.
    ! F90 conversion, by Nick Papior, feb/2018
    use precision, only: dp
    use sys, only: die

    integer, intent(in)  :: L
    real(dp), intent(in) :: X
    real(dp) :: B
    
    integer, parameter :: nterms = 100
    real(dp), parameter :: zero = 0.0_dp, one = 1.0_dp, tiny=1.0e-15_dp

    integer :: i, n
    real(dp) :: switch, term, x2, y, fnm1, fn, fnp1
    character(len=128):: msg

    if ( abs(x) < max(1, 2*L - 1) ) then
      ! Use power series
      TERM = ONE
      do I = 1 , L
        TERM = TERM * X/(2*I+1)
      end do
      X2 = X * X
      FN = ZERO
      do I = 1 , NTERMS
        FN = FN + TERM
        TERM = -TERM * X2 / (4._dp * I * (I + L + 0.5_dp))
        if ( abs(TERM) < tiny ) exit
      end do
      if ( abs(TERM) > tiny ) then
        write(msg,*) 'BESSPH: SERIES HAS NOT CONVERGED. L,X=',L,X
        call die(trim(msg))
      end if
      B = FN
      
    else
      
      ! Use explicit expressions or recurrence relation
      select case ( L )
      case ( 0 )
        
        B = sin(X) / X
        
      case ( 1 )
        
        B = (sin(X) / X - cos(X) ) / X
        
      case default
        
        Y = ONE / X
        fnm1 = sin(X) * Y
        fn = (fnm1 - cos(X)) * Y
        DO N = 1 , L-1
          FNP1 = (2 * N + 1) * Y * fn - fnm1
          FNM1 = FN
          FN = FNP1
        end do
        B = fn
        
      end select
      
    end if
            
  end function bessph
  
end module m_bessph

