!!@LICENSE

module cellSubs

  use precision, only: dp   ! Double precision real kind

  implicit none

  private

  public :: reclat ! Reciprocal unit cell vectors
  public :: volcel  ! Volume of the unit cell

contains

  subroutine reclat(A,B,IOPT)
    
    !  CALCULATES RECIPROCAL LATTICE VECTORS. THEIR PRODUCT WITH DIRECT
    !  LATTICE VECTORS IS 1 IF IOPT=0 OR 2*PI IF IOPT=1
    real(dp), intent(in) :: A(3,3)
    real(dp), intent(out) :: B(3,3)
    integer, intent(in) :: iopt

    real(dp), parameter :: PI = 3.14159265358979323846264338327950288419716939937510_dp

    integer :: i
    real(dp):: c , ci
    
    B(1,1) = A(2,2)*A(3,3)-A(3,2)*A(2,3)
    B(2,1) = A(3,2)*A(1,3)-A(1,2)*A(3,3)
    B(3,1) = A(1,2)*A(2,3)-A(2,2)*A(1,3)
    B(1,2) = A(2,3)*A(3,1)-A(3,3)*A(2,1)
    B(2,2) = A(3,3)*A(1,1)-A(1,3)*A(3,1)
    B(3,2) = A(1,3)*A(2,1)-A(2,3)*A(1,1)
    B(1,3) = A(2,1)*A(3,2)-A(3,1)*A(2,2)
    B(2,3) = A(3,1)*A(1,2)-A(1,1)*A(3,2)
    B(3,3) = A(1,1)*A(2,2)-A(2,1)*A(1,2)
    C=1.D0
    if ( IOPT == 1 ) C = 2.D0*PI
    do I = 1 , 3
      CI = C/(A(1,I)*B(1,I)+A(2,I)*B(2,I)+A(3,I)*B(3,I))
      B(1,I) = B(1,I)*CI
      B(2,I) = B(2,I)*CI
      B(3,I) = B(3,I)*CI
    end do
    
  end subroutine reclat

  !--------------------------------------------------------------------

  function volcel( C ) result(V)
    !  CALCULATES THE VOLUME OF THE UNIT CELL
    real(dp), intent(in) :: C(3,3)
    real(dp) :: V
    
    V = ( C(2,1)*C(3,2) - C(3,1)*C(2,2) ) * C(1,3) + &
        ( C(3,1)*C(1,2) - C(1,1)*C(3,2) ) * C(2,3) + &
        ( C(1,1)*C(2,2) - C(2,1)*C(1,2) ) * C(3,3)
    V = ABS( V )
    
  end function volcel

end module cellsubs

