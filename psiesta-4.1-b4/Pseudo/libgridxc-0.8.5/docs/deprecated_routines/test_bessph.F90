program test_bessph_update

  ! One should get the bessph.f file from r39 and older
  use m_bessph_old, only: bessph_old => bessph
  ! This is the new bessph module
  use m_bessph
  
  integer, parameter :: dp = kind(1.d0)
  integer :: l, i
  real(dp) :: x, dx, b1, b2
  
  dx = 0.00001_dp
  do l = 0, 10
    do i = 1, 10000000
      x = (i-1) * dx
      b1 = bessph_old(l, x)
      b2 = bessph(l, x)
      
      if ( abs(b1 - b2) > 1.e-17_dp ) then
        print *, 'Different: ', l, x, b1, b2
      end if
      
    end do
  end do
  
end program test_bessph_update
