module simpson_int
use accuracy

contains

  subroutine sim_int(func, a, b, N, k1, k, ans)

  implicit none
  external func
  real(dp), intent(in)  :: a, b
  integer, intent(in)   :: N
  real(dp), intent(out) :: ans

  real(dp)              :: delta, k1, k
  real(dp), allocatable :: xi(:), sol(:), val1(:), val2(:), val3(:)
  integer               :: ii
        
  integer               :: error_alloc


  allocate(val1(N), val2(N), val3(N), xi(N),sol(N), stat=error_alloc)

  IF (error_alloc /= 0) STOP 'simpson_int: nicht genug Speicher!'
  
    delta = (b - a)/real(N)

  DO ii = 1, N
  
    xi(ii) = a + delta*(real(ii)-0.5)
    
  END DO 

  call func(N, xi + delta/2., val1, k1, k)
  call func(N, xi           , val2, k1, k)
  call func(N, xi - delta/2., val3, k1, k)
  
  sol(:) = delta/6.*(val1(:) + 4.*val2(:) + val3(:))

  !sol(:) = delta * val2(:)

  ans = sum(sol(:))

  deallocate(val1, val2, val3, xi, sol)

  end subroutine sim_int
    
end module simpson_int
