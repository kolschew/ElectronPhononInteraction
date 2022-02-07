module kvector
use accuracy 

contains
  subroutine kcreator(kmax, length, kout)
  
    implicit none
    real(dp), intent(in)        :: kmax 
    integer, intent(in)         :: length
    real(dp), intent(inout)     :: kout(:)
    
    integer                     :: ii
    real(dp)                    :: delta 
    
    
    delta = kmax/real(length - 1) 
    
    DO ii = 1, length
    
      kout(ii) = delta * real(ii - 1./2.)

    END DO
    
  end subroutine kcreator
  
  subroutine startvec(kmax, f0, k0, sigma, length, yout)
  
    implicit none
    integer, intent(in)   :: length
    real(dp), intent(in)  :: kmax, f0, k0, sigma
    real(dp), intent(out) :: yout(:)
    
    integer               :: ii
    real(dp), allocatable :: kvec(:)
    
    
    allocate(kvec(length))
    
    call kcreator(kmax, length, kvec)
    
    DO ii = 1, length
    
      yout(ii) = f0 * exp(-(kvec(ii)-k0)**2 / (2. * sigma ** 2))
      
    END DO
    
  end subroutine startvec
  
end module kvector



