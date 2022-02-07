module interpolation
use accuracy

contains

  subroutine lin_interpolation(N, xint, yint, x, y) 
  
  implicit none
  real(dp), intent(in)          :: x
  real(dp), intent(in)          :: xint(:), yint(:)
  integer, intent(in)           :: N
  real(dp), intent(out)         :: y
  
  integer                       :: ii
  
  
  IF ( x .lt. xint(1) .or. x .gt. xint(N) ) THEN
    
    write(*,*) 'Wert liegt nicht im Interpolationsbereich'
    
  ELSE  
  
    DO ii = 1, N
    
        IF ( x .ge. xint(ii) ) THEN
        
          IF ( x .le. xint(ii + 1) ) THEN
          
            y = yint(ii) + (yint(ii + 1) - yint(ii))/(xint(ii + 1) - xint(ii))*(x - xint(ii))
            
          END IF
          
        END IF
        
      END DO
      
  END IF
  
  end subroutine lin_interpolation
  
end module interpolation

