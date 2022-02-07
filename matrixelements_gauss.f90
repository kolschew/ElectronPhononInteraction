module matrixelements_gauss
use accuracy
use parameters
use simpson_int
use kvector

contains

  subroutine calc_ME(length, kvec, w_ij)
  implicit none
  
  integer, intent(in)                   :: length                
  real(dp), intent(inout)               :: kvec(length) 
  real(dp), intent(inout)               :: w_ij(length,length)
  
  real(dp)                              :: q_sol(length)
  integer                               :: idx_k, idx_kp
  integer, parameter                    :: N_int = 50
    
  q_sol = 0._dp
  
  DO idx_k = 1, length
    
    DO idx_kp = 1, length
      
        IF (idx_k .ne. idx_kp) THEN
        
          call sim_int(elmp, 0._dp, 2.*pi, N_int, kvec(idx_kp), kvec(idx_k), q_sol(idx_k))
          w_ij(idx_k, idx_kp) = mlo2 / 2 * q_sol(idx_k)
        
        
        ELSE
        
          w_ij(idx_k, idx_k) = 0._dp
          
        END IF
         !write(*,*) w_ij(idx_k, idx_kp)
    END DO !idx_kp!
  
  END DO !idx_k!

  end subroutine calc_ME

    
    subroutine elmp(N_phi, phi, val, k1, k)
    
        implicit none
        
        integer, intent(in)      :: N_phi
        real(dp), intent(inout)  :: phi(N_phi), k1, k
        real(dp), intent(inout)  :: val(N_phi)
        
        integer                  :: idx_phi
        
        
        DO idx_phi = 1, N_phi
        
            val(idx_phi) = 1. / sqrt(k1**2 + k**2 - k1*k*cos(phi(idx_phi)))
         
        END DO
     
  end subroutine elmp
  
end module matrixelements_gauss