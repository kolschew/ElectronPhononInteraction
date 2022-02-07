module matrixelements
use accuracy
use parameters
use kvector
use simpson_int

contains
    subroutine kpm(length, kvec, kp, km, kp_exist, km_exist)
    
      implicit none   
      integer, intent(in)                  :: length
      real(dp),  intent(inout)             :: kp(length), km(length), kp_exist(length), km_exist(length)
      real(dp), intent(inout)              :: kvec(length)
      
      integer                              :: idx_k
      
      
      DO idx_k = 1, length
      
        IF (sqrt(kvec(idx_k)**2 + (2*meff*wlo)/hbar) .lt. kvec(length)) then
        
          kp(idx_k) = sqrt(kvec(idx_k)**2 + (2*meff*wlo)/hbar)
          
          kp_exist(idx_k) = 1.
        
        ELSE 
        
          kp(idx_k) = 0.
        
          kp_exist(idx_k) = 0.
          
        END IF  
        
        IF (kvec(idx_k)**2 - (2*meff*wlo)/hbar .gt. 1e-12) then
        
            km(idx_k) = sqrt(kvec(idx_k)**2 - (2*meff*wlo)/hbar)
            
            km_exist(idx_k) = 1.
            
        ELSE
        
            km(idx_k) = 0.
            
            km_exist(idx_k) = 0.
            
        END IF
        
      END DO
        
    end subroutine kpm
    
    subroutine elements(length, kvec, kp, km, kp_exist, km_exist, Mkp, Mkm)
    
      implicit none
      integer, intent(in)       :: length
      real(dp), intent(inout)   :: kvec(length), Mkp(length), Mkm(length) 
      real(dp), intent(inout)   :: kp_exist(length), km_exist(length),kp(length), km(length)
      
      real(dp)                  :: qparp(length), qparm(length)
      real(dp)                  :: kplus, kminus, k
      integer                   :: ii
      
      
      kp = 0.
      km = 0.
      kp_exist = 0.
      km_exist = 0.
      qparp    = 0.
      qparm    = 0.
      
      call kpm(length, kvec, kp, km, kp_exist, km_exist)
      
      DO ii = 1, length
        kplus = kp(ii)
        kminus = km(ii)
        k     = kvec(ii)
        
        IF (kp_exist(ii) == 1.) THEN
        
          call sim_int(elmp, 0._dp, 2.*pi, length, kplus, k, qparp(ii))

        END IF
        
        IF (km_exist(ii) == 1.) THEN
        
          call sim_int(elmp, 0._dp, 2.*pi, length, kminus, k, qparm(ii))
            
        END IF
        
      END DO
      
      DO ii = 1, length
        
        Mkp(ii) = mlo2*pi*qparp(ii)
        Mkm(ii) = mlo2*pi*qparm(ii)
        !write(*,*) Mkp(ii)
        
      END DO

    end subroutine elements
    
    subroutine elmp(N_phi, phi, val, k1, k)
    
        implicit none
        
        integer, intent(in)      :: N_phi
        real(dp), intent(inout)  :: phi(N_phi), k1, k
        real(dp), intent(inout)  :: val(N_phi)
        
        real(dp)                 :: q_int(N_phi)
        integer                  :: idx_phi
        
          
        DO idx_phi = 1, N_phi
            
            q_int(idx_phi) = sqrt(k1**2 + k**2 - k1*k*cos(phi(idx_phi)))
            val(idx_phi) = 1. / q_int(idx_phi) * (2. / (q_int(idx_phi) * L_z) + (q_int(idx_phi) * L_z)&
            / ((q_int(idx_phi) * L_z) ** 2 + 4. * pi ** 2) + 2. * (exp(-q_int(idx_phi) * L_z - 1.)) &
            * 1. / (q_int(idx_phi) * L_z) + (q_int(idx_phi) * L_z) / ((q_int(idx_phi) * L_z) ** 2 &
            + 4. * pi ** 2))
           
        END DO
     
      end subroutine elmp
   
end module matrixelements
