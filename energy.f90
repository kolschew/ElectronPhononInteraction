module energy
use accuracy
use parameters
use brent
contains
  subroutine calc_energy(y_in, k_vec, length, dens, inner_energy)
    implicit none
    real(dp), intent(in)        :: y_in(:), k_vec(:)    !y- und k-Werte, die extern aus jedem Zeitschritt kommen
    integer, intent(in)         :: length               
    real(dp), intent(out)       :: dens, inner_energy   !Elektronendichte und innere Energie aus jedem Zeitschritt, für das GLS     
    
    real(dp)                    :: const, energy, energy_dist, k_dist
    integer                     :: ii
    
    k_dist = k_vec(2) - k_vec(1)                        !Abstand zweier k-Punkte, äquidistant
    const = meff / (pi * hbar ** 2)                     
    energy_dist = (hbar ** 2 * k_dist **2) / (2 * meff) 
     
    dens = 0.
    inner_energy = 0.
    
    DO ii = 1, length
    
      energy = (hbar ** 2 * k_vec(ii) **2) / (2 * meff)
      dens = dens + const * energy_dist * y_in(ii) 
      inner_energy = inner_energy + const * energy_dist * energy * y_in(ii)
      
    END DO   
    
  end subroutine calc_energy

  subroutine fcn(n, x, f_vec, i_flag)
    implicit none 
    integer, intent(in)         :: n, i_flag
    real(dp), intent(inout)     :: x(n), f_vec(n)
    
    f_vec(1) = exp(x(2)) - x(1)
    f_vec(2) = sin(x(1)) + x(2) ** 2
    
  end subroutine fcn    
  
  call brent1(fcn, 2, [1, 1], f_vec, 1e-6, info, WA, lwa)
end module energy