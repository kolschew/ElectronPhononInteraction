module parameters
use accuracy

  implicit none
  !Parameter sind in den Einheiten ps, meV, nm, K und pA 
  
  
  real(dp), parameter :: pi = 3.14159265359
  real(dp), parameter :: temp = 77.
  real(dp), parameter :: kb = 0.08617
  real(dp), parameter :: echarge = 1.602e5
  real(dp), parameter :: mfree = 5.686e-3
  real(dp), parameter :: meff = 0.064*mfree
  real(dp), parameter :: hbar = 0.6582
  real(dp), parameter :: wlo = 54.71
  real(dp), parameter :: beta = 1./(kb * temp)
  real(dp), parameter :: nlo = 1./(exp(beta*hbar*wlo)-1.)
  real(dp), parameter :: eps = 12.9
  real(dp), parameter :: eps0 = 1.418e6
  real(dp), parameter :: epsinf = 10.9
  real(dp), parameter :: mlo2 = (hbar*wlo*echarge**2)/(2*eps0)*(1./epsinf - 1./eps)  
  real(dp), parameter :: gamma = 8_dp
  real(dp), parameter :: L_z = 6.5
  
end module parameters