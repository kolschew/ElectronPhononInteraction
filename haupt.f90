program haupt
use accuracy
use kvector
use dgl
use runge_kutta
use matrixelements
use matrixelements_gauss
use parameters

  implicit none
  integer, parameter    :: length = 1000, dgl_step = 5000
  real(dp), parameter   :: kmax = 1.0_dp, f0 = 0.1_dp, k0 = sqrt(2.*meff*11.*wlo/hbar) !11 E_LO Verschoben!
  real(dp), parameter   :: sigma = 0.005_dp, tmax = 5_dp
  real(dp), allocatable :: kvec(:), Mkp(:), Mkm(:), init(:), kp_exist(:), km(:), kp(:), km_exist(:), w_ij(:,:)

  integer               :: error_alloc
  
  
  allocate(kvec(length), Mkp(length), Mkm(length), init(length)&
  &, kp(length), km(length), kp_exist(length), km_exist(length), w_ij(length, length), stat=error_alloc)
    if (error_alloc /= 0) stop 'haupt: nicht genug Speicher!'
    
  Mkp = 0.
  Mkm = 0.
  kvec = 0.
  
  open(unit=101, status='replace',action='write',file='data/test.dat')
  
  call kcreator(kmax, length, kvec)
  
  call startvec(kmax, f0, k0, sigma, length, init)
  
  call calc_ME(length, kvec, w_ij)
  
  call elements(length, kvec, kp, km, kp_exist, km_exist, Mkp, Mkm)
  
  call rk4(func, 0.0_dp, tmax, length, dgl_step, init, kvec, kp, km, kp_exist, km_exist, Mkp, Mkm, w_ij) 

    close(101)
  
  deallocate(kvec, Mkp, Mkm)
  
end program haupt