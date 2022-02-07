module runge_kutta
use accuracy
use parameters
use energy  
contains

 subroutine RK4(func, xstart, xend, length, nstep, init, kvec, kp, km, kp_exist, km_exist, Wp, Wm, w_ij)
 
    implicit none 
    external func
    real(dp), intent(in)  :: xstart, xend, init(:), kvec(:), km(:), kp(:), kp_exist(:), km_exist(:), Wp(:), Wm(:), w_ij(:,:)
    integer, intent(in)   :: nstep, length
    
    integer               :: ii, jj
    real(dp), allocatable :: yhelp1(:), yhelp2(:), yhelp3(:), yhelp4(:), yn(:)
    real(dp)              :: xn, delta, inner_energy, dens
    real(dp), allocatable :: ysol(:,:)
    integer, parameter    :: Nwrite = 10
    
    
    open(unit=100, status='replace',action='write',file='data/solution.dat')
  
    allocate(yhelp1(length))
    allocate(yhelp2(length))
    allocate(yhelp3(length))
    allocate(yhelp4(length))
    allocate(yn(length))
    allocate(ysol(nstep + 1, length))
    
    delta = (xend - xstart)/nstep
    xn = xstart
    yn = init
      
!     ysol(1,1) = xstart
    ysol(1,:) = init
    
    DO jj = 1, length
          
      write(100,*) xn, kvec(jj), hbar ** 2 * kvec(jj) ** 2 / (2. * meff), ysol(1, jj)
      
    END DO
    
    write(100,*)
    
    DO ii = 1, nstep
    
      IF (mod(ii, nstep / 10) == 0) write(*,*) ii, ' / ', nstep 
    
       call func(length, yn,                  xn,           yhelp1, kvec, kp, km, kp_exist, km_exist, Wp, Wm, w_ij)
       call func(length, yn + delta/2*yhelp1, xn + delta/2, yhelp2, kvec, kp, km, kp_exist, km_exist, Wp, Wm, w_ij)
       call func(length, yn + delta/2*yhelp2, xn + delta/2, yhelp3, kvec, kp, km, kp_exist, km_exist, Wp, Wm, w_ij)
       call func(length, yn + delta*yhelp3,   xn + delta,   yhelp4, kvec, kp, km, kp_exist, km_exist, Wp, Wm, w_ij)
 
       xn = xn + delta
       yn = yn + delta/6*(yhelp1 + 2*yhelp2 + 2*yhelp3 + yhelp4)
       
       ysol(ii + 1, :) = yn
       
       call calc_energy(ysol(ii + 1, :), kvec(:), length, dens, inner_energy)
       
       IF (mod(ii, Nwrite) == 0) THEN  
              
          DO jj = 1, length
              
              write(100,*) xn, kvec(jj), hbar ** 2 * kvec(jj) ** 2 / (2. * meff), ysol(ii + 1, jj), inner_energy
                
          END DO
          
          write(100,*)
          
       END IF 
       
       
       
    END DO
    
  end subroutine RK4
 
end module runge_kutta
