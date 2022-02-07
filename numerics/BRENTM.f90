!************************************************************************
!* This program solves a nonlinear system of two variables:             *
!*                                                                      *
!*             f(x,y)  =  0                                             *
!*             g(x,y)  =  0                                             *
!*                                                                      *
!* using the Non_linear_system subroutine.                              *
!* -------------------------------------------------------------------- *
!* REFERENCE: "Mathematiques en Turbo-Pascal part 1, By M. Ducamp and   *
!*             A. Reverchon, Editions EYROLLES, Paris, 1987" [BIBLI 03] *
!*                                                                      *
!*                                   F90 version by J-P Moreau, Paris   *
!*                                          (www.jpmoreau.fr)           *
!* -------------------------------------------------------------------- *
!* SAMPLE RUN:                                                          *
!*                  (Solve the non linear system                        *
!*                  x^2 + x + y^2 - 2          = 0                      *
!*                  x^2 + y - y^2 - 1 + log(x) = 0 )                    *
!*                                                                      *
!* Desired precision........... : 1e-10                                 *
!* Maximal number of iterations : 30                                    *
!* Componants of starting vector:  1  1                                 *
!*                                                                      *
!* Results in file nlinsyst.lst:                                        *
!* (Error status: 0)                                                    *
!*                                                                      *
!* ----------------------------------------------                       *
!*  Method for nonlinear systems of two equations                       *
!* ----------------------------------------------                       *
!*                                                                      *  
!*  System to be solved:                                                *
!*  x^2 + x + y^2 - 2          = 0                                      *
!*  x^2 + y - y^2 - 1 + log(x) = 0                                      *
!*                                                                      *  
!*  Starting vector:                                                    *
!*    1.000000      1.000000                                            *
!*                                                                      *  
!*  Error bound =  0.100000E-09                                         * 
!*  Maximal number of iterations =   30                                 *
!*                                                                      *  
!*  Solution vector:                                                    *
!*    0.915554      0.496191                                            *
!*                                                                      *  
!*  Number of iterations:    7                                          *
!* ----------------------------------------------                       *
!*                                                                      *  
!*  End of file nlinsyst.lst.                                           *
!*                                                                      *
!************************************************************************
PROGRAM TEST_NLINSYST

integer  n,maxit,iter   ! n     size of system
                        ! maxit maximal number of iterations
                        ! iter  number of iterations performed

real*8   eps,x0,y0,x,y  ! eps   desired accuracy
                        ! x0,y0 starting point
                        ! x,y   approximate solution

integer  error          !error code:  0 = OK
                        !             1 = error in evluating f(x,y) or g(x,y)
                        !             2 = singular system

! begin main program
! -------------------- read input --------------------------------
  n=2    !size of system
  print *,' '
  write(*,"(' Desired precision..: ')",advance='no') 
  read  *, eps
  write(*,"(' Maximal number of iterations : ')",advance='no') 
  read  *, maxit
  print *, 'Componants of starting point (x0,y0):'
  write(*,"('    x0 = ')",advance='no') 
  read  *, x0
  write(*,"('    y0 = ')",advance='no') 
  read  *, y0
! ------------ print input for checking purposes -----------------
  OPEN(UNIT=1,FILE='nlinsyst.lst',STATUS='UNKNOWN')
  write(1,*) '----------------------------------------------'
  write(1,*) ' Method for nonlinear systems of two equations'
  write(1,*) '----------------------------------------------'
  write(1,*) ' '
  write(1,*) ' System to be solved:'
  write(1,*) ' x^2 + x + y^2 - 2          = 0'
  write(1,*) ' x^2 + y - y^2 - 1 + log(x) = 0'
  write(1,*) ' '
  write(1,*) ' Starting vector:'
  write(1,50)  x0, y0
  write(1,*) ' '
  write(1,60)  eps
  write(1,70)  maxit

! ------------ solve nonlinear system ----------------------------
  call Non_linear_system(n,x0,y0,eps,maxit,x,y,iter,error)

! --------------------- print solution ---------------------------
  write(1,*) ' '
  write(1,*) ' Solution vector:'
  write(1,50)  x, y
  write(1,*) ' '
  write(1,80)  iter
  write(1,*) '----------------------------------------------'
  write(1,*) ' '
  write(1,*) ' End of file nlinsyst.lst.'
  write(1,*) ' '
  close(unit=1)

  print *,' '
  print *,'Results in file nlinsyst.lst.'
  print *,'Error status:',error
  print *,' '

  stop

50 format(f12.6,'  ',f12.6)
60 format('  Error bound = ',e13.6)
70 format('  Maximal number of iterations = ',i4)
80 format('  Number of iterations: ',i4)

end !of main program


real*8 Function f(x,y)
  real*8 x,y
  f=x*x+x+y*y-2.d0
end

real*8 Function g(x,y,rc)
  real*8 x,y
  integer rc
  rc=0
  if (x>0.d0) then
    g=x*x+y-y*y-1.d0+DLOG(x)
  else
    rc=1; g=0.d0
  end if
end


Subroutine Non_linear_system(n,x0,y0,eps,maxit,x,y,iter,ierror)
  integer irc,n,maxit,iter
  real*8  x0,y0,eps,x,y
  real*8  a,b,c,d,t,xm,xn,p,q
  real*8  f, g
  h = 0.01d0
  ierror=1; iter=0
  x=x0; y=y0
100 iter=iter+1
  if (iter > maxit) return
  a=f(x+h,y)
  b=g(x+h,y,irc); if (irc.ne.0) return
  a=(a-f(x-h,y))/2.d0/h
  b=(b-g(x-h,y,irc))/2.d0/h; if (irc.ne.0) return
  c=f(x,y+h)
  d=g(x,y+h,irc); if (irc.ne.0) return  
  c=(c-f(x,y-h))/2.d0/h
  d=(d-g(x,y-h,irc))/2.0/h; if (irc.ne.0) return
  t=a*d-b*c
  if (dabs(t)<1.d-12) then
    ierror=2
    return
  end if
  xm=f(x,y)
  xn=g(x,y,irc); if (irc.ne.0) return
  p=(xm*d-xn*c)/t
  q=(xn*a-xm*b)/t
  x=x-p; y=y-q
  if (dabs(p)+dabs(q) > eps) goto 100
  ierror=0
  return
end

!End of file nlinsyst.f90
