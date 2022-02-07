module gts_gitter 
use parameters

contains

  subroutine gts(a, b, n, n0, x, w)
        use Constants
        ! "halbes" GTS-Gitter mit einem Häufungspunkt bei a
        implicit none

        integer :: j           ! Laufindex

        real :: a              ! untere Intervallgrenze
        real :: b              ! obere Intervallgrenze

        integer :: n           ! Anz. der St.-Stellen im Intervall
        integer :: n0          ! Offset fuer Index

        real :: x(1:*)         ! Stuetzstellen
        real :: w(1:*)         ! Gewichtungen

        real :: jf
        real :: pi2


        if (n == 0) return

        pi2 = 2. * Pi

        DO j = 1, n

            jf = float(j)
            x(j+n0) = (jf / (2 * n + 1.) - 1. / (pi2) * sin(pi2 * jf / (2 * n + 1.))) * 2. * (b - a) + a
            w(j+n0) = (1. - cos(pi2 * jf / (2 * n + 1))) * 2. * (b - a) / (2 * n + 1.)

        END DO


    end subroutine gts
    
end module gts_gitter
