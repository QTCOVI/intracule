!-----------------------------------------------------------------------
      subroutine init_yntphi_recursive (yntphi, nmax)
!-----------------------------------------------------------------------
!     Inicializa la tabla yntphi(a,b) = Phi(2a, 2b)
!     usando relaciones de recurrencia.
!
!     yntphi(a,b) = INT_0^{2pi} sin^{2a} cos^{2b} dphi
!
!     Relaciones de recurrencia:
!       yntphi(a+1, b) = (2a+1)/(2a+2b+2) * yntphi(a, b)
!       yntphi(a, b+1) = (2b+1)/(2a+2b+2) * yntphi(a, b)
!
!     NUMERICAL NOTES (version 2, 2024):
!     --------------------------------
!     The recurrence is stable and all values are in [0, 2*pi]. No
!     overflow/underflow is possible for reasonable nmax (<= 20).
!     For large a+b the values can become very small (~1e-20 for
!     a=b=8). We add a check on the minimum value and issue a warning
!     if it goes below 1e-250.
!
!     The interface is unchanged.
!
!-----------------------------------------------------------------------
      implicit none
      real(kind=8), intent(out) :: yntphi(0:nmax,0:nmax)
      integer, intent(in) :: nmax

      integer :: a, b
      real(kind=8), parameter :: pi = 3.141592653589793238462643d0
      real(kind=8) :: minval
!-----------------------------------------------------------------------
!
!-----Caso base: a = 0, b = 0
!     Phi(0,0) = 2 * pi
!
      yntphi(0,0) = 2.0d0 * pi
!
!-----Recurrencia hacia la derecha: aumentar a
!
      do a = 0, nmax - 1
        yntphi(a+1,0) = yntphi(a,0) * dble(2*a+1)/dble(2*a+2)
      enddo
!
!-----Recurrencia hacia arriba: aumentar b
!
      do a = 0, nmax
        do b = 0, nmax - 1
          yntphi(a,b+1) = yntphi(a,b)*dble(2*b+1)/dble(2*a+2*b+2)
        enddo
      enddo
!
!-----Check on the minimum value
!
      minval = huge(1.0d0)
      do a = 0, nmax
        do b = 0, nmax
          if (yntphi(a,b) .gt. 0.0d0 .and. &
     &        yntphi(a,b) .lt. minval) then
            minval = yntphi(a,b)
          endif
        enddo
      enddo
      if (minval .lt. 1.0d-250) then
        write (*,'(a,e12.4)') &
     &    '# init_yntphi_recursive: WARNING min Phi = ',minval
        write (*,'(a)') &
     &    '# init_yntphi_recursive: close to underflow limit.'
      endif
!
      return
      end subroutine init_yntphi_recursive
