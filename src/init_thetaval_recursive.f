!-----------------------------------------------------------------------
      subroutine init_thetaval_recursive (thetaval, nxymax, nzmax)
!-----------------------------------------------------------------------
!     Inicializa la tabla thetaval(nxy, nz) = Theta(nxy, nz)
!     donde nxy = nx + ny, nz = nz
!
!     Relaciones de recurrencia:
!       Theta(nxy+2, nz) = (nxy+2)/(nxy+nz+3) * Theta(nxy, nz)
!       Theta(nxy, nz+2) = (nz+1)/(nxy+nz+3) * Theta(nxy, nz)
!
!     NUMERICAL NOTES (version 2, 2024):
!     --------------------------------
!     The recurrence is stable and all values are in [0, 2]. No
!     overflow/underflow is possible. However, for large nxy+nz the
!     values can become very small (~1e-20 for nxy=nz=16). We add a
!     check on the minimum value and issue a warning if it goes below
!     1e-250, which is close to the underflow limit of double precision.
!
!     The interface is unchanged.
!
!-----------------------------------------------------------------------
      implicit none
      real(kind=8), intent(out) :: thetaval(0:nxymax,0:nzmax)
      integer, intent(in) :: nxymax, nzmax

      integer :: nxy, nz
      real(kind=8) :: minval
!-----------------------------------------------------------------------
!
!-----Inicializar todo a cero
!
      thetaval = 0.0d0
!
!-----Caso base: nxy = 0
!     Theta(0, nz) = 2/(nz+1) para nz par
!
      do nz = 0, nzmax, 2
        thetaval(0, nz) = 2.0d0 / dble(nz + 1)
      enddo
!
!-----Recurrencia en nz (para nxy = 0)
!
      do nz = 0, nzmax - 2, 2
        thetaval(0, nz+2) = thetaval(0, nz)
     &    * dble(nz + 1) / dble(nz + 3)
      enddo
!
!-----Recurrencia en nxy
!
      do nz = 0, nzmax, 2
        do nxy = 2, nxymax, 2
          thetaval(nxy, nz) = thetaval(nxy-2, nz)
     &      * dble(nxy) / dble(nxy + nz + 1)
        enddo
      enddo
!
!-----Check on the minimum value
!
      minval = huge(1.0d0)
      do nz = 0, nzmax, 2
        do nxy = 0, nxymax, 2
          if (thetaval(nxy,nz) .gt. 0.0d0 .and. 
     &        thetaval(nxy,nz) .lt. minval) then
            minval = thetaval(nxy,nz)
          endif
        enddo
      enddo
      if (minval .lt. 1.0d-250) then
        write (*,'(a,e12.4)') 
     &    '# init_thetaval_recursive: WARNING min Theta = ',minval
        write (*,'(a)') 
     &    '# init_thetaval_recursive: close to underflow limit.'
      endif
 
      return
      end subroutine init_thetaval_recursive
