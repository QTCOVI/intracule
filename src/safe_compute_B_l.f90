      !=======================================================================
      ! safe_compute_B_l.f90 - Versión con real*8
      !=======================================================================

      subroutine safe_compute_B_l (lmax,mu,R,T,f_l,log_f_l,success)
      !-------------------------------------------------------------------------
      ! Subrutina robusta para calcular B_l con gestión de errores
      !-------------------------------------------------------------------------
      use gaussian_scaled_bessel_mod, &
          only: compute_gaussian_scaled_bessel, LOG_ZERO
      implicit none
      
      integer, intent(in) :: lmax
      real*8, intent(in) :: mu, R, T
      real*8, intent(out) :: f_l(0:lmax)
      real*8, intent(out) :: log_f_l(0:lmax)
      logical, intent(out) :: success
      
      integer :: ier
      integer :: l
      integer :: countsig
      
      ! Inicializar
      success = .false.
      countsig = 0
      f_l = 0.0d0
      log_f_l = LOG_ZERO
      
      ! Llamar a la rutina principal
      call compute_gaussian_scaled_bessel (lmax,mu,R,T,f_l,log_f_l,ier)
      
      ! Verificar si GSL funcionó correctamente
      if (ier == 0) then
        success = .true.
        return
      end if
      
      ! GSL falló, la rutina intentó la aproximación
      ! Verificar cuántos términos son significativos
      do l = 0, lmax
        if (log_f_l(l) > -700.0d0) then
          countsig = countsig + 1
        end if
      end do
      
      ! Si al menos un término es significativo, considerar éxito parcial
      if (countsig > 0) then
        success = .true.
        write(*,*) 'Warning: GSL falló (ier=', ier, &
          '), pero se usó aproximación para',countsig, 'términos'
      else
        ! Todos los términos son despreciables
        f_l = 0.0d0
        log_f_l = LOG_ZERO
        success = .false.
        write(*,*) 'Error: GSL failed, No significant terms, ier=', ier
        write(*,*) '  mu=', mu, ' R=', R, ' T=', T
      end if
      
      end subroutine safe_compute_B_l
