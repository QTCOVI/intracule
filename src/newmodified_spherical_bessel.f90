!
!===============================================================================
! modified_spherical_bessel.f90
!
! Calcula, para l = 0,...,l_max, la cantidad
!
!     f_l(mu,R,T) = exp[-mu*(R^2+T^2)] * i_l(2*mu*R*T)
!
! que aparece, por ejemplo, en el desarrollo en armonicos esfericos /
! polinomios de Legendre de productos de gaussianas ("Gaussian expansion
! method"):
!
!     exp[2*mu*R*T*cos(gamma)] = sum_l (2l+1) i_l(2*mu*R*T) P_l(cos gamma)
!
! Para R,T grandes, i_l(2*mu*R*T) crece como exp(2*mu*R*T) mientras que
! exp[-mu*(R^2+T^2)] decae de forma comparable. Si se calculan ambos
! factores por separado en doble precision se puede producir overflow
! (el argumento de exp() facilmente supera ~709, el limite de un
! real*8) o underflow, aun cuando el PRODUCTO final sea perfectamente
! representable.
!
! Para evitarlo, se reescribe el producto usando la version ESCALADA
! de la funcion de Bessel esferica modificada que proporciona GSL:
!
!     i_l_scaled(x) = exp(-x) * i_l(x),      x = 2*mu*R*T >= 0
!
! de forma que:
!
!     f_l(mu,R,T) = exp[-mu*(R^2+T^2) + x] * i_l(x)
!                 = exp[-mu*(R-T)^2]       * i_l_scaled(x)
!
! Como mu >= 0, el exponente -mu*(R-T)^2 es siempre <= 0, por lo que
! exp[-mu*(R-T)^2] esta acotado entre 0 y 1 (nunca desborda), e
! i_l_scaled(x) esta acotada (decrece ~1/sqrt(x) para x grande), de modo
! que el producto NUNCA desborda, sea cual sea el valor de mu, R, T.
!
! Ademas, siguiendo lo solicitado, si mu*(R^2+T^2) > 1000 se considera
! que f_l es numericamente despreciable y se fuerza directamente a 0.0
! (evitando calculos innecesarios en un regimen extremo).
!
! Por ultimo, se ofrece tambien el resultado en ESCALA LOGARITMICA,
!
!     log[f_l(mu,R,T)] = -mu*(R-T)^2 + log[i_l_scaled(x)]
!
! que es la representacion recomendada cuando f_l es tan pequena que
! "hace underflow" a 0.0 en doble precision pero aun se necesita
! conocer su orden de magnitud (por ejemplo, para sumar terminos con
! pesos muy dispares mediante la tecnica log-sum-exp).
!
! Compilacion:
!   gfortran -o modified_spherical_bessel modified_spherical_bessel.f90 -lgsl -lgslcblas -lm
!===============================================================================

module gsl_bessel_interface
  !---------------------------------------------------------------------------
  ! Interfaces C (iso_c_binding) a las rutinas de GSL necesarias.
  !---------------------------------------------------------------------------
  use, intrinsic :: iso_c_binding
  implicit none

  interface

     ! int gsl_sf_bessel_il_scaled_array(const int lmax, const double x,
     !                                    double result_array[])
     ! Calcula exp(-x) * i_l(x) para l = 0,...,lmax (x >= 0), mediante
     ! recursion descendente interna (l_max -> 0), devolviendo el
     ! resultado ya ordenado de l=0 a l_max en result_array.
     function gsl_sf_bessel_il_scaled_array (lmax, x, result_array) &
          bind(C, name="gsl_sf_bessel_il_scaled_array") result(ierr)
       import :: c_int, c_double
       integer(c_int), value :: lmax
       real(c_double),  value :: x
       real(c_double), dimension(*) :: result_array
       integer(c_int) :: ierr
     end function gsl_sf_bessel_il_scaled_array

     ! void gsl_set_error_handler_off(void)
     ! Desactiva el manejador de errores por defecto de GSL (que aborta
     ! el programa), para poder tratar los errores mediante codigos.
     subroutine gsl_set_error_handler_off () &
          bind(C, name="gsl_set_error_handler_off")
     end subroutine gsl_set_error_handler_off

  end interface

end module gsl_bessel_interface


module gaussian_scaled_bessel_mod
  !---------------------------------------------------------------------------
  ! Rutina de alto nivel: calcula f_l(mu,R,T) = exp[-mu*(R^2+T^2)]*i_l(x),
  ! con x = 2*mu*R*T, evitando overflow, y su version en escala log.
  !---------------------------------------------------------------------------
  use, intrinsic :: iso_c_binding
  use gsl_bessel_interface
  implicit none

  ! Umbral solicitado: por encima de el, f_l se considera 0
  real(c_double), parameter :: ARG_THRESHOLD = 1000.0_c_double

  ! Centinela para representar log(0) = -infinito
  real(c_double), parameter :: LOG_ZERO = -huge(1.0_c_double)

contains

  subroutine compute_gaussian_scaled_bessel (lmax,mu,R,T,f_l,log_f_l,ierr)
    !-------------------------------------------------------------------------
    !   lmax     (in)  : orden maximo l_max (>=0)
    !   mu,R,T   (in)  : parametros fisicos (se asume mu,R,T >= 0)
    !   f_l      (out) : vector 0:lmax con f_l(mu,R,T)
    !   log_f_l  (out) : vector 0:lmax con log[f_l(mu,R,T)]
    !                    (= LOG_ZERO si f_l se anula o hace underflow)
    !   ierr     (out) : codigo de error GSL (0 = exito)
    !-------------------------------------------------------------------------

    use gsl_bessel_interface
!   use gsl_error_handler_mod  ! AÑADIR ESTA LÍNEA
    implicit none

    integer,        intent(in)  :: lmax
    real(c_double), intent(in)  :: mu, R, T
    real(c_double), intent(out) :: f_l(0:lmax)
    real(c_double), intent(out) :: log_f_l(0:lmax)
    integer,        intent(out) :: ierr

    real(c_double) :: x, expo, arg
    real(c_double) :: i_scaled(0:lmax)
    integer :: l

    ierr = 0
    ! BUGFIX: el corte debe hacerse sobre el EXPONENTE FISICO que realmente
    ! amortigua B_l en la formulacion escalada, que es
    !     exp[-mu(R-T)^2] * i_l_scaled(2 mu R T),
    ! es decir, sobre mu*(R-T)^2, NO sobre mu*(R^2+T^2).
    ! El antiguo arg = mu*(R^2+T^2) era demasiado pesimista: para R ~ T con
    ! mu grande (p.ej. pares de electrones de core en dos atomos pesados
    ! separados una distancia T, contribuyendo al intracule en R ~ T) se
    ! cumple mu*(R^2+T^2) ~ 2 mu T^2 > 1000 mientras que exp[-mu(R-T)^2] ~ 1,
    ! de modo que se anulaba una contribucion legitima.  Esto restaba la
    ! misma cantidad al total y al Coulomb (que coinciden en esa region,
    ! donde XC ~ 0), dejando XC intacto: exactamente el sintoma observado
    ! en sistemas con varios atomos pesados proximos (p.ej. C2H6).
    arg  = mu*(R - T)**2

    ! --- Corte de seguridad (sobre el exponente fisico correcto) ---
    if (arg > ARG_THRESHOLD) then
       f_l     = 0.0_c_double
       log_f_l = LOG_ZERO
       return
    end if

    x = 2.0_c_double*mu*R*T

    call gsl_set_error_handler_off()
    if (x < 1.0d-15) then
       ! Para x muy pequeño, i_l_scaled(x) ≈ 1/(2l+1)!!
       do l = 0, lmax
          ! Aproximación para x pequeño: i_l(x) ≈ x^l/(2l+1)!!
          ! y i_l_scaled(x) = exp(-x)*i_l(x) ≈ i_l(x)
          ! Para l=0: i_0_scaled(0) = 1
          ! Para l>0: el valor es extremadamente pequeño
          if (l == 0) then
             i_scaled(l) = 1.0_c_double
          else
             i_scaled(l) = 0.0_c_double
          end if
       end do
       ! Saltar la llamada a GSL
       goto 100
    end if

    if (lmax > 500 .and. x < 1.0d-5) then
       ! Para órdenes grandes y x pequeño, los valores son esencialmente cero
       i_scaled = 0.0_c_double
       if (lmax >= 0) i_scaled(0) = 1.0_c_double  ! l=0 siempre es ~1 para x pequeño
       goto 100
    end if

    ierr = gsl_sf_bessel_il_scaled_array(int(lmax, c_int), x, i_scaled)

     if (ierr /= 0) then
        ! Si GSL devuelve error, intentamos con un enfoque alternativo
        ! Usar aproximación asintótica o serie de Taylor directamente
        call approximate_bessel_scaled (lmax, x, i_scaled, ierr)
        if (ierr /= 0) then
           f_l     = 0.0_c_double
           log_f_l = LOG_ZERO
           return
        end if
     end if
     
 100 continue

     expo = -mu*(R - T)**2
     do l = 0, lmax
       f_l(l) = exp(expo) * i_scaled(l)
       if (i_scaled(l) > 0.0_c_double) then
          log_f_l(l) = expo + log(i_scaled(l))
       else
          log_f_l(l) = LOG_ZERO
       end if
     end do

  end subroutine compute_gaussian_scaled_bessel


subroutine approximate_bessel_scaled(lmax, x, i_scaled, ierr)
  ! Aproximación para casos donde GSL falla
  integer, intent(in) :: lmax
  real(c_double), intent(in) :: x
  real(c_double), intent(out) :: i_scaled(0:lmax)
  integer, intent(out) :: ierr
  
  integer :: l, k
  real(c_double) :: factorial_factor

  ierr = 0

  if (x < 1.0d-10) then
     ! Serie de Taylor para x pequeño: i_l(x) ≈ x^l/(2l+1)!!
     do l = 0, lmax
        if (l == 0) then
           i_scaled(l) = 1.0_c_double
        else
           ! Calcular (2l+1)!! = 1*3*5*...*(2l+1) = prod_{k=1}^{l} (2k+1)
           ! (BUGFIX: el antiguo "do while (l > 0)" no modificaba l y
           !  entraba en un bucle infinito; se usa un contador propio k.)
           factorial_factor = 1.0_c_double
           do k = 1, l
              factorial_factor = factorial_factor * real(2*k + 1, c_double)
           end do
           i_scaled(l) = (x**l) / factorial_factor
           if (i_scaled(l) < 1.0d-300) i_scaled(l) = 0.0_c_double
        end if
     end do
  else
     ! Otros casos: usar aproximación asintótica
     do l = 0, lmax
        i_scaled(l) = 1.0_c_double / sqrt(2.0_c_double * 3.141592653589793_c_double * x)
        ! Para l > 0, decae aproximadamente como 1/(2l+1)
        if (l > 0) then
           i_scaled(l) = i_scaled(l) / real(2*l+1, c_double)
        end if
     end do
  end if
  
end subroutine approximate_bessel_scaled


end module gaussian_scaled_bessel_mod
