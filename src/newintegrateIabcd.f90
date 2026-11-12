!===============================================================================
! MODULE: intracule_integral
! 
! DESCRIPTION: Provides a routine to compute the radial integral
!              I(n,l,mu,T) = ∫_0^∞ R^n B_l(R) dR
!
! MATHEMATICAL BACKGROUND:
!   The general integral is given by (CORRECTED 2026-09: GR 6.631.1-type
!   derivation, sympy+mpmath verified 40 dps, see
!   ground-truth/derive_iexact_correct.py):
!   
!   I(n,l,mu,T) = √π/4 e^{-μT²} μ^{(l-n-1)/2} T^l
!                 Γ((n+l+1)/2)/Γ(l+3/2)
!                 _1F_1((n+l+1)/2; l+3/2; μT²)
!   
!   The OLD formula (ported from integrateIabcd.f90 v1) had wrong exponents:
!   √π/2 · 2^{(l-n-2)/2} · μ^{(l-n-3)/2} — off by factor 2^{(n-l)/2} mu.
!   It contradicted the T=0 branch of the same function. Fixed 2026-09-10.
!   jaqg 2026-09-09: bug 1 of HANDOFF_evelio_code_bugs.md — wrong 2-power
!   and mu-power exponents in the closed form; net correction 2^{(n-l)/2}·mu.
!   
!   This formula is valid for ALL n >= 0, l >= 0, and T > 0.
!   There is NO selection rule n = l+1.
!   
!   For T = 0, the integral reduces to:
!     I(n,l,0) = 1/2 μ^{-(n+1)/2} Γ((n+1)/2)  for l = 0
!     I(n,l,0) = 0                            for l > 0
!   
!   Corollary (used by mod_analytic_norm): for n = l+2, a = b = l+3/2, so
!   _1F_1(a;a;x) = e^x cancels e^{-x} and
!   I(l+2,l,mu,T)/T^l = √π/(4 mu^{3/2}), independent of l and T.
!===============================================================================

module intracule_integral
  use, intrinsic :: iso_c_binding
  implicit none
  private

  ! Public interface
  public :: radial_integral, radial_integral_log
  public :: disable_gsl_error_handler
  public :: get_status_message
! ! Public interface
! public :: radial_integral, radial_integral_log
! public :: disable_gsl_error_handler
! public :: get_status_message

  ! Mathematical constants
  real(8), parameter :: PI = 3.141592653589793238462643383279502884d0
  real(8), parameter :: LOG_HUGE = log(huge(1.0d0))
  real(8), parameter :: LOG_TINY = log(tiny(1.0d0))

  ! Status codes
  integer, parameter :: STATUS_SUCCESS = 0
  integer, parameter :: STATUS_UNDERFLOW = 1
  integer, parameter :: STATUS_OVERFLOW = 2
  integer, parameter :: STATUS_ZERO_BY_RULE = 3
  integer, parameter :: STATUS_GSL_ERROR = 4
  integer, parameter :: STATUS_INVALID_PARAM = 5

  ! GSL function prototypes (C bindings)
  interface
     function gsl_sf_hyperg_1F1(a, b, x) bind(c, name='gsl_sf_hyperg_1F1')
       use, intrinsic :: iso_c_binding
       real(c_double), value :: a, b, x
       real(c_double) :: gsl_sf_hyperg_1F1
     end function gsl_sf_hyperg_1F1

     function gsl_sf_hyperg_1F1_e(a, b, x, result) bind(c, name='gsl_sf_hyperg_1F1_e')
       use, intrinsic :: iso_c_binding
       real(c_double), value :: a, b, x
       type(c_ptr), value :: result
       integer(c_int) :: gsl_sf_hyperg_1F1_e
     end function gsl_sf_hyperg_1F1_e

     function gsl_sf_gamma(x) bind(c, name='gsl_sf_gamma')
       use, intrinsic :: iso_c_binding
       real(c_double), value :: x
       real(c_double) :: gsl_sf_gamma
     end function gsl_sf_gamma

     function gsl_sf_lngamma(x) bind(c, name='gsl_sf_lngamma')
       use, intrinsic :: iso_c_binding
       real(c_double), value :: x
       real(c_double) :: gsl_sf_lngamma
     end function gsl_sf_lngamma

     subroutine gsl_set_error_handler_off() bind(c, name='gsl_set_error_handler_off')
     end subroutine gsl_set_error_handler_off
  end interface

  ! GSL result structure
  type, bind(c) :: gsl_sf_result
     real(c_double) :: val
     real(c_double) :: err
  end type gsl_sf_result

contains

  !============================================================================
  ! SUBROUTINE: disable_gsl_error_handler
  !============================================================================
  subroutine disable_gsl_error_handler()
    call gsl_set_error_handler_off()
  end subroutine disable_gsl_error_handler

  !============================================================================
  ! SUBROUTINE: get_status_message
  !============================================================================
  subroutine get_status_message(status, msg)
    integer, intent(in) :: status
    character(len=*), intent(out) :: msg
    
    select case(status)
    case(0)  ! STATUS_SUCCESS
       msg = 'SUCCESS'
    case(1)  ! STATUS_UNDERFLOW
       msg = 'UNDERFLOW (zero, but log valid)'
    case(2)  ! STATUS_OVERFLOW
       msg = 'OVERFLOW'
    case(3)  ! STATUS_ZERO_BY_RULE
       msg = 'ZERO (T=0,l>0)'
    case(4)  ! STATUS_GSL_ERROR
       msg = 'GSL_ERROR (used asymptotic)'
    case(5)  ! STATUS_INVALID_PARAM
       msg = 'INVALID_PARAMETER'
    case default
       msg = 'UNKNOWN'
    end select
  end subroutine get_status_message



  subroutine radial_integral(n, l, mu, T, integral, status)
    integer, intent(in) :: n, l
    real(8), intent(in) :: mu, T
    real(8), intent(out) :: integral
    integer, intent(out) :: status

    real(8) :: a, b, x, hyperg
    real(8) :: log_prefactor, log_result
    real(8) :: log_T_factor, log_2_factor, log_mu_factor
    real(8) :: log_gamma_num, log_gamma_den
    real(8) :: d1, d2 ! NEW
    type(gsl_sf_result), target :: gsl_result
    integer(c_int) :: gsl_status
    logical :: use_asymptotic

    integral = 0.0d0
    status = STATUS_SUCCESS
    use_asymptotic = .false.

    ! Input validation
    if (n < 0 .or. l < 0) then
       integral = 0.0d0
       status = STATUS_INVALID_PARAM
       return
    end if

    if (mu <= 0.0d0) then
       integral = 0.0d0
       status = STATUS_INVALID_PARAM
       return
    end if

    ! T = 0 case
    if (T <= 0.0d0) then
       if (l == 0) then
          integral = 0.5d0 * (mu ** (-(n+1)/2.0d0)) * gsl_sf_gamma((n+1)/2.0d0)
          status = STATUS_SUCCESS
       else
          integral = 0.0d0
          status = STATUS_ZERO_BY_RULE
       end if
       return
    end if

    ! General case: T > 0
    a = 0.5d0 * dble(n + l + 1)
    b = dble(l) + 1.5d0
    x = mu * T * T

    !-------------------------------------------------------------------------
    ! Compute log of prefactor (excluding hypergeometric)
    !-------------------------------------------------------------------------
    ! I = √π/4 * e^{-μT²} * μ^{(l-n-1)/2} * T^l
    !     * Γ((n+l+1)/2) / Γ(l+3/2) * _1F_1((n+l+1)/2; l+3/2; μT²)
    ! (CORRECTED 2026-09-10: was √π/2 · 2^{(l-n-2)/2} · μ^{(l-n-3)/2} — wrong
    !  exponents, contradicted T=0 branch; verified ground-truth/
    !  derive_iexact_correct.py)
    ! jaqg 2026-09-09: bug 1 — wrong exponents, off by 2^{(n-l)/2}·mu.
    !-------------------------------------------------------------------------
    log_prefactor = 0.5d0 * log(PI) - 2.0d0 * log(2.0d0)
    log_prefactor = log_prefactor - mu * T * T
    log_2_factor = 0.0d0
    log_mu_factor = 0.5d0 * dble(l - n - 1) * log(mu)
    log_T_factor = dble(l) * log(T)
    log_gamma_num = gsl_sf_lngamma(a)
    log_gamma_den = gsl_sf_lngamma(b)

    log_prefactor = log_prefactor + log_2_factor + log_mu_factor + &
                    log_T_factor + log_gamma_num - log_gamma_den

    !-------------------------------------------------------------------------
    ! Evaluate _1F_1
    !-------------------------------------------------------------------------
    ! For large x, use asymptotic approximation directly
!   if (x > 200.0d0) then   ! OLD
    if (x > 300.0d0) then   ! NEW
       use_asymptotic = .true.
    else
       ! Try GSL
       gsl_status = gsl_sf_hyperg_1F1_e(a, b, x, c_loc(gsl_result))
       if (gsl_status == 0) then
          hyperg = gsl_result%val
          if (hyperg > 0.0d0 .and. hyperg == hyperg) then
             log_result = log_prefactor + log(hyperg)
             
             if (log_result > LOG_HUGE) then
                integral = huge(1.0d0)
                status = STATUS_OVERFLOW
                return
             else if (log_result < LOG_TINY) then
                integral = 0.0d0
                status = STATUS_UNDERFLOW
                return
             else
                integral = exp(log_result)
                status = STATUS_SUCCESS
                return
             end if
          end if
       end if
       use_asymptotic = .true.
    end if

    !-------------------------------------------------------------------------
    ! Asymptotic approximation: _1F_1(a;b;x) ~ Γ(b)/Γ(a) e^x x^{a-b}
    ! This is the leading term of the asymptotic expansion.
    !-------------------------------------------------------------------------
    if (use_asymptotic) then
       ! Add correction for better accuracy at moderate x
       ! _1F_1(a;b;x) ~ Γ(b)/Γ(a) e^x x^{a-b} [1 + (a-1)(a-b)/x + ...]

       log_result = log_prefactor + gsl_sf_lngamma(b) - gsl_sf_lngamma(a) + &
                    x + (a - b) * log(x)
       
       if (x > 10.0d0) then                                                   ! NEW
         d1 = (a - 1.0d0) * (a - b)                                           ! NEW
         d2 = 0.5d0 * (a - 1.0d0) * (a - b) * (b - a + 1.0d0) * (2.0d0 - a)   ! NEW
         ! log(1 + d1/x + d2/x^2)                                             ! NEW
         log_result = log_result + log(1.0d0 + d1/x + d2/(x*x))               ! NEW
       end if                                                                 ! NEW


!      ! Add first correction term for better accuracy                        !
!      if (x > 10.0d0 .and. x <= 200.0d0) then                                !
!         ! log(1 + (a-1)*(a-b)/x) ≈ (a-1)*(a-b)/x                            ! OLD
!         log_result = log_result + log(1.0d0 + (a - 1.0d0) * (a - b) / x)    !
!      end if                                                                 !
       
       if (log_result > LOG_HUGE) then
          integral = huge(1.0d0)
          status = STATUS_OVERFLOW
          return
       else if (log_result < LOG_TINY) then
          integral = 0.0d0
          status = STATUS_UNDERFLOW
          return
       else
          integral = exp(log_result)
!         if (x <= 200.0d0) then ! OLD
          if (x  >  10.0d0) then ! NEW
             status = STATUS_GSL_ERROR
          else
             status = STATUS_SUCCESS
          end if
          return
       end if
    end if

    ! Should never reach here
    integral = 0.0d0
    status = STATUS_GSL_ERROR

  end subroutine radial_integral




  !============================================================================
  ! SUBROUTINE: radial_integral_log
  !============================================================================
  subroutine radial_integral_log(n, l, mu, T, log_integral, status)
    integer, intent(in) :: n, l
    real(8), intent(in) :: mu, T
    real(8), intent(out) :: log_integral
    integer, intent(out) :: status

    real(8) :: a, b, x, log_hyperg, hyperg
    real(8) :: log_prefactor, log_2_factor
    real(8) :: log_mu_factor, log_T_factor, log_gamma_num, log_gamma_den
    real(8) :: d1, d2  ! NEW
    type(gsl_sf_result), target :: gsl_result
    integer(c_int) :: gsl_status
    logical :: use_asymptotic

    log_integral = -huge(1.0d0)
    status = STATUS_SUCCESS
    use_asymptotic = .false.

    ! Input validation
    if (n < 0 .or. l < 0) then
       log_integral = -huge(1.0d0)
       status = STATUS_INVALID_PARAM
       return
    end if

    if (mu <= 0.0d0) then
       log_integral = -huge(1.0d0)
       status = STATUS_INVALID_PARAM
       return
    end if

    ! T = 0 case
    if (T <= 0.0d0) then
       if (l == 0) then
          log_integral = -log(2.0d0) - (dble(n+1)/2.0d0) * log(mu) + &
                         gsl_sf_lngamma((n+1)/2.0d0)
          status = STATUS_SUCCESS
       else
          log_integral = -huge(1.0d0)
          status = STATUS_ZERO_BY_RULE
       end if
       return
    end if

    ! General case: T > 0
    a = 0.5d0 * dble(n + l + 1)
    b = dble(l) + 1.5d0
    x = mu * T * T

    ! Compute log of prefactor (excluding hypergeometric)
    ! CORRECTED 2026-09-10: √π/4, mu^{(l-n-1)/2}, no 2-power factor
    ! jaqg 2026-09-09: bug 1 — wrong exponents (see radial_integral), same fix.
    log_prefactor = 0.5d0 * log(PI) - 2.0d0 * log(2.0d0) - mu * T * T
    log_2_factor = 0.0d0
    log_mu_factor = 0.5d0 * dble(l - n - 1) * log(mu)
    log_T_factor = dble(l) * log(T)
    log_gamma_num = gsl_sf_lngamma(a)
    log_gamma_den = gsl_sf_lngamma(b)

    log_prefactor = log_prefactor + log_2_factor + log_mu_factor + &
                    log_T_factor + log_gamma_num - log_gamma_den

    ! Check if hypergeometric argument is large
!   if (x > 100.0d0) then ! OLD
    if (x > 300.0d0) then ! NEW
       use_asymptotic = .true.
    else
       ! Try GSL for moderate x
       gsl_status = gsl_sf_hyperg_1F1_e(a, b, x, c_loc(gsl_result))
       if (gsl_status == 0) then
          hyperg = gsl_result%val
          if (hyperg > 0.0d0 .and. hyperg == hyperg) then
             log_integral = log_prefactor + log(hyperg)
             if (log_integral < -1.0d300) then
                log_integral = -huge(1.0d0)
                status = STATUS_UNDERFLOW
             else
                status = STATUS_SUCCESS
             end if
             return
          end if
       end if
       use_asymptotic = .true.
    end if

    ! Asymptotic approximation
    if (use_asymptotic) then
!      log_hyperg = gsl_sf_lngamma(b) - gsl_sf_lngamma(a) + x + (a - b) * log(x)   ! OLD
!      log_integral = log_prefactor + log_hyperg                                   ! OLD

       log_hyperg = gsl_sf_lngamma(b) - gsl_sf_lngamma(a) + x + (a - b) * log(x) ! NEW
       if (x > 10.0d0) then                                                      ! NEW
         d1 = (a - 1.0d0) * (a - b)                                              ! NEW
         d2 = 0.5d0 * (a - 1.0d0) * (a - b) * (b - a + 1.0d0) * (2.0d0 - a)      ! NEW
         log_hyperg = log_hyperg + log(1.0d0 + d1/x + d2/(x*x))                  ! NEW
       end if                                                                    ! NEW
       log_integral = log_prefactor + log_hyperg                                 ! NEW
       
       if (log_integral < -1.0d300) then
          log_integral = -huge(1.0d0)
          status = STATUS_UNDERFLOW
       else if (x <= 100.0d0) then
          status = STATUS_GSL_ERROR
       else
          status = STATUS_SUCCESS
       end if
       return
    end if

    ! Should never reach here
    log_integral = -huge(1.0d0)
    status = STATUS_GSL_ERROR

  end subroutine radial_integral_log

end module intracule_integral
