! ============================================================
!  MODULE: kummer_pos
!  DESCRIPTION:
!     Numerically stable evaluation of the confluent
!     hypergeometric function 1F1(a,b,x) for REAL POSITIVE
!     a, b, x.
!
!     Strategy:
!       - For 0 < x < X_SERIES (=30): use the DEFINING SERIES,
!         which converges rapidly and gives full precision.
!       - For X_SERIES <= x < X_SWITCH (=100): use the
!         ASYMPTOTIC expansion in ordinary arithmetic.
!       - For x >= X_SWITCH: use the ASYMPTOTIC expansion
!         in logarithmic scale.
!
!     In all cases, 1F1(a,b,x) > 0 for a,b,x > 0.
!
!  REFERENCES:
!     - DLMF section 13.2 (defining series)
!     - DLMF section 13.7 (asymptotic expansion)
!
!  COMPILATION:
!     gfortran -O2 -Wall -o kummer_pos kummer_pos.f90
! ============================================================
module kummer_pos
    implicit none
    private

    integer, parameter, public :: dp = kind(1.0d0)

    real(dp), parameter :: PI = 3.141592653589793238462643383279502884197_dp

    ! --------------------------------------------------------
    ! Thresholds
    ! --------------------------------------------------------
    ! Below this, use the defining series
    real(dp), parameter, public :: X_SERIES   = 30.0_dp

    ! Above this, use the logarithmic asymptotic version
    real(dp), parameter, public :: X_SWITCH   = 100.0_dp

    ! Above this, overflow would occur in ordinary arithmetic
    real(dp), parameter, public :: X_OVERFLOW = 700.0_dp

    ! --------------------------------------------------------
    ! Series truncation parameters
    ! --------------------------------------------------------
    integer, parameter, public :: N_MAX_SERIES = 500
    integer, parameter, public :: N_MAX_ASYM   = 200
    integer, parameter, public :: N_MAX_LOG    = 2000
    real(dp), parameter, public :: TOL         = 1.0e-16_dp

    ! --------------------------------------------------------
    ! Public interfaces
    ! --------------------------------------------------------
    public :: kummer_1f1_dispatch
    public :: kummer_1f1_series
    public :: kummer_1f1_asym
    public :: kummer_1f1_log
    public :: log_gamma_pos

contains

    ! ========================================================
    !  log_gamma_pos: real log-Gamma for positive argument.
    ! ========================================================
    recursive function log_gamma_pos(z) result(res)
        real(dp), intent(in) :: z
        real(dp) :: res
        real(dp) :: zz, s
        real(dp) :: g, c(0:8)
        integer  :: i

        g = 7.0_dp
        c(0) =  0.99999999999980993_dp
        c(1) =  676.5203681218851_dp
        c(2) = -1259.1392167224028_dp
        c(3) =  771.32342877765313_dp
        c(4) = -176.61502916214059_dp
        c(5) =  12.507343278686905_dp
        c(6) =  -0.13857109526572012_dp
        c(7) =   9.9843695780195716e-6_dp
        c(8) =   1.5056327351493116e-7_dp

        zz = z - 1.0_dp
        s  = c(0)
        do i = 1, 8
            s = s + c(i) / (zz + real(i,dp))
        end do
        res = 0.5_dp*log(2.0_dp*PI) + (zz + 0.5_dp)*log(zz + g + 0.5_dp) &
              - (zz + g + 0.5_dp) + log(s)
    end function log_gamma_pos

    ! ========================================================
    !  kummer_1f1_series: DEFINING SERIES for a,b,x > 0.
    !
    !   1F1(a,b,x) = sum_{k=0}^inf (a)_k / (b)_k / k! * x^k
    !
    !  Converges for all x, but is most useful for moderate x.
    !  Returns the real value directly.
    ! ========================================================
    function kummer_1f1_series(a, b, x) result(res)
        real(dp), intent(in) :: a, b, x
        real(dp) :: res

        real(dp) :: term, sum_val
        integer  :: k

        term    = 1.0_dp
        sum_val = 1.0_dp
        do k = 0, N_MAX_SERIES-1
            term = term * (a + real(k,dp)) * x &
                        / ((b + real(k,dp)) * real(k+1,dp))
            sum_val = sum_val + term
            if (abs(term) < TOL * abs(sum_val)) exit
        end do
        res = sum_val

    end function kummer_1f1_series

    ! ========================================================
    !  kummer_1f1_asym: ASYMPTOTIC expansion for a,b,x > 0.
    !  Valid for x large (say x > 20-30).
    ! ========================================================
    function kummer_1f1_asym(a, b, x) result(res)
        real(dp), intent(in) :: a, b, x
        real(dp) :: res

        real(dp) :: s2, t
        real(dp) :: log_prefactor
        real(dp) :: mag_t_prev, mag_t_cur
        integer  :: k

        ! -------- Series S2 (associated with exp(x) x^(a-b)) --------
        s2 = 1.0_dp
        t  = 1.0_dp
        mag_t_prev = 1.0_dp
        do k = 0, N_MAX_ASYM-1
            t = t * (b - a + real(k,dp)) * (1.0_dp - a + real(k,dp)) &
                  / (real(k+1,dp) * x)
            mag_t_cur = abs(t)
            s2 = s2 + t
            if (abs(t) < TOL * abs(s2)) exit
            if (mag_t_cur > mag_t_prev .and. k > 1) exit
            mag_t_prev = mag_t_cur
        end do

        ! -------- Prefactor in log scale --------
        log_prefactor = log_gamma_pos(b) - log_gamma_pos(a) &
                        + x + (a - b) * log(x)

        res = exp(log_prefactor) * s2

    end function kummer_1f1_asym

    ! ========================================================
    !  kummer_1f1_log: LOGARITHMIC asymptotic expansion.
    !  Valid for any x > 0.
    ! ========================================================
    subroutine kummer_1f1_log(a, b, x, log_abs)
        real(dp), intent(in)  :: a, b, x
        real(dp), intent(out) :: log_abs

        real(dp) :: s2, t
        real(dp) :: log_prefactor
        real(dp) :: mag_t_prev, mag_t_cur
        integer  :: k

        s2 = 1.0_dp
        t  = 1.0_dp
        mag_t_prev = 1.0_dp
        do k = 0, N_MAX_LOG-1
            t = t * (b - a + real(k,dp)) * (1.0_dp - a + real(k,dp)) &
                  / (real(k+1,dp) * x)
            mag_t_cur = abs(t)
            s2 = s2 + t
            if (mag_t_cur > mag_t_prev .and. k > 1) exit
            mag_t_prev = mag_t_cur
        end do

        log_prefactor = log_gamma_pos(b) - log_gamma_pos(a) &
                        + x + (a - b) * log(x)

        log_abs = log_prefactor + log(s2)

    end subroutine kummer_1f1_log

    ! ========================================================
    !  kummer_1f1_dispatch: UNIFIED INTERFACE.
    !
    !  Uses the appropriate method based on x:
    !    - 0 < x < X_SERIES   -> defining series
    !    - X_SERIES <= x < X_SWITCH -> asymptotic (ordinary)
    !    - x >= X_SWITCH      -> asymptotic (logarithmic)
    ! ========================================================
    subroutine kummer_1f1_dispatch (n, l, mu, T, method, value, log_abs)
        real(dp), intent(in)  :: mu, T
        integer,  intent(in)  :: n, l
        integer,  intent(out) :: method    ! 1=series, 2=asym, 3=log
        real(dp), intent(out) :: value
        real(dp), intent(out) :: log_abs

        real(dp) :: a, b, x

        ! Mathematical constants
        real(8), parameter :: PI = 3.141592653589793238462643383279502884d0

        real(dp) :: log_prefactor, log_mu_factor, log_T_factor
        real(dp) :: log_gamma_num, log_gamma_den

        ! Input validation
        if (n < 0 .or. l < 0) then
           value = 0.0d0
           method = 0
           return
        end if
    
        if (mu <= 0.0d0) then
           value = 0.0d0
           method = 0
           return
        end if
    
        ! T = 0 case
        if (T <= 0.0d0) then
           if (l == 0) then
              value = 0.5d0 * (mu ** (-(n+1)/2.0d0)) * exp(log_gamma_pos((n+1)/2.0d0))
              method = 1
           else
              value = 0.0d0
              method = 1
           end if
           return
        end if


        a = 0.5d0 * dble(n + l + 1)
        b = dble(l) + 1.5d0
        x = mu * T * T

        log_prefactor = 0.5d0 * log(PI) - 2.0d0 * log(2.0d0)
        log_prefactor = log_prefactor - mu * T * T
        log_mu_factor = 0.5d0 * dble(l - n - 1) * log(mu)
        log_T_factor = dble(l) * log(T)

        log_gamma_num = log_gamma_pos(a)
        log_gamma_den = log_gamma_pos(b)

        log_prefactor = log_prefactor + log_mu_factor + log_T_factor
        log_prefactor = log_prefactor + log_gamma_num - log_gamma_den

        if (x < X_SERIES) then
            method  = 1
            value   = kummer_1f1_series(a, b, x)
            log_abs = log(value)
            value = exp(log_abs + log_prefactor)
        else if (x < X_SWITCH) then
            method  = 2
            value   = kummer_1f1_asym(a, b, x)
            log_abs = log(value)
            value = exp(log_abs + log_prefactor)
        else
            method  = 3
            call kummer_1f1_log(a, b, x, log_abs)
            value   = 0.0_dp   ! not used
            value = exp(log_abs + log_prefactor)
        end if

    end subroutine kummer_1f1_dispatch

end module kummer_pos
