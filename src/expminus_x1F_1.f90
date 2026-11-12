! ============================================================
!  MODULE: kummer_exp_neg
!  DESCRIPTION:
!     Numerically stable evaluation of the product
!
!         G(a,b,x) = exp(-x) * 1F1(a,b,x)
!
!     for REAL POSITIVE a, b, x.
!
!     The factor exp(-x) cancels the exponential growth of
!     1F1(a,b,x), so G(a,b,x) is an ALGEBRAIC function of x
!     in the asymptotic regime:
!
!         G(a,b,x) ~ Gamma(b)/Gamma(a) * x^(a-b) * S2(x)
!
!     This eliminates the overflow problem that affected
!     1F1(a,b,x) itself, and reduces the implementation to
!     just two methods:
!
!       - For 0 < x < X_SERIES (=30): DEFINING SERIES for G
!         (exponentially damped series, full precision)
!       - For x >= X_SERIES: ASYMPTOTIC expansion for G
!         (algebraic, no overflow, error ~ exp(-x))
!
!     In all cases, G(a,b,x) > 0 for a,b,x > 0.
!
!  REFERENCES:
!     - DLMF section 13.2 (defining series)
!     - DLMF section 13.7 (asymptotic expansion)
!     - Abramowitz & Stegun, Handbook of Mathematical
!       Functions, Chapter 13
!
!  COMPILATION:
!     gfortran -O2 -Wall -o kummer_exp_neg kummer_exp_neg.f90
! ============================================================
module kummer_exp_neg
    implicit none
    private

    ! Default precision (double precision)
    integer, parameter, public :: dp = kind(1.0d0)

    ! Mathematical constants
    real(dp), parameter :: PI = 3.141592653589793238462643383279502884197_dp

    ! --------------------------------------------------------
    ! Threshold
    ! --------------------------------------------------------
    ! Below this, use the defining series for G.
    ! Above this, use the asymptotic expansion for G.
    real(dp), parameter, public :: X_SERIES = 30.0_dp

    ! --------------------------------------------------------
    ! Series truncation parameters
    ! --------------------------------------------------------
    integer, parameter, public :: N_MAX_SERIES = 500
    integer, parameter, public :: N_MAX_ASYM   = 200
    real(dp), parameter, public :: TOL         = 1.0e-16_dp

    ! --------------------------------------------------------
    ! Public interfaces
    ! --------------------------------------------------------
    public :: kummer_G_dispatch
    public :: kummer_G_series
    public :: kummer_G_asym
    public :: log_gamma_pos

contains

    ! ========================================================
    !  log_gamma_pos: real log-Gamma for positive argument.
    !  Lanczos approximation with g=7, n=9 (precision ~1e-15).
    !  No reflection formula needed (z > 0 always).
    !  Declared RECURSIVE for consistency, though not needed.
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
    !  kummer_G_series: DEFINING SERIES for G(a,b,x)
    !  for a,b,x > 0 and moderate x.
    !
    !     G(a,b,x) = exp(-x) * 1F1(a,b,x)
    !              = exp(-x) * sum_{k=0}^inf (a)_k/(b)_k/k! * x^k
    !
    !  Convergence: rapid for x < 30. Full double precision
    !  reached with N ~ 50 terms.
    !
    !  Recurrence:
    !   t_0 = 1
    !   t_{k+1} = t_k * (a+k) * x / ((b+k) * (k+1))
    ! ========================================================
    function kummer_G_series(a, b, x) result(res)
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

        res = exp(-x) * sum_val

    end function kummer_G_series

    ! ========================================================
    !  kummer_G_asym: ASYMPTOTIC expansion for G(a,b,x)
    !  for a,b,x > 0 and large x.
    !
    !     G(a,b,x) ~ Gamma(b)/Gamma(a) * x^(a-b) * S2(x)
    !
    !  with
    !   S2(x) = sum_{k=0}^{N-1} (b-a)_k (1-a)_k / k! * x^(-k)
    !
    !  NO exponential factor remains: the exp(-x) has
    !  exactly cancelled the exp(x) of Kummer's expansion.
    !
    !  Recurrence:
    !   t_0 = 1
    !   t_{k+1} = t_k * (b-a+k)(1-a+k) / ((k+1) x)
    !
    !  Adaptive truncation: stop when terms start growing.
    ! ========================================================
    function kummer_G_asym(a, b, x) result(res)
        real(dp), intent(in) :: a, b, x
        real(dp) :: res

        real(dp) :: s2, t
        real(dp) :: log_prefactor
        real(dp) :: mag_t_prev, mag_t_cur
        integer  :: k

        ! -------- Series S2 (the asymptotic series) --------
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
        ! log( Gamma(b)/Gamma(a) * x^(a-b) )
        log_prefactor = log_gamma_pos(b) - log_gamma_pos(a) &
                        + (a - b) * log(x)

        ! -------- Assemble result --------
        res = exp(log_prefactor) * s2

    end function kummer_G_asym

    ! ========================================================
    !  kummer_G_dispatch: UNIFIED INTERFACE.
    !
    !  Uses the appropriate method based on x:
    !    - 0 < x < X_SERIES  -> defining series for G (method=1)
    !    - x >= X_SERIES     -> asymptotic expansion  (method=2)
    !
    !  Outputs:
    !    method : 1 or 2
    !    value  : real value of G(a,b,x)
    ! ========================================================
    subroutine kummer_G_dispatch (n, l, mu, T, method, value)
        real(dp), intent(in)  :: mu, T
        integer,  intent(in)  :: n, l
        integer,  intent(out) :: method    ! 1=series, 2=asym
        real(dp), intent(out) :: value

        real(dp) :: a, b, x

        ! Mathematical constants
        real(dp), parameter :: PI = 3.141592653589793238462643383279502884d0

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
        log_mu_factor = 0.5d0 * dble(l - n - 1) * log(mu)
        log_T_factor = dble(l) * log(T)

        log_gamma_num = log_gamma_pos(a)
        log_gamma_den = log_gamma_pos(b)

        log_prefactor = log_prefactor + log_mu_factor + log_T_factor
        log_prefactor = log_prefactor + log_gamma_num - log_gamma_den

        if (x < X_SERIES) then
            method = 1
            value  = kummer_G_series(a, b, x)
        else
            method = 2
            value  = kummer_G_asym(a, b, x)
        end if
        value = value * exp(log_prefactor)

    end subroutine kummer_G_dispatch

end module kummer_exp_neg
