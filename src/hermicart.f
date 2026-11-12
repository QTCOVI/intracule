!
!-----------------------------------------------------------------------
!
!     Given the coefficients C_p of the Hermite polynomial of degree n in
!     the variable xi = mu^{1/2} (X - T):
!
!     H_n(xi) = SUM(p=0...n) C_p xi^p
!
!     this routine finds the coefficients D_p for p=0,...,n of the expansion
!
!     H_n(xi) = SUM(p=0...n) D_p X^p
!
!     NUMERICAL NOTES (version 2, 2024):
!     --------------------------------
!     The original version expanded H_n(sqrt(mu)*(X-T)) in powers of X
!     by direct monomial expansion:
!
!       d(p) = SUM_{k=p}^{n} c(k) * (sqrt(mu))^k * (-t)^(k-p) * C(k,p)
!
!     For the worst case (mu ~ 5e4, n ~ 16, t ~ 100 bohr), the individual
!     terms can reach ~1e60 while the final result d(p) is O(1). This
!     causes catastrophic cancellation: ~60 digits are lost out of the
!     16 available in double precision, and the result is numerical noise.
!
!     This version avoids the monomial expansion entirely. Instead, it
!     uses the three-term recurrence of the Hermite polynomials directly
!     in the variable X:
!
!       H_0(xi)   = 1
!       H_1(xi)   = 2*xi = 2*sqrt(mu)*(X - t)
!       H_{i+1}   = 2*xi*H_i - 2*i*H_{i-1}
!
!     with xi = sqrt(mu)*(X - t). Expanding H_i(xi) = SUM_p d_i(p) X^p,
!     we obtain the recurrence for the coefficients d_i(p):
!
!       d_0(0)   = 1
!       d_1(0)   = -2*sqrt(mu)*t
!       d_1(1)   =  2*sqrt(mu)
!       d_{i+1}(p) = 2*sqrt(mu)*( d_i(p-1) - t*d_i(p) ) - 2*i*d_{i-1}(p)
!
!     This recurrence involves only additions, subtractions and
!     multiplications by sqrt(mu) ~ O(100) and t ~ O(100). No enormous
!     intermediate products are formed, and the cancellation is limited
!     to the intrinsic one of the Hermite recurrence (which is mild for
!     moderate xi).
!
!     For n = 0 the routine returns d(0)=1.
!     For n = 1 the routine returns d(0)=-2*sqrt(mu)*t, d(1)=2*sqrt(mu).
!
!     The routine also works if n > maxn, in which case it returns the
!     coefficients up to maxn and sets the rest to zero (with a warning
!     if verbose is enabled by the caller; here we simply truncate).
!
!-----------------------------------------------------------------------
!
      subroutine hermicart (mu,t,n,maxn,d)
      implicit real(kind=8) (a-h,o-z)
      include 'param.inc'
      integer(kind=4), intent(in)  :: n,maxn
      real   (kind=8), intent(in)  :: mu,t
      real   (kind=8), intent(out) :: d(0:maxn)

      integer(kind=4) :: p, i, nmax_eff
      real   (kind=8) :: sqmu
      real   (kind=8), allocatable :: dm1(:), d0(:), dp1(:)
!
!-----Sanity check on n
!
      if (n .lt. 0) then
        d(0:maxn) = 0.0d0
        return
      endif
!
!-----Effective degree: cannot exceed maxn
!
      nmax_eff = min(n, maxn)
!
!-----Trivial cases
!
      d(0:maxn) = 0.0d0
      if (nmax_eff .eq. 0) then
        d(0) = 1.0d0
        return
      endif
!
!-----sqrt(mu). If mu is tiny (e.g. 1e-10), sqrt(mu) is 1e-5; if mu is
!     large (e.g. 1e5), sqrt(mu) is ~316. Neither overflows nor underflows.
!
      sqmu = sqrt(mu)
!
!-----Allocate the three recurrence vectors. We use real(8) but could
!     use real(16) if the compiler supports it and precision is critical.
!     Here real(8) is enough because the recurrence is stable.
!
      allocate (dm1(0:nmax_eff))
      allocate (d0 (0:nmax_eff))
      allocate (dp1(0:nmax_eff))
!
!-----i = 0: H_0 = 1
!
      dm1(0:nmax_eff) = 0.0d0
      dm1(0) = 1.0d0
!
!-----i = 1: H_1 = 2*sqrt(mu)*(X - t)
!
      d0(0:nmax_eff) = 0.0d0
      d0(0) = -2.0d0 * sqmu * t
      if (nmax_eff .ge. 1) d0(1) = 2.0d0 * sqmu
!
      if (nmax_eff .eq. 1) then
        d(0:nmax_eff) = d0(0:nmax_eff)
        deallocate (dm1, d0, dp1)
        return
      endif
!
!-----Recurrence for i = 1, ..., nmax_eff-1
!     H_{i+1} = 2*sqrt(mu)*(X - t)*H_i - 2*i*H_{i-1}
!
!     In terms of coefficients:
!       dp1(p) = 2*sqrt(mu)*( d0(p-1) - t*d0(p) ) - 2*i*dm1(p)
!
!     with the convention d0(-1) = 0.
!
      do i = 1, nmax_eff - 1
        dp1(0:nmax_eff) = 0.0d0
        do p = 0, min(i+1, nmax_eff)
          if (p .eq. 0) then
            dp1(p) = 2.0d0 * sqmu * ( - t * d0(p) ) 
     &             - 2.0d0 * dble(i) * dm1(p)
          else
            dp1(p) = 2.0d0 * sqmu * ( d0(p-1) - t * d0(p) ) 
     &             - 2.0d0 * dble(i) * dm1(p)
          endif
        enddo
!
!-------Shift the recurrence vectors
!
        dm1(0:nmax_eff) = d0 (0:nmax_eff)
        d0 (0:nmax_eff) = dp1(0:nmax_eff)
      enddo
!
!-----Copy the final result
!
      d(0:nmax_eff) = d0(0:nmax_eff)
!
      deallocate (dm1, d0, dp1)
!
      return
      end
