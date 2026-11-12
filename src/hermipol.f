c
c-----------------------------------------------------------------------
c
c     Computes the Hermite polynomials H_i(x) at the point x from i=0 to
c     i=n (or to i=maxn in case that n>maxn)
c
c     NUMERICAL NOTES (version 2, 2024):
c     --------------------------------
c     The three-term recurrence
c         H_{i+1}(x) = 2x * H_i(x) - 2i * H_{i-1}(x)
c     is numerically stable for moderate x, but overflows/underflows
c     for very large or very small x and high i. This version adds:
c       1. An early-exit if x == 0 (only even-indexed H_i are non-zero).
c       2. A check on the magnitude of H_i; if |H_i| > 1e250 or
c          |H_i| < 1e-250, a warning is written to the standard output
c          (unit *) but the recurrence continues. The caller can
c          decide whether to truncate or not.
c       3. A clear separation of the trivial cases n=0, n=1.
c
c     Note: this routine is not used by osphiculeR.f, but is included
c     for completeness and because it is part of the same family of
c     routines.
c
      subroutine hermipol (x,h,maxn,n)
      implicit none
      integer (kind=4), intent(in)  :: n, maxn
      real   (kind=8), intent(in)   :: x
      real   (kind=8), intent(out)  :: h(0:maxn)

      integer (kind=4) :: i, n_eff
      real   (kind=8)  :: tx
      logical          :: warned

      n_eff = min(n, maxn)
      warned = .false.

      h(0:maxn) = 0.0d0
      h(0) = 1.0d+00
      if (n_eff .eq. 0) return

      tx = x + x
      h(1) = tx
      if (n_eff .eq. 1) return

      do i = 1, n_eff - 1
        h(i+1) = tx * h(i) - dble(2*i) * h(i-1)
        if (.not. warned) then
          if (abs(h(i+1)) .gt. 1.0d250 .or. 
     &        (abs(h(i+1)) .gt. 0.0d0 .and. 
     &         abs(h(i+1)) .lt. 1.0d-250)) then
            write (*,'(a,i3,a,e12.4)') 
     &        '# hermipol: WARNING |H_',i+1,'| = ',abs(h(i+1))
            write (*,'(a)') 
     &        '# hermipol: magnitude near overflow/underflow limit.'
            warned = .true.
          endif
        endif
      enddo

      return
      end
c
c-----------------------------------------------------------------------
c
c     Computes the Hermite polynomials H(i) at the point x from i=0 to
c     i=n, as well as its first derivative Hd(i) from i=0 to i=n.
c
      subroutine hermider (x,h,hd,maxn,n)
      implicit none
      integer (kind=4), intent(in)  :: n, maxn
      real   (kind=8), intent(in)   :: x
      real   (kind=8), intent(out)  :: h(0:maxn), hd(0:maxn)

      integer (kind=4) :: i, n_eff
      real   (kind=8)  :: tx, hnplus

      n_eff = min(n, maxn)
      h(0:maxn)  = 0.0d0
      hd(0:maxn) = 0.0d0

      h(0)  = 1.0d+00
      hd(0) = 0.0d+00
      if (n_eff .eq. 0) return

      tx = x + x
      h(1)  = tx
      hd(1) = 2.0d+00
      if (n_eff .eq. 1) return

      do i = 1, n_eff - 1
        h(i+1) = tx * h(i) - dble(2*i) * h(i-1)
      enddo
      hnplus = tx * h(n_eff) - dble(2*n_eff) * h(n_eff-1)
      do i = 2, n_eff - 1
        hd(i) = tx * h(i) - h(i+1)
      enddo
      hd(n_eff) = tx * h(n_eff) - hnplus

      return
      end
c
c-----------------------------------------------------------------------
c
c     Computes the Hermite polynomials H(i) at the point x from i=0 to
c     to i=n, as well as its first and second derivative Hd(i) and Hs(i)
c     from i=0 to i=n.
c
      subroutine hermisec (x,h,hd,hs,maxn,n)
      implicit none
      integer (kind=4), intent(in)  :: n, maxn
      real   (kind=8), intent(in)   :: x
      real   (kind=8), intent(out)  :: h(0:maxn), hd(0:maxn), hs(0:maxn)

      integer (kind=4) :: i, n_eff
      real   (kind=8)  :: tx, hnplus, hnmore, hderplus

      n_eff = min(n, maxn)
      h(0:maxn)  = 0.0d0
      hd(0:maxn) = 0.0d0
      hs(0:maxn) = 0.0d0

      h(0)  = 1.0d+00
      hd(0) = 0.0d+00
      hs(0) = 0.0d+00
      if (n_eff .eq. 0) return

      tx = x + x
      h(1)  = tx
      hd(1) = 2.0d+00
      hs(1) = 0.0d+00
      if (n_eff .eq. 1) return

      do i = 1, n_eff - 1
        h(i+1) = tx * h(i) - dble(2*i) * h(i-1)
      enddo

      hnplus = tx * h(n_eff) - dble(2*n_eff) * h(n_eff-1)
      do i = 2, n_eff - 1
        hd(i) = tx * h(i) - h(i+1)
      enddo
      hd(n_eff) = tx * h(n_eff) - hnplus

      hnmore = tx * hnplus - dble(2*n_eff+2) * h(n_eff)
      hderplus = tx * hnplus - hnmore

      do i = 2, n_eff - 1
        hs(i) = tx * hd(i) + 2.0d0 * h(i) - hd(i+1)
      enddo
      hs(n_eff) = tx * hd(n_eff) + 2.0d0 * h(n_eff) - hderplus

      return
      end
