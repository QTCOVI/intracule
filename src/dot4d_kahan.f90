!-----------------------------------------------------------------------
!     Compensated (Kahan) dot product of two 4D arrays of size
!     (n1,n2,n3,n4). Returns sum_{p,q,r,s} a(p,q,r,s)*b(p,q,r,s).
!-----------------------------------------------------------------------
      real(kind=8) function dot4d_kahan (a, b, n1, n2, n3, n4)
      implicit none
      integer(kind=4), intent(in) :: n1, n2, n3, n4
      real(kind=8),    intent(in) :: a(n1,n2,n3,n4)
      real(kind=8),    intent(in) :: b(n1,n2,n3,n4)

      integer(kind=4) :: p, q, r, s
      real(kind=8)    :: sumk, c, y, t

      sumk = 0.0d0
      c    = 0.0d0
      do s = 1,n4
        do r = 1,n3
          do q = 1,n2
            do p = 1,n1
              y    = a(p,q,r,s) * b(p,q,r,s) - c
              t    = sumk + y
              c    = (t - sumk) - y
              sumk = t
            enddo
          enddo
        enddo
      enddo
      dot4d_kahan = sumk
      return
      end
