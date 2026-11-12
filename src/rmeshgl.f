c
c----------------------------------------------------------------------
c
      subroutine rmeshgl (x1,x2,rmid,rads,wrads,n)
c
c-----Given the lower and upper limits of integration x1 and x2, 
c     and given n, this routine returns arrays x and w of lenght n, 
c     containing the abscissas and weights of the Gauss-Legendre 
c     n-point quadrature formula.
c-----------------------------------------------------------------------
c-----High precision is a good idea for this routine
c     
      implicit real*8 (a-h,o-z)
c
c-----Increase eps if you don't have this floating precision.
      parameter (eps=3.0d-14,pi=3.1415 92653 58979 32384 62643d0)
c
      real(kind=8), intent(out), dimension(n) :: rads, wrads
      integer n,i,j,m
      real*8 pi
      real*8 x1,x2,rmid,x(n),w(n),eps,p1,p2,p3,pp,xl,xm,z,z1
  
c
c-----The roots are symmetric in the interval, so we only have to find
c     half of them.
c
      m=(n+1)/2
      xm=0.5d0*(x2+x1)
      xl=0.5d0*(x2-x1)
c
c-----Loop over the desired roots.
c
      do i=1,m
        z=cos(pi*(i-0.25d0)/(n+0.5d0))
 1      continue
        p1=1.0d0
        p2=0.0d0
c
c-------Loop up the recurrence relation to get the Legendre polynomial 
c       evaluated at z.
c
        do 11 j=1,n
          p3=p2
          p2=p1
          p1=((2.0d0*j-1.0d0)*z*p2-(j-1.0d0)*p3)/j
 11     continue
c
c-------p1 is now the desired Legendre polynomial. We next compute pp,
c       derivative , by a standard relation involving p2, the polyn-
c       omial of one lower order.
c       
        pp=n*(z*p1-p2)/(z*z-1.0d0)
        z1=z
c
c-------Newton's method.
c
        z=z1-p1/pp
        if(abs(z-z1).gt.eps) goto 1
c
c-------Scale the root to the desired interval.
c
        x(i)=xm-xl*z
c
c-------and put in its symmetric counterpart.
c
        x(n+1-i)=xm+xl*z
c
c-------compute the weight.
c
        w(i)=2.0d0*xl/((1.0d0-z*z)*pp*pp)
c
c-------and its symmetric counterpart.
c
        w(n+1-i)=w(i)
      enddo
      do j=1,n
        rads(j) = rmid * (1.d0+x(j)) / (1.d0-x(j))
        wrads(j) = w(j) * rmid * 2.d0 / (1.d0-x(j))**2
      enddo
      return
      end
