      subroutine cmcq()
c
      use mod_prec, only : dp, ip
      use space_for_wfnbasis
c
      include 'implicit.inc'
      include 'error.inc'
      include 'cmcq.inc'
      include 'io.inc'
      include 'datatm.inc'
      include 'wfn.inc'
c
      real(kind=dp) :: wi,xx,yy,zz
      integer(kind=ip) :: nrot
c
      cm = 0d0
      cq = 0d0
      wt = 0d0
      zt = 0d0
      imatrix = 0d0
c
c     compute the center of mass and charge
c
      do i = 1, ncent
        wt = wt + wgatm(int(charget(i)))
        cm(1) = cm(1) + wgatm(int(charget(i))) * xyz(i,1)
        cm(2) = cm(2) + wgatm(int(charget(i))) * xyz(i,2)
        cm(3) = cm(3) + wgatm(int(charget(i))) * xyz(i,3)
        zt = zt + charget(i)             
        cq(1) = cq(1) + charget(i) * xyz(i,1) 
        cq(2) = cq(2) + charget(i) * xyz(i,2) 
        cq(3) = cq(3) + charget(i) * xyz(i,3) 
      enddo
      cm = cm / wt
      cq = cq / zt
c
      write(uout,'(" # Center of mass : ",3(F8.4))') cm
      write(uout,'(" # Center of charge : ",3(F8.4))') cq
c
c     compute the inertia matrix and diagonalize it
c
      do i = 1, ncent
        wi = wgatm(i)
        xx = xyz(i,1)
        yy = xyz(i,2)
        zz = xyz(i,3)
        imatrix(1,1) = imatrix(1,1) + wi * (yy*yy + zz*zz)
        imatrix(2,2) = imatrix(2,2) + wi * (zz*zz + xx*xx)
        imatrix(3,3) = imatrix(3,3) + wi * (xx*xx + yy*yy)
        imatrix(1,2) = imatrix(1,2) - wi * xx * yy
        imatrix(2,3) = imatrix(2,3) - wi * yy * zz
        imatrix(3,1) = imatrix(3,1) - wi * zz * xx
      enddo
      imatrix(2,1) = imatrix(1,2)
      imatrix(3,2) = imatrix(2,3)
      imatrix(1,3) = imatrix(3,1)
c
      call jacobi(imatrix,3,3,einertia,dimatrix,nrot)
      if (nrot .lt. 0) call error('cmcq','fail diag inertia',faterr)
c
      write (uout,'(" # Inertia eigenvalues (au) : ",3(1x,E15.8))') 
     &  (einertia(i),i=1,3)
c
      end subroutine
