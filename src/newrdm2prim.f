!-----------------------------------------------------------------------
!
!-----Compute the 2-RDM over Primitive Cartesian Gaussians (CG) from the
!     standard 2-RDM over Canonical MOs.
!
      subroutine newrdm2prim (rho2prim)
      USE space_for_wfncoef
      USE mod_rdm
      implicit real(kind=8) (a-h,o-z)
      include   'wfn.inc'
!
!-----The returned 2-RDM over CGs
!
      real (kind=8),intent(out) :: rho2prim(nprims,nprims,nprims,nprims)
!
!-----Working space arrays to perform the four-index transformation from
!     MOs to CGs
!
      real(kind=8), allocatable,dimension (:,:,:,:)  :: rho2ijkd
      real(kind=8), allocatable,dimension (:,:,:,:)  :: rho2ijcd
      real(kind=8), allocatable,dimension (:,:,:,:)  :: rho2ibcd
!
!-----Local variables
!
      integer (kind=4) :: a,b,c,d,i,j,k,l,ab,cd
      real    (kind=8) :: val,c1,c1c2,c1c2c3,prodc

      allocate (rho2ijkd(nmo,nmo,nmo,nprims),stat=ier)
      if (ier.ne.0) goto 1000
      allocate (rho2ijcd(nmo,nmo,nprims,nprims),stat=ier)
      if (ier.ne.0) goto 1000
      allocate (rho2ibcd(nmo,nprims,nprims,nprims),stat=ier)
      if (ier.ne.0) goto 1000

      rho2prim = 0.0d+00
      ab = 0
      do a=1,nprims
        do b=1,nprims
          ab = ab + 1
          cd = 0
          do c=1,nprims
            do d=1,nprims
              cd = cd + 1
              if (cd.le.ab) then
                do i=1,nmo
                  do j=1,nmo
                    do k=1,nmo
                      rho2ijkd(i,j,k,d) = dot_product(
     &                rdm2(i,j,k,1:nmo),coef(nmo+1:nmo+nmo,d))
                    enddo
                  enddo
                enddo
                do i=1,nmo
                  do j=1,nmo
                    rho2ijcd(i,j,c,d) = dot_product(
     &              rho2ijkd(i,j,1:nmo,d),coef(nmo+1:nmo+nmo,c))
                  enddo
                enddo
                do i=1,nmo
                  rho2ibcd(i,b,c,d) = dot_product(
     &            rho2ijcd(i,1:nmo,c,d),coef(nmo+1:nmo+nmo,b))
                enddo
                val = dot_product(
     &          rho2ibcd(1:nmo,b,c,d),coef(nmo+1:nmo+nmo,a))
              endif
              rho2prim(a,b,c,d) = val
              rho2prim(c,d,a,b) = val
              rho2prim(b,a,d,c) = val
              rho2prim(d,c,b,a) = val
            enddo
          enddo
        enddo
      enddo
      deallocate (rho2ijkd,stat=ier)
      if (ier.ne.0) stop '# newrdm2prim.f: Cannot deallocate rho2ijkd()'
      deallocate (rho2ijcd,stat=ier)
      if (ier.ne.0) stop '# newrdm2prim.f: Cannot deallocate rho2ijcd()'
      deallocate (rho2ibcd,stat=ier)
      if (ier.ne.0) stop '# newrdm2prim.f: Cannot deallocate rho2ibcd()'
      return

 1000 continue
!
!-----Calculation not involving the auxiliary arrays 
!     rho2ijkd(), rho2ijcd(), and rho2ibcd
!
      rho2prim = 0.0d+00
      ab = 0
      do a=1,nprims
        do b=1,nprims
          ab = ab + 1
          cd = 0
          do c=1,nprims
            do d=1,nprims
              cd = cd + 1
              val = 0.0d+00
              if (cd.le.ab) then
                val = 0.0d+00
                do i=1,nmo
                  c1 = coef(nmo+i,a)
                  do j=1,nmo
                    c1c2 = c1 * coef(nmo+j,b)
                    do k=1,nmo
                      c1c2c3 = c1c2 * coef(nmo+k,c)
                      do l=1,nmo
                        prodc = c1c2c3 * coef(nmo+l,d)
                        val = val + rdm2(i,j,k,l) * prodc
                      enddo
                    enddo
                  enddo
                enddo
              endif
              rho2prim(a,b,c,d) = val
              rho2prim(c,d,a,b) = val
              rho2prim(b,a,d,c) = val
              rho2prim(d,c,b,a) = val
            enddo
          enddo
        enddo
      enddo
      return
      end
