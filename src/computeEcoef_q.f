!
!-----------------------------------------------------------------------
!     binary128 twin of computeEcoef.f: E_t^ij coefficients, Gaussian
!     prefactors and pair data in real(16), stored in binary128_kernels.
!     Used only by the binary128 path of radial_integ.f.
!-----------------------------------------------------------------------
!
      subroutine computeEcoef_q ()
      USE        space_for_wfnbasis
      USE        space_for_othershells
      USE        binary128_kernels
      implicit real(kind=8) (a-h,o-z)
      include   'param.inc'
      include   'wfn.inc'
      integer(kind=4) p,q,pq,lx,npairshell,ica,icb,j,la,lb
      real(kind=16) ax(3),bx(3),za,zb,d2(3)
      real(kind=16), allocatable :: ce(:,:,:,:)

      lx = maxlprim
      npairshell = nshells*(nshells+1)/2
      if (allocated(etij_q)) deallocate (etij_q,psmall_q,xmu_q,
     &                                   prefac_q,pbig_q)
      allocate (etij_q(npairshell,0:2*lx,0:lx,0:lx,3))
      allocate (psmall_q(npairshell),xmu_q(npairshell))
      allocate (prefac_q(npairshell,4),pbig_q(npairshell,3))
      allocate (ce(-1:2*lx,-1:lx,-1:lx,3))
      etij_q = 0.0_16
      pq = 0
      do p=1,nshells
        ica = center(p)
        ax(1:3) = real(xyz(ica,1:3),16)
        la = lang(p)
        za = real(zeta(p),16)
        do q=1,p
          pq = pq + 1
          icb = center(q)
          bx(1:3) = real(xyz(icb,1:3),16)
          lb = lang(q)
          zb = real(zeta(q),16)
          psmall_q(pq) = za+zb
          pbig_q(pq,:) = (za*ax(:)+zb*bx(:))/(za+zb)
          xmu_q(pq) = za*zb/(za+zb)
          d2(:) = (ax(:)-bx(:))**2
          prefac_q(pq,1:3) = exp(-xmu_q(pq)*d2(1:3))
          prefac_q(pq,4) = prefac_q(pq,1)*prefac_q(pq,2)*prefac_q(pq,3)
          do j=1,3
            call etijcalc_q (j,lx,la,lb,ce,za,zb,ax(j),bx(j))
            etij_q(pq,0:2*lx,0:lx,0:lx,j) = ce(0:2*lx,0:lx,0:lx,j)
          enddo
        enddo
      enddo
      deallocate (ce)
      return
      end
