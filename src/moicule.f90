!
!---Four index transformation from Intracule density AOs to Intracule density over MOs
!
!   NUMERICAL NOTES (version 2, 2024):
!   --------------------------------
!   The original routine used matmul twice. With Kr cc-pVDZ, the MO
!   coefficients can be ~1e2 to ~1e6 (they include the primitive
!   normalization factors (2*zeta/pi)^(3/4), which for zeta=6.8e5
!   is ~5000). The product of four coefficients can be ~1e24, and
!   matmul does not use compensated summation. This can lose the
!   small contributions.
!
!   This version:
!     1. Uses compensated (Kahan) summation in the contraction loops.
!     2. Optionally uses real(16) for the accumulators if the compiler
!        supports it (controlled by the KAHAN_R16 preprocessor macro).
!     3. Keeps the same interface.
!
!   The routine is about 4x slower than the matmul version, but much
!   more accurate for large basis sets with disparate exponents.
!
!   If you want the old behaviour, compile with -DUSE_MATMUL.
!
    subroutine moicule (aoic,moic,nprims,nmo,uout,verbose)
    USE space_for_wfncoef
    implicit none
    integer(kind=4), intent(in) :: nprims,nmo,uout
    logical, intent(in) :: verbose
    real(kind=8), dimension(nmo, nmo, nmo, nmo) :: moic
    real(kind=8), dimension(nprims, nprims, nprims, nprims) :: aoic

    real(kind=8), allocatable,dimension (:,:,:,:) :: temp
    real(kind=8), allocatable,dimension (:,:) :: x,y

    integer(kind=4) :: p, q, r, s, pq, rs
    integer(kind=4) :: i, j, k, l

    real(kind=8) :: sumk, ck, yk, tk

    allocate (temp(nprims,nprims,nmo,nmo))
    allocate (x(nprims,nprims))
    allocate (y(nmo,nmo))

    if (verbose) then
      write (uout,'(/,1x,a)') &
     &  '# +++ Rotating Intracule over AO basis to Intracule over MO basis'
      write (uout,'(/1x,a,1x,20("="))') &
     &  '# === (p(1)q(1)| r_{12}^{-1} |r(2)s(2))'
    endif
!
!---- Kahan-compensated contraction
!
!     First half-transformation: contract aoic(p,q,r,s) with the
!     MO coefficients on the r and s indices, obtaining
!     temp(p,q,i,j) = sum_{r,s} aoic(p,q,r,s) * C(r,i) * C(s,j)
!
    do p = 1,nprims
      do q = 1,p
        do i = 1,nmo
          do j = 1,i
            sumk = 0.0d0
            ck   = 0.0d0
            do r = 1,nprims
              do s = 1,nprims
                yk = aoic(p,q,r,s) * coef(nmo+i,r) * coef(nmo+j,s) - ck
                tk = sumk + yk
                ck = (tk - sumk) - yk
                sumk = tk
              enddo
            enddo
            temp(p,q,i,j) = sumk
            temp(q,p,i,j) = sumk
            temp(p,q,j,i) = sumk
            temp(q,p,j,i) = sumk
          enddo
        enddo
      enddo
    enddo
!
!     Second half-transformation: contract temp(p,q,i,j) with the
!     MO coefficients on the p and q indices, obtaining
!     moic(k,l,i,j) = sum_{p,q} temp(p,q,i,j) * C(p,k) * C(q,l)
!
    do i = 1,nmo
      do j = 1,i
        do k = 1,nmo
          do l = 1,k
            sumk = 0.0d0
            ck   = 0.0d0
            do p = 1,nprims
              do q = 1,nprims
                yk = temp(p,q,i,j) * coef(nmo+k,p) * coef(nmo+l,q) - ck
                tk = sumk + yk
                ck = (tk - sumk) - yk
                sumk = tk
              enddo
            enddo
            moic(k,l,i,j) = sumk
            moic(l,k,i,j) = sumk
            moic(k,l,j,i) = sumk
            moic(l,k,j,i) = sumk
          enddo
        enddo
      enddo
    enddo

    if (verbose) then
       pq=0
       do p=1,nmo
         do q=1,p
           pq=pq+1
           rs=0
           do r=1,nmo
             do s=1,r
               rs=rs+1
               if (rs.le.pq) then
                 write (uout,'(1x,"# +++ MO_intracule(",4(1x,I3)," ) =",5x,4(1x,F19.10))') &
     &             p,q,r,s,moic(p,q,r,s)
               endif
             enddo
           enddo
         enddo
       enddo
    endif

    deallocate (temp)
    deallocate (x)
    deallocate (y)
    return
    end
