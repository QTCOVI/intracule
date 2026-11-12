!-----------------------------------------------------------------------
!  binary128.f90
!
!  IEEE binary128 (real(16)) path of the analytic radial integration of
!  the intracule (radial_integ.f, norm pass).  Ported from the binary128
!  chain of J. A. Quinonero (v0.6.3, evaluators/mod_eval_bessic.f90 and
!  bessic/*.f90) into this program.
!
!  For heavy atoms (e.g. Kr2 with d shells and large |T|) the Jtab
!  l-sum  SUM_l SUM_m S_lm(T^) G_lm(nx,ny,nz) I(n,l,mu,T)  and its
!  contraction with the Hermite-Cartesian coefficients A_nx A_ny A_nz
!  cancel by 12-15 orders of magnitude, and binary64 loses the norm.
!  Here every object of that chain is rebuilt in binary128:
!
!    E_t^{ij} coefficients and Gaussian prefactors (computeEcoef_q.f)
!    Hermite -> Cartesian coefficients            (hermicart_q)
!    phi/theta power integrals and G_lm table      (build_gtab_q)
!    real spherical harmonics S_lm(T^)             (rsh_q)
!    radial integrals, terminating Kummer polynomial (radial_integral_q)
!
!  For the norm n = nx+ny+nz+2 and only l with (nx+ny+nz-l) even
!  contribute, so n-l = 2m with m >= 1 and
!     e^{-x} 1F1(a;b;x) = 1F1(1-m;b;-x),  x = mu T^2,
!  a polynomial of degree m-1 with positive terms; no library 1F1.
!
!  Mode (input keyword BINARY128):
!     b128_mode = 0  NEVER   binary64 only (previous behaviour)
!     b128_mode = 1  AUTO    binary64 per quartet with a rounding-error
!                            estimate; quartets whose estimated relative
!                            error exceeds b128_tol are recomputed here
!     b128_mode = 2  ALWAYS  every quartet in binary128
!-----------------------------------------------------------------------
      module binary128_kernels
      implicit none
      save
      integer, parameter :: qp = 16
      real(qp), parameter :: QPI = &
        3.141592653589793238462643383279502884197_qp

      integer  :: b128_mode = 1
      real(8)  :: b128_tol  = 1.0d-10
      integer(8) :: b128_nquad = 0, b128_ntot = 0
!
!-----binary128 copies of the shell-pair tables (computeEcoef_q.f)
!
      real(qp), allocatable :: etij_q(:,:,:,:,:)
      real(qp), allocatable :: psmall_q(:), xmu_q(:)
      real(qp), allocatable :: prefac_q(:,:), pbig_q(:,:)
!
!-----binary128 angular table G_lm(nx,ny,nz) and its dimensions
!
      real(qp), allocatable :: gtab_q(:,:,:,:,:)
      integer :: lmaxg_q = -1, nmaxc_q = -1
!
!-----Legendre / real-SH normalisation tables (binary128)
!
      real(qp), allocatable :: pcoef_q(:,:), qcoef_q(:,:,:), anlm_q(:,:)
      integer :: lmax_slm_q = -1

      contains
!
!-----------------------------------------------------------------------
!     INT_0^2pi sin^a cos^b dphi  (twin of phiintegrals.f)
!-----------------------------------------------------------------------
      subroutine phiintegrals_q (nxmax, nymax, xi)
      integer, intent(in) :: nxmax, nymax
      real(qp), intent(out) :: xi(0:nxmax,0:nymax)
      real(qp) :: xj(0:nxmax,0:nymax)
      integer :: nx, ny
      xj = 0.0_qp
      xj(0,0) = QPI/2.0_qp
      if (nxmax.ge.1) xj(1,0) = 1.0_qp
      if (nymax.ge.1) xj(0,1) = 1.0_qp
      do ny = 0, nymax
        do nx = 2, nxmax
          xj(nx,ny) = (real(nx-1,qp)/real(nx+ny,qp)) * xj(nx-2,ny)
        enddo
      enddo
      do nx = 0, nxmax
        do ny = 2, nymax
          xj(nx,ny) = (real(ny-1,qp)/real(nx+ny,qp)) * xj(nx,ny-2)
        enddo
      enddo
      do nx = 0, nxmax
        do ny = 0, nymax
          if (mod(nx,2).eq.0 .and. mod(ny,2).eq.0) then
            xi(nx,ny) = 4.0_qp * xj(nx,ny)
          else
            xi(nx,ny) = 0.0_qp
          endif
        enddo
      enddo
      end subroutine phiintegrals_q
!
!-----------------------------------------------------------------------
!     INT_0^pi sin^{p+1} cos^q dtheta  (twin of thetaintegrals.f)
!-----------------------------------------------------------------------
      subroutine thetaintegrals_q (maxnxy, maxnz, xi)
      integer, intent(in) :: maxnxy, maxnz
      real(qp), intent(out) :: xi(0:maxnxy,0:maxnz)
      real(qp) :: xj(0:maxnxy+1,0:maxnz)
      integer :: p, q, nxy, nz
      xj = 0.0_qp
      xj(0,0) = QPI/2.0_qp
      xj(1,0) = 1.0_qp
      if (maxnz.ge.1) xj(0,1) = 1.0_qp
      do q = 0, maxnz
        do p = 2, maxnxy+1
          xj(p,q) = (real(p-1,qp)/real(p+q,qp)) * xj(p-2,q)
        enddo
      enddo
      do p = 0, maxnxy+1
        do q = 2, maxnz
          xj(p,q) = (real(q-1,qp)/real(p+q,qp)) * xj(p,q-2)
        enddo
      enddo
      do nxy = 0, maxnxy
        do nz = 0, maxnz
          if (mod(nz,2).eq.0) then
            xi(nxy,nz) = 2.0_qp * xj(nxy+1,nz)
          else
            xi(nxy,nz) = 0.0_qp
          endif
        enddo
      enddo
      end subroutine thetaintegrals_q
!
!-----------------------------------------------------------------------
!     Legendre-derivative coefficients and real-SH normalisations
!     (twin of init_slmangint in slmangint.f90)
!-----------------------------------------------------------------------
      subroutine init_slmangint_q (lmax)
      integer, intent(in) :: lmax
      integer :: l, m, c, j
      real(qp) :: fac, pref, ratio
      if (lmax_slm_q .ge. lmax) return
      if (allocated(pcoef_q)) deallocate (pcoef_q)
      if (allocated(qcoef_q)) deallocate (qcoef_q)
      if (allocated(anlm_q))  deallocate (anlm_q)
      allocate (pcoef_q(0:lmax,0:lmax))
      allocate (qcoef_q(0:lmax,0:lmax,0:lmax))
      allocate (anlm_q(0:lmax,0:lmax))
      pcoef_q = 0.0_qp
      qcoef_q = 0.0_qp
      anlm_q  = 0.0_qp
      pcoef_q(0,0) = 1.0_qp
      if (lmax.ge.1) pcoef_q(1,1) = 1.0_qp
      do l = 2, lmax
        do c = 0, l
          pref = 0.0_qp
          if (c.ge.1) pref = real(2*l-1,qp) * pcoef_q(l-1,c-1)
          pcoef_q(l,c) = (pref - real(l-1,qp)*pcoef_q(l-2,c)) / real(l,qp)
        enddo
      enddo
      do l = 0, lmax
        do m = 0, l
          do c = 0, l-m
            fac = 1.0_qp
            do j = 1, m
              fac = fac * real(c+j,qp)
            enddo
            qcoef_q(l,m,c) = pcoef_q(l,c+m) * fac
          enddo
        enddo
      enddo
      do l = 0, lmax
        anlm_q(l,0) = sqrt(real(2*l+1,qp)/(4.0_qp*QPI))
        do m = 1, l
          ratio = 1.0_qp
          do j = l-m+1, l+m
            ratio = ratio / real(j,qp)
          enddo
          anlm_q(l,m) = sqrt(real(2*l+1,qp)/(2.0_qp*QPI) * ratio)
        enddo
      enddo
      lmax_slm_q = lmax
      end subroutine init_slmangint_q
!
!-----------------------------------------------------------------------
!     G_lm(nx,ny,nz) in binary128 (twin of slmangint in slmangint.f90)
!-----------------------------------------------------------------------
      function slmangint_q (l,m,nx,ny,nz,np,nt,yntphi,thetaval) result (g)
      integer, intent(in) :: l,m,nx,ny,nz,np,nt
      real(qp), intent(in) :: yntphi(0:np,0:np), thetaval(0:nt,0:nt)
      real(qp) :: g, gk, phipart, thetasum
      integer :: mm, k, kstart, c, aphi, bphi, pth, j
      g = 0.0_qp
      mm = abs(m)
      if (l .gt. mm+nx+ny+nz) return
      pth = mm + nx + ny
      if (pth .gt. nt) return
      kstart = 0
      if (m .lt. 0) kstart = 1
      do k = kstart, mm, 2
        aphi = (mm-k) + nx
        bphi = k + ny
        if (aphi.gt.np .or. bphi.gt.np) cycle
        phipart = yntphi(aphi,bphi)
        if (phipart .eq. 0.0_qp) cycle
        gk = 1.0_qp
        do j = 1, k
          gk = gk * real(mm-j+1,qp) / real(j,qp)
        enddo
        if (mod(k,2).eq.0) then
          if (mod(k/2,2).eq.1)     gk = -gk
        else
          if (mod((k-1)/2,2).eq.1) gk = -gk
        endif
        thetasum = 0.0_qp
        do c = l-mm, 0, -2
          if (c+nz .gt. nt) cycle
          thetasum = thetasum + qcoef_q(l,mm,c) * thetaval(pth,c+nz)
        enddo
        g = g + gk * phipart * thetasum
      enddo
      g = anlm_q(l,mm) * g
      end function slmangint_q
!
!-----------------------------------------------------------------------
!     Build gtab_q(0:lmax,-lmax:lmax,0:nmaxc,0:nmaxc,0:nmaxc) once
!     (twin of build_gtab; same index bounds as radial_integ.f)
!-----------------------------------------------------------------------
      subroutine build_gtab_q (lmax, nmaxc)
      integer, intent(in) :: lmax, nmaxc
      integer :: ntab, l, m, nx, ny, nz
      real(qp), allocatable :: yntphi(:,:), thetaval(:,:)
      if (allocated(gtab_q) .and. lmaxg_q.eq.lmax .and. &
          nmaxc_q.eq.nmaxc) return
      if (allocated(gtab_q)) deallocate (gtab_q)
      ntab = 2*lmax
      allocate (yntphi(0:ntab,0:ntab), thetaval(0:ntab,0:ntab))
      call phiintegrals_q (ntab, ntab, yntphi)
      call thetaintegrals_q (ntab, ntab, thetaval)
      call init_slmangint_q (lmax)
      allocate (gtab_q(0:lmax,-lmax:lmax,0:nmaxc,0:nmaxc,0:nmaxc))
      gtab_q = 0.0_qp
      !$omp parallel do collapse(3) default(shared) &
      !$omp   private(nx,ny,nz,l,m) schedule(dynamic)
      do nx = 0, nmaxc
        do ny = 0, nmaxc
          do nz = 0, nmaxc
            do l = 0, min(lmax, nx+ny+nz)
              do m = -l, l
                gtab_q(l,m,nx,ny,nz) = &
                  slmangint_q(l,m,nx,ny,nz,ntab,ntab,yntphi,thetaval)
              enddo
            enddo
          enddo
        enddo
      enddo
      !$omp end parallel do
      deallocate (yntphi, thetaval)
      lmaxg_q = lmax
      nmaxc_q = nmaxc
      end subroutine build_gtab_q
!
!-----------------------------------------------------------------------
!     Real spherical harmonics, no Condon-Shortley phase, binary128
!     (twin of compute_real_spherical_harmonics; S(l*(l+1)+m))
!-----------------------------------------------------------------------
      subroutine rsh_q (lmax, theta, phi, s)
      integer, intent(in) :: lmax
      real(qp), intent(in) :: theta, phi
      real(qp), intent(out) :: s(0:(lmax+1)**2-1)
      real(qp) :: nbar(0:lmax,0:lmax), x, sinth, sqrt2, a, b
      integer :: l, m
      x = cos(theta)
      sinth = sin(theta)
      sqrt2 = sqrt(2.0_qp)
      nbar = 0.0_qp
      nbar(0,0) = 1.0_qp/sqrt(4.0_qp*QPI)
      do m = 1, lmax
        nbar(m,m) = sqrt(real(2*m+1,qp)/real(2*m,qp))*sinth*nbar(m-1,m-1)
      enddo
      do m = 0, lmax
        if (m+1 .le. lmax) nbar(m+1,m) = sqrt(real(2*m+3,qp))*x*nbar(m,m)
        do l = m+2, lmax
          a = sqrt( (real(2*l-1,qp)*real(2*l+1,qp)) / &
                    (real(l-m,qp)*real(l+m,qp)) )
          b = sqrt( (real(2*l+1,qp)*real(l+m-1,qp)*real(l-m-1,qp)) / &
                    (real(2*l-3,qp)*real(l-m,qp)*real(l+m,qp)) )
          nbar(l,m) = a*x*nbar(l-1,m) - b*nbar(l-2,m)
        enddo
      enddo
      do l = 0, lmax
        s(l*(l+1)) = nbar(l,0)
        do m = 1, l
          s(l*(l+1)+m) = sqrt2 * nbar(l,m) * cos(real(m,qp)*phi)
          s(l*(l+1)-m) = sqrt2 * nbar(l,m) * sin(real(m,qp)*phi)
        enddo
      enddo
      end subroutine rsh_q
!
!-----------------------------------------------------------------------
!     Hermite -> Cartesian coefficients by the three-term recurrence
!     (twin of hermicart.f):  H_n(sqrt(mu)(X-t)) = SUM_p d(p) X^p
!-----------------------------------------------------------------------
      subroutine hermicart_q (mu, t, n, maxn, d)
      integer, intent(in) :: n, maxn
      real(qp), intent(in) :: mu, t
      real(qp), intent(out) :: d(0:maxn)
      real(qp) :: sqmu, dm1(0:maxn), d0(0:maxn), dp1(0:maxn)
      integer :: i, p, ne
      d = 0.0_qp
      if (n .lt. 0) return
      ne = min(n, maxn)
      if (ne .eq. 0) then
        d(0) = 1.0_qp
        return
      endif
      sqmu = sqrt(mu)
      dm1 = 0.0_qp
      dm1(0) = 1.0_qp
      d0 = 0.0_qp
      d0(0) = -2.0_qp*sqmu*t
      d0(1) =  2.0_qp*sqmu
      do i = 1, ne-1
        dp1 = 0.0_qp
        dp1(0) = -2.0_qp*sqmu*t*d0(0) - 2.0_qp*real(i,qp)*dm1(0)
        do p = 1, min(i+1, ne)
          dp1(p) = 2.0_qp*sqmu*(d0(p-1) - t*d0(p)) &
                 - 2.0_qp*real(i,qp)*dm1(p)
        enddo
        dm1 = d0
        d0  = dp1
      enddo
      d(0:ne) = d0(0:ne)
      end subroutine hermicart_q
!
!-----------------------------------------------------------------------
!     E_t^{ij} by the Hermite recursion (twin of etijcalc.f)
!-----------------------------------------------------------------------
      subroutine etijcalc_q (m, lamx, la, lb, ce, a, b, ax, bx)
      integer, intent(in) :: m, lamx, la, lb
      real(qp), intent(out) :: ce(-1:2*lamx,-1:lamx,-1:lamx,3)
      real(qp), intent(in) :: a, b, ax, bx
      integer :: i, j, t
      real(qp) :: p, pa, pb, tp
      ce(:,:,:,m) = 0.0_qp
      ce(0,0,0,m) = 1.0_qp
      if (la+lb .eq. 0) return
      p  = a + b
      pa = -b*(ax-bx)/p
      pb =  a*(ax-bx)/p
      tp = 1.0_qp/(2.0_qp*p)
      do i = 0, la
        do j = 0, lb
          do t = 1, i+j
            ce(t,i,j,m) = tp*(i*ce(t-1,i-1,j,m) + j*ce(t-1,i,j-1,m)) &
                          / real(t,qp)
          enddo
          if (i.lt.la) ce(0,i+1,j,m) = pa*ce(0,i,j,m) + ce(1,i,j,m)
          if (j.lt.lb) ce(0,i,j+1,m) = pb*ce(0,i,j,m) + ce(1,i,j,m)
        enddo
      enddo
      end subroutine etijcalc_q
!
!-----------------------------------------------------------------------
!     I(n,l,mu,T) = INT_0^inf R^n B_l(R) dR for n-l = 2m, m >= 1:
!       sqrt(pi)/4 mu^{(l-n-1)/2} T^l (b)_{m-1} 1F1(1-m;b;-mu T^2),
!     b = l+3/2, evaluated in log space (all terms positive).
!     T = 0: (1/2) mu^{-(n+1)/2} Gamma((n+1)/2) for l = 0, else 0.
!     status: 0 ok, 3 zero by rule, 5 invalid (odd n-l or m < 1).
!-----------------------------------------------------------------------
      subroutine radial_integral_q (n, l, mu, t, val, status)
      integer, intent(in) :: n, l
      real(qp), intent(in) :: mu, t
      real(qp), intent(out) :: val
      integer, intent(out) :: status
      real(qp) :: b, x, logx, lpre, lt, ls
      integer :: k, m
      val = 0.0_qp
      status = 0
      if (n.lt.0 .or. l.lt.0 .or. mu.le.0.0_qp) then
        status = 5
        return
      endif
      if (t .le. 0.0_qp) then
        if (l .eq. 0) then
          val = 0.5_qp * mu**(-real(n+1,qp)/2.0_qp) &
                * gamma(real(n+1,qp)/2.0_qp)
        else
          status = 3
        endif
        return
      endif
      if (mod(n-l,2).ne.0 .or. n-l.lt.2) then
        status = 5
        return
      endif
      m = (n-l)/2
      b = real(l,qp) + 1.5_qp
      x = mu*t*t
      logx = log(x)
      lpre = 0.5_qp*log(QPI) - 2.0_qp*log(2.0_qp) &
           + real(l-n-1,qp)/2.0_qp*log(mu) + real(l,qp)*log(t)
      do k = 0, m-2
        lpre = lpre + log(b + real(k,qp))
      enddo
      lt = 0.0_qp
      ls = 0.0_qp
      do k = 1, m-1
        lt = lt + log(real(m-k,qp)) - log(b+real(k-1,qp)) &
                + logx - log(real(k,qp))
        if (ls .gt. lt) then
          ls = ls + log(1.0_qp + exp(lt-ls))
        else
          ls = lt + log(1.0_qp + exp(ls-lt))
        endif
      enddo
      val = exp(lpre + ls)
      end subroutine radial_integral_q
!
!-----------------------------------------------------------------------
!     One shell quartet (pq|rs) of the norm pass in binary128.
!     ang?(3,d?) are the Cartesian exponents of the components of each
!     shell; abcd(di,dj,dk,dl) returns INT_0^inf I_abcd(R) dR / (4 pi)
!     exactly as intra_pairnorm (radial_integ.f) does in binary64.
!-----------------------------------------------------------------------
      subroutine quartet_norm_q (pq, rs, di, dj, dk, dl, angp, angq, &
                                 angr, angs, lmaxp, abcd)
      integer, intent(in) :: pq, rs, di, dj, dk, dl, lmaxp
      integer, intent(in) :: angp(3,di), angq(3,dj), angr(3,dk), angs(3,dl)
      real(8), intent(out) :: abcd(di,dj,dk,dl)
      real(qp) :: tvec(3), mu, sqrtmu, piz, tis, thetat, phit
      real(qp) :: coefe1, coefe2, eedpi, pw, tm, rad, iabcd
      real(qp), allocatable :: aac(:,:,:,:,:,:)
      real(qp) :: jtab(0:4*lmaxp,0:4*lmaxp,0:4*lmaxp)
      real(qp) :: dd(0:4*lmaxp), rsh(0:(12*lmaxp+1)**2-1)
      integer  :: maxterm(di,dj,dk,dl,3)
      integer  :: ka, kb, kc, kd, k, i, j, ip, jp, it, itp, iv
      integer  :: nx, ny, nz, nmx, nmy, nmz, lq, lup, l, m, nn, st
      real(qp), parameter :: cut = 1.0e-30_qp

      allocate (aac(di,dj,dk,dl,0:4*lmaxp,3))
      tvec(:) = pbig_q(rs,:) - pbig_q(pq,:)
      mu  = psmall_q(pq)*psmall_q(rs)/(psmall_q(pq)+psmall_q(rs))
      sqrtmu = sqrt(mu)
      piz = sqrt(QPI/(psmall_q(pq)+psmall_q(rs)))
      tis = sqrt(sum(tvec*tvec))
!
!-----A_n^{abcd} coefficients, per Cartesian component
!
      aac = 0.0_qp
      do k = 1, 3
        do ka = 1, di
          i = angp(k,ka)
          do kb = 1, dj
            j = angq(k,kb)
            do kc = 1, dk
              ip = angr(k,kc)
              do kd = 1, dl
                jp = angs(k,kd)
                maxterm(ka,kb,kc,kd,k) = i+j+ip+jp
                do it = 0, i+j
                  pw = 1.0_qp
                  if (mod(it,2).eq.1) pw = -1.0_qp
                  coefe1 = etij_q(pq,it,i,j,k)*prefac_q(pq,k)
                  if (coefe1 .eq. 0.0_qp) cycle
                  do itp = 0, ip+jp
                    coefe2 = etij_q(rs,itp,ip,jp,k)*prefac_q(rs,k)
                    if (coefe2 .eq. 0.0_qp) cycle
                    eedpi = coefe1*coefe2*piz*sqrtmu**(it+itp)*pw
                    call hermicart_q (mu, tvec(k), it+itp, 4*lmaxp, dd)
                    do iv = 0, it+itp
                      aac(ka,kb,kc,kd,iv,k) = aac(ka,kb,kc,kd,iv,k) &
                                            + eedpi*dd(iv)
                    enddo
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
      nmx = maxval(maxterm(:,:,:,:,1))
      nmy = maxval(maxterm(:,:,:,:,2))
      nmz = maxval(maxterm(:,:,:,:,3))
      lq  = nmx + nmy + nmz
      if (tis .gt. 1.0e-12_qp) then
        thetat = acos(max(-1.0_qp, min(1.0_qp, tvec(3)/tis)))
        phit   = atan2(tvec(2), tvec(1))
      else
        thetat = 0.0_qp
        phit   = 0.0_qp
      endif
      call rsh_q (lq, thetat, phit, rsh)
!
!-----Jtab(nx,ny,nz) = SUM_l I(n,l) SUM_m S_lm(T^) G_lm(nx,ny,nz)
!
      jtab = 0.0_qp
      do nx = 0, nmx
        do ny = 0, nmy
          do nz = 0, nmz
            lup = nx + ny + nz
            nn  = lup + 2
            tm  = 0.0_qp
            do l = 0, lup
              if (mod(lup-l,2) .ne. 0) cycle
              if (tis .le. 1.0e-12_qp .and. l .gt. 0) exit
              call radial_integral_q (nn, l, mu, tis, rad, st)
              if (st .ne. 0 .or. rad .eq. 0.0_qp) cycle
              do m = -l, l
                tm = tm + rsh(l*(l+1)+m)*gtab_q(l,m,nx,ny,nz)*rad
              enddo
            enddo
            jtab(nx,ny,nz) = tm
          enddo
        enddo
      enddo
!
!-----Contraction
!
      do ka = 1, di
        do kb = 1, dj
          do kc = 1, dk
            do kd = 1, dl
              iabcd = 0.0_qp
              do nx = 0, maxterm(ka,kb,kc,kd,1)
                if (abs(aac(ka,kb,kc,kd,nx,1)) .lt. cut) cycle
                do ny = 0, maxterm(ka,kb,kc,kd,2)
                  if (abs(aac(ka,kb,kc,kd,ny,2)) .lt. cut) cycle
                  do nz = 0, maxterm(ka,kb,kc,kd,3)
                    if (abs(aac(ka,kb,kc,kd,nz,3)) .lt. cut) cycle
                    iabcd = iabcd + aac(ka,kb,kc,kd,nx,1) &
                          * aac(ka,kb,kc,kd,ny,2)*aac(ka,kb,kc,kd,nz,3) &
                          * jtab(nx,ny,nz)
                  enddo
                enddo
              enddo
              abcd(ka,kb,kc,kd) = real(iabcd,8)
            enddo
          enddo
        enddo
      enddo
      deallocate (aac)
      end subroutine quartet_norm_q

      end module binary128_kernels
