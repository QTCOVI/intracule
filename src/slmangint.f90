!-----------------------------------------------------------------------
!  slmangint.f90
!
!  Angular integrals of a real spherical harmonic against the intracule
!  monomial weight, for the Bessel / spherical-harmonic replacement of the
!  Taylor expansion of  exp[2 mu R.T]  (see next.txt and lebedev.pdf).
!
!  For the code's angular convention  R^ = (sinT sinP, sinT cosP, cosT)
!  (i.e. the x axis carries sin(phi), the y axis carries cos(phi) -- this
!  matches how newtayloric.f indexes yntphi(nx,ny) and thetaval(nx+ny,nz))
!  this module returns
!
!     G_lm(nx,ny,nz) = INT_0^pi INT_0^2pi  S_lm(theta,phi)
!                       * sin^{nx+ny+1}theta cos^{nz}theta
!                       * sin^{nx}phi cos^{ny}phi  dtheta dphi
!
!  With the orthonormal real spherical harmonics S_lm and the addition
!  theorem  SUM_m S_lm(R^) S_lm(T^) = (2l+1)/(4pi) P_l(cos gamma), the
!  spherically averaged coupling table of newtayloric.f becomes the SHORT,
!  FINITE sum (note G_lm = 0 for l > nx+ny+nz):
!
!     J(nx,ny,nz) = SUM_{l=0}^{nx+ny+nz} B_l SUM_{m=-l}^{l} S_lm(T^) G_lm
!     B_l         = exp[-mu(|R|-|T|)^2] * i_l,scaled(2 mu |R| |T|)
!
!  Both B_l and S_lm(T^) are supplied by the caller's existing routines.
!  This module implements ONLY the geometry-independent G_lm, exactly as
!  requested, reusing the existing yntphi / thetaval tables.
!
!  NORMALISATION / PHASE CONVENTION (must match the caller's S_lm(T^)!):
!    Orthonormal real spherical harmonics, NO Condon-Shortley phase,
!    m>0 -> cos(m phi), m<0 -> sin(|m| phi):
!       m =0 : S_l0  = sqrt[(2l+1)/(4pi)]              P_l(cosT)
!       m>0 : S_lm  = sqrt[(2l+1)/(2pi) (l-m)!/(l+m)!] Pl^m(cosT) cos(m phi)
!       m<0 : S_lm  = sqrt[(2l+1)/(2pi) (l-|m|)!/(l+|m|)!] Pl^|m|(cosT) sin(|m| phi)
!    with  Pl^m(u) = (1-u^2)^{m/2} d^m/du^m P_l(u)     (no (-1)^m).
!    A per-(l,m) sign difference from the caller's routine would break the
!    addition theorem: if the parent S_lm uses the Condon-Shortley phase or
!    a different real form, flip the single marked line in slm_setnorm.
!-----------------------------------------------------------------------
      module space_for_slmangint
      implicit none
      save
      integer :: lmax_slm = -1
!     pcoef(l,c)   : coefficient of u^c in the Legendre polynomial P_l(u)
!     qcoef(l,m,c) : coefficient of u^c in d^m/du^m P_l(u)
!     anlm (l,m)   : real-SH normalisation N_l^m   (m = |m| >= 0)
      real(kind=8), allocatable :: pcoef(:,:)
      real(kind=8), allocatable :: qcoef(:,:,:)
      real(kind=8), allocatable :: anlm(:,:)
      end module space_for_slmangint
!
!-----------------------------------------------------------------------
!     Build the Legendre-derivative coefficient tables and the real-SH
!     normalisations up to degree lmax.  Call ONCE (like phiintegrals /
!     thetaintegrals), before any call to slmangint / build_gtab.
!-----------------------------------------------------------------------
!
      subroutine init_slmangint (lmax)
      use space_for_slmangint
      implicit none
      integer, intent(in) :: lmax
      integer :: l, m, c
      real(kind=8) :: fac, pref, ratio
      real(kind=8), parameter :: pi = 3.14159265358979323846d0

      if (lmax_slm .ge. lmax) return          ! already built large enough
      if (allocated(pcoef)) deallocate (pcoef)
      if (allocated(qcoef)) deallocate (qcoef)
      if (allocated(anlm))  deallocate (anlm)
      allocate (pcoef(0:lmax,0:lmax))
      allocate (qcoef(0:lmax,0:lmax,0:lmax))
      allocate (anlm (0:lmax,0:lmax))
      pcoef = 0.0d0
      qcoef = 0.0d0
      anlm  = 0.0d0
!
!-----Legendre polynomials by the standard recurrence
!        l P_l(u) = (2l-1) u P_{l-1}(u) - (l-1) P_{l-2}(u)
!     stored as coefficients of u^c (only c with parity of l are nonzero).
!
      pcoef(0,0) = 1.0d0
      if (lmax .ge. 1) pcoef(1,1) = 1.0d0
      do l = 2, lmax
        do c = 0, l
          pref = 0.0d0
          if (c .ge. 1) pref = dble(2*l-1) * pcoef(l-1,c-1)
          pcoef(l,c) = ( pref - dble(l-1)*pcoef(l-2,c) ) / dble(l)
        enddo
      enddo
!
!-----q^{l,m}(u) = d^m/du^m P_l(u).  Coefficient of u^c is
!        pcoef(l,c+m) * (c+m)!/c!  =  pcoef(l,c+m) * PROD_{j=1}^{m}(c+j)
!
      do l = 0, lmax
        do m = 0, l
          do c = 0, l-m
            fac = 1.0d0                       ! product (c+1)(c+2)...(c+m)
            block
              integer :: j
              do j = 1, m
                fac = fac * dble(c+j)
              enddo
            end block
            qcoef(l,m,c) = pcoef(l,c+m) * fac
          enddo
        enddo
      enddo
!
!-----real-SH normalisations
!        m =0 : sqrt[(2l+1)/(4pi)]
!        m>0 : sqrt[(2l+1)/(2pi) * (l-m)!/(l+m)!]
!     the factorial ratio is built as PROD_{j=l-m+1}^{l+m} 1/j  (overflow safe)
!
      do l = 0, lmax
        anlm(l,0) = sqrt( dble(2*l+1)/(4.0d0*pi) )
        do m = 1, l
          ratio = 1.0d0
          block
            integer :: j
            do j = l-m+1, l+m
              ratio = ratio / dble(j)
            enddo
          end block
          anlm(l,m) = sqrt( dble(2*l+1)/(2.0d0*pi) * ratio )
!         ---- Condon-Shortley knob (see header): to ADD the (-1)^m phase,
!         ---- uncomment the next line so it matches the parent S_lm routine.
!         if (mod(m,2).eq.1) anlm(l,m) = -anlm(l,m)
        enddo
      enddo
      lmax_slm = lmax
      end subroutine init_slmangint
!
!-----------------------------------------------------------------------
!     One angular integral  G_lm(nx,ny,nz)  (see module header).
!     yntphi(0:np,0:np)   : INT_0^2pi sin^a phi cos^b phi dphi = yntphi(a,b)
!     thetaval(0:nt,0:nt) : INT_0^pi sin^{p+1}theta cos^q theta = thetaval(p,q)
!     The caller must size np >= |m|+nx  and (|m|+ny) and  nt >= |m|+nx+ny
!     and >= (l+nz).  init_slmangint(>=l) must have run.
!-----------------------------------------------------------------------
!
      function slmangint (l,m,nx,ny,nz,np,nt,yntphi,thetaval) result (g)
      use space_for_slmangint
      implicit none
      integer, intent(in) :: l,m,nx,ny,nz,np,nt
      real(kind=8), intent(in) :: yntphi(0:np,0:np)
      real(kind=8), intent(in) :: thetaval(0:nt,0:nt)
      real(kind=8) :: g
      integer :: mm, k, kstart, c, aphi, bphi, pth
      real(kind=8) :: gk, phipart, thetasum

      g = 0.0d0
      mm = abs(m)
      if (l .gt. mm+nx+ny+nz) return       ! (also caught by zero tables)
      pth = mm + nx + ny                    ! sin-power index of thetaval
      if (pth .gt. nt) return
!
!-----Fourier expansion of cos(mm*phi) (m>=0) or sin(mm*phi) (m<0):
!        cos(mm p) = SUM_{k=0,2,..} C(mm,k)(-1)^{k/2}  cos^{mm-k} sin^k
!        sin(mm p) = SUM_{k=1,3,..} C(mm,k)(-1)^{(k-1)/2} cos^{mm-k} sin^k
!     Under R^=(sinT sinP,sinT cosP,cosT) the x-power (mm-k) becomes a
!     sin(phi) power and the y-power (k) a cos(phi) power, hence:
!        aphi = (mm-k)+nx   (sin phi),   bphi = k+ny   (cos phi)
!
      if (m .ge. 0) then
        kstart = 0
      else
        kstart = 1
      endif
      do k = kstart, mm, 2
        aphi = (mm-k) + nx
        bphi = k + ny
        if (aphi .gt. np .or. bphi .gt. np) cycle
        phipart = yntphi(aphi,bphi)
        if (phipart .eq. 0.0d0) cycle            ! parity: both must be even
        gk = binomd(mm,k)
        if (mod(k,2).eq.0) then
          if (mod(k/2,2).eq.1)     gk = -gk       ! (-1)^{k/2}
        else
          if (mod((k-1)/2,2).eq.1) gk = -gk       ! (-1)^{(k-1)/2}
        endif
!
!-------theta part: SUM_c q^{l,mm}_c * thetaval(pth, c+nz)   (c same parity l-mm)
!
        thetasum = 0.0d0
        do c = l-mm, 0, -2
          if (c+nz .gt. nt) cycle
          thetasum = thetasum + qcoef(l,mm,c) * thetaval(pth, c+nz)
        enddo
        g = g + gk * phipart * thetasum
      enddo
      g = anlm(l,mm) * g

      contains
        real(kind=8) function binomd (n,i)
        integer, intent(in) :: n,i
        integer :: j
        binomd = 1.0d0
        do j = 1, i
          binomd = binomd * dble(n-j+1) / dble(j)
        enddo
        end function binomd
      end function slmangint
!
!-----------------------------------------------------------------------
!     Fill the full table  gtab(l,m,nx,ny,nz)  once, for all
!     0<=l<=lmax, -l<=m<=l, 0<=nx<=nmx, 0<=ny<=nmy, 0<=nz<=nmz.
!     lmax should be nmx+nmy+nmz (higher l give exactly zero).
!     This is the geometry-independent object reused for every R and T.
!-----------------------------------------------------------------------
!
      subroutine build_gtab (gtab,lmax,nmx,nmy,nmz,np,nt,yntphi,thetaval)
      use space_for_slmangint
      implicit none
      integer, intent(in) :: lmax,nmx,nmy,nmz,np,nt
      real(kind=8), intent(in)  :: yntphi(0:np,0:np)
      real(kind=8), intent(in)  :: thetaval(0:nt,0:nt)
      real(kind=8), intent(out) :: gtab(0:lmax,-lmax:lmax,0:nmx,0:nmy,0:nmz)
      integer :: l,m,nx,ny,nz
      real(kind=8), external :: slmangint

      call init_slmangint (lmax)
      gtab = 0.0d0
      !$omp parallel do collapse(3) default(shared) &
      !$omp   private(nx,ny,nz,l,m) schedule(dynamic)
      do nx = 0, nmx
        do ny = 0, nmy
          do nz = 0, nmz
            do l = 0, min(lmax, nx+ny+nz)
              do m = -l, l
                gtab(l,m,nx,ny,nz) = &
                  slmangint(l,m,nx,ny,nz,np,nt,yntphi,thetaval)
              enddo
            enddo
          enddo
        enddo
      enddo
      !$omp end parallel do
      end subroutine build_gtab
