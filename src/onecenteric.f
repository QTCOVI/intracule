
!-----------------------------------------------------------------------
!
      subroutine onecenteric (verbose,wfnfile,rval,xicMO,xicMOCoul,
     &  xicXC,xicMODF,nrval,kmaxcap,thresh)
!
!------------------------------------------------------------------------
!
      USE space_for_wfnbasis
      USE space_for_wfncoef
      USE space_for_primgto
      USE space_for_othershells
      USE space_for_inexcule
      USE mod_rdm
      implicit real(kind=8) (a-h,o-z)
      include   'param.inc'
      include   'wfn.inc'
      include   'corr.inc'
      include   'io.inc'
      include   'primgto.inc'
!
!-----Input variables
!
      logical,           intent(in)  :: verbose
      character (len=*), intent(in)  :: wfnfile
      integer (kind=4),  intent(in)  :: nrval
      real (kind=8),     intent(in)  :: rval(nrval)
      real (kind=8),     intent(out) :: xicMO (nrval)
      real (kind=8),     intent(out) :: xicMOCoul(nrval)
      real (kind=8),     intent(out) :: xicXC(nrval)
      real (kind=8),     intent(out) :: xicMODF(nrval)
      real (kind=8),     intent(inout)  :: thresh
      integer (kind=4),  intent(in)  :: kmaxcap

!
!-----Shared data.  The 2-RDMs and the angular-integral tables are built
!     once in the host below and only read inside the parallel quartet loop;
!     aoicule is shared but written in disjoint blocks (one (p,q) bra-block
!     per task); all per-pair scratch is declared inside intra_pair.
!
      real(kind=8), allocatable,dimension (:,:,:,:)     :: rdm2coul
      real(kind=8), allocatable,dimension (:,:,:,:)     :: rdm2df
      real(kind=8), allocatable,dimension (:,:,:,:)     :: aoicule
      real(kind=8), allocatable,dimension (:,:,:,:)     :: moic
      real(kind=8), allocatable,dimension (:,:)         :: thetaval
!
!-----Shell-offset table and (p,q) bra-pair map for the parallel quartet loop
!
      integer(kind=4), allocatable,dimension (:) :: shoff,pmap,qmap
!
!-----Angular integrals over phi
!
      real(kind=8), allocatable,dimension (:,:) :: yntphi

      integer (kind=4)  p,q,r,s
!
!-----array to store the coefficiente of the Hermite polynomials up to n=20
!
      integer(kind=4), parameter :: maxn = 20
      real   (kind=8) Chermite(0:maxn,0:maxn)
      real   (kind=8) logfac(0:kmaxcap)

      write (uout,'(a)')' # ++++++++++++++++++++++++++++++++++++++++++'
      write (uout,'(a,I8)')' # KMAX in Taylor expansion = ',kmaxcap
      write (uout,'(a,1PE15.8)')' # Requested convergence = ',thresh
      write (uout,'(a)')' # ++++++++++++++++++++++++++++++++++++++++++'
      if (maxn.lt.4*maxlprim) then
        write (uerr,1) maxn,4*maxlprim
        stop
      else
        ixg = ngtox()
        Chermite(0:maxn,0:maxn) = 0.0d+00
        call hermitec (Chermite,maxn)
      endif
!
!-----Compute 2-RDM
!
      write (uout,50) 'Canonical Orbitals'
      call timer (2,irdm2,'_rdm2MOs____________',-1)
      call compute_rdm (ndets,nelact,ncore,nmo,uout,uerr,nrecwfn, 
     &    cciqa,mp2iqa,icorr,verbose,wfnfile)
      call timer (4,irdm2,'_rdm2MOs____________',-1)
!
!-----Compute the Coulombic part of the RDM2 and the Dirac-Fock 
!     part of the XC RDM2
!
      allocate (rdm2coul(nmo,nmo,nmo,nmo))
      allocate (rdm2df(nmo,nmo,nmo,nmo))
      rdm2coul = 0.0d+00
      rdm2df = 0.0d+00
      do p=1,nmo
        do q=1,nmo
          do r=1,nmo
            do s=1,nmo
              rdm2coul(p,q,r,s) = rdm1(p,q) * rdm1(r,s)
              rdm2df(p,s,q,r) = rdm1(p,q) * rdm1(r,s)
            enddo
          enddo
        enddo
      enddo
!
!-----Classification of GTFs in shells
!
      call timer (2,ionecenteric,'_onecenteric________',-1)
      nshells = sum(ngroup(1:ncent))
      call allocate_space_for_othershells (nshells,ngtoG)
      write (uout,'(a)')' # Classification of CGFs into shells'
      i = 0
      do ic=1,ncent
        do m=1,ngroup(ic)
          i = i+1
          center(i) = ic
          nzshell(i) = nzexp(ic,m)
          ii = nuexp(ic,m,1)
          lang(i) = nlm(ityp(ii),1)+nlm(ityp(ii),2)+nlm(ityp(ii),3)
          zeta(i) = oexp(ii)
          ordp(i,1:nzshell(i)) = nuexp(ic,m,1:nzexp(ic,m))
        enddo
      enddo 
!
!-----We compute the INT_0^2\pi sin(x)^(2k) cos(x)^(2l) dx integrals
!
      write (uout,'(a)')' # Computing Phi(nx,ny) angular integrals'
      nxmaxphi = 4*maxlprim + kmaxcap + 2
      allocate (yntphi(0:nxmaxphi,0:nxmaxphi))
C-----call init_yntphi_recursive (yntphi,nxmaxphi)
      call phiintegrals (nxmaxphi,nxmaxphi,yntphi)
!
!-----------------------------------------------------------------------
!
!-----Compute the E_t^ij, E_u^kl and E_v^mn coefficients
!
      write (uout,51)
      call computeEcoef ()

      write (uout,'(a)')' # Computing Theta(nx,ny,nz) angular integrals'
      nxmaxtheta = 4*maxlprim + kmaxcap
      nxy_max = 2 * nxmaxtheta   ! porque nx+ny puede llegar hasta 2*nxmaxtheta
      allocate(thetaval(0:nxy_max, 0:nxmaxtheta+2))
C-----call init_thetaval_recursive (thetaval,nxy_max,nxmaxtheta+2)
      call thetaintegrals (nxy_max,nxmaxtheta+2,thetaval)
!
!-----Table of log(k!) for the overflow-safe Taylor summation below
!
      write (uout,'(a,I8)')' # Obtain Logs of integers up to ',kmaxcap
      logfac(0) = 0.0d+00
      do k = 1, kmaxcap
        logfac(k) = logfac(k-1) + log(dble(k))
      enddo
!
!-----Print the r values of the grid
!
      write (uout,20)
      write (uout,21) (rval(i),i=1,nrval)
      write (uout,20)
!
!-----aoicule(a,b,c,d) is the spherically averaged I_abcd(R) intracule density
!     between CGTFs.  It is shared across the parallel bra-pair loop (each
!     task writes a disjoint block); the per-pair scratch is private and
!     lives inside intra_pair.
!
      write (uout,9)
!
!-----Precompute the shell-offset table (the AO index at which each shell
!     starts) and the (p,q) bra-pair map, so the parallel quartet loop below
!     can recover all of its indices from a single flat task id, with no
!     sequential counters to serialise.
!
      npairs = nshells*(nshells+1)/2
      allocate (shoff(nshells))
      allocate (pmap(npairs))
      allocate (qmap(npairs))
      shoff(1) = 0
      do p = 2,nshells
        shoff(p) = shoff(p-1) + nzshell(p-1)
      enddo
      ipq = 0
      do p = 1,nshells
        do q = 1,p
          ipq = ipq + 1
          pmap(ipq) = p
          qmap(ipq) = q
        enddo
      enddo
      allocate (aoicule(nprims,nprims,nprims,nprims))
      allocate (moic(nmo,nmo,nmo,nmo))
!
!-----Loop over R values (serial).  The work is dominated by the shell-
!     quartet integrals, which are parallelised over the bra-pair index
!     ipq -> (p,q): each task scatters its quartets into a DISJOINT (p,q)
!     block of the shared aoicule (the within-pair symmetry only permutes
!     indices inside that block), so the writes are race-free.  Per-pair
!     scratch lives in intra_pair and is thread-private.  R is kept serial
!     to avoid nested parallelism; the quartet loop already exposes
!     O(nshells^2) tasks, so this scales even for a single R point.
!
      do irval=1,nrval
        rr = rval(irval)
        r2 = rr*rr
        aoicule = 0.0d+00
!$omp parallel do schedule(dynamic) default(shared) private(ipq)
        do ipq = 1,npairs
          call intra_pair (ipq,rr,r2,aoicule)
        enddo
!$omp end parallel do
        call moicule (aoicule,moic,nprims,nmo,uout,verbose)
        xicMO(irval) = dot4d_kahan(rdm2,moic,nmo,nmo,nmo,nmo)
        xicMOCoul(irval) = dot4d_kahan(rdm2coul,moic,nmo,nmo,nmo,nmo)
        xicXC(irval) = xicMO(irval) - xicMOCoul(irval)
        xicMODF(irval) =-dot4d_kahan(rdm2df,moic,nmo,nmo,nmo,nmo)/2.0d0
        write (uout,100) rval(irval),xicMO(irval),xicMOCoul(irval),
     &    xicXC(irval),xicMODF(irval)
      enddo
      deallocate (aoicule)
      deallocate (moic)
      deallocate (pmap,qmap,shoff)
      deallocate (rdm2coul)
      deallocate (rdm2df)
      call timer (4,ionecenteric,'_onecenteric________',-1)
      return
!
!-----Formats (host scope; referenced by the writes above)
!
 1    format (/,' #',/,
     &  ' # The routine hermitec.f must be changed increasing the',/,
     &  ' # current value of maxn = ',I2,' up to at least a value of ',
     & I4,/,' #')
 9    format (//' # Spherical Intracule density',1x,
     & '(Tot <- 2RDM, Coul <- 2RDM Coul, XC <- 2RDM(XC) ,',
     & 'XC-DF <- 2RDM(XC-Dirac-Fock))',/,' # ',
     & 109('-'),/,' # ',6x,'R',14x,'IDensity (Tot)',9x,
     & 'IDensity (Coul)',10x,'IDensity (XC)',10x,'Density (XC-DF)',/,
     & ' # ',109('-'))
 20   format (' # ',80('+'))
 21   format (' # Spherically averaged intracule density at R = ',/,
     &   ' # ',80('+'),/,1000(5(1x,F15.8),/))
 50   format (' # Computing 1- and 2-RDM in the basis of ',a)
 51   format (' # Computing E_ij^t coefficients')
 100  format (1x,F15.8,8(2x,E22.15))

      contains
!
!-----------------------------------------------------------------------
!     Per bra-pair worker for one R.  For bra-pair ipq -> (p,q) it loops
!     over all ket pairs (r,s), evaluates the quartet intracule integrals
!     and scatters them (with the (pq|rs) within-pair permutational
!     symmetry) into the shared aoicule.  Every array declared here is
!     local, hence thread-private; aoicule is shared but each ipq writes a
!     disjoint (p,q) block.  The read-only tables (etij, prefac, pbig,
!     psmall, yntphi, thetaval, Chermite, logfac, shell data) come in by
!     host association.
!-----------------------------------------------------------------------
!
      subroutine intra_pair (ipq,rr,r2,aoicule)
      implicit real(kind=8) (a-h,o-z)
      integer (kind=4), intent(in) :: ipq
      real    (kind=8), intent(in) :: rr,r2
      real    (kind=8) aoicule(nprims,nprims,nprims,nprims)
      integer (kind=4) p,q,pq,r,s,rs,di,dj,dk,dl,idx,jdx,kdx,ldx
      integer (kind=4) labcd,i,j,k,ka,kb,kc,kd,ip,jp,it,itp,ival,kk
      integer (kind=4) nx,ny,nz,nxk,nyk,nzk,kx,ky,kz,k0,kup,kwmax
      integer (kind=4) nmaxx,nmaxy,nmaxz,ntermx,ntermy,ntermz
      integer (kind=4) i1,i2,j1,j2,k1,k2,l1,l2,ia,ja,la
      integer (kind=4) lab,lcd
      real    (kind=8) mu,Iabcd
      real    (kind=8) wx(0:kmaxcap),wy(0:kmaxcap),wz(0:kmaxcap)
      real    (kind=8) dd(0:maxn)
      real(kind=8),allocatable,dimension(:,:,:,:)     :: abcd
      real(kind=8),allocatable,dimension(:,:,:,:,:,:) :: aac
      real(kind=8),allocatable,dimension(:,:,:)       :: Jtab
      integer(kind=4),allocatable,dimension(:,:,:,:,:):: maxterm

      allocate (abcd(ixg,ixg,ixg,ixg))
      allocate (aac(ixg,ixg,ixg,ixg,0:4*maxlprim,1:3))
      allocate (maxterm(ixg,ixg,ixg,ixg,1:3))
      allocate (Jtab(0:4*maxlprim,0:4*maxlprim,0:4*maxlprim))
!
!-----recover the bra-pair indices and AO offsets for this task
!
      p   = pmap(ipq)
      q   = qmap(ipq)
      pq  = ipq
      di  = nzshell(p)
      dj  = nzshell(q)
      idx = shoff(p)
      jdx = shoff(q)
      pf1 = prefac(pq,1)*prefac(pq,2)*prefac(pq,3)
      lab = lang(p) + lang(q)
      do r = 1,nshells
        dk = nzshell(r)
        kdx = shoff(r)
        do s = 1,r
          rs = r*(r-1)/2 + s
          dl = nzshell(s)
          ldx = shoff(s)
          lcd  = lang(r) + lang(s)
          prodz = psmall(pq)*psmall(rs)
          addz = psmall(pq)+psmall(rs)
          mu = prodz/addz
          twomurr = 2.0d+00 * mu * rr
          sqrtmu = sqrt(mu)
          piz = sqrt(pi/addz)
!
!---------this is the quantity T of the 'nintracule.pdf' document
!
          pf2 = prefac(rs,1)*prefac(rs,2)*prefac(rs,3)
          xk = pf1 * pf2 * piz**3
          xmur2 = mu*r2
          i1 = idx + 1
          i2 = idx + di
          j1 = jdx + 1
          j2 = jdx + dj
          k1 = kdx + 1
          k2 = kdx + dk
          l1 = ldx + 1
          l2 = ldx + dl
          if (xmur2 .gt. 1000.0D+00) then
            abcd(1:di,1:dj,1:dk,1:dl) = 0.0d+00
            goto 222
          endif
          xdamp = -xmur2 - log(4.0d+00*pi)
!
!---------Sum of angular momenta of the four shells of CGs
!
          labcd = lab + lcd
 
!---------aac = The A_nx^{abcd},  A_ny^{abcd},  A_nz^{abcd} 
!         quantities of the 'nuevointracule.pdf' document
!
          aac(:,:,:,:,:,:) = 0.0d+00
          do k=1,3 ! the three components x,y,z
            do ka=1,di
              i = nlm(ityp(ordp(p,ka)),k)
              do kb=1,dj
                j = nlm(ityp(ordp(q,kb)),k)
                do kc=1,dk
                  ip = nlm(ityp(ordp(r,kc)),k)
                  do kd=1,dl
                    jp = nlm(ityp(ordp(s,kd)),k)
                    maxterm(ka,kb,kc,kd,k) = i+j+ip+jp
                    do it=0,i+j
                      pow=(-1.0d+00)**it
                      coefe1 = etij(pq,it,i,j,k)*prefac(pq,k)
                      if (abs(coefe1).lt.1.0d-14) cycle
                      do itp=0,ip+jp
                        coefe2 = etij(rs,itp,ip,jp,k)*prefac(rs,k)
                        if (abs(coefe2).lt.1.0d-14) cycle
                        EtijEtij = coefe1 * coefe2
                        xmupow = sqrtmu**(it+itp) 
                        eedpi = EtijEtij * piz * xmupow * pow
                        if (abs(eedpi).lt.1.0d-14) cycle
                        call hermicart (mu,0.0d+00,it+itp,maxn,dd)
                        do ival=0,it+itp
                          aac(ka,kb,kc,kd,ival,k) = eedpi * 
     &                    dd(ival) + aac(ka,kb,kc,kd,ival,k) 
                        enddo
                      enddo
                    enddo
                  enddo
                enddo
              enddo
            enddo
          enddo
!
!---------Maximum Cartesian-Hermite orders actually needed by any of
!         the di*dj*dk*dl Cartesian components of this quartet
!
          nmaxx = maxval(maxterm(1:di,1:dj,1:dk,1:dl,1))
          nmaxy = maxval(maxterm(1:di,1:dj,1:dk,1:dl,2))
          nmaxz = maxval(maxterm(1:di,1:dj,1:dk,1:dl,3))
!
!---------Damped Taylor weights (2*mu*R*T_i)^k/k! * exp(xdamp/3),
!         built in log space (so neither (2muRT)^k nor 1/k!
!         overflows) and reused for every (nx,ny,nz) and every
!         Cartesian component below, so the exp() calls leave the
!         innermost loops.  They are filled lazily up to the highest
!         order kwmax the adaptive loop below actually reaches (which
!         for s/p / loosely-coupled quartets is small), instead of
!         eagerly up to kmaxcap.
!
          xdamp3 = xdamp / 3.0d+00
          kwmax = -1
!
!---------Taylor expansion of the coupling factor e^{2*mu*R.T}
!         (Eq. 52 of fusion3.pdf).  The Gaussian parities fix
!         kx==nx, ky==ny, kz==nz (mod 2), so only total orders k of
!         the same parity as nx+ny+nz contribute; we step k by 2
!         from k0 and sum the multinomial terms p+q+r=k, increasing
!         k adaptively until the contribution falls below thresh.
!
!         The coupling integral depends only on (nx,ny,nz) and the
!         quartet geometry, NOT on the Cartesian indices, so it is
!         tabulated once in Jtab(nx,ny,nz) and reused below.
!
          do nx=0,nmaxx
            do ny=0,nmaxy
              do nz=0,nmaxz
                k0 = mod(nx+ny+nz,2)
                kup = k0
                tm = 0.0d+00
                do k = k0, kup, 2
!
!-----------------extend the weight tables on demand (each order is
!                 built exactly once per quartet)
!
                  if (k.gt.kwmax) then
                    do kk = kwmax+1, k
                      wx(kk) = exp(xdamp3 - logfac(kk))
                      wy(kk) = exp(xdamp3 - logfac(kk))
                      wz(kk) = exp(xdamp3 - logfac(kk))
                    enddo
                    kwmax = k
                  endif
                  contrib = 0.0d+00
                  do kx = 0, k
                    if (kx.gt.0) cycle
                    nxk = nx+kx
                    if (mod(nxk,2).ne.0) cycle
                    do ky = 0, k-kx
                      if (ky.gt.0) cycle
                      kz = k-kx-ky
                      if (kz.gt.0) cycle
                      nyk = ny+ky
                      nzk = nz+kz
                      if (mod(nyk,2).ne.0) cycle
                      if (mod(nzk,2).ne.0) cycle
                      thetaint = thetaval(nxk+nyk,nzk)
                      if (thetaint.eq.0.0d+00) cycle
                      phiint = yntphi(nxk,nyk)
                      if (phiint.eq.0.0d+00) cycle
                      wxyz = wx(kx)*wy(ky)*wz(kz)
                      contrib = contrib + wxyz * thetaint * phiint
                    enddo  ! ky
                  enddo  ! kx
                  tm = tm + contrib
                  acontrib = abs(contrib)
                  atm = abs(tm)
                  if (k.ge.0 .and. acontrib.le.thresh*atm) exit
                enddo  ! k (adaptive Taylor order)
                Jtab(nx,ny,nz) = rr**(nx+ny+nz+2) * tm
                otm = Jtab(nx,ny,nz)
              enddo ! nz
            enddo ! ny
          enddo ! nx
!
!---------Evaluate the I_abcd(R) integral: contract the A_nx/A_ny/A_nz
!         coefficients with the precomputed coupling table Jtab
!
          do ka=1,di
            do kb=1,dj
              do kc=1,dk
                do kd=1,dl
                  ntermx = maxterm(ka,kb,kc,kd,1)
                  ntermy = maxterm(ka,kb,kc,kd,2)
                  ntermz = maxterm(ka,kb,kc,kd,3)
                  Iabcd = 0.0d+00
                  do nx=0,ntermx
                    axis = aac(ka,kb,kc,kd,nx,1)
                    if (abs(axis).lt.1.0d-08) cycle
                    do ny=0,ntermy
                      ayis = aac(ka,kb,kc,kd,ny,2)
                      if (abs(ayis).lt.1.0d-08) cycle
                      do nz=0,ntermz
                        azis = aac(ka,kb,kc,kd,nz,3)
                        if (abs(azis).lt.1.0d-08) cycle
                        if (Jtab(nx,ny,nz).eq.0.0d+00) cycle
                        axyz = axis*ayis*azis
                        Iabcd = Iabcd + axyz * Jtab(nx,ny,nz)
                      enddo ! nz
                    enddo ! ny
                  enddo  ! nx
                  abcd(ka,kb,kc,kd) = Iabcd
                enddo ! kd
              enddo ! kc
            enddo ! kb
          enddo ! ka
 222      continue
          aoicule(i1:i2,j1:j2,k1:k2,l1:l2) =
     &      abcd(1:di,1:dj,1:dk,1:dl)
          do ia=i1,i2
            do ja=j1,j2
              do ka=k1,k2
                do la=l1,l2
                  aoicule(ja,ia,ka,la) = aoicule(ia,ja,ka,la)
                  aoicule(ia,ja,la,ka) = aoicule(ia,ja,ka,la)
                  aoicule(ja,ia,la,ka) = aoicule(ia,ja,ka,la)
                enddo
              enddo
            enddo
          enddo
        end do !s
      end do !r
!
!-----Free this task's per-pair scratch
!
      deallocate (abcd)
      deallocate (aac)
      deallocate (Jtab)
      deallocate (maxterm)
      return
      end subroutine intra_pair
      end
