!
!-----------------------------------------------------------------------
!
      subroutine radial_integ (verbose,wfnfile,gsl,kummer,expkummer)
!
!------------------------------------------------------------------------
!
      USE space_for_wfnbasis
      USE space_for_wfncoef
      USE space_for_primgto
      USE space_for_othershells
      USE space_for_inexcule
      USE space_for_slmangint
      USE mod_rdm
      USE intracule_integral
      USE binary128_kernels
      USE kummer_pos, ONLY: kummer_1f1_asym, kummer_1f1_series,
     &    kummer_1f1_log, kummer_1f1_dispatch

      USE kummer_exp_neg, ONLY: kummer_G_dispatch, kummer_G_series, 
     &    kummer_G_asym
!
!-----Needed to use the bessel_scaled_gauss routines
!
      use, intrinsic :: iso_c_binding
      use gsl_bessel_interface
      use gaussian_scaled_bessel_mod
!     real(c_double), allocatable :: f_l(:), log_f_l(:)
!
!-----Needed to the real_spherical_harmonics routines
!
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use real_spherical_harmonics_mod

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

!
!-----Shared data.  The 2-RDMs and the angular-integral tables are built
!     once in the host below and only read inside the parallel quartet loop;
!     aoicule is shared but written in disjoint blocks (one (p,q) bra-block
!     per task); all per-pair scratch is declared inside intra_pair.
!
      real(kind=8), allocatable,dimension (:,:,:,:) :: rdm2coul
      real(kind=8), allocatable,dimension (:,:,:,:) :: rdm2diracfock
      real(kind=8), allocatable,dimension (:,:,:,:) :: aoicule
      real(kind=8), allocatable,dimension (:,:,:,:) :: moic
      real(kind=8), allocatable,dimension (:,:)     :: thetaval
      real(kind=8), allocatable,dimension (:,:)     :: yntphi
!
!-----Shell-offset table and (p,q) bra-pair map for the parallel quartet loop
!
      integer(kind=4), allocatable,dimension (:) :: shoff,pmap,qmap

      integer (kind=4)  p,q,r,s
!
!-----array to store the coefficiente of the Hermite polynomials up to n=20
!
      integer(kind=4), parameter :: maxn = 20
      real   (kind=8) Chermite(0:maxn,0:maxn)
!
!-----Geometry-independent angular integrals  gtab(l,m,nx,ny,nz) =
!     INT S_lm f dtheta dphi  (slmangint.f90).  Built ONCE below and shared
!     read-only across every R and every quartet.  nmaxc = 4*maxlprim is the
!     safe upper bound on nx,ny,nz (matches the aac index range); the l/m
!     ranges then run up to Lmaxg = 3*nmaxc = nx+ny+nz|max.
!
      real(kind=8), allocatable,dimension (:,:,:,:,:) :: gtab
      integer(kind=4) nmaxc,Lmaxg,ntab
 
      logical integIabcd
      logical gsl,kummer,expkummer

      if (maxn.lt.4*maxlprim) then
        write (uerr,1) maxn,4*maxlprim
        stop
      else
        ixg = ngtox()
        Chermite(0:maxn,0:maxn) = 0.0d+00
        call hermitec (Chermite,maxn)
      endif
      call compute_rdm (ndets,nelact,ncore,nmo,uout,uerr,nrecwfn, 
     &    cciqa,mp2iqa,icorr,verbose,wfnfile)
!
!-----Compute the Coulombic part of the RDM2 and the Dirac-Fock 
!     part of the XC RDM2
!
      allocate (rdm2coul(nmo,nmo,nmo,nmo))
      allocate (rdm2diracfock(nmo,nmo,nmo,nmo))
      rdm2coul = 0.0d+00
      rdm2diracfock = 0.0d+00
      do p=1,nmo
        do q=1,nmo
          do r=1,nmo
            do s=1,nmo
              rdm2coul(p,q,r,s) = rdm1(p,q) * rdm1(r,s)
              rdm2diracfock(p,s,q,r) = rdm1(p,q) * rdm1(r,s)
            enddo
          enddo
        enddo
      enddo
!
!-----Classification of GTFs in shells
!
      call timer (2,iradial_integ,'_radial_integ_______',-1)
      nshells = sum(ngroup(1:ncent))
      call allocate_space_for_othershells (nshells,ngtoG)
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
!-----Angular index bounds.  nx,ny,nz <= nmaxc = 4*maxlprim (aac range);
!     l,m <= Lmaxg = 3*nmaxc; the phi/theta power tables need indices up to
!     ntab = 2*Lmaxg (aphi=(mm-k)+nx, pth=mm+nx+ny, etc. -- see slmangint).
!
      nmaxc = 4*maxlprim
      Lmaxg = 3*nmaxc
      ntab  = 2*Lmaxg
!
!-----INT_0^2pi sin^a cos^b dphi  and  INT_0^pi sin^{p+1} cos^q dtheta,
!     now only up to ntab (was kmaxcap-sized in the Taylor path).
!
      allocate (yntphi(0:ntab,0:ntab))
      call phiintegrals (ntab,ntab,yntphi)
      allocate (thetaval(0:ntab,0:ntab))
      call thetaintegrals (ntab,ntab,thetaval)
!
!-----------------------------------------------------------------------
!
!-----Compute the E_t^ij, E_u^kl and E_v^mn coefficients
!
!     write (uout,51)
      call computeEcoef ()
!
!-----Build the geometry-independent angular integrals G_lm(nx,ny,nz) ONCE.
!     This is the object that replaces the per-quartet Taylor summation.
!
      allocate (gtab(0:Lmaxg,-Lmaxg:Lmaxg,0:nmaxc,0:nmaxc,0:nmaxc))
      call build_gtab (gtab,Lmaxg,nmaxc,nmaxc,nmaxc,ntab,ntab,
     &  yntphi,thetaval)
!
!-----binary128 tables for the norm pass (BINARY128 AUTO or ALWAYS)
!
      b128_nquad = 0
      b128_ntot  = 0
      if (b128_mode.gt.0) then
        call computeEcoef_q ()
        call build_gtab_q (Lmaxg,nmaxc)
      endif
!
!-----aoicule(a,b,c,d) is the spherically averaged I_abcd(R) intracule density
!     between CGTFs.  It is shared across the parallel bra-pair loop (each
!     task writes a disjoint block); the per-pair scratch is private and
!     lives inside intra_pair.
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
!-----Integration in the interval r€[0,infinity] of the radial intracule 
!
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
      integIabcd = .true.
      aoicule = 0.0d+00
!$omp parallel do schedule(dynamic) default(shared) private(ipq)
      do ipq = 1,npairs
        call intra_pairnorm (ipq,aoicule,integIabcd)
      enddo
!$omp end parallel do
      moic(:,:,:,:) = 0.0d+00
      call moicule (aoicule,moic,nprims,nmo,uout,verbose)
      f42 = 4.0d0 * pi
!     xeee = f42 * sum(rdm2 * moic(:,:,:,:))
      xeee = f42 * dot4d_kahan(rdm2,moic,nmo,nmo,nmo,nmo)
!     xeco = f42 * sum(rdm2coul * moic(:,:,:,:))
      xeco = f42 * dot4d_kahan(rdm2coul,moic,nmo,nmo,nmo,nmo)
      xexc = xeee - xeco
!     xeDF = f42 * sum(rdm2diracfock * moic(:,:,:,:)) * 0.5d+00
      xeDF = f42 * dot4d_kahan(rdm2diracfock,moic,
     &             nmo,nmo,nmo,nmo) * 0.5d+00

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
      integIabcd = .false.
      aoicule = 0.0d+00
!$omp parallel do schedule(dynamic) default(shared) private(ipq)
      do ipq = 1,npairs
        call intra_pairnorm (ipq,aoicule,integIabcd)
      enddo
!$omp end parallel do
      moic(:,:,:,:) = 0.0d+00
      call moicule (aoicule,moic,nprims,nmo,uout,verbose)
      f42 = 2.0d0 * pi
!     repee = f42 * sum(rdm2 * moic(:,:,:,:)) 
      repee = f42 * dot4d_kahan(rdm2,moic,nmo,nmo,nmo,nmo) 
!     repco = f42 * sum(rdm2coul * moic(:,:,:,:)) 
      repco = f42 * dot4d_kahan(rdm2coul,moic,nmo,nmo,nmo,nmo) 
      repxc = repee - repco
!     repDF = - f42 * sum(rdm2diracfock * moic(:,:,:,:)) * 0.5d+00 
      repDF = - f42 * dot4d_kahan(rdm2diracfock,moic,
     &                nmo,nmo,nmo,nmo) * 0.5d+00 
      ntotel = nint(totel)
      ntt = ntotel*(ntotel-1)
      ntc = ntotel*ntotel
      ntxc = ntotel
      write (uout,81) repee,repco,repxc,repDF,xeee,
     &        ntt,xeco,ntc,-xexc,ntxc,xeDF
      if (b128_mode.eq.0) then
        write (uout,82) 'NEVER (binary64 only)'
      else
        if (b128_mode.eq.1) write (uout,82) 'AUTO'
        if (b128_mode.eq.2) write (uout,82) 'ALWAYS'
        write (uout,83) b128_nquad,b128_ntot
      endif

      deallocate (aoicule)
      deallocate (moic)
      deallocate (pmap,qmap,shoff)
      deallocate (gtab)
      deallocate (yntphi,thetaval)
      deallocate (rdm2coul)
      deallocate (rdm2diracfock)
      call deallocate_space_for_inexcule  ()
      call timer (4,iradial_integ,'_radial_integ_______',-1)
      return
!
!-----Formats (host scope; referenced by the writes above)
!
 1    format (/,' #',/,
     &  ' # The routine hermitec.f must be changed increasing the',/,
     &  ' # current value of maxn = ',I2,' up to at least a value of ',
     & I4,/,' #')
 81   format (' # ',84('-')/,
     & ' # Exact (no numerical) Integrated Quantities'/,' # ',84('-')/,
     & ' # Electron-Electron repulsion energy          = ',F15.8,/,
     & ' # Coulomb repulsion energy                    = ',F15.8,/,
     & ' # Exchange-Correlation (XC) repulsion energy  = ',F15.8,/,
     & ' # Dirac-Fock XC repulsion energy              = ',F15.8,/,
     & ' # Total Ordered pairs                         = ',F15.8,
     &    ' (it should be ',I4,')',/,
     & ' # Coulomb Ordered Pairs                       = ',F15.8,
     &    ' (it should be ',I4,')',/,
     & ' # XC Ordered pairs                            = ',F15.8,
     &    ' (it should be ',I4,')',/,
     & ' # Dirac-Fock Ordered Pairs                    = ',F15.8,/,
     & ' # ',61('+'))
 82   format (' # Ordered pairs (norm pass): BINARY128 ',a)
 83   format (' # Quartets in binary128 = ',I12,' of ',I12,/,
     & ' # (the repulsion energies are always evaluated in binary64)')


      contains

      subroutine intra_pairnorm (ipq,aoicule,integIabcd)
      implicit real(kind=8) (a-h,o-z)
      integer (kind=4), intent(in) :: ipq
      logical,          intent(in) :: integIabcd
      real    (kind=8) aoicule(nprims,nprims,nprims,nprims)
      integer (kind=4) p,q,pq,r,s,rs,di,dj,dk,dl,idx,jdx,kdx,ldx
      integer (kind=4) i,j,k,ka,kb,kc,kd,ip,jp,it,itp,ival
      integer (kind=4) nx,ny,nz
      integer (kind=4) nmaxx,nmaxy,nmaxz,ntermx,ntermy,ntermz
      integer (kind=4) i1,i2,j1,j2,k1,k2,l1,l2,ia,ja,la
      real    (kind=8) mu,Iabcd,Tvec(3)
      integer (kind=4) l,Lq,lup
      real    (kind=8) almT(0:Lmaxg,-Lmaxg:Lmaxg)
      real    (kind=8) dd(0:maxn)
      real(kind=8),allocatable,dimension(:,:,:,:)     :: abcd
      real(kind=8),allocatable,dimension(:,:,:,:,:,:) :: aac
      real(kind=8),allocatable,dimension(:,:,:)       :: Jtab
      real(kind=8),allocatable,dimension(:,:,:,:,:,:) :: aabs
      real(kind=8),allocatable,dimension(:,:,:)       :: Jabs
      integer(kind=4),allocatable,dimension(:,:)      :: angp,angq
      integer(kind=4),allocatable,dimension(:,:)      :: angr,angs
      logical doquad,estim
      real(kind=8) errI,tma,skipI
      real(kind=8), parameter :: eps64 = 2.220446049250313d-16
      integer(kind=4),allocatable,dimension(:,:,:,:,:):: maxterm

      real    (kind=8), allocatable :: rsh(:)
      real(kind=8) :: intg, log_abs
      integer :: status, method

      ! Disable GSL's default error handler
      call disable_gsl_error_handler()

      allocate (rsh(0:(Lmaxg+1)**2-1))
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
      do r = 1,nshells
        dk = nzshell(r)
        kdx = shoff(r)
        do s = 1,r
          rs = r*(r-1)/2 + s
          dl = nzshell(s)

          allocate (abcd(di,dj,dk,dl))
          allocate (aac(di,dj,dk,dl,0:4*maxlprim,1:3))
          allocate (maxterm(di,dj,dk,dl,1:3))
          allocate (Jtab(0:4*maxlprim,0:4*maxlprim,0:4*maxlprim))
!
!---------binary128: ALWAYS -> straight to quartet_norm_q; AUTO -> run
!         binary64 with a rounding-error estimate and redo if needed
!
          estim  = integIabcd .and. b128_mode.eq.1
          doquad = integIabcd .and. b128_mode.eq.2
          if (integIabcd .and. b128_mode.gt.0) then
            allocate (angp(3,di),angq(3,dj),angr(3,dk),angs(3,dl))
            do k=1,3
              angp(k,1:di) = nlm(ityp(ordp(p,1:di)),k)
              angq(k,1:dj) = nlm(ityp(ordp(q,1:dj)),k)
              angr(k,1:dk) = nlm(ityp(ordp(r,1:dk)),k)
              angs(k,1:dl) = nlm(ityp(ordp(s,1:dl)),k)
            enddo
          endif
          if (estim) then
            allocate (aabs(di,dj,dk,dl,0:4*maxlprim,1:3))
            allocate (Jabs(0:4*maxlprim,0:4*maxlprim,0:4*maxlprim))
            aabs = 0.0d+00
            Jabs = 0.0d+00
          endif

          ldx = shoff(s)
          Tvec(:) = pbig(rs,:)-pbig(pq,:) 
          prodz = psmall(pq)*psmall(rs)
          addz = psmall(pq)+psmall(rs)
          mu = prodz/addz
          sqrtmu = sqrt(mu)
          piz = sqrt(pi/addz)
!
!---------this is the quantity T of the 'nintracule.pdf' document
!
          tis = norm2(Tvec)
          i1 = idx + 1
          i2 = idx + di
          j1 = jdx + 1
          j2 = jdx + dj
          k1 = kdx + 1
          k2 = kdx + dk
          l1 = ldx + 1
          l2 = ldx + dl
          if (doquad) goto 500
 
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
                        tt = Tvec(k) 
                        call hermicart (mu,tt,it+itp,maxn,dd)
                        do ival=0,it+itp
                          aac(ka,kb,kc,kd,ival,k) = eedpi * 
     &                    dd(ival) + aac(ka,kb,kc,kd,ival,k) 
                        enddo
                        if (estim) then
                          do ival=0,it+itp
                            aabs(ka,kb,kc,kd,ival,k) = abs(eedpi *
     &                      dd(ival)) + aabs(ka,kb,kc,kd,ival,k)
                          enddo
                        endif
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
          Lq = nmaxx + nmaxy + nmaxz 
!
!---------polar angles of T for the S_lm(T^) routine (standard
!         convention).  T=0 -> only l=0 survives, direction irrelevant.
!
          if (tis .gt. 1.0d-12) then
            thetaT = acos( max(-1.0d0,min(1.0d0,Tvec(3)/tis)) )
            phiT   = atan2( Tvec(2), Tvec(1) )
          else
            thetaT = 0.0d+00
            phiT   = 0.0d+00
          endif
          call compute_real_spherical_harmonics 
     &           (Lq,thetaT,phiT,rsh,.false.)
!
!---------precompute A_lm = S_lm(T^)  (independent of nx,ny,nz)
!
          do l = 0,Lq
            do m = -l,l
              almT(l,m) = rsh(l*(l+1)+m) 
            enddo
          enddo
!
          do nx=0,nmaxx
            do ny=0,nmaxy
              do nz=0,nmaxz
                lup = nx + ny + nz
                nn = lup + 2
                if (.not.integIabcd) nn = nn-1
                tm = 0.0d+00
                if (tis <= 1.0d-12) then
                  ! Para T=0, solo l=0 contribuye
                  l = 0
                  if (gsl) then
                    call radial_integral (nn,l,mu,tis,intg,status)
                  else if (expkummer) then
                    call kummer_G_dispatch (nn,l,mu,tis,method,intg)
                    if (method == 0) then
                      stop 'radial_integ.f: Illegal '//
     &                     'parameters passed to kummer_G_dispatch'
                    endif
                  elseif (kummer) then
                    call kummer_1f1_dispatch 
     &              (nn,l,mu,tis,method,intg,log_abs)
                    if (method == 0) then
                      stop 'radial_integ.f: Illegal '//
     &                     'parameters passed to kummer_G_dispatch'
                    endif
                  else
                    stop 'radial_integ.f: Unknown integration !!!'
                  endif
                  if (intg /= 0.0d0) then
                    do m = -l, l
                      tm = tm + almT(l,m) * gtab(l,m,nx,ny,nz) * intg
                    enddo
                  end if
                  if (estim) Jabs(nx,ny,nz) = abs(tm)
                else
                  ! Caso general T>0
                  do l = 0,lup
                    if (gsl) then
                      call radial_integral (nn,l,mu,tis,intg,status)
                    else if (expkummer) then
                      call kummer_G_dispatch (nn,l,mu,tis,method,intg)
                      if (method == 0) then
                        stop 'radial_integ.f: Illegal '//
     &                       'parameters passed to kummer_G_dispatch'
                      endif
                    elseif (kummer) then
                      call kummer_1f1_dispatch 
     &                     (nn,l,mu,tis,method,intg,log_abs)
                      if (method == 0) then
                        stop 'radial_integ.f: Illegal '//
     &                       'parameters passed to kummer_1f1_dispatch'
                      endif
                    else
                      stop 'radial_integ.f: Unknown integration !!!'
                    endif              
                    if (intg == 0.0d0) cycle
                    do m = -l, l
                     tm = tm + almT(l,m) * gtab(l,m,nx,ny,nz) * intg
                    enddo
                    if (estim) then
                      tma = 0.0d+00
                      do m = -l, l
                        tma = tma + abs(almT(l,m)*gtab(l,m,nx,ny,nz))
                      enddo
                      Jabs(nx,ny,nz) = Jabs(nx,ny,nz) + tma*abs(intg)
                    endif
                  enddo
                end if
                Jtab(nx,ny,nz) = tm
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
                        axyzis = axis * ayis * azis
                        Iabcd = Iabcd + axyzis * Jtab(nx,ny,nz)
                      enddo ! nz
                    enddo ! ny
                  enddo  ! nx
                  abcd(ka,kb,kc,kd) = Iabcd
                  if (estim .and. .not.doquad) then
!
!-------------------rounding bound of the binary64 chain, plus the terms
!                   dropped by the 1e-8 cutoffs of the contraction above
!
                    errI = 0.0d+00
                    skipI = 0.0d+00
                    do nx=0,ntermx
                      axis = aac(ka,kb,kc,kd,nx,1)
                      do ny=0,ntermy
                        ayis = aac(ka,kb,kc,kd,ny,2)
                        do nz=0,ntermz
                          azis = aac(ka,kb,kc,kd,nz,3)
                          errI = errI + aabs(ka,kb,kc,kd,nx,1) *
     &                      aabs(ka,kb,kc,kd,ny,2) *
     &                      aabs(ka,kb,kc,kd,nz,3) * Jabs(nx,ny,nz)
                          if (abs(axis).lt.1.0d-08 .or.
     &                        abs(ayis).lt.1.0d-08 .or.
     &                        abs(azis).lt.1.0d-08) skipI = skipI +
     &                      abs(axis*ayis*azis*Jtab(nx,ny,nz))
                        enddo
                      enddo
                    enddo
                    errI = eps64 * errI + skipI
                    if (errI.gt.b128_tol*abs(Iabcd) .and.
     &                  errI.gt.1.0d-20) doquad = .true.
                  endif
                enddo ! kd
              enddo ! kc
            enddo ! kb
          enddo ! ka
!
!---------binary128 recomputation of the whole quartet
!
 500      continue
          if (doquad) then
            call quartet_norm_q (pq,rs,di,dj,dk,dl,angp,angq,angr,angs,
     &                           maxlprim,abcd)
!$omp atomic
            b128_nquad = b128_nquad + 1
          endif
          if (integIabcd .and. b128_mode.gt.0) then
!$omp atomic
            b128_ntot = b128_ntot + 1
            deallocate (angp,angq,angr,angs)
          endif
          if (estim) deallocate (aabs,Jabs)
          aoicule(i1:i2,j1:j2,k1:k2,l1:l2) =
     &          abcd(1:di,1:dj,1:dk,1:dl)
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
          deallocate (abcd)
          deallocate (aac)
          deallocate (maxterm)
          deallocate (Jtab)
        end do !s
      end do !r
      return
      end subroutine intra_pairnorm

      end
