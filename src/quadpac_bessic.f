!
!-----------------------------------------------------------------------
!
      subroutine quadpac_bessic (wfnfile,rval,xtot,xcoul,xxc,xdf)
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
      common /quadarea/ lwfun,verbose, usemuller, first_r_value, wfnaux
      logical verbose
      logical usemuller
      logical first_r_value
      character (len=132) wfnaux
!
!-----Input variables
!
      character (len=*), intent(in)  :: wfnfile
      real (kind=8),     intent(in)  :: rval
      real (kind=8),     intent(out) :: xtot, xcoul, xxc, xdf
!
!-----Shared data.  The 2-RDMs and the angular-integral tables are built
!     once in the host below and only read inside the parallel quartet loop;
!     aoicule is shared but written in disjoint blocks (one (p,q) bra-block
!     per task); all per-pair scratch is declared inside intra_pair.
!
      real(kind=8), allocatable,dimension(:,:,:,:), save :: rdm2coul
      real(kind=8), allocatable,dimension(:,:,:,:), save :: rdm2df
      real(kind=8), allocatable,dimension (:,:), save    :: thetaval

      real(kind=8), allocatable,dimension (:,:,:,:)     :: aoicule
      real(kind=8), allocatable,dimension (:,:,:,:)     :: moic
!
!-----Shell-offset table and (p,q) bra-pair map for the parallel quartet loop
!
      integer(kind=4), allocatable,dimension (:),save :: shoff,pmap,qmap
!
!-----Angular integrals over phi
!
      real(kind=8), allocatable, save,dimension (:,:) :: yntphi

      integer (kind=4)  p,q,r,s
      integer (kind=4), save :: ncapas
      integer (kind=4), save :: npairs
!
!-----array to store the coefficiente of the Hermite polynomials up to n=20
!
      integer(kind=4), parameter :: maxn = 20
      real   (kind=8), save :: Chermite(0:maxn,0:maxn)
!
!-----Geometry-independent angular integrals  gtab(l,m,nx,ny,nz) =
!     INT S_lm f dtheta dphi  (slmangint.f90).  Built ONCE below and shared
!     read-only across every R and every quartet.  nmaxc = 4*maxlprim is the
!     safe upper bound on nx,ny,nz (matches the aac index range); the l/m
!     ranges then run up to Lmaxg = 3*nmaxc = nx+ny+nz|max.
!
      real(kind=8), allocatable,save,dimension (:,:,:,:,:) :: gtab
      integer(kind=4) nmaxc,ntab
      integer(kind=4), save :: Lmaxg
 
      logical success

      call timer (2,iquadpac_bessic,'_quadpac_bessic_____',-1)
      if (maxn.lt.4*maxlprim) then
        write (uerr,1) maxn,4*maxlprim
        stop
      else
        ixg = ngtox()
        Chermite(0:maxn,0:maxn) = 0.0d+00
        call hermitec (Chermite,maxn)
      endif

      if (first_r_value) then
        call compute_rdm (ndets,nelact,ncore,nmo,uout,uerr,nrecwfn, 
     &      cciqa,mp2iqa,icorr,verbose,wfnfile)
!
!-------Compute the Coulombic part of the RDM2 and the Dirac-Fock 
!       part of the XC RDM2
!
        if (.not.allocated(rdm2coul)) then
          allocate (rdm2coul(nmo,nmo,nmo,nmo))
        endif
        if (.not.allocated(rdm2df)) then
          allocate (rdm2df(nmo,nmo,nmo,nmo))
        endif
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
!-------Classification of GTFs in shells
!
        ncapas = sum(ngroup(1:ncent))
        call allocate_space_for_othershells (ncapas,ngtoG)
!       write (uout,'(a)')' # Classification of CGFs into shells'
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
!-------Angular index bounds.  nx,ny,nz <= nmaxc = 4*maxlprim (aac range);
!       l,m <= Lmaxg = 3*nmaxc; the phi/theta power tables need indices up to
!       ntab = 2*Lmaxg (aphi=(mm-k)+nx, pth=mm+nx+ny, etc. -- see slmangint).
!
        nmaxc = 4*maxlprim
        Lmaxg = 3*nmaxc
        ntab  = 2*Lmaxg
!
!-------INT_0^2pi sin^a cos^b dphi  and  INT_0^pi sin^{p+1} cos^q dtheta,
!       now only up to ntab (was kmaxcap-sized in the Taylor path).
!
        if (.not.allocated(yntphi)) then
          allocate (yntphi(0:ntab,0:ntab))
        endif
        call phiintegrals (ntab,ntab,yntphi)
        if (.not.allocated(thetaval)) then
          allocate (thetaval(0:ntab,0:ntab))
        endif
        call thetaintegrals (ntab,ntab,thetaval)
      endif 
!
!-----------------------------------------------------------------------
      if (first_r_value) then
!
!-------Compute the E_t^ij, E_u^kl and E_v^mn coefficients
!
        call computeEcoef ()
!
!-------Build the geometry-independent angular integrals G_lm(nx,ny,nz) ONCE.
!       This is the object that replaces the per-quartet Taylor summation.
!
        if (.not.allocated(gtab)) then
          allocate (gtab(0:Lmaxg,-Lmaxg:Lmaxg,0:nmaxc,0:nmaxc,0:nmaxc))
        endif
        call build_gtab (gtab,Lmaxg,nmaxc,nmaxc,nmaxc,ntab,ntab,
     &    yntphi,thetaval)
        npairs = ncapas*(ncapas+1)/2
        if (.not.allocated(shoff)) then
          allocate (shoff(ncapas))
        endif
        if (.not.allocated(pmap)) then
          allocate (pmap(npairs))
        endif
        if (.not.allocated(qmap)) then
          allocate (qmap(npairs))
        endif
        shoff(1) = 0
        do p = 2,ncapas
          shoff(p) = shoff(p-1) + nzshell(p-1)
        enddo
        ipq = 0
        do p = 1,ncapas
          do q = 1,p
            ipq = ipq + 1
            pmap(ipq) = p
            qmap(ipq) = q
          enddo
        enddo
      endif
!
!-----aoicule(a,b,c,d) is the spherically averaged I_abcd(R) intracule density
!     between CGTFs.  It is shared across the parallel bra-pair loop (each
!     task writes a disjoint block); the per-pair scratch is private and
!     lives inside intra_pair.
!
      allocate (aoicule(nprims,nprims,nprims,nprims))
      allocate (moic(nmo,nmo,nmo,nmo))
      rr = rval
      r2 = rr*rr
      aoicule = 0.0d+00
!$omp parallel do schedule(dynamic) default(shared) private(ipq)
      do ipq = 1,npairs
        call intra_pair (ipq,rr,r2,aoicule)
      enddo
!$omp end parallel do
      call moicule (aoicule,moic,nprims,nmo,uout,verbose)
      xtot     = dot4d_kahan(rdm2,moic,nmo,nmo,nmo,nmo)
      xcoul = dot4d_kahan(rdm2coul,moic,nmo,nmo,nmo,nmo)
      xxc     = xtot - xcoul
      xdf   = - dot4d_kahan(rdm2df,moic,nmo,nmo,nmo,nmo)*0.5d+00

      deallocate (aoicule)
      deallocate (moic)
      call timer (4,iquadpac_bessic,'_quadpac_bessic_____',-1)
      return

 1    format (/,' #',/,
     &  ' # The routine hermitec.f must be changed increasing the',/,
     &  ' # current value of maxn = ',I2,' up to at least a value of ',
     & I4,/,' #')

      contains
!
!-----------------------------------------------------------------------
!     Per bra-pair worker for one R.  For bra-pair ipq -> (p,q) it loops
!     over all ket pairs (r,s), evaluates the quartet intracule integrals
!     and scatters them (with the (pq|rs) within-pair permutational
!     symmetry) into the shared aoicule.  Every array declared here is
!     local, hence thread-private; aoicule is shared but each ipq writes a
!     disjoint (p,q) block.  The read-only tables (etij, prefac, pbig,
!     psmall, Chermite, gtab, Lmaxg, shell data) come in by host association.
!-----------------------------------------------------------------------
!
        subroutine intra_pair (ipq,rr,r2,aoicule)
        implicit real(kind=8) (a-h,o-z)
        integer (kind=4), intent(in) :: ipq
        real    (kind=8), intent(in) :: rr,r2
        real    (kind=8) aoicule(nprims,nprims,nprims,nprims)
        integer (kind=4) p,q,pq,r,s,rs,di,dj,dk,dl,idx,jdx,kdx,ldx
        integer (kind=4) i,j,k,ka,kb,kc,kd,ip,jp,it,itp,ival
        integer (kind=4) nx,ny,nz
        integer (kind=4) nmaxx,nmaxy,nmaxz,ntermx,ntermy,ntermz
        integer (kind=4) i1,i2,j1,j2,k1,k2,l1,l2,ia,ja,la
        real    (kind=8) mu,Iabcd,Tvec(3)
        integer (kind=4) l,Lq,lup
        real    (kind=8) blvec(0:Lmaxg)
        real    (kind=8) logblvec(0:Lmaxg)
        real    (kind=8) almT(0:Lmaxg,-Lmaxg:Lmaxg)
        real    (kind=8) dd(0:maxn)
        real(kind=8),allocatable,dimension(:,:,:,:)     :: abcd
        real(kind=8),allocatable,dimension(:,:,:,:,:,:) :: aac
        real(kind=8),allocatable,dimension(:,:,:)       :: Jtab
        integer(kind=4),allocatable,dimension(:,:,:,:,:):: maxterm
  
        real    (kind=8), allocatable :: rsh(:)
  
        allocate (rsh(0:(Lmaxg+1)**2-1))
!
!-------recover the bra-pair indices and AO offsets for this task
!
        p   = pmap(ipq)
        q   = qmap(ipq)
        pq  = ipq
        di  = nzshell(p)
        dj  = nzshell(q)
        idx = shoff(p)
        jdx = shoff(q)
        pf1 = prefac(pq,1)*prefac(pq,2)*prefac(pq,3)
        do r = 1,ncapas
          dk = nzshell(r)
          kdx = shoff(r)
          do s = 1,r
            rs = r*(r-1)/2 + s
            dl = nzshell(s)
  
            allocate (abcd(di,dj,dk,dl))
            allocate (aac(di,dj,dk,dl,0:4*maxlprim,1:3))
            allocate (maxterm(di,dj,dk,dl,1:3))
            allocate (Jtab(0:4*maxlprim,0:4*maxlprim,0:4*maxlprim))
  
            ldx = shoff(s)
            Tvec(:) = pbig(rs,:)-pbig(pq,:) 
            prodz = psmall(pq)*psmall(rs)
            addz = psmall(pq)+psmall(rs)
            mu = prodz/addz
            sqrtmu = sqrt(mu)
            piz = sqrt(pi/addz)
!
!-----------this is the quantity T of the 'nintracule.pdf' document
!
            tis = norm2(Tvec)
            T2 = tis * tis
            pf2 = prefac(rs,1)*prefac(rs,2)*prefac(rs,3)
            xk = pf1 * pf2 * piz**3
            alpha = 2.0d+00 * mu * rr * tis
            xmur2t2 = mu*(r2+T2)
            i1 = idx + 1
            i2 = idx + di
            j1 = jdx + 1
            j2 = jdx + dj
            k1 = kdx + 1
            k2 = kdx + dk
            l1 = ldx + 1
            l2 = ldx + dl
!
!-----------Screen: the whole quartet is damped by exp[-mu(|R|-|T|)^2]
!           (the stable Gaussian factor of B_l), so it underflows to
!           zero when mu(|R|-|T|)^2 is large.  (This replaces the old
!           mu(R^2+T^2) screen, which was too pessimistic for the
!           numerically-stable Bessel form.)
!
            if (mu*(rr-tis)**2 .gt. 700.0D+00) then
              abcd(1:di,1:dj,1:dk,1:dl) = 0.0d+00
              goto 222
            endif
 
!-----------aac = The A_nx^{abcd},  A_ny^{abcd},  A_nz^{abcd} 
!           quantities of the 'nuevointracule.pdf' document
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
     &                      dd(ival) + aac(ka,kb,kc,kd,ival,k) 
                          enddo
                        enddo
                      enddo
                    enddo
                  enddo
                enddo
              enddo
            enddo
!
!-----------Maximum Cartesian-Hermite orders actually needed by any of
!           the di*dj*dk*dl Cartesian components of this quartet
!
            nmaxx = maxval(maxterm(1:di,1:dj,1:dk,1:dl,1))
            nmaxy = maxval(maxterm(1:di,1:dj,1:dk,1:dl,2))
            nmaxz = maxval(maxterm(1:di,1:dj,1:dk,1:dl,3))
!
!-----------Bessel / spherical-harmonic coupling table.
!
!           J(nx,ny,nz) = SUM_{l=0}^{nx+ny+nz} B_l
!                         SUM_{m=-l}^{l} S_lm(T^) G_lm(nx,ny,nz)
!
!           with  B_l = exp[-mu(|R|-|T|)^2] i_l,scaled(2 mu |R| |T|)
!           and   G_lm precomputed in gtab (host association).  The
!           l-sum is finite: G_lm = 0 for l > nx+ny+nz.  Nothing here
!           depends on the Cartesian indices, so J is tabulated once in
!           Jtab(nx,ny,nz) and contracted with aac below (unchanged).
!
            Lq = nmaxx + nmaxy + nmaxz 
!
!-----------polar angles of T for the S_lm(T^) routine (standard
!           convention).  T=0 -> only l=0 survives, direction irrelevant.
!
            if (tis .gt. 1.0d-12) then
              thetaT = acos( max(-1.0d0,min(1.0d0,Tvec(3)/tis)) )
              phiT   = atan2( Tvec(2), Tvec(1) )
            else
              thetaT = 0.0d+00
              phiT   = 0.0d+00
            endif
            xarg = 2.0d+00 * mu * rr * tis

!           !>>>> USER ROUTINE (1): scaled modified spherical Bessel
!           times the Gaussian,  blvec(l) = exp[-mu(|R|-|T|)^2] i_l^s(x).
!           Adapt this CALL to the parent package's actual routine.
!
            call safe_compute_B_l (Lq,mu,rr,tis,blvec,logblvec,success)
            call compute_real_spherical_harmonics 
     &        (Lq,thetaT,phiT,rsh,.false.)
!
!-----------precompute A_lm = B_l S_lm(T^)  (independent of nx,ny,nz)
!
            do l = 0,Lq
              do m = -l,l
                almT(l,m) = blvec(l) * rsh(l*(l+1)+m)
              enddo
            enddo
!
            do nx=0,nmaxx
              do ny=0,nmaxy
                do nz=0,nmaxz
                  lup = nx + ny + nz
                  tm = 0.0d+00
                  do l = 0,lup
                    do m = -l,l
                      tm = tm + almT(l,m)*gtab(l,m,nx,ny,nz)
                    enddo
                  enddo
                  Jtab(nx,ny,nz) = rr**(nx+ny+nz+2) * tm
                  otm = Jtab(nx,ny,nz)
                enddo ! nz
              enddo ! ny
            enddo ! nx
!
!-----------Evaluate the I_abcd(R) integral: contract the A_nx/A_ny/A_nz
!           coefficients with the precomputed coupling table Jtab
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
                          Iabcd = Iabcd
     &                      + axis*ayis*azis * Jtab(nx,ny,nz)
                        enddo ! nz
                      enddo ! ny
                    enddo  ! nx
                    abcd(ka,kb,kc,kd) = Iabcd
                  enddo ! kd
                enddo ! kc
              enddo ! kb
            enddo ! ka
 222        continue
            aoicule(i1:i2,j1:j2,k1:k2,l1:l2) =
     &            abcd(1:di,1:dj,1:dk,1:dl)
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
        end subroutine intra_pair
      end subroutine quadpac_bessic
