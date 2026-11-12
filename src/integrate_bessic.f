!
!-----------------------------------------------------------------------
!
!     BESSEL / SPHERICAL-HARMONIC variant of newtayloric.f
!     --------------------------------------------------------------------
!     The angular coupling table J(nx,ny,nz) is built here from the
!     addition-theorem expansion of  exp[2 mu R.T]  (see next.txt,
!     lebedev.pdf), NOT from the >5000-term Taylor expansion:
!
!        J(nx,ny,nz) = SUM_{l=0}^{nx+ny+nz} B_l SUM_{m=-l}^{l} S_lm(T^) G_lm
!        B_l   = exp[-mu(|R|-|T|)^2] * i_l,scaled(2 mu |R| |T|)   (routine 1)
!        S_lm(T^)                                                 (routine 2)
!        G_lm  = INT S_lm f dtheta dphi   <- slmangint.f90 (built once, gtab)
!
!     The l-sum is FINITE (G_lm = 0 for l > nx+ny+nz) and everything is kept
!     O(1) in magnitude by the scaled Bessel function, so kmaxcap/thresh and
!     the whole log-space damping machinery of newtayloric.f are gone.
!     kmaxcap/thresh remain in the argument list only for a drop-in swap of
!     the CALL name; they are UNUSED here.
!
!     Two external routines (marked  !>>>> USER ROUTINE  below) are assumed
!     to exist in the parent package, exactly as stated in next.txt:
!        (1) the scaled modified spherical Bessel times the Gaussian, B_l
!        (2) the real spherical harmonics S_lm at the polar angles of T
!     Their names/interfaces here are PLACEHOLDERS -- adapt the two CALL
!     lines to the parent's actual routines.  Their normalisation/phase MUST
!     match slmangint.f90 (orthonormal real SH, no Condon-Shortley phase,
!     m>0 -> cos, m<0 -> sin); see the slmangint header.
!-----------------------------------------------------------------------
!
      subroutine integrate_bessic (verbose,quts,wfnfile,rval,xicMO,
     &  xicMOCoul,xicXC,xicMODF,nrval,iquts,thresh,usemuller,
     &  gsl,kummer,expkummer)
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
      logical,           intent(in)  :: verbose,quts
      character (len=*), intent(in)  :: wfnfile
      integer (kind=4),  intent(in)  :: nrval,iquts
      real (kind=8),     intent(in)  :: rval(nrval)
      real (kind=8),     intent(out) :: xicMO (nrval)
      real (kind=8),     intent(out) :: xicMOCoul(nrval)
      real (kind=8),     intent(out) :: xicXC(nrval)
      real (kind=8),     intent(out) :: xicMODF(nrval)
      real (kind=8),     intent(inout)  :: thresh

!
!-----Shared data.  The 2-RDMs and the angular-integral tables are built
!     once in the host below and only read inside the parallel quartet loop;
!     aoicule is shared but written in disjoint blocks (one (p,q) bra-block
!     per task); all per-pair scratch is declared inside intra_pair.
!
      real(kind=8), allocatable,dimension (:,:,:,:)     :: rdm2coul
      real(kind=8), allocatable,dimension (:,:,:,:)     :: rdm2diracfock
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
!
!-----Geometry-independent angular integrals  gtab(l,m,nx,ny,nz) =
!     INT S_lm f dtheta dphi  (slmangint.f90).  Built ONCE below and shared
!     read-only across every R and every quartet.  nmaxc = 4*maxlprim is the
!     safe upper bound on nx,ny,nz (matches the aac index range); the l/m
!     ranges then run up to Lmaxg = 3*nmaxc = nx+ny+nz|max.
!
      real(kind=8), allocatable,dimension (:,:,:,:,:) :: gtab
      integer(kind=4) nmaxc,Lmaxg,ntab
 
      real(kind=8), allocatable,dimension (:,:)  :: sprim,smomo,coeft

      logical success
      logical integIabcd
      logical usemuller
      logical gsl,kummer,expkummer

      write (uout,'(a)')' # ++++++++++++++++++++++++++++++++++++++++++'
      write (uout,'(a)')' # Bessel / spherical-harmonic coupling table'
      write (uout,'(a,1PE15.8)')' # (kmaxcap/thresh unused) thresh = ',
     &  thresh
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
      if (usemuller) then
        write (uout,50) 'Natural Orbitals'
        write (uout,'(a)') ' # Müller approximation for XC DM is used'
      else
        write (uout,50) 'Canonical Orbitals'
      endif
      call timer (2,irdm2,'_rdm2MOs____________',-1)
      call compute_rdm (ndets,nelact,ncore,nmo,uout,uerr,nrecwfn, 
     &    cciqa,mp2iqa,icorr,verbose,wfnfile)
      call timer (4,irdm2,'_rdm2MOs____________',-1)
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
!-----check orthonormality of natural MOs
!
      allocate (sprim(nprims,nprims))
      allocate (smomo(nmo,nmo))
      allocate (coeft(nmo,nprims))
      call gtogto (sprim)
      coeft = matmul(coef(1:nmo,1:nprims),sprim(:,:))
      smomo = matmul(coeft(:,:),transpose(coef(1:nmo,1:nprims)))
!
!     Si prefieres utilizar LAPACK comenta las dos lineas anteriores 
!     y sustituyelas por las siguientes 
!
!     call dgemm('N','N',nmo,nprims,nprims,1.0d0,coef(1,1),nmo, 
!    &     sprim(1,1),nprims,0.0d0,coeft(1,1),nmo)
!
!     call dgemm('N','T',nmo,nmo,nprims,1.0d0,coeft(1,1),nmo, 
!    &     coef(1,1),nmo,0.0d0,smomo(1,1),nmo)

c
c-----Orthogonality between natural MOs
c
      write (uout,*)
      write (uout,'(1x,a)') '# Testing orthonormality of natural MOs'
      do i=1,nmo
        do j=i,nmo
          write (uout,223) i,j,smomo(i,j)
        enddo
      enddo
      coeft = matmul(coef(1+nmo:nmo+nmo,1:nprims),sprim(:,:))
      smomo = matmul(coeft(:,:),transpose(coef(1+nmo:nmo+nmo,1:nprims)))
!
!     Si prefieres utilizar LAPACK comenta las dos lineas anteriores 
!     y sustituyelas por las siguientes 
!
!     call dgemm('N','N',nmo,nprims,nprims,1.0d0,coef(1+nmo,1),nmo,&
!    &     sprim(1,1),nprims,0.0d0,coeft(1,1),nmo)
!
!     call dgemm('N','T',nmo,nmo,nprims,1.0d0,coeft(1,1),nmo,
!    &     coef(1+nmo,1),nmo,0.0d0,smomo(1,1),nmo)
c
c-----Orthogonality between natural MOs
c
      write (uout,*)
      write (uout,'(1x,a)') '# Testing orthonormality of canonical MOs'
      do i=1,nmo
        do j=i,nmo
          write (uout,224) i,j,smomo(i,j)
        enddo
      enddo

      deallocate (sprim,smomo,coeft)
!
!-----Checking the fulfillment of Int rho d1 = N
!
      write (uout,*)
      write (uout,'(1x,a)') '# Testing fulfillment of INT rho d1 = N'
      rho1int = 0.0d+00
      do p=1,nmo
        rho1int = rho1int + rdm1(p,p)
      enddo
      reler = (rho1int-totel)/totel*100.0d+00
      write (uout,221) rho1int,totel,reler
!
!-----Checking the fulfillment of Int rho_2 d1 d2 = N(N-1)
!
      write (uout,*)
      write (uout,'(1x,a)') '# Idem INT rho2 d1 d2 = N(N-1)'
      rho2int = 0.0d+00
      ntotel = nint(totel)
      ntt = ntotel*(ntotel-1)
      do p=1,nmo
        do r=1,nmo
          rho2int = rho2int + rdm2(p,p,r,r)
        enddo
      enddo
      reler = (rho2int-ntt)/dble(ntt)*100.0d+00
      write (uout,220) rho2int,dble(ntt),reler
!
!-----Classification of GTFs in shells
!
      call timer (2,iintegrate_bessic,'_integrate_bessic___',-1)
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
      write (uout,'(a)')' # Computing Phi / Theta angular integrals'
      allocate (yntphi(0:ntab,0:ntab))
      call phiintegrals (ntab,ntab,yntphi)
      allocate (thetaval(0:ntab,0:ntab))
      call thetaintegrals (ntab,ntab,thetaval)
!
!-----------------------------------------------------------------------
!
!-----Compute the E_t^ij, E_u^kl and E_v^mn coefficients
!
      write (uout,51)
      call computeEcoef ()
!
!-----Build the geometry-independent angular integrals G_lm(nx,ny,nz) ONCE.
!     This is the object that replaces the per-quartet Taylor summation.
!
      write (uout,'(a,I4)')
     &  ' # Building S_lm angular integrals up to l = ',Lmaxg
      allocate (gtab(0:Lmaxg,-Lmaxg:Lmaxg,0:nmaxc,0:nmaxc,0:nmaxc))
      call build_gtab (gtab,Lmaxg,nmaxc,nmaxc,nmaxc,ntab,ntab,
     &  yntphi,thetaval)
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
!       xicMO(irval)     = sum(rdm2 * moic(:,:,:,:))
        xicMO(irval)     = dot4d_kahan(rdm2,moic,nmo,nmo,nmo,nmo)
!       xicMOCoul(irval) = sum(rdm2coul * moic(:,:,:,:))
        xicMOCoul(irval) = dot4d_kahan(rdm2coul,moic,
     &                     nmo,nmo,nmo,nmo)   
        xicXC(irval)     = xicMO(irval) - xicMOCoul(irval)
!       xicMODF(irval)   = - sum(rdm2diracfock * moic(:,:,:,:))*0.5d+00
        xicMODF(irval)   = - dot4d_kahan(rdm2diracfock,moic,
     &                       nmo,nmo,nmo,nmo)*0.5d+00
        if (quts) write (iquts,35)
!
!-------Print the densities
!
        write (uout,100) rval(irval),xicMO(irval),xicMOCoul(irval),
     &    xicXC(irval),xicMODF(irval)
      enddo
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

      deallocate (aoicule)
      deallocate (moic)
      deallocate (pmap,qmap,shoff)
      deallocate (gtab)
      deallocate (yntphi,thetaval)
      deallocate (rdm2coul)
      deallocate (rdm2diracfock)
      call deallocate_space_for_inexcule  ()
      call timer (4,iintegrate_bessic,'_integrate_bessic___',-1)
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
 35   format (' # ',54('-+-'))
 223  format (' # Overlap between NATMOs ',I3,' and ',I2,' = ',F25.18)
 224  format (' # Overlap between CANMOs ',I3,' and ',I2,' = ',F25.18)
 221  format (' # INT RHO1(1) d1      = ',F15.8,/,
     &        ' # N                   = ',F15.8,/,
     &        ' # Relative error      = ',F15.8,' %',/)
 220  format ( ' # INT RHO2(1,2) d1 d2 = ',F15.8,/,
     &         ' # N(N-1)              = ',F15.8,/,
     &         ' # Relative error      = ',F15.8,' %',/)
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
!---------Screen: the whole quartet is damped by exp[-mu(|R|-|T|)^2]
!         (the stable Gaussian factor of B_l), so it underflows to
!         zero when mu(|R|-|T|)^2 is large.  (This replaces the old
!         mu(R^2+T^2) screen, which was too pessimistic for the
!         numerically-stable Bessel form.)
!
          if (mu*(rr-tis)**2 .gt. 700.0D+00) then
            abcd(1:di,1:dj,1:dk,1:dl) = 0.0d+00
            if (quts) write (iquts,34) rr,p,q,r,s
            goto 222
          endif
 
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
!---------Bessel / spherical-harmonic coupling table.
!
!         J(nx,ny,nz) = SUM_{l=0}^{nx+ny+nz} B_l
!                       SUM_{m=-l}^{l} S_lm(T^) G_lm(nx,ny,nz)
!
!         with  B_l = exp[-mu(|R|-|T|)^2] i_l,scaled(2 mu |R| |T|)
!         and   G_lm precomputed in gtab (host association).  The
!         l-sum is finite: G_lm = 0 for l > nx+ny+nz.  Nothing here
!         depends on the Cartesian indices, so J is tabulated once in
!         Jtab(nx,ny,nz) and contracted with aac below (unchanged).
!
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
          xarg = 2.0d+00 * mu * rr * tis

!         !>>>> USER ROUTINE (1): scaled modified spherical Bessel
!         times the Gaussian,  blvec(l) = exp[-mu(|R|-|T|)^2] i_l^s(x).
!         Adapt this CALL to the parent package's actual routine.
!
          call safe_compute_B_l (Lq,mu,rr,tis,blvec,logblvec,success)
          call compute_real_spherical_harmonics 
     &           (Lq,thetaT,phiT,rsh,.false.)
!
!---------precompute A_lm = B_l S_lm(T^)  (independent of nx,ny,nz)
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
                if (quts) then
                  write (iquts,33) rr,p,q,r,s,nx,ny,nz,tm,otm
                endif
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
                        Iabcd = Iabcd
     &                    + axis*ayis*azis * Jtab(nx,ny,nz)
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
 33   format (' # R = ',F12.5,'  Quartet ',4I4,3x,'(px,py,pz)=',3I2,3x,
     & 'J_{px,py,pz} = ',1PE15.8,3x,
     &  'R^{px+py+pz+2) J_{px,py,pz} = ',1PE15.8)
 34   format (' # R = ',F12.5,'  Quartet ',4I4,'   mu(R^2+T^2) > 1000')
      end subroutine intra_pair




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
                      stop 'integrate_bessic.f: Illegal '//
     &                     'parameters passed to kummer_G_dispatch'
                    endif
                  elseif (kummer) then
                    call kummer_1f1_dispatch 
     &              (nn,l,mu,tis,method,intg,log_abs)
                    if (method == 0) then
                      stop 'integrate_bessic.f: Illegal '//
     &                     'parameters passed to kummer_G_dispatch'
                    endif
                  else
                    stop 'integrate_bessic.f: Unknown integration !!!'
                  endif
                  if (intg /= 0.0d0) then
                    do m = -l, l
                      tm = tm + almT(l,m) * gtab(l,m,nx,ny,nz) * intg
                    enddo
                  end if
                else
                  ! Caso general T>0
                  do l = 0,lup
                    if (gsl) then
                      call radial_integral (nn,l,mu,tis,intg,status)
                    else if (expkummer) then
                      call kummer_G_dispatch (nn,l,mu,tis,method,intg)
                      if (method == 0) then
                        stop 'integrate_bessic.f: Illegal '//
     &                       'parameters passed to kummer_G_dispatch'
                      endif
                    elseif (kummer) then
                      call kummer_1f1_dispatch 
     &                     (nn,l,mu,tis,method,intg,log_abs)
                      if (method == 0) then
                        stop 'integrate_bessic.f: Illegal '//
     &                       'parameters passed to kummer_1f1_dispatch'
                      endif
                    else
                      stop 'integrate_bessic.f: Unknown integration !!!'
                    endif              
                    if (intg == 0.0d0) cycle
                    do m = -l, l
                     tm = tm + almT(l,m) * gtab(l,m,nx,ny,nz) * intg
                    enddo
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
                enddo ! kd
              enddo ! kc
            enddo ! kb
          enddo ! ka
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
