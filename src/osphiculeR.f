!
!-----------------------------------------------------------------------
!
      subroutine osphiculeR (verbose,wfnfile,rval,
     &  xicMO,xicMOCoul,xicXC,xicMODF,nrval,nptsgj)
!
!------------------------------------------------------------------------
!
!-----Computes the spherically averaged intracule density of an atom
!
!     NUMERICAL NOTES (version 2, 2024):
!     --------------------------------
!     Major changes with respect to version 1:
!       1. The product EtijEtij * piz * xmupow is computed in
!          log-scale to avoid underflow/overflow.
!       2. The threshold on beta is relative (beta > 700) instead of
!          absolute (beta > 1000).
!       3. The accumulation in sumxyz uses Kahan summation.
!       4. The accumulation in xicMO, xicMOCoul, xicMODF uses Kahan
!          summation.
!       5. hermicart is called with the same interface, but the new
!          version of hermicart does not use the c() argument.
!       6. The local variable mu has been renamed to muval to avoid
!          a name clash with the intrinsic function argument list.
!
!-----------------------------------------------------------------------
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
      integer(kind=4), parameter :: nxpx = 15
      logical,           intent(in)  :: verbose
      character (len=*), intent(in)  :: wfnfile
      integer (kind=4),  intent(in)  :: nrval
      real (kind=8),     intent(in)  :: rval(nrval)
      real (kind=8),     intent(out) :: xicMO (nrval)
      real (kind=8),     intent(out) :: xicMOCoul(nrval)
      real (kind=8),     intent(out) :: xicXC(nrval)
      real (kind=8),     intent(out) :: xicMODF(nrval)

      real(kind=8), parameter :: fourpi   = 4.0d+00*pi
      real(kind=8), parameter :: beta_max = 700.0d0

      real(kind=8), allocatable,dimension (:,:,:,:,:,:) :: aac
      real(kind=8), allocatable,dimension (:,:,:,:)     :: aoicule
      real(kind=8), allocatable,dimension (:,:,:,:)     :: rdm2coul
      real(kind=8), allocatable,dimension (:,:,:,:)     :: rdm2diracfock
      real(kind=8), allocatable,dimension (:,:,:,:)     :: moic
      real(kind=8), allocatable,dimension (:,:,:,:)     :: abcd
      real(kind=8), allocatable,dimension (:)           :: dd
      real(kind=8), allocatable,dimension (:,:)         :: tmp,tmpt
      real(kind=8), allocatable,dimension (:,:)         :: intphi
      real(kind=8), allocatable,dimension (:)           :: rrpow
      integer (kind=4), allocatable, dimension (:,:,:,:,:) :: maxterm
      integer (kind=4), allocatable, dimension (:,:,:,:,:) :: maxt
      real(kind=8),allocatable,dimension (:,:)  :: thetaval
!
!-----Abscissas and weights of the Gauss-Jacobi quadrature
!
      real(kind=8), allocatable,dimension (:,:) :: xgj,wgj

      real    (kind=8)  muval
      real    (kind=8)  sqrtmu
      real    (kind=8)  piz
      real    (kind=8)  beta
      real    (kind=8)  fx
      real    (kind=8)  prodz
      real    (kind=8)  addz
      real    (kind=8)  rr
      real    (kind=8)  r2
      real    (kind=8)  tt
      real    (kind=8)  pow
      real    (kind=8)  coefe1
      real    (kind=8)  coefe2
      real    (kind=8)  EtijEtij
      real    (kind=8)  eedpi
      real    (kind=8)  axis
      real    (kind=8)  ayis
      real    (kind=8)  azis
      real    (kind=8)  annn
      real    (kind=8)  yntphi
      real    (kind=8)  ynttheta
      real    (kind=8)  ynt
      real    (kind=8)  alim
      real    (kind=8)  blim
      real    (kind=8)  dbi
      real    (kind=8)  lab
      real    (kind=8)  lcd
      real    (kind=8)  labcd
      integer (kind=4)  p,q,pq,r,s,rs,di,dj,dk,dl,idx,jdx,kdx,ldx
      integer (kind=4)  a,b,c,d
      integer (kind=4)  ia,ja,ka,la
      integer (kind=4)  i,j,k
      integer (kind=4)  ic,m,ii
      integer (kind=4)  irval
      integer (kind=4)  ixg
      integer (kind=4)  nshells
      integer (kind=4)  kindis
      integer (kind=4)  nxmaxtheta
      integer (kind=4)  nxymax
      integer (kind=4)  it,itp,ij,ijp,ip,jp
      integer (kind=4)  ival
      integer (kind=4)  nx,ny,nz,nxy,nxyz
      integer (kind=4)  ntermx,ntermy,ntermz
      integer (kind=4)  i1,i2,j1,j2,k1,k2,l1,l2
!
!-----array to store the coefficiente of the Hermite polynomials up to n=20
!
      integer(kind=4), parameter :: maxn = 20
      real   (kind=8) Chermite(0:maxn,0:maxn)
!
!-----Kahan summation variables
!
      real(kind=8) :: sumxyz, csum, yk, tk
      real(kind=8) :: xicMO_k, xicMO_c
      real(kind=8) :: xicMOCoul_k, xicMOCoul_c
      real(kind=8) :: xicMODF_k, xicMODF_c
!
!-----Log-scale variables
!
      real(kind=8) :: log_eedpi, log_EtijEtij, log_piz, log_xmupow
      real(kind=8) :: sign_eedpi
!
!-----This routine can only be used in atoms
!
      if (ncent.ne.1) write (uerr,200)
      if (ncent.ne.1) stop
!
      if (maxn.lt.4*maxlprim) then
        write (uerr,1) maxn,4*maxlprim
        stop
      else
        allocate (dd(0:maxn))
        ixg = ngtox()
        allocate (abcd(ixg,ixg,ixg,ixg))
        allocate (tmp(ixg,ixg))
        allocate (tmpt(ixg,ixg))
        allocate (aac(ixg,ixg,ixg,ixg,0:4*maxlprim,1:3))
        allocate (maxterm(ixg,ixg,ixg,ixg,1:3))
        allocate (maxt(nprims,nprims,nprims,nprims,1:3))
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
      write (uout,*) '# Done'
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
      call timer (2,isphiculeR,'_osphiculeR_________',-1)
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
!-----Allocate space for the abscissas and weights of the Gauss-Jacobi
!     quadratures that will be used to compute the theta integrals,
!     and compute them for nxy from 0 to 2*maxlprim
!
      allocate (xgj(0:2*maxlprim,nptsgj))
      allocate (wgj(0:2*maxlprim,nptsgj))
      alim = -1.0d0
      blim = +1.0d0
      kindis = 4
      do i=0,2*maxlprim
        dbi = dble(i)
        call cdgqf (nptsgj, kindis, dbi, dbi, xgj(i,:), wgj(i,:))
      enddo
!
!-----We compute the INT_0^2\pi sin(x)^(2k) cos(x)^(2l) dx integrals
!
      allocate (intphi(0:2*maxlprim,0:2*maxlprim))
      call init_yntphi_recursive (intphi,2*maxlprim)
      nxmaxtheta = 4*maxlprim
      nxymax = 2 * nxmaxtheta
      allocate(thetaval(0:nxymax, 0:nxmaxtheta+2))
      call init_thetaval_recursive (thetaval,nxymax,nxmaxtheta+2)

!
!-----Print the r values of the grid
!
      write (uout,20)
      write (uout,21) (rval(i),i=1,nrval)
      write (uout,20)
!
!-----------------------------------------------------------------------
!
!-----Compute the E_t^ij, E_u^kl and E_v^mn coefficients
!
      write (uout,51)
      call computeEcoef ()
      write (uout,*) '# Done'
!
!-----aoicule(a,b,c,d) is the spherically averaged I_abcd(R) intracule density
!     between CGTFs
!
      allocate (aoicule(nprims,nprims,nprims,nprims))
      allocate (moic(nmo,nmo,nmo,nmo))
      allocate (rrpow(0:4*maxlprim))
!
      write (uout,9)
!
!-----Loop over the different values of R
!
      do irval=1,nrval
!
!-------The distance at which the intracule density is computed
!
        rr = rval(irval)
        rrpow(0) = 1.0d+00
        do j=1,4*maxlprim
          rrpow(j) = rrpow(j-1) * rr
        enddo
        r2 = rr * rr
        aoicule = 0.0d+00
!
!-------Run over p,q,r,s shells of primitive Cartesian Gaussians (CG)
!
        pq = 0
        idx = 0
        do p = 1,nshells
          di = nzshell(p)
          jdx = 0
          do q = 1,p
            dj = nzshell(q)
            pq = pq + 1
            lab  = lang(p) + lang(q)
            kdx = 0
            rs = 0
            do r = 1,nshells
              dk = nzshell(r)
              ldx = 0
              do s = 1,r
                rs = rs + 1
                dl = nzshell(s)
                lcd  = lang(r) + lang(s)
                prodz = psmall(pq)*psmall(rs)
                addz  = psmall(pq)+psmall(rs)
                muval = prodz/addz
                sqrtmu = sqrt(muval)
                piz = sqrt(pi/addz)
                beta = muval * r2
                i1 = idx + 1
                i2 = idx + di
                j1 = jdx + 1
                j2 = jdx + dj
                k1 = kdx + 1
                k2 = kdx + dk
                l1 = ldx + 1
                l2 = ldx + dl
!
!---------------Relative threshold on beta (instead of absolute 1000)
!
                if (beta .gt. beta_max) then
                  abcd(1:di,1:dj,1:dk,1:dl) = 0.0d+00
                  goto 222
                endif
                fx = exp(-beta)
!
!---------------Sum of angular momenta of the four shells of CGs
!
                labcd = lab + lcd
!
!---------------aac = The A_nx^{abcd},  A_ny^{abcd},  A_nz^{abcd} quantities
!
                do k=1,3 ! the three components x,y,z
                  do ka=1,di
                    a = ordp(p,ka)
                    itipa = ityp(a)
                    i = nlm(itipa,k)
                    do kb=1,dj
                      b = ordp(q,kb)
                      itipb = ityp(b)
                      j = nlm(itipb,k)
                      ij = i+j
                      do kc=1,dk
                        c = ordp(r,kc)
                        itipc = ityp(c)
                        ip = nlm(itipc,k)
                        do kd=1,dl
                          d = ordp(s,kd)
                          itipd = ityp(d)
                          jp = nlm(itipd,k)
                          ijp = ip+jp
                          aac(ka,kb,kc,kd,:,k) = 0.0d+00
                          pow = 1.0d+00
                          maxterm(ka,kb,kc,kd,k) = ij+ijp
                          do it=0,ij
                            if (mod(it+ij,2).ne.0) then
                              pow = -pow
                              cycle
                            endif
                            coefe1 = etij(pq,it,i,j,k)*prefac(pq,k)
                            do itp=0,ijp
                              if (mod(itp+ijp,2).ne.0) cycle
                              coefe2 = etij(rs,itp,ip,jp,k)*prefac(rs,k)
                              EtijEtij = coefe1 * coefe2
!
!---------------------------Log-scale computation of eedpi
!
                              if (EtijEtij .ne. 0.0d0) then
                                log_EtijEtij = log(abs(EtijEtij))
                              else
                                log_EtijEtij = -huge(1.0d0)
                              endif
                              log_piz    = log(piz)
                              log_xmupow = 0.5d0 * dble(it+itp) *
     &                                             log(muval)
                              log_eedpi = log_EtijEtij + log_piz + 
     &                                    log_xmupow
!
!-----------------------------If log_eedpi is too negative, the term
!                             underflows; skip it.
!
                              if (log_eedpi .lt. -700.0d0) then
                                cycle
                              endif
!
                              sign_eedpi = sign(1.0d0, EtijEtij * pow)
                              eedpi = sign_eedpi * exp(log_eedpi)
!
                              tt = 0.0d+00
                              ittp = it+itp
                              call hermicart (muval,tt,ittp,maxn,dd)
                              do ival=0,it+itp
                                aac(ka,kb,kc,kd,ival,k) = eedpi *
     &                          dd(ival) + aac(ka,kb,kc,kd,ival,k)
                              enddo
                            enddo
                            pow = -pow
                          enddo
                        enddo
                      enddo
                    enddo
                  enddo
                enddo
!
!---------------Evaluate the I_abcd(R) integral
!
                do ka=1,di
                  do kb=1,dj
                    do kc=1,dk
                      do kd=1,dl
                        ntermx = maxterm(ka,kb,kc,kd,1)
                        ntermy = maxterm(ka,kb,kc,kd,2)
                        ntermz = maxterm(ka,kb,kc,kd,3)
!
!-----------------------Only even values of nx and ny are neccessary.
!                       Otherwise the PHI integral (yntphi) is zero.
!
                        sumxyz = 0.0d0
                        csum   = 0.0d0
                        do nx=0,ntermx,2
                          axis = aac(ka,kb,kc,kd,nx,1)
                          if (abs(axis).lt.1.0d-08) cycle
                          do ny=0,ntermy,2
                            ayis = aac(ka,kb,kc,kd,ny,2)
                            if (abs(ayis).lt.1.0d-08) cycle
                            ayis = axis * ayis
                            nxy = nx+ny
!
!---------------------------Evaluate the angular integral of
!                           exp[-2\mu*vec{R}*vec{T}] * X^nx * Y^ny * Z^nz
!
                            yntphi = intphi(ny/2,nx/2)
                            do nz=0,ntermz
                              nxyz = nx+ny+nz
                              if (mod(nxyz,2) .ne. 0) cycle
!
!-----------------------------Final value of the A_{nx,ny,nz) coefficient
!
                              azis = aac(ka,kb,kc,kd,nz,3)
                              if (abs(azis).lt.1.0d-08) cycle
                              annn = ayis * azis
                              ynttheta = thetaval(nx+ny,nz) * fx
                              ynt = ynttheta * yntphi * rrpow(nxyz)
!
!-----------------------------Kahan summation in sumxyz
!
                              yk = annn * ynt - csum
                              tk = sumxyz + yk
                              csum = (tk - sumxyz) - yk
                              sumxyz = tk
                            enddo ! nz
                          enddo ! ny
                        enddo  ! nx
                        abcd(ka,kb,kc,kd) = sumxyz * r2
                      enddo ! kd
                    enddo ! kc
                  enddo ! kb
                enddo ! ka
 222            continue
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
                ldx = ldx + dl
              end do !s
              kdx = kdx + dk
            end do !r
            jdx = jdx + dj
          end do !q
          idx = idx + di
        end do !p

        call moicule (aoicule,moic,nprims,nmo,uout,verbose)
!
!-------Kahan summation in the final sums
!
        xicMO_k = 0.0d0
        xicMO_c = 0.0d0
        do p=1,nmo
          do q=1,nmo
            do r=1,nmo
              do s=1,nmo
                yk = rdm2(p,q,r,s) * moic(p,q,r,s) - xicMO_c
                tk = xicMO_k + yk
                xicMO_c = (tk - xicMO_k) - yk
                xicMO_k = tk
              enddo
            enddo
          enddo
        enddo
        xicMO(irval) = xicMO_k / fourpi
!
        xicMOCoul_k = 0.0d0
        xicMOCoul_c = 0.0d0
        do p=1,nmo
          do q=1,nmo
            do r=1,nmo
              do s=1,nmo
                yk = rdm2coul(p,q,r,s) * moic(p,q,r,s) - xicMOCoul_c
                tk = xicMOCoul_k + yk
                xicMOCoul_c = (tk - xicMOCoul_k) - yk
                xicMOCoul_k = tk
              enddo
            enddo
          enddo
        enddo
        xicMOCoul(irval) = xicMOCoul_k / fourpi
!
        xicXC(irval) = xicMO(irval) - xicMOCoul(irval)
!
        xicMODF_k = 0.0d0
        xicMODF_c = 0.0d0
        do p=1,nmo
          do q=1,nmo
            do r=1,nmo
              do s=1,nmo
                yk = rdm2diracfock(p,q,r,s) * moic(p,q,r,s) - xicMODF_c
                tk = xicMODF_k + yk
                xicMODF_c = (tk - xicMODF_k) - yk
                xicMODF_k = tk
              enddo
            enddo
          enddo
        enddo
        xicMODF(irval) = - xicMODF_k * 0.5d+00 / fourpi
!
        write (uout,100) rval(irval),xicMO(irval),xicMOCoul(irval),
     &    xicXC(irval),xicMODF(irval)

      enddo ! R values
!
!-----Output
!
 100  format (1x,F15.8,8(2x,E22.15))

      deallocate (aoicule)
      deallocate (moic)
      deallocate (dd)
      deallocate (abcd)
      deallocate (aac)
      deallocate (rdm2coul)
      deallocate (rdm2diracfock)

      call timer (4,isphiculeR,'_osphiculeR_________',-1)

 200  format (' #',/,
     &  ' # osphiculeR.f: This routine can only be used in atoms',
     &  /,' #')
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
 50   format (' # Computing 2RDM over ',a)
 51   format (' # Computing E_ij^t coefficients')
      return
      end
