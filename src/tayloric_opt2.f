!
!-----------------------------------------------------------------------
!
      subroutine tayloric_opt2 (verbose,wfnfile,rval,xicMO,xicMOCoul,
     &  xicXC,xicMODF,nrval,tmask,thresh)
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
      integer (kind=4),  intent(inout)  :: tmask

!
!-----Shared (read-only during the parallel R loop) data.  The 2-RDMs and
!     the angular-integral tables are built once in the host below and only
!     read inside intra_oneR; all per-R scratch is declared in intra_oneR.
!
      real(kind=8), allocatable,dimension (:,:,:,:)     :: rdm2coul
      real(kind=8), allocatable,dimension (:,:,:,:)     :: rdm2diracfock
      real(kind=8), allocatable,dimension (:,:)         :: thetaval
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
!-----maximum Taylor order for the e^{2*mu*R.T} coupling expansion and
!     table of log(k!) used by the adaptive, overflow-safe summation
!
      integer(kind=4), parameter :: kmaxcap = 400
      real   (kind=8) logfac(0:kmaxcap)


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
      call timer (2,itayloric,'_tayloric_opt2______',-1)
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
!-----We compute the INT_0^2\pi sin(x)^(2k) cos(x)^(2l) dx integrals
!
      nxmaxphi = 4*maxlprim + kmaxcap + 2
      allocate (yntphi(0:nxmaxphi,0:nxmaxphi))
      call phiintegrals (nxmaxphi,nxmaxphi,yntphi)
!
!-----------------------------------------------------------------------
!
!-----Compute the E_t^ij, E_u^kl and E_v^mn coefficients
!
      write (uout,51)
      call computeEcoef ()

      nxmaxtheta = 4*maxlprim + kmaxcap
      nxy_max = 2 * nxmaxtheta   ! porque nx+ny puede llegar hasta 2*nxmaxtheta
      allocate(thetaval(0:nxy_max, 0:nxmaxtheta+2))
      call thetaintegrals (nxy_max,nxmaxtheta+2,thetaval)
!
!-----Table of log(k!) for the overflow-safe Taylor summation below
!
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
!     between CGTFs; it (and the rest of the per-R scratch) is allocated
!     per worker inside intra_oneR so each thread has its own copy.
!
      write (uout,9)
!
!-----Loop over the different values of R.  Each R is independent (aoicule
!     is rebuilt from scratch and the densities are stored by index), so the
!     loop is parallelised over irval with OpenMP.  All per-R scratch lives
!     in the internal routine intra_oneR and is therefore thread-private;
!     the read-only tables come in by host association and are shared.
!
!$omp parallel do default(shared) schedule(dynamic)
      do irval=1,nrval
        call intra_oneR (irval)
      enddo
!$omp end parallel do
!
!-----Print the densities in grid order (kept out of the parallel region so
!     the output is ordered and free of inter-thread interleaving)
!
      do irval=1,nrval
        write (uout,100) rval(irval),xicMO(irval),xicMOCoul(irval),
     &    xicXC(irval),xicMODF(irval)
      enddo

      deallocate (rdm2coul)
      deallocate (rdm2diracfock)
      call timer (4,itayloric,'_tayloric_opt2______',-1)
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
 50   format (' # Computing 2RDM over ',a)
 51   format (' # Computing E_ij^t coefficients')
 100  format (1x,F15.8,8(2x,E22.15))

      contains
!
!-----------------------------------------------------------------------
!     Per-R worker: builds aoicule(R), transforms it to the MO basis, and
!     stores the four intracule densities for grid point irval.  Every
!     array declared here is local, hence thread-private under the OpenMP
!     loop above; the read-only tables (etij, prefac, pbig, psmall, yntphi,
!     thetaval, Chermite, logfac, the shell data and rdm2*) come in by host
!     association and are only read, so they are safely shared.
!-----------------------------------------------------------------------
!
      subroutine intra_oneR (irval)
      implicit real(kind=8) (a-h,o-z)
      integer (kind=4), intent(in) :: irval
      integer (kind=4) p,q,pq,r,s,rs,di,dj,dk,dl,idx,jdx,kdx,ldx
      integer (kind=4) i,j,k,ka,kb,kc,kd,ip,jp,it,itp,ival,kk
      integer (kind=4) nx,ny,nz,nxk,nyk,nzk,kx,ky,kz,k0,kup,kwmax
      integer (kind=4) nmaxx,nmaxy,nmaxz,ntermx,ntermy,ntermz
      integer (kind=4) i1,i2,j1,j2,k1,k2,l1,l2,ia,ja,la
      integer (kind=4) lab,lcd,kxmax,kymax,kzmax
      real    (kind=8) mu,Iabcd,Tvec(3)
      real    (kind=8) wx(0:kmaxcap),wy(0:kmaxcap),wz(0:kmaxcap)
      real    (kind=8) dd(0:maxn)
      logical ejex,ejey,ejez
      real(kind=8),allocatable,dimension(:,:,:,:)     :: aoicule,moic
      real(kind=8),allocatable,dimension(:,:,:,:)     :: abcd
      real(kind=8),allocatable,dimension(:,:,:,:,:,:) :: aac
      real(kind=8),allocatable,dimension(:,:,:)       :: Jtab
      integer(kind=4),allocatable,dimension(:,:,:,:,:):: maxterm

      allocate (abcd(ixg,ixg,ixg,ixg))
      allocate (aac(ixg,ixg,ixg,ixg,0:4*maxlprim,1:3))
      allocate (maxterm(ixg,ixg,ixg,ixg,1:3))
      allocate (Jtab(0:4*maxlprim,0:4*maxlprim,0:4*maxlprim))
      allocate (aoicule(nprims,nprims,nprims,nprims))
      allocate (moic(nmo,nmo,nmo,nmo))
      moic = 0.0d+00
!
!-------The distance at which the intracule density is computed
!
        rr = rval(irval)
*       rr = 0.82384615
        r2 = rr*rr
        aoicule = 0.0d+00
        pq = 0
        idx = 0
        do p = 1,nshells
          di = nzshell(p)
          jdx = 0
          do q = 1,p
            dj = nzshell(q)
            pq = pq + 1
            pf1 = prefac(pq,1)*prefac(pq,2)*prefac(pq,3)
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
                Tvec(:) = pbig(rs,:)-pbig(pq,:) 
                prodz = psmall(pq)*psmall(rs)
                addz = psmall(pq)+psmall(rs)
                mu = prodz/addz
                twomurr = 2.0d+00 * mu * rr
                twox = twomurr*Tvec(1)
                twoy = twomurr*Tvec(2)
                twoz = twomurr*Tvec(3)

                kxmax = tmask
                ejex = .false.
                if (twox.gt.0.0d+00) then
                  xtwox = log(+twox)
                  signox = +1.0d+00
                elseif (twox.lt.0.0d+00) then
                  xtwox = log(-twox)
                  signox = -1.0d+00
                else
                  xtwox = 0.0d+00
                  signox = +1.0d+00
                  kxmax = 0
                  ejex = .true.
                endif

                kymax = tmask
                ejey = .false.
                if (twoy.gt.0.0d+00) then
                  xtwoy = log(+twoy)
                  signoy = +1.0d+00
                elseif (twoy.lt.0.0d+00) then
                  xtwoy = log(-twoy)
                  signoy = -1.0d+00
                else
                  xtwoy = 0.0d+00
                  signoy = +1.0d+00
                  kymax = 0
                  ejey = .true.
                endif

                ejez = .false.
                kzmax = tmask
                if (twoz.gt.0.0d+00) then
                  xtwoz = log(+twoz)
                  signoz = +1.0d+00
                elseif (twoz.lt.0.0d+00) then
                  xtwoz = log(-twoz)
                  signoz = -1.0d+00
                else
                  xtwoz = 0.0d+00
                  signoz = +1.0d+00
                  kzmax = 0
                  ejez = .true.
                endif

                sqrtmu = sqrt(mu)
                piz = sqrt(pi/addz)
!
!---------------this is the quantity T of the 'nintracule.pdf' document
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
                if (xmur2t2 .gt. 1000.0+00) then
                  abcd(1:di,1:dj,1:dk,1:dl) = 0.0d+00
                  goto 222
                endif
                xdamp = -xmur2t2 - log(4.0d+00*pi)
!
!---------------aac = The A_nx^{abcd},  A_ny^{abcd},  A_nz^{abcd} 
!               quantities of the 'nuevointracule.pdf' document
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
     &                          dd(ival) + aac(ka,kb,kc,kd,ival,k) 
                              enddo
                            enddo
                          enddo
                        enddo
                      enddo
                    enddo
                  enddo
                enddo
!
!---------------Maximum Cartesian-Hermite orders actually needed by any of
!               the di*dj*dk*dl Cartesian components of this quartet
!
                nmaxx = maxval(maxterm(1:di,1:dj,1:dk,1:dl,1))
                nmaxy = maxval(maxterm(1:di,1:dj,1:dk,1:dl,2))
                nmaxz = maxval(maxterm(1:di,1:dj,1:dk,1:dl,3))
!
!---------------Damped Taylor weights (2*mu*R*T_i)^k/k! * exp(xdamp/3),
!               built in log space (so neither (2muRT)^k nor 1/k!
!               overflows) and reused for every (nx,ny,nz) and every
!               Cartesian component below, so the exp() calls leave the
!               innermost loops.  They are filled lazily up to the highest
!               order kwmax the adaptive loop below actually reaches (which
!               for s/p / loosely-coupled quartets is small), instead of
!               eagerly up to kmaxcap.
!
                xdamp3 = xdamp / 3.0d+00
                base = abs(twox)+abs(twoy)+abs(twoz)
                kwmax = -1
!
!---------------Taylor expansion of the coupling factor e^{2*mu*R.T}
!               (Eq. 52 of fusion3.pdf).  The Gaussian parities fix
!               kx==nx, ky==ny, kz==nz (mod 2), so only total orders k of
!               the same parity as nx+ny+nz contribute; we step k by 2
!               from k0 and sum the multinomial terms p+q+r=k, increasing
!               k adaptively until the contribution falls below thresh.
!
!               The coupling integral depends only on (nx,ny,nz) and the
!               quartet geometry, NOT on the Cartesian indices, so it is
!               tabulated once in Jtab(nx,ny,nz) and reused below.
!
                do nx=0,nmaxx
                  do ny=0,nmaxy
                    do nz=0,nmaxz
                      k0 = mod(nx+ny+nz,2)
                      if (base.eq.0.0d+00) then
                        kup = k0
                      else
                        kup = kmaxcap
                      endif
                      term = 0.0d+00
                      do k = k0, kup, 2
!
!-----------------------extend the weight tables on demand (each order is
!                       built exactly once per quartet)
!
                        if (k.gt.kwmax) then
                          do kk = kwmax+1, k
                            wx(kk) = exp(xdamp3 + kk*xtwox - logfac(kk))
                            if (signox.lt.0.0d0.and.mod(kk,2).eq.1)
     &                        wx(kk) = -wx(kk)
                            wy(kk) = exp(xdamp3 + kk*xtwoy - logfac(kk))
                            if (signoy.lt.0.0d0.and.mod(kk,2).eq.1)
     &                        wy(kk) = -wy(kk)
                            wz(kk) = exp(xdamp3 + kk*xtwoz - logfac(kk))
                            if (signoz.lt.0.0d0.and.mod(kk,2).eq.1)
     &                        wz(kk) = -wz(kk)
                          enddo
                          kwmax = k
                        endif
                        contrib = 0.0d+00
                        do kx = 0, k
                          if (ejex .and. kx.gt.0) cycle
                          nxk = nx+kx
                          if (mod(nxk,2).ne.0) cycle
                          do ky = 0, k-kx
                            if (ejey .and. ky.gt.0) cycle
                            kz = k-kx-ky
                            if (ejez .and. kz.gt.0) cycle
                            nyk = ny+ky
                            nzk = nz+kz
                            if (mod(nyk,2).ne.0) cycle
                            if (mod(nzk,2).ne.0) cycle
                            thetaint = thetaval(nxk+nyk,nzk)
                            if (thetaint.eq.0.0d+00) cycle
                            phiint = yntphi(nxk,nyk)
                            if (phiint.eq.0.0d+00) cycle
                            contrib = contrib + wx(kx)*wy(ky)*wz(kz)
     &                                 * thetaint * phiint
                          enddo  ! ky
                        enddo  ! kx
                        term = term + contrib
                        if (dble(k).ge.base .and.
     &                    abs(contrib).le.thresh*abs(term)) exit
                      enddo  ! k (adaptive Taylor order)
                      Jtab(nx,ny,nz) = rr**(nx+ny+nz+2) * term
                    enddo ! nz
                  enddo ! ny
                enddo ! nx
!
!---------------Evaluate the I_abcd(R) integral: contract the A_nx/A_ny/A_nz
!               coefficients with the precomputed coupling table Jtab
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
                              Iabcd = Iabcd
     &                          + axis*ayis*azis * Jtab(nx,ny,nz)
                            enddo ! nz
                          enddo ! ny
                        enddo  ! nx
                        abcd(ka,kb,kc,kd) = Iabcd
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
!       xicMO(irval) = sum(rdm2 * moic(:,:,:,:))
        xicMO(irval) = dot4d_kahan(rdm2,moic,nmo,nmo,nmo,nmo)
!       xicMOCoul(irval) = sum(rdm2coul * moic(:,:,:,:))
        xicMOCoul(irval) = dot4d_kahan(rdm2coul,moic,nmo,nmo,nmo,nmo)
        xicXC(irval) = xicMO(irval) - xicMOCoul(irval)
!       xicMODF(irval) = - sum(rdm2diracfock * moic(:,:,:,:)) * 0.5d+00
        xicMODF(irval) = - dot4d_kahan(rdm2diracfock,moic,
     &                     nmo,nmo,nmo,nmo) * 0.5d+00
*       stop 'aqui'
!
!-------Free this worker's per-R scratch
!
      deallocate (aoicule)
      deallocate (moic)
      deallocate (abcd)
      deallocate (aac)
      deallocate (Jtab)
      deallocate (maxterm)
      return
      end subroutine intra_oneR
      end
