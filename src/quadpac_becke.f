!
!-----------------------------------------------------------------------
!
      subroutine quadpac_becke (wfnfile,rval,xtot,xcoul,xxc,xdf)
!
!------------------------------------------------------------------------
!
      USE space_for_wfnbasis
      USE space_for_wfncoef
      USE space_for_primgto
      USE space_for_othershells
      USE space_for_inexcule
      USE space_for_angles
      USE mod_rdm
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
!-----Shared data.  The 2-RDMs and the (read-only) Lebedev grid in module
!     space_for_angles are used inside the parallel quartet loop; aoicule
!     is shared but written in disjoint blocks (one (p,q) bra-block per
!     task); all per-pair scratch is declared inside intra_pair.
!
      real(kind=8), allocatable,dimension (:,:,:,:), save  :: rdm2coul
      real(kind=8), allocatable,dimension (:,:,:,:), save  :: rdm2df
      real(kind=8), allocatable,dimension (:,:,:,:)     :: aoicule
      real(kind=8), allocatable,dimension (:,:,:,:)     :: moic
!
!-----Shell-offset table and (p,q) bra-pair map for the parallel quartet loop
!
      integer(kind=4), allocatable,dimension (:),save :: shoff,pmap,qmap

      integer (kind=4)  p,q,r,s
      integer (kind=4), save :: nang
      integer (kind=4), save :: ncapas
      integer (kind=4), save :: npairs
!
!-----array to store the coefficients of the Hermite polynomials up to n=20
!
      integer(kind=4), parameter :: maxn = 20
      real   (kind=8), save :: Chermite(0:maxn,0:maxn)

      call timer (2,iquadpac_becke,'_quadpac_becke______',-1)
      if (maxn.lt.4*maxlprim) then
        write (uerr,1) maxn,4*maxlprim
        stop
      else
        ixg = ngtox()
        Chermite(0:maxn,0:maxn) = 0.0d+00
        call hermitec (Chermite,maxn)
      endif
      nang = size(angw)
!
!-----Compute 2-RDM
!
      if (first_r_value) then
!       write (uout,50) 'Canonical Orbitals'
        call timer (2,irdm2,'_rdm2MOs____________',-1)
        call compute_rdm (ndets,nelact,ncore,nmo,uout,uerr,nrecwfn,
     &      cciqa,mp2iqa,icorr,verbose,wfnfile)
        call timer (4,irdm2,'_rdm2MOs____________',-1)
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
      endif
!
!-----------------------------------------------------------------------

      if (first_r_value) then
        call computeEcoef ()   ! creo que yo hace falta la primera vez
!
!-------Precompute the shell-offset table (the AO index at which each shell
!       starts) and the (p,q) bra-pair map, so the parallel quartet loop
!       below can recover all of its indices from a single flat task id,
!       with no sequential counters to serialise.
!
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
!       allocate (shoff(ncapas))
!       allocate (pmap(npairs))
!       allocate (qmap(npairs))
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


!     if (first_r_value) then
!      write (*,*) 'quadpac_becke: FIRST call, ncapas = ', ncapas, 
!    &        ' npairs = ', npairs, ' nang = ', nang
!     else
!      write (*,*) 'quadpac_becke: subsequent call, npairs = ', npairs,
!    &        ' allocated(pmap) = ', allocated(pmap)
!      endif
!      flush(6)


      allocate (aoicule(nprims,nprims,nprims,nprims))
      allocate (moic(nmo,nmo,nmo,nmo))
      rr = rval
      r2 = rr*rr
      aoicule = 0.0d+00
!$omp parallel do schedule(dynamic) default(shared) private(ipq)
      do ipq = 1,npairs
        call intra_pair (ipq,rr,r2,nang,aoicule)
      enddo
!$omp end parallel do
      call moicule (aoicule,moic,nprims,nmo,uout,verbose)
      xtot  = + dot4d_kahan(rdm2,moic,nmo,nmo,nmo,nmo)
      xcoul = + dot4d_kahan(rdm2coul,moic,nmo,nmo,nmo,nmo)
      xxc   = + xtot - xcoul
      xdf   = - dot4d_kahan(rdm2df,moic,nmo,nmo,nmo,nmo)*0.5d+00

      deallocate (aoicule)
      deallocate (moic)
      call timer (4,iquadpac_becke,'_quadpac_becke______',-1)
      return
!
!-----Formats (host scope; referenced by the writes above)
!
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
!     psmall, Chermite, the Lebedev grid in space_for_angles, shell data)
!     come in by host / module association.
!
!     This is the Lebedev-quadrature worker: the angular average
!     J_{nxnynz}(mu,R,T) of Eq. 46 of lebedev.pdf is evaluated by summing
!     over the angular grid instead of using the analytical theta/phi
!     integral tables.
!-----------------------------------------------------------------------
!
      subroutine intra_pair (ipq,rr,r2,nang,aoicule)
      implicit real(kind=8) (a-h,o-z)
      integer (kind=4), intent(in) :: ipq,nang
      real    (kind=8), intent(in) :: rr,r2
      real    (kind=8) aoicule(nprims,nprims,nprims,nprims)
      integer (kind=4) p,q,pq,r,s,rs,di,dj,dk,dl,idx,jdx,kdx,ldx
      integer (kind=4) i,j,k,ka,kb,kc,kd,ip,jp,it,itp,ival
      integer (kind=4) nx,ny,nz,ig,nn
      integer (kind=4) nmaxx,nmaxy,nmaxz,ntermx,ntermy,ntermz
      integer (kind=4) i1,i2,j1,j2,k1,k2,l1,l2,ia,ja,la
      integer (kind=4) lab,lcd
      real    (kind=8) mu,Iabcd,Tvec(3)
      real    (kind=8) dd(0:maxn)
      real(kind=8),allocatable,dimension(:,:,:,:)     :: abcd
      real(kind=8),allocatable,dimension(:,:,:,:,:,:) :: aac
      real(kind=8),allocatable,dimension(:,:,:)       :: Jtab
      integer(kind=4),allocatable,dimension(:,:,:,:,:):: maxterm
      real(kind=8),allocatable,dimension(:)           :: wexp
      real(kind=8),allocatable,dimension(:,:)         :: xpow,ypow,zpow

      allocate (abcd(ixg,ixg,ixg,ixg))
      allocate (aac(ixg,ixg,ixg,ixg,0:4*maxlprim,1:3))
      allocate (maxterm(ixg,ixg,ixg,ixg,1:3))
      allocate (Jtab(0:4*maxlprim,0:4*maxlprim,0:4*maxlprim))
      allocate (wexp(nang))
      allocate (xpow(nang,0:4*maxlprim))
      allocate (ypow(nang,0:4*maxlprim))
      allocate (zpow(nang,0:4*maxlprim))

      fourpiinv = 1.0d+00 / (4.0d+00*pi)
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
      lab = lang(p) + lang(q)
      do r = 1,ncapas
        dk = nzshell(r)
        kdx = shoff(r)
        do s = 1,r
          rs = r*(r-1)/2 + s
          dl = nzshell(s)
          ldx = shoff(s)
                lcd  = lang(r) + lang(s)
                Tvec(:) = pbig(rs,:)-pbig(pq,:)
                prodz = psmall(pq)*psmall(rs)
                addz = psmall(pq)+psmall(rs)
                mu = prodz/addz
                twomurr = 2.0d+00 * mu * rr
                twox = twomurr*Tvec(1)
                twoy = twomurr*Tvec(2)
                twoz = twomurr*Tvec(3)

                sqrtmu = sqrt(mu)
                piz = sqrt(pi/addz)
!
!---------------this is the quantity T of the 'nintracule.pdf' document
!
                tis = norm2(Tvec)
                T2 = tis * tis
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
!
!---------------aac = The A_nx^{abcd},  A_ny^{abcd},  A_nz^{abcd}
!               quantities of the 'nuevointracule.pdf' document
!               (the D^{abcd} coefficients of Eq. 45 of lebedev.pdf)
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
!---------------Angular average J_{nxnynz}(mu,R,T) of Eq. 46 of
!               lebedev.pdf by NUMERICAL Lebedev-Laikov quadrature.
!
!               Writing X=R*st*cp, Y=R*st*sp, Z=R*ct, the integrand of
!               Eq. 46 becomes (st*cp)^nx (st*sp)^ny ct^nz times the
!               Boltzmann factor exp[-mu(R-T)^2], so
!
!                 J_{nxnynz} = (1/4pi) SUM_g angw(g) * wexp(g)
!                              * (st*cp)^nx (st*sp)^ny ct^nz ,
!
!               with wexp(g) = exp(-mu(R^2+T^2) + 2muR.T_g), and
!               2muR.T_g = stcp(g)*twox + stsp(g)*twoy + ct(g)*twoz.
!               The exponent equals -mu*(R*u_g - T)^2 <= 0, so no
!               overflow protection is needed.  angw is normalised so
!               that SUM(angw)=4pi.  We build per-point power tables in
!               (st*cp), (st*sp) and ct once, then accumulate.
!
                do ig=1,nang
                  dotp = stcp(ig)*twox + stsp(ig)*twoy + ct(ig)*twoz
                  wexp(ig) = angw(ig)*fourpiinv*exp(-xmur2t2+dotp)
                  xpow(ig,0) = 1.0d+00
                  do nn=1,nmaxx
                    xpow(ig,nn) = xpow(ig,nn-1)*stcp(ig)
                  enddo
                  ypow(ig,0) = 1.0d+00
                  do nn=1,nmaxy
                    ypow(ig,nn) = ypow(ig,nn-1)*stsp(ig)
                  enddo
                  zpow(ig,0) = 1.0d+00
                  do nn=1,nmaxz
                    zpow(ig,nn) = zpow(ig,nn-1)*ct(ig)
                  enddo
                enddo
                do nx=0,nmaxx
                  do ny=0,nmaxy
                    do nz=0,nmaxz
                      acc = 0.0d+00
                      do ig=1,nang
                        acc = acc + wexp(ig)*xpow(ig,nx)
     &                            * ypow(ig,ny)*zpow(ig,nz)
                      enddo
                      Jtab(nx,ny,nz) = rr**(nx+ny+nz+2) * acc
                    enddo ! nz
                  enddo ! ny
                enddo ! nx
!
!---------------Evaluate the I_abcd(R) integral: contract the A_nx/A_ny/A_nz
!               coefficients with the angular table Jtab (Eq. 45)
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
              end do !s
            end do !r
!
!-----Free this task's per-pair scratch
!
      deallocate (abcd)
      deallocate (aac)
      deallocate (Jtab)
      deallocate (maxterm)
      deallocate (wexp)
      deallocate (xpow,ypow,zpow)
      return
      end subroutine intra_pair
      end
