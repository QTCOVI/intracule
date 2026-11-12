!
!-----------------------------------------------------------------------
!
      subroutine lebicule (verbose,wfnfile,rval,xicMO,xicMOCoul,
     &  xicXC,xicMODF,nrval,nang)
!
!------------------------------------------------------------------------
!
      USE space_for_wfnbasis
      USE space_for_wfncoef
      USE space_for_primgto
      USE space_for_othershells
      USE space_for_inexcule
      USE mod_rdm
      USE space_for_angles
      implicit real(kind=8) (a-h,o-z)
      include   'param.inc'
      include   'wfn.inc'
      include   'corr.inc'
      include   'io.inc'
      include   'primgto.inc'
      include   'parallel.inc'
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

      real(kind=8), allocatable,dimension (:,:,:,:,:,:) :: aac
      real(kind=8), allocatable,dimension (:,:,:,:)     :: aoicule
      real(kind=8), allocatable,dimension (:,:,:,:)     :: rdm2coul
      real(kind=8), allocatable,dimension (:,:,:,:)     :: rdm2diracfock
      real(kind=8), allocatable,dimension (:,:,:,:)     :: moic
      real(kind=8), allocatable,dimension (:,:,:,:)     :: abcd
      real(kind=8), allocatable,dimension (:)           :: dd
      real(kind=8), allocatable,dimension (:)           :: rp
      real(kind=8), allocatable,dimension (:,:)         :: stcppow
      real(kind=8), allocatable,dimension (:,:)         :: stsppow
      real(kind=8), allocatable,dimension (:,:)         :: ctpow
      integer (kind=4), allocatable, dimension (:,:,:,:,:) :: maxterm
      real(kind=8), parameter :: fourpi = 4.0d+00*pi

      real    (kind=8)  mu
      real    (kind=8)  angint
      real    (kind=8)  Iabcd
      real    (kind=8)  Tvec(3), RmT(3)
      integer (kind=4)  p,q,pq,r,s,rs,di,dj,dk,dl,idx,jdx,kdx,ldx
      integer (kind=4)  a
!
!-----array to store the coefficiente of the Hermite polynomials up to n=20
!
      integer(kind=4), parameter :: maxn = 20
      real   (kind=8) Chermite(0:maxn,0:maxn)
!
!     Precalcular potencias de stcp, stsp, ct
!
      allocate(stcppow(0:4*maxlprim,nang))
      allocate(stsppow(0:4*maxlprim,nang))
      allocate(ctpow(0:4*maxlprim,nang))
      do ia = 1, nang
        stcppow(0,ia) = 1.0d0
        stsppow(0,ia) = 1.0d0
        ctpow(0,ia) = 1.0d0
        do n = 1, 4*maxlprim
          stcppow(n,ia) = stcppow(n-1,ia) * stcp(ia)
          stsppow(n,ia) = stsppow(n-1,ia) * stsp(ia)
          ctpow(n,ia) = ctpow(n-1,ia) * ct(ia)
        end do
      end do

      if (maxn.lt.4*maxlprim) then
        write (uerr,1) maxn,4*maxlprim
        stop
      else
        allocate (dd(0:maxn))
        ixg = ngtox()
        allocate (abcd(ixg,ixg,ixg,ixg))
        allocate (aac(ixg,ixg,ixg,ixg,0:4*maxlprim,1:3))
        allocate (maxterm(ixg,ixg,ixg,ixg,1:3))
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
      call timer (2,ilebicule,'_lebicule___________',-1)
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
!-----------------------------------------------------------------------
!
!-----Compute the E_t^ij, E_u^kl and E_v^mn coefficients
!
      write (uout,51)
      call computeEcoef ()
c
c-----Powers of r
c
      allocate (rp(0:4*maxlprim))
!
!-----Print the r values of the grid
!
      write (uout,20)
      write (uout,21) (rval(i),i=1,nrval)
      write (uout,20)
!
!-----aoicule(a,b,c,d) is the spherically averaged I_abcd(R) intracule density 
!     between CGTFs
!
      if (.not.allocated(aoicule)) then
        allocate (aoicule(nprims,nprims,nprims,nprims),stat=ier)
        if (ier.ne.0) then
          stop '# lebicule.f: Cannot allocate aoicule() array'
        endif
      else
        stop '# lebicule.f: aoicule() array is already allocated'
      endif
      allocate (moic(nmo,nmo,nmo,nmo))
      moic = 0.0d+00
      write (uout,9)
!
!-----Loop over the different values of R
!
      do irval=1,nrval
!
!-------The distance at which the intracule density is computed
!
        rp(0) = 1.0d+00
        rr = rval(irval)
        do j=1,4*maxlprim 
          rp(j) = rr * rp(j-1)
        enddo
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
                sqrtmu = sqrt(mu)
                piz = sqrt(pi/addz)
!
!---------------this is the quantity T of the 'nintracule.pdf' document
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
c
c---------------exp{-mu*[rr-tis]**2} is the maximun value of 
c               exp{-mu*[\bm{r}-\bm{T}]**2}
c
                if (mu*(rr-tis)**2 .gt. 1000.0d+00) then
                  abcd(1:di,1:dj,1:dk,1:dl) = 0.0d+00
                  goto 222
                endif
! 
!---------------aac = The A_nx^{abcd},  A_ny^{abcd},  A_nz^{abcd} 
!               quantities of the 'nuevointracule.pdf' document
!               PARALLELIZED VERSION
!
                aac(:,:,:,:,:,:) = 0.0d+00
                do k = 1, 3
                  do ka = 1, di
                    i = nlm(ityp(ordp(p,ka)),k)
                    do kb = 1, dj
                      j = nlm(ityp(ordp(q,kb)),k)
                      do kc = 1, dk
                        ip = nlm(ityp(ordp(r,kc)),k)
                        do kd = 1, dl
                          jp = nlm(ityp(ordp(s,kd)),k)
                          maxterm(ka,kb,kc,kd,k) = i+j+ip+jp
                          do it = 0, i+j
                            pow = (-1.0d0)**it
                            coefe1 = etij(pq,it,i,j,k)*prefac(pq,k)
                            if (abs(coefe1).lt.1.0d-14) cycle
                            do itp = 0, ip+jp
                              coefe2 = etij(rs,itp,ip,jp,k)*prefac(rs,
     &                           k)
                              if (abs(coefe2).lt.1.0d-14) cycle
                              EtijEtij = coefe1 * coefe2
                              xmupow = sqrtmu**(it+itp) 
                              eedpi = EtijEtij * piz * xmupow * pow
                              if (abs(eedpi).lt.1.0d-14) cycle
                              tt = Tvec(k) 
                              call hermicart (mu,tt,it+itp,maxn,dd)
                              do ival = 0, it+itp
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
!---------------Evaluate the I_abcd(R) integral
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
                          do ny=0,ntermy
                            ayis = aac(ka,kb,kc,kd,ny,2)
                            if (abs(ayis).lt.1.0d-08) cycle
                            nxy = nx+ny
                            do nz=0,ntermz
                              nxyz = nx+ny+nz
                              azis = aac(ka,kb,kc,kd,nz,3)
                              Dnxnynz = axis * ayis * azis
!
! Bucle de integración angular
!
                              angint = 0.0d+00
!$OMP PARALLEL DO PRIVATE(ia,rp1,rp2,rp3,RmT,rpow_x,rpow_y,
!$OMP& rpow_z,xmuF2,eF2,xyzpow) REDUCTION(+:angint)
                              do ia = 1, nang
                                rp1 = rr * stcp(ia)
                                rp2 = rr * stsp(ia)
                                rp3 = rr * ct(ia)
                                RmT(1) = rp1 - Tvec(1)
                                RmT(2) = rp2 - Tvec(2)
                                RmT(3) = rp3 - Tvec(3)
                                
                                rpow_x = rp(nx) * stcppow(nx,ia)
                                rpow_y = rp(ny) * stsppow(ny,ia)
                                rpow_z = rp(nz) * ctpow(nz,ia)
                                
                                xmuF2 = mu 
     &                                * (RmT(1)*RmT(1)
     &                                + RmT(2)*RmT(2)
     &                                + RmT(3)*RmT(3))
                                eF2 = exp(-xmuF2)
                                xyzpow = eF2 * rpow_x * rpow_y * rpow_z
                                angint = angint + xyzpow * angw(ia)
                              end do
!$OMP END PARALLEL DO
                              
                              angint = angint / (4.0d+00*pi)
                              Iabcd = Iabcd + Dnxnynz * r2 * angint
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
        xicMOCoul(irval) = dot4d_kahan(rdm2coul,moic,
     &                     nmo,nmo,nmo,nmo)
        xicXC(irval) = xicMO(irval) - xicMOCoul(irval)
!       xicMODF(irval) = - sum(rdm2diracfock * moic(:,:,:,:)) * 0.5d+00
        xicMODF(irval) = - dot4d_kahan(rdm2diracfock,moic,
     &                     nmo,nmo,nmo,nmo) * 0.5d+00
        write (uout,100) rval(irval),xicMO(irval),xicMOCoul(irval),
     &    xicXC(irval),xicMODF(irval)
      enddo ! R values
!
!-----Output
!
 100  format (1x,F15.8,8(2x,E22.15))

      if (allocated(aoicule)) then
        deallocate (aoicule,stat=ier)
        if (ier.ne.0) then
          stop ' # lebicule.f: Cannot deallocate aoicule()'
        endif
      else
        stop ' # lebicule.f: aoicule() array is not allocated'
      endif
      deallocate (moic)
      deallocate (dd)
      deallocate (abcd)
      deallocate (aac)
      deallocate (rdm2coul)
      deallocate (rdm2diracfock)
      call timer (4,ilebicule,'_lebicule___________',-1)

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
