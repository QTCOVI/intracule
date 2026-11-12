      subroutine icule1_mem (ider,verbose,screening,wfnfile,filegrid)
!
!------------------------------------------------------------------------
!
!-----This routine computes the intracule density, its gradient and its
!     Laplacian at the (x,y,z) points read in from file 'filegrid'.
!
!     For this purpose the following tasks a are performed.
!
!     1) The second order reduced density matrix in the canonical MO
!     basis (2-RDM) is computed in the 'compute_rdm' routine.
!
!     2) The (x,y,z) points read in from the 'filegrid' file are classified 
!     such that three (new) different lists of points (one for x, one for y, 
!     and one for z) are obtained.
!
!     3) The expansion coefficients of each product of two cartesian 
!     primitive Gaussian functions in terms of Hermite Gaussians are
!     computed in the 'computeEcoef' routine.
!
!     4) The intracule density is computed with a recently developed
!     algorithm.
!
!-----------------------------------------------------------------------
!
      USE space_for_wfnbasis
      USE space_for_wfncoef
      USE space_for_primgto
      USE space_for_othershells
      USE space_for_inexcule
      USE mod_rdm
      USE xyzgrid
      implicit real(kind=8) (a-h,o-z)
      include   'param.inc'
      include   'wfn.inc'
      include   'corr.inc'
      include   'io.inc'
      include   'primgto.inc'
!
!-----Input variables
!
      logical,           intent(in) :: verbose
      character (len=*), intent(in) :: wfnfile

      real(kind=8), parameter :: xlog10 = 2.30258509299405D+00
      integer(kind=4), parameter :: nexpmax = 15
      real(kind=8), parameter :: xnlog = xlog10 * nexpmax

      real(kind=8)  abcd  (ngtoG,ngtoG,ngtoG,ngtoG)
      real(kind=8)  abcdx (ngtoG,ngtoG,ngtoG,ngtoG)

      real(kind=8), allocatable,dimension (:,:,:,:,:,:) :: aoicule
      real(kind=8), allocatable,dimension (:,:,:,:,:,:) :: aoiculex
      real(kind=8), allocatable,dimension (:,:,:,:)     :: rho2prim
      real(kind=8), allocatable,dimension (:)           :: h
      real(kind=8), allocatable,dimension (:,:)         :: tmp 
  
      real    (kind=8)  pxicule,pxiculex(3),grid(3)
      real    (kind=8)  muis,musqrt,chi,expchi,sqpi
      real    (kind=8)  xintra0(3),xintra1(3)
      real    (kind=8)  efexyz(3)
      integer (kind=4)  p,q,pq,r,s,rs,di,dj,dk,dl,idx,jdx,kdx,ldx
      integer (kind=4)  a,b,c,d

      character (len = *) :: filegrid
      logical             :: screening
!
!-----Maximum L value in primitive cartesian Gaussians (CG) and its square
!
      lamx = maxlprim
      lamxsq = (lamx+1)*(lamx+2)/2
!
!-----Hermite polynomials
!
      allocate (h(0:4*lamx+ider))
!
!-----rho2prim is the 2-RDM matrix between primitive CGs
!
      allocate (rho2prim(nprims,nprims,nprims,nprims),stat=ier)
      if (ier.ne.0) then
        stop '# icule1_mem.f: Cannot allocate rho2prim() array'
      endif

      allocate (tmp(lamxsq,lamxsq),stat=ier)
      if (ier.ne.0) stop '# icule1_mem.f: Cannot allocate tmp() array'
!
      write (uout,'(a)') ' # +++++++++++++++++++++++++++++++++++++++'
      write (uout,'(a)') ' # Intracule densities in cubic-like grids'
      write (uout,'(a)') ' # GRADIENT = YES      LAPLACIAN = NO  '
      write (uout,'(a)') ' # Memory version is used'
      write (uout,'(a)') ' # +++++++++++++++++++++++++++++++++++++++'
!
!-----Compute 2-RDM
!
      call timer (2,irdm2,'_rdm2MOs____________',-1)
      call compute_rdm (ndets,nelact,ncore,nmo,uout,uerr,nrecwfn, 
     &    cciqa,mp2iqa,icorr,verbose,wfnfile)
      call timer (4,irdm2,'_rdm2MOs____________',-1)
c
c-----Compute 2-RDM in the basis of primitive Cartesian Gaussians
c
      call timer (2,irho2prim,'_rho2prim___________',-1)
      ione = 1
      npp = nprims
      call newrdm2prim (rho2prim)
      call timer (4,irho2prim,'_rho2prim___________',-1)

      call timer (2,kicule1_mem,'_icule1_mem_________',-1)
!
!-----Classify grid points
!
      call clasgrid (filegrid,verbose,uout)
!
!-----Classification of GTFs in shells
!
      nshells = sum(ngroup(1:ncent))
      call allocate_space_for_othershells (nshells,ngtoG)
      i=0
      do ic=1,ncent
        do m=1,ngroup(ic)
          i=i+1
          center(i)=ic
          nzshell(i)=nzexp(ic,m)
          ii=nuexp(ic,m,1)
          lang(i)=nlm(ityp(ii),1)+nlm(ityp(ii),2)+nlm(ityp(ii),3)
          zeta(i)=oexp(ii)
          ordp(i,1:nzshell(i))=nuexp(ic,m,1:nzexp(ic,m))
        enddo
      enddo 
!
!========================================================================
!
!-----Compute the E_t^ij, E_u^kl and E_v^mn coefficients
!
      call computeEcoef ()
!
!-----Loop over grid points X (k=1), Y(k=2), and Z (k=3)
!
      nxiculedim = max(npts(1),npts(2),npts(3))
!
!-----aoicule(a,b,c,d) is the I_abcd(R) intracule density between CGs
!
      npp = nprims
      allocate (aoicule(npp,npp,npp,npp,3,nxiculedim),stat=ier)
      if (ier.ne.0) then
        stop '# icule1_mem.f: Cannot allocate aoicule() array'
      endif
      allocate (aoiculex(npp,npp,npp,npp,3,nxiculedim),stat=ier)
      if (ier.ne.0) then
        stop '# icule1_mem.f: Cannot allocate aoiculex() array'
      endif
      aoicule(:,:,:,:,:,:)   = 0.0d+00
      aoiculex(:,:,:,:,:,:)  = 0.0d+00
      do k=1,3
        do l=1,npts(k)
          xxx = gridp(l,k)
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
                  dl = nzshell(s)
                  rs = rs + 1
!
!-----------------Screening of relevant integrals
!
                  abcd   = 0.0d+00
                  abcdx  = 0.0d+00
                  xtest1 = eprefac(pq) + eprefac(rs)
                  if (xtest1.gt.xnlog.and.screening) goto 1000
                  efexyz(:) = pbig(pq,:)-pbig(rs,:) + gridp(l,:)
                  prodz = psmall(pq)*psmall(rs)
                  addz = psmall(pq)+psmall(rs)
                  muis = prodz/addz
                  xtest2 = muis * dot_product(efexyz(:),efexyz(:))
                  if (xtest1+xtest2.gt.xnlog.and.screening) goto 1000

                  lcd  = lang(r) + lang(s)
                  labcd  = lab + lcd
                  efe = efexyz(k)
                  musqrt = sqrt(muis)
                  chi = musqrt * efe
                  expchi = exp(-chi*chi)
                  sqpi = sqrt(pi/addz)
                  sqpiexp = sqpi * expchi
                  ppp = prefac(pq,k) * prefac(rs,k)
                  prodboth = sqpiexp * ppp

                  call hermipol (chi,h,4*lamx+ider,labcd+ider)

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
                          ipow   = +1
                          ipowx  = -1
                          xintra0(k) = 0.0d+00
                          xintra1(k) = 0.0d+00
                          do it=0,ij
                            coefe1 = etij(pq,it,i,j,k)
                            do itp=0,ijp
                              coefe2 = etij(rs,itp,ip,jp,k)
                              p3 = coefe1 * coefe2 * prodboth
                              itt0 = it+itp
                              itt1 = itt0+1
                              sqt0 = musqrt**itt0 * ipow
                              sqt1 = musqrt**itt1 * ipowx
                              xintra0(k) = xintra0(k) + p3*sqt0*h(itt0)
                              xintra1(k) = xintra1(k) + p3*sqt1*h(itt1)
                            enddo
                            ipow   = -ipow
                            ipowx  = -ipowx
                          enddo
                          abcd (ka,kb,kc,kd) = xintra0(k)
                          abcdx(ka,kb,kc,kd) = xintra1(k)
                        enddo
                      enddo
                    enddo
                  enddo
 1000             continue
                  i1 = idx + 1
                  i2 = idx + di
                  j1 = jdx + 1
                  j2 = jdx + dj
                  k1 = kdx + 1
                  k2 = kdx + dk
                  l1 = ldx + 1
                  l2 = ldx + dl
                  aoicule  (i1:i2,j1:j2,k1:k2,l1:l2,k,l) 
     &              = abcd  (1:di,1:dj,1:dk,1:dl)
                  aoiculex (i1:i2,j1:j2,k1:k2,l1:l2,k,l) 
     &              = abcdx (1:di,1:dj,1:dk,1:dl)
                  do ia=i1,i2
                   do ja=j1,j2
                    do ka=k1,k2
                     do la=l1,l2
                      aoicule  (ja,ia,ka,la,k,l) 
     &                = aoicule  (ia,ja,ka,la,k,l)
                      aoicule  (ia,ja,la,ka,k,l) 
     &                = aoicule  (ia,ja,ka,la,k,l)
                      aoicule  (ja,ia,la,ka,k,l) 
     &                = aoicule  (ia,ja,ka,la,k,l)
                      aoiculex (ja,ia,ka,la,k,l) 
     &                = aoiculex (ia,ja,ka,la,k,l)
                      aoiculex (ia,ja,la,ka,k,l) 
     &                = aoiculex (ia,ja,ka,la,k,l)
                      aoiculex (ja,ia,la,ka,k,l) 
     &                = aoiculex (ia,ja,ka,la,k,l)
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
        enddo ! l
      enddo ! k

      write (uout,9)
      do n=1,ngridxyz
        n1 = indgrid(n,1)
        n2 = indgrid(n,2)
        n3 = indgrid(n,3)
        grid(1) = gridp(n1,1) 
        grid(2) = gridp(n2,2) 
        grid(3) = gridp(n3,3) 
        pxicule = sum(rho2prim 
     &    * aoicule (:,:,:,:,1,n1) 
     &    * aoicule (:,:,:,:,2,n2) 
     &    * aoicule (:,:,:,:,3,n3))
        pxiculex(1) = sum(rho2prim 
     &    * aoiculex(:,:,:,:,1,n1) 
     &    * aoicule (:,:,:,:,2,n2) 
     &    * aoicule (:,:,:,:,3,n3))
        pxiculex(2) = sum(rho2prim 
     &    * aoicule (:,:,:,:,1,n1) 
     &    * aoiculex(:,:,:,:,2,n2) 
     &    * aoicule (:,:,:,:,3,n3))
        pxiculex(3) = sum(rho2prim 
     &    * aoicule (:,:,:,:,1,n1) 
     &    * aoicule (:,:,:,:,2,n2) 
     &    * aoiculex(:,:,:,:,3,n3))
        write (uout,10) grid(:),pxicule,pxiculex(1:3)
      enddo

      deallocate (h,stat=ier)
      if (ier.ne.0) stop '# icule1_mem.f: Cannot deallocate h()'
      deallocate (aoicule,stat=ier)
      if (ier.ne.0) stop '# icule1_mem.f: Cannot deallocate aoicule()'
      deallocate (aoiculex,stat=ier)
      if (ier.ne.0) stop '# icule1_mem.f: Cannot deallocate aoiculex()'
      deallocate (tmp,stat=ier)
      if (ier.ne.0) stop '# icule1_mem.f: Cannot deallocate tmp()'
      deallocate (rho2prim,stat=ier)
      if (ier.ne.0) stop '# icule1_mem.f: Cannot deallocate rho2prim()'

      call timer (4,kicule1_mem,'_icule1_mem_________',-1)
 9    format (7x,'X',16x,'Y',17x,'Z',16x,'Intracule',7x,'Grad_x',
     & 10x,'Grad_y',10x,'Grad_z',/,1x,114('+'))
 10   format (3(1x,E17.10),10(1x,E15.8))
      return
      end
