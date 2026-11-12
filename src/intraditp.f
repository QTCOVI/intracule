      subroutine intraditp (rmax,ngrid,verbose,screening,wfnfile)
!
!------------------------------------------------------------------------
!
      USE space_for_wfnbasis
      USE space_for_wfncoef
      USE space_for_primgto
      USE space_for_othershells
      USE space_for_inexcule
      USE mod_rdm
      USE grids
      implicit real(kind=8) (a-h,o-z)
      include   'param.inc'
      include   'wfn.inc'
      include   'corr.inc'
      include   'io.inc'
      include   'primgto.inc'
      include   'lengrec.inc'
!
!-----Input variables
!
      real    (kind=8),  intent(in) :: rmax
      integer (kind=4),  intent(in) :: ngrid
      logical,           intent(in) :: verbose
      character (len=*), intent(in) :: wfnfile

      real(kind=8), parameter :: xlog10 = 2.30258509299405D+00
      integer(kind=4), parameter :: nexpmax = 15
      real(kind=8), parameter :: xnlog = xlog10 * nexpmax

      real(kind=8)  abcd  (ngtoG,ngtoG,ngtoG,ngtoG)

      real(kind=8), allocatable,dimension (:,:,:,:,:) :: aoicule
      real(kind=8), allocatable,dimension (:,:,:,:)   :: moic
      real(kind=8), allocatable,dimension (:,:,:,:)   :: rdm2coul
      real(kind=8), allocatable,dimension (:,:,:,:)   :: rdm2df
      real(kind=8), allocatable,dimension (:,:,:,:)   :: aoictot
      real(kind=8), allocatable,dimension (:)         :: h
  
      real    (kind=8)  pxicule,grid(3)
      real    (kind=8)  muis,musqrt,chi,expchi,sqpi
      real    (kind=8)  xintra0(3)
      real    (kind=8)  efexyz(3)
      integer (kind=4)  p,q,pq,r,s,rs,di,dj,dk,dl,idx,jdx,kdx,ldx
      integer (kind=4)  a,b,c,d

      character (len=132) :: filei
      logical             :: exfil,screening
!
!-----Maximum L value in primitive cartesian Gaussians (CG) and its square
!
      lamx = maxlprim
      lamxsq = (lamx+1)*(lamx+2)/2
!
!-----Hermite polynomials
!
      allocate (h(0:4*lamx))
      write (uout,30)
 30   format (' # ',80('+'),/,' # Intracule densities in the defined',
     & ' cubic grid:  GRADIENT = NO   LAPLACIAN = NO',/,' # TOT=Full',
     & '  C=Coulomb   XC=Exc-Corr  XC-DF=DiracFock XC',/,' # ',80('+'))
!
!-----Compute 2-RDM in the basis of canonical MOs
!
      call timer (2,irdm2,'_rdm2MOs____________',-1)
      call compute_rdm (ndets,nelact,ncore,nmo,uout,uerr,nrecwfn,
     &    cciqa,mp2iqa,icorr,verbose,wfnfile)
      call timer (4,irdm2,'_rdm2MOs____________',-1)
!
!-----Compute the Coulombic part of the RDM2 and the Dirac-Fock part
!     of the XC RDM2
!
      allocate (rdm2coul(nmo,nmo,nmo,nmo),stat=ier)
      if (ier.ne.0) stop '# intraditp.f: Cannot allocate rdm2coul()'
      allocate (rdm2df(nmo,nmo,nmo,nmo),stat=ier)
      if (ier.ne.0) stop '# intraditp.f: Cannot allocate rdm2df()'
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

      call timer (2,kintraditp,'_intraditp__________',-1)
!
!-----Classify grid points
!
      call raditp (rmax,ngrid,verbose,uout)
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
!
!-----Open direct access file to store the intracule densities and derivatives
!
      lu = 1
      exfil = .true.
      filei = trim(wfnfile)//'.iculedat'
      lrec = npp*npp*npp*npp*RecLength
      do while (exfil)
        inquire (unit=lu,opened=exfil)
        if (.not.exfil) then
          open (unit    = lu,
     &          file    = trim(filei),
     &          access  = 'direct',
     &          recl    = lrec,
     &          form    = 'unformatted',
     &          iostat  = ier)
          if (ier.ne.0) then
            write (0,*) 'intraditp.f: Cannot open '//trim(filei)
            stop
          endif
        else
          lu=lu+1
        endif
      enddo

      allocate (aoicule(npp,npp,npp,npp,3),stat=ier)
      if (ier.ne.0) stop '# intraditp.f: Cannot allocate aoicule()'
      allocate (moic(nmo,nmo,nmo,nmo),stat=ier)
      if (ier.ne.0) stop '# intraditp.f: Cannot allocate moic() array'
      allocate (aoictot(npp,npp,npp,npp),stat=ier)
      if (ier.ne.0) stop '# intraditp.f: Cannot allocate aoictot()'

      irec = 0
      do k=1,3
        do l=1,npts(k)
          aoicule  (:,:,:,:,k)  = 0.0d+00
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

                  call hermipol (chi,h,4*lamx,labcd)

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
                          xintra0(k) = 0.0d+00
                          do it=0,ij
                            coefe1 = etij(pq,it,i,j,k)
                            do itp=0,ijp
                              coefe2 = etij(rs,itp,ip,jp,k)
                              p3 = coefe1 * coefe2 * prodboth
                              itt0 = it+itp
                              itt1 = itt0+1
                              sqt0 = musqrt**itt0 * ipow
                              xintra0(k) = xintra0(k) + p3*sqt0*h(itt0)
                            enddo
                            ipow   = -ipow
                          enddo
                          abcd  (ka,kb,kc,kd) = xintra0(k)
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
                  aoicule  (i1:i2,j1:j2,k1:k2,l1:l2,k) 
     &              = abcd  (1:di,1:dj,1:dk,1:dl)
                  do ia=i1,i2
                   do ja=j1,j2
                    do ka=k1,k2
                     do la=l1,l2
                      aoicule  (ja,ia,ka,la,k)=aoicule  (ia,ja,ka,la,k)
                      aoicule  (ia,ja,la,ka,k)=aoicule  (ia,ja,ka,la,k)
                      aoicule  (ja,ia,la,ka,k)=aoicule  (ia,ja,ka,la,k)
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
          write (lu,rec=irec+1) aoicule   (:,:,:,:,k)
          irec = irec + 1
        enddo ! l
      enddo ! k
!     write (uout,*) '# ',irec,' records written in "',trim(filei),'"'

      write (uout,9)
      do n=1,ngridxyz
        l1=indgrid(n,1)
        l2=indgrid(n,2)
        l3=indgrid(n,3)
        irecx = l1-1
        read (lu,rec=irecx+1) aoicule  (:,:,:,:,1)
        irecy = l2-1 + npts(1)
        read (lu,rec=irecy+1) aoicule  (:,:,:,:,2)
        irecz = l3-1 + npts(1) + npts(2)
        read (lu,rec=irecz+1) aoicule  (:,:,:,:,3)

        aoictot=aoicule(:,:,:,:,1)*aoicule(:,:,:,:,2)*aoicule(:,:,:,:,3)
c
c-------Transform intracule density matrix from the primitive Cartesian
c       Gaussians to the Canonical MO basis
c
        call moicule (aoictot,moic,npp,nmo,uout,verbose)
!       pxicule   = + sum(rdm2 * moic)
        pxicule   = + dot4d_kahan(rdm2,moic,nmo,nmo,nmo,nmo)
!       pxiculeC  = + sum(rdm2coul * moic)
        pxiculeC  = + dot4d_kahan(rdm2coul,moic,nmo,nmo,nmo,nmo)
        pxiculeXC = + pxicule - pxiculeC
!       pxiculeDF = - sum(rdm2df * moic) * 0.5d+00
        pxiculeDF = - dot4d_kahan(rdm2df,moic,
     &                nmo,nmo,nmo,nmo) * 0.5d+00
        grid(1) = gridp(l1,1) 
        grid(2) = gridp(l2,2) 
        grid(3) = gridp(l3,3) 
        write (uout,10) grid(1:3),pxicule,pxiculeC,pxiculeXC,pxiculeDF
      enddo
!
!-----Delete file containing the intracule densities and derivatives.
!
      close (lu,status='delete')

      deallocate (h,stat=ier)
      if (ier.ne.0) stop '# intraditp.f: Cannot deallocate h'
      deallocate (aoicule,stat=ier)
      if (ier.ne.0) stop '# intraditp.f: Cannot deallocate aoicule'
      deallocate (aoictot,stat=ier)
      if (ier.ne.0) stop '# intraditp.f: Cannot deallocate aoicule()'
      deallocate (moic,stat=ier)
      if (ier.ne.0) stop '# intraditp.f: Cannot deallocate moic()'
      deallocate (rdm2coul,stat=ier)
      if (ier.ne.0) stop '# intraditp.f: Cannot deallocate rdm2coul()'
      deallocate (rdm2df,stat=ier)
      if (ier.ne.0) stop '# intraditp.f: Cannot deallocate rdm2df()'

      call timer (4,kintraditp,'_intraditp__________',-1)

 9    format (' # ',4x,'X',16x,'Y',17x,'Z',17x,'TOT',13x,'C',15x,'XC',
     & 14x,'XC-DF',/,' # ',116('+'))
 10   format (1x,3(1x,F17.10),10(1x,E15.8))
      return
      end
