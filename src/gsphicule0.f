!
!-----------------------------------------------------------------------
!
      subroutine gsphicule0 (verbose,wfnfile,rval,
     &  xicMO,xicMOCoul,xicXC,xicMODF,nrval)
!
!------------------------------------------------------------------------
!
!    Computes the intracular density when only 1s-type GTFs are involved
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
      real (kind=8),     intent(out) :: xicMO(nrval)
      real (kind=8),     intent(out) :: xicMOCoul(nrval)
      real (kind=8),     intent(out) :: xicXC(nrval)
      real (kind=8),     intent(out) :: xicMODF(nrval)

      real(kind=8), parameter :: exp_limit = 700.0d0
      real(kind=8), parameter :: alpha_limit = 1.0d-12

      real(kind=8), allocatable,dimension (:,:,:,:)     :: aoicule
      real(kind=8), allocatable,dimension (:,:,:,:)     :: rdm2coul
      real(kind=8), allocatable,dimension (:,:,:,:)     :: rdm2diracfock
      real(kind=8), allocatable,dimension (:,:,:,:)     :: moic
      real    (kind=8)  mu
      real    (kind=8)  Iabcd
      real    (kind=8)  txyz(3)
      integer (kind=4)  p,q,pq,r,s,rs
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
      call timer (2,isphiculeR,'_gsphicule0_________',-1)
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
      call computeEcoef ()
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
      allocate (aoicule(nprims,nprims,nprims,nprims))
      allocate (moic(nmo,nmo,nmo,nmo))

      write (uout,9)
!
!-----Loop over the different values of R
!
      do irval=1,nrval
        rr = rval(irval)
        r2 = rr * rr
        aoicule = 0.0d+00
        pq = 0
        do p = 1,nshells
          do q = 1,p
            pq = pq + 1
            pf1 = prefac(pq,1)*prefac(pq,2)*prefac(pq,3)
            rs = 0
            do r = 1,nshells
              do s = 1,r
                rs = rs + 1
                prodz = psmall(pq)*psmall(rs)
                addz = psmall(pq)+psmall(rs)
                mu = prodz/addz
                sqrtmu = sqrt(mu)
                piz = sqrt(pi/addz)
                txyz(:) = pbig(rs,:)-pbig(pq,:) 
                tis = norm2(txyz)
                t2 = tis * tis
                pf2 = prefac(rs,1)*prefac(rs,2)*prefac(rs,3)
                xk = pf1 * pf2 * piz**3
                alpha = 2.0d+00 * mu * rr * tis
                twoalpha = alpha+alpha
                beta = mu * (r2+t2)
                ebeta = exp(-beta)

                if (abs(alpha) .lt. alpha_limit) then
                  ! Caso alpha ≈ 0: usar límite sinh(alpha)/alpha → 1
                  Iabcd = xk * r2 * exp(-mu*(r2+t2))
                else
                  ! Calcular los dos términos exponenciales
                  term1 = -mu*(rr-tis)**2
                  term2 = -mu*(rr+tis)**2
                  if (term1 .gt. exp_limit) then
                    Iabcd = 0.0d0
                  else if (term1 .lt. -exp_limit) then
                    Iabcd = 0.0d0
                  else
                    Iabcd = xk * r2 * (exp(term1)-exp(term2))/twoalpha
                  endif
                endif
                aoicule(p,q,r,s) = Iabcd
                aoicule(p,q,s,r) = Iabcd
                aoicule(q,p,r,s) = Iabcd
                aoicule(q,p,s,r) = Iabcd
              enddo
            enddo
          enddo
        enddo
        call moicule (aoicule,moic,nprims,nmo,uout,verbose)
!       xicMO(irval) = sum(rdm2 * moic(:,:,:,:))
        xicMO(irval) = dot4d_kahan(rdm2,moic,nmo,nmo,nmo,nmo)
!       xicMOCoul(irval) = sum(rdm2coul * moic(:,:,:,:))
        xicMOCoul(irval) = dot4d_kahan(rdm2coul,moic,nmo,nmo,nmo,nmo)
        xicXC(irval) = xicMO(irval) - xicMOCoul(irval)
!       xicMODF(irval) = -sum(rdm2diracfock*moic(:,:,:,:))/2.0D0
        xicMODF(irval) = -dot4d_kahan(rdm2diracfock,moic,
     &                   nmo,nmo,nmo,nmo)/2.0D0
        write (uout,100) rval(irval),xicMO(irval),xicMOCoul(irval),
     &    xicXC(irval),xicMODF(irval)
      enddo ! R values
!
!-----Output
!
 100  format (1x,F15.8,8(2x,E22.15))

      deallocate (aoicule)
      deallocate (moic)
      deallocate (rdm2coul)
      deallocate (rdm2diracfock)
      call timer (4,isphiculeR,'_gsphicule0_________',-1)

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
      return
      end
