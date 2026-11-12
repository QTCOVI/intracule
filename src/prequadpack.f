! 
!-----------------------------------------------------------------------
!
      subroutine prequadpack ()
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
!-----Classification of GTFs in shells
!
      call timer (2,iprequadpack,'_prequadpack________',-1)
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
!-----Compute the E_t^ij, E_u^kl and E_v^mn coefficients
!
      write (uout,51)
      call computeEcoef ()
      call deallocate_space_for_inexcule  ()
      call timer (4,iprequadpack,'_prequadpack________',-1)
      return
!
!-----Formats (host scope; referenced by the writes above)
!
 51   format (' # Computing E_ij^t coefficients')
      return
      end subroutine prequadpack
