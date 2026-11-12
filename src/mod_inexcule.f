      module space_for_inexcule
      save

      real(kind=8), allocatable,dimension (:,:,:,:) :: ecoef
      real(kind=8), allocatable,dimension (:,:,:,:,:) :: etij
      real(kind=8), allocatable,dimension (:) :: psmall
      real(kind=8), allocatable,dimension (:) :: xmu
      real(kind=8), allocatable,dimension (:,:) :: prefac
      real(kind=8), allocatable,dimension (:) :: eprefac
      real(kind=8), allocatable,dimension (:,:) :: pbig

      end module space_for_inexcule
c
c-----------------------------------------------------------------------
c
      subroutine allocate_space_for_inexcule (npairshell,lx)
      USE        space_for_inexcule

      if (.not.allocated(ecoef)) then
        allocate (ecoef(-1:2*lx,-1:lx,-1:lx,3),stat=ier)
        if (ier.ne.0) then
          stop 'mod_inexcule.f: Cannot allocate ecoef()'
        endif
      endif

      if (.not.allocated(etij)) then
        allocate (etij(npairshell,0:lx+lx,0:lx,0:lx,3),stat=ier)
        if (ier.ne.0) then
          stop 'mod_inexcule.f: Cannot allocate etij()'
        endif
      endif

      if (.not.allocated(psmall)) then
        allocate (psmall(npairshell),stat=ier)
        if (ier.ne.0) then
          stop 'mod_inexcule.f: Cannot allocate psmall()'
        endif
      endif

      if (.not.allocated(xmu)) then
        allocate (xmu(npairshell),stat=ier)
        if (ier.ne.0) then
          stop 'mod_inexcule.f: Cannot allocate xmu()'
        endif
      endif

      if (.not.allocated(prefac)) then
        allocate (prefac(npairshell,4),stat=ier)
        if (ier.ne.0) then
          stop 'mod_inexcule.f: Cannot allocate prefac()'
        endif
      endif

      if (.not.allocated(eprefac)) then
        allocate (eprefac(npairshell),stat=ier)
        if (ier.ne.0) stop 'mod_inexcule.f: Cannot allocate eprefac()'
      endif
      if (.not.allocated(pbig)) then
        allocate (pbig(npairshell,3),stat=ier)
        if (ier.ne.0) then
          stop 'mod_inexcule.f: Cannot allocate pbig()'
        endif
      endif
c
      end subroutine allocate_space_for_inexcule
c
c-----------------------------------------------------------------------
c
      subroutine deallocate_space_for_inexcule ()
      USE        space_for_inexcule

      deallocate (ecoef,stat=ier)
      if (ier.ne.0) then
        stop 'mod_inexcule.f: Cannot deallocate ecoef()'
      endif
      deallocate (etij,stat=ier)
      if (ier.ne.0) then
        stop 'mod_inexcule.f: Cannot deallocate etij()'
      endif
      deallocate (psmall,stat=ier)
      if (ier.ne.0) then
        stop 'mod_inexcule.f: Cannot deallocate psmall()'
      endif
      deallocate (xmu,stat=ier)
      if (ier.ne.0) then
        stop 'mod_inexcule.f: Cannot deallocate xmu()'
      endif
      deallocate (prefac,stat=ier)
      if (ier.ne.0) then
        stop 'mod_inexcule.f: Cannot deallocate prefac()'
      endif
      deallocate (eprefac,stat=ier)
      if (ier.ne.0) stop 'mod_inexcule.f: Cannot deallocate eprefac()'
      deallocate (pbig,stat=ier)
      if (ier.ne.0) then
        stop 'mod_inexcule.f: Cannot deallocate pbig()'
      endif
c
      end subroutine deallocate_space_for_inexcule
