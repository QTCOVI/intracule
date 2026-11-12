      module space_for_othershells
      save
      integer(kind=4),allocatable,dimension (:)   :: center,lang,nzshell
      integer(kind=4),allocatable,dimension (:,:) :: ordp
      real(kind=8),allocatable,dimension (:)      :: zeta
      real(kind=8),allocatable,dimension (:,:)    :: prefactor
      logical, allocatable,dimension (:,:)        :: lshell
      end module space_for_othershells
c
c-----------------------------------------------------------------------
c
      subroutine allocate_space_for_othershells (nshells,ngtoG)
      USE        space_for_othershells
c
      if (.not.allocated(center)) then
        allocate (center(nshells),stat=ier) 
        if (ier.ne.0) stop 'othershell.f: Cannot allocate center()'
      endif
      if (.not.allocated(lang)) then
        allocate (lang(nshells),stat=ier) 
        if (ier.ne.0) stop 'othershell.f: Cannot allocate lang()'
      endif
      if (.not.allocated(nzshell)) then
        allocate (nzshell(nshells),stat=ier) 
        if (ier.ne.0) stop 'othershell.f: Cannot allocate nzshell()'
      endif
      if (.not.allocated(zeta)) then
        allocate (zeta(nshells),stat=ier) 
        if (ier.ne.0) stop 'othershell.f: Cannot allocate zeta()'
      endif
      if (.not.allocated(ordp)) then
        allocate (ordp(nshells,ngtoG),stat=ier) 
        if (ier.ne.0) stop 'othershell.f: Cannot allocate ordp()'
      endif
      if (.not.allocated(prefactor)) then
        allocate (prefactor(nshells,nshells),stat=ier) 
        if (ier.ne.0) stop 'othershell.f: Cannot allocate prefactor()'
      endif
      if (.not.allocated(lshell)) then
        allocate (lshell(nshells,nshells),stat=ier) 
        if (ier.ne.0) stop 'othershell.f: Cannot allocate lshell()'
      endif
      end subroutine allocate_space_for_othershells
c
c-----------------------------------------------------------------------
c
      subroutine deallocate_space_for_othershells ()
      USE        space_for_othershells
c
      if (allocated(center)) then
        deallocate (center,stat=ier) 
        if (ier.ne.0) stop 'othershell.f: Cannot deallocate center()'
      endif
      if (allocated(lang)) then
        deallocate (lang,stat=ier) 
        if (ier.ne.0) stop 'othershell.f: Cannot deallocate lang()'
      endif
      if (allocated(nzshell)) then
        deallocate (nzshell,stat=ier) 
        if (ier.ne.0) stop 'othershell.f: Cannot deallocate nzshell()'
      endif
      if (allocated(zeta)) then
        deallocate (zeta,stat=ier) 
        if (ier.ne.0) stop 'othershell.f: Cannot deallocate zeta()'
      endif
      if (allocated(ordp)) then
        deallocate (ordp,stat=ier) 
        if (ier.ne.0) stop 'othershell.f: Cannot deallocate ordp()'
      endif
      if (allocated(prefactor)) then
        deallocate (prefactor,stat=ier) 
        if (ier.ne.0) stop 'othershell.f: Cannot deallocate prefactor()'
      endif
      if (allocated(lshell)) then
        deallocate (lshell,stat=ier) 
        if (ier.ne.0) stop 'othershell.f: Cannot deallocate lshell()'
      endif
      end subroutine deallocate_space_for_othershells
c
c-----------------------------------------------------------------------
c
