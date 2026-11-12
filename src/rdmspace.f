
      module space_for_rdm1
      save
      real(kind=8),allocatable,dimension (:,:)  :: c1ea,c1eb,c1et
      real(kind=8),allocatable,dimension (:,:)  :: onerdm
      end module space_for_rdm1
c-----------------------------------------------------------------------
      subroutine allocate_space_for_rdm1 (nmo)
      USE        space_for_rdm1
c
      if (.not.allocated(c1ea)) then
        allocate (c1ea(nmo,nmo),stat=ier)
        if (ier.ne.0) stop 'rdmspace.f: Cannot allocate c1ea()'
      endif
      if (.not.allocated(c1eb)) then
        allocate (c1eb(nmo,nmo),stat=ier)
        if (ier.ne.0) stop 'rdmspace.f: Cannot allocate c1eb()'
      endif
      if (.not.allocated(c1et)) then
        allocate (c1et(nmo,nmo),stat=ier)
        if (ier.ne.0) stop 'rdmspace.f: Cannot allocate c1et()'
      endif
      if (.not.allocated(onerdm)) then
        allocate (onerdm(nmo,nmo),stat=ier)
        if (ier.ne.0) stop 'rdmspace.f: Cannot allocate onerdm()'
      endif
c
      end subroutine allocate_space_for_rdm1
c-----------------------------------------------------------------------
      subroutine deallocate_space_for_rdm1
      USE        space_for_rdm1
c
      if (allocated(c1ea)) then
        deallocate (c1ea,stat=ier)
        if (ier.ne.0) stop 'rdmspace.f: Cannot deallocate c1ea()'
      endif
      if (allocated(c1eb)) then
        deallocate (c1eb,stat=ier)
        if (ier.ne.0) stop 'rdmspace.f: Cannot deallocate c1eb()'
      endif
      if (allocated(c1et)) then
        deallocate (c1et,stat=ier)
        if (ier.ne.0) stop 'rdmspace.f: Cannot deallocate c1et()'
      endif
      if (allocated(onerdm)) then
        deallocate (onerdm,stat=ier)
        if (ier.ne.0) stop 'rdmspace.f: Cannot deallocate onerdm()'
      endif
c
      end subroutine deallocate_space_for_rdm1

      module space_for_rdm2
      save
      real(kind=8),allocatable,dimension (:,:)  :: c2e
      end module space_for_rdm2
c-----------------------------------------------------------------------
      subroutine allocate_space_for_rdm2 (nmo)
      USE        space_for_rdm2
c
      if (.not.allocated(c2e)) then
        allocate (c2e(nmo*(nmo+1)/2,nmo*(nmo+1)/2),stat=ier)
        if (ier.ne.0) stop 'rdmspace.f: Cannot allocate c2e()'
      endif
c
      end subroutine allocate_space_for_rdm2
c-----------------------------------------------------------------------
      subroutine deallocate_space_for_rdm2
      USE        space_for_rdm2
c
      if (allocated(c2e)) then
        deallocate (c2e,stat=ier)
        if (ier.ne.0) stop 'rdmspace.f: Cannot deallocate c2e()'
      endif
c
      end subroutine deallocate_space_for_rdm2
