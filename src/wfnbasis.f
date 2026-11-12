      module space_for_wfnbasis
      save
      real    (kind=8)                              :: totalenergy,gamma
      real    (kind=8),allocatable,dimension (:,:)  :: xyz
      real    (kind=8),allocatable,dimension (:)    :: oexp
      real    (kind=8),allocatable,dimension (:)    :: occ,eorb,charge
      real    (kind=8),allocatable,dimension (:)    :: zcore,charget
      integer (kind=4),allocatable,dimension (:)    :: icen,ityp
      character(len=8),allocatable,dimension (:)    :: atnam
      logical         ,allocatable,dimension (:)    :: ghost
      end module space_for_wfnbasis
c-----------------------------------------------------------------------
      subroutine allocate_space_for_wfnbasis (nprims,nmo,ncent)
      USE        space_for_wfnbasis
c
      if (.not.allocated(xyz)) then
        allocate (xyz(ncent,3),stat=ier) 
        if (ier.ne.0) then
          stop 'wfnbasis.f: Already allocated xyz() array ?'
        endif
      endif
c
      if (.not.allocated(oexp)) then
        allocate (oexp(nprims),stat=ier) 
        if (ier.ne.0) then
          stop 'wfnbasis.f: Already allocated oexp() array ?'
        endif
      endif
c
      if (.not.allocated(occ)) then
        allocate (occ(nmo),stat=ier) 
        if (ier.ne.0) then
          stop 'wfnbasis.f: Already allocated occ() array ?'
        endif
      endif
c
      if (.not.allocated(eorb)) then
        allocate (eorb(nmo),stat=ier) 
        if (ier.ne.0) then
          stop 'wfnbasis.f: Already allocated eorb() array ?'
        endif
      endif
c
      if (.not.allocated(charge)) then
        allocate (charge(ncent),stat=ier) 
        if (ier.ne.0) then
          stop 'wfnbasis.f: Already allocated charge() array ?'
        endif
      endif
c
      if (.not.allocated(ghost)) then
        allocate (ghost(ncent),stat=ier) 
        if (ier.ne.0) then
          stop 'wfnbasis.f: Already allocated ghost() array ?'
        endif
      endif
c
      if (.not.allocated(zcore)) then
        allocate (zcore(ncent),stat=ier) 
        if (ier.ne.0) then
          stop 'wfnbasis.f: Already allocated zcore() array ?'
        endif
      endif
c
      if (.not.allocated(charget)) then
        allocate (charget(ncent),stat=ier) 
        if (ier.ne.0) then
          stop 'wfnbasis.f: Already allocated charget() array ?'
        endif
      endif
c
      if (.not.allocated(icen)) then
        allocate (icen(nprims),stat=ier) 
        if (ier.ne.0) then
          stop 'wfnbasis.f: Already allocated icen() array ?'
        endif
      endif
c
      if (.not.allocated(ityp)) then
        allocate (ityp(nprims),stat=ier) 
        if (ier.ne.0) then
          stop 'wfnbasis.f: Already allocated ityp() array ?'
        endif
      endif
c
      if (.not.allocated(atnam)) then
        allocate (atnam(ncent),stat=ier) 
        if (ier.ne.0) then
          stop 'wfnbasis.f: Already allocated atnam() array ?'
        endif
      endif
c
      end subroutine allocate_space_for_wfnbasis
c-----------------------------------------------------------------------
      subroutine deallocate_space_for_wfnbasis
      USE        space_for_wfnbasis
c
      if (allocated(xyz)) then
        deallocate (xyz,stat=ier) 
        if (ier.ne.0) then
          stop 'wfnbasis.f: Already deallocated xyz() array ?'
        endif
      endif
c
      if (allocated(oexp)) then
        deallocate (oexp,stat=ier) 
        if (ier.ne.0) then
          stop 'wfnbasis.f: Already deallocated oexp() array ?'
        endif
      endif
c
      if (allocated(occ)) then
        deallocate (occ,stat=ier) 
        if (ier.ne.0) then
          stop 'wfnbasis.f: Already deallocated occ() array ?'
        endif
      endif
c
      if (allocated(eorb)) then
        deallocate (eorb,stat=ier) 
        if (ier.ne.0) then
          stop 'wfnbasis.f: Already deallocated eorb() array ?'
        endif
      endif
c
      if (allocated(charge)) then
        deallocate (charge,stat=ier) 
        if (ier.ne.0) then
          stop 'wfnbasis.f: Already deallocated charge() array ?'
        endif
      endif
c
      if (allocated(ghost)) then
        deallocate (ghost,stat=ier) 
        if (ier.ne.0) then
          stop 'wfnbasis.f: Already deallocated ghost() array ?'
        endif
      endif
c
      if (allocated(zcore)) then
        deallocate (zcore,stat=ier) 
        if (ier.ne.0) then
          stop 'wfnbasis.f: Already deallocated zcore() array ?'
        endif
      endif
c
      if (allocated(charget)) then
        deallocate (charget,stat=ier) 
        if (ier.ne.0) then
          stop 'wfnbasis.f: Already deallocated charget() array ?'
        endif
      endif
c
      if (allocated(icen)) then
        deallocate (icen,stat=ier) 
        if (ier.ne.0) then
          stop 'wfnbasis.f: Already deallocated icen() array ?'
        endif
      endif
c
      if (allocated(ityp)) then
        deallocate (ityp,stat=ier) 
        if (ier.ne.0) then
          stop 'wfnbasis.f: Already deallocated ityp() array ?'
        endif
      endif
c
      if (allocated(atnam)) then
        deallocate (atnam,stat=ier) 
        if (ier.ne.0) then
          stop 'wfnbasis.f: Already allocated atnam() array ?'
        endif
      endif
c
      end subroutine deallocate_space_for_wfnbasis
