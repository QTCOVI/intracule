      module space_for_atdata
      save
      real(kind=8), allocatable, dimension (:)    :: ctenor
      real(kind=8), allocatable, dimension (:)    :: rhoat
      real(kind=8), allocatable, dimension (:)    :: xelect
      real(kind=8), allocatable, dimension (:)    :: elpop
      end module space_for_atdata
c-----------------------------------------------------------------------
      subroutine allocate_space_for_atdata (ncent)
      USE        space_for_atdata
      if (.not.allocated(ctenor)) then
        allocate (ctenor(ncent),stat=ier) 
        if (ier.ne.0) then
          write (0,*) 'atdata.f.f: Already allocated ctenor ?'
        endif
      endif
      if (.not.allocated(rhoat)) then
        allocate (rhoat(ncent),stat=ier) 
        if (ier.ne.0) then
          write (0,*) 'atdata.f.f: Already allocated rhoat ?'
        endif
      endif
      if (.not.allocated(xelect)) then
        allocate (xelect(ncent),stat=ier) 
        if (ier.ne.0) then
          write (0,*) 'atdata.f.f: Already allocated xelect ?'
        endif
      endif
      if (.not.allocated(elpop)) then
        allocate (elpop(ncent),stat=ier) 
        if (ier.ne.0) then
          write (0,*) 'atdata.f.f: Already allocated elpop ?'
        endif
      endif
      end subroutine allocate_space_for_atdata
c-----------------------------------------------------------------------
      subroutine deallocate_space_for_atdata
      USE        space_for_atdata
      if (allocated(ctenor)) then
        deallocate (ctenor,stat=ier)
        if (ier.ne.0) then
          write (0,*) 'atdata.f.f: Already deallocated ctenor ?'
        endif
      endif
      if (allocated(rhoat)) then
        deallocate (rhoat,stat=ier)
        if (ier.ne.0) then
          write (0,*) 'atdata.f.f: Already deallocated rhoat ?'
        endif
      endif
      if (allocated(xelect)) then
        deallocate (xelect,stat=ier)
        if (ier.ne.0) then
          write (0,*) 'atdata.f.f: Already deallocated xelect ?'
        endif
      endif
      if (allocated(elpop)) then
        deallocate (elpop,stat=ier)
        if (ier.ne.0) then
          write (0,*) 'atdata.f.f: Already deallocated elpop ?'
        endif
      endif
      end subroutine deallocate_space_for_atdata
