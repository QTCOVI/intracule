      module space_for_corr
      save
      integer,     allocatable, dimension (:)   :: ialpha,ibeta,pind
      integer,     allocatable, dimension (:,:) :: ipak
      real(kind=8),allocatable, dimension (:)   :: cidet
      end module space_for_corr
c-----------------------------------------------------------------------
      subroutine allocate_space_for_corr (nmo)
      USE        space_for_corr
****  parameter  (maxactive=200)
      integer(kind=4) i,j,ij,nmo
c
      if (.not.allocated(ipak)) then
        allocate (ipak(nmo*(nmo+2),2),stat=ier) 
        if (ier.ne.0) then
          stop 'corrarea.f: Cannot allocate ipak()'
        else
          ij=0
          do i=1,nmo
            do j=1,i
              ij=ij+1
              ipak(ij,1)=i
              ipak(ij,2)=j
            enddo
          enddo
        endif
      endif
c
      if (.not.allocated(ialpha)) then
        allocate (ialpha(nmo),stat=ier) 
        if (ier.ne.0) stop 'corrarea.f: Cannot allocate ialpha()'
      endif
c
      if (.not.allocated(ibeta)) then
        allocate (ibeta(nmo),stat=ier) 
        if (ier.ne.0) stop 'corrarea.f: Cannot allocate ibeta()'
      endif
c
      if (.not.allocated(pind)) then
        allocate (pind(nmo*(nmo+2)),stat=ier) 
        if (ier.ne.0) stop 'corrarea.f: Cannot allocate pind()'
      endif
c
      if (.not.allocated(cidet)) then
        allocate (cidet(0:nmo+nmo),stat=ier) 
        if (ier.ne.0) stop 'corrarea.f: Cannot allocate cidet()'
      endif
c
      end subroutine allocate_space_for_corr
c-----------------------------------------------------------------------
      subroutine deallocate_space_for_corr
      USE        space_for_corr
c
      if (allocated(ipak)) then
        deallocate (ipak,stat=ier) 
        if (ier.ne.0) stop 'corrarea.f: Cannot deallocate ipak()'
      endif
c
      if (allocated(ialpha)) then
        deallocate (ialpha,stat=ier) 
        if (ier.ne.0) stop 'corrarea.f: Cannot deallocate ialpha()'
      endif
c
      if (allocated(ibeta)) then
        deallocate (ibeta,stat=ier) 
        if (ier.ne.0) stop 'corrarea.f: Cannot deallocate ibeta()'
      endif
c
      if (allocated(pind)) then
        deallocate (pind,stat=ier) 
        if (ier.ne.0) stop 'corrarea.f: Cannot deallocate pind()'
      endif
c
      if (allocated(cidet)) then
        deallocate (cidet,stat=ier) 
        if (ier.ne.0) stop 'corrarea.f: Cannot deallocate cidet()'
      endif
c
      end subroutine deallocate_space_for_corr
