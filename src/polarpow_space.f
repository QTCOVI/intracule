      module space_for_polarpow
      save
      real(kind=8), allocatable, dimension (:,:) :: stp,ctp
      real(kind=8), allocatable, dimension (:,:) :: spp,cpp
      end module space_for_polarpow
c-----------------------------------------------------------------------
      subroutine allocate_space_for_polarpow (nang,maxpow)
      USE space_for_polarpow
      allocate (stp(nang,0:maxpow))
      allocate (ctp(nang,0:maxpow))
      allocate (spp(nang,0:maxpow))
      allocate (cpp(nang,0:maxpow))
      end subroutine allocate_space_for_polarpow
c-----------------------------------------------------------------------
      subroutine deallocate_space_for_polarpow ()
      USE space_for_polarpow
      deallocate (stp)
      deallocate (ctp)
      deallocate (spp)
      deallocate (cpp)
      end subroutine deallocate_space_for_polarpow
c-----------------------------------------------------------------------
