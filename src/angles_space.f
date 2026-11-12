      module space_for_angles
      save
      real(kind=8), allocatable, dimension (:) :: ct,st,cp,sp,angw
      real(kind=8), allocatable, dimension (:) :: stsp,stcp
      end module space_for_angles
c-----------------------------------------------------------------------
      subroutine allocate_space_for_angles (nang)
      USE        space_for_angles
      allocate (ct(nang))
      allocate (st(nang))
      allocate (cp(nang))
      allocate (sp(nang))
      allocate (angw(nang))
      allocate (stsp(nang))
      allocate (stcp(nang))
      end subroutine allocate_space_for_angles
c-----------------------------------------------------------------------
      subroutine deallocate_space_for_angles ()
      USE        space_for_angles
      deallocate (ct,st,cp,sp,angw,stsp,stcp)
      end subroutine deallocate_space_for_angles
