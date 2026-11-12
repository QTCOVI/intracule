module acumula_intracule
   implicit none

   ! Arrays accesibles desde cualquier rutina que haga 'use acumula_intracule'
   real(kind=8), allocatable :: rr(:)
   real(kind=8), allocatable :: valt(:), valcoul(:), valxc(:), valxcdf(:)
   integer(kind=4) :: nquad = 0
   integer(kind=4) :: maxquad = 100

contains

   subroutine acumula_init()
      allocate (rr(maxquad))
      allocate (valt(maxquad))
      allocate (valcoul(maxquad))
      allocate (valxc(maxquad))
      allocate (valxcdf(maxquad))
      nquad = 0
   end subroutine

   subroutine acumula_push (r,r1,r2,r3,r4)
      real(kind=8), intent(in) :: r,r1,r2,r3,r4
      real(kind=8), allocatable :: t1(:),t2(:),t3(:),t4(:),t5(:)
      if (nquad == maxquad) then
         allocate(t1(maxquad+100))
         allocate(t2(maxquad+100))
         allocate(t3(maxquad+100))
         allocate(t4(maxquad+100))
         allocate(t5(maxquad+100))
         t1(1:nquad) = rr(1:nquad)
         t2(1:nquad) = valt(1:nquad)
         t3(1:nquad) = valcoul(1:nquad)
         t4(1:nquad) = valxc(1:nquad)
         t5(1:nquad) = valxcdf(1:nquad)
         call move_alloc (t1, rr)
         call move_alloc (t2, valt)
         call move_alloc (t3, valcoul)
         call move_alloc (t4, valxc)
         call move_alloc (t5, valxcdf)
         maxquad = maxquad + 100
      end if
      nquad = nquad + 1
      rr(nquad)      = r
      valt(nquad)    = r1
      valcoul(nquad) = r2
      valxc(nquad)   = r3
      valxcdf(nquad) = r4
   end subroutine

   subroutine acumula_free()
      if (allocated(rr)) deallocate(rr, valt, valcoul, valxc, valxcdf)
      nquad = 0
   end subroutine

end module acumula_intracule
