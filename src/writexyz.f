      subroutine writexyz(filedat)

      use mod_prec, only : dp, ip
      use space_for_wfnbasis, only : xyz, charge
      include 'implicit.inc'
      include 'io.inc'
      include 'datatm.inc'
      include 'wfn.inc'

      character(len=*), intent(in) :: filedat

      integer(kind=ip), dimension(:), allocatable :: z      
      integer(kind=ip) :: i, ixyz, istat
      integer(kind=ip) :: find_free_unit
  
      allocate (z(ncent),stat=istat)
      z(:) = int(charge(:))

      ixyz = find_free_unit()
      open (ixyz,file=filedat)         
      write (ixyz,*) ncent
      write (ixyz,*)
      do i = 1,ncent
        write (ixyz,'(a2,x,3(f20.10,x))') namatm(z(i)), 
     &          xyz(i,:)*0.52917720859d0
      end do
      close (ixyz)

      deallocate (z,stat=istat)

      end subroutine writexyz
