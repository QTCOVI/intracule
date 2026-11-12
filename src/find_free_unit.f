      function find_free_unit ()

        use mod_prec, only : ip
        include 'implicit.inc'
        include 'error.inc'
     
        integer(kind=ip), parameter :: max_unit = 240
        integer(kind=ip) :: find_free_unit
        integer(kind=ip) :: iunit
        logical :: opnd
     
        find_free_unit = -1    ! <-- valor por defecto si no se encuentra
        
        do iunit = 100, max_unit
          !inquire (unit=iunit,opened=opnd,exist=t_exist)
          inquire (unit=iunit,opened=opnd)
          !if (.not.opnd .and. t_exist) then
          if (.not.opnd) then
            find_free_unit = iunit
            return
          end if
        end do 
     
        call error ('find_free_unit','free unit not found ?!?',faterr)
     
        return
     
      end function find_free_unit
