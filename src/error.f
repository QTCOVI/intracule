      subroutine error (routine,message,errortype)
        use mod_prec, only : ip
        implicit none
        include 'io.inc'
        include 'error.inc'
        include 'parallel.inc'

        character(len=*), intent(in) :: routine !< Routine calling the error
        character(len=*), intent(in) :: message !< The message
        integer(kind=ip), intent(in) :: errortype !< Fatal, warning or info

        character(len=20) :: chtype

        if (errortype.eq.faterr) then
          chtype='ERROR'
        else if (errortype.eq.warning) then
          chtype='WARNING'
        else if (errortype.eq.noerr) then
          chtype='COMMENT'
        else
          chtype='UNKNOWN'
        endif
    
        if (errortype.eq.faterr) then
!$omp critical (io)
          write (uerr,100) trim(chtype),trim(routine),trim(message)
!$omp end critical (io)
           stop 1
        else if (errortype.eq.warning) then
!$omp critical (io)
          write (uout,101) trim(chtype),trim(routine),trim(message)
!$omp end critical (io)
!$omp atomic
          nwarms=nwarms+1
        else if (errortype.eq.noerr) then
!$omp critical (io)
          write (uout,101) trim(chtype),trim(routine),trim(message)
!$omp end critical (io)
!$omp atomic
          ncomms=ncomms+1
        endif

100     format(A,"(",A,"): ",A)
101     format(1x,"# ",A,"(",A,"): ",A)

      end subroutine error
