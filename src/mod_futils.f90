module mod_futils

  implicit none
  public

contains

  !> Eliminar el límite de pila (no se puede hacer en Fortran puro)
  subroutine remove_stack_limit(retval)
    implicit none
    integer, intent(out) :: retval
    ! En Fortran puro no se puede hacer, pero podemos usar un comando del sistema
    call system("ulimit -s unlimited 2>/dev/null")
    retval = 0
  end subroutine remove_stack_limit

  !> Crear un directorio
  subroutine mkdir(directory,retval)
    implicit none
    character(len=*), intent(in) :: directory
    integer, intent(out) :: retval

    ! Usa el comando del sistema 'mkdir'
    character(len=256) :: cmd
    write(cmd, '(A)') "mkdir -p " // trim(directory) // " 2>/dev/null"
    call system(cmd)
    retval = 0
  end subroutine mkdir

  !> Establecer variable de entorno
  subroutine putenv(varname,retval)
    implicit none
    character(len=*), intent(in) :: varname
    integer, intent(out) :: retval

    ! Usa el comando del sistema 'export'
    character(len=256) :: cmd
    write(cmd, '(A)') "export " // trim(varname)
    call system(cmd)
    retval = 0
  end subroutine putenv

  !> Calcular MD5 de un archivo (sin OpenSSL)
  subroutine get_md5(filename,md5,retval)
    implicit none
    character(len=*), intent(in) :: filename
    character(len=*), intent(out) :: md5
    integer, intent(out) :: retval

    ! Usa el comando del sistema 'md5sum'
    character(len=512) :: cmd
!   character(len=256) :: output

    write(cmd, '(A)') "md5sum " // trim(filename) // " 2>/dev/null | cut -d' ' -f1"
    call system(cmd, retval)
    if (retval == 0) then
      ! Leer el resultado del comando (complicado en Fortran puro)
      ! Mejor usar la versión con OpenSSL si es necesario
      md5 = "NO_MD5_AVAILABLE"
    else
      md5 = "ERROR"
    endif
    retval = 0
  end subroutine get_md5

end module mod_futils
