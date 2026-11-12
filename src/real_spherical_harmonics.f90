!===============================================================================
! real_spherical_harmonics.f90
!
! Calcula los ARMONICOS ESFERICOS REALES S_{l,m}(theta,phi) para
! l = 0,...,l_max, m = -l,...,l, y angulos theta, phi arbitrarios,
! almacenandolos en un array de UN SOLO INDICE.
!
! ---------------------------------------------------------------------------
! IMPLEMENTACION 100% FORTRAN (sin dependencia de GSL)
! ---------------------------------------------------------------------------
! La version anterior de este fichero llamaba a gsl_sf_legendre_array_e.
! Esa funcion produjo un error de corrupcion de memoria (malloc corrupted)
! en el entorno del usuario, muy probablemente por una incompatibilidad
! de la interfaz C (size_t/enum) con la version de GSL instalada, que no
! ha podido reproducirse ni depurarse por no disponer aqui de gfortran/GSL.
!
! Para eliminar por completo ese riesgo, esta version implementa
! DIRECTAMENTE EN FORTRAN la misma recursion normalizada y estable que usa
! GSL internamente (sin ningun binding C), de modo que:
!   - No hay riesgo de incompatibilidad de ABI con ninguna libreria externa.
!   - Se mantiene la fiabilidad numerica para l,m grandes (no se calculan
!     nunca factoriales ni P_l^m sin normalizar).
!
! ---------------------------------------------------------------------------
! FIABILIDAD NUMERICA PARA l, m GRANDES
! ---------------------------------------------------------------------------
! Se trabaja siempre con la funcion YA normalizada
!
!     N_l^m(x) = sqrt[ (2l+1)/(4 pi) * (l-m)!/(l+m)! ] * P_l^m(x)
!
! usando una recursion de tres terminos en la que los factoriales se han
! cancelado algebraicamente de antemano, de modo que solo aparecen
! cocientes de numeros de tamano moderado (nunca (l+m)! o (l-m)! por
! separado). Esto evita el overflow/underflow que se produciria si se
! calculase P_l^m sin normalizar (que puede alcanzar ~10^300) y el
! prefactor factorial por separado.
!
! Formulas usadas (deduccion estandar; ver p.ej. Holmes & Featherstone 2002):
!
!   N_0^0 = 1/sqrt(4 pi)
!
!   N_m^m(x)   = -sqrt[(2m+1)/(2m)] * sin(theta) * N_{m-1}^{m-1}(x)     (m>=1)
!   N_{m+1}^m(x) = sqrt(2m+3) * x * N_m^m(x)
!   N_l^m(x)   = a_l^m * x * N_{l-1}^m(x) - b_l^m * N_{l-2}^m(x)        (l>=m+2)
!
!   a_l^m = sqrt[ (2l-1)(2l+1) / ((l-m)(l+m)) ]
!   b_l^m = sqrt[ (2l+1)(l+m-1)(l-m-1) / ((2l-3)(l-m)(l+m)) ]
!
! (el signo "-" en la recursion sectorial es la fase de Condon-Shortley
!  (-1)^m; es opcional, ver el argumento condon_shortley mas abajo)
!
! ---------------------------------------------------------------------------
! DEFINICION DE LOS ARMONICOS ESFERICOS REALES
! ---------------------------------------------------------------------------
!     S_{l,0}  =              N_l^0(cos theta)
!     S_{l,m}  = sqrt(2) * N_l^m(cos theta) * cos(m*phi)   , m > 0
!     S_{l,-m} = sqrt(2) * N_l^m(cos theta) * sin(m*phi)   , m > 0
!
! normalizados a integral_(4pi) |S_lm|^2 dOmega = 1.
!
! ---------------------------------------------------------------------------
! ALMACENAMIENTO EN UN UNICO INDICE
! ---------------------------------------------------------------------------
!     idx_ylm(l,m) = l*(l+1) + m ,   m = -l,...,l
!
! El bloque de un l dado ocupa los indices l^2,...,(l+1)^2-1; el array
! completo tiene tamano (lmax+1)^2 (indices 0..(lmax+1)^2-1).
!
! NOTA: si l_max es muy grande (por encima de ~30000-40000), cambia los
! indices enteros de "integer" (32 bits) a "integer(kind=8)", ya que el
! tamano del array y algunos indices crecen como l_max^2.
!
! Compilacion (ya NO requiere GSL para este fichero):
!   gfortran -o real_spherical_harmonics real_spherical_harmonics.f90
!===============================================================================

module real_spherical_harmonics_mod
  use, intrinsic :: iso_fortran_env, only: dp => real64
  implicit none

  real(dp), parameter :: PIVAL = 3.14159265358979323846_dp

contains

  pure function idx_ylm(l, m) result(idx)
    ! Indice de un solo entero para S_{l,m}: bloque de l en [l^2,(l+1)^2-1]
    integer, intent(in) :: l, m
    integer :: idx
    idx = l*(l+1) + m
  end function idx_ylm

  subroutine compute_real_spherical_harmonics (lmax, theta, phi, S, condon_shortley)
    !-------------------------------------------------------------------------
    !   lmax     (in)  : l maximo (>=0)
    !   theta,phi(in)  : angulos en RADIANES (theta en [0,pi], phi libre)
    !   S        (out) : array 0:(lmax+1)**2-1 con los armonicos reales,
    !                    S(idx_ylm(l,m)) = S_{l,m}
    !   condon_shortley (in, opcional) : .true. para incluir la fase
    !             (-1)^m, .false. (por defecto) para excluirla. Solo
    !             cambia el signo de los terminos con m impar.
    !-------------------------------------------------------------------------
    integer,  intent(in)  :: lmax
    real(dp), intent(in)  :: theta, phi
    real(dp), intent(out) :: S(0:(lmax+1)**2-1)
    logical,  intent(in), optional :: condon_shortley

    real(dp), allocatable :: Nbar(:,:)   ! Nbar(l,m), 0<=m<=l<=lmax
    real(dp) :: x, sinth, sqrt2, cs_sign, a, b, cmphi, smphi
    integer :: l, m
    logical :: use_cs

    use_cs = .false.
    if (present(condon_shortley)) use_cs = condon_shortley
    cs_sign = 1.0_dp
    if (use_cs) cs_sign = -1.0_dp

    x     = cos(theta)
    sinth = sin(theta)
    sqrt2 = sqrt(2.0_dp)

    allocate(Nbar(0:lmax,0:lmax))
    Nbar = 0.0_dp

    ! --- Valores sectoriales N_m^m ---
    Nbar(0,0) = 1.0_dp/sqrt(4.0_dp*PIVAL)
    do m = 1, lmax
       Nbar(m,m) = cs_sign * sqrt(real(2*m+1,dp)/real(2*m,dp)) * sinth * Nbar(m-1,m-1)
    end do

    ! --- Recursion en l, para cada m fijo ---
    do m = 0, lmax
       if (m+1 <= lmax) Nbar(m+1,m) = sqrt(real(2*m+3,dp)) * x * Nbar(m,m)

       do l = m+2, lmax
          a = sqrt( (real(2*l-1,dp)*real(2*l+1,dp)) / (real(l-m,dp)*real(l+m,dp)) )
          b = sqrt( (real(2*l+1,dp)*real(l+m-1,dp)*real(l-m-1,dp)) &
                    / (real(2*l-3,dp)*real(l-m,dp)*real(l+m,dp)) )
          Nbar(l,m) = a*x*Nbar(l-1,m) - b*Nbar(l-2,m)
       end do
    end do

    ! --- Construccion de los armonicos esfericos reales ---
    do l = 0, lmax
       S(idx_ylm(l,0)) = Nbar(l,0)
       do m = 1, l
          cmphi = cos(real(m,dp)*phi)
          smphi = sin(real(m,dp)*phi)
          S(idx_ylm(l, m)) = sqrt2 * Nbar(l,m) * cmphi
          S(idx_ylm(l,-m)) = sqrt2 * Nbar(l,m) * smphi
       end do
    end do

    deallocate(Nbar)

  end subroutine compute_real_spherical_harmonics

end module real_spherical_harmonics_mod


!program main
!  use, intrinsic :: iso_fortran_env, only: dp => real64
!  use real_spherical_harmonics_mod
!  implicit none
!
!  integer :: lmax, l, m, lm
!  real(dp) :: x,y,z,theta,phi,r
!  real(dp), allocatable :: rsh(:)
!
!  print *, "Introduzca l_max (entero >= 0): "
!  read(*,*) lmax
!  if (lmax < 0) then
!     print *, "Error: l_max debe ser >= 0."
!     stop 1
!  end if
!
!  print *, "Introduzca x,y,z"
!  read (*,*) x,y,z
!
!  r = sqrt(x*x + y*y + z*z)
!  if (r > 0.0_dp) then
!       ! Se acota z/rr a [-1,1] para evitar que un redondeo de punto
!       ! flotante (p.ej. 1.0000000000000002) haga que acos reciba un
!       ! argumento fuera de su dominio.
!       theta = acos( max(-1.0_dp, min(1.0_dp, z/r)) )
!  else
!     theta = 0.0_dp   ! convenio para el origen, donde theta no esta definido
!  end if
!
!  phi = atan2(y, x)          ! devuelve un valor en (-pi, pi]
!  if (phi < 0.0_dp) phi = phi + 2.0_dp*PI   ! se desplaza a [0, 2*pi)
!
!  allocate (rsh(0:(lmax+1)**2-1))
!
!! condon_shortley=.false. (valor por defecto): sin la fase (-1)^m
!  call compute_real_spherical_harmonics (lmax,theta,phi,rsh,condon_shortley=.false.)
!
!  write (6,1) 
!! print *
!! print '(A)', "  l     m         (x,y,z) ==> (theta,phi)/rad S_lm(theta,phi)              indice"
!! print '(A)', "---------------------------------------------------------------------------------"
!  lm = 0 
!  do l = 0, lmax
!     do m = -l, l
!!       write (*, '(I4,2X,I4,6X,E20.12,8X,I8)') l,m,x,y,z,theta,phi,rsh(lm),idx_ylm(l,m)
!        write (6,2) l,m,x,y,z,theta,phi,rsh(lm),l*(l+1)+m
!        lm = lm + 1
!     end do
!  end do
!1 format ('   l     m',17x,'x',15x,'y',15x,'z',12x,'theta',12x,'phi',11x,'S_lm(theta,phi)'14x,'Indice'/,138('-'))
!2 format (2(I4,2X),4X,5(1x,F15.8),3x,E20.12,8X,I8)
!
!  deallocate (rsh)
!
!end program main
