! Copyright (c) 2018 
! Jose Luis Casals Sainz <jluiscasalssainz@gmail.com>,
! Ángel Martín Pendás <angel@fluor.quimica.uniovi.es> and 
! Evelio Francisco Miguelez <evelio@uniovi.es>. 
!
! integ is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or (at
! your option) any later version.
! 
! integ is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
module mod_prec

  implicit none
  save
  private

  public :: dp, dp_size
  public :: ip, ip_size
  public :: size_t, size_t_size
  public :: smallest
  public :: print_kind_info

  integer, parameter :: real_4 = selected_real_kind(6, 30)
  integer, parameter :: real_8 = selected_real_kind(14, 200)
  integer, parameter :: real_4_size = 4
  integer, parameter :: real_8_size = 8

  integer, parameter :: int_size = bit_size(0)/8

  integer, parameter :: int_1  = selected_int_kind(2)
  integer, parameter :: int_1_size = bit_size(int(0, int_1))/8

  integer, parameter :: int_2  = selected_int_kind(4)
  integer, parameter :: int_2_size = bit_size(int(0, int_2))/8

  integer, parameter :: int_4 = selected_int_kind(5)
  integer, parameter :: int_4_size = bit_size(int(0, int_4))/8

  integer, parameter :: int_8 = selected_int_kind(10)
  integer, parameter :: int_8_size = bit_size(int(0, int_8))/8

  integer, parameter :: ip = int_4
  integer, parameter :: ip_size = int_4_size
  integer, parameter :: sp = real_4
  integer, parameter :: sp_size = real_4_size
  integer, parameter :: dp = real_8
  integer, parameter :: dp_size = real_8_size
  integer, parameter :: size_t = int_8
  integer, parameter :: size_t_size = int_8_size
  
  real(kind=dp), parameter :: smallest = epsilon(0.0_dp)

contains

  subroutine print_kind_info(uout)

    implicit none
    integer(kind=ip), intent(in) :: uout
 
    write (uout,'(t2,a)') '# '
    write (uout,'(t2,a)') '# Data type information '
    write (uout,'(t2,a,/,t2,a,2(/,t2,a,i2),3(/,t2,a,E15.8))') &
      '# REAL ', &
      '#  Data type name : dp', &
      '#  Kind value : ', kind(0.0_dp), &
      '#  Precision : ', precision(0.0_dp), &
      '#  Smallest nonnegligible quantity relative to 1 : ', epsilon(0.0_dp), &
      '#  Smallest positive number : ', tiny(0.0_dp), &
      '#  Largest representable number : ', huge(0.0_dp)
    write (uout,'(t2,a,/,t2,a,2(/,t2,a,i2),t2,a,i20)') &
      '# INTEGER ', &
      '#  Data type name : ip', &
      '#  Kind value : ', kind(0_ip), &
      '#  Bit size : ', bit_size(0_ip), &
      '#  Largest representable number : ', huge(0_ip)
    write (uout,'(t2,a,/,t2,a,2(/,t2,a,i2),t2,a,i20)') &
      '# INTEGER ', &
      '#  Data type name : size_t', &
      '#  Kind value : ', kind(0_size_t), &
      '#  Bit size : ', bit_size(0_size_t), &
      '#  Largest representable number : ', huge(0_size_t)
    write (uout,'(t2,a,/,t2,a,/,t2,a,i2)') &
      '# LOGICAL ', &
      '#  Data type name : default', &
      '#  Kind value : ', kind(.true.)
    write (uout,'(t2,a,/,t2,a,/,t2,a,i2)') &
      '# CHARACTER ', &
      '#  Data type name : default', &
      '#  Kind value : ', kind('c')
 
  end subroutine print_kind_info

end module mod_prec
