real*8 function funbess(r) result(val)
   use acumula_intracule
   implicit none
   real(kind=8), intent(in) :: r
   real(kind=8) :: r1, r2, r3, r4
   logical :: verx, usemx
   common /quadarea/ lwfun, verbose, usemuller, first_r_value, wfnaux
   integer lwfun
   logical verbose
   logical usemuller
   logical first_r_value
   character(len=132) :: wfnaux
   character(len=132) :: wfnfile

   wfnfile = trim(wfnaux)
   verx = verbose
   usemx = usemuller
   call quadpac_bessic (wfnfile,r,r1,r2,r3,r4)
   val = r1
!
!--Temporal mofification for using Vee as objective function in the integration
!
!  val = r1/r
   call acumula_push (r,r1,r2,r3,r4)
   first_r_value = .false.
   write (lwfun,1) r,r1,r2,r3,r4
1  format (1x,F15.8,8(2x,E22.15))
end function funbess
