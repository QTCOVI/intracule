!
!     Computes the coefficients of the first NMAX physicist's
!     Hermite polynomials.
!
!     Chermite (i,j) gives the coefficient of x^j of Hermite polynomial
!     of degree i.
!
!     NUMERICAL NOTES (version 2, 2024):
!     --------------------------------
!     This routine is only used to fill the table Chermite(i,j). The
!     version 2 of hermicart does NOT use this table, so hermitec is
!     kept for backward compatibility and for any other routine that
!     might need the coefficients.
!
!     The table is filled with EXACT integer values (converted to
!     double precision). For n <= 20 the largest coefficient is
!     ~10^15, which fits well within double precision (max ~10^308).
!
!     A sanity check is performed at the end: if any coefficient is
!     not exactly representable as an integer in double precision, a
!     warning is issued. This is not expected for n <= 20.
!
!-----------------------------------------------------------------------
      subroutine hermitec (Chermite,nmax)
      implicit none
      integer(kind=4), intent(in)  :: nmax
      real   (kind=8), intent(out) :: Chermite(0:nmax,0:nmax)

      integer(kind=4) :: i, j
      real   (kind=8) :: maxc

      Chermite = 0.0D+00
!
!     Hermite polynomial of degree  0
!
        Chermite( 0, 0) =                  1.0D+00
!
!     Hermite polynomial of degree  1
!
        Chermite( 1, 1) =                  2.0D+00
!
!     Hermite polynomial of degree  2
!
        Chermite( 2, 0) =                 -2.0D+00
        Chermite( 2, 2) =                  4.0D+00
!
!     Hermite polynomial of degree  3
!
        Chermite( 3, 1) =                -12.0D+00
        Chermite( 3, 3) =                  8.0D+00
!
!     Hermite polynomial of degree  4
!
        Chermite( 4, 0) =                 12.0D+00
        Chermite( 4, 2) =                -48.0D+00
        Chermite( 4, 4) =                 16.0D+00
!
!     Hermite polynomial of degree  5
!
        Chermite( 5, 1) =                120.0D+00
        Chermite( 5, 3) =               -160.0D+00
        Chermite( 5, 5) =                 32.0D+00
!
!     Hermite polynomial of degree  6
!
        Chermite( 6, 0) =               -120.0D+00
        Chermite( 6, 2) =                720.0D+00
        Chermite( 6, 4) =               -480.0D+00
        Chermite( 6, 6) =                 64.0D+00
!
!     Hermite polynomial of degree  7
!
        Chermite( 7, 1) =              -1680.0D+00
        Chermite( 7, 3) =               3360.0D+00
        Chermite( 7, 5) =              -1344.0D+00
        Chermite( 7, 7) =                128.0D+00
!
!     Hermite polynomial of degree  8
!
        Chermite( 8, 0) =               1680.0D+00
        Chermite( 8, 2) =             -13440.0D+00
        Chermite( 8, 4) =              13440.0D+00
        Chermite( 8, 6) =              -3584.0D+00
        Chermite( 8, 8) =                256.0D+00
!
!     Hermite polynomial of degree  9
!
        Chermite( 9, 1) =              30240.0D+00
        Chermite( 9, 3) =             -80640.0D+00
        Chermite( 9, 5) =              48384.0D+00
        Chermite( 9, 7) =              -9216.0D+00
        Chermite( 9, 9) =                512.0D+00
!
!     Hermite polynomial of degree 10
!
        Chermite(10, 0) =             -30240.0D+00
        Chermite(10, 2) =             302400.0D+00
        Chermite(10, 4) =            -403200.0D+00
        Chermite(10, 6) =             161280.0D+00
        Chermite(10, 8) =             -23040.0D+00
        Chermite(10,10) =               1024.0D+00
!
!     Hermite polynomial of degree 11
!
        Chermite(11, 1) =            -665280.0D+00
        Chermite(11, 3) =            2217600.0D+00
        Chermite(11, 5) =           -1774080.0D+00
        Chermite(11, 7) =             506880.0D+00
        Chermite(11, 9) =             -56320.0D+00
        Chermite(11,11) =               2048.0D+00
!
!     Hermite polynomial of degree 12
!
        Chermite(12, 0) =             665280.0D+00
        Chermite(12, 2) =           -7983360.0D+00
        Chermite(12, 4) =           13305600.0D+00
        Chermite(12, 6) =           -7096320.0D+00
        Chermite(12, 8) =            1520640.0D+00
        Chermite(12,10) =            -135168.0D+00
        Chermite(12,12) =               4096.0D+00
!
!     Hermite polynomial of degree 13
!
        Chermite(13, 1) =           17297280.0D+00
        Chermite(13, 3) =          -69189120.0D+00
        Chermite(13, 5) =           69189120.0D+00
        Chermite(13, 7) =          -26357760.0D+00
        Chermite(13, 9) =            4392960.0D+00
        Chermite(13,11) =            -319488.0D+00
        Chermite(13,13) =               8192.0D+00
!
!     Hermite polynomial of degree 14
!
        Chermite(14, 0) =          -17297280.0D+00
        Chermite(14, 2) =          242161920.0D+00
        Chermite(14, 4) =         -484323840.0D+00
        Chermite(14, 6) =          322882560.0D+00
        Chermite(14, 8) =          -92252160.0D+00
        Chermite(14,10) =           12300288.0D+00
        Chermite(14,12) =            -745472.0D+00
        Chermite(14,14) =              16384.0D+00
!
!     Hermite polynomial of degree 15
!
        Chermite(15, 1) =         -518918400.0D+00
        Chermite(15, 3) =         2421619200.0D+00
        Chermite(15, 5) =        -2905943040.0D+00
        Chermite(15, 7) =         1383782400.0D+00
        Chermite(15, 9) =         -307507200.0D+00
        Chermite(15,11) =           33546240.0D+00
        Chermite(15,13) =           -1720320.0D+00
        Chermite(15,15) =              32768.0D+00
!
!     Hermite polynomial of degree 16
!
        Chermite(16, 0) =          518918400.0D+00
        Chermite(16, 2) =        -8302694400.0D+00
        Chermite(16, 4) =        19372953600.0D+00
        Chermite(16, 6) =       -15498362880.0D+00
        Chermite(16, 8) =         5535129600.0D+00
        Chermite(16,10) =         -984023040.0D+00
        Chermite(16,12) =           89456640.0D+00
        Chermite(16,14) =           -3932160.0D+00
        Chermite(16,16) =              65536.0D+00
!
!     Hermite polynomial of degree 17
!
        Chermite(17, 1) =        17643225088.0D+00
        Chermite(17, 3) =       -94097203200.0D+00
        Chermite(17, 5) =       131736084480.0D+00
        Chermite(17, 7) =       -75277762560.0D+00
        Chermite(17, 9) =        20910489600.0D+00
        Chermite(17,11) =        -3041525760.0D+00
        Chermite(17,13) =          233963520.0D+00
        Chermite(17,15) =           -8912896.0D+00
        Chermite(17,17) =             131072.0D+00
!
!     Hermite polynomial of degree 18
!
        Chermite(18, 0) =       -17643225088.0D+00
        Chermite(18, 2) =       317578051584.0D+00
        Chermite(18, 4) =      -846874804224.0D+00
        Chermite(18, 6) =       790416523264.0D+00
        Chermite(18, 8) =      -338749947904.0D+00
        Chermite(18,10) =        75277762560.0D+00
        Chermite(18,12) =        -9124577280.0D+00
        Chermite(18,14) =          601620480.0D+00
        Chermite(18,16) =          -20054016.0D+00
        Chermite(18,18) =             262144.0D+00
!
!     Hermite polynomial of degree 19
!
        Chermite(19, 1) =      -670442586112.0D+00
        Chermite(19, 3) =      4022655254528.0D+00
        Chermite(19, 5) =     -6436248616960.0D+00
        Chermite(19, 7) =      4290832760832.0D+00
        Chermite(19, 9) =     -1430277586944.0D+00
        Chermite(19,11) =       260050452480.0D+00
        Chermite(19,13) =       -26671841280.0D+00
        Chermite(19,15) =         1524105216.0D+00
        Chermite(19,17) =          -44826624.0D+00
        Chermite(19,19) =             524288.0D+00
!
!     Hermite polynomial of degree 20
!
        Chermite(20, 0) =       670442586112.0D+00
        Chermite(20, 2) =    -13408851197952.0D+00
        Chermite(20, 4) =     40226554642432.0D+00
        Chermite(20, 6) =    -42908325511168.0D+00
        Chermite(20, 8) =     21454162755584.0D+00
        Chermite(20,10) =     -5721110347776.0D+00
        Chermite(20,12) =       866834841600.0D+00
        Chermite(20,14) =       -76205260800.0D+00
        Chermite(20,16) =         3810263040.0D+00
        Chermite(20,18) =          -99614720.0D+00
        Chermite(20,20) =            1048576.0D+00
!
!-----Sanity check: all coefficients should be exactly representable
!     as integers in double precision for n <= 20. If any coefficient
!     exceeds 2^53 (~9.0e15), it is no longer exactly representable,
!     and the table cannot be trusted for exact arithmetic.
!
      maxc = 0.0d0
      do i = 0, nmax
        do j = 0, nmax
          if (abs(Chermite(i,j)) .gt. maxc) maxc = abs(Chermite(i,j))
        enddo
      enddo
      if (maxc .gt. 9.0d15) then
        write (*,'(a,i3,a,e12.4)') 
     &    '# hermitec: WARNING max |Chermite| = ',maxc, 
     &    ' for nmax = ',nmax
        write (*,'(a)') 
     &    '# hermitec: coefficients are no longer exactly representable'
        write (*,'(a)') 
     &    '# hermitec: in double precision. The table may be inaccurate'
      endif

!
      return
      end
