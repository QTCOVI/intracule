!
!.....stdio.inc - standard input, output, and error units.
!
      integer, parameter :: stdin  = 5
      integer, parameter :: stdout = 6
      integer, parameter :: stderr = 0
      integer            uin, uout, udat, uerr, udatt
      common  /idat/     uin, uout, udat, uerr, udatt
