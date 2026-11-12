      integer (kind=4) function ngtox ()
      USE space_for_wfnbasis
      implicit real*8 (a-h,o-z)
      include 'param.inc'
      include 'wfn.inc'
      integer(kind=4) lp1lp2

      lp1lp2 = 0
      do i=1,nprims
       lp1lp2 = max (lp1lp2,sum(nlm(ityp(i),1:3)))
      enddo
      lp1lp2 = (lp1lp2+1)*(lp1lp2+2)/2
      ngtox = lp1lp2
      return
      end
