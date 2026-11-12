      subroutine initial()
c
c.....initialize some variables
c
      include  'implicit.inc'
      include  'param.inc'

      sqpi = sqrt(pi)
c
c.....nlm keeps the nlm values of x^n y^l z^m gaussian
c
      do j=1,3
        do i=1,35
          nlm(i,j)=0
        enddo
      enddo
c
c-----Definition of Cartesian Gaussian Functions
c
c.....p's
c
      nlm(2,1)=1      !px
      nlm(3,2)=1      !py
      nlm(4,3)=1      !pz
c
c.....d's
c
      nlm(5,1)=2      !xx
      nlm(6,2)=2      !yy
      nlm(7,3)=2      !zz
      nlm(8,1)=1      !xy
      nlm(8,2)=1
      nlm(9,1)=1      !xz
      nlm(9,3)=1
      nlm(10,2)=1     !yz
      nlm(10,3)=1
c
c.....f's
c
      nlm(11,1)=3     !xxx
      nlm(12,2)=3     !yyy
      nlm(13,3)=3     !zzz    
      nlm(14,1)=2     !xxy
      nlm(14,2)=1
      nlm(15,1)=2     !xxz
      nlm(15,3)=1
      nlm(16,2)=2     !yyz
      nlm(16,3)=1 
      nlm(17,1)=1     !xyy
      nlm(17,2)=2 
      nlm(18,1)=1     !xzz
      nlm(18,3)=2     
      nlm(19,2)=1     !yzz
      nlm(19,3)=2
      nlm(20,1)=1     !xyz
      nlm(20,2)=1 
      nlm(20,3)=1 
c
c.....g's
c
      nlm(23,3)=4     !zzzz
      nlm(29,2)=1     !yzzz
      nlm(29,3)=3
      nlm(32,2)=2     !yyzz 
      nlm(32,3)=2 
      nlm(27,2)=3     !yyyz
      nlm(27,3)=1 
      nlm(22,2)=4     !yyyy
      nlm(28,1)=1     !xzzz 
      nlm(28,3)=3 
      nlm(35,1)=1     !xyzz
      nlm(35,2)=1
      nlm(35,3)=2
      nlm(34,1)=1     !xyyz
      nlm(34,2)=2 
      nlm(34,3)=1
      nlm(26,1)=1     !xyyy
      nlm(26,2)=3 
      nlm(31,1)=2     !xxzz
      nlm(31,3)=2
      nlm(33,1)=2     !xxyz 
      nlm(33,2)=1
      nlm(33,3)=1
      nlm(30,1)=2     !xxyy
      nlm(30,2)=2 
      nlm(25,1)=3     !xxxz
      nlm(25,3)=1 
      nlm(24,1)=3     !xxxy
      nlm(24,2)=1
      nlm(21,1)=4     !xxxx
c
c     fact(i) = i! 
c     facts(i) = (i+1/2)! 
c     facd(i) = i!!
c----------------------------------------------------------------------
c
      fact(0) = 1.0d0
      facts(-1) = sqpi
      facd(-1) = 1d0
      facd(0) = 1d0
      facts(0) = facts(-1) * .5d0
      do i=1,mxfact
        fact(i) = fact(i-1) * dfloat(i)
        facts(i) = facts(i-1) * (2d0*dble(i)+1d0) * .5d0
        facd(i) = dble(i)*facd(i-2)
      end do
c
c.....negative double factorials
c
      mf21 = mxfact/2 - 1
      do i = 1, mf21
        facd(-i-i-1) = (-1)**i / facd(i+i-1)
      enddo
c----------------------------------------------------------------------
c               ind(i) = i+(i-1)+(i-2)+...0 
c----------------------------------------------------------------------
!     ind(0) = 0
!     do i = 1, mxind
!        ind(i) = ind(i-1) + i
!     end do
c----------------------------------------------------------------------
c            combinatorial numbers com(i,j)=(i over j)
c----------------------------------------------------------------------
      do i=0, mcomb
        do j=0,i/2
          ij=i-j
          comb(i,j)=fact(i)/fact(j)/fact(ij)
          if (j.ne.ij) comb(i,ij)=comb(i,j)
        enddo
      enddo
      return
      end 
