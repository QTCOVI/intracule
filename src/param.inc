!
!.....global parameters
!
      integer, parameter :: maxtype=35
      real(kind=8), parameter :: pi = 3.141592653589793238462643d0
      common       /datai/ nlm(35,3)
      integer, parameter :: mxfact=100
      integer, parameter :: mxind=200
      integer, parameter :: mcomb=20
      common  / const / fact(0:mxfact), facts(-1:mxfact), 
     &                  comb(0:mcomb,0:mcomb), facd(-mxfact:mxfact)
      logical icorr,rhf,rohf,onedet,wfnat,cciqa,mp2iqa,onedcou,dmftxc
      common  /logic / icorr,rhf,rohf,onedet,wfnat,cciqa,mp2iqa,onedcou,
     &                 dmftxc
      !
!.....Integer parameter that defines the type of Density Matrix
!     Functional Theory (DMFT) exchange-correlation that is used.
!
      integer idmft
      common     /typexc/ idmft
!
!.....MAXLPRIM = Maximun l value in a primitive cartesian GTO
!
      common     / maxlvalue / maxlprim
!
!-----The variable NRECWFN controls the number of records of the logical
!     unit IWFN in the 'rdwfn.f' routine before starting to read the
!     determinants of the calculation (in a CAS calculation) or the 
!     1-RDM and 2-RDM (in a MP2 or CCIQA calculation).
!
      common /readwfn/ nrecwfn
