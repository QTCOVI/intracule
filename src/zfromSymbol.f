      real*8 function ZfromSymbol (symbol)
      include 'datatm.inc'
      integer(kind=4) i
      parameter (mline=132)
      character*(*) symbol
      character*(mline) uppcase

      do i=1,atminstall
        if (uppcase(symbol(1:2)).eq.uppcase(namatm(i)(1:2))) then
          ZfromSymbol=dble(i)
          return
        endif
      enddo
      ZfromSymbol = 0.0d+00
      return
      call error ('ZfromSymbol.f','Unknown atomic symbol in WFN',-1)
      end
