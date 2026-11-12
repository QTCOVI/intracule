c
c-----------------------------------------------------------------------
c
      logical function setdble (dblevar, line, lp)
c
c.....setdble - get a real number from line and sets dblevar to it.
c
c     If a valid real number is not found, setdble returns .false.
c
      real(kind=8)      dblevar
      character*(*)     line
      integer           lp
c
      logical           isdigit, hasdigits
      integer           dsign, dexpo, desig
      real(kind=8)      dval, dpow, ten
      parameter (ten=1.0d+01)
c
      integer           istate, none, sign, digit, dot, decim, expo
      integer           esign, edigit, done
      parameter (none=0, sign=1, digit=2, dot=3, decim=4, expo=5,
     &           esign=6, edigit=7, done=8)
      character*(1)     ch, blank, newline, null
      isdigit(ch) = ch.ge.'0' .and. ch.le.'9'
      blank = char(32)
      newline = char(10)
      null = char(0)
c
c.....we don't know for sure ...
c
      setdble = .false.
c
c.....test for end of string
c
      ll = len (line)
      if (lp.gt.ll) return
c
c.....skip blanks
c
      istate = none
      do while (lp.le.ll .and. istate.eq.none)
         if (line(lp:lp).ne.blank) istate = done
         lp = lp + 1
      enddo
      if (istate.eq.none) return
      lp = lp - 1
c
c.....finite authomaton
c
      i = lp
      istate = none
      hasdigits = .false.
      dsign = 1
      dval  = 0d0
      dpow  = 1d0
      dexpo = 0
      desig = 1
      do while (istate.ne.done .and. i.le.ll)
         if (line(i:i).eq.'+' .or.
     &       line(i:i).eq.'-') then
            if (istate.eq.none) then
               istate = sign
               if (line(i:i).eq.'-') dsign = -1
            else if (istate.eq.expo) then
               istate = esign
               if (line(i:i).eq.'-') desig = -1
            else
               return
            endif
         else if (line(i:i).eq.'.') then
            if (istate.eq.none .or.
     &          istate.eq.sign .or.
     &          istate.eq.digit) then
               istate = dot
            else
               return
            endif
         else if (index('eEdD',line(i:i)).gt.0) then
            if (istate.eq.digit .or.
     &          istate.eq.decim .or.
     &         (istate.eq.dot .and. hasdigits)) then
               istate = expo
            else
               return
            endif
         else if (isdigit(line(i:i))) then
            if (istate.eq.none .or.
     &          istate.eq.sign .or.
     &          istate.eq.digit) then
               istate = digit
               hasdigits = .true.
               dval = ten*dval + ichar(line(i:i)) - ichar('0')
            else if (istate.eq.dot .or.
     &               istate.eq.decim) then
               istate = decim
               dval = ten*dval + ichar(line(i:i)) - ichar('0')
               dpow = ten*dpow
            else if (istate.eq.expo  .or.
     &               istate.eq.esign .or.
     &               istate.eq.edigit) then
               istate = edigit
               dexpo = int(ten)*dexpo + ichar(line(i:i)) - ichar('0')
            else
               return
            endif
         else if (index(blank//','//null//newline,line(i:i)).gt.0) then
            if (istate.eq.digit  .or.
     &          istate.eq.decim  .or.
     &          istate.eq.edigit .or.
     &         (istate.eq.dot .and. hasdigits)) then
               istate = done
            else
               return
            endif
         else
            return
         endif
         i = i + 1
      enddo
c
c.....check for done when the end is reached
c
      if (istate.ne.done .and. i.eq.ll+1) then
         if (istate.eq.digit  .or.
     &       istate.eq.decim  .or.
     &       istate.eq.edigit .or.
     &      (istate.eq.dot .and. hasdigits)) then
            istate = done
         else
            return
         endif
      else if (istate.eq.done) then
         i = i - 1
      endif
c
c.....set the variable if done
c
      if (istate.eq.done) then
         dblevar = (dsign*dval/dpow) * ten**(desig*dexpo)
         lp = i
         setdble = .true.
      endif
c
c.....end
c
      return
      end
c
c-----------------------------------------------------------------------
c
c-----Analizes a line of text with the structure 
c
c     'var1=value  var2=value .... ' 
c
c     and returns in the array var() the values of var1, var2, ..., i.e.
c
c     var(1) = value1
c     var(2) = value2
c            :
c
c     An arbitrary number of blank characters can be found before and 
c     after the '=' symbol of the analized line.
c
c-----------------------------------------------------------------------
c
      subroutine anline (line,var,n,maxdata) 
c
      implicit none
      character(len=*)  line
      real   (kind=8) var(maxdata)
      real   (kind=8) value
      integer (kind=4) maxdata
      integer (kind=4) n,lp,ll,leng
!     character(len=1)  ch
      logical lineleft,setdble

      var(:) = 0.0d+00
      ll = leng(line)
      lineleft = .true.
      n = 0
      lp = 1
      do while (lineleft)
        do while (line(lp:lp).ne.'=')
          lp = lp + 1
        enddo
        lp = lp + 1
        if (setdble(value,line,lp)) then
          n = n + 1
          var(n) = value
        else
          return
        endif
        if (n .eq. maxdata .or. lp .ge. ll) then
          lineleft = .false.
          if (n.eq.maxdata) then
            write (0,2) '# anline.f: There are probably un-read data'
          endif
        endif
      enddo
      return
!1    format (' # analine.f: Current line has not the expected format')
 2    format (' # anline.f: There are probably un-read data')
      end
