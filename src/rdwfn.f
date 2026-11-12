c
c-----------------------------------------------------------------------
c
      subroutine rdwfn (wfnfile)
c
c.....read wavefunction with slightly modified AIMPAC routine.
c
      USE        space_for_wfnbasis
      USE        space_for_primgto
      USE        space_for_wfncoef
      USE        space_for_corr
      USE        space_for_rdm1
      USE        space_for_rdm2
      include   'implicit.inc'
      include   'param.inc'
      include   'io.inc'
      include   'wfn.inc'
      include   'corr.inc'
      include   'primgto.inc'

      parameter (mline=132)
      real(kind=8),    allocatable,dimension (:)      :: oexpa
      real(kind=8),    allocatable,dimension (:,:)    :: coefa,coefat
      character(len=mline), allocatable,dimension (:)  :: eorbocc
      integer,   allocatable,dimension (:)      :: icena,itypa
      integer(kind=4) :: iord(ngtoG),ordtyp(ngtoG)
      real   (kind=8) :: cord(ngtoG)
      integer    leng
      logical    nonum
      logical    ok
      character*(mline) line,symbol,cutsym
      logical setdble,setword
      character(len=1)  ch
c
c-----Arrays needed in case of a WFN file based on a CCSD calculation
c
      real(kind=8) , allocatable,dimension(:,:) :: v1mata
      real(kind=8) , allocatable,dimension(:)   :: d1mata
c
      integer        ntfu(0:4)
      character * 1  iotipe(0:4)
      character * 1  blank
      character*(*)  wfnfile
      data  iotipe(0) /'S'/, iotipe(1) /'P'/, iotipe(2) /'D'/, 
     &      iotipe(3) /'F'/, iotipe(4) /'G'/
      data  blank /' '/
c
c.....local variables
c
      CHARACTER*80 wfnttl
      CHARACTER*4  mode
      character*17 label
      character*8  check
c
      iwfn      =  udat
      icorr     = .false.
      rhf       = .false.
      rohf      = .false.
      uhf       = .false.
      wfnat     = .false.
      onedcou   = .false.
      cciqa     = .false.  ! By default WFN is not a Coupled Cluster one
      mp2iqa    = .false.  ! By default WFN is not a MP2 one
      ndou      =  0
      nsin      =  0
      nalpha    =  0
      nbeta     =  0
      occmin    =  2.0d+00
      occmax    = -2.0d+00
      totel     =  0.0d+00
c
c-----Centers with a pseudo-potencial or effective core potential
c
      ncentcore = 0
c
      read (iwfn,101) wfnttl
      read (iwfn,102) mode,nmo,nprims,ncent
      nprimsorig = nprims
      nmotot = nmo+nmo
      npor   = nprims
c
c.....Allocate memory for the wfn coefficients.
c
      call allocate_space_for_wfncoef (nprims,nmo)
c
c.....Allocate memory for the wfn basis, geometry, ...
c
      call allocate_space_for_wfnbasis (nprims,nmo,ncent)
c
c.....Allocate memory for primitive gtos related arrays
c
      call allocate_space_for_primgto (nprims,ncent,mgrp,ngtoG)
c
c.....Allocate ialpha(), ibeta(), pind(), ipak(), and cidet().
c
      call allocate_space_for_corr (nmo)
c
c-----Allocate space for 1-RDM.
c
      call allocate_space_for_rdm1 (nmo)
c
c-----Allocate space for lines that contain (except in CCIQA and MP2
c     calculations) the orbital energies and occupations.
c
      allocate (eorbocc(nmo))

      ghost = .false.
      do i = 1,ncent
        read (iwfn,'(a)') line
        read (line,'(a8)') atnam(i)
        line = line(24:)
        lp = 1
        ok = .true.
        ok = setdble(xyz(i,1),line,lp) .and. ok
        ok = setdble(xyz(i,2),line,lp) .and. ok
        ok = setdble(xyz(i,3),line,lp) .and. ok
        if (.not.ok) then
          backspace (iwfn)
          read (iwfn,'(a)') line
          line = line(25:)
          read (line,'(3F12.8,10x,F5.1)') (xyz(i,k),k=1,3),charge(i)
        else
          ch=line(lp:lp)
          nonum = ch.eq.'+'.or.ch.eq.'-'.or.ch.eq.'.'.or.
     &           (ch.ge.'0'.and.ch.le.'9')
          do while ((.not.nonum).and.lp.lt.leng(line))
            lp = lp+1
            ch = line(lp:lp)
            nonum = ch.eq.'+' .or. ch.eq.'-' .or. ch.eq.'.'
     &        .or. (ch.ge.'0'.and.ch.le.'9')
          enddo
          ok = setdble(charge(i),line,lp)
          if (.not.ok) write (0,182) i+2
          if (.not.ok) stop
        endif

        if (abs(charge(i)).lt.1.0D-6) ghost(i)=.true.
        charget(i)=charge(i)
        if (.not.ghost(i)) then
          backspace (iwfn)
          read (iwfn,'(a)') line
          lp = 1
          ok = setword(cutsym,line,lp)
          symbol=''
          k=0
          do l=1,leng(cutsym)
            if (cutsym(l:l).ge.'0'.and.cutsym(l:l).le.'9') cycle
            k=k+1
            symbol(k:k)=cutsym(l:l)
          enddo
          if (leng(symbol).gt.2) then
            call error ('rdwfn.f','WFN: Non-valid atomic symbol',+1)
          elseif (leng(symbol).eq.1) then
            symbol(2:2)=blank
          endif
          zw = ZfromSymbol(symbol)
          diffz = zw - charge(i)
          if (abs(diffz).gt.1.0D-2) then
            ncentcore=ncentcore+1
!           thereiscore(ncentcore) = i
            zcore(ncentcore) = diffz
            charget(i) = charge(i) + diffz
          endif
        endif
      enddo
c
      read (iwfn,104) (icen(i),i=1,nprims)
      read (iwfn,104) (ityp(i),i=1,nprims)
      read (iwfn,105) (oexp(i),i=1,nprims)

      if (.not.allocated(oexpa)) then
        allocate (oexpa(nprims),stat=ier)
        if (ier.ne.0) stop 'rdwfn.f: Cannot allocate oexpa()'
      endif
      if (.not.allocated(icena)) then
        allocate (icena(nprims),stat=ier)
        if (ier.ne.0) stop 'rdwfn.f: Cannot allocate icena()'
      endif
      if (.not.allocated(itypa)) then
        allocate (itypa(nprims),stat=ier)
        if (ier.ne.0) stop 'rdwfn.f: Cannot allocate itypa()'
      endif

      do i=1,nprims
        if (ityp(i).gt.maxtype) then
          stop '# rdwfn.f: promolden cannot work with L > 4 primitives'
        endif
      enddo
      do i=1,nmo
        read (iwfn,'(a)') eorbocc(i)
        eorbocc(i) = trim (eorbocc(i))
        read (iwfn,107) (coef(i,j),j=1,nprims)
        coef(i+nmo,1:nprims)=coef(i,1:nprims)
c
c-------Extract occ(i) and eorb(i) from eorbocc(i)
c
        occ(i) = 0.0d+00
        eorb(i) = 0.0d+00
        lp = 1
        do while (eorbocc(i)(lp:lp).ne.'=' .and. lp.lt.mline)
          lp = lp + 1
        enddo
        if (lp.ge.mline) then
          occ(i) = 0.0d+00
          eorb(i) = 0.0d+00
        else
          lp = lp + 1
          if (.not.setdble(occ(i),eorbocc(i),lp)) then
            write (uerr,873) 'occ()'
          endif
          do while (eorbocc(i)(lp:lp).ne.'=' .and. lp.lt.mline)
            lp = lp + 1
          enddo
          if (lp.ge.mline) then
            eorb(i) = 0.0d+00
          else
            lp = lp + 1
            if (.not.setdble(eorb(i),eorbocc(i),lp)) then
              write (uerr,873) 'eorb()'
            endif
          endif
        endif
        totel = totel + occ(i)
      enddo
      read (iwfn,108) check
      if (check .ne. 'END DATA') stop '# rdwfn.f: End card not found 1'
C
C.....READ IN TOTAL SCF ENERGY AND -V/T
C
      read (iwfn,'(a)') line
c
c-----recover total energy and virial
c
      read (line,'(a10)') label
      line=line(11:)
      lp = 1
      do while (line(lp:lp).ne.'=' .and. lp.lt.mline)
        lp = lp + 1
      enddo
      if (lp.ge.mline) then
        totalenergy = 0.0d+00
        gamma = 0.0d+00
      else
        lp = lp + 1
        if (.not.setdble(totalenergy,line,lp)) then
          write (uerr,873) 'totalenergy'
        endif
        do while (line(lp:lp).ne.'=' .and. lp.lt.mline)
          lp = lp + 1
        enddo
        if (lp.ge.mline) then
          gamma = 0.0d+00
        else
          lp = lp + 1
          if (.not.setdble(gamma,line,lp)) then
            write (uerr,873) 'virial'
          endif
        endif
      endif

      nrecwfn = 0
      nrrec = 0
      if (label(1:5).eq.'MCSCF' .or.
     &    label(1:5).eq.'ALDET' .or.
     &    label(1:3).eq.'CIS' .or.
     &    label(1:5).eq.'ORMAS' .or.
     &    label(1:5).eq.'GENCI') then
c
c-------Determine nalpha, nbeta, ...
c
        do i=1,nmo
          if (occ(i).gt.occmax) occmax=occ(i)
          if (occ(i).lt.occmin) occmin=occ(i)
          if (occ(i).lt.0.0d+00) then
            if (abs(abs(occ(i))-1.0d+00).lt.1.0d-6) then
              occ(i) = -occ(i)
              nbeta = nbeta+1
              ibeta(nbeta) = i
              nsin = nsin+1
            else
              wfnat = .true.
            endif
          else
            if (abs(occ(i)-1.0d+00).lt.1.0d-6) then
              nalpha = nalpha+1
              ialpha(nalpha) = i
              nsin = nsin+1
            else if (abs(occ(i)-2.0d+00).lt.1.0d-6) then
              nalpha = nalpha+1
              nbeta = nbeta+1
              ialpha(nalpha) = i
              ibeta(nbeta) = i
              ndou = ndou+1
            else
              wfnat = .true.
            endif
          endif
        enddo
c
c.......correlated case. Equal number of NMO
c
        read (iwfn,*,end=9999) label
        icorr = .true.
        do i=nmo+1,nmo+nmo
          read (iwfn,*) label
          read (iwfn,107) (coef(i,j),j=1,nprims)
        enddo
        read (iwfn,108) check
        if (check .ne. 'END DATA') stop 'rdwfn.f: End card not found 2'
        read (iwfn,*) label
        read (iwfn,*) nelact, ndets, ncore, nact
        nel = nelact+2*ncore
c
c-------Set wfnat = .false. since, not only natural MOs, but also
c       canonical MOs have read in.
c
        wfnat = .false.
c
c-------Recompute the number of ALPHA and BETA MOs.
c
        nalpha = 0
        nbeta = 0
        read (iwfn,'(a)') line
        read (iwfn,*) (cidet(k),k=0,nelact)
        if (ncore.gt.0) then
          nalpha = ncore
          nbeta = ncore
          do i=1,ncore
            ialpha(i) = i
            ibeta(i) = i
          enddo
        endif
        do k=1,nelact
          mm=int(cidet(k))
          if (mm.gt.0) then
            nalpha = nalpha+1
            ialpha (nalpha) = ncore+int(cidet(k))
          else
            nbeta = nbeta+1
            ibeta (nbeta) = ncore-int(cidet(k))
          endif
        enddo
        totel = dble (nel)
c
c-------Number of records plus 1 read from iwfn to this point.
c
        nrecwfn = ncent+2*(nprims/20)+(2*nmo+1)*(nprims/5)+2*nmo+9+nrrec
        if (mod(nprims,5 ).ne.0) nrecwfn = nrecwfn + 2*nmo + 1
        if (mod(nprims,20).ne.0) nrecwfn = nrecwfn + 2
        close (iwfn)
      elseif (label(1:5).eq.'CCIQA'. or. label(1:3).eq.'MP2') then
c
c-------Number of records plus 1 read from iwfn to this point.
c
        nrecwfn = ncent+2*(nprims/20)+(nmo+1)*(nprims/5)+nmo+5+nrrec
        if (mod(nprims,20).ne.0) nrecwfn = nrecwfn + 2
        if (mod(nprims,5 ).ne.0) nrecwfn = nrecwfn + nmo + 1
        if (label(1:5).eq.'CCIQA') cciqa=.true.
        if (label(1:3).eq.'MP2')  mp2iqa=.true.
        ndets=1
        cidet(0)=1.0d+00
c
c.......Test that the number of electrons is integer
c
        nel = int(anint(totel))
        if (abs(totel-nel).gt.1.0d-6) write (uerr,331) totel
        nelact = nel
        ncore  = 0
        nact   = nmo
        close (iwfn)
      else
c
c-------Determine nalpha, nbeta, ...
c
        totel = 0.0d+00
        do i=1,nmo
          totel = totel + occ(i)
          if (occ(i).gt.occmax) occmax=occ(i)
          if (occ(i).lt.occmin) occmin=occ(i)
          if (occ(i).lt.0.0d+00) then
            if (abs(abs(occ(i))-1.0d+00).lt.1.0d-6) then
              occ(i) = -occ(i)
              nbeta = nbeta+1
              ibeta(nbeta) = i
              nsin = nsin+1
            else
              wfnat = .true.
            endif
          else
            if (abs(occ(i)-1.0d+00).lt.1.0d-6) then
              nalpha = nalpha+1
              ialpha(nalpha) = i
              nsin = nsin+1
            else if (abs(occ(i)-2.0d+00).lt.1.0d-6) then
              nalpha = nalpha+1
              nbeta = nbeta+1
              ialpha(nalpha) = i
              ibeta(nbeta) = i
              ndou = ndou+1
            else
              wfnat = .true.
            endif
          endif
        enddo
        if (abs(occmin-2.0d+00).lt.1.0d-6) then
          rhf = .true.
        elseif (abs(abs(occmin)-1.0d+00).lt.1.0d-6  .and.
     &          abs(occmax-1.0d+00).lt.1.0d-6) then
          uhf = .true.
        elseif (abs(occmin-1.0d+00).lt.1.0d-6 .and.
     &          abs(occmax-2.0d+00).lt.1.0d-6) then
          rohf = .true.
        endif
c
        cidet(:) = 0.0d+00
        ndets = 1
        onedcou = .true.
        cidet(0) = 1.0d+00
c
c-------Test the total number of electrons and print a warning if
c       it is not an integer number.
c
        nelaux = int(anint(totel))
        epsocc = 1.0d-06
        if (abs(totel-dble(nelaux)).gt.epsocc) then
          write (uerr,338) totel,nelaux,epsocc
        endif
        nel = int(anint(totel))
        nelact = nel
        ncore  = 0
        nact   = nmo
        if (.not.(rhf.or.rohf.or.uhf)) then
          write (uerr,332)
        endif
        if (nalpha.gt.0) then
          cidet(1:nalpha) = +ialpha(1:nalpha)
        endif
        if (nbeta.gt.0) then
          cidet(nalpha+1:nalpha+nbeta) = -ibeta(1:nbeta)
        endif
        close (iwfn)
      endif
 9999 continue
c
c-----MOs in the WFN file are expressed in terms of unnormalized carte-
c     sian gaussians (CG). It could happens that a given CG is repeated.
c     In that case, we condensate all the coefficients in a single one.
c     
      if (.not.allocated(oexpa)) then
        allocate (oexpa(nprims),stat=ier)
        if (ier.ne.0) stop 'rdwfn.f: Cannot allocate oexpa()'
      endif
      if (.not.allocated(coefa)) then
        allocate (coefa(nmo+nmo,nprims),stat=ier)
        if (ier.ne.0) stop 'rdwfn.f: Cannot allocate coefa()'
      endif
      if (.not.allocated(coefat)) then
        allocate (coefat(nmofull,nprims),stat=ier)
        if (ier.ne.0) stop 'rdwfn.f: Cannot allocate coefat()'
      endif
      if (.not.allocated(icena)) then
        allocate (icena(nprims),stat=ier)
        if (ier.ne.0) stop 'rdwfn.f: Cannot allocate icena()'
      endif
      if (.not.allocated(itypa)) then
        allocate (itypa(nprims),stat=ier)
        if (ier.ne.0) stop 'rdwfn.f: Cannot allocate itypa()'
      endif
      nprimsold=0
      do 10 j=1,nprims
        if (j.eq.1) then
          nprimsold=nprimsold+1
          icena(nprimsold)=icen(j)
          oexpa(nprimsold)=oexp(j)
          itypa(nprimsold)=ityp(j)
          do i=1,nmotot
            coefa(i,nprimsold)=coef(i,j)
          enddo
        else
          do m=1,nprimsold
            if (icen(j).eq.icena(m) .and.
     &          ityp(j).eq.itypa(m) .and.
     &          abs(oexp(j)-oexpa(m)).le.1.0d-10)  then
                do i=1,nmotot
                  coefa(i,m)=coefa(i,m)+coef(i,j)
                enddo
                goto 10
            endif
          enddo
          nprimsold=nprimsold+1
          icena(nprimsold)=icen(j)
          oexpa(nprimsold)=oexp(j)
          itypa(nprimsold)=ityp(j)
          do i=1,nmotot
            coefa(i,nprimsold)=coef(i,j)
          enddo
        endif
 10   continue
c
c-----Recompute the original variables
c
      write (uout,*)'# '
      write (uout,311) nprims,nprimsold
 311  format (1x,'# Input number of Primitives ',I8,' reduced to ',I8)
      nprims=nprimsold
c
c-----Deallocate and allocate again those arrays depending on 'nprims',
c     whose dimensions might have changed with respect to their initial 
c     values due to the above line 'nprims=nprimsold'
c
      deallocate (oexp,icen,ityp,coef,xnorm,icenat)
      allocate (oexp(nprims),stat=ier)
      if (ier.ne.0) stop 'rdwfn.f: Cannot allocate oexp() array ?'
      allocate (icen(nprims),stat=ier)
      if (ier.ne.0) stop 'rdwfn.f: Cannot allocate icen() array ?'
      allocate (ityp(nprims),stat=ier)
      if (ier.ne.0) stop 'rdwfn.f: Cannot allocate ityp() array ?'
      allocate (coef(nmo+nmo,nprims),stat=ier)
      if (ier.ne.0) stop 'rdwfn.f: Cannot allocate coef() array ?'
      allocate (xnorm(nprims),stat=ier)
      if (ier.ne.0) stop 'rdwfn.f: Cannot allocate xnorm() array ?'
      allocate (icenat(nprims,ncent),stat=ier)
      if (ier.ne.0) stop 'rdwfn.f: Can not allocate icenat() array ?'

      do j=1,nprims
        icen(j)=icena(j)
        oexp(j)=oexpa(j)
        ityp(j)=itypa(j)
        do i=1,nmotot
          coef(i,j)=coefa(i,j)
        enddo
      enddo
c
c.....Determine primitives corresponding to each center.
c
      do ic=1,ncent
        npc(ic)=0
      enddo
      do j=1,nprims
        ic=icen(j)
        npc(ic)=npc(ic)+1
        inda=npc(ic)
        icenat(inda,ic)=j
      enddo
c
c.....Classify primitives in each center by types and similar exponents.
c
      ngroup(1:ncent)=0
      do ic=1,ncent
        if (npc(ic).eq.0) cycle
        do 33 j=1,npc(ic)
          k=icenat(j,ic)
          itk=ityp(k)
          isuk=nlm(itk,1)+nlm(itk,2)+nlm(itk,3)
          if (j.eq.1) then
            ngroup(ic)=1
            nzexp(ic,1)=1
            nuexp(ic,1,1)=k
          else
            do m=1,ngroup(ic)
              inda=nuexp(ic,m,1)
              itd=ityp(inda)
              isud=nlm(itd,1)+nlm(itd,2)+nlm(itd,3)
              if (abs(oexp(k)-oexp(inda)).lt.1.0d-8) then
                if (itk.eq.1.and.itd.eq.1) then
                  write (uerr,*) 
     &            'Warning TWO S primitives with equal exponents'
                   stop
                else
                  if (isuk.eq.isud) then
                    nzexp(ic,m)=nzexp(ic,m)+1
                    nuexp(ic,m,nzexp(ic,m))=k
                    goto 33
                  endif
                endif
              endif
            enddo
            ngroup(ic)=ngroup(ic)+1
            if (ngroup(ic).gt.mgrp) then
              write (uerr,*) 'Faterr: Increase MGRP in primgto.inc file'
              stop
            endif
            nzexp(ic,ngroup(ic))=1
            nuexp(ic,ngroup(ic),1)=k
          endif   
 33     continue
      enddo
c
c.....Reconstruct the values of ityp(),icen(),oexp(), and coef()
c
      i = 0
      do ic=1,ncent
        if (ngroup(ic).eq.0) cycle
        do m=1,ngroup(ic)
          do k=1,nzexp(ic,m)
            j = nuexp(ic,m,k)
            alph = oexp(j)
            itip = ityp(j)
            icentro = icen(j)
            i = i + 1
            itypa (i) = itip
            oexpa (i) = alph
            icena (i) = icentro
            do n=1,nmotot
              coefa(n,i) = coef(n,j)
            enddo
          enddo
        enddo
      enddo
c 
      do i = 1, nprims
        ityp(i) = itypa (i)
        oexp(i) = oexpa (i)
        icen(i) = icena (i)
        coef(1:nmotot,i) = coefa(1:nmotot,i)
c
c.......Normalization constants of Cartesian Gaussians.
c
        it1=nlm(ityp(i),1)
        it2=nlm(ityp(i),2)
        it3=nlm(ityp(i),3)
        it1d=it1+it1-1
        it2d=it2+it2-1
        it3d=it3+it3-1
        isu=it1+it2+it3
        dblisu2=(isu+isu+3)*0.25d0
        xnorm(i)=(2.0d+00/pi)**0.75d0*2.0d+00**isu*oexp(i)**dblisu2
        xnorm(i)=xnorm(i)/sqrt(facd(it1d)*facd(it2d)*facd(it3d))
      enddo
c
      npcant = 0
      do ic=1,ncent
        if (npc(ic).eq.0) cycle
        do k=1,npc(ic)
          icenat(k,ic) = k + npcant
        enddo
        npcant = npcant + npc(ic)
      enddo
c
c.....Reconstruct the values of nuexp(). Now, they are ordered.
c.....Determine also which is the maximum value of ngroup(ic)
c
      i = 0
      numshells = 0
      maxgrp=0
      do ic=1,ncent
        numshells = numshells + ngroup(ic)
        if (ngroup(ic).gt.maxgrp) maxgrp=ngroup(ic)
        do m=1,ngroup(ic)
           do k=1,nzexp(ic,m)
             i = i + 1
             nuexp(ic,m,k) = i
           enddo
        enddo
      enddo
c
      write (uout,*)'# '
      write (uout,*)'# Description of the Primitive Basis Set'
      write (uout,604) nprims
      do ic=1,ncent
        ntfu(0:4) = 0
        do m=1,ngroup(ic)
          zz = oexp(nuexp(ic,m,1))
          ii = nuexp(ic,m,1)
          itip = ityp(ii)
          it1 = nlm(itip,1)
          it2 = nlm(itip,2)
          it3 = nlm(itip,3)
          isu = it1+it2+it3
          ntfu(isu) = ntfu(isu)+1
          nzicm = nzexp(ic,m)
c
c---------Re-ordering primitives in each shell such their ityp() values are 
c         ordered by increasing value; i.e. 2 3 4 instead of 4 2 3, 
c         5 6 7 8 9 10 instead of 8 5 7 9 10 6, etc.
c
          if (nzicm.gt.1) then  ! i.e. only p,d,f,.. primitives
            forall (i=1:nzicm) iord(i)=i
            ordtyp(1:nzicm) = ityp(nuexp(ic,m,1:nzicm))
            call iqcksort (ordtyp,iord,nzicm,1,nzicm)
            ityp(nuexp(ic,m,1:nzicm))=ordtyp(iord(1:nzicm))
            do k=1,nmotot
              cord(1:nzicm) = coef(k,nuexp(ic,m,1:nzicm))
              forall (i=1:nzicm) coef(k,nuexp(ic,m,i)) = cord(iord(i))
            enddo
          endif
        enddo
        write (uout,2120) ic,(ntfu(nn),nn=0,4)
      enddo
c
c.....construct supervector density matrix
c
!     call allocate_space_for_pmat (nprims)
c
c-----Be careful if the WFN file comes from a MOLPRO CCSD calculation 
c
      if (cciqa) then
        if (.not.allocated(v1mata)) then
          allocate (v1mata(nmo,nmo),stat=ier)
          if (ier.ne.0) stop 'rdwfn.f: Cannot allocate v1mata()'
        endif
        if (.not.allocated(d1mata)) then
          allocate (d1mata(nmo),stat=ier)
          if (ier.ne.0) stop 'rdwfn.f: Cannot allocate d1mata()'
        endif
        c1et(1:nmo,1:nmo)=0.0d+00
c
c-------Re-read lines up to the definition of 1-RDM and 2-RDM
c
        open (iwfn,file=trim(wfnfile),status='old')         
        do i=1,nrecwfn
          read (udat,'(a)') line
        enddo
 1250   read (udat,*,err=1230) ifil,icol,rho1val
        if (ifil.gt.nmo.or.icol.gt.nmo) then
          stop 'ccwfn.f: Row and/or column index > NMO in input file'
        else
          c1et(ifil,icol)=rho1val
        endif
        goto 1250
 1230   close (iwfn)
c
c-------Test the symmetric character of c1et().
c
        do i=2,nmo
          do j=1,i-1
            cij=c1et(i,j)
            cji=c1et(j,i)
            if (abs(cij).gt.1.0D-8) then
              if (abs(cji).gt.1.0D-8) then
                 if (abs(cij-cji).gt.1.0D-6) then
                   write (uerr,600) i,j,cij,j,i,cji
                 endif
                 c1et(i,j) = 0.5d0 * (cij+cji)
                 c1et(j,i) = c1et(i,j)
              else
                 c1et(j,i) = c1et(i,j)
              endif
            else
              if (abs(cji).gt.1.0D-8) then
                 c1et(i,j) = c1et(j,i)
              else
                 c1et(i,j) = 0.5d0 * (cij+cji)
                 c1et(j,i) = c1et(i,j)
              endif
            endif
          enddo
        enddo
        totel = 0.0d+00
        gamma = 0.0d+00
        do i=1,nmo
          totel = totel + c1et(i,i)
        enddo
        nel=int(anint(totel))
c
c-------diagonalize total 1st-order matrix
c
c-------Store also c1et() in onerdm() to not destroy it.
c
        onerdm(1:nmo,1:nmo)=c1et(1:nmo,1:nmo)
        call jacobi (c1et,nmo,nmo,d1mata,v1mata,nrot)
        occ(1:nmo)=d1mata(1:nmo)
        do i=1,nmo
          do j=1,nprims
            tmp = 0.0d+00
            do k=1,nmo
              tmp = tmp + v1mata(k,i) * coef(k+nmo,j)
            enddo
            coef(i,j) = tmp
          enddo
        enddo
        if (allocated(v1mata)) then
          deallocate (v1mata,stat=ier)
          if (ier.ne.0) stop 'rdwfn.f: Cannot deallocate v1mata()'
        endif
        if (allocated(d1mata)) then
          deallocate (d1mata,stat=ier)
          if (ier.ne.0) stop 'rdwfn.f: Cannot deallocate d1mata()'
        endif
      else
      endif
c
c.....Deallocate arrays.
c
      if (allocated(oexpa)) then
        deallocate (oexpa,stat=ier)
        if (ier.ne.0) stop 'rdwfn.f: Cannot deallocate oexpa()'
      endif
      if (allocated(coefa)) then
        deallocate (coefa,stat=ier)
        if (ier.ne.0) stop 'rdwfn.f: Cannot deallocate coefa()'
      endif
      if (allocated(coefat)) then
        deallocate (coefat,stat=ier)
        if (ier.ne.0) stop 'rdwfn.f: Cannot deallocate coefat()'
      endif
      if (allocated(icena)) then
        deallocate (icena,stat=ier)
        if (ier.ne.0) stop 'rdwfn.f: Cannot deallocate icena()'
      endif
      if (allocated(itypa)) then
        deallocate (itypa,stat=ier)
        if (ier.ne.0) stop 'rdwfn.f: Cannot deallocate itypa()'
      endif

      return
101   format (A80)
102   format (4X,A4,10X,3(I5,15X))
104   format (20X,20I3)
105   format (10X,5E14.7)
107   format (5E16.8)
108   format(A8)
2120  format (1x,'# CENTER ',i4,' : [ ',
     &   I2,'s | ',I2,'p |',I2,'d | ',I2,'f | ',I2,'g ] basis')
604   format (' # Total number of Primitive Gaussians: ',I6)
331   format (/,' # rdwfn.f: !!! Number of e- is not integer',/,
     & ' # rdwfn.f: Number of electrons = ',F15.8,/,' #')
332   format (/,
     & ' # rdwfn.f: SDW which is not of RHF, ROHF, UHF type',/,
     & ' # Some of the modules of PROMOLDEN will not properly work',
     & ' in case that canonical MOs are necessary')
182   format (1x,'# rdwfn.f: Error reading charge() of line ',
     &  I2,' of WFN file')
600   format (' # ¡¡¡ Warning 1-RDM is not symmetric !!!',/,
     &    ' # 1-RDM(',I3,',',I3,' ) = ',E17.10,' but ',/,
     &    ' # 1-RDM(',I3,',',I3,' ) = ',E17.10)
873   format(' # rdwfn.f: Bad format of line before reading ',a)
338   format (' #',/,' # ',6('WARNING -- '),/,' #',54(' '),/,
     & ' # Number of electrons from the Occupation Numbers =',F15.8,/
     & ' # The closest integer to the above number is      =',I15,/,
     & ' # These two numbers differ in a quantity         >=',E15.8,/   
     & ' # ',53(' '),/,' # ',6('WARNING -- '))
      end
