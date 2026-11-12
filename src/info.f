      subroutine info (filedat)
c
c.....initialize some variables
c
      USE                 space_for_wfnbasis
      USE                 space_for_corr
      include            'implicit.inc'
      include            'param.inc'
      include            'wfn.inc'
      include            'datatm.inc'
      include            'io.inc'
      include            'corr.inc'
c
      character*(132)  filedat
      integer            nform(0:106), nn
      integer,allocatable,dimension (:)   :: ialps,ibets
c
c.....write some info about molecule
c
      nl=131
      do while (filedat(nl:nl).eq.' ')
         nl=nl-1
      enddo
      write (uout,10) filedat(1:nl), ncent
      call cmcq ()
      do i=0,106
         nform(i)=0
      enddo
      wmol=0d0
c
      write (uout,11, advance='no')
      do i=1,ncent
         nz=int(charget(i))
         if (nz.le.103) wmol=wmol+wgatm(nz)
         nform(nz)=nform(nz)+1
      enddo
      nn=0
      do i=0,106
         nn=nform(i)
         if (nn.ne.0) then 
           write (uout,12,advance='no') namatm(i)(1:2),'(',nn,')'     
         end if
         nn=0
      enddo
      write (uout,13) wmol
c
      write (uout,20)
      do i=1,ncent
         write (uout,21) i, atnam(i)(1:4),(xyz(i,j),j=1,3)
      enddo
c
c-----Maximum angular momentum in the basis of primitive Gaussians
c
      maxlprim = 0
      do i=1,nprims
         ii = ityp(i)
         ni = nlm(ii,1)+nlm(ii,2)+nlm(ii,3)
         maxlprim = max(ni,maxlprim)
      enddo
      write (uout,30) npor, nprims, nmo, maxlprim
c
      allocate (ialps(nmo),stat=ier)
      if (ier.ne.0) stop 'info.f: Cannot allocate ialps'
      allocate (ibets(nmo),stat=ier)
      if (ier.ne.0) stop 'info.f: Cannot allocate ibets'
c
      if (icorr) then
         write (uout,41) 
         chtot=anint(sum(charge(1:ncent))-anint(totel))
         write (uout,31) totalenergy,gamma,dble(anint(totel)),
     &              chtot,ndets,nelact,nact,ncore
         write (uout,270)
         write (uout,'(a)') ' # OCCs from file "'//trim(filedat)//'"'
         do i=1,nmo
           write (uout,271) i,occ(i),eorb(i)
         enddo
      else
        if (rhf) then
          write (uout,42) 
        elseif (rohf) then
          write (uout,43) 
        elseif (uhf) then
          write (uout,44) 
        else
          write (uout,45) 
        endif
        chtot=anint(sum(charge(1:ncent))-anint(totel))
        if (cciqa.or.mp2iqa) then
          write (uout,460) totalenergy,gamma,dble(anint(totel)),chtot
        elseif (rhf.or.uhf.or.rohf) then
          write (uout,46) totalenergy,gamma,dble(anint(totel)),
     &      chtot,nsin,ndou,nalpha-ndou,nbeta-ndou
        else
          write (uout,461) totalenergy,gamma,dble(anint(totel)),chtot
        endif
        write (uout,270)
        if (cciqa.or.mp2iqa) then
          write (uout,'(a)') ' # OCCs Computed in the "rdwfn.f" routine'
        else
          write (uout,'(a)') ' # OCCs from file "'//trim(filedat)//'"'
        endif
        do i=1,nmo
          write (uout,271) i,occ(i),eorb(i)
        enddo
        if (uhf.or.rohf) then
          nsal=0
          do i=1,nalpha
            ia=ialpha(i)
            do j=1,nbeta
              if (ia.eq.ibeta(j)) goto 49
            enddo
            nsal=nsal+1
            ialps(nsal)=ia
 49         continue
          enddo
          write (uout,47) 'ALPHA',(ialps(i),i=1,nsal)
          nsbe=0
          do i=1,nbeta
            ib=ibeta(i)
            do j=1,nalpha
              if (ib.eq.ialpha(j)) goto 59
            enddo
            nsbe=nsbe+1
            ibets(nsbe)=ib
 59         continue
          enddo
          write (uout,47) 'BETA ',(ibets(i),i=1,nsbe)
        endif
      endif
      deallocate (ialps,stat=ier)
      if (ier.ne.0) stop 'info.f: Cannot deallocate ialps'
      deallocate (ibets,stat=ier)
      if (ier.ne.0) stop 'info.f: Cannot deallocate ibets'
      return
c
c----------------------------------------------------------------
 10   format (1x, '# WAVEFUNCTION ACCOUNT:',/,
     &        1x, '# File: ', a,/
     &        1x, '# Number of Centers:' ,i4)
 11   format (1x, '# Molecule: ')
 13   format (/1x, '# Molecular Weight:', f14.6)
 12   format (a2,a1,i4,a1,1x)
 20   format (1x,'# Cartesian coordinates of centers:'/)
 21   format (6x,i3,1x,a4,1x,3f14.8)
 30   format (/,1x,'# Number of Primitives         : ',i4,4x,
     &             'reduced to ',I4,/
     &        1x,'# Number of Molecular Orbitals : ', i4,/
     &        1x,'# Maximum l in Basis Set       : ', i4,/
     &        1x,'#',72('-'))
 41   format (1x, '# MULTIDETERMINANTAL WAVEFUNCTION')
 42   format (1x, '# RESTRICTED CLOSED SHELL HF WAVEFUNCTION (RHF)')
 43   format (1x, '# GENERAL MONODETERMINANTAL WAVEFUNCTION')
 44   format (1x, '# UNRESTRICTED HF WAVEFUNCTION (UHF)')
 45   format (1x, '# WAVEFUNCTION IS NOT RHF, ROHF, NOR UHF')
 31   format (1x, '# Total energy                 :', F20.10,/
     &        1x, '# Virial                       :', F13.8,/
     &        1x, '# Number of Electrons          :', F16.10,/
     &        1x, '# Total charge of the molecule :', F16.10,/
     &        1x, '# Number of Determinants       :', i9,/
     &        1x, '# Number of Active Electrons   :', i9,/
     &        1x, '# Number of Active Orbitals    :', i9,/
     &        1x, '# Number of Core Orbitals      :', i9,/
     &        1x,'#',72('-'))
 460  format (1x, '# Total energy                          :', F20.10,/
     &        1x, '# Virial                                :', F13.8,/
     &        1x, '# Number of Electrons                   :', F16.10,/
     &        1x, '# Total charge of the molecule          :',F16.10,/,
     &        1x,'#',72('-'))
 46   format (1x, '# Total energy                          :', F20.10,/
     &        1x, '# Virial                                :', F13.8,/
     &        1x, '# Number of Electrons                   :', F16.10,/
     &        1x, '# Total charge of the molecule          :', F16.10,/
     &        1x, '# Number of Singly occupied MOs         :', i9,/
     &        1x, '# Number of Double occupied MOs         :', i9,/
     &        1x, '# Number of Open-Shell (ALPHA/BETA) MOs :',2i9,/
     &        1x,'#',72('-'))
 461  format (1x, '# Total energy                          :', F20.10,/
     &        1x, '# Virial                                :', F13.8,/
     &        1x, '# Number of Electrons                   :', F16.10,/
     &        1x, '# Total charge of the molecule          :', F16.10,/
     &        1x,'#',72('-'))
 47   format (1x, '# Open Shell ',a,' MOs are: ',50I3,/,
     &   1000(25x,50I3,/))
 270  format (' # Natural MOs occupations (OCC) & Orbital energies')
 271  format (' # MO ',I4,3x,'OCC = ',F15.8,3x,'ORB ENERGY = ',f15.8)
c----------------------------------------------------------------
      end
