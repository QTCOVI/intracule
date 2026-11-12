module mod_rdm
  save
  real   (kind=8), allocatable, dimension(:,:)      :: rdm1
  real   (kind=8), allocatable, dimension(:,:,:,:)  :: rdm2
  real   (kind=8), allocatable, dimension(:,:)      :: rdm1a,rdm1b
  real   (kind=8), allocatable, dimension(:,:,:,:)  :: rdm2aa,rdm2bb,rdm2ab,rdm2ba
  real   (kind=8), allocatable, dimension(:,:)      :: c2e,ex2mat

contains

  subroutine compute_rdm (ndets,nelact,ncore,nmo,uout,uerr,nrecwfn,&
      cciqa,mp2iqa,icorr,verbose,filename)

    use space_for_corr
    implicit none

    integer(kind=4),  intent(in) :: ndets,nelact,ncore,nmo
    integer(kind=4),  intent(in) :: uout,uerr
    integer(kind=4),  intent(in) :: nrecwfn
    logical,          intent(in) :: cciqa,mp2iqa,icorr,verbose
    character(len=*), intent(in) :: filename

    real   (kind=8) :: cd(ndets)
    integer(kind=8) :: det(2,ndets)
    integer(kind=4) :: nela,nelb
    integer(kind=4) :: i, j, k, l, i1, j1, i2, j2, ij, kl
    integer(kind=4) :: maxij, minij, maxkl, minkl, ii, jj
    integer(kind=4) :: ier
    real   (kind=8) :: xcrho2,aklmn

    integer(kind=4) :: udat,ifil,icol
    real   (kind=8) :: value,cij,cji,epsrho
    real   (kind=8) :: xnorm
    real   (kind=8) :: dkl,dmn
    character*(132) :: line
!
!---dmft.inc
!
    logical dmftxc
    common /typedmftxc/ dmftxc
    integer idmft
    common /typexc/ idmft
!
!---1-RDM and 2-RDM obtained with the Müller approximation
!
    if (dmftxc) then
      allocate (rdm1 (nmo,nmo))
      allocate (rdm2 (nmo,nmo,nmo,nmo))
      call muller (nmo, rdm1, rdm2)
      return
    endif

    if (.not.allocated(rdm1)) then
      allocate (rdm1 (nmo,nmo),stat=ier)
      if (ier.ne.0) stop '# mod_rdm.f90: Cannot allocate rdm1()'
    endif
    if (.not.allocated(rdm1a)) then
      allocate (rdm1a(nmo,nmo),stat=ier)
      if (ier.ne.0) stop '# mod_rdm.f90: Cannot allocate rdm1a()'
    endif
    if (.not.allocated(rdm1b)) then
      allocate (rdm1b(nmo,nmo),stat=ier)
      if (ier.ne.0) stop '# mod_rdm.f90: Cannot allocate rdm1b()'
    endif
    if (.not.allocated(rdm2)) then
      allocate (rdm2  (nmo,nmo,nmo,nmo),stat=ier)
      if (ier.ne.0) stop '# mod_rdm.f90: Cannot allocate rdm2()'
    endif
    if (.not.allocated(rdm2aa)) then
      allocate (rdm2aa(nmo,nmo,nmo,nmo),stat=ier)
      if (ier.ne.0) stop '# mod_rdm.f90: Cannot allocate rdm2aa()'
    endif
    if (.not.allocated(rdm2bb)) then
      allocate (rdm2bb(nmo,nmo,nmo,nmo),stat=ier)
      if (ier.ne.0) stop '# mod_rdm.f90: Cannot allocate rdm2bb()'
    endif
    if (.not.allocated(rdm2ab)) then
      allocate (rdm2ab(nmo,nmo,nmo,nmo),stat=ier)
      if (ier.ne.0) stop '# mod_rdm.f90: Cannot allocate rdm2ab()'
    endif
    if (.not.allocated(rdm2ba)) then
      allocate (rdm2ba(nmo,nmo,nmo,nmo),stat=ier)
      if (ier.ne.0) stop '# mod_rdm.f90: Cannot allocate rdm2ba()'
    endif
    if (.not.allocated(c2e)) then
      allocate (c2e(nmo*(nmo+1)/2,nmo*(nmo+1)/2),stat=ier)
      if (ier.ne.0) stop '# mod_rdm.f90: Cannot allocate c2e()'
    endif
    if (.not.allocated(ex2mat)) then
      allocate (ex2mat(nmo*(nmo+1)/2,nmo*(nmo+1)/2),stat=ier)
      if (ier.ne.0) stop '# mod_rdm.f90: Cannot allocate ex2mat()'
    endif

    rdm1   = 0.0d+00
    rdm1a  = 0.0d+00
    rdm1b  = 0.0d+00
    rdm2   = 0.0d+00
    rdm2aa = 0.0d+00
    rdm2bb = 0.0d+00
    rdm2ab = 0.0d+00
    rdm2ba = 0.0d+00
    c2e    = 0.0d+00
    ex2mat = 0.0d+00

    udat = 1
    if (cciqa.or.mp2iqa) then
      open (udat,file=trim(filename),status='old')
      do i=1,nrecwfn
        read (udat,'(a)') line
      enddo
 125  read (udat,*,err=123) ifil,icol,value
      if (ifil.gt.nmo.or.icol.gt.nmo) then
        stop 'mod_rdm.f90: Row and/or column index > NMO in the input file'
      else
        rdm1(ifil,icol)=value
      endif
      goto 125
 123  read (udat,*,end=124) i,j,k,l,value
      if (i.gt.nmo.or.j.gt.nmo.or.k.gt.nmo.or.l.gt.nmo) then
        stop 'mod_rdm.f90: Row and/or column index > NMO in the input file'
      else
        rdm2 (i,j,k,l) = value
      endif
      goto 123
 124  close (udat)
!
!-----Symmetrize 1-RDM
!
      do i=2,nmo
        do j=1,i-1
          cij=rdm1(i,j)
          cji=rdm1(j,i)
          if (abs(cij).gt.1.0D-8) then
            if (abs(cji).gt.1.0D-8) then
               if (abs(cij-cji).gt.1.0D-6) then
                 write (uerr,600) i,j,cij,j,i,cji
               endif
               rdm1(i,j) = 0.5d0 * (cij+cji)
               rdm1(j,i) = rdm1(i,j)
            else
               rdm1(j,i) = rdm1(i,j)
            endif
          else
            if (abs(cji).gt.1.0D-8) then
               rdm1(i,j) = rdm1(j,i)
            else
               rdm1(i,j) = 0.5d0 * (cij+cji)
               rdm1(j,i) = rdm1(i,j)
            endif
          endif
        enddo
      enddo
    else
      if (icorr) then
!
!-------The maximum number of allowed MOs in CAS calculations is 64
!
        if (nmo.ge.64) then
          write (uerr,*) '# Maximum number of allowed MOs is 64'
          write (uerr,*) '# The current value is ',nmo
          stop
        endif

        open (udat,file=trim(filename),status='old')
        do i=1,nrecwfn
          read (udat,'(a)') line
        enddo
        xnorm = 0.0d+00
        nela = ncore
        nelb = ncore
        det(:,:) = 0
        do i=1,ndets
          read (udat,*) cd(i),(cidet(j),j=1,nelact)
          if (i.eq.1) then
            do j = 1,nelact
              if (cidet(j).gt.0) nela = nela + 1
              if (cidet(j).lt.0) nelb = nelb + 1
            end do
          end if
          xnorm = xnorm + cd(i)*cd(i)
          if (ncore.gt.0) then
            do j=1,ncore
              det(1,i) = ibset(det(1,i),j-1)
              det(2,i) = ibset(det(2,i),j-1)
            enddo
          endif
          do j = 1,nelact
            if (cidet(j).gt.0) then
              det(1,i) = ibset(det(1,i),ncore+int(cidet(j))-1)
            else
              det(2,i) = ibset(det(2,i),ncore+abs(int(cidet(j)))-1)
            end if
          end do
        enddo
        xnorm = 1.0d+00/sqrt(xnorm)
        cd(1:ndets) = cd(1:ndets) * xnorm
        close (udat)
      else
        cd(1) = 1.0d+00
        nela = 0
        nelb = 0
        do j=1,nelact
          if (int(cidet(j)).gt.0) then
            nela = nela + 1
            ialpha(nela)=+int(cidet(j))
          endif
          if (int(cidet(j)).lt.0) then
            nelb = nelb + 1
            ibeta(nelb)=-int(cidet(j))
          endif
        enddo
        det(:,:) = 0
        do j=1,nela
          det(1,1) = ibset(det(1,1),ialpha(j)-1)
        enddo
        do j=1,nelb
          det(2,1) = ibset(det(2,1),ibeta(j)-1)
        enddo
      endif
      call compute_rdm2 (cd,det,ndets,nela,nelb)
      call compute_rdm1 (cd,det,ndets,nela,nelb)
    endif
!
!---Print results
!
    if (verbose) then
      epsrho = 1.0d-14
      write (uout,*) '1-RDM'
      write (uout,*) '    e1    e1'
      write (uout,*) '     i     j     rdm1(i,j)'
      write (uout,*) '--------------------------------------'
      do i=1,nmo
        do j=1,nmo
          if (abs(rdm1(i,j)).gt.epsrho) then
             write (uout,'(1x,2(1x,i5),1x,e25.18)') i,j,rdm1(i,j)
          endif
        enddo
      enddo
      write (uout,*) '2-RDM'
      write (uout,*) '    e1    e1    e2    e2'
      write (uout,*) '     i     j     k     l      rdm2(i,j,k,l)'
      write (uout,*) '--------------------------------------------------'
      do i=1,nmo
        do j=1,nmo
          do k=1,nmo
            do l=1,nmo
              if (abs(rdm2(i,j,k,l)).gt.epsrho) then
                 write (uout,'(1x,4(1x,i5),1x,e25.18)') i,j,k,l,rdm2(i,j,k,l)
              endif
            enddo
          enddo
        enddo
      enddo
    endif
!
!---Condense rdm2() ---> c2e()
!
    do i = 1,nmo
      do j = 1,nmo
        do k = 1,nmo
          do l = 1,nmo
            maxij=max(i,j)
            maxkl=max(k,l)
            minij=min(i,j)
            minkl=min(k,l)
            ii=maxij*(maxij-1)/2+minij
            jj=maxkl*(maxkl-1)/2+minkl
            c2e(ii,jj)=c2e(ii,jj)+rdm2(i,j,k,l)
          end do
        end do
      end do
    end do
!
!-----Symmetrize c2e()
!
      do i=2,nmo*(nmo+1)/2
        do j=1,nmo*(nmo+1)/2-1
          cij=c2e(i,j)
          cji=c2e(j,i)
          if (abs(cij-cji).gt.1.0d-6) then
            write (uerr,*) i,j,cij,cji
            write (uerr,*) 'mod_rdm.f90: WARNING Non symmetric c2e() array'
          else
            c2e(i,j)=0.5d0*(cij+cji)
            c2e(j,i)=0.5d0*(cij+cji)
          endif
        enddo
      enddo
!
      if (verbose) then
        write (uout,*) '# '
        write (uout,*) '# Condensed 2-RDM in the Canonical MO basis'
        write (uout,54) 
        ij = 0
        do i1=1,nmo
          do j1=1,i1
             ij = ij + 1
             dkl = 2.0d+00
             if (i1.eq.j1) dkl=1.0d+00
             kl = 0
             do i2=1,nmo
               do j2=1,i2
                  kl = kl +1 
                  if (kl.gt.ij) cycle
                  dmn=2.0d+00
                  if (i2.eq.j2) dmn=1.0d+00
                  aklmn=rdm1(i1,j1) * rdm1(i2,j2) * dkl * dmn
                  xcrho2=c2e(ij,kl)-aklmn
                  ex2mat(ij,kl) = xcrho2
                  ex2mat(kl,ij) = xcrho2
                  if (abs(c2e(ij,kl)).gt.epsrho .or. abs(aklmn) .gt.epsrho .or. &
                      abs(xcrho2)  .gt.epsrho) then
                      write (uout,55) i1,j1,i2,j2,c2e(kl,ij),xcrho2,aklmn
                  endif
               enddo
             enddo
          enddo
        enddo
      endif
 54 format (' #',12x,'e1    e1    e2    e2',10x,'TOT',23x,'XC',24x,'Coul')
 55 format (' #  rho_2 ',4(1x,i5),3x,3(1x,E25.18))
600 format (' # ¡¡¡ Warning 1-RDM is not symmetric !!!',/,&
     ' # 1-RDM(',I3,',',I3,' ) = ',E17.10,' but ',/,&
     ' # 1-RDM(',I3,',',I3,' ) = ',E17.10)

  end subroutine compute_rdm

  subroutine compute_rdm1 (cd,det,ndets,nela,nelb)

    implicit none

    integer(kind=4), intent(in) :: ndets,nela,nelb
    real   (kind=8), intent(in) :: cd(ndets)
    integer(kind=8), intent(in) :: det(2,ndets)

    integer(kind=4) :: dega, degb, ia(2), i, j, k, l
    integer(kind=4) :: p, occa(nela), occb(nelb)
    real(kind=8) :: phase, c
    
    rdm1  = 0.0d+00
    rdm1a = 0.0d+00
    rdm1b = 0.0d+00
    do k = 1,ndets
      call occ_list (det(1,k), nela, occa)
      call occ_list (det(2,k), nelb, occb)
      c = cd(k)*cd(k)
      do p = 1,nela
        j = occa(p)
        rdm1(j,j) = rdm1(j,j) + c
        rdm1a(j,j) = rdm1a(j,j) + c
      end do
      do p = 1,nelb
        j = occb(p)
        rdm1(j,j) = rdm1(j,j) + c
        rdm1b(j,j) = rdm1b(j,j) + c
      end do
      ! Non diagonal terms
      do l = 1,k-1
        dega = n_excitations(det(1,k),det(1,l))
        degb = n_excitations(det(2,k),det(2,l))
        if (dega+degb > 1) then
          cycle
        else if (degb == 0) then
          call get_single_excitation (det(1,k),det(1,l),ia)
          i = ia(1)
          j = ia(2)
          phase = cre_des_sign(j,i,det(1,k))
          c = phase*cd(k)*cd(l)
          rdm1(j,i) = rdm1(j,i) + c 
          rdm1(i,j) = rdm1(i,j) + c 
          rdm1a(j,i) = rdm1a(j,i) + c 
          rdm1a(i,j) = rdm1a(i,j) + c 
        else if (dega == 0) then
          call get_single_excitation (det(2,k),det(2,l),ia)
          i = ia(1)
          j = ia(2)
          phase = cre_des_sign(j,i,det(2,k))
          c = phase*cd(k)*cd(l)
          rdm1(j,i) = rdm1(j,i) + c 
          rdm1(i,j) = rdm1(i,j) + c 
          rdm1b(j,i) = rdm1b(j,i) + c 
          rdm1b(i,j) = rdm1b(i,j) + c 
        end if
      end do
    end do

  end subroutine compute_rdm1

  subroutine compute_rdm2 (cd,det,ndets,nela,nelb)

    implicit none

    integer(kind=4), intent(in) :: ndets,nela,nelb
    real   (kind=8), intent(in) :: cd(ndets)
    integer(kind=8), intent(in) :: det(2,ndets)

    integer(kind=4) :: dega, degb, ia(2), i, j, k, l, iajb(4), jb(2)
    integer(kind=4) :: occa(nela), occb(nelb), p, pp, q, qq, a, b
    real(kind=8) :: phase, c
    

    rdm2   = 0.0d+00
    rdm2aa = 0.0d+00
    rdm2bb = 0.0d+00
    rdm2ab = 0.0d+00
    rdm2ba = 0.0d+00
    do k = 1,ndets
      call occ_list (det(1,k), nela, occa)
      call occ_list (det(2,k), nelb, occb)
      c = cd(k)*cd(k)
      do p = 1,nela
        pp = occa(p)
        do q = 1,nela
          qq = occa(q)
          rdm2(pp,pp,qq,qq) = rdm2(pp,pp,qq,qq) + c
          rdm2(pp,qq,qq,pp) = rdm2(pp,qq,qq,pp) - c
          rdm2aa(pp,pp,qq,qq) = rdm2aa(pp,pp,qq,qq) + c
          rdm2aa(pp,qq,qq,pp) = rdm2aa(pp,qq,qq,pp) - c
        end do
        do q = 1,nelb
          qq = occb(q)
          rdm2(pp,pp,qq,qq) = rdm2(pp,pp,qq,qq) + c
          rdm2ab(pp,pp,qq,qq) = rdm2ab(pp,pp,qq,qq) + c
        end do
      end do
      do p = 1,nelb
        pp = occb(p)
        do q = 1,nelb
          qq = occb(q)
          rdm2(pp,pp,qq,qq) = rdm2(pp,pp,qq,qq) + c
          rdm2(pp,qq,qq,pp) = rdm2(pp,qq,qq,pp) - c
          rdm2bb(pp,pp,qq,qq) = rdm2bb(pp,pp,qq,qq) + c
          rdm2bb(pp,qq,qq,pp) = rdm2bb(pp,qq,qq,pp) - c
        end do
        do q = 1,nela
          qq = occa(q)
          rdm2(pp,pp,qq,qq) = rdm2(pp,pp,qq,qq) + c
          rdm2ba(pp,pp,qq,qq) = rdm2ba(pp,pp,qq,qq) + c
        end do
      end do
      ! Non diagonal terms
      !do l = 1,k-1
      do l = 1,ndets
        dega = n_excitations (det(1,k),det(1,l))
        degb = n_excitations (det(2,k),det(2,l))
        if (dega+degb > 2) then
          cycle
        else if (dega+degb == 1) then
          if (degb == 0) then
            call get_single_excitation (det(1,k),det(1,l),ia)
            call occ_list (det(1,k), nela, occa)
            call occ_list (det(2,k), nelb, occb)
            i = ia(1)
            a = ia(2)
            phase = cre_des_sign(a,i,det(1,k))
            c = phase*cd(k)*cd(l)
            do q = 1,nela
              qq = occa(q)
              rdm2(a,i,qq,qq) = rdm2(a,i,qq,qq) + c
              rdm2(a,qq,qq,i) = rdm2(a,qq,qq,i) - c
              rdm2(qq,i,a,qq) = rdm2(qq,i,a,qq) - c
              rdm2(qq,qq,a,i) = rdm2(qq,qq,a,i) + c
              rdm2aa(a,i,qq,qq) = rdm2aa(a,i,qq,qq) + c
              rdm2aa(a,qq,qq,i) = rdm2aa(a,qq,qq,i) - c
              rdm2aa(qq,i,a,qq) = rdm2aa(qq,i,a,qq) - c
              rdm2aa(qq,qq,a,i) = rdm2aa(qq,qq,a,i) + c
            end do
            do q = 1,nelb
              qq = occb(q)
              rdm2(a,i,qq,qq) = rdm2(a,i,qq,qq) + c
              rdm2(qq,qq,a,i) = rdm2(qq,qq,a,i) + c
              rdm2ab(a,i,qq,qq) = rdm2ab(a,i,qq,qq) + c
              rdm2ab(qq,qq,a,i) = rdm2ab(qq,qq,a,i) + c
            end do
          else if (dega == 0) then
            call get_single_excitation (det(2,k),det(2,l),ia)
            call occ_list (det(1,k), nela, occa)
            call occ_list (det(2,k), nelb, occb)
            i = ia(1)
            a = ia(2)
            phase = cre_des_sign(a,i,det(2,k))
            c = phase*cd(k)*cd(l)
            do q = 1,nelb
              qq = occb(q)
              rdm2(a,i,qq,qq) = rdm2(a,i,qq,qq) + c
              rdm2(a,qq,qq,i) = rdm2(a,qq,qq,i) - c
              rdm2(qq,i,a,qq) = rdm2(qq,i,a,qq) - c
              rdm2(qq,qq,a,i) = rdm2(qq,qq,a,i) + c
              rdm2bb(a,i,qq,qq) = rdm2bb(a,i,qq,qq) + c
              rdm2bb(a,qq,qq,i) = rdm2bb(a,qq,qq,i) - c
              rdm2bb(qq,i,a,qq) = rdm2bb(qq,i,a,qq) - c
              rdm2bb(qq,qq,a,i) = rdm2bb(qq,qq,a,i) + c
            end do
            do q = 1,nela
              qq = occa(q)
              rdm2(a,i,qq,qq) = rdm2(a,i,qq,qq) + c
              rdm2(qq,qq,a,i) = rdm2(qq,qq,a,i) + c
              rdm2ba(a,i,qq,qq) = rdm2ba(a,i,qq,qq) + c
              rdm2ba(qq,qq,a,i) = rdm2ba(qq,qq,a,i) + c
            end do
          end if
        else if (dega+degb == 2) then
          if (degb == 0) then
            call get_double_excitation (det(1,k),det(1,l),iajb)
            i = iajb(1)
            a = iajb(2)
            j = iajb(3)
            b = iajb(4)
            if (a > j .or. i > b) then
              phase = cre_des_sign(a,j,det(1,k))
              phase = phase*cre_des_sign(b,i,det(1,k))
              c = phase*cd(k)*cd(l)
              rdm2(b,i,a,j) = rdm2(b,i,a,j) + c
              rdm2(b,j,a,i) = rdm2(b,j,a,i) - c
              rdm2(a,i,b,j) = rdm2(a,i,b,j) - c
              rdm2(a,j,b,i) = rdm2(a,j,b,i) + c
              rdm2aa(b,i,a,j) = rdm2aa(b,i,a,j) + c
              rdm2aa(b,j,a,i) = rdm2aa(b,j,a,i) - c
              rdm2aa(a,i,b,j) = rdm2aa(a,i,b,j) - c
              rdm2aa(a,j,b,i) = rdm2aa(a,j,b,i) + c
            else
              phase = cre_des_sign(a,i,det(1,k))
              phase = phase*cre_des_sign(b,j,det(1,k))
              c = phase*cd(k)*cd(l)
              rdm2(b,i,a,j) = rdm2(b,i,a,j) - c
              rdm2(b,j,a,i) = rdm2(b,j,a,i) + c
              rdm2(a,i,b,j) = rdm2(a,i,b,j) + c
              rdm2(a,j,b,i) = rdm2(a,j,b,i) - c
              rdm2aa(b,i,a,j) = rdm2aa(b,i,a,j) - c
              rdm2aa(b,j,a,i) = rdm2aa(b,j,a,i) + c
              rdm2aa(a,i,b,j) = rdm2aa(a,i,b,j) + c
              rdm2aa(a,j,b,i) = rdm2aa(a,j,b,i) - c
            end if
          else if (dega == 0) then
            call get_double_excitation (det(2,k),det(2,l),iajb)
            i = iajb(1)
            a = iajb(2)
            j = iajb(3)
            b = iajb(4)
            if (a > j .or. i > b) then
              phase = cre_des_sign(a,j,det(2,k))
              phase = phase*cre_des_sign(b,i,det(2,k))
              c = phase*cd(k)*cd(l)
              rdm2(b,i,a,j) = rdm2(b,i,a,j) + c
              rdm2(b,j,a,i) = rdm2(b,j,a,i) - c
              rdm2(a,i,b,j) = rdm2(a,i,b,j) - c
              rdm2(a,j,b,i) = rdm2(a,j,b,i) + c
              rdm2bb(b,i,a,j) = rdm2bb(b,i,a,j) + c
              rdm2bb(b,j,a,i) = rdm2bb(b,j,a,i) - c
              rdm2bb(a,i,b,j) = rdm2bb(a,i,b,j) - c
              rdm2bb(a,j,b,i) = rdm2bb(a,j,b,i) + c
            else
              phase = cre_des_sign(a,i,det(2,k))
              phase = phase*cre_des_sign(b,j,det(2,k))
              c = phase*cd(k)*cd(l)
              rdm2(b,i,a,j) = rdm2(b,i,a,j) - c
              rdm2(b,j,a,i) = rdm2(b,j,a,i) + c
              rdm2(a,i,b,j) = rdm2(a,i,b,j) + c
              rdm2(a,j,b,i) = rdm2(a,j,b,i) - c
              rdm2bb(b,i,a,j) = rdm2bb(b,i,a,j) - c
              rdm2bb(b,j,a,i) = rdm2bb(b,j,a,i) + c
              rdm2bb(a,i,b,j) = rdm2bb(a,i,b,j) + c
              rdm2bb(a,j,b,i) = rdm2bb(a,j,b,i) - c
            end if
          else
            call get_single_excitation (det(1,k),det(1,l),ia)
            call get_single_excitation (det(2,k),det(2,l),jb)
            i = ia(1)
            a = ia(2)
            j = jb(1)
            b = jb(2)
            phase = cre_des_sign(a,i,det(1,k))
            phase = phase*cre_des_sign(b,j,det(2,k))
            c = phase*cd(k)*cd(l)
            rdm2(a,i,b,j) = rdm2(a,i,b,j) + c
            rdm2(b,j,a,i) = rdm2(b,j,a,i) + c
            rdm2ab(a,i,b,j) = rdm2ab(a,i,b,j) + c
            rdm2ba(b,j,a,i) = rdm2ba(b,j,a,i) + c
          end if
        end if
      end do
    end do

  end subroutine compute_rdm2

  integer(kind=4) function n_excitations (str1, str2)

    implicit none
    integer(kind=8), intent(in)  :: str1, str2

    n_excitations = popcnt(ieor(str1, str2))
    n_excitations = ishft(n_excitations,-1)

  end function

  subroutine get_single_excitation (str1, str2, ia)

    implicit none
    integer(kind=4), intent(out) :: ia(2)
    integer(kind=8), intent(in) :: str1, str2

    integer(kind=8) :: tmp, particle, hole

    tmp = ieor(str1, str2)
    particle = iand(tmp, str2)
    hole = iand(tmp, str1)

    ia(2) = trailz(particle) + 1
    ia(1) = trailz(hole) + 1

  end subroutine get_single_excitation

  subroutine get_double_excitation (str1, str2, iajb)

    implicit none
    integer(kind=4), intent(out) :: iajb(4)
    integer(kind=8), intent(in) :: str1, str2

    integer(kind=8) :: tmp, particle, hole

    tmp = ieor(str1, str2)
    particle = iand(tmp, str2)
    hole = iand(tmp, str1)

    iajb(2) = trailz(particle) + 1
    particle = iand(particle, not(ishft(1_8,iajb(2)-1)))
    iajb(4) = trailz(particle) + 1

    iajb(1) = trailz(hole) + 1
    hole = iand(hole, not(ishft(1_8,iajb(1)-1)))
    iajb(3) = trailz(hole) + 1

  end subroutine

  real(kind=8) function cre_des_sign (p,q,str) result(phase)

    implicit none
    integer(kind=4), intent(in) :: p, q
    integer(kind=8), intent(in) :: str

    integer(kind=8) :: mask, pp, qq

    pp = p - 1_8
    qq = q - 1_8
    phase = 0.0d+00

    if (pp > qq) then
      mask = ishft(1_8,pp) - ishft(1_8,qq+1)
      phase = (-1.0d+00)**popcnt(iand(str,mask))
      return
    else !if (pp < qq) then
      mask = ishft(1_8,qq) - ishft(1_8,pp+1)
      phase = (-1.0d+00)**popcnt(iand(str,mask))
      return
    end if

  end function

  subroutine occ_list (str,nelec,occ)

    implicit none
    integer(kind=8), intent(in) :: str
    integer(kind=4), intent(in) :: nelec
    integer(kind=4), intent(out), dimension(nelec) :: occ
  
    integer(kind=4) :: idx, i, i_occ
  
    idx = 1
    do i = 0,63
      i_occ = int(iand(ishft(str,-i),1_8),4)
      if (i_occ.eq.1) then
        occ(idx) = i + 1
        idx = idx + 1
      end if
    end do

  end subroutine

end module mod_rdm
