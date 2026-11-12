      module space_for_slm
c
c     Real spherical harmonic (RSH) S_{l,m} will be stored as
c     element (l*(l+1)+m) of a supervector slm(0:maxl(maxl+2))
c
c     rshnorm() ....... normalization constant of the RSH
c     il0() ........... index of the S_{l,0} RSH (l*(l+1))
c     pil0() .......... pil0(k) = Number of RSH with l < k
c
      save
      real(kind=8),  allocatable, dimension (:)   :: rshnorm,sqrin
      real(kind=8),  allocatable, dimension (:,:) :: deltam
      integer, allocatable, dimension (:)   :: il0,pil0
      integer, allocatable, dimension (:,:) :: jlm,islm
      end module space_for_slm
c-----------------------------------------------------------------------
      subroutine allocate_space_for_slm (maxl)
      USE        space_for_slm
      maxl2=maxl+maxl
      allocate (rshnorm(0:maxl2*(maxl2+2)),stat=ier) 
      if (ier.ne.0) stop 'slmarea.f: Cannot allocate rshnorm()'
      allocate (sqrin(0:2*maxl2+1),stat=ier) 
      if (ier.ne.0) stop 'slmarea.f: Cannot allocate sqrin()'
      allocate (deltam(0:maxl,0:maxl),stat=ier) 
      if (ier.ne.0) stop 'slmarea.f: Cannot allocate deltam()'
      allocate (il0(0:maxl2),stat=ier) 
      if (ier.ne.0) stop 'slmarea.f: Cannot allocate il0()'
      allocate (pil0(0:maxl2),stat=ier)
      if (ier.ne.0) stop 'slmarea.f: Cannot allocate pil0()'
      allocate (jlm(0:maxl*(maxl+2),2),stat=ier) 
      if (ier.ne.0) stop 'slmarea.f: Cannot allocate jlm()'
      allocate (islm(0:maxl*(maxl+2),0:maxl*(maxl+2)),stat=ier) 
      if (ier.ne.0) stop 'slmarea.f: Cannot allocate islm()'
      end subroutine allocate_space_for_slm
c-----------------------------------------------------------------------
      subroutine deallocate_space_for_slm
      USE        space_for_slm
      deallocate (rshnorm,stat=ier) 
      if (ier.ne.0) stop 'slmarea.f: Cannot deallocate rshnorm()'
      deallocate (sqrin,stat=ier) 
      if (ier.ne.0) stop 'slmarea.f: Cannot deallocate sqrin()'
      deallocate (deltam,stat=ier) 
      if (ier.ne.0) stop 'slmarea.f: Cannot deallocate deltam()'
      deallocate (il0,stat=ier) 
      if (ier.ne.0) stop 'slmarea.f: Cannot deallocate il0()'
      deallocate (pil0,stat=ier)
      if (ier.ne.0) stop 'slmarea.f: Cannot deallocate pil0()'
      deallocate (jlm,stat=ier) 
      if (ier.ne.0) stop 'slmarea.f: Cannot deallocate jlm()'
      deallocate (islm,stat=ier) 
      if (ier.ne.0) stop 'slmarea.f: Cannot deallocate islm()'
      end subroutine deallocate_space_for_slm
