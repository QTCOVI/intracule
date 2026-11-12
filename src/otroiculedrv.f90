!
!-----------------------------------------------------------------------
!
      subroutine otroiculedrv (totel,udat,uin,uout,uerr,ncent,verb,wfnfile)
      USE space_for_angles
      USE mod_rdm
      USE acumula_intracule
      USE binary128_kernels, only: b128_mode, b128_tol
      implicit none
!
!-----------------------------------------------------------------------
!
!     Obtains INTRACULE-DENSITY PROPERTIES OF A MOLECULE
!
!-----Input variables
!
!     totel         ---> Total number of electrons of the molecule
!     udat          ---> Logical unit of the WFN file
!     uin,uout,uerr ---> logical input, output, and error units
!     ncent         ---> Number of atoms of the molecule
!     verb          ---> If verb=.true. a verbose output is produced
!     wfnfile       ---> Name of the wavefunction file
!
!-----The orders than can be used between the INTRACULE and END keywords
!     are the following:
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! VERBOSE
!
!     ==> Large output is invoked. It can possibly change the value of
!         the input parameter 'verb'
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! INTEGRATE  [  { ANALYTIC | NUMERICAL } ]
!
!     ==> Request the computation of the radial Intracule density either
!         analytically (ANALYTIC) or numerically (NUMERICAL}. The default
!         is ANALYTIC.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! BINARY128 [ AUTO [tol] | ALWAYS | NEVER ]
!
!     ==> IEEE binary128 (real(16)) evaluation of the analytic norm pass
!         (INTEGRATE ANALYTIC: N(N-1), N^2, -N).  AUTO (the default when
!         the keyword is absent) runs each shell quartet in binary64 with
!         a running rounding-error estimate and recomputes in binary128
!         the quartets whose estimated relative error exceeds tol
!         (default 1.0d-10).  BINARY128 alone, or ALWAYS, uses binary128
!         for every quartet; NEVER restores the binary64-only path.
!         The repulsion energies (R^-1 moments) stay in binary64.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! BESSEL
!
!     ==> Request the Bessel-like angular integration mode in computing
!         the spherically averaged (or radial) Intracule density. 
!         The default is .FALSE.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! RMESH rmesh
!
!     ==> Defines the radial quadrature and r-mapping in performing
!     the radial integration of the radial Intracule density.
!
!     rmesh = 1 => The radial grid is rmeshb:
!      - Becke radial grid
!      - r = rmid * (1+q)/(1-q)
!      - the q-mesh is uniform on the interval (-1,+1)
!
!     rmesh = 2 => The radial grid is rmeshh1
!      - Eulec-Maclaurin radial grid with n=1
!      - r = rmid * q/(1-q)
!      - the q-mesh is uniform on the interval (0,+1).
!
!     rmesh = 3 => The radial grid is rmeshh2
!      - Eulec-Maclaurin radial grid with n=2
!      - r = rmid * q^2/(1-q)^2
!      - the q-mesh is uniform on the interval (0,+1).
!
!     The default is rmesh = 2
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! NRAD
!
!     ==> Defines the number of r-points in radial integrations
!     radial Intracule density. The default value is 101.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! RMID
!
!     ==> Defines the value of the 'rmid' parameter in r-mappings
!     (See the RMESH keyword). The default value is 0.5
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! RADIALGRID file_with_r_values [  { ANALYTIC | NUMERICAL } ] 
!
!     ==> This keyword activates the computation of the radial 
!     Intracule density at the r points given in 'file_with_r_values'. 
!     Each line in this file contains a single 'r' value. 
!     The meaning of ANALYTIC and NUMERICAL is the same as the above
!     in the INTEGRATE order.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! ABSQUADERR epsabs
!
!     ==> This keyword redefines the absolute error requested in the
!     QUADPACK integration of the intracule density. The default value
!     for epsabs is 1E-5.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! RELQUADERR epsrel
!
!     ==> This keyword redefines the relative error requested in the
!     QUADPACK integration of the intracule density. The default value
!     for epsrel is 1E-5.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! DERIV { deriv }
!
!     ==> Request the computation of the Intracule density only 
!     (deriv = 0), the Intracule density and gradient (deriv = 1),
!     or the Intracule density, gradient and Laplacian (deriv = 2)
!     in the following order.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! GENGRID file_with_xyz_values
!
!     ==> This keyword activates the computation of the Intracule
!     density at the (x,y,z) points given in the 'file_with_xyz_values'
!     file. If DERIV = 1 the gradient of the Intracule density is also
!     computed. If DERIV = 2, the gradient and Laplacian are also 
!     obtained.
!
! FILEGRID file_with_xyz_values
! NOCUBIC file_with_xyz_values
!
!     ==> These two orders are equivalent to the GENGRID above
!
!     When IDER = 0  (set up using 'DERIV 0' in the input), only the intracule density
!     is computed. When IDER = 1  (set up using 'DERIV 1' in the input), the intracule
!     density and its gradient are computed. Finally, when IDER = 2  (set up using
!     'DERIV 1' in the input), the intracule density, its gradient and its Laplacian
!     are computed.
!
!     These three computations are carried out in the iculegen0, iculegen1, and
!     iculegen2 routines
!
!     In these routines the following steps are followed.
!
!     1) The second order reduced density matrix in the canonical MO
!     basis (2-RDM) is computed in the 'compute_rdm' routine.
!
!     2) The (x,y,z) points read in from the 'filegrid' file are classified
!     such that three (new) different lists of points (one for x, one for y,
!     and one for z) are obtained.
!
!     3) The expansion coefficients of each product of two cartesian
!     primitive Gaussian functions in terms of Hermite Gaussians are
!     computed in the 'computeEcoef' routine.
!
!     4) The intracule density is computed
!
!     The routines iculegen0, iculegen1, and iculegen2 are used when
!     DERIV = 0, 1, and 2, respectively
!
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! CUBIC  [ { MEM | DISK } [ ngrid [ rmax ] ] ]
!
!     ==> Activates the computation of the intracule density on a cubic grid
!         with NGRID points in each direction with the limits  -RMAX <= X,Y,Z <= RMAX
!         If neither MEM nor DISK keywords are used, the DISK version is used.
!         The default values of NGRID and RMAX are 81 and 10.0, respectively.
!         The calculations with MEM  are carried out in the 'icule0grid' routine.
!         The calculations with DISK are carried out in the 'intraditp'  routine.
!
!         If the order 'DERIV 2' is given, the calculation also includes the gradient
!         and the Laplacian and is performed with the MEM version in 'icule2grid'.
!
!         One can change from DISK to MEM or viceversa by using the KEYWORDS:
!
!         MEMORY   ==> MEM = .true.   and  DISK = .false.
!         NODISK   ==> MEM = .true.   and  DISK = .false.
!         DISK     ==> MEM = .false.  and  DISK = .true.
!         NOMEMORY ==> MEM = .false.  and  DISK = .true.
!
! CUBICGRID ...
! NOGENGRID ...
!
!     ==> These two orders aque equivalent to the above CUBIC order
!
!     The routines
!
!       icule0_disk   (DERIV = 0)
!       icule0_mem    (DERIV = 0)
!       icule1_disk   (DERIV = 1)
!       icule1_mem    (DERIV = 1)
!       icule2_disk   (DERIV = 2)
!       icule2_mem    (DERIV = 2)
!
!      are used
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! QUARTETS
!
!     ==> Activates writing of Intracule contributions due to all 
!     quartets of primitive CGF shells
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! GAUSS_JACOBI_POINTS nptsgj
!
!     ==> Defines the number of points to be used in the Gauss-Jacobi 
!     angular integration that is performed in some of the variants  
!     of the radial Intracule density computation. The default 
!     is nptsgj = 40
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! LEBEDEV [ nangular ]
!
!     ==> Activates the selection of Lebedev grid with 'nangular' points
!     in the NUMERICAL angular integration of Intracule density.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! TDESIGN [ nangular ]
!
!     ==> Activates the selection of the spherical T-design grid with 
!     'nangular' points in the NUMERICAL angular integration of 
!     Intracule density.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! SCREENING
!
!     ==> If this keyword is given quartets of primitive Cartesian Gaussian
!     functions with a SMALL coefficient are neglected in some routines,
!     but not all. This is the default.   
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! NOSCREENING
!
!     ==> If this keyword is given in the input file, the above screening 
!     is obviated; i.e. this keywords sets the 'screening' variable
!     to .false.
!
!-----------------------------------------------------------------------
!
      integer (kind=4),   intent (in) :: udat,uin,uout,uerr,ncent
      real    (kind=8),   intent (in) :: totel
      character (len=*),  intent (in) :: wfnfile
      logical,  intent (inout) :: verb

      include 'dmft.inc'
!
!-----Local variables
!
      real(kind=8), parameter :: pi = 3.141592653589793238462643d0
      real(kind=8), parameter :: fourpi = 4.0d0*pi
      real(kind=8), parameter :: third = 1.0d0/3.0d0
      integer(kind=4), parameter :: maxrval = 2000
      character (len=132) word,uppcase,filegrid,line,oldline,qutsfile
      integer(kind=4) ::  lp,lpold,ider,leng,lword,ios,luf,irval,iquts,lwfun
      integer(kind=4) ::  irmesh,nrad,nptsgj
      real   (kind=8) ::  rmid,r,dr,f4r,f42
      real   (kind=8) ::  repee,repxc,repco,repDF,xeee,xeco,xexc,xeDF
      integer(kind=4) ::  ntotel,ntt,ntc,ntxc
      integer(kind=4) ::  i,nangular
      real   (kind=8) ::  rvalue,rval(maxrval) 
      logical ok,setword,setint,setdble,disk,readfile,gengrid,cubicgrid
      logical goon,screening,exfil,radialgrid,noclaude,quartets
      logical claudenoopt,claudeopt1,claudeopt2,bessel
      logical dointeg,lebedev,tdesign,analytic,numerical,interpol
      logical :: usequadpack
      real   (kind=8) ::  rmax
      real   (kind=8) ::  thresh,threshdef
      integer(kind=4) ::  ngrid
      integer(kind=4) ::  kmaxval,kmaxdef

      real(kind=8), allocatable,dimension (:) :: rads,wrads
      real(kind=8), allocatable,dimension (:) :: xicMO
      real(kind=8), allocatable,dimension (:) :: xicou
      real(kind=8), allocatable,dimension (:) :: xicXC
      real(kind=8), allocatable,dimension (:) :: xicDF

      !
      ! Variables used in the DQAGI routine of QUADPACK for radial integration
      !
      real(kind=8) :: bound, result, aberr
      integer(kind=4) :: inf, neval, ier, last, j
      integer(kind=4), parameter :: limit = 100
      real(kind=8) :: epsabs
      real(kind=8) :: epsrel 
      integer(kind=4) :: iwork(limit)
      real(kind=8) :: work(limit*4)
      real(kind=8) :: absquaderr,relquaderr
      integer(kind=4), allocatable,dimension (:) :: iord
      common /quadarea/ lwfun,verbose,usemuller,first_r_value,wfnaux
      real (kind=8) :: coul, xcval, xcdf
      logical verbose
      logical usemuller
      logical first_r_value
      logical gsl, kummer, expkummer
      character (len=132) wfnaux
      external funbess,funbecke
      
      lwfun       = uout
      goon        = .true.
      ider        = 0
      readfile    = .false.
      screening   = .true.
      gengrid     = .true.
      cubicgrid   = .false.
      verb        = .false.
      disk        = .true.
      radialgrid  = .false.
      dointeg     = .false.
      analytic    = .true.
      numerical   = .false.
      interpol    = .false.
      noclaude    = .false.
      claudenoopt = .false.
      bessel      = .true.
      claudeopt1  = .true.
      claudeopt2  = .false.
      quartets    = .false.
      irmesh      = 2
      nrad        = 101
      rmid        = 0.5d+00
      rmax        = 10.0d+00
      ngrid       = 81
      luf         = 1
      kmaxdef     = 500
      threshdef   = 1.0D-14
      usemuller   = .false.
      gsl         = .true.
      kummer      = .false.
      expkummer   = .false.
      epsabs = 1.0d-5
      epsrel = 1.0d-5
      usequadpack = .false.
!
!-----Max value of k in intracule Taylor series
!
      kmaxval = kmaxdef
      thresh  = threshdef
!
!-----Angular integration features
!
      lebedev    = .true.  ! Lebedev angular quadrature is used by default
      tdesign    = .false. ! Spherical T-design angular quadrature is not used by default
      nangular   = 600
      call drvleb (nangular)
      if (allocated(ct  )) deallocate (ct  )
      if (allocated(st  )) deallocate (st  )
      if (allocated(cp  )) deallocate (cp  )
      if (allocated(sp  )) deallocate (sp  )
      if (allocated(stsp)) deallocate (stsp)
      if (allocated(stcp)) deallocate (stcp)
      if (allocated(angw)) deallocate (angw)
      call allocate_space_for_angles (nangular)
      call lebgrid (ct,st,cp,sp,angw,nangular)
      forall (i=1:nangular) stsp(i) = st(i) * sp(i)
      forall (i=1:nangular) stcp(i) = st(i) * cp(i)
!
!-----Driver of the input
!
      nptsgj     = 40    ! Default number points of the Gauss-Jacobi quadrature
      dmftxc = .false.
      idmft  = 2  ! Default XC functional when dmftxc becomes .true.
      do while (goon) 
        read (uin,'(a)',end=110) line
        oldline = trim(line)
        line = uppcase(line)
        lp = 1
        ok = setword(word,line,lp)
        lword = leng(word)
        if (lword.eq.0.or.word(1:1).eq.'#') then
        elseif (word(1:lword).eq.'END') then
          goon=.false.
        elseif (word(:lword).eq.'DMFTXC') then
           usemuller = .true.
        elseif (word(1:lword).eq.'VERBOSE') then
          verb = .true.
        elseif (word(1:lword).eq.'BINARY128') then
          b128_mode = 2
          if (setword(word,line,lp)) then
            if (word.eq.'AUTO') then
              b128_mode = 1
              ok = setdble (b128_tol,line,lp)
            elseif (word.eq.'NEVER' .or. word.eq.'OFF') then
              b128_mode = 0
            elseif (word.eq.'ALWAYS' .or. word.eq.'ON') then
              b128_mode = 2
            else
              stop 'otroiculedrv: BINARY128 [ AUTO [tol] | ALWAYS | NEVER ]'
            endif
          endif
        elseif (word(1:lword).eq.'BESSEL') then
          bessel = .true.
        elseif (word(1:lword).eq.'NOBESSEL') then
          bessel = .false.
        elseif (word(1:lword).eq.'NOQUADPACK') then
          usequadpack = .false.
        elseif (word(1:lword).eq.'QUADPACK') then
          usequadpack = .true.
        elseif (word(1:lword).eq.'GSL') then
          gsl         = .true.
          kummer      = .false.
          expkummer   = .false.
        elseif (word(1:lword).eq.'KUMMER') then
          gsl         = .false.
          kummer      = .true.
          expkummer   = .false.
        elseif (word(1:lword).eq.'EXPKUMMER') then
          gsl         = .false.
          kummer      = .false.
          expkummer   = .true.
        elseif (word(1:lword).eq.'NOCLAUDE') then
          noclaude    = .true.
          claudenoopt = .false.
          claudeopt1  = .false.
          claudeopt2  = .false.
        elseif (word(1:lword).eq.'CLAUDEOPT1') then
          noclaude = .false.
          claudenoopt = .false.
          claudeopt1 = .true.
          claudeopt2 = .false.
        elseif (word(1:lword).eq.'CLAUDEOPT2') then
          noclaude = .false.
          claudenoopt = .false.
          claudeopt1 = .false.
          claudeopt2 = .true.
        elseif (word(1:lword).eq.'CLAUDENOOPT') then
          noclaude = .false.
          claudenoopt = .true.
          claudeopt1 = .false.
          claudeopt2 = .false.
        elseif (word(1:lword).eq.'SCREENING') then
          screening = .true.
        elseif (word(1:lword).eq.'NOSCREENING') then
          screening = .false.
        elseif (word(1:lword).eq.'QUARTETS') then
          quartets = .true.
        elseif (word(1:lword).eq.'MAXK_TAYLOR') then
          ok = setint (kmaxval,line,lp)
          if (ok) kmaxval = abs(kmaxval)
        elseif (word(1:lword).eq.'THRESH_TAYLOR') then
          ok = setdble (thresh,line,lp)
          if (ok) thresh = abs(thresh)
        elseif (word(1:lword).eq.'ABSQUADERR') then
          ok = setdble (absquaderr,line,lp)
          if (ok) epsabs = abs(absquaderr)
        elseif (word(1:lword).eq.'RELQUADERR') then
          ok = setdble (relquaderr,line,lp)
          if (ok) epsrel = abs(relquaderr)
        elseif (word(1:lword).eq.'RMESH') then
          ok = setint(irmesh,line,lp)
        elseif (word(1:lword).eq.'NRAD') then
          ok = setint(nrad,line,lp)
          if (ok) nrad = max(2,min(1000,abs(nrad)))
        elseif (word(1:lword).eq.'RMID') then
          ok = setdble (rmid,line,lp)
          if (ok) rmid = max(0.05d+00,min(10.0d+00,abs(rmid)))
        elseif (word(1:lword).eq.'LEBEDEV') then
          lebedev = .true.
          tdesign = .false.
          ok = setint(nangular,line,lp)
          if (ok) nangular = max(6,min(20000,abs(nangular)))
          call drvleb (nangular)
          if (allocated(ct)) deallocate (ct)
          if (allocated(st)) deallocate (st)
          if (allocated(cp)) deallocate (cp)
          if (allocated(sp)) deallocate (sp)
          if (allocated(stsp)) deallocate (stsp)
          if (allocated(stcp)) deallocate (stcp)
          if (allocated(angw)) deallocate (angw)
          call allocate_space_for_angles (nangular)
          call lebgrid (ct,st,cp,sp,angw,nangular)
          forall (i=1:nangular) stsp(i) = st(i) * sp(i)
          forall (i=1:nangular) stcp(i) = st(i) * cp(i)
        elseif (word(1:lword).eq.'TDESIGN') then
          lebedev = .false.
          tdesign = .true.
          ok = setint(nangular,line,lp)
          if (ok) nangular = max(6,min(20000,abs(nangular)))
          call drvtdsg(nangular)
          if (allocated(ct)) deallocate (ct)
          if (allocated(st)) deallocate (st)
          if (allocated(cp)) deallocate (cp)
          if (allocated(sp)) deallocate (sp)
          if (allocated(stsp)) deallocate (stsp)
          if (allocated(stcp)) deallocate (stcp)
          if (allocated(angw)) deallocate (angw)
          call allocate_space_for_angles (nangular)
          call seltdsg (ct,st,cp,sp,angw,nangular)
          forall (i=1:nangular) stsp(i) = st(i) * sp(i)
          forall (i=1:nangular) stcp(i) = st(i) * cp(i)
        elseif (word(1:lword).eq.'INTEGRATE') then
          dointeg = .true.
          if (setword(word,line,lp)) then
            if (uppcase(word).eq.'ANALYTIC') then
              analytic = .true.
              numerical = .false.
            else if (uppcase(word).eq.'NUMERICAL') then
              numerical = .true.
              analytic = .false.
            else
            endif
          endif
        elseif (word(1:lword).eq.'CUBICGRID' .or. &
                word(1:lword).eq.'CUBIC' .or.     &
                word(1:lword).eq.'NOGENGRID') then
!
!---------A cubic grid is used. Example 'CUBICGRID MEM 51 8.0'
!
          cubicgrid  = .true.
          dointeg    = .false.
          radialgrid = .false.
          gengrid    = .false.
          lpold = lp
          if (setword (word,line,lp)) then
            if (uppcase(word(1:3)).eq.'MEM') then
              disk = .false.
            elseif (uppcase(word(1:4)).eq.'DISK') then
              disk = .true.
            else
              lp = lpold
            endif
          endif
          ok = setint (ngrid,line,lp)
          if (ok) then
            ngrid = abs(ngrid)
            ok = setdble (rmax,line,lp)
            rmax = abs(rmax)
          endif
        elseif (word(1:lword).eq.'MEMORY' .or. &
                word(1:lword).eq.'NODISK') then
          disk = .false.
        elseif (word(1:lword).eq.'DISK' .or. &
                word(1:lword).eq.'NOMEMORY') then
          disk = .true.
        elseif (word(1:lword).eq.'RADIALGRID') then
          radialgrid = .true.
          cubicgrid  = .false.
          dointeg    = .false.
          gengrid    = .false.
          if (setword(filegrid,oldline,lp)) then
            exfil = .true.
            luf = 1
            do while (exfil)
              luf = luf+1
              inquire (unit=luf,opened=exfil)
            enddo
            open (luf,status='unknown',file=trim(filegrid),iostat=ios)
            if (ios.ne.0) then
              stop '# otroiculedrv.f: Error oppening radialgrid file'
            else
              readfile  = .true.
            endif
            if (setword(word,line,lp)) then
              if (uppcase(word).eq.'ANALYTIC') then
                analytic = .true.
                numerical = .false.
              else if (uppcase(word).eq.'NUMERICAL') then
                numerical = .true.
                analytic = .false.
              else
              endif
            endif
          else
            write (uerr,1000) 'RADIALGRID' 
            stop
          endif
        elseif (word(1:lword).eq.'FILEGRID' .or. &
                word(1:lword).eq.'GENGRID'  .or. &
                word(1:lword).eq.'NOCUBIC') then
          radialgrid = .false.
          cubicgrid  = .false.
          dointeg    = .false.
          gengrid    = .true.
          if (setword(filegrid,oldline,lp)) then
            exfil = .true.
            luf = 1
            do while (exfil)
              luf = luf+1
              inquire (unit=luf,opened=exfil)
            enddo
            open (luf,status='unknown',file=trim(filegrid),iostat=ios)
            if (ios.ne.0) then
              stop '# otroiculedrv.f: Error oppening filegrid file'
            else
              readfile  = .false.
              close (luf)
            endif
          else
            write (uerr,1000)
            stop
          endif
        elseif (word(1:lword).eq.'GAUSS_JACOBI_POINTS') then
          ok = setint (nptsgj,line,lp)
          if (ok) nptsgj = abs(nptsgj)
        elseif (word(1:lword).eq.'DERIV') then
          ok = setint (ider,line,lp)
          if (ok.and.(ider.gt.2.or.ider.lt.0)) ider = 0
        endif
      enddo
 110  continue

      if (.not.readfile .and. radialgrid) then
        stop '# otroiculedrv.f: Error oppening filegrid file'
      endif
!-----------------------------------------------------------------------
!
!-----Obtain a radial grid of NRAD points with the grid irmesh
!
!     irmesh = 1 => The radial grid is rmeshb:
!       - Becke radial grid
!       - r = rmid * (1+q)/(1-q)
!       - the q-mesh is uniform on the interval (-1,+1)

!      irmesh = 2 => The radial grid is rmeshh1
!       - Eulec-Maclaurin radial grid with n=1
!       - r = rmid * q/(1-q)
!       - the q-mesh is uniform on the interval (0,+1).

!      irmesh = 3 => The radial grid is rmeshh2
!       - Eulec-Maclaurin radial grid with n=2
!       - r = rmid * q^2/(1-q)^2
!       - the q-mesh is uniform on the interval (0,+1).
!-----------------------------------------------------------------------
!
!
!     RUN the part of the routine associated to cubicgrid = .true.
!
!
!-----------------------------------------------------------------------
!
      if (cubicgrid) then
        if (ider.eq.2) then
          call icule2grid (rmax,ngrid,ider,.false.,.false.,wfnfile)
        else
          if (disk) then
            call intraditp (rmax,ngrid,.false.,.false.,wfnfile)
          else
            call icule0grid (rmax,ngrid,.false.,.false.,wfnfile)
          endif
        endif
        goto 1001
      endif
!
!-----------------------------------------------------------------------
!
!
!     RUN the part of the routine associated to dointeg = .true.
!
!
!-----------------------------------------------------------------------
!
      if (dointeg .or. radialgrid) then
        if (dointeg) then
          if (.not.allocated(rads)) allocate (rads(nrad))
          if (.not.allocated(wrads)) allocate (wrads(nrad))
          write (uout,10)
          if (.not.usequadpack) then
            if (irmesh.eq.1) then
              write (uout,11)
              call rmeshgl (-1.0d0,+1.0d0,rmid,rads,wrads,nrad)
            elseif (irmesh.eq.2) then
              write (uout,12)
              call rmeshh1 (nrad,rmid,rads,wrads)
            else
              write (uout,13)
              call rmeshh2 (nrad,rmid,rads,wrads)
            endif
            write (uout,*) '# Number of radial points is ',nrad
            write (uout,*) '# The value of RMID is       ',rmid
            allocate (xicMO(nrad))
            allocate (xicou(nrad))
            allocate (xicXC(nrad))
            allocate (xicDF(nrad))
!
!-----------IF requested, open file to write quartets contribucions
!
            if (quartets) then
              exfil = .true.
              iquts = 16
              do while (exfil)
                iquts = iquts+1
                inquire (unit=iquts,opened=exfil)
              enddo
              qutsfile = trim(wfnfile)//'.quartets'
              open (iquts,status='unknown',file=qutsfile,iostat=ios)
            endif
          endif
        else
          usequadpack = .false.
          irval = 0
          do
            read (luf,*,iostat=ios) rvalue
            if (ios .ne. 0) exit
            if (irval.eq.maxrval) then
              stop '# otroiculedrv.f: Increase the value of maxrval'
            endif
            irval = irval+1
            rval(irval) = rvalue
          enddo
          if (irval == 0) then
            write (0,'(a)') 'File "'//trim(filegrid)//'" is empty'
            stop
          endif
          nrad = irval
          allocate (rads(nrad))
          allocate (wrads(nrad))
          allocate (xicMO(nrad))
          allocate (xicou(nrad))
          allocate (xicXC(nrad))
          allocate (xicDF(nrad))
          rads(1:nrad) = rval(1:nrad)
          wrads(1:nrad) = 0.0d+00
        endif

        if (usequadpack) then
          call prequadpack ()
          call acumula_init()
          close (udat)
          bound = 0.0d0   ! Límite inferior
          inf = 1         ! 1 = (bound, +infinito)
          wfnaux = trim(wfnfile)
          nquad = 0
          first_r_value = .true.
          write (uout,8) 'Using QUADPACK radial integration'
          if (analytic) then
            write (uout,8) 'Using ANALYTIC (Bessel) angular integration'
            write (uout,9)
            call dqagi (funbess,bound,inf,epsabs,epsrel,result,&
             aberr,neval,ier,limit,limit*4,last,iwork,work)
          else
            write (uout,8) 'Using NUMERICAL (Becke) angular integration'
            write (uout,9)
            call dqagi (funbecke,bound,inf,epsabs,epsrel,result,&
             aberr,neval,ier,limit,limit*4,last,iwork,work)
          endif
          result = result * fourpi
          aberr = aberr * fourpi
          write (uout,29)
          allocate (iord(nquad))
          forall (i=1:nquad) iord(i)=i
          call qcksort (rr, iord, 1, nquad) ! Ordering I(R) value by increasing value of R
          do i = 1, nquad
            j = iord(i)
            coul = valcoul(j)
            xcval = valxc(j)
            xcdf = valxcdf(j)
            write (uout,99) rr(j),valt(j),coul,xcval,xcdf
          end do
          call acumula_free()
          ntotel = nint(totel)
          ntt = ntotel*(ntotel-1)
          write (uout,*)   '# ---------------------------------------'
          write (uout,'(a,6x,F15.8,a,I8,a)') ' # Total Ordered pairs', &
                  result,'    (It should be',ntt,')'
          write (uout,'(a,F15.8)')' # Estimated absolute error ',aberr
          write (uout,'(a,I15)') ' # Number of evaluations    ',neval
          write (uout,'(a,I15)') ' # Error code               ',ier
          write (uout,*) '# ----------------------------------------'
!
!         Analytic (in principle, exact) integradial of radial intracule 
!
          call radial_integ (verbose,wfnfile,gsl,kummer,expkummer)
        else
          write (uout,8) 'Using STANDARD radial integration'
          if (ncent.eq.1) then
            call onecenteric (verb,wfnfile,rads,xicMO,xicou,xicXC,xicDF,nrad,1,thresh)
            goto 1002
          endif
          if (analytic) then
            write (uout,8) 'Using ANALYTIC (Bessel) angular integration'
            call integrate_bessic (verb,quartets,wfnfile,rads,xicMO,xicou,&
               xicXC,xicDF,nrad,iquts,thresh,usemuller,gsl,kummer,expkummer)
          else
            write (uout,8) 'Using NUMERICAL (Becke) angular integration'
            call lebicule_claude1 (verb,wfnfile,rads,xicMO,xicou,xicXC,&
               xicDF,nrad,kmaxval,thresh)
            call radial_integ (verbose,wfnfile,gsl,kummer,expkummer)
          endif
 1002     continue          
        endif
        if (quartets) close (iquts)
!
!-------Integration of radial intracule density
!
        if (dointeg .and. (.not.usequadpack)) then
!
!---------Integration of the intracule densities are here
!
          repee  = 0.0d+00
          repxc  = 0.0d+00
          repco  = 0.0d+00
          repdf  = 0.0d+00
          xeee   = 0.0d+00
          xeco   = 0.0d+00
          xexc   = 0.0d+00
          xeDF   = 0.0d+00
          ntotel = nint(totel)
          ntt = ntotel*(ntotel-1)
          ntc = ntotel*ntotel
          ntxc = ntotel
          do irval=1,nrad
            r = rads(irval)
            dr = r*2.0
            f42 = fourpi * wrads(irval) 
            f4r = f42 / dr
            repee = repee +  xicMO(irval) * f4r
            repxc = repxc +  xicXC(irval) * f4r
            repco = repco +  xicou(irval) * f4r
            repDF = repDF +  xicDF(irval) * f4r
            xeee  = xeee  +  xicMO(irval) * f42 
            xeco  = xeco  +  xicou(irval) * f42 
            xexc  = xexc  -  xicXC(irval) * f42 
            xeDF  = xeDF  -  xicDF(irval) * f42 
          enddo
          write (uout,80)
          write (uout,81) repee,repco,repxc,repDF,xeee,ntt,xeco,ntc,xexc,ntxc,xeDF
          write (uerr,*) '# End of Intracule Density Integration Ok!'
        endif
!
!        Analytic (in principle, exact) integradial of radial intracule 
!
        if (.not.bessel) then
          call radial_integ (verbose,wfnfile,gsl,kummer,expkummer)
        endif
        goto 1001
      endif

      if (gengrid) then
        if (ider.eq.0) then
          call iculegen0 (ider,verb,screening,wfnfile,filegrid)
        elseif (ider.eq.1) then
          call iculegen1 (ider,verb,screening,wfnfile,filegrid)
        elseif (ider.eq.2) then
          call iculegen2 (ider,verb,screening,wfnfile,filegrid)
        else
        endif
      else
        if (ider.eq.0) then
          if (disk) then
            call icule0_disk (ider,verb,screening,wfnfile,filegrid)
          else
            call icule0_mem (ider,verb,screening,wfnfile,filegrid)
          endif
        elseif (ider.eq.1) then
          if (disk) then
            call icule1_disk (ider,verb,screening,wfnfile,filegrid)
          else
            call icule1_mem (ider,verb,screening,wfnfile,filegrid)
          endif
        elseif (ider.eq.2) then
          if (disk) then
            call icule2_disk (ider,verb,screening,wfnfile,filegrid)
          else
            call icule2_mem (ider,verb,screening,wfnfile,filegrid)
          endif
        else
        endif
      endif

 1001 continue
      
      if (allocated(rads)) deallocate (rads)
      if (allocated(wrads)) deallocate (wrads)

      close (udat)
      return

 1000 format(' # otroiculedrv.f: NAME of the file grid not read in',&
         /,' # after the "',a,'" keyword')
 80   format (' # ',61('+'),/,' # Energetic and Electron Population from Intracule',/,' # ',61('+'))
 81   format (&
       ' # Electron-Electron repulsion energy          = ',F15.8,/, &
       ' # Coulomb repulsion energy                    = ',F15.8,/, &
       ' # Exchange-Correlation (XC) repulsion energy  = ',F15.8,/, &
       ' # Dirac-Fock XC repulsion energy              = ',F15.8,/, &
       ' # Total Ordered pairs',24x,' = ',F15.8, ' (it should be ',I4,')',/, &
       ' # Coulomb Ordered Pairs',22x,' = ',F15.8, ' (it should be ',I4,')',/, &
       ' # XC Ordered pairs',27x,' = ',F15.8, ' (it should be ',I4,')',/, &
       ' # Dirac-Fock Ordered Pairs',19x,' = ',F15.8,/,' # ',61('+'))
 9    format (//,' # Q U A D P A C K   I N T E G R A T I O N',/,   &
       ' # Spherical Intracule density',1x,                        &
       '(Tot <- 2RDM, Coul <- 2RDM Coul, XC <- 2RDM(XC) ,',        &
       'XC-DF <- 2RDM(XC-Dirac-Fock))',/,' # ',                    &
       109('-'),/,' # ',6x,'R',14x,'IDensity (Tot)',9x,            &
       'IDensity (Coul)',10x,'IDensity (XC)',10x,'Density (XC-DF)',&
        /,' # ',109('-'),/,' # (First writing --> Not ordered values of R)') 
 10   format (' #',/,' # Radial Intracule Density',/,' # ',41('-'))
 11   format (' # Becke radial grid used',/,' # r = rmid * (1+q)/(1-q)', &
         /,' # q-mesh = uniform on the interval (-1,+1)')
 12   format (' # Eulec-Maclaurin radial grid with n=1 used',/,          &
       ' # r = rmid * q/(1-q)',/,' # q-mesh = uniform on the interval (0,+1)')
 13   format (' # Eulec-Maclaurin radial grid with n=2 used',/,          &
         ' # r = rmid * q^2/(1-q^2)',/,' # q-mesh = uniform on the interval (0,+1)')
 29   format (' # ',109('+'),/,' # (Second writing --> Ordered values of R)',/,' # ',109('+'))
 99   format (1x,F15.8,8(2x,E22.15))
 8    format (/,' # ',a)
      end
