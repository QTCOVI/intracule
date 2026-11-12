module xyzgrid
  save
  integer (kind=4), parameter :: maxpts = 501
  integer (kind=4)            :: ngridxyz 
  real    (kind=8), allocatable, dimension(:,:) :: gridp
  integer (kind=4), allocatable, dimension(:,:) :: indgrid
  integer (kind=4), allocatable, dimension(:)   :: npts

contains

  subroutine clasgrid (filegrid,verbose,uout)

  implicit none
  character (len=*), intent(in)    :: filegrid
  integer (kind=4),  intent(in)    :: uout
  logical,           intent(in)    :: verbose

  real    (kind=8)  :: value(3)
  integer (kind=4)  :: lu,ios,i,j,k,n
  logical           :: exfil
  character(len=3)  :: coor(3)

  coor(1) = 'X()'
  coor(2) = 'Y()'
  coor(3) = 'Z()'

  lu = 0
  exfil = .true.
  do while (exfil) 
    lu = lu + 1
    inquire (unit=lu,opened=exfil)
  enddo
  open (lu,status='unknown',file=trim(filegrid))
!
!-Count the number lines of the grid file
!
  n = 0
  do 
    read (lu,*,iostat=ios) (value(i),i=1,3)
    if (ios.ne.0) exit
    n = n + 1
  enddo
  rewind (lu)

  ngridxyz = n
  allocate (npts(3))
  allocate (gridp(maxpts,3))
  allocate (indgrid(ngridxyz,3),stat=ios)
  if (ios.ne.0) then
    stop '# classify_grid.f: Cannot allocate the indgrid() array'
  endif

  if (verbose) then
    write (uout,'(a)') ' # ++++++++++++++++++++++++++++'
    write (uout,'(a)') ' # Input list of (x,y,z) points'
    write (uout,'(a)') ' # ++++++++++++++++++++++++++++'
    write (uout,'(6x,a,17x,a,17x,a)') 'X','Y','Z'
  endif
  npts(:) = 0
  n = 0
  do 
    read (lu,*,iostat=ios) (value(k),k=1,3)
    if (ios.ne.0) exit
    if (verbose) write (uout,'(3(1x,E17.10))') value(1:3)
    
    n = n + 1

    loopk: do k=1,3
      do i=1,npts(k)
        if (value(k).eq.gridp(i,k)) then
          indgrid(n,k) = i
          cycle loopk
        endif
      enddo
      npts(k) = npts(k) + 1
      if (npts(k).gt.maxpts) then
        if (verbose) write (uout,1) coor(k)
        return
      endif
      indgrid(n,k) = npts(k) 
      gridp(npts(k),k) = value(k)
    enddo loopk
  enddo
  
  if (verbose) then
     write (uout,*) '# +++++++++++++++++++++++++++++'
     write (uout,*) '# Classification of grid points'
     write (uout,*) '# +++++++++++++++++++++++++++++'
  endif
  do k=1,3
    if (verbose) write (uout,*) npts(k),coor(k)(1:3)//' points'
    do j=1,npts(k)
      if (verbose) write (uout,'(F17.10)') gridp(j,k)
    enddo
  enddo
 
  if (verbose) then
    write (uout,*) '# +++++++++++++++++++++++++++++++++++++'
    write (uout,*) '# (X,Y,Z) Order Number in the list'
    write (uout,*) '# +++++++++++++++++++++++++++++++++++++'
  endif
  do i=1,n
    if (verbose) write (uout,'(" # (OrdX,OrdY,OrdZ) = ",3I6)') indgrid(i,:) 
  enddo

1 format (' # classifygrid.f: Too many ',a,' points. Increase MAXPTS in clasgrid.f90 routine')

  end subroutine clasgrid
end module xyzgrid
