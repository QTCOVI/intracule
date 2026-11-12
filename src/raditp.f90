module grids
  save
! integer (kind=4), parameter :: maxpts = 1000
  integer (kind=4)            :: ngridxyz 
  real    (kind=8), allocatable, dimension(:,:) :: gridp
  integer (kind=4), allocatable, dimension(:,:) :: indgrid
  integer (kind=4), allocatable, dimension(:)   :: npts

contains

  subroutine raditp (rmax,ngrid,verbose,uout)

  implicit none
  real    (kind=8),  intent(in)    :: rmax
  integer (kind=4),  intent(in)    :: uout
  integer (kind=4),  intent(in)    :: ngrid
  logical,           intent(in)    :: verbose

  real    (kind=8)  :: delta,x(3)
  integer (kind=4)  :: i,j,k,n
  character(len=3)  :: coor(3)

  coor(1) = 'X()'
  coor(2) = 'Y()'
  coor(3) = 'Z()'

  ngridxyz = ngrid * ngrid * ngrid
  delta = (rmax+rmax)/(ngrid-1)
  allocate (gridp(ngridxyz,3))
  allocate (indgrid(ngridxyz,3))

  n = 0
  x(1) = -rmax
  do i=1,ngrid
    gridp(i,1) = x(1)
    x(2) = -rmax
    do j=1,ngrid
      gridp(j,2) = x(2)
      x(3) = -rmax
      do k=1,ngrid
        gridp(k,3) = x(3)
        n = n + 1
        indgrid(n,1) = i
        indgrid(n,2) = j
        indgrid(n,3) = k
        gridp(n,:) = x(:)
        if (verbose) write (6,'(3(1x,F12.5))') x(:)
        x(3) = x(3) + delta
      enddo
      x(2) = x(2) + delta
    enddo
    x(1) = x(1) + delta
  enddo

  allocate (npts(3))
  npts(1:3) = ngrid

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
  end subroutine raditp
end module grids
