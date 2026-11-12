c
c----------------------------------------------------------------------
c
      subroutine getdate (line)
c
c----------------------------------------------------------------------
c     Extracts Time and date from the operating system.
c     Presently, this routine probably only works well in Linux.
c----------------------------------------------------------------------
*     character*(*) line
*     integer leng
c-----------------------------------------------------------------------
*     call system ('date > file_with_date')
*     open (unit=99,file='file_with_date')
*     read (99,'(a)') line
*     close (99)
*     call system ('rm file_with_date')
*     line = line(1:leng(line))
*     return
        character(len=*), intent(inout) :: line
        character(len=132) :: sdate
        integer(kind=4) :: values(8)

        call date_and_time(values=values)
        write(sdate,'(I4.4,".",I2.2,".",I2.2,", ",
     & I2.2,":",I2.2,":",I2.2,".",I3.3)') values(1:3),values(5:8)
        line = trim(adjustl(sdate))
      end
