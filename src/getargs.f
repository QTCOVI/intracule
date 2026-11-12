c
c-----------------------------------------------------------------------
c
      subroutine getargs (argc, argv)
c
c.....getargs - get command line arguments and argument count
c     on return, argc will contain the number of arguments and
c     argv will be a character array containing the arguments.
c
      parameter         (mstr=132)
      integer           argc
      character*(*)     argv(*)
      argc=iargc()
      do i=1,argc
         call getarg(i,argv(i))
      enddo
      return
      end
c
c-----------------------------------------------------------------------
c
      integer function leng (string)
c
c.....leng - obtains the leng of the string, assuming that the blancs at
c     the end are dummy.
c
      character*(*)     string
      character*1       blank, null, newline
     
      blank=' '
      null=char(0)
      newline=char(10)
c
      do 10 i = len(string), 1, -1
         if (string(i:i) .ne. blank .and.
     &       string(i:i) .ne. null  .and.
     &       string(i:i) .ne. newline      ) then
            leng = i
            return
         endif
 10   continue
      leng = 0
      return
      end
