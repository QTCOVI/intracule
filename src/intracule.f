C----------------------------------------------------------------------
      program          intracule
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c (c) Evelio Francisco, Jose Antonio Quinonero, and Angel Martin Pendas,
c     University of Oviedo, July (2026)
c
c-----------------------------------------------------------------------------------
c
      USE               space_for_wfnbasis
      USE               space_for_wfncoef
      USE               space_for_primgto
      USE               space_for_atdata
      USE               mod_rdm
      USE               mod_futils, only : remove_stack_limit
      USE               iso_c_binding, only : c_int
      include          'implicit.inc'
      include          'io.inc'
      include          'param.inc'
      include          'wfn.inc'
      include          'corr.inc'
      include          'datatm.inc'
      include          'parallel.inc'
      include          'error.inc'
c
      parameter         (mline=132)
      integer           argc
      character*(mline) argv(2),filein,fileout,filedat,line,word,
     &                  uppcase, oline, wdate
      logical    goon,largwr,setword,verb,test,show_help,wide_help
c
c-----Interfaz para isatty (funcion C de la libc)
c
      interface
         function c_isatty(fd) bind(C, name='isatty') result(r)
            import :: c_int
            integer(c_int), value :: fd
            integer(c_int) :: r
         end function c_isatty
      end interface
      integer(c_int) :: tty0
c
c------------------------------------------------------------------
c
c.....Initialize MPI and remove stack limit
c
      call init_parallel ()
c
c     This remove the stack limit, specially on intel compiler but
c     only in the main thread. This is needed for MPI and OMP
c
      call remove_stack_limit (ier)
      call cpu_time (tiempo1)
!$    tiempo1=omp_get_wtime()
      call timer (0,ipid,'_intracule__________',uout)
      call timer (1,ipid,'_intracule__________',uout)
c
c.....Default is short output in two-electron integration routines.
c
      largwr = .false.
c
c.....Parallel options
c
      ichunk = 1
      nthreads = 1
c
c-----------------------------------------------------------------------
c
c-----Logical input/output units
c
      uin  = stdin
      uout = stdout
      uerr = stderr
      udat = 1
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c-----Command line processing
c
      call getargs(argc,argv)
c
!-----Determinar si hay que mostrar la ayuda y parar:
!     - "-h" o "--help" como primer argumento, o
!     - sin argumentos y stdin es un terminal (no redirigido)
c
      show_help = .false.
      wide_help = .false.
      if (trim(argv(1)).eq.'--HELP') wide_help = .true.

      if (argc.ge.1) then
         if (trim(argv(1)).eq.'-h' .or.
     &       trim(argv(1)).eq.'--help') show_help = .true.
      else if (argc.eq.0) then
         tty0 = c_isatty(0_c_int)
         if (tty0 .ne. 0_c_int) show_help = .true.
      endif

      if (show_help) then
         write (6,'(a)') 'Usage: "intracule [ file.in [ file.out ] ]"'//
     &   ' or "intracule [ < file.in [ > file.out ] ]"'
         write (6,'(a)') ''
         write (6,'(a)') 'Options:'
         write (6,'(a)') '  -h, --help   show this help and exit'
         write (6,'(a)') '  --HELP --> print wide INFO and exit'
         write (6,'(a)') ''
         write (6,'(a)') 'Reads input from stdin, writes to stdout.'
         write (6,'(a)') 
     &    'Use shell redirection: "intracule < water.inp > water.out"'
         stop
      endif
c
      if (wide_help) then
        call usage (6)
        call usage_intracule (6)
        stop
      endif
c
!-----Mas de dos argumentos: error
c
      if (argc.ge.3 .or. argc.lt.0) then
       write (uerr,'(a)')'Usage: "intracule [ file.in [ file.out ] ]"'//
     &   ' or "intracule [ < file.in [ > file.out ] ]"'
       write (6,'(a)') 'Usage: "intracule [ file.in [ file.out ] ]"'//
     &   ' or "intracule [ < file.in [ > file.out ] ]"'
       stop
      endif
c
!-----Modo posicional: ./intracule file.in [file.out]
c
      if (argc.ge.1) then
         filein=argv(1)
         open (uin, file=filein, status='old')
         if (argc.eq.2) then
            fileout=argv(2)
            open (uout,file=fileout,status='unknown')
         endif
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c-----Determine the date of the calculation and print banner.
c
      call getdate (wdate)
      call banner (uout,trim(wdate))
c
c-----First reading of the input file just to known if largwr output
c     will be requested in the 'rdwfn.f' routine
c
      write (uout,'(1x,a)') '#'
      write (uout,'(1x,a)') '# Begin reading the input file of the run'
      write (uout,'(1x,a)') '#'
      iln = 1
      do
        read (uin,'(a)',end=36) line
        write (uout,'(4x,a,I3,2x,a,2x,a)') 'Line ',iln,'-->',trim(line)
        lp=1
        if (setword(word,line,lp)) then
          word = uppcase(word)
          if (trim(word).eq.'VERBOSE'   .or.
     &        trim(word).eq.'BIGOUTPUT' .or.
     &        trim(word).eq.'LARGE' ) then
              largwr = .true.
          endif
        endif
        iln = iln+1
      enddo
 36   write (uout,'(1x,a)') '#'
      write (uout,'(1x,a)') '# End reading the input file of the run'
      write (uout,'(1x,a)') '#'
      rewind (uin)
c
c-----True reading of the input
c
      read (uin,'(a)') line
      lp=1
      if (setword(word,line,lp)) then
        filedat = trim(word)
      else
        stop '# intracule.f: Fatal error. WFN file name not read in'
      endif
c
c-----Open the standard WFN file
c
      open (unit=udat,file=trim(filedat),status='old',iostat=ios)
c
c----------------------------------------------------------------------
c     Read the input data from the wfn file and initialize
c----------------------------------------------------------------------
c
c-----Compute factorials, define the types of the Cartesian Gaussians,etc
c
      call initial()
      call info_parallel ()
c
c-----Read the wavefunction file
c
      call rdwfn (filedat)
      call writexyz (trim(filedat)//".xyz")
      call info (filedat)
      goon = .true.
      verb = .false.
      test = .false.
      do while (goon)
 13      read (uin,'(a)',end=11) line
         if (leng(line).eq.0) goto 13
         oline=line
         line=uppcase(line)
c
c........easy analyze line
c
         lline=len(line)
         lp=1
         do while (line(lp:lp) .eq. ' ')
            lp=lp+1
         enddo
         np=lp
         do while (lp.le.lline .and. line(lp:lp).ne.' ')
            lp=lp+1
         enddo
         lword=lp-np
         word=line(np:lp-1)
         if (word(1:1).eq.'#') then
            continue
         else if (word(:lword).eq.'THREADS') then
           line=line(lp:)
           read (line,*) nthreads
!$         nthreads = iabs(nthreads)
!$         write (uout,'(/,1x,a,i3,/,1x,a,I3)')
!$   &     '# Number of OMP THREADS changed to ', nthreads,
!$   &     '# Maximum Number of OMP THREADS is ', OMP_GET_MAX_THREADS()
!$         call omp_set_num_threads (nthreads)
         elseif (word(:lword).eq.'VERBOSE') then
           verb = .true.
         elseif (word(:lword).eq.'TEST') then
           test = .true.
         elseif (word(:lword).eq.'INTRACULE') then
            oline = oline(lp:)
            mlp = maxlprim
            nshells = sum(ngroup(1:ncent))
            if (test) then
              call newiculedrv (totel,udat,uin,uout,uerr,ncent,mlp,verb,
     &           filedat)
            else
              call otroiculedrv (totel,udat,uin,uout,uerr,ncent,verb,
     &           filedat)
            endif
         endif
      enddo
 11   continue
      call timer (4,ipid,'_intracule__________',-1)
      call timer (6,ipid,'_intracule__________',uout)
      call cpu_time (tiempo2)
!$    tiempo2=omp_get_wtime()
      write (uout,1011) tiempo2-tiempo1
      call end_parallel ()

 1011 format (' # Total elapsed time =',F16.6,' seconds')
      end
