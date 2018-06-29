! This software is distributed under the
! terms of the GNU General Public License
! v3 or any later version.
! Copyright Stefan Gary, 2018.
!----------------------------------------

      program tcdf2xy

!----------------------------------------
! This program will load the trajectory
! data (in cdf format) and write out a
! subset of the trajectory lon,lat positions
! to an ASCII .xy file for plotting.
!
! Comand line information must contain:
! the input file name, the first traj
! number (ti), the last traj number (tf),
! and the skip between ti and tf (ts).
!
! One may also inclue the number of points
! for the output trajectory as an option
! for the last command line option.
!
!--------Example of invokation----------
!
! tcdfsort infile.cdf <ti> <tf> <ts> [<np>]
!
!---------------------------------------

        implicit none

        ! Declare external functions for
        ! NetCDF.
        integer, external :: ncopn
        integer, external :: ncdid
        integer, external :: ncvid

        ! Other parameters
        integer :: ncfloat = 5
        integer :: ncnowrite = 0
        integer :: ncclobber = 0

        ! File, dimension, and variable IDs
        integer :: ncfid, xyfid
        integer :: trajdid, timedid
        integer :: lamvid, phivid, depvid
        integer :: tempvid, saltvid, rhovid
        integer :: uvid, vvid, wvid
        integer :: readstart(2)
        integer :: readcount(2)
        integer :: vdims(2)

        ! Dimension, variable names
        character(len=4) :: trajdnam = 'traj'
        character(len=4) :: timednam = 'time'

        character(len=3) :: lamvnam = 'lam'
        character(len=3) :: phivnam = 'phi'
        character(len=3) :: depvnam = 'dep'
        character(len=4) :: tempvnam = 'temp'
        character(len=4) :: saltvnam = 'salt'
        character(len=3) :: rhovnam = 'rho'
        character(len=1) :: uvnam = 'u'
        character(len=1) :: vvnam = 'v'
        character(len=1) :: wvnam = 'w'
        
        character(len=50) :: dummy

        ! Exitcode
        integer :: exitcode

        ! Declare sizes of variables
        ! npts = initial number of points
        integer :: ntraj, npts

        ! Counters
        ! p = points(time), t = traj
        integer :: t, p, skipcount, trajcount

        ! Trajectory variables
        real, allocatable :: lam(:,:)
        real, allocatable :: phi(:,:)
        real, allocatable :: dep(:,:)
        real, allocatable :: temp(:,:)
        real, allocatable :: salt(:,:)
        real, allocatable :: rho(:,:)
        real, allocatable :: u(:,:)
        real, allocatable :: v(:,:)
        real, allocatable :: w(:,:)

        ! Command line args, File names.
        integer :: num_command_arg
        integer :: arg_len
        character(len=50) :: infname
        character(len=50) :: inti,intf,ints,innp
        integer :: ti,tf,ts,np
        character(len=7) :: xyfname = 'traj.xy'

        !-----------------------------------------
        write(*,*) ' Starting tcdf2xy...'

        !------Get command line information------
        ! First argument: Name of file to restart.
        ! Second argument: Final number of points.
        ! Third argument: Trim mode.
        num_command_arg = command_argument_count()
       
        if(num_command_arg .lt. 4) then
           write(*,*) ' Error: Too few command line arguments.'
           stop

        elseif(num_command_arg .gt. 5) then
           write(*,*) ' Warning: Too many command line arguments.'
           stop

        else

           ! Get input file name.
           call get_command_argument(1, infname, arg_len, exitcode)

           ! Test validity of input file name by opening input file.
           write(*,*) 'Open input netcdf file...'
           ncfid = ncopn(trim(infname),ncnowrite,exitcode)

           write(*,*) 'Get dimension IDs...'
           timedid = ncdid(ncfid,timednam,exitcode)
           trajdid = ncdid(ncfid,trajdnam,exitcode)

           write(*,*) 'Get dimension sizes...'
           call ncdinq(ncfid,timedid,dummy,npts,exitcode)
           call ncdinq(ncfid,trajdid,dummy,ntraj,exitcode)

           write(*,*) 'There are ',npts,' points in each traj.'
           write(*,*) 'There are ',ntraj,' trajectories.'

           write(*,*) 'Get variable IDs...'
           lamvid = ncvid(ncfid,lamvnam,exitcode)
           phivid = ncvid(ncfid,phivnam,exitcode)
           depvid = ncvid(ncfid,depvnam,exitcode)

           write(*,*) 'Allocating space...' 
           allocate(lam(npts,ntraj))
           allocate(phi(npts,ntraj))
           allocate(dep(npts,ntraj))
           lam = 0.0
           phi = 0.0
           dep = 0.0
              
           write(*,*) 'Specify the read/write vectors...'
           readstart(1) = 1
           readstart(2) = 1

           readcount(1) = npts
           readcount(2) = ntraj
           
           write(*,*) 'Load position data...'
           call ncvgt(ncfid,lamvid,readstart,readcount,lam,exitcode)
           call ncvgt(ncfid,phivid,readstart,readcount,phi,exitcode)
           call ncvgt(ncfid,depvid,readstart,readcount,dep,exitcode)

           write(*,*) 'Close input file...'
           call ncclos(ncfid,exitcode)

           write(*,*) 'Get the operation information...'
           call get_command_argument(2, inti, arg_len, exitcode)
           call get_command_argument(3, intf, arg_len, exitcode)
           call get_command_argument(4, ints, arg_len, exitcode)

           write(*,*) 'Convert strings to reals...'
           read(inti,'(i10)') ti
           read(intf,'(i10)') tf
           read(ints,'(i10)') ts

           write(*,*) ' Printing out traj: [',ti,':',ts,':',tf,']'

           if ( num_command_arg .gt. 4 ) then
              ! We have included the optional number of points.

              ! Read final number of points from command line.
              call get_command_argument(5, innp, arg_len, exitcode)
              read(innp,'(i10)') np

              if ( (np .gt. npts) .or. (np .lt. 1) ) then
                 write(*,*) 'Error: Nonsense value for np!'
                 stop
              endif

           else
              ! Otherwise, preserve the length of the traj.
              np = npts

           endif

           write(*,*) ' Cropping traj from ',npts,' to ',np,' points.'

        endif

        !------Done reading from the command line------

        write(*,*) 'Opening traj.xy...'
        xyfid = 42
        open(unit=xyfid,file=xyfname,status='new',form='formatted')

        ! Initialize skipcounter
        skipcount = 0
        trajcount = 1

        ! Loop over each trajectory
        do t = 1,ntraj

           ! Check that we are within the traj # range requested
           if( (t.ge.ti) .and. (t.le.tf) ) then

              if ( (trajcount .eq. 1) .or. (skipcount .ge. ts) ) then
                 
                 write(*,*) 'Writing traj ',t,'...'

                 do p = 1,np
                    write(xyfid,'(f10.6,x,f10.6)') lam(p,t),phi(p,t)
                 enddo

                 write(xyfid,'(a1)')'>'

                 trajcount = trajcount + 1

                 skipcount = 0

              else

                 skipcount = skipcount + 1

              endif

           endif

        enddo
     
        close(xyfid)

        write(*,*) 'Wrote ',trajcount,' trajectories.'

        write(*,*) 'Opening launch.xy...'
        xyfid = 43
        open(unit=xyfid,file='launch.xy',status='new',form='formatted')

        ! Write the start position for each trajectory
        do t = 1,ntraj
           write(xyfid,'(f10.6,x,f20.6)') lam(1,t),dep(1,t)
        enddo

        close(xyfid)

        write(*,*) 'Deallocate input trajectories...'
        deallocate(lam,phi,dep)

        write(*,*) 'End of tcdf2xy.'

      end program tcdf2xy

!----------------------------------------
