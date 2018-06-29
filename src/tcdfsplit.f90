! This software is distributed under the
! terms of the GNU General Public License
! v3 or any later version.
! Copyright Stefan Gary, 2018.
!----------------------------------------

      program tcdfsplit

!----------------------------------------
! This program will load the trajectory
! data (in cdf format) and split the file
! into many subensembles.  The basic goal
! is to split trajectory launches up from
! a big run containing multiple launches.
!
! As it is now, temporal variability in the
! launches is ignored, so blocks of every
! n trajectories are split off into a ntraj/n
! output files.
!
! More generality can be added later, but I
! don't need it now.  See flags in print_help.
!
!--------Example of invokation----------
!
! tcdfsplit -n 10 -I infile.nc -T -S -H -s -v
!
! will separate infile.nc into however many files
! contain 10 trajectories.
!
!---------------------------------------
!
!----------------------------------------

        use tcdfio
        use load_tcdf
        use netcdfio
        use basicfun

        implicit none

        ! File, dimension, and variable IDs
        ! _i = input file
        ! _n = no (reject file)
        ! _y = yes (accepted file)
        ! _j = indef. term. (with -J opt.)
        integer :: ncid_i, ncid_y, ncid_n, ncid_j, vid

        integer :: num_split_file, split_file_io
        integer, allocatable :: split_by_first_index(:)
        
        ! Trajectory variables
        real, allocatable :: lam(:)
        real, allocatable :: phi(:)
        real, allocatable :: dep(:)
        real, allocatable :: temp(:)
        real, allocatable :: salt(:)
        real, allocatable :: rho(:)
        real, allocatable :: u(:)
        real, allocatable :: v(:)
        real, allocatable :: w(:)
        real, allocatable :: q(:)
        real, allocatable :: dtdz(:)
        real, allocatable :: drdz(:)
        real, allocatable :: year(:)
        real, allocatable :: month(:)
        real, allocatable :: day(:)
        real, allocatable :: age(:)
        real, allocatable :: top(:)
        real, allocatable :: bot(:)
        real, allocatable :: age_in(:)
        real, allocatable :: age_out(:)
        real, allocatable :: bdep(:)
        real, allocatable :: iit(:)
        real, allocatable :: jjt(:)
        
        ! Counters
        ! p = points(time), t = traj
        ! y = num traj in yes file
        ! n = num traj in no file
        ! y_j = num traj in end.nc file
        integer :: t, p, y, n, pp, ppp, y_j

        ! Command line args, File names.
        integer :: num_command_arg
        integer :: arg_len, arg_count
        character(len=50) :: arg_string
        character(len=10) :: arg_int
        character(len=5) :: arg_flag
        character(len=10) :: arg_real

        character(len=6) :: outname_start = 'split_'
        character(len=3) :: outname_end = '.nc'
        integer :: outname_base_num, outname_num
        character(len=14) :: outname
        character(len=5) :: outnname = 'no.nc'
        character(len=6) :: outjname = 'end.nc'
        integer :: nkeep_io

        ! Flags for possible operations.
        logical :: l_num, l_file

        ! Number of floats in each file
        integer :: nsplit

        ! Accumulating number of floats in file
        integer :: current_number_traj_written

        ! Accumulating number of files that have been opened
        integer :: current_number_output_file

        ! Flags for presence of variables on command line.
        logical :: li, lt, ls, lh, lj, lu, lv, lw, lr
        logical :: l_stamp_id, l_verbose

        ! Vector of 1,0 to mark presence of traj in
        ! yes file (1) or in no file (0).
        integer, allocatable :: inyes(:)

        ! Vector to hold list of traj to extract.
        integer, allocatable :: nkeep(:)

        ! Counters of traj in yes and no files
        ! (to double check).
        integer :: yntraj, nntraj, nslow

        ! Spatial and temporal unit conversions
        real :: pi, deg2rad, deg2m
        real, parameter :: radius = 6370.0e03   ! [meters]
        real :: tstep

        !-----------------------------------------
        write(*,*) ' Starting tcdfsplit...'

        ! Initialize operation flag (only two)
        l_num = .false.
        l_file = .false.

        ! Default value is zero
        nsplit = 0
        
        ! Initialize other flags
        lt = .false.
        ls = .false.
        lh = .false.
        lj = .false.
        lu = .false.
        lv = .false.
        lw = .false.
        lr = .false.
        l_verbose = .false.
        l_stamp_id = .false.
        li = .false.

        !------Get command line information------
        ! First argument: operation
        ! Other arguments can vary in order
        num_command_arg = command_argument_count()
       
        if(num_command_arg .lt. 4) then
           write(*,*) ' Error: Too few command line arguments.'
           call print_help()
           stop
        else

           arg_count = 1

           ! Get operation
           call get_command_argument(arg_count, arg_string, arg_len, exitcode)
           call checkexit(exitcode)
           arg_count = arg_count + 1

           ! Check that operation is valid
           if( index(trim(arg_string),'-N') .ne. 0 ) then
              write(*,*) 'Sorting traj start positions...'
              l_num = .true.

              ! Read in the number of floats desired in each file:
              call get_command_argument(arg_count, arg_int, arg_len, exitcode)
              call checkexit(exitcode)
              arg_count = arg_count + 1
              read(arg_int,'(i10)') nsplit
              write(*,*) 'Will split the main file into files with ',nsplit,' floats.'

           elseif( index(trim(arg_string),'-F') .ne. 0 ) then
              write(*,*) 'Using split file to subdivide this ensemble...'
              l_file = .true.

              ! Get the name of the file
              call get_command_argument(arg_count, arg_string, arg_len, exitcode)
              call checkexit(exitcode)
              arg_count = arg_count + 1

              ! Open the file
              split_file_io = 42
              open(unit=split_file_io,file=trim(arg_string),status='old',form='formatted')
              write(*,*) 'Opened ',trim(arg_string)

              ! Count number of lines
              do t=1,9999999
                 read(split_file_io,'(i10)',end=17) 
              enddo
17            num_split_file = t - 1
              write(*,*)' Number of files to split to: ',num_split_file

              ! Clean up and close.
              rewind(split_file_io)
              close(split_file_io)

              ! Allocate space to hold indices for splitting
              allocate(split_by_first_index(num_split_file+1))

              ! Read and store indices.
              open(unit=split_file_io,file=trim(arg_string),status='old',form='formatted')
              do t = 1,num_split_file
                 read(split_file_io,'(i10)') split_by_first_index(t)
                 write(*,*) ' File: ',t,' starts with traj index:',split_by_first_index(t)
              enddo
              rewind(split_file_io)
              close(split_file_io)
              ! Put essentially a dummy ending in the index list
              ! that will never be triggered.
              split_by_first_index(num_split_file+1) = ntraj+1
              
           else

              write(*,*) 'Unrecognized operation.'
              call print_help()
              stop

           endif

           ! Loop through all the other command line flags
           do

              ! Test that we are still reading all the flags.
              if ( arg_count .gt. num_command_arg ) exit

              call get_command_argument(arg_count, arg_flag, arg_len, exitcode)
              call checkexit(exitcode)
              arg_count = arg_count + 1
              
              if ( index(trim(arg_flag),'-T') .ne. 0 ) then
                 ! Copy temperature
                 lt = .true.

              elseif ( index(trim(arg_flag),'-S') .ne. 0 ) then
                 ! Copy salinity
                 ls = .true.

              elseif ( index(trim(arg_flag),'-H') .ne. 0 ) then
                 ! Copy bottom depth
                 lh = .true.

              elseif ( index(trim(arg_flag),'-J') .ne. 0 ) then
                 ! Copy iit,jjt
                 lj = .true.
                 
              elseif ( index(trim(arg_flag),'-s') .ne. 0 ) then
                 ! Stamp each traj with ID
                 l_stamp_id = .true.

              elseif ( index(trim(arg_flag),'-v') .ne. 0 ) then
                 ! Verbose mode
                 l_verbose = .true.

              elseif ( index(trim(arg_flag),'-h') .ne. 0 ) then
                 ! Print help and quit
                 call print_help()
                 stop

              elseif ( index(trim(arg_flag),'-I') .ne. 0 ) then
                 ! Read input file name
                 li = .true.

                 call get_command_argument(arg_count, arg_string, arg_len, exitcode)
                 arg_count = arg_count + 1
                 call checkexit(exitcode)

                 ! Test validity of input file name by opening input file.
                 write(*,*) 'Open input netcdf file...'
                 ncid_i = ncopn(trim(arg_string),ncnowrit,exitcode)

              else
                 write(*,*) 'ERROR: Unrecognized command line flag.'
                 call print_help()
                 stop

              endif

           enddo

        endif   ! Done reading command line options.
        !=======================================================

        !=======================================================
        !                  INITIALIZATIONS
        !=======================================================

        ! Check that we got an input file name.
        if ( .not. li ) then
           write(*,*) 'Error: No input file specified with -I.'
           call print_help()
           stop
        endif

        ! Gather information about the input file.
        if ( l_verbose ) write(*,*) 'Get dimension IDs...'
        timedid = ncdid(ncid_i,timednam,exitcode)
        trajdid = ncdid(ncid_i,trajdnam,exitcode)

        if ( l_verbose ) write(*,*) 'Get dimension sizes...'
        call ncdinq(ncid_i,timedid,dummy,npts,exitcode)
        call ncdinq(ncid_i,trajdid,dummy,ntraj,exitcode)

        write(*,*) 'There are ',npts,' points in each traj.'
        write(*,*) 'There are ',ntraj,' trajectories.'

        ! Check input -N
        if ( (l_num) .and. (nsplit .ge. ntraj) ) then
           write(*,*) 'ERROR: number of traj per file'
           write(*,*) 'is either the same as or bigger'
           write(*,*) 'than the number of input traj.'
           call print_help()
           stop
        endif

        ! Check input -F
        if ( l_file .and. (split_by_first_index(1) .ne. 1) ) then
           write(*,*) 'ERROR: first line of split file not'
           write(*,*) 'equal to 1.'
           call print_help()
           stop
        endif
        do t=2,num_split_file
           ! Check that this number is always more than the
           ! preceeding line and that total number of traj
           ! is not exceeded in the split file.
           if ( split_by_first_index(t) .le. split_by_first_index(t-1) ) then
              write(*,*) 'ERROR: split file must have increasing'
              write(*,*) 'trajectory index values with each line.'
              call print_help()
              stop
           endif
           if ( split_by_first_index(t) .gt. ntraj ) then
              write(*,*) 'ERROR: split file index exceeds maximum'
              write(*,*) 'number of input trajectories.'
              call print_help()
              stop
           endif
        enddo

        if ( l_verbose ) write(*,*) 'Allocating space for one trajectory...'
        allocate(lam(npts))
        allocate(phi(npts))
        allocate(dep(npts))
        if (lu) then
           allocate(u(npts))
        endif
        if (lv) then
           allocate(v(npts))
        endif
        if (lw) then
           allocate(w(npts))
        endif
        if (lt) then
           allocate(temp(npts))
        endif
        if (ls) then
           allocate(salt(npts))
        endif
        if (lr) then
           allocate(rho(npts))
        endif
        allocate(q(npts))
        allocate(dtdz(npts))
        allocate(drdz(npts))
        allocate(age(npts))
        allocate(top(npts))
        allocate(bot(npts))
        allocate(year(npts))
        allocate(month(npts))
        allocate(day(npts))
        if (lh) then
           allocate(bdep(npts))
        endif
        if (lj) then
           allocate(iit(npts))
           allocate(jjt(npts))
        endif
        ! Specify read vectors
        lag_readst2d(1) = 1
        lag_readct2d(1) = npts
        lag_readct2d(2) = 1

        ! Specify write vectors
        lag_writest2d(1) = 1
        lag_writect2d(1) = npts
        lag_writect2d(2) = 1

        ! Initialize counters
        current_number_traj_written = 0
        current_number_output_file = 0
        outname_base_num = 10000

        if ( l_verbose ) write(*,*) 'Splitting: Loop over each trajectory...'

        ! Split input into files set by -N flag.
        do t = 1,ntraj
           
           if ( l_verbose ) write(*,*) 'For traj ',t,' out of ',ntraj
           
           !write(*,*) 'Specify the read vectors...'
           lag_readst2d(2) = t

           !--------------------------------------------------------
           ! Load position information
           if ( l_verbose ) write(*,*) 'Load lam...'
           call get_lag_var(ncid_i,lamvnam,lam)
           if ( l_verbose ) write(*,*) 'Load phi...'
           call get_lag_var(ncid_i,phivnam,phi)
           if ( l_verbose ) write(*,*) 'Load dep...'
           call get_lag_var(ncid_i,depvnam,dep)
              
           ! Load other variables
           if (lt) then
              if ( l_verbose ) write(*,*) 'Load temp...'
              call get_lag_var(ncid_i,tempvnam,temp)
           endif
           if (ls) then
              if ( l_verbose ) write(*,*) 'Load salt...'
              call get_lag_var(ncid_i,saltvnam,salt)
           endif
           if (lh) then
              if ( l_verbose ) write(*,*) 'Load bdep...'
              call get_lag_var(ncid_i,bdepvnam,bdep)
           endif
           if (lj) then
              if ( l_verbose ) write(*,*) 'Load iit...'
              call get_lag_var(ncid_i,iitvnam,iit)
              if ( l_verbose ) write(*,*) 'Load jjt...'
              call get_lag_var(ncid_i,jjtvnam,jjt)
           endif

           !--------------------------------------------------------
           ! Check for whether to open a new file or not.
           if ( (current_number_output_file .eq. 0) .or. &
                (l_num .and. (current_number_traj_written .ge. nsplit)) .or. &
                (l_file .and. (t .eq. split_by_first_index(current_number_output_file+1))) ) then
                
                ! Start of loop, must initialize output file
                ! OR
                ! Done writing to that file, reset the counter
                ! and create a new file and start copying to that file.

                if ( current_number_output_file .eq. 0 ) then
                   ! Do nothing because beginning of loop
                else
                   ! Need to close previous file
                   call ncclos(ncid_y,exitcode)
                endif

                !------Create, define (one of many) output netcdf files------   
                current_number_output_file = current_number_output_file + 1
                outname_num = outname_base_num + current_number_output_file
                write(outname,'(a,i5,a)') outname_start,outname_num,outname_end
                write(*,*) 'Creating output tcdf file ',outname,'...'
                ncid_y = nccre(outname,ncclobber,exitcode)

                if ( l_verbose ) write(*,*) 'Creating dimensions...'
                ! We set the traj dimension to unlimited because
                ! this is the last index in the variable.  Also,
                ! we are more likely to have more trajectories
                ! than numbers of points (at 3 day intervals).
                ! For most case, we let the length of the traj
                ! stay the same (npts).  However, for worms,
                ! traj will be cut to length of 2*nw.
                timedid = ncddef(ncid_y,timednam,npts,exitcode)
                trajdid = ncddef(ncid_y,trajdnam,0,exitcode)

                vdims_out(1) = timedid
                vdims_out(2) = trajdid

                call dup_lag_var(ncid_i,ncid_y,lamvnam)
                if ( l_verbose ) write(*,*) 'Created longitude variable.'
                call dup_lag_var(ncid_i,ncid_y,phivnam)
                if ( l_verbose ) write(*,*) 'Created latitude variable.'
                call dup_lag_var(ncid_i,ncid_y,depvnam)
                if ( l_verbose ) write(*,*) 'Created depth variable.'
                if ( l_stamp_id ) then
                   ! Need to include ID information in output
                   vid = ncvdef(ncid_y,id_vnam,ncint,1,trajdid,exitcode)
                   if ( l_verbose ) write(*,*) 'Created ID variable'
                endif
                if (lt) then
                   call dup_lag_var(ncid_i,ncid_y,tempvnam)
                   if ( l_verbose ) write(*,*) 'Created temperature variable.'
                endif
                if (ls) then
                   call dup_lag_var(ncid_i,ncid_y,saltvnam)
                   if ( l_verbose ) write(*,*) 'Created salinity variable.'
                endif
                if (lh) then
                   call dup_lag_var(ncid_i,ncid_y,bdepvnam)
                   if ( l_verbose ) write(*,*) 'Created bottom depth variable.'
                endif
                if (lj) then
                   call dup_lag_var(ncid_i,ncid_y,iitvnam)
                   if ( l_verbose ) write(*,*) 'Created i-index position variable.'
                   call dup_lag_var(ncid_i,ncid_y,jjtvnam)
                   if ( l_verbose ) write(*,*) 'Created j-index position variable.'
                endif
                
                if ( l_verbose ) write(*,*) 'Done creating file...'
                call ncendf(ncid_y,exitcode)

                ! Reset the counter of trajectories written to
                ! the current file
                current_number_traj_written = 0

             endif   ! Done with opening new files
             
             !--------------------------------------------------------
             if ( l_verbose ) write(*,*) 'Copy input to output file...'
             current_number_traj_written = current_number_traj_written + 1
             lag_writest2d(2) = current_number_traj_written
             call put_lag_var(ncid_y,lamvnam,lam,npts)
             call put_lag_var(ncid_y,phivnam,phi,npts)
             call put_lag_var(ncid_y,depvnam,dep,npts)
             if (l_stamp_id) then
                ! Get stamp variable ID and then write ID
                vid = ncvid(ncid_y,id_vnam,exitcode)
                call ncvpt1(ncid_y,vid,current_number_traj_written,1,exitcode)
             endif
             if (lt) call put_lag_var(ncid_y,tempvnam,temp,npts)
             if (ls) call put_lag_var(ncid_y,saltvnam,salt,npts)
             if (lh) call put_lag_var(ncid_y,bdepvnam,bdep,npts)
             if (lj) then
                call put_lag_var(ncid_y,iitvnam,iit,npts)
                call put_lag_var(ncid_y,jjtvnam,jjt,npts)
             endif
          enddo   ! End of looping over traj points

        !==================================================
        !                    CLEAN UP
        !==================================================
        if ( l_verbose ) write(*,*) 'Close input file...'
        call ncclos(ncid_i,exitcode)

        if ( l_verbose ) write(*,*) 'Close output file...'
        call ncclos(ncid_y,exitcode)

        if ( l_verbose ) write(*,*) 'Deallocate input trajectories...'
        deallocate(lam,phi,dep)
        if ( lt ) deallocate(temp)
        if ( ls ) deallocate(salt)
        if ( lh ) deallocate(bdep)

        write(*,*) 'End of tcdfsplit.'

      end program tcdfsplit

!---------------------------------------
!---------------------------------------

      subroutine print_help()

!---------------------------------------
! This subroutine prints help info
! about this program.
!---------------------------------------

        write(*,*) '-------------------------------------------'
        write(*,*) ' tcdfsplit <-N <N> | -F <splitfile> > -I infile.nc'
        write(*,*) '          [-s] [-v] [-h] [-T] [-S] [-H] [-J]'
        write(*,*) ' '
        write(*,*) '-------------------------------------------'
        write(*,*) 'This program will split a big traj ensemble'
        write(*,*) 'specified by the -I flag into:'
        write(*,*) 'Option 1) many files with N trajectories as'
        write(*,*) '          specified by the -N flag, or'
        write(*,*) 'Option 2) many files with different numbers'
        write(*,*) '          of traj as specified by the given'
        write(*,*) '          splitfile with format:'
        write(*,*) '          index_start_first_file = 1'
        write(*,*) '          index_start_second_file'
        write(*,*) '          ...'
        write(*,*) '          index_start_last_file'
        write(*,*) '          With this format, ALL traj. will'
        write(*,*) '          be written to output since the'
        write(*,*) '          first line must be 1 and the last'
        write(*,*) '          file will assume read to the end.'
        write(*,*) ' THE FIRST FLAG MUST BE THE OPERATION FLAG'
        write(*,*) ' (-N or -F) FOLLOWED BY THE -I FLAG. Other,'
        write(*,*) ' subsequent, flags can go in any order.'
        write(*,*) '-------------------------------------------'
        write(*,*) 'Other optional flags are:'
        write(*,*) ' -v = verbose mode'
        write(*,*) ' -h = prints this message'
        write(*,*) ' -s = stamp each traj with an ID before'
        write(*,*) '      chopping traj up into worms or lines.'
        write(*,*) '      This is useful for going back to look'
        write(*,*) '      at the original trajectories if you find'
        write(*,*) '      something interesting after later sortings.'
        write(*,*) '      ID numbers are written to the output file.'
        write(*,*) ' -T = copy trajectory temperature (temp)'
        write(*,*) ' -S = copy trajectory salinity (salt)'
        write(*,*) ' -H = copy trajectory bottom depth (bdep)'
        write(*,*) ' -J = copy trajectory i and j index locations'
        write(*,*) '      on the original model grid (iit,jjt)'
        write(*,*) '--------------------------------------------------'

        return

      end subroutine print_help

!---------------------------------------
