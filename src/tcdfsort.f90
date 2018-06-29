! This software is distributed under the
! terms of the GNU General Public License
! v3 or any later version.
! Copyright Stefan Gary, 2018.
!----------------------------------------

      program tcdfsort

!----------------------------------------
! This program will load the trajectory
! data (in cdf format) and sort the
! trajectories as specified by the
! operation flag and optional sorting
! requirements.
!
! Possible operation flags (w/ info) are:
!
! -start         traj that start in domain
! -through       traj that go through domain
! -fastest <n>   n fastest traj to domain
! -number <file> trajectories with the indicies
!                listed in the single column
!                file specified.
! -lexit         traj that exit the domain,
!                  store everything after the
!                  LAST exit point
! -fexit         traj that exit the domain,
!                  store everything after the
!                  FIRST exit point (even if
!                  the traj goes back into
!                  the domain).
! -wormsin  <n>
! -wormsout <n>
! -linesin  <n>
! -linesout <n>
!
! One of the operation flags must be
! specified.  The operation flag must be
! the first item on the command line
!
! The sorting domain is set by the optional
! variable flags:
!
!         Definition                    Example
!
! -X = zonal limits, in deg.        -X -100.0 -20.0
! -Y = meridional limits, in deg.   -Y  -90.0  90.0
! -Z = depth limits, [m] down       -Z    0.0 500.0
! -U = zonal velocity [m/s]         -U   -1.0   1.0
! -V = meridional velocity [m/s]    -V   -1.0   1.0
! -W = vertical velocity [m/s]      -W   -1.0   1.0
! -T = temperature [oC]             -T   -2.0  42.0
! -S = salinity [PSU]               -S    0.0  40.0
! -R = density (Need to include     -R  36.97 36.99 2.0
!      refp, in standard sigma
!      notation
! -Q = potential vorticity          -Q    0.0   3.0
!      [ (ms)^-1 ]
! -A = vertical gradient of temp    -A    0.0   0.01
!      [ oC/m ]
! -B = vertical gradient of dens    -B    0.0   0.01
!      [ kg/m^4 ]
! -P = probability for keeping      -P    0.1
!      a float in the data set.
! -C = copy date information
! -E = age of float                 -E    0.0 100.0
!
! If a variable/domain flag is followed by "all"
! (.e.g. -Tall) then that variable will be written
! to output but default (i.e. very wide) ranges
! will be applied to the sort.  This is done so that
! even if you don't want to sort by a Lagrangian
! variable, it is still written to the output file.
! If the variable is not on the command line domain
! list, then it will not be included in the output
! file(s).
!
! Floats from the input data set that are
! in the doman as specified by the variable
! flags and the selected operation are written
! to yes.nc and floats that are not in the
! selection are written to no.nc. The optional
! flag:
!
! -N = if present, write a no.nc file
!
! allows the user to specify whether or not
! we need to retain the rejected trajectories.
!
! The -M flag allows the user to mask out all
! traj after they have passed through the sort
! box.  This is useful if the user would like
! to make a point to point plot and remove
! recirculations after first arrival.
!
! The -I flag specifies the name of the input
! file.  This flag is required.
!
! The -h flag prints help
!
! The -J flag allows the user to output
! lines with indefinite terminations.
!
!--------Example of invokation----------
!
! tcdfsort -fastest 10 -I infile.nc -X -100.0 -20.0 -Y 14.5 15.5 -Tall
!
! will choose the 10 fastest traj
! that make it to a box [-100.0, -20.0]
! by [14.5,15.5].  The variables to be written
! to output are temperature as well as
! position.
!
!---------------------------------------
! Version 2:           Loads trajectories
!                      one-by-one instead
!                      of all at once to
!                      optimize memory
!                      use and file writing.
!                      Also, we remove
!                      no.nc from mode -Ft
!
! Version 3:            Implements sorting
!                       by all possible traj
!                       variables, not just
!                       by mode.  This is a
!                       direct complement to
!                       tcdftrim3.
!
! Version 4:            Adds selection by
!                       trajectory index
!                       ( -number operation).
!
! Version 5:            Adds selection by
!                       exiting the domain
!                       ( -exit operation)
!
! Version 6:            Change -exit to two
!                       options: -fexit and
!                       -lexit for the first
!                       and last exits.  Also,
!                       add support for opening
!                       dtdz and drdt vars.
!                       Also added support for
!                       packed variables.  This
!                       version also uses the
!                       _lag_var routines for
!                       modularized writing and
!                       reading of FLAME traj.
!
! Version 7:            Added four new modes of
!                       operation:
!                         wormsin, wormsout,
!                         linesin, linesout
!                       Also, changed the masking
!                       to work with data that can
!                       be generally packed.
!
! Version 8:            Added layer sorting filters
!                       Minimum line length can be
!                       specified for lines in/out.
!                       Added -J option for writing
!                       out the indefinite
!                       terminations of lines.
!
! Version 9:            Restructured the layer options
!                       and added the -end mode/operator.
!                       Also added -c, -f, and -k options
!                       to compute velocity from the
!                       particle positions.
!
!                       Added -l option (for other layers)
!                       Changed -L option to be a required
!                       thickness range instead of a minimum
!                       thickness.
!
! Version 10:           Changed the in/out of box check into
!                       a single function (like layer_check).
!                       This is done so that we can easily
!                       test changes in .gt. or .ge. (etc)
!                       in the in/out box check.  The .ge.
!                       may be adding some minor numerical
!                       noise to the exit/entrance budgets
!                       based on worms.
!
! Version 11:           Added the stamp_ID and keep_ID options
!                       that allow for the assignment of an ID
!                       number to each trajectory before they
!                       are broken up into lines/worms and then
!                       keep the ID all the way through any
!                       sorting steps that specify keep_ID.
!
!                       Added -d option to allow for sorting by
!                       day in month.  This was added so that
!                       a specific date can be searched for
!                       (in particular the temporal discontinuity
!                       in recycled model fields).
!
! Version 12:           Added the bottom depth variable, use
!                       with flag -H
!                       Merged this code to be compiled with
!                       tcdfio.f90.
!
! Lots of flag letters are used, the remaining ones are:
! g, j, o, p, q, r, u, v, w, x, z
!
!---------------WISH LIST----------------
! Full modularization of tcdf programs
! ( *_lag_var subroutines and all the
! default and mask values should be in
! an external file.  Also, they should be
! the identical programs as used by the
! traj calc code and other tcdf programs).
! Then, any changes in tcdf format can be
! implemented universally.  To use tcdf,
! One will need to always convert ARIANE
! output to tcdf format.
!
! More careful implementation of data
! masking/fill values etc.  This needs to
! be listed in the netcdf file as an
! attribute, etc, as well as each mask
! value must work for packed and unpacked
! variables.  Need to experiment with the
! the actual masking values.
!
! Careful checks of the presence or
! absence of variables and then using
! available information to reconstruct
! missing information.  For example, if
! traj data has T and S only and we want
! to sort by density, use the T and S to
! build up density after checking for
! the presence of density in input file.
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
        
        ! npts = initial number of points
        integer :: lines_min_npts

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
        character(len=4) :: arg_int
        character(len=5) :: arg_flag
        character(len=10) :: arg_real

        ! Values of domain limiting properties
        ! Set default values here (very wide ranges).
        real :: x1, x2, y1, y2, z1, z2
        real :: u1, u2, v1, v2, w1, w2
        real :: t1, t2, s1, s2, r1, r2
        real :: q1, q2, a1, a2, b1, b2
        real :: e1, e2, k1, k2, m1, m2
        real :: d1, d2, yr1, yr2, h1, h2
        real :: percent_keep = 1.0
        real, allocatable :: current_random(:)
        integer :: nf, nw, current_id

        ! The number of possible ranges to
        ! check is defined here.
        integer, parameter :: n_limits = 18

        ! The domain limiting properties will be
        ! loaded into the two column vector
        real :: all_min_max(2,n_limits)

        ! ...which is then passed to box_check along
        ! with a single column vector storing all
        ! of the various data at that point.
        real :: all_pt_data(n_limits)

        character(len=50) :: inmode
        character(len=6) :: outyname = 'yes.nc'
        character(len=5) :: outnname = 'no.nc'
        character(len=6) :: outjname = 'end.nc'
        integer :: nkeep_io

        ! Flags for possible operations.
        logical :: lstart, lthrough, lfastest
        logical :: lnumber, lno, lmask, lnone
        logical :: lfexit, llexit, lend
        logical :: lwormsin, lwormsout
        logical :: llinesin, llinesout
        logical :: inbox_old, inbox_now, inbox
        integer :: p_beg, p_end

        ! Flags for presence of variables on command line.
        logical :: lx, ly, lz, lh
        logical :: lu, lv, lw
        logical :: lt, ls, lr, ldump
        logical :: lq, lp, li, lj, lk, lday, lmon, lyear
        logical :: la, lb, lc, le, lage_const, l_equality
        logical :: l_stamp_id, l_keep_id
        logical :: lo, ll, lf, lg, llayer
        logical :: inlayer, lfound
        logical :: ltt, lbb, lll
        logical :: lfore, lback, lcent

        ! Vector of 1,0 to mark presence of traj in
        ! yes file (1) or in no file (0).
        integer, allocatable :: inyes(:)

        ! Vector counting number of steps to reaching
        ! the box we search for.
        integer, allocatable :: n2box(:)

        ! Vector holding traj num of the n fastest
        ! trajectories.
        integer, allocatable :: nbest(:)

        ! Vector to hold list of traj to extract.
        integer, allocatable :: nkeep(:)

        ! Counters of traj in yes and no files
        ! (to double check).
        integer :: yntraj, nntraj, nslow

        ! Layer thickness variables
        real :: min_thick, max_thick
        real :: thick_frac
        real :: frac_tmp

        ! Start location logging variables
        ! Column order is lon, lat, value, ID 
        real, allocatable :: xyzymdvi(:,:)

        ! File ID for trajectory dump information
        integer :: tdo
        integer :: subseg_per_traj

        ! Function for checking layer
        logical, external :: layer_check

        ! Function for checking in/out of box
        logical, external :: box_check

        ! Spatial and temporal unit conversions
        real :: pi, deg2rad, deg2m
        real, parameter :: radius = 6370.0e03   ! [meters]
        real :: tstep

        !-----------------------------------------
        write(*,*) ' Starting trajsort...'

        ! Initialize operation flags
        lstart = .false.
        lend = .false.
        lthrough = .false.
        lfastest = .false.
        lnumber = .false.
        lno = .false.
        lnone = .false.
        lmask = .false.
        lfexit = .false.
        llexit = .false.
        lwormsin = .false.
        lwormsout = .false.
        llinesin = .false.
        llinesout = .false.
        inbox = .false.
        lcent = .false.
        lfore = .false.
        lback = .false.

        ! Initialize domain flags
        lx = .false.
        ly = .false.
        lz = .false.
        lu = .false.
        lv = .false.
        lw = .false.
        lt = .false.
        ls = .false.
        lr = .false.
        lq = .false.
        li = .false.
        lj = .false.
        lk = .false.
        lday = .false.
        lmon = .false.
        lyear = .false.
        lp = .false.
        la = .false.
        lb = .false.
        lc = .false.
        le = .false.
        lo = .false.
        ll = .false.
        lf = .false.
        lg = .false.
        lh = .false.
        ltt = .false.
        lbb = .false.
        lll = .false.
        lage_const = .false.
        llayer = .false.
        ldump = .false.
        l_equality = .true.
        l_stamp_id = .false.
        l_keep_id = .false.

        ! Initialize domain size
        x1 = lam_min
        x2 = lam_max
        y1 = phi_min
        y2 = phi_max
        z1 = dep_min
        z2 = dep_max
        u1 = u_min
        u2 = u_max
        v1 = v_min
        v2 = v_max
        w1 = w_min
        w2 = w_max
        t1 = t_min
        t2 = t_max
        s1 = s_min
        s2 = s_max
        r1 = r_min
        r2 = r_max
        q1 = q_min
        q2 = q_max
        a1 = dtdz_min
        a2 = dtdz_max
        b1 = drdz_min
        b2 = drdz_max
        e1 = age_min
        e2 = age_max
        k1 = rld_min
        k2 = rld_max
        h1 = bdep_min
        h2 = bdep_max
        d1 = 1.0
        d2 = 31.0
        m1 = 1.0
        m2 = 12.0
        yr1 = 0.0
        yr2 = 5000.0
        lines_min_npts = 0
        min_thick = 0.0
        max_thick = 10000.0
        thick_frac = 0.5

        ! Start random seed sequence and set
        ! size for one random number
        call init_random_seed()
        allocate(current_random(1))

        !------Get command line information------
        ! First argument: operation
        ! Other arguments can vary in order
        num_command_arg = command_argument_count()
       
        if(num_command_arg .lt. 3) then
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
           if( index(trim(arg_string),'-start') .ne. 0 ) then
              write(*,*) 'Sorting traj start positions...'
              lstart = .true.

           elseif( index(trim(arg_string),'-end') .ne. 0 ) then
              write(*,*) 'Sorting traj end positions...'
              lend = .true.

           elseif ( index(trim(arg_string),'-through') .ne. 0 ) then
              write(*,*) 'Sorting for traj that pass through a box...'
              lthrough = .true.

           elseif ( index(trim(arg_string),'-fexit') .ne. 0 ) then
              write(*,*) 'Sorting for traj at first exit of box...'
              lfexit = .true.

              if ( llexit ) then
                 write(*,*) ' ERROR: Cannot run with both -fexit and -lexit.'
                 stop
              endif

           elseif ( index(trim(arg_string),'-lexit') .ne. 0 ) then
              write(*,*) 'Sorting for traj at last exit of box...'
              llexit = .true.

              if ( lfexit ) then
                 write(*,*) ' ERROR: Cannot run with both -fexit and -lexit.'
                 stop
              endif

           elseif ( index(trim(arg_string),'-fastest') .ne. 0 ) then
              write(*,*) 'Sorting for the fastest traj to a box...'
              lfastest = .true.

              ! Read in the number of floats desired:
              call get_command_argument(arg_count, arg_int, arg_len, exitcode)
              call checkexit(exitcode)
              arg_count = arg_count + 1
              read(arg_int,'(i4)') nf
              write(*,*) 'Will find the ',nf,' fastest floats.'

           elseif ( index(trim(arg_string),'-wormsin') .ne. 0 ) then
              write(*,*) 'Sorting for worms going into box...'
              lwormsin = .true.

              ! Read in the desired worm half length:
              call get_command_argument(arg_count, arg_int, arg_len, exitcode)
              call checkexit(exitcode)
              arg_count = arg_count + 1
              read(arg_int,'(i4)') nw
              write(*,*) 'Will find worms of half length ',nw,'.'

           elseif ( index(trim(arg_string),'-wormsout') .ne. 0 ) then
              write(*,*) 'Sorting for worms going out of box...'
              lwormsout = .true.

              ! Read in the desired worm half length:
              call get_command_argument(arg_count, arg_int, arg_len, exitcode)
              call checkexit(exitcode)
              arg_count = arg_count + 1
              read(arg_int,'(i4)') nw
              write(*,*) 'Will find worms of half length ',nw,'.'

           elseif ( index(trim(arg_string),'-linesin') .ne. 0 ) then
              write(*,*) 'Sorting for lines inside the box...'
              llinesin = .true.

              ! Read in the desired min length:
              call get_command_argument(arg_count, arg_int, arg_len, exitcode)
              call checkexit(exitcode)
              arg_count = arg_count + 1
              read(arg_int,'(i4)') lines_min_npts
              write(*,*) 'Will find lines with minimum length ',lines_min_npts,'.'

           elseif ( index(trim(arg_string),'-linesout') .ne. 0 ) then
              write(*,*) 'Sorting for lines outside the box...'
              llinesout = .true.

              ! Read in the desired min length:
              call get_command_argument(arg_count, arg_int, arg_len, exitcode)
              call checkexit(exitcode)
              arg_count = arg_count + 1
              read(arg_int,'(i4)') lines_min_npts
              write(*,*) 'Will find lines with minimum length ',lines_min_npts,'.'

           elseif ( index(trim(arg_string),'-number') .ne. 0 ) then
              write(*,*) 'Sorting for traj with index numbers...'
              lnumber = .true.

              ! Get the name of the file
              call get_command_argument(arg_count, arg_string, arg_len, exitcode)
              call checkexit(exitcode)
              arg_count = arg_count + 1

              ! Open the file
              nkeep_io = 42
              open(unit=nkeep_io,file=trim(arg_string),status='old',form='formatted')
              write(*,*) 'Opened ',trim(arg_string)

              ! Count number of lines
              do t=1,9999999
                 read(nkeep_io,'(i10)',end=17) 
              enddo
17            nf = t - 1
              write(*,*)' Number of floats to extract: ',nf

              ! Clean up and close.
              rewind(nkeep_io)
              close(nkeep_io)

              ! Allocate space to hold indices
              allocate(nkeep(nf))

              ! Read and store indices.
              open(unit=nkeep_io,file=trim(arg_string),status='old',form='formatted')
              do t = 1,nf
                 read(nkeep_io,'(i10)') nkeep(t)
                 write(*,*) nkeep(t)
              enddo
              rewind(nkeep_io)
              close(nkeep_io)

           else
              write(*,*) 'Unrecognized operation.'
              call print_help
              stop
           endif

           ! Check input.
           if ( lwormsin ) then
              if ( lwormsout ) then
                 write(*,*) ' Cannot run with both worms in and out!'
                 stop
              endif
           endif

           if ( llinesin ) then
              if ( llinesout ) then
                 write(*,*) ' Cannot run with both lines in and out!'
                 stop
              endif
           endif

           if ( lwormsout ) then
              if ( lwormsin ) then
                 write(*,*) ' Cannot run with both worms in and out!'
                 stop
              endif
           endif

           if ( llinesout ) then
              if ( llinesin ) then
                 write(*,*) ' Cannot run with both lines in and out!'
                 stop
              endif
           endif

           ! Loop through all the other command line flags
           do

              ! Test that we are still reading all the flags.
              if ( arg_count .gt. num_command_arg ) exit

              call get_command_argument(arg_count, arg_flag, arg_len, exitcode)
              call checkexit(exitcode)
              arg_count = arg_count + 1
              
              if ( index(trim(arg_flag),'-X') .ne. 0 ) then
                 ! Read in x limits on domain
                 lx = .true.

                 if (index(trim(arg_flag),'all') .ne. 0 ) then
                    ! We do not want to modify range of this variable,
                    ! but we do want to write it to output.
                 else
                    ! Get max and min values.
                    call get_command_argument(arg_count, arg_real, arg_len, exitcode)
                    arg_count = arg_count + 1
                    call checkexit(exitcode)
                    read(arg_real,'(f10.4)') x1

                    call get_command_argument(arg_count, arg_real, arg_len, exitcode)
                    arg_count = arg_count + 1
                    call checkexit(exitcode)
                    read(arg_real,'(f10.4)') x2

                    call check_limits(x1,x2)

                    write(*,'(a10,a2,a7,f10.4,x,f10.4)') 'Property: ',trim(arg_flag),' spans ',x1,x2
                 endif

              elseif ( index(trim(arg_flag),'-Y') .ne. 0 ) then
                 ! Read in y limits on domain
                 ly = .true.

                 if (index(trim(arg_flag),'all') .ne. 0 ) then
                    ! We do not want to modify range of this variable,
                    ! but we do want to write it to output.
                 else
                    ! Get max and min values.
                    call get_command_argument(arg_count, arg_real, arg_len, exitcode)
                    arg_count = arg_count + 1
                    call checkexit(exitcode)
                    read(arg_real,'(f10.4)') y1
                    
                    call get_command_argument(arg_count, arg_real, arg_len, exitcode)
                    arg_count = arg_count + 1
                    call checkexit(exitcode)
                    read(arg_real,'(f10.4)') y2
                 
                    call check_limits(y1,y2)

                    write(*,'(a10,a2,a7,f10.4,x,f10.4)') 'Property: ',trim(arg_flag),' spans ',y1,y2
                 endif

              elseif ( index(trim(arg_flag),'-Z') .ne. 0 ) then
                 ! Read in z limits on domain
                 lz = .true.

                 if (index(trim(arg_flag),'all') .ne. 0 ) then
                    ! We do not want to modify range of this variable,
                    ! but we do want to write it to output.
                 else
                    ! Get max and min values.
                    call get_command_argument(arg_count, arg_real, arg_len, exitcode)
                    arg_count = arg_count + 1
                    call checkexit(exitcode)
                    read(arg_real,'(f10.4)') z1
                 
                    call get_command_argument(arg_count, arg_real, arg_len, exitcode)
                    arg_count = arg_count + 1
                    call checkexit(exitcode)
                    read(arg_real,'(f10.4)') z2
                 
                    call check_limits(z1,z2)

                    write(*,'(a10,a2,a7,f10.4,x,f10.4)') 'Property: ',trim(arg_flag),' spans ',z1,z2
                 endif

              elseif ( index(trim(arg_flag),'-U') .ne. 0 ) then
                 ! Read in u limits on domain
                 lu = .true.
                 
                 if (index(trim(arg_flag),'all') .ne. 0 ) then
                    ! We do not want to modify range of this variable,
                    ! but we do want to write it to output.
                 else
                    ! Get max and min values.
                    call get_command_argument(arg_count, arg_real, arg_len, exitcode)
                    arg_count = arg_count + 1
                    call checkexit(exitcode)
                    read(arg_real,'(f10.4)') u1
                 
                    call get_command_argument(arg_count, arg_real, arg_len, exitcode)
                    arg_count = arg_count + 1
                    call checkexit(exitcode)
                    read(arg_real,'(f10.4)') u2
                 
                    call check_limits(u1,u2)

                    write(*,'(a10,a2,a7,f10.4,x,f10.4)') 'Property: ',trim(arg_flag),' spans ',u1,u2
                 endif

              elseif ( index(trim(arg_flag),'-V') .ne. 0 ) then
                 ! Read in v limits on domain
                 lv = .true.

                 if (index(trim(arg_flag),'all') .ne. 0 ) then
                    ! We do not want to modify range of this variable,
                    ! but we do want to write it to output.
                 else
                    ! Get max and min values.
                    call get_command_argument(arg_count, arg_real, arg_len, exitcode)
                    arg_count = arg_count + 1
                    call checkexit(exitcode)
                    read(arg_real,'(f10.4)') v1
                 
                    call get_command_argument(arg_count, arg_real, arg_len, exitcode)
                    arg_count = arg_count + 1
                    call checkexit(exitcode)
                    read(arg_real,'(f10.4)') v2
                 
                    call check_limits(v1,v2)

                    write(*,'(a10,a2,a7,f10.4,x,f10.4)') 'Property: ',trim(arg_flag),' spans ',v1,v2
                 endif

              elseif ( index(trim(arg_flag),'-W') .ne. 0 ) then
                 ! Read in w limits on domain
                 lw = .true.
                 
                 if (index(trim(arg_flag),'all') .ne. 0 ) then
                    ! We do not want to modify range of this variable,
                    ! but we do want to write it to output.
                 else
                    ! Get max and min values.
                    call get_command_argument(arg_count, arg_real, arg_len, exitcode)
                    arg_count = arg_count + 1
                    call checkexit(exitcode)
                    read(arg_real,'(f10.4)') w1
                 
                    call get_command_argument(arg_count, arg_real, arg_len, exitcode)
                    arg_count = arg_count + 1
                    call checkexit(exitcode)
                    read(arg_real,'(f10.4)') w2
                 
                    call check_limits(w1,w2)

                    write(*,'(a10,a2,a7,f10.4,x,f10.4)') 'Property: ',trim(arg_flag),' spans ',w1,w2
                 endif

              elseif ( index(trim(arg_flag),'-T') .ne. 0 ) then
                 ! Read in t limits on domain
                 lt = .true.
                 
                 if (index(trim(arg_flag),'all') .ne. 0 ) then
                    ! We do not want to modify range of this variable,
                    ! but we do want to write it to output.
                 else
                    ! Get max and min values.
                    call get_command_argument(arg_count, arg_real, arg_len, exitcode)
                    arg_count = arg_count + 1
                    call checkexit(exitcode)
                    read(arg_real,'(f10.4)') t1
                 
                    call get_command_argument(arg_count, arg_real, arg_len, exitcode)
                    arg_count = arg_count + 1
                    call checkexit(exitcode)
                    read(arg_real,'(f10.4)') t2
                 
                    call check_limits(t1,t2)

                    write(*,'(a10,a2,a7,f10.4,x,f10.4)') 'Property: ',trim(arg_flag),' spans ',t1,t2
                 endif

              elseif ( index(trim(arg_flag),'-S') .ne. 0 ) then
                 ! Read in s limits on domain
                 ls = .true.
                 
                 if (index(trim(arg_flag),'all') .ne. 0 ) then
                    ! We do not want to modify range of this variable,
                    ! but we do want to write it to output.
                 else
                    ! Get max and min values.
                    call get_command_argument(arg_count, arg_real, arg_len, exitcode)
                    arg_count = arg_count + 1
                    call checkexit(exitcode)
                    read(arg_real,'(f10.4)') s1
                 
                    call get_command_argument(arg_count, arg_real, arg_len, exitcode)
                    arg_count = arg_count + 1
                    call checkexit(exitcode)
                    read(arg_real,'(f10.4)') s2
                 
                    call check_limits(s1,s2)

                    write(*,'(a10,a2,a7,f10.4,x,f10.4)') 'Property: ',trim(arg_flag),' spans ',s1,s2
                 endif

              elseif ( index(trim(arg_flag),'-R') .ne. 0 ) then
                 ! Read in d limits on domain
                 lr = .true.
                 
                 if (index(trim(arg_flag),'all') .ne. 0 ) then
                    ! We do not want to modify range of this variable,
                    ! but we do want to write it to output.
                 else
                    ! Get max and min values.
                    call get_command_argument(arg_count, arg_real, arg_len, exitcode)
                    arg_count = arg_count + 1
                    call checkexit(exitcode)
                    read(arg_real,'(f10.4)') r1
                 
                    call get_command_argument(arg_count, arg_real, arg_len, exitcode)
                    arg_count = arg_count + 1
                    call checkexit(exitcode)
                    read(arg_real,'(f10.4)') r2

                    call check_limits(r1,r2)

                    write(*,'(a10,a2,a7,f10.4,x,f10.4,x,f5.2)') 'Property: ',trim(arg_flag),' spans ',r1,r2
                 endif

              elseif ( index(trim(arg_flag),'-Q') .ne. 0 ) then
                 ! Read in q limits on domain
                 lq = .true.
                 
                 if (index(trim(arg_flag),'all') .ne. 0 ) then
                    ! We do not want to modify range of this variable,
                    ! but we do want to write it to output.
                 else
                    ! Get max and min values.
                    call get_command_argument(arg_count, arg_real, arg_len, exitcode)
                    arg_count = arg_count + 1
                    call checkexit(exitcode)
                    read(arg_real,'(f10.4)') q1
                 
                    call get_command_argument(arg_count, arg_real, arg_len, exitcode)
                    arg_count = arg_count + 1
                    call checkexit(exitcode)
                    read(arg_real,'(f10.4)') q2

                    call check_limits(q1,q2)

                    write(*,'(a10,a2,a7,f10.4,x,f10.4)') 'Property: ',trim(arg_flag),' spans ',q1,q2
                 endif

              elseif ( index(trim(arg_flag),'-A') .ne. 0 ) then
                 ! Read in dtdz limits on domain
                 la = .true.
                 
                 if (index(trim(arg_flag),'all') .ne. 0 ) then
                    ! We do not want to modify range of this variable,
                    ! but we do want to write it to output.
                 else
                    ! Get max and min values.
                    call get_command_argument(arg_count, arg_real, arg_len, exitcode)
                    arg_count = arg_count + 1
                    call checkexit(exitcode)
                    read(arg_real,'(f10.4)') a1
                 
                    call get_command_argument(arg_count, arg_real, arg_len, exitcode)
                    arg_count = arg_count + 1
                    call checkexit(exitcode)
                    read(arg_real,'(f10.4)') a2

                    call check_limits(a1,a2)

                    write(*,'(a10,a2,a7,f10.4,x,f10.4)') 'Property: ',trim(arg_flag),' spans ',a1,a2
                 endif

              elseif ( index(trim(arg_flag),'-B') .ne. 0 ) then
                 ! Read in drdz limits on domain
                 lb = .true.
                 
                 if (index(trim(arg_flag),'all') .ne. 0 ) then
                    ! We do not want to modify range of this variable,
                    ! but we do want to write it to output.
                 else
                    ! Get max and min values.
                    call get_command_argument(arg_count, arg_real, arg_len, exitcode)
                    arg_count = arg_count + 1
                    call checkexit(exitcode)
                    read(arg_real,'(f10.4)') b1
                 
                    call get_command_argument(arg_count, arg_real, arg_len, exitcode)
                    arg_count = arg_count + 1
                    call checkexit(exitcode)
                    read(arg_real,'(f10.4)') b2

                    call check_limits(b1,b2)

                    write(*,'(a10,a2,a7,f10.4,x,f10.4)') 'Property: ',trim(arg_flag),' spans ',b1,b2
                 endif

              elseif ( index(trim(arg_flag),'-E') .ne. 0 ) then
                 ! Read in age limits on domain
                 le = .true.
                 
                 if (index(trim(arg_flag),'all') .ne. 0 ) then
                    ! We do not want to modify range of this variable,
                    ! but we do want to write it to output.
                 else
                    ! Get max and min values.
                    call get_command_argument(arg_count, arg_real, arg_len, exitcode)
                    arg_count = arg_count + 1
                    call checkexit(exitcode)
                    read(arg_real,'(f10.4)') e1
                 
                    call get_command_argument(arg_count, arg_real, arg_len, exitcode)
                    arg_count = arg_count + 1
                    call checkexit(exitcode)
                    read(arg_real,'(f10.4)') e2

                    call check_limits(e1,e2)

                    write(*,'(a10,a2,a7,f10.4,x,f10.4)') 'Property: ',trim(arg_flag),' spans ',e1,e2
                 endif

              elseif ( index(trim(arg_flag),'-C') .ne. 0 ) then
                 ! Copy date info at end
                 lc = .true.

              elseif ( index(trim(arg_flag),'-D') .ne. 0 ) then
                 ! Dump trajectory chopping info
                 if ( lwormsin .or. lwormsout .or. &
                      llinesin .or. llinesout ) then
                    ldump = .true.

                    tdo = 44
                    open(unit=tdo,file='chop_traj.dump',status='new',form='formatted')
                 else
                    write(*,*) 'WARNING: -D option ignored b/c'
                    write(*,*) 'we are not in lines or worms mode.'
                 endif
              elseif ( index(trim(arg_flag),'-e') .ne. 0 ) then
                 ! Construct age information from traj len.
                 write(*,*) 'Will construct age from trajectory length.'
                 lage_const = .true.

              elseif ( index(trim(arg_flag),'-a') .ne. 0 ) then
                 ! Compute ranges based on non-equality.
                 write(*,*) 'Will compute in or out with < and >.'
                 l_equality = .false.

              elseif ( index(trim(arg_flag),'-s') .ne. 0 ) then
                 ! Stamp each traj (and all lines/worms cut from
                 ! this traj) with an ID
                 write(*,*) 'Will stamp ID on worms or lines.'
                 l_stamp_id = .true.

              elseif ( index(trim(arg_flag),'-i') .ne. 0 ) then
                 ! Carry over the ID for each traj to output.
                 write(*,*) 'Will keep IDs for input file.'
                 l_keep_id = .true.

              elseif ( index(trim(arg_flag),'-P') .ne. 0 ) then
                 ! Read in percent_keep
                 lp = .true.

                 call get_command_argument(arg_count, arg_real, arg_len, exitcode)
                 arg_count = arg_count + 1
                 call checkexit(exitcode)
                 read(arg_real,'(f10.4)') percent_keep
                 
                 write(*,'(a10,a2,a4,f10.4)') 'Property: ',trim(arg_flag),' is ',percent_keep
                 
              elseif ( index(trim(arg_flag),'-I') .ne. 0 ) then
                 ! Read input file name
                 li = .true.

                 call get_command_argument(arg_count, arg_string, arg_len, exitcode)
                 arg_count = arg_count + 1
                 call checkexit(exitcode)

                 ! Test validity of input file name by opening input file.
                 write(*,*) 'Open input netcdf file...'
                 ncid_i = ncopn(trim(arg_string),ncnowrit,exitcode)

              elseif ( index(trim(arg_flag),'-h') .ne. 0 ) then
                 ! Print help information and exit
                 call print_help()
                 stop

              elseif ( index(trim(arg_flag),'-N') .ne. 0 ) then
                 ! We do not write no.nc
                 lno = .true.

              elseif ( index(trim(arg_flag),'-n') .ne. 0 ) then
                 ! We do not write anything to netcdf,
                 ! only output the totals in and/or out.
                 lnone = .true.

              elseif ( index(trim(arg_flag),'-M') .ne. 0 ) then
                 ! We do mask out traj after arriving at the box.
                 lmask = .true.

              elseif ( index(trim(arg_flag),'-t') .ne. 0 ) then
                 ! Require floats to be below top of layer
                 ltt = .true.
                 llayer = .true.

              elseif ( index(trim(arg_flag),'-b') .ne. 0 ) then
                 ! Require floats to be above bottom of layer
                 lbb = .true.
                 llayer = .true.

              elseif ( index(trim(arg_flag),'-l') .ne. 0 ) then
                 ! Change layer sorting to look at
                 ! top2 and bot2 instead of top1 and
                 ! bot1 (default).
                 lll = .true.
                 llayer = .true.

                 layer_top_vnam = 'top2' 
                 layer_bot_vnam = 'bot2'

              elseif ( index(trim(arg_flag),'-O') .ne. 0 ) then
                 ! Require floats to be outcropped
                 lo = .true.
                 llayer = .true.
              
              elseif ( index(trim(arg_flag),'-L') .ne. 0 ) then
                 llayer = .true.

                 if (index(trim(arg_flag),'all') .ne. 0 ) then
                    ! We simply want to copy the layer information,
                    ! not apply it to sorting.
                 else
                    ! Require floats to be in layer of min thickness
                    ll = .true.

                    call get_command_argument(arg_count, arg_real, arg_len, exitcode)
                    arg_count = arg_count + 1
                    call checkexit(exitcode)
                    read(arg_real,'(f10.4)') min_thick

                    call get_command_argument(arg_count, arg_real, arg_len, exitcode)
                    arg_count = arg_count + 1
                    call checkexit(exitcode)
                    read(arg_real,'(f10.4)') max_thick

                    write(*,'(a10,a2,a7,f10.4,x,f10.4)') 'Property: ',trim(arg_flag),' spans ',min_thick,max_thick
                 endif

              elseif ( index(trim(arg_flag),'-F') .ne. 0 ) then
                 ! Choose the float closest to the
                 ! specified fraction of depth of
                 ! water column at launch.
                 lf = .true.
                 lc = .true.      ! We must have date info!
                 !ltt = .true.     ! We must always be in the layer
                 !lbb = .true.     ! so automatically set ltt, lbb
                 llayer = .true.

                 ! Check that we're looking at the starting
                 ! positions of the floats
                 if ( lstart ) then
                    ! We're OK, read in layer thickness fraction
                    call get_command_argument(arg_count, arg_real, arg_len, exitcode)
                    arg_count = arg_count + 1
                    call checkexit(exitcode)
                    read(arg_real,'(f10.4)') thick_frac
                 else
                    write(*,*) 'ERROR: -F option only works for start mode.'
                    stop
                 endif

              elseif ( index(trim(arg_flag),'-G') .ne. 0 ) then
                 ! Choose the float with the minimum
                 ! stratification at launch.
                 lg = .true.
                 lc = .true.      ! We must have date info!
                 !ltt = .true.     ! We must always be in the layer
                 !lbb = .true.     ! so automatically set ltt, lbb
                 llayer = .true.

                 ! Check that we're looking at starting
                 ! location of floats
                 if ( lstart ) then
                    ! We're OK, move on.
                 else
                    write(*,*) 'ERROR: -G option only works for start mode.'
                    stop
                 endif

              elseif ( index(trim(arg_flag),'-H') .ne. 0 ) then
                 ! Use bottom depth (bdep) as a sort criterium
                 lh = .true.

                 if (index(trim(arg_flag),'all') .ne. 0 ) then
                    ! We do not want to modify range of this variable,
                    ! but we do want to write it to output.
                 else
                    ! Get max and min values.
                    call get_command_argument(arg_count, arg_real, arg_len, exitcode)
                    arg_count = arg_count + 1
                    call checkexit(exitcode)
                    read(arg_real,'(f10.4)') h1
                 
                    call get_command_argument(arg_count, arg_real, arg_len, exitcode)
                    arg_count = arg_count + 1
                    call checkexit(exitcode)
                    read(arg_real,'(f10.4)') h2

                    call check_limits(h1,h2)

                    write(*,'(a10,a2,a7,f10.4,x,f10.4)') 'Property: ',trim(arg_flag),' spans ',h1,h2
                 endif

              elseif ( index(trim(arg_flag),'-J') .ne. 0 ) then
                 ! Keep all indefinite linesin/linesout.
                 lj = .true.
                 if ( (llinesin .or. llinesout) ) then
                    ! This choice makes sense
                 else
                    ! This flag does nothing here!
                    write(*,*) 'WARNING: -J option ignored b/c not in lines mode.'
                    lj = .false.
                 endif

              elseif ( index(trim(arg_flag),'-K') .ne. 0 ) then
                 ! Read in limits to fractional depth within layer
                 lk = .true.
                 llayer = .true.

                 ! Get max and min values.
                 call get_command_argument(arg_count, arg_real, arg_len, exitcode)
                 arg_count = arg_count + 1
                 call checkexit(exitcode)
                 read(arg_real,'(f10.4)') k1

                 call get_command_argument(arg_count, arg_real, arg_len, exitcode)
                 arg_count = arg_count + 1
                 call checkexit(exitcode)
                 read(arg_real,'(f10.4)') k2
                 
                 call check_limits(k1,k2)

                 write(*,'(a10,a2,a7,f10.4,x,f10.4)') 'Property: ',trim(arg_flag),' spans ',k1,k2
              
              elseif ( index(trim(arg_flag),'-d') .ne. 0 ) then
                 ! Read in limits to day
                 lday = .true.

                 ! Need to activate date reading
                 lc = .true.

                 ! Get max and min values.
                 call get_command_argument(arg_count, arg_real, arg_len, exitcode)
                 arg_count = arg_count + 1
                 call checkexit(exitcode)
                 read(arg_real,'(f10.4)') d1

                 call get_command_argument(arg_count, arg_real, arg_len, exitcode)
                 arg_count = arg_count + 1
                 call checkexit(exitcode)
                 read(arg_real,'(f10.4)') d2
                 
                 call check_limits(d1,d2)

                 write(*,'(a10,a2,a7,f10.4,x,f10.4)') 'Property: ',trim(arg_flag),' spans ',d1,d2

              elseif ( index(trim(arg_flag),'-m') .ne. 0 ) then
                 ! Read in limits to month
                 lmon = .true.

                 ! Need to activate date reading
                 lc = .true.

                 ! Get max and min values.
                 call get_command_argument(arg_count, arg_real, arg_len, exitcode)
                 arg_count = arg_count + 1
                 call checkexit(exitcode)
                 read(arg_real,'(f10.4)') m1

                 call get_command_argument(arg_count, arg_real, arg_len, exitcode)
                 arg_count = arg_count + 1
                 call checkexit(exitcode)
                 read(arg_real,'(f10.4)') m2
                 
                 call check_limits(m1,m2)

                 write(*,'(a10,a2,a7,f10.4,x,f10.4)') 'Property: ',trim(arg_flag),' spans ',m1,m2

              elseif ( index(trim(arg_flag),'-y') .ne. 0 ) then
                 ! Read in limits to month
                 lyear = .true.

                 ! Need to activate date reading
                 lc = .true.

                 ! Get max and min values.
                 call get_command_argument(arg_count, arg_real, arg_len, exitcode)
                 arg_count = arg_count + 1
                 call checkexit(exitcode)
                 read(arg_real,'(f10.4)') yr1

                 call get_command_argument(arg_count, arg_real, arg_len, exitcode)
                 arg_count = arg_count + 1
                 call checkexit(exitcode)
                 read(arg_real,'(f10.4)') yr2
                 
                 call check_limits(yr1,yr2)

                 write(*,'(a10,a2,a7,f10.4,x,f10.4)') 'Property: ',trim(arg_flag),' spans ',yr1,yr2

              elseif ( index(trim(arg_flag),'-c') .ne. 0 ) then
                 ! Centered difference for velocity.
                 lcent = .true.

                 ! Read number of days between time steps
                 call get_command_argument(arg_count, arg_real, arg_len, exitcode)
                 arg_count = arg_count + 1
                 call checkexit(exitcode)
                 read(arg_real,'(f10.4)') tstep
                 
                 write(*,'(a10,a2,a4,f10.4)') 'Property: ',trim(arg_flag),' is ',tstep

              elseif ( index(trim(arg_flag),'-f') .ne. 0 ) then
                 ! Forward difference for velocity.
                 lfore = .true.

                 ! Read number of days between time steps
                 call get_command_argument(arg_count, arg_real, arg_len, exitcode)
                 arg_count = arg_count + 1
                 call checkexit(exitcode)
                 read(arg_real,'(f10.4)') tstep
                 
                 write(*,'(a10,a2,a4,f10.4)') 'Property: ',trim(arg_flag),' is ',tstep

              elseif ( index(trim(arg_flag),'-k') .ne. 0 ) then
                 ! Backward difference for velocity.
                 lback = .true.

                 ! Read number of days between time steps
                 call get_command_argument(arg_count, arg_real, arg_len, exitcode)
                 arg_count = arg_count + 1
                 call checkexit(exitcode)
                 read(arg_real,'(f10.4)') tstep
                 
                 write(*,'(a10,a2,a4,f10.4)') 'Property: ',trim(arg_flag),' is ',tstep

              else
                 write(*,*) 'Error: Unexpected flag.'
                 stop
              endif

           enddo

        endif   ! Done reading command line options.
        !=======================================================

        !=======================================================
        ! Load up the domain limits vector for box_check
        !=======================================================
        all_min_max(1,1) = x1
        all_min_max(2,1) = x2

        all_min_max(1,2) = y1
        all_min_max(2,2) = y2

        all_min_max(1,3) = z1
        all_min_max(2,3) = z2

        all_min_max(1,4) = t1
        all_min_max(2,4) = t2

        all_min_max(1,5) = s1
        all_min_max(2,5) = s2

        all_min_max(1,6) = r1
        all_min_max(2,6) = r2

        all_min_max(1,7) = u1
        all_min_max(2,7) = u2

        all_min_max(1,8) = v1
        all_min_max(2,8) = v2

        all_min_max(1,9) = w1
        all_min_max(2,9) = w2

        all_min_max(1,10) = a1
        all_min_max(2,10) = a2

        all_min_max(1,11) = b1
        all_min_max(2,11) = b2

        all_min_max(1,12) = e1
        all_min_max(2,12) = e2

        all_min_max(1,13) = q1
        all_min_max(2,13) = q2

        all_min_max(1,14) = k1
        all_min_max(2,14) = k2

        all_min_max(1,15) = d1
        all_min_max(2,15) = d2

        all_min_max(1,16) = m1
        all_min_max(2,16) = m2

        all_min_max(1,17) = yr1
        all_min_max(2,17) = yr2

        all_min_max(1,18) = h1
        all_min_max(2,18) = h2

        !=======================================================
        ! Compute velocities if requested
        if ( lcent .and. lfore ) then
           write(*,*) 'ERROR: Cannot do both -c and -f!'
           stop
        endif

        if ( lcent .and. lback ) then
           write(*,*) 'ERROR: Cannot do both -c and -k!'
           stop
        endif

        if ( lback .and. lfore ) then
           write(*,*) 'ERROR: Cannot do both -k and -f!'
           stop
        endif

        if ( lcent .or. lfore .or. lback ) then
           write(*,*) 'Computing velocities from particle positions.'
           ! Initialize the unit conversions
           pi = 4.0*atan(1.0)
           deg2rad = pi/180.0
           deg2m = deg2rad*radius
           tstep = tstep*24.0*60.0*60.0
        endif

        if ( lnone .and. (lwormsin .or. lwormsout .or. &
                          llinesin .or. llinesout) ) then
           write(*,*) 'ERROR: Cannot specify -n with worms or lines!'
           stop
        endif

        ! Actually, one could stamp the ID on each traj.  Could be
        ! useful later, so keep it.
        !if ( l_stamp_id .and. .not. (lwormsin .or. lwormsout .or. &
        !                             llinesin .or. llinesout) ) then
        !   write(*,*) 'ERROR: Does not make sense to have -s without'
        !   write(*,*) 'linesin, linesout, wormsin, or wormsout.'
        !   stop
        !endif

        if ( l_stamp_id .and. l_keep_id ) then
           write(*,*) 'ERROR: Cannot specify stamp ID (-s) and'
           write(*,*) 'keep ID (-i) at the same time.'
           stop
        endif

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
        write(*,*) 'Get dimension IDs...'
        timedid = ncdid(ncid_i,timednam,exitcode)
        trajdid = ncdid(ncid_i,trajdnam,exitcode)

        write(*,*) 'Get dimension sizes...'
        call ncdinq(ncid_i,timedid,dummy,npts,exitcode)
        call ncdinq(ncid_i,trajdid,dummy,ntraj,exitcode)

        write(*,*) 'There are ',npts,' points in each traj.'
        write(*,*) 'There are ',ntraj,' trajectories.'

        ! Check input ntraj and fastest ntraj
        if ( (lfastest) .and. (nf .ge. ntraj) ) then
           write(*,*) 'Error: Too many fast traj. requested.'
           stop
        endif

        if ( lfastest .or. lthrough .or. lfexit .or. llexit ) then
           ! Allocate space for vector storing the number of
           ! steps until reaching the box.  Initialize for
           ! impossible to reach the box for all traj.
           !
           ! In the case for lfexit, n2box represents the
           ! number of time steps before the first time
           ! the traj reaches the edge of the box.
           !
           ! In the case for llexit, n2box represents the
           ! number of time steps before the last time
           ! the traj reaches the edge of the box (it is
           ! about to exit).
           allocate(n2box(ntraj))
           !if ( lfexit .or. llexit ) then
           !   n2box = 0
           !else
           n2box = npts + 1
           !endif
        endif

        if ( lnumber ) then
           ! Check that the largest index is within the
           ! number of available trajectories.
           if ( maxval(nkeep) .gt. ntraj ) then
              write(*,*) 'Error: max index value in index list is bigger than ntraj.'
              stop
           endif
           
           ! Check that we do not request more trajectories
           ! than are available
           if ( nf .gt. ntraj ) then
              write(*,*) 'Error: too many traj indices requested.'
              stop
           endif
        endif

        allocate(inyes(ntraj))
        inyes = 0
        yntraj = 0
        nntraj = 0
        
        write(*,*) 'Allocating space for one trajectory...'
        allocate(lam(npts))
        allocate(phi(npts))
        allocate(dep(npts))
        allocate(u(npts))
        allocate(v(npts))
        allocate(w(npts))
        allocate(temp(npts))
        allocate(salt(npts))
        allocate(rho(npts))
        allocate(q(npts))
        allocate(dtdz(npts))
        allocate(drdz(npts))
        allocate(age(npts))
        allocate(top(npts))
        allocate(bot(npts))
        allocate(year(npts))
        allocate(month(npts))
        allocate(day(npts))
        allocate(bdep(npts))
        if (lg .or. lf) then
           allocate(xyzymdvi(8,ntraj))
           xyzymdvi = -999.0
        endif
        if ( lwormsin .or. lwormsout .or. &
             llinesin .or. llinesout ) then
           allocate(age_in(npts))
           allocate(age_out(npts))
        endif

        ! Initialize values to values within the default
        ! ranges.  If values are not read in from the
        ! input file, then these values will never reject
        ! trajectories.
        lam = lam_def
        phi = phi_def
        dep = dep_def
        u = u_def
        v = v_def
        w = w_def
        temp = t_def
        salt = s_def
        rho = r_def
        q = q_def
        dtdz = dtdz_def
        drdz = drdz_def
        age = age_def
        top = layer_top_def
        bot = layer_bot_def
        year = year_def
        month = month_def
        day = day_def
        bdep = bdep_def

        ! Specify read vectors
        lag_readst2d(1) = 1
        lag_readct2d(1) = npts
        lag_readct2d(2) = 1

        ! Specify write vectors
        lag_writest2d(1) = 1
        if ( lwormsin .or. lwormsout ) then
           lag_writect2d(1) = 2*nw
        else
           lag_writect2d(1) = npts
        endif
        lag_writect2d(2) = 1

        if ( .not. lnone ) then
           !------Create, define the output netcdf yes file-------        
           write(*,*) 'Creating output cdf yes file...'
           ncid_y = nccre(outyname,ncclobber,exitcode)

           write(*,*) 'Creating dimensions...'
           ! We set the traj dimension to unlimited because
           ! this is the last index in the variable.  Also,
           ! we are more likely to have more trajectories
           ! than numbers of points (at 3 day intervals).
           ! For most case, we let the length of the traj
           ! stay the same (npts).  However, for worms,
           ! traj will be cut to length of 2*nw.
           if ( lwormsin .or. lwormsout ) then
              timedid = ncddef(ncid_y,timednam,(2*nw),exitcode)
           else
              timedid = ncddef(ncid_y,timednam,npts,exitcode)
           endif
           trajdid = ncddef(ncid_y,trajdnam,0,exitcode)

           vdims_out(1) = timedid
           vdims_out(2) = trajdid

           call dup_lag_var(ncid_i,ncid_y,lamvnam)
           write(*,*) 'Created longitude variable.'
           call dup_lag_var(ncid_i,ncid_y,phivnam)
           write(*,*) 'Created latitude variable.'
           call dup_lag_var(ncid_i,ncid_y,depvnam)
           write(*,*) 'Created depth variable.'
           if (lu .and. ( .not. (lcent .or. lfore .or. lback) ) ) then
              ! Assume U is available from input file
              call dup_lag_var(ncid_i,ncid_y,uvnam)
              write(*,*) 'Created zonal velocity variable.'
           elseif ( lu .and. (lcent .or. lfore .or. lback) ) then
              ! Assume U is not available, define fresh
              vid = ncvdef(ncid_y,uvnam,ncshort,2,vdims_out,exitcode)
              call ncapt(ncid_y,vid,'add_offset',ncfloat,1,0.0,exitcode)
              call ncapt(ncid_y,vid,'scale_factor',ncfloat,1,0.001,exitcode)
              write(*,*) 'Created NEW zonal velocity variable.'
           endif
           if (lv .and. ( .not. (lcent .or. lfore .or. lback) ) ) then
              call dup_lag_var(ncid_i,ncid_y,vvnam)
              write(*,*) 'Created meridional velocity variable.'
           elseif (lv .and. (lcent .or. lfore .or. lback) ) then
              vid = ncvdef(ncid_y,vvnam,ncshort,2,vdims_out,exitcode)
              call ncapt(ncid_y,vid,'add_offset',ncfloat,1,0.0,exitcode)
              call ncapt(ncid_y,vid,'scale_factor',ncfloat,1,0.001,exitcode)
              write(*,*) 'Created NEW meridional velocity variable.'
           endif
           if (lw .and. ( .not. (lcent .or. lfore .or. lback) ) ) then
              call dup_lag_var(ncid_i,ncid_y,wvnam)
              write(*,*) 'Created vertical velocity variable.'
           elseif (lw .and. ( lcent .or. lfore .or. lback ) ) then
              vid = ncvdef(ncid_y,wvnam,ncshort,2,vdims_out,exitcode)
              call ncapt(ncid_y,vid,'add_offset',ncfloat,1,0.0,exitcode)
              call ncapt(ncid_y,vid,'scale_factor',ncfloat,1,0.001,exitcode)
              write(*,*) 'Created NEW vertical velocity variable.'
           endif
           if ( l_stamp_id .or. l_keep_id ) then
              ! Need to include ID information in output
              vid = ncvdef(ncid_y,id_vnam,ncint,1,trajdid,exitcode)
           endif
           if (lt) then
              call dup_lag_var(ncid_i,ncid_y,tempvnam)
              write(*,*) 'Created temperature variable.'
           endif
           if (ls) then
              call dup_lag_var(ncid_i,ncid_y,saltvnam)
              write(*,*) 'Created salinity variable.'
           endif
           if (lr) then
              call dup_lag_var(ncid_i,ncid_y,rhovnam)
              write(*,*) 'Created density variable.'
           endif
           if (lq) then
              call dup_lag_var(ncid_i,ncid_y,qvnam)
              write(*,*) 'Created potential vorticity variable.'
           endif
           if (la) then
              ! SMALL OUTPUT
              call dup_lag_var(ncid_i,ncid_y,dtdz_vnam)
              write(*,*) 'Created dtdz variable.'
           endif
           if (lb) then
              call dup_lag_var(ncid_i,ncid_y,drdz_vnam)
              write(*,*) 'Created drdz variable.'
           endif
           if (lh) then
              call dup_lag_var(ncid_i,ncid_y,bdepvnam)
              write(*,*) 'Created bdep variable.'
           endif
           if ( lwormsin .or. lwormsout ) then
              ! Do not copy age from file - create
              ! locally.  Worms get age
              ! reset or constructed from files that
              ! did not originally have an age var.
              ! Therefore, age in input file is
              ! automatically, and always, ignored.
              !
              ! Lines get their age reset simply
              ! by being constructed as lines,
              ! but a separate age variable is not
              ! needed.
              
              ! Create variable in file
              write(*,*) 'Creating age variable...'
              vid = ncvdef(ncid_y,agevnam,ncshort,2,vdims_out,exitcode)
              
              ! Create add_offset and scale_factor attributes
              call ncapt(ncid_y,vid,'add_offset',ncfloat,1,32000.0,exitcode)
              call ncapt(ncid_y,vid,'scale_factor',ncfloat,1,1.0,exitcode)
              
           elseif (le .and. (.not. lage_const) ) then
              ! The user has still specified age
              ! sorting, so it must already be
              ! in the file (and if not, crash).
              ! However, if the user knowingly
              ! adds -e to the command line (saying
              ! "I know there is no age variable in this
              ! input file, but I still want to sort
              ! by length", do not carry over age info.
              ! Will be loaded later.
              call dup_lag_var(ncid_i,ncid_y,agevnam)
              write(*,*) 'Created age variable.'
           endif
           if (lc) then
              call dup_lag_var(ncid_i,ncid_y,year_vnam)
              write(*,*) 'Created year variable.'
              call dup_lag_var(ncid_i,ncid_y,month_vnam)
              write(*,*) 'Created month variable.'
              call dup_lag_var(ncid_i,ncid_y,day_vnam)
              write(*,*) 'Created day variable.'
           endif
           if (llayer) then
              call dup_lag_var(ncid_i,ncid_y,layer_top_vnam)
              write(*,*) 'Created layer top variable.'
              call dup_lag_var(ncid_i,ncid_y,layer_bot_vnam)
              write(*,*) 'Created layer bottom variable.'
           endif
           write(*,*) 'Done creating file...'
           call ncendf(ncid_y,exitcode)
        endif

        if (lj) then
           !------Create, define the indef. term. file-----------
           ! During command line check, this value would
           ! have been set to false if any mode other than
           ! lines in/out were selected.
           write(*,*) 'Creating output cdf indef. term. file...'
           ncid_j = nccre(outjname,ncclobber,exitcode)

           write(*,*) 'Creating dimensions...'
           timedid = ncddef(ncid_j,timednam,npts,exitcode)
           trajdid = ncddef(ncid_j,trajdnam,0,exitcode)

           vdims_out(1) = timedid
           vdims_out(2) = trajdid

           call dup_lag_var(ncid_i,ncid_j,lamvnam)
           write(*,*) 'Created longitude variable.'
           call dup_lag_var(ncid_i,ncid_j,phivnam)
           write(*,*) 'Created latitude variable.'
           call dup_lag_var(ncid_i,ncid_j,depvnam)
           write(*,*) 'Created depth variable.'
           if (lu .and. ( .not. (lcent .or. lfore .or. lback) ) ) then
              call dup_lag_var(ncid_i,ncid_j,uvnam)
              write(*,*) 'Created zonal velocity variable.'
           elseif (lu .and. (lcent .or. lfore .or. lback) ) then
              vid = ncvdef(ncid_j,uvnam,ncshort,2,vdims_out,exitcode)
              call ncapt(ncid_j,vid,'add_offset',ncfloat,1,0.0,exitcode)
              call ncapt(ncid_j,vid,'scale_factor',ncfloat,1,0.001,exitcode)
              write(*,*) 'Created NEW zonal velocity variable.'
           endif
           if (lv .and. ( .not. (lcent .or. lfore .or. lback) ) ) then
              call dup_lag_var(ncid_i,ncid_j,vvnam)
              write(*,*) 'Created meridional velocity variable.'
           elseif (lv .and. (lcent .or. lfore .or. lback) ) then
              vid = ncvdef(ncid_j,vvnam,ncshort,2,vdims_out,exitcode)
              call ncapt(ncid_j,vid,'add_offset',ncfloat,1,0.0,exitcode)
              call ncapt(ncid_j,vid,'scale_factor',ncfloat,1,0.001,exitcode)
              write(*,*) 'Created NEW meridional velocity variable.'
           endif
           if (lw .and. ( .not. (lcent .or. lfore .or. lback) ) ) then
              call dup_lag_var(ncid_i,ncid_j,wvnam)
              write(*,*) 'Created vertical velocity variable.'
           elseif (lw .and. (lcent .or. lfore .or. lback) ) then
              vid = ncvdef(ncid_j,wvnam,ncshort,2,vdims_out,exitcode)
              call ncapt(ncid_j,vid,'add_offset',ncfloat,1,0.0,exitcode)
              call ncapt(ncid_j,vid,'scale_factor',ncfloat,1,0.001,exitcode)
              write(*,*) 'Created NEW vertical velocity variable.'
           endif
           if ( l_stamp_id .or. l_keep_id ) then
              ! Need to include ID information in output
              vid = ncvdef(ncid_j,id_vnam,ncint,1,trajdid,exitcode)
           endif
           if (lt) then
              call dup_lag_var(ncid_i,ncid_j,tempvnam)
              write(*,*) 'Created temperature variable.'
           endif
           if (ls) then
              call dup_lag_var(ncid_i,ncid_j,saltvnam)
              write(*,*) 'Created salinity variable.'
           endif
           if (lr) then
              call dup_lag_var(ncid_i,ncid_j,rhovnam)
              write(*,*) 'Created density variable.'
           endif
           if (lq) then
              call dup_lag_var(ncid_i,ncid_j,qvnam)
              write(*,*) 'Created potential vorticity variable.'
           endif
           if (la) then
              call dup_lag_var(ncid_i,ncid_j,dtdz_vnam)
              write(*,*) 'Created dtdz variable.'
           endif
           if (lb) then
              call dup_lag_var(ncid_i,ncid_j,drdz_vnam)
              write(*,*) 'Created drdz variable.'
           endif
           if (le .and. (.not. lage_const) ) then
              ! The user has still specified age
              ! sorting, so it must already be
              ! in the file (and if not, crash).
              ! However, if the user knowingly
              ! adds -e to the command line (saying
              ! "I know there is no age variable in this
              ! input file, but I still want to sort
              ! by length", do not carry over age info.
              ! Will be loaded later.
              call dup_lag_var(ncid_i,ncid_j,agevnam)
              write(*,*) 'Created age variable.'
           endif
           if (lc) then
              call dup_lag_var(ncid_i,ncid_j,year_vnam)
              write(*,*) 'Created year variable.'
              call dup_lag_var(ncid_i,ncid_j,month_vnam)
              write(*,*) 'Created month variable.'
              call dup_lag_var(ncid_i,ncid_j,day_vnam)
              write(*,*) 'Created day variable.'
           endif
           if (llayer) then
              call dup_lag_var(ncid_i,ncid_j,layer_top_vnam)
              write(*,*) 'Created layer top variable.'
              call dup_lag_var(ncid_i,ncid_j,layer_bot_vnam)
              write(*,*) 'Created layer bottom variable.'
           endif
           if ( lh ) then
              call dup_lag_var(ncid_i,ncid_j,bdepvnam)
              write(*,*) 'Created bottom depth variable.'
           endif
           write(*,*) 'Done creating file...'
           call ncendf(ncid_j,exitcode)
        endif

        !------Create, define the output netcdf no file-------
        if (lno .and. ( lstart .or. lthrough .or. lfastest .or.&
                        lfexit .or. llexit .or. lnumber .or. lend ) .and. &
                        .not. lnone ) then
           ! Rejected trajectory file is only created for modes:
           ! start, through, fastest, lexit, fexit, number, and end.
           ! It does not make sense to have a no file
           ! for worms or lines since we're chopping up
           ! all traj into subsegments in those modes.

           write(*,*) 'Creating output cdf no file...'
           ncid_n = nccre(outnname,ncclobber,exitcode)

           write(*,*) 'Creating dimensions...'
           ! We set the traj dimension to unlimited because
           ! this is the last index in the variable.  Also,
           ! we are more likely to have more trajectories
           ! than numbers of points (at 3 day intervals).
           timedid = ncddef(ncid_n,timednam,npts,exitcode)
           trajdid = ncddef(ncid_n,trajdnam,0,exitcode)
           
           vdims_out(1) = timedid
           vdims_out(2) = trajdid
           
           write(*,*) 'Creating variables...'
           call dup_lag_var(ncid_i,ncid_n,lamvnam)
           call dup_lag_var(ncid_i,ncid_n,phivnam)
           call dup_lag_var(ncid_i,ncid_n,depvnam)
           if (lu .and. ( .not. (lcent .or. lfore .or. lback) ) ) then
              call dup_lag_var(ncid_i,ncid_n,uvnam)
           elseif (lu .and. (lcent .or. lfore .or. lback) ) then
              vid = ncvdef(ncid_n,uvnam,ncshort,2,vdims_out,exitcode)
              call ncapt(ncid_n,vid,'add_offset',ncfloat,1,0.0,exitcode)
              call ncapt(ncid_n,vid,'scale_factor',ncfloat,1,0.001,exitcode)
           endif
           if (lv .and. ( .not. (lcent .or. lfore .or. lback) ) ) then
              call dup_lag_var(ncid_i,ncid_n,vvnam)
           elseif (lv .and. (lcent .or. lfore .or. lback) ) then
              vid = ncvdef(ncid_n,vvnam,ncshort,2,vdims_out,exitcode)
              call ncapt(ncid_n,vid,'add_offset',ncfloat,1,0.0,exitcode)
              call ncapt(ncid_n,vid,'scale_factor',ncfloat,1,0.001,exitcode)
           endif
           if (lw .and. ( .not. (lcent .or. lfore .or. lback) ) ) then
              call dup_lag_var(ncid_i,ncid_n,wvnam)
           elseif (lw .and. (lcent .or. lfore .or. lback) ) then
              vid = ncvdef(ncid_n,wvnam,ncshort,2,vdims_out,exitcode)
              call ncapt(ncid_n,vid,'add_offset',ncfloat,1,0.0,exitcode)
              call ncapt(ncid_n,vid,'scale_factor',ncfloat,1,0.001,exitcode)
           endif
           if ( l_stamp_id .or. l_keep_id ) then
              ! Need to include ID information in output
              vid = ncvdef(ncid_n,id_vnam,ncint,1,trajdid,exitcode)
           endif
           if (lt) call dup_lag_var(ncid_i,ncid_n,tempvnam)
           if (ls) call dup_lag_var(ncid_i,ncid_n,saltvnam)
           if (lr) call dup_lag_var(ncid_i,ncid_n,rhovnam)
           if (lq) call dup_lag_var(ncid_i,ncid_n,qvnam)
           if (la) call dup_lag_var(ncid_i,ncid_n,dtdz_vnam)
           if (lb) call dup_lag_var(ncid_i,ncid_n,drdz_vnam)
           if (le) call dup_lag_var(ncid_i,ncid_n,agevnam)
           if (lc) then
              call dup_lag_var(ncid_i,ncid_n,year_vnam)
              call dup_lag_var(ncid_i,ncid_n,month_vnam)
              call dup_lag_var(ncid_i,ncid_n,day_vnam)
           endif
           if (llayer) then
              call dup_lag_var(ncid_i,ncid_n,layer_top_vnam)
              call dup_lag_var(ncid_i,ncid_n,layer_bot_vnam)
           endif
           if (lh) call dup_lag_var(ncid_i,ncid_n,bdepvnam)
           write(*,*) 'Done creating file...'
           call ncendf(ncid_n,exitcode)
        endif

        !=======================================================
        !                  SORT TRAJECTORIES
        !=======================================================
        if ( lnumber ) then
           ! We don't need to sort, instead,
           ! update the inyes and yntraj look
           ! up tables.
           yntraj = nf
           do t = 1,nf
              inyes(nkeep(t)) = 1
           enddo

           ! NOTE: Currently the whole file is loaded
           ! when running in -number mode - a special
           ! check could be made for each traj to
           ! not load data if in number mode for traj
           ! we don't want.  However, that would come
           ! at the cost of slowing down the other
           ! operators.
        else
           ! Other operators require sorting.

           ! Initialize the counter of worms or lines
           if ( lwormsin .or. lwormsout .or. &
                llinesin .or. llinesout ) then
              y = 0
              y_j = 0
           endif

           write(*,*) 'Sorting: Loop over each trajectory...'
           do t = 1,ntraj
           
              !write(*,*) 'For traj ',t,' out of ',ntraj
        
              !write(*,*) 'Specify the read vectors...'
              lag_readst2d(2) = t

              !write(*,*) 'Load traj data...'
              inlayer = .true.
              !--------------------------------------------------------
              ! Load position information
              call get_lag_var(ncid_i,lamvnam,lam)
              call get_lag_var(ncid_i,phivnam,phi)
              call get_lag_var(ncid_i,depvnam,dep)
              !--------------------------------------------------------
              ! Load velocity information
              if (lu .and. (.not. (lcent .or. lfore .or. lback))) then
                 call get_lag_var(ncid_i,uvnam,u)
              endif
              if (lv .and. (.not. (lcent .or. lfore .or. lback))) then
                 call get_lag_var(ncid_i,vvnam,v)
              endif
              if (lw .and. (.not. (lcent .or. lfore .or. lback))) then
                 call get_lag_var(ncid_i,wvnam,w)
              endif
              !--------------------------------------------------------
              ! Compute velocity from the position information with
              ! choice of three different differencing schemes.
              if (lcent) then
                 do p = 2,npts-1
                    ! Use the avg. latitude between the centered
                    ! difference endpoints to determine the latitude
                    ! spacing degrees to meter conversion.  This is
                    ! instead of the center point because the center
                    ! point is not necessarily in the center of the
                    ! two centered difference endpoints.
                    !        change in lon      to m        account for sphere
                    u(p) = ((lam(p+1)-lam(p-1))*deg2m*cos(((phi(p+1)+phi(p-1))/2.0)*deg2rad))/(2.0*tstep)
                    v(p) = ((phi(p+1)-phi(p-1))*deg2m)/(2.0*tstep)
                    w(p) = (dep(p+1)-dep(p-1))/(2.0*tstep)
                 enddo
                 ! Endpoints
                 u(1) = u(2)
                 v(1) = v(2)
                 w(1) = w(2)
                 u(npts) = u(npts-1)
                 v(npts) = v(npts-1)
                 w(npts) = w(npts-1)
              elseif (lfore) then
                 do p = 1,npts-1
                    u(p) = ((lam(p+1)-lam(p))*deg2m*cos(((phi(p+1)+phi(p))/2.0)*deg2rad))/tstep
                    v(p) = ((phi(p+1)-phi(p))*deg2m)/tstep
                    w(p) = (dep(p+1)-dep(p))/tstep
                 enddo
                 u(npts) = u(npts-1)
                 v(npts) = v(npts-1)
                 w(npts) = w(npts-1)
              elseif (lback) then
                 do p = 2,npts
                    u(p) = ((lam(p)-lam(p-1))*deg2m*cos(((phi(p)+phi(p-1))/2.0)*deg2rad))/tstep
                    v(p) = ((phi(p)-phi(p-1))*deg2m)/tstep
                    w(p) = (dep(p)-dep(p-1))/tstep
                 enddo
                 u(1) = u(2)
                 v(1) = v(2)
                 w(1) = w(2)
              endif
              !--------------------------------------------------------
              if ( l_keep_id ) then
                 ! Get the vid of the id variable
                 vid = ncvid(ncid_i,id_vnam,exitcode)
                 
                 ! Load the ID value into current_id
                 call ncvgt1(ncid_i,vid,t,current_id,exitcode)
              else
                 current_id = t
              endif

              ! Load scalar variables
              if (lt) call get_lag_var(ncid_i,tempvnam,temp)
              if (ls) call get_lag_var(ncid_i,saltvnam,salt)
              if (lr) call get_lag_var(ncid_i,rhovnam,rho)
              if (lq) call get_lag_var(ncid_i,qvnam,q)
              if (la) call get_lag_var(ncid_i,dtdz_vnam,dtdz)
              if (lb) call get_lag_var(ncid_i,drdz_vnam,drdz)
              if (le .and. .not. (lwormsin.or.lwormsout.or.&
                                  llinesin.or.llinesout.or.&
                                  lage_const )) then
                 call get_lag_var(ncid_i,agevnam,age)
              else
                 ! Worms and lines operations will reset
                 ! the age by default, so do not load
                 ! any age information for them.
                 ! Rather, set a default age by the
                 ! time step counter.
                 ! "blind" aging - don't check for
                 ! valid length.
                 do p = 1,npts
                    age(p) = p
                 enddo
              endif
              ! Any age information loaded above will be
              ! overwritten by valid age information
              ! as determined by the valid points in
              ! the trajectory.
              if ( lage_const ) then
                 ! Check for valid length, and if
                 ! not valid, flag with mask.
                 do p = 1,npts
                    if ( (lam(p) .eq. lam_mask) .or.&
                         (phi(p) .eq. phi_mask) .or.&
                         (dep(p) .eq. dep_mask) ) then
                       age(p) = age_mask
                    else
                       age(p) = p
                    endif
                 enddo
              endif
              if (lc) call get_lag_var(ncid_i,year_vnam,year)
              if (lc) call get_lag_var(ncid_i,month_vnam,month)
              if (lc) call get_lag_var(ncid_i,day_vnam,day)
              if (llayer) call get_lag_var(ncid_i,layer_top_vnam,top)
              if (llayer) call get_lag_var(ncid_i,layer_bot_vnam,bot)
              if (lh) call get_lag_var(ncid_i,bdepvnam,bdep)

              !===============================================
              ! SORT STARTING POSITIONS
              !===============================================
              if ( lstart ) then
                 ! Test for position within layer
                 inlayer = layer_check(dep(1),top(1),bot(1),&
                                       min_thick,max_thick,ll,lo,ltt,lbb)

                 frac_tmp = 0.0
                 if ( lk ) then
                    ! It is worth computing the fractional
                    ! layer location.
                    frac_tmp = ((bot(1) - dep(1))/(bot(1) - top(1)))
                 endif

                 ! Load the point vector
                 all_pt_data(1) = lam(1)
                 all_pt_data(2) = phi(1)
                 all_pt_data(3) = dep(1)
                 all_pt_data(4) = temp(1)
                 all_pt_data(5) = salt(1)
                 all_pt_data(6) = rho(1)
                 all_pt_data(7) = u(1)
                 all_pt_data(8) = v(1)
                 all_pt_data(9) = w(1)
                 all_pt_data(10) = dtdz(1)
                 all_pt_data(11) = drdz(1)
                 all_pt_data(12) = age(1)
                 all_pt_data(13) = q(1)
                 all_pt_data(14) = frac_tmp
                 all_pt_data(15) = day(1)
                 all_pt_data(16) = month(1)
                 all_pt_data(17) = year(1)
                 all_pt_data(18) = bdep(1)

                 ! Test for start position.
                 !if ( (lam(1) .ge. x1) .and. (lam(1) .le. x2) .and. &
                 !     (phi(1) .ge. y1) .and. (phi(1) .le. y2) .and. &
                 !     (dep(1) .ge. z1) .and. (dep(1) .le. z2) .and. &
                 !     (temp(1) .ge. t1) .and. (temp(1) .le. t2) .and. &
                 !     (salt(1) .ge. s1) .and. (salt(1) .le. s2) .and. &
                 !     (rho(1) .ge. r1) .and. (rho(1) .le. r2) .and. &
                 !     (u(1) .ge. u1) .and. (u(1) .le. u2) .and. &
                 !     (v(1) .ge. v1) .and. (v(1) .le. v2) .and. &
                 !     (w(1) .ge. w1) .and. (w(1) .le. w2) .and. &
                 !     (dtdz(1) .ge. a1) .and. (dtdz(1) .le. a2) .and. &
                 !     (drdz(1) .ge. b1) .and. (drdz(1) .le. b2) .and. &
                 !     (age(1) .ge. e1) .and. (age(1) .le. e2) .and. &
                 !     (q(1) .ge. q1) .and. (q(1) .le. q2) .and. &
                 !     (frac_tmp .ge. k1) .and. (frac_tmp .le. k2) .and. &
                 !     (month(1) .ge. m1) .and. (month(1) .le. m2) .and. &
                 !     (year(1) .ge. yr1) .and. (year(1) .le. yr2) .and. inlayer ) then
                 if ( box_check(all_min_max,all_pt_data,n_limits,l_equality) .and. inlayer ) then
                    ! This start point is within the start zone and we
                    ! are sorting only the start of each trajectory.

                    if (lg .or. lf) then
                       ! This float could be a float that we're
                       ! looking for, but we need to make certain
                       ! it's the best one of its launch column.

                       ! For every possible entry in the list of
                       ! launch locations, xyzymdvi, see if the
                       ! current float launch location matches
                       ! an already detected launch location.
                       ! If yntraj is 0, then this loop will
                       ! not execute and a new float will be
                       ! added to the xyzymdvi list.
                       lfound = .false.
                       do p = 1,yntraj
                          if ( (lam(1) .eq. xyzymdvi(1,p)) .and. &
                               (phi(1) .eq. xyzymdvi(2,p)) .and. &
                               (year(1) .eq. xyzymdvi(4,p)) .and. &
                               (month(1) .eq. xyzymdvi(5,p)) .and. &
                               (day(1) .eq. xyzymdvi(6,p)) ) then
                             lfound = .true.
                             exit
                          endif
                       enddo

                       if (lfound) then
                          ! A duplicate launch location was found.
                          if (lg) then
                             ! Pick the float with the lowest dtdz
                             ! or if the dtdz is the same, pick the
                             ! deeper float.
                             if ( (xyzymdvi(7,p) .gt. dtdz(1)) .or. &
                                 ((xyzymdvi(7,p) .eq. dtdz(1)) .and. xyzymdvi(3,p) .lt. dep(1))) then
                                ! This newly discovered float is a
                                ! better match.  Replace the previous
                                ! float with this one.
                                inyes(int(xyzymdvi(8,p))) = 0
                                inyes(t) = 1

                                xyzymdvi(3,p) = dep(1)
                                xyzymdvi(7,p) = dtdz(1)
                                xyzymdvi(8,p) = t
                             else
                                ! The previously discovered float at
                                ! this launch location is a better match.
                                ! Do nothing to change.
                             endif
                          elseif (lf) then
                             ! Pick the float closest to the given
                             ! fraction of the column thickness.
                             frac_tmp = (dep(1)-top(1))/(bot(1)-top(1))
                             if ( abs(xyzymdvi(7,p)-thick_frac) .gt. &
                                  abs(frac_tmp - thick_frac) ) then
                                ! This newly discovered float is a
                                ! better match.  Replace the previous
                                ! float with this one.
                                inyes(int(xyzymdvi(8,p))) = 0
                                inyes(t) = 1

                                xyzymdvi(3,p) = dep(1)
                                xyzymdvi(7,p) = frac_tmp
                                xyzymdvi(8,p) = t
                             else
                                ! The previously discovered float at
                                ! this launch location is a better match.
                                ! Do nothing to change.
                             endif
                          endif
                       else
                          ! Add a new entry into the launch location
                          ! list because no duplicate (x,y) was found.
                          ! The only time that the number of yes traj
                          ! can be augmented in lf or lg mode is when
                          ! we find a new launch location in xyt space.
                          inyes(t) = 1
                          yntraj = yntraj + 1
                          xyzymdvi(1,yntraj) = lam(1)
                          xyzymdvi(2,yntraj) = phi(1)
                          xyzymdvi(3,yntraj) = dep(1)
                          xyzymdvi(4,yntraj) = year(1)
                          xyzymdvi(5,yntraj) = month(1)
                          xyzymdvi(6,yntraj) = day(1)
                          if (lg) then
                             ! Comparison value is statification
                             xyzymdvi(7,yntraj) = dtdz(1)
                          elseif (lf) then
                             ! Comparison value is frac depth in layer.
                             xyzymdvi(7,yntraj) = (dep(1)-top(1))/(bot(1)-top(1))
                          endif
                          xyzymdvi(8,yntraj) = real(t)
                       endif
                    else
                       ! Just log the locations of the valid
                       ! floats and do not check launch locations
                       ! because no special lf or lg mode.
                       inyes(t) = 1
                       yntraj = yntraj + 1
                    endif
                 else
                    ! Do nothing here because the trajectory is not
                    ! in the selected launched box.
                 endif

              !===============================================
              ! SORT ENDING POSITIONS
              !===============================================
              elseif ( lend ) then
                 ! Test for position within layer
                 inlayer = layer_check(dep(npts),top(npts),bot(npts),&
                                       min_thick,max_thick,ll,lo,ltt,lbb)

                 frac_tmp = 0.0
                 if ( lk ) then
                    ! It is worth computing the fractional
                    ! layer location.
                    frac_tmp = ((bot(npts) - dep(npts))/(bot(npts) - top(npts)))
                 endif

                 ! Load the point vector
                 all_pt_data(1) = lam(npts)
                 all_pt_data(2) = phi(npts)
                 all_pt_data(3) = dep(npts)
                 all_pt_data(4) = temp(npts)
                 all_pt_data(5) = salt(npts)
                 all_pt_data(6) = rho(npts)
                 all_pt_data(7) = u(npts)
                 all_pt_data(8) = v(npts)
                 all_pt_data(9) = w(npts)
                 all_pt_data(10) = dtdz(npts)
                 all_pt_data(11) = drdz(npts)
                 all_pt_data(12) = age(npts)
                 all_pt_data(13) = q(npts)
                 all_pt_data(14) = frac_tmp
                 all_pt_data(15) = day(npts)
                 all_pt_data(16) = month(npts)
                 all_pt_data(17) = year(npts)
                 all_pt_data(18) = bdep(npts)

                 ! Test end position.
                 !if ( (lam(npts) .ge. x1) .and. (lam(npts) .le. x2) .and. &
                 !     (phi(npts) .ge. y1) .and. (phi(npts) .le. y2) .and. &
                 !     (dep(npts) .ge. z1) .and. (dep(npts) .le. z2) .and. &
                 !     (temp(npts) .ge. t1) .and. (temp(npts) .le. t2) .and. &
                 !     (salt(npts) .ge. s1) .and. (salt(npts) .le. s2) .and. &
                 !     (rho(npts) .ge. r1) .and. (rho(npts) .le. r2) .and. &
                 !     (u(npts) .ge. u1) .and. (u(npts) .le. u2) .and. &
                 !     (v(npts) .ge. v1) .and. (v(npts) .le. v2) .and. &
                 !     (w(npts) .ge. w1) .and. (w(npts) .le. w2) .and. &
                 !     (dtdz(npts) .ge. a1) .and. (dtdz(npts) .le. a2) .and. &
                 !     (drdz(npts) .ge. b1) .and. (drdz(npts) .le. b2) .and. &
                 !     (age(npts) .ge. e1) .and. (age(npts) .le. e2) .and. &
                 !     (q(npts) .ge. q1) .and. (q(npts) .le. q2) .and. &
                 !     (frac_tmp .ge. k1) .and. (frac_tmp .le. k2) .and. &
                 !     (month(npts) .ge. m1) .and. (month(npts) .le. m2) .and. &
                 !     (year(npts) .ge. yr1) .and. (year(npts) .le. yr2) .and. inlayer ) then
                 if ( box_check(all_min_max,all_pt_data,n_limits,l_equality) .and. inlayer ) then
                    ! This ending point is within the zone and we
                    ! are sorting only the end of each trajectory.

                    ! Just log the locations of the valid
                    ! floats and do not check launch locations
                    ! because no special lf or lg mode.
                    inyes(t) = 1
                    yntraj = yntraj + 1
                 else
                    ! Do nothing here because the trajectory is not
                    ! in the selected launched box.
                 endif

              !===============================================
              ! SEARCH FOR TRAJ THAT GO THROUGH THE BOX
              !===============================================
              elseif ( lthrough .or. lfastest ) then
                 ! Test that trajectory has passed through a box (at any time)
                 ! Count the number of points until the first arrival at the
                 ! box.
                 
                 ! Loop over each point in the trajectory
                 do p = 1,npts
                    
                    ! Test for position within layer
                    inlayer = layer_check(dep(p),top(p),bot(p),&
                                          min_thick,max_thick,ll,lo,ltt,lbb)

                    frac_tmp = 0.0
                    if ( lk ) then
                       ! It is worth computing the fractional
                       ! layer location.
                       frac_tmp = ((bot(p) - dep(p))/(bot(p) - top(p)))
                    endif

                    ! Load the point vector
                    all_pt_data(1) = lam(p)
                    all_pt_data(2) = phi(p)
                    all_pt_data(3) = dep(p)
                    all_pt_data(4) = temp(p)
                    all_pt_data(5) = salt(p)
                    all_pt_data(6) = rho(p)
                    all_pt_data(7) = u(p)
                    all_pt_data(8) = v(p)
                    all_pt_data(9) = w(p)
                    all_pt_data(10) = dtdz(p)
                    all_pt_data(11) = drdz(p)
                    all_pt_data(12) = age(p)
                    all_pt_data(13) = q(p)
                    all_pt_data(14) = frac_tmp
                    all_pt_data(15) = day(p)
                    all_pt_data(16) = month(p)
                    all_pt_data(17) = year(p)
                    all_pt_data(18) = bdep(p)

                    ! Test that we have gone through the box.
                    !if ( (lam(p) .ge. x1) .and. (lam(p) .le. x2) .and. &
                    !     (phi(p) .ge. y1) .and. (phi(p) .le. y2) .and. &
                    !     (dep(p) .ge. z1) .and. (dep(p) .le. z2) .and. &
                    !     (temp(p) .ge. t1) .and. (temp(p) .le. t2) .and. &
                    !     (salt(p) .ge. s1) .and. (salt(p) .le. s2) .and. &
                    !     (rho(p) .ge. r1) .and. (rho(p) .le. r2) .and. &
                    !     (u(p) .ge. u1) .and. (u(p) .le. u2) .and. &
                    !     (v(p) .ge. v1) .and. (v(p) .le. v2) .and. &
                    !     (w(p) .ge. w1) .and. (w(p) .le. w2) .and. &
                    !     (dtdz(p) .ge. a1) .and. (dtdz(p) .le. a2) .and. &
                    !     (drdz(p) .ge. b1) .and. (drdz(p) .le. b2) .and. &
                    !     (age(p) .ge. e1) .and. (age(p) .le. e2) .and. &
                    !     (q(p) .ge. q1) .and. (q(p) .le. q2) .and. &
                    !     (month(p) .ge. m1) .and. (month(p) .le. m2) .and. inlayer ) then
                    if ( box_check(all_min_max,all_pt_data,n_limits,l_equality) .and. inlayer ) then   
                       ! We have just entered/touched the box.
                       inyes(t) = 1
                       yntraj = yntraj + 1
                       
                       ! Exit this loop
                       exit
                       
                    else
                       ! Do nothing since the trajectory has not yet
                       ! reached the box.
                    endif
                    
                 enddo
                 
                 ! Store the number of points it takes to first
                 ! reach the box.
                 if ( inyes(t) .eq. 1 ) then
                    ! We have indeed found the point
                    n2box(t) = p
                 endif

              !===============================================
              ! SEARCH FOR LAST EXIT OF TRAJ
              !===============================================
              elseif ( llexit ) then
                 ! Test whether this traj exits the box.  Keep track
                 ! of the number of time steps that it took to exit
                 ! the box.

                 ! Loop over each point in the trajectory
                 do p = 1,npts
                    
                    ! Test for position within layer
                    inlayer = layer_check(dep(p),top(p),bot(p),&
                                          min_thick,max_thick,ll,lo,ltt,lbb)

                    frac_tmp = 0.0
                    if ( lk ) then
                       ! It is worth computing the fractional
                       ! layer location.
                       frac_tmp = ((bot(p) - dep(p))/(bot(p) - top(p)))
                    endif

                    ! Load the point vector
                    all_pt_data(1) = lam(p)
                    all_pt_data(2) = phi(p)
                    all_pt_data(3) = dep(p)
                    all_pt_data(4) = temp(p)
                    all_pt_data(5) = salt(p)
                    all_pt_data(6) = rho(p)
                    all_pt_data(7) = u(p)
                    all_pt_data(8) = v(p)
                    all_pt_data(9) = w(p)
                    all_pt_data(10) = dtdz(p)
                    all_pt_data(11) = drdz(p)
                    all_pt_data(12) = age(p)
                    all_pt_data(13) = q(p)
                    all_pt_data(14) = frac_tmp
                    all_pt_data(15) = day(p)
                    all_pt_data(16) = month(p)
                    all_pt_data(17) = year(p)
                    all_pt_data(18) = bdep(p)

                    ! Test that we are in the box.
                    !if ( (lam(p) .ge. x1) .and. (lam(p) .le. x2) .and. &
                    !     (phi(p) .ge. y1) .and. (phi(p) .le. y2) .and. &
                    !     (dep(p) .ge. z1) .and. (dep(p) .le. z2) .and. &
                    !     (temp(p) .ge. t1) .and. (temp(p) .le. t2) .and. &
                    !     (salt(p) .ge. s1) .and. (salt(p) .le. s2) .and. &
                    !     (rho(p) .ge. r1) .and. (rho(p) .le. r2) .and. &
                    !     (u(p) .ge. u1) .and. (u(p) .le. u2) .and. &
                    !     (v(p) .ge. v1) .and. (v(p) .le. v2) .and. &
                    !     (w(p) .ge. w1) .and. (w(p) .le. w2) .and. &
                    !     (dtdz(p) .ge. a1) .and. (dtdz(p) .le. a2) .and. &
                    !     (drdz(p) .ge. b1) .and. (drdz(p) .le. b2) .and. &
                    !     (age(p) .ge. e1) .and. (age(p) .le. e2) .and. &
                    !     (q(p) .ge. q1) .and. (q(p) .le. q2) .and. &
                    !     (month(p) .ge. m1) .and. (month(p) .le. m2) .and. inlayer ) then
                    if ( box_check(all_min_max,all_pt_data,n_limits,l_equality) .and. inlayer ) then  
                       ! This point is in the box.
                       inbox = .true.

                       ! Move to the next point to see if we exit...

                    else
                       ! We are not in the box, so this
                       ! traj will eventually be written to
                       ! the yes file.
                       inyes(t) = 1

                       ! However, we need to wait until the
                       ! very last time a traj exits the box.
                       ! Do not augment the number of yes floats!
                       
                       ! Test that that previous point was in
                       ! the box.
                       if ( inbox ) then
                          ! We have discovered a recent exit!
                          ! Log the number of the previous time
                          ! step when the particle was in the box,
                          ! just before exiting.
                          n2box(t) = p - 1
                          !write(*,*) 'traj ',t,' at point ',p,' went OUT!'
                       else
                          ! This is not a recent exit, do
                          ! nothing.
                       endif

                       ! Set the inbox flag to false because
                       ! this particle is no longer in the box.
                       ! The only time this flag can be reset to
                       ! true is if the particle re-enters the
                       ! box.
                       inbox = .false.

                    endif
                 enddo

              !===============================================
              ! SEARCH FOR THE FIRST EXIT OF THE TRAJ FROM BOX
              !===============================================
              elseif ( lfexit ) then
                 ! Test whether this traj exits the box.  Keep track
                 ! of the number of time steps that it took to get
                 ! to the first exit of the box.

                 ! Loop over each point in the trajectory
                 do p = 1,npts
                    
                    ! Test for position within layer
                    inlayer = layer_check(dep(p),top(p),bot(p),&
                                          min_thick,max_thick,ll,lo,ltt,lbb)

                    frac_tmp = 0.0
                    if ( lk ) then
                       ! It is worth computing the fractional
                       ! layer location.
                       frac_tmp = ((bot(p) - dep(p))/(bot(p) - top(p)))
                    endif

                    ! Load the point vector
                    all_pt_data(1) = lam(p)
                    all_pt_data(2) = phi(p)
                    all_pt_data(3) = dep(p)
                    all_pt_data(4) = temp(p)
                    all_pt_data(5) = salt(p)
                    all_pt_data(6) = rho(p)
                    all_pt_data(7) = u(p)
                    all_pt_data(8) = v(p)
                    all_pt_data(9) = w(p)
                    all_pt_data(10) = dtdz(p)
                    all_pt_data(11) = drdz(p)
                    all_pt_data(12) = age(p)
                    all_pt_data(13) = q(p)
                    all_pt_data(14) = frac_tmp
                    all_pt_data(15) = day(p)
                    all_pt_data(16) = month(p)
                    all_pt_data(17) = year(p)
                    all_pt_data(18) = bdep(p)

                    ! Test that we are in the box.
                    !if ( (lam(p) .ge. x1) .and. (lam(p) .le. x2) .and. &
                    !     (phi(p) .ge. y1) .and. (phi(p) .le. y2) .and. &
                    !     (dep(p) .ge. z1) .and. (dep(p) .le. z2) .and. &
                    !     (temp(p) .ge. t1) .and. (temp(p) .le. t2) .and. &
                    !     (salt(p) .ge. s1) .and. (salt(p) .le. s2) .and. &
                    !     (rho(p) .ge. r1) .and. (rho(p) .le. r2) .and. &
                    !     (u(p) .ge. u1) .and. (u(p) .le. u2) .and. &
                    !     (v(p) .ge. v1) .and. (v(p) .le. v2) .and. &
                    !     (w(p) .ge. w1) .and. (w(p) .le. w2) .and. &
                    !     (dtdz(p) .ge. a1) .and. (dtdz(p) .le. a2) .and. &
                    !     (drdz(p) .ge. b1) .and. (drdz(p) .le. b2) .and. &
                    !     (age(p) .ge. e1) .and. (age(p) .le. e2) .and. &
                    !     (q(p) .ge. q1) .and. (q(p) .le. q2) .and. &
                    !     (month(p) .ge. m1) .and. (month(p) .le. m2) .and. inlayer ) then
                    if ( box_check(all_min_max,all_pt_data,n_limits,l_equality) .and. inlayer ) then
                       ! We are in the box.
                       ! Move to the next point to see if we exit...
                    else
                       ! We are not in the box, so this
                       ! traj will eventually be written to
                       ! the yes file.
                       inyes(t) = 1

                       ! Report the position just before exit.
                       ! Note that if we are not in the box at
                       ! launch, then n2box = 0.
                       n2box(t) = p - 1

                       ! Augment counter
                       yntraj = yntraj + 1

                       ! Break the loop when the first exit
                       ! is detected.
                       exit

                    endif
                 enddo

              !==============================================
              ! CHOP TRAJ UP INTO WORMS OR LINES
              !==============================================
              elseif ( lwormsout .or. lwormsin .or. &
                       llinesin .or. llinesout ) then
                 ! Until now, all we have done is load
                 ! the information for the current
                 ! trajectory into memory.

                 ! The entry point and exit point each
                 ! have an age of 1.
                 age_in = 0
                 age_out = 0

                 ! Set the beginning point for a line (mask)
                 p_beg = 0

                 subseg_per_traj = 0

                 ! Loop over each point in the trajectory
                 do p = 1,npts

                    ! Test for position within layer
                    inlayer = layer_check(dep(p),top(p),bot(p),&
                                          min_thick,max_thick,ll,lo,ltt,lbb)

                    frac_tmp = 0.0
                    if ( lk ) then
                       ! It is worth computing the fractional
                       ! layer location.
                       frac_tmp = ((bot(p) - dep(p))/(bot(p) - top(p)))
                    endif

                    ! Load the point vector
                    all_pt_data(1) = lam(p)
                    all_pt_data(2) = phi(p)
                    all_pt_data(3) = dep(p)
                    all_pt_data(4) = temp(p)
                    all_pt_data(5) = salt(p)
                    all_pt_data(6) = rho(p)
                    all_pt_data(7) = u(p)
                    all_pt_data(8) = v(p)
                    all_pt_data(9) = w(p)
                    all_pt_data(10) = dtdz(p)
                    all_pt_data(11) = drdz(p)
                    all_pt_data(12) = age(p)
                    all_pt_data(13) = q(p)
                    all_pt_data(14) = frac_tmp
                    all_pt_data(15) = day(p)
                    all_pt_data(16) = month(p)
                    all_pt_data(17) = year(p)
                    all_pt_data(18) = bdep(p)

                    ! Test that we are in the box.
                    !if ( (lam(p) .ge. x1) .and. (lam(p) .le. x2) .and. &
                    !     (phi(p) .ge. y1) .and. (phi(p) .le. y2) .and. &
                    !     (dep(p) .ge. z1) .and. (dep(p) .le. z2) .and. &
                    !     (temp(p) .ge. t1) .and. (temp(p) .le. t2) .and. &
                    !     (salt(p) .ge. s1) .and. (salt(p) .le. s2) .and. &
                    !     (rho(p) .ge. r1) .and. (rho(p) .le. r2) .and. &
                    !     (u(p) .ge. u1) .and. (u(p) .le. u2) .and. &
                    !     (v(p) .ge. v1) .and. (v(p) .le. v2) .and. &
                    !     (w(p) .ge. w1) .and. (w(p) .le. w2) .and. &
                    !     (dtdz(p) .ge. a1) .and. (dtdz(p) .le. a2) .and. &
                    !     (drdz(p) .ge. b1) .and. (drdz(p) .le. b2) .and. &
                    !     (age(p) .ge. e1) .and. (age(p) .le. e2) .and. &
                    !     (q(p) .ge. q1) .and. (q(p) .le. q2) .and. &
                    !     (month(p) .ge. m1) .and. (month(p) .le. m2) .and. inlayer ) then
                    if ( box_check(all_min_max,all_pt_data,n_limits,l_equality) .and. inlayer ) then
                       ! We are in the box.
                       inbox_now = .true.

                       if ( p .eq. 1 ) then
                          ! If it is the first point, set equal
                          ! previous time step flag too.
                          inbox_old = inbox_now

                          ! If we are also searching for linesin,
                          ! then this first point, which is inside
                          ! the domain, is the beginning of the
                          ! first line in.
                          if ( llinesin ) then
                             p_beg = 1
                             subseg_per_traj = subseg_per_traj + 1
                          endif
                       endif

                       if ( inbox_now .eqv. inbox_old ) then
                          ! Current and old positions are
                          ! both in the box.  No entering
                          ! detected, so skip to next point.

                          ! Age the float in the domain.
                          if ( p .eq. 1 ) then
                             age_in(p) = 1
                          elseif ( (p .eq. npts) .and. lj .and. llinesin ) then
                             ! Line with an indefinite ending
                             ! terminates within the domain.
                             ! User has requested that it be
                             ! added to the other output.
                             p_end = npts
                             age_in(p) = age_in(p-1) + 1

                             ! The beginning of this line was noted
                             ! earlier in p_beg.  However, check for
                             ! consistency with the age vector:
                             if ( age_in(p_beg) .ne. 1 ) then
                                write(*,*) 'ERROR: age_in(p_beg) is not 1, rather =',age_in(p_beg)
                                write(*,*) 'p_beg =',p_beg
                                write(*,*) 'p_end =',p_end
                                write(*,*) 'traj  =',t
                                write(*,*) 'point =',p
                             endif

                             ! Check that the line is long enough
                             ! to be written, otherwise, skip it.
                             if ( (p_end-p_beg+1) .ge. lines_min_npts ) then
                                ! Line is long enough.

                                ! Augment counter of number of lines.
                                y_j = y_j + 1

                                ! Write this line to output
                                lag_writest2d(2) = y_j

                                !write(*,*) '--------------------------------'
                                !write(*,*) 'Writing in line #',y
                                !write(*,*) 'from ',p_beg,' to ',p_end
                                call put_lag_var_seg_pad(ncid_j,lamvnam,lam,npts,p_beg,p_end,lam_mask)
                                call put_lag_var_seg_pad(ncid_j,phivnam,phi,npts,p_beg,p_end,phi_mask)
                                call put_lag_var_seg_pad(ncid_j,depvnam,dep,npts,p_beg,p_end,dep_mask)
                                ! Lines do not need age since the index on lag vars is age.
                                if (lu) call put_lag_var_seg_pad(ncid_j,uvnam,u,npts,p_beg,p_end,u_mask)
                                if (lv) call put_lag_var_seg_pad(ncid_j,vvnam,v,npts,p_beg,p_end,v_mask)
                                if (lw) call put_lag_var_seg_pad(ncid_j,wvnam,w,npts,p_beg,p_end,w_mask)
                                if (l_stamp_id) then
                                   ! Get stamp variable ID and then write ID
                                   vid = ncvid(ncid_j,id_vnam,exitcode)
                                   call ncvpt1(ncid_j,vid,y_j,current_id,exitcode)
                                endif
                                if (lt) call put_lag_var_seg_pad(ncid_j,tempvnam,temp,npts,p_beg,p_end,t_mask)
                                if (ls) call put_lag_var_seg_pad(ncid_j,saltvnam,salt,npts,p_beg,p_end,s_mask)
                                if (lr) call put_lag_var_seg_pad(ncid_j,rhovnam,rho,npts,p_beg,p_end,r_mask)
                                if (lq) call put_lag_var_seg_pad(ncid_j,qvnam,q,npts,p_beg,p_end,q_mask)
                                if (la) call put_lag_var_seg_pad(ncid_j,dtdz_vnam,dtdz,npts,p_beg,p_end,dtdz_mask)
                                if (lb) call put_lag_var_seg_pad(ncid_j,drdz_vnam,drdz,npts,p_beg,p_end,drdz_mask)
                                if (lc) then
                                   call put_lag_var_seg_pad(ncid_j,year_vnam,year,npts,p_beg,p_end,year_mask)
                                   call put_lag_var_seg_pad(ncid_j,month_vnam,month,npts,p_beg,p_end,month_mask)
                                   call put_lag_var_seg_pad(ncid_j,day_vnam,day,npts,p_beg,p_end,day_mask)
                                endif
                                if (llayer) then
                                   call put_lag_var_seg_pad(ncid_j,layer_top_vnam,top,npts,p_beg,p_end,layer_top_mask)
                                   call put_lag_var_seg_pad(ncid_j,layer_bot_vnam,bot,npts,p_beg,p_end,layer_bot_mask)
                                endif
                                if (lh) call put_lag_var_seg_pad(ncid_j,bdepvnam,bdep,npts,p_beg,p_end,bdep_mask)
                             endif  ! End of long line check
                          else
                             age_in(p) = age_in(p-1) + 1
                          endif

                          ! Since we are not exiting or entering,
                          ! we do not set p_beg here.

                       else
                          !=============ENTRY DETECTED===============

                          ! Reset the age_in to entry point,
                          ! regardless of the sorting type.
                          ! (Since age_in for points outside
                          ! is zero, this should be the same
                          ! as age_in(p) = age_in(p-1) + 1.
                          ! However, for clarity, just this:
                          age_in(p) = 1

                          if ( lwormsin ) then

                             ! Augment counter of number of worms.
                             y = y + 1
                             subseg_per_traj = subseg_per_traj + 1

                             ! Write this worm to output
                             lag_writest2d(2) = y

                             !write(*,*) '--------------------------------'
                             !write(*,*) 'Writing in worm #',y
                             !write(*,*) 'from ',(p-nw),' to ',(p+nw-1)
                             call put_lag_var_seg(ncid_y,lamvnam,lam,npts,(p-nw),(p+nw-1),lam_mask)
                             call put_lag_var_seg(ncid_y,phivnam,phi,npts,(p-nw),(p+nw-1),phi_mask)
                             call put_lag_var_seg(ncid_y,depvnam,dep,npts,(p-nw),(p+nw-1),dep_mask)
                             call put_lag_var_seg(ncid_y,agevnam,age_out,npts,(p-nw),(p+nw-1),age_mask)
                             if (lu) call put_lag_var_seg(ncid_y,uvnam,u,npts,(p-nw),(p+nw-1),u_mask)
                             if (lv) call put_lag_var_seg(ncid_y,vvnam,v,npts,(p-nw),(p+nw-1),v_mask)
                             if (lw) call put_lag_var_seg(ncid_y,wvnam,w,npts,(p-nw),(p+nw-1),w_mask)
                             if (l_stamp_id) then
                                ! Get stamp variable ID and then write ID to newly created worm.
                                vid = ncvid(ncid_y,id_vnam,exitcode)
                                call ncvpt1(ncid_y,vid,y,current_id,exitcode)
                             endif
                             if (lt) call put_lag_var_seg(ncid_y,tempvnam,temp,npts,(p-nw),(p+nw-1),t_mask)
                             if (ls) call put_lag_var_seg(ncid_y,saltvnam,salt,npts,(p-nw),(p+nw-1),s_mask)
                             if (lr) call put_lag_var_seg(ncid_y,rhovnam,rho,npts,(p-nw),(p+nw-1),r_mask)
                             if (lq) call put_lag_var_seg(ncid_y,qvnam,q,npts,(p-nw),(p+nw-1),q_mask)
                             if (la) call put_lag_var_seg(ncid_y,dtdz_vnam,dtdz,npts,(p-nw),(p+nw-1),dtdz_mask)
                             if (lb) call put_lag_var_seg(ncid_y,drdz_vnam,drdz,npts,(p-nw),(p+nw-1),drdz_mask)
                             if (lc) then
                                call put_lag_var_seg(ncid_y,year_vnam,year,npts,(p-nw),(p+nw-1),year_mask)
                                call put_lag_var_seg(ncid_y,month_vnam,month,npts,(p-nw),(p+nw-1),month_mask)
                                call put_lag_var_seg(ncid_y,day_vnam,day,npts,(p-nw),(p+nw-1),day_mask)
                             endif
                             if (llayer) then
                                call put_lag_var_seg(ncid_y,layer_top_vnam,top,npts,(p-nw),(p+nw-1),layer_top_mask)
                                call put_lag_var_seg(ncid_y,layer_bot_vnam,bot,npts,(p-nw),(p+nw-1),layer_bot_mask)
                             endif
                             if (lh) call put_lag_var_seg(ncid_y,bdepvnam,bdep,npts,(p-nw),(p+nw-1),bdep_mask)
                          endif

                          if ( llinesout ) then
                             ! We just entered the domain, so now
                             ! we know the end of a line that was
                             ! once outside the domain.
                             p_end = p - 1

                             ! The beginning of this line was noted
                             ! earlier in p_beg.  However, check for
                             ! consistency with the age vector:
                             if ( age_out(p_beg) .ne. 1 ) then
                                write(*,*) 'ERROR: age_out(p_beg) is not 1, rather =',age_out(p_beg)
                                write(*,*) 'p_beg =',p_beg
                                write(*,*) 'p_end =',p_end
                                write(*,*) 'traj  =',t
                                write(*,*) 'point =',p
                             endif

                             ! Check that lines are long enough,
                             ! otherwise skip line and do not output.
                             if ( (p_end-p_beg+1) .ge. lines_min_npts ) then
                                ! Long enough!
                                ! Augment counter of number of lines.
                                y = y + 1

                                ! Write this line to output
                                lag_writest2d(2) = y

                                !write(*,*) '--------------------------------'
                                !write(*,*) 'Writing out line #',y
                                !write(*,*) 'from ',p_beg,' to ',p_end
                                call put_lag_var_seg_pad(ncid_y,lamvnam,lam,npts,p_beg,p_end,lam_mask)
                                call put_lag_var_seg_pad(ncid_y,phivnam,phi,npts,p_beg,p_end,phi_mask)
                                call put_lag_var_seg_pad(ncid_y,depvnam,dep,npts,p_beg,p_end,dep_mask)
                                ! Lines do not need age since traj age is index value on the lag_vars
                                if (lu) call put_lag_var_seg_pad(ncid_y,uvnam,u,npts,p_beg,p_end,u_mask)
                                if (lv) call put_lag_var_seg_pad(ncid_y,vvnam,v,npts,p_beg,p_end,v_mask)
                                if (lw) call put_lag_var_seg_pad(ncid_y,wvnam,w,npts,p_beg,p_end,w_mask)
                                if (l_stamp_id) then
                                   ! Get stamp variable ID and then write ID to newly created line.
                                   vid = ncvid(ncid_y,id_vnam,exitcode)
                                   call ncvpt1(ncid_y,vid,y,current_id,exitcode)
                                endif
                                if (lt) call put_lag_var_seg_pad(ncid_y,tempvnam,temp,npts,p_beg,p_end,t_mask)
                                if (ls) call put_lag_var_seg_pad(ncid_y,saltvnam,salt,npts,p_beg,p_end,s_mask)
                                if (lr) call put_lag_var_seg_pad(ncid_y,rhovnam,rho,npts,p_beg,p_end,r_mask)
                                if (lq) call put_lag_var_seg_pad(ncid_y,qvnam,q,npts,p_beg,p_end,q_mask)
                                if (la) call put_lag_var_seg_pad(ncid_y,dtdz_vnam,dtdz,npts,p_beg,p_end,dtdz_mask)
                                if (lb) call put_lag_var_seg_pad(ncid_y,drdz_vnam,drdz,npts,p_beg,p_end,drdz_mask)
                                if (lc) then
                                   call put_lag_var_seg_pad(ncid_y,year_vnam,year,npts,p_beg,p_end,year_mask)
                                   call put_lag_var_seg_pad(ncid_y,month_vnam,month,npts,p_beg,p_end,month_mask)
                                   call put_lag_var_seg_pad(ncid_y,day_vnam,day,npts,p_beg,p_end,day_mask)
                                endif
                                if (llayer) then
                                   call put_lag_var_seg_pad(ncid_y,layer_top_vnam,top,npts,p_beg,p_end,layer_top_mask)
                                   call put_lag_var_seg_pad(ncid_y,layer_bot_vnam,bot,npts,p_beg,p_end,layer_bot_mask)
                                endif
                                if (lh) call put_lag_var_seg_pad(ncid_y,bdepvnam,bdep,npts,p_beg,p_end,bdep_mask)
                             endif  ! End of long line check
                          endif

                          if ( llinesin ) then
                             ! We just entered the domain, so
                             ! here we set the beginning of a
                             ! line in.  We must wait until the
                             ! detection of the end of this traj
                             ! at a later point.
                             p_beg = p
                             subseg_per_traj = subseg_per_traj + 1
                          endif
                       endif
                    else
                       ! We are out of the box.
                       inbox_now = .false.

                       if ( p .eq. 1 ) then
                          ! If it is the first point, set equal
                          ! previous time step flag too.
                          inbox_old = inbox_now

                          ! If we are searching for linesout,
                          ! then this point, which is
                          ! outside of the domain, is the
                          ! beginning of the first lineout.
                          if ( llinesout ) then
                             p_beg = 1
                             subseg_per_traj = subseg_per_traj + 1
                          endif
                       endif

                       if ( inbox_now .eqv. inbox_old ) then
                          ! Current and old positions are
                          ! the same.  No exiting detected,
                          ! so skip to the next point.

                          ! Augment age out of box.
                          if ( p .eq. 1 ) then
                             age_out(p) = 1
                          elseif ( (p .eq. npts) .and. lj .and. llinesout ) then
                             ! This is an indefinite termination
                             ! of a trajectory outside of the box
                             ! because the traj ends, not because
                             ! of actual exiting.  If the user 
                             ! specified -J, then record this 
                             ! line out.
                             p_end = npts
                             age_out(p) = age_out(p-1) + 1

                             ! The beginning of this line was noted
                             ! earlier in p_beg.  However, check for
                             ! consistency with the age vector:
                             if ( age_out(p_beg) .ne. 1 ) then
                                write(*,*) 'ERROR: age_out(p_beg) is not 1, rather =',age_out(p_beg)
                                write(*,*) 'p_beg =',p_beg
                                write(*,*) 'p_end =',p_end
                                write(*,*) 'traj  =',t
                                write(*,*) 'point =',p
                             endif

                             ! Check that the line is long enough
                             ! to be written, otherwise, skip it.
                             if ( (p_end-p_beg+1) .ge. lines_min_npts ) then
                                ! Line is long enough.

                                ! Augment counter of number of lines.
                                y_j = y_j + 1

                                ! Write this line to output
                                lag_writest2d(2) = y_j

                                !write(*,*) '--------------------------------'
                                !write(*,*) 'Writing in line #',y
                                !write(*,*) 'from ',p_beg,' to ',p_end
                                call put_lag_var_seg_pad(ncid_j,lamvnam,lam,npts,p_beg,p_end,lam_mask)
                                call put_lag_var_seg_pad(ncid_j,phivnam,phi,npts,p_beg,p_end,phi_mask)
                                call put_lag_var_seg_pad(ncid_j,depvnam,dep,npts,p_beg,p_end,dep_mask)
                                ! Lines do not need age since the index on lag vars is age.
                                if (lu) call put_lag_var_seg_pad(ncid_j,uvnam,u,npts,p_beg,p_end,u_mask)
                                if (lv) call put_lag_var_seg_pad(ncid_j,vvnam,v,npts,p_beg,p_end,v_mask)
                                if (lw) call put_lag_var_seg_pad(ncid_j,wvnam,w,npts,p_beg,p_end,w_mask)
                                if (l_stamp_id) then
                                   ! Get stamp variable ID and then write ID to newly created line.
                                   vid = ncvid(ncid_j,id_vnam,exitcode)
                                   call ncvpt1(ncid_j,vid,y_j,current_id,exitcode)
                                endif
                                if (lt) call put_lag_var_seg_pad(ncid_j,tempvnam,temp,npts,p_beg,p_end,t_mask)
                                if (ls) call put_lag_var_seg_pad(ncid_j,saltvnam,salt,npts,p_beg,p_end,s_mask)
                                if (lr) call put_lag_var_seg_pad(ncid_j,rhovnam,rho,npts,p_beg,p_end,r_mask)
                                if (lq) call put_lag_var_seg_pad(ncid_j,qvnam,q,npts,p_beg,p_end,q_mask)
                                if (la) call put_lag_var_seg_pad(ncid_j,dtdz_vnam,dtdz,npts,p_beg,p_end,dtdz_mask)
                                if (lb) call put_lag_var_seg_pad(ncid_j,drdz_vnam,drdz,npts,p_beg,p_end,drdz_mask)
                                if (lc) then
                                   call put_lag_var_seg_pad(ncid_j,year_vnam,year,npts,p_beg,p_end,year_mask)
                                   call put_lag_var_seg_pad(ncid_j,month_vnam,month,npts,p_beg,p_end,month_mask)
                                   call put_lag_var_seg_pad(ncid_j,day_vnam,day,npts,p_beg,p_end,day_mask)
                                endif
                                if (llayer) then
                                   call put_lag_var_seg_pad(ncid_j,layer_top_vnam,top,npts,p_beg,p_end,layer_top_mask)
                                   call put_lag_var_seg_pad(ncid_j,layer_bot_vnam,bot,npts,p_beg,p_end,layer_bot_mask)
                                endif
                                if (lh) call put_lag_var_seg_pad(ncid_j,bdepvnam,bdep,npts,p_beg,p_end,bdep_mask)
                             endif  ! End of long line check
                          else
                             age_out(p) = age_out(p-1) + 1
                          endif

                          ! Since we are not exiting or entering,
                          ! we do not set p_beg here.

                       else
                          !=============EXIT DETECTED=================

                          ! This is the exit point, so
                          ! set the age_out, regardless
                          ! of the sorting type.
                          age_out(p) = 1

                          if ( lwormsout ) then             

                             ! Write to wormsout

                             ! Augment counter of number of worms.
                             y = y + 1
                             subseg_per_traj = subseg_per_traj + 1

                             ! Write this worm to output
                             lag_writest2d(2) = y

                             !write(*,*) '--------------------------------'
                             !write(*,*) 'Writing out worm #',y
                             !write(*,*) 'from ',(p-nw),' to ',(p+nw-1)
                             call put_lag_var_seg(ncid_y,lamvnam,lam,npts,(p-nw),(p+nw-1),lam_mask)
                             call put_lag_var_seg(ncid_y,phivnam,phi,npts,(p-nw),(p+nw-1),phi_mask)
                             call put_lag_var_seg(ncid_y,depvnam,dep,npts,(p-nw),(p+nw-1),dep_mask)
                             call put_lag_var_seg(ncid_y,agevnam,age_in,npts,(p-nw),(p+nw-1),age_mask)
                             if (lu) call put_lag_var_seg(ncid_y,uvnam,u,npts,(p-nw),(p+nw-1),u_mask)
                             if (lv) call put_lag_var_seg(ncid_y,vvnam,v,npts,(p-nw),(p+nw-1),v_mask)
                             if (lw) call put_lag_var_seg(ncid_y,wvnam,w,npts,(p-nw),(p+nw-1),w_mask)
                             if (l_stamp_id) then
                                ! Get stamp variable ID and then write ID to newly created line.
                                vid = ncvid(ncid_y,id_vnam,exitcode)
                                call ncvpt1(ncid_y,vid,y,current_id,exitcode)
                             endif
                             if (lt) call put_lag_var_seg(ncid_y,tempvnam,temp,npts,(p-nw),(p+nw-1),t_mask)
                             if (ls) call put_lag_var_seg(ncid_y,saltvnam,salt,npts,(p-nw),(p+nw-1),s_mask)
                             if (lr) call put_lag_var_seg(ncid_y,rhovnam,rho,npts,(p-nw),(p+nw-1),r_mask)
                             if (lq) call put_lag_var_seg(ncid_y,qvnam,q,npts,(p-nw),(p+nw-1),q_mask)
                             if (la) call put_lag_var_seg(ncid_y,dtdz_vnam,dtdz,npts,(p-nw),(p+nw-1),dtdz_mask)
                             if (lb) call put_lag_var_seg(ncid_y,drdz_vnam,drdz,npts,(p-nw),(p+nw-1),drdz_mask)
                             if (lc) then
                                call put_lag_var_seg(ncid_y,year_vnam,year,npts,(p-nw),(p+nw-1),year_mask)
                                call put_lag_var_seg(ncid_y,month_vnam,month,npts,(p-nw),(p+nw-1),month_mask)
                                call put_lag_var_seg(ncid_y,day_vnam,day,npts,(p-nw),(p+nw-1),day_mask)
                             endif
                             if (llayer) then
                                call put_lag_var_seg(ncid_y,layer_top_vnam,top,npts,(p-nw),(p+nw-1),layer_top_mask)
                                call put_lag_var_seg(ncid_y,layer_bot_vnam,bot,npts,(p-nw),(p+nw-1),layer_bot_mask)
                             endif
                             if (lh) call put_lag_var_seg(ncid_y,bdepvnam,bdep,npts,(p-nw),(p+nw-1),bdep_mask)
                          endif

                          if ( llinesout ) then
                             ! We just exited the domain, so
                             ! here we set the beginning of a
                             ! lineOUT.  We must wait until the
                             ! detection of the end of this line
                             ! at a later point.
                             p_beg = p
                             subseg_per_traj = subseg_per_traj + 1
                          endif

                          if ( llinesin ) then
                             ! We have now found the end of the
                             ! current lineIN.
                             p_end = p - 1

                             ! The beginning of this line was noted
                             ! earlier in p_beg.  However, check for
                             ! consistency with the age vector:
                             if ( age_in(p_beg) .ne. 1 ) then
                                write(*,*) 'ERROR: age_in(p_beg) is not 1, rather =',age_in(p_beg)
                                write(*,*) 'p_beg =',p_beg
                                write(*,*) 'p_end =',p_end
                                write(*,*) 'traj  =',t
                                write(*,*) 'point =',p
                             endif

                             ! Check that the line is long enough
                             ! to be written, otherwise, skip it.
                             if ( (p_end-p_beg+1) .ge. lines_min_npts ) then
                                ! Line is long enough.

                                ! Augment counter of number of lines.
                                y = y + 1

                                ! Write this line to output
                                lag_writest2d(2) = y

                                !write(*,*) '--------------------------------'
                                !write(*,*) 'Writing in line #',y
                                !write(*,*) 'from ',p_beg,' to ',p_end
                                call put_lag_var_seg_pad(ncid_y,lamvnam,lam,npts,p_beg,p_end,lam_mask)
                                call put_lag_var_seg_pad(ncid_y,phivnam,phi,npts,p_beg,p_end,phi_mask)
                                call put_lag_var_seg_pad(ncid_y,depvnam,dep,npts,p_beg,p_end,dep_mask)
                                ! Lines do not need age since the index on lag vars is age.
                                if (lu) call put_lag_var_seg_pad(ncid_y,uvnam,u,npts,p_beg,p_end,u_mask)
                                if (lv) call put_lag_var_seg_pad(ncid_y,vvnam,v,npts,p_beg,p_end,v_mask)
                                if (lw) call put_lag_var_seg_pad(ncid_y,wvnam,w,npts,p_beg,p_end,w_mask)
                                if (l_stamp_id) then
                                   ! Get stamp variable ID and then write ID to newly created line.
                                   vid = ncvid(ncid_y,id_vnam,exitcode)
                                   call ncvpt1(ncid_y,vid,y,current_id,exitcode)
                                endif
                                if (lt) call put_lag_var_seg_pad(ncid_y,tempvnam,temp,npts,p_beg,p_end,t_mask)
                                if (ls) call put_lag_var_seg_pad(ncid_y,saltvnam,salt,npts,p_beg,p_end,s_mask)
                                if (lr) call put_lag_var_seg_pad(ncid_y,rhovnam,rho,npts,p_beg,p_end,r_mask)
                                if (lq) call put_lag_var_seg_pad(ncid_y,qvnam,q,npts,p_beg,p_end,q_mask)
                                if (la) call put_lag_var_seg_pad(ncid_y,dtdz_vnam,dtdz,npts,p_beg,p_end,dtdz_mask)
                                if (lb) call put_lag_var_seg_pad(ncid_y,drdz_vnam,drdz,npts,p_beg,p_end,drdz_mask)
                                if (lc) then
                                   call put_lag_var_seg_pad(ncid_y,year_vnam,year,npts,p_beg,p_end,year_mask)
                                   call put_lag_var_seg_pad(ncid_y,month_vnam,month,npts,p_beg,p_end,month_mask)
                                   call put_lag_var_seg_pad(ncid_y,day_vnam,day,npts,p_beg,p_end,day_mask)
                                endif
                                if (llayer) then
                                   call put_lag_var_seg_pad(ncid_y,layer_top_vnam,top,npts,p_beg,p_end,layer_top_mask)
                                   call put_lag_var_seg_pad(ncid_y,layer_bot_vnam,bot,npts,p_beg,p_end,layer_bot_mask)
                                endif
                                if (lh) call put_lag_var_seg_pad(ncid_y,bdepvnam,bdep,npts,p_beg,p_end,bdep_mask)
                             endif  ! End of long line check
                          endif
                       endif  ! End of recent exit check

                    endif   ! End of in box check

                    ! Moving onto the next time step,
                    ! so update the inbox flag
                    inbox_old = inbox_now

                 enddo   ! End of looping over traj points

                 ! Check that the age_in and age_out
                 ! vectors are consistent with each
                 ! other.  Cannot be zero at the same
                 ! time because the float is either in
                 ! or out of the domain.
                 do p = 1,npts
                    if ( (age_in(p) .eq. 0) .and. (age_out(p) .eq. 0) ) then
                       write(*,*) 'ERROR: age_in and age_out mismatch!'
                       write(*,*) 'age_in =',age_in(p)
                       write(*,*) 'age_out=',age_out(p)
                       write(*,*) 'traj   =',t
                       write(*,*) 'point  =',p
                    endif
                    if ( (age_in(p) .gt. 0) .and. (age_out(p) .gt. 0) ) then
                       write(*,*) 'ERROR: age_in and age_out mismatch!'
                       write(*,*) 'age_in =',age_in(p)
                       write(*,*) 'age_out=',age_out(p)
                       write(*,*) 'traj   =',t
                       write(*,*) 'point  =',p
                    endif
                 enddo

                 if (ldump) write(tdo,'(i10,x,i10)') t,subseg_per_traj

              else
                 write(*,*) ' Unexpected operation code!  Cannot sort.'
                 stop
              endif  ! End of start, through, fastest, or exit check
           enddo   ! End of looping over traj

           ! Compute the number of yes traj
           ! for only the last exit case because
           ! we could not augment the counter
           ! in the loop.
           if ( llexit ) then
              yntraj = 0
              do t = 1,ntraj
                 yntraj = yntraj + inyes(t)
              enddo
           endif
        endif   ! End of number check
        
        !==================================================
        !               POST-SORTING WRITING
        !==================================================
        !
        ! For the modes listed below, write to output
        ! after sorting.  This is not strictly necessary
        ! for any of the modes except fastest, but is
        ! retained here because of the orginal structure
        ! of the code.  Room for making things more
        ! efficient here b/c traj are read a second time
        ! before being written.
        !
        ! The worms and lines modes execute all writing
        ! from within the sorting loops.  The start,
        ! through, fexit and lexit modes could do this as
        ! well, but they are not set up to do this now.
        !
        !--------------------------------------------------
        if ( lstart .or. lthrough .or. lfastest .or.&
             lfexit .or. llexit .or. lnumber .or. lend ) then

           nntraj = ntraj - yntraj
           
           ! Double check
           if ( nntraj .lt. 0 ) then
              write(*,*) 'Error: Found too many yes traj!'
              stop
           endif

           ! If we are running in search for fastest trajectory
           ! mode, modify the inyes vector, and yntraj, and
           ! nntraj accordingly.  Also, kill out any points
           ! that are beyond the arrival point.
           if ( lfastest ) then

              ! Allocate space for a vector to hold the nbest
              ! (traj num of the n fastest trajectories).
              allocate(nbest(nf))
              
              ! Initial nbest with the first nf trajectories.
              do t = 1,nf
                 nbest(t) = t
              enddo

              ! For all the remaining trajectories,
              do t = (nf+1),ntraj
                 
                 ! Find the slowest traj currently in
                 ! nbest.
                 nslow = 1
                 do p = 1,nf
                    if ( n2box(nbest(p)) .gt. n2box(nbest(nslow)) ) then
                       nslow = p
                    endif
                 enddo
                 !write(*,*) 'Slowest traj: ',nslow,' w/ len: ',n2box(nbest(nslow))

                 ! Test if the current trajectory
                 ! is faster than this slowest traj.
                 if ( n2box(nbest(nslow)) .gt. n2box(t) ) then
                    ! If the current traj is faster,
                    ! replace the slowest traj in nbest
                    ! with the current trajectory.
                    nbest(nslow) = t
                    !write(*,*) 'Faster traj: ',n2box(t)
                 endif
              enddo

              ! Update inyes with the trajectory numbers in
              ! nbest only.  Then set yntraj and nntraj to
              ! the corresponding values.
              inyes = 0
              do t = 1,nf
                 inyes(nbest(t)) = 1
              enddo
              
              yntraj = nf
              nntraj = ntraj - nf
              
           endif
        
           write(*,*) ' Started with ',ntraj,' trajectories.'
           write(*,*) ' Found ',yntraj, ' trajectories for yes.nc'
           write(*,*) ' Found ',nntraj, ' trajectories for no.nc'
           
           !------Perform the random culling of traj here------
           if (lp) then
              write(*,*) ' Randomly selecting ',percent_keep,' percent of traj.'
              
              ! Cycle through each of the floats
              do t = 1,ntraj

                 if ( inyes(t) .eq. 1 ) then
                    ! We test whether or not to throw this float out.

                    ! Get a random number [0,1]
                    call random_number(current_random)
                    
                    if ( current_random(1) .gt. percent_keep ) then
                       ! We remove this float from inyes
                       ! by modifing inyes, yntraj, nntraj
                       inyes(t) = 0
                       yntraj = yntraj - 1
                       nntraj = nntraj + 1
                    endif
                 endif
              enddo

              write(*,*) ' Retained ',yntraj,' trajectories for yes.nc.'
              write(*,*) ' Remaining ',nntraj,' trajectories for no.nc.'
           endif

           !---------------------------------------------------
           if ( .not. lnone ) then
              write(*,*) ' Writing traj to output files...'
              y = 1
              n = 1
              do t = 1,ntraj
              
                 !write(*,*) 'For traj ',t,' out of ',ntraj
           
                 ! Load traj data again
                 !write(*,*) 'Specify the read vectors...'
                 lag_readst2d(2) = t
           
                 !write(*,*) 'Load traj data...'
                 !--------------------------------------------------------
                 ! Load position variables
                 call get_lag_var(ncid_i,lamvnam,lam)
                 call get_lag_var(ncid_i,phivnam,phi)
                 call get_lag_var(ncid_i,depvnam,dep)
                 !--------------------------------------------------------
                 ! Load velocity variables
                 if (lu .and. (.not. (lcent .or. lfore .or. lback))) then
                    call get_lag_var(ncid_i,uvnam,u)
                 endif
                 if (lv .and. (.not. (lcent .or. lfore .or. lback))) then
                    call get_lag_var(ncid_i,vvnam,v)
                 endif
                 if (lw .and. (.not. (lcent .or. lfore .or. lback))) then
                    call get_lag_var(ncid_i,wvnam,w)
                 endif
                 !--------------------------------------------------------
                 ! Compute velocity from the position information with
                 ! choice of three different differencing schemes.
                 if (lcent) then
                    do p = 2,npts-1
                       ! Use the avg. latitude between the centered
                       ! difference endpoints to determine the latitude
                       ! spacing degrees to meter conversion.  This is
                       ! instead of the center point because the center
                       ! point is not necessarily in the center of the
                       ! two centered difference endpoints.
                       !        change in lon      to m        account for sphere
                       u(p) = ((lam(p+1)-lam(p-1))*deg2m*cos(((phi(p+1)+phi(p-1))/2.0)*deg2rad))/(2.0*tstep)
                       v(p) = ((phi(p+1)-phi(p-1))*deg2m)/(2.0*tstep)
                       w(p) = (dep(p+1)-dep(p-1))/(2.0*tstep)
                    enddo
                    ! Endpoints
                    u(1) = u(2)
                    v(1) = v(2)
                    w(1) = w(2)
                    u(npts) = u(npts-1)
                    v(npts) = v(npts-1)
                    w(npts) = w(npts-1)
                 elseif (lfore) then
                    do p = 1,npts-1
                       u(p) = ((lam(p+1)-lam(p))*deg2m*cos(((phi(p+1)+phi(p))/2.0)*deg2rad))/tstep
                       v(p) = ((phi(p+1)-phi(p))*deg2m)/tstep
                       w(p) = (dep(p+1)-dep(p))/tstep
                    enddo
                    u(npts) = u(npts-1)
                    v(npts) = v(npts-1)
                    w(npts) = w(npts-1)
                 elseif (lback) then
                    do p = 2,npts
                       u(p) = ((lam(p)-lam(p-1))*deg2m*cos(((phi(p)+phi(p-1))/2.0)*deg2rad))/tstep
                       v(p) = ((phi(p)-phi(p-1))*deg2m)/tstep
                       w(p) = (dep(p)-dep(p-1))/tstep
                    enddo
                    u(1) = u(2)
                    v(1) = v(2)
                    w(1) = w(2)
                 endif
                 !--------------------------------------------------------
                 if ( l_keep_id ) then
                    ! Get the vid of the id variable
                    vid = ncvid(ncid_i,id_vnam,exitcode)

                    ! Load the ID value into current_id
                    call ncvgt1(ncid_i,vid,t,current_id,exitcode)
                 else
                    current_id = t
                 endif
                 ! Load tracer variables
                 if (lt) call get_lag_var(ncid_i,tempvnam,temp)
                 if (ls) call get_lag_var(ncid_i,saltvnam,salt)
                 if (lr) call get_lag_var(ncid_i,rhovnam,rho)
                 if (lq) call get_lag_var(ncid_i,qvnam,q)
                 if (la) call get_lag_var(ncid_i,dtdz_vnam,dtdz)
                 if (lb) call get_lag_var(ncid_i,drdz_vnam,drdz)
                 if (le .and. (.not. lage_const) ) call get_lag_var(ncid_i,agevnam,age)
                 if (lc) then
                    call get_lag_var(ncid_i,year_vnam,year)
                    call get_lag_var(ncid_i,month_vnam,month)
                    call get_lag_var(ncid_i,day_vnam,day)
                 endif
                 if (llayer) then
                    call get_lag_var(ncid_i,layer_top_vnam,top)
                    call get_lag_var(ncid_i,layer_bot_vnam,bot)
                 endif
                 if (lh) call get_lag_var(ncid_i,bdepvnam,bdep)

                 if ( inyes(t) .eq. 1 ) then

                    ! Check for overrun
                    if ( y .gt. yntraj ) then
                       write(*,*) ' Error: Too many yes traj.'
                       stop
                    endif
                    
                    if ( (lfastest .or. lthrough) .and. lmask ) then
                       ! For all trajectories, mask out any
                       ! points after n2box.
                       do p = 1,npts
                          if ( p .gt. n2box(t) ) then
                             lam(p) = lam_mask
                             phi(p) = phi_mask
                             dep(p) = dep_mask
                             u(p) = u_mask
                             v(p) = v_mask
                             w(p) = w_mask
                             temp(p) = t_mask
                             salt(p) = s_mask
                             rho(p) = r_mask
                             q(p) = q_mask
                             dtdz(p) = dtdz_mask
                             drdz(p) = drdz_mask
                             age(p) = age_mask
                             bdep(p) = bdep_mask
                             if (lc) then
                                year(p) = year_mask
                                month(p) = month_mask
                                day(p) = day_mask
                             endif
                             if (llayer) then
                                top(p) = layer_top_mask
                                bot(p) = layer_bot_mask
                             endif
                          endif
                       enddo
                    endif

                    ! For the case of exit mode, mask
                    ! out all traj which are before the
                    ! last exiting of the box for each
                    ! particle by making the last exiting
                    ! event the first point in the traj.
                    ! This will reset the age of the traj.
                    if ( llexit .or. lfexit ) then
                       if ( lmask ) then
                          ! MASKING = KEEP THE INITIAL POSITIONS
                          ! and mask out any points after n2box.
                          do p = 1,npts
                             if ( p .gt. n2box(t) ) then
                                lam(p) = lam_mask
                                phi(p) = phi_mask
                                dep(p) = dep_mask
                                u(p) = u_mask
                                v(p) = v_mask
                                w(p) = w_mask
                                temp(p) = t_mask
                                salt(p) = s_mask
                                rho(p) = r_mask
                                q(p) = q_mask
                                dtdz(p) = dtdz_mask
                                drdz(p) = drdz_mask
                                age(p) = age_mask
                                bdep(p) = bdep_mask
                                if (lc) then
                                   year(p) = year_mask
                                   month(p) = month_mask
                                   day(p) = day_mask
                                endif
                                if (llayer) then
                                   top(p) = layer_top_mask
                                   bot(p) = layer_bot_mask
                                endif
                             endif
                          enddo
                       else
                          ! NO MASKING = KEEP THE EXIT POSITIONS
                          ! If the last point in the box
                          ! is n2box, and we wish to trim
                          ! all traj to this point in time,
                          ! then the length of the traj
                          ! will become npts-n2box+1.
                          
                          ! It is possible for n2box to be
                          ! zero if the particle was not in
                          ! the box to begin with. To avoid
                          ! a segmentation fault, set n2box
                          ! to the first point in that case
                          ! so the whole out trajectory is
                          ! written to output.
                          !=================================
                          ! WORKING HERE! NEED TO TEST THIS
                          !=================================
                          if (n2box(t) .eq. 0 ) n2box(t) = 1
                          
                          do p = 1,(npts-n2box(t)+1)
                             
                             pp = p + n2box(t) - 1
                             
                             lam(p) = lam(pp)
                             phi(p) = phi(pp)
                             dep(p) = dep(pp)
                             u(p) = u(pp)
                             v(p) = v(pp)
                             w(p) = w(pp)
                             temp(p) = temp(pp)
                             salt(p) = salt(pp)
                             rho(p) = rho(pp)
                             q(p) = q(pp)
                             dtdz(p) = dtdz(pp)
                             drdz(p) = drdz(pp)
                             age(p) = age(pp)
                             bdep(p) = bdep(pp)
                             if (lc) then
                                year(p) = year(pp)
                                month(p) = month(pp)
                                day(p) = day(pp)
                             endif
                             if (llayer) then
                                top(p) = top(pp)
                                bot(p) = bot(pp)
                             endif
                          enddo
                          ! Mask out the tail end of
                          ! the trajectory after the
                          ! useful data has already been
                          ! shifted.
                          do p = (npts-n2box(t)+2),npts
                             lam(p) = lam_mask
                             phi(p) = phi_mask
                             dep(p) = dep_mask
                             u(p) = u_mask
                             v(p) = v_mask
                             w(p) = w_mask
                             temp(p) = t_mask
                             salt(p) = s_mask
                             rho(p) = r_mask
                             q(p) = q_mask
                             dtdz(p) = dtdz_mask
                             drdz(p) = drdz_mask
                             age(p) = age_mask
                             bdep(p) = bdep_mask
                             if (lc) then
                                year(p) = year_mask
                                month(p) = month_mask
                                day(p) = day_mask
                             endif
                             if (llayer) then
                                top(p) = layer_top_mask
                                bot(p) = layer_bot_mask
                             endif
                          enddo
                       endif
                    endif

                    !write(*,*) 'Specify write vectors...'
                    lag_writest2d(2) = y

                    !write(*,*) 'Write traj #',y,' to yes.nc...'
                    call put_lag_var(ncid_y,lamvnam,lam,npts)
                    call put_lag_var(ncid_y,phivnam,phi,npts)
                    call put_lag_var(ncid_y,depvnam,dep,npts)
                    if (lu) call put_lag_var(ncid_y,uvnam,u,npts)
                    if (lv) call put_lag_var(ncid_y,vvnam,v,npts)
                    if (lw) call put_lag_var(ncid_y,wvnam,w,npts)
                    if (l_stamp_id .or. l_keep_id) then
                       ! Get stamp variable ID and then write ID
                       vid = ncvid(ncid_y,id_vnam,exitcode)
                       call ncvpt1(ncid_y,vid,y,current_id,exitcode)
                    endif
                    if (lt) call put_lag_var(ncid_y,tempvnam,temp,npts)
                    if (ls) call put_lag_var(ncid_y,saltvnam,salt,npts)
                    if (lr) call put_lag_var(ncid_y,rhovnam,rho,npts)
                    if (lq) call put_lag_var(ncid_y,qvnam,q,npts)
                    if (la) call put_lag_var(ncid_y,dtdz_vnam,dtdz,npts)
                    if (lb) call put_lag_var(ncid_y,drdz_vnam,drdz,npts)
                    if (le .and. (.not. lage_const) ) call put_lag_var(ncid_y,agevnam,age,npts)
                    if (lc) then
                       call put_lag_var(ncid_y,year_vnam,year,npts)
                       call put_lag_var(ncid_y,month_vnam,month,npts)
                       call put_lag_var(ncid_y,day_vnam,day,npts)
                    endif
                    if (llayer) then
                       call put_lag_var(ncid_y,layer_top_vnam,top,npts)
                       call put_lag_var(ncid_y,layer_bot_vnam,bot,npts)
                    endif
                    if (lh) call put_lag_var(ncid_y,bdepvnam,bdep,npts)
                    y = y + 1
                    
                 elseif ( inyes(t) .eq. 0 ) then
                    
                    ! Check for overrun
                    if ( n .gt. nntraj ) then
                       write(*,*) ' Error: Too many no traj.'
                       stop
                    endif
                    
                    if (lno .and. ( lstart .or. lthrough .or. lfastest .or.&
                         lfexit .or. llexit .or. lnumber .or. lend ) ) then
                       !write(*,*) 'Specify write vectors...'
                       lag_writest2d(2) = n
                       
                       !write(*,*) 'Write data to no.nc...'
                       call put_lag_var(ncid_n,lamvnam,lam,npts)
                       call put_lag_var(ncid_n,phivnam,phi,npts)
                       call put_lag_var(ncid_n,depvnam,dep,npts)
                       if (lu) call put_lag_var(ncid_n,uvnam,u,npts)
                       if (lv) call put_lag_var(ncid_n,vvnam,v,npts)
                       if (lw) call put_lag_var(ncid_n,wvnam,w,npts)
                       if (l_stamp_id .or. l_keep_id) then
                          ! Get stamp variable ID and then write ID
                          vid = ncvid(ncid_n,id_vnam,exitcode)
                          call ncvpt1(ncid_n,vid,n,current_id,exitcode)
                       endif
                       if (lt) call put_lag_var(ncid_n,tempvnam,temp,npts)
                       if (ls) call put_lag_var(ncid_n,saltvnam,salt,npts)
                       if (lr) call put_lag_var(ncid_n,rhovnam,rho,npts)
                       if (lq) call put_lag_var(ncid_n,qvnam,q,npts)
                       if (la) call put_lag_var(ncid_n,dtdz_vnam,dtdz,npts)
                       if (lb) call put_lag_var(ncid_n,drdz_vnam,drdz,npts)
                       if (le .and. (.not. lage_const) ) call put_lag_var(ncid_n,agevnam,age,npts)
                       if (lc) then
                          call put_lag_var(ncid_n,year_vnam,year,npts)
                          call put_lag_var(ncid_n,month_vnam,month,npts)
                          call put_lag_var(ncid_n,day_vnam,day,npts)
                       endif
                       if (llayer) then
                          call put_lag_var(ncid_n,layer_top_vnam,top,npts)
                          call put_lag_var(ncid_n,layer_bot_vnam,bot,npts)
                       endif
                       if (lh) call put_lag_var(ncid_n,bdepvnam,bdep,npts)
                    endif
                    n = n + 1
                 else
                    write(*,*) 'Error in inyes array.'
                    stop
                 endif
              enddo
           endif
        endif   ! End of test of post-sorting traj writing

        !==================================================
        !                    CLEAN UP
        !==================================================
        write(*,*) 'Close input file...'
        call ncclos(ncid_i,exitcode)

        if ( .not. lnone ) then
           write(*,*) 'Close yes.nc...'
           call ncclos(ncid_y,exitcode)
        endif

        if (lj) then
           write(*,*) 'Close termination file...'
           call ncclos(ncid_j,exitcode)
        endif
        
        if (lno .and. ( lstart .or. lthrough .or. lfastest .or.&
                        lfexit .or. llexit .or. lnumber .or. lend ) .and.&
                        .not. lnone ) then
           write(*,*) 'Close no.nc...'
           call ncclos(ncid_n,exitcode)
        endif

        write(*,*) 'Deallocate input trajectories...'
        deallocate(lam,phi,dep,temp,salt,rho,u,v,w,q,drdz,dtdz,age,bdep)
        if (lc) deallocate(year,month,day)
        if (llayer) deallocate(top,bot)

        if ( lthrough .or. lfastest .or.&
             llexit .or. lfexit .or.&
             lstart .or. lnumber ) then
           write(*,*) 'Writing inyes to file...'
           ncid_i = 41
           open(unit=ncid_i,file='inyes.txt',status='new',form='formatted')
           do t = 1,ntraj
              write(ncid_i,'(i10)') inyes(t)
           enddo
           close(ncid_i)
        endif

        if( lthrough .or. lfastest .or. llexit .or. lfexit ) then
           write(*,*) 'Writing n2box to file...'
           ncid_i = 43
           open(unit=ncid_i,file='n2box.txt',status='new',form='formatted')
           do t = 1,ntraj
              write(ncid_i,'(i10)') n2box(t)
           enddo
           close(ncid_i)
        endif

        if ( ldump ) close(tdo)

        write(*,*) 'End of tcdfsort.'

      end program tcdfsort

!---------------------------------------
!---------------------------------------

      subroutine print_help()

!---------------------------------------
! This subroutine prints help info
! about this program.
!---------------------------------------

        write(*,*) '-------------------------------------------'
        write(*,*) ' tcdfsort -<start|end|through|fastest N|...'
        write(*,*) '         ...fexit|lexit|wormsin N|wormsout N|...'
        write(*,*) '         ...linesin|linesout|number>'
        write(*,*) '          -I infile.nc'
        write(*,*) '          [-(X,Y,Z,U,V,W,T,S,R,Q,E,A,B,K,d,m,y,H) min max]'
        write(*,*) '          [-P (0.0,1.0)] [-h] [-e] [-N] [-M]'
        write(*,*) '          [-O] [-<F frac|G>] [-L thk] [-C] [-l]'
        write(*,*) '          [-D] [-J] [-t] [-b] [-a] [-<c|f|k> dt]'
        write(*,*) '          [-s] [-i]'
        write(*,*) ' '
        write(*,*) '-------------------------------------------'
        write(*,*) 'This program will sort traj as specified'
        write(*,*) 'by the operation flag and optional sorting'
        write(*,*) 'domain.                         Stefan Gary'
        write(*,*) '-------------------------------------------'
        write(*,*) 'Possible operation flags are:'
        write(*,*) ' -start         sort by traj start position'
        write(*,*) ' '
        write(*,*) ' -end           sort by traj end position'
        write(*,*) '                Note that if the traj is'
        write(*,*) '                masked at the end, it will'
        write(*,*) '                not be written to yes.nc'
        write(*,*) ' '
        write(*,*) ' -through       traj that go through domain'
        write(*,*) ' '
        write(*,*) ' -fastest <N>   N fastest traj to domain'
        write(*,*) ' '
        write(*,*) ' -number <file> traj w/ index # listed in'
        write(*,*) '                the given file.  Any sorting'
        write(*,*) '                domain info is ignored.'
        write(*,*) ' '
        write(*,*) ' -lexit         traj that exit the box will'
        write(*,*) '                be kept.  Only positions after'
        write(*,*) '                the LAST exit attempt will be'
        write(*,*) '                retained and the traj will be'
        write(*,*) '                shifted so the age is reset.'
        write(*,*) '                NOTE: If the -M flag is selected,'
        write(*,*) '                then the opposite happens: all'
        write(*,*) '                traj points before the last'
        write(*,*) '                exit attempt are written out'
        write(*,*) '                and the rest of the traj is masked.'
        write(*,*) ' '
        write(*,*) ' -fexit         traj that exit the box will'
        write(*,*) '                be kept.  Only positions after'
        write(*,*) '                the FIRST exit attempt will be'
        write(*,*) '                retained and the traj will be'
        write(*,*) '                shifted so the age is reset.'
        write(*,*) '                NOTE: If the -M flag is selected,'
        write(*,*) '                then the opposite happens: all'
        write(*,*) '                traj points before the first'
        write(*,*) '                exit attempt are written out'
        write(*,*) '                and the rest of the traj is masked.'
        write(*,*) ' '
        write(*,*) ' -worms<in|out> <N>'
        write(*,*) '                will chop up traj into'
        write(*,*) '                subsegments, "worms," of'
        write(*,*) '                length 2N. The time step'
        write(*,*) '                at N+1 is the point a traj'
        write(*,*) '                is out (for wormsout) or'
        write(*,*) '                is in (for wormsin) the'
        write(*,*) '                user defined domain.  Thus,'
        write(*,*) '                point N is the point just'
        write(*,*) '                before exit (for worms out)'
        write(*,*) '                or entry (for worms in).'
        write(*,*) ' '
        write(*,*) ' -lines<in|out> <N>'
        write(*,*) '                will chop up traj into'
        write(*,*) '                subsegments, "lines," of'
        write(*,*) '                arbitrary minimum length N'
        write(*,*) '                that are outside (linesout),'
        write(*,*) '                inside (linesin).  The age'
        write(*,*) '                will be reset with each line'
        write(*,*) '                and the first point in the'
        write(*,*) '                line will be the last time'
        write(*,*) '                the traj was outside (linesin)'
        write(*,*) '                or inside (linesout) the domain.'
        write(*,*) ' '
        write(*,*) 'One of the operation flags must be specified.'
        write(*,*) 'The operation flag must be the first item on'
        write(*,*) 'the command line.'
        write(*,*) '------------------------------------------------'
        write(*,*) 'The sorting domain is set by the optional'
        write(*,*) 'variable flags (include as many as desired):'
        write(*,*) ' -X = zonal limits, in deg.'
        write(*,*) ' -Y = meridional limits, in deg.'
        write(*,*) ' -Z = depth limits, [m] down'
        write(*,*) ' -U = zonal velocity [m/s]'
        write(*,*) ' -V = meridional velocity [m/s]'
        write(*,*) ' -W = vertical velocity [m/s]'
        write(*,*) ' -T = temperature [oC]'
        write(*,*) ' -S = salinity [PSU]'
        write(*,*) ' -R = density'
        write(*,*) ' -Q = potential vorticity [(ms)^-1]'
        write(*,*) ' -A = vertical grad. of temp [oC m^-1]'
        write(*,*) ' -B = vertical grad. of dens [kg m^-4]'
        write(*,*) ' -E = age of float [time_steps]'
        write(*,*) ' -K = fractional depth within layer'
        write(*,*) ' -d = day (in a month)'
        write(*,*) ' -m = month'
        write(*,*) ' -y = year'
        write(*,*) ' -H = bottom bepth'
        write(*,*) ' '
!        write(*,*) ' NOTE: -K, -y only applies in -start'
!        write(*,*) '       and -end filtering modes!  Even if'
!        write(*,*) '       this flag is specified for the'
!        write(*,*) '       other modes, it will be ignored.'
!        write(*,*) ' '
        write(*,*) 'The variable flags also are used to indicate'
        write(*,*) 'which variables to write to the output file.'
        write(*,*) 'If you wish to write out a variable but not'
        write(*,*) 'use it in the sorting process, append all to'
        write(*,*) 'the flag (e.g. -Tall).'
        write(*,*) ' '
        write(*,*) 'Note that (U,V,W) are, by default, loaded'
        write(*,*) 'from the input file (if not in the file,'
        write(*,*) 'the program will crash).  However, if one'
        write(*,*) 'of -c, -f, or -k is specified, then (U,V,W)'
        write(*,*) 'are computed from the particle positions.'
        write(*,*) ' '
        write(*,*) 'The layer sorting domain flags are similar,'
        write(*,*) 'but owing to the unique behaviors of layers,'
        write(*,*) 'there are some special differences.'
        write(*,*) ' '
        write(*,*) ' -O = the user domain requires that the'
        write(*,*) '      layer be outcropped (at surf) but does'
        write(*,*) '      not require that the particle is in'
        write(*,*) '      the layer.'
        write(*,*) ' '
        write(*,*) ' -L thk = requires a minimum layer thickness'
        write(*,*) '      but does not require particle in layer.'
        write(*,*) '      Thicknesses less than or equal to zero'
        write(*,*) '      are ignored.'
        write(*,*) ' '
        write(*,*) ' -F frac = of all particles launched in'
        write(*,*) '      this column, select the particle that'
        write(*,*) '      is the closest to the frac percentage'
        write(*,*) '      of the layer depth.  E.g. -F 0.25 is'
        write(*,*) '      particle in upper quarter of column.'
        write(*,*) '      Requires particle to be in the layer,'
        write(*,*) '      so automatically sets -t, -b flags.'
        write(*,*) ' '
        write(*,*) ' -G = of all particles launched in this'
        write(*,*) '      column, get the one with the min dtdz.'
        write(*,*) '      Requires particle to be in the layer,'
        write(*,*) '      so automatically sets -t, -b flags.'
        write(*,*) ' '
        write(*,*) ' -b = particle must be at least under the top'
        write(*,*) '      boundary of the layer.'
        write(*,*) ' '
        write(*,*) ' -t = particle must be at least above the'
        write(*,*) '      bottom boundary of the layer.'
        write(*,*) ' '
        write(*,*) ' -l = by default, all layer operations use'
        write(*,*) '      variables top1 and bot1.  If this flag'
        write(*,*) '      is selected, variables top2 and bot2'
        write(*,*) '      are used instead.'
        write(*,*) ' '
        write(*,*) 'Options -F and -G only work in -start mode.'
        write(*,*) 'If one simply wants to carry through the'
        write(*,*) 'layer information without applying a layer'
        write(*,*) 'filter, use -Lall. Use -b and -t together to'
        write(*,*) 'require a particle always within the layer.'
        write(*,*) '--------------------------------------------------'
        write(*,*) 'Other optional flags are:'
        write(*,*) ' -P = probability [0,1] to keep a float (random)'
        write(*,*) ' -h = prints this message'
        write(*,*) ' -N = write the rejected floats to no.nc in'
        write(*,*) '      addition to writing the accepted floats'
        write(*,*) '      to yes.nc.  DOES NOT apply to -worms or'
        write(*,*) '      -lines modes.'
        write(*,*) ' -n = no output written to netcdf, only total'
        write(*,*) '      numbers of floats in and out.  This is'
        write(*,*) '      useful if only looking at the number of'
        write(*,*) '      traj. that satisfy a particular criteria'
        write(*,*) '      without having to wait a long time for'
        write(*,*) '      the ensembles to be written.  This flag'
        write(*,*) '      DOES NOT apply of -worms or -lines modes.' 
        write(*,*) ' -M = mask out trajectories after they have'
        write(*,*) '      passed through the sort box.  This mode'
        write(*,*) '      is only useful for -through, -fastest,'
        write(*,*) '      -lexit, and -fexit.  NOTE: The -M flag'
        write(*,*) '      does more than mask for the exit operations,'
        write(*,*) '      see above.'
        write(*,*) ' -C = copy date information once traj are stored.'
        write(*,*) ' -e = Age data is not provided in the input file,'
        write(*,*) '      so build it up from valid trajectory length.'
        write(*,*) '      Age data will not be copied to output.'
        write(*,*) ' -J = Detect all linesin/linesout, even the ones'
        write(*,*) '      that terminate within their domain (i.e.'
        write(*,*) '      the lines with indefinite terminations).'
        write(*,*) '      Lines with indefinite terminations will be'
        write(*,*) '      written to a separate output file, end.nc,'
        write(*,*) '      which be merged with the main file later.'
        write(*,*) ' -D = If in worms (lines) mode, optional file'
        write(*,*) '      listing number of exits/entrances (suseg)'
        write(*,*) '      found per traj is written. Format is traj'
        write(*,*) '      index number followed by number of subseg.'
        write(*,*) ' -c dt = centered difference to compute'
        write(*,*) '      particle velocity from positions.'
        write(*,*) ' -f dt = forward difference to compute'
        write(*,*) '      particle velocity from positions.'
        write(*,*) ' -k dt = backward difference to compute'
        write(*,*) '      particle velocity from positions.'
        write(*,*) 'NOTE: Only one of c, f, or k can be chosen. Also,'
        write(*,*) 'c, f, and k are only useful if U, V, or W are'
        write(*,*) 'selected, above.  Otherwise, they are ignored.'
        write(*,*) 'The required input for c, f, or k options is the'
        write(*,*) 'time step, in days, between trajectory position'
        write(*,*) 'updates.'
        write(*,*) ' '
        write(*,*) ' -a = ranges are computed with .lt. and .gt.'
        write(*,*) '      By default (without the flag), in or out'
        write(*,*) '      of the sorting domain is computed with'
        write(*,*) '      .le. and .ge.'
        write(*,*) ' '
        write(*,*) ' -s = stamp each traj with an ID before'
        write(*,*) '      chopping traj up into worms or lines.'
        write(*,*) '      This is useful for going back to look'
        write(*,*) '      at the original trajectories if you find'
        write(*,*) '      something interesting after later sortings.'
        write(*,*) '      ID numbers are written to the output file.'
        write(*,*) ' '
        write(*,*) ' -i = read and carry over the IDs stamped with'
        write(*,*) '      the -s option.  Without this option, ID'
        write(*,*) '      info from the input file will be ignored'
        write(*,*) '      and lost during the sorting.'
        write(*,*) ' '
        write(*,*) '--------------------------------------------------'

        return

      end subroutine print_help

!---------------------------------------
!------------------------------------------------
     logical function layer_check (dep_z,top_z,bot_z,thk_min,thk_max,l_l,l_o,l_tt,l_bb)
!------------------------------------------------     
! This function will output a T/F logical value
! depending on whether or not the current position
! of the particle is within the layer and the
! user specified requirements about the layer
! definition.
!
! dep_z = the depth of the particle at an instant
! top_z = the top of the layer at this instant
! bot_z = the bottom of the layer at this instant
! thk_min = the minimum thickness of the layer
! thk_max = the maximum thickness of the layer
! l_l   = if the user has specified the -L flag on
!         command line (min thickness)
! l_o   = if the user has specified the -O flag on
!         command line (outcropping)
! l_tt  = if the user has specified the -t flag on
!         command line (below top of layer)
! l_bb  = if the user has specified the -b flag on
!         command line (above bottom of layer)
!------------------------------------------------

       use tcdfio

       implicit none

       real, intent(in) :: dep_z
       real, intent(in) :: top_z
       real, intent(in) :: bot_z
       real, intent(in) :: thk_min
       real, intent(in) :: thk_max

       logical, intent(in) :: l_l
       logical, intent(in) :: l_o
       logical, intent(in) :: l_tt
       logical, intent(in) :: l_bb

       real :: thickness

       ! By default, set layer_check to true
       ! and this can only be set to false
       ! if one of the flags is active.
       layer_check = .true.

       !---------------------------------------
       ! Mask check: We only want to detect
       ! exits if we have valid layer data.
       ! Otherwise, no possibility of exiting.
       !---------------------------------------  
       if ( (top_z .eq. layer_top_mask) .or. &
            (bot_z .eq. layer_bot_mask) ) then
          ! Do nothing here - no chance to be out,
          ! err on the side of being in.
       else
          ! Check for the layer properties.
          if (l_l) then
             ! Check for minimum or maximum layer thickness.
             thickness = bot_z - top_z
             if ( thickness .lt. 0.0 ) then
                write(*,*) 'ERROR: Layer thickness less than 0!'
                write(*,*) 'bottom: ',bot_z
                write(*,*) 'top:    ',top_z
                stop
             endif
             if ( (thickness .lt. thk_min) .or. &
                  (thickness .gt. thk_max) ) then
                ! The layer is too thin or thick
                layer_check = .false.
             endif
          endif
          
          if (l_o) then
             ! Check for layer outcropping.
             ! (Use 10m instead of 0m as a check
             ! since a layer that close to the
             ! surface is effectively outcropped.)
             if ( top_z .gt. 10.0 ) then
                ! The layer is not outcropping,
                layer_check = .false.
             endif
          endif
          
          if (l_tt) then
             ! Check for particle position wrt.
             ! the top of layer.
             if ( dep_z .lt. top_z ) then
                ! The particle is above the top
                ! of the layer, so outside
                layer_check = .false.
             endif
          endif
          
          if (l_bb) then
             ! Check for particle position wrt.
             ! the bottom of the layer.  Also
             ! require that the bottom of the
             ! layer be greater than zero to
             ! ensure that we have a valid layer
             ! and not a masked out layer.  If
             ! we remove the bottom > 0 restriction,
             ! then we will have particles out of
             ! the layer for masked layers (-9 or 0).
             ! This is what should happen as masked
             ! layers occur when we have not found
             ! the layer, so the particle is out
             ! of the layer because there is no layer
             ! in this location.
             !if ( (dep_z .gt. bot_z) .and. (bot_z .gt. 0) ) then
             if ( dep_z .gt. bot_z ) then
                ! The particle is below the
                ! bottom of the layer, so outside
                layer_check = .false.
             endif
          endif
       endif

       return

     end function layer_check

!------------------------------------------------

!------------------------------------------------
     logical function box_check (lims,data,lens,lequ)
!------------------------------------------------     
! This function will output a T/F logical value
! depending on whether or not the current position
! and data of the particle (in data) is within the
! domain specified by the list of limits in lims.
! The variable lens is the length of lims and data.
! The flag lequ is for when equality is used or not
! used in the domain checks.
!------------------------------------------------

       implicit none

       logical, intent(in) :: lequ
       integer, intent(in) :: lens
       real, intent(in) :: lims(2,lens)
       real, intent(in) :: data(lens)

       ! By default, set box_check to true
       ! and this can only be set to false
       ! if we are outside one of the limits.
       box_check = .true.

       ! Check that we are inside the limits.
       if (lequ) then
          ! Use equality
          if ( (data(1) .ge. lims(1,1)) .and. (data(1) .le. lims(2,1)) .and. &
               (data(2) .ge. lims(1,2)) .and. (data(2) .le. lims(2,2)) .and. &
               (data(3) .ge. lims(1,3)) .and. (data(3) .le. lims(2,3)) .and. &
               (data(4) .ge. lims(1,4)) .and. (data(4) .le. lims(2,4)) .and. &
               (data(5) .ge. lims(1,5)) .and. (data(5) .le. lims(2,5)) .and. &
               (data(6) .ge. lims(1,6)) .and. (data(6) .le. lims(2,6)) .and. &
               (data(7) .ge. lims(1,7)) .and. (data(7) .le. lims(2,7)) .and. &
               (data(8) .ge. lims(1,8)) .and. (data(8) .le. lims(2,8)) .and. &
               (data(9) .ge. lims(1,9)) .and. (data(9) .le. lims(2,9)) .and. &
               (data(10) .ge. lims(1,10)) .and. (data(10) .le. lims(2,10)) .and. &
               (data(11) .ge. lims(1,11)) .and. (data(11) .le. lims(2,11)) .and. &
               (data(12) .ge. lims(1,12)) .and. (data(12) .le. lims(2,12)) .and. &
               (data(13) .ge. lims(1,13)) .and. (data(13) .le. lims(2,13)) .and. &
               (data(14) .ge. lims(1,14)) .and. (data(14) .le. lims(2,14)) .and. &
               (data(15) .ge. lims(1,15)) .and. (data(15) .le. lims(2,15)) .and. &
               (data(16) .ge. lims(1,16)) .and. (data(16) .le. lims(2,16)) .and. &
               (data(17) .ge. lims(1,17)) .and. (data(17) .le. lims(2,17)) .and. &
               (data(18) .ge. lims(1,18)) .and. (data(18) .le. lims(2,18)) ) then
             box_check = .true.
          else
             box_check = .false.
          endif
       else
          if ( (data(1) .gt. lims(1,1)) .and. (data(1) .lt. lims(2,1)) .and. &
               (data(2) .gt. lims(1,2)) .and. (data(2) .lt. lims(2,2)) .and. &
               (data(3) .gt. lims(1,3)) .and. (data(3) .lt. lims(2,3)) .and. &
               (data(4) .gt. lims(1,4)) .and. (data(4) .lt. lims(2,4)) .and. &
               (data(5) .gt. lims(1,5)) .and. (data(5) .lt. lims(2,5)) .and. &
               (data(6) .gt. lims(1,6)) .and. (data(6) .lt. lims(2,6)) .and. &
               (data(7) .gt. lims(1,7)) .and. (data(7) .lt. lims(2,7)) .and. &
               (data(8) .gt. lims(1,8)) .and. (data(8) .lt. lims(2,8)) .and. &
               (data(9) .gt. lims(1,9)) .and. (data(9) .lt. lims(2,9)) .and. &
               (data(10) .gt. lims(1,10)) .and. (data(10) .lt. lims(2,10)) .and. &
               (data(11) .gt. lims(1,11)) .and. (data(11) .lt. lims(2,11)) .and. &
               (data(12) .gt. lims(1,12)) .and. (data(12) .lt. lims(2,12)) .and. &
               (data(13) .gt. lims(1,13)) .and. (data(13) .lt. lims(2,13)) .and. &
               (data(14) .gt. lims(1,14)) .and. (data(14) .lt. lims(2,14)) .and. &
               (data(15) .gt. lims(1,15)) .and. (data(15) .lt. lims(2,15)) .and. &
               (data(16) .gt. lims(1,16)) .and. (data(16) .lt. lims(2,16)) .and. &
               (data(17) .gt. lims(1,17)) .and. (data(17) .lt. lims(2,17)) ) then
             ! In the box.
             box_check = .true.
          else
             ! Out of the box.
             box_check = .false.
          endif
       endif
       return

     end function box_check

!------------------------------------------------

