! This software is distributed under the
! terms of the GNU General Public License
! v3 or any later version.
! Copyright Stefan Gary, 2018.
!----------------------------------------

      program tcdfhist

!----------------------------------------
! This program will load the trajectory
! data (in cdf format) and grid the
! Lagrangian data onto an Eulerian grid.
!
! Version 1: Is based on tcdfprob6.f90.
!            However, the difference is
!            that geographic as well as
!            property axes can be selected
!            whereas tcdfprob6 can only
!            bin in geographic axes.  This
!            program only bins in 2 
!            dimensions, so only two axes
!            can be selected.  Some of
!            this code is generalized to
!            work with an arbitary number
!            of dimensions, but this has
!            not been completed yet.
!
!            Added -l option for using
!            top2/bot2 instead of top1/bot1.
!
! Usage: tcdfhist2d -<axis1> min max del
!                   -<axis2> min max del
!                   -I tcdf_input_file.nc
!                   [-P pts_min pts_max p_skip]
!                   [-F ] (include average and variance values)
!                   [-O ] (verbose mode)
!                   [-H ] (print help and exit)
!
! where the axes can be two of:
! -A dT/dz
! -B drho/dz
! -C age
! -D day
! -E year
! -K percentage depth within layer
! -L points along the axis (age proxy)
! -M month
! -N displacement
! -Q potential vorticity
! -R density
! -S salinity
! -T temp
! -U zonal velocity
! -V meridonal velocity
! -W vertical velocity
! -X longitude
! -Y latitude
! -Z depth
!
! WISH LIST: Allow for automatic detection
! of 1, 2, 3, or 4 sorting axes and set
! program running in 1, 2, 3, or 4 dimension
! modes.  I don't think anyone would want
! more than 4 dims because it will be hard
! to visualize.  Some possible useful combos:
! X,Y,1,1
! X,Y,Z,1
! X,Y,Z,time
! T,S,1,1
! T,S,1,time
! T,S,Z,time
! Note that time is always the last dimension.
! For compatibility with the contour finder,
! any 4D hist will be disqualified and only
! 2D (time flattened) or 3D (2D + time) will
! be allowed.  In the first case, just one
! contour will be created.  In the second case,
! many contour slices will be created.
! Also allow for 1D and 1D+time cases.
!
! Unused -G, -J
!
! added -K option.
! added C to list of avg. variables.
!   also added a scale factor to convert time
!   variable, -t option.
!----------------------------------------

        use netcdfio
        use tcdfio
        use load_tcdf
        use basicfun

        implicit none

        ! Dimensionality of the program
        ! Note output variable will always
        ! have 4 dimensions.
        ! Some dimensions can
        ! just have a length of 1, and if so,
        ! they are labeled as dummy dims.
        integer, parameter :: nax_max = 4
        integer :: nax_pt_found

        ! File, dimension, and variable IDs

        ! Input File
        integer :: ncfid

        !----Dimension, variable names-----

        ! Output CDF file
        integer :: ncoid
        
        character(len=4) :: hvnam = 'hist'

        integer :: writstart(nax_max)
        integer :: writcount(nax_max)
        integer :: ax_did(nax_max)
        integer :: ax_vid(nax_max)
        character(len=4) :: ax_vnam(nax_max)

        integer :: xavid, xvvid
        integer :: yavid, yvvid
        integer :: zavid, zvvid
        integer :: tavid, tvvid
        integer :: aavid, avvid
        integer :: cavid, cvvid
        integer :: ravid, rvvid

        ! Declare sizes of variables
        integer :: nout

        ! Counters
        ! i = axis1, j = axis2, k = axis3, l = axis4
        ! p = points, t = traj
        integer :: i, j, k, l, t, p

        ! Lagrangian variables (axis, npts)
        real, allocatable :: lag_var(:,:)
        real, allocatable :: temp_hold(:,:)
        real, allocatable :: salt_hold(:,:)
        real, allocatable :: top_hold(:,:)
        real, allocatable :: bot_hold(:,:)
        real, allocatable :: dep_hold(:,:)
        real, allocatable :: lam(:,:)
        real, allocatable :: phi(:,:)
        real, allocatable :: dep(:,:)
        real, allocatable :: temp(:,:)
        real, allocatable :: dtdz(:,:)
        real, allocatable :: rho(:,:)

        ! Axis limits
        real :: ax_min(nax_max)
        real :: ax_max(nax_max)
        real :: ax_del(nax_max)
        integer :: ax_ind(nax_max)

        ! Command line args, File names.
        integer :: num_command_arg
        integer :: arg_len, arg_count
        character(len=50) :: arg_string
        character(len=10) :: arg_int
        character(len=5) :: arg_flag
        character(len=10) :: arg_real
        character(len=7) :: histfname = 'hist.nc'

        ! Flags for keeping track of command
        ! line options
        logical :: la, lb, lc, ld, le
        logical :: lf, lg, lh, li, lk, ll
        logical :: lm, ln, lo, lp, lq
        logical :: lr, ls, lt
        logical :: lu, lv, lw
        logical :: lx, ly, lz, lll, ltt

        logical :: lxo, lyo, lzo, lto, lao, lco, lro
        logical :: l_quick
        
        ! Default values for time steps to use
        integer :: pmin = 1
        integer :: pmax = 1
        integer :: pskip = 0
        integer :: pskipcount

        ! Numbers of histoboxes
        ! (to be calculated from the
        ! domain values specified above)
        integer :: ne(nax_max)
        integer :: nc(nax_max)
        integer :: nax
        integer :: nc_max

        ! Histobox edges (axis,length)
        real, allocatable :: edges(:,:)

        ! Histobox centers (axis,length)
        real, allocatable :: centers(:,:)

        ! Histogram
        integer, allocatable :: hist(:,:,:,:)

        ! Avg and var of requested variables
        real, allocatable :: xavg(:,:,:,:), xvar(:,:,:,:)
        real, allocatable :: yavg(:,:,:,:), yvar(:,:,:,:)
        real, allocatable :: zavg(:,:,:,:), zvar(:,:,:,:)
        real, allocatable :: tavg(:,:,:,:), tvar(:,:,:,:)
        real, allocatable :: aavg(:,:,:,:), avar(:,:,:,:)
        real, allocatable :: cavg(:,:,:,:), cvar(:,:,:,:)
        real, allocatable :: ravg(:,:,:,:), rvar(:,:,:,:)

        integer, allocatable :: xcnt(:,:,:,:), ycnt(:,:,:,:)
        integer, allocatable :: zcnt(:,:,:,:), tcnt(:,:,:,:)
        integer, allocatable :: acnt(:,:,:,:), ccnt(:,:,:,:)
        integer, allocatable :: rcnt(:,:,:,:)

        real :: time_scale_factor

        !------Read command line------

        ! Initialize the flags
        la = .false.
        lb = .false.
        lc = .false.
        ld = .false.
        le = .false.
        lf = .false.
        lh = .false.
        li = .false.
        ll = .false.
        lk = .false.

        lm = .false.
        ln = .false.
        lo = .false.
        l_quick = .false.
        lp = .false.
        lq = .false.
        lr = .false.
        ls = .false.
        lt = .false.
        lu = .false.
        lv = .false.
        lw = .false.
        lx = .false.
        ly = .false.
        lz = .false.
        lll = .false.
        lxo = .false.
        lyo = .false.
        lzo = .false.
        lto = .false.
        lao = .false.
        lco = .false.
        lro = .false.
        ltt = .false.

        ! Initialize number of axis counter
        nax = 0

        ! Initialize the axis names and sizes to defaults
        ! if the number of requested axes is less than the
        ! total available axes (4).
        ! WORKING HERE
        writstart = 1 
        writcount = 1
        ax_vnam(1) = 'xxx1'
        ax_vnam(2) = 'xxx2'
        ax_vnam(3) = 'xxx3'
        ax_vnam(4) = 'xxx4'
        
        ! Axis limits
        ax_min = 0.0
        ax_max = 0.0
        ax_del = 0.0
        ! Note that for axes that are not used, this
        ! index will stay at 1 the whole time.
        ax_ind = 1

        
        !------Get command line information------
        ! First argument: Name of file to restart.
        num_command_arg = command_argument_count()
       
        if(num_command_arg .lt. 6) then
           write(*,*) ' Error: Must specify at least 1 axis and 1 input file!'
           call print_help()
           stop
        else

           arg_count = 1

           ! Loop through all the command line flags
           do

              ! Test that we are still reading all the flags.
              if ( arg_count .gt. num_command_arg ) exit

              call get_command_argument(arg_count, arg_flag, arg_len, exitcode)
              call checkexit(exitcode)
              arg_count = arg_count + 1
              
              !==============NON-AXIS COMMAND LINE FLAGS=================

              if ( index(trim(arg_flag),'-F') .ne. 0 ) then
                 ! Include average variables (listed)
                 lf = .true.

                 ! Read the requested variables from the
                 ! list following the flag
                 call get_command_argument(arg_count, arg_string, arg_len, exitcode)
                 arg_count = arg_count + 1
                 call checkexit(exitcode)

                 if ( index(trim(arg_string),'-') .ne. 0 ) then
                    ! The user has not specified a list
                    ! (this is the start of a new flag)
                    write(*,*) 'ERROR: Must include list of variables after -F!'
                    call print_help()
                    stop
                 else
                    if ( index(trim(arg_string),'X') .ne. 0 ) then
                       ! Output x variable.
                       write(*,*) 'Will output avg/var X...'
                       lxo = .true.
                    endif
                    if ( index(trim(arg_string),'Y') .ne. 0 ) then
                       ! Output y variable.
                       write(*,*) 'Will output avg/var Y...'
                       lyo = .true.
                    endif
                    if ( index(trim(arg_string),'Z') .ne. 0 ) then
                       ! Output z variable.
                       write(*,*) 'Will output avg/var Z...'
                       lzo = .true.
                    endif
                    if ( index(trim(arg_string),'T') .ne. 0 ) then
                       ! Output t variable.
                       write(*,*) 'Will output avg/var T...'
                       lto = .true.
                    endif
                    if ( index(trim(arg_string),'A') .ne. 0 ) then
                       ! Output dtdz variable.
                       write(*,*) 'Will output avg/var dtdz...'
                       lao = .true.
                    endif
                    if ( index(trim(arg_string),'C') .ne. 0 ) then
                       ! Output time variable.
                       write(*,*) 'Will output avg/var time step...'
                       lco = .true.
                    endif
                    if ( index(trim(arg_string),'R') .ne. 0 ) then
                       ! Output time variable.
                       write(*,*) 'Will output avg/var rho...'
                       lro = .true.
                    endif
                    ! ADD MORE VARIABLES HERE FOR GENERALITY!
                 endif

              elseif ( index(trim(arg_flag),'-H') .ne. 0 ) then
                 lh = .true.
                 call print_help()
                 stop

              elseif ( index(trim(arg_flag),'-O') .ne. 0 ) then
                 ! Verbose mode
                 lo = .true.

              elseif ( index(trim(arg_flag),'-q') .ne. 0 ) then
                 ! Quick time binning mode
                 l_quick = .true.

              elseif ( index(trim(arg_flag),'-l') .ne. 0 ) then
                 ! Use top2/bot2 instead of top1/bot1.
                 lll = .true.

                 layer_top_vnam = 'top2'
                 layer_bot_vnam = 'bot2'

              elseif ( index(trim(arg_flag),'-t') .ne. 0 ) then
                 ! Convert time units by the real number scale
                 ! factor before writing to output.
                 ltt = .true.

                 ! Read the scale factor from command line
                 call get_command_argument(arg_count, arg_real, arg_len, exitcode)
                 arg_count = arg_count + 1
                 call checkexit(exitcode)
                 read(arg_real,'(f10.4)') time_scale_factor

                 write(*,*)'Prop: ',trim(arg_flag),' is ',time_scale_factor

              elseif ( index(trim(arg_flag),'-P') .ne. 0 ) then
                 ! Read in time limits on domain
                 lp = .true.

                 ! Get max and min values.
                 call get_command_argument(arg_count, arg_int, arg_len, exitcode)
                 arg_count = arg_count + 1
                 call checkexit(exitcode)
                 read(arg_int,'(i10)') pmin

                 call get_command_argument(arg_count, arg_int, arg_len, exitcode)
                 arg_count = arg_count + 1
                 call checkexit(exitcode)
                 read(arg_int,'(i10)') pmax

                 call get_command_argument(arg_count, arg_int, arg_len, exitcode)
                 arg_count = arg_count + 1
                 call checkexit(exitcode)
                 read(arg_int,'(i10)') pskip

                 write(*,*)'Prop: ',trim(arg_flag),' spans ',pmin,pmax,' skip ',pskip

              elseif ( index(trim(arg_flag),'-I') .ne. 0 ) then
                 li = .true.

                 ! Get input file name
                 call get_command_argument(arg_count, arg_string, arg_len, exitcode)
                 arg_count = arg_count + 1
                 call checkexit(exitcode)

                 write(*,*) 'Open input netcdf file...'
                 ncfid = ncopn(trim(arg_string),ncnowrit,exitcode)

                 write(*,*) 'Get dimension IDs...'
                 timedid = ncdid(ncfid,timednam,exitcode)
                 trajdid = ncdid(ncfid,trajdnam,exitcode)

                 write(*,*) 'Get dimension sizes...'
                 call ncdinq(ncfid,timedid,dummy,npts,exitcode)
                 call ncdinq(ncfid,trajdid,dummy,ntraj,exitcode)
                 nall = npts*ntraj
                 nout = 0

                 write(*,*) 'There are ',npts,' points in each traj.'
                 write(*,*) 'There are ',ntraj,' trajectories.'
                 write(*,*) 'There are ',nall,' total points to sort.'
                 
              else

                 !=============AXIS DEFINING COMMAND LINE FLAGS=============

                 ! Axis detected, so augment number of axes.
                 nax = nax + 1

                 ! Check that number of axes is ok
                 if (nax .gt. nax_max) then
                    write(*,*) 'ERROR: More than ',nax_max,' axes are specified!'
                    stop
                 endif

                 if ( index(trim(arg_flag),'-A') .ne. 0 ) then
                    ! Proceed by setting flag and axis name
                    la = .true.
                    ax_vnam(nax) = dtdz_vnam

                 elseif ( index(trim(arg_flag),'-B') .ne. 0 ) then
                    lb = .true.
                    ax_vnam(nax) = drdz_vnam

                 elseif ( index(trim(arg_flag),'-C') .ne. 0 ) then
                    lc = .true.
                    ax_vnam(nax) = agevnam

                 elseif ( index(trim(arg_flag),'-D') .ne. 0 ) then
                    ld = .true.
                    ax_vnam(nax) = day_vnam

                 elseif ( index(trim(arg_flag),'-E') .ne. 0 ) then
                    le = .true.
                    ax_vnam(nax) = year_vnam

                 elseif ( index(trim(arg_flag),'-L') .ne. 0 ) then
                    ll = .true.
                    ax_vnam(nax) = 'lpts'

                 elseif ( index(trim(arg_flag),'-K') .ne. 0 ) then
                    lk = .true.
                    ax_vnam(nax) = 'perc'

                 elseif ( index(trim(arg_flag),'-M') .ne. 0 ) then
                    lm = .true.
                    ax_vnam(nax) = month_vnam

                 elseif ( index(trim(arg_flag),'-N') .ne. 0 ) then
                    ln = .true.
                    ax_vnam(nax) = 'disp'

                 elseif ( index(trim(arg_flag),'-Q') .ne. 0 ) then
                    lq = .true.
                    ax_vnam(nax) = qvnam

                 elseif ( index(trim(arg_flag),'-R') .ne. 0 ) then
                    lr = .true.
                    ax_vnam(nax) = rhovnam

                 elseif ( index(trim(arg_flag),'-S') .ne. 0 ) then
                    ls = .true.
                    ax_vnam(nax) = saltvnam

                 elseif ( index(trim(arg_flag),'-T') .ne. 0 ) then
                    lt = .true.
                    ax_vnam(nax) = tempvnam

                 elseif ( index(trim(arg_flag),'-U') .ne. 0 ) then
                    lu = .true.
                    ax_vnam(nax) = uvnam

                 elseif ( index(trim(arg_flag),'-V') .ne. 0 ) then
                    lv = .true.
                    ax_vnam(nax) = vvnam

                 elseif ( index(trim(arg_flag),'-W') .ne. 0 ) then
                    lw = .true.
                    ax_vnam(nax) = wvnam

                 elseif ( index(trim(arg_flag),'-X') .ne. 0 ) then
                    lx = .true.
                    ax_vnam(nax) = lamvnam

                 elseif ( index(trim(arg_flag),'-Y') .ne. 0 ) then
                    ly = .true.
                    ax_vnam(nax) = phivnam

                 elseif ( index(trim(arg_flag),'-Z') .ne. 0 ) then
                    lx = .true.
                    ax_vnam(nax) = depvnam
                    
                 else
                    write(*,*) 'ERROR: Unrecognized axis type! ',arg_flag
                    stop
                 endif
                    
                 ! All valid axes require min, max, and delta values.
                 call get_command_argument(arg_count, arg_real, arg_len, exitcode)
                 arg_count = arg_count + 1
                 call checkexit(exitcode)
                 read(arg_real,'(f10.4)') ax_min(nax)

                 call get_command_argument(arg_count, arg_real, arg_len, exitcode)
                 arg_count = arg_count + 1
                 call checkexit(exitcode)
                 read(arg_real,'(f10.4)') ax_max(nax)
                 
                 call check_limits(ax_min(nax),ax_max(nax))

                 call get_command_argument(arg_count, arg_real, arg_len, exitcode)
                 arg_count = arg_count + 1
                 call checkexit(exitcode)
                 read(arg_real,'(f10.4)') ax_del(nax)

                 write(*,4242) 'Prop: ',trim(arg_flag),' spans ',ax_min(nax),ax_max(nax),' step ',ax_del(nax)

              endif   ! End of axis definition check
           enddo   ! End of looping over all command line args
        endif   ! End of command line length check

4242    format(a6,a2,a7,f10.4,x,f10.4,a6,f10.4)

        ! Input error checking
        if ( li ) then
           ! We specified an input file -> good.
        else
           write(*,*) ' ERROR: Did not specifiy input file with -I' 
           call print_help()
           stop
        endif

        ! Data loading and/or set up
        if (lo) write(*,*) 'Allocating space...'
        allocate(lag_var(nax,npts))
        ! By default, the number of boxes is 1 (a collapsed axis)
        ! and the number of edges is 2 (2 edges around a single box).
        ne = 2
        nc = 1
        lag_var = 0.0
        nc_max = 0

        if ( lr .or. lro ) then
           ! Density axis specified, so
           ! need to load temp and salt
           ! to compute density.
           allocate(temp_hold(npts,1))
           allocate(salt_hold(npts,1))
           temp_hold = 0.0
           salt_hold = 0.0
        endif

        if ( lk ) then
           ! Percentage depth specifed
           ! so need to load layer top
           ! and bottom and depth.
           allocate(top_hold(npts,1))
           allocate(bot_hold(npts,1))
           allocate(dep_hold(npts,1))
           top_hold = 0.0
           bot_hold = 0.0
           dep_hold = 0.0
        endif

        if ( lp ) then
           ! If the user specified time limits, then
           ! keep them as they are.
        else
           ! The default time limits are not particularly
           ! useful, so if the user did not specify
           ! time limits, make something more sane
           ! (all time steps by default).
           write(*,*)'Auto reset time domain for all time steps!'
           pmin = 1
           pmax = npts
           pskip = 0
           write(*,*)'Prop: -P spans ',pmin,pmax,' skip ',pskip
        endif
        
        if ( lxo ) allocate(lam(npts,1))
        if ( lyo ) allocate(phi(npts,1))
        if ( lzo ) allocate(dep(npts,1))
        if ( lto ) allocate(temp(npts,1))
        if ( lao ) allocate(dtdz(npts,1))
        ! if ( lco ) ! No need to allocate anything for timestep.
        if ( lro ) allocate(rho(npts,1))

        !------Grid Initializations------

        write(*,*) '------------------------------------'
        write(*,*) ' Grid (min,max,delta,npts)'
        ! Loop over just the axes that are specified.
        ! The remaining axes will have nc = 1, ne = 2.
        do k = 1,nax
           ! Number of boxes
           nc(k) = aint((ax_max(k) - ax_min(k))/ax_del(k))

           ! Keep track of the longest axis dimension
           if ( nc(k) .gt. nc_max ) nc_max = nc(k)

           ! Number of box edges
           ne(k) = nc(k) + 1

           ! Print results to screen for error checking
           write(*,*) 'Axis ',k,': ',ax_min(k),ax_max(k),ax_del(k),nc(k)
        enddo
        write(*,*) '------------------------------------'

        write(*,*) 'Allocate space for histogram...'
        allocate(hist(nc(1),nc(2),nc(3),nc(4)))
        if (lxo) then
           allocate(xavg(nc(1),nc(2),nc(3),nc(4)))
           allocate(xvar(nc(1),nc(2),nc(3),nc(4)))
           allocate(xcnt(nc(1),nc(2),nc(3),nc(4)))
           xavg = lam_mask
           xvar = lam_mask
           xcnt = 0
        endif
        if (lyo) then
           allocate(yavg(nc(1),nc(2),nc(3),nc(4)))
           allocate(yvar(nc(1),nc(2),nc(3),nc(4)))
           allocate(ycnt(nc(1),nc(2),nc(3),nc(4)))
           yavg = phi_mask
           yvar = phi_mask
           ycnt = 0
        endif
        if (lzo) then
           allocate(zavg(nc(1),nc(2),nc(3),nc(4)))
           allocate(zvar(nc(1),nc(2),nc(3),nc(4)))
           allocate(zcnt(nc(1),nc(2),nc(3),nc(4)))
           zavg = dep_mask
           zvar = dep_mask
           zcnt = 0
        endif
        if (lto) then
           allocate(tavg(nc(1),nc(2),nc(3),nc(4)))
           allocate(tvar(nc(1),nc(2),nc(3),nc(4)))
           allocate(tcnt(nc(1),nc(2),nc(3),nc(4)))
           tavg = t_mask
           tvar = t_mask
           tcnt = 0
        endif
        if (lao) then
           allocate(aavg(nc(1),nc(2),nc(3),nc(4)))
           allocate(avar(nc(1),nc(2),nc(3),nc(4)))
           allocate(acnt(nc(1),nc(2),nc(3),nc(4)))
           aavg = dtdz_mask
           avar = dtdz_mask
           acnt = 0
        endif
        if (lco) then
           allocate(cavg(nc(1),nc(2),nc(3),nc(4)))
           allocate(cvar(nc(1),nc(2),nc(3),nc(4)))
           allocate(ccnt(nc(1),nc(2),nc(3),nc(4)))
           cavg = age_mask
           cvar = age_mask
           ccnt = 0
        endif
        if (lro) then
           allocate(ravg(nc(1),nc(2),nc(3),nc(4)))
           allocate(rvar(nc(1),nc(2),nc(3),nc(4)))
           allocate(rcnt(nc(1),nc(2),nc(3),nc(4)))
           ravg = r_mask
           rvar = r_mask
           rcnt = 0
        endif

        write(*,*) 'Initialize histograms to zero...'
        hist = 0

        ! Use the longest dimension to allocate space
        ! for the center and edge locations arrays.
        allocate(edges(nax,nc_max))
        allocate(centers(nax,nc_max))

        write(*,*) 'Compute histobox edges and centers...'
        ! Only loop over the active axes.
        do k = 1,nax

           ! Set minimum edge value
           edges(k,1) = ax_min(k)

           ! March along the axis adding all
           ! remaining box edges based on increment.
           do p = 2,ne(k)
              edges(k,p) = edges(k,p-1) + ax_del(k) 
           enddo

           ! Centers are the mean value of the
           ! surrouding edges.
           do p = 1,nc(k)
              centers(k,p) = (edges(k,p) + edges(k,p+1))/2.0
           enddo

           !write(*,*) '======================================='
           !write(*,*) 'Axis ',k,' centers:'
           !do p = 1,nc(k)
           !   write(*,*) centers(k,p)
           !enddo
        enddo
        write(*,*) '======================================='

        !------Create, define the output netcdf file-------
        
        write(*,*) 'Creating output cdf file...'
        ncoid = nccre(histfname,ncclobber,exitcode)

        write(*,*) 'Creating dimensions...'
        do k = 1,nax_max
           ax_did(k) = ncddef(ncoid,ax_vnam(k),nc(k),exitcode)
        enddo

        write(*,*) 'Creating variables...'
        do k = 1,nax_max
           ax_vid(k) = ncvdef(ncoid,ax_vnam(k),ncfloat,1,ax_did(k),exitcode)
           if ( index(trim(ax_vnam(k)),'dep') .ne. 0 ) then
              ! Indicate in netcdf file that this axis is
              ! positive down
              call ncaptc(ncoid,ax_vid(k),'positive',ncchar,4,'down',exitcode)
           endif
           if ( index(trim(ax_vnam(k)),'rho') .ne. 0 ) then
              ! Density axis is also positive down
              call ncaptc(ncoid,ax_vid(k),'positive',ncchar,4,'down',exitcode)
           endif
        enddo
        hvid = ncvdef(ncoid,hvnam,ncint,nax_max,ax_did,exitcode)
        ! The fill value complicates matters for plotting in
        ! Ferret.  Just remove it here and use ignore0 in
        ! Ferret.
        !call ncapt(ncoid,hvid,'_FillValue',ncint,1,0,exitcode)

        if (lxo) then
           xavid = ncvdef(ncoid,'xavg',ncfloat,nax_max,ax_did,exitcode)
           xvvid = ncvdef(ncoid,'xstd',ncfloat,nax_max,ax_did,exitcode)
        endif
        if (lyo) then
           yavid = ncvdef(ncoid,'yavg',ncfloat,nax_max,ax_did,exitcode)
           yvvid = ncvdef(ncoid,'ystd',ncfloat,nax_max,ax_did,exitcode)
        endif
        if (lzo) then
           zavid = ncvdef(ncoid,'zavg',ncfloat,nax_max,ax_did,exitcode)
           zvvid = ncvdef(ncoid,'zstd',ncfloat,nax_max,ax_did,exitcode)
        endif
        if (lto) then
           tavid = ncvdef(ncoid,'tavg',ncfloat,nax_max,ax_did,exitcode)
           tvvid = ncvdef(ncoid,'tstd',ncfloat,nax_max,ax_did,exitcode)
        endif
        if (lao) then
           aavid = ncvdef(ncoid,'aavg',ncfloat,nax_max,ax_did,exitcode)
           avvid = ncvdef(ncoid,'astd',ncfloat,nax_max,ax_did,exitcode)
        endif
        if (lco) then
           cavid = ncvdef(ncoid,'cavg',ncfloat,nax_max,ax_did,exitcode)
           cvvid = ncvdef(ncoid,'cstd',ncfloat,nax_max,ax_did,exitcode)
        endif
        if (lro) then
           ravid = ncvdef(ncoid,'ravg',ncfloat,nax_max,ax_did,exitcode)
           rvvid = ncvdef(ncoid,'rstd',ncfloat,nax_max,ax_did,exitcode)
        endif

        write(*,*) 'Done creating file...'
        call ncendf(ncoid,exitcode)

        write(*,*) 'Writing node locations to output file...'
        do k = 1,nax_max
           if ( ltt .and. ((index(trim(ax_vnam(k)),agevnam) .ne. 0) .or.&
                           (index(trim(ax_vnam(k)),'lpts') .ne. 0)) ) then
              write(*,*) 'Dividing the time axis by ',time_scale_factor
              call ncvpt(ncoid,ax_vid(k),1,nc(k), &
                   centers(k,1:nc(k))/time_scale_factor,exitcode)
           elseif ( index(trim(ax_vnam(k)),'xxx') .ne. 0 ) then
              ! Empty axis, 1 unit long
              call ncvpt(ncoid,ax_vid(k),1,1,1.0,exitcode)
           else
              ! Proceed as normal
              call ncvpt(ncoid,ax_vid(k),1,nc(k),centers(k,1:nc(k)),exitcode)
           endif
        enddo

        ! Writing vectors are not affected by input file.
        do k = 1,nax_max
           writstart(k) = 1
           writcount(k) = nc(k)
        enddo

!=====================Sort Points===========================

        write(*,*) 'Populating the histogram counts...'

        ! Initialize valid point counter
        nall = 0

        do t = 1,ntraj
           if (lo) write(*,*) 'For traj ',t,' of ',ntraj

           ! Specify the read/write vectors...
           lag_readst2d(1) = 1
           lag_readst2d(2) = t
           lag_readct2d(1) = npts
           lag_readct2d(2) = 1

           ! For each axis (binning dimension), get
           ! the Lagrangian data that correspond to
           ! that data.  The lag_var array is thus
           ! a series of Lagrangian time series,
           ! one time series corresponding to each
           ! of the binning dimensions.
           do k = 1,nax
              if ( index(trim(ax_vnam(k)),'xxx') .ne. 0 ) then
                 ! This axis is a dummy axis, so do nothing.
                 ! This will keep lag_var(k,:) = 0.0
              elseif ( index(trim(ax_vnam(k)),'lpts') .ne. 0 ) then
                 ! We construct age from the index value
                 ! of time along the trajectory.
                 do i = 1,npts
                    lag_var(k,i) = i
                 enddo
              elseif ( index(trim(ax_vnam(k)),'disp') .ne. 0 ) then
                 ! We need to load x, y, and z Lagrangian
                 ! variables and then compute the displacement
                 ! between each point.
                 ! NOT IMPLEMENTED YET !
              elseif ( index(trim(ax_vnam(k)),'rho') .ne. 0 ) then
                 ! The user has requested density binning.
                 ! For now, density is simply potential
                 ! density referenced to the surface.  Also,
                 ! we assume that we compute density from
                 ! available salinity and temperature.
                 call get_lag_var(ncfid,tempvnam,temp_hold(:,1))
                 call get_lag_var(ncfid,saltvnam,salt_hold(:,1))
                 do i = 1,npts
                    if ( (temp_hold(i,1) .eq. t_mask).or. &
                         (temp_hold(i,1) .gt. 50.0) .or. &
                         (temp_hold(i,1) .lt. -2.0) .or. &
                         (salt_hold(i,1) .eq. s_mask).or. &
                         (salt_hold(i,1) .gt. 50.0) .or. &
                         (salt_hold(i,1) .lt. 0.0) ) then
                       ! Flag with a fill value
                       lag_var(k,i) = r_mask
                    else
                       ! Good value!
                       ! DISABLED FOR NOW
                       !lag_var(k,i) = sigthet(0.0,temp_hold(i,1),&
                       !     salt_hold(i,1))
                       lag_var(k,i) = r_mask
                    endif
                 enddo
              elseif ( index(trim(ax_vnam(k)),'perc') .ne. 0 ) then
                 ! The user has requested binning by the relative
                 ! depth of the float within the given layer limits.
                 call get_lag_var(ncfid,layer_top_vnam,top_hold(:,1))
                 call get_lag_var(ncfid,layer_bot_vnam,bot_hold(:,1))
                 call get_lag_var(ncfid,depvnam,dep_hold(:,1))
                 do i = 1,npts
                    if ( bot_hold(i,1) .lt. top_hold(i,1) ) then
                       ! Sanity check unhappy, so print warning!
                       write(*,*) 'WARNING: Layer bottom above layer top!'
                    endif
                    if ( bot_hold(i,1) .eq. 0 ) then
                       ! No layer found at this location, so return
                       ! a mask value (otherwise, we get -INF...)
                       ! (even for the thinnest possible layer ~10m,
                       ! the particle is ~1000x above that layer,
                       ! which is beyond the depth of the ocean).
                       lag_var(k,i) = 999.9
                    else
                       lag_var(k,i) = ((bot_hold(i,1) - dep_hold(i,1))/&
                                         (bot_hold(i,1) - top_hold(i,1)))
                    endif
                 enddo
              else
                 ! All other normal Lagrangian variables
                 ! are read directly from the input file.
                 call get_lag_var(ncfid,ax_vnam(k),lag_var(k,:))
              endif
           enddo

           ! If we are to be averaging variables (in
           ! addition to simple position binning) then
           ! we also need to explicitly load the
           ! Lagrangian time series for each of the
           ! requested variables.  It is possible that
           ! the same variable may be read more than
           ! once.  If so, oh well.
           if ( lxo ) then
              call get_lag_var(ncfid,lamvnam,lam)
           endif
           if ( lyo ) then
              call get_lag_var(ncfid,phivnam,phi)
           endif
           if ( lzo ) then
              call get_lag_var(ncfid,depvnam,dep)
           endif
           if ( lto ) then
              call get_lag_var(ncfid,tempvnam,temp)
           endif
           if ( lao ) then
              call get_lag_var(ncfid,dtdz_vnam,dtdz)
           endif
           if ( lro ) then
              ! The user has requested density averaging.
              ! For now, density is simply potential
              ! density referenced to the surface.  Also,
              ! we assume that we compute density from
              ! available salinity and temperature.
              call get_lag_var(ncfid,tempvnam,temp_hold(:,1))
              call get_lag_var(ncfid,saltvnam,salt_hold(:,1))
              do i = 1,npts
                 if ( (temp_hold(i,1) .eq. t_mask).or. &
                      (temp_hold(i,1) .gt. 50.0) .or. &
                      (temp_hold(i,1) .lt. -2.0) .or. &
                      (salt_hold(i,1) .eq. s_mask).or. &
                      (salt_hold(i,1) .gt. 50.0) .or. &
                      (salt_hold(i,1) .lt. 0.0) ) then
                    ! Flag with a fill value
                    rho(i,1) = r_mask
                 else
                    ! Good value!
                    ! DISABLED FOR NOW
                    !rho(i,1) = sigthet(0.0,temp_hold(i,1),&
                    !     salt_hold(i,1))
                    rho(i,1) = r_mask
                 endif
              enddo
           endif

           ! Initialize the time skipping counter
           pskipcount = pskip

           do p = pmin,pmax
              ! For each point...

              if (pskipcount .ge. pskip) then
                 ! We have reached the end of the skip-stride.
                 ! We will use this timestep.

                 ! Reset the skipcounter for future skipping.
                 pskipcount = 0

                 ! Find the indices along each dimension
                 ! NOTE: WE COULD OPTIMIZE THIS TO RUN FASTER
                 ! BY CHECKING FOR A MASK VALUE AT EACH POINT
                 ! AND THEN SKIPPING TO THE NEXT TRAJ WHEN THE
                 ! FIRST MASK POINT IS REACHED (esp. for lines).
                 nax_pt_found = 0
                 do k = 1,nax

                    ! WORKING HERE
                    !write(*,*) k, p, lag_var(k,p)
                    if ( (index(trim(ax_vnam(k)),'lpts') .ne. 0) .and. l_quick ) then
                       ! Don't need to search for the time bin since
                       ! time is monotonically increasing.  Note, we
                       ! still need to ensure that ax_ind falls within
                       ! the binning range requested by the user.
                       if ( p .lt. edges(k,1) .or. p .gt. edges(k,ne(k)) ) then
                          ! Outside binning domain
                          ax_ind(k) = 0
                       else
                          ! Within binning domain
                          ax_ind(k) = p
                       endif
                    else
                       ! Search for a bin
                       ax_ind(k) = binsearch(edges(k,:),lag_var(k,p),1,ne(k))
                    endif
                    
                    ! Check that the point was found in the domain
                    if (ax_ind(k) .gt. 0) nax_pt_found = nax_pt_found + 1
                 enddo

                 ! Check that the point was found in the domain
                 if (nax_pt_found .eq. nax) then

                    ! The number of dimensions with valid points
                    ! equals the number of binning dimensions, so
                    ! so we have found a valid point.

                    ! Add one to the histogram
                    hist(ax_ind(1),ax_ind(2),ax_ind(3),ax_ind(4)) = &
                         hist(ax_ind(1),ax_ind(2),ax_ind(3),ax_ind(4)) + 1
                      
                    ! Add one to the total number of valid points
                    nall = nall + 1
 
                    ! Push the stack for the running mean and
                    ! variance of requested properties at this
                    ! location.
                    if ( lxo ) then
                       call olavgvar(lam(p,1),xcnt(ax_ind(1),ax_ind(2),&
                                                   ax_ind(3),ax_ind(4)),&
                                              xavg(ax_ind(1),ax_ind(2),&
                                                   ax_ind(3),ax_ind(4)),&
                                              xvar(ax_ind(1),ax_ind(2),&
                                                   ax_ind(3),ax_ind(4)))
                    endif
                    if ( lyo ) then
                       call olavgvar(phi(p,1),ycnt(ax_ind(1),ax_ind(2),&
                                                   ax_ind(3),ax_ind(4)),&
                                              yavg(ax_ind(1),ax_ind(2),&
                                                   ax_ind(3),ax_ind(4)),&
                                              yvar(ax_ind(1),ax_ind(2),&
                                                   ax_ind(3),ax_ind(4)))
                    endif
                    if ( lzo ) then
                       call olavgvar(dep(p,1),zcnt(ax_ind(1),ax_ind(2),&
                                                   ax_ind(3),ax_ind(4)),&
                                              zavg(ax_ind(1),ax_ind(2),&
                                                   ax_ind(3),ax_ind(4)),&
                                              zvar(ax_ind(1),ax_ind(2),&
                                                   ax_ind(3),ax_ind(4)))
                    endif
                    if ( lto ) then
                       call olavgvar(temp(p,1),tcnt(ax_ind(1),ax_ind(2),&
                                                    ax_ind(3),ax_ind(4)),&
                                               tavg(ax_ind(1),ax_ind(2),&
                                                    ax_ind(3),ax_ind(4)),&
                                               tvar(ax_ind(1),ax_ind(2),&
                                                    ax_ind(3),ax_ind(4)))
                    endif
                    if ( lao ) then
                       call olavgvar(dtdz(p,1),acnt(ax_ind(1),ax_ind(2),&
                                                    ax_ind(3),ax_ind(4)),&
                                               aavg(ax_ind(1),ax_ind(2),&
                                                    ax_ind(3),ax_ind(4)),&
                                               avar(ax_ind(1),ax_ind(2),&
                                                    ax_ind(3),ax_ind(4)))
                    endif
                    if ( lco ) then
                       if ( ltt ) then
                          ! Divide by the time_scale_factor
                          call olavgvar(real(p)/time_scale_factor,ccnt(ax_ind(1),ax_ind(2),&
                                                                       ax_ind(3),ax_ind(4)),&
                                                                  cavg(ax_ind(1),ax_ind(2),&
                                                                       ax_ind(3),ax_ind(4)),&
                                                                  cvar(ax_ind(1),ax_ind(2),&
                                                                       ax_ind(3),ax_ind(4)))
                       else
                          call olavgvar(real(p),ccnt(ax_ind(1),ax_ind(2),&
                                                     ax_ind(3),ax_ind(4)),&
                                                cavg(ax_ind(1),ax_ind(2),&
                                                     ax_ind(3),ax_ind(4)),&
                                                cvar(ax_ind(1),ax_ind(2),&
                                                     ax_ind(3),ax_ind(4)))
                       endif
                    endif
                    if ( lro ) then
                       call olavgvar(rho(p,1),rcnt(ax_ind(1),ax_ind(2),&
                                                   ax_ind(3),ax_ind(4)),&
                                              ravg(ax_ind(1),ax_ind(2),&
                                                   ax_ind(3),ax_ind(4)),&
                                              rvar(ax_ind(1),ax_ind(2),&
                                                   ax_ind(3),ax_ind(4)))
                    endif

                 else
                    ! This point is not in domain and
                    ! not sorted
                    nout = nout + 1

                 endif   !---End of point found check---
              else
                 ! We skip this time step and augment
                 ! the temporal skip counter.
                 pskipcount = pskipcount + 1
              endif   !---End of point skip check---
           enddo   !---End of looping over each point---
        enddo   !---End of looping over each trajectory---

        write(*,*) 'Close input file...'
        call ncclos(ncfid,exitcode)

        write(*,*) 'Deallocate trajectory information...'
        deallocate(lag_var)
        if (lxo) deallocate(lam)
        if (lyo) deallocate(phi)
        if (lzo) deallocate(dep)
        if (lto) deallocate(temp)
        if (lao) deallocate(dtdz)
        ! if (lco) ! Nothing to deallocate
        if (lro) deallocate(rho)

        ! Final computations for the standard deviation
        ! of the variables.
        do l = 1,nc(4)
           do k = 1,nc(3)
              do j = 1,nc(2)
                 do i = 1,nc(1)
                    if (lxo) then
                       ! Check that all boxes are good.
                       if ((xcnt(i,j,k,l).eq.0).or.(xavg(i,j,k,l).eq.lam_mask)) then
                          ! Bad point
                          xvar(i,j,k,l) = 0.0
                       else
                          ! Good point
                          xvar(i,j,k,l) = sqrt(xvar(i,j,k,l)/float(xcnt(i,j,k,l)))
                       endif
                    endif
                    if (lyo) then
                       if ((ycnt(i,j,k,l).eq.0).or.(yavg(i,j,k,l).eq.phi_mask)) then
                          yvar(i,j,k,l) = 0.0
                       else
                          yvar(i,j,k,l) = sqrt(yvar(i,j,k,l)/float(ycnt(i,j,k,l)))
                       endif
                    endif
                    if (lzo) then
                       if ((zcnt(i,j,k,l).eq.0).or.(zavg(i,j,k,l).eq.dep_mask)) then
                          zvar(i,j,k,l) = 0.0
                       else
                          zvar(i,j,k,l) = sqrt(zvar(i,j,k,l)/float(zcnt(i,j,k,l)))
                       endif
                    endif
                    if (lto) then
                       if ((tcnt(i,j,k,l).eq.0).or.(tavg(i,j,k,l).eq.t_mask)) then
                          tvar(i,j,k,l) = 0.0
                       else
                          tvar(i,j,k,l) = sqrt(tvar(i,j,k,l)/float(tcnt(i,j,k,l)))
                       endif
                    endif
                    if (lao) then
                       if ((acnt(i,j,k,l).eq.0).or.(aavg(i,j,k,l).eq.dtdz_mask)) then
                          avar(i,j,k,l) = 0.0
                       else
                          avar(i,j,k,l) = sqrt(avar(i,j,k,l)/float(acnt(i,j,k,l)))
                       endif
                    endif
                    if (lco) then
                       if ((ccnt(i,j,k,l).eq.0).or.(cavg(i,j,k,l).eq.age_mask)) then
                          cvar(i,j,k,l) = 0.0
                       else
                          cvar(i,j,k,l) = sqrt(cvar(i,j,k,l)/float(ccnt(i,j,k,l)))
                       endif
                    endif
                    if (lro) then
                       if ((rcnt(i,j,k,l).eq.0).or.(ravg(i,j,k,l).eq.r_mask)) then
                          rvar(i,j,k,l) = 0.0
                       else
                          rvar(i,j,k,l) = sqrt(rvar(i,j,k,l)/float(rcnt(i,j,k,l)))
                       endif
                    endif
                 enddo
              enddo
           enddo
        enddo

        write(*,*) 'Write histogram to output file...'
        call ncvpt(ncoid,hvid,writstart,writcount,hist,exitcode)
        if ( lxo ) then
           write(*,*) 'Write avg/var lam to output file...'
           call ncvpt(ncoid,xavid,writstart,writcount,xavg,exitcode)
           call ncvpt(ncoid,xvvid,writstart,writcount,xvar,exitcode)
        endif
        if ( lyo ) then
           write(*,*) 'Write avg/var phi to output file...'
           call ncvpt(ncoid,yavid,writstart,writcount,yavg,exitcode)
           call ncvpt(ncoid,yvvid,writstart,writcount,yvar,exitcode)
        endif
        if ( lzo ) then
           write(*,*) 'Write avg/var dep to output file...'
           call ncvpt(ncoid,zavid,writstart,writcount,zavg,exitcode)
           call ncvpt(ncoid,zvvid,writstart,writcount,zvar,exitcode)
        endif
        if ( lto ) then
           write(*,*) 'Write avg/var temp to output file...'
           call ncvpt(ncoid,tavid,writstart,writcount,tavg,exitcode)
           call ncvpt(ncoid,tvvid,writstart,writcount,tvar,exitcode)
        endif
        if ( lao ) then
           write(*,*) 'Write avg/var dtdz to output file...'
           call ncvpt(ncoid,aavid,writstart,writcount,aavg,exitcode)
           call ncvpt(ncoid,avvid,writstart,writcount,avar,exitcode)
        endif
        if ( lco ) then
           write(*,*) 'Write avg/var time to output file...'
           call ncvpt(ncoid,cavid,writstart,writcount,cavg,exitcode)
           call ncvpt(ncoid,cvvid,writstart,writcount,cvar,exitcode)
        endif
        if ( lro ) then
           write(*,*) 'Write avg/var rho to output file...'
           call ncvpt(ncoid,ravid,writstart,writcount,ravg,exitcode)
           call ncvpt(ncoid,rvvid,writstart,writcount,rvar,exitcode)
        endif
        write(*,*) '-----------------------------------------------------'
        write(*,*) 'Started with ntraj*npts points to sort: ',ntraj*npts
        write(*,*) 'Number of points used in sorting: ',nall
        write(*,*) 'Number of points outside the domain: ',nout
        write(*,*) '-----------------------------------------------------'

        write(*,*) 'Close hist.nc...'
        call ncclos(ncoid,exitcode)

!        write(*,*) 'Deallocate histogram space...'
!        deallocate(hist)

        write(*,*) 'End of tcdfhist.'

      end program tcdfhist

!----------------------------------------

!---------------------------------------

      subroutine print_help()

!---------------------------------------
! This subroutine prints help info
! about this program.
!---------------------------------------

        write(*,*) '-------------------------------------------'
        write(*,*) ' tcdfhist   -<axis1> min max del'
        write(*,*) '           [-<axis2> min max del]'
        write(*,*) '           [-<axis3> min max del]'
        write(*,*) '           [-<axis4> min max del]'
        write(*,*) '            -I tcdf_input_file.nc'
        write(*,*) '           [-P pts_min pts_max p_skip]'
        write(*,*) '           [-F list_of_vars] [-l]'
        write(*,*) '           [-O ] [-H ] [-t scale]'
        write(*,*) ' '
        write(*,*) '-------------------------------------------'
        write(*,*) ' This program will compute a 2D histogram'
        write(*,*) ' of the trajectory ensemble in the input'
        write(*,*) ' file in the axes specifed by the user.'
        write(*,*) ' '
        write(*,*) '-------------------------------------------'
        write(*,*) ' '
        write(*,*) ' The axes can be any of:'
        write(*,*) ' -A dT/dz   [oC/m]'
        write(*,*) ' -B drho/dz [kg/m^4]'
        write(*,*) ' -C age     [time steps]'
        write(*,*) ' -D day     [day]'
        write(*,*) ' -E year    [year]'
        write(*,*) ' -K % depth in layer [0,1]'
        write(*,*) ' -L points along the axis (age proxy)'
        write(*,*) ' -M month   [month]'
        write(*,*) ' -N displacement [m]'
        write(*,*) ' -Q potential vorticity [x10^12 m^-1 s^-1]'
        write(*,*) ' -R density [kg/m^3] pref = 0.0.  DISABLED!'
        write(*,*) ' -S salinity[PSU]'
        write(*,*) ' -T temp    [oC]'
        write(*,*) ' -U zonal velocity [m/s]'
        write(*,*) ' -V meridonal velocity [m/s]'
        write(*,*) ' -W vertical velocity [m/s]'
        write(*,*) ' -X longitude [deg]'
        write(*,*) ' -Y latitude  [deg]'
        write(*,*) ' -Z depth     [m]'
        write(*,*) ' '
        write(*,*) '-------------------------------------------'
        write(*,*) ' '
        write(*,*) ' The -l option is included to switch layer'
        write(*,*) ' operations from using top1/bot1 to'
        write(*,*) ' top2/bot2.  This allows for the inclusion' 
        write(*,*) ' of the local EDW layer in addition to the'
        write(*,*) ' thickest layer at that position.'
        write(*,*) ' '
        write(*,*) ' The -F option allows for computing the'
        write(*,*) ' bin-averaged values in each bin for the'
        write(*,*) ' given list of variables.  Also includes'
        write(*,*) ' the variance in addition to average. Ex:'
        write(*,*) ' -F XYZ'
        write(*,*) ' '
        write(*,*) ' -t scale = will divide the time step'
        write(*,*) '            increments by the real number'
        write(*,*) '            specified in scale to convert'
        write(*,*) '            from lag. time steps to actual'
        write(*,*) '            time units.'
        write(*,*) ' '
        write(*,*) ' -O will print out traj number currently'
        write(*,*) '    working on so you can see progress.'
        write(*,*) ' '
        write(*,*) ' -P will work over pts_min [1] to pts_max [all]'
        write(*,*) '    skipping every p_skip [0] points.'
        write(*,*) ' '
        write(*,*) ' -q will speed up time bin sorting only on'
        write(*,*) '    the time axis by assuming time is'
        write(*,*) '    monotonically increasing.'
        write(*,*) '-------------------------------------------'
        return

      end subroutine print_help

!---------------------------------------
