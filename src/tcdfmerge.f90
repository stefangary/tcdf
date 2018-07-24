!----------------------------------------
! This software is distributed under the
! terms of the GNU GPL 3.0 and any later
! version.
! Copyright Stefan Gary, 2018
!----------------------------------------

      program tcdfmerge

!----------------------------------------
! This program will load two trajectory
! ensembles (in cdf format) and append
! the first data set to the second data
! set.
!
! Which variables to append are specified
! by the command line flags:
!
! -U = zonal velocity
! -V = meridional velocity
! -W = vertical velocity
! -T = temperature
! -S = salinity
! -R = density
! -Q = potential vorticity
! -A = vertical gradient of temp
! -B = vertical gradient of dens
! -C = date information
! -E = age of float
!
! The -h flag prints help
!
! (x,y,z) position variables are
! automatically copied and do not
! need separate flags.
!----------------------------------------

        use tcdfio
        use load_tcdf
        use netcdfio
        use basicfun

        implicit none

        ! File, dimension, and variable IDs
        ! _i = input file
        ! _a = input file to append to.
        integer :: ncid_i, ncid_a, vid

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
        real, allocatable :: topl(:)
        real, allocatable :: botl(:)
        real, allocatable :: bdep(:)
        real, allocatable :: iit(:)
        real, allocatable :: jjt(:)

        ! Counters
        ! p = points(time), t = traj
        integer :: t, p

        ! Command line args, File names.
        integer :: num_command_arg
        integer :: arg_len, arg_count
        character(len=100) :: arg_string
        character(len=4) :: arg_int
        character(len=5) :: arg_flag
        character(len=10) :: arg_real

        ! Flags for presence of variables on command line.
        !logical :: lx, ly, lz
        logical :: lu, lv, lw
        logical :: lt, ls, lr
        logical :: lq, ll, lh, lj
        logical :: la, lb, lc, le

        ! Local numbers of input and append trajectories
        ! and points.
        integer :: ntraj_i, ntraj_a
        integer :: npts_i, npts_a

        !-----------------------------------------
        write(*,*) ' Starting tcdfmerge...'

        ! Initialize domain flags
        !lx = .false.
        !ly = .false.
        !lz = .false.
        lu = .false.
        lv = .false.
        lw = .false.
        lt = .false.
        ls = .false.
        lr = .false.
        lq = .false.
        la = .false.
        lb = .false.
        lc = .false.
        le = .false.
        ll = .false.
        lh = .false.
        lj = .false.
        
        !------Get command line information------
        ! First argument: input file name
        ! Second argument: append file name
        ! Other arguments are single flags and
        ! can vary in any order
        num_command_arg = command_argument_count()
       
        if(num_command_arg .lt. 2) then
           write(*,*) ' Error: Too few command line arguments.'
           call print_help()
           stop
        else

           ! Get first file name
           arg_count = 1
           call get_command_argument(arg_count, arg_string, arg_len, exitcode)
           call checkexit(exitcode)
           arg_count = arg_count + 1

           ! Check that name is valid by
           ! opening file.
           ncid_i = ncopn(trim(arg_string),ncnowrit,exitcode)

           ! Get append file name
           call get_command_argument(arg_count, arg_string, arg_len, exitcode)
           call checkexit(exitcode)
           arg_count = arg_count + 1

           ! Check that name is valid by
           ! opening file.
           ncid_a = ncopn(trim(arg_string),ncwrite,exitcode)

           ! Loop through all the other command line flags
           do

              ! Test that we are still reading all the flags.
              if ( arg_count .gt. num_command_arg ) exit

              call get_command_argument(arg_count, arg_flag, arg_len, exitcode)
              call checkexit(exitcode)
              arg_count = arg_count + 1
              
              if ( index(trim(arg_flag),'-U') .ne. 0 ) then
                 lu = .true.

              elseif ( index(trim(arg_flag),'-V') .ne. 0 ) then
                 lv = .true.

              elseif ( index(trim(arg_flag),'-W') .ne. 0 ) then
                 lw = .true.

              elseif ( index(trim(arg_flag),'-T') .ne. 0 ) then
                 lt = .true.

              elseif ( index(trim(arg_flag),'-S') .ne. 0 ) then
                 ls = .true.
                 
              elseif ( index(trim(arg_flag),'-R') .ne. 0 ) then
                 lr = .true.
                 
              elseif ( index(trim(arg_flag),'-Q') .ne. 0 ) then
                 lq = .true.

              elseif ( index(trim(arg_flag),'-A') .ne. 0 ) then
                 la = .true.

              elseif ( index(trim(arg_flag),'-B') .ne. 0 ) then
                 lb = .true.
                 
              elseif ( index(trim(arg_flag),'-E') .ne. 0 ) then
                 le = .true.

              elseif ( index(trim(arg_flag),'-C') .ne. 0 ) then
                 lc = .true.

              elseif ( index(trim(arg_flag),'-L') .ne. 0 ) then
                 ll = .true.

              elseif ( index(trim(arg_flag),'-H') .ne. 0 ) then
                 lh = .true.

              elseif ( index(trim(arg_flag),'-J') .ne. 0 ) then
                 lj = .true.

              elseif ( index(trim(arg_flag),'-l') .ne. 0 ) then
                 ll = .true.
                 layer_top_vnam = 'top2'
                 layer_bot_vnam = 'bot2'
                 
              elseif ( index(trim(arg_flag),'-h') .ne. 0 ) then
                 ! Print help information and exit
                 call print_help()
                 stop

              else
                 write(*,*) 'Error: Unexpected flag.'
                 stop
              endif

           enddo

        endif   ! Done reading command line options.
        !=======================================================

        !=======================================================
        !                  INITIALIZATIONS
        !=======================================================

        ! Gather information about the input files.

        write(*,*) 'Get dimension IDs...'
        timedid = ncdid(ncid_i,timednam,exitcode)
        trajdid = ncdid(ncid_i,trajdnam,exitcode)

        write(*,*) 'Get dimension sizes...'
        call ncdinq(ncid_i,timedid,dummy,npts_i,exitcode)
        call ncdinq(ncid_i,trajdid,dummy,ntraj_i,exitcode)

        write(*,*) '-----------File 1----------------------'
        write(*,*) 'There are ',npts_i,' points in each traj.'
        write(*,*) 'There are ',ntraj_i,' trajectories.'

        write(*,*) 'Get dimension IDs...'
        timedid = ncdid(ncid_a,timednam,exitcode)
        trajdid = ncdid(ncid_a,trajdnam,exitcode)

        write(*,*) 'Get dimension sizes...'
        call ncdinq(ncid_a,timedid,dummy,npts_a,exitcode)
        call ncdinq(ncid_a,trajdid,dummy,ntraj_a,exitcode)

        write(*,*) '-----------File 2----------------------'
        write(*,*) 'There are ',npts_a,' points in each traj.'
        write(*,*) 'There are ',ntraj_a,' trajectories.'

        if ( npts_i .ne. npts_a ) then
           write(*,*) 'ERROR: Trajectories have different lengths!'
           call print_help()
           stop
        endif

        npts = npts_a

        write(*,*) 'Allocating space for one trajectory...'
        allocate(lam(npts))
        allocate(phi(npts))
        allocate(dep(npts))
        if (lu) allocate(u(npts))
        if (lv) allocate(v(npts))
        if (lw) allocate(w(npts))
        if (lt) allocate(temp(npts))
        if (ls) allocate(salt(npts))
        if (lr) allocate(rho(npts))
        if (lq) allocate(q(npts))
        if (la) allocate(dtdz(npts))
        if (lb) allocate(drdz(npts))
        if (le) allocate(age(npts))
        if (lc) then
           allocate(year(npts))
           allocate(month(npts))
           allocate(day(npts))
        endif
        if (ll) then
           allocate(topl(npts))
           allocate(botl(npts))
        endif
        if (lh) allocate(bdep(npts))
        if (lj) then
           allocate(iit(npts))
           allocate(jjt(npts))
        endif
        
        ! Specify read vectors
        lag_readst2d(1) = 1
        lag_readct2d(1) = npts
        lag_readct2d(2) = 1

        lag_writest2d(1) = 1
        lag_writest2d(2) = ntraj_a

        lag_writect2d(1) = npts
        lag_writect2d(2) = 1

        write(*,*) 'Appending: Loop over each trajectory in file 1...'
        do t = 1,ntraj_i
           
           !write(*,*) 'For traj ',t,' out of ',ntraj_i
        
           !write(*,*) 'Specify the read vectors...'
           lag_readst2d(2) = t

           !write(*,*) 'Load traj data...'
           call get_lag_var(ncid_i,lamvnam,lam)
           call get_lag_var(ncid_i,phivnam,phi)
           call get_lag_var(ncid_i,depvnam,dep)
           if (lu) call get_lag_var(ncid_i,uvnam,u)
           if (lv) call get_lag_var(ncid_i,vvnam,v)
           if (lw) call get_lag_var(ncid_i,wvnam,w)
           if (lt) call get_lag_var(ncid_i,tempvnam,temp)
           if (ls) call get_lag_var(ncid_i,saltvnam,salt)
           if (lr) call get_lag_var(ncid_i,rhovnam,rho)
           if (lq) call get_lag_var(ncid_i,qvnam,q)
           if (la) call get_lag_var(ncid_i,dtdz_vnam,dtdz)
           if (lb) call get_lag_var(ncid_i,drdz_vnam,drdz)
           if (le) call get_lag_var(ncid_i,agevnam,age)
           if (lc) then
              call get_lag_var(ncid_i,year_vnam,year)
              call get_lag_var(ncid_i,month_vnam,month)
              call get_lag_var(ncid_i,day_vnam,day)
           endif
           if (ll) then
              call get_lag_var(ncid_i,layer_top_vnam,topl)
              call get_lag_var(ncid_i,layer_bot_vnam,botl)
           endif
           if (lh) call get_lag_var(ncid_i,bdepvnam,bdep)
           if (lj) then
              call get_lag_var(ncid_i,iitvnam,iit)
              call get_lag_var(ncid_i,jjtvnam,jjt)
           endif

           !write(*,*) ' Done loading data for traj ',t,' out of ',ntraj

           ! Write this traj to output
           lag_writest2d(2) = lag_writest2d(2) + 1

           call put_lag_var(ncid_a,lamvnam,lam,npts)
           call put_lag_var(ncid_a,phivnam,phi,npts)
           call put_lag_var(ncid_a,depvnam,dep,npts)
           if (lu) call put_lag_var(ncid_a,uvnam,u,npts)
           if (lv) call put_lag_var(ncid_a,vvnam,v,npts)
           if (lw) call put_lag_var(ncid_a,wvnam,w,npts)
           if (lt) call put_lag_var(ncid_a,tempvnam,temp,npts)
           if (ls) call put_lag_var(ncid_a,saltvnam,salt,npts)
           if (lr) call put_lag_var(ncid_a,rhovnam,rho,npts)
           if (lq) call put_lag_var(ncid_a,qvnam,q,npts)
           if (la) call put_lag_var(ncid_a,dtdz_vnam,dtdz,npts)
           if (lb) call put_lag_var(ncid_a,drdz_vnam,drdz,npts)
           if (le) call put_lag_var(ncid_a,agevnam,age,npts)
           if (lc) then
              call put_lag_var(ncid_a,year_vnam,year,npts)
              call put_lag_var(ncid_a,month_vnam,month,npts)
              call put_lag_var(ncid_a,day_vnam,day,npts)
           endif
           if (ll) then
              call put_lag_var(ncid_a,layer_top_vnam,topl,npts)
              call put_lag_var(ncid_a,layer_bot_vnam,botl,npts)
           endif
           if (lh) call put_lag_var(ncid_a,bdepvnam,bdep,npts)
           if (lj) then
              call put_lag_var(ncid_a,iitvnam,iit,npts)
              call put_lag_var(ncid_a,jjtvnam,jjt,npts)
           endif
           
           ! Check for overrun
           if ( lag_writest2d(2) .gt. (ntraj_a + ntraj_i) ) then
              write(*,*) 'ERROR: Writing too many traj!'
              write(*,*) 'lag_writest2d(2) = ',lag_writest2d(2)
              write(*,*) 'ntraj_a = ',ntraj_a
              write(*,*) 'ntraj_i = ',ntraj_i
           endif

        enddo

        !==================================================
        !                    CLEAN UP
        !==================================================
        write(*,*) 'Close input file...'
        call ncclos(ncid_i,exitcode)
        
        write(*,*) 'Close append file...'
        call ncclos(ncid_a,exitcode)

        write(*,*) 'Deallocate input trajectories...'
        deallocate(lam,phi,dep)
        if (lt) deallocate(temp)
        if (ls) deallocate(salt)
        if (lr) deallocate(rho)
        if (lu) deallocate(u)
        if (lv) deallocate(v)
        if (lw) deallocate(w)
        if (lq) deallocate(q)
        if (lb) deallocate(drdz)
        if (la) deallocate(dtdz)
        if (le) deallocate(age)
        if (lc) deallocate(year,month,day)
        if (ll) deallocate(topl,botl)
        if (lh) deallocate(bdep)
        if (lj) deallocate(iit,jjt)
        write(*,*) 'End of tcdfmerge.'

      end program tcdfmerge

!---------------------------------------

!---------------------------------------

      subroutine print_help()

!---------------------------------------
! This subroutine prints help info
! about this program.
!---------------------------------------

        write(*,*) '-------------------------------------------'
        write(*,*) ' tcdfmerge infile1.nc infile2.nc'
        write(*,*) '          [-<U|V|W|T|S|R|Q|E|A|B|C|L|l|J>] [-h]'
        write(*,*) '-------------------------------------------'
        write(*,*) 'This program will append the traj in the'
        write(*,*) 'first file to the second file.  The'
        write(*,*) 'copied are set by the flags and (x,y,z)'
        write(*,*) 'are always copied.'
        write(*,*) '-------------------------------------------'
        write(*,*) 'Variable flags:'
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
        write(*,*) ' -C = year, month, and day of time step.'
        write(*,*) ' -L = layer top and bottom (top1, bot1)'
        write(*,*) ' -l = use top2 and bot2 instead of *1.'
        write(*,*) ' -H = include bottom depth.'
        write(*,*) ' -J = include i,j index node locations'
        write(*,*) '      on the orginal model grid (iit,jjt).'
        write(*,*) ' '
        write(*,*) '-------------------------------------------'
        write(*,*) 'Other optional flags are:'
        write(*,*) ' -h = prints this message'
        write(*,*) '-------------------------------------------'

        return

      end subroutine print_help

!---------------------------------------
