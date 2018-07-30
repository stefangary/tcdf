! This software is distributed under the terms
! of the GNU General Public License v3 or any
! later version.
! Copyright Stefan Gary, 2018.
!----------------------------------------

      program ariane2tcdf

!----------------------------------------
! This program will load ARIANE traj
! data (in cdf format) and write that
! data to tcdf format (the same as the
! FLAME trajectory output).
!
! Version 2: Reads input
! traj-by-traj and writes output
! traj-by-traj to minimize memory use
! and optimize writing to the unlimited
! dimension (ntraj) in the tcdf files.
!
! Version 3: Merged with the
! tcdfio.f90 library.  Packs variables
! as well as finds the bottom depth
! based on traj_iU|jV|kW in the input
! file.  Bottom depth is added with -H
! on the command line.
!
! Added the -T, -S flags for
! converting over temperature
! and salinity from ARIANE output.
!
! BY DEFAULT, traj_lon|lat|depth are
! always converted.
!
! Added -J flag for
! storing index values.     
!----------------------------------------

        use params
        use netcdfio
        use grids
        use load_data
        use tcdfio
        use load_tcdf
        use basicfun

        implicit none

        ! File, dimension, and variable IDs

        ! Input and Output Files
        integer :: ncfid, ncoid
        integer :: itrajdid, itimedid

        integer :: ilamvid, iphivid, idepvid, iuvid, jvvid
        integer :: itempvid, isaltvid
        integer :: olamvid, ophivid, odepvid

        !----Dimension, variable names-----

        ! Input files are in ARIANE format,
        ! output files are in FLAME format.
        ! The FLAME format is different from
        ! ARIANE format because I wanted to
        ! make the unlimited dimension equal
        ! to the number of trajectories
        ! because there are often many more
        ! trajectories than time points.

        ! ARIANE traj calc netcdf output
        character(len=5) :: itrajdnam = 'ntraj'
        character(len=9) :: itimednam = 'nb_output'
        character(len=8) :: ilamvnam = 'traj_lon'
        character(len=8) :: iphivnam = 'traj_lat'
        character(len=10) :: idepvnam = 'traj_depth'
        character(len=7) :: iuvnam = 'traj_iU'
        character(len=7) :: jvvnam = 'traj_jV'
        character(len=7) :: kwvnam = 'traj_kW'
        character(len=9) :: itempvnam = 'traj_temp'
        character(len=9) :: isaltvnam = 'traj_salt'

        ! Counters
        ! p = points(time), t = traj
        integer :: t, p, ii, jj, kk

        ! Bottom depth
        real, allocatable :: bdep_map(:,:)

        ! FLAME trajectory variables
        real, allocatable :: olam(:,:)
        real, allocatable :: ophi(:,:)
        real, allocatable :: odep(:,:)
        real, allocatable :: bdep(:,:)
        real, allocatable :: otemp(:,:)
        real, allocatable :: osalt(:,:)
        real, allocatable :: oiit(:,:)
        real, allocatable :: ojjt(:,:)
        real(kind=8) :: lam_missing(1), phi_missing(1), dep_missing(1)
        real(kind=8) :: temp_missing(1), salt_missing(1)
        real(kind=8) :: iit_missing(1),jjt_missing(1)

        ! ARIANE trajectory variables
        ! Note, I changed ARIANE to output files
        ! with floats rather than doubles.
#if defined ( ARIANE_standard )
        real(kind=8), allocatable :: ilam(:,:)
        real(kind=8), allocatable :: iphi(:,:)
        real(kind=8), allocatable :: idep(:,:)
        real(kind=8), allocatable :: iu(:,:)
        real(kind=8), allocatable :: jv(:,:)
        real(kind=8), allocatable :: kw(:,:)
#else
        real, allocatable :: ilam(:,:)
        real, allocatable :: iphi(:,:)
        real, allocatable :: idep(:,:)
        real, allocatable :: iu(:,:)
        real, allocatable :: jv(:,:)
        real, allocatable :: kw(:,:)
#endif
        real, allocatable :: itemp(:,:)
        real, allocatable :: isalt(:,:)

        ! Command line args, File names.
        integer :: num_command_arg
        integer :: arg_len, arg_count
        character(len=50) :: in_file_name
        character(len=5) :: out_file_name = 't.cdf'
        character(len=50) :: arg_string
        character(len=4) :: arg_int
        character(len=5) :: arg_flag
        character(len=10) :: arg_real

        ! Flags
        logical :: li, lh, lp, l_help
        logical :: lm, lq, ls, lt, lj
        logical :: l_verbose
        
        ! Default values for flags
        li = .false.
        lj = .false.
        lh = .false.
        lp = .false.
        lm = .false.
        lq = .false.
        ls = .false.
        lt = .false.
        l_help = .false.
        l_verbose = .false.

        write(*,*) 'Starting ariane2tcdf...'
        
        !------Get command line information------
        ! First argument: Name of file to restart.
        num_command_arg = command_argument_count()
       
        if(num_command_arg .lt. 2) then
           write(*,*) ' Too few command line arguments.'
           write(*,*) ' Must specify at least filename with'
           write(*,*) ' ariane2tcdf -I filename'
           call print_help()
           stop
        else

           arg_count = 1

           ! Loop through all the other command line flags
           do

              ! Test that we are still reading all the flags.
              if ( arg_count .gt. num_command_arg ) exit

              call get_command_argument(arg_count, arg_flag, arg_len, exitcode)
              call checkexit(exitcode)
              arg_count = arg_count + 1
              
              if ( index(trim(arg_flag),'-I') .ne. 0 ) then
                 ! Read input file name
                 li = .true.

                 call get_command_argument(arg_count, arg_string, arg_len, exitcode)
                 arg_count = arg_count + 1
                 call checkexit(exitcode)

                 ! Test validity of input file name by opening input file.
                 if(l_verbose) write(*,*) 'Open input netcdf file...'
                 ncfid = ncopn(trim(arg_string),ncnowrit,exitcode)

                 if(l_verbose) write(*,*) 'Get dimension IDs...'
                 itimedid = ncdid(ncfid,itimednam,exitcode)
                 itrajdid = ncdid(ncfid,itrajdnam,exitcode)

                 if(l_verbose) write(*,*) 'Get dimension sizes...'
                 call ncdinq(ncfid,itimedid,dummy,npts,exitcode)
                 call ncdinq(ncfid,itrajdid,dummy,ntraj,exitcode)
                 nall = npts*ntraj

                 if(l_verbose) then
                    write(*,*) 'There are ',npts,' points in each traj.'
                    write(*,*) 'There are ',ntraj,' trajectories.'
                    write(*,*) 'There are ',nall,' total points to convert.'
                 endif
                 
              elseif ( index(trim(arg_flag),'-J') .ne. 0 ) then
                 ! We will add iit,jjt - grid indeces.
                 lj = .true.
              elseif ( index(trim(arg_flag),'-H') .ne. 0 ) then
                 ! We will add bottom depth
                 lh = .true.
              elseif ( index(trim(arg_flag),'-T') .ne. 0 ) then
                 ! We will add temperature
                 lt = .true.
              elseif ( index(trim(arg_flag),'-S') .ne. 0 ) then
                 ! We will add salinity
                 ls = .true.
              elseif ( index(trim(arg_flag),'-P') .ne. 0 ) then
                 ! We will pack output
                 lp = .true.
              elseif ( index(trim(arg_flag),'-M') .ne. 0 ) then
                 ! We will check for misssing values
                 lm = .true.
              elseif ( index(trim(arg_flag),'-q') .ne. 0 ) then
                 ! We will check for misssing values quickly
                 lq = .true.
              elseif ( index(trim(arg_flag),'-h') .ne. 0 ) then
                 ! Pring help and quit
                 l_help = .true.
                 call print_help()
              elseif ( index(trim(arg_flag),'-v') .ne. 0 ) then
                 ! Verbose output
                 l_verbose = .true.
              else
                 write(*,*) 'ERROR: Unknown command line flag!'
                 stop
              endif
           enddo
        endif

        if ( (.not. lm) .and. lq ) then
           write(*,*) 'WARNING: Cannot use -q without -M.'
           write(*,*) 'Ignoring -q...'
           lq = .false.
        endif

        !------Done reading from the command line------

        ! Error check at end of command line read
        if ( .not. li ) then
           write(*,*) 'ERROR: Cannot proceed without -I filename'
           stop
        endif

        if ( l_verbose ) write(*,*) 'Get variable IDs...'
        ilamvid = ncvid(ncfid,ilamvnam,exitcode)
        iphivid = ncvid(ncfid,iphivnam,exitcode)
        idepvid = ncvid(ncfid,idepvnam,exitcode)
        if ( lh .or. lj ) then
           iuvid = ncvid(ncfid,iuvnam,exitcode)
           jvvid = ncvid(ncfid,jvvnam,exitcode)
        endif
        if ( lt ) then
           itempvid = ncvid(ncfid,itempvnam,exitcode)
        endif
        if ( ls ) then
           isaltvid = ncvid(ncfid,isaltvnam,exitcode)
        endif

        !------Create and define the tcdf output file------

        if ( l_verbose ) write(*,*) ' Creating output tcdf file...'
        ncoid = nccre(out_file_name,ncclobber,exitcode)

        ! Define dimensions.  Note that we set the traj dimension
        ! to be unlimited because this is the second (last) index
        ! in the variable.  Also, we are more likely to have more
        ! trajectories than time steps (~25,000 traj for 5,000
        ! steps - in 3 days, this is ~30 years).  We want to
        ! minimize the record offset between variables, so make
        ! the longer dimension unlimited with a 0 length.  To
        ! convert back to normal lengths, change the "0" to
        ! "number".
        if ( l_verbose ) write(*,*) ' Creating Traj and time dimensions...'
        timedid = ncddef(ncoid,timednam,npts,exitcode)
        trajdid = ncddef(ncoid,trajdnam,0,exitcode)

        vdims_out(1) = timedid
        vdims_out(2) = trajdid
         
        if ( l_verbose ) write(*,*) ' Creating Lam, phi, dep variables...'
        call cre_lag_var(ncoid,lamvnam,lp)
        call cre_lag_var(ncoid,phivnam,lp)
        call cre_lag_var(ncoid,depvnam,lp)
        if ( lj ) then
           call cre_lag_var(ncoid,iitvnam,lp)
           call cre_lag_var(ncoid,jjtvnam,lp)
        endif
        if ( lh ) then
           call cre_lag_var(ncoid,bdepvnam,lp)
        endif
        if ( lt ) then
           call cre_lag_var(ncoid,tempvnam,lp)
        endif
        if ( ls ) then
           call cre_lag_var(ncoid,saltvnam,lp)
        endif

        if ( l_verbose ) write(*,*) ' Ending output file definition stage...'
        call ncendf(ncoid,exitcode)

        if ( l_verbose ) write(*,*) ' Allocating space for one trajectory...'   
        allocate(ilam(1,npts))
        allocate(olam(npts,1))
        ilam = 0.0
        olam = 0.0   

        allocate(iphi(1,npts))
        allocate(ophi(npts,1))
        iphi = 0.0
        ophi = 0.0

        allocate(idep(1,npts))
        allocate(odep(npts,1))
        idep = 0.0
        odep = 0.0

        if ( lt ) then
           allocate(itemp(1,npts))
           allocate(otemp(npts,1))
           itemp = 0.0
           otemp = 0.0
        endif
        if ( ls ) then
           allocate(isalt(1,npts))
           allocate(osalt(npts,1))
           isalt = 0.0
           osalt = 0.0
        endif

        if ( lj ) then
           allocate(oiit(npts,1))
           allocate(ojjt(npts,1))
           oiit = 0.0
           ojjt = 0.0
        endif
           
        if ( lh .or. lj ) then
           ! Input is i,j coordinates along each particle
           allocate(iu(1,npts))
           allocate(jv(1,npts))
           iu = 0.0
           jv = 0.0
        endif

        if ( lh ) then
           ! Output
           allocate(bdep(npts,1))
           bdep = 0.0

           ! Load the model grid and build up a bottom depth map
           call load_orca_mesh_mask

           ! Build a bottom depth map, for ORCA,
           ! sum up the partial depths of the cells.
           allocate(bdep_map(imt,jmt))
           do jj = 1,jmt
              do ii = 1,imt
                 do kk = 1,kmt(ii,jj)
                    bdep_map(ii,jj) = bdep_map(ii,jj) + dzt(ii,jj,kk)
                 enddo
                 ! Check reasonable limits on bathymetry
                 if ((bdep_map(ii,jj).gt.10000.0).or.(bdep_map(ii,jj).lt.0.0)) bdep_map(ii,jj) = 0.0
              enddo
           enddo

#ifdef clip_seamounts
           write(*,*) '!!!!!!!!!!!!!!!!!!!!CAUTION!!!!!!!!!!!!!!!!!!!!'
           write(*,*) '!!!!!!!!!!!!!!!!!!!!CAUTION!!!!!!!!!!!!!!!!!!!!'
           write(*,*) '!!!!!!!!!!!!!!!!!!!!CAUTION!!!!!!!!!!!!!!!!!!!!'
           write(*,*) 'All points in vicinity of Anton Dohrn Seamount'
           ! ONLY WORKS FOR VIKING20 CUT DOMAIN!!!
           do jj = 720,740
              do ii = 920,940
                 bdep_map(ii,jj) = 2222.0
              enddo
           enddo
           write(*,*) 'and the Hebrides Seamount Terrace are set to'
           ! ONLY WORKS FOR VIKING20 CUT DOMAIN!!!
           do jj = 690,710
              do ii = 940,955
                 bdep_map(ii,jj) = 2222.0
              enddo
           enddo
           write(*,*) '2222 m depth to filter them out of the bathy'
           write(*,*) 'sorting routines but recover later if you wish.'
           write(*,*) '!!!!!!!!!!!!!!!!!!!!CAUTION!!!!!!!!!!!!!!!!!!!!'
           write(*,*) '!!!!!!!!!!!!!!!!!!!!CAUTION!!!!!!!!!!!!!!!!!!!!'
           write(*,*) '!!!!!!!!!!!!!!!!!!!!CAUTION!!!!!!!!!!!!!!!!!!!!'
#endif
        endif

        ! Get missing values
        call ncagt(ncfid,ilamvid,'missing_value',lam_missing,exitcode)
        if(l_verbose) write(*,*) 'Lon missing_value = ',lam_missing(1)
        
        call ncagt(ncfid,iphivid,'missing_value',phi_missing,exitcode)
        if(l_verbose) write(*,*) 'Lat missing_value = ',phi_missing(1)
        
        call ncagt(ncfid,idepvid,'missing_value',dep_missing,exitcode)
        if(l_verbose) write(*,*) 'Dep missing_value = ',dep_missing(1)
        
        if ( lt ) then
           call ncagt(ncfid,itempvid,'missing_value',temp_missing,exitcode)
           write(*,*) 'Temp missing_value = ',temp_missing(1)
        endif
        if ( ls ) then
           call ncagt(ncfid,isaltvid,'missing_value',salt_missing,exitcode)
           write(*,*) 'Salt missing_value = ',salt_missing(1)
        endif
        
        !========Loop over each trajectory from infile==========
        do t = 1,ntraj

           if ( l_verbose ) write(*,*) 'For traj ',t,' of ',ntraj

           !write(*,*) 'Specify the read/write vectors...'
           lag_readst2d(1) = t
           lag_readst2d(2) = 1
           lag_readct2d(1) = 1
           lag_readct2d(2) = npts

           lag_writest2d(1) = 1
           lag_writest2d(2) = t
           lag_writect2d(1) = npts
           lag_writect2d(2) = 1

           !------Read variables from infile-----
           call ncvgt(ncfid,ilamvid,lag_readst2d,lag_readct2d,ilam,exitcode)
           call ncvgt(ncfid,iphivid,lag_readst2d,lag_readct2d,iphi,exitcode)
           call ncvgt(ncfid,idepvid,lag_readst2d,lag_readct2d,idep,exitcode)

           !------If bottom depth or indeces are called for------
           if ( lh .or. lj ) then
              ! Read extra variables from infile
              call ncvgt(ncfid,iuvid,lag_readst2d,lag_readct2d,iu,exitcode)
              call ncvgt(ncfid,jvvid,lag_readst2d,lag_readct2d,jv,exitcode)
           endif

           !------If temperature is called for------
           if ( lt ) then
              ! Read extra variable from infile
              call ncvgt(ncfid,itempvid,lag_readst2d,lag_readct2d,itemp,exitcode)
           endif

           !------If salinity is called for------
           if ( ls ) then
              ! Read extra variable from infile
              call ncvgt(ncfid,isaltvid,lag_readst2d,lag_readct2d,isalt,exitcode)
           endif

           !------Transfer variables pt-by-pt------
           ! Note that here, if we decide to pack the
           ! variables, the ARIANE masking values are
           ! waaaaaay outside of the packed variables
           ! range limits.  So, need to explicitly check
           ! for range limits and assign masking to be
           ! compatible with masking.  This will be very
           ! expensive computationally because of the
           ! repeated if statements.  However, if you
           ! don't do this, then the 1e+20 (typically)
           ! missing values will be set to zero. 
           do p = 1,npts
              if ( ilam(1,p) .eq. real(lam_missing(1)) ) then
                 olam(p,1) = lam_mask
              else
                 olam(p,1) = real(ilam(1,p))
              endif

              if ( iphi(1,p) .eq. real(phi_missing(1)) ) then
                 ophi(p,1) = phi_mask
              else
                 ophi(p,1) = real(iphi(1,p))
              endif

              ! Note that for ARIANE input, we need to
              ! negate the the depth coordinate.
              if ( idep(1,p) .eq. real(dep_missing(1)) ) then
                 odep(p,1) = dep_mask
              else
                 odep(p,1) = real(-1.0*idep(1,p))
              endif
           enddo

           ! If temperature is called for:
           if ( lt ) then
              do p = 1,npts
                 if ( itemp(1,p) .eq. real(temp_missing(1)) ) then
                    otemp(p,1) = t_mask
                 else
                    otemp(p,1) = real(itemp(1,p))
                 endif
              enddo
           endif

           ! If salinity is called for:
           if ( ls ) then
              do p = 1,npts
                 if ( isalt(1,p) .eq. real(salt_missing(1)) ) then
                    osalt(p,1) = s_mask
                 else
                    osalt(p,1) = real(isalt(1,p))
                 endif
              enddo
           endif

           ! Transfer data from iu,jv to oiit,ojjt
           if ( lj ) then
              do p = 1,npts
                 ! Convert from iu,jv grid to t-grid
                 ! Any rounding is taken care of in
                 ! the writing of trajectories.
                 oiit(p,1) = iu(1,p)+0.5
                 ojjt(p,1) = jv(1,p)+0.5
              enddo
           endif
           
           ! Make bottom depth in its own loop to avoid multiple
           ! if executions.
           if ( lh ) then
              do p = 1,npts
                 ! Find the depth of the water at this point
                 ! i is on the U-grid and j is on the V-grid,
                 ! so need to shift to T-grid and finally
                 ! convert to integers.  (The shift is the
                 ! opposite to get_launch_points.sh.)
                 ! Chances are pretty good that these values
                 ! should be nearly exactly on whole T-grid
                 ! nodes for launch only so use nint to round
                 ! rather than int, just in case there are
                 ! small numerical differences from exact grid
                 ! nodes.  If the result is just slightly less
                 ! than a whole grid node, then int would move
                 ! over by a whole grid node.
                 ! In general, if you use int, then any i from
                 ! x.00001 to x.99999 would end up at x.  But
                 ! the t-grid node essentially goes from x-0.5
                 ! to x+0.5, and the rounding of nint will take
                 ! you to exactly that requirement.
                 ! (Alternatively, if you didn't convert from i
                 ! and j grids, then you could just use int because
                 ! they are a half step behind.)
                 ii = nint(iu(1,p)+0.5)
                 jj = nint(jv(1,p)+0.5)

                 ! Sanity check for particles at the boundaries,
                 ! throw them back onto the edge.
                 if ( ii .lt. 1 ) ii = 1
                 if ( ii .gt. imt ) ii = imt
                 if ( jj .lt. 1 ) jj = 1
                 if ( jj .gt. jmt ) jj = jmt

                 ! Determine bottom depth based on location
                 bdep(p,1) = bdep_map(ii,jj)
              enddo
           endif

           !------Write to output file------
           call put_lag_var(ncoid,lamvnam,olam,npts)
           call put_lag_var(ncoid,phivnam,ophi,npts)
           call put_lag_var(ncoid,depvnam,odep,npts)
           if ( lj ) then
              call put_lag_var(ncoid,iitvnam,oiit,npts)
              call put_lag_var(ncoid,jjtvnam,ojjt,npts)
           endif
           if ( lh ) then
              call put_lag_var(ncoid,bdepvnam,bdep,npts)
           endif
           if ( lt ) then
              call put_lag_var(ncoid,tempvnam,otemp,npts)
           endif
           if ( ls ) then
              call put_lag_var(ncoid,saltvnam,osalt,npts)
           endif

        enddo
        !---------------------------------------------------
        if ( l_verbose ) write(*,*) 'Cleaning up...'
        deallocate(ilam,olam)
        deallocate(iphi,ophi)
        deallocate(idep,odep)
        if ( lj ) then
           deallocate(oiit,ojjt)
        endif
        if ( lj .or. lh ) then
           deallocate(iu,jv)
        endif
        if ( lh ) then
           deallocate(bdep)
        endif
        if ( lt ) then
           deallocate(itemp,otemp)
        endif
        if ( ls ) then
           deallocate(isalt,osalt)
        endif

        if ( l_verbose ) write(*,*) ' Closing input file...'
        call ncclos(ncfid,exitcode)

        if ( l_verbose ) write(*,*) ' Closing output file...'
        call ncclos(ncoid,exitcode)

        if ( l_verbose ) write(*,*) 'End of ariane2tcdf.'

      end program ariane2tcdf

!----------------------------------------

!----------------------------------------

subroutine print_help()

  implicit none

  write(*,*) 'ariane2tcdf version 3'
  write(*,*) 'Usage:'
  write(*,*) 'ariane2tcdf -I filename'
  write(*,*) ' [-P] [-H] [-M] [-q] [-h] [-T] [-S] [-v]'
  write(*,*) ' '
  write(*,*) ' -I specifies input.  Output is t.cdf'
  write(*,*) ' -P activates packing of data into short ints'
  write(*,*) ' -H determines bottom depth from traj_iU|jV'
  write(*,*) ' -T loads temperature from traj_temp'
  write(*,*) ' -S loads salinity from traj_salt'
  write(*,*) ' -M checks for missing_values'
  write(*,*) ' -q only used with -M, checks for missing lon'
  write(*,*) '    only, not the other variables.  Quick check.'
  write(*,*) ' -J loads traj_iU|jV and saves to output file'
  write(*,*) '    NOTE: if packing is selected, -P, then the'
  write(*,*) '    the output is automatically converted to'
  write(*,*) '    integers with int() and fractional locations'
  write(*,*) '    within boxes are lost.  If packing is not'
  write(*,*) '    selected, then they are stored as reals'
  write(*,*) '    with fractional positions within boxes'
  write(*,*) '    retained.'
  write(*,*) ' -h prints this message and quits'
  write(*,*) ' -v verbose mode'
  write(*,*) ' '
  write(*,*) ' Eventually, would be nice to add'
  write(*,*) ' -R flag for density.'
#if defined ( ARIANE_standard )
  write(*,*) ' '
  write(*,*) ' NOTE: This executable was compiled with'
  write(*,*) ' the preprocessor flag ARIANE_standard'
  write(*,*) ' which allows for using ARIANE output'
  write(*,*) ' with double precision traj_lon|lat|depth|time'
  write(*,*) ' |iU|jV|kW.'
  write(*,*) ' '
  write(*,*) ' If you want to use ARIANE output where'
  write(*,*) ' c_rltype has been changed to c_rstype'
  write(*,*) ' recompile without this flag.'
  write(*,*) ' '
#else
  write(*,*) ' '
  write(*,*) ' NOTE: This executable was compiled without'
  write(*,*) ' the ARIANE_standard preprocessor flag so'
  write(*,*) ' traj_lon|lat|depth|time|iU|jV|kW are all single'
  write(*,*) ' precision rather than double precision.'
  write(*,*) ' This is achieved by changing c_rltype to'
  write(*,*) ' c_rstype in ARIANE mod_save_netcdf.f90'
  write(*,*) ' for just those four variables.  Temperature'
  write(*,*) ' and salinity are always output at single'
  write(*,*) ' precision.'
#endif

end subroutine print_help

!----------------------------------------

