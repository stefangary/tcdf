!=======================================
! This software is distributed under the
! terms of the GNU GPL 3.0 and any later
! version.
! Copyright, Stefan Gary 2018
!=======================================
  program tcdf_flatten_hist
!=======================================

    use basicfun
    use netcdfio

    implicit none

    ! Command line args, File names.
    integer :: num_command_arg
    integer :: arg_len, arg_count
    character(len=50) :: arg_string
    character(len=4) :: arg_int
    character(len=5) :: arg_flag
    character(len=10) :: arg_real

    ! Flags for presence of variables on command line.
    logical :: lf, li, l_verbose

    ! Flatten dimension
    integer :: flat_dim = 1

    ! File ID
    integer :: hist_fid
    integer :: flat_hist_fid

    integer :: hist_vid
    integer :: i_vid
    integer :: j_vid
    integer :: d_vid
    integer :: dim_id
    
    ! Histogram numbers of dimensions
    integer, parameter :: nax_max = 4

    ! Histogram number of bin centers for each dimension
    integer :: nc(nax_max)
    integer :: nc_flat(nax_max)
    integer :: nc_max
    character(len=10) :: dnam(nax_max)
    integer :: dim_var_id(nax_max)
    
    ! Other NetCDF temporary variables
    integer :: vtype
    integer :: nvdims
    integer :: nvatts
    integer :: dimid
    integer :: dimsiz
    integer :: max_iterations
    integer :: mindex(2)

    ! Data holders
    integer, allocatable :: hist(:,:,:,:)
    integer, allocatable :: flat_hist(:,:,:,:)
    integer, allocatable :: tmp_flat_hist(:,:,:)
    
    ! dim_var_val contains the dimension variables.
    ! Accessed by dim,value_along_dimension
    real, allocatable :: dim_var_val(:,:)
    integer :: ii, jj, kk, ll, sum_hist

    !-----------------------------------------
    write(*,*) ' Starting tcdf_flatten_hist...'

    ! Initialize flags
    li = .false.
    lf = .false.
    l_verbose = .false.
    nc = 0
    nc_max = 0
    
    !------Get command line information------
    ! Arguments can vary in order
    num_command_arg = command_argument_count()
       
    if(num_command_arg .lt. 2) then
       write(*,*) ' Error: Too few command line arguments.'
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

          ! Probability flag
          if( index(trim(arg_flag),'-F') .ne. 0 ) then
             lf = .true.

             ! Read in the probability number:
             call get_command_argument(arg_count, arg_int, arg_len, exitcode)
             call checkexit(exitcode)
             arg_count = arg_count + 1
             read(arg_int,'(i4)') flat_dim
             write(*,*) 'Will flatten along dim ',flat_dim
              
          elseif ( index(trim(arg_flag),'-v') .ne. 0 ) then
             ! Verbose mode
             l_verbose = .true.

          elseif ( index(trim(arg_flag),'-h') .ne. 0 ) then
             ! print help and exit
             call print_help()
             stop

          elseif ( index(trim(arg_flag),'-I') .ne. 0 ) then
             ! Read input file name
             li = .true.

             call get_command_argument(arg_count, arg_string, arg_len, exitcode)
             arg_count = arg_count + 1
             call checkexit(exitcode)

             ! Test validity of input file name by opening input file.
             if ( l_verbose ) write(*,*) 'Open input netcdf file...'
             hist_fid = ncopn(trim(arg_string),ncnowrit,exitcode)

          else
             write(*,*) 'ERROR: Unrecognized command line flag.'
             write(*,*) arg_flag
             call print_help()
             stop   
          endif
       enddo   ! Done looping over command line args.
    endif   ! Done reading command line options.
    !=======================================================

    !=======================================================
    ! ERROR CHECKING
    !=======================================================

    if ( .not. li ) then
       write(*,*) 'ERROR: no input file specified.'
       stop
    endif

    if ( l_verbose ) then
       write(*,*) 'Will flatten histogram along dim = ',flat_dim
    endif

    if ( flat_dim .gt. nax_max ) then
       write(*,*) 'ERROR: Given dimension number ',flat_dim,' exceeds max allowed dimensions.'
       stop
    endif
    
    !=======================================================
    ! LOAD HISTOGRAM
    !=======================================================

    ! Get variable ID
    if ( l_verbose ) write(*,*) 'Getting variable ID...'
    hist_vid = ncvid(hist_fid,'hist',exitcode)

    ! Get variable dimensions
    if ( l_verbose ) write(*,*) 'Getting variable dimensions...'
    call ncvinq(hist_fid,hist_vid,arg_string,vtype,nvdims,vdims4d,nvatts,exitcode)

    ! Get dimension sizes
    if ( l_verbose ) write(*,*) 'Getting dimension sizes...'
    do ii = 1,nax_max
       ! For each axis
       call ncdinq(hist_fid,vdims4d(ii),dnam(ii),nc(ii),exitcode)
       if ( l_verbose ) write(*,*) 'Found dimension ',trim(dnam(ii)),&
            ' with length ',nc(ii),'.'

       if ( nc(ii) .gt. nc_max ) nc_max = nc(ii)
       
       ! Set reading limits to get histogram, later
       readst4d(ii) = 1
       readct4d(ii) = nc(ii)

       ! If the dimension is the one to be flattened,
       ! assign it to the number of centers for the
       ! output, flattened histogram.
       if ( ii .eq. flat_dim ) then
          nc_flat(ii) = 1
       else
          nc_flat(ii) = nc(ii)
       endif

       ! Set writing limits for flattened histogram, used later
       writest4d(ii) = 1
       writect4d(ii) = nc_flat(ii)
       
    enddo

    ! Check for sane request.
    if ( nc(flat_dim) .eq. 1 ) then
       write(*,*) 'ERROR: dimension ',flat_dim,' already length 1.'
       write(*,*) 'Quit.'
       stop
    endif
    
    ! Allocate space for hist, flat_hist, tmp_flat_hist
    if ( l_verbose ) write(*,*) 'Allocating work space...'
    allocate(hist(nc(1),nc(2),nc(3),nc(4)))
    allocate(flat_hist(nc_flat(1),nc_flat(2),nc_flat(3),nc_flat(4)))
    if( flat_dim .eq. 1 ) then
       allocate(tmp_flat_hist(nc_flat(2),nc_flat(3),nc_flat(4)))
    elseif ( flat_dim .eq. 2 ) then
       allocate(tmp_flat_hist(nc_flat(1),nc_flat(3),nc_flat(4)))
    elseif ( flat_dim .eq. 3 ) then
       allocate(tmp_flat_hist(nc_flat(1),nc_flat(2),nc_flat(4)))
    elseif ( flat_dim .eq. 4 ) then
       allocate(tmp_flat_hist(nc_flat(1),nc_flat(2),nc_flat(3)))
    else
       write(*,*) 'ERROR: Cannot work with flat_dim > 4.'
       stop
    endif
    allocate(dim_var_val(nax_max,nc_max))
    dim_var_val = 0.0
    
    ! Get histogram
    if ( l_verbose ) write(*,*) 'Loading histogram...'
    call ncvgt(hist_fid,hist_vid,readst4d,readct4d,hist,exitcode)

    !=======================================================
    ! INITIALIZATIONS
    !=======================================================

    ! Always need dimension information
    if ( l_verbose ) write(*,*) 'Read dimensions...'
    do ii = 1,nax_max
       d_vid = ncvid(hist_fid,trim(dnam(ii)),exitcode)
       call ncvgt(hist_fid,d_vid,1,nc(ii),dim_var_val(ii,1:nc(ii)),exitcode)
    enddo
    
    ! Create output file
    if ( l_verbose ) write(*,*) 'Create output file...'
    flat_hist_fid = nccre('flat_hist.nc',ncclobber,exitcode)

    do ii = 1,nax_max
       ! Create dimensions
       if ( l_verbose ) write(*,*) 'Create output dimensions...'
       vdims4d(ii) = ncddef(flat_hist_fid,dnam(ii),nc_flat(ii),exitcode)

       ! Create dimension variables
       if ( l_verbose ) write(*,*) 'Create output dimension variables...'
       dim_var_id(ii) = ncvdef(flat_hist_fid,dnam(ii),ncfloat,1,vdims4d(ii),exitcode)
    enddo
       
    ! Create hist variable in flat_hist.nc
    if ( l_verbose ) write(*,*) 'Create output flat_hist variable...'
    hist_vid = ncvdef(flat_hist_fid,'hist',ncint,nax_max,vdims4d,exitcode)
       
    ! Exit define mode
    if ( l_verbose ) write(*,*) 'Exit define mode...'
    call ncendf(flat_hist_fid,exitcode)
       
    ! Copy dimension variables over
    do ii = 1,nax_max
       writest1d(1) = 1
       writect1d(1) = nc_flat(ii)
       if ( ii .eq. flat_dim ) then
          ! Put the average value of the dimension variable
          ! for the flat_dim
          call ncvpt(flat_hist_fid,dim_var_id(ii),writest1d,writect1d,&
               sum(dim_var_val(ii,1:nc(flat_dim)))/nc(flat_dim),exitcode)
       else
          ! Simple copy of dimension variables.
          call ncvpt(flat_hist_fid,dim_var_id(ii),writest1d,writect1d,&
               dim_var_val(ii,writest1d(1):writect1d(1)),exitcode)
       endif
    enddo

    !=======================================================
    ! FLATTEN THE DIMENSION
    !=======================================================

    if ( l_verbose ) write(*,*) 'Flattening...'
    
    ! Sum over the dimension
    tmp_flat_hist = sum(hist,flat_dim)

    if ( l_verbose ) write(*,*) 'Copying flattened variable to output...'
    
    ! Copy back to output variable with the dummy dimension.
    if ( flat_dim .eq. 1 ) then
       do jj = 1,nc_flat(2)
          do kk = 1,nc_flat(3)
             do ll = 1,nc_flat(4)
                flat_hist(1,jj,kk,ll) = tmp_flat_hist(jj,kk,ll)
             enddo
          enddo
       enddo
    elseif ( flat_dim .eq. 2 ) then
       do ii = 1,nc_flat(1)
          do kk = 1,nc_flat(3)
             do ll = 1,nc_flat(4)
                flat_hist(ii,1,kk,ll) = tmp_flat_hist(ii,kk,ll)
             enddo
          enddo
       enddo
    elseif ( flat_dim .eq. 3 ) then
       do ii = 1,nc_flat(1)
          do jj = 1,nc_flat(2)
             do ll = 1,nc_flat(4)
                flat_hist(ii,jj,1,ll) = tmp_flat_hist(ii,jj,ll)
             enddo
          enddo
       enddo
    elseif ( flat_dim .eq. 4 ) then
       do ii = 1,nc_flat(1)
          do jj = 1,nc_flat(2)
             do kk = 1,nc_flat(3)
                flat_hist(ii,jj,kk,1) = tmp_flat_hist(ii,jj,kk)
             enddo
          enddo
       enddo
    else
       write(*,*) 'ERROR: Cannot work with flat_dim > 4.'
       stop
    endif

    !=======================================================
    ! WRITE OUTPUT
    !=======================================================

    if ( l_verbose ) write(*,*) 'Writing output...'
    
    call ncvpt(flat_hist_fid,hist_vid,writest4d,writect4d,flat_hist,exitcode)
    call ncclos(flat_hist_fid,exitcode)

    !=======================================================
    ! CLEAN UP
    !=======================================================
    if ( l_verbose ) write(*,*) 'Closing input file...'
    call ncclos(hist_fid,exitcode)

    deallocate(dim_var_val)
    deallocate(hist)
    deallocate(flat_hist)
    deallocate(tmp_flat_hist)
    write(*,*) 'tcdf_flatten_hist done.'

!=======================================
  end program tcdf_flatten_hist
!=======================================

  subroutine print_help()
    write(*,*) ' '
    write(*,*) 'tcdf_flatten_hist -F <dimension_number> -I <file.nc>'
    write(*,*) '                 [-v] [-h]'
    write(*,*) ' '
    write(*,*) 'This module will flatten a histogram'
    write(*,*) 'by adding up all the counts along the'
    write(*,*) 'specified dimension.'
    write(*,*) ' '
    write(*,*) 'Why would one want to do this when it'
    write(*,*) 'is possible to specify flattened dims'
    write(*,*) 'in tcdfhist? To avoid having to redo'
    write(*,*) 'the binning process once it has been'
    write(*,*) 'done across many dimensions.'
    write(*,*) ' '
    write(*,*) 'Currently, the histogram is called'
    write(*,*) 'flat_hist in a netCDF file, more flexibility'
    write(*,*) 'can be incorporated later.'
    write(*,*) ' '
    write(*,*) ' -F = specify the dimension number to flatten.'
    write(*,*) '      e.g. -F 3 will add up all the values along'
    write(*,*) '      the third dimension of the histogram.'
    write(*,*) '      DEFAULT: 1'
    write(*,*) ' '
    write(*,*) ' -I = specify the name of the input file.'
    write(*,*) ' '
    write(*,*) ' -v = verbose mode.'
    write(*,*) ' '
    write(*,*) ' -h = print this message and quit.'
    write(*,*) ' '
  end subroutine print_help
