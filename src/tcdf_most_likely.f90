!=======================================
! This software is distributed under the
! terms of the GNU GPL 3.0 and any later
! version.
! Copyright, Stefan Gary 2018
!=======================================
  program tcdf_most_likely
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
    logical :: lp, li, l_verbose, lx

    ! Probability contour level
    real :: prob_contour = 0.95
    integer :: iteration_limit = 100

    ! File ID
    integer :: hist_fid
    integer :: mask_fid
    integer :: xyz_fid

    integer :: hist_vid
    integer :: mask_vid
    integer :: i_vid
    integer :: j_vid
    integer :: d_vid
    integer :: dim_id
    
    ! Histogram numbers of dimensions
    integer, parameter :: nax_max = 4

    ! Histogram number of bin centers for each dimension
    integer :: nc(nax_max)
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
    integer, allocatable :: mask(:,:,:,:)
    logical, allocatable :: hist_mask(:,:)
    ! dim_var_val contains the dimension variables.
    ! Accessed by dim,value_along_dimension
    real, allocatable :: dim_var_val(:,:)
    integer :: ii, jj, ll, sum_hist
    integer :: ij(2)

    ! Accumulator for probability.
    ! Keep at double precision to allow
    ! for not loosing very small contributions
    ! to the total sum.
    real(kind=8) :: total_prob
    real(kind=8) :: delta_prob

    !-----------------------------------------
    write(*,*) ' Starting tcdf_most_likely...'

    ! Initialize flags
    lp = .false.
    li = .false.
    l_verbose = .false.
    lx = .false.
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
          if( index(trim(arg_flag),'-P') .ne. 0 ) then
             lp = .true.

             ! Read in the probability number:
             call get_command_argument(arg_count, arg_real, arg_len, exitcode)
             call checkexit(exitcode)
             arg_count = arg_count + 1
             read(arg_real,'(f10.4)') prob_contour
             write(*,*) 'Will search for ',prob_contour,' most likely contour.'
              
          elseif ( index(trim(arg_flag),'-v') .ne. 0 ) then
             ! Verbose mode
             l_verbose = .true.

          elseif ( index(trim(arg_flag),'-X') .ne. 0 ) then
             ! xyz output
             lx = .true.

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
             write(*,*) 'Open input netcdf file...'
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
       write(*,*) 'Will search for ',prob_contour,' contour'
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
       
       ! Set reading and writing limits to get histogram, later
       readst4d(ii) = 1
       readct4d(ii) = nc(ii)

       writest4d(ii) = 1
       writect4d(ii) = nc(ii)
    enddo

    ! Although it is possible get more generality later,
    ! for our purposes here, we want to allow for two cases:
    ! 1) a 2D only histogram -> first two dimensions in tcdfhist.
    ! 2) a 3D only histogram with the last dimension in time
    ! Check that the fourth dimension size is 1.  If not, panic.
    if ( nc(4) .ne. 1 ) then
       write(*,*) 'ERROR: 4th input dimension must be length 1.'
       write(*,*) 'Quit.'
       stop
    endif

    ! Allocate space for hist, mask, and prob
    if ( l_verbose ) write(*,*) 'Allocating work space...'
    allocate(hist(nc(1),nc(2),nc(3),nc(4)))
    allocate(mask(nc(1),nc(2),nc(3),nc(4)))
    allocate(hist_mask(nc(1),nc(2)))
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
    
    if ( lx ) then
       ! No need to create netcdf output file
    else
       ! Create output file
       if ( l_verbose ) write(*,*) 'Create output file...'
       mask_fid = nccre('ml_mask.nc',ncclobber,exitcode)

       do ii = 1,nax_max
          ! Create dimensions
          if ( l_verbose ) write(*,*) 'Create output dimensions...'
          vdims4d(ii) = ncddef(mask_fid,dnam(ii),nc(ii),exitcode)

          ! Create dimension variables
          if ( l_verbose ) write(*,*) 'Create output dimension variables...'
          dim_var_id(ii) = ncvdef(mask_fid,dnam(ii),ncfloat,1,vdims4d(ii),exitcode)
       enddo
       
       ! Create mask variable
       if ( l_verbose ) write(*,*) 'Create output mask variable...'
       mask_vid = ncvdef(mask_fid,'mask',ncint,4,vdims4d,exitcode)
       
       ! Exit define mode
       if ( l_verbose ) write(*,*) 'Exit define mode...'
       call ncendf(mask_fid,exitcode)
       
       ! Copy dimension variables over
       !do ii = 1,imax
       !   write(*,*) ivd(ii)
       !enddo
       do ii = 1,nax_max
          writest1d(1) = 1
          writect1d(1) = nc(ii)
          call ncvpt(mask_fid,dim_var_id(ii),writest1d,writect1d,&
                     dim_var_val(ii,writest1d(1):writect1d(1)),exitcode)
       enddo

    endif

    !=======================================================
    ! COMPUTE CONTOUR
    !=======================================================
    ! Set mask all to zeros.
    mask = 0

    ! Determine the maximum number of iterations (picks
    ! every box)
    max_iterations = nc(1)*nc(2)
    
    ! Loop over all the instances of the 3D dimension
    ! (this is the time dimension).
    do ll = 1,nc(3)

       ! Determine the total number of counts in the histogram
       sum_hist = sum(hist(:,:,ll,1))

       ! If there are no data points, no need to compute
       ! mask.
       if ( sum_hist .eq. 0 ) then
          if ( l_verbose ) write(*,*) 'No points in this histogram.'
       else
          ! Compute
             
          ! Set the histogram search mask to all .true.
          hist_mask = .true.
          
          ! Initialize the accumulating probability
          total_prob = 0.0
          
          ! Loop over iterations
          do ii = 1,max_iterations

             !if ( l_verbose ) write(*,*) 'Iteration ',ii,' of ',max_iterations
             
             ! Find the location of the current biggest
             ! value in the histogram.
             ij = maxloc(hist(:,:,ll,1),hist_mask)
             
             ! Determine the incremental addition of this
             ! probability.
             delta_prob = dble(hist(ij(1),ij(2),ll,1))/dble(sum_hist)
             
             ! If the sum is less than the desired probability,
             if ( total_prob + delta_prob .lt. dble(prob_contour) ) then
                
                ! Add this value to the total probability
                total_prob = total_prob + delta_prob
                
                ! Remove this box from the search mask
                hist_mask(ij(1),ij(2)) = .false.
                
                ! And flag this box on the output mask
                mask(ij(1),ij(2),ll,1) = 1
                
                ! And keep on searching for the next highest value...
                
             else
                
                ! Here, the new total probability will be at or
                ! exceed the probability contour, so we want to stop.
                !
                ! However, the question is whether it is best to
                ! include or exclude this extra little bit of
                ! probability.
                !
                ! Here, we decide to include or exclude based on
                ! which operation will bring the final total probability
                ! closest to the prob_contour that is desired.
                
                if ( total_prob + delta_prob .eq. dble(prob_contour) ) then
                   ! Exact equality, no brainer, keep this bit and exit.
                   total_prob = total_prob + delta_prob
                   mask(ij(1),ij(2),ll,1) = 1
                   exit
                   
                elseif ( (dble(prob_contour) - total_prob) .lt. (total_prob + delta_prob - dble(prob_contour)) ) then
                   ! We are closer to prob_contour if we DO NOT include this
                   ! last extra little bit, so do not include it.
                   exit
                   
                else
                   ! We are closer to prob_contour if we DO include this
                   ! last extra little bit, so do include it.
                   total_prob = total_prob + delta_prob
                   mask(ij(1),ij(2),ll,1) = 1
                   exit
                   
                endif
             endif
          enddo   ! End of looping over each possible point.
          
          ! Check result
          if ( total_prob .lt. dble(prob_contour) .and. ii .eq. max_iterations ) then
             write(*,*) 'WARNING: Did not exceed prob_contour even'
             write(*,*) 'WARNING: after going through all boxes!!!'
             write(*,*) 'WARNING: total_prob = ',total_prob
             write(*,*) 'WARNING: prob_contour = ',prob_contour
          endif
          
          if ( l_verbose ) write(*,*) 'Total_prob = ',total_prob
          if ( l_verbose ) write(*,*) 'Num flagged boxes = ',sum(mask(:,:,ll,1))
          if ( l_verbose ) write(*,*) 'Mask number ',ll,' after ',ii,' iterations.'
       endif
    enddo
    
    !=======================================================
    ! WRITE OUTPUT
    !=======================================================

    ! If the user requests xyz output
    if ( lx ) then
       xyz_fid = 42
       open(xyz_fid,file='ml_mask.xyz',status='new',form='formatted')       
       do ii = 1,nc(1)
          do jj = 1,nc(2)
             write(xyz_fid,'(f20.10,x,f20.10,x,i1)') dim_var_val(1,ii),dim_var_val(2,jj),mask(ii,jj,1,1)
          enddo
       enddo
    else
       call ncvpt(mask_fid,mask_vid,writest4d,writect4d,mask,exitcode)
       call ncclos(mask_fid,exitcode)
    endif

    !=======================================================
    ! CLEAN UP
    !=======================================================
    if ( l_verbose ) write(*,*) 'Closing input file...'
    call ncclos(hist_fid,exitcode)

    deallocate(dim_var_val)
    deallocate(hist)
    deallocate(mask)
    deallocate(hist_mask)

    write(*,*) 'tcdf_most_likely done.'

!=======================================
  end program tcdf_most_likely
!=======================================

  subroutine print_help()
    write(*,*) ' '
    write(*,*) 'tcdf_most_likely [-P <(0.0,1.0)>]  -I <file.nc>'
    write(*,*) '                 [-v] [-h] [-X]'
    write(*,*) ' '
    write(*,*) 'This module will compute the most'
    write(*,*) 'likely probability contour given'
    write(*,*) 'a histogram.'
    write(*,*) ' '
    write(*,*) 'Currently, the histogram is called'
    write(*,*) 'hist in a netCDF file, more flexibility'
    write(*,*) 'can be incorporated later.'
    write(*,*) ' '
    write(*,*) ' -P = specify the probability level of the'
    write(*,*) '      most likely contour, e.g. -P 95 will'
    write(*,*) '      result in the 95% most likely contour.'
    write(*,*) '      DEFAULT: 0.95'
    write(*,*) ' '
    write(*,*) ' -I = specify the name of the input file.'
    write(*,*) ' '
    write(*,*) ' -v = verbose mode.'
    write(*,*) ' '
    write(*,*) ' -h = print this message and quit.'
    write(*,*) ' '
    write(*,*) ' -X = output in .xyz in stead of .nc'
    write(*,*) ' '
    write(*,*) ' TO ADD LATER: -L <iteration_limit>, '
    write(*,*) '               -E <tolerance>,'
    write(*,*) '               -V <variable_name>' 
    write(*,*) '               -O <output_file_name>'
    write(*,*) ' '
  end subroutine print_help
