!--------------------------------------
! Merge many histogram files by adding
! them together.
!--------------------------------------
! This software is distributed under the
! terms of the GNU GPL v3 or any later
! version.
!--------------------------------------
program tcdf_hist_merge

  use basicfun
  use netcdfio
  
  implicit none

  ! Command line args, File names.
  integer :: num_command_arg
  integer :: arg_len, arg_count, file_count
  character(len=50) :: arg_string
  character(len=10) :: arg_int
  character(len=5) :: arg_flag
  character(len=10) :: arg_real

  ! Default output file name
  character(len=50) :: merged_hist_file_name

  ! Whether or not to get khist variable too
  logical :: lk, l_verbose

  ! File information
  integer :: hist_vid, khist_vid, hist_fid
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
  integer, allocatable :: hist_add(:,:,:,:)
  integer, allocatable :: hist_out(:,:,:,:)
  integer, allocatable :: khist_add(:,:,:,:)
  integer, allocatable :: khist_out(:,:,:,:)
  ! dim_var_val contains the dimension variables.
  ! Accessed by dim,value_along_dimension
  real, allocatable :: dim_var_val(:,:)
  integer :: ii, jj, ll, sum_hist
  integer :: ij(2)
  
  ! Read the command line and figure out
  ! how many files are listed.
  num_command_arg = command_argument_count()
  lk = .false.
  merged_hist_file_name = 'merged_hist.nc'
  
  if(num_command_arg .lt. 1) then
     write(*,*) ' Error: Must specify at least 1 flag!'
     call print_help()
     stop
  else
     arg_count = 1
     file_count = 1

     ! Loop through all the command line flags.
     do
        ! Test that we are still reading all the flags.
        if ( arg_count .gt. num_command_arg ) exit
        
        call get_command_argument(arg_count, arg_flag, arg_len, exitcode)
        call checkexit(exitcode)
        arg_count = arg_count + 1

        if ( index(trim(arg_flag),'-F') .ne. 0 ) then
           ! Start to read the file list

           ! Loop through all the input files
           do
              ! Test that we are still reading all the files.
              if ( arg_count .gt. num_command_arg ) exit
              call get_command_argument(arg_count, arg_string, arg_len, exitcode)
              call checkexit(exitcode)
              arg_count = arg_count + 1
        
              if ( l_verbose) write(*,*) 'Working on file: ',trim(arg_string),' which is arg_count ',arg_count

              ! Open the file
              if ( l_verbose) write(*,*) 'Open input netcdf file...'
              hist_fid = ncopn(trim(arg_string),ncnowrit,exitcode)

              ! Get variable ID
              if ( l_verbose ) write(*,*) 'Getting variable ID...'
              hist_vid = ncvid(hist_fid,'hist',exitcode)
              
              ! For the first file, allocate space for
              ! the storage variables and also create
              ! the output file.
              if ( file_count .eq. 1 ) then
                 if ( l_verbose ) write(*,*) 'Reading first file to init storage and output...'

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
                 
                 if ( l_verbose ) write(*,*) 'Allocate space...'
                 ! Initialize storage to zeros, see below for why.
                 allocate(hist_out(nc(1),nc(2),nc(3),nc(4)))
                 allocate(hist_add(nc(1),nc(2),nc(3),nc(4)))
                 hist_out = 0.0
                 hist_add = 0.0
                 allocate(dim_var_val(nax_max,nc_max))
                 dim_var_val = 0.0
                 if ( lk ) then
                    allocate(khist_out(nc(1),nc(2),nc(3),nc(4)))
                    allocate(khist_add(nc(1),nc(2),nc(3),nc(4)))
                    khist_out = 0.0
                    khist_add = 0.0
                 endif

                 if ( l_verbose ) write(*,*) 'Read dimension variables...'
                 do ii = 1,nax_max
                    d_vid = ncvid(hist_fid,trim(dnam(ii)),exitcode)
                    call ncvgt(hist_fid,d_vid,1,nc(ii),dim_var_val(ii,1:nc(ii)),exitcode)
                 enddo
              else
                 ! No set up to be done for other input files.
              endif

              ! Get histogram to add to the output
              if ( l_verbose ) write(*,*) 'Loading histogram...'
              call ncvgt(hist_fid,hist_vid,readst4d,readct4d,hist_add,exitcode)

              ! For all input files, add histogram counts
              ! into the storage histograms.  Since the first
              ! file is added to zeros in the storage, this
              ! means that if only one file is specified as input,
              ! the output will be the same as the input file.
              ! While this may not be the most efficient way of
              ! initializing an aggregation, it works well if
              ! you are limited on space and can only work with
              ! a small subset of uncompressed files.

              ! Add hist_add to hist_out
              if ( l_verbose ) write(*,*) 'Adding hists...'
              hist_out = hist_out + hist_add

              ! For the case where there is a khist, we also
              ! need to load and add that data to the khist_out
              if ( lk ) then
                 ! Get variable ID
                 if ( l_verbose ) write(*,*) 'Getting KHIST variable ID...'
                 hist_vid = ncvid(hist_fid,'khist',exitcode)

                 ! Get histogram to add to the output
                 if ( l_verbose ) write(*,*) 'Loading KHIST histogram...'
                 call ncvgt(hist_fid,hist_vid,readst4d,readct4d,khist_add,exitcode)

                 ! Add khist_add to khist_out
                 if ( l_verbose ) write(*,*) 'Adding KHISTS...'
                 khist_out = khist_out + khist_add
                 
              endif

              ! Close current input file.
              call ncclos(hist_fid,exitcode)

              ! Move to the next file
              file_count = file_count + 1
              
           enddo
        elseif ( index(trim(arg_flag),'-h') .ne. 0 ) then
           ! Print help and exit
           call print_help()
           stop
        elseif ( index(trim(arg_flag),'-O') .ne. 0 ) then
           ! Specify the name of the output file and
           ! change from default.
           call get_command_argument(arg_count, arg_string, arg_len, exitcode)
           call checkexit(exitcode)
           arg_count = arg_count + 1
           write(*,*) 'Changing default file name: ',merged_hist_file_name
           merged_hist_file_name = arg_string
           write(*,*) 'to: ',merged_hist_file_name
        elseif ( index(trim(arg_flag),'-K') .ne. 0 ) then
           ! Also sum up khist variables
           lk = .true.
        elseif ( index(trim(arg_flag),'-V') .ne. 0 ) then
           ! Print lots of output
           l_verbose = .true.
        else
           write(*,*) 'ERROR: tcdf_hist_merge unrecognized flag'
        endif
     enddo
  endif

  ! We're basically done - data from files were added
  ! to storage as their names were read from the
  ! command line.  All we need to do is create an
  ! output file and write the data to that file.
  if ( l_verbose ) write(*,*) 'Creating output file...'
  hist_fid = nccre(trim(merged_hist_file_name),ncclobber,exitcode)

  if ( l_verbose ) write(*,*) 'Create output dimensions and variables...'
  do ii = 1,nax_max
     ! Create dimensions
     vdims4d(ii) = ncddef(hist_fid,dnam(ii),nc(ii),exitcode)
     
     ! Create dimension variables
     dim_var_id(ii) = ncvdef(hist_fid,dnam(ii),ncfloat,1,vdims4d(ii),exitcode)
     if ( index(trim(dnam(ii)),'dep') .ne. 0 ) then
        ! Axis is positive down
        call ncaptc(hist_fid,dim_var_id(ii),'positive',ncchar,4,'down',exitcode)
     endif
     if ( index(trim(dnam(ii)),'rho') .ne. 0 ) then
        ! Axis is positive down
        call ncaptc(hist_fid,dim_var_id(ii),'positive',ncchar,4,'down',exitcode)
     endif
  enddo
       
  ! Create hist variable
  if ( l_verbose ) write(*,*) 'Create output hist variable...'
  hist_vid = ncvdef(hist_fid,'hist',ncint,4,vdims4d,exitcode)

  if ( lk ) then
     if ( l_verbose ) write(*,*) 'Creates output KHIST variable...'
     khist_vid = ncvdef(hist_fid,'khist',ncint,4,vdims4d,exitcode)
  endif
  
  ! Exit define mode
  if ( l_verbose ) write(*,*) 'Exit define mode...'
  call ncendf(hist_fid,exitcode)
       
  ! Copy dimension variables over
  if ( l_verbose ) write(*,*) 'Copy dimension variables over...'
  do ii = 1,nax_max
     writest1d(1) = 1
     writect1d(1) = nc(ii)
     call ncvpt(hist_fid,dim_var_id(ii),writest1d,writect1d,&
          dim_var_val(ii,writest1d(1):writect1d(1)),exitcode)
  enddo

  if ( l_verbose ) write(*,*) 'Writing histogram to output file...'
  call ncvpt(hist_fid,hist_vid,writest4d,writect4d,hist_out,exitcode)

  if ( lk ) then
     if ( l_verbose ) write(*,*) 'Writing KHIST to output file...'
     call ncvpt(hist_fid,khist_vid,writest4d,writect4d,khist_out,exitcode)
  endif

  if ( l_verbose ) write(*,*) 'Closing output file...'
  call ncclos(hist_fid,exitcode)
  
  if ( l_verbose ) write(*,*) 'tcdf_hist_merge done.'
end program tcdf_hist_merge

subroutine print_help()

  write(*,*) '---------------------------------'
  write(*,*) 'tcdf_merge_hist'
  write(*,*) '---------------------------------'
  write(*,*) 'usage:'
  write(*,*) 'tcdf_merge_hist [-O] [-h] [-K] [-V] -F file1.nc file2.nc ...'
  write(*,*) '---------------------------------'
  write(*,*) 'Histogram variables within each'
  write(*,*) 'file in this list will be added up'
  write(*,*) 'into a single histogram, a new'
  write(*,*) 'output file [merged_hist.nc] or'
  write(*,*) 'name specified by -O. -h prints'
  write(*,*) 'this help message. -F simply tells'
  write(*,*) 'the command line interpreter that'
  write(*,*) 'a list of files is coming and there'
  write(*,*) 'are no more flags/options after'
  write(*,*) 'this list.'
  write(*,*) ' '
  write(*,*) 'This program assumes all histograms'
  write(*,*) 'in all input files will have exactly'
  write(*,*) 'the same dimensions and variable names.'
  write(*,*) 'If not, then it simply crashes.'
  write(*,*) 'This was written to work with default'
  write(*,*) 'output from tcdfhist.'
  write(*,*) ' '
  write(*,*) 'Stefan Gary, 2018, distributed under'
  write(*,*) 'the terms of the GNU GPL v3 or later.'
  write(*,*) '----------------------------------'
  
end subroutine print_help
