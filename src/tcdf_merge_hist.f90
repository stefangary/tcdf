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
  integer :: arg_len, arg_count
  character(len=50) :: arg_string
  character(len=10) :: arg_int
  character(len=5) :: arg_flag
  character(len=10) :: arg_real

  ! Default output file name
  character(len=50) :: merged_hist_file_name
  merged_hist_file_name = 'merged_hist.nc'
  
  ! Read the command line and figure out
  ! how many files are listed.
  num_command_arg = command_argument_count()
       
  if(num_command_arg .lt. 1) then
     write(*,*) ' Error: Must specify at least 1 flag!'
     call print_help()
     stop
  else
     arg_count = 1

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
        
              write(*,*) 'Working on file: ',trim(arg_string),' which is arg_count ',arg_count

              ! For the first file, allocate space for
              ! the storage variables and also create
              ! the output file.
              if ( arg_count .eq. 1 ) then
                 write(*,*) 'Reading first file to init storage and output...'
                 ! WORKING HERE
                 ! Initialize storage to zeros, see below for why.
              else
                 ! Nothing to be done for all subsequent
                 ! files.
              endif

              ! For all input files, add histogram counts
              ! into the storage histograms.  Since the first
              ! file is added to zeros in the storage, this
              ! means that if only one file is specified as input,
              ! the output will be the same as the input file.
              ! While this may not be the most efficient way of
              ! initializing an aggregation, it works well if
              ! you are limited on space and can only work with
              ! a small subset of uncompressed files.
              
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
           
        else
           write(*,*) 'ERROR: tcdf_hist_merge unrecognized flag'
        endif
     enddo
  endif
end program tcdf_hist_merge

subroutine print_help()

  write(*,*) '---------------------------------'
  write(*,*) 'tcdf_merge_hist'
  write(*,*) '---------------------------------'
  write(*,*) 'usage:'
  write(*,*) 'tcdf_merge_hist [-O] [-h] -F file1.nc file2.nc ...'
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
