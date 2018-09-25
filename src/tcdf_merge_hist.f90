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

  ! Read the command line and figure out
  ! how many files are listed.
  num_command_arg = command_argument_count()
       
  if(num_command_arg .lt. 1) then
     write(*,*) ' Error: Must specify at least 1 input file!'
     call print_help()
     stop
  else
     arg_count = 1

     ! Loop through all the command line files.
     do
        ! Test that we are still reading all the files.
        if ( arg_count .gt. num_command_arg ) exit

        call get_command_argument(arg_count, arg_string, arg_len, exitcode)
        call checkexit(exitcode)
        arg_count = arg_count + 1

        write(*,*) arg_string
     enddo
  endif
end program tcdf_hist_merge

subroutine print_help()

  write(*,*) 'Placeholder for printing help.'
  
end subroutine print_help
