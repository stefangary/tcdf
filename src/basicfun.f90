! This software is distributed under the terms
! of the GNU General Public License v3 or any
! later version.
! Copyright Stefan Gary, 2018.
!---------------------------------------------

      module basicfun

        implicit none

      contains

!---------------------------------------------
!
! olavg, olavgvar, checkexit, check_limits
!
! init_random_seed, binsearch, binsearch2
!
!---------------------------------------------

!---------------------------------------------------------
    subroutine olavg(new_val,counter,store_avg)

!---------------------------------------------------------
! This subroutine implements an online averaging algorithm
! by incorporating each new_val into the updated variable
! store_avg.  The variable counter is needed to keep track
! of the number of points that are not fill values.  See
! the Technical Notes in subroutine olavgvar for more info.
!---------------------------------------------------------
          
        use params

        implicit none
          
        ! Declare local variables
        real, intent(in) :: new_val
        integer, intent(inout) :: counter
        real, intent(inout) :: store_avg
        real :: delta
        
        integer :: new_count
        if(new_val.eq.fill_real) then
           ! New value flagged blank, do not augment counter
           ! and skip in the averaging process.
        else
           ! Value is ok, compute the means and augment the counter
           new_count = counter + 1
             
           ! Original, old method for computing online average:
           !delta = new_val - store_avg
           !store_avg = store_avg + delta/real(counter)
         
           ! New method, eliminates any subtraction errors
           ! Note that initial value of store_avg is not
           ! important - the zero multiplication will
           ! eliminate it on the first iteration.
           store_avg = (real(counter)*store_avg + new_val)/real(new_count)
           counter = new_count
        endif
          
        return
      end subroutine olavg

!---------------------------------------------------------

!---------------------------------------------------------

    subroutine olavgvar(new_val,counter,store_avg,store_var)

!---------------------------------------------------------
! This subroutine implements an online averaging algorithm
! by incorporating each new_val into the updated variable
! store_avg.  The variance is also computed in the variable
! store_var.  The variable counter is needed to keep track
! of the number of points that are incorporated into the
! average (i.e. not fill values).
!
!--------------------Technical Details--------------------
! Use an online algorithm based on Knuth, D. E. (1998).
! _The_Art_of_Computer_Programming_, vol. 2:
! Seminumerical Algorithms, 3rd Ed. p232.
! Boston: Addison-Wesley.  ...which I found out
! about on wikipedia,
! http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
! See research notes for verification of algorithm.
! In the testing of this code, it is clear that the
! subtraction in the averaging does not cause the
! loss of digits for the numbers I'm using.
!
! After looping (and pushing the stack), the last step is
! to normalize the variance (store_var) by n (otherwise,
! our var values represent the sum of the squares of the
! differences between each value and its mean), but not
! the variance.  This must be done separately.
!
! In addition to checking for fill_values, we also check
! for NaN's.  This is a non-negligable increase in the
! computational cost.
!---------------------------------------------------------

          use params

          implicit none

          ! Declare local variables
          real, intent(in) :: new_val
          integer, intent(inout) :: counter
          real, intent(inout) :: store_avg
          real, intent(inout) :: store_var
          real(kind=8) :: delta
          
          integer :: new_count
          real(kind=8) :: dnew_val
          real(kind=8) :: dstore_avg
          real(kind=8) :: dstore_var
        
          if(new_val.eq.fill_real) then
             ! New value flagged blank, do not augment counter
             ! and skip in the averaging process.
          else
             ! Value is ok, compute the means and augment the counter
             new_count = counter + 1
             
             ! Convert to doubles
             dnew_val = dble(new_val)
             dstore_avg = dble(store_avg)
             dstore_var = dble(store_var)
             
             delta = dnew_val - dstore_avg
             
             store_avg = (real(counter)*store_avg + new_val)/real(new_count)
             counter = new_count
             
             dstore_avg = dble(store_avg)
             
             store_var = sngl(dstore_var + delta*(dnew_val - dstore_avg))
          endif
          
          return
        end subroutine olavgvar
!---------------------------------------------------------

!------------------------------------------------

     subroutine checkexit(code)

!------------------------------------------------
! This subroutine checks the output of the
! standard get_command_argument subroutine
! and will stop the code if something is
! wrong.
!------------------------------------------------

       integer :: code

       if (code .gt. 0) then
          write(*,*) ' Error in retrieving argument!'
          stop
       elseif(code .eq. -1) then
          write(*,*) ' Warning: argument truncated while passed!'
          stop
       elseif(code .eq. 0) then
#ifdef verbose
          write(*,*) ' Successful read from command line.'
#endif
       else
          write(*,*) ' Unknown exitcode from get_command_argument!'
          stop
       endif

       return

     end subroutine checkexit

!------------------------------------------------

!------------------------------------------------

     subroutine check_limits(low,high)

!------------------------------------------------
! This subroutine will check that low is less
! then high.  If not, then the values will be
! switched.
!------------------------------------------------

       implicit none

       real :: low,high,hold

       if ( low .lt. high ) then
          ! This is the way things should be
          ! so do nothing.  Note that nothing
          ! changes for two identical values.
       else
          hold = high
          high = low
          low = hold
       endif

       return

     end subroutine check_limits

!------------------------------------------------

!------------------------------------------------

     subroutine check_limits_i(low,high)

!------------------------------------------------
! This subroutine will check that low is less
! then high.  If not, then the values will be
! switched.  This is for integer inputs!
!------------------------------------------------

       implicit none

       integer :: low,high,hold

       if ( low .lt. high ) then
          ! This is the way things should be
          ! so do nothing.  Note that nothing
          ! changes for two identical values.
       else
          hold = high
          high = low
          low = hold
       endif

       return

     end subroutine check_limits_i

!------------------------------------------------

!---------------------------------------

      subroutine init_random_seed()

!---------------------------------------
! This subroutine initializes the
! random seed based on the system clock
! state.  Based on the gfortran manual.
!---------------------------------------

        integer :: i, n, clock
        integer, dimension(:), allocatable :: seed

        call random_seed(size = n)
        allocate(seed(n))

        call system_clock(count=clock)

        seed = clock + 37 * (/ (i-1, i=1, n) /)
        call random_seed(put=seed)

        deallocate(seed)

      end subroutine init_random_seed

!---------------------------------------

!----------------------------------------

      recursive function binsearch(rarray,rvalue,ileft,iright) result(ibox)

!----------------------------------------
! This function will return the integer
! index of the real number, rvalue, that
! lies within the bounds (edges) of the
! elements within the monotonically
! increasing array rarray.
!
! The subroutine uses a recursive binary
! (divide and conquor) search algorithm.
!
! Since this is a recursive algorithm,
! ileft and iright must be specified as
! the left and right bounds.
!
! This function returns -1 if the value
! is not within the array.
!
!-----------------------------------------

        implicit none

        real :: rarray(*)
        real :: rvalue
        integer :: ileft, iright, ibox, imid
        
        ! Check that we have not exhausted the search.
        if (iright .le. ileft) then
           write(*,*) 'WARNING: Binary search unsucessful!'
           write(*,*) 'R: ',iright,' L: ',ileft
           ibox = -1
           return
        endif

        ! Check that we have found the item.
        if ( (iright - ileft) .eq. 1 ) then
           if ( (rarray(iright).ge.rvalue).and.(rarray(ileft).le.rvalue) ) then
              ibox = ileft
              return

           else
              !write(*,*) 'WARNING: Value not found.'
              ibox = -1
              return

           endif
        endif

        ! Continue to check a subset of the array.
        ! Compute the index of midpoint.
        ! NOTE: INTEGER DIVISION => TRUNCATION
        imid = ileft + ((iright - ileft) / 2) 

        ! Check on which side we fall.
        if ( rvalue .gt. rarray(imid) ) then
           ! We fall to the right
           ibox = binsearch(rarray,rvalue,imid,iright)
           return

        elseif ( rvalue .le. rarray(imid) ) then
           ! We fall to the left
           ibox = binsearch(rarray,rvalue,ileft,imid)
           return
           
        else
           ! We have exhausted the search somehow.
           write(*,*) 'WARNING: Binary search unsucessful!'
           write(*,*) 'R: ',iright,' L: ',ileft, ' M: ',imid
           ibox = -1
           return
        endif
        
      end function binsearch
       
!----------------------------------------

!----------------------------------------

      recursive function binsearch2(rarray,dvalue,ileft,iright) result(ibox)

!----------------------------------------
! This function will return the integer
! index of the double prec. number, dvalue, that
! lies within the bounds (edges) of the
! elements within the monotonically
! increasing array rarray.
!
! The subroutine uses a recursive binary
! (divide and conquor) search algorithm.
!
! Since this is a recursive algorithm,
! ileft and iright must be specified as
! the left and right bounds.
!-----------------------------------------

        implicit none

        real :: rarray(*)
        real(kind=8) :: dvalue
        integer :: ileft, iright, ibox, imid
        
        ! Check that we have not exhausted the search.
        if (iright .le. ileft) then
           write(*,*) 'WARNING: Binary search unsucessful!'
           write(*,*) 'R: ',iright,' L: ',ileft
           ibox = -1
           return
        endif

        ! Check that we have found the item.
        if ( (iright - ileft) .eq. 1 ) then
           if ( (dble(rarray(iright)).ge.dvalue).and.(dble(rarray(ileft)).le.dvalue) ) then
              ibox = ileft
              return

           else
              !write(*,*) 'WARNING: Value not found.'
              ibox = -1
              return

           endif
        endif

        ! Continue to check a subset of the array.
        ! Compute the index of midpoint.
        ! NOTE: INTEGER DIVISION => TRUNCATION
        imid = ileft + ((iright - ileft) / 2) 

        ! Check on which side we fall.
        if ( dvalue .gt. dble(rarray(imid)) ) then
           ! We fall to the right
           ibox = binsearch2(rarray,dvalue,imid,iright)
           return

        elseif ( dvalue .le. dble(rarray(imid)) ) then
           ! We fall to the left
           ibox = binsearch2(rarray,dvalue,ileft,imid)
           return
           
        else
           ! We have exhausted the search somehow.
           write(*,*) 'WARNING: Binary search unsucessful!'
           write(*,*) 'R: ',iright,' L: ',ileft, ' M: ',imid
           ibox = -1
           return
        endif
        
      end function binsearch2
       
!----------------------------------------

      end module basicfun
