! This software is distributed under the
! terms of the GNU General Public License
! v3 or any later version.
! Copyright Stefan Gary, 2018.
!------------------------------------------
! Variable names and default values for
! Lagrangian variables.
!
! Should be used in conjunction with
! gridsio.
!
! Need to do some cleaning up and
! organization with this file and gridsio.
!------------------------------------------

      module tcdfio

        implicit none

        ! Dimension ID and names
        integer :: trajdid,timedid
        integer :: vtype, nvdims, nvatts
        integer :: vdims(2)
        integer :: vdims_out(2)

        ! First (i) dimension is time and
        ! Second (j) dimension is traj (unlim) 
        character(len=4) :: trajdnam = 'traj'
        character(len=4) :: timednam = 'time'
        character(len=50) :: dummy

        ! Declare sizes of variables
        ! npts = initial number of points
        integer :: ntraj, npts, nall

        ! Read and write index limits
        integer :: lag_readst2d(2)
        integer :: lag_readct2d(2)
        integer :: lag_writest2d(2)
        integer :: lag_writect2d(2)

        ! Specific values for scale_factors and
        ! add_offsets.
        ! Scale is by default 1
        real :: scale_factor(1) = 1.0
        real, parameter :: xy_scale(1) = 0.01
        real, parameter :: ts_scale(1) = 0.001
        real, parameter :: z_scale(1) = 0.1
        ! Add offset is by default 0.
        real :: add_offset(1) = 0.0
        real, parameter :: zero_offset(1) = 0.0
        real, parameter :: xy_offset(1) = 0.0
        real, parameter :: z_offset(1) = 3266.0
        real, parameter :: s_offset(1) = 30.0
        real, parameter :: year_offset(1) = 1900.0

        ! Variable ID and names
        character(len=4) :: lamvnam = 'lam '
        character(len=4) :: phivnam = 'phi '
        character(len=4) :: depvnam = 'dep '
        character(len=4) :: tempvnam = 'temp'
        character(len=4) :: saltvnam = 'salt'
        character(len=4) :: rhovnam = 'rho '
        character(len=4) :: uvnam = 'u   '
        character(len=4) :: vvnam = 'v   '
        character(len=4) :: wvnam = 'w   '
        character(len=4) :: qvnam = 'q   '
        character(len=4) :: l_yvnam = 'l_y '
        character(len=4) :: l_mvnam = 'l_m '
        character(len=4) :: l_dvnam = 'l_d '
        character(len=4) :: drdz_vnam = 'drdz'
        character(len=4) :: dtdz_vnam = 'dtdz'
        character(len=4) :: year_vnam = 'year'
        character(len=4) :: month_vnam = 'mon '
        character(len=4) :: day_vnam = 'day '
        character(len=4) :: agevnam = 'age '
        character(len=4) :: id_vnam = 'id  '
        character(len=4) :: rldvnam = 'rld '
        character(len=4) :: bdepvnam = 'bdep'
        ! iit, jjt, kkt are decimal index values for
        ! location of particle in i,j,k T-grid.
        character(len=4) :: iitvnam = 'iit '
        character(len=4) :: jjtvnam = 'jjt '
        character(len=4) :: kktvnam = 'kkt '
        character(len=4) :: kbotvnam = 'kbot'
        
        ! Layer names are set to defaults here.
        ! Can have multiple layer names in one
        ! file (top1, top2, top3, etc.) and those
        ! will be changed on the fly depending
        ! on command line flags.
        character(len=4) :: layer_top_vnam = 'top1'
        character(len=4) :: layer_bot_vnam = 'bot1'
        
        ! Variable mask values (mask values
        ! should not trip the range limits?)
        ! Note that you will need to be careful
        ! that the mask values fall within the
        ! possible range of the packed variables
        ! This is not so difficult with all variables
        ! EXCEPT FOR DEPTH (AND BOTTOM DEPTH) whose
        ! packing needs an add_offset to allow for
        ! sufficient range in the variable to be
        ! contained in a 5-digit short integer.
        ! Packed variables are stored as signed shorts,
        ! i.e. -32,768 to +32,768
        ! For lon and lat, this means that a scale factor
        ! of 0.01 and add_offset of 0 will give you a range
        ! from -327.68 to +327.68.
        ! Because out = stored_value*scale_factor + add_offset
        ! For depth and bottom depth, the range is from
        ! -3276.8+z_offset to +3276.8 + z_offset
        ! = -10.8 m to 6542.8 m
        ! which will cover the global ocean.  If you want
        ! to include the Mariana Trench, then you'll need
        ! to set the z_scale = 1.0 and z_offset = 32660.0
        ! and TEST THAT IT WORKS!
        real, parameter :: lam_mask = -300.0
        real, parameter :: phi_mask = -300.0        
        real, parameter :: dep_mask = -10.5
        real, parameter :: u_mask = -9.0
        real, parameter :: v_mask = -9.0
        real, parameter :: w_mask = -9.0
        real, parameter :: t_mask = -9.0
        real, parameter :: s_mask = -9.0
        real, parameter :: r_mask = -9.0
        real, parameter :: dtdz_mask = -999.0
        real, parameter :: drdz_mask = -999.0
        real, parameter :: q_mask = -999.0
        real, parameter :: year_mask = -9.0
        real, parameter :: month_mask = -9.0
        real, parameter :: day_mask = -9.0
        real, parameter :: age_mask = -9.0
        real, parameter :: layer_top_mask = -9.0
        real, parameter :: layer_bot_mask = -9.0
        real, parameter :: rld_mask = -9.0
        real, parameter :: lay_mask = -9.0
        real, parameter :: bdep_mask = -9.0
        real, parameter :: kbot_mask = -9
        real, parameter :: iit_mask = -9
        real, parameter :: jjt_mask = -9
        real, parameter :: kkt_mask = -9
        
        ! Variable default values (when
        ! initialized, to allow for never
        ! being edited out when the default
        ! variable limits are being used.
        real, parameter :: lam_def = 0.0
        real, parameter :: phi_def = 0.0        
        real, parameter :: dep_def = 0.0
        real, parameter :: u_def = 0.0
        real, parameter :: v_def = 0.0
        real, parameter :: w_def = 0.0
        real, parameter :: t_def = 1.0
        real, parameter :: s_def = 1.0
        real, parameter :: r_def = 1.0
        real, parameter :: dtdz_def = 0.0
        real, parameter :: drdz_def = 0.0
        real, parameter :: q_def = 0.0
        real, parameter :: year_def = 1.0
        real, parameter :: month_def = 1.0
        real, parameter :: day_def = 1.0
        real, parameter :: age_def = 0.0
        real, parameter :: layer_top_def = 0.0
        real, parameter :: layer_bot_def = 0.0
        real, parameter :: rld_def = 0.0
        real, parameter :: bdep_def = 0.0
        real, parameter :: kbot_def = 0
        real, parameter :: iit_def = 0
        real, parameter :: jjt_def = 0
        real, parameter :: kkt_def = 0
        
        ! Variable default ranges (values
        ! that the _1 and _2 values are
        ! set to before command line
        ! options are read).  The _def
        ! variables, above, should fit
        ! in each corresponding range.
        real, parameter :: lam_min = -180.0
        real, parameter :: lam_max =  180.0
        real, parameter :: phi_min = -90.0
        real, parameter :: phi_max =  90.0
        real, parameter :: dep_min = -10.0
        real, parameter :: dep_max =  9000.0
        real, parameter :: u_min =  -10.0
        real, parameter :: u_max =   10.0
        real, parameter :: v_min =  -10.0
        real, parameter :: v_max =   10.0
        real, parameter :: w_min =  -10.0
        real, parameter :: w_max =   10.0
        real, parameter :: t_min =  -10.0
        real, parameter :: t_max =   40.0
        real, parameter :: s_min =  -10.0
        real, parameter :: s_max =   40.0
        real, parameter :: r_min =  -10.0
        real, parameter :: r_max =   50.0
        real, parameter :: q_min = -50000.0
        real, parameter :: q_max =  50000.0
        real, parameter :: dtdz_min = -50000.0
        real, parameter :: dtdz_max =  50000.0
        real, parameter :: drdz_min = -50000.0
        real, parameter :: drdz_max =  50000.0
        real, parameter :: year_min = -10.0
        real, parameter :: year_max =  50000.0
        real, parameter :: mon_min = -10.0
        real, parameter :: mon_max =  12.0
        real, parameter :: day_min = -10.0
        real, parameter :: day_max =  31.0
        real, parameter :: age_min = -10.0
        real, parameter :: age_max = 32000.0
        real, parameter :: layer_top_max = 10000.0
        real, parameter :: layer_top_min = 0.0
        real, parameter :: layer_bot_max = 10000.0
        real, parameter :: layer_bot_min = 0.0
        real, parameter :: rld_min = 0.0
        real, parameter :: rld_max = 1.0
        real, parameter :: bdep_min = -1.0
        real, parameter :: bdep_max = 9000.0
        real, parameter :: kbot_min = 0
        real, parameter :: kbot_max = 9000
        real, parameter :: iit_min = 0
        real, parameter :: iit_max = 9000
        real, parameter :: jjt_min = 0
        real, parameter :: jjt_max = 9000
        real, parameter :: kkt_min = 0
        real, parameter :: kkt_max = 9000
        
      end module tcdfio

!------------------------------------------

!------------------------------------------

      module load_tcdf

        implicit none

      contains

        ! All functions for reading and
        ! writing tcdf (trajectory
        ! ensembles in netcdf) files.

        ! get_lag_var
        ! get_lag_var_seg
        ! put_lag_var
        ! dup_lag_var
        ! put_lag_var_seg
        ! put_lag_var_seg_pad
        ! cre_lag_var
!------------------------------------------------

     subroutine get_lag_var(fid,vn,ro)

!------------------------------------------------
! This subroutine will load a Lagrangian variable
! from a tcdf formatted file into a the real
! variable ro given,
! fid   = file ID
! vn    = variable name
! ro    = real output variable
!
! This subroutine will automatically search for
! the type of variable in the source file and
! convert that source variable (including unpack)
! to a real number array.
!
! Quite a bit of additional information is being
! passed through the lag_readst2d, lag_readct2d, npts,
! variables which are in tcdfio.
! This is done to make call lines simple and
! keep things that are constant out of the way.
!------------------------------------------------

       use tcdfio
       use netcdfio

       implicit none

       !==============INPUTS and OUTPUTS===============
       integer, intent(in) :: fid
       character(len=4), intent(in) :: vn
       real, dimension(npts), intent(out) :: ro

       !==============LOCAL VARIABLES==================
       integer :: vid
       character(len=50) :: dn   ! Dummy Name
       integer*1,allocatable :: byte_array(:)
       integer*2,allocatable :: short_array(:)
       integer :: o

       !===============================================
       
       !write(*,*) ' Getting variable ',trim(vn),' from file ',fid

       vid = ncvid(fid,trim(vn),exitcode)

       call ncvinq(fid,vid,dn,vtype,nvdims,vdims,nvatts,exitcode)

       if ( nvatts .eq. 0 ) then

          if ( vtype .eq. ncfloat ) then
             ! Data are real numbers, unpacked
             ! so read directly
             call ncvgt(fid,vid,lag_readst2d,lag_readct2d,ro,exitcode)

          elseif ( vtype .eq. ncbyte ) then
             ! Data are bytes, month and day
             allocate(byte_array(npts))
             call ncvgt(fid,vid,lag_readst2d,lag_readct2d,byte_array,exitcode)
             do o = 1,npts
                ro(o) = real(byte_array(o))
             enddo
             deallocate(byte_array)

          elseif ( vtype .eq. ncshort ) then
             ! Data are unpacked integers, convert to real
             ! This is reserved for iit,jjt,kkt,kbot.
             allocate(short_array(npts))
             call ncvgt(fid,vid,lag_readst2d,lag_readct2d,short_array,exitcode)
             do o = 1,npts
                ro(o) = real(short_array(o))
             enddo
             deallocate(short_array)
             
          else
             write(*,*) ' ERROR: Unrecognized variable format!'
          endif

       elseif ( (nvatts .eq. 1) .and. (vtype .eq. ncbyte) ) then
          ! Data have only an add_offset,
          ! must be years stored in bytes
          call ncagt(fid,vid,'add_offset',add_offset,exitcode)
          allocate(byte_array(npts))
          call ncvgt(fid,vid,lag_readst2d,lag_readct2d,byte_array,exitcode)
          do o = 1,npts
             ro(o) = real(byte_array(o)) + add_offset(1)
          enddo
          deallocate(byte_array)

       elseif ( (nvatts .eq. 2) .and. (vtype .eq. ncshort) ) then
          ! Data have both a scale_factor and an add_offset
          ! and are short integers.
          call ncagt(fid,vid,'add_offset',add_offset,exitcode)
          call ncagt(fid,vid,'scale_factor',scale_factor,exitcode)
          allocate(short_array(npts))
          call ncvgt(fid,vid,lag_readst2d,lag_readct2d,short_array,exitcode)
          do o = 1,npts
             ro(o) = real(short_array(o))*scale_factor(1) + add_offset(1)
          enddo
          deallocate(short_array)

       else
          write(*,*) ' ERROR: Incorrect combination of var atts and var type!'
          stop
       endif

       return

     end subroutine get_lag_var

!------------------------------------------------

!------------------------------------------------

     subroutine get_lag_var_seg(fid,vn,ro,ps,pc)

!------------------------------------------------
! This subroutine will load a Lagrangian variable
! from a tcdf formatted file into a the real
! variable ro given,
! fid   = file ID
! vn    = variable name
! ro    = real output variable
! ps    = the start index in traj (like readstart)
! pc    = the count index in traj (like readcount)
!
! This is different from get_lag_var because only
! a time-delimited segment within the whole traj.
! is pulled out.  It is assumed that ro has the
! correct length corresponding to pc.
!
! This subroutine will automatically search for
! the type of variable in the source file and
! convert that source variable (including unpack)
! to a real number array.
!
! Quite a bit of additional information is being
! passed through the lag_readst2d, lag_readct2d, npts,
! variables which are in tcdfio.
! This is done to make call lines simple and
! keep things that are constant out of the way.
!------------------------------------------------

       use tcdfio
       use netcdfio

       implicit none

       !==============INPUTS and OUTPUTS===============
       integer, intent(in) :: fid
       character(len=4), intent(in) :: vn
       integer, intent(in) :: ps, pc
       real, dimension(pc), intent(out) :: ro

       !==============LOCAL VARIABLES==================
       integer :: vid
       character(len=50) :: dn   ! Dummy Name
       integer*1,allocatable :: byte_array(:)
       integer*2,allocatable :: short_array(:)
       integer :: o
       integer :: local_readst(2)
       integer :: local_readct(2)
       !===============================================
       
       !write(*,*) ' Getting variable ',trim(vn),' from file ',fid

       vid = ncvid(fid,trim(vn),exitcode)

       call ncvinq(fid,vid,dn,vtype,nvdims,vdims,nvatts,exitcode)

       local_readst(1) = ps
       local_readst(2) = lag_readst2d(2)

       local_readct(1) = pc
       local_readct(2) = lag_readct2d(2)

       if ( nvatts .eq. 0 ) then

          if ( vtype .eq. ncfloat ) then
             ! Data are real numbers, unpacked
             ! so read directly
             call ncvgt(fid,vid,local_readst,local_readct,ro,exitcode)

          elseif ( vtype .eq. ncbyte ) then
             ! Data are bytes, month and day
             allocate(byte_array(pc))
             call ncvgt(fid,vid,local_readst,local_readct,byte_array,exitcode)
             do o = 1,pc
                ro(o) = real(byte_array(o))
             enddo
             deallocate(byte_array)

          elseif ( vtype .eq. ncshort ) then
             ! Data are iit,jjt,kkt,kbot
             allocate(short_array(pc))
             call ncvgt(fid,vid,local_readst,local_readct,short_array,exitcode)
             do o = 1,pc
                ro(o) = real(short_array(o))
             enddo
             deallocate(short_array)

          else
             write(*,*) ' ERROR: Unrecognized variable format!'
          endif

       elseif ( (nvatts .eq. 1) .and. (vtype .eq. ncbyte) ) then
          ! Data have only an add_offset,
          ! must be years stored in bytes
          call ncagt(fid,vid,'add_offset',add_offset,exitcode)
          allocate(byte_array(pc))
          call ncvgt(fid,vid,local_readst,local_readct,byte_array,exitcode)
          do o = 1,pc
             ro(o) = real(byte_array(o)) + add_offset(1)
          enddo
          deallocate(byte_array)

       elseif ( (nvatts .eq. 2) .and. (vtype .eq. ncshort) ) then
          ! Data have both a scale_factor and an add_offset
          ! and are short integers.
          call ncagt(fid,vid,'add_offset',add_offset,exitcode)
          call ncagt(fid,vid,'scale_factor',scale_factor,exitcode)
          allocate(short_array(pc))
          call ncvgt(fid,vid,local_readst,local_readct,short_array,exitcode)
          do o = 1,pc
             ro(o) = real(short_array(o))*scale_factor(1) + add_offset(1)
          enddo
          deallocate(short_array)

       else
          write(*,*) ' ERROR: Incorrect combination of var atts and var type!'
          stop
       endif

       return

     end subroutine get_lag_var_seg

!------------------------------------------------

!------------------------------------------------

     subroutine put_lag_var(fid,vn,ri,np)

!------------------------------------------------
! This subroutine will write a Lagrangian variable
! from real variable into a file.  It may be packed
! or unpacked, depending on the format of the
! variable in the original file.  The variable being
! put must first be created by dup_lag_var.
!
! fid   = file ID
! vn    = variable name
! ri    = real input variable
! np    = number of time steps in lag var
!
! This subroutine will automatically search for
! the type of variable to be written and
! convert the input real array into the correct
! format (including packing).
!
! Quite a bit of additional information is being
! passed through the lag_writest2d and lag_writect2d
! variables which are in tcdfio.  This is done to
! make call lines simple and keep things that
! are constant out of the way.
!------------------------------------------------

       use tcdfio
       use netcdfio

       implicit none

       !==============INPUTS and OUTPUTS===============
       integer, intent(in) :: fid
       real, intent(in) :: ri(*)
       character(len=4), intent(in) :: vn
       integer, intent(in) :: np

       !==============LOCAL VARIABLES==================
       integer :: vid
       character(len=50) :: dn   ! Dummy Name
       integer*1,allocatable :: byte_array(:)
       integer*2,allocatable :: short_array(:)
       integer :: o
       !===============================================

       vid = ncvid(fid,trim(vn),exitcode)

       call ncvinq(fid,vid,dn,vtype,nvdims,vdims,nvatts,exitcode)

       if ( nvatts .eq. 0 ) then

          if ( vtype .eq. ncfloat ) then
             ! Data are real numbers, unpacked
             ! so write directly
             call ncvpt(fid,vid,lag_writest2d,lag_writect2d,ri,exitcode)

          elseif ( vtype .eq. ncbyte ) then
             ! Data are bytes, month and day
             allocate(byte_array(np))
             do o = 1,np
                byte_array(o) = int(ri(o),1)
             enddo
             call ncvpt(fid,vid,lag_writest2d,lag_writect2d,byte_array,exitcode)
             deallocate(byte_array)

          elseif ( vtype .eq. ncshort ) then
             ! Data are iit,jjt,kkt,kbot
             allocate(short_array(np))
             do o = 1,np
                ! Note that for index values we use nint to
                ! round to the nearest integer because that
                ! will place particles within their respective
                ! grid boxes.  In the other data, e.g. x,y,z,t,s
                ! we are truncating extra digits after dividing
                ! by the scale factor, so int is more appropriate.
                short_array(o) = nint(ri(o),2)
             enddo
             call ncvpt(fid,vid,lag_writest2d,lag_writect2d,short_array,exitcode)
             deallocate(short_array)

          else
             write(*,*) ' ERROR: Unrecognized variable format!'
          endif

       elseif ( (nvatts .eq. 1) .and. (vtype .eq. ncbyte) ) then
          ! Data have only an add_offset,
          ! must be years stored in bytes
          call ncagt(fid,vid,'add_offset',add_offset,exitcode)
          allocate(byte_array(np))
          do o = 1,np
             byte_array(o) = int(ri(o) - add_offset(1),1)
          enddo
          call ncvpt(fid,vid,lag_writest2d,lag_writect2d,byte_array,exitcode)
          deallocate(byte_array)

       elseif ( (nvatts .eq. 2) .and. (vtype .eq. ncshort) ) then
          ! Data have both a scale_factor and an add_offset
          ! and are short integers.
          call ncagt(fid,vid,'add_offset',add_offset,exitcode)
          call ncagt(fid,vid,'scale_factor',scale_factor,exitcode)
          allocate(short_array(np))
          do o = 1,np
             short_array(o) = int((ri(o) - add_offset(1))/scale_factor(1),2)
          enddo
          call ncvpt(fid,vid,lag_writest2d,lag_writect2d,short_array,exitcode)
          deallocate(short_array)

       else
          write(*,*) ' ERROR: Incorrect combination of var atts and var type!'
          stop
       endif

       return

     end subroutine put_lag_var

!------------------------------------------------

!------------------------------------------------

     subroutine put_lag_var_seg2(fid,vn,ri,ps,pc)

!------------------------------------------------
! This subroutine will write a Lagrangian variable
! from real variable into a file.  It may be packed
! or unpacked, depending on the format of the
! variable in the original file.  The variable being
! put must first be created by dup_lag_var.
!
! fid   = file ID
! vn    = variable name
! ri    = real input variable
! ps    = start point of writing (like readstart)
! pc    = number of time steps to put (like readcount)
!
! This is different from put_lag_var_seg because
! we use readstart and readcount style notation
! here while put_lag_var_seg uses other notation.
! This is ultimately the most general version.
! 
! This subroutine will automatically search for
! the type of variable to be written and
! convert the input real array into the correct
! format (including packing).
!
! Quite a bit of additional information is being
! passed through the lag_writest2d and lag_writect2d
! variables which are in tcdfio.  This is done to
! make call lines simple and keep things that
! are constant out of the way.
!------------------------------------------------

       use tcdfio
       use netcdfio

       implicit none

       !==============INPUTS and OUTPUTS===============
       integer, intent(in) :: fid
       integer, intent(in) :: ps
       integer, intent(in) :: pc
       real, intent(in) :: ri(pc)
       character(len=4), intent(in) :: vn

       !==============LOCAL VARIABLES==================
       integer :: vid
       character(len=50) :: dn   ! Dummy Name
       integer*1,allocatable :: byte_array(:)
       integer*2,allocatable :: short_array(:)
       integer :: o
       integer :: local_writest(2)
       integer :: local_writect(2)
       !===============================================

       vid = ncvid(fid,trim(vn),exitcode)

       call ncvinq(fid,vid,dn,vtype,nvdims,vdims,nvatts,exitcode)

       local_writest(1) = ps
       local_writest(2) = lag_writest2d(2)

       local_writect(1) = pc
       local_writect(2) = lag_writect2d(2)

       if ( nvatts .eq. 0 ) then

          if ( vtype .eq. ncfloat ) then
             ! Data are real numbers, unpacked
             ! so write directly
             call ncvpt(fid,vid,local_writest,local_writect,ri,exitcode)

          elseif ( vtype .eq. ncbyte ) then
             ! Data are bytes, month and day
             allocate(byte_array(pc))
             do o = 1,pc
                byte_array(o) = int(ri(o),1)
             enddo
             call ncvpt(fid,vid,local_writest,local_writect,byte_array,exitcode)
             deallocate(byte_array)

          elseif ( vtype .eq. ncshort ) then
             ! Data are iit,jjt,kkt,kbot
             allocate(short_array(pc))
             do o = 1,pc
                short_array(o) = int(ri(o),2)
             enddo
             call ncvpt(fid,vid,local_writest,local_writect,short_array,exitcode)
             deallocate(short_array)

          else
             write(*,*) ' ERROR: Unrecognized variable format!'
          endif

       elseif ( (nvatts .eq. 1) .and. (vtype .eq. ncbyte) ) then
          ! Data have only an add_offset,
          ! must be years stored in bytes
          call ncagt(fid,vid,'add_offset',add_offset,exitcode)
          allocate(byte_array(pc))
          do o = 1,pc
             byte_array(o) = int(ri(o) - add_offset(1),1)
          enddo
          call ncvpt(fid,vid,local_writest,local_writect,byte_array,exitcode)
          deallocate(byte_array)

       elseif ( (nvatts .eq. 2) .and. (vtype .eq. ncshort) ) then
          ! Data have both a scale_factor and an add_offset
          ! and are short integers.
          call ncagt(fid,vid,'add_offset',add_offset,exitcode)
          call ncagt(fid,vid,'scale_factor',scale_factor,exitcode)
          allocate(short_array(pc))
          do o = 1,pc
             short_array(o) = int((ri(o) - add_offset(1))/scale_factor(1),2)
          enddo
          call ncvpt(fid,vid,local_writest,local_writect,short_array,exitcode)
          deallocate(short_array)

       else
          write(*,*) ' ERROR: Incorrect combination of var atts and var type!'
          stop
       endif

       return

     end subroutine put_lag_var_seg2

!------------------------------------------------

!------------------------------------------------

     subroutine dup_lag_var(fid1,fid2,vn)

!------------------------------------------------
! This subroutine will duplicate the root of
! the var in fid1 to fid2 without writing values
! into the data storage.
!
! fid1   = source file ID
! fid2   = destination file ID
! vn     = variable name
!
! This subroutine will automatically search for
! the type of variable to be written and
! convert the input real array into the correct
! format (including packing).
!------------------------------------------------

       use tcdfio
       use netcdfio

       implicit none

       !==============INPUTS and OUTPUTS===============
       integer, intent(in) :: fid1, fid2
       character(len=4), intent(in) :: vn

       !==============LOCAL VARIABLES==================
       integer :: vid1, vid2
       character(len=50) :: dn   ! Dummy Name

       !===============================================

       ! Get variable ID and info for source variable
       vid1 = ncvid(fid1,trim(vn),exitcode)
       call ncvinq(fid1,vid1,dn,vtype,nvdims,vdims,nvatts,exitcode)

       ! Create variable in other file
       vid2 = ncvdef(fid2,trim(vn),vtype,nvdims,vdims_out,exitcode)

       ! Copy over attributes if present
       if ( nvatts .eq. 0 ) then
          ! No attributes to copy over.
       elseif ( nvatts .eq. 1 ) then
          ! A single attribute means only an add_offset
          call ncagt(fid1,vid1,'add_offset',add_offset,exitcode)
          call ncapt(fid2,vid2,'add_offset',ncfloat,1,add_offset,exitcode)

       elseif ( nvatts .eq. 2 ) then
          ! Copy both an add_offset and scale_factor
          call ncagt(fid1,vid1,'add_offset',add_offset,exitcode)
          call ncapt(fid2,vid2,'add_offset',ncfloat,1,add_offset,exitcode)

          call ncagt(fid1,vid1,'scale_factor',scale_factor,exitcode)
          call ncapt(fid2,vid2,'scale_factor',ncfloat,1,scale_factor,exitcode)

       else
          write(*,*) ' ERROR: Incorrect number of var atts!'
          stop
       endif

       return

     end subroutine dup_lag_var

!------------------------------------------------

!------------------------------------------------

     subroutine put_lag_var_seg(fid,vn,ri,np,lo,hi,mask)

!------------------------------------------------
! This subroutine will put the segment (defined
! by index locations lo,hi) of the Lagrangian
! variable vn (values in ri) into the file fid.
! The length of ri is specified by np.
!
! Variable masking happens here as well, given
! by the mask value.
!------------------------------------------------

       use tcdfio
       use netcdfio

       implicit none

       !==============INPUTS and OUTPUTS===============
       integer, intent(in) :: fid
       integer, intent(in) :: np
       real, intent(in) :: ri(np)
       character(len=4), intent(in) :: vn
       integer, intent(in) :: lo
       integer, intent(in) :: hi
       real, intent(in) :: mask

       !==============LOCAL VARIABLES==================
       integer :: vid
       character(len=50) :: dn   ! Dummy Name
       integer*1,allocatable :: byte_array(:)
       integer*2,allocatable :: short_array(:)
       real, allocatable :: float_array(:)
       integer :: o, oo, len_seg, len_traj
       !===============================================

       ! Get variable information

       vid = ncvid(fid,trim(vn),exitcode)

       call ncvinq(fid,vid,dn,vtype,nvdims,vdims,nvatts,exitcode)

       !==============================================

       ! Initial error checking
       ! We allow lo to be less than 1
       ! hi to be greater than size(ri)
       ! because those overlapping values
       ! will be masked out.  What does not
       ! make any sense is for lo to be
       ! beyond size(ri) and for hi to be
       ! less than 1.

       if ( lo .gt. hi ) then
          write(*,*) 'ERROR in put_lag_var_seg: lo > hi'
          write(*,*) 'hi = ',hi
          write(*,*) 'lo = ',lo
          stop
       endif

       if ( hi .lt. 1 ) then
          write(*,*) 'ERROR in put_lag_var_seg: hi < 1'
          write(*,*) 'hi = ',hi
          stop
       endif

       len_traj = size(ri)
       if ( lo .gt. len_traj ) then
          write(*,*) 'ERROR in put_lag_var_seg: lo > len'
          write(*,*) 'lo = ',lo
          write(*,*) 'len= ',len_traj
          stop
       endif

       ! Get length of subsegment
       len_seg = hi - lo + 1

       ! Check valid length
       if ( len_seg .lt. 1 ) then
          write(*,*) 'ERROR in put_lag_var_seg: segment length < 1'
          write(*,*) 'hi = ',hi
          write(*,*) 'lo = ',lo
          write(*,*) 'len= ',len_seg
          stop
       endif

       !==============================================

       ! Write to file for each case

       if ( nvatts .eq. 0 ) then

          if ( vtype .eq. ncfloat ) then
             ! Data are real numbers, unpacked
             allocate(float_array(len_seg))
             do o = 1,len_seg
                oo = lo + o - 1
                ! For o = 1, oo = lo
                ! For o = len, oo = lo + len - 1 = hi
                if ( (oo .lt. 1) .or. (oo .gt. len_traj) ) then
                   float_array(o) = mask
                else
                   float_array(o) = ri(oo)
                endif
             enddo
             call ncvpt(fid,vid,lag_writest2d,lag_writect2d,float_array,exitcode)
             deallocate(float_array)

          elseif ( vtype .eq. ncbyte ) then
             ! Data are bytes, month and day, unpacked
             allocate(byte_array(len_seg))
             do o = 1,len_seg
                oo = lo + o - 1
                if ( (oo .lt. 1) .or. (oo .gt. len_traj) ) then
                   byte_array(o) = int(mask,1)
                else
                   byte_array(o) = int(ri(oo),1)
                endif
             enddo
             call ncvpt(fid,vid,lag_writest2d,lag_writect2d,byte_array,exitcode)
             deallocate(byte_array)

          elseif ( vtype .eq. ncshort ) then
             ! Data are iit,jjt,kkt,kbot
             allocate(short_array(len_seg))
             do o = 1,len_seg
                oo = lo + o - 1
                if ( (oo .lt. 1) .or. (oo .gt. len_traj) ) then
                   short_array(o) = int(mask,2)
                else
                   short_array(o) = int(ri(oo),2)
                endif
             enddo
             call ncvpt(fid,vid,lag_writest2d,lag_writect2d,short_array,exitcode)
             deallocate(short_array)
             
          else
             write(*,*) ' ERROR: Unrecognized variable format!'
          endif

       elseif ( (nvatts .eq. 1) .and. (vtype .eq. ncbyte) ) then
          ! Data have only an add_offset,
          ! must be years stored in bytes
          call ncagt(fid,vid,'add_offset',add_offset,exitcode)
          allocate(byte_array(len_seg))
          do o = 1,len_seg
             oo = lo + o - 1
             if ( (oo .lt. 1) .or. (oo .gt. len_traj) ) then
                byte_array(o) = int(mask - add_offset(1),1)
             else
                byte_array(o) = int(ri(oo) - add_offset(1),1)
             endif
          enddo
          call ncvpt(fid,vid,lag_writest2d,lag_writect2d,byte_array,exitcode)
          deallocate(byte_array)

       elseif ( (nvatts .eq. 2) .and. (vtype .eq. ncshort) ) then
          ! Data have both a scale_factor and an add_offset
          ! and are short integers.
          call ncagt(fid,vid,'add_offset',add_offset,exitcode)
          call ncagt(fid,vid,'scale_factor',scale_factor,exitcode)
          allocate(short_array(len_seg))
          do o = 1,len_seg
             oo = lo + o - 1
             if ( (oo .lt. 1) .or. (oo .gt. len_traj) ) then
                short_array(o) = int((mask - add_offset(1))/scale_factor(1),2)
             else
                short_array(o) = int((ri(oo) - add_offset(1))/scale_factor(1),2)
             endif
          enddo
          call ncvpt(fid,vid,lag_writest2d,lag_writect2d,short_array,exitcode)
          deallocate(short_array)

       else
          write(*,*) ' ERROR: Incorrect combination of var atts and var type!'
          stop
       endif

       return

     end subroutine put_lag_var_seg

!------------------------------------------------

!------------------------------------------------

     subroutine put_lag_var_seg_pad(fid,vn,ri,np,lo,hi,mask)

!------------------------------------------------
! This subroutine will put the segment (defined
! by index locations lo,hi) of the Lagrangian
! variable vn (values in ri) into the file fid.
! The length of ri is specified by np.
!
! Variable masking happens here as well, given
! by the mask value.
!
! The _pad variation on put_lag_var_seg is
! identical to the original version except that
! segments can have unequal lengths (i.e. for
! lines) so the unused portion of the traj is
! padded with the given mask values at the end.
!
! Therefore, put_lag_var_seg is best used with
! worms (because all worms have the same length)
! and put_lag_var_seg_pad is best used with
! lines (because lines have unequal lengths).
!------------------------------------------------

       use tcdfio
       use netcdfio

       implicit none

       !==============INPUTS and OUTPUTS===============
       integer, intent(in) :: fid
       integer, intent(in) :: np
       real, intent(in) :: ri(np)
       character(len=4), intent(in) :: vn
       integer, intent(in) :: lo
       integer, intent(in) :: hi
       real, intent(in) :: mask

       !==============LOCAL VARIABLES==================
       integer :: vid
       character(len=50) :: dn   ! Dummy Name
       integer*1,allocatable :: byte_array(:)
       integer*2,allocatable :: short_array(:)
       real, allocatable :: float_array(:)
       integer :: o, oo, len_seg, len_traj
       !===============================================

       ! Get variable information

       vid = ncvid(fid,trim(vn),exitcode)

       call ncvinq(fid,vid,dn,vtype,nvdims,vdims,nvatts,exitcode)

       !==============================================

       ! Initial error checking
       ! We allow lo to be less than 1
       ! hi to be greater than size(ri)
       ! because those overlapping values
       ! will be masked out.  What does not
       ! make any sense is for lo to be
       ! beyond size(ri) and for hi to be
       ! less than 1.

       if ( lo .gt. hi ) then
          write(*,*) 'ERROR in put_lag_var_seg: lo > hi'
          write(*,*) 'hi = ',hi
          write(*,*) 'lo = ',lo
          stop
       endif

       if ( hi .lt. 1 ) then
          write(*,*) 'ERROR in put_lag_var_seg: hi < 1'
          write(*,*) 'hi = ',hi
          stop
       endif

       len_traj = size(ri)
       if ( lo .gt. len_traj ) then
          write(*,*) 'ERROR in put_lag_var_seg: lo > len'
          write(*,*) 'lo = ',lo
          write(*,*) 'len= ',len_traj
          stop
       endif

       ! Get length of subsegment
       len_seg = hi - lo + 1

       ! Check valid length
       if ( len_seg .lt. 1 ) then
          write(*,*) 'ERROR in put_lag_var_seg: segment length < 1'
          write(*,*) 'hi = ',hi
          write(*,*) 'lo = ',lo
          write(*,*) 'len= ',len_seg
          stop
       endif

       !==============================================

       ! Write to file for each case

       if ( nvatts .eq. 0 ) then

          if ( vtype .eq. ncfloat ) then
             ! Data are real numbers, unpacked
             allocate(float_array(lag_writect2d(1)))
             float_array = mask
             do o = 1,len_seg
                oo = lo + o - 1
                ! For o = 1, oo = lo
                ! For o = len, oo = lo + len - 1 = hi
                if ( (oo .lt. 1) .or. (oo .gt. len_traj) ) then
                   float_array(o) = mask
                else
                   float_array(o) = ri(oo)
                endif
             enddo
             call ncvpt(fid,vid,lag_writest2d,lag_writect2d,float_array,exitcode)
             deallocate(float_array)

          elseif ( vtype .eq. ncbyte ) then
             ! Data are bytes, month and day, unpacked
             allocate(byte_array(lag_writect2d(1)))
             byte_array = int(mask,1)
             do o = 1,len_seg
                oo = lo + o - 1
                if ( (oo .lt. 1) .or. (oo .gt. len_traj) ) then
                   byte_array(o) = int(mask,1)
                else
                   byte_array(o) = int(ri(oo),1)
                endif
             enddo
             call ncvpt(fid,vid,lag_writest2d,lag_writect2d,byte_array,exitcode)
             deallocate(byte_array)

          elseif ( vtype .eq. ncshort ) then
             ! Data are iit,jjt,kkt,kbot
             allocate(short_array(lag_writect2d(1)))
             short_array = int(mask,2)
             do o = 1,len_seg
                oo = lo + o - 1
                if ( (oo .lt. 1) .or. (oo .gt. len_traj) ) then
                   short_array(o) = int(mask,2)
                else
                   short_array(o) = int(ri(oo),2)
                endif
             enddo
             call ncvpt(fid,vid,lag_writest2d,lag_writect2d,short_array,exitcode)
             deallocate(short_array)

          else
             write(*,*) ' ERROR: Unrecognized variable format!'
          endif

       elseif ( (nvatts .eq. 1) .and. (vtype .eq. ncbyte) ) then
          ! Data have only an add_offset,
          ! must be years stored in bytes
          call ncagt(fid,vid,'add_offset',add_offset,exitcode)
          allocate(byte_array(lag_writect2d(1)))
          byte_array = int(mask - add_offset(1),1)
          do o = 1,len_seg
             oo = lo + o - 1
             if ( (oo .lt. 1) .or. (oo .gt. len_traj) ) then
                byte_array(o) = int(mask - add_offset(1),1)
             else
                byte_array(o) = int(ri(oo) - add_offset(1),1)
             endif
          enddo
          call ncvpt(fid,vid,lag_writest2d,lag_writect2d,byte_array,exitcode)
          deallocate(byte_array)

       elseif ( (nvatts .eq. 2) .and. (vtype .eq. ncshort) ) then
          ! Data have both a scale_factor and an add_offset
          ! and are short integers.
          call ncagt(fid,vid,'add_offset',add_offset,exitcode)
          call ncagt(fid,vid,'scale_factor',scale_factor,exitcode)
          allocate(short_array(lag_writect2d(1)))
          short_array = int((mask - add_offset(1))/scale_factor(1),2)
          do o = 1,len_seg
             oo = lo + o - 1
             if ( (oo .lt. 1) .or. (oo .gt. len_traj) ) then
                short_array(o) = int((mask - add_offset(1))/scale_factor(1),2)
             else
                short_array(o) = int((ri(oo) - add_offset(1))/scale_factor(1),2)
             endif
          enddo
          call ncvpt(fid,vid,lag_writest2d,lag_writect2d,short_array,exitcode)
          deallocate(short_array)

       else
          write(*,*) ' ERROR: Incorrect combination of var atts and var type!'
          stop
       endif

       return

     end subroutine put_lag_var_seg_pad

!------------------------------------------------

!------------------------------------------------
     subroutine cre_lag_var(fid1,vn,lp)
!------------------------------------------------
! This subroutine creates a new lagrangian variable
! in the output file.  Very similar to dup_lag_var
! but useful for creating new variables when you
! don't have an input file as a template.
!
! All you need is a file id for an open netcdf
! tcdf file with dimensions already created and
! vdims_out already assigned.
!
! Also added the lp logical to set whether packing
! is activated or not.
!------------------------------------------------
       use tcdfio
       use netcdfio

       implicit none

       !==============INPUTS and OUTPUTS===============
       integer, intent(in) :: fid1
       character(len=4), intent(in) :: vn
       logical, intent(in) :: lp

       ! Local
       integer :: vid

       ! Create variable and set packing attributes
       if( (index(trim(vn),trim(lamvnam)) .ne. 0) .or. &
           (index(trim(vn),trim(phivnam)) .ne. 0) ) then
          ! Longitude or Latitude in x,y
          if ( lp ) then
             vid = ncvdef(fid1,vn,ncshort,2,vdims_out,exitcode)

             ! Write packing attributes
             call ncapt(fid1,vid,'scale_factor',ncfloat,1,xy_scale(1),exitcode)
             call ncapt(fid1,vid,'add_offset',ncfloat,1,xy_offset(1),exitcode)
          else
             vid = ncvdef(fid1,vn,ncfloat,2,vdims_out,exitcode)
          endif

       elseif ( (index(trim(vn),trim(iitvnam)) .ne. 0) .or. &
                (index(trim(vn),trim(jjtvnam)) .ne. 0) .or. &
                (index(trim(vn),trim(kktvnam)) .ne. 0) .or. &
                (index(trim(vn),trim(kbotvnam)) .ne. 0) ) then
          ! Longitude or latitude or depth in i,j,k
          if ( lp ) then
             vid = ncvdef(fid1,vn,ncshort,2,vdims_out,exitcode)
             ! No other packing attributes needed.
          else
             vid = ncvdef(fid1,vn,ncfloat,2,vdims_out,exitcode)
          endif
          
       elseif ( (index(trim(vn),trim(depvnam)) .ne. 0) .or. &
                (index(trim(vn),trim(bdepvnam)) .ne. 0) ) then
          ! Depth or bottom depth
          if ( lp ) then
             vid = ncvdef(fid1,vn,ncshort,2,vdims_out,exitcode)

             ! Write packing attributes
             call ncapt(fid1,vid,'scale_factor',ncfloat,1,z_scale(1),exitcode)
             call ncapt(fid1,vid,'add_offset',ncfloat,1,z_offset(1),exitcode)
          else
             vid = ncvdef(fid1,vn,ncfloat,2,vdims_out,exitcode)
          endif

       elseif (index(trim(vn),trim(tempvnam)) .ne. 0) then
          ! Temperature
          if ( lp ) then
             vid = ncvdef(fid1,vn,ncshort,2,vdims_out,exitcode)

             ! Write packing attributes
             call ncapt(fid1,vid,'scale_factor',ncfloat,1,ts_scale(1),exitcode)
             call ncapt(fid1,vid,'add_offset',ncfloat,1,zero_offset(1),exitcode)
          else
             vid = ncvdef(fid1,vn,ncfloat,2,vdims_out,exitcode)
          endif

       elseif (index(trim(vn),trim(saltvnam)) .ne. 0) then
          ! Salinity
          if ( lp ) then
             vid = ncvdef(fid1,vn,ncshort,2,vdims_out,exitcode)

             ! Write packing attributes
             call ncapt(fid1,vid,'scale_factor',ncfloat,1,ts_scale(1),exitcode)
             call ncapt(fid1,vid,'add_offset',ncfloat,1,s_offset(1),exitcode)
          else
             vid = ncvdef(fid1,vn,ncfloat,2,vdims_out,exitcode)
          endif

       else
          ! Missing Rho, A, B, layer, Q
          write(*,*) 'ERROR in tcdfio::cre_lag_var'
          write(*,*) 'Not yet setup to deal with variable:'
          write(*,*) vn
          return
       endif
       return

     end subroutine cre_lag_var
!------------------------------------------------
   end module load_tcdf
