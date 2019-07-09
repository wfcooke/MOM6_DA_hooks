!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This module contains a set of dummy interfaces for compiling the MOM6 DA
! driver code. These interfaces are not finalized and will be replaced by supported
! interfaces at some later date.
!
! 3/22/18
! matthew.harrison@noaa.gov
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module ocean_da_core_mod


  use mpp_domains_mod, only : domain2d
  use fms_mod, only : file_exist, read_data
  use fms_mod, only : open_namelist_file, check_nml_error, close_file
  use fms_mod, only : error_mesg, FATAL, NOTE
#ifdef INTERNAL_FILE_NML
  USE mpp_mod, ONLY: input_nml_file
#endif
  use mpp_mod, only : mpp_sum, stdout, stdlog, mpp_sync_self
  use mpp_mod, only : mpp_pe, mpp_root_pe
  use mpp_io_mod, only : mpp_open, mpp_close, MPP_ASCII, MPP_RDONLY, MPP_MULTI, MPP_SINGLE, MPP_NETCDF
  use mpp_io_mod, only : mpp_get_atts, mpp_get_info, mpp_get_fields, mpp_read, axistype, fieldtype, mpp_get_axes
  use mpp_io_mod, only : mpp_get_times
  use mpp_io_mod, only : mpp_get_axis_data, mpp_get_field_name
  use mpp_domains_mod, only : mpp_get_compute_domain, mpp_get_data_domain
  use mpp_domains_mod, only : domain2d, mpp_get_global_domain
  use time_manager_mod, only : time_type, set_time, get_date
  use time_manager_mod, only : decrement_time, increment_time
  use time_manager_mod, only : operator( <= ), operator( >= ), operator(*)
  use time_manager_mod, only : operator( - ), operator( > ), operator ( < )
  use get_cal_time_mod, only : get_cal_time
  use axis_utils_mod, only : frac_index
  use horiz_interp_type_mod, only: horiz_interp_type
  use horiz_interp_bilinear_mod, only : horiz_interp_bilinear_new
  use constants_mod, only : DEG_TO_RAD
  ! ODA_tools modules
  use ocean_da_types_mod, only : ocean_profile_type, grid_type
  use ocean_da_types_mod, only : forward_operator_type
  use ocean_da_types_mod, only : TEMP_ID, SALT_ID, MISSING_VALUE
  use ocean_da_types_mod, only : ODA_PFL, ODA_XBT, ODA_MRB, ODA_OISST
  use ocean_da_types_mod, only : UNKNOWN, MAX_LEVELS_FILE, MAX_LINKS
  use kdtree, only : kd_root, kd_search_nnearest, kd_init
  use loc_and_dist_mod, only : within_domain

  implicit none
  private
  public :: ocean_da_core_init
  public :: get_profiles

  ! Parameters
  integer, parameter :: PROFILE_FILE = 1
  integer, parameter :: SURFACE_FILE = 2
  integer, parameter :: ARGO_FILE = 3
! time window for DROP, MOORING and SATELLITE data respectively
  type(time_type) , dimension(10) :: time_window
  integer, dimension(10) :: type_count = 0

  integer, allocatable, dimension(:) :: lon1d, lat1d
  real, allocatable, dimension(:) :: glon1d, glat1d
  type(kd_root), allocatable :: kdroot

  ! ocean_obs_nml variables
  integer :: max_levels = 1000 !< maximium number of levels for a single profile
  real, dimension(10) :: obs_sbound = -89.0 !< set obs domain
  real, dimension(10) :: obs_nbound = 89.0 !< set obs domain
  real :: depth_cut = 2000.0
  real, dimension(10) :: data_window = 24.0
  integer, dimension(10) :: sec_offset = 0
  integer, dimension(10) :: day_offset = 0
  real, dimension(10) :: temp_error = 1.0
  real, dimension(10) :: salt_error = 0.2
  integer, dimension(10) :: impact_levels = 3
  real :: sst_vimpact_temp = -1.75
  integer :: sst_vimpact_levels = 5
  real, dimension(10) :: temp_dist = 200.0e3
  real, dimension(10) :: salt_dist = 200.0e3
  logical, dimension(10) :: temp_to_salt = .false.
  logical, dimension(10) :: salt_to_temp = .false.
  integer :: obs_days_plus, obs_days_minus
  logical :: temp_obs, salt_obs
  integer :: max_files = 30
  real :: shelf_depth = 500.0
  namelist /ocean_obs_nml/ max_levels, obs_sbound, obs_nbound, depth_cut, &
          data_window, sec_offset, day_offset, temp_error, salt_error, impact_levels, &
          temp_dist, salt_dist, temp_to_salt, salt_to_temp, &
          sst_vimpact_temp, sst_vimpact_levels, &
          temp_obs, salt_obs, max_files, obs_days_minus, obs_days_plus

contains

  subroutine ocean_da_core_init(Domain, T_grid, Profiles, model_time)
    type(domain2d), pointer, intent(in) :: Domain
    type(grid_type), pointer, intent(in) :: T_grid
    type(ocean_profile_type), pointer :: Profiles
    type(time_type), intent(in) :: model_time

    type obs_entry_type
       character(len=128) :: filename
       character(len=16)  :: file_type
    end type obs_entry_type

    type(time_type)  :: time_s, time_e
    integer :: data_seconds
    integer :: i, j, n, obs_variable, ni, nj
    integer :: ioun, io_status, ierr
    integer :: stdout_unit, stdlog_unit
    integer :: nfiles, nrecs, unit
    integer, dimension(:), allocatable :: filetype

    character(len=256) :: record
    character(len=128), dimension(:), allocatable :: input_files

    type(obs_entry_type) :: tbl_entry

    stdout_unit = stdout()
    stdlog_unit = stdlog()

    ioun = open_namelist_file()
    read(UNIT=ioun, NML=ocean_obs_nml, IOSTAT=io_status)
    ierr = check_nml_error(io_status,'ocean_obs_nml')
    call close_file(ioun)

    write (UNIT=stdlog_unit, NML=ocean_obs_nml)
    if ( mpp_pe() == mpp_root_pe() ) then
       write (UNIT=stdout_unit, NML=ocean_obs_nml)
    end if

    ! Allocate filetype* and input_files* variables
    allocate(filetype(max_files), input_files(max_files))

    filetype = -1
    input_files = ''

    time_s = decrement_time(model_time, 0, obs_days_minus)
    time_e = increment_time(model_time, 0, obs_days_plus)

    ni = T_grid%ni; nj = T_grid%nj
    allocate(lon1d(ni*nj))
    allocate(lat1d(ni*nj))
    allocate(glon1d(ni*nj))
    allocate(glat1d(ni*nj))
    allocate(kdroot)
    do i = 1, ni; do j = 1, nj
      lon1d((j-1)*ni+i) = i
      lat1d((j-1)*ni+i) = j
      glon1d((j-1)*ni+i) = T_grid%x(i,j)
      glat1d((j-1)*ni+i) = T_grid%y(i,j)
    enddo; enddo
    call kd_init(kdroot, glon1d, glat1d)

    ! time window for DROP, MOORING and SATELLITE data respectively
    ! will be available from namelist
    do i=1,10
      data_seconds = data_window(i) * 3600
      time_window(i) = set_time(data_seconds,0)
    enddo
    
    nfiles = 0
    nrecs=0
    call mpp_open(unit, 'ocean_obs_table', action=MPP_RDONLY)
    read_obs: do while ( nfiles <= max_files )
       read (UNIT=unit, FMT='(A)', IOSTAT=io_status) record
       if ( io_status < 0 ) then
          exit read_obs
       else if ( io_status > 0 ) then
          cycle read_obs
       else
          nrecs = nrecs + 1
          if ( record(1:1) == '#' ) cycle read_obs
          read ( UNIT=record, FMT=*, IOSTAT=io_status ) tbl_entry
          if ( io_status < 0 ) then
             exit read_obs
          else if ( io_status > 0 ) then
             cycle read_obs
          else
             nfiles = nfiles + 1
             write (UNIT=stdout_unit, FMT='("Obs filename:",A," type:",A)') tbl_entry%filename,tbl_entry%file_type
             input_files(nfiles) = tbl_entry%filename
             select case ( trim(tbl_entry%file_type) )
             case ('profiles')
                filetype(nfiles) = PROFILE_FILE
             case ('surface')
                filetype(nfiles) = SURFACE_FILE
             case ('argo')
                filetype(nfiles) = ARGO_FILE
             case default
                call error_mesg('ocean_da_core_mod::init_observations', 'error in obs_table entry format', FATAL)
             end select
          end if
       end if
    end do read_obs
    if ( nfiles > max_files ) then
       call error_mesg('ocean_da_core_mod::init_observations', 'number of obs files exceeds max_files parameter', FATAL)
    end if
    CALL mpp_close(unit)

    if( .not. associated(Profiles) ) allocate(Profiles)

    if ( temp_obs ) then
       obs_variable = TEMP_ID
       if ( mpp_pe() == mpp_root_pe() ) then
          write (UNIT=stdout_unit, FMT='("TEMP_ID = ",I5)') TEMP_ID
       end if
       do n=1, nfiles
          select case ( filetype(n) )
          case (PROFILE_FILE)
             call open_profile_dataset(Profiles, Domain, T_grid, &
                     trim(input_files(n)), time_s, time_e, obs_variable)
          case (SURFACE_FILE)
             call open_oisst_dataset(Profiles, Domain, T_grid, &
                     trim(input_files(n)), time_s, time_e, obs_variable)
          case (ARGO_FILE)
             call open_argo_dataset(Profiles, Domain, T_grid, &
                     trim(input_files(n)), time_s, time_e, obs_variable)
          case default
             call error_mesg('ocean_da_core_mod::init_observations', 'filetype not currently supported for temp_obs', FATAL)
          end select
       end do
    end if

    if ( salt_obs ) then
       obs_variable = SALT_ID
       if ( mpp_pe() == mpp_root_pe() ) then
          write (UNIT=stdout_unit, FMT='("SALT_ID = ",I5)') SALT_ID
       end if
       do n=1, nfiles
          select case ( filetype(n) )
          case (PROFILE_FILE)
             call open_profile_dataset(Profiles, Domain, T_grid, &
                     trim(input_files(n)), time_s, time_e, obs_variable)
          case (SURFACE_FILE)
          case (ARGO_FILE)
             call open_argo_dataset(Profiles, Domain, T_grid, &
                     trim(input_files(n)), time_s, time_e, obs_variable)
          case default
             call error_mesg('ocean_da_core_mod::init_observations', 'filetype not currently supported for salt_obs', FATAL)
          end select
       end do
    end if

    ! Deallocate before exiting routine
    deallocate(kdroot)
    deallocate(lon1d, lat1d, glat1d, glon1d)
    deallocate(filetype, input_files)
  end subroutine ocean_da_core_init

  subroutine open_profile_dataset(Profiles, Domain, T_grid, &
                  filename, time_start, time_end, obs_variable, localize)
    type(ocean_profile_type), pointer :: Profiles
    !< This is an unstructured recursive list of profiles
    !< which are either within the localized domain corresponding
    !< to the Domain argument, or the global profile list
    type(domain2d), pointer, intent(in) :: Domain !< MOM grid type for the local domain
    type(grid_type), pointer, intent(in) :: T_grid !< MOM grid type for the local domain
    character(len=*), intent(in) :: filename !< filename containing profile data
    type(time_type), intent(in) :: time_start, time_end !< start and end times for the analysis
    integer, intent(in), optional :: obs_variable !< If present, then extract corresponding data
    !< from file, otherwise, extract all available data which.
    logical, intent(in), optional :: localize !< Localize the observations to the current computational domain

    real :: lon, lat, time, profile_error, rlink, flag_t, flag_s
    integer :: ni, nj, nk
    real :: ri0, rj0
    real, dimension(MAX_LEVELS_FILE) :: depth, data, t_flag, s_flag
    real, dimension(MAX_LINKS, MAX_LEVELS_FILE) :: data_bfr, depth_bfr, t_flag_bfr, s_flag_bfr
    type(ocean_profile_type), pointer :: Prof

    integer :: unit, ndim, nvar, natt, nstation, max_profiles
    integer :: stdout_unit
    integer :: inst_type, var_id
    integer :: station_count, station_link
    integer :: num_levs, k, kk, i, j, i0, j0, k0, nlevs, a, nn, nlinks
    integer :: yr, mon, day, hr, min, sec
    integer :: ii, jj

    logical :: data_is_local, localize_data, cont
    logical :: data_in_period, data_in_compute
    logical, dimension(MAX_LEVELS_FILE) :: flag
    logical, dimension(MAX_LINKS, MAX_LEVELS_FILE) :: flag_bfr

    character(len=32) :: fldname, axisname, time_units
    character(len=138) :: emsg_local

    type(time_type) :: obs_time, profile_time
    type(axistype), pointer :: depth_axis, station_axis
    type(axistype), allocatable, dimension(:), target :: axes
    type(fieldtype), allocatable, dimension(:), target :: fields
    type(fieldtype), pointer :: field_lon, field_lat, field_flag, field_flag_s, field_time, field_depth, field_t, field_s
    type(fieldtype), pointer :: field_t_error, field_s_error, field_link, field_t_flag, field_s_flag

    integer :: isc, iec, jsc, jec, isd, ied, jsd, jed
    integer :: isg, ieg, jsg, jeg, halox, haloy, lon_len, blk
    integer :: profile_count = 0
    type(horiz_interp_type) :: Interp
    real :: lon_out(1, 1), lat_out(1, 1)
    real :: lat_bound = 59.0
    integer :: inds(1), r_num
    real :: dist(1), frac_lon, frac_lat, frac_k
    real, dimension(6) :: coef
    integer, dimension(8) :: state_index

    Prof=>Profiles
    do while (associated(Prof%next))
      Prof=>Prof%next
    end do

    if ( PRESENT(localize) ) then
       localize_data = localize
    else
       localize_data = .true.
    end if

    ni = T_grid%ni; nj = T_grid%nj; nk = T_grid%nk
    call mpp_get_compute_domain(Domain, isc, iec, jsc, jec)
    call mpp_get_data_domain(Domain, isd, ied, jsd, jed)
    call mpp_get_global_domain(Domain, isg, ieg, jsg, jeg)
    lon_len = ied-isd+1
    blk = (jed-jsd+1)*lon_len
    stdout_unit = stdout()

    select case ( trim(filename(13:15)) )
      case ('PFL')
        inst_type = ODA_PFL
      case ('XBT')
        inst_type = ODA_XBT
      case ('MRB')
        inst_type = ODA_MRB
      case default
    end select

    var_id=-1
    if ( obs_variable == TEMP_ID .or. obs_variable == SALT_ID) then
       var_id = obs_variable
    end if

    call mpp_open(unit, filename, form=MPP_NETCDF, fileset=MPP_SINGLE, threading=MPP_MULTI, action=MPP_RDONLY)
    call mpp_get_info(unit, ndim, nvar, natt, nstation)

    write (UNIT=stdout_unit, FMT='("Opened profile dataset: ",A)') trim(filename)
    if (nstation .EQ. 0) then
      write(UNIT=stdout_unit, FMT='("There are ZERO records in this dataset.")')
      call mpp_close(unit)
      return
    end if

    ! get axis information
    allocate(axes(ndim))
    call mpp_get_axes(unit, axes)
    do i=1, ndim
       call mpp_get_atts(axes(i), name=axisname)
       select case ( trim(axisname) )
       case ('depth_index')
          depth_axis => axes(i)
       case ('station_index')
          station_axis => axes(i)
       end select
    end do

    ! get field information
    allocate(fields(nvar))
    call mpp_get_fields(unit, fields)
    field_lon=>NULL();field_lat=>NULL()
    field_flag=>NULL();field_flag_s=>NULL();field_time=>NULL()
    field_t=>NULL();field_t_error=>NULL();field_t_flag=>NULL()
    field_s=>NULL();field_s_error=>NULL();field_s_flag=>NULL()
    field_depth=>NULL();field_link=>NULL()
    do i=1, nvar
      call mpp_get_atts(fields(i), name=fldname)
      select case (trim(fldname))
      case ('longitude')
        field_lon => fields(i)
      case ('latitude')
        field_lat => fields(i)
      case ('profile_flag')
        field_flag => fields(i)
      case ('profile_flag_s')
        field_flag_s => fields(i)
      case ('time')
        field_time => fields(i)
      case ('temp')
        field_t => fields(i)
      case ('salt')
        field_s => fields(i)
      case ('depth')
        field_depth => fields(i)
      case ('link')
        field_link => fields(i)
      case ('temp_error')
        field_t_error => fields(i)
      case ('temp_flag')
        field_t_flag => fields(i)
      case ('salt_error')
        field_s_error => fields(i)
      case ('salt_flag')
        field_s_flag => fields(i)
      end select
    end do

    call mpp_get_atts(depth_axis, len=nlevs)

    if ( nlevs > MAX_LEVELS_FILE ) then
       call error_mesg('ocean_da_core_mod::open_profile_dataset', 'increase parameter MAX_LEVELS_FILE', FATAL)
    else if (nlevs < 1) then
       call error_mesg('ocean_da_core_mod::open_profile_dataset', 'Value of nlevs is less than 1.', FATAL)
    end if

    if ( .NOT.ASSOCIATED(field_t) .and. .not. ASSOCIATED(field_s) ) then
       call error_mesg('ocean_da_core_mod::open_profile_dataset',&
            & 'profile dataset not used because data not needed for Analysis', NOTE)
       return
    end if

    write(UNIT=stdout_unit, FMT='("There are ",I8," records in this dataset.")') nstation
    write(UNIT=stdout_unit, FMT='("Searching for profiles . . .")')

    call mpp_get_atts(field_time, units=time_units)

    station_count = 1
    cont = .true.

    do while ( cont )
       depth = MISSING_VALUE  ! snz add
       data = MISSING_VALUE   ! snz add

       call mpp_read(unit, field_lon, lon, tindex=station_count)
       call mpp_read(unit, field_lat, lat, tindex=station_count)
       call mpp_read(unit, field_time, time, tindex=station_count)
       if ( var_id == TEMP_ID ) then
         call mpp_read(unit, field_flag, flag_t, tindex=station_count)
       else if ( var_id == SALT_ID ) then
         call mpp_read(unit, field_flag_s, flag_s, tindex=station_count)
       end if

       data_is_local = .false.
       data_in_period = .false.

       if ( lon .lt. 0.0 ) lon = lon + 360.0
       if ( lon .gt. 360.0 ) lon = lon - 360.0
       if ( lon .gt. 60.0 ) lon = lon - 360.0

       if ( lat < obs_sbound(inst_type) .or. lat > obs_nbound(inst_type) ) then
         station_count = station_count + 1
         if ( station_count .gt. nstation ) cont = .false.
         cycle
       end if

       obs_time = get_cal_time(time, time_units, 'julian')
       profile_time = increment_time(obs_time, sec_offset(inst_type),day_offset(inst_type))
       if ( profile_time >= time_start .and. profile_time <= time_end ) data_in_period = .true.
       if ( .not. data_in_period ) then
         station_count = station_count + 1
         if ( station_count .gt. nstation ) cont = .false.
         cycle
       end if

       if ( localize_data ) then
         call kd_search_nnearest(kdroot, lon, lat, &
                 1, inds, dist, r_num, .false.)
         data_is_local = within_domain(lon1d(inds(1)), lat1d(inds(1)), isd+1, ied-1, jsd+1, jed-1, ni, nj)
         data_in_compute = within_domain(lon1d(inds(1)), lat1d(inds(1)), isc, iec, jsc, jec, ni, nj)
       else
         data_is_local = .true.
       end if

       if (.not. data_is_local) then
          station_count = station_count + 1
          if ( station_count .gt. nstation ) cont = .false.
          cycle
       end if

       profile_count = profile_count + 1
       type_count(inst_type) = type_count(inst_type)+1

       call mpp_read(unit, field_depth, depth(1:nlevs), tindex=station_count)
       if ( var_id == TEMP_ID ) then
         call mpp_read(unit, field_t, data(1:nlevs), tindex=station_count)
         call mpp_read(unit, field_t_error, profile_error, tindex=station_count)
         call mpp_read(unit, field_t_flag, t_flag(1:nlevs), tindex=station_count)
       else if ( var_id == SALT_ID ) then
         call mpp_read(unit, field_s, data(1:nlevs), tindex=station_count)
         call mpp_read(unit, field_s_error, profile_error, tindex=station_count)
         call mpp_read(unit, field_s_flag, s_flag(1:nlevs), tindex=station_count)
       end if
       call mpp_read(unit, field_link, rlink, tindex=station_count)

       num_levs = 0
       do k=1, MAX_LEVELS_FILE
          flag(k) = .true.

          if ( depth(k) > depth_cut ) depth(k) = MISSING_VALUE  ! snz add for rdat-hybn

          if ( var_id == TEMP_ID ) then
             if ( data(k) .eq. MISSING_VALUE .or.&
                  & depth(k) .eq. MISSING_VALUE .or. t_flag(k) .ne. 0.0 ) then
                flag(k) = .false.
             else
                num_levs = num_levs + 1
             end if
          else if ( var_id == SALT_ID ) then
             if ( data(k) .eq. MISSING_VALUE .or.&
                  & depth(k) .eq. MISSING_VALUE .or. s_flag(k) .ne. 0.0 ) then
                flag(k) = .false.
             else
                num_levs = num_levs+1
             end if
          end if
       end do

       ! large profile are stored externally in separate records
       ! read linked records and combine profile
       station_link = station_count + 1
       nlinks = 0
       do while ( rlink > 0.0 )
          nlinks = nlinks + 1

          if ( nlinks > MAX_LINKS ) then
             print *,'nlinks=',nlinks,'in ',filename
             call error_mesg('ocean_da_core_mod::open_profile_dataset', 'increase parameter MAX_LINKS', FATAL)
          end if

          depth_bfr(nlinks,:) = MISSING_VALUE
          data_bfr(nlinks,:) = MISSING_VALUE
          call mpp_read(unit, field_depth, depth_bfr(nlinks,1:nlevs), tindex=station_link)
          if ( var_id == TEMP_ID ) then
             call mpp_read(unit, field_t, data_bfr(nlinks,1:nlevs), tindex=station_link)
             call mpp_read(unit, field_t_flag, t_flag_bfr(nlinks,1:nlevs), tindex=station_link)
          else if ( var_id == SALT_ID ) then
             call mpp_read(unit, field_s, data_bfr(nlinks,1:nlevs), tindex=station_link)
             call mpp_read(unit, field_s_flag, s_flag_bfr(nlinks,1:nlevs), tindex=station_link)
          end if
          call mpp_read(unit, field_link, rlink, tindex=station_link)
          station_link = station_link + 1
       end do
       station_count = station_link ! set record counter to start of next profile

       if ( nlinks > 0 ) then
          do nn=1, nlinks
             do k=1, MAX_LEVELS_FILE
                flag_bfr(nn,k) = .true.

                if ( depth_bfr(nn,k) > depth_cut ) depth_bfr(nn,k) = MISSING_VALUE  ! snz add for rdat-hybn

                if ( var_id == TEMP_ID ) then
                   if ( data_bfr(nn,k) .eq. MISSING_VALUE .or.&
                        & depth_bfr(nn,k) .eq. MISSING_VALUE .or.&
                        & t_flag_bfr(nn,k) .ne. 0.0 ) then
                      flag_bfr(nn,k) = .false.
                   else
                      num_levs = num_levs+1
                   end if
                else if (var_id == SALT_ID) then
                   if ( data_bfr(nn,k) .eq. MISSING_VALUE .or.&
                        & depth_bfr(nn,k) .eq. MISSING_VALUE .or.&
                        & s_flag_bfr(nn,k) .ne. 0.0 ) then
                      flag_bfr(nn,k) = .false.
                   else
                      num_levs = num_levs+1
                   end if
                end if
             end do
          end do
       end if

       if ( num_levs == 0 ) then
          if ( station_count .gt. nstation ) cont = .false.
          cycle
       end if

       ! allocate profile structure content and put in data
       allocate(Prof%depth(num_levs));Prof%depth(:)=MISSING_VALUE
       allocate(Prof%data(num_levs));Prof%data(:)=MISSING_VALUE
       allocate(Prof%flag(num_levs));Prof%flag(:)=.false.
       Prof%variable = var_id
       Prof%inst_type = inst_type
       Prof%levels = num_levs
       Prof%lat = lat; Prof%lon = lon
       Prof%nbr_xi = lon1d(inds(1)); Prof%nbr_yi = lat1d(inds(1))
       Prof%nbr_dist = dist(1)
       Prof%compute = data_in_compute
       Prof%time_window = time_window(inst_type)
       Prof%impact_levels = impact_levels(inst_type)
       Prof%temp_to_salt = temp_to_salt(inst_type)
       Prof%salt_to_temp = salt_to_temp(inst_type)
       if ( var_id == TEMP_ID ) then
          Prof%obs_error = temp_error(inst_type)
          Prof%loc_dist = temp_dist(inst_type)
       else if ( var_id == SALT_ID ) then
          Prof%obs_error = salt_error(inst_type)
          Prof%loc_dist = salt_dist(inst_type)
       end if
       kk = 1
       do k=1, MAX_LEVELS_FILE
          if ( flag(k) ) then
             if ( kk > Prof%levels ) then
                call error_mesg('ocean_da_core_mod::open_profile_dataset',&
                     & 'Loop value "kk" is greater than profile levels', FATAL)
             end if
             Prof%depth(kk) = depth(k)
             Prof%data(kk) = data(k)
             Prof%flag(kk) = flag(k)
             kk = kk + 1
          end if
       end do

       do nn=1, nlinks
          do k=1, MAX_LEVELS_FILE
             if ( flag_bfr(nn,k) ) then
                if ( kk > Prof%levels ) then
                   call error_mesg('ocean_da_core_mod::open_profile_dataset',&
                        & 'Loop value "kk" is greater than profile levels (bfr loop)', FATAL)
                end if
                Prof%depth(kk) = depth_bfr(nn,k)
                Prof%data(kk) = data_bfr(nn,k)
                Prof%flag(kk) = flag_bfr(nn,k)
                kk = kk + 1
             end if
          end do
       end do

       Prof%time = profile_time


       if ( lat < lat_bound ) then ! calculate interpolation coefficients
          ri0 = frac_index(lon, T_grid%x(:,jsg))
          rj0 = frac_index(lat, T_grid%y(isg,:))
          i0 = floor(ri0)
          j0 = floor(rj0)
          if ( i0 > ieg .or. j0 > jeg ) then
             write (UNIT=emsg_local, FMT='("i0 = ",I8,", j0 = ",I8)') mpp_pe(), i0, j0
             call error_mesg('ocean_da_core_mod::open_profile_dataset',&
                  & 'For regular grids, either i0 > ieg or j0 > jeg.  '//trim(emsg_local), FATAL)
          end if
          Prof%i_index = ri0
          Prof%j_index = rj0
       else ! tripolar grids
          lon_out(1,1) = lon*DEG_TO_RAD
          lat_out(1,1) = lat*DEG_TO_RAD
          call horiz_interp_bilinear_new (Interp, T_grid%x*DEG_TO_RAD, T_grid%y*DEG_TO_RAD,&
               & lon_out, lat_out, new_search=.true., no_crash_when_not_found=.true.)

          if ( Interp%wti(1,1,2) < 1.0 ) then
             i0 = Interp%i_lon(1,1,1)
          else
             i0 = Interp%i_lon(1,1,2)
          end if
          if ( Interp%wtj(1,1,2) < 1.0 ) then
             j0 = Interp%j_lat(1,1,1)
          else
             j0 = Interp%j_lat(1,1,2)
          end if
          if ( i0 > ieg .or. j0 > jeg ) then
             write (UNIT=emsg_local, FMT='("i0 = ",I6,", j0 = ",I6)') mpp_pe(), i0, j0
             call error_mesg('ocean_da_core_mod::open_profile_dataset',&
                  & 'For tripolar grids, either i0 > ieg or j0 > jeg', FATAL)
          end if
          if ( Interp%wti(1,1,2) < 1.0 ) then
             Prof%i_index =Interp%i_lon(1,1,1) + Interp%wti(1,1,2)
          else
             Prof%i_index =Interp%i_lon(1,1,2)
          end if
          if (Interp%wtj(1,1,2) < 1.0) then
             Prof%j_index =Interp%j_lat(1,1,1) + Interp%wtj(1,1,2)
          else
             Prof%j_index =Interp%j_lat(1,1,2)
          end if
       end if ! interpolation coefficients

       Prof%accepted = .true.

       if ( var_id == TEMP_ID .and. flag_t /= 0.0 ) Prof%accepted = .false.
       if ( var_id == SALT_ID .and. flag_s /= 0.0 ) Prof%accepted = .false.
       ! removing some profiles due to bad placement at 0/0 lat/lon
       if ( abs(Prof%lat) < 0.001 .and. abs(Prof%lon) < 0.1 ) Prof%accepted = .false.

       call get_date(Prof%time, yr, mon, day, hr, min, sec)
       if(yr.eq.1991 .and. mon.eq.12 .and. day.ge.5 .and. day.le.8) Prof%accepted=.false.

       if (i0 < 1 .or. j0 < 1) then
          Prof%accepted = .false.
       else
          Prof%basin_mask = T_grid%basin_mask(lon1d(inds(1)),lat1d(inds(1)))
       end if

       if ( Prof%accepted ) then ! check surface land-sea mask and depth of ocean around profile location
          if ( i0 /= ieg .and. j0 /= jeg ) then
             if (T_grid%mask(i0,j0,1) == 0.0 .or.&
                  & T_grid%mask(i0+1,j0,1) == 0.0 .or.&
                  & T_grid%mask(i0,j0+1,1) == 0.0 .or.&
                  & T_grid%mask(i0+1,j0+1,1) == 0.0 ) then
                Prof%accepted = .false.
             end if
             if (T_grid%bathyT(i0,j0) < shelf_depth .or.&
                  & T_grid%bathyT(i0+1,j0) < shelf_depth .or.&
                  & T_grid%bathyT(i0,j0+1) < shelf_depth .or.&
                  & T_grid%bathyT(i0+1,j0+1) < shelf_depth ) then
                Prof%accepted = .false.
             end if
          else if ( i0 == ieg .and. j0 /= jeg ) then
             if (T_grid%mask(i0,j0,1) == 0.0 .or.&
                  & T_grid%mask(1,j0,1) == 0.0 .or.&
                  & T_grid%mask(i0,j0+1,1) == 0.0 .or.&
                  & T_grid%mask(1,j0+1,1) == 0.0 ) then
                Prof%accepted = .false.
             end if
             if (T_grid%bathyT(i0,j0) < shelf_depth .or.&
                  & T_grid%bathyT(1,j0) < shelf_depth .or.&
                  & T_grid%bathyT(i0,j0+1) < shelf_depth .or.&
                  & T_grid%bathyT(1,j0+1) < shelf_depth ) then
                Prof%accepted = .false.
             end if
          else if ( i0 /= ieg .and. j0 == jeg ) then
             if ( T_grid%mask(i0,j0,1) == 0.0 .or. T_grid%mask(i0+1,j0,1) == 0.0 ) then
                Prof%accepted = .false.
             end if
             if ( T_grid%bathyT(i0,j0) < shelf_depth .or.&
                   & T_grid%bathyT(i0+1,j0) < shelf_depth ) then
                Prof%accepted = .false.
             end if
          else
             if ( T_grid%mask(i0,j0,1) == 0.0 ) then
                Prof%accepted = .false.
             end if
             if ( T_grid%bathyT(i0,j0) < shelf_depth ) then
                Prof%accepted = .false.
             end if
          end if
       end if ! check surface land-sea mask and depth of ocean

       if ( Prof%accepted ) then ! determine vertical position and check mask at depth
          allocate(Prof%k_index(Prof%levels))
          do k=1, Prof%levels
             Prof%k_index(k) = frac_index(Prof%depth(k), (/T_grid%z(i0,j0,:)/))
             if ( Prof%k_index(k) < 1.0 ) then
                if ( Prof%depth(k) < T_grid%z(i0,j0,1) ) then
                   Prof%k_index(k) = 0.0
                else if ( Prof%depth(k) > T_grid%z(i0,j0,nk) ) then
                    Prof%k_index(k) = real(nk)
                    Prof%flag(k)=.false.
                end if
             end if
             if ( k > 3 ) then ! thinning the profile observations to a maximum of 3 within each layer
                if (floor(Prof%k_index(k)) == floor(Prof%k_index(k-3))) then
                    Prof%flag(k)=.false.
                end if
             end if
             if ( Prof%k_index(k) > real(nk) ) then
                call error_mesg('ocean_da_core_mod::open_profile_dataset', 'Profile k_index is greater than nk', FATAL)
             else if ( Prof%k_index(k) < 0.0 ) then
                call error_mesg('ocean_da_core_mod::open_profile_dataset', 'Profile k_index is less than 0', FATAL)
             end if
             k0 = floor(Prof%k_index(k))

             if ( k0 >= 1 ) THEN ! snz add
                if ( Prof%flag(k) ) then ! flag
                   if ( i0 /= ieg .and. j0 /= jeg ) then
                      if ( T_grid%mask(i0,j0,k0) == 0.0 .or.&
                           & T_grid%mask(i0+1,j0,k0) == 0.0 .or.&
                           & T_grid%mask(i0,j0+1,k0) == 0.0 .or.&
                           & T_grid%mask(i0+1,j0+1,k0) == 0.0 ) then
                         Prof%flag(k) = .false.
                      end if
                   else if ( i0 == ieg .and. j0 /= jeg ) then
                      if ( T_grid%mask(i0,j0,k0) == 0.0 .or.&
                           & T_grid%mask(1,j0,k0) == 0.0 .or.&
                           & T_grid%mask(i0,j0+1,k0) == 0.0 .or.&
                           & T_grid%mask(1,j0+1,k0) == 0.0) then
                         Prof%flag(k) = .false.
                      end if
                   else if ( i0 /= ieg .and. j0 == jeg ) then
                      if ( T_grid%mask(i0,j0,k0) == 0.0 .or.&
                           & T_grid%mask(i0+1,j0,k0) == 0.0) then
                         Prof%flag(k) = .false.
                      end if
                   else
                      if ( T_grid%mask(i0,j0,k0) == 0.0 ) then
                         Prof%flag(k) = .false.
                      end if
                   end if

                   if ( i0 /= ieg .and. j0 /= jeg) then
                      if ( T_grid%mask(i0,j0,k0+1) == 0.0 .or.&
                           & T_grid%mask(i0+1,j0,k0+1) == 0.0 .or.&
                           & T_grid%mask(i0,j0+1,k0+1) == 0.0 .or.&
                           & T_grid%mask(i0+1,j0+1,k0+1) == 0.0 ) then
                         Prof%flag(k) = .false.
                      end if
                   else if ( i0 == ieg .and. j0 /= jeg ) then
                      if ( T_grid%mask(i0,j0,k0+1) == 0.0 .or.&
                           & T_grid%mask(1,j0,k0+1) == 0.0 .or.&
                           & T_grid%mask(i0,j0+1,k0+1) == 0.0 .or.&
                           & T_grid%mask(1,j0+1,k0+1) == 0.0) then
                         Prof%flag(k) = .false.
                      end if
                   else if ( i0 /= ieg .and. j0 == jeg ) then
                      if ( T_grid%mask(i0,j0,k0+1) == 0.0 .or.&
                           & T_grid%mask(i0+1,j0,k0+1) == 0.0) then
                         Prof%flag(k) = .false.
                      end if
                   else
                      if ( T_grid%mask(i0,j0,k0+1) == 0.0 ) then
                         Prof%flag(k) = .false.
                      end if
                   end if

                   if ( Prof%data(k) == MISSING_VALUE &
                           .or. Prof%depth(k) == MISSING_VALUE ) then
                      Prof%flag(k) = .false.
                   end if
                end if ! flag
             end if ! snz add
          end do
       end if ! determine vertical position and check mask at depth

       if ( Prof%accepted ) then ! calculate forward operator indices and weights
         allocate(Prof%obs_def(Prof%levels))
         ii = i0; jj = j0
         frac_lat = Prof%j_index - jj
         frac_lon = Prof%i_index - ii

         coef(1) = (1.0 - frac_lon) * (1.0 - frac_lat)
         coef(2) = frac_lon * (1.0 - frac_lat)
         coef(3) = (1.0 - frac_lon) * frac_lat
         coef(4) = frac_lon * frac_lat

         if ( ied > ni .and. ii < isd ) ii = ii + ni
         if ( isd < 1 .and. ii > ied ) ii = ii - ni

         do k=1, Prof%levels
           k0 = floor(Prof%k_index(k))
           frac_k = Prof%k_index(k) - k0

           if ( k0 == 0 ) then
             state_index(1) = (jj-jsd)*lon_len + ii-isd + 1
             state_index(2) = (jj-jsd)*lon_len + ii-isd + 2
             state_index(3) = (jj-jsd+1)*lon_len + ii-isd + 1
             state_index(4) = (jj-jsd+1)*lon_len + ii-isd + 2
             state_index(5) = state_index(1)
             state_index(6) = state_index(2)
             state_index(7) = state_index(3)
             state_index(8) = state_index(4)
           else if (k0 == nk ) then
             state_index(1) = (k0-1)*blk + (jj-jsd)*lon_len+ii-isd+1
             state_index(2) = (k0-1)*blk + (jj-jsd)*lon_len+ii-isd+2
             state_index(3) = (k0-1)*blk + (jj-jsd+1)*lon_len+ii-isd+1
             state_index(4) = (k0-1)*blk + (jj-jsd+1)*lon_len+ii-isd+2
             state_index(5) = state_index(1)
             state_index(6) = state_index(2)
             state_index(7) = state_index(3)
             state_index(8) = state_index(4)
           else
             state_index(1) = (k0-1)*blk + (jj-jsd)*lon_len + ii-isd + 1
             state_index(2) = (k0-1)*blk + (jj-jsd)*lon_len + ii-isd + 2
             state_index(3) = (k0-1)*blk + (jj-jsd+1)*lon_len + ii-isd + 1
             state_index(4) = (k0-1)*blk + (jj-jsd+1)*lon_len + ii-isd + 2
             state_index(5) = k0*blk + (jj-jsd)*lon_len + ii-isd + 1
             state_index(6) = k0*blk + (jj-jsd)*lon_len + ii-isd + 2
             state_index(7) = k0*blk + (jj-jsd+1)*lon_len + ii-isd + 1
             state_index(8) = k0*blk + (jj-jsd+1)*lon_len + ii-isd + 2
           end if


           do i = 1, 8
             if ( state_index( i ) < 0 ) then
               write (UNIT=emsg_local, FMT='("state_index(",I1,") = ",I8," < 0 at &
                       [ii,jj] = [",I5,",",I5,"] within [isc,iec] = [",I5,",",I5,"] &
                       and [jsc,jec] = [",I5,",",I5,"], with halox = ",I5,", haloy = ",I5,", &
                       k0 = ",I5,", blk = ",I5,", nk = ",I5)') &
                       i, state_index( i ), ii, jj, isc, iec, jsc, jec, halox, haloy, k0, blk, nk
               call error_mesg('ocean_da_core_mod::open_profile_dataset', trim(emsg_local), FATAL)
             end if
           end do

           if ( frac_lon == 0.0 ) then
             state_index(2) = state_index(1)
             state_index(4) = state_index(3)
             state_index(6) = state_index(5)
             state_index(8) = state_index(7)
           end if

           if ( frac_lat == 0.0 ) then
             state_index(3) = state_index(1)
             state_index(4) = state_index(2)
             state_index(7) = state_index(5)
             state_index(8) = state_index(6)
           end if

           coef(5) = 1.0 - frac_k
           coef(6) = frac_k

           if ( frac_k == 0.0 ) then
             state_index(5) = state_index(1)
             state_index(6) = state_index(2)
             state_index(7) = state_index(3)
             state_index(8) = state_index(4)
           end if

           call def_forward_operator(8, state_index(1:8), coef(1:6), Prof%obs_def(k))
         end do
       endif ! calculate forward operator indices and weights

      !if ( var_id == TEMP_ID .and. profile_count > 0 ) call xbt_drop_rate_adjust(profiles(profile_count))

       if ( station_count .gt. nstation ) cont = .false.
       allocate(Prof%next) ! allocate next profile and link it to current one
       Prof%next%prev=>Prof
       Prof=>Prof%next
    end do

    call mpp_sync_self()
    call mpp_close(unit)
  end subroutine open_profile_dataset

  subroutine open_argo_dataset(Profiles, Domain, T_grid, &
                  filename, time_start, time_end, obs_variable, localize)
    type(ocean_profile_type), pointer :: Profiles
    !< This is an unstructured recursive list of profiles
    !< which are either within the localized domain corresponding
    !< to the Domain argument, or the global profile list
    type(domain2d), pointer, intent(in) :: Domain !< MOM grid type for the local domain
    type(grid_type), pointer, intent(in) :: T_grid !< MOM grid type for the local domain
    character(len=*), intent(in) :: filename !< filename containing profile data
    type(time_type), intent(in) :: time_start, time_end !< start and end times for the analysis
    integer, intent(in), optional :: obs_variable !< If present, then extract corresponding data
    !< from file, otherwise, extract all available data which.
    logical, intent(in), optional :: localize !< Localize the observations to the current computational domain

    real :: lon, lat, time, rlink, var_type
    integer :: ni, nj, nk
    real :: ri0, rj0
    real, dimension(MAX_LEVELS_FILE) :: depth, data
    real, dimension(MAX_LINKS, MAX_LEVELS_FILE) :: data_bfr, depth_bfr
    type(ocean_profile_type), pointer :: Prof

    integer :: unit, ndim, nvar, natt, nstation, max_profiles
    integer :: stdout_unit
    integer :: inst_type, var_id
    integer :: station_count, station_link
    integer :: num_levs, k, kk, i, j, i0, j0, k0, nlevs, a, nn, nlinks
    integer :: yr, mon, day, hr, min, sec
    integer :: ii, jj

    logical :: data_is_local, localize_data, cont
    logical :: data_in_period, data_in_compute
    logical, dimension(MAX_LEVELS_FILE) :: flag
    logical, dimension(MAX_LINKS, MAX_LEVELS_FILE) :: flag_bfr

    character(len=32) :: fldname, axisname, time_units
    character(len=138) :: emsg_local

    type(time_type) :: obs_time, profile_time
    type(axistype), pointer :: depth_axis, station_axis
    type(axistype), allocatable, dimension(:), target :: axes
    type(fieldtype), allocatable, dimension(:), target :: fields
    type(fieldtype), pointer :: field_lon, field_lat, field_time, field_depth
    type(fieldtype), pointer :: field_t, field_s, field_link, field_var_type

    integer :: isc, iec, jsc, jec, isd, ied, jsd, jed
    integer :: isg, ieg, jsg, jeg, halox, haloy, lon_len, blk
    integer :: profile_count = 0, temp_count = 0
    type(horiz_interp_type) :: Interp
    real :: lon_out(1, 1), lat_out(1, 1)
    real :: lat_bound = 59.0
    integer :: inds(1), r_num
    real :: dist(1), frac_lon, frac_lat, frac_k
    real, dimension(6) :: coef
    integer, dimension(8) :: state_index

    Prof=>Profiles
    do while (associated(Prof%next))
      Prof=>Prof%next
    end do

    if ( PRESENT(localize) ) then
       localize_data = localize
    else
       localize_data = .true.
    end if

    ni = T_grid%ni; nj = T_grid%nj; nk = T_grid%nk
    call mpp_get_compute_domain(Domain, isc, iec, jsc, jec)
    call mpp_get_data_domain(Domain, isd, ied, jsd, jed)
    call mpp_get_global_domain(Domain, isg, ieg, jsg, jeg)
    lon_len = ied-isd+1
    blk = (jed-jsd+1)*lon_len
    stdout_unit = stdout()

    inst_type = ODA_PFL

    var_id=-1
    if ( obs_variable == TEMP_ID .or. obs_variable == SALT_ID) then
       var_id = obs_variable
    end if

    call mpp_open(unit, filename, form=MPP_NETCDF, fileset=MPP_SINGLE, threading=MPP_MULTI, action=MPP_RDONLY)
    call mpp_get_info(unit, ndim, nvar, natt, nstation)

    write (UNIT=stdout_unit, FMT='("Opened argo dataset: ",A)') trim(filename)
    if (nstation .EQ. 0) then
      write(UNIT=stdout_unit, FMT='("There are ZERO records in this dataset.")')
      call mpp_close(unit)
      return
    end if

    ! get axis information
    allocate(axes(ndim))
    call mpp_get_axes(unit, axes)
    do i=1, ndim
       call mpp_get_atts(axes(i), name=axisname)
       select case ( trim(axisname) )
       case ('depth_index')
          depth_axis => axes(i)
       case ('station_index')
          station_axis => axes(i)
       end select
    end do

    ! get field information
    allocate(fields(nvar))
    call mpp_get_fields(unit, fields)
    field_lon=>NULL();field_lat=>NULL();field_time=>NULL()
    field_t=>NULL();field_s=>NULL();field_var_type=>NULL()
    field_depth=>NULL();field_link=>NULL()
    do i=1, nvar
      call mpp_get_atts(fields(i), name=fldname)
      select case (trim(fldname))
      case ('longitude')
        field_lon => fields(i)
      case ('latitude')
        field_lat => fields(i)
      case ('var_type')
        field_var_type => fields(i)
      case ('time')
        field_time => fields(i)
      case ('temp')
        field_t => fields(i)
      case ('salt')
        field_s => fields(i)
      case ('depth')
        field_depth => fields(i)
      case ('link')
        field_link => fields(i)
      end select
    end do

    call mpp_get_atts(depth_axis, len=nlevs)

    if ( nlevs > MAX_LEVELS_FILE ) then
       call error_mesg('ocean_da_core_mod::open_profile_dataset', 'increase parameter MAX_LEVELS_FILE', FATAL)
    else if (nlevs < 1) then
       call error_mesg('ocean_da_core_mod::open_profile_dataset', 'Value of nlevs is less than 1.', FATAL)
    end if

    if ( .NOT.ASSOCIATED(field_t) .and. .not. ASSOCIATED(field_s) ) then
       call error_mesg('ocean_da_core_mod::open_profile_dataset',&
            & 'profile dataset not used because data not needed for Analysis', NOTE)
       return
    end if

    write(UNIT=stdout_unit, FMT='("There are ",I8," records in this dataset.")') nstation
    write(UNIT=stdout_unit, FMT='("Searching for profiles . . .")')

    call mpp_get_atts(field_time, units=time_units)

    station_count = 1
    cont = .true.

    do while ( cont )
       depth = MISSING_VALUE  ! snz add
       data = MISSING_VALUE   ! snz add

       call mpp_read(unit, field_lon, lon, tindex=station_count)
       call mpp_read(unit, field_lat, lat, tindex=station_count)
       call mpp_read(unit, field_time, time, tindex=station_count)
       call mpp_read(unit, field_var_type, var_type, tindex=station_count)
       if (var_id==SALT_ID .and. var_type==1) then
          temp_count = temp_count + 1
          station_count = station_count + 1
          if ( station_count .gt. nstation ) cont = .false.
          cycle
       end if

       data_is_local = .false.
       data_in_period = .false.

       if ( lon .lt. 0.0 ) lon = lon + 360.0
       if ( lon .gt. 360.0 ) lon = lon - 360.0
       if ( lon .gt. 60.0 ) lon = lon - 360.0

       if ( lat < obs_sbound(inst_type) .or. lat > obs_nbound(inst_type) ) then
         station_count = station_count + 1
         if ( station_count .gt. nstation ) cont = .false.
         cycle
       end if

       obs_time = get_cal_time(time, time_units, 'julian')
       profile_time = increment_time(obs_time, sec_offset(inst_type),day_offset(inst_type))
       if ( profile_time >= time_start .and. profile_time <= time_end ) data_in_period = .true.
       if ( .not. data_in_period ) then
         station_count = station_count + 1
         if ( station_count .gt. nstation ) cont = .false.
         cycle
       end if

       if ( localize_data ) then
         call kd_search_nnearest(kdroot, lon, lat, &
                 1, inds, dist, r_num, .false.)
         data_is_local = within_domain(lon1d(inds(1)), lat1d(inds(1)), isd+1, ied-1, jsd+1, jed-1, ni, nj)
         data_in_compute = within_domain(lon1d(inds(1)), lat1d(inds(1)), isc, iec, jsc, jec, ni, nj)
       else
         data_is_local = .true.
       end if

       if (.not. data_is_local) then
          station_count = station_count + 1
          if ( station_count .gt. nstation ) cont = .false.
          cycle
       end if

       profile_count = profile_count + 1
       type_count(inst_type) = type_count(inst_type)+1

       call mpp_read(unit, field_depth, depth(1:nlevs), tindex=station_count)
       if ( var_id == TEMP_ID ) then
         call mpp_read(unit, field_t, data(1:nlevs), tindex=station_count)
       else if ( var_id == SALT_ID ) then
         call mpp_read(unit, field_s, data(1:nlevs), tindex=station_count)
       end if
       call mpp_read(unit, field_link, rlink, tindex=station_count)

       num_levs = 0
       do k=1, MAX_LEVELS_FILE
          flag(k) = .true.
          if ( depth(k) > depth_cut ) depth(k) = MISSING_VALUE  ! snz add for rdat-hybn
          if ( data(k) .eq. MISSING_VALUE .or. depth(k) .eq. MISSING_VALUE) then
             flag(k) = .false.
          else
             num_levs = num_levs + 1
          end if
       end do

       ! large profile are stored externally in separate records
       ! read linked records and combine profile
       station_link = station_count + 1
       nlinks = 0
       do while ( rlink > 0.0 )
          nlinks = nlinks + 1

          if ( nlinks > MAX_LINKS ) then
             print *,'nlinks=',nlinks,'in ',filename
             call error_mesg('ocean_da_core_mod::open_profile_dataset', 'increase parameter MAX_LINKS', FATAL)
          end if

          depth_bfr(nlinks,:) = MISSING_VALUE
          data_bfr(nlinks,:) = MISSING_VALUE
          call mpp_read(unit, field_depth, depth_bfr(nlinks,1:nlevs), tindex=station_link)
          if ( var_id == TEMP_ID ) then
             call mpp_read(unit, field_t, data_bfr(nlinks,1:nlevs), tindex=station_link)
          else if ( var_id == SALT_ID ) then
             call mpp_read(unit, field_s, data_bfr(nlinks,1:nlevs), tindex=station_link)
          end if
          call mpp_read(unit, field_link, rlink, tindex=station_link)
          station_link = station_link + 1
       end do
       station_count = station_link ! set record counter to start of next profile

       if ( nlinks > 0 ) then
          do nn=1, nlinks
             do k=1, MAX_LEVELS_FILE
                flag_bfr(nn,k) = .true.

                if ( depth_bfr(nn,k) > depth_cut ) depth_bfr(nn,k) = MISSING_VALUE

                if ( data_bfr(nn,k) .eq. MISSING_VALUE .or. &
                        depth_bfr(nn,k) .eq. MISSING_VALUE) then
                   flag_bfr(nn,k) = .false.
                else
                   num_levs = num_levs+1
                end if
             end do
          end do
       end if

       if ( num_levs == 0 ) then
          if ( station_count .gt. nstation ) cont = .false.
          cycle
       end if

       ! allocate profile structure content and put in data
       allocate(Prof%depth(num_levs));Prof%depth(:)=MISSING_VALUE
       allocate(Prof%data(num_levs));Prof%data(:)=MISSING_VALUE
       allocate(Prof%flag(num_levs));Prof%flag(:)=.false.
       Prof%variable = var_id
       Prof%inst_type = inst_type
       Prof%levels = num_levs
       Prof%lat = lat; Prof%lon = lon
       Prof%nbr_xi = lon1d(inds(1)); Prof%nbr_yi = lat1d(inds(1))
       Prof%nbr_dist = dist(1)
       Prof%compute = data_in_compute
       Prof%time_window = time_window(inst_type)
       Prof%impact_levels = impact_levels(inst_type)
       Prof%temp_to_salt = temp_to_salt(inst_type)
       Prof%salt_to_temp = salt_to_temp(inst_type)
       if ( var_id == TEMP_ID ) then
          Prof%obs_error = temp_error(inst_type)
          Prof%loc_dist = temp_dist(inst_type)
       else if ( var_id == SALT_ID ) then
          Prof%obs_error = salt_error(inst_type)
          Prof%loc_dist = salt_dist(inst_type)
       end if
       kk = 1
       do k=1, MAX_LEVELS_FILE
          if ( flag(k) ) then
             if ( kk > Prof%levels ) then
                call error_mesg('ocean_da_core_mod::open_profile_dataset',&
                     & 'Loop value "kk" is greater than profile levels', FATAL)
             end if
             Prof%depth(kk) = depth(k)
             Prof%data(kk) = data(k)
             Prof%flag(kk) = flag(k)
             kk = kk + 1
          end if
       end do

       do nn=1, nlinks
          do k=1, MAX_LEVELS_FILE
             if ( flag_bfr(nn,k) ) then
                if ( kk > Prof%levels ) then
                   call error_mesg('ocean_da_core_mod::open_profile_dataset',&
                        & 'Loop value "kk" is greater than profile levels (bfr loop)', FATAL)
                end if
                Prof%depth(kk) = depth_bfr(nn,k)
                Prof%data(kk) = data_bfr(nn,k)
                Prof%flag(kk) = flag_bfr(nn,k)
                kk = kk + 1
             end if
          end do
       end do

       Prof%time = profile_time

       if ( lat < lat_bound ) then ! calculate interpolation coefficients
          ri0 = frac_index(lon, T_grid%x(:,jsg))
          rj0 = frac_index(lat, T_grid%y(isg,:))
          i0 = floor(ri0)
          j0 = floor(rj0)
          if ( i0 > ieg .or. j0 > jeg ) then
             write (UNIT=emsg_local, FMT='("i0 = ",I8,", j0 = ",I8)') mpp_pe(), i0, j0
             call error_mesg('ocean_da_core_mod::open_profile_dataset',&
                  & 'For regular grids, either i0 > ieg or j0 > jeg.  '//trim(emsg_local), FATAL)
          end if
          Prof%i_index = ri0
          Prof%j_index = rj0
       else ! tripolar grids
          lon_out(1,1) = lon*DEG_TO_RAD
          lat_out(1,1) = lat*DEG_TO_RAD
          call horiz_interp_bilinear_new (Interp, T_grid%x*DEG_TO_RAD, T_grid%y*DEG_TO_RAD,&
               & lon_out, lat_out, new_search=.true., no_crash_when_not_found=.true.)

          if ( Interp%wti(1,1,2) < 1.0 ) then
             i0 = Interp%i_lon(1,1,1)
          else
             i0 = Interp%i_lon(1,1,2)
          end if
          if ( Interp%wtj(1,1,2) < 1.0 ) then
             j0 = Interp%j_lat(1,1,1)
          else
             j0 = Interp%j_lat(1,1,2)
          end if
          if ( i0 > ieg .or. j0 > jeg ) then
             write (UNIT=emsg_local, FMT='("i0 = ",I6,", j0 = ",I6)') mpp_pe(), i0, j0
             call error_mesg('ocean_da_core_mod::open_profile_dataset',&
                  & 'For tripolar grids, either i0 > ieg or j0 > jeg', FATAL)
          end if
          if ( Interp%wti(1,1,2) < 1.0 ) then
             Prof%i_index =Interp%i_lon(1,1,1) + Interp%wti(1,1,2)
          else
             Prof%i_index =Interp%i_lon(1,1,2)
          end if
          if (Interp%wtj(1,1,2) < 1.0) then
             Prof%j_index =Interp%j_lat(1,1,1) + Interp%wtj(1,1,2)
          else
             Prof%j_index =Interp%j_lat(1,1,2)
          end if
       end if ! interpolation coefficients

       Prof%accepted = .true.

       if (i0 < 1 .or. j0 < 1) then
          Prof%accepted = .false.
       else
          Prof%basin_mask = T_grid%basin_mask(lon1d(inds(1)),lat1d(inds(1)))
       end if

       if ( Prof%accepted ) then ! check surface land-sea mask and depth of ocean around profile location
          if ( i0 /= ieg .and. j0 /= jeg ) then
             if (T_grid%mask(i0,j0,1) == 0.0 .or.&
                  & T_grid%mask(i0+1,j0,1) == 0.0 .or.&
                  & T_grid%mask(i0,j0+1,1) == 0.0 .or.&
                  & T_grid%mask(i0+1,j0+1,1) == 0.0 ) then
                Prof%accepted = .false.
             end if
             if (T_grid%bathyT(i0,j0) < shelf_depth .or.&
                  & T_grid%bathyT(i0+1,j0) < shelf_depth .or.&
                  & T_grid%bathyT(i0,j0+1) < shelf_depth .or.&
                  & T_grid%bathyT(i0+1,j0+1) < shelf_depth ) then
                Prof%accepted = .false.
             end if
          else if ( i0 == ieg .and. j0 /= jeg ) then
             if (T_grid%mask(i0,j0,1) == 0.0 .or.&
                  & T_grid%mask(1,j0,1) == 0.0 .or.&
                  & T_grid%mask(i0,j0+1,1) == 0.0 .or.&
                  & T_grid%mask(1,j0+1,1) == 0.0 ) then
                Prof%accepted = .false.
             end if
             if (T_grid%bathyT(i0,j0) < shelf_depth .or.&
                  & T_grid%bathyT(1,j0) < shelf_depth .or.&
                  & T_grid%bathyT(i0,j0+1) < shelf_depth .or.&
                  & T_grid%bathyT(1,j0+1) < shelf_depth ) then
                Prof%accepted = .false.
             end if
          else if ( i0 /= ieg .and. j0 == jeg ) then
             if ( T_grid%mask(i0,j0,1) == 0.0 .or. T_grid%mask(i0+1,j0,1) == 0.0 ) then
                Prof%accepted = .false.
             end if
             if ( T_grid%bathyT(i0,j0) < shelf_depth .or.&
                   & T_grid%bathyT(i0+1,j0) < shelf_depth ) then
                Prof%accepted = .false.
             end if
          else
             if ( T_grid%mask(i0,j0,1) == 0.0 ) then
                Prof%accepted = .false.
             end if
             if ( T_grid%bathyT(i0,j0) < shelf_depth ) then
                Prof%accepted = .false.
             end if
          end if
       end if ! check surface land-sea mask and depth of ocean

       if ( Prof%accepted ) then ! determine vertical position and check mask at depth
          allocate(Prof%k_index(Prof%levels))
          do k=1, Prof%levels
             Prof%k_index(k) = frac_index(Prof%depth(k), (/T_grid%z(i0,j0,:)/))
             if ( Prof%k_index(k) < 1.0 ) then
                if ( Prof%depth(k) < T_grid%z(i0,j0,1) ) then
                   Prof%k_index(k) = 0.0
                else if ( Prof%depth(k) > T_grid%z(i0,j0,nk) ) then
                    Prof%k_index(k) = real(nk)
                    Prof%flag(k)=.false.
                end if
             end if
             if ( k > 3 ) then ! thinning the profile observations to a maximum of 3 within each layer
                if (floor(Prof%k_index(k)) == floor(Prof%k_index(k-3))) then
                    Prof%flag(k)=.false.
                end if
             end if
             if ( Prof%k_index(k) > real(nk) ) then
                call error_mesg('ocean_da_core_mod::open_profile_dataset', 'Profile k_index is greater than nk', FATAL)
             else if ( Prof%k_index(k) < 0.0 ) then
                call error_mesg('ocean_da_core_mod::open_profile_dataset', 'Profile k_index is less than 0', FATAL)
             end if
             k0 = floor(Prof%k_index(k))

             if ( k0 >= 1 ) THEN ! snz add
                if ( Prof%flag(k) ) then ! flag
                   if ( i0 /= ieg .and. j0 /= jeg ) then
                      if ( T_grid%mask(i0,j0,k0) == 0.0 .or.&
                           & T_grid%mask(i0+1,j0,k0) == 0.0 .or.&
                           & T_grid%mask(i0,j0+1,k0) == 0.0 .or.&
                           & T_grid%mask(i0+1,j0+1,k0) == 0.0 ) then
                         Prof%flag(k) = .false.
                      end if
                   else if ( i0 == ieg .and. j0 /= jeg ) then
                      if ( T_grid%mask(i0,j0,k0) == 0.0 .or.&
                           & T_grid%mask(1,j0,k0) == 0.0 .or.&
                           & T_grid%mask(i0,j0+1,k0) == 0.0 .or.&
                           & T_grid%mask(1,j0+1,k0) == 0.0) then
                         Prof%flag(k) = .false.
                      end if
                   else if ( i0 /= ieg .and. j0 == jeg ) then
                      if ( T_grid%mask(i0,j0,k0) == 0.0 .or.&
                           & T_grid%mask(i0+1,j0,k0) == 0.0) then
                         Prof%flag(k) = .false.
                      end if
                   else
                      if ( T_grid%mask(i0,j0,k0) == 0.0 ) then
                         Prof%flag(k) = .false.
                      end if
                   end if

                   if ( i0 /= ieg .and. j0 /= jeg) then
                      if ( T_grid%mask(i0,j0,k0+1) == 0.0 .or.&
                           & T_grid%mask(i0+1,j0,k0+1) == 0.0 .or.&
                           & T_grid%mask(i0,j0+1,k0+1) == 0.0 .or.&
                           & T_grid%mask(i0+1,j0+1,k0+1) == 0.0 ) then
                         Prof%flag(k) = .false.
                      end if
                   else if ( i0 == ieg .and. j0 /= jeg ) then
                      if ( T_grid%mask(i0,j0,k0+1) == 0.0 .or.&
                           & T_grid%mask(1,j0,k0+1) == 0.0 .or.&
                           & T_grid%mask(i0,j0+1,k0+1) == 0.0 .or.&
                           & T_grid%mask(1,j0+1,k0+1) == 0.0) then
                         Prof%flag(k) = .false.
                      end if
                   else if ( i0 /= ieg .and. j0 == jeg ) then
                      if ( T_grid%mask(i0,j0,k0+1) == 0.0 .or.&
                           & T_grid%mask(i0+1,j0,k0+1) == 0.0) then
                         Prof%flag(k) = .false.
                      end if
                   else
                      if ( T_grid%mask(i0,j0,k0+1) == 0.0 ) then
                         Prof%flag(k) = .false.
                      end if
                   end if

                   if ( Prof%data(k) == MISSING_VALUE &
                           .or. Prof%depth(k) == MISSING_VALUE ) then
                      Prof%flag(k) = .false.
                   end if
                end if ! flag
             end if ! snz add
          end do
       end if ! determine vertical position and check mask at depth

       if ( Prof%accepted ) then ! calculate forward operator indices and weights
         allocate(Prof%obs_def(Prof%levels))
         ii = i0; jj = j0
         frac_lat = Prof%j_index - jj
         frac_lon = Prof%i_index - ii

         coef(1) = (1.0 - frac_lon) * (1.0 - frac_lat)
         coef(2) = frac_lon * (1.0 - frac_lat)
         coef(3) = (1.0 - frac_lon) * frac_lat
         coef(4) = frac_lon * frac_lat

         if ( ied > ni .and. ii < isd ) ii = ii + ni
         if ( isd < 1 .and. ii > ied ) ii = ii - ni

         do k=1, Prof%levels
           k0 = floor(Prof%k_index(k))
           frac_k = Prof%k_index(k) - k0

           if ( k0 == 0 ) then
             state_index(1) = (jj-jsd)*lon_len + ii-isd + 1
             state_index(2) = (jj-jsd)*lon_len + ii-isd + 2
             state_index(3) = (jj-jsd+1)*lon_len + ii-isd + 1
             state_index(4) = (jj-jsd+1)*lon_len + ii-isd + 2
             state_index(5) = state_index(1)
             state_index(6) = state_index(2)
             state_index(7) = state_index(3)
             state_index(8) = state_index(4)
           else if (k0 == nk ) then
             state_index(1) = (k0-1)*blk + (jj-jsd)*lon_len+ii-isd+1
             state_index(2) = (k0-1)*blk + (jj-jsd)*lon_len+ii-isd+2
             state_index(3) = (k0-1)*blk + (jj-jsd+1)*lon_len+ii-isd+1
             state_index(4) = (k0-1)*blk + (jj-jsd+1)*lon_len+ii-isd+2
             state_index(5) = state_index(1)
             state_index(6) = state_index(2)
             state_index(7) = state_index(3)
             state_index(8) = state_index(4)
           else
             state_index(1) = (k0-1)*blk + (jj-jsd)*lon_len + ii-isd + 1
             state_index(2) = (k0-1)*blk + (jj-jsd)*lon_len + ii-isd + 2
             state_index(3) = (k0-1)*blk + (jj-jsd+1)*lon_len + ii-isd + 1
             state_index(4) = (k0-1)*blk + (jj-jsd+1)*lon_len + ii-isd + 2
             state_index(5) = k0*blk + (jj-jsd)*lon_len + ii-isd + 1
             state_index(6) = k0*blk + (jj-jsd)*lon_len + ii-isd + 2
             state_index(7) = k0*blk + (jj-jsd+1)*lon_len + ii-isd + 1
             state_index(8) = k0*blk + (jj-jsd+1)*lon_len + ii-isd + 2
           end if


           do i = 1, 8
             if ( state_index( i ) < 0 ) then
               write (UNIT=emsg_local, FMT='("state_index(",I1,") = ",I8," < 0 at &
                       [ii,jj] = [",I5,",",I5,"] within [isc,iec] = [",I5,",",I5,"] &
                       and [jsc,jec] = [",I5,",",I5,"], with halox = ",I5,", haloy = ",I5,", &
                       k0 = ",I5,", blk = ",I5,", nk = ",I5)') &
                       i, state_index( i ), ii, jj, isc, iec, jsc, jec, halox, haloy, k0, blk, nk
               call error_mesg('ocean_da_core_mod::open_profile_dataset', trim(emsg_local), FATAL)
             end if
           end do

           if ( frac_lon == 0.0 ) then
             state_index(2) = state_index(1)
             state_index(4) = state_index(3)
             state_index(6) = state_index(5)
             state_index(8) = state_index(7)
           end if

           if ( frac_lat == 0.0 ) then
             state_index(3) = state_index(1)
             state_index(4) = state_index(2)
             state_index(7) = state_index(5)
             state_index(8) = state_index(6)
           end if

           coef(5) = 1.0 - frac_k
           coef(6) = frac_k

           if ( frac_k == 0.0 ) then
             state_index(5) = state_index(1)
             state_index(6) = state_index(2)
             state_index(7) = state_index(3)
             state_index(8) = state_index(4)
           end if

           call def_forward_operator(8, state_index(1:8), coef(1:6), Prof%obs_def(k))
         end do
       endif ! calculate forward operator indices and weights

      !if ( var_id == TEMP_ID .and. profile_count > 0 ) call xbt_drop_rate_adjust(profiles(profile_count))

       if ( station_count .gt. nstation ) cont = .false.
       allocate(Prof%next) ! allocate next profile and link it to current one
       Prof%next%prev=>Prof
       Prof=>Prof%next
    end do

    call mpp_sync_self()
    call mpp_close(unit)
  end subroutine open_argo_dataset

  subroutine open_oisst_dataset(Profiles, Domain, T_grid, &
                  filename, time_start, time_end, obs_variable, localize)
    type(ocean_profile_type), pointer :: Profiles
    !< This is an unstructured recursive list of profiles
    !< which are either within the localized domain corresponding
    !< to the Domain argument, or the global profile list
    type(domain2d), pointer, intent(in) :: Domain !< MOM grid type for the local domain
    type(grid_type), pointer, intent(in) :: T_grid !< MOM grid type for the local domain
    character(len=*), intent(in) :: filename !< filename containing profile data
    type(time_type), intent(in) :: time_start, time_end !< start and end times for the analysis
    integer, intent(in), optional :: obs_variable !< If present, then extract corresponding data
    !< from file, otherwise, extract all available data which.
    logical, intent(in), optional :: localize !< Localize the observations to the current computational domain

    real :: lon, lat, time
    integer :: nlat, nlon, ntime
    integer :: ni, nj, nk
    real :: ri0, rj0
    real :: depth, data
    logical :: flag
    type(ocean_profile_type), pointer :: Prof

    integer :: unit, ndim, nvar, natt, max_profiles
    integer :: stdout_unit
    integer :: inst_type, var_id
    integer :: num_levs, k, kk, i, j, i0, j0, k0, nlevs, a, nn, nlinks
    integer :: yr, mon, day, hr, min, sec
    integer :: ii, jj

    logical :: data_is_local, localize_data, data_in_period

    character(len=32) :: fldname, axisname, time_units
    character(len=138) :: emsg_local

    type(time_type) :: surface_time, obs_time
    type(axistype), pointer :: lon_axis, lat_axis, time_axis
    type(axistype), allocatable, dimension(:), target :: axes
    type(fieldtype), allocatable, dimension(:), target :: fields
    type(fieldtype), pointer :: field_sst

    real, allocatable, dimension(:) :: lons, lats, times
    real, allocatable, dimension(:,:) :: sfc_obs
    integer :: sfc_size(3)

    integer :: isc, iec, jsc, jec, isd, ied, jsd, jed
    integer :: isg, ieg, jsg, jeg, halox, haloy, lon_len, blk
    integer :: surface_count = 0
    type(horiz_interp_type) :: Interp
    real :: lon_out(1, 1), lat_out(1, 1)
    real :: lat_bound = 59.0
    integer :: inds(1), r_num
    real :: dist(1), frac_lon, frac_lat, frac_k
    real, dimension(6) :: coef
    integer, dimension(8) :: state_index

    Prof=>Profiles
    do while (associated(Prof%next))
      Prof=>Prof%next
    end do

    if ( PRESENT(localize) ) then
       localize_data = localize
    else
       localize_data = .true.
    end if

    ni = T_grid%ni; nj = T_grid%nj; nk = T_grid%nk
    call mpp_get_compute_domain(Domain, isc, iec, jsc, jec)
    call mpp_get_data_domain(Domain, isd, ied, jsd, jed)
    call mpp_get_global_domain(Domain, isg, ieg, jsg, jeg)
    lon_len = ied-isd+1
    blk = (jed-jsd+1)*lon_len
    stdout_unit = stdout()

    inst_type = ODA_OISST
    var_id = obs_variable

    call mpp_open(unit, filename, form=MPP_NETCDF, fileset=MPP_SINGLE, threading=MPP_MULTI, action=MPP_RDONLY)
    call mpp_get_info(unit, ndim, nvar, natt, ntime)

    write (UNIT=stdout_unit, FMT='("Opened surface dataset: ",A)') trim(filename)

    ! get axis information
    allocate(axes(ndim))
    call mpp_get_axes(unit, axes)
    do i=1, ndim
       call mpp_get_atts(axes(i), name=axisname)
       select case ( trim(axisname) )
       case ('lon')
          lon_axis => axes(i)
       case ('lat')
          lat_axis => axes(i)
       case ('time')
          time_axis => axes(i)
       end select
    end do

    call mpp_get_atts(lon_axis,len=nlon)
    call mpp_get_atts(lat_axis,len=nlat)
    call mpp_get_atts(time_axis,len=ntime)
    call mpp_get_atts(time_axis, units=time_units)

    allocate(lons(nlon), lats(nlat), times(ntime))
    allocate(sfc_obs(nlon,nlat))

    call mpp_get_axis_data(lon_axis, lons)
    call mpp_get_axis_data(lat_axis, lats)
    call mpp_get_times(unit, times)

    ! get field information
    allocate(fields(nvar))
    call mpp_get_fields(unit, fields)
    field_sst=>NULL()
    do i=1, nvar
      call mpp_get_atts(fields(i), name=fldname)
      select case (trim(fldname))
      case ('sst')
        field_sst => fields(i)
      end select
    end do

    call mpp_get_atts(field_sst, siz=sfc_size)
    write(UNIT=stdout_unit, FMT='("Searching for surface obs . . .")')

    num_levs = 1
    do k=1, ntime
      data_in_period = .false.
      time = times(k)
      obs_time = get_cal_time(time, time_units, 'julian')
      ! Weekly OISST is timed at beginning of the 7-day period, so increase time by 3.5 days
      surface_time = increment_time(obs_time, sec_offset(inst_type),day_offset(inst_type))

      if ( surface_time >= time_start .and. surface_time <= time_end ) data_in_period = .true.
      if ( .not. data_in_period ) cycle

      call mpp_read(unit, field_sst, sfc_obs, tindex=k)
      do j=1, nlat
        do i=1, nlon
          lon = lons(i)
          lat = lats(j)
          data = MISSING_VALUE   ! snz add
          data_is_local = .false.

          if ( lon .lt. 0.0 ) lon = lon + 360.0
          if ( lon .gt. 360.0 ) lon = lon - 360.0
          if ( lon .gt. 60.0 ) lon = lon - 360.0

          if ( lat < obs_sbound(inst_type) .or. lat > obs_nbound(inst_type) ) cycle

          if ( localize_data ) then
            call kd_search_nnearest(kdroot, lon, lat, &
                    1, inds, dist, r_num, .false.)
            data_is_local = within_domain(lon1d(inds(1)), lat1d(inds(1)), &
                    isd+1, ied-1, jsd+1, jed-1, ni, nj)
          else
            data_is_local = .true.
          end if

          if (.not. data_is_local) cycle

          surface_count = surface_count + 1
          type_count(inst_type) = type_count(inst_type)+1

          data = sfc_obs(i,j)
          depth = 0.5
          flag = .true.

          if ( data .gt. 50 .or. data .lt. -5 ) then
            flag = .false.
          end if

          ! allocate profile structure content and put in data
          allocate(Prof%depth(1));Prof%depth=depth
          allocate(Prof%data(1));Prof%data=data
          allocate(Prof%flag(1));Prof%flag=flag
          Prof%variable = var_id
          Prof%inst_type = inst_type
          Prof%levels = num_levs
          Prof%lat = lat; Prof%lon = lon
          Prof%nbr_xi = lon1d(inds(1)); Prof%nbr_yi = lat1d(inds(1))
          Prof%nbr_dist = dist(1)
          Prof%time_window = time_window(inst_type)
          Prof%impact_levels = impact_levels(inst_type)
          Prof%temp_to_salt = temp_to_salt(inst_type)
          Prof%salt_to_temp = salt_to_temp(inst_type)
          Prof%obs_error = temp_error(inst_type)
          Prof%loc_dist = temp_dist(inst_type)
          Prof%time = surface_time

          if ( data .lt. sst_vimpact_temp ) Prof%impact_levels = sst_vimpact_levels

          if ( lat < lat_bound ) then ! calculate interpolation coefficients
             ri0 = frac_index(lon, T_grid%x(:,jsg))
             rj0 = frac_index(lat, T_grid%y(isg,:))
             i0 = floor(ri0)
             j0 = floor(rj0)
             if ( i0 > ieg .or. j0 > jeg ) then
                write (UNIT=emsg_local, FMT='("i0 = ",I8,", j0 = ",I8)') mpp_pe(), i0, j0
                call error_mesg('ocean_da_core_mod::open_profile_dataset',&
                     & 'For regular grids, either i0 > ieg or j0 > jeg.  '//trim(emsg_local), FATAL)
             end if
             Prof%i_index = ri0
             Prof%j_index = rj0
          else ! tripolar grids
             lon_out(1,1) = lon*DEG_TO_RAD
             lat_out(1,1) = lat*DEG_TO_RAD
             call horiz_interp_bilinear_new (Interp, T_grid%x*DEG_TO_RAD, T_grid%y*DEG_TO_RAD,&
                  & lon_out, lat_out, new_search=.true., no_crash_when_not_found=.true.)

             if ( Interp%wti(1,1,2) < 1.0 ) then
                i0 = Interp%i_lon(1,1,1)
             else
                i0 = Interp%i_lon(1,1,2)
             end if
             if ( Interp%wtj(1,1,2) < 1.0 ) then
                j0 = Interp%j_lat(1,1,1)
             else
                j0 = Interp%j_lat(1,1,2)
             end if
             if ( i0 > ieg .or. j0 > jeg ) then
                write (UNIT=emsg_local, FMT='("i0 = ",I6,", j0 = ",I6)') mpp_pe(), i0, j0
                call error_mesg('ocean_da_core_mod::open_profile_dataset',&
                     & 'For tripolar grids, either i0 > ieg or j0 > jeg', FATAL)
             end if
             if ( Interp%wti(1,1,2) < 1.0 ) then
                Prof%i_index =Interp%i_lon(1,1,1) + Interp%wti(1,1,2)
             else
                Prof%i_index =Interp%i_lon(1,1,2)
             end if
             if (Interp%wtj(1,1,2) < 1.0) then
                Prof%j_index =Interp%j_lat(1,1,1) + Interp%wtj(1,1,2)
             else
                Prof%j_index =Interp%j_lat(1,1,2)
             end if
          end if ! interpolation coefficients

          Prof%accepted = .true.

          if (i0 < 1 .or. j0 < 1) then
             Prof%accepted = .false.
          else
             Prof%basin_mask = T_grid%basin_mask(lon1d(inds(1)),lat1d(inds(1)))
          end if

          if ( Prof%accepted ) then ! check surface land-sea mask
             if ( i0 /= ieg .and. j0 /= jeg ) then
                if (T_grid%mask(i0,j0,1) == 0.0 .or.&
                     & T_grid%mask(i0+1,j0,1) == 0.0 .or.&
                     & T_grid%mask(i0,j0+1,1) == 0.0 .or.&
                     & T_grid%mask(i0+1,j0+1,1) == 0.0 ) then
                   Prof%accepted = .false.
                end if
             else if ( i0 == ieg .and. j0 /= jeg ) then
                if (T_grid%mask(i0,j0,1) == 0.0 .or.&
                     & T_grid%mask(1,j0,1) == 0.0 .or.&
                     & T_grid%mask(i0,j0+1,1) == 0.0 .or.&
                     & T_grid%mask(1,j0+1,1) == 0.0 ) then
                   Prof%accepted = .false.
                end if
             else if ( i0 /= ieg .and. j0 == jeg ) then
                if ( T_grid%mask(i0,j0,1) == 0.0 .or. T_grid%mask(i0+1,j0,1) == 0.0 ) then
                   Prof%accepted = .false.
                end if
             else
                if ( T_grid%mask(i0,j0,1) == 0.0 ) then
                   Prof%accepted = .false.
                end if
             end if
          end if ! check surface land-sea mask
   
          if ( Prof%accepted ) then ! determine vertical position and check mask at depth
             allocate(Prof%k_index(Prof%levels))
             do kk=1, Prof%levels
                Prof%k_index(kk) = 0.0
             end do
          end if ! determine vertical position and check mask at depth
   
          if ( Prof%accepted ) then ! calculate forward operator indices and weights
            allocate(Prof%obs_def(Prof%levels))
            ii = i0; jj = j0
            frac_lat = Prof%j_index - jj
            frac_lon = Prof%i_index - ii
   
            coef(1) = (1.0 - frac_lon) * (1.0 - frac_lat)
            coef(2) = frac_lon * (1.0 - frac_lat)
            coef(3) = (1.0 - frac_lon) * frac_lat
            coef(4) = frac_lon * frac_lat
   
            if ( ied > ni .and. ii < isd ) ii = ii + ni
            if ( isd < 1 .and. ii > ied ) ii = ii - ni
   
            do kk=1, Prof%levels
              k0 = floor(Prof%k_index(kk))
              frac_k = Prof%k_index(kk) - k0
   
              if ( k0 == 0 ) then
                state_index(1) = (jj-jsd)*lon_len + ii-isd + 1
                state_index(2) = (jj-jsd)*lon_len + ii-isd + 2
                state_index(3) = (jj-jsd+1)*lon_len + ii-isd + 1
                state_index(4) = (jj-jsd+1)*lon_len + ii-isd + 2
                state_index(5) = state_index(1)
                state_index(6) = state_index(2)
                state_index(7) = state_index(3)
                state_index(8) = state_index(4)
              else if (k0 == nk ) then
                state_index(1) = (k0-1)*blk + (jj-jsd)*lon_len+ii-isd+1
                state_index(2) = (k0-1)*blk + (jj-jsd)*lon_len+ii-isd+2
                state_index(3) = (k0-1)*blk + (jj-jsd+1)*lon_len+ii-isd+1
                state_index(4) = (k0-1)*blk + (jj-jsd+1)*lon_len+ii-isd+2
                state_index(5) = state_index(1)
                state_index(6) = state_index(2)
                state_index(7) = state_index(3)
                state_index(8) = state_index(4)
              else
                state_index(1) = (k0-1)*blk + (jj-jsd)*lon_len + ii-isd + 1
                state_index(2) = (k0-1)*blk + (jj-jsd)*lon_len + ii-isd + 2
                state_index(3) = (k0-1)*blk + (jj-jsd+1)*lon_len + ii-isd + 1
                state_index(4) = (k0-1)*blk + (jj-jsd+1)*lon_len + ii-isd + 2
                state_index(5) = k0*blk + (jj-jsd)*lon_len + ii-isd + 1
                state_index(6) = k0*blk + (jj-jsd)*lon_len + ii-isd + 2
                state_index(7) = k0*blk + (jj-jsd+1)*lon_len + ii-isd + 1
                state_index(8) = k0*blk + (jj-jsd+1)*lon_len + ii-isd + 2
              end if
   
              if ( frac_lon == 0.0 ) then
                state_index(2) = state_index(1)
                state_index(4) = state_index(3)
                state_index(6) = state_index(5)
                state_index(8) = state_index(7)
              end if
   
              if ( frac_lat == 0.0 ) then
                state_index(3) = state_index(1)
                state_index(4) = state_index(2)
                state_index(7) = state_index(5)
                state_index(8) = state_index(6)
              end if
   
              coef(5) = 1.0 - frac_k
              coef(6) = frac_k
   
              if ( frac_k == 0.0 ) then
                state_index(5) = state_index(1)
                state_index(6) = state_index(2)
                state_index(7) = state_index(3)
                state_index(8) = state_index(4)
              end if
   
              call def_forward_operator(8, state_index(1:8), coef(1:6), Prof%obs_def(kk))
            end do
          endif ! calculate forward operator indices and weights
   
          allocate(Prof%next) ! allocate next profile and link it to current one
          Prof%next%prev=>Prof
          Prof=>Prof%next
        end do
      end do
    end do

    call mpp_sync_self()
    call mpp_close(unit)
  end subroutine open_oisst_dataset

  ! get profiles obs relevant to current analysis interval
  subroutine get_profiles(model_time, Profiles, Current_profiles)
    type(time_type), intent(in) :: model_time
    type(ocean_profile_type), pointer :: Profiles
    type(ocean_profile_type), pointer :: Current_profiles

    type(ocean_profile_type), pointer :: Prof
    integer :: nprof
    integer :: i, yr, mon, day, hr, min, sec
    integer :: stdout_unit
    type(ocean_profile_type), pointer :: PREVIOUS=>NULL()
    !integer, dimension(10) :: current_type

    type(time_type) :: tdiff

    !current_type = 0
    nprof = 0
    stdout_unit = stdout()

    write (UNIT=stdout_unit, FMT='("Gathering profiles for current analysis time")')
    call get_date(model_time, yr, mon, day, hr, min, sec)
    write (UNIT=stdout_unit, FMT='("Current YYYY/MM/DD = ",I4,"/",I2,"/",I2,",",I2,":",I2)') yr, mon, day, hr, min

    Prof=>Profiles
    Current_profiles=>NULL()
    do while(associated(Prof))
       Prof%cnext=>NULL()
       Prof%cprev=>NULL()
       if ( Prof%time <= model_time ) then
          tdiff = model_time - Prof%time
       else
          tdiff = Prof%time - model_time
       end if
       ! no tdiff criteria for monthly mean data like
       ! but tdiff criteria has to be set for daily data
       if ( tdiff <= 2*Prof%time_window .and. Prof%accepted ) then
          nprof = nprof + 1
          !current_type(Prof%inst_type) = current_type(Prof%inst_type) + 1
          Prof%tdiff = tdiff
          if (.not.associated(Current_profiles)) then
              Current_profiles=>Prof
          else
              Prof%cprev=>PREVIOUS
              Prof%cprev%cnext=>Prof
          endif
          PREVIOUS=>Prof
      end if
      if (associated(Prof%next)) then
          Prof=>Prof%next
      else
          Prof=>NULL()
      endif
    end do

    return
  end subroutine get_profiles

  !=======================================================================
  ! Puts FO interpolation coefficients and state variable indices into an FO_type data structure.
  subroutine def_forward_operator(num_state, state_ind, coef, obs_def)
    integer, intent(in) :: num_state
    integer, dimension(num_state), intent(in) :: state_ind
    real, dimension(num_state-2), intent(in) :: coef
    type(forward_operator_type), intent(inout) :: obs_def

    integer :: i

    ! Set aside storage for defining this ob
    obs_def%num = num_state
    allocate(obs_def%state_var_index(num_state), obs_def%coef(num_state-2))

    ! Load the state variable index and coefficient for each state variable
    do i = 1, num_state
       obs_def%state_var_index(i) = state_ind(i)
    end do

    do i=1, num_state-2
       obs_def%coef(i) = coef(i)
    end do
  end subroutine def_forward_operator

end module ocean_da_core_mod

