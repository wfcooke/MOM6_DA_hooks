module write_ocean_obs_mod

use mpp_io_mod, only : fieldtype, axistype, mpp_open, mpp_write_meta, mpp_write, mpp_close
use mpp_io_mod, only : MPP_OVERWR, MPP_APPEND, MPP_WRONLY, MPP_NETCDF, MPP_MULTI, MPP_SINGLE
use mpp_mod, only : mpp_error, FATAL
use ocean_da_types_mod, only : MISSING_VALUE
use ocean_da_types_mod, only : ocean_profile_type, MAX_LEVELS_FILE
use time_manager_mod, only : time_type, get_time, set_date, operator ( - )

implicit none
private

type(fieldtype), save :: lon_field, lat_field, time_field, data_field, &
        obs_err_field, flag_field, forecast_field, analysis_field, depth_field, &
        lon_index_field, lat_index_field, link_field

integer, parameter :: ref_yr=1900, ref_mon=1, ref_day=1, &
        ref_hr=0, ref_min=0, ref_sec=0,max_files=1000

integer :: ref_seconds, ref_days, chid, wmo_id
integer, save :: station_num(max_files), unit_num(max_files), nfiles
type(time_type) :: ref_time, time
logical :: module_is_initialized=.false.

public :: open_profile_file, write_profile, close_profile_file, write_ocean_obs_init

#include <netcdf.inc>

contains

function open_profile_file(name, grid_lon, grid_lat,thread,fset)

  character(len=*), intent(in) :: name
  real, dimension(:), optional, intent(in) :: grid_lon, grid_lat
  integer, intent(in), optional :: thread, fset

  integer :: i, open_profile_file, unit
  integer :: threading, fileset
  character(len=128) :: units, time_units
  real, dimension(MAX_LEVELS_FILE) :: array

  type(axistype) :: depth_axis, station_axis, lon_axis, lat_axis
  
  threading=MPP_MULTI
  fileset=MPP_MULTI
  
  if (PRESENT(thread)) threading=thread
  if (PRESENT(fset)) fileset=fset
  
  ref_time = set_date(ref_yr, ref_mon, ref_day, ref_hr, ref_min, ref_sec)
  call get_time(ref_time, ref_seconds, ref_days)
  call mpp_open(unit, trim(name), action=MPP_OVERWR, form=MPP_NETCDF,&
                threading=threading, fileset=fileset)
  
  open_profile_file = unit
  
  nfiles=nfiles+1
  if (nfiles > max_files) call mpp_error(FATAL,'max number of profiles exceeded&
       &in module write_ocean_data, increase param : max_files')
  
  unit_num(nfiles) = unit
  
  if (PRESENT(grid_lon) .and. PRESENT(grid_lat)) then
     call mpp_write_meta(unit, lon_axis, 'grid_longitude','degrees_E',&
          'observational grid longitude',cartesian='X',sense=1,data=grid_lon)
  
     call mpp_write_meta(unit, lat_axis, 'grid_latitude','degrees_N',&
          'observational grid latitude', cartesian='Y',sense=1,data=grid_lat)
  endif
  
  !call mpp_write_meta(unit,depth_axis,'depth_index','none','depth index',&
  !                  cartesian='Z',sense=-1)!,data=(/(float(i),i=1,MAX_LEVELS_FILE)/))
  !pgf90 complains about the above. This is a compiler bug. Workaround:
  array = (/(float(i),i=1,MAX_LEVELS_FILE)/)
  call mpp_write_meta(unit,depth_axis,'depth_index','none','depth index',&
                      cartesian='Z',sense=-1,data=array)
  
  call mpp_write_meta(unit,station_axis,'station_index','none',&
                      'station index', cartesian='T',sense=1)
  
  call mpp_write_meta(unit,lon_field,(/station_axis/),&
          'longitude','degrees_E','longitude',min=-301.0,max=61.0)
  call mpp_write_meta(unit,lat_field,(/station_axis/),&
          'latitude','degrees_N','latitude',min=-91.0,max=91.0)
  
  write(time_units,'(a,i4.4,a,i2.2,a,i2.2,a)')  'days since ',ref_yr,'-',ref_mon,'-',ref_day,' 00:00:00'
  
  call mpp_write_meta(unit,time_field,(/station_axis/),'time',trim(time_units),'time')
  
  units='none'
  call mpp_write_meta(unit,obs_err_field,(/station_axis/),&
          'obs_error',trim(units),'assimilated observation error',missing=MISSING_VALUE)
  
  call mpp_write_meta(unit,link_field,(/station_axis/),&
          'link','none','partial_profile flag')
  
  call mpp_write_meta(unit,data_field,(/depth_axis,station_axis/),&
          'data',trim(units),'in-situ observation',min=-10.0,max=50.0,missing=MISSING_VALUE)
 
  call mpp_write_meta(unit,forecast_field,(/depth_axis,station_axis/),&
          'obs_prior',trim(units),'first-guess in-situ values',min=-10.0,max=50.0,missing=MISSING_VALUE)
 
  call mpp_write_meta(unit,analysis_field,(/depth_axis,station_axis/),&
          'obs_posterior',trim(units),'analysis in-situ values',min=-10.0,max=50.0,missing=MISSING_VALUE)
 
  call mpp_write_meta(unit,depth_field,(/depth_axis,station_axis/),&
          'depth','meters','depth of obs',min=0.0,max=7000.0,missing=MISSING_VALUE)
  
  call mpp_write_meta(unit,flag_field,(/depth_axis,station_axis/),&  
          'obs_flag','none','each level flag',missing=MISSING_VALUE)
  
  call mpp_write(unit, depth_axis)
  
  if (PRESENT(grid_lon).and.PRESENT(grid_lat)) then
     call mpp_write(unit, lon_axis)
     call mpp_write(unit, lat_axis)
  endif

end function open_profile_file

subroutine write_profile(unit,profile)
  integer, intent(in) :: unit
  type(ocean_profile_type), intent(in) :: profile
  
  real, dimension(MAX_LEVELS_FILE) :: obs_data, depth
  real, dimension(MAX_LEVELS_FILE) :: forecast, analysis
  integer :: levels, secs, days, i, j, nlinks
  real :: days_since, station
  real, dimension(MAX_LEVELS_FILE) :: flag
  integer :: findex
  logical :: debug=.false.
  logical :: f_avail=.false. , a_avail=.false.
  ! find file index from file unit list
  
  findex=-1
  do i=1,nfiles
     if (unit_num(i) .eq. unit) then
         findex=i
         exit
     endif
  enddo
  
  if (findex .eq. -1) call mpp_error(FATAL,'Attempt write to unopened file in&
       &write_ocean_data_mod:write_profile_data')
  
  station_num(findex)=station_num(findex)+1
  station=station_num(findex)
  
  if(associated(profile%forecast)) f_avail = .true.
  if(associated(profile%analysis)) a_avail = .true.

  levels = min(profile%levels,MAX_LEVELS_FILE)
  obs_data=MISSING_VALUE; depth=MISSING_VALUE
  forecast=MISSING_VALUE; analysis=MISSING_VALUE
  flag=MISSING_VALUE
  
  obs_data(1:levels)=profile%data(1:levels)
  if(f_avail) forecast(1:levels)=profile%forecast(1:levels)
  if(a_avail) analysis(1:levels)=profile%analysis(1:levels)
  flag(1:levels)=profile%flag(1:levels)
  depth(1:levels)=profile%depth(1:levels)
  time = profile%time - ref_time
  call get_time(time, secs, days)
  days_since = days + secs/86400.

  call mpp_write(unit,data_field,obs_data,station)
  if(f_avail) call mpp_write(unit,forecast_field,forecast,station)
  if(a_avail) call mpp_write(unit,analysis_field,analysis,station)
  call mpp_write(unit,depth_field,depth,station)
  call mpp_write(unit,lon_field,profile%lon,station)
  call mpp_write(unit,lat_field,profile%lat,station)
  call mpp_write(unit,time_field,days_since,station)
  call mpp_write(unit,obs_err_field,profile%obs_error,station)
  nlinks = 0
  if (profile%levels .gt. MAX_LEVELS_FILE) then
      nlinks = ceiling(float(profile%levels)/float(MAX_LEVELS_FILE)) - 1
  endif
 
  if (nlinks .gt. 0) then
      call mpp_write(unit,link_field,1.,station)
  else
      call mpp_write(unit,link_field,0.,station)
  endif
 
  do i = 1, nlinks
     station_num(findex)=station_num(findex)+1
     station=station_num(findex)
     if (i.eq.nlinks) then
         levels = mod(profile%levels,MAX_LEVELS_FILE)
         if (levels .eq. 0) levels = MAX_LEVELS_FILE
     else
         levels = MAX_LEVELS_FILE
     endif
     obs_data = MISSING_VALUE; depth = MISSING_VALUE
     forecast=MISSING_VALUE; analysis=MISSING_VALUE
     flag = MISSING_VALUE
  
     obs_data(1:levels)=profile%data((MAX_LEVELS_FILE*i)+1:(MAX_LEVELS_FILE*i)+levels)
     if(f_avail) forecast(1:levels)=profile%forecast((MAX_LEVELS_FILE*i)+1:(MAX_LEVELS_FILE*i)+levels)
     if(a_avail) analysis(1:levels)=profile%analysis((MAX_LEVELS_FILE*i)+1:(MAX_LEVELS_FILE*i)+levels)
     flag(1:levels)=profile%flag((MAX_LEVELS_FILE*i)+1:(MAX_LEVELS_FILE*i)+levels)
     depth(1:levels)=profile%depth((MAX_LEVELS_FILE*i)+1:(MAX_LEVELS_FILE*i)+levels)
     call mpp_write(unit,data_field,obs_data,station)
     if(f_avail) call mpp_write(unit,forecast_field,forecast,station)
     if(a_avail) call mpp_write(unit,analysis_field,analysis,station)
     call mpp_write(unit,depth_field,depth,station)
     call mpp_write(unit,lon_field,profile%lon,station)
     call mpp_write(unit,lat_field,profile%lat,station)
     call mpp_write(unit,time_field,days_since,station)
     call mpp_write(unit,obs_err_field,profile%obs_error,station)
  
     if (i .lt. nlinks) then
         call mpp_write(unit,link_field,1.,station)
     else
         call mpp_write(unit,link_field,0.,station)
     endif
  
  enddo

end subroutine write_profile

subroutine close_profile_file(unit)

  integer, intent(in) :: unit

  call mpp_close(unit)

end subroutine close_profile_file

subroutine write_ocean_obs_init()

  module_is_initialized=.true.

  station_num=0;unit_num=0;nfiles=0

  return

end subroutine write_ocean_obs_init

end module write_ocean_obs_mod
