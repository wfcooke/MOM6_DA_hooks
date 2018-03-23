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
module kdtree


  implicit none
  private
  public :: kd_root

  type kd_root
     integer :: blank !< just a dumb type
  end type kd_root


end module kdtree
