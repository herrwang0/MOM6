! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> A template of a user to code up customized initial conditions.
module user_initialization

use MOM_error_handler, only : MOM_mesg, MOM_error, FATAL, is_root_pe
use MOM_dyn_horgrid, only : dyn_horgrid_type
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_get_input, only : directories
use MOM_grid, only : ocean_grid_type
use MOM_open_boundary, only : ocean_OBC_type, OBC_NONE, OBC_segment_type
use MOM_open_boundary, only : OBC_DIRECTION_E, OBC_DIRECTION_W, OBC_DIRECTION_N
use MOM_open_boundary, only : OBC_DIRECTION_S
use MOM_sponge, only : set_up_sponge_field, initialize_sponge, sponge_CS
use MOM_tracer_registry, only : tracer_registry_type
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : thermo_var_ptrs
use MOM_verticalGrid, only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

public USER_set_coord, USER_initialize_topography, USER_initialize_thickness
public USER_initialize_velocity, USER_init_temperature_salinity
public USER_initialize_sponges, USER_set_OBC_data, USER_set_rotation

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> A module variable that should not be used.
!! \todo Move this module variable into a control structure.
logical :: first_call = .true.

contains

!> Set vertical coordinates.
subroutine USER_set_coord(Rlay, g_prime, GV, US, param_file)
  type(verticalGrid_type),  intent(in)  :: GV      !< The ocean's vertical grid structure
  real, dimension(GV%ke),   intent(out) :: Rlay    !< Layer potential density [R ~> kg m-3].
  real, dimension(GV%ke+1), intent(out) :: g_prime !< The reduced gravity at each
                                                   !! interface [L2 Z-1 T-2 ~> m s-2].
  type(unit_scale_type),    intent(in)  :: US      !< A dimensional unit scaling type
  type(param_file_type),    intent(in)  :: param_file !< A structure indicating the
                                                   !! open file to parse for model
                                                   !! parameter values.

  call MOM_error(FATAL, &
    "USER_initialization.F90, USER_set_coord: " // &
    "Unmodified user routine called - you must edit the routine to use it")
  Rlay(:) = 0.0
  g_prime(:) = 0.0

  if (first_call) call write_user_log(param_file)

end subroutine USER_set_coord

!> Initialize topography.
subroutine USER_initialize_topography(D, G, param_file, max_depth, US)
  type(dyn_horgrid_type),          intent(in)  :: G !< The dynamic horizontal grid type
  real, dimension(G%isd:G%ied,G%jsd:G%jed), &
                                   intent(out) :: D !< Ocean bottom depth [Z ~> m]
  type(param_file_type),           intent(in)  :: param_file !< Parameter file structure
  real,                            intent(in)  :: max_depth !< Maximum model depth [Z ~> m]
  type(unit_scale_type),           intent(in)  :: US !< A dimensional unit scaling type

  call MOM_error(FATAL, &
    "USER_initialization.F90, USER_initialize_topography: " // &
    "Unmodified user routine called - you must edit the routine to use it")

  D(:,:) = 0.0

  if (first_call) call write_user_log(param_file)

end subroutine USER_initialize_topography

!> Initialize thicknesses in depth units.  These will be converted to thickness units later.
subroutine USER_initialize_thickness(h, G, GV, param_file, just_read)
  type(ocean_grid_type),   intent(in)  :: G  !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)  :: GV !< The ocean's vertical grid structure.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(out) :: h  !< The thicknesses being initialized [Z ~> m]
  type(param_file_type),   intent(in)  :: param_file !< A structure indicating the open
                                             !! file to parse for model parameter values.
  logical,                 intent(in)  :: just_read !< If true, this call will
                                             !! only read parameters without changing h.

  call MOM_error(FATAL, &
    "USER_initialization.F90, USER_initialize_thickness: " // &
    "Unmodified user routine called - you must edit the routine to use it")

  if (just_read) return ! All run-time parameters have been read, so return.

  h(:,:,1:GV%ke) = 0.0 ! h should be set in [Z ~> m].  It will be converted to thickness units
                       ! [H ~> m or kg m-2] once the temperatures and salinities are known.

  if (first_call) call write_user_log(param_file)

end subroutine USER_initialize_thickness

!> initialize velocities.
subroutine USER_initialize_velocity(u, v, G, GV, US, param_file, just_read)
  type(ocean_grid_type),                       intent(in)  :: G !< Ocean grid structure.
  type(verticalGrid_type),                     intent(in)  :: GV !< The ocean's vertical grid structure.
  real, dimension(SZIB_(G), SZJ_(G),SZK_(GV)), intent(out) :: u !< i-component of velocity [L T-1 ~> m s-1]
  real, dimension(SZI_(G), SZJB_(G),SZK_(GV)), intent(out) :: v !< j-component of velocity [L T-1 ~> m s-1]
  type(unit_scale_type),                       intent(in)  :: US !< A dimensional unit scaling type
  type(param_file_type),                       intent(in)  :: param_file !< A structure indicating the
                                                            !! open file to parse for model
                                                            !! parameter values.
  logical,                                     intent(in)  :: just_read !< If true, this call will
                                                      !! only read parameters without changing u & v.

  call MOM_error(FATAL, &
    "USER_initialization.F90, USER_initialize_velocity: " // &
    "Unmodified user routine called - you must edit the routine to use it")

  if (just_read) return ! All run-time parameters have been read, so return.

  u(:,:,1) = 0.0
  v(:,:,1) = 0.0

  if (first_call) call write_user_log(param_file)

end subroutine USER_initialize_velocity

!> This function puts the initial layer temperatures and salinities
!! into T(:,:,:) and S(:,:,:).
subroutine USER_init_temperature_salinity(T, S, G, GV, param_file, just_read)
  type(ocean_grid_type),                     intent(in)  :: G !< Ocean grid structure.
  type(verticalGrid_type),                   intent(in)  :: GV !< The ocean's vertical grid structure.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(out) :: T !< Potential temperature [C ~> degC].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(out) :: S !< Salinity [S ~> ppt].
  type(param_file_type),                     intent(in)  :: param_file !< A structure indicating the
                                                            !! open file to parse for model
                                                            !! parameter values.
  logical,                                   intent(in)  :: just_read !< If true, this call will only
                                                           !! read parameters without changing T & S.

  call MOM_error(FATAL, &
    "USER_initialization.F90, USER_init_temperature_salinity: " // &
    "Unmodified user routine called - you must edit the routine to use it")

  if (just_read) return ! All run-time parameters have been read, so return.

  T(:,:,1) = 0.0
  S(:,:,1) = 0.0

  if (first_call) call write_user_log(param_file)

end subroutine USER_init_temperature_salinity

!> Set up the sponges.
subroutine USER_initialize_sponges(G, GV, use_temp, tv, param_file, CSp, h)
  type(ocean_grid_type),   intent(in) :: G             !< Ocean grid structure.
  type(verticalGrid_type), intent(in) :: GV            !< The ocean's vertical grid structure.
  logical,                 intent(in) :: use_temp      !< If true, temperature and salinity are state variables.
  type(thermo_var_ptrs),   intent(in) :: tv            !< A structure containing pointers
                                                       !! to any available thermodynamic
                                                       !! fields, potential temperature and
                                                       !! salinity or mixed layer density.
                                                       !! Absent fields have NULL ptrs.
  type(param_file_type),   intent(in) :: param_file    !< A structure indicating the
                                                       !! open file to parse for model
                                                       !! parameter values.
  type(sponge_CS),         pointer    :: CSp           !< A pointer to the sponge control structure.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in) :: h             !< Layer thicknesses [H ~> m or kg m-2].
  call MOM_error(FATAL, &
    "USER_initialization.F90, USER_initialize_sponges: " // &
    "Unmodified user routine called - you must edit the routine to use it")

  if (first_call) call write_user_log(param_file)

end subroutine USER_initialize_sponges

!> This subroutine sets the properties of flow at open boundary conditions.
subroutine USER_set_OBC_data(OBC, tv, G, GV, US, param_file, tr_Reg)
  type(ocean_OBC_type),       pointer    :: OBC   !< This open boundary condition type specifies
                                                  !! whether, where, and what open boundary
                                                  !! conditions are used.
  type(thermo_var_ptrs),      intent(in) :: tv    !< A structure containing pointers to any
                                       !! available thermodynamic fields, including potential
                                       !! temperature and salinity or mixed layer density. Absent
                                       !! fields have NULL ptrs.
  type(ocean_grid_type),      intent(in) :: G     !< The ocean's grid structure.
  type(verticalGrid_type),    intent(in) :: GV    !< The ocean's vertical grid structure.
  type(unit_scale_type),      intent(in) :: US    !< A dimensional unit scaling type.
  type(param_file_type),      intent(in) :: param_file !< A structure indicating the
                                                  !! open file to parse for model
                                                  !! parameter values.
  type(tracer_registry_type), pointer    :: tr_Reg !< Tracer registry.
!  call MOM_error(FATAL, &
!   "USER_initialization.F90, USER_set_OBC_data: " // &
!   "Unmodified user routine called - you must edit the routine to use it")

  ! Local variables
  type(OBC_segment_type), pointer :: segment => NULL() ! A pointer to an open boundary segment
  real, dimension(SZK_(GV)) :: trans_dflt ! The layer transports normal to any open boundary segment
                                          ! that does not specify its own values, in mks units [m3 s-1]
  real, dimension(SZK_(GV)) :: vel_dflt   ! The layer velocities normal to any open boundary segment
                                          ! that does not specify its own values, in mks units [m s-1]
  real, dimension(SZK_(GV)) :: trans      ! The layer transports normal to one open boundary
                                          ! segment [H L2 T-1 ~> m3 s-1 or kg s-1]
  real, dimension(SZK_(GV)) :: vel        ! The layer velocities normal to one open boundary
                                          ! segment [L T-1 ~> m s-1]
  real :: ssh        ! The sea surface height imposed outside of the outflow segment [Z ~> m]
  character(len=48) :: seg_param_str ! The name of a segment-specific run-time parameter
  character(len=40) :: mdl = "user_initialization" ! This module's name.
  integer :: ssh_seg ! The number of the segment where the sea surface height is imposed, or 0
                     ! if there is no such segment
  integer :: is, ie, js, je, i, j, k, n, nz

  nz = GV%ke

  ! The first two parameters set the values that are used at any segment that does not specify
  ! its own values with the segment-specific parameters that are read below.
  call get_param(param_file, mdl, "OBC_USER_TRANSPORT", trans_dflt, &
                 "The transports normal to the open boundary segments, from the top layer to the "//&
                 "bottom layer.  This is used at any segment that does not set its own "//&
                 "OBC_SEGMENT_%%%_USER_TRANSPORT.", &
                 units="m3 s-1", default=0.0)
  call get_param(param_file, mdl, "OBC_USER_VELOCITY", vel_dflt, &
                 "The velocities normal to the open boundary segments, from the top layer to the "//&
                 "bottom layer.  This is used at any segment that does not set its own "//&
                 "OBC_SEGMENT_%%%_USER_VELOCITY.", &
                 units="m s-1", default=0.0)
  call get_param(param_file, mdl, "OBC_USER_FLATHER_SSH", ssh, &
                 "The sea surface height imposed outside of the outflow segment, taken to be the "//&
                 "last segment without a specified transport or velocity.  This is only used "//&
                 "with Flather open boundary conditions.", &
                 units="m", default=0.0, scale=US%m_to_Z)

  ssh_seg = 0

  do n=1,OBC%number_of_segments
    ! The segment-specific parameters are read on all PEs, including those where this segment is
    ! not present, so that they are logged consistently.
    write(seg_param_str, "('OBC_SEGMENT_',i3.3,'_USER_TRANSPORT')") n
    call get_param(param_file, mdl, trim(seg_param_str), trans, &
                   "The transports normal to this open boundary segment, from the top layer to "//&
                   "the bottom layer.", &
                   units="m3 s-1", defaults=trans_dflt, scale=(GV%m_to_H*US%m_to_L**2)*US%T_to_s)
    write(seg_param_str, "('OBC_SEGMENT_',i3.3,'_USER_VELOCITY')") n
    call get_param(param_file, mdl, trim(seg_param_str), vel, &
                   "The velocities normal to this open boundary segment, from the top layer to "//&
                   "the bottom layer.", &
                   units="m s-1", defaults=vel_dflt, scale=US%m_s_to_L_T)

    ! The sea surface height is imposed at the last segment where no transport or velocity has
    ! been specified, which is taken to be the outflow segment.
    if (all(trans(:) == 0.0) .and. all(vel(:) == 0.0)) ssh_seg = n

    segment => OBC%segment(n)
    if (.not. segment%on_pe) cycle

    if (segment%is_E_or_W) then
      is = segment%HI%IsdB ; ie = segment%HI%IedB
      js = segment%HI%jsd ; je = segment%HI%jed
      do k=1,nz ; do j=js,je ; do I=is,ie ; if (abs(OBC%segnum_u(I,j)) == n) then
        segment%normal_trans(I,j,k) = trans(k)
        segment%normal_vel(I,j,k) = vel(k)
      endif ; enddo ; enddo ; enddo
    else
      is = segment%HI%isd ; ie = segment%HI%ied
      js = segment%HI%JsdB ; je = segment%HI%JedB
      do k=1,nz ; do J=js,je ; do i=is,ie ; if (abs(OBC%segnum_v(i,J)) == n) then
        segment%normal_trans(i,J,k) = trans(k)
        segment%normal_vel(i,J,k) = vel(k)
      endif ; enddo ; enddo ; enddo
    endif
  enddo ! n-loop over the open boundary segments

  ! The barotropic solver reads the Flather sea surface height everywhere within the segment's
  ! index range, so it is set here without excluding the land points along the segment.
  if (ssh_seg > 0) then
    segment => OBC%segment(ssh_seg)
    if (segment%on_pe) then
      if (allocated(segment%SSH)) segment%SSH(:,:) = ssh
    endif
  endif

  if (first_call) call write_user_log(param_file)

end subroutine USER_set_OBC_data

subroutine USER_set_rotation(G, param_file)
  type(ocean_grid_type), intent(inout) :: G    !< The ocean's grid structure
  type(param_file_type), intent(in)    :: param_file !< A structure to parse for run-time parameters
  call MOM_error(FATAL, &
    "USER_initialization.F90, USER_set_rotation: " // &
    "Unmodified user routine called - you must edit the routine to use it")

  if (first_call) call write_user_log(param_file)

end subroutine USER_set_rotation

!> Write output about the parameter values being used.
subroutine write_user_log(param_file)
  type(param_file_type), intent(in) :: param_file !< A structure indicating the
                                                  !! open file to parse for model
                                                  !! parameter values.

  ! This include declares and sets the variable "version".
# include "version_variable.h"
  character(len=40)  :: mdl = "user_initialization" ! This module's name.

  call log_version(param_file, mdl, version)
  first_call = .false.

end subroutine write_user_log

!> \namespace user_initialization
!!
!!  This subroutine initializes the fields for the simulations.
!!  The one argument passed to initialize, Time, is set to the
!!  current time of the simulation.  The fields which might be initialized
!!  here are:
!!  - u - Zonal velocity [Z T-1 ~> m s-1].
!!  - v - Meridional velocity [Z T-1 ~> m s-1].
!!  - h - Layer thickness [H ~> m or kg m-2].  (Must be positive.)
!!  - G%bathyT - Basin depth [Z ~> m].
!!  - G%CoriolisBu - The Coriolis parameter [T-1 ~> s-1].
!!  - GV%g_prime - The reduced gravity at each interface [L2 Z-1 T-2 ~> m s-2].
!!  - GV%Rlay - Layer potential density (coordinate variable) [R ~> kg m-3].
!!  If ENABLE_THERMODYNAMICS is defined:
!!  - T - Temperature [C ~> degC].
!!  - S - Salinity [S ~> ppt].
!!  If BULKMIXEDLAYER is defined:
!!  - Rml - Mixed layer and buffer layer potential densities [R ~> kg m-3].
!!  If SPONGE is defined:
!!  - A series of subroutine calls are made to set up the damping
!!    rates and reference profiles for all variables that are damped
!!    in the sponge.
!!
!!  Any user provided tracer code is also first linked through this
!!  subroutine.
!!
!!  These variables are all set in the set of subroutines (in this
!!  file) USER_initialize_bottom_depth, USER_initialize_thickness,
!!  USER_initialize_velocity,  USER_initialize_temperature_salinity,
!!  USER_initialize_mixed_layer_density, USER_initialize_sponges,
!!  USER_set_coord, and USER_set_ref_profile.
!!
!!  The names of these subroutines should be self-explanatory. They
!!  start with "USER_" to indicate that they will likely have to be
!!  modified for each simulation to set the initial conditions and
!!  boundary conditions.  Most of these take two arguments: an integer
!!  argument specifying whether the fields are to be calculated
!!  internally or read from a NetCDF file; and a string giving the
!!  path to that file.  If the field is initialized internally, the
!!  path is ignored.

end module user_initialization
