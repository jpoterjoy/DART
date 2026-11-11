! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

module cov_cutoff_mod


use     types_mod, only : r8
use utilities_mod, only : error_handler, E_ERR, E_MSG, &
                          do_output, do_nml_file, do_nml_term, nmlfileunit, &
                          find_namelist_in_file, check_namelist_read
use location_mod,  only : location_type, query_location, is_vertical
use obs_kind_mod,  only : get_index_for_type_of_obs, get_num_types_of_obs

implicit none
private

public :: comp_cov_factor, vloc_has_table

character(len=*), parameter :: source = 'cov_cutoff_mod.f90'


!============================================================================

!---- namelist with default values
logical :: namelist_initialized = .false.

integer :: select_localization = 1
! Value 1 selects default Gaspari-Cohn cutoff
! Value 2 selects boxcar
! Value 3 selects ramped boxcar

! Vertical localization tables for obs types without a vertical coordinate
! (e.g. satellite-like layer observations; models with levels 1..vloc_nlev):
! for an ob of a type in vloc_types located with VERTISUNDEF, the factor is
! multiplied by vloc_factors(k, t) at the target's model level k (linear in
! between; vloc_surface_level for VERTISSURFACE targets; 1 for targets
! without a vertical). Types not listed, or located with a vertical, are
! unchanged.
integer, parameter :: max_vloc_types = 20, max_vloc_lev = 100
character(len=64)  :: vloc_types(max_vloc_types) = ''
real(r8)           :: vloc_factors(max_vloc_lev, max_vloc_types) = 1.0_r8
integer            :: vloc_nlev = 0
integer            :: vloc_surface_level = 1
integer, allocatable :: vloc_index(:)          ! obs type -> table (0: none)

namelist / cov_cutoff_nml / select_localization, vloc_types, vloc_factors, vloc_nlev, vloc_surface_level

!============================================================================

contains

!======================================================================



function comp_cov_factor(z_in, c, obs_loc, obs_type, target_loc, target_kind, &
   localization_override)
!----------------------------------------------------------------------
! function comp_cov_factor(z_in, c)
!
! Computes a covariance cutoff function from Gaspari and Cohn
! QJRMS, 125, 723-757.  (their eqn. 4.10)
!
! z_in is the distance while c is the cutoff distance. 
! For distances greater than 2c, the cov_factor returned goes to 0.

! Other ramping shapes are also available and can be selected by a namelist
! parameter. At present, these include a boxcar with the given halfwidth
! and a ramp in which the weight is set at 1.0 up to the half-width 
! distance and then decreases linearly to 0 at twice the half-width 
! distance.

! Additional information is passed in about the location and specific type of the
! observation and the location and generic kind of the variable being targeted for
! increments. These can be used for more refined algorithms that want to 
! make the cutoff a function of these additional arguments. 

implicit none

real(r8),                      intent(in) :: z_in, c
type(location_type), optional, intent(in) :: obs_loc, target_loc
integer,             optional, intent(in) :: obs_type, target_kind
integer,             optional, intent(in) :: localization_override
real(r8)                                  :: comp_cov_factor

real(r8) :: z, r
integer  :: iunit, io
integer  :: localization_type

!--------------------------------------------------------
! Initialize namelist if not already done
if(.not. namelist_initialized) then


   namelist_initialized = .true.

   ! Read the namelist entry
   call find_namelist_in_file("input.nml", "cov_cutoff_nml", iunit)
   read(iunit, nml = cov_cutoff_nml, iostat = io)
   call check_namelist_read(iunit, io, "cov_cutoff_nml")

   if (do_nml_file()) write(nmlfileunit,nml=cov_cutoff_nml)
   if (do_nml_term()) write(     *     ,nml=cov_cutoff_nml)


   call init_vloc()

   if (do_output()) then
      select case (select_localization)
         case (1)
            call error_handler(E_MSG,'comp_cov_factor:', &
               'Standard Gaspari Cohn localization selected')
         case (2)
            call error_handler(E_MSG,'comp_cov_factor:', &
               'Boxcar localization selected')
         case (3)
            call error_handler(E_MSG,'comp_cov_factor:', &
               'Ramped localization selected')
         case default
            call error_handler(E_ERR,'comp_cov_factor', &
               'Illegal value of "select_localization" in cov_cutoff_mod namelist', source)
      end select
   endif

endif
!---------------------------------------------------------

if(present(localization_override)) then
   localization_type = localization_override
else
   localization_type = select_localization
endif

z = abs(z_in)

!----------------------------------------------------------

if(localization_type == 1) then ! Standard Gaspari Cohn localization

   if( z >= c*2.0_r8 ) then

      comp_cov_factor = 0.0_r8

   else if( z <= c ) then
      r = z / c
      comp_cov_factor = &
           ( ( ( -0.25_r8*r +0.5_r8 )*r +0.625_r8 )*r -5.0_r8/3.0_r8 )*r**2 + 1.0_r8
!!$           r**5 * (-0.25_r8 ) + &
!!$           r**4 / 2.0_r8 +              &
!!$           r**3 * 5.0_r8/8.0_r8 -       &
!!$           r**2 * 5.0_r8/3.0_r8 + 1.0_r8
   else

      r = z / c
      comp_cov_factor = &
           ( ( ( ( r/12.0_r8 -0.5_r8 )*r +0.625_r8 )*r +5.0_r8/3.0_r8 )*r -5.0_r8 )*r &
!!$           r**5 / 12.0_r8  -  &
!!$           r**4 / 2.0_r8   +  &
!!$           r**3 * 5.0_r8 / 8.0_r8 + &
!!$           r**2 * 5.0_r8 / 3.0_r8 - 5.0_r8*r &
           + 4.0_r8 - 2.0_r8 / (3.0_r8 * r) 
   endif

else if(localization_type == 2) then ! BOXCAR localization

   if(z < 2.0_r8 * c) then
      comp_cov_factor = 1.0_r8
   else
      comp_cov_factor = 0.0_r8
   endif

else if(localization_type == 3) then ! Ramped localization

   if(z >= 2.0_r8 * c) then
      comp_cov_factor = 0.0_r8
   else if(z >= c .and. z < 2.0_r8 * c) then
      comp_cov_factor = (2.0_r8 * c - z) / c
   else
      comp_cov_factor = 1.0_r8
   endif

else ! Otherwise namelist parameter is illegal; this is an error

     call error_handler(E_ERR,'comp_cov_factor', &
              'Illegal value of "localization" in cov_cutoff_mod namelist', source)

endif

! vertical table of an ob type without a vertical coordinate
if (present(obs_type) .and. present(obs_loc) .and. present(target_loc)) then
   if (vloc_type_of(obs_type) > 0) then
      if (is_vertical(obs_loc, 'UNDEFINED')) &
         comp_cov_factor = comp_cov_factor*vloc_factor(vloc_type_of(obs_type), target_loc)
   endif
endif

end function comp_cov_factor

!----------------------------------------------------------------------
!> Map the namelist's vloc_types to obs type indices (once)

subroutine init_vloc()

integer :: t, i

allocate(vloc_index(max(1, get_num_types_of_obs())))
vloc_index = 0
do t = 1, max_vloc_types
   if (len_trim(vloc_types(t)) == 0) cycle
   if (vloc_nlev < 1 .or. vloc_nlev > max_vloc_lev) call error_handler(E_ERR, 'comp_cov_factor', &
      'vloc_types needs 1 <= vloc_nlev <= 100', source)
   i = get_index_for_type_of_obs(vloc_types(t))
   if (i < 1) call error_handler(E_ERR, 'comp_cov_factor', &
      'unknown obs type in vloc_types: '//trim(vloc_types(t)), source)
   vloc_index(i) = t
end do

end subroutine init_vloc

!----------------------------------------------------------------------
!> Table of obs type itype (0: none)

function vloc_type_of(itype)

integer, intent(in) :: itype
integer             :: vloc_type_of

vloc_type_of = 0
if (.not. allocated(vloc_index)) return
if (itype < 1 .or. itype > size(vloc_index)) return
vloc_type_of = vloc_index(itype)

end function vloc_type_of

!----------------------------------------------------------------------
!> True if obs type itype has a vertical table (for the LPF's obs bins)

function vloc_has_table(itype)

integer, intent(in) :: itype
logical             :: vloc_has_table
real(r8)            :: dummy

if (.not. namelist_initialized) dummy = comp_cov_factor(0.0_r8, 1.0_r8)
vloc_has_table = vloc_type_of(itype) > 0

end function vloc_has_table

!----------------------------------------------------------------------
!> Vertical factor of table t at the target's vertical

function vloc_factor(t, target_loc)

integer,             intent(in) :: t
type(location_type), intent(in) :: target_loc
real(r8)                        :: vloc_factor

real(r8) :: lv, f
integer  :: k

vloc_factor = 1.0_r8
if (is_vertical(target_loc, 'LEVEL')) then
   lv = min(max(query_location(target_loc, 'VLOC'), 1.0_r8), real(vloc_nlev, r8))
   k  = min(int(lv), vloc_nlev - 1)
   if (vloc_nlev == 1) then
      vloc_factor = vloc_factors(1, t)
   else
      f = lv - k
      vloc_factor = (1.0_r8 - f)*vloc_factors(k, t) + f*vloc_factors(k + 1, t)
   endif
else if (is_vertical(target_loc, 'SURFACE')) then
   vloc_factor = vloc_factors(vloc_surface_level, t)
endif

end function vloc_factor

end module cov_cutoff_mod

