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
! (e.g. satellite channels, layer observations): for an ob of a type in
! vloc_types located with VERTISUNDEF, the factor is multiplied by the type's
! table, vloc_factors(1:vloc_nlev, t), at the target's vertical:
!   vloc_coord(t) = 'level' (default): row k is model level k (targets located
!                   by level; linear in between);
!                   'pressure', 'height', 'scale_height': row k is at the
!                   vertical coordinate vloc_points(k, t) (Pa, m, or scale
!                   height; strictly monotonic), targets located in that
!                   coordinate (linear in log(p) for pressure, linear in the
!                   coordinate otherwise; beyond the ends, the end values).
! Targets at the surface (VERTISSURFACE) use row vloc_surface_level; targets
! without a vertical (VERTISUNDEF) get 1. A target located in any other
! vertical than the table's is an error (vloc_mismatch = 'error', default) or
! gets 1 with a warning once per table and vertical ('warn'). Types not
! listed, or located with a vertical, are unchanged.
integer, parameter :: max_vloc_types = 20, max_vloc_lev = 200
character(len=64)  :: vloc_types(max_vloc_types) = ''
real(r8)           :: vloc_factors(max_vloc_lev, max_vloc_types) = 1.0_r8
character(len=16)  :: vloc_coord(max_vloc_types) = 'level'
real(r8)           :: vloc_points(max_vloc_lev, max_vloc_types) = -huge(1.0_r8)
integer            :: vloc_nlev = 0
integer            :: vloc_surface_level = 1
character(len=8)   :: vloc_mismatch = 'error'
integer, allocatable :: vloc_index(:)          ! obs type -> table (0: none)
integer            :: vloc_ic(max_vloc_types) = 1                ! table coordinate: 1 level, 2 pressure, 3 height, 4 scale height
logical            :: vloc_warned(max_vloc_types, -2:4) = .false. ! vloc_mismatch = 'warn': warned for (table, target vertical)

namelist / cov_cutoff_nml / select_localization, vloc_types, vloc_factors, vloc_nlev, vloc_surface_level, &
                            vloc_coord, vloc_points, vloc_mismatch

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
real(r8) :: d
character(len=256) :: msg

if (trim(vloc_mismatch) /= 'error' .and. trim(vloc_mismatch) /= 'warn') call error_handler(E_ERR, &
   'comp_cov_factor', 'vloc_mismatch must be error or warn', source)
allocate(vloc_index(max(1, get_num_types_of_obs())))
vloc_index = 0
do t = 1, max_vloc_types
   if (len_trim(vloc_types(t)) == 0) cycle
   if (vloc_nlev < 1 .or. vloc_nlev > max_vloc_lev) then
      write(msg, '(A,I0)') 'vloc_types needs 1 <= vloc_nlev <= ', max_vloc_lev
      call error_handler(E_ERR, 'comp_cov_factor', msg, source)
   endif
   if (vloc_surface_level < 1 .or. vloc_surface_level > vloc_nlev) call error_handler(E_ERR, &
      'comp_cov_factor', 'vloc_surface_level must be a row of the tables (1..vloc_nlev)', source)
   select case (trim(vloc_coord(t)))
      case ('level')
         vloc_ic(t) = 1
      case ('pressure')
         vloc_ic(t) = 2
      case ('height')
         vloc_ic(t) = 3
      case ('scale_height')
         vloc_ic(t) = 4
      case default
         call error_handler(E_ERR, 'comp_cov_factor', 'vloc_coord must be level, pressure, height '// &
            'or scale_height; got '//trim(vloc_coord(t))//' for '//trim(vloc_types(t)), source)
   end select
   if (vloc_ic(t) > 1) then
      ! the rows' coordinates: all set, strictly monotonic, positive for pressure
      if (any(vloc_points(1:vloc_nlev, t) < -0.5_r8*huge(1.0_r8))) call error_handler(E_ERR, &
         'comp_cov_factor', 'vloc_points(1:vloc_nlev) must be set for '//trim(vloc_types(t)), source)
      if (vloc_ic(t) == 2 .and. any(vloc_points(1:vloc_nlev, t) <= 0.0_r8)) call error_handler(E_ERR, &
         'comp_cov_factor', 'vloc_points must be positive pressures (Pa) for '//trim(vloc_types(t)), source)
      if (vloc_nlev > 1) then
         d = vloc_points(2, t) - vloc_points(1, t)
         if (d == 0.0_r8 .or. any((vloc_points(2:vloc_nlev, t) - vloc_points(1:vloc_nlev-1, t))*d <= 0.0_r8)) &
            call error_handler(E_ERR, 'comp_cov_factor', &
               'vloc_points must be strictly monotonic for '//trim(vloc_types(t)), source)
      endif
   endif
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
if (vloc_ic(t) == 1 .and. is_vertical(target_loc, 'LEVEL')) then
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
else if (is_vertical(target_loc, 'UNDEFINED')) then
   vloc_factor = 1.0_r8
else if ((vloc_ic(t) == 2 .and. is_vertical(target_loc, 'PRESSURE')) .or. &
         (vloc_ic(t) == 3 .and. is_vertical(target_loc, 'HEIGHT')) .or. &
         (vloc_ic(t) == 4 .and. is_vertical(target_loc, 'SCALE_HEIGHT'))) then
   vloc_factor = vloc_interp(t, query_location(target_loc, 'VLOC'))
else
   call vloc_mismatch_found(t, target_loc)
endif

end function vloc_factor

!----------------------------------------------------------------------
!> Factor of a coordinate table (vloc_coord pressure, height or scale height)
!> at vertical value v: linear between rows (in log(p) for pressure), the end
!> rows' factors beyond the ends

function vloc_interp(t, v)

integer,  intent(in) :: t
real(r8), intent(in) :: v
real(r8)             :: vloc_interp

real(r8) :: x, x1, x2
integer  :: k, n
logical  :: up

n = vloc_nlev
if (n == 1) then
   vloc_interp = vloc_factors(1, t)
   return
endif
if (vloc_ic(t) == 2) then
   x = log(max(v, tiny(1.0_r8)))
else
   x = v
endif
up = vloc_points(n, t) > vloc_points(1, t)
do k = 1, n - 1
   x1 = vloc_points(k, t); x2 = vloc_points(k + 1, t)
   if (vloc_ic(t) == 2) then
      x1 = log(x1); x2 = log(x2)
   endif
   if (k == 1 .and. ((up .and. x <= x1) .or. (.not. up .and. x >= x1))) then
      vloc_interp = vloc_factors(1, t)
      return
   endif
   if ((up .and. x <= x2) .or. (.not. up .and. x >= x2)) then
      vloc_interp = vloc_factors(k, t) + (x - x1)/(x2 - x1)*(vloc_factors(k + 1, t) - vloc_factors(k, t))
      return
   endif
end do
vloc_interp = vloc_factors(n, t)

end function vloc_interp

!----------------------------------------------------------------------
!> A table met a target located in a vertical other than its own

subroutine vloc_mismatch_found(t, target_loc)

integer,             intent(in) :: t
type(location_type), intent(in) :: target_loc

character(len=16), parameter :: vname(-2:4) = [character(len=16) :: 'undefined', 'surface', &
                                               '?', 'level', 'pressure', 'height', 'scale_height']
character(len=256) :: msg
integer :: wv

wv = 0
if (is_vertical(target_loc, 'LEVEL'))        wv = 1
if (is_vertical(target_loc, 'PRESSURE'))     wv = 2
if (is_vertical(target_loc, 'HEIGHT'))       wv = 3
if (is_vertical(target_loc, 'SCALE_HEIGHT')) wv = 4
write(msg, '(7A)') 'vertical localization table of ', trim(vloc_types(t)), ' is in ', &
   trim(vloc_coord(t)), ' but a target is located by ', trim(vname(wv)), &
   ' (set vloc_coord to match, convert the state verticals, or vloc_mismatch = ''warn'')'
if (trim(vloc_mismatch) == 'error') call error_handler(E_ERR, 'comp_cov_factor', msg, source)
if (.not. vloc_warned(t, wv)) then
   vloc_warned(t, wv) = .true.
   call error_handler(E_MSG, 'comp_cov_factor', 'WARNING: '//trim(msg)//'; factor 1 used', source)
endif

end subroutine vloc_mismatch_found

end module cov_cutoff_mod

