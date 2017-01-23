module ocean_nphysicsC_mod
#define COMP isc:iec,jsc:jec
#define COMPXL isc-1:iec,jsc:jec
#define COMPYL isc:iec,jsc-1:jec
#define COMPXLYL isc-1:iec,jsc-1:jec
! 
!<CONTACT EMAIL="GFDL.Climate.Model.Info@noaa.gov"> Stephen M. Griffies
!</CONTACT>
!
!<REVIEWER EMAIL="tim.leslie@gmail.com"> Tim Leslie
!</REVIEWER>
!  
!<OVERVIEW>
! Thickness weighted and density weighted time tendency for tracer 
! from Laplacian neutral diffusion + Laplacian skew-diffusion.
!</OVERVIEW>
!
!<DESCRIPTION>
! This module computes the cell thickness weighted and density 
! weighted tracer tendency from small angle Laplacian neutral diffusion
! plus Laplacian skew-diffusion.  The algorithms for neutral diffusion
! are based on mom4p0d methods.  The algorithm for neutral skewsion 
! are based on a projection onto a few of the lowest baroclinic 
! modes. This module is experimental, and should be used with caution. 
!</DESCRIPTION>
!
! <INFO>
!
! <REFERENCE>
! S.M. Griffies, A. Gnanadesikan, R.C. Pacanowski, V. Larichev, 
! J.K. Dukowicz,  and R.D. Smith
! Isoneutral diffusion in a z-coordinate ocean model
! Journal of Physical Oceanography (1998) vol 28 pages 805-830
! </REFERENCE>
!
! <REFERENCE>
! S.M. Griffies 
! The Gent-McWilliams Skew-flux 
! Journal of Physical Oceanography (1998) vol 28 pages 831-841
! </REFERENCE>
!
! <REFERENCE>
! R. Ferrari, S.M. Griffies, A.J.G. Nurser, and G.K. Vallis  
! A boundary value problem for the parameterized mesoscale eddy transport
! Ocean Modelling, Volume 32, 2010, Pages 143-156.
! </REFERENCE>
!
! <REFERENCE>
! S.M. Griffies 
! Fundamentals of Ocean Climate Models (2004)
! Princeton University Press 
! </REFERENCE>
!
! <REFERENCE>
! S.M. Griffies: Elements of MOM (2012)
! </REFERENCE>
!
! <REFERENCE>
! D.B. Chelton,  R.A. deSzoeke, M.G. Schlax, K.E. Naggar, N. Siwertz
! Geographical Variability of the First Baroclinic Rossby Radius of Deformation
! Journal of Physical Oceanography (1998) vol 28 pages 433-460 
! </REFERENCE>
!
! <REFERENCE>
! G. Danabasoglu and J. C. McWilliams
! Sensitivity of the global ocean circulation to 
! parameterizations of mesoscale tracer transports
! Journal of Climate (1995) vol 8 pages 2967--2987 
! </REFERENCE>
!
! <REFERENCE>
! Gerdes, Koberle, and Willebrand
! The influence of numerical advection schemes on the results of ocean
! general circulation models, Climate Dynamics (1991), vol. 5, 
! pages 211--226. 
! </REFERENCE>
!
! <NOTE>
! Numerical implementation of the flux components follows the triad 
! approach documented in the references and implemented in MOM2 and MOM3.  
! The MOM algorithm accounts for partial bottom cells and generalized
! orthogonal horizontal coordinates and general vertical levels.  
! </NOTE> 
!
! <NOTE> 
! In steep neutral slope regions, neutral diffusive fluxes are tapered
! to zero with the tanh taper of Danabasoglu and McWilliams (1995) or the 
! quadratic scheme of Gerdes, Koberle, and Willebrand.  
!
! Traditional tapering is not required for the skew fluxes computed in 
! this module.  
! </NOTE> 
!
! </INFO>
!
!<NAMELIST NAME="ocean_nphysicsC_nml">
!
!  <DATA NAME="use_this_module" TYPE="logical">
!  Must be true to use this module. Default is false.
!  </DATA> 
!  <DATA NAME="debug_this_module" TYPE="logical">
!  For printing starting and ending checksums for restarts
!  </DATA> 
!
!  <DATA NAME="epsln_bv_freq" UNITS="kg/m4" TYPE="real">
!  Minimum buoyancy frequency accepted for the computation of 
!  baroclinic modes. Default epsln_bv_freq=1e-10.  Note there 
!  is also a minimum drhodz set in ocean_density.F90 via the
!  nml epsln_drhodz in that module.  We provide yet another minimum
!  here in case we need an extra regularization for the amplitude
!  of the baroclinic modes.  
!  </DATA>
!
!  <DATA NAME="do_neutral_diffusion" TYPE="logical">
!  To compute tendency from neutral diffusion.
!  Default do_neutral_diffusion=.true. 
!  </DATA> 
!  <DATA NAME="do_gm_skewsion" TYPE="logical">
!  To compute tendency from GM skewsion. Default do_gm_skewsion=.true. 
!  </DATA> 
!  <DATA NAME="gm_skewsion_modes" TYPE="logical">
!  To compute tendency from GM skewsion using streamfunction established
!  by baroclinic modes. Default gm_skewsion_modes=.false.
!  </DATA> 
!  <DATA NAME="gm_skewsion_bvproblem" TYPE="logical">
!  To compute tendency from GM skewsion using streamfunction established
!  by a boundary value problem. Default gm_skewsion_bvproblem=.true.
!  </DATA> 
!
!  <DATA NAME="number_bc_modes" TYPE="integer">
!  The number of baroclinic modes used to construct the eddy induced 
!  streamfunction. Default number_bc_modes=1.
!  </DATA> 
!  <DATA NAME="bvp_bc_mode" TYPE="integer">
!  The particular baroclinic mode used to construct the BVP streamfunction.
!  If bvp_bc_mode=0, then will set bc_speed=0 when computing the BVP streamfunction.
!  Default bvp_bc_mode=1.  
!  </DATA> 
!  <DATA NAME="bvp_constant_speed" TYPE="logical">
!  For taking a constant speed to be used for the calculation 
!  of the BVP streamfunction. Default bvp_constant_speed=.false.  
!  </DATA> 
!  <DATA NAME="bvp_speed" UNITS="m/s" TYPE="real">
!  For setting the speed weighting the second order derivative operator 
!  in the BVP streamfunction method: 
!  c^2 = max[bvp_min_speed, (bvp_speed-c_mode)^2].
!  If bvp_constant_speed, then  c^2 = bvp_speed^2. 
!  Default bvp_speed=0.0, in which case c^2 = c_mode^2.  
!  </DATA> 
!  <DATA NAME="bvp_min_speed" UNITS="m/s" TYPE="real">
!  For setting a minimum speed for use with the calculation 
!  of the BVP streamfunction. We need  bvp_min_speed>0 to ensure
!  that the second order derivative operator contributes to the 
!  calculation of the streamfunction.  
!  Default bvp_min_speed=0.1.  
!  </DATA> 
! 
!  <DATA NAME="bv_freq_smooth_vert" TYPE="logical">
!  To smooth the buoyancy frequency for use in 
!  computing the baroclinic modes. Generally this field has already 
!  been smooted in ocean_density_mod, but we maintain the possibility of 
!  further smoothing here.  Default bv_freq_smooth_vert=.false.
!  </DATA>
!  <DATA NAME="num_121_passes" TYPE="integer">
!  The number of 121 passes used to smooth buoyancy frequency when 
!  bv_freq_smooth_vert=.true.  Default num_121_passes=1.
!  </DATA>
!  <DATA NAME="min_bc_speed"  UNITS="m/s"  TYPE="real">
!  The minimum speed used for computing the baroclinic modes. 
!  Default min_bc_speed=1e-6
!  </DATA> 
!
!  <DATA NAME="smooth_bc_modes" TYPE="logical">
!  For doing a vertical 1-2-1 smoothing on the baroclinic modes
!  prior to normalization.  This is useful to reduce noise.
!  Default smooth_bc_modes=.false.
!  </DATA> 
!
!  <DATA NAME="smooth_psi" TYPE="logical">
!  For doing a horizontal 1-2-1 smoothing on the psix and psiy fields. 
!  This is useful to reduce noise. Default smooth_psi=.true.
!  </DATA> 
!
!  <DATA NAME="regularize_psi" TYPE="logical">
!  To reduce the magnitude of psi in regions of weak stratification, 
!  using the slope = smax_psi to set the overall scale of the max allowed
!  for psi. Default regularize_psi=.true.
!  </DATA> 
!  <DATA NAME="smax_modes" TYPE="real">
!  Maximum slope used for setting the overall scale of a modal 
!  contribution to the parameterized transport.   
!  Default smax_psi=0.1.  
!  </DATA> 
! 
!  <DATA NAME="diffusion_all_explicit" TYPE="logical">
!  To compute all contributions from neutral diffusion explicitly in time, including
!  the K33 diagonal piece.  This approach is available only when have small time 
!  steps and/or running with just a single tracer.  It is for testing purposes. 
!  </DATA> 
!
!  <DATA NAME="neutral_physics_limit" TYPE="logical">
!  When tracer falls outside a specified range, revert to horizontal 
!  diffusive fluxes at this cell. This is an ad hoc and incomplete attempt
!  to maintain monotonicity with the neutral physics scheme.  
!  Default neutral_physics_limit=.true.
!  </DATA> 
!  <DATA NAME="tmask_neutral_on" TYPE="logical">
!  If .true. then this logical reduces the neutral diffusive fluxes to 
!  horizontal/vertical diffusion next to boundaries.  
!  This approach has been found to reduce spurious 
!  extrema resulting from truncation of triads used to compute 
!  a neutral flux component.   
!  Default tmask_neutral_on=.false.
!  </DATA> 
!
!  <DATA NAME="dm_taper" TYPE="logical">
!  Set to true to use the tanh tapering scheme of Danabasoglu and McWilliams.
!  Default is true. 
!  </DATA> 
!  <DATA NAME="gkw_taper" TYPE="logical">
!  Set to true to use the quadradic tapering scheme of Gerdes, Koberle, and Willebrand.
!  Default is false. 
!  </DATA> 
!
!  <DATA NAME="neutral_eddy_depth" TYPE="logical">
!  Compute eddy_depth according to depth over which eddies feel the ocean surface.
!  Default neutral_eddy_depth=.true. 
!  </DATA> 
!
!  <DATA NAME="turb_blayer_min" TYPE="real">
!  Minimum depth of a surface turbulent boundary layer
!  used in the transition of the neutral diffusion fluxes
!  to the surface.  Note that in mom4p0, 
!  turb_blayer_min was always set to zero. 
!  </DATA> 
!
!  <DATA NAME="use_neutral_slopes_potrho" TYPE="logical">
!  To compute the neutral slopes based on globally referenced potential 
!  density rather than locally referenced potential density. This approach
!  is meant solely for sensitivity studies; it is not meant for realistic 
!  simulations.  
!  Default use_neutral_slopes_potrho=.false. 
!  </DATA> 
!  <DATA NAME="neutral_slopes_potrho_press" UNITS="dbar" TYPE="real">
!  The reference pressure used to compute neutral slopes when setting
!  use_neutral_slopes_potrho=.true. 
!  Default neutral_slopes_potrho_press=2000.0
!  </DATA> 
!
!  <DATA NAME="smooth_advect_transport" TYPE="logical">
!  For doing a horizontal 1-2-1 smoothing on the diagnosed  
!  uhrho_et_gm and vhrho_nt_gm fields.  
!  Default smooth_advect_transport=.true.
!  </DATA> 
!  <DATA NAME="smooth_advect_transport_num" TYPE="integer">
!  Number of iterations for the smooothing of horizontal transport. 
!  Default smooth_advect_transport_num=2.
!  </DATA> 
!
!</NAMELIST>

use constants_mod,           only: epsln, pi
use diag_manager_mod,        only: register_diag_field, register_static_field, need_data, send_data
use fms_mod,                 only: FATAL, WARNING, NOTE
use fms_mod,                 only: open_namelist_file, check_nml_error, close_file, write_version_number
use fms_io_mod,              only: string 
use mpp_domains_mod,         only: mpp_update_domains
use mpp_domains_mod,         only: CGRID_NE, WUPDATE, SUPDATE, EUPDATE, NUPDATE
use mpp_domains_mod,         only: cyclic_global_domain, global_data_domain 
use mpp_mod,                 only: input_nml_file, mpp_error, stdout, stdlog
use mpp_mod,                 only: mpp_clock_id, mpp_clock_begin, mpp_clock_end, CLOCK_ROUTINE
use time_manager_mod,        only: set_time, time_type, increment_time, operator ( + )

use ocean_density_mod,       only: density_derivs 
use ocean_domains_mod,       only: get_local_indices, set_ocean_domain
use ocean_nphysics_util_mod, only: ocean_nphysics_coeff_init, ocean_nphysics_coeff_end
use ocean_nphysics_util_mod, only: neutral_slopes, tracer_derivs 
use ocean_nphysics_util_mod, only: compute_eady_rate, compute_baroclinicity
use ocean_nphysics_util_mod, only: compute_rossby_radius, compute_bczone_radius
use ocean_nphysics_util_mod, only: compute_diffusivity, ocean_nphysics_util_restart
use ocean_nphysics_util_mod, only: transport_on_nrho_gm, transport_on_rho_gm, transport_on_theta_gm
use ocean_nphysics_util_mod, only: cabbeling_thermob_tendency
use ocean_nphysics_util_mod, only: compute_eta_tend_gm90
use ocean_nphysics_util_mod, only: watermass_diag_init, watermass_diag_ndiffuse, watermass_diag_sdiffuse
use ocean_nphysics_new_mod,  only: neutral_slopes=>gradrho
use ocean_nphysics_new_mod,  only: neutral_slopes=>calc_neutral_slope_vector
use ocean_operators_mod,     only: FAX, FAY, FMX, FMY, BDX_ET, BDY_NT
use ocean_parameters_mod,    only: missing_value, onehalf, onefourth, oneeigth, DEPTH_BASED
use ocean_parameters_mod,    only: rho0r, rho0, grav
use ocean_tracer_diag_mod,   only: diagnose_eta_tend_3dflux 
use ocean_types_mod,         only: ocean_grid_type, ocean_domain_type, ocean_density_type
use ocean_types_mod,         only: ocean_prog_tracer_type, ocean_thickness_type
use ocean_types_mod,         only: ocean_time_type, ocean_time_steps_type
use ocean_types_mod,         only: tracer_2d_type, tracer_3d_0_nk_type, tracer_3d_1_nk_type 
use ocean_util_mod,          only: write_note, write_line, write_warning
use ocean_util_mod,          only: diagnose_2d, diagnose_3d
use ocean_workspace_mod,     only: wrk1_2d, wrk1_v2d, wrk2_v2d
use ocean_workspace_mod,     only: wrk1, wrk2, wrk3, wrk4
use ocean_workspace_mod,     only: wrk1_v, wrk2_v, wrk3_v, wrk4_v

implicit none

public ocean_nphysicsC_init
public ocean_nphysicsC_end
public nphysicsC
public ocean_nphysicsC_restart

private 

type(ocean_grid_type), pointer   :: Grd => NULL()
type(ocean_domain_type), pointer :: Dom => NULL()
type(ocean_domain_type), save    :: Dom_flux

integer :: num_prog_tracers = 0

#include <ocean_memory.h>

#ifdef MOM_STATIC_ARRAYS

#else

#endif

! for diagnostics 
real, dimension(:,:,:), allocatable :: uhrho_et_hrm  ! i-component of advective mass transport
real, dimension(:,:,:), allocatable :: vhrho_nt_hrm  ! j-component of advective mass transport

! introduce following derived types so that do not need to know num_prog_tracers at compile time 
type(tracer_3d_1_nk_type), dimension(:), allocatable  :: dTdx    ! tracer partial derivative (tracer/m)
type(tracer_3d_1_nk_type), dimension(:), allocatable  :: dTdy    ! tracer partial derivative (tracer/m)
type(tracer_3d_0_nk_type), dimension(:), allocatable  :: dTdz    ! tracer partial derivative (tracer/m)

integer :: index_temp
integer :: index_salt
integer :: neutralrho_nk

character(len=*), parameter :: FILENAME=&
     __FILE__

logical :: module_is_initialized = .FALSE.

! time step settings 
real    :: dtime
real    :: two_dtime_inv

! constants
real :: pi_r
real :: grav_rho0_r
real :: grav_rho0r
real :: sqrt_grav

! vertical coordinate 
integer :: vert_coordinate_class

! lower and upper depth for vertically averaging ocean properties.
! read into ocean_nphysics_util_nml
real :: agm_closure_upper_depth
real :: agm_closure_lower_depth

! for eta_tend diagnostics (internally set)
logical :: diagnose_eta_tend_ndiff_flx=.false. 
logical :: diagnose_eta_tend_gm_flx   =.false. 

! for diagnosing advective GM transport 
logical :: diag_advect_transport       = .false. 
logical :: smooth_advect_transport     = .true. 
integer :: smooth_advect_transport_num = 2


!**************nml settings**************

! FIXME: explain this option
logical :: do_hrm                = .true.

! for the module as a whole 
logical :: use_this_module   = .false.
logical :: debug_this_module = .false.

!**************end of nml settings**************

namelist /ocean_nphysicsC_nml/ use_this_module, debug_this_module, &
          do_hrm             &

contains

!#######################################################################
! <SUBROUTINE NAME="ocean_nphysicsC_init">
!
! <DESCRIPTION>
! Initialize the neutral physics module by registering fields for 
! diagnostic output and performing some numerical checks to see 
! that namelist settings are appropriate. 
! </DESCRIPTION>
!
subroutine ocean_hrm_init(Grid, Domain, Time, Time_steps, Thickness, Dens, T_prog,&
           velocity, debug)

  type(ocean_grid_type),        intent(in), target   :: Grid
  type(ocean_domain_type),      intent(in), target   :: Domain
  type(ocean_time_type),        intent(in)           :: Time
  type(ocean_time_steps_type),  intent(in)           :: Time_steps
  type(ocean_thickness_type),   intent(in)           :: Thickness
  type(ocean_density_type),     intent(in)           :: Dens
  type(ocean_prog_tracer_type), intent(inout)        :: T_prog(:)
  type(ocean_velocity_type),    intent(inout)        :: Velocity

end subroutine ocean_hrm_init
! </SUBROUTINE>  NAME="ocean_hrm_init"

!#######################################################################
! <SUBROUTINE NAME="nphysicsC">
!
! <DESCRIPTION>
! This function computes the thickness weighted and density weighted
! time tendency for tracer from neutral physics.  Full discussion
! and details are provided by Griffies (2008). 
!
! Here is a brief summary.  
!
!---How the neutral diffusive flux components are computed:
!
! The vertical flux component is split into diagonal (3,3) and 
! off-diagonal (3,1) and (3,2) terms. The off-diagonal (3,1) and (3,2) 
! terms are included explicitly in time. The main contribution from the 
! (3,3) term to the time tendency is included implicitly in time 
! along with the usual contribution from diapycnal processes 
! (vertical mixing schemes).  This is the K33_implicit term.
! This approach is necessary with high vertical resolution, as 
! noted by Cox (1987).  However, splitting the vertical flux into 
! an implicit and explicit piece compromises the 
! integrity of the vertical flux component (see Griffies et al. 1998).
! So to minimize the disparity engendered by this split, the portion of 
! K33 that can be stably included explicitly in time is computed along 
! with the (3,1) and (3,2) terms. 
! 
! All other terms in the mixing tensor are included explicitly in time
! using a forward time step as required for temporal stability of 
! numerical diffusive processes.  
!
! The off-diagonal terms in the horizontal flux components, and all terms
! in the vertical flux component, are tapered in regions of steep neutral
! slope according to the requirements of linear stability.  MOM allows for 
! choice of two tapering schemes:
!
! (a) the tanh taper of Danabasoglu and McWilliams (1995)
! (b) the quadratic scheme of Gerdes, Koberle, and Willebrand (1991)
!
! Linear stability is far less stringent on the diagonal (1,1) and (2,2)
! part of the horizontal flux.  Indeed, these terms in practice need
! not be tapered in steep sloped regions. 
!
!---How the skew diffusive flux components are computed:
!
! The skew flux components are purely off-diagonal.  
! They are computed based on a vector streamfunction which 
! is built from a sum of baroclinic modes. 
! It is this part of the calculation that differs from 
! ocean_nphysicsA and ocean_nphysicsB.
! </DESCRIPTION>
!
subroutine horizontal_residual_mean (Time, Thickness, Dens, rho, T_prog, &
                      Velocity)

  type(ocean_time_type),        intent(in)    :: Time
  type(ocean_thickness_type),   intent(in)    :: Thickness
  type(ocean_density_type),     intent(in)    :: Dens
  real, dimension(isd:,jsd:,:), intent(in)    :: rho
  real, dimension(isd:,jsd:),   intent(in)    :: surf_blthick
  type(ocean_prog_tracer_type), intent(inout) :: T_prog(:)
  type(ocean_velocity_type),    intent(inout)        :: Velocity

  if (.not. use_this_module) return

  if ( .not. module_is_initialized ) then 
     call mpp_error(FATAL, &
     '==>Error from ocean_nphysicsC (neutral_physics): needs initialization')
  endif 

  if (size(T_prog(:)) /= num_prog_tracers) then 
    call mpp_error(FATAL, &
    '==>Error from ocean_nphysicsC (neutral_physics): inconsistent size of T_prog')
  endif 

  taum1 = Time%taum1
  tau   = Time%tau

  ! time dependent delqc geometric factor 
  do k=1,nk
    delqc(:,:,k,0) = Grd%fracdz(k,0)*Thickness%rho_dzt(:,:,k,tau)
    delqc(:,:,k,1) = Grd%fracdz(k,1)*Thickness%rho_dzt(:,:,k,tau)
  enddo

end subroutine horizontal_residual_mean
! </SUBROUTINE> NAME="nphysicsC"


!#######################################################################
! <SUBROUTINE NAME="compute_hrm_transport">
!
! <DESCRIPTION>
! Diagnose advective mass transport from HRM. 
!
! </DESCRIPTION>
!
subroutine compute_hrm_transport(Time, Dens, Thickness, Grid, T_prog, Velocity)

  type(ocean_time_type),      intent(in) :: Time
  type(ocean_density_type),   intent(in) :: Dens
  type(ocean_thickness_type), intent(in) :: Thickness
  type(ocean_velocity_type), intent(in) :: Velocity

  real, dimension(isd:,jsd:,:,:),  intent(in) :: dTdx
  real, dimension(isd:,jsd:,:,:),  intent(in) :: dTdy
  real, dimension(isd:,jsd:,0:,:), intent(in) :: dTdz

  real, dimension(isd:,jsd:,:,0:,0:), intent(out) :: slopex_xz
  real, dimension(isd:,jsd:,:,0:,0:), intent(out) :: slopey_yz
  real, dimension(isd:,jsd:,:,0:),    intent(out) :: absslope_z
  real, dimension(isd:,jsd:,:),       intent(out) :: absslope

  ! neutral density derivatives 
  real, dimension(isd:ied,jsd:jed,nk,0:1) :: drhodx_x
  real, dimension(isd:ied,jsd:jed,nk,0:1) :: drhody_y
  real, dimension(isd:ied,jsd:jed,nk,0:1) :: drhodz_z 

  integer :: stdoutunit
  stdoutunit=stdout()

  print*, 'CALLING compute_hrm_transport'

  call tracer_gradients(Time, Grid, Time%taum1, T_prog, Thickness, &
         dTdx, dTdy, dTdz)

  ! Calculate density gradients
  call gradrho(Dens, dTdx, dTdy, dTdz, drhodx_x, drhody_y, drhodz_z)

  ! Calculate neutral slopes
  call neutral_slopes(Time, Grid, drhodx_x, drhody_y, drhodz_z, &
         slopex_xz, slopey_yz, absslope_z, absslope)

  ! i-component of hrm mass transport. units ((kg/m^3)*m^2/sec)
  uhrho_et_hrm(:,:,:) = 0.0
  ! j-component of hrm mass transport. units ((kg/m^3)*m^2/sec)
  vhrho_nt_hrm(:,:,:) = 0.0

  wrk1_v(:,:,:,:) = 0.0 ! (rho*psix,rho*psiy) = 0.0
  do k=1,nk
     do j=jsc,jec
        do i=isc,iec
           uhrho_et_hrm(i,j,k) = 0.0
           vhrho_nt_hrm(i,j,k) = 0.0
        enddo
     enddo
  enddo

  ! update domains as per C-grid fluxes
  call mpp_update_domains(uhrho_et_gm(:,:,:), vhrho_nt_gm(:,:,:),&
                          Dom_flux%domain2d, gridtype=CGRID_NE)

  call diagnose_3d(Time, Grd, id_uhrho_et_hrm, uhrho_et_hrm(:,:,:))
  call diagnose_3d(Time, Grd, id_vhrho_nt_hrm, vhrho_nt_hrm(:,:,:))

end subroutine compute_hrm_transport
! </SUBROUTINE> NAME="compute_advect_transport"


!#######################################################################
! <SUBROUTINE NAME="ocean_nphysicsC_end">
!
! <DESCRIPTION>
! Write to restart.
! </DESCRIPTION>
!
subroutine ocean_hrm_end(Time)

  type(ocean_time_type), intent(in) :: Time
  
  if(.not. use_this_module) return

  call ocean_nphysics_coeff_end(Time, agm_array, aredi_array, rossby_radius, rossby_radius_raw, bczone_radius)

end subroutine ocean_nphysicsC_end
! </SUBROUTINE> NAME="ocean_nphysicsC_end"


end module ocean_nphysicsC_mod
