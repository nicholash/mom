module ocean_hrm_mod
#define COMP isc:iec,jsc:jec
#define COMPXL isc-1:iec,jsc:jec
#define COMPYL isc:iec,jsc-1:jec
#define COMPXLYL isc-1:iec,jsc-1:jec
#define COMPXLL isc-1:iec-1,jsc:jec
#define COMPYLL isc:iec,jsc-1:jec-1
! 
!<CONTACT EMAIL="yuehua.li@unsw.edu.au"> Yuehua Li
!</CONTACT>
!
!<REVIEWER EMAIL>
!</REVIEWER>
!  
!<OVERVIEW>
!
! This module is for diagnostics of Horizontal Residual Mean (HRM) method. 
! As an analogue of TRM, HRM aims to incoporate the unresolved spatial correlations into MOM.
! The HRM method is based on Taylor Series.
!
!</OVERVIEW>
!
!<DESCRIPTION>
!
! The HRM method doesn't need parameterization, but uses the density and velocity fields.
! It can use either B-grid or C-grid.
! This module is based on B-grid, but it is easy to adapt to C-grid.
!
!</DESCRIPTION>
!
! <INFO>
!
! <REFERENCE>
! Yuehua Li, Trevor McDougall, Shane Keating, Casimir de Lavergne and Gurvan Madec
! Horizontal Residual Mean(2017)
! </REFERENCE>
!
! <NOTE>
! </NOTE>
!
! </INFO>
!
!<NAMELIST NAME="ocean_hrm_nml">
!
!  <DATA NAME="use_this_module" TYPE="logical">
!  Must be true to use this module. Default is false.
!  </DATA> 
!  <DATA NAME="debug_this_module" TYPE="logical">
!  For printing starting and ending checksums for restarts
!  </DATA> 
!  <DATA NAME="vhrho_nt_hrm" TYPE="float">
!  The meridional HRM transport.
!  </DATA> 
!  <DATA NAME="uhrho_et_hrm" TYPE="float">
!  The zonal HRM transport.
!  </DATA> 
!
!</NAMELIST>

use constants_mod,   only: epsln, pi
use diag_manager_mod,        only: register_diag_field, register_static_field, need_data, send_data
use fms_mod,                 only: FATAL, WARNING, NOTE
use fms_mod,                 only: open_namelist_file, check_nml_error, close_file, write_version_number
use mpp_domains_mod,         only: mpp_update_domains
use mpp_domains_mod,         only: CGRID_NE, WUPDATE, SUPDATE, EUPDATE, NUPDATE
use mpp_mod,                 only:  mpp_pe
use ocean_operators_mod,     only: FDX_T, FDY_T, FMY, FMX
use ocean_parameters_mod,    only: missing_value
use mpp_domains_mod,         only: cyclic_global_domain, global_data_domain
use mpp_mod,                 only: input_nml_file, mpp_error, stdout, stdlog
use mpp_mod,                 only: mpp_clock_id, mpp_clock_begin, mpp_clock_end, CLOCK_ROUTINE
use ocean_domains_mod,       only: get_local_indices

use ocean_tracer_diag_mod,   only: diagnose_eta_tend_3dflux
use ocean_types_mod,         only: ocean_grid_type, ocean_domain_type, ocean_density_type
use ocean_types_mod,         only: ocean_prog_tracer_type, ocean_thickness_type
use ocean_types_mod,         only: ocean_velocity_type
!use ocean_types_mod,         only: ocean_external_mode_type
use ocean_types_mod,         only: ocean_time_type, ocean_time_steps_type
use ocean_types_mod,         only: tracer_3d_0_nk_type, tracer_3d_1_nk_type
use ocean_util_mod,          only: diagnose_2d, diagnose_3d
use ocean_parameters_mod, only: CONSERVATIVE_TEMP, POTENTIAL_TEMP, PRACTICAL_SALT, PREFORMED_SALT
use ocean_parameters_mod, only: rho0

implicit none

public ocean_hrm_init
public compute_hrm_transport
public ocean_hrm_end

private 

integer :: id_uhrho_et_hrm      =-1
integer :: id_vhrho_nt_hrm      =-1
integer :: id_T_trans_hrm       =-1
integer :: id_hrm_slopexN       =-1
integer :: id_hrm_slopexS       =-1
integer :: id_hrm_slopeyE       =-1
integer :: id_hrm_slopeyW       =-1
integer :: id_delt_rho_z        =-1
integer :: id_rhodz             =-1
integer :: id_rhodz_avg         =-1
integer :: id_rhodx_E_avg       =-1
integer :: id_rhodx_W_avg       =-1
integer :: id_rhodz_avg_E       =-1
integer :: id_rhodz_avg_W       =-1
integer :: id_delt_x            =-1
integer :: id_delt_y            =-1
integer :: id_delt_z            =-1
integer :: id_uy                =-1
integer :: id_vx                =-1
#include <ocean_memory.h>

! water mass transport and heat transport (for diagnostics) 
real, dimension(:,:,:), allocatable :: uhrho_et_hrm  ! i-component of advective mass transport
real, dimension(:,:,:), allocatable :: vhrho_nt_hrm  ! j-component of advective mass transport
real, dimension(:,:), allocatable :: T_trans       ! vertical sum of heat transport


integer :: index_temp
integer :: index_salt

character(len=*), parameter :: FILENAME=&
     __FILE__

logical :: module_is_initialized = .FALSE.


!**************nml settings**************

! for the module as a whole 
logical :: use_this_module   = .false.
logical :: debug_this_module = .false.

!**************end of nml settings**************

namelist /ocean_hrm_nml/ use_this_module, debug_this_module

! nml settings 

! for linear EOS 
logical :: eos_linear=.false.    
real    :: alpha_linear_eos=0.255
real    :: beta_linear_eos =0.0

! for teos10 or preteos10 eos
logical :: eos_preteos10 = .false.
logical :: eos_teos10    = .false.


contains

!#######################################################################
! <SUBROUTINE NAME="ocean_hrm_init">
!
! <DESCRIPTION>
! Initialize the HRM module by registering fields for 
! diagnostic output and performing some numerical checks to see 
! that namelist settings are appropriate. 
! </DESCRIPTION>
!
subroutine ocean_hrm_init(Grid, Domain, Time, Dens, T_prog,&
           velocity)

  type(ocean_grid_type),          intent(in), target   :: Grid
  type(ocean_domain_type),        intent(in), target   :: Domain
  type(ocean_time_type),          intent(in)           :: Time
  type(ocean_density_type),       intent(in)           :: Dens
  type(ocean_prog_tracer_type),   intent(inout)        :: T_prog(:)
  type(ocean_velocity_type),      intent(inout)        :: Velocity

  integer :: num_prog_tracers, n

  module_is_initialized = .TRUE.

  num_prog_tracers = size(T_prog)

  call get_local_indices(Domain,isd,ied,jsd,jed,isc,iec,jsc,jec)
  nk = Grid%nk

  index_temp=-1;index_salt=-1
  do n=1,num_prog_tracers
     if (T_prog(n)%name == 'temp') index_temp = n
     if (T_prog(n)%name == 'salt') index_salt = n
  enddo
  if (index_temp == -1 .or. index_salt == -1) then
     call mpp_error(FATAL, &
          '==>Error: temp and/or salt not identified in call to ocean_nphysics_util_init')
  endif


  allocate(vhrho_nt_hrm(isd:ied,jsd:jed,nk))
  allocate(uhrho_et_hrm(isd:ied,jsd:jed,nk))
  allocate(T_trans(isd:ied,jsd:jed))


  id_uhrho_et_hrm = register_diag_field ('ocean_model', 'uhrho_et_hrm', &
        Grid%tracer_axes_flux_x(1:3), Time%model_time,                  &
       'i-component of hrm mass transport',                            &
       '(kg/m^3)*m^2/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))

  id_vhrho_nt_hrm = register_diag_field ('ocean_model', 'vhrho_nt_hrm', &
        Grid%tracer_axes_flux_y(1:3), Time%model_time,                  &
       'j-component of hrm mass transport',                            &
       '(kg/m^3)*m^2/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))

  id_T_trans_hrm = register_diag_field ('ocean_model', 'T_trans_hrm', &
        Grid%tracer_axes(1:2), Time%model_time,                        &
       'j-component of hrm mass transport',                           &
       '(kg/m^3)*m^2/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))

  id_hrm_slopexS = register_diag_field ('ocean_model', 'hrm_slopexS', &
        Grid%tracer_axes_flux_x(1:3), Time%model_time,                        &
       'test diagnostic',                           &
       'unknown', missing_value=missing_value, range=(/-1.e10,1.e10/))

id_hrm_slopexN = register_diag_field ('ocean_model', 'hrm_slopexN', &
        Grid%tracer_axes_flux_x(1:3), Time%model_time,                        &
       'test diagnostic',                           &
       'unknown', missing_value=missing_value, range=(/-1.e10,1.e10/))

id_hrm_slopeyE = register_diag_field ('ocean_model', 'hrm_slopeyE', &
Grid%tracer_axes_flux_x(1:3), Time%model_time,                        &
'test diagnostic',                           &
'unknown', missing_value=missing_value, range=(/-1.e10,1.e10/))

id_hrm_slopeyW = register_diag_field ('ocean_model', 'hrm_slopeyW', &
Grid%tracer_axes_flux_x(1:3), Time%model_time,                        &
'test diagnostic',                           &
'unknown', missing_value=missing_value, range=(/-1.e10,1.e10/))

id_delt_rho_z = register_diag_field ('ocean_model', 'delt_rho_z', &
Grid%tracer_axes_flux_x(1:3), Time%model_time,                        &
'test diagnostic',                           &
'unknown', missing_value=missing_value, range=(/-1.e10,1.e10/))

id_rhodz = register_diag_field ('ocean_model', 'rho_dz', &
Grid%tracer_axes_flux_x(1:3), Time%model_time,                        &
'test diagnostic',                           &
'unknown', missing_value=missing_value, range=(/-1.e10,1.e10/))

id_rhodz_avg = register_diag_field ('ocean_model', 'rhodz_avg', &
Grid%tracer_axes_flux_x(1:3), Time%model_time,                        &
'test diagnostic',                           &
'unknown', missing_value=missing_value, range=(/-1.e10,1.e10/))

id_rhodx_E_avg = register_diag_field ('ocean_model', 'rhodx_E_avg', &
Grid%tracer_axes_flux_x(1:3), Time%model_time,                        &
'test diagnostic',                           &
'unknown', missing_value=missing_value, range=(/-1.e10,1.e10/))

id_rhodx_W_avg = register_diag_field ('ocean_model', 'rhodx_W_avg', &
Grid%tracer_axes_flux_x(1:3), Time%model_time,                        &
'test diagnostic',                           &
'unknown', missing_value=missing_value, range=(/-1.e10,1.e10/))

id_rhodz_avg_E = register_diag_field ('ocean_model', 'rhodz_avg_E', &
Grid%tracer_axes_flux_x(1:3), Time%model_time,                        &
'test diagnostic',                           &
'unknown', missing_value=missing_value, range=(/-1.e10,1.e10/))

id_rhodz_avg_W = register_diag_field ('ocean_model', 'rhodz_avg_W', &
Grid%tracer_axes_flux_x(1:3), Time%model_time,                        &
'test diagnostic',                           &
'unknown', missing_value=missing_value, range=(/-1.e10,1.e10/))

id_delt_x = register_diag_field ('ocean_model', 'delt_x', &
Grid%tracer_axes_flux_x(1:2), Time%model_time,                        &
'test diagnostic',                           &
'unknown', missing_value=missing_value, range=(/-1.e10,1.e10/))

id_delt_y = register_diag_field ('ocean_model', 'delt_y', &
Grid%tracer_axes_flux_x(1:2), Time%model_time,                        &
'test diagnostic',                           &
'unknown', missing_value=missing_value, range=(/-1.e10,1.e10/))

id_delt_z = register_diag_field ('ocean_model', 'delt_z', &
Grid%tracer_axes_flux_x(1:3), Time%model_time,                        &
'test diagnostic',                           &
'unknown', missing_value=missing_value, range=(/-1.e10,1.e10/))

id_uy = register_diag_field ('ocean_model', 'uy', &
Grid%tracer_axes_flux_x(1:3), Time%model_time,                        &
'test diagnostic',                           &
'unknown', missing_value=missing_value, range=(/-1.e10,1.e10/))

id_vx = register_diag_field ('ocean_model', 'vx', &
Grid%tracer_axes_flux_x(1:3), Time%model_time,                        &
'test diagnostic',                           &
'unknown', missing_value=missing_value, range=(/-1.e10,1.e10/))

end subroutine ocean_hrm_init
! </SUBROUTINE>  NAME="ocean_hrm_init"


!#######################################################################
! <SUBROUTINE NAME="compute_hrm_transport">
!
! <DESCRIPTION>
! Diagnose advective mass transport from HRM. 
!
! </DESCRIPTION>
!

subroutine compute_hrm_transport(Time, Dens, Grid, T_prog, Velocity)
!subroutine compute_hrm_transport(Time, Dens, Grid, T_prog, Velocity, Ext_mode)

  type(ocean_time_type),          intent(in) :: Time
  type(ocean_density_type),       intent(in) :: Dens
  type(ocean_grid_type),          intent(in) :: Grid
  type(ocean_prog_tracer_type),   intent(in) :: T_prog(:) 
  type(ocean_velocity_type),      intent(in) :: Velocity
  !type(ocean_external_mode_type), intent(in) :: Ext_mode


 real,dimension(isc:iec,jsc:jec,nk) :: salinity,temperature,rho,v,u
 real,dimension(isc:iec,jsc:jec-1,nk-1) :: UN,US,Uu,Ul
 real,dimension(isc:iec-1,jsc:jec,nk-1) :: VE,VW,Vu,Vl
 real,dimension(nk) :: h,delt_z
 real,dimension(isc:iec-1,jsc:jec) :: delt_x
 real,dimension(isc:iec,jsc:jec-1) :: delt_y
 real,dimension(isc:iec-1,jsc:jec,nk) :: rhodx, delt_rho_x,rhodx_W,rhodx_E,sa_e,temp_e
 real,dimension(isc:iec-1,jsc:jec-1,nk) :: rhodx_W_avg, rhodx_E_avg,rhody_N_avg,rhody_S_avg
 real,dimension(isc:iec,jsc:jec-1,nk) :: delt_rho_y,rhody, rhody_N,rhody_S,sa_n,temp_n
 real,dimension(isc:iec, jsc:jec, nk-2) ::  delt_rho_z, rhodz
 real,dimension(isc:iec, jsc:jec, nk) :: pressure
 real,dimension(isc:iec-1, jsc:jec-1, nk-2) ::slope_xz_faceE,slope_xz_faceW,slope_yz_faceN,slope_yz_faceS

 real,dimension(isd:ied,jsd:jed,nk) :: test_array, test_array1
 real,dimension(isd:ied,jsd:jed) :: test_array2
 real,dimension(nk,jsd:jed,isd:ied) :: test_array3
 real,dimension(isc:iec,jsc:jec-2,nk-1) :: uy
 real,dimension(isc:iec-2,jsc:jec,nk-1) :: vx

  integer :: tau
  real    :: press_standard,g
  integer :: i, j, k, m

  integer :: stdoutunit
  stdoutunit=stdout()

  if ( .not. module_is_initialized ) then 
     call mpp_error(FATAL, &
     '==>Error from ocean_hrm: needs initialization')
  endif 
  

  tau = Time%tau
  press_standard=10.1325
  g = 9.8

  ! i-component of hrm mass transport. units ((kg/m^3)*m^2/sec)
  uhrho_et_hrm(:,:,:) = 0.0
  ! j-component of hrm mass transport. units:kg/sec     units ((kg/m^3)*m^2/sec)
  vhrho_nt_hrm(:,:,:) = 0.0

  v = velocity%u(isc:iec,jsc:jec,:,2,tau) !The j component of velocity at t time level


  VE = (v(isc+1:,:,1:nk-1) + v(isc+1:,:,2:))/2.0 !Average the Eastern velocity onto depth z
  VW = (v(isc:iec-1,:,1:nk-1) + v(isc:iec-1,:,2:))/2.0 !Average the Western velocity onto depth z
  Vu = (v(isc:iec-1,:,1:nk-1) + v(isc+1:,:,1:nk-1))/2.0 !!Velocity at upper face Averaged to the center
  Vl = (v(isc:iec-1,:,2:) + v(isc+1:,:,2:))/2.0 !Velocity at lower face Averaged to the center


  u = velocity%u(isc:iec,jsc:jec,:,1,tau) !The i component of velocity at t time level


  UN = (u(:,jsc+1:,1:nk-1) + u(:,jsc+1:,2:))/2.0 !Average the Northern velocity onto depth z
  US = (u(:,jsc:jec-1,1:nk-1) + u(:,jsc:jec-1,2:))/2.0 !Average the Southern velocity onto depth z
  Uu = (u(:,jsc:jec-1,1:nk-1) + u(:,jsc+1:,1:nk-1))/2.0 !Average to the center
  Ul = (u(:,jsc:jec-1,2:) + u(:,jsc+1:,2:))/2.0 !Average to the center


  
  h = Grid%zt  ! distance from surface to grid point in level k (m)
  delt_z = Grid%dzt  ! thickness (m) of T cell at time tau/taup1 in ocean_types
  delt_x = Grid%dxte ! long-width between grid points at i+1 and i in T-cells (m)
  delt_y = Grid%dytn ! lat-width between grid points at j+1 and j in T-cells (m)

  ! Calculate neutral slopes
  rho = Dens%rho(isc:iec,jsc:jec,:,tau)
  
  !Calculate zonal slopes
  delt_rho_x =  rho(isc+1:iec,:,:) - rho (isc:iec-1,:,:)
  do k = 1,nk
    do j = jsc,jec
        do i = isc,iec-1
            if (delt_rho_x(i,j,k) > 10.0**(5.0)) then
                delt_rho_x(i,j,k) = 0.0
            endif
		enddo
	enddo
  enddo

  do k = 1,nk
		rhodx(:,:,k) = delt_rho_x(:,:,k)/ delt_x
  enddo
  
  rhodx_W = rhodx
  rhodx_E(isc:iec-2,:,:) = rhodx(isc+1:,:,:)
  rhodx_E(iec-1,:,:) = rhodx(iec-1,:,:)
  
  rhodx_W_avg = ( rhodx_W(:,jsc+1:,:) + rhodx_W(:,jsc:jec-1,:) )/2.0
  rhodx_E_avg = ( rhodx_E(:,jsc+1:,:) + rhodx_E(:,jsc:jec-1,:) )/2.0
  
  !!! From which variable should I read salinity and temperature fields?
  !!! What type of salinity and temperature are they?
  salinity    = Dens%rho_salinity(isc:iec,jsc:jec,:,tau)
  temperature = T_prog(index_temp)%field(isc:iec,jsc:jec,:,tau)
  
  ! The HRM calculation is on the northern and eastern faces of the grid box.
  sa_n = ( salinity(isc:iec,jsc:jec-1,:) + salinity(isc:iec,jsc+1:jec,:) )/2.0
  temp_n = ( temperature(isc:iec,jsc:jec-1,:) + temperature(isc:iec,jsc+1:jec,:) )/2.0
  
  sa_e = ( salinity(isc:iec-1,jsc:jec,:) + salinity(isc+1:iec,jsc:jec,:) )/2.0
  temp_e = ( temperature(isc:iec-1,jsc:jec,:) + temperature(isc+1:iec,jsc:jec,:) )/2.0
  
  !!! Getting the right pressure ?
  do k = 1,nk
  !pressure(:,:,k) = Dens%pressure_at_depth(isc:iec,jsc:jec,k+1) - Ext_mode%patm_t(isc:iec,jsc:jec,tau)
  !pressure(:,:,k) = Dens%pressure_at_depth(isc:iec,jsc:jec,k) - Ext_mode%patm_t(isc:iec,jsc:jec,tau)
  !pressure(:,:,k) = Dens%pressure_at_depth(isc:iec,jsc:jec,k) - press_standard
  pressure(:,:,k) = Dens%pressure_at_depth(isc:iec,jsc:jec,k)
  !pressure(:,:,k) = rho(isc:iec,jsc:jec,k)*g*h(k)*10 !unit is dbar
  enddo
  
  ! The density difference in vertical direction is the difference between two potential densities with the same reference pressure (at the target depth)
  do k = 2,nk-1
        do j = jsc,jec-1
                do i = isc,iec
                        delt_rho_z(i,j,k-1) = potential_density(sa_n(i,j, k+1 ),temp_n(i,j, k+1),pressure(i,j,k),T_prog) - potential_density(sa_n(i,j,k-1),temp_n(i,j,k-1),pressure(i,j,k),T_prog)
  				enddo
  		enddo
  enddo
  
  ! Calculating drho/dz
  do k =1,nk-2
		rhodz(:,:,k) = delt_rho_z(:,:,k) / (delt_z(k)/2 + delt_z(k+1) + delt_z(k+2)/2) !delt_z has one layer above and one under delt_rho_z
  enddo 
  
! Two slopes (Eastern and Western) are emanated from the center of the northern face
do k = 1,nk-2
        do j = jsc,jec-1
                do i = isc,iec-1
                        if (rhodz(i, j, k) == 0.0) then
                                slope_xz_faceE(i, j, k) = 0.0
                                slope_xz_faceW(i, j, k) = 0.0
                        else
                                slope_xz_faceE(i, j, k) = rhodx_E_avg(i,j,k+1)/rhodz(i, j, k)
                                slope_xz_faceW(i, j, k) = rhodx_W_avg(i,j,k+1)/rhodz(i, j, k)
                        endif
                enddo
        enddo
enddo


!Setting the limit of slope
do k = 1,nk-2
    do j = jsc,jec-1
        do i = isc,iec-1
            if ( slope_xz_faceW(i,j,k) > 10.0**(-2.0)) then
                slope_xz_faceW(i,j,k) = 10.0**(-2.0)
            else if (slope_xz_faceW(i,j,k) < - 10.0**(-2.0)) then
                slope_xz_faceW(i,j,k) = - 10.0**(-2.0)
            end if

            if ( slope_xz_faceE(i,j,k) > 10.0**(-2.0)) then
                slope_xz_faceE(i,j,k) = 10.0**(-2.0)
            else if (slope_xz_faceE(i,j,k) < - 10**(-2.0)) then
                slope_xz_faceE(i,j,k) = - 10.0**(-2.0)
            end if
        enddo
    enddo
enddo

  !Calculate meridional slopes
  delt_rho_y = rho(:,jsc+1:jec,:) - rho(:,jsc:jec-1,:)

do k = 1,nk
do j = jsc,jec-1
do i = isc,iec
if (delt_rho_y(i,j,k) > 10.0**(5.0)) then
delt_rho_y(i,j,k) = 0.0
endif
enddo
enddo
enddo

    do k = 1,nk
rhody(:,:,k) = delt_rho_y(:,:,k)/ delt_y
  enddo
  rhody_N (:,jsc:jec-2,:) = rhody(:,jsc+1:,:)
  rhody_N (:,jec-1,:) = rhody(:,jec-1,:)
  rhody_S  = rhody
  
  rhody_N_avg = ( rhody_N (isc+1:,:,:) + rhody_N(isc:iec-1,:,:) )/2.0
  rhody_S_avg = ( rhody_S (isc+1:,:,:) + rhody_S(isc:iec-1,:,:) )/2.0
  
  
    ! The density difference in vertical direction is the difference between two potential densities with the same reference pressure (at the target depth)
  do k = 2,nk-1
        do j = jsc,jec
                do i = isc,iec-1
                        delt_rho_z(i,j,k-1) = potential_density(sa_e(i,j, k+1 ),temp_e(i,j, k+1),pressure(i,j,k),T_prog) - potential_density(sa_e(i,j,k-1),temp_e(i,j,k-1),pressure(i,j,k),T_prog)
  				enddo
  		enddo
  enddo
  
  ! Calculating drho/dz
  do k =1,nk-2
		rhodz(:,:,k) = delt_rho_z(:,:,k) / (delt_z(k)/2 + delt_z(k+1) + delt_z(k+2)/2) !delt_z has one layer above and one under delt_rho_z
  enddo



  do k = 1,nk-2
        do j = jsc,jec-1
                do i = isc,iec-1
                        if (rhodz(i, j, k) == 0) then
                                slope_yz_faceN(i, j, k) = 0.0
                                slope_yz_faceS(i, j, k) = 0.0
                        else
                                slope_yz_faceN(i, j, k) =  rhody_N_avg(i,j,k+1)/rhodz(i, j, k)
                                slope_yz_faceS(i, j, k) =  rhody_S_avg(i,j,k+1)/rhodz(i, j, k)
                        endif
                enddo
        enddo
enddo

!Setting the limit of slope
do k = 1,nk-2
    do j = jsc,jec-1
        do i = isc,iec-1
            if ( slope_yz_faceN(i,j,k) > 10.0**(-2.0)) then
                slope_yz_faceN(i,j,k) = 10.0**(-2.0)
            else if (slope_yz_faceN(i,j,k) < - 10.0**(-2.0)) then
                slope_yz_faceN(i,j,k) = - 10.0**(-2.0)
            end if

            if ( slope_yz_faceS(i,j,k) > 10.0**(-2.0)) then
                slope_yz_faceS(i,j,k) = 10.0**(-2.0)
            else if (slope_yz_faceS(i,j,k) < - 10.0**(-2.0)) then
                slope_yz_faceS(i,j,k) = - 10.0**(-2.0)
            end if
        enddo
    enddo
enddo

                      
  !Calculate HRM from second depths
  do k=1,nk-2
     do j=jsc,jec-1
        do i=isc,iec-1
           uhrho_et_hrm(i,j,k) = 1.0/24.0*(UN(i+1, j, k)-US(i+1, j, k))*(slope_yz_faceN(i,j,k)+slope_yz_faceS(i,j,k))*(delt_y(i+1,j))**(2.0) + &
           1.0/48.0*(Uu(i+1, j, k)-Ul(i+1, j, k))/delt_z(k+1)*((slope_yz_faceN(i,j,k))**(2.0)+(slope_yz_faceS(i,j,k))**(2.0))*(delt_y(i+1, j))**(3.0) + &
           1.0/128.0*(Uu(i+1,j,k)-Ul(i+1,j,k))/delt_z(k+1)*(slope_yz_faceN(i,j,k)-slope_yz_faceS(i,j,k))**(2.0)*(delt_y(i+1, j))**(3.0)


           vhrho_nt_hrm(i,j,k) = 1.0/24.0*(VE(i, j+1, k)-VW(i, j+1, k))*(slope_xz_faceE(i,j,k)+slope_xz_faceW(i,j,k))*(delt_x(i,j+1))**(2.0) + &
           1.0/48.0*(Vu(i, j+1, k)-Vl(i, j+1, k))/delt_z(k+1)*((slope_xz_faceW(i,j,k))**(2.0)+(slope_xz_faceE(i,j,k))**(2.0))*(delt_x(i,j+1))**(3.0) + &
           1.0/128.0*(Vu(i,j+1,k)-Vl(i,j+1,k))/delt_z(k+1)*(slope_xz_faceE(i,j,k)-slope_xz_faceW(i,j,k))**(2.0)*(delt_x(i, j+1))**(3.0)
        enddo
     enddo
  enddo

  call compute_hrm_Ttransport(T_prog, T_trans, tau)

  if (id_uhrho_et_hrm > 0) then
    call diagnose_3d(Time, Grid, id_uhrho_et_hrm, uhrho_et_hrm(:,:,:))
  endif
  if (id_vhrho_nt_hrm > 0) then
    call diagnose_3d(Time, Grid, id_vhrho_nt_hrm, vhrho_nt_hrm(:,:,:))
  endif
  if (id_T_trans_hrm > 0) then
    call diagnose_2d(Time, Grid, id_T_trans_hrm, T_trans(:,:))
  endif

if (id_delt_rho_z > 0) then
test_array(:,:,:) = 0.0
test_array(isc:iec,jsc:jec,:nk-2) = delt_rho_z(:,:,:)
call diagnose_3d(Time, Grid, id_delt_rho_z, test_array(:, :, :))
endif

if (id_rhodz > 0) then
test_array(:,:,:) = 0.0
test_array(isc:iec,jsc:jec,:nk-2) = rhodz(:,:,:)
call diagnose_3d(Time, Grid, id_rhodz, test_array(:, :, :))
endif

!if (id_rhodz_avg > 0) then
!test_array(:,:,:) = 0.0
!test_array(isc:iec-1,jsc:jec-1,:nk-2) = rhodz_avg(:,:,:)
!call diagnose_3d(Time, Grid, id_rhodz_avg, test_array(:, :, :))
!endif

if (id_rhodx_E_avg > 0) then
test_array(:,:,:) = 0.0
test_array(isc:iec-1,jsc:jec-1,:nk) = rhodx_E_avg(:,:,:)
call diagnose_3d(Time, Grid, id_rhodx_E_avg, test_array(:, :, :))
endif

if (id_rhodx_W_avg > 0) then
test_array(:,:,:) = 0.0
test_array(isc:iec-1,jsc:jec-1,:nk) = rhodx_W_avg(:,:,:)
call diagnose_3d(Time, Grid, id_rhodx_W_avg, test_array(:, :, :))
endif

!if (id_rhodz_avg_E > 0) then
!test_array(:,:,:) = 0.0
!test_array(isc:iec-1,jsc:jec-1,:nk-2) = rhodz_avg_E(:,:,:)
!call diagnose_3d(Time, Grid, id_rhodz_avg_E, test_array(:, :, :))
!endif

!if (id_rhodz_avg_W > 0) then
!test_array(:,:,:) = 0.0
!test_array(isc:iec-1,jsc:jec-1,:nk-2) = rhodz_avg_W(:,:,:)
!call diagnose_3d(Time, Grid, id_rhodz_avg_W, test_array(:, :, :))
!endif

if (id_hrm_slopexS > 0) then
test_array1(:,:,:) = 0.0
test_array1(isc:iec-1,jsc:jec-1,:nk-2) = slope_yz_faceS(:,:,:)
call diagnose_3d(Time, Grid, id_hrm_slopexS, test_array1(:, :, :))
endif

if (id_hrm_slopexN > 0) then
test_array1(:,:,:) = 0.0
test_array1(isc:iec-1,jsc:jec-1,:nk-2) = slope_yz_faceN(:,:,:)
call diagnose_3d(Time, Grid, id_hrm_slopexN, test_array1(:, :, :))
endif

if (id_hrm_slopeyE > 0) then
test_array(:,:,:) = 0.0
test_array(isc:iec-1,jsc:jec-1,:nk-2) = slope_xz_faceE(:,:,:)
call diagnose_3d(Time, Grid, id_hrm_slopeyE, test_array(:, :, :))
endif

if (id_hrm_slopeyW > 0) then
test_array1(:,:,:) = 0.0
test_array1(isc:iec-1,jsc:jec-1,:nk-2) = slope_xz_faceW(:,:,:)
call diagnose_3d(Time, Grid, id_hrm_slopeyW, test_array1(:, :, :))
endif

if (id_delt_x> 0) then
test_array2(:,:) = 0.0
test_array2(isc:iec-1,jsc:jec) = delt_x(:,:)
call diagnose_2d(Time, Grid, id_delt_x, test_array2(:, :))
endif

if (id_delt_y> 0) then
test_array2(:,:) = 0.0
test_array2(isc:iec,jsc:jec-1) = delt_y(:,:)
call diagnose_2d(Time, Grid, id_delt_y, test_array2(:, :))
endif

!if (id_delt_z> 0) then
!test_array3(:,:) = 0.0
!test_array3 = spread(spread(delt_z,2,jec),3,iec)
!call diagnose_3d(Time, Grid, id_delt_z, test_array3)
!endif
!uy(:,:,:) = 10.0
!if (id_uy> 0) then
!uy(:,:,:) = 10.0
!uy(:,:,:)=UN(:,jsc:jec-2,:)  !-US
!call diagnose_3d(Time, Grid, id_uy, uy(:,:,:))
!endif

!vx(:,:,:) = 10.0
!if (id_vx> 0) then
!vx(:,:,:) = 10.0
!vx(:,:,:)=VE(isc:iec-2,:,:)  !-VW
!call diagnose_3d(Time, Grid, id_vx, vx(:,:,:))
!endif

end subroutine compute_hrm_transport

!#######################################################################
! <SUBROUTINE NAME="compute_hrm_Ttransport">
!
! <DESCRIPTION>
! Diagnose the vertically integrated meridional heat flux from HRM. 
! 
! </DESCRIPTION>
!
subroutine compute_hrm_Ttransport (T_prog, T_trans, tau)

! Potential temperature
type(ocean_prog_tracer_type), dimension(:), intent(in)    :: T_prog
real,dimension(isd:ied,jsd:jed), intent(inout) :: T_trans
integer, intent(in) :: tau

real,dimension(isd:ied,jsd:jed,nk-1) :: del_psi
real,dimension(isd:ied,jsd:jed,nk) :: avrg_psi

real :: cp0,rho0

integer :: i, j, k

cp0 = 3991.9
rho0 = 1000.0

avrg_psi(:,:,1) = 0.0 !Force the top layer to be zero
do k = 2,nk
     do j = jsc,jed
       	do i = isc,ied
        	avrg_psi(i,j,k) = (vhrho_nt_hrm(i,j,k-1) + vhrho_nt_hrm(i,j,k))/2.0
		enddo
	enddo
enddo
avrg_psi(:,:,nk) = 0.0 ! Force the bottom layer to be zero

del_psi = avrg_psi(:,:,2:nk) - avrg_psi(:,:,1:nk-1)

!Vertically integrated meridional heat flux
T_trans = rho0*cp0*sum(del_psi*T_prog(index_temp)%field(:, :, :, tau),DIM = 3)
T_trans = T_trans*10.0**(-15.0) ! Unit PW

end subroutine compute_hrm_Ttransport


!#######################################################################
! <FUNCTION NAME="potential_density">
!
! <DESCRIPTION>
! Compute potential density referenced to some given sea pressure. 
!
! Note that potential density referenced to the surface (i.e., sigma_0)
! has a zero sea pressure, so pressure=0.0 should be the argument
! to the function. 
!
! Note that pressure here is 
! sea pressure = absolute pressure - press_standard  (dbars)
!
! input pressure < 0 is an error, and model is brought down.  
! 
! </DESCRIPTION>
!
  function potential_density (salinity, theta, pressure,T_prog)

    !real, dimension(isd:,jsd:,:), intent(in) :: salinity
    !real, dimension(isd:,jsd:,:), intent(in) :: theta
    real, intent(in) :: salinity
    real, intent(in) :: theta
    real,                         intent(in) :: pressure
    type(ocean_prog_tracer_type), intent(in)  :: T_prog(:)

    !real, dimension(isd:ied,jsd:jed,nk) :: potential_density
	real    :: potential_density
    integer :: i, j, k
    real    :: t1, t2, s1, sp5, p1, p1t1
    real    :: num, den,n
      integer :: num_prog_tracers
  integer :: salt_variable
integer :: temp_variable
    
    ! polynomial coefficients in the preTEOS10 EOS
real :: a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11      
real :: b0, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12
! polynomial coefficients in the TEOS10 eos 
real :: v01,v02,v03,v04,v05,v06,v07,v08,v09,v10
real :: v11,v12,v13,v14,v15,v16,v17,v18,v19,v20
real :: v21,v22,v23,v24,v25,v26,v27,v28,v29,v30
real :: v31,v32,v33,v34,v35,v36,v37,v38,v39,v40
real :: v41,v42,v43,v44,v45,v46,v47,v48

num_prog_tracers = size(T_prog(:))
    do n=1,num_prog_tracers
       !if (T_prog(n)%name == 'temp')   index_temp   = n
       !if (T_prog(n)%name == 'salt')   index_salt   = n
       !if (T_prog(n)%name == 'fdelta') index_Fdelta = n
       if (T_prog(n)%longname == 'Conservative temperature') temp_variable = CONSERVATIVE_TEMP
       if (T_prog(n)%longname == 'Potential temperature')    temp_variable = POTENTIAL_TEMP
       if (T_prog(n)%longname == 'Practical Salinity')       salt_variable = PRACTICAL_SALT
       if (T_prog(n)%longname == 'Preformed Salinity')       salt_variable = PREFORMED_SALT
    enddo
    
	v01 =  9.998420897506056d+2
    v02 =  2.839940833161907
    v03 = -3.147759265588511d-2
    v04 =  1.181805545074306d-3
    v05 = -6.698001071123802
    v06 = -2.986498947203215d-2
    v07 =  2.327859407479162d-4
    v08 = -3.988822378968490d-2
    v09 =  5.095422573880500d-4
    v10 = -1.426984671633621d-5
    v11 =  1.645039373682922d-7
    v12 = -2.233269627352527d-2
    v13 = -3.436090079851880d-4
    v14 =  3.726050720345733d-6
    v15 = -1.806789763745328d-4
    v16 =  6.876837219536232d-7
    v17 = -3.087032500374211d-7
    v18 = -1.988366587925593d-8
    v19 = -1.061519070296458d-11
    v20 =  1.550932729220080d-10
    v21 =  1.0
    v22 =  2.775927747785646d-3
    v23 = -2.349607444135925d-5
    v24 =  1.119513357486743d-6
    v25 =  6.743689325042773d-10
    v26 = -7.521448093615448d-3
    v27 = -2.764306979894411d-5
    v28 =  1.262937315098546d-7
    v29 =  9.527875081696435d-10
    v30 = -1.811147201949891d-11
    v31 = -3.303308871386421d-5
    v32 =  3.801564588876298d-7
    v33 = -7.672876869259043d-9
    v34 = -4.634182341116144d-11
    v35 =  2.681097235569143d-12
    v36 =  5.419326551148740d-6
    v37 = -2.742185394906099d-5
    v38 = -3.212746477974189d-7
    v39 =  3.191413910561627d-9
    v40 = -1.931012931541776d-12
    v41 = -1.105097577149576d-7
    v42 =  6.211426728363857d-10
    v43 = -1.119011592875110d-10
    v44 = -1.941660213148725d-11
    v45 = -1.864826425365600d-14
    v46 =  1.119522344879478d-14
    v47 = -1.200507748551599d-15
    v48 =  6.057902487546866d-17
! 25 coefficients in the preTEOS10 equation of state 
    if(temp_variable==CONSERVATIVE_TEMP) then

        !jmfwg_rho   = 1017.842890411975d0
        !jmfwg_alpha = 2.436057013634649d-4
        !jmfwg_beta  = 7.314818108935248d-4

        a0  =  9.9983912878771446d+02
        a1  =  7.0687133522652896d+00
        a2  = -2.2746841916232965d-02
        a3  =  5.6569114861400121d-04
        a4  =  2.3849975952593345d+00
        a5  =  3.1761924314867009d-04
        a6  =  1.7459053010547962d-03
        a7  =  1.2192536310173776d-02
        a8  =  2.4643435731663949d-07
        a9  =  4.0525405332794888d-06
        a10 = -2.3890831309113187d-08
        a11 = -5.9016182471196891d-12

        b0  =  1.0000000000000000d+00 
        b1  =  7.0051665739672298d-03
        b2  = -1.5040804107377016d-05 
        b3  =  5.3943915288426715d-07
        b4  =  3.3811600427083414d-10
        b5  =  1.5599507046153769d-03
        b6  = -1.8137352466500517d-06
        b7  = -3.3580158763335367d-10
        b8  =  5.7149997597561099d-06
        b9  =  7.8025873978107375d-10
        b10 =  7.1038052872522844d-06
        b11 = -2.1692301739460094d-17
        b12 = -8.2564080016458560d-18

        ! Coefficients for neutral density based on McDougall/Jackett (2005).
        ! To be replaced by Klocker/McDougall approach in near future. 
        !rho_neutralrho=1024.43863927763d0
        
        !a0n =  1.0022048243661291d+003
        !a1n =  2.0634684367767725d-001
        !a2n =  8.0483030880783291d-002
        !a3n = -3.6670094757260206d-004
        !a4n = -1.4602011474139313d-003
        !a5n = -2.5860953752447594d-003
        !a6n = -3.0498135030851449d-007

        !b0n =  1.0000000000000000d+000 
        !b1n =  4.4946117492521496d-005
        !b2n =  7.9275128750339643d-005
        !b3n = -1.2358702241599250d-007
        !b4n = -4.1775515358142458d-009
        !b5n = -4.3024523119324234d-004
        !b6n =  6.3377762448794933d-006
        !b7n = -7.2640466666916413d-010
        !b8n = -5.1075068249838284d-005
        !b9n = -5.8104725917890170d-009


    elseif(temp_variable==POTENTIAL_TEMP) then 
    
        !jmfwg_rho   =  1017.728868019642d0
        !jmfwg_alpha = 2.525481286927133d-4
        !jmfwg_beta  = 7.379638527217575d-4

        a0  =  9.9984085444849347d+02
        a1  =  7.3471625860981584d+00
        a2  = -5.3211231792841769d-02
        a3  =  3.6492439109814549d-04
        a4  =  2.5880571023991390d+00
        a5  = -6.7168282786692355d-03
        a6  =  1.9203202055760151d-03
        a7  =  1.1798263740430364d-02
        a8  =  9.8920219266399117d-08
        a9  =  4.6996642771754730d-06
        a10 = -2.5862187075154352d-08
        a11 = -3.2921414007960662d-12

        b0  =  1.0000000000000000d+00 
        b1  =  7.2815210113327091d-03
        b2  = -4.4787265461983921d-05 
        b3  =  3.3851002965802430d-07
        b4  =  1.3651202389758572d-10
        b5  =  1.7632126669040377d-03
        b6  = -8.8066583251206474d-06
        b7  = -1.8832689434804897d-10
        b8  =  5.7463776745432097d-06
        b9  =  1.4716275472242334d-09
        b10 =  6.7103246285651894d-06
        b11 = -2.4461698007024582d-17
        b12 = -9.1534417604289062d-18

        ! Coefficients for neutral density based on McDougall/Jackett (2005).
        ! To be replaced by Klocker/McDougall approach in near future. 
        !rho_neutralrho=1024.59416751197d0

        !a0n =  1.0023063688892480d+003
        !a1n =  2.2280832068441331d-001
        !a2n =  8.1157118782170051d-002
        !a3n = -4.3159255086706703d-004
        !a4n = -1.0304537539692924d-004
        !a5n = -3.1710675488863952d-003
        !a6n = -1.7052298331414675d-007

        !b0n =  1.0000000000000000d+000 
        !b1n =  4.3907692647825900d-005
        !b2n =  7.8717799560577725d-005
        !b3n = -1.6212552470310961d-007
        !b4n = -2.3850178558212048d-009
        !b5n = -5.1268124398160734d-004
        !b6n =  6.0399864718597388d-006
        !b7n = -2.2744455733317707d-009
        !b8n = -3.6138532339703262d-005
        !b9n = -1.3409379420216683d-009
       
    endif 
    if ( .not. module_is_initialized ) then 
      call mpp_error(FATAL, &
      '==>Error in ocean_hrm_mod (potential_density): module must be initialized')
    endif 

    if(pressure < 0.0) then 
       call mpp_error(FATAL, &
       '==>Error in ocean_hrm_mod: potential density at pressure < 0 is not defined')
    endif 


    if(eos_linear) then

        !do k=1,nk
         !  do j=jsd,jed
          !    do i=isd,ied
           !      potential_density(i,j,k) = rho0 - alpha_linear_eos*theta(i,j,k) + beta_linear_eos*salinity(i,j,k)
            !  enddo
           !enddo
        !enddo
    potential_density = rho0 - alpha_linear_eos*theta + beta_linear_eos*salinity
    elseif(eos_preteos10) then 

        p1 = pressure
        if(p1 > 0.0) then 
					 t1  = theta
                     t2  = t1*t1

                     s1  = salinity
                     sp5 = sqrt(s1) 

                     p1t1 = p1*t1

                     num = a0 + t1*(a1 + t1*(a2+a3*t1) )    &
                          + s1*(a4 + a5*t1  + a6*s1)        &
                          + p1*(a7 + a8*t2 + a9*s1 + p1*(a10+a11*t2))

                     den = b0 + t1*(b1 + t1*(b2 + t1*(b3 + t1*b4)))      &
                          + s1*(b5 + t1*(b6 + b7*t2) + sp5*(b8 + b9*t2)) &  
                          + p1*(b10 + p1t1*(b11*t2 + b12*p1))

                     potential_density = num/(epsln+den)
                     
            !do k=1,nk
             !  do j=jsd,jed
             !     do i=isd,ied
              !       t1  = theta(i,j,k)
               !      t2  = t1*t1

                 !    s1  = salinity(i,j,k)
                  !   sp5 = sqrt(s1) 

!                     p1t1 = p1*t1

 !                    num = a0 + t1*(a1 + t1*(a2+a3*t1) )    &
  !                        + s1*(a4 + a5*t1  + a6*s1)        &
  !                        + p1*(a7 + a8*t2 + a9*s1 + p1*(a10+a11*t2))

   !                  den = b0 + t1*(b1 + t1*(b2 + t1*(b3 + t1*b4)))      &
!                          + s1*(b5 + t1*(b6 + b7*t2) + sp5*(b8 + b9*t2)) &  
 !                         + p1*(b10 + p1t1*(b11*t2 + b12*p1))

!                     potential_density(i,j,k) = num/(epsln+den)
 !                 enddo
  !             enddo
   !         enddo

        elseif(p1==0.0) then ! simplification for sigma_0
        t1  = theta
                     t2  = t1*t1

                     s1  = salinity
                     sp5 = sqrt(s1) 

                     num = a0 + t1*(a1 + t1*(a2+a3*t1) )  &
                          + s1*(a4 + a5*t1  + a6*s1)        

                     den = b0 + t1*(b1 + t1*(b2 + t1*(b3 + t1*b4)))      &
                          + s1*(b5 + t1*(b6 + b7*t2) + sp5*(b8 + b9*t2)) 

                     potential_density = num/(epsln+den)

            !do k=1,nk
             !  do j=jsd,jed
              !    do i=isd,ied
               !      t1  = theta(i,j,k)
                !     t2  = t1*t1

                 !    s1  = salinity(i,j,k)
                  !   sp5 = sqrt(s1) 

                   !  num = a0 + t1*(a1 + t1*(a2+a3*t1) )  &
                    !      + s1*(a4 + a5*t1  + a6*s1)        

                     !den = b0 + t1*(b1 + t1*(b2 + t1*(b3 + t1*b4)))      &
                      !    + s1*(b5 + t1*(b6 + b7*t2) + sp5*(b8 + b9*t2)) 

                     !potential_density(i,j,k) = num/(epsln+den)
                  !enddo
               !enddo
            !enddo

        endif  ! endif for value of pressure 

    else   ! eos_teos10

        p1 = pressure
        if(p1 > 0.0) then 
     				 t1  = theta
                     t2  = t1*t1

                     s1  = salinity
                     sp5 = sqrt(s1) 

                     p1t1 = p1*t1

                     num = v01 + t1*(v02 + t1*(v03 + v04*t1))               &
                          + s1*(v05 + t1*(v06 + v07*t1)                     &
                          + sp5*(v08 + t1*(v09 + t1*(v10 + v11*t1))))       &
                          + p1*(v12 + t1*(v13 + v14*t1) + s1*(v15 + v16*t1) &
                          + p1*(v17 + t1*(v18 + v19*t1) + v20*s1))

                     den =  v21 + t1*(v22 + t1*(v23 + t1*(v24 + v25*t1)))              &
                          + s1*(v26 + t1*(v27 + t1*(v28 + t1*(v29 + v30*t1))) + v36*s1 &
                          + sp5*(v31 + t1*(v32 + t1*(v33 + t1*(v34 + v35*t1)))))       &
                          + p1*(v37 + t1*(v38 + t1*(v39 + v40*t1))                     &
                          +     s1*(v41 + v42*t1)                                      &
                          +     p1*(v43 + t1*(v44 + v45*t1 + v46*s1)                   &
                          +     p1*(v47 + v48*t1)))

                     potential_density = num/(epsln+den)
           ! do k=1,nk
            !   do j=jsd,jed
             !     do i=isd,ied
              !       t1  = theta(i,j,k)
               !      t2  = t1*t1

                !     s1  = salinity(i,j,k)
                 !    sp5 = sqrt(s1) 

                  !   p1t1 = p1*t1

                   !  num = v01 + t1*(v02 + t1*(v03 + v04*t1))               &
                    !      + s1*(v05 + t1*(v06 + v07*t1)                     &
                     !     + sp5*(v08 + t1*(v09 + t1*(v10 + v11*t1))))       &
                      !    + p1*(v12 + t1*(v13 + v14*t1) + s1*(v15 + v16*t1) &
                       !   + p1*(v17 + t1*(v18 + v19*t1) + v20*s1))

                     !den =  v21 + t1*(v22 + t1*(v23 + t1*(v24 + v25*t1)))              &
                      !    + s1*(v26 + t1*(v27 + t1*(v28 + t1*(v29 + v30*t1))) + v36*s1 &
                       !   + sp5*(v31 + t1*(v32 + t1*(v33 + t1*(v34 + v35*t1)))))       &
                        !  + p1*(v37 + t1*(v38 + t1*(v39 + v40*t1))                     &
                         ! +     s1*(v41 + v42*t1)                                      &
                          !+     p1*(v43 + t1*(v44 + v45*t1 + v46*s1)                   &
                          !+     p1*(v47 + v48*t1)))

              !       potential_density(i,j,k) = num/(epsln+den)
             !     enddo
            !   enddo
           ! enddo

        elseif(p1==0.0) then ! for sigma_0
                     t1  = theta
                     t2  = t1*t1

                     s1  = salinity
                     sp5 = sqrt(s1) 

                     num = v01 + t1*(v02 + t1*(v03 + v04*t1))         &
                          + s1*(v05 + t1*(v06 + v07*t1)               &
                          + sp5*(v08 + t1*(v09 + t1*(v10 + v11*t1)))) 

                     den =  v21 + t1*(v22 + t1*(v23 + t1*(v24 + v25*t1)))              &
                          + s1*(v26 + t1*(v27 + t1*(v28 + t1*(v29 + v30*t1))) + v36*s1 &
                          + sp5*(v31 + t1*(v32 + t1*(v33 + t1*(v34 + v35*t1)))))


                     potential_density = num/(epsln+den)
            !do k=1,nk
             !  do j=jsd,jed
              !    do i=isd,ied
               !      t1  = theta(i,j,k)
                !     t2  = t1*t1

!                     s1  = salinity(i,j,k)
 !                    sp5 = sqrt(s1) 

  !                   num = v01 + t1*(v02 + t1*(v03 + v04*t1))         &
   !                       + s1*(v05 + t1*(v06 + v07*t1)               &
    !                      + sp5*(v08 + t1*(v09 + t1*(v10 + v11*t1)))) 

     !                den =  v21 + t1*(v22 + t1*(v23 + t1*(v24 + v25*t1)))              &
      !                    + s1*(v26 + t1*(v27 + t1*(v28 + t1*(v29 + v30*t1))) + v36*s1 &
       !                   + sp5*(v31 + t1*(v32 + t1*(v33 + t1*(v34 + v35*t1)))))


        !             potential_density(i,j,k) = num/(epsln+den)
         !         enddo
          !     enddo
           ! enddo

        endif  ! endif for pressure 

    endif  ! endif for eos choice 


  end function potential_density
! </FUNCTION> NAME="potential_density"

!#######################################################################
! <SUBROUTINE NAME="ocean_hrm_end">
!
! <DESCRIPTION>
! Write to restart.
! </DESCRIPTION>
!
subroutine ocean_hrm_end(Time)

  type(ocean_time_type), intent(in) :: Time

  if(.not. use_this_module) return

  deallocate(vhrho_nt_hrm)
  deallocate(uhrho_et_hrm)
  deallocate(T_trans)

end subroutine ocean_hrm_end
! </SUBROUTINE> NAME="ocean_hrm_end"


end module ocean_hrm_mod
