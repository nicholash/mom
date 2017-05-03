module ocean_hrm_mod
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
!</OVERVIEW>
!
!<DESCRIPTION>
!</DESCRIPTION>
!
! <INFO>
!
! <REFERENCE>
! </REFERENCE>
!
! <REFERENCE>
! </REFERENCE>
!
! <REFERENCE>
! </REFERENCE>
!
! <NOTE>
! </NOTE> 
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
!
!</NAMELIST>

use diag_manager_mod,        only: register_diag_field, register_static_field, need_data, send_data
use fms_mod,                 only: FATAL, WARNING, NOTE
use fms_mod,                 only: open_namelist_file, check_nml_error, close_file, write_version_number
use mpp_domains_mod,         only: mpp_update_domains
use mpp_domains_mod,         only: CGRID_NE, WUPDATE, SUPDATE, EUPDATE, NUPDATE
use ocean_parameters_mod,    only: missing_value
use mpp_domains_mod,         only: cyclic_global_domain, global_data_domain
use mpp_mod,                 only: input_nml_file, mpp_error, stdout, stdlog
use mpp_mod,                 only: mpp_clock_id, mpp_clock_begin, mpp_clock_end, CLOCK_ROUTINE

use ocean_nphysics_new_mod,  only: neutral_slopes, gradrho, tracer_gradients
use ocean_tracer_diag_mod,   only: diagnose_eta_tend_3dflux
use ocean_types_mod,         only: ocean_grid_type, ocean_domain_type, ocean_density_type
use ocean_types_mod,         only: ocean_prog_tracer_type, ocean_thickness_type
use ocean_types_mod,         only: ocean_velocity_type
use ocean_types_mod,         only: ocean_time_type, ocean_time_steps_type
use ocean_types_mod,         only: tracer_3d_0_nk_type, tracer_3d_1_nk_type
use ocean_util_mod,          only: diagnose_2d, diagnose_3d

implicit none

public ocean_hrm_init
public compute_hrm_transport
public ocean_hrm_end

private 

type(ocean_grid_type), pointer   :: Grd => NULL()
type(ocean_domain_type), pointer :: Dom => NULL()
type(ocean_domain_type), save    :: Dom_flux

integer :: id_uhrho_et_hrm      =-1
integer :: id_vhrho_nt_hrm      =-1

#include <ocean_memory.h>

#ifdef MOM_STATIC_ARRAYS

#else

#endif

! for diagnostics 
real, dimension(:,:,:), allocatable :: uhrho_et_hrm  ! i-component of advective mass transport
real, dimension(:,:,:), allocatable :: vhrho_nt_hrm  ! j-component of advective mass transport
real, dimension(:,:,:), allocatable :: T_trans       ! vertical sum of T transport

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

! for the module as a whole 
logical :: use_this_module   = .false.
logical :: debug_this_module = .false.

!**************end of nml settings**************

namelist /ocean_hrm_nml/ use_this_module, debug_this_module

contains

!#######################################################################
! <SUBROUTINE NAME="ocean_hrm_init">
!
! <DESCRIPTION>
! Initialize the neutral physics module by registering fields for 
! diagnostic output and performing some numerical checks to see 
! that namelist settings are appropriate. 
! </DESCRIPTION>
!
subroutine ocean_hrm_init(Grid, Domain, Time, Time_steps, Thickness, Dens, T_prog,&
           velocity)

  type(ocean_grid_type),        intent(in), target   :: Grid
  type(ocean_domain_type),      intent(in), target   :: Domain
  type(ocean_time_type),        intent(in)           :: Time
  type(ocean_time_steps_type),  intent(in)           :: Time_steps
  type(ocean_thickness_type),   intent(in)           :: Thickness
  type(ocean_density_type),     intent(in)           :: Dens
  type(ocean_prog_tracer_type), intent(inout)        :: T_prog(:)
  type(ocean_velocity_type),    intent(inout)        :: Velocity

  allocate (vhrho_nt_hrm(isd:ied,jsd:jed,nk))
  allocate (uhrho_et_hrm(isd:ied,jsd:jed,nk))

  id_uhrho_et_hrm = register_diag_field ('ocean_model', 'uhrho_et_hrm', &
        Grd%tracer_axes_flux_x(1:3), Time%model_time,                  &
       'i-component of hrm mass transport',                            &
       '(kg/m^3)*m^2/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))

  id_vhrho_nt_hrm = register_diag_field ('ocean_model', 'vhrho_nt_hrm', &
        Grd%tracer_axes_flux_y(1:3), Time%model_time,                  &
       'j-component of hrm mass transport',                            &
       '(kg/m^3)*m^2/sec', missing_value=missing_value, range=(/-1.e10,1.e10/))

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
subroutine compute_hrm_transport(Time, Dens, Thickness, Grid, T_prog, Velocity)

  type(ocean_time_type),      intent(in) :: Time
  type(ocean_density_type),   intent(in) :: Dens
  type(ocean_thickness_type), intent(in) :: Thickness
  type(ocean_grid_type),      intent(in) :: Grid
  type(ocean_prog_tracer_type), intent(in)        :: T_prog(:)
  type(ocean_velocity_type),  intent(in) :: Velocity

  real, dimension(isd:ied,jsd:jed,nk,size(T_prog)) :: dTdx, dTdy
  real, dimension(isd:ied,jsd:jed,0:nk,size(T_prog)) :: dTdz

  real, dimension(isd:ied,jsd:jed,nk,0:1,0:1) :: slopex_xz
  real, dimension(isd:ied,jsd:jed,nk,0:1,0:1) :: slopey_yz

  ! neutral density derivatives 
  real, dimension(isd:ied,jsd:jed,nk,0:1) :: drhodx_x
  real, dimension(isd:ied,jsd:jed,nk,0:1) :: drhody_y
  real, dimension(isd:ied,jsd:jed,nk,0:1) :: drhodz_z 

  real, dimension(isd:ied,jsd:jed,nk,0:1) :: absslope_z
  real, dimension(isd:ied,jsd:jed,nk) :: absslope

  real,dimension(isd:ied,jsd:jed,nk) :: vE,vW,uN,uS
  real,dimension(isd:ied+1,jsd:jed,nk-1) :: vu,vl
  real,dimension(isd:ied,jsd:jed+1,nk-1) :: uu,ul
  real,dimension(nk) :: z, delt_z
  real,dimension(isd:ied,jsd:jed) :: delt_x
  real,dimension(isd:ied, jsd:jed, nk, 0:1) :: slope_xz_avg, slope_yz_avg
  real,dimension(isd:ied+1, jsd:jed, nk) :: slope_xz_zonal,slope_xz_face,v
  real,dimension(isd:ied, jsd:jed+1, nk) :: slope_yz_meridional,slope_yz_face,u

  integer :: tau

  integer :: i, j, k

  integer :: stdoutunit
  stdoutunit=stdout()

  print*, 'CALLING compute_hrm_transport'

  tau = Time%tau

  call tracer_gradients(Time, Grid, Time%taum1, T_prog, Thickness, &
         dTdx, dTdy, dTdz)

  ! Calculate density gradients
  call gradrho(Dens, dTdx, dTdy, dTdz, drhodx_x, drhody_y, drhodz_z)

  ! Calculate neutral slopes
  call neutral_slopes(Time, Grid, drhodx_x, drhody_y, drhodz_z, &
         slopex_xz, slopey_yz, absslope_z, absslope)

  ! i-component of hrm mass transport. units ((kg/m^3)*m^2/sec)
  uhrho_et_hrm(:,:,:) = 0.0
  ! j-component of hrm mass transport. units:kg/sec     units ((kg/m^3)*m^2/sec)
  vhrho_nt_hrm(:,:,:) = 0.0

  v(2:,:,:) = velocity%u(:,:,:,2,tau) !The j component of velocity at t time level
  v(1,:,:) = 0.0
  ! ??? Is velocity located at the northeastern lower vertex of a tracer cell? Need to check the index for following operations
  vE = v(2:,:,:) !The velocity on the eastern side of a box
  vW = v(1:ied-1,:,:) !The velocity on the western side of a box
  VE(:,:,2:) = (vE(:,:,1:ied-1) + vE(:,:,2:))/2 !Average the velocity onto depth z
  VW(:,:,2:) = (vW(:,:,1:ied-1) + vW(:,:,2:))/2 !Average the velocity onto depth z
  vu = v(:,:,1:ied-1) !Velocity at upper face
  vl = v(:,:,2:)   !Velocity at lower face
  Vu(1:ied-1,:,:) = (vu(1:ied-1,:,:) + vu(2:,:,:))/2 !Average to the center
  Vl(1:ied-1,:,:) = (vl(1:ied-1,:,:) + vl(2:,:,:))/2 !Average to the center
  
  u(:,2:,:) = velocity%u(:,:,:,1,tau) !The i component of velocity at t time level
  u(:,1,:) = 0.0
  uN = u(:,2:,:) !The velocity on the northern side of a box
  uS = v(:,1:jed-1,:) !The velocity on the southern side of a box
  UN(:,:,2:) = (uN(:,:,1:ied-1) + uN(:,:,2:))/2 !Average the velocity onto depth z
  US(:,:,2:) = (uS(:,:,1:ied-1) + uS(:,:,2:))/2 !Average the velocity onto depth z
  uu = u(:,:,1:ied-1) !Velocity at upper face
  ul = u(:,:,2:)   !Velocity at lower face
  Uu(:,1:jed-1,:) = (vu(:,1:jed-1,:) + vu(:,2:,:))/2 !Average to the center
  Ul(:,1:jed-1,:) = (vl(:,1:jed-1,:) + vl(:,2:,:))/2 !Average to the center
  
  z = Grid%zt  ! distance from surface to grid point in level k (m) 
  delt_z = Grid%dzt  ! initial vertical resolution of T or U grid cells (m)
  delt_x = Grid%dxtn ! i-width of north face of T-cells (m)

  
  !Take the average of the upper and lower quater cells
  slope_xz_avg = (slopex_xz(:,:,:,:,0) + slopex_xz(:,:,:,:,1))/2.0
  
  slope_xz_zonal(2:ied,:,:) = (slope_xz_avg(:ied-1,:,:,1) + slope_xz_avg(2:,:,:,0))/2.0
  slope_xz_zonal(1,:,:) = slope_xz_avg(1,:,:,0)
  slope_xz_zonal(ied+1,:,:) = slope_xz_avg(ied,:,:,1)
  ! Kept the first face but not the last
  slope_xz_face(:,2:,:) = (slope_xz_zonal(:,2:,:) + slope_xz_zonal(:,1:jed-1,:))/2.0
  slope_xz_face(:,1,:) = slope_xz_zonal(:,1,:)
  
  slope_yz_avg = (slopey_yz(:,:,:,:,0) + slopey_yz(:,:,:,:,1))/2  
  slope_yz_meridional(:,2:jed,:) = (slope_yz_avg(:,:jed-1,:,1) + slope_yz_avg(:,2:,:,0))/2.0
  slope_yz_meridional(:,1,:) = slope_yz_avg(:,1,:,0)
  slope_yz_meridional(:,jed+1,:) = slope_yz_avg(:,jed,:,1)
  slope_yz_face(2:,:,:) = (slope_yz_meridional(2:,:,:) + slope_yz_meridional(1:jed-1,:,:))/2.0
  slope_yz_face(1,:,:) = slope_yz_meridional(1,:,:)
  
  
  !print *,size(slopexz) !To check the dimension of slopexz
  
 !Set a limit to the neutral slopes
do k = 1,size(slope_xz_face,3)
     do j = 1,size(slope_xz_face,2)
       do i = 1,size(slope_xz_face,1)
         if (slope_xz_face(i,j,k) > 10**-2) then
  					slope_xz_face(i,j,k) = 10**-2
  		 end if
  	  enddo
  	enddo
 enddo
 
  !Set a limit to the neutral slopes
do k = 1,size(slope_yz_face,3)
     do j = 1,size(slope_yz_face,2)
       do i = 1,size(slope_yz_face,1)
         if (slope_yz_face(i,j,k) > 10**-2) then
  					slope_yz_face(i,j,k) = 10**-2
  		 end if
  	  enddo
  	enddo
 enddo

  !Calculate HRM from second depths
  do k=2,nk
     do j=jsc,jec
        do i=isc,iec
           uhrho_et_hrm(i,j,k) = 1/24*(UN(i, j, k)-US(i, j, k))*(slope_yz_face(i,j,k)+slope_yz_face(i,j+1,k))*(delt_x(i,j))**2 + &
           1/48*(Uu(i, j, k-1)-Ul(i, j, k-1))/delt_z(k)*((slope_yz_face(i,j,k))**2+(slope_yz_face(i,j+1,k))**2)*(delt_x(i, j))**3

           vhrho_nt_hrm(i,j,k) = 1/24*(VE(i, j, k)-VW(i, j, k))*(slope_xz_face(i,j,k)+slope_xz_face(i+1,j,k))*(delt_x(i,j))**2 + &
           1/48*(Vu(i, j, k-1)-Vl(i, j, k-1))/delt_z(k)*((slope_xz_face(i,j,k))**2+(slope_xz_face(i+1,j,k))**2)*(delt_x(i, j))**3
        enddo
     enddo
  enddo

  call compute_hrm_Ttransport(T_prog, T_trans, tau)

  call diagnose_3d(Time, Grd, id_uhrho_et_hrm, uhrho_et_hrm(:,:,:))
  call diagnose_3d(Time, Grd, id_vhrho_nt_hrm, vhrho_nt_hrm(:,:,:))

end subroutine compute_hrm_transport

!#######################################################################
! <SUBROUTINE NAME="compute_hrm_Ttransport">
!
! <DESCRIPTION>
! Diagnose advective potential temperature transport from HRM. 
!
subroutine compute_hrm_Ttransport (T_prog, T_trans, tau)

! Potential temperature
type(ocean_prog_tracer_type), dimension(:), intent(in)    :: T_prog
real,dimension(isd:ied,jsd:jed), intent(inout) :: T_trans
integer, intent(in) :: tau

real,dimension(isd:ied,jsd:jed,nk) :: del_psi
real,dimension(isd:ied,jsd:jed,nk+1) :: avrg_psi
integer :: num_prog_tracers, n, index_temp

integer :: i, j, k

num_prog_tracers = size(T_prog)

! read the the potential temperature
index_temp = -1
do n=1,num_prog_tracers
    if (T_prog(n)%name == 'temp') then
        index_temp = n
    endif
enddo
if (index_temp == -1) then
    call mpp_error(FATAL, &
            '==>Error: temp not identified in call to compute_hrm_Ttransport')
endif
avrg_psi(:,:,1) = 0.0
   do k = 2,nk
     do j = jsd,jed
       do i = isd,ied
        avrg_psi(i,j,k) = (vhrho_nt_hrm(i,j,k) + vhrho_nt_hrm(i,j,k))/2
enddo
enddo
enddo
avrg_psi(:,:,nk+1) = 0.0

del_psi = avrg_psi(:,:,2:nk) - avrg_psi(:,:,1:nk-1)

T_trans = sum(del_psi*T_prog(index_temp)%field(:, :, :, tau),DIM = 3)

end subroutine compute_hrm_Ttransport


! </SUBROUTINE> NAME="compute_advect_transport"


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

  ! FIXME I don't know if we need this.

end subroutine ocean_hrm_end
! </SUBROUTINE> NAME="ocean_hrm_end"


end module ocean_hrm_mod
