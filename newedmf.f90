! Modified EDMF
module newedmf_mod

use           fms_mod, only: mpp_pe, mpp_root_pe, file_exist, error_mesg, open_file, &
                             read_data, write_data, write_version_number, check_nml_error, & 
                             open_namelist_file, stdlog, FATAL, close_file, nullify_domain
                             
use           time_manager_mod, only: time_type
use            diffusivity_mod, only: diffusivity, pbl_depth
use              vert_diff_mod, only: gcm_vert_diff_down_ed_new, gcm_vert_diff_up_ed_new, & 
                                      gcm_vert_diff_tke, surf_diff_type
use            mixed_layer_mod, only: mixed_layer, mixed_layer_noupdate
use  simple_sat_vapor_pres_mod, only: escomp, descomp
use              constants_mod, only: HLv,Cp_air,Grav,rdgas,rvgas,kappa,pstd_mks
use           diag_manager_mod, only: register_diag_field, send_data
 
use             transforms_mod, only: grid_domain  ! Added ZTAN 04/08/2012 for new version

implicit none
private

public newedmf_init, newedmf, newedmf_end

!-----------------------------------------------------------------------

character(len=128) :: version = '$Id: newedmf.F90,v 0.3.0 2018/04/10 $'
character(len=128) :: tagname = '$Name:  $'
character(len=10), parameter :: mod_name='newedmf'

integer :: id_thl_e, id_qt_e, id_ql_e, id_tke_e, id_a_e, & 
           id_thl_u, id_qt_u, id_ql_u, id_w_u, id_a_u, id_m_u, id_b_u, &
           id_k_t, id_k_m, id_ls, id_ustar, id_wstar, id_zstar, id_bstar, &
           id_flux_thl_ed, id_flux_qt_ed, id_flux_ql_ed, id_flux_tke_ed, &
           id_flux_thl_mf, id_flux_qt_mf, id_flux_ql_mf, id_flux_tke_mf, &
           id_rh, id_s, id_sigmas, id_ccov_ed, id_ccov_mf, &
           id_precip_ed, id_precip_mf, id_rain2d_ed, id_rain2d_mf, &
           id_entT, id_detT, &
           id_tend_tke_buoy, id_tend_tke_shear, id_tend_tke_detr, id_tend_tke_entr, &
           id_tend_tke_diss, id_tend_tke_corr, id_tend_tke_mf, id_tend_tke_ed, id_tend_tke_tot

integer :: id_t_g_sl, id_q_g_sl, id_ql_g_sl, id_u_g_sl, id_v_g_sl, &
           id_thl_g_sl, id_qt_g_sl, id_thv_g_sl, id_rho_g_sl, id_tke_g_sl, &
           id_zfull_sl, id_zhalf_sl, id_pfull_sl, id_phalf_sl, id_rh_sl, &
           id_ccov_ed_sl, id_ccov_mf_sl, id_rain2d_ed_sl, id_rain2d_mf_sl, &
           id_flux_t_sl, id_flux_q_sl, &  ! Grid-mean output above
           id_thl_e_sl, id_qt_e_sl, id_ql_e_sl, id_tke_e_sl, id_a_e_sl, & 
           id_thl_u_sl, id_qt_u_sl, id_ql_u_sl, id_w_u_sl, id_a_u_sl, id_m_u_sl, id_b_u_sl, &
           id_entr_u_sl, id_detr_u_sl, id_entT_sl, id_detT_sl, &
           id_k_t_sl, id_k_m_sl, id_ls_sl, id_ustar_sl, id_wstar_sl, id_zstar_sl, id_bstar_sl, &
           id_flux_thl_ed_sl, id_flux_qt_ed_sl, id_flux_ql_ed_sl, &
           id_flux_thl_mf_sl, id_flux_qt_mf_sl, id_flux_ql_mf_sl
           
! Debug only: slice output for tendencies
integer :: id_thldel_sl, id_qtdel_sl, id_tdel_sl, id_qdel_sl, id_liqdel_sl, &  
           id_udel_sl, id_vdel_sl, id_tkedel_sl, &
           id_thldel_diag_sl, id_qtdel_diag_sl, id_tdel_diag_sl, id_qdel_diag_sl, id_liqdel_diag_sl, &
           id_udel_diag_sl, id_vdel_diag_sl, id_tkedel_diag_sl, &
           id_tend_tke_buoy_sl, id_tend_tke_shear_sl, id_tend_tke_detr_sl, id_tend_tke_entr_sl, &
           id_tend_tke_diss_sl, id_tend_tke_corr_sl, id_tend_tke_mf_sl, id_tend_tke_ed_sl, id_tend_tke_tot_sl
           
! --- local constants --- !
real, parameter :: d622 = rdgas/rvgas
real, parameter :: d378 = 1.-d622
real, parameter :: d608 = (rvgas-rdgas)/rdgas
real, parameter :: ce   = 0.16 ! "Turbulent mixing coefficient"
real, parameter :: vonk = 0.41 ! Von Karman's Constant
real, parameter :: small_val = 1.e-10

logical :: updraft_initialized = .false.

! --- local options --- !
integer :: maxiter          = 2     ! Maximal iteration times for saturation adjustment
real    :: updraft_rescale  = 1.0   ! Rescaling of updraft prognostic term (1 = fully prognostic, <<1 = diagnostic)
integer :: updraft_number   = 1     ! Number of updrafts (multiple updrafts may not work yet) 
real    :: updraft_fraction = 0.1   ! Surface updraft fraction
real    :: updraft_exponent = 1.0   ! Scaling factor of multiple updrafts (not used for single updraft)

real    :: var_surf_fac     = 1.0   ! rescale the initial condition for updrafts

integer :: pbl_depth_opt    = 1     ! See subroutine compute_zstar_wstar_obl
real    :: pbl_fixed_depth  = 1000. ! fixed PBL depth (m) (only when pbl_depth_opt == 0)
real    :: pbl_depth_res    = 1.0   ! rescaling factor for PBL depth (only for k-profile scalar mixing)

integer :: updr_iter        = 1     ! Iteration for updraft
integer :: correct_buoy     = 12    ! Option for correcting buoyancy in updraft
integer :: ER               = 0     ! Option for entrainment rate
real    :: ER0_zmin         = 1.0   ! lower limit for z (upper limit for entrainment rate) (only when ER == 0)
real    :: c_tau_ER5        = 450.0 ! used in tau/wu entrainment/detrainment only (ER == 5)
real    :: ER_frac          = 1.0   ! Scaling factor for entrainment rate
real    :: ER_KF_min        = 0.0   ! Minimal fraction of KF mixture to be entrained (only when DR == 10)
integer :: DR               = 0     ! Option for detrainment rate
real    :: DR_frac          = 1.0   ! Scaling factor for detrainment rate

integer :: wu_opt           = 1     ! Option for wu scheme
real    :: wu_min           = 0.0   ! lower limit for updraft velocity

integer :: au_optL          = 1     ! Option for lateral entrainment/detrainment scheme
integer :: au_optB          = 1     ! Option for entrainment/detrainment control on updraft fraction
real    :: au_optB1_frac    = 1.0   ! Ratio of max updraft fraction to surface updraft fraction (only when au_optB == 1)
integer :: au_optB_wu       = 1     ! Option for top detrainment
integer :: au_optB_srf      = 1     ! Option for bottom entrainment adjustment

real    :: updr_surfht      = 100.0 ! Below this height, the updraft entrains from the buoyant tail of envr distribution.
                                    ! = 0.0: only entrain from the buoyant tail at the first level

logical :: mf_full_uw         = .true.
logical :: do_mf_env_implicit = .true.  ! Do implicit fix for environment part of MF tendency (updraft part is still explicit)

integer :: precip_upd_opt   = 1
real    :: precip_upd_fac   = 0.02      ! Precip relative threshold (only when precip_upd_opt == 1)
real    :: precip_upd_val   = 0.0005    ! Precip absolute threshold (only when precip_upd_opt == 2)



! Eddy Diffusivity options
integer :: diff_opt         = 2         ! Option for eddy diffusivity 
real    :: diff_min         = 0.0       ! Lower limit on eddy diffusivity
real    :: diff_rescale     = 1.0       ! Rescaling factor for eddy diffusivity
real    :: c_diff           = 0.5       ! Scaling parameter for TKE-based diffusivity
integer :: diff1_fix        = 0         ! Fixing k-prof diffusivity with TKE-based (when stable) (only when diff_opt == 1)

! Mixing Length options
integer :: tau_opt          = 2         ! Option for turbulent overturning timescale
integer :: ml_opt           = 1         ! Option for mixing length
logical :: ml_limiter       = .false.   ! Option for height-dependent limiter of ML
real    :: lmax             = 300.0     ! Upper limit for ML (only when ml_opt == 1 or 2)
real    :: dls_factor       = 2.5       ! dls = ls / dls_factor  ! Added ZTAN 11/17/2017

! Surface Flux Options
logical :: surf_explicit    = .false.   ! Use explicit scheme for updating SST 

! Cloud Scheme Options
integer :: sigma_opt        = 1
logical :: sigma_qonly      = .false.   ! Explicit qs
logical :: sigma_tonly      = .false. 
integer :: precip_env_opt   = 1
real    :: precip_env_fac   = 0.02      ! Precip relative threshold (only when precip_env_opt == 1)
real    :: precip_env_val   = 0.0005    ! Precip absolute threshold (only when precip_env_opt == 2)

! TKE Scheme Options
logical :: buoy_freeatm     = .true.    ! Allow buoyancy production of TKE in the free troposphere
logical :: pbuoy_corr       = .false.   ! Correct for pbuoy in MF term
integer :: ed_tke_opt       = 1         ! 


logical :: do_trans_dphi    = .false.   ! Transform thl_u, qt_u, ql_u to anomalies
                                        ! (thl_u - thl_e), (qt_u - qt_e), (ql_u - ql_e) ?
                                        
integer :: phi_u_corr_opt   = 0         ! Option for do_phiu_corr
                                        ! 0 := No correction, i.e., fill in environment value.
                                        ! 1 := Fill in updraft value from the last timestep.
                                        ! 2 := Fill in phi_u(ilev+1)*.7 + phi_e(ilev)*.3 (empirical).

! Use the new au/wu solver (do both together). If false, do wu first, au thereafter.
logical :: do_wu_first      = .false.
integer :: entr_opt         = 0  ! 0 := implicit entrainment; 
                                 ! 1 := explicit entrainment from last step;  
                                 ! 2 := explicit entrainment from lower level.  
logical :: do_buoy_avg      = .true.  ! Use the average buoyancy from the last step and the current step
logical :: do_compute_updr_lbl = .true.  ! Use the level-by-level updraft solver...

! Updraft variables
real, allocatable, dimension(:,:,:,:) :: a_u, w_u, thl_u, qt_u, ql_u
real, allocatable, dimension(:,:,:,:) :: buoy_o_u ! Added for two-step average... 10/10/2017
real, allocatable, dimension(:) :: init_a_u, init_scl_u 

logical :: debug_point = .false.  ! For debug ZTAN 09/12/2017
integer :: idb, jdb               ! For debug ZTAN 09/12/2017
integer :: sl_ind = 1             ! Which point to output for slide output. (Future: -1 for average)

real    :: gaussian_std           ! Gaussian mean for upper 10% as standard

namelist /newedmf_nml/ maxiter, & 
    updraft_rescale, updraft_number, updraft_fraction, updraft_exponent, & 
    pbl_depth_res, diff_min, diff_rescale, ml_limiter, var_surf_fac, & 
    pbl_depth_opt, pbl_fixed_depth, pbl_depth_res, updr_iter, correct_buoy,   &
    ER, ER0_zmin, ER_frac, c_tau_ER5, DR, ER_KF_min, DR_frac, wu_opt, wu_min, & 
    au_optL, au_optB, au_optB1_frac, au_optB_wu, au_optB_srf, updr_surfht,    &
    mf_full_uw, do_mf_env_implicit, diff_opt, diff_min, diff_rescale, c_diff, diff1_fix, &
    tau_opt, ml_opt, ml_limiter, lmax, surf_explicit, sigma_opt, sigma_qonly, sigma_tonly, &
    precip_upd_opt, precip_upd_fac, precip_upd_val,  &
    precip_env_opt, precip_env_fac, precip_env_val,  &
    buoy_freeatm, pbuoy_corr, ed_tke_opt, do_trans_dphi, phi_u_corr_opt, sl_ind, &
    do_wu_first, entr_opt, do_buoy_avg, dls_factor, do_compute_updr_lbl ! Add more, ed_leap_frog is removed
    
real :: missing_value = -1.e-10

contains

! main subroutine: newedmf

! param subroutines:
! compute_tend_tke, compute_tend_grid, compute_tend_cloud, compute_precip, cal_sigmas_new, 
! compute_tend_ed, compute_eddy_diffusivity, compute_mixing_length, compute_tau, 
! update_surf_flux_mf, update_grid_mean_mf, compute_tend_mf_first, compute_tend_mf_second, 
! compute_tend_mf, initialize_updr, compute_updr, compute_bulk_updraft, compute_phiu, 
! compute_precip_upd, compute_au, compute_wu, compute_upd_buoy, compute_ent_det, 
! compute_zstar_wstar_obl, cal_wstar

! utility subroutines: 
! cloud_mixing, cloud_decompose, compute_thv, compute_thv_checkcond, compute_ql, 
! envupd_decompose, area_sum_updr

! subroutines from other modules:
! vert_diff_mod: gcm_vert_diff_down_ed (compute_mu, compute_nu, uv_vert_diff, 
!    vert_diff_down_1, vert_diff_down_2, vert_diff_down_3, diff_surface, vert_diff_up, 
!    explicit_tend, compute_e, compute_f, gcm_vert_diff_up_ed_new)
! diffusivity_mod: pbl_depth

! ZTAN 11/23/2017: Starting to rewrite the compute_updr subroutine into a level-by-level version.


subroutine newedmf_init(idim, jdim, kdim, axes, Time, cloud_cover, mfcloud_cover)
! TO_DO 082016: Add more fields (au, thlu, qu, ...)? 

!                         ---  input arguments ---

      integer, intent(in) :: axes(4)
      integer, intent(in) :: idim, jdim, kdim
      integer, dimension(3) :: half = (/1,2,4/)
      integer, dimension(2) :: half_slice = (/2,4/)
      integer, dimension(1) :: sfc_slice = (/2/)
      type(time_type), intent(in) :: Time
      real, intent(out), dimension(:,:,:) :: cloud_cover, mfcloud_cover
!                         ---  local variables ---

      integer :: unit, ierr, io
      character(len=64) :: file ! Added by ZTAN 04/08/2012
      
!                         ---  executable code ---
      integer :: iupd
      real :: limpart_low, limpart_upp, limpart_tot, lower_lim, upper_lim, res_fac ! ,gaussian_std

! ADDED BY ZTAN 05/01/2011

! if(module_is_initialized) return

call write_version_number(version, tagname)

! added by cw 12/04/03:
unit = open_namelist_file ()
ierr=1
do while (ierr /= 0)
  read  (unit, nml=newedmf_nml, iostat=io, end=10)
  ierr = check_nml_error (io, 'newedmf_nml')
enddo
10 call close_file (unit)

if ( mpp_pe() == mpp_root_pe() )   write (stdlog(), nml=newedmf_nml)
! end CW addition

! END ZTAN ADDITION 05/01/2011

! Initialize updraft (initialize_updr - only initialization)
allocate(a_u   (idim, jdim, kdim, updraft_number))
allocate(w_u   (idim, jdim, kdim, updraft_number))
allocate(thl_u (idim, jdim, kdim, updraft_number))
allocate(qt_u  (idim, jdim, kdim, updraft_number))
allocate(ql_u  (idim, jdim, kdim, updraft_number))

allocate(buoy_o_u  (idim, jdim, kdim, updraft_number))  ! Added for two-step average... 10/10/2017

allocate(init_a_u  (updraft_number))
allocate(init_scl_u(updraft_number))

init_a_u = 0. 
init_scl_u = 0.

buoy_o_u = 0. !  10/10/2017

call gaussian_mean(0.9, 1.0, res_fac)
gaussian_std = res_fac ! This is a constant = 1.75498366

limpart_tot = 0.0
do iupd = 1, updraft_number
    limpart_tot = limpart_tot + updraft_exponent ** iupd
end do

limpart_low = 0.0; limpart_upp = 0.0
do iupd = 1, updraft_number
    limpart_upp = limpart_upp + updraft_exponent ** iupd
    lower_lim = 1.0 - updraft_fraction * (1.0 - limpart_low/limpart_tot)
    upper_lim = 1.0 - updraft_fraction * (1.0 - limpart_upp/limpart_tot)
    init_a_u (iupd) = upper_lim - lower_lim
    call gaussian_mean(lower_lim, upper_lim, res_fac)
    init_scl_u (iupd) = res_fac/gaussian_std  ! Note by ZTAN 09/25/2017: init_scl_u is normalized by gaussian_std
    limpart_low = limpart_upp
end do

file = 'INPUT/newedmf.res.nc'
if(file_exist(trim(file))) then
   call nullify_domain()
   call read_data(trim(file), 'cloud_cover', cloud_cover, grid_domain)
   call read_data(trim(file), 'mfcloud_cover', mfcloud_cover, grid_domain)
   call read_data(trim(file), 'a_u',   a_u,   grid_domain)
   call read_data(trim(file), 'w_u',   w_u,   grid_domain)
   call read_data(trim(file), 'thl_u', thl_u, grid_domain)
   call read_data(trim(file), 'qt_u',  qt_u,  grid_domain)
   call read_data(trim(file), 'ql_u',  ql_u,  grid_domain)
   updraft_initialized = .true.
else
   cloud_cover = 0.0
   mfcloud_cover = 0.0
   updraft_initialized = .false.
endif


! TO_DO 082016: SAVE UPDRAFT CONDITIONS for restart !!

! Output Parameters below are added by ZTAN: 04/19/2011
! TO_DO 082016: Add more outputs
      id_thl_e    = register_diag_field ( mod_name, 'thl_e', axes(1:3), Time, &
           'thl_e', 'K', missing_value=missing_value)  
      id_qt_e    = register_diag_field ( mod_name, 'qt_e', axes(1:3), Time, &
           'qt_e', 'kg/kg', missing_value=missing_value)
      id_ql_e    = register_diag_field ( mod_name, 'ql_e', axes(1:3), Time, &
           'ql_e', 'kg/kg', missing_value=missing_value)  
      id_tke_e    = register_diag_field ( mod_name, 'tke_e', axes(1:3), Time, &
           'tke_e', 'm^2/s^2', missing_value=missing_value)
      id_a_e    = register_diag_field ( mod_name, 'a_e', axes(1:3), Time, &
           'a_e', '1', missing_value=missing_value)  

      id_thl_u    = register_diag_field ( mod_name, 'thl_u', axes(1:3), Time, &
           'thl_u', 'K', missing_value=missing_value)  
      id_qt_u    = register_diag_field ( mod_name, 'qt_u', axes(1:3), Time, &
           'qt_u', 'kg/kg', missing_value=missing_value)
      id_ql_u    = register_diag_field ( mod_name, 'ql_u', axes(1:3), Time, &
           'ql_u', 'kg/kg', missing_value=missing_value)  
      id_w_u    = register_diag_field ( mod_name, 'w_u', axes(1:3), Time, &
           'w_u', 'm/s', missing_value=missing_value)
      id_a_u    = register_diag_field ( mod_name, 'a_u', axes(1:3), Time, &
           'a_u', '1', missing_value=missing_value)  
      id_m_u    = register_diag_field ( mod_name, 'm_u', axes(1:3), Time, &
           'm_u', 'm/s', missing_value=missing_value)  
      id_b_u    = register_diag_field ( mod_name, 'b_u', axes(1:3), Time, &
           'b_u', 'm/s^2', missing_value=missing_value)  
    
      id_k_m    = register_diag_field ( mod_name, 'k_m', axes(1:3), Time, &
           'k_m', 'm^2/s', missing_value=missing_value) 
      id_k_t    = register_diag_field ( mod_name, 'k_t', axes(1:3), Time, &
           'k_t', 'm^2/s', missing_value=missing_value)  
      id_ls    = register_diag_field ( mod_name, 'ls', axes(1:3), Time, &
           'ls', 'm', missing_value=missing_value)    
      id_ustar    = register_diag_field ( mod_name, 'ustar', axes(1:2), Time, &
           'ustar', 'm/s', missing_value=missing_value)  
      id_wstar    = register_diag_field ( mod_name, 'wstar', axes(1:2), Time, &
           'wstar', 'm/s', missing_value=missing_value) 
      id_zstar    = register_diag_field ( mod_name, 'zstar', axes(1:2), Time, &
           'zstar', 'm', missing_value=missing_value)  
      id_bstar    = register_diag_field ( mod_name, 'bstar', axes(1:2), Time, &
           'bstar', 'm/s^2', missing_value=missing_value) 

      id_flux_thl_ed = register_diag_field ( mod_name, 'flux_thl_ed', axes(half), Time, &
               'flux_thl_ed', 'K*kg/m^2/s', missing_value=missing_value)
      id_flux_qt_ed = register_diag_field ( mod_name, 'flux_qt_ed', axes(half), Time, &
               'flux_qt_ed', '1*kg/m^2/s', missing_value=missing_value)
      id_flux_ql_ed = register_diag_field ( mod_name, 'flux_ql_ed', axes(half), Time, &
               'flux_ql_ed', '1*kg/m^2/s', missing_value=missing_value)
      id_flux_tke_ed = register_diag_field ( mod_name, 'flux_tke_ed', axes(half), Time, &
               'flux_tke_ed', 'm^2/s^2*kg/m^2/s', missing_value=missing_value)
      id_flux_thl_mf = register_diag_field ( mod_name, 'flux_thl_mf', axes(half), Time, &
               'flux_thl_mf', 'K*kg/m^2/s', missing_value=missing_value)
      id_flux_qt_mf = register_diag_field ( mod_name, 'flux_qt_mf', axes(half), Time, &
               'flux_qt_mf', '1*kg/m^2/s', missing_value=missing_value)
      id_flux_ql_mf = register_diag_field ( mod_name, 'flux_ql_mf', axes(half), Time, &
               'flux_ql_mf', '1*kg/m^2/s', missing_value=missing_value)
      id_flux_tke_mf = register_diag_field ( mod_name, 'flux_tke_mf', axes(half), Time, &
               'flux_tke_mf', 'm^2/s^2*kg/m^2/s', missing_value=missing_value)
           
      id_rh    = register_diag_field ( mod_name, 'rh', axes(1:3), Time, &
           'rh', '1', missing_value=missing_value) 
      id_s    = register_diag_field ( mod_name, 's', axes(1:3), Time, &
           's', '1', missing_value=missing_value) 
      id_sigmas    = register_diag_field ( mod_name, 'sigmas', axes(1:3), Time, &
           'sigmas', '1', missing_value=missing_value) 
      id_ccov_ed    = register_diag_field ( mod_name, 'ccov_ed', axes(1:3), Time, &
           'ccov_ed', '1', missing_value=missing_value)
      id_ccov_mf    = register_diag_field ( mod_name, 'ccov_mf', axes(1:3), Time, &
           'ccov_mf', '1', missing_value=missing_value)  
      id_precip_ed    = register_diag_field ( mod_name, 'precip_ed', axes(1:3), Time, &
           'precip_ed', '1', missing_value=missing_value) 
      id_precip_mf    = register_diag_field ( mod_name, 'precip_mf', axes(1:3), Time, &
           'precip_mf', '1', missing_value=missing_value) 
      id_rain2d_ed    = register_diag_field ( mod_name, 'rain2d_ed', axes(1:2), Time, &
           'rain2d_ed', '1', missing_value=missing_value) 
      id_rain2d_mf    = register_diag_field ( mod_name, 'rain2d_mf', axes(1:2), Time, &
           'rain2d_mf', '1', missing_value=missing_value) 
           
      id_entT    = register_diag_field ( mod_name, 'entT', axes(1:3), Time, &
           'entT', 'kg/m^2/s', missing_value=missing_value) 
      id_detT    = register_diag_field ( mod_name, 'detT', axes(1:3), Time, &
           'detT', 'kg/m^2/s', missing_value=missing_value) 
           
      ! ZTAN 11/15/2017: TKE source/sink terms output
      id_tend_tke_buoy  = register_diag_field ( mod_name, 'tend_tke_buoy',  axes(1:3), Time, &
           'tend_tke_buoy',  'm^2/s^3', missing_value=missing_value) 
      id_tend_tke_shear = register_diag_field ( mod_name, 'tend_tke_shear', axes(1:3), Time, &
           'tend_tke_shear', 'm^2/s^3', missing_value=missing_value) 
      id_tend_tke_detr  = register_diag_field ( mod_name, 'tend_tke_detr',  axes(1:3), Time, &
           'tend_tke_detr',  'm^2/s^3', missing_value=missing_value) 
      id_tend_tke_entr  = register_diag_field ( mod_name, 'tend_tke_entr',  axes(1:3), Time, &
           'tend_tke_entr',  'm^2/s^3', missing_value=missing_value) 
      id_tend_tke_diss  = register_diag_field ( mod_name, 'tend_tke_diss',  axes(1:3), Time, &
           'tend_tke_diss',  'm^2/s^3', missing_value=missing_value) 
      id_tend_tke_corr  = register_diag_field ( mod_name, 'tend_tke_corr',  axes(1:3), Time, &
           'tend_tke_corr',  'm^2/s^3', missing_value=missing_value) 
      id_tend_tke_mf    = register_diag_field ( mod_name, 'tend_tke_mf',    axes(1:3), Time, &
           'tend_tke_mf',    'm^2/s^3', missing_value=missing_value) 
      id_tend_tke_ed    = register_diag_field ( mod_name, 'tend_tke_ed',    axes(1:3), Time, &
           'tend_tke_ed',    'm^2/s^3', missing_value=missing_value) 
      id_tend_tke_tot   = register_diag_field ( mod_name, 'tend_tke_tot',   axes(1:3), Time, &
           'tend_tke_tot',   'm^2/s^3', missing_value=missing_value) 
           
      ! Latitude slice
      id_t_g_sl = register_diag_field ( mod_name, 't_g_sl', axes(2:3), Time, &
           't_g_sl', 'K', missing_value=missing_value) 
      id_q_g_sl= register_diag_field ( mod_name, 'q_g_sl', axes(2:3), Time, &
           'q_g_sl', 'kg/kg', missing_value=missing_value) 
      id_ql_g_sl= register_diag_field ( mod_name, 'ql_g_sl', axes(2:3), Time, &
           'ql_g_sl', 'kg/kg', missing_value=missing_value) 
      id_u_g_sl= register_diag_field ( mod_name, 'u_g_sl', axes(2:3), Time, &
           'u_g_sl', 'm/s', missing_value=missing_value) 
      id_v_g_sl= register_diag_field ( mod_name, 'v_g_sl', axes(2:3), Time, &
           'v_g_sl', 'm/s', missing_value=missing_value) 
      id_thl_g_sl= register_diag_field ( mod_name, 'thl_g_sl', axes(2:3), Time, &
           'thl_g_sl', 'K', missing_value=missing_value) 
      id_qt_g_sl= register_diag_field ( mod_name, 'qt_g_sl', axes(2:3), Time, &
           'qt_g_sl', 'kg/kg', missing_value=missing_value) 
      id_thv_g_sl= register_diag_field ( mod_name, 'thv_g_sl', axes(2:3), Time, &
           'thv_g_sl', 'K', missing_value=missing_value) 
      id_rho_g_sl= register_diag_field ( mod_name, 'rho_g_sl', axes(2:3), Time, &
           'rho_g_sl', 'kg/m^3', missing_value=missing_value) 
      id_tke_g_sl= register_diag_field ( mod_name, 'tke_g_sl', axes(2:3), Time, &
           'tke_g_sl', 'm^2/s^2', missing_value=missing_value) 
      
      id_zfull_sl= register_diag_field ( mod_name, 'zfull_sl', axes(2:3), Time, &
           'zfull_sl', 'm', missing_value=missing_value) 
      id_zhalf_sl= register_diag_field ( mod_name, 'zhalf_sl', axes(2:3), Time, &
           'zhalf_sl', 'm', missing_value=missing_value) 
      id_pfull_sl= register_diag_field ( mod_name, 'pfull_sl', axes(2:3), Time, &
           'pfull_sl', 'Pa', missing_value=missing_value) 
      id_phalf_sl= register_diag_field ( mod_name, 'phalf_sl', axes(2:3), Time, &
           'phalf_sl', 'Pa', missing_value=missing_value) 
      id_rh_sl= register_diag_field ( mod_name, 'rh_sl', axes(2:3), Time, &
           'rh_sl', '1', missing_value=missing_value) 
      id_ccov_ed_sl= register_diag_field ( mod_name, 'ccov_ed_sl', axes(2:3), Time, &
           'ccov_ed_sl', '1', missing_value=missing_value) 
      id_ccov_mf_sl= register_diag_field ( mod_name, 'ccov_mf_sl', axes(2:3), Time, &
           'ccov_mf_sl', '1', missing_value=missing_value) 
       
      id_rain2d_ed_sl= register_diag_field ( mod_name, 'rain2d_ed_sl', axes(sfc_slice), Time, &
           'rain2d_ed_sl', 'N/A', missing_value=missing_value) 
      id_rain2d_mf_sl= register_diag_field ( mod_name, 'rain2d_mf_sl', axes(sfc_slice), Time, &
           'rain2d_mf_sl', 'N/A', missing_value=missing_value) 
      id_flux_t_sl= register_diag_field ( mod_name, 'flux_t_sl', axes(sfc_slice), Time, &
           'flux_t_sl', 'K*m/s', missing_value=missing_value) 
      id_flux_q_sl= register_diag_field ( mod_name, 'flux_q_sl', axes(sfc_slice), Time, &
           'flux_q_sl', 'kg/kg*m/s', missing_value=missing_value) 
      
      id_thl_e_sl    = register_diag_field ( mod_name, 'thl_e_sl', axes(2:3), Time, &
           'thl_e_sl', 'K', missing_value=missing_value)  
      id_qt_e_sl    = register_diag_field ( mod_name, 'qt_e_sl', axes(2:3), Time, &
           'qt_e_sl', 'kg/kg', missing_value=missing_value)
      id_ql_e_sl    = register_diag_field ( mod_name, 'ql_e_sl', axes(2:3), Time, &
           'ql_e_sl', 'kg/kg', missing_value=missing_value)  
      id_tke_e_sl    = register_diag_field ( mod_name, 'tke_e_sl', axes(2:3), Time, &
           'tke_e_sl', 'm^2/s^2', missing_value=missing_value)
      id_a_e_sl    = register_diag_field ( mod_name, 'a_e_sl', axes(2:3), Time, &
           'a_e_sl', '1', missing_value=missing_value)  

      id_thl_u_sl    = register_diag_field ( mod_name, 'thl_u_sl', axes(2:3), Time, &
           'thl_u_sl', 'K', missing_value=missing_value)  
      id_qt_u_sl    = register_diag_field ( mod_name, 'qt_u_sl', axes(2:3), Time, &
           'qt_u_sl', 'kg/kg', missing_value=missing_value)
      id_ql_u_sl    = register_diag_field ( mod_name, 'ql_u_sl', axes(2:3), Time, &
           'ql_u_sl', 'kg/kg', missing_value=missing_value)  
      id_w_u_sl    = register_diag_field ( mod_name, 'w_u_sl', axes(2:3), Time, &
           'w_u_sl', 'm/s', missing_value=missing_value)
      id_a_u_sl    = register_diag_field ( mod_name, 'a_u_sl', axes(2:3), Time, &
           'a_u_sl', '1', missing_value=missing_value)  
      id_m_u_sl    = register_diag_field ( mod_name, 'm_u_sl', axes(2:3), Time, &
           'm_u_sl', 'm/s', missing_value=missing_value)  
      id_b_u_sl    = register_diag_field ( mod_name, 'b_u_sl', axes(2:3), Time, &
           'b_u_sl', 'm/s^2', missing_value=missing_value) 
           
           
      id_entr_u_sl    = register_diag_field ( mod_name, 'entr_u_sl', axes(2:3), Time, &
           'entr_u_sl', 'm/s^2', missing_value=missing_value) 
      id_detr_u_sl    = register_diag_field ( mod_name, 'detr_u_sl', axes(2:3), Time, &
           'detr_u_sl', 'm/s^2', missing_value=missing_value) 
      id_entT_sl    = register_diag_field ( mod_name, 'entT_sl', axes(2:3), Time, &
           'entT_sl', 'kg/m^2/s', missing_value=missing_value) 
      id_detT_sl    = register_diag_field ( mod_name, 'detT_sl', axes(2:3), Time, &
           'detT_sl', 'kg/m^2/s', missing_value=missing_value) 
    
      id_k_m_sl    = register_diag_field ( mod_name, 'k_m_sl', axes(2:3), Time, &
           'k_m_sl', 'm^2/s', missing_value=missing_value) 
      id_k_t_sl    = register_diag_field ( mod_name, 'k_t_sl', axes(2:3), Time, &
           'k_t_sl', 'm^2/s', missing_value=missing_value)  
      id_ls_sl    = register_diag_field ( mod_name, 'ls_sl', axes(2:3), Time, &
           'ls_sl', 'm', missing_value=missing_value)    
      id_ustar_sl    = register_diag_field ( mod_name, 'ustar_sl', axes(sfc_slice), Time, &
           'ustar_sl', 'm/s', missing_value=missing_value)  
      id_wstar_sl    = register_diag_field ( mod_name, 'wstar_sl', axes(sfc_slice), Time, &
           'wstar_sl', 'm/s', missing_value=missing_value) 
      id_zstar_sl    = register_diag_field ( mod_name, 'zstar_sl', axes(sfc_slice), Time, &
           'zstar_sl', 'm', missing_value=missing_value)  
      id_bstar_sl    = register_diag_field ( mod_name, 'bstar_sl', axes(sfc_slice), Time, &
           'bstar_sl', 'm/s^2', missing_value=missing_value) 

      id_flux_thl_ed_sl = register_diag_field ( mod_name, 'flux_thl_ed_sl', axes(half_slice), Time, &
               'flux_thl_ed_sl', 'K*kg/m^2/s', missing_value=missing_value)
      id_flux_qt_ed_sl = register_diag_field ( mod_name, 'flux_qt_ed_sl', axes(half_slice), Time, &
               'flux_qt_ed_sl', '1*kg/m^2/s', missing_value=missing_value)
      id_flux_ql_ed_sl = register_diag_field ( mod_name, 'flux_ql_ed_sl', axes(half_slice), Time, &
               'flux_ql_ed_sl', '1*kg/m^2/s', missing_value=missing_value)
      id_flux_thl_mf_sl = register_diag_field ( mod_name, 'flux_thl_mf_sl', axes(half_slice), Time, &
               'flux_thl_mf_sl', 'K*kg/m^2/s', missing_value=missing_value)
      id_flux_qt_mf_sl = register_diag_field ( mod_name, 'flux_qt_mf_sl', axes(half_slice), Time, &
               'flux_qt_mf_sl', '1*kg/m^2/s', missing_value=missing_value)
      id_flux_ql_mf_sl = register_diag_field ( mod_name, 'flux_ql_mf_sl', axes(half_slice), Time, &
               'flux_ql_mf_sl', '1*kg/m^2/s', missing_value=missing_value)

      ! Debug only for tendency output...
      id_thldel_sl    = register_diag_field ( mod_name, 'thldel_sl', axes(2:3), Time, &
           'thldel_sl', 'K/s', missing_value=missing_value)  
      id_qtdel_sl    = register_diag_field ( mod_name, 'qtdel_sl', axes(2:3), Time, &
           'qtdel_sl', 'kg/kg/s', missing_value=missing_value)  
      id_tdel_sl    = register_diag_field ( mod_name, 'tdel_sl', axes(2:3), Time, &
           'tdel_sl', 'K/s', missing_value=missing_value)  
      id_qdel_sl    = register_diag_field ( mod_name, 'qdel_sl', axes(2:3), Time, &
           'qdel_sl', 'kg/kg/s', missing_value=missing_value)  
      id_liqdel_sl    = register_diag_field ( mod_name, 'liqdel_sl', axes(2:3), Time, &
           'liqdel_sl', 'kg/kg/s', missing_value=missing_value) 
      id_udel_sl    = register_diag_field ( mod_name, 'udel_sl', axes(2:3), Time, &
           'udel_sl', 'm/s^2', missing_value=missing_value)   
      id_vdel_sl    = register_diag_field ( mod_name, 'vdel_sl', axes(2:3), Time, &
           'vdel_sl', 'm/s^2', missing_value=missing_value)   
      id_tkedel_sl    = register_diag_field ( mod_name, 'tkedel_sl', axes(2:3), Time, &
           'tkedel_sl', 'm^2/s^3', missing_value=missing_value)  
          
      id_thldel_diag_sl    = register_diag_field ( mod_name, 'thldel_diag_sl', axes(2:3), Time, &
           'thldel_diag_sl', 'K/s', missing_value=missing_value)  
      id_qtdel_diag_sl    = register_diag_field ( mod_name, 'qtdel_diag_sl', axes(2:3), Time, &
           'qtdel_diag_sl', 'kg/kg/s', missing_value=missing_value)  
      id_tdel_diag_sl    = register_diag_field ( mod_name, 'tdel_diag_sl', axes(2:3), Time, &
           'tdel_diag_sl', 'K/s', missing_value=missing_value)  
      id_qdel_diag_sl    = register_diag_field ( mod_name, 'qdel_diag_sl', axes(2:3), Time, &
           'qdel_diag_sl', 'kg/kg/s', missing_value=missing_value)  
      id_liqdel_diag_sl    = register_diag_field ( mod_name, 'liqdel_diag_sl', axes(2:3), Time, &
           'liqdel_diag_sl', 'kg/kg/s', missing_value=missing_value) 
      id_udel_diag_sl    = register_diag_field ( mod_name, 'udel_diag_sl', axes(2:3), Time, &
           'udel_diag_sl', 'm/s^2', missing_value=missing_value)   
      id_vdel_diag_sl    = register_diag_field ( mod_name, 'vdel_diag_sl', axes(2:3), Time, &
           'vdel_diag_sl', 'm/s^2', missing_value=missing_value)   
      id_tkedel_diag_sl    = register_diag_field ( mod_name, 'tkedel_diag_sl', axes(2:3), Time, &
           'tkedel_diag_sl', 'm^2/s^3', missing_value=missing_value)  

      ! ZTAN 11/15/2017: TKE source/sink terms output
      id_tend_tke_buoy_sl  = register_diag_field ( mod_name, 'tend_tke_buoy_sl',  axes(2:3), Time, &
           'tend_tke_buoy_sl',  'm^2/s^3', missing_value=missing_value) 
      id_tend_tke_shear_sl = register_diag_field ( mod_name, 'tend_tke_shear_sl', axes(2:3), Time, &
           'tend_tke_shear_sl', 'm^2/s^3', missing_value=missing_value) 
      id_tend_tke_detr_sl  = register_diag_field ( mod_name, 'tend_tke_detr_sl',  axes(2:3), Time, &
           'tend_tke_detr_sl',  'm^2/s^3', missing_value=missing_value) 
      id_tend_tke_entr_sl  = register_diag_field ( mod_name, 'tend_tke_entr_sl',  axes(2:3), Time, &
           'tend_tke_entr_sl',  'm^2/s^3', missing_value=missing_value) 
      id_tend_tke_diss_sl  = register_diag_field ( mod_name, 'tend_tke_diss_sl',  axes(2:3), Time, &
           'tend_tke_diss_sl',  'm^2/s^3', missing_value=missing_value) 
      id_tend_tke_corr_sl  = register_diag_field ( mod_name, 'tend_tke_corr_sl',  axes(2:3), Time, &
           'tend_tke_corr_sl',  'm^2/s^3', missing_value=missing_value) 
      id_tend_tke_mf_sl    = register_diag_field ( mod_name, 'tend_tke_mf_sl',    axes(2:3), Time, &
           'tend_tke_mf_sl',    'm^2/s^3', missing_value=missing_value) 
      id_tend_tke_ed_sl    = register_diag_field ( mod_name, 'tend_tke_ed_sl',    axes(2:3), Time, &
           'tend_tke_ed_sl',    'm^2/s^3', missing_value=missing_value) 
      id_tend_tke_tot_sl   = register_diag_field ( mod_name, 'tend_tke_tot_sl',   axes(2:3), Time, &
           'tend_tke_tot_sl',   'm^2/s^3', missing_value=missing_value) 
               
end subroutine newedmf_init


subroutine newedmf (is, ie, js, je, Time,                     & 
                    delta_t, dt_real, tin, qin,               &
                    liqin, tkein, t_surf, uin, vin,           &
                    ustar, bstar, k_m, k_t,                   &
                    pfull, phalf, zfull, zhalf,               &
                    flux_t,flux_q, flux_u,flux_v, flux_r,     &
                    dtaudv_atm, dhdt_surf, dedt_surf,         &
                    dedq_surf, drdt_surf, dhdt_atm, dedq_atm, &
                    net_surf_sw_down, surf_lw_down,           &
                    Tri_surf, diss_heat,                      &
                    tdel, qdel, liqdel, tkedel, udel, vdel,   &
                    ccov_ed, ccov_mf, rain2d_ed, rain2d_mf,   &
                    edmf_update_ml, istep, deg_lat, deg_lon)


    integer, intent(in)                      :: is, ie, js, je
    type(time_type), intent(in)              :: Time
    real   , intent(in)                      :: delta_t, dt_real
    ! Note by ZTAN 09/10/2017: 
    ! delta_t is the GCM timestep used for GS tendencies; 
    ! dt_real is the MF timestep used for updraft variables, and for SST.
    ! After each call to the EDMF scheme: 
    ! (1) The updraft variables a_u, w_u, thl_u, qt_u, ql_u are only advanced for dt_real;
    ! (2) The GS tendencies tdel, qdel, liqdel, tkedel, udel, vdel are advanced for delta_t. 
    ! (3) The updraft variables after delta_t are REQUIRED for calculating the GS 
    !     tendencies, but they are NOT USED by the next call to EDMF and are not stored.
    
    real   , intent(in),    dimension(:,:,:) :: tin, qin, liqin, tkein, uin, vin, pfull, phalf, zfull, zhalf
    real   , intent(out),   dimension(:,:)   :: rain2d_ed, rain2d_mf
    real   , intent(out),   dimension(:,:,:) :: ccov_ed, ccov_mf
    real   , intent(out),   dimension(:,:,:) :: tdel, qdel, liqdel, tkedel, udel, vdel
    real   , intent(in),    dimension(:,:)   :: ustar, bstar, flux_r, net_surf_sw_down, surf_lw_down
    real   , intent(inout), dimension(:,:)   :: t_surf, flux_t, flux_q, flux_u, flux_v, &
                                                dtaudv_atm, dhdt_surf, dedt_surf, dedq_surf, drdt_surf, & 
                                                dhdt_atm, dedq_atm 
    real   , intent(inout), dimension(:,:,:) :: k_m, k_t, diss_heat
    logical, intent(in)                      :: edmf_update_ml
    integer, intent(in)                      :: istep
    real   , intent(in),    dimension(:)     :: deg_lat, deg_lon   ! For debug... ZTAN 09/2017
    type(surf_diff_type), intent(inout)      :: Tri_surf
    
    real  , dimension(size(tin,1),size(tin,2),size(tin,3)+1) ::  zhalf_tmp
    ! Grid scale variables (subscript 'g' does NOT mean geostrophic)
    real  , dimension(size(tin,1),size(tin,2),size(tin,3)) ::  t_g, q_g, ql_g, u_g, v_g, w_g, &
                                                               thl_g, qt_g, thv_g, rho_g, tke_g  
    real  , dimension(size(tin,2),size(tin,3))             ::  sl_tmpout
    ! Diagnosis: final GS tendencies                                                           
    real  , dimension(size(tin,1),size(tin,2),size(tin,3)) ::  & 
            tdel_diag, qdel_diag, liqdel_diag, tkedel_diag, udel_diag, vdel_diag, &
            t_g_diag, q_g_diag, ql_g_diag, tke_g_diag, u_g_diag, v_g_diag, thl_g_diag, qt_g_diag, w_g_diag, &
            thldel_diag, qtdel_diag, thldel, qtdel   ! Added for debug
    ! TO_DO 082016: Add input for omega??
    ! real  , dimension(size(tin,1),size(tin,2),size(tin,3)) ::  thl, qt, ql, thv, rho, ls, dls, rho_half
    real  , dimension(size(tin,1),size(tin,2),size(tin,3)) ::  ls, dls, rho_half
    real  , dimension(size(tin,1),size(tin,2)) :: & 
        zstar, thetav, wthv_surf, wthl_surf, wq_surf, wstar, oblength, tke_surf, &
        flux_t_in, flux_q_in ! For debug
    
    real  , dimension(size(tin,1),size(tin,2),size(tin,3)) ::  &   ! Environment values
        thl_e, qt_e, ql_e, w_e, a_e, u_e, v_e, tke_e  
    
    real  , dimension(size(tin,1),size(tin,2),size(tin,3)) ::  &   ! Bulk updraft values after dt_real
        thl_ut1, qt_ut1, ql_ut1, w_ut1, a_ut1, m_ut1, buoy_ut1, chic_ut1, &
        entL_ut1, entB_ut1, detL_ut1, detB_ut1, tend_thl_ut1, tend_qt_ut1, tend_ql_ut1
        
    real  , dimension(size(tin,1),size(tin,2),size(tin,3)) ::  &   ! Bulk updraft values after delta_t
        thl_ut2, qt_ut2, ql_ut2, w_ut2, a_ut2, m_ut2, buoy_ut2, chic_ut2, &
        entL_ut2, entB_ut2, detL_ut2, detB_ut2, tend_thl_ut2, tend_qt_ut2, tend_ql_ut2
        
    real  , dimension(size(tin,1),size(tin,2),size(tin,3)) ::  &   ! Bulk updraft values -- Total
        thl_ut, qt_ut, ql_ut, w_ut, a_ut, m_ut, buoy_ut, chic_ut, &
        entL_ut, entB_ut, detL_ut, detB_ut, tend_thl_ut, tend_qt_ut, tend_ql_ut
                
    real  , dimension(size(tin,1),size(tin,2),size(tin,3), updraft_number) ::  &   
        thl_u1,  qt_u1,  ql_u1,  w_u1,  a_u1, &  ! Updraft values after dt_real
        thl_u2,  qt_u2,  ql_u2,  w_u2,  a_u2, &  ! Updraft values after delta_t
        thl_uo,  qt_uo,  ql_uo,  w_uo,  a_uo     ! Updraft values from the old step -- used in do_phiu_corr
        
    real  , dimension(size(tin,1),size(tin,2),size(tin,3), updraft_number) ::  & 
        entL_u1, entB_u1, detL_u1, detB_u1, &   ! Individual updraft values after dt_real
        entL_u2, entB_u2, detL_u2, detB_u2, entr_u2, detr_u2, buoy_u2, &   ! Individual updraft values after delta_t
        entL_u, entB_u, detL_u, detB_u, entr_u, detr_u, buoy_u             ! Individual updraft values seen by diagnosis
        ! buoy_u2 and buoy_u added 10/10/2017
                
    real  , dimension(size(tin,1),size(tin,2),size(tin,3)) ::  &   ! MF tendencies for dt_real
        tend_thl_mf1, tend_qt_mf1, tend_ql_mf1
    real  , dimension(size(tin,1),size(tin,2),size(tin,3)+1) ::  & ! MF fluxes for dt_real
        flux_thl_mf1, flux_qt_mf1, flux_ql_mf1, net_mf_tot1
        
    real  , dimension(size(tin,1),size(tin,2),size(tin,3)) ::  &   ! MF tendencies for delta_t - dt_real
        tend_thl_mf2, tend_qt_mf2, tend_ql_mf2
    real  , dimension(size(tin,1),size(tin,2),size(tin,3)+1) ::  & ! MF fluxes for delta_t - dt_real
        flux_thl_mf2, flux_qt_mf2, flux_ql_mf2, net_mf_tot2
        
    real  , dimension(size(tin,1),size(tin,2),size(tin,3)) ::  &   ! MF tendencies -- Total
        tend_thl_mf, tend_qt_mf, tend_ql_mf
    real  , dimension(size(tin,1),size(tin,2),size(tin,3)+1) ::  & ! MF fluxes -- Total
        flux_thl_mf, flux_qt_mf, flux_ql_mf, net_mf_tot
    
    
    real  , dimension(size(tin,1),size(tin,2),size(tin,3)) ::  &   ! Grid-mean values after MF update - dt_real (st11)
        u_st11, v_st11, w_st11, tke_st11, thl_st11, qt_st11, ql_st11, t_st11, q_st11, thv_st11, rho_st11  
    real  , dimension(size(tin,1),size(tin,2),size(tin,3)) ::  &   ! Environment values after MF update - dt_real (en11)
        thl_en11, qt_en11, ql_en11, w_en11, a_en11, u_en11, v_en11, tke_en11
        
    real  , dimension(size(tin,1),size(tin,2),size(tin,3)) ::  &   ! Grid-mean values after MF update - delta_t (st12)
        u_st12, v_st12, w_st12, tke_st12, thl_st12, qt_st12, ql_st12, t_st12, q_st12, thv_st12, rho_st12  
    real  , dimension(size(tin,1),size(tin,2),size(tin,3)) ::  &   ! Environment values after MF update - delta_t (en12)
        thl_en12, qt_en12, ql_en12, w_en12, a_en12, u_en12, v_en12, tke_en12    
        
    real  , dimension(size(tin,1),size(tin,2),size(tin,3)) ::  &   ! Grid-mean values after MF update (st1)
        u_st1, v_st1, w_st1, tke_st1, thl_st1, qt_st1, ql_st1, t_st1, q_st1, thv_st1, rho_st1  
    real  , dimension(size(tin,1),size(tin,2),size(tin,3)) ::  &   ! Environment values after MF update (en1)
        thl_en1, qt_en1, ql_en1, w_en1, a_en1, u_en1, v_en1, tke_en1
    
    real  , dimension(size(tin,1),size(tin,2),size(tin,3)) ::  tend_qt_clip1, tend_ql_clip1, tend_ql_clip2 
                ! GS tendency of qt and ql due to clipping (#1: initial clipping; #2: ql-limiter so that ql_env >= 0)
    real  , dimension(size(tin,1),size(tin,2),size(tin,3)) ::  &   ! ED tendencies
        tend_thl_ed, tend_qt_ed, tend_ql_ed, tend_u_ed, tend_v_ed
    real  , dimension(size(tin,1),size(tin,2),size(tin,3)+1) ::  & ! ED fluxes
        flux_thl_ed, flux_qt_ed, flux_ql_ed, flux_u_ed, flux_v_ed  
    real  , dimension(size(tin,1),size(tin,2),size(tin,3)) ::  &   ! Environment values after ED update (en2)
        thl_en2, qt_en2, ql_en2, w_en2, a_en2, u_en2, v_en2, tke_en2
        
    real  , dimension(size(tin,1),size(tin,2),size(tin,3)) ::  &   ! Environment values after Cloud/Precip (en3)
        thl_en3, qt_en3, ql_en3, w_en3, a_en3, u_en3, v_en3, tke_en3
        
    real  , dimension(size(tin,1),size(tin,2),size(tin,3)) :: tke_en4   ! Environment TKE after update (en4)
        
    real  , dimension(size(tin,1),size(tin,2),size(tin,3)) ::  &   ! Variables and Tendencies for Cloud/Precip scheme
        rh, s, sigmas, ccov, wthl, wqt, tend_thl_cld, tend_qt_cld, tend_ql_cld
        
    real  , dimension(size(tin,1),size(tin,2),size(tin,3)) ::  &   ! Total tendencies
        tend_thl_tot, tend_qt_tot, tend_ql_tot, tend_t_tot, tend_q_tot, tend_u_tot, tend_v_tot
        
    real  , dimension(size(tin,1),size(tin,2),size(tin,3)) ::  &   ! 3D cloud and Rain
        precip_ed, precip_mf
        
    real  , dimension(size(tin,1),size(tin,2),size(tin,3)+1) :: flux_tke_mf, flux_tke_ed    ! TKE fluxes
    real  , dimension(size(tin,1),size(tin,2),size(tin,3)+1) :: flux_tke_mf1, flux_tke_mf2  ! TKE fluxes
    
    real  , dimension(size(tin,1),size(tin,2),size(tin,3)) :: &
        tend_tke_buoy, tend_tke_buoyD, tend_tke_buoyW, tend_tke_shear, tend_tke_detr, tend_tke_entr, &
        tend_tke_diss, tend_tke_corr, tend_tke_source, tend_tke_mf, tend_tke_ed, tend_tke_tot   ! TKE Tendencies

    ! real  , dimension(size(tin,1),size(tin,2)) :: rain2d_ed, rain2d_mf ! 2D Rain
    

    integer :: nlev, iupd, i, j
    real    :: delt_2nd, delt_fac1, delt_fac2

    logical :: used

    nlev = size(tin,3)
    
    if (istep == 1) then      
        ! i = mpp_pe(); j = mpp_root_pe()
        ! write(*,*) i, j, is, ie, js, je, deg_lat(1), deg_lat(je-js+1), deg_lon(1), deg_lon(ie-is+1)
        idb = 0; jdb = 0
        do i = 1, ie-is+1
          do j = 1, je-js+1
             ! if ((deg_lat(j) .le. 3.5) .and. (deg_lat(j) .ge. 0.5) .and. (deg_lon(i) .le. 2.5)) then
             ! if ((deg_lat(j) .le. -22.5) .and. (deg_lat(j) .ge. -24.5) .and. (deg_lon(i) .le. 2.5)) then
             ! if ((deg_lat(j) .le. -25.5) .and. (deg_lat(j) .ge. -27.5) .and. (deg_lon(i) .le. 2.5)) then
             if ((deg_lat(j) .le. 6.0) .and. (deg_lat(j) .ge. 3.0) .and. (deg_lon(i) .le. -2.5)) then 
                 write(*,"(A, 3F9.2)") 'DEBUG point: Lat/Lon/Pfull: ', deg_lat(j), deg_lon(i)
                 idb = i; jdb = j; debug_point = .true. 
             endif
          enddo
        enddo
    endif
    
    ! Make local copies of GS variables
    t_g   (:,:,:)  = tin  (:,:,:)
    q_g   (:,:,:)  = max(qin  (:,:,:), 0.0)   ! Clipping.. Added ZTAN 09/13/2017 
    ql_g  (:,:,:)  = max(liqin(:,:,:), 0.0)   ! Clipping.. Added ZTAN 09/13/2017
    u_g   (:,:,:)  = uin  (:,:,:)
    v_g   (:,:,:)  = vin  (:,:,:)
    w_g   (:,:,:)  = 0.0          ! TO_DO 082016: Add input for omega??
    tke_g (:,:,:)  = tkein(:,:,:)
    
    flux_t_in(:,:) = flux_t(:,:)
    flux_q_in(:,:) = flux_q(:,:)
    
    ! CLIPPING #1: keeps T unchanged while clipping q and ql to above zero 
    !              -> 'clipping' tendencies for q and ql (but NOT for T).
    tend_qt_clip1(:,:,:) = (q_g (:,:,:) + ql_g (:,:,:) - qin(:,:,:) - liqin(:,:,:))/delta_t
    tend_ql_clip1(:,:,:) = (ql_g (:,:,:) - liqin(:,:,:))/delta_t
    
    call cloud_mixing (t_g, q_g, ql_g, pfull, thl_g, qt_g)
    call compute_thv (thl_g, qt_g, ql_g, pfull, thv_g, rho_g) 
    ! NOTE by ZTAN 09/25/2017: the compute_thv subroutine computes grid-mean theta_v directly 
    !     from grid-mean theta_l, q_t, etc. This is inaccurate, because theta_v is not linear, 
    !     so theta_v_g == sum(a_i theta_v(thl_i, qt_i)) .NEQ. theta_v(sum(a_i thl_i), sum(a_i qt_i))

    if (updraft_initialized == .false.) then
       if (mpp_pe() == mpp_root_pe()) then
          write(*,*) 'Initialize updraft from env value'
       endif
       a_u (:,:,:,:) = 0. 
       w_u (:,:,:,:) = 0.
       if (do_trans_dphi) then  ! updraft anomalies from env are kept instead of updraft values.
                                ! thus, they are all initialized as zero (no anomaly from environment).
          do iupd = 1, updraft_number
             thl_u(:,:,:,iupd) = 0.
             qt_u (:,:,:,iupd) = 0.
             ql_u (:,:,:,iupd) = 0.
             a_u (:,:,nlev,iupd) = init_a_u(iupd)
          end do
       else                     ! initialize the updraft values directly with grid-mean values
          do iupd = 1, updraft_number
             thl_u(:,:,:,iupd) = thl_g (:,:,:)
             qt_u (:,:,:,iupd) = qt_g  (:,:,:)
             ql_u (:,:,:,iupd) = ql_g  (:,:,:)
             a_u (:,:,nlev,iupd) = init_a_u(iupd)
          end do
       end if
       updraft_initialized = .true.
    else
       if (istep == 1) write(*,*) 'Updraft already initialized'
    end if

    ! Modify the first (top) zhalf to be consistent... (the native value from the GCM is zero...)
    zhalf_tmp (:,:,2:nlev+1) = zhalf(:,:,2:nlev+1)
    zhalf_tmp (:,:,1) = zhalf(:,:,2) - (phalf(:,:,1) - phalf(:,:,2))/Grav/rho_g(:,:,1)

    ! initialize_updr - setting up boundary normalized anomalies of updrafts.
    !                   this is now already done in newedmf_init subroutine
    
    if (do_trans_dphi) then  ! transform back to actual updraft values from the anomalies 
                             ! by calling the subroutine 'trans_dphi_to_phiu'
                             
       ! Note 09/14/2017: ql_g may need to be adjusted up to guarantee a positive ql_e, 
       !                  and thus t_g and q_g may need to be updated.
       
       ! CLIPPING #2: keeps thl unchanged, changes T, q and ql 
       !              -> 'clipping' tendencies for q, ql AND T.
       tend_ql_clip2 (:,:,:) = ql_g (:,:,:) ! temporarily store the old value of ql_g
       call trans_dphi_to_phiu (thl_u, qt_u, ql_u, a_u, thl_g, qt_g, ql_g)
       tend_ql_clip2(:,:,:) = (ql_g (:,:,:) - tend_ql_clip2(:,:,:))/delta_t
       call cloud_decompose (thl_g, qt_g, ql_g, pfull, t_g, q_g)  ! update t_g and q_g
    else
       tend_ql_clip2 (:,:,:) = 0.0
    end if
    ! Compute clipping tendencies

    ! compute_zstar_wstar_obl: compute surface scales and obukhov length
    call compute_zstar_wstar_obl(Time, t_g, q_g, ql_g, thl_g, qt_g, thv_g, u_g, v_g,  &  ! grid_mean
                                 ustar, bstar, flux_t, flux_q,         &  ! surf_flux
                                 pfull, zfull,  & 
                                 zstar, wthv_surf, wthl_surf, wq_surf, wstar, oblength)                        
                                 
    ! --- START MF COMPUTATION --- !
    ! updraft/environment decomposition and surface condition is moved into compute_updr 
    ! -- Should make this into a subroutine later.
    
    if (do_compute_updr_lbl) then ! use the level-by-level updraft solver
       if (delta_t .gt. dt_real) then 
          ! The GCM uses leapfrog timestepping, so it is necessary to compute updraft 
          ! twice with dt_real and 2*dt_real (= delta_t), respectively.
          ! * The dt_real values are used only for updraft memory, 
          ! * The 2*dt_real (= delta_t) values are used for the main GCM timestepping.
          ! * In case of simple forward stepping, delta_t (= dt_real) values are used for both purposes.
          
          ! if ((mpp_pe() == mpp_root_pe()) .and. (istep == 1)) then
          !    write(*,*) 'LEAPFROG -- Calculate updraft twice, dt_real loop is only for memory'
          ! endif
       
          call compute_updr_lbl ( thl_u,  qt_u,  ql_u,  w_u,  a_u, &  ! updr (in) 
                             thl_u2, qt_u2, ql_u2, w_u2, a_u2, & ! updr (out)  <- updraft value at t+dt_real (for memory)
                             thl_e,  qt_e,  ql_e,  w_e,  a_e,  u_e,  v_e,  tke_e, & ! envr, now as output
                             thl_g,  qt_g,  ql_g,  w_g,        u_g,  v_g,  tke_g, & ! grid-mean as input
                             rho_g, zstar, ustar, wthl_surf, wq_surf, wstar, oblength, &  ! oblength added 09/25/2017
                             pfull, phalf, zfull, zhalf_tmp, dt_real,     & 
                             thl_ut2, qt_ut2, ql_ut2, w_ut2, a_ut2, & ! updr_out_tot (bulk updraft values)
                             m_ut2, buoy_ut2, chic_ut2, entL_ut2, entB_ut2, detL_ut2, detB_ut2, & ! updr_out_tot (bulk updraft values)
                             entL_u2, entB_u2, detL_u2, detB_u2, entr_u2, detr_u2, buoy_u2, &  ! individual updraft values
                             tend_thl_ut2, tend_qt_ut2, tend_ql_ut2, istep)
       endif
    
       call compute_updr_lbl ( thl_u,  qt_u,  ql_u,  w_u,  a_u, &  ! updr (in) 
                          thl_u1, qt_u1, ql_u1, w_u1, a_u1, & ! updr (out)  <- updraft value at t+delta_t (for GCM timestepping)
                          thl_e,  qt_e,  ql_e,  w_e,  a_e,  u_e,  v_e,  tke_e, & ! envr, now as output
                          thl_g,  qt_g,  ql_g,  w_g,        u_g,  v_g,  tke_g, & ! grid-mean as input
                          rho_g, zstar, ustar, wthl_surf, wq_surf, wstar, oblength, &  ! oblength added 09/25/2017
                          pfull, phalf, zfull, zhalf_tmp, delta_t,     & 
                          thl_ut, qt_ut, ql_ut, w_ut, a_ut, & ! updr_out_tot (bulk updraft values)
                          m_ut, buoy_ut, chic_ut, entL_ut, entB_ut, detL_ut, detB_ut, & ! updr_out_tot (bulk updraft values)
                          entL_u, entB_u, detL_u, detB_u, entr_u, detr_u, buoy_u, &  ! individual updraft values
                          tend_thl_ut, tend_qt_ut, tend_ql_ut, istep) 
                          
    else ! Use the old subroutine that is not level-by-level...
       if (delta_t .gt. dt_real) then 
          ! The GCM uses leapfrog timestepping, so it is necessary to compute updraft 
          ! twice with dt_real and 2*dt_real (= delta_t), respectively.
          ! * The dt_real values are used only for updraft memory, 
          ! * The 2*dt_real (= delta_t) values are used for the main GCM timestepping.
          ! * In case of simple forward stepping, delta_t (= dt_real) values are used for both purposes.
                   
          ! if ((mpp_pe() == mpp_root_pe()) .and. (istep == 1)) then
          !    write(*,*) 'LEAPFROG -- Calculate updraft twice, dt_real loop is only for memory'
          ! endif
       
          call compute_updr( thl_u,  qt_u,  ql_u,  w_u,  a_u, &  ! updr (in) 
                             thl_u2, qt_u2, ql_u2, w_u2, a_u2, & ! updr (out)  <- updraft value at t+dt_real (for memory)
                             thl_e,  qt_e,  ql_e,  w_e,  a_e,  u_e,  v_e,  tke_e, & ! envr, now as output
                             thl_g,  qt_g,  ql_g,  w_g,        u_g,  v_g,  tke_g, & ! grid-mean as input
                             rho_g, zstar, ustar, wthl_surf, wq_surf, wstar, oblength, &  ! oblength added 09/25/2017
                             pfull, phalf, zfull, zhalf_tmp, dt_real,     & 
                             thl_ut2, qt_ut2, ql_ut2, w_ut2, a_ut2, & ! updr_out_tot (bulk updraft values)
                             m_ut2, buoy_ut2, chic_ut2, entL_ut2, entB_ut2, detL_ut2, detB_ut2, & ! updr_out_tot (bulk updraft values)
                             entL_u2, entB_u2, detL_u2, detB_u2, entr_u2, detr_u2, buoy_u2, &  ! individual updraft values
                             tend_thl_ut2, tend_qt_ut2, tend_ql_ut2, istep)
       endif
    
       call compute_updr( thl_u,  qt_u,  ql_u,  w_u,  a_u, &  ! updr (in) 
                          thl_u1, qt_u1, ql_u1, w_u1, a_u1, & ! updr (out)  <- updraft value at t+delta_t (for GCM timestepping)
                          thl_e,  qt_e,  ql_e,  w_e,  a_e,  u_e,  v_e,  tke_e, & ! envr, now as output
                          thl_g,  qt_g,  ql_g,  w_g,        u_g,  v_g,  tke_g, & ! grid-mean as input
                          rho_g, zstar, ustar, wthl_surf, wq_surf, wstar, oblength, &  ! oblength added 09/25/2017
                          pfull, phalf, zfull, zhalf_tmp, delta_t,     & 
                          thl_ut, qt_ut, ql_ut, w_ut, a_ut, & ! updr_out_tot (bulk updraft values)
                          m_ut, buoy_ut, chic_ut, entL_ut, entB_ut, detL_ut, detB_ut, & ! updr_out_tot (bulk updraft values)
                          entL_u, entB_u, detL_u, detB_u, entr_u, detr_u, buoy_u, &  ! individual updraft values
                          tend_thl_ut, tend_qt_ut, tend_ql_ut, istep) 
    endif ! (do_compute_updr_lbl)
      
                       
    ! Compute the mass-flux part of the EDMF tendency
    call compute_tend_mf (thl_u1,  qt_u1,  ql_u1,  w_u1,  a_u1, &
                          thl_e,  qt_e,  ql_e,  w_e,  rho_g,  w_g, zhalf_tmp, &
                          tend_thl_mf, tend_qt_mf, tend_ql_mf, flux_thl_mf, flux_qt_mf, flux_ql_mf, net_mf_tot) !
    
    ! Update the grid-mean variables wT, fT
    call update_grid_mean_mf (u_g, v_g, w_g, tke_g, thl_g, qt_g, ql_g, pfull, delta_t, &
                              tend_thl_ut, tend_qt_ut, tend_ql_ut, tend_thl_mf, tend_qt_mf, tend_ql_mf, &
                              u_st1, v_st1, w_st1, tke_st1, thl_st1, qt_st1, ql_st1, &
                              t_st1, q_st1, thv_st1, rho_st1)
                          
    ! Compute decomp again to get updated environmental parameters fn, wn, an
    call envupd_decompose(thl_st1, qt_st1, ql_st1, w_st1,        u_st1, v_st1, tke_st1,  &
                          thl_u1,   qt_u1,  ql_u1,  w_u1,  a_u1,                         & 
                          thl_en1, qt_en1, ql_en1, w_en1, a_en1, u_en1, v_en1, tke_en1)
    
    if (do_mf_env_implicit) then
        call mf_env_implicit (thl_e, qt_e, ql_e, pfull, zhalf_tmp, rho_g, delta_t,       &
                              thl_en1, qt_en1, ql_en1, a_en1,                          &
                              thl_st1, qt_st1, ql_st1, t_st1, q_st1, thv_st1, rho_st1, &
                              tend_thl_mf, tend_qt_mf, tend_ql_mf, flux_thl_mf, flux_qt_mf, flux_ql_mf, net_mf_tot)
    end if    
    
    ! Modify surface flux to account for the updated grid-mean surface air condition                  
    call update_surf_flux_mf (t_g, q_g, t_st1, q_st1, dhdt_atm, dedq_atm, flux_t, flux_q)
    
    ! Added ZTAN 11/21/2017: keep a copy of previous updraft values for use in do_phiu_corr
    thl_uo = thl_u; qt_uo = qt_u; ql_uo = ql_u; w_uo = w_u; a_uo = a_u
    
    
    if (delta_t .gt. dt_real) then 
       ! Leapfrog: use updraft value at t+dt_real for memory ('u2' terms)
       thl_u = thl_u2; qt_u = qt_u2; ql_u = ql_u2; w_u = w_u2; a_u = a_u2; buoy_o_u = buoy_u2
    else
       ! Simple forward: use updraft value at t+delta_t for memory ('u1' terms)
       thl_u = thl_u1; qt_u = qt_u1; ql_u = ql_u1; w_u = w_u1; a_u = a_u1; buoy_o_u = buoy_u
    endif
    
    ! --- END MF COMPUTATION --- !
    ! 'en1', 'st1' values are evaluated after delta_t (= 2*dt_real if leapfrog, = dt_real if NOT leapfrog)
    ! dt_real SHOULD NOT BE USED afterwards (except for mixed_layer SST update)! 

    
    ! Compute ED part of the EDMF tendency
    ! 'en2' values are evaluated after delta_t
    call compute_tend_ed (  is, ie, js, je, Time,                                         & 
                            thl_en1, qt_en1, ql_en1, w_en1, a_en1, u_en1, v_en1, tke_en1, &
                            pfull, zfull, phalf, dt_real, delta_t,                        &
                            ustar, wstar, zstar, oblength,                                &
                            k_m, k_t, t_surf, flux_t, flux_q, flux_u, flux_v, flux_r,     &
                            dtaudv_atm, dhdt_surf, dedt_surf, dedq_surf, drdt_surf,       &
                            dhdt_atm, dedq_atm, net_surf_sw_down, surf_lw_down, Tri_surf, &
                            ls, dls, diss_heat, rho_half,                                 &
                            thl_en2, qt_en2, ql_en2, w_en2, a_en2, u_en2, v_en2, tke_en2, &
                            flux_thl_ed, flux_qt_ed, flux_ql_ed, flux_u_ed, flux_v_ed,    &
                            tend_thl_ed, tend_qt_ed, tend_ql_ed, tend_u_ed, tend_v_ed,    &
                            edmf_update_ml, istep )
    
    ! Compute Cloud/Precip tendency
    ! 'en3' values are also evaluated after delta_t (as 'en2')
    call compute_tend_cloud (thl_en2, qt_en2, ql_en2, w_en2, a_en2, u_en2, v_en2, tke_en2,  &
                             thl_en3, qt_en3, ql_en3, w_en3, a_en3, u_en3, v_en3, tke_en3,  &
                             delta_t, &  ! dt_real is removed
                             pfull, zfull, k_t, ls, dls, rho_half, flux_thl_ed, flux_qt_ed, &
                             rh, s, sigmas, ccov, wthl, wqt, tend_thl_cld, tend_qt_cld, tend_ql_cld)
    
    ! Sum up all tendencies, and generate diagnostic output
    call compute_tend_grid (              tend_qt_clip1, tend_ql_clip1, tend_ql_clip2,     &
                            tend_thl_mf,  tend_qt_mf,   tend_ql_mf,          &
                            tend_thl_ut,  tend_qt_ut,   tend_ql_ut,          &
                            tend_thl_ed,  tend_qt_ed,   tend_ql_ed, tend_u_ed,  tend_v_ed, &
                            tend_thl_cld, tend_qt_cld,  tend_ql_cld,         &                   
                            pfull, phalf, ccov, a_en3, ql_u1, a_u1,         &
                            tend_thl_tot, tend_qt_tot, tend_ql_tot,         & 
                            tend_u_tot, tend_v_tot, tend_t_tot, tend_q_tot, &
                            ccov_ed, ccov_mf, precip_ed, precip_mf, rain2d_ed, rain2d_mf )
    
    ! Compute TKE tendency, new TKE is given by tke_en4
    
    call compute_surface_tke(ustar, wstar, tke_surf) 
    ! tke_surf(:,:) = max( (3.75*(ustar(:,:)**2.) + 0.2*(wstar(:,:)**2.)), 1.0e-8)
    call compute_tend_tke  (thl_en3, qt_en3, ql_en3, w_en3, a_en3, u_en3, v_en3, tke_en3,      &
                            w_st1, tke_st1, flux_thl_ed, flux_qt_ed, flux_u_ed, flux_v_ed,     & 
                            pfull, phalf, zfull, zhalf_tmp, tke_surf, rho_g, ccov, dls, k_m, entL_ut, entB_ut, &
                            detL_u, detB_u, w_u1, net_mf_tot, zstar, delta_t,                 & ! dt_real is removed
                            flux_tke_mf, flux_tke_ed, &
                            tend_tke_buoy, tend_tke_buoyD, tend_tke_buoyW, tend_tke_shear,   &
                            tend_tke_detr, tend_tke_entr, tend_tke_diss, tend_tke_corr,      &
                            tend_tke_source, tend_tke_mf, tend_tke_ed, tend_tke_tot, tke_en4 )
    
    ! ed_leap_frog removed by ZTAN -- 09/10/2017
    thldel = tend_thl_tot * delta_t
     qtdel = tend_qt_tot  * delta_t     
      tdel = tend_t_tot   * delta_t
      qdel = tend_q_tot   * delta_t
    liqdel = tend_ql_tot  * delta_t
    tkedel = tend_tke_tot * delta_t
      udel = tend_u_tot   * delta_t
      vdel = tend_v_tot   * delta_t
    
    call envupd_combine  (thl_u1,   qt_u1,   ql_u1,   w_u1,   a_u1,                            & 
                          thl_en3,  qt_en3,  ql_en3,  w_en3,         u_en3,  v_en3,  tke_en4,  &
                          thl_g_diag, qt_g_diag, ql_g_diag, w_g_diag, u_g_diag, v_g_diag, tke_g_diag)
    
    call cloud_decompose (thl_g_diag, qt_g_diag, ql_g_diag, pfull, t_g_diag, q_g_diag)  ! NOTE ZTAN 2017/09/11: may not be exact!

    thldel_diag = thl_g_diag - thl_g
    qtdel_diag  = qt_g_diag  - qin - liqin
    tdel_diag   = t_g_diag   - tin
    qdel_diag   = q_g_diag   - qin
    liqdel_diag = ql_g_diag  - liqin
    tkedel_diag = tke_g_diag - tkein
    udel_diag   = u_g_diag   - uin
    vdel_diag   = v_g_diag   - vin
    

      
      ! DIAG OUTPUT
      if ( id_thl_e > 0 ) used = send_data ( id_thl_e, thl_en3, Time)
      if ( id_qt_e > 0 )  used = send_data ( id_qt_e, qt_en3, Time)
      if ( id_ql_e > 0 )  used = send_data ( id_ql_e, ql_en3, Time)
      if ( id_tke_e > 0 )  used = send_data ( id_tke_e, tke_en3, Time)
      if ( id_a_e > 0 )  used = send_data ( id_a_e, a_en3, Time)
      
      if ( id_thl_u > 0 ) used = send_data ( id_thl_u, thl_ut, Time)
      if ( id_qt_u > 0 )  used = send_data ( id_qt_u, qt_ut, Time)
      if ( id_ql_u > 0 )  used = send_data ( id_ql_u, ql_ut, Time)
      if ( id_w_u > 0 )  used = send_data ( id_w_u, w_ut, Time)
      if ( id_a_u > 0 )  used = send_data ( id_a_u, a_ut, Time)
      if ( id_m_u > 0 )  used = send_data ( id_m_u, m_ut, Time)
      if ( id_b_u > 0 )  used = send_data ( id_b_u, buoy_ut, Time)
      
      if ( id_k_t > 0 )  used = send_data ( id_k_t, k_t, Time)
      if ( id_k_m > 0 )  used = send_data ( id_k_m, k_m, Time)
      if ( id_ls > 0 )  used = send_data ( id_ls, ls, Time) 
      if ( id_ustar > 0 )  used = send_data ( id_ustar, ustar, Time) 
      if ( id_wstar > 0 )  used = send_data ( id_wstar, wstar, Time) 
      if ( id_bstar > 0 )  used = send_data ( id_bstar, bstar, Time) 
      if ( id_zstar > 0 )  used = send_data ( id_zstar, zstar, Time)
          
      if ( id_flux_thl_ed > 0 )  used = send_data ( id_flux_thl_ed, flux_thl_ed, Time) 
      if ( id_flux_qt_ed > 0 )  used = send_data ( id_flux_qt_ed, flux_qt_ed, Time) 
      if ( id_flux_ql_ed > 0 )  used = send_data ( id_flux_ql_ed, flux_ql_ed, Time) 
      if ( id_flux_tke_ed > 0 )  used = send_data ( id_flux_tke_ed, flux_tke_ed, Time) 
      if ( id_flux_thl_mf > 0 )  used = send_data ( id_flux_thl_mf, flux_thl_mf, Time) 
      if ( id_flux_qt_mf > 0 )  used = send_data ( id_flux_qt_mf, flux_qt_mf, Time) 
      if ( id_flux_ql_mf > 0 )  used = send_data ( id_flux_ql_mf, flux_ql_mf, Time) 
      if ( id_flux_tke_mf > 0 )  used = send_data ( id_flux_tke_mf, flux_tke_mf, Time) 
      
      if ( id_rh > 0 )  used = send_data ( id_rh, rh, Time)
      if ( id_s > 0 )  used = send_data ( id_s, s, Time)
      if ( id_sigmas > 0 )  used = send_data ( id_sigmas, sigmas, Time)
      if ( id_ccov_ed > 0 )  used = send_data ( id_ccov_ed, ccov_ed, Time)
      if ( id_ccov_mf > 0 )  used = send_data ( id_ccov_mf, ccov_mf, Time)
      if ( id_precip_ed > 0 )  used = send_data ( id_precip_ed, precip_ed, Time)
      if ( id_precip_mf > 0 )  used = send_data ( id_precip_mf, precip_mf, Time)
      if ( id_rain2d_ed > 0 )  used = send_data ( id_rain2d_ed, rain2d_ed, Time)
      if ( id_rain2d_mf > 0 )  used = send_data ( id_rain2d_mf, rain2d_mf, Time)
      
      if ( id_entT > 0) used = send_data( id_entT, entL_ut + entB_ut, Time)
      if ( id_detT > 0) used = send_data( id_detT, detL_ut + detB_ut, Time)

      ! ZTAN 11/15/2017: TKE source/sink terms output
      if ( id_tend_tke_buoy  > 0 ) used = send_data(id_tend_tke_buoy,  tend_tke_buoy,  Time)
      if ( id_tend_tke_shear > 0 ) used = send_data(id_tend_tke_shear, tend_tke_shear, Time)
      if ( id_tend_tke_detr  > 0 ) used = send_data(id_tend_tke_detr,  tend_tke_detr,  Time)
      if ( id_tend_tke_entr  > 0 ) used = send_data(id_tend_tke_entr,  tend_tke_entr,  Time)
      if ( id_tend_tke_diss  > 0 ) used = send_data(id_tend_tke_diss,  tend_tke_diss,  Time)
      if ( id_tend_tke_corr  > 0 ) used = send_data(id_tend_tke_corr,  tend_tke_corr,  Time)
      if ( id_tend_tke_mf    > 0 ) used = send_data(id_tend_tke_mf,    tend_tke_mf,    Time)
      if ( id_tend_tke_ed    > 0 ) used = send_data(id_tend_tke_ed,    tend_tke_ed,    Time)
      if ( id_tend_tke_tot   > 0 ) used = send_data(id_tend_tke_tot,   tend_tke_tot,   Time)

      ! SLICE OUTPUT -- After One Step, 
      ! ** THIS ONLY WORKS FOR SINGLE-CPU RUNS **

      if ( id_t_g_sl > 0 ) used = send_data ( id_t_g_sl, t_g(1,:,:), Time) 
      if ( id_q_g_sl > 0 ) used = send_data ( id_q_g_sl, q_g(1,:,:), Time) 
      if ( id_ql_g_sl > 0 ) used = send_data ( id_ql_g_sl, ql_g(1,:,:), Time) 
      if ( id_u_g_sl > 0 ) used = send_data ( id_u_g_sl, u_g(1,:,:), Time) 
      if ( id_v_g_sl > 0 ) used = send_data ( id_v_g_sl, v_g(1,:,:), Time) 
      if ( id_thl_g_sl > 0 ) used = send_data ( id_thl_g_sl, thl_g(1,:,:), Time) 
      if ( id_qt_g_sl > 0 ) used = send_data ( id_qt_g_sl, qt_g(1,:,:), Time) 
      if ( id_thv_g_sl > 0 ) used = send_data ( id_thv_g_sl, thv_g(1,:,:), Time) 
      if ( id_rho_g_sl > 0 ) used = send_data ( id_rho_g_sl, rho_g(1,:,:), Time) 
      if ( id_tke_g_sl > 0 ) used = send_data ( id_tke_g_sl, tke_g(1,:,:), Time) 
      if ( id_zfull_sl > 0 ) used = send_data ( id_zfull_sl, zfull(1,:,:), Time) 
      if ( id_zhalf_sl > 0 ) used = send_data ( id_zhalf_sl, zhalf_tmp(1,:,:), Time) 
      if ( id_pfull_sl > 0 ) used = send_data ( id_pfull_sl, pfull(1,:,:), Time) 
      if ( id_phalf_sl > 0 ) used = send_data ( id_phalf_sl, phalf(1,:,:), Time) 
       
      if ( id_rh_sl > 0 ) used = send_data ( id_rh_sl, rh(1,:,:), Time) 
      if ( id_ccov_ed_sl > 0 )  used = send_data ( id_ccov_ed_sl, ccov_ed(1,:,:), Time)
      if ( id_ccov_mf_sl > 0 )  used = send_data ( id_ccov_mf_sl, ccov_mf(1,:,:), Time)
      if ( id_rain2d_ed_sl > 0 )  used = send_data ( id_rain2d_ed_sl, rain2d_ed(1,:), Time)
      if ( id_rain2d_mf_sl > 0 )  used = send_data ( id_rain2d_mf_sl, rain2d_mf(1,:), Time)
      if ( id_flux_t_sl > 0 )  used = send_data ( id_flux_t_sl, flux_t(1,:), Time)
      if ( id_flux_q_sl > 0 )  used = send_data ( id_flux_q_sl, flux_q(1,:), Time)
      
      if ( id_thl_e_sl > 0 ) used = send_data ( id_thl_e_sl, thl_en3(1,:,:), Time)  ! <- Do transpose?? ZTAN 09/13/2017
      if ( id_qt_e_sl > 0 )  used = send_data ( id_qt_e_sl, qt_en3(1,:,:), Time)
      if ( id_ql_e_sl > 0 )  used = send_data ( id_ql_e_sl, ql_en3(1,:,:), Time)
      if ( id_tke_e_sl > 0 )  used = send_data ( id_tke_e_sl, tke_en3(1,:,:), Time)
      if ( id_a_e_sl > 0 )  used = send_data ( id_a_e_sl, a_en3(1,:,:), Time)
      
      if ( id_thl_u_sl > 0 ) used = send_data ( id_thl_u_sl, thl_ut(1,:,:), Time)
      if ( id_qt_u_sl > 0 )  used = send_data ( id_qt_u_sl, qt_ut(1,:,:), Time)
      if ( id_ql_u_sl > 0 )  used = send_data ( id_ql_u_sl, ql_ut(1,:,:), Time)
      if ( id_w_u_sl > 0 )  used = send_data ( id_w_u_sl, w_ut(1,:,:), Time)
      if ( id_a_u_sl > 0 )  used = send_data ( id_a_u_sl, a_ut(1,:,:), Time)
      if ( id_m_u_sl > 0 )  used = send_data ( id_m_u_sl, m_ut(1,:,:), Time)
      if ( id_b_u_sl > 0 )  used = send_data ( id_b_u_sl, buoy_ut(1,:,:), Time)
      
      if ( id_entr_u_sl > 0 )  used = send_data ( id_entr_u_sl, entr_u(1,:,:,1), Time) ! only send the first updraft
      if ( id_detr_u_sl > 0 )  used = send_data ( id_detr_u_sl, detr_u(1,:,:,1), Time) ! only send the first updraft
      if ( id_entT_sl > 0 )  used = send_data ( id_entT_sl, entL_ut(1,:,:) + entB_ut(1,:,:), Time)
      if ( id_detT_sl > 0 )  used = send_data ( id_detT_sl, detL_ut(1,:,:) + detB_ut(1,:,:), Time)
      
      if ( id_k_t_sl > 0 )  used = send_data ( id_k_t_sl, k_t(1,:,:), Time)
      if ( id_k_m_sl > 0 )  used = send_data ( id_k_m_sl, k_m(1,:,:), Time)
      if ( id_ls_sl > 0 )  used = send_data ( id_ls_sl, ls(1,:,:), Time) 
      if ( id_ustar_sl > 0 )  used = send_data ( id_ustar_sl, ustar(1,:), Time) 
      if ( id_wstar_sl > 0 )  used = send_data ( id_wstar_sl, wstar(1,:), Time) 
      if ( id_bstar_sl > 0 )  used = send_data ( id_bstar_sl, bstar(1,:), Time) 
      if ( id_zstar_sl > 0 )  used = send_data ( id_zstar_sl, zstar(1,:), Time)
          
      if ( id_flux_thl_ed_sl > 0 )  used = send_data ( id_flux_thl_ed_sl, flux_thl_ed(1,:,:), Time) 
      if ( id_flux_qt_ed_sl > 0 )  used = send_data ( id_flux_qt_ed_sl, flux_qt_ed(1,:,:), Time) 
      if ( id_flux_ql_ed_sl > 0 )  used = send_data ( id_flux_ql_ed_sl, flux_ql_ed(1,:,:), Time) 
      if ( id_flux_thl_mf_sl > 0 )  used = send_data ( id_flux_thl_mf_sl, flux_thl_mf(1,:,:), Time) 
      if ( id_flux_qt_mf_sl > 0 )  used = send_data ( id_flux_qt_mf_sl, flux_qt_mf(1,:,:), Time) 
      if ( id_flux_ql_mf_sl > 0 )  used = send_data ( id_flux_ql_mf_sl, flux_ql_mf(1,:,:), Time) 

      ! ZTAN: debug output -- diagnose tendencies. 09/14/2017 -> Changed to 1/dt tendencies
      if ( id_thldel_sl > 0 ) used = send_data ( id_thldel_sl, tend_thl_tot(1,:,:), Time)
      if ( id_qtdel_sl > 0 )  used = send_data ( id_qtdel_sl,  tend_qt_tot(1,:,:),  Time)
      if ( id_tdel_sl > 0 )   used = send_data ( id_tdel_sl,   tend_t_tot(1,:,:),   Time) 
      if ( id_qdel_sl > 0 )   used = send_data ( id_qdel_sl,   tend_q_tot(1,:,:),   Time) 
      if ( id_liqdel_sl > 0 ) used = send_data ( id_liqdel_sl, tend_ql_tot(1,:,:),  Time) 
      if ( id_udel_sl > 0 )   used = send_data ( id_udel_sl,   tend_u_tot(1,:,:),   Time) 
      if ( id_vdel_sl > 0 )   used = send_data ( id_vdel_sl,   tend_v_tot(1,:,:),   Time) 
      if ( id_tkedel_sl > 0 ) used = send_data ( id_tkedel_sl, tend_tke_tot(1,:,:), Time)

      ! ZTAN: DEBUG output -- diagnose ERROR IN tendencies... 09/14/2017
      if ( id_thldel_diag_sl > 0 ) used = send_data ( id_thldel_diag_sl, thldel(1,:,:) - thldel_diag(1,:,:), Time)
      if ( id_qtdel_diag_sl > 0 )  used = send_data ( id_qtdel_diag_sl,  qtdel(1,:,:)  - qtdel_diag(1,:,:),  Time)
      if ( id_tdel_diag_sl > 0 )   used = send_data ( id_tdel_diag_sl,   tdel(1,:,:)   - tdel_diag(1,:,:),   Time) 
      if ( id_qdel_diag_sl > 0 )   used = send_data ( id_qdel_diag_sl,   qdel(1,:,:)   - qdel_diag(1,:,:),   Time) 
      if ( id_liqdel_diag_sl > 0 ) used = send_data ( id_liqdel_diag_sl, liqdel(1,:,:) - liqdel_diag(1,:,:), Time) 
      if ( id_udel_diag_sl > 0 )   used = send_data ( id_udel_diag_sl,   udel(1,:,:)   - udel_diag(1,:,:),   Time) 
      if ( id_vdel_diag_sl > 0 )   used = send_data ( id_vdel_diag_sl,   vdel(1,:,:)   - vdel_diag(1,:,:),   Time) 
      if ( id_tkedel_diag_sl > 0 ) used = send_data ( id_tkedel_diag_sl, tkedel(1,:,:) - tkedel_diag(1,:,:), Time)
      
      ! ZTAN 11/15/2017: TKE source/sink terms output
      if ( id_tend_tke_buoy_sl  > 0 ) used = send_data(id_tend_tke_buoy_sl,  tend_tke_buoy(1,:,:),  Time)
      if ( id_tend_tke_shear_sl > 0 ) used = send_data(id_tend_tke_shear_sl, tend_tke_shear(1,:,:), Time)
      if ( id_tend_tke_detr_sl  > 0 ) used = send_data(id_tend_tke_detr_sl,  tend_tke_detr(1,:,:),  Time)
      if ( id_tend_tke_entr_sl  > 0 ) used = send_data(id_tend_tke_entr_sl,  tend_tke_entr(1,:,:),  Time)
      if ( id_tend_tke_diss_sl  > 0 ) used = send_data(id_tend_tke_diss_sl,  tend_tke_diss(1,:,:),  Time)
      if ( id_tend_tke_corr_sl  > 0 ) used = send_data(id_tend_tke_corr_sl,  tend_tke_corr(1,:,:),  Time)
      if ( id_tend_tke_mf_sl    > 0 ) used = send_data(id_tend_tke_mf_sl,    tend_tke_mf(1,:,:),    Time)
      if ( id_tend_tke_ed_sl    > 0 ) used = send_data(id_tend_tke_ed_sl,    tend_tke_ed(1,:,:),    Time)
      if ( id_tend_tke_tot_sl   > 0 ) used = send_data(id_tend_tke_tot_sl,   tend_tke_tot(1,:,:),   Time)

      ! END OF SLICE OUTPUT


      ! Added ZTAN 11/21/2017: redefine thl_e, qt_e, ql_e to be at t + dt_real (instead of t + delta_t)
      !                        ** Only for do_phiu_corr **
      thl_e(:,:,:) = thl_en3(:,:,:) + &
                     (delta_t - dt_real)/delta_t *(thl_e(:,:,:) - thl_en3(:,:,:))
      qt_e(:,:,:)  = qt_en3(:,:,:)  + &
                     (delta_t - dt_real)/delta_t *(qt_e(:,:,:)  - qt_en3(:,:,:))
      ql_e(:,:,:)  = ql_en3(:,:,:)  + &
                     (delta_t - dt_real)/delta_t *(ql_e(:,:,:)  - ql_en3(:,:,:))

      ! Added ZTAN 11/21/2017: Define 'updraft values' where a_u is small and the updraft 
      !                        does not exist,which serve as predicted updraft values. 
      !                        In case that updrafts reach these levels at the next step, 
      !                        these 'updraft values' will be used to compute updraft buoyancy. 
      call do_phiu_corr (thl_u,  qt_u,  ql_u,  a_u, &
                         thl_uo, qt_uo, ql_uo,      &
                         thl_e,  qt_e,  ql_e )
      
      if (do_trans_dphi) then ! transform back to the departure of updraft from env
          do iupd = 1, updraft_number
              thl_u(:,:,:,iupd) = thl_u(:,:,:,iupd) - thl_e(:,:,:)
              qt_u (:,:,:,iupd) = qt_u (:,:,:,iupd) - qt_e (:,:,:)
              ql_u (:,:,:,iupd) = ql_u (:,:,:,iupd) - ql_e (:,:,:)
          end do
      end if
         
end subroutine newedmf


subroutine newedmf_end(cloud_cover, mfcloud_cover)
real, intent(inout), dimension(:,:,:) :: cloud_cover, mfcloud_cover
character(len=64) :: file

file='RESTART/newedmf.res'
call nullify_domain()

call write_data(trim(file), 'cloud_cover', cloud_cover, grid_domain)
call write_data(trim(file), 'mfcloud_cover', mfcloud_cover, grid_domain)
! a_u, w_u, thl_u, qt_u, and ql_u are stored variables in the module, 
! not the input of this subroutine.
call write_data(trim(file), 'a_u',   a_u,   grid_domain)
call write_data(trim(file), 'w_u',   w_u,   grid_domain)
call write_data(trim(file), 'thl_u', thl_u, grid_domain)
call write_data(trim(file), 'qt_u',  qt_u,  grid_domain)
call write_data(trim(file), 'ql_u',  ql_u,  grid_domain)


end subroutine newedmf_end




! ------------------------- !
!        UPDATE TKE         !
! ------------------------- !

subroutine compute_tend_tke (thl_env, qt_env, ql_env, w_env, a_env, u_env, v_env, tke_env,    &
                             w_grid, tke_grid, flux_thl_ed, flux_qt_ed, flux_u_ed, flux_v_ed, & 
                             pfull, phalf, zfull, zhalf, tke_surf, rhot, ccov, dls, k_m, entL, entB,    &
                             detL_u, detB_u, w_upd, net_mf_tot, zstar, delta_t,               &
                             flux_tke_mf, flux_tke_ed, &
                             tend_tke_buoy, tend_tke_buoyD, tend_tke_buoyW, tend_tke_shear,   &
                             tend_tke_detr, tend_tke_entr, tend_tke_diss, tend_tke_corr,      &
                             tend_tke_source, tend_tke_mf, tend_tke_ed, tend_tke_tot, tke_envn )

    real   , intent(in) ,   dimension(:,:,:)   :: thl_env, qt_env, ql_env, w_env, a_env, u_env, v_env, tke_env
    real   , intent(in) ,   dimension(:,:,:)   :: w_grid, tke_grid
    real   , intent(in) ,   dimension(:,:,:)   :: flux_thl_ed, flux_qt_ed, flux_u_ed, flux_v_ed
    real   , intent(in) ,   dimension(:,:,:)   :: pfull, phalf, zfull, zhalf, rhot
    real   , intent(in) ,   dimension(:,:)     :: tke_surf
    real   , intent(in) ,   dimension(:,:,:)   :: ccov, dls, k_m, entL, entB, net_mf_tot
    real   , intent(in) ,   dimension(:,:,:,:) :: detL_u, detB_u, w_upd !, a_upd
    real   , intent(in) ,   dimension(:,:)     :: zstar
    real   , intent(in)                        :: delta_t  ! dt_real is removed
    real   , intent(out),   dimension(:,:,:)   :: flux_tke_mf, flux_tke_ed, tke_envn ! tke_envn Added to output
    real   , intent(out),   dimension(:,:,:)   :: tend_tke_buoy, tend_tke_buoyD, tend_tke_buoyW, tend_tke_shear, &
                                                  tend_tke_detr, tend_tke_entr, tend_tke_diss, tend_tke_corr, &
                                                  tend_tke_source, tend_tke_mf, tend_tke_ed, tend_tke_tot 
    
    ! And: total_tend, envr_new, updr_new, updr_out, edmf_flux, grid_mean
    
    real   , dimension(size(thl_env,1),size(thl_env,2),size(thl_env,3)) :: t_env, q_env, th_env, &
            alphaD, betaD, alphaW, betaW, qsat_env, dqsat_env, dudz, dvdz, dz
    
    real   , dimension(size(thl_env,1),size(thl_env,2),size(thl_env,3)) :: & 
            rwthl_f, rwqt_f, rwu_f, rwv_f  ! Mass and area weighted ED fluxes at full level
            
    real   , dimension(size(thl_env,1),size(thl_env,2),size(thl_env,3)+1) :: & 
            u_half, v_half  ! U and V at half levels
    
    real   , dimension(size(thl_env,1),size(thl_env,2),size(thl_env,3))   :: &  
            pbuoy, pbuoyD, pbuoyW, pshear, pentr, pdetr, pdiss, psource, pmf, pcorr, ped ! TKE source terms
    
    integer :: nlev, ilev, nupd, iupd
    nlev = size(thl_env,3)
    nupd = size(w_upd, 4)
    
    dz(:,:,:) = zhalf(:,:,2:nlev+1) - zhalf(:,:,1:nlev)
    
    ! NOTE1: edmf_flux are all given in density-weighted forms -->
    ! pbuoy, pshear in (kg/m^3)*(m^2/s^2)/s: a_n already included in flux wT, wq, wu, wv
    
    call cloud_decompose(thl_env, qt_env, ql_env, pfull, t_env, q_env)
    th_env(:,:,:) = t_env(:,:,:) /(pfull(:,:,:)/pstd_mks)**kappa
    
    ! --- Buoyancy production of TKE --- ! 
    rwthl_f(:,:,1:nlev) = (flux_thl_ed(:,:,1:nlev) + flux_thl_ed(:,:,2:nlev+1))*0.5
    rwqt_f (:,:,1:nlev) = (flux_qt_ed (:,:,1:nlev) + flux_qt_ed (:,:,2:nlev+1))*0.5
    
    alphaD(:,:,:) = 1. + d608 * qt_env(:,:,:)
    betaD (:,:,:) =      d608 * th_env(:,:,:)

    call compute_qsat(t_env, pfull, qsat_env, dqsat_env)
    
    alphaW(:,:,:) = (1.-qt_env + (1.+d608) * qsat_env * (1. + d622*HLv/(rdgas*t_env))) / & 
        (1. + d622*HLv*HLv*qsat_env/(rdgas*Cp_air*t_env*t_env))
    betaW (:,:,:) = th_env * (HLv/Cp_air/t_env*alphaW - 1.)
    
    ! pbuoy (unit: kg/m^3 * m^2/s^2 * s^-1)
    pbuoyD = (alphaD * rwthl_f + betaD * rwqt_f) * Grav / th_env
    pbuoyW = (alphaW * rwthl_f + betaW * rwqt_f) * Grav / th_env
    pbuoy = (1.-ccov) * pbuoyD + ccov * pbuoyW
    
    if (buoy_freeatm == .false.) then
        do ilev = 1, nlev
            where (zfull(:,:,ilev) .gt. zstar(:,:))
                pbuoy(:,:,ilev) = 0.
            end where
        end do
    end if
    
    ! --- Shear production of TKE --- ! 
    rwu_f(:,:,1:nlev) = (flux_u_ed(:,:,1:nlev) + flux_u_ed(:,:,2:nlev+1))*0.5
    rwv_f(:,:,1:nlev) = (flux_v_ed(:,:,1:nlev) + flux_v_ed(:,:,2:nlev+1))*0.5
    
    u_half(:,:,1) = u_env(:,:,1)
    v_half(:,:,1) = v_env(:,:,1)
    
    u_half(:,:,2:nlev) = (u_env(:,:,1:nlev-1) + u_env(:,:,2:nlev))*0.5
    v_half(:,:,2:nlev) = (v_env(:,:,1:nlev-1) + v_env(:,:,2:nlev))*0.5
     
    u_half(:,:,nlev+1) = 0.
    v_half(:,:,nlev+1) = 0. 
    
    dudz = (u_half(:,:,2:nlev+1) - u_half(:,:,1:nlev)) / dz(:,:,:)
    dvdz = (v_half(:,:,2:nlev+1) - v_half(:,:,1:nlev)) / dz(:,:,:)
    
    ! pshear (unit: kg/m^3 * m^2/s^2 * s^-1)
    pshear(:,:,:) = - dudz(:,:,:) * rwu_f(:,:,:) - dvdz(:,:,:) * rwv_f(:,:,:)
    
    ! --- Entrainment loss of TKE --- ! 
    pentr(:,:,:) = -(entL(:,:,:) + entB(:,:,:))*tke_env(:,:,:)/(-dz(:,:,:))
    
    ! --- Detrainment gain of TKE --- !  ! NOTE ZTAN 2017/09/11 -- 
                                         ! NEED TO CALCULATE IN UPDRAFT PART 
                                         !    use w_upd(endtime) may underestimate
    pdetr(:,:,:) = 0.
    do iupd = 1, nupd
        pdetr(:,:,:) = pdetr(:,:,:) + (detL_u(:,:,:,iupd) + detB_u(:,:,:,iupd)) * &
                       (0.5 * (w_upd(:,:,:,iupd) - w_env(:,:,:))**2.) / (-dz(:,:,:))
    end do
    
    ! --- Dissipation loss of TKE --- !  ! Changed to delta_t: ZTAN 09/11/2017 
    pdiss(:,:,:) = max(tke_env(:,:,:), 0.) * & 
                  (-1.+ (1.+ ce*delta_t/2./max(dls(:,:,:), 1.e-5) * max(tke_env(:,:,:), 0.)**0.5)**(-2.))

    pdiss(:,:,:) = pdiss(:,:,:) * a_env(:,:,:) * rhot(:,:,:) / delta_t
    
    ! --- Total non-EDMF gain/loss of TKE --- !
    psource = pbuoy + pshear + pentr + pdetr + pdiss
    ! All 'p' values include the rho*a_e factor !
    
    ! --- MF (compensating subsidence), using implicit upwind -- ! 
    ! call compute_pmf (w_grid, tke_grid, rhot, a_env, dz, w_upd, a_upd, pbuoy, dt_real, &
    !                   psource, pmf, pcorr, flux_tke_mf)
    call compute_pmf (tke_grid, rhot, a_env, dz, net_mf_tot, pbuoy, delta_t, &  ! dt_real Changed to delta_t: ZTAN 09/11/2017
                      psource, pmf, pcorr, flux_tke_mf)
                          
    ! --- ED, using k_m --- !
    call compute_ped (t_env, q_env, a_env, tke_env, psource, pmf, rhot, k_m, &
                      pfull, phalf, zfull, delta_t, tke_surf, ped, flux_tke_ed, tke_envn)
    
    ! --- Generate outputs --- !
    tend_tke_buoy  = pbuoy  / rhot
    tend_tke_buoyD = pbuoyD / rhot
    tend_tke_buoyW = pbuoyW / rhot
    tend_tke_shear = pshear / rhot
    tend_tke_detr  = pdetr  / rhot
    tend_tke_entr  = pentr  / rhot
    tend_tke_diss  = pdiss  / rhot
    tend_tke_corr  = pcorr  / rhot
    tend_tke_source = psource / rhot
    tend_tke_mf    = pmf    / rhot
    tend_tke_ed    = ped    / rhot
    tend_tke_tot   = (psource + pmf + ped) / rhot
    
    ! Note: tend_tke_tot = tend_tke_source + tend_tke_mf + tend_tke_ed
    !       tend_tke_source = tend_tke_buoy + tend_tke_shear + tend_tke_entr + 
    !                         tend_tke_detr + tend_tke_diss  + tend_tke_corr    (EXCEPT at the SFC)
    
end subroutine compute_tend_tke


subroutine compute_ped (t_env, q_env, a_env, tke_env, psource, pmf, rhot, k_m, &
                        pfull, phalf, zfull, delta_t, tke_surf, ped, flux_tke_ed, tke_n)

    real   , intent(in) ,   dimension(:,:,:)   :: t_env, q_env, a_env, tke_env
    real   , intent(in) ,   dimension(:,:,:)   :: pmf, rhot, k_m
    real   , intent(inout), dimension(:,:,:)   :: psource
    real   , intent(in) ,   dimension(:,:,:)   :: pfull, phalf, zfull
    real   , intent(in) ,   dimension(:,:)     :: tke_surf
    real   , intent(in)                        :: delta_t  ! dt_real is removed
    real   , intent(out),   dimension(:,:,:)   :: ped, flux_tke_ed, tke_n ! Added to output
    
    real   , dimension(size(t_env,1),size(t_env,2),size(t_env,3)) :: tke_0, dt_tke, nu ! , &
                                                                     ! dt_tke_0_pre_ed
    integer :: nlev
    
    nlev = size(t_env, 3)
       
    ! ed_leap_frog removed by ZTAN -- 09/10/2017 
        
    dt_tke(:,:,:) = (psource(:,:,:) + pmf(:,:,:))/rhot(:,:,:)/a_env(:,:,:) 
    tke_0(:,:,:) = tke_env(:,:,:)
    
    ! Added ZTAN 06/16/2017: prescribing TKE at surface
    ! dt_tke_0_pre_ed(:,:,:) = (psource(:,:,:) + pmf(:,:,:))/rhot(:,:,:)
    
    ! tke_0(:,:,nlev) = tke_surf(:,:) 
    dt_tke(:,:,nlev) = (tke_surf(:,:)  - tke_0(:,:,nlev))/delta_t
    ! dt_tke_0_pre_ed(:,:,nlev) = dt_tke(:,:,nlev)*a_env(:,:,nlev) 
    psource(:,:,nlev) = dt_tke(:,:,nlev)*a_env(:,:,nlev)*rhot(:,:,nlev) - pmf(:,:,nlev)
    
    if (ed_tke_opt == 1) then
        tke_0(:,:,:) = tke_0(:,:,:) + dt_tke(:,:,:) * delta_t
        dt_tke(:,:,:) = 0.
    end if
    
    call gcm_vert_diff_tke(1, 1, delta_t, t_env, q_env, tke_0, k_m, a_env, &
                           phalf, pfull, zfull, dt_tke, nu)
    
    tke_n(:,:,:) = tke_0(:,:,:) + dt_tke(:,:,:) * delta_t
    
    flux_tke_ed(:,:,1) = 0.
    flux_tke_ed(:,:,2:nlev) = nu(:,:,2:nlev) * (tke_n (:,:,2:nlev) - tke_n (:,:,1:nlev-1))
    flux_tke_ed(:,:,nlev+1) = 0.
    
    if (ed_tke_opt == 1) then
        ! ped(:,:,:) = dt_tke(:,:,:)/a_env(:,:,:)*rhot(:,:,:)
        ped(:,:,:) = dt_tke(:,:,:)*a_env(:,:,:)*rhot(:,:,:)     ! Corrected ZTAN 09/14/2017
    else
        ! ped(:,:,:) = dt_tke(:,:,:)/a_env(:,:,:)*rhot(:,:,:) - (psource(:,:,:) + pmf(:,:,:))
        ped(:,:,:) = dt_tke(:,:,:)*a_env(:,:,:)*rhot(:,:,:) - (psource(:,:,:) + pmf(:,:,:))    ! Corrected ZTAN 09/14/2017
        ! ped(:,:,:) = dt_tke(:,:,:)*a_env(:,:,:)*rhot(:,:,:) - dt_tke_0_pre_ed(:,:,:)*rhot(:,:,:) ! Added ZTAN 09/20/2017
    end if
    
    
end subroutine compute_ped


! CHANGED ZTAN 2017/09/11: eliminate net_mf, a_upd, w_upd, and change net_mf_tot to input
subroutine compute_pmf (tke_grid, rhot, a_env, dz, net_mf_tot, pbuoy, delta_t, &  ! dt_real Changed to delta_t: ZTAN 09/11/2017
                        psource, pmf, pcorr, flux_tke_mf)

    real   , intent(in) ,   dimension(:,:,:)   :: tke_grid, rhot, a_env, dz, net_mf_tot, pbuoy ! w_grid
    ! real   , intent(in) ,   dimension(:,:,:,:) :: w_upd, a_upd
    real   , intent(inout), dimension(:,:,:)   :: psource
    real   , intent(in)                        :: delta_t  ! dt_real
    real   , intent(out),   dimension(:,:,:)   :: pmf, pcorr, flux_tke_mf
    
    ! real   , dimension(size(a_upd,1),size(a_upd,2),size(a_upd,3),size(a_upd,4)) :: net_mf
    real   , dimension(size(a_env,1),size(a_env,2),size(a_env,3)) :: tke_grid_new
    ! real   , dimension(size(a_upd,1),size(a_upd,2),size(a_upd,3)+1) :: net_mf_tot
    real   , dimension(size(a_env,1),size(a_env,2)) :: c1, c2, ct, coefb, coefc, tke_grid_corr
    
    integer :: nlev, ilev !, nupd, iupd
    nlev = size(a_env,3)
    ! nupd = size(a_upd,4) 
    
    ! net_mf_tot (:,:,:) = 0.
    ! do iupd = 1, nupd
    !     net_mf(:,:,:,iupd) = rhot(:,:,:) * a_upd(:,:,:,iupd) * ( w_upd(:,:,:,iupd) - w_grid(:,:,:))
    !     net_mf_tot (:,:,1:nlev) = net_mf_tot (:,:,1:nlev) + net_mf(:,:,:,iupd)
    ! end do
    
    pmf(:,:,:) = 0.; pcorr(:,:,:) = 0.; flux_tke_mf(:,:,:) = 0.
    
    ! HIGHEST LEVEL (ILEV = 1)
    ilev = 1
    ct(:,:) = rhot(:,:,ilev)/delta_t + net_mf_tot(:,:,ilev+1)/a_env(:,:,ilev)/(-dz(:,:,ilev))  ! dt_real changed to delta_t
    c1(:,:) = rhot(:,:,ilev)/delta_t
    tke_grid_new(:,:,ilev) = & 
        (c1(:,:)*tke_grid(:,:,ilev) + psource(:,:,ilev))/ct
    
    pcorr(:,:,ilev) = 0.
    
    where (tke_grid_new(:,:,ilev) .lt. 0.)        ! Correct for negative TKE (all cases)
        pcorr(:,:,ilev) = - tke_grid_new(:,:,ilev) * ct(:,:)
    end where
    
    if ((diff_opt == 2) .and. (pbuoy_corr)) then  ! Correct for positive TKE (only when pbuoy_corr)
        where (tke_grid(:,:,ilev) .gt. 0.)
            coefb(:,:) = -pbuoy(:,:,ilev)/ct/sqrt(max(tke_grid(:,:,ilev), 1.e-4))
            coefc(:,:) =  pbuoy(:,:,ilev)/ct - tke_grid_new(:,:,ilev)
            tke_grid_corr(:,:) = ( -coefb(:,:) + sqrt(coefb(:,:)**2.-4*coefc(:,:)) )/2.
            pcorr(:,:,ilev) = pbuoy(:,:,ilev)*(sqrt(tke_grid_corr(:,:))/sqrt(max(tke_grid(:,:,ilev), 1.e-4)) - 1.)
        end where
    end if
    
    psource(:,:,ilev) = psource(:,:,ilev) + pcorr(:,:,ilev)
    tke_grid_new(:,:,ilev) = & 
        (c1(:,:)*tke_grid(:,:,ilev) + psource(:,:,ilev))/ct
    flux_tke_mf(:,:,ilev) = - net_mf_tot(:,:,ilev+1)*tke_grid_new(:,:,ilev)/a_env(:,:,ilev)
    pmf(:,:,ilev) = flux_tke_mf(:,:,ilev)/(-dz(:,:,ilev))
    
    ! Debug output of tke_grid_new1 is omitted
    
    
    ! OTHER LEVELS (ILEV = 2 : NLEV)
    do ilev = 2, nlev
        ct(:,:) = rhot(:,:,ilev)/delta_t + net_mf_tot(:,:,ilev+1)/a_env(:,:,ilev)/(-dz(:,:,ilev))  ! dt_real changed to delta_t
        c1(:,:) = rhot(:,:,ilev)/delta_t
        c2(:,:) = net_mf_tot(:,:,ilev)/a_env(:,:,ilev-1)/(-dz(:,:,ilev))
        tke_grid_new(:,:,ilev) = & 
            (c1(:,:)*tke_grid(:,:,ilev) + c2(:,:)*tke_grid_new(:,:,ilev-1) + psource(:,:,ilev))/ct

        pcorr(:,:,ilev) = 0.
    
        where (tke_grid_new(:,:,ilev) .lt. 0.)        ! Correct for negative TKE (all cases)
            pcorr(:,:,ilev) = - tke_grid_new(:,:,ilev) * ct(:,:)
        end where
        
        if ((diff_opt == 2) .and. (pbuoy_corr)) then  ! Correct for positive TKE (only when pbuoy_corr)
            where (tke_grid(:,:,ilev) .gt. 0.)
                coefb(:,:) = -pbuoy(:,:,ilev)/ct/sqrt(max(tke_grid(:,:,ilev), 1.e-4))
                coefc(:,:) =  pbuoy(:,:,ilev)/ct - tke_grid_new(:,:,ilev)
                tke_grid_corr(:,:) = ( -coefb(:,:) + sqrt(coefb(:,:)**2.-4*coefc(:,:)) )/2.
                pcorr(:,:,ilev) = pbuoy(:,:,ilev)*(sqrt(tke_grid_corr(:,:))/sqrt(max(tke_grid(:,:,ilev), 1.e-4)) - 1.)
            end where
        end if
    
        psource(:,:,ilev) = psource(:,:,ilev) + pcorr(:,:,ilev)
        tke_grid_new(:,:,ilev) = & 
            (c1(:,:)*tke_grid(:,:,ilev) + c2(:,:)*tke_grid_new(:,:,ilev-1) + psource(:,:,ilev))/ct
        flux_tke_mf(:,:,ilev) = - net_mf_tot(:,:,ilev+1)*tke_grid_new(:,:,ilev)/a_env(:,:,ilev)
        pmf(:,:,ilev) = (flux_tke_mf(:,:,ilev) - flux_tke_mf(:,:,ilev-1))/(-dz(:,:,ilev))
    
    end do
    

end subroutine compute_pmf


! ------------------------- !
!  CLOUD/PRECIP TENDENCIES  !
! ------------------------- !

subroutine compute_tend_grid (              tend_qt_clip1, tend_ql_clip1, tend_ql_clip2,     &
                              tend_thl_mf,  tend_qt_mf,   tend_ql_mf,          &
                              tend_thl_ut,  tend_qt_ut,   tend_ql_ut,          &
                              tend_thl_ed,  tend_qt_ed,   tend_ql_ed, tend_u_ed,  tend_v_ed, &
                              tend_thl_cld, tend_qt_cld,  tend_ql_cld,         &                   
                              pfull, phalf, ccov, a_env,  ql_upd, a_upd,       &
                              tend_thl_tot, tend_qt_tot,  tend_ql_tot,         & 
                              tend_u_tot, tend_v_tot,  tend_t_tot, tend_q_tot, &
                              ccov_ed, ccov_mf, precip_ed, precip_mf, rain2d_ed, rain2d_mf )
                              
    real   , intent(in) ,   dimension(:,:,:)   :: tend_qt_clip1, tend_ql_clip1, tend_ql_clip2  ! Clipping of qt and ql
    real   , intent(in) ,   dimension(:,:,:)   :: tend_thl_mf,  tend_qt_mf,   tend_ql_mf   ! MF Transport
    real   , intent(in) ,   dimension(:,:,:)   :: tend_thl_ut,  tend_qt_ut,   tend_ql_ut   ! MF Rain
    real   , intent(in) ,   dimension(:,:,:)   :: tend_thl_ed,  tend_qt_ed,   tend_ql_ed,  &
                                                  tend_u_ed,  tend_v_ed                    ! ED transport
    real   , intent(in) ,   dimension(:,:,:)   :: tend_thl_cld, tend_qt_cld,  tend_ql_cld  ! ED Rain
    real   , intent(in) ,   dimension(:,:,:)   :: pfull, phalf, ccov, a_env
    real   , intent(in) ,   dimension(:,:,:,:) :: ql_upd, a_upd
        
    real   , intent(out),   dimension(:,:,:)   :: tend_thl_tot, tend_qt_tot, tend_ql_tot, & 
                                                  tend_u_tot, tend_v_tot, tend_t_tot, tend_q_tot
    real   , intent(out),   dimension(:,:,:)   :: ccov_ed, ccov_mf, precip_ed, precip_mf
    real   , intent(out),   dimension(:,:)     :: rain2d_ed, rain2d_mf 
    
    integer :: nlev, ilev, nupd, iupd
    real :: small = 1.0e-10
    
    nlev = size(a_upd, 3)
    nupd = size(a_upd, 4)
    
    tend_thl_tot = tend_thl_mf + tend_thl_ut + tend_thl_ed + tend_thl_cld
    tend_qt_tot  = tend_qt_mf  + tend_qt_ut  + tend_qt_ed  + tend_qt_cld + tend_qt_clip1
    tend_ql_tot  = tend_ql_mf  + tend_ql_ut  + tend_ql_ed  + tend_ql_cld + tend_ql_clip1 + tend_ql_clip2
    tend_u_tot   = tend_u_ed
    tend_v_tot   = tend_v_ed
    ! tend_t_tot   = tend_thl_tot * ((pfull/pstd_mks)**kappa) + tend_ql_tot * HLv/Cp_air 
    ! The clipping tendency #1 for qt and ql should not enter T tendency
    tend_t_tot   = tend_thl_tot * ((pfull/pstd_mks)**kappa) + (tend_ql_tot - tend_ql_clip1) * HLv/Cp_air 
    tend_q_tot   = tend_qt_tot - tend_ql_tot
    
    precip_ed = - tend_qt_cld
    precip_mf = - tend_qt_ut
    ccov_ed   = ccov * a_env
    
    ccov_mf   = 0.
    do iupd = 1, nupd
        where (ql_upd(:,:,:,iupd) .gt. small)
            ccov_mf(:,:,:) = ccov_mf(:,:,:) + a_upd(:,:,:,iupd)
        end where
    end do
    
    do ilev = 1, nlev
        rain2d_ed(:,:) = rain2d_ed(:,:) + (phalf(:,:,ilev+1) - phalf(:,:,ilev)) * precip_ed(:,:,ilev)/Grav
        rain2d_mf(:,:) = rain2d_mf(:,:) + (phalf(:,:,ilev+1) - phalf(:,:,ilev)) * precip_mf(:,:,ilev)/Grav
    end do
    
end subroutine compute_tend_grid

subroutine compute_tend_cloud (thl_env, qt_env, ql_env, w_env, a_env, u_env, v_env, tke_env, &
                               thl_en3, qt_en3, ql_en3, w_en3, a_en3, u_en3, v_en3, tke_en3, &
                               delta_t, &  ! dt_real is removed
                               pfull, zfull, k_t, ls, dls, rho_half, flux_thl_ed, flux_qt_ed ,&
                               rh, s, sigmas, ccov, wthl, wqt, tend_thl_cld, tend_qt_cld, tend_ql_cld)

    real   , intent(in),    dimension(:,:,:)   :: thl_env, qt_env, ql_env, w_env, a_env, u_env, v_env, tke_env
    real   , intent(out),   dimension(:,:,:)   :: thl_en3, qt_en3, ql_en3, w_en3, a_en3, u_en3, v_en3, tke_en3
    real   , intent(in)                        :: delta_t ! dt_real is removed
    real   , intent(in),    dimension(:,:,:)   :: pfull, zfull, k_t, ls, dls, rho_half, flux_thl_ed, flux_qt_ed
    real   , intent(out),   dimension(:,:,:)   :: rh, s, sigmas, ccov, wthl, wqt
    real   , intent(out),   dimension(:,:,:)   :: tend_thl_cld, tend_qt_cld, tend_ql_cld

    real   , dimension(size(thl_env,1),size(thl_env,2),size(thl_env,3)) :: &
                                        pip, tl_env, qln_env, qsat_env, dqsat_env, &
                                        gamma, a, b, qone, tmpliq, liqnew, qnew, tnew, prec
    real :: tau_q = 600.
    integer :: nlev    
    
    nlev = size(thl_env, 3)    
    
    ! ed_leap_frog removed by ZTAN -- 09/10/2017
    
    pip(:,:,:) = (pfull(:,:,:)/pstd_mks)**kappa
    tl_env(:,:,:) = thl_env(:,:,:)* pip(:,:,:)
    qln_env(:,:,:) = ql_env(:,:,:)
    
    call compute_ql(tl_env, qt_env, qln_env, pfull, qsat_env, dqsat_env)
    rh(:,:,:) = qt_env(:,:,:)/qsat_env(:,:,:)
    gamma(:,:,:) = dqsat_env(:,:,:) * HLv / Cp_air
    s(:,:,:) = qt_env(:,:,:) - qsat_env(:,:,:)
    
    a(:,:,:) = 1./(1.+gamma(:,:,:))
    b(:,:,:) = a(:,:,:)*pip(:,:,:)*dqsat_env(:,:,:)
    
    if (sigma_qonly) then
        a(:,:,:) = 1.; b(:,:,:) = 0.
    elseif (sigma_tonly) then
        a(:,:,:) = 0.; b(:,:,:) = pip(:,:,:)*Cp_air/HLv
    end if
    
    
    call cal_sigmas_new(a_env, rho_half, a, b, k_t, thl_env, qt_env,           &
                        flux_thl_ed, flux_qt_ed, ls, dls, tau_q, tke_env, zfull, & 
                        sigmas, wthl, wqt)
                        
    qone(:,:,:) = s(:,:,:)/sigmas(:,:,:)
    ccov(:,:,:) = max(min(0.5 + 0.36 * atan(1.55 * qone(:,:,:)), 1.), 0.)
    ccov(:,:,1) = 0. ! Top level is 0 cloud cover, but WHY? 08/12/2011
    
    where (qone(:,:,:) .gt. 0.)
        tmpliq(:,:,:) = sigmas(:,:,:) * (ccov(:,:,:)*qone(:,:,:)+exp(-1.-(qone(:,:,:)**2.)/2.))
    elsewhere
        tmpliq(:,:,:) = sigmas(:,:,:) * exp(-1.+1.2*qone(:,:,:))
    end where
    
    tmpliq(:,:,:) = min(max(tmpliq(:,:,:), 0.), qt_env(:,:,:))
    liqnew(:,:,:) = tmpliq(:,:,:); qnew(:,:,:) = qt_env(:,:,:) - liqnew(:,:,:)
    tnew(:,:,:) = thl_env(:,:,:)*pip(:,:,:) + tmpliq(:,:,:) * HLv/Cp_air
    
    call compute_precip_env(qsat_env, liqnew, ccov, delta_t, prec) 
    
    thl_en3(:,:,:) =   (tnew(:,:,:) - liqnew(:,:,:) * HLv/Cp_air)/ pip(:,:,:)
     qt_en3(:,:,:) =    qnew(:,:,:) + liqnew(:,:,:)
     ql_en3(:,:,:) =  liqnew(:,:,:)
      w_en3(:,:,:) =   w_env(:,:,:) 
      a_en3(:,:,:) =   a_env(:,:,:) 
      u_en3(:,:,:) =   u_env(:,:,:) 
      v_en3(:,:,:) =   v_env(:,:,:) 
    tke_en3(:,:,:) = tke_env(:,:,:) 
    
    tend_thl_cld(:,:,:) = ( thl_en3(:,:,:) - thl_env(:,:,:) ) * a_env(:,:,:)/delta_t
    tend_qt_cld (:,:,:) = (  qt_en3(:,:,:) -  qt_env(:,:,:) ) * a_env(:,:,:)/delta_t 
    tend_ql_cld (:,:,:) = (  ql_en3(:,:,:) -  ql_env(:,:,:) ) * a_env(:,:,:)/delta_t
    

end subroutine compute_tend_cloud

    

subroutine compute_precip_env(qsat_env, liqnew, ccov, delt, prec) 

    real   , intent(in),    dimension(:,:,:)   :: qsat_env, ccov
    real   , intent(inout), dimension(:,:,:)   :: liqnew
    real   , intent(in)                        :: delt
    real   , intent(out)  , dimension(:,:,:)   :: prec
    
    real :: ct = 1.0e-4
    real :: cw = 5.0e-4
    
    prec(:,:,:) = 0.
    
    if (precip_env_opt == 1) then      ! Precip threshold with fixed fraction (precip_env_fac)
        prec(:,:,:) = max(0., liqnew(:,:,:) - precip_env_fac * qsat_env(:,:,:) * ccov(:,:,:))
    elseif (precip_env_opt == 2) then  ! Precip threshold with fixed value    (precip_env_val)
        prec(:,:,:) = max(0., liqnew(:,:,:) - precip_env_val * ccov(:,:,:))
    elseif (precip_env_opt == 3) then  ! New precip scheme 
        prec(:,:,:) = - delt * liqnew(:,:,:) *   &
                     (ct * (1. - exp(-(liqnew(:,:,:)/max(ccov(:,:,:), 1.e-10)/cw))) )
        prec(:,:,:) = min(liqnew(:,:,:), prec(:,:,:))
    end if
    
    liqnew(:,:,:) = liqnew(:,:,:) - prec(:,:,:)

end subroutine compute_precip_env



subroutine cal_sigmas_new(a_env, rho_half, a, b, k_t, thl, qt,           &
                          flux_thl_ed, flux_qt_ed, ls, dls, tau_q, tke_env, zfull, & 
                          sigmas, wthl, wqt)

    real   , intent(in),    dimension(:,:,:)   :: thl, qt, a_env, tke_env
    real   , intent(in),    dimension(:,:,:)   :: rho_half, zfull, a, b, k_t, ls, dls
    real   , intent(in),    dimension(:,:,:)   :: flux_thl_ed, flux_qt_ed
    real   , intent(in)                        :: tau_q
    real   , intent(out),   dimension(:,:,:)   :: sigmas, wthl, wqt

    real   , dimension(size(thl,1),size(thl,2),size(thl,3)) :: dqdz, dthldz, ra_half, pre_sigmas
    real   , dimension(size(thl,1),size(thl,2)) :: ha, hb
    integer :: nlev, ilev
    
    real :: small = 1.0e-10

    nlev = size(thl, 3) 
    
    ! --- Interpolation for the gradients of thl and qt --- !
    
    dqdz(:,:,:) = 0.; dthldz(:,:,:) = 0.
    
    ha(:,:) = zfull(:,:,2) - zfull(:,:,1)
    hb(:,:) = zfull(:,:,3) - zfull(:,:,1)
    dqdz  (:,:,1) = (ha * ha * ( qt(:,:,3) -  qt(:,:,1)) - hb * hb * ( qt(:,:,2) -  qt(:,:,1))) / &
                    (ha * hb * (ha - hb))
    dthldz(:,:,1) = (ha * ha * (thl(:,:,3) - thl(:,:,1)) - hb * hb * (thl(:,:,2) - thl(:,:,1))) / &
                    (ha * hb * (ha - hb))
                    
    do ilev = 2, nlev-1
        ha(:,:) = zfull(:,:,ilev-1) - zfull(:,:,ilev)
        hb(:,:) = zfull(:,:,ilev+1) - zfull(:,:,ilev)
        dqdz  (:,:,ilev) = (ha * ha * ( qt(:,:,ilev+1) -  qt(:,:,ilev)) - hb * hb * ( qt(:,:,ilev-1) -  qt(:,:,ilev))) / &
                           (ha * hb * (ha - hb))
        dthldz(:,:,ilev) = (ha * ha * (thl(:,:,ilev+1) - thl(:,:,ilev)) - hb * hb * (thl(:,:,ilev-1) - thl(:,:,ilev))) / &
                           (ha * hb * (ha - hb))    
    end do
    
    ha(:,:) = zfull(:,:,nlev-1) - zfull(:,:,nlev)
    hb(:,:) = zfull(:,:,nlev-2) - zfull(:,:,nlev)
    dqdz  (:,:,nlev) = (ha * ha * ( qt(:,:,nlev-2) -  qt(:,:,nlev)) - hb * hb * ( qt(:,:,nlev-1) -  qt(:,:,nlev))) / &
                       (ha * hb * (ha - hb))
    dthldz(:,:,nlev) = (ha * ha * (thl(:,:,nlev-2) - thl(:,:,nlev)) - hb * hb * (thl(:,:,nlev-1) - thl(:,:,nlev))) / &
                       (ha * hb * (ha - hb))
    
    ! --- Interpolation of wthl and wqt from half levels to full levels --- !
    ra_half(:,:,1) = 0.; 
    ra_half(:,:,2:nlev) = (a_env(:,:,1:nlev-1)+a_env(:,:,2:nlev))/2. * rho_half(:,:,2:nlev)
    
    wthl(:,:,1) = 0.; wqt(:,:,1) = 0.
    wthl(:,:,2:nlev) = flux_thl_ed(:,:,2:nlev)/ra_half(:,:,2:nlev)  ! This is at half level
    wqt (:,:,2:nlev) = flux_qt_ed (:,:,2:nlev)/ra_half(:,:,2:nlev)  ! This is at half level
    
    wthl(:,:,1:nlev-1) = (wthl(:,:,1:nlev-1)+wthl(:,:,2:nlev))/2.   ! Interpolate to full level
    wqt (:,:,1:nlev-1) = (wqt (:,:,1:nlev-1)+wqt (:,:,2:nlev))/2.   ! Interpolate to full level
    
    ! Note the bottom level of interpolation is not that accurate... but
    ! perhaps the probablistic cloud scheme is not that important there.
    
    
    ! --- Actual calculation for Sigma --- !
    
    if (sigma_opt == 1) then     ! SCM: oldsigma == 0
        pre_sigmas(:,:,:) = a*b * (wthl*dqdz + wqt*dthldz) - a*a * (wqt*dqdz) - b*b * (wthl*dthldz)
        pre_sigmas(:,:,:) = (max(pre_sigmas(:,:,:), 0.))**.5
        where (tke_env(:,:,:) .le. 0.)
            sigmas(:,:,:) = tau_q**.5/((2*ce)**.5)*pre_sigmas(:,:,:)
        elsewhere
            sigmas(:,:,:) = dls(:,:,:)**.5/(tke_env(:,:,:)**.25 * (2*ce)**.5)*pre_sigmas(:,:,:)
        end where
        sigmas(:,:,:) = max(sigmas(:,:,:), small)
        
    elseif (sigma_opt == 2) then ! SCM: oldsigma == -1
        pre_sigmas(:,:,:) = a*b * (wthl*dqdz + wqt*dthldz) - a*a * (wqt*dqdz) - b*b * (wthl*dthldz)
        pre_sigmas(:,:,:) = (max(pre_sigmas(:,:,:), 0.))**.5
        sigmas(:,:,:) = tau_q**.5/((2*ce)**.5)*pre_sigmas(:,:,:)
        sigmas(:,:,:) = max(sigmas(:,:,:), small)
        
    elseif (sigma_opt == 3) then ! SCM: oldsigma == 1
        pre_sigmas(:,:,:) = (a*dqdz)**2. + (b*dthldz)**2. - 2.*a*b*dqdz*dthldz
        pre_sigmas(:,:,:) = (max(pre_sigmas(:,:,:), 0.))**.5
        sigmas(:,:,:) = (2.0*tau_q*k_t(:,:,:))**.5 * pre_sigmas(:,:,:)
        sigmas(:,:,:) = max(sigmas(:,:,:), small)
        
    else
        sigmas = small
    end if
    

end subroutine cal_sigmas_new


! ------------------------- !
!      ED TENDENCIES        !
! ------------------------- !

! NOTE: kbot is not used yet...
subroutine compute_tend_ed (is, ie, js, je, Time,                                         & 
                            thl_env, qt_env, ql_env, w_env, a_env, u_env, v_env, tke_env, &
                            pfull, zfull, phalf, dt_real, delta_t,                        &
                            ustar, wstar, zstar, oblength,                                &
                            k_m, k_t, t_surf, flux_t, flux_q, flux_u, flux_v, flux_r,     &
                            dtaudv_atm, dhdt_surf, dedt_surf, dedq_surf, drdt_surf,       &
                            dhdt_atm, dedq_atm, net_surf_sw_down, surf_lw_down, Tri_surf, &
                            ls, dls, diss_heat, rho_half,                                 &
                            thl_en2, qt_en2, ql_en2, w_en2, a_en2, u_en2, v_en2, tke_en2, &
                            flux_thl_ed, flux_qt_ed, flux_ql_ed, flux_u_ed, flux_v_ed,    &
                            tend_thl_ed, tend_qt_ed, tend_ql_ed, tend_u_ed, tend_v_ed,    &
                            edmf_update_ml, istep )

    integer, intent(in)                        :: is, ie, js, je
    type(time_type), intent(in)                :: Time
    real   , intent(in) ,   dimension(:,:,:)   :: thl_env, qt_env, ql_env, w_env, a_env, u_env, v_env, tke_env
    real   , intent(in) ,   dimension(:,:,:)   :: pfull, zfull, phalf
    real   , intent(in) ,   dimension(:,:)     :: ustar, wstar, zstar, oblength
    real   , intent(in)                        :: dt_real, delta_t  ! delta_t for main ED; dt_real ONLY for SST
    
    real   , intent(inout), dimension(:,:,:)   :: k_m, k_t
    real   , intent(inout), dimension(:,:)     :: t_surf, flux_t, flux_q, flux_u, flux_v, &
                                                  dtaudv_atm, dhdt_surf, dedt_surf, dedq_surf, drdt_surf, &
                                                  dhdt_atm, dedq_atm
    real   , intent(in),    dimension(:,:)     :: flux_r, net_surf_sw_down, surf_lw_down
    type(surf_diff_type),   intent(inout)      :: Tri_surf
    logical, intent(in)                        :: edmf_update_ml
    integer, intent(in)                        :: istep
    
    real   , intent(out),   dimension(:,:,:)   :: thl_en2, qt_en2, ql_en2, w_en2, a_en2, u_en2, v_en2, tke_en2
    real   , intent(out),   dimension(:,:,:)   :: flux_thl_ed, flux_qt_ed, flux_ql_ed, flux_u_ed, flux_v_ed
    real   , intent(out),   dimension(:,:,:)   :: tend_thl_ed, tend_qt_ed, tend_ql_ed, tend_u_ed, tend_v_ed
    real   , intent(out),   dimension(:,:,:)   :: ls, dls, diss_heat, rho_half    
    
    real   , dimension(size(thl_env,1),size(thl_env,2),size(thl_env,3)) :: t_env, q_env, & 
                  dt_u, dt_v, dt_thl, dt_qt, dt_ql, nu, nu_m ! thl_enn, qt_enn, ql_enn, u_enn, v_enn are removed
    integer :: nlev
    
    
    nlev = size(thl_env, 3)
    
    call compute_eddy_diffusivity (thl_env, qt_env, ql_env, w_env, a_env, u_env, v_env, tke_env, &
         ustar, wstar, zstar, oblength, k_m, k_t, pfull, zfull, ls, dls)
    
    call cloud_decompose (thl_env, qt_env, ql_env, pfull, t_env, q_env)
    
    dt_u = 0.; dt_v = 0.; dt_thl = 0.; dt_qt = 0.; dt_ql = 0.
    
    ! ed_leap_frog removed by ZTAN -- 09/10/2017
    
    if (surf_explicit) then
        !! downward sweep  <-- Write in vert_diff.f90
        call gcm_vert_diff_down_ed_new (1, 1, delta_t,               &
             u_env, v_env, t_env, q_env, thl_env, qt_env, ql_env, &
             k_m, k_t, a_env, phalf, pfull, zfull,                & 
             flux_u, flux_v, dtaudv_atm * 0.,                     &
             dt_u, dt_v, dt_thl, dt_qt, dt_ql, diss_heat,         &
             Tri_surf, nu_m, nu, rho_half)
    
        !! update surface
        if (edmf_update_ml) then  ! Note: dt_real is used here, because SST is advanced simply forward
            ! ZTAN: 09/12/2017 -- HERE MAY BE WRONG -- May need to use 1, ie-is+1, 1, je-js+1 instead.
            ! call mixed_layer (is, ie, js, je, Time,                            &
            call mixed_layer (1, ie-is+1, 1, je-js+1, Time,                            &
                              t_surf, flux_t, flux_q, flux_r, flux_u, dt_real, &
                              net_surf_sw_down, surf_lw_down, Tri_surf,        &
                              dhdt_surf * 0., dedt_surf * 0., dedq_surf * 0., drdt_surf * 0.,      &
                              dhdt_atm * 0.,  dedq_atm * 0.)
        else
            ! call mixed_layer_noupdate (is, ie, js, je, Time,                   &
            call mixed_layer_noupdate (1, ie-is+1, 1, je-js+1, Time,                   &
                              t_surf, flux_t, flux_q, flux_r, flux_u, dt_real, &
                              net_surf_sw_down, surf_lw_down, Tri_surf,        &
                              dhdt_surf * 0., dedt_surf * 0., dedq_surf * 0., drdt_surf * 0.,      &
                              dhdt_atm * 0.,  dedq_atm * 0.)
        end if
        
        
    else
        
        call gcm_vert_diff_down_ed_new (1, 1, delta_t,            &
             u_env, v_env, t_env, q_env, thl_env, qt_env, ql_env, &
             k_m, k_t, a_env, phalf, pfull, zfull,                & 
             flux_u, flux_v, dtaudv_atm,                          &
             dt_u, dt_v, dt_thl, dt_qt, dt_ql, diss_heat,         &
             Tri_surf, nu_m, nu, rho_half)
         
        if (edmf_update_ml) then  ! Note: dt_real is used here, because SST is advanced simply forward
            ! call mixed_layer (is, ie, js, je, Time,                            &
            call mixed_layer (1, ie-is+1, 1, je-js+1, Time,                            &
                              t_surf, flux_t, flux_q, flux_r, flux_u, dt_real, &
                              net_surf_sw_down, surf_lw_down, Tri_surf,        &
                              dhdt_surf, dedt_surf, dedq_surf, drdt_surf,      &
                              dhdt_atm * a_env(:,:,nlev),  dedq_atm * a_env(:,:,nlev))
        else
            ! call mixed_layer_noupdate (is, ie, js, je, Time,                   &
            call mixed_layer_noupdate (1, ie-is+1, 1, je-js+1, Time,                   &
                              t_surf, flux_t, flux_q, flux_r, flux_u, dt_real, &
                              net_surf_sw_down, surf_lw_down, Tri_surf,        &
                              dhdt_surf, dedt_surf, dedq_surf, drdt_surf,      &
                              dhdt_atm * a_env(:,:,nlev),  dedq_atm * a_env(:,:,nlev))
        end if        
        
    end if

    !! upward sweep
    call gcm_vert_diff_up_ed_new (1, 1, delta_t, Tri_surf, dt_thl, dt_qt, dt_ql)
             
    ! New env values after delta_t (NOT dt_real! CHANGED 09/11/2017)
    thl_en2(:,:,:) = thl_env(:,:,:) + dt_thl(:,:,:) * delta_t ! dt_real
     qt_en2(:,:,:) =  qt_env(:,:,:) + dt_qt (:,:,:) * delta_t ! dt_real 
     ql_en2(:,:,:) =  ql_env(:,:,:) + dt_ql (:,:,:) * delta_t ! dt_real 
      w_en2(:,:,:) =   w_env(:,:,:)
      a_en2(:,:,:) =   a_env(:,:,:) 
      u_en2(:,:,:) =   u_env(:,:,:) + dt_u  (:,:,:) * delta_t ! dt_real
      v_en2(:,:,:) =   v_env(:,:,:) + dt_v  (:,:,:) * delta_t ! dt_real
    tke_en2(:,:,:) = tke_env(:,:,:)
    
    
    ! New env values after delt (the 'implicit' values) <- NOT correct
    ! thl_enn(:,:,:) = thl_env(:,:,:) + dt_thl(:,:,:) * delt
    !  qt_enn(:,:,:) =  qt_env(:,:,:) + dt_qt (:,:,:) * delt 
    !  ql_enn(:,:,:) =  ql_env(:,:,:) + dt_ql (:,:,:) * delt 
    !   u_enn(:,:,:) =   u_env(:,:,:) + dt_u  (:,:,:) * delt
    !   v_enn(:,:,:) =   v_env(:,:,:) + dt_v  (:,:,:) * delt
    
    flux_thl_ed(:,:,1) = 0.; flux_qt_ed (:,:,1) = 0.; flux_ql_ed (:,:,1) = 0.
    flux_u_ed  (:,:,1) = 0.; flux_v_ed  (:,:,1) = 0.
    
    ! flux_thl_ed(:,:,2:nlev) = nu  (:,:,2:nlev) * (thl_enn(:,:,2:nlev) - thl_enn(:,:,1:nlev-1))
    ! flux_qt_ed (:,:,2:nlev) = nu  (:,:,2:nlev) * (qt_enn (:,:,2:nlev) - qt_enn (:,:,1:nlev-1))
    ! flux_ql_ed (:,:,2:nlev) = nu  (:,:,2:nlev) * (ql_enn (:,:,2:nlev) - ql_enn (:,:,1:nlev-1))
    ! flux_u_ed  (:,:,2:nlev) = nu_m(:,:,2:nlev) * (u_enn  (:,:,2:nlev) - u_enn  (:,:,1:nlev-1))
    ! flux_v_ed  (:,:,2:nlev) = nu_m(:,:,2:nlev) * (v_enn  (:,:,2:nlev) - v_enn  (:,:,1:nlev-1))

    flux_thl_ed(:,:,2:nlev) = nu  (:,:,2:nlev) * (thl_en2(:,:,2:nlev) - thl_en2(:,:,1:nlev-1))
    flux_qt_ed (:,:,2:nlev) = nu  (:,:,2:nlev) * (qt_en2 (:,:,2:nlev) - qt_en2 (:,:,1:nlev-1))
    flux_ql_ed (:,:,2:nlev) = nu  (:,:,2:nlev) * (ql_en2 (:,:,2:nlev) - ql_en2 (:,:,1:nlev-1))
    flux_u_ed  (:,:,2:nlev) = nu_m(:,:,2:nlev) * (u_en2  (:,:,2:nlev) - u_en2  (:,:,1:nlev-1))
    flux_v_ed  (:,:,2:nlev) = nu_m(:,:,2:nlev) * (v_en2  (:,:,2:nlev) - v_en2  (:,:,1:nlev-1))

    flux_thl_ed(:,:,nlev+1) = flux_t(:,:) / Cp_air * ((pstd_mks/phalf(:,:,nlev+1))**kappa)
    flux_qt_ed (:,:,nlev+1) = flux_q(:,:)
    flux_ql_ed (:,:,nlev+1) = 0.
    flux_u_ed  (:,:,nlev+1) = flux_u(:,:)
    flux_v_ed  (:,:,nlev+1) = flux_v(:,:)
    ! 12:02
    
    tend_thl_ed(:,:,:) = dt_thl(:,:,:) * a_env(:,:,:) 
    tend_qt_ed (:,:,:) = dt_qt (:,:,:) * a_env(:,:,:) 
    tend_ql_ed (:,:,:) = dt_ql (:,:,:) * a_env(:,:,:) 
    tend_u_ed  (:,:,:) = dt_u  (:,:,:) 
    tend_v_ed  (:,:,:) = dt_v  (:,:,:) 
    
end subroutine compute_tend_ed



subroutine compute_eddy_diffusivity (thl_env, qt_env, ql_env, w_env, a_env, u_env, v_env, tke_env, &
           ustar, wstar, zstar, oblength, k_m, k_t, pfull, zfull, ls, dls)

    real   , intent(in) ,   dimension(:,:,:)   :: thl_env, qt_env, ql_env, w_env, a_env, u_env, v_env, tke_env
    real   , intent(in) ,   dimension(:,:)     :: ustar, wstar, zstar, oblength
    real   , intent(in) ,   dimension(:,:,:)   :: pfull, zfull
    real   , intent(inout), dimension(:,:,:)   :: k_m, k_t
    real   , intent(out),   dimension(:,:,:)   :: ls, dls
    
    real   , dimension(size(thl_env,1),size(thl_env,2),size(thl_env,3)) :: k_tke
    real   , dimension(size(thl_env,1),size(thl_env,2)) :: zstar_mom
    integer :: nlev, ilev
    
    nlev = size(thl_env,3)

    call compute_mixing_length (thl_env, qt_env, ql_env, tke_env, & 
                                wstar, zstar, oblength, pfull, zfull, ls)
    
    k_tke(:,:,:) = c_diff * ls(:,:,:) * (max(0., tke_env(:,:,:))**.5)
    
    if (diff_opt == 1) then      ! k-profile diffusivity <- WRITE THIS !!
            
        zstar_mom(:,:) = zstar(:,:) / pbl_depth_res
        do ilev = 1, nlev
            where (zfull(:,:,ilev) .lt. zstar(:,:))  ! Use rescaled zstar for scalar mixing
                k_t(:,:,ilev) = vonk * ((ustar(:,:)/wstar(:,:))**3.+ &
                                        (39. * vonk * zfull(:,:,ilev)/zstar(:,:)))**(1./3.)
                k_t(:,:,ilev) = k_t(:,:,ilev) * (zfull(:,:,ilev)/zstar(:,:)) * (1.-zfull(:,:,ilev)/zstar(:,:))**2.
                k_t(:,:,ilev) = k_t(:,:,ilev) * zstar(:,:) * wstar(:,:)
            elsewhere
                k_t(:,:,ilev) = 0.0
            end where
            
            where (zfull(:,:,ilev) .lt. zstar_mom(:,:))  ! Use original zstar for momentum mixing
                k_m(:,:,ilev) = vonk * ((ustar(:,:)/wstar(:,:))**3.+ &
                                        (39. * vonk * zfull(:,:,ilev)/zstar_mom(:,:)))**(1./3.)
                k_m(:,:,ilev) = k_m(:,:,ilev) * (zfull(:,:,ilev)/zstar_mom(:,:)) * (1.-zfull(:,:,ilev)/zstar_mom(:,:))**2.
                k_m(:,:,ilev) = k_m(:,:,ilev) * zstar_mom(:,:) * wstar(:,:)
            elsewhere
                k_m(:,:,ilev) = 0.0
            end where
            
            ! Fix diffusivity with TKE-formulation
            if (diff1_fix == 1) then
                where (wstar(:,:) .le. 0.0) ! Stable condition, reduce diffusivity to 10%
                    where (zfull(:,:,ilev) .lt. zstar(:,:)) 
                        k_t(:,:,ilev) = k_t(:,:,ilev) * 0.1
                    end where
                    where (zfull(:,:,ilev) .lt. zstar_mom(:,:)) 
                        k_m(:,:,ilev) = k_m(:,:,ilev) * 0.1
                    end where
                end where
            elseif (diff1_fix == 2) then
                where (wstar(:,:) .le. 0.0) ! Stable condition, use TKE-based diffusivity instead
                    where (zfull(:,:,ilev) .lt. zstar(:,:)) 
                        k_t(:,:,ilev) = k_tke(:,:,ilev)
                    end where
                    where (zfull(:,:,ilev) .lt. zstar_mom(:,:)) 
                        k_m(:,:,ilev) = k_tke(:,:,ilev)
                    end where
                end where
            
            elseif (diff1_fix == 3) then
                where (wstar(:,:) .le. 0.0) ! Stable condition, use TKE-based diffusivity, and reduce to 1% 
                    where (zfull(:,:,ilev) .lt. zstar(:,:)) 
                        k_t(:,:,ilev) = k_tke(:,:,ilev) * .01
                    end where
                    where (zfull(:,:,ilev) .lt. zstar_mom(:,:)) 
                        k_m(:,:,ilev) = k_tke(:,:,ilev) * .01
                    end where
                end where
                
            end if
            
        end do
        
    elseif (diff_opt == 2) then  ! TKE-based diffusivity
        k_m(:,:,:) = k_tke(:,:,:)
        k_t(:,:,:) = k_tke(:,:,:)
    end if
    
    k_m = max(diff_min, k_m) * diff_rescale
    k_t = max(diff_min, k_t) * diff_rescale
    
    dls(:,:,:) = ls(:,:,:) / dls_factor !2.5
           
end subroutine compute_eddy_diffusivity


subroutine compute_mixing_length (thl_env, qt_env, ql_env, tke_env, & 
                                  wstar, zstar, oblength, pfull, zfull, ls)

    real   , intent(in) ,   dimension(:,:,:)   :: thl_env, qt_env, ql_env, tke_env
    real   , intent(in) ,   dimension(:,:)     :: wstar, zstar, oblength
    real   , intent(in) ,   dimension(:,:,:)   :: pfull, zfull
    real   , intent(out),   dimension(:,:,:)   :: ls
    
    real   , dimension(size(thl_env,1),size(thl_env,2),size(thl_env,3)) :: thv_env, rho_env, N
    real   , dimension(size(thl_env,1),size(thl_env,2)) :: tau, L1, L2, L3, L23
    integer :: nlev, ilev
    
    nlev = size(thl_env,3)
    ls = 0.; L1 = 0.; L2 = 0.; L3 = 0.; L23 = 0.
    
    call compute_tau(zstar, wstar, tau)
    call compute_thv(thl_env, qt_env, ql_env, pfull, thv_env, rho_env)
    
    if (ml_opt == 0) then ! mixing length ~ k*z
        ls(:,:,:) = vonk * zfull(:,:,:)
        if (ml_limiter) then
            do ilev = 1, nlev
                ls(:,:,ilev) = ls(:,:,ilev) * exp(-zfull(:,:,ilev)/5./zstar(:,:))
            end do
        end if
        
    elseif (ml_opt == 1) then
        do ilev = 1, nlev
            L1(:,:) = tau(:,:) * (max(0., tke_env(:,:,ilev))**.5)
            if (ml_limiter) then
                L1(:,:) = L1(:,:) * exp(-zfull(:,:,ilev)/5./zstar(:,:))
            end if
            
            where (oblength(:,:) .lt. 0.)
                L2(:,:) = vonk * zfull(:,:,ilev) * ( (1. - 100.*zfull(:,:,ilev)/oblength(:,:))**0.2 )
            end where
            
            where (oblength(:,:) .gt. 0.)
                L2(:,:) = vonk * zfull(:,:,ilev) * ( (1. + 2.7 *zfull(:,:,ilev)/oblength(:,:))**(-1.))
            end where
            
            where (oblength(:,:) == 0.)  
            ! Neutral <- wthv = 0.0, Actual oblength = infty. (see compute_zstar_wstar_obl subroutine)
                L2(:,:) = vonk * zfull(:,:,ilev)
            end where
            
            where (L1(:,:) .gt. 0.0)
                ls(:,:,ilev) = 1./(1./L1(:,:) + 1./L2(:,:) + 1./lmax)
            end where
        end do
        
        
    elseif (ml_opt == 2) then
        N(:,:,:) = 0.0;
        N(:,:,1:nlev-1) = (grav/300. * max((thv_env(:,:,2:nlev) - thv_env(:,:,1:nlev-1)) / &
                                             (zfull(:,:,2:nlev) -   zfull(:,:,1:nlev-1)), 0. ))**0.5
        N(:,:,2:nlev) = N(:,:,2:nlev) + N(:,:,1:nlev-1)
        N(:,:,2:nlev-1) = N(:,:,2:nlev-1)/2.  ! This is runmean(N,2), different from SCM (runmean(N,3))
        
        do ilev = 1, nlev
            L1(:,:) = vonk * zfull(:,:,ilev)
            L2(:,:) = tau(:,:) * (max(0., tke_env(:,:,ilev))**.5)
            if (ml_limiter) then
                L2(:,:) = L2(:,:) * exp(-zfull(:,:,ilev)/5./zstar(:,:))
            end if
            
            where (N(:,:,ilev) .gt. 0.0)
                L3(:,:) = 0.7 * (max(0., tke_env(:,:,ilev))**.5)/N(:,:,ilev)
                ! L23(:,:) = 1./(1./L2(:,:) + 1./max(L3(:,:), (zhalf(:,:,ilev) - zhalf(:,:,ilev+1))/2.)) ! Suselj's SCM
                L23(:,:) = 1./(1./L2(:,:) + 1./L3(:,:))  ! Formulation in Suselj et al. [2012]
            elsewhere
                L23(:,:) = L2(:,:)
            end where
            
            ls(:,:,ilev) = L23(:,:) + (L1(:,:) - L23(:,:)) * exp(-zfull(:,:,ilev)/0.1/zstar(:,:))
        end do
        
    elseif (ml_opt == 3) then
        do ilev = 1, nlev
            ls(:,:,ilev) = tau(:,:) * (max(0., tke_env(:,:,ilev))**.5) + &
                (vonk * zfull(:,:,ilev) - tau(:,:) * & 
                (max(0., tke_env(:,:,ilev))**.5)) * exp(-zfull(:,:,ilev)/100.)
        end do
    end if
    
end subroutine compute_mixing_length


subroutine compute_tau (zstar, wstar, tau)

    real   , intent(in) ,   dimension(:,:)     :: wstar, zstar
    real   , intent(out),   dimension(:,:)     :: tau
    
    real :: small = 1.0e-10
    
    tau(:,:) = small
    
    if (tau_opt == 1) then
        tau(:,:) = 0.5*zstar(:,:)/max(wstar(:,:), small) 
        ! Different from SCM: tau = infty if wstar <= 0 (SCM gives tau = 0.0)
    elseif (tau_opt == 2) then
        tau(:,:) = 600.0
    elseif (tau_opt == 3) then
        tau(:,:) = 400.0
    end if

end subroutine compute_tau

! ------------------------- !
!      MF TENDENCIES        !
! ------------------------- !

subroutine update_surf_flux_mf (t_in, q_in, t_st1, q_st1, dhdt_atm, dedq_atm, flux_t, flux_q)
    real  , intent(in), dimension(:,:,:)    :: t_in, q_in, t_st1, q_st1   
    real  , intent(in), dimension(:,:)      :: dhdt_atm, dedq_atm
    real  , intent(inout), dimension(:,:)   :: flux_t, flux_q

    integer :: nlev
    nlev = size(t_in,3)
    
    flux_t(:,:) = flux_t(:,:) + dhdt_atm(:,:) * (t_st1(:,:,nlev) - t_in (:,:,nlev))
    flux_q(:,:) = flux_q(:,:) + dedq_atm(:,:) * (q_st1(:,:,nlev) - q_in (:,:,nlev))

end subroutine update_surf_flux_mf


subroutine mf_env_implicit (thl_env, qt_env, ql_env, pfull, zhalf, rhot, dt_real, &
                            thl_en1, qt_en1, ql_en1, a_en1,              &
                            thl_st1, qt_st1, ql_st1, t_st1, q_st1, thv_st1, rho_st1, &
                            tend_thl_mf, tend_qt_mf, tend_ql_mf, flux_thl_mf, flux_qt_mf, flux_ql_mf, net_mf_tot)
            
    real   , intent(in) ,   dimension(:,:,:)   :: thl_env,  qt_env,  ql_env, pfull, zhalf, rhot
    real   , intent(inout), dimension(:,:,:)   :: thl_en1,  qt_en1,  ql_en1
    real   , intent(in) ,   dimension(:,:,:)   :: a_en1
    real   , intent(inout), dimension(:,:,:)   :: thl_st1, qt_st1, ql_st1, &
                                                  t_st1, q_st1, thv_st1, rho_st1
    real   , intent(in)                        :: dt_real
    real   , intent(in)   , dimension(:,:,:)   :: net_mf_tot
    real   , intent(inout), dimension(:,:,:)   :: tend_thl_mf, tend_qt_mf, tend_ql_mf, &
                                                  flux_thl_mf, flux_qt_mf, flux_ql_mf
                                                  
    real   , dimension(size(thl_env,1),size(thl_env,2),size(thl_env,3)) :: dthl0, dqt0, dql0, dthl1, dqt1, dql1
    real   , dimension(size(thl_env,1),size(thl_env,2)) :: dz, c2, ct
    
    integer :: nlev, ilev
    nlev = size(thl_env,3)
    
    dthl0(:,:,:) = thl_en1(:,:,:) - thl_env(:,:,:)
    dqt0(:,:,:) = qt_en1(:,:,:) - qt_env(:,:,:)
    dql0(:,:,:) = ql_en1(:,:,:) - ql_env(:,:,:)
    
    ilev = 1
    dz(:,:) = zhalf(:,:,ilev+1) - zhalf(:,:,ilev)  ! DOUBLE-CHECK, because it is opposite to the other formulation !
    c2(:,:) = 0.; ct(:,:) = 1. - dt_real * net_mf_tot(:,:,ilev+1)/(rhot(:,:,ilev)*a_en1(:,:,ilev)*dz(:,:))
    dthl1(:,:,ilev) = dthl0(:,:,ilev) * (1-ct(:,:))/ct(:,:)
    dqt1 (:,:,ilev) = dqt0 (:,:,ilev) * (1-ct(:,:))/ct(:,:)
    dql1 (:,:,ilev) = dql0 (:,:,ilev) * (1-ct(:,:))/ct(:,:)
    
    do ilev = 2, nlev
        dz(:,:) = zhalf(:,:,ilev+1) - zhalf(:,:,ilev)  ! DOUBLE-CHECK, because it is opposite to the other formulation !
        c2(:,:) =    - dt_real * net_mf_tot(:,:,ilev)  /(rhot(:,:,ilev)*a_en1(:,:,ilev)*dz(:,:))
        ct(:,:) = 1. - dt_real * net_mf_tot(:,:,ilev+1)/(rhot(:,:,ilev)*a_en1(:,:,ilev)*dz(:,:))
        dthl1(:,:,ilev) = (dthl0(:,:,ilev) * (1-ct(:,:)) + (dthl1(:,:,ilev-1) + dthl0(:,:,ilev-1)) * c2(:,:) )/ct(:,:)
        dqt1 (:,:,ilev) = (dqt0 (:,:,ilev) * (1-ct(:,:)) + (dqt1 (:,:,ilev-1) + dqt0 (:,:,ilev-1)) * c2(:,:) )/ct(:,:)
        dql1 (:,:,ilev) = (dql0 (:,:,ilev) * (1-ct(:,:)) + (dql1 (:,:,ilev-1) + dql0 (:,:,ilev-1)) * c2(:,:) )/ct(:,:)
    end do
    
    thl_en1(:,:,:) = thl_en1(:,:,:) + dthl1(:,:,:)
    qt_en1 (:,:,:) = qt_en1 (:,:,:) + dqt1 (:,:,:)
    ql_en1 (:,:,:) = ql_en1 (:,:,:) + dql1 (:,:,:)
    
    thl_st1(:,:,:) = thl_st1(:,:,:) + dthl1(:,:,:) * a_en1(:,:,:)
    qt_st1 (:,:,:) = qt_st1 (:,:,:) + dqt1 (:,:,:) * a_en1(:,:,:)
    ql_st1 (:,:,:) = ql_st1 (:,:,:) + dql1 (:,:,:) * a_en1(:,:,:)
    
    call cloud_decompose(thl_st1, qt_st1, ql_st1, pfull, t_st1, q_st1)
    call compute_thv(thl_st1, qt_st1, ql_st1, pfull, thv_st1, rho_st1)
    
    tend_thl_mf(:,:,:) = tend_thl_mf(:,:,:) + dthl1(:,:,:) * a_en1(:,:,:) / dt_real
    tend_qt_mf (:,:,:) = tend_qt_mf (:,:,:) + dqt1 (:,:,:) * a_en1(:,:,:) / dt_real
    tend_ql_mf (:,:,:) = tend_ql_mf (:,:,:) + dql1 (:,:,:) * a_en1(:,:,:) / dt_real
    
    flux_thl_mf(:,:,2:nlev) = flux_thl_mf(:,:,2:nlev) - net_mf_tot(:,:,2:nlev) * (thl_en1(:,:,1:nlev-1) - thl_env(:,:,1:nlev-1))
    flux_qt_mf (:,:,2:nlev) = flux_qt_mf (:,:,2:nlev) - net_mf_tot(:,:,2:nlev) * (qt_en1 (:,:,1:nlev-1) - qt_env (:,:,1:nlev-1))
    flux_ql_mf (:,:,2:nlev) = flux_ql_mf (:,:,2:nlev) - net_mf_tot(:,:,2:nlev) * (ql_en1 (:,:,1:nlev-1) - ql_env (:,:,1:nlev-1))

end subroutine mf_env_implicit


subroutine update_grid_mean_mf (u_in,  v_in,  w_in,  tke_in,  thl_in,  qt_in,  ql_in, pfull, dt_real, &
                                tend_thl_ut, tend_qt_ut, tend_ql_ut, tend_thl_mf, tend_qt_mf, tend_ql_mf, &
                                u_st1, v_st1, w_st1, tke_st1, thl_st1, qt_st1, ql_st1, &
                                t_st1, q_st1, thv_st1, rho_st1)
        
    real   , intent(in) ,   dimension(:,:,:)   :: u_in,  v_in,  w_in,  tke_in,  thl_in,  qt_in,  ql_in,  pfull
    real   , intent(in)                        :: dt_real
    real   , intent(in) ,   dimension(:,:,:)   :: tend_thl_ut, tend_qt_ut, tend_ql_ut, tend_thl_mf, tend_qt_mf, tend_ql_mf
    real   , intent(out),   dimension(:,:,:)   :: u_st1, v_st1, w_st1, tke_st1, thl_st1, qt_st1, ql_st1      
    real   , intent(out),   dimension(:,:,:)   :: t_st1, q_st1, thv_st1, rho_st1
    
    u_st1   (:,:,:) = u_in   (:,:,:)
    v_st1   (:,:,:) = v_in   (:,:,:)
    w_st1   (:,:,:) = w_in   (:,:,:)
    tke_st1 (:,:,:) = tke_in (:,:,:)
    
    thl_st1 (:,:,:) = thl_in (:,:,:) + (tend_thl_ut(:,:,:) + tend_thl_mf(:,:,:)) * dt_real
    qt_st1  (:,:,:) = qt_in  (:,:,:) + (tend_qt_ut (:,:,:) + tend_qt_mf (:,:,:)) * dt_real 
    ql_st1  (:,:,:) = ql_in  (:,:,:) + (tend_ql_ut (:,:,:) + tend_ql_mf (:,:,:)) * dt_real
    
    call cloud_decompose(thl_st1, qt_st1, ql_st1, pfull, t_st1, q_st1)
    call compute_thv(thl_st1, qt_st1, ql_st1, pfull, thv_st1, rho_st1)

end subroutine update_grid_mean_mf



subroutine compute_tend_mf (thl_upd,  qt_upd,  ql_upd,  w_upd,  a_upd, &
                            thl_env,  qt_env,  ql_env,  w_env,  rhot, wt, zhalf, &
                            tend_thl_mf, tend_qt_mf, tend_ql_mf, flux_thl_mf, flux_qt_mf, flux_ql_mf, net_mf_tot)
                            
    real   , intent(in) ,   dimension(:,:,:,:) :: thl_upd,  qt_upd,  ql_upd,  w_upd,  a_upd                 
    real   , intent(in) ,   dimension(:,:,:)   :: thl_env,  qt_env,  ql_env,  w_env,  rhot, wt, zhalf
    real   , intent(out),   dimension(:,:,:)   :: tend_thl_mf, tend_qt_mf, tend_ql_mf, &
                                                  flux_thl_mf, flux_qt_mf, flux_ql_mf, net_mf_tot

    real   , dimension(size(thl_upd,1),size(thl_upd,2),size(thl_upd,3),size(thl_upd,4)) :: net_mf
    real   , dimension(size(thl_upd,1),size(thl_upd,2)) :: dz
    integer :: nlev, nupd, ilev, iupd
    
    nlev = size(thl_upd,3)
    nupd = size(thl_upd,4) 


    net_mf_tot (:,:,:) = 0.
    do iupd = 1, nupd
        ! net_mf(:,:,:,iupd) = rhot(:,:,:) * a_upd(:,:,:,iupd) * ( w_upd(:,:,:,iupd) - w_env(:,:,:)) ! ZTAN 2017/09/11: seems wrong!
        net_mf(:,:,:,iupd) = rhot(:,:,:) * a_upd(:,:,:,iupd) * ( w_upd(:,:,:,iupd) - wt(:,:,:))
        net_mf_tot (:,:,1:nlev) = net_mf_tot (:,:,1:nlev) + net_mf(:,:,:,iupd)
    end do
    
    flux_thl_mf = 0.; flux_qt_mf = 0.; flux_ql_mf = 0.
    
    if (mf_full_uw) then
        do iupd = 1, nupd
            flux_thl_mf (:,:,1:nlev) = flux_thl_mf (:,:,1:nlev) + net_mf(:,:,1:nlev,iupd) * thl_upd(:,:,1:nlev,iupd)
            flux_qt_mf  (:,:,1:nlev) = flux_qt_mf  (:,:,1:nlev) + net_mf(:,:,1:nlev,iupd) * qt_upd (:,:,1:nlev,iupd)
            flux_ql_mf  (:,:,1:nlev) = flux_ql_mf  (:,:,1:nlev) + net_mf(:,:,1:nlev,iupd) * ql_upd (:,:,1:nlev,iupd)
        end do
        flux_thl_mf (:,:,2:nlev) = flux_thl_mf (:,:,2:nlev) - net_mf_tot(:,:,2:nlev) * thl_env(:,:,1:nlev-1)
        flux_qt_mf  (:,:,2:nlev) = flux_qt_mf  (:,:,2:nlev) - net_mf_tot(:,:,2:nlev) * qt_env (:,:,1:nlev-1)
        flux_ql_mf  (:,:,2:nlev) = flux_ql_mf  (:,:,2:nlev) - net_mf_tot(:,:,2:nlev) * ql_env (:,:,1:nlev-1)
    else
        do iupd = 1, nupd
            flux_thl_mf (:,:,:) = flux_thl_mf (:,:,:) + net_mf(:,:,:,iupd) * (thl_upd(:,:,:,iupd) - thl_env(:,:,:))
            flux_qt_mf  (:,:,:) = flux_qt_mf  (:,:,:) + net_mf(:,:,:,iupd) * (qt_upd (:,:,:,iupd) - qt_env (:,:,:))
            flux_ql_mf  (:,:,:) = flux_ql_mf  (:,:,:) + net_mf(:,:,:,iupd) * (ql_upd (:,:,:,iupd) - ql_env (:,:,:))
        end do
    end if
    
    do ilev = nlev, 1, -1  ! (Does not work if flux_phi_mf(:,:,1) != 0)
        dz(:,:) = zhalf(:,:,ilev) - zhalf(:,:,ilev+1)
        tend_thl_mf(:,:,ilev) = - (flux_thl_mf(:,:,ilev) - flux_thl_mf(:,:,ilev+1)) / rhot(:,:,ilev) / dz(:,:) 
        tend_qt_mf (:,:,ilev) = - (flux_qt_mf (:,:,ilev) - flux_qt_mf (:,:,ilev+1)) / rhot(:,:,ilev) / dz(:,:)
        tend_ql_mf (:,:,ilev) = - (flux_ql_mf (:,:,ilev) - flux_ql_mf (:,:,ilev+1)) / rhot(:,:,ilev) / dz(:,:)
    end do
    
end subroutine compute_tend_mf




! ------------------------- !
!      Updraft model LBL    !
! ------------------------- !
! compute_updr_lbl: main subroutine for level-by-level updraft computation
subroutine compute_updr_lbl ( thl_upd,  qt_upd,  ql_upd,  w_upd,  a_upd,  & ! updr (in)
                              thl_unew, qt_unew, ql_unew, w_unew, a_unew, & ! updr (out), Added 09/11/2017 by ZTAN
                              thl_env,  qt_env,  ql_env,  w_env,  a_env,  u_env,  v_env,  tke_env,  & ! envr
                              thl_grid, qt_grid, ql_grid, w_grid,         u_grid, v_grid, tke_grid, & ! grid-mean
                              rhot, zstar, ustar, wthl_surf, wq_surf, wstar, oblength, &  ! oblength added 09/25/2017
                              pfull, phalf, zfull, zhalf, dt,     & 
                              thl_ut, qt_ut, ql_ut, w_ut, a_ut, & ! updr_out_tot (bulk updraft values)
                              m_ut, buoy_ut, chic_ut, entL_ut, entB_ut, detL_ut, detB_ut, & ! updr_out_tot (bulk updraft values)
                              entL_u, entB_u, detL_u, detB_u, entr_u, detr_u, buoy_u, &  ! individual updraft values
                              tend_thl_ut, tend_qt_ut, tend_ql_ut, istep)   ! For debug
                             
    real   , intent(in)                        :: dt
    real   , intent(in) ,   dimension(:,:,:,:) :: thl_upd,  qt_upd,  ql_upd,  w_upd,  a_upd 
    real   , intent(out) ,  dimension(:,:,:,:) :: thl_unew, qt_unew, ql_unew, w_unew, a_unew                
    real   , intent(out) ,  dimension(:,:,:)   :: thl_env,  qt_env,  ql_env,  w_env,  a_env,  u_env,  v_env,  tke_env               
    real   , intent(in) ,   dimension(:,:,:)   :: thl_grid, qt_grid, ql_grid, w_grid,         u_grid, v_grid, tke_grid
    real   , intent(in) ,   dimension(:,:,:)   :: rhot, pfull, phalf, zfull, zhalf     
    real   , intent(in) ,   dimension(:,:)     :: zstar, ustar, wthl_surf, wq_surf, wstar, oblength ! oblength added 09/25/2017    
    real   , intent(out),   dimension(:,:,:)   :: thl_ut, qt_ut, ql_ut, w_ut, a_ut, &
                                                  m_ut, buoy_ut, chic_ut, entL_ut, entB_ut, detL_ut, detB_ut, &
                                                  tend_thl_ut, tend_qt_ut, tend_ql_ut   
    real   , intent(out),   dimension(:,:,:,:) :: entL_u, entB_u, detL_u, detB_u, entr_u, detr_u, buoy_u  
                                                  ! entr_u and detr_u moved to output; buoy_u added 10/10/2017
    integer, intent(in)                        :: istep
    
    
    real   , dimension(size(thl_upd,1),size(thl_upd,2),size(thl_upd,3),size(thl_upd,4)) :: &  ! updr_old
                                                  thl_uin,  qt_uin,  ql_uin,  w_uin,  a_uin, &
                                                  thv_uin, rho_uin, chic_uin, buoy_uin, & ! chic_u is not actually used
                                                  tend_thl, tend_qt, tend_ql
                                          
    real   , dimension(size(thl_upd,1),size(thl_upd,2),size(thl_upd,3),size(thl_upd,4)) :: &  
                                                  thv_unew, rho_unew, m_unew
    real   , dimension(size(thl_upd,1),size(thl_upd,2),size(thl_upd,3)) :: thv_env, rho_env, thv_grid
            
    real   , dimension(size(thl_upd,1),size(thl_upd,2)) :: thl_std, qt_std, dz,      &
                                                           tl_2d, qsat_2d, dqsat_2d, &
                                                           thv_uinsum, thv_unewsum,  &
                                                           zero_2d, adj_fac2d, c3_lim
    real   , dimension(size(thl_upd,1),size(thl_upd,2),size(thl_upd,4)) :: thl_anom, qt_anom, c1, c2, c3, d4

    ! The general algorithm is: 
    ! 1. Compute updraft buoyancy, using updraft values from the LAST step;
    ! 2. Compute entrainment and detrainment rates, using updraft buoyancy and velocity from the LAST step;
    ! 3. Update updraft scalar values by combining three mass sources:
    !    i.   Updraft scalar values from the LAST step;
    !    ii.  Incoming mass flux and updraft scalar values from the lower level of the CURRENT step; 
    !    iii. Environmental scalar values from the LAST step, and entrainment rate computed in (2).
    ! 4. Update updraft buoyancy, using updraft scalar values from the CURRENT step; 
    ! 5. Compute updraft velocity and mass flux, using updraft buoyancy from the CURRENT step.
    
    ! This algorithm helps prevent the on-off oscillation of updraft and the resulting oscillation
    ! of updraft buoyancy.
    
    ! c1: mass of the last step divided by dt (= rho * a_u(old) * dz/dt)
    ! c2: mass flux from the lower level into the current level (= rho(low) * a_u(low) * w_u(low)).
    ! c3: entrained mass per unit time ('explicit' with mass flux of last step)
    !     (= rho * a_u(old) * w_u(old) * entr * dz)
    ! d4: detrained mass per unit time per unit mass flux (= rho * detr * dz);
    !     the actual detrained mass per unit time is 'implicit' with mass flux of this step (= d4*a_u(new)*w_u(new)).
    ! Units for c1, c2, c3, d4*a_u*w_u are all kg/m2/s.
    ! The explicit treatment for entrainment and the implicit treatment for detrainment help stabilize the scheme.
    
                                                  
    integer :: nlev, nupd, iupd, ilev
    
    real :: alpha_b  = 1./3.
    real :: alpha_d  = 0.075
    real :: radius_d = 1500.
    
    nlev = size(thl_upd,3)
    nupd = size(thl_upd,4)
    zero_2d(:,:) = 0.
    chic_uin(:,:,:,:) = 0.
    
    ! 0. Make a local copy of updraft variables
    thl_uin = thl_upd; qt_uin = qt_upd; ql_uin = ql_upd; w_uin = w_upd; a_uin = a_upd
    
    ! 1. Do the surface level
    !    * au is fixed at the surface, and phiu is prescribed: 
    !    * Thus: compute surface std; update surface phiu; compute entT and wu.
    
    ! 1.1. Adjust thl_uin, qt_uin, and ql_uin at the surface to be the Gaussian tail value (ZTAN 09/25/2017)
    !    * Now only adjust the lowest level and not any higher levels.
    call compute_surface_std(wthl_surf, ustar, zfull(:,:,nlev), oblength, thl_std)
    call compute_surface_std(wq_surf,   ustar, zfull(:,:,nlev), oblength, qt_std)  
    
    do iupd = 1, nupd 
       thl_anom(:,:,     iupd)  = thl_std(:,:) * init_scl_u(iupd) * gaussian_std * var_surf_fac
       qt_anom (:,:,     iupd)  = qt_std(:,:)  * init_scl_u(iupd) * gaussian_std * var_surf_fac
       thl_uin (:,:,nlev,iupd)  = thl_grid(:,:,nlev) + thl_anom(:,:, iupd)
       qt_uin  (:,:,nlev,iupd)  = qt_grid (:,:,nlev) +  qt_anom(:,:, iupd)
       
       tl_2d = thl_uin (:,:,nlev,iupd)/((pstd_mks/pfull(:,:,nlev))**kappa)
       call compute_ql_2d(tl_2d, qt_uin   (:,:,nlev,iupd), ql_uin   (:,:,nlev,iupd), & 
                          pfull(:,:,nlev), qsat_2d, dqsat_2d)
    end do
    
    ! 1.2. Compute the environmental values for ALL levels.
    call envupd_decompose(thl_grid, qt_grid, ql_grid, w_grid,         u_grid, v_grid, tke_grid, &
                          thl_uin,  qt_uin,  ql_uin,  w_uin,  a_uin,                            & 
                          thl_env,  qt_env,  ql_env,  w_env,  a_env,  u_env,  v_env,  tke_env)
                          
    ! 1.3. Compute thl_u, qt_u, ql_u at the lowest level
    dz(:,:) = zhalf(:,:,nlev) - zhalf(:,:,nlev+1)
    
    ! 1.3.1. entr_u, detr_u (turbulent ent/det rates) are simply set to zero at the lowest level.
    entr_u(:,:,nlev,:) = 0.; detr_u(:,:,nlev,:) = 0.
    
    do iupd = 1, nupd
       ! Trivial combination for prescribed thl_u, qt_u, ql_u 
       ! --> the obtained new updraft values are simply: thl_unew = thl_uin, qt_unew = qt_uin, ql_unew = ql_uin
       c1(:,:,iupd) = rhot(:,:,nlev)*a_uin(:,:,nlev,iupd)*dz(:,:)*updraft_rescale/dt
       c2(:,:,iupd) = 0.
       c3(:,:,iupd) = rhot(:,:,nlev)*a_uin (:,:,nlev,nupd)*w_uin (:,:,nlev,nupd)*entr_u(:,:,nlev,nupd)*dz(:,:) ! = 0
       d4(:,:,iupd) = rhot(:,:,nlev)*detr_u(:,:,nlev,nupd)*dz(:,:) ! = 0
       
       call compute_phiu_2d ( thl_uin (:,:,nlev,iupd), qt_uin (:,:,nlev,iupd), ql_uin (:,:,nlev,iupd), &  ! old values
                              thl_uin (:,:,nlev,iupd), qt_uin (:,:,nlev,iupd), ql_uin (:,:,nlev,iupd), &  ! lower values 
                              thl_uin (:,:,nlev,iupd), qt_uin (:,:,nlev,iupd), ql_uin (:,:,nlev,iupd), &  ! entrained values
                              c1(:,:,iupd), c2(:,:,iupd), c3(:,:,iupd), rhot(:,:,nlev), pfull (:,:,nlev), dz(:,:), dt, &
                              thl_unew(:,:,nlev,iupd), qt_unew(:,:,nlev,iupd), ql_unew(:,:,nlev,iupd), &
                              tend_thl(:,:,nlev,iupd), tend_qt(:,:,nlev,iupd), tend_ql(:,:,nlev,iupd))
    end do
    
    ! 1.4. Compute buoy_u at the lowest level, using a_u at the OLD step 
    !      (= a_u at the NEW step, since updraft fraction is fixed at the lowest level)
    thv_unewsum(:,:) = 0.
    do iupd = 1, nupd 
       call compute_thv_2d(thl_unew(:,:,nlev,iupd),  qt_unew(:,:,nlev,iupd), &
                            ql_unew(:,:,nlev,iupd),    pfull(:,:,nlev),      & 
                           thv_unew(:,:,nlev,iupd), rho_unew(:,:,nlev,iupd))
                           
       thv_unewsum(:,:) = thv_unewsum(:,:) + a_uin(:,:,nlev,iupd) * thv_unew(:,:,nlev,iupd) * updraft_rescale
    end do
    call compute_thv_2d(thl_env(:,:,nlev),  qt_env(:,:,nlev), &
                         ql_env(:,:,nlev),   pfull(:,:,nlev), & 
                        thv_env(:,:,nlev), rho_env(:,:,nlev))
    thv_grid(:,:,nlev) = thv_env(:,:,nlev) * a_env(:,:,nlev) + thv_unewsum(:,:)
    
    do iupd = 1, nupd
        buoy_u(:,:,nlev,iupd) = Grav * (thv_unew(:,:,nlev,iupd)/thv_grid(:,:,nlev)-1.)
    end do
    buoy_uin(:,:,nlev,:) = buoy_u(:,:,nlev,:)
    
    ! 1.5. Compute a_unew, w_unew, entB_u, detB_u, entL_u, detL_u at the lowest level
    do iupd = 1, nupd
       call compute_au_wu_lowest( a_uin(:,:,nlev,iupd), w_uin(:,:,nlev,iupd), &
                                  zero_2d(:,:), w_env(:,:,nlev), buoy_u(:,:,nlev,iupd), &
                                  c1(:,:,iupd), c2(:,:,iupd), c3(:,:,iupd), d4(:,:,iupd), &
                                  rhot (:,:,nlev), pfull (:,:,nlev), dz(:,:), dt, &
                                  alpha_b, alpha_d, radius_d, &
                                  a_unew(:,:,nlev,iupd), w_unew(:,:,nlev,iupd), m_unew(:,:,nlev,iupd), &
                                  entB_u(:,:,nlev,iupd), detB_u(:,:,nlev,iupd), &
                                  entL_u(:,:,nlev,iupd), detL_u(:,:,nlev,iupd))
                                  
       ! 1.5.1. Rescale tend_thl, tend_qt, tend_ql to account for the boundary massive entrainment
       !        (The actual 'c3' is entB_u, so the tendencies need to be adjusted by a factor of:
       !          (c1+entB_u)/(c1) -- Note that c2 and c3 are both zero)
       adj_fac2d(:,:) = 1. + entB_u(:,:,nlev,iupd)/(c1(:,:,iupd) + c2(:,:,iupd) + c3(:,:,iupd))
       tend_thl(:,:,nlev,iupd) = tend_thl(:,:,nlev,iupd) * adj_fac2d(:,:)
       tend_qt (:,:,nlev,iupd) = tend_qt (:,:,nlev,iupd) * adj_fac2d(:,:)
       tend_ql (:,:,nlev,iupd) = tend_ql (:,:,nlev,iupd) * adj_fac2d(:,:)
    end do
    
    ! 2. Upward loop for all non-surface levels
    do ilev = nlev-1, 1, -1
       dz(:,:) = zhalf(:,:,ilev) - zhalf(:,:,ilev+1)
       
       ! 2.1. Compute buoy_u using OLD thl_u, qt_u, ql_u, a_u (Only used for buoy_u-based entrainment!)
       call compute_thv_2d(thl_env(:,:,ilev),  qt_env(:,:,ilev), &
                            ql_env(:,:,ilev),   pfull(:,:,ilev), & 
                           thv_env(:,:,ilev), rho_env(:,:,ilev))
 
       thv_uinsum(:,:) = 0.
       do iupd = 1, nupd 
          call compute_thv_2d(thl_uin(:,:,ilev,iupd),  qt_uin(:,:,ilev,iupd), &
                               ql_uin(:,:,ilev,iupd),   pfull(:,:,ilev),      & 
                              thv_uin(:,:,ilev,iupd), rho_uin(:,:,ilev,iupd))
                           
          thv_uinsum(:,:) = thv_uinsum(:,:) + a_uin(:,:,ilev,iupd) * thv_uin(:,:,ilev,iupd) * updraft_rescale
       end do
       thv_grid(:,:,ilev) = thv_env(:,:,ilev) * a_env(:,:,ilev) + thv_uinsum(:,:)  ! This will be overwritten later.
    
       do iupd = 1, nupd
          buoy_uin(:,:,ilev,iupd) = Grav * (thv_uin(:,:,ilev,iupd)/thv_grid(:,:,ilev)-1.)
       end do      
       
       do iupd = 1, nupd
          ! 2.2. Use buoy_o, w_o to calculate entr_u and detr_u
          call compute_ent_det_2d (w_uin(:,:,ilev,iupd), buoy_uin(:,:,ilev,iupd), zstar(:,:), &
                                   zfull(:,:,ilev), entr_u(:,:,ilev,iupd), detr_u(:,:,ilev,iupd))  
          
          ! 2.2.1. Calculate c1, c2, c3, d4
          c1(:,:,iupd) = rhot(:,:,ilev)  *a_uin (:,:,ilev,   iupd)*dz(:,:)*updraft_rescale/dt
          c2(:,:,iupd) = rhot(:,:,ilev+1)*a_unew(:,:,ilev+1, iupd)*w_unew(:,:,ilev+1, iupd)
          c3(:,:,iupd) = rhot(:,:,ilev)  *a_uin (:,:,ilev,   iupd)*w_uin (:,:,ilev,   iupd)*entr_u(:,:,ilev,iupd)*dz(:,:)
          d4(:,:,iupd) = rhot(:,:,ilev)  *detr_u(:,:,ilev,   iupd)*dz(:,:)
          
                     
          ! 2.2.2. Entrainment limiter (for c3) -- entrained mass cannot exceed the environmental mass.
          call compute_c3_lim( a_uin(:,:,ilev, iupd), w_uin (:,:,ilev, iupd), thv_uin (:,:,ilev, iupd), &
                               w_unew(:,:,ilev+1, iupd), thv_unew(:,:,ilev+1, iupd), &
                               w_env (:,:,ilev), a_env(:,:,ilev), a_unew(:,:,nlev, iupd), & ! a_unew(nlev) is surface area fraction
                               thv_grid(:,:,ilev), &
                               c1(:,:,iupd), c2(:,:,iupd), d4(:,:,iupd), rhot(:,:,ilev), pfull (:,:,ilev), dz(:,:), dt, &
                               alpha_b, alpha_d, radius_d, &
                               c3_lim(:,:))
          c3(:,:,iupd) = min(c3(:,:,iupd), c3_lim(:,:)) 
          
          ! 2.3. Calculate thl_u, qt_u, ql_u at Level ilev
          call compute_phiu_2d (   &
                   thl_uin (:,:,ilev,  iupd), qt_uin (:,:,ilev,  iupd), ql_uin (:,:,ilev,  iupd), &  ! old values
                   thl_unew(:,:,ilev+1,iupd), qt_unew(:,:,ilev+1,iupd), ql_unew(:,:,ilev+1,iupd), &  ! lower values (new) 
                   thl_env (:,:,ilev),        qt_env (:,:,ilev),        ql_env (:,:,ilev),        &  ! entrained values
                   c1(:,:,iupd), c2(:,:,iupd), c3(:,:,iupd), rhot(:,:,ilev), pfull (:,:,ilev), dz(:,:), dt, &
                   thl_unew(:,:,ilev,  iupd), qt_unew(:,:,ilev,  iupd), ql_unew(:,:,ilev,  iupd), &
                   tend_thl(:,:,ilev,  iupd), tend_qt(:,:,ilev,  iupd), tend_ql(:,:,ilev,  iupd)  )
       end do
       
       ! 2.3. Recompute buoy_u at Level ilev, using thl_u/qt_u/ql_u at the NEW step (but a_u at the OLD step)
       thv_unewsum(:,:) = 0.
       do iupd = 1, nupd 
          call compute_thv_2d(thl_unew(:,:,ilev,iupd),  qt_unew(:,:,ilev,iupd), &
                               ql_unew(:,:,ilev,iupd),    pfull(:,:,ilev),      & 
                              thv_unew(:,:,ilev,iupd), rho_unew(:,:,ilev,iupd))
                           
          thv_unewsum(:,:) = thv_unewsum(:,:) + a_uin(:,:,ilev,iupd) * thv_unew(:,:,ilev,iupd) * updraft_rescale
       end do
       thv_grid(:,:,ilev) = thv_env(:,:,ilev) * a_env(:,:,ilev) + thv_unewsum(:,:)
    
       do iupd = 1, nupd
          buoy_u(:,:,ilev,iupd) = Grav * (thv_unew(:,:,ilev,iupd)/thv_grid(:,:,ilev)-1.)
       end do       

       ! 2.4. Compute a_unew, w_unew, entB_u (= 0), detB_u, entL_u, detL_u at Level ilev
       do iupd = 1, nupd
          call compute_au_wu_2d ( a_uin(:,:,ilev,iupd), w_uin(:,:,ilev,iupd), &
                                  w_unew(:,:,ilev+1, iupd), w_env(:,:,ilev), buoy_u(:,:,ilev,iupd), &
                                  a_unew(:,:,nlev, iupd), & ! a_unew(nlev) is surface area fraction for au-limiter
                                  c1(:,:,iupd), c2(:,:,iupd), c3(:,:,iupd), d4(:,:,iupd), &
                                  rhot (:,:,ilev), pfull (:,:,ilev), dz(:,:), dt, &
                                  alpha_b, alpha_d, radius_d, &
                                  a_unew(:,:,ilev,iupd), w_unew(:,:,ilev,iupd), m_unew(:,:,ilev,iupd), &
                                  entB_u(:,:,ilev,iupd), detB_u(:,:,ilev,iupd), &
                                  entL_u(:,:,ilev,iupd), detL_u(:,:,ilev,iupd))
          
       ! 2.4.1. Since boundary entrainment (entB_u) = 0, no need to rescale tend_thl, tend_qt, tend_ql 
       end do
    end do
   
    call compute_bulk_updraft(thl_unew, qt_unew, ql_unew, a_unew, tend_thl, tend_qt, tend_ql, &
                              buoy_u, chic_uin, m_unew, entL_u, entB_u, detL_u, detB_u, & 
                              thl_env, qt_env, ql_env, rhot, zhalf, phalf, nlev, nupd, &
                              thl_ut, qt_ut, ql_ut, w_ut, a_ut, & 
                              m_ut, buoy_ut, chic_ut, entL_ut, entB_ut, detL_ut, detB_ut, &
                              tend_thl_ut, tend_qt_ut, tend_ql_ut)
    ! DEBUG
    ! write(*,*) istep                          
end subroutine compute_updr_lbl
  
! compute_wu_lowest: compute w_u at the lowest level, assuming a_u and buoy_u are given
subroutine compute_au_wu_lowest (a_o, w1, w2, w3, buoy, &
                                 c1, c2, c3, d4, rhot, pfull, dz, dt, &
                                 alpha_b, alpha_d, radius_d, &
                                 a, w, m, entB, detB, entL, detL )
                                 
    real   , intent(in),    dimension(:,:) :: a_o, w1, & ! old values
                                                   w2, & ! lower values 
                                                   w3, & ! entrained values
                                                 buoy, & ! estimated buoyancy
                                              c1, c2, c3, d4, rhot, pfull, dz
    real   , intent(in)                    :: dt
    real   , intent(in)                    :: alpha_b, alpha_d, radius_d
    real   , intent(out),   dimension(:,:) :: a, w, m, entB, detB, entL, detL
    
    real   , dimension(size(dz,1),size(dz,2))  :: f_1, f_2, coef_2, coef_1, coef_0, db, beta_d
    real :: beta_b
    
    beta_b = 1 - alpha_b
    beta_d(:,:) = alpha_d/radius_d/sqrt(max(a_o(:,:), 1.e-3))
    db(:,:) = rhot(:,:)*dz(:,:)
    
	! Fix SFC fraction by:
	! (1) If buoy > 0 => w > 0, aN < a0 because entL = 0 => entB > 0.
	! (2) If buoy <= 0. Swap airmass, use detB = c1, calculate entB accordingly.

	! Calculate entB for target a(ilev, iupd) = a_o(ilev, iupd)
	! final w-equation: (c1 + c2 + c3 + {entB}) * {w} = ...
	!                   (c1 * w1) + (c2 * w2) + (c3 + {entB}) * w3 + ...
	!                   (db * beta_b * a) * buoy - db * a * beta_d * ({w} - w3)^2
	! final a-equation: (c1 + c2 + c3 - d4*a*{w} + {entB} - {detB}) = ...
	!                    rhot*a *dz*updraft_rescale/dt + rhot*a*{w}
	!                  -> f_2 * {w} + f_1 = {entB}
	! Thus     (c1 + c2 + c3 + f_2 * {w} + f_1) * {w} = ...
	!                   (c1 * w1) + (c2 * w2) + (c3 + f_2 * {w} + f_1) * w3 + ...
	!                   (db * beta_b * a) * buoy - db * a * beta_d * ({w} - w3)^2
	! ** Now entL and detL are set to zero at SFC.
	
	detB(:,:) = 0.0
		
	a(:,:) = a_o(:,:)
	f_1(:,:) = rhot(:,:)*a(:,:)*dz(:,:)*updraft_rescale/dt - &
		       (c1(:,:)+c2(:,:)+c3(:,:) - detB(:,:))
	f_2(:,:) = (rhot(:,:) + d4(:,:))*a(:,:)

	coef_2(:,:) = f_2(:,:) + db(:,:)*a(:,:)* beta_d(:,:)    ! coef_2 > 0 by definition
	coef_1(:,:) = c1(:,:)+c2(:,:)+c3(:,:) + f_1(:,:) - w3(:,:)*f_2(:,:) &
				  - 2. * db(:,:) * beta_d(:,:) * a(:,:) * w3(:,:)
	coef_0(:,:) = -(c1(:,:)*w1(:,:) + c2(:,:)*w2(:,:) + c3(:,:)*w3(:,:) + f_1(:,:)*w3(:,:) &
					+ db(:,:) * beta_b * a(:,:)*buoy(:,:) &
					- db(:,:) * beta_d(:,:) * a(:,:) * w3(:,:) * w3(:,:) )
				   
	where (coef_0(:,:) >= 0.)  ! Updraft at SFC is not buoyant
		w(:,:) = 0.0
		entB(:,:) = 0.
		
		! Note 11/26/2017: The old definition below doubly-counts... 
		!      The entrainment/detrainment at the lowest level is done at step 1.1 of the 
		!      compute_updr_lbl subroutine, and the additional entrainment is only needed
		!      to account for the mass flux of the outgoing updraft.
		! entB(:,:) = rhot(:,:)*a(:,:)*dz(:,:)*updraft_rescale/dt - c2(:,:) - c3(:,:)
		! detB(:,:) = c1(:,:)
		
	elsewhere
		w (:,:) = (-coef_1(:,:) + &
			sqrt(coef_1(:,:)**2. - 4.*coef_2(:,:)*coef_0(:,:)))/(2.*coef_2(:,:))
		entB(:,:) = f_1(:,:) + f_2(:,:) * w(:,:)
	end where ! (coef_0 >= 0.)
    
    entL(:,:) = c3(:,:) 
    detL(:,:) = d4(:,:) * a(:,:) * w(:,:)
    m(:,:) = rhot(:,:)*a(:,:)*w(:,:)
    
end subroutine compute_au_wu_lowest


! compute_au_wu: compute updated updraft fraction au and updraft velocity wu, with other auxiliary variables
! Note: based on the solver 'compute_au_wu_expentr_newp': Use explicit entrainment (but still implicit detrainment)
subroutine compute_au_wu_2d (a_o, w1, w2, w3, buoy, a_s, &
                             c1, c2, c3, d4, rhot, pfull, dz, dt, &
                             alpha_b, alpha_d, radius_d, &
                             a, w, m, entB, detB, entL, detL)
                                 
    real   , intent(in),    dimension(:,:) :: a_o, w1, & ! old values
                                                   w2, & ! lower values 
                                                   w3, & ! entrained values
                                                 buoy, & ! estimated buoyancy
                                              a_s,     & ! surface area fraction 
                                              c1, c2, c3, d4, rhot, pfull, dz
    real   , intent(in)                    :: dt
    real   , intent(in)                    :: alpha_b, alpha_d, radius_d
    real   , intent(out),   dimension(:,:) :: a, w, m, entB, detB, entL, detL

    real   , dimension(size(dz,1),size(dz,2))  :: f_1, f_2, coef_2, coef_1, coef_0, db, beta_d, &
                                                  w_adv, a_denom, c1s, c1d, atd
    integer, dimension(size(dz,1),size(dz,2))  :: adjust_flag
    
    real :: beta_b
    
    beta_b = 1 - alpha_b
    beta_d(:,:) = alpha_d/radius_d/sqrt(max(a_o(:,:), 1.e-3))
    db(:,:) = rhot(:,:)*dz(:,:)
    
    ! ------ Step 1 ------
    ! Assume entB = 0 and detB = 0, and calculate {w},{a},{entL},{detL}.
    ! w-equation: (c1 + c2 + c3) * {w} = ...
    !             (c1 * w1) + (c2 * w2) + (c3 * w3) + ...
    !             (db * beta_b) *{a}* buoy - db * {a} * beta_d * ({w} - w3)^2
    ! a-equation: (c1 + c2 + c3 - d4*{a}*{w}) = ...
    !              rhot*{a} *dz*updraft_rescale/dt + rhot*{a}*{w}
    ! Solve {w} first, then {a}, then compute entL and detL if {a} satisfies condition.
    ! Note: {a} = (c1+c2+c3) / [(d4+rhot)*{w} + rhot*dz*updraft_rescale/dt]
    ! Thus: [{w} - w_adv] * [(d4+rhot)*{w} + rhot*dz*updraft_rescale/dt]
    !           = (db * beta_b) * buoy - db *beta_d * ({w} - w3)^2 
    ! where w_adv =  = [(c1 * w1) + (c2 * w2) + (c3 * w3)]/(c1+c2+c3)

    where (c1(:,:) + c2(:,:) + c3(:,:) <= 0.)  !   no prev updraft and no incoming updraft
        ! ==> no current updraft and no massive entrainment. a_o should be zero here.
        w(:,:) = 0.; a(:,:) = a_o(:,:)
        entL(:,:) = 0.; entB(:,:) = 0.
        detL(:,:) = 0.; detB(:,:) = 0.
    elsewhere
        w_adv (:,:) = (c1(:,:)*w1(:,:)+c2(:,:)*w2(:,:)+c3(:,:)*w3(:,:))/(c1(:,:)+c2(:,:)+c3(:,:))
        coef_2(:,:) = d4(:,:) + rhot(:,:) + db(:,:) * beta_d(:,:) ! > 0 because d4 >= 0
        coef_1(:,:) =   rhot(:,:)*dz(:,:)*updraft_rescale/dt &
                      - (rhot(:,:)+d4(:,:))* w_adv(:,:) & 
                      - 2. * db(:,:) * beta_d(:,:) * w3(:,:) 
        coef_0(:,:) = - w_adv(:,:) * rhot(:,:)*dz(:,:)*updraft_rescale/dt & 
                      - db(:,:) * beta_b * buoy(:,:) &
                      + db(:,:) * beta_d(:,:) * w3(:,:) * w3(:,:)
                              
        where (coef_0(:,:) >= 0) 
        ! Both w_adv and buoy <= 0. Thus coef_1 >= 0, 
        ! Note that coef_2 > 0. Thus, both roots are non-positive.
        ! Updraft stops, boundary detrainment occurs (Step 3).
            w   (:,:) = 0.; a(:,:) = a_o(:,:) ! a will be adjusted by Step 3
            entL(:,:) = 0.; detL(:,:) = 0.
            entB(:,:) = 0.; detB(:,:) = 0.
        elsewhere
            w (:,:) = (-coef_1(:,:) + &
                  sqrt(coef_1(:,:)**2. - 4.*coef_2(:,:)*coef_0(:,:)))/(2.*coef_2(:,:))
                    
            a_denom (:,:) =  rhot(:,:)*dz(:,:)*updraft_rescale/dt + &
                            (rhot(:,:)+d4(:,:))*w(:,:)
            ! This denominator is always positive.
                    
            a(:,:) = (c1(:,:) + c2(:,:) + c3(:,:))/a_denom (:,:)
                    
            entB(:,:) = 0.; detB(:,:) = 0.
            entL(:,:) = c3(:,:) 
            detL(:,:) = d4(:,:) * a(:,:) * w(:,:)
        end where ! (coef_0(:,:) >= 0) 
    end where ! (c1(:,:) + c2(:,:) <= 0.)

    ! ------ Step 3 ------
    ! Adjust area fraction by tuning up detB.
    ! Note: this algorithm only adjusts down (since adjusting up requires modification of entB).
                
    ! Limitation on w_u
    adjust_flag(:,:) = 0
    if (au_optB_wu == 1) then ! detraining parcels with zero velocity
        where ( (w(:,:) <= 1.0e-6) .or. (a(:,:) <= 1.0e-6) )
            adjust_flag(:,:) = 2  ! Kill the updraft
        end where
    elseif (au_optB_wu == 11) then ! detraining parcels with zero velocity or zero area
        where ( (a(:,:) <= 1.0e-6))
            adjust_flag(:,:) = 2  ! Kill the updraft
        elsewhere
            where (w(:,:) <= 1.0e-6)
                adjust_flag(:,:) = 3  ! Kill the updraft gradually
            end where
        end where
    endif ! (au_optB_wu == 1) (au_optB_wu == 11)
            
    ! Limitation on a_u is top priority -- WILL OVERWRITE adjust_flag due to w_u.
    if (au_optB == 1) then    ! constrain updraft fraction < SFC fraction
        where (a(:,:)> a_s(:,:))
            a(:,:) = a_s(:,:);  adjust_flag(:,:) = 1
        end where
    elseif (au_optB == 2) then    ! constrain updraft fraction < 2*SFC fraction
        where (a(:,:)> 2.0*a_s(:,:))
            a(:,:) = 2.0*a_s(:,:);  adjust_flag(:,:) = 1
        end where
    ! au_optB == 100: constrain updraft fraction < 0.5/nupd, not coded yet
    endif ! (au_optB == 1, 2)
                
    where (adjust_flag(:,:) == 3) ! Adjust w to zero directly, 
                                  ! but detrain a over a 15 minute timescale.
        w   (:,:) = 0. 
        a   (:,:) = a_o(:,:) * exp(-dt/900.0)   
        entL(:,:) = c3(:,:); detL(:,:) = 0.
        entB(:,:) = 0. 
        detB(:,:) = c1(:,:)+c2(:,:)+c3(:,:) - &
                    rhot(:,:)*a(:,:)*dz(:,:)*updraft_rescale/dt
    end where ! (adjust_flag(:,:) == 3) 
                
    where (adjust_flag(:,:) == 2) ! Adjust to zero directly.
        w   (:,:) = 0.  
        a   (:,:) = 0.  
        entL(:,:) = c3(:,:); detL(:,:) = 0.
        entB(:,:) = 0.
        detB(:,:) = c1(:,:)+c2(:,:)+c3(:,:)
    end where ! (adjust_flag(:,:) == 2) 
                
    ! c1d is only used where adjust_flag == 1 (but should be prevailed by flag 2 or 3)
    if (au_optB_wu == 1) then   ! Detrainment everything. atd is target a.
        c1d(:,:) = c1(:,:)+c2(:,:)+c3(:,:)
        atd(:,:) = 0.
    elseif (au_optB_wu == 11) then ! Detrainment to target a (= atd).
        c1d(:,:) =   c1(:,:)+c2(:,:)+c3(:,:) &
                   - rhot(:,:)*a(:,:)*dz(:,:)*updraft_rescale/dt
        atd(:,:) = a(:,:)
    endif       
                                
    where (adjust_flag(:,:) == 1)
        ! Calculate detB for target a = a_o, with entB = 0
        ! final w-equation: (c1 + c2 + c3) * {w} = ...
        !                   (c1 * w1) + (c2 * w2) + c3 * w3  + ...
        !                   (db * beta_b) *a* buoy - db * a * beta_d * ({w} - w3)^2
        ! final a-equation: (c1 + c2 + c3 - d4*a{w} - {detB}) = ...
        !                    rhot*a *dz*updraft_rescale/dt + rhot*a*{w}
                    
                    
        coef_2(:,:) = db(:,:) * beta_d(:,:) *a(:,:)  ! coef_2 > 0 by definition
        coef_1(:,:) = c1(:,:)+c2(:,:) + c3(:,:) &
                      - 2.* db(:,:) * beta_d(:,:) *a(:,:) * w3(:,:)
        coef_0(:,:) = -(c1(:,:)*w1(:,:) + c2(:,:)*w2(:,:) + c3(:,:)*w3(:,:) + &
                        + db(:,:) * beta_b * a(:,:)*buoy(:,:) &
                        - db(:,:) * beta_d(:,:)* a(:,:) * w3(:,:) * w3(:,:))
                                                 
        where (coef_0(:,:) >= 0.) ! Updraft stops . <- Note coef_1 > 0.
                                  ! Should be prevailed by adjust_flag = 2 or 3.
            w   (:,:) = 0.
            a   (:,:) = atd(:,:)
            entL(:,:) = c3(:,:); detL(:,:) = 0.
            entB(:,:) = 0.
            detB(:,:) = c1d(:,:)
        elsewhere
            where (coef_2(:,:) == 0.) 
                ! d3 = 0, no entrainment, degenerate equation.
                ! This implies that c1+c2 > 0, thus coef_1 > 0, coef_0 < 0 ==> w > 0.
                w(:,:) = - coef_0(:,:) / coef_1(:,:) 
            elsewhere  
                ! Now coef_2 > 0 and coef_0 < 0 ==> w > 0
                w(:,:) = (-coef_1(:,:) + &
                    sqrt(coef_1(:,:)**2. - 4.*coef_2(:,:)*coef_0(:,:)))/(2.*coef_2(:,:))
            end where 
                        
            entL(:,:) = c3(:,:)
            detL(:,:) = d4(:,:) * a(:,:) * w(:,:)
            entB(:,:) = 0.
            detB(:,:) = c1(:,:)+c2(:,:)+c3(:,:)  - d4(:,:) * a(:,:) * w(:,:) &
                        - rhot(:,:)*a(:,:)*dz(:,:)*updraft_rescale/dt &
                        - rhot(:,:)*a(:,:) * w(:,:)
        end where ! (coef_0(:,:) >= 0) 
    end where ! (adjust_flag(:,:) == 1) 
             
    m(:,:) = rhot(:,:)*a(:,:)*w(:,:)
    
end subroutine compute_au_wu_2d



! compute_phiu_2d: compute updated updraft scalar values phiu, 
!                  given mixing ratio of old mass, lower mass, and entrained mass
subroutine compute_phiu_2d (thl_o, qt_o, ql_o, thl_l, qt_l, ql_l, thl_e, qt_e, ql_e, &
                            c1, c2, c3, rhot, pfull, dz, dt, &
                            thl, qt, ql, s_thl, s_qt, s_ql)

    real   , intent(in),    dimension(:,:) :: thl_o, qt_o, ql_o, & ! old values
                                              thl_l, qt_l, ql_l, & ! lower values 
                                              thl_e, qt_e, ql_e, & ! entrained values
                                              c1, c2, c3, rhot, pfull, dz
    real   , intent(in)                    :: dt
    real   , intent(out),   dimension(:,:) :: thl, qt, ql, s_thl, s_qt, s_ql

    real   , dimension(size(dz,1),size(dz,2))  :: d_thl, d_qt, d_ql, c1r, c2r, c3r, ct
    
    thl = 0.; qt = 0.; ql = 0.; s_thl = 0.; s_qt = 0.; s_ql = 0.; d_thl = 0.; d_qt = 0.; d_ql = 0.
            
    ct (:,:) = c1(:,:) + c2(:,:) + c3(:,:)
    where (ct(:,:) > 0.0)
        c1r(:,:) = c1(:,:)/ct(:,:); c2r(:,:) = c2(:,:)/ct(:,:); c3r(:,:) = c3(:,:)/ct(:,:)
    elsewhere
        c1r(:,:) = 0.; c2r(:,:) = 0.; c3r(:,:) = 1.  ! Fill in environment values for now
    end where
    
    thl(:,:) = c1r(:,:)*thl_o(:,:) + c2r(:,:)*thl_l(:,:) + c3r(:,:)*thl_e(:,:)
    qt (:,:) = c1r(:,:)*qt_o (:,:) + c2r(:,:)*qt_l (:,:) + c3r(:,:)*qt_e (:,:)
    ql (:,:) = c1r(:,:)*ql_o (:,:) + c2r(:,:)*ql_l (:,:) + c3r(:,:)*ql_e (:,:)
            
    ! Call saturation adjustment in updraft  <- 'where' does not work... Redundant calculation including environment?
    if (maxval(ct(:,:)) > 0.0) then 
        call compute_precip_upd(thl(:,:), qt(:,:), ql(:,:), pfull(:,:), dt, &
                                d_thl(:,:), d_qt(:,:), d_ql(:,:)) ! d_thl, d_qt, d_ql are in K, kg/kg, kg/kg
    else
        d_thl = 0.; d_qt = 0.; d_ql = 0.
    end if
                
    ! ct is in kg/m^3 * m/s; thus s_thl is in K/s, and s_qt and s_ql are in kg/kg/s.          
    where (ct(:,:) > 0.0)   
        s_thl(:,:) = ct(:,:)/dz(:,:)/rhot(:,:) * d_thl(:,:) ! fill in d_thl, d_qt, d_ql
        s_qt (:,:) = ct(:,:)/dz(:,:)/rhot(:,:) * d_qt (:,:)
        s_ql (:,:) = ct(:,:)/dz(:,:)/rhot(:,:) * d_ql (:,:)
                
        thl(:,:) = thl(:,:)+d_thl(:,:)
        qt (:,:) = qt (:,:)+d_qt (:,:)
        ql (:,:) = ql (:,:)+d_ql (:,:)
                
    elsewhere
        s_thl(:,:) = 0.0
        s_qt (:,:) = 0.0
        s_ql (:,:) = 0.0
                
        thl(:,:) = thl_e(:,:) 
        qt (:,:) = qt_e (:,:)
        ql (:,:) = ql_e (:,:)
    end where

end subroutine compute_phiu_2d                       
    

subroutine compute_c3_lim (a_o, w1, thv_o, w2, thv_l, w3, a_e, a_s, thv_g, &
                           c1, c2, d4, rhot, pfull, dz, dt, &
                           alpha_b, alpha_d, radius_d, &
                           c3_lim)
                           
    real   , intent(in),    dimension(:,:) :: a_o, w1, thv_o, & ! old values
                                              w2, thv_l, & ! lower values
                                              w3, a_e,   & ! entrained values
                                              a_s, thv_g, & ! env, surf, and grid values
                                              c1, c2, d4, rhot, pfull, dz
    real   , intent(in)                    :: dt
    real   , intent(in)                    :: alpha_b, alpha_d, radius_d
    real   , intent(out),   dimension(:,:) :: c3_lim
                                              
    real , dimension(size(dz,1),size(dz,2)) :: b_nent, ct_nent, au_lim, wu_lim, &
                                               f_1, f_2, coef_2, coef_1, coef_0, db, beta_d
    real :: beta_b
    
    beta_b = 1 - alpha_b
    beta_d(:,:) = alpha_d/radius_d/sqrt(max(a_o(:,:), 1.e-3))
    db(:,:) = rhot(:,:)*dz(:,:)
    
    ! Limit #1: entrained mass < environment mass
    c3_lim(:,:) = rhot(:,:)*a_e(:,:)*dz(:,:)/dt
    
    ! Limit #2: entrained mass should keep a_u limited below the max allowable fraction
    if (au_optB == 1) then       ! constrain updraft fraction < SFC fraction
        au_lim(:,:) = a_s(:,:)
    elseif (au_optB == 2) then   ! constrain updraft fraction < 2*SFC fraction
        au_lim(:,:) = 2.0*a_s(:,:)
    else ! no constraint
        au_lim(:,:) = 1.
    endif
    
    b_nent(:,:) = 0.; wu_lim(:,:) = 0.
    f_1(:,:) = 0.; f_2(:,:) = 0.; coef_2(:,:) = 0.; coef_1(:,:) = 0.; coef_0(:,:) = 0.
    ct_nent(:,:) = c1(:,:) + c2(:,:)
    
    where ( ct_nent(:,:) .gt. 0.) ! otherwise, c1=c2=c3=0 => no updraft at this level
    
       ! ROUGH ESTIMATE for thv, NON-ENTRAINING and NOT DOING SATURATION ADJUSTMENT.
       b_nent(:,:) = (c1(:,:)*thv_o(:,:) + c2(:,:)*thv_l(:,:))/ct_nent(:,:)
       b_nent(:,:) = Grav * (b_nent(:,:)/thv_g(:,:)-1.) ! Convert to buoyancy
       
       ! Calculate c3_lim for target a(ilev, iupd) = au_lim(ilev, iupd)
	   ! final w-equation: (c1 + c2 + {c3_lim}) * {w} = ...
	   !                   (c1 * w1) + (c2 * w2) + {c3_lim} * w3 + ...
	   !                   (db * beta_b * a_lim) * buoy - db * a_lim * beta_d * ({w} - w3)^2
	   ! final a-equation: (c1 + c2 - d4*a*{w} + {c3_lim}) = ...
	   !                    rhot*a_lim *dz*updraft_rescale/dt + rhot*a*{w}
	   !                  -> f_2 * {w} + f_1 = {c3_lim}
	   ! Thus     (c1 + c2 + f_2 * {w} + f_1) * {w} = ...
	   !                   (c1 * w1) + (c2 * w2) + (f_2 * {w} + f_1) * w3 + ...
	   !                   (db * beta_b * a_lim) * buoy - db * a_lim * beta_d * ({w} - w3)^2
       
       f_1(:,:) = rhot(:,:)*au_lim(:,:)*dz(:,:)*updraft_rescale/dt - (c1(:,:)+c2(:,:))
	   f_2(:,:) = (rhot(:,:) + d4(:,:))*au_lim(:,:)

	   coef_2(:,:) = f_2(:,:) + db(:,:)*au_lim(:,:)* beta_d(:,:)    ! coef_2 > 0 by definition
	   coef_1(:,:) = c1(:,:)+c2(:,:) + f_1(:,:) - w3(:,:)*f_2(:,:) &
			         - 2. * db(:,:) * beta_d(:,:) * au_lim(:,:) * w3(:,:)
	   coef_0(:,:) = -(c1(:,:)*w1(:,:) + c2(:,:)*w2(:,:) + f_1(:,:)*w3(:,:) &
					 + db(:,:) * beta_b * au_lim(:,:)*b_nent(:,:) &
					 - db(:,:) * beta_d(:,:) * au_lim(:,:) * w3(:,:) * w3(:,:) )
					       
       where (coef_0(:,:) >= 0.)  ! The non-entraining updraft is already non-buoyant <- Don't entrain!
		   wu_lim(:,:) = 0.0
		   c3_lim(:,:) = 0.
	   elsewhere  ! {c3_lim} now gives an estimate for max allowable entrainment 
	              ! -- IT IS A ROUGH ESTIMATE because buoy is estimated.
		   wu_lim (:,:) = (-coef_1(:,:) + &
			   sqrt(coef_1(:,:)**2. - 4.*coef_2(:,:)*coef_0(:,:)))/(2.*coef_2(:,:))
		   c3_lim(:,:) = max(min(c3_lim(:,:), f_1(:,:) + f_2(:,:) * wu_lim(:,:)), 0.)
	   end where           
    end where   
           
end subroutine compute_c3_lim


subroutine compute_ent_det_2d (w_upd, buoy, zstar, zfull, entr, detr) 

    real   , intent(in) ,   dimension(:,:) :: w_upd, buoy
    real   , intent(in) ,   dimension(:,:) :: zstar
    real   , intent(in) ,   dimension(:,:) :: zfull
    real   , intent(out),   dimension(:,:) :: entr, detr
    
    real , dimension(size(w_upd,1),size(w_upd,2)) :: zt 
    real :: c1  = 0.4
    real :: a1  = 0.15
    real :: a3  = 2.0
    real :: tau = 600.0
    
    ! used in b/wu^2 entrainment/detrainment only
    real :: c_eps = 0.12 ! 0.33 previously
    real :: c_del = 0.12 
    real :: del_b = 0.004
    
    entr = 0.
    detr = 0.
    
    if (ER == 0) then  ! Simple scaling: k/z
        entr(:,:) = c1 / max(zfull(:,:), ER0_zmin)
    elseif (ER == 1) then
        where (zfull(:,:) .ge. zstar(:,:))
            entr(:,:) = 0.0
        elsewhere
            entr(:,:) = c1 * (1./zfull(:,:) + 1./(zstar(:,:) - zfull(:,:)) )
        end where
    ! ER = 1.33 and ER = 2 and ER = 3 in SCM is not coded
    
    elseif (ER == 4) then
        entr(:,:) = c_eps*max(buoy(:,:),0.0)/(max(w_upd(:,:),3.0e-2)**2.0)
    elseif (ER == 5) then
        where (buoy(:,:) >= 0.0)
            entr(:,:) = 1.0/c_tau_ER5/max(w_upd(:,:),3.0e-2)
        end where
    end if
    entr = ER_frac * entr
    
    if (DR == 0) then  ! proportional to entrainment rate
        detr(:,:) = entr(:,:)
    ! DR = 3 in SCM is not coded

    elseif (DR == 4) then
        ! ZTAN 09/26/2017: staggering is eliminated; lateral detrainment formulation is disabled for the lowest level
        where (zfull(:,:) .le. zstar(:,:))
            detr(:,:) = 0.0
        elsewhere
            detr(:,:) = del_b + &
                        c_del*max(-buoy(:,:),0.0)/(max(w_upd(:,:),3.0e-2)**2.0)
        end where
    elseif (DR == 5) then
        where (buoy(:,:) <= 0.0)
            detr(:,:) = 1.0/c_tau_ER5/max(w_upd(:,:),3.0e-2)
        end where
    ! DR = 10 (KF scheme) is not coded
    
    end if
    detr = DR_frac * detr
    
end subroutine compute_ent_det_2d

! ------------------------- !
!      Updraft model OLD    !
! ------------------------- !

! compute_updr: main subroutine for updraft computation
subroutine compute_updr( thl_upd,  qt_upd,  ql_upd,  w_upd,  a_upd,  & ! updr (in)
                         thl_unew, qt_unew, ql_unew, w_unew, a_unew, & ! updr (out), Added 09/11/2017 by ZTAN
                         thl_env,  qt_env,  ql_env,  w_env,  a_env,  u_env,  v_env,  tke_env,  & ! envr
                         thl_grid, qt_grid, ql_grid, w_grid,         u_grid, v_grid, tke_grid, & ! grid-mean
                         rhot, zstar, ustar, wthl_surf, wq_surf, wstar, oblength, &  ! oblength added 09/25/2017
                         pfull, phalf, zfull, zhalf, dt_real,     & 
                         thl_ut, qt_ut, ql_ut, w_ut, a_ut, & ! updr_out_tot (bulk updraft values)
                         m_ut, buoy_ut, chic_ut, entL_ut, entB_ut, detL_ut, detB_ut, & ! updr_out_tot (bulk updraft values)
                         entL_u, entB_u, detL_u, detB_u, entr_u, detr_u, buoy_u, &  ! individual updraft values
                         tend_thl_ut, tend_qt_ut, tend_ql_ut, istep)   ! For debug
                        
    ! Updraft is not evolved with the main dynamics, thus dt_real should be
    ! used instead of delta_t (for leapfrog)
                             
    real   , intent(in)                        :: dt_real
    real   , intent(in) ,   dimension(:,:,:,:) :: thl_upd,  qt_upd,  ql_upd,  w_upd,  a_upd 
    real   , intent(out) ,  dimension(:,:,:,:) :: thl_unew, qt_unew, ql_unew, w_unew, a_unew                
    real   , intent(out) ,  dimension(:,:,:)   :: thl_env,  qt_env,  ql_env,  w_env,  a_env,  u_env,  v_env,  tke_env               
    real   , intent(in) ,   dimension(:,:,:)   :: thl_grid, qt_grid, ql_grid, w_grid,         u_grid, v_grid, tke_grid
    real   , intent(in) ,   dimension(:,:,:)   :: rhot, pfull, phalf, zfull, zhalf     
    real   , intent(in) ,   dimension(:,:)     :: zstar, ustar, wthl_surf, wq_surf, wstar, oblength ! oblength added 09/25/2017    
    real   , intent(out),   dimension(:,:,:)   :: thl_ut, qt_ut, ql_ut, w_ut, a_ut, &
                                                  m_ut, buoy_ut, chic_ut, entL_ut, entB_ut, detL_ut, detB_ut, &
                                                  tend_thl_ut, tend_qt_ut, tend_ql_ut   
    real   , intent(out),   dimension(:,:,:,:) :: entL_u, entB_u, detL_u, detB_u, entr_u, detr_u, buoy_u  
                                                  ! entr_u and detr_u moved to output; buoy_u added 10/10/2017
    integer, intent(in)                        :: istep
    
    real   , dimension(size(thl_upd,1),size(thl_upd,2),size(thl_upd,3),size(thl_upd,4)) :: &  ! updr_new
                                                  thl_utmp, qt_utmp, ql_utmp, w_utmp, a_utmp, thv_utmp, &
                                                  thl_uin,  qt_uin,  ql_uin,  w_uin,  a_uin,  &   
                                                  m_u, chic_u, tend_thl, tend_qt, tend_ql, buoy_prev
    
    real   , dimension(size(thl_upd,1),size(thl_upd,2),size(thl_upd,3)) :: thv_env
    real   , dimension(size(thl_upd,1),size(thl_upd,2)) :: thl_std, qt_std
    real   , dimension(size(thl_upd,1),size(thl_upd,2)) :: tl_2d, qsat_2d, dqsat_2d ! only for surf updraft adjustment
    
    integer :: nlev, nupd, iupd, ilev
    integer :: iter
    
    nlev = size(thl_upd,3)
    nupd = size(thl_upd,4)
    
    ! NOTE BY ZTAN 10/10/2017: Need to use the updated thl_uin and qt_uin for computation (instead of thl_upd and qt_upd) !!!
    thl_uin = thl_upd; qt_uin = qt_upd; ql_uin = ql_upd; w_uin = w_upd; a_uin = a_upd  ! Make a local copy
    
    ! Added by ZTAN 09/25/2017: adjust thl_upd and qt_upd at the surface to be the Gaussian tail value
    call compute_surface_std(wthl_surf, ustar, zfull(:,:,nlev), oblength, thl_std)
    call compute_surface_std(wq_surf,   ustar, zfull(:,:,nlev), oblength, qt_std)
    
    ! do iupd = 1, nupd 
    !     thl_uin (:,:,nlev,iupd)  = thl_grid(:,:,nlev) + &
    !                                thl_std(:,:) * init_scl_u(iupd) * gaussian_std * var_surf_fac
    !     qt_uin  (:,:,nlev,iupd)  = qt_grid (:,:,nlev)  + &
    !                                qt_std(:,:)  * init_scl_u(iupd) * gaussian_std * var_surf_fac
    !     tl_2d = thl_uin (:,:,nlev,iupd)/((pstd_mks/pfull(:,:,nlev))**kappa)
    !     call compute_ql_2d(tl_2d, qt_uin   (:,:,nlev,iupd), ql_uin   (:,:,nlev,iupd), & 
    !                        pfull(:,:,nlev), qsat_2d, dqsat_2d)
    ! end do
    
    do iupd = 1, nupd 
        do ilev = nlev, 1, -1
            if (updr_surfht == 0.0 .and. ilev == nlev) then
                thl_uin (:,:,ilev,iupd)  = thl_grid(:,:,ilev) + &
                                   thl_std(:,:) * init_scl_u(iupd) * gaussian_std * var_surf_fac
                qt_uin  (:,:,ilev,iupd)  = qt_grid (:,:,ilev)  + &
                                   qt_std(:,:)  * init_scl_u(iupd) * gaussian_std * var_surf_fac
                tl_2d = thl_uin (:,:,ilev,iupd)/((pstd_mks/pfull(:,:,ilev))**kappa)
                call compute_ql_2d(tl_2d, qt_uin   (:,:,ilev,iupd), ql_uin   (:,:,ilev,iupd), & 
                               pfull(:,:,ilev), qsat_2d, dqsat_2d)
            else
                where (zfull(:,:,ilev) .lt. updr_surfht)
                
                    thl_uin (:,:,ilev,iupd)  = thl_grid(:,:,ilev) + &
                                   thl_std(:,:) * init_scl_u(iupd) * gaussian_std * var_surf_fac
                    qt_uin  (:,:,ilev,iupd)  = qt_grid (:,:,ilev)  + &
                                   qt_std(:,:)  * init_scl_u(iupd) * gaussian_std * var_surf_fac
                ! elsewhere, don't do 'swap'
                end where
                
                if (minval(zfull(:,:,ilev)) .lt. updr_surfht) then
                    tl_2d = thl_uin (:,:,ilev,iupd)/((pstd_mks/pfull(:,:,ilev))**kappa)
                    call compute_ql_2d(tl_2d, qt_uin   (:,:,ilev,iupd), ql_uin   (:,:,ilev,iupd), & 
                               pfull(:,:,ilev), qsat_2d, dqsat_2d)
                end if
            endif
            
            
        end do
    end do
        
    ! Added by ZTAN 09/25/2017: Adjust thl, qt, ql for the environment
    ! Compute decomp to get environmental parameters ?n, wn, an for the previous model step
    call envupd_decompose(thl_grid, qt_grid, ql_grid, w_grid,         u_grid, v_grid, tke_grid, &
                          thl_uin,  qt_uin,  ql_uin,  w_uin,  a_uin,                            & 
                          thl_env,  qt_env,  ql_env,  w_env,  a_env,  u_env,  v_env,  tke_env)
    
    thl_unew = thl_uin; qt_unew = qt_uin; ql_unew = ql_uin; w_unew = w_uin; a_unew = a_uin
    
    do iter = 1, updr_iter
        thl_utmp = thl_unew; qt_utmp = qt_unew; ql_utmp = ql_unew; w_utmp = w_unew; a_utmp = a_unew 
        
        call compute_upd_buoy (  thl_utmp, qt_utmp, ql_utmp, a_utmp,  &   ! Use updraft (thl, qt, ql, a) from utmp
                                 thl_env,  qt_env,  ql_env,  a_env,   & 
                                 pfull, buoy_u, thv_utmp, thv_env  )
        
        if (do_buoy_avg) then
            if (iter == 1) then
                if (istep .ne. 1) then
                    buoy_u(:,:,:,:) = (buoy_u(:,:,:,:) + buoy_o_u(:,:,:,:)) * .5
                    ! for istep = 1, no averaging is done without iteration
                endif
            else
                buoy_u(:,:,:,:) = (buoy_u(:,:,:,:) + buoy_prev(:,:,:,:)) * .5
            endif
            buoy_prev(:,:,:,:) = buoy_u(:,:,:,:)
        endif
                                 
        call compute_chi_c     ( thl_utmp, qt_utmp, ql_utmp, thv_utmp, &
                                 thl_env,  qt_env,  ql_env,  thv_env,  & 
                                 pfull, chic_u)   
                                 
        call compute_ent_det (w_utmp, buoy_u, chic_u, zstar, zfull, nlev, nupd, entr_u, detr_u)  
            ! Use updraft (w) from utmp
        
        if (do_wu_first) then                                                                   
            call compute_wu    (w_uin,  w_env, buoy_u, entr_u, zfull, zhalf, dt_real, nlev, nupd, w_unew) 
            ! Use updraft (w) from PREV step - 'uin'
                                                                           
            call compute_au    (w_unew, a_uin, rhot, entr_u, detr_u, pfull, zfull, zhalf, dt_real, nlev, nupd, &
                           a_unew, m_u, entL_u, entB_u, detL_u, detB_u)    
            ! Use updraft (a) from PREV step - 'upd'
        else
            if (entr_opt == 0) then
                call compute_au_wu ( a_uin, w_uin, w_env, rhot, buoy_u, entr_u, detr_u, &
                                     pfull, zfull, zhalf, dt_real, nlev, nupd, &
                                     a_unew, w_unew, m_u, entL_u, entB_u, detL_u, detB_u)
            elseif (entr_opt == 1) then
                call compute_au_wu_expentr ( a_uin, w_uin, w_env, rhot, buoy_u, entr_u, detr_u, &
                                     pfull, zfull, zhalf, dt_real, nlev, nupd, &
                                     a_unew, w_unew, m_u, entL_u, entB_u, detL_u, detB_u, .false.)
            elseif (entr_opt == 2) then
                call compute_au_wu_expentr ( a_uin, w_uin, w_env, rhot, buoy_u, entr_u, detr_u, &
                                     pfull, zfull, zhalf, dt_real, nlev, nupd, &
                                     a_unew, w_unew, m_u, entL_u, entB_u, detL_u, detB_u, .true.)
            
            elseif (entr_opt == 10) then
                call compute_au_wu_newp ( a_uin, w_uin, w_env, rhot, buoy_u, entr_u, detr_u, &
                                     pfull, zfull, zhalf, dt_real, nlev, nupd, &
                                     a_unew, w_unew, m_u, entL_u, entB_u, detL_u, detB_u)
            elseif (entr_opt == 11) then
                call compute_au_wu_expentr_newp ( a_uin, w_uin, w_env, rhot, buoy_u, entr_u, detr_u, &
                                     pfull, zfull, zhalf, dt_real, nlev, nupd, &
                                     a_unew, w_unew, m_u, entL_u, entB_u, detL_u, detB_u, .false.)
            else
                call error_mesg('compute_au_wu', 'Invalid entr_opt', FATAL)
            end if
        endif
                           
        call compute_phiu (w_unew, a_unew, rhot, entL_u, entB_u, detL_u, detB_u, &
                           a_uin, thl_uin, qt_uin, ql_uin, thl_env, qt_env, ql_env, thl_grid, qt_grid, & 
                           ustar, wthl_surf, wq_surf, wstar, oblength, & ! oblength added 09/25/2017 
                           pfull, zfull, zhalf, dt_real, nlev, nupd, &
                           thl_unew, qt_unew, ql_unew, tend_thl, tend_qt, tend_ql)  
            ! Use updraft (a, thl, qt, ql) from PREV step - 'upd'
            ! tend_thl, tend_qt, and tend_ql are grid-mean tendencies.
    end do
    
    ! thl_upd = thl_unew; qt_upd = qt_unew; ql_upd = ql_unew; w_upd = w_unew; a_upd = a_unew 
    
    call compute_bulk_updraft(thl_unew, qt_unew, ql_unew, a_unew, tend_thl, tend_qt, tend_ql, &
                              buoy_u, chic_u, m_u, entL_u, entB_u, detL_u, detB_u, & 
                              thl_env, qt_env, ql_env, rhot, zhalf, phalf, nlev, nupd, &
                              thl_ut, qt_ut, ql_ut, w_ut, a_ut, & 
                              m_ut, buoy_ut, chic_ut, entL_ut, entB_ut, detL_ut, detB_ut, &
                              tend_thl_ut, tend_qt_ut, tend_ql_ut)
    
    ! Debug
    if (debug_point .and. dt_real .gt. 200.0) then 
         if (istep == 1) then
            write(*,*) init_scl_u(1), gaussian_std, var_surf_fac
         endif
         
         write(*,"(I6, 6F8.3)")   istep, wthl_surf(idb,jdb)*Cp_air, wq_surf(idb,jdb)*HLv, &
                                  thl_std(idb,jdb), qt_std(idb,jdb)*1.e3, oblength(idb,jdb), wstar(idb,jdb)
         write(*,"(A, 15F8.3)") '   L30', thl_grid(idb,jdb,30),      thl_utmp(idb,jdb,30,1),      &
                                          thl_env (idb,jdb,30),      thl_unew(idb,jdb,30,1),      &
                                          qt_grid (idb,jdb,30)*1.e3, qt_utmp (idb,jdb,30,1)*1.e3, &
                                          qt_env  (idb,jdb,30)*1.e3, qt_unew (idb,jdb,30,1)*1.e3, &
                                                                     thv_utmp(idb,jdb,30,1),      &
                                          thv_env (idb,jdb,30),      buoy_u  (idb,jdb,30,1)*1.e3, &
                                          a_utmp  (idb,jdb,30,1),    w_utmp  (idb,jdb,30,1),      &
                                          a_unew  (idb,jdb,30,1),    w_unew  (idb,jdb,30,1)
                                                         
         write(*,"(A, 15F8.3)") '   L29', thl_grid(idb,jdb,29),      thl_utmp(idb,jdb,29,1),      &
                                          thl_env (idb,jdb,29),      thl_unew(idb,jdb,29,1),      &
                                          qt_grid (idb,jdb,29)*1.e3, qt_utmp (idb,jdb,29,1)*1.e3, &
                                          qt_env  (idb,jdb,29)*1.e3, qt_unew (idb,jdb,29,1)*1.e3, &
                                                                     thv_utmp(idb,jdb,29,1),      &
                                          thv_env (idb,jdb,29),      buoy_u  (idb,jdb,29,1)*1.e3, &
                                          a_utmp  (idb,jdb,29,1),    w_utmp  (idb,jdb,29,1),      &
                                          a_unew  (idb,jdb,29,1),    w_unew  (idb,jdb,29,1)
                                                         
         write(*,"(A, 15F8.3)") '   L28', thl_grid(idb,jdb,28),      thl_utmp(idb,jdb,28,1),      &
                                          thl_env (idb,jdb,28),      thl_unew(idb,jdb,28,1),      &
                                          qt_grid (idb,jdb,28)*1.e3, qt_utmp (idb,jdb,28,1)*1.e3, &
                                          qt_env  (idb,jdb,28)*1.e3, qt_unew (idb,jdb,28,1)*1.e3, &
                                                                     thv_utmp(idb,jdb,28,1),      &
                                          thv_env (idb,jdb,28),      buoy_u  (idb,jdb,28,1)*1.e3, &
                                          a_utmp  (idb,jdb,28,1),    w_utmp  (idb,jdb,28,1),      &
                                          a_unew  (idb,jdb,28,1),    w_unew  (idb,jdb,28,1)
                                          
         write(*,"(A, 15F8.3)") '   L27', thl_grid(idb,jdb,27),      thl_utmp(idb,jdb,27,1),      &
                                          thl_env (idb,jdb,27),      thl_unew(idb,jdb,27,1),      &
                                          qt_grid (idb,jdb,27)*1.e3, qt_utmp (idb,jdb,27,1)*1.e3, &
                                          qt_env  (idb,jdb,27)*1.e3, qt_unew (idb,jdb,27,1)*1.e3, &
                                                                     thv_utmp(idb,jdb,27,1),      &
                                          thv_env (idb,jdb,27),      buoy_u  (idb,jdb,27,1)*1.e3, &
                                          a_utmp  (idb,jdb,27,1),    w_utmp  (idb,jdb,27,1),      &
                                          a_unew  (idb,jdb,28,1),    w_unew  (idb,jdb,28,1)
    end if
    
                      
end subroutine compute_updr


! compute_bulk_updraft: sum across multiple updrafts to get bulk properties
subroutine compute_bulk_updraft(thl_upd, qt_upd, ql_upd, a_upd, tend_thl, tend_qt, tend_ql, &
           buoy_u, chic_u, m_u, entL_u, entB_u, detL_u, detB_u, &
           thl_env, qt_env, ql_env, rhot, zhalf, phalf, nlev, nupd, &
           thl_ut, qt_ut, ql_ut, w_ut, a_ut, &
           m_ut, buoy_ut, chic_ut, entL_ut, entB_ut, detL_ut, detB_ut, &
           tend_thl_ut, tend_qt_ut, tend_ql_ut)
           
    real   , intent(in), dimension(:,:,:,:) :: thl_upd, qt_upd, ql_upd, a_upd, tend_thl, tend_qt, tend_ql
    real   , intent(in), dimension(:,:,:,:) :: buoy_u, chic_u, m_u, entL_u, entB_u, detL_u, detB_u
    real   , intent(in), dimension(:,:,:)   :: thl_env, qt_env, ql_env, rhot, zhalf, phalf
    integer, intent(in)                     :: nlev, nupd
    real   , intent(out),dimension(:,:,:)   :: thl_ut, qt_ut, ql_ut, w_ut, a_ut, &
                                               m_ut, buoy_ut, chic_ut, entL_ut, entB_ut, detL_ut, detB_ut, &
                                               tend_thl_ut, tend_qt_ut, tend_ql_ut 
    integer :: iupd, ilev
    
    a_ut = 0.; w_ut = 0.; m_ut = 0.; entL_ut = 0.; entB_ut = 0.; detL_ut = 0.; detB_ut = 0.
    thl_ut = 0.; qt_ut = 0.; ql_ut = 0.; buoy_ut = 0.; chic_ut = 0.
    tend_thl_ut = 0.; tend_qt_ut = 0.; tend_ql_ut = 0.
    
    do iupd = 1, nupd
    ! Simple-summation variables: 
           a_ut(:,:,:) =    a_ut(:,:,:) +   a_upd(:,:,:,iupd)
           m_ut(:,:,:) =    m_ut(:,:,:) +     m_u(:,:,:,iupd)
        entL_ut(:,:,:) = entL_ut(:,:,:) +  entL_u(:,:,:,iupd)
        entB_ut(:,:,:) = entB_ut(:,:,:) +  entB_u(:,:,:,iupd)
        detL_ut(:,:,:) = detL_ut(:,:,:) +  detL_u(:,:,:,iupd)
        detB_ut(:,:,:) = detB_ut(:,:,:) +  detB_u(:,:,:,iupd)
        
    
    ! Summed variables: tend_thl, tend_qt, tend_ql
        tend_thl_ut(:,:,:) = tend_thl_ut(:,:,:) + tend_thl(:,:,:,iupd)
         tend_qt_ut(:,:,:) =  tend_qt_ut(:,:,:) +  tend_qt(:,:,:,iupd)
         tend_ql_ut(:,:,:) =  tend_ql_ut(:,:,:) +  tend_ql(:,:,:,iupd)
        
    ! MF-weighted variables: thl, qt, ql, buoy, chic (weighted sum)
         thl_ut(:,:,:) =  thl_ut(:,:,:) + thl_upd(:,:,:,iupd) * m_u(:,:,:,iupd)
          qt_ut(:,:,:) =   qt_ut(:,:,:) +  qt_upd(:,:,:,iupd) * m_u(:,:,:,iupd)
          ql_ut(:,:,:) =   ql_ut(:,:,:) +  ql_upd(:,:,:,iupd) * m_u(:,:,:,iupd)
        buoy_ut(:,:,:) = buoy_ut(:,:,:) +  buoy_u(:,:,:,iupd) * m_u(:,:,:,iupd)
        chic_ut(:,:,:) = chic_ut(:,:,:) +  chic_u(:,:,:,iupd) * m_u(:,:,:,iupd)
        
    end do  
    
    ! Area-weighted variables: w
    where (a_ut(:,:,:) .gt. 0.0)
        w_ut(:,:,:) = m_ut(:,:,:)/a_ut(:,:,:)/rhot(:,:,:)        
    end where
    
    ! MF-weighted variables: thl, qt, ql, buoy, chic (divided by weight)
    where (m_ut(:,:,:) .gt. 0.0)
         thl_ut(:,:,:) =  thl_ut(:,:,:) / m_ut(:,:,:) 
          qt_ut(:,:,:) =   qt_ut(:,:,:) / m_ut(:,:,:)  
          ql_ut(:,:,:) =   ql_ut(:,:,:) / m_ut(:,:,:)  
        buoy_ut(:,:,:) = buoy_ut(:,:,:) / m_ut(:,:,:)  
        chic_ut(:,:,:) = chic_ut(:,:,:) / m_ut(:,:,:) 
    elsewhere ! No mass flux, fill in environment values
         thl_ut(:,:,:) = thl_env(:,:,:)
          qt_ut(:,:,:) =  qt_env(:,:,:)
          ql_ut(:,:,:) =  ql_env(:,:,:)
        ! chic_ut and buoy_ut are filled with values of the 1st updraft   
        buoy_ut(:,:,:) = buoy_u(:,:,:,1)  
        chic_ut(:,:,:) = chic_u(:,:,:,1)  
    end where

end subroutine compute_bulk_updraft


! compute_phiu: compute updated updraft scalar values phiu
subroutine compute_phiu (w, a, rhot, entL, entB, detL, detB, &
           a_o, thl_o, qt_o, ql_o, thl_e, qt_e, ql_e, thl_g, qt_g, &
           ustar, wthl_s, wq_s, wstar, oblength, pfull, zfull, zhalf, dt, nlev, nupd, &
           thl, qt, ql, s_thl, s_qt, s_ql)
! init_scl_u(nupd) is now a global array

    real   , intent(in),    dimension(:,:,:,:) :: w, a, entL, entB, detL, detB, a_o, thl_o, qt_o, ql_o
    real   , intent(in),    dimension(:,:,:)   :: rhot, thl_e, qt_e, ql_e, thl_g, qt_g, pfull, zfull, zhalf
    real   , intent(in) ,   dimension(:,:)     :: ustar, wthl_s, wq_s, wstar, oblength
    real   , intent(in)                        :: dt
    integer, intent(in)                        :: nlev, nupd
    real   , intent(out),   dimension(:,:,:,:) :: thl, qt, ql, s_thl, s_qt, s_ql

    real   , dimension(size(zfull,1),size(zfull,2))  :: dz, thl_ent, qt_ent, & ! tke_s, &
                                 c1, c2, c3, ct, d_thl, d_qt, d_ql, thl_std, qt_std
    integer :: iupd, ilev
    ! Note: phiu is not equal to phie when au first reaches 0, but at the next step it will be zero.
    
    thl = 0.; qt = 0.; ql = 0.; s_thl = 0.; s_qt = 0.; s_ql = 0.; d_thl = 0.; d_qt = 0.; d_ql = 0.
    
    ! call compute_surface_tke(ustar, wstar, tke_s)  ! This is no longer used by the new formulation
    ! tke_s(:,:) = max( (3.75*(ustar(:,:)**2.) + 0.2*(wstar(:,:)**2.)), 1.0e-8)
    
    ! Added ZTAN 09/25/2017: compute surface covariance according to manuscript
    ! Note: wthl_s and wqt_s are NOT density-weighted
    
    call compute_surface_std(wthl_s(:,:), ustar(:,:), zfull(:,:,nlev), oblength(:,:), thl_std(:,:))
    call compute_surface_std(wq_s(:,:),   ustar(:,:), zfull(:,:,nlev), oblength(:,:), qt_std(:,:)) 
    
    do iupd = 1, nupd
        do ilev = nlev, 1, -1
            dz(:,:) = zhalf(:,:,ilev) - zhalf(:,:,ilev+1)
            if (updr_surfht == 0.0 .and. ilev == nlev) then
                
                ! NOTE by ZTAN 09/25/2017: thl_g and qt_g should be used in place of thl_e and qt_e 
                !                          when calculating updraft boundary conditions !!!
                ! thl_ent(:,:) = thl_e(:,:,ilev) + & 
                !     1.8*wthl_s(:,:)/max(tke_s(:,:)**0.5, wstar(:,:)) * init_scl_u(iupd)* var_surf_fac
                ! qt_ent(:,:)  = qt_e (:,:,ilev) + & 
                !     1.8*wq_s(:,:)  /max(tke_s(:,:)**0.5, wstar(:,:)) * init_scl_u(iupd)* var_surf_fac
                
                thl_ent(:,:) = thl_g(:,:,ilev) + &
                               thl_std(:,:) * init_scl_u(iupd) * gaussian_std * var_surf_fac
                qt_ent(:,:)  = qt_g(:,:,ilev)  + &
                               qt_std(:,:)  * init_scl_u(iupd) * gaussian_std * var_surf_fac
            else 
                where (zfull(:,:,ilev) .lt. updr_surfht)
                    ! thl_ent(:,:) = thl_e(:,:,ilev) + & 
                    !     1.8*wthl_s(:,:)/max(tke_s(:,:)**0.5, wstar(:,:)) * init_scl_u(iupd)* var_surf_fac
                    ! qt_ent(:,:)  = qt_e (:,:,ilev) + & 
                    !     1.8*wq_s(:,:)  /max(tke_s(:,:)**0.5, wstar(:,:)) * init_scl_u(iupd)* var_surf_fac
                    
                    thl_ent(:,:) = thl_g(:,:,ilev) + &
                                   thl_std(:,:) * init_scl_u(iupd) * gaussian_std * var_surf_fac
                    qt_ent(:,:)  = qt_g(:,:,ilev)  + &
                                   qt_std(:,:)  * init_scl_u(iupd) * gaussian_std * var_surf_fac
                elsewhere
                    thl_ent(:,:) = thl_e(:,:,ilev)
                    qt_ent(:,:)  = qt_e (:,:,ilev)
                end where
            end if
            

            c1(:,:) = rhot(:,:,ilev)*a_o(:,:,ilev,iupd)*dz(:,:)*updraft_rescale/dt
            c3(:,:) = entL(:,:,ilev, iupd)+entB(:,:,ilev, iupd)
            if (ilev == nlev) then
                c2(:,:) = 0.0
            else
                c2(:,:) = rhot(:,:,ilev+1)*a(:,:,ilev+1, iupd)*w(:,:,ilev+1, iupd)
            end if
            ct(:,:) = c1(:,:) + c2(:,:) + c3(:,:)
            
            
            if (ilev == nlev) then
                where (ct(:,:) > 0.0)
                    c1(:,:) = c1(:,:)/ct(:,:); c2(:,:) = c2(:,:)/ct(:,:); c3(:,:) = c3(:,:)/ct(:,:)
                    thl(:,:,ilev, iupd) = c1(:,:)*thl_o(:,:,ilev,iupd) + c3*thl_ent(:,:)
                    qt (:,:,ilev, iupd) = c1(:,:)*qt_o (:,:,ilev,iupd) + c3*qt_ent (:,:)
                    ql (:,:,ilev, iupd) = c1(:,:)*ql_o (:,:,ilev,iupd) ! No entrainment??
                end where
                
                ! ZTAN 09/25/2017: changed surface condition -- detrain everything and 
                ! entrain the most 'buoyant' airmass
                ! Note: The effective entrainment rate is different from the w-equation
                ! c1(:,:) = 0.0; c2(:,:) = 0.0; c3(:,:) = 1.0
                ! ct(:,:) = rhot(:,:,ilev)*a(:,:,ilev,iupd)*dz(:,:)*updraft_rescale/dt
                ! thl(:,:,ilev, iupd) = thl_ent(:,:)
                ! qt (:,:,ilev, iupd) = qt_ent (:,:)
                ! ql (:,:,ilev, iupd) = 0.0  ! ql will be later calculated by compute_precip_upd 
                
            else
                where (ct(:,:) > 0.0)
                    c1(:,:) = c1(:,:)/ct(:,:); c2(:,:) = c2(:,:)/ct(:,:); c3(:,:) = c3(:,:)/ct(:,:)
                    thl(:,:,ilev, iupd) = c1(:,:)*thl_o(:,:,ilev,iupd) + c2(:,:)*thl(:,:,ilev+1,iupd) + c3(:,:)*thl_ent(:,:)
                    qt (:,:,ilev, iupd) = c1(:,:)*qt_o (:,:,ilev,iupd) + c2(:,:)*qt (:,:,ilev+1,iupd) + c3(:,:)*qt_ent (:,:)
                    ql (:,:,ilev, iupd) = c1(:,:)*ql_o (:,:,ilev,iupd) + c2(:,:)*ql (:,:,ilev+1,iupd) + c3(:,:)*ql_e(:,:,ilev)
                end where
            end if
            
            ! Call saturation adjustment in updraft  <- 'where' does not work... Redundant calculation including environment?
            if (maxval(ct(:,:)) > 0.0) then 
                call compute_precip_upd(thl(:,:,ilev,iupd), qt(:,:,ilev,iupd), ql(:,:,ilev,iupd), pfull(:,:,ilev), dt, &
                                d_thl(:,:), d_qt(:,:), d_ql(:,:)) ! d_thl, d_qt, d_ql are in K, kg/kg, kg/kg
            else
                d_thl = 0.; d_qt = 0.; d_ql = 0.
            end if
                
            ! ct is in kg/m^3 * m/s; thus s_thl is in K/s, and s_qt and s_ql are in kg/kg/s.          
            where (ct(:,:) > 0.0)   
                s_thl(:,:,ilev, iupd) = ct(:,:)/dz(:,:)/rhot(:,:,ilev) * d_thl(:,:) ! fill in d_thl, d_qt, d_ql
                s_qt (:,:,ilev, iupd) = ct(:,:)/dz(:,:)/rhot(:,:,ilev) * d_qt (:,:)
                s_ql (:,:,ilev, iupd) = ct(:,:)/dz(:,:)/rhot(:,:,ilev) * d_ql (:,:)
                
                thl(:,:,ilev, iupd) = thl(:,:,ilev, iupd)+d_thl(:,:)
                qt (:,:,ilev, iupd) = qt (:,:,ilev, iupd)+d_qt (:,:)
                ql (:,:,ilev, iupd) = ql (:,:,ilev, iupd)+d_ql (:,:)
                
            elsewhere
                s_thl(:,:,ilev, iupd) = 0.0
                s_qt (:,:,ilev, iupd) = 0.0
                s_ql (:,:,ilev, iupd) = 0.0
                
                thl(:,:,ilev, iupd) = thl_e(:,:,ilev) 
                qt (:,:,ilev, iupd) = qt_e (:,:,ilev)
                ql (:,:,ilev, iupd) = ql_e (:,:,ilev)
            end where
            
        end do
    end do     

end subroutine compute_phiu


! compute_precip_upd: do simple saturation adjustment in updraft <- should change to more sophisticated microphysics in the future
subroutine compute_precip_upd (thl_o, qt_o, ql_o, p, dt, d_thl, d_qt, d_ql)

    real   , intent(in),    dimension(:,:)   :: thl_o, qt_o, ql_o, p
    real   , intent(out),   dimension(:,:)   :: d_thl, d_qt, d_ql
    real   , intent(in)                      :: dt
    
    real   , dimension(size(thl_o,1),size(thl_o,2))  :: Tl, thl, qt, ql, qs, dqs, prec
    real :: ct = 1.0e-4
    real :: cw = 5.0e-4
    
    Tl(:,:) = thl_o(:,:)*(p(:,:)/pstd_mks)**kappa
    thl(:,:) = thl_o(:,:); qt(:,:) = qt_o(:,:); ql(:,:) = ql_o(:,:)
    
    call compute_ql_2d(Tl(:,:), qt(:,:), ql(:,:), p(:,:), qs(:,:), dqs(:,:))
    
    prec(:,:) = 0.
    
    if (precip_upd_opt == 1) then      ! Precip threshold with fixed fraction (precip_upd_fac)
        prec(:,:) = max(0., ql(:,:) - precip_upd_fac * qs(:,:))
    elseif (precip_upd_opt == 2) then  ! Precip threshold with fixed value    (precip_upd_val)
        prec(:,:) = max(0., ql(:,:) - precip_upd_val)
    elseif (precip_upd_opt == 3) then  ! New precip scheme 
        prec(:,:) = - dt * ql(:,:) * (ct * (1. - exp(-(ql(:,:)/cw))) )
        prec(:,:) = min(ql(:,:), prec(:,:))
    end if
    
     thl(:,:) = thl(:,:) * (1. + HLv * prec(:,:)/(Cp_air*Tl(:,:))) 
      qt(:,:) =  qt(:,:) - prec(:,:)
      ql(:,:) =  ql(:,:) - prec(:,:)
    
    d_thl(:,:) = thl(:,:) - thl_o(:,:)
    d_qt(:,:)  =  qt(:,:) -  qt_o(:,:)
    d_ql(:,:)  =  ql(:,:) -  ql_o(:,:)

end subroutine compute_precip_upd

! Added ZTAN: 09/20/2017
! compute_au_wu: compute updated updraft fraction au and updraft velocity wu, with other auxiliary variables
subroutine compute_au_wu (a_old, w_old, w_env, rhot, buoy, entr, detr, &
                          pfull, zfull, zhalf, dt, nlev, nupd, &
                          a, w, m, entL, entB, detL, detB)

    real   , intent(in),    dimension(:,:,:,:) :: a_old, w_old, buoy, entr, detr
    real   , intent(in),    dimension(:,:,:)   :: w_env, rhot, pfull, zfull, zhalf
    real   , intent(in)                        :: dt
    integer, intent(in)                        :: nlev, nupd
    real   , intent(out),   dimension(:,:,:,:) :: a, w, m, entL, entB, detL, detB
    
    real   , dimension(size(a_old,1),size(a_old,2))  :: dz, c1, c2, d3, d4, db, w1, w2, w3, &
                                                        w_adv, coef_0, coef_1, coef_2, a_denom, &
                                                        c1s, c1d, atd, f_1, f_2
    integer, dimension(size(a_old,1),size(a_old,2))  :: adjust_flag
    
    ! c1: old mass; c2: mass from below; d3: (entrained mass)/wu; d4: (entrained - detrained mass)/wu.
    ! w1: old wu; w2: wu from below; w3: w of entrained mass.
    ! w_adv: wu without buoyancy effect, (coef_0, coef_1, coef_2): coefficients of the quadratic equation.
    ! adjust_flag: flag for updraft fraction limiter.
    
    
    integer :: iupd, ilev
    real :: b1  = 2.5
    real :: b2  = 2.0
    real :: limit_fac = 100.

    a = 0.; w = 0.; m = 0.; entL = 0.; entB = 0.; detL = 0.; detB = 0.
    dz = 0.; c1 = 0.; c2 = 0.; d3 = 0.; d4 = 0.; db = 0.; w1 = 0.; w2 = 0.; w3 = 0.
    w_adv = 0.; coef_0 = 0.; coef_1 = 0.; coef_2 = 0.; a_denom = 0.
    c1s = 0.; c1d = 0.; atd = 0.; f_1 = 0.; f_2 = 0.
    
    do iupd = 1, nupd
        do ilev = nlev, 1, -1
            dz(:,:) = zhalf(:,:,ilev) - zhalf(:,:,ilev+1)
            c1(:,:) = rhot(:,:,ilev)*a_old(:,:,ilev,iupd)*dz(:,:)*updraft_rescale/dt
            w1(:,:) = w_old(:,:,ilev,iupd)
            w3(:,:) = w_env(:,:,ilev)
            
            if (ilev /= nlev) then
                w2(:,:) = w(:,:,ilev+1,iupd)
                c2(:,:) = rhot(:,:,ilev+1)*a(:,:,ilev+1, iupd)*w(:,:,ilev+1, iupd)
                
                ! Compute entrained and detrained mass coefficients
                if (au_optL == 1) then  ! use exponential formula (119)
                    ! Note: Substituting (entr-detr)*dz by exp[(entr-detr)*dz]-1.
                    ! This may be unstable if (entr-detr)*dz is not small (e.g., near SFC). 
                    where (entr(:,:,ilev,iupd)-detr(:,:,ilev, iupd) == 0.0)
                        d4(:,:) = 0.
                        d3(:,:) = rhot(:,:,ilev)*entr(:,:,ilev, iupd)*dz(:,:)
                    elsewhere
                        d4(:,:) = rhot(:,:,ilev)*(exp((entr(:,:,ilev, iupd)-detr(:,:,ilev, iupd))*dz(:,:))-1.)
                        d3(:,:) = entr(:,:,ilev, iupd)/(entr(:,:,ilev, iupd)-detr(:,:,ilev, iupd))*d4(:,:) 
                    end where
                else  ! use simple formula
                    d4(:,:) = rhot(:,:,ilev)*(entr(:,:,ilev,iupd)-detr(:,:,ilev, iupd))*dz(:,:)
                    d3(:,:) = rhot(:,:,ilev)* entr(:,:,ilev,iupd)*dz(:,:)
                endif
                
                where (d4(:,:) > rhot(:,:,ilev))  
                    ! MAY NEED TO LIMIT d4 < rhot because of implicit scheme.
                    d3(:,:) = d3(:,:)/d4(:,:)* rhot(:,:,ilev) ! now d3 = 
                    d4(:,:) = rhot(:,:,ilev)
                end where
                where (d4(:,:) < - limit_fac * rhot(:,:,ilev)) 
                    ! LIMIT d4 > - limit_fac * rhot? May be unnecessary... ZTAN 09/26/2017
                    d3(:,:) = d3(:,:)/d4(:,:)*(- limit_fac * rhot(:,:,ilev))
                    d4(:,:) = - limit_fac * rhot(:,:,ilev)
                end where
            else
                c2(:,:) = 0.; w2(:,:) = 0.; d3(:,:) = 0.; d4(:,:) = 0.
            endif  ! (ilev /= nlev)
            
            ! write(*,*), ilev, maxval(d3(:,:)), minval(d3(:,:)), maxval(d4(:,:)), minval(d4(:,:))
            db(:,:) = rhot(:,:,ilev)*dz(:,:)

            ! ------ Step 1 ------
            ! Assume entB = 0 and detB = 0, and calculate {w},{a},{entL},{detL}.
            ! w-equation: (c1 + c2 + d3*{a}*{w} * b1) * {w} = ...
            !             (c1 * w1) + (c2 * w2) + (d3*{a}*{w} * b1 * w3) + ...
            !             (db * b2) *{a}* buoy
            ! a-equation: (c1 + c2 + d4*{a}*{w}) = ...
            !              rhot*{a} *dz*updraft_rescale/dt + rhot*{a}*{w}
            ! Solve {w} first, then {a}, then compute entL and detL if {a} satisfies condition.
            
            where (c1(:,:) + c2(:,:) <= 0.)  !   no prev updraft and no incoming updraft
                ! ==> no current updraft without massive entrainment. a_old should be zero here.
                w(:,:,ilev,iupd) = 0.; a(:,:,ilev,iupd) = a_old(:,:,ilev,iupd)
                entL(:,:,ilev,iupd) = 0.; entB(:,:,ilev,iupd) = 0.
                detL(:,:,ilev,iupd) = 0.; detB(:,:,ilev,iupd) = 0.
            elsewhere
                w_adv (:,:) = (c1(:,:)*w1(:,:)+c2(:,:)*w2(:,:))/(c1(:,:)+c2(:,:))
                coef_2(:,:) = b1*d3(:,:) - d4(:,:) + rhot(:,:,ilev) ! >= 0 because d3 >= 0 and we limit d4 <= rhot
                coef_1(:,:) = rhot(:,:,ilev)*dz(:,:)*updraft_rescale/dt - &
                              b1*d3(:,:)*w3(:,:) - (rhot(:,:,ilev)-d4(:,:))* w_adv(:,:)
                coef_0(:,:) = - w_adv(:,:) * rhot(:,:,ilev)*dz(:,:)*updraft_rescale/dt - &
                              b2*db(:,:)*buoy(:,:,ilev, iupd)
                              
                ! where coef_2 < 0 <- removed with the above discussion
                
                where (coef_0(:,:) >= 0) 
                ! The non-diluted w_nent is negative, i.e., buoy << 0. 
                ! Thus updraft stops, boundary detrainment occurs (Step 3).
                    w   (:,:,ilev, iupd) = 0.; a(:,:,ilev, iupd) = a_old(:,:,ilev,iupd) ! a will be adjusted by Step 3
                    entL(:,:,ilev, iupd) = 0.; detL(:,:,ilev, iupd) = 0.
                    entB(:,:,ilev, iupd) = 0.; detB(:,:,ilev, iupd) = 0.
                elsewhere
                    where (coef_2(:,:) == 0)  
                        ! degenerate equation, implies d3 = 0 and d4 = rhot, and thus coef_1 > 0; note coef_0 < 0
                        w(:,:,ilev,iupd) = - coef_0(:,:) / coef_1(:,:) 
                    elsewhere
                        w (:,:,ilev, iupd) = (-coef_1(:,:) + &
                             sqrt(coef_1(:,:)**2. - 4.*coef_2(:,:)*coef_0(:,:)))/(2.*coef_2(:,:))
                    end where
                    a_denom (:,:) =  rhot(:,:,ilev)*dz(:,:)*updraft_rescale/dt + &
                                    (rhot(:,:,ilev)-d4(:,:))*w(:,:,ilev, iupd)
                    where (a_denom(:,:) < 0.)
                        a(:,:,ilev, iupd) = 1.e5 ! Temporarily filling in a large number,
                                                 ! Will be fixed by increased detrainment (detB).
                    elsewhere
                        a(:,:,ilev, iupd) = (c1(:,:) + c2(:,:))/a_denom (:,:)
                        ! At the limit of large entrainment (e.g., b is large but w is small)
                        ! -> a_denom = rho*dz/dt, c1 = rho*a_o*dz/dt, c2 = rho- * a- * w-
                        ! -> a = a_o + ( rho- /rho)*(w- /(dz/dt))* a- . If dz/dt << w-, we will have a >> a_o
                        ! This is expected !
                        ! And: d3 = d4 = rho => 
                        ! coef_2 = (b1-1)*rho, coef_1 = rho*dz/dt - b1*rho*w3, 
                        ! coef_0 = - w_adv * rho *dz/dt - b2*rho*dz*buoy
                        ! At the limit of dz/dt -> 0: 
                        ! coef_2 = (b1-1)*rho, coef_1 = -b1*rho*w3 ~ 0, coef_0 = -b2*rho*dz*buoy
                        ! => w = sqrt(b2*dz*buoy / (b1-1)) ~ sqrt(1.33*dz*buoy)
                        ! The dominant balance is: 
                        
                        ! However, if entrainment is slightly smaller (d4 = 0.9*rho), 
                        ! then a_denom = rho*(dz/dt +0.1w)
                        ! a = (a_o*rho*dz/dt + a- * rho- * w-)/[rho*(dz/dt +0.1w)]
                        ! If dz/dt << w-, w => a ~ a- * 10(w-/w) is still huge...
                        
                    end where
                    entB(:,:,ilev, iupd) = 0.; detB(:,:,ilev, iupd) = 0.
                    entL(:,:,ilev, iupd) =  d3(:,:)          * a(:,:,ilev, iupd) * w(:,:,ilev, iupd)
                    detL(:,:,ilev, iupd) = (d3(:,:)-d4(:,:)) * a(:,:,ilev, iupd) * w(:,:,ilev, iupd)
                end where ! (coef_0(:,:) >= 0) 
            end where ! (c1(:,:) + c2(:,:) <= 0.)


            ! write(*,*), ilev, maxval(entL(:,:)), minval(entL(:,:)), maxval(detL(:,:)), minval(detL(:,:))
            
            if (ilev == nlev) then
            ! ------ Step 2 ------
            ! Fix SFC fraction by:
            ! (1) If buoy > 0 => w > 0, aN < a0 because entL = 0 => entB > 0.
            ! (2) If buoy <= 0. Swap airmass, use detB = c1, calculate entB accordingly.
                ! if (au_optB_srf == 1) then ! Disabled ZTAN 09/26/2017
                !     c1s(:,:) = c1(:,:)
                ! else
                !     c1s(:,:) = 0.0
                ! endif
                
                ! where (w (:,:,ilev, iupd) < 1.0e-6)
                !     a(:,:,ilev, iupd) = a_old(:,:,ilev, iupd)
                !     w(:,:,ilev, iupd) = 0.0
                !     entL(:,:,ilev, iupd) = 0.0; detL(:,:,ilev, iupd) = 0.0
                !     entB(:,:,ilev, iupd) = rhot(:,:,ilev)*a(:,:,ilev, iupd)*dz(:,:)*updraft_rescale/dt - c2(:,:)
                !     detB(:,:,ilev, iupd) = c1(:,:)
                    
                ! elsewhere
                ! Calculate entB for target a(ilev, iupd) = a_o(ilev, iupd)
                ! final w-equation: (c1 + c2 + d3*a*{w} * b1 + {entB} * b1) * {w} = ...
                !                   (c1 * w1) + (c2 * w2) + (d3*a*{w} + {entB}) * b1 * w3 + (db * b2 * {a}) * buoy
                ! final a-equation: (c1 + c2 + d4*a{w} + {entB} - {detB}) = ...
                !                    rhot*a *dz*updraft_rescale/dt + rhot*a*{w}
                !                  -> f_2 * {w} + f_1 = {entB}
                ! 
                ! SFC w-equation at steady state: {w} = w1, w2 = w3 = 0, au = const; 
                !     also: c1 = rho*a*dz/dt, c2 = 0, d3 = d4 = 0, db = rho*dz
                ! final w-equation: {w} = (rho * dz * b2 * {a}) * buoy / ({entB} * b1) 
                ! (1) if the updraft is always swapped -> entB = c1 + rho*{a}*{w}; detB = c1
                !     -> this gives (rho * b1 * {a}) * {w}^2 + c1 * b1 * {w} = (rho * dz * b2 * {a}) * buoy
                !     -> i.e., b1*{w}*(dz/dt+{w}) = dz*b2*buoy (THIS IS PROBLEMATIC, especially if dt is small)
                ! (2) if the updraft is NOT swapped -> entB = rho*a*w; detB = 0.
                !     -> this gives (rho * b1 * {a}) * {w}^2 = (rho * dz * b2 * {a}) * buoy
                !     -> i.e., b1*{w}^2 = dz*b2*buoy
                
                detB(:,:,ilev, iupd) = 0.0 ! c1s(:,:)  ! Changed ZTAN 09/26/2017: Do not swap 
                    
                a(:,:,ilev, iupd) = a_old(:,:,ilev, iupd)
                f_1(:,:) = rhot(:,:,ilev)*a(:,:,ilev, iupd)*dz(:,:)*updraft_rescale/dt - &
                      (c1(:,:)+c2(:,:) - detB(:,:,ilev, iupd))
                f_2(:,:) = (rhot(:,:,ilev) - d4(:,:))*a(:,:,ilev, iupd)

                coef_2(:,:) = b1*(f_2(:,:) + d3(:,:) *a(:,:,ilev, iupd))   ! coef_2 >= 0 by definition
                coef_1(:,:) = c1(:,:)+c2(:,:) + b1*f_1(:,:) - b1*w3(:,:)*(f_2(:,:) + d3(:,:) *a(:,:,ilev, iupd))
                coef_0(:,:) = -(c1(:,:)*w1(:,:) + c2(:,:)*w2(:,:) + b1*f_1(:,:)*w3(:,:) + &
                                    b2*db(:,:)*a(:,:,ilev, iupd)*buoy(:,:,ilev, iupd))
                               
                where (coef_0(:,:) >= 0.)  ! Updraft at SFC is not buoyant
                    w(:,:,ilev, iupd) = 0.0
                    entL(:,:,ilev, iupd) = 0.0; detL(:,:,ilev, iupd) = 0.0
                    entB(:,:,ilev, iupd) = rhot(:,:,ilev)*a(:,:,ilev, iupd)*dz(:,:)*updraft_rescale/dt - c2(:,:)
                    detB(:,:,ilev, iupd) = c1(:,:)
                elsewhere
                    where (coef_2(:,:) == 0)  
                        ! degenerate equation, implies d3 = 0 and d4 = rhot, 
                        ! and thus coef_1 = c1 + b1*f_1 > 0, since f_1 = detB if a_u is fixed; note coef_0 < 0.
                        w(:,:,ilev,iupd) = - coef_0(:,:) / coef_1(:,:) 
                    elsewhere
                        w (:,:,ilev, iupd) = (-coef_1(:,:) + &
                             sqrt(coef_1(:,:)**2. - 4.*coef_2(:,:)*coef_0(:,:)))/(2.*coef_2(:,:))
                    end where
                    entB(:,:,ilev, iupd) = f_1(:,:) + f_2(:,:) * w(:,:,ilev, iupd)
                    entL(:,:,ilev, iupd) =  d3(:,:)          * a(:,:,ilev, iupd) * w(:,:,ilev, iupd)
                    detL(:,:,ilev, iupd) = (d3(:,:)-d4(:,:)) * a(:,:,ilev, iupd) * w(:,:,ilev, iupd)
                end where ! (coef_0 >= 0.)
                    
                ! end where ! (w (:,:,ilev, iupd) < 1.0e-6)
                
                ! Adjust entB and detB for scalars -- switch updraft thl/qt with environment buoyant tail
                ! ZTAN: 09/25/2017: This is disabled due to inconsistency. The 'switch' is now done at the beginning.
                ! detB(:,:,ilev, iupd) = detB(:,:,ilev, iupd) + c1(:,:)*99. ! 99
                ! entB(:,:,ilev, iupd) = entB(:,:,ilev, iupd) + c1(:,:)*99. ! 99
            else
            ! ------ Step 3 ------
            ! Adjust area fraction at other levels by tuning up detL.
            ! Note: this algorithm is applicable to adjusting down only.
                
                ! Limitation on w_u
                adjust_flag(:,:) = 0
                if (au_optB_wu == 1) then ! detraining parcels with zero velocity
                    where ( (w(:,:,ilev, iupd) <= 1.0e-6) .or. (a(:,:,ilev, iupd) <= 1.0e-6) )
                        adjust_flag(:,:) = 2  ! Kill the updraft
                    end where
                elseif (au_optB_wu == 11) then ! detraining parcels with zero velocity or zero area
                    where ( (a(:,:,ilev, iupd) <= 1.0e-6))
                        adjust_flag(:,:) = 2  ! Kill the updraft
                    elsewhere
                        where (w(:,:,ilev, iupd) <= 1.0e-6)
                            adjust_flag(:,:) = 3  ! Kill the updraft gradually
                        end where
                    end where
                endif ! (au_optB_wu == 1) (au_optB_wu == 11)
            
            
                ! Limitation on a_u is top priority -- will overwrite adjust_flag due to w_u.
                if (au_optB == 1) then    ! constrain updraft fraction < SFC fraction
                    where (a(:,:,ilev, iupd)> a(:,:,nlev, iupd))
                        a(:,:,ilev,iupd) = a(:,:,nlev, iupd);  adjust_flag(:,:) = 1
                    end where
                elseif (au_optB == 2) then    ! constrain updraft fraction < 2*SFC fraction
                    where (a(:,:,ilev, iupd)> 2.0*a(:,:,nlev, iupd))
                        a(:,:,ilev,iupd) = 2.0*a(:,:,nlev, iupd);  adjust_flag(:,:) = 1
                    end where
                elseif (au_optB == 100) then  ! constrain updraft fraction < 0.5/nupd
                    where (a(:,:,ilev, iupd)> .5/nupd)
                        a(:,:,ilev,iupd) = .5/nupd;  adjust_flag(:,:) = 1
                    end where
                endif ! (au_optB == 1, 2, 100)
                
                where (adjust_flag(:,:) == 3) ! Adjust w to zero directly, 
                                              ! but detrain a over a 15 minute timescale.
                    w   (:,:,ilev, iupd) = 0. 
                    a   (:,:,ilev, iupd) = a_old(:,:,ilev,iupd) * exp(-dt/900.0)   
                    entL(:,:,ilev, iupd) = 0.; detL(:,:,ilev, iupd) = 0.
                    entB(:,:,ilev, iupd) = 0. 
                    detB(:,:,ilev, iupd) = c1(:,:)+c2(:,:) - &
                         rhot(:,:,ilev)*a(:,:,ilev, iupd)*dz(:,:)*updraft_rescale/dt
                end where ! (adjust_flag(:,:) == 3) 
                
                where (adjust_flag(:,:) == 2) ! Adjust to zero directly.
                    w   (:,:,ilev, iupd) = 0.  
                    a   (:,:,ilev, iupd) = 0.  
                    entL(:,:,ilev, iupd) = 0.; detL(:,:,ilev, iupd) = 0.
                    entB(:,:,ilev, iupd) = 0.
                    detB(:,:,ilev, iupd) = c1(:,:)+c2(:,:)
                end where ! (adjust_flag(:,:) == 2) 
                
                ! c1d is only used where adjust_flag == 1
                if (au_optB_wu == 1) then   ! Detrainment everything. atd is target a.
                    c1d(:,:) = c1(:,:)+c2(:,:)
                    atd(:,:) = 0.
                elseif (au_optB_wu == 11) then ! Detrainment to target a (= atd).
                    c1d(:,:) = c1(:,:)+c2(:,:) &
                               - rhot(:,:,ilev)*a(:,:,ilev, iupd)*dz(:,:)*updraft_rescale/dt
                    atd(:,:) = a(:,:,ilev, iupd)
                endif
                
                where (adjust_flag(:,:) == 1)
                    ! Calculate detB for target a(ilev, iupd) = a_o(ilev, iupd), with entB = 0
                    ! final w-equation: (c1 + c2 + d3*a*{w} * b1) * {w} = ...
                    !                   (c1 * w1) + (c2 * w2) + d3*a*{w} * b1 * w3 + (db*a * b2) * buoy
                    ! final a-equation: (c1 + c2 + d4*a{w} - {detB}) = rhot*a *dz*updraft_rescale/dt + rhot*a*{w}
                    
                    coef_2(:,:) = b1* d3(:,:) *a(:,:,ilev, iupd)    ! coef_2 >= 0 by definition
                    coef_1(:,:) = c1(:,:)+c2(:,:) - b1*w3(:,:)* d3(:,:) *a(:,:,ilev, iupd)
                    coef_0(:,:) = -(c1(:,:)*w1(:,:) + c2(:,:)*w2(:,:) + &
                                    b2*db(:,:)*a(:,:,ilev, iupd)*buoy(:,:,ilev, iupd))
                                                 
                    where (coef_0(:,:) >= 0.) ! Updraft stops .
                                              ! Should be prevailed by adjust_flag = 2 or 3.
                        w   (:,:,ilev, iupd) = 0.
                        a   (:,:,ilev, iupd) = atd(:,:)
                        entL(:,:,ilev, iupd) = 0.; detL(:,:,ilev, iupd) = 0.
                        entB(:,:,ilev, iupd) = 0.
                        detB(:,:,ilev, iupd) = c1d(:,:)
                    elsewhere
                        where (coef_2(:,:) == 0.) 
                            ! d3 = 0, no entrainment, degenerate equation.
                            ! This implies that c1+c2 > 0, thus coef_1 > 0, coef_0 < 0 ==> w > 0.
                            w(:,:,ilev,iupd) = - coef_0(:,:) / coef_1(:,:) 
                        elsewhere  
                            ! Now coef_2 > 0 and coef_0 < 0 ==> w > 0
                            w(:,:,ilev,iupd) = (-coef_1(:,:) + &
                                sqrt(coef_1(:,:)**2. - 4.*coef_2(:,:)*coef_0(:,:)))/(2.*coef_2(:,:))
                        end where 
                        entL(:,:,ilev, iupd) =  d3(:,:)          * a(:,:,ilev, iupd) * w(:,:,ilev, iupd)
                        detL(:,:,ilev, iupd) = (d3(:,:)-d4(:,:)) * a(:,:,ilev, iupd) * w(:,:,ilev, iupd)
                        entB(:,:,ilev, iupd) = 0.
                        detB(:,:,ilev, iupd) = c1(:,:)+c2(:,:) + d4(:,:) * a(:,:,ilev, iupd) * w(:,:,ilev, iupd) &
                            - rhot(:,:,ilev)*a(:,:,ilev, iupd)*dz(:,:)*updraft_rescale/dt &
                            - rhot(:,:,ilev)*a(:,:,ilev, iupd) * w(:,:,ilev, iupd);
                    end where ! (coef_0(:,:) >= 0) 
                end where ! (adjust_flag(:,:) == 1) 
             
            endif ! (ilev == nlev)
            m(:,:,ilev, iupd) = rhot(:,:,ilev)*a(:,:,ilev,iupd)*w(:,:,ilev, iupd)
            
        enddo ! (ilev = nlev, 1, -1)
    enddo ! (iupd = 1, nupd)
            
end subroutine compute_au_wu


! Added ZTAN: 10/09/2017 -- New Pressure
! compute_au_wu: compute updated updraft fraction au and updraft velocity wu, with other auxiliary variables
subroutine compute_au_wu_newp (a_old, w_old, w_env, rhot, buoy, entr, detr, &
                               pfull, zfull, zhalf, dt, nlev, nupd, &
                               a, w, m, entL, entB, detL, detB)

    real   , intent(in),    dimension(:,:,:,:) :: a_old, w_old, buoy, entr, detr
    real   , intent(in),    dimension(:,:,:)   :: w_env, rhot, pfull, zfull, zhalf
    real   , intent(in)                        :: dt
    integer, intent(in)                        :: nlev, nupd
    real   , intent(out),   dimension(:,:,:,:) :: a, w, m, entL, entB, detL, detB
    
    real   , dimension(size(a_old,1),size(a_old,2))  :: dz, c1, c2, d3, d4, db, w1, w2, w3, &
                                                        w_adv, coef_0, coef_1, coef_2, a_denom, &
                                                        c1s, c1d, atd, f_1, f_2
    integer, dimension(size(a_old,1),size(a_old,2))  :: adjust_flag
    
    ! c1: old mass; c2: mass from below; d3: (entrained mass)/wu; d4: (entrained - detrained mass)/wu.
    ! w1: old wu; w2: wu from below; w3: w of entrained mass.
    ! w_adv: wu without buoyancy effect, (coef_0, coef_1, coef_2): coefficients of the quadratic equation.
    ! adjust_flag: flag for updraft fraction limiter.
    
    
    integer :: iupd, ilev
    real :: alpha_b  = 1./3.
    real :: beta_b
    real :: alpha_d  = 0.075
    real :: radius_d = 1500.
    real   , dimension(size(a_old,1),size(a_old,2))  :: beta_d
    
    real :: limit_fac = 100.
    
    beta_b = 1 - alpha_b
    
    a = 0.; w = 0.; m = 0.; entL = 0.; entB = 0.; detL = 0.; detB = 0.
    dz = 0.; c1 = 0.; c2 = 0.; d3 = 0.; d4 = 0.; db = 0.; w1 = 0.; w2 = 0.; w3 = 0.
    w_adv = 0.; coef_0 = 0.; coef_1 = 0.; coef_2 = 0.; a_denom = 0.
    c1s = 0.; c1d = 0.; atd = 0.; f_1 = 0.; f_2 = 0.
    
    do iupd = 1, nupd
        do ilev = nlev, 1, -1
            dz(:,:) = zhalf(:,:,ilev) - zhalf(:,:,ilev+1)
            c1(:,:) = rhot(:,:,ilev)*a_old(:,:,ilev,iupd)*dz(:,:)*updraft_rescale/dt
            w1(:,:) = w_old(:,:,ilev,iupd)
            w3(:,:) = w_env(:,:,ilev)
            
            if (ilev /= nlev) then
                w2(:,:) = w(:,:,ilev+1,iupd)
                c2(:,:) = rhot(:,:,ilev+1)*a(:,:,ilev+1, iupd)*w(:,:,ilev+1, iupd)
                
                ! Compute entrained and detrained mass coefficients
                if (au_optL == 1) then  ! use exponential formula (119)
                    ! Note: Substituting (entr-detr)*dz by exp[(entr-detr)*dz]-1.
                    ! This may be unstable if (entr-detr)*dz is not small (e.g., near SFC). 
                    where (entr(:,:,ilev,iupd)-detr(:,:,ilev, iupd) == 0.0)
                        d4(:,:) = 0.
                        d3(:,:) = rhot(:,:,ilev)*entr(:,:,ilev, iupd)*dz(:,:)
                    elsewhere
                        d4(:,:) = rhot(:,:,ilev)*(exp((entr(:,:,ilev, iupd)-detr(:,:,ilev, iupd))*dz(:,:))-1.)
                        d3(:,:) = entr(:,:,ilev, iupd)/(entr(:,:,ilev, iupd)-detr(:,:,ilev, iupd))*d4(:,:) 
                    end where
                else  ! use simple formula
                    d4(:,:) = rhot(:,:,ilev)*(entr(:,:,ilev,iupd)-detr(:,:,ilev, iupd))*dz(:,:)
                    d3(:,:) = rhot(:,:,ilev)* entr(:,:,ilev,iupd)*dz(:,:)
                endif
                
                where (d4(:,:) > rhot(:,:,ilev))  
                    ! MAY NEED TO LIMIT d4 < rhot because of implicit scheme.
                    d3(:,:) = d3(:,:)/d4(:,:)* rhot(:,:,ilev) ! now d3 = 
                    d4(:,:) = rhot(:,:,ilev)
                end where
                where (d4(:,:) < - limit_fac * rhot(:,:,ilev)) 
                    ! LIMIT d4 > - limit_fac * rhot? May be unnecessary... ZTAN 09/26/2017
                    d3(:,:) = d3(:,:)/d4(:,:)*(- limit_fac * rhot(:,:,ilev))
                    d4(:,:) = - limit_fac * rhot(:,:,ilev)
                end where
            else
                c2(:,:) = 0.; w2(:,:) = 0.; d3(:,:) = 0.; d4(:,:) = 0.
            endif  ! (ilev /= nlev)
            
            ! write(*,*), ilev, maxval(d3(:,:)), minval(d3(:,:)), maxval(d4(:,:)), minval(d4(:,:))
            db(:,:) = rhot(:,:,ilev)*dz(:,:)
            beta_d(:,:) = alpha_d/radius_d/sqrt(max(a_old(:,:,ilev,iupd), 1.e-3))
            
            ! ------ Step 1 ------
            ! Assume entB = 0 and detB = 0, and calculate {w},{a},{entL},{detL}.
            ! w-equation: (c1 + c2 + d3*{a}*{w}) * {w} = ...
            !             (c1 * w1) + (c2 * w2) + (d3*{a}*{w} * w3) + ...
            !             (db * beta_b) *{a}* buoy - db * {a} * beta_d * ({w} - w3)^2
            ! a-equation: (c1 + c2 + d4*{a}*{w}) = ...
            !              rhot*{a} *dz*updraft_rescale/dt + rhot*{a}*{w}
            ! Solve {w} first, then {a}, then compute entL and detL if {a} satisfies condition.
            
            where (c1(:,:) + c2(:,:) <= 0.)  !   no prev updraft and no incoming updraft
                ! ==> no current updraft without massive entrainment. a_old should be zero here.
                w(:,:,ilev,iupd) = 0.; a(:,:,ilev,iupd) = a_old(:,:,ilev,iupd)
                entL(:,:,ilev,iupd) = 0.; entB(:,:,ilev,iupd) = 0.
                detL(:,:,ilev,iupd) = 0.; detB(:,:,ilev,iupd) = 0.
            elsewhere
                w_adv (:,:) = (c1(:,:)*w1(:,:)+c2(:,:)*w2(:,:))/(c1(:,:)+c2(:,:))
                coef_2(:,:) = d3(:,:) - d4(:,:) + rhot(:,:,ilev) + db(:,:) * beta_d(:,:)
                ! > 0 because d3 >= 0 and we limit d4 <= rhot
                coef_1(:,:) = rhot(:,:,ilev)*dz(:,:)*updraft_rescale/dt &
                              - d3(:,:)*w3(:,:) - (rhot(:,:,ilev)-d4(:,:))* w_adv(:,:) &
                              - 2. * db(:,:) * beta_d(:,:) * w3(:,:)
                              
                coef_0(:,:) = - w_adv(:,:) * rhot(:,:,ilev)*dz(:,:)*updraft_rescale/dt & 
                              - db(:,:) * beta_b * buoy(:,:,ilev, iupd) &
                              + db(:,:) * beta_d(:,:) * w3(:,:) * w3(:,:)
                              
                ! where coef_2 <= 0 <- removed with the above discussion
                
                where (coef_0(:,:) >= 0) 
                ! The non-diluted w_nent is negative, i.e., buoy << 0. 
                ! Thus updraft stops, boundary detrainment occurs (Step 3).
                    w   (:,:,ilev, iupd) = 0.; a(:,:,ilev, iupd) = a_old(:,:,ilev,iupd) ! a will be adjusted by Step 3
                    entL(:,:,ilev, iupd) = 0.; detL(:,:,ilev, iupd) = 0.
                    entB(:,:,ilev, iupd) = 0.; detB(:,:,ilev, iupd) = 0.
                elsewhere
                    w (:,:,ilev, iupd) = (-coef_1(:,:) + &
                        sqrt(coef_1(:,:)**2. - 4.*coef_2(:,:)*coef_0(:,:)))/(2.*coef_2(:,:))
                    a_denom (:,:) =  rhot(:,:,ilev)*dz(:,:)*updraft_rescale/dt + &
                                    (rhot(:,:,ilev)-d4(:,:))*w(:,:,ilev, iupd)
                    where (a_denom(:,:) < 0.)
                        a(:,:,ilev, iupd) = 1.e5 ! Temporarily filling in a large number,
                                                 ! Will be fixed by increased detrainment (detB).
                    elsewhere
                        a(:,:,ilev, iupd) = (c1(:,:) + c2(:,:))/a_denom (:,:)
                        ! At the limit of large entrainment (e.g., b is large but w is small)
                        ! -> a_denom = rho*dz/dt, c1 = rho*a_o*dz/dt, c2 = rho- * a- * w-
                        ! -> a = a_o + ( rho- /rho)*(w- /(dz/dt))* a- . If dz/dt << w-, we will have a >> a_o
                        ! This is expected !
                        ! And: d3 = d4 = rho => 
                        ! coef_2 = (b1-1)*rho, coef_1 = rho*dz/dt - b1*rho*w3, 
                        ! coef_0 = - w_adv * rho *dz/dt - b2*rho*dz*buoy
                        ! At the limit of dz/dt -> 0: 
                        ! coef_2 = (b1-1)*rho, coef_1 = -b1*rho*w3 ~ 0, coef_0 = -b2*rho*dz*buoy
                        ! => w = sqrt(b2*dz*buoy / (b1-1)) ~ sqrt(1.33*dz*buoy)
                        ! The dominant balance is: 
                        
                        ! However, if entrainment is slightly smaller (d4 = 0.9*rho), 
                        ! then a_denom = rho*(dz/dt +0.1w)
                        ! a = (a_o*rho*dz/dt + a- * rho- * w-)/[rho*(dz/dt +0.1w)]
                        ! If dz/dt << w-, w => a ~ a- * 10(w-/w) is still huge...
                        
                    end where
                    entB(:,:,ilev, iupd) = 0.; detB(:,:,ilev, iupd) = 0.
                    entL(:,:,ilev, iupd) =  d3(:,:)          * a(:,:,ilev, iupd) * w(:,:,ilev, iupd)
                    detL(:,:,ilev, iupd) = (d3(:,:)-d4(:,:)) * a(:,:,ilev, iupd) * w(:,:,ilev, iupd)
                end where ! (coef_0(:,:) >= 0) 
            end where ! (c1(:,:) + c2(:,:) <= 0.)


            ! write(*,*), ilev, maxval(entL(:,:)), minval(entL(:,:)), maxval(detL(:,:)), minval(detL(:,:))
            
            if (ilev == nlev) then
            ! ------ Step 2 ------
            ! Fix SFC fraction by:
            ! (1) If buoy > 0 => w > 0, aN < a0 because entL = 0 => entB > 0.
            ! (2) If buoy <= 0. Swap airmass, use detB = c1, calculate entB accordingly.
                ! if (au_optB_srf == 1) then ! Disabled ZTAN 09/26/2017
                !     c1s(:,:) = c1(:,:)
                ! else
                !     c1s(:,:) = 0.0
                ! endif
                
                ! where (w (:,:,ilev, iupd) < 1.0e-6)
                !     a(:,:,ilev, iupd) = a_old(:,:,ilev, iupd)
                !     w(:,:,ilev, iupd) = 0.0
                !     entL(:,:,ilev, iupd) = 0.0; detL(:,:,ilev, iupd) = 0.0
                !     entB(:,:,ilev, iupd) = rhot(:,:,ilev)*a(:,:,ilev, iupd)*dz(:,:)*updraft_rescale/dt - c2(:,:)
                !     detB(:,:,ilev, iupd) = c1(:,:)
                    
                ! elsewhere
                ! Calculate entB for target a(ilev, iupd) = a_o(ilev, iupd)
                ! final w-equation: (c1 + c2 + d3*a*{w} + {entB}) * {w} = ...
                !                   (c1 * w1) + (c2 * w2) + (d3*a*{w} + {entB}) * w3 + ...
                !                   (db * beta_b * a) * buoy - db * a * beta_d * ({w} - w3)^2
                ! final a-equation: (c1 + c2 + d4*a{w} + {entB} - {detB}) = ...
                !                    rhot*a *dz*updraft_rescale/dt + rhot*a*{w}
                !                  -> f_2 * {w} + f_1 = {entB}
                ! 
                ! SFC w-equation at steady state: {w} = w1, w2 = w3 = 0, au = const; 
                !     also: c1 = rho*a*dz/dt, c2 = 0, d3 = d4 = 0, db = rho*dz
                ! final w-equation: {w} = (rho * dz * b2 * {a}) * buoy / ({entB} * b1) 
                ! (1) if the updraft is always swapped -> entB = c1 + rho*{a}*{w}; detB = c1
                !     -> this gives (rho * b1 * {a}) * {w}^2 + c1 * b1 * {w} = (rho * dz * b2 * {a}) * buoy
                !     -> i.e., b1*{w}*(dz/dt+{w}) = dz*b2*buoy (THIS IS PROBLEMATIC, especially if dt is small)
                ! (2) if the updraft is NOT swapped -> entB = rho*a*w; detB = 0.
                !     -> this gives (rho * b1 * {a}) * {w}^2 = (rho * dz * b2 * {a}) * buoy
                !     -> i.e., b1*{w}^2 = dz*b2*buoy
                
                detB(:,:,ilev, iupd) = 0.0 ! c1s(:,:)  ! Changed ZTAN 09/26/2017: Do not swap 
                    
                a(:,:,ilev, iupd) = a_old(:,:,ilev, iupd)
                f_1(:,:) = rhot(:,:,ilev)*a(:,:,ilev, iupd)*dz(:,:)*updraft_rescale/dt - &
                      (c1(:,:)+c2(:,:) - detB(:,:,ilev, iupd))
                f_2(:,:) = (rhot(:,:,ilev) - d4(:,:))*a(:,:,ilev, iupd)

                ! Change below
                coef_2(:,:) = f_2(:,:) + d3(:,:) *a(:,:,ilev, iupd) + db(:,:)*a(:,:,ilev, iupd)* beta_d(:,:)  ! coef_2 > 0 by definition
                coef_1(:,:) = c1(:,:)+c2(:,:) + f_1(:,:) - w3(:,:)*(f_2(:,:) + d3(:,:) *a(:,:,ilev, iupd)) &
                              - 2. * db(:,:) * beta_d(:,:) * a(:,:,ilev, iupd) * w3(:,:)
                coef_0(:,:) = -(c1(:,:)*w1(:,:) + c2(:,:)*w2(:,:) + f_1(:,:)*w3(:,:) &
                                + db(:,:) * beta_b * a(:,:,ilev, iupd)*buoy(:,:,ilev, iupd) &
                                - db(:,:) * beta_d(:,:) * a(:,:,ilev, iupd) * w3(:,:) * w3(:,:) )
                               
                where (coef_0(:,:) >= 0.)  ! Updraft at SFC is not buoyant
                    w(:,:,ilev, iupd) = 0.0
                    entL(:,:,ilev, iupd) = 0.0; detL(:,:,ilev, iupd) = 0.0
                    entB(:,:,ilev, iupd) = rhot(:,:,ilev)*a(:,:,ilev, iupd)*dz(:,:)*updraft_rescale/dt - c2(:,:)
                    detB(:,:,ilev, iupd) = c1(:,:)
                elsewhere
                    w (:,:,ilev, iupd) = (-coef_1(:,:) + &
                        sqrt(coef_1(:,:)**2. - 4.*coef_2(:,:)*coef_0(:,:)))/(2.*coef_2(:,:))
                    entB(:,:,ilev, iupd) = f_1(:,:) + f_2(:,:) * w(:,:,ilev, iupd)
                    entL(:,:,ilev, iupd) =  d3(:,:)          * a(:,:,ilev, iupd) * w(:,:,ilev, iupd)
                    detL(:,:,ilev, iupd) = (d3(:,:)-d4(:,:)) * a(:,:,ilev, iupd) * w(:,:,ilev, iupd)
                end where ! (coef_0 >= 0.)
                    
                ! end where ! (w (:,:,ilev, iupd) < 1.0e-6)
                
                ! Adjust entB and detB for scalars -- switch updraft thl/qt with environment buoyant tail
                ! ZTAN: 09/25/2017: This is disabled due to inconsistency. The 'switch' is now done at the beginning.
                ! detB(:,:,ilev, iupd) = detB(:,:,ilev, iupd) + c1(:,:)*99. ! 99
                ! entB(:,:,ilev, iupd) = entB(:,:,ilev, iupd) + c1(:,:)*99. ! 99
            else
            ! ------ Step 3 ------
            ! Adjust area fraction at other levels by tuning up detL.
            ! Note: this algorithm is applicable to adjusting down only.
                
                ! Limitation on w_u
                adjust_flag(:,:) = 0
                if (au_optB_wu == 1) then ! detraining parcels with zero velocity
                    where ( (w(:,:,ilev, iupd) <= 1.0e-6) .or. (a(:,:,ilev, iupd) <= 1.0e-6) )
                        adjust_flag(:,:) = 2  ! Kill the updraft
                    end where
                elseif (au_optB_wu == 11) then ! detraining parcels with zero velocity or zero area
                    where ( (a(:,:,ilev, iupd) <= 1.0e-6))
                        adjust_flag(:,:) = 2  ! Kill the updraft
                    elsewhere
                        where (w(:,:,ilev, iupd) <= 1.0e-6)
                            adjust_flag(:,:) = 3  ! Kill the updraft gradually
                        end where
                    end where
                endif ! (au_optB_wu == 1) (au_optB_wu == 11)
            
            
                ! Limitation on a_u is top priority -- will overwrite adjust_flag due to w_u.
                if (au_optB == 1) then    ! constrain updraft fraction < SFC fraction
                    where (a(:,:,ilev, iupd)> a(:,:,nlev, iupd))
                        a(:,:,ilev,iupd) = a(:,:,nlev, iupd);  adjust_flag(:,:) = 1
                    end where
                elseif (au_optB == 2) then    ! constrain updraft fraction < 2*SFC fraction
                    where (a(:,:,ilev, iupd)> 2.0*a(:,:,nlev, iupd))
                        a(:,:,ilev,iupd) = 2.0*a(:,:,nlev, iupd);  adjust_flag(:,:) = 1
                    end where
                elseif (au_optB == 100) then  ! constrain updraft fraction < 0.5/nupd
                    where (a(:,:,ilev, iupd)> .5/nupd)
                        a(:,:,ilev,iupd) = .5/nupd;  adjust_flag(:,:) = 1
                    end where
                endif ! (au_optB == 1, 2, 100)
                
                where (adjust_flag(:,:) == 3) ! Adjust w to zero directly, 
                                              ! but detrain a over a 15 minute timescale.
                    w   (:,:,ilev, iupd) = 0. 
                    a   (:,:,ilev, iupd) = a_old(:,:,ilev,iupd) * exp(-dt/900.0)   
                    entL(:,:,ilev, iupd) = 0.; detL(:,:,ilev, iupd) = 0.
                    entB(:,:,ilev, iupd) = 0. 
                    detB(:,:,ilev, iupd) = c1(:,:)+c2(:,:) - &
                         rhot(:,:,ilev)*a(:,:,ilev, iupd)*dz(:,:)*updraft_rescale/dt
                end where ! (adjust_flag(:,:) == 3) 
                
                where (adjust_flag(:,:) == 2) ! Adjust to zero directly.
                    w   (:,:,ilev, iupd) = 0.  
                    a   (:,:,ilev, iupd) = 0.  
                    entL(:,:,ilev, iupd) = 0.; detL(:,:,ilev, iupd) = 0.
                    entB(:,:,ilev, iupd) = 0.
                    detB(:,:,ilev, iupd) = c1(:,:)+c2(:,:)
                end where ! (adjust_flag(:,:) == 2) 
                
                ! c1d is only used where adjust_flag == 1
                if (au_optB_wu == 1) then   ! Detrainment everything. atd is target a.
                    c1d(:,:) = c1(:,:)+c2(:,:)
                    atd(:,:) = 0.
                elseif (au_optB_wu == 11) then ! Detrainment to target a (= atd).
                    c1d(:,:) = c1(:,:)+c2(:,:) &
                               - rhot(:,:,ilev)*a(:,:,ilev, iupd)*dz(:,:)*updraft_rescale/dt
                    atd(:,:) = a(:,:,ilev, iupd)
                endif
                
                where (adjust_flag(:,:) == 1)
                    ! Calculate detB for target a(ilev, iupd) = a_o(ilev, iupd), with entB = 0
                    ! final w-equation: (c1 + c2 + d3*a*{w}) * {w} = ...
                    !                   (c1 * w1) + (c2 * w2) + (d3*a*{w} * w3) + ...
                    !                   (db * beta_b) *a* buoy - db * a * beta_d * ({w} - w3)^2
                    ! final a-equation: (c1 + c2 + d4*a*{w} - {detB}) = ...
                    !                   rhot*a *dz*updraft_rescale/dt + rhot*a*{w}
                    coef_2(:,:) = d3(:,:) *a(:,:,ilev, iupd) &
                                  + db(:,:) * beta_d(:,:) *a(:,:,ilev, iupd)  ! coef_2 > 0 by definition
                    coef_1(:,:) = c1(:,:)+c2(:,:) - w3(:,:)* d3(:,:) *a(:,:,ilev, iupd) &
                                  - 2. * db(:,:) * beta_d(:,:) *a(:,:,ilev, iupd) * w3(:,:)
                    coef_0(:,:) = -(c1(:,:)*w1(:,:) + c2(:,:)*w2(:,:) &
                                    + db(:,:) * beta_b * a(:,:,ilev, iupd)*buoy(:,:,ilev, iupd) &
                                    - db(:,:) * beta_d(:,:)* a(:,:,ilev, iupd) * w3(:,:) * w3(:,:))
                                    
                    !coef_2(:,:) = b1* d3(:,:) *a(:,:,ilev, iupd)    ! coef_2 >= 0 by definition
                    !coef_1(:,:) = c1(:,:)+c2(:,:) - b1*w3(:,:)* d3(:,:) *a(:,:,ilev, iupd)
                    !coef_0(:,:) = -(c1(:,:)*w1(:,:) + c2(:,:)*w2(:,:) + &
                    !                b2*db(:,:)*a(:,:,ilev, iupd)*buoy(:,:,ilev, iupd))
                                                 
                    where (coef_0(:,:) >= 0.) ! Updraft stops .
                                              ! Should be prevailed by adjust_flag = 2 or 3.
                        w   (:,:,ilev, iupd) = 0.
                        a   (:,:,ilev, iupd) = atd(:,:)
                        entL(:,:,ilev, iupd) = 0.; detL(:,:,ilev, iupd) = 0.
                        entB(:,:,ilev, iupd) = 0.
                        detB(:,:,ilev, iupd) = c1d(:,:)
                    elsewhere
                        where (coef_2(:,:) == 0.) 
                            ! d3 = 0, no entrainment, degenerate equation.
                            ! This implies that c1+c2 > 0, thus coef_1 > 0, coef_0 < 0 ==> w > 0.
                            w(:,:,ilev,iupd) = - coef_0(:,:) / coef_1(:,:) 
                        elsewhere  
                            ! Now coef_2 > 0 and coef_0 < 0 ==> w > 0
                            w(:,:,ilev,iupd) = (-coef_1(:,:) + &
                                sqrt(coef_1(:,:)**2. - 4.*coef_2(:,:)*coef_0(:,:)))/(2.*coef_2(:,:))
                        end where 
                        entL(:,:,ilev, iupd) =  d3(:,:)          * a(:,:,ilev, iupd) * w(:,:,ilev, iupd)
                        detL(:,:,ilev, iupd) = (d3(:,:)-d4(:,:)) * a(:,:,ilev, iupd) * w(:,:,ilev, iupd)
                        entB(:,:,ilev, iupd) = 0.
                        detB(:,:,ilev, iupd) = c1(:,:)+c2(:,:) + d4(:,:) * a(:,:,ilev, iupd) * w(:,:,ilev, iupd) &
                            - rhot(:,:,ilev)*a(:,:,ilev, iupd)*dz(:,:)*updraft_rescale/dt &
                            - rhot(:,:,ilev)*a(:,:,ilev, iupd) * w(:,:,ilev, iupd);
                    end where ! (coef_0(:,:) >= 0) 
                end where ! (adjust_flag(:,:) == 1) 
             
            endif ! (ilev == nlev)
            m(:,:,ilev, iupd) = rhot(:,:,ilev)*a(:,:,ilev,iupd)*w(:,:,ilev, iupd)
            
        enddo ! (ilev = nlev, 1, -1)
    enddo ! (iupd = 1, nupd)
            
end subroutine compute_au_wu_newp


! Added ZTAN: 09/28/2017  --  Use explicit entrainment (but still implicit detrainment)
! compute_au_wu: compute updated updraft fraction au and updraft velocity wu, with other auxiliary variables
subroutine compute_au_wu_expentr (a_old, w_old, w_env, rhot, buoy, entr, detr, &
                          pfull, zfull, zhalf, dt, nlev, nupd, &
                          a, w, m, entL, entB, detL, detB, use_lower)

    real   , intent(in),    dimension(:,:,:,:) :: a_old, w_old, buoy, entr, detr
    real   , intent(in),    dimension(:,:,:)   :: w_env, rhot, pfull, zfull, zhalf
    real   , intent(in)                        :: dt
    integer, intent(in)                        :: nlev, nupd
    logical, intent(in)                        :: use_lower
    real   , intent(out),   dimension(:,:,:,:) :: a, w, m, entL, entB, detL, detB
    
    real   , dimension(size(a_old,1),size(a_old,2))  :: dz, c1, c2, c3, d4, db, w1, w2, w3, &
                                                        w_adv, b_adv, coef_0, coef_1, coef_2, a_denom, &
                                                        c1s, c1d, atd, f_1, f_2
    integer, dimension(size(a_old,1),size(a_old,2))  :: adjust_flag
    
    ! c1: old mass; c2: mass from below; c3: (entrained mass); d4: (detrained mass)/wu.
    ! NOTE: d4 is of opposite sign against the compute_au_wu subroutine.
    ! w1: old wu; w2: wu from below; w3: w of entrained mass.
    ! w_adv: wu without buoyancy effect, (coef_0, coef_1, coef_2): coefficients of the quadratic equation.
    ! adjust_flag: flag for updraft fraction limiter.
    
    
    integer :: iupd, ilev
    real :: b1  = 2.5
    real :: b2  = 2./3. ! Soares 04: 2.0
    real :: limit_fac = 100000.
    ! logical :: use_lower = .true.
    logical :: limit_entr = .true.

    a = 0.; w = 0.; m = 0.; entL = 0.; entB = 0.; detL = 0.; detB = 0.
    dz = 0.; c1 = 0.; c2 = 0.; c3 = 0.; d4 = 0.; db = 0.; w1 = 0.; w2 = 0.; w3 = 0.
    w_adv = 0.; b_adv = 0.; coef_0 = 0.; coef_1 = 0.; coef_2 = 0.; a_denom = 0.
    c1s = 0.; c1d = 0.; atd = 0.; f_1 = 0.; f_2 = 0.
    
    do iupd = 1, nupd
        do ilev = nlev, 1, -1
            dz(:,:) = zhalf(:,:,ilev) - zhalf(:,:,ilev+1)
            c1(:,:) = rhot(:,:,ilev)*a_old(:,:,ilev,iupd)*dz(:,:)*updraft_rescale/dt
            w1(:,:) = w_old(:,:,ilev,iupd)
            w3(:,:) = w_env(:,:,ilev)
            
            if (ilev /= nlev) then
                w2(:,:) = w(:,:,ilev+1,iupd)
                c2(:,:) = rhot(:,:,ilev+1)*a(:,:,ilev+1, iupd)*w(:,:,ilev+1, iupd)
                if (use_lower) then
                    c3(:,:) = rhot(:,:,ilev+1)*a(:,:,ilev+1, iupd)*w(:,:,ilev+1, iupd)*min(entr(:,:,ilev,iupd)*dz(:,:), limit_fac)
                else ! use old MF
                    c3(:,:) = rhot(:,:,ilev)*a_old(:,:,ilev,iupd)*w_old(:,:,ilev,iupd)*min(entr(:,:,ilev,iupd)*dz(:,:), limit_fac)
                end if
                ! For diagnostic scheme: c3(:,:) = c2(:,:)*entr(:,:,ilev,iupd)*dz(:,:)
                
                ! Note 10/10/2017: For stability, entrainment should be smaller than environmental airmass
                if (limit_entr) then
                    c3(:,:) = min(c3(:,:), rhot(:,:,ilev)*(1-a_old(:,:,ilev,iupd)*updraft_rescale)*dz(:,:)/dt)
                    ! This implies: a*w*entr*dz < (1-a)dz/dt --> entr < (1-a)/a /(w*dt)
                    ! If: entr = 1/tau/w --> dt < (1-a)/a * tau :: This is satisfied...
                end if

                ! Compute detrained mass coefficients with simple formula
                d4(:,:) = rhot(:,:,ilev)*min(detr(:,:,ilev, iupd)*dz(:,:), limit_fac)
            else
                c2(:,:) = 0.; w2(:,:) = 0.; c3(:,:) = 0.; d4(:,:) = 0.
            endif  ! (ilev /= nlev)
            
            ! write(*,*), ilev, maxval(d3(:,:)), minval(d3(:,:)), maxval(d4(:,:)), minval(d4(:,:))
            db(:,:) = rhot(:,:,ilev)*dz(:,:)

            ! ------ Step 1 ------
            ! Assume entB = 0 and detB = 0, and calculate {w},{a},{entL},{detL}.
            ! w-equation: (c1 + c2 + c3 * b1) * {w} = ...
            !             (c1 * w1) + (c2 * w2) + (c3 * b1 * w3) + ...
            !             (db * b2) *{a}* buoy
            ! a-equation: (c1 + c2 + c3 - d4*{a}*{w}) = ...
            !              rhot*{a} *dz*updraft_rescale/dt + rhot*{a}*{w}
            ! Solve {w} first, then {a}, then compute entL and detL if {a} satisfies condition.
            
            where (c1(:,:) + c2(:,:) + c3(:,:) <= 0.)  !   no prev updraft and no incoming updraft
                ! ==> no current updraft and no massive entrainment. a_old should be zero here.
                w(:,:,ilev,iupd) = 0.; a(:,:,ilev,iupd) = a_old(:,:,ilev,iupd)
                entL(:,:,ilev,iupd) = 0.; entB(:,:,ilev,iupd) = 0.
                detL(:,:,ilev,iupd) = 0.; detB(:,:,ilev,iupd) = 0.
            elsewhere
                w_adv (:,:) = (c1(:,:)*w1(:,:)+c2(:,:)*w2(:,:)+b1*c3(:,:)*w3(:,:))/(c1(:,:)+c2(:,:)+b1*c3(:,:))
                b_adv (:,:) = (c1(:,:)+c2(:,:)+c3(:,:))*buoy(:,:,ilev,iupd)/(c1(:,:)+c2(:,:)+b1*c3(:,:))
                coef_2(:,:) = d4(:,:) + rhot(:,:,ilev) ! > 0 because d4 >= 0
                coef_1(:,:) = rhot(:,:,ilev)*dz(:,:)*updraft_rescale/dt - &
                              (rhot(:,:,ilev)+d4(:,:))* w_adv(:,:)
                coef_0(:,:) = - w_adv(:,:) * rhot(:,:,ilev)*dz(:,:)*updraft_rescale/dt - &
                              b2*db(:,:)*b_adv(:,:)
                              
                where (coef_0(:,:) >= 0) 
                ! Both w_adv and buoy <= 0. Thus coef_1 >= 0, 
                ! Note that coef_2 > 0. Thus, both roots are non-positive.
                ! Updraft stops, boundary detrainment occurs (Step 3).
                    w   (:,:,ilev, iupd) = 0.; a(:,:,ilev, iupd) = a_old(:,:,ilev,iupd) ! a will be adjusted by Step 3
                    entL(:,:,ilev, iupd) = 0.; detL(:,:,ilev, iupd) = 0.
                    entB(:,:,ilev, iupd) = 0.; detB(:,:,ilev, iupd) = 0.
                elsewhere
                    w (:,:,ilev, iupd) = (-coef_1(:,:) + &
                        sqrt(coef_1(:,:)**2. - 4.*coef_2(:,:)*coef_0(:,:)))/(2.*coef_2(:,:))
                    
                    a_denom (:,:) =  rhot(:,:,ilev)*dz(:,:)*updraft_rescale/dt + &
                                    (rhot(:,:,ilev)+d4(:,:))*w(:,:,ilev, iupd)
                    ! This denominator is always positive.
                    
                    a(:,:,ilev, iupd) = (c1(:,:) + c2(:,:) + c3(:,:))/a_denom (:,:)
                    
                    entB(:,:,ilev, iupd) = 0.; detB(:,:,ilev, iupd) = 0.
                    entL(:,:,ilev, iupd) = c3(:,:) 
                    detL(:,:,ilev, iupd) = d4(:,:) * a(:,:,ilev, iupd) * w(:,:,ilev, iupd)
                end where ! (coef_0(:,:) >= 0) 
            end where ! (c1(:,:) + c2(:,:) <= 0.)


            ! write(*,*), ilev, maxval(entL(:,:)), minval(entL(:,:)), maxval(detL(:,:)), minval(detL(:,:))
            
            if (ilev == nlev) then
            ! ------ Step 2 ------
            ! Fix SFC fraction by:
            ! (1) If buoy > 0 => w > 0, aN < a0 because entL = 0 => entB > 0.
            ! (2) If buoy <= 0. Swap airmass, use detB = c1, calculate entB accordingly.

                ! Calculate entB for target a(ilev, iupd) = a_o(ilev, iupd)
                ! final w-equation: (c1 + c2 + c3 * b1 + {entB} * b1) * {w} = ...
                !                   (c1 * w1) + (c2 * w2) + (c3 + {entB}) * b1 * w3 + (db * b2 * a) * buoy
                ! final a-equation: (c1 + c2 + c3 - d4*a*{w} + {entB} - {detB}) = ...
                !                    rhot*a *dz*updraft_rescale/dt + rhot*a*{w}
                !                  -> f_2 * {w} + f_1 = {entB}
                
                detB(:,:,ilev, iupd) = 0.0 ! c1s(:,:)  ! Changed ZTAN 09/26/2017: Do not swap 
                    
                a(:,:,ilev, iupd) = a_old(:,:,ilev, iupd)
                f_1(:,:) = rhot(:,:,ilev)*a(:,:,ilev, iupd)*dz(:,:)*updraft_rescale/dt - &
                      (c1(:,:)+c2(:,:)+c3(:,:) - detB(:,:,ilev, iupd))
                f_2(:,:) = (rhot(:,:,ilev) + d4(:,:))*a(:,:,ilev, iupd)

                coef_2(:,:) = b1*f_2(:,:)   ! coef_2 > 0 by definition
                coef_1(:,:) = c1(:,:)+c2(:,:)+ b1*c3(:,:) + b1*f_1(:,:) - b1*w3(:,:)*f_2(:,:)
                coef_0(:,:) = -(c1(:,:)*w1(:,:) + c2(:,:)*w2(:,:) + b1*c3(:,:)*w3(:,:) + b1*f_1(:,:)*w3(:,:) + &
                                    b2*db(:,:)*a(:,:,ilev, iupd)*buoy(:,:,ilev, iupd))
                               
                where (coef_0(:,:) >= 0.)  ! Updraft at SFC is not buoyant
                    w(:,:,ilev, iupd) = 0.0
                    entL(:,:,ilev, iupd) = 0.0; detL(:,:,ilev, iupd) = 0.0
                    entB(:,:,ilev, iupd) = rhot(:,:,ilev)*a(:,:,ilev, iupd)*dz(:,:)*updraft_rescale/dt - c2(:,:) - c3(:,:)
                    detB(:,:,ilev, iupd) = c1(:,:)
                elsewhere
                    w (:,:,ilev, iupd) = (-coef_1(:,:) + &
                        sqrt(coef_1(:,:)**2. - 4.*coef_2(:,:)*coef_0(:,:)))/(2.*coef_2(:,:))
                    entB(:,:,ilev, iupd) = f_1(:,:) + f_2(:,:) * w(:,:,ilev, iupd)
                    entL(:,:,ilev, iupd) = c3(:,:) 
                    detL(:,:,ilev, iupd) = d4(:,:) * a(:,:,ilev, iupd) * w(:,:,ilev, iupd)
                end where ! (coef_0 >= 0.)
                
                ! Adjust entB and detB for scalars -- switch updraft thl/qt with environment buoyant tail
                ! ZTAN: 09/25/2017: This is disabled due to inconsistency. The 'switch' is now done at the beginning.
                ! detB(:,:,ilev, iupd) = detB(:,:,ilev, iupd) + c1(:,:)*99. ! 99
                ! entB(:,:,ilev, iupd) = entB(:,:,ilev, iupd) + c1(:,:)*99. ! 99
            else
            ! ------ Step 3 ------
            ! Adjust area fraction at other levels by tuning up detL.
            ! Note: this algorithm is applicable to adjusting down only.
                
                ! Limitation on w_u
                adjust_flag(:,:) = 0
                if (au_optB_wu == 1) then ! detraining parcels with zero velocity
                    where ( (w(:,:,ilev, iupd) <= 1.0e-6) .or. (a(:,:,ilev, iupd) <= 1.0e-6) )
                        adjust_flag(:,:) = 2  ! Kill the updraft
                    end where
                elseif (au_optB_wu == 11) then ! detraining parcels with zero velocity or zero area
                    where ( (a(:,:,ilev, iupd) <= 1.0e-6))
                        adjust_flag(:,:) = 2  ! Kill the updraft
                    elsewhere
                        where (w(:,:,ilev, iupd) <= 1.0e-6)
                            adjust_flag(:,:) = 3  ! Kill the updraft gradually
                        end where
                    end where
                endif ! (au_optB_wu == 1) (au_optB_wu == 11)
            
            
                ! Limitation on a_u is top priority -- will overwrite adjust_flag due to w_u.
                if (au_optB == 1) then    ! constrain updraft fraction < SFC fraction
                    where (a(:,:,ilev, iupd)> a(:,:,nlev, iupd))
                        a(:,:,ilev,iupd) = a(:,:,nlev, iupd);  adjust_flag(:,:) = 1
                    end where
                elseif (au_optB == 2) then    ! constrain updraft fraction < 2*SFC fraction
                    where (a(:,:,ilev, iupd)> 2.0*a(:,:,nlev, iupd))
                        a(:,:,ilev,iupd) = 2.0*a(:,:,nlev, iupd);  adjust_flag(:,:) = 1
                    end where
                elseif (au_optB == 100) then  ! constrain updraft fraction < 0.5/nupd
                    where (a(:,:,ilev, iupd)> .5/nupd)
                        a(:,:,ilev,iupd) = .5/nupd;  adjust_flag(:,:) = 1
                    end where
                endif ! (au_optB == 1, 2, 100)
                
                where (adjust_flag(:,:) == 3) ! Adjust w to zero directly, 
                                              ! but detrain a over a 15 minute timescale.
                    w   (:,:,ilev, iupd) = 0. 
                    a   (:,:,ilev, iupd) = a_old(:,:,ilev,iupd) * exp(-dt/900.0)   
                    entL(:,:,ilev, iupd) = c3(:,:); detL(:,:,ilev, iupd) = 0.
                    entB(:,:,ilev, iupd) = 0. 
                    detB(:,:,ilev, iupd) = c1(:,:)+c2(:,:)+c3(:,:) - &
                         rhot(:,:,ilev)*a(:,:,ilev, iupd)*dz(:,:)*updraft_rescale/dt
                end where ! (adjust_flag(:,:) == 3) 
                
                where (adjust_flag(:,:) == 2) ! Adjust to zero directly.
                    w   (:,:,ilev, iupd) = 0.  
                    a   (:,:,ilev, iupd) = 0.  
                    entL(:,:,ilev, iupd) = c3(:,:); detL(:,:,ilev, iupd) = 0.
                    entB(:,:,ilev, iupd) = 0.
                    detB(:,:,ilev, iupd) = c1(:,:)+c2(:,:)+c3(:,:)
                end where ! (adjust_flag(:,:) == 2) 
                
                ! c1d is only used where adjust_flag == 1
                if (au_optB_wu == 1) then   ! Detrainment everything. atd is target a.
                    c1d(:,:) = c1(:,:)+c2(:,:)+c3(:,:)
                    atd(:,:) = 0.
                elseif (au_optB_wu == 11) then ! Detrainment to target a (= atd).
                    c1d(:,:) = c1(:,:)+c2(:,:)+c3(:,:) &
                               - rhot(:,:,ilev)*a(:,:,ilev, iupd)*dz(:,:)*updraft_rescale/dt
                    atd(:,:) = a(:,:,ilev, iupd)
                endif
                
                where (adjust_flag(:,:) == 1)
                    ! Calculate detB for target a(ilev, iupd) = a_o(ilev, iupd), with entB = 0
                    ! final w-equation: (c1 + c2 + c3 * b1) * {w} = ...
                    !                   (c1 * w1) + (c2 * w2) + c3 * b1 * w3 + (db*a * b2) * buoy
                    ! final a-equation: (c1 + c2 + c3 - d4*a{w} - {detB}) = rhot*a *dz*updraft_rescale/dt + rhot*a*{w}
                    
                    coef_2(:,:) = 0.
                    coef_1(:,:) = c1(:,:)+c2(:,:) + b1*c3(:,:)
                    coef_0(:,:) = -(c1(:,:)*w1(:,:) + c2(:,:)*w2(:,:) + b1*c3(:,:)*w3(:,:) + &
                                    b2*db(:,:)*a(:,:,ilev, iupd)*buoy(:,:,ilev, iupd))
                                                 
                    where (coef_0(:,:) >= 0.) ! Updraft stops . <- Note coef_1 > 0.
                                              ! Should be prevailed by adjust_flag = 2 or 3.
                        w   (:,:,ilev, iupd) = 0.
                        a   (:,:,ilev, iupd) = atd(:,:)
                        entL(:,:,ilev, iupd) = c3(:,:); detL(:,:,ilev, iupd) = 0.
                        entB(:,:,ilev, iupd) = 0.
                        detB(:,:,ilev, iupd) = c1d(:,:)
                    elsewhere
                        w(:,:,ilev,iupd) = - coef_0(:,:) / coef_1(:,:) 
                        
                        entL(:,:,ilev, iupd) = c3(:,:)
                        detL(:,:,ilev, iupd) = d4(:,:) * a(:,:,ilev, iupd) * w(:,:,ilev, iupd)
                        entB(:,:,ilev, iupd) = 0.
                        detB(:,:,ilev, iupd) = c1(:,:)+c2(:,:)+c3(:,:)  - d4(:,:) * a(:,:,ilev, iupd) * w(:,:,ilev, iupd) &
                            - rhot(:,:,ilev)*a(:,:,ilev, iupd)*dz(:,:)*updraft_rescale/dt &
                            - rhot(:,:,ilev)*a(:,:,ilev, iupd) * w(:,:,ilev, iupd)
                    end where ! (coef_0(:,:) >= 0) 
                end where ! (adjust_flag(:,:) == 1) 
             
            endif ! (ilev == nlev)
            m(:,:,ilev, iupd) = rhot(:,:,ilev)*a(:,:,ilev,iupd)*w(:,:,ilev, iupd)
            
        enddo ! (ilev = nlev, 1, -1)
    enddo ! (iupd = 1, nupd)
            
end subroutine compute_au_wu_expentr

! Added ZTAN: 09/28/2017  --  Use explicit entrainment (but still implicit detrainment)
! compute_au_wu: compute updated updraft fraction au and updraft velocity wu, with other auxiliary variables
subroutine compute_au_wu_expentr_newp (a_old, w_old, w_env, rhot, buoy, entr, detr, &
                          pfull, zfull, zhalf, dt, nlev, nupd, &
                          a, w, m, entL, entB, detL, detB, use_lower)

    real   , intent(in),    dimension(:,:,:,:) :: a_old, w_old, buoy, entr, detr
    real   , intent(in),    dimension(:,:,:)   :: w_env, rhot, pfull, zfull, zhalf
    real   , intent(in)                        :: dt
    integer, intent(in)                        :: nlev, nupd
    logical, intent(in)                        :: use_lower
    real   , intent(out),   dimension(:,:,:,:) :: a, w, m, entL, entB, detL, detB
    
    real   , dimension(size(a_old,1),size(a_old,2))  :: dz, c1, c2, c3, d4, db, w1, w2, w3, &
                                                        w_adv, b_adv, coef_0, coef_1, coef_2, a_denom, &
                                                        c1s, c1d, atd, f_1, f_2
    integer, dimension(size(a_old,1),size(a_old,2))  :: adjust_flag
    
    ! c1: old mass; c2: mass from below; c3: (entrained mass); d4: (detrained mass)/wu.
    ! NOTE: d4 is of opposite sign against the compute_au_wu subroutine.
    ! w1: old wu; w2: wu from below; w3: w of entrained mass.
    ! w_adv: wu without buoyancy effect, (coef_0, coef_1, coef_2): coefficients of the quadratic equation.
    ! adjust_flag: flag for updraft fraction limiter.
    
    
    integer :: iupd, ilev
    real :: alpha_b  = 1./3.
    real :: beta_b
    real :: alpha_d  = 0.075
    real :: radius_d = 1500.
    real   , dimension(size(a_old,1),size(a_old,2))  :: beta_d
    
    real :: limit_fac = 100000.
    ! logical :: use_lower = .true.
    logical :: limit_entr = .true.

    beta_b = 1 - alpha_b
    
    a = 0.; w = 0.; m = 0.; entL = 0.; entB = 0.; detL = 0.; detB = 0.
    dz = 0.; c1 = 0.; c2 = 0.; c3 = 0.; d4 = 0.; db = 0.; w1 = 0.; w2 = 0.; w3 = 0.
    w_adv = 0.; b_adv = 0.; coef_0 = 0.; coef_1 = 0.; coef_2 = 0.; a_denom = 0.
    c1s = 0.; c1d = 0.; atd = 0.; f_1 = 0.; f_2 = 0.
    
    do iupd = 1, nupd
        do ilev = nlev, 1, -1
            dz(:,:) = zhalf(:,:,ilev) - zhalf(:,:,ilev+1)
            c1(:,:) = rhot(:,:,ilev)*a_old(:,:,ilev,iupd)*dz(:,:)*updraft_rescale/dt
            w1(:,:) = w_old(:,:,ilev,iupd)
            w3(:,:) = w_env(:,:,ilev)
            
            if (ilev /= nlev) then
                w2(:,:) = w(:,:,ilev+1,iupd)
                c2(:,:) = rhot(:,:,ilev+1)*a(:,:,ilev+1, iupd)*w(:,:,ilev+1, iupd)
                if (use_lower) then
                    c3(:,:) = rhot(:,:,ilev+1)*a(:,:,ilev+1, iupd)*w(:,:,ilev+1, iupd)*min(entr(:,:,ilev,iupd)*dz(:,:), limit_fac)
                else ! use old MF
                    c3(:,:) = rhot(:,:,ilev)*a_old(:,:,ilev,iupd)*w_old(:,:,ilev,iupd)*min(entr(:,:,ilev,iupd)*dz(:,:), limit_fac)
                end if
                ! For diagnostic scheme: c3(:,:) = c2(:,:)*entr(:,:,ilev,iupd)*dz(:,:)
                
                ! Note 10/10/2017: For stability, entrainment should be smaller than environmental airmass
                if (limit_entr) then
                    c3(:,:) = min(c3(:,:), rhot(:,:,ilev)*(1-a_old(:,:,ilev,iupd)*updraft_rescale)*dz(:,:)/dt)
                    ! This implies: a*w*entr*dz < (1-a)dz/dt --> entr < (1-a)/a /(w*dt)
                    ! If: entr = 1/tau/w --> dt < (1-a)/a * tau :: This is satisfied...
                end if

                ! Compute detrained mass coefficients with simple formula
                d4(:,:) = rhot(:,:,ilev)*min(detr(:,:,ilev, iupd)*dz(:,:), limit_fac)
            else
                c2(:,:) = 0.; w2(:,:) = 0.; c3(:,:) = 0.; d4(:,:) = 0.
            endif  ! (ilev /= nlev)
            
            ! write(*,*), ilev, maxval(d3(:,:)), minval(d3(:,:)), maxval(d4(:,:)), minval(d4(:,:))
            db(:,:) = rhot(:,:,ilev)*dz(:,:)
            beta_d(:,:) = alpha_d/radius_d/sqrt(max(a_old(:,:,ilev,iupd), 1.e-3))

            ! ------ Step 1 ------
            ! Assume entB = 0 and detB = 0, and calculate {w},{a},{entL},{detL}.
            ! w-equation: (c1 + c2 + c3) * {w} = ...
            !             (c1 * w1) + (c2 * w2) + (c3 * w3) + ...
            !             (db * beta_b) *{a}* buoy - db * {a} * beta_d * ({w} - w3)^2
            ! a-equation: (c1 + c2 + c3 - d4*{a}*{w}) = ...
            !              rhot*{a} *dz*updraft_rescale/dt + rhot*{a}*{w}
            ! Solve {w} first, then {a}, then compute entL and detL if {a} satisfies condition.
            ! Note: {a} = (c1+c2+c3) / [(d4+rhot)*{w} + rhot*dz*updraft_rescale/dt]
            ! Thus: [{w} - w_adv] * [(d4+rhot)*{w} + rhot*dz*updraft_rescale/dt]
            !           = (db * beta_b) * buoy - db *beta_d * ({w} - w3)^2 
            ! where w_adv =  = [(c1 * w1) + (c2 * w2) + (c3 * w3)]/(c1+c2+c3)
            
            where (c1(:,:) + c2(:,:) + c3(:,:) <= 0.)  !   no prev updraft and no incoming updraft
                ! ==> no current updraft and no massive entrainment. a_old should be zero here.
                w(:,:,ilev,iupd) = 0.; a(:,:,ilev,iupd) = a_old(:,:,ilev,iupd)
                entL(:,:,ilev,iupd) = 0.; entB(:,:,ilev,iupd) = 0.
                detL(:,:,ilev,iupd) = 0.; detB(:,:,ilev,iupd) = 0.
            elsewhere
                w_adv (:,:) = (c1(:,:)*w1(:,:)+c2(:,:)*w2(:,:)+c3(:,:)*w3(:,:))/(c1(:,:)+c2(:,:)+c3(:,:))
                coef_2(:,:) = d4(:,:) + rhot(:,:,ilev) + db(:,:) * beta_d(:,:) ! > 0 because d4 >= 0
                coef_1(:,:) = rhot(:,:,ilev)*dz(:,:)*updraft_rescale/dt &
                              - (rhot(:,:,ilev)+d4(:,:))* w_adv(:,:) & 
                              - 2. * db(:,:) * beta_d(:,:) * w3(:,:) 
                coef_0(:,:) = - w_adv(:,:) * rhot(:,:,ilev)*dz(:,:)*updraft_rescale/dt & 
                              - db(:,:) * beta_b * buoy(:,:,ilev, iupd) &
                              + db(:,:) * beta_d(:,:) * w3(:,:) * w3(:,:)
                              
                where (coef_0(:,:) >= 0) 
                ! Both w_adv and buoy <= 0. Thus coef_1 >= 0, 
                ! Note that coef_2 > 0. Thus, both roots are non-positive.
                ! Updraft stops, boundary detrainment occurs (Step 3).
                    w   (:,:,ilev, iupd) = 0.; a(:,:,ilev, iupd) = a_old(:,:,ilev,iupd) ! a will be adjusted by Step 3
                    entL(:,:,ilev, iupd) = 0.; detL(:,:,ilev, iupd) = 0.
                    entB(:,:,ilev, iupd) = 0.; detB(:,:,ilev, iupd) = 0.
                elsewhere
                    w (:,:,ilev, iupd) = (-coef_1(:,:) + &
                        sqrt(coef_1(:,:)**2. - 4.*coef_2(:,:)*coef_0(:,:)))/(2.*coef_2(:,:))
                    
                    a_denom (:,:) =  rhot(:,:,ilev)*dz(:,:)*updraft_rescale/dt + &
                                    (rhot(:,:,ilev)+d4(:,:))*w(:,:,ilev, iupd)
                    ! This denominator is always positive.
                    
                    a(:,:,ilev, iupd) = (c1(:,:) + c2(:,:) + c3(:,:))/a_denom (:,:)
                    
                    entB(:,:,ilev, iupd) = 0.; detB(:,:,ilev, iupd) = 0.
                    entL(:,:,ilev, iupd) = c3(:,:) 
                    detL(:,:,ilev, iupd) = d4(:,:) * a(:,:,ilev, iupd) * w(:,:,ilev, iupd)
                end where ! (coef_0(:,:) >= 0) 
            end where ! (c1(:,:) + c2(:,:) <= 0.)


            ! write(*,*), ilev, maxval(entL(:,:)), minval(entL(:,:)), maxval(detL(:,:)), minval(detL(:,:))
            
            if (ilev == nlev) then
            ! ------ Step 2 ------
            ! Fix SFC fraction by:
            ! (1) If buoy > 0 => w > 0, aN < a0 because entL = 0 => entB > 0.
            ! (2) If buoy <= 0. Swap airmass, use detB = c1, calculate entB accordingly.

                ! Calculate entB for target a(ilev, iupd) = a_o(ilev, iupd)
                ! final w-equation: (c1 + c2 + c3 + {entB}) * {w} = ...
                !                   (c1 * w1) + (c2 * w2) + (c3 + {entB}) * w3 + ...
                !                   (db * beta_b * a) * buoy - db * a * beta_d * ({w} - w3)^2
                ! final a-equation: (c1 + c2 + c3 - d4*a*{w} + {entB} - {detB}) = ...
                !                    rhot*a *dz*updraft_rescale/dt + rhot*a*{w}
                !                  -> f_2 * {w} + f_1 = {entB}
                ! Thus (c1 + c2 + c3 + f_2 * {w} + f_1) * {w} = ...
                !                   (c1 * w1) + (c2 * w2) + (c3 + f_2 * {w} + f_1) * w3 + ...
                !                   (db * beta_b * a) * buoy - db * a * beta_d * ({w} - w3)^2
                
                detB(:,:,ilev, iupd) = 0.0 ! c1s(:,:)  ! Changed ZTAN 09/26/2017: Do not swap 
                    
                a(:,:,ilev, iupd) = a_old(:,:,ilev, iupd)
                f_1(:,:) = rhot(:,:,ilev)*a(:,:,ilev, iupd)*dz(:,:)*updraft_rescale/dt - &
                      (c1(:,:)+c2(:,:)+c3(:,:) - detB(:,:,ilev, iupd))
                f_2(:,:) = (rhot(:,:,ilev) + d4(:,:))*a(:,:,ilev, iupd)

                coef_2(:,:) = f_2(:,:) + db(:,:)*a(:,:,ilev, iupd)* beta_d(:,:)    ! coef_2 > 0 by definition
                coef_1(:,:) = c1(:,:)+c2(:,:)+ c3(:,:) + f_1(:,:) - w3(:,:)*f_2(:,:) &
                              - 2. * db(:,:) * beta_d(:,:) * a(:,:,ilev, iupd) * w3(:,:)
                coef_0(:,:) = -(c1(:,:)*w1(:,:) + c2(:,:)*w2(:,:) + c3(:,:)*w3(:,:) + f_1(:,:)*w3(:,:) &
                                + db(:,:) * beta_b * a(:,:,ilev, iupd)*buoy(:,:,ilev, iupd) &
                                - db(:,:) * beta_d(:,:) * a(:,:,ilev, iupd) * w3(:,:) * w3(:,:) )
                               
                where (coef_0(:,:) >= 0.)  ! Updraft at SFC is not buoyant
                    w(:,:,ilev, iupd) = 0.0
                    entL(:,:,ilev, iupd) = 0.0; detL(:,:,ilev, iupd) = 0.0
                    entB(:,:,ilev, iupd) = rhot(:,:,ilev)*a(:,:,ilev, iupd)*dz(:,:)*updraft_rescale/dt - c2(:,:) - c3(:,:)
                    detB(:,:,ilev, iupd) = c1(:,:)
                elsewhere
                    w (:,:,ilev, iupd) = (-coef_1(:,:) + &
                        sqrt(coef_1(:,:)**2. - 4.*coef_2(:,:)*coef_0(:,:)))/(2.*coef_2(:,:))
                    entB(:,:,ilev, iupd) = f_1(:,:) + f_2(:,:) * w(:,:,ilev, iupd)
                    entL(:,:,ilev, iupd) = c3(:,:) 
                    detL(:,:,ilev, iupd) = d4(:,:) * a(:,:,ilev, iupd) * w(:,:,ilev, iupd)
                end where ! (coef_0 >= 0.)
                
                ! Adjust entB and detB for scalars -- switch updraft thl/qt with environment buoyant tail
                ! ZTAN: 09/25/2017: This is disabled due to inconsistency. The 'switch' is now done at the beginning.
                ! detB(:,:,ilev, iupd) = detB(:,:,ilev, iupd) + c1(:,:)*99. ! 99
                ! entB(:,:,ilev, iupd) = entB(:,:,ilev, iupd) + c1(:,:)*99. ! 99
            else
            ! ------ Step 3 ------
            ! Adjust area fraction at other levels by tuning up detL.
            ! Note: this algorithm is applicable to adjusting down only.
                
                ! Limitation on w_u
                adjust_flag(:,:) = 0
                if (au_optB_wu == 1) then ! detraining parcels with zero velocity
                    where ( (w(:,:,ilev, iupd) <= 1.0e-6) .or. (a(:,:,ilev, iupd) <= 1.0e-6) )
                        adjust_flag(:,:) = 2  ! Kill the updraft
                    end where
                elseif (au_optB_wu == 11) then ! detraining parcels with zero velocity or zero area
                    where ( (a(:,:,ilev, iupd) <= 1.0e-6))
                        adjust_flag(:,:) = 2  ! Kill the updraft
                    elsewhere
                        where (w(:,:,ilev, iupd) <= 1.0e-6)
                            adjust_flag(:,:) = 3  ! Kill the updraft gradually
                        end where
                    end where
                endif ! (au_optB_wu == 1) (au_optB_wu == 11)
            
            
                ! Limitation on a_u is top priority -- will overwrite adjust_flag due to w_u.
                if (au_optB == 1) then    ! constrain updraft fraction < SFC fraction
                    where (a(:,:,ilev, iupd)> a(:,:,nlev, iupd))
                        a(:,:,ilev,iupd) = a(:,:,nlev, iupd);  adjust_flag(:,:) = 1
                    end where
                elseif (au_optB == 2) then    ! constrain updraft fraction < 2*SFC fraction
                    where (a(:,:,ilev, iupd)> 2.0*a(:,:,nlev, iupd))
                        a(:,:,ilev,iupd) = 2.0*a(:,:,nlev, iupd);  adjust_flag(:,:) = 1
                    end where
                elseif (au_optB == 100) then  ! constrain updraft fraction < 0.5/nupd
                    where (a(:,:,ilev, iupd)> .5/nupd)
                        a(:,:,ilev,iupd) = .5/nupd;  adjust_flag(:,:) = 1
                    end where
                endif ! (au_optB == 1, 2, 100)
                
                where (adjust_flag(:,:) == 3) ! Adjust w to zero directly, 
                                              ! but detrain a over a 15 minute timescale.
                    w   (:,:,ilev, iupd) = 0. 
                    a   (:,:,ilev, iupd) = a_old(:,:,ilev,iupd) * exp(-dt/900.0)   
                    entL(:,:,ilev, iupd) = c3(:,:); detL(:,:,ilev, iupd) = 0.
                    entB(:,:,ilev, iupd) = 0. 
                    detB(:,:,ilev, iupd) = c1(:,:)+c2(:,:)+c3(:,:) - &
                         rhot(:,:,ilev)*a(:,:,ilev, iupd)*dz(:,:)*updraft_rescale/dt
                end where ! (adjust_flag(:,:) == 3) 
                
                where (adjust_flag(:,:) == 2) ! Adjust to zero directly.
                    w   (:,:,ilev, iupd) = 0.  
                    a   (:,:,ilev, iupd) = 0.  
                    entL(:,:,ilev, iupd) = c3(:,:); detL(:,:,ilev, iupd) = 0.
                    entB(:,:,ilev, iupd) = 0.
                    detB(:,:,ilev, iupd) = c1(:,:)+c2(:,:)+c3(:,:)
                end where ! (adjust_flag(:,:) == 2) 
                
                ! c1d is only used where adjust_flag == 1
                if (au_optB_wu == 1) then   ! Detrainment everything. atd is target a.
                    c1d(:,:) = c1(:,:)+c2(:,:)+c3(:,:)
                    atd(:,:) = 0.
                elseif (au_optB_wu == 11) then ! Detrainment to target a (= atd).
                    c1d(:,:) = c1(:,:)+c2(:,:)+c3(:,:) &
                               - rhot(:,:,ilev)*a(:,:,ilev, iupd)*dz(:,:)*updraft_rescale/dt
                    atd(:,:) = a(:,:,ilev, iupd)
                endif
                
                where (adjust_flag(:,:) == 1)
                    ! Calculate detB for target a(ilev, iupd) = a_o(ilev, iupd), with entB = 0
                    ! final w-equation: (c1 + c2 + c3) * {w} = ...
                    !                   (c1 * w1) + (c2 * w2) + c3 * w3  + ...
                    !                   (db * beta_b) *a* buoy - db * a * beta_d * ({w} - w3)^2
                    ! final a-equation: (c1 + c2 + c3 - d4*a{w} - {detB}) = ...
                    !                    rhot*a *dz*updraft_rescale/dt + rhot*a*{w}
                    
                    coef_2(:,:) = db(:,:) * beta_d(:,:) *a(:,:,ilev, iupd)  ! coef_2 > 0 by definition
                    coef_1(:,:) = c1(:,:)+c2(:,:) + c3(:,:) &
                                  - 2.* db(:,:) * beta_d(:,:) *a(:,:,ilev, iupd) * w3(:,:)
                    coef_0(:,:) = -(c1(:,:)*w1(:,:) + c2(:,:)*w2(:,:) + c3(:,:)*w3(:,:) + &
                                    + db(:,:) * beta_b * a(:,:,ilev, iupd)*buoy(:,:,ilev, iupd) &
                                    - db(:,:) * beta_d(:,:)* a(:,:,ilev, iupd) * w3(:,:) * w3(:,:))
                                                 
                    where (coef_0(:,:) >= 0.) ! Updraft stops . <- Note coef_1 > 0.
                                              ! Should be prevailed by adjust_flag = 2 or 3.
                        w   (:,:,ilev, iupd) = 0.
                        a   (:,:,ilev, iupd) = atd(:,:)
                        entL(:,:,ilev, iupd) = c3(:,:); detL(:,:,ilev, iupd) = 0.
                        entB(:,:,ilev, iupd) = 0.
                        detB(:,:,ilev, iupd) = c1d(:,:)
                    elsewhere
                        where (coef_2(:,:) == 0.) 
                            ! d3 = 0, no entrainment, degenerate equation.
                            ! This implies that c1+c2 > 0, thus coef_1 > 0, coef_0 < 0 ==> w > 0.
                            w(:,:,ilev,iupd) = - coef_0(:,:) / coef_1(:,:) 
                        elsewhere  
                            ! Now coef_2 > 0 and coef_0 < 0 ==> w > 0
                            w(:,:,ilev,iupd) = (-coef_1(:,:) + &
                                sqrt(coef_1(:,:)**2. - 4.*coef_2(:,:)*coef_0(:,:)))/(2.*coef_2(:,:))
                        end where 
                        
                        entL(:,:,ilev, iupd) = c3(:,:)
                        detL(:,:,ilev, iupd) = d4(:,:) * a(:,:,ilev, iupd) * w(:,:,ilev, iupd)
                        entB(:,:,ilev, iupd) = 0.
                        detB(:,:,ilev, iupd) = c1(:,:)+c2(:,:)+c3(:,:)  - d4(:,:) * a(:,:,ilev, iupd) * w(:,:,ilev, iupd) &
                            - rhot(:,:,ilev)*a(:,:,ilev, iupd)*dz(:,:)*updraft_rescale/dt &
                            - rhot(:,:,ilev)*a(:,:,ilev, iupd) * w(:,:,ilev, iupd)
                    end where ! (coef_0(:,:) >= 0) 
                end where ! (adjust_flag(:,:) == 1) 
             
            endif ! (ilev == nlev)
            m(:,:,ilev, iupd) = rhot(:,:,ilev)*a(:,:,ilev,iupd)*w(:,:,ilev, iupd)
            
        enddo ! (ilev = nlev, 1, -1)
    enddo ! (iupd = 1, nupd)
            
end subroutine compute_au_wu_expentr_newp


! compute_wu: compute updated updraft fraction au, and entrainment/detrainment masses (entL, entB, detL, detB)
subroutine compute_au (w, a_old, rhot, entr, detr, pfull, zfull, zhalf, dt, nlev, nupd, &
                       a, m, entL, entB, detL, detB)

    real   , intent(in),    dimension(:,:,:,:) :: w, a_old, entr, detr
    real   , intent(in),    dimension(:,:,:)   :: rhot, pfull, zfull, zhalf
    real   , intent(in)                        :: dt
    integer, intent(in)                        :: nlev, nupd
    real   , intent(out),   dimension(:,:,:,:) :: a, m, entL, entB, detL, detB

    real   , dimension(size(a_old,1),size(a_old,2),size(a_old,3),size(a_old,4)) :: ra_old, ra_new
    real   , dimension(size(a_old,1),size(a_old,2))  ::   dz, fact, a_new
    integer, dimension(size(a_old,1),size(a_old,2))  ::   adj_au_flag
    
    integer :: iupd, ilev
        
! options: au_optL, au_optB, au_optB_wu, au_optB_srf
    
    do iupd = 1, nupd
        ra_old(:,:,:,iupd) = a_old(:,:,:,iupd) * rhot(:,:,:)
    end do
    ra_new = 0.; a = 0.; m = 0.; entL = 0.; entB = 0.; detL = 0.; detB = 0.
    
    do iupd = 1, nupd
        do ilev = nlev, 1, -1
            dz(:,:) = zhalf(:,:,ilev) - zhalf(:,:,ilev+1)
            
            ! (1) Compute entL and detL
            if (au_optL == 1) then !  equation (119): z=z+
                if (ilev == nlev) then
                    ra_new(:,:,ilev,iupd) = ra_old(:,:,ilev,iupd)
                    entL(:,:,ilev,iupd) = 0.0
                    detL(:,:,ilev,iupd) = 0.0
                    entB(:,:,ilev,iupd) = ra_new(:,:,ilev,iupd) * w(:,:,ilev,iupd)
                else
                    ra_new(:,:,ilev,iupd) = ra_old(:,:,ilev,iupd) * dz(:,:)*updraft_rescale/dt + & !  Numerator
                        ra_new(:,:,ilev+1,iupd) * w(:,:,ilev+1,iupd)*exp(dz(:,:)*(entr(:,:,ilev,iupd)-detr(:,:,ilev, iupd)))
                    ra_new(:,:,ilev,iupd) = ra_new(:,:,ilev,iupd) / (dz(:,:)*updraft_rescale/dt + w(:,:,ilev,iupd)) ! Denominator
                    
                    where (entr(:,:,ilev,iupd)-detr(:,:,ilev, iupd) == 0.0)
                        fact(:,:) = dz(:,:)
                    elsewhere
                        fact(:,:) = exp(dz(:,:)*(entr(:,:,ilev, iupd)-detr(:,:,ilev, iupd)))-1.
                        fact(:,:) = fact(:,:)/(entr(:,:,ilev, iupd)-detr(:,:,ilev, iupd))
                    end where
                    entL(:,:,ilev, iupd) = entr(:,:,ilev, iupd)*ra_new(:,:,ilev+1, iupd)*w(:,:,ilev+1,iupd)*fact(:,:)
                    detL(:,:,ilev, iupd) = detr(:,:,ilev, iupd)*ra_new(:,:,ilev+1, iupd)*w(:,:,ilev+1,iupd)*fact(:,:)
                    entB(:,:,ilev, iupd) = 0.
                end if
                
            elseif (au_optL == 2) then  !  equation (120): z=z0
                if (ilev == nlev) then
                    ra_new(:,:,ilev,iupd) = ra_old(:,:,ilev,iupd)
                    where (entr(:,:,ilev,iupd)-detr(:,:,ilev, iupd) == 0.0)
                        fact(:,:) = dz(:,:)
                    elsewhere
                        fact(:,:) = 1.-exp(-dz(:,:)*(entr(:,:,ilev, iupd)-detr(:,:,ilev, iupd)))
                        fact(:,:) = fact(:,:)/(entr(:,:,ilev, iupd)-detr(:,:,ilev, iupd))       
                    end where
                    entL(:,:,ilev, iupd) = entr(:,:,ilev, iupd)*ra_new(:,:,ilev, iupd)*w(:,:,ilev,iupd)*fact(:,:)
                    detL(:,:,ilev, iupd) = detr(:,:,ilev, iupd)*ra_new(:,:,ilev, iupd)*w(:,:,ilev,iupd)*fact(:,:)
                    entB(:,:,ilev, iupd) = ra_new(:,:,ilev, iupd)*w(:,:,ilev, iupd)* & 
                                           exp(-dz(:,:)*(entr(:,:,ilev, iupd)-detr(:,:,ilev, iupd)))
                    ! Assumption: at the surface, lateral detrainment rate < lateral entrainment rate
                else
                    ra_new(:,:,ilev, iupd) = ra_old(:,:,ilev, iupd)*dz(:,:)*updraft_rescale/dt + &  ! Numerator
                         ra_new(:,:,ilev+1, iupd)*w(:,:,ilev+1,iupd)
                    ra_new(:,:,ilev, iupd) = ra_new(:,:,ilev, iupd)/ &        ! Denominator
                         (dz(:,:)*updraft_rescale/dt + w(:,:,ilev,iupd)* & 
                          exp(-dz(:,:)*(entr(:,:,ilev, iupd)-detr(:,:,ilev, iupd))))
                    where (entr(:,:,ilev,iupd)-detr(:,:,ilev, iupd) == 0.0)
                        fact(:,:) = dz(:,:)
                    elsewhere
                        fact(:,:) = 1.-exp(-dz(:,:)*(entr(:,:,ilev, iupd)-detr(:,:,ilev, iupd)))
                        fact(:,:) = fact(:,:)/(entr(:,:,ilev, iupd)-detr(:,:,ilev, iupd))
                    end where
                    entL(:,:,ilev, iupd) = entr(:,:,ilev, iupd)*ra_new(:,:,ilev, iupd)*w(:,:,ilev,iupd)*fact(:,:)
                    detL(:,:,ilev, iupd) = detr(:,:,ilev, iupd)*ra_new(:,:,ilev, iupd)*w(:,:,ilev,iupd)*fact(:,:)        
                    entB(:,:,ilev, iupd) = 0.
                end if
            
            end if
            
            a(:,:,ilev,iupd) = ra_new(:,:,ilev, iupd)/rhot(:,:,ilev)
            adj_au_flag (:,:) = 0
            a_new(:,:) = a(:,:,ilev,iupd)
            
            ! (2) Compute entB, detB (TO MEET a_u REQUIREMENTS), and update a_u
            
            if (au_optB == 1) then  ! DEFAULT: fixed surface area fraction, and constrain
                                    ! updraft fraction to be less than some multiples of surface fraction
                if (au_optL == 1) then ! Do nothing
                    where (a(:,:,ilev, iupd) .gt. a(:,:,nlev, iupd) * au_optB1_frac) ! Adjust: note that ilev < nlev is satisfied
                        a_new(:,:) = a(:,:, nlev, iupd) * au_optB1_frac              ! Note that a_new(:,:) < a(:,:,ilev, iupd)
                        adj_au_flag(:,:) = 1
                    end where
                elseif (au_optL == 2) then ! need to update (rescale) entL, detL
                    where (a(:,:,ilev, iupd) .gt. a(:,:,nlev, iupd) * au_optB1_frac) ! Adjust: note that ilev < nlev is satisfied
                        a_new(:,:) = a(:,:, nlev, iupd) * au_optB1_frac              ! Note that a_new(:,:) < a(:,:,ilev, iupd)
                        entL(:,:,ilev, iupd) = a_new(:,:)/a(:,:,ilev, iupd)*entL(:,:,ilev, iupd)
                        detL(:,:,ilev, iupd) = a_new(:,:)/a(:,:,ilev, iupd)*detL(:,:,ilev, iupd)
                        adj_au_flag(:,:) = 1
                    end where
                else
                    where (a(:,:,ilev, iupd) .gt. a(:,:,nlev, iupd) * au_optB1_frac) ! Adjust: note that ilev < nlev is satisfied
                        a_new(:,:) = a(:,:, nlev, iupd) * au_optB1_frac              ! Note that a_new(:,:) < a(:,:,ilev, iupd)
                        adj_au_flag(:,:) = 1
                    end where
                end if
                            
            elseif (au_optB == 2) then ! fixed surface area fraction, and constrain
                                       ! updraft fraction to be non-increasing                       
                if (ilev .lt. nlev) then
                    if (au_optL == 1) then ! Do nothing
                        where (a(:,:,ilev, iupd) .gt. a(:,:,ilev+1, iupd))   ! Adjust: note that ilev < nlev is satisfied
                            a_new(:,:) = a(:,:,ilev+1, iupd)     ! Note that a_new < a(ilev, iupd)
                            adj_au_flag(:,:) = 1
                        end where
                    elseif(au_optL == 2) then ! need to update (rescale) entL, detL
                        where (a(:,:,ilev, iupd) .gt. a(:,:,ilev+1, iupd))   ! Adjust: note that ilev < nlev is satisfied
                            a_new(:,:) = a(:,:,ilev+1, iupd)     ! Note that a_new < a(ilev, iupd)
                            entL(:,:,ilev, iupd) = a_new(:,:)/a(:,:,ilev, iupd)*entL(:,:,ilev, iupd)
                            detL(:,:,ilev, iupd) = a_new(:,:)/a(:,:,ilev, iupd)*detL(:,:,ilev, iupd)
                            adj_au_flag(:,:) = 1
                        end where
                    else
                        where (a(:,:,ilev, iupd) .gt. a(:,:,ilev+1, iupd))   ! Adjust: note that ilev < nlev is satisfied
                            a_new(:,:) = a(:,:,ilev+1, iupd)     ! Note that a_new < a(ilev, iupd)
                            adj_au_flag(:,:) = 1
                        end where
                    end if
                end if
            
            elseif (au_optB == 10) then ! fixed area fraction everywhere
                a_new(:,:) = a(:,:, nlev, iupd)
                if    (au_optL == 1) then ! Do nothing
                elseif(au_optL == 2) then ! need to update (rescale) entL, detL
                    entL(:,:,ilev, iupd) = a_new(:,:)/a(:,:,ilev, iupd)*entL(:,:,ilev, iupd)
                    detL(:,:,ilev, iupd) = a_new(:,:)/a(:,:,ilev, iupd)*detL(:,:,ilev, iupd)
                end if
                adj_au_flag(:,:) = 1
            end if
            
            if (au_optB_wu == 1)  then ! detraining parcels with zero velocity
                if (ilev == nlev) then ! adjust entrainment
                    where (w(:,:,ilev,iupd) .le. wu_min)
                        a_new(:,:) = a(:,:,nlev, iupd)
                        entB(:,:,ilev, iupd)  = entB(:,:,ilev, iupd) + &
                            ra_old(:,:,ilev, iupd)*dz(:,:)*updraft_rescale/dt
                        adj_au_flag(:,:) = 1
                    end where
                else  ! adjust area fraction: detrain everything
                    where (w(:,:,ilev,iupd) .le. wu_min)
                        a_new(:,:) = 0.0
                        adj_au_flag(:,:) = 1
                    end where
                end if
                 
            elseif (au_optB_wu == 2) then  ! detraining parcels with zero velocity naturally
                where (w(:,:,ilev,iupd) .le. wu_min)
                    a_new(:,:) = a(:,:,nlev, iupd)
                    entB(:,:,ilev, iupd)  = entB(:,:,ilev, iupd) + & 
                        ra_old(:,:,ilev, iupd)*dz(:,:)*updraft_rescale/dt
                    adj_au_flag(:,:) = 1
                end where
            end if
            
            if (au_optB_srf == 1) then
                 if (ilev == nlev) then
                     where (adj_au_flag(:,:) == 0)  ! adjust entrainment
                         a_new(:,:) = a(:,:,nlev, iupd)
                         entB(:,:,ilev, iupd)  = entB(:,:,ilev, iupd) + &
                             ra_old(:,:,ilev, iupd)*dz(:,:)*updraft_rescale/dt
                         adj_au_flag(:,:) = 1
                     end where
                 end if
            end if

            ! Diagnose detB that meets the requirement of a_new
            
            if (ilev == nlev) then
                where (adj_au_flag(:,:) == 1)
                    detB(:,:,ilev, iupd) = -rhot(:,:,ilev)*a_new(:,:)*(dz(:,:)*updraft_rescale/dt+w(:,:,ilev, iupd)) + &
                        ra_old(:,:,ilev, iupd)*dz(:,:)*updraft_rescale/dt +   &
                        (entL(:,:,ilev, iupd) - detL(:,:,ilev, iupd) + entB(:,:,ilev, iupd))
                    where (detB(:,:,ilev, iupd) .lt. 0.0)
                        entB(:,:,ilev, iupd) = entB(:,:,ilev, iupd) - detB(:,:,ilev, iupd)
                        detB(:,:,ilev, iupd) = 0.0
                    end where
                    a(:,:,ilev, iupd) = a_new(:,:)
                    ra_new(:,:,ilev, iupd) = rhot(:,:,ilev)*a_new(:,:)
                end where
            else
                where (adj_au_flag(:,:) == 1)
                    detB(:,:,ilev, iupd) = -rhot(:,:,ilev)*a_new(:,:)*(dz(:,:)*updraft_rescale/dt+w(:,:,ilev, iupd)) + &
                        ra_old(:,:,ilev, iupd)*dz(:,:)*updraft_rescale/dt + ra_new(:,:,ilev+1, iupd)*w(:,:,ilev+1, iupd) + &
                        (entL(:,:,ilev, iupd) - detL(:,:,ilev, iupd) + entB(:,:,ilev, iupd))
                    where (detB(:,:,ilev, iupd) .lt. 0.0)
                        entB(:,:,ilev, iupd) = entB(:,:,ilev, iupd) - detB(:,:,ilev, iupd)
                        detB(:,:,ilev, iupd) = 0.0
                    end where
                    a(:,:,ilev, iupd) = a_new(:,:)
                    ra_new(:,:,ilev, iupd) = rhot(:,:,ilev)*a_new(:,:)
                end where
            end if
                
            
            
            m(:,:,ilev, iupd) = rhot(:,:,ilev)*a(:,:,ilev, iupd)*w(:,:,ilev, iupd)
                       
        end do
    end do
    
end subroutine compute_au


! compute_wu: compute updated updraft velocity wu
subroutine compute_wu ( w_old, w_env, buoy, entr, zfull, zhalf, dt, nlev, nupd, w)

    real   , intent(in),    dimension(:,:,:,:) :: w_old, buoy, entr
    real   , intent(in),    dimension(:,:,:)   :: w_env, zfull, zhalf
    real   , intent(in)                        :: dt
    integer, intent(in)                        :: nlev, nupd
    real   , intent(out),   dimension(:,:,:,:) :: w

    real   , dimension(size(zfull,1),size(zfull,2))  :: dz, wa, wb2, wbl
    integer :: iupd, ilev
    real :: b1  = 2.5 ! = b in thesis. Bretherton 04 -- b = 2.
                      ! The value 2.5 is from a version of EDMF by Chung.
                      ! Soares 2004 and Witek 2011 use b1 = 1 (which means ZERO damping!!).
    real :: b2  = 2.0 ! = a in thesis. Bretherton 04 -- a = 1; 
                      ! Soares 2004 and Witek 2011 use b2 = 2 (which means NEGATIVE virtual mass!!).
    
    w = 0.0
    
    if (wu_opt == 0) then   !  simple explicit scheme, not coded yet
    
    elseif (wu_opt == 1) then  !  simple implicit scheme
        do iupd = 1, nupd
            do ilev = nlev, 1, -1
                dz(:,:) = zhalf(:,:,ilev) - zhalf(:,:,ilev+1)
                ! Limit wa to be positive (which may not hold if w_env >> 0, or more specifically, w_env >~ dz/dt)
                wa(:,:) = max(dz(:,:)*updraft_rescale*(1./dt - b1*w_env(:,:,ilev)*entr(:,:,ilev, iupd)), 0.0)
                if (ilev == nlev) then
                    wb2(:,:) = (1.+2.*b1*entr(:,:,ilev, iupd)*dz(:,:)) * &
                        ( 2.*dz(:,:)*updraft_rescale/dt * w_old(:,:,ilev,iupd) + & 
                          2.*b2*buoy(:,:,ilev,iupd)*dz(:,:))
                else
                    wb2(:,:) = (1.+2.*b1*entr(:,:,ilev, iupd)*dz(:,:)) * &
                        ((w(:,:,ilev+1,iupd))**2.+2.*dz(:,:)*updraft_rescale/dt * w_old(:,:,ilev,iupd) + & 
                          2.*b2*buoy(:,:,ilev,iupd)*dz(:,:))
                end if
                wb2(:,:) = max(wb2(:,:), 0.0)                
                w(:,:,ilev,iupd) = max( (-wa(:,:)+sqrt(wa(:,:)*wa(:,:) + wb2(:,:))) / & 
                                        (1.+2.*b1*entr(:,:,ilev,iupd)*dz(:,:)), wu_min) 
            end do
        end do
    elseif (wu_opt == 2) then  !  integrated implicit scheme 
        do iupd = 1, nupd
            do ilev = nlev, 1, -1
                dz(:,:) = zhalf(:,:,ilev) - zhalf(:,:,ilev+1)
                wa(:,:) = max(dz(:,:)*updraft_rescale/dt, 0.0)
                where (entr(:,:,ilev,iupd) .gt. 0.)
                    wbl(:,:) = b2 * buoy(:,:,ilev,iupd)/b1/entr(:,:,ilev,iupd)*(1.-exp(-2.*b1*entr(:,:,ilev,iupd)*dz(:,:)))
                elsewhere
                    wbl(:,:) = 2. * b2 * buoy(:,:,ilev,iupd)*dz(:,:)
                end where
                
                if (ilev == nlev) then
                    wb2(:,:) =  (w(:,:,ilev+1,iupd))**2.*exp(-2.*b1*entr(:,:,ilev,iupd)*dz(:,:))+ &
                                  2.*dz(:,:)*updraft_rescale/dt* w_old(:,:,ilev,iupd) + wbl(:,:)
                else
                    wb2(:,:) = 2.*dz(:,:)*updraft_rescale/dt* w_old(:,:,ilev,iupd) + wbl(:,:)
                end if
                wb2(:,:) = max(wb2(:,:), 0.0)
                w(:,:,ilev,iupd) = max( -wa(:,:)+sqrt(wa(:,:)*wa(:,:) + wb2(:,:)), wu_min)
                                        
            end do
        end do
    end if
    
end subroutine compute_wu


! compute_ent_det: compute entrainment and detrainment rates
subroutine compute_ent_det (w_upd, buoy, chi_c, zstar, zfull, nlev, nupd, entr, detr) 

    real   , intent(in) ,   dimension(:,:,:,:) :: w_upd, buoy, chi_c
    real   , intent(in) ,   dimension(:,:)     :: zstar
    real   , intent(in) ,   dimension(:,:,:)   :: zfull
    integer, intent(in)                        :: nlev, nupd
    real   , intent(out),   dimension(:,:,:,:) :: entr, detr
    
    real , dimension(size(w_upd,1),size(w_upd,2),size(w_upd,4)) :: zt 
    real :: c1  = 0.4
    real :: a1  = 0.15
    real :: a3  = 2.0
    real :: tau = 600.0
    
    ! used in b/wu^2 entrainment/detrainment only
    real :: c_eps = 0.12 ! 0.33 previously
    real :: c_del = 0.12 
    real :: del_b = 0.004
    
    
    integer :: iupd, ilev
    
    entr = 0.
    detr = 0.
    
    if (ER == 0) then  ! Simple scaling: k/z
        do iupd = 1, nupd
            entr(:,:,:,iupd) = c1 / max(zfull(:,:,:), ER0_zmin)
        end do
    elseif (ER == 1) then
        do iupd = 1, nupd
            do ilev = 1, nlev
                where (zfull(:,:,ilev) .ge. zstar(:,:))
                    entr(:,:,ilev,iupd) = 0.0
                elsewhere
                    entr(:,:,ilev,iupd) = c1 * (1./zfull(:,:,ilev) + 1./(zstar(:,:) - zfull(:,:,ilev)) )
                end where
            end do
        end do
        ! ER = 1.33 and ER = 2 in SCM is not coded
    elseif (ER == 3) then
        if (tau > 0.0) then
            do iupd = 1, nupd
                do ilev = 1, nlev-1
                    where (w_upd(:,:,ilev+1, iupd) .gt. 0.0)
                        entr(:,:,ilev,iupd) = a3/(tau*w_upd(:,:,ilev+1, iupd))
                    end where
                end do
            end do
        end if
    elseif (ER == 4) then
        do iupd = 1, nupd
            ! do ilev = 1, nlev-1
            !     where (w_upd(:,:,ilev+1, iupd) .gt. 0.0)
            !         entr(:,:,ilev,iupd) = max(c_eps*buoy(:,:,ilev+1, iupd)/(w_upd(:,:,ilev+1, iupd)**2.0), 0.0)
            !     end where
            ! end do
            ! 
            ! ZTAN 09/26/2017: staggering is eliminated; lateral entrainment formulation is disabled for the lowest level
            do ilev = 1, nlev-1
                entr(:,:,ilev,iupd) = c_eps*max(buoy(:,:,ilev, iupd),0.0)/(max(w_upd(:,:,ilev, iupd),3.0e-2)**2.0)
            end do
        end do 
    elseif (ER == 5) then
        do iupd = 1, nupd
            do ilev = 1, nlev-1
                where (buoy(:,:,ilev, iupd) >= 0.0)
                    entr(:,:,ilev,iupd) = 1.0/c_tau_ER5/max(w_upd(:,:,ilev, iupd),3.0e-2)
                end where
            end do
         end do
    end if
    entr = ER_frac * entr
    
    if (DR == 0) then  ! proportional to entrainment rate
        detr(:,:,:,:) = entr(:,:,:,:)
    elseif (DR == 1)  then  ! TO_DO: 082016 -- THIS MAY NOT BE CORRECT !!!
        call compute_zt_upd(buoy, zfull, nlev, nupd, zt)   
        ! zt: the highest level where buoy > 0.0 (different from SCM !!)
        do iupd = 1, nupd
            do ilev = 1, nlev
                where (zfull(:,:,ilev) .ge. zt(:,:,iupd))
                    entr(:,:,ilev,iupd) = 0.0
                elsewhere
                    entr(:,:,ilev,iupd) = c1 * (1./(zt(:,:,iupd) - zfull(:,:,ilev))) 
                end where
            end do
        end do
    elseif (DR == 4) then
        do iupd = 1, nupd
            ! ZTAN 09/26/2017: staggering is eliminated; lateral detrainment formulation is disabled for the lowest level
            do ilev = 1, nlev-1
                where (zfull(:,:,ilev) .le. zstar(:,:))
                    detr(:,:,ilev,iupd) = 0.0
                elsewhere
                    detr(:,:,ilev,iupd) = del_b + &
                        c_del*max(-buoy(:,:,ilev, iupd),0.0)/(max(w_upd(:,:,ilev, iupd),3.0e-2)**2.0)
                end where
            end do
        end do
    elseif (DR == 5) then
        do iupd = 1, nupd
            do ilev = 1, nlev-1
                where (buoy(:,:,ilev, iupd) <= 0.0)
                    detr(:,:,ilev,iupd) = 1.0/c_tau_ER5/max(w_upd(:,:,ilev, iupd),3.0e-2)
                end where
            end do
         end do
    elseif (DR == 10) then  ! Use KF scheme <-- NOTE: THIS MODIFIES entr too!!
        detr = entr * (1-chi_c)**2.
        entr = entr * max(ER_KF_min, chi_c**2.)
    end if
    detr = DR_frac * detr
    
end subroutine compute_ent_det


! compute_zt_upd: find the highest level where buoy > 0.0, for detrainment profile (different from SCM !!)
subroutine compute_zt_upd (buoy, zfull, nlev, nupd, zt) 

    real   , intent(in) ,   dimension(:,:,:,:) :: buoy
    real   , intent(in) ,   dimension(:,:,:)   :: zfull
    integer, intent(in)                        :: nlev, nupd
    real   , intent(out),   dimension(:,:,:)   :: zt
    
    integer :: ilev,iupd, nx, ny
    logical, dimension(size(zt,1),size(zt,2),size(zt,3)) :: zt_found
    nx = size(zt,1)
    ny = size(zt,2)
    
    zt = 0.0
    zt_found = .false.
    
    do iupd = 1, nupd
        do ilev = 1, nlev
            where ( (zt_found(:,:,iupd) == .false.) .and. ( buoy(:,:,iupd,ilev) .gt. 0.0))
               zt_found(:,:,iupd) = .true.
               zt(:,:,iupd) = zfull(:,:,ilev) + 100.0
            end where
        end do
    end do

end subroutine compute_zt_upd


! compute_upd_buoy: Compute buoyancy of updraft
subroutine compute_upd_buoy (  thl_upd,  qt_upd,  ql_upd,  a_upd,   &
                               thl_env,  qt_env,  ql_env,  a_env,   & 
                               pfull, buoy, thv_upd, thv_env)  ! Saturation adjustment ... May be time-consuming

    real   , intent(in) ,   dimension(:,:,:,:) :: thl_upd,  qt_upd,  ql_upd,  a_upd                 
    real   , intent(in) ,   dimension(:,:,:)   :: thl_env,  qt_env,  ql_env,  a_env
    real   , intent(in) ,   dimension(:,:,:)   :: pfull
    real   , intent(out),   dimension(:,:,:,:) :: buoy, thv_upd
    real   , intent(out),   dimension(:,:,:)   :: thv_env
    
    real  , dimension(size(thl_upd,1),size(thl_upd,2),size(thl_upd,3),size(thl_upd,4)) :: rho_upd, thv_upd_tmp
    real  , dimension(size(thl_upd,1),size(thl_upd,2),size(thl_upd,3)) :: rho_env, thv_updsum, thv_grid
    
    
    
    integer :: nupd, nlev, iupd, ilev
    nlev = size(thl_upd,3)
    nupd = size(thl_upd,4)
    
    call compute_thv   (thl_env, qt_env, ql_env, pfull, thv_env, rho_env)     ! NOTE: rho_env is dummy and should not be used!
    call compute_thv_4d(thl_upd, qt_upd, ql_upd, pfull, thv_upd_tmp, rho_upd) ! NOTE: rho_upd is dummy and should not be used!
    
    call compute_upd_thv_correction(thl_upd,  qt_upd,  ql_upd,  a_upd,  thv_upd_tmp, & 
                                     pfull, thv_upd)
                                     
    ! CHANGED by ZTAN 09/25/2017: thv_g is now used to calculate buoyancy.
    thv_updsum = 0.
    do iupd = 1, nupd
        thv_updsum  = thv_updsum + a_upd(:,:,:,iupd) * thv_upd(:,:,:,iupd) * updraft_rescale
        ! Note: rescale the updraft fraction, but not rescaling total mass flux
        ! i.e., au -> au * upd_res wu -> wu / upd_res
    end do
    thv_grid = thv_env * a_env + thv_updsum
    
    do iupd = 1, nupd
        buoy(:,:,:,iupd) = Grav * (thv_upd(:,:,:,iupd)/thv_grid(:,:,:)-1.)
    end do
    
    ! NOTE by ZTAN 09/25/2017: Old formulation incorrectly used thv_e.
    ! do iupd = 1, nupd
    !     buoy(:,:,:,iupd) = Grav * (thv_upd(:,:,:,iupd)/thv_env(:,:,:)-1.)
    ! end do
    
end subroutine compute_upd_buoy


! compute_upd_thv_correction: correcting the updraft buoyancy to allow for overshooting its previous top
! Note by ZTAN 09/25/2017 -- Switch to an implicit formulation (b of the three mixing components) ?? 
subroutine compute_upd_thv_correction(thl_upd,  qt_upd,  ql_upd,  a_upd,  thv_upd_tmp, & 
                               pfull, thv_upd)
                               
    real   , intent(in) ,   dimension(:,:,:,:) :: thl_upd,  qt_upd,  ql_upd,  a_upd, thv_upd_tmp  
    real   , intent(in) ,   dimension(:,:,:)   :: pfull
    real   , intent(out),   dimension(:,:,:,:) :: thv_upd
    
    real  , dimension(size(thl_upd,1),size(thl_upd,2),size(thl_upd,3),size(thl_upd,4)) :: rho_upd, qln_upd
    real  , dimension(size(thl_upd,1),size(thl_upd,2),size(thl_upd,3)) :: a_updsum
    integer :: nupd, nlev, iupd, ilev
    
    nlev = size(thl_upd,3)
    nupd = size(thl_upd,4)
    a_updsum = 0.
    
    if (correct_buoy == 1) then  ! use buoyancy at one level below, when none of the updrafts exist
        do iupd = 1, nupd
            a_updsum = a_updsum + a_upd(:,:,:,iupd)
            thv_upd (:,:,nlev,iupd) = thv_upd_tmp (:,:,nlev,iupd)
            thv_upd (:,:,1:nlev-1,iupd) = thv_upd_tmp(:,:,2:nlev,iupd) ! use 1-level below by default
        end do
        
        do iupd = 1, nupd 
            do ilev = nlev-1, 1, -1
                where (a_updsum(:,:,ilev) > 0.0)   ! At least one updraft at this level, switch back to actual value             
                    thv_upd(:,:,ilev,iupd) = thv_upd_tmp(:,:,ilev,iupd)
                end where
            end do
        end do
        
    elseif (correct_buoy == 2) then  ! use buoyancy at one level below and with saturation adjustment, when none of the updrafts exist
        do iupd = 1, nupd
            a_updsum = a_updsum + a_upd(:,:,:,iupd)
            thv_upd (:,:,nlev,iupd) = thv_upd_tmp (:,:,nlev,iupd)
            qln_upd(:,:,:,iupd) = ql_upd(:,:,:,iupd)
            call compute_thv_checkcond( &
                thl_upd(:,:,2:nlev,iupd), qt_upd(:,:,2:nlev,iupd), qln_upd(:,:,2:nlev,iupd), pfull(:,:,1:nlev-1), &
                thv_upd(:,:,1:nlev-1,iupd), rho_upd(:,:,1:nlev-1,iupd))
        end do
        
        do iupd = 1, nupd 
            do ilev = nlev-1, 1, -1
                where (a_updsum(:,:,ilev) > 0.0)   ! At least one updraft at this level, switch back to actual value             
                    thv_upd(:,:,ilev,iupd) = thv_upd_tmp(:,:,ilev,iupd)
                end where
            end do
        end do
        
    elseif (correct_buoy == 11) then  ! use buoyancy at one level below, when THIS updraft does not exist
        do iupd = 1, nupd
            thv_upd (:,:,nlev,iupd) = thv_upd_tmp (:,:,nlev,iupd)
            thv_upd (:,:,1:nlev-1,iupd) = thv_upd_tmp(:,:,2:nlev,iupd) ! use 1-level below by default
        end do
        
        where (a_upd (:,:,:,:) > 0.0) ! Updraft exists at this level, switch back to actual value   
            thv_upd (:,:,:,:) = thv_upd_tmp(:,:,:,:)
        end where
    
    elseif (correct_buoy == 12) then  ! use buoyancy at one level below and with saturation adjustment, when THIS updraft does not exist
        do iupd = 1, nupd
            thv_upd (:,:,nlev,iupd) = thv_upd_tmp (:,:,nlev,iupd)
            qln_upd(:,:,:,iupd) = ql_upd(:,:,:,iupd)
            call compute_thv_checkcond( &
                thl_upd(:,:,2:nlev,iupd), qt_upd(:,:,2:nlev,iupd), qln_upd(:,:,2:nlev,iupd), pfull(:,:,1:nlev-1), &
                thv_upd(:,:,1:nlev-1,iupd), rho_upd(:,:,1:nlev-1,iupd))
        end do

        where (a_upd (:,:,:,:) > 0.0) ! Updraft exists at this level, switch back to actual value   
            thv_upd (:,:,:,:) = thv_upd_tmp(:,:,:,:)
        end where
        
    else
        thv_upd(:,:,:,:) = thv_upd_tmp(:,:,:,:)
    end if

end subroutine compute_upd_thv_correction


! compute_chi_c: compute the critical mixing ratio between updraft and environment air for buoyancy reversal:
! chi_c = 0 => all mixture will sink; chi_c = 1 => all mixture will rise
! TO_DO_082016: may not work in current form -- NEED TO CHECK THE FORMULA AGAIN !!! 
subroutine compute_chi_c(thl_upd,  qt_upd,  ql_upd,  thv_upd, thl_env,  qt_env,  ql_env,  thv_env, pfull, chi_c)

    real   , intent(in) ,   dimension(:,:,:,:) :: thl_upd,  qt_upd,  ql_upd,  thv_upd             
    real   , intent(in) ,   dimension(:,:,:)   :: thl_env,  qt_env,  ql_env,  thv_env
    real   , intent(in) ,   dimension(:,:,:)   :: pfull
    real   , intent(out),   dimension(:,:,:,:) :: chi_c   
    
    real  , dimension(size(thl_upd,1),size(thl_upd,2),size(thl_upd,3)) :: pip, qsat_env, dqsat_env, rh_env, qln_env
    real  , dimension(size(thl_upd,1),size(thl_upd,2),size(thl_upd,3),size(thl_upd,4)) :: &
            qsat_upd, dqsat_upd, alp_upd, bet_upd, gam_upd, del_upd, dif_thv, qln_upd
    integer :: nupd, nlev, iupd, ilev, ix, iy, nx, ny
    
    nx = size(thl_upd,1)
    ny = size(thl_upd,2)
    nlev = size(thl_upd,3)
    nupd = size(thl_upd,4)
    
    pip = (pfull(:,:,:)/pstd_mks)**kappa
    qln_env(:,:,:) = ql_env(:,:,:)
    call compute_ql(thl_env*pip, qt_env, qln_env, pfull, qsat_env, dqsat_env)
    do iupd = 1, nupd
        qln_upd(:,:,:,iupd) = ql_upd(:,:,:,iupd)
        call compute_ql(thl_upd(:,:,:,iupd)*pip, qt_upd(:,:,:,iupd), qln_upd(:,:,:,iupd), &
                        pfull, qsat_upd(:,:,:,iupd), dqsat_upd(:,:,:,iupd))
        alp_upd(:,:,:,iupd) = Cp_air/HLv * pip * thl_upd(:,:,:,iupd)
    end do
    rh_env = qt_env/qsat_env
    gam_upd = HLv/Cp_air * dqsat_upd
    bet_upd = 1/(1+gam_upd)*(1+(1+d608)*gam_upd*alp_upd)
    
    do iupd = 1, nupd
        del_upd(:,:,:,iupd) = qsat_env*(bet_upd(:,:,:,iupd)-alp_upd(:,:,:,iupd))*(1-rh_env) & 
                            - alp_upd(:,:,:,iupd)*ql_upd(:,:,:,iupd)
        dif_thv(:,:,:,iupd) = thv_upd(:,:,:,iupd) - thv_env(:,:,:)
    end do
    
    chi_c = 0.
    do iupd = 1, nupd
        do ilev = 1, nlev
            do iy = 1, ny
                do ix = 1, nx
                    if     (dif_thv(ix,iy,ilev,iupd) .le. 0.0) then
                        chi_c(ix,iy,ilev,iupd) = 0.0
                    elseif (del_upd(ix,iy,ilev,iupd) .le. 0.0) then
                        chi_c(ix,iy,ilev,iupd) = 1.0
                    elseif ( ql_upd(ix,iy,ilev,iupd) .eq. 0.0) then
                        chi_c(ix,iy,ilev,iupd) = 1.0
                    else
                        chi_c(ix,iy,ilev,iupd) = Cp_air/HLv*pip(ix,iy,ilev)*dif_thv(ix,iy,ilev,iupd)/del_upd(ix,iy,ilev,iupd)  
                    end if
                end do
            end do
        end do
    end do
    chi_c = max(0.0, min(chi_c, 1.0))

end subroutine compute_chi_c


! ------------------------- !
!  Parametrizations (COEF)  !
! ------------------------- !

! compute_surface_tke: compute surface TKE based on wstar and zstar
subroutine compute_surface_tke(ustar, wstar, tke_surf)
    real  , intent(in) , dimension(:,:)   :: ustar, wstar
    real  , intent(out), dimension(:,:)   :: tke_surf
    real  :: small_val
    
    small_val = 1.0e-8
    
    tke_surf(:,:) = max( (3.75*(ustar(:,:)**2.) + 0.2*(wstar(:,:)**2.)), small_val)

end subroutine compute_surface_tke

! compute_surface_var: compute surface scalar covariance for scalars #1 and #2
subroutine compute_surface_var(scl_flux_1, scl_flux_2, ustar, z_ll, oblength, scl_var)
    real  , intent(in) , dimension(:,:)   :: scl_flux_1, scl_flux_2, ustar, z_ll, oblength
    real  , intent(out), dimension(:,:)   :: scl_var
    
    real  , dimension(size(ustar,1),size(ustar,2)) :: scl_star_1, scl_star_2
    real  :: small_val
    
    small_val = 1.0e-8
    
    scl_star_1(:,:) = scl_flux_1(:,:) / max(ustar(:,:), small_val) ! Question by ZTAN 09/25/2017: why a negative sign in SCAMPy ??
    scl_star_2(:,:) = scl_flux_2(:,:) / max(ustar(:,:), small_val) 
    
    where (oblength(:,:) < 0.0)
        scl_var(:,:) = 4.0 * scl_star_1(:,:) * scl_star_2(:,:) * &
                       (1.0 - 8.3 * z_ll(:,:)/oblength(:,:))**(-2.0/3.0)
    elsewhere
        scl_var(:,:) = 4.0 * scl_star_1(:,:) * scl_star_2(:,:)
    end where
    
end subroutine compute_surface_var

! compute_surface_std: compute signed surface scalar standard deviation for a scalar, where sign is same as flux
subroutine compute_surface_std(scl_flux, ustar, z_ll, oblength, scl_std)
    real  , intent(in) , dimension(:,:)   :: scl_flux, ustar, z_ll, oblength
    real  , intent(out), dimension(:,:)   :: scl_std
    
    real  , dimension(size(ustar,1),size(ustar,2)) :: scl_star
    real  :: small_val
    
    small_val = 1.0e-8
    
    scl_star(:,:) = scl_flux(:,:) / max(ustar(:,:), small_val)
    
    where (oblength(:,:) < 0.0)
        scl_std(:,:) = 2.0 * scl_star(:,:) * (1.0 - 8.3 * z_ll(:,:)/oblength(:,:))**(-1.0/3.0)
    elsewhere
        scl_std(:,:) = 2.0 * scl_star(:,:)
    end where
    
end subroutine compute_surface_std

! compute_zstar_wstar_obl: compute surface conditions and BL depth
subroutine compute_zstar_wstar_obl(Time, tin, qin, liqin, thl, qt, thv, u, v,  &  ! grid_mean
        ustar, bstar, flux_t, flux_q,                             &  ! surf_flux
        pfull, zfull,  & 
        zstar, wthv_surf, wthl_surf, wq_surf, wstar, oblength) ! ADD DIMENSIONS !!
        
    type(time_type), intent(in)           :: Time
    real  , intent(in) , dimension(:,:,:) :: tin, qin, liqin, thl, qt, thv, u, v, pfull, zfull
    real  , intent(in) , dimension(:,:)   :: ustar, bstar, flux_q, flux_t
    real  , intent(out), dimension(:,:)   :: zstar, wthv_surf, wthl_surf, wq_surf, wstar, oblength

    real, dimension(size(tin,1),size(tin,2),size(tin,3)) :: svcp
    integer :: num_levels
    
    num_levels = size(tin,3)
    
    call compute_svcp(tin, qin, liqin, thl, pfull, zfull, svcp)
    if  (pbl_depth_opt == 0)  then   ! Use fixed PBL depth
        zstar = pbl_fixed_depth
    else 
        call pbl_depth(Time, svcp, u,v, zfull, ustar, bstar, zstar)  
    end if
    
    zstar = zstar * pbl_depth_res
    
    call cal_wstar(Time, tin(:,:, num_levels), qin(:,:, num_levels), liqin(:,:, num_levels),  &
                    pfull(:,:, num_levels), zstar, flux_t, flux_q,  &
                    wthv_surf, wthl_surf, wq_surf, wstar)   ! PENDING...

    where (wthv_surf == 0.0)
          oblength = 0.0!  If calculated wthv =0, give oblength =0.
    elsewhere
          oblength = -ustar**3.0 * thv(:,:, num_levels)/vonk/Grav/wthv_surf ! This is Eqn (8) in 2010 paper
    end where    

end subroutine compute_zstar_wstar_obl       


! cal_wstar: calculate w_star and surface scalar fluxes (only for use in compute_zstar_wstar_obl)
subroutine cal_wstar(Time, t_atm, q_atm, liq_atm, p_atm, zstar, flux_t, flux_q, wthv, wthl, wq, wstar)
   type(time_type), intent(in) :: Time
   real   , intent(in) , dimension(:,:)   :: t_atm, q_atm,liq_atm, p_atm, zstar, flux_q, flux_t
   real   , intent(out), dimension(:,:)   :: wthv, wthl, wq, wstar
   
   real  , dimension(size(zstar,1),size(zstar,2)) ::  rho, tv_atm, th_atm,thv_atm, thv_surf, p_ratio

     p_ratio = (p_atm / pstd_mks ) ** (-kappa)
     tv_atm  = t_atm  * (1.0 + d608*q_atm - liq_atm)     ! virtual temperature
     th_atm  = t_atm  * p_ratio  ! potential T, using pstd_mks as refernce
     thv_atm = tv_atm * p_ratio  ! virt. potential T, using pstd_mks as reference 

     rho = p_atm / (rdgas * tv_atm)   ! density  
     wthv = flux_t/cp_air/rho * (1 + d608 * q_atm) + flux_q/rho*d608*th_atm
     wthl = flux_t/cp_air/rho
     wq = flux_q/rho
     
     where (wthv .ge. 0.0)
         wstar = (zstar* Grav *wthv / thv_atm)**(1.0/3.0)  ! Refer to Garratt book
         ! This expression was problematic: change by ZTAN 02/07/11:  1/3-> 1.0/3.0
     elsewhere
         wstar = -(-zstar* Grav *wthv / thv_atm)**(1.0/3.0)
     end where

end subroutine cal_wstar


! compute_svcp: criteria for convective boundary layer depth (only for use in compute_zstar_wstar_obl)
subroutine compute_svcp(tin, qin, liqin, thl, pfull, zfull, svcp)
    real  , intent(in) , dimension(:,:,:) :: tin, qin, liqin, thl, pfull, zfull
    real  , intent(out), dimension(:,:,:) :: svcp
    real, dimension(size(tin,1),size(tin,2),size(tin,3)) :: tl, th


    ! ?1: virtual with water loading; ?2: virtual without water loading; ?3: dry
    ! 0?: dry static energy; 1?: liquid static energy; 2?: pot temp; 3?: liquid pot temp
        
    if     (pbl_depth_opt == 1)  then  ! virtual static energy (w/ water loading effect)
        svcp = tin * (1.0 + d608*qin - liqin) + grav/cp_air*zfull
        
    elseif (pbl_depth_opt == 2)  then  ! virtual static energy (w/o water loading effect)
        svcp = tin * (1.0 + d608*qin) + grav/cp_air*zfull 
        
    elseif (pbl_depth_opt == 3)  then  ! dry static energy
        svcp = tin + grav/cp_air*zfull 
        
    elseif (pbl_depth_opt == 11) then  ! liquid virtual static energy (w/ water loading effect)
        tl = thl * ((pfull/pstd_mks)**kappa)
        svcp = tl  * (1.0 + d608*qin - liqin) + grav/cp_air*zfull
        
    elseif (pbl_depth_opt == 12) then  ! liquid virtual static energy (w/o water loading effect)
        tl = thl * ((pfull/pstd_mks)**kappa)
        svcp = tl  * (1.0 + d608*qin) + grav/cp_air*zfull
        
    elseif (pbl_depth_opt == 13) then  ! liquid static energy
        tl = thl * ((pfull/pstd_mks)**kappa)
        svcp = tl + grav/cp_air*zfull 
        
    elseif (pbl_depth_opt == 21) then  ! virtual pot temp  (w/ water loading effect)
        th = tin * ((pstd_mks/pfull)**kappa)
        svcp = th * (1.0 + d608*qin - liqin)
        
    elseif (pbl_depth_opt == 22) then  ! virtual pot temp  (w/o water loading effect)
        th = tin * ((pstd_mks/pfull)**kappa)
        svcp = th * (1.0 + d608*qin)
        
    elseif (pbl_depth_opt == 23) then  ! dry pot temp
        th = tin * ((pstd_mks/pfull)**kappa)
        svcp = th
            
    elseif (pbl_depth_opt == 31) then  ! liquid virtual pot temp (w/ water loading effect)
        svcp = thl * (1.0 + d608*qin - liqin)
         
    elseif (pbl_depth_opt == 32) then  ! liquid virtual pot temp (w/o water loading effect)
        svcp = thl * (1.0 + d608*qin) 
    
    elseif (pbl_depth_opt == 33) then  ! liquid pot temp
        svcp = thl
    
    else ! Other options
        svcp = 0.0
    end if

end subroutine compute_svcp


! ------------------------- !
!  Utilities (3D/4D fields) !
! ------------------------- !


! cloud_mixing: compute thl and qt from t, q, ql
subroutine cloud_mixing(t, q, ql, pfull, thl, qt)

    real   , intent(in) , dimension(:,:,:) :: t, q, ql, pfull
    real   , intent(out), dimension(:,:,:) :: thl, qt

    thl(:,:,:)= (t(:,:,:)-HLv/Cp_air*ql(:,:,:))*((pstd_mks/pfull(:,:,:))**kappa)
    qt(:,:,:) =  q(:,:,:)+ql(:,:,:)

end subroutine cloud_mixing


! cloud_decompose: compute t, q, ql from thl and qt
subroutine cloud_decompose(thl, qt, ql, pfull, t, q)

    real   , intent(in) , dimension(:,:,:) :: thl, qt, ql, pfull
    real   , intent(out), dimension(:,:,:) :: t, q

    q(:,:,:) = qt(:,:,:) - ql(:,:,:)
    t(:,:,:) = thl(:,:,:)/((pstd_mks/pfull(:,:,:))**kappa)+HLv/Cp_air*ql(:,:,:)

end subroutine cloud_decompose


! cloud_decompose_2d: compute t, q, ql from thl and qt
subroutine cloud_decompose_2d(thl, qt, ql, pfull, t, q)

    real   , intent(in) , dimension(:,:) :: thl, qt, ql, pfull
    real   , intent(out), dimension(:,:) :: t, q

    q(:,:) = qt(:,:) - ql(:,:)
    t(:,:) = thl(:,:)/((pstd_mks/pfull(:,:))**kappa)+HLv/Cp_air*ql(:,:)

end subroutine cloud_decompose_2d

! compute_thv: compute thv from thl, qt, ql (without saturation adjustment)
subroutine compute_thv(thl, qt, ql, pfull, thv, rho)

    real   , intent(in) , dimension(:,:,:) :: thl, qt, ql, pfull
    real   , intent(out), dimension(:,:,:) :: thv, rho

    real  , dimension(size(thl,1),size(thl,2),size(thl,3)) :: t, q
   
    call cloud_decompose(thl, qt, ql, pfull, t, q)
    thv(:,:,:) = t(:,:,:) * (1.0 + d608*q(:,:,:) - ql(:,:,:)) * ((pstd_mks/pfull(:,:,:))**kappa)
    rho(:,:,:) = pfull(:,:,:)/rdgas/(t(:,:,:) * (1.0 + d608*q(:,:,:) - ql(:,:,:)))
   
end subroutine compute_thv

! compute_thv_2d: compute thv from thl, qt, ql (without saturation adjustment)
subroutine compute_thv_2d(thl, qt, ql, pfull, thv, rho)

    real   , intent(in) , dimension(:,:) :: thl, qt, ql, pfull
    real   , intent(out), dimension(:,:) :: thv, rho

    real  , dimension(size(thl,1),size(thl,2)) :: t, q
   
    call cloud_decompose_2d(thl, qt, ql, pfull, t, q)
    thv(:,:) = t(:,:) * (1.0 + d608*q(:,:) - ql(:,:)) * ((pstd_mks/pfull(:,:))**kappa)
    rho(:,:) = pfull(:,:)/rdgas/(t(:,:) * (1.0 + d608*q(:,:) - ql(:,:)))
   
end subroutine compute_thv_2d

! compute_thv_4d: compute thv from thl, qt, ql (without saturation adjustment - 4d fields)
subroutine compute_thv_4d(thl, qt, ql, pfull, thv, rho)   ! similar to the trick in mo_profile_2d

    real   , intent(in) , dimension(:,:,:,:) :: thl, qt, ql
    real   , intent(in) , dimension(:,:,:)   :: pfull
    real   , intent(out), dimension(:,:,:,:) :: thv, rho

    integer :: nupd, iupd
    nupd = size(thl,4)
    do iupd = 1, nupd
        call compute_thv(thl(:,:,:,iupd), qt(:,:,:,iupd), ql(:,:,:,iupd), pfull, & 
                         thv(:,:,:,iupd), rho(:,:,:,iupd))
    end do
   
end subroutine compute_thv_4d



! compute_thv_checkcond: compute thv from thl, qt, ql (with saturation adjustment)
subroutine compute_thv_checkcond(thl, qt, ql, pfull, thv, rho)

    real   , intent(in) ,    dimension(:,:,:) :: thl, qt, pfull
    real   , intent(inout) , dimension(:,:,:) :: ql
    real   , intent(out),    dimension(:,:,:) :: thv, rho

    real  , dimension(size(thl,1),size(thl,2),size(thl,3)) :: t, q, tl, qsat, dqsat
   
    tl(:,:,:) = thl(:,:,:)/((pstd_mks/pfull(:,:,:))**kappa)
    call compute_ql(tl, qt, ql, pfull, qsat, dqsat)
    call cloud_decompose(thl, qt, ql, pfull, t, q)
    thv(:,:,:) = t(:,:,:) * (1.0 + d608*q(:,:,:) - ql(:,:,:)) * ((pstd_mks/pfull(:,:,:))**kappa)
    rho(:,:,:) = pfull(:,:,:)/rdgas/(t(:,:,:) * (1.0 + d608*q(:,:,:) - ql(:,:,:)))
   
end subroutine compute_thv_checkcond



! ! compute_thv_checkcond_0d: compute thv from thl, qt, ql (with saturation adjustment, 0d)
!subroutine compute_thv_checkcond_0d(thl, qt, ql, pfull, thv, rho)
!
!    real   , intent(in) ,  :: thl, qt, pfull
!    real   , intent(inout) :: ql
!    real   , intent(out),  :: thv, rho
!
!    real  , dimension(1,1,1) :: thl_3d, qt_3d, ql_3d, thv_3d, rho_3d
!    
!    thl_3d(1,1,1) = thl
!    qt_3d(1,1,1)  = qt
!    ql_3d(1,1,1)  = ql
!    pfull_3d(1,1,1) = pfull
!    
!    call compute_thv_checkcond(thl_3d, qt_3d, ql_3d, pfull_3d, thv_3d, rho_3d)
!    
!    thv = thv_3d(1,1,1)
!    rho = rho_3d(1,1,1)
!    ql  = ql_3d(1,1,1)
!   
!end subroutine compute_thv_checkcond_0d

! compute_qsat: similar to compute_ql, but from T, and no adjustment
subroutine compute_qsat(t, pfull, qsat, dqsat)

    real   , intent(in) ,    dimension(:,:,:) :: t, pfull
    real   , intent(out),    dimension(:,:,:) :: qsat, dqsat

    real  , dimension(size(t,1),size(t,2),size(t,3)) :: esat, desat
    real :: small = 1.0e-10
    
    call  escomp (t, esat)
    call descomp (t, desat)
       
        ! Saturation specific humidity
     qsat(:,:,:) = d622*esat(:,:,:)/max(pfull(:,:,:)-d378*esat(:,:,:), small)
    dqsat(:,:,:) = d622*pfull(:,:,:)*desat(:,:,:)/(max(pfull(:,:,:)-d378*esat(:,:,:), small)**2.0)
   
end subroutine compute_qsat


! compute_ql: compute ql from tl, qt (i.e., saturation adjustment)
subroutine compute_ql(tl, qt, ql, pfull, qsat, dqsat)

    real   , intent(in) ,    dimension(:,:,:) :: tl, qt, pfull
    real   , intent(inout) , dimension(:,:,:) :: ql
    real   , intent(out),    dimension(:,:,:) :: qsat, dqsat

    real  , dimension(size(tl,1),size(tl,2),size(tl,3)) :: esat, desat, gammaU
    real :: small = 1.0e-10
    integer :: iter
   
    do iter = 1, maxiter
        call  escomp (tl+ HLv*ql/Cp_air, esat)
        call descomp (tl+ HLv*ql/Cp_air, desat)
       
        ! Saturation specific humidity
         qsat(:,:,:) = d622*esat(:,:,:)/max(pfull(:,:,:)-d378*esat(:,:,:), small)
        dqsat(:,:,:) = d622*pfull(:,:,:)*desat(:,:,:)/(max(pfull(:,:,:)-d378*esat(:,:,:), small)**2.0)
        ! when p<es, it's impossible to saturate -> qs = 1 ! How to define dqs when qs = inf?
       
        gammaU(:,:,:) = HLv*HLv/(Cp_air*rvgas*((tl(:,:,:)+ HLv*ql(:,:,:)/Cp_air)**2.0))*qsat(:,:,:)
        ! Saturation excess (Cheinet and Teixeira 2003)
       
        ql(:,:,:) = (qt(:,:,:)-qsat(:,:,:)-ql(:,:,:))/(1.0+gammaU(:,:,:)) + ql(:,:,:)
    end do
   
    ql(:,:,:) = max(ql(:,:,:), 0.0)
   
end subroutine compute_ql

! compute_ql_4d: compute ql from tl, qt (i.e., saturation adjustment - 4d fields)
subroutine compute_ql_4d(tl, qt, ql, pfull, qsat, dqsat)   ! similar to the trick in mo_profile_2d


    real   , intent(in) ,    dimension(:,:,:,:) :: tl, qt
    real   , intent(in) ,    dimension(:,:,:)   :: pfull
    real   , intent(inout) , dimension(:,:,:,:) :: ql
    real   , intent(out),    dimension(:,:,:,:) :: qsat, dqsat
    

    integer :: nupd, iupd
    nupd = size(tl,4)
    do iupd = 1, nupd
        call compute_ql(tl(:,:,:,iupd), qt(:,:,:,iupd), ql(:,:,:,iupd), pfull, & 
                        qsat(:,:,:,iupd), dqsat(:,:,:,iupd))
    end do
   
end subroutine compute_ql_4d


! compute_ql_2d: compute ql from tl, qt (i.e., saturation adjustment - 2d fields)
subroutine compute_ql_2d(tl, qt, ql, pfull, qsat, dqsat)   ! similar to the trick in mo_profile_0d


    real   , intent(in) ,    dimension(:,:) :: tl, qt
    real   , intent(in) ,    dimension(:,:) :: pfull
    real   , intent(inout) , dimension(:,:) :: ql
    real   , intent(out),    dimension(:,:) :: qsat, dqsat
    
    real  , dimension(size(tl,1),size(tl,2),1) :: tl_3d, qt_3d, ql_3d, pfull_3d, qsat_3d, dqsat_3d
    
    tl_3d(:,:,1) = tl(:,:)
    qt_3d(:,:,1) = qt(:,:)
    ql_3d(:,:,1) = ql(:,:)
    pfull_3d(:,:,1) = pfull(:,:)
    
    call compute_ql(tl_3d, qt_3d, ql_3d, pfull_3d, qsat_3d, dqsat_3d)
    
    qsat(:,:) = qsat_3d(:,:,1)
    dqsat(:,:) = dqsat_3d(:,:,1)
    ql(:,:) = ql_3d(:,:,1)
   
end subroutine compute_ql_2d



! envupd_decompose: compute environment condition from grid-mean and updraft conditions
subroutine envupd_decompose(thl_grid, qt_grid, ql_grid, w_grid,        u_grid, v_grid, tke_grid, &
                            thl_upd,  qt_upd,  ql_upd,  w_upd,  a_upd,                           & 
                            thl_env,  qt_env,  ql_env,  w_env,  a_env, u_env,  v_env,  tke_env)

    real   , intent(in) , dimension(:,:,:)   :: thl_grid, qt_grid, ql_grid, w_grid, &
                                                u_grid, v_grid, tke_grid
    real   , intent(in) , dimension(:,:,:,:)    :: thl_upd, w_upd, a_upd
    real   , intent(inout) , dimension(:,:,:,:) :: qt_upd, ql_upd
    real   , intent(out), dimension(:,:,:)   :: thl_env, qt_env, ql_env, w_env, a_env, &
                                                u_env, v_env, tke_env

    real  , dimension(size(thl_grid,1),size(thl_grid,2),size(thl_grid,3)) :: &
        thl_updsum, qt_updsum, ql_updsum, w_updsum, a_updsum, qt_updsum_corr, ql_updsum_corr
        ! Note these values are the weighted sum of (au_i * fu_i) and not divided by sum(au_i)
    
    integer :: nupd, iupd
    nupd = size(thl_upd,4)
    
    call area_sum_updr(thl_upd, qt_upd, ql_upd, w_upd, a_upd, &
         thl_updsum, qt_updsum, ql_updsum, w_updsum, a_updsum)

    a_env   = 1.0 - a_updsum
    thl_env = (thl_grid - thl_updsum )/a_env
    w_env   = (w_grid   - w_updsum   )/a_env

    ! ADDED ZTAN 09/12/2017: Limit qt_env and ql_env to be both greater than zero; 
    !                        qt_upd and ql_upd may be adjusted down to satisfy this constraint.
    ! NOTE 1: Are the tendency terms still exact after this correction?
    ! NOTE 2: this correction satisfies: qt_env*a_env + qt_upd*a_upd(*res_fac) = qt_grid, 
    !                                and ql_env*a_env + ql_upd*a_upd(*res_fac) = ql_grid.
    !         Therefore, this correction does not induce a grid-mean tendency.
    
    qt_updsum_corr = min(qt_updsum, qt_grid)
    ql_updsum_corr = min(ql_updsum, ql_grid)
    qt_env  = (qt_grid  - qt_updsum_corr  )/a_env
    ql_env  = (ql_grid  - ql_updsum_corr  )/a_env
    
    do iupd = 1, nupd
       where (qt_updsum(:,:,:) .gt. 0.0)
           qt_upd(:,:,:,iupd) = qt_upd(:,:,:,iupd)*qt_updsum_corr(:,:,:)/qt_updsum(:,:,:)
       endwhere
       where (ql_updsum(:,:,:) .gt. 0.0)
           ql_upd(:,:,:,iupd) = ql_upd(:,:,:,iupd)*ql_updsum_corr(:,:,:)/ql_updsum(:,:,:)
       endwhere
    end do
    
    ! No decomposition for u and v
    u_env = u_grid
    v_env = v_grid
    
    ! unsure about TKE decomposition
    tke_env = tke_grid / a_env

end subroutine envupd_decompose

subroutine envupd_combine  (thl_upd,  qt_upd,  ql_upd,  w_upd,  a_upd,                           & 
                            thl_env,  qt_env,  ql_env,  w_env,         u_env,  v_env,  tke_env,  &
                            thl_grid, qt_grid, ql_grid, w_grid,        u_grid, v_grid, tke_grid)

    real   , intent(out), dimension(:,:,:)   :: thl_grid, qt_grid, ql_grid, w_grid, &
                                                u_grid, v_grid, tke_grid
    real   , intent(in) , dimension(:,:,:,:) :: thl_upd, qt_upd, ql_upd, w_upd, a_upd
    real   , intent(in) , dimension(:,:,:)   :: thl_env, qt_env, ql_env, w_env, &
                                                u_env, v_env, tke_env

    real  , dimension(size(thl_grid,1),size(thl_grid,2),size(thl_grid,3)) :: &
        thl_updsum, qt_updsum, ql_updsum, w_updsum, a_updsum, a_env
        ! Note these values are the weighted sum of (au_i * fu_i) and not divided by sum(au_i)
      
    call area_sum_updr(thl_upd, qt_upd, ql_upd, w_upd, a_upd, &
         thl_updsum, qt_updsum, ql_updsum, w_updsum, a_updsum)

    a_env   = 1.0 - a_updsum
    thl_grid = thl_env * a_env + thl_updsum
    qt_grid  = qt_env  * a_env + qt_updsum
    ql_grid  = ql_env  * a_env + ql_updsum
    w_grid   = w_env   * a_env + w_updsum

    ! No decomposition for u and v
    u_grid = u_env
    v_grid = v_env
    
    ! unsure about TKE decomposition
    tke_grid = tke_env * a_env

end subroutine envupd_combine

subroutine do_phiu_corr( thl_upd,  qt_upd,  ql_upd,  a_upd,  & ! updr (inout)
                         thl_uold, qt_uold, ql_uold,         & ! updr_old (in)
                         thl_env,  qt_env,  ql_env )           ! envr (in) -- Use value at dt_real

    real   , intent(in)   , dimension(:,:,:,:) :: thl_uold, qt_uold, ql_uold
    real   , intent(in)   , dimension(:,:,:)   :: thl_env,  qt_env,  ql_env 
    real   , intent(inout), dimension(:,:,:,:) :: thl_upd,  qt_upd,  ql_upd, a_upd

    integer :: nupd, iupd, nlev, ilev
    
    nlev = size(thl_upd, 3)
    nupd = size(thl_upd, 4)
    
    if (phi_u_corr_opt == 0) then
        do iupd = 1, nupd
           where (a_upd (:,:,:,iupd) < 1.e-4)
               thl_upd(:,:,:,iupd) = thl_env(:,:,:)
               qt_upd (:,:,:,iupd) = qt_env (:,:,:)
               ql_upd (:,:,:,iupd) = ql_env (:,:,:)
               a_upd  (:,:,:,iupd) = 0.0
           end where
        end do
        
    elseif (phi_u_corr_opt == 1) then
        do iupd = 1, nupd
           where (a_upd (:,:,:,iupd) < 1.e-4)
               thl_upd(:,:,:,iupd) = thl_uold(:,:,:,iupd)
               qt_upd (:,:,:,iupd) = qt_uold (:,:,:,iupd)
               ql_upd (:,:,:,iupd) = ql_uold (:,:,:,iupd)
               a_upd  (:,:,:,iupd) = 0.0
           end where
        end do
        
    elseif (phi_u_corr_opt == 2) then
        do iupd = 1, nupd
           ! Fill in the lowest level with env value 
           ! This is overwritten in the next step anyway, and it should not occur with fixed au_sfc.
           where (a_upd (:,:,nlev,iupd) < 1.e-4)
               thl_upd(:,:,nlev,iupd) = thl_env(:,:,nlev)
               qt_upd (:,:,nlev,iupd) = qt_env (:,:,nlev)
               ql_upd (:,:,nlev,iupd) = ql_env (:,:,nlev)
               a_upd  (:,:,nlev,iupd) = 0.0
           end where
           
           ! For other levels, fill in a blend of lower-level updraft and environment values 
           do ilev = nlev-1, 1, -1
               where (a_upd (:,:,ilev,iupd) < 1.e-4)
                   thl_upd(:,:,ilev,iupd) = thl_env(:,:,ilev) *.3 + thl_upd(:,:,ilev+1,iupd) *.7
                   qt_upd (:,:,ilev,iupd) = qt_env (:,:,ilev) *.3 + qt_upd (:,:,ilev+1,iupd) *.7
                   ql_upd (:,:,ilev,iupd) = ql_env (:,:,ilev) *.3 + ql_upd (:,:,ilev+1,iupd) *.7
                   a_upd  (:,:,ilev,iupd) = 0.0
               end where
            end do
        end do
    
    endif
    
end subroutine do_phiu_corr


subroutine trans_dphi_to_phiu (thl_upd, qt_upd, ql_upd, a_upd, thl_grid, qt_grid, ql_grid)

    real   , intent(in)   , dimension(:,:,:,:) :: a_upd
    real   , intent(in)   , dimension(:,:,:)   :: thl_grid, qt_grid
    real   , intent(inout), dimension(:,:,:)   :: ql_grid       ! ZTAN 09/11/2017 -- changed to inout
    real   , intent(inout), dimension(:,:,:,:) :: thl_upd, qt_upd, ql_upd
    
    integer :: nupd, iupd
    real  , dimension(size(thl_grid,1),size(thl_grid,2),size(thl_grid,3)) :: &
            a_updsum, thl_updsum, qt_updsum, ql_updsum  ! Sum of anomalies, not actual values
    
    nupd = size(thl_upd, 4)
    thl_updsum = 0.; qt_updsum = 0.; ql_updsum = 0.; a_updsum= 0.
    
    do iupd = 1, nupd
          a_updsum  =   a_updsum + a_upd(:,:,:,iupd)                       * updraft_rescale
        thl_updsum  = thl_updsum + a_upd(:,:,:,iupd) * thl_upd(:,:,:,iupd) * updraft_rescale 
         qt_updsum  =  qt_updsum + a_upd(:,:,:,iupd) *  qt_upd(:,:,:,iupd) * updraft_rescale 
         ql_updsum  =  ql_updsum + a_upd(:,:,:,iupd) *  ql_upd(:,:,:,iupd) * updraft_rescale 
    end do 
    
    ! CHANGED by ZTAN 09/10/2017: Add limiter to ql_grid 
    !   Force ql_grid >= ql_u * a_u, while keeping thl_grid and qt_grid unchanged
    ql_grid = max(ql_grid, ql_updsum)
    
    ! CHANGED by ZTAN 09/10/2017: Previous formula does not work for multiple updrafts
    ! do iupd = 1, nupd
    !    thl_upd(:,:,:,iupd) = thl_grid (:,:,:) + (1. - a_updsum(:,:,:)) * thl_upd(:,:,:,iupd) 
    !    qt_upd (:,:,:,iupd) = qt_grid  (:,:,:) + (1. - a_updsum(:,:,:)) * qt_upd (:,:,:,iupd) 
    !    ql_upd (:,:,:,iupd) = ql_grid  (:,:,:) + (1. - a_updsum(:,:,:)) * ql_upd (:,:,:,iupd) 
    ! end do 
    
    do iupd = 1, nupd
       thl_upd(:,:,:,iupd) = thl_grid (:,:,:) - thl_updsum(:,:,:) + thl_upd(:,:,:,iupd) 
       qt_upd (:,:,:,iupd) = qt_grid  (:,:,:) - qt_updsum (:,:,:) + qt_upd (:,:,:,iupd) 
       ql_upd (:,:,:,iupd) = ql_grid  (:,:,:) - ql_updsum (:,:,:) + ql_upd (:,:,:,iupd) 
       ! Note: Sum of the first two terms on the r.h.s. are thl_e, qt_e, ql_e, respectively.
       ! Because: f_grid - f_updsum = f_e*(1-a_u*res) + f_u*a_u*res - (f_u - f_e)*a_u*res = f_e.
    end do 
    
end subroutine trans_dphi_to_phiu


subroutine area_sum_updr(thl_upd, qt_upd, ql_upd, w_upd, a_upd, &
         thl_updsum, qt_updsum, ql_updsum, w_updsum, a_updsum)
         
    real   , intent(in) , dimension(:,:,:,:) :: thl_upd, qt_upd, ql_upd, w_upd, a_upd
    real   , intent(out), dimension(:,:,:)   :: &
                                    thl_updsum, qt_updsum, ql_updsum, w_updsum, a_updsum

    integer :: nupd, iupd
    
    nupd = size(thl_upd, 4)
    
    thl_updsum = 0.; qt_updsum = 0.; ql_updsum = 0.; w_updsum = 0.; a_updsum= 0.
    do iupd = 1, nupd
          a_updsum  =   a_updsum + a_upd(:,:,:,iupd)                       * updraft_rescale
        thl_updsum  = thl_updsum + a_upd(:,:,:,iupd) * thl_upd(:,:,:,iupd) * updraft_rescale 
         qt_updsum  =  qt_updsum + a_upd(:,:,:,iupd) *  qt_upd(:,:,:,iupd) * updraft_rescale 
         ql_updsum  =  ql_updsum + a_upd(:,:,:,iupd) *  ql_upd(:,:,:,iupd) * updraft_rescale 
          w_updsum  =   w_updsum + a_upd(:,:,:,iupd) *   w_upd(:,:,:,iupd)  
        ! Note: rescale the updraft fraction, but not rescaling total mass flux
        ! i.e., au -> au * upd_res wu -> wu / upd_res
    end do
    
end subroutine area_sum_updr


subroutine gaussian_mean(lower_lim, upper_lim, res_fac)

    real   , intent(in)  :: lower_lim, upper_lim
    real   , intent(out) :: res_fac
    !! NOT WRITTEN !!
    
    real :: upper_x, lower_x, upper_int, lower_int, pi_val
    
    pi_val = 4.E0 * atan(1.E0)
    
    call erfinv_calc(upper_lim, upper_x)
    call erfinv_calc(lower_lim, lower_x)
    
    upper_int = -exp(-upper_x*upper_x/2.0)/sqrt(2.0*pi_val)
    lower_int = -exp(-lower_x*lower_x/2.0)/sqrt(2.0*pi_val)
   
    res_fac = (upper_int - lower_int)/(upper_lim - lower_lim)

end subroutine gaussian_mean


! erfinv_calc: calculation of the inverse error function
!   Based on Wichura, Michael J. "Algorithm AS 241: The Percentage Points 
!   of the Normal Distribution." Journal of the Royal Statistical Society. 
!   Series C (Applied Statistics) 37, no. 3 (1988): 477-84.
subroutine erfinv_calc(P, PPND16)

    real   , intent(in)  :: P
    real   , intent(out) :: PPND16
    
    real :: Q, R
    
    ! real, parameter :: ZERO = 0.0E0
    real, parameter :: SMALL = 1.0E-15
    real, parameter :: ONE = 1.0E0
    real, parameter :: HALF = ONE/2.0E0
    real, parameter :: SPLIT1 = 0.425E0 
    real, parameter :: SPLIT2 = 5.0E0
    real, parameter :: CONST1 = 0.180625E0
    real, parameter :: CONST2 = 1.6E0

    real, parameter :: A0 = 3.3871328727963666080E0
    real, parameter :: A1 = 1.3314166789178437745E2
    real, parameter :: A2 = 1.9715909503065514427E3
    real, parameter :: A3 = 1.3731693765509461125E4
    real, parameter :: A4 = 4.5921953931549871457E4
    real, parameter :: A5 = 6.7265770927008700853E4
    real, parameter :: A6 = 3.3430575583588128105E4
    real, parameter :: A7 = 2.5090809287301226727E3

    real, parameter :: B1 = 4.2313330701600911252E1
    real, parameter :: B2 = 6.8718700749205790830E2
    real, parameter :: B3 = 5.3941960214247511077E3
    real, parameter :: B4 = 2.1213794301586595867E4
    real, parameter :: B5 = 3.9307895800092710610E4
    real, parameter :: B6 = 2.8729085735721942674E4
    real, parameter :: B7 = 5.2264952788528545610E3

    real, parameter :: C0 = 1.42343711074968357734E0
    real, parameter :: C1 = 4.63033784615654529590E0
    real, parameter :: C2 = 5.76949722146069140550E0
    real, parameter :: C3 = 3.64784832476320460504E0
    real, parameter :: C4 = 1.27045825245236838258E0
    real, parameter :: C5 = 2.41780725177450611770E-1
    real, parameter :: C6 = 2.27238449892691845833E-2
    real, parameter :: C7 = 7.74545014278341407640E-4

    real, parameter :: D1 = 2.05319162663775882187E0
    real, parameter :: D2 = 1.67638483018380384940E0
    real, parameter :: D3 = 6.89767334985100004550E-1
    real, parameter :: D4 = 1.48103976427480074590E-1
    real, parameter :: D5 = 1.51986665636164571966E-2
    real, parameter :: D6 = 5.47593808499534494600E-4
    real, parameter :: D7 = 1.05075007164441684324E-9

    real, parameter :: E0 = 6.65790464350110377720E0
    real, parameter :: E1 = 5.46378491116411436990E0
    real, parameter :: E2 = 1.78482653991729133580E0
    real, parameter :: E3 = 2.96560571828504891230E-1
    real, parameter :: E4 = 2.65321895265761230930E-2
    real, parameter :: E5 = 1.24266094738807843860E-3
    real, parameter :: E6 = 2.71155556874348757815E-5
    real, parameter :: E7 = 2.01033439929228813265E-7

    real, parameter :: F1 = 5.99832206555887937690E-1
    real, parameter :: F2 = 1.36929880922735805310E-1
    real, parameter :: F3 = 1.48753612908506148525E-2
    real, parameter :: F4 = 7.86869131145613259100E-4
    real, parameter :: F5 = 1.84631831751005468180E-5
    real, parameter :: F6 = 1.42151175831644588870E-7
    real, parameter :: F7 = 2.04426310338993978564E-15

    ! Calculation

    Q = P - HALF
    if (abs(Q) .le. SPLIT1) then
        R = CONST1 - Q * Q
        PPND16 = Q * (((((((A7 * R + A6) * R + A5) * R + A4) * R + A3)   &
          * R + A2) * R + A1) * R + A0) / (((((((B7 * R + B6) * R + B5)  &
          * R + B4) * R + B3) * R + B2) * R + B1) * R + ONE)
    else
        if (Q .lt. 0) then
            R = P
        else
            R = ONE - P
        endif
        
        if (R .le. SMALL) then
            ! PPND16 = ZERO
            R = SMALL
        endif
        
        R = sqrt(-log(R))
        if (R <= SPLIT2) then
             R = R - CONST2
             PPND16 = (((((((C7 * R + C6) * R + C5) * R + C4) * R     &
                 + C3) * R + C2) * R + C1) * R + C0) / (((((((D7 * R  &
                 + D6) * R + D5) * R + D4) * R + D3) * R + D2) * R    &
                 + D1) * R + ONE)
        else
             R = R - SPLIT2
             PPND16 = (((((((E7 * R + E6) * R + E5) * R + E4) * R     &
                 + E3) * R + E2) * R + E1) * R + E0) / (((((((F7 * R  &
                 + F6) * R + F5) * R + F4) * R + F3) * R + F2) * R    &
                 + F1) * R + ONE)
        endif
         
        if (Q .lt. 0) then
             PPND16 = -PPND16
        endif
         
    endif

end subroutine erfinv_calc

end module newedmf_mod
