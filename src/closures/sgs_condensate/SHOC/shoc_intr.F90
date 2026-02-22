module shoc_intr

  !------------------------------------------------------------------- !
  ! Module to interface E3SM with Simplified Higher                    !
  ! Order Closure (SHOC), developed                                    !
  !    by Peter Bogenschutz (Bogenschutz and Krueger 2013).            !
  !                                                                    !
  ! SHOC replaces the existing turbulence, shallow convection, and     !
  !   macrophysics in E3SM                                             !
  !                                                                    !
  !                                                                    !
  !---------------------------Code history---------------------------- !
  ! Authors:  P. Bogenschutz                                           !
  !                                                                    !
  !------------------------------------------------------------------- !

  use shr_kind_mod,  only: r8=>shr_kind_r8
  use shr_const_mod, only: SHR_CONST_PI
  use ppgrid,        only: pver, pverp
  use phys_control,  only: phys_getopts
  use physconst,     only: rair, cpair, gravit, latvap, latice, zvir, &
                           rh2o, karman, tms_orocnst, tms_z0fac
  use constituents,  only: pcnst, cnst_add
  use pbl_utils,     only: calc_ustar, calc_obklen
  use perf_mod,      only: t_startf, t_stopf
  use spmd_utils,    only: masterproc
  use cam_logfile,   only: iulog
  use shoc,          only: linear_interp, largeneg
  use spmd_utils,    only: masterproc
  use edmf,          only: do_edmf

  implicit none

  public :: shoc_init_cnst, shoc_implements_cnst

  ! define physics buffer indicies here
  integer :: tke_idx, &     ! turbulent kinetic energy
             tkh_idx, &
             tk_idx, &
             wthv_idx, &       ! buoyancy flux
             cld_idx, &          ! Cloud fraction
             concld_idx, &       ! Convective cloud fraction
             ast_idx, &          ! Stratiform cloud fraction
             alst_idx, &         ! Liquid stratiform cloud fraction
             aist_idx, &         ! Ice stratiform cloud fraction
             qlst_idx, &         ! Physical in-cloud LWC
             qist_idx, &         ! Physical in-cloud IWC
             dp_frac_idx, &      ! deep convection cloud fraction
	     cmeliq_idx, &       ! cmeliq_idx index in physics buffer
	     icwmrdp_idx, &      ! In cloud mixing ratio for deep convection
             sh_frac_idx, &      ! shallow convection cloud fraction
             relvar_idx, &       ! relative cloud water variance
	     kvh_idx, &          ! SHOC eddy diffusivity on thermo levels
             kvm_idx, &          ! SHOC eddy diffusivity on mom levels
             pblh_idx, &         ! PBL pbuf
             accre_enhan_idx, &  ! optional accretion enhancement factor for MG
             naai_idx, &         ! ice number concentration
             prer_evap_idx, &    ! rain evaporation rate
             qrl_idx, &          ! longwave cooling rate
             radf_idx, &
             tpert_idx, &
             fice_idx, &
	     vmag_gust_idx

  integer, public :: &
    ixtke = 0

  integer :: cmfmc_sh_idx = 0

  real(r8), parameter :: tke_tol = 0.0004_r8

  real(r8), parameter :: &
      theta0   = 300._r8, &             ! Reference temperature                     [K]
      ts_nudge = 86400._r8, &           ! Time scale for u/v nudging (not used)     [s]
      p0_shoc = 100000._r8, &
      shoc_tk1 = 268.15_r8, &
      shoc_tk2 = 238.15_r8, &
      shoc_liq_deep = 8.e-6, &
      shoc_liq_sh = 10.e-6, &
      shoc_ice_deep = 25.e-6, &
      shoc_ice_sh = 50.e-6

  logical      :: lq(pcnst)
  logical      :: lq2(pcnst)

  logical            :: history_budget
  integer            :: history_budget_histfile_num
  logical            :: micro_do_icesupersat

  character(len=16)  :: eddy_scheme      ! Default set in phys_control.F90
  character(len=16)  :: deep_scheme      ! Default set in phys_control.F90

  real(r8), parameter :: unset_r8 = huge(1.0_r8)

  real(r8) :: shoc_timestep = unset_r8  ! Default SHOC timestep set in namelist
  real(r8) :: dp1
  
  real(r8) :: shoc_length_fac = unset_r8
  
  integer :: edsclr_dim

  logical      :: prog_modal_aero
  real(r8) :: micro_mg_accre_enhan_fac = huge(1.0_r8) !Accretion enhancement factor from namelist

  integer, parameter :: ncnst=1
  character(len=8) :: cnst_names(ncnst)
  logical :: do_cnst=.true.

  logical :: liqcf_fix = .FALSE.  ! HW for liquid cloud fraction fix
  logical :: relvar_fix = .FALSE. !PMA for relvar fix

  contains

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

  subroutine shoc_register_e3sm()

#ifdef SHOC_SGS
    ! Add SHOC fields to pbuf
    use physics_buffer,  only: pbuf_add_field, dtype_r8, dyn_time_lvls
    use ppgrid,          only: pver, pverp, pcols

    call phys_getopts( eddy_scheme_out                 = eddy_scheme, &
                       deep_scheme_out                 = deep_scheme, &
                       history_budget_out              = history_budget, &
                       history_budget_histfile_num_out = history_budget_histfile_num, &
                       micro_do_icesupersat_out        = micro_do_icesupersat, &
                       micro_mg_accre_enhan_fac_out    = micro_mg_accre_enhan_fac)

    cnst_names=(/'TKE   '/)

    ! TKE is prognostic in SHOC and should be advected by dynamics
    call cnst_add('SHOC_TKE',0._r8,0._r8,0._r8,ixtke,longname='turbulent kinetic energy',cam_outfld=.false.)

    ! Fields that are not prognostic should be added to PBUF
    call pbuf_add_field('WTHV', 'global', dtype_r8, (/pcols,pver,dyn_time_lvls/), wthv_idx)
    call pbuf_add_field('TKH', 'global', dtype_r8, (/pcols,pver,dyn_time_lvls/), tkh_idx)
    call pbuf_add_field('TK', 'global', dtype_r8, (/pcols,pver,dyn_time_lvls/), tk_idx)

    call pbuf_add_field('pblh',       'global', dtype_r8, (/pcols/), pblh_idx)
    call pbuf_add_field('tke',        'global', dtype_r8, (/pcols, pverp/), tke_idx)
    call pbuf_add_field('kvh',        'global', dtype_r8, (/pcols, pverp/), kvh_idx)
    call pbuf_add_field('kvm',        'global', dtype_r8, (/pcols, pverp/), kvm_idx)
    call pbuf_add_field('tpert',      'global', dtype_r8, (/pcols/), tpert_idx)
    call pbuf_add_field('AST',        'global', dtype_r8, (/pcols,pver,dyn_time_lvls/),    ast_idx)
    call pbuf_add_field('AIST',       'global', dtype_r8, (/pcols,pver,dyn_time_lvls/),    aist_idx)
    call pbuf_add_field('ALST',       'global', dtype_r8, (/pcols,pver,dyn_time_lvls/),    alst_idx)
    call pbuf_add_field('QIST',       'global', dtype_r8, (/pcols,pver,dyn_time_lvls/),    qist_idx)
    call pbuf_add_field('QLST',       'global', dtype_r8, (/pcols,pver,dyn_time_lvls/),    qlst_idx)
    call pbuf_add_field('CONCLD',     'global', dtype_r8, (/pcols,pver,dyn_time_lvls/),    concld_idx)
    call pbuf_add_field('CLD',        'global', dtype_r8, (/pcols,pver,dyn_time_lvls/),    cld_idx)
    call pbuf_add_field('FICE',       'physpkg',dtype_r8, (/pcols,pver/), fice_idx)
    call pbuf_add_field('RAD_CLUBB',  'global', dtype_r8, (/pcols,pver/), radf_idx)
    call pbuf_add_field('CMELIQ',     'physpkg',dtype_r8, (/pcols,pver/), cmeliq_idx)

    call pbuf_add_field('vmag_gust',  'global', dtype_r8, (/pcols/),      vmag_gust_idx)

#endif

  end subroutine shoc_register_e3sm
  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

function shoc_implements_cnst(name)

  !--------------------------------------------------------------------
  ! Return true if specified constituent is implemented by this package
  !--------------------------------------------------------------------

  character(len=*), intent(in) :: name
  logical :: shoc_implements_cnst

  shoc_implements_cnst = (do_cnst .and. any(name == cnst_names))

end function shoc_implements_cnst

  subroutine shoc_init_cnst(name, q, gcid)

  !------------------------------------------------------------------- !
  ! Initialize the state for SHOC's prognostic variable                !
  !------------------------------------------------------------------- !

    character(len=*), intent(in)  :: name     ! constituent name
    real(r8),         intent(out) :: q(:,:)   ! mass mixing ratio (gcol, plev)
    integer,          intent(in)  :: gcid(:)  ! global column id

#ifdef SHOC_SGS
    if (trim(name) == trim('SHOC_TKE')) q = tke_tol
#endif

  end subroutine shoc_init_cnst

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

  subroutine shoc_readnl(nlfile)

  !------------------------------------------------------------------- !
  ! Read in any namelist parameters here                               !
  !   (currently none)                                                 !
  !------------------------------------------------------------------- !

    use units,           only: getunit, freeunit
    use namelist_utils,  only: find_group_name
    use cam_abortutils,  only: endrun
    use mpishorthand
    use edmf,            only: mf_readnl

    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

    integer :: iunit, read_status

    namelist /shocpbl_diff_nl/ shoc_timestep, shoc_length_fac
    
    !  Call SHOC+MF namelist
    call mf_readnl(nlfile)
    
    !  Read namelist to determine if SHOC history should be called
    if (masterproc) then
      iunit = getunit()
      open( iunit, file=trim(nlfile), status='old' )

      call find_group_name(iunit, 'shocpbl_diff_nl', status=read_status)
      if (read_status == 0) then
         read(unit=iunit, nml=shocpbl_diff_nl, iostat=read_status)
         if (read_status /= 0) then
            call endrun('shoc_readnl:  error reading namelist')
         end if
      end if

      close(unit=iunit)
      call freeunit(iunit)
    end if

#ifdef SPMD
! Broadcast namelist variables
      call mpibcast(shoc_timestep,           1,   mpir8,   0, mpicom)
      call mpibcast(shoc_length_fac,         1,   mpir8,   0, mpicom)
#endif

  end subroutine shoc_readnl

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

  subroutine shoc_init_e3sm(pbuf2d, dp1_in)

  !------------------------------------------------------------------- !
  ! Initialize SHOC for E3SM                                           !
  !------------------------------------------------------------------- !

    use physics_types,          only: physics_state, physics_ptend
    use ppgrid,                 only: pver, pverp, pcols
    use ref_pres,               only: pref_mid
    use time_manager,           only: is_first_step
    use hb_diff,                only: init_hb_diff
    use physics_buffer,         only: pbuf_get_index, pbuf_set_field, &
                                      physics_buffer_desc
    use rad_constituents,       only: rad_cnst_get_info, rad_cnst_get_mode_num_idx, &
                                      rad_cnst_get_mam_mmr_idx
    use constituents,           only: cnst_get_ind
    use shoc,                   only: shoc_init
    use cam_history,            only: horiz_only, addfld, add_default
    use error_messages,         only: handle_errmsg
    use trb_mtn_stress,         only: init_tms

    implicit none
    !  Input Variables
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    real(r8) :: dp1_in

    integer :: lptr
    integer :: nmodes, nspec, m, l
    integer :: ixnumliq
    integer :: ntop_shoc
    integer :: nbot_shoc
    character(len=128) :: errstring

    logical :: history_amwg

    lq(1:pcnst) = .true.
    edsclr_dim = pcnst

    !----- Begin Code -----

    ! ----------------------------------------------------------------- !
    ! Determine how many constituents SHOC will transport.  Note that
    ! SHOC does not transport aerosol consituents.  Therefore, need to
    ! determine how many aerosols constituents there are and subtract that
    ! off of pcnst (the total consituents)
    ! ----------------------------------------------------------------- !

    call phys_getopts(prog_modal_aero_out=prog_modal_aero, &
                      history_amwg_out = history_amwg, &
                      liqcf_fix_out   = liqcf_fix)

    ! Define physics buffers indexes
    cld_idx     = pbuf_get_index('CLD')         ! Cloud fraction
    concld_idx  = pbuf_get_index('CONCLD')      ! Convective cloud cover
    ast_idx     = pbuf_get_index('AST')         ! Stratiform cloud fraction
    alst_idx    = pbuf_get_index('ALST')        ! Liquid stratiform cloud fraction
    aist_idx    = pbuf_get_index('AIST')        ! Ice stratiform cloud fraction
    qlst_idx    = pbuf_get_index('QLST')        ! Physical in-stratus LWC
    qist_idx    = pbuf_get_index('QIST')        ! Physical in-stratus IWC
    dp_frac_idx = pbuf_get_index('DP_FRAC')     ! Deep convection cloud fraction
    icwmrdp_idx = pbuf_get_index('ICWMRDP')     ! In-cloud deep convective mixing ratio
    sh_frac_idx = pbuf_get_index('SH_FRAC')     ! Shallow convection cloud fraction
    relvar_idx      = pbuf_get_index('RELVAR')      ! Relative cloud water variance
    accre_enhan_idx = pbuf_get_index('ACCRE_ENHAN') ! accretion enhancement for MG
    prer_evap_idx   = pbuf_get_index('PRER_EVAP')
    qrl_idx         = pbuf_get_index('QRL')
    cmfmc_sh_idx    = pbuf_get_index('CMFMC_SH')
    tke_idx         = pbuf_get_index('tke')
    vmag_gust_idx = pbuf_get_index('vmag_gust')

    if (is_first_step()) then
      call pbuf_set_field(pbuf2d, wthv_idx, 0.0_r8)
      call pbuf_set_field(pbuf2d, tkh_idx, 0.0_r8)
      call pbuf_set_field(pbuf2d, tk_idx, 0.0_r8)
      call pbuf_set_field(pbuf2d, fice_idx, 0.0_r8)
      call pbuf_set_field(pbuf2d, tke_idx, tke_tol)
      call pbuf_set_field(pbuf2d, alst_idx, 0.0_r8)
      call pbuf_set_field(pbuf2d, aist_idx, 0.0_r8)

      call pbuf_set_field(pbuf2d, vmag_gust_idx,    1.0_r8)

    endif

    if (prog_modal_aero) then
       ! Turn off modal aerosols and decrement edsclr_dim accordingly
       call rad_cnst_get_info(0, nmodes=nmodes)

       do m = 1, nmodes
          call rad_cnst_get_mode_num_idx(m, lptr)
          lq(lptr)=.false.
          edsclr_dim = edsclr_dim-1

          call rad_cnst_get_info(0, m, nspec=nspec)
          do l = 1, nspec
             call rad_cnst_get_mam_mmr_idx(m, l, lptr)
             lq(lptr)=.false.
             edsclr_dim = edsclr_dim-1
          end do
       end do

       !  In addition, if running with MAM, droplet number is transported
       !  in dropmixnuc, therefore we do NOT want SHOC to apply transport
       !  tendencies to avoid double counted.  Else, we apply tendencies.
       call cnst_get_ind('NUMLIQ',ixnumliq)
       lq(ixnumliq) = .false.
       edsclr_dim = edsclr_dim-1
    endif

    ! Add SHOC fields
    call addfld('SHOC_TKE', (/'lev'/), 'A', 'm2/s2', 'TKE')
    call addfld('WTHV_SEC', (/'lev'/), 'A', 'K m/s', 'Buoyancy Flux')
    call addfld('SHOC_MIX', (/'lev'/), 'A', 'm', 'SHOC length scale')
    call addfld('TK',(/'lev'/), 'A', 'm2/s','Eddy viscosity for momentum')
    call addfld('TKH', (/'lev'/), 'A', 'm2/s', 'Eddy viscosity for heat')
    call addfld('W_SEC', (/'lev'/), 'A', 'm2/s2', 'Vertical velocity variance')
    call addfld('THL_SEC',(/'ilev'/), 'A', 'K2', 'Temperature variance')
    call addfld('QW_SEC',(/'ilev'/), 'A', 'kg2/kg2', 'Moisture variance')
    call addfld('QWTHL_SEC',(/'ilev'/), 'A', 'K kg/kg', 'Temperature and moisture correlation')
    call addfld('WTHL_SEC',(/'ilev'/), 'A', 'W/m2', 'Heat flux')
    call addfld('WQW_SEC',(/'ilev'/), 'A', 'W/m2', 'Moisture flux')
    call addfld('WTKE_SEC',(/'ilev'/), 'A', 'm3/s3', 'Vertical flux of turbulence')
    call addfld('UW_SEC',(/'ilev'/), 'A', 'm2/s2', 'Momentum flux')
    call addfld('VW_SEC',(/'ilev'/), 'A', 'm2/s2', 'Momentum flux')
    call addfld('W3',(/'ilev'/), 'A', 'm3/s3', 'Third moment vertical velocity')
    call addfld('WQL_SEC',(/'lev'/),'A', 'W/m2', 'Liquid water flux')
    call addfld('ISOTROPY',(/'lev'/),'A', 's', 'timescale')
    call addfld('CONCLD',(/'lev'/),  'A',        'fraction', 'Convective cloud cover')
    call addfld('BRUNT',(/'lev'/), 'A', 's-1', 'Brunt frequency')
    call addfld('RELVAR',(/'lev'/), 'A', 'kg/kg', 'SHOC cloud liquid relative variance')
    call addfld('TKE_DISSIP',(/'lev'/),'A', 'm2/s3', 'TKE budget: dissipation term')
    call addfld('TKE_PROD_B',(/'lev'/),'A', 'm2/s3', 'TKE budget: B production term')
    call addfld('TKE_PROD_SH',(/'lev'/),'A', 'm2/s3', 'TKE budget: shear production term')


    call add_default('SHOC_TKE', 1, ' ')
    call add_default('WTHV_SEC', 1, ' ')
    call add_default('SHOC_MIX', 1, ' ')
    call add_default('TK', 1, ' ')
    call add_default('TKH', 1, ' ')
    call add_default('W_SEC', 1, ' ')
    call add_default('THL_SEC', 1, ' ')
    call add_default('QW_SEC', 1, ' ')
    call add_default('QWTHL_SEC', 1, ' ')
    call add_default('WTHL_SEC', 1, ' ')
    call add_default('WQW_SEC', 1, ' ')
    call add_default('WTKE_SEC', 1, ' ')
    call add_default('UW_SEC', 1, ' ')
    call add_default('VW_SEC', 1, ' ')
    call add_default('W3', 1, ' ')
    call add_default('WQL_SEC', 1, ' ')
    call add_default('ISOTROPY',1,' ')
    call add_default('CONCLD',1,' ')
    call add_default('BRUNT',1,' ')
    call add_default('RELVAR',1,' ')
    call add_default('TKE_DISSIP',1,' ')
    call add_default('TKE_PROD_B',1,' ')
    call add_default('TKE_PROD_SH',1,' ')

    ! Extra output variables:
    call addfld('TEMPERATURE',(/'lev'/), 'A', 'K', 'Temperature')
    call addfld('THETA',(/'lev'/), 'A', 'K', 'Potential temperature')
    call addfld('THETAL',(/'lev'/), 'A', 'K', 'Liq. water potential temperature')
    call addfld('THETAV',(/'lev'/), 'A', 'K', 'Virtual potential temperature')
    call addfld('QW', (/'lev'/), 'A', 'kg/kg', 'Total water mixing ratio')
    call addfld('QV', (/'lev'/), 'A', 'kg/kg', 'Water vapour mixing ratio')
    call addfld('QL', (/'lev'/), 'A', 'kg/kg', 'Lquid water mixing ratio')
    call addfld('ZM',(/'ilev'/), 'A', 'm', 'Momentum heights')
    call addfld('ZT',(/'lev'/), 'A', 'm', 'Thermodynamic heights')
    call addfld('PBL_H', horiz_only, 'A', 'm', 'Planetary boundary layer depth')
    call addfld('rho', (/ 'lev' /), 'A', 'kg/m^3' , 'Density full levels')
    call addfld('rho_i', (/ 'ilev' /), 'A', 'kg/m^3' , 'Density mid levels')
    call addfld('CLD_SHOC', (/ 'lev' /), 'A', 'fraction' , 'cld from SHOC')
    call addfld('AST', (/ 'lev' /), 'A', 'fraction' , 'Stratiform cloud fraction from shoc_intr.F90')
    call addfld('DEEPCU', (/ 'lev' /), 'A', 'fraction' , 'Deep convection cloud fraction from shoc_intr.F90')
    call addfld('RI_B', (/ 'lev' /), 'A', 'm', 'Bulk Richardson number')
    call addfld('PBL_WTHL', horiz_only, 'A', 'm', 'JPL PBL as the first level where wthl = wthl_surf*0.01')
    call addfld('PBL_TKE', horiz_only, 'A', 'm', 'JPL PBL as the first level where tke = tke_surf*0.01')

    call add_default('TEMPERATURE',1,' ')
    call add_default('THETA',1,' ')
    call add_default('THETAL',1,' ')
    call add_default('THETAV',1,' ')
    call add_default('QW', 1, ' ')
    call add_default('QV', 1, ' ')
    call add_default('QL', 1, ' ')
    call add_default('ZM',1,' ')
    call add_default('ZT',1,' ')
    call add_default('PBL_H', 1, ' ')
    call add_default('rho', 1, ' ')
    call add_default('rho_i', 1, ' ')
    call add_default('CLD_SHOC', 1, ' ')
    call add_default('AST', 1, ' ')
    call add_default('DEEPCU', 1, ' ')
    call add_default('RI_B', 1, ' ')
    call add_default('PBL_WTHL', 1, ' ')
    call add_default('PBL_TKE', 1, ' ')

    if (do_edmf) then
      ! Add EDMF fields
      call addfld ( 'mf_dry_a'    , (/ 'ilev' /), 'A', 'fraction'    , 'Dry updraft area fraction (EDMF)' )
      call addfld ( 'mf_moist_a'  , (/ 'ilev' /), 'A', 'fraction'    , 'Moist updraft area fraction (EDMF)' )
      call addfld ( 'mf_moist_a_zt' , (/ 'lev' /), 'A', 'fraction'   , 'Moist updraft area fraction (EDMF) interpolated to the full levels' )
      call addfld ( 'mf_dry_w'    , (/ 'ilev' /), 'A', 'm/s'         , 'Dry updraft vertical velocity (EDMF)' )
      call addfld ( 'mf_moist_w'  , (/ 'ilev' /), 'A', 'm/s'         , 'Moist updraft vertical velocity (EDMF)' )
      call addfld ( 'mf_dry_qt'   , (/ 'ilev' /), 'A', 'kg/kg'       , 'Dry updraft total water mixing ratio (EDMF)' )
      call addfld ( 'mf_moist_qt' , (/ 'ilev' /), 'A', 'kg/kg'       , 'Moist updraft total water mixing ratio (EDMF)' )
      call addfld ( 'mf_dry_thl'  , (/ 'ilev' /), 'A', 'K'           , 'Dry updraft liquid-ice potential temperature (EDMF)' )
      call addfld ( 'mf_moist_thl', (/ 'ilev' /), 'A', 'K'           , 'Moist updraft liquid-ice potential temperature (EDMF)' )
      call addfld ( 'mf_dry_u'    , (/ 'ilev' /), 'A', 'm/s'         , 'Dry updraft zonal velocity (EDMF)' )
      call addfld ( 'mf_moist_u'  , (/ 'ilev' /), 'A', 'm/s'         , 'Moist updraft zonal velocity (EDMF)' )
      call addfld ( 'mf_dry_v'    , (/ 'ilev' /), 'A', 'm/s'         , 'Dry updraft meridional velocity (EDMF)' )
      call addfld ( 'mf_moist_v'  , (/ 'ilev' /), 'A', 'm/s'         , 'Moist updraft meridional velocity (EDMF)' )
      call addfld ( 'mf_moist_qc' , (/ 'ilev' /), 'A', 'kg/kg'       , 'Moist updraft condensate mixing ratio (EDMF)' )
      call addfld ( 'mf_moist_qc_zt' , (/ 'lev' /), 'A', 'kg/kg'     , 'Moist updraft condensate mixing ratio (EDMF) interpolated to the full levels' )
      call addfld ( 'mf_ae'       , (/ 'ilev' /), 'A', 'fraction'    , 'Environmental area fraction (EDMF)' )
      call addfld ( 'mf_aw'       , (/ 'ilev' /), 'A', 'm/s'         , 'Sum of a_i*w_i (EDMF)' )
      call addfld ( 'mf_awthv'    , (/ 'ilev' /), 'A', 'K m/s'       , 'Sum of a_i*w_i*thv_i (EDMF)' )
      call addfld ( 'mf_awthl'    , (/ 'ilev' /), 'A', 'K m/s'       , 'Sum of a_i*w_i*thl_i (EDMF)' )
      call addfld ( 'mf_awqt'     , (/ 'ilev' /), 'A', 'kg/kg m/s'   , 'Sum of a_i*w_i*q_ti (EDMF)' )
      call addfld ( 'mf_awql'     , (/ 'ilev' /), 'A', 'kg/kg m/s'   , 'Sum of a_i*w_i*q_li (EDMF)' )
      call addfld ( 'mf_awu'      , (/ 'ilev' /), 'A', 'm2/s2'       , 'Sum of a_i*w_i*u_i (EDMF)' )
      call addfld ( 'mf_awv'      , (/ 'ilev' /), 'A', 'm2/s2'       , 'Sum of a_i*w_i*v_i (EDMF)' )
      call addfld ( 'mf_thlflx'   , (/ 'ilev' /), 'A', 'K m/s'       , 'Kinematic thl flux by edmf (EDMF)' )
      call addfld ( 'mf_qtflx'    , (/ 'ilev' /), 'A', 'kg/kg m/s'   , 'Kinematic qt flux by edmf (EDMF)' )
      call addfld ( 'mf_thlflx_watts'    , (/ 'ilev' /), 'A', 'W/m2' , 'Dynamic thl flux by edmf (EDMF)' )
      call addfld ( 'mf_qtflx_watts'     , (/ 'ilev' /), 'A', 'W/m2' , 'Dynamic qt flux by edmf (EDMF)' )
      call addfld ( 'mf_dry_freq'  , horiz_only , 'A', 'fraction'    , 'Frequency of dry plume activation' )
      call addfld ( 'mf_moist_freq', horiz_only , 'A', 'fraction'    , 'Frequency of moist plume activation' )
      call addfld ( 'plumeheight', horiz_only , 'A', 'm'    , 'Plume height according to the one plume test' )
      call addfld ( 'plume_dry_height', horiz_only , 'A', 'm', 'Plume height according to the one dry plume test' )
      call addfld ( 'ztop_B', horiz_only , 'A', 'm', 'Mixed layer height according to min(wthv_sec)' )
      call addfld ( 'cfl_mf'      , (/ 'ilev' /), 'A', 'm/s'       , 'Mass-flux CFL number' )

     call addfld ( 'wthl_ed'    , (/ 'ilev' /), 'A', 'W/m2' , '(EDMF)' )
     call addfld ( 'wthl_mf'    , (/ 'ilev' /), 'A', 'W/m2' , '(EDMF)' )
     call addfld ( 'wqw_ed'    , (/ 'ilev' /), 'A', 'W/m2' , '(EDMF)' )
     call addfld ( 'wqw_mf'    , (/ 'ilev' /), 'A', 'W/m2' , '(EDMF)' )

      call addfld ( 'mf_nup', horiz_only, 'A', '[-]', 'Number of EDMF plumes')
      call addfld ( 'mf_L0', horiz_only, 'A', 'm', 'Distance a plume needs to travel to entrain once')
      call addfld ( 'mf_ent0', horiz_only, 'A', '[-]', 'Fractional mass of air entrained in a single entrainment event')
     
     
      call add_default( 'mf_dry_a'           , 1, ' ')
      call add_default( 'mf_moist_a'         , 1, ' ')
      call add_default( 'mf_moist_a_zt'      , 1, ' ')
      call add_default( 'mf_dry_w'           , 1, ' ')
      call add_default( 'mf_moist_w'         , 1, ' ')
      call add_default( 'mf_dry_qt'          , 1, ' ')
      call add_default( 'mf_moist_qt'        , 1, ' ')
      call add_default( 'mf_dry_thl'         , 1, ' ')
      call add_default( 'mf_moist_thl'       , 1, ' ')
      call add_default( 'mf_dry_u'           , 1, ' ')
      call add_default( 'mf_moist_u'         , 1, ' ')
      call add_default( 'mf_dry_v'           , 1, ' ')
      call add_default( 'mf_moist_v'         , 1, ' ')
      call add_default( 'mf_moist_qc'        , 1, ' ')
      call add_default( 'mf_moist_qc_zt'     , 1, ' ')
      call add_default( 'mf_ae'              , 1, ' ')
      call add_default( 'mf_aw'              , 1, ' ')
      call add_default( 'mf_awthv'           , 1, ' ')
      call add_default( 'mf_awthl'           , 1, ' ')
      call add_default( 'mf_awqt'            , 1, ' ')
      call add_default( 'mf_awql'            , 1, ' ')
      call add_default( 'mf_awu'             , 1, ' ')
      call add_default( 'mf_awv'             , 1, ' ')
      call add_default( 'mf_thlflx'          , 1, ' ')
      call add_default( 'mf_qtflx'           , 1, ' ')
      call add_default( 'mf_thlflx_watts'    , 1, ' ')
      call add_default( 'mf_qtflx_watts'     , 1, ' ')
      call add_default( 'mf_dry_freq'        , 1, ' ')
      call add_default( 'mf_moist_freq'      , 1, ' ')
      call add_default( 'plumeheight'        , 1, ' ')
      call add_default( 'plume_dry_height'   , 1, ' ')
      call add_default( 'ztop_B'             , 1, ' ')
      call add_default( 'cfl_mf'             , 1, ' ')

      call add_default( 'wthl_ed'          , 1, ' ')
      call add_default( 'wthl_mf'          , 1, ' ')
      call add_default( 'wqw_ed'          , 1, ' ')
      call add_default( 'wqw_mf'          , 1, ' ')
      
      call add_default( 'mf_nup' , 1, ' ')
      call add_default( 'mf_L0'  , 1, ' ')
      call add_default( 'mf_ent0', 1, ' ')
      
    endif


    ! ---------------------------------------------------------------!
    ! Initialize SHOC                                                !
    ! ---------------------------------------------------------------!

    ntop_shoc = 1    ! if >1, must be <= nbot_molec
    nbot_shoc = pver ! currently always pver

    call shoc_init( &
          pver, gravit, rair, rh2o, cpair, &
          zvir, latvap, latice, karman, &
          pref_mid, nbot_shoc, ntop_shoc, &
          shoc_length_fac )

    ! --------------- !
    ! End             !
    ! Initialization  !
    ! --------------- !

    dp1 = dp1_in

  end subroutine shoc_init_e3sm

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

  subroutine shoc_tend_e3sm( &
                             state, ptend_all, pbuf, hdtime, &
			     cmfmc, cam_in, sgh30, &
			     macmic_it, cld_macmic_num_steps, &
			     dlf, det_s, det_ice, alst_o)

  !------------------------------------------------------------------- !
  ! Provide tendencies of shallow convection , turbulence, and         !
  !   macrophysics from SHOC to E3SM                                   !
  !------------------------------------------------------------------- !

    use physics_types,  only: physics_state, physics_ptend, &
                              physics_state_copy, physics_ptend_init, &
                              physics_ptend_sum

    use physics_update_mod, only: physics_update

    use physics_buffer, only: pbuf_get_index, pbuf_old_tim_idx, pbuf_get_field, &
                              pbuf_set_field, physics_buffer_desc

    use ppgrid,         only: pver, pverp, pcols
    use constituents,   only: cnst_get_ind
    use camsrfexch,     only: cam_in_t
    use ref_pres,       only: top_lev => trop_cloud_top_lev
    use time_manager,   only: is_first_step
    use cam_abortutils, only: endrun
    use wv_saturation,  only: qsat
    use micro_mg_cam,   only: micro_mg_version
    use cldfrc2m,                  only: aist_vector
    use hb_diff,                   only: pblintd
    use trb_mtn_stress,            only: compute_tms
    use shoc,           only: shoc_main
    use cam_history,    only: outfld
    use scamMod,        only: single_column, iop_mode
    use edmf,           only: mf_nup, mf_L0, mf_ent0
    
    implicit none

   ! --------------- !
   ! Input Auguments !
   ! --------------- !

   type(physics_state), intent(in)    :: state                    ! Physics state variables                 [vary]
   type(cam_in_t),      intent(in)    :: cam_in
   real(r8),            intent(in)    :: hdtime                   ! Host model timestep                     [s]
   real(r8),            intent(in)    :: dlf(pcols,pver)          ! Detraining cld H20 from deep convection [kg/ks/s]
   real(r8),            intent(in)    :: cmfmc(pcols,pverp)       ! convective mass flux--m sub c           [kg/m2/s]
   real(r8),            intent(in)    :: sgh30(pcols)             ! std deviation of orography              [m]
   integer,             intent(in)    :: cld_macmic_num_steps     ! number of mac-mic iterations
   integer,             intent(in)    :: macmic_it                ! number of mac-mic iterations
   ! ---------------------- !
   ! Input-Output Auguments !
   ! ---------------------- !

   type(physics_buffer_desc), pointer :: pbuf(:)

   ! ---------------------- !
   ! Output Auguments !
   ! ---------------------- !

   type(physics_ptend), intent(out)   :: ptend_all                 ! package tendencies

   ! These two variables are needed for energy check
   real(r8),            intent(out)   :: det_s(pcols)              ! Integral of detrained static energy from ice
   real(r8),            intent(out)   :: det_ice(pcols)            ! Integral of detrained ice for energy check

   real(r8), intent(out) :: alst_o(pcols,pver)  ! H. Wang: for old liquid status fraction

   ! --------------- !
   ! Local Variables !
   ! --------------- !

   integer :: shoctop(pcols)

#ifdef SHOC_SGS

   type(physics_state) :: state1                ! Local copy of state variable
   type(physics_ptend) :: ptend_loc             ! Local tendency from processes, added up to return as ptend_all

   integer :: i, j, k, t, ixind, nadv
   integer :: ixcldice, ixcldliq, ixnumliq, ixnumice, ixq
   integer :: itim_old
   integer :: ncol, lchnk                       ! # of columns, and chunk identifier
   integer :: err_code                          ! Diagnostic, for if some calculation goes amiss.
   integer :: begin_height, end_height
   integer :: icnt

   real(r8) :: dtime                            ! SHOC time step                              [s]
   real(r8) :: edsclr_in(pcols,pver,edsclr_dim)      ! Scalars to be diffused through SHOC         [units vary]
   real(r8) :: edsclr_out(pcols,pver,edsclr_dim)
   real(r8) :: rcm_in(pcols,pver)
   real(r8) :: cloudfrac_shoc(pcols,pver)
   real(r8) :: newfice(pcols,pver)              ! fraction of ice in cloud at CLUBB start       [-]
   real(r8) :: exner(pcols,pver)
   real(r8) :: thlm(pcols,pver)
   real(r8) :: um(pcols,pver)
   real(r8) :: vm(pcols,pver)
   real(r8) :: rvm(pcols,pver)
   real(r8) :: rtm(pcols,pver)
   real(r8) :: rcm(pcols,pver)
   real(r8) :: rcm2(pcols,pver)                 ! cloud liquid variance                         [kg/kg]
   real(r8) :: ksrftms(pcols)                   ! Turbulent mountain stress surface drag        [kg/s/m2]
   real(r8) :: tautmsx(pcols)                   ! U component of turbulent mountain stress      [N/m2]
   real(r8) :: tautmsy(pcols)                   ! V component of turbulent mountain stress      [N/m2]
   real(r8) :: wm_zt(pcols,pver)
   real(r8) :: zt_g(pcols,pver)
   real(r8) :: dz_g(pcols,pver)
   real(r8) :: zi_g(pcols,pverp)
   real(r8) :: cloud_frac(pcols,pver)          ! CLUBB cloud fraction                          [fraction]
   real(r8) :: dlf2(pcols,pver)
   real(r8) :: isotropy(pcols,pver)
   real(r8) :: host_dx, host_dy
   real(r8) :: host_temp(pcols,pver)
   real(r8) :: host_dx_in(pcols), host_dy_in(pcols)
   real(r8) :: shoc_mix_out(pcols,pver), tk_in(pcols,pver), tkh_in(pcols,pver)
   real(r8) :: isotropy_out(pcols,pver), tke_zt(pcols,pver)
   real(r8) :: a_diss_out(pcols,pver), a_prod_bu_out(pcols,pver), a_prod_sh_out(pcols,pver)
   real(r8) :: w_sec_out(pcols,pver), thl_sec_out(pcols,pverp)
   real(r8) :: qw_sec_out(pcols,pverp), qwthl_sec_out(pcols,pverp)
   real(r8) :: wthl_sec_out(pcols,pverp), wqw_sec_out(pcols,pverp)
   real(r8) :: wtke_sec_out(pcols,pverp), uw_sec_out(pcols,pverp)
   real(r8) :: vw_sec_out(pcols,pverp), w3_out(pcols,pverp)
   real(r8) :: wqls_out(pcols,pver), brunt_out(pcols,pver)
   real(r8) :: rino(pcols,pver)
   real(r8) :: temperature(pcols,pver)  
   real(r8) :: theta(pcols,pver) 
   
   real(r8) :: wthl_output(pcols,pverp)
   real(r8) :: wqw_output(pcols,pverp)
   real(r8) :: wthv_output(pcols,pver)
   real(r8) :: wql_output(pcols,pver)
   
   ! MF, ED
   real(r8) :: wthl_sec_ed_out(pcols,pverp),wthl_sec_mf_out(pcols,pverp)
   real(r8) :: wthl_ed_out(pcols,pverp),wthl_mf_out(pcols,pverp)
   real(r8) :: wqw_sec_ed_out(pcols,pverp),wqw_sec_mf_out(pcols,pverp)
   real(r8) :: wqw_ed_out(pcols,pverp),wqw_mf_out(pcols,pverp)

   real(r8) :: obklen(pcols), ustar2(pcols), kinheat(pcols), kinwat(pcols)
   real(r8) :: dummy2(pcols), dummy3(pcols), kbfs(pcols), th(pcols,pver), thv(pcols,pver)
   real(r8) :: thv2(pcols,pver)

   real(r8) :: minqn, rrho(pcols,pver), rrho_i(pcols,pverp)    ! minimum total cloud liquid + ice threshold    [kg/kg]
   real(r8) :: cldthresh, frac_limit
   real(r8) :: ic_limit, dum1

   real(r8) :: wpthlp_sfc(pcols), wprtp_sfc(pcols), upwp_sfc(pcols), vpwp_sfc(pcols)
   real(r8) :: wtracer_sfc(pcols,edsclr_dim)

   ! For EDMF output
   real(r8), dimension(pcols,pverp) :: mf_dry_a, mf_moist_a, &
             mf_dry_w, mf_moist_w, &
             mf_dry_qt, mf_moist_qt, &
             mf_dry_thl, mf_moist_thl, &
             mf_dry_u, mf_moist_u, &
             mf_dry_v, mf_moist_v, &
             mf_moist_qc, cfl_mf, &
             mf_ae, mf_aw, mf_awthv, mf_awthl, mf_awqt, mf_awql, mf_awqi, mf_awu, mf_awv, &
             mf_thlflx, mf_qtflx, &         ! Kinematic fluxes
             mf_thlflx_out, mf_qtflx_out    ! Dynamic fluxes (W/m2)
   real(r8), dimension(pcols,pver) :: mf_ae_zt, mf_moist_a_zt, mf_moist_qc_zt
   real(r8), dimension(pcols)      :: mf_dry_freq, mf_moist_freq 
   real(r8), dimension(pcols)      :: plumeheight, plume_dry_height, ztop_B
   real(r8), dimension(pcols) :: pblh_wthl  ! PBL height for EDMF: first level where wthl = wthl_surf*0.01 (Common method in SBLs) [m]
   real(r8), dimension(pcols) :: pblh_tke   ! PBL height for EDMF: level where TKE vanishes (Common method in SBLs) [m]
   real(r8), dimension(pcols):: edmf_L0, edmf_ent0, edmf_nup
    
   ! Variables below are needed to compute energy integrals for conservation
   real(r8) :: ke_a(pcols), ke_b(pcols), te_a(pcols), te_b(pcols)
   real(r8) :: wv_a(pcols), wv_b(pcols), wl_b(pcols), wl_a(pcols)
   real(r8) :: se_dis(pcols), se_a(pcols), se_b(pcols), shoc_s(pcols,pver)
   real(r8) :: shoc_t(pcols,pver)

   ! --------------- !
   ! Pointers        !
   ! --------------- !

   real(r8), pointer, dimension(:,:) :: tke_zi  ! turbulent kinetic energy, interface
   real(r8), pointer, dimension(:,:) :: wthv    ! buoyancy flux
   real(r8), pointer, dimension(:,:) :: tkh
   real(r8), pointer, dimension(:,:) :: tk
   real(r8), pointer, dimension(:,:) :: cld      ! cloud fraction                               [fraction]
   real(r8), pointer, dimension(:,:) :: concld   ! convective cloud fraction                    [fraction]
   real(r8), pointer, dimension(:,:) :: ast      ! stratiform cloud fraction                    [fraction]
   real(r8), pointer, dimension(:,:) :: alst     ! liquid stratiform cloud fraction             [fraction]
   real(r8), pointer, dimension(:,:) :: aist     ! ice stratiform cloud fraction                [fraction]
   real(r8), pointer, dimension(:,:) :: cmeliq

   real(r8), pointer, dimension(:,:) :: qlst     ! Physical in-stratus LWC                      [kg/kg]
   real(r8), pointer, dimension(:,:) :: qist     ! Physical in-stratus IWC                      [kg/kg]
   real(r8), pointer, dimension(:,:) :: deepcu   ! deep convection cloud fraction               [fraction]
   real(r8), pointer, dimension(:,:) :: shalcu   ! shallow convection cloud fraction            [fraction]
   real(r8), pointer, dimension(:,:) :: khzt     ! eddy diffusivity on thermo levels            [m^2/s]
   real(r8), pointer, dimension(:,:) :: khzm     ! eddy diffusivity on momentum levels          [m^2/s]
   real(r8), pointer, dimension(:) :: pblh       ! planetary boundary layer height                [m]
   real(r8), pointer, dimension(:,:) :: dp_icwmr ! deep convection in cloud mixing ratio        [kg/kg]
   real(r8), pointer, dimension(:,:) :: cmfmc_sh ! Shallow convective mass flux--m subc (pcols,pverp) [kg/m2/s/]

   real(r8), pointer, dimension(:,:) :: prer_evap
   real(r8), pointer, dimension(:,:) :: accre_enhan
   real(r8), pointer, dimension(:,:) :: relvar

   logical :: lqice(pcnst)
   real(r8) :: relvarmax

   !------------------------------------------------------------------!
   !------------------------------------------------------------------!
   !------------------------------------------------------------------!
   !       MAIN COMPUTATION BEGINS HERE                               !
   !------------------------------------------------------------------!
   !------------------------------------------------------------------!
   !------------------------------------------------------------------!

 !  Get indicees for cloud and ice mass and cloud and ice number
   ic_limit   = 1.e-12_r8
   frac_limit = 0.01_r8

   call cnst_get_ind('Q',ixq)
   call cnst_get_ind('CLDLIQ',ixcldliq)
   call cnst_get_ind('CLDICE',ixcldice)
   call cnst_get_ind('NUMLIQ',ixnumliq)
   call cnst_get_ind('NUMICE',ixnumice)

   call physics_ptend_init(ptend_loc,state%psetcols, 'shoc', ls=.true., lu=.true., lv=.true., lq=lq)

   call physics_state_copy(state,state1)

   !  Determine number of columns and which chunk computation is to be performed on
   ncol = state%ncol
   lchnk = state%lchnk

   !  Determine time step of physics buffer
   itim_old = pbuf_old_tim_idx()

   !  Establish associations between pointers and physics buffer fields
   call pbuf_get_field(pbuf, tke_idx,     tke_zi)
   call pbuf_get_field(pbuf, wthv_idx,     wthv,     start=(/1,1,itim_old/), kount=(/pcols,pver,1/))
   call pbuf_get_field(pbuf, tkh_idx,      tkh,     start=(/1,1,itim_old/), kount=(/pcols,pver,1/))
   call pbuf_get_field(pbuf, tk_idx,       tk,     start=(/1,1,itim_old/), kount=(/pcols,pver,1/))
   call pbuf_get_field(pbuf, cld_idx,     cld,     start=(/1,1,itim_old/), kount=(/pcols,pver,1/))
   call pbuf_get_field(pbuf, concld_idx,  concld,  start=(/1,1,itim_old/), kount=(/pcols,pver,1/))
   call pbuf_get_field(pbuf, ast_idx,     ast,     start=(/1,1,itim_old/), kount=(/pcols,pver,1/))
   call pbuf_get_field(pbuf, alst_idx,    alst,    start=(/1,1,itim_old/), kount=(/pcols,pver,1/))
   call pbuf_get_field(pbuf, aist_idx,    aist,    start=(/1,1,itim_old/), kount=(/pcols,pver,1/))
   call pbuf_get_field(pbuf, qlst_idx,    qlst,    start=(/1,1,itim_old/), kount=(/pcols,pver,1/))
   call pbuf_get_field(pbuf, qist_idx,    qist,    start=(/1,1,itim_old/), kount=(/pcols,pver,1/))

   call pbuf_get_field(pbuf, prer_evap_idx, prer_evap)
   call pbuf_get_field(pbuf, accre_enhan_idx, accre_enhan)
   call pbuf_get_field(pbuf, cmeliq_idx,  cmeliq)
   call pbuf_get_field(pbuf, relvar_idx,  relvar)
   call pbuf_get_field(pbuf, dp_frac_idx, deepcu)
   call pbuf_get_field(pbuf, sh_frac_idx, shalcu)
   call pbuf_get_field(pbuf, kvm_idx,     khzm)
   call pbuf_get_field(pbuf, kvh_idx,     khzt)
   call pbuf_get_field(pbuf, pblh_idx,    pblh)
   call pbuf_get_field(pbuf, icwmrdp_idx, dp_icwmr)
   call pbuf_get_field(pbuf, cmfmc_sh_idx, cmfmc_sh)

   !  Determine SHOC time step and make it sub-step friendly
   !  For now we want SHOC time step to be 5 min.  However, there are certain
   !  instances when a 5 min time step will not be possible (based on
   !  host model time step or on macro-micro sub-stepping

   dtime = shoc_timestep

   !  Now check to see if dtime is greater than the host model
   !    (or sub stepped) time step.  If it is, then simply
   !    set it equal to the host (or sub step) time step.
   !    This section is mostly to deal with small host model
   !    time steps (or small sub-steps)

   if (dtime .gt. hdtime) then
     dtime = hdtime
   endif

   !  Now check to see if SHOC time step divides evenly into
   !    the host model time step.  If not, force it to divide evenly.
   !    We also want it to be 5 minutes or less.  This section is
   !    mainly for host model time steps that are not evenly divisible
   !    by 5 minutes

   if (mod(hdtime,dtime) .ne. 0) then
     dtime = hdtime/2._r8
     do while (dtime .gt. 300._r8)
       dtime = dtime/2._r8
     end do
   endif

   !  If resulting host model time step and SHOC time step do not divide evenly
   !    into each other, have model throw a fit.

   if (mod(hdtime,dtime) .ne. 0) then
     call endrun('shoc_tend_e3sm:  SHOC time step and HOST time step NOT compatible')
   endif

   !  determine number of timesteps SHOC core should be advanced,
   !  host time step divided by SHOC time step
   nadv = max(hdtime/dtime,1._r8)

   ! Set grid space, in meters. If SCM, set to a grid size representative
   !  of a typical GCM.  Otherwise, compute locally.
   print*,'IOP_MODE = ',iop_mode
   if (single_column .and. .not. iop_mode) then
     host_dx_in(:) = 100000._r8
     host_dy_in(:) = 100000._r8
   else if (iop_mode) then
     call grid_size_uniform(host_dx, host_dy)
     host_dx_in(:) = host_dx
     host_dy_in(:) = host_dy
   else
     call grid_size(state1, host_dx_in, host_dy_in)
   endif

   minqn = 0._r8
   newfice(:,:) = 0._r8
   where(state1%q(:ncol,:pver,ixcldice) .gt. minqn) &
       newfice(:ncol,:pver) = state1%q(:ncol,:pver,ixcldice)/(state1%q(:ncol,:pver,ixcldliq)+state1%q(:ncol,:pver,ixcldice))

   do k=1,pver
     do i=1,ncol
       exner(i,k) = 1._r8/((state1%pmid(i,k)/p0_shoc)**(rair/cpair))
     enddo
   enddo

   !  At each SHOC call, initialize mean momentum  and thermo SHOC state
   !  from the E3SM state

   do k=1,pver   ! loop over levels
     do i=1,ncol ! loop over columns

       rvm(i,k) = state1%q(i,k,ixq)
       rcm(i,k) = state1%q(i,k,ixcldliq)
       rtm(i,k) = rvm(i,k) + rcm(i,k)
       
       !if (um(i,k) .EQ. 0._r8) then
       !   um(i,k) = state1%u(i,k)
       !   !print*,'um(i,k) = state1%u(i,k)'
       !endif   
       um(i,k) = state1%u(i,k)
       vm(i,k) = state1%v(i,k)

       temperature(i,k) = state1%t(i,k)
       theta(i,k) = state1%t(i,k)*exner(i,k)
       thlm(i,k) = state1%t(i,k)*exner(i,k)-(latvap/cpair)*state1%q(i,k,ixcldliq)
       thv(i,k) = state1%t(i,k)*exner(i,k)*(1.0_r8+zvir*state1%q(i,k,ixq)-state1%q(i,k,ixcldliq))

       tke_zt(i,k) = max(tke_tol,state1%q(i,k,ixtke))

       ! Cloud fraction needs to be initialized for first
       !  PBL height calculation call
       cloud_frac(i,k) = alst(i,k)

     enddo
   enddo

   ! ------------------------------------------------- !
   ! Prepare inputs for SHOC call                      !
   ! ------------------------------------------------- !

   do k=1,pver
     do i=1,ncol
       dz_g(i,k) = state1%zi(i,k)-state1%zi(i,k+1)  ! compute thickness
     enddo
   enddo

   !  Define the SHOC thermodynamic grid (in units of m)
   wm_zt(:,pver) = 0._r8
   do k=1,pver
     do i=1,ncol
       zt_g(i,k) = state1%zm(i,k)-state1%zi(i,pver+1)
       rrho(i,k)=(1._r8/gravit)*(state1%pdel(i,k)/dz_g(i,k))
       wm_zt(i,k) = -1._r8*state1%omega(i,k)/(rrho(i,k)*gravit)
       shoc_s(i,k) = state1%s(i,k)
     enddo
   enddo

   do k=1,pverp
     do i=1,ncol
       zi_g(i,k) = state1%zi(i,k)-state1%zi(i,pver+1)
     enddo
   enddo

   ! Get Density on the interface grid
   call linear_interp(zt_g(:ncol,:pver),zi_g(:ncol,:pverp),rrho(:ncol,:pver),&
                      rrho_i(:ncol,:pverp),pver,pverp,ncol,0._r8)

   do i=1,ncol
      !  Surface fluxes provided by host model

      wpthlp_sfc(i) = cam_in%shf(i)/(cpair*rrho_i(i,pverp))   ! Sensible heat flux
      wprtp_sfc(i)  = cam_in%cflx(i,1)/(rrho_i(i,pverp))      ! Latent heat flux
      upwp_sfc(i)   = cam_in%wsx(i)/rrho_i(i,pverp)           ! Surface meridional momentum flux
      vpwp_sfc(i)   = cam_in%wsy(i)/rrho_i(i,pverp)           ! Surface zonal momentum flux
      wtracer_sfc(i,:) = 0._r8  ! in E3SM tracer fluxes are done elsewhere
   enddo

   !  Do the same for tracers
   icnt=0
   do ixind=1,pcnst
     if (lq(ixind))  then
       icnt=icnt+1
       do k=1,pver
         do i=1,ncol
           edsclr_in(i,k,icnt) = state1%q(i,k,ixind)
         enddo
       enddo
     end if
   enddo

   ! ------------------------------------------------- !
   ! Actually call SHOC                                !
   ! ------------------------------------------------- !

   call shoc_main( &
        ncol, pver, pverp, dtime, nadv, &                                                         ! Input
	    host_dx_in(:ncol), host_dy_in(:ncol), thv(:ncol,:),&                                      ! Input
        zt_g(:ncol,:), zi_g(:ncol,:), state%pmid(:ncol,:pver), state%pint(:ncol,:pverp), &        ! Input
        state1%pdel(:ncol,:pver), &                                                               ! Input
	    wpthlp_sfc(:ncol), wprtp_sfc(:ncol), upwp_sfc(:ncol), vpwp_sfc(:ncol), &                  ! Input
	    wtracer_sfc(:ncol,:), edsclr_dim, wm_zt(:ncol,:), &                                       ! Input
	    exner(:ncol,:),state1%phis(:ncol), &                                                      ! Input
	    shoc_s(:ncol,:), tke_zt(:ncol,:), thlm(:ncol,:), rtm(:ncol,:), &                          ! Input/Ouput
	    um(:ncol,:), vm(:ncol,:), edsclr_in(:ncol,:,:), &                                         ! Input/Output
	    wthv(:ncol,:),tkh(:ncol,:),tk(:ncol,:), &                                                 ! Input/Output
	    rcm(:ncol,:),cloud_frac(:ncol,:), &                                                       ! Input/Output
        pblh(:ncol),pblh_wthl(:ncol),pblh_tke(:ncol),&                                            ! Output
        rino(:ncol,:), &                                                                          ! Output
        shoc_mix_out(:ncol,:), isotropy_out(:ncol,:), &                                           ! Output (diagnostic)
        a_diss_out(:ncol,:), a_prod_bu_out(:ncol,:), a_prod_sh_out(:ncol,:), &                    ! Output (diagnostic)
        w_sec_out(:ncol,:), thl_sec_out(:ncol,:), qw_sec_out(:ncol,:), qwthl_sec_out(:ncol,:), &  ! Output (diagnostic)
        wthl_sec_out(:ncol,:), wqw_sec_out(:ncol,:), wtke_sec_out(:ncol,:), &                     ! Output (diagnostic)
        uw_sec_out(:ncol,:), vw_sec_out(:ncol,:), w3_out(:ncol,:), &                              ! Output (diagnostic)
        wqls_out(:ncol,:),brunt_out(:ncol,:),rcm2(:ncol,:), &                                     ! Output (diagnostic)
        mf_dry_a(:ncol,:), mf_moist_a(:ncol,:), &                                                 ! Output (EDMF diagnostic)
        mf_dry_w(:ncol,:), mf_moist_w(:ncol,:), &                                                 ! Output (EDMF diagnostic)
        mf_dry_qt(:ncol,:), mf_moist_qt(:ncol,:), &                                               ! Output (EDMF diagnostic)
        mf_dry_thl(:ncol,:), mf_moist_thl(:ncol,:), &                                             ! Output (EDMF diagnostic)
        mf_dry_u(:ncol,:),   mf_moist_u(:ncol,:), &                                               ! Output (EDMF diagnostic)
        mf_dry_v(:ncol,:),   mf_moist_v(:ncol,:), &                                               ! Output (EDMF diagnostic)
                             mf_moist_qc(:ncol,:), &                                              ! Output (EDMF diagnostic)
        mf_thlflx(:ncol,:), mf_qtflx(:ncol,:), &                                                  ! Output (EDMF diagnostic)
        mf_ae(:ncol,:), mf_aw(:ncol,:), &                                                         ! Output (EDMF diagnostic)
        mf_awthv(:ncol,:), &                                                                      ! Output (EDMF diagnostic)
        mf_awthl(:ncol,:), mf_awqt(:ncol,:), &                                                    ! Output (EDMF diagnostic)
        mf_awql(:ncol,:), mf_awqi(:ncol,:), &                                                     ! Output (EDMF diagnostic)
        mf_awu(:ncol,:), mf_awv(:ncol,:), cfl_mf(:ncol,:), &                                      ! Output (EDMF diagnostic)
        mf_dry_freq(:), mf_moist_freq(:), plumeheight(:), plume_dry_height(:), ztop_B(:), &       ! Output (EDMF diagnostic)
        wthl_sec_ed_out(:ncol,:), wthl_sec_mf_out(:ncol,:), &                                     ! Output (EDMF diagnostic)
        wqw_sec_ed_out(:ncol,:), wqw_sec_mf_out(:ncol,:))                   					  ! Output (EDMF diagnostic)  
                    
   ! Transfer back to pbuf variables
   do k=1,pver
     do i=1,ncol
     
       if (do_edmf) then
         call linear_interp(zi_g(:ncol,:pverp),zt_g(:ncol,:pver),mf_ae(:ncol,:pverp),&
                        mf_ae_zt(:ncol,:pver),pverp,pver,ncol,0._r8)
         cloud_frac(i,k) = min(mf_ae_zt(i,k)*cloud_frac(i,k),1._r8) 
       else
         cloud_frac(i,k) = min(cloud_frac(i,k),1._r8) 
       endif

       do ixind=1,edsclr_dim
         edsclr_out(i,k,ixind) = edsclr_in(i,k,ixind)
       enddo

     enddo
   enddo

   ! Eddy diffusivities and TKE are needed for aerosol activation code.
   !   Linearly interpolate from midpoint grid and onto the interface grid.
   !   The output variables for these routines (khzm, khzt, and tke_zi)
   !   are PBUF pointers
   call linear_interp(state%zm(:ncol,:pver),state%zi(:ncol,:pverp),&
                tk(:ncol,:pver),khzm(:ncol,:pverp),pver,pverp,ncol,0._r8)
   call linear_interp(state%zm(:ncol,:pver),state%zi(:ncol,:pverp),&
                tkh(:ncol,:pver),khzt(:ncol,:pverp),pver,pverp,ncol,0._r8)
   call linear_interp(state%zm(:ncol,:pver),state%zi(:ncol,:pverp),&
                tke_zt(:ncol,:pver),tke_zi(:ncol,:pverp),pver,pverp,ncol,tke_tol)

   !  Now compute the tendencies of SHOC to E3SM
   do k=1,pver
     do i=1,ncol

       ptend_loc%u(i,k) = (um(i,k)-state1%u(i,k))/hdtime
       ptend_loc%v(i,k) = (vm(i,k)-state1%v(i,k))/hdtime
       ptend_loc%q(i,k,ixq) = (rtm(i,k)-rcm(i,k)-state1%q(i,k,ixq))/hdtime ! water vapor
       ptend_loc%q(i,k,ixcldliq) = (rcm(i,k)-state1%q(i,k,ixcldliq))/hdtime   ! Tendency of liquid water
       ptend_loc%s(i,k) = (shoc_s(i,k)-state1%s(i,k))/hdtime

       ptend_loc%q(i,k,ixtke)=(tke_zt(i,k)-state1%q(i,k,ixtke))/hdtime ! TKE

   !  Apply tendencies to ice mixing ratio, liquid and ice number, and aerosol constituents.
   !  Loading up this array doesn't mean the tendencies are applied.
   !  edsclr_out is compressed with just the constituents being used, ptend and state are not compressed

       icnt=0
       do ixind=1,pcnst
         if (lq(ixind)) then
           icnt=icnt+1
           if ((ixind /= ixq) .and. (ixind /= ixcldliq) .and. (ixind /= ixtke)) then
             ptend_loc%q(i,k,ixind) = (edsclr_out(i,k,icnt)-state1%q(i,k,ixind))/hdtime ! transported constituents
           end if
         end if
       enddo

     enddo
   enddo

   cmeliq(:,:) = ptend_loc%q(:,:,ixcldliq)

   ! Update physics tendencies
   call physics_ptend_init(ptend_all, state%psetcols, 'shoc')
   call physics_ptend_sum(ptend_loc,ptend_all,ncol)
   call physics_update(state1,ptend_loc,hdtime)

   ! ------------------------------------------------------------ !
   ! ------------------------------------------------------------ !
   ! ------------------------------------------------------------ !
   ! The rest of the code deals with diagnosing variables         !
   ! for microphysics/radiation computation and macrophysics      !
   ! ------------------------------------------------------------ !
   ! ------------------------------------------------------------ !
   ! ------------------------------------------------------------ !

   ! --------------------------------------------------------------------------------- !
   !  COMPUTE THE ICE CLOUD DETRAINMENT                                                !
   !  Detrainment of convective condensate into the environment or stratiform cloud    !
   ! --------------------------------------------------------------------------------- !

   !  Initialize the shallow convective detrainment rate, will always be zero
   dlf2(:,:) = 0.0_r8

   lqice(:)        = .false.
   lqice(ixcldliq) = .true.
   lqice(ixcldice) = .true.
   lqice(ixnumliq) = .true.
   lqice(ixnumice) = .true.

   call physics_ptend_init(ptend_loc,state%psetcols, 'clubb_det', ls=.true., lq=lqice)
   do k=1,pver
      do i=1,ncol
         if( state1%t(i,k) > shoc_tk1 ) then
            dum1 = 0.0_r8
         elseif ( state1%t(i,k) < shoc_tk2 ) then
            dum1 = 1.0_r8
         else
            !Note: Denominator is changed from 30.0_r8 to (shoc_tk1 - shoc_tk2),
            !(clubb_tk1 - clubb_tk2) is also 30.0 but it introduced a non-bfb change
            dum1 = ( shoc_tk1 - state1%t(i,k) ) /(shoc_tk1 - shoc_tk2)
         endif

         ptend_loc%q(i,k,ixcldliq) = dlf(i,k) * ( 1._r8 - dum1 )
         ptend_loc%q(i,k,ixcldice) = dlf(i,k) * dum1
         ptend_loc%q(i,k,ixnumliq) = 3._r8 * ( max(0._r8, ( dlf(i,k) - dlf2(i,k) )) * ( 1._r8 - dum1 ) ) &
                                     / (4._r8*3.14_r8* shoc_liq_deep**3*997._r8) + & ! Deep    Convection
                                     3._r8 * (                         dlf2(i,k)    * ( 1._r8 - dum1 ) ) &
                                     / (4._r8*3.14_r8*shoc_liq_sh**3*997._r8)     ! Shallow Convection
         ptend_loc%q(i,k,ixnumice) = 3._r8 * ( max(0._r8, ( dlf(i,k) - dlf2(i,k) )) *  dum1 ) &
                                     / (4._r8*3.14_r8*shoc_ice_deep**3*500._r8) + & ! Deep    Convection
                                     3._r8 * (                         dlf2(i,k)    *  dum1 ) &
                                     / (4._r8*3.14_r8*shoc_ice_sh**3*500._r8)     ! Shallow Convection
         ptend_loc%s(i,k)          = dlf(i,k) * dum1 * latice

         ! Only rliq is saved from deep convection, which is the reserved liquid.  We need to keep
         !   track of the integrals of ice and static energy that is effected from conversion to ice
         !   so that the energy checker doesn't complain.
         det_s(i)                  = det_s(i) + ptend_loc%s(i,k)*state1%pdel(i,k)/gravit
         det_ice(i)                = det_ice(i) - ptend_loc%q(i,k,ixcldice)*state1%pdel(i,k)/gravit

      enddo
    enddo

    det_ice(:ncol) = det_ice(:ncol)/1000._r8  ! divide by density of water

    call physics_ptend_sum(ptend_loc,ptend_all,ncol)
    call physics_update(state1,ptend_loc,hdtime)

    ! For purposes of this implementation, just set relvar and accre_enhan to 1
    relvar(:,:) = 1.0_r8
    accre_enhan(:,:) = 1._r8

! +++ JShpund: add relative cloud liquid variance (a vectorized version based on CLUBB)
!     TODO: double check the hardcoded values ('relvarmax', '0.001_r8')
    relvarmax = 10.0_r8
    where (rcm(:ncol,:pver) /= 0.0 .and. rcm2(:ncol,:pver) /= 0.0) &
           relvar(:ncol,:pver) = min(relvarmax,max(0.001_r8,rcm(:ncol,:pver)**2.0/rcm2(:ncol,:pver)))

    ! --------------------------------------------------------------------------------- !
    !  Diagnose some quantities that are computed in macrop_tend here.                  !
    !  These are inputs required for the microphysics calculation.                      !
    !                                                                                   !
    !  FIRST PART COMPUTES THE STRATIFORM CLOUD FRACTION FROM SHOC CLOUD FRACTION       !
    ! --------------------------------------------------------------------------------- !

    ! HW: set alst to alst_o before getting updated
    if(liqcf_fix) then
      if(.not.is_first_step()) alst_o(:ncol,:pver) = alst(:ncol,:pver)
    endif

    !  initialize variables
    alst(:,:) = 0.0_r8
    qlst(:,:) = 0.0_r8

    do k=1,pver
      do i=1,ncol
        alst(i,k) = cloud_frac(i,k)
        qlst(i,k) = rcm(i,k)/max(0.01_r8,alst(i,k))  ! Incloud stratus condensate mixing ratio
      enddo
    enddo

    ! HW
    if(liqcf_fix) then
      if(is_first_step()) alst_o(:ncol,:pver) = alst(:ncol,:pver)
    endif

    ! --------------------------------------------------------------------------------- !
    !  THIS PART COMPUTES CONVECTIVE AND DEEP CONVECTIVE CLOUD FRACTION                 !
    ! --------------------------------------------------------------------------------- !
    deepcu(:,pver) = 0.0_r8
    shalcu(:,pver) = 0.0_r8
    
    do k=1,pver-1
      do i=1,ncol
        !  diagnose the deep convective cloud fraction, as done in macrophysics based on the
        !  deep convective mass flux, read in from pbuf.  Since shallow convection is never
        !  called, the shallow convective mass flux will ALWAYS be zero, ensuring that this cloud
        !  fraction is purely from deep convection scheme.
        deepcu(i,k) = max(0.0_r8,min(dp1*log(1.0_r8+500.0_r8*(cmfmc(i,k+1)-cmfmc_sh(i,k+1))),0.6_r8))
        shalcu(i,k) = 0._r8
        
        if (do_edmf) then
          call linear_interp(zi_g(:ncol,:pverp),zt_g(:ncol,:pver),mf_moist_a(:ncol,:pverp),&
                        mf_moist_a_zt(:ncol,:pver),pverp,pver,ncol,0._r8)
          call linear_interp(zi_g(:ncol,:pverp),zt_g(:ncol,:pver),mf_moist_qc(:ncol,:pverp),&
                        mf_moist_qc_zt(:ncol,:pver),pverp,pver,ncol,0._r8)              
        !  shalcu(i,k) = mf_moist_a_zt(i,k)              
        !else
        !  shalcu(i,k) = 0._r8
        endif

        if (deepcu(i,k) <= frac_limit .or. dp_icwmr(i,k) < ic_limit) then
          deepcu(i,k) = 0._r8
        endif

        !  using the deep convective cloud fraction, and SHOC cloud fraction (variable
        !  "cloud_frac"), compute the convective cloud fraction.  This follows the formulation
        !  found in macrophysics code.  Assumes that convective cloud is all nonstratiform cloud
        !  from SHOC plus the deep convective cloud fraction
        concld(i,k) = min(cloud_frac(i,k)-alst(i,k)+shalcu(i,k)+deepcu(i,k),0.80_r8)
      enddo
    enddo

    ! --------------------------------------------------------------------------------- !
    !  COMPUTE THE ICE CLOUD FRACTION PORTION                                           !
    !  use the aist_vector function to compute the ice cloud fraction                   !
    ! --------------------------------------------------------------------------------- !

    do k=1,pver
      call aist_vector(state1%q(:,k,ixq),state1%t(:,k),state1%pmid(:,k),state1%q(:,k,ixcldice), &
           state1%q(:,k,ixnumice),cam_in%landfrac(:),cam_in%snowhland(:),aist(:,k),ncol)
    enddo

    ! --------------------------------------------------------------------------------- !
    !  THIS PART COMPUTES THE LIQUID STRATUS FRACTION                                   !
    !                                                                                   !
    !  For now leave the computation of ice stratus fraction from macrop_driver intact  !
    !  because SHOC does nothing with ice.  Here I simply overwrite the liquid stratus !
    !  fraction that was coded in macrop_driver                                         !
    ! --------------------------------------------------------------------------------- !

    !  Recompute net stratus fraction using maximum over-lapping assumption, as done
    !  in macrophysics code, using alst computed above and aist read in from physics buffer

    cldthresh=1.e-18_r8

    do k=1,pver
      do i=1,ncol

        ast(i,k) = max(alst(i,k),aist(i,k))

        qist(i,k) = state1%q(i,k,ixcldice)/max(0.01_r8,aist(i,k))
      enddo
    enddo

    !  Probably need to add deepcu cloud fraction to the cloud fraction array, else would just
    !  be outputting the shallow convective cloud fraction

    do k=1,pver
      do i=1,ncol
        cloud_frac(i,k) = min(ast(i,k)+shalcu(i,k)+deepcu(i,k),1.0_r8)
      enddo
    enddo

    cld(:,1:pver) = cloud_frac(:,1:pver)

    ! --------------------------------------------------------!
    ! Output fields
    !---------------------------------------------------------!
    do k=1,pverp
      do i=1,ncol
        wthl_output(i,k) = wthl_sec_out(i,k) * rrho_i(i,k) * cpair
        wqw_output(i,k) = wqw_sec_out(i,k) * rrho_i(i,k) * latvap
        
        if (do_edmf) then
          ! Convert EDMF kinematic fluxes into dynamic fluxes
          mf_thlflx_out(i,k) = mf_thlflx(i,k) * rrho_i(i,k) * cpair
          mf_qtflx_out(i,k)  = mf_qtflx(i,k) * rrho_i(i,k) * latvap
          
          wthl_ed_out(i,k) = wthl_sec_ed_out(i,k) * rrho_i(i,k) * cpair
          wthl_mf_out(i,k) = wthl_sec_mf_out(i,k) * rrho_i(i,k) * cpair
          
          wqw_ed_out(i,k) = wqw_sec_ed_out(i,k) * rrho_i(i,k) * latvap
          wqw_mf_out(i,k) = wqw_sec_mf_out(i,k) * rrho_i(i,k) * latvap
        end if
      enddo
    enddo

    do k=1,pver
      do i=1,ncol
        wql_output(i,k) = wqls_out(i,k) * rrho(i,k) * latvap
        wthv_output(i,k) = wthv(i,k) * rrho(i,k) * cpair
      enddo
    enddo
    
    call outfld('SHOC_TKE', tke_zt, pcols, lchnk)
    call outfld('WTHV_SEC', wthv_output, pcols, lchnk)
    call outfld('SHOC_MIX', shoc_mix_out, pcols, lchnk)
    call outfld('TK', tk, pcols, lchnk)
    call outfld('TKH', tkh, pcols, lchnk)
    call outfld('W_SEC', w_sec_out, pcols, lchnk)
    call outfld('THL_SEC', thl_sec_out, pcols, lchnk)
    call outfld('QW_SEC', qw_sec_out, pcols, lchnk)
    call outfld('QWTHL_SEC', qwthl_sec_out, pcols, lchnk)
    call outfld('WTHL_SEC', wthl_output, pcols, lchnk)
    call outfld('WQW_SEC', wqw_output, pcols, lchnk)
    call outfld('WTKE_SEC', wtke_sec_out, pcols, lchnk)
    call outfld('UW_SEC', uw_sec_out, pcols, lchnk)
    call outfld('VW_SEC', vw_sec_out, pcols, lchnk)
    call outfld('W3', w3_out, pcols, lchnk)
    call outfld('WQL_SEC',wql_output, pcols, lchnk)
    call outfld('ISOTROPY',isotropy_out, pcols,lchnk)
    call outfld('CONCLD',concld,pcols,lchnk)
    call outfld('BRUNT',brunt_out,pcols,lchnk)
    call outfld('RELVAR',relvar,pcols,lchnk)
    call outfld('TKE_DISSIP',a_diss_out, pcols,lchnk)
    call outfld('TKE_PROD_B',a_prod_bu_out, pcols,lchnk)
    call outfld('TKE_PROD_SH',a_prod_sh_out, pcols,lchnk)

    !Extra output variables
    call outfld('TEMPERATURE',temperature,pcols,lchnk)
    call outfld('THETA',theta,pcols,lchnk)
    call outfld('THETAL',thlm,pcols,lchnk)
    call outfld('THETAV',thv,pcols,lchnk)
    call outfld('QW', rtm, pcols, lchnk)
    call outfld('QV', rvm, pcols, lchnk)
    call outfld('QL', rcm, pcols, lchnk)
    call outfld('ZM',zi_g,pcols,lchnk)
    call outfld('ZT',zt_g,pcols,lchnk)
    call outfld('PBL_H', pblh, pcols, lchnk)
    call outfld('rho_i', rrho_i, pcols, lchnk)
    call outfld('rho', rrho, pcols, lchnk)
    call outfld('CLD_SHOC', cld, pcols, lchnk)
    call outfld('AST', ast, pcols, lchnk)
    call outfld('DEEPCU', deepcu, pcols, lchnk)
    call outfld('RI_B', rino, pcols, lchnk)
    call outfld('PBL_WTHL', pblh_wthl, pcols, lchnk)
    call outfld('PBL_TKE', pblh_tke, pcols, lchnk)
   
    if (do_edmf) then
      ! EDMF outputs
      call outfld( 'mf_dry_a'    , mf_dry_a,     pcols, lchnk )
      call outfld( 'mf_moist_a'  , mf_moist_a,   pcols, lchnk )
      call outfld( 'mf_moist_a_zt', mf_moist_a_zt, pcols, lchnk )
      call outfld( 'mf_dry_w'    , mf_dry_w,     pcols, lchnk )
      call outfld( 'mf_moist_w'  , mf_moist_w,   pcols, lchnk )
      call outfld( 'mf_dry_qt'   , mf_dry_qt,    pcols, lchnk )
      call outfld( 'mf_moist_qt' , mf_moist_qt,  pcols, lchnk )
      call outfld( 'mf_dry_thl'  , mf_dry_thl,   pcols, lchnk )
      call outfld( 'mf_moist_thl', mf_moist_thl, pcols, lchnk )
      call outfld( 'mf_dry_u'    , mf_dry_u,     pcols, lchnk )
      call outfld( 'mf_moist_u'  , mf_moist_u,   pcols, lchnk )
      call outfld( 'mf_dry_v'    , mf_dry_v,     pcols, lchnk )
      call outfld( 'mf_moist_v'  , mf_moist_v,   pcols, lchnk )
      call outfld( 'mf_moist_qc' , mf_moist_qc,  pcols, lchnk )
      call outfld( 'mf_moist_qc_zt', mf_moist_qc_zt,  pcols, lchnk )
      call outfld( 'mf_ae'       , mf_ae,        pcols, lchnk )
      call outfld( 'mf_aw'       , mf_aw,        pcols, lchnk )
      call outfld( 'mf_awthv'    , mf_awthv,     pcols, lchnk )
      call outfld( 'mf_awthl'    , mf_awthl,     pcols, lchnk )
      call outfld( 'mf_awqt'     , mf_awqt,      pcols, lchnk )
      call outfld( 'mf_awql'     , mf_awql,      pcols, lchnk )
      call outfld( 'mf_awu'      , mf_awu,       pcols, lchnk )
      call outfld( 'mf_awv'      , mf_awv,       pcols, lchnk )
      call outfld( 'mf_thlflx'   , mf_thlflx,    pcols, lchnk )
      call outfld( 'mf_qtflx'    , mf_qtflx,     pcols, lchnk )
      call outfld( 'mf_qtflx_watts'    , mf_qtflx_out,     pcols, lchnk )
      call outfld( 'mf_thlflx_watts'   , mf_thlflx_out,    pcols, lchnk )
      call outfld( 'plumeheight', plumeheight, pcols, lchnk)
      call outfld( 'plume_dry_height', plume_dry_height, pcols, lchnk)
      call outfld( 'ztop_B', ztop_B, pcols, lchnk)
      call outfld( 'cfl_mf'      , cfl_mf,       pcols, lchnk )

      call outfld( 'wthl_ed', wthl_ed_out, pcols, lchnk )
      call outfld( 'wthl_mf', wthl_mf_out, pcols, lchnk )
      call outfld( 'wqw_ed', wqw_ed_out, pcols, lchnk )
      call outfld( 'wqw_mf', wqw_mf_out, pcols, lchnk )

      edmf_nup = real(mf_nup,r8)
      edmf_L0 = mf_L0
      edmf_ent0 = mf_ent0
      call outfld( 'mf_nup', edmf_nup, pcols, lchnk)
      call outfld( 'mf_L0', edmf_L0, pcols, lchnk)
      call outfld( 'mf_ent0', edmf_ent0, pcols, lchnk)
     
    end if

#endif
    return
  end subroutine shoc_tend_e3sm

  subroutine grid_size(state, grid_dx, grid_dy)
  ! Determine the size of the grid for each of the columns in state

  use phys_grid,       only: get_area_p
  use shr_const_mod,   only: shr_const_pi
  use physics_types,   only: physics_state
  use ppgrid,          only: pver, pverp, pcols

  type(physics_state), intent(in) :: state
  real(r8), intent(out)           :: grid_dx(pcols), grid_dy(pcols)   ! E3SM grid [m]

  real(r8), parameter :: earth_ellipsoid1 = 111132.92_r8 ! World Geodetic System 1984 (WGS84)
                                                         ! first coefficient, meters per degree longitude at equator
  real(r8), parameter :: earth_ellipsoid2 = 559.82_r8 ! second expansion coefficient for WGS84 ellipsoid
  real(r8), parameter :: earth_ellipsoid3 = 1.175_r8 ! third expansion coefficient for WGS84 ellipsoid

  real(r8) :: mpdeglat, column_area, degree
  integer  :: i

  do i=1,state%ncol
      ! determine the column area in radians
      column_area = get_area_p(state%lchnk,i)
      ! convert to degrees
      degree = sqrt(column_area)*(180._r8/shr_const_pi)

      ! Now find meters per degree latitude
      ! Below equation finds distance between two points on an ellipsoid, derived from expansion
      !  taking into account ellipsoid using World Geodetic System (WGS84) reference
      mpdeglat = earth_ellipsoid1 - earth_ellipsoid2 * cos(2._r8*state%lat(i)) + earth_ellipsoid3 * cos(4._r8*state%lat(i))
      grid_dx(i) = mpdeglat * degree
      grid_dy(i) = grid_dx(i) ! Assume these are the same
  enddo

  end subroutine grid_size

  subroutine grid_size_uniform(grid_dx, grid_dy)

    ! Estimate grid box size at equator using
    !  the earth radius set for this case.  This assumes
    !  that all grid points are uniform, which is
    !  reasonable for IOP mode (not currently compatible with RRM)

    use physical_constants, only: rearth, dd_pi
    use dimensions_mod, only: np, ne

    real(r8), intent(out) :: grid_dx, grid_dy

    grid_dx = dd_pi*rearth/(2000.d0*dble(ne*(np-1)))
    grid_dx = grid_dx*1000._r8
    grid_dy = grid_dx

  end subroutine grid_size_uniform

end module shoc_intr
