
module cloud_diagnostics

!---------------------------------------------------------------------------------
! Purpose:
!
! Put cloud physical specifications on the history tape
!  Modified from code that computed cloud optics
!
! Author: Byron Boville  Sept 06, 2002
!  Modified Oct 15, 2008
!    
!
!---------------------------------------------------------------------------------

   use shr_kind_mod,  only: r8=>shr_kind_r8
   use ppgrid,        only: pcols, pver,pverp
   use physconst,     only: gravit
   use cam_history,   only: outfld
   use cam_history,   only: addfld, horiz_only, add_default

   ! YQIN 10/25/22
   use spmd_utils,          only: masterproc
   use cam_logfile,         only: iulog
   use cam_history_support, only: fillvalue
   use infnan,              only: isnan, isinf
   use cam_abortutils,      only: endrun
   use wv_saturation,       only: qsat
   use physconst,           only: cappa

   ! YQIN 01/14/22
   use cdiag_pdf,    only: pdf1d_regime, pdf2d_regime, pdf3d_regime, N_REGIME, regime, N_CAT, cats,&
   ncfs, &
   mcfs, &
   bcfs_1d, &
   bcfs, &
   nlnticlwps, & 
   mlnticlwps, & 
   blnticlwps_1d, &
   blnticlwps, &
   nlncdncs, &
   mlncdncs, &
   blncdncs_1d, &
   blncdncs, &
   nrels, &
   mrels, &
   brels_1d, &
   brels, &
   nlnccns, &
   mlnccns, &
   blnccns_1d, &
   blnccns, &
   ncods, &
   mcods, &
   bcods_1d, &
   bcods, &
   n3ds, &
   m3ds, &
   b3ds_1d, &
   b3ds, &
   minccn, mincdnc, minticlwp, mincf, minrel, mincod,&
   LTS_threshold, &
   EIS_threshold_dry, EIS_threshold_wet, &
   RH750_threshold


   implicit none
   private
   save

   public :: cloud_diagnostics_init
   public :: cloud_diagnostics_calc
   public :: cloud_diagnostics_register


! Local variables
   integer :: dei_idx, mu_idx, lambda_idx, iciwp_idx, iclwp_idx, cld_idx  ! index into pbuf for cloud fields
   integer :: ixcldice, ixcldliq, rei_idx, rel_idx


   ! YQIN 10/23/22
   integer :: ccn02_idx,icwnc_idx,cdnumc_avg_idx,cdnumc_925_idx,cdnumc_wavg_idx
   integer :: ticlwp_idx = -1 ! ticlwp
   integer :: cllow_idx = -1 ! cllow
   integer :: lts_idx ! LTS
   integer :: rh750_idx ! RH750
   integer :: eis_idx ! EIS 
!   integer :: t_cldtop_idx = -1 ! cloud-top temperature
   integer :: ticiwp_idx = -1 ! ticiwp

   logical :: do_cld_diag, mg_clouds, rk_clouds
   integer :: conv_water_in_rad
   
   integer :: cicewp_idx = -1
   integer :: cliqwp_idx = -1
   integer :: cldemis_idx = -1
   integer :: cldtau_idx = -1
   integer :: nmxrgn_idx = -1
   integer :: pmxrgn_idx = -1

   ! Index fields for precipitation efficiency.
   integer :: acpr_idx, acgcme_idx, acnum_idx



contains

!===============================================================================
  subroutine cloud_diagnostics_register

    use phys_control,  only: phys_getopts
    use physics_buffer,only: pbuf_add_field, dtype_r8, dtype_i4


    character(len=16) :: rad_pkg, microp_pgk

    call phys_getopts(radiation_scheme_out=rad_pkg,microp_scheme_out=microp_pgk)
    rk_clouds = microp_pgk == 'RK'
    mg_clouds = microp_pgk == 'MG'

    if (rk_clouds) then
       call pbuf_add_field('CLDEMIS','physpkg', dtype_r8,(/pcols,pver/), cldemis_idx)
       call pbuf_add_field('CLDTAU', 'physpkg', dtype_r8,(/pcols,pver/), cldtau_idx)

       call pbuf_add_field('CICEWP', 'physpkg', dtype_r8,(/pcols,pver/), cicewp_idx)
       call pbuf_add_field('CLIQWP', 'physpkg', dtype_r8,(/pcols,pver/), cliqwp_idx)

       call pbuf_add_field('PMXRGN', 'physpkg', dtype_r8,(/pcols,pverp/), pmxrgn_idx)
       call pbuf_add_field('NMXRGN', 'physpkg', dtype_i4,(/pcols /),      nmxrgn_idx)
    else if (mg_clouds) then
       ! In cloud ice water path for radiation
       call pbuf_add_field('ICIWP',      'global', dtype_r8,(/pcols,pver/), iciwp_idx)
       ! In cloud liquid water path for radiation
       call pbuf_add_field('ICLWP',      'global', dtype_r8,(/pcols,pver/), iclwp_idx)

       ! YQIN 11/10/22 add into the pbuf
       call pbuf_add_field('TICLWP',     'global', dtype_r8, (/pcols/),       ticlwp_idx)
       call pbuf_add_field('CLLOW',      'global', dtype_r8, (/pcols/),       cllow_idx)
       call pbuf_add_field('CLIQWP',     'global', dtype_r8, (/pcols,pver/),  cliqwp_idx)
       call pbuf_add_field('TICIWP',     'global', dtype_r8, (/pcols/),       ticiwp_idx)

    endif
  end subroutine cloud_diagnostics_register

!===============================================================================
  subroutine cloud_diagnostics_init()
!-----------------------------------------------------------------------
    use physics_buffer,only: pbuf_get_index
    use phys_control,  only: phys_getopts
    use constituents,  only: cnst_get_ind
    use cloud_cover_diags, only: cloud_cover_diags_init

    implicit none

!-----------------------------------------------------------------------

    character(len=16) :: wpunits, sampling_seq
    logical           :: history_amwg          ! output the variables used by the AMWG diag package
    logical           :: history_verbose       ! produce verbose history output

    !-----------------------------------------------------------------------

    cld_idx    = pbuf_get_index('CLD')

    ! ----------------------------
    ! determine default variables
    ! ----------------------------
    call phys_getopts( history_amwg_out = history_amwg, &
                       history_verbose_out = history_verbose )

    if (mg_clouds) then

       call addfld ('ICWMR', (/ 'lev' /), 'A', 'kg/kg', 'Prognostic in-cloud water mixing ratio'                  )
       call addfld ('ICIMR', (/ 'lev' /), 'A', 'kg/kg', 'Prognostic in-cloud ice mixing ratio'                    )
       call addfld ('IWC', (/ 'lev' /), 'A', 'kg/m3', 'Grid box average ice water content'                      )
       call addfld ('LWC', (/ 'lev' /), 'A', 'kg/m3', 'Grid box average liquid water content'                   )


       if (history_amwg) then
          call add_default ('ICWMR', 1, ' ')
          call add_default ('ICIMR', 1, ' ')
          call add_default ('IWC      ', 1, ' ')
       end if

       dei_idx    = pbuf_get_index('DEI')
       mu_idx     = pbuf_get_index('MU')
       lambda_idx = pbuf_get_index('LAMBDAC')

       ! YQIN 10/23/22
       ccn02_idx  = pbuf_get_index('CCN02')
       icwnc_idx = pbuf_get_index('ICWNC')
       cdnumc_avg_idx = pbuf_get_index('CDNUMC_AVG')
       cdnumc_925_idx = pbuf_get_index('CDNUMC_925')
       cdnumc_wavg_idx = pbuf_get_index('CDNUMC_WAVG')
       eis_idx = pbuf_get_index('EIS')
       lts_idx = pbuf_get_index('LTS')
       rh750_idx = pbuf_get_index('RH750')
       rel_idx = pbuf_get_index('REL')


    elseif (rk_clouds) then

       rei_idx    = pbuf_get_index('REI')
       rel_idx    = pbuf_get_index('REL')

    endif

    call cnst_get_ind('CLDICE', ixcldice)
    call cnst_get_ind('CLDLIQ', ixcldliq)

    do_cld_diag = rk_clouds .or. mg_clouds

    if (.not.do_cld_diag) return
    
    call phys_getopts(conv_water_in_rad_out=conv_water_in_rad)

    if (rk_clouds) then 
       wpunits = 'gram/m2'
       sampling_seq='rad_lwsw'
    else if (mg_clouds) then 
       wpunits = 'kg/m2'
       sampling_seq=''
    endif

    call addfld ('ICLDIWP', (/ 'lev' /), 'A', wpunits,'In-cloud ice water path'               , sampling_seq=sampling_seq)
    call addfld ('ICLDTWP', (/ 'lev' /), 'A',wpunits,'In-cloud cloud total water path (liquid and ice)', &
         sampling_seq=sampling_seq)

    call addfld ('GCLDLWP',(/ 'lev' /), 'A',wpunits,'Grid-box cloud water path'             , &
         sampling_seq=sampling_seq)
    call addfld ('TGCLDCWP',horiz_only,    'A',wpunits,'Total grid-box cloud water path (liquid and ice)', &
         sampling_seq=sampling_seq, standard_name='atmosphere_mass_content_of_cloud_condensed_water')

    call addfld ('TGCLDLWP',horiz_only,    'A',wpunits,'Total grid-box cloud liquid water path', &
         sampling_seq=sampling_seq, standard_name='atmosphere_mass_content_of_cloud_liquid_water')

    call addfld ('TGCLDIWP',horiz_only,    'A',wpunits,'Total grid-box cloud ice water path'   , &
         sampling_seq=sampling_seq, standard_name='atmosphere_mass_content_of_cloud_ice')
    
    if(mg_clouds) then
       call addfld ('lambda_cloud',(/ 'lev' /),'I','1/meter','lambda in cloud')
       call addfld ('mu_cloud',(/ 'lev' /),'I','1','mu in cloud')
       call addfld ('dei_cloud',(/ 'lev' /),'I','micrometers','ice radiative effective diameter in cloud')

       ! YQIN 10/23/22
       call addfld ('ICLWP', horiz_only, 'A', 'g/m2', 'vertically-integrated in-cloud liquid water path', sampling_seq=sampling_seq)
       call addfld ('ICIWP', horiz_only, 'A', 'g/m2', 'vertically-integrated in-cloud ice water path', sampling_seq=sampling_seq)
       call addfld ('TREL', horiz_only, 'A', 'micron', 'cloud layer mean REL', sampling_seq=sampling_seq)
!       call addfld ('T_CLDTOP', horiz_only, 'A', 'K', 'cloud top temperature', sampling_seq=sampling_seq, flag_xyfill=.true., fill_value = fillvalue)

       ! YQIN 12/19/22
       call add_default ('ICLWP', 1, ' ')
       call add_default ('ICIWP', 1, ' ')
       call add_default ('TREL', 1, ' ')
!       call add_default ('T_CLDTOP', 1, ' ')


    endif

    if(rk_clouds) then
       call addfld ('rel_cloud',(/ 'lev' /),'I','1/meter','effective radius of liq in cloud', sampling_seq=sampling_seq)
       call addfld ('rei_cloud',(/ 'lev' /),'I','1','effective radius of ice in cloud', sampling_seq=sampling_seq)
    endif

    call addfld ('SETLWP',(/ 'lev' /), 'A','gram/m2','Prescribed liquid water path'          , sampling_seq=sampling_seq)
    call addfld ('LWSH',horiz_only,    'A','m','Liquid water scale height'             , sampling_seq=sampling_seq)
    call addfld ('EFFCLD',(/ 'lev' /), 'A','fraction','Effective cloud fraction'              , sampling_seq=sampling_seq)
    call addfld ('EMISCLD', (/ 'lev' /), 'A', '1','cloud emissivity'                      , sampling_seq=sampling_seq)

    call cloud_cover_diags_init(sampling_seq)


    if (history_amwg) then
       call add_default ('TGCLDLWP', 1, ' ')
       call add_default ('TGCLDIWP', 1, ' ')
       call add_default ('TGCLDCWP', 1, ' ')
       if (history_verbose) call add_default ('EMISCLD', 1, ' ')
    endif

    return
  end subroutine cloud_diagnostics_init

subroutine cloud_diagnostics_calc(state,  pbuf)
!===============================================================================
!
! Compute (liquid+ice) water path and cloud water/ice diagnostics
! *** soon this code will compute liquid and ice paths from input liquid and ice mixing ratios
! 
! **** mixes interface and physics code temporarily
!-----------------------------------------------------------------------
    use physics_types, only: physics_state    
    use physics_buffer,only: physics_buffer_desc, pbuf_get_field, pbuf_old_tim_idx
    use pkg_cldoptics, only: cldovrlap, cldclw,  cldems
    use conv_water,    only: conv_water_4rad
    use radiation,     only: radiation_do
    use cloud_cover_diags, only: cloud_cover_diags_out

    use ref_pres,       only: top_lev=>trop_cloud_top_lev

    ! YQIN 10/25/22
    use micro_mg_utils, only: qsmall,mincld
    use interpolate_data,only: vertinterp

    implicit none

! Arguments
    type(physics_state), intent(in)    :: state        ! state variables
    type(physics_buffer_desc), pointer :: pbuf(:)

! Local variables

    real(r8), pointer :: cld(:,:)       ! cloud fraction
    real(r8), pointer :: iciwp(:,:)   ! in-cloud cloud ice water path
    real(r8), pointer :: iclwp(:,:)   ! in-cloud cloud liquid water path
    real(r8), pointer :: dei(:,:)       ! effective radiative diameter of ice
    real(r8), pointer :: mu(:,:)        ! gamma distribution for liq clouds
    real(r8), pointer :: lambda(:,:)    ! gamma distribution for liq clouds
    real(r8), pointer :: rei(:,:)       ! effective radiative radius of ice
    real(r8), pointer :: rel(:,:)       ! effective radiative radius of liq

    ! YQIN 10/23/22
    real(r8), pointer :: icwnc(:,:)     ! in-cloud CDNC
    real(r8), pointer :: ccn02(:,:)     ! CCN at S = 0.2%
    real(r8), pointer :: cdnumc_avg(:)  ! cloud layer mean CDNC unit: #/cm3
    real(r8), pointer :: cdnumc_925(:)  ! in-cloud CDNC at 925 hPa
    real(r8), pointer :: cdnumc_wavg(:) ! cloud layer pressure-weighted mean CDNC unit: #/cm3


    real(r8), pointer :: cldemis(:,:)   ! cloud emissivity
    real(r8), pointer :: cldtau(:,:)    ! cloud optical depth
    real(r8), pointer :: cicewp(:,:)    ! in-cloud cloud ice water path
    real(r8), pointer :: cliqwp(:,:)    ! in-cloud cloud liquid water path

    integer,  pointer :: nmxrgn(:)      ! Number of maximally overlapped regions
    real(r8), pointer :: pmxrgn(:,:)    ! Maximum values of pressure for each

    integer :: itim_old

    real(r8) :: cwp   (pcols,pver)      ! in-cloud cloud (total) water path
    real(r8) :: gicewp(pcols,pver)      ! grid-box cloud ice water path
    real(r8) :: gliqwp(pcols,pver)      ! grid-box cloud liquid water path
    real(r8) :: gwp   (pcols,pver)      ! grid-box cloud (total) water path
    real(r8) :: tgicewp(pcols)          ! Vertically integrated ice water path
    real(r8) :: tgliqwp(pcols)          ! Vertically integrated liquid water path
    real(r8) :: tgwp   (pcols)          ! Vertically integrated (total) cloud water path

    ! YQIN 10/23/22
    real(r8), pointer :: ticlwp(:)      ! vertically integrated in-cloud liquid water path
    real(r8), pointer :: ticiwp(:)      ! vertically integrated in-cloud ice water path 
    real(r8) :: trel(pcols)             ! vertically integrated in-cloud effective radius (REL)
    real(r8), pointer :: cllow(:)       ! low cloud fraction 
    real(r8) :: cltot(pcols)            ! total cloud fraction


    real(r8) :: ficemr (pcols,pver)     ! Ice fraction from ice and liquid mixing ratios

    real(r8) :: icimr(pcols,pver)       ! In cloud ice mixing ratio
    real(r8) :: icwmr(pcols,pver)       ! In cloud water mixing ratio
    real(r8) :: iwc(pcols,pver)         ! Grid box average ice water content
    real(r8) :: lwc(pcols,pver)         ! Grid box average liquid water content

! old data
    real(r8) :: tpw    (pcols)          ! total precipitable water
    real(r8) :: clwpold(pcols,pver)     ! Presribed cloud liq. h2o path
    real(r8) :: hl     (pcols)          ! Liquid water scale height

    integer :: i,k                      ! loop indexes
    integer :: ncol, lchnk
    real(r8) :: rgrav

    real(r8) :: allcld_ice (pcols,pver) ! Convective cloud ice
    real(r8) :: allcld_liq (pcols,pver) ! Convective cloud liquid

    real(r8) :: effcld(pcols,pver)      ! effective cloud=cld*emis

    logical :: dosw,dolw
  
    ! YQIN 10/23/22
    real(r8) :: pdf_cf           (pcols,N_REGIME,ncfs)
    real(r8) :: pdf_rel          (pcols,N_REGIME,nrels) 
    real(r8) :: pdf_lnccn        (pcols,N_REGIME,nlnccns)
    real(r8) :: pdf_lncdnc       (pcols,N_REGIME,nlncdncs)
    real(r8) :: pdf_lncdnc925    (pcols,N_REGIME,nlncdncs)
    real(r8) :: pdf_lncdncWAVG   (pcols,N_REGIME,nlncdncs)
    real(r8) :: pdf_lnticlwp     (pcols,N_REGIME,nlnticlwps)

    real(r8) :: pdf_lnccn_lncdnc      (pcols,N_REGIME,nlnccns,nlncdncs)
    real(r8) :: pdf_lncdnc_rel        (pcols,N_REGIME,nlncdncs,nrels)

    real(r8) :: pdf_NDLWPCF     (pcols,N_REGIME,nlncdncs,nlnticlwps,ncfs)
    real(r8) :: pdf_NDLWPCF_out (pcols,n3ds)

    real(r8) :: pdf_relp             (pcols,N_REGIME,nrels) 
    real(r8) :: pdf_cfp              (pcols,N_REGIME,ncfs)
    real(r8) :: pdf_lnccnp           (pcols,N_REGIME,nlnccns)
    real(r8) :: pdf_lncdncp          (pcols,N_REGIME,nlncdncs)
    real(r8) :: pdf_lnticlwpp        (pcols,N_REGIME,nlnticlwps)

    real(r8) :: pdf_lnccnp_lncdncp   (pcols,N_REGIME,nlnccns,nlncdncs) ! PDF of CCN and CDNC
    real(r8) :: pdf_lncdncp_relp     (pcols,N_REGIME,nlncdncs,nrels)

    real(r8) :: pdf_NDLWPCFp     (pcols,N_REGIME,nlncdncs,nlnticlwps,ncfs)
    real(r8) :: pdf_NDLWPCFp_out (pcols,n3ds)

    ! YQIN 12/15/22 -- add LTS and RH750 to define regime
    real(r8), pointer :: LTS(:) ! K
    real(r8), pointer :: EIS(:) ! K
    real(r8), pointer :: RH750(:) 

    real(r8) :: LTSin(pcols), LTSin_threshold_dry,LTSin_threshold_wet

    integer :: j,l,lkj
    integer :: ireg
    integer :: p,w

    real(r8) :: cld_h(pcols,pver), cllow_h(pcols)
    real(r8) :: rel_h(pcols,pver), trel_h(pcols)
    real(r8) :: ccn02_sfc_h(pcols), ccn02_h(pcols,pver)
    real(r8) :: cdnumc_avg_h(pcols), icwnc_h(pcols,pver)
    real(r8) :: cdnumc_925_h(pcols)
    real(r8) :: ticlwp_h(pcols), cliqwp_h(pcols,pver)
    integer :: itype 


    real(r8) :: p_surf_t1(pcols)    ! data interpolated to a pressure surface
    real(r8) :: p_surf_t2(pcols)    ! data interpolated to a pressure surface
    real(r8) :: tem2(pcols,pver)    ! temporary workspace
    real(r8) :: ftem(pcols,pver)    ! temporary workspace

    real(r8) :: TS(pcols), SLP(pcols), T700(pcols), Z700(pcols) 
    real(r8) :: z3(pcols,pver)   ! geo-potential height

!    real(r8), pointer :: t_cldtop(:) ! cloud-top temperature estimate: K
    real(r8) :: cf_accum  ! accumulated CF downward 


  
!-----------------------------------------------------------------------
    if (.not.do_cld_diag) return


    if(rk_clouds) then
       dosw     = radiation_do('sw')      ! do shortwave heating calc this timestep?
       dolw     = radiation_do('lw')      ! do longwave heating calc this timestep?
    else
       dosw     = .true.
       dolw     = .true.
    endif

    if (.not.(dosw .or. dolw)) return

    ncol  = state%ncol
    lchnk = state%lchnk

    itim_old = pbuf_old_tim_idx()
    call pbuf_get_field(pbuf, cld_idx, cld, start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )

    if(mg_clouds)then

       call pbuf_get_field(pbuf, iclwp_idx, iclwp )
       call pbuf_get_field(pbuf, iciwp_idx, iciwp )
       call pbuf_get_field(pbuf, dei_idx, dei )
       call pbuf_get_field(pbuf, mu_idx, mu )
       call pbuf_get_field(pbuf, lambda_idx, lambda )

       ! YQIN 10/23/22
       call pbuf_get_field(pbuf, icwnc_idx, icwnc)
       call pbuf_get_field(pbuf, ccn02_idx,  ccn02)
       call pbuf_get_field(pbuf, cdnumc_avg_idx, cdnumc_avg)
       call pbuf_get_field(pbuf, cdnumc_925_idx, cdnumc_925)
       call pbuf_get_field(pbuf, cdnumc_wavg_idx, cdnumc_wavg)
       call pbuf_get_field(pbuf, rel_idx, rel)
       call pbuf_get_field(pbuf, lts_idx,   LTS)
       call pbuf_get_field(pbuf, rh750_idx, RH750)
       call pbuf_get_field(pbuf, eis_idx,   EIS)
!       call pbuf_get_field(pbuf, t_cldtop_idx, t_cldtop)


       call pbuf_get_field(pbuf, cllow_idx, cllow)
       call pbuf_get_field(pbuf, ticlwp_idx, ticlwp)
       call pbuf_get_field(pbuf, cliqwp_idx, cliqwp)
       call pbuf_get_field(pbuf, ticiwp_idx, ticiwp)

       call outfld('dei_cloud',dei(:,:),pcols,lchnk)
       call outfld('mu_cloud',mu(:,:),pcols,lchnk)
       call outfld('lambda_cloud',lambda(:,:),pcols,lchnk)

    elseif(rk_clouds) then

       call pbuf_get_field(pbuf, rei_idx, rei )
       call pbuf_get_field(pbuf, rel_idx, rel )

       call outfld('rel_cloud', rel, pcols, lchnk)
       call outfld('rei_cloud', rei, pcols, lchnk)

       if (cldemis_idx>0) then
          call pbuf_get_field(pbuf, cldemis_idx, cldemis )
       else
          allocate(cldemis(pcols,pver))
       endif
       if (cldtau_idx>0) then
          call pbuf_get_field(pbuf, cldtau_idx, cldtau )
       else
          allocate(cldtau(pcols,pver))
       endif

    endif

    if (cicewp_idx>0) then
       call pbuf_get_field(pbuf, cicewp_idx, cicewp )
    else
       allocate(cicewp(pcols,pver))
    endif
    ! YQIN 12/07/22
    !if (cliqwp_idx>0) then
    !   call pbuf_get_field(pbuf, cliqwp_idx, cliqwp )
    !else
    !   allocate(cliqwp(pcols,pver))
    !endif

    if (nmxrgn_idx>0) then
       call pbuf_get_field(pbuf, nmxrgn_idx, nmxrgn )
    else
       allocate(nmxrgn(pcols))
    endif

    if (pmxrgn_idx>0) then
       call pbuf_get_field(pbuf, pmxrgn_idx, pmxrgn )
    else
       allocate(pmxrgn(pcols,pverp))
    endif

! Compute liquid and ice water paths
    if(mg_clouds) then

       ! ----------------------------------------------------------- !
       ! Adjust in-cloud water values to take account of convective  !
       ! in-cloud water. It is used to calculate the values of       !
       ! iclwp and iciwp to pass to the radiation.                   !
       ! ----------------------------------------------------------- !
       if( conv_water_in_rad /= 0 ) then
          allcld_ice(:ncol,:) = 0._r8 ! Grid-avg all cloud liquid
          allcld_liq(:ncol,:) = 0._r8 ! Grid-avg all cloud ice
    
          call conv_water_4rad( state, pbuf, conv_water_in_rad, allcld_liq, allcld_ice )
       else
          allcld_liq(:ncol,top_lev:pver) = state%q(:ncol,top_lev:pver,ixcldliq)  ! Grid-ave all cloud liquid
          allcld_ice(:ncol,top_lev:pver) = state%q(:ncol,top_lev:pver,ixcldice)  !           "        ice
       end if

       ! ------------------------------------------------------------ !
       ! Compute in cloud ice and liquid mixing ratios                !
       ! Note that 'iclwp, iciwp' are used for radiation computation. !
       ! ------------------------------------------------------------ !


       iciwp = 0._r8
       iclwp = 0._r8
       icimr = 0._r8
       icwmr = 0._r8
       iwc = 0._r8
       lwc = 0._r8

       do k = top_lev, pver
          do i = 1, ncol
             ! Limits for in-cloud mixing ratios consistent with MG microphysics
             ! in-cloud mixing ratio maximum limit of 0.005 kg/kg
             icimr(i,k)     = min( allcld_ice(i,k) / max(0.0001_r8,cld(i,k)),0.005_r8 )
             icwmr(i,k)     = min( allcld_liq(i,k) / max(0.0001_r8,cld(i,k)),0.005_r8 )
             iwc(i,k)       = allcld_ice(i,k) * state%pmid(i,k) / (287.15_r8*state%t(i,k))
             lwc(i,k)       = allcld_liq(i,k) * state%pmid(i,k) / (287.15_r8*state%t(i,k))
             ! Calculate total cloud water paths in each layer
             iciwp(i,k)     = icimr(i,k) * state%pdel(i,k) / gravit
             iclwp(i,k)     = icwmr(i,k) * state%pdel(i,k) / gravit
          end do
       end do

       do k=1,pver
          do i = 1,ncol
             gicewp(i,k) = iciwp(i,k)*cld(i,k)
             gliqwp(i,k) = iclwp(i,k)*cld(i,k)
             cicewp(i,k) = iciwp(i,k)
             cliqwp(i,k) = iclwp(i,k)
          end do
       end do

    elseif(rk_clouds) then

       if (conv_water_in_rad /= 0) then
          call conv_water_4rad(state,pbuf,conv_water_in_rad,allcld_liq,allcld_ice)
       else
          allcld_liq = state%q(:,:,ixcldliq)
          allcld_ice = state%q(:,:,ixcldice)
       end if
    
       do k=1,pver
          do i = 1,ncol
             gicewp(i,k) = allcld_ice(i,k)*state%pdel(i,k)/gravit*1000.0_r8  ! Grid box ice water path.
             gliqwp(i,k) = allcld_liq(i,k)*state%pdel(i,k)/gravit*1000.0_r8  ! Grid box liquid water path.
             cicewp(i,k) = gicewp(i,k) / max(0.01_r8,cld(i,k))               ! In-cloud ice water path.
             cliqwp(i,k) = gliqwp(i,k) / max(0.01_r8,cld(i,k))               ! In-cloud liquid water path.
             ficemr(i,k) = allcld_ice(i,k) / max(1.e-10_r8,(allcld_ice(i,k) + allcld_liq(i,k)))
          end do
       end do
    endif

! Determine parameters for maximum/random overlap
    call cldovrlap(lchnk, ncol, state%pint, cld, nmxrgn, pmxrgn)

! Cloud cover diagnostics
! YQIN 11/15/22 add cltot and cllow as output variables 
    call cloud_cover_diags_out(lchnk, ncol, cld, state%pmid, nmxrgn, pmxrgn, cltot, cllow )
    
    tgicewp(:ncol) = 0._r8
    tgliqwp(:ncol) = 0._r8

    do k=1,pver
       tgicewp(:ncol)  = tgicewp(:ncol) + gicewp(:ncol,k)
       tgliqwp(:ncol)  = tgliqwp(:ncol) + gliqwp(:ncol,k)
    end do

    ! YQIN 10/23/22
    ! calculate vertically-integrated in-cloud liquid/ice water path
    ticlwp(:ncol) = 0._r8
    ticiwp(:ncol) = 0._r8
    do k=1,pver
       ticlwp(:ncol) = ticlwp(:ncol) + cliqwp(:ncol,k)
       ticiwp(:ncol) = ticiwp(:ncol) + cicewp(:ncol,k)
    end do 
    ! convert from kg/m2 to g/m2
    ticlwp(:ncol) = ticlwp(:ncol)*1.0e3
    ticiwp(:ncol) = ticiwp(:ncol)*1.0e3

    ! calculate cloud layer mean in-cloud REL
    trel(:ncol) = 0._r8
    do i=1,ncol
        j = 0
        do k=1,pver
            if (state%pmid(i,k) > 70000._r8 .and. rel(i,k) > 0._r8 .and. cliqwp(i,k) > qsmall)then
                trel(i) = trel(i) + rel(i,k)
                j = j + 1
            end if 
        end do 

        if (j > 0)then
            trel(i) = trel(i)/j
        else
            trel(i) = 0._r8
        end if

    end do

    call outfld( 'ICLWP'   , ticlwp,      pcols, lchnk )
    call outfld( 'ICIWP'   , ticiwp,      pcols, lchnk )
    call outfld( 'TREL'     , trel,        pcols, lchnk)

!    ! YQIN 12/22/22 derive cloud-top temperature
!    t_cldtop = 0._r8
!    do i = 1,ncol
!        cf_accum = 0._r8 ! accumulated CF downward
!        do k = 1, pver
!            if (icwmr(i,k) > 1.e-7_r8) then ! only consider liquid
!                cf_accum = cf_accum + cld(i,k)
!            end if 
!
!            if (cf_accum .ge. 1.0_r8) then
!                do j = 1, k
!                    if (icwmr(i,j) > 1.e-7_r8) then
!                        t_cldtop(i) = t_cldtop(i) + cld(i,j) * state%t(i,j)/cf_accum
!                    end if 
!                end do ! j
!                !write(*,*) 'k=',k,'cf_accum=',cf_accum,'cld=',cld(i,k),'t(k)=',state%t(i,k),'t(pver)=',state%t(i,pver),'t_cldtop=',t_cldtop(i)
!                exit
!            end if 
!
!            if ((k.eq.pver) .and. (cf_accum.lt.1.0_r8)) then
!                t_cldtop(i) = fillvalue
!            end if 
!        end do ! k
!    end do ! i
!
!    call outfld('T_CLDTOP', t_cldtop, pcols, lchnk)



    do itype=1,N_CAT
        if (itype.eq.1)then ! all clouds
            cld_h = cld
            cllow_h = cllow
            rel_h = rel
            trel_h = trel 
            ccn02_h = ccn02
            ccn02_sfc_h = ccn02(:,pver)
            icwnc_h = icwnc
            cdnumc_avg_h = cdnumc_avg 
            cdnumc_925_h = cdnumc_925
            cliqwp_h = cliqwp
            ticlwp_h = ticlwp
        else
            do i=1,pcols
                if (ticiwp(i)>qsmall*1.e3_r8) then ! columns with ice water 
                    cld_h(i,:) = 0._r8
                    cllow_h(i) = 0._r8
                    rel_h(i,:) = 0._r8
                    trel_h(i)  = 0._r8
                    ccn02_h(i,:) = 0._r8
                    ccn02_sfc_h(i) = 0._r8
                    icwnc_h(i,:) = 0._r8
                    cdnumc_avg_h(i) = 0._r8
                    cdnumc_925_h(i) = 0._r8
                    cliqwp_h(i,:) = 0._r8
                    ticlwp_h(i) = 0._r8
                else
                    cld_h(i,:)  = cld(i,:)
                    cllow_h(i)  = cllow(i)
                    rel_h(i,:)   = rel(i,:)
                    trel_h(i)   = trel(i)
                    ccn02_h(i,:) = ccn02(i,:)
                    ccn02_sfc_h(i) = ccn02(i,pver)
                    icwnc_h(i,:) = icwnc(i,:)
                    cdnumc_avg_h(i) = cdnumc_avg(i)
                    cdnumc_925_h(i) = cdnumc_925(i)
                    cliqwp_h(i,:) = cliqwp(i,:)
                    ticlwp_h(i) = ticlwp(i)
                end if
            end do ! i
        end if

        LTSin = EIS
        LTSin_threshold_dry = EIS_threshold_dry
        LTSin_threshold_wet = EIS_threshold_wet

        ! YQIN 10/23/22
        ! CF PDF 
        call pdf1d_regime(fillvalue,ncfs,bcfs_1d,N_REGIME,&
                    LTSin,RH750,LTSin_threshold_dry,LTSin_threshold_wet,RH750_threshold,mincf,&
                    .false.,&
                    cllow_h,&
                    pdf_cf, &
                    cld_h,&
                    pdf_cfp)

        do ireg = 1, N_REGIME
            call outfld( 'PDF_CF'//cats(itype)//regime(ireg),          pdf_cf(:,ireg,:),   pcols, lchnk)
!            call outfld('PDF_CFp'//cats(itype)//regime(ireg),          pdf_cfp(:,ireg,:),  pcols, lchnk)
        end do 

        ! REL PDF 
        call pdf1d_regime(fillvalue,nrels,brels_1d,N_REGIME,&
                    LTSin,RH750,LTSin_threshold_dry,LTSin_threshold_wet,RH750_threshold,minrel,&
                    .false.,&
                    trel_h,&
                    pdf_rel, &
                    rel_h,&
                    pdf_relp)
    
        do ireg = 1, N_REGIME
            call outfld( 'PDF_REL'//cats(itype)//regime(ireg),          pdf_rel(:,ireg,:),   pcols, lchnk)
!            call outfld('PDF_RELp'//cats(itype)//regime(ireg),          pdf_relp(:,ireg,:),  pcols, lchnk)
        end do 

        ! lnCCN PDF
        call pdf1d_regime(fillvalue,nlnccns,blnccns_1d,N_REGIME,&
                    LTSin,RH750,LTSin_threshold_dry,LTSin_threshold_wet,RH750_threshold,minccn,&
                    .true.,&
                    ccn02_sfc_h,&
                    pdf_lnccn, &
                    ccn02_h,&
                    pdf_lnccnp)
    
        do ireg = 1, N_REGIME
            call outfld( 'PDF_lnCCN'//cats(itype)//regime(ireg),          pdf_lnccn(:,ireg,:),   pcols, lchnk)
!            call outfld('PDF_lnCCNp'//cats(itype)//regime(ireg),          pdf_lnccnp(:,ireg,:),  pcols, lchnk)
        end do 
    
        ! lnCDNC PDF
        call pdf1d_regime(fillvalue,nlncdncs,blncdncs_1d,N_REGIME,&
                    LTSin,RH750,LTSin_threshold_dry,LTSin_threshold_wet,RH750_threshold,mincdnc,&
                    .true.,&
                    cdnumc_avg_h,&
                    pdf_lncdnc,&
                    icwnc_h,&
                    pdf_lncdncp)
    
        do ireg = 1, N_REGIME
            call outfld( 'PDF_lnCDNC'//cats(itype)//regime(ireg),          pdf_lncdnc(:,ireg,:),  pcols, lchnk)
!            call outfld('PDF_lnCDNCp'//cats(itype)//regime(ireg),          pdf_lncdncp(:,ireg,:), pcols, lchnk)
        end do 
    
        ! lnCDNC925 PDF
        call pdf1d_regime(fillvalue,nlncdncs,blncdncs_1d,N_REGIME,&
                    LTSin,RH750,LTSin_threshold_dry,LTSin_threshold_wet,RH750_threshold,mincdnc,&
                    .true.,&
                    cdnumc_925_h,&
                    pdf_lncdnc925 &
                    )
    
        do ireg = 1, N_REGIME
            call outfld('PDF_lnCDNC925'//cats(itype)//regime(ireg),          pdf_lncdnc925(:,ireg,:), pcols, lchnk)
        end do 
    
        ! lnTICLWP PDF 
        call pdf1d_regime(fillvalue,nlnticlwps,blnticlwps_1d,N_REGIME,&
                    LTSin,RH750,LTSin_threshold_dry,LTSin_threshold_wet,RH750_threshold,minticlwp,&
                    .true.,&
                    ticlwp_h,&
                    pdf_lnticlwp,&
                    cliqwp_h*1.e3_r8,&
                    pdf_lnticlwpp)
    
        do ireg = 1, N_REGIME
            call outfld( 'PDF_lnICLWP'//cats(itype)//regime(ireg),        pdf_lnticlwp(:,ireg,:), pcols, lchnk)
!            call outfld('PDF_lnICLWPp'//cats(itype)//regime(ireg),        pdf_lnticlwpp(:,ireg,:), pcols, lchnk)
        end do 
    
        ! joint PDF of CCN and CDNC
        call pdf2d_regime(fillvalue,nlnccns,blnccns_1d,nlncdncs,blncdncs_1d,N_REGIME,&
                    LTSin,RH750,LTSin_threshold_dry,LTSin_threshold_wet,RH750_threshold,minccn,mincdnc,&
                    .true.,.true.,&
                    ccn02_sfc_h,cdnumc_avg_h,&
                    pdf_lnccn_lncdnc)
    
        do ireg = 1, N_REGIME
            call outfld('PDF_lnCCN_lnCDNC'//cats(itype)//regime(ireg),          pdf_lnccn_lncdnc(:,ireg,:,:), pcols, lchnk)
        end do 

   
        ! joint PDF of CDNC and REL
        call pdf2d_regime(fillvalue,nlncdncs,blncdncs_1d,nrels,brels_1d,N_REGIME,&
                    LTSin,RH750,LTSin_threshold_dry,LTSin_threshold_wet,RH750_threshold,mincdnc,minrel,&
                    .true.,.false.,&
                    cdnumc_avg_h,trel_h,&
                    pdf_lncdnc_rel)
    
        do ireg = 1, N_REGIME
            call outfld('PDF_lnCDNC_REL'//cats(itype)//regime(ireg),          pdf_lncdnc_rel(:,ireg,:,:), pcols, lchnk)
        end do 
    

        ! 3-D PDF of CDNC, LWP and CF
        call pdf3d_regime(fillvalue,nlncdncs,blncdncs_1d,nlnticlwps,blnticlwps_1d,ncfs,bcfs_1d,N_REGIME,&
                    LTSin,RH750,LTSin_threshold_dry,LTSin_threshold_wet,RH750_threshold,mincdnc,minticlwp,mincf,&
                    .true.,.true.,.false.,&
                    cdnumc_avg_h,ticlwp_h,cllow_h,&
                    pdf_NDLWPCF &
        )
    
        pdf_NDLWPCF_out = 0._r8
        do ireg = 1,1 !N_REGIME
            do i = 1,ncol
                do l = 1,ncfs
                do k = 1,nlnticlwps
                do j = 1,nlncdncs
                    lkj = (l-1)*nlnticlwps*nlncdncs+(k-1)*nlncdncs+j
                    pdf_NDLWPCF_out(i,lkj) = pdf_NDLWPCF(i,ireg,j,k,l)
                end do
                end do
                end do
            end do
    
            !call outfld('PDF_NDLWPCF'//regime(ireg),          pdf_NDLWPCF_out, pcols, lchnk)
        end do 
    
!        ! 2D PDF of CCN and CDNC with level dimensions
!        call pdf2dp_regime(fillvalue,nlnccns,blnccns_1d,nlncdncs,blncdncs_1d,N_REGIME,&
!                     LTSin,RH750,LTSin_threshold_dry,LTSin_threshold_wet,RH750_threshold,minccn,mincdnc,&
!                     .true.,.true.,&
!                     ccn02, icwnc,&
!                     pdf_lnccnp_lncdncp &
!                     )
!        do ireg = 1,N_REGIME
!            call outfld('PDF_lnCCNp_lnCDNCp'//regime(ireg), pdf_lnccnp_lncdncp(:,ireg,:,:), pcols, lchnk)
!        end do 
!    
!        ! 2D PDF of CDNC and REL with level dimensions
!        call pdf2dp_regime(fillvalue,nlncdncs,blncdncs_1d,nrels,brels_1d,N_REGIME,&
!                     LTSin,RH750,LTSin_threshold_dry,LTSin_threshold_wet,RH750_threshold,minccn,mincdnc,&
!                     .true.,.false.,&
!                     icwnc,rel,&
!                     pdf_lncdncp_relp &
!                     )
!        do ireg = 1,N_REGIME
!            call outfld('PDF_lnCDNCp_RELp'//regime(ireg), pdf_lncdncp_relp(:,ireg,:,:), pcols, lchnk)
!        end do 
    
    end do ! itype

    ! ===================================================================== !


    tgwp(:ncol) = tgicewp(:ncol) + tgliqwp(:ncol)
    gwp(:ncol,:pver) = gicewp(:ncol,:pver) + gliqwp(:ncol,:pver)
    cwp(:ncol,:pver) = cicewp(:ncol,:pver) + cliqwp(:ncol,:pver)

    if(rk_clouds) then

       ! Cloud emissivity.
       call cldems(lchnk, ncol, cwp, ficemr, rei, cldemis, cldtau)
       
       ! Effective cloud cover
       do k=1,pver
          do i=1,ncol
             effcld(i,k) = cld(i,k)*cldemis(i,k)
          end do
       end do
       
       call outfld('EFFCLD'  ,effcld , pcols,lchnk)
       call outfld('EMISCLD' ,cldemis, pcols,lchnk)

    else if (mg_clouds) then

       ! --------------------------------------------- !
       ! General outfield calls for microphysics       !
       ! --------------------------------------------- !

       call outfld( 'IWC'      , iwc,         pcols, lchnk )
       call outfld( 'LWC'      , lwc,         pcols, lchnk )
       call outfld( 'ICIMR'    , icimr,       pcols, lchnk )
       call outfld( 'ICWMR'    , icwmr,       pcols, lchnk )

    endif

    call outfld('GCLDLWP' ,gwp    , pcols,lchnk)
    call outfld('TGCLDCWP',tgwp   , pcols,lchnk)
    call outfld('TGCLDLWP',tgliqwp, pcols,lchnk)
    call outfld('TGCLDIWP',tgicewp, pcols,lchnk)
    call outfld('ICLDTWP' ,cwp    , pcols,lchnk)
    call outfld('ICLDIWP' ,cicewp , pcols,lchnk)

! Compute total preciptable water in column (in mm)
    tpw(:ncol) = 0.0_r8
    rgrav = 1.0_r8/gravit
    do k=1,pver
       do i=1,ncol
          tpw(i) = tpw(i) + state%pdel(i,k)*state%q(i,k,1)*rgrav
       end do
    end do

! Diagnostic liquid water path (old specified form)

    call cldclw(lchnk, ncol, state%zi, clwpold, tpw, hl)
    call outfld('SETLWP'  ,clwpold, pcols,lchnk)
    call outfld('LWSH'    ,hl     , pcols,lchnk)
    
    if(rk_clouds) then
       if (cldemis_idx<0) deallocate(cldemis)
       if (cldtau_idx<0) deallocate(cldtau)
    endif
    if (cicewp_idx<0) deallocate(cicewp)
    ! YQIN 12/07/22
    !if (cliqwp_idx<0) deallocate(cliqwp)
    if (pmxrgn_idx<0) deallocate(pmxrgn)
    if (nmxrgn_idx<0) deallocate(nmxrgn)

    return
end subroutine cloud_diagnostics_calc


end module cloud_diagnostics
