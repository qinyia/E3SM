
module cdiag_pdf 

! ----------------------------------------
! Purpuse: diagnose N-D PDF and heatmap

! Author: Yi Qin (yi.qin@pnnl.gov)

! Jan 14, 2022
! ----------------------------------------
   use shr_kind_mod,  only: r8=>shr_kind_r8
   use ppgrid,        only: pcols, pver,pverp
   use cam_history,   only: addfld, horiz_only, add_default
   use cam_history_support, only: fillvalue
   use infnan,              only: isnan, isinf
   use cam_abortutils,      only: endrun



   implicit none
   private
   save

   public :: cdiag_pdf_init
   public :: cdiag_pdf_register 

   ! YQIN 11/23/22
   public :: pdf1d_regime, pdf2d_regime, pdf3d_regime, pdf2dp_regime
   public :: check_nan_inf

   ! YQIN 12/13/22
   integer, public, parameter :: N_REGIME = 5

   ! All; d1: LTS<15K,RH750<50%; d2: LTS<15K,RH750>50%; d3: LTS>15K,RH750<50%; d4: LTS>15K,RH750>50%
   character(len=3), public :: regime(1:N_REGIME) =(/'   ','_d1','_d2','_d3','_d4'/)

   integer, public, parameter :: N_CAT = 2 ! N_CAT=1: all clouds; N_CAT=2: only liquid clouds
   character(len=2), public :: cats(1:N_CAT) = (/"_A","_L"/)


   ! YQIN 10/23/22 --- define bins
   ! cloud fraction: every 0.1 fraction from 0 to 1  
   integer,  parameter, public :: ncfs = 13           ! number of bins
   real(r8)           , public :: bcfs_1d(ncfs+1)     ! bin bounds
   real(r8), target   , public :: bcfs(2,ncfs)        ! left and right bin bounds 
   real(r8), target   , public :: mcfs(ncfs)          ! bin center values

   ! effective radius: every 1 micron from 0 to 30
   integer, parameter , public :: nrels = 30
   real(r8)           , public :: brels_1d(nrels+1) 
   real(r8), target   , public :: brels(2,nrels)
   real(r8), target   , public :: mrels(nrels)

   ! lnCCN
   integer, parameter , public :: nlnccns = 30
   real(r8)           , public :: blnccns_1d(nlnccns+1)
   real(r8), target   , public :: blnccns(2,nlnccns)
   real(r8), target   , public :: mlnccns(nlnccns)

   ! lnCDNC
   integer, parameter , public :: nlncdncs = 30
   real(r8)           , public :: blncdncs_1d(nlncdncs+1)
   real(r8), target   , public :: blncdncs(2,nlncdncs)
   real(r8), target   , public :: mlncdncs(nlncdncs)

   ! ln cloud LWP
   integer, parameter , public :: nlnticlwps = 33  !12/20/22: change from 40 to 33
   real(r8)           , public :: blnticlwps_1d(nlnticlwps+1)
   real(r8), target   , public :: blnticlwps(2,nlnticlwps)
   real(r8), target   , public :: mlnticlwps(nlnticlwps)

   ! unified coordniate for CDNC-LWP-CF 
   integer, parameter , public :: n3ds = nlncdncs*nlnticlwps*ncfs
   real(r8),            public :: b3ds_1d(n3ds+1)
   real(r8), target   , public :: b3ds(2,n3ds)
   real(r8), target   , public :: m3ds(n3ds)

   ! COD 
   integer, parameter , public :: ncods = 30
   real(r8),            public :: bcods_1d(ncods+1)
   real(r8), target   , public :: bcods(2,ncods)
   real(r8), target   , public :: mcods(ncods)


   ! YQIN 
   real(r8), parameter, public :: minccn = 1._r8
   real(r8), parameter, public :: mincdnc = 1._r8
   real(r8), parameter, public :: mincf = 0._r8
   real(r8), parameter, public :: minticlwp = 1._r8
   real(r8), parameter, public :: minrel = 0._r8
   real(r8), parameter, public :: mincod = 0._r8

   real(r8), parameter :: limiter = -1.e6_r8

   real(r8), parameter, public :: LTS_threshold = 18.5 !15 ! K
   real(r8), parameter, public :: RH750_threshold = 50 ! %
   real(r8), parameter, public :: EIS_threshold = 1 !4 !5 !7 ! K

contains

! ================================================================
subroutine cdiag_pdf_register 

    ! YQIN 10/23/22
    use cam_history_support, only: add_hist_coord

       ! YQIN 10/23/22 add new coordinates 
       call add_hist_coord('hist_rel',      nrels,      'REL for histogram',     'micron',   mrels,     bounds_name='hist_rel_bnds',       bounds=brels)
       call add_hist_coord('hist_cf',       ncfs,       'CF for histogram',      'fraction', mcfs,      bounds_name='hist_cf_bnds',        bounds=bcfs)
       call add_hist_coord('hist_lnccn',    nlnccns,    'lnCCN for histogram',   'cm-3',     mlnccns,   bounds_name='hist_lnccn_bnds',     bounds=blnccns)
       call add_hist_coord('hist_lncdnc',   nlncdncs,   'lnCDNC for histogram',  'cm-3',     mlncdncs,  bounds_name='hist_lncdnc_bnds',    bounds=blncdncs)
       call add_hist_coord('hist_lnticlwp', nlnticlwps, 'lnICLWP for histogram', 'kg/m2',    mlnticlwps,bounds_name='hist_lnticlwp_bnds',  bounds=blnticlwps)
       ! YQIN 12/17/22 define a unified coordinate for 3-var (CDNC,LWP,CF) histogram/heatmap
       call add_hist_coord('hist_3d',       n3ds,       '3d coord for histogram', '',        m3ds,      bounds_name='hist_3d_bnds',        bounds=b3ds)
       call add_hist_coord('hist_cod',      ncods,      'COD coord for histogram', '',       mcods,     bounds_name='hist_cod_bnds',       bounds=bcods)

end subroutine cdiag_pdf_register 

! -----------------------------------------------------------------------

subroutine cdiag_pdf_init()

    ! YQIN 10/23/22
    integer  :: k
    real(r8) :: dbin
    integer  :: ireg
    integer  :: itype

       ! YQIN 12/19/22 -- output PDF histogram
       do itype=1,N_CAT
           do ireg=1,N_REGIME
               call addfld(       'PDF_CF'//cats(itype)//regime(ireg),       (/'hist_cf'/),      'A', '', 'PDF of CF',        flag_xyfill=.true., fill_value=fillvalue)
               call addfld(      'PDF_REL'//cats(itype)//regime(ireg),       (/'hist_rel'/),     'A', '', 'PDF of REL',       flag_xyfill=.true., fill_value=fillvalue)
               call addfld(    'PDF_lnCCN'//cats(itype)//regime(ireg),       (/'hist_lnccn'/),   'A', '', 'PDF of lnCCN',     flag_xyfill=.true., fill_value=fillvalue)
               call addfld(   'PDF_lnCDNC'//cats(itype)//regime(ireg),       (/'hist_lncdnc'/),  'A', '', 'PDF of lnCDNC',    flag_xyfill=.true., fill_value=fillvalue)
               call addfld('PDF_lnCDNC925'//cats(itype)//regime(ireg),       (/'hist_lncdnc'/),  'A', '', 'PDF of lnCDNC925', flag_xyfill=.true., fill_value=fillvalue)
               !call addfld('PDF_lnCDNCWAVG'//cats(itype)//regime(ireg),      (/'hist_lncdnc'/),  'A', '', 'PDF of lnCDNCWAVG',flag_xyfill=.true., fill_value=fillvalue)
               call addfld(  'PDF_lnICLWP'//cats(itype)//regime(ireg),       (/'hist_lnticlwp'/),'A', '', 'PDF of lnICLWP',   flag_xyfill=.true., fill_value=fillvalue)

               call addfld('PDF_lnCCN_lnCDNC'//cats(itype)//regime(ireg),    (/'hist_lnccn','hist_lncdnc'/),    'A', '', 'joint histogram of lnCCN and lnCDNC', flag_xyfill=.true., fill_value=fillvalue)
               call addfld(  'PDF_lnCDNC_REL'//cats(itype)//regime(ireg),    (/'hist_lncdnc', 'hist_rel'/),     'A', '', 'joint histogram of lnCDNC and REL', flag_xyfill=.true., fill_value=fillvalue)
               call addfld( 'PDFR_lnCCN_ALBA'//cats(itype)//regime(ireg),    (/'hist_lnccn','hist_cf'/),        'A', '', 'joint histogram of ALBA and lnCCN', flag_xyfill=.true., fill_value=fillvalue)
               call addfld( 'PDFR_lnCCN_ALBC'//cats(itype)//regime(ireg),    (/'hist_lnccn','hist_cf'/),        'A', '', 'joint histogram of ALBC and lnCCN', flag_xyfill=.true., fill_value=fillvalue)
               call addfld( 'PDFR_lnCCN_CODA'//cats(itype)//regime(ireg),    (/'hist_lnccn','hist_cod'/),        'A', '', 'joint histogram of CODA and lnCCN', flag_xyfill=.true., fill_value=fillvalue)


               !call addfld(  'PDF_NDLWPCF'//cats(itype)//regime(ireg),      (/'hist_3d'/),      'A', '', 'PDF of CDNC, LWP and CF')
               call addfld(   'PDFR_NDLWPCF'//cats(itype)//regime(ireg),     (/'hist_3d'/),      'A', '', 'PDF of CDNC, LWP and CF in rad step', flag_xyfill=.true., fill_value=fillvalue)
               call addfld(   'ALBA_NDLWPCF'//cats(itype)//regime(ireg),     (/'hist_3d'/),      'A', '', 'ALBA in CDNC, LWP and CF space in rad step', flag_xyfill=.true., fill_value=fillvalue)
               call addfld(   'ALBC_NDLWPCF'//cats(itype)//regime(ireg),     (/'hist_3d'/),      'A', '', 'ALBC in CDNC, LWP and CF space in rad step', flag_xyfill=.true., fill_value=fillvalue)
               call addfld(   'CODA_NDLWPCF'//cats(itype)//regime(ireg),     (/'hist_3d'/),      'A', '', 'grid-mean COD in CDNC, LWP and CF space in rad step', flag_xyfill=.true., fill_value=fillvalue)

!               call addfld(     'PDF_CFp'//cats(itype)//regime(ireg),       (/'hist_cf'/) ,                    'A', '', 'PDF of CFp')
!               call addfld(    'PDF_RELp'//cats(itype)//regime(ireg),       (/'hist_rel'/) ,                   'A', '', 'PDF of RELp')
!               call addfld(  'PDF_lnCCNp'//cats(itype)//regime(ireg),       (/'hist_lnccn'/) ,                 'A', '', 'PDF of lnCCNp')
!               call addfld( 'PDF_lnCDNCp'//cats(itype)//regime(ireg),       (/'hist_lncdnc'/) ,                'A', '', 'PDF of lnCDNCp')
!               call addfld('PDF_lnICLWPp'//cats(itype)//regime(ireg),       (/'hist_lnticlwp'/),               'A', '', 'PDF of lnICLWPp')
!    
!               call addfld('PDF_lnCCNp_lnCDNCp'//cats(itype)//regime(ireg), (/'hist_lnccn','hist_lncdnc'/),    'A', '', 'joint histogram of lnCCNp and lnCDNCp')
!               call addfld('PDF_lnCDNCp_RELp'//cats(itype)//regime(ireg),   (/'hist_lncdnc', 'hist_rel'/),     'A', '', 'joint histogram of lnCDNCp and RELp')
!    
!               call addfld('PDFR_NDLWPCFp'//cats(itype)//regime(ireg),      (/'hist_3d'/),                     'A', '', 'PDF of CDNCp, LWPp and CFp in rad step')
!               call addfld('CODA_NDLWPCFp'//cats(itype)//regime(ireg),      (/'hist_3d'/),                     'A', '', 'grid-mean COD in CDNCp, LWPp and CFp space in rad step', flag_xyfill=.true., fill_value=fillvalue)

           end do ! ireg
       end do ! itype

       do itype = 1,N_CAT
       do ireg = 1, N_REGIME
           call add_default (          'PDF_CF'//cats(itype)//regime(ireg),      1, ' ')
           call add_default (         'PDF_REL'//cats(itype)//regime(ireg),      1, ' ')
           call add_default (       'PDF_lnCCN'//cats(itype)//regime(ireg),      1, ' ')
           call add_default (      'PDF_lnCDNC'//cats(itype)//regime(ireg),      1, ' ')
           call add_default (   'PDF_lnCDNC925'//cats(itype)//regime(ireg),      1, ' ')
!           call add_default (  'PDF_lnCDNCWAVG'//cats(itype)//regime(ireg),      1, ' ')
           call add_default (     'PDF_lnICLWP'//cats(itype)//regime(ireg),      1, ' ')

           call add_default ('PDF_lnCCN_lnCDNC'//cats(itype)//regime(ireg),      1, ' ')
           call add_default (  'PDF_lnCDNC_REL'//cats(itype)//regime(ireg),      1, ' ')
           call add_default ( 'PDFR_lnCCN_ALBA'//cats(itype)//regime(ireg),      1, ' ')
           call add_default ( 'PDFR_lnCCN_ALBC'//cats(itype)//regime(ireg),      1, ' ')
           call add_default ( 'PDFR_lnCCN_CODA'//cats(itype)//regime(ireg),      1, ' ')

           !call add_default (   'PDF_NDLWPCF'//cats(itype)//regime(ireg),       1, ' ')
           call add_default (   'PDFR_NDLWPCF'//cats(itype)//regime(ireg),       1, ' ')
           call add_default (   'ALBA_NDLWPCF'//cats(itype)//regime(ireg),       1, ' ')
           call add_default (   'ALBC_NDLWPCF'//cats(itype)//regime(ireg),       1, ' ')
           call add_default (   'CODA_NDLWPCF'//cats(itype)//regime(ireg),       1, ' ')

!           call add_default (      'PDF_CFp'//cats(itype)//regime(ireg),        1, ' ')
!           call add_default (     'PDF_RELp'//cats(itype)//regime(ireg),        1, ' ')
!           call add_default (   'PDF_lnCCNp'//cats(itype)//regime(ireg),        1, ' ')
!           call add_default (  'PDF_lnCDNCp'//cats(itype)//regime(ireg),        1, ' ')
!           call add_default ( 'PDF_lnICLWPp'//cats(itype)//regime(ireg),        1, ' ')

!           call add_default ('PDF_lnCCNp_lnCDNCp'//cats(itype)//regime(ireg),   1, ' ')
!           call add_default ('PDF_lnCDNCp_RELp'//cats(itype)//regime(ireg),     1, ' ')

!           call add_default ('PDFR_NDLWPCFp'//cats(itype)//regime(ireg),       1, ' ')
!           call add_default ('CODA_NDLWPCFp'//cats(itype)//regime(ireg),       1, ' ')
       end do
       end do

       ! YQIN -- define bins 
       ! cloud fraction 
       ! use -1~0 as the first bin, 0-0.01 as the second bin, 0.01-0.10 as the third bin, and the
       ! following bins with 0.10 as the bin size.
       ! [-1.0, 0, 0.01, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00] - 13 bin edges
       !mcfs = [-0.5         0.005       0.050       0.150       0.250       0.350
       !         0.450       0.550       0.650       0.750       0.850       0.950 ] - 12 bins
       dbin = 0.10_r8
       bcfs_1d(1) = -1.0_r8
       bcfs_1d(2) = 0.00_r8
       bcfs_1d(3) = 0.01_r8
       bcfs_1d(4) = 0.10_r8
       mcfs(1) = -0.5_r8
       mcfs(2) = 0.005_r8
       mcfs(3) = 0.050_r8
       mcfs(4) = 0.150_r8
       do k=5, ncfs
           mcfs(k) = mcfs(k-1) + dbin
       end do 
       do k=5,ncfs+1
             bcfs_1d(k) = bcfs_1d(k-1) + dbin 
       end do
       do k=1,ncfs
             bcfs(1,k)=bcfs_1d(k)
             bcfs(2,k)=bcfs_1d(k+1)
       end do    

       ! REL
       dbin = 1.0_r8
       brels_1d(1) = -1.0_r8
       brels_1d(2) = 0.0_r8
       mrels(1) = -0.5_r8
       mrels(2) = dbin/2._r8
       do k=3, nrels
           mrels(k) = mrels(k-1) + dbin
       end do 
       do k=3,nrels+1
             brels_1d(k) = brels_1d(k-1) + dbin 
       end do
       do k=1,nrels
             brels(1,k)=brels_1d(k)
             brels(2,k)=brels_1d(k+1)
       end do    

       ! lnCCN
       dbin=log(10._r8)/8._r8
       blnccns_1d(1) = -10._r8
       blnccns_1d(2) = 0._r8
       mlnccns(1) = -5._r8
       mlnccns(2) = dbin/2.0_r8
       do k=3, nlnccns
           mlnccns(k) = mlnccns(k-1) + dbin
       end do 
       do k=3,nlnccns+1
             blnccns_1d(k) = blnccns_1d(k-1) + dbin
       end do
       do k=1,nlnccns
             blnccns(1,k)=blnccns_1d(k)
             blnccns(2,k)=blnccns_1d(k+1)
       end do    

       ! lnCDNC
       dbin = log(10._r8)/8._r8
       blncdncs_1d(1) = -10._r8
       blncdncs_1d(2) = 0._r8
       mlncdncs(1) = -5._r8
       mlncdncs(2) = dbin/2._r8
       do k=3, nlncdncs
           mlncdncs(k) = mlncdncs(k-1) + dbin
       end do 
       do k=3,nlncdncs+1
             blncdncs_1d(k) = blncdncs_1d(k-1) + dbin
       end do
       do k=1,nlncdncs
             blncdncs(1,k)=blncdncs_1d(k)
             blncdncs(2,k)=blncdncs_1d(k+1)
       end do    

       ! lnLWP
       dbin = log(10._r8)/8._r8
       blnticlwps_1d(1) = -10._r8
       blnticlwps_1d(2) = 0._r8
       mlnticlwps(1) = -5._r8
       mlnticlwps(2) = dbin/2._r8
       do k=3, nlnticlwps
           mlnticlwps(k) = mlnticlwps(k-1) + dbin
       end do 
       do k=3,nlnticlwps+1
             blnticlwps_1d(k) = blnticlwps_1d(k-1) + dbin
       end do
       do k=1,nlnticlwps
             blnticlwps(1,k)=blnticlwps_1d(k)
             blnticlwps(2,k)=blnticlwps_1d(k+1)
       end do    

       ! unified 3d coordinate
       dbin = 1.0_r8
       b3ds_1d(1) = 0._r8
       m3ds(1) = b3ds_1d(1)+dbin/2._r8
       do k=2, n3ds
           m3ds(k) = m3ds(k-1) + dbin
       end do 
       do k=2,n3ds+1
             b3ds_1d(k) = b3ds_1d(k-1) + dbin
       end do
       do k=1,n3ds
             b3ds(1,k)=b3ds_1d(k)
             b3ds(2,k)=b3ds_1d(k+1)
       end do    

       ! COD coordinate
       dbin = 2.0_r8
       bcods_1d(1) = 0._r8
       mcods(1) = bcods_1d(1)+dbin/2._r8
       do k=2, ncods
           mcods(k) = mcods(k-1) + dbin
       end do 
       do k=2,ncods+1
             bcods_1d(k) = bcods_1d(k-1) + dbin
       end do
       do k=1,ncods
             bcods(1,k)=bcods_1d(k)
             bcods(2,k)=bcods_1d(k+1)
       end do    

end subroutine cdiag_pdf_init

! -----------------------------------------------------------------------
! YQIN 12/13/22
subroutine pdf1d_regime(undefvar,nbins,bvars_1d,N_REGIME,&
                 LTS,RH750,LTS_threshold,RH750_threshold,varmin,&
                 logvar,&
                 var,&
                 pdf_var,&
                 varp, &
                 pdf_varp)

real(r8), intent(in) :: undefvar ! missing value 
integer, intent(in) :: nbins ! number of bins
real(r8), intent(in) :: bvars_1d(nbins+1) ! bounds of bins
integer,  intent(in) :: N_REGIME
real(r8), intent(in) :: LTS(pcols)
real(r8), intent(in) :: RH750(pcols)
real(r8), intent(in) :: LTS_threshold
real(r8), intent(in) :: RH750_threshold
real(r8), intent(in) :: varmin  ! set the minimum of input variable
logical, intent(in)  :: logvar ! flag to denote whether use the log or the linear bins

real(r8), intent(in) :: var(pcols) ! input variable
real(r8), intent(out) :: pdf_var(pcols,N_REGIME,nbins) 

real(r8), optional, intent(in) :: varp(pcols,pver) ! input variable with level 
real(r8), optional, intent(out) :: pdf_varp(pcols,N_REGIME,nbins)


integer i, k, p, w
integer idx,kk

real(r8) :: vartmp
real(r8) :: varptmp

pdf_var = 0._r8
if (present(varp)) then
pdf_varp = 0._r8
end if

do i=1,pcols

    if (var(i) .eq. undefvar) then
        pdf_var(i,:,:) = undefvar 
        cycle 
    end if 

    if (logvar) then
        if (var(i) > varmin) then
            vartmp = log(var(i))
        end if
    else
        vartmp = var(i)
    end if
    
do k=1,nbins

    ! the condition with limiter=-1.e6_r8 is to avoid the negative missing value from
    ! cosp simulator, which has a fillvalue as -1.e30
    if ((var(i) .le. varmin) .and. (var(i).gt.limiter)) then ! set the first bin
        kk = 1
        idx = 1
        pdf_var(i,idx,kk) = 1._r8
        if (LTS(i) .le. LTS_threshold) then
            if (RH750(i) .le. RH750_threshold) then
                idx = 2
            else
                idx = 3
            end if
        else
            if (RH750(i) .le. RH750_threshold) then
                idx = 4
            else
                idx = 5
            end if
        end if
        pdf_var(i,idx,kk) = 1._r8

    else if (vartmp.gt.bvars_1d(k) .and. vartmp.le.bvars_1d(k+1)) then
        kk = k
        idx = 1
        pdf_var(i,idx,kk) = 1._r8
        if (LTS(i) .le. LTS_threshold) then
            if (RH750(i) .le. RH750_threshold) then
                idx = 2
            else
                idx = 3
            end if
        else
            if (RH750(i) .le. RH750_threshold) then
                idx = 4
            else
                idx = 5
            end if
        end if
        pdf_var(i,idx,kk) = 1._r8
    end if

end do ! k
end do ! i

if (present(varp))then
    do i=1,pcols
    
    if (varp(i,0) .eq. undefvar) then
        pdf_varp(i,:,:) = undefvar
        cycle
    end if

    ! include lev dimension
    do p=1,pver

        if (logvar) then
            varptmp = log(varp(i,p))
        else
            varptmp = varp(i,p)
        end if
    
        if ((varp(i,p) .le. varmin) .and. (varp(i,p).gt.limiter)) then
            kk = 1
            idx = 1
            pdf_varp(i,idx,kk) = pdf_varp(i,idx,kk) + 1._r8
           
            if (LTS(i) .le. LTS_threshold) then
                if (RH750(i) .le. RH750_threshold) then
                    idx = 2
                else
                    idx = 3
                end if
            else
                if (RH750(i) .le. RH750_threshold) then
                    idx = 4
                else
                    idx = 5
                end if
            end if
            pdf_varp(i,idx,kk) = pdf_varp(i,idx,kk) + 1._r8
        end if
    
        do k=2,nbins
            if (varptmp.gt.bvars_1d(k) .and. varptmp.le.bvars_1d(k+1)) then
                kk = k
                idx = 1
                pdf_varp(i,idx,kk) = pdf_varp(i,idx,kk) + 1._r8

                if (LTS(i) .le. LTS_threshold) then
                    if (RH750(i) .le. RH750_threshold) then
                        idx = 2
                    else
                        idx = 3
                    end if
                else
                    if (RH750(i) .le. RH750_threshold) then
                        idx = 4
                    else
                        idx = 5
                    end if
                end if
                pdf_varp(i,idx,kk) = pdf_varp(i,idx,kk) + 1._r8
            end if
        end do ! k
    
    end do ! p
    end do ! i
end if

end subroutine pdf1d_regime

! ---------------------------------------------------------------
subroutine pdf2d_regime(undefvar,nbins1,bvars1_1d,nbins2,bvars2_1d,N_REGIME,&
                 LTS,RH750,LTS_threshold,RH750_threshold,var1min,var2min,&
                 logvar1,logvar2,&
                 var1,var2,&
                 pdf_var,&
                 var3,heatmap_var3 &
                 )
! calculate the joint PDF and heatmap 

real(r8), intent(in) :: undefvar ! missing value 
integer, intent(in) :: nbins1 ! number of bins
real(r8), intent(in) :: bvars1_1d(nbins1+1) ! bounds of bins
integer, intent(in) :: nbins2 ! number of bins
real(r8), intent(in) :: bvars2_1d(nbins2+1) ! bounds of bins

integer, intent(in) :: N_REGIME
real(r8), intent(in) :: LTS(pcols) 
real(r8), intent(in) :: RH750(pcols) 
real(r8), intent(in) :: LTS_threshold
real(r8), intent(in) :: RH750_threshold

real(r8), intent(in) :: var1min  ! set the minimum of input variable
real(r8), intent(in) :: var2min  ! set the minimum of input variable

logical,  intent(in) :: logvar1
logical,  intent(in) :: logvar2   ! .true.: log(var2); .false.: var2

real(r8), intent(in) :: var1(pcols) ! input variable
real(r8), intent(in) :: var2(pcols) ! input variable
real(r8), intent(out) :: pdf_var(pcols,N_REGIME,nbins1,nbins2) 


real(r8), optional, intent(in) :: var3(pcols) 
real(r8), optional, intent(out) :: heatmap_var3(pcols,N_REGIME,nbins1,nbins2)

integer i, j, k, p, w
integer idx,jj,kk

real(r8) :: var2tmp, var1tmp

pdf_var = 0._r8

if (present(var3)) then
    heatmap_var3 = undefvar 
end if

do i=1,pcols

    if (var1(i).eq.undefvar .or. var2(i).eq.undefvar) then
        pdf_var(i,:,:,:) = undefvar 
        cycle 
    end if 

    if (present(var3)) then
        if (var3(i).eq.undefvar) then
            heatmap_var3(i,:,:,:) = undefvar
            cycle
        end if
    end if 
        
    if (logvar1) then
        if (var1(i)>0._r8)then
            var1tmp = log(var1(i))
        end if
    else
        var1tmp = var1(i)
    end if

    if (logvar2) then
        if (var2(i)>0._r8)then
            var2tmp = log(var2(i))
        end if
    else
        var2tmp = var2(i)
    end if

    ! check var3 to ensure no nan and inf
    if (present(var3))then
        call check_nan_inf(var3(i),'pdf2d')
    end if

    do j=1,nbins1
    do k=1,nbins2
        if (var1(i) .le. var1min .and. (var1(i).gt.limiter)) then
            if (var2(i) .le. var2min .and. (var2(i).gt.limiter)) then
                ! ===========================================
                idx = 1
                jj = 1
                kk = 1
                pdf_var(i,idx,jj,kk) = 1._r8
                if (present(var3)) then
                    heatmap_var3(i,idx,jj,kk) = var3(i)
                end if

                if (LTS(i) .le. LTS_threshold) then
                    if (RH750(i) .le. RH750_threshold) then
                        idx = 2
                    else
                        idx = 3
                    end if
                else
                    if (RH750(i) .le. RH750_threshold) then
                        idx = 4
                    else
                        idx = 5
                    end if
                end if

                pdf_var(i,idx,jj,kk) = 1._r8
                if (present(var3)) then
                    heatmap_var3(i,idx,jj,kk) = var3(i)
                end if
                ! ===========================================

            else if (var2tmp.gt.bvars2_1d(k) .and. var2tmp.le.bvars2_1d(k+1)) then
                ! ===========================================
                idx = 1
                jj = 1
                kk = k
                pdf_var(i,idx,jj,kk) = 1._r8
                if (present(var3)) then
                    heatmap_var3(i,idx,jj,kk) = var3(i)
                end if

                if (LTS(i) .le. LTS_threshold) then
                    if (RH750(i) .le. RH750_threshold) then
                        idx = 2
                    else
                        idx = 3
                    end if
                else
                    if (RH750(i) .le. RH750_threshold) then
                        idx = 4
                    else
                        idx = 5
                    end if
                end if

                pdf_var(i,idx,jj,kk) = 1._r8
                if (present(var3)) then
                    heatmap_var3(i,idx,jj,kk) = var3(i)
                end if
                ! ===========================================

            end if
        else if (var1tmp.gt.bvars1_1d(j).and.var1tmp.le.bvars1_1d(j+1)) then
            if (var2(i) .le. var2min .and. (var2(i).gt.limiter)) then
                ! ===========================================
                idx = 1
                jj = j
                kk = 1
                pdf_var(i,idx,jj,kk) = 1._r8
                if (present(var3)) then
                    heatmap_var3(i,idx,jj,kk) = var3(i)
                end if
                if (LTS(i) .le. LTS_threshold) then
                    if (RH750(i) .le. RH750_threshold) then
                        idx = 2
                    else
                        idx = 3
                    end if
                else
                    if (RH750(i) .le. RH750_threshold) then
                        idx = 4
                    else
                        idx = 5
                    end if
                end if

                pdf_var(i,idx,jj,kk) = 1._r8
                if (present(var3)) then
                    heatmap_var3(i,idx,jj,kk) = var3(i)
                end if
                ! ===========================================

            else if (var2tmp.gt.bvars2_1d(k) .and. var2tmp.le.bvars2_1d(k+1)) then
                ! ===========================================

                idx = 1
                jj = j
                kk = k
                pdf_var(i,idx,jj,kk) = 1._r8
                if (present(var3)) then
                    heatmap_var3(i,idx,jj,kk) = var3(i)
                end if

                if (LTS(i) .le. LTS_threshold) then
                    if (RH750(i) .le. RH750_threshold) then
                        idx = 2
                    else
                        idx = 3
                    end if
                else
                    if (RH750(i) .le. RH750_threshold) then
                        idx = 4
                    else
                        idx = 5
                    end if
                end if

                pdf_var(i,idx,jj,kk) = 1._r8
                if (present(var3)) then
                    heatmap_var3(i,idx,jj,kk) = var3(i)
                end if
                ! ===========================================

            end if
        end if 
           
    end do ! k
    end do ! j
end do ! i

end subroutine pdf2d_regime

! YQIN 12/17/22
! ---------------------------------------------------------------
subroutine pdf3d_regime(undefvar,nbins1,bvars1_1d,nbins2,bvars2_1d,nbins3,bvars3_1d,N_REGIME,&
                 LTS,RH750,LTS_threshold,RH750_threshold,var1min,var2min,var3min,&
                 logvar1,logvar2,logvar3,&
                 var1,var2,var3,&
                 pdf_var,&
                 varout,heatmap_varout &
                 )
! calculate the joint PDF and heatmap 

real(r8), intent(in) :: undefvar            ! missing value 
integer, intent(in) :: nbins1               ! number of bins
real(r8), intent(in) :: bvars1_1d(nbins1+1) ! bounds of bins
integer, intent(in) :: nbins2               ! number of bins
real(r8), intent(in) :: bvars2_1d(nbins2+1) ! bounds of bins
integer, intent(in) :: nbins3               ! number of bins
real(r8), intent(in) :: bvars3_1d(nbins3+1) ! bounds of bins

integer, intent(in) :: N_REGIME
real(r8), intent(in) :: LTS(pcols) 
real(r8), intent(in) :: RH750(pcols) 
real(r8), intent(in) :: LTS_threshold
real(r8), intent(in) :: RH750_threshold

real(r8), intent(in) :: var1min  ! set the minimum of input variable
real(r8), intent(in) :: var2min  ! set the minimum of input variable
real(r8), intent(in) :: var3min  ! set the minimum of input variable

logical,  intent(in) :: logvar1
logical,  intent(in) :: logvar2   ! .true.: log(var2); .false.: var2
logical,  intent(in) :: logvar3   ! .true.: log(var3); .false.: var3

real(r8), intent(in) :: var1(pcols) ! input variable
real(r8), intent(in) :: var2(pcols) ! input variable
real(r8), intent(in) :: var3(pcols) ! input variable

real(r8), intent(out) :: pdf_var(pcols,N_REGIME,nbins1,nbins2,nbins3) 

real(r8), optional, intent(in) :: varout(pcols) 
real(r8), optional, intent(out) :: heatmap_varout(pcols,N_REGIME,nbins1,nbins2,nbins3)

integer i, j, k, l, p, w
integer idx,jj,kk,ll

real(r8) :: var1tmp,var2tmp,var3tmp

pdf_var = 0._r8

if (present(varout)) then
    heatmap_varout = undefvar 
end if

do i=1,pcols

    if (var1(i).eq.undefvar .or. var2(i).eq.undefvar .or. var3(i).eq.undefvar) then
        pdf_var(i,:,:,:,:) = undefvar 
        cycle 
    end if 

    if (present(varout)) then
        if (varout(i).eq.undefvar) then
            heatmap_varout(i,:,:,:,:) = undefvar 
            cycle 
        end if
    end if

    if (logvar1) then
        if (var1(i)>0._r8)then
            var1tmp = log(var1(i))
        end if
    else
        var1tmp = var1(i)
    end if

    if (logvar2) then
        if (var2(i)>0._r8)then
            var2tmp = log(var2(i))
        end if
    else
        var2tmp = var2(i)
    end if

    if (logvar3) then
        if (var3(i)>0._r8)then
            var3tmp = log(var3(i))
        end if
    else
        var3tmp = var3(i)
    end if

    ! check varout to ensure no nan and inf
    if (present(varout))then
        call check_nan_inf(varout(i),'pdf3d_regime')
    end if

    ! ==================
    do j=1,nbins1
    do k=1,nbins2
    do l=1,nbins3

        if (var1(i) .le. var1min .and. (var1(i).gt.limiter)) then
            if (var2(i) .le. var2min .and. (var2(i).gt.limiter)) then
                if (var3(i) .le. var3min .and. (var3(i).gt.limiter)) then
                    ! ===========================================
                    jj = 1
                    kk = 1
                    ll = 1
                    idx = 1
                    pdf_var(i,idx,jj,kk,ll) = 1._r8
                    if (present(varout)) then
                        heatmap_varout(i,idx,jj,kk,ll) = varout(i)
                    end if

                    if (LTS(i) .le. LTS_threshold) then
                        if (RH750(i) .le. RH750_threshold) then
                            idx = 2
                        else
                            idx = 3
                        end if
                    else
                        if (RH750(i) .le. RH750_threshold) then
                            idx = 4
                        else
                            idx = 5
                        end if ! RH
                    end if ! LTS

                    pdf_var(i,idx,jj,kk,ll) = 1._r8
                    if (present(varout)) then
                        heatmap_varout(i,idx,jj,kk,ll) = varout(i)
                    end if
                    ! ===========================================
                else if (var3tmp.gt.bvars3_1d(l) .and. var3tmp.le.bvars3_1d(l+1)) then
                    ! ===========================================
                    jj = 1
                    kk = 1
                    ll = l
                    idx = 1
                    pdf_var(i,idx,jj,kk,ll) = 1._r8
                    if (present(varout)) then
                        heatmap_varout(i,idx,jj,kk,ll) = varout(i)
                    end if

                    if (LTS(i) .le. LTS_threshold) then
                        if (RH750(i) .le. RH750_threshold) then
                            idx = 2
                        else
                            idx = 3
                        end if ! RH750
                    else
                        if (RH750(i) .le. RH750_threshold) then
                            idx = 4
                        else
                            idx = 5
                        end if ! RH
                    end if ! LTS

                    pdf_var(i,idx,jj,kk,ll) = 1._r8
                    if (present(varout)) then
                        heatmap_varout(i,idx,jj,kk,ll) = varout(i)
                    end if
                    ! ===========================================
                end if ! var3tmp 

            else if (var2tmp.gt.bvars2_1d(k) .and. var2tmp.le.bvars2_1d(k+1)) then
                if (var3(i) .le. var3min .and. (var3(i).gt.limiter)) then
                    ! ===========================================
                    jj = 1
                    kk = k
                    ll = 1
                    idx = 1
                    pdf_var(i,idx,jj,kk,ll) = 1._r8
                    if (present(varout)) then
                        heatmap_varout(i,idx,jj,kk,ll) = varout(i)
                    end if

                    if (LTS(i) .le. LTS_threshold) then
                        if (RH750(i) .le. RH750_threshold) then
                            idx = 2
                        else
                            idx = 3
                        end if
                    else
                        if (RH750(i) .le. RH750_threshold) then
                            idx = 4
                        else
                            idx = 5
                        end if ! RH
                    end if ! LTS

                    pdf_var(i,idx,jj,kk,ll) = 1._r8
                    if (present(varout)) then
                        heatmap_varout(i,idx,jj,kk,ll) = varout(i)
                    end if

                    ! ===========================================
                else if (var3tmp.gt.bvars3_1d(l) .and. var3tmp.le.bvars3_1d(l+1)) then
                    ! ===========================================
                    jj = 1
                    kk = k
                    ll = l
                    idx = 1
                    pdf_var(i,idx,jj,kk,ll) = 1._r8
                    if (present(varout)) then
                        heatmap_varout(i,idx,jj,kk,ll) = varout(i)
                    end if

                    if (LTS(i) .le. LTS_threshold) then
                        if (RH750(i) .le. RH750_threshold) then
                            idx = 2
                        else
                            idx = 3
                        end if ! RH750
                    else
                        if (RH750(i) .le. RH750_threshold) then
                            idx = 4
                        else
                            idx = 5
                        end if ! RH
                    end if ! LTS

                    pdf_var(i,idx,jj,kk,ll) = 1._r8
                    if (present(varout)) then
                        heatmap_varout(i,idx,jj,kk,ll) = varout(i)
                    end if
                    ! ===========================================
                end if ! var3tmp 
            end if ! var2tmp

        else if (var1tmp.gt.bvars1_1d(j) .and. var1tmp.le.bvars1_1d(j+1)) then
            if (var2(i) .le. var2min .and. (var2(i).gt.limiter)) then
                if (var3(i) .le. var3min .and. (var3(i).gt.limiter)) then
                    ! ===========================================
                    jj = j
                    kk = 1
                    ll = 1
                    idx = 1
                    pdf_var(i,idx,jj,kk,ll) = 1._r8
                    if (present(varout)) then
                        heatmap_varout(i,idx,jj,kk,ll) = varout(i)
                    end if

                    if (LTS(i) .le. LTS_threshold) then
                        if (RH750(i) .le. RH750_threshold) then
                            idx = 2
                        else
                            idx = 3
                        end if
                    else
                        if (RH750(i) .le. RH750_threshold) then
                            idx = 4
                        else
                            idx = 5
                        end if ! RH
                    end if ! LTS

                    pdf_var(i,idx,jj,kk,ll) = 1._r8
                    if (present(varout)) then
                        heatmap_varout(i,idx,jj,kk,ll) = varout(i)
                    end if

                    ! ===========================================
                else if (var3tmp.gt.bvars3_1d(l) .and. var3tmp.le.bvars3_1d(l+1)) then
                    ! ===========================================
                    jj = j
                    kk = 1
                    ll = l
                    idx = 1
                    pdf_var(i,idx,jj,kk,ll) = 1._r8
                    if (present(varout)) then
                        heatmap_varout(i,idx,jj,kk,ll) = varout(i)
                    end if

                    if (LTS(i) .le. LTS_threshold) then
                        if (RH750(i) .le. RH750_threshold) then
                            idx = 2
                        else
                            idx = 3
                        end if ! RH750
                    else
                        if (RH750(i) .le. RH750_threshold) then
                            idx = 4
                        else
                            idx = 5
                        end if ! RH
                    end if ! LTS

                    pdf_var(i,idx,jj,kk,ll) = 1._r8
                    if (present(varout)) then
                        heatmap_varout(i,idx,jj,kk,ll) = varout(i)
                    end if
                    ! ===========================================
                end if ! var3tmp 

            else if (var2tmp.gt.bvars2_1d(k) .and. var2tmp.le.bvars2_1d(k+1)) then
                if (var3(i) .le. var3min .and. (var3(i).gt.limiter)) then
                    ! ===========================================
                    jj = j
                    kk = k
                    ll = 1
                    idx = 1
                    pdf_var(i,idx,jj,kk,ll) = 1._r8
                    if (present(varout)) then
                        heatmap_varout(i,idx,jj,kk,ll) = varout(i)
                    end if

                    if (LTS(i) .le. LTS_threshold) then
                        if (RH750(i) .le. RH750_threshold) then
                            idx = 2
                        else
                            idx = 3
                        end if
                    else
                        if (RH750(i) .le. RH750_threshold) then
                            idx = 4
                        else
                            idx = 5
                        end if ! RH
                    end if ! LTS

                    pdf_var(i,idx,jj,kk,ll) = 1._r8
                    if (present(varout)) then
                        heatmap_varout(i,idx,jj,kk,ll) = varout(i)
                    end if

                    ! ===========================================
                else if (var3tmp.gt.bvars3_1d(l) .and. var3tmp.le.bvars3_1d(l+1)) then
                    ! ===========================================
                    jj = j
                    kk = k
                    ll = l
                    idx = 1
                    pdf_var(i,idx,jj,kk,ll) = 1._r8
                    if (present(varout)) then
                        heatmap_varout(i,idx,jj,kk,ll) = varout(i)
                    end if

                    if (LTS(i) .le. LTS_threshold) then
                        if (RH750(i) .le. RH750_threshold) then
                            idx = 2
                        else
                            idx = 3
                        end if ! RH750
                    else
                        if (RH750(i) .le. RH750_threshold) then
                            idx = 4
                        else
                            idx = 5
                        end if ! RH
                    end if ! LTS

                    pdf_var(i,idx,jj,kk,ll) = 1._r8
                    if (present(varout)) then
                        heatmap_varout(i,idx,jj,kk,ll) = varout(i)
                    end if
                    ! ===========================================
                end if ! var3tmp 
            end if ! var2tmp
        end if ! var1tmp
           
    end do ! l
    end do ! k
    end do ! j
end do ! i

end subroutine pdf3d_regime

! ---------------------------------------------------------------
subroutine pdf2dp_regime(undefvar,nbins1,bvars1_1d,nbins2,bvars2_1d,N_REGIME,&
                 LTS,RH750,LTS_threshold,RH750_threshold,var1min,var2min,&
                 logvar1,logvar2,&
                 var1,var2,&
                 pdf_var,&
                 var3,heatmap_var3 &
                 )

! calculate the joint PDF with variables having lev dimension

real(r8), intent(in) :: undefvar  ! missing value 
integer,  intent(in) :: nbins1 ! number of bins
real(r8), intent(in) :: bvars1_1d(nbins1+1) ! bounds of bins
integer,  intent(in) :: nbins2 ! number of bins
real(r8), intent(in) :: bvars2_1d(nbins2+1) ! bounds of bins

integer,  intent(in) :: N_REGIME
real(r8), intent(in) :: LTS(pcols) 
real(r8), intent(in) :: RH750(pcols) 
real(r8), intent(in) :: LTS_threshold
real(r8), intent(in) :: RH750_threshold

real(r8), intent(in) :: var1min  ! set the minimum of input variable
real(r8), intent(in) :: var2min  ! set the minimum of input variable

logical,  intent(in) :: logvar1
logical,  intent(in) :: logvar2   ! .true.: log(var2); .false.: var2

real(r8), intent(in) :: var1(pcols,pver) ! input variable
real(r8), intent(in) :: var2(pcols,pver) ! input variable
real(r8), intent(out) :: pdf_var(pcols,N_REGIME,nbins1,nbins2) 

real(r8), optional, intent(in) :: var3(pcols,pver) 
real(r8), optional, intent(out) :: heatmap_var3(pcols,N_REGIME,nbins1,nbins2)

integer i, j, k, p, w
integer idx, jj, kk
integer ireg

real(r8) :: var1tmp, var2tmp 
real(r8) :: freq_var3(pcols,N_REGIME,nbins1,nbins2)


pdf_var = 0._r8

if (present(var3)) then
    heatmap_var3 = 0._r8
    freq_var3 = 0._r8
end if

do i=1,pcols

if (var1(i,0).eq.undefvar .or. var2(i,0).eq.undefvar) then
    pdf_var(i,:,:,:) = undefvar     
    cycle
end if 

if (present(var3))then
    if (var3(i,0).eq.undefvar)then
        heatmap_var3(i,:,:,:) = undefvar 
        cycle
    end if 
end if 

do p=1,pver

    if (logvar1) then
        if (var1(i,p)>0._r8)then
            var1tmp = log(var1(i,p))
        end if
    else
        var1tmp = var1(i,p)
    end if

    if (logvar2) then
        if (var2(i,p)>0._r8)then
            var2tmp = log(var2(i,p))
        end if
    else
        var2tmp = var2(i,p)
    end if

    ! check var3 to ensure no nan and inf
    if (present(var3))then
        call check_nan_inf(var3(i,p),'pdf2dp_regime')
    end if

    if (var1(i,p) .le. var1min .and. (var1(i,p).gt.limiter)) then
        if (var2(i,p) .le. var2min .and. (var2(i,p).gt.limiter)) then
                ! ===========================================
                idx = 1
                jj = 1
                kk = 1
                pdf_var(i,idx,jj,kk) = pdf_var(i,idx,jj,kk) + 1._r8
                if (present(var3)) then
                    heatmap_var3(i,idx,jj,kk) = heatmap_var3(i,idx,jj,kk) + var3(i,p)
                    freq_var3(i,idx,jj,kk) = freq_var3(i,idx,jj,kk) + 1.0_r8
                end if

                if (LTS(i) .le. LTS_threshold) then
                    if (RH750(i) .le. RH750_threshold) then
                        idx = 2
                    else
                        idx = 3
                    end if
                else
                    if (RH750(i) .le. RH750_threshold) then
                        idx = 4
                    else
                        idx = 5
                    end if
                end if

                pdf_var(i,idx,jj,kk) = pdf_var(i,idx,jj,kk) + 1._r8
                if (present(var3)) then
                    heatmap_var3(i,idx,jj,kk) = heatmap_var3(i,idx,jj,kk) + var3(i,p)
                    freq_var3(i,idx,jj,kk) = freq_var3(i,idx,jj,kk) + 1.0_r8
                end if

    end if 
    end if 

    if (var1(i,p) .le. var1min .and. (var1(i,p).gt.limiter)) then
        do k=1,nbins2
            if (var2tmp.gt.bvars2_1d(k) .and. var2tmp.le.bvars2_1d(k+1)) then
                idx = 1
                jj = 1
                kk = k
                pdf_var(i,idx,jj,kk) = pdf_var(i,idx,jj,kk) + 1._r8
                if (present(var3)) then
                    heatmap_var3(i,idx,jj,kk) = heatmap_var3(i,idx,jj,kk) + var3(i,p)
                    freq_var3(i,idx,jj,kk) = freq_var3(i,idx,jj,kk) + 1.0_r8
                end if

                if (LTS(i) .le. LTS_threshold) then
                    if (RH750(i) .le. RH750_threshold) then
                        idx = 2
                    else
                        idx = 3
                    end if
                else
                    if (RH750(i) .le. RH750_threshold) then
                        idx = 4
                    else
                        idx = 5
                    end if
                end if

                pdf_var(i,idx,jj,kk) = pdf_var(i,idx,jj,kk) + 1._r8
                if (present(var3)) then
                    heatmap_var3(i,idx,jj,kk) = heatmap_var3(i,idx,jj,kk) + var3(i,p)
                    freq_var3(i,idx,jj,kk) = freq_var3(i,idx,jj,kk) + 1.0_r8
                end if
                ! ===========================================
            end if
        end do ! k
    end if

    if (var2(i,p) .le. var2min .and. (var2(i,p).gt.limiter)) then
        do j=1,nbins1
            if (var1tmp.gt.bvars1_1d(j) .and. var1tmp.le.bvars1_1d(j+1)) then
                ! ===========================================
                idx = 1
                jj = j
                kk = 1
                pdf_var(i,idx,jj,kk) = pdf_var(i,idx,jj,kk) + 1._r8
                if (present(var3)) then
                    heatmap_var3(i,idx,jj,kk) = heatmap_var3(i,idx,jj,kk) + var3(i,p)
                    freq_var3(i,idx,jj,kk) = freq_var3(i,idx,jj,kk) + 1.0_r8
                end if

                if (LTS(i) .le. LTS_threshold) then
                    if (RH750(i) .le. RH750_threshold) then
                        idx = 2
                    else
                        idx = 3
                    end if
                else
                    if (RH750(i) .le. RH750_threshold) then
                        idx = 4
                    else
                        idx = 5
                    end if
                end if

                pdf_var(i,idx,jj,kk) = pdf_var(i,idx,jj,kk) + 1._r8
                if (present(var3)) then
                    heatmap_var3(i,idx,jj,kk) = heatmap_var3(i,idx,jj,kk) + var3(i,p)
                    freq_var3(i,idx,jj,kk) = freq_var3(i,idx,jj,kk) + 1.0_r8
                end if
                ! ===========================================


            end if
        end do ! j
    end if 

    do j=1,nbins1
    do k=1,nbins2
        if (var1tmp.gt.bvars1_1d(j) .and. var1tmp.le.bvars1_1d(j+1)) then
            if (var2tmp.gt.bvars2_1d(k) .and. var2tmp.le.bvars2_1d(k+1)) then
                ! ===========================================
                idx = 1
                jj = j
                kk = k
                pdf_var(i,idx,jj,kk) = pdf_var(i,idx,jj,kk) + 1._r8
                if (present(var3)) then
                    heatmap_var3(i,idx,jj,kk) = heatmap_var3(i,idx,jj,kk) + var3(i,p)
                    freq_var3(i,idx,jj,kk) = freq_var3(i,idx,jj,kk) + 1.0_r8
                end if

                if (LTS(i) .le. LTS_threshold) then
                    if (RH750(i) .le. RH750_threshold) then
                        idx = 2
                    else
                        idx = 3
                    end if
                else
                    if (RH750(i) .le. RH750_threshold) then
                        idx = 4
                    else
                        idx = 5
                    end if
                end if

                pdf_var(i,idx,jj,kk) = pdf_var(i,idx,jj,kk) + 1._r8
                if (present(var3)) then
                    heatmap_var3(i,idx,jj,kk) = heatmap_var3(i,idx,jj,kk) + var3(i,p)
                    freq_var3(i,idx,jj,kk) = freq_var3(i,idx,jj,kk) + 1.0_r8
                end if
                ! ===========================================
            end if
        end if
    end do ! k
    end do ! j

end do ! p
end do ! i

! get the value mean in specific bin
if (present(var3)) then
    do i=1,pcols
        do ireg=1,N_REGIME
        do j=1,nbins1
        do k=1,nbins2
            if (freq_var3(i,ireg,j,k) .ne. 0._r8) then
                heatmap_var3(i,ireg,j,k) = heatmap_var3(i,ireg,j,k)/freq_var3(i,ireg,j,k)
            else
                heatmap_var3(i,ireg,j,k) = undefvar
            end if
        end do ! k
        end do ! j
        end do ! ireg
    end do ! i
end if

end subroutine pdf2dp_regime

! ---------------------------------------------------------------

subroutine check_nan_inf(var,subname)
    real(r8) var
    character(len=*) :: subname 
    if (isnan(var)) then
       write(*,*) var
       call endrun(subname // ':: subroutine failed (NaN)')
    else if (isinf(var)) then
       call endrun(subname // ':: subroutine failed (inf)')
    else if (var>huge(1._r8) .or. var<-huge(1._r8)) then
       call endrun(subname // ':: subroutine failed (overflow)')
    end if
end subroutine check_nan_inf
! ---------------------------------------------------------------


end module cdiag_pdf
