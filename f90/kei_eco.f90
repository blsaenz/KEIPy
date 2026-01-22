module kei_eco
  !-----------------------------------------------------------------------------
  !   CVS:$Id: ecosys_parms.F90,v 1.1 2004/05/12 18:27:14 ivan Exp $
  !   CVS:$Name:  $
  !-----------------------------------------------------------------------------
  !   Modified to include parameters for diazotrophs, JKM  4/2002
  !-----------------------------------------------------------------------------
  !   variables/subroutines/function used from other modules
  !   The following are used extensively in this ecosys, so are used at
  !   the module level. The use statements for variables that are only needed
  !   locally are located at the module subprogram level.
  !-----------------------------------------------------------------------------
  use kei_ecocommon

  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  !   public/private declarations
  !   The default is private - public parameters should be in kei_ecocommon module
  !-----------------------------------------------------------------------------
  Public :: ecosys_step, ecosys_init, ecosys_set_interior, ecosys_set_sflux

  Private
  Save

  !-----------------------------------------------------------------------------
  !   floating point constants used across ecosystem module
  !-----------------------------------------------------------------------------

  REAL(KIND=dbl_kind), PARAMETER :: &
       c1 = 1.0_dbl_kind,             &
       spd = 86400.0_dbl_kind,        & ! number of seconds in a day
       dps = c1 / spd,                & ! number of days in a second
       yps = c1 / (365_dbl_kind*spd), & ! number of years in a second
       c0 = 0.0_dbl_kind,             &
       c2 = 2.0_dbl_kind,              &
       c3 = 3.0_dbl_kind,              &
       c10 = 10.0_dbl_kind,            &
       p5 = 0.5_dbl_kind,              &
       c100 = 100.0_dbl_kind,        &
       c1000 = 1000.0_dbl_kind,        &
       rho_sw = 4.1_dbl_kind/3.996_dbl_kind,      & ! density of salt water (g/cm^3)
       T0_Kelvin = 273.15_dbl_kind
  !-----------------------------------------------------------------------------
  !   Redfield Ratios, dissolved & particulate
  !-----------------------------------------------------------------------------

  REAL(KIND=dbl_kind), PARAMETER :: &
       parm_Red_D_C_P  = 117.0_dbl_kind,                 & ! carbon:phosphorus
       parm_Red_D_N_P  =  16.0_dbl_kind,                 & ! nitrogen:phosphorus
       parm_Red_D_O2_P = 170.0_dbl_kind,                 & ! oxygen:phosphorus
       parm_Red_P_C_P  = parm_Red_D_C_P,                 & ! carbon:phosphorus
       parm_Red_D_C_N  = parm_Red_D_C_P/parm_Red_D_N_P,  & ! carbon:nitrogen
       parm_Red_P_C_N  = parm_Red_D_C_N,                 & ! carbon:nitrogen
       parm_Red_D_C_O2 = parm_Red_D_C_P/parm_Red_D_O2_P, & ! carbon:oxygen
       parm_Red_P_C_O2 = parm_Red_D_C_O2,                & ! carbon:oxygen
       parm_Red_Fe_C   = 3.0e-6_dbl_kind,                & ! iron:carbon
       parm_Red_D_C_O2_diaz = parm_Red_D_C_P/150.0_dbl_kind! carbon:oxygen
                                                           ! for diazotrophs

  !----------------------------------------------------------------------------
  !   ecosystem parameters accessible via input file
  !----------------------------------------------------------------------------

  REAL(KIND=dbl_kind) :: &
       parm_Fe_bioavail,      & ! fraction of Fe flux that is bioavailable
       parm_prod_dissolve,    & ! frac. of prod -> DOC
       parm_o2_min,           & ! min O2 needed for prod & consump. (nmol/cm^3)
       parm_no3_min,          & ! min NO3 needed for denitrification (mmol/m^3)
       parm_Rain_CaCO3,       & ! Rain ratio for CaCO3
       parm_Rain_SiO2,        & ! Rain ratio for SiO2
       parm_kappa_nitrif,     & ! nitrification inverse time constant (1/sec)
       parm_nitrif_par_lim,   & ! PAR limit for nitrif. (W/m^2)
       parm_POC_flux_ref,     & ! reference POC flux (nmol C/cm^2/sec)
       parm_rest_prod_tau,    & ! time-scale for restoring prod (sec)
       parm_rest_prod_z_c,    & ! depth-limit for restoring (cm)
       parm_z_umax_0,         & ! max. zoo growth rate on sphyto at tref (1/sec)
       parm_diat_umax_0,      & ! max. zoo growth rate on diat at tref (1/sec)
       parm_z_mort_0,         & ! zoo linear mort rate (1/sec)
       parm_z_mort2_0,        & ! zoo quad mort rate (1/sec/((mmol C/m3))
       parm_sd_remin_0,       & ! small detrital remineralization rate (1/sec)
       parm_sp_kNO3,          & ! sphyto nitrate half sat. coef. (mmol N/m3)
       parm_diat_kNO3,        & ! diatom nitrate half sat. coef. (mmol N/m3)
       parm_sp_kNH4,          & ! sphyto ammonium half sat. coef. (mmol N/m3)
       parm_diat_kNH4,        & ! diatom ammonium half sat. coef. (mmol N/m3)
       parm_sp_kFe,           & ! sphyto iron half sat. coef. (nmol Fe/m3)
       parm_diat_kFe,         & ! diatom iron half sat. coef. (nmol Fe/m3)
       parm_diat_kSiO3,       & ! diatom si half sat. coef. (mmol SiO3/m3)
       parm_sp_kPO4,          & ! sphyto PO4 uptake (mmol P/m^3)
       parm_diat_kPO4,        & ! diatom PO4 uptate (mmol P/m^3)
       parm_z_grz,            & ! grazing coef. for small phyto (mmol C/m^3)
       parm_alphaChl,         & ! init slope of P_I curve (GD98) (mmol C m^2/(mg Chl W sec))
       parm_alphaChlsp,         & ! init slope of P_I curve (GD98) (mmol C m^2/(mg Chl W sec))
       parm_alphaChldiat,         & ! init slope of P_I curve (GD98) (mmol C m^2/(mg Chl W sec))
       parm_alphaChlphaeo,         & ! init slope of P_I curve (GD98) (mmol C m^2/(mg Chl W sec))
       parm_labile_ratio,     & ! fraction of loss to DOC that routed directly to DIC (non-dimensional)
       parm_alphaDiaz,        & ! chl. spec. init. slope of diaz. P_I curve
       parm_diaz_umax_0         ! max. zoo growth rate on diazotrophs at tref (1/sec)

  !----------------------------------------------------------------------------
  !   ecosystem parameters not (yet?) accessible via input file
  !----------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  !     Parameters previously hardcoded in main body of bio_subs code
  !---------------------------------------------------------------------------

  real(kind=dbl_kind), parameter ::      &
      PCref = 3.0_dbl_kind * dps, & !max phyto C-spec. grth rate at tref (1/sec)
      PCrefSp = 4.5_dbl_kind * dps, & !max phyto C-spec. grth rate at tref (1/sec)
      PCrefDiat = 4.5_dbl_kind * dps, & !max phyto C-spec. grth rate at tref (1/sec)
      sp_mort    = 0.1_dbl_kind   * dps, & !sphyto mort rate (1/sec)
!      sp_mort2   = 0.009_dbl_kind * dps, & !sphyto quad. mort rate (1/sec/((mmol C/m3))
      diat_mort  = 0.1_dbl_kind   * dps, & !diatom mort rate (1/sec)
!      diat_mort2 = 0.009_dbl_kind * dps, & !diatom quad mort rate (1/sec/((mmol C/m3))

      sp_mort2   = 0.0009_dbl_kind * dps, & !sphyto quad. mort rate (1/sec/((mmol C/m3))
      diat_mort2   = 0.0009_dbl_kind * dps, & !sphyto quad. mort rate (1/sec/((mmol C/m3))


!      sp_mort    = 0.17_dbl_kind   * dps, & !sphyto mort rate (1/sec)
!      sp_mort2   = 0.0035_dbl_kind * dps, & !sphyto quad. mort rate (1/sec/((mmol C/m3))
!      diat_mort  = 0.17_dbl_kind   * dps, & !diatom mort rate (1/sec)
!      diat_mort2 = 0.0035_dbl_kind * dps, & !diatom quad mort rate (1/sec/((mmol C/m3))
      PCrefDiaz  = 0.4_dbl_kind  * dps,  & !max Diaz C-specific growth rate at tref (1/sec)
      diaz_mort  = 0.16_dbl_kind * dps,  & !diaz mort rate (1/sec)
      diaz_kPO4  = 0.005_dbl_kind,      & !diaz half-sat. const. for P (diatom value)
      diaz_kFe   = 0.1e-3_dbl_kind        !diaz half-sat. const. for Fe

  !---------------------------------------------------------------------
  !     Misc. Rate constants
  !---------------------------------------------------------------------
  real(kind=dbl_kind), parameter :: &
!       sp_agg_rate_max   = 0.2_dbl_kind, & !max agg. rate for small phyto (1/d)
!       diat_agg_rate_max = 0.2_dbl_kind, & !max agg. rate for diatoms (1/d)
       sp_agg_rate_max   = 0.75_dbl_kind, & !max agg. rate for small phyto (1/d)
       diat_agg_rate_max = 0.75_dbl_kind, & !max agg. rate for diatoms (1/d)
       diat_agg_rate_min = 0.01_dbl_kind,& !min agg. rate for diatoms (1/d)
       fe_scavenge_rate0 = 0.12_dbl_kind,& !init Fe scaveng. rate (% of ambient)
       fe_scavenge_thres1 = 0.6e-3_dbl_kind, & !upper thres. for Fe scavenging
       fe_scavenge_thres2 = 0.5e-3_dbl_kind, & !lower thres. for Fe scavenging
       dust_fescav_scale  = 0.833e8_dbl_kind, & !dust scavenging scale factor
       thres_fe           = 1.0e5_dbl_kind,   & !thres. depth for Fe diff. flux
       fe_max_scale1      = 3.0_dbl_kind,     & !unitless scaling coeff.
       fe_max_scale2      = 6.0_dbl_kind/1.4e-3_dbl_kind,& !unitless scaling coeff.
       fe_diff_rate       = 2.3148e-6_dbl_kind,&!fe diffusion rate
                                                !   (nmolFe/cm2/sec)
       f_fescav_P_iron    = 0.1_dbl_kind        !fraction of Fe scavenging
                                                !        to particulate Fe

  !---------------------------------------------------------------------
  !     Compute iron remineralization and flux out.
  !     dust remin gDust = 0.035 gFe      mol Fe     1e9 nmolFe
  !                        --------- *  ---------- * ----------
  !          gDust       55.847 gFe     molFe
  !
  !     dust_to_Fe          conversion - dust to iron (nmol Fe/g Dust)
  !---------------------------------------------------------------------
  real(kind=dbl_kind), parameter :: &
       dust_to_Fe=0.035_dbl_kind/55.847_dbl_kind*1.0e9_dbl_kind

  !----------------------------------------------------------------------------
  !     Partitioning of phytoplankton growth, grazing and losses
  !
  !     All f_* variables are fractions and are non-dimensional
  !----------------------------------------------------------------------------

  real(kind=dbl_kind), parameter ::     &
      z_ingest         = 0.15_dbl_kind,  & !zoo ingestion coefficient (non-dim)
      caco3_poc_min    = 0.4_dbl_kind,  & !minimum proportionality between
                                          !   QCaCO3 and grazing losses to POC
                                          !   (mmol C/mmol CaCO3)
      spc_poc_fac      = 0.22_dbl_kind, & !small phyto grazing factor (1/mmolC)
      f_graze_sp_poc_lim = 0.24_dbl_kind, &
      f_prod_sp_CaCO3  = 0.026_dbl_kind, & !fraction of sp prod. as CaCO3 prod.
      f_photosp_CaCO3  = 0.4_dbl_kind,  & !proportionality between small phyto
                                          !    production and CaCO3 production
      f_graze_sp_doc   = 0.34_dbl_kind, & !fraction sm. phyto. grazing to DOC
      f_graze_sp_dic   = c1 - z_ingest - f_graze_sp_doc, & !fraction to DIC
      f_z_grz_sqr_diat = 0.81_dbl_kind, &   ! original
!      f_z_grz_sqr_diat = 3.5_dbl_kind, &     ! sevrine recommendation
      f_graze_diat_poc = 0.26_dbl_kind, & !fraction diatom grazing to POC
      f_graze_diat_doc = 0.13_dbl_kind, & !fraction diatom grazing to DOC
      f_graze_diat_dic = c1 - z_ingest - f_graze_diat_poc &
                          - f_graze_diat_doc, & !fraction diatom grazing to DIC
      f_diat_loss_poc  = 0.05_dbl_kind, &  !fraction diatom loss to POC
      f_diat_loss_dc   = c1-f_diat_loss_poc, & !fraction diatom loss to DOC
      f_graze_diaz_zoo = 0.21_dbl_kind, & !fraction diaz. grazing to zoo
      f_graze_diaz_poc = 0.0_dbl_kind, &  !fraction diaz grazing to POC
      f_graze_diaz_doc = 0.24_dbl_kind, & !fraction diaz grazing to DOC
      f_graze_diaz_dic = c1-f_graze_diaz_zoo-f_graze_diaz_poc &
                         - f_graze_diaz_doc, & !fraction diaz grazing to DIC
      f_sp_zoo_detr   = 0.06666_dbl_kind,& !fraction of zoo losses to detrital
                                          !  pool when eating sphyto
      f_diat_zoo_detr = 0.1333_dbl_kind,& !fraction of zoo losses to detrital
                                          !  pool when eating diatoms
      f_diaz_zoo_detr = 0.03333_dbl_kind,& !fraction of zoo losses to detrital
                                          !  pool when eating diaz
      f_graze_CaCO3_remin = 0.33_dbl_kind, & !fraction of spCaCO3 grazing
                                             !          which is remin
      f_graze_si_remin    = 0.5_dbl_kind      !fraction of diatom Si grazing
                                             !          which is remin

  !----------------------------------------------------------------------------
  !     fixed ratios
  !----------------------------------------------------------------------------
  real(kind=dbl_kind), parameter :: &
       r_Nfix_photo=1.43_dbl_kind         ! N fix relative to C fix (non-dim)

  !-----------------------------------------------------------------------
  !     SET FIXED RATIOS for N/C, P/C, SiO3/C, Fe/C
  !     assumes C/N/P of 117/16/1 based on Anderson and Sarmiento, 1994
  !     for diazotrophs a N/P of 45 is assumed based on Letelier & Karl, 1998
  !-----------------------------------------------------------------------

  real(kind=dbl_kind), parameter ::  &
      Q             = 0.137_dbl_kind,  & !N/C ratio (mmol/mmol) of phyto & zoo
      Qp            = 0.00855_dbl_kind,& !P/C ratio (mmol/mmol) sphyto,diat,zoo
      Qp_diaz       = 0.002735_dbl_kind,& !diazotroph P/C ratio
      Qfe_zoo       = 2.5e-6_dbl_kind, & !zooplankton fe/C ratio
      gQsi_0        = 0.137_dbl_kind,  & !initial diatom Si/C ratio
      gQfe_diat_0   = 6.0e-6_dbl_kind, & !initial diatom fe/C ratio
      gQfe_sp_0     = 6.0e-6_dbl_kind, & !initial sphyto fe/C ratio
      gQfe_diaz_0   = 42.0e-6_dbl_kind,& !initial diaz. fe/C ratio
      gQfe_diat_min = 2.5e-6_dbl_kind, & !min diatom fe/C ratio
      gQsi_max      = 0.685_dbl_kind,  & !max diatom Si/C ratio
      gQsi_min      = 0.0685_dbl_kind, & !min diatom Si/C ratio
      gQsi_coef     = 2.5_dbl_kind, &
      gQfe_sp_min   = 2.5e-6_dbl_kind, & !min sphyto fe/C ratio
      gQfe_diaz_min = 14.0e-6_dbl_kind,& !min diaz fe/C ratio
      QCaCO3_max    = 0.4_dbl_kind,    & !max QCaCO3
      thetaN_max_sp   = 2.5_dbl_kind, & !sp max thetaN (Chl/N) (mg Chl/mmol N)
      thetaN_max_diat = 4.0_dbl_kind, & !diat max thetaN (Chl/N) (mg Chl/mmol N)
!      thetaN_max_sp   = 2.3_dbl_kind, & !sp max thetaN (Chl/N) (mg Chl/mmol N)
!      thetaN_max_diat = 3.0_dbl_kind, & !diat max thetaN (Chl/N) (mg Chl/mmol N)
      thetaN_max_diaz = 3.4_dbl_kind, & !diaz max thetaN (Chl/N) (mg Chl/mmol N)
      ! carbon:nitrogen ratio for denitrification
      ! net removal of 120 mols NO3 for 117 mols C (136 = 120 + 16)
      denitrif_C_N  = parm_Red_D_C_P/136.0_dbl_kind

  !----------------------------------------------------------------------------
  !     loss term threshold parameters, chl:c ratios
  !----------------------------------------------------------------------------

  real(kind=dbl_kind), parameter ::    &
      thres_z1          = 100.0e2_dbl_kind, & !threshold = C_loss_thres for z shallower than this (cm)
      thres_z2          = 200.0e2_dbl_kind, & !threshold = 0 for z deeper than this (cm)
      loss_thres_sp     = 0.003_dbl_kind, & !small phyto conc. where losses go to zero
!      loss_thres_sp     = 0.01_dbl_kind, & !small phyto conc. where losses go to zero - ben increased
      loss_thres_diat   = 0.03_dbl_kind, & !diat conc. where losses go to zero
!      loss_thres_diat   = 0.1_dbl_kind, & !diat conc. where losses go to zero - ben increased
      loss_thres_zoo    = 0.03_dbl_kind,  & !zoo conc. where losses go to zero
      loss_thres_diaz   = 0.01_dbl_kind,  & !diaz conc. where losses go to zero
      loss_thres_diaz2  = 0.001_dbl_kind, & !diaz conc. thres at low temp
      diaz_temp_thres   = 15.0_dbl_kind,  & !Temp. where diaz conc thres drops
      CaCO3_temp_thres1 = 1.0_dbl_kind,   & !upper temp threshold for CaCO3 prod
      CaCO3_temp_thres2 = -2.0_dbl_kind,  & !lower temp threshold
      CaCO3_sp_thres    = 3.0_dbl_kind      ! bloom condition thres (mmolC/m3)

  !---------------------------------------------------------------------
  !     attenuation coefficients for PAR and related parameters
  !---------------------------------------------------------------------
  real(kind=dbl_kind), parameter :: &
       k_chl = 0.03e-2_dbl_kind, & ! Chl atten. coeff. (1/cm/(mg Chl/m^3))
       k_h2o = 0.04e-2_dbl_kind, & ! water atten. coeff (1/cm)
       f_qsw_par = 0.45_dbl_kind   ! PAR fraction

  !---------------------------------------------------------------------
  !     Temperature parameters
  !---------------------------------------------------------------------
  real(kind=dbl_kind), parameter :: &
       Tref = 30.0_dbl_kind, & ! reference temperature (C)
       Q_10 = 2.0_dbl_kind     ! factor for temperature dependence (non-dim)

!-----------------------------------------------------------------------
!     derived type for implicit handling of sinking particulate matter
!-----------------------------------------------------------------------

      type sinking_particle
        real(kind=dbl_kind) :: &
            diss,        & ! dissolution length for soft subclass
            gamma,       & ! fraction of production -> hard subclass
            mass,        & ! mass of 1e9 base units in g
            rho            ! QA mass ratio of POC to this particle class
        real(kind=dbl_kind), dimension(imt,jmt) :: &
            sflux_in,    & ! incoming flux of soft subclass (base units/cm^2/sec)
            hflux_in,    & ! incoming flux of hard subclass (base units/cm^2/sec)
            prod,        & ! production term (base units/cm^3/sec)
            sflux_out,   & ! outgoing flux of soft subclass (base units/cm^2/sec)
            hflux_out,   & ! outgoing flux of hard subclass (base units/cm^2/sec)
            remin          ! remineralization term (base units/cm^3/sec)
        end type sinking_particle

      logical, dimension(imt,jmt) :: LAND_MASK = .true.

      real(kind=dbl_kind), dimension(km) :: &
        nutr_rest_time_inv ! inverse restoring time scale for nutrients (1/secs)

  !-----------------------------------------------------------------------------
  !   The current setting of xacc, a tolerance critera, will result in co2star
  !   being accurate to 3 significant figures (xx.y). Making xacc bigger will
  !   result in faster convergence also, but this is not recommended (xacc of
  !   10**-9 drops precision to 2 significant figures).
  !-----------------------------------------------------------------------------

  REAL(KIND=dbl_kind), PARAMETER :: xacc = 1e-10_dbl_kind
  INTEGER(KIND=int_kind), PARAMETER :: maxit = 100

  !-----------------------------------------------------------------------------
  !   declarations for function coefficients & species concentrations
  !-----------------------------------------------------------------------------

  REAL(KIND=dbl_kind), DIMENSION(imt) :: &
       k0, k1, k2, kw, kb, ks, kf, k1p, k2p, k3p, ksi, ff, &
       bt, st, ft, dic, ta, pt, sit

  REAL(KIND=dbl_kind), DIMENSION(imt,jmt) :: &
       PH_PREV

  !-----------------------------------------------------------------------------
  !   holding variables to account for different tracer precision
  !   in kpp and ecosys_mod
  !-----------------------------------------------------------------------------
  REAL(KIND=dbl_kind) :: &
       dz_eco(km), &
       dzr_eco(km), &
       zt_eco(km)

  INTEGER(KIND=int_kind), PARAMETER :: wavb = 31  ! PAR/PUR wavelength bins
  ! default surface spectrum, if only using single PAR band in physical model
  real(KIND=dbl_kind), PARAMETER, DIMENSION(wavb) :: &
     srfspec = (/ 2.3421E-02, 2.3241E-02, 2.4949E-02, 2.3446E-02, &
      2.8255E-02, 3.2075E-02, 3.3944E-02, 3.4012E-02, 3.5491E-02, &
      3.3682E-02, 3.4819E-02, 3.5322E-02, 3.3823E-02, 3.5034E-02, &
      3.5200E-02, 3.5458E-02, 3.5001E-02, 3.4680E-02, 3.4925E-02, &
      3.3849E-02, 3.4200E-02, 3.4444E-02, 3.3892E-02, 3.3023E-02, &
      3.3942E-02, 3.3048E-02, 3.1756E-02, 3.3287E-02, 3.2811E-02, &
      2.8603E-02, 3.0369E-02 /)

CONTAINS

!***********************************************************************

      subroutine ecosys_step(kforce,X,fice,albocn,Sref,dt_sec,qsw_ice,hmix,nt,start_year, &
                             par_phyto,absorp_in,absorp_out)

        USE kei_hacks

        implicit none

        ! subroutine arguments
        ! --------------------------------------------------------------
        real :: kforce(forcing_var_cnt)
        real(KIND=real_kind) :: X(NZP1,NSCLR)    ! tracers
        real(KIND=real_kind) :: fice ! fractional ice coverage
        real(KIND=real_kind) :: albocn ! ocean albedo - surf irradiance stuff should be outside/in one place!
        real, dimension(NZ) :: &
          qsw                 ! penetrative (absorbed?) solar heat flux (W/m^2)
        real(KIND=real_kind) :: Sref ! reference salinity to derive true salinity
         real(KIND=dbl_kind) :: dt_sec ! timestep in seconds
         real(KIND=real_kind) :: qsw_ice ! shortwave irradiance transmitted through ice
         real(KIND=real_kind) :: hmix ! mixed layer (+m)
         integer(KIND=int_kind) :: nt,start_year


        ! local variables
        ! --------------------------------------------------------------
        real(KIND=dbl_kind), dimension(imt,jmt,ecosys_tracer_cnt) :: &
          TRACER_MODULE, &    !  tracers
          DTRACER_MODULE      !  tracer changes
        real(KIND=dbl_kind), dimension(imt,jmt) :: &
             fice_eco, &
             ap_eco, &
             winds_SQR_eco, &
             t_eco, &
             s_eco, &
             qsw_eco
        integer(KIND=int_kind) :: i,j,k,tr_end
        real(kind=dbl_kind), dimension (ecosys_tracer_cnt) :: &
          eco_inject
        real(KIND=dbl_kind), intent(out) :: par_phyto(km) ! shortwave irradiance transmitted
        real(KIND=dbl_kind), intent(in) :: absorp_in(km) ! external sources shortwave irradiance absorption (1/m)
        real(KIND=dbl_kind), intent(out) :: absorp_out(km) ! total ecosystem shortwave irradiance absorption (1/m)

         real(KIND=dbl_kind) :: PAR_out,PAR_in,PAR_avg,dz_ml, KPARdz, dz1
         integer(KIND=int_kind) :: k_ml


!-----------------------------------------------------------------------
!     equate forcings in kpp model and ecosys module
!-----------------------------------------------------------------------

        qsw_eco = kforce(qswins_f_ind)*(c1-fice)*(c1-albocn) + &
          qsw_ice*fice

         dust_flux = kforce(dustf_f_ind) ! surface dust flux
!       iron_flux,         ! iron component of surface dust flux
         winds_SQR = (kforce(taux_f_ind)*100)**2 + &
           (kforce(tauy_f_ind)*100)**2
         atm_co2 =  atm_co2_const  ! atmospheric CO2 concentration (
         ap = ap_const*1013.25e+3_dbl_kind      ! dyne / cm^2

!-----------------------------------------------------------------------
!     hacks
!-----------------------------------------------------------------------

        CALL do_eco_hacks(X,nt,start_year)

!-----------------------------------------------------------------------
!     surface fluxes of dust/iron, O2, CO2
!-----------------------------------------------------------------------

        tr_end = ecosys_tracer_cnt+2
        do i=1,imt
          do j=1,jmt
            ! upgrade precision/shape of tracers for sflux calculations
            TRACER_MODULE(i,j,:) = X(1,3:tr_end)
            fice_eco(i,j) = fice
            ap_eco(i,j) = ap
            winds_SQR_eco(i,j) = winds_SQR
            t_eco(i,j) = X(1,1)
            s_eco(i,j) = X(1,2) + Sref

            ! call sflux
            call ecosys_set_sflux( &
              fice_eco,        & ! IFRAC = sea ice fraction (non-dimensional)
              ap_eco,          & ! ATM_PRESS = sea level atmospheric pressure (dyne/cm^2)
              winds_SQR_eco,  &  ! U10_SQR = 10m wind speed squared (cm^2/s^2)
              t_eco,          & ! SST = sea surface temperature (C)
              s_eco,          & ! SSS = sea surface salinity (psu)
              TRACER_MODULE,  & ! SURF_VALS = surface ecosys tracers
              DTRACER_MODULE ) ! STF_MODULE = computed source/sink ecosys tracers (output)

          enddo
        enddo


           ! update surface tracers, just once since kpp is not 3D
          TRACER_MODULE(1,1,:) = TRACER_MODULE(1,1,:) + &
            DTRACER_MODULE(1,1,:)*dt_sec
          !TRACER_MODULE(1,1,:) = max(c0,TRACER_MODULE(1,1,:))
                ! O2 flux can go too fast, and oscillate out of range - limit to saturation
        if (TRACER_MODULE(1,1,o2_ind) > O2SAT_tavg(1,1)) then
                    TRACER_MODULE(1,1,o2_ind) = O2SAT_tavg(1,1)
                endif
           X(1,3:tr_end) = TRACER_MODULE(1,1,:)

!-----------------------------------------------------------------------
!     find k-level light that phytoplankton see, averaged over mixed layer
!-----------------------------------------------------------------------
        PAR_out = max(c0, f_qsw_par * qsw_eco(1,1))
        par_phyto = c0
        PAR_avg = c0
        dz_ml = c0
        k_ml = 0
        do k=1,km
          PAR_in = PAR_out
          dz1 = dz_eco(k)/c100

          absorp_out(k) = absorp_in(k)/100._dbl_kind  ! 1/m -> 1/cm
          absorp_out(k) = absorp_out(k) &
                    + k_chl * (max(c0, X(k,2+diazChl_ind)) &
                    + max(c0, X(k,2+diatChl_ind)) &
                    + max(c0, X(k,2+spChl_ind))) &
                    + k_h2o

          if (PAR_in > 1.0e-6) then
            !KPARdz = (k_chl * (max(c0, X(k,2+diazChl_ind)) &
            !        + max(c0, X(k,2+diatChl_ind)) &
            !        + max(c0, X(k,2+spChl_ind))) &
            !        + k_h2o) * dz_eco(k)

            KPARdz =  absorp_out(k) * dz_eco(k)
            PAR_out = PAR_in * exp(-KPARdz)
            par_phyto(k) = (PAR_out + PAR_in)/c2
          else
            par_phyto(k) = c0
            PAR_out = c0
          endif
          if (dz_ml <= hmix) then
              PAR_avg = PAR_avg + par_phyto(k)*dz1
              dz_ml = dz_ml + dz1
              k_ml = k
          endif

          absorp_out(k) = absorp_in(k)/100._dbl_kind ! convert back to 1/m for output

        enddo
        if (k_ml > 1) then
          par_phyto(1:k_ml) = PAR_avg / dz_ml
        endif
        !print *,'hmix:', hmix, 'k_ml',k_ml, 'dz_m',dz_ml
        !print *,'par_phyto:', par_phyto


!-----------------------------------------------------------------------
!     ecosystem step
!-----------------------------------------------------------------------
        do k=1,km
          ! upgrade precision/shape of tracers for ecosys calculations
          if (k.eq.1) then
            eco_inject  = ice_to_ocean_eflux

            if (nt > 8664 .and. nt < 8792) then
              eco_inject(Fe_ind) = eco_inject(Fe_ind) + meltwater_fe
            endif

          else
            eco_inject = c0
          endif


          do i=1,imt
            do j=1,jmt
              TRACER_MODULE(i,j,:) = X(k,3:tr_end)
              t_eco(i,j) = X(k,1)

              ! call ecosys
              call ecosys_set_interior( &
                k,                  & ! k = vertical level index
                t_eco,              & ! TEMP = potential temperature (C)
                qsw_eco,             & ! SHF_QSW_Wpm2 = penetrative solar heat flux (W/m^2)
                TRACER_MODULE,      & ! TRACER_MODULE = current tracer values
                DTRACER_MODULE,     & ! DTRACER_MODULE = computed source/sink tracers (output)
                NZ,                 & ! KMT = lowest valid layer (<=km)
                dz_eco,              &  ! dz = layer thickness
                dzr_eco,            & ! dzr = reciprocal of dz
                zt_eco,             & ! zt = vert dist from sfc to midpoint of layer
                eco_inject , par_phyto(k))

                ! write out ecosys data
                tot_prod(k) = tot_prod_tavg(i,j)
                diat_Fe_lim(k) = diat_Fe_lim_tavg(i,j)
                diat_N_lim(k) = diat_N_lim_tavg(i,j)
                diat_P_lim(k) = diat_PO4_lim_tavg(i,j)
                diat_Si_lim(k) = diat_SiO3_lim_tavg(i,j)
                diat_light_lim(k) = diat_light_lim_tavg(i,j)
                sp_Fe_lim(k) = sp_Fe_lim_tavg(i,j)
                sp_N_lim(k) = sp_N_lim_tavg(i,j)
                sp_P_lim(k) = sp_PO4_lim_tavg(i,j)
                sp_light_lim(k) = sp_light_lim_tavg(i,j)
                graze_diat(k) = graze_diat_tavg(i,j)
                graze_sp(k) = graze_sp_tavg(i,j)
                graze_tot(k) = graze_tot_tavg(i,j)

              enddo
          enddo



          ! update ecosys tracers
          TRACER_MODULE(1,1,:) = TRACER_MODULE(1,1,:) + &
            DTRACER_MODULE(1,1,:)*dt_sec
          !TRACER_MODULE(1,1,:) = max(c0,TRACER_MODULE(1,1,:))
          X(k,3:tr_end) = TRACER_MODULE(1,1,:)


        enddo

      end subroutine ecosys_step


!***********************************************************************

      subroutine ecosys_init(dz,zt)

!-----------------------------------------------------------------------
!     subroutine arguments
!-----------------------------------------------------------------------

      real(kind=real_kind), dimension(km) :: &
          dz,         & ! layer thickness (m)
          zt            ! layer midpoint position (negative,m)

!-----------------------------------------------------------------------
!     local variable declarations
!-----------------------------------------------------------------------

      integer (kind=int_kind) :: &
       n,                        & ! tracer index
       k,                        & ! vertical level index
       ind,                      & ! tracer index for tracer name from namelist
       cnt                         ! count of tavg vars

!-----------------------------------------------------------------------
!     setup grid variables with correct precision for ecosys
!-----------------------------------------------------------------------
      dz_eco = dz*c100
      dzr_eco = c1/dz_eco
      zt_eco = zt*c100

!-----------------------------------------------------------------------
!     initialize ecosystem parameters
!-----------------------------------------------------------------------

      parm_Fe_bioavail    = 0.02_dbl_kind
      parm_prod_dissolve  = 0.67_dbl_kind
      parm_o2_min         = 4.0_dbl_kind
      parm_no3_min        = 32.0_dbl_kind
      parm_Rain_CaCO3     = 0.07_dbl_kind
      parm_Rain_SiO2      = 0.03_dbl_kind
      parm_kappa_nitrif   = 0.06_dbl_kind * dps       ! (= 1/( days))
      parm_nitrif_par_lim = 5.0_dbl_kind
      parm_POC_flux_ref   = 2.0e-3_dbl_kind
      parm_rest_prod_tau  = 30.0_dbl_kind * spd       ! (= 30 days)
      parm_rest_prod_z_c  = 7500_dbl_kind
      parm_z_umax_0       = 2.75_dbl_kind * dps ! original
      parm_diat_umax_0    = 2.07_dbl_kind * dps ! original
!      parm_z_umax_0       = 4.0_dbl_kind * dps
!      parm_z_umax_0       = 3.1_dbl_kind * dps  ! per wang and moore - standard
!      parm_z_umax_0       = 1.5_dbl_kind * dps  ! sevrine recommendation
!      parm_diat_umax_0    = 3.06_dbl_kind * dps

!      parm_z_umax_0    = 0.2_dbl_kind * dps
!      parm_diat_umax_0    = 0.2_dbl_kind * dps

      parm_z_umax_0    = 1.5_dbl_kind * dps
      parm_diat_umax_0    = 1.5_dbl_kind * dps


!      parm_z_mort_0       = 0.1_dbl_kind * dps
!      parm_z_mort2_0      = 0.45_dbl_kind * dps
      parm_z_mort_0       = 0.17_dbl_kind * dps ! wang and moore
      parm_z_mort2_0      = 0.0035_dbl_kind * dps ! wang and moore
      parm_sd_remin_0     = 0.01_dbl_kind * dps       ! (= 1/(100 days))
      parm_sp_kNO3        = 0.5_dbl_kind
      parm_diat_kNO3      = 2.5_dbl_kind
!      parm_sp_kNH4        = 0.005_dbl_kind
!      parm_diat_kNH4      = 0.08_dbl_kind
      parm_sp_kNH4        = 0.01_dbl_kind
      parm_diat_kNH4      = 0.1_dbl_kind
!      parm_sp_kFe         = 0.06e-3_dbl_kind
!      parm_diat_kFe       = 0.15e-3_dbl_kind
      parm_sp_kFe         = 0.035e-3_dbl_kind
      parm_diat_kFe       = 0.08e-3_dbl_kind
      parm_diat_kSiO3     = 1.0_dbl_kind
!      parm_sp_kPO4        = 0.0003125_dbl_kind
!      parm_diat_kPO4      = 0.005_dbl_kind
      parm_sp_kPO4        = 0.01_dbl_kind
      parm_diat_kPO4      = 0.1_dbl_kind
      parm_z_grz          = 1.05_dbl_kind     ! original
!      parm_z_grz          = 0.80_dbl_kind
!      parm_alphaChl       = 0.3_dbl_kind * dps
      parm_alphaChl       = 0.25_dbl_kind * dps
!      parm_alphaChlsp     = 0.60_dbl_kind * dps  ! ben playing
      parm_alphaChlsp     = 0.28_dbl_kind * dps
!      parm_alphaChldiat   = 0.51_dbl_kind * dps  ! ben playing
      parm_alphaChldiat   = 0.25_dbl_kind * dps
      parm_alphaChlphaeo   = 0.68_dbl_kind * dps
      parm_labile_ratio   = 0.70_dbl_kind
      parm_alphaDiaz      = 0.036_dbl_kind * dps
      parm_diaz_umax_0    = 1.2_dbl_kind * dps

      PH_PREV = 0.0_dbl_kind

!-----------------------------------------------------------------------
!     initialize restoring timescale (if required)
!-----------------------------------------------------------------------

      if (lrest_po4 .or. lrest_no3 .or. lrest_sio3) then
        do k = 1,km
          if (zt(k) < rest_z0) then
             nutr_rest_time_inv(k) = rest_time_inv_surf
          else if (zt(k) > rest_z1) then
             nutr_rest_time_inv(k) = rest_time_inv_deep
          else if (rest_z1 == rest_z0) then
             nutr_rest_time_inv(k) = rest_time_inv_surf + p5 * &
                (rest_time_inv_deep - rest_time_inv_surf)
          else
             nutr_rest_time_inv(k) = rest_time_inv_surf +  &
               (zt(k) - rest_z0) / (rest_z1 - rest_z0) *  &
               (rest_time_inv_deep - rest_time_inv_surf)
          endif
        enddo
      endif
!-----------------------------------------------------------------------


      allocate( &
         SCHMIDT_O2_tavg(imt,jmt), &
         XKW_tavg(imt,jmt),             AP_tavg(imt,jmt), &
         PV_CO2_tavg(imt,jmt),          PV_O2_tavg(imt,jmt), &
         O2SAT_tavg(imt,jmt),           FG_O2_tavg(imt,jmt), &
         SCHMIDT_CO2_tavg(imt,jmt),     PH_tavg(imt,jmt), &
         CO2STAR_tavg(imt,jmt),         DCO2STAR_tavg(imt,jmt), &
         pCO2SURF_tavg(imt,jmt),        DpCO2_tavg(imt,jmt), &
         FG_CO2_tavg(imt,jmt),          IRON_FLUX_tavg(imt,jmt), &
         PROD_tavg(imt,jmt),            PO4_RESTORE_tavg(imt,jmt), &
         NO3_RESTORE_tavg(imt,jmt),     SiO3_RESTORE_tavg(imt,jmt), &
         PAR_avg_tavg(imt,jmt),         PO4STAR_tavg(imt,jmt), &
         POC_FLUX_IN_tavg(imt,jmt),     POC_PROD_tavg(imt,jmt), &
         POC_REMIN_tavg(imt,jmt),       CaCO3_FLUX_IN_tavg(imt,jmt), &
         CaCO3_PROD_tavg(imt,jmt),      CaCO3_REMIN_tavg(imt,jmt), &
         SiO2_FLUX_IN_tavg(imt,jmt),    SiO2_PROD_tavg(imt,jmt), &
         SiO2_REMIN_tavg(imt,jmt),      dust_FLUX_IN_tavg(imt,jmt), &
         dust_REMIN_tavg(imt,jmt),      P_iron_FLUX_IN_tavg(imt,jmt), &
         P_iron_PROD_tavg(imt,jmt),     P_iron_REMIN_tavg(imt,jmt), &
         graze_sp_tavg(imt,jmt),        graze_diat_tavg(imt,jmt), &
         graze_tot_tavg(imt,jmt),       sp_loss_tavg(imt,jmt), &
         diat_loss_tavg(imt,jmt),       zoo_loss_tavg(imt,jmt), &
         sp_agg_tavg(imt,jmt),          diat_agg_tavg(imt,jmt), &
         photoC_sp_tavg(imt,jmt),       photoC_diat_tavg(imt,jmt), &
         tot_prod_tavg(imt,jmt),        DOC_prod_tavg(imt,jmt), &
         DOC_remin_tavg(imt,jmt),       Fe_scavenge_tavg(imt,jmt), &
         sp_N_lim_tavg(imt,jmt),        sp_Fe_lim_tavg(imt,jmt), &
         sp_PO4_lim_tavg(imt,jmt),      sp_light_lim_tavg(imt,jmt), &
         diat_N_lim_tavg(imt,jmt),      diat_Fe_lim_tavg(imt,jmt), &
         diat_PO4_lim_tavg(imt,jmt),    diat_SiO3_lim_tavg(imt,jmt), &
         diat_light_lim_tavg(imt,jmt),  CaCO3_form_tavg(imt,jmt), &
         diaz_Nfix_tavg(imt,jmt),       graze_diaz_tavg(imt,jmt), &
         diaz_loss_tavg(imt,jmt),       photoC_diaz_tavg(imt,jmt), &
         diaz_P_lim_tavg(imt,jmt),      diaz_Fe_lim_tavg(imt,jmt), &
         diaz_light_lim_tavg(imt,jmt),  Fe_scavenge_rate_tavg(imt,jmt), &
         DON_prod_tavg(imt,jmt),        DON_remin_tavg(imt,jmt), &
         DOFe_prod_tavg(imt,jmt),       DOFe_remin_tavg(imt,jmt), &
         DOP_prod_tavg(imt,jmt),        DOP_remin_tavg(imt,jmt), &
         bSi_form_tavg(imt,jmt),        photoFe_diaz_tavg(imt,jmt), &
         photoFe_diat_tavg(imt,jmt),    photoFe_sp_tavg(imt,jmt), &
         FvPE_DIC_tavg(imt,jmt),        FvPE_ALK_tavg(imt,jmt), &
         NITRIF_tavg(imt,jmt),          DENITRIF_tavg(imt,jmt))

      cnt = 0

      XKW_tavg             = c0    ; cnt = cnt + 1
      PV_CO2_tavg          = c0    ; cnt = cnt + 1
      PV_O2_tavg           = c0    ; cnt = cnt + 1
      AP_tavg              = c0    ; cnt = cnt + 1
      SCHMIDT_O2_tavg      = c0    ; cnt = cnt + 1
      O2SAT_tavg           = c0    ; cnt = cnt + 1
      FG_O2_tavg           = c0    ; cnt = cnt + 1
      SCHMIDT_CO2_tavg     = c0    ; cnt = cnt + 1
      PH_tavg              = c0    ; cnt = cnt + 1
      CO2STAR_tavg         = c0    ; cnt = cnt + 1
      DCO2STAR_tavg        = c0    ; cnt = cnt + 1
      pCO2SURF_tavg        = c0    ; cnt = cnt + 1
      DpCO2_tavg           = c0    ; cnt = cnt + 1
      FG_CO2_tavg          = c0    ; cnt = cnt + 1
      IRON_FLUX_tavg       = c0    ; cnt = cnt + 1
      PROD_tavg            = c0    ; cnt = cnt + 1
      PO4_RESTORE_tavg     = c0    ; cnt = cnt + 1
      NO3_RESTORE_tavg     = c0    ; cnt = cnt + 1
      SiO3_RESTORE_tavg    = c0    ; cnt = cnt + 1
      PAR_avg_tavg         = c0    ; cnt = cnt + 1
      PO4STAR_tavg         = c0    ; cnt = cnt + 1
      POC_FLUX_IN_tavg     = c0    ; cnt = cnt + 1
      POC_PROD_tavg        = c0    ; cnt = cnt + 1
      POC_REMIN_tavg       = c0    ; cnt = cnt + 1
      CaCO3_FLUX_IN_tavg   = c0    ; cnt = cnt + 1
      CaCO3_PROD_tavg      = c0    ; cnt = cnt + 1
      CaCO3_REMIN_tavg     = c0    ; cnt = cnt + 1
      SiO2_FLUX_IN_tavg    = c0    ; cnt = cnt + 1
      SiO2_PROD_tavg       = c0    ; cnt = cnt + 1
      SiO2_REMIN_tavg      = c0    ; cnt = cnt + 1
      dust_FLUX_IN_tavg    = c0    ; cnt = cnt + 1
      dust_REMIN_tavg      = c0    ; cnt = cnt + 1
      P_iron_FLUX_IN_tavg  = c0    ; cnt = cnt + 1
      P_iron_PROD_tavg     = c0    ; cnt = cnt + 1
      P_iron_REMIN_tavg    = c0    ; cnt = cnt + 1
      graze_sp_tavg        = c0    ; cnt = cnt + 1
      graze_diat_tavg      = c0    ; cnt = cnt + 1
      graze_tot_tavg       = c0    ; cnt = cnt + 1
      sp_loss_tavg         = c0    ; cnt = cnt + 1
      diat_loss_tavg       = c0    ; cnt = cnt + 1
      zoo_loss_tavg        = c0    ; cnt = cnt + 1
      sp_agg_tavg          = c0    ; cnt = cnt + 1
      diat_agg_tavg        = c0    ; cnt = cnt + 1
      photoC_sp_tavg       = c0    ; cnt = cnt + 1
      photoC_diat_tavg     = c0    ; cnt = cnt + 1
      tot_prod_tavg        = c0    ; cnt = cnt + 1
      DOC_prod_tavg        = c0    ; cnt = cnt + 1
      DOC_remin_tavg       = c0    ; cnt = cnt + 1
      Fe_scavenge_tavg     = c0    ; cnt = cnt + 1
      sp_N_lim_tavg        = c0    ; cnt = cnt + 1
      sp_Fe_lim_tavg       = c0    ; cnt = cnt + 1
      sp_PO4_lim_tavg      = c0    ; cnt = cnt + 1
      sp_light_lim_tavg    = c0    ; cnt = cnt + 1
      diat_N_lim_tavg      = c0    ; cnt = cnt + 1
      diat_Fe_lim_tavg     = c0    ; cnt = cnt + 1
      diat_PO4_lim_tavg    = c0    ; cnt = cnt + 1
      diat_SiO3_lim_tavg   = c0    ; cnt = cnt + 1
      diat_light_lim_tavg  = c0    ; cnt = cnt + 1
      CaCO3_form_tavg      = c0    ; cnt = cnt + 1
      diaz_Nfix_tavg       = c0    ; cnt = cnt + 1
      graze_diaz_tavg      = c0    ; cnt = cnt + 1
      diaz_loss_tavg       = c0    ; cnt = cnt + 1
      photoC_diaz_tavg     = c0    ; cnt = cnt + 1
      diaz_P_lim_tavg      = c0    ; cnt = cnt + 1
      diaz_Fe_lim_tavg     = c0    ; cnt = cnt + 1
      diaz_light_lim_tavg  = c0    ; cnt = cnt + 1
      Fe_scavenge_rate_tavg = c0   ; cnt = cnt + 1
      DON_prod_tavg        = c0    ; cnt = cnt + 1
      DON_remin_tavg       = c0    ; cnt = cnt + 1
      DOFe_prod_tavg       = c0    ; cnt = cnt + 1
      DOFe_remin_tavg      = c0    ; cnt = cnt + 1
      DOP_prod_tavg        = c0    ; cnt = cnt + 1
      DOP_remin_tavg       = c0    ; cnt = cnt + 1
      bSi_form_tavg        = c0    ; cnt = cnt + 1
      photoFe_diat_tavg    = c0    ; cnt = cnt + 1
      photoFe_diaz_tavg    = c0    ; cnt = cnt + 1
      photoFe_sp_tavg      = c0    ; cnt = cnt + 1
      FvPE_DIC_tavg        = c0    ; cnt = cnt + 1
      FvPE_ALK_tavg        = c0    ; cnt = cnt + 1
      NITRIF_tavg          = c0    ; cnt = cnt + 1
      DENITRIF_tavg        = c0    ; cnt = cnt + 1

      end subroutine ecosys_init

!***********************************************************************

      subroutine ecosys_set_interior(k,TEMP,SHF_QSW_Wpm2,TRACER_MODULE, &
          DTRACER_MODULE, KMT, dz, dzr, zt, eco_inject, par_phyto)

!-----------------------------------------------------------------------
!     arguments
!-----------------------------------------------------------------------

      integer(kind=int_kind), intent(in) :: &
        k                   ! vertical level index

      real(kind=dbl_kind), dimension(imt,jmt), intent(in) ::  &
        TEMP,               & ! potential temperature (C)
        SHF_QSW_Wpm2             ! penetrative solar heat flux (W/m^2)

      real(kind=dbl_kind), dimension(imt,jmt,ecosys_tracer_cnt), &
        intent(in) :: TRACER_MODULE    ! current tracer values

      real(kind=dbl_kind), dimension(imt,jmt,ecosys_tracer_cnt), &
        intent(out) :: DTRACER_MODULE  ! computed source/sink terms

      integer(kind=int_kind) :: KMT ! maximum valid depth layer

      real(kind=dbl_kind), dimension(km) :: &
          dz,        & ! layer thickness
          dzr,       & ! inverse later thickness
          zt          ! layer midpoint position
      real(kind=dbl_kind), dimension (ecosys_tracer_cnt) :: &
          eco_inject
       real(KIND=dbl_kind) :: par_phyto ! PAR that this layer experiences

!-----------------------------------------------------------------------
!     local variables
!-----------------------------------------------------------------------

      character(len=*), parameter :: &
        sub_name = 'ecosys_mod:ecosys_set_interior'

      real(kind=dbl_kind), parameter :: &
        epsC      = 1.00e-8_dbl_kind, & !small C concentration (mmol C/m^3)
        epsTinv   = 3.17e-8_dbl_kind, & !small inverse time scale (1/year) (1/sec)
        epsnondim = 1.00e-6_dbl_kind  !small non-dimensional number (non-dim)

      type(sinking_particle), save :: &
        POC,          & ! base units = nmol C
        P_CaCO3,      & ! base units = nmol CaCO3
        P_SiO2,       & ! base units = nmol SiO2
        dust,         & ! base units = g
        P_iron        ! base units = nmol Fe

      real(kind=dbl_kind), dimension(imt,jmt), save :: &
        QA_dust_def,  & ! incoming deficit in the QA(dust) POC flux
        PAR_out       ! photosynthetically available radiation (W/m^2)

      real(kind=dbl_kind), dimension(imt,jmt) :: &
        PO4_loc,      & ! local copy of model PO4
        NO3_loc,      & ! local copy of model NO3
        SiO3_loc,     & ! local copy of model SiO3
        NH4_loc,      & ! local copy of model NH4
        Fe_loc,       & ! local copy of model Fe
        O2_loc,       & ! local copy of model O2
        DOC_loc,      & ! local copy of model DOC
        spC_loc,      & ! local copy of model spC
        spChl_loc,    & ! local copy of model spChl
        spCaCO3_loc,  & ! local copy of model spCaCO3
        diatC_loc,    & ! local copy of model diatC
        diatChl_loc,  & ! local copy of model diatChl
        zooC_loc,     & ! local copy of model zooC
        spFe_loc,     & ! local copy of model spFe
        diatSi_loc,   & ! local copy of model diatSi
        diatFe_loc,   & ! local copy of model diatFe
        diazC_loc,    & ! local copy of model diazC
        diazChl_loc,  & ! local copy of model diazChl
        diazFe_loc,   & ! local copy of model diazFe
        DON_loc,      & ! local copy of model DON
        DOFe_loc,     & ! local copy of model DOFe
        DOP_loc       ! local copy of model DOP

      real(kind=dbl_kind) :: &
        z_grz_sqr,    & ! square of parm_z_grz (mmol C/m^3)^2
        C_loss_thres, & ! bio-C threshold at which losses go to zero (mmol C/m^3)
        f_loss_thres  ! fraction of grazing loss reduction at depth

      real(kind=dbl_kind), dimension(imt,jmt) :: &
        PAR_in,       & ! photosynthetically available radiation (W/m^2)
        KPARdz,       & ! PAR adsorption coefficient (non-dim)
        PAR_avg,      & ! average PAR over mixed layer depth (W/m^2)
        DOC_prod,     & ! production of DOC (mmol C/m^3/sec)
        DOC_remin,    & ! remineralization of DOC (mmol C/m^3/sec)
        NITRIF,       & ! nitrification (NH4 -> NO3) (mmol N/m^3/sec)
        DENITRIF,     & ! denitrification (NO3 -> N2) (mmol N/m^3/sec)
        RESTORE      ! restoring terms for nutrients (mmol ./m^3/sec)


      real(kind=dbl_kind), dimension(imt,jmt) :: &
        z_umax,       & ! max. zoo growth rate on sp at local T (1/sec)
        diat_umax,    & ! max. zoo growth rate on diatoms at local T (1/sec)
        z_mort,       & ! zoo respiration loss, (1/sec/((mmol C/m3))
        C_loss_diaz,  & ! bio-C threshold at which losses go to zero (mmol C/m^3)
        z_mort2,      & ! zoo quad mort rate, tohigherlevels (1/sec/((mmol C/m3))
        diaz_umax     ! max. zoo growth rate on diazotrophs at local T (1/sec)

      real(kind=dbl_kind), dimension(imt,jmt) :: &
        thetaC_sp,    & ! local Chl/C ratio in small phyto (mg Chl/mmol C)
        thetaC_diat,  & ! local Chl/C ratio in diatoms (mg Chl/mmol C)
        QCaCO3,       & ! small phyto CaCO3/C ratio (mmol CaCO3/mmol C)
        Tfunc,        & ! temp response function GD98 (non-dim)
        VNO3_sp,      & ! small phyto NO3 uptake rate (non-dim)
        VNH4_sp,      & ! small phyto NH4 uptake rate (non-dim)
        VNtot_sp,     & ! small phyto total N uptake rate (non-dim)
        VFeC_sp,      & ! ??? small phyto C-specific iron uptake (non-dim)
        VPO4_sp,      & ! ??? (non-dim)
        f_nut,        & ! nut limitation factor, modifies C fixation (non-dim)
        PCmax,        & ! max value of PCphoto at temperature TEMP (1/sec)
        PCphoto_sp,   & ! small C-specific rate of photosynth. (1/sec)
        photoC_sp,    & ! small phyto C-fixation (mmol C/m^3/sec)
        NO3_V_sp,     & ! nitrate uptake by small phyto (mmol NO3/m^3/sec)
        NH4_V_sp,     & ! ammonium uptake by small phyto (mmol NH4/m^3/sec)
        VNC_sp,       & ! small phyto C-specific N uptake rate (mmol N/mmol C/sec)
        pChl,         & ! Chl synth. regulation term (mg Chl/mmol N)
        photoacc_sp,  & ! Chl synth. term in photoadapt. (GD98) (mg Chl/m^3/sec)
        CaCO3_prod,   & ! prod. of CaCO3 by small phyto (mmol CaCO3/m^3/sec)
        VNO3_diat,    & ! diatom nitrate uptake rate (non-dim)
        VNH4_diat,    & ! diatom ammonium uptake rate (non-dim)
        VNtot_diat,   & ! diatom total N uptake rate (non-dim)
        VFeC_diat,    & ! diatom C-specific iron uptake (non-dim)
        VPO4_diat,    & ! diatom C-specific PO4 uptake (non-dim)
        VSiO3_diat,   & ! C-specific SiO3 uptake for diatoms (non-dim)
        PCphoto_diat, & ! diatom C-specific rate of photosynth. (1/sec)
        photoC_diat,  & ! diatom C-fixation (mmol C/m^3/sec)
        NO3_V_diat,   & ! nitrate uptake by diatoms (mmol NO3/m^3/sec)
        NH4_V_diat,   & ! ammonium uptake by diatoms (mmol NH4/m^3/sec)
        VNC_diat,     & ! diatom C-specific N uptake rate (mmol N/mmol C/sec)
        photoacc_diat,& ! Chl synth. term in photoadapt. (GD98) (mg Chl/m^3/sec)
        reduceV,      & ! factor in nutrient uptake (mmol C/m^3)^2
        graze_sp,     & ! grazing rate on small phyto (mmol C/m^3/sec)
        graze_sp_zoo, & ! graze_sp routed to zoo (mmol C/m^3/sec)
        graze_sp_poc, & ! graze_sp routed to poc (mmol C/m^3/sec)
        graze_sp_doc, & ! graze_sp routed to doc (mmol C/m^3/sec)
        graze_sp_dic  ! graze_sp routed to dic (mmol C/m^3/sec)

      real(kind=dbl_kind), dimension(imt,jmt) :: & ! max of 39 continuation lines
        graze_diat,    & ! grazing rate on diatoms (mmol C/m^3/sec)
        graze_diat_zoo,& ! graze_diat routed to zoo (mmol C/m^3/sec)
        graze_diat_poc,& ! graze_diat routed to poc (mmol C/m^3/sec)
        graze_diat_doc,& ! graze_diat routed to doc (mmol C/m^3/sec)
        graze_diat_dic,& ! graze_diat routed to dic (mmol C/m^3/sec)
        Pprime,        & ! used to limit phyto mort at low biomass (mmol C/m^3)
        sp_loss,       & ! small phyto non-grazing mort (mmol C/m^3/sec)
        sp_loss_poc,   & ! sp_loss routed to poc (mmol C/m^3/sec)
        sp_loss_doc,   & ! sp_loss routed to doc (mmol C/m^3/sec)
        sp_loss_dic,   & ! sp_loss routed to dic (mmol C/m^3/sec)
        sp_agg,        & ! small phyto agg loss (mmol C/m^3/sec)
        diat_loss,     & ! diatom non-grazing mort (mmol C/m^3/sec)
        diat_loss_poc, & ! diat_loss routed to poc (mmol C/m^3/sec)
        diat_loss_doc, & ! diat_loss routed to doc (mmol C/m^3/sec)
        diat_loss_dic, & ! diat_loss routed to dic (mmol C/m^3/sec)
        diat_agg,      & ! diatom agg (mmol C/m^3/sec)
        f_zoo_detr,    & ! frac of zoo losses into large detrital pool (non-dim)
        Fe_scavenge,   & ! loss of dissolved iron, scavenging (mmol Fe/m^3/sec)
        Zprime,        & ! used to limit zoo mort at low biomass (mmol C/m^3)
        zoo_loss,      & ! mortality on zooplankton (mmol C/m^3/sec)
        zoo_loss_doc,  & ! zoo_loss routed to doc (mmol C/m^3/sec)
        zoo_loss_dic,  & ! zoo_loss routed to dic (mmol C/m^3/sec)
        WORK,          & ! intermediate value in photsyntheis computation (1/sec)
        light_lim,     & ! light limitation factor
        Qsi,           & ! Diatom initial Si/C ratio (mmol Si/mmol C)
        gQsi,          & ! diatom Si/C ratio for growth (new biomass)
        Qfe_sp,        & ! small phyto init fe/C ratio (mmolFe/mmolC)
        gQfe_sp,       & ! small phyto fe/C for growth
        Qfe_diat,      & ! diatom init fe/C ratio
        gQfe_diat,     & ! diatom fe/C ratio for growth
        Qfe_diaz,      & ! diazotrophs init fe/C ratio
        gQfe_diaz      ! diazotroph fe/C ratio for new growth

      real(kind=dbl_kind), dimension(imt,jmt) :: & ! max of 39 continuation lines
        PCphoto_diaz,  & ! diazotroph C-specific rate of photosynth. (1/sec)
        photoC_diaz,   & ! diazotroph C-fixation (mmol C/m^3/sec)
        Vfec_diaz,     & ! diazotroph C-specific iron uptake (non-dim)
        Vpo4_diaz,     & ! diazotroph C-specific po4 uptake (non-dim)
        photoacc_diaz, & ! Chl synth. term in photoadapt. (GD98) (mg Chl/m^3/sec)
        Vnc_diaz,      & ! diazotroph C-specific N uptake rate (mmol N/mmol C/sec)
        diaz_Nexcrete, & ! diazotroph fixed N excretion
        diaz_Nfix,     & ! diazotroph total Nitrogen fixation (mmol N/m^3/sec)
        thetaC_diaz,   & ! local Chl/C ratio in diazotrophs (mg Chl/mmol C)
        photoFe_diaz,  & ! iron uptake by diazotrophs (mmolFe/m^3/sec)
        photoFe_diat,  & ! iron uptake by diatoms
        photoFe_sp,    & ! iron uptake by small phyto
        photoN_diaz,   & ! nitrogen for Nfixation added to diaz biomass
        photoSi_diat,  & ! silicon uptake by diatoms (mmol Si/m^3/sec)
        remaining_diazP,& ! used in routing P from diazotroph losses
        diaz_loss,     & ! diazotroph non-grazing mort (mmol C/m^3/sec)
        diaz_loss_doc, & ! mortality routed to DOM pool
        diaz_loss_dic, & ! mortality routed to remin
        diaz_loss_dop, & ! P from mort routed to DOP pool
        diaz_loss_dip, & ! P from mort routed to remin
             ! NOTE: P from diaz losses must be routed differently than
                       ! other elements to ensure that sinking detritus and
                       ! zooplankton pools get their fixed P/C ratios, the
                       ! remaining P is split evenly between DOP and PO4
        graze_diaz,    & ! grazing rate on diazotrophs (mmol C/m^3/sec)
        graze_diaz_poc,& ! grazing routed to sinking detr (mmol C/m^3/sec)
        graze_diaz_doc,& ! grazing routed to DOC (mmol C/m^3/sec)
        graze_diaz_dic,& ! grazing routed to remin (mmol C/m^3/sec)
        graze_diaz_zoo,& ! grazing routed to new zoo biomass
        DON_remin,     & ! portion of DON remineralized
        DOFe_remin,    & ! portion of DOFe remineralized
        DOP_remin,     & ! portion of DOP remineralized
        DOM_remin,     & ! fraction of DOM remineralized at current TEMP
        Fe_scavenge_rate ! annual scavenging rate of iron as % of ambient

      real(kind=dbl_kind), dimension(imt,jmt) :: & ! max of 39 continuation lines
        DON_prod,      & ! production of dissolved organic N
        DOFe_prod,     & ! produciton of dissolved organic Fe
        DOP_prod       ! production of dissolved organic P

!-----------------------------------------------------------------------

!      call timer_start(timer_interior)

!-----------------------------------------------------------------------

      DTRACER_MODULE = c0

!-----------------------------------------------------------------------
!     exit immediately if computations are not to be performed
!-----------------------------------------------------------------------

      if (.not. lsource_sink) then
 !       call timer_stop(timer_interior)
        return
      endif

!-----------------------------------------------------------------------
!     create local copies of model tracers
!     treat negative values as zero
!     apply mask to local copies
!-----------------------------------------------------------------------

      PO4_loc      = max(c0, TRACER_MODULE(:,:,po4_ind))
      NO3_loc      = max(c0, TRACER_MODULE(:,:,no3_ind))
      SiO3_loc     = max(c0, TRACER_MODULE(:,:,sio3_ind))
      NH4_loc      = max(c0, TRACER_MODULE(:,:,nh4_ind))
      Fe_loc       = max(c0, TRACER_MODULE(:,:,fe_ind))
      O2_loc       = max(c0, TRACER_MODULE(:,:,o2_ind))
      DOC_loc      = max(c0, TRACER_MODULE(:,:,doc_ind))
      spC_loc      = max(c0, TRACER_MODULE(:,:,spC_ind))
      spChl_loc    = max(c0, TRACER_MODULE(:,:,spChl_ind))
      spCaCO3_loc  = max(c0, TRACER_MODULE(:,:,spCaCO3_ind))
      diatC_loc    = max(c0, TRACER_MODULE(:,:,diatC_ind))
      diatChl_loc  = max(c0, TRACER_MODULE(:,:,diatChl_ind))
      zooC_loc     = max(c0, TRACER_MODULE(:,:,zooC_ind))

      !zooC_loc = c0

      spFe_loc     = max(c0, TRACER_MODULE(:,:,spFe_ind))
      diatSi_loc   = max(c0, TRACER_MODULE(:,:,diatSi_ind))
      diatFe_loc   = max(c0, TRACER_MODULE(:,:,diatFe_ind))
      diazC_loc    = max(c0, TRACER_MODULE(:,:,diazC_ind))
      diazChl_loc  = max(c0, TRACER_MODULE(:,:,diazChl_ind))
      diazFe_loc   = max(c0, TRACER_MODULE(:,:,diazFe_ind))
      DON_loc      = max(c0, TRACER_MODULE(:,:,don_ind))
      DOFe_loc     = max(c0, TRACER_MODULE(:,:,dofe_ind))
      DOP_loc      = max(c0, TRACER_MODULE(:,:,dop_ind))

      where (.not. LAND_MASK .or. k > KMT)
        PO4_loc      = c0
        NO3_loc      = c0
        SiO3_loc     = c0
        NH4_loc      = c0
        Fe_loc       = c0
        O2_loc       = c0
        DOC_loc      = c0
        spC_loc      = c0
        spChl_loc    = c0
        spCaCO3_loc  = c0
        diatC_loc    = c0
        diatChl_loc  = c0
        zooC_loc     = c0
        spFe_loc     = c0
        diatSi_loc   = c0
        diatFe_loc   = c0
        diazC_loc    = c0
        diazChl_loc  = c0
        diazFe_loc   = c0
        DON_loc      = c0
        DOFe_loc     = c0
        DOP_loc      = c0
      endwhere

!-----------------------------------------------------------------------
!     If any phyto box are zero, set others to zeros.
!-----------------------------------------------------------------------

      where (spC_loc == c0 .or. spChl_loc == c0 .or. spFe_loc == c0)
        spC_loc = c0
        spChl_loc = c0
        spCaCO3_loc = c0
        spFe_loc = c0
      endwhere

      where (diatC_loc == c0 .or. diatChl_loc == c0 .or. &
        diatFe_loc == c0 .or. diatSi_loc == c0)
        diatC_loc = c0
        diatChl_loc = c0
        diatFe_loc = c0
        diatSi_loc = c0
      endwhere

      where (diazC_loc == c0 .or. diazChl_loc == c0 .or. &
        diazFe_loc == c0)
        diazC_loc = c0
        diazChl_loc = c0
        diazFe_loc = c0
      endwhere

!-----------------------------------------------------------------------
!     set local variables, with incoming ratios
!-----------------------------------------------------------------------

      thetaC_sp   = spChl_loc   / (spC_loc + epsC)
      thetaC_diat = diatChl_loc / (diatC_loc + epsC)
      thetaC_diaz = diazChl_loc / (diazC_loc + epsC)
      Qsi         = diatSi_loc  / (diatC_loc + epsC)
      Qfe_diat    = diatFe_loc  / (diatC_loc + epsC)
      Qfe_sp      = spFe_loc    / (spC_loc + epsC)
      Qfe_diaz    = diazFe_loc  / (diazC_loc + epsC)
      where (Qsi > gQsi_max) Qsi = gQsi_max

!-----------------------------------------------------------------------
!     determine new elemental ratios for growth (new biomass)
!-----------------------------------------------------------------------

      gQsi      = gQsi_0
      gQfe_diat = gQfe_diat_0
      gQfe_sp   = gQfe_sp_0
      gQfe_diaz = gQfe_diaz_0

!-----------------------------------------------------------------------
!     modify these initial ratios under low ambient iron conditions
!-----------------------------------------------------------------------

      where (Fe_loc .lt. c2 * parm_diat_kfe)
        gQfe_diat = max((gQfe_diat * Fe_loc /(c2 * parm_diat_kfe)), &
            gQfe_diat_min)
      endwhere

      where ((Fe_loc .lt. c2 * parm_diat_kfe) .and. (Fe_loc > c0) .and. &
        (SiO3_loc .gt. (c2 * parm_diat_kSiO3)))
            gQsi = min(((gQsi*gQsi_coef*c2*parm_diat_kfe/Fe_loc) &
            - (gQsi_coef-c1)*gQsi_0), gQsi_max)
      endwhere

      where (Fe_loc == c0)
        gQsi = gQsi_max
      endwhere

      where (Fe_loc .lt. c2 * parm_sp_kfe)
        gQfe_sp = max((gQfe_sp*Fe_loc/(c2 * parm_sp_kfe)), &
            gQfe_sp_min)
      endwhere

      where (Fe_loc .lt. c2 * diaz_kFe)
        gQfe_diaz = max((gQfe_diaz*Fe_loc/(c2 * diaz_kFe)), &
            gQfe_diaz_min)
      endwhere

!-----------------------------------------------------------------------
!     Modify the initial si/C ratio under low ambient Si conditions
!-----------------------------------------------------------------------

      where (SiO3_loc .lt. (c2 * parm_diat_kSiO3))
        gQsi = max((gQsi*SiO3_loc/(c2 * parm_diat_kSiO3)), &
            gQsi_min)
      endwhere

!-----------------------------------------------------------------------
!     various k==1 initializations
!
!     0.45   fraction of incoming SW -> PAR (non-dim)
!-----------------------------------------------------------------------

      if (k == 1) then
        where (LAND_MASK)
           PAR_out = max(c0, f_qsw_par * SHF_QSW_Wpm2)
        elsewhere
           PAR_out = c0
        endwhere
        call init_particulate_terms(POC, P_CaCO3, P_SiO2, dust, P_iron, &
            QA_dust_def)
      endif

!-----------------------------------------------------------------------
!     QCaCO3 is the percentage of sp organic matter which is associated
!     with coccolithophores
!-----------------------------------------------------------------------

      QCaCO3 = spCaCO3_loc / (spC_loc + epsC)
      where (QCaCO3 > QCaCO3_max) QCaCO3 = QCaCO3_max

!-----------------------------------------------------------------------
!     compute PAR related quantities
!
!   0.03e-2   atten. coeff. per unit chlorophyll (1/cm/(mg Chl/m^3))
!   0.04e-2   atten. coeff. for water (1/cm)
!-----------------------------------------------------------------------

      PAR_in = PAR_out
      where (.not. LAND_MASK .or. k > KMT .or. PAR_out < 1.0e-6)
        PAR_in = c0
      end where ! added PAR_out test to prevent floating-point under-run - Ben Saenz 9/2011

      KPARdz = (k_chl * (spChl_loc + diatChl_loc &
        + diazChl_loc) + k_h2o) * dz(k)

      PAR_out = PAR_in * exp(-KPARdz)
      PAR_avg = PAR_in * (c1 - exp(-KPARdz)) / KPARdz

      !print *, k,'- PAR_out', PAR_out, 'PAR_avg', PAR_avg
      PAR_avg = par_phyto                  ! ben added 6/2013 b/c above was not averaging PAR over mixed layer

!-----------------------------------------------------------------------
!   Tref = 30.0 reference temperature (deg. C)
!
!   Using q10 formulation with Q10 value of 2.0 (Doney et al., 1996).
!-----------------------------------------------------------------------

      Tfunc = Q_10**(((TEMP + T0_Kelvin) &
        - (Tref + T0_Kelvin)) / c10)

!-----------------------------------------------------------------------
!   modify growth mort rates by Tfunc
!-----------------------------------------------------------------------

      z_umax    = parm_z_umax_0    * Tfunc
      diat_umax = parm_diat_umax_0 * Tfunc
      z_mort2   = parm_z_mort2_0   * Tfunc
      z_mort    = parm_z_mort_0    * Tfunc
      diaz_umax = parm_diaz_umax_0 * Tfunc

      DOM_remin = parm_sd_remin_0

!-----------------------------------------------------------------------
!     Get relative nutrient uptake rates for phytoplankton,
!     min. relative uptake rate modifies C fixation in the manner
!     that the min. cell quota does in GD98.
!-----------------------------------------------------------------------

      VNO3_sp = (NO3_loc / parm_sp_kNO3) &
        / (c1 + (NO3_loc / parm_sp_kNO3) + (NH4_loc / parm_sp_kNH4))

      VNH4_sp = (NH4_loc / parm_sp_kNH4) &
        / (c1 + (NO3_loc / parm_sp_kNO3) + (NH4_loc / parm_sp_kNH4))

      VNtot_sp = VNO3_sp + VNH4_sp
      sp_N_lim_tavg = VNtot_sp

!-----------------------------------------------------------------------
!     get relative Fe uptake by phytoplankton
!     get relative P uptake rates for phytoplankton
!-----------------------------------------------------------------------

      VFeC_sp = Fe_loc / (Fe_loc + parm_sp_kFe)
      sp_Fe_lim_tavg = VFeC_sp
      VPO4_sp = PO4_loc / (PO4_loc + parm_sp_kPO4)
      sp_PO4_lim_tavg = VPO4_sp

!-----------------------------------------------------------------------
!   Small Phytoplankton C-fixation - given light and relative uptake rates
!   determine small phyto nutrient limitation factor for carbon fixation
!-----------------------------------------------------------------------

      f_nut = min(VNtot_sp, VFeC_sp)
      f_nut = min(f_nut, VPO4_sp)

!-----------------------------------------------------------------------
!   get small phyto photosynth. rate, phyto C biomass change, photoadapt
!-----------------------------------------------------------------------

      PCmax = PCrefSp * f_nut * Tfunc

      light_lim = (c1 - exp((-c1 * parm_alphaChlsp * thetaC_sp * PAR_avg) &
        / (PCmax + epsTinv)))
      PCphoto_sp = PCmax * light_lim
      sp_light_lim_tavg = light_lim

      photoC_sp = PCphoto_sp * spC_loc

!-----------------------------------------------------------------------
!     Get nutrient uptakes by small phyto based on calculated C fixation
!     total N uptake (VNC_sp) is used in photoadaption
!-----------------------------------------------------------------------

      where (VNtot_sp > c0)
        NO3_V_sp = (VNO3_sp / VNtot_sp) * photoC_sp * Q
        NH4_V_sp = (VNH4_sp / VNtot_sp) * photoC_sp * Q
        VNC_sp = PCphoto_sp * Q
      elsewhere
        NO3_V_sp = c0
        NH4_V_sp = c0
        VNC_sp = c0
        photoC_sp = c0
      endwhere

      photoFe_sp = photoC_sp * gQfe_sp
      photoFe_sp_tavg = photoFe_sp

!-----------------------------------------------------------------------
!     calculate pChl, (used in photoadapt., GD98)
!     2.3   max value of thetaN (Chl/N ratio) (mg Chl/mmol N)
!     GD 98 Chl. synth. term
!-----------------------------------------------------------------------

      WORK = parm_alphaChlsp * thetaC_sp * PAR_avg
      where (WORK > c0)
        pChl = thetaN_max_sp * PCphoto_sp / WORK
        photoacc_sp = (pChl * VNC_sp / thetaC_sp) * spChl_loc
      elsewhere
        photoacc_sp = c0
      endwhere

!-----------------------------------------------------------------------
!     CaCO3 Production, parameterized as function of small phyto production
!     decrease CaCO3 as function of nutrient limitation decrease CaCO3 prod
!     at low temperatures increase CaCO3 prod under bloom conditions
!     maximum calcification rate is 40% of primary production
!-----------------------------------------------------------------------

      CaCO3_prod = f_prod_sp_CaCO3 * photoC_sp
      CaCO3_prod = CaCO3_prod * f_nut

      where (TEMP < CaCO3_temp_thres1) &
        CaCO3_prod = CaCO3_prod * max((TEMP-CaCO3_temp_thres2), c0) &
        / (CaCO3_temp_thres1-CaCO3_temp_thres2)

      where (spC_loc > CaCO3_sp_thres) &
        CaCO3_prod = min((CaCO3_prod*spC_loc/CaCO3_sp_thres), &
        (f_photosp_CaCO3*photoC_sp))

      CaCO3_form_tavg = CaCO3_prod

!-----------------------------------------------------------------------
!   Relative uptake rates for diatoms nitrate is VNO3, ammonium is VNH4
!-----------------------------------------------------------------------

      VNO3_diat = (NO3_loc / parm_diat_kNO3) &
        / (c1 + (NO3_loc / parm_diat_kNO3) + (NH4_loc / parm_diat_kNH4))

      VNH4_diat = (NH4_loc / parm_diat_kNH4) &
        / (c1 + (NO3_loc / parm_diat_kNO3) + (NH4_loc / parm_diat_kNH4))

      VNtot_diat = VNO3_diat + VNH4_diat
      diat_N_lim_tavg = VNtot_diat

!-----------------------------------------------------------------------
!     get relative Fe uptake by diatoms
!     get relative P uptake rates for diatoms
!     get relative SiO3 uptake rate for diatoms
!-----------------------------------------------------------------------

      VFeC_diat = Fe_loc / (Fe_loc + parm_diat_kFe)
      diat_Fe_lim_tavg = VFeC_diat
      VPO4_diat = PO4_loc / (PO4_loc + parm_diat_kPO4)
      diat_PO4_lim_tavg = VPO4_diat
      VSiO3_diat = SiO3_loc / (SiO3_loc + parm_diat_kSiO3)
      diat_SiO3_lim_tavg = VSiO3_diat

!-----------------------------------------------------------------------
!     Diatom carbon fixation and photoadapt.
!     determine diatom nutrient limitation factor for carbon fixation
!-----------------------------------------------------------------------

      f_nut = min(VNtot_diat, VFeC_diat)
      f_nut = min(f_nut, VSiO3_diat)
      f_nut = min(f_nut, VPO4_diat)

!-----------------------------------------------------------------------
!     get diatom photosynth. rate, phyto C biomass change, photoadapt
!-----------------------------------------------------------------------

      PCmax = PCrefDiat * f_nut * Tfunc

      light_lim = &
        (c1 - exp((-c1 * parm_alphaChldiat * thetaC_diat * PAR_avg) / &
        (PCmax + epsTinv)))
      PCphoto_diat = PCmax * light_lim
      diat_light_lim_tavg = light_lim

      photoC_diat = PCphoto_diat * diatC_loc

!-----------------------------------------------------------------------
!     Get nutrient uptake by diatoms based on C fixation
!-----------------------------------------------------------------------

      where (VNtot_diat > c0)
        NO3_V_diat = (VNO3_diat / VNtot_diat) * photoC_diat * Q
        NH4_V_diat = (VNH4_diat / VNtot_diat) * photoC_diat * Q
        VNC_diat = PCphoto_diat * Q
      elsewhere
        NO3_V_diat = c0
        NH4_V_diat = c0
        VNC_diat = c0
        photoC_diat = c0
      endwhere

      photoFe_diat = photoC_diat * gQfe_diat
      photoSi_diat = photoC_diat * gQsi

      photoFe_diat_tavg = photoFe_diat
      bSi_form_tavg = photoSi_diat

!-----------------------------------------------------------------------
!     calculate pChl, (used in photoadapt., GD98)
!     3.0   max value of thetaN (Chl/N ratio) (mg Chl/mmol N)
!     GD 98 Chl. synth. term
!-----------------------------------------------------------------------

      WORK = parm_alphaChldiat * thetaC_diat * PAR_avg  ! (mmol C m^2/(mg Chl W sec)) * mg Chl/(mmol C) * W --> m^2/sec
      where (WORK > c0)
        pChl = thetaN_max_diat * PCphoto_diat / WORK    ! mg Chl/mmol N * mmol C / (m^2/s) = mmol C * sec / (mmol N m^2)
        photoacc_diat = (pChl * VNC_diat / thetaC_diat) * diatChl_loc  ! mmol C * sec / (mmol N m^2) * mmol N * / (mg Chl / mmol C) * mmol C =

        ! Chl/N * C / (alpha * Chl/C * PAR) * C*N/C   / Chl/C * tot_C
        ! Chl/N * C / alpha * C/Chl / PAR * C*N/C  * C/Chl * tot_C
        !

      elsewhere
        photoacc_diat = c0
      endwhere

!-----------------------------------------------------------------------
!     get relative Fe uptake by diazotrophs
!     get relative P uptake rates for diazotrophs
!-----------------------------------------------------------------------

      Vfec_diaz = Fe_loc / (Fe_loc + diaz_kFe)
      diaz_Fe_lim_tavg = Vfec_diaz

      Vpo4_diaz = PO4_loc / (PO4_loc + diaz_kPO4)
      diaz_P_lim_tavg = Vpo4_diaz

      f_nut = min(Vpo4_diaz, Vfec_diaz)

!-----------------------------------------------------------------------
!   get diazotroph photosynth. rate, phyto C biomass change
!-----------------------------------------------------------------------

      PCmax = PCrefDiaz * f_nut * Tfunc

      light_lim = &
        (c1 - exp((-c1 * parm_alphaDiaz * thetaC_diaz * PAR_avg) / &
        (PCmax + epsTinv)))
      PCphoto_diaz = PCmax * light_lim
      diaz_light_lim_tavg = light_lim

      photoC_diaz = PCphoto_diaz * diazC_loc

!-----------------------------------------------------------------------
!     Get N fixation by diazotrophs based on C fixation,
!     Diazotrophs fix more than they need then 30% is excreted
!-----------------------------------------------------------------------

      diaz_Nfix = photoC_diaz * Q * r_Nfix_photo
      photoN_diaz = photoC_diaz * Q
      diaz_Nexcrete = diaz_Nfix - photoN_diaz
      diaz_Nfix_tavg = diaz_Nfix

      Vnc_diaz = PCphoto_diaz * Q

!-----------------------------------------------------------------------
!     Get Fe and po4 uptake by diazotrophs based on C fixation
!-----------------------------------------------------------------------

      photoFe_diaz = photoC_diaz * gQfe_diaz
      photoFe_diaz_tavg = photoFe_diaz

!-----------------------------------------------------------------------
!     calculate pChl, (used in photoadapt., GD98)
!     3.4   max value of thetaN (Chl/N ratio) (mg Chl/mmol N)
!     GD 98 Chl. synth. term
!-----------------------------------------------------------------------

      WORK = parm_alphaDiaz * thetaC_diaz * PAR_avg
      where (WORK > c0)
        pChl = thetaN_max_diaz * PCphoto_diaz / WORK
        photoacc_diaz = (pChl * Vnc_diaz / thetaC_diaz) * diazChl_loc
      elsewhere
        photoacc_diaz = c0
      endwhere

!-----------------------------------------------------------------------
!     CALCULATE GRAZING AND OTHER MORT
!-----------------------------------------------------------------------

!---------------------------------------------------------------------
!     calculate the loss threshold interpolation factor
!---------------------------------------------------------------------
         if (zt(k) > thres_z1) then
            if (zt(k) < thres_z2) then
               f_loss_thres = (thres_z2 - zt(k))/(thres_z2 - thres_z1)
            else
               f_loss_thres = c0
            endif
         else
            f_loss_thres = c1
         endif

!-----------------------------------------------------------------------
!     0.001 small phytoplankton threshold C concentration (mmol C/m^3)
!     get small phyto loss(in C units)
!     small phyto agg loss
!     get grazing rate (graze_sp) on small phyto  (in C units)
!-----------------------------------------------------------------------

      C_loss_thres = f_loss_thres * loss_thres_sp

      Pprime = max(spC_loc - C_loss_thres, c0)

      sp_loss = sp_mort * Pprime

      sp_agg = min((sp_agg_rate_max * dps) * Pprime, sp_mort2 * Pprime  &
         * Pprime)

      reduceV = Pprime * Pprime
      z_grz_sqr = parm_z_grz * parm_z_grz
      graze_sp = z_umax * zooC_loc * (reduceV / (reduceV + z_grz_sqr))

!-----------------------------------------------------------------------
!     routing of graze_sp & sp_loss
!     sp_agg all goes to POC
! currently assumes that 33% of grazed caco3 is remineralized
! if z_ingest ever changes, coefficients on routing grazed sp must change!
! min.%C routed to POC from grazing for ballast requirements = 0.4 * Qcaco3
! min.%C routed from sp_loss = 0.59 * QCaCO3, or P_CaCO3%rho
!-----------------------------------------------------------------------

      graze_sp_zoo = z_ingest * graze_sp
      graze_sp_poc = graze_sp * max((caco3_poc_min * QCaCO3), &
        min((spc_poc_fac * Pprime),f_graze_sp_poc_lim))

      graze_sp_doc = f_graze_sp_doc * graze_sp - graze_sp_poc
      graze_sp_dic = f_graze_sp_dic * graze_sp

      sp_loss_poc = QCaCO3 * sp_loss
      sp_loss_doc = (c1 - parm_labile_ratio)*(sp_loss-sp_loss_poc)
      sp_loss_dic = parm_labile_ratio * (sp_loss - sp_loss_poc)

!-----------------------------------------------------------------------
!     0.01 small diatom threshold C concentration (mmol C/m^3)
!     get diatom loss(in C units)
!     Diatom agg loss, min. 1%/day
!     get grazing rate (graze_diat) on diatoms  (in C units)
!-----------------------------------------------------------------------

      C_loss_thres = f_loss_thres * loss_thres_diat

      Pprime = max(diatC_loc - C_loss_thres, c0)

      diat_loss = diat_mort * Pprime

      diat_agg = min((diat_agg_rate_max * dps) * Pprime, &
        diat_mort2 * Pprime * Pprime)
      diat_agg = max((diat_agg_rate_min * dps) * Pprime, diat_agg)

!------------------------------------------------------------------------
!   Lower z_grz term for diatoms and diazotrophs, larger, more mobile
!   predators
!------------------------------------------------------------------------

      reduceV = Pprime * Pprime
      graze_diat = diat_umax * zooC_loc * &
        (reduceV / (reduceV + z_grz_sqr * f_z_grz_sqr_diat))

!-----------------------------------------------------------------------
!     routing of graze_diat & diat_loss
!     diat_agg all goes to POC
!     NOTE: if z_ingest is changed, coeff.s for poc,doc and dic must change!
!-----------------------------------------------------------------------

      graze_diat_zoo = z_ingest * graze_diat
      graze_diat_poc = f_graze_diat_poc * graze_diat
      graze_diat_doc = f_graze_diat_doc * graze_diat
      graze_diat_dic = f_graze_diat_dic * graze_diat

      diat_loss_poc = f_diat_loss_poc * diat_loss
      diat_loss_doc = (c1-parm_labile_ratio) * f_diat_loss_dc  &
        * diat_loss
      diat_loss_dic = parm_labile_ratio * f_diat_loss_dc * diat_loss

!-----------------------------------------------------------------------
!     0.03 small diazotroph threshold C concentration (mmol C/m^3)
!     Lower value used at temperatures < 16 deg. C, negligible biomass
!     get diazotroph loss(in C units)
!     get grazing rate (graze_diaz) on diazotrophs  (in C units)
!     no aggregation loss for diazotrophs
!-----------------------------------------------------------------------

      C_loss_diaz = f_loss_thres * loss_thres_diaz

      where (TEMP .lt. diaz_temp_thres) C_loss_diaz =  &
                f_loss_thres * loss_thres_diaz2

      Pprime = max(diazC_loc - C_loss_diaz, c0)

      diaz_loss = diaz_mort * Pprime

      reduceV = Pprime * Pprime
      graze_diaz = diaz_umax * zooC_loc &
        * (reduceV / (reduceV + z_grz_sqr))

!-----------------------------------------------------------------------
!     routing of graze_diaz & diaz_loss
!     NOTE: if z_ingest is changed, coeff.s for poc,doc and dic must change!
!     z_ingest for diaz = 0.21 based on O'Neil (1998)
!-----------------------------------------------------------------------

      graze_diaz_zoo = f_graze_diaz_zoo * graze_diaz
      graze_diaz_poc = f_graze_diaz_poc * graze_diaz
      graze_diaz_doc = f_graze_diaz_doc * graze_diaz
      graze_diaz_dic = f_graze_diaz_dic * graze_diaz

      diaz_loss_doc = (c1 - parm_labile_ratio) * diaz_loss
      diaz_loss_dic = parm_labile_ratio * diaz_loss

!--------------------------------------------------------------------------
!     Note as diazotrophs have different Qp, we must route
!       enough P into zoopl and sinking detritus pools to fill
!       their fixed p/C ratios.
!     The left over P (remaining_diazP) is split between DOP and DIP pools
!--------------------------------------------------------------------------

      remaining_diazP =  ((graze_diaz + diaz_loss) * Qp_diaz) &
        - ((graze_diaz_poc+graze_diaz_zoo) * Qp)
      diaz_loss_dop = (c1 - parm_labile_ratio) * remaining_diazP
      diaz_loss_dip = parm_labile_ratio * remaining_diazP

!-----------------------------------------------------------------------
!   get fractional factor for routing of zoo losses, based on food supply
!   more material is routed to large detrital pool when diatoms eaten
!-----------------------------------------------------------------------

      f_zoo_detr = (f_diat_zoo_detr * (graze_diat + epsC * epsTinv) &
        + f_sp_zoo_detr * (graze_sp + epsC * epsTinv) &
        + f_diaz_zoo_detr * (graze_diaz + epsC * epsTinv)) &
        / (graze_diat + graze_sp + graze_diaz + c3 &
        * epsC * epsTinv)

!-----------------------------------------------------------------------
!   compute zoo threshold C concentration (mmol C/m^3)
!-----------------------------------------------------------------------

      C_loss_thres = f_loss_thres * loss_thres_zoo

      Zprime = max(zooC_loc - C_loss_thres, c0)

      zoo_loss = z_mort2 * Zprime * Zprime + z_mort * Zprime

      zoo_loss_doc = (c1 - parm_labile_ratio) * (c1 - f_zoo_detr) &
        * zoo_loss
      zoo_loss_dic = parm_labile_ratio * (c1 - f_zoo_detr) * zoo_loss

!-----------------------------------------------------------------------
!   compute terms for DOM
!-----------------------------------------------------------------------
      DOC_prod = sp_loss_doc + graze_sp_doc + zoo_loss_doc &
        + diat_loss_doc + graze_diat_doc + diaz_loss_doc  &
        + graze_diaz_doc

      DON_prod = (DOC_prod * Q) + diaz_Nexcrete
      DOP_prod = (sp_loss_doc + graze_sp_doc + zoo_loss_doc &
        + diat_loss_doc + graze_diat_doc) * Qp + diaz_loss_dop
      DOFe_prod = (zoo_loss_doc * Qfe_zoo) &
        + (Qfe_sp * (graze_sp_doc + sp_loss_doc)) &
        + (Qfe_diat * (graze_diat_doc + diat_loss_doc)) &
        + (Qfe_diaz * (graze_diaz_doc + diaz_loss_doc))

      DOC_remin  = DOC_loc  * DOM_remin
      DON_remin  = DON_loc  * DOM_remin
      DOFe_remin = DOFe_loc * DOM_remin
      DOP_remin  = DOP_loc  * DOM_remin

!-----------------------------------------------------------------------
!     large detritus C
!-----------------------------------------------------------------------

      POC%prod = sp_agg + graze_sp_poc + sp_loss_poc + f_zoo_detr &
        * zoo_loss + diat_loss_poc + diat_agg + graze_diat_poc &
        + graze_diaz_poc

!-----------------------------------------------------------------------
!     large detrital CaCO3
!     33% of CaCO3 is remin when phyto are grazed
!-----------------------------------------------------------------------

      P_CaCO3%prod = ((c1 - f_graze_CaCO3_remin) * graze_sp + sp_loss  &
        + sp_agg) * QCaCO3

!-----------------------------------------------------------------------
!     large detritus SiO2
!     grazed diatom SiO2, 60% is remineralized
!-----------------------------------------------------------------------

      P_SiO2%prod = ((c1 - f_graze_si_remin) * graze_diat + diat_agg &
        + f_diat_loss_poc * diat_loss) * Qsi

      dust%prod = c0

!-----------------------------------------------------------------------
!  Compute iron scavenging :
!  1) compute in terms of loss per year per unit iron (%/year/fe)
!  2) scale by sinking POC/Dust flux
!  3) increase scavenging at higher iron (>0.6nM)
!  4) decrease scavenging rates at low iron
!  5) convert to net loss per second
!-----------------------------------------------------------------------

      Fe_scavenge_rate = fe_scavenge_rate0

      Fe_scavenge_rate = Fe_scavenge_rate  &
        * min(((POC%sflux_out + POC%hflux_out  &
        + ((dust%sflux_out + dust%hflux_out) * dust_fescav_scale)) &
        / parm_POC_flux_ref), fe_max_scale1)

      where (Fe_loc > fe_scavenge_thres1) &
        Fe_scavenge_rate = Fe_scavenge_rate  &
        + (Fe_loc - fe_scavenge_thres1) &
        * fe_max_scale2

      where (Fe_loc < fe_scavenge_thres2) &
        Fe_scavenge_rate = Fe_scavenge_rate  &
        * (Fe_loc / fe_scavenge_thres2)

      Fe_scavenge = yps * Fe_loc * Fe_scavenge_rate

      P_iron%prod = ((sp_agg + graze_sp_poc + sp_loss_poc) * Qfe_sp) &
        + (zoo_loss * f_zoo_detr * Qfe_zoo) &
        + ((diat_agg + graze_diat_poc + diat_loss_poc) * Qfe_diat) &
        + (graze_diaz_poc * Qfe_diaz) + (f_fescav_P_iron * Fe_scavenge)

      !call compute_particulate_terms(k, POC, P_CaCO3, P_SiO2, dust, &
      !  P_iron, QA_dust_def, TEMP, O2_loc, KMT, dz, dzr, zt)

      ! particulate remineralization hacket here to prevent total remin
      ! in bottom cell (i.e. we want particles to fall out, b/c we are
      ! not simulating the entire water column
      call compute_particulate_terms(k, POC, P_CaCO3, P_SiO2, dust, &
        P_iron, QA_dust_def, TEMP, O2_loc, KMT+2, dz, dzr, zt)

!-----------------------------------------------------------------------
!     nitrate & ammonium
!     nitrification in low light
!     use exponential decay of PAR across model level to compute taper factor
!-----------------------------------------------------------------------

      if (lrest_no3) then
         RESTORE = (NO3_CLIM(:,:,k) - NO3_loc) * nutr_rest_time_inv(k)
      else
         RESTORE = c0
      endif

      NO3_RESTORE_tavg = RESTORE

      where (PAR_out < parm_nitrif_par_lim)
         NITRIF = parm_kappa_nitrif * NH4_loc
         where (PAR_in > parm_nitrif_par_lim)
            NITRIF = NITRIF * log(PAR_out / parm_nitrif_par_lim) / &
                (-KPARdz)
         endwhere
      elsewhere
         NITRIF = c0
      endwhere

!-----------------------------------------------------------------------
! Compute denitrification under low O2 conditions
!-----------------------------------------------------------------------

      where ((O2_loc.le.parm_o2_min) .and. (NO3_loc.gt.parm_no3_min))
         DENITRIF = (DOC_remin + POC%remin) / denitrif_C_N
      elsewhere
         DENITRIF = c0
      endwhere

!-----------------------------------------------------------------------
!     nitrate and ammonium
!-----------------------------------------------------------------------

      DTRACER_MODULE(:,:,no3_ind) = RESTORE + NITRIF - NO3_V_diat &
        - NO3_V_sp - DENITRIF

      DTRACER_MODULE(:,:,nh4_ind) = -NH4_V_diat - NH4_V_sp - NITRIF &
        + Q * (zoo_loss_dic + sp_loss_dic + graze_sp_dic + diat_loss_dic &
        + graze_diat_dic + POC%remin + diaz_loss_dic &
        + graze_diaz_dic) + DON_remin

!-----------------------------------------------------------------------
!     dissolved iron
!-----------------------------------------------------------------------

      DTRACER_MODULE(:,:,fe_ind) = P_iron%remin  &
        + (Qfe_zoo * zoo_loss_dic) + DOFe_remin  &
        + (Qfe_sp * (sp_loss_dic + graze_sp_dic)) &
        + (Qfe_diat * (diat_loss_dic + graze_diat_dic)) &
        + (Qfe_diaz * (diaz_loss_dic + graze_diaz_dic)) &
        + graze_diaz_zoo * (Qfe_diaz-Qfe_zoo) &
        + graze_diat_zoo * (Qfe_diat-Qfe_zoo) &
        + graze_sp_zoo * (Qfe_sp-Qfe_zoo) &
        + eco_inject(Fe_ind) &
        - photoFe_sp - photoFe_diat - photoFe_diaz - Fe_scavenge

!-----------------------------------------------------------------------
!     dissolved SiO3
!-----------------------------------------------------------------------

      if (lrest_sio3) then
         RESTORE = (SiO3_CLIM(:,:,k) - SiO3_loc) * nutr_rest_time_inv(k)
      else
         RESTORE = c0
      endif

      SiO3_RESTORE_tavg = RESTORE

      DTRACER_MODULE(:,:,sio3_ind) = RESTORE + P_SiO2%remin &
        + Qsi * (f_graze_si_remin * graze_diat + f_diat_loss_dc  &
        * diat_loss) - photoSi_diat

!-----------------------------------------------------------------------
!     phosphate
!-----------------------------------------------------------------------

      if (lrest_po4) then
         RESTORE = (PO4_CLIM(:,:,k) - PO4_loc) * nutr_rest_time_inv(k)
      else
         RESTORE = c0
      endif

      PO4_RESTORE_tavg = RESTORE

      DTRACER_MODULE(:,:,po4_ind) = RESTORE + (Qp * (POC%remin &
        + zoo_loss_dic + sp_loss_dic + graze_sp_dic + diat_loss_dic &
        + graze_diat_dic - photoC_sp - photoC_diat)) + DOP_remin &
        + diaz_loss_dip - (photoC_diaz * Qp_diaz)

!-----------------------------------------------------------------------
!     small phyto Carbon
!-----------------------------------------------------------------------

      DTRACER_MODULE(:,:,spC_ind) = photoC_sp - graze_sp - sp_loss &
        - sp_agg

!-----------------------------------------------------------------------
!     small phyto Chlorophyll
!-----------------------------------------------------------------------

      DTRACER_MODULE(:,:,spChl_ind) = photoacc_sp &
        - thetaC_sp * (graze_sp + sp_loss + sp_agg)

!-----------------------------------------------------------------------
!     small phytoplankton CaCO3
!-----------------------------------------------------------------------

      DTRACER_MODULE(:,:,spCaCO3_ind) = CaCO3_prod &
        - (graze_sp + sp_loss + sp_agg) * QCaCO3

!-----------------------------------------------------------------------
!     diatom Carbon
!-----------------------------------------------------------------------

      !DTRACER_MODULE(:,:,diatC_ind) = photoC_diat - graze_diat &
      !  - diat_loss - diat_agg
      DTRACER_MODULE(:,:,diatC_ind) = photoC_diat - graze_diat &
        + eco_inject(diatC_ind) &
        - diat_loss - diat_agg

!-----------------------------------------------------------------------
!     diatom Chlorophyll
!-----------------------------------------------------------------------

      DTRACER_MODULE(:,:,diatChl_ind) = photoacc_diat &
        + eco_inject(diatChl_ind) &
        - thetaC_diat * (graze_diat + diat_loss + diat_agg)

!-----------------------------------------------------------------------
!     zoo Carbon
!-----------------------------------------------------------------------

      DTRACER_MODULE(:,:,zooC_ind) = graze_sp_zoo + graze_diat_zoo &
        + graze_diaz_zoo - zoo_loss

!-----------------------------------------------------------------------
!     dissolved organic Matter
!-----------------------------------------------------------------------

      DTRACER_MODULE(:,:,doc_ind)  = DOC_prod  - DOC_remin

      DTRACER_MODULE(:,:,don_ind)  = DON_prod  - DON_remin

      DTRACER_MODULE(:,:,dop_ind)  = DOP_prod  - DOP_remin

      DTRACER_MODULE(:,:,dofe_ind) = DOFe_prod - DOFe_remin

!-----------------------------------------------------------------------
!     small phyto Fe
!-----------------------------------------------------------------------

      DTRACER_MODULE(:,:,spFe_ind) =  photoFe_sp &
        - (Qfe_sp * (graze_sp+sp_loss+sp_agg))

!-----------------------------------------------------------------------
!     Diatom Fe
!-----------------------------------------------------------------------

      DTRACER_MODULE(:,:,diatFe_ind) =  photoFe_diat &
       + eco_inject(diatFe_ind) &
       - (Qfe_diat * (graze_diat+diat_loss+diat_agg))

!-----------------------------------------------------------------------
!     Diatom Si
!-----------------------------------------------------------------------

      DTRACER_MODULE(:,:,diatSi_ind) =  photoSi_diat &
        + eco_inject(diatSi_ind) &
        - (Qsi * (graze_diat+diat_loss+diat_agg))

!-----------------------------------------------------------------------
!     Diazotroph C
!-----------------------------------------------------------------------

      DTRACER_MODULE(:,:,diazC_ind) =  photoC_diaz - graze_diaz &
        - diaz_loss

!-----------------------------------------------------------------------
!     diazotroph Chlorophyll
!-----------------------------------------------------------------------

      DTRACER_MODULE(:,:,diazChl_ind) = photoacc_diaz &
        - thetaC_diaz * (graze_diaz + diaz_loss)

!-----------------------------------------------------------------------
!     Diazotroph Fe
!-----------------------------------------------------------------------

      DTRACER_MODULE(:,:,diazFe_ind) =  photoFe_diaz &
        - (Qfe_diaz * (graze_diaz + diaz_loss))

!-----------------------------------------------------------------------
!     dissolved inorganic Carbon
!-----------------------------------------------------------------------

      DTRACER_MODULE(:,:,dic_ind) = DOC_remin + POC%remin &
        + P_CaCO3%remin + f_graze_CaCO3_remin * graze_sp * QCaCO3 &
        + zoo_loss_dic + sp_loss_dic + graze_sp_dic + diat_loss_dic &
        + graze_diat_dic - photoC_sp - photoC_diat - CaCO3_prod &
        + graze_diaz_dic + diaz_loss_dic - photoC_diaz

!-----------------------------------------------------------------------
!     alkalinity
!-----------------------------------------------------------------------

      DTRACER_MODULE(:,:,alk_ind) = -DTRACER_MODULE(:,:,no3_ind) &
        + DTRACER_MODULE(:,:,nh4_ind) + c2 &
        * (P_CaCO3%remin + f_graze_CaCO3_remin * graze_sp  &
        * QCaCO3-CaCO3_prod)

!-----------------------------------------------------------------------
!     oxygen
!-----------------------------------------------------------------------

      DTRACER_MODULE(:,:,o2_ind) = (photoC_sp + photoC_diat)  &
        / parm_Red_D_C_O2 + photoC_diaz/parm_Red_D_C_O2_diaz

      where (O2_loc > parm_o2_min)
        DTRACER_MODULE(:,:,o2_ind) = DTRACER_MODULE(:,:,o2_ind) &
            + ((- POC%remin - DOC_remin - zoo_loss_dic - sp_loss_dic &
            - graze_sp_dic - diat_loss_dic - graze_diat_dic &
            - graze_diaz_dic - diaz_loss_dic) / parm_Red_P_C_O2)
      endwhere

!-----------------------------------------------------------------------
!     various tavg/history variables
!-----------------------------------------------------------------------

      PAR_avg_tavg          = PAR_avg
      graze_sp_tavg         = graze_sp
      graze_diat_tavg       = graze_diat
      graze_diaz_tavg       = graze_diaz
      graze_tot_tavg        = graze_sp + graze_diat + graze_diaz
      sp_loss_tavg          = sp_loss
      diat_loss_tavg        = diat_loss
      diaz_loss_tavg        = diaz_loss
      zoo_loss_tavg         = zoo_loss
      sp_agg_tavg           = sp_agg
      diat_agg_tavg         = diat_agg
      photoC_sp_tavg        = photoC_sp
      photoC_diat_tavg      = photoC_diat
      photoC_diaz_tavg      = photoC_diaz
      tot_prod_tavg         = photoC_sp + photoC_diat + photoC_diaz
      DOC_prod_tavg         = DOC_prod
      DOC_remin_tavg        = DOC_remin
      DON_prod_tavg         = DON_prod
      DON_remin_tavg        = DON_remin
      DOP_prod_tavg         = DOP_prod
      DOP_remin_tavg        = DOP_remin
      DOFe_prod_tavg        = DOFe_prod
      DOFe_remin_tavg       = DOFe_remin
      Fe_scavenge_tavg      = Fe_scavenge
      Fe_scavenge_rate_tavg = Fe_scavenge_rate
      NITRIF_tavg           = NITRIF
      DENITRIF_tavg         = DENITRIF

!-----------------------------------------------------------------------

      end subroutine ecosys_set_interior

!***********************************************************************

      subroutine init_particulate_terms(POC, P_CaCO3, P_SiO2, dust, &
        P_iron, QA_dust_def)

!-----------------------------------------------------------------------
!     Set incoming fluxes (put into outgoing flux for first level usage).
!     Set dissolution length, production fraction and mass terms.
!-----------------------------------------------------------------------

!      use time_management, only : thour00
!      use forcing_tools, only : update_forcing_data, interpolate_forcing

!-----------------------------------------------------------------------
!     arguments
!-----------------------------------------------------------------------

      type(sinking_particle), intent(out) :: &
        POC,          & ! base units = nmol C
        P_CaCO3,      & ! base units = nmol CaCO3
        P_SiO2,       & ! base units = nmol SiO2
        dust,         & ! base units = g
        P_iron        ! base units = nmol Fe

      real(kind=dbl_kind), dimension(imt,jmt), intent(out) :: &
        QA_dust_def     ! incoming deficit in the QA(dust) POC flux

!-----------------------------------------------------------------------
!     local variables
!-----------------------------------------------------------------------

      REAL(KIND=dbl_kind), DIMENSION(imt,jmt) :: &
        net_dust_in        ! net incoming dust flux

!-----------------------------------------------------------------------
!    based on mineral ballast model from Armstrong et al. 2000
!
!    July 2002, length scale for excess POC and bSI modified by temperature
!    Value given here is at Tref of 30 deg. C, JKM
!
!-----------------------------------------------------------------------

      POC%diss      = 13000.0_dbl_kind  ! diss. length (cm), modified by TEMP
      POC%gamma     = c0                ! not used
      POC%mass      = 12.01_dbl_kind   ! molecular weight of POC
      POC%rho       = c0               ! not used

      P_CaCO3%diss  = 60000.0_dbl_kind ! diss. length (cm)
      P_CaCO3%gamma = 0.55_dbl_kind     ! prod frac -> hard subclass
      P_CaCO3%mass  = 100.09_dbl_kind  ! molecular weight of CaCO3
      P_CaCO3%rho   = 0.07_dbl_kind * P_CaCO3%mass / POC%mass
                                       ! QA mass ratio for CaCO3
                                       ! This ratio is used in ecos_set_interior

      P_SiO2%diss   = 2200.0_dbl_kind  ! diss. length (cm), modified by TEMP
      P_SiO2%gamma  = 0.37_dbl_kind    ! prod frac -> hard subclass
      P_SiO2%mass   = 60.08_dbl_kind   ! molecular weight of SiO2
      P_SiO2%rho    = 0.035_dbl_kind * P_SiO2%mass / POC%mass
                                       ! QA mass ratio for SiO2

      dust%diss     = 60000.0_dbl_kind ! diss. length (cm)
      dust%gamma    = 0.97_dbl_kind    ! prod frac -> hard subclass
      dust%mass     = 1.0e9_dbl_kind   ! base units are already grams
      dust%rho      = 0.07_dbl_kind * dust%mass / POC%mass
                                       ! QA mass ratio for dust

      P_iron%diss   = 60000.0_dbl_kind ! diss. length (cm) - not used
      P_iron%gamma  = c0               ! prod frac -> hard subclass - not used
      P_iron%mass   = c0               ! not used
      P_iron%rho    = c0               ! not used

!-----------------------------------------------------------------------
!     Set incoming fluxes
!-----------------------------------------------------------------------

      P_CaCO3%sflux_out = c0
      P_CaCO3%hflux_out = c0

      P_SiO2%sflux_out = c0
      P_SiO2%hflux_out = c0

!------------------------------------------------------------------------
!     Reduce surface dust flux due to assumed instant surface dissolution
!------------------------------------------------------------------------
!      if (...) then
        net_dust_in = dust_flux * (c1 - parm_fe_bioavail)
        dust%sflux_out = (c1 - dust%gamma) * net_dust_in
        dust%hflux_out = dust%gamma * net_dust_in
!      else
!        net_dust_in    = c0
!        dust%sflux_out = c0
!        dust%hflux_out = c0
!      endif

      P_iron%sflux_out = c0
      P_iron%hflux_out = c0

!------------------------------------------------------------------------
!   Hard POC is QA flux and soft POC is excess POC.
!------------------------------------------------------------------------

      POC%sflux_out = c0
      POC%hflux_out = c0

!------------------------------------------------------------------------
!   Compute initial QA(dust) POC flux deficit.
!------------------------------------------------------------------------

      QA_dust_def = dust%rho * (dust%sflux_out + dust%hflux_out)

!-----------------------------------------------------------------------

      end subroutine init_particulate_terms

!***********************************************************************

      subroutine compute_particulate_terms(k, POC, P_CaCO3, P_SiO2, &
        dust, P_iron, QA_dust_def, TEMP, O2_loc, KMT, dz, dzr, zt)

!-----------------------------------------------------------------------
!   Compute outgoing fluxes and remineralization terms. Assumes that
!   production terms have been set. Incoming fluxes are assumed to be the
!   outgoing fluxes from the previous level.
!
!   It is assumed that there is no production of dust.
!
!   Instantaneous remineralization in the bottom cell is implemented by
!   setting the outgoing flux to zero.
!
!   For POC, the hard subclass is the POC flux qualitatively associated
!   with the ballast flux. The soft subclass is the excess POC flux.
!
!   Remineralization for the non-iron particulate pools is computing
!   by first computing the outgoing flux and then computing the
!   remineralization from conservation, i.e.
!      flux_in - flux_out + prod * dz - remin * dz == 0.
!
!   For iron, remineralization is first computed from POC remineralization
!   and then flux_out is computed from conservation. If the resulting
!   flux_out is negative or should be zero because of the sea floor, the
!   remineralization is adjusted.
!   Note: all the sinking iron is in the P_iron%sflux pool, hflux Fe not
!      explicitly tracked, it is assumed that total iron remin is
!            proportional to total POC remin.
!
!   Based upon Armstrong et al. 2000
!
!   July 2002, added temperature effect on remin length scale of
!       excess POC (all soft POM& Iron) and on SiO2.
!   new variable passed into ballast, Tfunc, main Temperature function
!   computed in ecosystem routine.  scaling factor for dissolution
!   of excess POC, Fe, and Bsi now varies with location (f(temperature)).
!
!   Added diffusive iron flux from sediments at depths < 1100m,
!   based on Johnson et al., 1999, value of 5 umolFe/m2/day,
!       this value too high, using 2 umolFe/m2/day here
!
!   Allow hard fraction of ballast to remin with long length scale 40,000m
!         thus ~ 10% of hard ballast remins over 4000m water column.
!
!   Sinking dust flux is decreased by assumed instant solubility/dissolution
!         at ocean surface from the parm_fe_bioavail.
!
!   Modified to allow different Q10 factors for soft POM and bSI remin,
!   water TEMP is now passed in instead of Tfunc (1/2005, JKM)
!-----------------------------------------------------------------------

      !use constants, only : T0_Kelvin
      !use grid, only : zt, dz, dzr, KMT
      !use global_reductions, only : global_count
      !use shr_sys_mod, only: shr_sys_abort

!-----------------------------------------------------------------------
!     arguments
!-----------------------------------------------------------------------

      integer(kind=int_kind), intent(in) :: k ! vertical model level
      type(sinking_particle), intent(inout) :: &
           POC,          & ! base units = nmol C
           P_CaCO3,      & ! base units = nmol CaCO3
           P_SiO2,       & ! base units = nmol SiO2
           dust,         & ! base units = g
           P_iron        ! base units = nmol Fe

      real(kind=dbl_kind), dimension(imt,jmt), intent(inout) :: &
           QA_dust_def     ! incoming deficit in the QA(dust) POC flux

      real(kind=dbl_kind), dimension(imt,jmt), intent(in) :: &
           TEMP          ! temperature for scaling functions

      ! used to change/increase POC%diss
      real(kind=dbl_kind), dimension(imt,jmt), intent(in) :: &
           O2_loc        ! dissolved oxygen

      integer(kind=int_kind) :: kmt ! maximum valid depth layer

      real(kind=dbl_kind), dimension(km) :: &
          dz,        & ! layer thickness
          dzr,       & ! inverse later thickness
          zt          ! layer midpoint position (negative)

!-----------------------------------------------------------------------
!     local variables
!-----------------------------------------------------------------------

      real(kind=dbl_kind) :: poc_diss ! diss. length used (cm)

      character(len=*), parameter :: &
           sub_name = 'ecosys_mod:compute_particulate_terms'

      real(kind=dbl_kind) :: &
           decay_CaCO3,      & ! scaling factor for dissolution of CaCO3
           decay_dust,       & ! scaling factor for dissolution of dust
           decay_POC_E,      & ! scaling factor for dissolution of excess POC
           decay_SiO2,       & ! scaling factor for dissolution of SiO2
           decay_Hard,       & ! scaling factor for dissolution of Hard Ballast
           POC_prod_avail,   & ! POC production available for excess POC flux
           new_QA_dust_def   ! outgoing deficit in the QA(dust) POC flux

      integer(kind=int_kind) :: &
           i, j                  ! loop indices

      real(kind=dbl_kind), dimension(imt,jmt) :: &
           TfuncP,        & ! temperature scaling from soft POM remin
           TfuncS         ! temperature scaling from soft POM remin

!-----------------------------------------------------------------------
!   incoming fluxes are outgoing fluxes from previous level
!-----------------------------------------------------------------------

      P_CaCO3%sflux_in = P_CaCO3%sflux_out
      P_CaCO3%hflux_in = P_CaCO3%hflux_out

      P_SiO2%sflux_in = P_SiO2%sflux_out
      P_SiO2%hflux_in = P_SiO2%hflux_out

      dust%sflux_in = dust%sflux_out
      dust%hflux_in = dust%hflux_out

      POC%sflux_in = POC%sflux_out
      POC%hflux_in = POC%hflux_out

      P_iron%sflux_in = P_iron%sflux_out
      P_iron%hflux_in = P_iron%hflux_out

!-----------------------------------------------------------------------
!     compute decay factors
!-----------------------------------------------------------------------

      decay_CaCO3 = exp(-dz(k) / P_CaCO3%diss)
      decay_dust  = exp(-dz(k) / dust%diss)
      decay_Hard  = exp(-dz(k) / 4.0e6_dbl_kind)

!----------------------------------------------------------------------
!   Tref = 30.0 reference temperature (deg. C)
!
!   Using q10 formulation with Q10 value of 1.12 soft POM (TfuncP) and
!       a Q10 value of 3.5 soft bSi (TfuncS)
!
!-----------------------------------------------------------------------

      TfuncP = 1.12_dbl_kind**(((TEMP + T0_Kelvin) &
        - (Tref + T0_Kelvin)) / c10)

      TfuncS = 4.0_dbl_kind**(((TEMP + T0_Kelvin) &
        - (Tref + T0_Kelvin)) / c10)

      do j = 1,jmt
         do i = 1,imt

            ! foo1:if (LAND_MASK(i,j) .and. k <= KMT(i,j)) then
            foo1:if (LAND_MASK(i,j) .and. k <= KMT) then

            !----------------------------------------------------------------
            ! decay_POC_E and decay_SiO2 set locally, modified by Tfunc
            !----------------------------------------------------------------

                ! increase POC diss where there is denitrification
                if (O2_loc(i,j).lt.parm_o2_min) then
                    poc_diss = 26000.0_dbl_kind  ! diss. length (cm)
                else
                    poc_diss = POC%diss
                endif

                decay_POC_E = exp(-dz(k) / (poc_diss / TfuncP(i,j)))
                decay_SiO2  = exp(-dz(k) / (P_SiO2%diss / TfuncS(i,j)))

            !-----------------------------------------------------------------
            !   Set outgoing fluxes for non-iron pools.
            !   The outoing fluxes for ballast materials are from the
            !   solution of the coresponding continuous ODE across the model
            !   level. The ODE has a constant source term and linear decay.
            !   It is assumed that there is no sub-surface dust production.
            !-----------------------------------------------------------------

                P_CaCO3%sflux_out(i,j) = P_CaCO3%sflux_in(i,j) &
                    * decay_CaCO3 + P_CaCO3%prod(i,j) &
                    * ((c1 - P_CaCO3%gamma) * (c1 - decay_CaCO3) &
                    * P_CaCO3%diss)

                P_CaCO3%hflux_out(i,j) = P_CaCO3%hflux_in(i,j) &
                    * decay_Hard + P_CaCO3%prod(i,j) &
                    * (P_CaCO3%gamma * dz(k))

                P_SiO2%sflux_out(i,j) = P_SiO2%sflux_in(i,j) &
                    * decay_SiO2 + P_SiO2%prod(i,j) &
                    * ((c1 - P_SiO2%gamma) * (c1 - decay_SiO2) &
                    * (P_SiO2%diss / TfuncS(i,j)))

                P_SiO2%hflux_out(i,j) = P_SiO2%hflux_in(i,j) &
                    * decay_Hard + P_SiO2%prod(i,j) &
                    * (P_SiO2%gamma * dz(k))

                dust%sflux_out(i,j) = dust%sflux_in(i,j) * decay_dust

                dust%hflux_out(i,j) = dust%hflux_in(i,j) * decay_Hard

            !-----------------------------------------------------------------
            !   Compute how much POC_prod is available for deficit reduction
            !   and excess POC flux after subtracting off fraction of non-dust
            !   ballast production from net POC_prod.
            !-----------------------------------------------------------------

                 POC_prod_avail = POC%prod(i,j) - P_CaCO3%rho &
                     * P_CaCO3%prod(i,j) - P_SiO2%rho * P_SiO2%prod(i,j)

            !-----------------------------------------------------------------
            !   Check for POC production bounds violations
            !-----------------------------------------------------------------

                 foo2: if (POC_prod_avail < c0) then
                     print *, sub_name,  &
                         ': mass ratio of ballast ', &
                         'production exceeds POC production'
                 endif foo2

            !-----------------------------------------------------------------
            !   Compute 1st approximation to new QA_dust_def, the QA_dust
            !   deficit leaving the cell. Ignore POC_prod_avail at this stage.
            !-----------------------------------------------------------------

                 foo3: if (QA_dust_def(i,j) > 0) then
                     new_QA_dust_def = QA_dust_def(i,j)   &
                     * (dust%sflux_out(i,j) + dust%hflux_out(i,j)) &
                     / (dust%sflux_in(i,j) + dust%hflux_in(i,j))
                 else
                     new_QA_dust_def = c0
                 endif foo3

            !-----------------------------------------------------------------
            !   Use POC_prod_avail to reduce new_QA_dust_def.
            !-----------------------------------------------------------------

                 foo4: if (new_QA_dust_def > c0) then
                    new_QA_dust_def = new_QA_dust_def - POC_prod_avail &
                        * dz(k)
                    if (new_QA_dust_def < c0) then
                       POC_prod_avail = -new_QA_dust_def / dz(k)
                       new_QA_dust_def = c0
                    else
                       POC_prod_avail = c0
                    endif
                 endif foo4

                 QA_dust_def(i,j) = new_QA_dust_def

            !-----------------------------------------------------------------
            !   Compute outgoing POC fluxes. QA POC flux is computing using
            !   ballast fluxes and new_QA_dust_def. If no QA POC flux came in
            !   and no production occured, then no QA POC flux goes out. This
            !   shortcut is present to avoid roundoff cancellation errors from
            !   the dust%rho * dust_flux_out - QA_dust_def computation.
            !   Any POC_prod_avail still remaining goes into excess POC flux.
            !-----------------------------------------------------------------

                 foo5: if (POC%hflux_in(i,j) == c0 .and.  &
                    POC%prod(i,j) == c0) then
                    POC%hflux_out(i,j) = c0
                 else
                    POC%hflux_out(i,j) = P_CaCO3%rho &
                         * (P_CaCO3%sflux_out(i,j) &
                         + P_CaCO3%hflux_out(i,j)) &
                         + P_SiO2%rho * (P_SiO2%sflux_out(i,j) &
                         + P_SiO2%hflux_out(i,j)) + dust%rho &
                         * (dust%sflux_out(i,j) + dust%hflux_out(i,j)) &
                         - new_QA_dust_def
                    POC%hflux_out(i,j) = max(POC%hflux_out(i,j), c0)
                 endif foo5

                 POC%sflux_out(i,j) = POC%sflux_in(i,j) * decay_POC_E &
                     + POC_prod_avail *((c1 - decay_POC_E) &
                     * (poc_diss / TfuncP(i,j)))

            !-----------------------------------------------------------------
            !   Compute remineralization terms. It is assumed that there is no
            !   sub-surface dust production.
            !-----------------------------------------------------------------

                 P_CaCO3%remin(i,j) = P_CaCO3%prod(i,j) &
                     + ((P_CaCO3%sflux_in(i,j) - P_CaCO3%sflux_out(i,j)) &
                     + (P_CaCO3%hflux_in(i,j) - P_CaCO3%hflux_out(i,j))) &
                     * dzr(k)

                 P_SiO2%remin(i,j) = P_SiO2%prod(i,j) &
                     + ((P_SiO2%sflux_in(i,j) - P_SiO2%sflux_out(i,j)) &
                     + (P_SiO2%hflux_in(i,j) - P_SiO2%hflux_out(i,j))) &
                     * dzr(k)

                 POC%remin(i,j) = POC%prod(i,j) + ((POC%sflux_in(i,j) &
                     - POC%sflux_out(i,j)) + (POC%hflux_in(i,j) &
                     - POC%hflux_out(i,j))) * dzr(k)

                 dust%remin(i,j) = &
                     ((dust%sflux_in(i,j) - dust%sflux_out(i,j)) &
                     + (dust%hflux_in(i,j) - dust%hflux_out(i,j))) &
                     * dzr(k)

            !-----------------------------------------------------------------
            !   Compute iron remineralization and flux out.
            !-----------------------------------------------------------------

                 foo6: if (POC%sflux_in(i,j) + POC%hflux_in(i,j) == c0)  &
                       then
                    P_iron%remin(i,j) = (POC%remin(i,j) * parm_Red_Fe_C)
                 else
                    P_iron%remin(i,j) = (POC%remin(i,j) &
                        * (P_iron%sflux_in(i,j) + P_iron%hflux_in(i,j)) &
                        / (POC%sflux_in(i,j) + POC%hflux_in(i,j)))
                 endif foo6

                 P_iron%sflux_out(i,j) = P_iron%sflux_in(i,j) + dz(k) &
                     * ((c1 - P_iron%gamma) * P_iron%prod(i,j) &
                     - P_iron%remin(i,j))

                 foo7: if (P_iron%sflux_out(i,j) < c0) then
                    P_iron%sflux_out(i,j) = c0
                    P_iron%remin(i,j) = P_iron%sflux_in(i,j) * dzr(k) &
                        + (c1 - P_iron%gamma) * P_iron%prod(i,j)
                 endif foo7

            !----------------------------------------------------------------
            !  Compute iron release from dust remin/dissolution
            !
            !   dust remin gDust = 0.035 / 55.847 * 1.0e9 = 626712.0 nmolFe
            !                       gFe     molFe     nmolFe
            !----------------------------------------------------------------

                 P_iron%remin(i,j) = P_iron%remin(i,j) + dust%remin(i,j) &
                     * dust_to_Fe

            else
                 P_CaCO3%sflux_out(i,j) = c0
                 P_CaCO3%hflux_out(i,j) = c0
                 P_CaCO3%remin(i,j) = c0

                 P_SiO2%sflux_out(i,j) = c0
                 P_SiO2%hflux_out(i,j) = c0
                 P_SiO2%remin(i,j) = c0

                 dust%sflux_out(i,j) = c0
                 dust%hflux_out(i,j) = c0
                 dust%remin(i,j) = c0

                 POC%sflux_out(i,j) = c0
                 POC%hflux_out(i,j) = c0
                 POC%remin(i,j) = c0

                 P_iron%sflux_out(i,j) = c0
                 P_iron%hflux_out(i,j) = c0
                 P_iron%remin(i,j) = c0
            endif foo1

          !---------------------------------------------------------------------
          !   Remineralize everything in bottom cell.
          !---------------------------------------------------------------------

            ! bar: if (LAND_MASK(i,j) .and. k == KMT(i,j)) then
            bar: if (LAND_MASK(i,j) .and. k == KMT) then
               P_CaCO3%remin(i,j) = P_CaCO3%remin(i,j) &
                   + (P_CaCO3%sflux_out(i,j) + P_CaCO3%hflux_out(i,j)) &
                   * dzr(k)
               P_CaCO3%sflux_out(i,j) = c0
               P_CaCO3%hflux_out(i,j) = c0

               P_SiO2%remin(i,j) = P_SiO2%remin(i,j) &
                   + (P_SiO2%sflux_out(i,j) + P_SiO2%hflux_out(i,j)) &
                   * dzr(k)
               P_SiO2%sflux_out(i,j) = c0
               P_SiO2%hflux_out(i,j) = c0

               dust%remin(i,j) = dust%remin(i,j) + (dust%sflux_out(i,j) &
                   + dust%hflux_out(i,j)) * dzr(k)
               dust%sflux_out(i,j) = c0
               dust%hflux_out(i,j) = c0

               POC%remin(i,j) = POC%remin(i,j) + (POC%sflux_out(i,j) &
                   + POC%hflux_out(i,j)) * dzr(k)
               POC%sflux_out(i,j) = c0
               POC%hflux_out(i,j) = c0

               P_iron%remin(i,j) = P_iron%remin(i,j) &
                   + (P_iron%sflux_out(i,j) + P_iron%hflux_out(i,j)) &
                   * dzr(k)
               P_iron%sflux_out(i,j) = c0
               P_iron%hflux_out(i,j) = c0

          !-------------------------------------------------------------------
          ! Add diffusive iron flux if depth < 1100.0m, based on
          ! Johnson et al.1999, value of 5.0 umolFe/m2/day.
          !
          ! eco2.42 - eco2.43, 0.5 umolFe/m2/day, 5.78704e-7 nmolFe/cm2/sec.
          ! begin run eco2.45, 1.0 umolFe/m2/day,  1.1574e-6 nmolFe/cm2/sec.
          ! begin run eco2.06, 2.0 umolFe/m2/day,  2.3148e-6 nmolFe/cm2/sec
          !-------------------------------------------------------------------

               if (zt(k) < thres_fe) then
                   P_iron%remin(i,j) = P_iron%remin(i,j) &
                       + (fe_diff_rate * dzr(k))
               endif

            endif bar

         enddo
      enddo

!-----------------------------------------------------------------------
!     Set history variables.
!-----------------------------------------------------------------------

!      POC_FLUX_IN_tavg     = POC%sflux_in + POC%hflux_in
!      POC_PROD_tavg        = POC%prod
!      POC_REMIN_tavg       = POC%remin
!      CaCO3_FLUX_IN_tavg   = P_CaCO3%sflux_in + P_CaCO3%hflux_in
!      CaCO3_PROD_tavg      = P_CaCO3%prod
!      CaCO3_REMIN_tavg     = P_CaCO3%remin
!      SiO2_FLUX_IN_tavg    = P_SiO2%sflux_in + P_SiO2%hflux_in
!      SiO2_PROD_tavg       = P_SiO2%prod
!      SiO2_REMIN_tavg      = P_SiO2%remin
!      dust_FLUX_IN_tavg    = dust%sflux_in + dust%hflux_in
!      dust_REMIN_tavg      = dust%remin
!      P_iron_FLUX_IN_tavg  = P_iron%sflux_in + P_iron%hflux_in
!      P_iron_PROD_tavg     = P_iron%prod
!      P_iron_REMIN_tavg    = P_iron%remin

!-----------------------------------------------------------------------

      end subroutine compute_particulate_terms

!***********************************************************************

      subroutine ecosys_set_sflux(IFRAC,ATM_PRESS,U10_SQR,SST,SSS, &
        SURF_VALS,STF_MODULE)

!-----------------------------------------------------------------------
!     compute surface fluxes for ecosystem BGC tracers
!-----------------------------------------------------------------------

      !use time_management, only : thour00, nsteps_run, check_time_flag
      !use forcing_tools, only : update_forcing_data, interpolate_forcing
      !use constants, only : rho_sw, c1000, p5
      !use co2calc, only : co2calc_row
      !use global_reductions, only : broadcast_scalar
      !use registry, only : named_field_get, named_field_set

!-----------------------------------------------------------------------
!     argument declarations
!-----------------------------------------------------------------------

      real (kind=dbl_kind), dimension(imt,jmt), intent(in) :: &
        IFRAC          & ! sea ice fraction (non-dimensional)
      , ATM_PRESS      & ! sea level atmospheric pressure (dyne/cm^2)
      , U10_SQR        & ! 10m wind speed squared (cm^2/s^2)
      , SST            & ! sea surface temperature (C)
      , SSS            ! sea surface salinity (psu)

      real (kind=dbl_kind), dimension(imt,jmt,ecosys_tracer_cnt), &
        intent(in) :: SURF_VALS ! module tracers (fmol/cm^3)

      real (kind=dbl_kind), dimension(imt,jmt,ecosys_tracer_cnt), &
        intent(inout) :: STF_MODULE ! surface fluxes (fmol/cm^2/s)

!-----------------------------------------------------------------------
!     local parameters
!-----------------------------------------------------------------------

      real(kind=dbl_kind), parameter :: &
          a         = 6.97e-9_dbl_kind, & ! in s/cm, from a = 0.251 cm/hr s^2/m^2 in Wannikhof 2014
          phlo_init = 5.0_dbl_kind,   & ! low bound for ph for no prev soln
          phhi_init = 9.0_dbl_kind,   & ! high bound for ph for no prev soln
          del_ph    = 0.25_dbl_kind   ! delta-ph for prev soln
          !a         = 8.6e-9_dbl_kind,& ! a = 0.31 cm/hr s^2/m^2 in (s/cm)    -- older value
!-----------------------------------------------------------------------
!     local variable declarations
!-----------------------------------------------------------------------

      real (kind=dbl_kind), dimension(imt,jmt) :: &
        IFRAC_USED      & ! used ice fraction (m^2/m^2)
      , AP_USED         & ! used atmospheric pressure (atm)
      , XKW             & ! a * U10_SQR (cm/s)
      , XKW_ICE         & ! common portion of piston vel., a*(1-fice1)*u**2 (cm/s)
      , SCHMIDT_USED    & ! Schmidt no. of O2 or CO2 (non-dimensional)
      , PV              & ! piston velocity of O2 or CO2 (cm/s)
      , O2SAT_1atm      & ! O2 saturation @ 1 atm (mmol/m^3)
      , O2SAT_USED      & ! used O2 saturation (mmol/m^3)
      , FLUX            & ! gas flux of O2 or CO2 (nmol/cm^2/s)
      , XCO2            ! atmospheric CO2 (2D)

      real (kind=dbl_kind), dimension(imt) :: &
        PHLO_row        & ! lower bound for pH solver
      , PHHI_row        & ! upper bound for pH solver
      , PH_NEW          & ! computed PH from solver
      , CO2STAR_ROW     & ! CO2STAR from solver
      , DCO2STAR_ROW    & ! DCO2STAR from solver
      , pCO2SURF_ROW    & ! pCO2SURF from solver
      , DpCO2_ROW       ! DpCO2 from solver

      logical (kind=log_kind), dimension(imt) :: &
        MASK_row        ! mask for pH solver

      integer (kind=int_kind) :: &
        j               ! 'latitudinal' index


!-----------------------------------------------------------------------

      STF_MODULE = c0

!-----------------------------------------------------------------------

      IFRAC_USED = IFRAC

      XKW = a * U10_SQR

!-----------------------------------------------------------------------
!        convert ATM_PRESS to atm
!-----------------------------------------------------------------------

      AP_USED = ATM_PRESS / 1013.25e+3_dbl_kind ! (dyne/cm^2 -> atm)

      XKW_tavg = XKW
      AP_tavg  = AP_USED

!-----------------------------------------------------------------------
!     Compute XKW_ICE. XKW is zero over land, so XKW_ICE is too.
!-----------------------------------------------------------------------

      if (lflux_gas_o2 .or. lflux_gas_co2) then
         XKW_ICE = XKW
         where (IFRAC_USED > c0 .and. IFRAC_USED < c1)
            XKW_ICE = (c1 - IFRAC_USED) * XKW_ICE
         endwhere
         where (IFRAC_USED >= c1)
            XKW_ICE = c0
         endwhere
      endif

!-----------------------------------------------------------------------
!     compute O2 flux
!-----------------------------------------------------------------------

      if (lflux_gas_o2) then
         SCHMIDT_USED = SCHMIDT_O2(SST)
         SCHMIDT_O2_tavg = SCHMIDT_USED
         O2SAT_1atm = O2SAT(SST,SSS)
         where (LAND_MASK)
            PV = XKW_ICE * sqrt(660.0_dbl_kind / SCHMIDT_USED)
            O2SAT_USED = AP_USED * O2SAT_1atm
            O2SAT_tavg = O2SAT_USED
            FLUX = PV * (O2SAT_USED - SURF_VALS(:,:,o2_ind))
            STF_MODULE(:,:,o2_ind) = FLUX
            FG_O2_tavg = FLUX
            PV_O2_tavg = PV
         else where
            O2SAT_tavg = c0
            FG_O2_tavg = c0
            PV_O2_tavg = c0
         endwhere
      endif

!-----------------------------------------------------------------------
!     Set XCO2
!-----------------------------------------------------------------------

!      select case (atm_co2_iopt)
!      case (atm_co2_iopt_const)
        XCO2 = atm_co2_const
!      case (atm_co2_iopt_model)
!        call named_field_get(atm_co2_nf_ind, XCO2)
!      end select

!-----------------------------------------------------------------------
!     compute CO2 flux, computing disequilibrium one row at a time
!-----------------------------------------------------------------------

      if (lflux_gas_co2) then
         SCHMIDT_USED = SCHMIDT_CO2(SST)
         SCHMIDT_CO2_tavg = SCHMIDT_USED
         where (LAND_MASK)
            PV = XKW_ICE * sqrt(660.0_dbl_kind / SCHMIDT_USED)
         else where
            PV = c0
         endwhere

         PV_CO2_tavg = PV

         do j = 1,jmt
            !where (PH_PREV(:,j) /= c0)
            !   PHLO_row = PH_PREV(:,j) - del_ph
            !   PHHI_row = PH_PREV(:,j) + del_ph
            !else where
               PHLO_row = 4. !phlo_init
               PHHI_row = 10. !phhi_init
            !endwhere


            call co2calc_row(LAND_MASK(:,j), SST(:,j), SSS(:,j), &
                 SURF_VALS(:,j,dic_ind), SURF_VALS(:,j,alk_ind), &
                 SURF_VALS(:,j,po4_ind), SURF_VALS(:,j,sio3_ind), &
                 PHLO_row, PHHI_row, PH_NEW, XCO2(:,j), &
                 AP_USED(:,j), CO2STAR_ROW, &
                 DCO2STAR_ROW, pCO2SURF_ROW, DpCO2_ROW)

            PH_PREV(:,j) = PH_NEW
            PH_tavg(:,j) = PH_NEW

            CO2STAR_tavg(:,j)  = CO2STAR_ROW
            DCO2STAR_tavg(:,j) = DCO2STAR_ROW
            pCO2SURF_tavg(:,j) = pCO2SURF_ROW
            DpCO2_tavg(:,j)    = DpCO2_ROW

            FLUX(:,j) = PV(:,j) * DCO2STAR_ROW
         enddo

         STF_MODULE(:,:,dic_ind) = STF_MODULE(:,:,dic_ind) + FLUX
         FG_CO2_tavg = FLUX

!-----------------------------------------------------------------------
!        set air-sea co2 gas flux named field, converting units from
!        nmol/cm^2/s (positive down) to kg CO2/m^2/s (positive down)
!-----------------------------------------------------------------------

         !call named_field_set(sflux_co2_nf_ind, 44.0e-8_dbl_kind * FLUX)

      endif

!      if (iron_flux%has_data) then
         FLUX = 0.
!      else
!         FLUX = c0
!      endif

      FLUX = FLUX * parm_Fe_bioavail

      STF_MODULE(:,:,fe_ind) = FLUX
      IRON_FLUX_tavg = FLUX


!-----------------------------------------------------------------------

      end subroutine ecosys_set_sflux

  !***********************************************************************

  SUBROUTINE co2calc_row(mask, t, s, dic_in, ta_in, pt_in, sit_in, &
       phlo, phhi, ph, xco2_in, atmpres, co2star, dco2star, pCO2surf, dpco2)

    !---------------------------------------------------------------------------
    !   SUBROUTINE CO2CALC
    !
    !   PURPOSE : Calculate delta co2*, etc. from total alkalinity, total CO2,
    !             temp, salinity (s), etc.
    !---------------------------------------------------------------------------

    ! USE constants, ONLY : c0, c1, c10, c1000, rho_sw, T0_Kelvin

    !---------------------------------------------------------------------------
    !   input arguments
    !---------------------------------------------------------------------------

    LOGICAL(KIND=log_kind), DIMENSION(imt), INTENT(IN) :: mask
    REAL(KIND=dbl_kind), DIMENSION(imt), INTENT(IN) :: &
         t,        & ! temperature (degrees C)
         s,        & ! salinity (PSU)
         dic_in,   & ! total inorganic carbon (mmol/m^3)
         ta_in,    & ! total alkalinity (meq/m^3)
         pt_in,    & ! inorganic phosphate (mmol/m^3)
         sit_in,   & ! inorganic silicate (mmol/m^3)
         phlo,     & ! lower limit of pH range
         phhi,     & ! upper limit of pH range
         xco2_in,  & ! atmospheric mole fraction CO2 in dry air (ppmv)
         atmpres     ! atmospheric pressure (atmosphere)

    !---------------------------------------------------------------------------
    !   output arguments
    !---------------------------------------------------------------------------

    REAL(KIND=dbl_kind), DIMENSION(imt), INTENT(OUT) :: &
         ph,       & ! computed ph values, for initial guess on next time step
         co2star,  & ! CO2*water (mmol/m^3)
         dco2star, & ! delta CO2 (mmol/m^3)
         pco2surf, & ! oceanic pCO2 (ppmv)
         dpco2       ! Delta pCO2, i.e, pCO2ocn - pCO2atm (ppmv)

    !---------------------------------------------------------------------------
    !   local variable declarations
    !---------------------------------------------------------------------------

    INTEGER(KIND=int_kind) :: i

    REAL(KIND=dbl_kind) :: &
         mass_to_vol,  & ! (mol/kg) -> (mmol/m^3)
         vol_to_mass,  & ! (mmol/m^3) -> (mol/kg)
         tk,           & ! temperature (K)
         is,           & ! ionic strength
         scl,          & ! chlorinity
         co2starair,   & ! co2star saturation
         tk100, tk1002, invtk, dlogtk, is2, sqrtis, &
         s2, sqrts, s15, htotal2

    REAL(KIND=dbl_kind), DIMENSION(imt) :: &
         xco2,         & ! atmospheric CO2 (atm)
         htotal,       & ! free concentration of H ion
         x1, x2          ! bounds on htotal for solver

    !---------------------------------------------------------------------------
    !   check for existence of ocean points
    !---------------------------------------------------------------------------

    IF (COUNT(mask) == 0) THEN
       ph          = c0
       co2star     = c0
       dco2star    = c0
       pCO2surf    = c0
       dpCO2       = c0
       RETURN
    END IF

    !---------------------------------------------------------------------------
    !   set unit conversion factors
    !---------------------------------------------------------------------------

    ! changed for use with tracers in mmol/m^3
    !mass_to_vol = 1.0e3_dbl_kind * rho_sw
    mass_to_vol = 1.0e6_dbl_kind * rho_sw
    vol_to_mass = c1 / mass_to_vol

    !---------------------------------------------------------------------------
    !   convert tracer units to per mass & xco2 from uatm to atm
    !---------------------------------------------------------------------------

    DO i = 1,imt
       IF (mask(i)) THEN
          dic(i)  = dic_in(i)  * vol_to_mass
          ta(i)   = ta_in(i)   * vol_to_mass
          pt(i)   = pt_in(i)   * vol_to_mass
          sit(i)  = sit_in(i)  * vol_to_mass
          xco2(i) = xco2_in(i) * 1e-6_dbl_kind

          !---------------------------------------------------------------------
          !   Calculate all constants needed to convert between various
          !   measured carbon species. References for each equation are
          !   noted in the code.  Once calculated, the constants are stored
          !   and passed in the common block "const". The original version
          !   of this code was based on the code by Dickson in Version 2 of
          !   "Handbook of Methods for the Analysis of the Various Parameters
          !   of the Carbon Dioxide System in Seawater", DOE, 1994 (SOP No. 3,
          !   p25-26).
          !   Derive simple terms used more than once
          !---------------------------------------------------------------------

          tk       = T0_Kelvin + t(i)
          tk100    = tk * 1e-2_dbl_kind
          tk1002   = tk100 * tk100
          invtk    = c1 / tk
          dlogtk   = LOG(tk)

          is       = 19.924 * s(i) / (c1000 - 1.005 * s(i))
          is2      = is * is
          sqrtis   = SQRT(is)
          sqrts    = SQRT(s(i))
          s15      = s(i) ** 1.5
          s2       = s(i) ** 2
          scl      = s(i) / 1.80655

          !---------------------------------------------------------------------
          !   f = k0(1-pH2O)*correction term for non-ideality
          !   Weiss & Price (1980, Mar. Chem., 8, 347-359;
          !                 Eq 13 with table 6 values)
          !---------------------------------------------------------------------

          ff(i) = EXP(-162.8301 + 218.2968/tk100 + 90.9241*LOG(tk100) - &
               1.47696*tk1002 + s(i)*(.025695 - .025225*tk100 + &
               0.0049867*tk1002))

          !---------------------------------------------------------------------
          !   K0 from Weiss 1974
          !---------------------------------------------------------------------

          k0(i) = EXP(93.4517/tk100 - 60.2409 + 23.3585*LOG(tk100) + &
               s(i)*(.023517 - 0.023656 * tk100 + 0.0047036 * tk1002))

          !---------------------------------------------------------------------
          !   k1 = [H][HCO3]/[H2CO3]
          !   k2 = [H][CO3]/[HCO3]
          !   Millero p.664 (1995) using Mehrbach et al. data on seawater scale
          !---------------------------------------------------------------------

          k1(i) = 10**(-1*(3670.7*invtk - 62.008 + 9.7944*dlogtk - &
               0.0118*s(i) + 0.000116*s2))

          k2(i) = 10**(-1*(1394.7*invtk + 4.777 - 0.0184*s(i) + 0.000118*s2))

          !---------------------------------------------------------------------
          !   kb = [H][BO2]/[HBO2]
          !   Millero p.669 (1995) using data from Dickson (1990)
          !---------------------------------------------------------------------

          kb(i) = EXP((-8966.90 - 2890.53*sqrts - 77.942*s(i) + &
               1.728*s15 - 0.0996*s2)*invtk + &
               (148.0248 + 137.1942*sqrts + 1.62142*s(i)) + &
               (-24.4344 - 25.085*sqrts - 0.2474*s(i)) * &
               dlogtk + 0.053105*sqrts*tk)

          !---------------------------------------------------------------------
          !   k1p = [H][H2PO4]/[H3PO4]
          !   DOE(1994) eq 7.2.20 with footnote using data from Millero (1974)
          !---------------------------------------------------------------------

          k1p(i) = EXP(-4576.752*invtk + 115.525 - 18.453 * dlogtk + &
               (-106.736*invtk + 0.69171) * sqrts + &
               (-0.65643*invtk - 0.01844) * s(i))

          !---------------------------------------------------------------------
          !   k2p = [H][HPO4]/[H2PO4]
          !   DOE(1994) eq 7.2.23 with footnote using data from Millero (1974))
          !---------------------------------------------------------------------

          k2p(i) = EXP(-8814.715*invtk + 172.0883 - 27.927 * dlogtk + &
               (-160.340*invtk + 1.3566) * sqrts + &
               (0.37335*invtk - 0.05778) * s(i))

          !---------------------------------------------------------------------
          !   k3p = [H][PO4]/[HPO4]
          !   DOE(1994) eq 7.2.26 with footnote using data from Millero (1974)
          !---------------------------------------------------------------------

          k3p(i) = EXP(-3070.75*invtk - 18.141 +  &
               (17.27039*invtk + 2.81197) * sqrts + &
               (-44.99486*invtk - 0.09984) * s(i))

          !---------------------------------------------------------------------
          !   ksi = [H][SiO(OH)3]/[Si(OH)4]
          !   Millero p.671 (1995) using data from Yao and Millero (1995)
          !---------------------------------------------------------------------

          ksi(i) = EXP(-8904.2*invtk + 117.385 - 19.334 * dlogtk + &
               (-458.79*invtk + 3.5913) * sqrtis + &
               (188.74*invtk - 1.5998) * is + &
               (-12.1652*invtk + 0.07871) * is2 + &
               LOG(1.0-0.001005*s(i)))

          !---------------------------------------------------------------------
          !   kw = [H][OH]
          !   Millero p.670 (1995) using composite data
          !---------------------------------------------------------------------

          kw(i) = EXP(-13847.26*invtk + 148.9652 - 23.6521 * dlogtk + &
               (118.67*invtk - 5.977 + 1.0495 * dlogtk) * &
               sqrts - 0.01615 * s(i))

          !---------------------------------------------------------------------
          !   ks = [H][SO4]/[HSO4]
          !   Dickson (1990, J. chem. Thermodynamics 22, 113)
          !---------------------------------------------------------------------

          ks(i) = EXP(-4276.1*invtk + 141.328 - 23.093*dlogtk + &
               (-13856*invtk + 324.57 - 47.986*dlogtk) * &
               sqrtis + &
               (35474*invtk - 771.54 + 114.723*dlogtk) * is - &
               2698*invtk*is**1.5 + 1776*invtk*is2 + &
               LOG(1.0 - 0.001005*s(i)))

          !---------------------------------------------------------------------
          !   kf = [H][F]/[HF]
          !   Dickson and Riley (1979) -- change pH scale to total
          !---------------------------------------------------------------------

          kf(i) = EXP(1590.2*invtk - 12.641 + 1.525*sqrtis + &
               LOG(1.0 - 0.001005*s(i)) +  &
               LOG(1.0 + (0.1400/96.062)*(scl)/ks(i)))

          !---------------------------------------------------------------------
          !   Calculate concentrations for borate, sulfate, and fluoride
          !   bt : Uppstrom (1974)
          !   st : Morris & Riley (1966)
          !   ft : Riley (1965)
          !---------------------------------------------------------------------

          bt(i) = 0.000232 * scl/10.811
          st(i) = 0.14 * scl/96.062
          ft(i) = 0.000067 * scl/18.9984

          x1(i) = c10 ** (-phhi(i))
          x2(i) = c10 ** (-phlo(i))

       END IF ! if mask
    END DO ! i loop

    !---------------------------------------------------------------------------
    !   If DIC and TA are known then either a root finding or iterative
    !   method must be used to calculate htotal. In this case we use
    !   the Newton-Raphson "safe" method taken from "Numerical Recipes"
    !   (function "rtsafe.f" with error trapping removed).
    !
    !   As currently set, this procedure iterates about 12 times. The
    !   x1 and x2 values set below will accomodate ANY oceanographic
    !   values. If an initial guess of the pH is known, then the
    !   number of iterations can be reduced to about 5 by narrowing
    !   the gap between x1 and x2. It is recommended that the first
    !   few time steps be run with x1 and x2 set as below. After that,
    !   set x1 and x2 to the previous value of the pH +/- ~0.5.
    !---------------------------------------------------------------------------

    CALL drtsafe_row(mask, x1, x2, xacc, htotal)

    !---------------------------------------------------------------------------
    !   Calculate [CO2*] as defined in DOE Methods Handbook 1994 Ver.2,
    !   ORNL/CDIAC-74, Dickson and Goyet, eds. (Ch 2 p 10, Eq A.49)
    !---------------------------------------------------------------------------

    DO i = 1,imt
       IF (mask(i)) THEN

          htotal2 = htotal(i) ** 2
          co2star(i) = dic(i) * htotal2 / &
               (htotal2 + k1(i)*htotal(i) + k1(i)*k2(i))
          co2starair = xco2(i) * ff(i) * atmpres(i)
          dco2star(i) = co2starair - co2star(i)
          ph(i) = -LOG10(htotal(i))

          !---------------------------------------------------------------------
          !   Add two output arguments for storing pCO2surf
          !   Should we be using K0 or ff for the solubility here?
          !---------------------------------------------------------------------

          pCO2surf(i) = co2star(i) / ff(i)
          dpCO2(i)    = pCO2surf(i) - xco2(i) * atmpres(i)

          !---------------------------------------------------------------------
          !   Convert units of output arguments
          !   Note: pCO2surf and dpCO2 are calculated in atm above.
          !---------------------------------------------------------------------

          co2star(i)  = co2star(i) * mass_to_vol
          dco2star(i) = dco2star(i) * mass_to_vol

          pCO2surf(i) = pCO2surf(i) * 1e6_dbl_kind
          dpCO2(i)    = dpCO2(i) * 1e6_dbl_kind

       ELSE ! if mask

          ph(i)       = c0
          co2star(i)  = c0
          dco2star(i) = c0
          pCO2surf(i) = c0
          dpCO2(i)    = c0

       END IF ! if mask
    END DO ! i loop

  END SUBROUTINE co2calc_row

  !*****************************************************************************

  SUBROUTINE talk_row(mask, x, fn, df)

    !---------------------------------------------------------------------------
    !   This routine computes TA as a function of DIC, htotal and constants.
    !   It also calculates the derivative of this function with respect to
    !   htotal. It is used in the iterative solution for htotal. In the call
    !   "x" is the input value for htotal, "fn" is the calculated value for
    !   TA and "df" is the value for dTA/dhtotal.
    !---------------------------------------------------------------------------

    !USE constants, ONLY : c1, c2, c3

    !---------------------------------------------------------------------------
    !   input arguments
    !---------------------------------------------------------------------------

    LOGICAL(KIND=log_kind), DIMENSION(imt), INTENT(IN) :: mask
    REAL(KIND=dbl_kind), DIMENSION(imt), INTENT(IN) :: x

    !---------------------------------------------------------------------------
    !   output arguments
    !---------------------------------------------------------------------------

    REAL(KIND=dbl_kind), DIMENSION(imt), INTENT(OUT) :: fn, df

    !---------------------------------------------------------------------------
    !   local variable declarations
    !---------------------------------------------------------------------------

    INTEGER(KIND=int_kind) :: i

    REAL(KIND=dbl_kind) :: &
         x1, x2, x3, k12, k12p, k123p, a, a2, da, b, b2, db, c

    !---------------------------------------------------------------------------

    DO i = 1,imt
       IF (mask(i)) THEN
          x1 = x(i)
          x2 = x1 * x1
          x3 = x2 * x1
          k12 = k1(i) * k2(i)
          k12p = k1p(i) * k2p(i)
          k123p = k12p * k3p(i)
          a = x3 + k1p(i) * x2 + k12p * x1 + k123p
          a2 = a * a
          da = c3 * x2 + c2 * k1p(i) * x1 + k12p
          b = x2 + k1(i) * x1 + k12
          b2 = b * b
          db = c2 * x1 + k1(i)
          c = c1 + st(i)/ks(i)

          !---------------------------------------------------------------------
          !   fn = hco3+co3+borate+oh+hpo4+2*po4+silicate-hfree-hso4-hf-h3po4-ta
          !---------------------------------------------------------------------

          fn(i) = k1(i) * x1 * dic(i)/b + &
               c2 * dic(i) * k12/b + &
               bt(i)/(c1 + x1/kb(i)) + &
               kw(i)/x1 + &
               pt(i) * k12p * x1/a + &
               c2 * pt(i) * k123p/a + &
               sit(i)/(c1 + x1/ksi(i)) - &
               x1/c - &
               st(i)/(c1 + ks(i)/x1/c) - &
               ft(i)/(c1 + kf(i)/x1) - &
               pt(i) * x3/a - &
               ta(i)

          !---------------------------------------------------------------------
          !   df = d(fn)/dx
          !---------------------------------------------------------------------

          df(i) = ((k1(i)*dic(i)*b) - k1(i)*x1*dic(i)*db)/b2 - &
               c2 * dic(i) * k12 * db/b2 - &
               bt(i)/kb(i)/(c1+x1/kb(i)) ** 2 - &
               kw(i)/x2 + &
               (pt(i) * k12p * (a - x1 * da))/a2 - &
               c2 * pt(i) * k123p * da/a2 - &
               sit(i)/ksi(i)/(c1+x1/ksi(i)) ** 2 - &
               c1/c + &
               st(i) * (c1 + ks(i)/x1/c)**(-2) * (ks(i)/c/x2) + &
               ft(i) * (c1 + kf(i)/x1)**(-2) * kf(i)/x2 - &
               pt(i) * x2 * (c3 * a - x1 * da)/a2

       END IF ! if mask
    END DO ! i loop

  END SUBROUTINE talk_row

  !*****************************************************************************

  SUBROUTINE drtsafe_row(mask_in, x1, x2, xacc, soln)

    !---------------------------------------------------------------------------
    !   Vectorized version of drtsafe, which was a modified version of
    !      Numerical Recipes algorithm.
    !   Keith Lindsay, Oct 1999
    !
    !   Algorithm comment :
    !      Iteration from Newton's method is used unless it leaves
    !      bracketing interval or the dx is > 0.5 the previous dx.
    !      In that case, bisection method is used.
    !---------------------------------------------------------------------------

    !USE constants, ONLY : c0, c2
    !USE shr_sys_mod, ONLY : shr_sys_abort

    !---------------------------------------------------------------------------
    !   input arguments
    !---------------------------------------------------------------------------

    LOGICAL(KIND=log_kind), DIMENSION(imt), INTENT(IN) :: mask_in
    REAL(KIND=dbl_kind), DIMENSION(imt), INTENT(IN) :: x1, x2
    REAL(KIND=dbl_kind), INTENT(IN) :: xacc

    !---------------------------------------------------------------------------
    !   output arguments
    !---------------------------------------------------------------------------

    REAL(KIND=dbl_kind), DIMENSION(imt), INTENT(OUT) :: soln

    !---------------------------------------------------------------------------
    !   local variable declarations
    !---------------------------------------------------------------------------

    LOGICAL(KIND=log_kind) :: leave_bracket, dx_decrease
    LOGICAL(KIND=log_kind), DIMENSION(imt) :: mask
    INTEGER(KIND=int_kind) ::  i, it
    REAL(KIND=dbl_kind) :: temp
    REAL(KIND=dbl_kind), DIMENSION(imt) :: xlo, xhi, flo, fhi, f, df, dxold, dx

    !---------------------------------------------------------------------------
    !   bracket root at each location and set up first iteration
    !---------------------------------------------------------------------------

    mask = mask_in

    CALL talk_row(mask, x1, flo, df)
    CALL talk_row(mask, x2, fhi, df)

    DO i = 1,imt
       IF (mask(i)) THEN
          IF (flo(i) .LT. c0) THEN
             xlo(i) = x1(i)
             xhi(i) = x2(i)
          ELSE
             xlo(i) = x2(i)
             xhi(i) = x1(i)
             temp = flo(i)
             flo(i) = fhi(i)
             fhi(i) = temp
          END IF
          soln(i) = 0.5 * (xlo(i) + xhi(i))
          dxold(i) = ABS(xlo(i) - xhi(i))
          dx(i) = dxold(i)
       END IF
    END DO

    CALL talk_row(mask, soln, f, df)

    !---------------------------------------------------------------------------
    !   perform iterations, zeroing mask when a location has converged
    !---------------------------------------------------------------------------

    DO it = 1,maxit
       DO i = 1,imt
          IF (mask(i)) THEN
             leave_bracket = ((soln(i)-xhi(i))*df(i)-f(i)) * &
                  ((soln(i)-xlo(i))*df(i)-f(i)) .GE. 0
             dx_decrease = ABS(c2 * f(i)) .LE. ABS(dxold(i) * df(i))
             IF (leave_bracket .OR. .NOT. dx_decrease) THEN
                dxold(i) = dx(i)
                dx(i) = 0.5 * (xhi(i) - xlo(i))
                soln(i) = xlo(i) + dx(i)
                IF (xlo(i) .EQ. soln(i)) mask(i) = .FALSE.
             ELSE
                dxold(i) = dx(i)
                dx(i) = -f(i) / df(i)
                temp = soln(i)
                soln(i) = soln(i) + dx(i)
                IF (temp .EQ. soln(i)) mask(i) = .FALSE.
             END IF
             IF (ABS(dx(i)) .LT. xacc) mask(i) = .FALSE.
          END IF
       END DO

       IF (.NOT. ANY(mask)) RETURN

       CALL talk_row(mask, soln, f, df)

       DO i = 1,imt
          IF (mask(i)) THEN
             IF (f(i) .LT. c0) THEN
                xlo(i) = soln(i)
                flo(i) = f(i)
             ELSE
                xhi(i) = soln(i)
                fhi(i) = f(i)
             END IF
          END IF
       END DO

    END DO ! iteration loop

    print *, 'lack of convergence in drtsafe_row, exiting'
    CALL exit(0)

  END SUBROUTINE drtsafe_row

  !*****************************************************************************

      function SCHMIDT_O2(SST)

!-----------------------------------------------------------------------
!     Compute Schmidt number of O2 in seawater as function of SST.
!
!     ref : Keeling et al, Global Biogeochem. Cycles, Vol. 12,
!     No. 1, pp. 141-163, March 1998
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     result & argument declarations
!-----------------------------------------------------------------------

      real(kind=dbl_kind), dimension(imt,jmt) :: SCHMIDT_O2

      real(kind=dbl_kind), dimension(imt,jmt), intent(in) :: SST

!-----------------------------------------------------------------------
!     coefficients in expansion
!-----------------------------------------------------------------------

      real(kind=dbl_kind), parameter :: &
        a = 1638.0_dbl_kind &
      , b = 81.83_dbl_kind &
      , c = 1.483_dbl_kind &
      , d = 0.008004_dbl_kind

!-----------------------------------------------------------------------

      SCHMIDT_O2 = a + SST * (-b + SST * (c + SST * (-d)))

!-----------------------------------------------------------------------

      end function SCHMIDT_O2

  !***********************************************************************

      function O2SAT(SST, SSS)

!-----------------------------------------------------------------------
!     Computes oxygen saturation concentration at 1 atm total pressure
!     in mol/m^3 given the temperature (t, in deg C) and the salinity (s,
!     in permil).
!
!     FROM GARCIA AND GORDON (1992), LIMNOLOGY and OCEANOGRAPHY.
!     THE FORMULA USED IS FROM PAGE 1310, EQUATION (8).
!
!     *** NOTE: THE "A_3*TS^2" TERM (IN THE PAPER) IS INCORRECT. ***
!     *** IT SHOULDN'T BE THERE.                                ***
!
!     O2SAT IS DEFINED BETWEEN T(freezing) <= T <= 40(deg C) AND
!     0 permil <= S <= 42 permil
!     CHECK VALUE:  T = 10.0 deg C, S = 35.0 permil,
!     O2SAT = 282.015 mmol/m^3
!-----------------------------------------------------------------------

      !use constants, only : T0_Kelvin

!-----------------------------------------------------------------------
!     result & argument declarations
!-----------------------------------------------------------------------

      real(kind=dbl_kind), dimension(imt,jmt) :: O2SAT

      real(kind=dbl_kind), dimension(imt,jmt), intent(in) :: &
        SST    & ! sea surface temperature (C)
      , SSS    ! sea surface salinity (psu)

!-----------------------------------------------------------------------
!     local variables
!-----------------------------------------------------------------------

      real(kind=dbl_kind), dimension(imt,jmt) :: TS

!-----------------------------------------------------------------------
!     coefficients in expansion
!-----------------------------------------------------------------------

      real(kind=dbl_kind), parameter :: &
        a_0 = 2.00907_dbl_kind &
      , a_1 = 3.22014_dbl_kind &
      , a_2 = 4.05010_dbl_kind &
      , a_3 = 4.94457_dbl_kind &
      , a_4 = -2.56847E-1_dbl_kind &
      , a_5 = 3.88767_dbl_kind &
      , b_0 = -6.24523E-3_dbl_kind &
      , b_1 = -7.37614E-3_dbl_kind &
      , b_2 = -1.03410E-2_dbl_kind &
      , b_3 = -8.17083E-3_dbl_kind &
      , c_0 = -4.88682E-7_dbl_kind

!-----------------------------------------------------------------------

      TS = LOG( ((T0_Kelvin+25.0_dbl_kind) - SST) / (T0_Kelvin + SST) )

      O2SAT = EXP(a_0+TS*(a_1+TS*(a_2+TS*(a_3+TS*(a_4+TS*a_5)))) + &
        SSS*( (b_0+TS*(b_1+TS*(b_2+TS*b_3))) + SSS*c_0 ))

!---------------------------------------------------------------------------
!   Convert from ml/l to mmol/m^3
!---------------------------------------------------------------------------

      O2SAT = O2SAT / 0.0223916_dbl_kind

!-----------------------------------------------------------------------

      end function O2SAT

  !***********************************************************************

      function SCHMIDT_CO2(SST)

!-----------------------------------------------------------------------
!     Compute Schmidt number of CO2 in seawater as function of SST.
!
!     ref : Wanninkhof, J. Geophys. Res, Vol. 97, No. C5,
!     pp. 7373-7382, May 15, 1992
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     result & argument declarations
!-----------------------------------------------------------------------

      real(kind=dbl_kind), dimension(imt,jmt) :: SCHMIDT_CO2

      real(kind=dbl_kind), dimension(imt,jmt), intent(in) :: SST

!-----------------------------------------------------------------------
!     coefficients in expansion
!-----------------------------------------------------------------------

      real(kind=dbl_kind), parameter :: &
        a = 2073.1_dbl_kind &
      , b = 125.62_dbl_kind &
      , c = 3.6276_dbl_kind &
      , d = 0.043219_dbl_kind

!-----------------------------------------------------------------------

      SCHMIDT_CO2 = a + SST * (-b + SST * (c + SST * (-d)))

!-----------------------------------------------------------------------

      end function SCHMIDT_CO2


end module kei_eco
