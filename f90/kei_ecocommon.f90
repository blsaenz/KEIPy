module kei_ecocommon

  use kei_parameters
	implicit none

	Public
	Save

	  integer, parameter :: &
          log_kind = kind(.true.);

	  integer, parameter :: &
          int_kind = 4;

	  integer, parameter :: &
          real_kind = 4;

	  integer, parameter :: &
          dbl_kind = 8;

     integer (kind=int_kind), parameter :: &
          ecosys_tracer_cnt = 24

     ! number of vertical layers
     ! sourced using equivalence to kei_parameters, for code compatibility
     integer (kind=int_kind), parameter :: &
          km = NZ
     !equivalence (nz,km)

      integer (kind=int_kind), parameter :: &
        imt = 1,  &   ! x (longitudinal) dimension of ecosys variables
        jmt = 1       ! y (latitudinal) dimension of ecosys variables

!-----------------------------------------------------------------------
!     tracer indices
!-----------------------------------------------------------------------

    integer(kind=int_kind), parameter :: &
      po4_ind     =  1,  & ! dissolved inorganic phosphate
      no3_ind     =  2,  & ! dissolved inorganic nitrate
      sio3_ind    =  3,  & ! dissolved inorganic silicate
      nh4_ind     =  4,  & ! dissolved ammonia
      fe_ind      =  5,  & ! dissolved inorganic iron
      o2_ind      =  6,  & ! dissolved oxygen
      dic_ind     =  7,  & ! dissolved inorganic carbon
      alk_ind     =  8,  & ! alkalinity
      doc_ind     =  9,  & ! dissolved organic carbon
      spC_ind     = 10,  & ! small phytoplankton carbon
      spChl_ind   = 11,  & ! small phytoplankton chlorophyll
      spCaCO3_ind = 12,  & ! small phytoplankton caco3
      diatC_ind   = 13,  & ! diatom carbon
      diatChl_ind = 14,  & ! diatom chlorophyll
      zooC_ind    = 15,  & ! zooplankton carbon
      spFe_ind    = 16,  & ! small phytoplankton iron
      diatSi_ind  = 17,  & ! diatom silicon
      diatFe_ind  = 18,  & ! diatom iron
      diazC_ind   = 19,  & ! diazotroph carbon
      diazChl_ind = 20,  & ! diazotroph Chlorophyll
      diazFe_ind  = 21,  & ! diazotroph iron
      don_ind     = 22,  & ! dissolved organic nitrogen
      dofe_ind    = 23,  & ! dissolved organic iron
      dop_ind     = 24     ! dissolved organic phosphorus

    character (len = 8), dimension(ecosys_tracer_cnt), &
      parameter :: eco_tracer_name = (/ &
        'PO4     ', &
        'NO3     ', &
        'SiO3    ', &
        'NH4     ', &
        'Fe      ', &
        'O2      ', &
        'DIC     ', &
        'ALK     ', &
        'DOC     ', &
        'spC     ', &
        'spChl   ', &
        'spCaCO3 ', &
        'diatC   ', &
        'diatChl ', &
        'zooC    ', &
        'spFe    ', &
        'diatSi  ', &
        'diatFe  ', &
        'diazC   ', &
        'diazChl ', &
        'diazFe  ', &
        'DON     ', &
        'DOFe    ', &
        'DOP     ' /)

!-----------------------------------------------------------------------
!     terms that may be of interest to output
!-----------------------------------------------------------------------

       real (kind=dbl_kind), dimension(:,:), allocatable, target :: &
        XKW_tavg, AP_tavg, PV_CO2_tavg, PV_O2_tavg, &
        SCHMIDT_O2_tavg, O2SAT_tavg, &
        FG_O2_tavg, SCHMIDT_CO2_tavg, PH_tavg, CO2STAR_tavg, &
        DCO2STAR_tavg, pCO2SURF_tavg, DpCO2_tavg, FG_CO2_tavg, &
        IRON_FLUX_tavg, PROD_tavg, PO4_RESTORE_tavg, NO3_RESTORE_tavg, &
        SiO3_RESTORE_tavg, PAR_avg_tavg, PO4STAR_tavg, POC_FLUX_IN_tavg, &
        POC_PROD_tavg, POC_REMIN_tavg, CaCO3_FLUX_IN_tavg, &
        CaCO3_PROD_tavg, CaCO3_REMIN_tavg, SiO2_FLUX_IN_tavg, &
        SiO2_PROD_tavg, SiO2_REMIN_tavg, dust_FLUX_IN_tavg, &
        dust_REMIN_tavg, P_iron_FLUX_IN_tavg, P_iron_PROD_tavg, &
        P_iron_REMIN_tavg, graze_sp_tavg, graze_diat_tavg, &
        graze_tot_tavg, sp_loss_tavg, diat_loss_tavg, zoo_loss_tavg, &
        sp_agg_tavg, diat_agg_tavg, photoC_sp_tavg, photoC_diat_tavg, &
        tot_prod_tavg, DOC_prod_tavg, DOC_remin_tavg, Fe_scavenge_tavg, &
        sp_N_lim_tavg, sp_Fe_lim_tavg, sp_PO4_lim_tavg, &
        sp_light_lim_tavg, diat_N_lim_tavg, diat_Fe_lim_tavg, &
        diat_PO4_lim_tavg, diat_SiO3_lim_tavg, diat_light_lim_tavg, &
        CaCO3_form_tavg, diaz_Nfix_tavg, graze_diaz_tavg,diaz_loss_tavg, &
        photoC_diaz_tavg, diaz_P_lim_tavg, diaz_Fe_lim_tavg, &
        diaz_light_lim_tavg, Fe_scavenge_rate_tavg, DON_prod_tavg, &
        DON_remin_tavg, DOFe_prod_tavg, DOFe_remin_tavg, DOP_prod_tavg, &
        DOP_remin_tavg, bSi_form_tavg, photoFe_diaz_tavg, &
        photoFe_diat_tavg, photoFe_sp_tavg, FvPE_DIC_tavg,  &
        FvPE_ALK_tavg, NITRIF_tavg, DENITRIF_tavg

       ! output parameters
       real (kind=real_kind), dimension(km) :: &
          tot_prod, diat_Fe_lim, diat_light_lim, graze_diat, graze_tot, &
          diat_N_lim, diat_P_lim, diat_Si_lim, sp_N_lim, sp_P_lim, sp_Fe_lim, &
          graze_sp, sp_light_lim

!-----------------------------------------------------------------------
!     forcing variables required for ecosys module
!-----------------------------------------------------------------------

		real (kind=dbl_kind), parameter :: &
		 atm_co2_const = 400.0_dbl_kind, &  ! default constant CO2
		 ap_const = 1.0_dbl_kind             ! default constant air pressure (atm)

		real(kind=dbl_kind) :: &
		 dust_flux,          & ! surface dust flux
!       iron_flux,         & ! iron component of surface dust flux
		 winds_SQR,          & ! wind-speed ** 2
		 xkw,                & ! a * U10_SQR
		 atm_co2,            & ! atmospheric CO2 concentration
		 ap                    ! atmospheric pressure

		real(kind=dbl_kind), dimension (ecosys_tracer_cnt) :: &
      ice_to_ocean_eflux

		! restoring nutrient climatology switches
		logical, parameter :: &
			lrest_po4 = .false. , &        ! po4 climatological restoring switch
			lrest_no3 = .false. , &        ! no3 climatological restoring switch
			lrest_sio3 = .false. , &       ! sio3 climatological restoring switch
			lflux_gas_co2 = .true., &      ! atmospheric co2 flux switch
		  lflux_gas_o2 = .true., &       ! atmospheric o2 flux switch
		  lsource_sink = .true.          ! T = compute ecosys, F = inorganic only

		! restoring nutrient climatology profiles
		real(kind=dbl_kind), dimension(imt,jmt,km) :: &
			PO4_CLIM, NO3_CLIM, SiO3_CLIM

		! timescales and depths for nutrient restoring
		real (kind=dbl_kind), parameter :: &
			rest_time_inv_surf = 0.0_dbl_kind, &  ! 0 = instant, inverse ?
			rest_time_inv_deep = 0.0_dbl_kind, &  ! 0 = instant, inverse ?
			rest_z0            = 1000.0_dbl_kind, &  ! m
			rest_z1            = 2000.0_dbl_kind     ! m


end module kei_ecocommon
