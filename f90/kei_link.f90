!!!     KPP-Ecosystem-Ice (KEI, pronounced 'key') Model
!!!     ==================================================================
!!!
!!!     This model derives from Large et al. [1994], Doney et al. [1996](KPP mixing),
!!!     Ukita and Martinson [2001](Mixed layer - ice interactions), Saenz and Arrigo
!!!     [2012 and 2014] (SIESTA sea ice model), Hunke and Lipscomb [2008] (CICE v4),
!!!     Moore et al. [2002,2004] (CESM ecosystem model, with some parameterizations
!!!     derived from Wang and Moore [2011]).
!!!
!!!     DEPENDENCIES:  1) A FORTRAN 90 compiler (gfortran and ifort have been tested)
!!!                    2) CPP pre-processor
!!!                    3) NetCDF 3.6.1 or higher, with FORTRAN 90 interface built
!!!                    4) LAPACK (used by SIESTA for solving heat transfer)
!!!                    5) MATLAB (not required to run compile or run KEI, but
!!!                        input data generation, run scripts, and analysis tools
!!!                        that greatly facilite use of KEI are MATLAB-based)!!!
!!!
!!!     Much of the original input/output code and options from the 1990s have been
!!!     removed in favor of standardizing on netcdf input/output files.  Some of the
!!!     code is still probably adaptable - contact Ben Saenz (blsaenz@gmail.com)
!!!     for older copies and more information.
!!!
!!!     As of Dec 2013, the atmospheric transfer model and coupling functionality
!!!     are unknown/likely broken. Currently the model is only being forced
!!!     at the atm/ocean and atm/ice interface by atmospheric parameters
!!!     derived from climatologies.  There is also much other legacy code
!!!     that should probably be deleted...
!!!
!!!     AS MODIFIED APR 18 1991 , WGL
!!!                 Jul 19      , jan
!!!                 Aug 12      , wgl
!!!                 Sep 18      , wgl
!!!                 Mar 17      , jan
!!!                 Apr 17      , jan : implicit integration
!!!                 Apr 27      , jan : include biology
!!!                 Jul 20      , jan : new restart from memory dump
!!!                 Nov 94      , wgl : for new KPP numerics
!!!                 Jun 04      , wgl : dgm ice model
!!!                 Dec 13      , BLS : Converted to Fortran 90,
!!!                                     SIESTA (Saenz and Arrigo 2012,2014) ice model,
!!!                                     netcdf data input/output, Los Alamos
!!!                                     CICE atmospheric/ocean interaction,
!!!                                     Moore et al. (CESM) Ecosystem model

MODULE link
  use kei_kinds, only: i4, r4, r8, log_kind
  !!! Use Statements and Variables/Globals/Parameters
  !!! --------------------------------------------------------------------
  USE kei_parameters
  USE kei_common
  USE kei_icecommon
  USE kei_ice
  USE kei_ecocommon
  USE kei_ocn, ONLY: init_ocn, ocnstep
  USE kei_eco
  USE kei_fluxes, ONLY: init_flx, atmflx, calflx, o2iflx, topflx
  USE kei_hacks
  USE kei_sw


  IMPLICIT NONE

  PUBLIC

  !!! Internal Variables
  !!! --------------------------------------------------------------------
  real(r4), save :: &
    U(NZP1,NVEL),  &    !!! momentum
    X(NZP1,NSCLR),  &   !!! scalars
    aflx(NSFLXSP2), &   !!!
    flx(11),        &   !!! diagnostic flux data structure
    Xprev(NZP1,2)       !!! T&S before (potential) data assimilation

  REAL(r8), save :: &
    timed, &            !!! model time
    dtday, &            !!! time step in days
    deltat, &           !!! time step
    HC(0:1), &          !!! heat content past/present
    SC(0:1), &          !!! salt content past/present
    HCocn, &            !!! heat content ocean
    SCocn, &            !!! salt content ocean
    HCice, &            !!! heat content ice
    SCice, &            !!! salt content ice
    FCice, &            !!!
    par_phyto(NZ), &    ! par (Watts) at depth from ecosys that may be needed by seaweed modules
    absorp_in(NZ), &    ! external sources shortwave irradiance absorption (1/m)
    absorp_out(NZ)      ! total ecosystem shortwave irradiance absorption (1/m)

  ! maybe can get rid of these?
  real(r4),save :: &
    atmtemp(NPLEV),  & !!! atmospheric model level temperature (K)
    atmshum(NPLEV),  & !!! atmospheric model level specific humidity(g/g)
    atmsave(NPSAVE), & !!! misc output from atm model to be stored
    atmsaveav(NPSAVE)  !!! averaged over storage time interval

  real(r4), save :: kforce(forcing_var_cnt)   !!! forcing data structure - revised for f2py
  real(r4), save :: absorp_ocn_step(NZ)       !!! cast for ``ocnstep`` (expects r4; ``absorp_out`` is r8)

  !!! counters and helper variables
  integer(i4), save :: i,ii,jptr,nt,ns_calls,nii,iaccum,nisteps,no, &
    ntlast,ntimelast,estep,top_f,top_unstable
  real(r4), save :: DHCdt,DSCdt,zml,DHSdt, &
    HCold,DHC,rho_top,rho_bot, &
    alpha,beta,exppr,v10


CONTAINS

!!! ************************************************************************
  SUBROUTINE set_tracers(U_in,X_in)
    real(r4), INTENT(IN) :: &
      U_in(NZ,NVEL),  &
      X_in(NZ,NSCLR)
    integer(i4) :: i

      U(1:NZ,:) = U_in
      X(1:NZ,:) = X_in

      ! copy last value the boundary cell, which is now invisible to python
      do i = 1, NVEL
        U(NZP1,i) = U(NZ,i)
      enddo
      do i = 1, NSCLR
        X(NZP1,i) = X(NZ,i)
      enddo

  END SUBROUTINE set_tracers

!!! ************************************************************************
  SUBROUTINE get_tracers(U_out,X_out)
    real(r4), INTENT(OUT) :: &
      U_out(NZ,NVEL),  &
      X_out(NZ,NSCLR)

      U_out = U(1:NZ,:)
      X_out = X(1:NZ,:)
      ! Internal storage keeps salinity anomaly (S - Sref); Python/output expects total salinity.
      X_out(:,2) = X_out(:,2) + Sref

  END SUBROUTINE get_tracers

!!! ************************************************************************
  SUBROUTINE get_sw(U_out,X_out)
    real(r4), INTENT(OUT) :: &
      U_out(NZ,NVEL),  &
      X_out(NZ,NSCLR)

      U_out = U(1:NZ,:)
      X_out = X(1:NZ,:)

  END SUBROUTINE get_sw

!!! ************************************************************************
  SUBROUTINE get_fluxes(flux_out)
    real(r4), INTENT(OUT) :: &
      flux_out(NSFLXS,5)

      flux_out = sflux(:,:,0)

  END SUBROUTINE get_fluxes

!!! ************************************************************************
  SUBROUTINE set_grid(dm_in,hm_in,zm_in)
    real(r4), INTENT(IN) :: &
      dm_in(NZ), hm_in(NZ), zm_in(NZ)
    real(r4), PARAMETER :: dzp1 = 0.1

      dm(0:NZM1) = dm_in
      hm(1:NZ) = hm_in
      zm(1:NZ) = zm_in

      hm(NZP1) = dzp1
      dm(NZ) = dm(NZM1) + dzp1 ! this is wrong? why does dm not have one more bound?
      zm(NZP1) = -1.0 * dm(NZM1) - dzp1/2.0

      zmp = abs(zm)

  END SUBROUTINE set_grid

!!! ************************************************************************
  SUBROUTINE set_forcing(kforce_in)
    real(r4), INTENT(IN) :: kforce_in(forcing_var_cnt)

      kforce = kforce_in

  END SUBROUTINE set_forcing


!!! ************************************************************************
  SUBROUTINE set_param_real(param,value)
    CHARACTER (len=*), intent(in) :: param
    real(r4), intent(in) :: value

    if (param == 'dlon') then
      dlon = value
    elseif (param == 'dlat') then
      dlat = value
    elseif (param == 'sal_correction_rate') then
      sal_correction_rate = value
    elseif (param == 'sw_scale_factor') then
      sw_scale_factor = value
    elseif (param == 'absorp_baseline') then
      absorp_baseline = value
    else
      print *,'set_param_real: param not understood or available!'
    endif

  END SUBROUTINE set_param_real

!!! ************************************************************************
  SUBROUTINE set_param_int(param,value)
    CHARACTER (len=*), intent(in) :: param
    integer(i4), intent(in) :: value
    if (param == 'leco') then
      leco = value
    elseif (param == 'lsw') then
      lsw = value
    elseif (param == 'lice') then
      lice = value
    elseif (param == 'nend') then
      nend = value
    elseif (param == 'jerlov') then
      jerlov = value
    else
      print *,'set_param_real: param not understood or available!'
    endif

  END SUBROUTINE set_param_int


!!! ************************************************************************
  SUBROUTINE get_data_real(param,value)
    CHARACTER (len=*), intent(in) :: param
    real(r4), intent(out) :: value

    if (param == 'time') then
      value = time
    elseif (param == 'dtsec') then
      value = REAL(dtsec, r4)
    !elseif (param == 'day') then
    !  value = day
    !elseif (param == 'hour') then
    !  value = hour
    elseif (param == 'hmx') then
      value = hmix
    elseif (param == 'zml') then
      value = zml

    elseif (param == 'hi') then
      value = hice(0)
    elseif (param == 'hs') then
      value = hsn(0)
    elseif (param == 'fice') then
      value = fice

    elseif (param == 'atm_flux_to_ocn_surface') then
      value = atm_flux_to_ocn_surface
    elseif (param == 'atm_flux_to_ice_surface') then
      value = atm_flux_to_ice_surface
    elseif (param == 'ice_ocean_bottom_flux') then
      value = ice_ocean_bottom_flux
    elseif (param == 'ice_ocean_bottom_flux_potential') then
      value = ice_ocean_bottom_flux_potential
    elseif (param == 'total_ice_melt') then
      value = total_ice_melt
    elseif (param == 'total_ice_freeze') then
      value = total_ice_freeze
    elseif (param == 'frazil_ice_volume') then
      value = frazil_ice_volume
    elseif (param == 'congelation_ice_volume') then
      value = congelation_ice_volume
    elseif (param == 'snow_ice_volume') then
      value = snow_ice_volume
    elseif (param == 'snow_precip_mass') then
      value = snow_precip_mass


    elseif (param == 'sw_sst') then
      value = sw_sst(1,1)
    elseif (param == 'sw_chl') then
      value = sw_chl(1,1)
    elseif (param == 'sw_par') then
      value = sw_par(1,1)
    elseif (param == 'sw_no3') then
      value = sw_no3_(1,1)
    elseif (param == 'sw_nh4') then
      value = sw_nh4_(1,1)
    elseif (param == 'sw_po4') then
      value = sw_po4_(1,1)
    elseif (param == 'sw_fe') then
      value = sw_fe_(1,1)
    elseif (param == 'sw_i') then
      value = DBLE(sw_i)
    elseif (param == 'sw_absorp') then
      value = sw_absorp

    else
        value = -99999.0
    endif

  END SUBROUTINE get_data_real

!!! ************************************************************************
  SUBROUTINE get_data_int(param,value)
    CHARACTER (len=*), intent(in) :: param
    integer(i4), intent(out) :: value

    if (param == 'ni') then
      value = ni_cur
    elseif (param == 'ns') then
      value = ns_cur

    else
        value = -99999
    endif

  END SUBROUTINE get_data_int


!!! ************************************************************************
  SUBROUTINE get_nz_data(param,nzp1_data)
    CHARACTER (len=*), intent(in) :: param
    real(r4), INTENT(OUT) :: &
      nzp1_data(NZ)

    ! nzp1 length
    if (param == 'wU') then
      nzp1_data = wU(1:NZ,1)
    elseif (param == 'wV') then
      nzp1_data = wU(1:NZ,2)
    elseif (param == 'wW') then
      nzp1_data = wU(1:NZ,3)
    elseif (param == 'wT') then
      nzp1_data = wX(1:NZ,1)
    elseif (param == 'wS') then
      nzp1_data = wX(1:NZ,2)
    elseif (param == 'wB') then
      nzp1_data = wX(1:NZ,3)
    elseif (param == 'Tprev') then
      nzp1_data = Xprev(1:NZ,1)
    elseif (param == 'Sprev') then
      nzp1_data = Xprev(1:NZ,2) + sref

    ! nz length
    elseif (param == 'par') then
      nzp1_data = par_phyto(1:NZ)
    elseif (param == 'tot_prod') then
      nzp1_data = tot_prod
    elseif (param == 'sp_Fe_lim') then
      nzp1_data = sp_Fe_lim
    elseif (param == 'sp_N_lim') then
      nzp1_data = sp_N_lim
    elseif (param == 'sp_P_lim') then
      nzp1_data = sp_P_lim
    elseif (param == 'sp_light_lim') then
      nzp1_data = sp_light_lim
    elseif (param == 'diat_Fe_lim') then
      nzp1_data = diat_Fe_lim
    elseif (param == 'diat_N_lim') then
      nzp1_data = diat_N_lim
    elseif (param == 'diat_P_lim') then
      nzp1_data = diat_P_lim
    elseif (param == 'diat_Si_lim') then
      nzp1_data = diat_Si_lim
    elseif (param == 'diat_light_lim') then
      nzp1_data = diat_light_lim
    elseif (param == 'graze_sp') then
      nzp1_data = graze_sp
    elseif (param == 'graze_diat') then
      nzp1_data = graze_diat
    elseif (param == 'graze_tot') then
      nzp1_data = graze_tot
    elseif (param == 'km') then
      nzp1_data = difm(1:NZ)
    elseif (param == 'ks') then
      nzp1_data = difs(1:NZ)
    elseif (param == 'kt') then
      nzp1_data = dift(1:NZ)
    elseif (param == 'ghat') then
      nzp1_data = ghat(1:NZ)

    else
      nzp1_data = -99999.0
    endif

  END SUBROUTINE get_nz_data


!!! ************************************************************************
  SUBROUTINE get_ice_data(param,ice_data)
    CHARACTER (len=*), intent(in) :: param
    real(r4), INTENT(OUT) :: &
      ice_data(nni)

    if (param == 'dzi') then
      ice_data = dzi
    elseif (param == 'Ti') then
      ice_data = Ti
    elseif (param == 'Si') then
      ice_data = Si

    else
      ice_data = -99999.0
    endif

  END SUBROUTINE get_ice_data


!!! ************************************************************************
  SUBROUTINE get_snow_data(param,snow_data)
    CHARACTER (len=*), intent(in) :: param
    real(r4), INTENT(OUT) :: &
      snow_data(nns)

    if (param == 'dzs') then
      snow_data = dzs
    elseif (param == 'Ts') then
      snow_data = Ts

    else
      snow_data = -99999.0
    endif

  END SUBROUTINE get_snow_data


!!! ************************************************************************
  SUBROUTINE get_sw_data(sw_tracer_block)
    REAL(r8), INTENT(OUT) :: &
      sw_tracer_block(n_sw_outputs)

    CALL sw_get_outputs(sw_tracer_block)

  END SUBROUTINE get_sw_data


!!! ************************************************************************
! Call this BEFORE setting parameters through the link interface ...
SUBROUTINE KEI_param_init

  !!! Physical & other constants
  CALL init_constants_params

  !!! init hacks, if needed
  CALL hacks_init()

END SUBROUTINE KEI_param_init



!!! ************************************************************************

!!! Must be called after loading forcing and initial conditions!!!
SUBROUTINE KEI_compute_init


  !!! Set time-method parameters (incl. ndtflx)
  CALL init_tm

  !!! Set ocean vertical grid, t-grid dimensions,
  !!!     rlat/rlon from dlat/dlon, and coriolis parameter
  !!!     load inital profiles
  CALL init_env(U,X)

  !!! initialize all kinds of stuff for the kei_eco module
  IF (leco) then
    CALL ecosys_init(hm,zmp)
    if (lsw) then
      CALL sw_init()
    endif
  endif

  !!! Atmosphere init "load VAF"
  !!!     timed=startt + dtsec/spd * float(ndtflx)/2.   !!!startt
  timed = time
  !CALL init_atm(timed,X,atmtemp,atmshum,kforce)

  !!! Initialize ice model:
  IF (lice)  THEN
    CALL init_ice(X(1,1),X(1,2)+Sref,HCice,FCice,SCice,startyear)
  ENDIF

  !!! initialize store/output with ns_calls=0, don't store + :return with ns=1
  !ns_calls = 0
  !IF (lstore) THEN
  !  CALL store (ns_calls,U,X,atmtemp,atmshum,atmsave,sflux,dout)
  !ENDIF

  ! subroutine prtprof is in "radinit.f" or when there is no
  !     rad/conv code provided as a dummy routine in "atm.f"
  !IF(LAFLX.eq.2 .or. LAFLX.eq.3 .or. LAFLX.eq.5) CALL prtprof

  ! Print input data to screen
  !CALL pinfo(6,U,X)

  !!! Initialize ocean model: set up temporary grids (if needed), and
  !!! compute initial hmix and Tref
  CALL init_ocn(U,X)

  last_wct_ii = 999999  !!! set index for previous assimilation depth to undefined

  !!! Initialize surface fluxes
  CALL init_flx(X,kforce,jptr,timed)
  flx=0.  !!! zero diagnostic flux structure

  !!! Initialize conservation tracking
  !CALL init_cons(U,X,HCice,FCice,SCice,HCocn,SCocn,HC,HCold,SC,DHCdt,DSCdt)

  !init_dout(X,zml,DHCdt,DSCdt)

  CALL mixedl(X(1,1),Tref,zml) ! X(:,1) ??
  Xprev = X(:,1:2)
  dtday = dtsec / 86400.D0
  nt    = nstart
  !!! perform initial store/write
  !IF(lstore) THEN
  !  CALL init_dout(X,zml,DHCdt,DSCdt)
  !  nt    = nstart
  !dtday = dtsec / spd
    !!! --- netcdf output routine ---- added 8/2011 Ben Saenz
  !  Xprev = X(:,1:2)
 !   CALL store_nc(nt,ns_calls,U,X,Xprev,sflux,zml,DHCdt,DSCdt,REAL(dtday),flx)
 !   CALL store(ns_calls,U,X,atmtemp,atmshum,atmsave,sflux,dout)
  !ENDIF

  !!! Close doc file (now ready for analysis while model is running)
  !CLOSE(nudoc)


  !!! Prepare main loop


END SUBROUTINE KEI_compute_init



SUBROUTINE KEI_compute_step(nt_in,doy)

  !!! TIME INTEGRATION main loop: do 1000 nt= nstart, nend-ndtflx, ndtflx
  !WRITE(6,*)
  !WRITE(6,*) '*******************************'
  !WRITE(6,*)     ' BEGIN  INTEGRATION LOOP'
  !WRITE(6,*)
  integer(i4), intent(in) :: nt_in
  REAL(r8), intent(in) :: doy ! day-of-year, 1-1 Jan
  REAL(r8) :: swh_out,mwp_out,cmag_out ! MACMODS expects doubles

    nt = nt_in

  ! DO WHILE (nt .le. nend-ndtflx)
    !!! time at midpoint of next flux interval
    timed = startt + dtday * ( nt + float(ndtflx) / 2. )

    !!! Get all fluxes on levels 1, 2, and 3 at timed
    CALL atm_step(kforce)
    CALL atmflx(jptr,timed)      !!! compute level 2, 3 and 1 fluxes
    CALL calflx(jptr)            !!! from past fluxes get jptr=0 fluxes
    flx(3) = sflux(6,2,jptr)
    flx(4) = sflux(5,2,jptr)
    flx(1) = sflux(3,2,jptr)
    flx(2) = sflux(4,2,jptr)
    flx(5) = 0.

    !!! ice sub-loop
    IF(lice) then

      IF(nt < zero_ice_before) THEN
        kforce(ic_f_ind) = 0.
        kforce(ain_f_ind) = 0.
        kforce(aout_f_ind) = 0.
        DivU = 0.
      ENDIF

      !!!load ocn-ice frazil ice fluxes + Tfrz at sfc
      CALL o2iflx(X,0)
      flx(5) = -sflux(7,4,0)*Fl + sflux(8,4,0) !!! new frazil growth + melt potential
      !!! Integrate ice model with this forcing, for nisteps, each ndtice long
      nisteps = 0
      DO nii = 1, ndtflx, ndtice
        nisteps = nisteps + 1
        !write(6,*) 'ice flux and step',nisteps,Tfrz
        !!! increment hice and hsn get new fice and TI0
        !!! load n=1,6 level 4 fluxes
        HCold = HCice
        v10 = sqrt(uz**2+vz**2)

        CALL icestep(nt,nisteps,timed, &                !!! timing
              X(1,1),X(1,2)+Sref,rho(0),TZ,Qz,hum,v10,msl, &
              DivU,kforce(ic_f_ind),kforce(ain_f_ind), & !!! fluxes/forcing
              kforce(aout_f_ind),sflux(:,:,0),NSFLXS, &    !!! fluxes/forcing
              HCice,FCice,SCice,flx,X)  !!!out:  heat,h2o & salt contents of ice
      ENDDO
    ENDIF  !!! end ice sub-loop


    !WRITE(6,*) 'end ice get ocn top'
    !WRITE(6,*) '*** focn: ',focn

    IF (eddy_hack) THEN
      CALL do_eddy_hack(nt,X)
    ENDIF

    !!! Get top of ocean fluxes then update focn
    atm_flux_to_ocn_surface = 0.
    CALL topflx(0,kforce)

    !CALL sprint(6,0)

    CALL assimilate_ts(X,nt)
    !!!CALL assimilate_ts_profile(X,kforce,nt)

    !!! record before-assimilation T & S
    Xprev(:,1) = X(:,1)
    Xprev(:,2) = X(:,2)

    !!! Integrate ocean model with this forcing over ndtflux interval
    DO 200 no = 1, ndtflx, ndtocn
         ntime = nt + no
         time = startt + dtday *  ntime
         absorp_ocn_step = real(absorp_out, kind=r4)
         CALL ocnstep(U,X,kforce,absorp_ocn_step) ! ``ocnstep`` dummy is r4; ecosystem keeps r8 in ``absorp_out``

    200  CONTINUE

      CALL HSC(U,X,HCice,FCice,SCice,HCocn,SCocn,HC,SC,DHCdt,DSCdt)
      !WRITE(6,*) 'Ocean output at  ',ntime,HC,SC
      !WRITE(6,*) HCocn,SCocn,HCice,SCice
      !CALL pprint(6,U,X)
      DHC =  HCice -  HCold
      !!!write(1,*) ntime,HCold,HCice,DHC,DHC/1800.
      !!!write(1,*) ntime,HC(iold),HC(inew),HCocn,HCice,DHCdt
      !!!write(2,*) ntime,SC(iold),SC(inew),SCocn,SCice,DSCdt

      !!! Set dout variables -- maybe
      !CALL set_dout(X,zml,DHCdt,DSCdt)
      CALL mixedl(X(1,1),Tref,zml) ! X(:,1) ??


      !!! Accumulate ocean model variables and fluxes, and compute averages
!       IF(lsaveaverages) THEN
!         CALL compav(ntime,inct,U,Uav,NZP1*NVEL)
!         CALL compav(ntime,inct,X,Xav,NZP1*NSCLR)
!         CALL compav(ntime,inct,dout(2),doutav(2),NDOUT-1)
!         doutav(1) = dout(1)  !!! do not average timestep
!         CALL compav(ntime,inct,sflux,sfluxav,5*NSFLXS)
!         CALL compav(ntime,inct,atmsave,atmsaveav,NPSAVE)
!       ENDIF

    !!! Accumulate biology fluxes
      absorp_in = absorp_baseline
      if (lsw .and. sw_i > 0) THEN
            absorp_in(sw_i) = absorp_in(sw_i) + sw_absorp
      endif
      IF(leco) THEN
        !!!do estep=1,60
          CALL ecosys_step(kforce,X,fice,albocn,Sref,dtsec,sflux(3,4,0),hmix,nt,startyear, &
                           par_phyto,absorp_in,absorp_out)
        !!!enddo
          if (lsw) THEN

              ! kei_sw_step(U,X,doy,dtday,nt,par_phyto,swh_in,mwp_in,cmag_in)
              swh_out = DBLE(kforce(swh_f_ind))
              mwp_out = DBLE(kforce(mwp_f_ind))
              cmag_out = DBLE(kforce(cmag_f_ind))
              call kei_sw_step(U,X,doy,dtday,nt,par_phyto,swh_out,mwp_out,cmag_out)

          endif


        ENDIF

    !!! Store and/or Print ?
!      IF( mod(ntime,inct).eq.0 ) THEN
!        IF( lstore ) THEN
!          !!! --- netcdf output routine ---- added 8/2011 Ben Saenz
!          CALL store_nc(nt,ns_calls,U,X,Xprev,sflux,zml,DHCdt,DSCdt,REAL(dtday),flx)
!          !!! -----------------------------------------------------
!          IF(lsaveaverages) THEN
!             CALL store (ns_calls,Uav,Xav,atmtemp,atmshum,atmsaveav, &
!                           sfluxav,doutav)
!          else
!             CALL store (ns_calls,U,X,atmtemp,atmshum,atmsave,sflux,dout)
!          ENDIF
!        else
!          !!! print out results at this timestep
!          write(6,*)
!          write(6,*) 'time step',ntime,dout(1),(dout(i),i=8,10)
!          CALL pprint(6,U,X)
!          CALL sprint(6,0)
! !!!             CALL kprint(6)
!          CALL fprint(6)
!          IF(lice) CALL iprint(6)
!          IF(LAFLX.eq.2 .or. LAFLX.eq.3 .or. LAFLX.eq.5) &
!            CALL atmprint(6,atmtemp,atmshum,atmsave)
!        ENDIF
!      ENDIF

      ntlast = nt

      !!! CONTINUE TIME INTEGRATION main loop
      focn = 1. - fice - flnd
      iold = inew
      inew = 1 - iold
      nt = nt + ndtflx

    ! ENDDO  !!! end of 1000 / main loop
!!!

!!! Store everything in memory for restart
!!!
!  CALL storerestart(U,X,aflx,atmtemp,atmshum,atmsave,atmsaveav, &
!                     jptr,ntlast,ntimelast,'.rst      ')
!!!      CALL writerestart(U,X,aflx,atmtemp,atmshum,atmsave,atmsaveav, &
!!!                       jptr,nt)
!!!
END SUBROUTINE KEI_compute_step


!!! **********************************************************************
!!! subroutine init_cons
!!! Purpose: Initialize conservation tracking variables
!!! ----------------------------------------------------------------------
  SUBROUTINE init_cons(U,X,HCice,FCice,SCice,HCocn,SCocn,HC,HCold,SC,DHCdt,DSCdt)

    USE kei_parameters
    USE kei_common
    USE kei_icecommon

    IMPLICIT NONE

    real(r4) :: U(NZP1,NVEL), X(NZP1,NSCLR),DHCdt,DSCdt,HCold
    REAL(r8) :: HC(0:1),SC(0:1),HCocn,SCocn,HCice,FCice,SCice

    HC(iold) =  0.0
    SC(iold) =  0.0
    CALL HSC(U,X,HCice,FCice,SCice,HCocn,SCocn,HC,SC,DHCdt,DSCdt)
    HC(iold) = HC(inew)
    SC(iold) = SC(inew)
    DHCdt    = 0.0
    DSCdt    = 0.0
    write(6,*)
    write(6,*) ' initial conditions',HCocn,HCice,HC(iold), &
                                    SCocn,SCice,SC(iold)
    HCold = 0.0


  END SUBROUTINE init_cons
!!! **********************************************************************



!!! **********************************************************************
!!! subroutine mixedl
!!! Purpose: calculate mixed-layer depth based on delta T
!!! ----------------------------------------------------------------------
  SUBROUTINE mixedl(t,tsfc,zml)
!!!     compute mixed layer depth zml: zml<=0.

    USE kei_parameters !!! include 'parameter.inc'
    USE kei_common !!!include 'common.inc'

    implicit none

    real(r4), intent(in) :: t(nzp1),tsfc
    real(r4), intent(out) :: zml
    real(r4) :: delt
    integer(i4) :: kz,k

    delt=0.35   !!! temperature difference for mixed layer depth
    kz=1

    IF( t(1).gt.(tsfc+delt) ) THEN
       zml = zm(1)* (-delt) / (tsfc-t(1))
       GOTO 151
    ENDIF

    IF( t(1).lt.(tsfc-delt) ) THEN
       zml = zm(1)* (delt) / (tsfc-t(1))
       k=1
       GOTO 151
    ENDIF

    DO 15 k=2,nz+1
       IF( t(k).gt.(tsfc+delt) ) THEN
          zml = zm(k-1) + (zm(k)-zm(k-1)) &
                    * (t(k-1)-tsfc-delt) / (t(k-1)-t(k))
          GOTO 151
       ENDIF
       IF( t(k).lt.(tsfc-delt) ) THEN
          zml = zm(k-1) + (zm(k)-zm(k-1)) &
                    * (t(k-1)-tsfc+delt) / (t(k-1)-t(k))
          GOTO 151
       ENDIF
15   CONTINUE
    k   = nzp1
    zml = zm(k)          !!! max value
151  CONTINUE
!!!      write(6,*) ' k=',k,zm(k),zm(k-1),zml
!!!      write(6,*) '       ',t(k),t(k-1),tsfc

    RETURN
  END SUBROUTINE mixedl

!!! ************************************************************************


!***************************************************************

    Subroutine HSC(U,X,HCice,FCice,SCice,     & !  inputs
    HCocn,SCocn,HC,SC,DHCdt,DSCdt ) ! outputs

!     Subroutine HSC(U,X,HCice,FCine,SCice,HCocn,SCocn,HC,SC,DHCdt)

!      Compute the Heat, Salt and Freshwater Content of a profile
!      Iputs from ice.f are heat, h20 and salt  kg/m2 of the ice.

    use kei_parameters
    use kei_common
    use kei_icecommon

    implicit none

		! inputs
    real(r4) :: U(NZP1,NVEL), X(NZP1,NSCLR)
    REAL(r8) :: HCice,FCice,SCice

		! outputs
    REAL(r8) :: HC(0:1),SC(0:1),HCocn,SCocn
    real(r4) :: DHCdt,DSCdt

		! local
    integer(i4) :: n
    real(r4) :: rhoCP,Szero,rhof,rhoocn,Sal,Saltice

!                                 Ocean Heat Content
    HCocn = 0.0
    SCocn = 0.0

    Do 10 n = 1,NZ
    !                           Allow for conservation
        if(LAFLX == 0) then
            rhoCP = rhosw * CPo
            Szero = Sref / 1000.
            rhof  = rhofrsh
            rhoocn= rhosw
        else
            rhoCP = rho(n)*CP(n)
            Szero = (Sref + X(n,2)) / 1000.
            rhof  = rhoh2o
            rhoocn= rho(n)
        endif

        HCocn = HCocn + rhoCP * (X(n,1)) * hm(n)          ! J/m2
        Sal   = (Sref + X(n,2))  /1000.
        SCocn = SCocn +  hm(n) * Sal    ! kgS/kgsw * m
    10 END DO

!        Use (Ice and Snow) Heat,h2o and Salt  Contents from icestep

    HC(inew) = HCocn + fice * HCice
    Saltice  = -FCice * Szero / rhof + SCice / rhoocn ! kgS/kgsw *m
    SCice    = fice  * Saltice
    SC(inew) = SCocn + SCice

    DHCdt    = (HC(inew) - HC(iold) ) / (dtsec * ndtflx)
    DSCdt    = (SC(inew) - SC(iold) ) / (dtsec * ndtflx)

    return
    end Subroutine HSC


END MODULE link





