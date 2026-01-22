
!Notes on fortran building (gfortran):
!f2py -c mag_kinds_mod.f90 mag_parameters_mod.f90 mag_calc.f90 -m mag_calc_fortran --debug --f90flags="-m64 -g -c -O3 -fopenmp -fdec-math -fopenmp -march=native -funroll-loops" -lgomp
!
MODULE macmods_calc

    !use iso_c_binding, only: c_int, c_float, c_double
    use macmods_kinds_mod
    use macmods_param_mod
    use macmods_util_mod
    implicit none

    PUBLIC

 CONTAINS
!
!
!   SUBROUTINE c_mag_calc(lat,lon,sst,par,swh,cmag,no3,nflux,mask, &
!                       params, &
!                       Q,B,d_B,d_Q,Growth2,d_Be,d_Bm,d_Ns,harv,GRate,B_N)
!     real(c_float),   dimension, intent(in) :: lat,lon,sst,par,swh,cmag,no3,nflux
!     integer(c_int),  dimension, intent(in) :: mask
!     real(c_float),   dimension(npar), intent(in) :: params
!     real(c_float),   dimension, intent(inout) :: Q,B,d_B,d_Q,Growth2,d_Be,d_Bm,d_Ns,harv,GRate,B_N
!     CALL mag_calc(lat,lon,sst,par,swh,cmag,no3,nflux,mask, &
!                   params, &
!                   Q,B,d_B,d_Q,Growth2,d_Be,d_Bm,d_Ns,harv,GRate,B_N)
!   END SUBROUTINE c_mag_calc

  SUBROUTINE mag_calc(lat,lon,sst,par,chl,swh,mwp,cmag,nflux, &
                      no3,nh4,po4,fe, &
                      do_harvest,seed_now, mask, GD_count, B_0000, n_harv, t_harv, &
                      params, &
                      Q,QP,QFe,B,d_B,d_Q,Growth2,d_Be,d_Bm,d_Bm_wave,harv,GRate,B_N,Gave,Dave, &
                      d_NO3,d_NH4,d_PO4,d_Fe, &
                      d_DIC,d_O2, &
                      d_DOC,d_DON,d_DOP,d_DOFe, &
                      d_POC,d_PON,d_POP,d_POFe, &
                      min_lim,gQout,gTout,gEout,gHout, &
                      turn_of_day,dt_mag)

    !use macmods_kinds_mod
    !use macmods_parameters_mod
    !use OMP_LIB, only : omp_set_num_threads

    implicit none

    ! input/output variables
    !----------------------------------------------------
    !integer(kind=i4), parameter :: imt = 2160  ! apparently can't pull this from module :\
    !integer(kind=i4), parameter :: jmt = 4320
    !integer(kind=i4), parameter :: npar1 = 70

    real(kind=r8), INTENT(IN), DIMENSION(imt,jmt) :: &
    lat, &
    lon, &
    sst, &
    par, &
    chl, &
    swh, &
    mwp, &
    cmag, &
    nflux

    real(kind=r8), INTENT(INOUT), DIMENSION(imt,jmt) :: &
    no3, &
    nh4, &
    po4, &
    fe 

    integer(kind=i4), INTENT(IN), DIMENSION(imt,jmt) :: &
    do_harvest

    integer(kind=i4), INTENT(INOUT), DIMENSION(imt,jmt) :: &
    seed_now, &
    mask, &
    GD_count
    !n_harv, &
    !t_harv, &
    !min_lim

    !integer(kind=i4), INTENT(IN) :: &
    !ny,ny,npar

    real(kind=r8), INTENT(IN) :: &
      params(npar), &
      dt_mag

    real(kind=r8), INTENT(INOUT), DIMENSION(imt,jmt) :: &
    Q, &
    QP, &
    QFe, &
    B, &
    B_0000, &
    d_B, &
    d_Q, &
    Growth2, &
    d_Be, &
    d_Bm, &
    d_Bm_wave, &
    d_NO3, &
    d_NH4, &
    d_PO4, &
    d_Fe, &
    d_DIC, &
    d_O2, &
    d_DOC, &
    d_DON, &
    d_DOP, &
    d_DOFe, &
    d_POC, &
    d_PON, &
    d_POP, &
    d_POFe, &
    harv, &
    GRate, &
    B_N, &
    Gave, &
    Dave, &
    gQout, &
    gTout, &
    gEout, &
    gHout, &
    n_harv, &
    t_harv, &
    min_lim

    logical, INTENT(IN) :: &
    turn_of_day

    ! shared internal variables
    !----------------------------------------------------
    real(kind=r8) ::  & ! real params
      KsNO3, & 
      PARs, &
      PARc, &
      Qmin,   &
      Qmax,   &
      Qmin_P,   &
      Qmax_P,   &
      Qmin_Fe,   &
      Qmax_Fe, &
      spp_NtoC_mass, &
      spp_NtoP_mass, &
      spp_NtoFe_mass 

    integer(kind=i4) :: &
      uptake_type, use_NH4, use_P, use_Fe, harvest_type, harvest_nmax, & ! integer params
      growth_lim_type, Q_lim_type ! integer params

    ! private internal variables
    !----------------------------------------------------
    logical :: bad_forcing,is_valid

    integer(kind=i4) :: i,j,thread

    real(kind=r8) :: &
      day_h,         &
      vQ_N,            &
      vQ_P,            &
      vQ_Fe,            &
      Growth,     &
      WP,     &
      M_Wave,     &
      M,     &
      par_watts,    &
      atten,    &
      chlmin,    &
      A,    &
      N0,    &
      harv1,  &
      mlim,   &
      Tlim,   &
      Llim,   &
      Qlim,   &
      kcap_slope, &
      B_calc, &
      G_calc, &
      Gmax_crowding, &
      dBdt, &
      dNdt, &
      dPdt, &
      dFedt, &
      Q_new, &
      Uptake_NO3,     &
      Uptake_NH4, &
      Uptake_PO4, &
      Uptake_Fe, &
      UptakeN,     &
      UptakeP, &
      UptakeFe, &
      gQ, &
      gQP, &
      gQFe, &
      N_exude, &
      C_exude, &
      e_scale, &
      exude, &
      Bm_loss, &
      Bm_wave, &
      Bm_POM, &
      Bm_DOM, &
      Bm_remin, &
      d0_no3, &
      d0_nh4, &
      d0_po4, &
      d0_fe


    ! whole model domain calculations for current timestep
    ! ------------------------------------------------------------
    KsNO3 = params(mp_spp_Ks_NO3)
    PARs = params(mp_spp_PARs)
    PARc = params(mp_spp_PARc)
    Qmin = params(mp_spp_Qmin)
    Qmax = params(mp_spp_Qmax)
    Qmin_P = params(mp_spp_Qmin_P)
    Qmax_P = params(mp_spp_Qmax_P)
    Qmin_Fe = params(mp_spp_Qmin_Fe)
    Qmax_Fe = params(mp_spp_Qmax_Fe)
    uptake_type = int(params(mp_N_uptake_type))
    growth_lim_type = int(params(mp_growth_lim_type))
    Q_lim_type = int(params(mp_Q_lim_type))
    use_NH4  = int(params(mp_use_NH4))
    use_P = int(params(mp_use_P))
    use_Fe = int(params(mp_use_Fe))
    harvest_type = int(params(mp_harvest_type))
    harvest_nmax = int(params(mp_harvest_nmax))
    spp_NtoC_mass = mp(mp_spp_NtoC) * mw_n/mw_c ! mw_n/mw_c
    spp_NtoP_mass = mp(mp_spp_NtoP) * mw_n/mw_p ! mw_n/mw_p
    spp_NtoFe_mass = mp(mp_spp_NtoFe) * mw_n/mw_Fe ! mw_n/mw_Fe

  !CALL omp_set_num_threads(10)

  !$OMP PARALLEL &
  !$OMP DEFAULT(SHARED) &
  !$OMP PRIVATE(bad_forcing,is_valid,i,j,thread, &
  !$OMPday_h,vQ_N,vQ_P,vQ_Fe,Growth,WP,M_Wave,M,par_watts,atten,chlmin,A,N0,harv1,mlim,Tlim,Llim,Qlim, &
  !$OMPkcap_slope,B_calc,G_calc,Gmax_crowding,dBdt,dNdt,dPdt,dFedt,Q_new,Uptake_NO3,Uptake_NH4, &
  !$OMPUptake_PO4,Uptake_Fe,UptakeN,UptakeP,UptakeFe,gQ,gQP,gQFe,N_exude,C_exude,e_scale,exude,Bm_loss, &
  !$OMPBm_POM,Bm_DOM,Bm_remin,d0_no3,d0_nh4,d0_po4,d0_fe)
  !$OMP DO SCHEDULE(DYNAMIC,chunk_size)

    ! ----------------------------------------------------------------
    !   do i=1708,1710
    !do j=1001,1001
    !do j=425,425
    do j=1,jmt
      !do i=1001,1001
      !do i=375,375
      do i=1,imt
        !print *,"Dude",i,j
      if (mask(i,j) > 0) then

        if (seed_now(i,j) > 0) then
            B(i,j) = params(mp_spp_seed)
            Q(i,j) = Qmin + no3(i,j)*(Qmax-Qmin)/35.0_r8
            QP(i,j) = MIN(Q(i,j) / spp_NtoP_mass, QMax_P)
            QP(i,j) = MAX(QP(i,j), QMin_P)
            QFe(i,j) = MIN(Q(i,j) / spp_NtoFe_mass, QMax_Fe)
            QFe(i,j) = MAX(QFe(i,j), QMin_Fe)
            t_harv(i,j) = 0
        endif

        ! Nutrient Uptake
        ! ----------------------------------------------------------(cmag,mwp,CD,Vmax,Ks,Nconc,vQ)
        vQ_N = uptake_velocity(Q(i,j), Qmin, Qmax, uptake_type)
        Uptake_NO3 = nutrient_uptake(cmag(i,j),mwp(i,j),params(mp_spp_CD), &
                                     params(mp_spp_Vmax),params(mp_spp_Ks_NO3), &
                                     no3(i,j),vQ_N)
        Uptake_NO3 = Uptake_NO3 * mw_n / 1.e3_r8
        UptakeN = Uptake_NO3 / params(mp_spp_dry_sa)
        !print *,lambda,vQ,vNuTw_NO3,UptakeN,Q(i,j)

        gQ = Q_limitation(Q(i,j),Qmax,Qmin,Q_lim_type)
        gQP = c1
        gQFe = c1
        if (use_NH4 == 1) then
            ! ammonia uptake, add to total N uptake
            Uptake_NH4 = nutrient_uptake(cmag(i,j),mwp(i,j),params(mp_spp_CD), &
                                         params(mp_spp_Vmax_NH4),params(mp_spp_Ks_NH4), &
                                         nh4(i,j),vQ_N)
            Uptake_NH4 = Uptake_NH4 * mw_n / 1.e3_r8
            UptakeN = UptakeN + Uptake_NH4 / params(mp_spp_dry_sa)
        endif

        if (use_P == 1) then
            ! Phosphorus uptake and quotient-based growth limitation
            vQ_P = uptake_velocity(QP(i,j), Qmin_P, Qmax_P, uptake_type)
            Uptake_PO4 = nutrient_uptake(cmag(i,j),mwp(i,j),params(mp_spp_CD), &
                                         params(mp_spp_Vmax_PO4),params(mp_spp_Ks_PO4), &
                                         po4(i,j),vQ_P)
            Uptake_PO4 = Uptake_PO4 * mw_p / 1.e3_r8
            UptakeP = Uptake_PO4 / params(mp_spp_dry_sa)

            gQP = Q_limitation(QP(i,j),Qmax_P,Qmin_P,Q_lim_type)
        endif

        if (use_Fe == 1) then
            ! Iron uptake and quotient-based growth limitation
            vQ_Fe = uptake_velocity(QFe(i,j), Qmin_Fe, Qmax_Fe, uptake_type)
            Uptake_Fe = nutrient_uptake(cmag(i,j),mwp(i,j),params(mp_spp_CD), &
                                         params(mp_spp_Vmax_Fe),params(mp_spp_Ks_Fe), &
                                         fe(i,j),vQ_Fe)
            Uptake_Fe = Uptake_Fe * mw_Fe / 1.e3_r8
            UptakeFe = Uptake_Fe / params(mp_spp_dry_sa)

            gQFe = Q_limitation(QFe(i,j),Qmax_Fe,Qmin_Fe,Q_lim_type)
        endif

        ! find DIC limitation ?


        ! Growth
        ! ----------------------------------------------------------
        ! Growth, nitrogen movement from Ns to Nf = umax*gQ*gT*gE*gH; [per day]
        ! Output:
        !   Growth, [h-1]
        !   gQ, quota-limited growth
        !       from Wheeler and North 1980 Fig. 2
        !   gT, temperature-limited growth
        !       piecewise approach taken from Broch and Slagstad 2012 (for sugar
        !       kelp) and optimized for Macrocystis pyrifera
        !   gE, light-limited growth
        !       from Dean and Jacobsen 1984
        !   gH, carrying capacity-limited growth


        ! temperature limitation
        Tlim = temp_lim(sst(i,j),params(mp_spp_Topt1),params(mp_spp_K1), &
                      params(mp_spp_Topt2),params(mp_spp_K2))

        ! light limitation
        par_watts = par(i,j) ! *2.515376387217542_r8 <-- no longer internally converting par - should be in W/m2
        ! clh < 0.0 means par-at-depth (in Watts) is already passed in
        if ( chl(i,j) > c0 ) then
            ! attentuation according to MARBL
            chlmin = max(0.02_r8,chl(i,j))
            chlmin = min(30.0_r8,chlmin)
            if (chlmin < 0.13224_r8) then
                atten = -0.000919_r8*(chlmin**0.3536_r8) ! 1/cm
            else
                atten = -0.001131_r8*(chlmin**0.4562_r8) ! 1/cm
            endif
            par_watts = par_watts * exp(atten*params(mp_farm_depth)*100.0_r8)
        endif

        if (par_watts < PARc) then
          Llim = c0
        elseif (par_watts > PARs) then
          Llim = c1
        else
          Llim = (par_watts-PARc)/(PARs-PARc)*exp(-(par_watts-PARc)/(PARs-PARc)+c1)
        endif

        ! consider daylength if timestep is > 1/2 a day
        if (dt_mag > 0.51) then
            day_h = daylength(lat(i,j),lon(i,j),c0,params(mp_dte)) ! daylength in h
            Llim = min(day_h,Llim*day_h)
        endif

        ! Carrying capacity
        ! ----------------------------------------------------------
          ! Gmax_crowding -> density-limited growth (ranges from the max growth rate to the death rate, excluding wave mortality)
          ! This expression follows Xiao et al (2019 and ignores wave mortality when
        ! thinking about the death rate

        B_calc = B(i,j) !*params(mp_spp_line_sep) ! converting from g/m2 to g/m

        kcap_slope = -0.75_r8

        !A = params(mp_spp_kcap_rate)/(params(mp_spp_kcap)**(-1.44_r8))
        !gH = A*B_calc**(-1.44_r8)
        A = params(mp_spp_kcap_rate)/(params(mp_spp_kcap)**(kcap_slope))
        Gmax_crowding = A * B_calc**kcap_slope
        Gmax_crowding = max(c0,Gmax_crowding)
        Gmax_crowding = min(c1,Gmax_crowding)


        ! Mortality
        ! ----------------------------------------------------------
        ! d_wave = frond loss due to waves; dependent on Hs, significant
        ! wave height [m]; Rodrigues et al. 2013 demonstrates linear relationship
        ! between Hs and frond loss rate in Macrocystis [d-1] (continuous)
        ! Duarte and Ferreira (1993) find a linear relationship between wave power and mortality in red seaweed.
        if (params(mp_breakage_type) == breakage_Duarte_Ferreira) then
          !WP = rho.*g.^2/(64*pi)*swh.^2.*Tw
          WP = 1025._r8*9.8_r8**2 / (64._r8*pi)*swh(i,j)**2 * mwp(i,j) /1.e3_r8 ! [kW]
          ! [Duarte and Ferreira (1993), in daily percentage]
          M_wave = 2.3_r8*1e-4_r8*WP + 2.2_r8*1e-3_r8 * params(mp_wave_mort_factor)
        else !if (params(mp_breakage_type) == breakage_Rodrigues) then
          M_wave  = params(mp_spp_death) * swh(i,j) * params(mp_wave_mort_factor)
        endif

        ! protect against bad values here -- this should make nans zero?:
        is_valid = valid_float(M_wave)
        if (is_valid) then
            M_wave = max(0.0,M_wave) ! returns 0 if M_wave=NAN (does it in fortran?)
            M_wave = min(0.999,M_wave)
        else
            M_wave = c0
        endif

        ! M = M_wave + general Mortality rate; [d-1]
        M = params(mp_spp_death) + M_wave


        ! Forward Euler growth
        ! ----------------------------------------------------------

        ! minimum limiting nutrient
        Qlim = min(gQ,gQP,gQFe)

        B_calc = B(i,j)
        if (growth_lim_type == 0) then
            G_calc = min(params(mp_spp_Gmax_cap),Gmax_crowding) * Tlim * Llim * Qlim
        else
            G_calc = min(params(mp_spp_Gmax_cap),Gmax_crowding) * Tlim * min(Llim,Qlim)
        endif
        dBdt = B_calc * (G_calc - M) * dt_mag

        GRate(i,j) = G_calc
        Growth = G_calc * B_calc * dt_mag
        Growth2(i,j) = Growth2(i,j) + Growth

        d_B(i,j) = dBdt
        d_Bm(i,j) = d_Bm(i,j) + B_calc * M * dt_mag
        Bm_loss = B_calc * params(mp_spp_death) * dt_mag
        Bm_wave = B_calc * M_wave * dt_mag
        d_Bm_wave(i,j) = d_Bm_wave(i,j) + Bm_wave

        ! debug terms
        !UptakeN = Upc1*(1.0 - 1.0/(1.0 + (Qmax - N/B) * Upc2))
        !dQdt = Q * (1.0 / (1.0 + G * dt) - 1.0) + (UptakeN - (Q - Qmin) * E) * dt
        !alt_dQdt = (N+dNdt)/(B+dBdt)-Q #dQdt

        ! Particulate accounting - someday maybe pass back POM with molar fractions instead?
        Bm_POM = Bm_loss * mp(mp_spp_POM_from_death) + Bm_wave * mp(mp_spp_POM_from_breakage)
        d_POC(i,j) = d_POC(i,j) + Bm_POM * mp(mp_spp_BtoC) / mw_c ! mmol/m3
        d_PON(i,j) = d_PON(i,j) + Bm_POM * Q(i,j) / mw_n ! mmol/m3

        ! Exudation
        ! Exudation could be modified from params(mp_spp_DOM_Growth_ratio) to be a function of Q + light
        exude = params(mp_spp_E) + G_calc * params(mp_spp_E_Growth_ratio) * params(mp_spp_E_N_fraction) ! fractional N-loss (1/day)
        N0 = B(i,j)*Q(i,j)
        N_exude = (N0 - Qmin*B_calc)*exude*dt_mag ! mg N/m3
        C_exude = N_exude * mp(mp_spp_E_C_fraction) / mp(mp_spp_E_N_fraction)! mg C/m3
        d_Be(i,j) = d_Be(i,j) + N_exude

        ! Dissolved organics accounting
        Bm_DOM = Bm_loss * (c1 - mp(mp_spp_POM_from_death)) + d_Bm_wave(i,j) * (c1 - mp(mp_spp_POM_from_breakage))
        Bm_remin = mp(mp_spp_labile_ratio) * Bm_DOM ! instant remin of labile DOM
        Bm_DOM = Bm_DOM - Bm_remin
        d_DON(i,j) = d_DON(i,j) + ( Bm_DOM*Q(i,j) + &
                                    N_exude * (c1-mp(mp_spp_labile_ratio)) &
                                  ) / mw_n ! mg/m3 -> mmol/m3
        d_DOC(i,j) = d_DOC(i,j) + ( Bm_DOM * mp(mp_spp_BtoC) + &
                                    C_exude * (c1-mp(mp_spp_labile_ratio)) &
                                  ) / mw_c ! mmol C/m3

        ! DIC & O2 accounting
        d_DIC(i,j) = d_DIC(i,j) + ( ((Bm_remin - Growth) * mp(mp_spp_BtoC)) + &
                                     C_exude * mp(mp_spp_labile_ratio) &
                                  ) / mw_c ! mmol C/m3
        d_O2(i,j) = d_O2(i,j) + ( (Growth - Bm_remin) * mp(mp_spp_BtoC) &
                                   - C_exude * mp(mp_spp_labile_ratio) &
                                ) / mw_c * mp(mp_spp_CtoO) / c2 ! mmol O2/m3


        ! N output terms
        ! ----------------------------------------------------------
        dNdt = d_N_internal_new_vQ(N0, B(i,j), Qmin, Qmax, UptakeN, vQ_N, M, exude, dt_mag, 1.0_r8)
        Q_new = (N0+dNdt) / (B(i,j)+dBdt)
        d_Q(i,j) = Q_new - Q(i,j)
        B_N(i,j) = 1.0e3_r8/Q_new

        if (use_NH4==1) then
          d0_nh4 =  ( N_exude * mp(mp_spp_labile_ratio) & ! labile remin to NH4
                      - Uptake_NH4 * dt_mag * B_calc / params(mp_spp_dry_sa) & ! Uptake_NH4 is per surface area, have to convert
                    ) / mw_n
          d0_no3 = - Uptake_NO3 * dt_mag * B_calc / params(mp_spp_dry_sa) / mw_n  ! Uptake_NO3 is per surface area, have to convert - confusing, change?
          d_NH4(i,j) = d_NH4(i,j) + d0_nh4  ! external record
          nh4(i,j) = nh4(i,j) + d0_nh4      ! update forcing for potential sub-stepping

        else

          d0_no3 = ( N_exude * mp(mp_spp_labile_ratio) & ! labile remin to NO3
                     - Uptake_NO3 * dt_mag * B_calc / params(mp_spp_dry_sa) & ! Uptake_NO3 is per surface area, have to convert - confusing, change?
                   ) / mw_n
        endif
        d_NO3(i,j) = d_NO3(i,j) + d0_no3  ! external record
        no3(i,j) = no3(i,j) + d0_no3      ! update forcing for potential sub-stepping


        ! P output terms
        ! ----------------------------------------------------------
        ! this setup, with d_N_internal_new_vQ and exudation based on N might cause QP to drop below Qmin_P a little?
        dPdt = c0
        if (use_P == 1) then
          N0 = B(i,j)*QP(i,j)
          e_scale = mp(mp_spp_E_P_fraction) / mp(mp_spp_E_N_fraction)
          dPdt = d_N_internal_new_vQ(N0, B(i,j), Qmin_P, Qmax_P, UptakeP, vQ_P, M, exude, dt_mag, e_scale)
          Q_new = (N0+dPdt) / (B(i,j)+dBdt)

          !d_QP(i,j) = Q_new - QP(i,j)
          d_POP(i,j) = d_POP(i,j) + Bm_POM * QP(i,j) / mw_p ! mmol/m3
          d_DOP(i,j) = d_DOP(i,j) + ( Bm_DOM * QP(i,j)  &
                                      + N_exude * e_scale * (c1 - mp(mp_spp_labile_ratio)) &
                                    ) / mw_p

          d0_po4 =  ( N_exude * e_scale * mp(mp_spp_labile_ratio) &
                      - UptakeP * dt_mag * B_calc &
                    ) / mw_p  ! mmol P per m3
          d_PO4(i,j) = d_PO4(i,j) + d0_po4  ! external record
          po4(i,j) = po4(i,j) + d0_po4      ! update forcing for potential sub-stepping

          QP(i,j) = Q_new
        endif

        ! Fe output terms
        ! ----------------------------------------------------------
        ! this setup, with d_N_internal_new_vQ and exudation based on N might cause QFe to drop below Qmin_Fe a little?
        dFedt = c0
        if (use_Fe == 1) then
          N0 = B(i,j)*QFe(i,j)
          e_scale = mp(mp_spp_E_Fe_fraction) / mp(mp_spp_E_N_fraction)
          dFedt = d_N_internal_new_vQ(N0, B(i,j), Qmin_Fe, Qmax_Fe, UptakeFe, vQ_Fe, M, exude, dt_mag, e_scale)
          Q_new = (N0+dFedt) / (B(i,j)+dBdt)

          !d_QFe(i,j) = Q_new - QFe(i,j)
          d_POFe(i,j) = d_POFe(i,j) + Bm_POM * QFe(i,j) / mw_Fe ! mmol/m3
          d_DOFe(i,j) = d_DOFe(i,j) + ( Bm_DOM * QFe(i,j)  & ! mg Fe/m3
                                      + N_exude * e_scale * (c1 - mp(mp_spp_labile_ratio)) & ! mg Fe/m3
                                    ) / mw_Fe ! mmol Fe/m3

          d0_fe = ( N_exude * e_scale * mp(mp_spp_labile_ratio) &
                    - UptakeFe * dt_mag * B_calc &
                  ) / mw_Fe  ! mmol P per m3
          d_Fe(i,j) = d_Fe(i,j) + d0_fe  ! external record
          fe(i,j) = fe(i,j) + d0_fe      ! update forcing for potential sub-stepping

          QFe(i,j) = Q_new
        endif


        ! Additional Output terms
        ! ----------------------------------------------------------

        !print *,day_h,Growth,gQ,gT ,gE , gH

        gQout(i,j) = Qlim
        gTout(i,j) = Tlim
        gEout(i,j) = Llim
        gHout(i,j) = Gmax_crowding

        !           1, 2, 3, 4
        mlim = min(Qlim,Tlim,Llim,Gmax_crowding)
        if (Tlim == mlim) then  ! depending on floating point type, this may not work in fortran
            min_lim(i,j) = 2
        elseif (Llim == mlim) then
            min_lim(i,j) = 3
        elseif (Qlim == mlim) then
            min_lim(i,j) = 1
        elseif (Gmax_crowding == mlim) then
            min_lim(i,j) = 4
        endif


        ! Update State Variables - other state variables updated earlier
        ! ----------------------------------------------------------
        B(i,j) = B(i,j) + d_B(i,j)
        Q(i,j) = Q(i,j) + d_Q(i,j)


        ! Harvest
        ! ----------------------------------------------------------

        ! increment growth/death running averages
        if (seed_now(i,j) > 0) then
           Gave(i,j) = Growth ! growth [1/timestep]
           Dave(i,j) = M    ! death [1/timestep]
        endif
        Gave(i,j) = Gave(i,j) + (Growth-Gave(i,j)) / (params(mp_harvest_avg_period) / dt_mag )
        Dave(i,j) = Dave(i,j) + (M-Dave(i,j)) / (params(mp_harvest_avg_period) / dt_mag )

        seed_now(i,j) = 0

        harv1 = c0 ! init this to zero harvest (fraction)

        if (turn_of_day) then

            ! counter for if death rate exceeds growth rate
            if (B_0000(i,j) > B(i,j)) then
                GD_count(i,j) = GD_count(i,j) + 1
            else
                GD_count(i,j) = 0
            endif
            B_0000(i,j) = B(i,j)

            !print *, do_harvest(i,j)
            !print *, n_harv(i,j),t_harv(i,j),params(mp_harvest_nmax)
            !print *, B(i,j),params(mp_harvest_kg)*params(mp_spp_line_sep)*1.0e3_r8
            !print *, GD_count(i,j),params(mp_harvest_avg_period)

            ! check for, and perform harvest

            if (params(mp_harvest_schedule) == 0) then
                ! fixed harvest
                if (do_harvest(i,j) == 1) then
                    if (params(mp_harvest_type) == 0) then
                        ! harvest to seed weight
                        if (B(i,j) > c0) then
                            harv1 = max(c0,(B(i,j)-params(mp_spp_seed)/params(mp_spp_line_sep))/B(i,j))
                        endif
                        n_harv(i,j) = n_harv(i,j) + 1
                        t_harv(i,j) = t_harv(i,j) + 1

                    else ! params[mp_harvest_type] == 1:
                        ! fractionally harvest, but not if below mp_harvest_kg
                        !if (B(i,j) >= params(mp_harvest_kg)*params(mp_spp_line_sep)*1.0e3_r8) then
                        !    harv1 = c0
                        !else
                            harv1 = params(mp_harvest_f)
                            n_harv(i,j) = n_harv(i,j) + 1
                            t_harv(i,j) = t_harv(i,j) + 1
                        !endif

                    endif
                endif
            else
                !  conditional harvest
                if (do_harvest(i,j) == 2) then
                    if (params(mp_harvest_type) == 1) then
                        ! we are done for the season, if not already harvested
                        harv1 = params(mp_harvest_f)
                    else
                        harv1 = 0.99_r8
                    endif
                    mask(i,j) = 0
                    n_harv(i,j) = n_harv(i,j) + 1
                    t_harv(i,j) = t_harv(i,j) + 1
                elseif (do_harvest(i,j) == 1) then

                    if (params(mp_harvest_schedule) == 1) then

                        ! within harvest span - check for declining growth
                        !if Gave[i,j]/Dave[i,j] < 1.0:  ! test of declining biomass over time
                        if (GD_count(i,j) > params(mp_harvest_avg_period)) then ! test of negative growth over time
                            if (params(mp_harvest_type) == 1) then
                                ! we are done for the season, if not already harvested
                                harv1 = params(mp_harvest_f)
                            else
                                harv1 = 0.99_r8
                                !harv1 = max(c0,(B(i,j)-params(mp_spp_seed)/params(mp_spp_line_sep))/B(i,j))
                            endif
                            mask(i,j) = 0
                            n_harv(i,j) = n_harv(i,j) + 1
                            t_harv(i,j) = t_harv(i,j) + 1
                        endif

                        if (harv1 < 1.0e-4_r8) then
                            if (t_harv(i,j) < (params(mp_harvest_nmax)-1)) then   !!!!!!!!!!!!!!!!! add -1
                                ! if not already harvested, check to see if conditions are OK for incremental harvest
                                if (B(i,j) >= params(mp_harvest_kg)*params(mp_spp_line_sep)*1.0e3_r8) then
                                    harv1 = params(mp_harvest_f)
                                    n_harv(i,j) = n_harv(i,j) + 1
                                    t_harv(i,j) = t_harv(i,j) + 1
                                endif
                            endif
                        endif

                    else ! period harvest -- params[mp_harvest_schedule] == 2:

                        ! do harvest according to "type"
                        if (params(mp_harvest_type) == 0) then
                            ! harvest to seed weight
                            harv1 = max(c0,(B(i,j)-params(mp_spp_seed)/params(mp_spp_line_sep))/B(i,j))
                            n_harv(i,j) = n_harv(i,j) + 1
                            t_harv(i,j) = t_harv(i,j) + 1
                        else ! params[mp_harvest_type] == 1:
                            ! fractionally harvest, but not if below mp_harvest_kg
                            if (B(i,j) > params(mp_harvest_kg)*params(mp_spp_line_sep)*1.0e3_r8) then
                                harv1 = params(mp_harvest_f)
                                n_harv(i,j) = n_harv(i,j) + 1
                                t_harv(i,j) = t_harv(i,j) + 1
                            endif
                        endif

                    endif

                endif
            endif
        endif
        !if (harv1 > c0) then  ! reset mean biomass growth/death records so they don't swing around?
        !   Gave(i,j) = Growth ! growth [1/timestep]
        !   Dave(i,j) = M    ! death [1/timestep]
        !endif
        !print *, harv1, n_harv(i,j), t_harv(i,j)

        harv1 = harv1 * B(i,j)  ! biomass harvested
        B(i,j) = B(i,j) - harv1
        harv(i,j) = harv(i,j) + harv1

        ! ## sanity checks, at end of calc to prevent output of bad data
        if (mask(i,j) > 0) then
            bad_forcing = .false.
            ! FORTRAN test for NaN - not equal to itself!
            if (sst(i,j) /= sst(i,j)) then
              bad_forcing = .true.
            endif
            if (par(i,j) /= par(i,j)) then
              bad_forcing = .true.
            endif
            if (no3(i,j) /= no3(i,j)) then
              bad_forcing = .true.
            endif
            if (swh(i,j) /= swh(i,j)) then
              bad_forcing = .true.
            endif
            if (cmag(i,j) /= cmag(i,j)) then
              bad_forcing = .true.
            endif
            if (nflux(i,j) /= nflux(i,j)) then
              bad_forcing = .true.
            endif
            if (B(i,j) /= B(i,j)) then
              bad_forcing = .true.
            endif
            if (Q(i,j) /= Q(i,j)) then
              bad_forcing = .true.
            endif
            ! not sure what to set things too - can't use nan in fortran?  maybe using ieee stuff
            if (bad_forcing) then
              mask(i,j) = 0
              B(i,j) = -999.9
              Q(i,j) = -999.9
              harv(i,j) = c0
              d_B(i,j) = -999.9
              d_Q(i,j) = -999.9
              GRate(i,j) = -999.9
              Growth2(i,j) = c0
              d_Be(i,j) = -999.9
              d_Bm(i,j) = -999.9
              !d_Ns(i,j) = -999.9
              B_N(i,j) = -999.9
              n_harv(i,j) = 0
              Gave(i,j) = -999.9
              Dave(i,j) = -999.9
              min_lim(i,j) = 0
              gQout(i,j) = -999.9
              gTout(i,j) = -999.9
              gEout(i,j) = -999.9
              gHout(i,j) = -999.9
            endif

        endif

      endif
      enddo
    enddo

  !$OMP END DO
  !$OMP END PARALLEL


  end SUBROUTINE mag_calc


  real(kind=r8) PURE FUNCTION d_N_internal_new_vQ(N, B, Qmin, Qmax, UptakeN, vQ, D, E, dt, E_scaler)

    ! This finds the new internal nutrient concentration, in mmol (as opposed to Q units)
    ! using the new vQ formulation (as opposed to the old linear uptake velocity)
    ! E_scaler : scale factor applied to the exudation rate for a particular nutrient, since the exudate may have altered nutrient ratios compared to biomass

    use macmods_kinds_mod, only : r8
    implicit none

    real(kind=r8), INTENT(IN) :: N, B, Qmin, Qmax, UptakeN, vQ, D, E, dt, E_scaler

    ! internal variables
    !----------------------------------------------------
    real(kind=r8) :: Upc1, Upc2
    real(kind=r8), parameter :: c0 = 0.0_r8
    real(kind=r8), parameter :: c1 = 1.0_r8

    ! calculate some intermediate terms, protecting for divide by zero
    Upc1 = c0
    if (vQ > c0) then
      Upc1 = UptakeN/vQ
    endif
    Upc2 = 50.0_r8/(Qmax-Qmin)

    ! below is the forward Euler accounting for nutrient uptake, and nutrient loss
    ! due to death/poc loss (D) and exudation (E)

    d_N_internal_new_vQ = (B*Upc1*(c1 - c1/(c1 + (Qmax - N/B) * Upc2)) - N*D - (N - Qmin*B)*E*E_scaler)*dt

  END FUNCTION d_N_internal_new_vQ


  PURE SUBROUTINE growth_forward_euler(N, B, Gmax, A, Kslope, Qmin, Qmax, E, D, &
                    Llim, Tlim, growth_lim_type, Upc1, Upc2, dt, &
                    B_calc, N_calc, G_calc, dBdt, dNdt, E_scaler)
!
!     Generates the forward euler method (delta) biomass at time t+1, incorporating
!     growth rate dependencies on Q (modified Droop version) and crowding, of the
!     form growth_rate_max = A*Biomass**(Kslope).
!
!     Parameters
!     ----------
!     N : step initial Nutrient mass (mg / m-2; tracer, what we are solving for)
!     B : step initial biomass (g / m-2; tracer, what we are solving for)
!     Gmax : maximum growth rate (1/day)
!     A : crowding constant, dependent on species parameters
!     Kslope : crowding constant, dependent on species parameters
!     Qmin : minimum species Nutrient storage (mg Nut/g biomass)
!     Qmax : maximum Nutrient storage (mg Nut/g biomass)
!     E : exudation rate (1/day)
!     D : death rate (1/day)
!     Llim : light limitation, considered constant [0-1]
!     Tlim : temperature limitation, considered constant  [0-1]
!     growth_lim_type : growth rate calc switch
!     Upc1: Uptake constant 1 (environment+mechalis menton limited)
!     Upc2: Uptake constant 2 (50./(Qmax-Qmin), part of the modified Droop Q-based uptake limtation formulation
!     dt : time step (days)
!     E_scaler : scale factor applied to the exudation rate for a particular nutrient, since the exudate may have altered nutrient ratios compared to biomass
!
!     Returns
!     -------
!     B_calc : best estimate of biomass used for step calc
!     N_calc : best estimate of nutrient mass used for step calc
!     G_mid = best estimate of biological growth rate used for step calc
!     dBdt = biomass change at t+dt
!     dNdt = Q change at t+dt
!
    use macmods_kinds_mod, only : i4,r8
    implicit none        

    real(kind=r8), INTENT(IN) :: &
    N, B, Gmax, A, Kslope, Qmin, Qmax, E, D, &
    Llim, Tlim, Upc1, Upc2, dt, E_scaler

    integer(kind=i4), INTENT(IN) ::growth_lim_type

    real(kind=r8), INTENT(OUT) :: &
    B_calc, N_calc, G_calc, dBdt, dNdt

    ! internal variables
    !----------------------------------------------------
    real(kind=r8) :: &
    Gmax_crowding, Q, Qlim

    Gmax_crowding = A * B**Kslope
    Gmax_crowding = max(0.0_r8,Gmax_crowding)
    Gmax_crowding = min(1.0_r8,Gmax_crowding)

    ! Calculate G and Qlim based on current values
    Q = N/B
    Qlim = (Q - Qmin) / Q * Qmax/(Qmax-Qmin)
    if (growth_lim_type == 0) then
        G_calc = min(Gmax,Gmax_crowding) * Tlim * Llim * Qlim
    else
        G_calc = min(Gmax,Gmax_crowding) * Tlim * min(Llim,Qlim)
    endif

    ! debug terms
    !UptakeN = Upc1*(1.0 - 1.0/(1.0 + (Qmax - N/B) * Upc2))
    !dQdt = Q * (1.0 / (1.0 + G * dt) - 1.0) + (UptakeN - (Q - Qmin) * E) * dt
    !alt_dQdt = (N+dNdt)/(B+dBdt)-Q #dQdt

    dBdt = B * (G_calc - D) * dt
    dNdt = (B*Upc1*(1.0_r8 - 1.0_r8/(1.0_r8 + (Qmax - N/B) * Upc2)) - N*D - (N - Qmin*B)*E*E_scaler)*dt

    ! return B,Q,G in addition to dBdt and dQdt for calculation of some output terms
    B_calc = B
    N_calc = N

  end SUBROUTINE growth_forward_euler



  real(kind=r8) PURE FUNCTION uptake_velocity(Q, Qmin, Qmax, uptake_type)

    use macmods_kinds_mod, only : r8,i4
    implicit none

    ! internal variables
    !----------------------------------------------------
    real(kind=r8), parameter :: c0 = 0.0_r8
    real(kind=r8), parameter :: c1 = 1.0_r8

    real(kind=r8), INTENT(IN) :: Qmax, Qmin, Q
    integer(kind=i4), INTENT(IN) :: uptake_type

    !print *,j,lambda
    ! Quota-limited uptake: maximum uptake when Q is minimum and
    ! approaches zero as Q increases towards maximum; Possible that Q
    ! is greater than Qmax. Set any negative values to zero.
    if (uptake_type == 0) then
        uptake_velocity = (Qmax-Q)/(Qmax-Qmin)

    ! new vQ formulation that permits high uptake even with high Q, so there is not negative feedback
    ! between high growth rate/high Q usage and reduced growth rate
    ! The reasoning behind this uptake curve is that seaweed would not reduce uptake under high-growth,
    ! high-nutrient conditions, just because stores are full - it's really a time stepping issue. High Q->
    ! leads to nutrient limitation in one time-step.  If at moderate Q, seaweeds are also prevented from
    ! uptaking short-pulsed nutrients, like closer to the scale of the timestep.  Basically the linear
    ! function above does not allow response to changing nutrients at rates that seaweeds are
    ! capable of.
    elseif (uptake_type == 1) then
        uptake_velocity = c1 - c1/(c1+(max(c0,Qmax-Q))*50.0_r8/(Qmax-Qmin))
    endif

    uptake_velocity = max(c0,uptake_velocity)
    uptake_velocity = min(c1,uptake_velocity)

  end FUNCTION uptake_velocity




  real(kind=r8) PURE FUNCTION nutrient_uptake(cmag,mwp,CD,Vmax,Ks,Nconc,vQ)

    use macmods_kinds_mod, only : r8
    implicit none

    real(kind=r8), INTENT(IN) :: cmag, mwp, CD, Vmax, Ks, Nconc, vQ

    ! internal variables
    !----------------------------------------------------
    real(kind=r8), parameter :: c0 = 0.0_r8
    real(kind=r8), parameter :: c1 = 1.0_r8

    real(kind=r8) :: lambda, N_u, vNuTw

    lambda = lambda_NO3(cmag, mwp, CD, Vmax, Ks, Nconc)

    ! Below is what we call "Uptake Factor." It varies betwen 0
    ! and 1 and includes kinetically limited uptake and
    ! mass-transfer-limited uptake (oscillatory + uni-directional flow)
    N_u = Nconc*1000.0_r8
    vNuTw = N_u / (Ks * ((N_u/Ks)  + 0.5_r8 * (lambda+sqrt(lambda**2 + 4.0_r8 * (N_u/Ks)))))
    vNuTw = max(c0,vNuTw)
    vNuTw = min(c1,vNuTw)

    ! Uptake Rate [mg N/g(dry)/d]
    ! Nutrient Uptake Rate = Max Uptake * v[Ci,u,Tw] * vQ
    ! converted from umol N/m2/d -> mg N/g(dry)/d by 14.0067 / 1e3
    nutrient_uptake = Vmax * vNuTw * vQ ! [umol/m2/d]


  end FUNCTION nutrient_uptake


  real(kind=r8) PURE FUNCTION Q_limitation(Q,Qmax,Qmin,Q_lim_type)

    use macmods_kinds_mod, only : r8,i4
    implicit none

    real(kind=r8), INTENT(IN) :: Q,Qmax,Qmin
    integer(kind=i4), INTENT(IN) :: Q_lim_type


    if (Q_lim_type == 0) then
        Q_limitation = (Q - Qmin) / Q * Qmax/(Qmax-Qmin) ! Droop scaled from 0-1
    elseif (Q_lim_type == 1) then
        Q_limitation = (Q - Qmin) / (Qmax-Qmin) ! Freider et al.
    elseif (Q_lim_type == 2) then
        Q_limitation = 1.0_r8 - exp(-11.0_r8 * (Q - Qmin) / (Qmax-Qmin)) ! more permissive than Droop scaled
    endif

    Q_limitation = max(0.0_r8,Q_limitation)
    Q_limitation = min(1.0_r8,Q_limitation)

  end FUNCTION Q_limitation

  ! a wrapper that references macmods output variables from block storage and
  ! calls mag_calc
  SUBROUTINE mag_calc_block(lat,lon,sst,par,chl,swh,mwp,cmag,nflux, &
    no3,nh4,po4,fe, &
    do_harvest,seed_now, mask, GD_count, B_0000, &
    params, &
    MO, &
    turn_of_day,dt_mag)

    implicit none

    real(kind=r8), INTENT(IN), DIMENSION(imt,jmt) :: &
      lat, &
      lon, &
      sst, &
      par, &
      chl, &
      swh, &
      mwp, &
      cmag, &
      nflux

    real(kind=r8), INTENT(INOUT), DIMENSION(imt,jmt) :: &
      no3, &
      nh4, &
      po4, &
      fe, &
      B_0000

    integer(kind=i4), INTENT(IN), DIMENSION(imt,jmt) :: &
      do_harvest

    integer(kind=i4), INTENT(INOUT), DIMENSION(imt,jmt) :: &
      seed_now, &
      mask, &
      GD_count

    real(kind=r8), INTENT(INOUT), DIMENSION(imt,jmt,n_outputs) :: &
        MO

    real(kind=r8), INTENT(IN) :: &
      params(npar), &
      dt_mag

    logical, INTENT(IN) :: &
      turn_of_day

     CALL mag_calc(lat,lon,sst,par,chl,swh,mwp,cmag,nflux, &
        no3,nh4,po4,fe, &
        do_harvest,seed_now, mask, GD_count, B_0000, MO(:,:,mo_n_harv), MO(:,:,mo_t_harv), &
        params, &
        MO(:,:,mo_QN),MO(:,:,mo_QP),MO(:,:,mo_QFe),MO(:,:,mo_B),MO(:,:,mo_d_B), &
        MO(:,:,mo_d_QN),MO(:,:,mo_Growth),MO(:,:,mo_d_Be),MO(:,:,mo_d_Bm), &
        MO(:,:,mo_d_Bm_wave),MO(:,:,mo_harv),MO(:,:,mo_GRate),MO(:,:,mo_B_N), &
        MO(:,:,mo_Gave),MO(:,:,mo_Dave),MO(:,:,mo_d_NO3),MO(:,:,mo_d_NH4), &
        MO(:,:,mo_d_PO4),MO(:,:,mo_d_Fe),MO(:,:,mo_d_DIC),MO(:,:,mo_d_O2), &
        MO(:,:,mo_d_DOC),MO(:,:,mo_d_DON),MO(:,:,mo_d_DOP),MO(:,:,mo_d_DOFe), &
        MO(:,:,mo_d_POC),MO(:,:,mo_d_PON),MO(:,:,mo_d_POP),MO(:,:,mo_d_POFe), &
        MO(:,:,mo_min_lim),MO(:,:,mo_gQ),MO(:,:,mo_gT),MO(:,:,mo_gE),MO(:,:,mo_gH), &
        turn_of_day,dt_mag)

end SUBROUTINE mag_calc_block



end MODULE macmods_calc




