module macmods_mod

    !use marbl_kinds_mod, only : int_kind
    !use marbl_kinds_mod, only : r8
    use macmods_kinds_mod, only : r8,i4,log_kind
    use macmods_param_mod
    use macmods_util_mod
    use macmods_calc, only: mag_calc_block

    implicit none

    Private

    ! The general plan for creating a driver for MACMODS:
    ! 
    ! 1) Build the "macmods_parameters.txt" file (see necessary param order in the parameters module)
    ! 2) call macmods_init(nCalcSteps,dtDays)
    ! 3) call macmods_set_forcing() on static data like lat/lon arrays, etc., if they don't change
    ! 4) Inside the model time step loop:
    !       call macmods_set_forcing() on everything needed for the time step
    !       call macmods_tendency_compute()
    !       use the output tendencies to update model state variables
    !       write out public state a output variables in whatever way you want


    ! We allocate storage for the whole grid, assuming memory is available,
    ! then restrict usage of particular grid cells in dynamically input
    ! parameters.
    ! --------------------------------------------------------------------


    ! public tracer and output storage

    ! Now all outputs/tracers are in a block of reals.  This makes it easier to
    ! link libraries without having all these variables in the driver namespace.
    real(kind=r8), dimension(imt,jmt,n_outputs) :: MO


!     real(kind=r8), dimension(imt,jmt), public :: & ! (currently only support one active depth, no km dimension)
!       B,      & ! macroalgae biomass [g/m2] dry weight, currently fixed biomass:C ratio
!       QN,     & ! nitrogen quotient [mg N/g B]
!       QP,     & ! phosphorus quotient [mg P/g B]
!       QFe,    & ! iron quotient [mg P/g B]
!       Gave,   & ! average growth rate, sort of dynamic and used for harvest or senescence calcs [1/d]
!       Dave      ! average death rate, sort of dynamic and used for harvest or senescence calcs [1/d]
!
!
!     ! public output storage, could be moved to elsewhere
!     real(kind=r8), dimension(imt,jmt), public :: & ! (currently only support one active depth, no km dimension)
!       Growth,  & ! macroalgae biomass [g/m2] dry weight, currently fixed biomass:C ratio
!       d_Be,    & ! nitrogen lost to exudation [mg N/m2], cumulative since last write
!       d_Bm,    & ! biomass that died [g/m2] dry weight, cumulative since last write
!       d_Bm_wave,    & ! biomass that broke off [g/m2] dry weight, part of d_Bm, cumulative since last write
!       d_B,   & ! change in B [g/m2] dry weight, instantaneous for subtimestep
!       d_Q,   & ! change in nitrogen quotient [mg N/g B], instantaneous for subtimestep
!       d_NO3,   & ! seaweed update of NO3 [mmol/m2], cumulative since last write
!       d_NH4,   & ! seaweed update of NH4 [mmol/m2], cumulative since last write
!       d_PO4,   & ! seaweed update of P [mmol/m2], cumulative since last write
!       d_Fe,    & ! seaweed update of Fe [mmol/m2], cumulative since last write
!       d_DIC,   & ! seaweed update of DIC [mmol/m2], cumulative since last write
!       d_O2,    & ! seaweed update of dissolved O2 [mmol/m2], cumulative since last write
!       d_DOC,   & ! seaweed contribution to DOC pool [mmol/m3]
!       d_DON,   & ! seaweed contribution to DON pool [mmol/m3]
!       d_DOP,   & ! seaweed contribution to DOP pool [mmol/m3]
!       d_DOFe,  & ! seaweed contribution to DOFe pool [mmol/m3]
!       d_POC,   & ! seaweed contribution to POC pool [mmol/m3]
!       d_PON,   & ! seaweed contribution to PON pool [mmol/m3]
!       d_POP,   & ! seaweed contribution to POP pool [mmol/m3]
!       d_POFe,  & ! seaweed contribution to POFe pool [mmol/m3]
!       harv,    & ! harvest biomass [g/m2], cumulative since last write
!       Grate,   & ! average growth rate, sort of dynamic and used for harvest or senescence calcs [1/d]
!       B_N,     & ! hmm...
!
!       gQ,      & ! nutrient limitation term [0-1]
!       gT,      & ! temperature limitation [0-1]
!       gE,      & ! light limitation [0-1]
!       gH         ! crowding limitation [0-1]
!
!     integer(kind=i4), dimension(imt,jmt) :: &
!       n_harv,  & ! number of harvests since last output/write [#]
!       t_harv,  & ! number of harvests since seeding [#]
!       min_lim    ! miniumum limitation term on growth [0-1]    seed_now(nx1,ny1), &

    ! public forcing storage, could be moved to elsewhere
    real(kind=r8), dimension(imt,jmt), public :: &
      lat,  &  ! degrees latitude
      lon,  &  ! degrees longitude
      sst,  &  ! sea (surface) temperature - legacy name, will be at whatever depth assigned
      par,  &  ! photosyntheticaly active radiation (W)
      chl,  &  ! phytoplankton chloropyll-a, used for attenuating light at depth, if enabled
      swh,  &  ! swell height (m)
      mwp,  &  ! mean wave period (s)
      cmag,  &  ! current magnitude (m/s)
      nflux, &  ! maximum vertical NO3-flux at 100m, if used
      no3,  &  ! ambient Nitrate (mmol/m3)
      nh4,  &  ! ambient ammonia (mml/m3)
      po4,  &  ! ambient photphate (mmol/m3)
      fe       ! ambient bioavailable iron (mmol/m3)

    ! internal data arrays
    real(kind=r8), dimension(imt,jmt) :: &
      doy_seeding,  & ! day-of-year (Jan 1 = 1) date of seeding
      doy_diff,     & ! nitrogen lost to exudation [mg N/m2], cumulative since last write
      doy_mod

    integer(kind=i4), dimension(imt,jmt) :: &
      do_harvest, & ! stage of the harvest cycle(s)
      seed,       & ! (re)seed flag
      seed_month, & ! optional array specifying which month of the year to seed in
      mask,       & ! active cell flag
      GD_count      ! increment indicating consecutive days of no-growth


    real(kind=r8), dimension(imt,jmt) :: &
      B_0000      ! biomass from previous day [g/m2] dry weight


    ! FORTRAN-compatibility data arrays
    integer(kind=i4), dimension(imt,jmt) :: &
      days_since_seeding ! integer days since seeding/last harvest


    ! public tendencies, could be moved to elsewhere
    real(kind=r8), dimension(imt,jmt) :: &
      tend_DIC, &
      tend_NH4, &
      tend_NO3, &
      tend_PO4, &
      tend_Fe, &
      tend_DOC, &
      tend_DON, &
      tend_DOP, &
      tend_DOFe, &
      tend_POC, &
      tend_PON, &
      tend_POP, &
      tend_POFe, &
      tend_O2, &
      absorp  ! output extra shortwave absorption



    integer(kind=i4) :: &
      use_NH4,   &   ! Is a parameter but an integer copy is helpful
      use_PO4,   &   ! Is a parameter but an integer copy is helpful
      use_Fe,    &   ! Is a parameter but an integer copy is helpful
      mag_idx        ! depth index into forcing/external model data

    ! timing variables
    integer(kind=i4) :: calc_steps

    real(kind=r8), ALLOCATABLE :: step_doy(:),step_year(:)

    PUBLIC :: macmods_get_output_1, macmods_set_tracer_name, macmods_set_tracer_1_name, &
      macmods_set_tracer_1, macmods_get_output_block, macmods_get_output_block_1, &
      macmods_zero_outputs, macmods_init, macmods_get_depth, macmods_set_forcing, &
      macmods_tendency_compute

  CONTAINS

    !---------------------------------------------------------------------
    !  Dynamically set all tracers by name.  I wonder if we should use masking...
    !---------------------------------------------------------------------
     SUBROUTINE macmods_set_tracer_name(tracer_name,tracer_values)

      CHARACTER(LEN=*), INTENT(IN) :: tracer_name
      REAL(KIND=r8), INTENT(IN), DIMENSION(IMT,JMT) :: tracer_values

      ! There are only a few tracers, we don't want to be able to set
      ! output variables in the block, so just list them here
      SELECT CASE (TRIM(tracer_name))
          CASE ('B')
              MO(:,:,mo_B) = tracer_values
          CASE ('QN')
              MO(:,:,mo_QN) = tracer_values
          CASE ('QP')
              MO(:,:,mo_QP) = tracer_values
          CASE ('QFe')
              MO(:,:,mo_QFe) = tracer_values
          CASE DEFAULT
            print *,'set_tracer: tracer name not known: '//tracer_name
            call exit(-1)
      END SELECT

    END SUBROUTINE macmods_set_tracer_name


    !---------------------------------------------------------------------
    !  Dynamically set one tracer by name
    !---------------------------------------------------------------------
    SUBROUTINE macmods_set_tracer_1_name(tracer_name,tracer_value,i,j)

      CHARACTER(LEN=*), INTENT(IN) :: tracer_name
      REAL(KIND=r8), INTENT(IN) :: tracer_value
      INTEGER(KIND=i4), INTENT(IN) :: i,j

      ! There are only a few tracers, we don't want to be able to set
      ! output variables in the block, so just list them here
      SELECT CASE (TRIM(tracer_name))
          CASE ('B')
              MO(i,j,mo_B) = tracer_value
          CASE ('QN')
              MO(i,j,mo_QN) = tracer_value
          CASE ('QP')
              MO(i,j,mo_QP) = tracer_value
          CASE ('QFe')
              MO(i,j,mo_QFe) = tracer_value
          CASE DEFAULT
            print *,'set_tracer_1: tracer name not known: '//tracer_name
            call exit(-1)
      END SELECT

    END SUBROUTINE macmods_set_tracer_1_name

    !---------------------------------------------------------------------
    !  Dynamically set one tracer by index
    !---------------------------------------------------------------------
    SUBROUTINE macmods_set_tracer_1(tracer_idx,tracer_value,i,j)

      INTEGER(KIND=i4), INTENT(IN) :: tracer_idx
      REAL(KIND=r8), INTENT(IN) :: tracer_value
      INTEGER(KIND=i4), INTENT(IN) :: i,j

      ! There are only a few tracers, we don't want to be able to set
      ! output variables in the block, so just list them here
      IF (tracer_idx <= n_tracers) then
          MO(i,j,tracer_idx) = tracer_value
      ELSE
          print *,'macmods_set_tracer_1: tracer_idx not known: ',tracer_idx
          call exit(-1)
      END IF

    END SUBROUTINE macmods_set_tracer_1

    !---------------------------------------------------------------------
    !  Get the whole output block
    !---------------------------------------------------------------------
    SUBROUTINE macmods_get_output_block(mo_out)
      REAL(KIND=r8), INTENT(OUT), DIMENSION(IMT,JMT,n_outputs) :: mo_out

      mo_out = mo(:,:,:)

    END SUBROUTINE macmods_get_output_block

    !---------------------------------------------------------------------
    !  Get the single location output block
    !---------------------------------------------------------------------
    SUBROUTINE macmods_get_output_block_1(mo_out,i,j)
      REAL(KIND=r8), INTENT(OUT), DIMENSION(n_outputs) :: mo_out
      INTEGER(KIND=i4), INTENT(IN) :: i,j

      mo_out = mo(i,j,:)

    END SUBROUTINE macmods_get_output_block_1

    !---------------------------------------------------------------------
    !  Get a single output variable by name
    !---------------------------------------------------------------------
    SUBROUTINE macmods_get_output(mo_name,mo_out)

      CHARACTER(LEN=*), INTENT(IN) :: mo_name
      REAL(KIND=r8), INTENT(OUT), DIMENSION(IMT,JMT) :: mo_out
      INTEGER(KIND=i4) :: mo_idx

      mo_idx = get_mo_idx(mo_name)
      mo_out = mo(:,:,mo_idx)

    END SUBROUTINE macmods_get_output


    !---------------------------------------------------------------------
    !  Get a single output variable by name and i,j index location - function
    !---------------------------------------------------------------------
    FUNCTION macmods_get_output_1(mo_idx,i,j) RESULT(mo_out)

      INTEGER(KIND=i4), INTENT(IN) :: mo_idx,i,j
      REAL(KIND=r8) :: mo_out
      
      mo_out = -c0
      IF (mo_idx >=1 .and. mo_idx <= n_outputs) then
          mo_out = mo(i,j,mo_idx)
      ELSE
        print *,'macmods_set_tracer_1: mo_idx out of range: ',mo_idx
        call exit(-1)
      END IF
    
    END FUNCTION macmods_get_output_1

    !---------------------------------------------------------------------
    !  Zero cummulcative/instantaneous outputs (e.g. after tracer 
    !  accounting and/or writing)
    !---------------------------------------------------------------------
    SUBROUTINE macmods_zero_outputs()

      mo(:,:,mo_i_cummulative:n_outputs) = c0

    END SUBROUTINE macmods_zero_outputs    

    subroutine macmods_init(read_txt_file)

      logical(kind=log_kind), intent(IN) :: read_txt_file

      ! does anything need to be passed into here from the model/driver running MACMODS, or do we let the driver take care of things?

      call macmods_load_parameters(read_txt_file)

      call macmods_init_timing() !nCalcSteps,dtDays)
      call macmods_init_memory()
      call macmods_init_seeding_harvest()

    end subroutine macmods_init



    subroutine macmods_load_parameters(read_txt_file)

      logical(kind=log_kind) :: read_txt_file

      call set_macmods_params_defaults()

      ! read a dumb, ordered list of params, perhaps?  Borrow PTM code for reading params?
      if (read_txt_file) then
        call read_macmods_params()
      endif

      use_NH4 = NINT(mp(mp_use_NH4))
      use_PO4 = NINT(mp(mp_use_P))
      use_Fe = NINT(mp(mp_use_Fe))


    end subroutine macmods_load_parameters

    ! don't seem to need this - we rely on external driver to supply the
    ! day-of-year, which is the only thing we need at this point.
    subroutine macmods_init_timing() !nCalcSteps,dtDays)

    end subroutine macmods_init_timing


    ! Don't seem to need this, memory is static with parameters for now
    subroutine macmods_init_memory()

    end subroutine

    ! init surface and other forcing data structures
    ! Don't seem to need this, memory is static with parameters for now
    subroutine macmods_init_forcing()

    end subroutine



    ! load data needed for mag_calc
    subroutine macmods_set_forcing(f_var_name,f_values)

      CHARACTER(LEN=*), INTENT(IN) :: f_var_name
      REAL(KIND=r8), DIMENSION(IMT,JMT), INTENT(IN) :: f_values

      SELECT CASE (TRIM(f_var_name))
          CASE ('lat')
            lat = f_values
          CASE ('lon')
            lon = f_values
          CASE ('sst')
            sst = f_values
          CASE ('par')
            par = f_values
          CASE ('chl')
            chl = f_values
          CASE ('swh')
            swh = f_values
          CASE ('mwp')
            mwp = f_values
          CASE ('cmag')
            cmag = f_values
          CASE ('nflux')
            nflux = f_values
          CASE ('no3')
            no3 = f_values
          CASE ('nh4')
            nh4 = f_values
          CASE ('po4')
            po4 = f_values
          CASE ('fe')
            fe = f_values

      END SELECT
    end subroutine macmods_set_forcing

    ! use run parameters to construct seeding and harvest setup
    subroutine macmods_init_seeding_harvest(init_mask)


      integer(kind=i4), dimension(imt,jmt), optional :: init_mask

      seeding_type = NINT(mp(mp_seeding_type))

      if ((seeding_type == init_seeding_type) .or. &
          (seeding_type == seed_on_harvest_type)) then

          IF (PRESENT(init_mask)) THEN
              WHERE (init_mask >= 1)
                  MO(:,:,mo_B) = mp(mp_spp_seed)
                  MO(:,:,mo_QN) = mp(mp_spp_initQ)
                  MO(:,:,mo_QP) = mp(mp_spp_initQP)
                  MO(:,:,mo_QFe) = mp(mp_spp_initQFe)
                  mask = 1
              END WHERE
          ELSE
              MO(:,:,mo_B) = mp(mp_spp_seed)
              MO(:,:,mo_QN) = mp(mp_spp_initQ)
              MO(:,:,mo_QP) = mp(mp_spp_initQP)
              MO(:,:,mo_QFe) = mp(mp_spp_initQFe)
              mask = 1
          ENDIF
      else
          ! mask everything initially, to skip computation until 1st seeding
          mask = 0
      endif

      do_harvest = 0


    end subroutine

    ! calc seeding per grid cell
    subroutine macmods_calc_seeding(doy_in)

        real(kind=r8), intent(IN) :: doy_in
        integer(kind=i4) :: month,day_of_month

!     def find_seeding(self,nstep):
!         """If using seeding file, check if month change and set seeding array & mask."""
!         if self.p['seeding_type'] == 1:
!             if self.p['seeding_type'] == self.monthly_seeding_type:    # monthly
!                 if self.forcing.day[nstep] == 1: # turn of the month
!                     print('Seeding month:',self.forcing.month[nstep])
!                     seed_mask = self.seed_mask(self.forcing.month[nstep])
!                     self.partial_seed[...] = 0
!                     self.partial_seed[seed_mask] = 1 # tell mag_calc to seed
!                     self.mask[seed_mask] = 1 # update computation mask
!                     return self.partial_seed
!
!         return self.do_not_seed


        seed = 0

        if (NINT(mp(mp_seeding_type)) == monthly_seeding_type) then
            month = get_month_from_doy(doy_in)
            day_of_month = get_dom_from_doy(doy_in,1999) ! not worrying about leap year yet
            if (day_of_month == 1) then
                WHERE (seed_month == month)
                    seed = 1
                    mask = 1
                ENDWHERE
            endif
        endif

    end subroutine

    ! calc harvest function input per grid cell
    subroutine macmods_check_harvest(doy_in)


      real(kind=r8), intent(IN) :: doy_in
      real(kind=r8) :: this_period
      integer(kind=i4) :: harv_period, harv_sched, i

!     def check_harvest(self, nstep):
!
!         # find how may days away from seeding+harvest_freq we currently are
!         # why does harvest get determined by seeding? These names could be changed for clarity
!         doy_diff = self.dt_doy[nstep] - self.seeding_doy
!         doy_diff[doy_diff < 0] = doy_diff[
!                                      doy_diff < 0] + 365  # convert so that doy_diff is positive relative to seeding
!         doy_diff2 = np.copy(doy_diff)
!         doy_diff[doy_diff == 0] = 1  # don't harvest on seed doy...
!         doy_mod = np.mod(doy_diff, self.p['mp_harvest_freq'])
!
!         self.do_harvest[...] = 0
!
!         if self.p['mp_harvest_schedule'] == 0:
!             # this could be done ahead of time ... but maybe not for multi-year runs?
!             self.do_harvest[doy_mod == 0] = 1
!
!         elif self.p['mp_harvest_schedule'] == 1:
!
!             if self.p['mp_harvest_span'] < 0:
!                 # no harvest span -- harvest anytime kg trigger is reached
!                 self.do_harvest[...] = 1
!             else:
!                 # if harvest schedule is conditional/flexible, find if we are in harvest span
!                 harvest_mask = doy_mod >= self.p['mp_harvest_freq'] - self.p['mp_harvest_span']
!                 self.do_harvest[harvest_mask] = 1  # within harvest span
!                 self.do_harvest[doy_mod == 0] = 2  # done - final harvest, cell will be masked from further calcs
!
!         elif self.p['mp_harvest_schedule'] == 2:  # fixed within period
!             harv_period = np.int32(np.float32(self.p['mp_harvest_span']) / self.p['mp_harvest_nmax'])
!             self.do_harvest[doy_mod == 0] = 2  # done - final harvest, cell will be masked from further calcs
!             for i in range(1, self.p['mp_harvest_nmax']):
!                 self.do_harvest[doy_diff2 == (self.p['mp_harvest_freq'] - harv_period * i)] = 1

      harv_sched = NINT(mp(mp_harvest_schedule)) ! conveniece variable

      doy_diff = doy_in - doy_seeding
      WHERE (doy_diff < c0) doy_diff = doy_diff + 365

      do_harvest = 0   ! init to no harvest


      if ((harv_sched==0) .or. (harv_sched==0)) then
          WHERE (doy_diff == c0) doy_diff = 1.0 ! don't harvest on seed day
          doy_mod = MOD(doy_diff,mp(mp_harvest_freq))
          call round_to_real_if_close_array2(doy_mod,0.0001_r8) ! tidy up potential slop in calc, rounding to 0.0001 of a day

          if (harv_sched == 0) then

              WHERE (doy_diff == c0) doy_diff = 1.0 ! don't harvest on seed day


              WHERE (doy_mod == c0) do_harvest = 1 ! works because doy_mod and do_harvest are same size?

          elseif (harv_sched == 1) then

              WHERE (doy_diff == c0) doy_diff = 1.0 ! don't harvest on seed day
              if (NINT(mp(mp_harvest_span)) < 0) then
                  ! no harvest span -- harvest anytime kg trigger is reached
                  do_harvest = 1
              else
                  WHERE (doy_mod >= mp(mp_harvest_freq) - mp(mp_harvest_span)) do_harvest = 1
                  WHERE (doy_mod == c0) do_harvest = 2
              endif
          endif

      elseif (harv_sched == 2) then

        doy_mod = MOD(doy_diff,mp(mp_harvest_freq))
          call round_to_real_if_close_array2(doy_mod,0.0001_r8) ! tidy up potential slop in calc, rounding to 0.0001 of a day

          harv_period = NINT(mp(mp_harvest_span) / NINT(mp(mp_harvest_nmax)))
          WHERE (doy_mod == c0) do_harvest = 2  ! done - final harvest, cell will be masked from further calcs
          do i = 1,NINT(mp(mp_harvest_nmax))
              this_period = round_to_real_if_close(mp(mp_harvest_freq) - harv_period * i, 0.0001_r8)
              WHERE(doy_diff == this_period) do_harvest = 1
          enddo

      endif

    end subroutine macmods_check_harvest

    ! ---------------------------------------------------------
    ! compute step
    ! ---------------------------------------------------------
    subroutine macmods_tendency_compute(dt_external,t_doy_in)

      real(kind=r8), intent(IN) :: dt_external,t_doy_in   ! days
      real(kind=r8) :: dt_internal, t_doy, t_doy_next ! days
      integer(kind=i4) :: nSubCalcs, nStep
      logical(kind=log_kind) :: turn_of_day

      ! zero tendencies for ecosys
!       tend_DIC = c0
!       tend_NH4 = c0
!       tend_NO3 = c0
!       tend_PO4 = c0
!       tend_Fe  = c0
!       tend_POC = c0
!       tend_DOC = c0
!       tend_DON = c0
!       tend_DOFe = c0
!       tend_O2  = c0

      call macmods_check_harvest(t_doy_in)

      call macmods_calc_seeding(t_doy_in)

      ! sub-step growth
      if (dt_external == mp(mp_dt_mag)) then
          dt_internal = dt_external
          nSubCalcs = 1
      else
          nSubCalcs = NINT(dt_external/mp(mp_dt_mag)) ! rNINT does rounding, as opposed to FLOOR
          nSubCalcs = MIN(1,nSubCalcs)
          dt_internal = dt_external/nSubCalcs ! hopefully not too far off mp(mp_dt_mag)
      endif

      t_doy = t_doy_in

      do nStep=1,nSubCalcs

        ! used to determine thing that only get checked daily
        turn_of_day = .false.
        t_doy_next = t_doy_next + dt_internal
        if (INT(t_doy_next) .ne. INT(t_doy)) then
            turn_of_day = .true.
        endif

!         call mag_calc(lat,lon,sst,par,chl,swh,mwp,cmag,nflux, &
!                         no3,nh4,po4,fe, &
!                         do_harvest,seed, mask, GD_count, B_0000, n_harv, t_harv, &
!                         mp, & ! from parameters module
!                         QN,QP,QFe,B, &
!                         d_B,d_Q,Growth,d_Be,d_Bm,d_Bm_wave,harv,GRate,B_N,Gave,Dave, &
!                         d_NO3,d_NH4,d_PO4,d_Fe, &
!                         d_DIC,d_O2, &
!                         d_DOC,d_DON,d_DOP,d_DOFe, &
!                         d_POC,d_PON,d_POP,d_POFe, &
!                         min_lim,gQ,gT,gE,gH, &
!                         turn_of_day,dt_internal)

        call mag_calc_block(lat,lon,sst,par,chl,swh,mwp,cmag,nflux, &
                        no3,nh4,po4,fe, &
                        do_harvest,seed, mask, GD_count, B_0000, &
                        mp, & ! from parameters module
                        MO, & ! tracer/output block
                        turn_of_day,dt_internal)

        t_doy = t_doy_next

      enddo

      ! update tendencies, which are different than d_Nutrient, in case we can do some internal reporting from macmods.


!       tend_DIC = tend_DIC + d_DIC
!       tend_O2 = tend_O2 + d_O2
!
!       tend_NH4 = tend_NH4 + d_NH4
!       tend_NO3 = tend_NO3 + d_NO3
!       tend_PO4 = tend_PO4 + d_PO4
!       tend_Fe = tend_Fe + d_Fe
!
!       tend_DOC = tend_DOC + d_DOC
!       tend_DON = tend_DON + d_DON
!       tend_DOP = tend_DOP + d_DOP
!       tend_DOFe = tend_DOFe + d_DOFe
!
!       tend_DOC = tend_DOC + d_POC
!       tend_DON = tend_DON + d_PON
!       tend_DOP = tend_DOP + d_POP
!       tend_DOFe = tend_DOFe + d_POFe

      ! update other outputs

      ! absorp[mag_idx] = ???


    end subroutine

    ! write out state - for diagnostics or whatever
    ! let's rely on python for this for now
    subroutine macmods_report(time)

       real(kind=r8), intent(IN) :: time   ! days

    end subroutine macmods_report

    ! helper routine to return depth (if depth cycling) given t_doy_in
    subroutine macmods_get_depth(t_doy_in,depth)

      real(kind=r8), intent(IN) :: t_doy_in   ! days
      real(kind=r8), intent(OUT) :: depth   ! m
      real(kind=r8) :: hour, t_up, t_dn ! hour of day (decimal)

      t_up = mp(mp_depth_cycle_t_photic)
      t_dn = mp(mp_depth_cycle_t_deep)

      hour = floor(t_doy_in)*24.0_r8

      depth = mp(mp_depth_cycle_z_deep)
      if (hour >= t_up) then
          if (t_dn > t_up) then
              if (hour < t_dn) then
                  depth = mp(mp_depth_cycle_z_photic)
              endif
          else
            depth = mp(mp_depth_cycle_z_photic)
          endif
      endif

    end subroutine macmods_get_depth


end module macmods_mod