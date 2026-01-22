! ======================================================================
! MACMODS Paramaters, Types, & Variables
! ======================================================================
module macmods_param_mod

  use macmods_kinds_mod, only : r8,i4,log_kind

  implicit none
  public
  save

  ! Dimensions
  integer, parameter, public ::                 &
       imt             = 1                    , &   ! grid is fixed
       jmt             = 1                    , &    ! grid is fixed
       chunk_size     = 40                           !


  !real, parameter :: d_a(7) = (/ -1.36471e-1, 4.68181e-2, 8.07004e-1, -7.45353e-3, -2.94418e-3, 3.43570e-5, 3.48658e-5/)

  !---------------------------------------------------------------------
  !  Parameters Enumeration
  !---------------------------------------------------------------------
  ! Parameters are passed in to functions as a double array, so we need to
  ! know which is which.
  ! To a new parameter (sorry, FORTRAN is hard):
  !    1) Add +1 to npar
  !    2) Add the parameter name to the enumeratation, preferably at the end
  !    3) Add the WRITE lines, in the correct order of the enum, to 'write_mp_parameters_block'
  !    4) Add a default value assigment in 'set_macmods_params_defaults'
  !    5) Add a CASE with the new param in 'get_mp_idx'

  integer(kind=i4), parameter, public :: &
      npar = 76

  integer(kind=i4), parameter, public :: &

      mp_spp_kcap            =  1  , &
      mp_spp_Gmax_cap        =  2  , &
      mp_spp_PARs            =  3  , &
      mp_spp_PARc            =  4  , &
      mp_spp_BtoSA           =  5  , &
      mp_spp_line_sep        =  6  , &
      mp_spp_kcap_rate       =  7  , &
      mp_spp_Topt1           =  8  , &
      mp_spp_K1              =  9  , &
      mp_spp_Topt2           = 10  , &
      mp_spp_K2              = 11  , &
      mp_spp_CD              = 12  , &
      mp_spp_dry_sa          = 13  , &
      mp_spp_dry_wet         = 14  , &
      mp_spp_E               = 15  , &
      mp_spp_seed            = 16  , &
      mp_spp_death           = 17  , &

      mp_spp_BtoC            = 18  , &
      mp_spp_NtoC            = 19  , &
      mp_spp_NtoP            = 20  , &
      mp_spp_NtoFe           = 21  , &
      mp_spp_CtoO            = 22  , &

      mp_spp_labile_ratio    = 23  , &
      mp_spp_spare           = 24  , &
      mp_spp_E_Growth_ratio  = 25  , &
      mp_spp_E_C_fraction    = 26  , &
      mp_spp_E_N_fraction    = 27  , &
      mp_spp_E_P_fraction    = 28  , &
      mp_spp_E_Fe_fraction   = 29  , &

      mp_spp_POM_from_breakage = 30  , &
      mp_spp_POM_from_grazing  = 31  , &
      mp_spp_POM_from_death    = 32  , &

      mp_spp_Vmax            = 33  , &
      mp_spp_Ks_NO3          = 34  , &
      mp_spp_Vmax_NH4        = 35  , &
      mp_spp_Ks_NH4          = 36  , &
      mp_spp_Q0              = 37  , &
      mp_spp_Qmax            =  38 , &
      mp_spp_Qmin            = 39  , &
      mp_spp_initQ           = 40  , &

      mp_spp_Vmax_PO4        = 41  , &
      mp_spp_Ks_PO4          = 42  , &
      mp_spp_Qmin_P          = 43  , &
      mp_spp_Qmax_P          = 44  , &
      mp_spp_initQP          = 45  , &

      mp_spp_Vmax_Fe         = 46  , &
      mp_spp_Ks_Fe           = 47  , &
      mp_spp_Qmin_Fe         = 48  , &
      mp_spp_Qmax_Fe         = 49  , &
      mp_spp_initQFe         = 50  , &

      mp_use_P           = 51  , &
      mp_use_Fe          = 52  , &
      mp_use_C           = 53  , &
      mp_use_NH4         = 54  , &

      mp_seeding_type        = 55  , &
      mp_harvest_type        = 56  , &
      mp_harvest_schedule    = 57  , &
      mp_harvest_avg_period  = 58  , &
      mp_harvest_kg          = 59  , &
      mp_harvest_f           = 60  , &
      mp_harvest_freq        = 61  , &
      mp_harvest_span        = 62  , &
      mp_harvest_nmax        = 63  , &
      mp_breakage_type       = 64  , &
      mp_dte                 = 65  , &
      mp_dt_mag              = 66  , &
      mp_N_flux_limit        = 67  , &
      mp_farm_depth          = 68  , &
      mp_wave_mort_factor    = 69  , &
      mp_N_uptake_type   = 70  , &
      mp_growth_lim_type = 71  , &
      mp_Q_lim_type      = 72  , &

      mp_depth_cycle_z_photic = 73  , &
      mp_depth_cycle_z_deep   = 74  , &
      mp_depth_cycle_t_photic = 75  , &
      mp_depth_cycle_t_deep   = 76

  ! storage for the parameters (maybe this should not be public)
  real(kind=r8), dimension(npar), public :: mp


  ! Tracers/Outputs array enumeration
  integer(kind=i4), parameter, public :: &
      n_tracers = 4, &
      n_outputs = 36, &
      mo_i_cummulative = 8 ! start of cummulative (or instantaneous) outpiut vars

  integer(kind=i4), parameter, public :: &
    mo_B         =1, & !units'#','long_name'KEI flux structure count'}, -- Tracer
    mo_QN        =2, & !units'mg N/g B','long_name'nitrogen quotient'}, -- Tracer
    mo_QP        =3, & !units'mg P/g B','long_name'phosphorus quotient'}, -- Tracer
    mo_QFe       =4, & !units'mg Fe/g B','long_name'iron quotient'}, -- Tracer
    mo_Gave      =5, & !units'1/day','long_name'average growth rate, sort of dynamic and used for harvest or senescence calcs'},
    mo_Dave      =6, & !units'1/day','long_name'average death rate, sort of dynamic and used for harvest or senescence calcs'},
    mo_t_harv    =7, & !units'#','long_name'number of harvests since seeding'},

    ! below are cummulative of instantanous and should be reset after accounting and/or writing

    mo_Growth    =8, & !units'g/m2','long_name'macroalgae biomass dry weight, currently fixed biomass:C ratio'},
    mo_n_harv    =9, & !units'#','long_name'number of harvests since last output/write'},
    mo_harv      =10, & !units'g/m2','long_name'harvest biomass, cumulative since last write'},
    mo_d_Be      =11, & !units'mg N/m2','long_name'nitrogen lost to exudation, cumulative since last write'},
    mo_d_Bm      =12, & !units'g/m2','long_name'biomass that died dry weight, cumulative since last write'},
    mo_d_Bm_wave =13, & !units'g/m2','long_name'biomass that broke off dry weight, part of d_Bm, cumulative since last write'},
    mo_d_B       =14, & !units'g/m2','long_name'change in B dry weight, instantaneous for subtimestep'},
    mo_d_QN      =15, & !units'mg N/g B','long_name'change in nitrogen quotient, instantaneous for subtimestep'},
    mo_d_NO3     =16, & !units'mmol/m2','long_name'seaweed update of NO3, cumulative since last write'},
    mo_d_NH4     =17, & !units'mmol/m2','long_name'seaweed update of NH4, cumulative since last write'},
    mo_d_PO4     =18, & !units'mmol/m2','long_name'seaweed update of PO4, cumulative since last write'},
    mo_d_Fe      =19, & !units'mmol/m2','long_name'seaweed update of Fe, cumulative since last write'},
    mo_d_DIC     =20, & !units'mmol/m2','long_name'seaweed update of DIC [mmol/m2], cumulative since last write'},
    mo_d_O2      =21, & !units'mmol/m2','long_name'seaweed update of dissolved O2, cumulative since last write'},
    mo_d_DOC     =22, & !units'mmol/m2','long_name'seaweed contribution to DOC pool'},
    mo_d_DON     =23, & !units'mmol/m2','long_name'seaweed contribution to DON pool'},
    mo_d_DOP     =24, & !units'mmol/m2','long_name'seaweed contribution to DOP pool'},
    mo_d_DOFe    =25, & !units'mmol/m2','long_name'seaweed contribution to DOFe pool'},
    mo_d_POC     =26, & !units'mmol/m2','long_name'seaweed contribution to POC pool'},
    mo_d_PON     =27, & !units'mmol/m2','long_name'seaweed contribution to PON pool'},
    mo_d_POP     =28, & !units'mmol/m2','long_name'seaweed contribution to POP pool'},
    mo_d_POFe    =29, & !units'mmol Fe/m2','long_name'seaweed contribution to DOFe pool'},
    mo_Grate     =30, & !units'1/day','long_name'average growth rate, sort of dynamic and used for harvest or senescence calcs'},
    mo_B_N       =31, & !units'?','long_name'cant remember what this is'},
    mo_gQ        =32, & !units'fractional','long_name'nutrient limitation term'},
    mo_gT        =33, & !units'fractional','long_name'temperature limitation'},
    mo_gE        =34, & !units'fractional','long_name'light limitation '},
    mo_gH        =35, & !units'fractional','long_name'crowding limitation'},
    mo_min_lim   =36    !units'factional','long_name'miniumum limitation term on growth'},


  !---------------------------------------------------------------------
  !  FORTRAN-compatibility parameters -- do we still need this?
  !---------------------------------------------------------------------
  ! These parameters help replace Python functionality
  integer(kind=i4), public :: &

    seeding_type

  CONTAINS


  !---------------------------------------------------------------------
  !  Dynamically set a parameter by name
  !---------------------------------------------------------------------
  SUBROUTINE set_mp(mp_name,mp_value)

    CHARACTER(LEN=*), INTENT(IN) :: mp_name
    REAL(KIND=r8), INTENT(IN) :: mp_value
    INTEGER(KIND=i4) :: mp_idx

    mp_idx = get_mp_idx(mp_name)
    mp(mp_idx) = mp_value
  END SUBROUTINE set_mp

  !---------------------------------------------------------------------
  !  Get a parameter value by name
  !---------------------------------------------------------------------
  FUNCTION get_mp(mp_name) RESULT(mp_value)

    CHARACTER(LEN=*), INTENT(IN) :: mp_name
    INTEGER(KIND=i4) :: mp_idx
    REAL(KIND=r8) :: mp_value

    mp_idx = get_mp_idx(mp_name)
    mp_value = mp(mp_idx)
  END FUNCTION get_mp

  !---------------------------------------------------------------------
  !  Read parameters, in exact order, from a simple flat file
  !  Only reads the first number; anything after the number is ignored
  !  such that the exact output from 'write_mp_parameters_block' can
  !  be read in to set all parameters
  !---------------------------------------------------------------------
  SUBROUTINE read_macmods_params()

    implicit none

    integer :: i

    open(unit=21, file='macmods_parameters.txt', form='FORMATTED', access='sequential')
    do i=1,npar
      read(21,*) mp(i)
    enddo

  end SUBROUTINE read_macmods_params


  !---------------------------------------------------------------------
  !  Maintain a working set of defauly params, for testing or whatever
  !---------------------------------------------------------------------
  SUBROUTINE set_macmods_params_defaults()

    implicit none

    real(kind=r8) :: spp_NtoP_mass,spp_NtoFe_mass,spp_NtoC_mass

    ! Macrocystis
    ! ------------------------------------------------------------------

    mp(mp_spp_kcap) = 2000.0_r8
    mp(mp_spp_Gmax_cap) = 0.2_r8
    mp(mp_spp_PARs) = 212.4_r8/4.57_r8
    mp(mp_spp_PARc) = 20.45_r8/4.57_r8
    mp(mp_spp_Q0) = 32.0_r8
    mp(mp_spp_BtoSA) = 1.0_r8
    mp(mp_spp_line_sep) = 1.0_r8
    mp(mp_spp_kcap_rate) = 0.05_r8
    mp(mp_spp_Topt1) = 13.0_r8 !13.0_r8
    mp(mp_spp_K1) = 0.04_r8
    mp(mp_spp_Topt2) = 18.0_r8
    mp(mp_spp_K2) = 0.05_r8
    mp(mp_spp_CD) = 0.5_r8
    mp(mp_spp_dry_sa) = 58.0_r8
    mp(mp_spp_dry_wet) = 0.094_r8
    mp(mp_spp_E) = 0.01_r8
    mp(mp_spp_seed) = 200.0_r8
    mp(mp_spp_death) = 0.01_r8

    mp(mp_spp_BtoC) = 0.3_r8

    ! these are mole ratios
    ! ----------------------------------------------------
    mp(mp_spp_NtoC) = 1.0_r8 / 16.0_r8 ! Paine et al. 2023
    mp(mp_spp_NtoP) = 16.0_r8          ! Redfield
   ! ~200 mg-Fe/Kg-Biomass = 3.56 umol/g-Biomass, Smith et al. 2010, an average of the two kelps in the paper
   ! math: 0.3 gC/g-Biomass / (12.011 gC/molC) * 1/16 molN/molC * 1e6 umol/mol = 1561 umolN/g-Biomass
   ! N:Fe = 1561 umolN/g-Biomass / 3.56 umol/g-Biomass = 438.5
    mp(mp_spp_NtoFe) = mp(mp_spp_BtoC) / 12.011_r8 * mp(mp_spp_NtoC) * 1.e6_r8 / 3.56 ! mw_c = 12.011
    mp(mp_spp_CtoO) = 1.0_r8

    ! we deal with mass a lot in macmods, handy to have these for calcs
    ! can't use these constants in this module though, b/c namespace conflict in f2py/python
    ! ----------------------------------------------------
    ! mw_c      = 12.011_r8,  & ! molecular weight C
    ! mw_n      = 14.00672_r8,  & ! molecular weight N
    ! mw_p      = 30.973762_r8, & ! molecular weight P
    ! mw_Fe     = 55.845_r8,    & ! molecular weight Fe
    spp_NtoC_mass = mp(mp_spp_NtoC) * 14.00672_r8 / 12.011_r8 ! mw_n/mw_c
    spp_NtoP_mass = mp(mp_spp_NtoP) * 14.00672_r8 / 30.973762_r8 ! mw_n/mw_p
    spp_NtoFe_mass = mp(mp_spp_NtoFe) * 14.00672_r8 / 55.845_r8 ! mw_n/mw_Fe

    mp(mp_spp_labile_ratio) = 0.7_r8  ! BEC
    mp(mp_spp_spare) = 999.0_r8
    mp(mp_spp_E_Growth_ratio) = 0.2_r8
    mp(mp_spp_E_N_fraction) = 0.01_r8/0.2_r8 ! fractional mass - 5% N? Probally not
    mp(mp_spp_E_P_fraction) = mp(mp_spp_E_N_fraction) / spp_NtoP_mass ! fractional mass
    mp(mp_spp_E_Fe_fraction) = mp(mp_spp_E_N_fraction) / spp_NtoFe_mass ! fractional mass
    ! Below could use re-thinking. Logic now is that it 1/2 leftover is C, the other 1/2 Oxygen
    ! Remember it's in mass!!
    mp(mp_spp_E_C_fraction) = (1.0_r8 - mp(mp_spp_E_N_fraction) - mp(mp_spp_E_P_fraction) - mp(mp_spp_E_Fe_fraction))*0.5_r8 ! fractional mass

    mp(mp_spp_POM_from_breakage) = 0.95  ! fractional POC (vs DOC)
    mp(mp_spp_POM_from_grazing) = 0.80
    mp(mp_spp_POM_from_death) = 0.80

    mp(mp_spp_Vmax) = 12.8_r8 * 24.0_r8 * mp(mp_spp_dry_sa)
    mp(mp_spp_Ks_NO3) = 10130.0_r8
    mp(mp_spp_Vmax_NH4) = mp(mp_spp_Vmax)
    mp(mp_spp_Ks_NH4) = mp(mp_spp_Ks_NO3)
    mp(mp_spp_Qmin) = 10.18_r8
    mp(mp_spp_Qmax) = 54.0_r8

    mp(mp_spp_Vmax_PO4) = mp(mp_spp_Vmax) / mp(mp_spp_NtoP)
    mp(mp_spp_Ks_PO4) = mp(mp_spp_Ks_NO3) / mp(mp_spp_NtoP)
    mp(mp_spp_Qmin_P) = mp(mp_spp_Qmin) / spp_NtoP_mass
    mp(mp_spp_Qmax_P) =  mp(mp_spp_Qmax) / spp_NtoP_mass

    mp(mp_spp_Vmax_Fe) = mp(mp_spp_Vmax) / mp(mp_spp_NtoFe)
    mp(mp_spp_Ks_Fe) = mp(mp_spp_Ks_NO3) / mp(mp_spp_NtoFe)
    mp(mp_spp_Qmin_Fe) = mp(mp_spp_Qmin) / spp_NtoFe_mass
    mp(mp_spp_Qmax_Fe) = mp(mp_spp_Qmax) / spp_NtoFe_mass

    mp(mp_spp_initQ) = (mp(mp_spp_Qmax) + mp(mp_spp_Qmin)) / 2.0_r8
    mp(mp_spp_initQP) = (mp(mp_spp_Qmax_P) + mp(mp_spp_Qmin_P)) / 2.0_r8
    mp(mp_spp_initQFe) = (mp(mp_spp_Qmax_Fe) + mp(mp_spp_Qmin_Fe)) / 2.0_r8

    ! Controls and options
    ! ------------------------------------------------------------------

    mp(mp_use_P) = 1.0_r8
    mp(mp_use_Fe) = 1.0_r8
    mp(mp_use_C) = 1.0_r8
    mp(mp_use_NH4) = 1.0_r8

    mp(mp_seeding_type) = 0._r8  ! init_seeding_type = 0.

    mp(mp_harvest_type) = 1._r8
    mp(mp_harvest_schedule) = 1._r8
    mp(mp_harvest_avg_period) = 1.35_r8
    mp(mp_harvest_kg) = 1.35_r8
    mp(mp_harvest_f) = 0.8_r8
    mp(mp_harvest_freq) = 220._r8
    mp(mp_harvest_span) = 150._r8
    mp(mp_harvest_nmax) = 2._r8
    mp(mp_breakage_type) = 0._r8 ! breakage_Duarte_Ferreira = 0

    mp(mp_dte) = 1.0_r8/24.0_r8
    mp(mp_dt_mag) = 1.0_r8/24.0_r8
    mp(mp_N_flux_limit) = 0.0_r8
    mp(mp_farm_depth) = 2.0_r8
    mp(mp_wave_mort_factor) = 1.0_r8 ! scalar for wave/breakage mortality calculation
    mp(mp_N_uptake_type) = 1.0_r8    ! 0=gMACMODS linear, 1=non-linear
    mp(mp_growth_lim_type) =  1.0_r8 ! 0=multiplicative/interacting, 1=minimum/non-interacting
    mp(mp_Q_lim_type) = 0.0_r8       ! 0=Droop scaled 0-1 (gMACMODS), 1=Frieder et al./linear, 2=power law

    mp(mp_depth_cycle_z_photic) = 2.0_r8
    mp(mp_depth_cycle_z_deep) = 40.0_r8
    mp(mp_depth_cycle_t_photic) = 6.0_r8
    mp(mp_depth_cycle_t_deep) = 20.0_r8


  end SUBROUTINE set_macmods_params_defaults



  SUBROUTINE write_mp_parameters_block(text_block)
    !> Writes parameter name, index, and value from the 'mp' array
    !> into a long character text block.
    !>
    !> Args:
    !>   text_block (CHARACTER(LEN=*), INTENT(OUT)): The long character
    !>                                               variable to write to.
    !>
    !> Notes:
    !>   - Assumes the parameter enumeration (mp_spp_...) and the
    !>     'mp' array are accessible in the calling scope.
    !>   - 'wp' should be defined for the desired working precision of 'mp'.
    !>   - The length of 'text_block' should be sufficient to hold the output.

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(OUT) :: text_block
    CHARACTER(LEN=200) :: line

    ! Initialize the text block to empty
    text_block = ""



    ! Write each parameter information to the text block
    WRITE(line, '(ES20.8, A, I4, A, A)') mp(mp_spp_kcap), ',', mp_spp_kcap, ',', 'mp_spp_kcap'
    text_block = TRIM(text_block) // TRIM(line) // NEW_LINE(' ')

    WRITE(line, '(ES20.8, A, I4, A, A)') mp(mp_spp_Gmax_cap), ',', mp_spp_Gmax_cap, ',', 'mp_spp_Gmax_cap'
    text_block = TRIM(text_block) // TRIM(line) // NEW_LINE(' ')

    WRITE(line, '(ES20.8, A, I4, A, A)') mp(mp_spp_PARs), ',', mp_spp_PARs, ',', 'mp_spp_PARs'
    text_block = TRIM(text_block) // TRIM(line) // NEW_LINE(' ')

    WRITE(line, '(ES20.8, A, I4, A, A)') mp(mp_spp_PARc), ',', mp_spp_PARc, ',', 'mp_spp_PARc'
    text_block = TRIM(text_block) // TRIM(line) // NEW_LINE(' ')

    WRITE(line, '(ES20.8, A, I4, A, A)') mp(mp_spp_BtoSA), ',', mp_spp_BtoSA, ',', 'mp_spp_BtoSA'
    text_block = TRIM(text_block) // TRIM(line) // NEW_LINE(' ')

    WRITE(line, '(ES20.8, A, I4, A, A)') mp(mp_spp_line_sep), ',', mp_spp_line_sep, ',', 'mp_spp_line_sep'
    text_block = TRIM(text_block) // TRIM(line) // NEW_LINE(' ')

    WRITE(line, '(ES20.8, A, I4, A, A)') mp(mp_spp_kcap_rate), ',', mp_spp_kcap_rate, ',', 'mp_spp_kcap_rate'
    text_block = TRIM(text_block) // TRIM(line) // NEW_LINE(' ')

    WRITE(line, '(ES20.8, A, I4, A, A)') mp(mp_spp_Topt1), ',', mp_spp_Topt1, ',', 'mp_spp_Topt1'
    text_block = TRIM(text_block) // TRIM(line) // NEW_LINE(' ')

    WRITE(line, '(ES20.8, A, I4, A, A)') mp(mp_spp_K1), ',', mp_spp_K1, ',', 'mp_spp_K1'
    text_block = TRIM(text_block) // TRIM(line) // NEW_LINE(' ')

    WRITE(line, '(ES20.8, A, I4, A, A)') mp(mp_spp_Topt2), ',', mp_spp_Topt2, ',', 'mp_spp_Topt2'
    text_block = TRIM(text_block) // TRIM(line) // NEW_LINE(' ')

    WRITE(line, '(ES20.8, A, I4, A, A)') mp(mp_spp_K2), ',', mp_spp_K2, ',', 'mp_spp_K2'
    text_block = TRIM(text_block) // TRIM(line) // NEW_LINE(' ')

    WRITE(line, '(ES20.8, A, I4, A, A)') mp(mp_spp_CD), ',', mp_spp_CD, ',', 'mp_spp_CD'
    text_block = TRIM(text_block) // TRIM(line) // NEW_LINE(' ')

    WRITE(line, '(ES20.8, A, I4, A, A)') mp(mp_spp_dry_sa), ',', mp_spp_dry_sa, ',', 'mp_spp_dry_sa'
    text_block = TRIM(text_block) // TRIM(line) // NEW_LINE(' ')

    WRITE(line, '(ES20.8, A, I4, A, A)') mp(mp_spp_dry_wet), ',', mp_spp_dry_wet, ',', 'mp_spp_dry_wet'
    text_block = TRIM(text_block) // TRIM(line) // NEW_LINE(' ')

    WRITE(line, '(ES20.8, A, I4, A, A)') mp(mp_spp_E), ',', mp_spp_E, ',', 'mp_spp_E'
    text_block = TRIM(text_block) // TRIM(line) // NEW_LINE(' ')

    WRITE(line, '(ES20.8, A, I4, A, A)') mp(mp_spp_seed), ',', mp_spp_seed, ',', 'mp_spp_seed'
    text_block = TRIM(text_block) // TRIM(line) // NEW_LINE(' ')

    WRITE(line, '(ES20.8, A, I4, A, A)') mp(mp_spp_death), ',', mp_spp_death, ',', 'mp_spp_death'
    text_block = TRIM(text_block) // TRIM(line) // NEW_LINE(' ')

    WRITE(line, '(ES20.8, A, I4, A, A)') mp(mp_spp_BtoC), ',', mp_spp_BtoC, ',', 'mp_spp_BtoC'
    text_block = TRIM(text_block) // TRIM(line) // NEW_LINE(' ')

    WRITE(line, '(ES20.8, A, I4, A, A)') mp(mp_spp_NtoC), ',', mp_spp_NtoC, ',', 'mp_spp_NtoC'
    text_block = TRIM(text_block) // TRIM(line) // NEW_LINE(' ')

    WRITE(line, '(ES20.8, A, I4, A, A)') mp(mp_spp_NtoP), ',', mp_spp_NtoP, ',', 'mp_spp_NtoP'
    text_block = TRIM(text_block) // TRIM(line) // NEW_LINE(' ')

    WRITE(line, '(ES20.8, A, I4, A, A)') mp(mp_spp_NtoFe), ',', mp_spp_NtoFe, ',', 'mp_spp_NtoFe'
    text_block = TRIM(text_block) // TRIM(line) // NEW_LINE(' ')

    WRITE(line, '(ES20.8, A, I4, A, A)') mp(mp_spp_CtoO), ',', mp_spp_CtoO, ',', 'mp_spp_CtoO'
    text_block = TRIM(text_block) // TRIM(line) // NEW_LINE(' ')

    WRITE(line, '(ES20.8, A, I4, A, A)') mp(mp_spp_labile_ratio), ',', mp_spp_labile_ratio, ',', 'mp_spp_labile_ratio'
    text_block = TRIM(text_block) // TRIM(line) // NEW_LINE(' ')

    WRITE(line, '(ES20.8, A, I4, A, A)') mp(mp_spp_spare), ',', mp_spp_spare, ',', 'mp_spp_spare'
    text_block = TRIM(text_block) // TRIM(line) // NEW_LINE(' ')

    WRITE(line, '(ES20.8, A, I4, A, A)') mp(mp_spp_E_Growth_ratio), ',', mp_spp_E_Growth_ratio, ',', 'mp_spp_E_Growth_ratio'
    text_block = TRIM(text_block) // TRIM(line) // NEW_LINE(' ')

    WRITE(line, '(ES20.8, A, I4, A, A)') mp(mp_spp_E_C_fraction), ',', mp_spp_E_C_fraction, ',', 'mp_spp_E_C_fraction'
    text_block = TRIM(text_block) // TRIM(line) // NEW_LINE(' ')

    WRITE(line, '(ES20.8, A, I4, A, A)') mp(mp_spp_E_N_fraction), ',', mp_spp_E_N_fraction, ',', 'mp_spp_E_N_fraction'
    text_block = TRIM(text_block) // TRIM(line) // NEW_LINE(' ')

    WRITE(line, '(ES20.8, A, I4, A, A)') mp(mp_spp_E_P_fraction), ',', mp_spp_E_P_fraction, ',', 'mp_spp_E_P_fraction'
    text_block = TRIM(text_block) // TRIM(line) // NEW_LINE(' ')

    WRITE(line, '(ES20.8, A, I4, A, A)') mp(mp_spp_E_Fe_fraction), ',', mp_spp_E_Fe_fraction, ',', 'mp_spp_E_Fe_fraction'
    text_block = TRIM(text_block) // TRIM(line) // NEW_LINE(' ')

    WRITE(line, '(ES20.8, A, I4, A, A)') mp(mp_spp_POM_from_breakage), ',', mp_spp_POM_from_breakage, ',', 'mp_spp_POM_from_breakage'
    text_block = TRIM(text_block) // TRIM(line) // NEW_LINE(' ')

    WRITE(line, '(ES20.8, A, I4, A, A)') mp(mp_spp_POM_from_grazing), ',', mp_spp_POM_from_grazing, ',', 'mp_spp_POM_from_grazing'
    text_block = TRIM(text_block) // TRIM(line) // NEW_LINE(' ')

    WRITE(line, '(ES20.8, A, I4, A, A)') mp(mp_spp_POM_from_death), ',', mp_spp_POM_from_death, ',', 'mp_spp_POM_from_death'
    text_block = TRIM(text_block) // TRIM(line) // NEW_LINE(' ')

    WRITE(line, '(ES20.8, A, I4, A, A)') mp(mp_spp_Vmax), ',', mp_spp_Vmax, ',', 'mp_spp_Vmax'
    text_block = TRIM(text_block) // TRIM(line) // NEW_LINE(' ')

    WRITE(line, '(ES20.8, A, I4, A, A)') mp(mp_spp_Ks_NO3), ',', mp_spp_Ks_NO3, ',', 'mp_spp_Ks_NO3'
    text_block = TRIM(text_block) // TRIM(line) // NEW_LINE(' ')

    WRITE(line, '(ES20.8, A, I4, A, A)') mp(mp_spp_Vmax_NH4), ',', mp_spp_Vmax_NH4, ',', 'mp_spp_Vmax_NH4'
    text_block = TRIM(text_block) // TRIM(line) // NEW_LINE(' ')

    WRITE(line, '(ES20.8, A, I4, A, A)') mp(mp_spp_Ks_NH4), ',', mp_spp_Ks_NH4, ',', 'mp_spp_Ks_NH4'
    text_block = TRIM(text_block) // TRIM(line) // NEW_LINE(' ')

    WRITE(line, '(ES20.8, A, I4, A, A)') mp(mp_spp_Q0), ',', mp_spp_Q0, ',', 'mp_spp_Q0'
    text_block = TRIM(text_block) // TRIM(line) // NEW_LINE(' ')

    WRITE(line, '(ES20.8, A, I4, A, A)') mp(mp_spp_Qmax), ',', mp_spp_Qmax, ',', 'mp_spp_Qmax'
    text_block = TRIM(text_block) // TRIM(line) // NEW_LINE(' ')

    WRITE(line, '(ES20.8, A, I4, A, A)') mp(mp_spp_Qmin), ',', mp_spp_Qmin, ',', 'mp_spp_Qmin'
    text_block = TRIM(text_block) // TRIM(line) // NEW_LINE(' ')

    WRITE(line, '(ES20.8, A, I4, A, A)') mp(mp_spp_initQ), ',', mp_spp_initQ, ',', 'mp_spp_initQ'
    text_block = TRIM(text_block) // TRIM(line) // NEW_LINE(' ')

    WRITE(line, '(ES20.8, A, I4, A, A)') mp(mp_spp_Vmax_PO4), ',', mp_spp_Vmax_PO4, ',', 'mp_spp_Vmax_PO4'
    text_block = TRIM(text_block) // TRIM(line) // NEW_LINE(' ')

    WRITE(line, '(ES20.8, A, I4, A, A)') mp(mp_spp_Ks_PO4), ',', mp_spp_Ks_PO4, ',', 'mp_spp_Ks_PO4'
    text_block = TRIM(text_block) // TRIM(line) // NEW_LINE(' ')

    WRITE(line, '(ES20.8, A, I4, A, A)') mp(mp_spp_Qmin_P), ',', mp_spp_Qmin_P, ',', 'mp_spp_Qmin_P'
    text_block = TRIM(text_block) // TRIM(line) // NEW_LINE(' ')

    WRITE(line, '(ES20.8, A, I4, A, A)') mp(mp_spp_Qmax_P), ',', mp_spp_Qmax_P, ',', 'mp_spp_Qmax_P'
    text_block = TRIM(text_block) // TRIM(line) // NEW_LINE(' ')

    WRITE(line, '(ES20.8, A, I4, A, A)') mp(mp_spp_initQP), ',', mp_spp_initQP, ',', 'mp_spp_initQP'
    text_block = TRIM(text_block) // TRIM(line) // NEW_LINE(' ')

    WRITE(line, '(ES20.8, A, I4, A, A)') mp(mp_spp_Vmax_Fe), ',', mp_spp_Vmax_Fe, ',', 'mp_spp_Vmax_Fe'
    text_block = TRIM(text_block) // TRIM(line) // NEW_LINE(' ')

    WRITE(line, '(ES20.8, A, I4, A, A)') mp(mp_spp_Ks_Fe), ',', mp_spp_Ks_Fe, ',', 'mp_spp_Ks_Fe'
    text_block = TRIM(text_block) // TRIM(line) // NEW_LINE(' ')

    WRITE(line, '(ES20.8, A, I4, A, A)') mp(mp_spp_Qmin_Fe), ',', mp_spp_Qmin_Fe, ',', 'mp_spp_Qmin_Fe'
    text_block = TRIM(text_block) // TRIM(line) // NEW_LINE(' ')

    WRITE(line, '(ES20.8, A, I4, A, A)') mp(mp_spp_Qmax_Fe), ',', mp_spp_Qmax_Fe, ',', 'mp_spp_Qmax_Fe'
    text_block = TRIM(text_block) // TRIM(line) // NEW_LINE(' ')

    WRITE(line, '(ES20.8, A, I4, A, A)') mp(mp_spp_initQFe), ',', mp_spp_initQFe, ',', 'mp_spp_initQFe'
    text_block = TRIM(text_block) // TRIM(line) // NEW_LINE(' ')

    WRITE(line, '(ES20.8, A, I4, A, A)') mp(mp_use_P), ',', mp_use_P, ',', 'mp_use_P'
    text_block = TRIM(text_block) // TRIM(line) // NEW_LINE(' ')

    WRITE(line, '(ES20.8, A, I4, A, A)') mp(mp_use_Fe), ',', mp_use_Fe, ',', 'mp_use_Fe'
    text_block = TRIM(text_block) // TRIM(line) // NEW_LINE(' ')

    WRITE(line, '(ES20.8, A, I4, A, A)') mp(mp_use_C), ',', mp_use_C, ',', 'mp_use_C'
    text_block = TRIM(text_block) // TRIM(line) // NEW_LINE(' ')

    WRITE(line, '(ES20.8, A, I4, A, A)') mp(mp_use_NH4), ',', mp_use_NH4, ',', 'mp_use_NH4'
    text_block = TRIM(text_block) // TRIM(line) // NEW_LINE(' ')

    WRITE(line, '(ES20.8, A, I4, A, A)') mp(mp_seeding_type), ',', mp_seeding_type, ',', 'mp_seeding_type'
    text_block = TRIM(text_block) // TRIM(line) // NEW_LINE(' ')

    WRITE(line, '(ES20.8, A, I4, A, A)') mp(mp_harvest_type), ',', mp_harvest_type, ',', 'mp_harvest_type'
    text_block = TRIM(text_block) // TRIM(line) // NEW_LINE(' ')

    WRITE(line, '(ES20.8, A, I4, A, A)') mp(mp_harvest_schedule), ',', mp_harvest_schedule, ',', 'mp_harvest_schedule'
    text_block = TRIM(text_block) // TRIM(line) // NEW_LINE(' ')

    WRITE(line, '(ES20.8, A, I4, A, A)') mp(mp_harvest_avg_period), ',', mp_harvest_avg_period, ',', 'mp_harvest_avg_period'
    text_block = TRIM(text_block) // TRIM(line) // NEW_LINE(' ')

    WRITE(line, '(ES20.8, A, I4, A, A)') mp(mp_harvest_kg), ',', mp_harvest_kg, ',', 'mp_harvest_kg'
    text_block = TRIM(text_block) // TRIM(line) // NEW_LINE(' ')

    WRITE(line, '(ES20.8, A, I4, A, A)') mp(mp_harvest_f), ',', mp_harvest_f, ',', 'mp_harvest_f'
    text_block = TRIM(text_block) // TRIM(line) // NEW_LINE(' ')

    WRITE(line, '(ES20.8, A, I4, A, A)') mp(mp_harvest_freq), ',', mp_harvest_freq, ',', 'mp_harvest_freq'
    text_block = TRIM(text_block) // TRIM(line) // NEW_LINE(' ')

    WRITE(line, '(ES20.8, A, I4, A, A)') mp(mp_harvest_span), ',', mp_harvest_span, ',', 'mp_harvest_span'
    text_block = TRIM(text_block) // TRIM(line) // NEW_LINE(' ')

    WRITE(line, '(ES20.8, A, I4, A, A)') mp(mp_harvest_nmax), ',', mp_harvest_nmax, ',', 'mp_harvest_nmax'
    text_block = TRIM(text_block) // TRIM(line) // NEW_LINE(' ')

    WRITE(line, '(ES20.8, A, I4, A, A)') mp(mp_breakage_type), ',', mp_breakage_type, ',', 'mp_breakage_type'
    text_block = TRIM(text_block) // TRIM(line) // NEW_LINE(' ')

    WRITE(line, '(ES20.8, A, I4, A, A)') mp(mp_dte), ',', mp_dte, ',', 'mp_dte'
    text_block = TRIM(text_block) // TRIM(line) // NEW_LINE(' ')

    WRITE(line, '(ES20.8, A, I4, A, A)') mp(mp_dt_mag), ',', mp_dt_mag, ',', 'mp_dt_mag'
    text_block = TRIM(text_block) // TRIM(line) // NEW_LINE(' ')

    WRITE(line, '(ES20.8, A, I4, A, A)') mp(mp_N_flux_limit), ',', mp_N_flux_limit, ',', 'mp_N_flux_limit'
    text_block = TRIM(text_block) // TRIM(line) // NEW_LINE(' ')

    WRITE(line, '(ES20.8, A, I4, A, A)') mp(mp_farm_depth), ',', mp_farm_depth, ',', 'mp_farm_depth'
    text_block = TRIM(text_block) // TRIM(line) // NEW_LINE(' ')

    WRITE(line, '(ES20.8, A, I4, A, A)') mp(mp_wave_mort_factor), ',', mp_wave_mort_factor, ',', 'mp_wave_mort_factor'
    text_block = TRIM(text_block) // TRIM(line) // NEW_LINE(' ')

    WRITE(line, '(ES20.8, A, I4, A, A)') mp(mp_N_uptake_type), ',', mp_N_uptake_type, ',', 'mp_N_uptake_type'
    text_block = TRIM(text_block) // TRIM(line) // NEW_LINE(' ')

    WRITE(line, '(ES20.8, A, I4, A, A)') mp(mp_growth_lim_type), ',', mp_growth_lim_type, ',', 'mp_growth_lim_type'
    text_block = TRIM(text_block) // TRIM(line) // NEW_LINE(' ')

    WRITE(line, '(ES20.8, A, I4, A, A)') mp(mp_Q_lim_type), ',', mp_Q_lim_type, ',', 'mp_Q_lim_type'
    text_block = TRIM(text_block) // TRIM(line) // NEW_LINE(' ')

    WRITE(line, '(ES20.8, A, I4, A, A)') mp(mp_depth_cycle_z_photic), ',', mp_depth_cycle_z_photic, ',', 'mp_depth_cycle_z_photic'
    text_block = TRIM(text_block) // TRIM(line) // NEW_LINE(' ')

    WRITE(line, '(ES20.8, A, I4, A, A)') mp(mp_depth_cycle_z_deep), ',', mp_depth_cycle_z_deep, ',', 'mp_depth_cycle_z_deep'
    text_block = TRIM(text_block) // TRIM(line) // NEW_LINE(' ')

    WRITE(line, '(ES20.8, A, I4, A, A)') mp(mp_depth_cycle_t_photic), ',', mp_depth_cycle_t_photic, ',', 'mp_depth_cycle_t_photic'
    text_block = TRIM(text_block) // TRIM(line) // NEW_LINE(' ')

    WRITE(line, '(ES20.8, A, I4, A, A)') mp(mp_depth_cycle_t_deep), ',', mp_depth_cycle_t_deep, ',', 'mp_depth_cycle_t_deep'
    text_block = TRIM(text_block) // TRIM(line) // NEW_LINE(' ')


  END SUBROUTINE write_mp_parameters_block


  FUNCTION get_mp_idx(mp_name) RESULT(mp_idx)
    !> Returns the enumeration number for a given parameter name string.
    !>
    !> Args:
    !>     mp_name (CHARACTER(LEN=*), INTENT(IN)): The string
    !>                                                    representing the parameter name
    !>                                                    (e.g., 'mp_spp_kcap').
    !>
    !> Returns:
    !>     mp_idx (INTEGER(KIND=I4)): The corresponding enumeration
    !>                                          number, or 0 if not found.
    !>
    !> Notes:
    !>     - Assumes the parameter enumeration (mp_spp_...) is accessible
    !>       in the calling scope.

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) :: mp_name
    INTEGER(KIND=i4) :: mp_idx

    SELECT CASE (TRIM(mp_name))
        CASE ('mp_spp_kcap')
            mp_idx = mp_spp_kcap
        CASE ('mp_spp_Gmax_cap')
            mp_idx = mp_spp_Gmax_cap
        CASE ('mp_spp_PARs')
            mp_idx = mp_spp_PARs
        CASE ('mp_spp_PARc')
            mp_idx = mp_spp_PARc
        CASE ('mp_spp_BtoSA')
            mp_idx = mp_spp_BtoSA
        CASE ('mp_spp_line_sep')
            mp_idx = mp_spp_line_sep
        CASE ('mp_spp_kcap_rate')
            mp_idx = mp_spp_kcap_rate
        CASE ('mp_spp_Topt1')
            mp_idx = mp_spp_Topt1
        CASE ('mp_spp_K1')
            mp_idx = mp_spp_K1
        CASE ('mp_spp_Topt2')
            mp_idx = mp_spp_Topt2
        CASE ('mp_spp_K2')
            mp_idx = mp_spp_K2
        CASE ('mp_spp_CD')
            mp_idx = mp_spp_CD
        CASE ('mp_spp_dry_sa')
            mp_idx = mp_spp_dry_sa
        CASE ('mp_spp_dry_wet')
            mp_idx = mp_spp_dry_wet
        CASE ('mp_spp_E')
            mp_idx = mp_spp_E
        CASE ('mp_spp_seed')
            mp_idx = mp_spp_seed
        CASE ('mp_spp_death')
            mp_idx = mp_spp_death
        CASE ('mp_spp_BtoC')
            mp_idx = mp_spp_BtoC
        CASE ('mp_spp_NtoC')
            mp_idx = mp_spp_NtoC
        CASE ('mp_spp_NtoP')
            mp_idx = mp_spp_NtoP
        CASE ('mp_spp_NtoFe')
            mp_idx = mp_spp_NtoFe
        CASE ('mp_spp_CtoO')
            mp_idx = mp_spp_CtoO
        CASE ('mp_spp_labile_ratio')
            mp_idx = mp_spp_labile_ratio
        CASE ('mp_spp_spare')
            mp_idx = mp_spp_spare
        CASE ('mp_spp_E_Growth_ratio')
            mp_idx = mp_spp_E_Growth_ratio
        CASE ('mp_spp_E_C_fraction')
            mp_idx = mp_spp_E_C_fraction
        CASE ('mp_spp_E_N_fraction')
            mp_idx = mp_spp_E_N_fraction
        CASE ('mp_spp_E_P_fraction')
            mp_idx = mp_spp_E_P_fraction
        CASE ('mp_spp_E_Fe_fraction')
            mp_idx = mp_spp_E_Fe_fraction
        CASE ('mp_spp_POM_from_breakage')
            mp_idx = mp_spp_POM_from_breakage
        CASE ('mp_spp_POM_from_grazing')
            mp_idx = mp_spp_POM_from_grazing
        CASE ('mp_spp_POM_from_death')
            mp_idx = mp_spp_POM_from_death
        CASE ('mp_spp_Vmax')
            mp_idx = mp_spp_Vmax
        CASE ('mp_spp_Ks_NO3')
            mp_idx = mp_spp_Ks_NO3
        CASE ('mp_spp_Vmax_NH4')
            mp_idx = mp_spp_Vmax_NH4
        CASE ('mp_spp_Ks_NH4')
            mp_idx = mp_spp_Ks_NH4
        CASE ('mp_spp_Q0')
            mp_idx = mp_spp_Q0
        CASE ('mp_spp_Qmax')
            mp_idx = mp_spp_Qmax
        CASE ('mp_spp_Qmin')
            mp_idx = mp_spp_Qmin
        CASE ('mp_spp_initQ')
            mp_idx = mp_spp_initQ
        CASE ('mp_spp_Vmax_PO4')
            mp_idx = mp_spp_Vmax_PO4
        CASE ('mp_spp_Ks_PO4')
            mp_idx = mp_spp_Ks_PO4
        CASE ('mp_spp_Qmin_P')
            mp_idx = mp_spp_Qmin_P
        CASE ('mp_spp_Qmax_P')
            mp_idx = mp_spp_Qmax_P
        CASE ('mp_spp_initQP')
            mp_idx = mp_spp_initQP
        CASE ('mp_spp_Vmax_Fe')
            mp_idx = mp_spp_Vmax_Fe
        CASE ('mp_spp_Ks_Fe')
            mp_idx = mp_spp_Ks_Fe
        CASE ('mp_spp_Qmin_Fe')
            mp_idx = mp_spp_Qmin_Fe
        CASE ('mp_spp_Qmax_Fe')
            mp_idx = mp_spp_Qmax_Fe
        CASE ('mp_spp_initQFe')
            mp_idx = mp_spp_initQFe
        CASE ('mp_use_P')
            mp_idx = mp_use_P
        CASE ('mp_use_Fe')
            mp_idx = mp_use_Fe
        CASE ('mp_use_C')
            mp_idx = mp_use_C
        CASE ('mp_use_NH4')
            mp_idx = mp_use_NH4
        CASE ('mp_harvest_type')
            mp_idx = mp_harvest_type
        CASE ('mp_harvest_schedule')
            mp_idx = mp_harvest_schedule
        CASE ('mp_harvest_avg_period')
            mp_idx = mp_harvest_avg_period
        CASE ('mp_harvest_kg')
            mp_idx = mp_harvest_kg
        CASE ('mp_harvest_f')
            mp_idx = mp_harvest_f
        CASE ('mp_harvest_freq')
            mp_idx = mp_harvest_freq
        CASE ('mp_harvest_span')
            mp_idx = mp_harvest_span
        CASE ('mp_harvest_nmax')
            mp_idx = mp_harvest_nmax
        CASE ('mp_breakage_type')
            mp_idx = mp_breakage_type
        CASE ('mp_dte')
            mp_idx = mp_dte
        CASE ('mp_dt_mag')
            mp_idx = mp_dt_mag
        CASE ('mp_N_flux_limit')
            mp_idx = mp_N_flux_limit
        CASE ('mp_farm_depth')
            mp_idx = mp_farm_depth
        CASE ('mp_wave_mort_factor')
            mp_idx = mp_wave_mort_factor
        CASE ('mp_N_uptake_type')
            mp_idx = mp_N_uptake_type
        CASE ('mp_growth_lim_type')
          mp_idx = mp_growth_lim_type
        CASE ('mp_Q_lim_type')
          mp_idx = mp_Q_lim_type
        CASE ('mp_depth_cycle_z_photic')
          mp_idx = mp_depth_cycle_z_photic
        CASE ('mp_depth_cycle_z_deep')
          mp_idx = mp_depth_cycle_z_deep
        CASE ('mp_depth_cycle_t_photic')
          mp_idx = mp_depth_cycle_t_photic
        CASE ('mp_depth_cycle_t_deep')
          mp_idx = mp_depth_cycle_t_deep
        CASE DEFAULT
          print *,'mp param name not known: '//mp_name
          call exit(-1)
    END SELECT
  END FUNCTION get_mp_idx


  FUNCTION get_mo_idx(mo_name) RESULT(mo_idx)
    !> Returns the enumeration number for a given parameter name string.
    !>
    !> Args:
    !>     mp_name (CHARACTER(LEN=*), INTENT(IN)): The string
    !>                                                    representing the parameter name
    !>                                                    (e.g., 'mp_spp_kcap').
    !>
    !> Returns:
    !>     mp_idx (INTEGER(KIND=I4)): The corresponding enumeration
    !>                                          number, or 0 if not found.
    !>
    !> Notes:
    !>     - Assumes the parameter enumeration (mp_spp_...) is accessible
    !>       in the calling scope.

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) :: mo_name
    INTEGER(KIND=i4) :: mo_idx

    SELECT CASE (TRIM(mo_name))
      case ('mo_B')
        mo_idx = mo_B
      case ('mo_QN')
        mo_idx = mo_QN
      case ('mo_QP')
        mo_idx = mo_QP
      case ('mo_QFe')
        mo_idx = mo_QFe
      case ('mo_Gave')
        mo_idx = mo_Gave
      case ('mo_Dave')
        mo_idx = mo_Dave
      case ('mo_Growth')
        mo_idx = mo_Growth
      case ('mo_d_Be')
        mo_idx = mo_d_Be
      case ('mo_d_Bm')
        mo_idx = mo_d_Bm
      case ('mo_d_Bm_wave')
        mo_idx = mo_d_Bm_wave
      case ('mo_d_B')
        mo_idx = mo_d_B
      case ('mo_d_QN')
        mo_idx = mo_d_QN
      case ('mo_d_NO3')
        mo_idx = mo_d_NO3
      case ('mo_d_NH4')
        mo_idx = mo_d_NH4
      case ('mo_d_PO4')
        mo_idx = mo_d_PO4
      case ('mo_d_Fe')
        mo_idx = mo_d_Fe
      case ('mo_d_DIC')
        mo_idx = mo_d_DIC
      case ('mo_d_O2')
        mo_idx = mo_d_O2
      case ('mo_d_DOC')
        mo_idx = mo_d_DOC
      case ('mo_d_DON')
        mo_idx = mo_d_DON
      case ('mo_d_DOP')
        mo_idx = mo_d_DOP
      case ('mo_d_DOFe')
        mo_idx = mo_d_DOFe
      case ('mo_d_POC')
        mo_idx = mo_d_POC
      case ('mo_d_PON')
        mo_idx = mo_d_PON
      case ('mo_d_POP')
        mo_idx = mo_d_POP
      case ('mo_d_POFe')
        mo_idx = mo_d_POFe
      case ('mo_harv')
        mo_idx = mo_harv
      case ('mo_Grate')
        mo_idx = mo_Grate
      case ('mo_B_N')
        mo_idx = mo_B_N
      case ('mo_gQ')
        mo_idx = mo_gQ
      case ('mo_gT')
        mo_idx = mo_gT
      case ('mo_gE')
        mo_idx = mo_gE
      case ('mo_gH')
        mo_idx = mo_gH
      case ('mo_n_harv')
        mo_idx = mo_n_harv
      case ('mo_t_harv')
        mo_idx = mo_t_harv
      case ('mo_min_lim')
        mo_idx = mo_min_lim
      CASE DEFAULT
        print *,'mp param name not known: '//mo_name
        call exit(-1)
    END SELECT
  END FUNCTION get_mo_idx


end module macmods_param_mod