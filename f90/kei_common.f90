	module kei_common
  use kei_kinds, only: i4, r4, r8, log_kind

		use kei_parameters

    implicit none

		logical(kind=log_kind), parameter :: write_output = .true.

		! constants
		real(r4), save :: &
			grav,vonk,sbc,twopi,onepi,TK0,spd,dpy,CPo

		! times — declare ``time`` and ``dtsec`` on separate lines (NumPy f2py can map
		! consecutive ``real(r8)`` module scalars incorrectly for Python ``kei.kei_common``).
		REAL(r8), save :: time
		REAL(r8), save :: dtsec
		real(r4), save :: &
			startt,finalt
		integer(i4), save :: &
			ntime,nstart,nend,startyear

		! timocn
		real(r4), save :: &
			dto
		integer(i4), save :: &
			ndtocn

		! timatm
		real(r4), save :: &
			dta
		integer(i4), save :: &
			ndtatm

		! vert pgrid
		real(r4), save :: &
			DMAX, &
			zm(NZP1), &    ! midpoint layer depth (negative below srf)
			zmp(NZP1), &   ! positive distance from srf to layer midpoint
			hm(NZP1), &    ! layer thickness
			dm(0:NZ)       ! layer bounardies (positive from srf)

		! vert tgrid
		integer(i4), save :: &
			NZDIV,NZT,NZP1t,last_wct_ii

		! location
		real(r4), save :: &
			dlat,dlon,rlat,rlon,f

		! proc switches
		logical(kind=log_kind), save :: &
			LKPP  , LRI, LDD, LICE, LECO, LNBFLX, LTGRID, &
      LSW, & ! seaweed switch
    	lsaveaverages, lstretchgrid, lrepeatdat, &
    	lradtq, lfluxSSTdat, &
    	lLWupSSTdat, lrhdat, lclddat, &
    	LKPROF     ! added 7/2011 by ben saenz, to enable compilation, does not appear to be used
    logical(kind=log_kind), save :: dlLWupSSTdat

    ! proc pars
    integer(i4), save :: &
    	jerlov

    ! ecosys additions
    REAL(r8), save :: &
      absorp_baseline           ! shortwave absorption to be added to ecosys absorption

    ! ocn inital cond
    integer(i4), save :: &
    		mu,mx
    character (len=40) :: ncf_file
    integer(i4), save :: ncf_id
  	real(r4), save :: &
			hu,SSU(NVEL),hx,SSX(NSCLR)

		! ocn state
		real(r4), save :: &
			focn,Tref,SSref,uref,vref,Qs,ug0,vg0,albocn

		! ocn paras
		real(r4), save :: &
			CP(0:NZP1tmax),rho(0:NZP1tmax), &
    	rhofrsh,rhoh2o,rhob,rhosw,Sref,epsw, &
    	talpha(0:NZP1tmax),sbeta(0:NZP1tmax)
    real(r4), save :: &
      sal_correction_rate, & ! (kg/s)
      sw_scale_factor

    ! ocn energy
    real(r4), save :: &
    	eflx,esnk,Tmke,Ptke,rmke(NZP1)
    ! ``qsw_absorbed`` lives in module ``kei_ocn`` (merged from former ``kei_ocncommon``).
    ! atm advec
    integer(i4), save :: &
    	modatm
    real(r4), save :: &
      atmad

    ! atm state
    integer(i4), save :: &
    	LAFLX
    real(r4), save :: &
      u10a

    ! atm paras
    real(r4), save :: &
      az,azeta,psima,psisa,ustara,bp,cpa,rhoa

    ! lnd state
    real(r4), save :: &
      flnd,alblnd,alblndlnd,alblndsnw

    ! flx calc
    integer(i4), save :: &
      ndtflx,ndtld
    real(r4), save :: &
      deltax(NJDT),xbar,denom

		! flx sfc
    real(r4), save :: &
      sflux(NSFLXS,5,0:NJDT),VAF(NFDATAP1)

		! flx profs
    real(r4), save :: &
      wU(0:NZtmax,NVP1),wX(0:NZtmax,NSP1), &
    	wXNT(0:NZtmax,NSCLR), &
    	wtI(0:NZ,0:NJDT),wsI(0:NZ,0:NJDT)

		! x coeff
    real(r4), save :: &
      cdw,ctw,cqw,cdi,cti,cqi

		! kprof in
    real(r4), save :: &
      ustar,B0,buoy(NZP1tmax)

		! kprof out
    integer(i4), save :: &
      kmix
    real(r4), save :: &
      hmix,difm(0:NZtmax),difs(0:NZtmax),ghat(NZtmax),V2(NZP1tmax)

		! Ri mixing
    logical(kind=log_kind), save :: &
    	firstF
    real(r4), save :: &
      deltaRi, &
    	Ri(NZtmax), &
    	Fint(NZtmax), &
    	F_of_Ri(0:MRp1) ! changed from weird 'F of Ri' variable, not sure how this ever compiled! - saenz 7/2001

  	! dble diff
    real(r4), save :: &
      dift(0:NZtmax),Rp(NZtmax)

    ! forcing
    integer(i4), save :: &
      nharm(NSFLXS)
    real(r4), save :: &
      harm0(NSFLXS),areg(NSFLXS), &
    	ampharm(0:NDHARM,NSFLXS), &
    	freharm(NDHARM,NSFLXS), &
    	phharm(NDHARM,NSFLXS),xnoise(NSFLXS)

		! forcread
    integer(i4), save :: &
    	idataold,idatanew
    real(r4), save :: &
      fdata(nfdata,0:1),delfdata(nfdata)
 !   character (len=40) :: fdatafile

		! forcrepeat
    REAL(r8), save :: &
      repeatdatdays,offdatday,startdatday

		! forc data
    real(r4), save :: &
      TZdat,Tdewdat,clddat,SSTdat

		! in/output
		!integer, save :: &
    !  inct,incz,incx,incy
    !real, save :: &
    ! buffer(NBUFF),dout(NDOUT)
    logical(kind=log_kind), save :: &
    	LSTORE

		! run info
    integer(i4), save :: &
			ntitle
    real(r4), save :: &
      runseq
    character (len=52) :: dirname
    character (len=52) :: outdirname
    character (len=50), dimension(4) :: title

		! averaging of output arrays
    real(r4), save :: &
      Uav(NZP1,NVEL),Xav(NZP1,NSCLR), & !doutav(NDOUT), &
    	sfluxav(NSFLXS,5),rmkeav(NZ),eflxav,esnkav, &
    	Ptkeav

    integer(i4), SAVE :: &
      nmodeadv(2), &  ! number of advection types : temp(1), salinity(2)
      modeadv(1:maxmodeadv,1:2) ! modes for tracers additions/subtractions
    real(r4), SAVE :: &
      advection(1:maxmodeadv,2) ! modes for tracer advection

		real(r4), save :: &
			uZ, vZ, QSWins, QLWdwn, TZ, qZ, Prain, Psnow, DivU, SSTvaf, &
			QSWup,QLWup, msl, hum

    real(r4), save :: atm_flux_to_ocn_surface

    equivalence &
			(vaf(1),uZ), &   ! U windspeed (m/s)
			(vaf(2),vZ), &  ! V windspeed (m/s)
			(vaf(3),QSWins), &  ! short wave insolation flux (W/m^2 ??)
			(vaf(4),QLWdwn), &  ! long wave radiation flux down (W/m^2 ??)
			(vaf(5),TZ), &    ! air temp (degC)
			(vaf(6),qZ) , &  ! humidity (kg/m^3)
			(vaf(7),Prain), & ! rain precip (kg/m^2/s)
			(vaf(8),Psnow), &  ! snow precip (kg/m^2/s)
			(vaf(9),DivU) , &  ! ice divergence (units??)
			(vaf(10),SSTvaf), & ! forced sea surface temp (degC)
			(vaf(11),msl), & ! mean sea level pressure (mbar)
			(vaf(12),hum), & ! air specific humidity (kg/kg)
			(vaf(nfdata),QSWup), &
			(vaf(nfdatap1),QLWup)

    integer(i4), save :: f_wct   ! switch for water column temperature data assimilation
    real(r4), save :: wct_interp(NZP1) ! water column assimilation temperatures


	CONTAINS

    SUBROUTINE init_constants_params

      USE kei_icecommon

      IMPLICIT NONE

      ! run options
      startyear = 2000
      startt    = 0.0    ! (days) - from 0:00z 1 Jan
      nstart    = 0      ! pretty sure this will always be zero - no need to restart runs with modern computers
      nend      = 13140  ! number total time steps
      dtsec     = 3600.0_r8 ! time step in seconds
      ndtatm    = 1      ! sub-dt atm steps
      ndtice    = 1      ! sub-dt ice steps, likely not used b/c SIESTA does internal stepping
      ndtocn    = 1      ! sub-dt ocean/kpp steps

      LKPP      = .true.
      LRI       = .true.
      LDD       = .true.
      LICE      = .false.
      LECO      = .true.
      LSW       = .true.
      LTGRID    = .true.
      NZDIV     = 4
      LNBFLX    = .true.

      !DMAX      = 200.  ! dunno, hopefully get rid of me
      dlon      = -66.0  ! degrees E
      dlat      = -63.0  ! degrees N

      fice      = 0.0
      hice(0)   = 0.0
      hsn(0)    = 0.0
      TI0       = 0.0
      hfsnow    = 0.0
      !!! hopefully everything after this in rdinput is not needed

      ! numerical constants
      spd  = 24. * 3600.  ! seconds per day                    s/day
      dpy  = 365.00       ! days per year                      day/yr
      twopi= 8.* atan(1.) ! 2 * pie
      onepi= twopi/2.     ! pie
      ! physics
      grav = 9.816        ! Gravitational Constant             m/s2
      vonk = 0.4          ! Von Karman's constant
      TK0  = 273.15       ! Kelvin of 0 degree Celsius         K
      C2K  = TK0
      sbc  = 5.67E-8      ! Stefan-Boltzmann-Constant          W/m2/K4
      ! atmos
      az   = 10.0         ! height of atmospheric state variables  m
      ! ocean
      rhofrsh = 1000.      ! density of freshwater
      rhosw= 1000.26      ! reference ocean density for exact budgets
      CPo  = 4.1e6/rhosw  ! reference heat capacity for exact budgets
      EL   = 4187.*(597.31-0.56525*1.0) ! latent heat evap at 1C in j/kg
      epsw = 1.0          ! cor.fac. for departure of h2o from blackbody
      albocn = 0.07       ! albedo for seawater, only used when running
      !                           without rad/conv model:
      nsice = 1           ! bndry between deep and shallow frazil

      ! ice
      Tmlt = 0.0          ! melting  temperature of ice        C
      Tfrz = TfrzC(34.,0.)! freezing temperature of sea water  C (dynamic)
      Sfrazil = 4.0       ! salinity of deep frazil ice        o/ooo
      Sice =   4.0        ! salinity of ocean 2 ice frazil     o/ooo
      Sslush =   18.0        ! salinity of slushy sea ice     o/ooo
      rhoice = 897.       ! density of ice                     kg/m3
      CPice  = 1.883e6 / rhoice  ! ice heat capacity  j/kg/C
      epsi = 1.0          ! emmissivity of ice
      FL   =  334000.     ! specific heat of fusion for ice    J/kg
      SL   = EL + FL      ! latent of sublimation for ice + snow  J/kg
      FLSN =  FL          ! specific heat of fusion for snow   J/kg
      albice = 0.80       ! initial albedo (dynamic)

      ! snow
      rhosn     = 330.    ! density of snow                    kg/m3
      rhofsnow  = 944.    ! density of flooded snow             kg/m3
      ! rhofsnow  = 600.    ! density of frozen snow             kg/m3
      CPsn      = CPice
      CPfsnow   = CPice
      ! land
      alblndlnd = 0.1     ! albedo for bare landsurface
      alblndsnw = 0.85    ! albedo for snowcovered land surface

      ! other stuff

      lradtq = .false. ! use Tair and q(humidity) from rad/conv model
      lfluxSSTdat = .false. ! use "mobs" data SST to compute fluxes
      jerlov = 3       ! jerlov water type for SW absorption:
      !                  0(seasonal variation, see "fluxes.f")
      !                  1=I, 2=IA, 3=IB, 4=II, 5=III
      LAFLX = 6  ! use data - this should go away at some point, as atm is simplified


      ! ocean advection (through modification of right hand side: "rhsmod")
      nmodeadv(1) = 0    ! number of advection types : temp
      nmodeadv(2) = 0    ! ..     .. ..        ..    : salinity
      !                                  (0 = no advection/rhs modification,
      !                                  maximum = maxmodeadv = 6)
      modeadv(1:maxmodeadv,1) = (/4,7,0,0,0,0/)   ! set modes for temperature
      modeadv(1:maxmodeadv,2) = (/4,0,0,0,0,0/) ! set modes for salinity/tracers
      advection(1:maxmodeadv,1) = (/-2.00, -23.0 ,0.,0.,0.,0./) ! Temp advection (W/m2)
      advection(1:maxmodeadv,2) = (/1.00E-6, 0.,0.,0.,0.,0./)  ! Salt advection (kg/m2s)
      !     mode (km = kmixe = index of gridpoint just below hmix,
      !           dm = depth at d(kmixe)
      !     1 : Steady upper layer horizontal advection
      !     2 : Steady mixed layer horizontal advection to km-1
      !     3 : Steady horizontal advection throughout the entire column
      !     4 : Steady vert. advection (=deep horizontal) 120m to bottom
      !     5 : Steady bottom diffusion (lowest layer)
      !     6 : Seasonal boundary layer horizontal advection to dm
      !     7 : Seasonal thermocline horizontal advection km to 100m

      ! numerics
      lsaveaverages = .true. ! saving averaged variables,fluxes,dout

      ! grid
      lstretchgrid = .false. ! stretched grid, or (false) constant

      ! forcing data
      lrepeatdat = .false. ! repeat forcing data
      repeatdatdays = 365.0  !        after "repeatdatdays" (days)

      ! rad/conv model
      dlLWupSSTdat = .false. ! use data SST for LWup for rad/conv model
      lrhdat = .true.   ! using data (not .77) for q in rad/conv
      lclddat = .true.  ! using data for cloud fraction
      modatm = 3        ! mode for atm advection (in"updatm")
      !                                  0: no, 1: lower 18&17&16,2:lower 1/2,
      !                                  3: entire col., 4: combination
      atmad = 50.0     ! advection to atm (W/m2)

      sal_correction_rate = 0.
      sw_scale_factor = 1.0  ! 1 = use full shortwave

      f_wct = 0 ! default no data assimilation
      wct_interp = -999.

      ! ecosys extras
      absorp_baseline = 0.

      RETURN
    END SUBROUTINE init_constants_params

!**********************************************************************
    ! Literal kind(4) so f2py-generated wrappers do not reference r4 from kei_kinds.
    real(4) FUNCTION TfrzC(S,Db)

			implicit none

			real(4), intent(in) :: S, Db
!     Freezing point of water in degrees C at salinity S in PSU
!                                         and pressure Db in decibars
    !TfrzC = (-0.0575 +1.710523e-3 *sqrt(S) -2.154996e-4 *S) *S &
    !- 7.53e-4 *Db

		TfrzC = -0.054*S - 7.53e-4 *Db


    return
    end FUNCTION TfrzC

! ******************************************************************

    real(4) FUNCTION CPSW(S,T1,P0)

! ******************************************************************
! UNITS:
!       PRESSURE        P0       DECIBARS
!       TEMPERATURE     T        DEG CELSIUS (IPTS-68)
!       SALINITY        S        (IPSS-78)
!       SPECIFIC HEAT   CPSW     J/(KG DEG C)
! ***
! REF: MILLERO ET AL,1973,JGR,78,4499-4507
!       MILLERO ET AL, UNESCO REPORT NO. 38 1981 PP. 99-188.
! PRESSURE VARIATION FROM LEAST SQUARES POLYNOMIAL
! DEVELOPED BY FOFONOFF 1980.
! ***
! CHECK VALUE: CPSW = 3849.500 J/(KG DEG. C) FOR S = 40 (IPSS-78),
! T = 40 DEG C, P0= 10000 DECIBARS

		! inputs
    real(4) :: S,T1,P0

    ! local
    real(4) :: T,P,SR,A,B,C,CP0,CP1,CP2

!   check that temperature is above -2
    T = T1
    if(T < -2.) T = -2.

!   SCALE PRESSURE TO BARS
    P=P0/10.
! ***
! SQRT SALINITY FOR FRACTIONAL TERMS
    SR = SQRT(ABS(S))
! SPECIFIC HEAT CP0 FOR P=0 (MILLERO ET AL ,UNESCO 1981)
    A = (-1.38385E-3*T+0.1072763)*T-7.643575
    B = (5.148E-5*T-4.07718E-3)*T+0.1770383
    C = (((2.093236E-5*T-2.654387E-3)*T+0.1412855)*T &
    -3.720283)*T+4217.4
    CP0 = (B*SR + A)*S + C
! CP1 PRESSURE AND TEMPERATURE TERMS FOR S = 0
    A = (((1.7168E-8*T+2.0357E-6)*T-3.13885E-4)*T+1.45747E-2)*T &
    -0.49592
    B = (((2.2956E-11*T-4.0027E-9)*T+2.87533E-7)*T-1.08645E-5)*T &
    +2.4931E-4
    C = ((6.136E-13*T-6.5637E-11)*T+2.6380E-9)*T-5.422E-8
    CP1 = ((C*P+B)*P+A)*P
! CP2 PRESSURE AND TEMPERATURE TERMS FOR S > 0
    A = (((-2.9179E-10*T+2.5941E-8)*T+9.802E-7)*T-1.28315E-4)*T &
    +4.9247E-3
    B = (3.122E-8*T-1.517E-6)*T-1.2331E-4
    A = (A+B*SR)*S
    B = ((1.8448E-11*T-2.3905E-9)*T+1.17054E-7)*T-2.9558E-6
    B = (B+9.971E-8*SR)*S
    C = (3.513E-13*T-1.7682E-11)*T+5.540E-10
    C = (C-1.4300E-12*T*SR)*S
    CP2 = ((C*P+B)*P+A)*P
! SPECIFIC HEAT RETURN
    CPSW = CP0 + CP1 + CP2
    RETURN
    END FUNCTION CPSW

!*******************************************************************
!    subroutine acopy(a,b,ndim)
!!                                 copy array a to b
!    dimension a(ndim),b(ndim)
!
!    do 5 n=1,ndim
!        b(n) = a(n)
!    5 END DO
!    return
!    end subroutine acopy
!!*****************************************************
!    subroutine aset(array,ndim,value)
!!                         set array to value
!    dimension array(ndim)
!    do 5 n=1,ndim
!        array(n) = value
!    5 END DO
!    return
!    end subroutine aset
!*****************************************************



	end module kei_common

