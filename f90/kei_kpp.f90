!********************************************************************
    subroutine vmix(U,X,hmixn,kmixn)
!  Interface between 1-d model and vertical mixing

    use kei_parameters
    use kei_common
		use kei_icecommon
		use kei_subs1D

		implicit none

    integer, parameter :: imt = NX*NY


! inputs including those from common.inc and parameter.inc
    real :: U(nzp1,nvel),X(nzp1,nsclr)
!     real zm(nzp1) ! vertical layer grid                  (m)
!     real hm(nzp1) ! layer thicknesses                    (m)
!     real sflux(NSFLXS,5,0:NJDT)  !  surface flux array
!     real f        ! local inertial frequecy
!     real time     !  ocean time in days
!     integer jerlov ! water type; =0 for seasonal cycle
!     integer nz,nzp1  ! number of vertical levels
!     logical LRI   ! Ri_iw interior mixing switch
!     logical LDD   ! double diffusion switch

! outputs including those to common.inc
!     real difm(0:kmp1)  ! vertical viscosity coefficient  (m^2/s)
!     real difs(0:kmp1)  ! vertical scalar diffusivity     (m^2/s)
!     real dift(0:kmp1)  ! vertical temperature diffusivity(m^2/s)
!     real ghat(km)     ! nonlocal transport              (s/m^2)
    real :: hmixn          ! boundary layer depth (m)
    integer :: kmixn
!     real talpha(0:nzp1)   ! alpha
!     real sbeta(0:nzp1)    ! beta
!     real wU(0,1),wU(0,2)  ! kinematic surface momintum fluxes
!     real wX(0,1),wX(0,2)  ! kinematic surface heat and salt fluxes
!     real CP(0:nzp1)       ! specific heat             (j  / kg / deg)
!     real rho(0:nzp1)      ! density                   (kg / m^3)
!     real rhoh2o, rhob     ! density of freshwater and of brine (kg/m^3)

! local
    real :: Shsq(nzp1)    ! (local velocity shear)^2       (m/s)^2
    real :: dVsq(nzp1)    ! (velocity shear re sfc)^2      (m/s)^2
    real :: dbloc(nz)     ! local delta buoyancy            (m/s^2)
    real :: Ritop(nz)     ! numerator of bulk Richardson Number (m/s)^2
!          Ritop = (-z - -zref)* delta buoyancy w/ respect to sfc(m/s^2)
    real :: alphaDT(nz)   ! alpha * DT across interfaces
    real :: betaDS(nz)    ! beta  * DS across interfaces
!     logical LKPP       ! kpp boundary layer mixing switch
!          month  1   2   3   4   5   6   7   8   9   10  11  12
    integer :: jerl(12) = &
    	(/ 2 , 2 , 2 , 3 , 3 , 3 , 4 , 4 , 4 , 4 , 3 , 2 /)

		integer :: n,k,kl,mon,jwtype(imt),kmixn1(imt)
		real :: epsilon, epsln,swtime,alpha,beta,exppr,sigma,sigma0,tau, &
			Bo,Bosol,zref,wz,bref,del,dlimit,vlimit,ustar1(imt),Bo1(imt),&
			Bosol1(imt),hmixn1(imt)

    real :: Shsq1(imt,nzp1)    ! (local velocity shear)^2       (m/s)^2
    real :: dVsq1(imt,nzp1)    ! (velocity shear re sfc)^2      (m/s)^2
    real :: dbloc1(imt,nz)     ! local delta buoyancy            (m/s^2)
    real :: Ritop1(imt,nz)     ! numerator of bulk Richardson Number (m/s)^2
!          Ritop = (-z - -zref)* delta buoyancy w/ respect to sfc(m/s^2)
    real :: alphaDT1(imt,nzp1)   ! alpha * DT across interfaces
    real :: betaDS1(imt,nzp1)    ! beta  * DS across interfaces

     real :: difm1(imt,0:nzp1)  ! vertical viscosity coefficient  (m^2/s)
     real :: difs1(imt,0:nzp1)  ! vertical scalar diffusivity     (m^2/s)
     real :: dift1(imt,0:nzp1)  ! vertical temperature diffusivity(m^2/s)
     real :: ghat1(imt,nz)     ! nonlocal transport              (s/m^2)

! ------------------------------------------------------

    epsilon = 0.1
    epsln   = 1.e-20

!                        ensure that bottom layer isn't zero thickness
    hm(nzp1) = AMAX1(hm(nzp1),epsln)
!                        find the jerlov water type for this month
    if(jerlov < 1) then
        swtime = time
        11 if(swtime >= dpy) then
            swtime = swtime-dpy
            goto 11
        endif
        mon = int(dble(swtime/dpy)*12.) + 1 !(swtime/30.147) with dpy=365.
        mon = MIN0(mon,12)
        jwtype = jerl(mon)
    else
        jwtype = jerlov
    endif

! calculate density of fresh water and brine in surface layer
    alpha = 1.
    beta  = 1.
    exppr = 0.0
    call ABK80(0.0,X(1,1),-zm(1),alpha,beta,exppr,sigma0,sigma)
    rhoh2o    = 1000. + sigma0
    call ABK80(Sice,X(1,1),-zm(1),alpha,beta,exppr,sigma0,sigma)
    rhob      = 1000. + sigma0

!     calculate temperature and salt contributions of buoyancy gradients
!      calculate buoyancy profile (m/s**2) on gridlevels
    do 10 k=1,nzp1
        call ABK80(X(k,2)+Sref,X(k,1),-zm(k),alpha,beta,exppr, &
        sigma0,sigma)
        rho(k)    = 1000. + sigma0
        CP(k)     = CPSW(X(k,2)+Sref,X(k,1),-zm(k))
        talpha(k) = alpha
        sbeta(k)  = beta
        buoy(k)   = - grav * sigma0 / 1000.
    10 END DO

    if(LAFLX /= 0) then
        rho(0)    = rho(1)
        CP(0)     = CP(1)
        talpha(0) = talpha(1)
        sbeta(0)  = sbeta(1)
    else
    !                CONSERVE
        rho(0) = rhosw
        CP(0)  = CPo
        rhoh2o = rhofrsh
        talpha(0) = talpha(1)
        sbeta(0)  = sbeta(1)
    !             write(6,*) 'conserve', rhoh2o,rhofrsh
    endif


! calculate kinematic surface momentum fluxes
    wU(0,1) = -sflux(1,5,0) / rho(0)
    wU(0,2) = -sflux(2,5,0) / rho(0)
    tau     = sqrt( sflux(1,5,0)**2 + sflux(2,5,0)**2 )
    ustar   = sqrt( tau / rho(0) )
    ustar1 = ustar

! total turbulent kinematic temperature flux (C m/s)
    wX(0,1)  = -(sflux(4,5,0)+sflux(8,5,0)) /rho(0) /CP(0) ! as+ i2o+frazil

! total turbulent kinematic salinity flux (o/oo m/s)
!    frazil brine = added to deep  + fresh shallow + shallow salt

    wX(0,2) = Sref   *sflux(6,5,0)/rhoh2o  &  !  i2o+as  freshwater>0
    - 1000.  *sflux(5,5,0)/rho(0)  &  !  i2o salt <0
    - Sref   *sflux(8,5,0)/rhoh2o/FL & !  o2i shallow frazil frsh
    - 1000.  *sflux(9,5,0)/rho(0)    ! o2i shallow frazil salt

!       write(6,*) 'kpp',wX(0,1),wX(0,2),rhoh2o,rho(0),FL
!       call sprint(6,0)

!           deep frazil heat and salt are in NT fluxes.

! calculate total kinematic surface buoyancy flux (m**2/s**3)
    Bo = - grav*(talpha(0)*wX(0,1) - sbeta(0)*wX(0,2) )
    wX(0,NSP1) =  - Bo
    Bosol = grav * talpha(0) * sflux(3,5,0) / (rho(0)*CP(0))

!     calculate temperature and salt contributions of buoyancy gradients
!               on interfaces for double diffusion
    do 105  n = 1,nz
        alphaDT(n) =0.5 *(talpha(n)+talpha(n+1)) *(X(n,1) -X(n+1,1))
        betaDS(n)  =0.5 *(sbeta(n) + sbeta(n+1)) *(X(n,2) -X(n+1,2))
    105 END DO

!                      compute buoyancy and shear profiles
    do 115  n = 1,nz
        zref =  epsilon * zm(n)

    !          compute reference buoyancy and velocity
        wz    = AMAX1(zm(1),zref)
        uref  = U(1,1) * wz / zref
        vref  = U(1,2) * wz / zref
        bref  = buoy(1)* wz / zref
        do 125 kl = 1,nz
            IF(zref >= zm(kl)) go to 126
            wz = AMIN1(zm(kl)-zm(kl+1),zm(kl)-zref)
            del = 0.5 * wz / (zm(kl) - zm(kl+1))
            uref =uref -wz*( U(kl,1) + del *( U(kl+1,1)- U(kl,1))) /zref
            vref =vref -wz*( U(kl,2) + del *( U(kl+1,2)- U(kl,2))) /zref
            bref =bref -wz*(buoy(kl) + del *(buoy(kl+1)-buoy(kl))) /zref
        125 END DO
        126 continue

        Ritop(n) = (zref - zm(n)) * (bref - buoy(n))
        dbloc(n) = buoy(n) - buoy(n+1)
        dVsq(n)  = (Uref - U(n,1))**2   + ( Vref  - U(n,2))**2
        Shsq(n)  = (U(n,1)-U(n+1,1))**2 + (U(n,2)-U(n+1,2))**2
    115 END DO

    Bo1 = Bo
    Bosol1 = Bosol

    do n = 1,nz
    dbloc1(1,n) = dbloc(n)
    Ritop1(1,n) = Ritop(n)
    ghat1(1,n) = ghat(n)
    alphaDT1(1,n) = alphaDT(n)  ! npz1 long, but apparently last element is not used?
    betaDS1(1,n) = betaDS(n)    ! npz1 long, but apparently last element is not used?
    enddo
    do n = 1,nzp1
    Shsq1(1,n) = Shsq(n)
    dVsq1(1,n) = dVsq(n)

    difm1(1,n) = difm(n)
    difm1(1,n) = dift(n)
    difs1(1,n) = difs(n)
    enddo

    difm1(1,0) = difm(0)
    difm1(1,0) = dift(0)
    difs1(1,0) = difs(0)

    call kppmix( LRI, LDD, LKPP, &
    zm, hm , Shsq1, dVsq1, &
    ustar1 , Bo1    , Bosol1 , alphaDT1, betaDS1 , &
    dbloc1 , Ritop1 , f , jwtype, &
    difm1  , difs1  , dift1  , ghat1 , hmixn1, kmixn1)

    hmixn = hmixn1(1)
    kmixn = kmixn1(1)
    do n = 0,nzp1
	  difm(n) = difm1(1,n)
	  dift(n) = dift1(1,n)
	  difs(n) = difs1(1,n)
    enddo

!      limit the bottom diffusity and viscosity
!      zero diffusivities for no bottom flux option
    if(LNBFLX) then
        dlimit = 0.0
        vlimit = 0.0
        do n=1,nsclr
            wxNT(nz,n) = 0.0
        enddo
    else
        dlimit = 0.00001
        vlimit = 0.0001
    endif
    do k=nz,nzp1
        difm(k) = vlimit
        difs(k) = dlimit
        dift(k) = dlimit
    enddo
    ghat(nz) = 0.0

    return
    end subroutine vmix
!*********************************************************************

    subroutine kppmix( LRI, LDD, LKPP, &
    zgrid , hwide , Shsq, dVsq, &
    ustar , Bo    , Bosol , alphaDT , betaDS , &
    dbloc , Ritop , Coriol, jwtype, &
    visc  , difs  , dift  , ghats , hbl , kbl )
!.......................................................................

!     Main driver subroutine for kpp vertical mixing scheme and
!     interface to greater ocean model

!     written by : bill large,   june,  6, 1994
!     modified by: jan morzel,   june, 30, 1994
!                  bill large, august, 11, 1994
!                  bill large, november 1994,   for 1d code

!.......................................................................

    use kei_parameters

		implicit none

    integer, parameter :: km = NZ
    integer, parameter :: kmp1 = nzp1
    integer, parameter :: imt = NX*NY
    integer, parameter :: mdiff = 3  ! number of diffusivities for local arrays

! input
    real :: zgrid(imt,kmp1)   ! vertical grid (<= 0)            (m)
    real :: hwide(imt,kmp1)   ! layer thicknesses               (m)
    real :: Shsq(imt,kmp1)    ! (local velocity shear)^2       (m/s)^2
    real :: dVsq(imt,kmp1)    ! (velocity shear re sfc)^2      (m/s)^2
    real :: ustar(imt)        ! surface friction velocity       (m/s)
    real :: Bo(imt)           ! surface turbulent buoy. forcing (m^2/s^3)
    real :: Bosol(imt)        ! radiative buoyancy forcing      (m^2/s^3)
    real :: alphaDT(imt,kmp1) ! alpha * DT  across interfaces
    real :: betaDS(imt,kmp1)  ! beta  * DS  across interfaces
    real :: dbloc(imt,km)     ! local delta buoyancy            (m/s^2)
    real :: Ritop(imt,km)     ! numerator of bulk Richardson Number (m/s)^2
!          Ritop = (-z - -zref)* delta buoyancy w/ respect to sfc(m/s^2)
    real :: Coriol            ! Coriolis parameter              (s^{-1})
    integer :: jwtype(imt)    ! Jerlov water type               (1 -- 5)
    logical :: LRI, LDD, LKPP ! mixing process switches

! output
    real :: visc(imt,0:kmp1)  ! vertical viscosity coefficient  (m^2/s)
    real :: difs(imt,0:kmp1)  ! vertical scalar diffusivity     (m^2/s)
    real :: dift(imt,0:kmp1)  ! vertical temperature diffusivity(m^2/s)
    real :: ghats(imt,km)     ! nonlocal transport              (s/m^2)
    real :: hbl(imt)          ! boundary layer depth (m)
    integer :: kbl(imt)       ! index of first grid level below hbl

! local
    real :: bfsfc(imt)        ! surface buoyancy forcing        (m^2/s^3)
    real :: ws(imt)           ! momentum velocity scale
    real :: wm(imt)           ! scalar   velocity scale
    real :: caseA(imt)        ! = 1 in case A; =0 in case B
    real :: stable(imt)       ! = 1 in stable forcing; =0 in unstable
    real :: dkm1(imt,mdiff)   ! boundary layer difs at kbl-1 level
    real :: gat1(imt,mdiff)   ! shape function at sigma=1
    real :: dat1(imt,mdiff)   ! derivative of shape function at sigma=1
    real :: blmc(imt,km,mdiff)! boundary layer mixing coefficients
    real :: sigma(imt)        ! normalized depth (d / hbl)
    real :: Rib(imt,2)        ! bulk Richardson number
		integer :: ki,i


! zero the mixing coefficients
    do 310 ki=0,km
        do 300 i =1,imt
            visc(i,ki) = 0.0
            difs(i,ki) = 0.0
            dift(i,ki) = 0.0
        300 END DO
    310 END DO

! compute RI and IW interior diffusivities everywhere
    IF(LRI) THEN
        call ri_iwmix ( km  , kmp1, imt   , &
        Shsq,  dbloc , zgrid , &
        visc, difs, dift  )
    ENDIF

! add double diffusion if desired
    IF(LDD) THEN
        call ddmix    ( km  , kmp1, imt   , &
        alphaDT, betaDS , zgrid , &
        visc, difs, dift  )
    ENDIF

! fill the bottom kmp1 coefficients for blmix
    do 100 i = 1,imt
        visc (i,kmp1) = visc(i,km)
        difs (i,kmp1) = difs(i,km)
        dift (i,kmp1) = dift(i,km)
    100 END DO


!  compute boundary layer mixing coefficients ??
    IF(LKPP) THEN
    ! diagnose the new boundary layer depth
        call  bldepth (km   , kmp1 , imt   , zgrid, hwide, dVsq, &
        dbloc, Ritop, ustar , Bo   , Bosol, Coriol, jwtype, &
        hbl  , bfsfc, stable, caseA, kbl  , &
        Rib  , sigma, wm    , ws  )

    ! compute boundary layer diffusivities
        call blmix   (km   , kmp1 , imt , mdiff , zgrid, hwide , &
        ustar, bfsfc, hbl , stable, caseA, &
        visc , difs , dift, kbl   , &
        gat1 , dat1 , dkm1, blmc  , ghats, &
        sigma, wm   , ws  )

    ! enhance diffusivity at interface kbl - 1
        call enhance (km   , kmp1 , imt , mdiff , dkm1  , visc , &
        difs , dift , hbl , kbl   , zgrid , caseA, &
        blmc ,ghats )

    ! combine interior and boundary layer coefficients and nonlocal term
        do 200 ki= 1,km
            do 190 i = 1,imt
                if(ki < kbl(i)) then
                    visc(i,ki)=blmc(i,ki,1)
                    difs(i,ki)=blmc(i,ki,2)
                    dift(i,ki)=blmc(i,ki,3)
                else
                    ghats(i,ki)=0.
                endif
            190 END DO
        200 END DO

    ENDIF              ! of LKPP

    return
    end subroutine kppmix

! ********************************************************************

    subroutine  bldepth ( &
    km   , kmp1 , imt   , zgrid, hwide, dVsq , &
    dbloc, Ritop, ustar , Bo   , Bosol, Coriol, jwtype, &
    hbl  , bfsfc, stable, caseA, kbl  , &
    Rib  , sigma, wm    , ws )

!     the oceanic planetray boundary layer depth, hbl, is determined as
!     the shallowest depth where the bulk richardson number is
!     equal to the critical value, Ricr.

!     bulk richardson numbers are evaluated by computing velocity and
!     buoyancy differences between values at zgrid(kl) < 0 and surface
!     reference values.
!     in this configuration, the reference values are equal to the
!     values in the surface layer.
!     when using a very fine vertical grid, these values should be
!     computed as the vertical average of velocity and buoyancy from
!     the surface down to epsilon*zgrid(kl).

!     when the bulk richardson number at k exceeds Ricr, hbl is
!     linearly interpolated between grid levels zgrid(k) and zgrid(k-1).

!     The water column and the surface forcing are diagnosed for
!     stable/ustable forcing conditions, and where hbl is relative
!     to grid points (caseA), so that conditional branches can be
!     avoided in later subroutines.

		implicit none

!  model
    integer :: km,kmp1      ! number of vertical levels
    integer :: imt          ! number of horizontal grid points
    real :: zgrid(imt,kmp1) ! vertical grid (<= 0)              (m)
    real :: hwide(imt,kmp1) ! layer thicknesses                 (m)

!  input
    real :: dVsq(imt,kmp1)  ! (velocity shear re sfc)^2      (m/s)^2
    real :: dbloc(imt,km)   ! local delta buoyancy              (m/s^2)
    real :: Ritop(imt,km)   ! numerator of bulk Richardson Number (m/s)^2
!          Ritop = (-z - -zref)* delta buoyancy w/ respect to sfc(m/s^2)
    real :: ustar(imt)      ! surface friction velocity         (m/s)
    real :: Bo(imt)         ! surface turbulent buoyancy forcing(m^2/s^3)
    real :: Bosol(imt)      ! radiative buoyancy forcing        (m^2/s^3)
    real :: Coriol          ! Coriolis parameter                (1/s)
    integer :: jwtype(imt)  ! Jerlov water type                 (1 to 5)

!  output
    real :: hbl(imt)        ! boundary layer depth              (m)
    real :: bfsfc(imt)      ! Bo+radiation absorbed to d=hbf*hbl(m^2/s^3)
    real :: stable(imt)     ! =1 in stable forcing; =0 unstable
    real :: caseA(imt)      ! =1 in case A, =0 in case B
    integer :: kbl(imt)     ! index of first grid level below hbl

!  local
    real :: Rib(imt,3)      ! Bulk Richardson number
    real :: sigma(imt)      ! normalized depth (d/hbl)
    real :: wm(imt),ws(imt) ! turbulent velocity scales         (m/s)

    real, save :: epsln =  1.e-20
    real, save :: DelVmin =  .005
    real, save :: Ricr =  0.50
    real, save :: epsilon =  0.1
    real, save :: cekman =  0.7
    real, save :: cmonob =  1.0
    real, save :: cs = 98.96
    real, save :: cv =  2.00
    real, save :: vonk =  0.4
    real, save :: hbf =  1.0

		integer :: i,kl,kupper,kup,kdn,ktmp
		double precision :: Vtc,z_upper,z_up,bvsq,Vtsq,slope_up,a_co,b_co,c_co, &
			sqrt_Arg,hekman,hmonob,hlimit

! find bulk Richardson number at every grid level until > Ric

! note: the reference depth is -epsilon/2.*zgrid(i,k), but the reference
!       u,v,t,s values are simply the surface layer values,
!       and not the averaged values from 0 to 2*ref.depth,
!       which is necessary for very fine grids(top layer < 2m thickness)
! note: max values when Ricr never satisfied are
!       kbl(i)=km and hbl(i) -zgrid(i,km)

    Vtc =  cv * sqrt(0.2/cs/epsilon) / vonk**2 / Ricr

!     indices for array Rib(i,k), the bulk Richardson number.
    kupper =1
    kup    =2
    kdn    =3
    z_upper = 0.0

!     initialize hbl and kbl to bottomed out values
    do 100 i = 1,imt
        Rib(i,kupper) = 0.0
        Rib(i,kup)    = 0.0
        z_up      = zgrid(i,1)
        kbl(i)    = km
        hbl(i)    = -zgrid(i,km)
    100 END DO

    do 200 kl = 2,km

    !        compute bfsfc = sw fraction at hbf * zgrid
        call swfrac_imt(imt,hbf,zgrid(1,kl),jwtype,bfsfc)

        do 190 i = 1,imt

        !           use caseA as temporary array
            caseA(i)  = -zgrid(i,kl)

        !           compute bfsfc= Bo + radiative contribution down to hbf * hbl
            bfsfc(i)  = Bo(i) &
            + Bosol(i) * (1. - bfsfc(i))
            stable(i) = 0.5 + SIGN( 0.5, bfsfc(i) )
            sigma(i)  = stable(i) * 1. + (1.-stable(i)) * epsilon

        190 END DO

    !        compute velocity scales at sigma, for hbl= caseA = -zgrid(i,kl)
        call wscale(imt, sigma, caseA, ustar, bfsfc,   wm, ws)

        do 180 i = 1,imt

        !           compute the turbulent shear contribution to Rib
        ! MODA.1
            bvsq = &
            dbloc(i,kl  )  / (zgrid(i,kl  )-zgrid(i,kl+1))
        !    $       0.5* ( dbloc(i,kl-1) / (zgrid(i,kl-1)-zgrid(i,kl  ))+
        !    $              dbloc(i,kl  ) / (zgrid(i,kl  )-zgrid(i,kl+1)) )
            Vtsq = - zgrid(i,kl) * ws(i) * sqrt(abs(bvsq)) * Vtc

        !           compute bulk Richardson number at new level, dunder
            Rib(i,kdn) = Ritop(i,kl) / (dVsq(i,kl)+Vtsq+epsln)
            if((kbl(i) == km) .AND. (Rib(i,kdn) >= Ricr)) then

            !              linear interpolate to find hbl where Rib = Ricr
            !              ka = kup
            !              ku = kdn
            !              hbl(i) = -zgrid(i,kl-1) + (zgrid(i,kl-1)-zgrid(i,kl)) *
            !    $                  (Ricr - Rib(i,ka)) / (Rib(i,ku)-Rib(i,ka))

            ! MODQ.1
            !               quadratic interpolation
                slope_up =  ( Rib(i,kupper) - Rib(i,kup) ) &
                / ( z_up - z_upper )
                a_co = ( Rib(i,kdn) - RIB(i,kup) &
                + slope_up * ( zgrid(i,kl) - z_up ) ) &
                / ( (z_up - zgrid(i,kl))**2 )
                b_co = slope_up + 2. * a_co * z_up
                c_co = RIB(i,kup) + z_up * (a_co*z_up + slope_up) &
                - Ricr
                sqrt_arg = b_co**2-4.*a_co*c_co
                if ( ( (abs(b_co) > epsln) .AND. &
                (abs(a_co)/abs(b_co) <= epsln) ) .OR. &
                (sqrt_arg            <= 0.0) ) then
                !              if(.true.) then
                    HBL(i) = -z_up + (z_up - zgrid(i,kl)) * &
                    (Ricr         - RIB(i,kup))/ &
                    (RIB(i,kdn) - RIB(i,kup))
                    if( .FALSE. ) then
                        write(6,*) 'kpp linear', z_upper, z_up,hbl(i),zgrid(i,kl)
                        write(6,*) Rib(i,kupper),Rib(i,kup),Rib(i,kdn)
                        write(6,*) slope_up,a_co,b_co,c_co,sqrt_arg
                        write(6,*) kupper,kup,kdn,kl,kbl(i)
                        stop
                    endif
                else
                    hbl(i) = (-b_co + sqrt(sqrt_arg)) / (2.*a_co)
                    if( .FALSE. ) then
                        write(6,*) 'kpp quad', z_upper, z_up,hbl(i),zgrid(i,kl)
                        write(6,*) Rib(i,kupper),Rib(i,kup),Rib(i,kdn)
                        write(6,*) slope_up,a_co,b_co,c_co,sqrt_arg,epsln
                        write(6,*) kupper,kup,kdn,kl,kbl(i)
                        stop
                    endif
                endif

                kbl(i) = kl
            endif
        180 END DO
    !-----------------------------------------------------------------------

    !       swap klevel indices and move to next level

    !-----------------------------------------------------------------------

        ktmp   = kupper
        kupper = kup
        kup    = kdn
        kdn    = ktmp
        z_upper = z_up
        z_up    = zgrid(1,kl)

    200 END DO

! compare hbl to limits
    call swfrac_imt(imt,-1.0,hbl,jwtype,bfsfc)

    do 300 i = 1,imt
        bfsfc(i)  = Bo(i) &
        + Bosol(i) * (1. - bfsfc(i))
        stable(i) = 0.5 + SIGN( 0.5, bfsfc(i) )
        bfsfc(i)  = bfsfc(i) + stable(i) * epsln !ensures bfsfc never=0
    300 END DO

!                          Option to ignore limits
    if( .TRUE. ) then

    !             check for HEKman or hmonob (11/22/94)
        do 400 i = 1,imt
            if(bfsfc(i) > 0.0) THEN
                hekman = cekman * ustar(i) / (abs(Coriol)+epsln)
            !         hekman = AMAX1(hekman,-zgrid(i,1))
                hmonob = cmonob * ustar(i)*ustar(i)*ustar(i) &
                /vonk / (bfsfc(i)+epsln)
            !         hmonob = AMAX1(hmonob,-zgrid(i,1))
                hlimit = stable(i)     * AMIN1(hekman,hmonob) + &
                (stable(i)-1.) * zgrid(i,km)
                hbl(i) = AMIN1(hbl(i),hlimit)
                hbl(i) = AMAX1(hbl(i),-zgrid(i,1))
            ENDIF
        !         hbl(i)= AMAX1(hbl(i),10.0)
            kbl(i) = km
        400 END DO
    !                find new kbl
        do 405 kl=2,km
            do 415 i =1,imt
                if((kbl(i) == km) .AND. (-zgrid(i,kl) > hbl(i))) then
                    kbl(i) = kl
                endif
            415 END DO
        405 END DO

    ! find stability and buoyancy forcing for final hbl values
        call swfrac_imt(imt,-1.0,hbl,jwtype,bfsfc)
        do 500 i = 1,imt
            bfsfc(i)  = Bo(i) &
            + Bosol(i) * (1. - bfsfc(i))
            stable(i) = 0.5 + SIGN( 0.5, bfsfc(i) )
            bfsfc(i)  = bfsfc(i) + stable(i) * epsln
        500 END DO

    endif          ! LIMITS

! determine caseA and caseB
    do 600 i = 1,imt
        caseA(i)  = 0.5 + &
        SIGN( 0.5,-zgrid(i,kbl(i)) -0.5*hwide(i,kbl(i)) -hbl(i))
    600 END DO

    return
    end subroutine bldepth

! *********************************************************************

    subroutine wscale(imt, sigma, hbl, ustar, bfsfc, &
    wm , ws   )

!     compute turbulent velocity scales.
!     use a 2D-lookup table for wm and ws as functions of ustar and
!     zetahat (=vonk*sigma*hbl*bfsfc).

		implicit none


! lookup table
    integer, parameter :: ni = 890            !  number of values for zehat
    integer, parameter :: nj = 48             ! number of values for ustar

    real, save :: wmt(0:ni+1,0:nj+1)           ! lookup table for wm
    real, save :: wst(0:ni+1,0:nj+1)           ! lookup table for ws
    real, save :: deltaz                       ! delta zehat in table
    real, save :: deltau                       ! delta ustar in table
    real, save :: zmin = -4.e-7 ! zehat limits for table
    real, save :: zmax = 0.0  ! m3/s3
    real, save :: umin = 0.   ! ustar limits for table
    real, save :: umax = .04  ! m/s
    logical, save :: firstf = .true.

!  model
    integer :: imt          ! number of horizontal grid points

!  input
    real :: sigma(imt)      ! normalized depth (d/hbl)
    real :: hbl(imt)        ! boundary layer depth (m)
    real :: ustar(imt)      ! surface friction velocity         (m/s)
    real :: bfsfc(imt)      ! total surface buoyancy flux       (m^2/s^3)

!  output
    real :: wm(imt),ws(imt) ! turbulent velocity scales at sigma

! local
    real :: zehat           ! = zeta *  ustar**3
    real :: zeta            ! = stability parameter d/L

    real, save :: epsln =  1.0e-20
    real, save :: c1 =  5.0
    real, save :: am = 1.257
    real, save :: cm = 8.380
    real, save :: c2 = 16.0
    real, save :: zetam = - 0.2
    real, save :: as = -28.86
    real, save :: cs = 98.96
    real, save :: c3 = 16.0
    real, save :: zetas = - 1.0
    real, save :: vonk = 0.40

    integer :: i,j,iz,ju,jup1,izp1
    real :: usta,zdiff,udiff,zfrac,ufrac,fzfrac,wam,wbm,was,wbs,ucube
    double precision :: temp_iz      ! add 4/24/2012, Saenz, to allow full column mixing with crash

! construct the wm and ws lookup tables

    if(firstf) then

        deltaz = (zmax-zmin)/(ni+1)
        deltau = (umax-umin)/(nj+1)

        do 100 i=0,ni+1
            zehat = deltaz*(i) + zmin
            do 90 j=0,nj+1
                usta = deltau*(j) + umin
                zeta = zehat/(usta**3+epsln)

                if(zehat >= 0.) then
                    wmt(i,j) = vonk*usta/(1.+c1*zeta)
                    wst(i,j) = wmt(i,j)
                else
                    if(zeta > zetam) then
                        wmt(i,j) = vonk* usta * (1.-c2*zeta)**(1./4.)
                    else
                        wmt(i,j) = vonk* (am*usta**3 - cm*zehat)**(1./3.)
                    endif
                    if(zeta > zetas) then
                        wst(i,j) = vonk* usta * (1.-c3*zeta)**(1./2.)
                    else
                        wst(i,j) = vonk* (as*usta**3 - cs*zehat)**(1./3.)
                    endif
                endif
            90 END DO
        100 END DO
        firstf=.false.
    endif

! use lookup table for zehat < zmax  ONLY;  otherwise use stable formulae
    do 200 i = 1,imt
        zehat = vonk * sigma(i) * hbl(i) * bfsfc(i)

        IF (zehat <= zmax) THEN

            zdiff  = zehat-zmin
            temp_iz = zdiff/deltaz
            if (abs(temp_iz) > ni) then
            	iz = ni
            else
            	iz = int( zdiff/deltaz )
            endif
            !iz = int( zdiff/deltaz ) ! replaced w/ if statement above, Saenz 4/12
            !iz = min( iz , ni )      ! replaced w/ if statement above, Saenz 4/12
            iz = max( iz , 0  )
            izp1=iz+1

            udiff  = ustar(i)-umin
            ju = int( udiff/deltau)
            ju = min( ju , nj )
            ju = max( ju , 0  )
            jup1=ju+1

            zfrac = zdiff/deltaz - float(iz)
            ufrac = udiff/deltau - float(ju)

            fzfrac= 1.-zfrac
            wam   = (fzfrac)  * wmt(iz,jup1) + zfrac*wmt(izp1,jup1)
            wbm   = (fzfrac)  * wmt(iz,ju  ) + zfrac*wmt(izp1,ju  )
            wm(i) = (1.-ufrac)* wbm          + ufrac*wam

            was   = (fzfrac)  * wst(iz,jup1) + zfrac*wst(izp1,jup1)
            wbs   = (fzfrac)  * wst(iz,ju  ) + zfrac*wst(izp1,ju  )
            ws(i) = (1.-ufrac)* wbs          + ufrac*was

        ELSE

            ucube = ustar(i)**3
            wm(i) = vonk * ustar(i) * ucube / (ucube + c1 * zehat)
            ws(i) = wm(i)

        ENDIF

    200 END DO

    return
    end subroutine wscale

! **********************************************************************

    subroutine ri_iwmix ( km  , kmp1, imt   , &
    Shsq, dbloc , zgrid , &
    visc, difs, dift  )

!     compute interior viscosity diffusivity coefficients due to
!     shear instability (dependent on a local richardson number)
!     and due to background internal wave activity.

		implicit none


!  input
    real :: Shsq(imt,kmp1)    ! (local velocity shear)^2          (m/s)^2
    real :: dbloc(imt,km)     ! local delta buoyancy              (m/s^2)
    real :: zgrid(imt,kmp1)   ! vertical grid (<= 0)              (m)
    integer :: km,kmp1        ! number of vertical levels
    integer :: imt            ! number of horizontal grid points

! output
    real :: visc(imt,0:kmp1)  ! vertical viscosivity coefficient  (m^2/s)
    real :: difs(imt,0:kmp1)  ! vertical scalar diffusivity       (m^2/s)
    real :: dift(imt,0:kmp1)  ! vertical temperature diffusivity  (m^2/s)

! local variables
		integer :: i,ki
    real :: Rig,Rigg          ! local richardson number
    real :: fri               ! function of Rig

    !save epsln,Riinfty,Ricon,difm0,difs0,difmcon,difscon, &
    !difmiw,difsiw,c1

    real, save ::   epsln = 1.e-16    ! a small number
    real, save ::   Riinfty = 0.8       ! default = 0.7
    real, save ::   Ricon = -0.2      ! note: exp was repl by multiplication
    real, save ::   difm0 = 0.005    ! max visc due to shear instability
    real, save ::   difs0 = 0.005    ! max diff ..  .. ..    ..
    real, save ::   difmiw = 0.0001    ! background/internal waves visc(m^2/s)
    real, save ::   difsiw = 0.00001   ! ..         ..       ..    diff(m^2/s)
    real, save ::   difmcon = 0.0000     ! max visc for convection  (m^2/s)
    real, save ::   difscon = 0.0000     ! max diff for convection  (m^2/s)
    real, save ::   c1 = 1.0
    real, save ::   c0 = 0.0
    integer, save ::   mRi = 2                ! number of vertical smoothing passes

    integer :: mr
    real :: ratio,fcon

!     compute interior gradient Ri at all interfaces, except surface

!-----------------------------------------------------------------------
!     compute interior gradient Ri at all interfaces ki=1,km, (not surface)
!       use visc(imt,ki=1,km) as temporary storage to be smoothed
!       use dift(imt,ki=1,km) as temporary storage of unsmoothed Ri
!       use difs(imt,ki=1,km) as dummy in smoothing call

    do 110 ki = 1, km
        do 100 i = 1,imt
            Rig  = dbloc(i,ki) * (zgrid(i,ki)-zgrid(i,ki+1)) / &
            ( Shsq(i,ki) + epsln)
            dift(i,ki) = Rig
            visc(i,ki) = dift(i,ki)
        100 END DO
    110 END DO

!-----------------------------------------------------------------------
!     vertically smooth Ri mRi times
    do mr = 1,mRi
        call z121(kmp1,imt,c0,Riinfty,visc,difs)
    enddo

!-----------------------------------------------------------------------
!                           after smoothing loop
    do 210 ki = 1, km
        do 200 i = 1,imt
        !  evaluate f of unsmooth Ri (fri) for convection        store in fcon
        !  evaluate f of   smooth Ri (fri) for shear instability store in fri

            Rigg  = AMAX1( dift(i,ki) , Ricon )
            ratio = AMIN1( (Ricon-Rigg)/Ricon , c1 )
            fcon  = (c1 - ratio*ratio)
            fcon  = fcon * fcon * fcon

            Rigg  = AMAX1( visc(i,ki) , c0 )
            ratio = AMIN1( Rigg/Riinfty , c1 )
            fri   = (c1 - ratio*ratio)
            fri   = fri * fri * fri

        !      if(i.eq.1)write(6,*)ki,dift(i,ki),visc(i,ki),Rigg,ratio,fri

        ! ************************   Overwrite with Gent's PP **********
        !           fcon = 0.0
        !           Rigg  = AMAX1( dift(i,ki) , c0 )
        !           fri   = c1 / (c1 + 10. * Rigg )
        !           difm0 = 0.1 * fri
        !           difs0 = 0.1 * fri * fri

        !  ************************   Overwrite with original PP
        !           fcon = 0.0
        !           Rigg  = AMAX1( dift(i,ki) , c0 )
        !           fri   = c1 / (c1 +  5. * Rigg )
        !           difm0 = 0.01 * fri
        !           difs0 = (difmiw + fri * difm0)

        ! ----------------------------------------------------------------------
        !            evaluate diffusivities and viscosity
        !    mixing due to internal waves, and shear and static instability

            visc(i,ki) = (difmiw + fcon * difmcon + fri * difm0)
            difs(i,ki) = (difsiw + fcon * difscon + fri * difs0)
            dift(i,ki) = difs(i,ki)

        200 END DO
    210 END DO

! ------------------------------------------------------------------------
!         set surface values to 0.0

    do i = 1,imt
        visc(i,0)    = c0
        dift(i,0)    = c0
        difs(i,0)    = c0
    enddo

    return
    end subroutine ri_iwmix

! *********************************************************************
    Subroutine z121 (kmp1,imt,vlo,vhi,V,w)

!    Apply 121 smoothing in k to 2-d array V(i,k=1,km)
!     top (0) value is used as a dummy
!     bottom (kmp1) value is set to input value from above.

!  input
    real :: V(imt,0:kmp1)  ! 2-D array to be smoothed in kmp1 direction
    real :: w(imt,0:kmp1)  ! 2-D array of internal weights to be computed

    km  = kmp1 - 1

    do i=1,imt
        w(i,0)    =   0.0
        w(i,kmp1) =   0.0
        V(i,0)    =   0.0
        V(i,kmp1) =   0.0
    enddo

    do k=1,km
        do i=1,imt
            if((V(i,k) < vlo) .OR. (V(i,k) > vhi)) then
                w(i,k) = 0.0
            !           w(i,k) = 1.0
            else
                w(i,k) = 1.0
            endif
        enddo
    enddo

    do k=1,km
        do i=1,imt
            tmp    = V(i,k)
            V(i,k) = w(i,k-1)*V(i,0)+2.*V(i,k)+w(i,k+1)*V(i,k+1)
            wait   = w(i,k-1) + 2.0 + w(i,k+1)
            V(i,k) = V(i,k) / wait
            V(i,0) = tmp
        enddo
    enddo

    return
    end Subroutine z121
! *********************************************************************

    subroutine ddmix (km  , kmp1, imt   , &
    alphaDT, betaDS, zgrid, &
    visc, difs, dift  )

!     Rrho dependent interior flux parameterization.
!     Add double-diffusion diffusivities to Ri-mix values at blending
!     interface and below.

		implicit none

! input
		integer :: km,kmp1,imt
    real :: alphaDT(imt,kmp1)  ! alpha * DT  across interfaces
    real :: betaDS(imt,kmp1)   ! beta  * DS  across interfaces
    real :: zgrid(imt,kmp1)

! output
    real :: visc(imt,0:kmp1)  ! interior viscosity           (m^2/s)
    real :: dift(imt,0:kmp1)  ! interior thermal diffusivity (m^2/s)
    real :: difs(imt,0:kmp1)  ! interior scalar  diffusivity (m^2/s)

! local
		integer :: ki,i
    real :: Rrho              ! dd parameter
    real :: diffdd            ! double diffusion diffusivity scale
    real :: prandtl           ! prandtl number

    !save Rrho0,dsfmax

    real, save :: Rrho0  =  2.55   ! Rp=(alpha*delT)/(beta*delS)
    real, save :: dsfmax = 1.0e-4  ! 0.0001 m2/s = 1 cm2/s

    do 100 ki= 1, km

        do 90 i = 1,imt

        !     salt fingering case
            if((alphaDT(i,ki) > betaDS(i,ki)) .AND. &
            (betaDS(i,ki) > 0.)) then
                Rrho  = MIN(alphaDT(i,ki) / betaDS(i,ki) , Rrho0)
                diffdd     =         1.0-((Rrho-1)/(Rrho0-1))
                diffdd     = dsfmax*diffdd*diffdd*diffdd
                dift(i,ki) = dift(i,ki) + 0.7*diffdd
                difs(i,ki) = difs(i,ki) + diffdd

            !     diffusive convection relative to molecular diffusivity 1.5e-6 m2/s
            else if ((betaDS(i,ki) < alphaDT(i,ki)) .AND. &
                (alphaDT(i,ki) < 0.0))  then
                Rrho    = alphaDT(i,ki) / betaDS(i,ki)
                diffdd  = 1.5e-6*9.0*0.101*exp(4.6*exp(-0.54*(1/Rrho-1)))
                prandtl = 0.15*Rrho
                if (Rrho > 0.5) prandtl = (1.85-0.85/Rrho)*Rrho
                dift(i,ki) = dift(i,ki) + diffdd
                difs(i,ki) = difs(i,ki) + prandtl*diffdd

            endif
        90 END DO
    100 END DO

    return
    end subroutine ddmix

! *********************************************************************

    subroutine blmix &
    (km   , kmp1 , imt , mdiff , zgrid, hwide , &
    ustar, bfsfc, hbl , stable, caseA, &
    visc , difs , dift, kbl   , &
    gat1 , dat1 , dkm1, blmc  , ghats, &
    sigma, wm   , ws  )
!     mixing coefficients within boundary layer depend on surface
!     forcing and the magnitude and gradient of interior mixing below
!     the boundary layer ("matching").

		implicit none

! UTION if mixing bottoms out at hbl = -zgrid(km) THEN
!   fictious layer kmp1 is needed with small but finite width (eg. 1.e-10)
! model
    integer :: km,kmp1        ! number of vertical levels
    integer :: imt            ! number of horizontal grid points
    integer :: mdiff          ! number of viscosities + diffusivities
    real :: zgrid(imt,kmp1)   ! vertical grid (<=0)               (m)
    real :: hwide(imt,kmp1)   ! layer thicknesses                 (m)

! input
    real :: ustar(imt)        ! surface friction velocity         (m/s)
    real :: bfsfc(imt)        ! surface buoyancy forcing        (m^2/s^3)
    real :: hbl(imt)          ! boundary layer depth              (m)
    real :: stable(imt)       ! = 1 in stable forcing
    real :: caseA(imt)        ! = 1 in case A
    real :: visc(imt,0:kmp1)  ! vertical viscosity coefficient    (m^2/s)
    real :: difs(imt,0:kmp1)  ! vertical scalar diffusivity       (m^2/s)
    real :: dift(imt,0:kmp1)  ! vertical temperature diffusivity  (m^2/s)
    integer :: kbl(imt)       ! index of first grid level below hbl

! output
    real :: gat1(imt,mdiff)
    real :: dat1(imt,mdiff)
    real :: dkm1(imt,mdiff)   ! boundary layer difs at kbl-1 level
    real :: blmc(imt,km,mdiff)! boundary layer mixing coefficients(m^2/s)
    real :: ghats(imt,km)     ! nonlocal scalar transport

!  local
    real :: sigma(imt)        ! normalized depth (d / hbl)
    real :: ws(imt), wm(imt)  ! turbulent velocity scales         (m/s)

    !save epsln,epsilon,c1,am,cm,c2,zetam,as,cs,c3,zetas, &
    !cstar,grav,vonk

    real, save :: epsln =  1.0e-20
    real, save :: epsilon =  0.1
    real, save :: c1 =  5.0
    real, save :: am = 1.257
    real, save :: cm = 8.380
    real, save :: c2 = 16.0
    real, save :: zetam = - 0.2
    real, save :: as = -28.86
    real, save :: cs = 98.96
    real, save :: c3 = 16.0
    real, save :: zetas = - 1.0
    real, save :: cstar = 5.
    real, save :: grav = 9.816
    real, save :: vonk = 0.40

		integer :: ki, i,kn
		real :: cg,delhat,R,dvdzup,dvdzdn,viscp,difsp,diftp,visch,difsh, &
			difth,f1,sig,a1,a2,a3,Gm,Gs,Gt


    cg = cstar * vonk * (cs * vonk * epsilon)**(1./3.)

! compute velocity scales at hbl
    do 100 i = 1,imt
        sigma(i) = stable(i) * 1.0 + (1.-stable(i)) * epsilon
    100 END DO

    call wscale(imt, sigma, hbl, ustar, bfsfc,   wm, ws)
    do 200 i = 1, imt
        kn    = ifix(caseA(i)+epsln) *(kbl(i) -1) + &
        (1-ifix(caseA(i)+epsln)) * kbl(i)

    ! find the interior viscosities and derivatives at hbl(i)
        delhat = 0.5*hwide(i,kn) - zgrid(i,kn) - hbl(i)
        R      = 1.0 - delhat / hwide(i,kn)
        dvdzup = (visc(i,kn-1) - visc(i,kn)) / hwide(i,kn)
        dvdzdn = (visc(i,kn)   - visc(i,kn+1)) / hwide(i,kn+1)
        viscp  = 0.5 * ( (1.-R) * (dvdzup + abs(dvdzup))+ &
        R  * (dvdzdn + abs(dvdzdn)) )

        dvdzup = (difs(i,kn-1) - difs(i,kn)) / hwide(i,kn)
        dvdzdn = (difs(i,kn)   - difs(i,kn+1)) / hwide(i,kn+1)
        difsp  = 0.5 * ( (1.-R) * (dvdzup + abs(dvdzup))+ &
        R  * (dvdzdn + abs(dvdzdn)) )

        dvdzup = (dift(i,kn-1) - dift(i,kn)) / hwide(i,kn)
        dvdzdn = (dift(i,kn)   - dift(i,kn+1)) / hwide(i,kn+1)
        diftp  = 0.5 * ( (1.-R) * (dvdzup + abs(dvdzup))+ &
        R  * (dvdzdn + abs(dvdzdn)) )

        visch  = visc(i,kn) + viscp * delhat
        difsh  = difs(i,kn) + difsp * delhat
        difth  = dift(i,kn) + diftp * delhat

        f1 = stable(i) * c1 * bfsfc(i) / (ustar(i)**4+epsln)

        gat1(i,1) = visch / hbl(i) / (wm(i)+epsln)
        dat1(i,1) = -viscp / (wm(i)+epsln) + f1 * visch
        dat1(i,1) = min(dat1(i,1),0.)

        gat1(i,2) = difsh  / hbl(i) / (ws(i)+epsln)
        dat1(i,2) = -difsp / (ws(i)+epsln) + f1 * difsh
        dat1(i,2) = min(dat1(i,2),0.)

        gat1(i,3) = difth /  hbl(i) / (ws(i)+epsln)
        dat1(i,3) = -diftp / (ws(i)+epsln) + f1 * difth
        dat1(i,3) = min(dat1(i,3),0.)

    !           Turn off interior matching here
    !       gat1(i,1) = 0.0001
    !       gat1(i,2) = 0.00001
    !       gat1(i,3) = 0.00001
    !       do m=1,3
    !       dat1(i,m) = 0.0
    !       enddo

    200 END DO

    do 300 ki = 1,km

    !     compute turbulent velocity scales on the interfaces

        do 290 i  = 1,imt
            sig     = (-zgrid(i,ki) + 0.5 * hwide(i,ki)) / hbl(i)
            sigma(i)= stable(i)*sig + (1.-stable(i))*AMIN1(sig,epsilon)
        290 END DO
        call wscale(imt, sigma, hbl, ustar, bfsfc,   wm,  ws)

    !     compute the dimensionless shape functions at the interfaces

        do 280 i = 1,imt
            sig = (-zgrid(i,ki) + 0.5 * hwide(i,ki)) / hbl(i)
            a1 = sig - 2.
            a2 = 3.-2.*sig
            a3 = sig - 1.

            Gm = a1 + a2 * gat1(i,1) + a3 * dat1(i,1)
            Gs = a1 + a2 * gat1(i,2) + a3 * dat1(i,2)
            Gt = a1 + a2 * gat1(i,3) + a3 * dat1(i,3)

        !     compute boundary layer diffusivities at the interfaces

            blmc(i,ki,1) = hbl(i) * wm(i) * sig * (1. + sig * Gm)
            blmc(i,ki,2) = hbl(i) * ws(i) * sig * (1. + sig * Gs)
            blmc(i,ki,3) = hbl(i) * ws(i) * sig * (1. + sig * Gt)

        !     nonlocal transport term = ghats * <ws>o
            ghats(i,ki) = (1.-stable(i)) * cg / (ws(i) * hbl(i) +epsln)

        280 END DO

    300 END DO

! find diffusivities at kbl-1 grid level
    do 400 i=1,imt
    	!print *,'zgrid - i: ',i,' kbl: ',kbl(i)
        sig      =  -zgrid(i,kbl(i)-1)  / hbl(i)
        sigma(i) =stable(i) * sig + (1.-stable(i)) * AMIN1(sig,epsilon)
    400 END DO

    call wscale(imt, sigma, hbl, ustar, bfsfc,   wm, ws)

    do 500 i = 1,imt
        sig = -zgrid(i,kbl(i)-1) / hbl(i)
        a1= sig - 2.
        a2 = 3.-2.*sig
        a3 = sig - 1.
        Gm = a1 + a2 * gat1(i,1) + a3 * dat1(i,1)
        Gs = a1 + a2 * gat1(i,2) + a3 * dat1(i,2)
        Gt = a1 + a2 * gat1(i,3) + a3 * dat1(i,3)
        dkm1(i,1) = hbl(i) * wm(i) * sig * (1. + sig * Gm)
        dkm1(i,2) = hbl(i) * ws(i) * sig * (1. + sig * Gs)
        dkm1(i,3) = hbl(i) * ws(i) * sig * (1. + sig * Gt)
    500 END DO

    return
    end subroutine blmix

! ******************************************************************

    subroutine enhance (km   , kmp1  , imt   , mdiff , dkm1  , visc , &
    difs , dift  , hbl   , kbl   , zgrid , caseA, &
    blmc , ghats )

! enhance the diffusivity at the kbl-.5 interface

		implicit none

! input
    integer :: km,kmp1        ! number of vertical levels
    integer :: imt            ! number of horizontal grid points
    integer :: mdiff          ! number of viscosities + diffusivities
    integer :: kbl(imt)       ! grid above hbl
    real :: hbl(imt)          ! boundary layer depth             (m)
    real :: dkm1(imt,mdiff)   ! bl diffusivity at kbl-1 grid level
    real :: zgrid(imt,kmp1)   ! vertical grid (<= 0)             (m)
    real :: visc(imt,0:kmp1)  ! enhanced viscosity               (m^2/s)
    real :: difs(imt,0:kmp1)  ! enhanced thermal diffusivity     (m^2/s)
    real :: dift(imt,0:kmp1)  ! enhanced scalar  diffusivity     (m^2/s)
    real :: caseA(imt)        ! = 1 in caseA, = 0 in case B

! input/output
    real :: ghats(imt,km)     ! nonlocal transport               (s/m**2)
!                              modified ghats at kbl(i)-1 interface
! output
    real :: blmc(imt,km,mdiff)! enhanced bound. layer mixing coeff.

! local
		integer :: ki,i
		real :: dkmp5,dstar
    real :: delta             ! fraction hbl lies beteen zgrid neighbors


    do 100 ki=1,km-1
        do 90 i = 1,imt

            if(ki == (kbl(i) - 1) ) then

                delta = (hbl(i)+zgrid(i,ki)) / (zgrid(i,ki)-zgrid(i,ki+1))

                dkmp5 = caseA(i) * visc(i,ki) + (1.-caseA(i)) * blmc(i,ki,1)
                dstar = (1.-delta)**2 * dkm1(i,1) + delta**2 * dkmp5
                blmc(i,ki,1) = (1.-delta) * visc(i,ki) + delta * dstar

                dkmp5 = caseA(i) * difs(i,ki) + (1.-caseA(i)) * blmc(i,ki,2)
                dstar = (1.-delta)**2 * dkm1(i,2) + delta**2 * dkmp5
                blmc(i,ki,2) = (1.-delta) * difs(i,ki) + delta * dstar

                dkmp5 = caseA(i) * dift(i,ki) + (1.-caseA(i)) * blmc(i,ki,3)
                dstar = (1.-delta)**2 * dkm1(i,3) + delta**2 * dkmp5
                blmc(i,ki,3) = (1.-delta) * dift(i,ki) + delta * dstar

                ghats(i,ki) = (1.-caseA(i)) * ghats(i,ki)

            endif

        90 END DO
    100 END DO

    return
    end subroutine enhance

!***********************************************************************
