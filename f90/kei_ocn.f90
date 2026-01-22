    SUBROUTINE  ocnstep (U,X,kforce,absorp_in)
!-----------------------------------------------------------------------
! Note in this version:
!   -  ADVECTIVE CORRECTIONS AVAILABLE
!   - read forcing, but compute SW from okta model   (fread in atmrad.f)
!   - use okta model for PAPA          (use cloudpapa in solar in bio.f)
!   - Jerlov water type II                            (SWDK in fluxes.f)
!   - albedo for ocean is set to 0.06, which is used for QSW from fcomp
!        or fread when rad/conv is not running:
!        albocn=0.06                             (init cnsts in input.f)
!   - no net fresh water flux into ocean when forcing with state
!        variables (laflx >= 1):
!        sflux(7,2,jptr) = - sflux(6,2,jptr)        (atmflx in fluxes.f)
!   - use psnow flux data to input observed SST,
!                                                    (fread in atmrad.f)
!-----------------------------------------------------------------------

!     Main driver for ocean module.
!     Integration is performed only on the permanent grid
!     Written   3 Mar 1991 - WGL
!     Modified  5 Jun 1992 - jan : implicit scheme
!              16 Nov      - jan : latest version
!              16 Nov 1994 - wgl : new KPP codes no temporary grid

    use kei_parameters
    use kei_common
		use kei_ocncommon

		implicit none

    !integer, parameter :: imt = NX*NY

! Input/Output
    real :: U(NZP1,NVEL), X(NZP1,NSCLR)
    real :: kforce(forcing_var_cnt)
    real :: absorp_in(NZ)

! Local Common Blocks -- relocated to module

!    real :: hmixd(0:1),         & !  storage arrays for extrapolations
!    Us(NZP1,NVEL ,0:1), & !  ..      ..     ..  ..
!    Xs(NZP1,NSCLR,0:1)  ! ..      ..     ..  ..
!    integer :: old,new          ! extrapolation index for Us,Xs,hmixd
!    common/ saveUXh / &
!    old,new,Us,Xs,hmixd

! Local

    real :: Un(NZP1,NVEL),      & !  new profiles
    Xn(NZP1,NSCLR),     & !  ..  ..
    hmixe,              & ! estimated hmix (integration input )
    hmixn,              & ! new comp. hmix (    ..      output)
    tol                   ! tolerance in hmix iteration
    integer :: &
    iter                ! number of iterations

	integer :: k,l,n,kmixn,kmixe,jwtype
	real :: deltaz,rhonot,sw_frac

! Estimate new profiles by extrapolation
    do 20 k=1,NZP1
        do 22 l=1,NVEL
            Un(k,l)=2.*Us(k,l,new)-Us(k,l,old)
        22 END DO
        do 24 l=1,NSCLR
            Xn(k,l)=2.*Xs(k,l,new)-Xs(k,l,old)
        24 END DO
    20 END DO

! Don't Estimate new profiles by extrapolation
!	Un = Us(:,:,new)
!	Xn = Xs(:,:,new)


!              Iteration loop for semi-implicit integration of KPP
! Reset iteration counter
    iter = 0
    call vmix(Un,Xn,hmixe,kmixe)
    call ocnint(NZ,zm,hm,dm,1,kmixe,U,X,Un,Xn)
    iter = 1
    IF (LKPP) THEN
        45 continue
        call vmix(Un,Xn,hmixn,kmixn)
        call ocnint(NZ,zm,hm,dm,1,kmixe,U,X,Un,Xn)
        iter = iter + 1

    !   check iteration for convergence
        tol = hmixtolfrac*hm(kmixn)
        if(kmixn == NZP1) tol = hmixtolfrac*hm(NZ)
        if(abs(hmixn-hmixe) > tol)then
            if (iter < itermax) then
                hmixe = hmixn
                kmixe = kmixn
                goto 45
            else
            !          use shallower hmix
                if(hmixn > hmixe) then
                    hmixe = hmixn                  ! comment out for hmix data
                    kmixe = kmixn                  ! ..      ..  ..  hmix data
                    goto 45                        ! ..      ..  ..  hmix data
                endif
            endif
        endif
        if( iter > (itermax+1) ) then
            write(6,1009) time,ntime,      & !  comment out for hmix data
            hmixe,hmixn,hmixn-hmixe,kmixn,iter
            1009 format('  long iteration at',f11.4,' days =',i6,' steps',/, &
            '  hmixest=',f7.2,' hmixnew=',f7.2,' diff=',f6.1, &
            ' kmixn=',i3,' iteration=',i3)
        endif
    ENDIF

!  Output  Results from permanent grid iterations to common.inc
! Compute diagnostic fluxes for writing to dat file
    do 40 k=1,NZ
        deltaz = 0.5*(hm(k)+hm(k+1))
        do 41 n=1,NSCLR
            wX(k,n)=-difs(k)*((Xn(k,n)-Xn(k+1,n))/deltaz- ghat(k)*wX(0,n))
        41 END DO
        if(LDD) &
        wX(k,1)= -dift(k)*((Xn(k,1)-Xn(k+1,1))/deltaz- ghat(k)*wX(0,1))
        wX(k,nsp1)= grav * ( talpha(k)*wX(k,1) - sbeta(k) * wX(k,2) )
        do 42 n=1,NVEL
            wU(k,n)= -difm(k)* (Un(k,n)-Un(k+1,n))/deltaz
        42 END DO
    40 END DO


! Compute energetics
    rhonot = 1026.
    Eflx = 0.5 * ( (U(1,1) + Un(1,1)) * sflux(1,5,0) + &
    (U(1,2) + Un(1,2)) * sflux(2,5,0) )
    Esnk = -0.5* rhonot  * ( (U(NZ,1) + Un(NZ,1)) * wU(NZ,1) + &
    (U(NZ,2) + Un(NZ,2)) * wU(NZ,2) )
    Ptke = 0.0
!     use "amax1" to prevent "underflow" in single precision
    do 120 k=1,NZ-1
        Ptke= Ptke- 0.5*( amax1(wU(k,1),1.E-10)* &
        (rhonot   * (U(k  ,1) + Un(k  ,1)) - &
        rhonot   * (U(k+1,1) + Un(k+1,1)) ) + &
        amax1(wU(k,2),1.E-10)* &
        (rhonot   * (U(k  ,2) + Un(k  ,2)) - &
        rhonot   * (U(k+1,2) + Un(k+1,2)) ) )
    120 END DO
    Tmke = 0.0
    do 130 k=1,NZP1
        rmke(k) = 0.5 * rhonot * (Un(k,1)**2 + Un(k,2)**2) * hm(k)
        Tmke    = Tmke + rmke(k)
    130 END DO

		! absorbed shortwave irradiance (W m-2) (assuming no scattering back out...)
		qsw_absorbed(1) = 1.0  ! init to total shortwave fraction
    jwtype = jerlov  ! ben added in 2021!
    do k=1,NZ
				call swfrac(1,1.0,zm(k),jwtype,sw_frac)
				!qsw_absorbed(k) = (qsw_absorbed(k) - sw_frac) * &
				!	kforce(qswins_f_ind)
				qsw_absorbed(k) = sw_frac * kforce(qswins_f_ind) ! not absorbed, just incident
				if (k .lt. NZ) then
					qsw_absorbed(k+1) = sw_frac
				endif
    enddo

!     check heat and salt budgets
!     call budget(X,Xn)

! Set new profiles
    do 30 k=1,NZP1         ! values at NZP1 only change for slab ocean
        do 31 n=1,NVEL
            U(k,n) = Un(k,n)
        31 END DO
        do 32 n=1,NSCLR
            X(k,n) = Xn(k,n)
        32 END DO
    30 END DO

! Set correct surface values, and CP and rho profiles for new profiles
! Get latest profiles
    call vmix(U,X,hmixn,kmixn)
    hmix = hmixn
    kmix = kmixn
    uref = U(1,1)
    vref = U(1,2)
    Tref = X(1,1)
    SSref= X(1,2)

! Save variables for next timestep
    old = new
    new = 1 - old
    hmixd(new) = hmix
    do 60 k=1,NZP1
        do 61 l=1,NVEL
            Us(k,l,new)=U(k,l)
        61 END DO
        do 62 l=1,NSCLR
            Xs(k,l,new)=X(k,l)
        62 END DO
    60 END DO

    return
    END SUBROUTINE

!**********************************************************************

    SUBROUTINE ocnint(nzi,z,h,d,intri,kmixe,Uo,Xo,Un,Xn)

!     Integrate the ocn model by backwards Euler(implicit)discretization
!     On input : Un,Xn are estimated profiles which are used
!                to estimate diffusivity profiles at new time.
!              : Updated diffusivities from Un Xn are in common
!     On output: Un,Xn are new profiles after integration.

!     Written  19 March 1991 - jan

    use kei_parameters
    use kei_common
		use kei_ocncommon

		implicit none

! Input
    real ::  Uo(nzi+1,NVEL), Xo(nzi+1,NSCLR), & !  old profiles
    z(nzi+1),h(nzi+1),d(0:nzi)
    integer :: intri,                         & !  index for tri.diag. coeff
    nzi, &
    kmixe
! Output
    real ::  Un(nzi+1,NVEL), Xn(nzi+1,NSCLR)  ! new profiles

! Common tridiagonal matrix factors (set in "init ocn") -- relocated to module
!    real :: tri(0:NZtmax,0:1,NGRID)    ! dt/dz/dz factors in trid. matrix
!    common/ trifac / tri
! Local
    real :: cu (NZtmax), & ! upper coeff for (k-1) on k line of trid.matrix
    cc (NZtmax), & !  central ...     (k  ) ..
    cl (NZtmax), & !  lower .....     (k-1) ..
    rhs(NZtmax)    ! right-hand-side terms
		integer :: i,n,imode
		real :: frot,ghatflux


! Compute biological fluxes: wXNT(0:nzi,3-nsclr)
    ! if ( LBIO ) call biomain(Xn)

! ********************************************************************
! U and V solution of tridiagonal matrix
!                               set coefficients of tridiagonal matrix
    call tridcof(difm,nzi,intri,cu,cc,cl)
!                               U right hand side and solution
    frot = f
    rhs(1)= Uo(1,1) + dto*(frot*.5*(Uo(1,2)+Un(1,2)) - wU(0,1)/h(1) )
    do 110 i=2,nzi-1
        rhs(i)= Uo(i,1) + dto*frot*.5*(Uo(i,2)+Un(i,2))
    110 END DO
    i=nzi   ! bottom
    rhs(i)= Uo(i,1) + dto*frot*.5*(Uo(i,2)+Un(i,2)) &
    + tri(i,1,intri)*difm(i)*Uo(i+1,1)
    call tridmat(cu,cc,cl,rhs,Uo(1,1),nzi,Un(1,1), &
    difm)

!                                 V rhs and solution
    rhs(1)= Uo(1,2) - dto*(frot*.5*(Uo(1,1)+Un(1,1)) + wU(0,2)/h(1) )
    do 120 i=2,nzi-1
        rhs(i)= Uo(i,2) - dto*frot*.5*(Uo(i,1)+Un(i,1))
    120 END DO
    i=nzi
    rhs(i)= Uo(i,2) - dto*frot*.5*(Uo(i,1)+Un(i,1)) &
    + tri(i,1,intri)*difm(i)*Uo(i+1,2)
    call tridmat(cu,cc,cl,rhs,Uo(1,2),nzi,Un(1,2),difm)

! *******************************************************************
! Scalar solutions of tridiagonal matrix
!     Temperature (different from other scalars because of ghat-term
!                  and double diffusion)
!     ghatflux = wX(0,1) - (1-SWDK(-hmixe,real(time)))
!    $                     * sflux(3,5,0) / rho(0) / CP(0)
!     ghatflux = wX(0,1) - (1-SWDK(-d(1) ,real(time)))
!    $                     * sflux(3,5,0) / rho(0) / CP(0)
    ghatflux = wX(0,1)
    call tridcof(dift,nzi,intri,cu,cc,cl)
    call tridrhs(h,Xo(1,1),wXNT(0,1),dift,ghat,wX(0,1),ghatflux, &
    dto,nzi,intri,rhs)
!                                     modify rhs for advection
    do imode=1,nmodeadv(1)
        call rhsmod(1,modeadv(imode,1),advection(imode,1), &
        time,dto,dpy,kmixe,d(kmixe),nzi,rho,cp,h,z,rhs)
    enddo
    call tridmat(cu,cc,cl,rhs,Xo(1,1),nzi,Xn(1,1),dift)

!     Salinity and other scalars
    call tridcof(difs,nzi,intri,cu,cc,cl)
    do 200 n=2,NSCLR
        ghatflux = wX(0,n)
        call tridrhs(h,Xo(1,n),wXNT(0,n),difs,ghat,wX(0,n),ghatflux, &
        dto,nzi,intri,rhs)
    !                                       modify rhs for advections
        do imode=1,nmodeadv(2)
            call rhsmod(2,modeadv(imode,2),advection(imode,2), &
            time,dto,dpy,kmixe,d(kmixe),nzi,rho,cp,h,z,rhs)
        enddo
        call tridmat(cu,cc,cl,rhs,Xo(1,n),nzi,Xn(1,n),difs)
    200 END DO

    return
    end SUBROUTINE ocnint

!***********************************************************************

    SUBROUTINE tridcof(diff,nzi,ind,cu,cc,cl)

!     Compute coefficients for tridiagonal matrix (dimension=nzi).
!     Note: cu(1) = 0. and cl(nzi) = 0. are necessary conditions.
!-----
    use kei_parameters
		use kei_ocncommon

		implicit none

! Input
    real :: diff(0:nzi)   ! diffusivity profile on interfaces
    integer :: nzi,     & ! dimension of field
    ind      ! index for tri-coefficients: = kmixo for t-grid,
!                                                    =     1 for p-grid.
! Output
    real :: cu(nzi),    & !  upper coeff. for (k-1) on k line of trid.matrix
    cc(nzi),    & !  central ...      (k  ) ..
    cl(nzi)     ! lower .....      (k-1) ..
! Common tridiagonal factors (set in "init ocn<") -- relocated to module
!    real :: tri(0:NZtmax,0:1,NGRID)    ! dt/dz/dz factors in trid. matrix
!    common/ trifac / tri

! Local
		integer :: i


! In the surface layer
    cu(1) = 0.
    cc(1) = 1. + tri(1,1,ind)*diff(1)   ! 1.+ dto/h(1)/dzb(1)*diff(1)
    cl(1) =    - tri(1,1,ind)*diff(1)   !   - dto/h(1)/dzb(1)*diff(1)
! Inside the domain
    do 10 i=2,nzi
        cu(i) =    - tri(i,0,ind)*diff(i-1)
        cc(i) = 1. + tri(i,1,ind)*diff(i)   + tri(i,0,ind)*diff(i-1)
        cl(i) =    - tri(i,1,ind)*diff(i)
    10 END DO
! In the bottom layer
    cl(nzi)= 0.
    return
    end SUBROUTINE tridcof

!**********************************************************************

    SUBROUTINE tridrhs(h,yo,ntflux,diff,ghat,sturflux,ghatflux, &
    dto,nzi,ind,rhs)

!     Compute right hand side of tridiagonal matrix for scalar fields:
!     =  yo (old field)
!      + flux-divergence of ghat
!      + flux-divergence of non-turbulant fluxes
!     Note: surface layer needs +dto/h(1) * surfaceflux
!           bottom  ..... ..... +dto/h(nzi)*diff(nzi)/dzb(nzi)*yo(nzi+1)

    use kei_parameters
		use kei_ocncommon

		implicit none


! Input
    real :: dto           ! timestep interval (seconds)
    real :: h(nzi+1),     & !  layer thickness
    yo(nzi+1),    & !  old profile
    ntflux(0:nzi),& !  non-turbulent flux = wXNT(0:nzi,1:2)
    diff(0:nzi),  & !  diffusivity profile on interfaces
    ghat(nzi),    & !  ghat turbulent flux
    sturflux,     & !  surface turbulent (kinematic) flux = wX(0,n)
    ghatflux      ! surface flux for ghat: includes solar flux
    integer :: nzi,       & !  dimension of field
    ind        ! index for tri-coefficients:=kmixo for t-grid,
!                                                     =    1 for p-grid
! Output
    real :: rhs(nzi)      ! right hand side
! Common tridiagonal factors (set in "init ocn") -- relocated to module
!    real :: tri(0:NZtmax,0:1,NGRID)    ! dt/dz/dz factors in trid. matrix
!    common/ trifac / tri

! Local
		integer :: i

! In the surface layer (dto/h(1)=tri(0,1,ind)
    rhs(1)= yo(1) + dto/h(1) * &
    ( ghatflux* diff(1)*ghat(1) - sturflux &
    +ntflux(1) - ntflux( 0 ) )
! Inside the domain
    do 10 i=2,nzi-1
        rhs(i)= yo(i) + dto/h(i) * &
        ( ghatflux*(diff(i)*ghat(i) - diff(i-1)*ghat(i-1)) &
        +ntflux(i) - ntflux(i-1) )
    10 END DO
! In the bottom layer
    if(nzi > 1) then   ! not for slab ocean
        i=nzi
        rhs(i)= yo(i) + dto/h(i) * &
        ( ghatflux*(diff(i)*ghat(i) - diff(i-1)*ghat(i-1)) &
        +ntflux(i) - ntflux(i-1) ) &
        + yo(i+1)*tri(i,1,ind)*diff(i)
    endif
    return
    end SUBROUTINE tridrhs

!**********************************************************************
    subroutine rhsmod(jsclr,mode,A,time, &
    dto,dpy,km,dm,nzi,rho,cp,h,z,rhs)

!     Modify rhs to correct scalar, jsclr,
!     for advection according to mode
! mode = 1 : Steady upper layer horizontal advection
!        2 : Steady mixed layer horizontal advection to km-1
!        3 : Steady horizontal advection throughout the entire column
!        4 : Steady vertical advection (= deep horizontal) below 100m
!            to bottom
!            (Change: start below 100m, instead of at layer 16, and
!            do not advect into bottom layer, 7-1-93)
!        5 : Steady bottom diffusion
!        6 : Seasonal mixed layer horizontal advection to dm
!        7 : Seasonal thermocline horizontal advection to 1.5 dm

! Input
    real :: h(nzi+1),     & !  layer thickness
    z(nzi+1),     & !  z grid levels (added as input on 7-1-93)
    rhs(nzi),     & !  right hand side from tridrhs
    rho(nzi),     & !  density
    cp (nzi)      ! specific heat
    DOUBLE PRECISION :: time        ! time in days from jan 1 of any year.
    real :: dto,          & !  ocean time step
    dpy,          & !  days per year (added as input on 7-1-93)
    dm,           & !  depth d(km+.5)
    A             ! advection of heat(W/m2) or Salt(PSU m/s)
    integer :: nzi,       & !  vertical dimension of field
    km,       & !  index of gridpoint just below h
    mode,       & !  type of advection
    jsclr        ! scalar

! Output
!     real rhs(nzi)      ! modified right hand side

! Internal
    real :: f(12) = &     ! monthly partion of annual advection
    	(/.05,.05,0.0,0.0,0.0,0.0,0.0,.05,.15,.20,.30,.20/)
    real :: xsA(21) =  &  ! yearly excess of heat
			(/48.26,21.73,29.02,56.59,19.94,15.96,18.28, &
    	40.52,37.06,29.83,29.47,15.77, 1.47,14.55, &
    	4.22,28.19,39.54,19.58,20.27,11.19,21.72/)
!     data f/.1,.1,6*0.0,.1,.3,.4,.2/
!    data f/.05,.05,5*0.0,.05,.15,.20,.30,.20/
!    data xsA/48.26,21.73,29.02,56.59,19.94,15.96,18.28, &
!    40.52,37.06,29.83,29.47,15.77, 1.47,14.55, &
!    4.22,28.19,39.54,19.58,20.27,11.19,21.72/

    if(mode <= 0) return
!               find ocean year
    iyr = 1 + idint((time-75.0)/dpy)
    day = time - dpy * (idint(time/dpy))  ! 365.25))
    month = 1 + int(12. * day / dpy)    ! 367.)
    if(month > 12) then
        write(6,*) 'STOP rhsmod (ocn.f):'
        write(6,*) '     rounding error, month gt 12 =',month
        stop 97
    endif

!       Am = -12. * f(month) * (xsA(iyr) - 0.0 )       ! Annual
    Am =  12. * f(month) * A                       ! Seasonal
!       Am = A                                          ! Steady

    if(mode == 1) then
    !                          correct upper layer advection
        if(jsclr == 1) fact = dto * Am / (rho(1)*cp(1))
        if(jsclr == 2) fact = dto * Am * 0.033
        rhs(1) = rhs(1) &
        + fact / h(1)

    else if(mode == 2) then
    !                          correct mixed layer advection
        delta = 0.0
        do 205 n=1,km-1
            delta = delta + h(n)
        205 END DO
        do 215 n=1,km-1
            if(jsclr == 1) fact = dto * Am / (rho(n)*cp(n))
            if(jsclr == 2) fact = dto * Am * 0.033
            rhs(n) = rhs(n) &
            + fact  / delta
        215 END DO

    else if (mode == 3) then
    !                               throughout whole water column
        delta = 0.0
        do 305 n=1,nzi
            delta = delta + h(n)
        305 END DO
        do 315 n=1,nzi
            if(jsclr == 1) fact = dto * Am / (rho(n)*cp(n))
            if(jsclr == 2) fact = dto * Am * 0.033
            rhs(n) = rhs(n) &
            + fact / delta
        315 END DO


    else if (mode == 4) then
    !                          vertical advection = deep horizontal
        nzend=nzi-1                          ! nzend=nzi (change:7-1-93)
        n1=0                                 ! n1=16     (change:7-1-93)
        401 n1=n1+1
        if(z(n1) >= -100.) goto 401
        delta = 0.0
        do 405 n=n1,nzend
            delta = delta + h(n)
        405 END DO
        do 415 n=n1,nzend
            if(jsclr == 1) fact = dto * Am / (rho(n)*cp(n))
            if(jsclr == 2) fact = dto * Am * 0.033
            rhs(n) = rhs(n) &
            + fact / delta
        415 END DO

    else if(mode == 5) then
    !                          correct bottom layer diffusion
        if(jsclr == 1) fact = dto * Am / (rho(nzi)*cp(nzi))
        if(jsclr == 2) fact = dto * Am * 0.033
        rhs(nzi) = rhs(nzi) &
        + fact / h(nzi)

    else

    !              seasonal mixed layer or thermocline advection
    !               find ocean year
    !       iyr = 1 + idint((time-75.0)/dpy)
    !       day = time - dpy * (idint(time/dpy))  ! 365.25))
    !       month = 1 + int(12. * day / dpy)    ! 367.)
    ! iag
    !       if(month.gt.12) then
    !          write(6,*) 'STOP rhsmod (ocn.f):'
    !          write(6,*) '     rounding error, month gt 12 =',month
    !          stop 97
    !       endif
    ! iag
    !       Am = -12. * f(month) * (xsA(iyr) - 0.0 )       ! Annual
    !       Am =  12. * f(month) * A                       ! Seasonal
    !       Am = A                                          ! Steady

        if(mode == 6) then
        !                            mixed layer to dm
            n1 = 1
            depth = h(1)
            dmax  = dm -  0.5 * (h(km) + h(km-1))
            delta = 0.0
            do 605 n =n1,nzi
                n2    = n
                delta = delta + h(n)
                depth = depth + h(n+1)
                if(depth >= dmax) go to 606
            605 END DO
            606 continue

        else if (mode == 7) then
        !                                 thermocline to 100m
            n1 = km - 1
            depth = dm - 0.5 * h(km)
            dmax = 100.
            delta = 0.0
            do 705 n=n1,nzi
                n2 = n
                delta = delta + h(n)
                depth = depth + h(n+1)
                if(depth >= dmax) go to 706
            705 END DO
            706 continue

        else
            write(6,*) 'STOP in rhsmod (ocn.f):'
            write(6,*) '      mode out of range, mode=',mode
            stop 96
        endif

    !                          Finish both 6 and 7 here
        do 615 n=n1,n2
            if(jsclr == 1) fact = dto * Am / (rho(n)*cp(n))
            if(jsclr == 2) fact = dto * Am * 0.033
            rhs(n) = rhs(n) + fact  / delta
        615 END DO

    endif

    return
    end subroutine rhsmod

!**********************************************************************

    SUBROUTINE tridmat(cu,cc,cl,rhs,yo,nzi,yn,diff)

!     Solve tridiagonal matrix for new vector yn, given right hand side
!     vector rhs. Note: yn(nzi+1) = yo(nzi+1).
!-----
    use kei_parameters
		!use kei_ocncommon

		implicit none

! Input
    real :: cu (nzi),    & !  upper coeff. for (k-1) on k line of tridmatrix
    cc (nzi),    & !  central ...      (k  ) ..
    cl (nzi),    & !  lower .....      (k-1) ..
    rhs(nzi),    & !  right hand side
    yo(nzi+1),   & !  old field
    diff(0:nzi)
    integer :: nzi       ! dimension of matrix
! Output
    real :: yn(nzi+1)    ! new field
! Local
    real :: gam(NZtmax), & !  temporary array for tridiagonal solver
    bet          ! ...
    integer :: i

! Solve tridiagonal matrix.
    bet   = cc(1)
!     yn(1) = (rhs(1) + tri(0,1,ind)*surflux) / bet    ! surface
    yn(1) =  rhs(1) / bet    ! surface
    do 21 i=2,nzi
        gam(i)= cl(i-1)/bet
        bet   = cc(i) - cu(i)*gam(i)
        if(bet == 0.) then
            write(6,*) '** algorithm for solving tridiagonal matrix fails'
            write(6,*) '** bet=',bet
            write(6,*) '** i=',i,' cc=',cc(i),' cu=',cu(i),' gam=',gam(i)
            bet=1.E-12
        !          Pause 3
        endif
    !        to avoid "Underflow" at single precision on the sun
        yn(i) =      (rhs(i)  - cu(i)  *yn(i-1)  )/bet
    !        yni   =      (rhs(i)  - cu(i)  *yn(i-1)  )/bet
    !        yn(i) = max( (rhs(i)  - cu(i)  *yn(i-1)  )/bet , 1.E-12 )
    !        if(yni.lt.0.)
    !    +   yn(i) = min( (rhs(i)  - cu(i)  *yn(i-1)  )/bet ,-1.E-12 )
    21 END DO

!     yn(nzi)  = (rhs(nzi)- cu(nzi)*yn(nzi-1)
!    +                    + tri(nzi,1,ind)*diff(nzi)*yo(nzi+1) )/bet
!                                                      ! bottom
    do 22 i=nzi-1,1,-1
        yn(i)  = yn(i) - gam(i+1)*yn(i+1)
    22 END DO
    yn(nzi+1) = yo(nzi+1)

    return
    end SUBROUTINE tridmat

!********************************************************************

    SUBROUTINE init_ocn(U,X)

!     Initialize ocean model:
!     Set coefficients for tridiagonal matrix solver.
!     Compute hmix and diffusivity profiles for initial profile.
!     Prepare for first time step.

    use kei_parameters
		use kei_common
		use kei_ocncommon

		implicit none

! Input
    real ::  U(NZP1,NVEL),X(NZP1,NSCLR)

! Common Blocks -- relocated to module
!    real :: tri(0:NZtmax,0:1,NGRID)! dt/dz/dz factors in trid. matrix
!    common/ trifac / tri

!    real :: hmixd(0:1),            & !  storage arrays for extrapolations
!    Us(NZP1,NVEL ,0:1),    & !  ..      ..     ..  ..
!    Xs(NZP1,NSCLR,0:1)     ! ..      ..     ..  ..
!    integer :: old,new             ! extrapolation index for Us,Xs,hmixd
!    common/ saveUXh / &
!    old,new,Us,Xs,hmixd

! Local variables
    real :: dzb(NZ)                ! diff. between grid-levels below z(j)
    integer :: k,n,l,kmix0
    real :: hmix0,deltaz


! Compute factors for coefficients of tridiagonal matrix elements.
!     tri(0     ,1,.........) : dt/h(1) factor for rhs flux
!     tri(k=1:NZ,0,.........) : dt/h(k)/ {dzb(k-1)=z(k-1)-z(k)=dzabove}
!     tri(k=1:NZ,1,.........) : dt/h(k)/ {dzb(k  )=z(k)-z(k+1)=dzbelow}

    do 10 k=1,NZ
        dzb(k)     = zm(k) - zm(k+1)
    10 END DO

    tri(0,1,1)    = dto/hm(1)
    tri(1,1,1)    = dto/hm(1)/dzb(1)
    do 20 k=2,NZ
        tri(k,1,1) = dto/hm(k)/dzb(k)
        tri(k,0,1) = dto/hm(k)/dzb(k-1)
    20 END DO

! Determine hmix for initial profile:
    call vmix(U,X,hmix0,kmix0)
    hmix = hmix0
    kmix = kmix0
    Tref = X(1,1)
    write(6,*) 'init ocn',hmix,kmix,Tref,f

! Evaluate initial fluxes (to write to output data file)

    do 30 k=1,NZ
        deltaz = 0.5*(hm(k)+hm(k+1))
        do 31 n=1,NSCLR
            wX(k,n)=-difs(k)*((X(k,n)-X(k+1,n))/deltaz- ghat(k)*wX(0,n))
        31 END DO
        if(LDD) &
        wX(k,1)= -dift(k)*((X(k,1)-X(k+1,1))/deltaz- ghat(k)*wX(0,1))
        wX(k,nsp1)= grav * ( talpha(k)*wX(k,1) - sbeta(k) * wX(k,2) )
        do 32 n=1,NVEL
            wU(k,n)= -difm(k)* (U(k,n)-U(k+1,n))/deltaz
        32 END DO
    30 END DO

! Prepare for first time step

!     indices for extrapolation
    old = 0
    new = 1
!     initialize array for extrapolating hmixd,Us,Xs
    hmixd(0) = hmix
    hmixd(1) = hmix
    do 60 k=1,NZP1
        do 62 l=1,NVEL
            Us(k,l,0)=U(k,l)
            Us(k,l,1)=U(k,l)
        62 END DO
        do 64 l=1,NSCLR
            Xs(k,l,0)=X(k,l)
            Xs(k,l,1)=X(k,l)
        64 END DO
    60 END DO

    return
    end SUBROUTINE init_ocn

!********************************************************************
    subroutine swfrac( imt, fact, z, jwtype, swdk )
!     compute fraction of solar short-wave flux penetrating to specified
!     depth (times fact) due to exponential decay in  Jerlov water type
!     reference : two band solar absorption model of simpson and
!     paulson (1977)

    parameter(nwtype=5) ! max number of different water types

!  model
    integer :: imt         ! number of horizontal grid points

!  input
    real :: fact           ! scale  factor to apply to depth array
    real :: z         ! vertical height ( <0.) for desired sw
!                           fraction                                 (m)
    integer :: jwtype ! index for jerlov water type

!  output
    real :: swdk      !  short wave (radiation) fractional decay

!  local
!     jerlov water type :  I       IA      IB      II      III
!                jwtype    1       2       3       4       5

    real, save :: rfac(nwtype) = (/  0.58 ,  0.62 ,  0.67 ,  0.77 ,  0.78 /)
    real, save :: a1(nwtype) = (/  0.35 ,  0.6  ,  1.0  ,  1.5  ,  1.4  /)
    real, save :: a2(nwtype) = (/ 23.0  , 20.0  , 17.0  , 14.0  ,  7.9  /)

    do 100 i = 1,imt

        swdk =      rfac(jwtype)  * exp(z*fact/a1(jwtype)) &
        + (1.-rfac(jwtype)) * exp(z*fact/a2(jwtype))

    100 END DO

    return
    end subroutine swfrac


!********************************************************************

!********************************************************************
    subroutine swfrac_imt( imt, fact, z, jwtype, swdk )
!     compute fraction of solar short-wave flux penetrating to specified
!     depth (times fact) due to exponential decay in  Jerlov water type
!     reference : two band solar absorption model of simpson and
!     paulson (1977)

    parameter(nwtype=5) ! max number of different water types

!  model
    integer :: imt         ! number of horizontal grid points

!  input
    real :: fact           ! scale  factor to apply to depth array
    real :: z(imt)         ! vertical height ( <0.) for desired sw
!                           fraction                                 (m)
    integer :: jwtype(imt) ! index for jerlov water type

!  output
    real :: swdk(imt)      !  short wave (radiation) fractional decay

!  local
!     jerlov water type :  I       IA      IB      II      III
!                jwtype    1       2       3       4       5

    real, save :: rfac(nwtype) = (/  0.58 ,  0.62 ,  0.67 ,  0.77 ,  0.78 /)
    real, save :: a1(nwtype) = (/  0.35 ,  0.6  ,  1.0  ,  1.5  ,  1.4  /)
    real, save :: a2(nwtype) = (/ 23.0  , 20.0  , 17.0  , 14.0  ,  7.9  /)

    do 100 i = 1,imt

        swdk(i) =      rfac(jwtype(i))  * exp(z(i)*fact/a1(jwtype(i))) &
        + (1.-rfac(jwtype(i))) * exp(z(i)*fact/a2(jwtype(i)))

    100 END DO

    return
    end subroutine swfrac_imt


!********************************************************************


    subroutine budget (X1,X2)

    use kei_parameters
		use kei_common

		implicit none

		! Inputs
		real :: X1(nzp1,nsclr), X2(nzp1,nsclr)

		! Local
		integer :: k
		real :: T1,T2,S1,S2,fltn,dhdt,dfdt,Qtop,Ftop,Qbot,Fbot,delt,fact, &
			rhs,diff

    T1 = 0.0
    T2 = 0.0
    S1 = 0.0
    S2 = 0.0
    fltn = 1. / float(nz)
    write(6,*)
    do 5 k=1,nz
        T1  =  T1 + X1(k,1)  * hm(k)
        T2  =  T2 + X2(k,1)  * hm(k)
        S1  =  S1 + X1(k,2)  * hm(k)
        S2  =  S2 + X2(k,2)  * hm(k)
        write(6,*) k,X1(k,1),X2(k,1),X1(k,2),X2(k,2)
    5 END DO
    dhdt = (T2 - T1) * rho(0) * CP(0) / dtsec
    dfdt = (S1 - S2) * rho(0) / Sref  / dtsec
    Qtop = sflux(3,5,0) + sflux(4,5,0)
    Ftop = sflux(6,5,0)
    Qbot = -rho(0) * CP(0) * (wX(nz,1) + wXNT(nz,1) )
    Fbot =  rhoh2o / Sref  *  wX(nz,2)

    write(6,*) 'heat ',Qtop,dhdt,Qbot
    write(6,*) 'salt ',Ftop,dfdt,Fbot
    write(6,*) dtsec,rho(0),CP(0),Sref,sflux(3,5,0),sflux(4,5,0)

    do 15 k=1,nz
        delt = (X2(k,1)-X1(k,1))
        fact = dtsec / hm(k)
        rhs  = fact * (wX(k,1)-wX(k-1,1) + wXNT(k,1)-wXNT(k-1,1) )
        diff = delt - rhs
        write(6,*) wX(k-1,1),wXNT(k-1,1),wX(k,1),wXNT(k,1)
        write(6,*) fact,delt,rhs,diff
    15 END DO

    return
    end subroutine budget

