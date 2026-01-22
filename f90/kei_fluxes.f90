 
    SUBROUTINE calflx(jptr)

!     Calculate atm fluxes at dtcal/2  into sflux(n,i,0), i=1,3
!     from flux computations at ndtld in sflux(n,i,j=1,NJDT)
! NB  for NJDT = 1, sflux(n,i,0) = sflux(n,i,1)

    use kei_parameters
    use kei_common

    implicit none

    ! inputs
    integer :: jptr

    ! local vars
    integer :: n,i,j,jx
    real :: bnum

!             OPTION I  Find linear regression extrapolation

    do n=1,nsflxs
        do i=1,3
            jx = jptr
            bnum  = 0.0
            sflux(n,i,0) = 0.0
            do j=1,NJDT
                sflux(n,i,0) = sflux(n,i,0) + sflux(n,i,jx)
                bnum = bnum + sflux(n,i,jx) * deltax(j)
                jx = jx - 1
                if(jx < 1) jx = NJDT
            ENDDO
            sflux(n,i,0) = sflux(n,i,0)/NJDT - xbar * bnum / denom
    ENDDO
    ENDDO

!                   OPTION II simple average
    if( .FALSE. ) then
        do n=1,nsflxs
            do i=1,3
                sflux(n,i,0) = 0.0
                do j=1,NJDT
                    sflux(n,i,0) = sflux(n,i,0) + sflux(n,i,j) / float(NJDT)
                ENDDO
            ENDDO
        ENDDO
    endif

    return
    end SUBROUTINE calflx

!*****************************************************************

    SUBROUTINE atmflx(jptr,timed)

!     Load all atm fluxes from present near surface state variables
!     i=1,3 ie. bottom of atm, atm-ocn, and atm-ice fluxes.

!     timed : time (days) to evaluate SW in "swokta" and Jerlov water
!             type in "SWDK"

    use kei_parameters
    use kei_common
    use kei_icecommon

    implicit none

    interface
     real pure function qsat(mode, TK)
     implicit none
     integer, intent(in) :: mode
     real, intent(in) :: TK
     end function qsat
    end interface

    ! Input
    integer :: jptr
    DOUBLE PRECISION :: timed

    ! local vars
    integer :: n
    real :: diurnal, Us, Vs, Ts1, du, dv, umag, Tk, To, psim10, psis10, &
      zoz0, zozt, zozq, rtcd, rtct, rtcq, cloudf, ea, Tdew, fofc, &
      tairk, QSWocn, QSWice, QSWdn, Tsst, QLWdown, DQDTice, Tice

    real :: ak = 0.0027
    real :: bk = 0.000142
    real :: ck1 = 0.0000764
    real :: f1 = 0.174
    real :: umin = 0.10
    real :: ustrmin = 1.e4

    !data ak,bk,ck1 / 0.0027  , 0.000142 , 0.0000764 /
    !data f1, umin, ustrmin / 0.174, 0.10, 1.E-4/

! Set flux index

    jptr = jptr + 1
    if(jptr > NJDT) jptr = 1

!         Enable diurnal cycle
    IF(LAFLX == 0) then
        diurnal = MAX(0.0,sin(3.1416*(2.*timed-0.5)))
        if(timed <= 1.0) then
        !       write(6,56) timed,diurnal
        ! 6     format(' diurnal',2f9.3)
        endif
    ELSE
        diurnal = 1.0
    ENDIF
    diurnal = 1.0


! When forcing specified by atmospheric state variables,
! calculate fluxes.

    IF (LAFLX >= 1) THEN
    !        atm to ocn fluxes  send weighted by focn to top of ocean (5)
        Us = 0.0
        Vs = 0.0
        if(lfluxSSTdat) then
            Ts1 = SSTvaf     ! Equivalenced in COMMON to a vaf
        else
            Ts1 = Tref
        endif
        du = uZ - Us
        dv = vZ - Vs
        umag = sqrt(du**2 + dv**2)
        BP = 100.                          ! kPa
        Qs = 0.98 * QSAT(0,TK0+Ts1)         ! kg/m3
    !         get humidity (kg/m3) according to COMMON Equivalences Qz or Tdew
        if((LAFLX /= 1) .AND. (LAFLX /= 6) ) then  ! use dewpoint
!            Qz = QSAT(0,TK0+Tdew)
            Qz = QSAT(0,TK0+Tdewdat)   ! changed Tdew->Tdewdat, since Tdew appear to be uninitialized/not correct - ben saenz 7/2011
        endif

        TK = TZ + TK0
        To = TK / (1. - TK * QZ * f1 / BP )
        rhoa = 1.29  * (tk0/To) * (bp/101.)   !kg/m3
        cpa  = 1004. * (1. + 9.0 * QZ/rhoa)
    !            use some present atm parameters to avoid iterating
    !           az is to be set to state variable height (m) in init cnsts
        ustara = sqrt(sqrt(sflux(1,2,0)**2 + sflux(2,2,0)**2)/rhoa)
        azeta = az*vonk*grav*(sflux(5,1,0)/rhoa/cpa/to + &
        sflux(6,1,0)*tk*f1/bp/1000.) / (ustara**3 + ustrmin**3)
        if(azeta > 1.0) azeta =  1.0
        if(azeta < -2.0) azeta = -2.0
        call fzol(azeta,psima,psisa)
        call fzol(azeta*10./az,psim10,psis10)
        u10a = umag + ustara * (alog(10./az) + psima - psim10) / vonk
        u10a = AMAX1(u10a,umin)

    !        compute ocean flux transfer coefficients
        zoz0 = az * exp(vonk/sqrt(ak/u10a + bk + ck1*u10a )) / 10.
        if(azeta > 0.0) then
            zozt = az / 2.2E-9
        else
            zozt = az / 4.9E-5
        endif
        zozq = az / 9.5E-5
        rtcd = vonk / ( alog(zoz0) - psima )
        rtct = vonk / ( alog(zozt) - psisa )
        rtcq = vonk / ( alog(zozq) - psisa )
        cdw  =  rtcd * rtcd
        ctw  =  rtcd * rtct
        cqw  =  rtcd * rtcq

        if(LAFLX == 4) then

        ! QSWins =  Shortwave Insolation
            call swokta(rlon,rlat,dpy,timed,cloudf,QSWins)

        ! Downwelling Longwave
        !     water vapor pressure at (tair) = saturation water vapor
        !     press at (tdew)
            ea = 10.**( 0.7859 + 0.03477*Tdew ) / ( 1. + 0.00412*Tdew )
        !     cloud correction (Bunker,1976, for 50N)
            fofc = 1. - 0.72*cloudf
            tairk = TZ + tk0
            QLWdwn = -1. * epsw * sbc * tairk**4 * &
            (0.39-0.05*sqrt(ea))*fofc &
            - 4.* epsw * sbc * tairk**3 * (Ts1 - TZ) &
            + epsw * sbc * (Ts1+TK0)**4

            QSWocn          = focn * (1.-albocn) * QSWins
            if(LICE) QSWice = fice * (1.-albice) * QSWins
        else

            QSWocn          = focn * (1.-albocn) * vaf(3)
            QLWdwn          = vaf(4)
            QSWice          = fice * (1.-albice) * vaf(3)

        endif ! LAFLX=4

        if(LAFLX == 6) then
            QSWup  =  -albocn * QSWins
            QLWup  =  -1. * epsw * sbc * (Ts1 +TK0)**4
        endif  ! LAFLX=6

    !                      load fraction weighted air-sea interface
        sflux(1,2,jptr) = focn * rhoa * cdw * umag * du
        sflux(2,2,jptr) = focn * rhoa * cdw * umag * dv
        sflux(3,2,jptr) = focn * (QSWins + QSWup)
        sflux(4,2,jptr) = focn * (QLWdwn + QLWup)
        sflux(5,2,jptr) = focn * rhoa * cpa * ctw * umag * (TZ-Ts1)
        sflux(6,2,jptr) = focn * cqw * umag * (qZ-Qs)
        sflux(7,2,jptr) = focn * Prain
        sflux(8,2,jptr) = focn * Psnow
        sflux(9,2,jptr) = focn * QLWup

        print *,sflux(1,2,jptr), "focn * rhoa * cdw * umag * du"
        print *,sflux(2,2,jptr), "focn * rhoa * cdw * umag * dv"
        print *,sflux(3,2,jptr), "focn * (QSWins + QSWup)"
        print *,sflux(4,2,jptr), "focn * (QLWdwn + QLWup)"
        print *,sflux(5,2,jptr), "focn * rhoa * cpa * ctw * umag * (TZ-Ts1)"
        print *,sflux(6,2,jptr), "focn * cqw * umag * (qZ-Qs)"
        print *,sflux(7,2,jptr), "focn * Prain"
        print *,sflux(8,2,jptr), "focn * Psnow"
        print *,sflux(9,2,jptr), "focn * QLWup"

    ! ice flux transfer coefficients
        if(LICE) then
        !               compute atm to ice fluxes unwieghted to top of ice
        !                              For now ignore stability dependencies !!!!
            du = uZ - uI
            dv = vZ - vI
            umag = sqrt(du**2 + dv**2 )
            qI0  =   QSAT(1,TK0 + TI0)
            cdi = 1.75E-3
            cti = cdi
            cqi = cdi

            sflux(1,3,jptr) =  cdi * rhoa * umag * du
            sflux(2,3,jptr) =  cdi * rhoa * umag * dv
            sflux(3,3,jptr) =  QSWins
            sflux(4,3,jptr) =  QLWdwn
            sflux(5,3,jptr) =  rhoa * cpa * cti * umag * (TZ-TI0)
            shs = rhoa * cpa * cti * umag  ! added Saenz 10/2011 - to incorporate temp dependence of sensible heat flux in heat balance
            sflux(6,3,jptr) =  cqi * umag * (qZ - qI0)
            sflux(7,3,jptr) =  Prain
            sflux(8,3,jptr) =  Psnow
            sflux(9,3,jptr) = -epsi * sbc * (TI0+TK0)**4

            DQDTice = qI0  * 5897.8 / (TI0+TK0)**2
            DHDTice = -umag * (rhoa * cpa * cti + SL * cqi * DQDTice)

        endif

    !        Total fluxes at bottom of atm.
    !        (Note: fluxes from land fractions are not included)

        do n=1,NSFLXS
            sflux(n,1,jptr) = sflux(n,2,jptr) + fice * sflux(n,3,jptr)
        ENDDO
        sflux(3,1,jptr) = QSWins          - fice * albice * QSWins &
        + focn * QSWup
        sflux(4,1,jptr) = QLWdwn +         sflux(9,2,jptr) &
        + fice *  sflux(9,3,jptr)

    ELSE  ! LAFLX<=0
    ! When forcing with specified fluxes:

    !      LAFLX = 0 : forcing with prescribed fixed fluxes
    !                  VAF(i),i=1,nsflxs  read in input
    !      LAFLX < 0 : forcing with time varying fluxes
    !                  VAF(i),i=1,nsflxs  computed in atm
    !         write(6,*) 'start atm',vaf
    !           vaf(3) is to be  net SW into ice ocn system
    !           vaf(4) is net LW into ice-ocn system
        QSWdn = diurnal*vaf(3)/(fice*(1.-albice)+focn*(1.-albocn))

    !       atm to ice and ocn stresses
        do n=1, 2
            sflux(n,1,jptr) =  VAF(n)
            sflux(n,2,jptr) =  focn * VAF(n)
            sflux(n,3,jptr) =  VAF(n)
        ENDDO

    !         redo but use appropriate fractions of SW flux
        sflux(3,1,jptr) = QSWdn * (1.-fice*albice -focn*albocn)
        sflux(3,2,jptr) = QSWdn * focn*(1.-albocn)
        sflux(3,3,jptr) = QSWdn

    !       set LW radiation terms    QLWnet = vaf(4)
        if(lfluxSSTdat) then
            Tsst = SSTdat
            Tice = SSTdat
        else
            Tsst = Tref
            Tice = Ti0
        endif

        sflux(9,2,jptr) = - focn * epsw * sbc * (Tsst +TK0)**4
        sflux(9,3,jptr) = -        epsi * sbc * (Tice +TK0)**4
        sflux(9,1,jptr) = sflux(9,2,jptr) + fice * sflux(9,3,jptr)

        QLWdown = vaf(4) -sflux(9,2,jptr) -fice *sflux(9,3,jptr)
        sflux(4,1,jptr) = QLWdown
        sflux(4,2,jptr) = focn * QLWdown
        sflux(4,3,jptr) = QLWdown
    !                                  zero both ocn and ice until ice sfc fixed
        sflux(4,2,jptr) = - sflux(9,2,jptr)
        sflux(4,3,jptr) = - sflux(9,3,jptr)
        sflux(4,1,jptr) = - sflux(9,1,jptr)

    !          put all sensible and latent heat flux directly into the ocean; none to sea-ice
        do n=5,6
            sflux(n,1,jptr) =  VAF(n)
            sflux(n,2,jptr) =  VAF(n)
            sflux(n,3,jptr) =  0.0
        ENDDO

    !           distribute precip between ocean and ice, bcs no effect on ice sfc energy balance
        do n=7,8
            sflux(n,1,jptr) =  VAF(n)
            sflux(n,2,jptr) =  focn * VAF(n)
            sflux(n,3,jptr) =  VAF(n)
        ENDDO

    ENDIF    ! (forcing with state variables or fluxes)

    return
    end SUBROUTINE atmflx

!****************************************************************

    subroutine o2iflx(X,jptr)

!      from the ocean temperature profile load the ocean-ice fluxes
!      sflux(7,5)= top of ocean non-turbulent frazil (S=Sfrazil) <=0
!      sflux(8,5)= top of ocean tubulent heat >= 0
!      sflux(9,5)= top of ocean turbulent salt <=0

!      sflux(7,4)= bottom of ice freshwater frazil
!      sflux(8,4)= Melt potential (W/m2) <= 0
!      sflux(9,4)= bottom of ice fazil salt flux (S=Sice) <=0

!      compute non-turbulent heat, wtI,
!                        and salt, wsI,  fluxes due to deep frazil ice

!      find frocn for ice fraction calculations

    use kei_parameters
    use kei_common
    use kei_icecommon

    implicit none

    interface
     FUNCTION SWDK(z,ftime)
     implicit none
     REAL, intent(IN) :: z,ftime
     REAL :: SWDK
     end function SWDK
    end interface

    ! inputs
    real :: X(NZP1,NSCLR)
    integer :: jptr

    ! locals var
    integer :: NS1,k
    real :: rhoCP, Szero, rhof, rhoocn, fk



!     EL = 4187. * (597.31 - 0.56525 * Tc  )  ! set in init cnsts

! ice fluxes ?

    if(Lice) then
    !      treat layers 1,nsice ice formation/melt as a surface turbulent fluxes
    !       nsice set in init_cnsts if <0 use boundary layer

        if(NSICE < 0) then
            NS1 = MAX(1,kmix-1)
        else
            NS1 = NSICE
        endif
        SWFACS = ( 1. - SWDK(-dm(NS1),real(time)) )

        wtI(nz,jptr) = 0.0
        wsI(nz,jptr) = 0.0
        sflux(7,5,jptr) = 0.0
        sflux(8,5,jptr) = 0.0

        do k=NZ,1,-1
        !                           Allow for conservation
            if(LAFLX == 0) then
                rhoCP = rhosw * CPo
                Szero = Sref
                rhof  = rhofrsh
                rhoocn= rhosw
            else
                rhoCP = rho(k)*CP(k)
                Tfrz  = TfrzC(X(k,2)+Sref, -zm(k) )
                Szero = Sref + X(k,2)
                rhof  = rhoh2o
                rhoocn= rho(k)
            endif

        !        fk is the potential amount of ice formation(fk>0) or
        !        melting(fk<0) in layer k  (kgIce / m2/ s)
            fk = rhoCP *hm(k) *(Tfrz-X(k,1)) / (FL * tfscl)


            IF (k <= NS1) THEN   ! shallow (turbulent) frazil
                if(fk >= 0.0) then   ! no non turb contribution to (7,5)
                    sflux(8,5,jptr) = sflux(8,5,jptr) + FL * fk
                    fk = 0.0          ! keeps wxI(k-1) = wxI(k)
                else  ! melt
                    if(fk <= sflux(7,5,jptr)) then
                    !            all deep ice melts, remaining heat can melt sea ice as (8,4)<0
                        sflux(8,5,jptr) =sflux(8,5,jptr) &
                          +FL*(fk-sflux(7,5,jptr))
                        fk = sflux(7,5,jptr)   !zero (7,5), heat &salt with wxI(k-1)
                    endif ! hot ocean (fk<0) melts (fk) deep ice; no melt potential
                endif
            ELSE                      ! deep non-turbulent frazil
                if(sflux(7,5,jptr) >= fk) then
                !                fk < 0 and melts all rising deep ice
                    fk = sflux(7,5,jptr)  !again (7,5)=0 ; heat &salt with wxI(k-1)
                endif
            ENDIF
        !        here fk is addition (or subtraction) to rising deep ice

            sflux(7,5,jptr) = sflux(7,5,jptr) - fk
            wtI(k-1,jptr) = wtI(k,jptr) - fk * FL / rhoCP
            wsI(k-1,jptr) = wsI(k,jptr) &
            - fk * Szero / rhof         &       !  freshwater loss
            + fk * Sfrazil/(1.-Sfrazil/1000.) / rhoocn  ! salt loss

        ENDDO
    !                    exit with sflux(7,5)
    !           set (8,4)<0; (8,5)>0; (7,4)<0; (9,4)<0; (9,5)<0

    !     ocean to ice fluxes  frazil ice at (Sice)
        sflux(8,4,jptr) =AMIN1(0.0, sflux(8,5,jptr))! Potential to melt (W/m2)<0
        sflux(8,5,jptr) =AMAX1(0.0, sflux(8,5,jptr)) ! sfc turb heat (W/m2)>0
        sflux(7,4,jptr) =sflux(7,5,jptr)  & !  kg of fresh deep frazil kg/m2/s <0
        -sflux(8,5,jptr)/FL ! kg of shallow frazil kg/m2/s <0
        sflux(9,4,jptr) =sflux(7,4,jptr) &
        * Sice /(1000.-Sice) ! salt S=Sice
        sflux(9,5,jptr) =                   & !   (9,5) = turbulent salt
        sflux(9,4,jptr)                     & !   total -
        - sflux(7,5,jptr) * Sfrazil / (1000.-Sfrazil) ! non-turbulent

    !     open water thawing rate for fice calculations
        FROCN =    sflux(4,2,jptr)      + sflux(5,2,jptr) &
        + sflux(6,2,jptr) * EL + sflux(3,2,jptr) * SWFACS &
        - sflux(8,2,jptr) * FL

    endif  ! Lice ice ocean frazil fluxes

    return
    end subroutine o2iflx

!**************************************
    subroutine TOPFLX(jptr,kforce)

    use kei_parameters
    use kei_common
    use kei_icecommon
    use kei_hacks, only: ic_conform

    implicit none

    ! inputs
    integer :: jptr
    real :: kforce(forcing_var_cnt)

    ! local vars
    integer :: n
    real :: aice,focn1

!           Once all i=1,4 levels are loaded in jptr=0
!       load top of ocean fluxes (sflux level 5 inputs)
!    NB  Deep Frazil brine flux (sflux(7,5)) accounted for in ntflx
!    Shallow Frazil turbulent brine flux sflux(8,5) is from o2iflux
!       complete the nonturbulent ocean fluxes call ntflx
!     EL = 4187. * (597.31 - 5.6525)
    aice = (1. - focn)          !   for compatable fluxes
    if (ic_conform .eq. 2) then
       if (aice .eq. 0.) then
          aice = kforce(ic_f_ind)
          focn1 = (1.- aice)/focn
      else
        focn1 = 1.
      endif
    else
        focn1 = focn
    endif

    do n=1,3
        sflux(n,5,jptr) = sflux(n,2,jptr)*focn1 + aice*sflux(n,4,jptr)
    END DO

!   total turbulent surface heat flux (W/m2)
    atm_flux_to_ocn_surface = sflux(4,2,jptr)*focn1        & ! ocn long net (down+up)
!    >                + sflux(9,2,jptr)        ! ocn long up
    + sflux(5,2,jptr)*focn1      &  ! sensible heat
    + sflux(6,2,jptr)*EL   &  ! latent
    - sflux(8,2,jptr)*FL  ! melt snow

    sflux(4,5,jptr) = atm_flux_to_ocn_surface &
    + sflux(4,4,jptr) !*aice    ! ice-ocean - aice partitioning now in kei_ice

    atm_flux_to_ocn_surface = atm_flux_to_ocn_surface + sflux(3,2,jptr)  ! add in shortwave

    print *, 'atm_flux_to_ocn_surface: ',atm_flux_to_ocn_surface

    ! removed - already done in sia2 melting routines ...
!    if (sflux(6,4,jptr) > 0.0) then
!        sflux(4,5,jptr) = sflux(4,5,jptr) &
!        +CPo*(Tmlt-Tref)* sflux(6,4,jptr) !*fice   ! change T of melted ice/snow
!    endif

!    total turbulent fresh water forcing (kg Fwater/m**2/s)
    sflux(6,5,jptr) = sflux(6,2,jptr)  &      ! evaporation
    + sflux(7,2,jptr)       & ! rainfall
    + sflux(8,2,jptr)       & ! snowfall
    + sflux(6,4,jptr) !*aice    ! ice-ocean - aice partitioning now in kei_ice

!     total turbulent salt  flux (kgS/m**2/s)
    sflux(5,5,jptr) = &
        + sflux(5,4,jptr) & !*aice - aice partitioning now in kei_ice
        + sal_correction_rate  ! kg/s

!  flux profiles through ocean (non-turbulent)
    call ntflx(jptr)


!    if(iabs(LAFLX) == 4) then
    if(abs(LAFLX) == 4) then
        do  n=1,3
            sflux(n,5,jptr) = sflux(n,5,jptr) + harm0(n)
        ENDDO
        sflux(4,5,jptr) = sflux(4,5,jptr) + harm0(4) + harm0(5)
        sflux(6,5,jptr) = sflux(6,5,jptr) + harm0(6) + harm0(7)
    endif

    RETURN
    end subroutine TOPFLX

!**********************************************************************

    SUBROUTINE ntflx(jptr)

!     compute non-turbulent heat/salt fluxes through ocean:
!     - the vertical profile of non-turb heat flux at interfaces equals
!       sum of solar flux  + deep ice formation.
!     - the vertical profile of non-turb salt flux at interfaces equals
!       brine rejection due to deep ice formation.

    use kei_parameters
    use kei_common

    implicit none

    interface
     FUNCTION SWDK(z,ftime)
     implicit none
     REAL, intent(IN) :: z,ftime
     REAL :: SWDK
     end function SWDK
    end interface

    ! inputs
    integer :: jptr

    ! local vars
    integer :: k
    real :: temp_s,temp_wxnt

!    calculate surface kinematic solar heat flux (C m/s)

    !print *,'shortwave ',sflux(3,5,jptr),'penetration stuff:'
    !print *,'layer decay xwnt'
    do k=0,NZ
        wXNT(k,1) = wtI(k,jptr) &
    !   - sflux(3,5,jptr) * SWDK(-dm(k),real(timed)) &
        - sflux(3,5,jptr) * SWDK(-dm(k),real(time )) &
        / (rho(k)*CP(k))
        temp_s = SWDK(-dm(k),real(time ))
        temp_wxnt = sflux(3,5,jptr) * temp_s / (rho(0)*CP(0))
        !print *,k,temp_s,temp_wxnt,time,real(time)

        wXNT(k,2) = wsI(k,jptr)
    ENDDO

    if(LNBFLX) then   !  no bottom fluxes
        wXNT(nz,1) = 0.0
        wXNT(nz,2) = 0.0
    endif

    return
    end SUBROUTINE ntflx

!*************************************************************

    REAL FUNCTION SWDK(z,ftime)
!     Compute SW exponential decay factor parameters
!     for given Jerlov water types.
!     Process parameter jerlov (in "common.inc") is set in "CAOIM.f":
!     it specifies the water type to be used.  If jerlov=0, local array
!     jerl is used to get different water type for each month.

    use kei_parameters
    use kei_common

    implicit none

    ! Input
    real, intent(IN) :: z, ftime               ! time in days for solar flux

    ! Local
    integer :: mon,j
    real :: swtime
    integer, parameter :: nmax=5
!         types =  I       IA      IB      II      III
!             j =  1       2       3       4       5
    real :: Rfac(nmax) = (/  0.58 ,  0.62 ,  0.67 ,  0.77 ,  0.78 /)
    real :: a1(nmax) = (/  0.35 ,  0.6  ,  1.0  ,  1.5  ,  1.4  /)
    real :: a2(nmax) = (/ 23.0  , 20.0  , 17.0  , 14.0  ,  7.9  /)
!          month  1   2   3   4   5   6   7   8   9   10  11  12
    integer :: jerl(12) = &
      (/ 2 , 2 , 2 , 3 , 3 , 3 , 4 , 4 , 4 , 4 , 3 , 2 /)

    swtime = ftime
    do while(swtime >= dpy)
        swtime = swtime-dpy
    enddo
    mon = int(dble(swtime/dpy)*12.) + 1 !(swtime/30.147) with dpy=365.
    mon = MIN0(mon,12)

    if ( jerlov < 1. ) then
        j = jerl(mon)
    else
        j = jerlov
    endif
!      write(6,*) 'time=',ftime,' mon=',mon,' j=',j

    SWDK =        Rfac(j)  * dexp(dble(z/a1(j))) &
    + (1.0-Rfac(j)) * dexp(dble(z/a2(j)))

    return
    END FUNCTION SWDK


!*************************************************************

    SUBROUTINE init_flx(X,kforce,jptr,timed)

!     Set up parameters for calculating fluxes and initialize fluxes.
!     intermediate values computed every ndtld
!     common deltax(NJDT), xbar, denom

    use kei_parameters
    use kei_common

    implicit none

! Input
    real ::   X(NZP1,NSCLR)
    real :: kforce(forcing_var_cnt)
    integer :: jptr

    DOUBLE PRECISION :: timed
! Local
    integer :: j,itemp
    real :: dt,bnum
    real ::   xv(NJDT)

! extrapolation parameters for fluxes
    if ( NJDT == 1 ) then
        ndtld = ndtflx
        xbar = -0.5 * dtsec * ndtflx
        deltax(1) = 0.0
        denom = 1.0
    else
        itemp = NJDT
        ndtld = ndtflx / (itemp-1)
        xbar = 0.0
        do j=1,NJDT
            xv(j) = -0.5 * dt * (ndtflx + 2 * ndtld * (j-1) )
            xbar = xbar + xv(j)
        ENDDO
        xbar = xbar / NJDT
        denom = 0.0
        bnum  = 0.0
        do j=1,NJDT
            deltax(j) = xv(j) - xbar
            denom = denom + deltax(j)**2
        ENDDO

    endif

!    Initialize flux arrays
    !call aset(  wtI,    nzp1*(njdt+1),0.0)
    !call aset(  wsI,    nzp1*(njdt+1),0.0)
    !call aset(sflux,nsflxs*5*(njdt+1),0.0)
    wtI = 0.0
    wsI = 0.0
    sflux = 0.0

    jptr = 0
!                            level 1,2,3
    do j=1,NJDT
        call atmflx(jptr,timed)
    ENDDO
    call calflx(jptr)

    ! ice - ocean fluxes iff lice
!     if(lice) then
!      ! ocean to ice frazil fluxes
!      call o2iflx(X,0)  ! frazil ocean to ice fluxes
!
!      ! use n=1,3 level  4 from iceflx
!      call iceflx(0)
!      ! zero ice ocean fluxes n=4,6 to complete level 4
!      do 220 n=4,6
!        sflux(n,4,0) = 0.0
!        20 continue
!     endif

!                          load top of ocean fluxes
    call topflx(0,kforce)

    return
    end SUBROUTINE init_flx

!***********************************************************************

    SUBROUTINE FZOL (ZOL,PSIM,PSIG)
!     ESTIMATE THE STABILITY FUNCTIONS
!     USING THE MOST UP TO DATE FLUX-PROFILE RELATIONSHIPS
    implicit none

    ! inputs
    real, intent(in) :: ZOL
    real, intent(out) :: PSIM,PSIG

    ! local
    real :: X

    IF(ZOL < 0.0) then
!           UNSTABLE
      X=(1.0-16.0*ZOL)**0.25
      PSIG=2.0*ALOG(0.5+0.5*X**2)
      PSIM=0.5*PSIG+2.0*ALOG(0.5*(1.0+X))-2.0*ATAN(X)+1.571
    else

!          STABLE
      PSIG=-5.0*ZOL
      PSIM=-5.0*ZOL
    ENDIF
    RETURN
    END SUBROUTINE FZOL

!***********************************************************************

    REAL PURE FUNCTION qsat(mode, TK)
!     calculate saturation humidity in (kg/m3) for given
!     temperature in (K), over either water mode = 0,
!                                or   ice   mode = 1.
    implicit none

    integer, intent(in) :: mode
    real, intent(in) :: TK
    real :: p_h2o

    if(mode == 1) THEN
        qsat  =  11637800./exp(5897.8/TK)
    else
        qsat  =  640380./exp(5107.4/TK)
    endif

    return
    end FUNCTION qsat

! ******************************************************************

    subroutine swokta(rlon,rlat,dpy,timed,cloudf,sw)

!     compute short wave solar radiation flux (W/m2)
!     with cloud correction based on Dobson & Smith (1988),
!     where clear sky values(jokta=0) are computed with coefficients
!     for jokta=1, because D&S clear sky transmissivity is too low.
!     Change: use "dpy" as input (7-1-93)

    implicit none

! input
    real :: rlat,     &  ! latitude                          (radians)
    rlon,              &  ! longitude (+=east)                (radians)
    cloudf,           &   ! cloud fraction from 0 to 1.0
    dpy
    DOUBLE PRECISION :: timed ! time in Julian days, G.M.T.    (days)
    integer :: jokta       ! cloud cover                       (oktas)

! output
    real :: sw             ! sw flux                           (W/m2)

! local
    real :: day,  & !   integer of time
    hour,         & !   fraction of day                   (0. to 1.)
    ror0,         & !   (r/r0)**2 eccentricity factor
    ds,           & !   declination of the earth          (radians)
    h,            & !   loc, hour angle since solar noon  (radians)
    sinphi,       & !   sine of solar elevation (phi)
    trans           !   atm. transmission factor
    integer :: mokta       ! cloud cover = MIN(1,jokta)        (oktas)

!     jokta   =  0    1    2    3    4    5    6    7    8    9
    real :: a(0:9) = (/ .400,.517,.474,.421,.380,.350,.304,.230,.106,.134 /) !   D&S linear regression coeffs:
    real :: b(0:9) = (/ .386,.317,.381,.413,.468,.457,.438,.384,.285,.295 /) !   depend on okta(9=sky obscured)
    real :: s0 = 1368.0 ! solar constant, = 1368            (W/m2)
    real :: pi = 3.141592653
    real :: twopi = 6.283185307

    real :: dss,cosphi,sws,phids,timelat,phid,transcl,swcl

    jokta = int(cloudf*8. + 0.5)
! set time variables

    day  = aint(timed)
    hour = timed-day

! compute sine of solar elevation (sinphi) as a function of time

    ror0   = 1. + 0.033*cos(twopi*day/dpy)      ! 365.25)
    ds     = 0.409 * cos(twopi*(day-173)/dpy)   ! 365.25)
    h      = twopi*hour + rlon
    sinphi = sin(rlat)*sin(ds) - cos(rlat)*cos(ds)* &
    cos(h)

! calculate transmittance (tr) and surface short wave

    if(sinphi >= 0) then
        mokta = max(1,jokta)
        trans = a(mokta) + b(mokta) * sinphi
        sw    = s0 * ror0 * sinphi * trans
    else
        sw    = 0.
    endif

! alternate formula
    if( .FALSE. ) then
        dss    = 23.45 * sin(twopi*(284+day)/dpy  )
        write(6,*) 'timed=',timed,' dec =',ds/twopi*360.,' decs=',dss
        dss    = dss/360.*twopi
        h      = twopi*hour + rlon - pi   ! hour angle since solar noon
        cosphi = sin(rlat)*sin(dss) + cos(rlat)*cos(dss)* &
        cos(h)
        sws    = s0 * ror0 * cosphi * trans
        phids  = 90. - acos(cosphi)/twopi*360.  ! solar elevation
        if(sws < 0.) sws= 0.

        timelat =24.+(h+pi)/pi*12 !time of day(Local Apparent Time=LAT)
        if(timelat >= 24.) timelat=timelat-24.
        phid    = asin(sinphi)/twopi*360. ! solar elevation in degrees
        transcl = a(1) + b(1) * sinphi
        swcl    = s0 * ror0 * sinphi * transcl     ! clear sky sw
        if(swcl < 0.) swcl=0.

        write(6,1000) jokta,phid,phids,hour*24.,timelat,swcl,sw,sws
        1000 format('okta=',i1,' sel=',f5.1,' sels=',f5.1, &
        ' hrGMT=',f4.1,' hrLAT=',f4.1, &
        ' SWcl=',f6.1,' SW=',f6.1,' SWs=',f6.1)
    endif

    return
    end subroutine swokta



