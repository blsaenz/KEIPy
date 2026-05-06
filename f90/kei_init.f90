
    SUBROUTINE init_tm
!     ==================
!     Modified  3 Mar  1991 - wgl
!              18 Sep  1991 - wgl , coupled mods
!              17 Mar  1992 - jan , implicit scheme

!     set parameters for the time integration procedures

      use kei_parameters ! include 'parameter.inc'
      use kei_common !include 'common.inc'
      use kei_icecommon ! include 'icecommon.inc'

      implicit none

! atmosphere

    dta    = dtsec * ndtatm

! ice

    dti    = dtsec * ndtice
    icemax =  8
    tfscl  =  1.0 * dtsec*MAX0(ndtatm,ndtocn,ndtice)
    inew   = 1
    iold   = 0

! ocean
!     time increment for ocean model (in seconds)
    dto    = dtsec * ndtocn

! set flux time step

    ndtflx = MAX0(ndtatm,ndtocn,ndtice)
    if((((ndtflx/ndtatm)*ndtatm) /= ndtflx) .or. &
    (((ndtflx/ndtocn)*ndtocn) /= ndtflx) .or. &
    (((ndtflx/ndtice)*ndtice) /= ndtflx)) THEN

      100 write(6,*) 'STOP in init tm (input.f):'
      write(6,*) '     time steps do not factor, ndtatm=',ndtatm, &
      ' ndtocn=',ndtocn,' ndtice=',ndtice,' ndtflx=',ndtflx
      stop 93

    ELSE

    ! set time parameters for main loop. Adjust nend so that run ends with
    ! storage call, so that a possible restart run can begin with the next
    ! storage call after inct timesteps, and the plotting can be continuous.

      !     nend   = (nend / ndtflx) * ndtflx
      !nend   = (nend / inct  ) * inct
      finalt = startt + dtsec / spd * nend
      time   = startt
      ntime  = nstart

    endif

    RETURN
  end SUBROUTINE init_tm

! ****************************************************************

SUBROUTINE init_env(U,X)
!     ===================
!     Modified  3 March 1991   -  WGL

!     the physical environment for the model:
!     vertical grid and geographic location

  use kei_parameters ! include 'parameter.inc'
  use kei_common !include 'common.inc'
  use kei_icecommon ! include 'icecommon.inc'
  use kei_ecocommon

  implicit none

  ! function arguments
  real(r4) :: &
    U(NZP1,NVEL), &
    X(NZP1,NSCLR)

  ! local
  integer(i4) :: i

  zmp = abs(zm)

  ! load initial ecosystem profile
  !if (leco)
  !  do i=1,ecosys_tracer_cnt
  !    call kei_init_real(ncf_file,trim(eco_tracer_name(i)),X(:,i+2))
  !  enddo
  !endif

  ! Calculate and remove reference salinity

  Sref=(X(1,2)+X(nzp1,2))/2.
  write(6,1010) Sref
  1010 format(/'salinity reference value =',f11.6/)
  do i=1,nzp1
      X(i,2)=X(i,2)-Sref
  ENDDO

  !    Initial surface temp set SSTdat
  SSTdat = X(1,1)
  Tref   = X(1,1)


  ! set (maximum) dimensions for temporary grids
  if(LTGRID) then
      NZt   = NZ + (NZL+NZU) * (NZDIV-1)
      NZP1t = NZt+1
  else
      NZt   = 1
      NZP1t = 1
  endif

  ! set fraction of land, and initial open water fraction
  if( .NOT. lice) fice = 0.0
  flnd = 0.
  focn = 1. - fice - flnd
  write(6,*) 'input fice' , lice, focn, fice, flnd

  ! compute geographic location
  rlat=dlat*twopi/360.
  rlon=dlon*twopi/360.

  !            Coriolis Parameter
  if(abs(dlat) < 0.5) then
      f = 2. * (twopi/86164.) * sin( 0.5*twopi/360.)
  else
      f = 2. * (twopi/86164.) * sin(dlat*twopi/360.)
  endif
  fCoriolis = f

  return
end SUBROUTINE init_env

