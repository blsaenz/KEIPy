! ======================================================================
! Global Paramaters, Types, & Variables
! ======================================================================
module macmods_util_mod

  use macmods_kinds_mod, only : r8,i4,log_kind

  use macmods_param_mod

  implicit none

  public

  ! Named enumerations
  integer, parameter, public ::                 &
       init_seeding_type        = 0             , &
       monthly_seeding_type     = 1             , &
       seed_on_harvest_type     = 2             , &
       breakage_Duarte_Ferreira = 0             , &   ! fancy breakage
       breakage_Rodrigues       = 1                ! death rate scales with swh

  real(kind=r8), parameter, public ::        &
      c0     =    0.0_r8                   , &
      c1     =    1.0_r8                   , &
      c2     =    2.0_r8                   , &
      c3     =    3.0_r8                   , &
      c4     =    4.0_r8                   , &
      c10    =   10.0_r8                   , &
      c1000  = 1000.0_r8                   , &
      p001   =    0.001_r8                 , &
      p5     =    0.5_r8                   , &
      pi     =    3.14159265358979323846_r8


  !---------------------------------------------------------------------
  !  Unit Conversion
  !---------------------------------------------------------------------

  real(kind=r8), parameter, public :: &
      mw_c      = 12.011_r8,  & ! molecular weight C
      mw_n      = 14.00672_r8,  & ! molecular weight N
      mw_p      = 30.973762_r8, & ! molecular weight P
      mw_Fe     = 55.845_r8,    & ! molecular weight Fe
      sphr      = 3600.0_r8,    & ! number of seconds in an hour
      spd       = 86400.0_r8,   & ! number of seconds in a day
      dpy       = 365.0_r8,     & ! number of days in a year
      spy       = dpy*spd,      & ! number of seconds in a year
      hrps      = c1 / sphr,    & ! number of hours in a second
      dps       = c1 / spd,     & ! number of days in a second
      ypd       = c1 / dpy,     & ! number of years in a day
      yps       = c1 / spy,     & ! number of years in a second
      cmperm    = 100.0_r8,     & ! cm per meter
      mpercm    = .01_r8          ! meters per cm

  !---------------------------------------------------------------------
  !  Physical Constants
  !---------------------------------------------------------------------

  real(kind=r8), parameter, public :: &
      vonkar    =   0.4_r8,            & ! von Karman constant
      T0_Kelvin = 273.15_r8,           & ! freezing T of fresh water (K)
      K_Boltz   =   8.617330350e-5_r8, & ! Boltzmann constant (eV/K)
      rho_sw    =   1.026_r8,          & ! density of salt water (g/cm^3)
      epsC      =   1.0e-8_r8,         & ! small C concentration (mmol C/m^3)
      epsTinv   =   3.17e-8_r8,        & ! small inverse time scale (1/year) (1/sec)
      molw_Fe   =  55.845_r8,          & ! molecular weight of iron (gFe / mol Fe)
      R13C_std  =   1.0_r8,            & ! actual 13C/12C PDB standard ratio (Craig, 1957) = 1123.72e-5_r8
      R14C_std =    1.0_r8               ! actual 14C/12C NOSAMS standard ratio = 11.76e-13_r8



  CONTAINS

!	real(kind=d8) PURE FUNCTION sind(deg)
!		sind = sin(deg*3.141592653589793_d8/180._d8)

!	real(kind=d8) PURE FUNCTION asind(n)
!		asind = asin(n)/3.141592653589793_d8*180._d8

!	real(kind=d8) PURE FUNCTION cosd(deg)
!		cosd = cos(deg*3.141592653589793_d8/180._d8)

!	real(kind=d8) PURE FUNCTION acosd(n)
!		acosd = acos(n)/3.141592653589793_d8*180._d8

    real(kind=r8) PURE FUNCTION daylength_r8(lat,lon,alt,dte)
!    real(kind=r8) FUNCTION daylength(lat,lon,alt,dte)

      use macmods_kinds_mod, only : r8
      implicit none

      real(kind=r8), INTENT(IN) :: lat,lon,alt,dte
      real(kind=r8) :: n2000, Js, M, C, lambda, delta, h, omega ! dte_2000_days, Jt

      ! main function that computes daylength and noon time
      ! https://en.wikipedia.org/wiki/Sunrise_equation

      ! number of days since Jan 1st, 2000 12:00 UT
      !dte_2000_days = 730490.0_r8
      !n2000 = dte - dte_2000_days + 68.184_r8/86400._r8
      n2000 = dte - 0.5 + 68.184_r8/86400._r8

      ! mean solar moon
      Js = n2000 - lon/360._r8

      ! solar mean anomaly
      M = mod(357.5291_r8 + 0.98560028_r8*Js,360._r8)

      ! center
      C = 1.9148_r8*sind(M) + 0.0200_r8*sind(2.0_r8*M) + 0.0003_r8*sind(3.0_r8*M)

      ! ecliptic longitude
      lambda = mod(M + C + 180._r8 + 102.9372_r8,360._r8)

      ! solar transit -- don't need for only daylength
      !Jt = 2451545.5_r8 + Js + 0.0053_r8*sind(M) - 0.0069_r8*sind(2.0_r8*lambda)

      ! Sun declination
      delta = asind(sind(lambda)*sind(23.44_r8))

      ! hour angle (day expressed in geometric degrees)
      h = (sind(-0.83_r8 - 2.076_r8*sqrt(alt)/60._r8) - sind(lat)*sind(delta))/(cosd(lat)*cosd(delta))


      ! to avoid meaningless complex angles: forces omega to 0 or 12h
      if (h < -1) then
        omega = 180._r8
      elseif (h > 1) then
        omega = 0._r8
      else
        omega = acosd(h)
      endif

      !print *,dte_2000_days,n2000,Js,M,C,lambda,delta,h,omega,lat,lon,dte

      daylength_r8 = omega/180._r8

    end FUNCTION daylength_r8

    real(kind=d8) PURE FUNCTION daylength(lat,lon,alt,dte)
!    real(kind=r8) FUNCTION daylength(lat,lon,alt,dte)

      use macmods_kinds_mod, only : r8,d8
      implicit none

      real(kind=r8), INTENT(IN) :: lat,lon,alt,dte
      real(kind=d8) :: n2000, Js, M, C, lambda, delta, h, omega ! dte_2000_days, Jt

      ! main function that computes daylength and noon time
      ! https://en.wikipedia.org/wiki/Sunrise_equation

      ! number of days since Jan 1st, 2000 12:00 UT
      !dte_2000_days = 730490.0_d8
      !n2000 = dte - dte_2000_days + 68.184_d8/86400._d8
      n2000 = dte - 0.5 + 68.184_d8/86400._d8

      ! mean solar moon
      Js = n2000 - lon/360._d8

      ! solar mean anomaly
      M = mod(357.5291_d8 + 0.98560028_d8*Js,360._d8)

      ! center
      C = 1.9148_d8*sind(M) + 0.0200_d8*sind(2.0_d8*M) + 0.0003_d8*sind(3.0_d8*M)

      ! ecliptic longitude
      lambda = mod(M + C + 180._d8 + 102.9372_d8,360._d8)

      ! solar transit -- don't need for only daylength
      !Jt = 2451545.5_d8 + Js + 0.0053_d8*sind(M) - 0.0069_d8*sind(2.0_d8*lambda)

      ! Sun declination
      delta = asind(sind(lambda)*sind(23.44_d8))

      ! hour angle (day expressed in geometric degrees)
      h = (sind(-0.83_d8 - 2.076_d8*sqrt(alt)/60._d8) - sind(lat)*sind(delta))/(cosd(lat)*cosd(delta))


      ! to avoid meaningless complex angles: forces omega to 0 or 12h
      if (h < -1) then
        omega = 180._d8
      elseif (h > 1) then
        omega = 0._d8
      else
        omega = acosd(h)
      endif

      !print *,dte_2000_days,n2000,Js,M,C,lambda,delta,h,omega,lat,lon,dte

      daylength = omega/180._d8

    end FUNCTION daylength


    real(kind=r8) PURE FUNCTION lambda_NO3(magu,Tw,CD,VmaxNO3,KsNO3,NO3)

      use macmods_kinds_mod, only : r8, i4
      implicit none

      real(kind=r8), INTENT(IN) :: magu,Tw,CD,VmaxNO3,KsNO3,NO3

      integer(kind=i4), parameter :: n_length = 25
      real(kind=r8), parameter :: pi = 3.14159265358979323846_r8
      real(kind=r8), parameter :: visc = 1.e-6_r8 * 86400._r8
      real(kind=r8), parameter :: Dm = (18_r8*3.65e-11_r8 + 9.72e-10_r8) * 86400._r8
      real(kind=r8), dimension(n_length) :: vval  ! two v's b/c val seems to be a keyword in fortran

      integer(kind=i4) :: n
      real(kind=r8) :: DBL,oscillatory,flow,beta, magu_m, Tw_s,NO3_u

      ! unit conversions, minima
      magu_m = max(1.0_r8,magu*86400.0_r8) ! m/day
      Tw_s = max(0.01_r8,Tw)/86400.0_r8  ! whaa...?
      NO3_u = NO3*1000.0_r8 ! converting from uM to umol/m3


      DBL = 10._r8 * (visc / (sqrt(CD) * abs(magu_m)))

      ! 1. Oscillatory Flow

      do n = 1,n_length
        vval(n) = (1-exp((-Dm * n**2 * pi**2 *Tw_s)/(2._r8*DBL**2)))/(n**2 * pi**2)
      enddo
      oscillatory = ((4._r8*DBL)/Tw_s) * sum(vval)

      ! 2. Uni-directional Flow
      flow = Dm / DBL

      beta = flow + oscillatory

      lambda_NO3 = c1 + (VmaxNO3 / (beta*KsNO3)) - (NO3_u/KsNO3);

    end FUNCTION lambda_NO3


    real(kind=r8) PURE FUNCTION temp_lim(sst,Topt1,K1,Topt2,K2)

      use macmods_kinds_mod, only : r8
      implicit none

      real(kind=r8), INTENT(IN) :: sst,Topt1,K1,Topt2,K2

      if (sst >= Topt1) then
        if (sst <= Topt2) then
          temp_lim = 1._r8
        else
          temp_lim = exp(-K2*(sst-Topt2)**2)
        endif
      else
        temp_lim = exp(-K1*(sst-Topt1)**2)
      endif
      temp_lim = max(0._r8,temp_lim)
      temp_lim = min(1._r8,temp_lim)

    end FUNCTION temp_lim


    logical(kind=log_kind) PURE FUNCTION valid_float(float_to_test)


      use macmods_kinds_mod, only : r8

      implicit none
      real(kind=r8), INTENT(IN) :: float_to_test
      real(kind=r8) :: infinity

      valid_float = .TRUE.

      infinity = HUGE(float_to_test)
      if(float_to_test > infinity) then
        valid_float = .FALSE.
      endif

      infinity = TINY(float_to_test)
      if(float_to_test < infinity) then
        valid_float = .FALSE.
      endif

      if (float_to_test /= float_to_test) then
        valid_float = .FALSE.
      endif

    end FUNCTION valid_float


    pure function round_to_int_if_close(x, tolerance) result(rounded_value)
      real(kind=r8), intent(in) :: x
      real(kind=r8), intent(in) :: tolerance
      real(kind=r8) :: fractional_part
      integer(kind=i4) :: rounded_value

      fractional_part = abs(x - floor(x + 0.5_r8)) ! Fractional part around nearest integer

      if (fractional_part <= tolerance) then
        rounded_value = nint(x) ! Round to nearest integer
      else
        rounded_value = floor(x + 0.5_r8) ! Effectively round to nearest integer for consistency
      endif

    end function round_to_int_if_close


    pure function round_to_real_if_close(x, tolerance) result(rounded_value_real)
      real(kind=r8), intent(in) :: x
      real(kind=r8), intent(in) :: tolerance
      real(kind=r8) :: fractional_part
      real(kind=r8) :: rounded_value_real

      fractional_part = abs(x - floor(x + 0.5_r8))

      if (fractional_part <= tolerance) then
        rounded_value_real = real(nint(x), kind(x)) ! Round to nearest integer and convert to real
      else
        rounded_value_real = floor(x + 0.5_r8) ! Effectively round to nearest integer as a real
      endif

    end function round_to_real_if_close


    subroutine round_to_real_if_close_array2(x, tolerance)
      real(kind=r8), dimension(:,:), intent(inout) :: x
      real(kind=r8), intent(in) :: tolerance

      WHERE (abs(x - floor(x + 0.5_r8)) <= tolerance)
        x = floor(x + 0.5_r8)
      END WHERE

    end subroutine round_to_real_if_close_array2


    subroutine round_to_real_if_close_array1(x, tolerance)
      real(kind=r8), dimension(:), intent(inout) :: x
      real(kind=r8), intent(in) :: tolerance

      WHERE (abs(x - floor(x + 0.5_r8)) <= tolerance)
        x = floor(x + 0.5_r8)
      END WHERE

    end subroutine round_to_real_if_close_array1

    pure function get_month_from_doy(day_of_year) result(month)
      real(kind=r8), intent(in) :: day_of_year
      integer(kind=i4) :: month

      integer(kind=i4), parameter :: days_in_month(12) = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
      real(kind=r8) :: cumulative_days

      cumulative_days = 0.0_r8
      if (day_of_year < 1.0 .or. day_of_year > 366.0) then
        month = 0 ! Indicate an invalid day of year
        return
      end if

      do month = 1, 12
        cumulative_days = cumulative_days + real(days_in_month(month))
        if (day_of_year <= cumulative_days) then
          return
        end if
      end do

      ! Should not reach here for a valid day of year
      month = 0

    end function get_month_from_doy


    pure function get_dom_from_doy(day_of_year, year) result(day_of_month)
      real(kind=r8), intent(in) :: day_of_year
      integer(kind=i4), intent(in) :: year
      integer(kind=i4) :: day_of_month
      logical(kind=log_kind) :: is_leap

      integer(kind=i4), parameter :: days_in_month_non_leap(12) = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
      integer(kind=i4), parameter :: days_in_month_leap(12) = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
      real(kind=r8) :: cumulative_days
      integer(kind=i4) :: month
      integer(kind=i4), dimension(12) :: current_days_in_month

      is_leap = (mod(year, 4) == 0 .and. mod(year, 100) /= 0) .or. (mod(year, 400) == 0)

      cumulative_days = 0.0_r8
      if (is_leap) then
        current_days_in_month = days_in_month_leap
        if (day_of_year < 1.0 .or. day_of_year > 366.0) then
          day_of_month = 0 ! Indicate an invalid day of year
          return
        end if
      else
        current_days_in_month = days_in_month_non_leap
        if (day_of_year < 1.0 .or. day_of_year > 365.0) then
          day_of_month = 0 ! Indicate an invalid day of year
          return
        end if
      end if

      do month = 1, 12
        if (day_of_year <= cumulative_days + real(current_days_in_month(month))) then
          day_of_month = int(day_of_year - cumulative_days)
          return
        end if
        cumulative_days = cumulative_days + real(current_days_in_month(month))
      end do

      day_of_month = 0 ! Should not reach here for a valid day of year

    end function get_dom_from_doy

end module macmods_util_mod