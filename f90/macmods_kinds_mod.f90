module macmods_kinds_mod

  implicit none

  public

  ! ------------------------------------------------------------
  ! Elemental data types/sizes
  ! ------------------------------------------------------------

  integer, parameter :: &
      i4              = 4             , & ! integers are 4 bytes
      log_kind        = kind(.true.)  , & ! logical kind
      r4              = 4             , & ! floats are 4 bytes in KPP
      r8              = 8             , & ! doubles are set to 4 for KPP
      d8              = 8                 ! there are few functions that require 8-byte precision

end module macmods_kinds_mod