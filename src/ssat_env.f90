module ssat_env
  use, intrinsic :: iso_fortran_env, only: error_unit
  public

  real, parameter :: g = 9.81, pi = 4. * atan(1.)
  integer, parameter :: ndim = 2

contains
  elemental function degrees(rad) result(deg)
    real, intent(in) :: rad
    real :: deg

    deg = rad * 180. / pi
  end function

  subroutine fatal(s)
    character(*), intent(in) :: s

    write(error_unit, *) trim(s)
    error stop
  end subroutine

  elemental function radians(deg) result(rad)
    real, intent(in) :: deg
    real :: rad

    rad = deg * pi / 180.
  end function
end module
