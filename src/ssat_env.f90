module ssat_env
  use, intrinsic :: iso_fortran_env, only: error_unit

  real, parameter :: pi = 4. * atan(1.)
  integer, parameter :: ndim = 2

contains
  subroutine fatal(s)
    character(*), intent(in) :: s

    write(error_unit, *) trim(s)
    error stop
  end subroutine
end module
