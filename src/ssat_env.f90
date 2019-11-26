module ssat_env
  use, intrinsic :: iso_fortran_env, only: error_unit
  public

  real, parameter :: g = 9.81

contains
  subroutine alert(s)
    character(*), intent(in) :: s

    write(error_unit, *) trim(s)
  end subroutine

  subroutine fatal(s)
    character(*), intent(in) :: s

    write(error_unit, *) trim(s)
    stop
  end subroutine
end module
