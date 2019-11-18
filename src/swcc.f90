program main
  use ssat_env, only: fatal, pi
  use rvcont
  implicit none

  character(len=255) :: arg, datafile, gsdfile, usage
  integer :: i, id, num
  real :: e, n, rho_p, t_res

  usage = 'usage: swcc [-nporosity] [-pdensity] -Ggrainsize [-snum file]'
  num = 0
  rho_p = 2.65
  do i = 1, command_argument_count()
    call get_command_argument(i, arg)
    if(arg(:1) == '-') then
      select case(arg(2:2))
      case('n')
        read(arg(3:), *) n
      case('p')
        read(arg(3:), *) rho_p
      case('G')
        gsdfile = arg(3:)
      case('s')
        read(arg(3:), *) num
      case default
        call fatal(usage)
      end select
    else
      datafile = arg
    end if
  end do
  if(len_trim(datafile) /= 0 .and. num == 0) call fatal(usage)
  e = n / (1. - n)

  stop

contains
  elemental function psi(r)
    real, intent(in) :: r
    real, parameter :: alpha = 1.38, cc = 130.
    real :: ptf

    ptf = cc / r * sqrt(3. / (2. * e) * (1. &
        / (2. * pi * rho_p) * mass(r * 1.0e-3) / r**2)**(alpha - 1.)
  end function

  elemental function theta(r)
    real, intent(in) :: r
    real :: theta

    theta = n * mascum(r * 1.0e-3)
  end function
end program
