program main
  use ieee_arithmetic
  use mesh, only: grid_t
  use ssat_env, only: fatal
  use scx
  implicit none

  character(len=255) :: arg, datafile, precfile, usage
  integer :: bins, i, id
  real :: dt, dx, dy, xlim(2)

  usage = 'usage: seep [-bnum] [-dx|yvalue] [-tstep] [-Pprecfile] [-x0|1value] file'
  bins = 100
  dt = 1.
  xlim(1) = ieee_value(xlim(1), ieee_negative_inf)
  xlim(2) = ieee_value(xlim(2), ieee_positive_inf)
  do i = 1, command_argument_count()
    call get_command_argument(i, arg)
    if(arg(:1) == '-') then
      select case(arg(2:2))
      case('b')
        read(arg(3:), *) bins
      case('d')
        if(arg(3:3) == 'x') then
          read(arg(4:), *) dx
        else if(arg(3:3) == 'y') then
          read(arg(4:), *) dy
        else
          call fatal(usage)
        end if
      case('t')
        read(arg(3:), *) dt
      case('P')
        precfile = arg(3:)
      case('x')
        if(arg(3:3) == '0') then
          read(arg(4:), *) xlim(1)
        else if(arg(3:3) == '1') then
          read(arg(4:), *) xlim(2)
        else
          call fatal(usage)
        end if
      case default
        call fatal(usage)
      end select
    else
      datafile = arg
    end if
  end do
  if(len_trim(datafile) == 0) call fatal(usage)

  call scxini(datafile, xlim)
  call scxdel
  stop

contains
end program
