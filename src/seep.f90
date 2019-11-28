program main
  use ieee_arithmetic
  use grid, only: grid_t
  use ssat_env, only: alert, fatal
  use vg, only: vg_t
  use wa_env, only: barom
  use scx
  implicit none

  character(len=255) :: arg, data_file, msg, pwp_file, prec_file, theta_file, usage
  integer :: bins, j, stat
  real :: beta(2), bound, dt, dx, dy, k, i, n, r, p, t, xlim(2), x, y
  type(grid_t) :: pwp, theta
  type(vg_t) :: swc

  usage = 'usage: seep [-bnum] [-dx|yvalue] [-tstep] [-Pprecipitation] [-IA|Bvalue] infile [-o[W]outfile]'
  bins = 100
  dt = 1.
  dx = 1.
  dy = 1.
  xlim(1) = ieee_value(xlim(1), ieee_negative_inf)
  xlim(2) = ieee_value(xlim(2), ieee_positive_inf)
  do j = 1, command_argument_count()
    call get_command_argument(j, arg)
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
      case('o')
        if(arg(3:3) == 'W') then
          theta_file = arg(4:)
        else
          pwp_file = arg(3:)
        end if
      case('t')
        read(arg(3:), *) dt
      case('P')
        prec_file = arg(3:)
      case('I')
        if(arg(3:3) == 'A') then
          read(arg(4:), *) xlim(1)
        else if(arg(3:3) == 'B') then
          read(arg(4:), *) xlim(2)
        else
          call fatal(usage)
        end if
      case default
        call fatal(usage)
      end select
    else
      data_file = arg
    end if
  end do
  if(len_trim(data_file) == 0) call fatal(usage)

  call scxini(data_file, xlim)
  if(len_trim(pwp_file) > 0) call pwp%init(scxdim, dx, dy)
  if(len_trim(theta_file) > 0) call theta%init(scxdim, dx, dy)

  bound = scxdim(1, 1) + scxdim(2, 1)
  x = scxdim(1, 1) + .5 * dx
  do while(x < bound)
    y = scxtop(x) - .5 * dy
    do while(y > scxdim(1, 2))
      call scxwa(x, y, beta, k, i, n, r)
      t = (i - r) / (n - r)
      call theta%set(x, y, t)
      swc%a = beta(1)
      swc%n = beta(2)
      p = barom(y) - t * swc%matsuc(t)
      call pwp%set(x, y, p)
      y = y - dy
    end do
    x = x + dx
  end do
  call scxdel
  call theta%dump(theta_file, msg, stat)
  if(stat /= 0) call alert(msg)
  call pwp%dump(pwp_file, msg, stat)
  if(stat /= 0) call alert(msg)
  stop
end program
