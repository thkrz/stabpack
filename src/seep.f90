program main
  use ieee_arithmetic
  use grid, only: grid_t
  use ssat_env, only: alert, fatal
  use vg, only: vg_t
  use wa_env, only: barom, hystp
  use scx
  implicit none

  character(len=255) :: arg, data_file, msg, pwp_file, prec_file, wa_file, usage
  integer :: bins, j, stat
  real :: beta(2), bound, dt, dx, dy, h, i, k, n, r, p, t, xlim(2), x, y
  type(grid_t) :: pwp, wa
  type(vg_t) :: swc

  usage = 'usage: seep [-bnum] [-dx|yvalue] [-tstep] [-Pprecipitation] [-IA|Bvalue] infile'
  bins = 100
  dt = 1.
  dx = 1.
  dy = 1.
  pwp_file = 'pwp.grd'
  wa_file = 'wa.grd'
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
  call pwp%init(scxdim, dx, dy)
  call wa%init(scxdim, dx, dy)

  if(len_trim(prec_file) == 0) then
    call waini
  else
    call wasim
  end if

  bound = scxdim(1, 1) + scxdim(2, 1)
  x = scxdim(1, 1) + .5 * dx
  do while(x < bound)
    y = scxtop(x) - .5 * dy
    h = 0.
    do while(y > scxdim(1, 2))
      call scxwa(x, y, beta, k, i, n, r)
      t = (i - r) / (n - r)
      call wa%set(x, y, t)
      if(t == 1) then
        h = h + dy
        p = hystp(h)
      else
        h = 0.
        swc%a = beta(1)
        swc%n = beta(2)
        p = barom(y) - t * swc%matsuc(t)
      end if
      call pwp%set(x, y, p)
      y = y - dy
    end do
    x = x + dx
  end do
  call scxdel
  call wa%dump(wa_file, msg, stat)
  if(stat /= 0) call alert(msg)
  call pwp%dump(pwp_file, msg, stat)
  if(stat /= 0) call alert(msg)
  stop

contains
  subroutine waini
    real :: beta(2), bound, i, k, n, r, x, y

    bound = scxdim(1, 1) + scxdim(2, 1)
    x = scxdim(1, 1) + .5 * dx
    do while(x < bound)
      y = scxtop(x) - .5 * dy
      do while(y > scxdim(1, 2))
        call scxwa(x, y, beta, k, i, n, r)
        t = (i - r) / (n - r)
        call wa%set(x, y, t)
        y = y - dy
      end do
      x = x + dx
    end do
  end subroutine
end program
