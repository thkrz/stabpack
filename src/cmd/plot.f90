program main
  use opt
  use ieee_arithmetic
  use ssat_env, only: fatal
  use scx
  implicit none

  character(len=255) :: arg, datafile, usage
  character :: c

  usage = 'usage: stab [-ddx] [-ffos] [-n[B]num] [-rratio] [-IA|Bvalue] file'
  xlim(1) = ieee_value(xlim(1), ieee_negative_inf)
  xlim(2) = ieee_value(xlim(2), ieee_positive_inf)
  do while(optget(arg, c))
    select case(c)
    case('d')
      read(arg, *) dx
    case('f')
      read(arg, *) fos
    case('n')
      if(arg(1:1) == 'B') then
        read(arg(2:), *) bezn
      else
        read(arg, *) num
      end if
    case('r')
      read(arg, *) rd
    case('I')
      call flgget(arg, c, .true.)
      select case(c)
      case('A')
        read(arg, *) xlim(1)
      case('B')
        read(arg, *) xlim(2)
      case default
        call fatal(usage)
      end select
    case default
      call fatal(usage)
    end select
  end do
  call flgopt(datafile, err)
  if(err == -1 .or. len_trim(datafile) == 0) call fatal(usage)

  call scxini(datafile, xlim)
  stop
end program
