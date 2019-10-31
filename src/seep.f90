program main
  use mesh, only: mesh_t
  use ssat_env, only: fatal
  use scx
  implicit none

  character(len=255) :: arg, datafile, precfile
  integer :: bins, i, id
  real :: dt, xlim(2)

  bins = 100
  d = 1.
  xlim = 0
  do i = 1, get_argument_count()
    call get_command_argument(i, arg)
    if(arg(:1) == '-') then
      select case(arg(2:2))
      case('b')
        read(arg(3:), *) bins
      case('d')
        read(arg(3:), *) dt
      case('p')
        precfile = arg(3:)
      case('x')
        if(arg(3:3) == '0') then
          read(arg(4:), *) xlim(1)
        else if(arg(3:3) == '1') then
          read(arg(4:), *) xlim(2)
        end if
      end select
    else
      datafile = arg
    end if
  end do
  if(len_trim(datafile) == 0) call fatal('file missing.')

  call scxini(datafile, xlim)
  call scxdel
  stop

contains
end program
