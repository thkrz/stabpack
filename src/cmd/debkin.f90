program main
  use flag, only: flgget
  use num_env, only: swap
  use ssat_env, only: alert, fatal
  implicit none

  character(len=255) :: input, msg, param, usage
  real :: drag, depth(2)
  real, allocatable :: arr(:, :)
  integer :: err, i, id, m, n
  character :: c

  usage = 'usage: debkin [-ddrag] -zT|Bdepth -Pprofile'
  drag = 0
  depth = 0
  do while(flgget(param, c))
    select case(c)
    case('d')
      read(param, *) drag
    case('z')
      if(flgget(param, c, .true.)) call fatal(usage)
      select case(c)
      case('T')
        read(param, *) depth(2)
      case('B')
        read(param, *) depth(1)
      case default
        call fatal(usage)
      end select
    case('P')
      input = param
    case default
      call fatal(usage)
    end select
  end do
  if(len_trim(input) == 0 .or. depth(1) <= depth(2)) call fatal(usage)

  open(newunit=id, file=input, status='old', iostat=err, iomsg=msg)
  if(err /= 0) call fatal(msg)
  read(id, *, iostat=err, iomsg=msg) m, n
  if(err /= 0) call fatal(msg)
  n = n + 1
  allocate(arr(2+n, m))
  read(id, *, iostat=err, iomsg=msg) (arr(:1+n, i), i = 1, m)
  if(err /= 0) call fatal(msg)
  close(id)
  call debkin(arr(1, :), arr(1+n, :), arr(2+n, :))
  print '(I5,1X,I5)', m, n
  do i = 1, m
    print '(*(F7.2))', arr(:, i)
  end do
  deallocate(arr)
  stop

contains
end program
