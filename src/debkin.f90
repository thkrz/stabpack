program main
  use ssat_env, only: fatal
  implicit none

  character(len=255) :: arg, input, msg
  real :: drag, depth(2)
  real, allocatable :: arr(:, :)
  integer :: err, i, id, m, n

  usage = 'usage: debkin [-cdrag] -d0|1depth -Pprofile'
  drag = .04
  depth = 0
  do i = 1, command_argument_count()
    call get_command_argument(i, arg)
    if(arg(:1) == '-') then
      select case(arg(2:2))
      case('c')
        read(arg(3:), *) drag
      case('d')
        if(arg(3:3) == '0') then
          read(arg(4:), *) depth(1)
        else if(arg(3:3) == '1') then
          read(arg(4:), *) depth(2)
        else
          call fatal(usage)
        end if
      case('P')
        input = arg(3:)
      case default
        call fatal(usage)
      end select
    else
      call fatal(usage)
    end if
  end do
  if(len_trim(input) == 0 .or. depth(2) <= depth(1)) call fatal(usage)

  open(newunit=id, file=input, status='old', iostat=err, iomsg=msg)
  if(err /= 0) call fatal(msg)
  read(id, *, iostat=err, iomsg=msg) m, n
  if(err /= 0) call fatal(msg)
  n = n + 1
  allocate(arr(2+n, m))
  read(id, *, iostat=err, iomsg=msg) (arr(:, i), i = 1, m)
  if(err /= 0) call fatal(msg)
  close(id)
  call debkin(arr(1, :), arr(1+n, :), arr(2+n, :))
  print '(I5,X,I5)', m, n
  print '(*(F6.2))', (arr(:, i), i = 1, m)
  stop

contains
  pure subroutine debkin(x, y, d)
    real, intent(in) :: x(:), y(:)
    real, intent(out) :: d(:)
    integer :: i, n

    n = size(x)
  end subroutine
end program
