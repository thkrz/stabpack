program main
  use num_env, only: swap
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
  deallocate(arr)
  stop

contains
  subroutine debkin(x, y1, y2)
    real, intent(in) :: x(:), y1(:)
    real, intent(out) :: y2(:)
    real :: k(size(y2))
    integer :: i, n

    n = size(x)
    i = minloc(y)
    k = 0
    call kin_(x, y, drag, 1, k)
    call kin_(x, y, drag, -1, k)
    k = k / maxval(k)
    if(k(i) /= 1) call fatal('invalid drag.')
    y2 = y - k * (depth(2) - depth(1)) + depth(1)
  end subroutine

  pure subroutine kin_(x, y, crr, dir, r)
    real, intent(in) :: x(:), y(:), crr
    integer, intent(in) :: dir
    real, intent(inout) :: r(:)
    real :: wa(size(r))
    integer :: e1, e2, e3, i, ii, j, jj, n

    n = size(x)
    e1 = 2
    e2 = n
    e3 = dir
    if(dir < 0) call swap(e1, e2)
    wa = 0
    do i = e1, e2, e3
      j = i - 1
      h = dir * (y(j) - y(i))
      s = sqrt((x(i) + x(j))**2 + h**2)
      alpha = asin(h / s)
      ii = i
      jj = j
      if(dir < 0) call swap(ii, jj)
      wa(ii) = max(0, h + wa(jj) - crr * cos(alpha) * s)
    end do
    r = r + wa
  end subroutine
end program
