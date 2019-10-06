module ssat_env
  use, intrinsic :: iso_fortran_env, only: error_unit
  public

  real, parameter :: g = 9.81, pi = 4. * atan(1.)

contains
  elemental function degrees(rad) result(deg)
    real, intent(in) :: rad
    real :: deg

    deg = rad * 180. / pi
  end function

  subroutine fatal(s)
    character(*), intent(in) :: s

    write(error_unit, *) trim(s)
    error stop
  end subroutine

  elemental function radians(deg) result(rad)
    real, intent(in) :: deg
    real :: rad

    rad = deg * pi / 180.
  end function
end module

module slope
  use interp1d, only: interp
  use ssat_env, only: fatal, radians
  implicit none
  private
  public slope_close
  public slope_load
  public slope_stratum
  public slope_top
  public stratum

  type stratum
    real h
    real w
    real phi
    real c
    real theta
    real res
    real sat
    real a
    real n
    real k
  end type

  real, allocatable :: xh(:), yh(:)
  type(stratum), allocatable :: strata(:)
  real :: tana
  logical :: is_init = .false.

contains
  subroutine slope_close
    deallocate(strata)
    deallocate(xh)
    deallocate(yh)
    is_init = .false.
  end subroutine

  subroutine slope_load(name)
    character(*), intent(in) :: name
    real, intent(in) :: alpha
    integer :: id, i, n, stat
    character(255) :: msg

    open(newunit=id, file=name, action='read', status='old', iostat=stat, iomsg=msg)
    if(stat /= 0) call fatal(msg)
    read(id, *, iostat=stat, iomsg=msg) alpha
    if(stat /= 0) call fatal(msg)
    tana = tan(radians(alpha))
    read(id, *, iostat=stat, iomsg=msg) n
    if(stat /= 0) call fatal(msg)
    allocate(xh(n))
    allocate(yh(n))
    read(id, *, iostat=stat, iomsg=msg) ((x(i), y(i)), i=1,n)
    if(stat /= 0) call fatal(msg)

    read(id, *, iostat=stat, iomsg=msg) n
    if(stat /= 0) call fatal(msg)
    allocate(strata(n))
    read(id, *, iostat=stat, iomsg=msg) (strata(i), i=1,n)
    if(stat /= 0) call fatal(msg)
    close(id)
    is_init = .true.
  end subroutine

  pure function slope_stratum(x, y) result(s)
    real, intent(in) :: x, y
    type(stratum), pointer :: s
    integer :: i, n
    real :: x0, y0, y1

    if(.not. is_init) error stop

    n = size(p, 1)
    y1 = slope_top(x)
    do i = 1, n
      associate(h => strata(i)%h)
        x0 = interp(h, yh, xh)
        y0 = slope_top(x0) + tana * (x - x0)
        if(y >= y0 .and. y < y1) then
          s => strata(i)
          return
        end if
      end associate
    end do
  end subroutine

  elemental function slope_top(x) result(y)
    real, intent(in) :: x
    real :: y

    if(.not. is_init) error stop

    y = interp(x, xh, yh)
  end function
end module

module critss
  use bez
  use fmin
  use razdol
  use slope
  implicit none
  private
  public critss_find

contains
  subroutine critss_find(num, n, a, b)
    integer, intent(in) :: num, n
    real, intent(in) :: a, b
    real :: omega(num)

    w = x(n) - x(1)
    dx = w / num
    omega = slope_top(
  end subroutine
end module

module pest
  implicit none
end module

program main
  use, intrinsic :: iso_fortran_env, only: error_unit
  use ssat_env, only: fatal
  use bez
  use critss
  use slope
  implicit none

  character(len=255) :: arg, input, msg
  character(len=4) :: mode
  integer :: id, stat
  real :: alpha

  namelist /CONFIG/ alpha, input, mode

  if(command_argument_count() /= 1) call fatal('control file missing.')
  call get_command_argument(1, arg)
  open(newunit=id, file=trim(arg), status='old', iostat=stat, iomsg=msg, action='read')
  if(stat /= 0) call fatal(msg)
  read(id, nml=CONFIG, iostat=stat, iomsg=msg)
  if(stat /= 0) call fatal(msg)
  close(id)

  call slope_init(input, alpha)
  if(mode == 'crit') call critfind(0, 0, 0., 0.)
  if(mode == 'pest') error stop
  call slope_fin
  stop
contains
  pure subroutine mos(p, slices, stat)
    real, intent(in) :: p(:, :)
    type(slice), intent(out) :: slices(:)
    integer, intent(in) :: stat
    real :: c(size(p, 2)), t(0:size(slices)), y(size(slices)+1)
    integer :: i, j, k, n
    type(stratum) :: s

    k = size(p, 1)
    n = size(slices)
    stat = 0
    do concurrent(i = 0:n)
      t(i) = real(i) / n
    end do

    d = rd * norm2(p(1, :) - p(k, :))
    xy = bezc(t, p)
    if (any(xy(2:, 1) - xy(:n, 1) < 0)) then
      stat = -1
      return
    end if
    y = slope_top(xy(:, 1))
    do i = 1, n
      j = i + 1
      h(1) = y(i) - xy(i, 2)
      h(2) = y(j) - xy(j, 2)
      if (any(h < 0) .or. any(h) > d) then
        stat = -1
        return
      end if
      b = xy(j, 1) - xy(i, 1)
      c = .5 * (xy(j, :) + xy(i, :))
      s = slope_stratum(c(1), c(2))
      ! usw
    end do
  end subroutine
end program
