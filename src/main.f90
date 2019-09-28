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
  public slope_fin
  public slope_init
  public slope_prop
  public slope_top

  real, allocatable :: h(:, :), p(:, :)
  real :: tana
  logical :: is_init = .false.

contains
  subroutine slope_fin
    deallocate(p)
    deallocate(h)
    is_init = .false.
  end subroutine

  subroutine slope_init(name, alpha)
    character(*), intent(in) :: name
    real, intent(in) :: alpha
    integer :: id, i, n, stat
    character(255) :: msg

    tana = tan(radians(alpha))

    open(newunit=id, file=name, action='read', status='old', iostat=stat, iomsg=msg)
    if(stat /= 0) call fatal(msg)
    read(id, *, iostat=stat, iomsg=msg) n
    if(stat /= 0) call fatal(msg)
    allocate(h(n, 2))
    read(id, *, iostat=stat, iomsg=msg) (h(i, :), i=1,n)
    if(stat /= 0) call fatal(msg)

    read(id, *, iostat=stat, iomsg=msg) n
    if(stat /= 0) call fatal(msg)
    allocate(p(n, 10))
    read(id, *, iostat=stat, iomsg=msg) (p(i, :), i=1,n)
    if(stat /= 0) call fatal(msg)
    close(id)
    is_init = .true.
  end subroutine

  pure subroutine slope_prop(x, y, w, phi, c) ! theta, a, n, k
    real, intent(in) :: x, y
    real, intent(out) :: w, phi, c
    integer :: i, n
    real :: x0, y0, y1

    if(.not. is_init) error stop

    n = size(p, 1)
    y1 = slope_top(x)
    do i = 1, n
      x0 = interp(p(i, 1), h(:, 2), h(:, 1))
      y0 = slope_top(x0) + tana * (x - x0)
      if(y >= y0 .and. y < y1) then
        w = p(i, 2)
        phi = radians(p(i, 3))
        c = p(i, 4)
        ! theta = (p(i, 5)- p(i, 6)) / (p(i, 7) - p(i, 6))
        ! a = p(i, 8)
        ! n = p(i, 9)
        ! k = p(i, 10)
        return
      end if
    end do
  end subroutine

  elemental function slope_top(x) result(y)
    real, intent(in) :: x
    real :: y

    if(.not. is_init) error stop

    y = interp(x, h(:, 1), h(:, 2))
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
  end subroutine
end module

module pest
  implicit none
end module

program main
  use, intrinsic :: iso_fortran_env, only: error_unit
  use ssat_env, only: fatal
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

  call slopeinit(input, alpha)
  if(mode == 'crit') call critfind(0, 0, 0., 0.)
  if(mode == 'pest') error stop
  call slopefin
  stop
end program
