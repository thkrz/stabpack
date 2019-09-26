module slope
  use interp1d, only: interp
  use ssat_env, only: fatal
  implicit none
  private
  public slopefin
  public slopeinit
  public slopeparam
  public slopesurf

  integer, parameter :: nparam = 10
  real, allocatable :: elev(:, :), prop(:, :), x0(:)
  real :: tana
  logical :: is_initialized = .false.

contains
  subroutine slopefin
    deallocate(x0)
    deallocate(prop)
    deallocate(elev)
    is_initialized = .false.
  end subroutine

  subroutine slopeinit(name, alpha)
    character(*), intent(in) :: name
    real, intent(in) :: alpha
    integer :: id, i, n, stat
    character(255) :: msg

    tana = tan(alpha)

    open(newunit=id, file=name, action='read', status='old', iostat=stat, iomsg=msg)
    if(stat /= 0) call fatal(msg)
    read(id, *, iostat=stat, iomsg=msg) n
    if(stat /= 0) call fatal(msg)
    allocate(elev(n, ndim))
    read(id, *, iostat=stat, iomsg=msg) elev
    if(stat /= 0) call fatal(msg)

    read(id, *, iostat=stat, iomsg=msg) n
    if(stat /= 0) call fatal(msg)
    allocate(prop(n, nparam))
    read(id, *, iostat=stat, iomsg=msg) prop
    if(stat /= 0) call fatal(msg)
    close(id)

    allocate(x0(n))
    do concurrent(i = 1, n)
      x0(i) = interp(prop(i, 1), elev(:, 2), elev(:, 1))
    end do
    is_initialized = .true.
  end subroutine

  pure subroutine slopeparam(x, y, w, c, phi) ! theta, res, sat, a, n, k
    real, intent(in) :: x, y
    real, intent(out) :: w, c, phi
    integer :: i, n
    real :: y0, y1

    if(.not. is_initialized) error stop

    n = size(prop, 1)
    y1 = slopesurf(x)
    do i = 1, n
      y0 = slopesurf(x0(i)) + tana * (x - x0(i))
      if(y < y1 .and. y >= y0) then
        w = prop(i, 2)
        phi = prop(i, 3)
        c = prop(i, 4)
        return
      end if
    end do
  end subroutine

  elemental function slopesurf(x) result(y)
    real, intent(in) :: x
    real :: y

    if(.not. is_initialized) error stop

    y = interp(x, elev(:, 1), elev(:, 2))
  end function
end module
