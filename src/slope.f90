module slope
  use interp1d
  use ssat_env
  implicit none
  private
  public slopefin
  public slopeinit
  public slopeparam
  public slopesurf

  integer, parameter :: nparam = 10
  real, allocatable :: elev(:, :), prop(:, :), x0(:)
  logical :: is_initialized = .false.

contains
  subroutine slopefin
    deallocate(x0)
    deallocate(prop)
    deallocate(elev)
    is_initialized = .false.
  end subroutine

  subroutine slopeinit(name)
    character(*), intent(in) :: name
    integer :: id, i, n

    open(newunit=id, file=name, action='read', status='old')
    read(id, *) n
    allocate(elev(n, ndim))
    read(id, *) elev

    read(id, *) n
    allocate(prop(n, nparam))
    read(id, *) prop
    close(id)

    allocate(x0(n))
    do concurrent(i = 1, n)
      x0(i) = interp(prop(i, 1), elev(:, 2), elev(:, 1))
    end do
    is_initialized = .true.
  end subroutine

  pure subroutine slopeparam(x, y, alpha, params) ! theta, res, sat, a, n, k
    real, intent(in) :: x, y, alpha
    real, intent(out) :: params(nparam)
    integer :: i, n
    real :: tana, x0, y0, y1

    if(.not. is_initialized) error stop

    tana = tan(alpha)
    n = size(prop, 1)
    y1 = slopesurf(x)
    do i = 1, n
      y0 = slopesurf(x0(i)) + tana * (x - x0(i))
      if(y < y1 .and. y >= y0) then
        params = prop(i, 2:)
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
