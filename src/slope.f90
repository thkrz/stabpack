module slope
  use interp1d
  use ssat_env
  implicit none
  private
  public slope_finalize
  public slope_init
  public slope_parameters
  public slope_surface

  integer, parameter :: nparam = 10
  real, dimension(:, :), allocatable :: elev, prop
  logical :: is_initialized = .false.

contains
  subroutine slope_finalize
    deallocate(prop)
    deallocate(elev)
    is_initialized = .false.
  end subroutine

  subroutine slope_init(name)
    character(*), intent(in) :: name
    integer :: id, n

    open(newunit=id, file=name, action='read', status='old')
    read(id, *) n
    allocate(elev(n, ndim))
    read(id, *) elev

    read(id, *) n
    allocate(prop(n, nparam))
    read(id, *) prop
    close(id)
    is_initialized = .true.
  end subroutine

  pure subroutine slope_parameters(x, y, alpha, gamma, phi, c) ! theta, res, sat, a, n, k
    real, intent(in) :: x, y, alpha
    real, intent(out) :: gamma, phi, c
    integer :: i, n
    real :: x0, y0, y1

    if(.not. is_initialized) error stop

    n = size(prop, 1)
    y1 = slope_surface(x)
    do i = 1, n
      x0 = interp(prop(i, 1), elev(:, 2), elev(:, 1))
      y0 = slope_surface(x0) + tan(alpha) * (x - x0)
      if(y < y1 .and. y >= y0) then
        gamma = prop(i, 2)
        phi = radians(prop(i, 3))
        c = prop(i, 4)
        return
      end if
    end do
  end subroutine

  elemental function slope_surface(x) result(y)
    real, intent(in) :: x
    real :: y

    if(.not. is_initialized) error stop

    y = interp(x, elev(:, 1), elev(:, 2))
  end function
end module
