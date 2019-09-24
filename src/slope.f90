module slope
  implicit none
  private
  public slope_finalize
  public slope_init
  public slope_parameters
  public slope_surface

  real, parameter :: pi = 4. * atan(1.)
  integer, parameter :: ndim = 2, nparam = 9
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

  pure subroutine slope_parameters(a, alpha, gamma, phi, c) ! a, n, res, sat, theta
    real, intent(in) :: a(ndim), alpha
    real, intent(out) :: gamma, phi, c
    integer :: i, n
    real :: y0(size(prop, 1)), y1

    n = size(prop, 1)
    y1 = slope_surface(a(1))
    if(alpha == 0) then
      y0 = prop(:, 1)
    end if
    do i = 1, n
      if(a(2) < y1 .and. a(2) >= y0(i)) then
        gamma = prop(i, 2)
        phi = prop(i, 3) * pi / 180.
        c = prop(i, 4)
        return
      end if
    end do
  end subroutine

  elemental function slope_surface(x) result(y)
    real, intent(in) :: x
    real :: y
    integer :: i, j, n

    if(x < elev(1, 1)) then
      y = elev(1, 2)
      return
    end if
    n = size(elev, 1)
    do i = 1, n - 1
      j = i + 1
      if(x >= elev(i, 1) .and. x < elev(j, 1)) then
        y = elev(i, 2) + (x - elev(i, 1)) * (elev(j, 2) - elev(i, 2)) / (elev(j, 1) - elev(i, 1))
        return
      end if
    end do
    y = elev(n, 2)
  end function
end module
