module slope
  use interp1d, only: interp
  use ssat_env, only: fatal, radians
  implicit none
  private
  public slopefin
  public slopeinit
  public slopeparam
  public sloperidge

  real, allocatable :: h(:, :), p(:, :), x0(:)
  real :: tana
  logical :: is_initialized = .false.

contains
  subroutine slopefin
    deallocate(x0)
    deallocate(p)
    deallocate(h)
    is_initialized = .false.
  end subroutine

  subroutine slopeinit(name, alpha)
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

    allocate(x0(n))
    do concurrent(i=1:n)
      x0(i) = interp(p(i, 1), h(:, 2), h(:, 1))
    end do
    is_initialized = .true.
  end subroutine

  pure subroutine slopeparam(x, y, w, phi, c) ! theta, a, n, k
    real, intent(in) :: x, y
    real, intent(out) :: w, phi, c
    integer :: i, n
    real :: y0, y1

    if(.not. is_initialized) error stop

    n = size(p, 1)
    y1 = sloperidge(x)
    do i = 1, n
      y0 = sloperidge(x0(i)) + tana * (x - x0(i))
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

  elemental function sloperidge(x) result(y)
    real, intent(in) :: x
    real :: y

    if(.not. is_initialized) error stop

    y = interp(x, h(:, 1), h(:, 2))
  end function
end module
