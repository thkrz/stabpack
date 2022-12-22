module scx
  use ieee_arithmetic, only: ieee_is_finite, ieee_is_nan, ieee_quiet_nan, ieee_value
  use bez, only: bezcrv
  use intp1d, only: interp
  implicit none
  private
  public scxcut
  public scxdel
  public scxdim
  public scxini
  public scxtop

  real, parameter :: rd = 1./3.
  integer, allocatable :: scxtop(:)
  real :: scxdim(2, 2)

contains
  pure subroutine scxcut(n, p, w, c, phi, u, alpha, b, h, stat)
    integer, intent(in) :: n
    real, intent(in) :: p(:, :)
    real, intent(out), dimension(n) :: w, c, phi, u, alpha, b
    real, intent(out) :: h(0:n)
    integer, intent(out) :: stat
    real :: d, z
    real, dimension(2) :: q, r
    integer :: i, j, k

    d = rd * norm2(p(:, 1) - p(:, size(p, 2)))
    h(0) = 0
    q = p(:, 1)
    stat = 0
    do k = 1, n
      r = bezcrv(real(k)/n, p)
      if(q(1) > r(1)) then
        stat = -k
        return
      end if

      b(k) = r(1) - q(1)
      h(k) = scxtop(r(1)) - r(2)
      if(h(k) < 0 .or. h(k) > d) then
        stat = k
        return
      end if

      alpha(k) = asin((q(2) - r(2)) / b(k))
      z = .5 * (h(k-1) + h(k))
      i = floor(r(1))
      j = floor(r(2) - z)
      u(k) = 0 ! get pwp from grid i, j
      c(k) = 0
      phi(k) = 0
      w(k) = 0

      q = r
    end do
  end subroutine

  subroutine scxdel
    deallocate(scxtop)
    scxdim = 0
  end subroutine

  subroutine scxini(x, y, dx, dy, msg, stat)
    real, intent(in) :: x(:), y(:)
    real, intent(in) :: dx, dy
    character(*), intent(out) :: msg
    integer, intent(out) :: stat
    integer :: m, n

    stat = 0
    n = size(x)
    if(n /= size(y)) then
      stat = -1
      msg = 'x and y must be the same size'
      return
    end if
    if any((x(2:) - x(:n-1)) <= 0) then
      stat = -1
      msg = 'x must be monotonically increasing'
      return
    end if

    scxdim(1, 1) = x(1)
    scxdim(1, 2) = minval(y)
    scxdim(2, 1) = x(n)
    scxdim(2, 2) = maxval(y)

    y = y - scxdim(1, 2)
    m = ceiling((x(n) - x(1)) / dx)
    allocate(scxtop(0:m))
    do concurrent i = 0:m
      scxtop(i) = nint(interp(i * dx, x, y) / dy)
    end do
  end subroutine

  elemental function scxtop(x) result(y)
  end function
end module
