module scx
  use air, only: barom
  use bez, only: bezc
  use interp1d, only: interp
  use vgmod, only: vgms
  implicit none
  private

  type stra_t
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

contains
  elemental function scxcrk(x) result(y)
    real, intent(in) :: x
    real :: y

    y = 0
    ! z = 2. * c / w * tan(atan(1.) + .4 * phi)
  end function

  subroutine scxdel
  end subroutine

  subroutine scxnew(name)
    character(*), intent(in) :: name
  end subroutine

  elemental function scxtop(x) result(y)
    real, intent(in) :: x
    real :: y

    y = interp(x, xe, ye)
  end function

  pure subroutine scxslc(n, rel, p, w, c, phi, u, alpha, b, h, stat)
    integer, intent(in) :: n
    real, intent(in) :: rel, p(:, :)
    real, intent(out), dimension(n) :: w, c, phi, u, alpha, b
    real, intent(out) :: h(0:n)
    integer, intent(out) :: stat
    real :: d
    real, dimension(2) :: m, q, r
    integer :: i, j

    stat = 0
    d = rel * norm2(p(1, :) - p(size(p, 1), :))
    q = p(1, :)
    do i = 1, n
      j = i - 1
      r = bezc(real(i)/n, p)
      if(q(1) > r(1)) then
        stat = -i
        return
      end if
      if(scxtop(r(1)) - s

      b(i) = r(1) - q(1)
      h(j) = scxtop(q(1)) - q(2)
      h(i) = scxtop(r(1)) - r(2)
      if(h(i) < 0 .or. h(i) < d) then
        stat = i
        return
      end if

      alpha(i) = asin((q(2) - r(2)) / b(i))
      m = .5 * (q + r)
      call scxvar(m(1), m(2), c(i), phi(i), w(i), u(i))

      q = r
    end do
  end subroutine

  pure subroutine scxvar(x, y, c, phi, w, u)
    real, intent(in) :: x, y
    real, intent(out) :: c, phi, w, u

    c = 0
    phi = 0
    w = 0
    u = 0
    ! w_i = sum(w(:i))
    ! t = (theta - res) / (sat - res)
    ! ua = barom(y)
    ! u = ua - t * vgms(t, a, n)
  end subroutine
end module
