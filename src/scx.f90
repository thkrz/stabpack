module scx
  use air, only: barom
  use bez, only: bezc
  use grndwt, only: hsp, piezom
  use interp1d, only: interp
  use vgmod, only: vgms
  implicit none
  private
  public scxcrk
  public scxdel
  public scxnew
  public scxslc
  public scxtop
  public scxvar

  type stra_t
    real xoc
    real yoc
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

  real :: tana

contains
  elemental function scxcrk(x) result(y)
    real, intent(in) :: x
    real :: y

    y = yoc + tana * (x - xoc)
    ! z = 2. * c / w * tan(atan(1.) + .4 * phi)
  end function

  subroutine scxdel
  end subroutine

  subroutine scxnew(name)
    character(*), intent(in) :: name
  end subroutine

  pure subroutine scxslc(n, rd, p, w, c, phi, u, alpha, b, h, stat)
    integer, intent(in) :: n
    real, intent(in) :: rd, p(:, :)
    real, intent(out), dimension(n) :: w, c, phi, u, alpha, b
    real, intent(out) :: h(0:n)
    integer, intent(out) :: stat
    real :: d
    real, dimension(2) :: m, q, r
    integer :: i, j

    stat = 0
    d = rd * norm2(p(1, :) - p(size(p, 1), :))
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

  elemental function scxtop(x) result(y)
    real, intent(in) :: x
    real :: y

    y = interp(x, xe, ye)
  end function

  pure subroutine scxvar(x, y, c, phi, w, u)
    real, intent(in) :: x, y
    real, intent(out) :: c, phi, w, u
    real :: ua, uw

    c = 0
    phi = 0
    w = 0
    u = 0
    ! w = sum(w(:i))
    ! t = min(max((theta - res) / (sat - res), 0.), 1.)
    ! ua = barom(y)
    ! uw = 0
    ! if(t == 1.) then
    !   y0 = scxtop(l) + piezom(x, h, l)
    !   uw = max(hsp(y0 - y), 0)
    ! else if(t > 0) then
    !   uw = t * vgms(t, a, n)
    ! end if
    ! u = ua - uw
  end subroutine
end module
