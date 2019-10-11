module scx
  use air, only: barom
  use bez, only: bezc
  use grndwt, only: hsp, piezom
  use interp1d, only: interp
  use vgmod, only: vgms
  implicit none
  private
  public scxcrk
  public scxcut
  public scxdim
  public scxdel
  public scxmat
  public scxnew
  public scxtop

  type aqua_t
    real h
    real l
    real y
  end type

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

  real, allocatable :: xs(:), ys(:)
  real :: scxdim(2), tana

contains
  elemental function scxcrk(x) result(y)
    real, intent(in) :: x
    real :: y

    y = 0
    ! y = yoc + tana * (x - xoc)
    ! z = 2. * c / w * tan(atan(1.) + .4 * phi)
  end function

  pure subroutine scxcut(n, rd, p, w, c, phi, u, alpha, b, h, stat)
    integer, intent(in) :: n
    real, intent(in) :: rd, p(:, :)
    real, intent(out), dimension(n) :: w, c, phi, u, alpha, b
    real, intent(out) :: h(0:n)
    integer, intent(out) :: stat
    real :: d
    real, dimension(2) :: m, q, r
    integer :: i

    d = rd * norm2(p(1, :) - p(size(p, 1), :))
    h(0) = 0
    q = p(1, :)
    stat = 0
    do i = 1, n
      r = bezc(real(i)/n, p)
      if(q(1) > r(1)) then
        stat = -i
        return
      end if

      b(i) = r(1) - q(1)
      h(i) = scxtop(r(1)) - r(2)
      if(h(i) < 0 .or. h(i) > d) then
        stat = i
        return
      end if

      alpha(i) = asin((q(2) - r(2)) / b(i))
      m = .5 * (q + r)
      call scxvar(m(1), m(2), c(i), phi(i), w(i), u(i))

      q = r
    end do
  end subroutine

  subroutine scxdel
  end subroutine

  pure subroutine scxmat(x, y, c, phi, w, u)
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

  subroutine scxnew(name)
    character(*), intent(in) :: name

    scxdim = (/ xs(size(xs)) - xs(1), maxval(ys) - minval(ys) /)
  end subroutine

  elemental function scxtop(x) result(y)
    real, intent(in) :: x
    real :: y

    y = interp(x, xs, ys)
  end function
end module
