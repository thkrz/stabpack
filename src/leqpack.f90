module gle
  implicit none
  private
  public gleslv

contains
  pure function gleslv(n, lambda, model)
    integer, intent(in) :: n
    real, intent(in) :: lambda
    character(2), intent(in) :: model
    real :: theta(n + 1)

    theta = lambda
    theta(1) = 0
    theta(n + 1) = 0
  end function

  pure subroutine ff(x, fval, fderiv)
    real, intent(in) :: x
    real, intent(out) :: fval, fderiv

  end subroutine
  pure function fos(x)
    integer :: i
    real :: fos, z(n + 1)

    z(1) = 0
    do i = 2, n + 1
      z(i) = zn(x, z(i - 1), alpha, theta(i - 1:i), k, c, phi, w, u)
    end do
    fos = z(n + 1)
  end function

  pure function zn(f, z0, alpha, theta, k, c, phi, w, u)
    real, intent(in) :: f, z0, alpha, theta(2), k(2), c, phi, w, u
    real :: cm, phim, sec, zn

    cm = c / f
    phim = phi / f
    sec = 1. / cos(alpha - theta(2) - phim)
    zn = sec * (-w * sin(alpha - phim) - u * sin(phim) + cos(phim) &
       * cm - w * cos(alpha - phim) * k(1) + w * sin(alpha - phim) &
       * k(2) + cos(alpha - theta(1) - phim) * z0)
  end function

  pure function dzn(f, dz0, z0, alpha, theta, k, c, phi, w, u)
    real, intent(in) :: f, dz0, z0, alpha, theta(2), k(2), c, phi, w, u
    real :: cm, dzn, f3, phim, phi2, sec

    cm = c / f
    phim = phi / f
    phi2 = 2. * phi
    f3 = 2. * f**3
    sec = 1. / cos(alpha - phim - theta(2))
    dzn = sec**2 / f3 * (-f * ((c - u * phi2) &
        * cos(alpha - theta(2)) + c * cos(alpha - 2. * phim - theta(1)) &
        + w * phi2 * cos(theta(2))) + c * phi2 * sin(alpha - theta(2)) &
        + f3 * cos(alpha - phim - theta(2)) &
        * cos(alpha - phim - theta(1)) * dz0 + f * phi2 &
        * (w * sin(theta(2)) * k(1) + w * cos(theta(2)) &
        * k(2) - sin(theta(2) - theta(1)) * z0))
  end function

  pure function hn(z, h0, alpha, theta, k, b, hc, w)
    real, intent(in) :: z(2), h0, alpha, theta(2), k(2), b, hc, w
    real :: bcos, hn, sec

    bcos = b / cos(alpha)
    sec = 1. / cos(theta(2))
    hn = -sec / (2. * z(2)) * (2. * w * hc * k(1) &
       + (bcos * sin(alpha + theta(1)) &
       - 2. * cos(theta(1)) * h0) * z(1) + bcos &
       * sin(alpha - theta(2)) * z(2))
  end function

  pure function dhn(z, dh0, h0, alpha, theta, k, b, hc, w)
    real, intent(in) :: z(2), dh0, h0, alpha, theta(2), k(2), b, hc, w
    real :: dhn, sec

    sec = 1. / cos(theta(2))
    dhn = sec / (2. * z(2)) * (2. * cos(theta(1)) * dh0 * z(1) &
        + b * sec * z(2) - (2. * w * hc * k(1) + (b / cos(alpha) &
        * sin(alpha + theta(1) - 2. * cos(theta(1)) * h0) * z(1)) &
        * tan(theta(2))
  end function
end module
