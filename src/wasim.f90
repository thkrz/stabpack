module grndwa
  implicit none
  private

contains
  elemental function piezom(x, dh, l) result(h)
    real, intent(in) :: x, dh, l
    real :: h

    h = dh * sqrt(1. - (x / l)**2)
  end function
end module

module wasim
  use intgrt, only: nsimp
  use rtfind only: rtnewt
  use vg
  implicit none
  private
  public wasga

  type(vg_t) :: swc

contains
  pure function wasga(k, psi, time, delta) result(f)
    real, intent(in) :: k, psi, time, delta
    real :: f, kt, pd

    kt = k * time
    pd = abs(psi) * delta
    call rtnewt(funcd, kt, f)

  contains
    pure subroutine funcd(x, fval, fderiv)
      real, intent(in) :: x
      real, intent(out) :: fval, fderiv

      fval = x - pd * log(abs(1. + x / pd)) - kt
      fderiv = x / (pd + x)
    end subroutine
  end subroutine

  subroutine wasini(name, dt, x, y)
  end subroutine

  pure subroutine wasim(ep, ts, beta, ksat, i0, n, r, h, z)
    real, intent(in) :: ep(:, :), ts, beta(:, :)
    real, intent(in), dimension(:) :: ksat, i0, n, r, h
    real, intent(out) :: z(:)
    real, dimension(size(z)) :: k, psi
    real :: eps, t(0:size(z)), theta, v, zd
    integer :: i, ii, j, jj, m, num

    m = size(h)
    num = size(z)
    eps = 1. / (2. * num)
    i = num
    t(0) = 0
    theta = (wci - wcr) / (wcs - wcr)
    do j = 1, num
      jj = j - 1
      t(j) = real(j) / num
      if(abs(t(jj) - theta) < eps) i = jj
    end do

    z(:i) = h
    v = ep(1, 1)
    zd = ga(ks, psi_(theta), ts, n - wci)
    where (z == 0) z = zd
    v = v - zd

    do concurrent(j = i:num-1)
      jj = j + 1
      k(j) = ksat * nsimp(swc%relhc, t(j), t(jj))
      psi(j) = nsimp(swc%matsuc, t(j), t(jj))
    end do
  contains
    pure subroutine seep
      integer, save :: d = num

      dz = c * (k(d) * psi(d) / z(j) + k(d))
      z(j) = z(j) + dz * ts
    end subroutine
  end subroutine
end module
