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

module wa
  use intgrt, only: nsimp
  use rtfind only: rtnewt
  use swc
  implicit none
  private

contains
  pure function ga(k, psi, t, delta)
    real, intent(in) :: k, psi, t, delta
    real :: ga, kt, pd

    kt = k * t
    pd = abs(psi) * delta
    call rtnewt(funcd, kt, ga)

  contains
    pure subroutine funcd(x, fval, fderiv)
      real, intent(in) :: x
      real, intent(out) :: fval, fderiv

      fval = x - pd * log(abs(1. + x / pd)) - kt
      fderiv = x / (pd + x)
    end subroutine
  end subroutine

  pure subroutine wasim(num, ts, ks, n, wap, wci, wcr, wcs)
    integer, intent(in) :: num
    real, intent(in) :: ts, ks, n, wap(:), wci, wcr, wcs
    real, dimension(num) :: k, psi, z
    real :: eps, t(0:num), theta
    integer :: i, j, jj

    eps = 1. / (2. * num)
    i = num
    t(0) = 0
    theta = (wci - wcr) / (wcs - wcr)
    z = ga(ks, psi_(1.), ts, n - wci)
    do j = 1, num
      jj = j - 1
      t(j) = real(j) / num
      if(abs(t(jj) - theta) < eps) i = jj
    end do

    do concurrent(j = i:num-1)
      jj = j + 1
      k(j) = nsimp(k_, t(j), t(jj))
      psi(j) = nsimp(psi_, t(j), t(jj))
    end do

  contains
    pure function k_(x)
      real, intent(in) :: x
      real :: k_

      k_ = ks * swcrhs(x, wap)
    end function

    pure function psi_(x)
      real, intent(in) :: x
      real :: psi_

      psi_ = swcms(x, wap)
    end function

    pure subroutine seep
      integer, save :: d = num

      dz = c * (k(d) * psi(d) / z(j) + k(d))
      z(j) = z(j) + dz * ts
    end subroutine
  end subroutine
end module
