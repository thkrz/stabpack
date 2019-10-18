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
  use rtfind only: rtnewt
  use spec, only: simpn
  use swc
  implicit none
  private

  real, allocatable :: p(:)

contains

  pure subroutine seep(n, t0, dt, z0)
    integer, intent(in) :: n
    real, intent(in) :: t0, z0
    integer :: d, i, j, jj
    real :: dv, dz, eps, k(n), psi(n), t(0:n), z(n), zd

    eps = 1. / (2. * n)
    t = 0
    do j = 1, n
      t(j) = real(i) / n
      if(abs(t(j) - t0) < eps) i = j - 1
    end do

    do concurrent(j = i:n-1)
      k(j) = simpn(k2, t(j), t(j+1))
      psi(j) = simpn(psi2, t(j), t(j+1))
    end do

    d = 0
    dv = 1. / n
    z = zd
    do while(v > 0)
      do j = i, n
        if(v <= 0) exit
        if(d < j) d = j
        if(d == i) then
          c = 1. / (1. - t(i))
        else
          c = 1. / (t(d) .- t(i))
        end if
        dz = c * (k(d) * psi(d) / z(j) + k(d))
        z(j) = z(j) + dz * dt
      end do
    end do
  end subroutine

  pure function zd(k, psi, t, delta)
    real, intent(in) :: k, psi, t, delta
    real :: kt, pd, zd

    kt = k * t
    pd = abs(psi) * delta
    call rtnewt(funcd, kt, zd)

  contains
    pure subroutine funcd(x, fval, fderiv)
      real, intent(in) :: x
      real, intent(out) :: fval, fderiv

      fval = x - pd * log(abs(1. + x / pd)) - kt
      fderiv = x / (pd + x)
    end subroutine
  end subroutine

end module
