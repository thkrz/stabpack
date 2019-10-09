module razdol
  implicit none
  private
  public razslv

contains
  pure subroutine razslv(n, w, c, phi, u, alpha, b, h, f, e, t, mu, p0)
    integer, intent(in) :: n
    real, dimension(n), intent(in) :: w, c, phi, u, alpha, b
    real, intent(in) :: h(n + 1), f
    real, dimension(n + 1), intent(out) :: e, t
    real, intent(out) :: mu
    real, optional, intent(in) :: p0
    integer :: i, j
    real :: phim, tana, tana2

    e = 0
    t = 0
    if(present(p0)) then
      e(1) = p0
      t(1) = p0
    end if

    do i = 2, n + 1
      j = i - 1
      tana = tan(alpha(j))
      tana2 = tana**2
      phim = tan(phi(j)) / f
      e(i) = 2. * t(j) * (tana - phim) + e(j) * (1. + 2. * tana * phim - tana2) &
           + w(j) * b(j) / 3. * (h(i) + 2. * h(j)) * (tana - phim) &
           - b(j) * (1. + tana2) * (c(j) / f - u(j) * phim)
      e(i) = e(i) / (1. + tana2)
      t(i) = (e(i) + e(j)) * tana - t(j) + w(j) * b(j) / 6. * (h(i) - h(j))
    end do
    mu = 1. - e(n + 1) / maxval(e)
  end subroutine

  pure function razsld(n, w, c, phi, u, alpha, b, h, f, e, t) result(mu)
    integer, intent(in) :: n
    real, dimension(n), intent(in) :: w, c, phi, u, alpha, b
    real, intent(in) :: h(n + 1), f
    real, dimension(n + 1), intent(in) :: e, t
    integer :: i, j
    real, dimension(n + 1) :: de, dt
    real :: cosa2, mu, tana, tanp

    de = 0
    dt = 0
    do i = 2, n + 1
      j = i - 1
      cosa2 = cos(alpha(j))**2
      tana = tan(alpha(j))
      tanp = tan(phi(j))
      de(i) = cosa2 * (6. * f * dt(j) * (f * tana - tanp) &
            + (2. * b(j) * w(j) * h(j) + b(j) * w(j) * h(i) &
            + 6. * t(j) - 6. * e(j) * tana) * tanp + 3. * b(j)/cosa2 &
            * (c(j) - u(j) * tanp) + 3. * f * de(j) &
            * (f - f * tana**2 + 2. * tana * tanp))
      de(i) = de(i) / (3. * f**2)
      dt(i) = (de(i) + de(j)) * tana + dt(j)
    end do
    i = maxloc(e, 1)
    mu = (de(i) * e(n + 1) - de(n + 1) * e(i)) / e(i)**2
  end function
end module
