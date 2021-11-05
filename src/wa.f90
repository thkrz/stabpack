module wa_env
  implicit none

  real, parameter :: g = 9.81, rho = 1.
  real, parameter :: p0 = 101.325, T0 = 288.15

contains
  elemental function barom(h)
    real, intent(in) :: h
    real, parameter :: dT = .0065
    real :: barom

    barom = p0 * (1. - dT * h / T0)**5.255
  end function

  elemental function hystp(h, p1)
    real, intent(in) :: h
    real, optional, intent(in) :: p1
    real :: hystp

    hystp = rho * g * h
    if(present(p1)) hystp = hystp + p1
  end function

  elemental function rwc(t, s, r)
    real :: intent(in) :: t, s, r

    rwc = (t - r) / (s - r)
  end function

  elemental function vwc(t, s, r)
    real :: intent(in) :: t, s, r

    vwc = t * (s - r) + r
  end function
end module

module swcc
  implicit none
  private
  public vg_se
  public vg_pm
  public vg_kr

contains
  elemental function vg_se(h, a, n) result(se)
    real, intent(in) :: h, a, n
    real :: m, se

    m = 1. - 1. / n
    se = (1. + (a * h)**n)**(-m)
  end function

  elemental function vg_pm(t, a, n) result(pm)
    real, intent(in) :: t, a, n
    real :: m, me, pm

    m = 1. - 1. / n
    me = a * (.046 * m + 2.07 * m**2 + 19.5 * m**3) &
       / (1. + 4.7 * m + 16. * m**2)
    pm = 1. / a * (t**(1 / m) - 1.)**(1 / n)
    pm = max(pm, me)
  end function

  elemental function vg_kr(t, a, n) result(kr)
    real, intent(in) :: t, a, n
    real :: m, kr

    m = 1. - 1. / n
    kr = sqrt(t) * (1. - (1. - t**(1. / m))**m)**2
  end function
end module

module popr
  use wa_env
  use swcc
  implicit none

contains
  pure function pwp(z, t, n, r, beta, h) result(u)
    real, intent(in) :: z, t, n, r, beta(2)
    real, intent(in), optional :: h

    theta = rwc(t, n, r)
    ua = barom(h-z)
    if(theta == 1) then
      pm = ua - hystp(z, barom(h))
    else
      pm = vg_pm(theta, beta(1), beta(2))
    end if
    u = ua - pm
  end function
end module

