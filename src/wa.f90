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

module vg
  implicit none
  private
  public vg_t

  type vg_t
    real a
    real n
  contains
    procedure :: se => vg_t_se
    procedure :: pm => vg_t_pm
    procedure :: kr => vg_t_kr
  end type

contains
  elemental function vg_t_se(self, h) result(se)
    class(vg_t), intent(in) :: self
    real, intent(in) :: h
    real :: m, se

    associate(a => self%a, n => self%n)
      m = 1. - 1. / n
      se = (1. + (a * h)**n)**(-m)
    end associate
  end function

  elemental function vg_t_pm(self, t) result(pm)
    class(vg_t), intent(in) :: self
    real, intent(in) :: t
    real :: m, me, pm

    associate(a => self%a, n => self%n)
      m = 1. - 1. / n
      me = a * (.046 * m + 2.07 * m**2 + 19.5 * m**3) &
         / (1. + 4.7 * m + 16. * m**2)
      pm = 1. / a * (t**(1 / m) - 1.)**(1 / n)
    end associate
    pm = max(pm, me)
  end function

  elemental function vg_t_kr(self, t) result(kr)
    class(vg_t), intent(in) :: self
    real, intent(in) :: t
    real :: m, kr

    associate(a => self%a, n => self%n)
      m = 1. - 1. / n
      kr = sqrt(t) * (1. - (1. - t**(1. / m))**m)**2
    end associate
  end function
end module
