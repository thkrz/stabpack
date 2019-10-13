module pwp
  implicit none
  private
  public barom
  public hystp

  real, parameter :: g = 9.81, rho = 1.
  real, parameter :: p0 = 101.325, T0 = 288.15

contains
  elemental function barom(h)
    real, parameter :: dT = .0065
    real, intent(in) :: h
    real :: barom

    barom = p0 * (1. - dT * h / T0)**5.255
  end function

  elemental function hystp(h, p)
    real, intent(in) :: h
    real, intent(in), optional :: p
    real :: hystp

    hystp = rho * g * h
    if(present(p)) hystp = hystp + p
  end function
end module

module swc
  implicit none
  private
  public swcms
  public swcrhc
  public swcewc
  public swcset

  abstract interface
    pure function fcn(x, p)
      real, intent(in) :: x, p(:)
      real :: fcn
    end function
  end interface

  procedure(fcn), pointer :: swcms => vgms
  procedure(fcn), pointer :: swcrhc => vgrhc
  procedure(fcn), pointer :: swcewc => vgewc

contains
  subroutine swcset(model)
    character(len=2), intent(in) :: model

    if(model == 'bc') then
      swcms => bcms
      swcrhc => bcrhc
      swcewc => bcewc
    else if(model == 'fx') then
      swcms => null()
      swcrhc => null()
      swcewc => null()
    else if(model == 'vg') then
      swcms => vgms
      swcrhc => vgrhc
      swcewc => vgewc
    end if
  end subroutine

  pure function bcms(h, p)
    real, intent(in) :: h, p(:)
    real :: bcms

    associate(hb => p(1), lambda => p(2))
      bcms = h
    end associate
  end function

  pure function bcrhc(h, p)
    real, intent(in) :: h, p(:)
    real :: bcrhc

    bcrhc = 1
    associate(hb => p(1), lambda => p(2))
      if(h >= hb) return
      bcrhc = (h / hb)**(-2. - 3. * lambda)
    end associate
  end function

  pure function bcewc(h, p)
    real, intent(in) :: h, p(:)
    real :: bcewc

    bcewc = 1
    associate(hb => p(1), lambda => p(2))
      if(h >= hb) return
      bcewc = (h / hb)**(-lambda)
    end associate
  end function

  pure function vgms(t, p)
    real, intent(in) :: t, p(:)
    real :: m, me, vgms

    associate(a => p(1), n => p(2))
      m = 1. - 1. / n
      me = a * (.046 * m + 2.07 * m**2 + 19.5 * m**3) &
         / (1. + 4.7 * m + 16. * m**2)
      vgms = max(1. / a * (t**(1 / m) - 1.)**(1 / n), me)
    end associate
  end function

  pure function vgrhc(t, p)
    real, intent(in) :: t, p(:)
    real :: m, vgrhc

    associate(n => p(2))
      m = 1. - 1. / n
      vgrhc = sqrt(t) * (1. - (1. - t**(1. / m))**m)**2
    end associate
  end function

  pure function vgewc(h, p)
    real, intent(in) :: h, p(:)
    real :: m, vgewc

    associate(a => p(1), n => p(2))
      m = 1. - 1. / n
      vgewc = (1. + (a * h)**n)**(-m)
    end associate
  end function
end module
