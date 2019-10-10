module grndwt
  implicit none
  private
  public piezom

contains
  elemental function piezom(x, h0, l) result(h)
    real, intent(in) :: x, h0, l
    real :: h

    h = h0 * sqrt(1. - (x / l)**2)
  end function
end module

module vgmod
  implicit none
  private
  public vgms
  public vgrhc
  public vgsm

contains
  elemental function vgms(t, a, n)
    real, intent(in) :: t, a, n
    real :: m, me, mt, vgms

    m = 1. - 1. / n
    me = a * (.046 * m + 2.07 * m**2 + 19.5 * m**3) &
       / (1. + 4.7 * m + 16. * m**2)
    mt = t**(1. / m)
    vgms = max(1. / a * ((1 - mt) / mt)**(1. / n), me)
  end function

  elemental function vgrhc(t, n)
    real, intent(in) :: t, n
    real :: m, vgrhc

    m = 1. - 1. / n
    vgrhc = sqrt(t) * (1. - (1. - t**(1. / m))**m)**2
  end function

  elemental function vgsm(h, a, n)
    real, intent(in) :: h, a, n
    real :: m, vgsm

    m = 1. - 1. / n
    vgsm = (1. + (a * h)**n)**(-m)
  end function
end module
