module vg
  implicit none
  private
  public vgms
  public vgrhc
  public vgsm

contains
  elemental function vgms(t, a, n)
    real, intent(in) :: t, a, n
    real :: m, mt, vgms

    m = 1. - 1. / n
    mt = t**(1. / m)
    vgms = 1. / a * ((1 - mt) / mt)**(1. / n)
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
