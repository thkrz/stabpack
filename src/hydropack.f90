module bc
  implicit none
  private
  public bcms
  public bcrhc
  public bcsm

contains
  elemental function bcms(t, a, n)
    real, intent(in) :: t, a, n
    real :: m, mt, bcms

    m = 1. - 1. / n
    mt = t**(1. / m)
    bcms = 1. / a * ((1 - mt) / mt)**(1. / n)
  end function

  elemental function bcrhc(t, n)
    real, intent(in) :: t, n
    real :: m, bcrhc

    m = 1. - 1. / n
    bcrhc = sqrt(t) * (1. - (1. - t**(1. / m))**m)**2
  end function

  elemental function bcsm(h, a, n)
    real, intent(in) :: h, a, n
    real :: m, bcsm

    m = 1. - 1. / n
    bcsm = (1. + (a * h)**n)**(-m)
  end function
end module

module fx
  implicit none
  private
  public fxms
  public fxrhc
  public fxsm

contains
  elemental function fxms(t, a, n)
    real, intent(in) :: t, a, n
    real :: m, mt, fxms

    m = 1. - 1. / n
    mt = t**(1. / m)
    fxms = 1. / a * ((1 - mt) / mt)**(1. / n)
  end function

  elemental function fxrhc(t, n)
    real, intent(in) :: t, n
    real :: m, fxrhc

    m = 1. - 1. / n
    fxrhc = sqrt(t) * (1. - (1. - t**(1. / m))**m)**2
  end function

  elemental function fxsm(h, a, n)
    real, intent(in) :: h, a, n
    real :: m, fxsm

    m = 1. - 1. / n
    fxsm = (1. + (a * h)**n)**(-m)
  end function
end module

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
