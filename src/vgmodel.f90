module vgmodel
  implicit none
  private
  public vg_rhc
  public vg_wc

contains
  elemental function vg_rhc(t, n)
    real, intent(in) :: t, n
    real :: m, vg_rhc

    m = 1. - 1. / n
    vg_rhc = sqrt(t) * (1. - (1. - t**(1. / m))**m)**2
  end function

  elemental function vg_wc(h, a, n)
    real, intent(in) :: h, a, n
    real :: m, vg_wc

    m = 1. - 1. / n
    vg_wc = (1. / (1. + (a * h)**n))**m
  end function
end module
