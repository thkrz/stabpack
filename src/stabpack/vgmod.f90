module vgmod
contains
  elemental function rhc(t, m)
    real, intent(in) :: t, m
    real :: rhc

    rhc = sqrt(t) * (1. - (1. - t**(1. / m))**m)**2
  end function

  elemental function wc(h, a, n)
    real, intent(in) :: h, a, n
    real :: m, wc

    m = 1. - 1. / n
    wc = (1. / (1. + (a * h)**n))**m
  end function
end module
