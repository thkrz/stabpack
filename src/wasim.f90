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
  use spec, only: hyp2f1
  implicit none
  private

contains
  pure subroutine seep
    ! dz = 1. / (td - ti) * (K(td)*psi(td)/z(i) + K(td))
  end subroutine

  pure function k2(a, b)
    real, intent(in) :: a, b
    real :: k2

    k2 = (k(b) - k(a)) / (b - a)
  contains
    pure function k(t)
      real, intent(in) :: t
      real :: k, nn, tt

      tt = t**(n / (n - 1.))
      nn = 1. / n
      k = -t / a * (tt - 1.)**(1. + nn) * hyp2f1(1., 2., 2. - nn, tt)
    end function
  end function

  pure function psi2(a, b)
    real, intent(in) :: a, b
    real :: psi2

    psi2 = (psi(b) - psi(a)) / (b - a)
  contains
    pure function psi(t)
      real, intent(in) :: t
      real :: c1, c2, psi, nn, n3, tt, t32
      real :: h1, h2, h3, h4, h5

      nn = 1. / (n - 1.)
      tt = t**nn
      t32 = t**(3. / 2.)
      c1 = 6. * (n - 1.) * t**(5. / 2 + nn) / (5. * n - 3.)
      c2 = 3. * (n - 1.) * t**(7. / 2. + 2. * nn) / (7. * n - 3.)
      n3 = 3. / (2. * n)
      a = 5. / 2. - n3
      b = 1. / n
      c = 7. / 2. - n3
      h1 = hyp2f1(a, b, c, tt)
      b = 2. * b
      h2 = hyp2f1(a, b, c, tt)
      a = c
      c = 9. / 2. - n3
      h3 = hyp2f1(a, b, c, tt)
      a = 1. / n
      b = n3 * (n - 1.)
      c = 5. / 2. - n3
      h4 = hyp2f1(a, b, c, tt)
      a = a * 2.
      h5 = hyp2f1(a, b, c, tt)
      psi = 2. / 3. * (t32 + c1 * (h1 - h2) + c2 * h3 &
          - 2. * t32 * h4 + t32 * h5
    end function
  end function
end module
