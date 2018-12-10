module hydro
contains
  subroutine hydest(k, t, s, x, r2)
    integer, parameter :: l = 4
    integer, intent(in) :: k
    real, intent(in) :: t(2, k), s
    real, intent(inout) :: x(l)
    real, intent(out) :: r2
    integer :: lwa, info, ipvt(l)
    real :: fvec(k), fjac(k, l), tol, wa(5 * l + k)
    real :: ss_res, ss_tot

    x(1) = .005
    x(2) = .5
    x(3) = 2.
    x(4) = t(2, k)
    tol = sqrt(epsilon(1.))
    lwa = size(wa)
    call lmder1(fcn, k, l, x, fvec, fjac, k, tol, info, ipvt, wa, lwa)

    ss_res = sum((t(2, :) - wcont(t(1, :), x(1), x(2), x(3), x(4), s))**2)
    ss_tot = sum((t(2, :) - sum(t(2, :))/size(t, 2))**2)
    r2 = 1. - ss_res / ss_tot

  contains
    subroutine fcn(k, l, x, fvec, fjac, ldfjac, iflag)
      integer, intent(in) :: k, l, ldfjac, iflag
      real, intent(in) :: x(l)
      real, intent(out) :: fvec(k), fjac(ldfjac, l)
      integer :: i

      if(iflag /= 2) then
        fvec = t(2, :) - wcont(t(1, :), x(1), x(2), x(3), x(4), s)
      else
        do concurrent(i = 1:k)
          fjac(i, :) = grad(t(1, i), x(1), x(2), x(3), x(4))
        end do
      end if
    end subroutine

    pure function grad(x, a, m, n, r)
      real, intent(in) :: x, a, m, n, r
      real :: grad(l)
      real :: c

      c = 1 / (1 + (a * x)**n)
      grad(1) = -m * n * (s - r) * x * (a * x)**(n - 1) * c**(m + 1)
      grad(2) = -(s - r) * c**m * log(c)
      grad(3) = m * (s - r) * (a * x)**n * c**(m + 1) * log(a * x)
      grad(4) = -1 + c**m
    end function
  end subroutine

  elemental function matsuc(t, a, m, n, r, s)
    real, intent(in) :: t, a, m, n, r, s
    real :: matsuc

    matsuc = (1 / a**n * (((s - r) / (t - r))**(1 / m) - 1))**(1 / n)
  end function

  elemental function relhc(t, m, r, s)
    real, intent(in) :: t, m, r, s
    real :: d, relhc

    d = (t - r) / (s - r)
    relhc = sqrt(d) * (1 - (1 - d**(1 / m))**m)**2
  end function

  elemental function wcont(h, a, m, n, r, s)
    real, intent(in) :: h, a, m, n, r, s
    real :: wcont

    wcont = r + (s - r) * (1 / (1 + (a * h)**n))**m
  end function
end module
