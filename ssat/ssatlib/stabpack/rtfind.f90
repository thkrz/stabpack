module rtfind
  private
  public :: zbrent, rtnewt

  integer, parameter :: MAXIT = 100

contains
  pure subroutine rtnewt(funcd, x1, x2, x, tol, stat)
    interface
      pure subroutine funcd(x, fval, fderiv)
        real, intent(in) :: x
        real, intent(out) :: fval, fderiv
      end subroutine
    end interface
    real, intent(in) :: x1, x2
    real, intent(out) :: x
    real, intent(in), optional :: tol
    integer, intent(out), optional :: stat
    real :: df, dx, f, xacc
    integer :: j

    if(present(stat)) stat = 0
    xacc = merge(tol, sqrt(epsilon(1.)), present(tol))
    x = .5 * (x1 + x2)
    do j = 1, MAXIT
      call funcd(x, f, df)
      dx = f / df
      x = x - dx
      if((x1 - x) * (x2 - x) < 0.) then
        if(present(stat)) stat = -1
        return
      end if
      if(abs(dx) < xacc) return
    end do
    if(present(stat)) stat = -2
  end subroutine

  pure subroutine zbrent(fcn, x1, x2, b, tol, stat)
    interface
      pure function fcn(x)
        real, intent(in) :: x
        real :: fcn
      end function
    end interface
    real, intent(in) :: x1, x2
    real, intent(out) :: b
    real, intent(in), optional :: tol
    integer, intent(out), optional :: stat
    integer :: i
    real :: a, c, d, e, fa, fb, fc, min1, min2, &
            p, q, r, s, xacc, xm

    if(present(stat)) stat = 0
    a = x1
    b = x2
    c = 0
    d = 0
    e = 0
    fa = fcn(a)
    fb = fcn(b)
    if(fa * fb > 0) then
      if(present(stat)) stat = -1
      return
    end if
    fc = fb
    do i = 1, MAXIT
      if(fb * fc > 0) then
        c = a
        fc = fa
        d = b - a
        e = d
      end if
      if(abs(fc) < abs(fa)) then
        a = b
        b = c
        c = a
        fa = fb
        fb = fc
        fc = fa
      end if
      xacc = 2. * epsilon(x1) * abs(b) + .5 * merge(tol, sqrt(epsilon(1.)), present(tol))
      xm = .5 * (c - b)
      if(abs(xm) <= xacc .or. fb == 0.) return
      if(abs(e) >= xacc .and. abs(fa) > abs(fb)) then
        s = fb / fa
        if(a == b) then
          p = 2. * xm * s
          q = 1. - s
        else
          q = fa / fc
          r = fb / fc
          p = s * (2. * xm * q * (q - r) - (b - a) * (r - 1.))
          q = (q - 1.) * (r - 1.) * (s - 1.)
        end if
        if(p > 0) q = -q
        p = abs(p)
        min1 = 3. * xm * q - abs(xacc*q)
        min2 = abs(e*q)
        if(2.*p < min(min1, min2)) then
          e = d
          d = p / q
        else
          d = xm
          e = d
        end if
      else
        d = xm
        e = d
      end if
      a = b
      fa = fb
      if(abs(d) > xacc) then
        b = b + d
      else
        b = b + sign(xacc, xm)
      end if
      fb = fcn(b)
    end do
    if(present(stat)) stat = -2
  end subroutine
end module
