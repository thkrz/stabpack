module bez
  use ssat_env, only: ndim, pi
  implicit none
  private
  public beza
  public bezc
  public bezl

contains
  pure subroutine beza(p)
    real, parameter :: k = 4. / 3. * (sqrt(2.) - 1.)
    real, intent(inout) :: p(:, :)
    integer :: i, m, n
    real, dimension(ndim) :: a, b, dx
    real :: beta, cosb, sinb, r, rot(ndim, ndim)

    n = size(p, 1)
    m = n / 2
    dx = p(n, :) - p(1, :)
    beta = .5 * pi - atan(dx(2) / dx(1))
    cosb = cos(beta)
    sinb = sin(beta)
    rot(1, 1) = cosb
    rot(1, 2) = -sinb
    rot(2, 1) = sinb
    rot(2, 2) = cosb
    a = p(1, 1)
    b = matmul(rot, dx) + a

    cosb = cos(.25 * pi)
    sinb = cosb
    r = abs(a(2) - b(2)) / (2. * sinb)
    do i = 2, m - 1
      p(i, :) = p(1, :)
    end do
    dx =  k * r * (/ sinb, cosb /)
    p(m, :) = matmul(inv(rot), dx) + a
    p(m + 1, :) = matmul(inv(rot), b + dx - a) + a
    do i = m + 2, n - 1
      p(i, :) = p(n, :)
    end do
  end subroutine

  pure function bezc(t, p) result(c)
    real, intent(in) :: t, p(:, :)
    real :: c(size(p, 1))
    integer :: i, n

    c = 0
    if(t < 0 .or. t > 1) return
    n = size(p, 2) - 1
    do i = 0, n
      c = c + b(t, i, n) * p(:, i+1)
    end do

  contains
    pure function b(t, i, n)
      real, intent(in) :: t
      integer, intent(in) :: i, n
      real :: b

      b = bico(n, i) * t**i * (1 - t)**(n - i)
    end function

    pure function bico(n, k)
      integer, intent(in) :: n, k
      integer :: bico

      bico = factln(n) / (factln(k) * factln(n-k))
    end function

    pure function factln(n)
      integer, intent(in) :: n
      integer :: i, factln

      factln = product([(i, i=1,n)])
    end function
  end function

  pure subroutine bezl(p)
    real, intent(inout) :: p(:, :)
    integer :: i, n

    n = size(p, 1)
    do i = 2, n - 1
      p(i, :) = (p(n, :) - p(1, :)) * real(i) / n + p(1, :)
    end do
  end subroutine

  pure function inv(a) result(b)
    real, intent(in) :: a(2, 2)
    real :: b(2, 2)
    real :: detinv

    detinv = 1. / (a(1, 1) * a(2, 2) - a(1, 2) * a(2, 1))

    b(1, 1) = +detinv * a(2, 2)
    b(2, 1) = -detinv * a(2, 1)
    b(1, 2) = -detinv * a(1, 2)
    b(2, 2) = +detinv * a(1, 1)
  end function
end module

module fmin
  implicit none
  private
  public amoeba

  integer, parameter :: MAXIT = 1000

contains
  pure subroutine amoeba(fcn, x, tol, fn, stat)
    interface
      pure function fcn(x)
        real, intent(in) :: x(:)
        real :: fcn
      end function
    end interface
    real, intent(inout) :: x(:)
    real, intent(in), optional :: tol
    real, intent(out), optional :: fn
    integer, intent(out), optional :: stat
    real :: c, alpha, beta, gamma, delta, &
            f(size(x)+1), fc, fe, fr, ftmp, &
            xx(size(x), size(x)+1), xc(size(x)), &
            xe(size(x)), xr(size(x)), xtmp(size(x)), &
            ident(size(x), size(x)), p(size(x)), p1, p2, facc
    integer :: i, j, k, n

    n = size(x)
    alpha = 1.
    beta = 1. + 2. / n
    gamma = .75 - 1. / (2 * n)
    delta = 1. - 1. / n
    c = 1.5
    p1 = c * (sqrt(n + 1.) + n - 1.) / (sqrt(2.) * n)
    p2 = c * (sqrt(n + 1.) - 1.) / (sqrt(2.) * n)
    p = p2
    ident = 0
    do concurrent(i = 1:n)
      ident(i, i) = 1
    end do
    xx(:, 1) = x
    do concurrent(i = 2:n+1)
      xx(:, i) = x + p + (p1 - p2) * ident(:, i-1)
    end do
    do concurrent(j = 1:n+1)
      f(j) = fcn(xx(:, j))
    end do
    facc = merge(tol, sqrt(epsilon(1.)), present(tol))
    if(present(stat)) stat = 0
    do i = 1, MAXIT
      do j = 2, n+1
        k = j - 1
        ftmp = f(j)
        xtmp = xx(:, j)
        do while(k >= 1 .and. f(k) > ftmp)
          f(k+1) = f(k)
          xx(:, k+1) = xx(:, k)
          k = k - 1
        end do
        f(k+1) = ftmp
        xx(:, k+1) = xtmp
      end do
      x = sum(xx(:, 1:n), 2) / n
      xr = x + alpha * (x - xx(:, n+1))
      fr = fcn(xr)
      if(fr >= f(1) .and. fr <= f(n)) then
        xx(:, n+1) = xr
        f(n+1) = fr
      else if(fr < f(1)) then
        xe = x + beta * (xr - x)
        fe = fcn(xe)
        if(fe < f(1)) then
          xx(:, n+1) = xe
          f(n+1) = fe
        else
          xx(:, n+1) = xr
          f(n+1) = fr
        end if
      else if(fr > f(n)) then
        if(fr >= f(n+1)) then
          xc = x + gamma * (xx(:, n+1) - x)
        else
          xc = x + gamma * (xr - x)
        end if
        fc = fcn(xc)
        if(fc < f(n+1)) then
          xx(:, n+1) = xc
          f(n+1) = fc
        else
          do j = 2, n+1
            xx(:, j) = delta * (xx(:, j) - xx(:, 1)) + xx(:, 1)
            f(j) = fcn(xx(:, j))
          end do
        end if
      end if
      if(present(fn)) fn = f(1)
      if(sqrt(sum(((f - fcn(x))**2) / n)) < facc) return
    end do
    if(present(stat)) stat = -1
  end subroutine
end module

module interp1d
  implicit none
  private
  public interp

contains
  pure function interp(x, xp, fp) result(y)
    real, intent(in) :: x, xp(:), fp(:)
    integer :: i, j, n
    real :: y

    n = size(xp)
    y = fp(1)
    if(x < xp(1)) return
    do i = 1, n - 1
      j = i + 1
      if(x >= xp(i) .and. x < xp(j)) then
        y = fp(i) + (x - xp(i)) * (fp(j) - fp(i)) / (xp(j) - xp(i))
        return
      end if
    end do
    if(x > xp(n)) y = fp(n)
  end function
end module

module rtfind
  implicit none
  private
  public rtnewt
  public zbrent

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
