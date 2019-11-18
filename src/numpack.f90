module bez
  implicit none
  private
  public bezarc
  public bezcrv
  public bezfit
  public bezlin
  public bezsim

contains
  pure subroutine bezarc(a, b, p)
    real, parameter :: k = 4. / 3. * (sqrt(2.) - 1.), &
                       pi = 4. * atan(1.)
    real, intent(in) :: a(:), b(:)
    real, intent(out) :: p(:, :)
    integer :: i, m, n
    real :: br(2), dx(2)
    real :: beta, cosb, sinb, r, rot(2, 2)

    n = size(p, 2)
    if (mod(n, 2) /= 0) error stop

    m = n / 2
    dx = b - a
    beta = .5 * pi - atan(dx(2) / dx(1))
    cosb = cos(beta)
    sinb = sin(beta)
    rot(1, 1) = cosb
    rot(1, 2) = -sinb
    rot(2, 1) = sinb
    rot(2, 2) = cosb
    br = matmul(rot, dx) + a

    cosb = cos(.25 * pi)
    sinb = cosb
    r = abs(a(2) - br(2)) / (2. * sinb)
    do i = 1, m - 1
      p(:, i) = a
    end do
    dx =  k * r * (/ sinb, cosb /)
    rot = inv(rot)
    p(:, m) = matmul(rot, dx) + a
    p(:, m + 1) = matmul(rot, br + dx - a) + a
    do i = m + 2, n
      p(:, i) = b
    end do
  end subroutine

  pure function bezcrv(t, p) result(c)
    real, intent(in) :: t, p(:, :)
    real :: c(size(p, 1))
    integer :: i, n

    c = 0
    if(t < 0 .or. t > 1) error stop
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

  pure recursive function bezfit(x, tol) result(p)
    real, intent(in) :: x(:, :), tol
    real :: a(size(x, 1)), a1, a2, a12, d, t, t1
    real, allocatable :: p(:, :), q(:, :), r(:, :)
    real, dimension(2) :: c1, c2, c12
    integer :: i, k, n

    n = size(x, 2)
    if (n < 3) error stop
    k = n - 1
    allocate(p(2, 4))
    p(:, 1) = x(:, 1)
    p(:, 4) = x(:, n)

    a1 = 0
    a2 = 0
    a12 = 0
    c1 = 0
    c2 = 0
    do i = 0, k
      t = real(i) / k
      t1 = 1. - t
      a1 = a1 + t**2 * t1**4
      a2 = a2 + t**4 * t1**2
      a12 = a12 + t**3 * t1**3
      c12 = x(:, i + 1) - t1**3 * p(1, :) - t**3 * p(:, 4)
      c1 = c1 + 3. * t * t1**2 * c12
      c2 = c2 + 3. * t**2 * t1 * c12
    end do
    a1 = 9. * a1
    a2 = 9. * a2
    a12 = 9. * a12
    d = (a1 * a2 - a12 * a12)
    p(:, 2) = (a2 * c1 - a12 * c2) / d
    p(:, 3) = (a1 * c2 - a12 * c1) / d

    do concurrent(i = 0:k)
      t = real(i) / k
      a(i + 1) = norm2(x(:, i + 1) - bezcrv(t, p))
    end do
    if (all(a < tol)) return

    k = maxloc(a, 1)
    if (k < 3 .or. n - k < 2) return

    q = bezfit(x(:, :k), tol)
    r = bezfit(x(:, k:), tol)
    n = size(q, 1)
    k = size(r, 1)
    deallocate(p)
    allocate(p(n + k - 1, 2))
    p(:, :n) = q
    p(:, n:) = r
    deallocate(q)
    deallocate(r)
  end function

  pure subroutine bezlin(a, b, p)
    real, intent(in) :: a(:), b(:)
    real, intent(out) :: p(:, :)
    real :: dx(size(a))
    integer :: i, j, n

    dx = b - a
    n = size(p, 2) - 1
    do i = 0, n
      j = i + 1
      p(:, j) = dx * real(i) / n + a
    end do
  end subroutine

  pure function bezsim(p, q) result(lambda)
    real, intent(in), dimension(:, :) :: p, q
    real :: lambda
    integer :: i, n

    n = size(p, 2)
    if(n /= size(q, 2)) error stop
    lambda = 0
    do i = 1, n
      lambda = lambda + norm2(p(:, i) - q(:, i))
    end do
    lambda = lambda / n
  end function

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
  public brent
  public mnbrak

contains
  subroutine amoeba(fcn, x, fn, tol, stat)
    interface
      function fcn(x)
        real, intent(in) :: x(:)
        real :: fcn
      end function
    end interface
    real, intent(inout) :: x(:)
    real, intent(out) :: fn
    real, intent(in), optional :: tol
    integer, intent(out), optional :: stat
    real :: c, alpha, beta, gamma, delta, &
            f(size(x)+1), fc, fe, fr, ftmp, &
            xx(size(x), size(x)+1), xc(size(x)), &
            xe(size(x)), xr(size(x)), xtmp(size(x)), &
            ident(size(x), size(x)), p(size(x)), p1, p2, facc
    integer :: i, j, k, itmax, n

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
    do j = 1, n+1
      f(j) = fcn(xx(:, j))
    end do
    facc = merge(tol, sqrt(epsilon(1.)), present(tol))
    if(present(stat)) stat = 0
    itmax = n * 200
    do i = 1, itmax
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
      fn = f(1)
      if(sqrt(sum(((f - fcn(x))**2) / n)) < facc) return
    end do
    if(present(stat)) stat = -1
  end subroutine

  function brent(ax, bx, cx, func, tol, xmin, stat)
    real, intent(in) :: ax, bx, cx, tol
    real, intent(out) :: xmin
    integer, intent(out), optional :: stat
    real :: brent
    interface
      function func(x)
        real, intent(in) :: x
        real :: func
      end function
    end interface
    real, parameter :: cgold = 0.3819660, zeps = 1.0e-3 * epsilon(ax)
    integer, parameter :: itmax = 100
    integer :: i
    real :: a, b, c, d, e, etemp, fu, fv, fw, fx, p, q, r, tol1, tol2, u, v, w, x, xm

    if(present(stat)) stat = 0
    a = min(ax, cx)
    b = max(ax, cx)
    v = bx
    w = v
    x = v
    e = 0
    fx = func(x)
    fv = fx
    fw = fx
    do i = 1, itmax
      xm = .5 * (a + b)
      tol1 = tol * abs(x) + zeps
      tol2 = 2. * tol1
      if(abs(x - xm) <= (tol2 - .5 * (b - a))) then
        xmin = x
        brent = fx
        return
      end if
      if(abs(e) > tol1) then
        r = (x - w) * (fx - fv)
        q = (x - v) * (fx - fw)
        p = (x - v) * q - (x - w) * r
        q = 2. * (q - r)
        if(q > 0) p = -p
        q = abs(q)
        etemp = e
        e = d
        if(abs(p) >= abs(.5 * q * etemp) .or. &
          p <= q * (a - x) .or. p >= q * (b - x)) then
          e = merge(a - x, b - x, x >= xm)
          d = cgold * e
        else
          d = p / q
          u = x + d
          if(u - a < tol2 .or. b - u < tol2) d = sign(tol1, xm - x)
        end if
      else
        e = merge(a - x, b - x, x >= xm)
        d = cgold * e
      end if
      u = merge(x + d, x + sign(tol1, d), abs(d) >= tol1)
      fu = func(u)
      if(fu <= fx) then
        if(u >= x) then
          a = x
        else
          b = x
        end if
        call shft(v, w, x, u)
        call shft(fv, fw, fx, fu)
      else
        if(u < x) then
          a = u
        else
          b = u
        end if
        if(fu <= fw .or. w == x) then
          v = w
          fv = fw
          w = u
          fw = fu
        else if(fu <= fv .or. v == x .or. v == w) then
          v = u
          fv = fu
        end if
    end do
    if(present(stat)) stat = -1
  end function

  subroutine mnbrak(ax, bx, cx, fa, fb, fc, func)
    real, intent(inout) :: ax, bx
    real, intent(out) :: cx, fa, fb, fc
    interface
      function func(x)
        real, intent(in) :: x
        real :: func
      end function
    end interface
    real, parameter :: gold = 1.618034, glimit = 100.0
    real :: fu, q, r, u, ulim

    fa = func(ax)
    fb = func(ab)
    if(fb > fa) then
      call swap(ax, bx)
      call swap(fa, fb)
    end if
    cx = bx + gold * (bx - ax)
    fc = func(cx)
    do
      if(fb < fc) return
      r = (bx - ax) * (fb - fc)
      q = (bx - cx) * (fb - fa)
      u = bx - ((bx - cx) * q - (bx - ax) * r) / (2. * sign(max(abs(q - r), tiny(1.)), q - r))
      ulim = bx + glimit * (cx - bx)
      if((bx - u) * (u - cx) > 0.) then
        fu = func(u)
        if(fu < fc) then
          ax = bx
          fx = fb
          bx = u
          fb = fu
          return
        else if(fu > fb) then
          cx = u
          fc = fu
          return
        end if
        u = cx + gold * (cx - bx)
        fu = func(u)
      else if((cx - u) * (u - ulim) > 0.) then
        fu = func(u)
        if(fu < fc) then
          bx = cx
          cx = u
          u = cx + gold * (cx - bx)
          call shft(fb, fc, fu, func(u))
        end if
      else if((u - ulim) * (ulim - cx) >= 0.) then
        u = ulim
        fu = func(u)
      else
        u = cx + gold * (cx - bx)
        fu = func(u)
      end if
      call shft(ax, bx, cx, u)
      call shft(fa, fb, fc, fu)
    end do
  end subroutine

  pure subroutine shft(a, b, c, d)
    real, intent(out) :: a
    real, intent(inout) :: b, c
    real, intent(in) :: d

    a = b
    b = c
    c = d
  end subroutine

  pure subroutine swap(a, b)
    real, intent(inout) :: a, b
    real :: c

    c = a
    a = b
    b = c
  end subroutine
end module

module intp1d
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

  integer, parameter :: itmax = 1000

contains
  pure subroutine rtnewt(funcd, p, x0, x, tol, stat)
    interface
      pure subroutine funcd(x, p, fval, fderiv)
        real, intent(in) :: x, p(:)
        real, intent(out) :: fval, fderiv
      end subroutine
    end interface
    real, intent(in) :: p(:), x0
    real, intent(out) :: x
    real, intent(in), optional :: tol
    integer, intent(out), optional :: stat
    real :: df, dx, f, xacc
    integer :: j

    if(present(stat)) stat = 0
    xacc = merge(tol, sqrt(epsilon(1.)), present(tol))
    x = x0
    do j = 1, itmax
      call funcd(x, p, f, df)
      dx = f / df
      x = x - dx
      if(abs(dx) < xacc) return
    end do
    if(present(stat)) stat = -1
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
    do i = 1, itmax
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

module spcfnc
  implicit none
  private
  public hyp2f1

contains
  elemental function hyp2f1(a, b, c, z) result(f1)
    integer, parameter :: itmax = 1000
    real, intent(in) :: a, b, c, z
    real :: aa, bb, cc, f1, fac, temp
    integer :: n

    fac = 1.
    temp = fac
    aa = a
    bb = b
    cc = c
    do n = 1, itmax
      fac = fac * ((aa * bb) / cc) * z / n
      f1 = temp + fac
      if(f1 == temp) return
      temp = f1
      aa = aa + 1.
      bb = bb + 1.
      cc = cc + 1.
    end do
  end function
end module
