module fmin
  use math, only: shft, swap
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
    real, optional, intent(in) :: tol
    integer, optional, intent(out) :: stat
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
    integer, optional, intent(out) :: stat
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
    real :: a, b, d, e, etemp, fu, fv, fw, fx, p, q, r, tol1, tol2, u, v, w, x, xm

    if(present(stat)) stat = 0
    brent = bx
    a = min(ax, cx)
    b = max(ax, cx)
    v = bx
    w = v
    x = v
    e = 0.
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
    fb = func(bx)
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
          fa = fb
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
end module
