module fmin
  implicit none
  private
  public amoeba

contains
  subroutine amoeba(fcn, x, fn, tol, stat)
    interface
      real function fcn(x)
        real, intent(in) :: x(:)
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
end module
