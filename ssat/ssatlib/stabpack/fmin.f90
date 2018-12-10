subroutine amoeba(fcn, x, n)
  integer, intent(in) :: n
  real, intent(inout) :: x(n)
  external :: fcn
  real :: c, alpha, beta, gamma, delta, &
          f(n+1), fc, fe, fr, ftmp, &
          xx(n, n+1), xc(n), &
          xe(n), xr(n), xtmp(n), &
          ident(n, n), p(n), p1, p2, facc
  integer :: i, j, k, maxit

  alpha = 1.
  beta = 1. + 2. / n
  gamma = .75 - 1. / (2 * n)
  delta = 1. - 1. / n
  c = 1.5
  p1 = c * (sqrt(n + 1.) + n - 1.) / (sqrt(2.) * n)
  p2 = c * (sqrt(n + 1.) - 1.) / (sqrt(2.) * n)
  p = p2
  ident = 0
  do i = 1, n
    ident(i, i) = 1
  end do
  xx(:, 1) = x
  do i = 2, n+1
    xx(:, i) = x + p + (p1 - p2) * ident(:, i-1)
  end do
  do j = 1, n+1
    f(j) = fcn(xx(:, j))
  end do
  facc = sqrt(epsilon(1.))
  maxit = 1000
  do i = 1, maxit
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
    ! if(present(fn)) fn = f(1)
    if(sqrt(sum(((f - fcn(x))**2) / n)) < facc) return
  end do
end subroutine
