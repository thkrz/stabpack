module clusan
  implicit none
  private
  public kmeans

  integer, parameter :: MAXIT = 10000

contains
  pure subroutine kmeans(x, m)
    real, intent(in) :: x(:, :)
    real, intent(inout) :: m(:, :)
    real :: tol, var, varmin
    integer :: c(size(m, 1), size(m, 2)), s(size(x, 2))
    integer :: i, j, k, num, n, z

    tol = sqrt(epsilon(1.))
    k = size(m, 2)
    n = size(x, 2)
    do z = 1, MAXIT
      do j = 1, n
        varmin = -1.
        do i = 1, k
          var = norm2(x(:, j) - m(:, i))
          if(varmin > var .or. varmin < 0) then
            varmin = var
            s(j) = i
          end if
        end do
      end do
      do concurrent(i = 1:k)
        num = count(s == i)
        do concurrent(j = 1:size(m, 1))
          c(j, i) = sum(x(j, :), s == i) / num
        end do
      end do
      if(norm2(c - m) < tol) exit
      m = c
    end do
  end subroutine
end module

module rvcont
  implicit none
  private

contains
  elemental function loglap_cdf(x, d, a, b) result(c)
    real, intent(in) :: x, d, a, b
    real :: ab, c, p

    c = 0
    if(x < 0) return

    ab = a / (a + b)
    if(x < d) then
      c = ab * (x / d)**b
    else
      c = 1. - ab * (d/ x)**a
    end if
  end function

  pure subroutine loglap_fit(y, d, a, b)
    real, intent(in), dimension(:) :: x, y
    real, intent(out) :: d, a, b
    real :: imv(3, 3)

  contains
    pure function h(d)
      real, intent(in) :: d
      real :: h

      call pest(d, a, b)
      h = log(a + b) + b + (a * b) / (a + b)
    end function

    pure subroutine pest(d, a, b)
      real, intent(in) :: d
      real, intent(out) :: a, b
      real :: pd, py, sqlp
      integer :: i, n

      n = size(y)
      pd = 1
      py = 1
      do i = 1, n
        pd = pd * merge(y(i) / d, 1., d /= 0)
        py = py * merge(d / y(i), 1., y(i) /= 0)
      end do
      sqlp = sqrt(log(pd) * log(py))
      a = n / (log(pd) + sqlp)
      b = n / (log(py) + sqlp)
    end subroutine
  end subroutine

  elemental function loglap_pdf(x, d, a, b) result(p)
    real, intent(in) :: x, d, a, b
    real :: c, p

    p = 0
    if(x <= 0) return

    c = 1. / d * (a * b) / (a + b)
    if(x < d) then
      p = c * (x / d)**(b - 1)
    else
      p = c * (d / x)**(a + 1)
    end if
  end function

  elemental function loglap_ppf(p, d, a, b) result(x)
    real, intent(in) :: p, d, a, b
    real :: ab, x

    x = 0
    if(p <= 0 .or. p >= 1) return

    ab = a / (a + b)
    if(p <= ab) then
      x = d * (p * (a + b) / a)**(1. / b)
    else
      x = d * ((1. - p) * (a + b) / b)**(-1. / a)
    end if
  end function
end module
