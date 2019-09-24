module bezier
  use ssat_env, only: ndim, pi
  implicit none
  private
  public bezier_arc
  public bezier_curve
  public bezier_line

contains
  pure subroutine bezier_arc(p)
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

  pure function bezier_curve(t, p) result(c)
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

  pure subroutine bezier_line(p)
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
