module bezier
  implicit none
  private
  public bezier_arc
  public bezier_curve

  integer, parameter :: ndim = 2
  real, parameter :: pi = 4. * atan(1.)

contains
  pure subroutine bezier_arc(p)
    real, parameter :: k = 4. / 3. * (sqrt(2.) - 1.)
    real, intent(inout) :: p(:, :)
    integer :: i, m, n
    real, dimension(2) :: a, b, dx
    real :: r, phi, theta, rot(2, 2)

    n = size(p, 1)
    m = n / 2
    dx = p(n, :) - p(1, :)
    theta = .5 * pi - atan(dx(2) / dx(1))
    rot(1, 1) = cos(theta)
    rot(1, 2) = -sin(theta)
    rot(2, 1) = sin(theta)
    rot(2, 2) = cos(theta)
    a = p(1)
    b = dot_product(rot, dx) + a

    phi = .25 * pi
    r = abs(a(2) - b(2)) / (2. * sin(phi))
    p(m, 1) = a(1) + k * r * sin(phi)
    p(m, 2) = a(2) + k * r * cos(phi)
    p(m + 1, 1) = b(1) + k * r * sin(phi)
    p(m + 1, 2) = b(2) + k * r * cos(phi)
    do i = 2, m - 1
      p(i, :) = p(1, :)
    end do
    do i = m, m + 1
      p(i, :) = dot_product(inv(rot), p(i, :) - a) + a
    end do
    do i = m + 2, n - 1
      p(i, :) = p(n, :)
    end do
  end subroutine

  elemental function bez(t, p)
    real, intent(in) :: t, p(:, :)
    real :: bez(size(p, 1))
    integer :: i, n

    bez = 0
    if(t < 0 .or. t > 1) return
    n = size(p, 2) - 1
    do i = 0, n
      bez = bez + b(t, i, n) * p(:, i+1)
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
