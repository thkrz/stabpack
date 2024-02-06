module bez
  use math, only: bico, fact, inv2, pi
  implicit none
  private
  public bezarc
  public bezcrv
  public bezfit
  public bezlin
  public bezsim

contains
  pure subroutine bezarc(a, b, p)
    real, intent(in) :: a(:), b(:)
    real, intent(out) :: p(:, :)
    real, parameter :: k = 4. / 3. * (sqrt(2.) - 1.)
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
    do concurrent(i = 1:m-1)
      p(:, i) = a
    end do
    dx =  k * r * [sinb, cosb]
    call inv2(rot)
    p(:, m) = matmul(rot, dx) + a
    p(:, m+1) = matmul(rot, br + dx - a) + a
    do concurrent(i = m+2:n)
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
      a(i+1) = norm2(x(:, i+1) - bezcrv(t, p))
    end do
    if (all(a < tol)) return

    k = maxloc(a, 1)
    if (k < 3 .or. n - k < 2) return

    q = bezfit(x(:, :k), tol)
    r = bezfit(x(:, k:), tol)
    n = size(q, 1)
    k = size(r, 1)
    deallocate(p)
    allocate(p(n+k-1, 2))
    p(:, :n) = q
    p(:, n:) = r
    deallocate(q)
    deallocate(r)
  end function

  pure subroutine bezlin(a, b, p)
    real, intent(in) :: a(:), b(:)
    real, intent(out) :: p(:, :)
    real :: dx(size(a))
    integer :: i, n

    dx = b - a
    n = size(p, 2) - 1
    do concurrent(i=0:n)
      p(:, i+1) = dx * real(i) / n + a
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
      lambda = lambda + hypot(p(:, i) - q(:, i))
    end do
    lambda = lambda / n
  end function
end module
