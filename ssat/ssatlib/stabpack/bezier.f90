module bezier
  implicit none

contains
  pure function arc(k, n, a, b, h) result(c)
    real(8), parameter :: kappa = 4.d0 / 3.d0 * (sqrt(2.d0) - 1.d0)
    integer, intent(in) :: k
    integer, intent(in) :: n
    real(8), intent(in) :: a(n), b(n)
    real(8), intent(in) :: h
    real(8) :: c(k+1, n)
    real(8) :: d, r, phi
    integer :: i

    d = norm2(a - b)
    r = h * 0.d05 + d**2 / (8.d0 * h)
    phi = acos(1.d0 - h / r)
    i = (k + 1) / 2 + 1
    c(:i, 1) = a(1)
    c(:i, 2) = a(2)
    c(i:, 1) = b(1)
    c(i:, 2) = b(2)
    c(2:k, 1) = c(2:k, 2) + kappa * r * cos(phi)
    c(2:k, 2) = c(2:k, 2) - kappa * r * cos(phi)
  end function

  pure function bez(p, t)
    real(8), intent(in) :: p(:, :), t
    real(8) :: bez(size(p, 1))
    integer :: i, n

    bez = 0
    if(t < 0 .or. t > 1) return
    n = size(p, 2)
    do i = 0, n
      bez = bez + b(t, i, n) * p(:, i+1)
    end do

  contains
    pure function b(t, i, n)
      real(8), intent(in) :: t
      integer, intent(in) :: i, n
      real(8) :: b

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
end module
