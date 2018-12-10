module bezier
contains
  pure subroutine bezarc(h, p)
    real, parameter :: k = 4. / 3. * (sqrt(2.) - 1.)
    real, intent(in) :: h
    real, intent(inout) :: p(:, :)
    integer :: i, n
    real :: a(2), b(2), c, r, phi
  
    n = size(p, 2)
    c = norm2(p(:, 1) - p(:, n))
    r = h * .5 + c**2 / (8. * h)
    phi = acos(1. - h / r)
    a(1) = p(1, 1) + k * r * sin(phi)
    a(2) = p(2, 1) - k * r * cos(phi)
    b(1) = p(1, n) + k * r * sin(phi)
    b(2) = p(2, n) - k * r * cos(phi)
    do concurrent(i = 2:n-1)
      p(:, i) = merge(b, a, i > n/2)
    end do
  end subroutine

  pure function bez(p, t)
    real, intent(in) :: p(:, :), t
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
end module
