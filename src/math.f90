module math
  implicit none
  private
  public pi
  public bico
  public deg
  public fact
  public inv2
  public rad
  public shft
  public swap

  real, parameter :: pi = 4. * atan(1.)

  interface swap
    procedure iswap, rswap
  end interface

contains
  pure function bico(n, k)
    integer, intent(in) :: n, k
    integer :: bico

    bico = fact(n) / (fact(k) * fact(n-k))
  end function

  elemental function deg(r) result(d)
    real, intent(in) :: r
    real :: d

    d = r * 180. / pi
  end function

  pure function fact(n)
    integer, intent(in) :: n
    integer :: i, fact

    fact = product([(i, i=1,n)])
  end function

  pure subroutine inv2(a)
    real, intent(inout) :: a(2, 2)
    real :: b(2, 2)
    real :: detinv

    detinv = 1. / (a(1, 1) * a(2, 2) - a(1, 2) * a(2, 1))

    b(1, 1) = +detinv * a(2, 2)
    b(2, 1) = -detinv * a(2, 1)
    b(1, 2) = -detinv * a(1, 2)
    b(2, 2) = +detinv * a(1, 1)
    a = b
  end subroutine

  elemental function rad(d) result(r)
    real, intent(in) :: d
    real :: r

    r = d * pi / 180.
  end function

  pure subroutine shft(a, b, c, d)
    real, intent(out) :: a
    real, intent(inout) :: b, c
    real, intent(in) :: d

    a = b
    b = c
    c = d
  end subroutine

  pure subroutine iswap(a, b)
    integer, intent(inout) :: a, b
    integer :: c

    c = a
    a = b
    b = c
  end subroutine

  pure subroutine rswap(a, b)
    real, intent(inout) :: a, b
    real :: c

    c = a
    a = b
    b = c
  end subroutine
end module
