module critss
  use bez
  use fmin, only: amoeba
  use razdol, only: razslv
  use slope
  implicit none
  private

contains
  pure subroutine optim(m, n, a, b)
    integer, intent(in) :: m, n
    real, intent(in) :: a(:), b(:)
    real :: p(n, size(a)), t(m + 1)
    real :: mu
    integer :: i

    p(1, :) = a
    p(n, :) = b
    beza(p)
    do concurrent(i = 0, m)
      t(i + 1) = real(i) / m
    end do
    !reshape
    amoeba(f, p(2:n-1, :), fn=mu)
  contains
    pure function f(x)
      real, intent(in) :: x
      real :: r(m + 1), xy(m + 1, 2)
      integer :: i

      p(2:n-1, :) = x
      do concurrent(i = 1, m + 1)
        xy(i, :) = bezc(t(i), p)
      end do
      r = xy(:m, :) - xy(
    end function
  end subroutine
end module
