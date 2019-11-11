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
    integer :: c(size(m)), s(size(x, 2))
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

module rcontd
  implicit none
  private

contains
end module

