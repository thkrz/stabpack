module clusan
  implicit none
  private
  public kmeans

  integer, parameter :: ITMAX = 10000

contains
  pure subroutine kmeans(x, m)
    real, intent(in) :: x(:, :)
    real, intent(inout) :: m(:, :)
    real :: c(size(m, 1), size(m, 2)), tol, v, vmin
    logical :: mask(size(x, 2))
    integer :: s(size(x, 2))
    integer :: i, j, k, num, n, z

    tol = sqrt(epsilon(1.))
    k = size(m, 2)
    n = size(x, 2)
    do z = 1, ITMAX
      do j = 1, n
        vmin = -1.
        do i = 1, k
          v = norm2(x(:, j) - m(:, i))
          if(vmin > v .or. vmin < 0) then
            vmin = v
            s(j) = i
          end if
        end do
      end do
      do i = 1, k
        mask = s == i
        num = count(mask)
        do concurrent(j = 1:size(m, 1))
          c(j, i) = sum(x(j, :), mask) / num
        end do
      end do
      if(norm2(c - m) < tol) exit
      m = c
    end do
  end subroutine
end module
