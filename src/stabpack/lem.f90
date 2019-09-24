module leq
  private

contains
  pure subroutine razdol
  end subroutine

  pure subroutine spence(n, alpha, b, theta, c, phi, w, u, fos)
    integer, intent(in) :: n
    real, dimension(n), intent(in) :: alpha, b, c, phi, w, u
    real, intent(in) :: theta
    real, intent(out) :: fos
    real, dimension(n+1) :: h, z
    real :: f
  contains
    pure subroutine force(x, fval, fderiv)
      real, intent(in) :: x
      real, intent(out) :: fval, fderiv
      integer :: i, j

      do i = 1, n
        z(i+1) = z(i) * alpha(i)
      end do
      fval = z(n+1)

      fderiv = z(1)
      do i = 1, n
        fderiv = fderiv * i
      end do
    end subroutine

  end subroutine
end module
