module rtfind
  implicit none
  private
  public rtnewt

  integer, parameter :: itmax = 1000

contains
  pure subroutine rtnewt(funcd, p, x0, x, tol, stat)
    interface
      pure subroutine funcd(x, p, fval, fderiv)
        real, intent(in) :: x, p(:)
        real, intent(out) :: fval, fderiv
      end subroutine
    end interface
    real, intent(in) :: p(:), x0
    real, intent(out) :: x
    real, optional, intent(in) :: tol
    integer, optional, intent(out) :: stat
    real :: df, dx, f, xacc
    integer :: j

    if(present(stat)) stat = 0
    xacc = merge(tol, sqrt(epsilon(1.)), present(tol))
    x = x0
    do j = 1, itmax
      call funcd(x, p, f, df)
      dx = f / df
      x = x - dx
      if(abs(dx) < xacc) return
    end do
    if(present(stat)) stat = -1
  end subroutine
end module
