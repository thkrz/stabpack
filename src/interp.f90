module interp1d
  implicit none
  private
  public interp

contains
  pure function interp(x, xp, fp) result(y)
    real, intent(in) :: x, xp(:), fp(:)
    integer :: i, j, n

    n = size(xp)
    if(x < xp(1)) then
      y = fp(1)
      return
    end if
    do i = 1, n - 1
      j = i + 1
      if(x >= xp(i) .and. x < xp(j)) then
        y = fp(i) + (x - xp(i)) * (fp(j) - fp(i)) / (xp(j) - xp(i))
        return
      end if
    end do
    if(x > xp(n)) y = fp(n)
  end function
end module
