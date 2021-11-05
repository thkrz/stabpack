module scx
  use ieee_arithmetic, only: ieee_is_finite, ieee_is_nan, ieee_quiet_nan, ieee_value
  use bez, only: bezcrv
  use intp1d, only: interp
  implicit none
  private
  public scxcut
  public scxdel
  public scxdim
  public scxini
  public scxmat
  public scxtop

  type stra_t
    real w
    real phi
    real c
    real E
    real nu
    real n
    real i
    real r
    real beta(2)
    real k
  end type

  type(stra_t) :: scxmat
  real, allocatable :: omega(:, :)
  real :: scxdim(2, 2)

contains
  pure subroutine scxcut(n, rd, p, w, c, phi, u, alpha, b, h, stat)
    integer, intent(in) :: n
    real, intent(in) :: rd, p(:, :)
    real, intent(out), dimension(n) :: w, c, phi, u, alpha, b
    real, intent(out) :: h(0:n)
    integer, intent(out) :: stat
    real :: d, z
    real, dimension(2) :: q, r
    integer :: i

    d = rd * norm2(p(:, 1) - p(:, size(p, 2)))
    h(0) = 0
    q = p(:, 1)
    stat = 0
    do i = 1, n
      r = bezcrv(real(i)/n, p)
      if(q(1) > r(1)) then
        stat = -i
        return
      end if

      b(i) = r(1) - q(1)
      h(i) = scxtop(r(1)) - r(2)
      if(h(i) < 0 .or. h(i) > d) then
        stat = i
        return
      end if

      alpha(i) = asin((q(2) - r(2)) / b(i))
      ! m = .5 * (q + r)
      ! z = scxtop(m(1)) - m(2)
      z = .5 * (h(i-1) + h(i))
      u(i) = 0 ! get pwp from grid z

      q = r
    end do

    c = scxmat%c
    phi = scxmat%phi
    w = scxmat%w
  end subroutine

  subroutine scxdel
    deallocate(omega)
    scxdim = 0
  end subroutine

  subroutine scxini(name, xlim, iomsg, iostat)
    character(*), intent(in) :: name
    real, intent(in) :: xlim(2)
    character(*), intent(out) :: iomsg
    integer, intent(out) :: iostat
    character(len=255) :: fnam
    integer :: fid
    namelist /MATERIAL/ scxmat

    open(newunit=fid, file=name, status='old', iostat=iostat, iomsg=iomsg)
    if(iostat /= 0) return
    read(fid, nml=MATERIAL, iostat=iostat, iomsg=iomsg)
    if(iostat /= 0) return
    close(fid)

    scxdim(1, 1) = merge(xlim(1), omega(1, 1), ieee_is_finite(xlim(1)))
    scxdim(1, 2) = minval(omega(2, :))
    scxdim(2, 1) = merge(xlim(2), omega(1, m), ieee_is_finite(xlim(2))) - scxdim(1, 1)
    scxdim(2, 2) = maxval(omega(2, :)) - scxdim(1, 2)
  end subroutine

  elemental function scxtop(x) result(y)
    real, intent(in) :: x
    real :: y

    y = interp(x, omega(1, :), omega(2, :))
  end function
end module
