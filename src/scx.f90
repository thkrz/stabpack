module scx
  use bez, only: bezc
  use intp1d, only: interp
  use mesh, only: mesh_t
  implicit none
  private
  protected scxdim
  protected scxp_m
  public scxcrk
  public scxcut
  public scxdel
  public scxini
  public scxmat
  public scxtop
  public scxwas

  type stra_t
    real w
    real phi
    real c
    real n
    real i
    real r
    real p(2)
    real k
    real x
    real y
    integer d
  contains
    procedure :: bot => stra_t_bot
  end type

  type(mesh_t) :: scxf_m, scxp_m
  type(stra_t), allocatable :: strata
  real, allocatable :: ridge(:, :)
  real :: scxdim(2), tana

contains
  elemental function scxcrk(x) result(y)
    real, intent(in) :: x
    real :: y
    integer :: i

    y = scxtop(x)
    i = 1
    do while(strata(i)%c == 0)
      y = strata(i)%bot(x)
      i = i + 1
    end do
    associate(c => strata(i)%c, w => strata(i)%w,&
              phi => strata(i)%phi)
      y = y - 2. * c / w * tan(atan(1.) + .4 * phi)
    end associate
  end function

  pure subroutine scxcut(n, rd, p, w, c, phi, u, alpha, b, h, stat)
    integer, intent(in) :: n
    real, intent(in) :: rd, p(:, :)
    real, intent(out), dimension(n) :: w, c, phi, u, alpha, b
    real, intent(out) :: h(0:n)
    integer, intent(out) :: stat
    real :: d
    real, dimension(2) :: m, q, r
    integer :: i

    d = rd * norm2(p(1, :) - p(size(p, 1), :))
    h(0) = 0
    q = p(1, :)
    stat = 0
    do i = 1, n
      r = bezc(real(i)/n, p)
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
      m = .5 * (q + r)
      u(i) = scxp_m%get(m(1), m(2))
      call scxmat(m(1), m(2), c(i), phi(i), w(i))

      q = r
    end do
  end subroutine

  subroutine scxdel
  end subroutine

  subroutine scxini(name)
    character(*), intent(in) :: name
    integer :: n

    n = size(ridge, 2)
    scxdim = (/ ridge(1, n) - ridge(1, 1), &
                maxval(ridge(2, :)) - minval(ridge(2, :)) /)
  end subroutine

  pure subroutine scxmat(x, y, c, phi, w)
    real, intent(in) :: x, y
    real, intent(out) :: c, phi, w
    real :: h0, h1, l
    integer :: i, n

    n = size(strata)
    h0 = scxtop(x)
    l = h0 - y
    i = 0
    w = 0
    do while(y < h0 .and. i < n)
      i = i + 1
      h1 = strata(i)%bot(x)
      w = w + (h0 - max(h1, y)) / l * strata(i)%w
      h0 = h1
    end do
    c = strata(i)%c
    phi = strata(i)%phi
  end subroutine

  elemental function scxtop(x) result(y)
    real, intent(in) :: x
    real :: y

    y = interp(x, ridge(1, :), ridge(2, :))
  end function

  elemental function stra_t_bot(self, x) result(y)
    class(stra_t), intent(in) :: self
    real, intent(in) :: x
    real :: y

    if(self%d > 0) then
      y = interp(x, ridge(1, :), ridge(2+self%d, :))
    else
      y = self%y + tana * (x - self%x)
    end if
  end function
end module
