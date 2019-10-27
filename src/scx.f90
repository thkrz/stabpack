module scx
  use bez, only: bezcrv
  use intp1d, only: interp
  use mesh, only: mesh_t
  use ssat_env
  implicit none
  private
  protected scxcrk
  protected scxdim
  protected scxnum
  protected scxp_m
  public scxcut
  public scxdel
  public scxini
  public scxmat
  public scxtop
  public scxwa

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

  type(mesh_t) :: scxp_m
  type(stra_t), allocatable :: strata
  real, allocatable, dimension(:, :) :: ridge, scxcrk
  real :: scxdim(2, 2), tana

contains
  ! elemental function scxcrk(x) result(y)
  !   real, intent(in) :: x
  !   real :: y
  !  integer :: i

  !   y = scxtop(x)
  !   i = 1
  !   do while(strata(i)%c == 0)
  !     y = strata(i)%bot(x)
  !     i = i + 1
  !   end do
  !   associate(c => strata(i)%c, w => strata(i)%w,&
  !             phi => strata(i)%phi)
  !     y = y - 2. * c / w * tan(atan(1.) + .4 * phi)
  !   end associate
  ! end function

  pure subroutine scxcut(n, rd, p, w, c, phi, u, alpha, b, h, stat)
    integer, intent(in) :: n
    real, intent(in) :: rd, p(:, :)
    real, intent(out), dimension(n) :: w, c, phi, u, alpha, b
    real, intent(out) :: h(0:n)
    integer, intent(out) :: stat
    real :: d
    real, dimension(2) :: m, q, r
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
      m = .5 * (q + r)
      u(i) = scxp_m%get(m(1), m(2))
      call scxmat(m(1), m(2), c(i), phi(i), w(i))

      q = r
    end do
  end subroutine

  subroutine scxdel
    deallocate(scxcrk)
    deallocate(strata)
    deallocate(ridge)
    scxdim = 0
    scxnum = 0
  end subroutine

  subroutine scxini(name)
    character(*), intent(in) :: name
    character(len=255) :: msg
    real :: alpha, h
    integer :: err, id, i, m, n

    open(newunit=id, file=name, status='old', iostat=err, iomsg=msg)
    if(stat /= 0) call fatal(msg)
    read(id, *, iostat=err, iomsg=msg) alpha
    if(stat /= 0) call fatal(msg)
    tana = tan(radians(alpha))
    read(id, *, iostat=err, iomsg=msg) m, n
    if(stat /= 0) call fatal(msg)
    allocate(ridge(2 + n, m))
    read(id, *, iostat=err, iomsg=msg) (ridge(:, i), i=1,m)
    if(stat /= 0) call fatal(msg)
    read(id, *, iostat=err, iomsg=msg) m
    if(stat /= 0) call fatal(msg)
    scxnum = m
    allocate(strata(m))
    do i = 1, m
      read(id, '10F', iostat=err, iomsg=msg) h, strata(i)
      if(stat /= 0) call fatal(msg)
      strata(i)%x = interp(h, ridge(2, :), ridge(1, :))
      strata(i)%y = scxtop(x)
      strata(i)%d = merge(i, 0, i <= n)
    end do
    read(id, *, iostat=err, iomsg=msg) m
    if(stat /= 0) call fatal(msg)
    allocate(scxcrk(2, m))
    read(id, *, iostat=err, iomsg=msg) (scxcrk(:, i), i=1,m)
    if(stat /= 0) call fatal(msg)

    m = size(ridge, 2)
    scxdim(1, 1) = ridge(1, 1)
    scxdim(1, 2) = minval(ridge(2, :))
    scxdim(2, 1) = ridge(1, m) - scxdim(1, 1)
    scxdim(2, 2) = maxval(ridge(2, :)) - scxdim(1, 2)
  end subroutine

  pure subroutine scxmat(x, y, c, phi, w)
    real, intent(in) :: x, y
    real, intent(out) :: c, phi, w
    real :: h0, h1, l
    integer :: i

    h0 = scxtop(x)
    l = h0 - y
    i = 0
    w = 0
    do while(y < h0 .and. i < scxnum)
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

  pure subroutine scxwa(x, p, k, n, i, r, h)
    real, intent(in) :: x
    real, intent(out) :: p(2, scxnum)
    real, intent(out), dimension(scxnum) :: k, n, i, r, h
    real :: y0, y1
    integer :: j

    y0 = scxtop(x)
    do j = 1, scxnum
      p(:, j) = strata(j)%p
      k(j) = strata(j)%k
      n(j) = strata(j)%n
      i(j) = strata(j)%i
      r(j) = strata(j)%r
      y1 = strata(j)%bot(x)
      h(j) = y0 - y1
      y0 = y1
    end do
  end subroutine

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
