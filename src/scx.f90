module scx
  use ieee_arithmetic, only: ieee_is_finite, ieee_is_nan, ieee_quiet_nan, ieee_value
  use bez, only: bezcrv
  use grid, only: grid_t
  use intp1d, only: interp
  use num_env, only: rad, swap
  implicit none
  private
  public scxcrk
  public scxcut
  public scxdel
  public scxdim
  public scxini
  public scxmat
  public scxnum
  public scxtop
  public scxwa

  type crk_t
    integer i
    real h
  end type

  type stra_t
    real w
    real phi
    real c
    real E
    real nu
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

  logical :: pwpini
  type(grid_t) :: pwp
  type(stra_t), allocatable :: strata(:)
  type(crk_t), allocatable :: scxcrk(:)
  real, allocatable :: omega(:, :)
  real :: scxdim(2, 2), tana
  integer :: scxnum

contains
  elemental function cdepth(x) result(y)
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

  subroutine ykin(d, x, y1, y2)
    real, intent(in) :: d, x(:), y1(:)
    real, intent(out) :: y2(:)
    real :: alpha, drag, h, k(size(y2)), s
    integer :: i, j, n

    n = size(x)
    i = minloc(y1, 1)

    j = maxloc(y1, 1)
    h = y1(j) - y1(i)
    s = sqrt((x(j) - x(i))**2 + h**2)
    alpha = asin(h / s)
    drag = h / (2. * cos(alpha) * s)

    k = kin_(x, y1, drag, 1) + kin_(x, y1, drag, -1)
    k = k / maxval(k)
    !if(k(i) /= 1) call fatal('invalid drag.')
    y2 = y1 - k * d

  contains
    pure function kin_(x, y, crr, dir) result(r)
      real, intent(in) :: x(:), y(:), crr
      integer, intent(in) :: dir
      real :: alpha, h, r(size(x)), s
      integer :: e1, e2, e3, i, ii, j, jj

      e1 = 2
      e2 = size(x)
      e3 = dir
      if(dir < 0) call swap(e1, e2)
      r = 0
      do i = e1, e2, e3
        j = i - 1
        h = dir * (y(j) - y(i))
        s = sqrt((x(i) - x(j))**2 + h**2)
        alpha = asin(h / s)
        ii = i
        jj = j
        if(dir < 0) call swap(ii, jj)
        r(ii) = max(0., h + r(jj) - crr * cos(alpha) * s)
      end do
    end function
  end subroutine

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
      u(i) = merge(pwp%get(m(1), m(2)), 0., pwpini)
      call scxmat(m(1), m(2), c(i), phi(i), w(i))

      q = r
    end do
  end subroutine

  subroutine scxdel
    deallocate(scxcrk)
    deallocate(strata)
    deallocate(omega)
    scxdim = 0
    scxnum = 0
  end subroutine

  subroutine scxini(name, xlim, pf, iomsg, iostat)
    character(*), intent(in) :: name
    real, intent(in) :: xlim(2)
    character(*), intent(in), optional :: pf
    character(*), intent(out) :: iomsg
    integer, intent(out) :: iostat
    character(len=255) :: set
    character(len=1) :: mode
    real :: alpha, h
    integer :: id, a, b, i, m, n, pos
    real, allocatable :: wa(:, :)
    integer, allocatable, dimension(:) :: wa_crk, wa_str

    pwpini = present(pf) .and. len_trim(pf) > 0
    if(pwpini) then
      call pwp%load(pf, iomsg, iostat)
      if(iostat /= 0) return
    end if

    open(newunit=id, file=name, status='old', iostat=iostat, iomsg=iomsg)
    if(iostat /= 0) return

    read(id, *, iostat=iostat, iomsg=iomsg) m
    if(iostat /= 0) return
    allocate(wa(2, m))
    allocate(wa_crk(m))
    allocate(wa_str(m))

    n = 0
    scxnum = 1
    do i = 1, m
      read(id, *, iostat=iostat, iomsg=iomsg) set
      if(iostat /= 0) return
      pos = index(set, 'auto')
      if(pos /= 0) then
        if(i == 1 .or. i == m) call fatal('invalid input file')
        wa(:, i) = ieee_value(1.0, ieee_quiet_nan)
        if(pos == 1) then
          read(set(5:), *) wa(2, i)
        else
          read(set(:pos), *) wa(1, i)
        end if
      else
        read(set, *) wa(:, i)
      end if
      if(index(set, 'CR', .true.) /= 0) then
        n = n + 1
        wa_crk(n) = i
      else if(index(set, 'OC', .true.) /= 0) then
        wa_str(scxnum) = i
        scxnum = scxnum + 1
      end if
    end do

    do i = 2, m - 1
      a = 0
      b = 0
      if(ieee_is_nan(wa(1, i))) then
        a = 1
        b = 2
      else if(ieee_is_nan(wa(2, i))) then
        a = 2
        b = 1
      end if
      if(a + b == 3) wa(a, i) = interp(wa(b, i),&
        [wa(b, i - 1), wa(b, i + 1)],&
        [wa(a, i - 1), wa(a, i + 1)])
    end do

    read(id, *, iostat=iostat, iomsg=iomsg) n
    if(iostat /= 0) return
    allocate(omega(2 + n, m))
    omega(:2, :) = wa
    deallocate(wa)
    scxnum = scxnum + n
    allocate(strata(scxnum))

    do i = 1, n
      read(id, *, iostat=iostat, iomsg=iomsg) mode, h
      if(iostat /= 0) return
      read(id, '(11(F5.2))', iostat=iostat, iomsg=iomsg) strata(i)
      if(iostat /= 0) return
      select case(mode)
        case('C')
          omega(2 + i, :) = omega(2 + i - 1, :) - h
        case('K')
          call ykin(h, omega(1, :), omega(2 + i - 1, :), omega(2 + i, :))
      end select
      strata(i)%d = i
    end do

    read(id, *, iostat=iostat, iomsg=iomsg) alpha
    if(iostat /= 0) return
    tana = tan(rad(alpha))
    do i = 1, scxnum - n
      read(id, '(11(F5.2))', iostat=iostat, iomsg=iomsg) strata(i)
      if(iostat /= 0) return
      if(i < scxnum - n) then
        strata(i + n)%x = omega(1, wa_str(i))
        strata(i + n)%y = omega(2, wa_str(i))
      end if
      strata(i)%d = 0
    end do
    deallocate(wa_str)
    close(id)

    allocate(scxcrk(n))
    do i = 1, n
      scxcrk(i)%i = wa_crk(i)
      scxcrk(i)%h = cdepth(omega(1, wa_crk(i)))
    end do
    deallocate(wa_crk)

    scxdim(1, 1) = merge(xlim(1), omega(1, 1), ieee_is_finite(xlim(1)))
    scxdim(1, 2) = minval(omega(2, :))
    scxdim(2, 1) = merge(xlim(2), omega(1, m), ieee_is_finite(xlim(2))) - scxdim(1, 1)
    scxdim(2, 2) = maxval(omega(2, :)) - scxdim(1, 2)
  end subroutine

  pure subroutine scxmat(x, y, c, phi, w)
    real, intent(in) :: x, y
    real, intent(out) :: c, phi, w
    real :: y0, y1, l
    integer :: i

    y0 = scxtop(x)
    l = y0 - y
    i = 0
    w = 0
    do while(y < y0 .and. i < scxnum)
      i = i + 1
      y1 = strata(i)%bot(x)
      w = w + (y0 - max(y1, y)) / l * strata(i)%w
      y0 = y1
    end do
    c = strata(i)%c
    phi = strata(i)%phi
  end subroutine

  elemental function scxtop(x) result(y)
    real, intent(in) :: x
    real :: y

    y = interp(x, omega(1, :), omega(2, :))
  end function

  pure subroutine scxwa(x, y, beta, k, i, n, r)
    real, intent(in) :: x, y
    real, intent(out) :: beta(2), k, i, n, r
    integer :: j

    j = 1
    do while(y < strata(j)%bot(x) .and. j <= scxnum)
      j = j + 1
    end do
    beta = strata(j)%p
    k = strata(j)%k
    i = strata(j)%i
    n = strata(j)%n
    r = strata(j)%r
  end subroutine

  elemental function stra_t_bot(self, x) result(y)
    class(stra_t), intent(in) :: self
    real, intent(in) :: x
    real :: y

    if(self%d > 0) then
      y = interp(x, omega(1, :), omega(2+self%d, :))
    else
      y = self%y + tana * (x - self%x)
    end if
  end function
end module
