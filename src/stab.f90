module stab
  use bez, only: bezarc, bezlin
  use fmin, only: amoeba
  use razdol, only: razslv
  use scx, only: scxcut
  implicit none
  private
  public incrit

  real, parameter :: maxmu = 10.

  subroutine incrit(num, q, r, p0, mode, mu, p, stat)
    integer, intent(in) :: num
    real, intent(in) :: q(2), r(2), p0
    character, intent(in) :: mode
    real, intent(out) :: mu, p(:, :)
    integer, intent(out) :: stat
    real, dimension(num) :: w, c, phi, u, alpha, b
    real, dimension(0:num) :: e, h, t
    real :: popt(2 * size(p, 2) - 4)
    integer :: bezn

    bezn = size(p, 2) - 1
    stat = 0
    if(mode == 'A') then
      call bezarc(q, r, p)
    else if(mode == 'L')
      call bezlin(q, r, p)
    else
      stat = -1
      return
    end if
    popt = pack(p(:, 2:bezn), .true.)
    call amoeba(mos, popt, fn=mu, stat=stat)
    if(mu == maxmu) stat = -1
  contains
    function mos(x)
      real, intent(in) :: x(:)
      integer :: err

      p(:, 2:bezn) = reshape(x, (/ 2, bezn - 1 /))
      scxcut(num, rd, p, w, c, phi, u, alpha, b, h, err)
      if(err /= 0) then
        mos = maxmu
        return
      end if
      call razslv(num, w, c, phi, u, alpha, b, h, fos, e, t, mu, p0)
      mos = mu
    end function
  end subroutine
end module

program main
  use ieee_arithmetic
  use ssat_env, only: fatal
  use scx
  implicit none

  type res_t
    real :: mu
    real, allocatable :: p(:, :)
    type(res_t), pointer :: next => null()
  end type

  character(len=255) :: arg, datafile, usage
  integer :: bezn, i, id, m, n, num
  real :: a(2), b(2), bound, dx, fos, rd, xlim(2)
  type(res_t) :: result

  usage = 'usage: stab [-ddx] [-ffos] [-n[B]num] [-rratio] [-x0|1value] file'
  bezn = 3
  dx = 5.
  fos = 1.3
  rd = .3
  xlim(1) = ieee_value(xlim(1), ieee_negative_inf)
  xlim(2) = ieee_value(xlim(2), ieee_positive_inf)
  do i = 1, command_argument_count()
    call get_command_argument(i, arg)
    if(arg(:1) == '-') then
      select case(arg(2:2))
      case('d')
        read(arg(3:), *) dx
      case('f')
        read(arg(3:), *) fos
      case('n')
        if(arg(3:3) == 'B') then
          read(arg(4:), *) bezn
        else
          read(arg(3:), *) num
        end if
      case('r')
        read(arg(3:), *) rd
      case('x')
        if(arg(3:3) == '0') then
          read(arg(4:), *) xlim(1)
        else if(arg(3:3) == '1') then
          read(arg(4:), *) xlim(2)
        else
          call fatal(usage)
        end if
      case default
        call fatal(usage)
      end select
    else
      datafile = arg
    end if
  end do
  if(len_trim(datafile) == 0) call fatal(usage)

  call scxini(datafile, xlim)
  bound = scxdim(1, 1) + scxdim(2, 1)
  m = size(scxcrk, 2)
  !$omp parallel do private(i, a, b)
  do i = 1, m
    a = scxcrk(:2, i)
    b(1) = a(1) + dx
    do while(b(1) < bound)
      b(2) = scxtop(b(1))
      call inloc(a, b, scxcrk(3, i))
      b(1) = b(1) + dx
    end do
  end do
  !$omp end parallel do

  n = ceiling(scxdim(2, 1) / dx)
  !$omp parallel do private(i, j, a, b)
  do i = 0, n-1
    a(1) = scxdim(1, 1) + i * dx
    a(2) = scxtop(a(1))
    do j = i+1, n
      b(1) = scxdim(1, 1) + j * dx
      b(2) = scxtop(b(1))
      call inloc(a, b, 0.)
    end do
  end do
  !$omp end parallel do
  call scxdel
  stop

contains
  subroutine inloc(a, b, p0)
    real, intent(in) :: a(2), b(2), p0
    use bez, only: bezarc, bezlin
    use stab, only: incrit
    real :: mu, p(2, bezn + 1)
    integer :: err

    call bezarc(a, b, p)
    call incrit(num, a, b, p0, mu, p, err)
    if(err == 0) call resadd(mu, p)

    call bezlin(a, b, p)
    call incrit(num, a, b, p0, mu, p, err)
    if(err == 0) call resadd(mu, p)
  end subroutine

  subroutine resadd(mu, p)
    real, intent(in) :: mu, p(:, :)
    type(res_t) :: i
    type(res_t), pointer :: n, r, temp

    i%mu = mu
    allocate(i%p, source=p)
    allocate(n, source=i)

    !$omp critical add
    if(.not. assoicated(result)) then
      result => n
      return
    end if
    if(mu < result%mu) then
      temp => result%next
      result => n
      result%next => temp
      return
    end if
    r => result
    do while(associated(r%next) .and. mu > r%next%mu)
      r => r%next
    end do
    temp => r%next
    r%next => n
    r%next%next => temp
    !$omp end critical add
  end subroutine
end program
