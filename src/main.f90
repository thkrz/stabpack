program main
  use bez, only: bezarc, bezlin
  use fmin, only: amoeba
  use razdol, only: razslv
  use ssat_env, only: fatal
  use scx
  implicit none

  type res_t
    real :: mu
    real, allocatable :: p(:, :)
    type(res_t), pointer :: next => null()
  end type

  character(len=255) :: arg, datafile, msg, precipitation
  character(len=4) :: mode = 'stab'
  integer :: bezn = 3, i, id, m, n, num = 100, stat
  real :: a(2), b(2), dx = 5., fos = 1.3, p0, rd = 0.3, xlim(2) = (/ -huge(1.), huge(1.) /)
  type(res_t) :: result

  namelist /CONFIG/ bezn, datafile, dx, fos, mode, num, precipitation, rd, xlim(2)

  if(command_argument_count() /= 1) call fatal('control file missing.')
  call get_command_argument(1, arg)
  open(newunit=id, file=trim(arg), status='old', iostat=stat, iomsg=msg, action='read')
  if(stat /= 0) call fatal(msg)
  read(id, nml=CONFIG, iostat=stat, iomsg=msg)
  if(stat /= 0) call fatal(msg)
  close(id)

  call scxini(datafile)
  m = size(scxcrk, 2)
  n = ceiling(scxdim(2, 1) / dx)
  !$omp parallel do private(i, j, a, b, p0)
  do i = 1, m
    a = scxcrk(:, i)
    p0 = hystp(scxtop(a(1)) - a(2))
    b(1) = a(1)
    do j = i+1, n
      b(1) = b(1) + dx
      b(2) = scxtop(b(1))
      call stab(a, b, p0)
    end do
  end do
  !$omp parallel do private(i, j, a, b)
  do i = 0, n-1
    a(1) = i * dx
    a(2) = scxtop(a(1))
    do j = i+1, n
      b(1) = j * dx
      b(2) = scxtop(b(2))
      call stab(a, b, 0.)
    end do
  end do
  call scxdel
  stop

contains
  subroutine stab(a, b, p0)
    real, intent(in) :: a(2), b(2), p0
    real :: h(0:num), p(2, bezn + 1)
    real, dimension(num) :: w, c, phi, u, alpha, b
    integer :: err

    call bezarc(a, b, p)
    call calc
    call bezlin(a, b, p)
    call calc
  contains
    subroutine calc
      real :: mu, popt(2 * bezn - 2)
      integer :: err

      popt = pack(p(:, 2:bezn), .true.)
      call amoeba(mos, popt, fn=mu, stat=err)
      if(err == 0 .and. mu < 1) then
        p(:, 2:bezn) = reshape(popt, (/ 2, bezn - 1 /))
        call resadd(mu, p)
      end if
    end subroutine

    pure function mos(x)
      real, intent(in) :: x(:)
      real :: pn(2, bezn + 1)
      real, dimension(0:num) :: e, t
      inetger :: err

      pn(:, 1) = p(:, 1)
      pn(:, 2:bezn) = reshape(x, (/ 2, bezn - 1 /))
      pn(:, bezn + 1) = p(:, bezn + 1)
      scxcut(num, rd, pn, w, c, phi, u, alpha, b, h, err)
      if(err /= 0) then
        mos = 10.
        return
      end if
      call razdslv(num, w, c, phi, u, alpha, b, h, fos, e, t, mu, p0)
      mos = mu
    end function
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
