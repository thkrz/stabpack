program main
  use ssat_env, only: fatal
  use scx
  implicit none

  type res_t
    real :: mu
    real, allocatable :: p(:, :)
    type(res_t), pointer :: next => null()
  end type

  character(len=255) :: arg, datafile, msg, precipitation
  character(len=4) :: mode = 'crit'
  integer :: bezn, i, id, num, stat, n
  real :: a(2), b(2), fos, dx, rd
  type(res_t) :: result

  namelist /CONFIG/ bezn, datafile, dx, fos, mode, num, precipitation, rd

  if(command_argument_count() /= 1) call fatal('control file missing.')
  call get_command_argument(1, arg)
  open(newunit=id, file=trim(arg), status='old', iostat=stat, iomsg=msg, action='read')
  if(stat /= 0) call fatal(msg)
  read(id, nml=CONFIG, iostat=stat, iomsg=msg)
  if(stat /= 0) call fatal(msg)
  close(id)

  call scxini(datafile)
  n = ceiling(scxdim(1) / dx)
  !$omp do private(i, j, a, b, x)
  do i = 0, n - 1
    x = i * dx
    a = (/ x, scxtop(x) /)
    do j = i, n
      x = x + dx
      b = (/ x, scxtop(x) /)
      call crstab(a, b)
    end do
  end do
  call scxdel
  stop

contains
  subroutine crstab(a, b)
    real, intent(in) :: a(2), b(2)
    real :: h(0:num), p(2, bezn + 1)
    real, dimension(num) :: w, c , phi, u, alpha, b
    integer :: stat

    call beza(a, b, p)
    call scxcut(num, rd, p, w, c, phi, u, alpha, b, h, stat)
    if(stat == 0) call crit(p)
    call bezl(a, b, p)
    call scxcut(num, rd, p, w, c, phi, u, alpha, b, h, stat)
    if(stat == 0) call crit(p)
  contains
    subroutine crit
      real :: mu
      integer :: stat

      if(stat == 0 .and. mu < 1) call resadd(mu, p)
    end subroutine
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
