module mesh
  implicit none
  private
  public mesh_t

  type mesh_t
    integer m
    integer n
    real dx
    real dy
    real x
    real y
    real nan
    real, allocatable :: map(:, :)
  contains
    procedure :: dump => mesh_t_dump
    procedure :: get => mesh_t_get
    procedure :: load => mesh_t_load
    procedure :: init => mesh_t_init
    procedure :: set => mesh_t_set
  end type

contains
  subroutine mesh_t_dump(self, name)
    class(mesh_t), intent(in) :: self
    character(*), intent(in) :: name
    integer :: id, i, j
  end subroutine

  pure function mesh_t_get(self, x, y) result(v)
    class(mesh_t), intent(in) :: self
    real, intent(in) :: x, y
    real :: v
    integer :: i, j

    i = floor(dim(x - self%x) / self%dx)
    j = n - floor(dim(y - self%y) / self%dy)
    if(i > m .or. j < 0) then
      v = 0
    else
      v = self%map(j, i)
    end if
  end function

  subroutine mesh_t_load(self, name)
    class(mesh_t), intent(in) :: self
    character(*), intent(in) :: name
    integer :: id, i, j
  end subroutine

  pure subroutine mesh_t_init(self, extent, dx, dy, nan)
    class(mesh_t), intent(inout) :: self
    real, intent(in) :: extent(2, 2), dx
    real, intent(in), optional :: dy, nan

    self%dx = dx
    self%dy = merge(dy, dx, present(dy))
    self%nan = merge(nan, -999., present(nan))
    self%m = ceilling(extent(2, 1) / self%dx)
    self%n = ceilling(extent(2, 2) / self%dy)
    self%x = extent(1, 1)
    self%y = extent(1, 2)
    allocate(self%map(self%n, self%m))
    self%map = self%nan
  end subroutine

  pure subroutine mesh_t_set(self, x, y, v)
    class(mesh_t), intent(in) :: self
    real, intent(in) :: x, y, v
    integer:: i, j

    i = floor(dim(x - self%x) / self%dx)
    j = n - floor(dim(y - self%y) / self%dy)
    if(i <= m .and j >= 0) self%map(j, i) = v
  end function
end module