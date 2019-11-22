module mesh
  implicit none
  private
  public grid_t

  type grid_t
    integer m
    integer n
    real dx
    real dy
    real x
    real y
    real nan
    real, allocatable :: map(:, :)
  contains
    procedure :: dump => grid_t_dump
    procedure :: get => grid_t_get
    procedure :: load => grid_t_load
    procedure :: init => grid_t_init
    procedure :: set => grid_t_set
  end type

contains
  subroutine grid_t_dump(self, name)
    class(grid_t), intent(in) :: self
    character(*), intent(in) :: name
    integer :: id, i, j
  end subroutine

  pure function grid_t_get(self, x, y) result(v)
    class(grid_t), intent(in) :: self
    real, intent(in) :: x, y
    real :: v
    integer :: i, j

    i = floor(dim(x, self%x) / self%dx)
    j = self%n - floor(dim(y, self%y) / self%dy)
    if(i > self%m .or. j < 0) then
      v = self%nan
    else
      v = self%map(j, i)
    end if
  end function

  subroutine grid_t_load(self, name)
    class(grid_t), intent(in) :: self
    character(*), intent(in) :: name
    integer :: id, i, j
  end subroutine

  pure subroutine grid_t_init(self, extent, dx, dy, nan)
    class(grid_t), intent(inout) :: self
    real, intent(in) :: extent(2, 2), dx
    real, intent(in), optional :: dy, nan

    self%dx = dx
    self%dy = merge(dy, dx, present(dy))
    self%nan = merge(nan, -999., present(nan))
    self%m = ceiling(extent(2, 1) / self%dx)
    self%n = ceiling(extent(2, 2) / self%dy)
    self%x = extent(1, 1)
    self%y = extent(1, 2)
    allocate(self%map(self%n, self%m))
    self%map = self%nan
  end subroutine

  pure subroutine grid_t_set(self, x, y, v)
    class(grid_t), intent(inout) :: self
    real, intent(in) :: x, y, v
    integer:: i, j

    i = floor(dim(x, self%x) / self%dx)
    j = self%n - floor(dim(y, self%y) / self%dy)
    if(i <= self%m .and. j >= 0) self%map(j, i) = v
  end subroutine
end module
