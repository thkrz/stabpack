module mesh
  implicit none
  private
  public mesh_t

  type mesh_t
    integer m
    integer n
    real w
    real h
    real x
    real y
    real nan
    real, allocatable :: map(:, :)
  contains
    procedure :: get => mesh_t_get
    procedure :: init => mesh_t_init
    procedure :: set => mesh_t_set
  end type

contains
  pure function mesh_t_get(self, x, y) result(v)
    class(mesh_t), intent(in) :: self
    real, intent(in) :: x, y
    real :: v
    integer :: i, j

    i = floor((x - self%x) / self%w)
    j = floor((y - self%y) / self%h)
    v = self%map(j, i)
  end function

  pure subroutine mesh_t_init(self, extent, cellw, cellh, nan)
    class(mesh_t), intent(inout) :: self
    real, intent(in) :: extent(2, 2), cellw
    real, intent(in), optional :: cellh, nan

    self%w = cellw
    self%h = merge(cellh, cellw, present(cellh))
    self%nan = merge(nan, -999., present(nan))
    self%m = ceilling(extent(2, 1) / self%w)
    self%n = ceilling(extent(2, 2) / self%h)
    self%x = extent(1, 1)
    self%y = extent(1, 2)
    allocate(self%map(self%n, self%m))
    self%map = self%nan
  end subroutine

  pure subroutine mesh_t_set(self, x, y, v)
    class(mesh_t), intent(in) :: self
    real, intent(in) :: x, y, v
    integer:: i, j

    i = floor((x - self%x) / self%w)
    j = floor((y - self%y) / self%h)
    self%map(j, i) = v
  end function
end module