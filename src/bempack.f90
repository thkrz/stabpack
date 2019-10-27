module mesh
  implicit none
  private
  public mesh_t

  type mesh_t
    integer m
    integer n
    real w
    real h
    real xll
    real yll
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


  end function

  pure subroutine mesh_t_init(self, m, n, w, h, x0, y0, nan)
    class(mesh_t), intent(inout) :: self
    integer, intent(in) :: m, n
    real, intent(in) :: w, h, x0, y0
    real, intent(in), optional :: nan

    allocate(self%map(n, m))
    self%map = merge(nan, 0, present(nan))
    self%m = m
    self%n = n
    self%w = w
    self%h = h
    self%xll = x0
    self%yll = y0
  end subroutine

  pure subroutine mesh_t_set(self, x, y, v)
    class(mesh_t), intent(in) :: self
    real, intent(in) :: x, y, v
  end function
end module