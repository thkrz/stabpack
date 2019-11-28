module grid
  use, intrinsic :: ieee_arithmetic
  use netcdf
  implicit none
  private
  public grid_t

  type grid_t
    private
    integer m
    integer n
    real dx
    real dy
    real x
    real y
    real, allocatable :: map(:, :)
  contains
    procedure :: dump => grid_t_dump
    procedure :: get => grid_t_get
    procedure :: load => grid_t_load
    procedure :: init => grid_t_init
    procedure :: set => grid_t_set
  end type

contains
  subroutine grid_t_dump(self, name, iomsg, iostat)
    class(grid_t), intent(in) :: self
    character(*), intent(in) :: name
    character(*), intent(out) :: iomsg
    integer, intent(out) :: iostat
    integer :: dimids(2), ncid, varid, void

    iostat = 0
    ncid = -1
    if(err(nf90_create(name, ior(NF90_NETCDF4, NF90_CLASSIC_MODEL), ncid))) return
    if(err(nf90_def_dim(ncid, 'm', self%m, dimids(1)))) return
    if(err(nf90_def_dim(ncid, 'n', self%n, dimids(2)))) return
    if(err(nf90_def_var(ncid, 'map', NF90_DOUBLE, dimids, varid))) return
    if(err(nf90_put_att(ncid, varid, 'dx', self%dx))) return
    if(err(nf90_put_att(ncid, varid, 'dy', self%dy))) return
    if(err(nf90_put_att(ncid, varid, 'origin_x', self%x))) return
    if(err(nf90_put_att(ncid, varid, 'origin_y', self%y))) return
    if(err(nf90_enddef(ncid))) return
    if(err(nf90_put_var(ncid, varid, self%map))) return
    void = nf90_close(ncid)

  contains
    function err(stat)
      integer, intent(in) :: stat
      logical :: err

      err = stat /= NF90_NOERR
      if(err) then
        if(ncid /= -1) void = nf90_close(ncid)
        iostat = -1
        iomsg = trim(nf90_strerror(stat))
      end if
    end function
  end subroutine

  pure function grid_t_get(self, x, y) result(v)
    class(grid_t), intent(in) :: self
    real, intent(in) :: x, y
    real :: v
    integer :: i, j

    i = floor(dim(x, self%x) / self%dx)
    j = self%n - floor(dim(y, self%y) / self%dy)
    if(i > self%m .or. j < 0) then
      v = ieee_value(v, ieee_quiet_nan)
    else
      v = self%map(j, i)
    end if
  end function

  subroutine grid_t_load(self, name, iomsg, iostat)
    class(grid_t), intent(inout) :: self
    character(*), intent(in) :: name
    character(*), intent(out) :: iomsg
    integer, intent(out) :: iostat
    integer :: dimid, ncid, varid, void

    iostat = 0
    ncid = -1
    if(err(nf90_open(name, NF90_NOWRITE, ncid))) return
    if(err(nf90_inq_dimid(ncid, 'm', dimid))) return
    if(err(nf90_inquire_dimension(ncid, dimid, len=self%m))) return
    if(err(nf90_inq_dimid(ncid, 'n', dimid))) return
    if(err(nf90_inquire_dimension(ncid, dimid, len=self%n))) return
    allocate(self%map(self%n, self%m))
    if(err(nf90_inq_varid(ncid, 'map', varid))) return
    if(err(nf90_get_att(ncid, varid, 'dx', self%dx))) return
    if(err(nf90_get_att(ncid, varid, 'dy', self%dy))) return
    if(err(nf90_get_att(ncid, varid, 'origin_x', self%x))) return
    if(err(nf90_get_att(ncid, varid, 'origin_y', self%y))) return
    if(err(nf90_get_var(ncid, varid, self%map))) return
    void = nf90_close(ncid)

  contains
    function err(stat)
      integer, intent(in) :: stat
      logical :: err

      err = stat /= NF90_NOERR
      if(err) then
        if(allocated(self%map)) deallocate(self%map)
        if(ncid /= -1) void = nf90_close(ncid)
        iostat = -1
        iomsg = trim(nf90_strerror(stat))
      end if
    end function
  end subroutine

  pure subroutine grid_t_init(self, extent, dx, dy)
    class(grid_t), intent(inout) :: self
    real, intent(in) :: extent(2, 2), dx
    real, intent(in), optional :: dy
    real :: nan

    self%dx = dx
    self%dy = merge(dy, dx, present(dy))
    self%m = ceiling(extent(2, 1) / self%dx)
    self%n = ceiling(extent(2, 2) / self%dy)
    self%x = extent(1, 1)
    self%y = extent(1, 2)
    allocate(self%map(self%n, self%m))
    nan = ieee_value(nan, ieee_quiet_nan)
    self%map = nan
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
