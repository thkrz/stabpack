module flag
  implicit none
  private
  public flgget
  public flgopt

contains
  logical function flgget(param, c, rep)
    character(*), intent(inout) :: param
    character, intent(out) :: c
    logical, intent(in), optional :: rep
    character(len=255) :: arg
    integer, save :: i = 1

    if(present(rep) .and. rep .eqv. .true.) then
      flgget = len_trim(param) < 2
      c = param(1:1)
      param = param(2:)
      return
      ! never reached
    end if
    flgget = .false.
    if(i > command_argument_count()) return
    call get_command_argument(i, arg)
    do while(arg(:1) /= '-')
      i = i + 1
      if(i > command_argument_count()) return
      call get_command_argument(i, arg)
    end do
    c = arg(2:2)
    param = arg(3:)
    i = i + 1
    flgget = .true.
  end function

  logical function flgopt(opt)
    character(*), intent(out) :: opt
    character(len=255) :: arg
    integer, save :: i = 1

    flgopt = .false.
    if(i > command_argument_count()) return
    call get_command_argument(i, arg)
    do while(arg(:1) == '-')
      i = i + 1
      if(i > command_argument_count()) return
      call get_command_argument(i, arg)
    end do
    i = i + 1
    opt = arg
    flgopt = .true.
  end function
end module
