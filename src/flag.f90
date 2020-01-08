module flag
  implicit none
  private
  public flgget
  public flgopt

contains
  subroutine flgget(param, c, rep)
    character(*), intent(inout) :: param
    character, intent(out) :: c
    logical, intent(in), optional :: rep
    character(len=255) :: arg
    integer, save :: i = 1

    if(present(rep) .and. rep .eqv. .true.) then
      c = param(1:1)
      param = param(2:)
      return
    end if
    c = ''
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
  end subroutine

  subroutine flgopt(opt, stat)
    character(*), intent(out) :: opt
    integer, intent(out) :: stat
    integer, save :: i = 1

    stat = -1
    if(i > command_argument_count()) return
    call get_command_argument(i, opt)
    do while(opt(:1) == '-')
      i = i + 1
      if(i > command_argument_count()) return
      call get_command_argument(i, opt)
    end do
    i = i + 1
    stat = 0
  end subroutine
end module
