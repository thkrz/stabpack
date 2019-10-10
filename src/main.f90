module crit
  implicit none
  private
end module

module pest
  implicit none
  private
end module

program main
  use ssat_env, only: fatal
  use scx
  implicit none

  character(len=255) :: arg, input, msg
  character(len=4) :: mode
  integer :: id, stat
  real :: rd

  namelist /CONFIG/ input, mode, rd

  if(command_argument_count() /= 1) call fatal('control file missing.')
  call get_command_argument(1, arg)
  open(newunit=id, file=trim(arg), status='old', iostat=stat, iomsg=msg, action='read')
  if(stat /= 0) call fatal(msg)
  read(id, nml=CONFIG, iostat=stat, iomsg=msg)
  if(stat /= 0) call fatal(msg)
  close(id)

  call scxnew(input)
  if(mode == 'crit') call fndsls
  if(mode == 'pest') call fndpar
  call scxdel
  stop

contains
  subroutine fndsls
    use crit
  end subroutine

  subroutine fndpar
    use pest
  end subroutine
end program
