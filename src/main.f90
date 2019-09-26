program main
  use, intrinsic :: iso_fortran_env, only: error_unit
  use ssat_env, only: fatal
  use slope
  implicit none

  character(len=255) :: arg, input, msg
  integer :: id, stat
  real :: alpha

  namelist /CONFIG/ alpha, input

  if(command_argument_count() /= 1) call fatal('control file missing.')
  call get_command_argument(1, arg)
  open(newunit=id, file=trim(arg), status='old', iostat=stat, iomsg=msg, action='read')
  if(stat /= 0) call fatal(msg)
  read(id, nml=CONFIG, iostat=stat, iomsg=msg)
  if(stat /= 0) call fatal(msg)
  close(id)

  call slopeinit(input, alpha)
  call slopefin
  stop
end program
