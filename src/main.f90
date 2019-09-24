program main
  use, intrinsic :: iso_fortran_env, only: error_unit
  use ssat_env, only: fatal
  use slope
  implicit none

  character(len=255) :: arg, input, msg
  integer :: id, stat

  namelist /CONFIG/ input

  if(command_argument_count() /= 1) call fatal("control file missing.")
  call get_command_argument(1, arg)
  open(newunit=id, file=trim(arg), status='old', iostat=stat, iomsg=msg, action='read')
  if(stat /= 0) call fatal(msg)
  read(id, nml=CONFIG, iostat=stat, iomsg=msg)
  if(stat /= 0) call fatal(msg)
  close(id)

  call slope_init(input)
  call slope_finalize
  stop
end program
