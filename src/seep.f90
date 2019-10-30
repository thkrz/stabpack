program main
  use mesh, only: mesh_t
  use ssat_env, only: fatal
  use scx
  implicit none

  character(len=255) :: arg, datafile, msg, precfile
  integer :: bins = 100, err, id
  real :: step = 1., xlim(2) = 0

  namelist /CONFIG/ bins, datafile, precfile, step, xlim

  if(command_argument_count() /= 1) call fatal('control file missing.')
  call get_command_argument(1, arg)
  open(newunit=id, file=trim(arg), status='old', iostat=err, iomsg=msg, action='read')
  if(err /= 0) call fatal(msg)
  read(id, nml=CONFIG, iostat=err, iomsg=msg)
  if(err /= 0) call fatal(msg)
  close(id)

  call scxini(datafile, xlim)
  call scxdel
  stop

contains
end program
