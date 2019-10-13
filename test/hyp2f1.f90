program hyp2f1_test
  use spec, only: hyp2f1
  implicit none

  character(len=255) :: arg
  real :: a, b, c, z

  call get_command_argument(1, arg)
  read(arg, *) a
  call get_command_argument(2, arg)
  read(arg, *) b
  call get_command_argument(3, arg)
  read(arg, *) c
  call get_command_argument(4, arg)
  read(arg, *) z
  print *, hyp2f1(a, b, c, z)
end program
