module lsff
  private

  integer, parameter :: iso14688 = 12

  enum, bind(c)
    public
    enumerator :: C=1
    enumerator :: PHI
    enumerator :: RHO_BULK
    enumerator :: RHO_PART
    enumerator :: E
    enumerator :: PR
    enumerator :: THETA
    enumerator :: THETA_RES
    enumerator :: THETA_SAT
    enumerator :: K_SAT
    enumerator :: A
    enumerator :: M
    enumerator :: N
    enumerator :: LAST
  end enum

  type, public :: group
    real :: attr(LAST)
    character(31) :: name
    real, dimension(iso14688) :: scacc
    real, allocatable :: verti(:, :)
  end type
end
