module grndwt
  implicit none
  private
  public piezom

contains
  elemental function piezom(x, dh, l) result(h)
    real, intent(in) :: x, dh, l
    real :: h

    h = dh * sqrt(1. - (x / l)**2)
  end function
end module

module soilwt
  implicit none
  private

! contains
!   pure subroutine hydrus(step, mm, time)
!     real :: dz

!     do i = 1, n
!       dz = 1./ (td - ti) * (K(td)*psi(td)/z(i) + K(td))
!     end do
!   end subroutine
end module
