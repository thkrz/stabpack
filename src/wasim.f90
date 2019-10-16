module grndwa
  implicit none
  private

contains
  elemental function piezom(x, dh, l) result(h)
    real, intent(in) :: x, dh, l
    real :: h

    h = dh * sqrt(1. - (x / l)**2)
  end function
end module

module wasim
  use spec, only: simpn
  use swc
  implicit none
  private

  real, allocatable :: p(:)

contains
  pure subroutine seep(n, t0)
    integer, intent(in) :: n
    real, intent(in) :: t0
    integer :: d, e, i, j, k
    real :: eps, k(n), psi(n), t
    
    d = n
    e = n
    eps = 1. / (2. * n)
    do j = 1, k
      i = j - 1
      t = real(i) / n
      if(abs(t - t0) < eps) exit
    end do    
    
    ! dz = 1. / (td - ti) * (K(td)*psi(td)/z(i) + K(td))
  end subroutine

  pure function k(a, b)
    real, intent(in) :: a, b
    real :: k

    k = simpn(qk, a, b) / (b - a)
  contains
    pure function qk(t)
      real, intent(in) :: t
      real :: qk

      qk = swcrhc(t, p)
    end function
  end function

  pure function psi(a, b)
    real, intent(in) :: a, b
    real :: psi

    psi = simpn(qpsi, a, b) / (b - a)
  contains
    pure function qpsi(t)
      real, intent(in) :: t
      real :: qpsi

      qpsi = swcms(t, p)
    end function
  end function
end module
