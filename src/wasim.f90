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
  pure subroutine seep(n, t0, dt)
    integer, intent(in) :: n
    real, intent(in) :: t0
    integer :: d, i, j
    real :: a, b, eps, k(n), psi(n), t, z(n), zd
   
    eps = 1. / (2. * n)
    do j = 1, n
      i = j - 1
      t = real(i) / n
      if(abs(t - t0) < eps) exit
    end do
    
    do concurrent(j = i:n-1)
      d = j + 1
      a = real(i) / n
      b = real(j) / n
      k(j) = simpn(k2, a, b)
      psi(j) = simpn(psi2, a, b)
    end do
    
    d = n
    z = zd
    do
      do j = i+1, d
        
      end do
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
