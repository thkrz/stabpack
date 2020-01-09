module grndwa
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

! module wasim
!   use intgrt, only: simpr
!   use rtfind, only: rtnewt
!   use vg
!   implicit none
!   private
!   public waga
!   public waseep

!   type(vg_t) :: swc

! contains
!   pure function waga(k, psi, t, delta) result(f)
!     real, intent(in) :: k, psi, t, delta
!     real :: f, kt, pd

!     kt = k * t
!     pd = abs(psi) * delta
!     call rtnewt(funcd, kt, f)

!   contains
!     pure subroutine funcd(x, fval, fderiv)
!       real, intent(in) :: x
!       real, intent(out) :: fval, fderiv

!       fval = x - pd * log(abs(1. + x / pd)) - kt
!       fderiv = x / (pd + x)
!     end subroutine
!   end function

!   subroutine waini(name, dt, x, y)
!   end subroutine

!   pure subroutine waseep(t, dt, v0, beta, ksat, i0, n, r, z)
!     real, intent(in) :: t, dt, v0, beta(2), ksat, i0, n, r
!     real, intent(out) :: z(:)
!     real, dimension(size(z)) :: k, psi
!     real :: eps, se(0:size(z)), se0, z0
!     integer :: i, ii, j, jj, m, num

!     num = size(z)
!     eps = 1. / (2. * num)
!     i = num
!     se0 = (i0 - r) / (n - r)
!     se(0) = 0
!     do j = 1, num
!       jj = j - 1
!       se(j) = real(j) / num
!       if(abs(se(jj) - se0) < eps) i = jj
!     end do

!     swc%a = beta(1)
!     swc%n = beta(2)
!     z(:i) = h
!     z0 = waga(ksat, swc%matsuc(se0), dt, n - i0)
!     where (z == 0) z = z0
!     do j = i, num - 1
!       jj = j + 1
!       k(j) = ksat * simpr(swc%relhc, se(j), se(jj))
!       psi(j) = simpr(swc%matsuc, se(j), se(jj))
!     end do
!   contains
!     pure subroutine seep
!       integer, save :: d = num

!       dz = c * (k(d) * psi(d) / z(j) + k(d))
!       z(j) = z(j) + dz * ts
!     end subroutine
!   end subroutine
! end module
