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
  use intgrt, only: nsimp
  use rtfind only: rtnewt
  use vg
  implicit none
  private
  public waga

  type(vg_t) :: swc

contains
  pure function waga(k, psi, time, delta) result(f)
    real, intent(in) :: k, psi, time, delta
    real :: f, kt, pd

    kt = k * time
    pd = abs(psi) * delta
    call rtnewt(funcd, kt, f)

  contains
    pure subroutine funcd(x, fval, fderiv)
      real, intent(in) :: x
      real, intent(out) :: fval, fderiv

      fval = x - pd * log(abs(1. + x / pd)) - kt
      fderiv = x / (pd + x)
    end subroutine
  end subroutine

  subroutine wasini(name, dt, x, y)
  end subroutine

  pure subroutine waseep(t, dt, v, h, beta, ksat, i0, n, r, z)
    real, intent(in), :: t, dt, v, h, beta, ksat, i0, n, r
    real, intent(inout) :: z(:)
    real, dimension(size(z)) :: k, psi
    real :: dn, dz, se(0:size(z)), se0, tt, vv, z0, ztemp
    integer :: d, i, ii, j, jj, m, num

    num = size(z)
    se(0) = 0
    se0 = (i0 - r) / (n - r)
    do j = 1, num
      jj = j - 1
      se(j) = real(j) / num
      k(j) = ksat * nsimp(swc%relhc, se(jj), se(j))
      psi(j) = nsimp(swc%matsuc, se(jj), se(j))
    end do

    i = nint(i0 * num / n)
    z(:i) = h
    z0 = waga(ksat, swc%matsuc(i0), dt, n - i0)
    d = num
    dn = 1. / num
    tt = t
    vv = v
    do while(tt > 0)
      do j = i + 1, d
        if(z(j) == 0) then
          z(j) = z0
          vv = vv - z0 * dn
        end if
        if(vv <= 0) exit
        ! TODO: relative water content ok?
        dz = 1. / (1. - se(d)) * (k(d) * psi(d) / z(j) + k(d))
        ztemp = dz * dt
        vv = vv - ztemp * dn
        z(j) = z(j) + ztemp
      end do
      tt = tt - dt
    end do
  contains
  ! contains
  !   pure subroutine seep
  !     integer, save :: d = num

  !     dz = c * (k(d) * psi(d) / z(j) + k(d))
  !     z(j) = z(j) + dz * ts
  !   end subroutine
  end subroutine
end module
