module mos
  use bezier

contains
  subroutine mossearch(rd, num, x0, xn, model, slice)
    interface
      pure function model(x)
        real, intent(in) :: x(:)
        real :: model
      end function
      pure function slice(x)
        real, intent(in) :: x(:)
        real, dimension(:) :: slice
      end function
    end interface

    p0 = (x0, omega(x0))
    pn = (xn, omega(xn))
    t = 0 .. 1
    p = asarc(p0, pn)
    fosmin = amoeba(fos(p))
    p = asline(p0, pn)
    fosmin = amoeba(fos(p))
  end subroutine

  pure function fos(x)
    xy = bez(t, p0, x, pn)

    xy(1:, 1) - xy(:n, 1) < 0
    omega(xy(:, 1)) - xy(:, 2) > 0
    xy(:, 2) - rd * norm(p0-pn) < 0
      return 'super high number'
    slices = slice(xy)
    fos = model(slices)
  end function
end module
