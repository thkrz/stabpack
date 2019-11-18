module rvcont
  implicit none
  private
  public loglap
  public lognrm

  type loglap_gen
    real d
    real a
    real b
  contains
    procedure :: cdf => loglap_cdf
    procedure :: fit => loglap_fit
    procedure :: pdf => loglap_pdf
    procedure :: ppf => loglap_ppf
  end type

  type lognrm_gen
    real m
    real s
  contains
    procedure :: cdf => lognrm_cdf
    procedure :: fit => lognrm_fit
    procedure :: pdf => lognrm_pdf
    procedure :: ppf => lognrm_ppf
  end type

  type(loglap_gen) :: loglap
  type(lognrm_gen) :: lognrm

contains
  elemental function loglap_cdf(self, x) result(p)
    class(loglap_gen), intent(in) :: self
    real, intent(in) :: x
    real :: ab, p

    p = 0
    if(x < 0) return

    ab = self%a / (self%a + self%b)
    if(x < self%d) then
      p = ab * (x / self%d)**self%b
    else
      p = 1. - ab * (self%d / x)**self%a
    end if
  end function

  subroutine loglap_fit(self, y)
    use fmin, only: brent
    class(loglap_gen), intent(inout) :: self
    real, intent(in), dimension(:) :: y
    real :: fmin, xmin
    integer :: n

    n = size(y)
    fmin = brent(0., 1., 2., h, 1.0e-4, xmin)
    self%d = xmin

  contains
    function h(d)
      real, intent(in) :: d
      real :: h

      call tail(d, self%a, self%b)
      h = log(self%a + self%b) + self%b + (self%a * self%b) &
        / (self%a + self%b)
    end function

    pure subroutine tail(d, a, b)
      real, intent(in) :: d
      real, intent(out) :: a, b
      real :: p, q, r
      integer :: i

      p = 1
      q = 1
      do i = 1, n
        p = p * merge(y(i) / d, 1., d /= 0)
        q = q * merge(d / y(i), 1., y(i) /= 0)
      end do
      r = sqrt(log(p) * log(q))
      a = n / (log(p) + r)
      b = n / (log(q) + r)
    end subroutine
  end subroutine

  elemental function loglap_pdf(self, x) result(p)
    class(loglap_gen), intent(in) :: self
    real, intent(in) :: x
    real :: c, p

    p = 0
    if(x <= 0) return

    c = 1. / self%d * (self%a * self%b) / (self%a + self%b)
    if(x < d) then
      p = c * (x / slef%d)**(self%b - 1)
    else
      p = c * (self%d / x)**(self%a + 1)
    end if
  end function

  elemental function loglap_ppf(p, d, a, b) result(x)
    class(loglap_gen), intent(in) :: self
    real, intent(in) :: p
    real :: ab, x

    x = 0
    if(p <= 0 .or. p >= 1) return

    ab = self%a / (self%a + self%b)
    if(p <= ab) then
      x = self%d * (p * (self%a + self%b) / self%a)**(1. / self%b)
    else
      x = self%d * ((1. - p) * (self%a + self%b) / self%b)**(-1. / self%a)
    end if
  end function

  elemental function lognrm_cdf(self, x) result(p)
    class(lognrm_gen), intent(in) :: self
    real, intent(in) :: x
    real :: p

    p = .5 + .5 * erf(log(x - self%m) / (sqrt(2.) * self%s))
  end function

  subroutine lognrm_fit(self, y)
    class(lognrm_gen), intent(inout) :: self
    real :: p
    integer :: n

    n = size(y)
    p = 0
    do i = 1, n
      p = p + log(y(i))
    end do
    self%m = p / n
    self%s = 0
    do i = 1, n
      self%s = self%s + (log(y(i)) - p)**2
    end do
    self%s = sqrt(self%s / n)
  end subroutine

  elemental function lognrm_pdf(self, x) result(p)
    class(lognrm_gen), intent(in) :: self
    real, intent(in) :: x
    real :: p

    p = 1. / (x * self%s * sqrt(2. * pi)) * exp(-(log(x) - self%m)**2 / (2. * self%s**2))
  end function

  elemental function lognrm_ppf(self, p) result(x)
    class(lognrm_gen), intent(in) :: self
    real, intent(in) :: p
    real :: x

    x = 1
  end function
end module
