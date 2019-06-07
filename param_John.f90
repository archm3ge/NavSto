! ==========================================================
! ----------------------------------------------------------
! ++++++++++  GLOBAL PARAMETERS AND VARIABLES ++++++++++++++
! ----------------------------------------------------------

! Note:
! "ind" dimensionless wavenum. "dim" dimensional. dim wavenum = 2*pi/L ind wavenum
! --- make sure Nx & Ny are even preferably powers of small primes (like 2,3)

module param_John
  use, intrinsic :: iso_c_binding
  implicit none

! ----------------------------------------------------------
! set the dimension of space
  integer, parameter       :: Nx = 256, Ny = 256            ! physical grid resolution in x and y
  integer, parameter       :: ktx = Nx/2, kty=Ny/2          ! max ind wavenumber
  integer, parameter       :: iktx = ktx+1,ikty=Ny          ! total number of wavenumbers
  real, parameter          :: fftnorm=float(Nx*Ny)          ! fft-size to normalize the transform

! ----------------------------------------------------------
! misc constants
  real, parameter          :: pi = 2.*asin(1.)
  complex, parameter       :: zi = cmplx(0.,1.)
  real, parameter          :: Lx = 2*pi, Ly = 2*pi          ! domain size in physical space

! ----------------------------------------------------------
! wavenumbers
  real, save               :: kxa(iktx), kya(ikty)          ! dim wavenumber
  integer, dimension(iktx,ikty),save :: L                   ! Mask for wavenumber truncation

! ----------------------------------------------------------
! time
  real, save               :: time
  real, save               :: delt = 0.01/256

end module param_John
