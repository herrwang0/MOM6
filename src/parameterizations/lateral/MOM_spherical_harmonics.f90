!> Laplace's spherical harmonics transforms
module MOM_spherical_harmonics
use MOM_cpu_clock,     only : cpu_clock_id, cpu_clock_begin, cpu_clock_end, &
  CLOCK_MODULE
use MOM_error_handler, only : MOM_error, MOM_mesg, FATAL, WARNING
use MOM_file_parser,   only : get_param, log_version, param_file_type
use MOM_grid,          only : ocean_grid_type
use MOM_unit_scaling,  only : unit_scale_type
use MOM_coms_infra, only : sum_across_PEs

implicit none ; private

public spherical_harmonics_init, spherical_harmonics_end, spherical_harmonics_forward, spherical_harmonics_inverse, SHOrderDegreeToIndex

#include <MOM_memory.h>

type, public :: sht_CS
  logical :: initialized = .False. !< True if this control structure has been initialized.
  integer :: nOrder !< Maximum degree of associated Legendre polynomials [nodim].
  integer :: lmax !< Total number of associated Legendre polynomials = (nOrder+1)*(nOrder+2)/2 [nodim].
  real, allocatable :: sinCoLatT(:,:), cosCoLatT(:,:) !< Precomputed sine and cosine of colatitude at the t-cells [nondim].
  real, allocatable :: complexFactorRe(:,:,:), complexFactorIm(:,:,:), complexExpRe(:,:,:), complexExpIm(:,:,:)
  real, allocatable :: aRecurrenceCoeff(:,:), bRecurrenceCoeff(:,:)
  real, allocatable :: &
    pmn(:,:),   & !< P(n,m)
    pmnm2(:,:), & !< P(n-2,m)
    pmnm1(:,:)    !< P(n-1,m)
  logical :: bfb
end type sht_CS

contains
subroutine spherical_harmonics_forward(G, CS, var, SnmRe, SnmIm)
  type(ocean_grid_type), intent(in) :: G
  type(sht_CS), intent(inout) :: CS
  real, intent(in)  :: var(:,:)
  real, intent(out) :: SnmRe(:), SnmIm(:)

  integer :: i, j
  integer :: is, ie, js, je
  integer :: m, n, l
  real, allocatable :: Snm_local(:), SnmRe_local(:), SnmIm_local(:)

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  if (CS%bfb) then
    ! allocate(Snm_local_reproSum(2*CS%lmax)); Snm_local_reproSum = 0.0
    ! allocate(SnmRe_local_reproSum(CS%lmax)); SnmRe_local_reproSum = 0.0
    ! allocate(SnmIm_local_reproSum(CS%lmax)); SnmIm_local_reproSum = 0.0
  else
    allocate(Snm_local(2*CS%lmax)); Snm_local = 0.0
    allocate(SnmRe_local(CS%lmax)); SnmRe_local = 0.0
    allocate(SnmIm_local(CS%lmax)); SnmIm_local = 0.0
  endif

  do m = 0, CS%nOrder
    !------------
    ! n = m
    !------------
    n = m
    ! Calculate associated Legendre polynomial for n=m (output pmnm2)
    call associatedLegendrePolynomials(n, m, l, CS, G)

    ! Compute local integral contribution
    if (CS%bfb) then
      ! do iCell = endIdx, startIdx, -1
      !     SnmRe_local_reproSum(iCell,l) = SnmRe_local_reproSum(iCell,l) + sshSmoothed(iCell)*pmnm2(iCell)*complexFactorRe(iCell,m+1)
      !     SnmIm_local_reproSum(iCell,l) = SnmIm_local_reproSum(iCell,l) + sshSmoothed(iCell)*pmnm2(iCell)*complexFactorIm(iCell,m+1)
      ! enddo
    else
      do j = js,je ; do i = is,ie
        SnmRe_local(l) = SnmRe_local(l) + var(i,j) * CS%pmnm2(i,j) * CS%complexFactorRe(i,j,m+1)
        SnmIm_local(l) = SnmIm_local(l) + var(i,j) * CS%pmnm2(i,j) * CS%complexFactorIm(i,j,m+1)
      enddo ; enddo
    endif

    !------------
    ! n = m+1
    !------------
    n = m+1
    if (n <= CS%nOrder) then
      ! Calculate associated Legendre polynomial for n = m+1 using recurrence relationship
      call associatedLegendrePolynomials(n, m, l, CS, G)

      ! Compute local integral contribution
      if (CS%bfb) then
        ! do iCell = endIdx, startIdx, -1
        !     SnmRe_local_reproSum(iCell,l) = SnmRe_local_reproSum(iCell,l) + sshSmoothed(iCell)*pmnm1(iCell)*complexFactorRe(iCell,m+1)
        !     SnmIm_local_reproSum(iCell,l) = SnmIm_local_reproSum(iCell,l) + sshSmoothed(iCell)*pmnm1(iCell)*complexFactorIm(iCell,m+1)
        ! enddo
      else
        do j = js,je ; do i = is,ie
          SnmRe_local(l) = SnmRe_local(l) + var(i,j) * CS%pmnm1(i,j) * CS%complexFactorRe(i,j,m+1)
          SnmIm_local(l) = SnmIm_local(l) + var(i,j) * CS%pmnm1(i,j) * CS%complexFactorIm(i,j,m+1)
        enddo ; enddo
      endif
    endif

    !------------
    ! n > m+1
    !------------
    do n = m+2,CS%nOrder
      ! Calculate associated Legendre polynomial using recurrence relationship
      ! call associatedLegendrePolynomials(n, m, l, CS%pmnm2, CS%pmnm1, CS%pmn)
      call associatedLegendrePolynomials(n, m, l, CS, G)

      ! Update associated Ledgendre polynomial values for next recurrence
      do j = js,je ; do i = is,ie
        CS%pmnm2(i,j) = CS%pmnm1(i,j)
        CS%pmnm1(i,j) = CS%pmn(i,j)
      enddo ; enddo

      ! Compute local integral contribution
      if (CS%bfb) then
        ! do iCell = endIdx, startIdx, -1
        !     SnmRe_local_reproSum(iCell,l) = SnmRe_local_reproSum(iCell,l) + sshSmoothed(iCell)*pmn(iCell)*complexFactorRe(iCell,m+1)
        !     SnmIm_local_reproSum(iCell,l) = SnmIm_local_reproSum(iCell,l) + sshSmoothed(iCell)*pmn(iCell)*complexFactorIm(iCell,m+1)
        ! enddo
      else
        do j = js,je ; do i = is,ie
          SnmRe_local(l) = SnmRe_local(l) + var(i,j) * CS%pmn(i,j) * CS%complexFactorRe(i,j,m+1)
          SnmIm_local(l) = SnmIm_local(l) + var(i,j) * CS%pmn(i,j) * CS%complexFactorIm(i,j,m+1)
        enddo ; enddo
      endif
    enddo ! n loop
  enddo ! m loop

  ! call mpas_timer_stop('Parallel SAL: Forward Transform')

  ! call mpas_timer_start('Parallel SAL: Communication')
  if (CS%bfb) then
    ! do m = 1,lmax
    !     do iCell = 1,nCellsOwned 
    !       Snm_local_reproSum(iCell,m) = SnmRe_local_reproSum(iCell,m)
    !       Snm_local_reproSum(iCell,lmax+m) = SnmIm_local_reproSum(iCell,m)
    !     enddo
    ! enddo
  else
    do m = 1,CS%lmax
      Snm_local(m) = SnmRe_local(m)
      Snm_local(CS%lmax+m) = SnmIm_local(m)
    enddo
  endif

  ! Compute global integral by summing local contributions
  if (CS%bfb) then
     !threadNum = mpas_threading_get_thread_num()
     !if ( threadNum == 0 ) then
        !  Snm = mpas_global_sum_nfld(Snm_local_reproSum,dminfo%comm)
     !endif
  else
    call sum_across_PEs(Snm_local, 2*CS%lmax)
  endif

  do m = 1,CS%lmax
    SnmRe(m) = Snm_local(m)
    SnmIm(m) = Snm_local(CS%lmax+m)
  enddo
end subroutine spherical_harmonics_forward

subroutine spherical_harmonics_inverse(G, CS, SnmRe, SnmIm, var)
  type(ocean_grid_type), intent(in) :: G
  type(sht_CS), intent(inout) :: CS
  real, intent(out)  :: var(:,:)
  real, intent(in) :: SnmRe(:), SnmIm(:)

  integer :: is, ie, js, je
  integer :: m, n, l
  integer :: i, j
  real :: mFac
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  var = 0.0
  do m = 0,CS%nOrder
    if (m>0) then
      mFac = 2.0
    else
      mFac = 1.0
    endif

    !------------
    ! n = m
    !------------
    n = m
    ! Calculate associated Legendre polynomial using recurrence relationship
    ! call associatedLegendrePolynomials(n, m, l, CS%pmnm2, CS%pmnm1, CS%pmn)
    call associatedLegendrePolynomials(n, m, l, CS, G)

    ! Sum together product of spherical harmonic functions and coefficients
    do j = js,je ; do i = is,ie
      var(i,j) = var(i,j) &
        + mFac * CS%pmnm2(i,j) * (SnmRe(l) * CS%complexExpRe(i,j,m+1) + SnmIm(l) * CS%complexExpIm(i,j,m+1))
    enddo ; enddo

    !------------
    ! n = m+1
    !------------
    n = m+1
    if (n <= CS%nOrder) then
      ! Calculate associated Legendre polynomial using recurrence relationship
      ! call associatedLegendrePolynomials(n, m, l, CS%pmnm2, CS%pmnm1, CS%pmn)
      call associatedLegendrePolynomials(n, m, l, CS, G)

      ! Sum together product of spherical harmonic functions and coefficients
      do j = js,je ; do i = is,ie
        var(i,j) = var(i,j) &
          + mFac * CS%pmnm1(i,j) * (SnmRe(l) * CS%complexExpRe(i,j,m+1) + SnmIm(l) * CS%complexExpIm(i,j,m+1))
      enddo ; enddo
    endif

    !------------
    ! n > m+1
    !------------
    do n = m+2,CS%nOrder
      ! Calculate associated Legendre polynomial using recurrence relationship
      ! call associatedLegendrePolynomials(n, m, l, CS%pmnm2, CS%pmnm1, CS%pmn)
      call associatedLegendrePolynomials(n, m, l, CS, G)

      ! Update associated Ledgendre polynomial values for next recurrence
      do j = js,je ; do i = is,ie
        CS%pmnm2(i,j) = CS%pmnm1(i,j)
        CS%pmnm1(i,j) = CS%pmn(i,j)
      enddo ; enddo

      ! Sum together product of spherical harmonic functions and coefficients
      do j = js,je ; do i = is,ie
        var(i,j) = var(i,j) &
          + mFac * CS%pmn(i,j) * (SnmRe(l) * CS%complexExpRe(i,j,m+1) + SnmIm(l) * CS%complexExpIm(i,j,m+1))
      enddo ; enddo
    enddo ! n loop
  enddo ! m loop
end subroutine spherical_harmonics_inverse

subroutine spherical_harmonics_init(G, param_file, CS)
  type(ocean_grid_type),   intent(in) :: G    !< The ocean's grid structure.
  type(sht_CS), intent(inout) :: CS
  type(param_file_type),   intent(in)    :: param_file !< A structure indicating
  real :: PI                   ! 3.1415926... calculated as 4*atan(1)
  real :: RADIAN
  integer :: is, ie, js, je
  integer :: i, j
  integer :: m, n
  ! This include declares and sets the variable "version".
# include "version_variable.h"
  character(len=40)  :: mdl = "MOM_spherical_harmonics" ! This module's name.


  PI = 4.0*atan(1.0)
  RADIAN = PI / 180.0

  call log_version(param_file, mdl, version, "")

  call get_param(param_file, '', "SHT_DEGREE", CS%nOrder, &
                 "Order of spherical harmonics transformation. ", default=0)
  CS%lmax = (CS%nOrder + 1) * (CS%nOrder + 2) / 2
  call get_param(param_file, '', "SHT_BFB", CS%bfb, &
                 "If true, use bfb sum. Default is False.", default=.False.)

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  allocate(CS%pmn(is:ie,js:je))  ; CS%pmn(:,:)   = 0.0
  allocate(CS%pmnm1(is:ie,js:je)); CS%pmnm1(:,:) = 0.0
  allocate(CS%pmnm2(is:ie,js:je)); CS%pmnm2(:,:) = 0.0

  ! Compute recurrence relationship coefficients
  allocate(CS%aRecurrenceCoeff(CS%nOrder+1,CS%nOrder+1)); CS%aRecurrenceCoeff(:,:) = 0.0
  allocate(CS%bRecurrenceCoeff(CS%nOrder+1,CS%nOrder+1)); CS%bRecurrenceCoeff(:,:) = 0.0

  do m = 0, CS%nOrder
    do n = m, CS%nOrder
      if (m /= n) then
        CS%aRecurrenceCoeff(n+1,m+1) = sqrt(real((2*n-1)*(2*n+1)) / real((n-m)*(n+m)))
        CS%bRecurrenceCoeff(n+1,m+1) = sqrt(real((2*n+1)*(n+m-1)*(n-m-1)) / real((n-m)*(n+m)*(2*n-3)))
      endif
    enddo
  enddo

  ! Precompute complex exponential factors
  allocate(CS%complexFactorRe(is:ie, js:je, CS%nOrder+1)); CS%complexFactorRe(:,:,:) = 0.0
  allocate(CS%complexFactorIm(is:ie, js:je, CS%nOrder+1)); CS%complexFactorIm(:,:,:) = 0.0
  allocate(CS%complexExpRe(is:ie, js:je, CS%nOrder+1)); CS%complexExpRe(:,:,:) = 0.0
  allocate(CS%complexExpIm(is:ie, js:je, CS%nOrder+1)); CS%complexExpIm(:,:,:) = 0.0

  do m = 0, CS%nOrder
    do j = js,je ; do i = is,ie
      CS%complexExpRe(i, j, m+1)    = cos(real(m) * (G%geolonT(i,j)*RADIAN))
      CS%complexExpIm(i, j, m+1)    = sin(real(m) * (G%geolonT(i,j)*RADIAN))
      CS%complexFactorRe(i, j, m+1) = CS%complexExpRe(i, j, m+1) * G%areaT(i,j) / G%Rad_Earth**2
      CS%complexFactorIm(i, j, m+1) = CS%complexExpIm(i, j, m+1) * G%areaT(i,j) / G%Rad_Earth**2
    enddo ; enddo
  enddo

  ! allocate(Snm_local(2*lmax),Snm(2*lmax))
  ! allocate(SnmRe_local(lmax),SnmRe(lmax))
  ! allocate(SnmIm_local(lmax),SnmIm(lmax))
  ! allocate(SnmIm_local_reproSum(nCellsOwned,lmax))
  ! allocate(SnmRe_local_reproSum(nCellsOwned,lmax))
  ! allocate(Snm_local_reproSum(nCellsOwned,2*lmax))


  ! Pre-compute sin and cos of latCell (co-latitude) values
  allocate(CS%sinCoLatT(is:ie,js:je)); CS%sinCoLatT(:,:) = 0.0
  allocate(CS%cosCoLatT(is:ie,js:je)); CS%cosCoLatT(:,:) = 0.0
  do j = js,je ; do i = is,ie
    CS%sinCoLatT(i,j) = sin(0.5*PI - G%geolatT(i,j)*RADIAN)
    CS%cosCoLatT(i,j) = cos(0.5*PI - G%geolatT(i,j)*RADIAN)
  enddo ; enddo
end subroutine spherical_harmonics_init

subroutine spherical_harmonics_end(CS)
  type(sht_CS), intent(inout) :: CS

  deallocate(CS%sinCoLatT, CS%cosCoLatT)
  deallocate(CS%complexFactorRe, CS%complexFactorIm, CS%complexExpRe, CS%complexExpIm)
  deallocate(CS%pmnm2, CS%pmnm1, CS%pmn)
  deallocate(CS%aRecurrenceCoeff, CS%bRecurrenceCoeff)
end subroutine spherical_harmonics_end

subroutine associatedLegendrePolynomials(n, m, l, CS, G)
  type(ocean_grid_type),   intent(in) :: G
  integer, intent(in) :: n
  integer, intent(in) :: m
  integer, intent(out) :: l
  type(sht_CS), intent(inout) :: CS
  ! real, dimension(:,:), intent(inout) :: pmnm2
  ! real, dimension(:,:), intent(inout) :: pmnm1
  ! real, dimension(:,:), intent(inout) :: pmn

  integer :: i, j, k
  real, parameter :: PI = 4.0*atan(1.0)                 ! 3.1415926... calculated as 4*atan(1)
  integer :: is, ie, js, je
 is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  l = SHOrderDegreeToIndex(n,m, CS%nOrder)

  if (n == m) then
    do j = js, je ; do i = is, ie
      CS%pmnm2(i,j) = sqrt(1.0/(4.0*PI)) * CS%sinCoLatT(i,j)**m
      do k = 1, m
        CS%pmnm2(i,j) = CS%pmnm2(i,j) * sqrt(real(2*k+1)/real(2*k))
      enddo
    enddo ; enddo
  else if (n == m+1) then
    do j = js, je ; do i = is, ie
      CS%pmnm1(i,j) = CS%aRecurrenceCoeff(n+1,m+1) * CS%cosCoLatT(i,j) * CS%pmnm2(i,j)
    enddo ; enddo
  else
    do j = js, je ; do i = is, ie
      CS%pmn(i,j) =  CS%aRecurrenceCoeff(n+1,m+1) * CS%cosCoLatT(i,j) * CS%pmnm1(i,j) &
                   - CS%bRecurrenceCoeff(n+1,m+1) * CS%pmnm2(i,j)
    enddo ; enddo
  endif
end subroutine associatedLegendrePolynomials

function SHOrderDegreeToIndex(n,m, nOrder) result(l)!{{{

  integer :: l
  integer :: n
  integer :: m
  integer :: nOrder

  l = (nOrder+1)*m - m*(m+1)/2 + n+1

end function SHOrderDegreeToIndex

end module MOM_spherical_harmonics