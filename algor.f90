!****m* algor/algor ===========================================================!
! NAME                                                                         !
!   algor                                                                      !
!------------------------------------------------------------------------------!
! DESCRIPTION                                                                  !
!------------------------------------------------------------------------------!
! AUTHOR                                                                       !
!   Aaron Hopkinson                                                            !
!****==========================================================================!
module algor
  use constants, only: dp
  implicit none
  private

  integer,            parameter                   ::  i32 = selected_int_kind(9)
  logical,                          save          ::  rand_initialized = .false.
  integer(kind=i32),  dimension(4), save, public  ::  rand_seed

  public  ::  algor_init_rand
  public  ::  algor_uniform_rand
  public  ::  algor_gauss_rand
  public  ::  algor_cross_product
  public  ::  algor_symmetric_matrix_index
  public  ::  algor_factorial
  public  ::  algor_3d_rotation
  public  ::  algor_exp
  public  ::  algor_calc_inverse_real
  public  ::  algor_kronecker_delta

contains

!****s* algor/algor_init_rand =================================================!
! NAME                                                                         !
!   algor_init_rand (PUBLIC)                                                   !
!------------------------------------------------------------------------------!
! DESCRIPTION                                                                  !
!   Initializes the random number generator.                                   !
!------------------------------------------------------------------------------!
! ARGUMENTS                                                                    !
!   integer,  optional,     intent(in)    ::  seed                             !
!------------------------------------------------------------------------------!
! AUTHOR                                                                       !
!   Based on Matt Probert's simple MD code                                     !
!==============================================================================!
subroutine algor_init_rand(seed)
  use hash_tables,  only: hash_function
  use io,           only: io_err, io_str_to_real
  implicit none
  integer,  optional,     intent(in)    ::  seed

  !f90 intrinsic time call
  character(len=10) ::  system_time             !length is crucial ...
  character(len=11) ::  str
  real(kind=dp)     ::  rtime
  integer           ::  i, istat

  if (present(seed)) then
    ! fixed seed
    rand_seed(1) = int(seed,i32)
  else
    ! get seed based on time
    call date_and_time(time=system_time)        !character string hhmmss.xxx
    rtime = io_str_to_real(system_time)         !convert to real
    rand_seed(1) = int((rtime * 1000.0_dp),i32) !0<rtime<235959999.0 which fits within huge(1)
  end if

  ! xorshift RNG requires 4 numbers..
  do i = 1, size(rand_seed,1)
    rand_seed(i) = rand_seed(1) + int(i-1,i32)*4_i32

    ! convert number to string..
    write(str,'(I11)',iostat=istat) rand_seed(i)
    if (istat .ne. 0) call io_err("algor_init_rand: Error converting integer to string")

    ! pass string to hash function and return an integer... set this integer to be the seed
    rand_seed(i) = hash_function(str)
  end do

  rand_initialized = .true.
end subroutine algor_init_rand


!****f* algor/algor_uniform_rand ==============================================!
! NAME                                                                         !
!   algor_uniform_rand (PUBLIC)                                                !
!------------------------------------------------------------------------------!
! DESCRIPTION                                                                  !
!------------------------------------------------------------------------------!
! ARGUMENTS                                                                    !
!   real(kind=dp),  optional, intent(in)  ::  rmin, rmax                       !
!==============================================================================!
function algor_uniform_rand(rmin, rmax)
  use io, only: io_err
  implicit none
  real(kind=dp),  optional, intent(in)  ::  rmin, rmax
  real(kind=dp)                         ::  algor_uniform_rand

  real(kind=dp)                         ::  int_min, int_max    ! internal min and max
  real(kind=dp)                         ::  rscale

  if (.not. rand_initialized) call algor_init_rand()

  algor_uniform_rand = algor_rand_xorshift()

  if (present(rmin) .and. present(rmax)) then
    ! make sure int_min is smallest, int_max is largest
    if (rmin < rmax) then
      ! what we expect/hope for
      int_min = rmin
      int_max = rmax
    else if (rmin > rmax) then
      ! swap
      int_min = rmax
      int_max = rmin
    else if (rmin == rmax) then
      call io_err("algor_uniform_rand: 2 distinct numbers expected")
    end if

    rscale = int_max-int_min
    algor_uniform_rand = rscale*algor_uniform_rand    ! rescale
    algor_uniform_rand = algor_uniform_rand + int_min ! shift
  else if (present(rmin) .or. present(rmax)) then
    call io_err("algor_uniform_rand: Only 1 argument provided - expect 2")
  end if
end function algor_uniform_rand

function algor_rand_xorshift()
  ! Marsaglia xorshift RNG - based on 'xor128' - should have period of 2^128 - 1
  ! (May be less because of mapping to positive numbers with addition)
  ! For ref: Paper initial seed: x = 123456789, y = 362436069, z = 521288629, w = 88675123
  implicit none
  real(kind=dp)                 ::  algor_rand_xorshift
  integer(kind=i32),  parameter ::  a = 11_i32, b = 8_i32, c = 19_i32
  integer(kind=i32),  parameter ::  huge_i32 = huge(1_i32)
  real(kind=dp),      parameter ::  r_huge_i32 = real(huge_i32,dp)
  ! local vars:
  integer(kind=i32) ::  x, y, z, w, t

  x = rand_seed(1)
  y = rand_seed(2)
  z = rand_seed(3)
  w = rand_seed(4)

  t = ieor(x, ishft(x,  a))
  t = ieor(t, ishft(t, -b))
  x = y
  y = z
  z = w
  w = ieor(w, ishft(w, -c))
  w = ieor(w, t)
  t = w

  rand_seed(:) = (/ x, y, z, w /)

  if (t .lt. 0_i32) t = t + huge_i32
  algor_rand_xorshift = real(t,dp)/r_huge_i32
end function algor_rand_xorshift

!****f* algor/algor_gauss_rand ================================================!
! NAME                                                                         !
!   algor_gauss_rand (PUBLIC)                                                  !
!------------------------------------------------------------------------------!
! DESCRIPTION                                                                  !
!------------------------------------------------------------------------------!
! ARGUMENTS                                                                    !
!   real(kind=dp),  optional, intent(in)  ::  mean, stddev                     !
!==============================================================================!
function algor_gauss_rand(mean, stddev)
  use io, only: io_err
  implicit none
  ! generate Gaussian distributed random numbers
  ! using the Marsaglia polar method
  ! default behaviour (no arguments) is a normally distributed number
  ! (ie: mean = 0, stddev = 1)
  real(kind=dp),  optional, intent(in)  ::  mean, stddev
  real(kind=dp)                         ::  algor_gauss_rand

  real(kind=dp)   ::  w, x, y
  real(kind=dp)   ::  int_mean, int_stddev
  real(kind=dp),  save  ::  extra
  logical,        save  ::  stored = .false.

  if (.not.(present(stddev) .and. present(mean))) then
    ! if neither do normal distribution
    int_stddev = 1.0_dp
    int_mean = 0.0_dp
  else if (present(stddev) .and. present(mean)) then
    ! use either normal or user specified stddev and mean
    int_stddev = stddev
    int_mean = mean
  else
    !if (present(stddev) .or. present(mean)) - ie: only one argument
    call io_err("algor_gauss_rand: Missing an argument")
  end if

  if (stored) then
    ! use the results from last time
    ! rescale/shift for new stddev and mean
    extra = extra*int_stddev + int_mean
    algor_gauss_rand = extra
    stored = .false.  ! reset
  else
    w = 2.0_dp ! > 1
    do while ((w .ge. 1.0_dp) .or. (w .eq. 0.0_dp))
      x = algor_uniform_rand(-1.0_dp, 1.0_dp)
      y = algor_uniform_rand(-1.0_dp, 1.0_dp)
      w = x*x + y*y
    end do

    w = sqrt((-2.0_dp*log(w))/w)
    algor_gauss_rand = int_stddev*x*w + int_mean
    extra = y*w     ! for next time
    stored = .true.
  end if
end function algor_gauss_rand


!****f* algor/ ================================================================!
! NAME                                                                         !
!------------------------------------------------------------------------------!
! DESCRIPTION                                                                  !
!------------------------------------------------------------------------------!
! ARGUMENTS                                                                    !
!==============================================================================!
function algor_cross_product(a,b)
  implicit none
  real(kind=dp),  dimension(3), intent(in)  ::  a, b
  real(kind=dp),  dimension(3)              ::  algor_cross_product

  algor_cross_product(1) = a(2)*b(3)-b(2)*a(3)
  algor_cross_product(2) = b(1)*a(3)-a(1)*b(3)
  algor_cross_product(3) = a(1)*b(2)-b(1)*a(2)
end function algor_cross_product


!****f* algor/ ================================================================!
! NAME                                                                         !
!------------------------------------------------------------------------------!
! DESCRIPTION                                                                  !
!------------------------------------------------------------------------------!
! ARGUMENTS                                                                    !
!==============================================================================!
function algor_symmetric_matrix_index(row,column,n)
  ! symmetric so it doesn't actually matter which is row and which is column
  ! these are just named to make things easier to compare with rest of code...
  implicit none
  integer,  intent(in)  ::  row, column, n
  integer               ::  i, j
  integer               ::  algor_symmetric_matrix_index

  ! function gives index of a 1D array that represents
  ! the (row, column) value of an N^2 symmetric array

  ! need to use right i and j
  if (row .le. column) then
    ! upper triangle
    i = row
    j = column
  else
    ! lower triangle
    i = column
    j = row
  end if

  ! for upper triangle: index = row*N - row(row-1)/2 - N + col
  algor_symmetric_matrix_index = i*n - (i*(i-1))/2 - n + j
end function algor_symmetric_matrix_index

!****f* algor/algor_factorial =================================================!
! NAME                                                                         !
!   Returns the factorial of the input.                                        !
!------------------------------------------------------------------------------!
! DESCRIPTION                                                                  !
!------------------------------------------------------------------------------!
! ARGUMENTS                                                                    !
!------------------------------------------------------------------------------!
! AUTHOR                                                                       !
!   Aaron Hopkinson/Edward Higgins?                                            !
!==============================================================================!
recursive function algor_factorial(i) result(res)
  implicit none
  integer               ::  res
  integer,  intent(in)  ::  i

  if (i .ne. 0) then
    res = i*algor_factorial(i-sign(1,i))
    if (i .lt. 0) res = sign(res,-1)
  else
    res = 1
  end if
end function algor_factorial

!****f* algor/algor_3d_rotation ===============================================!
! NAME                                                                         !
!   Performs 3d rotation of a vector                                           !
!------------------------------------------------------------------------------!
! DESCRIPTION                                                                  !
!------------------------------------------------------------------------------!
! ARGUMENTS                                                                    !
!------------------------------------------------------------------------------!
! AUTHOR                                                                       !
!   Aaron Hopkinson                                                            !
!==============================================================================!
function algor_3d_rotation(vector, angles)
  implicit none
  real(kind=dp),  dimension(3), intent(in)  ::  vector
  real(kind=dp),  dimension(3), intent(in)  ::  angles
  real(kind=dp),  dimension(3)              ::  algor_3d_rotation

  real(kind=dp),  dimension(3,3)  ::  rotation_matrix
  real(kind=dp) ::  cosalpha, sinalpha, cosbeta, sinbeta, cosgamma, singamma

  cosalpha = cos(angles(1))
  sinalpha = sin(angles(1))
  cosbeta  = cos(angles(2))
  sinbeta  = sin(angles(2))
  cosgamma = cos(angles(3))
  singamma = sin(angles(3))

  ! col 1:
  rotation_matrix(:,1) = (/ cosbeta*cosgamma, cosalpha*singamma + sinalpha*sinbeta*cosgamma, &
                                                                  & sinalpha*singamma - cosalpha*sinbeta*cosgamma /)
  ! col 2:
  rotation_matrix(:,2) = (/ -cosbeta*singamma, cosalpha*cosgamma - sinalpha*sinbeta*singamma, &
                                                                  & sinalpha*cosgamma + cosalpha*sinbeta*singamma /)
  ! col 3:
  rotation_matrix(:,3) = (/ sinbeta, -sinalpha*cosbeta, cosalpha*cosbeta /)

  ! perform rotation
  algor_3d_rotation = matmul(rotation_matrix, vector)
end function algor_3d_rotation

function algor_exp(x)
  implicit none
  real(kind=dp) , intent(in)  ::  x
  real(kind=dp)               ::  algor_exp

  if (x .gt. -500.0_dp) then
    algor_exp = exp(x)
  else
    algor_exp = 0.0_dp
  end if
end function algor_exp

subroutine algor_calc_inverse_real(array)
  use io, only: io_err
  implicit none
  real(kind=dp),  dimension(:,:), intent(inout) ::  array

  ! lapack vars:
  integer,        dimension(:), allocatable ::  ipiv
  real(kind=dp),  dimension(:), allocatable ::  work
  integer :: m, n, lwork, info

  m = size(array, 1)
  n = size(array, 2)

  allocate(ipiv(min(m,n)), stat=info)
  if (info .ne. 0) call io_err("algor_calc_inverse_real: Could not allocate ipiv array")

  call dgetrf(m, n, array, m, ipiv, info)
  if (info .ne. 0) call io_err("algor_calc_inverse_real: dgetrf failed")


  n = min(m,n)
  lwork = n

  ! can calc optimum lwork - ignore for now
  allocate(work(lwork), stat=info)
  if (info .ne. 0) call io_err("algor_calc_inverse_real: Could not allocate work array")

  call dgetri(n, array, n, ipiv, work, lwork, info)
  if (info .ne. 0) call io_err("algor_calc_inverse_real: dgetri failed")

  deallocate(ipiv, stat=info)
  if (info .ne. 0) call io_err("algor_calc_inverse_real: Could not deallocate ipiv array")

  deallocate(work, stat=info)
  if (info .ne. 0) call io_err("algor_calc_inverse_real: Could not deallocate work array")
end subroutine algor_calc_inverse_real

function algor_kronecker_delta(i,j)
    implicit none
    integer,  intent(in)  ::  i, j
    integer               ::  algor_kronecker_delta

    if (i .ne. j) then
      algor_kronecker_delta = 0
    else
      algor_kronecker_delta = 1
    end if
end function algor_kronecker_delta
end module algor
