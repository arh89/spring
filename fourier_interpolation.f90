module fourier_interpolation

  use constants,  only: dp, two_pi, cmplx_i
  implicit none

  private

  type, public  ::  interp_type
    integer                                       ::  ndims = 0
    integer,          dimension(2)                ::  grid_size = (/ 0, 0 /)
    real(kind=dp),    dimension(:,:), allocatable ::  grid_values
    complex(kind=dp), dimension(:,:), allocatable ::  fourier_coeffs
    logical                                       ::  ready = .false.
  end type interp_type

  ! likely won't need 1D routines, but might as well keep anyway
  interface dft_forward
    module procedure dft_forward_1d
    module procedure dft_forward_2d
  end interface dft_forward

  interface dft_backward
    module procedure dft_backward_1d
    module procedure dft_backward_2d
  end interface dft_backward

  interface fourier_interpolation_potential
    module procedure fourier_interpolation_potential_1d
    module procedure fourier_interpolation_potential_2d
  end interface fourier_interpolation_potential

  public  ::  fourier_interpolation_init
  public  ::  fourier_interpolation_read
  public  ::  fourier_interpolation_transform
  public  ::  fourier_interpolation_potential
  public  ::  fourier_interpolation_deallocate

contains

subroutine fourier_interpolation_init(keyword, interp, dealloc_real)
  use io, only: io_err
  implicit none
  character(len=*),             intent(in)      ::  keyword
  type(interp_type),            intent(inout)   ::  interp
  logical,            optional, intent(in)      ::  dealloc_real    ! default false
  integer ::  istat

  ! read real values at grid points from IO internal hash table
  call fourier_interpolation_read(keyword, interp)

  call fourier_interpolation_transform(interp)

  ! only deallocate if we explicitly ask to (safest option)
  if (present(dealloc_real) .and. dealloc_real) then
    deallocate(interp%grid_values, stat=istat)
    if (istat .ne. 0) call io_err("fourier_interpolation_init: Could not deallocate interp%grid_values")
  end if
end subroutine fourier_interpolation_init

subroutine fourier_interpolation_transform(interp)
  use io, only: io_err
  implicit none
  type(interp_type),  intent(inout)   ::  interp

  if (.not. allocated(interp%grid_values)) then
    call io_err("fourier_interpolation_transform: interp%grid_values array not allocated")
  end if

  ! forward transform - manually select between 1D and 2D and pass a slice if necessary
  if (interp%ndims .eq. 1) then
    call dft_forward_1d(interp%grid_values(:,0), interp%fourier_coeffs(:,0))
  else if (interp%ndims .eq. 2) then
    call dft_forward_2d(interp%grid_values, interp%fourier_coeffs)
  else
    call io_err("fourier_interpolation_transform: interp%ndims should be 1 or 2")
  end if
  interp%ready = .true.
end subroutine fourier_interpolation_transform

subroutine fourier_interpolation_allocate(interp)
  use io, only: io_err
  implicit none
  type(interp_type),            intent(inout)   ::  interp
  integer ::  istat

  if ((interp%grid_size(1) .lt. 1) .or. (interp%grid_size(2) .lt. 1)) then
    call io_err("fourier_interpolation_allocate: grid is not >= 1 along all dimensions")
  end if

  allocate(interp%grid_values(0:interp%grid_size(1)-1, 0:interp%grid_size(2)-1), stat=istat)
  if (istat .ne. 0) call io_err("fourier_interpolation_allocate: Could not allocate interp%grid_values")

  allocate(interp%fourier_coeffs(0:interp%grid_size(1)-1, 0:interp%grid_size(2)-1), stat=istat)
  if (istat .ne. 0) call io_err("fourier_interpolation_allocate: Could not allocate interp%fourier_coeffs")
end subroutine fourier_interpolation_allocate

subroutine fourier_interpolation_deallocate(interp)
  use io, only: io_err
  implicit none
  type(interp_type),            intent(inout)   ::  interp
  integer ::  istat

  ! might have deallocated this already..
  if (allocated(interp%grid_values)) then
    deallocate(interp%grid_values, stat=istat)
    if (istat .ne. 0) call io_err("fourier_interpolation_deallocate: Could not deallocate interp%grid_values")
  end if

  ! no reason why this should not be allocated if we are calling this, but let's be safe (could throw error)
  if (allocated(interp%fourier_coeffs)) then
    deallocate(interp%fourier_coeffs, stat=istat)
    if (istat .ne. 0) call io_err("fourier_interpolation_deallocate: Could not deallocate interp%fourier_coeffs")
  end if

  ! why not?
  interp%grid_size(:) = 0
  interp%ndims = 0
  interp%ready = .false.
end subroutine fourier_interpolation_deallocate

subroutine fourier_interpolation_read(keyword, interp)
  use io, only: max_line_len, io_input_get_data,  io_str_get_token, io_str_get_num_tokens, io_err
  implicit none
  character(len=*),             intent(in)    ::  keyword
  type(interp_type),            intent(inout) ::  interp

  character(len=max_line_len),  allocatable,  dimension(:)  ::  values
  character(len=max_line_len) ::  str
  character(len=10)           ::  line_num_str            ! any possible line num..
  logical                     ::  found
  integer,  dimension(2)      ::  grid_indices
  integer                     ::  itoken, iline, nlines
  integer                     ::  istat

  ! Error checking kept fairly basic here:
  !   - Don't check whether value supplied more than once
  !   - Don't check that all values are supplied.
  ! (Benefit: Don't enforce any particular ordering of points)

  ! Block should look like:
  !   Note: values are: int, [int,] real (for all lines after first)

  ! for 1D:
  ! npts
  ! 0           value(1)
  ! npts-1      value(npts)

  ! for 2D:
  ! na      nb
  ! 0       0       value(1)
  ! ...     ...     ...
  ! na-1    nb-1    value(na*nb)

  ! don't need to worry about length of keyword, io module will complain if too long (as does hash..)
  call io_input_get_data(keyword, values, found)
  if (.not. found) call io_err("fourier_interpolation_read: "//trim(keyword)//" not found in input hash")

  nlines = size(values, 1)

  ! must have at least two lines
  if (nlines .lt. 2) call io_err("fourier_interpolation_read: "//trim(keyword)//" block must be >= 2 lines long")

  ! check first line has correct number of tokens (1 for 1D,   2 for 2D)
  interp%ndims = io_str_get_num_tokens(values(1))
  if ((interp%ndims .ne. 1) .and. (interp%ndims .ne. 2))  &
    & call io_err("fourier_interpolation_read: First line of "//trim(keyword)//" block must have 1 or 2 values (ndims)")

  ! get number of grid points
  if (interp%ndims .eq. 1) then
    read(values(1), *, iostat=istat) interp%grid_size(1)
    interp%grid_size(2) = 1
  else
    ! interp%ndims .eq. 2
    read(values(1), *, iostat=istat) interp%grid_size(1), interp%grid_size(2)
  end if
  if (istat .ne. 0) call io_err("fourier_interpolation_read: "//trim(keyword)//" block: could not read number of grid points")


  ! check number of lines makes sense
  if (nlines-1 .ne. interp%grid_size(1)*interp%grid_size(2)) &
    & call io_err("fourier_interpolation_read: "//trim(keyword)//" block has unexpected number of grid points")

  call fourier_interpolation_allocate(interp) ! will check we have valid number of grid points along each dim

  interp%grid_values(:,:) = 0.0_dp            ! prevent doing DFT of uninitialized data if grid points missing
                                              ! (due setting one twice, for example) - result will likely be wrong though..

  ! values can be in any order, but must all be there for this to work properly
  ! Note, could create logical array of size Nx*Ny and check off each element when added, then look to see if any not set..
  do iline = 2, nlines
    ! store line num as string for useful error messages..
    write(line_num_str, '(I10)', iostat=istat) iline
    if (istat .ne. 0) call io_err("fourier_interpolation_read: Could not store line number string")
    line_num_str = adjustl(line_num_str)  ! remove lhs padding

    if (io_str_get_num_tokens(values(iline)) .ne. interp%ndims+1) then
      call io_err("fourier_interpolation_read: line "//trim(line_num_str) &
                  & //" in "//trim(keyword)//" block must have ndims+1 values")
    end if

    grid_indices(:) = 0
    ! check that grid indices are valid (ie: >=0, <na or <nb)
    do itoken = 1, interp%ndims    ! ia, ib

      str = io_str_get_token(values(iline), itoken)
      read(str, *, iostat=istat) grid_indices(itoken)
      if (istat .ne. 0) then
        call io_err("fourier_interpolation_read: "//trim(keyword)//&
          & " could not read grid point index on line "//trim(line_num_str)//" of block")
      end if

      if ((grid_indices(itoken) .lt. 0) .or. (grid_indices(itoken) .ge. interp%grid_size(itoken))) then
        call io_err("fourier_interpolation_read: "//trim(keyword)//" invalid grid point index on line "&
          &//trim(line_num_str)//" of block")
      end if
    end do

    str = io_str_get_token(values(iline), interp%ndims+1)
    read(str, *, iostat=istat) interp%grid_values(grid_indices(1), grid_indices(2))
    if (istat .ne. 0) then
      call io_err("fourier_interpolation_read: "//trim(keyword)//" could not read grid point on line "&
        &//trim(line_num_str)//" of block")
    end if
  end do
end subroutine fourier_interpolation_read

subroutine dft_forward_1d(rinput, coutput)
  use io, only: io_err
  implicit none
  real(kind=dp),    dimension(0:), intent(in)  ::  rinput
  complex(kind=dp), dimension(0:), intent(out) ::  coutput
  ! local vars:
  integer               ::  n, k, nelements

  ! computes the 1D DFT (see numerical recipes)

  if (size(rinput,1) .ne. size(coutput,1)) call io_err("dft_forward_1d: Size of two arrays should match")

  nelements = size(rinput,1)

  coutput(:) = cmplx(0.0_dp, 0.0_dp, dp)

  do n = 0, nelements-1
    do k = 0, nelements-1
      coutput(n) = coutput(n) + rinput(k) * exp( (two_pi*real(k*n,dp)/real(nelements,dp)) * cmplx_i )
    end do !k
  end do !n
end subroutine dft_forward_1d

subroutine dft_forward_2d(rinput, coutput)
  use io, only: io_err
  implicit none
  real(kind=dp),    dimension(0:,0:), intent(in)  ::  rinput
  complex(kind=dp), dimension(0:,0:), intent(out) ::  coutput
  ! local vars:
  complex(kind=dp)      ::  contrib
  integer               ::  n1, n2, k1, k2, nelements1, nelements2

  ! computes the 2D DFT (see numerical recipes)
  ! could do this with nested 1D transforms (to save code repetition) but clearer to have it all in one routine

  if ((size(rinput,1) .ne. size(coutput,1)) .or. &
  &   (size(rinput,2) .ne. size(coutput,2))) call io_err("dft_forward_2d: Size of two arrays should match")

  nelements1 = size(rinput,1)
  nelements2 = size(rinput,2)

  coutput(:,:) = cmplx(0.0_dp, 0.0_dp, dp)

  do n2 = 0, nelements2-1
    do n1 = 0, nelements1-1

      do k2 = 0, nelements2-1
        contrib = exp( (two_pi*real(k2*n2,dp)/real(nelements2,dp)) * cmplx_i )

        do k1 = 0, nelements1-1
          coutput(n1,n2) = coutput(n1,n2) + rinput(k1,k2) &
          &                    * contrib * exp( (two_pi*real(k1*n1,dp)/real(nelements1,dp)) * cmplx_i )
        end do !k1
      end do !k2

    end do !n1
  end do !n2
end subroutine dft_forward_2d

subroutine dft_backward_1d(cinput, routput)
  use io, only: io_err
  implicit none
  complex(kind=dp), dimension(0:), intent(in)   ::  cinput
  real(kind=dp),    dimension(0:), intent(out)  ::  routput
  ! local vars:
  complex(kind=dp)  ::  coutput
  integer           ::  n, k, nelements

  ! computes the 1D inverse DFT (see numerical recipes)

  if (size(cinput,1) .ne. size(routput,1)) call io_err("dft_backward_1d: Size of two arrays should match")

  nelements = size(cinput,1)

  do k = 0, nelements-1
    coutput = cmplx(0.0_dp, 0.0_dp, dp)

    do n = 0, nelements-1
      coutput = coutput + cinput(n) * exp( (-two_pi*real(k*n,dp)/real(nelements,dp)) * cmplx_i )
    end do !n

    routput(k) = real(coutput,dp)
  end do !k
  routput = routput/real(nelements,dp)
end subroutine dft_backward_1d

subroutine dft_backward_2d(cinput, routput)
  use io, only: io_err
  implicit none
  complex(kind=dp), dimension(0:,0:), intent(in)  ::  cinput
  real(kind=dp),    dimension(0:,0:), intent(out) ::  routput
  ! local vars:
  complex(kind=dp)      ::  coutput, contrib
  integer               ::  n1, n2, k1, k2, nelements1, nelements2

  ! computes the 2D inverse DFT (see numerical recipes)
  ! could do this with nested 1D transforms (to save code repetition) but clearer to have it all in one routine

  if ((size(cinput,1) .ne. size(routput,1)) .or. &
  &   (size(cinput,2) .ne. size(routput,2))) call io_err("dft_backward_2d: Size of two arrays should match")

  nelements1 = size(cinput,1)
  nelements2 = size(cinput,2)


  do k2 = 0, nelements2-1
    do k1 = 0, nelements1-1
      coutput = cmplx(0.0_dp, 0.0_dp, dp)

      do n2 = 0, nelements2-1
        contrib = exp( (-two_pi*real(k2*n2,dp)/real(nelements2,dp)) * cmplx_i )

        do n1 = 0, nelements1-1
          coutput = coutput + cinput(n1,n2) &
          &                    * contrib * exp( (-two_pi*real(k1*n1,dp)/real(nelements1,dp)) * cmplx_i )
        end do !n1
      end do !n2
      routput(k1,k2) = real(coutput,dp)

    end do !k1
  end do !k2
  routput(:,:) = routput(:,:)/real(nelements1*nelements2,dp)
end subroutine dft_backward_2d

! wrapper routines (automatically pass correct slice of array)
subroutine fourier_interpolation_potential_1d(interp, pos, routput, rgrad)
  use io, only: io_err
  implicit none
  type(interp_type),            intent(in)  ::  interp
  real(kind=dp),                intent(in)  ::  pos
  real(kind=dp),                intent(out) ::  routput
  real(kind=dp),      optional, intent(out) ::  rgrad

  if ((interp%ndims .ne. 1) .or. (.not. interp%ready)) then
    call io_err("fourier_interpolation_potential_1d: interp type does not contain 1D data/has not been transformed")
  end if

  if (present(rgrad)) then
    call fourier_interp_pot_1d(interp%fourier_coeffs(:,0), pos, routput, rgrad)
  else
    call fourier_interp_pot_1d(interp%fourier_coeffs(:,0), pos, routput)
  end if
end subroutine fourier_interpolation_potential_1d

subroutine fourier_interpolation_potential_2d(interp, pos, routput, rgrad)
  use io, only: io_err
  implicit none
  type(interp_type),                          intent(in)  ::  interp
  real(kind=dp),      dimension(2),           intent(in)  ::  pos
  real(kind=dp),                              intent(out) ::  routput
  real(kind=dp),      dimension(2), optional, intent(out) ::  rgrad

  if ((interp%ndims .ne. 2) .or. (.not. interp%ready)) then
    call io_err("fourier_interpolation_potential_2d: interp type does not contain 2D data/has not been transformed")
  end if

  if (present(rgrad)) then
    call fourier_interp_pot_2d(interp%fourier_coeffs, pos, routput, rgrad)
  else
    call fourier_interp_pot_2d(interp%fourier_coeffs, pos, routput)
  end if
end subroutine fourier_interpolation_potential_2d

! 'short name' routines are lower level... called by wrapper routines
subroutine fourier_interp_pot_1d(cinput, pos, routput, rgrad)
  implicit none
  complex(kind=dp), dimension(0:),  intent(in)  ::  cinput
  real(kind=dp),                    intent(in)  ::  pos
  real(kind=dp),                    intent(out) ::  routput
  real(kind=dp),    optional,       intent(out) ::  rgrad
  ! local vars:
  complex(kind=dp)  ::  cgrad
  complex(kind=dp)  ::  component, coutput
  integer           ::  n, nelements

  ! outputs inverse fourier transform of data set interpolated at a pos
  nelements = size(cinput,1)

  coutput = cmplx(0.0_dp, 0.0_dp, dp)
  cgrad = cmplx(0.0_dp, 0.0_dp, dp)

  ! first half of n
  do n = 0, nelements/2
    component = cinput(n)*exp(-two_pi*real(n,dp)*pos*cmplx_i)
    coutput = coutput + component
    cgrad = cgrad + component*real(n,dp)
  end do

  ! second half of n
  do n = nelements/2+1, nelements-1
    component = cinput(n)*exp(-two_pi*real(n-nelements,dp)*pos*cmplx_i)
    coutput = coutput + component
    cgrad = cgrad + component*real(n-nelements,dp)
  end do

  ! only do the rest of the work with cgrad if we need to
  if (present(rgrad)) then
    cgrad = cgrad*(-cmplx_i)
    rgrad = two_pi*real(cgrad,dp)/real(nelements,dp)
  end if

  routput = real(coutput,dp)/real(nelements,dp)
end subroutine fourier_interp_pot_1d

subroutine fourier_interp_pot_2d(cinput, pos, routput, rgrad)
  implicit none
  complex(kind=dp), dimension(0:,0:),           intent(in)  ::  cinput
  real(kind=dp),    dimension(2),               intent(in)  ::  pos
  real(kind=dp),                                intent(out) ::  routput
  real(kind=dp),    dimension(2),     optional, intent(out) ::  rgrad
  ! local vars:
  complex(kind=dp), dimension(2)  ::  cgrad
  complex(kind=dp)                ::  contrib
  complex(kind=dp)                ::  component, coutput
  integer   ::  n1, n2, nelements1, nelements2

  ! outputs inverse fourier transform of data set interpolated at a pos
  nelements1 = size(cinput,1)
  nelements2 = size(cinput,2)

  coutput = cmplx(0.0_dp, 0.0_dp, dp)
  cgrad = cmplx(0.0_dp, 0.0_dp, dp)

  ! first half of n2
  do n2 = 0, nelements2/2

    contrib = exp(-two_pi*real(n2,dp)*pos(2)*cmplx_i)

    ! first half of n1
    do n1 = 0, nelements1/2
      component = cinput(n1,n2)*exp(-two_pi*real(n1,dp)*pos(1)*cmplx_i)*contrib
      coutput = coutput + component
      cgrad(1) = cgrad(1) + component*real(n1,dp)
      cgrad(2) = cgrad(2) + component*real(n2,dp)
    end do !n1

    ! second half of n2
    do n1 = nelements1/2+1, nelements1-1
      component = cinput(n1,n2)*exp(-two_pi*real(n1-nelements1,dp)*pos(1)*cmplx_i)*contrib
      coutput = coutput + component
      cgrad(1) = cgrad(1) + component*real(n1-nelements1,dp)
      cgrad(2) = cgrad(2) + component*real(n2,dp)
    end do !n1
  end do !n2

  ! second half of n2
  do n2 = nelements2/2+1, nelements2-1

    contrib = exp(-two_pi*real(n2-nelements2,dp)*pos(2)*cmplx_i)

    ! first half of n1
    do n1 = 0, nelements1/2
      component = cinput(n1,n2)*exp(-two_pi*real(n1,dp)*pos(1)*cmplx_i)*contrib
      coutput = coutput + component
      cgrad(1) = cgrad(1) + component*real(n1,dp)
      cgrad(2) = cgrad(2) + component*real(n2-nelements2,dp)
    end do !n1

    ! second half of n2
    do n1 = nelements1/2+1, nelements1-1
      component = cinput(n1,n2)*exp(-two_pi*real(n1-nelements1,dp)*pos(1)*cmplx_i)*contrib
      coutput = coutput + component
      cgrad(1) = cgrad(1) + component*real(n1-nelements1,dp)
      cgrad(2) = cgrad(2) + component*real(n2-nelements2,dp)
    end do !n1
  end do !n2

  ! only do the rest of the work with cgrad if we need to
  if (present(rgrad)) then
    cgrad(:) = cgrad(:)*(-cmplx_i)
    rgrad(:) = two_pi*real(cgrad,dp)/real(nelements2*nelements1,dp)
  end if

  routput = real(coutput,dp)/real(nelements2*nelements1,dp)
end subroutine fourier_interp_pot_2d

end module fourier_interpolation
