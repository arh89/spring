module pes_scan
  ! Can either: 
  !     * Raster along a surface that is in the a,b plane to get minimum energy and corresponding height along c
  !     * Raster along a line across the a,b plane to get minimum energy and corresponding height along c
  !     * Get minimum energy and corresponding height along c at a point on the a,b plane

  use constants,  only: dp
  implicit none
  
  private

  ! TSPD parameters (for opt of z coord)
  integer,          parameter ::  default_max_iter = 100
  real(kind=dp),    parameter ::  default_initial_alpha = 1.0E-05_dp
  real(kind=dp),    parameter ::  default_max_disp  = 5.0E-1_dp       ! Ang
  real(kind=dp),    parameter ::  default_force_tol = 1.0E-5_dp       ! eV/Ang 
  real(kind=dp),    parameter ::  default_energy_tol = 1.0E-8_dp      ! eV
  real(kind=dp),    parameter ::  default_disp_tol = 1.0E-3_dp        ! Ang 

  ! line scan parameters:
  real(kind=dp),    parameter ::  default_ls_height = 3.0_dp          ! Ang
  real(kind=dp),    parameter ::  default_ls_dx = 1.0E-1_dp           ! Ang

  character(len=4), parameter ::  default_opt_method = 'tpsd'

  type, public  :: pes_scan_type
    integer                       ::  file_unit                       ! IO unit
    integer                       ::  probe_id                        ! id of probe atom within cell
    integer,        dimension(2)  ::  grid_points                     ! number of grid points along a,b
    real(kind=dp),  dimension(2)  ::  dx_real                         ! dx along the real length of a,b
    real(kind=dp),  dimension(2)  ::  dx_frac                         ! dx in the fractional space along a,b
    real(kind=dp),  dimension(3)  ::  lat_vec_len                     ! lengths of a,b,c vectors
    integer,        dimension(3)  ::  supercell                       ! supercell size

    ! TPSD:
    integer                       ::  max_iter
    real(kind=dp)                 ::  initial_alpha
    real(kind=dp)                 ::  force_tol
    real(kind=dp)                 ::  energy_tol
    real(kind=dp)                 ::  disp_tol
    real(kind=dp)                 ::  max_disp

    ! line scan:
    real(kind=dp)                 ::  line_scan_height
    real(kind=dp)                 ::  line_scan_dx

    ! optimization method (either 'tpsd' or 'scan')
    character(len=4)              ::  opt_method
  end type pes_scan_type

  public  ::  pes_scan_init
  public  ::  do_pes_scan

contains

subroutine pes_scan_init(pes, cell)
  use cell,       only: cell_type, cell_reallocate, cell_shift_to_unit_cell, cell_frac_to_cart, cell_supercell 
  use io,         only: io_input_get_single_value, max_line_len, io_input_get_data, io_str_get_num_tokens, &
                      & io_str_get_token, io_str_to_real, io_str_to_int, io_err
  use constants,  only: element_symbol_to_z, units_natural_to_atomic, element_mass

  implicit none
  type(pes_scan_type),  intent(inout) ::  pes
  type(cell_type),      intent(inout) ::  cell

  character(len=max_line_len)                                 ::  keyword
  character(len=max_line_len),  allocatable,  dimension(:)    ::  input_data
  real(kind=dp),                allocatable,  dimension(:,:)  ::  tmp_atom_frac_pos
  real(kind=dp),                allocatable,  dimension(:)    ::  tmp_atom_mass
  integer,                      allocatable,  dimension(:)    ::  tmp_atom_species
  real(kind=dp) ::  vacuum_size
  logical       ::  kw_found
  integer       ::  iatom, probe_species, itoken, ntokens, istat

  ! routine only does some very basic error checking..

  ! get num grid points from input
  keyword = 'pes_grid_points'
  call io_input_get_data(keyword, input_data, kw_found)
  if (.not. kw_found) call io_err('pes_scan_init: '//trim(keyword)//' not found in input file')

  if (size(input_data,1) .ne. 1) call io_err('pes_scan_init: '//trim(keyword)//' should be a single line')

  ntokens = io_str_get_num_tokens(input_data(1))
  if (ntokens .ne. 2) call io_err('pes_scan_init: '//trim(keyword)//' expects 2 values')

  read(unit=input_data(1), fmt=*, iostat=istat) pes%grid_points(1), pes%grid_points(2)
  if (istat .ne. 0) call io_err('pes_scan_init: Error reading '//trim(keyword)//' values')

  if ((pes%grid_points(1) .lt. 1) .or. (pes%grid_points(2) .lt. 1)) &
    &   call io_err('pes_scan_init: Number of grid points along each dimension must be >= 1')

  ! set max iter
  call io_input_get_single_value('pes_max_iter', pes%max_iter, default_max_iter)
  if (pes%max_iter .lt. 1) call io_err('pes_scan_init: pes_max_iter should be > 0')

  ! set force tol
  call io_input_get_single_value('pes_force_tol', pes%force_tol, default_force_tol)
  if (pes%force_tol .le. 0.0_dp) call io_err('pes_scan_init: pes_force_tol should be > 0.0')
  pes%force_tol = units_natural_to_atomic(force=pes%force_tol)

  ! set energy tol
  call io_input_get_single_value('pes_energy_tol', pes%energy_tol, default_energy_tol)
  if (pes%energy_tol .le. 0.0_dp) call io_err('pes_scan_init: pes_energy_tol should be > 0.0')
  pes%energy_tol = units_natural_to_atomic(energy=pes%energy_tol)

  ! set disp tol
  call io_input_get_single_value('pes_disp_tol', pes%disp_tol, default_disp_tol)
  if (pes%disp_tol .le. 0.0_dp) call io_err('pes_scan_init: pes_disp_tol should be > 0.0')
  pes%disp_tol = units_natural_to_atomic(length=pes%disp_tol)

  ! set max displacement
  call io_input_get_single_value('pes_max_disp', pes%max_disp, default_max_disp)
  if (pes%disp_tol .le. 0.0_dp) call io_err('pes_scan_init: pes_max_disp should be > 0.0')
  pes%max_disp = units_natural_to_atomic(length=pes%max_disp)

  ! set initial alpha
  call io_input_get_single_value('pes_initial_alpha', pes%initial_alpha, default_initial_alpha)
  if (pes%initial_alpha .le. 0.0_dp) call io_err('pes_scan_init: pes_initial_alpha should be > 0.0')

  ! set line_scan_height
  call io_input_get_single_value('pes_line_scan_height', pes%line_scan_height, default_ls_height)
  if (pes%line_scan_height .le. 0.0_dp) call io_err('pes_scan_init: pes_line_scan_height should be > 0.0')
  pes%line_scan_height = units_natural_to_atomic(length=pes%line_scan_height)

  ! set line_scan_dx
  call io_input_get_single_value('pes_line_scan_dx', pes%line_scan_dx, default_ls_dx)
  if (pes%line_scan_dx .lt. 0.0_dp) call io_err('pes_scan_init: pes_line_scan_dx should be >= 0.0')
  pes%line_scan_dx = units_natural_to_atomic(length=pes%line_scan_dx)

  ! set optimization method
  call io_input_get_single_value('pes_opt_method', pes%opt_method, default_opt_method)
  if ((pes%opt_method .ne. 'tpsd') .and. (pes%opt_method .ne. 'scan')) then
    call io_err('pes_scan_init: pes_opt_method must be either "tpsd" or "scan"')
  end if

  ! make supercell
  keyword = 'pes_supercell'
  call io_input_get_data(keyword, input_data, kw_found)
  if (.not. kw_found) then
    pes%supercell(:) = 1
  else
    if (size(input_data,1) .ne. 1) call io_err('pes_scan_init: '//trim(keyword)//' should be a single line')
    ntokens = io_str_get_num_tokens(input_data(1))
    if (ntokens .ne. 3) call io_err('pes_scan_init: '//trim(keyword)//' expects 3 values')

    do itoken = 1, 3
      pes%supercell(itoken) = io_str_to_int(io_str_get_token(input_data(1),itoken))
      if (pes%supercell(itoken) .lt. 1) call io_err('pes_scan_init: Supercell size should be >= 1')
    end do

    ! not 1x1x1 supercell
    if (sum(pes%supercell(:)) .ne. 3) then
      call cell_supercell(cell, pes%supercell(1), pes%supercell(2), pes%supercell(3))
    end if
  end if

  ! must be done *after* creating supercell if necessary
  ! get probe from input (frac coords)
  keyword = 'pes_probe'
  call io_input_get_data(trim(keyword), input_data, kw_found)
  if (.not. kw_found) call io_err('pes_scan_init: '//trim(keyword)//' missing in input file')

  if (size(input_data, 1) .ne. 1) call io_err('pes_scan_init: Input must contain *one* probe atom')

  ! allocate tmp arrays to store some cell info
  allocate(tmp_atom_frac_pos(3, cell%natoms), stat=istat)
  if (istat .ne. 0) call io_err('pes_scan_init: Could not allocate tmp_atom_frac_pos array')

  allocate(tmp_atom_mass(cell%natoms), stat=istat)
  if (istat .ne. 0) call io_err('pes_scan_init: Could not allocate tmp_atom_mass array')

  allocate(tmp_atom_species(cell%natoms), stat=istat)
  if (istat .ne. 0) call io_err('pes_scan_init: Could not allocate tmp_atom_species array')

  ! copy to tmp arrays
  tmp_atom_frac_pos(:,:) = cell%atom_frac_pos(:,:)
  tmp_atom_mass(:) = cell%atom_mass(:)
  tmp_atom_species(:) = cell%atom_species(:)

  ! reallocate cell arrays to new size
  cell%natoms = cell%natoms + 1
  pes%probe_id = cell%natoms
  call cell_reallocate(cell)

  ! copy old data back into cell
  do iatom = 1, cell%natoms-1
    cell%atom_frac_pos(:,iatom) = tmp_atom_frac_pos(:,iatom)
    cell%atom_mass(iatom) = tmp_atom_mass(iatom)
    cell%atom_species(iatom) = tmp_atom_species(iatom)
  end do

  ! deallocate tmp arrays
  deallocate(tmp_atom_frac_pos, stat=istat)
  if (istat .ne. 0) call io_err('pes_scan_init: Could not deallocate tmp_atom_frac_pos array')

  deallocate(tmp_atom_mass, stat=istat)
  if (istat .ne. 0) call io_err('pes_scan_init: Could not deallocate tmp_atom_mass array')

  deallocate(tmp_atom_species, stat=istat)
  if (istat .ne. 0) call io_err('pes_scan_init: Could not deallocate tmp_atom_species array')

  ! add probe
  ntokens = io_str_get_num_tokens(input_data(1))
  if (ntokens .ne. 4) call io_err("pes_scan_init: Expected 4 values per line for "//trim(keyword)//" block")

  ! first token is atom species
  probe_species = element_symbol_to_z(trim(io_str_get_token(input_data(1), 1)))
  if (probe_species .lt. 1) call io_err("pes_scan_init: Unexpected atom species")

  cell%atom_species(pes%probe_id) = probe_species

  do itoken = 2, ntokens !2,4
    ! itoken-1 is coordinate index, iatom = iline
    cell%atom_frac_pos(itoken-1, pes%probe_id) = io_str_to_real(io_str_get_token(input_data(1), itoken))
  end do

  ! set mass of the probe:
  cell%atom_mass(pes%probe_id) = units_natural_to_atomic(mass=element_mass(probe_species))

  ! apply pbcs, get cart positions
  call cell_shift_to_unit_cell(cell)
  call cell_frac_to_cart(cell)

  ! input_data should always be allocated:
  if (allocated(input_data)) then
    deallocate(input_data, stat=istat)
    if (istat .ne. 0) call io_err('pes_scan_init: Could not deallocate input_data array')
  end if

  ! lengths of a, b and c
  pes%lat_vec_len(1) = sqrt(dot_product(cell%lattice_vectors(:,1),cell%lattice_vectors(:,1)))
  pes%lat_vec_len(2) = sqrt(dot_product(cell%lattice_vectors(:,2),cell%lattice_vectors(:,2)))
  pes%lat_vec_len(3) = sqrt(dot_product(cell%lattice_vectors(:,3),cell%lattice_vectors(:,3)))

  ! calculate dx_frac and dx_real
  ! (smallest possible value of grid_points is 1, don't worry about division by zero)
  pes%dx_frac(1) = 1.0_dp/real(pes%grid_points(1), dp)
  pes%dx_frac(2) = 1.0_dp/real(pes%grid_points(2), dp)

  ! take into account supercell:
  pes%dx_frac(1) = pes%dx_frac(1)/real(pes%supercell(1),dp)
  pes%dx_frac(2) = pes%dx_frac(2)/real(pes%supercell(2),dp)

  pes%dx_real(1) = pes%dx_frac(1)*pes%lat_vec_len(1)
  pes%dx_real(2) = pes%dx_frac(2)*pes%lat_vec_len(2)

  ! check that line_scan_height doesn't cause probe to go past upper cell boundary
  ! (assume surface in xy plane and optimization is along z...)
  if (pes%opt_method .eq. 'scan') then
    vacuum_size = (1.0_dp - cell%atom_frac_pos(3,pes%probe_id))*pes%lat_vec_len(3)
    if (vacuum_size .lt. pes%line_scan_height) then
      call io_err("pes_scan_init: Probe height will exceed cell boundary -" &
      & //" reduce pes_line_scan_height or increase vacuum region")
    end if
  end if
end subroutine pes_scan_init

subroutine do_pes_scan(pes, cell)
  use cell,       only: cell_type, cell_shift_to_unit_cell, cell_frac_to_cart
  use io,         only: io_err, io_open_file, io_close_file, seedname
  use constants,  only: units_atomic_to_natural

  implicit none

  type(pes_scan_type),  intent(inout) ::  pes
  type(cell_type),      intent(inout) ::  cell

  real(kind=dp),  dimension(:,:), allocatable ::  forces
  real(kind=dp),  dimension(3)    ::  probe_initial_frac_pos
  real(kind=dp)                   ::  energy
  integer                         ::  ia, ib
  integer                         ::  istat
  integer                         ::  num_iters

  call io_open_file(trim(seedname)//'.pes', pes%file_unit, 'replace')

  probe_initial_frac_pos(:) = cell%atom_frac_pos(:, pes%probe_id)

  allocate(forces(3,cell%natoms), stat=istat)
  if (istat .ne. 0) call io_err("do_pes_scan: Could not allocate forces array")

!  write(pes%file_unit, *), 'grid points = ', pes%grid_points(1),  pes%grid_points(2)
!  write(pes%file_unit, *), 'dx_frac = ', pes%dx_frac(1),  pes%dx_frac(2)
!  write(pes%file_unit, *), 'dx_real = ', pes%dx_real(1),  pes%dx_real(2)

  ! raster over surface
  do ib = 0, pes%grid_points(2)-1

    ! probe position (b):
    cell%atom_frac_pos(2, pes%probe_id) = probe_initial_frac_pos(2) + real(ib,dp)*pes%dx_frac(2)

    do ia = 0, pes%grid_points(1)-1

      ! probe position (a and c):
      cell%atom_frac_pos(1, pes%probe_id) = probe_initial_frac_pos(1) + real(ia,dp)*pes%dx_frac(1)
      cell%atom_frac_pos(3, pes%probe_id) = probe_initial_frac_pos(3)

      ! resync cell after shifting probe
      call cell_shift_to_unit_cell(cell)
      call cell_frac_to_cart(cell)

      select case (pes%opt_method)
        case ('tpsd')
          call pes_tpsd_line_minimise(pes, cell, forces, energy, num_iters)
      
          ! write final iter to file
          write(pes%file_unit, *, iostat=istat) cell%atom_frac_pos(:,pes%probe_id), units_atomic_to_natural(energy=energy), &
          & units_atomic_to_natural(force=forces(:,pes%probe_id)), num_iters, '<-- final'
          if (istat .ne. 0) call io_err("do_pes_scan: Could not write to output file")

        case ('scan')
          call pes_line_scan(pes, cell, forces, energy)
        case default
          call io_err("do_pes_scan: Unexpected optimization method")
      end select
    end do
  end do

  call io_close_file(pes%file_unit)

  deallocate(forces, stat=istat)
  if (istat .ne. 0) call io_err("do_pes_scan: Could not deallocate forces array")
end subroutine do_pes_scan

subroutine pes_line_scan(pes, cell, forces, energy)
  use cell,       only: cell_type, cell_cart_to_frac, cell_frac_to_cart, cell_shift_to_unit_cell
  use potential,  only: pot_get_forces
  use constants,  only: units_atomic_to_natural
  use io,         only: io_err
  implicit none
  type(pes_scan_type),                        intent(in)    ::  pes
  type(cell_type),                            intent(inout) ::  cell
  real(kind=dp),  dimension(3, cell%natoms),  intent(out)   ::  forces
  real(kind=dp),                              intent(out)   ::  energy

  ! local variables:
  real(kind=dp) ::  initial_pos, pos
  integer       ::  istat

  ! initial positions and forces:
  initial_pos = cell%atom_cart_pos(3, pes%probe_id)
  pos = initial_pos
  call pot_get_forces(cell, forces, energy)

  write(pes%file_unit, *, iostat=istat) cell%atom_frac_pos(:,pes%probe_id), units_atomic_to_natural(energy=energy), &
  & units_atomic_to_natural(force=forces(:,pes%probe_id))
  if (istat .ne. 0) call io_err("do_pes_scan: Could not write to output file")

  do while (abs(pos-initial_pos) .lt. pes%line_scan_height)
    ! update cell:
    cell%atom_cart_pos(3, pes%probe_id) = pos + pes%line_scan_dx
    call cell_cart_to_frac(cell)                  ! keep everything in sync

    !  make sure PBCs are enforced correctly
    call cell_shift_to_unit_cell(cell)
    call cell_frac_to_cart(cell)
    pos = cell%atom_cart_pos(3, pes%probe_id)

    ! forces:
    call pot_get_forces(cell, forces, energy)

    ! write iter data:
    write(pes%file_unit, *, iostat=istat) cell%atom_frac_pos(:,pes%probe_id), units_atomic_to_natural(energy=energy), &
    & units_atomic_to_natural(force=forces(:,pes%probe_id))
    if (istat .ne. 0) call io_err("do_pes_scan: Could not write to output file")
  end do
end subroutine pes_line_scan

subroutine pes_tpsd_line_minimise(pes, cell, forces, energy, iter)
  use cell,       only: cell_type, cell_cart_to_frac, cell_frac_to_cart, cell_shift_to_unit_cell
  use potential,  only: pot_get_forces
  use constants,  only: units_atomic_to_natural
  use io,         only: io_err
  implicit none
  type(pes_scan_type),                        intent(in)    ::  pes
  type(cell_type),                            intent(inout) ::  cell
  real(kind=dp),  dimension(3, cell%natoms),  intent(out)   ::  forces
  real(kind=dp),                              intent(out)   ::  energy
  integer,                                    intent(out)   ::  iter

  ! local variables:
  logical       ::  force_conv, disp_conv, energy_conv, conv
  real(kind=dp) ::  delta_pos, delta_grad, pos, pos_old, energy_old, grad, grad_old, alpha
  integer       ::  istat

  ! Routine does 1D TPSD minimization along z (not necessarily c) for a probe atom on top of a fixed surface (in xy plane)
  ! TPSD algorithm:
  ! g_k = -force_k
  ! delta_x = x_k - x_k-1
  ! delta_g = g_k - g_k-1
  ! alpha_k = dot_product(delta_x, delta_g)/dot_product(delta_g, delta_g)
  ! x_k+1 = x_k - alpha_k*g_k

  ! We will "throw away" first step by picking small alpha - need to have a smart way of choosing this
  alpha = pes%initial_alpha
  iter = 0

  ! initial positions and forces:
  pos_old = cell%atom_cart_pos(3, pes%probe_id)
  call pot_get_forces(cell, forces, energy)
  grad_old = -forces(3, pes%probe_id)
  grad = grad_old
  energy_old = energy

  ! assume disp and energy are fine at first.. check only for small force at first iter
  force_conv = .false.
  disp_conv = .false.
  energy_conv = .false.
  conv = .false.

  if (abs(grad) .le. pes%force_tol) then
    force_conv =  .true.
    conv = .true.
  end if

  ! write initial info
  write(pes%file_unit, *, iostat=istat) cell%atom_frac_pos(:,pes%probe_id), units_atomic_to_natural(energy=energy), &
  & units_atomic_to_natural(force=forces(:,pes%probe_id)), alpha, iter
  if (istat .ne. 0) call io_err("do_pes_scan: Could not write to output file")

  ! only update position if we actually need to (initial guess might have been at minimum)
  do while (.not. conv)
    if (iter .eq. pes%max_iter) then
      ! print warning - could exit here..
      write(pes%file_unit, *, iostat=istat) '*** WARNING: Max number iterations reached ***'
      if (istat .ne. 0) call io_err("pes_tpsd_line_minimize: Could not write to output file")
      exit
    end if

    ! reset flags
    force_conv = .false.
    disp_conv = .false.
    energy_conv = .false.
    conv = .false.

    iter = iter + 1
    pos = pos_old - alpha*grad_old                ! update position

    ! allow for max disp
    if (abs(pos-pos_old) .gt. pes%max_disp) then
      write(pes%file_unit, *, iostat=istat) '*** WARNING: Displacement larger than pes%max_disp - moving by pes%max_disp ***'
      if (istat .ne. 0) call io_err("pes_tpsd_line_minimize: Could not write to output file")
      pos = pos_old + sign(pes%max_disp, (pos-pos_old))
    end if

    ! update cell:
    cell%atom_cart_pos(3, pes%probe_id) = pos
    call cell_cart_to_frac(cell)                  ! keep everything in sync

    !  make sure PBCs are enforced correctly
    call cell_shift_to_unit_cell(cell)
    call cell_frac_to_cart(cell)
    pos = cell%atom_cart_pos(3, pes%probe_id)

    ! forces:
    call pot_get_forces(cell, forces, energy)
    grad = -forces(3, pes%probe_id)

    if (abs(grad) .le. pes%force_tol) force_conv = .true.
    if (abs(pos-pos_old) .le. pes%disp_tol) disp_conv = .true.
    if (abs(energy-energy_old) .le. pes%energy_tol) energy_conv = .true.

    if (force_conv .and. disp_conv .and. energy_conv) conv = .true.

    ! write iter data:
    write(pes%file_unit, *, iostat=istat) cell%atom_frac_pos(:,pes%probe_id), units_atomic_to_natural(energy=energy), &
    & units_atomic_to_natural(force=forces(:,pes%probe_id)), alpha, iter
    if (istat .ne. 0) call io_err("do_pes_scan: Could not write to output file")

    ! new step size:
    delta_pos = pos - pos_old
    delta_grad = grad - grad_old
    alpha = (delta_pos*delta_grad)/(delta_grad*delta_grad)

    ! old values for next iter:
    pos_old = pos
    grad_old = grad
    energy_old = energy
  end do
end subroutine pes_tpsd_line_minimise

end module pes_scan
