module phonon

  use constants,  only: dp
  use cell,       only: cell_type
  implicit none

  private

  ! default phonon parameters:
  real(kind=dp),                parameter ::  default_dx = 0.005_dp                 ! Ang
  integer,        dimension(3), parameter ::  default_supercell = (/ 1, 1, 1 /)
  real(kind=dp),                parameter ::  default_qpoint_path_spacing = 0.1_dp  ! 1/Ang
  integer,        dimension(3), parameter ::  default_mp_grid_size = (/ 1, 1, 1 /)
  logical,                      parameter ::  default_phonon_geom_opt = .false.


  type, public  ::  phonon_type
    ! public interface variables 
    integer,      dimension(3)  ::  supercell_size
    real(kind=dp)               ::  dx                  ! perturbation size
    real(kind=dp)               ::  qpoint_path_spacing
    integer,      dimension(3)  ::  mp_grid_size
    logical                     ::  do_geom_opt

    ! internal variables
    logical               ::  grid_specified
    logical               ::  path_specified
    type(cell_type)       ::  calc_cell                 ! calculation supercell
    integer               ::  ncells

    ! mapping variables (from local cells to supercell)
    integer,          allocatable,  dimension(:,:)        ::  global_atom_index
    integer,          allocatable,  dimension(:)          ::  local_atom_index
    integer,          allocatable,  dimension(:)          ::  atom_cell_index

    real(kind=dp),    allocatable,  dimension(:,:)        ::  supercell_corner_vectors

    ! tmp work arrays
    real(kind=dp),    allocatable,  dimension(:,:)        ::  packed_force
    real(kind=dp),    allocatable,  dimension(:,:)        ::  packed_dforce

    real(kind=dp),    allocatable,  dimension(:,:,:,:,:)  ::  force_const_matrix
    complex(kind=dp), allocatable,  dimension(:,:,:,:)    ::  dynamical_matrix
    complex(kind=dp), allocatable,  dimension(:,:)        ::  packed_dynamical_matrix
    
    real(kind=dp),    allocatable,  dimension(:)          ::  frequencies

    ! qpoint variables
    real(kind=dp),    allocatable,  dimension(:,:)        ::  path_qpoint_nodes
    real(kind=dp),    allocatable,  dimension(:,:)        ::  qpoints
  end type

  public  ::  phonon_init
  public  ::  do_phonon_calculation
  public  ::  phonon_deallocate

contains

subroutine phonon_read_input(phonon)
  use io, only: io_input_get_single_value, max_line_len, io_input_get_data, io_str_get_num_tokens, &
              & io_str_get_token, io_str_to_int, io_err
  implicit none
  type(phonon_type),  intent(inout) ::  phonon

  character(len=max_line_len) ::  keyword
  character(len=max_line_len),  allocatable,  dimension(:)  ::  input_data
  logical ::  kw_found
  integer ::  iline, nlines
  integer ::  i, ntokens, istat

  keyword = 'phonon_dx'
  call io_input_get_single_value(keyword, phonon%dx, default_dx)


  keyword = 'phonon_supercell'
  call io_input_get_data(keyword, input_data, kw_found)
  if (kw_found) then
    ! check for single line
    if (size(input_data,1) .ne. 1) call io_err("phonon_read_input: "//trim(keyword)//" must be a single line")

    ! check for number of values - must be 1 or 3
    ntokens = io_str_get_num_tokens(input_data(1))
    if (ntokens .eq. 1) then
      ! set all values to equal single value.. (isotropic expansion)
      phonon%supercell_size(1) = io_str_to_int(io_str_get_token(input_data(1), 1))
      phonon%supercell_size(:) = phonon%supercell_size(1)
    else if (ntokens .eq. 3) then
      ! set all 3 independently
      read(input_data(1), *, iostat=istat) phonon%supercell_size
      if (istat .ne. 0) call io_err("phonon_read_input: Could not read "//trim(keyword)//" values")
    else
      call io_err("phonon_read_input: "//trim(keyword)//" must contain 1 or 3 values")
    end if
    do i = 1, 3
      if (phonon%supercell_size(i) .lt. 1) then
        call io_err("phonon_read_input: "//trim(keyword)//": Each component must be >= 1")
      end if
    end do
  else
    ! use default
    phonon%supercell_size(:) = default_supercell(:)
  end if


  keyword = 'phonon_qpoint_mp_grid'
  call io_input_get_data(keyword, input_data, kw_found)
  if (kw_found) then
    ! check for single line
    if (size(input_data,1) .ne. 1) call io_err("phonon_read_input: "//trim(keyword)//" must be a single line")

    ! check for number of values - must be 1 or 3
    ntokens = io_str_get_num_tokens(input_data(1))
    if (ntokens .eq. 1) then
      ! set all values to equal single value.. (isotropic expansion)
      phonon%mp_grid_size(1) = io_str_to_int(io_str_get_token(input_data(1), 1))
      phonon%mp_grid_size(:) = phonon%mp_grid_size(1)
    else if (ntokens .eq. 3) then
      ! set all 3 independently
      read(input_data(1), *, iostat=istat) phonon%mp_grid_size(:)
      if (istat .ne. 0) call io_err("phonon_read_input: Could not read "//trim(keyword)//" values")
    else
      call io_err("phonon_read_input: "//trim(keyword)//" must contain 1 or 3 values")
    end if

    do i = 1, 3
      if (phonon%mp_grid_size(i) .lt. 1) then
        call io_err("phonon_read_input: "//trim(keyword)//": Each component must be >= 1")
      end if
    end do
    ! remember that grid was specified
    phonon%grid_specified = .true.
  else
    ! use default
    phonon%mp_grid_size(:) = default_mp_grid_size(:)
  end if


  keyword = 'phonon_qpoint_path_spacing'
  call io_input_get_single_value(keyword, phonon%qpoint_path_spacing, default_qpoint_path_spacing)
  if (phonon%qpoint_path_spacing .le. 0.0_dp) call io_err("phonon_read_input: "//trim(keyword)//": Must be > 0.0")


  keyword = 'phonon_qpoint_path'
  call io_input_get_data(keyword, input_data, kw_found)
  ! have no default path.. use MP_grid instead..
  if (kw_found) then
    ! check we have at least 2 points (lines)
    nlines = size(input_data,1)
    if (nlines .lt. 2) call io_err("phonon_read_input: "//trim(keyword)//" must have at least 2 points")

    ! allocate path_qpoint_nodes array
    allocate(phonon%path_qpoint_nodes(3,nlines), stat=istat)
    if (istat .ne. 0) call io_err("phonon_read_input: Could not allocate phonon%path_qpoint_nodes array")

    do iline = 1, nlines
      ! check for number of values - must have 3 points..
      ntokens = io_str_get_num_tokens(input_data(iline))
      if (ntokens .ne. 3) call io_err("phonon_read_input: "//trim(keyword)//" must have 3 values on each line")

      read(input_data(iline), *, iostat=istat) phonon%path_qpoint_nodes(:,iline)
      if (istat .ne. 0) call io_err("phonon_read_input: "//trim(keyword)//": Could not read node data")
    end do
    ! remember that path was specified
    phonon%path_specified = .true.
  end if

  keyword = 'phonon_geom_opt'
  call io_input_get_single_value(keyword, phonon%do_geom_opt, default_phonon_geom_opt)

  if (allocated(input_data)) then
    deallocate(input_data, stat=istat)
    if (istat .ne. 0) call io_err("phonon_read_input: Could not deallocate input_data array")
  end if
end subroutine phonon_read_input

subroutine phonon_init(phonon, cell)
  use io,         only: io_err
  use constants,  only: units_natural_to_atomic
  implicit none
  type(phonon_type),  intent(inout) ::  phonon
  type(cell_type),    intent(in)    ::  cell
  real(kind=dp) ::  one_ang_in_bohr

  phonon%grid_specified = .false.
  phonon%path_specified = .false.

  call phonon_read_input(phonon)

  if (phonon%grid_specified .and. phonon%path_specified) then
    call io_err("phonon_init: Only *one* of a qpoint path, or qpoint grid must be specified at once")
  else if ((.not. phonon%grid_specified) .and. (.not. phonon%path_specified)) then
    ! if neither, then use default MP grid.. (value is set in phonon_read_input)
    ! could determine default MP grid size based on a default spacing (as with points along paths)
    phonon%grid_specified = .true.
  end if

  ! convert units
  phonon%dx = units_natural_to_atomic(length=phonon%dx)

  ! convert to 1/Ang to 1/Bohr..
  one_ang_in_bohr = units_natural_to_atomic(length=1.0_dp)
  phonon%qpoint_path_spacing = phonon%qpoint_path_spacing/one_ang_in_bohr
  
  ! number of unit cells..
  phonon%ncells = phonon%supercell_size(1)*phonon%supercell_size(2)*phonon%supercell_size(3)

  call phonon_allocate(phonon, cell)
end subroutine phonon_init

subroutine phonon_allocate(phonon, cell)
  use io, only: io_err
  implicit none
  type(phonon_type),  intent(inout) ::  phonon
  type(cell_type),    intent(in)    ::  cell
  integer ::  istat, supercell_natoms

  supercell_natoms = cell%natoms*phonon%ncells

  allocate(phonon%packed_force(3, supercell_natoms), stat=istat)
  if (istat .ne. 0) call io_err("phonon_allocate: Could not allocate phonon%packed_force array")

  allocate(phonon%packed_dforce(3, supercell_natoms), stat=istat)
  if (istat .ne. 0) call io_err("phonon_allocate: Could not allocate phonon%packed_dforce array")

  allocate(phonon%force_const_matrix(3, 3, cell%natoms, cell%natoms, phonon%ncells), stat=istat)
  if (istat .ne. 0) call io_err("phonon_allocate: Could not allocate phonon%force_const_matrix array")

  allocate(phonon%global_atom_index(cell%natoms, phonon%ncells), stat=istat)
  if (istat .ne. 0) call io_err("phonon_allocate: Could not allocate phonon%global_atom_index array")

  allocate(phonon%local_atom_index(supercell_natoms), stat=istat)
  if (istat .ne. 0) call io_err("phonon_allocate: Could not allocate phonon%local_atom_index array")

  allocate(phonon%atom_cell_index(supercell_natoms), stat=istat)
  if (istat .ne. 0) call io_err("phonon_allocate: Could not allocate phonon%atom_cell_index array")

  allocate(phonon%dynamical_matrix(3, 3, cell%natoms, cell%natoms), stat=istat)
  if (istat .ne. 0) call io_err("phonon_allocate: Could not allocate phonon%dynamical_matrix array")

  allocate(phonon%packed_dynamical_matrix(3*cell%natoms, 3*cell%natoms), stat=istat)
  if (istat .ne. 0) call io_err("phonon_allocate: Could not allocate phonon%packed_dynamical_matrix array")

  allocate(phonon%frequencies(3*cell%natoms), stat=istat)
  if (istat .ne. 0) call io_err("phonon_allocate: Could not allocate phonon%frequencies array")

  allocate(phonon%supercell_corner_vectors(3, phonon%ncells), stat=istat)
  if (istat .ne. 0) call io_err("phonon_allocate: Could not allocate phonon%supercell_corner_vectors array")
end subroutine phonon_allocate

subroutine phonon_deallocate(phonon)
  use cell, only: cell_deallocate
  use io,   only: io_err
  implicit none
  type(phonon_type),  intent(inout) ::  phonon
  integer ::  istat

  call cell_deallocate(phonon%calc_cell)

  deallocate(phonon%packed_force, stat=istat)
  if (istat .ne. 0) call io_err("phonon_deallocate: Could not deallocate phonon%packed_force array")

  deallocate(phonon%packed_dforce, stat=istat)
  if (istat .ne. 0) call io_err("phonon_deallocate: Could not deallocate phonon%packed_dforce array")

  deallocate(phonon%force_const_matrix, stat=istat)
  if (istat .ne. 0) call io_err("phonon_deallocate: Could not deallocate phonon%force_const_matrix array")

  deallocate(phonon%global_atom_index, stat=istat)
  if (istat .ne. 0) call io_err("phonon_deallocate: Could not deallocate phonon%global_atom_index array")

  deallocate(phonon%local_atom_index, stat=istat)
  if (istat .ne. 0) call io_err("phonon_deallocate: Could not deallocate phonon%local_atom_index array")

  deallocate(phonon%atom_cell_index, stat=istat)
  if (istat .ne. 0) call io_err("phonon_deallocate: Could not deallocate phonon%atom_cell_index array")

  deallocate(phonon%dynamical_matrix, stat=istat)
  if (istat .ne. 0) call io_err("phonon_deallocate: Could not deallocate phonon%dynamical_matrix array")

  deallocate(phonon%packed_dynamical_matrix, stat=istat)
  if (istat .ne. 0) call io_err("phonon_deallocate: Could not deallocate phonon%packed_dynamical_matrix array")

  deallocate(phonon%frequencies, stat=istat)
  if (istat .ne. 0) call io_err("phonon_deallocate: Could not deallocate phonon%frequencies array")

  deallocate(phonon%supercell_corner_vectors, stat=istat)
  if (istat .ne. 0) call io_err("phonon_deallocate: Could not deallocate phonon%supercell_corner_vectors array")
end subroutine phonon_deallocate

subroutine phonon_build_supercell(phonon, cell)
  use cell, only: cell_allocate, cell_calc_vol, cell_get_recip, cell_frac_to_cart, cell_min_img_init, cell_min_img
  implicit none
  type(phonon_type),  intent(inout) ::  phonon
  type(cell_type),    intent(in)    ::  cell
  ! local vars
  integer,        dimension(3)  ::  isc
  real(kind=dp),  dimension(3)  ::  dx_frac
  integer ::  icell, jcell, kcell
  integer ::  iatom, icomp, old_atom, new_atom, cell_index

  dx_frac(:) = (/ 1.0_dp/real(phonon%supercell_size(1),dp), &
                & 1.0_dp/real(phonon%supercell_size(2),dp), &
                & 1.0_dp/real(phonon%supercell_size(3),dp) /)

  phonon%calc_cell%natoms = cell%natoms*phonon%ncells

  call cell_allocate(phonon%calc_cell)

  ! expand lattice vectors - icomp runs over a, b, c
  do icomp = 1, 3
    phonon%calc_cell%lattice_vectors(:,icomp) = cell%lattice_vectors(:,icomp)*real(phonon%supercell_size(icomp),dp)
  end do

  ! shift everything down to 'corner' of cell
  do iatom = 1, cell%natoms
    do icomp = 1, 3
      phonon%calc_cell%atom_frac_pos(icomp, iatom) = cell%atom_frac_pos(icomp, iatom)*dx_frac(icomp)
    end do
  end do

  ! add in extra atoms to phonon%calc_cell
  new_atom = 1
  cell_index = 0
  do icell = 1, phonon%supercell_size(1)
    do jcell = 1, phonon%supercell_size(2)
      do kcell = 1, phonon%supercell_size(3)
        isc(:) = (/ icell-1, jcell-1, kcell-1 /)
        cell_index = cell_index + 1
        do old_atom = 1, cell%natoms
          do icomp = 1, 3
            phonon%calc_cell%atom_frac_pos(icomp, new_atom) = phonon%calc_cell%atom_frac_pos(icomp, old_atom) &
                                                              & + real(isc(icomp),dp)*dx_frac(icomp)
          end do
          phonon%calc_cell%atom_mass(new_atom) = cell%atom_mass(old_atom)
          phonon%calc_cell%atom_species(new_atom) = cell%atom_species(old_atom)

          ! get cell index, and store local and global atom numbers
          phonon%global_atom_index(old_atom, cell_index) = new_atom
          phonon%local_atom_index(new_atom) = old_atom
          phonon%atom_cell_index(new_atom) = cell_index
          new_atom = new_atom + 1
        end do
      end do
    end do
  end do

  ! update rest of cell - order important
  call cell_calc_vol(phonon%calc_cell)
  call cell_get_recip(phonon%calc_cell)
  call cell_frac_to_cart(phonon%calc_cell)
  call cell_min_img_init(phonon%calc_cell)

  ! vectors from unit cell to all cells in supercell
  old_atom = 1 ! first atom in each cell
  do icell = 1, phonon%ncells
    new_atom = phonon%global_atom_index(old_atom, icell)
    call cell_min_img(phonon%calc_cell, old_atom, new_atom, phonon%supercell_corner_vectors(:,icell))
  end do
end subroutine phonon_build_supercell

! use global atom index..
subroutine phonon_perturb_atom(phonon, iatom, direction)
  use potential,  only: pot_get_forces
  use cell,       only: cell_cart_to_frac, cell_frac_to_cart, cell_shift_to_unit_cell
  use io,         only: io_err
  implicit none
  type(phonon_type),  intent(inout) ::  phonon
  integer,            intent(in)    ::  iatom           ! atom we are moving...
  integer,            intent(in)    ::  direction       ! 1,2,3 is x,y,z

  ! local vars:
  integer ::  perturb_dir
  real(kind=dp) ::  energy
  real(kind=dp),  dimension(3)  ::  orig_pos

  ! zero dforce array..
  phonon%packed_dforce(:,:) = 0.0_dp

  ! store original cartesian position
  orig_pos(:) = phonon%calc_cell%atom_cart_pos(:,iatom)

  ! centered difference
  ! dy/dx = (y(x+h) - y(x-h)) / 2h

  do perturb_dir = -1, 1, 2

    ! new position
    phonon%calc_cell%atom_cart_pos(:, iatom) = orig_pos(:)
    phonon%calc_cell%atom_cart_pos(direction, iatom) = orig_pos(direction) + real(perturb_dir,dp)*phonon%dx

    ! sync cell coords
    call cell_cart_to_frac(phonon%calc_cell)
    call cell_shift_to_unit_cell(phonon%calc_cell)
    call cell_frac_to_cart(phonon%calc_cell)

    ! calc forces
    call pot_get_forces(phonon%calc_cell, phonon%packed_force, energy)

    ! dforce = f(y+h) - f(y-h)
    phonon%packed_dforce(:,:) = phonon%packed_dforce(:,:) + real(perturb_dir,dp)*phonon%packed_force(:,:)
  end do

  ! make sure cell is as we found it
  phonon%calc_cell%atom_cart_pos(:, iatom) = orig_pos(:)
  call cell_cart_to_frac(phonon%calc_cell)
  call cell_shift_to_unit_cell(phonon%calc_cell)
  call cell_frac_to_cart(phonon%calc_cell)
       
  ! dforce = (f(y+u) - f(y-u)) / 2u
  phonon%packed_dforce(:,:) = 0.5_dp*(phonon%packed_dforce(:,:))/phonon%dx
end subroutine phonon_perturb_atom

subroutine phonon_calc_force_const_matrix(phonon, cell)
  implicit none
  type(phonon_type),  intent(inout) ::  phonon
  type(cell_type),    intent(in)    ::  cell
  ! local vars:
  integer ::  local_iatom, global_iatom
  integer ::  local_jatom, global_jatom
  integer ::  local_cell_index
  integer ::  icomp

  ! iatom and comp are the perturbation atoms/directions
  ! jatom etc are the atoms which feel the response to the perturbation

  ! local_iatom is index in unit cell..
  ! this may not be the same number as the global atom index
  do local_iatom = 1, cell%natoms

    ! get global atom index - cell_index is always 1
    global_iatom = phonon%global_atom_index(local_iatom, 1)

    ! perturbation direction
    do icomp = 1, 3

      ! perturb global_iatom along icomp
      call phonon_perturb_atom(phonon, global_iatom, icomp)

      ! need to map packed dforce to unpacked FCM array
      do global_jatom = 1, phonon%calc_cell%natoms

        ! jatom is global atom index, need to get local atom index..
        local_jatom = phonon%local_atom_index(global_jatom)

        ! and local cell index...
        local_cell_index = phonon%atom_cell_index(global_jatom)

        ! do mapping...
        ! FCM(comp, comp', atom, atom', cell_num).. where atom' and comp' are the perturbations
        phonon%force_const_matrix(:, icomp, local_jatom, local_iatom, local_cell_index) = -phonon%packed_dforce(:, global_jatom)
      end do

    end do
  end do
end subroutine phonon_calc_force_const_matrix

subroutine phonon_calc_dynamical_matrix(phonon, cell, qpoint)
  use constants,  only: cmplx_i
  use io,         only: io_err
  implicit none
  type(phonon_type),                intent(inout) ::  phonon
  type(cell_type),                  intent(in)    ::  cell
  real(kind=dp),      dimension(3), intent(in)    ::  qpoint
  ! local variables
  real(kind=dp),  dimension(3)  ::  qvector
  complex(kind=dp)  ::  phase, contribution
  real(kind=dp)     ::  q_dot_R_cell
  integer           ::  icomp, jcomp, iatom, jatom, icell
  real(kind=dp)     ::  mass_iatom, mass_jatom

  ! qpoint is in fractions of recip lattice vectors, qvector is an absolute position
  qvector = matmul(transpose(phonon%calc_cell%recip_lattice), qpoint)

  phonon%dynamical_matrix(:,:,:,:) = (0.0_dp, 0.0_dp)
  do iatom = 1, cell%natoms

    ! just in case ordering is not preserved in each cell..
    mass_iatom = phonon%calc_cell%atom_mass(phonon%global_atom_index(iatom,1))

    do jatom = 1, cell%natoms

      ! use same cell as above, just to be safe
      mass_jatom = phonon%calc_cell%atom_mass(phonon%global_atom_index(jatom,1))

      do icell = 1, phonon%ncells
        q_dot_R_cell = dot_product(qvector(:), phonon%supercell_corner_vectors(:, icell))
        phase = exp(-cmplx_i*cmplx(q_dot_R_cell,kind=dp))

        do icomp = 1,3
          do jcomp = 1,3
            ! calculate contribution for this cell; jatom, iatom, jcomp, icomp
            contribution = cmplx(phonon%force_const_matrix(jcomp, icomp, jatom, iatom, icell),kind=dp)*phase
            phonon%dynamical_matrix(jcomp, icomp, jatom, iatom) = phonon%dynamical_matrix(jcomp, icomp, jatom, iatom) + contribution
          end do  !jcomp
        end do  !icomp
      end do  !icell

      ! divide through by sqrt(mass_i * mass_j)
      phonon%dynamical_matrix(:, :, jatom, iatom) = &
      &                 phonon%dynamical_matrix(:, :, jatom, iatom)/cmplx(sqrt(mass_iatom*mass_jatom),kind=dp)
    end do  !jatom
  end do  !iatom
end subroutine phonon_calc_dynamical_matrix

subroutine phonon_diagonalize_dynamical_matrix(phonon, cell)
  use io, only: io_err
  implicit none
  type(phonon_type),  intent(inout) ::  phonon
  type(cell_type),    intent(in)    ::  cell

  ! local vars:
  complex(kind=dp), dimension(:), allocatable ::  work
  real(kind=dp),    dimension(:), allocatable ::  rwork
  integer ::  k1, k2, a1, a2, m, n
  integer ::  lwork, istat

  ! put dynamical matrix in packed form
  do k1 = 1, cell%natoms
    do a1 = -2, 0
      n = 3*k1 + a1
      do k2 = 1, cell%natoms
        do a2 = -2, 0
          m = 3*k2 + a2
          phonon%packed_dynamical_matrix(m,n) = phonon%dynamical_matrix(a2, a1, k2, k1)
        end do
      end do
    end do
  end do

  n = 3*cell%natoms
  lwork = 2*n - 1

  allocate(work(lwork), stat=istat)
  if (istat .ne. 0) call io_err("phonon_diagonalize_dynamical_matrix: Could not allocate work array")

  allocate(rwork(3*n-2), stat=istat)
  if (istat .ne. 0) call io_err("phonon_diagonalize_dynamical_matrix: Could not allocate rwork array")

  call zheev('N', 'U', n, phonon%packed_dynamical_matrix, n, phonon%frequencies, work, lwork, rwork, istat)
  if (istat .ne. 0) call io_err("phonon_diagonalize_dynamical_matrix: Could not diagonalize dynamical matrix")

  deallocate(work, stat=istat)
  if (istat .ne. 0) call io_err("phonon_diagonalize_dynamical_matrix: Could not deallocate work array")

  deallocate(rwork, stat=istat)
  if (istat .ne. 0) call io_err("phonon_diagonalize_dynamical_matrix: Could not deallocate rwork array")
end subroutine phonon_diagonalize_dynamical_matrix

subroutine phonon_define_qpoints_path(phonon)
  use io, only: io_err
  implicit none
  type(phonon_type),  intent(inout) ::  phonon
  ! local vars
  integer ::  inode, ipt, jpt, jstart, npts, total_npts, istat
  real(kind=dp),  dimension(3)  ::  start_qpoint, end_qpoint
  real(kind=dp),  dimension(3)  ::  start_abs_qpoint, end_abs_qpoint, path_vec, path_dx
  real(kind=dp) ::  length_path

  ! first get total number of points.. (slightly slow way of doing this..)
  total_npts = 1      ! include start point..
  do inode = 1, size(phonon%path_qpoint_nodes, 2)-1
    start_qpoint(:) = phonon%path_qpoint_nodes(:,inode)
    end_qpoint(:) = phonon%path_qpoint_nodes(:,inode+1)

    start_abs_qpoint = matmul(transpose(phonon%calc_cell%recip_lattice), start_qpoint)
    end_abs_qpoint = matmul(transpose(phonon%calc_cell%recip_lattice), end_qpoint)
    path_vec = end_abs_qpoint - start_abs_qpoint

    length_path = sqrt(dot_product(path_vec, path_vec))

    ! round up number of points (so phonon_qpoint_path_spacing is always >= actual path spacing)
    npts = ceiling(length_path/phonon%qpoint_path_spacing)

    total_npts = total_npts + npts
  end do

  ! now we can allocate to correct number of points 
  allocate(phonon%qpoints(3,total_npts), stat=istat)
  if (istat .ne. 0) call io_err("phonon_define_qpoint_path: Could not allocate phonon%qpoints array")
 
  ! now need to set up the qpoints
  ipt = 0
  jstart = 0
  do inode = 1, size(phonon%path_qpoint_nodes, 2)-1
    start_qpoint(:) = phonon%path_qpoint_nodes(:,inode)
    end_qpoint(:) = phonon%path_qpoint_nodes(:,inode+1)

    start_abs_qpoint = matmul(transpose(phonon%calc_cell%recip_lattice), start_qpoint)
    end_abs_qpoint = matmul(transpose(phonon%calc_cell%recip_lattice), end_qpoint)
    path_vec = end_abs_qpoint - start_abs_qpoint

    length_path = sqrt(dot_product(path_vec, path_vec))

    ! round up number of points (so phonon_qpoint_path_spacing is always >= actual path spacing)
    npts = ceiling(length_path/phonon%qpoint_path_spacing)

    path_dx = (end_qpoint - start_qpoint)/real(npts,dp)

    ! skip first point if past first part of path
    if (inode .gt. 1) jstart = 1

    do jpt = jstart, npts
      ipt = ipt + 1
      phonon%qpoints(:,ipt) = start_qpoint(:) + real(jpt,dp)*path_dx(:)
    end do
  end do
end subroutine phonon_define_qpoints_path

subroutine phonon_define_mp_grid(phonon)
  use io, only: io_err
  implicit none
  type(phonon_type),  intent(inout) ::  phonon
  ! local variables
  real(kind=dp) ::  up, ur, us
  integer       ::  p, r, s, iq
  integer       ::  np, nr, ns, total_nq
  integer       ::  istat

  np = phonon%mp_grid_size(1)
  nr = phonon%mp_grid_size(2)
  ns = phonon%mp_grid_size(3)

  total_nq = np*nr*ns
  allocate(phonon%qpoints(3, total_nq), stat=istat)
  if (istat .ne. 0) call io_err("phonon_define_mp_grid: Could not allocate phonon%qpoints array")

  ! first define the evenly spaced grid points inside the BZ - using Monkhorst-Pack scheme
  iq = 0
  do p = 1, np
    up = real((2*p - np - 1),dp)/real((2*np),dp)

    do r = 1, nr
      ur = real((2*r - nr - 1),dp)/real((2*nr),dp)

      do s = 1, ns
        us = real((2*s - ns - 1),dp)/real((2*ns),dp)
        iq = iq + 1
        phonon%qpoints(:, iq) = (/ up, ur, us /)
      end do !s

    end do  !r
  end do  !p
end subroutine phonon_define_mp_grid

subroutine do_phonon_calculation(phonon, cell)
  use io, only: stdout, io_err
  implicit none
  type(phonon_type),  intent(inout) ::  phonon
  type(cell_type),    intent(in)    ::  cell
  ! local vars
  integer :: ifreq, ipoint, istat

  call phonon_build_supercell(phonon, cell)

  call phonon_calc_force_const_matrix(phonon, cell)

!  do a = 1, 3
!  do ap = 1, 3
!  do k = 1, cell%natoms
!  xsum = 0.0_dp
!  do i = 1, phonon%ncells
!    do kp = 1, cell%natoms
!      xsum = xsum + phonon%force_const_matrix(a, ap, k, kp, i)
!    end do
!  end do
!  print *, xsum
!  end do
!  end do
!  end do

  ! only one of these will be set at once (phonon_init ensures this)
  if (phonon%grid_specified) then
    call phonon_define_mp_grid(phonon)
  else if (phonon%path_specified) then
    call phonon_define_qpoints_path(phonon)
  end if

  ! for each qpoint, calc dynamical matrix and diagonalize
  do ipoint = 1, size(phonon%qpoints,2)
    call phonon_calc_dynamical_matrix(phonon, cell, phonon%qpoints(:,ipoint))
    call phonon_diagonalize_dynamical_matrix(phonon, cell)
!    do i = 1, size(phonon%packed_dynamical_matrix,1) 
!      do j = 1, size(phonon%packed_dynamical_matrix,2)
!        print *, phonon%packed_dynamical_matrix(i,j)-conjg(phonon%packed_dynamical_matrix(j,i))
      !if ((phonon%packed_dynamical_matrix(i,j)-conjg(phonon%packed_dynamical_matrix(j,i))) .ne. (0.0_dp, 0.0_dp)) then
      !  print *, i, j, phonon%packed_dynamical_matrix(i,j), phonon%packed_dynamical_matrix(j,i)
      !end if
!      end do
!    end do

    write(stdout, *) 'qpoint:'
    write(stdout, *) phonon%qpoints(:,ipoint)
    write(stdout, *)
    write(stdout, *) 'frequencies**2'
    do ifreq = 1, size(phonon%frequencies,1)
      write(stdout, *) phonon%frequencies(ifreq) 
    end do
    write(stdout, *)
  end do

  deallocate(phonon%qpoints, stat=istat)
  if (istat .ne. 0) call io_err("do_phonon_calculation: Could not deallocate phonon%qpoints array")
end subroutine do_phonon_calculation
end module phonon
