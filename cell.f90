!****m* cell/cell =============================================================!
! NAME                                                                         !
!   cell                                                                       !
!------------------------------------------------------------------------------!
! DESCRIPTION                                                                  !
!   Module containing the cell type and routines for manipulating cells.       !
!                                                                              !
!   The cell contains the atomic species, masses and positions (in cartesian   !
!   and fractional spaces), as well as the lattice and reciprocal lattice      !
!   vectors.                                                                   !
!------------------------------------------------------------------------------!
! AUTHOR                                                                       !
!   Aaron Hopkinson                                                            !
!****==========================================================================!
module cell
  use constants, only: dp
  implicit none

  private

  type, public  ::  cell_type
    integer                                     ::  natoms
    real(kind=dp)                               ::  volume
    real(kind=dp),  dimension(3,3)              ::  lattice_vectors                   ! columns are a, b, c
    real(kind=dp),  dimension(3,3)              ::  recip_lattice
    real(kind=dp),  dimension(:,:), allocatable ::  atom_frac_pos                     ! ndims, natoms
    real(kind=dp),  dimension(:,:), allocatable ::  atom_cart_pos                     ! ndims, natoms
    real(kind=dp),  dimension(:),   allocatable ::  atom_mass
    integer,        dimension(:),   allocatable ::  atom_species

    real(kind=dp)                               ::  min_img_cutoff
  end type cell_type

  ! cell independent min image vars:
  logical,  save                                ::  min_img_initialized = .false.
  integer,  save, dimension(3,27)               ::  min_img_c_vectors

  ! temporary cell perturbation vars (not required after init):
  real(kind=dp) ::  lattice_perturb, atom_perturb

  ! public routines:
  public  ::  cell_init
  public  ::  cell_write
  public  ::  cell_allocate
  public  ::  cell_deallocate
  public  ::  cell_copy
  public  ::  cell_read_input
  public  ::  cell_shift_to_unit_cell
  public  ::  cell_frac_to_cart
  public  ::  cell_cart_to_frac
  public  ::  cell_get_recip
  public  ::  cell_min_img_init
  public  ::  cell_min_img
  public  ::  cell_min_img_vec
  public  ::  cell_calc_vol
  public  ::  cell_supercell
  public  ::  cell_perturb
  public  ::  cell_reallocate
  public  ::  cell_read_checkpoint
  public  ::  cell_write_checkpoint

contains

! checkpointing things.. checkpoint_open must be called first! these are low level...
subroutine cell_read_checkpoint(cell)
  use checkpoint, only: checkpoint_unit
  implicit none
  type(cell_type),  intent(inout) ::  cell
  integer :: i, j

  read(checkpoint_unit)  cell%natoms
  call cell_reallocate(cell)
  read(checkpoint_unit)  cell%volume
  do i = 1, 3
    do j = 1, 3
      read(checkpoint_unit)  cell%lattice_vectors(j,i)
    end do
  end do
  do i = 1, 3
    do j = 1, 3
      read(checkpoint_unit)  cell%recip_lattice(j,i)
    end do
  end do
  do i = 1, cell%natoms
    do j = 1, 3
      read(checkpoint_unit)  cell%atom_frac_pos(j,i)
    end do
  end do
  do i = 1, cell%natoms
    do j = 1, 3
      read(checkpoint_unit)  cell%atom_cart_pos(j,i)
    end do
  end do
  do i = 1, cell%natoms
    read(checkpoint_unit)  cell%atom_mass(i)
  end do
  do i = 1, cell%natoms
    read(checkpoint_unit)  cell%atom_species(i)
  end do
  read(checkpoint_unit)  cell%min_img_cutoff
end subroutine cell_read_checkpoint

subroutine cell_write_checkpoint(cell)
  use checkpoint, only: checkpoint_unit
  implicit none
  type(cell_type),  intent(in) ::  cell
  integer :: i, j

  write(checkpoint_unit)  cell%natoms
  write(checkpoint_unit)  cell%volume
  do i = 1, 3
    do j = 1, 3
      write(checkpoint_unit)  cell%lattice_vectors(j,i)
    end do
  end do
  do i = 1, 3
    do j = 1, 3
      write(checkpoint_unit)  cell%recip_lattice(j,i)
    end do
  end do
  do i = 1, cell%natoms
    do j = 1, 3
      write(checkpoint_unit)  cell%atom_frac_pos(j,i)
    end do
  end do
  do i = 1, cell%natoms
    do j = 1, 3
      write(checkpoint_unit)  cell%atom_cart_pos(j,i)
    end do
  end do
  do i = 1, cell%natoms
    write(checkpoint_unit)  cell%atom_mass(i)
  end do
  do i = 1, cell%natoms
    write(checkpoint_unit)  cell%atom_species(i)
  end do
  write(checkpoint_unit)  cell%min_img_cutoff
end subroutine cell_write_checkpoint


!****s* cell/perturb ==========================================================!
! NAME                                                                         !
!   cell_perturb (PUBLIC)                                                      !
!------------------------------------------------------------------------------!
! ARGUMENTS                                                                    !
!   type(cell_type),    intent(inout)   ::  cell                               !
!   real(kind=dp),      intent(in)      ::  atomic_pert                        !
!   real(kind=dp),      intent(in)      ::  lattice_pert                       !
!------------------------------------------------------------------------------!
! DESCRPTION                                                                   !
!   A subroutine to perturb the atomic positions and/or lattice parameters of  !
!   a given cell.                                                              !
!------------------------------------------------------------------------------!
! AUTHOR                                                                       !
!   Edward Higgins                                                             !
!==============================================================================!
subroutine cell_perturb(cell, atomic_pert, lattice_pert)
  use constants,  only: pi
  use io,         only: io_err
  use algor,      only: algor_uniform_rand
  implicit none
  type(cell_type),            intent(inout) ::  cell
  real(kind=dp),    optional, intent(in)    ::  atomic_pert
  real(kind=dp),    optional, intent(in)    ::  lattice_pert

  real(kind=dp),    dimension(3,3)  ::  new_lattice
  real(kind=dp)                     ::  r, theta, phi
  integer                           ::  iatom

  if ((.not. present(atomic_pert)) .and. (.not. present(lattice_pert))) &
    call io_err("cell_perturb: No perturbation to apply")

  if (present(lattice_pert)) then
    if (lattice_pert .gt. 0.0_dp) then
      new_lattice(1,:) =  [ cell%lattice_vectors(1,1) + algor_uniform_rand(-lattice_pert/2, lattice_pert/2), &
                            cell%lattice_vectors(1,2) + algor_uniform_rand(-lattice_pert/2, lattice_pert/2), &
                            cell%lattice_vectors(1,3) + algor_uniform_rand(-lattice_pert/2, lattice_pert/2)  ]

      new_lattice(2,:) =  [ cell%lattice_vectors(2,1) + algor_uniform_rand(-lattice_pert/2, lattice_pert/2), &
                            cell%lattice_vectors(2,2) + algor_uniform_rand(-lattice_pert/2, lattice_pert/2), &
                            cell%lattice_vectors(2,3) + algor_uniform_rand(-lattice_pert/2, lattice_pert/2)  ]

      new_lattice(3,:) =  [ cell%lattice_vectors(3,1) + algor_uniform_rand(-lattice_pert/2, lattice_pert/2), &
                            cell%lattice_vectors(3,2) + algor_uniform_rand(-lattice_pert/2, lattice_pert/2), &
                            cell%lattice_vectors(3,3) + algor_uniform_rand(-lattice_pert/2, lattice_pert/2)  ]

      cell%lattice_vectors = new_lattice

      ! Sync the cell
      call cell_calc_vol(cell)
      call cell_get_recip(cell)
      call cell_frac_to_cart(cell)
      call cell_min_img_init(cell)
    end if
  end if

  if (present(atomic_pert)) then
    if (atomic_pert .gt. 0.0_dp) then
      do iatom = 1, cell%natoms
        r = algor_uniform_rand(0.0_dp, 1.0_dp) ** (1.0_dp/3.0_dp)
        r = r * atomic_pert
        theta = algor_uniform_rand(-pi, pi)
        phi = algor_uniform_rand(0.0_dp, 1.0_dp)
        phi = acos(2.0_dp*phi - 1.0_dp)
        cell%atom_cart_pos(1,iatom) = cell%atom_cart_pos(1,iatom) + r * cos(theta) * sin(phi)
        cell%atom_cart_pos(2,iatom) = cell%atom_cart_pos(2,iatom) + r * sin(theta) * sin(phi)
        cell%atom_cart_pos(3,iatom) = cell%atom_cart_pos(3,iatom) + r * cos(phi)
      end do

      ! Sync the cell
      call cell_cart_to_frac(cell)
      call cell_shift_to_unit_cell(cell)
      call cell_frac_to_cart(cell)
    end if
  end if
end subroutine cell_perturb

subroutine cell_supercell(cell, supercell_a, supercell_b, supercell_c)
  use io, only: io_err
  implicit none
  type(cell_type),  intent(inout) ::  cell
  integer,          intent(in)    ::  supercell_a
  integer,          intent(in)    ::  supercell_b
  integer,          intent(in)    ::  supercell_c

  real(kind=dp),  allocatable,  dimension(:,:)  ::  tmp_atom_frac_pos
  real(kind=dp),  allocatable,  dimension(:)    ::  tmp_atom_mass
  integer,        allocatable,  dimension(:)    ::  tmp_atom_species

  integer,                      dimension(3)    ::  nsc, isc
  real(kind=dp),                dimension(3)    ::  dx_frac
  integer ::  icell, jcell, kcell, tmp_natoms
  integer ::  iatom, icomp, old_atom, new_atom, istat

  ! error checking (valid number supercells)
  if ((supercell_a .lt. 1) .or. (supercell_b .lt. 1) .or. (supercell_c .lt. 1)) then
    call io_err("cell_supercell: Supercell size less than 1 along one or more dimensions")
  end if

  nsc(:) = (/ supercell_a, supercell_b, supercell_c /)
  dx_frac(:) = (/ 1.0_dp/real(nsc(1),dp), 1.0_dp/real(nsc(2),dp), 1.0_dp/real(nsc(3),dp) /)

  tmp_natoms = cell%natoms*nsc(1)*nsc(2)*nsc(3)

  ! allocate tmp arrays
  allocate(tmp_atom_frac_pos(3, tmp_natoms), stat=istat)
  if (istat .ne. 0) call io_err("cell_supercell: Could not allocate tmp_atom_frac_pos array")

  allocate(tmp_atom_mass(tmp_natoms), stat=istat)
  if (istat .ne. 0) call io_err("cell_supercell: Could not allocate tmp_atom_mass array")

  allocate(tmp_atom_species(tmp_natoms), stat=istat)
  if (istat .ne. 0) call io_err("cell_supercell: Could not allocate tmp_atom_species array")

  ! expand lattice vectors - icomp runs over a, b, c
  do icomp = 1, 3
    cell%lattice_vectors(:,icomp) = cell%lattice_vectors(:,icomp)*real(nsc(icomp),dp)
  end do

  ! shift everything down to 'corner' of cell
  do iatom = 1, cell%natoms
    do icomp = 1, 3
      cell%atom_frac_pos(icomp, iatom) = cell%atom_frac_pos(icomp, iatom)*dx_frac(icomp)
    end do
  end do

  ! add in extra atoms to tmp arrays
  new_atom = 1
  do icell = 1, nsc(1)
    do jcell = 1, nsc(2)
      do kcell = 1, nsc(3)
        isc(:) = (/ icell-1, jcell-1, kcell-1 /)
        do old_atom = 1, cell%natoms
          do icomp = 1, 3
            tmp_atom_frac_pos(icomp, new_atom) = cell%atom_frac_pos(icomp, old_atom) + real(isc(icomp),dp)*dx_frac(icomp)
          end do
          tmp_atom_mass(new_atom) = cell%atom_mass(old_atom)
          tmp_atom_species(new_atom) = cell%atom_species(old_atom)
          new_atom = new_atom + 1
        end do
      end do
    end do
  end do

  ! give cell more atoms and reallocate to correct size
  cell%natoms = tmp_natoms
  call cell_reallocate(cell)

  ! copy tmp data over - deep copy
  do iatom = 1, cell%natoms
    cell%atom_frac_pos(:, iatom) = tmp_atom_frac_pos(:, iatom)
    cell%atom_mass(iatom) = tmp_atom_mass(iatom)
    cell%atom_species(iatom) = tmp_atom_species(iatom)
  end do

  ! update rest of cell - order important
  call cell_calc_vol(cell)
  call cell_get_recip(cell)
  call cell_frac_to_cart(cell)
  call cell_min_img_init(cell)

  ! deallocate tmp arrays
  deallocate(tmp_atom_frac_pos, stat=istat)
  if (istat .ne. 0) call io_err("cell_supercell: Could not deallocate tmp_atom_frac_pos array")

  deallocate(tmp_atom_mass, stat=istat)
  if (istat .ne. 0) call io_err("cell_supercell: Could not deallocate tmp_atom_mass array")

  deallocate(tmp_atom_species, stat=istat)
  if (istat .ne. 0) call io_err("cell_supercell: Could not deallocate tmp_atom_species array")
end subroutine cell_supercell

!****s* cell/cell_write =======================================================!
! NAME                                                                         !
!   cell_write (PUBLIC)                                                        !
!------------------------------------------------------------------------------!
! ARGUMENTS                                                                    !
!   type(cell_type),    intent(in)  ::  cell                                   !
!                                                                              !
!   logical,  optional, intent(in)  ::  suppress_output_file                   !
!------------------------------------------------------------------------------!
! DESCRIPTION                                                                  !
!   Writes the contents of the cell to file seedname.cell (unless              !
!   suppress_output_file is present and .true., in which case it outputs to    !
!   stdout unit).                                                              !
!                                                                              !
!   If there is no seedname (very unlikely), we instead write to whatever      !
!   stdout is connected to.                                                    !
!==============================================================================!
subroutine cell_write(cell, suppress_output_file)
  use io,         only: seedname, stdout, io_open_file, io_close_file
  use constants,  only: element_symbol
  implicit none
  type(cell_type),    intent(in)  ::  cell
  logical,  optional, intent(in)  ::  suppress_output_file

  integer ::  out_unit, i, j

!+----------------------------------------------------------------------------------------------------------------------+
!|                                                         CELL                                                         |
!+-----------------------------------------------------------+----------------------------------------------------------+
!|                      lattice vectors                      |                reciprocal lattice vectors                |
!|           a                  b                  c         |          a*                 b*                 c*        |
!| fffff.ffffffffffff fffff.ffffffffffff fffff.ffffffffffff  | fffff.ffffffffffff fffff.ffffffffffff fffff.ffffffffffff |
!+-----------------------------------------------------------+----------------------------------------------------------+
!| volume fffffffff.ffffffff                                                                                            |
!+----------------------------------------------------------------------------------------------------------------------+
!| num atoms iiiiiiiiii                                                                                                 |
!+----------+---------------+---------------------------------------------+---------------------------------------------+
!| species  |     mass      |              cartesian position             |             fractional position             |
!|          |               |          x             y             z      |         u             v             w       |
!|   aaa    | ffffff.ffffff |  fffff.fffffff fffff.fffffff fffff.fffffff  |  fff.fffffffff fff.fffffffff fff.fffffffff  |
!+----------+---------------+---------------------------------------------+---------------------------------------------+

  out_unit = stdout

  if ((.not. present(suppress_output_file)) .or. (.not. suppress_output_file)) then
    ! only open file if we have a seedname, regardless of suppress_output_file
    if (len_trim(seedname) .ge. 1) then
      call io_open_file(trim(seedname)//'.cell', out_unit, 'replace')
    end if
  end if

  write(out_unit, "('+', 118('-'), '+')")
  write(out_unit, "('|', 57X, 'CELL', 57X, '|')")
  write(out_unit, "('+', 59('-'), '+', 58('-'), '+')")
  write(out_unit, "('|', 22X, 'lattice vectors', 22X, '|', 16X, 'reciprocal lattice vectors', 16X, '|')")
  write(out_unit, "('|', 11X, 'a', 18X, 'b', 18X, 'c', 9X, '|', 10X, 'a*', 17X, 'b*', 17X, 'c*', 8X, '|')")
  do i = 1, 3
    write(out_unit, "('|', 1X, 3(F18.12, 1X), 1X, '|', 1X, 3(F18.12, 1X), '|')") &
    & (cell%lattice_vectors(i,j), j=1,3), (cell%recip_lattice(i,j), j=1,3)  ! might need to transpose recip_lattice
  end do
  write(out_unit, "('+', 59('-'), '+', 58('-'), '+')")
  write(out_unit, "('|', 1X, 'volume', 1X, F18.8, 92X, '|')") cell%volume
  write(out_unit, "('+', 118('-'), '+')")
  write(out_unit, "('|', 1X, 'num atoms', 1X, I10, 97X, '|')") cell%natoms
  write(out_unit, "('+', 10('-'), '+', 15('-'), '+', 45('-'), '+', 45('-'), '+')")
  write(out_unit, "('|', 1X, 'species', 2X, '|', 5X, 'mass', 6X, '|', 14X, &
  & 'cartesian position', 13X, '|', 13X, 'fractional position', 13X, '|')")
  write(out_unit, "('|', 10X, '|', 15X, '|', 10X, 'x', 13X, 'y', 13X, 'z', 6X, '|', 9X, 'u', 13X, 'v', 13X, 'w', 7X, '|')")
  do i = 1, cell%natoms
    write(out_unit, "('|', 3X, A3, 4X, '| ', F13.6, 1X, '|', 2X, 3(F13.7, 1X), 1X, '|', 2X, 3(F13.9, 1X), 1X, '|')") &
    & element_symbol(cell%atom_species(i)), cell%atom_mass(i), (cell%atom_cart_pos(j,i), j=1,3), (cell%atom_frac_pos(j,i), j=1,3)
  end do
  write(out_unit, "('+', 10('-'), '+', 15('-'), '+', 45('-'), '+', 45('-'), '+')")

  ! only close file if we opened a new one
  if (out_unit .ne. stdout) call io_close_file(out_unit)
end subroutine cell_write


!****s* cell/cell_init ========================================================!
! NAME                                                                         !
!   cell_init (PUBLIC)                                                         !
!------------------------------------------------------------------------------!
! ARGUMENTS                                                                    !
!   type(cell_type),  intent(inout) ::  cell                                   !
!------------------------------------------------------------------------------!
! DESCRIPTION                                                                  !
!   Initializes cell ready for use.                                            !
!                                                                              !
!   Gets lattice vectors, atomic species and fractional positions from input   !
!   (via cell_read_input) then converts everything to atomic units and         !
!   populates the rest of the cell (reciprocal lattice vectors, cartesian      !
!   positions, masses etc).                                                    !
!------------------------------------------------------------------------------!
! NOTES                                                                        !
!   If we allow for atomic positions to be inputted in cartesian form, this    !
!   routine will need updating.                                                !
!==============================================================================!
subroutine cell_init(cell)
  use constants,  only: element_mass, units_natural_to_atomic
  implicit none
  type(cell_type),  intent(inout) ::  cell
  integer ::  i

  cell%natoms = 0

  ! get cell from input - this gets the lattice vectors, atomic species and positions (in fractional coordinates)
  call cell_read_input(cell)

  ! with known atomic species, can populate mass array
  do i = 1, cell%natoms
    ! mass units are in amu, we want to rescale such that electron_mass = 1
    cell%atom_mass(i) = units_natural_to_atomic(mass=element_mass(cell%atom_species(i)))
  end do

  ! rescale lattice vectors from "natural" units to atomic units
  cell%lattice_vectors = units_natural_to_atomic(length=cell%lattice_vectors)

  call cell_calc_vol(cell)

  ! must be done after lattice vectors have had units converted:
  call cell_get_recip(cell)             ! need vol first
  call cell_shift_to_unit_cell(cell)
  call cell_frac_to_cart(cell)
  call cell_min_img_init(cell)

  ! perform perturbation if necessary
  atom_perturb = units_natural_to_atomic(length=atom_perturb)
  lattice_perturb = units_natural_to_atomic(length=lattice_perturb)
  call cell_perturb(cell, atom_perturb, lattice_perturb)
end subroutine cell_init


!****s* cell/cell_allocate ====================================================!
! NAME                                                                         !
!   cell_allocate (PUBLIC)                                                     !
!------------------------------------------------------------------------------!
! ARGUMENTS                                                                    !
!   type(cell_type),  intent(inout) ::  cell                                   !
!------------------------------------------------------------------------------!
! DESCRIPTION                                                                  !
!   Allocates all of the allocatable arrays within the cell given the number   !
!   of atoms contained within the cell. (Not the lattice_vectors and           !
!   recip_lattice arrays, which are statically allocated.)                     !
!------------------------------------------------------------------------------!
! NOTES                                                                        !
!   cell%natoms must be > 0.                                                   !
!                                                                              !
!   None of the arrays must already be allocated.                              !
!==============================================================================!
subroutine cell_allocate(cell)
  use io,  only: io_err
  implicit none
  type(cell_type),  intent(inout) ::  cell
  integer ::  istat

  if (cell%natoms .lt. 1) call io_err("cell_allocate: Cell contains no atoms")

  ! only allocated if not already
  if (.not. allocated(cell%atom_frac_pos)) then
    allocate(cell%atom_frac_pos(3, cell%natoms), stat=istat)
    if (istat .ne. 0) call io_err("cell_allocate: Could not allocate cell%atom_frac_pos")
  else
    call io_err("cell_allocate: cell%atom_frac_pos already allocated")
  end if

  if (.not. allocated(cell%atom_cart_pos)) then
    allocate(cell%atom_cart_pos(3, cell%natoms), stat=istat)
    if (istat .ne. 0) call io_err("cell_allocate: Could not allocate cell%atom_cart_pos")
  else
    call io_err("cell_allocate: cell%atom_cart_pos already allocated")
  end if

  if (.not. allocated(cell%atom_species)) then
    allocate(cell%atom_species(cell%natoms), stat=istat)
    if (istat .ne. 0) call io_err("cell_allocate: Could not allocate cell%atom_species")
  else
    call io_err("cell_allocate: cell%atom_species already allocated")
  end if

  if (.not. allocated(cell%atom_mass)) then
    allocate(cell%atom_mass(cell%natoms), stat=istat)
    if (istat .ne. 0) call io_err("cell_allocate: Could not allocate cell%atom_mass")
  else
    call io_err("cell_allocate: cell%atom_mass already allocated")
  end if
end subroutine cell_allocate


!****s* cell/cell_deallocate ==================================================!
! NAME                                                                         !
!   cell_deallocate (PUBLIC)                                                   !
!------------------------------------------------------------------------------!
! ARGUMENTS                                                                    !
!   type(cell_type),  intent(inout) ::  cell                                   !
!------------------------------------------------------------------------------!
! DESCRIPTION                                                                  !
!   Deallocates all of the allocatable arrays within the cell.                 !
!   (Not the lattice_vectors and recip_lattice arrays, which are statically    !
!   allocated.)                                                                !
!------------------------------------------------------------------------------!
! NOTES                                                                        !
!   None of the arrays must already be unallocated.                            !
!==============================================================================!
subroutine cell_deallocate(cell)
  use io, only: io_err
  implicit none

  type(cell_type),  intent(inout) ::  cell
  integer ::  istat

  if (allocated(cell%atom_frac_pos)) then
    deallocate(cell%atom_frac_pos, stat=istat)
    if (istat .ne. 0) call io_err("cell_deallocate: Could not deallocate cell%atom_frac_pos")
  end if

  if (allocated(cell%atom_cart_pos)) then
    deallocate(cell%atom_cart_pos, stat=istat)
    if (istat .ne. 0) call io_err("cell_deallocate: Could not deallocate cell%atom_cart_pos")
  end if

  if (allocated(cell%atom_species)) then
    deallocate(cell%atom_species, stat=istat)
    if (istat .ne. 0) call io_err("cell_deallocate: Could not deallocate cell%atom_species")
  end if

  if (allocated(cell%atom_mass)) then
    deallocate(cell%atom_mass, stat=istat)
    if (istat .ne. 0) call io_err("cell_deallocate: Could not deallocate cell%atom_mass")
  end if

  cell%natoms = 0   ! might as well
end subroutine cell_deallocate

subroutine cell_reallocate(cell)
  use io, only: io_err
  implicit none

  type(cell_type),  intent(inout) ::  cell
  logical ::  do_alloc
  integer ::  istat

  do_alloc = .true.
  if (allocated(cell%atom_frac_pos)) then
    do_alloc = .false.

    if (size(cell%atom_frac_pos, 2) .ne. cell%natoms) then
      deallocate(cell%atom_frac_pos, stat=istat)
      if (istat .ne. 0) call io_err("cell_reallocate: Could not deallocate cell%atom_frac_pos")
      do_alloc = .true.
    end if
  end if
  if (do_alloc) then
    allocate(cell%atom_frac_pos(3, cell%natoms), stat=istat)
    if (istat .ne. 0) call io_err("cell_reallocate: Could not allocate cell%atom_frac_pos")
  end if


  do_alloc = .true.
  if (allocated(cell%atom_cart_pos)) then
    do_alloc = .false.

    if (size(cell%atom_cart_pos, 2) .ne. cell%natoms) then
      deallocate(cell%atom_cart_pos, stat=istat)
      if (istat .ne. 0) call io_err("cell_reallocate: Could not deallocate cell%atom_cart_pos")
      do_alloc = .true.
    end if
  end if
  if (do_alloc) then
    allocate(cell%atom_cart_pos(3, cell%natoms), stat=istat)
    if (istat .ne. 0) call io_err("cell_reallocate: Could not allocate cell%atom_cart_pos")
  end if


  do_alloc = .true.
  if (allocated(cell%atom_species)) then
    do_alloc = .false.

    if (size(cell%atom_species, 1) .ne. cell%natoms) then
      deallocate(cell%atom_species, stat=istat)
      if (istat .ne. 0) call io_err("cell_reallocate: Could not deallocate cell%atom_species")
      do_alloc = .true.
    end if
  end if
  if (do_alloc) then
    allocate(cell%atom_species(cell%natoms), stat=istat)
    if (istat .ne. 0) call io_err("cell_reallocate: Could not allocate cell%atom_species")
  end if


  do_alloc = .true.
  if (allocated(cell%atom_mass)) then
    do_alloc = .false.

    if (size(cell%atom_mass, 1) .ne. cell%natoms) then
      deallocate(cell%atom_mass, stat=istat)
      if (istat .ne. 0) call io_err("cell_reallocate: Could not deallocate cell%atom_mass")
      do_alloc = .true.
    end if
  end if
  if (do_alloc) then
    allocate(cell%atom_mass(cell%natoms), stat=istat)
    if (istat .ne. 0) call io_err("cell_reallocate: Could not allocate cell%atom_mass")
  end if
end subroutine cell_reallocate


!****s* cell/cell_copy ========================================================!
! NAME                                                                         !
!   cell_copy (PUBLIC)                                                         !
!------------------------------------------------------------------------------!
! ARGUMENTS                                                                    !
!   type(cell_type),  intent(in)    ::  incell                                 !
!   type(cell_type),  intent(inout) ::  outcell                                !
!------------------------------------------------------------------------------!
! DESCRIPTION                                                                  !
!   Copies the contents of incell to outcell on a component by component basis !
!   (a deep copy).                                                             !
!==============================================================================!
subroutine cell_copy(incell, outcell)
  use io, only: io_err
  implicit none
  type(cell_type),  intent(in)    ::  incell
  type(cell_type),  intent(inout) ::  outcell

  outcell%natoms            = incell%natoms
  outcell%volume            = incell%volume
  outcell%lattice_vectors   = incell%lattice_vectors
  outcell%recip_lattice     = incell%recip_lattice
  outcell%min_img_cutoff    = incell%min_img_cutoff

  ! can only allocate new cell (and copy these arrays if natoms > 0)
  if (outcell%natoms .gt. 0) then
    call cell_reallocate(outcell)

    outcell%atom_frac_pos   = incell%atom_frac_pos
    outcell%atom_cart_pos   = incell%atom_cart_pos
    outcell%atom_mass       = incell%atom_mass
    outcell%atom_species    = incell%atom_species
  end if
end subroutine cell_copy


!****s* cell/cell_read_input ==================================================!
! NAME                                                                         !
!   cell_read_input (PRIVATE)                                                  !
!------------------------------------------------------------------------------!
! ARGUMENTS                                                                    !
!   type(cell_type),  intent(inout) ::  cell                                   !
!------------------------------------------------------------------------------!
! DESCRIPTION                                                                  !
!   Sets the lattice vectors, allocates the cell and stores the atomic species !
!   and atomic positions (in fractional space) based on the contents of the    !
!   input hash table.                                                          !
!==============================================================================!
subroutine cell_read_input(cell)
  use io, only: max_line_len, io_input_get_data, io_str_get_num_tokens, io_str_get_token, io_str_to_real, &
    & io_input_get_single_value, io_err
  use constants,  only: element_symbol_to_z
  implicit none
  type(cell_type),  intent(inout) ::  cell

  character(len=max_line_len)                               ::  keyword
  character(len=max_line_len),  allocatable,  dimension(:)  ::  values
  logical ::  found
  integer ::  iline, itoken, species
  integer ::  nlines, ntokens
  integer ::  istat

  ! begin lattice vectors --------------------------------------------------------------------------
  keyword = 'lattice_vectors'
  call io_input_get_data(trim(keyword), values, found)
  if (.not. found) call io_err("cell_read_input: Did not find "//trim(keyword)//" in input file")

  nlines = size(values,1)
  if (nlines .ne. 3)  call io_err("cell_read_input: Expected "//trim(keyword)//" block to be 3 lines")

  do iline = 1, nlines
    ntokens = io_str_get_num_tokens(values(iline))

    if (ntokens .ne. 3) call io_err("cell_read_input: Expected 3 values per line for "//trim(keyword)//" block")
    do itoken = 1, ntokens
      cell%lattice_vectors(itoken, iline) = io_str_to_real(io_str_get_token(values(iline), itoken))
    end do
  end do
  ! end lattice vectors ----------------------------------------------------------------------------

  ! begin atomic postions --------------------------------------------------------------------------
  keyword = 'atomic_positions'
  call io_input_get_data(trim(keyword), values, found)
  if (.not. found) call io_err("cell_read_input: No "//trim(keyword)//" in input file")

  nlines = size(values,1)
  cell%natoms = nlines

  call cell_allocate(cell)

  do iline = 1, nlines
    ntokens = io_str_get_num_tokens(values(iline))
    if (ntokens .ne. 4) call io_err("cell_read_input: Expected 4 values per line for "//trim(keyword)//" block")

    ! first token is atom species
    species = element_symbol_to_z(trim(io_str_get_token(values(iline), 1)))
    if (species .lt. 1) call io_err("cell_read_input: Unexpected atom species")

    cell%atom_species(iline) = species  ! iatom = iline

    do itoken = 2, ntokens !2,4
      ! itoken-1 is coordinate index, iatom = iline
      cell%atom_frac_pos(itoken-1, iline) = io_str_to_real(io_str_get_token(values(iline), itoken))
    end do
  end do
  ! end atomic positions ---------------------------------------------------------------------------

  ! begin cell_perturb ----------------------------------------------------------------------------
  call io_input_get_single_value('cell_perturb_lattice', lattice_perturb, -1.0_dp)
  call io_input_get_single_value('cell_perturb_atoms', atom_perturb, -1.0_dp)
  ! end cell_perturb ------------------------------------------------------------------------------

  ! finalization - should always be allocated to get this far, but no harm in checking...
  if (allocated(values)) then
    deallocate(values, stat=istat)
    if (istat .ne. 0) call io_err("cell_read_input: Could not deallocate values array")
  end if
end subroutine cell_read_input


!****s* cell/cell_shift_to_unit_cell ==========================================!
! NAME                                                                         !
!   cell_shift_to_unit_cell (PUBLIC)                                           !
!------------------------------------------------------------------------------!
! ARGUMENTS                                                                    !
!   type(cell_type),  intent(inout) ::  cell                                   !
!------------------------------------------------------------------------------!
! DESCRIPTION                                                                  !
!   Enforces periodic boundary conditions such that anything leaving one side  !
!   of the cell will return through the other side.                            !
!------------------------------------------------------------------------------!
! NOTES                                                                        !
!   This operation is carried out in fractional space, therefore the           !
!   fractional positions must be correct. The result is that the fractional    !
!   positions are always shifted into the unit cell.                           !
!==============================================================================!
subroutine cell_shift_to_unit_cell(cell)
  implicit none
  type(cell_type),  intent(inout) ::  cell

  ! enforces periodic boundary conditions
  ! anything leaving from one side of the cell will return through the other side
  cell%atom_frac_pos(:,:) = cell%atom_frac_pos(:,:) - real(floor(cell%atom_frac_pos(:,:)),dp)
end subroutine cell_shift_to_unit_cell


!****s* cell/cell_frac_to_cart ================================================!
! NAME                                                                         !
!   cell_frac_to_cart (PUBLIC)                                                 !
!------------------------------------------------------------------------------!
! ARGUMENTS                                                                    !
!   type(cell_type),  intent(inout) ::  cell                                   !
!------------------------------------------------------------------------------!
! DESCRIPTION                                                                  !
!   Updates the cartesian atomic positions from the fractional positions.      !
!==============================================================================!
subroutine cell_frac_to_cart(cell)
  implicit none
  type(cell_type),  intent(inout) ::  cell

  cell%atom_cart_pos=matmul(cell%lattice_vectors, cell%atom_frac_pos)
end subroutine cell_frac_to_cart


!****s* cell/cell_cart_to_frac=================================================!
! NAME                                                                         !
!   cell_cart_to_frac (PUBLIC)                                                 !
!------------------------------------------------------------------------------!
! ARGUMENTS                                                                    !
!   type(cell_type),  intent(inout) ::  cell                                   !
!------------------------------------------------------------------------------!
! DESCRIPTION                                                                  !
!   Updates the fractional atomic positions from the cartesian positions.      !
!==============================================================================!
subroutine cell_cart_to_frac(cell)
  use constants, only: two_pi
  implicit none
  type(cell_type),  intent(inout) ::  cell

  cell%atom_frac_pos=matmul(cell%recip_lattice/two_pi, cell%atom_cart_pos) ! need divide by two_pi because of definition..
end subroutine cell_cart_to_frac


!****s* cell/cell_get_recip ===================================================!
! NAME                                                                         !
!   cell_get_recip (PUBLIC)                                                    !
!------------------------------------------------------------------------------!
! ARGUMENTS                                                                    !
!   type(cell_type),  intent(inout) ::  cell                                   !
!------------------------------------------------------------------------------!
! DESCRIPTION                                                                  !
!   Calculates the reciprocal lattice vectors from the lattice vectors.        !
!==============================================================================!
subroutine cell_get_recip(cell)
  use constants,  only: two_pi
  use algor,      only: algor_cross_product
  implicit none
  type(cell_type),  intent(inout) ::  cell
  real(kind=dp)                   ::  two_pi_over_vol

  ! calc recip_lattice:
  two_pi_over_vol = two_pi/cell%volume  ! might not want the two_pi here...

  cell%recip_lattice(:,1) = two_pi_over_vol*algor_cross_product(cell%lattice_vectors(:,2), cell%lattice_vectors(:,3))
  cell%recip_lattice(:,2) = two_pi_over_vol*algor_cross_product(cell%lattice_vectors(:,3), cell%lattice_vectors(:,1))
  cell%recip_lattice(:,3) = two_pi_over_vol*algor_cross_product(cell%lattice_vectors(:,1), cell%lattice_vectors(:,2))

  ! transpose matrix
  cell%recip_lattice(:,:) = transpose(cell%recip_lattice(:,:))
end subroutine cell_get_recip


!****s* cell/cell_min_img_init ================================================!
! NAME                                                                         !
!   cell_min_img_init (PUBLIC)                                                 !
!------------------------------------------------------------------------------!
! ARGUMENTS                                                                    !
!   type(cell_type),  intent(inout) ::  cell                                   !
!------------------------------------------------------------------------------!
! DESCRIPTION                                                                  !
!   Initializes cell_min_img - on first pass we need to calculate all of the   !
!   "c_vectors" (these are used for all possible corners of adjacent cells)    !
!   - this is indepedent of any one particular cell.                           !
!                                                                              !
!   Next we calculate the perpendicular cell widths and store the square of    !
!   half of minimum width, which is used as a cutoff for cell_min_img          !
!   - anything within this distance is guaranteed to be the minimum image.     !
!==============================================================================!
subroutine cell_min_img_init(cell)
  use algor, only: algor_cross_product
  implicit none
  type(cell_type),      intent(inout) ::  cell

  integer,  parameter,  dimension(3)  ::  ck_values = (/ -1, 0, 1 /)
  real(kind=dp),        dimension(3)  ::  axb, bxc, cxa
  real(kind=dp)                       ::  w_a, w_b, w_c
  integer                             ::  i, j, k, l

  ! this part is universal for all cells - only do this once, ever:
  if (.not. min_img_initialized) then
    ! fill up the matrix with all of the possible c_vectors
    ! there will be 27 possible c vectors
    l = 1
    do i = 1,3
      do j = 1,3
        do k = 1,3
          min_img_c_vectors(:,l) = (/ ck_values(i), ck_values(j), ck_values(k) /)
          l = l+1
        end do
      end do
    end do
    min_img_initialized = .true.
  end if

  ! get perpendicular cell widths
  bxc = algor_cross_product(cell%lattice_vectors(:,2), cell%lattice_vectors(:,3))
  cxa = algor_cross_product(cell%lattice_vectors(:,3), cell%lattice_vectors(:,1))
  axb = algor_cross_product(cell%lattice_vectors(:,1), cell%lattice_vectors(:,2))
  w_a = abs(dot_product(cell%lattice_vectors(:,1), bxc))/sqrt(dot_product(bxc,bxc))
  w_b = abs(dot_product(cell%lattice_vectors(:,2), cxa))/sqrt(dot_product(cxa,cxa))
  w_c = abs(dot_product(cell%lattice_vectors(:,3), axb))/sqrt(dot_product(axb,axb))
  cell%min_img_cutoff = 0.5_dp*min(w_a, w_b, w_c)
  cell%min_img_cutoff = cell%min_img_cutoff*cell%min_img_cutoff ! square to save a sqrt in cell_min_img
end subroutine cell_min_img_init


!****s* cell/cell_min_img =====================================================!
! NAME                                                                         !
!   cell_min_img (PUBLIC)                                                      !
!------------------------------------------------------------------------------!
! ARGUMENTS                                                                    !
!   type(cell_type),              intent(in)              ::  cell             !
!   integer,                      intent(in)              ::  iatom, jatom     !
!   real(kind=dp),  dimension(3), intent(out)             ::  rij_out          !
!                                                                              !
!   real(kind=dp),                intent(out),  optional  ::  distance         !
!   character(len=4),             intent(in),   optional  ::  space            !
!------------------------------------------------------------------------------!
! DESCRIPTION                                                                  !
!   Returns the vector between the closest image of atoms under periodic       !
!   boundary conditions.                                                       !
!                                                                              !
!   Optionally we can also return the length of this vector (which must be in  !
!   cartesian space) - and/or change the space that the vector between the     !
!   atoms is defined in (default is cartesian).                                !
!==============================================================================!
subroutine cell_min_img(cell, iatom, jatom, rij_out, distance, space)
  use io, only: io_err
  implicit none
  type(cell_type),              intent(in)              ::  cell
  integer,                      intent(in)              ::  iatom, jatom
  real(kind=dp),  dimension(3), intent(out)             ::  rij_out
  real(kind=dp),                intent(out),  optional  ::  distance
  character(len=4),             intent(in),   optional  ::  space         ! either 'frac' or 'cart', space to output rij_out

  character(len=4)                ::  int_space
  real(kind=dp)                   ::  int_dist
  real(kind=dp),    dimension(3)  ::  ri, rj

  if (iatom .eq. jatom) then
    rij_out = 0.0_dp
    if (present(distance)) distance = 0.0_dp
    return
  end if

  int_space = 'cart'
  if (present(space)) then
    int_space = space
  end if

  select case (int_space)
    case ('cart')
      ri = cell%atom_cart_pos(:,iatom)
      rj = cell%atom_cart_pos(:,jatom)
    case ('frac')
      ri = cell%atom_frac_pos(:,iatom)
      rj = cell%atom_frac_pos(:,jatom)
    case default
      call io_err("cell_min_img: Unrecognised space")
  end select

  call cell_min_img_vec(cell, ri, rj, rij_out, int_dist, int_space)
  if (present(distance)) distance = int_dist
end subroutine cell_min_img

subroutine cell_min_img_vec(cell, ri, rj, rij_out, distance, space)
  use io,         only: io_err
  use constants,  only: two_pi
  ! 'space' variable determines input space of ri and rj, as well as output space of rij_out
  ! these are in cartesian space by default (ie: unless space == 'frac')
  ! distance is optional, real space distance
  implicit none
  type(cell_type),              intent(in)              ::  cell
  real(kind=dp),  dimension(3), intent(in)              ::  ri, rj
  real(kind=dp),  dimension(3), intent(out)             ::  rij_out
  real(kind=dp),                intent(out),  optional  ::  distance
  character(len=4),             intent(in),   optional  ::  space         ! either 'frac' or 'cart'

  character(len=4)                          ::  int_space
  real(kind=dp),  dimension(3)              ::  int_ri, int_rj, rij, fij, rij_min
  real(kind=dp)                             ::  smallestdist, currentdist
  integer                                   ::  k, l, cmin

  ! for speedup - just care about 8 vectors, not all 27
  logical,        dimension(27)             ::  is_c_important
  integer,        dimension(3)              ::  disallowed_c

  int_space = 'cart'
  if (present(space)) then
    int_space = space
  end if

  ! don't check whether ri == rj in vector routine
  ! set up internal position variables
  select case (int_space)
    case ('cart')
      int_ri = ri
      int_rj = rj
      fij = ri - rj                                 ! fij = rij
      fij = matmul(cell%recip_lattice/two_pi, fij)  ! convert rij to fractionals
    case ('frac')
      ! convert ri, rj to cartesian
      int_ri = matmul(cell%lattice_vectors, ri)
      int_rj = matmul(cell%lattice_vectors, rj)
      fij = ri - rj
    case default
      call io_err("cell_min_img_vec: Unrecognised space")
  end select

  ! fij = fi - fj (frac coords)
  ! rij = ri - rj (cart coords)
  rij = int_ri - int_rj

  rij_min = rij

  ! calc distance**2 for within unit cell
  currentdist = dot_product(rij,rij)

  ! if sqrt of this is less than smallest half width then we have minimum image already
  ! (cell%min_img_cutoff is this width**2)
  if (currentdist .gt. cell%min_img_cutoff) then
    ! we need to do more work...

    ! first pick the allowed c_vectors from the list (can use all, but this saves comparing all 27 vectors)
    ! disallowed_c(k) is basically just the sign of the s_k
    ! (c_k can be either 0 or -SIGN(s_k) - doesn't matter what we choose if we have s_k = 0 due to symmetry)
    do k = 1,3
      if (fij(k) .ge. 0.0_dp) then
        disallowed_c(k) = 1
      else
        ! negative
        disallowed_c(k) = -1
      end if
    end do

    ! make an array with true and false to represent if we care about this vector..
    ! initialize to .true. and make .false. if we don't care..
    is_c_important = .true.
    do l = 1,27
      do k = 1,3
        if (min_img_c_vectors(k,l) .eq. disallowed_c(k)) then
          is_c_important(l) = .false.
          exit    ! entire vector useless, proceed to next vector...
        end if
      end do
    end do

    smallestdist = huge(1.0_dp) ! big..
    do l = 1,27
      if (is_c_important(l)) then

        ! set rij to be ri... (vary ri and use same rj to find min image)
        rij = int_ri

        do k = 1,3
          ! at first, rij is ri... so this will give
          ! rij = ri + sum_k(c_k * H_k)     {H_k is k-th column of lattice_vectors, ie: a, b, c}
          rij(:) = rij(:) + real(min_img_c_vectors(k,l),dp)*cell%lattice_vectors(:,k)
        end do
        ! rij = ri - rj + sum_k(c_k * H_k)
        rij = rij - int_rj
        currentdist = dot_product(rij,rij)

        ! store smallest distance
        if (currentdist .lt. smallestdist) then
          smallestdist = currentdist
          rij_min = rij
          cmin = l                                ! index of smallest c in min_img_c_vectors
        end if
      end if
    end do

    currentdist = smallestdist
    fij(:) = fij(:) + real(min_img_c_vectors(:,cmin),dp)  ! use the minimum c_vector
  end if

  select case (int_space)
    case ('cart')
      rij_out = rij_min
    case ('frac')
      rij_out = fij
  end select

  if (present(distance)) distance = sqrt(currentdist)
end subroutine cell_min_img_vec


subroutine cell_calc_vol(cell)
  ! could be a function
  use algor,  only: algor_cross_product
  implicit none
  type(cell_type),  intent(inout) ::  cell

  ! vol = a .dot. (b x c)
  cell%volume = dot_product(cell%lattice_vectors(:,1), algor_cross_product(cell%lattice_vectors(:,2),cell%lattice_vectors(:,3)))
end subroutine cell_calc_vol
end module cell
