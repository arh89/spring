module fourier_pes
  use constants,              only: dp
  use cell,                   only: cell_type
  use fourier_interpolation,  only: interp_type

  implicit none

  private

  integer, parameter :: default_supercell_size = 1 ! not a supercell..

  type  ::  pes_type
    type(interp_type) ::  energies
    type(interp_type) ::  spring_consts
    type(interp_type) ::  heights

    ! sometimes need supercell expansion of PES in order to reduce image interactions
    ! only isotropic expansions (in a/b) supported for now...
    integer ::  supercell_size

    ! associated with a particular cell.. store inv_metric here..
    ! cell vectors must not change or this needs recalculating..
    ! (no routine coded - variable cell PES does not make sense..)
    real(kind=dp),  dimension(3,3)  ::  cell_inv_metric
  end type pes_type

  type(pes_type),   save  ::  pes

  public  ::  fourier_pes_init
  public  ::  fourier_pes_get_potential
  public  ::  fourier_pes_get_forces
  public  ::  fourier_pes_finalize

contains

subroutine fourier_pes_init(cell)
  use constants,              only: units_natural_to_atomic
  use algor,                  only: algor_calc_inverse_real
  use io,                     only: io_err, io_input_get_single_value
  use fourier_interpolation,  only: fourier_interpolation_read, fourier_interpolation_transform
  implicit none
  type(cell_type),  intent(in)  ::  cell
  real(kind=dp) ::  ev_to_hartree, ang_to_bohr
  integer       ::  istat

  ! make sure cell contains at least one atom, otherwise PES is useless
  if (cell%natoms .lt. 1) call io_err("fourier_pes_init: Cell does not contain any atoms")

  ! check that the surface is in the xy plane
  if ((abs(cell%lattice_vectors(3,1)) .gt. epsilon(1.0_dp)) .or.  &         ! z component of a
    & (abs(cell%lattice_vectors(3,2)) .gt. epsilon(1.0_dp)))      &         ! z component of b
    & call io_err("fourier_pes_init: surface not in xy plane")

  ! and that c is along z
  if ((abs(cell%lattice_vectors(1,3)) .gt. epsilon(1.0_dp)) .or.  &         ! x component of c
    & (abs(cell%lattice_vectors(2,3)) .gt. epsilon(1.0_dp)))      &         ! y component of c
    & call io_err("fourier_pes_init: c is not aligned along z")

  ! read from input file - can't use init routine as need to convert units
  ! fourier_interpolation_read will exit if keyword does not exist
  call fourier_interpolation_read('fourier_2d_pes_energy', pes%energies)
  call fourier_interpolation_read('fourier_2d_pes_spring', pes%spring_consts)
  call fourier_interpolation_read('fourier_2d_pes_height', pes%heights)

  ! see if we have a supercell
  call io_input_get_single_value('pes_supercell_size', pes%supercell_size, default_supercell_size)
  if (pes%supercell_size .lt. 1) call io_err("fourier_pes_init: Supercell size must be >= 1")

  ! check all grids are 2D
  if ((pes%energies%ndims .ne. 2) .or.      &
  &   (pes%spring_consts%ndims .ne. 2) .or. &
  &   (pes%heights%ndims .ne. 2))           &
  & call io_err("fourier_pes_init: Fourier PES for energies, spring constants "//&
  & "and heights must all be defined on 2D grids")

  ! convert units from user to atomic units
  ev_to_hartree = units_natural_to_atomic(energy=1.0_dp)
  ang_to_bohr = units_natural_to_atomic(length=1.0_dp)

  ! energies -> eV to Hartrees
  pes%energies%grid_values(:,:) = pes%energies%grid_values(:,:)*ev_to_hartree

  ! spring consts -> eV Ang^-2 to Hartree Bohr^-2
  pes%spring_consts%grid_values(:,:) = pes%spring_consts%grid_values(:,:)*ev_to_hartree/(ang_to_bohr**2)

  ! heights -> Ang to Bohr
  pes%heights%grid_values(:,:) = pes%heights%grid_values(:,:)*ang_to_bohr

  ! make sure that no heights are negative
  if (minval(pes%heights%grid_values(:,:)) .lt. 0.0_dp) then
    call io_err("fourier_pes_init: pes%heights contains negative heights")
  end if

  ! check that max height is not larger than cell height
  ! (know that c is along z)
  if (maxval(pes%heights%grid_values(:,:)) .gt. cell%lattice_vectors(3,3)) then
    call io_err("fourier_pes_init: pes%heights contains heights greater than c(z)")
  end if

  call fourier_interpolation_transform(pes%energies)
  call fourier_interpolation_transform(pes%spring_consts)
  call fourier_interpolation_transform(pes%heights)

  ! can deallocate grid_values for each of the grids
  ! everything should be allocated at this point (else transform would not work)
  deallocate(pes%energies%grid_values, stat=istat)
  if (istat .ne. 0) call io_err("fourier_pes_init: Could not deallocate pes%energies%grid_values array")

  deallocate(pes%spring_consts%grid_values, stat=istat)
  if (istat .ne. 0) call io_err("fourier_pes_init: Could not deallocate pes%spring_consts%grid_values array")

  deallocate(pes%heights%grid_values, stat=istat)
  if (istat .ne. 0) call io_err("fourier_pes_init: Could not deallocate pes%heights%grid_values array")

  ! calculate inv_metric
  ! g = h^T h
  pes%cell_inv_metric = matmul(transpose(cell%lattice_vectors), cell%lattice_vectors)
  call algor_calc_inverse_real(pes%cell_inv_metric)
end subroutine fourier_pes_init

subroutine fourier_pes_get_potential(cell, potential)
  use fourier_interpolation,  only: fourier_interpolation_potential
  implicit none
  type(cell_type),  intent(in)  ::  cell
  real(kind=dp),    intent(out) ::  potential

  ! local vars:
  real(kind=dp),  dimension(3)  ::  atom_frac_pos
  real(kind=dp) ::  interaction_energy, equilib_height, displacement, spring_const
  real(kind=dp) ::  half_cell_height
  integer       ::  iatom, j

  ! assume potential energy is just the sum of all the interactions...
  ! and atoms do not interact with one another - only the surface

  half_cell_height = 0.5_dp*cell%lattice_vectors(3,3)   ! c along z only

  potential = 0.0_dp

  do iatom = 1, cell%natoms
    atom_frac_pos(:) = cell%atom_frac_pos(:,iatom)

    ! now rescale coords to map to smaller cell
    do j = 1, 2
      atom_frac_pos(j) = atom_frac_pos(j)*real(pes%supercell_size,dp)
      ! and map back in range 0..1 (shouldn't matter..)
      atom_frac_pos(j) = atom_frac_pos(j) - real(floor(atom_frac_pos(j)),dp)
    end do

    ! get interaction energy, equilibrium height and spring const at this point (a, b) - ie: fractional pos
    call fourier_interpolation_potential(pes%energies, atom_frac_pos(1:2), interaction_energy)
    call fourier_interpolation_potential(pes%heights, atom_frac_pos(1:2), equilib_height)
    call fourier_interpolation_potential(pes%spring_consts, atom_frac_pos(1:2), spring_const)

    ! next add the harmonic contribution..
    ! cell%atom_cart_pos is always mapped back to unit cell (so no negative numbers)
    displacement = cell%atom_cart_pos(3,iatom) - equilib_height     ! c along z

    ! make sure we interact with correct image of surface
    if (abs(displacement) .gt. half_cell_height)  then
      if (displacement .gt. 0.0_dp) then
        displacement = displacement - cell%lattice_vectors(3,3)
      else
        displacement = displacement + cell%lattice_vectors(3,3)
      end if
    end if

    interaction_energy = interaction_energy + 0.5_dp*spring_const*displacement**2
    potential = potential + interaction_energy
  end do
end subroutine fourier_pes_get_potential

subroutine fourier_pes_get_forces(cell, force, potential)
  use fourier_interpolation,  only: fourier_interpolation_potential
  implicit none
  type(cell_type),                              intent(in)  ::  cell
  real(kind=dp),    dimension(3, cell%natoms),  intent(out) ::  force
  real(kind=dp),                                intent(out) ::  potential

  ! local vars:
  real(kind=dp),  dimension(3)  ::  atom_frac_pos
  real(kind=dp),  dimension(2)  ::  energy_grad, spring_grad, height_grad
  real(kind=dp) ::  interaction_energy, equilib_height, displacement, spring_const
  real(kind=dp) ::  half_cell_height
  integer       ::  iatom, j

  ! assume potential energy is just the sum of all the interactions...
  ! and atoms do not interact with one another - only the surface

  half_cell_height = 0.5_dp*cell%lattice_vectors(3,3)   ! c along z only

  potential = 0.0_dp
  force(:,:) = 0.0_dp

  do iatom = 1, cell%natoms
    atom_frac_pos(:) = cell%atom_frac_pos(:,iatom)

    ! now rescale coords to map to smaller cell as necessary
    do j = 1, 2
      atom_frac_pos(j) = atom_frac_pos(j)*real(pes%supercell_size,dp)
      ! and map back in range 0..1 (shouldn't matter..)
      atom_frac_pos(j) = atom_frac_pos(j) - real(floor(atom_frac_pos(j)),dp)
    end do

    ! get interaction energy, equilibrium height and spring const at this point (a, b) - ie: fractional pos
    call fourier_interpolation_potential(pes%energies, atom_frac_pos(1:2), interaction_energy, energy_grad)
    call fourier_interpolation_potential(pes%heights, atom_frac_pos(1:2), equilib_height, height_grad)
    call fourier_interpolation_potential(pes%spring_consts, atom_frac_pos(1:2), spring_const, spring_grad)

    ! get displacement for the harmonic contribution..
    ! cell%atom_cart_pos is always mapped back to unit cell (so no negative numbers)
    displacement = cell%atom_cart_pos(3,iatom) - equilib_height     ! c along z

    ! make sure we interact with correct image of surface and with force in correct direction...
    if (abs(displacement) .gt. half_cell_height)  then
      if (displacement .gt. 0.0_dp) then
        displacement = displacement - cell%lattice_vectors(3,3)
      else
        displacement = displacement + cell%lattice_vectors(3,3)
      end if
    end if

    ! force due to change in min_energy_surface
    force(1:2,iatom) = -energy_grad

    ! now need to add in components due to spatial change in equilibrium height and spring const
    ! change in spring const term:
    force(1:2,iatom) = force(1:2,iatom) - 0.5_dp*(displacement**2)*spring_grad

    ! change in height term:
    force(1:2,iatom) = force(1:2,iatom) + spring_const*displacement*height_grad

    ! forces are along the a and b directions at this point
    ! now shift forces to be in cartesian coordinate system.. ! force(3,iatom) = 0
    force(:,iatom) = matmul(pes%cell_inv_metric, force(:,iatom))
    force(:,iatom) = matmul(cell%lattice_vectors, force(:,iatom))

    ! finally add in z component:
    force(3,iatom) = -spring_const*displacement

    interaction_energy = interaction_energy + 0.5_dp*spring_const*displacement**2
    potential = potential + interaction_energy
  end do
end subroutine fourier_pes_get_forces

subroutine fourier_pes_finalize
  use fourier_interpolation,  only: fourier_interpolation_deallocate
  implicit none

  call fourier_interpolation_deallocate(pes%energies)
  call fourier_interpolation_deallocate(pes%spring_consts)
  call fourier_interpolation_deallocate(pes%heights)
end subroutine fourier_pes_finalize

end module fourier_pes
