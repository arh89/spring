module lj
  use constants,  only: dp
  use io,         only: io_err
  use cell,       only: cell_type
  implicit none

  private

  ! Defaults from http://www.physics.buffalo.edu/phy411-506/topic1/topic1-lec2.pdf
  real(kind=dp),  parameter ::  lj_epsilon =   0.0003785_dp
  real(kind=dp),  parameter ::  lj_sigma   =   6.425_dp
  real(kind=dp),  parameter ::  rcut       = lj_sigma * 2.5_dp

  public  ::  lj_init
  public  ::  lj_get_potential
  public  ::  lj_get_forces
  public  ::  lj_finalize

contains

subroutine lj_init(cell)
  implicit none
  type(cell_type),  intent(in) ::  cell

  if(any(cell%atom_species /= 18)) call io_err("Lennard-Jones only works for Argon, sorry to be 'that guy'...")

end subroutine lj_init

subroutine lj_finalize()
  implicit none

  ! twiddle thumbs

end subroutine lj_finalize

subroutine lj_get_potential(cell, potential)
  implicit none
  type(cell_type),  intent(in)  ::  cell
  real(kind=dp),    intent(out) ::  potential

  real(kind=dp),  dimension(:,:), allocatable :: tmp_force
  integer :: ierr

  allocate(tmp_force(3,cell%natoms), stat=ierr)
  if(ierr /= 0) call io_err("Unable to allocate tmp_force")

  call lj_get_forces(cell, tmp_force, potential)

  deallocate(tmp_force, stat=ierr)
  if(ierr /= 0) call io_err("Unable to deallocate tmp_force")

end subroutine lj_get_potential

subroutine lj_get_forces(cell, force, potential)
  use cell, only: cell_min_img
  implicit none
  type(cell_type),                            intent(in)  ::  cell
  real(kind=dp),    dimension(3,cell%natoms), intent(out) ::  force
  real(kind=dp),                              intent(out) ::  potential

  real(kind=dp), dimension(4) ::  rij
  real(kind=dp), dimension(3) ::  force_contrib
  real(kind=dp)               ::  lj_raw
  real(kind=dp)               ::  pot_contrib
  real(kind=dp)               ::  force_shift
  real(kind=dp)               ::  pot_shift
  real(kind=dp)               ::  sigma_over_r_6, sigma_over_r_12
  real(kind=dp)               ::  sigma_over_rcut_6, sigma_over_rcut_12
  integer                     ::  iatom, jatom

  sigma_over_rcut_6 = (lj_sigma/rcut)**6
  sigma_over_rcut_12 = sigma_over_rcut_6**2

  ! v(r) = v_lj(r) + force_shift*r + pot_shift
  ! ie:
  ! force_shift is the shift required to make force 0 at cutoff (r=rcut)
  ! pot_shift is the shift required to make potential 0 at cutoff and includes a contribution due to presence of force_shift
  force_shift = 24.0_dp*lj_epsilon*(2.0_dp*sigma_over_rcut_12 - sigma_over_rcut_6)/rcut
  pot_shift = 4.0_dp*lj_epsilon*(7.0_dp*sigma_over_rcut_6 - 13.0_dp*sigma_over_rcut_12)

  potential = 0.0_dp
  force(:,:) = 0.0_dp !(3,natoms)

  do iatom = 1, cell%natoms-1
    do jatom = iatom+1, cell%natoms

      ! get distance between i, j
      call cell_min_img(cell, iatom, jatom, rij(1:3), distance=rij(4))

      ! defined such that force and potential are non-zero only for r<rcut
      if (rij(4) .lt. rcut) then
        sigma_over_r_6 = (lj_sigma/rij(4))**6
        sigma_over_r_12 = sigma_over_r_6**2

        lj_raw = 4.0_dp*lj_epsilon*(sigma_over_r_12 - sigma_over_r_6)
        pot_contrib = lj_raw + force_shift*rij(4) + pot_shift           ! contribution to potential due to ij pair
        potential = potential + pot_contrib

        ! make rij a unit vector
        rij(1:3) = rij(1:3)/rij(4)

        ! need to check direction here.. and sign of force_shift
        force_contrib(:) = (24.0_dp*lj_epsilon*(2.0_dp*sigma_over_r_12 - sigma_over_r_6)/rij(4) + force_shift)*rij(1:3)
        force(:,iatom) = force(:,iatom) + force_contrib(:)
        force(:,jatom) = force(:,jatom) - force_contrib(:)    ! Newton's third law..
      end if
    end do
  end do
end subroutine lj_get_forces
end module lj
