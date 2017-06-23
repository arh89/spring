module model_fourier_pes

  use constants,  only: dp
  implicit none

  private

  real(kind=dp),  parameter   ::  vtop = 0.3_dp
  real(kind=dp),  parameter   ::  vbridge = 0.1_dp
  real(kind=dp),  parameter   ::  vfcc = -0.1_dp

  real(kind=dp),  parameter   ::  spring_const = 1.0_dp

  ! Must supply 3 of vtop, vbridge, vfcc, vhcp:
  ! vtop    = 3*(a1 + a2)
  ! vbridge = 3*a2 - a1
  ! vfcc    = -(3/2)*(a1 + a2 + sqrt(3)*b1)
  ! vhcp    = -(3/2)*(a1 + a2 - sqrt(3)*b1)

  public :: model_fourier_pes_get_forces

contains

subroutine model_fourier_pes_get_forces(cell, potential)
  use constants,  only: two_pi
  use cell,       only: cell_type
  implicit none

  type(cell_type),  intent(inout) ::  cell
  real(kind=dp),    intent(out)   ::  potential

  real(kind=dp) ::  displacement
  real(kind=dp) ::  interaction
  integer       ::  iatom

  ! We assume that potential energy is just the sum of all the interactions
  ! and atoms do not interact with one another - only the surface.

  ! We also assume that lattice vectors are correct on entry to routine:
  !   That the surface is in xy plane (a and b define this plane)
  !   c is only along z (ie: no components in x and y)
  !   a and b are at 60 degrees to one another with correct length..
  !   atom positions are mapped back to unit cell

  potential = 0.0_dp
  cell%atom_force(:,:) = 0.0_dp

  do iatom = 1, cell%natoms

    ! first we deal with the xy components of the force - for this we need the energy surface:
    call model_fourier_potential_2d(cell, cell%atom_cart_pos(1:2, iatom, cell%atom_force(1:2, iatom))

    ! cell%atom_force(3, iatom) = 0

    ! next we need the harmonic contribution for force along z...
    ! assume optimal height is at base of cell (0, 0, 0)
    displacement = cell%atom_cart_pos(3,iatom)                        ! z - 0.0.. always positive

    ! if in top half of cell.. interact with top surface
    if (cell%atom_frac_pos(3, iatom) .gt. 0.5_dp) &
      & displacement = -(cell%lattice_vectors(3,3) - displacement)    ! make sure sign correct

    cell%atom_force(3,iatom) = -spring_const*displacement

    interaction = interaction + 0.5_dp*spring_const*displacement**2

    potential = potential + interaction
  end do
end subroutine fourier_pes_get_forces

subroutine model_fourier_potential_2d(cell, pos, potential, force)
  implicit none
  type(cell_type),                          intent(in)  ::  cell
  real(kind=dp),    dimension(2),           intent(in)  ::  pos       ! along x and y
  real(kind=dp),                            intent(out) ::  potential
  real(kind=dp),    dimension(2), optional, intent(out) ::  force     ! along x and y

  ! Fourier coefficients:
  real(kind=dp),  parameter   ::  a1 = (vtop - vbridge)/4.0_dp
  real(kind=dp),  parameter   ::  a2 = (3.0_dp*vbridge + vtop)/12.0_dp
  real(kind=dp),  parameter   ::  b1 = -(3.0_dp*vtop + 2.0_dp*vfcc/3.0_dp)/sqrt(3.0_dp)

  real(kind=dp) ::  c1, c2, c3, s1, s2, s3
  real(kind=dp) ::  g1, g2, g3

  ! Note: MUST have correct lattice vectors or potential is ill-defined
  ! gamma = 4.0_dp*pi/(sqrt(3.0_dp)*lattice_const)
  g1x = cell%recip_lattice(1,1) ! gamma * 1.0_dp
  g1y = cell%recip_lattice(1,2) ! gamma * 0.0_dp
  g2x = cell%recip_lattice(2,1) ! gamma * 0.5_dp                ! gamma*cos(pi/3.0_dp)
  g2y = cell%recip_lattice(2,2) ! gamma * 0.5_dp*sqrt(3.0_dp)   ! gamma*sin(pi/3.0_dp))

  ! g .dot. r
  g1 =  g1x*pos(1) + g1y*pos(2)
  g2 =  g2x*pos(1) + g2y*pos(2)
  g3 = -g2x*pos(1) + g2y*pos(2)

  ! save a few trig functions
  c1 = cos(g1)
  c2 = cos(g2)
  c3 = cos(g3)
  s1 = sin(g1)
  s2 = sin(g2)
  s3 = sin(g3)

  potential = a1*(c1 + c2 + c3) + 2.0_dp*a2*(c1**2 + c2**2 + c3**2 - 1.5_dp) + b1*(s1 - s2 + s3)

  if (present(force)) then
    force(1) = -a1*(-g1x*s1 - g2x*(s2 - s3)) - 4.0_dp*a2(-g1x*s1*c1 - g2x*(s2*c2 - s3*c3)) - b1*(g1x*c1 -g2x*(c2 + c3))
    force(2) = -a1*(-g1y*s1 - g2y*(s2 + s3)) - 4.0_dp*a2(-g1y*s1*c1 - g2y*(s2*c2 + s3*c3)) - b1*(g1y*c1 - g2y*(c2 - c3))
  end if
end subroutine model_fourier_potential_2d
end module model_fourier_pes
