module potential 
  ! EAM 6
  use algor
  use constants,  only: dp
  use cell,       only: cell_type
  implicit none

  public

  character(len=10),  parameter ::  potential_type  = 'eam'      ! eam or harmonic
  logical,            parameter ::  non_local_correction = .false.
  logical,            save      ::  pot_initialized = .false.

  ! EAM potential:
  ! Variable names are for first use - when no longer required, pointers (correct name) point to these arrays for storage
  ! For both Ni and H:
  real(kind=dp),  dimension(:,:,:), allocatable           ::  rij_vec
  real(kind=dp),  dimension(:,:),   allocatable,  target  ::  rij_mag       ! seperate from rij_vec array for efficiency...
  real(kind=dp),  dimension(:),     allocatable,  target  ::  rho_bar
  real(kind=dp),  dimension(:,:),   allocatable           ::  d_rho_atomic
  real(kind=dp),  dimension(:,:),                 pointer ::  d_phi
  real(kind=dp),  dimension(:),                   pointer ::  d_f_rb

  real(kind=dp),  dimension(:,:),   allocatable,  target  ::  rho_hat
  real(kind=dp),  dimension(:,:),                 pointer ::  d_f_rh
  logical                                                 ::  h_present

  integer,        parameter                     ::  h_num = 1, ni_num = 28  ! faster than lookup in consts


  ! Harmonic potential:
  real(kind=dp)                 ::  harm_period = 100.0_dp                            ! isotropic...
  real(kind=dp),  dimension(3)  ::  harm_spring_const = (/ 1.0_dp, 1.0_dp, 1.0_dp /)  ! can be anisotropic..

  public  ::  pot_get_potential
  public  ::  pot_get_forces

contains

subroutine pot_get_potential(cell, energy)
  use io, only: io_err
  implicit none
  type(cell_type),  intent(in)  ::  cell
  real(kind=dp),    intent(out) ::  energy

  select case (potential_type)
    case ('eam')
      if (pot_initialized .eqv. .false.) call pot_eam_init(cell)
      call pot_eam_potential(cell, energy)
    case default
      call io_err("pot_get_potential: unsupported potential type")
  end select
end subroutine pot_get_potential

subroutine pot_get_forces(cell, force, potential)
  use constants,  only: two_pi, units_natural_to_atomic
  use io,         only: io_err
  implicit none
  type(cell_type),                          intent(in)    ::  cell
  real(kind=dp),  dimension(3,cell%natoms), intent(out)   ::  force
  real(kind=dp),                            intent(out)   ::  potential
  ! local for SHO
  real(kind=dp) ::  omega

  select case (potential_type)
    case ('eam')
      if (pot_initialized .eqv. .false.) call pot_eam_init(cell)
      call pot_eam_forces(cell, force, potential)
    case default
      call io_err("Error in pot_get_force: unsupported potential type")
  end select
end subroutine pot_get_forces

subroutine pot_eam_init(cell)
  use io,   only: io_err
  implicit none
  type(cell_type),  intent(in)  ::  cell
  integer ::  iatom, istat, n_hydrogen

  h_present = .false.
  n_hydrogen = 0

  ! first check that all atoms are either hydrogen or nickel:
  do iatom = 1, cell%natoms
    if ((cell%atom_species(iatom) .ne. h_num) .and. (cell%atom_species(iatom) .ne. ni_num)) then
      call io_err("Error in potential - pot_eam_init: Unexpected atom type.")
    else
      if (cell%atom_species(iatom) .eq. h_num) n_hydrogen = n_hydrogen + 1
    end if
  end do

  if (n_hydrogen .gt. 0) h_present = .true.

  allocate(rij_vec(3, cell%natoms, cell%natoms), stat=istat)!arraysize), stat=istat)
  if (istat .ne. 0) call io_err("Error in pot_eam_init - failed to allocate rij_vec array.")

  allocate(rij_mag(cell%natoms, cell%natoms), stat=istat)!arraysize), stat=istat)
  if (istat .ne. 0) call io_err("Error in pot_eam_init - failed to allocate rij_mag array.")

  allocate(rho_bar(cell%natoms), stat=istat)
  if (istat .ne. 0) call io_err("Error in pot_eam_init - failed to allocate rho_bar array.")

  allocate(d_rho_atomic(cell%natoms, cell%natoms), stat=istat)
  if (istat .ne. 0) call io_err("Error in pot_eam_init - failed to allocate d_rho_atomic array.")

  ! initialize pointers here:
  ! once we are done with data on RHS, we will store LHS data in these arrays instead
  ! pointers save memory and make code easier to follow
  d_phi => rij_mag
  d_f_rb => rho_bar

  allocate(rho_hat(2, max(1,n_hydrogen)), stat=istat)
  if (istat .ne. 0) call io_err("pot_eam_init - could not allocate rho_hat array.")
  d_f_rh => rho_hat

  pot_initialized = .true.
end subroutine pot_eam_init

subroutine pot_eam_potential(cell, potential)
  use constants,  only: units_atomic_to_natural, units_natural_to_atomic
  use cell,       only: cell_min_img
  !use algor,      only: algor_symmetric_matrix_index
  implicit none
  type(cell_type),  intent(in)  ::  cell
  real(kind=dp),    intent(out) ::  potential
  ! local variables
  real(kind=dp) ::  rho_a
  integer ::  iatom, jatom, ihydrogen!, ielement, i

  potential = 0.0_dp  ! we have sums

  ! rij vectors and pairwise terms:
  do iatom = 1, cell%natoms!-1
    do jatom = 1, cell%natoms!iatom+1, cell%natoms
      if (iatom .ne. jatom) then
        !ielement = algor_symmetric_matrix_index(iatom, jatom, cell%natoms)
        ! could instead set ielement = 1 outside loop then do ielement = ielement+1 inside inner loop
        ! followed by the same thing in the outer loop (after the inner loop has completed) to match this.
        ! speed difference is expected to be small, so just call this function to keep things simple.

        ! rij(4, ielement) is magintude of vector
        call cell_min_img(cell, iatom, jatom, rij_vec(:,iatom,jatom), distance=rij_mag(iatom,jatom))!ielement), distance=rij(4,ielement))
        ! make unit vectors but keep original magnitude
        rij_vec(:,iatom,jatom) = rij_vec(:,iatom,jatom)/rij_mag(iatom,jatom)!ielement) = rij(1:3,ielement)/rij(4,ielement)

        ! scale magnitude: length from Bohr to Angstroms
        rij_mag(iatom,jatom) = units_atomic_to_natural(length=rij_mag(iatom,jatom))

        ! Calculate the pairwise potential energy contributions:
        ! Since we are only working on half of the matrix (excluding the diagonal), we don't need to x0.5
        ! Note: no need to store each of these terms as they only contribute to the potential energy
        potential = potential + 0.5_dp*eam_phi(rij_mag(iatom,jatom), cell%atom_species(iatom), cell%atom_species(jatom))
        !eam_phi(rij(4,ielement), cell%atom_species(iatom), cell%atom_species(jatom))
      end if
    end do  !jatom
  end do  !iatom


  ! Get density dependent terms
  rho_bar(:) = 0.0_dp

  ! extra terms for H
  rho_hat(:,:) = 0.0_dp
  ihydrogen = 1

  do iatom = 1, cell%natoms

    do jatom = 1, cell%natoms
      ! i .ne. j
      if (iatom .ne. jatom) then
        ! get index of reduced N^2 array
!        ielement = algor_symmetric_matrix_index(iatom, jatom, cell%natoms)

        rho_a = eam_rho(rij_mag(iatom,jatom), cell%atom_species(jatom))!ielement), cell%atom_species(jatom))
        rho_bar(iatom) = rho_bar(iatom) + rho_a

        if (h_present .and. (cell%atom_species(iatom) .eq. h_num)) then
          rho_hat(1,ihydrogen) = rho_hat(1,ihydrogen) + rho_a**2        ! to power Q1=2
          rho_hat(2,ihydrogen) = rho_hat(2,ihydrogen) + sqrt(rho_a)     ! to power Q2=0.5
        end if

      end if
    end do  !jatom

    ! rho_hat(:,1) = 0.0_dp if hydrogen doesn't exist.. won't affect eam_f_ni anyway
    potential = potential + eam_f(rho_bar(iatom), rho_hat(1,ihydrogen), rho_hat(2,ihydrogen), cell%atom_species(iatom))
    if (h_present .and. (cell%atom_species(iatom) .eq. h_num)) ihydrogen = ihydrogen + 1
  end do !iatom

  ! rescale units into atomic
  potential = units_natural_to_atomic(energy=potential)
end subroutine pot_eam_potential

subroutine pot_eam_forces(cell, force, potential)
  use constants,  only: units_natural_to_atomic
!  use algor,      only: algor_symmetric_matrix_index
  implicit none
  type(cell_type),                          intent(in)  ::  cell
  real(kind=dp),  dimension(3,cell%natoms), intent(out) ::  force
  real(kind=dp),                            intent(out) ::  potential
  ! local variables
  real(kind=dp) :: contrib
  integer ::  iatom, jatom, ihydrogen, jhydrogen, ihat

  ! need to fill out some of the matrices anyway, so potential is basically free...
  call pot_eam_potential(cell, potential)

  ! get everything we need for the forces in one go...
  do iatom = 1, cell%natoms
    ! have rho_bar, need dF/rho_bar         0 is rho_bar derivative
    d_f_rb(iatom) = eam_d_f(rho_bar(iatom), 0, cell%atom_species(iatom))

    do jatom = 1, cell%natoms
      if (iatom .ne. jatom) then
        ! MUST do d_rho_atomic first since d_phi points to the rij_mag array!
        d_rho_atomic(iatom,jatom) = eam_d_rho(rij_mag(iatom,jatom), cell%atom_species(jatom))
        d_phi(iatom, jatom) = eam_d_phi(rij_mag(iatom,jatom), cell%atom_species(iatom), cell%atom_species(jatom))
      end if
    end do
  end do

  if (h_present) then
    do ihydrogen = 1, size(rho_hat,2)
      do ihat = 1, 2
        d_f_rh(ihat,ihydrogen) = eam_d_f_h(rho_hat(ihat,ihydrogen), ihat)   ! always use eam_d_f_h
      end do
    end do
  end if

  ! Can actually calculate the forces now
  force(:,:) = 0.0_dp

  ihydrogen = 0
  do iatom = 1, cell%natoms

    if (h_present .and. (cell%atom_species(iatom) .eq. h_num)) ihydrogen = ihydrogen + 1
    jhydrogen = 0

    do jatom = 1, cell%natoms
      if (iatom .ne. jatom) then

        contrib = d_f_rb(jatom)*d_rho_atomic(jatom,iatom) + d_f_rb(iatom)*d_rho_atomic(iatom,jatom) + d_phi(jatom,iatom)

        if (h_present) then
          if (cell%atom_species(iatom) .eq. h_num) then
            if(d_rho_atomic(iatom,jatom) > 0.0_dp) &
            contrib = contrib + d_f_rh(1,ihydrogen)*2.0_dp*d_rho_atomic(iatom,jatom) &
                            & + d_f_rh(2,ihydrogen)*0.5_dp/sqrt(d_rho_atomic(iatom,jatom))
          end if

          if (cell%atom_species(jatom) .eq. h_num) then
            jhydrogen = jhydrogen + 1
            if(d_rho_atomic(jatom,iatom) > 0.0_dp) &
            contrib = contrib + d_f_rh(1,jhydrogen)*2.0_dp*d_rho_atomic(jatom,iatom) &
                            & + d_f_rh(2,jhydrogen)*0.5_dp/sqrt(d_rho_atomic(jatom,iatom))
          end if
        end if

        force(:,iatom) = force(:,iatom) - contrib*rij_vec(:, iatom, jatom)
      end if
    end do
  end do

  ! rescale units:
  force(:,:) = units_natural_to_atomic(force=force(:,:))
end subroutine pot_eam_forces

function eam_phi(rij, atomi, atomj)
  implicit none
  real(kind=dp),  intent(in)  ::  rij
  integer,        intent(in)  ::  atomi, atomj
  real(kind=dp)               ::  eam_phi
  ! associated constants:
  real(kind=dp),  parameter   ::  c = 14.3888_dp  ! eV Angstrom

  eam_phi = c*eam_z(rij, atomi)*eam_z(rij, atomj)/rij
end function eam_phi

function eam_d_phi(rij, atomi, atomj)
  implicit none
  real(kind=dp),  intent(in)  ::  rij
  integer,        intent(in)  ::  atomi, atomj
  real(kind=dp)               ::  eam_d_phi
  ! associated constants:
  real(kind=dp),  parameter   ::  c = 14.3888_dp  ! eV Angstrom

  eam_d_phi = c*((eam_d_z(rij,atomi)*eam_z(rij,atomj) + eam_d_z(rij,atomj)*eam_z(rij,atomi))/rij &
  &           - (eam_z(rij,atomi)*eam_z(rij,atomj))/(rij**2))
end function eam_d_phi

function eam_z(rij, atomtype)
  implicit none
  real(kind=dp),  intent(in)  ::  rij
  integer,        intent(in)  ::  atomtype
  real(kind=dp)               ::  eam_z

  ! do no error checking for atom type - we do this on initialization
  select case (atomtype)
    case (h_num)
      ! hydrogen
      eam_z = eam_z_h(rij)*eam_s(rij)    ! multiply by smoothing function
    case (ni_num)
      ! nickel
      eam_z = eam_z_ni(rij)*eam_s(rij)   ! multiply by smoothing function
  end select
end function eam_z

function eam_s(rij)
  use constants,  only: pi, two_pi
  implicit none
  real(kind=dp),  intent(in)  ::  rij
  real(kind=dp)               ::  eam_s
  ! parameters:
  real(kind=dp),  parameter   ::  r_c = 5.0_dp    ! Angstrom
  real(kind=dp),  parameter   ::  delta = 5.0_dp  ! Angstrom
  real(kind=dp),  parameter   ::  r_c_p_delta = r_c+delta

  if (rij .le. r_c) then
    eam_s = 1.0_dp
  else if ((rij .gt. r_c) .and. (rij .lt. r_c_p_delta)) then
    eam_s = (rij-r_c_p_delta)/(-delta) - sin(pi*((2.0_dp*rij - 2.0_dp*r_c - delta)/delta))/two_pi
    if (eam_s .lt. 0.0_dp) eam_s = 0.0_dp   ! fix for when the above goes negative (specific values of 9.99...)
  else
    ! rij .ge. r_c+delta
    eam_s = 0.0_dp
  end if
end function eam_s

function eam_d_s(rij)
  use constants,  only: pi
  implicit none
  real(kind=dp),  intent(in)  ::  rij
  real(kind=dp)               ::  eam_d_s
  ! parameters:
  real(kind=dp),  parameter   ::  r_c = 5.0_dp    ! Angstrom
  real(kind=dp),  parameter   ::  delta = 5.0_dp  ! Angstrom
  real(kind=dp),  parameter   ::  r_c_p_delta = r_c+delta

  if (rij .le. r_c) then
    eam_d_s = 0.0_dp
  else if ((rij .gt. r_c) .and. (rij .lt. r_c_p_delta)) then
    eam_d_s = -(1.0_dp + cos(pi*(2.0_dp*rij - 2.0_dp*r_c - delta)/delta))/delta
  else
    ! rij .ge. r_c+delta
    eam_d_s = 0.0_dp
  end if
end function eam_d_s

function eam_d_z(rij, atomtype)
  implicit none
  real(kind=dp),  intent(in)  ::  rij
  integer,        intent(in)  ::  atomtype
  real(kind=dp)               ::  eam_d_z

  ! no error checking, assume correct atom types
  select case (atomtype)
    case (h_num)
      eam_d_z = eam_d_z_h(rij)*eam_s(rij) + eam_z_h(rij)*eam_d_s(rij)
    case (ni_num)
      eam_d_z = eam_d_z_ni(rij)*eam_s(rij) + eam_z_ni(rij)*eam_d_s(rij)
  end select
end function eam_d_z

function eam_z_ni(rij)
  implicit none
  real(kind=dp),  intent(in)  ::  rij
  real(kind=dp)               ::  eam_z_ni
  ! associated constants:
  real(kind=dp),  parameter   ::  z0 = 10.0_dp  ! Dimensionless
  real(kind=dp),  parameter   ::  a = 0.537_dp  ! Angstrom
  real(kind=dp),  parameter   ::  b = 1.116_dp  ! Angstrom
! real(kind=dp),  parameter   ::  c = 1.0_dp    ! Dimensionless

  ! same as EAM5, but written as 1/a, 1/b (same constants)
! eam_z_ni = z0*(1.0_dp + (rij/b)**c)*algor_exp(-rij/a)
  eam_z_ni = z0*(1.0_dp + (rij/b))*algor_exp(-rij/a)
end function eam_z_ni

function eam_d_z_ni(rij)
  implicit none
  real(kind=dp),  intent(in)  ::  rij
  real(kind=dp)               ::  eam_d_z_ni
  ! associated constants:
  real(kind=dp),  parameter   ::  z0 = 10.0_dp  ! Dimensionless
  real(kind=dp),  parameter   ::  a = 0.537_dp  ! Angstrom
  real(kind=dp),  parameter   ::  b = 1.116_dp  ! Angstrom
! real(kind=dp),  parameter   ::  c = 1.0_dp    ! Dimensionless

  ! same as EAM5, but written as 1/a, 1/b (same constants)

  !eam_d_z_ni = z0*algor_exp(-rij/a)*( ((rij/b)**c)*( (c/rij) - (1.0_dp/a) ) - (1.0_dp/a) )
  eam_d_z_ni = z0*algor_exp(-rij/a)*( (-1.0_dp/a) + (1.0_dp/b) - (rij/(a*b)) )
end function eam_d_z_ni

function eam_z_h(rij)
  implicit none
  real(kind=dp),  intent(in)  ::  rij
  real(kind=dp)               ::  eam_z_h
  integer ::  j
  ! associated constants:
  real(kind=dp),  dimension(3), parameter   ::  a = (/ 0.374_dp, 0.010_dp, 0.011_dp /)    ! Angstrom
  real(kind=dp),  dimension(3), parameter   ::  b = (/ 0.932_dp, 0.549_dp, 0.377_dp /)    ! Angstrom
  real(kind=dp),  dimension(3), parameter   ::  c = (/ 2.775_dp, 142.678_dp, 95.455_dp /) ! Dimensionless
  real(kind=dp),  dimension(3), parameter   ::  d = (/ 0.00752_dp, 8.0036_dp, 2.401_dp /) ! Dimensionless, Angstrom^-1, Angstrom

  eam_z_h = 0.0_dp
  do j = 1, 3
    eam_z_h = eam_z_h + ((rij/b(j))**c(j))*algor_exp(-rij/a(j))
  end do
  eam_z_h = eam_z_h + d(1)*sech(d(2)*(rij-d(3)))
end function eam_z_h

function eam_d_z_h(rij)
  implicit none
  real(kind=dp),  intent(in)  ::  rij
  real(kind=dp)               ::  eam_d_z_h
  real(kind=dp)               ::  argument
  integer ::  j
  ! associated constants:
  real(kind=dp),  dimension(3), parameter   ::  a = (/ 0.374_dp, 0.010_dp, 0.011_dp /)    ! Angstrom
  real(kind=dp),  dimension(3), parameter   ::  b = (/ 0.932_dp, 0.549_dp, 0.377_dp /)    ! Angstrom
  real(kind=dp),  dimension(3), parameter   ::  c = (/ 2.775_dp, 142.678_dp, 95.455_dp /) ! Dimensionless
  real(kind=dp),  dimension(3), parameter   ::  d = (/ 0.00752_dp, 8.0036_dp, 2.401_dp /) ! Dimensionless, Angstrom^-1, Angstrom


  eam_d_z_h = 0.0_dp
  do j = 1, 3
    eam_d_z_h = eam_d_z_h + algor_exp(-rij/a(j))*((rij/b(j))**c(j))*(c(j)/rij - (1.0_dp/a(j)))
  end do
  argument = d(2)*(rij-d(3))
  eam_d_z_h = eam_d_z_h - d(1)*d(2)*tanh(argument)*sech(argument)
end function eam_d_z_h

function eam_f(rho_bar, rho_hat_1, rho_hat_2, atomtype)
  implicit none
  real(kind=dp) ::  eam_f
  real(kind=dp),  intent(in)  ::  rho_bar
  real(kind=dp),  intent(in)  ::  rho_hat_1
  real(kind=dp),  intent(in)  ::  rho_hat_2
  integer,        intent(in)  ::  atomtype

  select case (atomtype)
    case (h_num)
      ! hydrogen
      eam_f = eam_f_h(rho_bar, rho_hat_1, rho_hat_2)
    case (ni_num)
      ! nickel
      eam_f = eam_f_ni(rho_bar)
  end select
end function eam_f

function eam_d_f(rho_bar_or_hat, derivative, atomtype)
  implicit none
  real(kind=dp),    intent(in)  ::  rho_bar_or_hat
  integer,          intent(in)  ::  derivative
  integer,          intent(in)  ::  atomtype
  real(kind=dp)                 ::  eam_d_f

  select case (atomtype)
    case (h_num)
      eam_d_f = eam_d_f_h(rho_bar_or_hat, derivative)
    case (ni_num)
      eam_d_f = eam_d_f_ni(rho_bar_or_hat)
  end select
end function eam_d_f

function eam_f_ni(rho_bar)
  implicit none
  real(kind=dp),  intent(in)  ::  rho_bar
  real(kind=dp)               ::  eam_f_ni
  ! associated constants:
  real(kind=dp), parameter    ::  a_ni = -303.289_dp    ! eV Angstrom^3
  real(kind=dp), parameter    ::  alpha_ni = 7.647_dp   ! Angstrom^3
  real(kind=dp), parameter    ::  a_s = -2.015E10_dp    ! eV Angstrom^15
  real(kind=dp), parameter    ::  rho_c = 0.11_dp       ! Angstrom^-3
  real(kind=dp), parameter    ::  b_ni = 87.987_dp      ! eV^1/3 Angstrom^3
  real(kind=dp), parameter    ::  beta_ni = 75.075_dp   ! Angstrom^3
  real(kind=dp), parameter    ::  b_s = -5.384E8_dp     ! eV Angstrom^12
  real(kind=dp), parameter    ::  delta = 0.01_dp       ! Angstrom^-3
  real(kind=dp), parameter    ::  c_ni = -522.774_dp    ! eV Angstrom^3
  real(kind=dp), parameter    ::  gamma_ni = 373.379_dp ! Angstrom^3
  real(kind=dp), parameter    ::  c_s = -4.060E6_dp     ! eV Angstrom^9
  real(kind=dp), parameter    ::  d_ni = 39.421_dp      ! eV^1/5 Angstrom^3
  real(kind=dp), parameter    ::  delta_ni = 56.342_dp  ! Angstrom^3
  real(kind=dp), parameter    ::  d_s = -11.031_dp      ! eV
  real(kind=dp), parameter    ::  rho_c_m_delta = rho_c-delta
  real(kind=dp) ::  rho_bar_m_rho_c

  if ((rho_bar .ge. 0.0_dp) .and. (rho_bar .lt. rho_c_m_delta)) then
    ! EAM6 paper wrong here:
    eam_f_ni = a_ni*rho_bar*algor_exp(-alpha_ni*rho_bar) + ((b_ni*rho_bar)**3)*algor_exp(-beta_ni*rho_bar) &
    &         + c_ni*rho_bar*algor_exp(-gamma_ni*rho_bar) + ((d_ni*rho_bar)**5)*algor_exp(-delta_ni*rho_bar)
  else if ((rho_bar .ge. rho_c_m_delta) .and. (rho_bar .le. rho_c)) then
    rho_bar_m_rho_c = rho_bar - rho_c
    eam_f_ni = a_s*(rho_bar_m_rho_c)**5 + b_s*(rho_bar_m_rho_c)**4 + c_s*(rho_bar_m_rho_c)**3 + d_s
  else
    ! ie: rho_bar > rho_c
    eam_f_ni = d_s
  end if
end function eam_f_ni

function eam_d_f_ni(rho_bar)
  implicit none
  real(kind=dp),  intent(in)  ::  rho_bar
  real(kind=dp)               ::  eam_d_f_ni
  ! associated constants:
  real(kind=dp), parameter    ::  a_ni = -303.289_dp    ! eV Angstrom^3
  real(kind=dp), parameter    ::  alpha_ni = 7.647_dp   ! Angstrom^3
  real(kind=dp), parameter    ::  a_s = -2.015E10_dp    ! eV Angstrom^15
  real(kind=dp), parameter    ::  rho_c = 0.11_dp       ! Angstrom^-3
  real(kind=dp), parameter    ::  b_ni = 87.987_dp      ! eV^1/3 Angstrom^3
  real(kind=dp), parameter    ::  beta_ni = 75.075_dp   ! Angstrom^3
  real(kind=dp), parameter    ::  b_s = -5.384E8_dp     ! eV Angstrom^12
  real(kind=dp), parameter    ::  delta = 0.01_dp       ! Angstrom^-3
  real(kind=dp), parameter    ::  c_ni = -522.774_dp    ! eV Angstrom^3
  real(kind=dp), parameter    ::  gamma_ni = 373.379_dp ! Angstrom^3
  real(kind=dp), parameter    ::  c_s = -4.060E6_dp     ! eV Angstrom^9
  real(kind=dp), parameter    ::  d_ni = 39.421_dp      ! eV^1/5 Angstrom^3
  real(kind=dp), parameter    ::  delta_ni = 56.342_dp  ! Angstrom^3
  real(kind=dp), parameter    ::  rho_c_m_delta = rho_c-delta
  real(kind=dp) ::  rho_bar_m_rho_c

  if ((rho_bar .ge. 0.0_dp) .and. (rho_bar .lt. rho_c_m_delta)) then
    ! EAM6 paper wrong here:
    eam_d_f_ni = a_ni*algor_exp(-alpha_ni*rho_bar)*(1.0_dp - alpha_ni*rho_bar) &
    & + (b_ni**3)*algor_exp(-beta_ni*rho_bar)*(3.0_dp*rho_bar**2 - beta_ni*rho_bar**3) &
    & + c_ni*algor_exp(-gamma_ni*rho_bar)*(1.0_dp - gamma_ni*rho_bar)  &
    & + (d_ni**5)*algor_exp(-delta_ni*rho_bar)*(5.0_dp*rho_bar**4 - delta_ni*rho_bar**5)
  else if ((rho_bar .ge. rho_c_m_delta) .and. (rho_bar .le. rho_c)) then
    rho_bar_m_rho_c = rho_bar - rho_c
    eam_d_f_ni = 5.0_dp*a_s*(rho_bar_m_rho_c)**4 + 4.0_dp*b_s*(rho_bar_m_rho_c)**3 + 3.0_dp*c_s*(rho_bar_m_rho_c)**2
  else
    ! ie: rho_bar > rho_c
    eam_d_f_ni = 0.0_dp
  end if
end function eam_d_f_ni

function eam_f_h(rho_bar, rho_hat_1, rho_hat_2)
  implicit none
  real(kind=dp),  intent(in)  ::  rho_bar
  real(kind=dp),  intent(in)  ::  rho_hat_1
  real(kind=dp),  intent(in)  ::  rho_hat_2
  real(kind=dp)               ::  eam_f_h
  ! associated constants:
  real(kind=dp),  dimension(2),   parameter ::  e_h = (/ 476.121_dp, -543.394_dp /)   ! eV Angstrom^3
  real(kind=dp),  dimension(2),   parameter ::  epsilon_h = (/ 5.072_dp, 5.285_dp /)  ! Angstrom^3
  ! Fortran is column major (row, column)
  real(kind=dp),  dimension(3,2), parameter ::  &
  ! xi:       11 (eV),   12 (Ang^6), 13 (Ang-6),   21 (eV),   22 (Ang^3/2),  23 (Ang^-3/2)
  & xi = reshape((/ -0.04084_dp, 3.686E8_dp, 13563.0_dp, -0.00325_dp,   54.772_dp,     0.812_dp /), shape(xi))
  real(kind=dp),  dimension(2) :: rho_hat
  integer :: j, k

  eam_f_h = 0.0_dp
  do j = 1,2
    eam_f_h = eam_f_h + e_h(j)*rho_bar*algor_exp(epsilon_h(j)*rho_bar)
  end do

  if (non_local_correction) then
    rho_hat(:) = (/ rho_hat_1, rho_hat_2 /)
    do k = 1,2
      eam_f_h = eam_f_h + xi(1,k)*sech(xi(2,k)*(rho_hat(k)-xi(3,k))) 
    end do
  end if
end function eam_f_h

function eam_d_f_h(rho_bar_or_hat, derivative)
  use io, only: io_err
  implicit none
  real(kind=dp),    intent(in)  ::  rho_bar_or_hat
  integer,          intent(in)  ::  derivative
  real(kind=dp)                 ::  eam_d_f_h
  ! associated constants:
  real(kind=dp),  dimension(2),   parameter ::  e_h = (/ 476.121_dp, -543.394_dp /)   ! eV Angstrom^3
  real(kind=dp),  dimension(2),   parameter ::  epsilon_h = (/ 5.072_dp, 5.285_dp /)  ! Angstrom^3
  ! Fortran is column major (row, column)
  real(kind=dp),  dimension(3,2), parameter ::  &
  ! xi:       11 (eV),   12 (Ang^6), 13 (Ang-6),   21 (eV),   22 (Ang^3/2),  23 (Ang^-3/2)
  & xi = reshape((/ -0.04084_dp, 3.686E8_dp, 13563.0_dp, -0.00325_dp,   54.772_dp,     0.812_dp /), shape(xi))
  integer :: j, k
  real(kind=dp) :: argument

  k = 0
  select case (derivative)
    case (0)    ! rho_bar
      ! primary term:
      eam_d_f_h = 0.0_dp
      do j = 1,2
        eam_d_f_h = eam_d_f_h + e_h(j)*algor_exp(epsilon_h(j)*rho_bar_or_hat)*(1.0_dp + rho_bar_or_hat*epsilon_h(j))
      end do
    case (1,2)  ! rho_hat_1/rho_hat_2
      k = derivative
    case default
      call io_err("Error in potential - eam_d_f_h: Unexpected derivative")
  end select
  
  ! non-local term:
  if (k .ne. 0) then
    if (non_local_correction) then
      ! k is either 1, 2 or 0
      argument = xi(2,k)*(rho_bar_or_hat-xi(3,k))
      eam_d_f_h = -xi(1,k)*xi(2,k)*tanh(argument)*sech(argument)
    else
      eam_d_f_h = 0.0_dp
    end if
  end if
end function eam_d_f_h

function eam_rho(rij, atomtype)
  implicit none
  real(kind=dp)               ::  eam_rho
  real(kind=dp),  intent(in)  ::  rij
  integer,        intent(in)  ::  atomtype

  ! no error checking for atom type - we do this on initialization
  select case (atomtype)
    case (h_num)
      ! Hydrogen
      eam_rho = eam_rho_h(rij)*eam_s(rij)
    case (ni_num)
      ! Nickel
      eam_rho = eam_rho_ni(rij)*eam_s(rij)
  end select
end function eam_rho

function eam_d_rho(rij, atomtype)
  implicit none
  real(kind=dp)               ::  eam_d_rho
  real(kind=dp),  intent(in)  ::  rij
  integer,        intent(in)  ::  atomtype

  ! no error checking for atom type - we do this on initialization
  select case (atomtype)
    case (h_num)
      ! Hydrogen
      eam_d_rho = eam_d_rho_h(rij)*eam_s(rij) + eam_rho_h(rij)*eam_d_s(rij)
    case (ni_num)
      ! Nickel
      eam_d_rho = eam_d_rho_ni(rij)*eam_s(rij) + eam_rho_ni(rij)*eam_d_s(rij)
  end select
end function eam_d_rho

function eam_rho_h(rij)
  use constants,  only: pi, units_atomic_to_natural
  implicit none
  real(kind=dp),  intent(in)  ::  rij
  real(kind=dp)               ::  eam_rho_h
  real(kind=dp)               ::  bohr_radius

  bohr_radius = units_atomic_to_natural(length=1.0_dp)
  eam_rho_h = algor_exp(-2.0_dp*rij/bohr_radius)/(pi*bohr_radius**3)
end function eam_rho_h

function eam_d_rho_h(rij)
  use constants,  only: pi, units_atomic_to_natural 
  implicit none
  real(kind=dp),  intent(in)  ::  rij
  real(kind=dp)               ::  eam_d_rho_h
  real(kind=dp)               ::  bohr_radius

  bohr_radius = units_atomic_to_natural(length=1.0_dp)
  eam_d_rho_h = -2.0_dp*algor_exp(-2.0_dp*rij/bohr_radius)/(pi*bohr_radius**4)
end function eam_d_rho_h

function eam_rho_ni(rij)
  implicit none
  real(kind=dp)               ::  eam_rho_ni
  real(kind=dp),  intent(in)  ::  rij
  real(kind=dp),  parameter   ::  n_s = 2.0_dp
  real(kind=dp),  parameter   ::  n = 10.0_dp
  real(kind=dp),  parameter   ::  n_m_ns = n-n_s

  eam_rho_ni = n_s*eam_rho_shell_ni(rij,'s') + n_m_ns*eam_rho_shell_ni(rij,'d')
end function eam_rho_ni

function eam_d_rho_ni(rij)
  implicit none
  real(kind=dp)               ::  eam_d_rho_ni
  real(kind=dp),  intent(in)  ::  rij
  real(kind=dp),  parameter   ::  n_s = 2.0_dp
  real(kind=dp),  parameter   ::  n = 10.0_dp
  real(kind=dp),  parameter   ::  n_m_ns = n-n_s

  eam_d_rho_ni = n_s*eam_d_rho_shell_ni(rij,'s') + n_m_ns*eam_d_rho_shell_ni(rij,'d')
end function eam_d_rho_ni

function eam_rho_shell_ni(rij, shelltype)
  use constants,  only: pi
  use io,         only: io_err
  implicit none
  real(kind=dp)                 ::  eam_rho_shell_ni
  real(kind=dp),    intent(in)  ::  rij
  character(len=1), intent(in)  ::  shelltype
  integer ::  i
  ! parameters:
  real(kind=dp),  parameter   ::  over_4pi = 1.0_dp/(4.0_dp*pi)
  ! s shell
  real(kind=dp),  dimension(8), parameter ::  c_s = (/ -0.00389_dp, -0.02991_dp, -0.03189_dp, 0.15289_dp, -0.20048_dp, &
  & -0.05423_dp, 0.49292_dp, 0.61875_dp /)
  ! d shell
  real(kind=dp),  dimension(2), parameter ::  c_d = (/ 0.42120_dp, 0.70658_dp /)

  eam_rho_shell_ni = 0.0_dp

  select case (shelltype)
    case ('s')
      do i = 1,8
        eam_rho_shell_ni = eam_rho_shell_ni + c_s(i)*eam_r_shell_ni(rij,i,'s')
      end do
    case ('d')
      do i = 1,2
        eam_rho_shell_ni = eam_rho_shell_ni + c_d(i)*eam_r_shell_ni(rij,i,'d')
      end do
    case default
      call io_err("Error in potential: eam_rho_shell_ni - unexpected shell type")
  end select

  eam_rho_shell_ni = over_4pi*eam_rho_shell_ni**2
end function eam_rho_shell_ni

function eam_d_rho_shell_ni(rij, shelltype)
  use constants,  only: two_pi
  use io,         only: io_err
  implicit none
  real(kind=dp)                 ::  eam_d_rho_shell_ni
  real(kind=dp),    intent(in)  ::  rij
  character(len=1), intent(in)  ::  shelltype
  integer ::  i
  real(kind=dp) :: ciri
  ! parameters:
  real(kind=dp),  parameter   ::  over_2pi = 1.0_dp/two_pi
  ! s shell
  real(kind=dp),  dimension(8), parameter ::  c_s = (/ -0.00389_dp, -0.02991_dp, -0.03189_dp, 0.15289_dp, -0.20048_dp, &
  & -0.05423_dp, 0.49292_dp, 0.61875_dp /)
  ! d shell
  real(kind=dp),  dimension(2), parameter ::  c_d = (/ 0.42120_dp, 0.70658_dp /)

  eam_d_rho_shell_ni = 0.0_dp
  ciri = 0.0_dp

  select case (shelltype)
    case ('s')
      do i = 1,8
        ciri = ciri + c_s(i)*eam_r_shell_ni(rij,i,'s')
        eam_d_rho_shell_ni = eam_d_rho_shell_ni + c_s(i)*eam_d_r_shell_ni(rij,i,'s')
      end do
    case ('d')
      do i = 1,2
        ciri = ciri + c_d(i)*eam_r_shell_ni(rij,i,'d')
        eam_d_rho_shell_ni = eam_d_rho_shell_ni + c_d(i)*eam_d_r_shell_ni(rij,i,'d')
      end do
    case default
      call io_err("Error in potential: eam_d_rho_shell_ni - unexpected shell type")
  end select

  ! might need abs here
  eam_d_rho_shell_ni = over_2pi*eam_d_rho_shell_ni*ciri
end function eam_d_rho_shell_ni

function eam_r_shell_ni(rij, i, shelltype)
  use io,     only: io_err
  implicit none
  real(kind=dp)                 ::  eam_r_shell_ni
  real(kind=dp),    intent(in)  ::  rij
  integer,          intent(in)  ::  i
  character(len=1), intent(in)  ::  shelltype
  integer ::  n_i
  ! parameters:                               zeta is in Ang^-1
  ! s shell:
  real(kind=dp),  dimension(8), parameter ::  zeta_s = (/ 54.88885_dp, 38.48431_dp, 27.42703_dp, 20.88204_dp, 10.95707_dp, &
  & 7.31958_dp, 3.92650_dp, 2.15289_dp /)
  ! d shell:
  real(kind=dp),  dimension(2), parameter ::  zeta_d = (/ 12.67582_dp, 5.43253_dp /)
  real(kind=dp),  dimension(4), parameter ::  sqrt_factorial_2ni = (/ sqrt(2.0_dp),   sqrt(24.0_dp), &
                                                                  &   sqrt(720.0_dp), sqrt(40320.0_dp) /)

  select case (shelltype)
    case ('s')
      ! assume i is between 1 and 8 (inclusive) - don't do error checking here
      ! make use of integer division
      n_i = (i+1)/2

      eam_r_shell_ni = ((2.0_dp*zeta_s(i))**(real(n_i,dp)+0.5_dp))*(rij**(n_i-1))*(algor_exp(-zeta_s(i)*rij)) &
      & / sqrt_factorial_2ni(n_i) 
    case ('d')
      ! again, don't check i - assume everything is okay
      !n_i = 3 ! always
      ! same functional form as above, except we know some values
      eam_r_shell_ni = ((2.0_dp*zeta_d(i))**3.5_dp)*(rij**2)*(algor_exp(-zeta_d(i)*rij))/sqrt_factorial_2ni(3)
    case default
      call io_err("Error in potential: eam_r_shell_ni - unexpected shell type")
  end select
end function eam_r_shell_ni

function eam_d_r_shell_ni(rij, i, shelltype)
  use algor,  only: algor_factorial
  use io,     only: io_err
  implicit none
  real(kind=dp)                 ::  eam_d_r_shell_ni
  real(kind=dp),    intent(in)  ::  rij
  integer,          intent(in)  ::  i
  character(len=1), intent(in)  ::  shelltype
  integer ::  n_i
  ! parameters:                               zeta is in Ang^-1
  ! s shell:
  real(kind=dp),  dimension(8), parameter ::  zeta_s = (/ 54.88885_dp, 38.48431_dp, 27.42703_dp, 20.88204_dp, 10.95707_dp, &
  & 7.31958_dp, 3.92650_dp, 2.15289_dp /)
  ! d shell:
  real(kind=dp),  dimension(2), parameter ::  zeta_d = (/ 12.67582_dp, 5.43253_dp /)
  real(kind=dp),  dimension(4), parameter ::  sqrt_factorial_2ni = (/ sqrt(2.0_dp),   sqrt(24.0_dp), &
                                                                  &   sqrt(720.0_dp), sqrt(40320.0_dp) /)
  select case (shelltype)
    case ('s')
      ! assume i is between 1 and 8 (inclusive) - don't do error checking here
      ! make use of integer division
      n_i = (i+1)/2

      eam_d_r_shell_ni = ((2.0_dp*zeta_s(i))**(real(n_i,dp)+0.5_dp))* &
      & (algor_exp(-zeta_s(i)*rij)*(real((n_i-1),dp)*(rij**(n_i-2)) - zeta_s(i)*(rij**(n_i-1)))) &
      & / sqrt_factorial_2ni(n_i)
    case ('d')
      ! again, don't check i - assume everything is okay
      !n_i = 3 ! always
      ! same functional form as above, except we know some values
      eam_d_r_shell_ni = ((2.0_dp*zeta_d(i))**3.5_dp)* &
      & (algor_exp(-zeta_d(i)*rij)*(2.0_dp*rij - zeta_d(i)*(rij**2)))/sqrt_factorial_2ni(3)
    case default
      call io_err("Error in potential: eam_d_r_shell_ni - unexpected shell type")
  end select
end function eam_d_r_shell_ni

function sech(x)
  implicit none
  real(kind=dp),  intent(in)  ::  x
  real(kind=dp)               ::  sech
  if(x > -500) then
    sech = 1.0_dp/cosh(x)
  else
    sech = 0
  end if
end function sech
end module potential
