module eam5
  ! EAM 5 potential for H-Ni
  use constants,  only: dp
  use cell,       only: cell_type
  use algor,      only: algor_exp
  implicit none

  private

  logical,            save      ::  pot_initialized = .false.

  ! EAM potential variables:
  ! Uses module level variables - not suitable for multiple cells
  ! (with different natoms) without reinitializing
  real(kind=dp),  dimension(:,:,:), allocatable,  save  ::  rij_vec
  real(kind=dp),  dimension(:,:),   allocatable,  save  ::  d_rho_atomic
  real(kind=dp),  dimension(:,:),   allocatable,  save  ::  d_phi
  real(kind=dp),  dimension(:),     allocatable,  save  ::  d_f_rb
  integer,        parameter   ::  h_num = 1, ni_num = 28

  public  ::  eam_init
  public  ::  eam_get_potential
  public  ::  eam_get_forces
  public  ::  eam_finalize

contains

subroutine eam_get_potential(cell, potential)
  ! wrapper routine - initializes if necessary
  implicit none
  type(cell_type),  intent(in)  ::  cell
  real(kind=dp),    intent(out) ::  potential

  if (.not. pot_initialized) call eam_init(cell)
  call eam_potential(cell, potential)
end subroutine eam_get_potential

subroutine eam_get_forces(cell, force, potential)
  ! wrapper routine - initializes if necessary
  implicit none
  type(cell_type),                          intent(in)    ::  cell
  real(kind=dp),  dimension(3,cell%natoms), intent(out)   ::  force
  real(kind=dp),                            intent(out)   ::  potential

  ! don't check size of forces array - assume correct size
  if (.not. pot_initialized) call eam_init(cell)
  call eam_forces(cell, force, potential)
end subroutine eam_get_forces

subroutine eam_init(cell)
  use io,   only: io_err
  implicit none
  type(cell_type),  intent(in)  ::  cell
  integer ::  iatom, istat

  ! first check that all atoms are either hydrogen or nickel:
  do iatom = 1, cell%natoms
    if ((cell%atom_species(iatom) .ne. h_num) .and. (cell%atom_species(iatom) .ne. ni_num)) &
      & call io_err("eam_init: Unexpected atom type")
  end do

  ! assume nothing previously allocated:
  allocate(rij_vec(3, cell%natoms, cell%natoms), stat=istat)
  if (istat .ne. 0) call io_err("eam_init: Could not allocate rij_vec array")

  allocate(d_rho_atomic(cell%natoms, cell%natoms), stat=istat)
  if (istat .ne. 0) call io_err("eam_init: Could not allocate d_rho_atomic array")

  allocate(d_phi(cell%natoms, cell%natoms), stat=istat)
  if (istat .ne. 0) call io_err("eam_init: Could not allocate d_phi array")

  allocate(d_f_rb(cell%natoms), stat=istat)
  if (istat .ne. 0) call io_err("eam_init: Could not allocate d_f_rb array")

  pot_initialized = .true.
end subroutine eam_init

subroutine eam_finalize
  use io,   only: io_err
  implicit none
  integer ::  istat

  ! deallocate everything:
  ! if pot_initialized, then these *should* be allocated, but check anyway..
  if (allocated(rij_vec)) then
    deallocate(rij_vec, stat=istat)
    if (istat .ne. 0) call io_err("eam_finalize: Could not deallocate rij_vec array")
  end if

  if (allocated(d_rho_atomic)) then
    deallocate(d_rho_atomic, stat=istat)
    if (istat .ne. 0) call io_err("eam_finalize: Could not deallocate d_rho_atomic array")
  end if

  if (allocated(d_phi)) then
    deallocate(d_phi, stat=istat)
    if (istat .ne. 0) call io_err("eam_finalize: Could not deallocate d_phi array")
  end if

  if (allocated(d_f_rb)) then
    deallocate(d_f_rb, stat=istat)
    if (istat .ne. 0) call io_err("eam_finalize: Could not deallocate d_f_rb array")
  end if

  pot_initialized = .false.
end subroutine eam_finalize

!------------------------------------------
! private routines below here
!------------------------------------------

subroutine eam_potential(cell, potential)
  use constants,  only: units_atomic_to_natural, units_natural_to_atomic
  use cell,       only: cell_min_img
  implicit none
  type(cell_type),  intent(in)  ::  cell
  real(kind=dp),    intent(out) ::  potential
  ! local variables
  real(kind=dp) ::  rho_bar, rij_mag
  integer       ::  iatom, jatom, ispecies, jspecies

  ! possibly inefficient method - we don't take into account symmetric matrices etc - might be doing unnecessary work
  potential = 0.0_dp

  !$omp parallel do default(none), schedule(static), &
  !$omp & private(iatom, jatom, rho_bar, ispecies, jspecies, rij_mag), shared(rij_vec, cell), reduction(+:potential)
  do iatom = 1, cell%natoms
    rho_bar = 0.0_dp
    ispecies = cell%atom_species(iatom)
    do jatom = 1, cell%natoms
      if (iatom .ne. jatom) then
        jspecies = cell%atom_species(jatom)
        call cell_min_img(cell, iatom, jatom, rij_vec(:,iatom,jatom), distance=rij_mag)

        ! make unit vectors but keep original magnitude
        rij_vec(:,iatom,jatom) = rij_vec(:,iatom,jatom)/rij_mag

        ! scale magnitude: length from Bohr to Angstroms
        rij_mag = units_atomic_to_natural(length=rij_mag)

        ! densities:
        rho_bar = rho_bar + eam_rho(rij_mag, jspecies)

        ! pairwise potential energy contributions:
        potential = potential + 0.5_dp*eam_phi(rij_mag, ispecies, jspecies)
      end if
    end do !jatom

    ! density dependent potential energy contribution:
    potential = potential + eam_f(rho_bar, ispecies)
  end do !iatom
  !$omp end parallel do

  ! rescale units into atomic
  potential = units_natural_to_atomic(energy=potential)
end subroutine eam_potential

subroutine eam_forces(cell, force, potential)
  use constants,  only: units_atomic_to_natural, units_natural_to_atomic
  use cell,       only: cell_min_img
  implicit none
  type(cell_type),                          intent(in)  ::  cell
  real(kind=dp),  dimension(3,cell%natoms), intent(out) ::  force
  real(kind=dp),                            intent(out) ::  potential
  ! local variables
  real(kind=dp) :: rho_bar, contrib, rij_mag
  integer ::  iatom, jatom, ispecies, jspecies

  ! Merge potential and force calculations for efficiency:
  ! (inefficient method - we don't take into account symmetric matrices etc - might be doing unnecessary work)
  potential = 0.0_dp
  force(:,:) = 0.0_dp

  !$omp parallel default(none), &
  !$omp & private(iatom, jatom, rho_bar, ispecies, jspecies, contrib, rij_mag), &
  !$omp & shared(rij_vec, cell, d_f_rb, d_rho_atomic, d_phi, force, potential)
  !$omp do schedule(static), reduction(+:potential)
  do iatom = 1, cell%natoms
    rho_bar = 0.0_dp
    ispecies = cell%atom_species(iatom)
    do jatom = 1, cell%natoms
      if (iatom .ne. jatom) then
        jspecies = cell%atom_species(jatom)
        call cell_min_img(cell, iatom, jatom, rij_vec(:,iatom,jatom), distance=rij_mag)

        ! make unit vectors but keep original magnitude
        rij_vec(:,iatom,jatom) = rij_vec(:,iatom,jatom)/rij_mag

        ! scale magnitude: length from Bohr to Angstroms
        rij_mag = units_atomic_to_natural(length=rij_mag)

        rho_bar = rho_bar + eam_rho(rij_mag, jspecies)
        d_rho_atomic(iatom,jatom) = eam_d_rho(rij_mag, jspecies)
        d_phi(iatom, jatom) = eam_d_phi(rij_mag, ispecies, jspecies)

        ! pairwise potential energy contributions:
        potential = potential + 0.5_dp*eam_phi(rij_mag, ispecies, jspecies)
      end if
    end do !jatom

    ! have rho_bar, need dF/rho_bar
    d_f_rb(iatom) = eam_d_f(rho_bar, ispecies)

    ! density dependent potential energy contribution:
    potential = potential + eam_f(rho_bar, ispecies)
  end do !iatom
  !$omp end do


  !$omp do schedule(static)
  do iatom = 1, cell%natoms
    do jatom = 1, cell%natoms
      if (iatom .ne. jatom) then
        contrib = d_f_rb(jatom)*d_rho_atomic(jatom,iatom) + d_f_rb(iatom)*d_rho_atomic(iatom,jatom) + d_phi(jatom,iatom)

        force(:,iatom) = force(:,iatom) - contrib*rij_vec(:, iatom, jatom)
      end if
    end do
  end do
  !$omp end do
  !$omp end parallel

  ! rescale units
  potential = units_natural_to_atomic(energy=potential)
  force(:,:) = units_natural_to_atomic(force=force(:,:))
end subroutine eam_forces

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
  real(kind=dp),  parameter   ::  a = 1.8633_dp ! Angstrom^-1
  real(kind=dp),  parameter   ::  b = 0.8957_dp ! Angstrom^-1
  ! c = 1.0
  eam_z_ni = z0*(1.0_dp + b*rij)*algor_exp(-a*rij)
end function eam_z_ni

function eam_d_z_ni(rij)
  implicit none
  real(kind=dp),  intent(in)  ::  rij
  real(kind=dp)               ::  eam_d_z_ni
  ! associated constants:
  real(kind=dp),  parameter   ::  z0 = 10.0_dp  ! Dimensionless
  real(kind=dp),  parameter   ::  a = 1.8633_dp ! Angstrom^-1
  real(kind=dp),  parameter   ::  b = 0.8957_dp ! Angstrom^-1
  ! c = 1.0

  eam_d_z_ni = z0*algor_exp(-a*rij)*( b-a-a*b*rij )
end function eam_d_z_ni

function eam_z_h(rij)
  implicit none
  real(kind=dp),  intent(in)  ::  rij
  real(kind=dp)               ::  eam_z_h
  ! associated constants:
  real(kind=dp),  parameter   ::  z0 = 0.1959_dp  ! Dimensionless
  real(kind=dp),  parameter   ::  a = 1.7957_dp ! Angstrom^-1
  real(kind=dp),  parameter   ::  b = 3.2108_dp ! Angstrom^-1
  ! c = 1.0

  eam_z_h = z0*(1.0_dp + b*rij)*algor_exp(-a*rij)
end function eam_z_h

function eam_d_z_h(rij)
  implicit none
  real(kind=dp),  intent(in)  ::  rij
  real(kind=dp)               ::  eam_d_z_h
  ! associated constants:
  real(kind=dp),  parameter   ::  z0 = 0.1959_dp  ! Dimensionless
  real(kind=dp),  parameter   ::  a = 1.7957_dp ! Angstrom^-1
  real(kind=dp),  parameter   ::  b = 3.2108_dp ! Angstrom^-1
  ! c = 1.0

  eam_d_z_h = z0*algor_exp(-a*rij)*( b-a-a*b*rij )
end function eam_d_z_h

function eam_f(rho_bar, atomtype)
  implicit none
  real(kind=dp) ::  eam_f
  real(kind=dp),  intent(in)  ::  rho_bar
  integer,        intent(in)  ::  atomtype

  select case (atomtype)
    case (h_num)
      ! hydrogen
      eam_f = eam_f_h(rho_bar)
    case (ni_num)
      ! nickel
      eam_f = eam_f_ni(rho_bar)
  end select
end function eam_f

function eam_d_f(rho_bar, atomtype)
  implicit none
  real(kind=dp),    intent(in)  ::  rho_bar
  integer,          intent(in)  ::  atomtype
  real(kind=dp)                 ::  eam_d_f

  select case (atomtype)
    case (h_num)
      ! hydrogen
      eam_d_f = eam_d_f_h(rho_bar)
    case (ni_num)
      ! nickel
      eam_d_f = eam_d_f_ni(rho_bar)
  end select
end function eam_d_f

function eam_f_ni(rho_bar)
  implicit none
  real(kind=dp),  intent(in)  ::  rho_bar
  real(kind=dp)               ::  eam_f_ni
  ! associated constants:
  real(kind=dp), parameter    ::  a = -126.5009308_dp       ! eV Angstrom^3
  real(kind=dp), parameter    ::  alpha = 0.3362141252_dp   ! Angstrom^3
  real(kind=dp), parameter    ::  a_s = -1.657208422E10_dp  ! eV Angstrom^15
  real(kind=dp), parameter    ::  rho_c = 0.21_dp           ! Angstrom^-3
  real(kind=dp), parameter    ::  b = 6033.287109_dp        ! eV Angstrom^9
  real(kind=dp), parameter    ::  beta = 11.02211666_dp     ! Angstrom^3
  real(kind=dp), parameter    ::  b_s = -5.13226816E8_dp    ! eV Angstrom^12
  real(kind=dp), parameter    ::  delta = 0.01_dp           ! Angstrom^-3
  real(kind=dp), parameter    ::  c = -209.7682800_dp       ! eV Angstrom^3
  real(kind=dp), parameter    ::  gamma_ = 51.76818085_dp   ! Angstrom^3
  real(kind=dp), parameter    ::  c_s = -4.377938E6_dp      ! eV Angstrom^9
  real(kind=dp), parameter    ::  d_s = -19.23463058_dp     ! eV
  real(kind=dp), parameter    ::  rho_c_m_delta = rho_c-delta
  real(kind=dp) ::  rho_bar_m_rho_c

  if ((rho_bar .ge. 0.0_dp) .and. (rho_bar .le. rho_c_m_delta)) then
    eam_f_ni = a*rho_bar*algor_exp(-alpha*rho_bar) + b*rho_bar**3*algor_exp(-beta*rho_bar) &
    &         + c*rho_bar*algor_exp(-gamma_*rho_bar)
  else if ((rho_bar .gt. rho_c_m_delta) .and. (rho_bar .le. rho_c)) then
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
  real(kind=dp), parameter    ::  a = -126.5009308_dp       ! eV Angstrom^3
  real(kind=dp), parameter    ::  alpha = 0.3362141252_dp   ! Angstrom^3
  real(kind=dp), parameter    ::  a_s = -1.657208422E10_dp  ! eV Angstrom^15
  real(kind=dp), parameter    ::  rho_c = 0.21_dp           ! Angstrom^-3
  real(kind=dp), parameter    ::  b = 6033.287109_dp        ! eV Angstrom^9
  real(kind=dp), parameter    ::  beta = 11.02211666_dp     ! Angstrom^3
  real(kind=dp), parameter    ::  b_s = -5.13226816E8_dp    ! eV Angstrom^12
  real(kind=dp), parameter    ::  delta = 0.01_dp           ! Angstrom^-3
  real(kind=dp), parameter    ::  c = -209.7682800_dp       ! eV Angstrom^3
  real(kind=dp), parameter    ::  gamma_ = 51.76818085_dp   ! Angstrom^3
  real(kind=dp), parameter    ::  c_s = -4.377938E6_dp      ! eV Angstrom^9
  real(kind=dp), parameter    ::  d_s = -19.23463058_dp     ! eV
  real(kind=dp), parameter    ::  rho_c_m_delta = rho_c-delta
  real(kind=dp) ::  rho_bar_m_rho_c

  if ((rho_bar .ge. 0.0_dp) .and. (rho_bar .le. rho_c_m_delta)) then
    eam_d_f_ni = a*algor_exp(-alpha*rho_bar)*(1.0_dp - alpha*rho_bar) &
    & + b*rho_bar**2*algor_exp(-beta*rho_bar)*(3.0_dp - beta*rho_bar) &
    & + c*algor_exp(-gamma_*rho_bar)*(1.0_dp - gamma_*rho_bar)
  else if ((rho_bar .gt. rho_c_m_delta) .and. (rho_bar .le. rho_c)) then
    rho_bar_m_rho_c = rho_bar - rho_c
    eam_d_f_ni = 5.0_dp*a_s*(rho_bar_m_rho_c)**4 + 4.0_dp*b_s*(rho_bar_m_rho_c)**3 + 3.0_dp*c_s*(rho_bar_m_rho_c)**2
  else
    ! ie: rho_bar > rho_c
    eam_d_f_ni = 0.0_dp
  end if
end function eam_d_f_ni

function eam_f_h(rho_bar)
  implicit none
  real(kind=dp),  intent(in)  ::  rho_bar
  real(kind=dp)               ::  eam_f_h
  ! associated constants:
  real(kind=dp),  parameter   ::  alpha_h = -70.5461_dp     ! eV*Ang^3
  real(kind=dp),  parameter   ::  beta_h  = 6.9507_dp       ! Ang^3

  eam_f_h = alpha_h*rho_bar*algor_exp(-beta_h*rho_bar)
end function eam_f_h

function eam_d_f_h(rho_bar)
  implicit none
  real(kind=dp),    intent(in)  ::  rho_bar
  real(kind=dp)                 ::  eam_d_f_h
  ! associated constants:
  real(kind=dp),  parameter   ::  alpha_h = -70.5461_dp     ! eV*Ang^3
  real(kind=dp),  parameter   ::  beta_h  = 6.9507_dp       ! Ang^3

  eam_d_f_h = alpha_h*algor_exp(-beta_h*rho_bar)*(1.0_dp - beta_h*rho_bar)
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
      call io_err("eam_rho_shell_ni: Unexpected shell type")
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
      call io_err("eam_d_rho_shell_ni: Unexpected shell type")
  end select

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
      ! n_i = 3 ! always
      ! same functional form as above, except we know some values
      eam_r_shell_ni = ((2.0_dp*zeta_d(i))**3.5_dp)*(rij**2)*(algor_exp(-zeta_d(i)*rij))/sqrt_factorial_2ni(3)
    case default
      call io_err("eam_r_shell_ni: Unexpected shell type")
  end select
end function eam_r_shell_ni

function eam_d_r_shell_ni(rij, i, shelltype)
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
      ! n_i = 3 ! always
      ! same functional form as above, except we know some values
      eam_d_r_shell_ni = ((2.0_dp*zeta_d(i))**3.5_dp)* &
      & (algor_exp(-zeta_d(i)*rij)*(2.0_dp*rij - zeta_d(i)*(rij**2)))/sqrt_factorial_2ni(3)
    case default
      call io_err("eam_d_r_shell_ni : Unexpected shell type")
  end select
end function eam_d_r_shell_ni
end module eam5
