module potential 
  use constants,  only: dp
  use cell,       only: cell_type
  implicit none

  private

  ! IO default parameters:
  character(len=10),  parameter                   ::  default_potential = 'eam'

  character(len=20),                save, public  ::  potential_type = default_potential
  real(kind=dp),      dimension(3), save          ::  harmonic_spring

  public  ::  pot_init
  public  ::  pot_get_potential
  public  ::  pot_get_forces
  public  ::  pot_finalize

contains

subroutine pot_init(cell)
  use lj,           only: lj_init
  use eam5,         only: eam_init
  use fourier_pes,  only: fourier_pes_init
  use io,           only: io_err
  use constants,    only: units_natural_to_atomic, units_atomic_to_natural
  implicit none
  type(cell_type),  intent(in)  ::  cell

  real(kind=dp) ::  ev_to_hartree, ang_to_bohr, atomic_time_to_fs

  call pot_read_input
  
  ! check we have correct potential, call initialization routines if necessary
  select case (potential_type)
    case ('lj')
      call lj_init(cell)
    case ('eam')
      call eam_init(cell)
    case ('fourier_pes')
      call fourier_pes_init(cell)
    case ('harmonic')
      ! assume cell is cubic.. can check
      ! convert units
      ev_to_hartree = units_natural_to_atomic(energy=1.0_dp)
      ang_to_bohr = units_natural_to_atomic(length=1.0_dp)
      atomic_time_to_fs = units_atomic_to_natural(time=1.0_dp)
      ! eV Ang^-2 --> Hartree Bohr^-2
      harmonic_spring(:) = harmonic_spring(:)*ev_to_hartree/(ang_to_bohr**2)
      ! arh: probably testing..
      if ((cell%natoms .eq. 1) .and. (cell%atom_species(1) .eq. 1)) then
        print *, 'Harmonic potential:'
        print *, 'omega =', sqrt(harmonic_spring(:)/cell%atom_mass(1))/atomic_time_to_fs, &
                & 'sum:', sum(sqrt(harmonic_spring(:)/cell%atom_mass(1))/atomic_time_to_fs)
      end if
    case ('free')
    case default
      call io_err("pot_init: Unsupported potential type")
  end select
end subroutine pot_init

subroutine pot_read_input
  use io,           only: max_line_len, io_input_get_data, io_input_get_single_value, io_str_get_num_tokens, &
                        & io_str_to_real, io_err
  implicit none
  character(len=max_line_len),  dimension(:), allocatable ::  values
  logical ::  kw_found
  integer ::  ntokens, istat

  call io_input_get_single_value('potential', potential_type, default_potential)

  ! check for harmonic_spring_const if necessary:
  if (potential_type .eq. 'harmonic') then
    call io_input_get_data('harmonic_spring_const', values, kw_found)
    if (kw_found) then
      if (size(values,1) .ne. 1) call io_err("pot_init: harmonic_spring_const must be a single line")

      ntokens = io_str_get_num_tokens(values(1)) 
      if (ntokens .eq. 1) then
        harmonic_spring(:) = io_str_to_real(values(1)) 
      else if (ntokens .eq. 3) then
        read(unit=values(1), fmt=*, iostat=istat) harmonic_spring(1), harmonic_spring(2), harmonic_spring(3)
        if (istat .ne. 0) call io_err("pot_read_input: Could not read harmonic_spring_const")
      else
        call io_err("pot_read_input: harmonic_spring_const expects 1 or 3 values")
      end if
    else
      ! set default
      harmonic_spring(:) = 1.0_dp ! eV Ang^-2
    end if
  end if

  if (allocated(values)) then
    deallocate(values, stat=istat)
    if (istat .ne. 0) call io_err("pot init: Could not deallocate values array")
  end if
end subroutine pot_read_input

! wrapper routines for selecting between potentials
subroutine pot_get_potential(cell, potential)
  use lj,           only: lj_get_potential
  use eam5,         only: eam_get_potential
  use fourier_pes,  only: fourier_pes_get_potential
  use io,           only: io_err
  implicit none
  type(cell_type),  intent(in)    ::  cell
  real(kind=dp),    intent(out)   ::  potential
  real(kind=dp),    dimension(3)  ::  lat_xxyyzz, disp
  integer ::  iatom

  select case (potential_type)
    case ('lj')
      call lj_get_potential(cell, potential)
    case ('eam')
      call eam_get_potential(cell, potential)
    case ('fourier_pes')
      call fourier_pes_get_potential(cell, potential)
    case ('harmonic')
      ! assume cell is cubic..
      potential = 0.0_dp
      lat_xxyyzz(:) = (/ cell%lattice_vectors(1,1), cell%lattice_vectors(2,2), cell%lattice_vectors(3,3) /)
      do iatom = 1, cell%natoms
        disp(:) = cell%atom_cart_pos(:,iatom)-0.5_dp*lat_xxyyzz(:)
        potential = potential + 0.5_dp*dot_product(harmonic_spring(:), disp(:)**2)
      end do
    case ('free')
      potential = 0.0_dp
    case default
      call io_err("pot_get_potential: Unsupported potential type")
  end select
end subroutine pot_get_potential

subroutine pot_get_forces(cell, force, potential)
  use lj,           only: lj_get_forces
  use eam5,         only: eam_get_forces
  use fourier_pes,  only: fourier_pes_get_forces
  use io,           only: io_err
  implicit none
  type(cell_type),                          intent(in)    ::  cell
  real(kind=dp),  dimension(3,cell%natoms), intent(out)   ::  force
  real(kind=dp),                            intent(out)   ::  potential
  real(kind=dp),  dimension(3)  ::  lat_xxyyzz, disp
  integer ::  iatom

  select case (potential_type)
    case ('lj')
      call lj_get_forces(cell, force, potential)
    case ('eam')
      call eam_get_forces(cell, force, potential)
    case ('fourier_pes')
      call fourier_pes_get_forces(cell, force, potential)
    case ('harmonic')
      ! assume cell is cubic..
      potential = 0.0_dp
      lat_xxyyzz(:) = (/ cell%lattice_vectors(1,1), cell%lattice_vectors(2,2), cell%lattice_vectors(3,3) /)
      do iatom = 1, cell%natoms
        disp(:) = cell%atom_cart_pos(:,iatom)-0.5_dp*lat_xxyyzz(:)
        potential = potential + 0.5_dp*dot_product(harmonic_spring(:), disp(:)**2)
        force(:,iatom) = -harmonic_spring(:)*disp(:)
      end do
    case ('free')
      force(:,:) = 0.0_dp
      potential = 0.0_dp
    case default
      call io_err("pot_get_forces: Unsupported potential type")
  end select
end subroutine pot_get_forces

subroutine pot_finalize
  use lj,           only: lj_finalize
  use eam5,         only: eam_finalize
  use fourier_pes,  only: fourier_pes_finalize
  use io,           only: io_err
  implicit none

  select case (potential_type)
    case ('lj')
      call lj_finalize
    case ('eam')
      call eam_finalize
    case ('fourier_pes')
      call fourier_pes_finalize
    case ('harmonic')
    case ('free')
    case default
      call io_err("pot_finalize: Unsupported potential type")
  end select
end subroutine pot_finalize
end module potential
