module singlepoint
  use constants,  only: dp
  use cell,       only: cell_type
  implicit none

  private

  ! Input defaults:
  logical,  parameter ::  default_calc_stress = .false.
  logical,  parameter ::  default_calc_forces = .true.

  ! Work variables: (module level - may want to change to type)..
  logical ::  calc_stress, calc_forces
  real(kind=dp),  dimension(:,:), allocatable ::  forces
  real(kind=dp),  dimension(:,:), allocatable ::  stress


  public  ::  singlepoint_init
  public  ::  do_singlepoint
  public  ::  singlepoint_deallocate

contains

subroutine singlepoint_init(cell)
  use io, only: io_err
  implicit none
  type(cell_type),  intent(in)  ::  cell
  integer ::  istat

  call singlepoint_read_input

  ! only allocate arrays if necessary
  if (calc_forces) then
    allocate(forces(3, cell%natoms), stat=istat)
    if (istat .ne. 0) call io_err("singlepoint_init: Could not allocate forces array")
  end if

  if (calc_stress) then
    allocate(stress(3,3), stat=istat)
    if (istat .ne. 0) call io_err("singlepoint_init: Could not allocate stress array")
  end if
end subroutine singlepoint_init

subroutine singlepoint_read_input
  use io, only: io_input_get_single_value
  implicit none

  call io_input_get_single_value('singlepoint_calculate_forces', calc_forces, default_calc_forces)
  call io_input_get_single_value('singlepoint_calculate_stress', calc_stress, default_calc_stress)
end subroutine singlepoint_read_input

subroutine singlepoint_deallocate
  use io, only: io_err
  implicit none
  integer :: istat

  ! could check calc_forces / calc_stress variables.. but this is safer..
  if (allocated(forces)) then
    deallocate(forces, stat=istat)
    if (istat .ne. 0) call io_err("singlepoint_deallocate: Could not deallocate forces array")
  end if

  if (allocated(stress)) then
    deallocate(stress, stat=istat)
    if (istat .ne. 0) call io_err("singlepoint_deallocate: Could not deallocate stress array")
  end if

  ! why not?
  calc_forces = default_calc_forces
  calc_stress = default_calc_stress
end subroutine singlepoint_deallocate

subroutine do_singlepoint(cell)
  use constants,  only: units_atomic_to_natural, element_symbol
  use potential,  only: pot_get_potential, pot_get_forces
  use geometry,   only: geom_calc_stress_tensor
  use io,         only: io_err, stdout
  implicit none

  type(cell_type),  intent(in)  ::  cell
  ! local vars:
  real(kind=dp) ::  energy, max_abs_f_comp, max_abs_s_comp
  integer ::  iatom, icomp, jcomp, istat

  
  if (calc_forces) then
    call pot_get_forces(cell, forces, energy)
    forces(:,:) = units_atomic_to_natural(force=forces(:,:))
    max_abs_f_comp = 0.0_dp
  else
    call pot_get_potential(cell, energy)
  end if
  
  energy = units_atomic_to_natural(energy=energy)

  write(stdout, *, iostat=istat) "Total energy:", energy
  if (istat .ne. 0) call io_err("do_singlepoint: Could not write to output")

  write(stdout, *, iostat=istat) "Energy per atom:", energy/real(cell%natoms,dp)
  if (istat .ne. 0) call io_err("do_singlepoint: Could not write to output")
  


  if (calc_forces) then
    write(stdout, *, iostat=istat)
    if (istat .ne. 0) call io_err("do_singlepoint: Could not write to output")

    write(stdout, *, iostat=istat) "Forces:" 
    if (istat .ne. 0) call io_err("do_singlepoint: Could not write to output")

    do iatom = 1, cell%natoms
      write(stdout, *, iostat=istat) element_symbol(cell%atom_species(iatom)), iatom, (forces(icomp, iatom), icomp = 1,3)
      if (istat .ne. 0) call io_err("do_singlepoint: Could not write to output")

      do icomp = 1, 3
        if (abs(forces(icomp, iatom)) .gt. max_abs_f_comp) max_abs_f_comp = abs(forces(icomp,iatom))
      end do
    end do

    write(stdout, *, iostat=istat) 'Max |F| comp:', max_abs_f_comp
    if (istat .ne. 0) call io_err("do_singlepoint: Could not write to output")
  end if



  if (calc_stress) then
    call geom_calc_stress_tensor(stress, cell)
    stress(:,:) = units_atomic_to_natural(pressure=stress(:,:))
    max_abs_s_comp = 0.0_dp

    write(stdout, *, iostat=istat)
    if (istat .ne. 0) call io_err("do_singlepoint: Could not write to output")

    write(stdout, *, iostat=istat) 'Stress Tensor:'
    if (istat .ne. 0) call io_err("do_singlepoint: Could not write to output")

    do icomp = 1, 3
      write(stdout, *, iostat=istat) stress(icomp,:)
      if (istat .ne. 0) call io_err("do_singlepoint: Could not write to output")

      do jcomp = 1, 3
        if (abs(stress(icomp, jcomp)) .gt. max_abs_s_comp) max_abs_s_comp = abs(stress(icomp,jcomp))
      end do
    end do

    write(stdout, *, iostat=istat) 'Max |S| comp:', max_abs_s_comp 
    if (istat .ne. 0) call io_err("do_singlepoint: Could not write to output")
  end if
end subroutine do_singlepoint

end module singlepoint 
