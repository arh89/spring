module geometry

  use constants,  only: dp
  use cell,       only: cell_type
  implicit none

  private
  
  ! default geom opt parameters:
  integer,        parameter ::  default_max_iter = 100
  real(kind=dp),  parameter ::  default_tpsd_initial_lambda = 1.0E-05_dp
  real(kind=dp),  parameter ::  default_max_disp = 5.0E-01_dp             ! Ang
  real(kind=dp),  parameter ::  default_force_tol = 1.0E-05_dp            ! eV/Ang
  real(kind=dp),  parameter ::  default_energy_tol = 1.0E-08_dp           ! eV
  real(kind=dp),  parameter ::  default_disp_tol = 1.0E-03_dp             ! Ang
  real(kind=dp),  parameter ::  default_stress_tol = 1.0E-01_dp           ! GPa
  logical,        parameter ::  default_fix_cell = .false.
  logical,        parameter ::  default_fix_atoms = .false.

  ! hardcoded for now - possibly get from input?
  real(kind=dp),  parameter   ::  d_epsilon = 5.0E-07_dp

  type, public  ::  geom_type
    integer                         ::  output_unit_num

    real(kind=dp),  dimension(3,3)              ::  orig_lattice_vectors

    real(kind=dp),  dimension(3,3)              ::  strain_tensor
    real(kind=dp),  dimension(3,3)              ::  stress_tensor
    real(kind=dp),  dimension(:,:), allocatable ::  inv_hessian
    real(kind=dp),  dimension(:,:), allocatable ::  atom_forces
    real(kind=dp),  dimension(:),   allocatable ::  x_vec
    real(kind=dp),  dimension(:),   allocatable ::  f_vec
    real(kind=dp),  dimension(:),   allocatable ::  x_vec_old
    real(kind=dp),  dimension(:),   allocatable ::  f_vec_old
    real(kind=dp),  dimension(:),   allocatable ::  delta_vec
    real(kind=dp)                               ::  lambda
    real(kind=dp)                               ::  pressure
    real(kind=dp)                               ::  enthalpy
    real(kind=dp)                               ::  potential_energy

    ! set from input
    logical                         ::  fix_cell
    logical                         ::  fix_atoms
    integer                         ::  max_iter
    real(kind=dp)                   ::  force_tol
    real(kind=dp)                   ::  energy_tol
    real(kind=dp)                   ::  disp_tol
    real(kind=dp)                   ::  stress_tol
    real(kind=dp)                   ::  max_disp
    real(kind=dp)                   ::  tpsd_initial_lambda
    real(kind=dp),  dimension(3,3)  ::  ext_pressure_tensor
  end type geom_type
  
  public  ::  geom_init
  public  ::  do_geom_opt
  public  ::  geom_deallocate
  public  ::  geom_calc_stress_tensor

contains

subroutine geom_init(geom, cell)
  use io,         only: io_open_file, seedname
  use constants,  only: units_natural_to_atomic
  implicit none
  type(geom_type),  intent(inout) ::  geom
  type(cell_type),  intent(inout) ::  cell

  ! get parameters from input file
  call geom_read_input(geom)

  ! convert parameters into atomic units:
  geom%force_tol = units_natural_to_atomic(force=geom%force_tol)
  geom%energy_tol = units_natural_to_atomic(energy=geom%energy_tol)
  geom%disp_tol = units_natural_to_atomic(length=geom%disp_tol)
  geom%stress_tol = units_natural_to_atomic(pressure=geom%stress_tol)
  geom%max_disp = units_natural_to_atomic(length=geom%max_disp)

  geom%ext_pressure_tensor(:,:) = units_natural_to_atomic(pressure=geom%ext_pressure_tensor(:,:))

  ! allocate everything
  call geom_allocate(geom, cell)

  geom%orig_lattice_vectors(:,:) = cell%lattice_vectors(:,:)

  ! pressure = Tr(pressure_tensor)/3
  geom%pressure = (geom%ext_pressure_tensor(1,1) + geom%ext_pressure_tensor(2,2) + geom%ext_pressure_tensor(3,3))/3.0_dp

  geom%strain_tensor(:,:) = 0.0_dp

  call geom_inv_hessian_init(geom)

  call io_open_file(trim(seedname)//'.geom', geom%output_unit_num, 'replace')
end subroutine geom_init


subroutine geom_allocate(geom, cell)
  use io, only: io_err
  implicit none
  type(geom_type),  intent(inout) ::  geom
  type(cell_type),  intent(in)    ::  cell
  integer ::  array_size, istat

  array_size = 9+3*cell%natoms

  allocate(geom%x_vec(array_size), stat=istat)
  if (istat .ne. 0) call io_err("geom_allocate: Could not allocate x_vec array")

  allocate(geom%x_vec_old(array_size), stat=istat)
  if (istat .ne. 0) call io_err("geom_allocate: Could not allocate x_vec_old array")

  allocate(geom%f_vec(array_size), stat=istat)
  if (istat .ne. 0) call io_err("geom_allocate: Could not allocate f_vec array")

  allocate(geom%f_vec_old(array_size), stat=istat)
  if (istat .ne. 0) call io_err("geom_allocate: Could not allocate f_vec_old array")

  allocate(geom%inv_hessian(array_size, array_size), stat=istat)
  if (istat .ne. 0) call io_err("geom_allocate: Could not allocate inv_hessian array")

  allocate(geom%delta_vec(array_size), stat=istat)
  if (istat .ne. 0) call io_err("geom_allocate: Could not allocate delta_vec array")

  allocate(geom%atom_forces(3, cell%natoms), stat=istat)
  if (istat .ne. 0) call io_err("geom_allocate: Could not allocate atom_forces array")
end subroutine geom_allocate


subroutine geom_deallocate(geom)
  use io, only: io_err
  implicit none
  type(geom_type),  intent(inout) ::  geom
  integer ::  istat

  deallocate(geom%x_vec, stat=istat)
  if (istat .ne. 0) call io_err("geom_deallocate: Could not deallocate x_vec array")

  deallocate(geom%x_vec_old, stat=istat)
  if (istat .ne. 0) call io_err("geom_deallocate: Could not deallocate x_vec_old array")

  deallocate(geom%f_vec, stat=istat)
  if (istat .ne. 0) call io_err("geom_deallocate: Could not deallocate f_vec array")

  deallocate(geom%f_vec_old, stat=istat)
  if (istat .ne. 0) call io_err("geom_deallocate: Could not deallocate f_vec_old array")

  deallocate(geom%inv_hessian, stat=istat)
  if (istat .ne. 0) call io_err("geom_deallocate: Could not deallocate inv_hessian array")

  deallocate(geom%delta_vec, stat=istat)
  if (istat .ne. 0) call io_err("geom_deallocate: Could not deallocate delta_vec array")

  deallocate(geom%atom_forces, stat=istat)
  if (istat .ne. 0) call io_err("geom_deallocate: Could not deallocate atom_forces array")
end subroutine geom_deallocate


subroutine geom_inv_hessian_init(geom)
  implicit none
  type(geom_type),  intent(inout) ::  geom
  integer ::  i

  ! for now, let's just initialize to the identity.. can eventually use Pfrommer method
  geom%inv_hessian(:,:) = 0.0_dp

  do i = 1, size(geom%inv_hessian, 1)
    geom%inv_hessian(i,i) = 1.0_dp
  end do
end subroutine geom_inv_hessian_init

subroutine geom_read_input(geom)
  use io, only: io_input_get_single_value, max_line_len, io_input_get_data, io_str_get_num_tokens, io_err
  implicit none
  type(geom_type),  intent(inout) ::  geom

  character(len=max_line_len) ::  keyword
  character(len=max_line_len),  allocatable,  dimension(:)  ::  input_data
  logical ::  kw_found
  integer ::  iline, istat

  ! single value parameters:
  call io_input_get_single_value('geom_max_iter', geom%max_iter, default_max_iter)
  if (geom%max_iter .le. 0) call io_err("geom_read_input: geom_max_iter must be > 0")

  call io_input_get_single_value('geom_tpsd_initial_lambda', geom%tpsd_initial_lambda, default_tpsd_initial_lambda)
  if (geom%tpsd_initial_lambda .le. 0.0_dp) call io_err("geom_read_input: geom_tpsd_initial_lambda must be > 0.0")

  call io_input_get_single_value('geom_force_tol', geom%force_tol, default_force_tol)
  if (geom%force_tol .le. 0.0_dp) call io_err("geom_read_input: geom_force_tol must be > 0.0")

  call io_input_get_single_value('geom_energy_tol', geom%energy_tol, default_energy_tol)
  if (geom%energy_tol .le. 0.0_dp) call io_err("geom_read_input: geom_energy_tol must be > 0.0")

  call io_input_get_single_value('geom_disp_tol', geom%disp_tol, default_disp_tol)
  if (geom%disp_tol .le. 0.0_dp) call io_err("geom_read_input: geom_disp_tol must be > 0.0")

  call io_input_get_single_value('geom_stress_tol', geom%stress_tol, default_stress_tol)
  if (geom%stress_tol .le. 0.0_dp) call io_err("geom_read_input: geom_stress_tol must be > 0.0")
  
  call io_input_get_single_value('geom_max_disp', geom%max_disp,  default_max_disp)
  if (geom%max_disp .le. 0.0_dp) call io_err("geom_read_input: geom_max_disp must be > 0.0")

  call io_input_get_single_value('geom_fix_cell', geom%fix_cell,  default_fix_cell)

  call io_input_get_single_value('geom_fix_atoms', geom%fix_atoms,  default_fix_atoms)

  ! block parameters:
  keyword = 'external_pressure'
  call io_input_get_data(keyword, input_data, kw_found)
  if (.not. kw_found) then
    ! zero external pressure:
    geom%ext_pressure_tensor(:,:) = 0.0_dp
  else
    ! input_data must be allocated
    if (size(input_data, 1) .ne. 3) call io_err("geom_read_input: "//trim(keyword)//" block must contain 3 lines")

    ! don't do anything smart - require whole 3x3 array - don't check that it's valid..
    do iline = 1, 3
      if (io_str_get_num_tokens(input_data(iline)) .ne. 3) &
        & call io_err("geom_read_input: "//trim(keyword)//" block must contain 3 values on each line")

      read(input_data(iline), *, iostat=istat) geom%ext_pressure_tensor(iline, :)      
      if (istat .ne. 0) call io_err("geom_read_input: "//trim(keyword)//": Error reading line")
    end do
  end if

  ! clean up
  if (allocated(input_data)) then
    deallocate(input_data, stat=istat)
    if (istat .ne. 0) call io_err("geom_read_input: Could not deallocate input_data array")
  end if
end subroutine geom_read_input

subroutine do_geom_opt(geom, cell)
  use constants,  only: units_atomic_to_natural, element_symbol
  use cell,       only: cell_min_img_vec
  use potential,  only: pot_get_forces, pot_get_potential
  use io,         only: io_err, stdout
  implicit none
  type(geom_type),  intent(inout) ::  geom
  type(cell_type),  intent(inout) ::  cell

  real(kind=dp),  allocatable,  dimension(:,:)  ::  old_cart_pos
  real(kind=dp),  dimension(3)  ::  disp_vec
  real(kind=dp)   ::  old_enthalpy, delta_enthalpy
  real(kind=dp)   ::  max_force, max_stress_comp, max_disp
  real(kind=dp)   ::  current_force_sq, current_disp
  integer         ::  i, iter


  ! routine does TPSD optimization, fixed/variable cell, but no other constraints..

  ! need array to store old cart positions for calculating displacements if fixed cell
  if (geom%fix_cell) then
    allocate(old_cart_pos(3, cell%natoms), stat=i)
    if (i .ne. 0) call io_err("do_geom_opt: Could not allocate old_cart_pos array")

    old_cart_pos(:,:) = cell%atom_cart_pos(:,:)
  end if

  call geom_build_x_vec(geom, cell)

  if (geom%fix_atoms) then
    call pot_get_potential(cell, geom%potential_energy)
    geom%atom_forces = 0.0_dp
  else
    call pot_get_forces(cell, geom%atom_forces, geom%potential_energy)
  end if

  if (geom%fix_cell) then
    geom%stress_tensor(:,:) = 0.0_dp
  else
    call geom_calc_stress_tensor(geom%stress_tensor, cell)
  end if

  write(stdout,*) "Starting geometry optimisation,"
  write(stdout,*) "  Default energy_tol = ", units_atomic_to_natural(energy=geom%energy_tol)  , "eV"
  write(stdout,*) "  Default force_tol  = ", units_atomic_to_natural(force=geom%force_tol)    , "eV/Ang"
  write(stdout,*) "  Default disp_tol   = ", units_atomic_to_natural(length=geom%disp_tol)    , "Ang"
  if(.not. geom%fix_cell) write(stdout,*) "  Default stress_tol = ", units_atomic_to_natural(pressure=geom%stress_tol), "GPa"

  
  call geom_build_f_vec(geom, cell)
  ! geom%atom_forces are not in cartesian space at this point

  geom%enthalpy = geom%potential_energy + geom%pressure*cell%volume

  geom%lambda = geom%tpsd_initial_lambda

  write(stdout,*)
  write(stdout,*) "  Initial Lambda     =", geom%lambda
  write(stdout,*) "  Initial enthalpy   =", units_atomic_to_natural(energy  =geom%enthalpy), "eV"
  write(stdout,*) "  Initial max_force  =", units_atomic_to_natural(force   =maxval(geom%atom_forces)), "eV/Ang"
if(.not. geom%fix_cell) write(stdout,*)"  Initial max_stress =",units_atomic_to_natural(pressure=maxval(geom%stress_tensor)), "GPa"


  ! just to start the whole thing, make copies of x_vec and f_vec
  geom%x_vec_old = geom%x_vec
  geom%f_vec_old = geom%f_vec
  old_enthalpy = geom%enthalpy

  iter = 0
  call geom_write_output(geom, cell, iter)
  do !i = 1, geom%max_iter

    iter = iter + 1
    write(stdout,*)
    write(stdout,*) "Beginning TPSD iteration", iter

    ! X_(i+1) = X_i + lambda*H_i*F_i 
    geom%x_vec = geom%x_vec + geom%lambda*matmul(geom%inv_hessian, geom%f_vec)

    ! with new positions, can now update everything else
    call geom_x_vec_to_cell(geom, cell)
    ! geom_x_vec_to_cell also updates geom%strain_tensor..
    
    ! should strain tensor be reinitialized to zero at every step?

    if (geom%fix_atoms) then
      call pot_get_potential(cell, geom%potential_energy)
      geom%atom_forces = 0.0_dp
    else
      call pot_get_forces(cell, geom%atom_forces, geom%potential_energy)
    end if

    ! not sure which "pressure" to use here.. External pressure?
    geom%enthalpy = geom%potential_energy + geom%pressure*cell%volume

    if (geom%fix_cell) then
      geom%stress_tensor(:,:) = 0.0_dp
    else
      call geom_calc_stress_tensor(geom%stress_tensor, cell)
    end if

    ! get max abs(force) before we build new f_vec
    max_force = 0.0_dp
    do i = 1, cell%natoms
      current_force_sq = dot_product(geom%atom_forces(:,i), geom%atom_forces(:,i)) 
      if (current_force_sq .gt. max_force) max_force = current_force_sq
    end do
    max_force = sqrt(max_force)

    ! max abs(stress) component - only calc if variable cell
    if (.not. geom%fix_cell) then
      max_stress_comp = maxval(abs(geom%stress_tensor/cell%volume+geom%ext_pressure_tensor))
    end if

    ! calc max abs(disp).. Always use same cell, but must use cell_min_img
    if (geom%fix_cell) then
      max_disp = 0.0_dp
      do i = 1, cell%natoms
        call cell_min_img_vec(cell, old_cart_pos(:,i), cell%atom_cart_pos(:,i), disp_vec, current_disp)
        if (current_disp .gt. max_disp) max_disp = current_disp
      end do
    end if


    ! enforce max disp? Don't want to move too far in a single step..


    ! calc change in enthalpy/energy:
    delta_enthalpy = old_enthalpy - geom%enthalpy



!!!!!!!!!!! to do
    ! can now write output...
    call geom_write_output(geom, cell, iter)
    write(stdout,*) "  Lambda     =", geom%lambda
    write(stdout,*) "  Enthalpy   =", units_atomic_to_natural(energy  =geom%enthalpy), 'eV'
    if(abs(delta_enthalpy) <= geom%energy_tol) then
      write(stdout,*) "  Y Max_dH     =", units_atomic_to_natural(energy  =delta_enthalpy), 'eV'
    else
      write(stdout,*) "  N Max_dH     =", units_atomic_to_natural(energy  =delta_enthalpy), 'eV'
    end if

    if(max_force <= geom%force_tol) then
      write(stdout,*) "  Y Max_force  =", units_atomic_to_natural(force   =max_force), 'eV/Ang'
    else
      write(stdout,*) "  N Max_force  =", units_atomic_to_natural(force   =max_force), 'eV/Ang'
    end if
    if(geom%fix_cell) then
      if(max_disp <= geom%disp_tol) then
        write(stdout,*) "  Y  Max_disp   =", units_atomic_to_natural(length  =max_disp), "Ang"
      else
        write(stdout,*) "  N  Max_disp   =", units_atomic_to_natural(length  =max_disp), "Ang"
      end if
    else
      if(max_stress_comp <= geom%stress_tol) then
        write(stdout,*) "  Y Max_stress =", units_atomic_to_natural(pressure=max_stress_comp), 'GPa'
      else
        write(stdout,*) "  N Max_stress =", units_atomic_to_natural(pressure=max_stress_comp), 'GPa'
      end if
    end if

!    write(stdout,*) 'Stress:                                   |  Strain'
!    do i=1,3
!      write(stdout,*) geom%stress_tensor(i,:), '|', geom%strain_tensor(i,:)
!    end do
    ! exit conditions
    if( (max_force <= geom%force_tol) .and. (abs(delta_enthalpy) <= geom%energy_tol) &
      .and. (max_disp <= geom%disp_tol) .and. (geom%fix_cell .or. (max_stress_comp <= geom%stress_tol))) then
      write(stdout,*)
      write(stdout,*) "Converged in ", iter, "iterations"
      write(stdout,*)
      write(stdout,*) "====== Final cell ======"
      write(stdout,*) "begin lattice_vectors"
      do i = 1, 3 
        write(stdout,*) units_atomic_to_natural(length=cell%lattice_vectors(:,i))
      end do
      write(stdout,*) "end lattice_vectors"
      write(stdout,*) 
      write(stdout,*) "begin atomic_positions"
      do i = 1, cell%natoms
        write(stdout,*) element_symbol(cell%atom_species(i)), cell%atom_frac_pos(:,i)
      end do
      write(stdout,*) "end atomic_positions"
      exit
    else if (iter .eq. geom%max_iter) then
      write(stdout,*)
      write(stdout,*) 'Failed to converge in', iter, 'iterations'
      exit
    end if
!!!!!!!!!! end to do

    ! below here necessary for next step only:
    call geom_build_f_vec(geom, cell)
    ! geom%atom_forces are not in cartesian space at this point

    ! for TPSD, let's not update geom%inv_hessian - this will need to be changed for other methods (eg: BFGS)
    ! calculate lambda for next time
    call geom_tpsd_calc_step_size(geom)

    ! store things for next time..
    if (geom%fix_cell) old_cart_pos(:,:) = cell%atom_cart_pos(:,:)
    geom%x_vec_old = geom%x_vec
    geom%f_vec_old = geom%f_vec
    old_enthalpy = geom%enthalpy
  end do
end subroutine do_geom_opt


subroutine geom_tpsd_calc_step_size(geom)
  implicit none
  type(geom_type),  intent(inout) ::  geom

  real(kind=dp) ::  dx_dot_dg, dg_dot_dg
  
  ! TPSD:
  ! lambda = <delta_x, delta_g>/<delta_g, delta_g>
  ! (or lambda = <delta_x, delta_x>/<delta_x, delta_g>)

  ! and:
  ! delta_g = -delta_f = -( F_i - F_(i-1) )

  ! in Pfrommer this is: -(H_i*F_i - H_(i-1)*F_(i-1)) because we multiply by H
  ! if H_i = H_(i-1) = H_0     =>     delta_g = -H_0*( F_i - F_(i-1) )
  
  ! use geom%delta_vec as a work array..
  geom%delta_vec = geom%f_vec - geom%f_vec_old
  geom%delta_vec = -matmul(geom%inv_hessian, geom%delta_vec)
  
  dx_dot_dg = dot_product((geom%x_vec - geom%x_vec_old), geom%delta_vec)
  dg_dot_dg = dot_product(geom%delta_vec, geom%delta_vec)

  geom%lambda = dx_dot_dg/dg_dot_dg
end subroutine geom_tpsd_calc_step_size


subroutine geom_build_x_vec(geom, cell)
  implicit none
  type(geom_type),  intent(inout) ::  geom
  type(cell_type),  intent(in)    ::  cell

  integer ::  i, j, k

  do i = 1, 3 ! row
    do j = 1, 3 ! col
      geom%x_vec(3*(i-1)+j) = geom%strain_tensor(i,j)
    end do
  end do

  ! s is the fractional coordinates..
  ! start at element 10..
  k = 9
  do i = 1, cell%natoms
    do j = 1, 3
      k = k + 1
      geom%x_vec(k) = cell%atom_frac_pos(j,i)
    end do
  end do
end subroutine geom_build_x_vec


subroutine geom_x_vec_to_cell(geom, cell)
  implicit none
  type(geom_type),  intent(inout) ::  geom
  type(cell_type),  intent(inout) ::  cell

  integer :: i, j, k

  do i = 1, 3 ! row
    do j = 1, 3 ! col
      geom%strain_tensor(i,j) = geom%x_vec(3*(i-1)+j)
    end do
  end do

  ! s is the fractional coordinates..
  ! start at element 10..
  k = 9
  do i = 1, cell%natoms
    do j = 1, 3
      k = k + 1
      cell%atom_frac_pos(j,i) = geom%x_vec(k)
    end do
  end do

  call geom_apply_strain(geom%strain_tensor, geom%orig_lattice_vectors, cell)
end subroutine geom_x_vec_to_cell


subroutine geom_apply_strain(strain, orig_lattice, cell)
  implicit none
  real(kind=dp),  dimension(3,3), intent(in)    ::  strain
  real(kind=dp),  dimension(3,3), intent(in)    ::  orig_lattice
  type(cell_type),                intent(inout) ::  cell

  real(kind=dp),  dimension(3,3)  ::  identity
  integer ::  i

  ! apply strain to change lattice vectors - resync everything else (cart positions etc)
  identity = 0.0_dp
  do i = 1, 3
    identity(i,i) = 1.0_dp
  end do

  cell%lattice_vectors = matmul((identity + strain), orig_lattice)

  call geom_sync_cell(cell)
end subroutine geom_apply_strain


subroutine geom_sync_cell(cell)
  use cell, only: cell_calc_vol, cell_get_recip, cell_shift_to_unit_cell, cell_frac_to_cart, cell_min_img_init
  implicit none
  type(cell_type),  intent(inout) ::  cell

  ! assuming that we have either changed cell vectors (either directly or applying a strain)
  ! and/or changed fractional positions

  call cell_calc_vol(cell)
  call cell_get_recip(cell)             ! need vol first
  call cell_shift_to_unit_cell(cell)
  call cell_frac_to_cart(cell)
  call cell_min_img_init(cell)          ! cell size/shape may have changed
end subroutine geom_sync_cell


subroutine geom_build_f_vec(geom, cell)
  use constants,  only: two_pi
  use algor,      only: algor_calc_inverse_real
  implicit none

  type(geom_type),  intent(inout)   ::  geom
  type(cell_type),  intent(in)      ::  cell

  real(kind=dp),    dimension(3,3)  ::  metric_tensor
  real(kind=dp),    dimension(3,3)  ::  f_strain
  integer                           ::  i, j, k

  ! strain components of f_vec:
  ! f^(strain) = -(stress + pressure*vol)*(1 + strain^T)^-1

  ! first calculate rightmost term:
  ! f_strain is identity first..
  f_strain = 0.0_dp
  do i = 1, 3
    f_strain(i,i) = 1.0_dp
  end do

  f_strain = f_strain + transpose(geom%strain_tensor)
  call algor_calc_inverse_real(f_strain)

  f_strain = matmul(-(geom%stress_tensor + geom%ext_pressure_tensor*cell%volume), f_strain)

  ! now have strain components.. can put into F_vec
  ! same sort of form as the X_vec...
  do j = 1, 3
    do i = 1, 3
      geom%f_vec(3*(i-1)+j) = f_strain(i,j)
    end do
  end do

  ! force terms:
  ! "obtained by multiplying the forces on the atoms in the lattice coordinates with the metric tensor"

  ! forces are currently in cartesian space -> convert to fractional..
  geom%atom_forces = matmul(cell%recip_lattice/two_pi, geom%atom_forces)

  metric_tensor = matmul(transpose(cell%lattice_vectors), cell%lattice_vectors)

  ! multiply by metric tensor
  geom%atom_forces = matmul(metric_tensor, geom%atom_forces)

  ! start at element 10..
  k = 9
  do j = 1, cell%natoms
    do i = 1, 3
      k = k + 1
      geom%f_vec(k) = geom%atom_forces(i,j)
    end do
  end do
end subroutine geom_build_f_vec


subroutine geom_calc_stress_tensor(stress_tensor, cell)
  use cell,       only: cell_copy, cell_deallocate
  use potential,  only: pot_get_potential
  implicit none
  real(kind=dp),    dimension(3,3), intent(out) ::  stress_tensor
  type(cell_type),                  intent(in)  ::  cell

  type(cell_type)                 ::  tmp_cell
  real(kind=dp),  dimension(3,3)  ::  d_strain
  real(kind=dp)                   ::  energy
  integer ::  i, j, k

  stress_tensor(:,:) = 0.0_dp

  do i = 1, 3
    do j = i, 3
      ! centered difference requires perturbations in each direction:
      do k = -1, 1, 2 ! -1 = back,  +1 = forward
        call cell_copy(cell, tmp_cell)

        ! d_strain is epsilon' tensor
        d_strain = 0.0_dp
        d_strain(i,j) = real(k,dp)*d_epsilon

        call geom_apply_strain(d_strain, tmp_cell%lattice_vectors, tmp_cell)

        ! calc new energy - don't need forces
        call pot_get_potential(tmp_cell, energy)

        ! sum energies..
        stress_tensor(i,j) = stress_tensor(i,j) + real(k,dp)*energy
      end do

      ! use centered difference
      ! dy/dx = (y(x+h) - y(x-h)) / 2*h
      stress_tensor(i,j) = stress_tensor(i,j)/(2.0_dp*d_epsilon*cell%volume)

      if (i .ne. j) stress_tensor(j,i) = stress_tensor(i,j)
    end do
  end do

  call cell_deallocate(tmp_cell)
end subroutine geom_calc_stress_tensor

subroutine geom_write_output(geom, cell, iter)
  use constants,  only: element_symbol
  implicit none
  type(geom_type),  intent(in)  ::  geom
  type(cell_type),  intent(in)  ::  cell
  integer,          intent(in)  ::  iter
  ! local vars
  real(kind=dp),  dimension(3,3)  ::  stress
  integer :: i

  if (iter .eq. 0) then
    write (geom%output_unit_num, 99) 'BEGIN header'
    write (geom%output_unit_num, *)  cell%natoms
    write (geom%output_unit_num, 99) 'END header'
    write (geom%output_unit_num, *)  ' '
  end if

  write (geom%output_unit_num, 1) iter!, dE_str, Fmax_str, dRmax_str, Smax_str
  write (geom%output_unit_num, 2) geom%potential_energy, geom%enthalpy

  ! lattice vectors (lines are a, b, c)
  do i = 1, 3
     write (geom%output_unit_num, 4) cell%lattice_vectors(:,i)
  end do

  if (.not. geom%fix_cell) then
     !cell stresses
     stress(:,:) = geom%stress_tensor(:,:)/cell%volume + geom%ext_pressure_tensor(:,:)
     do i=1,3
        write (geom%output_unit_num,6) stress(i,:)
     end do
  end if

  ! positions
  do i = 1, cell%natoms
    write (geom%output_unit_num, 7) element_symbol(cell%atom_species(i)), i, cell%atom_cart_pos(:,i)
  end do

  ! forces
  do i = 1, cell%natoms
    write (geom%output_unit_num, 9) element_symbol(cell%atom_species(i)), i, geom%atom_forces(:,i)
  end do

  write (geom%output_unit_num, *) ' '

  1   format(21x,i18,34x,4(3x,a),10x,'  <-- c')
  2   format(18x,2(3x,es24.16e3),27x,'  <-- E')
  4   format(18x,3(3x,es24.16e3),    '  <-- h')
  6   format(18x,3(3x,es24.16e3),    '  <-- S')
  7   format(1x,a8,1x,i8,3(3x,es24.16e3),'  <-- R')
  9   format(1x,a8,1x,i8,3(3x,es24.16e3),'  <-- F')
  99  format(1x,a)
end subroutine geom_write_output

end module geometry
