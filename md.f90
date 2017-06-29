module md

  use constants,  only: dp
  use io,         only: max_line_len
  use cell,       only: cell_type
  implicit none

  private
  character(len=3),   parameter   ::  output_file_extension = '.md'

  ! IO defaults:
  integer,            parameter   ::  default_num_timesteps = 1000
  real(kind=dp),      parameter   ::  default_timestep = 1.0_dp
  character(len=3),   parameter   ::  default_ensemble = 'nve'
  character(len=20),  parameter   ::  default_thermostat = 'langevin'
  real(kind=dp),      parameter   ::  default_temperature = 300.0_dp
  real(kind=dp),      parameter   ::  default_langevin_time = 100.0_dp
  character(len=6),   parameter   ::  default_trajectory_units = 'atomic'
  integer,            parameter   ::  default_output_interval = 1
  integer,            parameter   ::  default_checkpoint_interval = 10000

  type, public  ::  md_type
    character(len=max_line_len+10)  ::  output_filename
    character(len=6)        ::  trajectory_units          ! 'atomic' or 'user'

    integer                 ::  timestep
    integer                 ::  ntimesteps
    integer                 ::  ndof
    real(kind=dp)           ::  dt
    real(kind=dp)           ::  kinetic_energy
    real(kind=dp)           ::  potential_energy
    real(kind=dp)           ::  inst_temp
    character(len=3)        ::  ensemble                  ! nve or nvt
    character(len=20)       ::  thermostat                ! langevin only for now
    real(kind=dp)           ::  temperature
    logical                 ::  fix_com
    integer                 ::  output_interval
    integer                 ::  checkpoint_interval

    ! langevin vars:
    real(kind=dp)           ::  langevin_time
    real(kind=dp)           ::  langevin_noise_stddev
    real(kind=dp)           ::  langevin_force_correction

    real(kind=dp),  allocatable,  dimension(:,:)    ::  pos
    real(kind=dp),  allocatable,  dimension(:,:)    ::  momentum
    real(kind=dp),  allocatable,  dimension(:,:)    ::  force
  end type md_type

  logical,  save  ::  temperature_specified ! in input file - initialize NVE w/temp

  public  ::  md_init
  public  ::  md_initial_momenta
  public  ::  md_copy
  public  ::  md_deallocate
  public  ::  do_md
  public  ::  md_sync_coords
  public  ::  md_velocity_verlet_p1
  public  ::  md_velocity_verlet_p2
  public  ::  md_get_forces
  public  ::  md_langevin_force
  public  ::  md_zero_net_force
  public  ::  md_calc_kinetic_energy
  public  ::  md_write_output
  public  ::  md_continuation
  public  ::  md_type_read_checkpoint
  public  ::  md_type_write_checkpoint

contains

subroutine md_init(md, cell)
  use io,         only: seedname, io_err
  use constants,  only: units_natural_to_atomic
  use checkpoint, only: continuation
  implicit none
  type(md_type),    intent(inout) ::  md
  type(cell_type),  intent(inout) ::  cell
  ! local vars
  integer       ::  iatom

  if (continuation) call io_err("md_init: md_init called rather than md_continuation")

  ! number of degrees of freedom
  ! ...first assume no constraints
  ! (this *must* be done before calling md_set_defaults)
  md%ndof = 3*cell%natoms

  call md_set_defaults(md)

  ! overwrite values by user - check for errors
  call md_read_input(md)

  ! now subtract 3 if linear momentum of COM is fixed (zero)
  if (md%fix_com) md%ndof = md%ndof - 3
  ! sanity check
  if (md%ndof .lt. 0) call io_err("md_init: Negative number of degrees for freedom")

  ! convert to atomic units
  md%dt = units_natural_to_atomic(time=md%dt)

  ! rescale temperature - always safe because we assign a default even with NVE...
  md%temperature = units_natural_to_atomic(temperature=md%temperature)
  md%langevin_time = units_natural_to_atomic(time=md%langevin_time)
  md%langevin_noise_stddev = sqrt((2.0_dp*md%temperature)/(md%langevin_time*md%dt)) ! sqrt(2kT/(dt*langevin_time))
  md%langevin_force_correction = 1.0_dp/(1.0_dp + 0.5_dp*md%dt/md%langevin_time)   ! 1/(1 + gamma*dt/2)

  call md_allocate(md, cell)

  ! copy positions
  do iatom = 1, cell%natoms
    md%pos(:,iatom) = cell%atom_cart_pos(:,iatom)
  end do

  ! intial momenta from Gaussian dist
  call md_initial_momenta(md, cell)

  ! initial forces and potential energy
  call md_get_forces(md, cell)

  ! doesn't matter if Langevin force called here or not...
  ! shouldn't need to zero net force, but do so anyway..
  if (md%fix_com) then
    call md_zero_net_force(md, cell)
  end if

  ! kinetic energy & instantaneous temperature
  call md_calc_kinetic_energy(md, cell)
  if (md%ndof .ne. 0) then
    md%inst_temp = 2.0_dp*md%kinetic_energy/real(md%ndof,dp)
  else
    md%inst_temp = 0.0_dp   ! always, because KE = 0, but avoids NaN
  end if

  ! create trajectory filename
  md%output_filename = trim(seedname)//trim(output_file_extension)

  ! write simulation params
  call md_echo_params(md)
end subroutine md_init

subroutine md_continuation(md, cell, close_checkpoint)
  use algor,      only: rand_seed
  use cell,       only: cell_read_checkpoint
  use checkpoint, only: checkpoint_open, checkpoint_unit, checkpoint_close, continuation
  use constants,  only: units_natural_to_atomic, units_atomic_to_natural
  use io,         only: stdout, io_err
  implicit none
  type(md_type),      intent(inout) ::  md
  type(cell_type),    intent(inout) ::  cell
  logical,  optional, intent(in)    ::  close_checkpoint
  ! local vars:
  type(md_type) ::  tmp_md
  logical ::  int_close
  integer ::  i

  if (.not. continuation) call io_err("md_continuation: md_continuation called rather than md_init")
  ! should we close checkpoint file after use? default yes..
  if (present(close_checkpoint) .and. (.not. close_checkpoint)) then
    int_close = .false.
  else
    int_close = .true.
  end if

  ! number of degrees of freedom
  ! ...first assume no constraints
  ! (this *must* be done before calling md_set_defaults)
  tmp_md%ndof = 3*cell%natoms

  ! read checkpoint data..
  call checkpoint_open('read')
  do i = 1, size(rand_seed)
    read(checkpoint_unit) rand_seed(i)
  end do
  call cell_read_checkpoint(cell)
  call md_type_read_checkpoint(md, cell)
  if (int_close) call checkpoint_close

  write(stdout,*) "Starting molecular dynamics continuation..."
  write(stdout,*) ""

  ! increment timestep by 1.. don't want to repeat last data point..
  md%timestep = md%timestep

  ! read in input params into temporary md type..
  call md_set_defaults(tmp_md)

  ! overwrite values by user - check for errors
  call md_read_input(tmp_md)

  ! now subtract 3 if linear momentum of COM is fixed (zero)
  if (tmp_md%fix_com) tmp_md%ndof = tmp_md%ndof - 3
  ! sanity check
  if (tmp_md%ndof .lt. 0) call io_err("md_continuation_init: Negative number of degrees for freedom")

  ! convert to atomic units
  tmp_md%dt = units_natural_to_atomic(time=tmp_md%dt)
  ! rescale temperature - always safe because we assign a default even with NVE...
  tmp_md%temperature = units_natural_to_atomic(temperature=tmp_md%temperature)
  tmp_md%langevin_time = units_natural_to_atomic(time=tmp_md%langevin_time)
  tmp_md%langevin_noise_stddev = sqrt((2.0_dp*tmp_md%temperature)/(tmp_md%langevin_time*tmp_md%dt)) ! sqrt(2kT/(dt*langevin_time))
  tmp_md%langevin_force_correction = 1.0_dp/(1.0_dp + 0.5_dp*tmp_md%dt/tmp_md%langevin_time)   ! 1/(1 + gamma*dt/2)

  ! now compare against what we already have..
  if (tmp_md%trajectory_units .ne. md%trajectory_units) then
    write(stdout,*) 'Warning: md_trajectory_units has changed on continuation'
    write(stdout,*) 'From:', md%trajectory_units, 'To:', tmp_md%trajectory_units
    write(stdout,*)
    md%trajectory_units = tmp_md%trajectory_units
  end if

  if (tmp_md%ntimesteps .ne. md%ntimesteps) then
    write(stdout,*) 'Warning: md_num_timesteps has changed on continuation'
    write(stdout,*) 'From:', md%ntimesteps, 'To:', tmp_md%ntimesteps
    write(stdout,*)
    md%ntimesteps = tmp_md%ntimesteps
  end if

  if (tmp_md%fix_com .neqv. md%fix_com) then
    write(stdout,*) 'Warning: md_fix_com has changed on continuation'
    write(stdout,*) 'From:', md%fix_com, 'To:', tmp_md%fix_com
    write(stdout,*)
    md%fix_com = tmp_md%fix_com
    if (abs(md%ndof - tmp_md%ndof) > 1) then
      call io_err("md_continuation_init: ndof should not change by more than 1 for change in fix_com")
    end if
    md%ndof = tmp_md%ndof
  else
    if (tmp_md%ndof .ne. md%ndof) then
      call io_err("md_continuation_init: ndof has somehow changed on continuation while fix_com is unchanged")
    end if
  end if

  if (tmp_md%dt .ne. md%dt) then
    write(stdout,*) 'Warning: md_timestep has changed on continuation'
    write(stdout,*) 'From:', units_atomic_to_natural(time=md%dt), &
                    & 'To:', units_atomic_to_natural(time=tmp_md%dt)
    write(stdout,*)
    md%dt = tmp_md%dt
  end if

  if (tmp_md%ensemble .ne. md%ensemble) then
    write(stdout,*) 'Warning: md_ensemble has changed on continuation'
    write(stdout,*) 'From:', md%ensemble, 'To:', tmp_md%ensemble
    write(stdout,*)
    md%ensemble = tmp_md%ensemble
  end if

  if (tmp_md%thermostat .ne. md%thermostat) then
    write(stdout,*) 'Warning: md_thermostat has changed on continuation'
    write(stdout,*) 'From:', md%thermostat, 'To:', tmp_md%thermostat
    write(stdout,*)
    md%thermostat = tmp_md%thermostat
  end if

  if (tmp_md%temperature .ne. md%temperature) then
    write(stdout,*) 'Warning: md_temperature has changed on continuation'
    write(stdout,*) 'From:', units_atomic_to_natural(temperature=md%temperature), &
                    & 'To:', units_atomic_to_natural(temperature=tmp_md%temperature)
    write(stdout,*)
    md%temperature = tmp_md%temperature
  end if

  if (tmp_md%output_interval .ne. md%output_interval) then
    write(stdout,*) 'Warning: md_output_interval has changed on continuation'
    write(stdout,*) 'From:', md%output_interval, 'To:', tmp_md%output_interval
    write(stdout,*)
    md%output_interval = tmp_md%output_interval
  end if

  if (tmp_md%checkpoint_interval .ne. md%checkpoint_interval) then
    write(stdout,*) 'Warning: md_checkpoint_interval has changed on continuation'
    write(stdout,*) 'From:', md%checkpoint_interval, 'To:', tmp_md%checkpoint_interval
    write(stdout,*)
    md%checkpoint_interval = tmp_md%checkpoint_interval
  end if

  if (tmp_md%langevin_time .ne. md%langevin_time) then
    write(stdout,*) 'Warning: md_langevin_time has changed on continuation'
    write(stdout,*) 'From:', units_atomic_to_natural(time=md%langevin_time), &
                    & 'To:', units_atomic_to_natural(time=tmp_md%langevin_time)
    write(stdout,*)
    md%trajectory_units = tmp_md%trajectory_units
    md%langevin_noise_stddev = tmp_md%langevin_noise_stddev
    md%langevin_force_correction = tmp_md%langevin_force_correction
  else
    if ((tmp_md%langevin_noise_stddev .ne. md%langevin_noise_stddev) .or. &
      & (tmp_md%langevin_force_correction .ne. md%langevin_force_correction)) then
      call io_err("md_continuation_init: Langevin derived parameters have somehow changed while time has remained constant")
    end if
  end if

  call md_echo_params(md)
end subroutine md_continuation

subroutine md_echo_params(md)
  use constants,  only: units_atomic_to_natural
  use io,         only: stdout
  implicit none
  type(md_type),  intent(in)  ::  md

  write(stdout,*) ""
  write(stdout,*) "Molecular dynamics parameters"
  write(stdout,*) "  Default md_num_timesteps = ", default_num_timesteps, " steps"
  write(stdout,*) "  Default md_timestep = ", default_timestep, " fs"
  write(stdout,*) "  Default md_ensemble = ", default_ensemble
  write(stdout,*) "  Default md_thermostat = ", default_thermostat
  write(stdout,*) "  Default md_temperature = ", default_temperature, " K"
  write(stdout,*) "  Default md_langevin_time = ", default_langevin_time, " fs"
  write(stdout,*) "  Default md_trajectory_units = ", default_trajectory_units
  write(stdout,*) "  Default md_output_interval = ", default_output_interval, " steps"
  write(stdout,*) "  Default md_checkpoint_interval = ", default_checkpoint_interval, " steps"
  write(stdout,*) ""
  write(stdout,*) "========== Simulation parameters =========="
  write(stdout,*) "  md_num_timesteps = ", md%ntimesteps, " steps"
  write(stdout,*) "  md_timestep = ", units_atomic_to_natural(time=md%dt), " fs"
  write(stdout,*) "  md_ensemble = ", md%ensemble
  write(stdout,*) "  md_thermostat = ", md%thermostat
  write(stdout,*) "  md_temperature = ", units_atomic_to_natural(temperature=md%temperature), " K"
  write(stdout,*) "  md_langevin_time = ", units_atomic_to_natural(time=md%langevin_time), " fs"
  write(stdout,*) "  md_trajectory_units = ", md%trajectory_units
  write(stdout,*) "  md_output_interval = ", md%output_interval, " steps"
  write(stdout,*) "  md_checkpoint_interval = ", md%checkpoint_interval, " steps"
  write(stdout,*) "  md_fix_com = ", md%fix_com
  write(stdout,*) ""
end subroutine md_echo_params

subroutine md_read_input(md)
  use io, only: io_input_get_single_value, io_query_keyword, io_err
  implicit none
  type(md_type),  intent(inout) ::  md

  call io_input_get_single_value('md_num_timesteps', md%ntimesteps)
  if (md%ntimesteps .lt. 1) call io_err("md_read_input: md_num_timesteps must be > 0")

  call io_input_get_single_value('md_timestep', md%dt)
  if (md%dt .le. 0.0_dp) call io_err("md_read_input: md_timestep must be > 0.0")

  call io_input_get_single_value('md_ensemble', md%ensemble)
  if ((md%ensemble .ne. 'nve') .and. (md%ensemble .ne. 'nvt')) then
    call io_err("md_read_input: md_ensemble must be nve or nvt")
  end if

  call io_input_get_single_value('md_thermostat', md%thermostat)
  if (md%thermostat .ne. 'langevin') call io_err("md_read_input: md_thermostat must be langevin")

  ! assume no temperature specified in input file..
  temperature_specified = .false.
  ! should only be 1... >1 is a block and will be caught by io_input_get_single_value
  if (io_query_keyword('md_temperature') .ne. 0) temperature_specified = .true.

  ! set default, or use value from input file..
  call io_input_get_single_value('md_temperature', md%temperature)
  if (md%temperature .le. 0.0_dp) call io_err("md_read_input: md_temperature must be > 0.0")

  call io_input_get_single_value('md_langevin_time', md%langevin_time)
  if (md%langevin_time .le. 0.0_dp) call io_err("md_read_input: md_langevin_time must be > 0.0")

  call io_input_get_single_value('md_fix_com', md%fix_com)

  call io_input_get_single_value('md_trajectory_units', md%trajectory_units)
  if ((md%trajectory_units .ne. 'atomic') .and. (md%trajectory_units .ne. 'user')) then
    call io_err("md_read_input: md_trajectory_units must be atomic or user")
  end if

  call io_input_get_single_value('md_output_interval', md%output_interval)
  if (md%output_interval .le. 0) call io_err("md_read_input: md_output_interval must be > 0")

  call io_input_get_single_value('md_checkpoint_interval', md%checkpoint_interval)
  if (md%output_interval .le. 0) call io_err("md_read_input: md_checkpoint_interval must be > 0")
end subroutine md_read_input

subroutine md_set_defaults(md)
  implicit none
  type(md_type),  intent(inout) ::  md
  logical ::  default_fix_com

  md%timestep = 0
  md%ntimesteps = default_num_timesteps
  md%dt = default_timestep
  md%ensemble = default_ensemble
  md%thermostat = default_thermostat
  md%temperature = default_temperature
  md%langevin_time = default_langevin_time

  ! md%ndof must be calculated before this routine called - at this stage, just 3*cell%natoms:
  if (md%ndof .eq. 3) then
    ! cell%natoms == 1:
    default_fix_com = .false.
  else
    ! cell%natoms > 1
    default_fix_com = .true.
  end if
  md%fix_com = default_fix_com

  md%trajectory_units = default_trajectory_units
  md%output_interval = default_output_interval
  md%checkpoint_interval = default_checkpoint_interval
end subroutine md_set_defaults

subroutine md_allocate(md, cell)
  use io, only: io_err
  implicit none
  type(md_type),    intent(inout) ::  md
  type(cell_type),  intent(in)    ::  cell
  integer                         ::  istat

  ! allocate the position, momentum and force arrays
  allocate(md%pos(3, cell%natoms), stat=istat)
  if (istat .ne. 0) call io_err("md_allocate: Could not allocate md%pos array")

  allocate(md%momentum(3, cell%natoms), stat=istat)
  if (istat .ne. 0) call io_err("md_allocate: Could not allocate md%momentum array")

  allocate(md%force(3, cell%natoms), stat=istat)
  if (istat .ne. 0) call io_err("md_allocate: Could not allocate md%force array")
end subroutine md_allocate

subroutine md_deallocate(md)
  use io, only: io_err
  implicit none
  type(md_type),    intent(inout) ::  md
  integer                         ::  istat

  ! deallocate the position, momentum and force arrays
  if (allocated(md%pos)) then
    deallocate(md%pos, stat=istat)
    if (istat .ne. 0) call io_err("md_deallocate: Could not deallocate md%pos array")
  end if

  if (allocated(md%momentum)) then
    deallocate(md%momentum, stat=istat)
    if (istat .ne. 0) call io_err("md_deallocate: Could not deallocate md%momentum array")
  end if

  if (allocated(md%force)) then
    deallocate(md%force, stat=istat)
    if (istat .ne. 0) call io_err("md_deallocate: Could not deallocate md%force array")
  end if
end subroutine md_deallocate

subroutine md_copy(md_in, cell, md_out)
  implicit none
  type(md_type),    intent(in)    ::  md_in
  type(cell_type),  intent(in)    ::  cell
  type(md_type),    intent(inout) ::  md_out

  ! assumes md_in is set up correctly (md_init called on this)
  ! md_allocate needs cell%natoms
  call md_allocate(md_out, cell)

  ! then do deep copy...
  md_out%output_filename = md_in%output_filename
  md_out%trajectory_units = md_in%trajectory_units
  md_out%timestep = md_in%timestep
  md_out%ntimesteps = md_in%ntimesteps
  md_out%ndof = md_in%ndof
  md_out%dt = md_in%dt
  md_out%kinetic_energy = md_in%kinetic_energy
  md_out%potential_energy = md_in%potential_energy
  md_out%inst_temp = md_in%inst_temp
  md_out%ensemble = md_in%ensemble
  md_out%thermostat = md_in%thermostat
  md_out%temperature = md_in%temperature
  md_out%fix_com = md_in%fix_com
  md_out%langevin_time = md_in%langevin_time
  md_out%langevin_noise_stddev = md_in%langevin_noise_stddev
  md_out%langevin_force_correction = md_in%langevin_force_correction

  md_out%pos(:,:) = md_in%pos(:,:)
  md_out%momentum(:,:) = md_in%momentum(:,:)
  md_out%force(:,:) = md_in%force(:,:)
end subroutine md_copy

subroutine md_initial_momenta(md, cell)
  use algor,  only: algor_gauss_rand
  implicit none
  type(md_type),    intent(inout) ::  md
  type(cell_type),  intent(in)    ::  cell
  ! local variables
  real(kind=dp),  dimension(3)  ::  com_momentum
  real(kind=dp)                 ::  sqrt_temp, temp_rescale
  real(kind=dp)                 ::  thermal_stddev
  integer                       ::  iatom, icomp

  if ((md%ensemble .eq. 'nvt') .or. ((md%ensemble .eq. 'nve') .and. (temperature_specified))) then
    ! if temperature was given in input file, initialize momenta at this temp, even for NVE

    sqrt_temp = sqrt(md%temperature)
    ! draw momenta from a Gaussian distribution for initialization
    ! mean: 0, stddev = sqrt(kT*mass) : kb = 1
    com_momentum(:) = 0.0_dp
    do iatom = 1, cell%natoms
      thermal_stddev = sqrt(cell%atom_mass(iatom))*sqrt_temp
      do icomp = 1,3
        md%momentum(icomp, iatom) = algor_gauss_rand(mean=0.0_dp, stddev=thermal_stddev)
      end do
      com_momentum(:) = com_momentum(:) + md%momentum(:,iatom)
    end do
    com_momentum(:) = com_momentum(:)/real(cell%natoms, dp)

    ! remove any centre of mass momentum..
    if (md%fix_com) then
      do iatom = 1, cell%natoms
        md%momentum(:,iatom) = md%momentum(:,iatom)-com_momentum(:)
      end do
    end if

    ! calculate actual temperature..
    call md_calc_kinetic_energy(md, cell)
    if (md%ndof .ne. 0) then
      md%inst_temp = 2.0_dp*md%kinetic_energy/real(md%ndof,dp)
    else
      md%inst_temp = 0.0_dp   ! always, because KE = 0, but avoids NaN
    end if

    ! now rescale to get correct temperature..
    if (md%inst_temp .ne. md%temperature) then
      ! make sure no divide by zero.. don't bother rescaling if inst_temp is too small..
      if (md%inst_temp .gt. epsilon(1.0_dp)) then
        temp_rescale = sqrt(md%temperature/md%inst_temp)

        do iatom = 1, cell%natoms
          md%momentum(:,iatom) = md%momentum(:,iatom)*temp_rescale
        end do
      end if
    end if

  else if (md%ensemble .eq. 'nve') then
    ! temperature not given in input file - md%temperature is meaningless (default)
    ! so assume all energy is potential due to atomic positions.. no initial momenta
    md%momentum(:,:) = 0.0_dp
  end if
end subroutine md_initial_momenta

subroutine do_md(md, cell)
  use checkpoint, only: continuation
  implicit none
  type(md_type),    intent(inout) ::  md
  type(cell_type),  intent(inout) ::  cell
  integer :: t, t_start

  if (.not. continuation) call md_write_output(md, cell)

  t_start = md%timestep + 1 ! when not continuation, md%timestep = 0

  do t = t_start, md%ntimesteps

    call md_velocity_verlet_p1(md, cell)

    ! keep coordinates in sync before calculating forces
    call md_sync_coords(md, cell)

    call md_get_forces(md, cell)

    ! add random noise to forces if using langevin
    if ((md%ensemble .eq. 'nvt') .and. (md%thermostat .eq. 'langevin')) then
      call md_langevin_force(md, cell)
    end if

    ! make sure no net force
    if (md%fix_com) then
      call md_zero_net_force(md, cell)
    end if

    call md_velocity_verlet_p2(md, cell)

    ! end of timestep - can calculate all other quantities
    call md_calc_kinetic_energy(md, cell)
    if (md%ndof .ne. 0) then
      md%inst_temp = 2.0_dp*md%kinetic_energy/real(md%ndof,dp)
    else
      md%inst_temp = 0.0_dp   ! always, because KE = 0, but avoids NaN
    end if

    md%timestep = t
    if (mod(t, md%output_interval) .eq. 0) call md_write_output(md, cell)
    if (mod(t, md%checkpoint_interval) .eq. 0) call md_write_checkpoint(md, cell)
  end do

  ! write final checkpoint information
  call md_write_checkpoint(md, cell)
end subroutine do_md

subroutine md_velocity_verlet_p1(md, cell)
  implicit none
  type(md_type),    intent(inout)   ::  md
  type(cell_type),  intent(inout)   ::  cell
  integer ::  iatom

  do iatom = 1, cell%natoms
    ! mv(t+dt/2) = mv(t) + F(t)*dt/2
    md%momentum(:,iatom) = md%momentum(:,iatom) + 0.5_dp*md%force(:,iatom)*md%dt
    ! r(t+dt) = r(t) + v(t+dt/2)*dt
    md%pos(:,iatom) = md%pos(:,iatom) + md%momentum(:,iatom)*md%dt/cell%atom_mass(iatom)
  end do
end subroutine md_velocity_verlet_p1

subroutine md_velocity_verlet_p2(md, cell)
  implicit none
  type(md_type),    intent(inout) ::  md
  type(cell_type),  intent(in)    ::  cell
  integer ::  iatom

  do iatom = 1, cell%natoms
    ! mv(t+dt) = mv(t+dt/2) + F(t+dt)*dt/2
    md%momentum(:,iatom) = md%momentum(:,iatom) + 0.5_dp*md%force(:,iatom)*md%dt
  end do
end subroutine md_velocity_verlet_p2

subroutine md_sync_coords(md, cell)
  use cell, only: cell_shift_to_unit_cell, cell_cart_to_frac, cell_frac_to_cart
  implicit none
  type(md_type),    intent(in)    ::  md
  type(cell_type),  intent(inout) ::  cell

  ! potentially slow because of copy
  ! (could instead create new cell routines which take in position array...)
  cell%atom_cart_pos(:,:) = md%pos(:,:)

  call cell_cart_to_frac(cell)        ! update fractionals
  call cell_shift_to_unit_cell(cell)  ! apply PBCs (only works on fractionals)
  call cell_frac_to_cart(cell)        ! update Cartesian after PBCs from corrected fractionals
end subroutine md_sync_coords

subroutine md_get_forces(md, cell)
  use potential,  only: pot_get_forces
  implicit none
  type(md_type),    intent(inout) ::  md
  type(cell_type),  intent(in)    ::  cell

  call pot_get_forces(cell, md%force, md%potential_energy)
end subroutine md_get_forces

subroutine md_zero_net_force(md, cell)
  implicit none
  type(md_type),    intent(inout) ::  md
  type(cell_type),  intent(in)    ::  cell
  ! local vars
  real(kind=dp),  dimension(3)    ::  net_force
  integer ::  iatom

  net_force(:) = 0.0_dp
  do iatom = 1, cell%natoms
    net_force(:) = net_force(:) + md%force(:,iatom)
  end do

  ! distribute among each atom..
  net_force(:) = net_force(:)/real(cell%natoms, dp)

  do iatom = 1, cell%natoms
    md%force(:,iatom) = md%force(:,iatom) - net_force(:)
  end do
end subroutine md_zero_net_force

subroutine md_langevin_force(md, cell)
  use algor,  only: algor_gauss_rand
  implicit none
  type(md_type),    intent(inout) ::  md
  type(cell_type),  intent(in)    ::  cell
  real(kind=dp) ::  langevin_gamma, sqrt_mass
  integer ::  iatom, j

  langevin_gamma = 1.0_dp/md%langevin_time

  do iatom = 1, cell%natoms
    sqrt_mass = sqrt(cell%atom_mass(iatom))
    do j = 1,3
      ! F = F - gamma*mv + sqrt(2mkT/(tau*dt))*N(0,1)
      md%force(j,iatom) = md%force(j,iatom) - langevin_gamma*md%momentum(j,iatom) &
      &                   + sqrt_mass*algor_gauss_rand(mean=0.0_dp, stddev=md%langevin_noise_stddev)

      ! correct for only having v(t + dt/2) multiplying by 1/(1 + gamma*dt/2)
      md%force(j,iatom) = md%langevin_force_correction*md%force(j,iatom)
    end do
  end do
end subroutine md_langevin_force

subroutine md_calc_kinetic_energy(md, cell)
  implicit none
  type(md_type),    intent(inout) ::  md
  type(cell_type),  intent(in)    ::  cell
  integer ::  iatom

  md%kinetic_energy = 0.0_dp
  do iatom = 1, cell%natoms
    md%kinetic_energy = md%kinetic_energy + 0.5_dp*dot_product(md%momentum(:,iatom), md%momentum(:,iatom))/cell%atom_mass(iatom)
  end do
end subroutine md_calc_kinetic_energy

subroutine  md_write_output(md, cell)
  use io,         only: io_open_file, io_close_file
  use constants,  only: element_symbol, units => units_atomic_to_natural
  implicit none
  type(md_type),    intent(in)  ::  md
  type(cell_type),  intent(in)  ::  cell
  real(kind=dp) ::  tot_energy
  integer ::  i, output_unit_num

  if (md%timestep .eq. 0) then
    call io_open_file(md%output_filename, output_unit_num, 'replace')
    write (output_unit_num, 99) 'BEGIN header'
    write (output_unit_num, *)   cell%natoms, md%ntimesteps+1
    write (output_unit_num, 99) 'END header'
    write (output_unit_num, *)  ' '
  else
    call io_open_file(md%output_filename, output_unit_num, 'append')
  end if

  tot_energy = md%kinetic_energy+md%potential_energy

  if (md%trajectory_units .eq. 'user') then
    ! rescale units
    write (output_unit_num, 1) units(time=real(md%timestep,dp)*md%dt)               ! time
    write (output_unit_num, 2) units(energy=md%potential_energy), units(energy=tot_energy), &
                                & units(energy=md%kinetic_energy)                   ! PE, TE, KE
    write (output_unit_num, 3) units(temperature=md%inst_temp)                      ! temperature

    ! lattice vectors (lines are a, b, c)
    do i = 1, 3
      write (output_unit_num, 4) units(length=cell%lattice_vectors(:,i))
    end do

    ! positions
    do i = 1, cell%natoms
      write (output_unit_num, 5) element_symbol(cell%atom_species(i)), i, units(length=md%pos(:,i))
    end do

    ! velocities
    do i = 1, cell%natoms
      write (output_unit_num, 6) element_symbol(cell%atom_species(i)), i, units(velocity=(md%momentum(:,i)/cell%atom_mass(i)))
    end do

    ! forces
    do i = 1, cell%natoms
      write (output_unit_num, 7) element_symbol(cell%atom_species(i)), i, units(force=md%force(:,i))
    end do


  ! atomic units - no change
  else
    write (output_unit_num, 1) real(md%timestep,dp)*md%dt                           ! time
    write (output_unit_num, 2) md%potential_energy, tot_energy, md%kinetic_energy   ! PE, TE, KE
    write (output_unit_num, 3) md%inst_temp                                         ! temperature

    ! lattice vectors (lines are a, b, c)
    do i = 1, 3
      write (output_unit_num, 4) cell%lattice_vectors(:,i)
    end do

    ! positions
    do i = 1, cell%natoms
      write (output_unit_num, 5) element_symbol(cell%atom_species(i)), i, md%pos(:,i)
    end do

    ! velocities
    do i = 1, cell%natoms
      write (output_unit_num, 6) element_symbol(cell%atom_species(i)), i, md%momentum(:,i)/cell%atom_mass(i)
    end do

    ! forces
    do i = 1, cell%natoms
      write (output_unit_num, 7) element_symbol(cell%atom_species(i)), i, md%force(:,i)
    end do
  end if

  write (output_unit_num, *) ' '

  call io_close_file(output_unit_num)

  1   format(21x,es24.16e3)                           ! time
  2   format(18x,3(3x,es24.16e3),'  <-- E')           ! energy
  3   format(21x,es24.16e3,T100, '  <-- T')           ! temperature
  4   format(18x,3(3x,es24.16e3),'  <-- h')           ! cell vectors (line is a vector)
  5   format(1x,a8,1x,i8,3(3x,es24.16e3),'  <-- R')   ! position
  6   format(1x,a8,1x,i8,3(3x,es24.16e3),'  <-- V')   ! velocity
  7   format(1x,a8,1x,i8,3(3x,es24.16e3),'  <-- F')   ! force
  99  format(1x,a)                                    ! headers
end subroutine md_write_output

subroutine md_write_checkpoint(md, cell)
  use checkpoint, only: checkpoint_open, checkpoint_unit, checkpoint_close
  use algor,      only: rand_seed
  use cell,       only: cell_write_checkpoint
  implicit none
  type(md_type),    intent(in)  ::  md
  type(cell_type),  intent(in)  ::  cell

  call checkpoint_open('write')
  write(checkpoint_unit) rand_seed
  call cell_write_checkpoint(cell)
  call md_type_write_checkpoint(md)
  call checkpoint_close
end subroutine md_write_checkpoint

subroutine md_type_read_checkpoint(md, cell)
  use checkpoint, only: checkpoint_unit
  implicit none
  type(md_type),    intent(inout) ::  md
  type(cell_type),  intent(in)    ::  cell
  integer :: i, j

  ! always safe to call md_deallocate... will only dealloc if necessary
  ! might be slow if not required though..
  call md_deallocate(md)
  call md_allocate(md, cell)
  read(checkpoint_unit)  md%output_filename
  read(checkpoint_unit)  md%trajectory_units
  read(checkpoint_unit)  md%timestep
  read(checkpoint_unit)  md%ntimesteps
  read(checkpoint_unit)  md%ndof
  read(checkpoint_unit)  md%dt
  read(checkpoint_unit)  md%kinetic_energy
  read(checkpoint_unit)  md%potential_energy
  read(checkpoint_unit)  md%inst_temp
  read(checkpoint_unit)  md%ensemble
  read(checkpoint_unit)  md%thermostat
  read(checkpoint_unit)  md%temperature
  read(checkpoint_unit)  md%fix_com
  read(checkpoint_unit)  md%output_interval
  read(checkpoint_unit)  md%checkpoint_interval
  read(checkpoint_unit)  md%langevin_time
  read(checkpoint_unit)  md%langevin_noise_stddev
  read(checkpoint_unit)  md%langevin_force_correction
  do i = 1, cell%natoms
    do j = 1, 3
      read(checkpoint_unit)  md%pos(j,i)
    end do
  end do
  do i = 1, cell%natoms
    do j = 1, 3
      read(checkpoint_unit)  md%momentum(j,i)
    end do
  end do
  do i = 1, cell%natoms
    do j = 1, 3
      read(checkpoint_unit)  md%force(j,i)
    end do
  end do
end subroutine md_type_read_checkpoint

subroutine md_type_write_checkpoint(md)
  use checkpoint, only: checkpoint_unit
  implicit none
  type(md_type),    intent(in)  ::  md
  integer :: i, j

  write(checkpoint_unit)  md%output_filename
  write(checkpoint_unit)  md%trajectory_units
  write(checkpoint_unit)  md%timestep
  write(checkpoint_unit)  md%ntimesteps
  write(checkpoint_unit)  md%ndof
  write(checkpoint_unit)  md%dt
  write(checkpoint_unit)  md%kinetic_energy
  write(checkpoint_unit)  md%potential_energy
  write(checkpoint_unit)  md%inst_temp
  write(checkpoint_unit)  md%ensemble
  write(checkpoint_unit)  md%thermostat
  write(checkpoint_unit)  md%temperature
  write(checkpoint_unit)  md%fix_com
  write(checkpoint_unit)  md%output_interval
  write(checkpoint_unit)  md%checkpoint_interval
  write(checkpoint_unit)  md%langevin_time
  write(checkpoint_unit)  md%langevin_noise_stddev
  write(checkpoint_unit)  md%langevin_force_correction
  do i = 1, size(md%pos,2)
    do j = 1, 3
      write(checkpoint_unit)  md%pos(j,i)
    end do
  end do
  do i = 1, size(md%momentum,2)
    do j = 1, 3
      write(checkpoint_unit)  md%momentum(j,i)
    end do
  end do
  do i = 1, size(md%force,2)
    do j = 1, 3
      write(checkpoint_unit)  md%force(j,i)
    end do
  end do
end subroutine md_type_write_checkpoint
end module md
