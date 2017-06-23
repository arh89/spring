module pathint

  use constants,  only: dp
  use cell,       only: cell_type
  use md,         only: md_type
  implicit none

  private

  ! fftw
  integer,  parameter ::  i64 = selected_int_kind(18)
  integer,  parameter ::  fftw_forward = -1
  integer,  parameter ::  fftw_backward = +1
  integer,  parameter ::  fftw_measure = 0

  integer,            parameter ::  default_num_beads = 1
  logical,            parameter ::  default_bead_trajectories = .false.
  logical,            parameter ::  default_bead_interactions = .true.
  character(len=11),  parameter ::  default_propagator_variables = 'primitive'
  real(kind=dp),      parameter ::  default_adiabaticity = 1.0_dp

  type, public  ::  bead_type
    integer           ::  prev_bead_index
    integer           ::  next_bead_index
    type(cell_type)   ::  cell
    type(md_type)     ::  md
  end type bead_type

  type, public  ::  ring_polymer_type
    logical                                           ::  bead_interactions
    real(kind=dp)                                     ::  adiabaticity
    character(len=11)                                 ::  propagator_variables  ! 'primitive' or 'normalmodes'
    logical                                           ::  output_bead_trajectories
    integer                                           ::  nbeads
    real(kind=dp)                                     ::  chain_frequency
    real(kind=dp),    dimension(:),     allocatable   ::  nm_eigenvalues
    type(bead_type),  dimension(:),     allocatable   ::  beads
    logical                                           ::  nm_thermostat_centroid

    ! work arrays
    integer(kind=i64)                                 ::  fftw_plan_forward
    integer(kind=i64)                                 ::  fftw_plan_backward
    real(kind=dp),    dimension(:),     allocatable   ::  fftw_real_n
    complex(kind=dp), dimension(:),     allocatable   ::  fftw_cmplx_in_n
    complex(kind=dp), dimension(:),     allocatable   ::  fftw_cmplx_out_n
    complex(kind=dp), dimension(:),     allocatable   ::  fftw_cmplx_out_no2p1
    real(kind=dp),    dimension(:,:),   allocatable   ::  beads_langevin_force

    ! thermodynamic estimators
    real(kind=dp)                                     ::  primitive_energy_est
    real(kind=dp)                                     ::  virial_energy_est
  end type ring_polymer_type

  public  ::  pathint_init
  public  ::  pathint_continuation
  public  ::  do_pimd
  public  ::  pathint_deallocate

contains

subroutine pathint_echo_params(ring_polymer)
  use io, only: stdout
  implicit none
  type(ring_polymer_type),  intent(in)  ::  ring_polymer

  write(stdout,*) ""
  write(stdout,*) "Path integral molecular dynamics parameters"
  write(stdout,*) "  Default pathint_num_beads = ", default_num_beads
  write(stdout,*) "  Default pathint_bead_trajectories = ", default_bead_trajectories
  write(stdout,*) "  Default pathint_bead_interactions = ", default_bead_interactions
  write(stdout,*) "  Default pathint_propagator_variables = ", default_propagator_variables
  write(stdout,*) "  Default pathint_adiabaticity = ", default_adiabaticity
  write(stdout,*) ""
  write(stdout,*) "========== Simulation parameters =========="
  write(stdout,*) "  pathint_num_beads = ", ring_polymer%nbeads
  write(stdout,*) "  pathint_bead_trajectories = ", ring_polymer%output_bead_trajectories
  write(stdout,*) "  pathint_bead_interactions = ", ring_polymer%bead_interactions
  write(stdout,*) "  pathint_propagator_variables = ", ring_polymer%propagator_variables
  write(stdout,*) "  pathint_adiabaticity = ", ring_polymer%adiabaticity
  if (ring_polymer%propagator_variables .eq. 'normalmodes') then
    write(stdout,*) "  pathint_nm_thermostat_centroid = ", ring_polymer%nm_thermostat_centroid
  end if
  write(stdout,*) ""
end subroutine pathint_echo_params

subroutine pathint_init(ring_polymer, centroid_cell, centroid_md)
  use cell,       only: cell_copy
  use md,         only: md_init, md_copy, md_initial_momenta, md_get_forces, md_zero_net_force, &
                      & md_calc_kinetic_energy
  use io,         only: io_err
  use checkpoint, only: continuation
  implicit none
  type(ring_polymer_type),  intent(inout) ::  ring_polymer
  type(cell_type),          intent(inout) ::  centroid_cell
  type(md_type),            intent(inout) ::  centroid_md
  character(len=12) ::  file_ext
  integer :: ibead, istat, ndof

  if (continuation) call io_err("pathint_init: pathint_init called rather than pathint_continuation")
  ! centroid_cell already initialized, let's initialize the centroid_md type here too
  ! (including reading in the md parameters etc)

  ! md_init will initialize (centroid) momenta as well as some other things
  ! ignore this for now - we recalc centroid properties anyway..
  ! md_init will continue centroid_md and centroid_cell if possible..
  call md_init(centroid_md, centroid_cell)

  ! might need to enfore NVT?..
  call pathint_set_defaults(ring_polymer)
  call pathint_read_input(ring_polymer)

  ! now echo pathint only params
  call pathint_echo_params(ring_polymer)

  ! allocate beads array and copy everything:
  allocate(ring_polymer%beads(ring_polymer%nbeads), stat=istat)
  if (istat .ne. 0) call io_err("pathint_init: Could not allocate beads array")

  do ibead = 1, ring_polymer%nbeads
    ! copy cells
    call cell_copy(centroid_cell, ring_polymer%beads(ibead)%cell)
    ! copy md
    call md_copy(centroid_md, centroid_cell, ring_polymer%beads(ibead)%md)

    ! set index of next and previous beads in chain
    ring_polymer%beads(ibead)%prev_bead_index = ibead - 1
    ring_polymer%beads(ibead)%next_bead_index = ibead + 1
  end do

  ! enforce boundary conditions on bead indices
  ring_polymer%beads(1)%prev_bead_index = ring_polymer%nbeads
  ring_polymer%beads(ring_polymer%nbeads)%next_bead_index = 1

  ! chain frequency: w_p = sqrt(P)/(hbar*beta) = sqrt(P)*kT/hbar
  ring_polymer%chain_frequency = sqrt(real(ring_polymer%nbeads, dp))*centroid_md%temperature

  ! change output filenames in each bead's md object (regardless of whether we output or not)
  do ibead = 1, ring_polymer%nbeads
    write(file_ext, '(I11)') ibead
    file_ext = '_'//adjustl(file_ext)
    ring_polymer%beads(ibead)%md%output_filename = trim(ring_polymer%beads(ibead)%md%output_filename)//trim(file_ext)
  end do

  ! Bookkeeping finished.. now can focus on physics..
  if (ring_polymer%propagator_variables .eq. 'normalmodes') then
    call pathint_normal_modes_init(ring_polymer)
  end if

  ! everything in primitive modes at the moment..
  ! set initial bead positions (spread out from centroid)
  call pathint_spread_beads(ring_polymer, centroid_cell, centroid_md)

  ! change masses (for momenta)
  call pathint_rescale_masses(ring_polymer, 'forward')

  ! initial bead momenta and forces
  do ibead = 1, ring_polymer%nbeads
    ! should work for NVT and NVE (based on same routines for md)
    call md_initial_momenta(ring_polymer%beads(ibead)%md, ring_polymer%beads(ibead)%cell)
    call md_get_forces(ring_polymer%beads(ibead)%md, ring_polymer%beads(ibead)%cell)
  end do


  ! forces currently on primitive modes, transform coordinates..
  if (ring_polymer%propagator_variables .eq. 'normalmodes') then
    call pathint_normal_modes_pos_forward(ring_polymer)
    call pathint_normal_modes_force_forward(ring_polymer)
  end if

  ! nearest neighbour forces (need to be using correct variables)
  ! should always be true except for testing..
  if (ring_polymer%bead_interactions) then
    ! have rescaled masses, so rescale back (physical masses)
    call pathint_rescale_masses(ring_polymer, 'backward')
    call pathint_bead_neighbour_interactions(ring_polymer)
  end if
 
  ! zero net force for each replica..
  do ibead = 1, ring_polymer%nbeads
    if (ring_polymer%beads(ibead)%md%fix_com) then
      call md_zero_net_force(ring_polymer%beads(ibead)%md, ring_polymer%beads(ibead)%cell)
    end if
  end do
  
  ! which masses do we need here? physical or rescaled? centroid mass = physical so won't matter too much..
  ! calc KE for each bead (not centroid, this is averaged from each bead)
  do ibead = 1, ring_polymer%nbeads
    call md_calc_kinetic_energy(ring_polymer%beads(ibead)%md, ring_polymer%beads(ibead)%cell)
    ndof = ring_polymer%beads(ibead)%md%ndof
    if (ndof .ne. 0) then
      ring_polymer%beads(ibead)%md%inst_temp = 2.0_dp*ring_polymer%beads(ibead)%md%kinetic_energy/real(ndof,dp)
    else
      ring_polymer%beads(ibead)%md%inst_temp = 0.0_dp ! because KE = 0, but avoids NaN
    end if
  end do

  ! transform back to primitive modes if necessary.. (for calculation of estimators)
  if (ring_polymer%propagator_variables .eq. 'normalmodes') then
    call pathint_normal_modes_pos_backward(ring_polymer)

    ! can finally calculate centroid properties..
    call pathint_calc_centroid(ring_polymer, centroid_cell, centroid_md)
  end if

  ! calc estimators
  call pathint_calc_estimators(ring_polymer, centroid_cell, centroid_md)
end subroutine pathint_init

subroutine pathint_read_input(ring_polymer)
  use io, only: io_input_get_single_value, io_err
  implicit none
  type(ring_polymer_type),  intent(inout) ::  ring_polymer

  call io_input_get_single_value('pathint_num_beads', ring_polymer%nbeads)
  if (ring_polymer%nbeads .lt. 1) call io_err("pathint_read_input: nbeads must be >= 1")

  call io_input_get_single_value('pathint_bead_trajectories', ring_polymer%output_bead_trajectories)

  call io_input_get_single_value('pathint_propagator_variables', ring_polymer%propagator_variables)
  if ((ring_polymer%propagator_variables .ne. 'primitive') .and. &
    & (ring_polymer%propagator_variables .ne. 'normalmodes')) then
      call io_err("pathint_read_input: Unsupported propagator variable transformation")
  end if

  call io_input_get_single_value('pathint_adiabaticity', ring_polymer%adiabaticity)
  if ((ring_polymer%adiabaticity .lt. 0.0_dp) .or. (ring_polymer%adiabaticity .gt. 1.0_dp)) then
    call io_err("pathint_read_input: Adiabaticity value must be 0 < gamma < 1")
  end if

  if (ring_polymer%adiabaticity .ne. default_adiabaticity) then
    ! if doing CMD, cannot use primitive variables..
    if (ring_polymer%propagator_variables .eq. 'primitive') then
        call io_err("pathint_read_input: Cannot set adiabaticity if pathint_propagator_variables is not 'normalmodes'")
    end if
  end if

  ! by default we thermostat the centroid in CMD...
  ! but can override just in case do/don't want to do this (testing)
  call io_input_get_single_value('pathint_nm_thermostat_centroid', ring_polymer%nm_thermostat_centroid, .true.)

  ! if primitive, we must thermostat the centroid..
  if ((ring_polymer%propagator_variables .eq. 'primitive') .and. (.not. ring_polymer%nm_thermostat_centroid)) then
    call io_err("pathint_read_input: If pathint_propagator_variables is 'primitive' the centroid must be thermostatted")
  end if

  call io_input_get_single_value('pathint_bead_interactions', ring_polymer%bead_interactions)
end subroutine pathint_read_input

subroutine pathint_set_defaults(ring_polymer)
  implicit none
  type(ring_polymer_type),  intent(inout) ::  ring_polymer

  ring_polymer%nbeads = default_num_beads
  ring_polymer%output_bead_trajectories = default_bead_trajectories
  ring_polymer%propagator_variables = default_propagator_variables
  ring_polymer%adiabaticity = default_adiabaticity
  ring_polymer%bead_interactions = default_bead_interactions
  ring_polymer%nm_thermostat_centroid = .true.
end subroutine pathint_set_defaults

subroutine pathint_spread_beads(ring_polymer, centroid_cell, centroid_md)
  use constants,  only: two_pi
  use md,         only: md_sync_coords
  use algor,      only: algor_uniform_rand, algor_3d_rotation
  implicit none
  type(ring_polymer_type),  intent(inout) ::  ring_polymer
  type(cell_type),          intent(in)    ::  centroid_cell
  type(md_type),            intent(in)    ::  centroid_md

  real(kind=dp),  dimension(3)  ::  centroid_pos, bead_pos, rotation
  real(kind=dp)                 ::  theta, dtheta, rgy
  integer                       ::  ibead, iatom

  ! only do if nbeads > 1
  if (ring_polymer%nbeads .eq. 1) return

  ! divide 2pi radians evenly among beads
  dtheta = two_pi/real(ring_polymer%nbeads,dp)

  ! for each atom, distribute beads
  do iatom = 1, centroid_cell%natoms

    centroid_pos(:) = centroid_cell%atom_cart_pos(:,iatom)

    ! free particle rgy = hbar/sqrt(12mkT)
    ! and bead-bead distance = sqrt(beta*hbar^2/Pm)
    rgy = 1.0_dp/sqrt(12.0_dp*centroid_cell%atom_mass(iatom)*centroid_md%temperature)

    ! choose 3 angles randomly between 0 -> 2*pi
    rotation(:) = (/ algor_uniform_rand(0.0_dp,two_pi), &
                   & algor_uniform_rand(0.0_dp,two_pi), &
                   & algor_uniform_rand(0.0_dp,two_pi) /)

    do ibead = 1, ring_polymer%nbeads

      ! evenly distribute beads around the xy plane
      ! polar coordinates:
      !     theta = 0..2*pi
      ! x = r*cos(theta)
      ! y = r*sin(theta)
      theta = real(ibead,dp)*dtheta

      bead_pos(1) = rgy*cos(theta)
      bead_pos(2) = rgy*sin(theta)
      bead_pos(3) = 0.0_dp

      ! then perform a rotation in 3D space (no need to rotate around z axis, but no harm in doing this)
      bead_pos(:) = centroid_pos + algor_3d_rotation(bead_pos(:), rotation(:))
      ring_polymer%beads(ibead)%md%pos(:,iatom) = bead_pos(:)
    end do
  end do

  ! cartesian positions should now be correct, resync cell and md positions array
  ! note: if we go beyond a cell boundary, md%pos will track this
  do ibead = 1, ring_polymer%nbeads
    call md_sync_coords(ring_polymer%beads(ibead)%md, ring_polymer%beads(ibead)%cell)
  end do
end subroutine pathint_spread_beads

subroutine pathint_normal_modes_init(ring_polymer)
  use constants,  only: two_pi
  use io,         only: io_err
  implicit none
  type(ring_polymer_type),    intent(inout) ::  ring_polymer
  ! local vars:
  real(kind=dp) ::  lambda
  integer ::  i, istat

  ! die if nbeads is not even..
  if (mod(ring_polymer%nbeads,2) .ne. 0) then
    call io_err("pathint_normal_modes_init: normal modes only supports even number of beads")
  end if

  ! allocate work arrays
  ! fftw
  allocate(ring_polymer%fftw_real_n(ring_polymer%nbeads), stat=istat)
  if (istat .ne. 0) call io_err("pathint_normal_modes_init: Could not allocate fftw_real_n array")

  allocate(ring_polymer%fftw_cmplx_in_n(ring_polymer%nbeads), stat=istat)
  if (istat .ne. 0) call io_err("pathint_normal_modes_init: Could not allocate fftw_cmplx_in_n array")
  
  allocate(ring_polymer%fftw_cmplx_out_n(ring_polymer%nbeads), stat=istat)
  if (istat .ne. 0) call io_err("pathint_normal_modes_init: Could not allocate fftw_cmplx_out_n array")

  allocate(ring_polymer%fftw_cmplx_out_no2p1(ring_polymer%nbeads/2 + 1), stat=istat)
  if (istat .ne. 0) call io_err("pathint_normal_modes_init: Could not allocate fftw_cmplx_out_no2p1 array")

  ! eigenvalues
  allocate(ring_polymer%nm_eigenvalues(ring_polymer%nbeads), stat=istat)
  if (istat .ne. 0) call io_err("pathint_normal_modes_init: Could not allocate nm_eigenvalues array")

  !-----------------------------------------------------------------------------
  ! From: Cao & Martyna, JCP 104, 2028 (1996), Appendix:
  !
  ! Free particle eigenvalue values of each normal modes:
  !
  ! 1. The eigenvalue of the (1)st mode is zero.
  !
  ! 2. The eigenvalue of the (P)th mode is 4P.
  !
  ! 3. The eigenvalues of the (2I-2)th and the (2I-1)th modes are 
  !     4P[1-cos(2pi[I-1]/P)]
  !-----------------------------------------------------------------------------

  ! store nm eigenvalues:
  ring_polymer%nm_eigenvalues(1) = 0.0_dp
  ring_polymer%nm_eigenvalues(ring_polymer%nbeads) = 4.0_dp*real(ring_polymer%nbeads,dp)

  do i = 2, (ring_polymer%nbeads/2)
    lambda = 4.0_dp*real(ring_polymer%nbeads,dp)*(1.0_dp - cos( (two_pi*real((i-1),dp))/real(ring_polymer%nbeads,dp) ))
    ring_polymer%nm_eigenvalues(2*i-2) = lambda
    ring_polymer%nm_eigenvalues(2*i-1) = lambda
  end do

  ! initialize fftw plans (use measure rather than estimate to optimize for performance)
  call dfftw_plan_dft_r2c_1d(ring_polymer%fftw_plan_forward, ring_polymer%nbeads, &
                          & ring_polymer%fftw_real_n, ring_polymer%fftw_cmplx_out_no2p1, fftw_measure)
  call dfftw_plan_dft_1d(ring_polymer%fftw_plan_backward, ring_polymer%nbeads, &
    & ring_polymer%fftw_cmplx_in_n, ring_polymer%fftw_cmplx_out_n, fftw_backward, fftw_measure)
end subroutine pathint_normal_modes_init

subroutine pathint_normal_modes_pos_forward(ring_polymer)
  use io, only: io_err
  implicit none

  type(ring_polymer_type),  intent(inout) ::  ring_polymer

  ! local vars
  integer ::  icomp, iatom, ibead
  integer ::  natoms, nbeads, i, cn

  nbeads = ring_polymer%nbeads
  natoms = ring_polymer%beads(1)%cell%natoms  ! same for each cell..

  cn = nbeads/2 + 1

  ! work goes here..
  do icomp = 1, 3
    do iatom = 1, natoms

      !---------------------------------------------------------------------------
      ! From: Cao & Martyna, JCP 104, 2028 (1996), Appendix:
      !
      ! To transform the Cartesian positions to normal modes coordinates
      !
      ! 1. Set the real part of the Ith element of the a complex vector equal to
      !      the Cartesian position of the (I)th bead  (I=1, P).
      !
      ! 2. Perform a scaled backward fast Fourier transform on the vector.
      !
      ! 3. Set the value of the (1)st mode (the centroid position) equal to the
      !      real part of the (1)st element of the transformed vector.
      !
      ! 4. Set the value of the Pth mode equal to the real part of the 
      !      (P/2 + 1)st element of the transformed vector.
      !
      ! 5. Set the value of the (2I-2)th mode equal to the real part of the
      !      (I)th element of the transformed vector (I=2, P/2).
      !
      ! 6. Set the value of the (2I-1)th mode equal to the imaginary part of the
      !      (I)th element of the transformed vector (I=2, P/2).
      !-----------------------------------------------------------------------------

      do ibead = 1, nbeads
        ! get positions (cartesian/primitive variables) from md type
        ! rather than cell type (to allow movement over cell boundary)

        ! step 1:
        ! no imaginary part, so just use real array
        ring_polymer%fftw_real_n(ibead) = ring_polymer%beads(ibead)%md%pos(icomp,iatom)
      end do ! ibead

      ! with correct position, can now do some magic..
      ! based on size of "the a complex vector" (vector called 'a'), and the fact that it only contains
      ! real data.. and we want complex results.. do r2c fft..

      ! don't need to create and destroy plan every time (since everything the same.. but let's just be safe..)
      ! step 2:
      ! zero cvec for safety..
      ring_polymer%fftw_cmplx_out_no2p1(:) = cmplx(0.0_dp, 0.0_dp, dp)

      ! do transform
      call dfftw_execute(ring_polymer%fftw_plan_forward)
      ! and scale factor..
      ring_polymer%fftw_cmplx_out_no2p1(:) = ring_polymer%fftw_cmplx_out_no2p1(:)/real(nbeads,dp)

      ! cvec should now contain everything we need..
      ! copy everything back into real array
      ! step 3 & 4:
      ring_polymer%fftw_real_n(1) = real(ring_polymer%fftw_cmplx_out_no2p1(1),dp)
      ring_polymer%fftw_real_n(nbeads) = real(ring_polymer%fftw_cmplx_out_no2p1(cn),dp)   ! cn = nbeads/2 + 1 

      ! step 5 & 6:
      do i = 2, (nbeads/2) 
        ! step 5:
        ring_polymer%fftw_real_n(2*i-2) = real(ring_polymer%fftw_cmplx_out_no2p1(i),dp)
        ! step 6:
        ring_polymer%fftw_real_n(2*i-1) = real(aimag(ring_polymer%fftw_cmplx_out_no2p1(i)), dp)
      end do

      ! now copy back into ring polymer type..
      do ibead = 1, nbeads
        ring_polymer%beads(ibead)%md%pos(icomp,iatom) = ring_polymer%fftw_real_n(ibead)
      end do

    ! should be done..
    end do !iatom
  end do !icomp
end subroutine pathint_normal_modes_pos_forward

subroutine pathint_normal_modes_pos_backward(ring_polymer)
  use io, only: io_err
  implicit none

  type(ring_polymer_type),  intent(inout) ::  ring_polymer

  ! local vars
  integer ::  icomp, iatom, ibead
  integer ::  ri3, ri4, ci
  integer ::  natoms, nbeads, i

  nbeads = ring_polymer%nbeads
  natoms = ring_polymer%beads(1)%cell%natoms  ! same for each cell..

  ! work goes here..
  do icomp = 1, 3
    do iatom = 1, natoms

      do ibead = 1, nbeads
        ! get normal modes from md type (note cell positions *must* be primitive/cartesian variables)
        ! copy to real array first..
        ring_polymer%fftw_real_n(ibead) = ring_polymer%beads(ibead)%md%pos(icomp,iatom)
      end do ! ibead

      !-----------------------------------------------------------------------------
      ! From: Cao & Martyna, JCP 104, 2028 (1996), Appendix:
      !
      ! To transform the normal modes coordinates to Cartesian positions
      !
      ! 1. Set the real part of the (1)st element of complex vector equal to the value of the
      !      (1)st mode (the centroid position).
      !
      ! 2. Set the real part of the (P/2+1)st element of the complex vector equal to the value
      !      of the Pth mode.
      !
      ! 3. Set the real part of the (I)th and the (P-I+2)th elements of the complex vector
      !      equal to the value of the (2I-2)th mode (I=1, P/2).
      !
      ! 4. Set the imaginary part of the (I)th and the (P-I+2)th elements of the complex vector
      !      equal to the plus and minus the value of the (2I-1)th mode (I=1, P/2), respectively.
      !
      ! 5. Perform an unscaled forward fast Fourier transform on the vector.
      !
      ! 6. The Cartesian position of the (I)th bead is equal to the real part of the Ith element
      !      of the transformed vector (I=1, P).
      !-----------------------------------------------------------------------------

      ! zero everything for safety:
      ring_polymer%fftw_cmplx_in_n(:) = cmplx(0.0_dp, 0.0_dp, dp)

      ! step 1 & 2:
      ring_polymer%fftw_cmplx_in_n(1) = cmplx(ring_polymer%fftw_real_n(1), 0.0_dp, dp)
      ring_polymer%fftw_cmplx_in_n((nbeads/2) + 1) = cmplx(ring_polymer%fftw_real_n(nbeads), 0.0_dp, dp)

      ! step 3 & 4:
      ! when i = 1, ri3 == 0, which is invalid.. this leads to ci == nbeads + 1, which is also invalid..
      do i = 2, (nbeads/2)
        ri3 = 2*i-2       ! real array index, step 3
        ci = nbeads-i+2   ! cmplx array index, steps 3 and 4.. other index is i
        ri4 = 2*i-1       ! real array index, step 4

        ring_polymer%fftw_cmplx_in_n(i) = cmplx(ring_polymer%fftw_real_n(ri3),   ring_polymer%fftw_real_n(ri4), dp)
        ring_polymer%fftw_cmplx_in_n(ci) = cmplx(ring_polymer%fftw_real_n(ri3), -ring_polymer%fftw_real_n(ri4), dp)
      end do

      ! step 5:
      ! now do unscaled fft (go from c2r, which is backwards, not forwards...)
      ! have a complex array which is size nbeads rather than nbeads/2 + 1 - need to use complex to complex
      ! zero for safety
      ring_polymer%fftw_cmplx_out_n(:) = cmplx(0.0_dp, 0.0_dp, dp)
      call dfftw_execute(ring_polymer%fftw_plan_backward)

      ! step 6:
      ! now copy everything from cvec_out into ring polymer type..
      do ibead = 1, nbeads
        ring_polymer%beads(ibead)%md%pos(icomp,iatom) = real(ring_polymer%fftw_cmplx_out_n(ibead), dp)
      end do

    ! should be done..
    end do !iatom
  end do !icomp
end subroutine pathint_normal_modes_pos_backward

subroutine pathint_normal_modes_force_forward(ring_polymer)
  use io, only: io_err
  implicit none

  type(ring_polymer_type),  intent(inout) ::  ring_polymer

  ! local vars
  integer ::  icomp, iatom, ibead
  integer ::  natoms, nbeads, i, cn

  nbeads = ring_polymer%nbeads
  natoms = ring_polymer%beads(1)%cell%natoms  ! same for each cell..
  cn = nbeads/2 + 1                           ! size of complex work array

  ! work goes here..
  do icomp = 1, 3
    do iatom = 1, natoms

      !-----------------------------------------------------------------------------
      ! From: Cao & Martyna, JCP 104, 2028 (1996), Appendix:
      !
      ! To get the forces on the normal modes:
      !
      ! 1. Set the real part of a complex vector equal to the Cartesian force (I=1, P).
      !
      ! 2. Set the imaginary part of the complex vector equal to zero (I=1, P).
      !
      ! 3. Perform an unscaled backward fast Fourier transform on the vector.
      !
      ! 4. The force on the (1)st mode is the real part of the 1st element of the
      !     transformed vector.
      !
      ! 5. The force on the (P)th mode is the real part of the (P/2+1)th element of
      !     the transformed vector.
      !
      ! 6. The force on the (1)st mode (the centroid) is the real part of the 1st
      !     element of the transformed vector.
      !
      ! 7. The force on the (2I-2)th mode is twice the real part of the Ith element
      !     of the transformed vector (I=2, P/2).
      !
      ! 8. The force on the (2I-1)th mode is twice the imaginary part of the Ith
      !     element of the transformed vector (I=2, P/2).
      !-----------------------------------------------------------------------------

      do ibead = 1, nbeads
        ! get forces (in cartesian/primitive variable space) from md type

        ! step 1 & 2:
        ! no need to use complex array, since imaginary part is zero
        ring_polymer%fftw_real_n(ibead) = ring_polymer%beads(ibead)%md%force(icomp,iatom)
      end do ! ibead

      ! since complex array does not contain any imaginary components, can just use real array
      ! and receive N/2 + 1 complex array..
      ! zero cvec for safety..
      ring_polymer%fftw_cmplx_out_no2p1(:) = cmplx(0.0_dp, 0.0_dp, dp)

      ! step 3:
      ! should be *unscaled* transform..
      ! don't need to create and destroy plan every time (since everything the same.. but let's just be safe..)
      call dfftw_execute(ring_polymer%fftw_plan_forward)

      ! cvec should now contain everything we need..
      ! step 4 & 5:
      ring_polymer%fftw_real_n(1) = real(ring_polymer%fftw_cmplx_out_no2p1(1),dp)
      ring_polymer%fftw_real_n(nbeads) = real(ring_polymer%fftw_cmplx_out_no2p1(cn),dp)   ! cn = nbeads/2 + 1 
      ! step 6 is redundant (same as step 4)

      do i = 2, (nbeads/2) 
        ! step 7:
        ring_polymer%fftw_real_n(2*i-2) = 2.0_dp*real(ring_polymer%fftw_cmplx_out_no2p1(i),dp)
        ! step 8:
        ring_polymer%fftw_real_n(2*i-1) = 2.0_dp*real(aimag(ring_polymer%fftw_cmplx_out_no2p1(i)), dp)
      end do

      ! now copy back into ring polymer type..
      do ibead = 1, nbeads
        ring_polymer%beads(ibead)%md%force(icomp,iatom) = ring_polymer%fftw_real_n(ibead)
      end do

    ! should be done..
    end do !iatom
  end do !icomp
end subroutine pathint_normal_modes_force_forward

subroutine pathint_normal_modes_finalize(ring_polymer)
  use io, only: io_err
  type(ring_polymer_type),  intent(inout) ::  ring_polymer
  integer :: istat

  deallocate(ring_polymer%fftw_real_n, stat=istat)
  if (istat .ne. 0) call io_err("pathint_normal_modes_finalize: Could not deallocate fftw_real_n")

  deallocate(ring_polymer%fftw_cmplx_in_n, stat=istat)
  if (istat .ne. 0) call io_err("pathint_normal_modes_finalize: Could not deallocate fftw_cmplx_in_n")

  deallocate(ring_polymer%fftw_cmplx_out_n, stat=istat)
  if (istat .ne. 0) call io_err("pathint_normal_modes_finalize: Could not deallocate fftw_cmplx_out_n")

  deallocate(ring_polymer%fftw_cmplx_out_no2p1, stat=istat)
  if (istat .ne. 0) call io_err("pathint_normal_modes_finalize: Could not deallocate fftw_cmplx_out_no2p1")

  deallocate(ring_polymer%nm_eigenvalues, stat=istat)
  if (istat .ne. 0) call io_err("pathint_deallocate: Could not deallocate nm_eigenvalues array")

  call dfftw_destroy_plan(ring_polymer%fftw_plan_forward)
  call dfftw_destroy_plan(ring_polymer%fftw_plan_backward)
end subroutine pathint_normal_modes_finalize

subroutine do_pimd(ring_polymer, centroid_cell, centroid_md)
  use md,         only: md_velocity_verlet_p1, md_velocity_verlet_p2, md_sync_coords, md_get_forces, &
                      & md_zero_net_force, md_calc_kinetic_energy, md_write_output
  use constants,  only: units_atomic_to_natural
  use io,         only: stdout, io_close_file
  use checkpoint, only: continuation
  implicit none
  type(ring_polymer_type),  intent(inout) ::  ring_polymer
  type(cell_type),          intent(inout) ::  centroid_cell
  type(md_type),            intent(inout) ::  centroid_md
  ! local vars:
  integer ::  itime, ibead, ndof, start_time

  start_time = centroid_md%timestep + 1

  ! don't write out duplicate data..
  if (.not. continuation) then
    ! write iter 0 data
    call md_write_output(centroid_md, centroid_cell)

    if (ring_polymer%output_bead_trajectories) then
      do ibead = 1, ring_polymer%nbeads
        call md_write_output(ring_polymer%beads(ibead)%md, ring_polymer%beads(ibead)%cell)
      end do
    end if

    ! estimators (calculated in pathint_init)
    write(stdout,*) units_atomic_to_natural(energy=ring_polymer%primitive_energy_est), '<-PTE'
    write(stdout,*) units_atomic_to_natural(energy=ring_polymer%virial_energy_est), '<-VTE'
  end if


  ! real time loop..
  do itime = start_time, centroid_md%ntimesteps

    ! transform into normal modes for propagation
    if (ring_polymer%propagator_variables .eq. 'normalmodes') then
      call pathint_normal_modes_pos_forward(ring_polymer)
    end if

    ! change each bead's mass so that momenta are rescaled during propagation...
    call pathint_rescale_masses(ring_polymer, 'forward')

    ! first part of propagation
    do ibead = 1, ring_polymer%nbeads
      call md_velocity_verlet_p1(ring_polymer%beads(ibead)%md, ring_polymer%beads(ibead)%cell)
    end do

    ! transform back to primitive modes if necessary.. (for force calculation)
    if (ring_polymer%propagator_variables .eq. 'normalmodes') then
      call pathint_normal_modes_pos_backward(ring_polymer)
    end if

    do ibead = 1, ring_polymer%nbeads
      ! sync coords
      call md_sync_coords(ring_polymer%beads(ibead)%md, ring_polymer%beads(ibead)%cell)

      ! get forces
      call md_get_forces(ring_polymer%beads(ibead)%md, ring_polymer%beads(ibead)%cell)
    end do

    ! transform back to NM variables and get NM forces..
    if (ring_polymer%propagator_variables .eq. 'normalmodes') then
      call pathint_normal_modes_pos_forward(ring_polymer)
      call pathint_normal_modes_force_forward(ring_polymer)
    end if

    ! bead-bead interactions (need to using correct variables here..)
    ! should always be true except for testing..
    if (ring_polymer%bead_interactions) then
      ! first rescale masses back to physical masses..
      call pathint_rescale_masses(ring_polymer, 'backward')
      call pathint_bead_neighbour_interactions(ring_polymer)
      ! rescale back so we can finish propagation..
      call pathint_rescale_masses(ring_polymer, 'forward')
    end if

    if ((centroid_md%ensemble .eq. 'nvt') .and. (centroid_md%thermostat .eq. 'langevin')) then
      ! might need to move this..
      ! Langevin on each bead (+ centroid, implicitly)
      call pathint_langevin_force(ring_polymer, centroid_cell, centroid_md)
    end if

    ! final part of propagation
    do ibead = 1, ring_polymer%nbeads
      ! zero net force for each replica..
      if (ring_polymer%beads(ibead)%md%fix_com) then
        call md_zero_net_force(ring_polymer%beads(ibead)%md, ring_polymer%beads(ibead)%cell)
      end if

      call md_velocity_verlet_p2(ring_polymer%beads(ibead)%md, ring_polymer%beads(ibead)%cell)
    end do

    ! finally rescale masses back to actual (physical) masses.. (for output)
    call pathint_rescale_masses(ring_polymer, 'backward')

    do ibead = 1, ring_polymer%nbeads

      ! calc KE for each bead (not centroid, this is averaged from each bead)
      call md_calc_kinetic_energy(ring_polymer%beads(ibead)%md, ring_polymer%beads(ibead)%cell)
      ndof = ring_polymer%beads(ibead)%md%ndof
      if (ndof .ne. 0) then
        ring_polymer%beads(ibead)%md%inst_temp = 2.0_dp*ring_polymer%beads(ibead)%md%kinetic_energy/real(ndof,dp)
      else
        ring_polymer%beads(ibead)%md%inst_temp = 0.0_dp ! because KE = 0, but avoids NaN
      end if

      ! copy current timestep back to each bead..
      ring_polymer%beads(ibead)%md%timestep = itime
    end do
    centroid_md%timestep = itime

    ! transform back to primitive modes if necessary.. (for calculation of estimators)
    if (ring_polymer%propagator_variables .eq. 'normalmodes') then
      call pathint_normal_modes_pos_backward(ring_polymer)
    end if

    ! write output
    if (mod(itime, centroid_md%output_interval) .eq. 0) then

      ! centroid properties must be calculated before estimators..
      call pathint_calc_centroid(ring_polymer, centroid_cell, centroid_md)

      ! now calc estimators
      call pathint_calc_estimators(ring_polymer, centroid_cell, centroid_md)
      write(stdout,*) units_atomic_to_natural(energy=ring_polymer%primitive_energy_est), '<-PTE'
      write(stdout,*) units_atomic_to_natural(energy=ring_polymer%virial_energy_est), '<-VTE'

      call pathint_calc_rgy(ring_polymer, centroid_cell)

      call md_write_output(centroid_md, centroid_cell)
      if (ring_polymer%output_bead_trajectories) then
        do ibead = 1, ring_polymer%nbeads
          call md_write_output(ring_polymer%beads(ibead)%md, ring_polymer%beads(ibead)%cell)
        end do
      end if
    end if

    ! write checkpoint if necessary..
    if (mod(itime, centroid_md%checkpoint_interval) .eq. 0) then
      call pathint_write_checkpoint(ring_polymer, centroid_cell, centroid_md)
    end if
  end do

  ! write final checkpoint data if finished
  call pathint_write_checkpoint(ring_polymer, centroid_cell, centroid_md)
end subroutine do_pimd

subroutine pathint_rescale_masses(ring_polymer, direction)
  use constants,  only: two_pi
  use io,         only: io_err
  implicit none
  type(ring_polymer_type),  intent(inout) ::  ring_polymer
  character(len=*),         intent(in)    ::  direction     ! 'forward' or 'backward'

  real(kind=dp) ::  scale_factor
  integer       ::  ibead, iatom

  if ((direction .ne. 'forward') .and. (direction .ne. 'backward')) then
    call io_err("pathint_rescale_masses: Direction must be forward or backward")
  end if

  select case (ring_polymer%propagator_variables)
    case ('primitive')
      ! m' = mP/((2*pi*hbar)**2)  =>  m' = mP/((2*pi)**2), (hbar = 1 in atomic units)
      scale_factor = real(ring_polymer%nbeads,dp)/(two_pi**2)
      if (direction .eq. 'backward') scale_factor = 1.0_dp/scale_factor

      do ibead = 1, ring_polymer%nbeads
        do iatom = 1, ring_polymer%beads(ibead)%cell%natoms ! should be same for each bead..
          ring_polymer%beads(ibead)%cell%atom_mass(iatom) = ring_polymer%beads(ibead)%cell%atom_mass(iatom)*scale_factor
        end do
      end do

    case ('normalmodes')
      ! m_1' = m
      ! m_k' = m*lambda_k*(gamma**2)
      do ibead = 2, ring_polymer%nbeads
        scale_factor = ring_polymer%nm_eigenvalues(ibead)*(ring_polymer%adiabaticity**2)
        if (direction .eq. 'backward') scale_factor = 1.0_dp/scale_factor

        do iatom = 1, ring_polymer%beads(ibead)%cell%natoms ! should be same for each bead..
          ring_polymer%beads(ibead)%cell%atom_mass(iatom) = ring_polymer%beads(ibead)%cell%atom_mass(iatom)*scale_factor
        end do
      end do

    case default
      call io_err("pathint_rescale_masses: Unsupported propagator variable transformation")
  end select
end subroutine pathint_rescale_masses

subroutine pathint_bead_neighbour_interactions(ring_polymer)
  use constants,  only: two_pi
  use cell,       only: cell_min_img_vec
  implicit none
  type(ring_polymer_type),  intent(inout) ::  ring_polymer

  real(kind=dp),  dimension(3)  ::  current_bead_pos, prev_bead_pos, next_bead_pos
  real(kind=dp),  dimension(3)  ::  bead_force
  real(kind=dp),  dimension(3)  ::  current_next_pos, current_prev_pos  ! distances between beads
  real(kind=dp),  dimension(3)  ::  pos
  real(kind=dp) ::  mass
  integer :: iatom, ibead, prev_bead, next_bead

  ! assume valid propagator_variables here...

  do ibead = 1, ring_polymer%nbeads
    ! don't actually need these for NM..
    prev_bead = ring_polymer%beads(ibead)%prev_bead_index
    next_bead = ring_polymer%beads(ibead)%next_bead_index

    do iatom = 1, ring_polymer%beads(ibead)%cell%natoms
      bead_force(:) = ring_polymer%beads(ibead)%md%force(:,iatom)
      mass = ring_polymer%beads(ibead)%cell%atom_mass(iatom)

      select case (ring_polymer%propagator_variables)
        case ('primitive')
          current_bead_pos(:) = ring_polymer%beads(ibead)%cell%atom_cart_pos(:,iatom)
          prev_bead_pos(:) = ring_polymer%beads(prev_bead)%cell%atom_cart_pos(:,iatom)
          next_bead_pos(:) = ring_polymer%beads(next_bead)%cell%atom_cart_pos(:,iatom)


          ! need:
          ! 2*x_k - x_(k+1) - x_(k-1)
          ! => x_k - x_(k+1)  +  x_k - x_(k-1)
          ! distance between next and current beads
          call cell_min_img_vec(ring_polymer%beads(ibead)%cell, &
                            & current_bead_pos(:), next_bead_pos(:), current_next_pos(:))

          ! distance between current and previous beads
          call cell_min_img_vec(ring_polymer%beads(ibead)%cell, &
                            & current_bead_pos(:), prev_bead_pos(:), current_prev_pos(:))

          pos(:) = current_next_pos(:)  + current_prev_pos(:)
        case ('normalmodes')
          ! NM transformation means we don't need to consider nearest neighbour interactions here..
          ! but do need to use NM variable which is inside the md type, not the cell!
          pos(:) = ring_polymer%beads(ibead)%md%pos(:,iatom)

          ! rescale mass as necessary for this part - have physical masses:
          mass = mass*ring_polymer%nm_eigenvalues(ibead)
      end select

      ! raw force/number beads
      bead_force(:) = bead_force(:)/real(ring_polymer%nbeads,dp)

      ! now add nearest neighbour contributions
      bead_force(:) = bead_force(:) - mass*(ring_polymer%chain_frequency**2)*pos(:)

      ! copy back to md type
      ring_polymer%beads(ibead)%md%force(:,iatom) = bead_force(:)
    end do
  end do
end subroutine pathint_bead_neighbour_interactions

subroutine pathint_langevin_force(ring_polymer, centroid_cell, centroid_md)
  use md,     only: md_langevin_force
  use io,     only: io_err
  use algor,  only: algor_gauss_rand
  implicit none
  type(ring_polymer_type),  intent(inout) ::  ring_polymer
  type(cell_type),          intent(in)    ::  centroid_cell
  type(md_type),            intent(in)    ::  centroid_md
  ! local vars:
!  real(kind=dp),  dimension(3)  ::  centroid_langevin_force
!  real(kind=dp),  dimension(3)  ::  net_beads_langevin_force
!  real(kind=dp) ::  langevin_gamma, sqrt_mass
  integer ::  ibead !iatom, icomp, ibead, istat

  ! should be able to thermostat each bead separately in NM trans.. - probably not as good for primitive vars
  if (ring_polymer%nm_thermostat_centroid) then
    call md_langevin_force(ring_polymer%beads(1)%md, ring_polymer%beads(1)%cell)
  end if

  do ibead = 2, ring_polymer%nbeads
    call md_langevin_force(ring_polymer%beads(ibead)%md, ring_polymer%beads(ibead)%cell)
  end do

  ! below is wrong, but good basis for fix...

!  ! allocate tmp beads langevin force array on first pass...
!  if (.not. allocated(ring_polymer%beads_langevin_force)) then
!    allocate(ring_polymer%beads_langevin_force(3, ring_polymer%nbeads), stat=istat)
!    if (istat .ne. 0) call io_err("pathint_langevin_force: Could not allocate beads_langevin_force array")
!  end if
!
!  langevin_gamma = 1.0_dp/centroid_md%langevin_time
!  do iatom = 1, centroid_cell%natoms
!
!    ! Langevin contribution to centroid force (assume no external potential for now)..
!    ! Don't actually update centroid force etc since we do this when we average over beads..
!    sqrt_mass = sqrt(centroid_cell%atom_mass(iatom))
!    centroid_langevin_force(:) = 0.0_dp
!    do icomp = 1, 3
!      centroid_langevin_force(icomp) = -langevin_gamma*centroid_md%momentum(icomp,iatom) &
!        & + sqrt_mass*algor_gauss_rand(mean=0.0_dp, stddev=centroid_md%langevin_noise_stddev)
!    end do
!
!    ! divide centroid langevin force between each bead
!    centroid_langevin_force(:) = centroid_langevin_force(:)/real(ring_polymer%nbeads,dp)
!
!    ! now get langevin force for each of the beads
!    ! can't use md_langevin_force as that includes force due to potential...
!    net_beads_langevin_force(:) = 0.0_dp
!    do ibead = 1, ring_polymer%nbeads
!
!      ! keep things general..
!      langevin_gamma = 1.0_dp/ring_polymer%beads(ibead)%md%langevin_time
!      sqrt_mass = sqrt(ring_polymer%beads(ibead)%cell%atom_mass(iatom))
!
!      ! langevin force on each bead..
!      do icomp = 1,3
!        ring_polymer%beads_langevin_force(icomp, ibead) = &
!          & -langevin_gamma*ring_polymer%beads(ibead)%md%momentum(icomp,iatom) &
!          & + sqrt_mass*algor_gauss_rand(mean=0.0_dp, &
!          & stddev=ring_polymer%beads(ibead)%md%langevin_noise_stddev)
!      end do
!
!      net_beads_langevin_force(:) = net_beads_langevin_force(:) + ring_polymer%beads_langevin_force(:,ibead)
!    end do
!
!    ! divide net force among all beads
!    net_beads_langevin_force(:) = net_beads_langevin_force(:)/real(ring_polymer%nbeads,dp)
!
!    do ibead = 1, ring_polymer%nbeads
!    ! zero net bead langevin force + add centroid langevin force..
!      ring_polymer%beads_langevin_force(:,ibead) = ring_polymer%beads_langevin_force(:,ibead) &
!        & - net_beads_langevin_force(:) + centroid_langevin_force(:)
!
!      ! copy back to md%force for each bead
!      ring_polymer%beads(ibead)%md%force(:,iatom) = &
!        & ring_polymer%beads(ibead)%md%force(:,iatom) + ring_polymer%beads_langevin_force(:,ibead)
!
!      ! correct for v(t + dt/2)
!      ring_polymer%beads(ibead)%md%force(:,iatom) = &
!        & ring_polymer%beads(ibead)%md%force(:,iatom)*ring_polymer%beads(ibead)%md%langevin_force_correction
!    end do
!  end do
end subroutine pathint_langevin_force

subroutine pathint_calc_centroid(ring_polymer, centroid_cell, centroid_md)
  use md, only: md_sync_coords
  use io, only: io_err
  implicit none
  type(ring_polymer_type),  intent(in)    ::  ring_polymer
  type(cell_type),          intent(inout) ::  centroid_cell
  type(md_type),            intent(inout) ::  centroid_md
  real(kind=dp),  dimension(3)  ::  centroid_pos
  real(kind=dp),  dimension(3)  ::  centroid_momentum
  real(kind=dp),  dimension(3)  ::  centroid_force
  integer ::  iatom, ibead

  ! if NM trans, the positions must be in primitive.., everything else in NMs (force and momentum)

    do iatom = 1, centroid_cell%natoms

      centroid_pos(:) = 0.0_dp
      centroid_momentum(:) = 0.0_dp
      centroid_force(:) = 0.0_dp

      do ibead = 1, ring_polymer%nbeads
        ! use md%pos, rather than cell%atom_cart_pos as we want to track movement beyond cell boundaries
        centroid_pos(:) = centroid_pos(:) + ring_polymer%beads(ibead)%md%pos(:,iatom)

        centroid_momentum(:) = centroid_momentum(:) + ring_polymer%beads(ibead)%md%momentum(:,iatom)
        centroid_force(:) = centroid_force(:) + ring_polymer%beads(ibead)%md%force(:,iatom)
      end do

      centroid_pos(:) = centroid_pos(:)/real(ring_polymer%nbeads,dp)
      centroid_momentum(:) = centroid_momentum(:)/real(ring_polymer%nbeads,dp)
      centroid_force(:) = centroid_force(:)/real(ring_polymer%nbeads,dp)

      ! use this instead..
      if (ring_polymer%propagator_variables .eq. 'normalmodes') then
        centroid_momentum(:) = ring_polymer%beads(1)%md%momentum(:,iatom)
        centroid_force(:) = ring_polymer%beads(1)%md%force(:,iatom)
      end if

      ! put into centroid_md%pos again, rather than centroid_cell%atom_cart_pos
      centroid_md%pos(:,iatom) = centroid_pos(:)
      centroid_md%momentum(:,iatom) = centroid_momentum(:)
      centroid_md%force(:,iatom) = centroid_force(:)
    end do

    ! KE and PE
    centroid_md%kinetic_energy = 0.0_dp
    centroid_md%potential_energy = 0.0_dp

    do ibead = 1, ring_polymer%nbeads
      centroid_md%potential_energy = centroid_md%potential_energy + ring_polymer%beads(ibead)%md%potential_energy
      centroid_md%kinetic_energy = centroid_md%kinetic_energy + ring_polymer%beads(ibead)%md%kinetic_energy
    end do
    centroid_md%potential_energy = centroid_md%potential_energy/real(ring_polymer%nbeads,dp)
    centroid_md%kinetic_energy = centroid_md%kinetic_energy/real(ring_polymer%nbeads,dp)

    ! temperature
    if (centroid_md%ndof .ne. 0) then
      centroid_md%inst_temp = 2.0_dp*centroid_md%kinetic_energy/real(centroid_md%ndof,dp)
    else
      centroid_md%inst_temp = 0.0_dp ! always, because KE = 0, but avoids NaN
    end if

    ! resync the centroid_cell positions based on centroid_md%pos, so centroid_cell positions
    ! are mapped back to unit cell
    call md_sync_coords(centroid_md, centroid_cell)

    ! things can be simplier if normal modes.. easier to flip back to primitive.. but only position makes sense...
!    case ('normalmodes')
!      ! position, force and momentum given by bead index 1 (zeroth Fourier mode)
!      do iatom = 1, centroid_cell%natoms
!        centroid_md%pos(:,iatom) = ring_polymer%beads(1)%md%pos(:,iatom)
!        centroid_md%momentum(:,iatom) = ring_polymer%beads(1)%md%momentum(:,iatom)
!        centroid_md%force(:,iatom) = ring_polymer%beads(1)%md%force(:,iatom)
!      end do
!
!      centroid_md%potential_energy = ring_polymer%beads(1)%md%potential_energy
!      centroid_md%kinetic_energy = ring_polymer%beads(1)%md%kinetic_energy
!      centroid_md%inst_temp = ring_polymer%beads(1)%md%inst_temp
!    case default
!      call io_err("pathint_calc_centroid: Unsupported propagator variable transformation")
!  end select

end subroutine pathint_calc_centroid

subroutine pathint_calc_estimators(ring_polymer, centroid_cell, centroid_md)
  use cell, only: cell_min_img_vec
  implicit none
  type(ring_polymer_type),  intent(inout) ::  ring_polymer
  type(cell_type),          intent(in)    ::  centroid_cell
  type(md_type),            intent(in)    ::  centroid_md
  ! local vars
  real(kind=dp),  dimension(3)  ::  curr_next, curr_centroid, rik, rikp1, ric
  real(kind=dp) ::  ext_pot
  integer ::  iatom, ibead, next_bead

  ! energy estimators
  ring_polymer%primitive_energy_est = 0.0_dp
  ring_polymer%virial_energy_est = 0.0_dp

  ext_pot = 0.0_dp
  do ibead = 1, ring_polymer%nbeads
    next_bead = ring_polymer%beads(ibead)%next_bead_index

    do iatom = 1, centroid_cell%natoms

      ric = centroid_cell%atom_cart_pos(:,iatom)
      rik = ring_polymer%beads(ibead)%cell%atom_cart_pos(:,iatom)
      rikp1 = ring_polymer%beads(next_bead)%cell%atom_cart_pos(:,iatom)

      ! calc:
      ! r_i_k - r_i_(k+1)
      ! r_i_k - r_i_centroid
      ! doesn't matter which cell we use - only need lattice vectors which should be the same..
      call cell_min_img_vec(centroid_cell, rik, rikp1, curr_next)
      call cell_min_img_vec(centroid_cell, rik, ric, curr_centroid)

      ! sum: -m_i*(r_i_k - r_i_(k+1))^2
      ring_polymer%primitive_energy_est = ring_polymer%primitive_energy_est &
      & -ring_polymer%beads(ibead)%cell%atom_mass(iatom)*dot_product(curr_next, curr_next)

      
      ! sum: (r_i_k - r_i_centroid) . dU/dr_i_k  =  sum: -(r_i_k - r_i_centroid) . F_i_k
      ring_polymer%virial_energy_est = ring_polymer%virial_energy_est &
      & - dot_product(curr_centroid, ring_polymer%beads(ibead)%md%force(:,iatom))
    end do

    ext_pot = ext_pot + ring_polymer%beads(ibead)%md%potential_energy
  end do

  ext_pot = ext_pot/real(ring_polymer%nbeads,dp)

  ! sum: -0.5*m_i*w^2*(r_i_k - r_i_(k+1))^2
  ring_polymer%primitive_energy_est = 0.5_dp*(ring_polymer%chain_frequency**2)*ring_polymer%primitive_energy_est

  ! sum: 0.5*(1/P)*(r_i_k - r_i_centroid) . dU/dr_i_k
  ring_polymer%virial_energy_est = 0.5_dp*ring_polymer%virial_energy_est/real(ring_polymer%nbeads,dp)

  ! add potential term:
  ring_polymer%primitive_energy_est = ring_polymer%primitive_energy_est + ext_pot
  ring_polymer%virial_energy_est = ring_polymer%virial_energy_est + ext_pot

  ! dN should be number of degrees of freedom (take into account fixed CoM)..
  ! primitive: dNPkT/2 + other terms
  ring_polymer%primitive_energy_est = 0.5_dp*real(centroid_md%ndof*ring_polymer%nbeads,dp)*centroid_md%temperature &
  & + ring_polymer%primitive_energy_est

  ! virial: dNkT/2 + other terms
  ring_polymer%virial_energy_est = 0.5_dp*real(centroid_md%ndof,dp)*centroid_md%temperature &
  & + ring_polymer%virial_energy_est
end subroutine pathint_calc_estimators

subroutine pathint_calc_rgy(ring_polymer, centroid_cell)
  use cell, only: cell_min_img_vec
  use io,   only: io_err, stdout
  implicit none
  type(ring_polymer_type),  intent(inout) ::  ring_polymer
  type(cell_type),          intent(in)    ::  centroid_cell
  ! local vars
!  real(kind=dp),  dimension(3,3)  ::  rgy_matrix
  real(kind=dp),  dimension(:,:), allocatable   ::  r_bead_centroid
  real(kind=dp) :: rgy
  integer ::  iatom, ibead, istat!, i, j

  ! slow to keep allocating this..
  allocate(r_bead_centroid(3,ring_polymer%nbeads), stat=istat)
  if (istat .ne. 0) call io_err('pathint_calc_rgy: Could not allocate r_bead_centroid array')

  do iatom = 1, centroid_cell%natoms
    rgy = 0.0_dp
    do ibead = 1, ring_polymer%nbeads
      call cell_min_img_vec(centroid_cell, ring_polymer%beads(ibead)%cell%atom_cart_pos(:,iatom), &
        & centroid_cell%atom_cart_pos(:,iatom), r_bead_centroid(:,ibead))

      rgy = rgy + dot_product(r_bead_centroid(:,ibead),r_bead_centroid(:,ibead))
    end do
    rgy = rgy/real(ring_polymer%nbeads,dp)
    rgy = sqrt(rgy)

    write(stdout, *) iatom, rgy, '<-- RGY'

!    rgy_matrix(:,:) = 0.0_dp
!
!    do i = 1, 3
!      do j = i, 3
!        do ibead = 1, ring_polymer%nbeads
!          rgy_matrix(i,j) = rgy_matrix(i,j) + r_bead_centroid(i,ibead)*r_bead_centroid(j,ibead)
!        end do
!        ! can output for each atom
!        rgy_matrix(i,j) = sqrt(rgy_matrix(i,j)/real(ring_polymer%nbeads,dp))
!        if (i .ne. j) rgy_matrix(j,i) = rgy_matrix(i,j)
!      end do
!    end do
  end do

  deallocate(r_bead_centroid, stat=istat)
  if (istat .ne. 0) call io_err('pathint_calc_rgy: Could not deallocate r_bead_centroid array')
end subroutine pathint_calc_rgy

subroutine pathint_deallocate(ring_polymer)
  use cell, only: cell_deallocate
  use md,   only: md_deallocate
  use io,   only: io_err
  implicit none
  type(ring_polymer_type),  intent(inout) ::  ring_polymer
  integer ::  ibead, istat

  ! clean up
  ! some of these optionally allocated
  if (ring_polymer%propagator_variables .eq. 'normalmodes') then
    call pathint_normal_modes_finalize(ring_polymer)
  end if

  if (allocated(ring_polymer%beads_langevin_force)) then
    deallocate(ring_polymer%beads_langevin_force, stat=istat)
    if (istat .ne. 0) call io_err("pathint_deallocate: Could not deallocate beads_langevin_force array")
  end if
  
  do ibead = 1, ring_polymer%nbeads
    call cell_deallocate(ring_polymer%beads(ibead)%cell)
    call md_deallocate(ring_polymer%beads(ibead)%md)
  end do

  deallocate(ring_polymer%beads, stat=istat)
  if (istat .ne. 0) call io_err("pathint_deallocate: Could not deallocate beads array")
end subroutine pathint_deallocate

subroutine pathint_continuation(ring_polymer, centroid_cell, centroid_md)
  use checkpoint, only: checkpoint_open, checkpoint_unit, checkpoint_close, continuation
  use md,         only: md_continuation, md_type_read_checkpoint
  use cell,       only: cell_read_checkpoint
  use io,         only: io_err
  implicit none
  type(ring_polymer_type),  intent(inout) ::  ring_polymer
  type(cell_type),          intent(inout) ::  centroid_cell
  type(md_type),            intent(inout) ::  centroid_md
  integer :: ibead
  ! continuation vars:
  integer :: tmp_nbeads
  logical :: tmp_thermostat_centroid
  logical :: tmp_output
  character(len=11) :: tmp_prop
  real(kind=dp) ::  tmp_gamma

  if (.not. continuation) call io_err("pathint_continuation: pathint_continuation called rather than pathint_init")

  ! this allows md params to change..
  call md_continuation(centroid_md, centroid_cell, .false.)

  ! might need to enfore NVT?..
  call pathint_set_defaults(ring_polymer)
  call pathint_read_input(ring_polymer)

  ! store input values:
  tmp_nbeads = ring_polymer%nbeads
  tmp_output = ring_polymer%output_bead_trajectories
  tmp_prop = ring_polymer%propagator_variables
  tmp_gamma = ring_polymer%adiabaticity
  tmp_thermostat_centroid = ring_polymer%nm_thermostat_centroid

  ! don't need to copy nm_eigenvalues
  read(checkpoint_unit)  ring_polymer%nm_thermostat_centroid
  read(checkpoint_unit)  ring_polymer%adiabaticity
  read(checkpoint_unit)  ring_polymer%propagator_variables
  read(checkpoint_unit)  ring_polymer%output_bead_trajectories
  read(checkpoint_unit)  ring_polymer%nbeads

  ! exit if any differ..
  if (tmp_nbeads .ne. ring_polymer%nbeads) then
    call io_err("pathint_continuation: pathint_num_beads changes on continuation")
  end if

  if (tmp_output .neqv. ring_polymer%output_bead_trajectories) then
    call io_err("pathint_continuation: pathint_bead_trajectories changes on continuation")
  end if

  if (tmp_prop .ne. ring_polymer%propagator_variables) then
    call io_err("pathint_continuation: pathint_propagator_variables changes on continuation")
  end if

  if (tmp_gamma .ne. ring_polymer%adiabaticity) then
    call io_err("pathint_continuation: pathint_adiabaticity changes on continuation")
  end if

  if (tmp_thermostat_centroid .neqv. ring_polymer%nm_thermostat_centroid) then
    call io_err("pathint_continuation: pathint_nm_thermostat_centroid changes on continuation")
  end if

  ! now echo pathint only params
  call pathint_echo_params(ring_polymer)

  read(checkpoint_unit)  ring_polymer%chain_frequency
  allocate(ring_polymer%beads(ring_polymer%nbeads))
  ! now copy bead contents
  do ibead = 1, ring_polymer%nbeads
    read(checkpoint_unit)  ring_polymer%beads(ibead)%prev_bead_index
    read(checkpoint_unit)  ring_polymer%beads(ibead)%next_bead_index
    call cell_read_checkpoint(ring_polymer%beads(ibead)%cell)
    call md_type_read_checkpoint(ring_polymer%beads(ibead)%md, ring_polymer%beads(ibead)%cell)
  end do

  ! need to copy over any params that changed in centroid md..
  ! just manually copy everything (..except things on per bead basis.. for now)..bit of a hack
  ! this will also increment each bead's timestep so we don't duplicate data..
  do ibead = 1, ring_polymer%nbeads
    ring_polymer%beads(ibead)%md%trajectory_units = centroid_md%trajectory_units
    ring_polymer%beads(ibead)%md%timestep = centroid_md%timestep
    ring_polymer%beads(ibead)%md%ntimesteps = centroid_md%ntimesteps
    ring_polymer%beads(ibead)%md%ndof = centroid_md%ndof
    ring_polymer%beads(ibead)%md%dt = centroid_md%dt
    ring_polymer%beads(ibead)%md%ensemble = centroid_md%ensemble
    ring_polymer%beads(ibead)%md%thermostat = centroid_md%thermostat
    ring_polymer%beads(ibead)%md%temperature = centroid_md%temperature
    ring_polymer%beads(ibead)%md%fix_com = centroid_md%fix_com
    ring_polymer%beads(ibead)%md%output_interval = centroid_md%output_interval
    ring_polymer%beads(ibead)%md%checkpoint_interval = centroid_md%checkpoint_interval
    ring_polymer%beads(ibead)%md%langevin_time = centroid_md%langevin_time
    ring_polymer%beads(ibead)%md%langevin_noise_stddev = centroid_md%langevin_noise_stddev
    ring_polymer%beads(ibead)%md%langevin_force_correction = centroid_md%langevin_force_correction
  end do

  read(checkpoint_unit)  ring_polymer%primitive_energy_est
  read(checkpoint_unit)  ring_polymer%virial_energy_est

  call checkpoint_close

  ! chain frequency: w_p = sqrt(P)/(hbar*beta) = sqrt(P)*kT/hbar
  ring_polymer%chain_frequency = sqrt(real(ring_polymer%nbeads, dp))*centroid_md%temperature

  if (ring_polymer%propagator_variables .eq. 'normalmodes') then
    call pathint_normal_modes_init(ring_polymer)
  end if
end subroutine pathint_continuation

subroutine pathint_write_checkpoint(ring_polymer, centroid_cell, centroid_md)
  use checkpoint, only: checkpoint_open, checkpoint_unit, checkpoint_close
  use cell,       only: cell_write_checkpoint 
  use md,         only: md_type_write_checkpoint
  use algor,      only: rand_seed
  implicit none
  type(ring_polymer_type),  intent(in)  ::  ring_polymer
  type(cell_type),          intent(in)  ::  centroid_cell
  type(md_type),            intent(in)  ::  centroid_md
  integer :: i, ibead

  call checkpoint_open('write')
  do i = 1, size(rand_seed)
    write(checkpoint_unit) rand_seed(i)
  end do
  call cell_write_checkpoint(centroid_cell)
  call md_type_write_checkpoint(centroid_md)

  ! don't need to copy nm_eigenvalues or beads_langevin_force
  write(checkpoint_unit)  ring_polymer%nm_thermostat_centroid
  write(checkpoint_unit)  ring_polymer%adiabaticity
  write(checkpoint_unit)  ring_polymer%propagator_variables
  write(checkpoint_unit)  ring_polymer%output_bead_trajectories
  write(checkpoint_unit)  ring_polymer%nbeads
  write(checkpoint_unit)  ring_polymer%chain_frequency

  ! now write bead contents
  do ibead = 1, ring_polymer%nbeads
    write(checkpoint_unit)  ring_polymer%beads(ibead)%prev_bead_index
    write(checkpoint_unit)  ring_polymer%beads(ibead)%next_bead_index
    call cell_write_checkpoint(ring_polymer%beads(ibead)%cell)
    call md_type_write_checkpoint(ring_polymer%beads(ibead)%md)
  end do

  write(checkpoint_unit)  ring_polymer%primitive_energy_est
  write(checkpoint_unit)  ring_polymer%virial_energy_est

  call checkpoint_close
end subroutine pathint_write_checkpoint

end module pathint
