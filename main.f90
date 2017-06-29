program main
  use constants,    only: dp
  use checkpoint,   only: continuation
  use io,           only: io_init, io_input_finalize, io_err, io_input_get_single_value, &
                        & io_query_keyword
  use cell,         only: cell_type, cell_init, cell_deallocate
  use md,           only: md_type, md_init, md_continuation, do_md, md_deallocate
  use algor,        only: algor_init_rand
  use potential,    only: pot_init, pot_finalize
  use singlepoint,  only: singlepoint_init, do_singlepoint, singlepoint_deallocate
  use geometry,     only: geom_type, geom_init, do_geom_opt, geom_deallocate
  use pathint,      only: ring_polymer_type, pathint_init, pathint_continuation, do_pimd, &
                        & pathint_deallocate
  use phonon,       only: phonon_type, phonon_init, do_phonon_calculation, phonon_deallocate
  use pes_scan,     only: pes_scan_type, pes_scan_init, do_pes_scan
  implicit none

  type(cell_type)         ::  cell
  type(md_type)           ::  md
  type(geom_type)         ::  geom
  type(ring_polymer_type) ::  ring_polymer
  type(phonon_type)       ::  phonon
  type(pes_scan_type)     ::  pes_scan

  character(len=32) ::  task
  integer           ::  rand_seed
  logical           ::  input_errors


  call io_init('replace')

  ! was random_seed included in input file?
  rand_seed = io_query_keyword('random_seed')
  if (rand_seed .eq. 0) then
    ! no.. call without initial seed
    call algor_init_rand()
  else if (rand_seed .eq. 1) then
    ! set required = .true. because it should be in the hash table!
    call io_input_get_single_value('random_seed', rand_seed, required=.true.)
    call algor_init_rand(rand_seed)
  else
    call io_err('main: random_seed expects a single value')
  end if

  call io_input_get_single_value('task', task, 'singlepoint')

  ! optionally turn off io_input_finalize (won't quit if input parameter entered incorrectly)
  call io_input_get_single_value('input_errors', input_errors, .true.)

  call io_input_get_single_value('continuation', continuation, .false.)
  ! if continuation, ignore what we did with the random seed above.. will get old value..

  ! done with main program input keywords..
  call cell_init(cell)

  call pot_init(cell)

  select case(task)
    case('sp','spt','singlepoint')
      call singlepoint_init(cell)
      if (input_errors) call io_input_finalize
      call do_singlepoint(cell)
      call singlepoint_deallocate

    case('md','moleculardynamics')
      if (continuation) then
        call md_continuation(md, cell)
      else
        call md_init(md, cell)
      end if
      if (input_errors) call io_input_finalize
      call do_md(md, cell)
      call md_deallocate(md)

    case('pimd','pathint','pathintegrals')
      if (continuation) then
        call pathint_continuation(ring_polymer, cell, md)
      else
        call pathint_init(ring_polymer, cell, md)
      end if
      if (input_errors) call io_input_finalize
      call do_pimd(ring_polymer, cell, md)
      call pathint_deallocate(ring_polymer)

    case('geomopt','geometryoptimization')
      call geom_init(geom, cell)
      if (input_errors) call io_input_finalize
      call do_geom_opt(geom, cell)
      call geom_deallocate(geom)

    case('phonon')
      call phonon_init(phonon, cell)
      if (phonon%do_geom_opt) then
        call geom_init(geom, cell)
        if (input_errors) call io_input_finalize
        call do_geom_opt(geom, cell)
        call geom_deallocate(geom)
      else
        if (input_errors) call io_input_finalize
      end if
      call do_phonon_calculation(phonon, cell)
      call phonon_deallocate(phonon)

    case('pes_scan')
      call pes_scan_init(pes_scan, cell)
      if (input_errors) call io_input_finalize
      call do_pes_scan(pes_scan, cell)

    case default
      call io_err('main: Unrecognised task '//task)
  end select

  call pot_finalize

  call cell_deallocate(cell)

end program main
