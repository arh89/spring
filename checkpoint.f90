module checkpoint
  use io, only: max_line_len
  implicit none

  private

  character(len=*), parameter       ::  check_ext = '.check'

  ! public vars
  logical,            save, public  ::  continuation
  integer,            save, public  ::  checkpoint_unit
  character(len=5),   save, public  ::  check_status = 'none' ! 'read', 'write'

  ! public routines:
  public  ::  checkpoint_open
  public  ::  checkpoint_close
contains

subroutine checkpoint_open(file_action)
  use io, only: seedname, io_err, io_find_unit
  implicit none
  character(len=*),  intent(in)  ::  file_action
  ! local vars
  character(len=max_line_len+len(check_ext))  :: filename
  logical  ::  file_exists
  integer  ::  istat

  if (check_status .eq. file_action) return

  ! check_status differs from file_action.. need to reopen file..
  if (check_status .ne. 'none') then
    ! only close if we have actually opened a file..
    call checkpoint_close
  end if

  filename = trim(seedname)//trim(check_ext)

  call io_find_unit(checkpoint_unit)

  select case (file_action)
  case ('read')
    ! to read, we require file to exist
    inquire(file=trim(filename), exist=file_exists)
    if (.not. file_exists) call io_err("checkpoint_open: File; "//trim(filename)//" does not exist")

    open(unit=checkpoint_unit, file=trim(filename), status='old', action='read', &
      & position='rewind', access='stream', form='unformatted', iostat=istat)

  case ('write')
    ! don't care if file already exists - will create new file if it doesn't exist
    open(unit=checkpoint_unit, file=trim(filename), status='replace', action='write', &
      & position='rewind', access='stream', form='unformatted', iostat=istat)

  case default
    call io_err("checkpoint_open: file_action must be 'read' or 'write'")
  end select

  check_status = file_action

  if (istat .ne. 0) call io_err("checkpoint_open: Failed to open file; "//trim(filename))
end subroutine checkpoint_open

subroutine checkpoint_close
  use io, only: io_close_file
  implicit none

  call io_close_file(checkpoint_unit)
  check_status = 'none'
end subroutine checkpoint_close
end module checkpoint
