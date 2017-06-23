!****m* io/io =================================================================!
! NAME                                                                         !
!   io                                                                         !
!------------------------------------------------------------------------------!
! DESCRIPTION                                                                  !
!   A collection of routines for file IO, including an implementation of       !
!   free-form input (which stores its data in a hash table) and error handling.!
!------------------------------------------------------------------------------!
! NOTES                                                                        !
!   io_input_get_single_value_real and io_str_to_real each use the module      !
!   level definition of a double precision data type. This is to break any     !
!   cyclic dependencies between this module and the constants module.          !
!------------------------------------------------------------------------------!
! AUTHOR                                                                       !
!   Aaron Hopkinson                                                            !
!****==========================================================================!
module io
  use iso_fortran_env,  only: output_unit
  use hash_tables,      only: hash_table, str_len
  implicit none

  private
                                               ! This must be the same as in consts...
  integer,          parameter                 ::  dp = selected_real_kind(15,300)
  integer,          parameter,        public  ::  max_line_len  = str_len

  ! public variables:
  character(len=max_line_len),  save, public  ::  seedname = ' '                      ! initialize to blank, so we can call io_err
  integer,                      save, public  ::  stdout = output_unit                ! ...without a seedname


  ! defines the user interface:
  character(len=*), parameter                 ::  input_file_extension = '.in'
  character(len=*), parameter                 ::  input_comment_marker = '#'
  character(len=*), parameter                 ::  input_begin_block_marker = 'begin'
  character(len=*), parameter                 ::  input_end_block_marker = 'end'
  character(len=*), parameter                 ::  error_file_extension = '.error'
  character(len=*), parameter                 ::  output_file_extension = '.out'

  ! internal variables:
  type(hash_table),             save          ::  input_table
  logical,                      save          ::  err_in_use = .false.                ! used for io_err recursion
  integer,                      save, public  ::  err_unit_num = output_unit          ! default is stdout


  interface io_input_get_single_value
    module procedure io_input_get_single_value_int
    module procedure io_input_get_single_value_real
    module procedure io_input_get_single_value_logical
    module procedure io_input_get_single_value_str
  end interface io_input_get_single_value

  ! public routines:
  public  ::  io_init
  public  ::  io_open_file
  public  ::  io_find_unit
  public  ::  io_close_file
  public  ::  io_err
  public  ::  io_read_input_file
  public  ::  io_input_finalize
  public  ::  io_query_keyword
  public  ::  io_input_get_data
  public  ::  io_input_get_single_value
  public  ::  io_input_hash_echo
  public  ::  io_str_get_num_tokens
  public  ::  io_str_get_token
  public  ::  io_str_to_lcase
  public  ::  io_str_to_ucase
  public  ::  io_str_to_int
  public  ::  io_str_to_real
  public  ::  io_str_to_logical

contains

!****s* io/io_init ============================================================!
! NAME                                                                         !
!   io_init (PUBLIC)                                                           !
!------------------------------------------------------------------------------!
! ARGUMENTS                                                                    !
!   character(len=*), optional, intent(in)  ::  stdout_action                  !
!------------------------------------------------------------------------------!
! DESCRIPTION                                                                  !
!   Sets up IO ready for use.                                                  !
!                                                                              !
!   Gets seedname from command line arguments and calls io_read_input_file,    !
!   which parses and stores the significant content of the input file:         !
!   seedname//input_file_extension into the default (module level) hash table. !
!                                                                              !
!   This routine also opens the main output file:                              !
!   seedname//output_file_extension and stores the unit number in the stdout   !
!   variable - unless stdout_action is 'stdout', in which case we use the OS   !
!   stdout.                                                                    !
!------------------------------------------------------------------------------!
! NOTES                                                                        !
!   Only the first command line argument is used at present.                   !
!                                                                              !
!   input_file_extension is trimmed away from the right hand side of the first !
!   command line argument (if present).                                        !
!==============================================================================!
subroutine io_init(stdout_action)
  implicit none

  character(len=*), optional, intent(in)  ::  stdout_action
  character(len=7)  ::  int_stdout_action
  integer           ::  extension_pos

  ! default stdout action
  int_stdout_action = 'append'

  if (present(stdout_action)) then
    if ((stdout_action .ne. 'stdout') .and. (stdout_action .ne. 'replace')  &
      & .and. (stdout_action .ne. 'append')) then
      write(*,*) 'io_init: Invalid stdout_action'
      stop
    end if
    int_stdout_action = stdout_action
  end if

  if (command_argument_count() .lt. 1) then
    ! give user help...
    write (*,*) 'io_init: No seedname'
    stop
  else
    ! only use first argument - if more, these are ignored (for now)
    call get_command_argument(1, seedname)
    seedname = adjustl(seedname)

    ! remove input file extension if included
    ! - go from the back to cut off only the final instance
    !   this way if this extension is part of the seedname, we can preserve it by typing entire filename
    extension_pos = index(seedname, input_file_extension, back=.true.)
    if (extension_pos .gt. 0) then
      seedname = seedname(1:extension_pos-1)
    end if

    ! read input file into default hash table
    call io_read_input_file(trim(seedname)//input_file_extension, input_table)

    ! open output file - only surpress this if flag set to true
    if (int_stdout_action .ne. 'stdout') then
      call io_open_file(trim(seedname)//output_file_extension, stdout, int_stdout_action)
    end if
  end if
end subroutine io_init
 

!****s* io/io_open_file =======================================================!
! NAME                                                                         !
!   io_open_file (PUBLIC)                                                      !
!------------------------------------------------------------------------------!
! ARGUMENTS                                                                    !
!   character(len=*),           intent(in)  ::  filename                       !
!   integer,                    intent(out) ::  unit_num                       !
!                                                                              !    
!   character(len=*), optional, intent(in)  ::  file_action                    !
!------------------------------------------------------------------------------!
! DESCRIPTION                                                                  !
!   Attempts to open a file for reading/writing.                               !
!   Do so by calling io_find_unit, which returns the first free unit number.   !
!                                                                              !
!   file_action is optional and can be either:                                 !
!     'read'    - file must exist and is marked as read only                   !
!     'replace' - file need not exist and is replaced if it does               !
!     'append'  - file need not exist, if it does, we append to existing file  !
!                                                                              !
!   In the case of file_action = 'replace' or 'append', if the file does not   !
!   exist then we create a new one. In both of these cases, we are able to     !
!   read from, or write to, the file, though this may involve repositioning    !
!   within the file.                                                           !
!                                                                              !
!   If file_action is not present, the default action is 'append'.             !
!------------------------------------------------------------------------------!
! NOTES                                                                        !
!   Routine must be recursive because we can call io_err, which can call this. !
!==============================================================================!
recursive subroutine io_open_file(filename, unit_num, file_action)
  implicit none

  character(len=*),           intent(in)  ::  filename
  integer,                    intent(out) ::  unit_num
  character(len=*), optional, intent(in)  ::  file_action

  character(len=7)    ::  int_action
  logical             ::  file_exists
  integer             ::  istat

  if (len_trim(filename) .lt. 1) call io_err("io_open_file: Blank filename")

  if (present(file_action)) then

    ! give error if invalid file action
    if  ((file_action .ne. 'read') .and. (file_action .ne. 'append') .and. (file_action .ne. 'replace')) then
      call io_err("io_open_file: Unsupported file action")
    end if

    ! then set internal file action:
    int_action = file_action
  else
    int_action = 'append' ! probably the safest
  end if

  ! search until we find suitable unit
  call io_find_unit(unit_num)

  ! check to see if file exists..
  inquire(file=trim(filename), exist=file_exists)

  select case (int_action)
    case ('read')
      ! to read, we require file to exist
      if (.not. file_exists) call io_err("io_open_file: File; "//trim(filename)//" does not exist")

      open(unit=unit_num, file=trim(filename), status='old', action='read', position='rewind', iostat=istat)

    case ('replace')
      ! don't care if file already exists - will create new file if it doesn't exist
      open(unit=unit_num, file=trim(filename), status='replace', action='readwrite', position='rewind', iostat=istat)

    case ('append')
      ! *probably* don't care if the file exists.. (in terms of flagging a warning to user)
      if (file_exists) then
        open(unit=unit_num, file=trim(filename), status='old', action='readwrite', position='append', iostat=istat)
      else
        open(unit=unit_num, file=trim(filename), status='new', action='readwrite', position='rewind', iostat=istat)
      end if

    case default
      ! should be impossible
      call io_err("io_open_file: Invalid internal file action")
  end select

  if (istat .ne. 0) call io_err("io_open_file: Failed to open file; "//trim(filename))
end subroutine io_open_file


!****s* io/io_find_unit =======================================================!
! NAME                                                                         !
!   io_find_unit (PUBLIC)                                                      !
!------------------------------------------------------------------------------!
! ARGUMENTS                                                                    !
!   integer,  intent(out) ::  unit_num                                         !
!------------------------------------------------------------------------------!
! DESCRIPTION                                                                  !
!   Attempts to find the first free unit for file IO. Does so by scanning      !
!   through a list of possible units (up to max_unit) to find one that both    !
!   exists and is not open.                                                    !
!------------------------------------------------------------------------------!
! NOTES                                                                        !
!   Subroutine could equally be a function.                                    !
!                                                                              !
!   We could consider storing a list of previously closed units and selecting  !
!   from this.                                                                 !
!                                                                              !
!   Routine must be recursive because we can call io_err, which can call       !
!   io_open_file, which calls this.                                            !
!==============================================================================!
recursive subroutine io_find_unit(unit_num)
  implicit none

  integer,  intent(out) ::  unit_num

  integer,  parameter   ::  min_unit = 10, max_unit = 99 
  integer               ::  iunit
  logical               ::  unit_exists, unit_open, success

  success  = .false.

  do iunit = min_unit, max_unit
    inquire(unit=iunit, exist=unit_exists, opened=unit_open)

    ! want unit to exist and not be in use
    if (unit_exists .and. (.not. unit_open)) then
      success = .true.
      unit_num = iunit      ! only set once we know it is safe
      exit
    end if
  end do

  if (.not. success) call io_err("io_find_unit: Could not open any file units") 
end subroutine io_find_unit


!****s* io/io_close_file ======================================================!
! NAME                                                                         !
!   io_close_file (PUBLIC)                                                     !
!------------------------------------------------------------------------------!
! ARGUMENTS                                                                    !
!   integer,  intent(in)  ::  unit_num                                         !
!------------------------------------------------------------------------------!
! DESCRIPTION                                                                  !
!   Closes the opened file associated with unit_num.                           !
!------------------------------------------------------------------------------!
! NOTES                                                                        !
!   Routine must be recursive because we can call io_err, which can call this. !
!==============================================================================!
recursive subroutine io_close_file(unit_num)
  implicit none

  integer,  intent(in)  ::  unit_num

  logical ::  unit_exists, unit_open
  integer ::  istat

  ! check that both the unit exists and is opened:
  ! no point in giving error if unit doesn't exist/isn't open --
  ! actually F95 standard says closing a unit that doesn't exist/isn't connected to a file is safe... 
  inquire(unit=unit_num, exist=unit_exists, opened=unit_open)
  if (unit_exists .and. unit_open) then
    close(unit_num, iostat=istat)
    if (istat .ne. 0) call io_err("io_close file: Failed to close unit")
  end if
end subroutine io_close_file


!****s* io/io_err =============================================================!
! NAME                                                                         !
!   io_err (PUBLIC)                                                            !
!------------------------------------------------------------------------------!
! ARGUMENTS                                                                    !
!   character(len=*), intent(in)  ::  string                                   !
!------------------------------------------------------------------------------!
! DESCRIPTION                                                                  !
!   Attempts to write an error to file: seedname//error_file_extension         !
!   if a seedname is present - if not, we write error to standard out, then    !
!   we abort.                                                                  !
!                                                                              !
!   Routine is recursive to deal with errors that arise from io_open_file,     !
!   io_close_file, or their children. This is a rare case but it allows us to  !
!   attempt to preserve errors in a relatively straightforward way.            !
!==============================================================================!
recursive subroutine io_err(string)
  use iso_fortran_env,  only: output_unit
  implicit none

  character(len=*), intent(in)  ::  string

  ! if first time calling this:
  if (.not. err_in_use) then
    call io_err_initialize              ! open error file if we can

    ! err_unit_num could be a file, or stdout
    write(err_unit_num, *)  string
    
    ! if we opened a file, then close it
    ! (doesn't actually matter if we close stdout, since this won't give an error and we're finished anyway)
    if (err_unit_num .ne. output_unit) call io_close_file(err_unit_num)
    stop
  else
    ! any further calls don't cause a stop - only the 'root' call
    ! err_unit_num could be a file, or stdout - although it might be safer if we always write to stdout here...
    write(err_unit_num, *) string
  end if
end subroutine io_err


!****s* io/io_read_input_file =================================================!
! NAME                                                                         !
!   io_read_input_file (PUBLIC)                                                !
!------------------------------------------------------------------------------!
! ARGUMENTS                                                                    !
!   character(len=*), intent(in)    ::  filename                               !
!   type(hash_table), intent(inout) ::  table                                  !
!------------------------------------------------------------------------------!
! DESCRIPTION                                                                  !
!   Opens, parses and stores the content of the input file into a hash table.  !
!   Allows for free-form input.                                                !
!------------------------------------------------------------------------------!
! NOTES                                                                        !
!   Input table must be uninitialized.                                         !
!                                                                              !
!   We discard the end block line without checking that the keywords match.    !
!                                                                              !
!   Nested blocks not allowed.                                                 !
!==============================================================================!
subroutine io_read_input_file(filename, table)
  use hash_tables,  only: hash_table_init
  implicit none

  character(len=*), intent(in)    ::  filename
  type(hash_table), intent(inout) ::  table


  character(len=max_line_len)                               ::  line_str
  character(len=max_line_len),  allocatable,  dimension(:)  ::  blockdata
  character(len=6)  ::  buffer_format               ! only for max_line_len up to 3 digits (999)
  character(len=1)  ::  str_type, block_line_type
  integer           ::  nlines, nblanklines
  integer           ::  nblocklines, nblockblanks
  integer           ::  iblock, iblockdata
  integer           ::  file_unit
  integer           ::  istat

  ! no need to check length of filename - this is done in io_open_file
  call io_open_file(trim(filename), file_unit, 'read')

  ! allow table to expand if we have many collisions but don't shrink because we will be
  ! removing (hopefully) all elements anyway (= lots of shrinking = expensive) 
  ! and then destroying entire table when we're done
  call hash_table_init(table, nbuckets=16, can_expand=.true., can_shrink=.false.)
  
  nlines = 0
  nblanklines = 0
  block_line_type = 'e' ! assume end of block unless we encounter one in the file

  ! make sure we read into the whole of the string length..
  write(buffer_format, '(I3)') max_line_len
  buffer_format = '('//trim(adjustl(buffer_format))//'A)'

  do
    ! read all lines..
    read(file_unit, fmt=buffer_format, iostat=istat, end=100) line_str
    if (istat .ne. 0) call io_err("io_read_input_file: Error reading from "//trim(filename))

    nlines = nlines + 1

    if (len_trim(line_str) .le. 0) then
      ! blank line
      nblanklines = nblanklines + 1
      cycle ! rather than having lots of branches
    end if

    ! not a blank line (len_trim(line_str) > 0)
    str_type = io_str_get_type(line_str)

    ! decide what to do with each type:
    select case (str_type)
      case ('n')  ! normal line
        call io_parse_store_str(line_str, table)

      case ('c')  ! comment
        nblanklines = nblanklines + 1

      case ('b')  ! begin block
        ! store the entire block in one go

        ! first loop over block to count number of lines (so we can allocate array)
        nblocklines = 1         ! start of block line
        nblockblanks = 0
        block_line_type = 'n'   ! we hope they're all "normal" lines... (don't allow for nested blocks)

        ! count various numbers of lines (and check for nested block)
        do while (block_line_type .ne. 'e') ! line isn't 'end block type'
          read(file_unit, fmt=buffer_format, iostat=istat, end=100) line_str
          if (istat .ne. 0) call io_err("io_read_input_file: Error reading from "//trim(filename))

          nblocklines = nblocklines + 1

          if (len_trim(line_str) .le. 0) then
            ! blank line inside block
            nblockblanks = nblockblanks + 1
            cycle
          end if

          ! not a blank line (len_trim(line_str) > 0)
          block_line_type = io_str_get_type(line_str)

          ! treat comments as blanks and give error if we have nested block
          if (block_line_type .eq. 'c') then
            nblockblanks = nblockblanks + 1
          else if (block_line_type .eq. 'b') then
            call io_err("io_read_input_file: Nested block in "//trim(filename))
          end if
        end do

        ! make sure block not empty (subtract two to take account for beginning and end block lines)
        if (nblocklines-2 .eq. nblockblanks) call io_err("io_read_input_file: Block is empty")

        ! block_line_type should be 'e' on the last line of a block (or if we haven't passed through one)
        ! now we backspace through the file to beginning of block
        do iblock = 1, nblocklines
          backspace(file_unit, iostat=istat)
          if (istat .ne. 0) call io_err("io_read_input file: Error rewinding through "//trim(filename))
        end do

        ! allocate array to store block (subtract one element for the end block line):
        allocate(blockdata(nblocklines-nblockblanks-1), stat=istat)
        if (istat .ne. 0) call io_err("io_read_input_file: Error allocating array to store block data")

        ! store the block data:
        iblockdata = 1
        do iblock = 1, nblocklines-1   ! we don't care about the last line (end block line) 
          read(file_unit, fmt=buffer_format, iostat=istat) line_str
          if (istat .ne. 0) call io_err("io_read_input_file: Error reading block data in "//trim(filename))

          if (len_trim(line_str) .gt. 0) then
            ! not blank - check it's not a comment - then store
            ! no need to check for nested blocks anymore
            block_line_type = io_str_get_type(line_str)
            if (block_line_type .ne. 'c') then
              blockdata(iblockdata) = line_str
              iblockdata = iblockdata + 1
            end if
          end if
        end do

        ! read the last line so we don't get an unmatched end block error later
        read(file_unit, fmt=buffer_format, iostat=istat) line_str
        if (istat .ne. 0) call io_err("io_read_input_file: Error reading block data in "//trim(filename))

        ! could check for mismatched begin and end block names here.. for now, let's not bother

        ! reset block_line_type so we don't get error thinking we left mid-way through a block
        block_line_type = 'e'

        ! finally, we have all the block stored and are ready to parse the data and store in hash table
        call io_parse_store_block(blockdata, table)

        ! done with block data
        deallocate(blockdata, stat=istat)
        if (istat .ne. 0) call io_err("io_read_input_file: Error deallocating block data")

      case ('e') ! unpaired end block
        call io_err("io_read_input_file: Error in "//trim(filename)//" unpaired end block command")

      case default
      ! should not ever get here
      call io_err("io_read_input_file: Error getting correct line type in "//trim(filename))
    end select
  end do
  100 continue

  ! give error if we left the file in the middle of a block
  if (block_line_type .ne. 'e') call io_err("io_read_input_file: Exited "//trim(filename)//" in the middle of a block.")

  ! Give error if file contains nothing
  if (nlines .eq. nblanklines) call io_err("io_read_input_file: Nothing in "//trim(filename))

  call io_close_file(file_unit)
end subroutine io_read_input_file


!****s* io/io_input_finalize ==================================================!
! NAME                                                                         !
!   io_input_finalize (PUBLIC)                                                 !
!------------------------------------------------------------------------------!
! ARGUMENTS                                                                    !
!   None.                                                                      !
!------------------------------------------------------------------------------!
! DESCRIPTION                                                                  !
!   Writes any keywords that remain in the input hash table to the error file  !
!   and then stops. Hopefully all keywords are removed by calls to             !
!   io_input_get_data.                                                         !
!                                                                              !
!   If no keywords remain, then the hash table is destroyed and the program    !
!   can continue as normal.                                                    !
!==============================================================================!
subroutine io_input_finalize
  use iso_fortran_env,  only: output_unit 
  use hash_tables,      only: hash_table_size, hash_table_get_keywords, hash_table_destroy
  implicit none

  character(len=max_line_len),  allocatable,  dimension(:)  ::  keywords
  integer ::  ikeyword, nkeywords
  integer ::  istat

  nkeywords = hash_table_size(input_table)

  if (nkeywords .gt. 0) then
    allocate(keywords(nkeywords), stat=istat)
    if (istat .ne. 0) call io_err("io_input_finalize: Could not allocate keyword array")

    call hash_table_get_keywords(input_table, keywords)

    ! if io_err was called before here, we won't get this far.. so no need to worry about the rest:
    call io_err_initialize        ! open error file if we can

    write(err_unit_num, *) "The following keywords were found in the input file but not used:"
    write(err_unit_num, *)
    do ikeyword = 1, nkeywords
      write(err_unit_num, *) trim(keywords(ikeyword))
    end do

    write(err_unit_num, *)
    write(err_unit_num, *) "It may be possible that these keywords are valid, but not recognised for the calculation " &
    & //"you are trying to do."
    write(err_unit_num, *) "Please comment out/delete the relevant lines and try again to continue."

    deallocate(keywords, stat=istat)
    if (istat .ne. 0) call io_err("io_input_finalize: Could not deallocate keyword array")
    
    ! did io_err_initialize open a file?
    if (err_unit_num .ne. output_unit) call io_close_file(err_unit_num)
  end if

  ! might not be worth destroying if we're going to stop anyway, but let's try and deallocate/free up as much memory as possible
  call hash_table_destroy(input_table)

  if (nkeywords .gt. 0) stop
end subroutine io_input_finalize


!****f* io/io_query_keyword ===================================================!
! NAME                                                                         !
!   io_query_keyword (PUBLIC)                                                  !
!------------------------------------------------------------------------------!
! ARGUMENTS                                                                    !
!   character(len=*), intent(in)  ::  keyword                                  !
!------------------------------------------------------------------------------!
! DESCRIPTION                                                                  !
!   Queries the default input hash table for a keyword (see if it existed in   !
!   the input file).                                                           !
!------------------------------------------------------------------------------!
! RESULT                                                                       !
!   Same as hash_table_query:                                                  !
!   Returns number of matches for a given keyword.                             !
!     - 0 if keyword doesn't exist                                             !
!     - 1 if keyword matches a single value (scalar)                           !
!     - N>1 if the keyword matches a block (number of lines)                   !
!------------------------------------------------------------------------------!
! NOTES                                                                        !
!   Since this is just a wrapper for hash_table_query, we do not allow for     !
!   the desired table to be passed as an input. For more flexibility, use that !
!   routine instead.                                                           !
!==============================================================================!
function io_query_keyword(keyword)
  use hash_tables,  only: hash_table_query
  implicit none

  character(len=*), intent(in)  ::  keyword
  integer                       ::  io_query_keyword

  ! although hash_table_query checks length of input, check it here so we can make use of io_err
  if (len_trim(keyword) .gt. max_line_len) call io_err("io_query_keyword: keyword too long")
  
  io_query_keyword = hash_table_query(input_table, trim(keyword))
end function io_query_keyword


!****s* io/io_input_get_data ==================================================!
! NAME                                                                         !
!   io_input_get_data (PUBLIC)                                                 !
!------------------------------------------------------------------------------!
! ARGUMENTS                                                                    !
!   character(len=*),                           intent(in)  ::  keyword        !
!   character(len=max_line_len),  allocatable,  intent(out) ::  values(:)      !
!   logical,                                    intent(out) ::  found          !
!------------------------------------------------------------------------------!
! DESCRIPTION                                                                  !
!   Essentially just a wrapper for hash_table_get_remove, but which always     !
!   acts on input_table (the default IO hash table).                           !
!                                                                              !
!   Gets the values (either single line or block) associated with the keyword  !
!   from input_table, and then removes it.                                     !
!                                                                              !
!   found is .true. if the keyword existed in the input file/hash table.       !
!------------------------------------------------------------------------------!
! NOTES                                                                        !
!   Calls io_query_keyword to return the number of lines associated with the   !
!   keyword in the hash table. We don't give an error if there is no match.    !
!                                                                              !
!   In order to simplify use, we make use of an array for the values (even if  !
!   we only have a single line) - this is slower, but we don't have to treat   !
!   blocks or single lines differently.                                        !
!                                                                              !
!   Note: values array will need deallocating after last call to this routine. !
!   To be safe, will have to check if (allocated(values)).                     !
!==============================================================================!
subroutine io_input_get_data(keyword, values, found)
  use hash_tables,  only: hash_table_get_remove
  implicit none

  character(len=*),                                         intent(in)  ::  keyword
  character(len=max_line_len),  allocatable,  dimension(:), intent(out) ::  values
  logical,                                                  intent(out) ::  found

  integer ::  nlines, istat

  nlines = io_query_keyword(keyword)

  if (nlines .gt. 0) then
    found = .true.

    ! make sure we do not need to deallocate block before each call (only after the last one)
    if (allocated(values) .and. (size(values,1) .ne. nlines)) then
      deallocate(values, stat=istat)
      if (istat .ne. 0) call io_err("io_input_get_data: Could not deallocate values array")
    end if

    allocate(values(nlines), stat=istat)
    if (istat .ne. 0) call io_err("io_input_get_data: Could not allocate values array")

    call hash_table_get_remove(input_table, keyword, values)
  else
    ! don't deallocate unnecessarily
    ! might be able to use previously allocated array in next call
    found = .false.
  end if
end subroutine io_input_get_data


subroutine io_input_get_single_value_int(keyword, int_return, int_default, required)
  implicit none

  character(len=*),           intent(in)  ::  keyword
  integer,                    intent(out) ::  int_return
  integer,          optional, intent(in)  ::  int_default
  logical,          optional, intent(in)  ::  required

  character(len=max_line_len),  allocatable,  dimension(:)  ::  input_data
  logical ::  kw_found
  logical ::  die
  integer ::  istat

  ! do we require that the routine find something in the input table?
  if (present(required) .and. required) then
    die = .true.
  else
    die = .false.
  end if

  call io_input_get_data(trim(keyword), input_data, kw_found)
  if (kw_found) then
    ! if keyword found - array must be allocated
    ! check number of lines
    if (size(input_data, 1) .ne. 1) &
      & call io_err("io_input_get_single_value_int: Expected single line for "//trim(keyword))

    ! check number of tokens
    if (io_str_get_num_tokens(input_data(1)) .ne. 1) &
      & call io_err("io_input_get_single_value_int: Expected single value for "//trim(keyword))

    ! set:
    int_return = io_str_to_int(io_str_get_token(input_data(1), 1))
  else
    if (die) call io_err("io_input_get_single_value_int: "//trim(keyword)//" must exist in input")
    ! set default value if we can..
    if (present(int_default)) int_return = int_default
  end if

  if (allocated(input_data)) then
    deallocate(input_data, stat=istat)
    if (istat .ne. 0) call io_err("io_input_get_single_value_int: Could not deallocate input_data array")
  end if
end subroutine io_input_get_single_value_int


subroutine io_input_get_single_value_real(keyword, real_return, real_default, required)
  implicit none

  character(len=*),           intent(in)  ::  keyword
  real(kind=dp),              intent(out) ::  real_return
  real(kind=dp),    optional, intent(in)  ::  real_default
  logical,          optional, intent(in)  ::  required

  character(len=max_line_len),  allocatable,  dimension(:)  ::  input_data
  logical ::  kw_found
  logical ::  die
  integer ::  istat

  ! do we require that the routine find something in the input table?
  if (present(required) .and. required) then
    die = .true.
  else
    die = .false.
  end if

  call io_input_get_data(trim(keyword), input_data, kw_found)
  if (kw_found) then
    ! if keyword found - array must be allocated
    ! check number of lines
    if (size(input_data, 1) .ne. 1) &
      & call io_err("io_input_get_single_value_real: Expected single line for "//trim(keyword))

    ! check number of tokens
    if (io_str_get_num_tokens(input_data(1)) .ne. 1) &
      & call io_err("io_input_get_single_value_real: Expected single value for "//trim(keyword))

    ! set:
    real_return = io_str_to_real(io_str_get_token(input_data(1), 1))
  else
    if (die) call io_err("io_input_get_single_value_real: "//trim(keyword)//" must exist in input")
    ! set default value if we can..
    if (present(real_default)) real_return = real_default
  end if

  if (allocated(input_data)) then
    deallocate(input_data, stat=istat)
    if (istat .ne. 0) call io_err("io_input_get_single_value_real: Could not deallocate input_data array")
  end if
end subroutine io_input_get_single_value_real


subroutine io_input_get_single_value_logical(keyword, logical_return, logical_default, required)
  implicit none

  character(len=*),           intent(in)  ::  keyword
  logical,                    intent(out) ::  logical_return
  logical,          optional, intent(in)  ::  logical_default
  logical,          optional, intent(in)  ::  required

  character(len=max_line_len),  allocatable,  dimension(:)  ::  input_data
  logical ::  kw_found
  logical ::  die
  integer ::  istat

  ! do we require that the routine find something in the input table?
  if (present(required) .and. required) then
    die = .true.
  else
    die = .false.
  end if

  call io_input_get_data(trim(keyword), input_data, kw_found)
  if (kw_found) then
    ! if keyword found - array must be allocated
    ! check number of lines
    if (size(input_data, 1) .ne. 1) &
      & call io_err("io_input_get_single_value_logical: Expected single line for "//trim(keyword))

    ! check number of tokens
    if (io_str_get_num_tokens(input_data(1)) .ne. 1) &
      & call io_err("io_input_get_single_value_logical: Expected single value for "//trim(keyword))

    ! set:
    logical_return = io_str_to_logical(io_str_get_token(input_data(1), 1))
  else
    if (die) call io_err("io_input_get_single_value_logical: "//trim(keyword)//" must exist in input")
    ! set default value if we can..
    if (present(logical_default)) logical_return = logical_default
  end if

  if (allocated(input_data)) then
    deallocate(input_data, stat=istat)
    if (istat .ne. 0) call io_err("io_input_get_single_value_logical: Could not deallocate input_data array")
  end if
end subroutine io_input_get_single_value_logical


subroutine io_input_get_single_value_str(keyword, str_return, str_default, required)
  implicit none

  character(len=*),           intent(in)  ::  keyword
  character(len=*),           intent(out) ::  str_return
  character(len=*), optional, intent(in)  ::  str_default
  logical,          optional, intent(in)  ::  required

  character(len=max_line_len),  allocatable,  dimension(:)  ::  input_data
  logical ::  kw_found
  logical ::  die
  integer ::  istat

  if (len_trim(str_default) .gt. max_line_len) &
    & call io_err("io_input_get_single_value_str: str_default length > max_line_len")

  if (len(str_return) .gt. max_line_len) &
    & call io_err("io_input_get_single_value_str: str_return length > max_line_len")

  ! do we require that the routine find something in the input table?
  if (present(required) .and. required) then
    die = .true.
  else
    die = .false.
  end if

  call io_input_get_data(trim(keyword), input_data, kw_found)
  if (kw_found) then
    ! if keyword found - array must be allocated
    ! check number of lines
    if (size(input_data, 1) .ne. 1) &
      & call io_err("io_input_get_single_value_str: Expected single line for "//trim(keyword))

    ! check number of tokens
    if (io_str_get_num_tokens(input_data(1)) .ne. 1) &
      & call io_err("io_input_get_single_value_str: Expected single value for "//trim(keyword))

    ! set:
    str_return = io_str_get_token(input_data(1), 1)
  else
    if (die) call io_err("io_input_get_single_value_str: "//trim(keyword)//" must exist in input")
    ! set default value if we can..
    if (present(str_default)) str_return = str_default
  end if

  if (allocated(input_data)) then
    deallocate(input_data, stat=istat)
    if (istat .ne. 0) call io_err("io_input_get_single_value_str: Could not deallocate input_data array")
  end if
end subroutine io_input_get_single_value_str


!****s* io/io_input_hash_echo =================================================!
! NAME                                                                         !
!   io_input_hash_echo (PUBLIC)                                                !
!------------------------------------------------------------------------------!
! ARGUMENTS                                                                    !
!   logical,  optional, intent(in)  ::  suppress_output_file                   !
!------------------------------------------------------------------------------!
! DESCRIPTION                                                                  !
!   Outputs the contents of the hash table. By default we write to a file:     !
!   seedname.hash - but we have the option of suppressing extra files and      !
!   instead write to whatever 'stdout' is connected to.                        !
!                                                                              !
!   If we have no seedname, we write to whatever stdout is connected to.       !
!                                                                              !
!   Probably just used for debugging, but public anyway.                       !
!==============================================================================!
subroutine io_input_hash_echo(suppress_output_file)
  use hash_tables,  only: hash_table_list
  implicit none

  logical,  optional, intent(in)  ::  suppress_output_file
  integer ::  out_unit

  out_unit = stdout

  if ((.not. present(suppress_output_file)) .or. (.not. suppress_output_file)) then
    ! only create output file if we have a seedname - regardless of suppress_output_file
    if (len_trim(seedname) .ge. 1) then
      call io_open_file(trim(seedname)//'.hash', out_unit, 'replace')
    end if
  end if

  call hash_table_list(input_table, out_unit)

  ! only close if we opened a new unit:
  if (out_unit .ne. stdout) call io_close_file(out_unit)
end subroutine io_input_hash_echo


!****f* io/io_str_get_num_tokens ==============================================!
! NAME                                                                         !
!   io_str_get_num_tokens (PUBLIC)                                             !
!------------------------------------------------------------------------------!
! ARGUMENTS                                                                    !
!   character(len=*), intent(in)  ::  str                                      !
!------------------------------------------------------------------------------!
! DESCRIPTION                                                                  !
!   Counts the number of 'tokens' (or 'words') in a string.                    !
!------------------------------------------------------------------------------!
! NOTES                                                                        !
!   Only works if the tokens are separated by single spaces. ie: We assume     !
!   that io_str_shrink has been called beforehand.                             !
!==============================================================================!
function io_str_get_num_tokens(str)
  implicit none

  integer                       ::  io_str_get_num_tokens
  character(len=*), intent(in)  ::  str

  character(len=1), parameter ::  space = ' '
  integer ::  i, j, strlen

  ! if we assume that we have used str_shrink before then we can just count the spaces..
  strlen = len_trim(str)
  io_str_get_num_tokens = 1 ! always N-1 spaces to N tokens
  i = 1 ! get into loop
  j = 1
  do while (i .gt. 0)
    i = index(str(j:strlen), space)
    j = j+i
    if (i .gt. 0) io_str_get_num_tokens = io_str_get_num_tokens + 1
  end do
end function io_str_get_num_tokens


!****f* io/io_str_get_token ===================================================!
! NAME                                                                         !
!   io_str_get_token (PUBLIC)                                                  !
!------------------------------------------------------------------------------!
! ARGUMENTS                                                                    !
!   character(len=*), intent(in)  ::  str                                      !
!   integer,          intent(in)  ::  ntoken                                   !
!------------------------------------------------------------------------------!
! DESCRIPTION                                                                  !
!   Returns the ntoken-th token/'word' of a string.                            !
!------------------------------------------------------------------------------!
! NOTES                                                                        !
!   Only works if the tokens are separated by single spaces. ie: We assume     !
!   that io_str_shrink has been called beforehand.                             !
!                                                                              !
!   If token does not exist, we give an error. Therefore the user should have  !
!   called io_str_get_num_tokens beforehand to determine how many tokens exist !
!   within a string.                                                           !
!==============================================================================!
function io_str_get_token(str, ntoken)
  implicit none

  character(len=*), intent(in)  ::  str
  integer,          intent(in)  ::  ntoken
  character(len=len(str))       ::  io_str_get_token

  character(len=1), parameter   ::  space = ' '
  integer ::  i, j, itoken
  integer ::  strlen

  strlen = len_trim(str)
  
  ! assume we have shrunk the string..
  j = 0

  do itoken = 1, ntoken
    if (j .ge. strlen) call io_err("io_str_get_token: could not find token")
    i = index(str(j+1:strlen)//space, space, back=.false.)
    i = i+j
    io_str_get_token = str(j+1:i-1)
    j = i
  end do
end function io_str_get_token


!****s* io/io_str_to_lcase ====================================================!
! NAME                                                                         !
!   io_str_to_lcase (PUBLIC)                                                   !
!------------------------------------------------------------------------------!
! ARGUMENTS                                                                    !
!   character(len=*), intent(inout) ::  str                                    !
!------------------------------------------------------------------------------!
! DESCRIPTION                                                                  !
!   Converts a string to lowercase.                                            !
!==============================================================================!
subroutine io_str_to_lcase(str)
  implicit none

  character(len=*), intent(inout) ::  str

  character(len=len(str))       ::  strcopy ! for manipulations
  character(len=26), parameter  ::  lcase='abcdefghijklmnopqrstuvwxyz'
  character(len=26), parameter  ::  ucase='ABCDEFGHIJKLMNOPQRSTUVWXYZ'
  integer :: i, imin, imax, ipos
  character(len=1) :: firstchar
  
  strcopy = adjustl(str)
  firstchar = strcopy(1:1)
  imin = index(str, firstchar, back=.false.)
  imax = len_trim(str)
  do i = imin, imax
    ipos = index(ucase, str(i:i), back=.false.)
    if (ipos .gt. 0) str(i:i) = lcase(ipos:ipos)
  end do
end subroutine io_str_to_lcase


!****s* io/io_str_to_ucase ====================================================!
! NAME                                                                         !
!   io_str_to_ucase (PUBLIC)                                                   !
!------------------------------------------------------------------------------!
! ARGUMENTS                                                                    !
!   character(len=*), intent(inout) ::  str                                    !
!------------------------------------------------------------------------------!
! DESCRIPTION                                                                  !
!   Converts a string to uppercase.                                            !
!==============================================================================!
subroutine io_str_to_ucase(str)
  implicit none

  character(len=*), intent(inout) ::  str

  character(len=len(str))       ::  strcopy ! for manipulations
  character(len=26), parameter  ::  lcase='abcdefghijklmnopqrstuvwxyz'
  character(len=26), parameter  ::  ucase='ABCDEFGHIJKLMNOPQRSTUVWXYZ'
  integer :: i, imin, imax, ipos
  character(len=1) :: firstchar
  
  strcopy = adjustl(str)
  firstchar = strcopy(1:1)
  imin = index(str, firstchar, back=.false.)
  imax = len_trim(str)
  do i = imin, imax
    ipos = index(lcase, str(i:i), back=.false.)
    if (ipos .gt. 0) str(i:i) = ucase(ipos:ipos)
  end do
end subroutine io_str_to_ucase


!****f* io/io_str_to_int ======================================================!
! NAME                                                                         !
!   io_str_to_int (PUBLIC)                                                     !
!------------------------------------------------------------------------------!
! ARGUMENTS                                                                    !
!   character(len=*), intent(in)  ::  str                                      !
!------------------------------------------------------------------------------!
! DESCRIPTION                                                                  !
!   Converts a string to an integer using internal read.                       !
!------------------------------------------------------------------------------!
! NOTES                                                                        !
!   Some extra (and probably unnecessary) checks done to confirm that the      !
!   string is a valid number before we attempt to read.                        !
!   iostat .ne. 0 would probably pick this up though..                         !
!==============================================================================!
function io_str_to_int(str)
  implicit none

  character(len=*), intent(in)  ::  str
  integer                       ::  io_str_to_int
  character(len=12),  parameter ::  validchars = '-+0123456789'
  integer ::  strlen
  integer ::  istat
  integer ::  i

  strlen = len_trim(str)
  do i = 1, strlen
    if (scan(str(i:i), validchars, back=.false.) .ne. 1) then
      call io_err("io_str_to_int: Not a valid integer")
    end if
  end do
  ! might want to use a format statement here to be sure:
  read(str(1:strlen), *, iostat=istat) io_str_to_int
  if (istat .ne. 0) call io_err("io_str_to_int: Error reading str")
end function io_str_to_int


!****f* io/io_str_to_real =====================================================!
! NAME                                                                         !
!   io_str_to_real (PUBLIC)                                                    !
!------------------------------------------------------------------------------!
! ARGUMENTS                                                                    !
!   character(len=*), intent(in)  ::  str                                      !
!------------------------------------------------------------------------------!
! DESCRIPTION                                                                  !
!   Converts a string to a (double precision) real using internal read.        !
!------------------------------------------------------------------------------!
! NOTES                                                                        !
!   Some extra (and probably unnecessary) checks done to confirm that the      !
!   string is a valid number before we attempt to read.                        !
!   iostat .ne. 0 would probably pick this up though..                         !
!==============================================================================!
function io_str_to_real(str)
  implicit none

  character(len=*), intent(in)  ::  str
  real(kind=dp)                 ::  io_str_to_real
  character(len=15),  parameter ::  validchars = '-+.Ee0123456789'
  integer ::  strlen
  integer ::  istat
  integer ::  i

  strlen = len_trim(str)
  do i = 1,strlen
    if(scan(str(i:i), validchars, back=.false.) .ne. 1) then
      call io_err("io_str_to_real: Not a valid real")
    end if
  end do
  ! might want to use a format statement here to be sure:
  read(str(1:strlen), *, iostat=istat) io_str_to_real
  if (istat .ne. 0) call io_err("io_str_to_int: Error reading str")
end function io_str_to_real


function io_str_to_logical(str)
  implicit none

  character(len=*), intent(in)  ::  str
  logical                       ::  io_str_to_logical

  ! assume everything lowercase... (io_str_to_lcase called first)
  select case (str)
    case ('1', 'true', 't', '.true.', '.t.', 'yes', 'y', 'on')
      io_str_to_logical = .true.
    case ('0', 'false', 'f', '.false.', '.f.', 'no', 'n', 'off')
      io_str_to_logical = .false.
    case default
      call io_err("io_str_to_logical: Invalid logical string")
  end select
end function io_str_to_logical


!==============================================================================!
! PRIVATE ROUTINES BELOW HERE:                                                 !
!==============================================================================!


!****s* io/io_err_initialize ==================================================!
! NAME                                                                         !
!   io_err_initialize (PRIVATE)                                                !
!------------------------------------------------------------------------------!
! ARGUMENTS                                                                    !
!   None.                                                                      !
!------------------------------------------------------------------------------!
! DESCRIPTION                                                                  !
!   Sets the err_in_use flag to .true. so that further errors do not cause a   !
!   stop (attempting to preserve errors). If we have a seedname, then we       !
!   attempt to open the error file and change err_unit_num.                    !
!------------------------------------------------------------------------------!
! NOTES                                                                        !
!   Any errors arising from the call to io_open_file will be written to stdout !
!   (which is the default value of err_unit_num).                              !
!==============================================================================!
subroutine io_err_initialize
  implicit none
  
  err_in_use = .true. ! any errors from now will be treated differently (no stop)
  
  ! if we have a seedname, open file (and change err_unit_num)
  if (len_trim(seedname) .gt. 0) then
    call io_open_file(trim(seedname)//error_file_extension, err_unit_num, 'append')
  end if
end subroutine io_err_initialize


!****s* io/io_parse_store_str =================================================!
! NAME                                                                         !
!   io_parse_store_str (PRIVATE)                                               !
!------------------------------------------------------------------------------!
! ARGUMENTS                                                                    !
!   character(len=*), intent(inout) ::  str                                    ! 
!   type(hash_table), intent(inout) ::  table                                  !
!------------------------------------------------------------------------------!
! DESCRIPTION                                                                  !
!   Parses a string and stores it in the hash table.                           !
!                                                                              !
!   String is converted to lowercase, endline comments removed, assignment     !
!   operators are converted to spaces, then all redundant spaces are removed,  !
!   such that a single space separates the keyword and each value in a line.   !
!==============================================================================!
subroutine io_parse_store_str(str, table)
  use hash_tables,  only: hash_table_add
  implicit none

  character(len=*), intent(inout) ::  str
  type(hash_table), intent(inout) ::  table

  character(len=len(str)) ::  keyword, values

  call io_str_to_lcase(str)                             ! make lowercase
  call io_str_remove_comments(str)                      ! remove comments
  call io_str_assignment_to_space(str)                  ! make assignment operators spaces
  call io_str_shrink(str)                               ! remove unnecessary spaces and shift to left
  keyword = io_str_get_token(str, 1)                    ! keyword is first token in string
  values = str(len_trim(keyword)+2:len(str))            ! values is the rest of the string (single spaces)

  ! add to hash table
  call hash_table_add(table, keyword, values)
end subroutine io_parse_store_str


!****s* io/io_parse_store_block ===============================================!
! NAME                                                                         !
!   io_parse_store_block (PRIVATE)                                             !
!------------------------------------------------------------------------------!
! ARGUMENTS                                                                    !
!   character(len=*), dimension(:), intent(inout) ::  blockdata                !
!   type(hash_table),               intent(inout) ::  table                    !
!------------------------------------------------------------------------------!
! DESCRIPTION                                                                  !
!   Parses a block and stores it in the hash table.                            !
!                                                                              !
!   Each line is converted to lowercase, endline comments removed, then all    !
!   redundant spaces are removed, such that a single space separates each      !
!   value in a line. The block keyword comes from the second token of the      !
!   first line, after these operations have been carried out.                  !
!==============================================================================!
subroutine io_parse_store_block(blockdata, table)
  use hash_tables,  only: hash_table_add
  implicit none

  character(len=*), dimension(:), intent(inout) ::  blockdata
  type(hash_table),               intent(inout) ::  table

  integer ::  iline, nlines

  nlines = size(blockdata,1)

  ! similar to treating a single line, except we don't remove assignment operators
  do iline = 1, nlines
    call io_str_to_lcase(blockdata(iline))                        ! make lowercase
    call io_str_remove_comments(blockdata(iline))                 ! remove endline comments
    call io_str_shrink(blockdata(iline))                          ! remove unnecessary spaces and shift to left
  end do

  ! get keyword:
  ! first let's check to see if we actually have a keyword (could do this sooner, but doing it here keeps io_read_input_file neat)
  if (io_str_get_num_tokens(blockdata(1)) .ne. 2) then
    call io_err("io_parse_store_block: Block does not have exactly one keyword associated with it")
  end if

  blockdata(1) = io_str_get_token(blockdata(1), 2)                ! keyword is second token in first line

  ! can now store in hash table
  call hash_table_add(table, blockdata(1), blockdata(2:nlines))
end subroutine io_parse_store_block


!****f* io/io_str_get_type ====================================================!
! NAME                                                                         !
!   io_str_get_type (PRIVATE)                                                  !
!------------------------------------------------------------------------------!
! ARGUMENTS                                                                    !
!   character(len=*), intent(in)  ::  str                                      !
!------------------------------------------------------------------------------!
! DESCRIPTION                                                                  !
!   Determines what type of line a string is in the input file.                !
!------------------------------------------------------------------------------!
! RESULT                                                                       !
!   Returns a character representing the type of string (in an input file).    !
!     - 'n' for a normal line                                                  !
!     - 'b' for a begin block line                                             !
!     - 'e' for an end block line                                              !
!     - 'c' for a comment line                                                 !
!------------------------------------------------------------------------------!
! NOTES                                                                        !
!   Input string must not be blank (check length before calling).              !
!==============================================================================!
function io_str_get_type(str)
  implicit none

  character(len=1)              ::  io_str_get_type
  character(len=*), intent(in)  ::  str

  character(len=len(str)) ::  strcopy

  ! marker lengths
  integer,  parameter ::  len_comment_marker = len(input_comment_marker)
  integer,  parameter ::  len_begin_block_marker = len(input_begin_block_marker)+1  ! +1 to make sure not part of keyword
  integer,  parameter ::  len_end_block_marker = len(input_end_block_marker)+1      ! ditto
  integer,  parameter ::  max_marker_len = max(len_comment_marker, len_begin_block_marker, len_end_block_marker)

  if(len_trim(str) .lt. 1) call io_err("io_str_get_type: Error: blank string")

  io_str_get_type = 'n'                           ! normal

  strcopy = adjustl(str)
  call io_str_to_lcase(strcopy(1:max_marker_len)) ! no need to convert entire string just yet...

  if (strcopy(1:len_comment_marker) .eq. input_comment_marker) then
    io_str_get_type = 'c'                         ! comment
  ! add space to end of markers to make sure they're not part of a keyword
  else if (strcopy(1:len_begin_block_marker) .eq. input_begin_block_marker//' ') then
    io_str_get_type = 'b'                         ! block begin
  else if (strcopy(1:len_end_block_marker) .eq. input_end_block_marker//' ') then
    io_str_get_type = 'e'                         ! block end
  end if
end function io_str_get_type


!****s* io/io_str_remove_comments =============================================!
! NAME                                                                         !
!   io_str_remove_comments (PRIVATE)                                           !
!------------------------------------------------------------------------------!
! ARGUMENTS                                                                    !
!   character(len=*), intent(inout) ::  str                                    !
!------------------------------------------------------------------------------!
! DESCRIPTION                                                                  !
!   Blanks out all characters of a string from the first comment character     !
!   onwards.                                                                   !
!------------------------------------------------------------------------------!
! NOTES                                                                        !
!   Does not take the comment character being inside quotation marks into      !
!   account. This should not be a problem because it is unlikely that tokens   !
!   will contain these.                                                        !
!==============================================================================!
subroutine io_str_remove_comments(str)
  implicit none

  character(len=*), intent(inout) ::  str

  character(len=1), parameter     ::  space = ' '
  integer                         ::  commentpos, strlen

  strlen = len(str)
  commentpos = index(str, input_comment_marker, back=.false.)
  if (commentpos .gt. 0) then
    str(commentpos:strlen) = repeat(space, strlen-commentpos+1)
  end if
end subroutine io_str_remove_comments


!****s* io/io_str_assignment_to_space =========================================!
! NAME                                                                         !
!   io_str_assigment_to_space (PRIVATE)                                        !
!------------------------------------------------------------------------------!
! ARGUMENTS                                                                    !
!   character(len=*), intent(inout) ::  str                                    !
!------------------------------------------------------------------------------!
! DESCRIPTION                                                                  !
!   Converts all characters that are assignment operators to spaces in a       !
!   string.                                                                    !
!==============================================================================!
subroutine io_str_assignment_to_space(str)
  implicit none

  character(len=*), intent(inout) ::  str

  ! could make str_assign global, but these are likely all we'll ever want
  character(len=*), parameter ::  str_assign = ';:='    ! valid assignment operators
  character(len=1), parameter ::  space = ' '
  integer ::  i

  i = scan(str, str_assign)
  do while (i .ne. 0)
    str(i:i) = space
    i = scan(str, str_assign)
  end do
end subroutine io_str_assignment_to_space


!****s* io/io_str_shrink ======================================================!
! NAME                                                                         !
!   io_str_shrink (PRIVATE)                                                    !
!------------------------------------------------------------------------------!
! ARGUMENTS                                                                    !
!   character(len=*), intent(inout) ::  str                                    !
!------------------------------------------------------------------------------!
! DESCRIPTION                                                                  !
!   Shifts a string to the left and removes all redundant spaces (ie: every    !
!   'word' is separated by only a single space).                               !
!==============================================================================!
subroutine io_str_shrink(str)
  implicit none

  character(len=*), intent(inout) ::  str

  character(len=1), parameter ::  space = ' '
  integer ::  i, j, strlen

  ! i is abs position in string, j is relative in substring
  i=0
  str = adjustl(str)
  strlen = len_trim(str)
  j = index(str, space)
  do while (j .gt. 0)
    i=j+i
    str(i+1:strlen) = adjustl(str(i+1:strlen))
    j = index(str(i+1:strlen), space)
  end do
end subroutine io_str_shrink
end module io
