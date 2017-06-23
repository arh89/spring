!****m* hash_tables/hash_tables ===============================================!
! NAME                                                                         !
!   hash_tables                                                                !
!------------------------------------------------------------------------------!
! DESCRIPTION                                                                  !
!   An implementation of hash tables (associative arrays).                     !
!                                                                              !
!   Internally the hash table is stored as an allocatable array of doubly      !
!   linked lists of elements - linked lists are used for collision management  !
!   (where multiple keywords hash to the same bucket).                         !
!   Each element of a list is used to store either a keyword-value pair, or a  !
!   keyword-block pair (where the values associated with the keyword are       !
!   stored in a 1D array of strings).                                          !
!                                                                              !
!   Each table is able to shrink/expand depending on the current load factor   !
!   of the table, to balance perfomance and memory usage.                      !
!   This can be set on a table by table basis.                                 !
!                                                                              !
!   Caching is also supported in order to increase the performance of repeated !
!   access to the same element of a table.                                     !
!   This is enabled at the module level.                                       !
!------------------------------------------------------------------------------!
! AUTHOR                                                                       !
!   Aaron Hopkinson                                                            !
!****==========================================================================!
module hash_tables
  implicit none

  private

  integer,  parameter,  public  ::  str_len = 128                 ! max length of strings stored in table

  integer,  parameter           ::  i32 = selected_int_kind(9)
  integer,  parameter           ::  i64 = selected_int_kind(18)

  ! Performance related parameters:
  integer,  parameter ::  default_num_buckets = 16
  logical,  parameter ::  default_can_expand = .true.
  logical,  parameter ::  default_can_shrink = .false.

  ! These values are module level and unable to be changed on a per table basis:
  real,     parameter ::  max_load_factor = 0.75
  real,     parameter ::  min_load_factor = 0.25
  real,     parameter ::  expand_scale    = 2.00
  real,     parameter ::  shrink_scale    = 1.0/expand_scale
  logical,  parameter ::  enable_cache = .true.
  
  ! Will need lists of these for collision management
  type :: hash_element
    integer                                           ::  hash_id
    character(len=str_len)                            ::  keyword
    character(len=str_len)                            ::  value
    character(len=str_len), dimension(:), allocatable ::  block
    type(hash_element),     pointer                   ::  next_element
    type(hash_element),     pointer                   ::  prev_element
  end type hash_element

  type :: hash_pointer
    type(hash_element), pointer ::  head
    type(hash_element), pointer ::  tail
  end type hash_pointer

  type, public  ::  hash_table
    logical                                       ::  initialized = .false.
    real                                          ::  load_factor
    logical                                       ::  can_expand
    logical                                       ::  can_shrink
    integer                                       ::  nvalues
    integer                                       ::  nbuckets
    type(hash_pointer), dimension(:), allocatable ::  buckets
    type(hash_element), pointer                   ::  last_access
    character(len=str_len)                        ::  last_keyword
  end type hash_table


  public  ::  hash_table_init
  public  ::  hash_table_size
  public  ::  hash_table_query
  public  ::  hash_table_add
  public  ::  hash_table_get
  public  ::  hash_table_remove
  public  ::  hash_table_get_remove
  public  ::  hash_table_get_keywords
  public  ::  hash_table_resize
  public  ::  hash_table_list
  public  ::  hash_table_destroy
  public  ::  hash_function

  interface hash_table_add
    module procedure hash_table_add_value
    module procedure hash_table_add_block
  end interface hash_table_add

  interface hash_table_get
    module procedure hash_table_get_value
    module procedure hash_table_get_block
  end interface hash_table_get

  interface hash_table_get_remove
    module procedure hash_table_get_remove_value
    module procedure hash_table_get_remove_block
  end interface hash_table_get_remove

contains

!****s* hash_tables/hash_table_init ===========================================!
! NAME                                                                         !
!   hash_table_init (PUBLIC)                                                   !
!------------------------------------------------------------------------------!
! DESCRIPTION                                                                  !
!   Sets up the hash table ready for use.                                      !
!   Buckets are allocated and internal state variables are reset.              !
!------------------------------------------------------------------------------!
! ARGUMENTS                                                                    !
!   type(hash_table),   intent(inout) ::  table                                !
!                                                                              !
!   (to override module defaults):                                             !
!   integer,  optional, intent(in)    ::  nbuckets                             !
!   logical,  optional, intent(in)    ::  can_expand                           !
!   logical,  optional, intent(in)    ::  can_shrink                           !
!------------------------------------------------------------------------------!
! NOTES                                                                        !
!   Table should not already be initialized.                                   ! 
!                                                                              !
!   If nbuckets is present and <1, then we do not give an error but instead    !
!   allocate the default number.                                               !
!                                                                              !
!   Optional arguments are necessary only to optimize memory usage/speed.      !
!==============================================================================!
subroutine hash_table_init(table, nbuckets, can_expand, can_shrink)
  implicit none
  type(hash_table),   intent(inout) ::  table
  integer,  optional, intent(in)    ::  nbuckets
  logical,  optional, intent(in)    ::  can_expand
  logical,  optional, intent(in)    ::  can_shrink
  integer ::  ibucket, istat

  if (table%initialized) stop 'Error in hash_table_init: table already allocated/initialized'

  table%nvalues = 0
  table%load_factor = 0.0
  table%can_shrink = default_can_shrink
  table%can_expand = default_can_expand

  if (present(can_shrink)) table%can_shrink = can_shrink
  if (present(can_expand)) table%can_expand = can_expand

  ! don't give error if nbuckets < 1, just assign default..
  table%nbuckets = default_num_buckets
  if (present(nbuckets) .and. (nbuckets .ge. 1)) table%nbuckets = nbuckets

  allocate(table%buckets(table%nbuckets), stat=istat)
  if (istat .ne. 0) stop 'Error in hash_table_init: could not allocate hash table'

  do ibucket = 1, table%nbuckets
    nullify(table%buckets(ibucket)%head)
    nullify(table%buckets(ibucket)%tail)
  end do

  if (enable_cache) call hash_table_reset_cache(table)

  table%initialized = .true.
end subroutine hash_table_init


!****f* hash_tables/hash_table_size ===========================================!
! NAME                                                                         !
!   hash_table_size (PUBLIC)                                                   !
!------------------------------------------------------------------------------!
! DESCRIPTION                                                                  !
!   Returns number of keyword-value pairs stored in the hash table (a block    !
!   is treated as a single pair). Note: This is not the number of buckets.     !
!------------------------------------------------------------------------------!
! ARGUMENTS                                                                    !
!   type(hash_table), intent(in)  ::  table                                    !
!------------------------------------------------------------------------------!
! NOTES                                                                        !
!   Table must be initialized.                                                 !
!==============================================================================!
function hash_table_size(table)
  type(hash_table), intent(in)  ::  table
  integer                       ::  hash_table_size

  if (.not. table%initialized) stop 'Error in hash_table_size: table not initialized'

  ! externally, we only care about the number of values - not the number of buckets
  hash_table_size = table%nvalues
end function hash_table_size


!****f* hash_tables/hash_table_query ==========================================!
! NAME                                                                         !
!   hash_table_query (PUBLIC)                                                  !
!------------------------------------------------------------------------------!
! DESCRIPTION                                                                  !
!   Queries a hash table for a given keyword - returns number of lines that    !
!   match.                                                                     !
!------------------------------------------------------------------------------!
! ARGUMENTS                                                                    !
!   type(hash_table), intent(inout) ::  table                                  !
!   character(len=*), intent(in)    ::  keyword                                !
!------------------------------------------------------------------------------!
! RESULT                                                                       !
!   Returns number of matches for a given keyword.                             !
!     - 0 if keyword doesn't exist                                             !
!     - 1 if keyword matches a single value (scalar)                           !
!     - N>1 if the keyword matches a block (number of lines)                   !
!------------------------------------------------------------------------------!
! NOTES                                                                        !
!   Table must be initialized.                                                 !
!==============================================================================!
function hash_table_query(table, keyword)
  implicit none
  type(hash_table), intent(inout) ::  table
  character(len=*), intent(in)    ::  keyword
  integer                         ::  hash_table_query
  ! local vars:
  type(hash_element), pointer ::  ptr
  integer ::  hash_id, bucket_id

  if (.not. table%initialized) stop 'Error in hash_table_query: table not initialized'

  call hash_table_find_element(table, keyword, hash_id, bucket_id, ptr)
  if (.not. associated(ptr)) then
    ! keyword does not exist
    hash_table_query = 0
  else
    if (allocated(ptr%block)) then
      hash_table_query = size(ptr%block,1)
    else
      hash_table_query = 1 ! scalar
    end if
  end if
end function hash_table_query


!****s* hash_tables/hash_table_add_value ======================================!
! NAME                                                                         !
!   hash_table_add_value (PUBLIC, as hash_table_add)                           !
!------------------------------------------------------------------------------!
! DESCRIPTION                                                                  !
!   Adds a keyword value pair to the table. In this case the value is a single !
!   line, rather than a block.                                                 !
!                                                                              !
!   If the keyword is already in the table, then its value will be overwritten !
!   by the new value.                                                          !
!------------------------------------------------------------------------------!
! ARGUMENTS                                                                    !
!   type(hash_table), intent(inout) ::  table                                  !
!   character(len=*), intent(in)    ::  keyword                                !
!   character(len=*), intent(in)    ::  value                                  !
!                                                                              !
!   Note that keyword and value must len_trim to no longer than the module     !
!   level str_len (as the hash table cannot store this).                       !
!------------------------------------------------------------------------------!
! NOTES                                                                        !
!   Table must be initialized.                                                 !
!                                                                              !
!   It is encouraged that the user calls the generic hash_table_add routine    !
!   instead. This eliminates having to worry about whether the 'value' is a    !
!   single line or a block.                                                    !
!==============================================================================!
recursive subroutine hash_table_add_value(table, keyword, value)
  ! must be recursive so that hash table can grow if necessary
  implicit none
  type(hash_table), intent(inout) ::  table
  character(len=*), intent(in)    ::  keyword
  character(len=*), intent(in)    ::  value
  ! local vars:
  type(hash_element), pointer ::  element
  type(hash_element), pointer ::  ptr
  integer ::  hash_id, bucket_id, istat

  if (.not. table%initialized) stop 'Error in hash_table_add_value: table not initialized'

  ! check we won't lose any information - keyword checked in hash_table_find_element
  if (len_trim(value) .gt. str_len) stop 'Error in hash_table_add_value: value too long'

  call hash_table_find_element(table, keyword, hash_id, bucket_id, ptr)

  ! check to see if keyword already exists (already have element with this)
  if (associated(ptr)) then

    ! deallocate if we have a block
    if (allocated(ptr%block)) then
      deallocate(ptr%block, stat=istat)
      if (istat .ne. 0) stop 'Error in hash_table_add_value: could not deallocate block'
    end if

    ! set value
    ptr%value = value
    nullify(ptr)                      ! no longer necessary
  else
    ! no, keyword doesn't exist.. create new element and append to list (could be head)
    ! we know hash_id from output of hash_table_find_element
    call hash_table_create_element_value(keyword, value, hash_id, element)
    if (enable_cache) call hash_table_update_cache(table, element)

    if (associated(table%buckets(bucket_id)%head)) then
      ! append to list
      table%buckets(bucket_id)%tail%next_element => element
      element%prev_element => table%buckets(bucket_id)%tail
    else
      ! create new list
      table%buckets(bucket_id)%head => element
      nullify(element%prev_element)
    end if

    table%buckets(bucket_id)%tail => element
    nullify(element%next_element)
    table%nvalues = table%nvalues + 1
    table%load_factor = real(table%nvalues)/real(table%nbuckets)
    nullify(element)                  ! no harm in doing this..
    if (table%can_expand .and. (table%load_factor .gt. max_load_factor)) call hash_table_resize(table, expand_scale)
  end if
end subroutine hash_table_add_value


!****s* hash_tables/hash_table_add_block ======================================!
! NAME                                                                         !
!   hash_table_add_block (PUBLIC, as hash_table_add)                           !
!------------------------------------------------------------------------------!
! DESCRIPTION                                                                  !
!   Adds a keyword-block pair to the table. In this case the value is a block  !
!   rather than a single line.                                                 !
!                                                                              !
!   If the keyword is already in the table, then its value will be overwritten !
!   by the new block.                                                          !
!------------------------------------------------------------------------------!
! ARGUMENTS                                                                    !
!   type(hash_table),               intent(inout) ::  table                    !
!   character(len=*),               intent(in)    ::  keyword                  !
!   character(len=*), dimension(:), intent(in)    ::  block                    !
!                                                                              !
!   Note that keyword and each line of the block must len_trim to no longer    !
!   than the module level str_len (as the hash table cannot store this).       !
!------------------------------------------------------------------------------!
! NOTES                                                                        !
!   Table must be initialized.                                                 !
!                                                                              !
!   It is encouraged that the user calls the generic hash_table_add routine    !
!   instead. This eliminates having to worry about whether the 'value' is a    !
!   single line or a block.                                                    !
!                                                                              !
!   Note that a block of length 1 (a single line) will be added to the table   !
!   in the same manner as a single line. (As if done by hash_table_add_value.) !
!==============================================================================!
recursive subroutine hash_table_add_block(table, keyword, block)
  ! must be recursive so that hash table can grow if necessary
  implicit none
  type(hash_table),               intent(inout) ::  table
  character(len=*),               intent(in)    ::  keyword
  character(len=*), dimension(:), intent(in)    ::  block
  ! local vars:
  type(hash_element), pointer ::  element
  type(hash_element), pointer ::  ptr
  integer ::  hash_id, bucket_id, iblock, blocklines, istat

  if (.not. table%initialized) stop 'Error in hash_table_add_block: table not initialized'

  blocklines = size(block,1)

  if (blocklines .lt. 1) stop 'Error in hash_table_add_block: block contains no values'

  ! check we won't lose any information - keyword checked in hash_table_find_element
  do iblock = 1, blocklines
    if (len_trim(block(iblock)) .gt. str_len) stop 'Error in hash_table_add_block: value too long in block'
  end do

  call hash_table_find_element(table, keyword, hash_id, bucket_id, ptr)

  ! check to see if keyword already exists (already have element with this)
  if (associated(ptr)) then

    ! if block is a single line long - just use the scalar 'value'
    if (blocklines .eq. 1)  then
      ptr%value = block(1)

      ! deallocate if necessary..
      if (allocated(ptr%block)) then
        deallocate(ptr%block, stat=istat)
        if (istat .ne. 0) stop 'Error in hash_table_add_block: cannot deallocate block'
      end if

    else
      ! block is muliple lines - reset 'value' and store in block

      ! can reuse block if it is exactly the right size
      if (allocated(ptr%block) .and. (size(ptr%block,1) .ne. blocklines)) then
        deallocate(ptr%block, stat=istat)
        if (istat .ne. 0) stop 'Error in hash_table_add_block: cannot deallocate block'
      end if

      if (.not. allocated(ptr%block)) then
        allocate(ptr%block(blocklines), stat=istat)
        if (istat .ne. 0) stop 'Error in hash_table_add_block: cannot allocate block'
      end if

      ptr%value = ' '
      ptr%block(:) = block(:)
    end if

    nullify(ptr)                      ! no longer necessary
  else
    ! no, keyword doesn't exist.. create new element and append to list (could be head)
    ! we know hash_id from output of hash_table_find_element
    if (blocklines .eq. 1) then
      ! again, use scalar 'value' if block is one line long
      call hash_table_create_element_value(keyword, block(1), hash_id, element)
    else
      call hash_table_create_element_block(keyword, block, hash_id, element)
    end if
    if (enable_cache) call hash_table_update_cache(table, element)

    if (associated(table%buckets(bucket_id)%head)) then
      ! append to list
      table%buckets(bucket_id)%tail%next_element => element
      element%prev_element => table%buckets(bucket_id)%tail
    else
      ! create new list
      table%buckets(bucket_id)%head => element
      nullify(element%prev_element)
    end if

    table%buckets(bucket_id)%tail => element
    nullify(element%next_element)
    table%nvalues = table%nvalues + 1
    table%load_factor = real(table%nvalues)/real(table%nbuckets)
    nullify(element)                  ! no harm in doing this..
    if (table%can_expand .and. (table%load_factor .gt. max_load_factor)) call hash_table_resize(table, expand_scale)
  end if
end subroutine hash_table_add_block


!****s* hash_tables/hash_table_get_value ======================================!
! NAME                                                                         !
!   hash_table_get_value (PUBLIC, as hash_table_get)                           !
!------------------------------------------------------------------------------!
! DESCRIPTION                                                                  !
!   Returns the value (single line) associated with the keyword in the hash    !
!   table.                                                                     !
!                                                                              !
!   If element_index is provided and is > 1, it is assumed that the keyword    !
!   is associated with a block and returns that line (line numbers starting    !
!   at 1).                                                                     !
!------------------------------------------------------------------------------!
! ARGUMENTS                                                                    !
!   type(hash_table),       intent(inout) ::  table                            !
!   character(len=*),       intent(in)    ::  keyword                          !
!   character(len=str_len), intent(out)   ::  value                            !
!                                                                              !
!   integer,  optional,     intent(in)    ::  element_index                    !
!------------------------------------------------------------------------------!
! NOTES                                                                        !
!   Table must be initialized.                                                 !
!                                                                              !
!   Keyword assumed to exist in table. Gives error if keyword not found.       !
!     ie: It is assumed that hash_table_query has already been called.         !
!==============================================================================!
subroutine hash_table_get_value(table, keyword, value, element_index)
  implicit none
  type(hash_table),       intent(inout) ::  table
  character(len=*),       intent(in)    ::  keyword
  character(len=str_len), intent(out)   ::  value
  integer,  optional,     intent(in)    ::  element_index
  ! local vars
  type(hash_element), pointer ::  ptr
  integer ::  hash_id, bucket_id

  if (.not. table%initialized) stop 'Error in hash_table_get_value: table not initialized'

  call hash_table_find_element(table, keyword, hash_id, bucket_id, ptr)

  ! hash_table_get_from_ptr_value will flag error if ptr not associated.. no need to deal with this here
  call hash_table_get_from_ptr_value(ptr, value, element_index)

  nullify(ptr) ! no harm in doing this..
end subroutine hash_table_get_value


!****s* hash_tables/hash_table_get_block ======================================!
! NAME                                                                         !
!   hash_table_get_block (PUBLIC, as hash_table_get)                           !
!------------------------------------------------------------------------------!
! DESCRIPTION                                                                  !
!   Returns the block associated with the keyword in the hash table.           !
!------------------------------------------------------------------------------!
! ARGUMENTS                                                                    !
!   type(hash_table),                     intent(inout) ::  table              !
!   character(len=*),                     intent(in)    ::  keyword            !
!   character(len=str_len), dimension(:), intent(out)   ::  block              !
!------------------------------------------------------------------------------!
! NOTES                                                                        !
!   Table must be initialized.                                                 !
!                                                                              !
!   Keyword assumed to exist in table. Gives error if keyword not found.       !
!     ie: It is assumed that hash_table_query has already been called.         !
!                                                                              !
!   The 'block' array *must* be allocated and large enough to store the entire !
!   block. (Size given by hash_table_query)                                    !
!==============================================================================!
subroutine hash_table_get_block(table, keyword, block)
  implicit none
  type(hash_table),                     intent(inout) ::  table
  character(len=*),                     intent(in)    ::  keyword
  character(len=str_len), dimension(:), intent(out)   ::  block
  ! local vars
  type(hash_element), pointer ::  ptr
  integer ::  hash_id, bucket_id
  
  if (.not. table%initialized) stop 'Error in hash_table_get_block: table not initialized'

  call hash_table_find_element(table, keyword, hash_id, bucket_id, ptr)

  ! hash_table_get_from_ptr_block will flag error if ptr not associated.. no need to deal with this here
  call hash_table_get_from_ptr_block(ptr, block)

  nullify(ptr) ! no harm in doing this..
end subroutine hash_table_get_block


!****s* hash_tables/hash_table_remove =========================================!
! NAME                                                                         !
!   hash_table_remove (PUBLIC)                                                 !
!------------------------------------------------------------------------------!
! DESCRIPTION                                                                  !
!   Removes a keyword-value pair from the hash table (the value associated     !
!   with a keyword can be a block).                                            !
!------------------------------------------------------------------------------!
! ARGUMENTS                                                                    !
!   type(hash_table), intent(inout) ::  table                                  !
!   character(len=*), intent(in)    ::  keyword                                !
!------------------------------------------------------------------------------!
! NOTES                                                                        !
!   Table must be initialized.                                                 !
!                                                                              !
!   Keyword must be in table (else we give an error) - ie: hash_table_query    !
!   assumed to be > 0.                                                         !
!==============================================================================!
subroutine hash_table_remove(table, keyword)
  implicit none
  type(hash_table), intent(inout) ::  table
  character(len=*), intent(in)    ::  keyword
  ! local vars:
  type(hash_element), pointer     ::  ptr
  integer :: hash_id, bucket_id

  if (.not. table%initialized) stop 'Error in hash_table_remove: table not initialized'

  call hash_table_find_element(table, keyword, hash_id, bucket_id, ptr)
  call hash_table_remove_from_ptr(table, ptr)
  if (enable_cache) call hash_table_reset_cache(table)
end subroutine hash_table_remove


!****s* hash_tables/hash_table_get_remove_value ===============================!
! NAME                                                                         !
!   hash_table_get_remove_value (PUBLIC, as hash_table_get_remove)             !
!------------------------------------------------------------------------------!
! DESCRIPTION                                                                  !
!   Gets a single line value associated with a keyword from a table and then   !
!   removes the keyword-value pair from the hash table.                        !
!------------------------------------------------------------------------------!
! ARGUMENTS                                                                    !
!   type(hash_table),       intent(inout) ::  table                            !
!   character(len=*),       intent(in)    ::  keyword                          !
!   character(len=str_len), intent(out)   ::  value                            !
!------------------------------------------------------------------------------!
! NOTES                                                                        !
!   Table must be initialized.                                                 !
!                                                                              !
!   The return/requested value must be stored as a single line in the hash     !
!   table.                                                                     !
!==============================================================================!
subroutine hash_table_get_remove_value(table, keyword, value)
  implicit none
  type(hash_table),       intent(inout) ::  table
  character(len=*),       intent(in)    ::  keyword
  character(len=str_len), intent(out)   ::  value
  ! local vars
  type(hash_element), pointer     ::  ptr
  integer ::  hash_id, bucket_id

  if (.not. table%initialized) stop 'Error in hash_table_get_remove_value: table not initialized'

  call hash_table_find_element(table, keyword, hash_id, bucket_id, ptr)
  call hash_table_get_from_ptr_value(ptr, value)
  call hash_table_remove_from_ptr(table, ptr)
  if (enable_cache) call hash_table_reset_cache(table)
end subroutine hash_table_get_remove_value


!****s* hash_tables/hash_table_get_remove_block ===============================!
! NAME                                                                         !
!   hash_table_get_remove_block (PUBLIC, as hash_table_get_remove)             !
!------------------------------------------------------------------------------!
! DESCRIPTION                                                                  !
!   Gets a block associated with a keyword from a table and then removes the   !
!   keyword-block pair from the hash table.                                    !
!------------------------------------------------------------------------------!
! ARGUMENTS                                                                    !
!   type(hash_table),                     intent(inout) ::  table              !
!   character(len=*),                     intent(in)    ::  keyword            !
!   character(len=str_len), dimension(:), intent(out)   ::  block              !
!------------------------------------------------------------------------------!
! NOTES                                                                        !
!   Table must be initialized.                                                 !
!                                                                              !
!   Output (block) array must be large enough to store the entire associated   !
!   block.                                                                     !
!==============================================================================!
subroutine hash_table_get_remove_block(table, keyword, block)
  implicit none
  type(hash_table),                     intent(inout) ::  table
  character(len=*),                     intent(in)    ::  keyword
  character(len=str_len), dimension(:), intent(out)   ::  block
  ! local vars
  type(hash_element), pointer     ::  ptr
  integer ::  hash_id, bucket_id

  if (.not. table%initialized) stop 'Error in hash_table_get_remove_block: table not initialized'

  call hash_table_find_element(table, keyword, hash_id, bucket_id, ptr)
  call hash_table_get_from_ptr_block(ptr, block)
  call hash_table_remove_from_ptr(table, ptr)
  if (enable_cache) call hash_table_reset_cache(table)
end subroutine hash_table_get_remove_block


!****s* hash_tables/hash_table_get_keywords ===================================!
! NAME                                                                         !
!   hash_table_get_keywords (PUBLIC)                                           !
!------------------------------------------------------------------------------!
! DESCRIPTION                                                                  !
!   Returns to an array all of the keywords that are associated with values    !
!   (single line or blocks) in the table.                                      !
!------------------------------------------------------------------------------!
! ARGUMENTS                                                                    !
!   type(hash_table),               intent(in)  ::  table                      !
!   character(len=*), dimension(:), intent(out) ::  keywords                   !
!------------------------------------------------------------------------------!
! NOTES                                                                        !
!   Table must be initialized.                                                 !
!                                                                              !
!   Output array must be allocated to a size large enough to store all of the  !
!   keywords in the table. (Assume hash_table_size has been called).           !
!==============================================================================!
subroutine hash_table_get_keywords(table, keywords)
  implicit none
  type(hash_table),               intent(in)  ::  table
  character(len=*), dimension(:), intent(out) ::  keywords
  ! local vars
  type(hash_element), pointer   ::  ptr
  integer                       ::  ibucket, ikeyword

  if (.not. table%initialized) stop 'Error in hash_table_get_keywords: table not initialized'

  if (size(keywords,1) .lt. table%nvalues) &
  & stop 'Error in hash_table_get_keywords: keyword array too small'

  if (len(keywords) .lt. str_len) &
  & stop 'Error in hash_table_get_keywords: keyword array character length too small'

  ikeyword = 1

  if (table%nvalues .gt. 0) then
    do ibucket = 1, table%nbuckets
      ! if bucket in use..
      if (associated(table%buckets(ibucket)%head)) then
        ptr => table%buckets(ibucket)%head

        ! loop over elements in bucket (until end of list)
        do
          keywords(ikeyword) = ptr%keyword
          ikeyword = ikeyword + 1
        
          if (.not. associated(ptr%next_element)) exit
          ptr => ptr%next_element
        end do
      end if
    end do
  end if
end subroutine hash_table_get_keywords


!****s* hash_tables/hash_table_resize =========================================!
! NAME                                                                         !
!   hash_table_resize (PUBLIC)                                                 !
!------------------------------------------------------------------------------!
! DESCRIPTION                                                                  !
!   Shrinks or expands a hash table (used if the enable_expand or              !
!   enable_shrink flags are set) by temporarily creating a new hash table of   !
!   same size as the existing one, copying each element into this temporary    !
!   table, resizing the initial table, copying all of the elements back into   !
!   the original table and then destroying the temporary table.                !
!------------------------------------------------------------------------------!
! ARGUMENTS                                                                    !
!   type(hash_table), intent(inout) ::  table                                  !
!   real,             intent(in)    ::  scale_factor                           !
!------------------------------------------------------------------------------!
! NOTES                                                                        !
!   Routine sufficiently low level that we assume the table is initialized.    !
!   It is expected that the end user won't call this routine themselves and    !
!   will instead make use of automatic table resizing. Hence we don't check    !
!   for table initialization.                                                  !
!==============================================================================!
subroutine hash_table_resize(table, scale_factor)
  implicit none
  type(hash_table), intent(inout) ::  table
  real,             intent(in)    ::  scale_factor
  ! local vars:
  type(hash_table)            ::  tmp_table
  type(hash_element), pointer ::  ptr
  integer ::  ibucket, istat

  ! create new table with same number of buckets
  call hash_table_init(tmp_table, table%nbuckets)

  ! copy pointers to head and tail for each bucket
  do ibucket = 1, table%nbuckets
    if (associated(table%buckets(ibucket)%head)) then
      tmp_table%buckets(ibucket)%head => table%buckets(ibucket)%head
      nullify(table%buckets(ibucket)%head)
    end if

    if (associated(table%buckets(ibucket)%tail)) then
      tmp_table%buckets(ibucket)%tail => table%buckets(ibucket)%tail
      nullify(table%buckets(ibucket)%tail)
    end if
  end do

  ! deallocate and reinitialize with new number of buckets
  deallocate(table%buckets, stat=istat)
  if (istat .ne. 0) stop 'Error in hash_table_resize: cannot deallocate buckets array'

  table%initialized = .false.   ! to allow for reinitialization without calling hash_table_destroy
  call hash_table_init(table, int(tmp_table%nbuckets*scale_factor))

  nullify(ptr)

  ! loop through old (tmp table), for each item in list, create new elements for new table
  do ibucket = 1, tmp_table%nbuckets

    ! if list in use
    if (associated(tmp_table%buckets(ibucket)%head)) then
      ptr => tmp_table%buckets(ibucket)%head

      ! loop through each element and create copy for new table
      do
        if (allocated(ptr%block)) then
          call hash_table_add_block(table, ptr%keyword, ptr%block)
        else
          call hash_table_add_value(table, ptr%keyword, ptr%value)
        end if

        ! exit at end of list
        if (.not. associated(ptr%next_element)) exit
        ptr => ptr%next_element
      end do
    end if

  end do

  ! just for safety, reset the cache of the new table (although it should be fine)
  if (enable_cache) call hash_table_reset_cache(table)

  ! free up memory and recalculate load factor for new table
  call hash_table_destroy(tmp_table)
  table%load_factor = real(table%nvalues)/real(table%nbuckets)
end subroutine hash_table_resize


!****s* hash_tables/hash_table_list ===========================================!
! NAME                                                                         !
!   hash_table_list (PUBLIC)                                                   !
!------------------------------------------------------------------------------!
! DESCRIPTION                                                                  !
!   Outputs (most of) the content of the hash table. Used for debugging.       !
!------------------------------------------------------------------------------!
! ARGUMENTS                                                                    !
!   type(hash_table),   intent(in)  ::  table                                  !
!                                                                              !
!   integer,  optional, intent(in)  ::  unit_num                               !
!------------------------------------------------------------------------------!
! NOTES                                                                        !
!   Table must be initialized.                                                 !
!                                                                              !
!   If unit_num is specified, can write to this unit - else, write to stdout.  !
!   Note that unit_num must be open/writable etc.. not done in here.           !
!                                                                              !
!   Formatting set up for an str_len of 100.. larger than this will look bad.  !
!==============================================================================!
subroutine hash_table_list(table, unit_num)
  use iso_fortran_env,  only: output_unit
  implicit none
  type(hash_table),   intent(in)  ::  table
  integer,  optional, intent(in)  ::  unit_num
  ! local vars
  type(hash_element), pointer   ::  ptr
  integer                       ::  ibucket, iblock
  integer                       ::  out_unit

  if (.not. table%initialized) stop 'Error in hash_table_list: table not initialized'

  !+----------------------------------------------------------------------------------------------------------------------+
  !| hash table                                                                                                           |
  !+----------------------------------------------------------------------------------------------------------------------+
  !| general info :                                                                                                       |
  !|     nbuckets : xxxxxxxxxx  | nvalues   : xxxxxxxxxx | load_factor : x.xxxx | can_expand : x | can_shrink : x         |
  !|                                                                                                                      |
  !|     warning: str_len is larger than this routine is able to print (str_len > 100)                                    |
  !+----------------------------------------------------------------------------------------------------------------------+
  !|     hash_id  : xxxxxxxxxxx                          |  bucket_id  : xxxxxxxxxx                                       |
  !|     keyword  : xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx  |
  !|     value    : xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx  |
  !+----------------------------------------------------------------------------------------------------------------------+

  if (.not. present(unit_num)) then
    out_unit = output_unit
  else
    out_unit = unit_num
  end if

  write(out_unit, fmt=100)
  write(out_unit, "('| hash table', T120, '|')")
  write(out_unit, fmt=100)
  write(out_unit, "('| general info :', T120, '|')")

  write(out_unit, fmt=300) table%nbuckets, table%nvalues, table%load_factor, table%can_expand, table%can_shrink

  ! explain that we may truncate output.. (output to fit in window of 120 characters)
  if (str_len > 100) then
    write(out_unit, fmt=200)
    write(out_unit, "('|', T7, A, T120, '|')") "warning: str_len is larger than this routine is able to print (str_len > 100)"
  end if

  write(out_unit, fmt=100)

  do ibucket = 1, table%nbuckets
    ! if bucket in use..
    if (associated(table%buckets(ibucket)%head)) then
      ptr => table%buckets(ibucket)%head

      ! loop over elements in bucket (until end of list)
      do
        write(out_unit, "('|', T7, 'hash_id  : ', I11, T55, '|  bucket_id  : ', I10, T120, '|')") ptr%hash_id, ibucket
        write(out_unit, "('|', T7, 'keyword  : ', A100, T120, '|')") ptr%keyword

        if (allocated(ptr%block)) then
          do iblock = 1, size(ptr%block,1)
            if (iblock .eq. 1) then 
              write(out_unit, "('|', T7, 'block    : ', A100, T120, '|')") ptr%block(iblock)
            else
              write(out_unit, "('|', T18, A100, T120, '|')") ptr%block(iblock)
            end if
          end do
        else
          write(out_unit, "('|', T7, 'value    : ', A100, T120, '|')") ptr%value
        end if

        write(out_unit, fmt=100)
        if (.not. associated(ptr%next_element)) exit
        ptr => ptr%next_element
      end do
    end if
  end do

100 format ('+', 118('-'), '+') ! top border
200 format ('|', T120, '|')     ! empty line
300 format &
& ('|     nbuckets : ',I10,'  | nvalues   : ',I10,' | load_factor : ',F6.4,' | can_expand :',L2,' | can_shrink :',L2,T120,'|')
end subroutine hash_table_list


!****s* hash_tables/hash_table_destroy ========================================!
! NAME                                                                         !
!   hash_table_destroy (PUBLIC)                                                !
!------------------------------------------------------------------------------!
! DESCRIPTION                                                                  !
!   Destroys a hash table, freeing up memory and allowing it to be reused.     !
!------------------------------------------------------------------------------!
! ARGUMENTS                                                                    !
!   type(hash_table), intent(inout) ::  table                                  !
!------------------------------------------------------------------------------!
! NOTES                                                                        !
!   Table must be initialized.                                                 !
!==============================================================================!
subroutine hash_table_destroy(table)
  implicit none
  type(hash_table), intent(inout) ::  table
  ! local variables:
  type(hash_element), pointer     ::  ptr, tmp_ptr
  integer                         ::  ibucket, istat

  ! for safety this routine nullifies every pointer even if the element is deallocated
  ! should hopefully be left with no dangling pointers anywhere..

  if (.not. table%initialized) stop 'Error in hash_table_destroy: table not initialized'

  nullify(ptr)

  do ibucket = 1, table%nbuckets

    ! if bucket in use..
    if (associated(table%buckets(ibucket)%tail)) then
      ptr => table%buckets(ibucket)%tail

      ! loop over all elements in bucket (from tail -> head)
      do
        ! point to previous element (if we can)
        nullify(tmp_ptr)
        if (associated(ptr%prev_element)) tmp_ptr => ptr%prev_element

        ! remove element
        nullify(ptr%next_element, ptr%prev_element)
        if (allocated(ptr%block)) then
          deallocate(ptr%block, stat=istat)
          if (istat .ne. 0) stop 'Error in hash_table_destroy: could not deallocate block'
        end if
        deallocate(ptr, stat=istat)
        if (istat .ne. 0) stop 'Error in hash_table_destroy: could not deallocate element'

        ! exit if no previous element, else move backwards
        if (.not. associated(tmp_ptr)) exit
        ptr => tmp_ptr
      end do
    end if

    nullify(table%buckets(ibucket)%head, table%buckets(ibucket)%tail)
  end do

  deallocate(table%buckets, stat=istat)
  if (istat .ne. 0) stop 'Error in hash_table_destroy: could not deallocate buckets'

  table%nvalues = 0
  table%nbuckets = 0
  table%load_factor = 0.0
  table%initialized = .false.
  if (enable_cache) call hash_table_reset_cache(table)
end subroutine hash_table_destroy

!==============================================================================!
! PRIVATE ROUTINES BELOW HERE:                                                 !
!==============================================================================!


!****s* hash_tables/hash_table_find_element ===================================!
! NAME                                                                         !
!   hash_table_find_element (PRIVATE)                                          !
!------------------------------------------------------------------------------!
! DESCRIPTION                                                                  !
!   Returns a pointer to an element of a hash table associated with the        !
!   keyword. If keyword does not exist in the table, then pointer is null.     !
!                                                                              !
!   The hash_id and bucket_id corresponding to this element are also returned. !
!------------------------------------------------------------------------------!
! ARGUMENTS                                                                    !
!   type(hash_table),   intent(inout) ::  table                                !
!   character(len=*),   intent(in)    ::  keyword                              !
!   integer,            intent(out)   ::  hash_id                              !
!   integer,            intent(out)   ::  bucket_id                            !
!   type(hash_element), pointer       ::  ptr                                  !
!------------------------------------------------------------------------------!
! NOTES                                                                        !
!   Since this is a low level routine (returns a pointer) we assume that the   !
!   table is initialized and everything is allocated as expected.              !
!==============================================================================!
subroutine hash_table_find_element(table, keyword, hash_id, bucket_id, ptr)
  implicit none
  type(hash_table),   intent(inout) ::  table
  character(len=*),   intent(in)    ::  keyword
  integer,            intent(out)   ::  hash_id
  integer,            intent(out)   ::  bucket_id
  type(hash_element), pointer       ::  ptr

  ! Returns pointer to element containing keyword if it exists - null otherwise (or bucket not in use)
  ! hash_id and bucket_id are also returned

  ! check the table can fit in the keyword..
  if (len_trim(keyword) .gt. str_len) stop 'Error in hash_table_find_element: keyword too long'

  hash_id = hash_function(trim(keyword))
  bucket_id = modulo(hash_id, table%nbuckets) + 1

  nullify(ptr)

  ! check to see if cache already has correct element..
  if (enable_cache) then
    if (associated(table%last_access) .and. (table%last_keyword .eq. keyword)) then
      ! sanity check:
      if (hash_id .ne. table%last_access%hash_id) stop 'Error in hash_table_find_element: contradictory hash_id'

      ! last_access is the correct element
      ptr => table%last_access
    end if
  end if

  ! if we used cache value, then won't go into here..
  if (.not. associated(ptr)) then
    ! check to see if bucket is in use
    if (associated(table%buckets(bucket_id)%head)) then

      ! start at head of bucket (must be associated)
      ptr => table%buckets(bucket_id)%head
      
      ! go through list until we have matching keyword or reach end
      do
        if (ptr%keyword .eq. keyword) then
          ! sanity check:
          if (hash_id .ne. ptr%hash_id) stop 'Error in hash_table_find_element: contradictory hash_id'
          exit
        end if

        ! won't go here if we found the keyword..
        if (.not. associated(ptr%next_element)) then
          nullify(ptr)
          ! ran out of list and didn't find keyword..
          exit
        end if
        ptr => ptr%next_element
      end do
    end if
    ! no need to nullify if head not associated, because we have done no work, so ptr still null

    ! update cache
    if (enable_cache) call hash_table_update_cache(table, ptr)
  end if
end subroutine hash_table_find_element


!****s* hash_tables/hash_table_create_element_value ===========================!
! NAME                                                                         !
!   hash_table_create_element_value (PRIVATE)                                  !
!------------------------------------------------------------------------------!
! DESCRIPTION                                                                  !
!   Creates an element (for use in a hash table) with the keyword, value       !
!   which is a single line only) and hash_id stored.                           !
!------------------------------------------------------------------------------!
! ARGUMENTS                                                                    !
!   character(len=*), intent(in)  ::  keyword                                  !
!   character(len=*), intent(in)  ::  value                                    !
!   integer,          intent(in)  ::  hash_id                                  !
!   type(hash_element), pointer   ::  element                                  !
!==============================================================================!
subroutine hash_table_create_element_value(keyword, value, hash_id, element)
  implicit none
  character(len=*), intent(in)  ::  keyword
  character(len=*), intent(in)  ::  value
  integer,          intent(in)  ::  hash_id
  type(hash_element), pointer   ::  element
  integer ::  istat

  ! store new value in memory
  allocate(element, stat=istat)
  if (istat .ne. 0) stop 'Error in hash_table_create_element_value: could not allocate element'
  element%hash_id = hash_id
  element%keyword = keyword
  element%value = value
end subroutine hash_table_create_element_value


!****s* hash_tables/hash_table_create_element_block ===========================!
! NAME                                                                         !
!   hash_table_create_element_block (PRIVATE)                                  !
!------------------------------------------------------------------------------!
! DESCRIPTION                                                                  !
!   Creates an element (for use in a hash table) with the keyword, block and   !
!   hash_id stored.                                                            !
!------------------------------------------------------------------------------!
! ARGUMENTS                                                                    !
!   character(len=*),               intent(in)  ::  keyword                    !
!   character(len=*), dimension(:), intent(in)  ::  block                      !
!   integer,                        intent(in)  ::  hash_id                    !
!   type(hash_element), pointer                 ::  element                    !
!==============================================================================!
subroutine hash_table_create_element_block(keyword, block, hash_id, element)
  implicit none
  character(len=*),               intent(in)  ::  keyword
  character(len=*), dimension(:), intent(in)  ::  block
  integer,                        intent(in)  ::  hash_id
  type(hash_element), pointer                 ::  element
  integer ::  istat

  ! creates block regardless of size - if block is one line long
  ! - we should use hash_table_create_element_value instead

  ! store new value in memory
  allocate(element, stat=istat)
  if (istat .ne. 0) stop 'Error in hash_table_create_element_block: could not allocate element'
  element%hash_id = hash_id
  element%keyword = keyword
  allocate(element%block(size(block,1)), stat=istat)
  if (istat .ne. 0) stop 'Error in hash_table_create_element block: could not allocate block'
  element%block(:) = block(:)
  element%value = ' '
end subroutine hash_table_create_element_block


!****s* hash_tables/hash_table_get_from_ptr_value =============================!
! NAME                                                                         !
!   hash_table_get_from_ptr_value (PRIVATE)                                    !
!------------------------------------------------------------------------------!
! DESCRIPTION                                                                  !
!   Returns the value (single line) associated with a pointer (which points    !
!   to an element of a hash table). If element_index is included, this can     !
!   be used to return a single line of a block.                                !
!------------------------------------------------------------------------------!
! ARGUMENTS                                                                    !
!   type(hash_element),     pointer     ::  ptr                                !
!   character(len=str_len), intent(out) ::  value                              !
!                                                                              !
!   integer,  optional,     intent(in)  ::  element_index                      !
!------------------------------------------------------------------------------!
! NOTES                                                                        !
!   If ptr is null, it is assumed that the keyword does not exist in the       !
!   table. (And therefore flags an error.)                                     !
!                                                                              !
!   If element_index is included and is > 1, then it is assumed that the block !
!   contains at least this many values (else we flag an error).                !
!   (Assume that the result of hash_table_query is known.)                     !
!==============================================================================!
subroutine hash_table_get_from_ptr_value(ptr, value, element_index)
  implicit none
  type(hash_element),     pointer     ::  ptr
  character(len=str_len), intent(out) ::  value
  integer,  optional,     intent(in)  ::  element_index

  ! if we don't have null pointer, then keyword must match at this stage
  if (.not. associated(ptr)) stop 'Error in hash_table_get_from_ptr_value: ptr not associated - keyword does not exit'

  if (present(element_index)) then
    select case (element_index)
      case (:0) ! element_index <= 0
        stop 'Error in hash_table_get_from_ptr_value: element_index must be >= 1'
      case (1)  ! element_index == 1
        ! could be either block or scalar
        if (allocated(ptr%block)) then
          value = ptr%block(1)
        else
          value = ptr%value
        end if
      case (2:) ! element_index >= 2
        ! must be inside block in this case
        if (.not. allocated(ptr%block)) then
          stop 'Error in hash_table_get_from_ptr_value: block not allocated'
        else
          if (size(ptr%block,1) .lt. element_index) then
            stop 'Error in hash_table_get_from_ptr_value: could not get requested element from block'
          end if
          value = ptr%block(element_index)
        end if
    end select
  else
    ! should just be a scalar value - if block, give error as we can't return this to a scalar
    if (allocated(ptr%block)) &
      & stop 'Error in hash_table_get_from_ptr_value: cannot return block to a scalar - consider getting with element_index'
    value = ptr%value
  end if
end subroutine hash_table_get_from_ptr_value


!****s* hash_tables/hash_table_get_from_ptr_block =============================!
! NAME                                                                         !
!   hash_table_get_from_ptr_block (PRIVATE)                                    !
!------------------------------------------------------------------------------!
! DESCRIPTION                                                                  !
!   Returns the block associated with a pointer (which points to an element of !
!   a hash table).                                                             !
!------------------------------------------------------------------------------!
! ARGUMENTS                                                                    !
!   type(hash_element),     pointer                   ::  ptr                  !
!   character(len=str_len), dimension(:), intent(out) ::  block                !
!------------------------------------------------------------------------------!
! NOTES                                                                        !
!   If ptr is null, it is assumed that the keyword does not exist in the       !
!   table. (And we flag an error.)                                             !
!                                                                              !
!   Output array must be large enough to store entire block. (Can be larger).  !
!                                                                              !
!   If keyword is associated with a single value (not a block), this is still  !
!   able to be returned to an array.                                           !
!==============================================================================!
subroutine hash_table_get_from_ptr_block(ptr, block)
  implicit none
  type(hash_element),     pointer                   ::  ptr
  character(len=str_len), dimension(:), intent(out) ::  block
  integer ::  iblock

  ! if we don't have null pointer, then keyword must match at this stage
  if (.not. associated(ptr)) stop 'Error in hash_table_get_from_ptr_block: ptr not associated - keyword does not exist'

  if (size(block,1) .lt. 1) stop 'Error in hash_table_get_from_ptr_block: block array not large enough'

  if (allocated(ptr%block)) then
    ! block allocated - must be here - so copy everything
    if (size(block,1) .lt. size(ptr%block,1)) stop 'Error in hash_table_get_from_ptr_block: block array not large enough'

    ! ptr%block can be smaller than block
    do iblock = 1, size(ptr%block,1)
      block(iblock) = ptr%block(iblock)
    end do
  else
    ! must be a single value - input block might have some empty values in this case...
    block(1) = ptr%value
  end if
end subroutine hash_table_get_from_ptr_block


!****s* hash_tables/hash_table_remove_from_ptr ================================!
! NAME                                                                         !
!   hash_table_remove_from_ptr (PRIVATE)                                       !
!------------------------------------------------------------------------------!
! DESCRIPTION                                                                  !
!   Removes the element of the hash table associated with the pointer (which   !
!   points to an element of the hash table), while preserving the integrity    !
!   of the associated linked list.                                             !
!------------------------------------------------------------------------------!
! ARGUMENTS                                                                    !
!   type(hash_table),   intent(inout) ::  table                                !
!   type(hash_element), pointer       ::  ptr                                  !
!------------------------------------------------------------------------------!
! NOTES                                                                        !
!   If ptr is null, it is assumed that the keyword does not exist in the       !
!   table. (We get an error.)                                                  !
!==============================================================================!
subroutine hash_table_remove_from_ptr(table, ptr)
  implicit none
  type(hash_table),   intent(inout) ::  table
  type(hash_element), pointer       ::  ptr
  ! local vars:
  type(hash_element), pointer       ::  tmp_ptr
  logical ::  is_head, is_tail
  integer ::  bucket_id, istat

  if (.not. associated(ptr)) stop 'Error in hash_table_remove_from_ptr: ptr not associated - keyword does not exist'

  bucket_id = modulo(ptr%hash_id, table%nbuckets) + 1

  nullify(tmp_ptr)

  is_head = (.not. associated(ptr%prev_element))
  is_tail = (.not. associated(ptr%next_element))

  if (is_head) then

    if (is_tail) then
      ! have nowhere to go, nullify head and tail
      nullify(table%buckets(bucket_id)%head)
      nullify(table%buckets(bucket_id)%tail)
    else
      ! is head, not tail - can move head to next element along
      tmp_ptr => ptr%next_element                 ! tmp_ptr is to right
      table%buckets(bucket_id)%head => tmp_ptr    ! move head to right
      nullify(tmp_ptr%prev_element)               ! make right/new head's pointer to left null
    end if

  else  ! not head

    if (is_tail) then
      ! not head, is tail
      tmp_ptr => ptr%prev_element               ! tmp_ptr is to left
      table%buckets(bucket_id)%tail => tmp_ptr  ! move tail to left
      nullify(tmp_ptr%next_element)             ! make left/new tail's pointer to right null
    else
      ! neither head nor tail - we need to move around a bit here

      ! move to left
      tmp_ptr => ptr%prev_element
      tmp_ptr%next_element => ptr%next_element  ! left's pointer to right points to right of original

      ! then to right (from original pointer)
      tmp_ptr => ptr%next_element
      tmp_ptr%prev_element => ptr%prev_element  ! right's pointer to left points to left or original
    end if
  end if

  nullify(tmp_ptr)  ! might as well

  ! should be able to remove the element now..
  nullify(ptr%prev_element, ptr%next_element)
  if (allocated(ptr%block)) then
    deallocate(ptr%block, stat=istat)
    if (istat .ne. 0) stop 'Error in hash_table_remove_from_ptr: error removing element - could not deallocate block'
  end if
  deallocate(ptr, stat=istat)
  if (istat .ne. 0) stop 'Error in hash_table_remove_from_ptr: could not remove element from table'

  table%nvalues = table%nvalues - 1
  table%load_factor = real(table%nvalues)/real(table%nbuckets)

  if (table%can_shrink .and. (table%load_factor .lt. min_load_factor)) call hash_table_resize(table, shrink_scale)
end subroutine hash_table_remove_from_ptr


!****s* hash_tables/hash_table_update_cache ===================================!
! NAME                                                                         !
!   hash_table_update_cache (PRIVATE)                                          !
!------------------------------------------------------------------------------!
! DESCRIPTION                                                                  !
!   Updates the cache of the hash table to store the keyword of the last       !
!   element accessed as well as a pointer to the element.                      !
!                                                                              !
!   This reduces the cost of repeated access to the same element (eg: a query  !
!   followed by a get/remove - especially useful for blocks).                  !
!------------------------------------------------------------------------------!
! ARGUMENTS                                                                    !
!   type(hash_table),   intent(inout) ::  table                                !
!   type(hash_element), pointer       ::  last_access                          !
!------------------------------------------------------------------------------!
! NOTES                                                                        !
!   Table must be initialized (note that we don't check)                       !
!                                                                              !
!   Only used if module level enable_cache is .TRUE.                           !
!==============================================================================!
subroutine hash_table_update_cache(table, last_access)
  implicit none
  type(hash_table),   intent(inout) ::  table
  type(hash_element), pointer       ::  last_access

  if (associated(last_access)) then
    table%last_access => last_access
    table%last_keyword = last_access%keyword
  else
    call hash_table_reset_cache(table)
  end if
end subroutine hash_table_update_cache


!****s* hash_tables/hash_table_reset_cache ====================================!
! NAME                                                                         !
!   hash_table_reset_cache (PRIVATE)                                           !
!------------------------------------------------------------------------------!
! DESCRIPTION                                                                  !
!   Clears the cache of a hash table.                                          !
!------------------------------------------------------------------------------!
! ARGUMENTS                                                                    !
!   type(hash_table), intent(inout) ::  table                                  !
!------------------------------------------------------------------------------!
! NOTES                                                                        !
!   Table must be initialized (note that we don't check)                       !
!                                                                              !
!   Only used if module level enable_cache is .TRUE.                           !
!==============================================================================!
subroutine hash_table_reset_cache(table)
  implicit none
  type(hash_table), intent(inout) ::  table

  nullify(table%last_access)
  table%last_keyword = ' '
end subroutine hash_table_reset_cache


!****f* hash_tables/hash_function =============================================!
! NAME                                                                         !
!   hash_function (PRIVATE)                                                    !
!------------------------------------------------------------------------------!
! DESCRIPTION                                                                  !
!   Hashes a string - basically a wrapper for fnv_1a_hash and all of the       !
!   conversions between 64 bit, 32 bit and default (probably 32 bit) integers. !
!------------------------------------------------------------------------------!
! ARGUMENTS                                                                    !
!   character(len=*), intent(in)  ::  str                                      !
!==============================================================================!
function hash_function(str)
  implicit none
  character(len=*), intent(in)  ::  str
  integer(kind=i32)             ::  tmp_hash
  integer                       ::  hash_function ! default kind

  ! this function is basically just a wrapper for fnv_1a_hash
  tmp_hash = i64_to_i32(fnv_1a_hash(str))
  hash_function = int(tmp_hash)           ! probably pointless - cast to default kind
end function hash_function


!****f* hash_tables/fnv_1a_hash ===============================================!
! NAME                                                                         !
!   fnv_1a_hash (PRIVATE)                                                      !
!------------------------------------------------------------------------------!
! DESCRIPTION                                                                  !
!   Hashes the input string using the FNV 1a algorithm, see:                   !
!   http://www.isthe.com/chongo/tech/comp/fnv/                                 !
!   for more details.                                                          !
!                                                                              !
!   This relies on the use of 32 bit unsigned integers - to get around this    !
!   we use a 64 bit int and then convert back (using only the first 32 bits).  !
!------------------------------------------------------------------------------!
! ARGUMENTS                                                                    !
!   character(len=*), intent(in)  ::  str                                      !
!==============================================================================!
function fnv_1a_hash(str)
  implicit none
  character(len=*), intent(in)  ::  str
  integer(kind=i64),  parameter ::  offset_basis = 2166136261_i64
  integer(kind=i64),  parameter ::  fnv_prime = 16777619_i64
  integer(kind=i64)             ::  octet
  integer(kind=i64)             ::  fnv_1a_hash
  integer ::  ioctet

  ! returns a 32-bit hash even though we use the i64 kind..
  ! we treat the ASCII code for each character of the string as an octet
  ! rather than building a very large integer
  fnv_1a_hash = offset_basis
  do ioctet = 1, len(str)
    octet = iachar(str(ioctet:ioctet), kind=i64)
    fnv_1a_hash = ieor(fnv_1a_hash, octet)
    fnv_1a_hash = fnv_1a_hash * fnv_prime
  end do
end function fnv_1a_hash


!****f* hash_tables/i64_to_i32 ================================================!
! NAME                                                                         !
!   i64_to_i32 (PRIVATE)                                                       !
!------------------------------------------------------------------------------!
! DESCRIPTION                                                                  !
!   Copies the first 32 bits of a 64 bit integer into a 32 bit integer.        !
!   Used to get around Fortran not having an unsigned int data type.           !
!------------------------------------------------------------------------------!
! ARGUMENTS                                                                    !
!   integer(kind=i64),  intent(in)  ::  x                                      !
!------------------------------------------------------------------------------!
! NOTES                                                                        !
!   Relies on the kinds i64 and i32 being at least 64 and 32 bits.             !
!                                                                              !
!   Because this copies the *first* 32 bits of the 64 bit number into the 32   !
!   bit number, this might be incorrect if we should instead be copying the    !
!   *last* 32 bits. For our purposes this should still be fine though.         !
!==============================================================================!
function i64_to_i32(x)
  implicit none
  integer(kind=i64),  intent(in)  ::  x
  integer(kind=i32)               ::  i64_to_i32
  integer :: ibit                 ! kind shouldn't matter, as long as it can hold 0->31

  do ibit = 0, bit_size(i64_to_i32)-1
    if (btest(x,ibit)) then
      i64_to_i32 = ibset(i64_to_i32,ibit)
    else
      ! could probably set i64_to_i32 = 0 then change only the bits that are 1 in x
      ! but this is probably safer (albeit slightly slower)
      i64_to_i32 = ibclr(i64_to_i32,ibit)
    end if
  end do
end function i64_to_i32
end module hash_tables
