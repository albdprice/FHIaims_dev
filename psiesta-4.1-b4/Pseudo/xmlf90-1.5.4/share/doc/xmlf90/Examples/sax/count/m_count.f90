MODULE m_count
!
! Contributed by Jon Wakelin
!

  USE xmlf90_sax

  IMPLICIT NONE

  private 
  PUBLIC :: begin_element_handler, end_element_handler, pcdata_chunk_handler

  TYPE, public :: hash
     CHARACTER(len=50) :: elm
     INTEGER           :: num
  END TYPE hash

  TYPE(hash), DIMENSION(50), public  :: element_hash
  INTEGER, public, save              ::  nhash = 0

  INTEGER, private, save             :: n = 0 

!--------------------------------------------------------
CONTAINS

  SUBROUTINE begin_element_handler(name,attributes)
    character(len=*), intent(in)   :: name
    TYPE(dictionary_t), INTENT(in) :: attributes

    LOGICAL :: match  
    INTEGER :: pmatch
    INTEGER :: i

    match = .false.

!!! First time through loop element must be unique...
    IF (n == 0) THEN
       element_hash(n+1)%elm = name
       element_hash(n+1)%num = 1
       nhash=nhash+1
    ELSE

!!! ...thereafter we will have to check if it is unique
       DO i=1,nhash
          IF (name == element_hash(i)%elm) THEN
             match = .true. ! set .true. if element already exists
             pmatch = i     ! and record the position at which the match occured
                            ! NB there can only ever be 1 or 0 matches
          ENDIF
       ENDDO

!!! If element already exists increment the counter for THIS element
       IF (match) THEN
          element_hash(pmatch)%num = element_hash(pmatch)%num + 1
       ELSE
!!! Otherwise make a new entry in the hash 
          element_hash(n+1)%elm = name
          element_hash(n+1)%num = 1
          nhash=nhash+1
       ENDIF
    ENDIF
    n=nhash

  END SUBROUTINE begin_element_handler


!--------------------------------------------------------------------------
! End tag handler
  SUBROUTINE end_element_handler(name)
    character(len=*), intent(in)     :: name
  END SUBROUTINE end_element_handler
  
  ! PCDATA handler
  SUBROUTINE pcdata_chunk_handler(chunk)
    CHARACTER(len=*), INTENT(in) :: chunk
  END SUBROUTINE pcdata_chunk_handler
  
END MODULE m_count
