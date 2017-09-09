! Module implementing an OO hash table (dictionary) in Fortran 2003.
! Compiles and runs with accompanying test program under the Intel 
! Fortran Compiler, version 11.1.046

! Copyright (c) Izaak Beekman 2010

    ! This program is free software: you can redistribute it and/or modify
    ! it under the terms of the GNU Lesser General Public License as published by
    ! the Free Software Foundation, either version 3 of the License, or
    ! (at your option) any later version.

    ! This program is distributed in the hope that it will be useful,
    ! but WITHOUT ANY WARRANTY; without even the implied warranty of
    ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    ! GNU Lesser General Public License for more details.

    ! You should have received a copy of the GNU Lesser General Public License
    ! along with this program.  If not, see <http://www.gnu.org/licenses/>.

MODULE hashtbl
  USE String_Mod
  IMPLICIT NONE ! Use strong typing
  INTEGER, PARAMETER :: tbl_size=50

  TYPE sllist
     TYPE(sllist), POINTER  :: child => NULL()
     CHARACTER, ALLOCATABLE :: key(:)
     INTEGER                :: val=0
   CONTAINS
     PROCEDURE :: put  => put_sll
     PROCEDURE :: get  => get_sll
     PROCEDURE :: free => free_sll
  END TYPE sllist

  TYPE hash_tbl_sll
     TYPE(sllist), ALLOCATABLE :: vec(:)
     INTEGER                   :: vec_len = 0
     LOGICAL                   :: is_init = .FALSE.
   CONTAINS
     PROCEDURE :: init => init_hash_tbl_sll
     PROCEDURE :: put  => put_hash_tbl_sll
     PROCEDURE :: get  => get_hash_tbl_sll
     PROCEDURE :: free => free_hash_tbl_sll
  END TYPE hash_tbl_sll

  PUBLIC :: hash_tbl_sll
  INTEGER, PARAMETER :: tbl_length = 100
  !TYPE(hash_tbl_sll), PRIVATE :: table
  
  CONTAINS

  RECURSIVE SUBROUTINE put_sll(list,key,val)
    CLASS(sllist), INTENT(inout) :: list
    CHARACTER(*),  INTENT(in)    :: key
    INTEGER,       INTENT(in)    :: val
    INTEGER                      :: keylen
    INTEGER   :: i
    CHARACTER :: keyVec(LEN_TRIM(key))

    keylen = LEN_TRIM(key)
    FORALL (i=1:LEN_TRIM(key))
       keyVec(i) = key(i:i)
    END FORALL
    IF (ALLOCATED(list%key)) THEN
       IF ( (list%key .NES. keyVec) .OR. SIZE(list%key)/=SIZE(keyVec)) THEN
          IF ( .NOT. ASSOCIATED(list%child) ) ALLOCATE(list%child)
          CALL put_sll(list%child,key,val)
       ELSE   
         list%val = val
       END IF
    ELSE
       IF (.NOT. ALLOCATED(list%key)) ALLOCATE(list%key(keylen))
       list%key = keyVec
       list%val = val
    END IF
  END SUBROUTINE put_sll


  RECURSIVE SUBROUTINE get_sll(list,key,val)
    CLASS(sllist),                 INTENT(in)    :: list
    CHARACTER(len=*),              INTENT(in)    :: key
    INTEGER, INTENT(out)   :: val
    INTEGER :: i
    CHARACTER :: keyVec(LEN(key))

    FORALL (i=1:LEN(key))
       keyVec(i) = key(i:i)
    END FORALL

    IF (ALLOCATED(list%key) .AND. SIZE(list%key)==SIZE(keyVec) .AND. (list%key .EQS. keyVec)) THEN
       val = list%val
    ELSE IF(ASSOCIATED(list%child)) THEN ! keep going
       CALL get_sll(list%child,key,val)
    ELSE ! At the end of the list, no key found
       val=0
       RETURN
    END IF
  END SUBROUTINE get_sll


  RECURSIVE SUBROUTINE free_sll(list)
    CLASS(sllist), INTENT(inout) :: list
    IF (ASSOCIATED(list%child)) THEN
       CALL free_sll(list%child)
       DEALLOCATE(list%child)
    END IF
    list%child => NULL()
    IF (ALLOCATED(list%key)) DEALLOCATE(list%key)
  END SUBROUTINE free_sll
  
  SUBROUTINE init_hash_tbl_sll(tbl,tbl_len)
    CLASS(hash_tbl_sll),   INTENT(inout) :: tbl
    INTEGER,     OPTIONAL, INTENT(in)    :: tbl_len

    IF (ALLOCATED(tbl%vec)) DEALLOCATE(tbl%vec)
    IF (PRESENT(tbl_len)) THEN
       ALLOCATE(tbl%vec(0:tbl_len-1))
       tbl%vec_len = tbl_len
    ELSE
       ALLOCATE(tbl%vec(0:tbl_size-1))
       tbl%vec_len = tbl_size
    END IF
    tbl%is_init = .TRUE.
  END SUBROUTINE init_hash_tbl_sll

  ! The first part of the hashing procedure using the string
  ! collating sequence
  ELEMENTAL FUNCTION sum_string(str) RESULT(sig)
    CHARACTER(len=*), INTENT(in)   :: str
    INTEGER                        :: sig
    CHARACTER, DIMENSION(LEN(str)) :: tmp
    INTEGER :: i

    FORALL (i=1:LEN(str))
       tmp(i) = str(i:i)
    END FORALL
    sig = SUM(ICHAR(tmp))
  END FUNCTION sum_string


  SUBROUTINE put_hash_tbl_sll(tbl,key,val)
    CLASS(hash_tbl_sll), INTENT(inout) :: tbl
    CHARACTER(len=*),    INTENT(in)    :: key
    INTEGER,    INTENT(in)    :: val
    INTEGER                            :: hash

    hash = MOD(sum_string(key),tbl%vec_len)
    CALL tbl%vec(hash)%put(key=key,val=val)
  END SUBROUTINE put_hash_tbl_sll


  SUBROUTINE get_hash_tbl_sll(tbl,key,val)
    CLASS(hash_tbl_sll), INTENT(in)  :: tbl
    CHARACTER(len=*),    INTENT(in)  :: key
    INTEGER,             INTENT(out) :: val
    INTEGER                          :: hash

    hash = MOD(sum_string(key),tbl%vec_len)
    CALL tbl%vec(hash)%get(key=key,val=val)
  END SUBROUTINE get_hash_tbl_sll


  SUBROUTINE free_hash_tbl_sll(tbl)
    CLASS(hash_tbl_sll), INTENT(inout) :: tbl    
    INTEGER     :: i, low, high

    low  = LBOUND(tbl%vec,dim=1)
    high = UBOUND(tbl%vec,dim=1) 
    IF (ALLOCATED(tbl%vec)) THEN
       DO i=low,high
          CALL tbl%vec(i)%free()
       END DO
       DEALLOCATE(tbl%vec)
    END IF
    tbl%is_init = .FALSE.
  END SUBROUTINE free_hash_tbl_sll

SUBROUTINE InsertHash(Table,Reactant,NumberSpecies)
  TYPE(hash_tbl_sll) :: Table
  CHARACTER(*) :: Reactant
  INTEGER :: NumberSpecies
  INTEGER :: out
  CALL table%get(key=Reactant,val=out)
  IF (out==0) THEN
    NumberSpecies=NumberSpecies+1
    CALL table%put(key=Reactant,val=NumberSpecies)
  END IF  
END SUBROUTINE InsertHash

FUNCTION GetHash(Table,Reactant)
  INTEGER :: GetHash
  TYPE(hash_tbl_sll) :: Table
  CHARACTER(*) :: Reactant
  CALL table%get(key=Reactant,val=GetHash)
END FUNCTION GetHash

SUBROUTINE RemoveBlank(String)
  CHARACTER(*) :: String

  INTEGER :: i,LenString,Pos
  LenString=LEN(String)
  Pos=1
  DO i=1,LenString
    IF (INDEX(String(Pos:Pos),' ')>0) THEN
      String(Pos:)=String(Pos+1:)
      Pos=Pos-1
    END IF
    Pos=Pos+1
  END DO   
END SUBROUTINE RemoveBlank

SUBROUTINE InitHashTable(Table,TableLength)
  TYPE(hash_tbl_sll) :: Table
  INTEGER :: TableLength
  CALL Table%init(TableLength)
END SUBROUTINE InitHashTable

FUNCTION NumElemHashTable(Table)
  INTEGER :: NumElemHashTable
  TYPE(hash_tbl_sll) :: Table

  INTEGER :: i
  TYPE(sllist), POINTER :: child => NULL()

  NumElemHashTable=0
  DO i=LBOUND(Table%vec,dim=1),UBOUND(Table%vec,dim=1)
    IF (ALLOCATED(table%vec(i)%key)) THEN 
      NumElemHashTable=NumElemHashTable+1
    END IF
    Child=>table%vec(i)%Child
    DO 
      IF (ASSOCIATED(Child)) THEN
        Child=>Child%Child
        NumElemHashTable=NumElemHashTable+1
      ELSE  
        EXIT
      END IF  
    END DO
  END DO
END FUNCTION NumElemHashTable

!SUBROUTINE ListToHashTable(List,Table)
!  CHARACTER(*) :: List(:)
!  TYPE(hash_tbl_sll) :: Table
!
!  INTEGER :: i
!  DO i=1,SIZE(List)
!    CALL Table%put(TRIM(ADJUSTL(List(i))),i)
!  END DO  
!END SUBROUTINE ListToHashTable

!SUBROUTINE HashTableToList(Table,List)
!  TYPE(hash_tbl_sll) :: Table
!  CHARACTER(*) :: List(:)
!
!  INTEGER :: i,j
!  TYPE(sllist), POINTER :: child => NULL()
!
!  DO i=LBOUND(table%vec,dim=1), UBOUND(table%vec,dim=1)
!    IF (ALLOCATED(table%vec(i)%key)) THEN 
!      DO j=1,SIZE(table%vec(i)%key)
!        List(table%vec(i)%Val)(j:j)=table%vec(i)%key(j)
!      END DO
!    END IF
!    Child=>table%vec(i)%Child
!    DO 
!      IF (ASSOCIATED(Child)) THEN
!        DO j=1,SIZE(Child%key)
!          List(Child%Val)(j:j)=Child%key(j)
!        END DO
!        Child=>Child%Child
!      ELSE  
!        EXIT
!      END IF  
!    END DO
!  END DO
!END SUBROUTINE HashTableToList

SUBROUTINE PrintHashTable(Table)
  TYPE(hash_tbl_sll) :: Table

  INTEGER :: i,Sum
  TYPE(sllist), POINTER :: child => NULL()

  sum = 0
  PRINT*, 'Indices of the hash table with content:'
  DO i=LBOUND(table%vec,dim=1), UBOUND(table%vec,dim=1)
    IF (ALLOCATED(table%vec(i)%key)) THEN 
      PRINT*, i,table%vec(i)%key,table%vec(i)%Val
      sum=sum+1
    END IF
    Child=>table%vec(i)%Child
    DO 
      IF (ASSOCIATED(Child)) THEN
        PRINT*, i,Child%key,Child%Val
        Child=>Child%Child
        sum=sum+1
      ELSE  
        EXIT
      END IF  
    END DO
  END DO
  PRINT*, 'Total used elements:', sum
END SUBROUTINE PrintHashTable

END MODULE hashtbl
