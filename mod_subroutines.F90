MODULE MOD_SUBROUTINES
   USE MOD_PRECISION
   IMPLICIT NONE

   CONTAINS

SUBROUTINE WRITE_CSV_HEADER(filename, headers, num_headers)
   CHARACTER(LEN=*), INTENT(IN) :: filename
   CHARACTER(LEN=*), ALLOCATABLE, DIMENSION(:), INTENT(IN) :: headers
   INTEGER, INTENT(IN) :: num_headers
   INTEGER :: i
   INTEGER, PARAMETER :: unit_number = 10

   OPEN(unit_number, file=filename, status='replace', action='write')

   DO i = 1, num_headers
      IF (i /= num_headers) THEN
         WRITE(unit_number, '(a,1x)', advance='no') TRIM(headers(i))
      ELSE
         WRITE(unit_number, '(a)', advance='no') TRIM(headers(i))
      END IF
   END DO

   WRITE(unit_number, *)

   CLOSE(unit_number)
END SUBROUTINE WRITE_CSV_HEADER

SUBROUTINE append_to_csv(filename, data)
      CHARACTER(LEN=*), INTENT(IN) :: filename
      REAL(KIND=wp), DIMENSION(:), INTENT(IN) :: data
      INTEGER :: i
      INTEGER, PARAMETER :: unit_number = 10

      OPEN(unit_number, file=filename, status='old', action='write', position='append')

      DO i = 1, SIZE(data)
         IF (i /= SIZE(data)) THEN
            WRITE(unit_number, '(f0.24,1x)', advance='no') data(i)
         ELSE
            WRITE(unit_number, '(f0.24)', advance='no') data(i)
         END IF
      END DO

      WRITE(unit_number, *)

      CLOSE(unit_number)
END SUBROUTINE append_to_csv

END MODULE MOD_SUBROUTINES
