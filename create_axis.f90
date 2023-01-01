MODULE domain_tools

  USE ISO_FORTRAN_ENV
  IMPLICIT NONE
  SAVE

  CONTAINS

  !This subroutine allocates an axis given a start and an
  !end and a number of elements. It puts the start and the
  !end of the domain at the edges of the cells. It optionally
  !takes an "nghosts" parameter that specifies ghost cells
  !at each end of the array. You can optionally specify the
  !"delta" parameter which returns the spacing between cell
  !centres that is consistent with the requested domain length
  !and number of cells
  SUBROUTINE create_axis(axis, nels, axis_range, nghosts, delta)

    REAL(REAL64), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: axis
    INTEGER, INTENT(IN) :: nels
    REAL(REAL64), DIMENSION(2), INTENT(IN) :: axis_range
    INTEGER, INTENT(IN), OPTIONAL :: nghosts
    REAL(REAL64), INTENT(OUT), OPTIONAL :: delta
    REAL(REAL64) :: adelta
    INTEGER :: astart, aend, i, nghosts_i

    nghosts_i = 0
    IF (PRESENT(nghosts)) nghosts_i = nghosts

    astart = 1 - nghosts_i
    aend = nels + nghosts_i
    adelta = (axis_range(2) - axis_range(1)) / REAL(nels, REAL64)

    ALLOCATE(axis(astart:aend))
    axis(astart) = axis_range(1) + adelta/2.0_REAL64 &
        - adelta * REAL(nghosts_i, REAL64)

    DO i = astart+1, aend
      axis(i) = axis(i-1) + adelta
    END DO

    IF (PRESENT(delta)) delta = adelta

  END SUBROUTINE create_axis


END MODULE domain_tools
