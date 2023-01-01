PROGRAM mini_project

    USE ISO_FORTRAN_ENV
    USE domain_tools
    USE command_line
    IMPLICIT NONE
    
    INTEGER:: i, j, int_val, nx, ny, k, count
    INTEGER, parameter :: nghost = 1
    REAL(REAL64), DIMENSION(:, :), ALLOCATABLE :: rho, phi
    REAL, parameter :: dt = 0.01
    REAL(REAL64):: delta_x, delta_y, delta, denom, var, error, xx, yy, etot, drms, et, d, sumd, dx2, dy2, coef, nxreal, rh
    REAL(REAL64), DIMENSION(2):: axis_range = [-1, 1]
    REAL(REAL64), DIMENSION(:), ALLOCATABLE:: x_ax, y_ax
    CHARACTER(LEN=6):: problem, str_val
    LOGICAL:: sucess1, sucess2, sucess3
   
    !read arguments from the command line, if they exist
    CALL parse_args
   
    sucess1 = get_arg("nx", int_val)
    IF (sucess1) THEN
      nx = int_val
    ELSE
      PRINT *, "Failed to get nx"
    END IF   
        
    sucess2 = get_arg("ny", int_val)
    IF (sucess2) THEN
      ny = int_val
    ELSE
      PRINT *, "Failed to get ny"
    END IF
    
    sucess3 = get_arg("problem", str_val)
    IF (sucess3) THEN
      problem = str_val
    ELSE
      PRINT *, "Failed to get problem"
    END IF
    
  
    ! allocate arrays
    
    ALLOCATE(rho(0:nx+1, 0:ny+1))
    ALLOCATE(phi(0:nx+1, 0:ny+1))

    
    ! create axes with domain range (-1, 1)
         
    CALL create_axis(x_ax, nx, axis_range, nghost, delta)
    delta_x = delta
    !print *, delta_x
    
    CALL create_axis(y_ax, ny, axis_range, nghost, delta)
    delta_y = delta
    !print *, delta_y
    
    ! set rho
    
    IF (problem=="null") THEN
        rho = rho*0
    END IF
    
    IF (problem=="single") THEN 
        DO i = 1, nx
            DO j = 1, ny
                rho(i,j) = EXP(-(x_ax(i)/0.1)**2 - (y_ax(j)/0.1)**2)
            END DO
        END DO
    END IF
    !print *, rho

    
     IF (problem=="double") THEN
        DO i = 1, nx
            DO j = 1, ny
               var =  EXP(-((x_ax(i)-0.75)/0.2)**2 - ((y_ax(j)-0.75)/0.2)**2)
               rho(i,j) = EXP(-((x_ax(i)+0.25)/0.1)**2 - ((y_ax(j)+0.25)/0.1)**2)  + var
            END DO
        END DO
    END IF
    
    
    ! Gauss-Seidel algorithm
    
    error = 1
    nxreal = REAL(nx, REAL64)
    coef = SQRT(REAL(1/nxreal, REAL64))
    
    
    dx2 = delta_x**2
    dy2 = delta_y**2
    denom = -(2/(dx2) + 2/(dy2))
    !print *, denom
    count = 0 
    
    DO WHILE (error >= 0.00001)
         etot = 0
         sumd = 0
         count = count + 1
         DO i = 1, nx
               DO j = 1, ny
                  
                  phi(i, j) = (rho(i,j) - (phi(i+1,j) + phi(i-1,j))/(dx2) - (phi(i,j+1) + phi(i,j-1))/(dy2) )/denom
                  xx = (phi(i-1,j) - 2*phi(i,j) + phi(i+1,j))/(dx2)
                  yy = (phi(i,j-1) - 2*phi(i,j) + phi(i,j+1))/(dy2) 
                  et = ABS( xx + yy - rho(i,j) )
                  !et(i, j) = ABS( (phi(i-1,j) - 2*phi(i,j) + phi(i+1,j))/(delta_x**2) + (phi(i,j-1) - 2*phi(i,j) + phi(i,j+1))/(delta_y**2) - rho(i,j) )
                  !d(i, j) = SUM( (phi(i-1,j) - 2*phi(i,j) + phi(i+1,j))/(delta_x**2) + (phi(i,j-1) - 2*phi(i,j) + phi(i,j+1))/(delta_y**2))
                  d = xx + yy
                  etot = etot + et
                  sumd = sumd + d
               END DO
         END DO
         drms = coef * SQRT(sumd)
         IF (drms /= 0) THEN
             error = etot/drms
             print *, "error", error
         END IF
    END DO
   ! print *, "number of iterations", count
END PROGRAM mini_project

