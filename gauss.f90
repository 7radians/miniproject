PROGRAM mini_project

    USE ISO_FORTRAN_ENV
    USE domain_tools
    USE command_line
    USE particle_mover_mod
    USE write_netcdf

    IMPLICIT NONE
    
    ! Precision
    INTEGER, PARAMETER :: dp = real64

    ! Variables for Gauss Seidel
    INTEGER:: i, j, k, count
    INTEGER, parameter :: nghost = 1
    REAL(dp), DIMENSION(:, :), ALLOCATABLE :: rho, phi, ex, ey
    REAL(dp):: delta_x, delta_y, delta, denom, var, error, xx, yy, etot, drms, et, d, sumd, dx2, dy2, coef, nxreal, rh
    REAL(dp), DIMENSION(2):: axis_range = [-1, 1]
    REAL(dp), DIMENSION(:), ALLOCATABLE:: x_ax, y_ax

    ! Variables for command line read in
    CHARACTER(LEN=6):: problem, str_val
    LOGICAL:: sucess1, sucess2, sucess3, exists
    INTEGER:: int_val, nx, ny

    ! Variables for particle mover
    REAL(dp), PARAMETER :: dt = 0.01 ! Timestep
    ! Store particle position, velocity, and acceleration in derived type
    TYPE(particle), DIMENSION(0:10) :: trajectory

    ! Variables for netcdf write
    TYPE(run_data_struct) :: run_data ! Derived type for run_data
    INTEGER :: ierr ! Store error code


    ! read arguments from the command line, if they exist
    CALL parse_args
   
    ! nx read in
    sucess1 = get_arg("nx", int_val, exists=exists)
    IF (sucess1) THEN
      nx = int_val
      IF (nx<=0) THEN
        PRINT*, "Argument 'nx' must be greater than zero."
        sucess1=.FALSE.
      END IF
    ELSE IF (exists) THEN
      PRINT*, "Argument 'nx' could not be parsed, requires integer input."
    ELSE
      PRINT*, "Argument 'nx' not found in input."
    END IF   

    ! ny read in    
    sucess2 = get_arg("ny", int_val, exists=exists)
    IF (sucess2) THEN
      ny = int_val
      IF (ny<=0) THEN
        PRINT*, "Argument 'ny' must be greater than zero."
        sucess2=.FALSE.
      END IF
    ELSE IF (exists) THEN
      PRINT*, "Argument 'ny' could not be parsed, requires integer input."
    ELSE
      PRINT*, "Argument 'ny' not found in input."
    END IF   
    
    ! problem read in
    sucess3 = get_arg("problem", str_val, exists=exists)
    IF (sucess3) THEN
      problem = str_val
      IF (.NOT.((problem=='null') .OR. (problem=='single') .OR. (problem=='double'))) THEN
        PRINT*, "Argument 'problem' must have value 'null', 'single', or 'double'."
        sucess3=.FALSE.
      END IF
    ELSE
      PRINT *, "Argument 'problem' not found in input."
    END IF
    
    ! Stop program if any variables fail to read or meet requirements.
    IF (.NOT.(sucess1 .AND. sucess2 .AND. sucess3)) THEN
      PRINT*, "Command line inputs incorrectly formatted."
      STOP 10 ! Generate non-zero exit code for bash scripting
    END IF
    
    
    ! GAUSS SEIDEL
    ! allocate arrays
    
    ALLOCATE(rho(0:nx+1, 0:ny+1))
    ALLOCATE(phi(0:nx+1, 0:ny+1))
    ALLOCATE(ex(1:nx, 1:ny))
    ALLOCATE(ey(1:nx, 1:ny))

    
    ! create axes with domain range (-1, 1)
         
    CALL create_axis(x_ax, nx, axis_range, nghost, delta)
    delta_x = delta

    CALL create_axis(y_ax, ny, axis_range, nghost, delta)
    delta_y = delta
    
    ! set rho
    
    IF (problem=="null") THEN
        rho = rho*0.0_dp
    END IF
    
    IF (problem=="single") THEN 
        DO i = 1, nx
            DO j = 1, ny
                rho(i,j) = EXP(-(x_ax(i)/0.1_dp)**2 - (y_ax(j)/0.1_dp)**2)
            END DO
        END DO
    END IF
    
     IF (problem=="double") THEN
        DO i = 1, nx
            DO j = 1, ny
               var =  EXP(-((x_ax(i)-0.75_dp)/0.2_dp)**2 - ((y_ax(j)-0.75_dp)/0.2_dp)**2)
               rho(i,j) = EXP(-((x_ax(i)+0.25_dp)/0.1_dp)**2 - ((y_ax(j)+0.25_dp)/0.1_dp)**2)  + var
            END DO
        END DO
    END IF
    
    print*, 'Axes and rho initialised. Running Gauss-Seidel algorithm...' ! Testing
    
    ! Gauss-Seidel algorithm
    
    error = 1.0_dp
    nxreal = REAL(nx, dp)
    coef = SQRT(REAL(1/nxreal, dp))
    
    
    dx2 = delta_x**2
    dy2 = delta_y**2
    denom = -(2.0_dp/(dx2) + 2.0_dp/(dy2))
    count = 0 
    
    DO WHILE (error >= 0.00001_dp)
         etot = 0.0_dp
         sumd = 0.0_dp
         count = count + 1
         DO i = 1, nx
               DO j = 1, ny
                  
                  phi(i, j) = (rho(i,j) - (phi(i+1,j) + phi(i-1,j))/(dx2) - (phi(i,j+1) + phi(i,j-1))/(dy2) )/denom
                  xx = (phi(i-1,j) - 2.0_dp*phi(i,j) + phi(i+1,j))/(dx2)
                  yy = (phi(i,j-1) - 2.0_dp*phi(i,j) + phi(i,j+1))/(dy2) 
                  et = ABS( xx + yy - rho(i,j) )
                  d = xx + yy
                  etot = etot + et
                  sumd = sumd + d
               END DO
         END DO
         drms = coef * SQRT(sumd)
         IF (drms /= 0) THEN
             error = etot/drms
         END IF
         !print*, 'while loop iteration' ! Testing
    END DO
    
    print*, 'Gauss Seidel completed. Calculating electric field...' ! Testing

    ! calculate the electric field
   
    DO i = 1, nx
        DO j = 1, ny
            ex(i,j) = (phi(i+1,j) - phi(i-1,j))/(2.0_dp*delta_x)
            ey(i,j) = (phi(i,j+1) - phi(i,j-1))/(2.0_dp*delta_y)
        END DO
    END DO

    print*, 'Electric field calculated. Running particle mover...' ! Testing
   
  
    ! PARTICLE MOVER
    CALL particle_mover(trajectory, problem, ex, ey, delta_x, delta_y, dt)
    
    print*, 'Particle moved. Writing to netCDF file...' ! Testing

    ! WRITE NETCDF
    ! Set up run data
    run_data%data_writer = 'gauss.f90'
    run_data%problem = problem
    run_data%nx = nx
    run_data%ny = ny

    CALL write_data('output_data', ierr, run_data, rho, phi, ex, ey, trajectory, delta_x, delta_y, dt)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, 'Write to netCDF file failed.'
      STOP 11 ! Generate non-zero exit code for bash scripting
    ELSE
      PRINT*, 'Write to netCDF file successful.'
    END IF


    ! deallocate arrays
    DEALLOCATE(rho)
    DEALLOCATE(phi)
    DEALLOCATE(ex)
    DEALLOCATE(ey)


   
END PROGRAM mini_project

