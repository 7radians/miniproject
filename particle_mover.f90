MODULE particle_mover_mod
    
    USE iso_fortran_env
    IMPLICIT NONE
    SAVE
    
    ! Declare derived type to store particle values
    TYPE :: particle
        REAL(REAL64) :: x, y, vx, vy, ax, ay 
    END TYPE
    
    CONTAINS

    SUBROUTINE particle_mover(traj, problem, E_field_x, E_field_y, dx, dy, dt)
        ! Precision
        INTEGER, PARAMETER :: dp = real64

        ! Dummy variables
        TYPE(particle), DIMENSION(0:), INTENT(OUT) :: traj ! store evolution of particle
        REAL(dp), DIMENSION(:,:), INTENT(IN) :: E_field_x, E_field_y
        CHARACTER(LEN=6), INTENT(IN) :: problem
        REAL(dp), INTENT(IN) :: dx, dy, dt

        ! Declare variables
        REAL(dp) :: m, q, mu ! physical constants
        INTEGER(INT32) :: i, iters ! Loop varible and Number of iterations
        INTEGER(INT32) :: cell_x, cell_y ! Indices of E_field array

        q = -1.0_dp
        m = 1.0_dp
        mu = 1.0_dp

        iters = SIZE(traj)-1

        ! Initialise x and v from command line read-in (problem)
        ! Validation already carried out so should be one of three values
        SELECT CASE (problem)
            CASE('null')
                traj(0)%x = 0.0_dp
                traj(0)%y = 0.0_dp
                traj(0)%vx = 0.1_dp
                traj(0)%vy = 0.1_dp 
            CASE('single')
                traj(0)%x = 0.1_dp
                traj(0)%y = 0.0_dp
                traj(0)%vx = 0.0_dp
                traj(0)%vy = 0.0_dp 
            CASE('double')
                traj(0)%x = 0.0_dp
                traj(0)%y = 0.5_dp
                traj(0)%vx = 0.0_dp
                traj(0)%vy = 0.0_dp 
        END SELECT

        ! Get grid cell of particle from x,y
        cell_x = 1 + FLOOR((traj(0)%x - 1.0_dp)/dx)
        cell_y = 1 + FLOOR((traj(0)%y - 1.0_dp)/dy)
    
        ! Initialise a
        traj(0)%ax = q*E_field_x(cell_x, cell_y)/m
        traj(0)%ay = q*E_field_y(cell_x, cell_y)/m

        ! Move particle
        DO i = 0, iters-1
            ! New position
            traj(i+1)%x = traj(i)%x + traj(i)%vx*dt + 0.5_dp*traj(i)%ax*(dt**2)
            traj(i+1)%y = traj(i)%y + traj(i)%vy*dt + 0.5_dp*traj(i)%ay*(dt**2)

            ! New grid point
            cell_x = 1 + FLOOR((traj(i+1)%x - 1.0_dp)/dx)
            cell_y = 1 + FLOOR((traj(i+1)%y - 1.0_dp)/dy)

            ! New acceleration
            traj(i+1)%ax = q*E_field_x(cell_x,cell_y)/m
            traj(i+1)%ay = q*E_field_y(cell_x,cell_y)/m

            ! New velocity
            traj(i+1)%vx = traj(i)%vx + dt*0.5_dp*(traj(i+1)%ax + traj(i)%ax)
            traj(i+1)%vy = traj(i)%vy + dt*0.5_dp*(traj(i+1)%ay + traj(i)%ay)

        END DO
        
    END SUBROUTINE
    
END MODULE particle_mover_mod