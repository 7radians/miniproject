MODULE write_netcdf

    ! Adapted from Workshop 7 netcdf example

    USE iso_fortran_env
    USE netcdf
    USE particle_mover_mod

    IMPLICIT NONE

    TYPE :: run_data_struct
        CHARACTER(LEN=30) :: data_writer ! store name of producing code
        INTEGER(INT32) :: nx, ny ! store command line arguments
        CHARACTER(LEN=6) :: problem ! store command line arguments
    END TYPE

    CONTAINS

    SUBROUTINE write_data(filename, ierr, run_data, rho, phi, E_field_x, E_field_y, trajectory, dx, dy, dt)

        ! Precision
        INTEGER, PARAMETER :: dp = real64

        ! Dummy arguments
        CHARACTER(LEN=*), INTENT(IN) :: filename
        INTEGER(INT32), INTENT(OUT) :: ierr
        TYPE(run_data_struct), INTENT(IN) :: run_data
        REAL(dp), DIMENSION(:,:), INTENT(IN) :: rho, phi, E_field_x, E_field_y
        TYPE(particle), DIMENSION(0:), INTENT(IN) :: trajectory
        REAL(dp), INTENT(IN) :: dx, dy, dt

        ! Id variables
        INTEGER :: file_id, rho_id, phi_id, Ex_id, Ey_id
        INTEGER :: traj_x_id, traj_y_id, traj_vx_id, traj_vy_id, traj_ax_id, traj_ay_id
        INTEGER :: i

        ! E field variables
        INTEGER, PARAMETER :: ndims_grid = 2
        INTEGER, DIMENSION(ndims_grid) :: size_grid, dims_grid_ids
        CHARACTER(LEN=6), DIMENSION(ndims_grid) :: dims_grid = (/'x_grid','y_grid'/)

        ! Trajectory variables
        INTEGER :: size_traj, dims_traj_id
        CHARACTER(LEN=1) :: dims_traj = 't'
        
        ! Axes
        REAL(dp), DIMENSION(:), ALLOCATABLE :: x_axis, y_axis, t_axis
        INTEGER :: x_id, y_id, t_id

        ! Units
        CHARACTER(LEN=15), PARAMETER :: unit = 'arbitrary units' 


        ! Get sizes of arrays 
        size_grid = SHAPE(rho)
        size_traj = SIZE(trajectory)
        
        ! Create axes
        ALLOCATE(x_axis(size_grid(1)))
        ALLOCATE(y_axis(size_grid(2)))
        ALLOCATE(t_axis(0:size_traj-1))

        ! Grid axis
        ! Values defined at the centre of grid cells, in domain -1 to 1
        ! x axis
        DO i = 1, size_grid(1)
            x_axis(i) = -1.0_dp + dx*(i-0.5_dp)
        END DO
        ! y axis
        DO i = 1, size_grid(2)
            y_axis(i) = -1.0_dp + dy*(i-0.5_dp)
        END DO
        ! t axis
        DO i = 0, size_traj-1
            t_axis(i) = dt*i
        END DO

        
        ! Create file with filename, overwrite if needed.
        ierr = nf90_create(filename, NF90_CLOBBER, file_id)
        ! If creating file fails, return to caller.
        IF (ierr /= nf90_noerr) THEN
            PRINT*, TRIM(nf90_strerror(ierr))
            RETURN
        END IF


        !!! Write metadata

        ! Run data
        ! Identifying producing code
        ierr = nf90_put_att(file_id, NF90_GLOBAL, "data_writer", run_data%data_writer)
        IF (ierr /= nf90_noerr) THEN
            PRINT*, TRIM(nf90_strerror(ierr))
            RETURN
        END IF

        ! Write command line inputs
        ierr = nf90_put_att(file_id, NF90_GLOBAL, "problem", run_data%problem)
        IF (ierr /= nf90_noerr) THEN
            PRINT*, TRIM(nf90_strerror(ierr))
            RETURN
        END IF
        
        ierr = nf90_put_att(file_id, NF90_GLOBAL, "nx", run_data%nx)
        IF (ierr /= nf90_noerr) THEN
            PRINT*, TRIM(nf90_strerror(ierr))
            RETURN
        END IF

        ierr = nf90_put_att(file_id, NF90_GLOBAL, "ny", run_data%ny)
        IF (ierr /= nf90_noerr) THEN
            PRINT*, TRIM(nf90_strerror(ierr))
            RETURN
        END IF
        

        ! Grid variables metadata
        ! Dimensions
        DO i = 1, ndims_grid
            ierr = nf90_def_dim(file_id, dims_grid(i), size_grid(i), dims_grid_ids(i))
            IF (ierr /= nf90_noerr) THEN
                PRINT*, TRIM(nf90_strerror(ierr))
                RETURN
            END IF
        END DO

        ! rho
        ierr = nf90_def_var(file_id, "rho", NF90_REAL, dims_grid_ids, rho_id)
        IF (ierr /= nf90_noerr) THEN
            PRINT*, TRIM(nf90_strerror(ierr))
            RETURN
        END IF
        ! rho units
        ierr = nf90_put_att(file_id, rho_id, 'unit', unit)
        IF (ierr /= nf90_noerr) THEN
            PRINT*, TRIM(nf90_strerror(ierr))
            RETURN
        END IF

        ! phi
        ierr = nf90_def_var(file_id, "phi", NF90_REAL, dims_grid_ids, phi_id)
        IF (ierr /= nf90_noerr) THEN
            PRINT*, TRIM(nf90_strerror(ierr))
            RETURN
        END IF
        ! phi units
        ierr = nf90_put_att(file_id, phi_id, 'unit', unit)
        IF (ierr /= nf90_noerr) THEN
            PRINT*, TRIM(nf90_strerror(ierr))
            RETURN
        END IF

        ! x electric field 
        ierr = nf90_def_var(file_id, "Ex", NF90_REAL, dims_grid_ids, Ex_id)
        IF (ierr /= nf90_noerr) THEN
            PRINT*, TRIM(nf90_strerror(ierr))
            RETURN
        END IF
        ! ex units
        ierr = nf90_put_att(file_id, Ex_id, 'unit', unit)
        IF (ierr /= nf90_noerr) THEN
            PRINT*, TRIM(nf90_strerror(ierr))
            RETURN
        END IF
        ! y electric field 
        ierr = nf90_def_var(file_id, "Ey", NF90_REAL, dims_grid_ids, Ey_id)
        IF (ierr /= nf90_noerr) THEN
            PRINT*, TRIM(nf90_strerror(ierr))
            RETURN
        END IF
        ! ey units
        ierr = nf90_put_att(file_id, Ey_id, 'unit', unit)
        IF (ierr /= nf90_noerr) THEN
            PRINT*, TRIM(nf90_strerror(ierr))
            RETURN
        END IF

        ! Axes metadata
        ! x axis of grids
        ierr = nf90_def_var(file_id, "x_axis", NF90_REAL, dims_grid_ids(1), x_id)
        IF (ierr /= nf90_noerr) THEN
            PRINT*, TRIM(nf90_strerror(ierr))
            RETURN
        END IF
        ! x axis units
        ierr = nf90_put_att(file_id, x_id, 'unit', unit)
        IF (ierr /= nf90_noerr) THEN
            PRINT*, TRIM(nf90_strerror(ierr))
            RETURN
        END IF
        ! y axis of grids
        ierr = nf90_def_var(file_id, "y_axis", NF90_REAL, dims_grid_ids(2), y_id)
        IF (ierr /= nf90_noerr) THEN
            PRINT*, TRIM(nf90_strerror(ierr))
            RETURN
        END IF
        ! y axis units
        ierr = nf90_put_att(file_id, y_id, 'unit', unit)
        IF (ierr /= nf90_noerr) THEN
            PRINT*, TRIM(nf90_strerror(ierr))
            RETURN
        END IF


        ! Trajectory metadata
        ! Dimension
        ierr = nf90_def_dim(file_id, dims_traj, size_traj, dims_traj_id)
        IF (ierr /= nf90_noerr) THEN
            PRINT*, TRIM(nf90_strerror(ierr))
            RETURN
        END IF

        ! x position
        ierr = nf90_def_var(file_id, "x", NF90_REAL, dims_traj_id, traj_x_id)
        IF (ierr /= nf90_noerr) THEN
            PRINT*, TRIM(nf90_strerror(ierr))
            RETURN
        END IF
        ! x units
        ierr = nf90_put_att(file_id, traj_x_id, 'unit', unit)
        IF (ierr /= nf90_noerr) THEN
            PRINT*, TRIM(nf90_strerror(ierr))
            RETURN
        END IF

        ! y position
        ierr = nf90_def_var(file_id, "y", NF90_REAL, dims_traj_id, traj_y_id)
        IF (ierr /= nf90_noerr) THEN
            PRINT*, TRIM(nf90_strerror(ierr))
            RETURN
        END IF
        ! y units
        ierr = nf90_put_att(file_id, traj_y_id, 'unit', unit)
        IF (ierr /= nf90_noerr) THEN
            PRINT*, TRIM(nf90_strerror(ierr))
            RETURN
        END IF

        ! x velocity
        ierr = nf90_def_var(file_id, "vx", NF90_REAL, dims_traj_id, traj_vx_id)
        IF (ierr /= nf90_noerr) THEN
            PRINT*, TRIM(nf90_strerror(ierr))
            RETURN
        END IF     
        ! vx units
        ierr = nf90_put_att(file_id, traj_vx_id, 'unit', unit)
        IF (ierr /= nf90_noerr) THEN
            PRINT*, TRIM(nf90_strerror(ierr))
            RETURN
        END IF

        ! y velocity
        ierr = nf90_def_var(file_id, "vy", NF90_REAL, dims_traj_id, traj_vy_id)
        IF (ierr /= nf90_noerr) THEN
            PRINT*, TRIM(nf90_strerror(ierr))
            RETURN
        END IF
        ! vy units
        ierr = nf90_put_att(file_id, traj_vy_id, 'unit', unit)
        IF (ierr /= nf90_noerr) THEN
            PRINT*, TRIM(nf90_strerror(ierr))
            RETURN
        END IF

        ! x acceleration
        ierr = nf90_def_var(file_id, "ax", NF90_REAL, dims_traj_id, traj_ax_id)
        IF (ierr /= nf90_noerr) THEN
            PRINT*, TRIM(nf90_strerror(ierr))
            RETURN
        END IF
        ! ax units
        ierr = nf90_put_att(file_id, traj_ax_id, 'unit', unit)
        IF (ierr /= nf90_noerr) THEN
            PRINT*, TRIM(nf90_strerror(ierr))
            RETURN
        END IF

        ! y acceleration
        ierr = nf90_def_var(file_id, "ay", NF90_REAL, dims_traj_id, traj_ay_id)
        IF (ierr /= nf90_noerr) THEN
            PRINT*, TRIM(nf90_strerror(ierr))
            RETURN
        END IF
        ! ay units
        ierr = nf90_put_att(file_id, traj_ay_id, 'unit', unit)
        IF (ierr /= nf90_noerr) THEN
            PRINT*, TRIM(nf90_strerror(ierr))
            RETURN
        END IF

        ! Axis metadata
        ! t axis of trajectory
        ierr = nf90_def_var(file_id, "t_axis", NF90_REAL, dims_traj_id, t_id)
        IF (ierr /= nf90_noerr) THEN
            PRINT*, TRIM(nf90_strerror(ierr))
            RETURN
        END IF
        ! t units
        ierr = nf90_put_att(file_id, t_id, 'unit', unit)
        IF (ierr /= nf90_noerr) THEN
            PRINT*, TRIM(nf90_strerror(ierr))
            RETURN
        END IF


        !!! End metadata
        ierr = nf90_enddef(file_id)
        IF (ierr /= nf90_noerr) THEN
            PRINT*, TRIM(nf90_strerror(ierr))
            RETURN
        END IF


        !!! Write data
        ! Grid data
        ! rho
        ierr = nf90_put_var(file_id, rho_id, rho)
        IF (ierr /= nf90_noerr) THEN
            PRINT*, TRIM(nf90_strerror(ierr))
            RETURN
        END IF

        ! phi
        ierr = nf90_put_var(file_id, phi_id, phi)
        IF (ierr /= nf90_noerr) THEN
            PRINT*, TRIM(nf90_strerror(ierr))
            RETURN
        END IF

        ! x electric field 
        ierr = nf90_put_var(file_id, Ex_id, E_field_x)
        IF (ierr /= nf90_noerr) THEN
            PRINT*, TRIM(nf90_strerror(ierr))
            RETURN
        END IF

        ! y electric field 
        ierr = nf90_put_var(file_id, Ey_id, E_field_y)
        IF (ierr /= nf90_noerr) THEN
            PRINT*, TRIM(nf90_strerror(ierr))
            RETURN
        END IF

        ! x axis
        ierr = nf90_put_var(file_id, x_id, x_axis)
        IF (ierr /= nf90_noerr) THEN
            PRINT*, TRIM(nf90_strerror(ierr))
            RETURN
        END IF

        ! y axis
        ierr = nf90_put_var(file_id, y_id, y_axis)
        IF (ierr /= nf90_noerr) THEN
            PRINT*, TRIM(nf90_strerror(ierr))
            RETURN
        END IF


        ! Trajectory data
        ! x position
        ierr = nf90_put_var(file_id, traj_x_id, trajectory%x)
        IF (ierr /= nf90_noerr) THEN
            PRINT*, TRIM(nf90_strerror(ierr))
            RETURN
        END IF

        ! y position
        ierr = nf90_put_var(file_id, traj_y_id, trajectory%y)
        IF (ierr /= nf90_noerr) THEN
            PRINT*, TRIM(nf90_strerror(ierr))
            RETURN
        END IF

        ! x velocity
        ierr = nf90_put_var(file_id, traj_vx_id, trajectory%vx)
        IF (ierr /= nf90_noerr) THEN
            PRINT*, TRIM(nf90_strerror(ierr))
            RETURN
        END IF

        ! y velocity
        ierr = nf90_put_var(file_id, traj_vy_id, trajectory%vy)
        IF (ierr /= nf90_noerr) THEN
            PRINT*, TRIM(nf90_strerror(ierr))
            RETURN
        END IF

        ! x acceleration
        ierr = nf90_put_var(file_id, traj_ax_id, trajectory%ax)
        IF (ierr /= nf90_noerr) THEN
            PRINT*, TRIM(nf90_strerror(ierr))
            RETURN
        END IF

        ! y acceleration
        ierr = nf90_put_var(file_id, traj_ay_id, trajectory%ay)
        IF (ierr /= nf90_noerr) THEN
            PRINT*, TRIM(nf90_strerror(ierr))
            RETURN
        END IF

        ! t axis
        ierr = nf90_put_var(file_id, t_id, t_axis)
        IF (ierr /= nf90_noerr) THEN
            PRINT*, TRIM(nf90_strerror(ierr))
            RETURN
        END IF
    

        ! Close the file
        ierr = nf90_close(file_id)
        IF (ierr /= nf90_noerr) THEN
            PRINT*, TRIM(nf90_strerror(ierr))
            RETURN
        END IF

        DEALLOCATE(x_axis)
        DEALLOCATE(y_axis)
        DEALLOCATE(t_axis)

    END SUBROUTINE


END MODULE write_netcdf