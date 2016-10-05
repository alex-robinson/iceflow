
module iceflow

    use nml 
    
    ! Add physics modules here 
    ! use streamice

    implicit none
    
    ! Define all parameters needed for the yelmo module
    type iceflow_param_class
        character (len=256) :: method

    end type

    type iceflow_state_class
        ! Model variables that the define the state of the domain

        ! 2D variables
        integer, allocatable :: mask(:,:) 
        real(4), allocatable :: H_ice(:,:), z_srf(:,:), z_bed(:,:)

        ! 3D variables
        real(4), allocatable :: vx(:,:,:), vy(:,:,:), vz(:,:,:) 

    end type

    type iceflow_class

        ! Parameters
        type(iceflow_param_class) :: par

        ! Variables
        type(iceflow_state_class) :: now

    end type

    private
    public :: iceflow_param_class, iceflow_state_class, iceflow_class
    public :: iceflow_init, iceflow_update, iceflow_end 

contains

    subroutine iceflow_init(flow,filename,nx,ny,nz)
        ! Initialize the iceflow_class object
        ! (load parameters, allocate arrays, define initial values)

        implicit none 

        type(iceflow_class) :: flow 
        character(len=*)    :: filename 
        integer :: nx, ny, nz 

        ! Load iceflow parameters for this domain 
        call iceflow_par_load(flow%par,filename)

        ! Allocate iceflow arrays 
        call iceflow_alloc(flow%now,nx,ny,nz)

        return 

    end subroutine iceflow_init

    subroutine iceflow_init_state(flow,filename)
        ! Initialize the state of the iceflow_class variables
        ! (Define inline, load from file, external values, etc)

        implicit none 

        type(iceflow_class) :: flow 
        character(len=*)    :: filename 

        ! Set initial values 
        flow%now%mask = 0.0 

        return 

    end subroutine iceflow_init_state

    subroutine iceflow_update(flow,dt)
        ! Update the state of the iceflow_class variables
        ! given new boundary conditions, time step, etc.

        implicit none 

        type(iceflow_class) :: flow 
        real(4) :: dt  

        ! Calculate material properties and viscosity 

        ! First calculate basal velocity, then 
        ! the horizontal velocity field, then the vertical velocity field
        ! to get full 3D velocity field [vx,vy,vz]

        if (trim(flow%par%method) .eq. "sia") then 
            ! Solver using sia 

            ! == TO DO == 


        else if (trim(flow%par%method) .eq. "ssa") then 
            ! Solve for velocity only using ssa 

            ! == TO DO == 

            
        else if (trim(flow%par%method) .eq. "sia-ssa") then 
            ! Solve using hybrid method 

            ! == TO DO == 

            
        else 

            write(*,*) "iceflow_udpate:: error: solver method not recognized: "//trim(flow%par%method)
            stop 

        end if 

        return 

    end subroutine iceflow_update

    subroutine iceflow_par_load(par,filename,init)
        ! Load parameters from namelist file 

        type(iceflow_param_class) :: par
        character(len=*)    :: filename 
        logical, optional :: init 
        logical :: init_pars 

        init_pars = .FALSE.
        if (present(init)) init_pars = .TRUE. 
 
        ! Load parameter values from file using nml library 
        call nml_read(filename,"iceflow","method",par%method,init=init_pars)
        
        return

    end subroutine iceflow_par_load

    subroutine iceflow_alloc(now,nx,ny,nz)

        implicit none 

        type(iceflow_state_class) :: now 
        integer :: nx, ny, nz  

        ! Allocate all arrays to proper size

        ! 2D arrays
        allocate(now%mask(nx,ny))
        allocate(now%H_ice(nx,ny))
        allocate(now%z_srf(nx,ny))
        allocate(now%z_bed(nx,ny))

        ! 3D arrays
        allocate(now%vx(nx,ny,nz))
        allocate(now%vy(nx,ny,nz))
        allocate(now%vz(nx,ny,nz))
 
        return 
    end subroutine iceflow_alloc 

    subroutine iceflow_dealloc(now)

        implicit none 

        type(iceflow_state_class) :: now

        ! Deallocate all allocatable arrays

        deallocate(now%mask,now%H_ice,now%z_srf,now%z_bed)
        deallocate(now%vx,now%vy,now%vz)

        return 

    end subroutine iceflow_dealloc 

end module iceflow
