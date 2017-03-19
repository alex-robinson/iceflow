
module iceflow

    use nml
    use solver_advdiff2D

    ! == TO DO ==   ! Add physics modules here
    ! use streamice
    
    implicit none
    
    ! Define all parameters needed for the yelmo module
    type iceflow_param_class
        character (len=256) :: method
        integer :: mix_method

        integer :: nx, ny, nz
        real(4) :: dx, dy, dz
        
        ! Internal parameters
        logical :: use_sia, use_ssa 
        
    end type

    type iceflow_state_class
        ! Model variables that the define the state of the domain

        ! Boundary variables
        integer, allocatable :: mask(:,:)
        real(4), allocatable :: H_ice(:,:), z_srf(:,:), z_bed(:,:)

        ! Dynamics variables
        real(4), allocatable :: ux_sia(:,:), uy_sia(:,:)
        real(4), allocatable :: ux_ssa(:,:), uy_ssa(:,:)
        real(4), allocatable :: f_ssa(:,:)
        real(4), allocatable :: ux(:,:,:), uy(:,:,:), uz(:,:,:)

    end type

    type iceflow_class

        ! Parameters
        type(iceflow_param_class) :: par

        ! Variables
        type(iceflow_state_class) :: now

    end type
    
    private
    public :: iceflow_param_class, iceflow_state_class, iceflow_class
    public :: iceflow_init, iceflow_init_state, iceflow_update, iceflow_end

contains

    subroutine iceflow_init(flow,filename,nx,ny,nz,dx,dy,dz)
        ! Initialize the iceflow_class object
        ! (load parameters, allocate arrays, define initial values)

        implicit none

        type(iceflow_class), intent(OUT) :: flow
        character(len=*), intent(IN)     :: filename
        integer, intent(IN) :: nx, ny, nz
        real(4), intent(IN) :: dx, dy, dz

        ! Load iceflow parameters for this domain
        call iceflow_par_load(flow%par,filename)

        ! Store grid information
        flow%par%dx = dx
        flow%par%dy = dy
        flow%par%dz = dz
        
        ! Allocate iceflow arrays
        call iceflow_alloc(flow%now,nx,ny,nz)
        
         ! Make sure solver choice makes sense
        if (trim(flow%par%method) .ne. "sia" .or. &
            trim(flow%par%method) .ne. "ssa" .or. &
            trim(flow%par%method) .ne. "sia-ssa") then
            write(*,*) "iceflow_udpate:: error: solver method not recognized: "//trim(flow%par%method)
            stop
        end if
        
        ! === Set internal parameters ==========
        
        ! Set switches to control whether sia/ssa calcs are performed
        flow%par%use_sia = .FALSE.
        flow%par%use_ssa = .TRUE. 
        if (index(trim(flow%par%method), "sia") .gt. 0) flow%par%use_sia = .TRUE.
        if (index(trim(flow%par%method), "ssa") .gt. 0) flow%par%use_ssa = .TRUE.
        
        write(*,*) "iceflow_init:: iceflow object initialized."

        return

    end subroutine iceflow_init

    subroutine iceflow_init_state(flow)
        ! Initialize the state of the iceflow_class variables
        ! (Define inline, load from file, external values, etc)

        implicit none

        type(iceflow_class) :: flow

        ! Set initial values
        flow%now%mask  =     1
        flow%now%H_ice = 100.0
        flow%now%z_srf = 100.0
        flow%now%z_bed =   0.0

        flow%now%ux = 0.0
        flow%now%uy = 0.0
        flow%now%uz = 0.0


        write(*,*) "iceflow_init_state:: iceflow state initialized."

        return

    end subroutine iceflow_init_state

    subroutine iceflow_update(flow,H_ice,z_srf,z_bed,f_grnd,T_ice,dt)
        ! Update the state of the iceflow_class variables
        ! given new boundary conditions, time step, etc.

        ! Note: input variables are provided on central grid, but
        ! velocity variables are calculated on staggered grid.
        ! Return velocity on central grid?

        implicit none

        type(iceflow_class), intent(INOUT) :: flow
        real(4), intent(IN) :: H_ice(:,:), z_srf(:,:), z_bed(:,:), f_grnd(:,:)
        real(4), intent(IN) :: T_ice(:,:,:)
        real(4), intent(IN) :: dt     ! External timestep to be matched by internal adaptive steps

       

        ! Calculate material properties and viscosity


        ! First calculate basal velocity, then
        ! the horizontal velocity field, then the vertical velocity field
        ! to get full 3D velocity field [vx,vy,vz]

        if (flow%par%use_sia) then
            ! Solver using sia this timestep

            ! == TO DO ==

        end if 
        
        if (flow%par%use_ssa) then
            ! Solve for velocity only using ssa

            ! == TO DO ==

        end if 
        
        if (trim(flow%par%method) .eq. "sia-ssa") then
            ! Hybrid method, update the mixing fraction

            ! 1. Update velocity mixing fraction f_ssa
            call determine_mixing_fraction(flow%par,flow%now%f_ssa)

            ! == TO DO ==

        end if
        
        ! Calculate the hybrid velocity fields (2D and 3D) 
        
        ! == TO DO ==
        
        
        
        
        ! Calculate some diagnostics
        
        ! Calculate the basal to surface velocity ratio, f_vbvs
        where ( ux(:,:,1)**2+uy(:,:,1)**2 .gt. 0.0) 
            vel%now%f_vbvs = ((ux(:,:,nz)**2+uy(:,:,nz)**2)**0.5 &
                                            / (ux(:,:,1)**2+uy(:,:,1)**2)**0.5)
        elsewhere 
            ! No ice (or no velocity)
            vel%now%f_vbvs = 1.0 
        end where 
        
        
        return

    end subroutine iceflow_update
    
    subroutine iceflow_update_ssa(vel,time)

        implicit none 

        type(grisli_vel_class), intent(INOUT) :: vel
        real(4), intent(IN) :: time 

        ! Local variables 
        real(4) :: dt 

        ! Get current time step 
        dt = time - vel%now%time 

        if (dt .ge. dt_ssa) then 
            ! Call ssa solver, update the ssa gradient fields 

            ! Store old solution
            vel%now%ux_ssa_old = vel%now%ux_ssa 
            vel%now%uy_ssa_old = vel%now%uy_ssa 
                
            ! Call the ssa solver
            call diagnoshelf(vel%now%ux_ssa,vel%now%uy_ssa,vel%now%f_ssa_mx,vel%now%f_ssa_my)

            ! Make an extra iteration at start to equilibrate   
            ! ajr: probably not needed...                           
            if ((restart.eq.0).and.(nt.lt.10)) then                                                          
                call diagnoshelf(vel%now%ux_ssa,vel%now%uy_ssa,vel%now%f_ssa_mx,vel%now%f_ssa_my)
            end if

            ! Calculate the velocity gradient 
            if (dt .gt. 0.0) then 
                vel%now%dux_ssa = (vel%now%ux_ssa - vel%now%ux_ssa_old) / dt
                vel%now%duy_ssa = (vel%now%uy_ssa - vel%now%uy_ssa_old) / dt
            else 
                vel%now%dux_ssa = 0.0 
                vel%now%duy_ssa = 0.0
            end if 

            ! =========
            ! To do: determine dt_ssa as a function of the ssa gradient 
            ! =========


            ! Update the current ssa time
            vel%now%time = time 

        else 
            ! Apply ssa gradient 

        end if 


        return 

    end subroutine iceflow_update_ssa
    
    subroutine iceflow_update_ssa_mixing(vel,z_srf,z_bed,H_ice,H_water,f_grnd,is_grz,f_pmp,time)
        ! Calculate the fraction of ssa that should be applied in the model
        ! Updates variables: f_ssa and ssa_active
        
        implicit none 

        type(grisli_vel_class), intent(INOUT) :: vel
        real(4), intent(IN)  :: z_srf(:,:), z_bed(:,:), H_ice(:,:), H_water(:,:), f_grnd(:,:), f_pmp(:,:) 
        logical, intent(IN)  :: is_grz(:,:)
        real(4), intent(IN)  :: time 

        !if (time .lt. -105000.0) then

        ! x-direction 
        call calc_ssa_fraction(vel%par,vel%now%f_ssa_mx,f_grnd,is_grz, &
                               vel%now%ux_bar,z_srf,H_water,f_pmp,time)

        ! y-direction
        call calc_ssa_fraction(vel%par,vel%now%f_ssa_my,f_grnd,is_grz, &
                               vel%now%uy_bar,z_srf,H_water,f_pmp,time)
        
        ! central mask
        call calc_ssa_fraction_centered(vel%par,vel%now%f_ssa,vel%now%f_ssa_mx,vel%now%f_ssa_my)

!         vel%now%f_ssa_mx = vel%now%f_ssa 
!         vel%now%f_ssa_my = vel%now%f_ssa 
        
        ! Determine potentially ssa active regions for computational efficiency
        vel%now%ssa_active = define_ssa_active(vel%par,z_srf,z_bed,H_ice,H_water,f_grnd,vel%now%f_ssa)

        !else 
    
        !   vel%now%f_ssa    =  vel%now%f_ssa
        !   vel%now%f_ssa_mx =  vel%now%f_ssa_mx
        !   vel%now%f_ssa_my =  vel%now%f_ssa_my

        !endif

        ! If desired, restrict ssa regions where necessary 
        if (vel%par%restrict_ssa) then 
            where ( .not. vel%now%ssa_active ) 
                vel%now%f_ssa    = 0.0 
                vel%now%f_ssa_mx = 0.0
                vel%now%f_ssa_my = 0.0 

            end where 
        end if 

        return 

    end subroutine iceflow_update_ssa_mixing


    subroutine iceflow_end(flow)
        ! Terminate the iceflow_class object
        ! (Finalize some calculations, deallocate arrays, etc.)

        implicit none

        type(iceflow_class) :: flow

        ! Deallocate iceflow arrays
        call iceflow_dealloc(flow%now)

        return

    end subroutine iceflow_end

    subroutine iceflow_par_load(par,filename,init)
        ! Load parameters from namelist file

        type(iceflow_param_class) :: par
        character(len=*)    :: filename
        logical, optional :: init
        logical :: init_pars

        init_pars = .FALSE.
        if (present(init)) init_pars = .TRUE.

        ! Load parameter values from file using nml library
        call nml_read(filename,"iceflow_pars","method",par%method,init=init_pars)



        write(*,*) "iceflow_par_load:: parameters loaded from: "//trim(filename)

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

        allocate(now%ux_sia(nx,ny))
        allocate(now%uy_sia(nx,ny))
        allocate(now%ux_ssa(nx,ny))
        allocate(now%uy_ssa(nx,ny))
        allocate(now%f_ssa(nx,ny))

        ! 3D arrays
        allocate(now%ux(nx,ny,nz))
        allocate(now%uy(nx,ny,nz))
        allocate(now%uz(nx,ny,nz))


        write(*,*) "iceflow_alloc:: arrays allocated (nx,ny,nz): ", nx, ny, nz

        return
    end subroutine iceflow_alloc

    subroutine iceflow_dealloc(now)

        implicit none

        type(iceflow_state_class) :: now

        ! Deallocate all allocatable arrays

        deallocate(now%mask,now%H_ice,now%z_srf,now%z_bed)
        deallocate(now%ux_sia,now%uy_sia)
        deallocate(now%ux_ssa,now%uy_ssa)
        deallocate(now%f_ssa)
        deallocate(now%ux,now%uy,now%uz)

        return

    end subroutine iceflow_dealloc

end module iceflow
