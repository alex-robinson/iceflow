
module iceflow

    use nml 
    
    ! Add physics modules here 
    ! use streamice

    implicit none
    
    ! Define all parameters needed for the yelmo module
    type iceflow_param_class
        character (len=256) :: name, boundary(30), method
        integer             :: nx, ny, npts
        double precision    :: p1, p2 
    end type

    type iceflow_state_class
        ! Model variables that the define the state of the domain
        integer,          allocatable, dimension(:,:) :: mask 
        double precision, allocatable, dimension(:,:) :: x2D, y2D, lon2D, lat2D
        double precision, allocatable, dimension(:,:) :: zs, zb, zb0 

    end type

    type iceflow_class

        ! Parameters
        type(iceflow_param_class) :: par

        ! Variables
        type(iceflow_state_class) :: now

    end type

    private
    public :: iceflow_param_class, iceflow_state_class, iceflow_class
    public :: iceflow_par_load
    public :: iceflow_init, iceflow_update, iceflow_end 
contains

    subroutine iceflow_init(flow,nx,ny,nz)

        implicit none 

        type(iceflow_class) :: flow 
        integer :: nx, ny, nz 
        
        call iceflow_alloc(now,nx,ny,nz)

        return 

    end subroutine iceflow_init


    subroutine calc_velocity(par,vx,vy,vz)

        implicit none
        
        type(iceflow_param_class) :: par

        double precision, dimension(:,:,:), pointer :: vx, vy, vz

        
        ! Calculate material properties and viscosity 


        ! First calculate basal velocity, then 
        ! the horizontal velocity field, then the vertical velocity field
        ! to get full 3D velocity field [vx,vy,vz]

!         call calc_vxy_b_sia(time, z_sl)
!         call calc_vxy_sia(dzeta_c, dzeta_t)
!         call calc_vz_sia(dxi, deta, dzeta_c, dzeta_t)

        ! Next calculate the full strain-rate tensor, 
        ! the full effective strain rate and the shear fraction
!         call calc_dxyz(dxi, deta, dzeta_c, dzeta_t)

        
        return

    end subroutine calc_velocity

    subroutine calc_vxy_b()

        implicit none


    end subroutine calc_vxy_b


    subroutine iceflow_alloc(now,nx,ny)

        implicit none 

        type(iceflow_state_class) :: now 
        integer :: nx, ny  

        allocate(now%mask(nx,ny))
        allocate(now%x2D(nx,ny))
        allocate(now%y2D(nx,ny))
        allocate(now%lon2D(nx,ny))
        allocate(now%lat2D(nx,ny))
        allocate(now%zs(nx,ny))
        allocate(now%zb(nx,ny))
        allocate(now%zb0(nx,ny))
 
        return 
    end subroutine iceflow_alloc 

    subroutine iceflow_dealloc(now)

        implicit none 

        type(iceflow_state_class) :: now

        deallocate(now%mask,now%x2D,now%y2D,now%lon2D,now%lat2D, &
                   now%zs,now%zb,now%zb0)

        return 

    end subroutine iceflow_dealloc 

    subroutine iceflow_par_load(par,filename,init)

        type(iceflow_param_class) :: par
        character(len=*)    :: filename 
        logical, optional :: init 
        logical :: init_pars 

        init_pars = .FALSE.
        if (present(init)) init_pars = .TRUE. 
 
        ! Store local parameter values in output object
        call nml_read(filename,"yelmo_velocity_greve","boundary",  par%boundary,  init=init_pars)
        call nml_read(filename,"yelmo_velocity_greve","method",    par%method,    init=init_pars)
        
        return

    end subroutine iceflow_par_load

end module iceflow
