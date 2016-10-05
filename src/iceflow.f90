
module yelmo_velocity

    use nml 
    
    implicit none
    
    ! Define all parameters needed for the yelmo module
    type yvel_param_class
        character (len=256) :: name, boundary(30), method
        integer             :: nx, ny, npts
        double precision    :: p1, p2 
    end type

    type yvel_state_class
        ! Model variables that the define the state of the domain
        integer,          allocatable, dimension(:,:) :: mask 
        double precision, allocatable, dimension(:,:) :: x2D, y2D, lon2D, lat2D
        double precision, allocatable, dimension(:,:) :: zs, zb, zb0 

    end type

    type boundary_opt_class 
        logical :: tsurf, smb 
    end type

    type yvel_class

        type(yvel_param_class)   :: par        ! physical parameters
        type(boundary_opt_class)  :: bnd, bnd0  ! boundary switches (bnd0 for equilibration)

        ! All variables
        type(yvel_state_class) :: now

    end type

    private
    public :: yvel_class, yvel_state_class, yvel_param_class
    public :: yvel_par_load, yvel_alloc, yvel_dealloc
    public :: yvel_boundary_define

contains

    subroutine calc_velocity(par,vx,vy,vz)

        implicit none
        
        type(yvel_param_class) :: par

        double precision, dimension(:,:,:), pointer :: vx, vy, vz

        

        ! First calculate basal velocity



        ! Next calculate the horizontal velocity field



        ! Now calculate the vertical velocity field

        
        return

    end subroutine calc_velocity

    subroutine calc_vxy_b()

        implicit none


    end subroutine calc_vxy_b


    subroutine yvel_alloc(now,nx,ny)

        implicit none 

        type(yvel_state_class) :: now 
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
    end subroutine yvel_alloc 

    subroutine yvel_dealloc(now)

        implicit none 

        type(yvel_state_class) :: now

        deallocate(now%mask,now%x2D,now%y2D,now%lon2D,now%lat2D, &
                   now%zs,now%zb,now%zb0)

        return 

    end subroutine yvel_dealloc 

    subroutine yvel_par_load(par,filename,init)

        type(yvel_param_class) :: par
        character(len=*)    :: filename 
        logical, optional :: init 
        logical :: init_pars 

        init_pars = .FALSE.
        if (present(init)) init_pars = .TRUE. 
 
        ! Store local parameter values in output object
        call nml_read(filename,"yelmo_velocity_greve","boundary",  par%boundary,  init=init_pars)
        call nml_read(filename,"yelmo_velocity_greve","method",    par%method,    init=init_pars)
        
        return

    end subroutine yvel_par_load

    subroutine yvel_boundary_define(bnd,boundary)

        implicit none 

        type(boundary_opt_class) :: bnd 
        character(len=256) :: boundary(:)
        integer :: q 

        ! First set all boundary fields to false
        bnd%tsurf   = .FALSE. 
        bnd%smb     = .FALSE. 

        ! Now find boundary fields 
        do q = 1,size(boundary)

            select case(trim(boundary(q)))

                case("tsurf")
                    bnd%tsurf   = .TRUE. 
                case("smb")
                    bnd%smb     = .TRUE. 
                case DEFAULT 
                    ! pass 
            end select 
        end do 

        return 

    end subroutine yvel_boundary_define


end module yelmo_velocity
