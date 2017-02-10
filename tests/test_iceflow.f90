program test_iceflow
    ! Program to test development of iceflow library

    use iceflow
    use ncio

    implicit none

    type(iceflow_class) :: flow1
    integer :: t, nt

    ! Input variables to iceflow_update
    real(4), allocatable :: H_ice(:,:), z_srf(:,:), z_bed(:,:), f_grnd(:,:)
    real(4), allocatable :: T_ice(:,:,:)

    ! Initialize our local flow1 object
    call iceflow_init(flow1,"Greenland.nml",nx=40,ny=40,nz=50,dx=40.0,dy=40.0,dz=1.0)

    ! Initialize the state of the variables
    call iceflow_init_state(flow1)

    ! == TO DO ==
    ! Load some precalculated state variables, internal temperature,
    ! ice thickness, etc, to test flow calculations

    ! Calculate several time steps
    nt = 10

    do t = 1, nt

        ! Update flow1
        call iceflow_update(flow1,dt=1.0,dx=dx,dy=dy)

    end do

    write(*,*)
    write(*,*) "test_iceflow complete."
    write(*,*)

end program test_iceflow
