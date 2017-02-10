module solver_advdiff2D

    implicit none

contains

    subroutine  advdiff2D(H,H_old,B_old,advx,advy,Dffx,Dffy,mdot,dt,dx)
        ! To solve the 2D adevection-diffusion equation:
        ! dH/dt =
        ! M H = Frelax

        implicit none

        real(4), intent(OUT) :: H(:,:)          ! new H returned from solver
        real(4), intent(IN)  :: H_old(:,:)      ! H from previous time step
        real(4), intent(IN)  :: B_old(:,:)      ! H from previous time step
        real(4), intent(IN)  :: Advx(:,:)       ! Advective velocity - x direction
        real(4), intent(IN)  :: Advy(:,:)       ! Advective velocity - y direction
        real(4), intent(IN)  :: Dffx(:,:)       ! Diffusive velocity - x direction
        real(4), intent(IN)  :: Dffy(:,:)       ! Diffusive velocity - y direction
        real(4), intent(IN)  :: mdot(:,:)       ! Total column mass balance
        real(4), intent(IN)  :: dt, dx          ! Timestep and resolution (assumes dx=dy)

        ! Local variables
        integer :: i, j, nx, ny
        real(4) :: dx1, dtdx, dtdx2

        real(4) :: frdx,frdy                    ! helper diffusion terms
        real(4) :: fraxw,fraxe,frays,frayn      ! helper advection terms
        real(4) :: reste, delh, testh
        logical :: stopp
        integer :: ntour

        real(4), allocatable :: crelax(:,:)      ! diagnonale de M
        real(4), allocatable :: arelax(:,:)      ! sous diagonale selon x
        real(4), allocatable :: brelax(:,:)      ! sur diagonale selon x
        real(4), allocatable :: drelax(:,:)      ! sous diagonale selon y
        real(4), allocatable :: erelax(:,:)      ! sur diagonale selon y
        real(4), allocatable :: frelax(:,:)      ! vecteur
        real(4), allocatable :: c_west(:,:)      ! sur demi mailles Ux
        real(4), allocatable :: c_east(:,:)      ! sur demi mailles Ux
        real(4), allocatable :: c_north(:,:)     ! sur demi mailles Uy
        real(4), allocatable :: c_south(:,:)     !sur demi mailles Ux

        real(4), allocatable :: bdx(:,:)         ! Bedrock slope - x direction
        real(4), allocatable :: bdy(:,:)         ! Bedrock slope - y direction
        real(4), allocatable :: hdx(:,:)         ! Thickness gradient - x direction
        real(4), allocatable :: hdy(:,:)         ! Thickness gradient - y direction

        real(4), allocatable :: deltaH(:,:)     ! Change in H

        real(4), parameter :: omega  = 2.5   ! parametre schema temporel de la resolution partie diffusion
        real(4), parameter :: mu_adv = 1.0   ! parametre schema temporel de la resolution advection
        real(4), parameter :: upwind = 1.0   ! schema spatial pour l'advection

        ! Determine array size
        nx = size(H,1)
        ny = size(H,2)

        ! Allocate local arrays
        allocate(crelax(nx,ny))
        allocate(arelax(nx,ny))
        allocate(brelax(nx,ny))
        allocate(drelax(nx,ny))
        allocate(erelax(nx,ny))
        allocate(frelax(nx,ny))
        allocate(c_west(nx,ny))
        allocate(c_east(nx,ny))
        allocate(c_north(nx,ny))
        allocate(c_south(nx,ny))

        allocate(bdx(nx,ny))
        allocate(bdy(nx,ny))
        allocate(hdx(nx,ny))
        allocate(hdy(nx,ny))

        allocate(deltaH(nx,ny))

        ! Define some helpful values
        dx1   = 1.0/dx
        dtdx2 = dt/(dx**2)
        dtdx  = dt/dx

        ! Caculate the gradients of H and bedrock
        hdx = dx1*(H_old-eoshift(H_old,shift=-1,boundary=0.0,dim=1))
        hdy = dx1*(H_old-eoshift(H_old,shift=-1,boundary=0.0,dim=2))
        bdx = dx1*(B_old-eoshift(B_old,shift=-1,boundary=0.0,dim=1))
        bdy = dx1*(B_old-eoshift(B_old,shift=-1,boundary=0.0,dim=2))

        ! Initialize relaxation arrays
        Arelax = 0.0
        Brelax = 0.0
        Drelax = 0.0
        Erelax = 0.0
        Crelax = 1.0
        Frelax = 0.0
        DeltaH = 0.0

        ! Set boundary value to zero
        H(1,:)  = 0.0  ! left  border
        H(nx,:) = 0.0  ! right border
        H(:,1)  = 0.0  ! lower border
        H(:,ny) = 0.0  ! upper border

        ! Modify coefficients depending on method (upwind, central)
        if (upwind .eq. 1.0) then
            ! Upwind method

            where (Advx.ge.0.0)
                c_west = 1.0
                c_east = 0.0
            elsewhere
                c_west = 0.0
                c_east = 1.0
            end where

            where (Advy.ge.0.0)
                c_south = 1.0
                c_north = 0.0
            elsewhere
                c_south = 0.0
                c_north = 1.0
            end where

        else if (upwind .lt. 1.0) then
            ! Central method

            c_west  = 0.5
            c_east  = 0.5
            c_south = 0.5
            c_north = 0.5

        end if

        ! attribution des elements des diagonales
        do j=2,ny-1
        do i=2,nx-1

            !  sous diagonale en x
            arelax(i,j) = -omega*dtdx2*Dffx(i,j)   &               ! partie diffusive en x
                    - mu_adv*dtdx*c_west(i,j)*Advx(i,j)            ! partie advective en x

            !  sur diagonale en x
            brelax(i,j) = -omega*dtdx2*Dffx(i+1,j)  &              ! partie diffusive
                    + mu_adv*dtdx*c_east(i+1,j)*Advx(i+1,j)        ! partie advective

            !  sous diagonale en y
            drelax(i,j) = -omega*dtdx2*Dffy(i,j)   &               ! partie diffusive en y
                    - mu_adv*dtdx*c_south(i,j)*Advy(i,j)           ! partie advective en y

            !  sur diagonale en y
            erelax(i,j) = -omega*dtdx2*Dffy(i,j+1)  &              ! partie diffusive
                    + mu_adv*dtdx*c_north(i,j+1)*Advy(i,j+1)       ! partie advective



            ! diagonale
            crelax(i,j) = omega*dtdx2*((Dffx(i+1,j)+Dffx(i,j))  &
                                     + (Dffy(i,j+1)+Dffy(i,j)))

            crelax(i,j) = crelax(i,j) + mu_adv*dtdx* &
                      ( (c_west(i+1,j)*Advx(i+1,j) - c_east(i,j)*Advx(i,j)) &
                      +(c_south(i,j+1)*Advy(i,j+1) - c_north(i,j)*Advy(i,j)))

            crelax(i,j) = 1.0 + crelax(i,j)


            ! terme du vecteur

            frdx = -Dffx(i,j) * (Bdx(i,j)  +(1.-omega)*Hdx(i,j))         &  ! partie diffusive en x
                 +  Dffx(i+1,j)*(Bdx(i+1,j)+(1.-omega)*Hdx(i+1,j))

            frdy = -Dffy(i,j) * (Bdy(i,j)  +(1.-omega)*Hdy(i,j))         &  ! partie diffusive en y
                 +  Dffy(i,j+1)*(Bdy(i,j+1)+(1.-omega)*Hdy(i,j+1))

            fraxw = -c_west(i,j)*  Advx(i,j) * H_old(i-1,j)           &  ! partie advective en x
                  +  c_west(i+1,j)*Advx(i+1,j)*H_old(i,j)                ! venant de l'west

            fraxe = -c_east(i,j) * Advx(i,j) * H_old(i,j)             &  ! partie advective en x
                  +  c_east(i+1,j)*Advx(i+1,j)*H_old(i+1,j)              ! venant de l'est

            frays = -c_south(i,j) * Advy(i,j) * H_old(i,j-1)          &  ! partie advective en y
                  +  c_south(i,j+1)*Advy(i,j+1)*H_old(i,j)               ! venant du sud

            frayn = -c_north(i,j) * Advy(i,j) * H_old(i,j)            &  ! partie advective en y
                  +  c_north(i,j+1)*Advy(i,j+1)*H_old(i,j+1)             ! venant du nord

            ! Combine all terms
            frelax(i,j) = H_old(i,j) + dt*mdot(i,j) + dtdx*(frdx+frdy) &
                     + (1.0-mu_adv)*dtdx*((fraxw+fraxe)+(frays+frayn))

        end do
        end do

        ! Relaxation loop
        stopp = .false.
        ntour=0

        do while(.not. stopp)

            ntour = ntour + 1

            ! calculate change in H
            do j=2,ny-1
            do i=2,nx-1

                reste = (((arelax(i,j)*H(i-1,j) + drelax(i,j)*H(i,j-1)) &
                        + (brelax(i,j)*H(i+1,j) + erelax(i,j)*H(i,j+1))) &
                        + crelax(i,j)*H(i,j))- frelax(i,j)

                deltaH(i,j) = reste/crelax(i,j)

            end do
            end do

            ! Adjust H to new value
            H = H - deltaH

            ! Check stopping criterion
            delh = 0

            do j=2,ny-1
            do i=2,nx-1
                delh=delh+deltaH(i,j)**2
            end do
            end do

            if (delh.gt.0.) then
                testh=sqrt(delh)/((nx-2)*(ny-2))
            else
                testh=0.0
            end if

            stopp = (testh.lt.1.e-4).or.(ntour.gt.100)

        end do ! End of relaxation loop

        ! write(*,*) "resol_adv_diff_2D: ", stopp, testh, ntour, maxval(deltaH)

        return

    end subroutine advdiff2D



end module solver_advdiff2D
