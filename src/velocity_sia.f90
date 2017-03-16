module velocity_sia

    implicit none
    
    
    
    
    private
    public :: calc_uv_sia
    
contains

    subroutine calc_uv_sia(par,z_srf,H_ice,At_int,dx,dy)
        ! Calculate the depth integrated horizontal velocity field
        ! using the SIA approximation
        
        implicit none
        
        implicit none

        type(grisli_vel_class), intent(INOUT) :: vel
        real(4), intent(IN) :: z_srf(:,:), H_ice(:,:)
        real(4), intent(IN) :: At_int(:,:)             ! Depth-integrated rate factor
        real(4), intent(IN) :: dx, dy 
        
        ! Local variables
        real(4) :: glenexp
        real(4), parameter :: INV_4DX = 1.0/(4.0*dx)
        real(4), parameter :: INV_4DY = 1.0/(4.0*dy)

        
        ! 1. Calculate surface slopes =============================
        
        do j=1,ny
        do i=2,nx
            sdx(i,j)=(s(i,j)-s(i-1,j))/dx
        end do
        end do

        do j=2,ny
        do i=1,nx
            sdy(i,j)=(s(i,j)-s(i,j-1))/dy
        end do
        end do

        do j=2,ny
        do i=2,nx-1
            sdxmy(i,j)= ((s(i+1,j)-s(i-1,j))+s(i+1,j-1)-s(i-1,j-1))*inv_4dx
        end do
        end do

        do j=2,ny-1
        do i=2,nx
            sdymx(i,j)= ((s(i,j+1)-s(i,j-1))+s(i-1,j+1)-s(i-1,j-1))*inv_4dy
        enddo
        enddo

        slope2mx = sdx**2 + sdymx**2
        slope2my = sdy**2 + sdxmy**2

        slope2mx(1,:)  = 0.0
        slope2mx(:,1)  = 0.0
        slope2mx(:,ny) = 0.0

        slope2my(1,:)  = 0.0
        slope2my(:,1)  = 0.0
        slope2my(nx,:) = 0.0
        
        
        ! 2. Calculate SIA sliding via sliding module ===============
        
        call sliding_sia_update(sliding_sia1,sed1%now%H,sdx,sdy)
        
        
        ! 3. Calculate the SIA diffusive solution ===================

        glenexp=max(0.0,(par%e_glen-1.0)/2.0)

        ddy = ((slope2my**glenexp)*(ROG)**par%e_glen) *HMY**(par%e_glen+1)
        ddx = ((slope2mx**glenexp)*(ROG)**par%e_glen) *HMX**(par%e_glen+1)


        ! 4. Calculate deformation speed (ux_sia, uy_sia) =========================

        ! --- y-deformation ------------------

        ! Calculate diffmy, for use with diffusion, multipy with H
        ! (initially without sliding component)
        diffmy = ddy*s2a_my(:,:,1)
        
        ! Calculate the depth-integrated SIA velocity (without sliding) and basal velocity
        vel%now%uy_sia = diffmy*(-sdy)
        uby            = ddby  *(-sdy)

        ! Add basal sliding component to diffmy 
        diffmy = diffmy + ddby   

        
        ! --- x-deformation ------------------

        ! Calculate diffmx, for use with diffusion, multipy with H
        ! (initially without sliding component)
        diffmx = ddx*s2a_mx(:,:,1)
        
        ! Calculate the depth-integrated SIA velocity (without sliding) and basal velocity
        vel%now%ux_sia = diffmx*(-sdx)
        ubx            = ddbx  *(-sdx)

        ! Add basal sliding component to diffmx 
        diffmx = diffmx + ddbx 
        
        
        ! 5. Limit the deformational velocity =====================
        ! ajr: this needs revision as diffmx/my will subsequently be inconsistent
        
        call limit_vel(vel%now%ux_sia,vel%par%ulim_sia)
        call limit_vel(vel%now%uy_sia,vel%par%ulim_sia)
        
        return
        
    end subroutine calc_uv_sia
    
    
    subroutine calc_vel3D_sia()
        ! Calculate the 3D velocity field
        ! using the SIA approximation
        
        implicit none
        
        
        
        
        return
        
    end subroutine calc_vel3D_sia




end module velocity_sia
