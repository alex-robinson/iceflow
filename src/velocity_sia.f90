module velocity_sia
    
    use iceflow_types
    
    implicit none
    
    
    
    
    private
    public :: calc_uv_sia
    
contains

    subroutine calc_diff_sia(ux_sia,uy_sia,diff_mx,diff_my,ddx,ddy,z_srf,H_mx,H_my, &
                                At_int_mx,At_int_my,is_float,dx,dy,e_glen)
        ! Calculate the depth integrated horizontal velocity field
        ! and intermediate variables using the SIA approximation
        
        implicit none

        real(4), intent(OUT) :: ux_sia(:,:), uy_sia(:,:)
        real(4), intent(OUT) :: diff_mx(:,:), diff_mx(:,:)
        real(4), intent(OUT) :: ddx(:,:), ddx(:,:)
        real(4), intent(IN) :: z_srf(:,:)
        real(4), intent(IN) :: H_mx(:,:), H_my(:,:)
        real(4), intent(IN) :: At_int_mx(:,:), At_int_my             ! Depth-integrated rate factor (ie, at surface) - called s2a_mx in grisli
        ! note: event. staggared values should be determined here locally
        logical, intent(IN) :: is_float(:,:)
        real(4), intent(IN) :: dx, dy, e_glen
        
        ! Local variables
        real(4), parameter :: INV_4DX = 1.0/(4.0*dx)
        real(4), parameter :: INV_4DY = 1.0/(4.0*dy)
        real(4), allocatable :: sdx(:,:), sdy(:,:), sdx_my(:,:), sdy_mx(:,:)
        real(4), allocatable :: slope_mx(:,:), slope_my(:,:)
        integer :: i, j, nx, ny
        
        nx = size(z_srf,1)
        ny = size(z_srf,2)
        
        allocate(sdx(nx,ny),sdy(nx,ny))
        allocate(sdx_my(nx,ny),sdy_mx(nx,ny))
        allocate(slope_mx(nx,ny),slope_my(nx,ny))
        
        
        ! 1. Calculate surface slopes =============================
        
        do j=1,ny
        do i=2,nx
            sdx(i,j)=(z_srf(i,j)-z_srf(i-1,j))/dx
        end do
        end do

        do j=2,ny
        do i=1,nx
            sdy(i,j)=(z_srf(i,j)-z_srf(i,j-1))/dy
        end do
        end do

        do j=2,ny
        do i=2,nx-1
            sdx_my(i,j)= ((z_srf(i+1,j)-z_srf(i-1,j))+z_srf(i+1,j-1)-z_srf(i-1,j-1))*inv_4dx
        end do
        end do

        do j=2,ny-1
        do i=2,nx
            sdy_mx(i,j)= ((z_srf(i,j+1)-z_srf(i,j-1))+z_srf(i-1,j+1)-z_srf(i-1,j-1))*inv_4dy
        enddo
        enddo

        slope_mx = (sdx**2 + sdy_mx**2)**0.5
        slope_my = (sdy**2 + sdx_my**2)**0.5

        slope_mx(1,:)  = 0.0
        slope_mx(:,1)  = 0.0
        slope_mx(:,ny) = 0.0

        slope_my(1,:)  = 0.0
        slope_my(:,1)  = 0.0
        slope_my(nx,:) = 0.0
        
        ! 3. Calculate the 2D SIA diffusive solution ===================
        
        ! --- x-deformation ------------------
        
        where (.not. is_float)
        
            ddx = slope_mx**(e_glen-1.0) * (ROG)**e_glen * H_mx**(e_glen+1.0)
            
            ! Calculate diffmx, for use with diffusion, multipy with H
            ! (initially without sliding component)
            diff_mx = ddx*At_int_mx
            
            ! Calculate the depth-integrated SIA velocity (without sliding) and basal velocity
            ux_sia = diff_mx*(-sdx)
            
        elsewhere
        
            ddx    = 0.0
            diffmx = 0.0
            ux_sia = 0.0
            
        end where
        
        ! --- y-deformation ------------------
        
        where (.not. is_float)
        
            ddy = slope_my**(e_glen-1.0) * (ROG)**e_glen * Hmy**(e_glen+1.0)
            
            ! Calculate diffmy, for use with diffusion, multipy with H
            ! (initially without sliding component)
            diff_my = ddy*At_int_my
            
            ! Calculate the depth-integrated SIA velocity (without sliding) and basal velocity
            uy_sia = diff_my*(-sdy)
            
        elsewhere
        
            ddy     = 0.0
            diff_my = 0.0
            uy_sia  = 0.0
            
        end where
        
        ! 5. Limit the deformational velocity =====================
        ! ajr: this needs revision as diffmx/my will subsequently be inconsistent
        ! Better: either limit slope in diff calculation, or move this outside of routine
        !call limit_vel(vel%now%ux_sia,vel%par%ulim_sia)
        !call limit_vel(vel%now%uy_sia,vel%par%ulim_sia)
        
        return
        
    end subroutine calc_uv_sia
    
    
    subroutine calc_vel3D_sia()
        ! Calculate the 3D velocity field
        ! using the SIA approximation
        
        implicit none
        
        
        
        
        return
        
    end subroutine calc_vel3D_sia




end module velocity_sia
