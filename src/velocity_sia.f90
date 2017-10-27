module velocity_sia
    
    !use iceflow_types
    
    implicit none
    
    real(4), parameter :: rho_ice = 910.0 
    real(4), parameter :: G       = 9.8 
    real(4), parameter :: rhog    = rho_ice*G 

    private
    public :: calc_vxy_sia
    
contains
    
    subroutine calc_vxy_sia(ux_sia,uy_sia,diff_mx,diff_my,ddx,ddy,z_srf,H_mx,H_my, &
                                ATT_int_srf_mx,ATT_int_srf_my,is_float,dx,dy,e_glen,ulim)
        ! Calculate the depth integrated horizontal velocity field
        ! and intermediate variables using the SIA approximation
        
        implicit none

        real(4), intent(OUT) :: ux_sia(:,:), uy_sia(:,:)
        real(4), intent(OUT) :: diff_mx(:,:), diff_my(:,:)
        real(4), intent(OUT) :: ddx(:,:), ddy(:,:)
        real(4), intent(IN) :: z_srf(:,:)
        real(4), intent(IN) :: H_mx(:,:), H_my(:,:)
        real(4), intent(IN) :: ATT_int_srf_mx(:,:), ATT_int_srf_my(:,:)     ! Depth-integrated rate factor (ie, at surface) - called s2a_mx in grisli
        ! note: event. staggared values should be determined here locally
        logical, intent(IN) :: is_float(:,:)
        real(4), intent(IN) :: dx, dy, e_glen, ulim 
        
        ! Local variables
        real(4), allocatable :: sdx(:,:), sdy(:,:), sdx_my(:,:), sdy_mx(:,:)
        real(4), allocatable :: slope_mx(:,:), slope_my(:,:)
        real(4) :: inv_4dx, inv_4dy 
        integer :: i, j, nx, ny
        
        inv_4dx = 1.0/(4.0*dx)
        inv_4dy = 1.0/(4.0*dy)
        
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
        
            ddx = slope_mx**(e_glen-1.0) * (rhog)**e_glen * H_mx**(e_glen+1.0)
            
            ! Calculate diffmx, for use with diffusion, multipy with H
            ! (initially without sliding component)
            diff_mx = ddx*ATT_int_srf_mx
            
            ! Calculate the depth-integrated SIA velocity (without sliding) and basal velocity
            ux_sia = diff_mx*(-sdx)
            
        elsewhere
        
            ddx     = 0.0
            diff_mx = 0.0
            ux_sia  = 0.0
            
        end where
        
        ! --- y-deformation ------------------
        
        where (.not. is_float)
        
            ddy = slope_my**(e_glen-1.0) * (rhog)**e_glen * H_my**(e_glen+1.0)
            
            ! Calculate diffmy, for use with diffusion, multipy with H
            ! (initially without sliding component)
            diff_my = ddy*ATT_int_srf_my
            
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
        call limit_vel(ux_sia,ulim)
        call limit_vel(uy_sia,ulim)
        
        return
        
    end subroutine calc_vxy_sia

    elemental subroutine limit_vel(u,u_lim)
        ! Apply a velocity limit to points with an SIA contribution
        ! (for stability)

        implicit none 

        real(4), intent(INOUT) :: u 
        real(4), intent(IN)    :: u_lim

        u = min(u, u_lim)
        u = max(u,-u_lim)

        return 

    end subroutine limit_vel


end module velocity_sia
