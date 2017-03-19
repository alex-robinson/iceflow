module velocity_3D

    implicit none
    
    private
    public :: calc_vel3D
    
contains
    
    
    subroutine calc_vel3D(ux,uy,uz)
        ! Calculate the 3D velocity field
        ! using sia, ssa or hybrid method
        
        implicit none
        
        real(4), intent(OUT) :: ux(:,:,:)
        real(4), intent(OUT) :: uy(:,:,:)
        real(4), intent(OUT) :: uz(:,:,:)
        
        type(grisli_vel_class), intent(INOUT) :: vel

        real :: hdd
        real(4) :: utmp_sia,  utmp_ssa
        real(4) :: sutmp_sia, sutmp_ssa
        real(4) :: f_ssa_now 

        do k=1,nz

            do j=2,ny-1
            do i=2,nx-1

                ! x-direction =============================================
                
                utmp_sia  = ddx(i,j)* sa_mx(i,j,k) *(-sdx(i,j))
                sutmp_sia = ddx(i,j)*s2a_mx(i,j,k) *sdx(i,j)*hmx(i,j)

                utmp_ssa  = vel%now%ux_ssa(i,j)
                sutmp_ssa = (nz-k)*1.0/(nz-1.0)*vel%now%ux_ssa(i,j)*hmx(i,j)

                ! Combine the calculations to obtain hybrid values
                ux(i,j,k)  = calc_ssa_mixed(utmp_sia,utmp_ssa,vel%now%f_ssa_mx(i,j),vel%par%mix_method)
!                 sux(i,j,k) = calc_ssa_mixed(sutmp_sia,sutmp_ssa,vel%now%f_ssa_mx(i,j),vel%par%mix_method)
                sux(i,j,k) = sutmp_sia 

                ! y-direction =============================================
                
                utmp_sia  = ddy(i,j)* sa_my(i,j,k) *(-sdy(i,j))
                sutmp_sia = ddy(i,j)*s2a_my(i,j,k) *sdy(i,j)*hmy(i,j)

                utmp_ssa  = vel%now%uy_ssa(i,j)
                sutmp_ssa = (nz-k)*1.0/(nz-1.0)*vel%now%uy_ssa(i,j)*hmy(i,j)

                ! Combine the calculations to obtain hybrid values
                uy(i,j,k)  = calc_ssa_mixed(utmp_sia,utmp_ssa,vel%now%f_ssa_my(i,j),vel%par%mix_method)
!                 suy(i,j,k) = calc_ssa_mixed(sutmp_sia,sutmp_ssa,vel%now%f_ssa_my(i,j),vel%par%mix_method)
                suy(i,j,k) = sutmp_sia 

            end do
            end do
        end do
        
        ! *************************** Z VELOCITIES ******************

        do j=2,ny-1
        do i=2,nx-1

            divu(i,j) = ((uxbar(i+1,j)*hmx(i+1,j)-uxbar(i,j)*hmx(i,j))  &
                       + (uybar(i,j+1)*hmy(i,j+1)-uybar(i,j)*hmy(i,j)))/dx

            uzr(i,j,nz) = bmelt(i,j)
            xx(i,j)     = uzr(i,j,1)
            hdd         = bm(i,j)-bmelt(i,j)-divu(i,j)

            ! -------------------------------------------------------------------
            ! attention uzr contient maintenant :
            ! uzr=uz +u ((1-e)*dH/dx+dB/dx)+ ((1-e)*dH/dt+dB/dt)

            do k=1,nz-1
                
                ! "new-mix-vert"
                utmp_sia = 0.0 
                utmp_ssa = 0.0

                ! NOTE: The "centered" f_ssa map is used here (and only here), since 
                ! the vertical velocity is a function of divergence which includes
                ! both the x- and y- components of velocity 

                utmp_sia = uzr(i,j,nz)-((sux(i+1,j,k)-sux(i,j,k)) &
                               + (suy(i,j+1,k)-suy(i,j,k)))/dx + hdd*cde(k)

!                 if (vel%now%f_ssa(i,j) .gt. 0.0) then 
!                     ! Calculate SSA component since point is not fully SIA 

                    utmp_ssa = bm(i,j)+(k-1.0)/(nz-1.0)*(bmelt(i,j)-bm(i,j))
                
!                 end if

                ! Combine the calculations to obtain hybrid values
                uzr(i,j,k)  = calc_ssa_mixed(utmp_sia,utmp_ssa,vel%now%f_ssa(i,j),vel%par%mix_method)
                
            end do
            
        end do
        end do

        ! === Calculate diagnostic quantities ====

        ! ajr: the below calculation is only diagnostic, but should be fixed to properly
        ! handle the dt_ssa timestep in place of dtt 
!         ! Determine rates of change
!         uzsdot = (uzr(:,:,1) - xx) / dtt
        uzk    = -bdot - hdot + bm - bmelt
        
        
        return
        
    end subroutine calc_vel3D


    

end module velocity_3D
