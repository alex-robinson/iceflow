module mixing
    ! Mix variables according to hybrid solver assumptions
    
    implicit none
    
    
    
    
    private
    public :: calc_ssa_mixed
    public :: calc_ssa_fraction
    public :: calc_ssa_fraction_centered
    public :: smooth_f_ssa
    public :: define_ssa_active
    
contains
    
    elemental function calc_ssa_mixed(var_sia,var_ssa,f_ssa,mix_method) result(var)
        ! Calculate the total velocity as a mixture of SIA and SSA velocities 
        ! Also used for mixing other variables (heat production, bmelt, etc)

        implicit none 

        real(4), intent(IN) :: var_sia, var_ssa, f_ssa
        integer, intent(IN) :: mix_method 
        real(4) :: var 

        select case(mix_method)

            case(1)
                ! Simple addition of sia and ssa velocities everywhere
                ! (No if-statement needed since var_ssa is zero where f_ssa==0)
                var = var_ssa + var_sia   

            case(-2,-1,0,2,3)
                ! Weighted average of sia and ssa velocities
                ! case -2, -1 and 0: binary cases (f_ssa eq 0 or 1), ie no real mixing
                ! case 2: f_ssa determined by v_ssa 
                ! case 3: f_ssa determined by Hwat 
                var = var_ssa*f_ssa+var_sia*(1.0-f_ssa)

            case DEFAULT 

                ! Pass - nothing happens 
                ! This case should be caught as an error during parameter loading

        end select 

        return 

    end function calc_ssa_mixed

    subroutine calc_ssa_fraction(par,f_ssa,f_grnd,is_grz,u_bar,z_srf,H_water,f_pmp,time)
        ! This fraction subroutine is called separately for the x- and y-direction
        ! fraction map for the respective staggered grids 

        implicit none 

        type(grisli_vel_param_class), intent(IN) :: par
        real(4), intent(OUT) :: f_ssa(:,:)
        logical, intent(IN)  :: is_grz(:,:)
        real(4), intent(IN)  :: f_grnd(:,:), u_bar(:,:), z_srf(:,:), H_water(:,:), f_pmp(:,:)
        real(4), intent(IN)  :: time 

        ! Local variables 
        real(4) :: z_srf_lim 

        ! First determine which areas should be streaming
        
        ! Calculate the fraction of sia or ssa that should be applied in the model
        select case(par%mix_method)

            case(-2)
                ! SIA is used everywhere inland, use ssa on floating ice shelves
                f_ssa = 0.0 
                where (f_grnd .eq. 0.0) f_ssa = 1.0 

            case(-1)
                ! SSA is used everywhere, so use ssa everywhere 
                f_ssa = 1.0 

            case(0) 
                ! Assign ssa points according to these conditions:
                ! 1. Either H_water is above minimum threshold,
                ! 2. or the point is floating or partially floating
                ! 3. Point is at the grounding line, or is a neighbor of the grounding line

                f_ssa = 0.0 

!                where (H_water .gt. par%H_water_ssa .or. f_grnd .lt. 1.0 .or. is_grz) ! Floating or not fully grounded, or grounding zone
!                    f_ssa = 1.0
!                end where 

!                 where (f_pmp .gt. 0.5 .or. f_grnd .lt. 1.0 .or. is_grz) ! Floating or not fully grounded, or grounding zone
!                     f_ssa = 1.0
!                 end where

                z_srf_lim = 2000.0 + 500.0*sin(2*pi*(time*1e-3)/10.0)

                where (z_srf .lt. z_srf_lim)
                    f_ssa = 1.0 
                elsewhere 
                    f_ssa = 0.0 
                end where 

                ! Smooth out the f_ssa field to avoid sharp transitions 
                !call smooth_f_ssa(f_ssa)

            case(1)
                ! SSA is used everywhere (for now), so calculate it everywhere 
                ! SSA and SIA will be summed via mixing.
                f_ssa = 1.0 

                ! To do: limit ssa to regions with effective pressure below a threshold 

!                 ! Smooth out the f_ssa field to avoid sharp transitions 
!                 call smooth_f_ssa(f_ssa)
                
            case(2)
                ! Fraction as a function of vertically-integrated velocity
                f_ssa = calc_fraction_v(par,f_grnd,is_grz,u_bar)

!                 ! Smooth out the f_ssa field to avoid sharp transitions 
!                 call smooth_f_ssa(f_ssa)
                
            case(3)
                ! Fraction as a function of basal water pressure
                f_ssa = calc_fraction_hwat(par,f_grnd,is_grz,H_water)
                
!                 ! Smooth out the f_ssa field to avoid sharp transitions 
!                 call smooth_f_ssa(f_ssa)
                
            case DEFAULT 
                write(*,*) "grisli_velocity_fraction:: error: &
                            &chosen mix_method not recognized: ", par%mix_method 
                stop 

        end select 

!         ! Limit fraction to values above a limit (eg, 1%)
!         ! to reduce computations
!         where(f_ssa .lt. 0.01) f_ssa = 0.0 

        return 

    end subroutine calc_ssa_fraction

    subroutine calc_ssa_fraction_centered(par,f_ssa,f_ssa_mx,f_ssa_my)
        ! This fraction subroutine is called separately for the x- and y-direction
        ! fraction map for the respective staggered grids 

        implicit none 

        type(grisli_vel_param_class), intent(IN) :: par
        real(4), intent(OUT) :: f_ssa(:,:)
        real(4), intent(IN)  :: f_ssa_mx(:,:), f_ssa_my(:,:)

        ! Local variables
        real(4), allocatable :: tmp1(:,:), tmp2(:,:)

        allocate(tmp1(size(f_ssa,1),size(f_ssa,2)))
        allocate(tmp2(size(f_ssa,1),size(f_ssa,2)))
        
        tmp1 = f_ssa_mx 
        tmp2 = f_ssa_my 

        ! Get mean of centered component matrices 
        f_ssa = (tmp1+tmp2) * 0.5 

!         ! Correct the central matrix (method specific)
!         select case(par%mix_method)

!             case(-2,-1,0,1)
!                 ! Binary metric, so ensure only values of 0,1 in central matrix
!                 where(f_ssa .gt. 0.0) f_ssa = 1.0 

!         end select 

        return 

    end subroutine calc_ssa_fraction_centered

    elemental function calc_fraction_v(par,f_grnd,is_grz,u_bar) result(f)
        ! Calculate the value of the fraction of velocity of SSA and SIA

        implicit none

        type(grisli_vel_param_class), intent(IN) :: par
        real(4), intent(IN)  :: f_grnd
        logical, intent(IN)  :: is_grz
        real(4), intent(IN)  :: u_bar   !ux_bar, uy_bar
        real(4) :: f,c,v,v1
        real(4), parameter :: pi = 3.1415927

        c = par%ssa_vref ! Calibration term (reference ssa velocity, ~100 m/a)

        ! Initially set fraction to zero (all SIA)
        f = 0.0
        v = 0.0

        if ( (f_grnd .lt. 1.0 .or. is_grz) ) then
            ! Floating ice shelf or grounding zone region
            f = 1.0

        else
            ! Inland ice stream
!             v = sqrt(ux_ssa**2 + uy_ssa**2)
!             v = u_ssa 
            v = u_bar

!             ! Bueler-Brown (2009) equation 
!             f = (2/pi)*atan((v/c)**2)

            ! New approach (explicitly reaches 0-1 values at extremes)
            v1 = min(abs(v)/c,1.0)
            f = 1.0-(0.5*cos(pi*v1)+0.5)

        end if

        return

    end function calc_fraction_v

    elemental function calc_fraction_hwat(par,f_grnd,is_grz,H_water) result(f)
        ! Calculate the value of the fraction of velocity of SSA and SIA

        implicit none

        type(grisli_vel_param_class), intent(IN) :: par
        real(4), intent(IN) :: H_water
        real(4), intent(IN)  :: f_grnd
        logical, intent(IN)  :: is_grz
        real(4) :: f
        
        if ( (f_grnd .lt. 1.0 .or. is_grz) ) then
            ! Floating ice shelf or grounding zone region
            f = 1.0
        
        else if (H_water .gt. par%Hmin) then
            ! Inland ice stream
            f = min(1.0,(H_water-par%Hmin)/(par%Hmax-par%Hmin))
        
        else 
            ! Not enough water available
            f = 0.0

        end if

        return

    end function calc_fraction_hwat

    subroutine smooth_f_ssa(f_ssa)
        ! Remove isolated non-ssa points from within an ssa region
        ! Also, eliminate sharp boundaries from f_ssa==1 to f_ssa==0

        implicit none 

        real(4), intent(INOUT) :: f_ssa(:,:)
        
        ! Local variables 
        integer :: i, j, nx, ny 
        real(4), allocatable :: f_ssa_old(:,:) 

        nx = size(f_ssa,1)
        ny = size(f_ssa,2)

        allocate(f_ssa_old(nx,ny))

        ! Note: step 2 should achieve both goals of eliminating isolated points
        ! and eliminating sharp boundaries, so only perform step 2
        
!         ! 1. Eliminate isolated points
!         f_ssa_old = f_ssa 

!         do i = 2, nx-1 
!         do j = 2, ny-1 

!             if (f_ssa_old(i,j) .eq. 0.0) then 
!                 ! Check inactive points

!                 if (count(f_ssa_old(i-1:i+1,j-1:j+1) .gt. 0.0) .gt. 5) then 
!                     ! Neighborhood is dominated by active ssa points,
!                     ! so treat it as ssa as well for now 
!                     f_ssa(i,j) = sum(f_ssa_old(i-1:i+1,j-1:j+1))/9.0

!                 end if 

!             end if 

!         end do 
!         end do 
        
        ! 2. Eliminate sharp boundaries 
        f_ssa_old = f_ssa 

        do i = 3, nx-2 
        do j = 3, ny-2 

            if (f_ssa_old(i,j) .eq. 0.0) then 
                ! Check inactive points

                if (count(f_ssa_old(i-2:i+2,j-2:j+2) .eq. 1.0) .gt. 0) then 
                    ! Neighborhood contains binary transition,
                    ! so smooth it out
                    f_ssa(i,j) = sum(f_ssa_old(i-2:i+2,j-2:j+2))/25.0

                end if 

            end if 

        end do 
        end do 

        return 

    end subroutine smooth_f_ssa

    elemental function define_ssa_active(par,z_srf,z_bed,H_ice,H_water,f_grnd,f_ssa) result(ssa_active)

        implicit none 

        type(grisli_vel_param_class), intent(IN) :: par 
        real(4), intent(IN) :: z_srf, z_bed, H_ice, H_water 
        real(4), intent(IN) :: f_grnd, f_ssa 
        logical :: ssa_active

        select case(par%ssa_domain) 

            case(0)

                ssa_active = .TRUE. 

            case(1)
                ssa_active = .TRUE. 

                ! Limit ssa regions based on Jorge's cool criteria 
                if (H_water .eq. 0.0 .and. z_srf .gt. 3000.0 .and. H_ice .gt. 1000.0) ssa_active = .FALSE. 

            case DEFAULT 

                ! SSA remains active by default 
                ssa_active = .TRUE. 

        end select 

        return 

    end function define_ssa_active

end module mixing
