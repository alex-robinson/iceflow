module grisli_velocity
    ! Module based on code in files `diffusiv-polyn-0.5.F90`, 
    ! `velocities-polyn-0.3.F90` and `mix-SIA-L1_mod.F90` 

    use grisli_global
    use deform_declar 
    use sliding_sia 

    use diagno_mod

    implicit none 

    type grisli_vel_param_class
        integer :: mix_method     ! Method for mixing sia and ssa velocity solutions
        integer :: ssa_domain     ! Where should ssa be calculated 
        logical :: restrict_ssa   ! Switch to apply ssa_domain limits or not
        real(4) :: H_water_ssa    ! [m] Water pressure limit above which ssa should be activated
        real(4) :: ssa_vref       ! Velocity fraction scalar
        real(4) :: Hmin, Hmax     ! Water pressure limits
        real(4) :: ulim_sia       ! Horizontal velocity limit 
    end type 

    type grisli_vel_state_class
        real(4), allocatable :: f_ssa(:,:) 
        real(4), allocatable :: f_ssa_mx(:,:), f_ssa_my(:,:)
        logical, allocatable :: ssa_active(:,:) 
        real(4), allocatable :: ux_bar(:,:), uy_bar(:,:) 
        real(4), allocatable :: ux_ssa(:,:), uy_ssa(:,:) 
        real(4), allocatable :: ux_sia(:,:), uy_sia(:,:)
        real(4), allocatable :: f_vbvs(:,:) 

        integer, allocatable :: uxbar_sign(:,:) 
        integer, allocatable :: uybar_sign(:,:)

        real(4) :: time 
        real(4), allocatable :: ux_ssa_old(:,:), uy_ssa_old(:,:)
        real(4), allocatable :: dux_ssa(:,:),    duy_ssa(:,:)
    end type 

    type grisli_vel_class 
        type(grisli_vel_param_class) :: par 

        ! State variables 
        type(grisli_vel_state_class) :: now 

    end type

    private 
    public :: grisli_vel_class
    public :: grisli_velocity_init
    public :: grisli_velocity_update_f_ssa
    public :: grisli_velocity_update
    public :: calc_ssa_mixed
    public :: limit_vel 

contains 

    subroutine grisli_velocity_init(vel,nml_in,nx,ny)

        implicit none 

        type(grisli_vel_class), intent(OUT) :: vel
        integer, intent(IN)  :: nml_in  
        integer, intent(IN)  :: nx, ny 

        ! Parameters 
        integer :: mix_method
        integer :: ssa_domain 
        logical :: restrict_ssa
        real(4) :: H_water_ssa
        real(4) :: ssa_vref, Hmin, Hmax, ulim_sia

        namelist /gvel_par/ mix_method, ssa_domain, restrict_ssa, H_water_ssa, ssa_vref, Hmin, Hmax, ulim_sia 

        ! Store initial values in local parameter values 
        mix_method   = vel%par%mix_method
        ssa_domain   = vel%par%ssa_domain
        restrict_ssa = vel%par%restrict_ssa  
        H_water_ssa  = vel%par%H_water_ssa
        ssa_vref     = vel%par%ssa_vref 
        Hmin         = vel%par%Hmin 
        Hmax         = vel%par%Hmax 
        ulim_sia     = vel%par%ulim_sia 

        ! grisli style!
        rewind(nml_in)
        read(nml_in,nml=gvel_par)

        vel%par%mix_method   = mix_method
        vel%par%ssa_domain   = ssa_domain 
        vel%par%restrict_ssa = restrict_ssa 
        vel%par%H_water_ssa  = H_water_ssa
        vel%par%ssa_vref     = ssa_vref 
        vel%par%Hmin         = Hmin 
        vel%par%Hmax         = Hmax 
        vel%par%ulim_sia     = ulim_sia 

        ! Allocate the velocity object
        call grisli_vel_allocate(vel%now,nx,ny)

        ! Set initial f_ssa to SIA only
        vel%now%f_ssa_mx = 0.0 
        vel%now%f_ssa_my = 0.0
        vel%now%f_ssa    = 0.0 

        vel%now%f_vbvs     = 0.0 

        ! Consistency checks 
        if (vel%par%mix_method .lt. -2 .or. &
            vel%par%mix_method .gt.  3) then 
            write(*,*) "grisli_velocity_init:: error: mix_method must be one of -2,-1,0,1,2,3"
            write(*,*) "mix_method = ", mix_method
            stop 
        end if 

        if (vel%par%mix_method .eq. 1) then 
            write(*,*) "grisli_velocity_init:: Note: for mix_method=1, &
                       &sliding should not be permitted in SIA (set k_weert=0!)."
        end if  

        return 

    end subroutine grisli_velocity_init

    subroutine grisli_velocity_update_f_ssa(vel,z_srf,z_bed,H_ice,H_water,f_grnd,is_grz,f_pmp,time)
        ! Calculate the fraction of ssa that should be applied in the model

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

    end subroutine grisli_velocity_update_f_ssa

    subroutine grisli_velocity_update(vel,H_sed,neff,diffusive,velocities)

        implicit none 

        type(grisli_vel_class), intent(INOUT) :: vel
        real(4), intent(IN) :: H_sed(:,:), neff(:,:)
        logical, intent(IN) :: diffusive, velocities 

        ! First call replacement for global function `diffusiv`
        if (diffusive) call calc_vel_diffusive(vel,H_sed,neff)

        ! Next call replacement for global function `velocities`
        if (velocities) call calc_velocities(vel)

        return 

    end subroutine grisli_velocity_update

    subroutine grisli_velocity_update_ssa(vel,time)

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

    end subroutine grisli_velocity_update_ssa

    subroutine calc_vel_diffusive(vel,H_sed,neff)
        ! From original grisli routine `diffusiv`
        ! Routine is used to estimate a new uxbar/uybar
        ! (ie, vertically integrated velocity fields)

        ! To do: Calculate SIA and SSA sections fully separately, store solutions,
        ! then decide what to do with them (mix them, only use sia, etc.)

        implicit none

        type(grisli_vel_class), intent(INOUT) :: vel
        real(4), intent(IN) :: H_sed(:,:), neff(:,:) 

        real(4) :: glenexp, INV_4DX, INV_4DY

        ! Original hard-coded limits for SIA velocities
!         ulgliss=5500.0
!         ultot=5899.0

        ! inverse de dx pour eviter les divisions =1/(4*dx), 1/(4*dy)
        INV_4DX=1.0/(4.0*dx)
        INV_4DY=1.0/(4.0*dy)

        ! Initialization of the diffusion
        diffmx = 0.0
        diffmy = 0.0

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

        ! === Calculate the diffusive solution (relevant for SIA) ===================
        
        ! Calcul de ddy  et ddx pour chaque valeur de iglen
        ddy = 0.0 
        ddx = 0.0 
        do iglen=n1poly,n2poly
            glenexp=max(0.0,(glen(iglen)-1.0)/2.0)
            DDY(:,:,iglen) = ((SLOPE2my**glenexp)*(ROG)**glen(iglen)) *HMY**(glen(iglen)+1)
            DDX(:,:,iglen) = ((SLOPE2mx**glenexp)*(ROG)**glen(iglen)) *HMX**(glen(iglen)+1)
        end do

        !      -------------------------------------------------
        !      GLISSEMENT
        !      DDBX et DDBY termes dus au glissement 
        !      relation avec la vitesse de glissement UXB et UYB 
        !      UXB=-DDBX*SDX et UYB=-DDBY*SDY
        !      -------------------------------------------------

        !ajr: now using sliding_sia module to determine uxb, uyb, DDBX, DDBY
        call sliding_sia_update(sliding_sia1,H_sed,neff,sdx,sdy)

        ! === Calculate deformation speed (ux_sia, uy_sia) =========================

        ! *In these calculations uxbar,uybar will be the final output for the model
        ! to use, which can be a hybrid combination of ssa, sia and sia-sliding solutions.
        ! However, below they are used as intermediate variables to calculate the individual
        ! components. 

        ! --- y-deformation ------------------

        ! Calculate diffmy, for use with diffusion, multipy with H
        ! (initially without sliding component)
        diffmy = 0.0
        do iglen = n1poly, n2poly
            diffmy = diffmy + ddy(:,:,iglen)*(s2a_my(:,:,1,iglen))
        end do
        
        ! Calculate the depth-integrated SIA velocity including sliding
        ! and separate it into basal velocity and deformational (sia) flow
        vel%now%uy_sia = diffmy*(-sdy)
        uby            = ddby  *(-sdy)

        ! Add basal sliding component to diffmy 
        diffmy = diffmy + ddby 

        ! --- x-deformation ------------------

        ! Calculate diffmx, for use with diffusion, multipy with H
        ! (initially without sliding component)
        diffmx = 0.0
        do iglen = n1poly, n2poly
            diffmx = diffmx + ddx(:,:,iglen)*(s2a_mx(:,:,1,iglen))
        end do
        
        ! Calculate the depth-integrated SIA velocity including sliding
        ! and separate it into basal velocity and deformational (sia) flow
        vel%now%ux_sia = diffmx*(-sdx)
        ubx            = ddbx  *(-sdx)

        ! Add basal sliding component to diffmx 
        diffmx = diffmx + ddbx  

        ! Limit the deformational velocity
        call limit_vel(vel%now%ux_sia,vel%par%ulim_sia)
        call limit_vel(vel%now%uy_sia,vel%par%ulim_sia)
            
!         ! ajr: test: kill all non-sia regions
!         where (vel%now%f_ssa_mx .gt. 0.0) 
!             vel%now%ux_sia = 0.0 
! !             diffmx         = 0.0      
!         end where 

!         where (vel%now%f_ssa_my .gt. 0.0) 
!             vel%now%uy_sia = 0.0 
! !             diffmy         = 0.0      
!         end where 
        
        return
        
    end subroutine calc_vel_diffusive

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

    subroutine calc_velocities(vel)
        ! Converted from velocities-polyn-0.3.F90
        ! moyenne sa_mx et sa_my
        ! Note: Attention la vitesse verticale est "super reduite"
        ! Note: With new mixing method, the above "super reduite" statement may not be true
        !       This should be checked. 

        implicit none

        type(grisli_vel_class), intent(INOUT) :: vel

        real :: hdd
        real(4) :: utmp_sia,  utmp_ssa
        real(4) :: sutmp_sia, sutmp_ssa
        real(4) :: f_ssa_now 

        do k=1,nz

            do j=2,ny-1
            do i=2,nx-1

                ! x-direction =============================================
                utmp_sia  = 0.0 
                utmp_ssa  = 0.0 
                sutmp_sia = 0.0 
                sutmp_ssa = 0.0 

                utmp_sia  = ddbx(i,j)
                sutmp_sia = ddbx(i,j)*cde(k)

                do iglen=n1poly,n2poly
                    utmp_sia  =  utmp_sia + ddx(i,j,iglen)* sa_mx(i,j,k,iglen)
                    sutmp_sia = sutmp_sia + ddx(i,j,iglen)*s2a_mx(i,j,k,iglen)
                end do
                utmp_sia  =  utmp_sia*(-sdx(i,j))
                sutmp_sia = sutmp_sia*sdx(i,j)*hmx(i,j)

                if (flotmx(i,j)) then 
                    utmp_sia   = 0.0 
                    sutmp_sia  = 0.0 
                end if 

                if (vel%now%f_ssa_mx(i,j) .gt. 0.0) then
                    ! Calculate SSA component since point is not fully SIA 

!                     utmp_ssa  = uxbar(i,j)
!                     sutmp_ssa = (nz-k)*1.0/(nz-1.0)*uxbar(i,j)*hmx(i,j)

                    utmp_ssa  = vel%now%ux_ssa(i,j)
                    sutmp_ssa = (nz-k)*1.0/(nz-1.0)*vel%now%ux_ssa(i,j)*hmx(i,j)
                    
                end if 

                ! Combine the calculations to obtain hybrid values
                ux(i,j,k)  = calc_ssa_mixed(utmp_sia,utmp_ssa,vel%now%f_ssa_mx(i,j),vel%par%mix_method)
!                 sux(i,j,k) = calc_ssa_mixed(sutmp_sia,sutmp_ssa,vel%now%f_ssa_mx(i,j),vel%par%mix_method)
                sux(i,j,k) = sutmp_sia 

                if (k .eq. nz) then 
                    ubx_sia(i,j) = utmp_sia 
                    ubx_ssa(i,j) = utmp_ssa 
                    ubx(i,j)     = ux(i,j,k) 
                end if 

                if (k .eq. 1) then 
                    if (utmp_sia .eq. 0.0 .or. utmp_ssa .eq. 0.0) then 
                        vel%now%uxbar_sign(i,j) = 1.0 
                    else if (utmp_sia .gt. 0.0 .and. utmp_ssa .gt. 0.0) then 
                        vel%now%uxbar_sign(i,j) = 1.0 
                    else if (utmp_sia .lt. 0.0 .and. utmp_ssa .lt. 0.0) then 
                        vel%now%uxbar_sign(i,j) = 1.0 
                    else 
                        vel%now%uxbar_sign(i,j) = 0.0 
                    end if 

!                     vel%now%uxbar_sign(i,j) = sign(1.0,utmp_sia)*sign(1.0,utmp_ssa)
                end if 

                ! y-direction =============================================
                utmp_sia  = 0.0 
                utmp_ssa  = 0.0 
                sutmp_sia = 0.0 
                sutmp_ssa = 0.0 

                utmp_sia  = ddby(i,j)
                sutmp_sia = ddby(i,j)*cde(k)

                do iglen=n1poly,n2poly
                    utmp_sia  =  utmp_sia + ddy(i,j,iglen)* sa_my(i,j,k,iglen)
                    sutmp_sia = sutmp_sia + ddy(i,j,iglen)*s2a_my(i,j,k,iglen)
                end do
                utmp_sia  =  utmp_sia*(-sdy(i,j))
                sutmp_sia = sutmp_sia*sdy(i,j)*hmy(i,j)

                if (flotmy(i,j)) then 
                    utmp_sia   = 0.0 
                    sutmp_sia  = 0.0 
                end if 

                if (vel%now%f_ssa_my(i,j) .gt. 0.0) then
                    ! Calculate SSA component since point is not fully SIA 

!                     utmp_ssa  = uybar(i,j)
!                     sutmp_ssa = (nz-k)*1.0/(nz-1.0)*uybar(i,j)*hmy(i,j)

                    utmp_ssa  = vel%now%uy_ssa(i,j)
                    sutmp_ssa = (nz-k)*1.0/(nz-1.0)*vel%now%uy_ssa(i,j)*hmy(i,j)
                    
                end if 

                ! Combine the calculations to obtain hybrid values
                uy(i,j,k)  = calc_ssa_mixed(utmp_sia,utmp_ssa,vel%now%f_ssa_my(i,j),vel%par%mix_method)
!                 suy(i,j,k) = calc_ssa_mixed(sutmp_sia,sutmp_ssa,vel%now%f_ssa_my(i,j),vel%par%mix_method)
                suy(i,j,k) = sutmp_sia 

                if (k .eq. nz) then 
                    uby_sia(i,j) = utmp_sia 
                    uby_ssa(i,j) = utmp_ssa 
                    uby(i,j)     = uy(i,j,k) 
                end if 

                if (k .eq. 1) then 
                    if (utmp_sia .eq. 0.0 .or. utmp_ssa .eq. 0.0) then 
                        vel%now%uybar_sign(i,j) = 1.0 
                    else if (utmp_sia .gt. 0.0 .and. utmp_ssa .gt. 0.0) then 
                        vel%now%uybar_sign(i,j) = 1.0 
                    else if (utmp_sia .lt. 0.0 .and. utmp_ssa .lt. 0.0) then 
                        vel%now%uybar_sign(i,j) = 1.0 
                    else 
                        vel%now%uybar_sign(i,j) = 0.0 
                    end if 

!                     vel%now%uybar_sign(i,j) = sign(1.0,utmp_sia)*sign(1.0,utmp_ssa)
                end if 

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

        ! Calculate the basal to surface velocity ratio, f_vbvs
        where ( ux(:,:,1)**2+uy(:,:,1)**2 .gt. 0.0) 
            vel%now%f_vbvs = ((ux(:,:,nz)**2+uy(:,:,nz)**2)**0.5 &
                                            / (ux(:,:,1)**2+uy(:,:,1)**2)**0.5)
        elsewhere 
            ! No ice (or no velocity)
            vel%now%f_vbvs = 1.0 
        end where 

        return 

    end subroutine calc_velocities

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
                if (H_water .eq. 0.0 .and. z_srf .gt. 3000.0 .and. H_ice .gt. 1000) ssa_active = .FALSE. 

            case DEFAULT 

                ! SSA remains active by default 
                ssa_active = .TRUE. 

        end select 

        return 

    end function define_ssa_active

    ! === OBJECT MANAGEMENT ==============================


    subroutine grisli_vel_allocate(now,nx,ny)

        implicit none 

        type(grisli_vel_state_class) :: now 
        integer :: nx, ny 

        call grisli_vel_deallocate(now)

        allocate(now%f_ssa(nx,ny))
        allocate(now%f_ssa_mx(nx,ny),now%f_ssa_my(nx,ny))
        allocate(now%ssa_active(nx,ny))
        allocate(now%ux_bar(nx,ny),now%uy_bar(nx,ny))
        allocate(now%ux_ssa(nx,ny),now%uy_ssa(nx,ny))
        allocate(now%ux_sia(nx,ny),now%uy_sia(nx,ny))

        allocate(now%f_vbvs(nx,ny))
        allocate(now%uxbar_sign(nx,ny),now%uybar_sign(nx,ny))

        allocate(now%ux_ssa_old(nx,ny),now%uy_ssa_old(nx,ny))
        allocate(now%dux_ssa(nx,ny),now%duy_ssa(nx,ny))
        

        ! Set all variables to zero intially 
        now%f_ssa    = 0.0
        now%f_ssa_mx = 0.0
        now%f_ssa_my = 0.0
        now%ux_bar   = 0.0 
        now%uy_bar   = 0.0 
        now%ux_ssa   = 0.0 
        now%uy_ssa   = 0.0 
        now%ux_sia   = 0.0 
        now%uy_sia   = 0.0 
        
        now%f_vbvs   = 0.0 

        now%ux_ssa_old = 0.0 
        now%uy_ssa_old = 0.0 
        now%dux_ssa    = 0.0 
        now%duy_ssa    = 0.0 
        
        now%ssa_active = .TRUE. 

        return 

    end subroutine grisli_vel_allocate 

    subroutine grisli_vel_deallocate(now)

        implicit none 

        type(grisli_vel_state_class) :: now 

        if (allocated(now%f_ssa))      deallocate(now%f_ssa)
        if (allocated(now%f_ssa_mx))   deallocate(now%f_ssa_mx)
        if (allocated(now%f_ssa_my))   deallocate(now%f_ssa_my)
        if (allocated(now%ssa_active)) deallocate(now%ssa_active)
        if (allocated(now%ux_bar))     deallocate(now%ux_bar)
        if (allocated(now%uy_bar))     deallocate(now%uy_bar)
        if (allocated(now%ux_ssa))     deallocate(now%ux_ssa)
        if (allocated(now%uy_ssa))     deallocate(now%uy_ssa)
        if (allocated(now%ux_sia))     deallocate(now%ux_sia)
        if (allocated(now%uy_sia))     deallocate(now%uy_sia)
        
        if (allocated(now%f_vbvs))     deallocate(now%f_vbvs)
        if (allocated(now%uxbar_sign))     deallocate(now%uxbar_sign)
        if (allocated(now%uybar_sign))     deallocate(now%uybar_sign)
        
        if (allocated(now%ux_ssa_old)) deallocate(now%ux_ssa_old)
        if (allocated(now%uy_ssa_old)) deallocate(now%uy_ssa_old)
        if (allocated(now%dux_ssa))    deallocate(now%dux_ssa)
        if (allocated(now%duy_ssa))    deallocate(now%duy_ssa)
        
        return 

    end subroutine grisli_vel_deallocate 
    
end module grisli_velocity 
