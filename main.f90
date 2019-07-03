!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Main program
!
! - Simulate emissions and chemical reactions of gases, aerosol processes as well as 
!   transport of gases and aerosol particles within the planetary boundary layer with a
!   column model.
! - Check Fortran conventions at http://www.fortran90.org/src/best-practices.html
! - Check code conventions at
!   http://www.cesm.ucar.edu/working_groups/Software/dev_guide/dev_guide/node7.html
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

program main

use chemistry_mod
use aerosol_mod

implicit none

!-----------------------------------------------------------------------------------------
! Control variables (can be moved to an input file in future)
!-----------------------------------------------------------------------------------------
logical :: use_emission   = .true.
logical :: use_chemistry  = .true.
logical :: use_deposition = .false.
logical :: use_aerosol    = .true.
character(len=255), parameter :: input_dir  = './input'
character(len=255), parameter :: output_dir = './output'

!-----------------------------------------------------------------------------------------
! Constants
!-----------------------------------------------------------------------------------------
! Double precision
! http://en.wikipedia.org/wiki/Double_precision_floating-point_format
!integer, parameter :: dp = selected_real_kind(15, 307)

! Physics constants
real(dp), parameter :: PI     = 2*asin(1.0_dp)                  ! the constant pi
real(dp), parameter :: grav   = 9.81_dp                         ! [m s-2], gravitation
real(dp), parameter :: Rgas   = 8.3144598_dp                    ! [J mol-1 K-1], universal gas constant
real(dp), parameter :: NA     = 6.022140857e23_dp               ! [molec mol-1], Avogadro's number 
real(dp), parameter :: mm_air = 28.96e-3_dp                     ! [kg mol-1], mean molar mass of air
real(dp), parameter :: kb     = 1.38064852e-23_dp               ! [m2 kg s-2 K-1], Boltzmann constant
real(dp), parameter :: Cp     = 1012.0_dp                       ! [J kg-1 K-1], air specific heat at constant pressure,
real(dp), parameter :: p00    = 1.01325e5_dp                    ! [Pa], reference pressure at surface
real(dp), parameter :: nu_air = 1.59e-5_dp                      ! [m2 s-1], kinematic viscosity of air
real(dp), parameter :: Omega  = 2*PI/(24.0_dp*60.0_dp*60.0_dp)  ! [rad s-1], Earth angular speed
real(dp), parameter :: lambda = 300.0_dp                        ! maximum mixing length, meters
real(dp), parameter :: vonk   = 0.4_dp                          ! von Karman constant, dimensionless

real(dp), parameter :: ug = 10.0d0, vg = 0.0d0  ! [m s-1], geostrophic wind

! Latitude and longitude of Hyytiala
!real(dp), parameter :: latitude_deg  = 61.8455d0  ! [degN]
!real(dp), parameter :: longitude_deg = 24.2833d0  ! [degE]
real(dp), parameter :: latitude_deg  = 56.1d0  ! [degN]
real(dp), parameter :: longitude_deg = 13.42d0  ! [degE]

real(dp), parameter :: latitude      = latitude_deg  * PI/180.0d0  ! [rad]
real(dp), parameter :: longitude     = longitude_deg * PI/180.0d0  ! [rad]

real(dp), parameter :: fcor = 2*Omega*sin(latitude)  ! Coriolis parameter at Hyytiala

real(dp), parameter :: C_T1 = 95 * 1000    ! [J mol-1]
real(dp), parameter :: C_T2 = 230 * 1000   ! [J mol-1]
real(dp), parameter :: T_S  = 303.15 ! [K]
real(dp), parameter :: T_M  = 314    ! [K]
real(dp), parameter :: alpha = 0.0027
real(dp), parameter :: beta  = 0.09  ! [K-1]
real(dp), parameter :: C_L1  = 1.006
real(dp), parameter :: D_m  = 0.0538    ! [g/cm2]
real(dp), parameter :: ee = 100    !  [ng/g (needle-dry-weight)/h]
real(dp), parameter :: Mm_iso = 68   ![g/mol]
real(dp), parameter :: Mm_mono = 136 ![g/mol]
real(dp), parameter :: ppb = 1e-9_dp

!-----------------------------------------------------------------------------------------
! Grid parameters
!-----------------------------------------------------------------------------------------
integer, parameter :: nz = 50  ! [-], number of height levels

! Model height levels, [m]
real(dp), parameter, dimension(nz) :: &
  hh = (/    0,   10,   20,   30,   40,   50,   60,   70,   80,   90, &
           100,  120,  140,  160,  180,  200,  230,  260,  300,  350, &
           400,  450,  500,  550,  600,  650,  700,  800,  900, 1000, &
          1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, &
          2100, 2200, 2300, 2400, 2500, 2600, 2700, 2800, 2900, 3000 /)

real(dp), parameter :: hc = 10.0_dp  ! [m], canopy height

!-----------------------------------------------------------------------------------------
! Time variables
!-----------------------------------------------------------------------------------------
integer, parameter :: one_hour = 60*60  ! [s], one hour in seconds

real(dp) :: time                  ! [s], current time
real(dp) :: time_start, time_end  ! [s], start and end time

real(dp) :: dt         ! [s], time step for main loop, usually is equal to meteorology time step
real(dp) :: dt_emis    ! [s], time step for emission calculation
real(dp) :: dt_chem    ! [s], time step for chemistry calculation
real(dp) :: dt_depo    ! [s], time step for deposition calculation
real(dp) :: dt_aero    ! [s], time step for aerosol calculation
real(dp) :: dt_output  ! [s], time step for output

real(dp) :: time_start_emission    ! [s], time to start calculating emission
real(dp) :: time_start_chemistry   ! [s], time to start calculating chemistry
real(dp) :: time_start_deposition  ! [s], time to start calculating deposition
real(dp) :: time_start_aerosol     ! [s], time to start calculating aerosol

integer :: daynumber_start  ! [day], start day of year
integer :: daynumber        ! [day], current day of year

integer :: counter  ! [-], counter of time steps

!-----------------------------------------------------------------------------------------
! Meteorology variables
!-----------------------------------------------------------------------------------------
real(dp), dimension(nz  ) :: uwind, &    ! [m s-1], u component of wind
                             uwind_new, & ! [m s-1], u component of wind for the loop tmp varible
                             vwind, &    ! [m s-1], v component of wind
                             vwind_new,& ! [m s-1], v component of wind for the loop tmp varible
                             theta,  &   ! [K], potential temperature
                             theta_new   ! [K], potential temperature for the loop tmp varible
real(dp), dimension(nz  ) :: temp, &     ! [K], air temperature
                             pres        ! [Pa], air pressure
real(dp), dimension(nz  ) :: km,   &     ! for momentum 
                             kh,   &     ! for heat and other scalars
                             mix_len,  & ! mixing length
                             Ri ,    &   ! Rechardson number
                             fm,    &    ! Dyer-Businger form for momentum
                             fh,    &     ! Dyer-Businger form for h
                             wind_shear  ! wind shear
integer :: i, j  ! used for loops
!-----------------------------------------------
! emission variables
!-----------------------------------------------
real(dp)  :: F_veg_iso           ! [(molecule number) m-3 s-1] surface emission flux of Isoprene 
real(dp)  :: F_veg_mono          ! [(molecule number)  m-3 s-1] surface emission flux of Monoterpene
real(dp)  :: Factor_iso
real(dp)  :: Factor_mono 
real(dp)  :: C_L                         ! Light
real(dp)  :: C_T                        ! Temerature T (leaf)
real(dp)  :: PAR                       ! rediation
real(dp)  :: T_env                      ! [K] environment temperature 


!----------------------------------------------
! chemistry varibles
!----------------------------------------------
real(dp),dimension(nz)  :: Mair        ! Air molecules concentration [molecules/cm3]
real(dp),dimension(nz)  :: O2   ! Oxygen concentration [molecules/cm3]
real(dp),dimension(nz)  :: N2   ! Nitrogen concentration [molecules/cm3]
real(dp)  :: H2O  = 1.0D16                               ! Water molecules [molecules/cm3]
integer   :: layer    ! used for loops
real(dp),dimension(neq,nz) :: Cons, cons_new, dcondt     ! concentration matrix with neq rows and nz columns
real(dp),dimension(neq)    :: Conc 
real(dp)  :: exp_coszen             ! the radiation related quantities 
real(dp), dimension(nr_bins,nz) :: particle_conc_1d, particle_conc_1d_new
real(dp), dimension(nz) :: PM_1d, PN_1d, PV_1d
real(dp), dimension(nr_cond,nz) :: cond_sink_1d = 0.001
!-----------------------------------------------------------------------------------------
! Initialization
!-----------------------------------------------------------------------------------------

call time_init()                 ! initialize time
call meteorology_init()          ! initialize meteorology
do i=1,nz
  call Aerosol_init(diameter, particle_mass, particle_volume, particle_conc_1d(:,i), &
  particle_density, nucleation_coef, molecular_mass, molar_mass, &
  molecular_volume, molecular_dia, mass_accomm)
enddo

call open_files()        ! open output files
write(10,'(50f16.4)') hh
write(11,'(50f16.4)') hh
write(12,'(50f16.4)') hh
write(13,'(50f16.4)') hh
write(14,'(50f16.4)') hh
write(15,'(50f16.4)') hh


call write_files(time)   ! write initial values



!-----------------------------------------------------------------------------------------
! Start main loop
!-----------------------------------------------------------------------------------------
do while (time <= time_end)
  !---------------------------------------------------------------------------------------
  ! Meteorology
  !---------------------------------------------------------------------------------------
  ! Set lower boundary condition
  call surface_values(theta(1), time+dt)  ! theta = temperature at the surface   
  
  do i=1,nz-1
    wind_shear(i) = sqrt( ( ( uwind(i+1) - uwind(i) ) / ( hh(i+1) - hh(i) ) )**2 +    &
                    ( ( vwind(i+1) - vwind(i) ) / ( hh(i+1) - hh(i) ) )**2 )
    mix_len(i) = ( vonk * ( hh(i+1) + hh(i) ) / 2.0 ) /                               & 
                 ( 1.0 + ( vonk * ( hh(i+1) + hh(i) )/2.0 ) / lambda)  
  
    Ri(i) = ( grav / ( ( theta(i+1) + theta(i) ) / 2 ) ) *                              &
            ( ( theta(i+1) - theta(i) ) / ( ( uwind(i+1) - uwind(i) )**2 + ( vwind(i+1) - vwind(i) )**2 ) ) * &
            ( hh(i+1) - hh(i) )
    if (Ri(i)<0) then
      fm(i)=(1-16*Ri(i))**(0.5)
      fh(i)=(1-16*Ri(i))**(0.75)
    else if (Ri(i)>=0 .and. Ri(i)<0.2) then
      fm(i) = max((1-5*Ri(i))**2,0.1)
      fh(i) = fm(i)
    else if (Ri(i)>=0.2) then
      fm(i) = 0.1
      fh(i) = 0.1
    endif
  enddo
  km = mix_len**2 * wind_shear * fm
  kh = mix_len**2 * wind_shear * fh

  do i=2,nz-1
    uwind_new(i) = uwind(i)               +      &
               dt*fcor*(vwind(i)-vg)      +      &
               dt*2/(hh(i+1)-hh(i-1))*(km(i)*(uwind(i+1)-uwind(i))/ &
               (hh(i+1)-hh(i))-(km(i-1)*(uwind(i)-uwind(i-1))/(hh(i)-hh(i-1))))  
    vwind_new(i) = vwind(i)               -      &
               dt*fcor*(uwind(i)-ug)      +      &
               dt*2/(hh(i+1)-hh(i-1))*(km(i)*(vwind(i+1)-vwind(i))/ &
               (hh(i+1)-hh(i))-(km(i-1)*(vwind(i)-vwind(i-1))/(hh(i)-hh(i-1))))  
    
    theta_new(i) = theta(i)               +      &
               dt*2/(hh(i+1)-hh(i-1))*(kh(i)*(theta(i+1)-theta(i))/ &
               (hh(i+1)-hh(i))-(kh(i-1)*(theta(i)-theta(i-1))/(hh(i)-hh(i-1)))) 
               
   
    enddo 
  uwind(2:nz-1) = uwind_new(2:nz-1)
  vwind(2:nz-1) = vwind_new(2:nz-1)
  theta(2:nz-1) = theta_new(2:nz-1)
  
  !print*, "before define temp and pres"
  temp = theta - (grav/Cp)*hh 
  pres = barometric_law(p00, temp, hh)
  !print*, "after define temp and pres"
  !print*,temp,pres
  
 



  
  ! Update meteorology

  !---------------------------------------------------------------------------------------
  ! Emission
  !---------------------------------------------------------------------------------------
  ! Start to calculate emission after time_start_emission
  ! Compute emission part every dt_emis, multiplying 1000 to convert s to ms to make mod easier
  if ( use_emission .and. time >= time_start_emission ) then
    if ( mod( nint((time - time_start_emission)*1000.0d0), nint(dt_emis*1000.0d0) ) == 0 ) then
      ! Calculate emission rates
      !print*,"start to do emission"
      T_env = theta(2) - 0.0098 * 10
      PAR = 1000 * get_exp_coszen(time,daynumber,latitude)
      C_L = alpha * C_L1 * PAR / sqrt( 1 + alpha**2 * PAR**2)
      C_T = exp( C_T1 * ( T_env - T_S ) / (Rgas * T_env * T_S ))/        &
            ( 1 + exp( C_T2 * ( T_env - T_M ) / (Rgas * T_env * T_S )))  
      
      Factor_iso  =    NA/Mm_iso /10 /1000000
      Factor_mono =  NA/Mm_mono /10 /1000000

      F_veg_iso  = D_m * ee * (0.00001/3600) * (C_L * C_T) * 1. * Factor_iso
      F_veg_mono = D_m * ee * (0.00001/3600)* exp( beta * (T_env - T_s)) * 1. * Factor_mono
      !print*, "iso = ",F_veg_iso, "  a-pinen = ",F_veg_mono
    end if
  end if

  if ( use_emission .and. (.not. use_chemistry) ) then
    ! Add emission to the number concentrations of compounds
  end if

  !---------------------------------------------------------------------------------------
  ! Deposition
  !---------------------------------------------------------------------------------------
  ! Start to calculate gas dry deposition velocity after time_start_deposition
  ! Compute deposition part every dt_depo, multiplying 1000 to convert s to ms to make mod easier
  if ( use_deposition .and. time >= time_start_deposition ) then
    if ( mod( nint((time - time_start_deposition)*1000.0d0), nint(dt_depo*1000.0d0) ) == 0 ) then
      ! Calculate deposition velocity

      ! Remove deposited concentration at level 2 which includes canopy and soil
    end if
  end if

  !---------------------------------------------------------------------------------------
  ! Chemistry
  !---------------------------------------------------------------------------------------
  ! Start to calculate chemical reactions only after some time to save the computation time
  ! Compute chemistry part every dt_chem, multiplying 1000 to convert s to ms to make mod easier
  if ( use_chemistry .and. time >= time_start_chemistry ) then
    if ( mod( nint((time - time_start_chemistry)*1000.0d0), nint(dt_chem*1000.0d0) ) == 0 ) then
      ! Solve chemical equations for each layer except boundaries
      
      ! chemisrty , dt = 10s
      !print*,"chemistry"
      exp_coszen = get_exp_coszen(time,daynumber,latitude)
      do layer = 2, nz-1
        Mair(layer) = pres(layer) * NA / ( Rgas * temp(layer) ) * 1e-6
        O2(layer) = 0.21 * Mair(layer)
        N2(layer) = 0.78 * Mair(layer)
        !Cons(1,layer)  = 24.0d0   * Mair(layer) * ppb     ! O3 concentration
        Cons(1,layer)  = 40.0d0   * Mair(layer) * ppb     ! O3 concentration for task
        Cons(5,layer)  = 0.2d0    * Mair(layer) * ppb     ! NO2
        Cons(6,layer)  = 0.07d0   * Mair(layer) * ppb     ! NO
        !Cons(9,layer)  = 100.0d0  * Mair(layer) * ppb     ! CO
        Cons(9,layer)  = 200.0d0  * Mair(layer) * ppb     ! CO for task
        Cons(11,layer) = 1759.0d0 * Mair(layer) * ppb     ! CH4
        !Cons(20,layer) = 0.5d0    * Mair(layer) * ppb     ! SO2
        Cons(20,layer) = 2.0d0    * Mair(layer) * ppb     ! SO2 for task
        !Cons(13,layer) = 2.2d0    * Mair(layer) * ppb     ! C5H8   ini  = 0
        !Cons(23,layer) = 2.2d0    * Mair(layer) * ppb     ! alpha-p ini=0
      

        if (layer == 2) then
          !print*, "Iso = ", F_veg_iso, "a-pinene = ", F_veg_mono
          call chemistry_step(Cons(1:neq,layer),time,time+dt_chem,&
              O2(layer), N2(layer), Mair(layer), H2O, temp(layer), exp_coszen, &
              F_veg_iso, F_veg_mono, Cond_sink_1d(:,layer) )
        else 
          call chemistry_step(Cons(1:neq,layer),time,time+dt_chem,&
                             O2(layer), N2(layer), Mair(layer), H2O, temp(layer), exp_coszen, &
                             0.0d0,0.0d0, Cond_sink_1d(:,layer))
        endif 
      enddo
      
    
    end if  ! every dt_chem
  end if

  ! Update concentrations of gas phase compounds if any of these processes are considered
  ! Deposition should not be used alone because it calculates nothing in that case
  if (use_emission .or. use_chemistry) then
    ! Trick to make bottom flux zero
    cons(1:neq,1) = cons(1:neq,2)
   ! cons_new(1:neq,50) = 0.

    ! Concentrations can not be lower than 0
    do i=1,neq
      do j=i,nz
        if(cons(i,j)<0.0d0) then
          cons(i,j) = 0.0d0
        endif
      enddo
    enddo

    !cons = max(cons,0.0d0)
    ! Mixing of chemical species, where i is layer
    do i = 2, nz-1
      cons_new(1:neq,i) = cons(1:neq, i)   +    &
                  dt*2/(hh(i+1)-hh(i-1))*(kh(i)*(cons(1:neq, i+1)-cons(1:neq, i))/ &
                  (hh(i+1)-hh(i))-(kh(i-1)*(cons(1:neq, i)-cons(1:neq, i-1))/(hh(i)-hh(i-1))))
    enddo
    cons(1:neq,1:nz-1) = cons_new(1:neq,1:nz-1)
    cons(1:neq,1) = cons(1:neq,2)
    ! Set the constraints above again for output
  end if

  !---------------------------------------------------------------------------------------
  ! Aerosol
  !---------------------------------------------------------------------------------------
  ! Start to calculate aerosol processes only after some time to save the computation time
  ! Compute aerosol part every dt_aero, multiplying 1000 to convert s to ms to make mod easier
  if ( use_aerosol .and. time >= time_start_aerosol ) then
    if ( mod( nint((time - time_start_aerosol)*1000.0d0), nint(dt_aero*1000.0d0) ) == 0 ) then
      ! Nucleation, condensation, coagulation and deposition of particles
      ! Nucleation loop
      do i=2,nz-1
        cond_vapour(1) = cons(21,i) * 1.D6   ! H2SO4
        cond_vapour(2) = cons(25,i) * 1.D6   ! ELVOC
        call Nucleation(particle_conc_1d(:,i), cond_vapour, dt_aero)  

        call Coagulation(dt_aero, particle_conc_1d(:,i), diameter, temp(i),pres(i),particle_mass)
        
        call Condensation(dt_aero, temp(i), pres(i), mass_accomm, molecular_mass, &
        molecular_volume, molar_mass, molecular_dia, particle_mass, particle_volume, &
        particle_conc_1d(:,i), diameter, cond_vapour, Cond_sink_1d(:,i))
      enddo
    end if


    ! Trick to make bottom flux zero
    particle_conc_1d(1:neq,1) = particle_conc_1d(1:neq,2)
    ! Concentrations can not be lower than 0 [molec m-3]
    do i=1,nr_bins
      do j=i,nz
        if(particle_conc_1d(i,j)<0.0d0) then
          particle_conc_1d = 0.0d0
        endif
      enddo
    enddo

    ! Mixing of aerosol particles
     do i = 2, nz-1
       particle_conc_1d_new(:,i) = particle_conc_1d(:,i)   +    &
                  dt*2/(hh(i+1)-hh(i-1))*(kh(i)*(particle_conc_1d(:, i+1)-particle_conc_1d(:, i))/ &
                  (hh(i+1)-hh(i))-(kh(i-1)*(particle_conc_1d(:,i)-particle_conc_1d(:,i-1))/(hh(i)-hh(i-1))))
     enddo
     particle_conc_1d(:,1:nz-1) = particle_conc_1d_new(:,1:nz-1)
     particle_conc_1d(:,1) = particle_conc_1d(:,2)

    ! Set the constraints above again for output
     do i=1,nr_bins
      do j=i,nz
        if(particle_conc_1d(i,j)<0.0d0) then
          particle_conc_1d = 0.0d0
        endif
      enddo
    enddo
    ! Update related values, e.g., total number concentration, total mass concentration
      
     do i=1,nz
       PN_1d(i) = SUM(particle_conc_1d(:,i))*1.D-6                 
       PM_1d(i) = SUM(particle_conc_1d(:,i)*particle_mass)*1.D9 
       PV_1d(i) = sum(particle_conc_1d(:,i)*particle_volume) * 1.D12
     enddo
  end if

  !---------------------------------------------------------------------------------------
  ! Ending loop actions
  !---------------------------------------------------------------------------------------
  ! Advance to next time step
  time = time + dt

  ! Write data every dt_output [s]
  if ( mod( nint((time - time_start)*1000.0d0), nint(dt_output*1000.0d0) ) == 0 ) then
    write(*, '(a8,f8.3,a8)') 'time = ', time/one_hour, '   hours'
    call write_files(time)
  end if

  ! Count loop number
  counter = counter + 1

end do

!-----------------------------------------------------------------------------------------
! Finalization
!-----------------------------------------------------------------------------------------
! Close all the opened files
call close_files()

! Count total time steps
write(*,*) counter,'time steps'


contains


!-----------------------------------------------------------------------------------------
! subroutine open_files()
!
! Open needed files
!-----------------------------------------------------------------------------------------
subroutine open_files()
  logical :: dir_exist

  ! Create a new directory if it does not exist
  inquire(file=trim(adjustl(output_dir)), exist=dir_exist)
  if (.not. dir_exist) then
    ! This line may change for different operating systems
    call system('mkdir ' // trim(adjustl(output_dir)))
  end if

  ! Open files to write output results
  open( 8,file=trim(adjustl(output_dir))//'/time.dat' ,status='replace',action='write')
  open( 9,file=trim(adjustl(output_dir))//'/hh.dat'   ,status='replace',action='write')
  open(10,file=trim(adjustl(output_dir))//'/uwind.dat',status='replace',action='write')
  open(11,file=trim(adjustl(output_dir))//'/vwind.dat',status='replace',action='write')
  open(12,file=trim(adjustl(output_dir))//'/theta.dat',status='replace',action='write')
  open(13,file=trim(adjustl(output_dir))//'/km.dat',status='replace',action='write')
  open(14,file=trim(adjustl(output_dir))//'/kh.dat',status='replace',action='write')
  open(15,file=trim(adjustl(output_dir))//'/Ri.dat',status='replace',action='write')
  open(16,file=trim(adjustl(output_dir))//'/F_veg_iso+mono.dat',status='replace',action='write')
  open(17,file=trim(adjustl(output_dir))//'/a-pinene',status='replace',action='write')
  open(18,file=trim(adjustl(output_dir))//'/H2SO4',status='replace',action='write')
  open(19,file=trim(adjustl(output_dir))//'/CO',status='replace',action='write')
  open(20,file=trim(adjustl(output_dir))//'/Isoprene',status='replace',action='write')
  open(21,file=trim(adjustl(output_dir))//'/OH',status='replace',action='write')
  open(22,file=trim(adjustl(output_dir))//'/HO2',status='replace',action='write')
  open(23,file=trim(adjustl(output_dir))//'/ELVOC',status='replace',action='write')
  open(24,file=trim(adjustl(output_dir))//'/PM_1d',status='replace',action='write')
  open(25,file=trim(adjustl(output_dir))//'/PN_1d',status='replace',action='write')
  open(26,file=trim(adjustl(output_dir))//'/PV_1d',status='replace',action='write')

  open(50,file=trim(adjustl(output_dir))//'/diameter',status='replace',action='write')
end subroutine open_files


!-----------------------------------------------------------------------------------------
! subroutine write_files(time)
!
! Write data to files at time
!-----------------------------------------------------------------------------------------
subroutine write_files(time)
  real(dp) :: time  ! current time
  character(255) :: outfmt_one_scalar, outfmt_two_scalar, outfmt_level, outfmt_mid_level

  ! Output real data with scientific notation with 16 decimal digits
  outfmt_one_scalar = '(es25.16e3)'                               ! for scalar
  write(outfmt_level     , '(a, i3, a)') '(', nz  , 'es25.16e3)'  ! for original levels
  write(outfmt_mid_level , '(a, i3, a)') '(', nz-1, 'es25.16e3)'  ! for middle levels
  write(outfmt_two_scalar, '(a, i3, a)') '(', 2   , 'es25.16e3)'  ! for two scalars

  ! Only output hh once
  if (time == time_start) then
    write(9, outfmt_level) hh
    write(50,outfmt_level) diameter
  end if

  ! Output every output time step
  write( 8, outfmt_one_scalar) time/(24*one_hour)  ! [day]
  write(10, outfmt_level     ) uwind
  write(11, outfmt_level     ) vwind
  write(12, outfmt_level     ) theta
  write(13, outfmt_level     ) km
  write(14, outfmt_level     ) kh
  write(15, outfmt_level     ) Ri
  write(16, outfmt_level     ) F_veg_iso,F_veg_mono
  write(17, outfmt_level     ) Cons(23,1:nz)
  write(18, outfmt_level     ) Cons(21,1:nz)
  write(19, outfmt_level     ) Cons(9,1:nz)
  write(20, outfmt_level     ) Cons(13,1:nz)
  write(21, outfmt_level     ) Cons(3,1:nz)
  write(22, outfmt_level     ) Cons(8,1:nz)
  write(23, outfmt_level     ) Cons(25,1:nz)
  write(24, outfmt_level     ) PM_1d
  write(25, outfmt_level     ) PN_1d
  write(26, outfmt_level     ) PV_1d
  
  

  
end subroutine write_files


!-----------------------------------------------------------------------------------------
! subroutine Close_Files()
!
! Close files
!-----------------------------------------------------------------------------------------
subroutine close_files()
  close(8)
  close(9)
  close(10)
  close(11)
  close(12)
  close(13)
  close(14)
  close(15)
  close(16)
  close(17)
  close(18)
  close(19)
  close(20)
  close(21)
  close(22)
  close(23)
  close(24)
  close(25)
  close(26)

  close(50)

end subroutine close_files


!-----------------------------------------------------------------------------------------
! subroutine time_init()
!
! Time initiation
!-----------------------------------------------------------------------------------------
subroutine time_init()
  ! Basic time variables
  time_start = 0.0d0
  time_end   = 5.0d0 * 24.0d0 * one_hour
  time       = time_start

  ! Time steps
  dt        = 0.5d0
  dt_emis   = 0.5d0
  dt_chem   = 10.0d0
  dt_depo   = 10.0d0
  dt_aero   = 10.0d0
  dt_output = 3600.0d0

  ! Day number
  !daynumber_start = 31+28+31+30+31+30+31+10  ! day is Aug. 10
  daynumber_start = 31+28+31+6 ! day is Aug. 10
  daynumber       = daynumber_start

  ! Start time for each process
  !time_start_emission   = 3*24*one_hour
  time_start_emission   = 4.625D0*24*one_hour

  time_start_chemistry  = 3*24*one_hour
  time_start_deposition = 3*24*one_hour
  !time_start_aerosol    = 3*24*one_hour
  time_start_aerosol    = 4.0D0*24*one_hour


  ! Loop number
  counter = 0
end subroutine time_init


!-----------------------------------------------------------------------------------------
! subroutine meteorology_init()
!
! Meteorology initiation
!-----------------------------------------------------------------------------------------
subroutine meteorology_init()
  ! Wind velocity
  uwind         = 0.0d0
  uwind(nz)     = ug
  uwind(2:nz-1) = uwind(nz) * hh(2:nz-1)/hh(nz)

  vwind = 0.0d0
  vwind(nz) = vg

  ! Potential temperature
  !theta     = 273.15d0 + 25.0d0
  !theta(nz) = 273.15d0 + 30.0d0
  theta     = 273.15d0 + 0.0d0
  theta(nz) = 273.15d0 + 5.0d0

  ! Air temperature and pressure
  temp = theta - (grav/Cp)*hh 
  pres = barometric_law(p00, temp, hh)

  ! km and kh 
   km = 5
   kh = 5

  ! initial the concentrations of 1-D model
   Cons(1:neq,1:nz)  = 0.0d0


end subroutine meteorology_init


!-----------------------------------------------------------------------------------------
! Get the surface values from the input data file
! Now only the temperature is used.
!-----------------------------------------------------------------------------------------
subroutine surface_values(temperature, time)

  ! (Note: can also get water concentrantion, in ppt, if modify this
  ! subroutine to also use column 8)
  !
  ! Data is taken from:
  ! http://avaa.tdata.fi/web/smart

  real(dp), intent(in)            :: time ! input, in seconds
  real(dp), intent(out)           :: temperature ! output, in Kelvin
  logical, save                   :: first_time = .true.
  real(dp), dimension(8,50), save :: surface_data
  real(dp), dimension(50), save   :: temperature_data
  real(dp), parameter             :: seconds_in_day = 24*60*60
  real(dp), parameter             :: seconds_in_30min = 30*60
  integer                         :: index
  real(dp) :: time24h, time30min, time24plus15, temp1, temp2, x

  ! Only when called for the first time, read in data from file
  ! With this trick, we don't need to open the file in the main program
  IF (first_time) THEN
     !open(30, file=trim(adjustl(input_dir))//'/hyytiala_2011_8_10_t_h2o.dat', status='old')
     open(30, file=trim(adjustl(input_dir))//'/hyltemossa_2018_4_06_t_h2o.dat', status='old')
     read(30, *) surface_data
     temperature_data(1:50) = surface_data(7,1:50) ! in Celcius
     first_time = .false.
  end IF

  time24h = modulo(time, seconds_in_day) ! time modulo 24 hours
  time24plus15 = time24h + 15*60 ! time since 23:45 previous day
  time30min = modulo(time24plus15, seconds_in_30min)
  index = 1 + floor(time24plus15/seconds_in_30min)

  temp1 = temperature_data(index)
  temp2 = temperature_data(index + 1)
  x = time30min/seconds_in_30min
  
  ! linear interpolation between previous and next temperature data value
  temperature = temp1 + x*(temp2 - temp1) + 273.15_dp  ! now in Kelvin
end subroutine surface_values


!-----------------------------------------------------------------------------------------
! Calculate the radiation related quantities
!-----------------------------------------------------------------------------------------
real(dp) function get_exp_coszen(time,daynumber,latitude)
  real(dp), intent(in) :: time,latitude
  INTEGER, intent(in) :: daynumber
  real(dp) :: hourangle,zenith,coszen
  hourangle = get_hourangle(time)
  zenith = solar_zenith_angle(hourangle,daynumber,latitude)
  coszen = cos(zenith)
  IF (coszen > 0) THEN  ! sun is above horizon
     get_exp_coszen = exp(-0.575_dp/coszen)
  ELSE
     get_exp_coszen = 0.0_dp
  endIF
end function get_exp_coszen


real(dp) function get_hourangle(time)
  real(dp), intent(in) :: time
  real(dp), parameter :: one_day = 24*one_hour
  get_hourangle = modulo(time,one_day)/one_day * 2 * pi - pi
end function get_hourangle


real(dp) function solar_zenith_angle(hourangle,daynumber,latitude)
  ! http://en.wikipedia.org/wiki/Solar_elevation_angle
  ! http://en.wikipedia.org/wiki/Position_of_the_Sun
  INTEGER, intent(in) :: daynumber
  real(dp), intent(in) :: hourangle,latitude
  real(dp) :: declination,elevation
  real(dp), parameter :: to_rad = pi/180.0_dp

  declination = -23.44_dp * to_rad * cos(2 * pi * (daynumber + 10)/365.0_dp)
  elevation = cos(hourangle)*cos(declination)*cos(latitude) &
       + sin(declination)*sin(latitude)
  solar_zenith_angle = pi/2.0_dp - elevation
  ! Notes:
  ! - Not tested near equador or on the southern hemisphere.
  ! - solar_zenith_angle can be larger than pi/2, it just means
  !   the sun is below horizon.
  ! - solar_zenith_angle assumes time is in local solar time, which
  !   is usually not exactly true
end function solar_zenith_angle


!-----------------------------------------------------------------------------------------
! Other functions
!-----------------------------------------------------------------------------------------
function barometric_law(p00, tempK, h) result(p)
  real(dp), intent(in) :: p00, tempK(nz), h(nz)
  real(dp) :: p(nz)
  real(dp) :: dh(nz)
  dh(2:nz) = h(2:nz) - h(1:nz-1)

  p(1) = p00
  do i=2, nz
    p(i) = p(i-1)*exp(-mm_air*grav/(Rgas*(tempK(i-1)+tempK(i))/2.0d0)*dh(i))
  end do
end function barometric_law

end program main


