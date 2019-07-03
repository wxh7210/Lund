MODULE aerosol_mod

IMPLICIT NONE

PRIVATE

PUBLIC :: nr_bins
PUBLIC :: DP, cond_vapour, diameter, particle_mass, particle_volume, particle_conc, &
  particle_density, nucleation_coef, molecular_mass, molar_mass, &
  molecular_volume, molecular_dia, mass_accomm, &
  PN, PM  , PV , cond_sink, nr_cond
PUBLIC :: Aerosol_init, Nucleation, Condensation, Coagulation, dry_dep_velocity

!! ====================== Definition of variables =====================================================================
INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15,300)
! so that numbers will be in 64bit floating point
! http://en.wikipedia.org/wiki/Double_precision_floating-point_format

REAL(DP), PARAMETER :: pi = 2D0*ASIN(1D0)
REAL(DP), PARAMETER :: ka = 0.4D0              ! von Karman constant, dimensionless
REAL(DP), PARAMETER :: g = 9.81D0              ! gravitation const. m s^-2
INTEGER,  PARAMETER ::  nr_bins = 100          ! Number of particle size bins
INTEGER, PARAMETER ::  nr_cond = 2             ! Number of condensable vapours
REAL(DP), PARAMETER :: Rg=8.3145D0             ! Universal gas constant J mol^-1 K^-1
REAL(DP), PARAMETER :: Na=6.022D23             ! Avogadro's number 
REAL(DP), PARAMETER :: Mair=28.96D-3           ! Mean molecular weight of air
REAL(DP), PARAMETER :: kb = 1.381d-23          ! Boltzmann constant [(m2*kg)/(s2*K)]
  
REAL(DP), DIMENSION(nr_bins) :: diameter  ,&    ! Diameter of each size bin
     particle_mass                        ,&    ! mass of one particle in each size bin
     particle_conc                        ,&    ! number concentration in each size bin
     particle_volume                      ,&    ! volume concentration in each size bin 
     coag_loss                            ,&    ! coagulation loss rate of particles in each size bin
     v_dep                                      ! Dry deposition velocity of particles
     
REAL(DP), DIMENSION(nr_cond) :: molecular_mass   ,&   ! molecular mass of the condensing vapours [kg/#]
     molecular_volume                            ,&   ! Molecule volume of condensable vapours [m^3]
     molecular_dia                               ,&   ! Molecule diameter of condensable vapours [m]
     molar_mass                                  ,&   ! Molar mass of condensable vapours [kg/m^3]
     cond_vapour                                      ! Concentration of condensable vapours [molec/m^3] 
                                      
REAL(DP), DIMENSION(nr_cond) :: Cond_sink = 0.001     ! Assumed initial condensation sink of vapours [s^-1]

REAL(DP) ::         PN, PM, PV               ! Total particle number and mass concentration [cm^-3] and Volume
     
REAL(DP) :: vd_SO2,vd_O3, vd_HNO3           ! Dry deposition velocity of SO2, O3 & HNO3     

REAL(DP)    particle_density              ,&    ! [Kg]
            nucleation_coef               ,&    ! Nucleation coefficient 
            mass_accomm                         ! mass accomodation coefficient 

REAL(dp) :: nucleation_rate                     ! #/(m3*s)
 

CONTAINS

SUBROUTINE Aerosol_init(diameter, particle_mass, particle_volume, particle_conc, &
  particle_density, nucleation_coef, molecular_mass, molar_mass, &
  molecular_volume, molecular_dia, mass_accomm)

    !! ====================== Definition of variables =====================================================================

    REAL(DP), DIMENSION(nr_bins), INTENT(OUT) :: diameter       , &    ! diamter of each size bin
                                                 particle_mass  , &    ! mass of one particle
                                                 particle_volume  , &  ! volume of one particle
                                                 particle_conc         ! number concentration

    REAL(DP), DIMENSION(nr_cond), INTENT(OUT) :: molecular_mass ,& ! molecular mass of the condensing vapours [kg/#]
                                           molecular_volume ,&     ! [m3]
                                           molecular_dia, &        ! [m]
                                           molar_mass              ! molar mass of the condensing vapours [kg/mol]

    REAL(DP), INTENT(OUT) :: nucleation_coef, mass_accomm

    REAL(DP), DIMENSION(nr_cond) :: density             ! Bulk density of condensing vapours [kg/m^3]
 
    REAL(DP)  ::  particle_density
    
    INTEGER :: i
    mass_accomm = 1D0   ! Mass accommodation coefficient 
    
    !nucleation_coef = 1D-20 
    nucleation_coef = 1D-21   ! for the pontus' task
    
    ! Particle diameters between 2D-9 and 2.5D-6 m:
    diameter(1)=2D-9 
    DO i=2,nr_bins
      diameter(i)=diameter(i-1)*(2.5D-6/diameter(1))**(1D0/(nr_bins-1))
    END DO
      
    particle_conc = 1D0 ! Assume an initial particle number concentration of 1 m^-3
    where((abs(diameter-2D-7)-MINVAL(abs(diameter-2D-7)))<1D-20)  particle_conc=2D8 ! add 200 cm^-3 200 nm sized accumulation mode particles
    
    particle_density = 1.4D3                                        ! Assumed fixed particle density [kg/m^3]
    particle_volume = 1D0/6D0 * pi * diameter**3                      ! Single particle volume (m^3)
    particle_mass=  1D0/6D0 * pi * diameter**3 * particle_density     ! [kg]

    density = (/1.84D3, 1.4D3/)                                     ! density of sulphuric acid and SOA
    molar_mass = (/0.098D0, 0.3D0/)                                 ! H2SO4 and ELVOC
    molecular_mass = molar_mass / Na                                ! molecular mass [kg]
    molecular_volume = molecular_mass / density                     ! molecular volume [m^3]
    molecular_dia = (6D0 * molecular_volume / pi )**(1D0/3D0)       ! molecular diameter [m]

END SUBROUTINE Aerosol_init
  



!----------------------------
SUBROUTINE Nucleation(particle_conc, cond_vapour, timestep)                            ! (Add input and output variables here)
  
  ! Consider how kinetic H2SO4 nucleation influence the number concentrations of particles 
  ! in the fist size bin particle_conc(1) within one model time step
  REAL(DP), dimension(nr_bins),INTENT(INOUT) :: particle_conc
  REAL(DP)             :: particle_conc_1_new
  real(dp), intent(in) :: timestep
  real(dp), dimension(2),intent(in) :: cond_vapour
  particle_conc_1_new = particle_conc(1)                           +   &
                        nucleation_coef * cond_vapour(1)**2 * timestep 
  particle_conc(1) = particle_conc_1_new

END SUBROUTINE Nucleation



!----------------------------------
!----------------------------------
SUBROUTINE Condensation(timestep, temperature, pressure, mass_accomm, molecular_mass, &
  molecular_volume, molar_mass, molecular_dia, particle_mass, particle_volume, &
  particle_conc, diameter, cond_vapour,Cond_sink) ! Add more variables if you need it
  
    REAL(DP), DIMENSION(nr_bins), INTENT(IN) :: diameter, particle_mass
    REAL(DP), DIMENSION(nr_cond), INTENT(IN) :: molecular_mass, molecular_dia, &
    molecular_volume, molar_mass
    REAL(DP), INTENT(IN) :: timestep, temperature, pressure, mass_accomm

    REAL(DP), DIMENSION(nr_bins), INTENT(INOUT) :: particle_conc

    REAL(DP), DIMENSION(2), INTENT(IN) :: cond_vapour  ! condensing vapour concentrations, which is H2SO4 and organics (ELVOC) [#/m^3]
    
    REAL(DP), DIMENSION(nr_bins), INTENT(IN)   :: particle_volume

    real(dp), dimension(nr_cond),intent(out) :: Cond_sink  ! condensation sink of H2SO4 and ELVOC

    REAL(DP), DIMENSION(nr_bins)   ::  slip_correction, diffusivity, speed_p, &
                                       particle_conc_new , particle_volume_new
    
    REAL(DP), DIMENSION(nr_cond)   ::  diffusivity_gas, speed_gas
    
    REAL(DP) :: dyn_visc, l_gas, dens_air , x1 , x2
    
    INTEGER :: i,j

    real(dp),dimension(nr_cond,nr_bins) :: mean_free_path,  &     ! mean free path of gas molecule
                                           kn            ,  &     ! Knudsen
                                           F_S_correction,  &     ! fuchs-sutugin correction factor
                                           CR                     ! collision rate between a prticle and a gas
    ! Add more variabels as you need it...

    dyn_visc = 1.8D-5*(temperature/298D0)**0.85D0  ! dynamic viscosity of air
    dens_air=Mair*pressure/(Rg*temperature)        ! Air density
    l_gas=2D0*dyn_visc/(pressure*SQRT(8D0*Mair/(pi*Rg*temperature))) ! Gas mean free path in air (m)

    slip_correction = 1D0+(2D0*l_gas/(diameter))*&
    (1.257D0+0.4D0*exp(-1.1D0/(2D0*l_gas/diameter))) ! Cunninghams slip correction factor (Seinfeld and Pandis eq 9.34) 
    
    diffusivity = slip_correction*kb*temperature/(3D0*pi*dyn_visc*diameter)   ! Diffusivity for the different particle sizes m^2/s
    speed_p = SQRT(8D0*kb*temperature/(pi*particle_mass))                     ! speed of particles (m/s)
  
    diffusivity_gas=5D0/(16D0*Na*molecular_dia**2D0*dens_air)*&
    SQRT(Rg*temperature*Mair/(2D0*pi)*(molar_mass+Mair)/molar_mass)            ! Diffusivity of condensable vapours (m^2 s^-1)

    ! Thermal velocity of vapour molecule
    speed_gas=SQRT(8D0*kb*temperature/(pi*molecular_mass)) ! speed of H2SO4 molecule
    
    ! Calculate the Fuchs-Sutugin correction factor:  (beta)
    do i=1,nr_cond
      do j=1,nr_bins
        mean_free_path(i,j) = 3 * ( diffusivity_gas(i) + diffusivity(j) ) /  &
                               ( (speed_gas(i)**2 + speed_p(j)**2)**0.5   ) 
        kn(i,j) = 2 * mean_free_path(i,j) / ( diameter(j)+molecular_dia(i) )
        F_S_correction(i,j) = ( 0.75 * ( 1 + kn(i,j) ) )  /   &
                         ( kn(i,j)**2 + kn(i,j) + 0.283 * kn(i,j) + 0.75 )
        CR(i,j) = 2 * pi * (diffusivity_gas(i) + diffusivity(j)) * (molecular_dia(i)+ diameter(j) ) &
                    * F_S_correction(i,j)
      enddo
    enddo
    ! Calculate the Collision rate (CR [m^3/2]) between gas molecules (H2SO4 and ELVOC) and the particles:
    
    ! Calculate the new single particle volume after condensation (particle_volume_new):
       !! wang: loop for all particle size
      do i=1, nr_bins
        particle_volume_new(i) =  particle_volume(i)            +  &
                                  CR(1,i) * cond_vapour(1) * molecular_volume(1) * timestep + &
                                  CR(2,i) * cond_vapour(2) * molecular_volume(2) * timestep   ! molucular_volume will change ?   
      enddo
      !particle_conc_new = particle_volume * particle_density

    ! Use the full-stationary method to divide the particles between the existing size bins (fixed diameter grid):
      particle_conc_new=0D0 ! Initialise a new vector with the new particle concentrations
      particle_conc_new(nr_bins)=particle_conc(nr_bins)  

    DO j = 1,nr_bins-1
    ! Add equations that redistributes the particle number concentration 
    !in size bin 1 to nr_bins-1 to the fixed volume (diameter) grid    
      x1 = (particle_volume(j+1)- particle_volume_new(j)) / (particle_volume(j+1) - particle_volume(j))
      x2 = 1D0 - x1
      particle_conc_new(j)   = particle_conc_new(j) + x1 * particle_conc(j)
      particle_conc_new(j+1) = particle_conc_new(j+1) + x2 * particle_conc(j) 
    END DO
    ! Update the particle concentration in the particle_conc vector:
    particle_conc = particle_conc_new   ! already done

    ! caculate the condensation sink of H2SO4 and ELVOC 
    do i=1,nr_cond
      Cond_sink(i) = sum(particle_conc * CR(i,:))   
    enddo
    
END SUBROUTINE Condensation


!-----------------------------------
!-----------------------------------
SUBROUTINE Coagulation(timestep, particle_conc, diameter, &
  temperature,pressure,particle_mass) ! Add more variables if you need it
  
    REAL(DP), DIMENSION(nr_bins), INTENT(IN) :: diameter
    REAL(DP), DIMENSION(nr_bins), INTENT(INOUT) :: particle_conc
    REAL(DP), INTENT(IN) :: timestep
    REAL(DP), DIMENSION(nr_bins), INTENT(IN) :: particle_mass       ! mass of one particle                                 
    REAL(DP), INTENT(IN) :: temperature, pressure
    
    REAL(DP), DIMENSION(nr_bins,nr_bins) :: coagulation_coef        ! coagulation coefficients [m^3/s]

    REAL(DP), DIMENSION(nr_bins) :: slip_correction, diffusivity, dist, speed_p, &
    Beta_Fuchs, free_path_p , loss1,  loss2,   loss

    REAL(DP) ::       dyn_visc, &                                   ! dynamic viscosity, kg/(m*s)
                      l_gas                                         ! Gas mean free path in air
    INTEGER  :: i,j

! The Coagulation coefficient is calculated according to formula 13.56 in Seinfield and Pandis (2006), Page 603
dyn_visc = 1.8D-5*(temperature/298.0d0)**0.85                                              ! Dynamic viscosity of air

l_gas=2D0*dyn_visc/(pressure*SQRT(8D0*Mair/(pi*Rg*temperature)))                        ! Gas mean free path in air (m)

slip_correction = 1D0+(2D0*l_gas/(diameter))*&
(1.257D0+0.4D0*exp(-1.1D0/(2D0*l_gas/diameter)))                                        ! Cunninghams slip correction factor (Seinfeld and Pandis eq 9.34)

diffusivity = slip_correction*kb*temperature/(3D0*pi*dyn_visc*diameter)                 ! Diffusivity for the different particle sizes m^2/s

speed_p = SQRT(8D0*kb*temperature/(pi*particle_mass))                                   ! Speed of particles (m/s)

free_path_p = 8D0*diffusivity/(pi*speed_p)                                              ! Particle mean free path (m)

dist = (1D0/(3D0*diameter*free_path_p))*((diameter+free_path_p)**3D0 &
-(diameter**2D0+free_path_p**2D0)**(3D0/2D0))-diameter                    ! mean distance from the center of a sphere reached by particles leaving the sphere's surface (m)

DO i = 1,nr_bins
   Beta_Fuchs = 1D0/((diameter+diameter(i))/(diameter+diameter(i)+&
   2D0*(dist**2D0+dist(i)**2D0)**0.5D0)+8D0*(diffusivity+diffusivity(i))/&
   (((speed_p**2D0+speed_p(i)**2D0)**0.5D0)*(diameter+diameter(i))))                    ! Fuchs correction factor from Seinfeld and Pandis, 2006, p. 600

   coagulation_coef(i,:) = 2D0*pi*Beta_Fuchs*(diameter*diffusivity(i)+&
   diameter*diffusivity+diameter(i)*diffusivity+diameter(i)*diffusivity(i))             ! coagulation rates between two particles of all size combinations  (m^3/s)    
END DO

! Write equations that considers how the particle number concentration in each size bin 
!(particle_conc) is influenced by the coagulation sink (loss of smaller particles when
! they collide with larger ones)
   ! double loop + if(i==j)  
loss2 = 0.
do i=1,nr_bins
  do j=i,nr_bins
    if (i==j) then
      loss1(i) = -0.5 * coagulation_coef(i,j) * particle_conc(i)**2 * timestep
    else
      loss2(i) = loss2(i) - particle_conc(i) * coagulation_coef(i,j)*particle_conc(j) * timestep
    endif  
  enddo
enddo
loss = loss1 + loss2
particle_conc = particle_conc + loss
! You can first calculate the loss (loss1) do to self-coagulation between particles in the same size bin
! and then calculate the loss (loss2) due to coagulation with larger particles
! Then add the two loss terms together loss = loss1 + loss2 
 




END SUBROUTINE Coagulation
  
SUBROUTINE dry_dep_velocity(diameter,particle_density,temperature,pressure,DSWF, & 
   Richards_nr10m,wind_speed10m) ! Add more variables if you need it
   
      REAL(dp), DIMENSION(nr_bins), INTENT(IN) :: diameter
      
      REAL(dp), INTENT(IN) :: temperature, pressure, Richards_nr10m, DSWF, &
      wind_speed10m, particle_density
            
      REAL(dp) :: z0m, r_coll, a_landuse, j_landuse, v_kinematic,dyn_visc,l_gas,Pr,beta,&
      gam,zr,u_friction,dens_air, L_Ob, raO3, raSO2, raHNO3, raisoprene, raapinene
      
      ! Specific parameters for the surface resistance of gases:
      REAL(dp) :: rj,rlu,rac, &
       DiffusivityH2O, D_ratio_SO2, D_ratio_O3, D_ratio_HNO3, D_ratio_isoprene, D_ratio_apinene, &
       DiffusivitySO2, DiffusivityO3, DiffusivityHNO3, Diffusivityisoprene, Diffusivityapinene,&
       z_roughSO2, z_roughO3, z_roughHNO3, z_roughisoprene, z_roughapinene, &
       ScSO2, ScO3, ScHNO3, Scisoprene, Scapinene, &
       rbSO2, rbO3, rbHNO3, rbisoprene, rbapinene, &
       H_effSO2, H_effO3, H_effHNO3, H_effisoprene, H_effapinene, &
       f0_SO2, f0_O3, f0_HNO3, f0_isoprene, f0_apinene, &
       rclSO2, rclO3, rgsSO2, rgsO3
            
       dens_air = Mair*pressure/(Rg*temperature)    ! Air density (kg/m^3)
       dyn_visc = 1.8D-5*(temperature/298.)**0.85   ! dynamic viscosity of air (kg/(m*s))
       v_kinematic = dyn_visc/dens_air              ! kinematic viscosity of air (m^2/s)
       
       zr=10D0                                      ! Reference height [m]
       L_Ob=zr/Richards_nr10m                       ! Monin-Obukhov length scale
       z0m = 0.9D0 ! Surface roughness length for momentum evergreen, needleleaf trees (m)     
       u_friction=ka*wind_speed10m/(log(zr/z0m))     ! Friction velocity (Eq. 16.67 from Seinfeld and Pandis, 2006)
 
       ! Land use category paramaters from Seinfeld and Pandis, 2006 Table 19.2: 
       r_coll = 2D-3 ! radius of collector evergreen, needleleaf trees
       
       ! coefficients based on land use categories (evergreen, needleleaf trees)
       a_landuse = 1D0
       j_landuse = 0.56D0
    
       Pr = 0.95D0        ! Turbulent Prandtl number (when ka = 0.4 (Hogstrom, 1988))
       beta = 7.8D0       ! When ka = 0.4 (Hogstrom, 1988)
       gam = 11.6D0       ! When ka = 0.4 (Hogstrom, 1988)

       ! Calculate the particle sedimentation velocity:
               
       ! Calculation of aerodynamic resistance for particles for:
       ! stable boundary layer (Ri>1D-6)
       ! neutral boundary layer (abs(Ri)<1D-6
       ! unstable boundary layer Ri<-1D-6
       
        
       ! Calculate the quasi-laminar resistance (rb) for particles:
       
       ! Calculate the dry deposition velocity for particles:
       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! Calculate the dry deposition velocity for O3, SO2, HNO3, isoprene and a-pinene: 
       
       ! Resistance components used when calculating the surface resistance for gases, 
       ! table 19.3 Seinfeld and Pandis, 2006: 
      
      ! The minimum, bulk canopy stomatal resistance for water vapor:
      rj = 130D0 ! (s/m) Summer, evergreen, needleleaf
      
      ! The resistance of the outer surfaces in the upper canopy
      rlu = 2000D0 ! (s/m) Summer, evergreen, needleleaf
      
      ! transfer resistance on the ground (that depends only on canopy height)
      rac = 2000D0 ! (s/m) Summer, evergreen, needleleaf
      
      ! resistance for uptake by soil, leaf litter, and so on at the ground SO2
      rgsSO2 = 500D0 ! (s/m) Summer, evergreen, needleleaf
      
      ! restistance for uptake by soil, leaf litter, and so on at the ground, O3
      rgsO3 = 200D0 ! (s/m) Summer, evergreen, needleleaf
      
      ! resistance for uptake by leaves,twigs, and other exposed surfaces, SO2
      rclSO2 = 2000D0 ! (s/m) Summer, evergreen, needleleaf  
      
      ! resistance for uptake by leaves, twigs, and other exposed surfaces, O3
      rclO3 = 1000D0 ! (s/m) Summer, evergreen, needleleaf    
        
      ! Diffusion coefficients of selected gases
      DiffusivityH2O = 0.234D-4 ! Diffusion coefficient of water vapor in air (m^2/s), table 16.2 Seinfeld and Pandis
      
      ! ratio between diffusivity of water vapor and SO2, O3 or HNO3 from table 19.4 Seinfeld and Pandis
      D_ratio_SO2 = 1.89D0
      D_ratio_O3 = 1.63D0
      D_ratio_HNO3 = 1.87D0
      D_ratio_isoprene = 2.7D0 ! Estimated 
      D_ratio_apinene = 4D0  ! Estimated 
      
      
      DiffusivitySO2 = DiffusivityH2O/D_ratio_SO2    ! Diffusivity of SO2 (m^2/s)
      DiffusivityO3 = DiffusivityH2O/D_ratio_O3      ! Diffusivity of O3 (m^2/s)
      DiffusivityHNO3 = DiffusivityH2O/D_ratio_HNO3  ! Diffusivity of HNO3 (m^2/s)
      Diffusivityisoprene = DiffusivityH2O/D_ratio_isoprene  ! Diffusivity of isoprene (m^2/s)
      Diffusivityapinene = DiffusivityH2O/D_ratio_apinene  ! Diffusivity of apinene (m^2/s)
      
      ! Calculate the aerodynamic resistance for O3, SO2, HNO3, isoprene & a-pinene (ra) in similar way as
      ! for particles: 
       		
      ! Calculate the quasi-laminar resistance for O3, SO2, HNO3, isoprene & a-pinene (rb):      
      
      ! Calculation of surface resistance for O3, SO2, HNO3, isoprene & a-pinene (rc):                         
      
       ! Effective Henry's lay const:
       H_effSO2 = 1D5   ! M atm^-1
       H_effO3 = 1D-2   ! M atm^-1
       H_effHNO3 = 1D14 ! M atm^-1
       H_effisoprene = 1.2D-2 ! M atm^-1
       H_effapinene = 3D-2 ! M atm^-1
       
       ! Noramlized reactivity, table 19.4 from Seinfeld and Pandis, 2006:
       f0_SO2 = 0D0
       f0_O3 = 1D0
       f0_HNO3 = 0D0
       f0_isoprene = 0D0
       f0_apinene = 0D0
      
      ! Calculate the bulk canopy stomatal resistance (rst)
      
      ! Calculate the combined stomatal and mesophyll resistance (rsm):
      
      ! Calculate the resistance of the outer surfaces in the upper canopy (rlu):
      
      ! Calculate the resistance to transfer by buoyant convection (rdc):
      
      ! Calculate the resistance of the exposed surfaces in the lower portions of 
      ! structures of the canopy (rcl): 
              
      ! Calculate the resistance of the exposed surfaces on the groud 
      !(soil,leaf litter, ground) (rgs):
      
      ! Combine all resistances in order to get the total surface resistance 
      ! for O3, SO2, HNO3, isoprene and a-pinene (rc):
       
      ! Finally calculate the dry deposition velocity of SO2, O3, HNO3, isoprene and a-pinene:                                        

END SUBROUTINE dry_dep_velocity

END MODULE aerosol_mod
