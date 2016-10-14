#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: 
!
! !INTERFACE:
   module fabm_niva_brom_partitioning
!
! !DESCRIPTION:
!
! !USES:

   use fabm_types

   implicit none

!  default: all is private.
   private
!
! !REVISION HISTORY:!
!  Original author(s): Evgeniy Yakushev
!

! !PUBLIC DERIVED TYPES:
   type,extends(type_base_model),public :: type_niva_brom_partitioning
!     Variable identifiers
!      type (type_state_variable_id)        :: id_Cu,id_CuS,id_Cu_prt
!      type (type_state_variable_id)        :: id_Ba, id_BaSO4, id_SO4   
      type (type_state_variable_id)        :: id_Subst_dis, id_Subst_biota, id_Subst_POM, id_Subst_DOM, id_Subst_miner , id_Subst_tot 
      type (type_state_variable_id)        :: id_Phy, id_Het, id_Baae, id_Bhae, id_Baan, id_Bhan, id_NH4, id_PON, id_DON, id_Sipart
      type (type_state_variable_id)        :: id_Mn4, id_FeS, id_FeS2
      type (type_dependency_id)            :: id_temp, id_par
      type (type_diagnostic_variable_id)   :: id_Subst_tot_diss, id_Subst_tot_part        
      
      
      !---- Ba---------!
      !real(rk) :: K_BaSO4=5.    ! BaSO4 equilibrium constant (Solubility Product Constant) (uM)=5  ( 5 umol/l, Wiki,09)  
      !real(rk) :: K_BaSO4_form=1.4e-6 ! Specific rate of precipitation of BaSO4 from Ba with SO4 (1/day)=1.4e-6 (5x10-4 uM/yr,  Arndt,09)
      !real(rk) :: K_BaSO4_diss=8.e-11 ! Specific rate of dissollution of BaSO4 to Ba and SO4  (1/day)=8.e-11 (3x10-8 1/yr, Arndt,09)
!     Model parameters
       !---- Sinking---!      
      real(rk) :: Wsed, Wphy, Whet, Wbact, Wm             
!      real(rk) :: Wsed= 5. !1Rate of sinking of detritus (POP, PON)d-1 !!  Wdetr=1.5 (Savchuk, Wulff,1996),!Wdetr= 3.5; 20. (Gregoire,2000)   
      
   contains
      procedure :: initialize
      procedure :: do
   end type
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the BROM equilibrium constant model
!
! !INTERFACE:
   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
! 
!
! !INPUT PARAMETERS:
   class (type_niva_brom_partitioning), intent(inout), target :: self
   integer,                     intent(in)            :: configunit
!
! !REVISION HISTORY:
!  Original author(s): Evgeniy Yakushev
!
!EOP
!-----------------------------------------------------------------------
!BOC
       !---- Sinking---------
    call self%get_parameter(self%Wsed,  'Wsed', '[1/day]',  'Rate of sinking of detritus (POP, PON)',           default=5.00_rk)     
    call self%get_parameter(self%Wphy,  'Wphy', '[m/day]',  'Rate of sinking of Phy',                           default=0.10_rk)    
    call self%get_parameter(self%Whet,  'Whet', '[m/day]',  'Rate of sinking of Het',                           default=1.00_rk)     
    call self%get_parameter(self%Wbact,  'Wbact','[1/day]',  'Rate of sinking of bacteria (Bhae,Baae,Bhan,Baan)',default=0.4_rk)
    call self%get_parameter(self%Wm,     'Wm',   '[1/day]',  'Rate of accelerated sinking of particles with settled Mn hydroxides', default=7.0_rk) 

    call self%register_state_variable(self%id_Subst_tot, 'Subst_tot', 'mol/m**3','Subst_tot', minimum=0.0_rk)

!    call self%register_state_dependency(self%id_SO4, 'SO4', 'mmol/m**3','SO4')
     
   call self%register_state_dependency(self%id_Phy, 'Phy', 'mmol/m**3','Phy')
   call self%register_state_dependency(self%id_Het, 'Het', 'mmol/m**3','Het')
   call self%register_state_dependency(self%id_Baae, 'Baae', 'mmol/m**3','aerobic autotrophic bacteria')
   call self%register_state_dependency(self%id_Bhae, 'Bhae', 'mmol/m**3','aerobic heterotrophic bacteria')
   call self%register_state_dependency(self%id_Baan, 'Baan', 'mmol/m**3','anaerobic aurotrophic bacteria')
   call self%register_state_dependency(self%id_Bhan, 'Bhan', 'mmol/m**3','anaerobic heterotrophic bacteria')
   call self%register_state_dependency(self%id_NH4,'NH4','mmol/m**3','ammonium')
   call self%register_state_dependency(self%id_PON,'PON','mmol/m**3','particulate organic nitrogen')
   call self%register_state_dependency(self%id_DON,'DON','mmol/m**3','dissolved organic nitrogen')
   call self%register_state_dependency(self%id_Sipart,   'Sipart', 'mmol/m**3', 'Si particulate')
   call self%register_state_dependency(self%id_Mn4,   'Mn4', 'mmol/m**3', 'Mn4')
   call self%register_state_dependency(self%id_FeS,   'FeS', 'mmol/m**3', 'FeS')
   call self%register_state_dependency(self%id_FeS2,  'FeS2', 'mmol/m**3', 'FeS2')
   call self%register_state_dependency(self%id_Subst_dis,  'Subst_dis', 'mol/m**3', 'Subst_dis') !,'Subst dissolved', minimum=0.0_rk)
   call self%register_state_dependency(self%id_Subst_biota, 'Subst_biota', 'mol/m**3','Subst_biota') !, minimum=0.0_rk,vertical_movement=-self%Wphy/86400._rk)
   call self%register_state_dependency(self%id_Subst_POM, 'Subst_POM', 'mol/m**3','Subst_POM') !, minimum=0.0_rk,vertical_movement=-self%Wsed/86400._rk)
   call self%register_state_dependency(self%id_Subst_DOM, 'Subst_DOM', 'mol/m**3','Subst_DOM') !, minimum=0.0_rk)
   call self%register_state_dependency(self%id_Subst_miner, 'Subst_miner', 'mol/m**3','Subst_miner') !, minimum=0.0_rk,vertical_movement=-self%Wm/86400._rk)


   call self%register_dependency(self%id_temp,standard_variables%temperature)
   call self%register_dependency(self%id_par, standard_variables%downwelling_photosynthetic_radiative_flux)  
   
   call self%register_diagnostic_variable(self%id_Subst_tot_diss,'Subst_tot_diss','mol/m**3','Subst_tot_diss')
   call self%register_diagnostic_variable(self%id_Subst_tot_part,'Subst_tot_part','mol/m**3','Subst_tot_part')   

   ! Specify that are rates computed in this module are per day (default: per second)
   self%dt = 86400.
   
   end subroutine initialize
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: 
!
! !INTERFACE:
   subroutine do(self,_ARGUMENTS_DO_)
!
! !DESCRIPTION:
! 
!
! !INPUT PARAMETERS:
   class (type_niva_brom_partitioning),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_
!
! !REVISION HISTORY:
!  Original author(s): 
!
! !LOCAL VARIABLES:
   real(rk) :: Subst_dis, Subst_biota, Subst_POM, Subst_DOM, Subst_miner , Subst_tot 
   real(rk) :: temp, SO4  , Iz
   real(rk) :: Subst_tot_diss,Subst_tot_part
   real(rk) :: Phy, Het, Baae, Bhae, Baan, Bhan, NH4, PON, DON, Sipart
   real(rk) :: FeS, FeS2, Mn4
   real(rk) :: dSubst_dis, dSubst_biota, dSubst_POM, dSubst_DOM, dSubst_miner, dSubst_tot
   real(rk) :: dSubst_tot_diss, dSubst_tot_part
   real(rk) :: pol_bio  ! pollutant in BIOTA,"ng?"/l
   real(rk) :: pol_dom  ! pollutant in DOM, "ng?"/l
   real(rk) :: pol_pom  ! pollutant in POM, "ng?"/l
   real(rk) :: pol_dis  ! pollutant total dissolved, "ng?"/l
   real(rk) :: pol_pa  ! pollutant total particulate, "ng?"/l
   real(rk) :: pol_miner  ! pollutant in clay, "ng?"/l
   
   real(rk) :: Ptotal   ! total pollutant, "ng?"/l
   real(rk) :: Pfree     ! pollutant in dissolved INORGANIC,"ng?"/l, i.e. not partitioned  
   
   real(rk) :: sha_bio ! % share of polutant in living organisms
   real(rk) :: sha_pom ! % share of polutant in POM
   real(rk) :: sha_dom ! % share of polutant in DOM
!   real(rk) :: sha_b1nd ! % share of polutant binded with metals
   real(rk) :: sha_miner ! % share of polutant in clay (Sipart)  
   real(rk) :: sha_free ! % share of 'free' polutant 
 !  real(rk) ::  Om_BaSO4, baso4_prec, baso4_disss,ba_flux 
   real(rk) ::  K_Ce137_deg=0.00006  !d-1; for 137Cs, decays at rate that is constant (137-Cs  = 0.0230 yr-1) in all compartments (Raoul)
   real(rk) :: uMn2lip=0.0084 !coeff.to transfer POM (uM N)->(g DryWeight/l) 
   real(rk) :: rho_FeS= 5.90E7 !    # Density of FeS [mmolFe/m3] (default = 5.90E7 mmolFe/m3)
   real(rk) :: rho_FeS2=4.17E7 !    # Density of FeS2 [mmolFe/m3] (default = 4.17E7 mmolFe/m3)
   real(rk) :: rho_Mn4= 5.78E7 !    # Density of Mn4 [mmolMn/m3] (default = 5.78E7 mmolMn/m3)   
!====================================================
!  Hg
   !real(rk) ::  Kow_bio = 199526. ! part.coeff. POM(in C)/water (Allisson05)
   !real(rk) ::  Kow_pom = 199526. ! part.coeff. POM(in C)/water (Allisson05)
   !real(rk) ::  Kow_dom = 199526. ! part.coeff. POM(in C)/water (Allisson05)
! /Hg
!  Ni
   real(rk) ::  Kow_bio = 39811. ! part.coeff. POM(in C)/water (Allisson05)
   real(rk) ::  Kow_pom = 39811. ! part.coeff. POM(in C)/water (Allisson05)
   real(rk) ::  Kow_dom = 125893. ! part.coeff. POM(in C)/water (Allisson05)
   real(rk) ::  Ksorb_min = 100000. ! part.coeff. POM(in C)/water (Allisson05)   
! /Ni   
!====================================================
!EOP 
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Environment
   _GET_(self%id_temp,temp)              ! temperature
   _GET_(self%id_par,Iz)              ! local photosynthetically active radiation
    
   _GET_(self%id_Subst_dis,Subst_dis)   
   _GET_(self%id_Subst_biota,Subst_biota)   
   _GET_(self%id_Subst_POM,Subst_POM)   
   _GET_(self%id_Subst_DOM,Subst_DOM)   
   _GET_(self%id_Subst_miner,Subst_miner)   
 !  _GET_(self%id_Subst_tot,Subst_tot)
   
   _GET_(self%id_Phy,Phy)    
   _GET_(self%id_Het,Het)    
   _GET_(self%id_Baae,Baae)
   _GET_(self%id_Bhae,Bhae)    
   _GET_(self%id_Baan,Baan)    
   _GET_(self%id_Bhan,Bhan)    
   _GET_(self%id_NH4,NH4)    
   _GET_(self%id_PON,PON)    
   _GET_(self%id_DON,DON)    
   _GET_(self%id_Sipart,Sipart)   

! Let's first assume that all the polutant is dissolved INORGANIC...
!        Ptotal = pol_di + pol_pa  !total amount of pollutant 
        Ptotal = Subst_dis+ Subst_biota+ Subst_POM+ Subst_DOM+ Subst_miner  !total amount of pollutant 

!! shares of pollutant:
       if((Phy+Het+Baae+Bhae+Baan+Bhan)<=0.) then 
        sha_bio=0. 
       else
        sha_bio=uMn2lip/1000.*(Phy+Het+Baae+Bhae+Baan+Bhan)  ! Volume(weight in kg, g->kg=/1000) of BIO
       endif 
       
       if(PON<=0.) then 
        sha_pom=0. 
       else
        sha_pom=uMn2lip/1000.*PON  ! Volume(weight in kg, g->kg=/1000) of BIO
       endif     
       
       if(DON<=0.) then 
        sha_dom=0. 
       else
        sha_dom=uMn2lip/1000.*DON  ! Volume(weight in kg, g->kg=/1000) of BIO
       endif    
              
       if(Sipart<=0.) then 
        sha_miner=0. 
       else
        sha_miner=(FeS/rho_FeS + FeS2/rho_FeS2 + Mn4/rho_Mn4)/1000. ! Volume( in kg of water??) of minerals
       endif  

     sha_free = 1.-sha_bio-sha_pom-sha_dom-sha_miner ! Volume(weight in kg) of 1l of water minus volumes of org. and part. forms 
!! The free poll.conc. left as INORGANIC
      Pfree = Ptotal* sha_free /(sha_free +Kow_bio*sha_bio              &
     &            +Kow_pom*sha_pom +Kow_dom*sha_dom +Ksorb_min*sha_miner) 
!      Pfree = Ptotal* sha_free /(sha_free +Kow_bio*sha_bio              &
!     &            +Kow_pom*sha_pom +Kow_dom*sha_dom) 
!! poll.conc. partitioned to biota
     pol_bio=max(0.,Kow_bio*Pfree*sha_bio/sha_free)
!! poll.conc. partitioning to POM
     pol_pom=max(0.,Kow_pom*Pfree*sha_pom/sha_free)
!! poll.conc. partitioning to DOM
     pol_dom=max(0.,Kow_dom*Pfree*sha_dom/sha_free)
!! poll.conc. partitioning to Minerals
     pol_miner=max(0.,Ksorb_min*Pfree*sha_miner/sha_free) 
!
!!NEW dissoved and particulate pollutants:
!!          pol_di = Pfree + pol_dom
!!          pol_pa = Ptotal-pol_di 
!          
!     if(sha_free.lt.0.001) then
!          pol_di= 0.001*Ptotal
!          pol_pa = 0.999*Ptotal
!          goto 12
!     endif
!
!     if((sha_free+sha_dom).ge.1.) then
!          pol_di=Ptotal
!          pol_pa = 0.
!     else
!          pol_di = Pfree + pol_dom
!            if (pol_di.ge.Ptotal) pol_di=Ptotal
!            if(pol_di.le.0.) pol_di = 0.
!          pol_pa = Ptotal-pol_di !pol_bio + pol_pom (to increase numeric.accuracy..)
!           if(pol_pa.le.0.) pol_pa = 0.
!     endif
!12     continue


    dSubst_dis   = -Subst_dis   +(Pfree)
    dSubst_biota = -Subst_biota +(pol_bio)
    dSubst_POM   = -Subst_POM   +(pol_pom)
    dSubst_DOM   = -Subst_DOM   +(pol_dom)
    dSubst_miner = -Subst_miner +(pol_miner)
!    dSubst_tot = -Subst_tot +(1.-K_Ce137_deg)*(Pfree +pol_dom +pol_bio +pol_pom)
    dSubst_tot_diss  = dSubst_dis+dSubst_DOM
    dSubst_tot_part  = dSubst_POM+dSubst_miner
    
   _SET_ODE_(self%id_Subst_dis,  dSubst_dis)
   _SET_ODE_(self%id_Subst_biota, dSubst_biota)
   _SET_ODE_(self%id_Subst_POM, dSubst_POM)
   _SET_ODE_(self%id_Subst_DOM, dSubst_DOM)
   _SET_ODE_(self%id_Subst_miner, dSubst_miner)
!   _SET_ODE_(self%id_Subst_tot, dSubst_tot)

   _SET_DIAGNOSTIC_(self%id_Subst_tot_diss, dSubst_tot_diss)   
   _SET_DIAGNOSTIC_(self%id_Subst_tot_part, dSubst_tot_part)      
   
! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine do
!EOC

end module