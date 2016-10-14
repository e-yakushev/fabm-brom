#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: 
!
! !INTERFACE:
   module fabm_niva_brom_metals_ni
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
!  Original author(s): Jorn Bruggeman, Svetlana Pakhomova, Evgeniy Yakushev
!

! !PUBLIC DERIVED TYPES:
   type,extends(type_base_model),public :: type_niva_brom_metals_ni
!     Variable identifiers
!      type (type_state_variable_id)        :: id_Cu,id_CuS,id_Cu_prt
!      type (type_state_variable_id)        :: id_Ba, id_BaSO4, id_SO4     
      type (type_state_variable_id)        ::  id_H2S, id_O2, id_Fe2,id_Fe3,id_Baae,id_Bhae,id_Baan,id_Bhan        
      type (type_dependency_id)            :: id_temp,id_par,id_depth
      type (type_state_variable_id)        :: id_Ni,id_Ni2S,id_Ni_biota, id_Ni_POM, id_Ni_DOM, id_Ni_miner
      

!      type (type_state_variable_id)        :: id_Hg0,id_Hg2,id_MeHg,id_HgS      
!      type (type_horizontal_dependency_id) :: id_windspeed, id_pCO2a
      !---- Ba---------!
!      real(rk) :: K_BaSO4, K_BaSO4_form, K_BaSO4_diss  
!     Model parameters
      !---- Hg---------!
      real(rk) :: K_hg2_mehg , K_mehg_hg2, K_hg2_hg0, K_hg0_hg2
      real(rk) :: K_mehg_irr_degr ! photo-degradation
      real(rk) :: K_hg0_irr_ox ! photo-oxidation of hg0    
      real(rk) :: K_hg2_irr_red ! photo-reduction of Hg2      
      real(rk) :: K_HgS, K_hgs_form, K_hgs_ox, K_hgs_diss, O2s_nf      
      !----Point source parameters----!  
!      real(rk) :: Q_source, R_layer, Ba_source, Fe2_source, O2_source, Hg_source     
       !---- Sinking---!      
      real(rk) :: Wsed, Wphy, Whet, Wbact, Wm            
      
   contains
      procedure :: initialize
      procedure :: do
!      procedure :: do_surface
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
   class (type_niva_brom_metals_ni), intent(inout), target :: self
   integer,                     intent(in)            :: configunit
!
! !REVISION HISTORY:
!  Original author(s): 
!
!EOP
!-----------------------------------------------------------------------
!BOC
!-----------------------------------------------------------------------   
! Parameters, i.e. rate constants  
    call self%get_parameter(self%Wsed,  'Wsed', '[1/day]',  'Rate of sinking of detritus (POP, PON)',           default=5.00_rk)     
    call self%get_parameter(self%Wphy,  'Wphy', '[m/day]',  'Rate of sinking of Phy',                           default=0.10_rk)    
    call self%get_parameter(self%Whet,  'Whet', '[m/day]',  'Rate of sinking of Het',                           default=1.00_rk)     
    call self%get_parameter(self%Wbact,  'Wbact','[1/day]',  'Rate of sinking of bacteria (Bhae,Baae,Bhan,Baan)',default=0.4_rk)
    call self%get_parameter(self%Wm,     'Wm',   '[1/day]',  'Rate of accelerated sinking of particles with settled Mn hydroxides', default=7.0_rk) 

   call self%get_parameter(self%O2s_nf, 'O2s_nf', '[uM O]','half saturation for nitrification',default=4.488_rk)

      !---- Ni---------!      
   call self%register_state_variable(self%id_Ni,    'Ni',    'mmol/m**3','Ni',         minimum=0.0_rk)
   call self%register_state_variable(self%id_Ni_biota,'Ni_biota','mmol/m**3','Ni_biota', minimum=0.0_rk,vertical_movement=-self%Wphy/86400._rk)
   call self%register_state_variable(self%id_Ni_POM,'Ni_POM','mmol/m**3','Ni_POM', minimum=0.0_rk,vertical_movement=-self%Wsed/86400._rk)
   call self%register_state_variable(self%id_Ni_DOM,'Ni_DOM','mmol/m**3','Ni_DOM', minimum=0.0_rk)
   call self%register_state_variable(self%id_Ni_miner,'Ni_miner','mmol/m**3','Ni_miner', minimum=0.0_rk,vertical_movement=-self%Wsed/86400._rk)
   call self%register_state_variable(self%id_Ni2S,  'Ni2S',  'mmol/m**3','Ni2S',       minimum=0.0_rk,vertical_movement=-self%Wm/86400._rk)

   call self%register_state_dependency(self%id_O2,  'O2',    'mmol/m**3', 'O2')
   call self%register_state_dependency(self%id_Fe2, 'Fe2',   'mmol/m**3', 'Fe2')
   call self%register_state_dependency(self%id_Fe3, 'Fe3',   'mmol/m**3', 'Fe3')
   call self%register_state_dependency(self%id_H2S, 'H2S',   'mmol/m**3',  'H2S')    
   call self%register_state_dependency(self%id_Baae, 'Baae', 'mmol/m**3','Aerobic Autotrophic Bacteria')
   call self%register_state_dependency(self%id_Bhae, 'Bhae', 'mmol/m**3','Aerobic Heterotrophic Bacteria')
   call self%register_state_dependency(self%id_Baan, 'Baan', 'mmol/m**3','Anaerobic Autotrophic Bacteria')
   call self%register_state_dependency(self%id_Bhan, 'Bhan', 'mmol/m**3','Anaerobic Heterotrophic Bacteria')
   
 !  call self%register_dependency(self%id_pco2,'pCO2','ppm','pCO2')
   
   call self%register_dependency(self%id_temp, standard_variables%temperature) 
   call self%register_dependency(self%id_depth,standard_variables%pressure)  
   call self%register_dependency(self%id_par, standard_variables%downwelling_photosynthetic_radiative_flux)
!   call self%register_dependency(self%id_Izt,'Izt','W/m2','downwelling_photosynthetic_radiative_flux')
   
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
   class (type_niva_brom_metals_ni),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_
!
! !REVISION HISTORY:
!  Original author(s): 
!
! !LOCAL VARIABLES:
   real(rk) ::  Ni, Ni2S, Ni_biota, Ni_POM, Ni_DOM, Ni_miner
   real(rk) ::  temp, O2, Fe2, Fe3, depth
!   real(rk) ::  Hg0, Hg2, MeHg, HgS
   real(rk) ::  Iz, Bhan, Baan, H2S !, Om_HgS, Hg2_flux 
!   real(rk) ::  hgs_form, hgs_diss, hgs_ox, hg2_hg0, hg0_hg2, hg0_irr_ox, hg2_irr_red, mehg_irr_degr, hg2_mehg, mehg_hg2
!EOP 
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Environment
   _GET_(self%id_temp,temp)              ! temperature
   _GET_(self%id_depth,depth)            ! local photosynthetically active radiation
   
   _GET_(self%id_par,Iz)              ! local photosynthetically active radiation
!     _GET_(self%id_Izt,Izt)              ! downwelling_photosynthetic_radiative_flux  
   
    ! Our own state variables
    !_GET_(self%id_Hg0,Hg0)
    !_GET_(self%id_Hg2,Hg2)    
    !_GET_(self%id_MeHg,MeHg)
    !_GET_(self%id_HgS,HgS)      

    _GET_(self%id_Baan,Baan)   
    _GET_(self%id_Bhan,Bhan)
    _GET_(self%id_H2S,H2S)
    _GET_(self%id_Fe2,Fe2)
    _GET_(self%id_Fe3,Fe3)
    _GET_(self%id_O2,O2)   
    _GET_(self%id_Ni,Ni) 
    _GET_(self%id_Ni2S,Ni2S) 
    _GET_(self%id_Ni_biota,Ni_biota)
    _GET_(self%id_Ni_POM,Ni_POM)
    _GET_(self%id_Ni_DOM,Ni_DOM)
    _GET_(self%id_Ni_miner,Ni_miner)

!    ! Hg species (Knigthes 2008)
!!% Hg0 bioreduction  Hg0 -> Hg2+  ()
!    hg0_hg2=self%K_hg0_hg2*Hg0         
!!% Hg2 biooxydation  Hg2+ + 0.5O2 + 2H+-> Hg0 + H2O   ()
!    hg2_hg0=self%K_hg2_hg0*Hg2 * 0.5*(1.+tanh(o2-self%O2s_nf))  
!!% Hg2 methylation Hg2+  -> MeHg   ()
!    hg2_mehg=self%K_hg2_mehg*Hg2*Bhan
!!% MeHg demethylation MeHg  -> Hg2+   ()
!    mehg_hg2=self%K_mehg_hg2*MeHg
!
!    Om_HgS=H2S*Hg2/(self%K_HgS) 
!!% HgS formation Hg2+ + H2S -> HgS + 2H+ ()
!    hgs_form=self%K_hgs_form*max(0._rk,(Om_HgS-1._rk))
!    if (Hg2<0.0000001.or.H2S<0.00001) hgs_form=0.
!!% HgS dissolution  HgS + 2H+ -> Hg2+ + H2S   ()
!    hgs_diss=self%K_hgs_diss*HgS*max(0._rk,(1._rk-Om_HgS))      
!!% HgS oxydation  HgS + 2O2 -> Hg2+ + SO4   ()
! !   hgs_diss=self%K_hgs_ox*HgS *O2     
!!% Hg2 biooxydation  Hg2+ + 0.5O2 + 2H+-> Hg0 + H2O   ()
!    hg2_hg0=self%K_hg2_hg0*Hg2*O2     
!!% Hg2 photo reduction  Hg2+ -> Hg0   ()
!    hg2_irr_red=self%K_hg2_irr_red*Hg2 *Iz/80.
!!% Hg0 photo oxydation  Hg0 -> Hg2+   ()
!    hg0_irr_ox=self%K_hg0_irr_ox*Hg0 *Iz/80.
!!% Hg2 methylation Hg2+  -> MeHg   ()
!    hg2_mehg=self%K_hg2_mehg*Hg2*Bhan
!!% MeHg demethylation MeHg  -> Hg2+   ()
!    mehg_hg2=self%K_mehg_hg2*MeHg
!!% MeHg photo degradation MeHg  -> Hg2+   
!    mehg_irr_degr=self%K_mehg_irr_degr*MeHg
    

   _SET_ODE_(self%id_Ni,0.0_rk)
   _SET_ODE_(self%id_Ni2S,0.0_rk)
   _SET_ODE_(self%id_Ni_biota,0.0_rk)
   _SET_ODE_(self%id_Ni_POM,0.0_rk)
   _SET_ODE_(self%id_Ni_DOM,0.0_rk)
   _SET_ODE_(self%id_Ni_miner,0.0_rk)
   !_SET_ODE_(self%id_Hg0, -hg0_hg2+hg2_hg0)
   !_SET_ODE_(self%id_Hg2, hg0_hg2-hg2_hg0-hg2_mehg+mehg_hg2-hgs_form+hgs_diss)    
   !_SET_ODE_(self%id_MeHg,hg2_mehg-mehg_hg2)
   !_SET_ODE_(self%id_HgS,+hgs_form-hgs_diss)
   _SET_ODE_(self%id_Baan,0.)   
   _SET_ODE_(self%id_Bhan,0.)
   _SET_ODE_(self%id_Fe3,0.)
   _SET_ODE_(self%id_H2S,0.)
   _SET_ODE_(self%id_O2,0.)      

! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine do
!EOC
!! !INTERFACE:
!   subroutine do_surface(self,_ARGUMENTS_DO_SURFACE_)
!!
!! !DESCRIPTION:
!! 
!! !INPUT PARAMETERS:
!   class (type_niva_brom_metals_ni),intent(in) :: self
!   _DECLARE_ARGUMENTS_DO_SURFACE_
!!
!! !LOCAL VARIABLES:
!   real(rk)                   :: Ba, Q_Ba
!
!   _HORIZONTAL_LOOP_BEGIN_
!      _GET_(self%id_Ba,Ba)
!   
!!/*---------------------------------------------------Ba flux from the surface
!                Q_Ba = 0.1
!
!   _SET_SURFACE_EXCHANGE_(self%id_Ba,Q_Ba)
!
!   _HORIZONTAL_LOOP_END_
!   
!   end subroutine 
!!-----------------------------------------------------------------------
!! !INTERFACE:
!   subroutine do_surface(self,_ARGUMENTS_DO_SURFACE_)
!!
!! !DESCRIPTION:
!! Sea water Hg(0) exchange.   
!   
!! !INPUT PARAMETERS:
!   class (type_niva_brom_metals_ni),intent(in) :: self
!   _DECLARE_ARGUMENTS_DO_SURFACE_
!!
!! !LOCAL VARIABLES:
!   real(rk) :: pCO2w, xk, Ox, Q_Hg0, Hg0
!   real(rk) :: temp, Kc0, salt
!   real(rk) :: Sc, TK, fwind !PML
!   real(rk) :: pCO2a
!   real(rk) :: windspeed
!
!   _HORIZONTAL_LOOP_BEGIN_
!      _GET_(self%id_temp,temp)              ! temperature
!!      _GET_(self%id_salt,salt)              ! salinity
!!      _GET_(self%id_pCO2w,pCO2w)
!      _GET_HORIZONTAL_(self%id_windspeed,windspeed)
!      _GET_HORIZONTAL_(self%id_pCO2a,pCO2a)
!!      _GET_HORIZONTAL_(self%id_Hg0,Hg0)
!      
!      TK=(temp + 273.15)
!
!! PML
!! calculate the Scmidt number and unit conversions
!          Sc=2073.1-125.62*temp+3.6276*temp**2.0-0.0432190*temp**3.0
!          fwind =  (0.222d0 * windspeed**2d0 + 0.333d0 * windspeed)*(Sc/660.d0)**(-0.5)
!          fwind=fwind*24.d0/100.d0   ! convert to m/day
!
!! flux depends on the difference in partial pressures, wind and henry
!! here it is rescaled to mmol/m2/d
!!          flux = fwind * HENRY * ( PCO2A - PCO2W ) * dcf      
!      
!      
!      Q_Hg0= fwind * ( pCO2a/100.- max(0e0,Hg0))
!
!      _SET_SURFACE_EXCHANGE_(self%id_Hg0,Q_Hg0)
!
!   _HORIZONTAL_LOOP_END_
!
!   end subroutine do_surface
end module