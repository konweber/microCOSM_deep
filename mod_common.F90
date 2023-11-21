! -*- f90 -*-
MODULE MOD_COMMON
!variable declarations etc for microCOSM model
USE MOD_PRECISION
USE MOD_BOXES
IMPLICIT NONE

! timestepping variables
INTEGER                  :: nstepmax
REAL(KIND=wp), PARAMETER :: dt           = 86400._wp !/ 24._wp
REAL(KIND=wp), PARAMETER :: sperd        = 86400._wp
REAL(KIND=wp), PARAMETER :: dperyr       = 365._wp
REAL(KIND=wp), PARAMETER :: speryr       = 31536000._wp
! carbon quota per cell of prokaryotic biomass (g C cell-1)
REAL(KIND=wp), PARAMETER :: qc_prokar    = 1.0e-11_wp
REAL(KIND=wp), PARAMETER :: weight_c     = 12._wp

! Conversion factors
REAL(KIND=wp), PARAMETER :: conv_molkg_molm3 = 1024.5_wp 
REAL(KIND=wp), PARAMETER :: conv_cmhr_ms     = 1._wp/3.6e5_wp

REAL(KIND=wp), PARAMETER :: permil       = 1._wp/conv_molkg_molm3
REAL(KIND=wp), PARAMETER :: umolkg2molm3 = conv_molkg_molm3 * 1.e-6_wp 
REAL(KIND=wp), PARAMETER :: nmolkg2molm3 = conv_molkg_molm3 * 1.e-9_wp 
REAL(KIND=wp), PARAMETER :: uatm2atm     = 1.e-6_wp 
REAL(KIND=wp), PARAMETER :: molps2gtcyr  = 106._wp * 12._wp * speryr * 1.e-15_wp
! NEED TO ADD CONVERSION FACTOR BETWEEN PROKARYOTIC BIOMASS AND ABUNDANCE
REAL(KIND=wp), PARAMETER :: cellsmuL2cellsm3 = 1.0e9_wp
REAL(KIND=wp), PARAMETER :: cellsm32molm3 = qc_prokar/weight_c

REAL(KIND=wp) :: zero, one, two, three, four, five,                      &
                    six, seven, eight, nine
REAL(KIND=wp) :: ten, hundred, thousand

REAL(KIND=wp), DIMENSION(nbox, nbox) :: K, R, P

REAL(KIND=wp) :: psi, dif
 
! biogeochemical tracers internal
! Modification: added LDOC and prokaryotic biomass
REAL(KIND=wp), DIMENSION(nbox) :: theta, salt, dic, alk, po4, no3, fet, lt, ldoc, pb    

! extra biogeochem...
REAL(KIND=wp), DIMENSION(nbox) :: ph, sit, carb, feprime, bioP
REAL(KIND=wp)                  :: pstar 

! some arrays for average accumulators
REAL(KIND=wp), DIMENSION(nbox) :: thetaM, saltM, exportM, pco2M,          &
                                  dicM, alkM, po4M, no3M, fetM, ltM, sitM 
! Addition of LDOC and prokaryotic biomass
REAL(KIND=wp), DIMENSION(nbox) :: ldocM, pbM                                  
REAL(KIND=wp)                  :: timeM, pco2A, pstarM

! time derivatives 
REAL(KIND=wp), DIMENSION(nbox) :: dthetadt, dsaltdt
REAL(KIND=wp), DIMENSION(nbox) :: ddicdt, dalkdt, dpo4dt, dno3dt,         &
                                  dfetdt, dltdt, dsitdt
! Addition of LDOC and prokaryotic biomass time derivatives
REAL(KIND=wp), DIMENSION(nbox) :: dldocdt, dpbdt  

! Redfield ratios      
REAL(KIND=wp) :: rCP, rNP, rPO2, rCN, rCO2, rCFe, rSIP, rCACO3

REAL(KIND=wp), DIMENSION(nbox) :: pco2ocean, fluxCO2  
REAL(KIND=wp)                  :: netco2flux, atmos_moles, atmos_carbon,  &
                                  pco2atmos
         
! DIC gas exchange piston velocity coefficient
REAL(KIND=wp)                  :: Kwexch_av
REAL(KIND=wp), DIMENSION(nbox) :: wind, fopen

! Iron cycle parameters 
! atomic weight of iron = 56
REAL(KIND=wp)                  :: weight_fe, fe_sol, beta, Kscav, relaxfe  
! iron input rate
REAL(KIND=wp), DIMENSION(nbox) :: fe_depo, fe_pptmask
! Dynamic Ligand variables
REAL(KIND=wp)                  :: gamma_lt, lt_lifetime, ligphi
REAL(KIND=wp), DIMENSION(nbox) :: dlambdadz, lambda, ligP

! export related
! half saturation constants 
REAL(KIND=wp)                  :: kfe, kpo4, kno3, klight
! Additon of half saturation constant for prokaryotic species
REAL(KIND=wp)                  :: kfe_p, kldoc_p
REAL(KIND=wp)                  :: phi ! The fraction of the total remineralized carbon that yields LDOC
REAL(KIND=wp), DIMENSION(nbox) :: export, eratio
REAL(KIND=wp)                  :: alpha 
REAL(KIND=wp), DIMENSION(nbox) :: light, ilimit, plimit, nlimit, flimit
! Addition of limitation codes for prokaryotic species
REAL(KIND=wp), DIMENSION(nbox) :: flimit_p, ldoclimit_p
! nutrient limitation codes
INTEGER                        :: lim

! prokaryotic parameters
REAL(KIND=wp)                  :: rFeC_pb, m_l, m_q, kappa, mu0, pge
REAL(KIND=wp), DIMENSION(nbox) :: pbp, pbl, ldocP, cadd_pbl

CONTAINS

!=======================================================================
SUBROUTINE COMMON_ASSIGNMENTS()
IMPLICIT NONE 

   one  =  1._wp
   two  =  2._wp
   three=  3._wp
   four =  4._wp
   five =  5._wp
   six  =  6._wp
   seven=  7._wp
   eight=  8._wp
   nine =  9._wp
   ten  = 10._wp
   hundred = 100._wp
   thousand= 1000._wp

! Redfield ratios      
   rCP  = 106._wp
   rNP  = 16._wp
   rPO2 = 170._wp     
   rCN  = rCP/rNP 
   rCO2 = rCP/rPO2 
   rCFe = rCP/1.e-3_wp 
   rSIP = 15._wp
   rCACO3 = 10.e-2_wp 

   ph   = eight

! half saturation constants 
   kfe    = 0.1e-9_wp*conv_molkg_molm3
   kpo4   = 0.1e-6_wp*conv_molkg_molm3
   kno3   = 0.1e-6_wp*conv_molkg_molm3*rNP
   klight = 30._wp

! Iron cycle parameters 
   weight_fe = 56._wp
!solubility of iron:
   fe_sol = 0.0025_wp
! conditional stability FeL: (mol kg-1)-1
   beta   = 1.0e9_wp 
! Free Fe scavenging rate: (s-1) 
   Kscav  = 1.0e-7_wp    
! relaxfe (s) 
   relaxfe = 0.01_wp * speryr  
! DIC gas exchange piston velocity coefficient
   Kwexch_av = 0.337_wp

! Prokaryotic parameters
! Prokaryotic biomass carbon to iron ratio
   rFeC_pb = 40._wp * 1.e-6_wp
   mu0 = 0.0_wp !0.1_wp/sperd ! in units of s-1
! Prokaryotic linear mortality rate
   m_l = 1.0e-10_wp ! in units of s-1
! Prokaryotic quadratic mortality rate
   m_q = 1.0e-17_wp ! treat pb in cells per m3, this is in units of m3 per cell per s
! fraction of dead prokaryotic biomass that released as LDOC
   kappa = 0.0_wp
! Prokaryotic half saturation constant for iron
   kfe_p = 0.1e-9_wp*conv_molkg_molm3
! Prokaryotic half saturation constant for LDOC
   kldoc_p =  1.0e-6_wp*conv_molkg_molm3
! Prokaryotic growth efficiency
   pge = 0.15_wp

! LDOC parameters
   phi = 0.0_wp

   
RETURN
END SUBROUTINE COMMON_ASSIGNMENTS
!=======================================================================

END MODULE MOD_COMMON
