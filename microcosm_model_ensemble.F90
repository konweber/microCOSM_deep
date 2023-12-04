! -*- f90 -*-
! atmosphere-ocean carbon cycle box model
! mick follows, march 2015
! convert matlab to f90 - march/june 2016
! significant work by jonathan lauderale june 2016-oct 2019
! refresh, parallelization, and expansion by jonathan lauderdale 2020-
!
!         --------------------------------
!         |          ATMOSPHERE          |
!         |                              |
!         --------------------------------
!         | SUR- |                       |
!         | FACE ->   SURFACE  2         |
!         |  1   |                   |   |
!         |---^----------------------v---|
!         |   |                          | 
!         |          DEEP OCEAN  3       |
!         |                              | 
!         |                              | 
!         --------------------------------
!=======================================================================
       PROGRAM MICROCOSM_MODEL_ENSEMBLE
!=======================================================================

       USE MOD_PRECISION
       USE MOD_BOXES
       USE MOD_MODELMAIN_ENSEMBLE
       USE MOD_SUBROUTINES

       IMPLICIT NONE

       INTEGER :: outstepmax, id, id0
       
       REAL(KIND=wp) ::                                                &
            maxyears,                                                  &
            outputyears,                                               &
            m2deg,                                                     &
!            gaovla_opt,                                                &
            ligphi_in,                                                 &
            lt_lifein,                                                 &
            alpha_yr,                                                  &
            atpco2in,                                                  &
            psi_in,                                                    &
            dif_in,                                                    &
            lt_lifetime_factor

! Input arrays (nbox dimensionesix)
       REAL(KIND=wp), dimension (nbox) ::                              & 
            dx,                                                        &
            dy,                                                        &
            dz,                                                        &
            depth,                                                     &
            latitude,                                                  &
            thin,                                                      & 
            sain,                                                      &
            cain,                                                      & 
            alin,                                                      & 
            phin,                                                      & 
            niin,                                                      & 
            fein,                                                      & 
            ltin,                                                      & 
            fe_input,                                                  &
            dldz_in,                                                   &
            wind_in,                                                   &
            fopen_in,                                                  &
            eratio_in,                                                 &
            pbin,                                                      &
            ldocin,                                                    &
            rCLig_in

       REAL(KIND=wp), dimension (nbox, nbox) ::                        & 
            Kin,                                                       &
            Rin,                                                       &
            Pin
            
! Output arrays (nbox, by timestep dimensionesix)
       REAL(KIND=wp), dimension (:,:), allocatable ::                  &
            thout,                                                     &
            sout,                                                      &
            cout,                                                      &
            aout,                                                      &
            pout,                                                      &
            nout,                                                      &
            fout,                                                      &
            lout,                                                      &
            expout,                                                    &
            ocpco2out,                                                   &
            pbout,                                                     &
            ldocout                                                  
            
       REAL(KIND=wp), dimension (:), allocatable   ::                  &
            tout,                                                      &
            psout,                                                     &
            atpco2out

       INTEGER, dimension (:), allocatable   ::                        &
            nlout


       ! Define the containers for the parameters that are cycled through
       ! in the different ensemble runs
       ! Dependent on the type of ensemble one wants to run, this must
       ! be adjusted
       REAL(KIND=wp), dimension (:), allocatable ::                    &
            rFeC_pb_ens,                                               &
            mu0_ens,                                                   &
            m_l_ens,                                                   &
            m_q_ens,                                                   &
            kappa_ens,                                                 &
            kfe_p_ens,                                                 &
            kldoc_p_ens,                                               &
            pge_ens,                                                   &
            phi_ens,                                                   &
            rCLig_ens,                                                 &
            lt_lifein_ens,                                             &
            ligphi_ens,                                                &
            lt_lifetime_factor_ens

       ! Define indices for the different ensemble parameters
       INTEGER ::                                                      &
            irFeC_pb,                                                  &
            imu0,                                                      &
            im_l,                                                      &
            im_q,                                                      &
            ikappa,                                                    &
            ikfe_p,                                                    &
            ikldoc_p,                                                  &
            ipge,                                                      &
            iphi,                                                      &
            irCLig,                                                    &
            ilt_lifein,                                                &
            iligphi,                                                   &
            ilt_lifetime_factor

       ! csv file parameters
       CHARACTER(len=255) :: filename_csv
       CHARACTER(len=10), allocatable :: headers(:)
       INTEGER :: i

       ALLOCATE(headers(45))

              headers(1)  = 'id       '
              headers(2)  = 'dt(s)    '
              headers(3)  = 't(yr)    '
              headers(4)  = 'rFeC_pb  '
              headers(5)  = 'mu0      '
              headers(6)  = 'm_l      '
              headers(7)  = 'm_q      '
              headers(8)  = 'kappa    '
              headers(9)  = 'kfe_p    '
              headers(10) = 'kldoc_p  '
              headers(11) = 'pge      '
              headers(12) = 'phi      '
              headers(13) = 'rCLig    '
              headers(14) = 'lt_lifet '
              headers(15) = 'ligphi   '
              headers(16) = 'lt_deepf '
              headers(17) = 'PB(1)    '
              headers(18) = 'PB(2)    '
              headers(19) = 'PB(3)    '
              headers(20) = 'LDOC(1)  '
              headers(21) = 'LDOC(2)  '
              headers(22) = 'LDOC(3)  '
              headers(23) = 'Fe(1)    '
              headers(24) = 'Fe(2)    '
              headers(25) = 'Fe(3)    '
              headers(26) = 'PO4(1)   '
              headers(27) = 'PO4(2)   '
              headers(28) = 'PO4(3)   '
              headers(29) = 'NO3(1)   '
              headers(30) = 'NO3(2)   '
              headers(31) = 'NO3(3)   '
              headers(32) = 'Lig(1)   '
              headers(33) = 'Lig(2)   '
              headers(34) = 'Lig(3)   '
              headers(35) = 'DIC(1)   '
              headers(36) = 'DIC(2)   '
              headers(37) = 'DIC(3)   '
              headers(38) = 'ALK(1)   '
              headers(39) = 'ALK(2)   '
              headers(40) = 'ALK(3)   '
              headers(41) = 'OCPCO2(1)'
              headers(42) = 'OCPCO2(2)'
              headers(43) = 'OCPCO2(3)'
              headers(44) = 'ATPCO2   '
              headers(45) = 'Limit    '

              
       ! Set the output filename
              filename_csv = 'ensemble18.csv'

!=======================================================================
! Time parameters
       ! Input some initial parameters
       maxyears   = 1.0e4_wp
       outputyears= 1.0e4_wp
       outstepmax = int((maxyears/outputyears)+1)
!=======================================================================
       
! allocate memory
       allocate ( tout      (outstepmax) )
       allocate ( nlout     (outstepmax) )
       allocate ( psout     (outstepmax) )
       allocate ( atpco2out (outstepmax) )
       allocate ( thout     (outstepmax,nbox) )
       allocate ( sout      (outstepmax,nbox) )
       allocate ( cout      (outstepmax,nbox) )
       allocate ( aout      (outstepmax,nbox) )
       allocate ( pout      (outstepmax,nbox) )
       allocate ( nout      (outstepmax,nbox) )
       allocate ( fout      (outstepmax,nbox) )
       allocate ( lout      (outstepmax,nbox) )
       allocate ( expout    (outstepmax,nbox) )
       allocate ( ocpco2out   (outstepmax,nbox) )
       allocate ( pbout     (outstepmax,nbox) )
       allocate ( ldocout   (outstepmax,nbox) )
       allocate ( rFeC_pb_ens  (1:1) )
       allocate ( mu0_ens      (1:3) )
       allocate ( m_l_ens      (1:4) )
       allocate ( m_q_ens      (1:1) )
       allocate ( kappa_ens    (1:1) )
       allocate ( kfe_p_ens    (1:2) )
       allocate ( kldoc_p_ens  (1:2) )
       allocate ( pge_ens      (1:1) )
       allocate ( phi_ens      (1:3) )
       allocate ( rCLig_ens    (1:1) )
       allocate ( lt_lifein_ens(1:1) )
       allocate ( ligphi_ens   (1:3) )
       allocate ( lt_lifetime_factor_ens(1:1) )
!=======================================================================

! Initialize input arguements
       thin     =      0._wp
       sain     =     34._wp
       cain     =   2150._wp
       alin     =   2350._wp
!       phin     =      1._wp
       phin     =      2._wp
!       niin     =     16._wp
       niin     =     36._wp
       fein     =      1._wp
       ltin     =      2._wp
       atpco2in =    280._wp
       pbin     =      100._wp ! in cells µL-1
       ldocin   =      2.0_wp ! in mumol kg-1
       
! Overturning and mixing rates (m3/s)
!       psi_in = 20.e6_wp
       dif_in =  0.e6_wp
       
! Wind speed (m/s)for CO2 gas fluxes
       wind_in      =   0._wp
       
! Open surface fraction in contact with atmoshpere 
!  1 => fully open, <1 => flux impeded (e.g. by sea ice)
       fopen_in     =  0._wp
       
! Gamma over lambda for ligands "optimum" value (Lauderdale et al 2020)
!       gaovla_opt   = 0._wp !4398._wp
!       gaovla_opt   = 4398._wp

! Gamma ligand production rate (in phosphate, not carbon, units)
!       gamma_in     = 5.e-5_wp*106._wp
!       gamma_in     = 5.e-5_wp*106._wp

       ligphi_in = 5e-8_wp
       rCLig = 30.0_wp ! number of carbon atoms per ligand molecule
       
! Lambda ligand lifetime (s)
       lt_lifein    = 0.0_wp !(10._wp * 365._wp * 24._wp * 3600._wp) ! 1._wp/((5.e-5_wp*106._wp/106._wp)/4398._wp)
!       lt_lifein    = 1._wp/((gamma_in/106._wp)/gaovla_opt)
       
! Dust deposition in g Fe m-2 year-1
       fe_input     =  0._wp

! Biological production maximum rate (mol P/yr)
       alpha_yr      = 6.e-6_wp

! Deep ocean box lifetime modifier to capture the gradient due to
! photodegradation near the surface and slower loss in the deep
       dldz_in       =   0._wp

! Export ratio (export production / total production)
       eratio_in    = 0._wp

!=======================================================================
! Parameters to cycle through

! Prokaryotic parameters
! Prokaryotic biomass carbon to iron ratio

       ! rFeC_pb_ens(1) = 5.0_wp * 1.e-6_wp
       rFeC_pb_ens(1) = 40._wp * 1.e-6_wp
       ! rFeC_pb_ens(3) = 80._wp * 1.e-6_wp

! Prokaryotic maximum growth rate

       ! mu0_ens(1)  = 0.0_wp
       mu0_ens(1)  = 0.01_wp * 1.0_wp/sperd ! 
       mu0_ens(2)  = 0.1_wp  * 1.0_wp/sperd ! 0.1 d-1, in units of s-1
       mu0_ens(3)  = 1.0_wp  * 1.0_wp/sperd !

       ! mu0_ens(3)  = 1.0_wp  * 1.0_wp/sperd ! 
       ! mu0_ens(4)  = 10._wp  * 1.0_wp/sperd ! 
       ! mu0_ens(5)  = 100._wp * 1.0_wp/sperd ! 
       


! Prokaryotic linear mortality rate

       ! already tried m_l = 1 d-1 which is too high, always gives wrong numbers

       ! m_l_ens(1)  = 1.0e-2_wp * 1.0_wp/sperd ! in units of s-1
       ! m_l_ens(2)  = 1.0e-3_wp * 1.0_wp/sperd ! in units of s-1
       m_l_ens(1)  = 1.0e-2_wp * 1.0_wp/sperd ! in units of s-1
       m_l_ens(2)  = 1.0e-3_wp * 1.0_wp/sperd ! in units of s-1
       m_l_ens(3)  = 1.0e-4_wp * 1.0_wp/sperd ! in units of s-1
       m_l_ens(4)  = 1.0e-5_wp * 1.0_wp/sperd ! in units of s-1
       ! m_l_ens(4)  = 1.0e-6_wp * 1.0_wp/sperd ! in units of s-1
       ! m_l_ens(3)  = 1.0e-7_wp * 1.0_wp/sperd ! in units of s-1
       ! m_l_ens(9) = 0                        ! in units of s-1

! Prokaryotic quadratic mortality rate

       m_q_ens(1)  = 1.0e-19_wp ! treat pb in cells per m3, this is in units of m3 per cell per s
       ! m_q_ens(2)  = 1.0e-20_wp ! treat pb in cells per m3, this is in units of m3 per cell per s
       ! m_q_ens(2)  = 1.0e-21_wp ! treat pb in cells per m3, this is in units of m3 per cell per s
       ! m_q_ens(4)  = 1.0e-22_wp ! treat pb in cells per m3, this is in units of m3 per cell per s
       ! m_q_ens(3)  = 1.0e-23_wp ! treat pb in cells per m3, this is in units of m3 per cell per s


! fraction of dead prokaryotic biomass that released as LDOC

       ! kappa_ens(1) = 0.0_wp
       kappa_ens(1) = 1.0_wp
       


! Prokaryotic half saturation constant for iron

       kfe_p_ens(1) = 1.0e-10_wp *conv_molkg_molm3
       kfe_p_ens(2) = 1.0e-9_wp  *conv_molkg_molm3
       ! kfe_p_ens(3) = 1.0e-8_wp  *conv_molkg_molm3

! Prokaryotic half saturation constant for LDOC

       ! kldoc_p_ens(1) = 1.0e-7_wp *conv_molkg_molm3
       kldoc_p_ens(1) = 1.0e-6_wp *conv_molkg_molm3
       kldoc_p_ens(2) = 1.0e-5_wp *conv_molkg_molm3

! Prokaryotic growth efficiency
       pge_ens(1) = 0.25_wp

! LDOC parameters
       ! phi_ens(1) = 0.0_wp
       ! phi_ens(1) = 1.0e-10_wp
       ! phi_ens(2) = 1.0e-9_wp
       ! phi_ens(3) = 1.0e-8_wp
       ! phi_ens(1) = 1.0e-7_wp
       ! phi_ens(2) = 1.0e-6_wp
       ! phi_ens(1) = 1.0e-5_wp
       ! phi_ens(2) = 1.0e-4_wp
       phi_ens(1) = 5.0e-3_wp
       phi_ens(2) = 1.0e-2_wp
       phi_ens(3) = 2.0e-2_wp
       

! Ligand parameters
       ! rCLig_ens(1) = 15.0_wp
       rCLig_ens(1) = 30.0_wp

! Ligand lifetimes for the surface ocean
       lt_lifein_ens(1) = 0.0_wp ! automatically sets the ad-hoc term to zero
       ! lt_lifein_ens(1) = 1._wp * 365._wp * 24._wp * 3600._wp
       ! lt_lifein_ens(2) = 5._wp * 365._wp * 24._wp * 3600._wp
       ! lt_lifein_ens(3) = 10._wp * 365._wp * 24._wp * 3600._wp
       ! lt_lifein_ens(4) = 50.0_wp * 365._wp * 24._wp * 3600._wp
       ! lt_lifein_ens(5) = 100.0_wp * 365._wp * 24._wp * 3600._wp
       ! lt_lifein_ens(6) = 1000.0_wp * 365._wp * 24._wp * 3600._wp

! Ligand production rates as fraction of the primary production
       ! ligphi_ens(1) = 1.0e-8_wp
       ! ligphi_ens(2) = 5.0e-8_wp
       ! ligphi_ens(3) = 1.0e-7_wp
       ligphi_ens(1) = 1.0e-5_wp
       ligphi_ens(2) = 5.0e-5_wp
       ligphi_ens(3) = 1.0e-4_wp
       ! ligphi_ens(3) = 1.0e-5_wp
       ! ligphi_ens(8) = 5.0e-5_wp
       ! ligphi_ens(9) = 1.0e-4_wp
       ! ligphi_ens(10) = 1.0e-3_wp

! Deep ocean box lifetime modifier, longer lifetime in the deep ocean
       lt_lifetime_factor_ens(1) = 1.0_wp
       ! lt_lifetime_factor_ens(2) = 5.0_wp
       ! lt_lifetime_factor_ens(3) = 10.0_wp
       ! lt_lifetime_factor_ens(4) = 50.0_wp
       ! lt_lifetime_factor_ens(5) = 100.0_wp
       ! lt_lifetime_factor_ens(6) = 1000.0_wp

! File number identifier start for the ensemble
       id0            =     57000

! Array inputs
#if defined(SIXBOX)
! For a 20SV AMOC, psi_in (i.e. Southern Ocean upwelling) needs to be 2x
       psi_in= 20.e6_wp

       dx    = [17.0e6_wp, 17.0e6_wp, 17.0e6_wp,                       &
                17.0e6_wp, 17.0e6_wp, 17.0e6_wp ]
       dy    = [ 4.0e6_wp,  4.0e6_wp,  2.0e6_wp,                       &
                 2.0e6_wp,  8.0e6_wp,  8.0e6_wp ]
       dz    = [ 100._wp, 3900._wp,  100._wp,                          & 
                3900._wp,  100._wp, 3900._wp ]
                
       depth = [       dz(1)/2._wp,                                    &
                 dz(1)+dz(2)/2._wp,                                    &
                       dz(3)/2._wp,                                    &
                 dz(3)+dz(4)/2._wp,                                    &
                       dz(5)/2._wp,                                    &
                 dz(5)+dz(6)/2._wp ]
                 
       m2deg    = 180._wp/(dy(1)+dy(3)+dy(5))  

       latitude = [((dy(1)/2._wp)+(dy(3)/2._wp)),                      &
                   ((dy(2)/2._wp)+(dy(4)/2._wp)),                      &
                   ((dy(3)/2._wp)              ),                      &
                   ((dy(4)/2._wp)              ),                      &
                   ((dy(5)/2._wp)+(dy(3)/2._wp)),                      &
                   ((dy(6)/2._wp)+(dy(4)/2._wp))                       &
                  ]
       latitude = -90._wp+(latitude*m2deg)

! define arrays (nbox*nbox long) of box connectivity for mixing and overturning (by rows)
! Box 1 mixes with box 2 and 3; 
! Box 2 mixes with box 1 and 4; 
! Box 3 mixes with box 1, 4 and 5; 
! Box 4 mixes with box 2, 3, and 6.
! Box 5 mixes with box 3 and 6.
! Box 4 mixes with box 4 and 5.
!                       Box1    Box2    Box3    Box4     Box5    Box6 
       Kin = RESHAPE([ 0.0_wp, 1.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, & ! Box1
                       1.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, & ! Box2
                       1.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 1.0_wp, 0.0_wp, & ! Box3
                       0.0_wp, 1.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, & ! Box4
                       0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, & ! Box5
                       0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 1.0_wp, 0.0_wp  & ! Box6
                     ], [ nbox, nbox ] )
       
! Box 1 is upstream of box 2; 
! Box 2 is upstream of box 1 and box 4; 
! Box 3 is upstream of box 1, 4, and 5; 
! Box 4 is upstream of box 2, 3, and 6;
! Box 5 is upstream of box 1 and box 6;
! Box 6 is upstream of box 4 and box 5.
!                       Box1     Box2     Box3     Box4     Box5     Box6  
       Pin = RESHAPE([ 0.00_wp, 1.00_wp, 0.00_wp, 0.00_wp, 0.00_wp, 0.00_wp, & ! Box1
                       0.35_wp, 0.00_wp, 0.00_wp, 0.90_wp, 0.00_wp, 0.00_wp, & ! Box2
                       0.15_wp, 0.00_wp, 0.00_wp, 0.55_wp, 0.30_wp, 0.00_wp, & ! Box3
                       0.00_wp, 0.25_wp, 1.00_wp, 0.00_wp, 0.00_wp, 1.10_wp, & ! Box4
                       0.50_wp, 0.00_wp, 0.00_wp, 0.00_wp, 0.00_wp, 0.50_wp, & ! Box5
                       0.00_wp, 0.00_wp, 0.00_wp, 0.90_wp, 0.70_wp, 0.00_wp  & ! Box6
                     ], [ nbox, nbox ] )

! define array of remineralization coefficients (by rows)
! -1 indicates all of export is lost from cell, while 
! +1 indicates all of export is remineralized (gained) by cell
! Box 1 loses export from Box 1, which is completely remineralized in Box 2
! Box 3 loses export from Box 3, which is completely remineralized in Box 4 
! Box 5 loses export from Box 5, which is completely remineralized in Box 6 
!                       Box1    Box2    Box3    Box4     Box5    Box6 
       Rin = RESHAPE([-1.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, & ! Box1
                       0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, & ! Box2
                       0.0_wp, 0.0_wp,-1.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, & ! Box3
                       0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, & ! Box4
                       0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp,-1.0_wp, 1.0_wp, & ! Box5
                       0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp  & ! Box6
                     ], [ nbox, nbox ] )
                     
! Initial conditions
       thin(1:6)= [   20.0_wp,    4.0_wp,   -1.00_wp,   -1.00_wp,   20.0_wp,    4.0_wp ]
       sain(1:6)= [   35.5_wp,   35.5_wp,   34.75_wp,   34.75_wp,   35.0_wp,   35.0_wp ]
       cain(1:6)= [ 2100.0_wp, 2400.0_wp, 2100.00_wp, 2400.00_wp, 2100.0_wp, 2400.0_wp ]
       alin(1:6)= [ 2350.0_wp, 2400.0_wp, 2300.00_wp, 2400.00_wp, 2300.0_wp, 2400.0_wp ]
       phin(1:6)= [    0.0_wp,    2.5_wp,    2.50_wp,    2.50_wp,    0.0_wp,    2.5_wp ]
       niin(1:6)= [    0.0_wp,   36.0_wp,   30.00_wp,   36.00_wp,    0.0_wp,   36.0_wp ]

! Wind speed (m/s)for CO2 gas fluxes
       wind_in(1:6) = [ 5._wp, 0._wp, 10._wp, 0._wp, 5._wp, 0._wp ]
       
! Open surface fraction in contact with atmoshpere 
!  1 => fully open, <1 => flux impeded (e.g. by sea ice)
       fopen_in(1:6)= [ 1._wp, 0._wp, 1._wp, 0._wp, 1._wp, 0._wp ]

! Export ratio is smaller than 1 for the surface boxes
! Export ratio is 1 for the deep boxes
       eratio_in(1:6)= [ 0.1_wp, 0.5_wp, 0.1_wp, 0.5_wp, 0.1_wp, 0.5_wp ]
! ==============================================================================
! Surface values need ajustment !!!
! ==============================================================================       
       
! Dust deposition in g Fe m-2 year-1
! Hydrothermal vent input of 1 Gmol/yr (Tagliabue et al., 2010)
! mol Fe/yr * g/mol * 1/area  == g Fe m-2 year-1....
!divide by 2.5e-3 because fe_sol is multiplied again within model.
       fe_input(1:6)= [ 1.5e-1_wp, (1.e9_wp*56._wp)/(2.5e-3_wp*dx(2)*dy(2)),   &
                        1.5e-3_wp, (1.e9_wp*56._wp)/(2.5e-3_wp*dx(4)*dy(4)),   &
                        7.5e-2_wp, (1.e9_wp*56._wp)/(2.5e-3_wp*dx(6)*dy(6)) ]

! Deep ocean box lifetime modifier to capture the gradient due to
! photodegradation near the surface and slower loss in the deep
! Modification: first order loss in the deep is set to 0
       dldz_in(1:6)  = [ 1._wp, 0._wp, 1._wp, 0._wp, 1._wp, 0._wp ]
       
#elif defined(FOURBOX)
! For a 20SV AMOC, psi_in (i.e. Southern Ocean upwelling) needs to be 2x
       psi_in= 40.e6_wp

       dx    = [17.0e6_wp, 17.0e6_wp, 17.0e6_wp, 17.0e6_wp]
       dy    = [ 1.0e6_wp,  3.0e6_wp, 12.0e6_wp, 16.0e6_wp]
       dz    = [50._wp,    50._wp,    50._wp,     5050._wp]
       depth = [       dz(1)/2._wp,                                    &
                       dz(2)/2._wp,                                    &
                       dz(3)/2._wp,                                    &
                 dz(1)+dz(4)/2._wp]
                 
       m2deg    = 180._wp/dy(4)  

       latitude = [(           +(dy(1)/2._wp)),                        &
                   (dy(1)      +(dy(2)/2._wp)),                        &
                   (dy(1)+dy(2)+(dy(3)/2._wp)),                        &
                   (           +(dy(4)/2._wp))                         &
                  ]
       latitude = -90._wp+(latitude*m2deg)

! define arrays (nbox*nbox long) of box connectivity for mixing and overturning (by rows)
! Box 1 mixes with box 2 and 4; 
! Box 2 mixes with box 1, 3 and 4; 
! Box 3 mixes with box 2 and 4; 
! Box 4 mixes with box 1, 2, and 3.
!                       Box1    Box2    Box3    Box4 
       Kin = RESHAPE([ 0.0_wp, 1.0_wp, 0.0_wp, 1.0_wp,                 & ! Box1
                       1.0_wp, 0.0_wp, 1.0_wp, 1.0_wp,                 & ! Box2
                       0.0_wp, 1.0_wp, 0.0_wp, 1.0_wp,                 & ! Box3
                       1.0_wp, 1.0_wp, 1.0_wp, 0.0_wp                  & ! Box4
                     ], [ nbox, nbox ] )
       
! Box 1 is upstream of box 4 (ie AABW downwelling); 
! Box 2 is upstream of box 1 and 3 (ie Antarctic divergence); 
! Box 3 is upstream of box 4 (ie NADW downwelling); 
! Box 4 is upstream of box 2 (ie SO upwelling).
!                       Box1    Box2    Box3    Box4 
       Pin = RESHAPE([ 0.0_wp, 0.0_wp, 0.0_wp, 0.5_wp,                 & ! Box1
                       0.5_wp, 0.0_wp, 0.5_wp, 0.0_wp,                 & ! Box2
                       0.0_wp, 0.0_wp, 0.0_wp, 0.5_wp,                 & ! Box3
                       0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp                  & ! Box4
                     ], [ nbox, nbox ] )

! define array of remineralization coefficients (by rows)
! -1 indicates all of export is lost from cell, while 
! +1 indicates all of export is remineralized (gained) by cell
! Box 1 loses export from Box 1, which is completely remineralized in Box 4 (Box 2 is adjacent)
! Box 2 loses export from Box 2, which is also completely remineralized in Box 4 (Box 1 is adjacent)
!                       Box1    Box2    Box3    Box4 
       Rin = RESHAPE([-1.0_wp, 0.0_wp, 0.0_wp, 1.0_wp,                 & ! Box1
                       0.0_wp,-1.0_wp, 0.0_wp, 1.0_wp,                 & ! Box2
                       0.0_wp, 0.0_wp,-1.0_wp, 1.0_wp,                 & ! Box3
                       0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp                  & ! Box4
                     ], [ nbox, nbox ] )
                     
! Initial conditions
       thin(1:4)= [   -1._wp,    2._wp,   20._wp  ,    4._wp   ]
       sain(1:4)= [   35._wp,   34._wp,   35.50_wp,   34.75_wp ]
       cain(1:4)= [ 2100._wp, 2100._wp, 2100._wp  , 2400._wp   ]
       alin(1:4)= [ 2350._wp, 2300._wp, 2300._wp  , 2400._wp   ]
       phin(1:4)= [    2._wp,    2._wp,    0.0_wp ,    2.5_wp  ]
       niin(1:4)= [   32._wp,   32._wp,    0._wp  ,   36._wp   ]

! Wind speed (m/s)for CO2 gas fluxes
       wind_in(1:4) = [ 10._wp, 10._wp, 5._wp, 0._wp ]
       
! Open surface fraction in contact with atmoshpere 
!  1 => fully open, <1 => flux impeded (e.g. by sea ice)
       fopen_in(1:4)= [ 0.5_wp, 1._wp, 1._wp, 0._wp ]

! Export ratio is smaller than 1 for the surface boxes
! Export ratio is 1 for the deep boxes
       eratio_in(1:4)= [ 0.25_wp, 0.5_wp, 0.5_wp, 0.1_wp ]       
       
! Dust deposition in g Fe m-2 year-1
       fe_input(1:4)= [ 1.5e-3_wp, 1.5e-3_wp, 1.5e-1_wp,               &
! Hydrothermal vent input of 1 Gmol/yr (Tagliabue et al., 2010)
! mol Fe/yr * g/mol * 1/area  == g Fe m-2 year-1....
!divide by 2.5e-3 because fe_sol is multiplied again within model.
        (1.e9_wp*56._wp)/(2.5e-3_wp*dx(4)*dy(4)) ]

! Deep ocean box lifetime modifier to capture the gradient due to
! photodegradation near the surface and slower loss in the deep
! Modification: first order loss in the deep is set to 0
       dldz_in(1:4)  = [ 1._wp, 1._wp, 1._wp, 0._wp ]
#else

! Default to the three box model
       psi_in = 20.e6_wp
       
       dx    = [ 17.0e6_wp, 17.0e6_wp,   17.0e6_wp ]
       dy    = [  4.0e6_wp, 12.0e6_wp,   16.0e6_wp ]  
       dz    = [ 50.0_wp  , 50.0_wp  , 5050.0_wp   ]

       depth = [ dz(1)/2._wp, dz(2)/2._wp, dz(1)+dz(3)/2._wp]

       m2deg    = 180._wp/dy(3)  
       latitude = [(dy(1)       /2._wp) ,                              &
                   (dy(1)+(dy(2)/2._wp)),                              &
                   (dy(3)       /2._wp)                                &
                   ]
       latitude = -90._wp+(latitude*m2deg)

! Mixing mask array
       Kin = RESHAPE( [ 0._wp, 1._wp, 1._wp,                           &
                        1._wp, 0._wp, 1._wp,                           &
                        1._wp, 1._wp, 0._wp ],                         &
                      [ nbox, nbox ] )
       
! Overturning mask array
       Pin = RESHAPE([ 0._wp, 1._wp, 0._wp,                            &
                       0._wp, 0._wp, 1._wp,                            &
                       1._wp, 0._wp, 0._wp ],                          &
                        [ nbox, nbox ] )
       
! Remineralization coefficients 
! -1 indicates all of export is lost from cell, while 
! +1 indicates all of export is remineralized (gained) by cell
!  the first box (column one) loses export from box 1,
!  the second box (col two) loses export from box 2,
!  and the third box (col three) gains export from boxes 1 and 2 
       Rin = RESHAPE([ -1._wp, 0._wp, 1._wp,                           &
                        0._wp,-1._wp, 1._wp,                           &
                        0._wp, 0._wp, 0._wp ],                         &
                      [ nbox, nbox ] )

! Initial conditions
       thin(1:3)= [    2._wp,   20._wp  ,    4._wp   ]
       sain(1:3)= [   34._wp,   35.50_wp,   34.75_wp ]
       cain(1:3)= [ 2200._wp, 2100._wp  , 2400._wp   ]
       alin(1:3)= [ 2350._wp, 2350._wp  , 2400._wp   ]
       phin(1:3)= [    2._wp,    0._wp  ,    2.5_wp  ]
       niin(1:3)= [   25._wp,    0._wp  ,   35._wp   ]

! Initial concentrationesix in (u/n)mol/kg
! Here are some equilibrated values run for 100,000 yrs (round-off error notwithstanding)
! Make sure to compile without -DFIXATMPCO2
!       cain(1:3)= [ 2263.27105_wp, 2104.02729_wp, 2358.11830_wp ]
!       alin(1:3)= [ 2396.24755_wp, 2388.24068_wp, 2399.60156_wp ]
!       phin(1:3)= [ 1.85304_wp   , 0.31325_wp   , 2.49804_wp    ]
!       niin(1:3)= [ 24.68043_wp  , 0.04392_wp   , 35.00046_wp   ]
!       fein(1:3)= [ 0.01007_wp   , 0.37382_wp   , 0.55782_wp    ]
!       ltin(1:3)= [ 1.62217_wp   , 1.58451_wp   , 1.57992_wp    ]

! Wind speed (m/s)for CO2 gas fluxes
       wind_in(1:3) = [ 10._wp, 5._wp, 0._wp ]
       
! Open surface fraction in contact with atmoshpere 
!  1 => fully open, <1 => flux impeded (e.g. by sea ice)
       fopen_in(1:3)= [ 1._wp, 1._wp, 0._wp ]

! Export ratio is smaller than 1 for the surface boxes
! Export ratio is 1 for the deep boxes
       eratio_in(1:3)= [ 0.25_wp, 0.1_wp, 0.5_wp ]   

! Dust deposition in g Fe m-2 year-1
       fe_input(1:3)= [ 1.5e-3_wp, 1.5e-1_wp,                          &
! Hydrothermal vent input of 1 Gmol/yr (Tagliabue et al., 2010)
! mol Fe/yr * g/mol * 1/area  == g Fe m-2 year-1....
!divide by 2.5e-3 because fe_sol is multiplied again within model.
        (1.e9_wp*56._wp)/(2.5e-3_wp*dx(3)*dy(3)) ]
! Deep ocean box lifetime modifier to capture the gradient due to
! photodegradation near the surface and slower loss in the deep
! Modification: first order loss in the deep is set to 0
       dldz_in(1:3)  = [ 1._wp, 1._wp, 0.0_wp ]
#endif

!=======================================================================
! Initialize csv file
!=======================================================================
       CALL WRITE_CSV_HEADER(filename_csv, headers, SIZE(headers))

       DEALLOCATE(headers)
       
!=======================================================================
! Start loop over ensemble
!=======================================================================

id = id0

DO irFeC_pb = 1, size(rFeC_pb_ens)
    rFeC_pb = rFeC_pb_ens(irFeC_pb)

    DO imu0 = 1, size(mu0_ens)
        mu0 = mu0_ens(imu0)

        DO im_l = 1, size(m_l_ens)
            m_l = m_l_ens(im_l)

            DO im_q = 1, size(m_q_ens)
                m_q = m_q_ens(im_q)

                DO ikappa = 1, size(kappa_ens)
                    kappa = kappa_ens(ikappa)

                    DO ikfe_p = 1, size(kfe_p_ens)
                        kfe_p = kfe_p_ens(ikfe_p)

                        DO ikldoc_p = 1, size(kldoc_p_ens)
                            kldoc_p = kldoc_p_ens(ikldoc_p)

                            DO ipge = 1, size(pge_ens)
                                pge = pge_ens(ipge)

                                DO iphi = 1, size(phi_ens)
                                    phi = phi_ens(iphi)

                                    DO irCLig = 1, size(rCLig_ens)
                                        rCLig = rCLig_ens(irCLig)

                                        DO ilt_lifein = 1, size(lt_lifein_ens)
                                            lt_lifein = lt_lifein_ens(ilt_lifein)

                                          DO iligphi = 1, size(ligphi_ens)
                                            ligphi_in = ligphi_ens(iligphi)

                                            DO ilt_lifetime_factor = 1, size(lt_lifetime_factor_ens)
                                                lt_lifetime_factor = lt_lifetime_factor_ens(ilt_lifetime_factor)



                                                 WRITE(*,*) 'id = ', id

                                                        dldz_in(1:3)  = [ 1._wp, 1._wp, 1.0_wp / lt_lifetime_factor ]

                                                        ! part of the deep box ligands are refractory, assumed 50%
                                                        rCLig_in(1:3) = [ rCLig, rCLig, rCLig / 10.0_wp ]

                                                 call model(                                                    &
                                                        id,                                                        &
                                                        filename_csv,                                              &
                                                        maxyears,                                                  &
                                                        outputyears,                                               &
                                                        outstepmax,                                                &
                                                        dx,                                                        &
                                                        dy,                                                        &
                                                        dz,                                                        &
                                                        depth,                                                     &
                                                        latitude,                                                  &
                                                        Kin,                                                       &
                                                        Rin,                                                       &
                                                        Pin,                                                       &
                                                        psi_in,                                                    &
                                                        dif_in,                                                    &    
                                                        alpha_yr,                                                  &
                                                        ligphi_in,                                                 &
                                                        lt_lifein,                                                 &
                                                        dldz_in,                                                   &
                                                        fe_input,                                                  &
                                                        wind_in,                                                   &
                                                        fopen_in,                                                  &
                                                        rCLig_in,                                                  &
                                                        thin,                                                      &
                                                        sain,                                                      &
                                                        cain,                                                      &
                                                        alin,                                                      &
                                                        phin,                                                      &
                                                        niin,                                                      &
                                                        fein,                                                      &
                                                        ltin,                                                      &
                                                        atpco2in,                                                  &
                                                        eratio_in,                                                 &
                                                        pbin,                                                      &
                                                        ldocin,                                                    &
                                                        tout,                                                      &            
                                                        thout,                                                     &
                                                        sout,                                                      &
                                                        cout,                                                      &
                                                        aout,                                                      &
                                                        pout,                                                      &
                                                        nout,                                                      &
                                                        fout,                                                      &
                                                        lout,                                                      &
                                                        expout,                                                    &
                                                        nlout,                                                     &
                                                        psout,                                                     &
                                                        ocpco2out,                                                 &
                                                        atpco2out,                                                 &
                                                        pbout,                                                     &
                                                        ldocout                                                    &
                                                 )

                                                 id = id + 1
                                              END DO ! ilt_lifetime_factor
                                          END DO ! ligphi_ens
                                        END DO ! lt_lifein
                                    END DO ! rCLig_ens
                                END DO ! phi_ens
                            END DO ! pge_ens
                        END DO ! kldoc_p_ens
                    END DO ! kfe_p_ens
                END DO ! ikappa
            END DO ! m_q_ens
        END DO ! m_l_ens
    END DO ! mu0_ens
END DO ! rFeC_pb_ens


!=======================================================================
       END PROGRAM MICROCOSM_MODEL_ENSEMBLE
!=======================================================================