

! $Id: ras.F90,v 1.24.2.8.2.1.4.2.2.1.14.1.2.1.2.1.12.1.34.1.32.1.2.1.22.1.4.1.20.1.2.3.6.5.24.3.4.3.6.2.2.2.2.6.4.3.30.1.12.1 2016-08-12 20:17:28 ltakacs Exp $

MODULE RAS_ADM

   use ESMF
   use GEOS_Mod
   use GEOS_UtilsMod, only : DQSAT=>GEOS_DQsat
!   use module_ras
   !use aer_cloud, only: AerProps, getINsubset
!   use aer_cloud
   use RASPARAMS
   use qsat_util

   IMPLICIT NONE

   PRIVATE
   PUBLIC RASE_BWD, RASE_FWD

   character(len=32), parameter, public :: RasCodes(-2:7) = (/ &
         "Default                         ", & ! -2
         "Invalid Code                    ", & ! -1
         "Cloud detraining here is active ", & ! 0
         "PBL h < layer's h*              ", & ! 1
         "No Valid lambda                 ", & ! 2
         "Lambda out of bounds            ", & ! 3
         "A < Acrit                       ", & ! 4
         "Negative Kernel                 ", & ! 5
         "Invalid Code                    ", & ! 6
         "RH Trigger not met              "  & ! 7
         /)

CONTAINS
!         pko, plo, phio, phie, qlo, qio,                  &
!FLX, FLXC,                  &
!         CNV_CVW,                                         &
!         CNV_QC,                                          &
!         ENTLAM,                                          &
!         CLAN,                                            &
!         HHO, HSO,PRECU,                                  &
!         AEROPROPS,                                   & !!!!!AER_CLOUD
!         CNV_FICE,                                        & 
!         CNV_NICE,                                        & 
!         CNV_NDROP,                                    &   !DONIF
  SUBROUTINE RASE(idim, irun, k0, icmin, dt, doras, cpo, alhlo, alhl1, &
&   tice, gravo, seedras, sige, kcbl, wgt0, wgt1, tpert, qpert, tho, qho&
&   , uho, vho, qss, dqs, cnv_fraction, rasal2_2d, co_auto, ple, pke, &
&   clw, flxd, cnv_prc3, cnv_updfrc, rasparams, itrcr, xho, &
&   triedlev_diag, fscav, disske)
    IMPLICIT NONE
!!!!!!!!!!!!!!!======================================     
!!!!!!!!!!   Subroutine ARGact: finds the activated droplet number from Abdul_Razzak and Ghan 2000.
!      ! Tailored for GOCART AEROSOL and works with AERO input 
!      !Written by Donifan Barahona
!      !!donifan.o.barahona@nasa.gov
!!!!!!!!!!!!!!!====================================
!
!      SUBROUTINE ARGact (TEMP, WX, NCPL_ACT, NCPL_AMB,  CDUST, CSOOT, LEV, ISBASE, DDUSTAMB, DSOOTAMB, ENT_PARAM)
!         !
!         Integer, intent(in)     ::  LEV    
!         LOGICAL,  intent(in)     ::   ISBASE
!
!         REAL(8), intent(inout)     ::   TEMP, WX, ENT_PARAM
!         REAL(8), intent(out)     ::   NCPL_ACT, NCPL_AMB, CSOOT, DSOOTAMB
!         REAL(8), DIMENSION(NDUSTMAX), INTENT(OUT) :: CDUST, DDUSTAMB
!         integer                 :: INDEX, NMODES, naux       
!
!         type(AerProps)  :: AER, auxaer      
!
!
!         REAL(8)     ::      kappa, alfa, beta, Akoh, G, T, smax, fi, gi, nui, &
!                       citai, ui, aux1, PACT,  Ntot, auxx, aux, auxconc, W, alph, aseasalt_aux, f_seasalt1
!         REAL(8), dimension (30) ::  SMI, TPI, SIGI 
!
!
!         SMI=0.0      
!         TPI = 0.0
!         SIGI =2.0
!         NCPL_ACT=0.0
!         NCPL_AMB=0.0
!         CDUST=0.0
!         CSOOT=0.0
!         DDUSTAMB =1.0e-9
!         DSOOTAMB= 1.0e-9
!         W=MIN(WX*(1.0-ENT_PARAM), 20.0)    
!         call init_Aer(AER)  
!         AER  =   AERO(LEV) 
!
!         PACT=0.0 !activation probability of entrained aerosol 
!   auxconc =0.0
!    aseasalt_aux  = 0.0
!!!!!!!!!!!activate aerosol transported from cloud base 
!             NMODES =  AER_BASE%nmods
!             TPI(1:nmodes) = AER_BASE%num(1:nmodes)
!             SIGI(1:nmodes) = AER_BASE%sig(1:nmodes)                          
!
!             
!             Ntot= 0.0
!              do index = 1, nmodes 
!	              if (AER_BASE%kap(index) .gt. 0.1) Ntot =  Ntot + TPI(index)  
!              end do
!         
!
!         if ((Ntot .lt. 1.0e4) .or. (TEMP .lt. 245.0) .or. (W .lt. 0.01)) then !no activation if aerosol < 1e-4 1/cm3           
! 
!            NCPL_ACT  = 0.0
!         else
!
!            ! Calculate constants. These fits were obtained from detailed correlations of physical properties. G is actually 1/G
!            T = min(max(TEMP, 243.0), 323.0)     
!            alfa=2.8915E-08*(T**2) - 2.1328E-05*T + 4.2523E-03
!            beta=exp(3.49996E-04*T**2 - 2.27938E-01*T + 4.20901E+01)
!            G=exp(-2.94362E-06*T**3 + 2.77941E-03*T**2 - 8.92889E-01*T + 1.18787E+02)
!            Akoh= 0.66e-6/T  !from Seinfeld and Pandis (1998)
!     
!            !=======================================================
!            !Activate droplets   
!            !=======================================================
!            !Calculate maximum supersaturation according to ARG2002
!
!            auxx=0.0 
!            
!            
!            DO INDEX = 1, NMODES            
!                
!               kappa=  max(AER_BASE%kap(INDEX), 0.001)
!             
!                  SMI (INDEX) = ((0.667*Akoh/AER_BASE%dpg(INDEX))**1.5)/SQRT(2.0*kappa)   ! Critical supersat for mode I   
!                  SMI=MAX(SMI, 1.0e-5)   
!                   
!              if ((TPI(INDEX) .gt. 1e4) .and.  (kappa .gt. 0.1)) then                       
!                  fi=0.5*exp(2.5*SIGI(INDEX)) !sigi is now log(sigi)
!                  gi=1.0+0.25*SIGI(INDEX)
!                  nui=((alfa*W*G)**1.5)/(2.0*MAPL_PI*980.0*beta*TPI(INDEX))
!                  citai = 0.667*Akoh*SQRT(alfa*W*G)
!                  aux1=fi*((citai/nui)**1.5) + gi*(SMI(INDEX)*SMI(INDEX)/(nui+(3.0*citai)))**0.75
!                  aux1=aux1/(SMI(INDEX)*SMI(INDEX))      
!                  auxx=auxx+aux1                  
!                end if
!            end do
!
!  !Calculate number of activated droplets
!            if (auxx .gt. 0.0) then
!               smax = 1/sqrt(auxx)
!               auxx=0.0
!
!                   DO INDEX = 1, NMODES
!                        if ((TPI(INDEX) .gt. 1e4) .and. (AER_BASE%kap(index) .gt. 0.1)) then
!                           ui=sqrt(2.0)*log(SMI(INDEX)/smax)/3.0
!                           aux1=0.5*TPI(INDEX)*(1.0-ERFAPP(ui))
!                           auxx=auxx+aux1
!                           AER_BASE%num(index) = max(TPI(INDEX) -aux1, 0.0) !remove already activated aerosol
!                        end if                    
!                   END DO
!                  NCPL_ACT=auxx             
!            else
!                  NCPL_ACT = 0.0             
!            end if
!
!     
! 
!       end if
!
!         !now filllup dust and soot number
!         NMODES =  AER%nmods
!
!         call getINsubset(1, AER,  auxaer)
!         CDUST(1:auxaer%nmods)= auxaer%num(1:auxaer%nmods)
!         DDUSTAMB(1:auxaer%nmods)= auxaer%dpg(1:auxaer%nmods)
!         call getINsubset(2, AER,  auxaer)
!         naux = max(auxaer%nmods, 1)
!         CSOOT= sum(auxaer%num) 
!         DSOOTAMB= sum(auxaer%dpg)/naux
!
!
!         PACT=1.0 ! fraction of entrained aerosol that is activated
!         auxconc =0.0
!         aseasalt_aux  = 0.0
!
!         do index = 1, nmodes 
!	       if (AER%kap(index) .gt. 0.8)  auxconc = AER%num(index) + auxconc
!           if (AER_BASE%kap(index) .gt. 0.8)    aseasalt_aux  = aseasalt_aux  + AER_BASE%num(index)*AER_BASE%dpg(index)*AER_BASE
!%dpg(index)*1.61*MAPL_PI !assumes a fixed sigma = 2.0
!
!         end do
!       aseasalt = max(aseasalt, aseasalt_aux)
!     
!	  NCPL_AMB=auxconc !Activate  entrained aerosol with kappa>0.8   
!
!
!
!      END SUBROUTINE ARGact
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      !=====Subroutine INfreezing=========
!      ! Freeze droplets in immersion and contact ice nucleation modes, according to Barahona et al. GMD (2014)
!!!!!!!!!!!!!!!!!!!!!!
!
!      subroutine INfreezing(QL, NL, NIN, NDUST, NSOOT, INDUST, INSOOT, TEMP, PRE, WSUB_, DDUST, DSOOT)  !NIN freezing tendency
!         REAL(8), INTENT( IN) :: TEMP, NSOOT, WSUB_, PRE, QL, NL, DSOOT
!         REAL(8), DIMENSION(NDUSTMAX), INTENT (IN) ::  NDUST, DDUST
!         REAL(8), DIMENSION(NDUSTMAX), INTENT (OUT) ::  INDUST
!         REAL(8), INTENT(OUT) ::  NIN, INSOOT 
!         REAL(8), DIMENSION(NDUSTMAX) ::  INDUST_C
!
!         REAL(8) :: a, b, c , d, Tx, n05, ui, aux, nssoot, nsdust, ninbkg, SI, acorr, &
!               dnsd, dnss, coolr, WSUB, ahet, INSOOT_C
!
!         REAL(8) :: nssoot_c, nsdust_c, mfp, nslip, ndfaer, viscosity, lam_r, taux, rho_a, &
!                fdust_drop, fsoot_drop, min_ns_dust, min_ns_soot, nsss, INsea, dnsss, min_ns_seasalt 
!
!         logical :: demott, Drop_mediated
!         integer :: ix
!
!         min_ns_dust= 3.75e6 !limits ice nucleation to -12 !new 02/10/14 
!         min_ns_soot= 3.75e9 !limits ice nucleation to -18
!         min_ns_seasalt = 4.0e2 !limits ice nucleation to -5
!         
!         demott=.false.
!         INDUST=0.0
!         INSOOT=0.0
!         INDUST_C=0.0
!         INSOOT_C=0.0   
!         NIN=0.0
!         Drop_mediated = .false.
!       INsea = 0.0 ! sea salt only in immersion
!   
!   
!   ! note for sea salt we just assume that it is equal to the current droplet concentration and take the area from the calculati
!on at cloud base
!
!         ! fraction of dust and soot within droplets
!         fdust_drop= FDROP_DUST
!         fsoot_drop = FDROP_SOOT
!
!
!         WSUB=MAX(MIN(WSUB_, 10.0), 0.8)
!         coolr=5.0e-3*WSUB  !Approximation to saturated cooling rate 
!         n05=sum(NDUST)+NSOOT
!
!         if (TEMP .gt. T_ICE_MAX) then 
!            return
!         end if
!
!         if (TEMP .lt. T_ICE_ALL   ) then 
!            return
!         end if
!
!
!         if ((QL .le. 1.0e-10) .or. (NL .le. 1.0)) then 
!            return
!         end if
!
!         !Background IN
!         ! SI at water saturation
!
!         rho_a = PRE*100.0/MAPL_RGAS/TEMP 
!         ninbkg=0.0
!         SI = -1.2379e-2+3.3595 !Ice supersat at water sat. Derived from Murphy and Koop 2005
!         !if (TEMP .lt. 260.0)  ninbkg=coolr*42.8*exp(3.88*si)*0.1 !tendency in IN from background IN. Derived from Phillips et 
!al. 2007
!
!
!         Tx = max(TEMP-273.16, -38.0 )
!
!         lam_r=min((MAPL_PI*950.0*NL/rho_a/QL)**(1./3.), 1.0e8)
!         viscosity=1.8e-5*(TEMP/298.0)**0.85    ! Viscosity (kg/m/s)
!         mfp=2.0*viscosity/(PRE  &                   ! Mean free path (m)
!               *sqrt(8.0*28.96e-3/(MAPL_PI*8.314409*TEMP)))        
!
!         if ((n05 .gt.1.0) .and. (TEMP .lt. 272.0)) then
!
!            nsdust=  max(exp(-0.517*Tx + 8.934)-min_ns_dust, 0.0) !From Niemand 2012
!            nssoot= max(1.0e4*exp(-0.0101*Tx*Tx - 0.8525*Tx + 0.7667)-min_ns_soot, 0.0) !Murray (review_ 2012)
!            dnsd  = 0.517*nsdust
!            dnss  = max(-(-2.0*0.0101*Tx -0.8525)*nssoot, 0.0)
!
!            !ns in  contact. It is assumed that in contact is T-3 immersion
!            taux=max(Tx-3.0, -35.0)
!            nsdust_c= max(exp(-0.517*taux + 8.934)-min_ns_dust, 0.0) !From Niemand 2012
!            nssoot_c= max(1.0e4*exp(-0.0101*taux*taux - 0.8525*taux + 0.7667)-min_ns_soot, 0.0) !Murray (review_ 2012)
!
!            aux=0.0        
!            acorr=2.7e7 !m2/m3 correction to the area due to non sphericity and aggregation. Assumes 10 m2/g (Murray 2011)
!
!
!            DO ix=1, NDUSTMAX
!               !Immersion
!               ahet=0.52*DDUST(ix)*DDUST(ix)*DDUST(ix)*acorr*exp(4.5*log(2.0)*log(2.0)*log(2.0))    !this needs to be improved
!
!               INDUST(ix) = NDUST(ix)*exp(-nsdust*ahet)* &
!                     dnsd*coolr*ahet*fdust_drop
!               !Contact   
!               nslip =1.0+(2.0*mfp/DDUST(ix))*(1.257+(0.4*exp(-(1.1*DDUST(ix)*0.5/mfp))))! Slip correction factor               
!     
!               ndfaer =1.381e-23*TEMP*nslip*(1.0-exp(-nsdust_c*ahet)) /(12.*MAPL_PI*viscosity*DDUST(ix))             
!               INDUST_C(ix) = 2.0*MAPL_PI*ndfaer*NDUST(ix)*NL/lam_r
!
!            END DO
!
!
!            acorr=8.0e7 !m2/m3 correction to the area due to non sphericity and aggregation  Assumes 50 m2/g (Popovicheva 2003) 
!      
!            ahet =0.52*DSOOT*DSOOT*DSOOT*acorr*exp(4.5*log(2.0)*log(2.0)*log(2.0))             
!            INSOOT=fsoot_drop*NSOOT*exp(-nssoot*ahet)*dnss*ahet*coolr !
!
!            nslip =1.0+(2.0*mfp/DSOOT)*(1.257+(0.4*exp(-(1.1*DSOOT*0.5/mfp))))! Slip correction factor                    
!            ndfaer =1.381e-23*TEMP*nslip*(1.0-exp(-nssoot_c*ahet)) /(12.*MAPL_PI*viscosity*DSOOT)             
!            INSOOT_c= 2.0*MAPL_PI*ndfaer*NSOOT*NL/lam_r
!
!         ! sea salt
!         nsss =  -0.459*TEMP +128.6235 ! from Demott et al. PNAS, 2015  
!         nsss=  max(exp(nsss)-min_ns_seasalt, 0.0)           
!    	 dnsss=  max(0.459*nsss, 0.0)
!         INsea= aseasalt*dnsss*coolr 
!
! 
!         end if
!
!	    NIN =ninbkg+ INSOOT + SUM(INDUST) + INSOOT_C + SUM(INDUST_C) + INsea!
!         INSOOT=INSOOT +INSOOT_C
!         INDUST =INDUST + INDUST_C
!
!	    
!      end subroutine INfreezing
!
!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      !===================Subroutine Qgrowth======================
!      !Partitions water and ice according to Korolev and Mazin 2005. Assume spheres by now.
!!!!!!!!!!!!!!!!!!!!!
!      subroutine Qgrowth(TEMP, PRE, QICE0, NIPRE, QT, NINUC, DQIG, RIM, FNDRIM) !freezing of IN according to Demott et al 2010 (
!everything SI)
!         REAL(8), INTENT(IN) :: TEMP, QICE0, NIPRE, NINUC,  PRE, QT
!         REAL(8), INTENT(INOUT) :: DQIG, RIM, FNDRIM
!
!         !REAL(8) :: A, denI, aux, Dco, SI, denA, DQold, DQnew 
!         REAL(8) :: DIFF, DENICE, DENAIR, K1, K2, K3, SI, AUX, DC, TEFF, PLo, TEo, TC, &
!               DQnew, DQold, rho_a, LWC, IWC, qmin_rim
!
!
!         if (TEMP .gt. 272.15) then
!            DQIG =0.0
!            return
!         end if
!
!         TC=TEMP-273.0 
!         PLo = max(PRE, 10.0) !limits  of the correlations 
!         TEo = max(190.0, TEMP)
!
!         qmin_rim = 1.0e-12
!
!
!         DENICE= 1000.0*(0.9167 - 1.75e-4*TC -5.0e-7*TC*TC) !From PK 97
!         DENAIR= PLo*100.0/MAPL_RGAS/TEMP
!         DIFF=(0.211*1013.25/(PLo+0.1))*(((TEo+0.1)/273.0)**1.94)*1e-4  !From Seinfeld and Pandis 2006
!
!         K1 = EXP(7.1170e-4*TEo*TEo-0.43563*TEo+78.744) 
!         K2 = EXP(-9.933e-3*TEo+25.26)
!         K3 = EXP(7.1772e-4*TEo*TEo-0.44055*TEo+73.996)
!
!
!         AUX= 210368.0 + 131.438*TEMP - (3.32373E6/TEMP)- (41729.1*LOG(TEMP)) !From Murphy and Koop 2005
!         SI=exp(-aux/8.314/TEMP)-1.0 !ratio of pw/pi-1
!         rho_a = PRE*100.0/MAPL_RGAS/TEMP 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!         if  ((NIPRE .gt. 1.0) .and. (QICE0 .gt. 1.0e-10)) then 
!            DC=max((QICE0/(NIPRE*500.0*MAPL_PI))**(0.333), 40.0e-6) !Assumme monodisperse size distribution about size distribut
!ion is made. 
!         else   
!            DC = 40.0e-6
!         end if
!
!         AUX=  NIPRE*DENICE*MAPL_PI*DC*DC
!         TEFF = DENAIR*2.0*((K1*DIFF+K2)*DC+(K3/0.1))
!
!         if  (AUX .gt. 1.0e-12) then 
!            TEFF=min(TEFF/AUX, 1.0e20)   
!            DQold= SI/TEFF  
!         else
!            DQold=0.0
!         end if
!
!         ! Calculate rimming fraction
!         IWC =  QICE0*rho_a
!         LWC =  max(QT-QICE0, 0.0)*rho_a
!         aux = DQold ! only due to deposition
!
!
!         !Account for rimming
!
!         if ((LWC .gt. qmin_rim)  .and. (IWC .gt. qmin_rim))  then 
!            RIM = 6.0e-5/(LWC*(IWC**0.17)) !Fom Lin and Colle, NRW, 2011
!            RIM = 1.0/(1.0+RIM)
!            RIM  = min (0.95, RIM)
!            DQold =  DQold*(1 + RIM/(1.0-RIM))
!            FNDrim =  max(min(rho_a*(DQold -aux)/LWC, 1.0), 0.0) !Fraction of liquid condensate removed due to riming   
!         END if
!
!
!
!
!!!!!!!!!!!!!!!!!!!!!!recently nucleated!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!         AUX=  NINUC*DENICE*MAPL_PI*20.0e-6*20.0e-6
!         TEFF = DENAIR*2.0*((K1*DIFF+K2)*DC+(K3/0.1))
!
!         if  (AUX .gt. 1.0e-12)  then 
!            TEFF=min(TEFF/AUX, 1.0e10)
!            DQnew= SI/TEFF
!         else
!            DQnew = 0.0
!         end if
!
!         DQIG = DQold+DQnew     
!
!      end subroutine Qgrowth
!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!
!
!      !*************************************************************
!      ! Function PDG07 (simplified background ice nucleation 
!      !                     spectra according to Phillips et. al. 2007).  
!      ! si is supersaturation wrt ice and T is in K 
!      !************************************************************  
!
!      subroutine PDG07_ice(si, Tx, N)     
!
!         REAL(8), intent(in) :: si, Tx
!         REAL(8), intent(out)  :: N 
!         N=0.0
!
!         !if (Tx .le. 243.0)then
!
!
!         N=1000.0*exp(-0.388)*(exp(3.88*si)-1.0)/0.76
!         !elseif (Tx .le. 260.0) then
!         !  N=60.0*exp(-0.639)*(exp(12.96*si)-1.0)/0.76   
!         !end if      
!      end subroutine PDG07_ice
!
!
!      !*********************************************************************
!*********************************************************************
!*********************************************************************
!******************** Relaxed Arakawa-Schubert ***********************
!************************ Parameterization ***************************
!********************** SCALAR RAS-1 VERSION  ************************
!************************* 31 DECEMBER 1999 **************************
!*********************************************************************
!************************** Developed By *****************************
!*********************************************************************
!************************ Shrinivas Moorthi **************************
!******************************* and *********************************
!************************** Max J. Suarez ****************************
!*********************************************************************
!******************** Laboratory for Atmospheres *********************
!****************** NASA/GSFC, Greenbelt, MD 20771 *******************
!*********************************************************************
!*********************************************************************
!  Input:
!  ------
! 
!     K0      : Number of vertical levels (increasing downwards)
!
!     DT      : Time step in seconds
!
!     RASAL   : Array of dimension K-1 containing relaxation parameters
!               for cloud-types detraining at those levels
!
!     CPO     : Specific heat at constant pressure (J/kg/K)
!
!     ALHLO   : Latent Heat of condensation (J/kg)
!
!     ALHL1   : Latent Heat of condensation + fusion (J/kg)
!
!     GRAVO   : Acceleration due to gravity (m/s^2)
!
!     PLE     : 2D array of dimension (IDIM,K0+1) containing pressure
!               in hPa at the interfaces of K-layers from top of the 
!               atmosphere to the bottom  (mb)
!
!     PKE     : 2D array of dimension (IDIM,K0+1) containing (PRS/P00) **
!               RKAP.  i.e. Exner function at layer edges.
!
!     PKL     : 2D array of dimension (IDIM,K0) ) containing the
!               Exner function at the layers.
!
!     QSS     : 2D array of dimension (IDIM,K0  ) containing the
!               saturation specific humidity at the layers. (kg/kg)
!
!     DQS     : 2D array of dimension (IDIM,K0  ) containing
!               d(qss)/dt at the layers.  (1/K)
!   
!     CNV_FRACTION    : 1D array of dimension (IDIM) containing
!               fraction of grid cell considered to be convective
!   
!  Update:
!  -------
!
!     THO     : 2D array of dimension (IDIM,K0) containing potential
!               temperature (K)
!
!     QHO     : 2D array of dimension (IDIM,K0) containing specific
!               humidity (kg/kg)
!
!     UHO     : 2D array of dimension (IDIM,K0) containing u-wind (m/s)
!
!     VHO     : 2D array of dimension (IDIM,K0) containing v-wind (m/s)
!
!  Output:
!  -------
!!
!     CLW     : 2D array of dimension (IDIM,K0) containing the
!               detrained cloud liquid water.  (kg/m^2/s) 
!
!     FLX     : 2D array of dimension (IDIM,K0) containing the
!               cloud-base mass flux for each cloud type ordered by
!               detrainment level.   (kg/m^2/s) 
!
!     FLXD    : 2D array of dimension (IDIM,K0) containing the
!               detrained  mass flux for each cloud type ordered by
!               detrainment level.   (kg/m^2/s) 
!
!     FLXC    : 2D array of dimension (IDIM,K0+1) containing the
!               total cloud mass flux for all cloud types through
!               the top of each level. (e.g., FLXC(K)=SUM(FLX(ICMIN:K))
!               and  FLXD(L) = FLXC(L+1)-FLSD(L) )
!                (kg/m^2/s) 
!
!     PRECU   : 1D (IDIM) Locally-handled convective precip
!               Zero if older version of RAS-1 is used. Nonzero if
!               RAS-2 is used.
!
!   AEROPROPS, Structure containing aerosol propoerties (in)
!   CNV_NICE, CNV_DROP.  Flux of ice crystals and droplet number at det level (1/m^2/s) (out)
!   CNV_FICE: Ice fraction in the detrained condensate (out)
!************************************************************************
!  ARGUMENTS
    INTEGER, INTENT(IN) :: idim, irun, k0, icmin
!, QLO, QIO, CLAN
    REAL*8, DIMENSION(idim, k0), INTENT(INOUT) :: tho, qho, uho, vho
!,PHIE
    REAL*8, DIMENSION(idim, k0 + 1), INTENT(IN) :: ple, pke
!,PLO,PKO,PHIO
    REAL*8, DIMENSION(idim, k0), INTENT(IN) :: qss, dqs
    REAL*8, DIMENSION(idim), INTENT(IN) :: cnv_fraction, rasal2_2d
    REAL*8, DIMENSION(k0 + 1), INTENT(IN) :: sige
!, FLX
    REAL*8, DIMENSION(idim, k0), INTENT(OUT) :: clw, flxd
!      REAL(8), DIMENSION (IDIM,K0+1), INTENT(  OUT) ::  FLXC
    REAL*8, DIMENSION(idim, k0), INTENT(OUT) :: cnv_prc3
!, CNV_CVW, CNV_QC
    REAL*8, DIMENSION(idim, k0), INTENT(OUT) :: cnv_updfrc
!      REAL(8), DIMENSION (IDIM,K0  ), INTENT(  OUT) ::  ENTLAM
!      REAL(8), DIMENSION (IDIM,K0  ), INTENT(  OUT) ::  HHO, HSO
    REAL*8, INTENT(IN) :: dt, cpo, alhlo, gravo
    REAL*8, INTENT(IN) :: alhl1, tice
    INTEGER, INTENT(IN) :: itrcr
    INTEGER, DIMENSION(idim, 2), INTENT(IN) :: seedras
    INTEGER, DIMENSION(idim), INTENT(IN) :: kcbl
    INTEGER, DIMENSION(idim), INTENT(IN) :: doras
    REAL*8, DIMENSION(idim), INTENT(IN) :: tpert, qpert
    REAL*8, DIMENSION(idim), INTENT(IN) :: co_auto
    REAL*8, DIMENSION(idim) :: mxdiam
!      REAL(8), DIMENSION (IDIM     ), INTENT(  OUT) ::  PRECU
    REAL*8, DIMENSION(idim, k0), INTENT(IN) :: wgt0, wgt1
!     REAL(8), DIMENSION(:),          INTENT(IN   ) ::  RASPARAMS
    TYPE(RASPARAM_TYPE), INTENT(IN) :: rasparams
!      INTEGER, DIMENSION(IDIM,K0), INTENT(  OUT) ::  IRC
    REAL*8, OPTIONAL, INTENT(INOUT) :: xho(idim, k0, itrcr)
    REAL*8, OPTIONAL, INTENT(OUT) :: triedlev_diag(idim, k0)
    REAL*8, OPTIONAL, INTENT(OUT) :: disske(idim, k0)
! Fraction scavenged per km
    REAL*8, OPTIONAL, INTENT(IN) :: fscav(itrcr)
! = 0 (no scav), = 1 (full scav)
!  LOCALS
    REAL*8, DIMENSION(k0) :: poi_sv, qoi_sv, uoi_sv, voi_sv
    REAL*8, DIMENSION(k0) :: poi, qoi, uoi, voi, dqq, bet, gam, cll, &
&   tmpt
    REAL*8, DIMENSION(k0) :: poi_c, qoi_c
    REAL*8, DIMENSION(k0) :: prh, pri, ght, dpt, dpb, pki, dissk0, &
&   dissk1, clantnd
    REAL*8, DIMENSION(k0) :: tcu, qcu, ucu, vcu, cln, rns, pol, dm
    REAL*8, DIMENSION(k0) :: qst, ssl, rmf, rnn, rn1, rmfc, rmfp
    REAL*8, DIMENSION(k0) :: gms, eta, gmh, eht, gm1, hcc, rmfd
    REAL*8, DIMENSION(k0) :: hol, hst, qol, zol, hcld, cll0, cllx, clli&
&   , cllb
    REAL*8, DIMENSION(k0) :: wsp, lambdsv, bke, cvw, updfrc
    REAL*8, DIMENSION(k0) :: rasal, mtkwi, updfrp, bk2, bk3, dll0, dllx
    REAL*8, DIMENSION(itrcr) :: xht
    REAL*8, DIMENSION(k0, itrcr) :: xoi, xcu, xoi_sv
    REAL*8, DIMENSION(k0 + 1) :: prj, prs, qht, sht, zet, xyd, xyd0
!      INTEGER,  DIMENSION(K0-1) :: RC
    INTEGER :: k, my_pe
    REAL*8, DIMENSION(idim, k0) :: lambdsv2
    REAL*8 :: tx2, tx3, uht, vht, akm, acr, alm, tth, qqh, shtrg, wspbl&
&   , dqx
!, BKE
    REAL*8 :: wfn, tem, trg, trgexp, evp, wlq, qcc, mtkw_max
    REAL*8 :: shtrg_fac, sige_minhol, wfnog
    INTEGER :: i, ic, l, icl, itr, icl_c, n_dtl
    INTEGER :: ndtlexpon
    INTEGER, DIMENSION(:), ALLOCATABLE :: icl_v
!  RASE GLOBAL CONSTANTS
    REAL*8 :: grav, cp, alhl, cpbg, alhi, cpi, gravi, ddt, lbcp, obg, &
&   afc
!!!!!!!!!
    REAL*8 :: fricfac, dpth_bl, wupdrft, pblfrac, autorampb, co_zdep
    REAL*8 :: rasal1, rasal2, rasal2i, co_t, rasncl, friclambda, sdqvt1&
&   , sdqv2
    REAL*8 :: lambda_fac, strapping, acritfac, hmintrigger, lldisaggxp
    REAL*8 :: lambmx_fac, diammn_min, rdtlexpon, cli_crit, sdqv3, &
&   maxdallowed_d, maxdallowed_s
    REAL*8 :: rhmn, rhmx, cldmicro, fdrop_dust, fdrop_soot
    INTEGER :: kstrap, rasal_exp
    REAL*8 :: cld_radius, areal_frac, spect_mflx, cvw_cbase
!!!!!!!!!
    REAL*8, PARAMETER :: onepkap=1.+2./7., daylen=86400.0
!      REAL(8), PARAMETER :: PBLFRAC = 0.5
    REAL*8, PARAMETER :: rhmax=0.9999
!  LAMBDA LIMITS
    REAL*8 :: lambda_min
    REAL*8 :: lambda_max
!  TRIGGER PARAMETERS
! Density of liquid water in kg/m^3
    REAL*8, PARAMETER :: rho_w=1.0e3
    LOGICAL :: dyna_strapping, do_tracers, smooth_hst
    CHARACTER(len=esmf_maxstr) :: cbl_style
    REAL*8, DIMENSION(k0) :: tcu8, qcu8, pcu, flx8
    REAL*8, DIMENSION(k0, itrcr + 2) :: rcu
!, dpd, tla
    REAL*8 :: cup
    LOGICAL :: revap, wrkfun, calkpb, crtfun, lprnt, dndrft
    REAL*8, DIMENSION(k0) :: toi8, qoi8, prsm8, phil8, qli8, qii8, &
&   trcfac
    REAL*8, DIMENSION(k0) :: alfind, alfint, alfinq, rhc_ls
    REAL*8, DIMENSION(k0 + 1) :: prs8, phih8
    REAL*8, DIMENSION(k0, itrcr + 2) :: roi8
    REAL*8 :: fracbl, dt8, rasalf
    INTEGER :: kpbl
! no inhibition for =1.0
    REAL*8, SAVE :: max_neg_bouy=1.0
!!real*8 :: ALFINT = 0.5
!!real*8 :: ALFINQ = 0.5
! not used
    REAL*8, SAVE :: rhfacl=0.0
! no inhibition
    REAL*8, SAVE :: rhfacs=0.0
! 1 degree resolution
    REAL*8, SAVE :: garea=1.e10
!!real*8 :: ALFIND = 1.0
!!real*8 :: RHC_LS = 0.80
    REAL*8, SAVE :: dsfc=0.001
    REAL*8, SAVE :: cd=1.e-3
    REAL*8, SAVE :: wfnc=0.0
    REAL*8, SAVE :: tla=-1.0
    REAL*8, SAVE :: dpd=300.
!  SCAVANGING RELATED PARAMETERS
! layer thickness in km
    REAL*8 :: delzkm
! fraction of tracer *not* scavenged
    REAL*8 :: fnoscav
! Fraction scavenged per km
    REAL*8 :: fscav_(itrcr)
    INTRINSIC PRESENT
    INTRINSIC INT
    INTRINSIC SQRT
    INTRINSIC SUM
    INTRINSIC ALLOCATED
! ************************AER_CLOUD *********************************************
!      TYPE(AerProps),    DIMENSION(IDIM, K0),    INTENT(IN)    :: AEROPROPS !DONIF
!      REAL(8), DIMENSION (IDIM,K0), INTENT( OUT) ::  CNV_NDROP, CNV_FICE, CNV_NICE
!      REAL(8),  DIMENSION(K0) :: CNVNDROP, CNVNICE, CNVFICE !DONIF
!
!      TYPE(AerProps), DIMENSION(K0) ::  AERO   !AEROSOL VERTICAl ARRAY FOR ALL SPECIES    
!      REAL(8) ::  T_ICE_ALL, T_ICE_MAX, ASEASALT, f_seasalt 
!      INTEGER, PARAMETER :: NDUSTMAX = 10
!
!      INTEGER :: INDEX        
!      TYPE(AerProps) :: AERAUX, AER_BASE
!
!      CNV_FICE  =0.0
!      CNV_NDROP =0.0
!      CNV_NICE  =0.0
!      CNVFICE  =0.0
!      CNVNDROP =0.0
!      CNVNICE  =0.0
!      T_ICE_ALL= 238.0
!      T_ICE_MAX= MAPL_TICE
!      CLDMICRO = 0.0
!      FDROP_DUST = 0.1
!      FDROP_SOOT = 0.01
! *********************************************************************
    IF (PRESENT(fscav)) THEN
      fscav_ = fscav
    ELSE
! NO SCAVENGING BY DEFAULT
      fscav_ = 0.0
    END IF
    IF (irun .LE. 0) THEN
      RETURN
    ELSE
      IF (PRESENT(triedlev_diag)) triedlev_diag = 0.
      cnv_prc3 = 0.
      cnv_updfrc = 0.
!      CNV_CVW   =0. 
!      CNV_QC    =0.
!      ENTLAM    =0.
      IF (PRESENT(disske)) disske = 0.
!!LAMBDSV2 = 0.
!      IRC = -2
!      SMOOTH_HST   = .TRUE. 
      smooth_hst = .false.
      fricfac = rasparams%cufricfac
      shtrg_fac = rasparams%shr_lambda_fac
! MAT CO_AUTO is now passed in from outside
!     in order to allow this code to run over im*jm 
!     columns
!CO_AUTO      = RASPARAMS(3)     !  ---  3
      cli_crit = rasparams%qc_crit_cn
      rasal1 = rasparams%rasal1
      rasal2 = rasparams%rasal2
      rasncl = rasparams%rasncl
      lambda_fac = rasparams%lambda_fac
      lambmx_fac = rasparams%lambmx_fac
      diammn_min = rasparams%min_diameter
      friclambda = rasparams%cufriclambda
      rdtlexpon = rasparams%rdtlexpon
      strapping = rasparams%strapping
      sdqv2 = rasparams%sdqv2
      sdqv3 = rasparams%sdqv3
      sdqvt1 = rasparams%sdqvt1
      acritfac = rasparams%acritfac
      hmintrigger = rasparams%hmintrigger
      lldisaggxp = rasparams%lldisaggxp
      pblfrac = rasparams%pblfrac
      autorampb = rasparams%rasautorampb
      co_zdep = rasparams%autoc_cn_zdep
      maxdallowed_s = rasparams%maxdallowed_s
      maxdallowed_d = rasparams%maxdallowed_d
      rasal_exp = rasparams%rasal_exp
      rhmn = rasparams%ras_rhmin
      rhmx = rasparams%ras_rhfull
      cldmicro = rasparams%cldmicro
      fdrop_dust = rasparams%fdrop_dust
      fdrop_soot = rasparams%fdrop_soot
      IF (strapping .LE. 0.0) THEN
        dyna_strapping = .true.
      ELSE
        dyna_strapping = .false.
        kstrap = INT(strapping)
      END IF
      do_tracers = PRESENT(xho) .AND. itrcr .GT. 0
      wupdrft = 2.500
      grav = gravo
      alhl = alhlo
      cp = cpo
      cpi = 1.0/cp
      alhi = 1.0/alhl
      gravi = 1.0/grav
      cpbg = cp*gravi
      ddt = daylen/dt
      afc = -(1.04e-4*SQRT(dt*113.84))
      lbcp = alhl*cpi
      obg = 100.*gravi
!      HHO = MAPL_UNDEF
!      HSO = MAPL_UNDEF
!#ifdef RAS2
!      dt8    =  dt
!
!      call set_ras_afc(dt8)
!      call ras_init(k0,1)
!
!      kpbl  = k0
!      lprnt = .false.
!      trcfac(:) = 1.0
!      revap  = .true.
!      wrkfun = .false.
!      calkpb = .true.
!      crtfun = .true.
!      dndrft = .true.
!      wfnc   = 0.0
!
!#endif
      DO i=1,irun
        IF (doras(i) .NE. 0) THEN
!===================AER_CLOUD
!            AERO =  AEROPROPS(I, :)
!            CNVFICE  =0.0
!            CNVNDROP =0.0
!            CNVNICE  =0.0
!#ifndef RAS2
!!CALL FINDBASE
          k = kcbl(i)
!         rc(icmin) = 0
          CALL FINDDTLS()
          IF (k .GT. 0) THEN
            CALL STRAP(0)
            CALL HTEST()
!            HHO(I,:) = HOL
!            HSO(I,:) = HST
            DO icl_c=1,n_dtl
              icl = icl_v(icl_c)
              IF (do_tracers) xcu(icmin:, :) = 0.
              IF (PRESENT(triedlev_diag)) triedlev_diag(i, icl) = 1.
! This change makes cumulus friction
              ucu(icmin:) = 0.
! correct.
              vcu(icmin:) = 0.
              IF (icl .GT. icmin) CALL CLOUDE(icl)
            END DO
!               ENTLAM(I,ICL) = ALM
            IF (SUM(rmf(icmin:k)) .GT. 0.0) THEN
              CALL RNEVP()
              CALL STRAP(1)
            ELSE
              CALL STRAP(2)
            END IF
          ELSE
            CALL STRAP(2)
          END IF
        END IF
      END DO
!         PRECU = 0.0  ! Zero out precip - TBD w/ in progno_cloud
!#else
!
!
!         K = K0
!         CALL FINDDTLS
!
!         !===>  K      INPUT   THE RISE & THE INDEX OF THE SUBCLOUD LAYER
!         !===>  KD     INPUT   DETRAINMENT LEVEL ( 1<= KD < K )
!         !===>  M      INPUT   NUMBER OF TRACERS. MAY BE ZERO.
!         !===>  RASALF INPUT
!         !===>  FRACBL INPUT   MASS FLUX LIMIT
!         !===>  MAX_NEG_BOUY  INPUT  FRACTION OF NEG ALLOWED
!         !===>  ALFINT(K) INPUT   INTERFACE UPWIND (=1) PARAMETER FOR T
!         !===>  ALFINQ(K) INPUT   INTERFACE UPWIND (=1) PARAMETER FOR Q
!         !===>  RHFACL INPUT   CRITICAL RH AT PBL TOP (LAND)
!         !===>  RHFACS INPUT   CRITICAL RH AT PBL TOP (SEA) ???
!         !===>  GAREA  INPUT   AREA OF GRID BOX
!         !===>  ALFIND(K) INPUT   INTERFACE UPWIND (=1) PARAMETER FOR DOWNDRAFT
!         !===>  RHC_LS(K) INPUT   CRITICAL RH FOR RE-EVAP
!
!         ALFINT = 0.5
!         ALFINQ = 0.5
!         ALFIND = 0.5
!         RHC_LS = 0.8
!
!         tcu8  = 0.0
!         qcu8  = 0.0
!         rcu  = 0.0
!         pcu  = 0.0
!         flx8 = 0.0
!         cup  = 0.0
!
!         cloudloop: DO ICL = ICMIN,K0-1
!            IF (DO_TRACERS) XCU(ICMIN:,:) = 0.
!
!
!            rasalf = rasal(ICL)
!            FRACBL = pblfrac
!
!            toi8 = pko(i,:)*tho(i,:) ! tbdone
!            qoi8 = qho(i,:)
!            roi8(:,1) = uho(i,:)
!            roi8(:,2) = vho(i,:)
!            if(present(xho)) roi8(:,3:size(xho,3)+2) = xho(i,:,:) 
!
!            prs8  = ple(i,:)
!            prsm8 = plo(i,:) ! tbdone
!            phil8 = phio(i,:) ! tbdone
!            phih8 = phie(i,:) ! tbdone
!            qli8  = qlo(i,:) ! tbdone
!            qii8  = qio(i,:) ! tbdone
!
!            !===>  TOI(K)     INOUT   TEMPERATURE             KELVIN
!            !===>  QOI(K)     INOUT   SPECIFIC HUMIDITY       NON-DIMENSIONAL
!            !===>  ROI(K,M)   INOUT   TRACER                  ARBITRARY
!            !===>  PRS(K+1)   INPUT   PRESSURE @ EDGES        MB
!            !===>  PRSM(K)    INPUT   PRESSURE @ LAYERS       MB
!            !===>  PHIL(K)    INPUT   GEOPOTENTIAL @ LAYERS IN MKS units
!            !===>  PHIH(K+1)  INPUT   GEOPOTENTIAL @ EDGES  IN MKS units
!            !===>  QLI(K)     INOUT   LIQUID WATER            NON-DIMENSIONAL
!            !===>  QII(K)     INOUT   ICE                     NON-DIMENSIONAL
!            !===>  KPBL       INOUT   PBL TOP
!            !===>  DSFC       INOUT   GUSTINESS (m/s)
!            !===>  CD         INOUT   DRAG COEF. (NON-DIMENSIONAL)
!            !===>  LPRNT      INOUT
!            !===>  TRCFAC     INOUT   tracer pg factor for u AND V 
!            !===>  TCU(K  )   UPDATE  TEMPERATURE TENDENCY       DEG
!            !===>  QCU(K  )   UPDATE  WATER VAPOR TENDENCY       (G/G)
!            !===>  RCU(K,M)   UPDATE  TRACER TENDENCIES          ND
!            !===>  PCU(K-1)   UPDATE  PRECIP @ BASE OF LAYER     KG/M^2
!            !===>  FLX(K  )   UPDATE  UPDRAFT MASS FLUX @ TOP OF LAYER   KG/M^2
!            !===>  CUP        UPDATE  PRECIPITATION AT THE SURFACE KG/M^2
!            !===>  REVAP      UPDATE  LOGICAL TO CONTROL REEVAP IN ENVIRONMENT
!            !===>  DT         UPDATE  TIME STEP USED FOR GUSTINESS
!            !===>  WFNC       UPDATE  WORK FUNCTION (in if CRTFUN = F, out if WRKFUN = T ) 
!            !===>  WRKFUN     UPDATE  LOGICAL TO DO WRKFUNC ONLY
!            !===>  CALKBP     UPDATE  LOGICAL TO COMPUTE KPBL INTERNALLY
!            !===>  CRTFUN     UPDATE  LOGICAL TO CONTROL CRITICAL WORK FUNCTION
!            !===>  TLA        UPDATE  TILTING ANGLE FOR DD (NEG USE TABLE)
!            !===>  DNDRFT     INPUT   LOGICAL TO CONTROL DOWNDRAFT
!            !===>  DPD        INPUT   Minumum Cloud Depth for DOWNDRFAT Computation hPa
!
!            call CLOUD(                                                       &
!                  &                  K0, ICL, size(xcu,2)                            &
!                  &,                 RASALF, FRACBL, MAX_NEG_BOUY                    &
!                  &,                 ALFINT, ALFINQ, RHFACL, RHFACS, garea           &
!                  &,                 ALFIND, RHC_LS                                  &
!                  
!                  &,                 TOI8, QOI8, ROI8, PRS8, PRSM8, phil8, phih8     &
!                  &,                 QLI8, QII8, KPBL, DSFC                            &
!                  &,                 CD,lprnt, trcfac                                &
!                  
!                  &,                 TCU8, QCU8, RCU, PCU, FLX8                         &
!                  &,                 CUP, REVAP, DT8                                 &
!                  &,                 WFNC, WRKFUN, CALKPB, CRTFUN, TLA, DnDRFT, DPD)  
!
!            tho(i,:)  = toi8/pko(i,:)
!            qho(i,:)  = qoi8
!            uho(i,:)  = roi8(:,1) 
!            vho(i,:)  = roi8(:,2)
!            qlo(i,:)  = qli8 ! tbdone
!            qio(i,:)  = qii8  ! tbdone
!
!
!            if(present(xho)) xho(i,:,:) =roi8(:,3:size(xho,3)+2)
!         ENDDO cloudloop
!
!         ! temporary zeroing of outputs needed by rest of moist
!         !-----------------------------------------------------
!         cnv_prc3(i,:)   = 0.0
!         clw(i,:)        = 0.0
!         HHO(i,:)        = 0.0
!         HSO(i,:)        = 0.0
!         cnv_updfrc(i,:) = 0.0001
!         cnv_qc(i,:)     = 0.0
!         flxc(i,1:k0)    = flx8 / DT
!         flx(i,:)        = 0.0
!
!         flxd(i,1:k0-1)  = MAX( flxc(i,2:k0)- flxc(i,1:k0-1) , 0. )
!         flxd(i,k0)      = 0.0
!
!         precu(i)        = cup / DT   ! 
!
!         !
!         ! Tiedtke-style prognostic anvil fraction. 
!         ! Normally done in progno_cloud. Moved here
!         ! since RAS-2 manages its own detraining 
!         ! ice and liquid condensate.
!         !---------------------------------------------
!         dm           = 100.*(prs8(2:k0+1)-prs8(1:k0) )/grav
!         CLANTND      = FLXD(i,:) / DM
!         CLAN(i,:)    = CLAN(i,:) + CLANTND*DT
!         CLAN(i,:)    = MIN( CLAN(i,:) , 0.99 )
!
!         ! Now detraining mass flux has to be zeroed
!         ! or clouds will be crated twice
!         !--------------------------------
!         flxd(i,:)    = 0.
!
!#endif
      IF (ALLOCATED(icl_v)) THEN
        DEALLOCATE(icl_v)
      END IF
      RETURN
    END IF

  CONTAINS
!*********************************************************************
    SUBROUTINE CLOUDE(ic)
      IMPLICIT NONE
!=======================================
      INTEGER, INTENT(IN) :: ic
      REAL*8 :: deep_fact, cu_diam, wscale
!, dQx
      REAL*8 :: cli, te_a, c00_x, cli_crit_x, pete, toki, gmhx, hstx
      REAL*8 :: dt_lyr, rate, cvw_x, closs, f2, f3, f4, f5
      INTEGER :: k700
      INTRINSIC MIN
      INTRINSIC MAX
      INTRINSIC SQRT
      INTRINSIC EXP
      REAL*8 :: min2
      REAL*8 :: min1
      REAL*8 :: x5
      REAL*8 :: x4
      REAL*8 :: x3
      REAL*8 :: x2
      REAL*8 :: x1
      LOGICAL :: mask(k-ic+1)
      REAL*8 :: max2
      REAL*8 :: max1
      REAL*8 :: y1
!         !=============================AER_CLOUD local variables ====================
!         REAL(8) :: WBASE, NDROP, NICE, FP_D, FF_A, FP_I, FICE, &
!               NDROP_AMB, NSOOT_AMB, NSOOT, NIN, INSOOT, dCVW2, QICE, &
!               dQICE, dQIG, FPICE, dNICE, dNDROP, DSOOT_AMB, DSOOT, QLIQ, dQLIQ, FPRECIP, AUX, QT, &
!               MAXNICE, MAXNDROP, MINNICE, MINNDROP, NDROP_ACT, RIMM, FNDRIM, TminusTa, Tparcel , &
!	      alph_e, beta_e, RH_AMB, ECRIT
!
!         REAL(8), DIMENSION (NDUSTMAX) :: NDUST, NDUST_AMB, INDUST, DDUST_AMB, DDUST 
!         INTEGER :: INX
!
!
!         T_ICE_ALL = 238.0 !A little higher so there is ice at the freezing level
!         WBASE=1.0 
!         FICE=0.0
!         NICE=0.0
!         NDROP=0.0
!         NDROP_AMB=0.0
!         NDROP_ACT=0.0
!         NIN = 0.0
!         NDUST=0.0
!         NSOOT=0.0
!         NDUST_AMB=0.0
!         NSOOT_AMB = 0.0
!         dCVW2=0.0
!         QICE=0.0
!         QLIQ=0.0
!         FPICE = 0.0
!         INDUST=0.0
!         INSOOT= 0.0
!         QT= 0.0
!         FPRECIP=0.0   
!         FNDRIM = 0.0
!         RIMM = 0.0   
!         TminusTa = 0.0
!
!         f_seasalt = 0.0
!         aseasalt = 0.0
!
!         call init_Aer(AER_BASE)
!
!         !AER_CLOUD=============================
      alm = 0.
      IF (1. .GT. (qoi(k)/qst(k)-rhmn)/(rhmx-rhmn)) THEN
        trg = (qoi(k)/qst(k)-rhmn)/(rhmx-rhmn)
      ELSE
        trg = 1.
      END IF
      IF (0.0 .LT. (autorampb-sige(ic))/0.2) THEN
        y1 = (autorampb-sige(ic))/0.2
      ELSE
        y1 = 0.0
      END IF
      IF (1.0 .GT. y1) THEN
        f4 = y1
      ELSE
        f4 = 1.0
      END IF
! to 1 at SIG=AUTORAMPB-0.2
      IF (sige(ic) .GE. 0.5) THEN
        f5 = 1.0
      ELSE
        f5 = 1.0 - 2.*co_zdep*(0.5-sige(ic))
        IF (f5 .LT. 0.0) THEN
          f5 = 0.0
        ELSE
          f5 = f5
        END IF
      END IF
      IF (trg .LE. 1.0e-5) THEN
! TRIGGER  =========>>
!            RC(IC) = 7
        RETURN
      ELSE
!  RECOMPUTE SOUNDING UP TO DETRAINMENT LEVEL
        poi_c = poi
        qoi_c = qoi
        poi_c(k) = poi_c(k) + tpert(i)
        qoi_c(k) = qoi_c(k) + qpert(i)
        zet(k+1) = 0.
        sht(k+1) = cp*poi_c(k)*prj(k+1)
        DO l=k,ic,-1
          IF (qst(l)*rhmax .GT. qoi_c(l)) THEN
            qol(l) = qoi_c(l)
          ELSE
            qol(l) = qst(l)*rhmax
          END IF
          IF (0.000 .LT. qol(l)) THEN
            qol(l) = qol(l)
          ELSE
            qol(l) = 0.000
          END IF
          ssl(l) = cp*prj(l+1)*poi_c(l) + grav*zet(l+1)
          hol(l) = ssl(l) + qol(l)*alhl
          hst(l) = ssl(l) + qst(l)*alhl
          tem = poi_c(l)*(prj(l+1)-prj(l))*cpbg
          zet(l) = zet(l+1) + tem
          zol(l) = zet(l+1) + (prj(l+1)-prh(l))*poi_c(l)*cpbg
        END DO
        DO l=ic+1,k
          tem = (prj(l)-prh(l-1))/(prh(l)-prh(l-1))
          sht(l) = ssl(l-1) + tem*(ssl(l)-ssl(l-1))
          qht(l) = .5*(qol(l)+qol(l-1))
        END DO
! SMOOTH HSTAR W/ 1-2-1 Filter
        IF (smooth_hst) THEN
! save for later
          hstx = hst(ic)
          DO l=k-1,ic+1,-1
            hst(l) = 0.25*(hst(l+1)+hst(l-1)) + 0.5*hst(l)
          END DO
          DO l=ic,ic
            hst(l) = 0.5*hst(l+1) + 0.5*hst(l)
          END DO
        END IF
!  CALCULATE LAMBDA, ETA, AND WORKFUNCTION
        lambda_min = .2/mxdiam(i)
        lambda_max = .2/200.
!     LAMBDA_MIN = .2/(LAMBDA_FAC*DPTH_BL)
!     LAMBDA_MAX = .2/( MAX( LAMBMX_FAC*DPTH_BL , DIAMMN_MIN ) )
        IF (hol(k) .LE. hst(ic)) THEN
! CANNOT REACH IC LEVEL  ======>>
!            RC(IC) = 1
          RETURN
        ELSE
!  LAMBDA CALCULATION: MS-A18
          tem = (hst(ic)-hol(ic))*(zol(ic)-zet(ic+1))
          DO l=ic+1,k-1
            tem = tem + (hst(ic)-hol(l))*(zet(l)-zet(l+1))
          END DO
          IF (tem .LE. 0.0) THEN
! NO VALID LAMBDA  ============>>
!            RC(IC) = 2
            RETURN
          ELSE
            alm = (hol(k)-hst(ic))/tem
            IF (alm .GT. lambda_max) THEN
!            RC(IC) = 3
              RETURN
            ELSE
              toki = 1.0
              IF (alm .LT. lambda_min) toki = (alm/lambda_min)**2
!we can probably replace this by a actual distribution based on grid cell size
!RC(IC) = 6
!RETURN
!LAMBDSV(IC) = ALM
!  ETA CALCULATION: MS-A2
              DO l=ic+1,k
                eta(l) = 1.0 + alm*(zet(l)-zet(k))
              END DO
              eta(ic) = 1.0 + alm*(zol(ic)-zet(k))
!  WORKFUNCTION CALCULATION:  MS-A22
              wfn = 0.0
              hcc(k) = hol(k)
              DO l=k-1,ic+1,-1
                hcc(l) = hcc(l+1) + (eta(l)-eta(l+1))*hol(l)
                tem = hcc(l+1)*dpb(l) + hcc(l)*dpt(l)
                eht(l) = eta(l+1)*dpb(l) + eta(l)*dpt(l)
                wfn = wfn + (tem-eht(l)*hst(l))*gam(l)
              END DO
              hcc(ic) = hst(ic)*eta(ic)
              wfn = wfn + (hcc(ic+1)-hst(ic)*eta(ic+1))*gam(ic)*dpb(ic)
!  VERTICAL VELOCITY/KE CALCULATION (ADDED 12/2001 JTB)
              bk3(k) = 0.0
              bk2(k) = 0.0
              bke(k) = 0.0
              hcld(k) = hol(k)
              DO l=k-1,ic,-1
                hcld(l) = (eta(l+1)*hcld(l+1)+(eta(l)-eta(l+1))*hol(l))/&
&                 eta(l)
                tem = (hcld(l)-hst(l))*(zet(l)-zet(l+1))/(1.0+lbcp*dqq(l&
&                 ))
                IF (cldmicro .LE. 0.0) THEN
                  bke(l) = bke(l+1) + grav*tem/(cp*prj(l+1)*poi(l))
                  IF (tem .LT. 0.0) THEN
                    max1 = 0.0
                  ELSE
                    max1 = tem
                  END IF
                  bk2(l) = bk2(l+1) + grav*max1/(cp*prj(l+1)*poi(l))
                  IF (tem .GT. 0.0) THEN
                    min1 = 0.0
                  ELSE
                    min1 = tem
                  END IF
                  bk3(l) = bk3(l+1) + grav*min1/(cp*prj(l+1)*poi(l))
                  IF (bk2(l) .LT. 0.0) THEN
                    max2 = 0.0
                  ELSE
                    max2 = bk2(l)
                  END IF
                  cvw(l) = SQRT(2.0*max2)
                END IF
              END DO
! 1.0 / ( 5.0*ALM )
              cu_diam = 1000.
!   ALPHA CALCULATION 
              rasal2i = rasal2_2d(i)
              IF (zet(ic) .LT. 2000.) rasal(ic) = rasal1
              IF (zet(ic) .GE. 2000.) THEN
                IF (1.0 .GT. (zet(ic)-2000.)/8000.) THEN
                  min2 = (zet(ic)-2000.)/8000.
                ELSE
                  min2 = 1.0
                END IF
!WMP     RASAL(IC) = RASAL1 + (RASAL2i-RASAL1)*(ZET(IC) - 2000.)/8000.
                rasal(ic) = rasal1 + (rasal2i-rasal1)*min2**rasal_exp
              END IF
!WMP  RASAL(IC) = MIN( RASAL(IC) , 1.0e5 )
              rasal(ic) = dt/rasal(ic)
              mask(1:k-ic+1) = cvw(ic:k) .LT. 1.00
              WHERE (mask(1:k-ic+1)) 
                cvw(ic:k) = 1.00
              ELSEWHERE
                cvw(ic:k) = cvw(ic:k)
              END WHERE
!  NOTE THIS "CENTRALIZES" A KLUGE PRESENT IN OTHER LOCATIONS.
!  CLEAN UP SOME TIME.      -JTB 12/04/03
!  TEST FOR CRITICAL WORK FUNCTION
              CALL ACRITN(pol(ic), prs(k), acr)
              IF (wfn .LE. acr) THEN
! SUB-CRITICAL WORK FUNCTION ======>>
!            RC(IC) = 4
                RETURN
              ELSE
!  CLOUD TOP WATER AND MOMENTUM (TIMES ETA(IC)) MS-A16
! Tracer scavenging
! RAS loops over a series of plumes all having common cloud base level K 
! and different detrainment levels IC.  The plumes operate sequentially
! on the grid box mean quantities (wind, moisture, tracer) and so each
! subsequent plume is seeing the effects of previous plumes.  We parameterize
! scavenging following Liu et al. [JGR, 2001], their equation 1:
!  AEROSOL FRACTION SCAVENGED = 1 - exp(-FSCAV*DZ)
! where FSCAV is a specified scavenging efficiency [km-1] and DZ is the
! distance [km] the tracer traverses in the plume from it's entrainment
! level to its detrainment level.  We write the aerosol fraction surviving as:
!  FNOSCAV = exp(- FSCAV_(ITR) * DZ)
! The total scavenging is proportional to the convective mass flux, which
! is not explicitly solved for at this point.
                IF (do_tracers) THEN
                  DO itr=1,itrcr
!           Scavenging of the below cloud tracer
                    delzkm = (zet(ic)-zet(k))/1000.
                    x4 = EXP(-(fscav_(itr)*delzkm))
                    IF (x4 .GT. 1.) THEN
                      x1 = 1.
                    ELSE
                      x1 = x4
                    END IF
                    IF (x1 .LT. 0.) THEN
                      fnoscav = 0.
                    ELSE
                      fnoscav = x1
                    END IF
                    xht(itr) = xoi(k, itr)*fnoscav
                  END DO
                END IF
                wlq = qol(k)
                uht = uoi(k)
                vht = voi(k)
                rnn(k) = 0.
                cll0(k) = 0.
!print *, '========================================='
                DO l=k-1,ic,-1
                  tem = eta(l) - eta(l+1)
                  wlq = wlq + tem*qol(l)
                  uht = uht + tem*uoi(l)
                  vht = vht + tem*voi(l)
                  IF (do_tracers) THEN
                    DO itr=1,itrcr
!         Scavenging of the entrained tracer.  Updates transported tracer mass.
                      delzkm = (zet(ic)-zet(l+1))/1000.
                      x5 = EXP(-(fscav_(itr)*delzkm))
                      IF (x5 .GT. 1.) THEN
                        x2 = 1.
                      ELSE
                        x2 = x5
                      END IF
                      IF (x2 .LT. 0.) THEN
                        fnoscav = 0.
                      ELSE
                        fnoscav = x2
                      END IF
                      xht(itr) = xht(itr) + tem*xoi(l, itr)*fnoscav
                    END DO
                  END IF
!!!! How much condensate (CLI) is present here? 
                  IF (l .GT. ic) THEN
                    tx2 = 0.5*(qst(l)+qst(l-1))*eta(l)
                    tx3 = 0.5*(hst(l)+hst(l-1))*eta(l)
                    qcc = tx2 + gm1(l)*(hcc(l)-tx3)
                    cll0(l) = wlq - qcc
                  ELSE
                    cll0(l) = wlq - qst(ic)*eta(ic)
                  END IF
                  IF (cll0(l) .LT. 0.00) THEN
                    cll0(l) = 0.00
                  ELSE
                    cll0(l) = cll0(l)
                  END IF
! condensate (kg/kg)
                  cli = cll0(l)/eta(l)
! Temperature (K)
                  te_a = poi(l)*prh(l)
!=====================================================================
                  IF (cldmicro .LE. 0.0) THEN
!AER_CLOUD MIcrophysics considering activation and nucleation 
!               !recompute vertical velocity
!
!               Tparcel = TE_A
!               CVW(K) = 0.8  ! Assume a below cloud base  W of 0.8 m s-1            
!               BK2(K)   = 0.0
!
!     
!               TEM     = (HCLD(L)-HST(L) )/ (1.0+LBCP*DQQ(L))  
!               TminusTa = max(min(TEM/CP, 5.0), 0.0) !limit DT to 5 K. According to Wei, JAS, 1998   
!	     TEM =0.33*TminusTa*CO_AUTO(I)/TE_A !Bouyancy term, effciency =0.5 mwr Roode et al    	     
!
!               BK2(L)  = BK2(L+1) + GRAV * TEM*(ZET(L)-ZET(L+1)) 
!               BK2(L) = BK2(L) - (ZET(L)-ZET(L+1))*(BK2(L+1)*ALM + CLI*GRAV)  !Account for drag from entrainment of stagnat air 
!and condesate loading
!               CVW(L) = max(SQRT(  2.0* MAX( BK2(L) , 0.0 )  ), 1.0) 
!
!
!	    CVW_X = MIN(CVW(L), 50.0)
!               DT_LYR  =  max(( ZET(L)-ZET(L+1) )/CVW_X, 1.0) !Sanity check 
!               TEM   = ETA(L) - ETA(L+1)
!
!               Tparcel  =  TE_A + TminusTa
!
!
!
!!!!!!!!!!account for entrainment effects on activation !!!!!!!!!!!
!               ! Barahona and Nenes, JGR, 2007
!               alph_e = 2.8915e-8*Tparcel*Tparcel -2.1328e-5*Tparcel+4.2523e-3
!               beta_e = MAPL_ALHL*TminusTa/MAPL_RVAP/Tparcel/Tparcel
!               RH_AMB=QOI(L)/QST(L)
!               ECRIT  = max(RH_AMB -beta_e, 1.0e-6) 
!               ECRIT =  alph_e/ECRIT
!               ! print *, L, Tparcel, RH_AMB, ECRIT, ALM
!	           ECRIT =  ALM/ECRIT
!       ! ECRIT = 0.0 ! do not use this for now
!               !Print *, ECRIT
!
!
!               if (L .eq. K-1) then
!
!                  FICE=0.0
!                  NICE=0.0
!                  NDROP=0.0
!                  NIN =0.0
!                  NDUST_AMB =0.0
!                  NSOOT_AMB = 0.0
!                  NSOOT=0.0
!                  NDUST= 0.0
!
!        
!                  AER_BASE =  AERO(L)
!                  RATE=0.0
!                  FPRECIP=0.0
!                  !initial conditions     
!                  call ARGact(Tparcel, CVW_X, NDROP_ACT, NDROP_AMB, NDUST_AMB, NSOOT_AMB, L,  .true., DDUST_AMB, DSOOT_AMB, ECRI
!T) !cloud droplet number and INsource at cloud base 
!                  NDUST=NDUST_AMB
!                  NSOOT=NSOOT_AMB
!                  DDUST=DDUST_AMB
!                  DSOOT=DSOOT_AMB                                     
!
!               else 
!                  call ARGact(Tparcel, CVW_X, NDROP_ACT, NDROP_AMB, NDUST_AMB, NSOOT_AMB, L, .false., DDUST_AMB, DSOOT_AMB, ECRI
!T) !cloud droplet number above cloud base  
!
!               end if
!
!               QT = CLI
!               RATE = 0.0
!               FPRECIP = 0.0
!
!               if (QT .gt. 0.0) then
!
!                  ! if FICE is already >= 1.0 then the cloud is glaciated and there is no need to do anymore partitioning
!
!                  if (FICE .ge. 1.0) then
!
!
!                     CALL  Qremoval(RATE, FICE, FP_D, FP_I, Tparcel,  & 
!                           POL(L), QT,  NICE, NDROP, CVW_X, FPICE, &
!                           DT_LYR, RIMM, CO_AUTO(I)) 
!
!
!
!                     dNICE = -NICE*FP_I 
!                     NICE  =  (NICE +dNICE*DT_LYR)*ETA(L+1)/ETA(L) !ice
!
!                     MINNICE = max(QICE*(1.0-RATE*DT_LYR), 0.0)/4.0e-8!assuming maximum vol radius 250 microns
!                     MAXNICE = max(QICE*(1.0-RATE*DT_LYR), 0.0)/2.51e-12 !assuming minimum vol radius 10 microns
!
!                     NICE=MIN(max(NICE, MINNICE), MAXNICE)
!
!                     FICE = 1.0
!
!                  else 
!
!                     ! Cloud is not completely glaciated do the whole thing
!                     ! ALL this subroutines return tendencies
!
!
!                     CALL  INfreezing(QLIQ, NDROP, NIN, NDUST, NSOOT, INDUST, INSOOT, Tparcel, POL(L), CVW_X, DDUST, DSOOT)  !ca
!lculate the freezing fraction of the aerosol at this level
!
!                     NIN = min(NIN, NDROP/DT_LYR)
!
!                     call Qgrowth(Tparcel, POL(L), QICE, NICE, QT, NIN, dQIG, RIMM, FNDRIM)
!
!                     CALL  Qremoval(RATE, FICE, FP_D, FP_I, Tparcel,  & 
!                           POL(L), QT,  NICE, NDROP, CVW_X, FPICE, &
!                           DT_LYR,  RIMM, CO_AUTO(I)) 
!
!
!
!                     !ice number tendency: -precip + freezin
!                     dNICE = -NICE*FP_I  + NIN     
!                     NICE  =  (NICE +dNICE*DT_LYR)*ETA(L+1)/ETA(L) !ice
!                     NICE =max(NICE, 0.0)
!
!
!                     !ice mass tendency: growth - precip
!                     dQICE = -QICE*FPICE + dQIG
!                     QICE  =  min((QICE + dQICE*DT_LYR)*ETA(L+1)/ETA(L), QT) !ice
!                     QICE=max(min(QT, QICE), 0.0)
!
!
!                     ! Liquid Tendency: source/evap -  precip 
!                     !dQLIQ = max((CLI-QICE), -QLIQ)/DT_LYR -QLIQ*max(RATE-FPICE, 0.0) 
!                     ! dQLIQ = CLI*(1.0-RATE*DT_LYR)/DT_LYR -dQICE - QLIQ*max(RATE-FPICE, 0.0)           
!                     !QLIQ  =  max((QLIQ + dQLIQ*DT_LYR)*ETA(L+1)/ETA(L), 0.0) !liquid. This is actually diagnostic
!                     QLIQ=max((QT-QICE), 0.0)
!
!
!                     !droplet number tendency: -precip - freezin + activation + activated entrained aerosol 
!
!
!                     dNDROP =-NDROP*FP_D - NIN -  FNDRIM*NDROP/DT_LYR + max(NDROP_ACT-NDROP, 0.0)/DT_LYR          
!
!                     !dNDROP =-NDROP*FP_D - NIN -  FNDRIM*NDROP/DT_LYR + NDROP_ACT/DT_LYR         
!
!                     NDROP =  (NDROP + dNDROP*DT_LYR)*ETA(L+1)/ETA(L) + &
!                           (ZET(L) - ZET(L+1))*ALM*MAX((NDROP_AMB-NDROP), 0.0)
!
!                     !Aerosol tendency: Entrainment - freezing 
!
!                     NDUST = (NDUST - INDUST*DT_LYR)*ETA(L+1)/ETA(L) + &
!                           (ZET(L) - ZET(L+1))*ALM*MAX(NDUST_AMB-NDUST, 0.0) 
!
!                     NSOOT =  (NSOOT - INSOOT*DT_LYR)*ETA(L+1)/ETA(L)  + &     
!                           (ZET(L) - ZET(L+1))*ALM*MAX(NSOOT_AMB-NSOOT, 0.0)  
!
!
!                           
!                     !Update FICE and perform Sanity checks
!
!
!                     MINNDROP = (1.0-FICE)*QLIQ*max(1.0-RATE*DT_LYR, 0.0)/2.e-10    !assuming maximum vol radius 36 microns
!                     MAXNDROP = (1.0-FICE)*QLIQ*max(1.0-RATE*DT_LYR, 0.0)/3.35e-14 !assuming minimum vol radius 2 microns
!                     MINNICE = QICE/4.0e-8!assuming maximum vol radius 250 microns
!                     MAXNICE = QICE/2.51e-12 !assuming minimum vol radius 10 microns
!
!                     IF ((NICE .gt. MAXNICE) .or. (NICE .lt. MINNICE))   then    
!                        !print *, 'nilim', NICE*1e-6, MINNICE*1e-6, MAXNICE*1e-6
!                     END IF
!
!                     IF ((NDROP .gt. MAXNDROP) .or. (NDROP .lt. MINNDROP))      then 
!                        !print *, 'ndroplim', NDROP*1e-6, MINNDROP*1e-6, MAXNDROP*1e-6
!                     end if
!
!
!                     NSOOT=MAX(NSOOT, 0.0)
!                     NDUST=MAX(NDUST, 0.0)              
!
!                     NDROP=MIN(max(NDROP, MINNDROP), MAXNDROP)
!                     NICE=MIN(max(NICE, MINNICE), MAXNICE)
!
!                     FICE=max(min(QICE/QT, 1.0), 0.0)
!
!                     IF (FICE .ge. 1.0) THEN !Complete glaciation 
!                        NICE=NICE+NDROP 
!                        NDROP = 0.0
!                        QICE  = QT
!                        QLIQ= 0.0
!                     END IF
!
!                     IF (Tparcel .LT. T_ICE_ALL) THEN !instantaneous freezing
!                        NICE=NICE+NDROP 
!                        NDROP = 0.0
!                        FICE  = 1.0
!                        QICE  = QT
!                        QLIQ=0.0
!                     END IF
!
!                     IF (Tparcel .GT. T_ICE_MAX) THEN !instantaneous melting
!                        NDROP=NICE+NDROP 
!                        NICE = 0.0
!                        FICE  = 0.0
!                        QICE  = 0.0
!                        QLIQ=QT
!                     END IF
!
!                  END IF
!
!               else 
!
!                  FICE =0.0 
!                  QICE = 0.0
!                  QLIQ = 0.0
!                  NICE= 0.0 
!                  NDROP = 0.0
!                  RATE =0.0
!               end if
!
!               FPRECIP= RATE*DT_LYR
!
!               !RATE=RATE*F4
!               ! NDROP=NDROP*F4
!               !NICE=NICE*(1.0-F4)
!
!               !print *, TE_A, FICE, 'NICE', NICE*1e-6, 'NDROP', NDROP*1e-6, L 
!               !print *, 'FPI', FP_I*DT_LYR, 'FPD', FP_D*DT_LYR, 'FPICE', FPICE, 'FPRE', FPRECIP, QT, QLIQ
!
!Bacmeister 2006 microphysics
                    CALL SUNDQ3_ICE(te_a, sdqv2, sdqv3, sdqvt1, f2, f3)
! * F5  ! F4 reduces AUTO for shallow clouds, F5 modifies auto for deep clouds
                    c00_x = co_auto(i)*f2*f3*f4
                    cli_crit_x = cli_crit/(f2*f3)
                    rate = c00_x*(1.0-EXP(-(cli**2/cli_crit_x**2)))
                  END IF
                  IF (cvw(l) .LT. 1.00) THEN
                    cvw_x = 1.00
                  ELSE
                    cvw_x = cvw(l)
                  END IF
! really trust it at low values
! l.h.s. DT_LYR => time in layer (L,L+1)
                  dt_lyr = (zet(l)-zet(l+1))/cvw_x
                  closs = cll0(l)*rate*dt_lyr
                  IF (closs .GT. cll0(l)) THEN
                    closs = cll0(l)
                  ELSE
                    closs = closs
                  END IF
                  cll0(l) = cll0(l) - closs
                  dll0(l) = closs
                  IF (closs .GT. 0.) THEN
                    wlq = wlq - closs
                    rnn(l) = closs
                  ELSE
                    rnn(l) = 0.
                  END IF
                END DO
!AER_CLOUD=======================================
!            CNVNDROP(IC)=NDROP
!            CNVNICE(IC)=NICE
!            CNVFICE(IC)=FICE
                wlq = wlq - qst(ic)*eta(ic)
!     CALCULATE GAMMAS AND KERNEL
! MS-A30 (W/O GRAV)
                gms(k) = (sht(k)-ssl(k))*pri(k)
! MS-A31 (W/O GRAV)
                gmh(k) = gms(k) + (qht(k)-qol(k))*pri(k)*alhl
! MS-A37 (W/O GRAV)
                akm = gmh(k)*gam(k-1)*dpb(k-1)
                tx2 = gmh(k)
                DO l=k-1,ic+1,-1
                  gms(l) = (eta(l)*(sht(l)-ssl(l))+eta(l+1)*(ssl(l)-sht(&
&                   l+1)))*pri(l)
                  gmh(l) = gms(l) + (eta(l)*(qht(l)-qol(l))+eta(l+1)*(&
&                   qol(l)-qht(l+1)))*alhl*pri(l)
                  tx2 = tx2 + (eta(l)-eta(l+1))*gmh(l)
                  akm = akm - gms(l)*eht(l)*pki(l) + tx2*ght(l)
                END DO
                gms(ic) = eta(ic+1)*(ssl(ic)-sht(ic+1))*pri(ic)
                akm = akm - gms(ic)*eta(ic+1)*dpb(ic)*pki(ic)
                gmh(ic) = gms(ic) + (eta(ic+1)*(qol(ic)-qht(ic+1))*alhl+&
&                 eta(ic)*(hst(ic)-hol(ic)))*pri(ic)
                IF (smooth_hst) gmhx = gms(ic) + (eta(ic+1)*(qol(ic)-qht&
&                   (ic+1))*alhl+eta(ic)*(hstx-hol(ic)))*pri(ic)
!    CLOUD BASE MASS FLUX
                IF (akm .GE. 0.0 .OR. wlq .LT. 0.0) THEN
!  =========>
!            RC(IC) = 5
                  RETURN
                ELSE
! MS-A39 MASS-FLUX IN Pa/step
                  wfn = -((wfn-acr)/akm)
                  x3 = rasal(ic)*trg*toki*wfn
                  IF (x3 .GT. (prs(k+1)-prs(k))*(100.*pblfrac)) THEN
                    wfn = (prs(k+1)-prs(k))*(100.*pblfrac)
                  ELSE
                    wfn = x3
                  END IF
!    CUMULATIVE PRECIP AND CLOUD-BASE MASS FLUX FOR OUTPUT
                  wfnog = wfn*gravi
                  tem = wfn*gravi
! (kg/m^2/step)
                  cll(ic) = cll(ic) + wlq*tem
! (kg/m^2/step)
                  rmf(ic) = rmf(ic) + tem
! (kg/m^2/step)
                  rmfd(ic) = rmfd(ic) + tem*eta(ic)
                  DO l=ic+1,k
! (kg/m^2/step)
                    rmfp(l) = tem*eta(l)
! (kg/m^2/step)
                    rmfc(l) = rmfc(l) + rmfp(l)
                    dllx(l) = dllx(l) + tem*dll0(l)
                    IF (cvw(l) .GT. 0.0) THEN
                      updfrp(l) = rmfp(l)*(ddt/daylen)*1000./(cvw(l)*prs&
&                       (l))
                    ELSE
                      updfrp(l) = 0.0
                    END IF
! current cloud; incloud condensate        
                    clli(l) = cll0(l)/eta(l)
!  cumulative grid mean convective condensate        
                    cllb(l) = cllb(l) + updfrp(l)*clli(l)
                    updfrc(l) = updfrc(l) + updfrp(l)
                  END DO
!    THETA AND Q CHANGE DUE TO CLOUD TYPE IC
                  DO l=ic,k
! (kg/m^2/step)
                    rns(l) = rns(l) + rnn(l)*tem
                    gmh(l) = gmh(l)*wfn
                    gms(l) = gms(l)*wfn
                    qoi(l) = qoi(l) + (gmh(l)-gms(l))*alhi
                    poi(l) = poi(l) + gms(l)*pki(l)*cpi
                    qst(l) = qst(l) + gms(l)*bet(l)*cpi
                  END DO
                  IF (smooth_hst) THEN
                    gmhx = gmhx*wfn
                    dqx = (gmhx-gmh(ic))*alhi
                    rns(ic) = rns(ic) + dqx/(pri(ic)*grav)
                  END IF
                  IF (do_tracers) THEN
!*FRICFAC*0.5
                    wfn = wfn*0.5*1.0
                    tem = wfn*pri(k)
                    DO itr=1,itrcr
                      xcu(k, itr) = xcu(k, itr) + tem*(xoi(k-1, itr)-xoi&
&                       (k, itr))
                    END DO
                    DO itr=1,itrcr
                      DO l=k-1,ic+1,-1
                        tem = wfn*pri(l)
                        xcu(l, itr) = xcu(l, itr) + tem*((xoi(l-1, itr)-&
&                         xoi(l, itr))*eta(l)+(xoi(l, itr)-xoi(l+1, itr)&
&                         )*eta(l+1))
                      END DO
                    END DO
                    tem = wfn*pri(ic)
                    DO itr=1,itrcr
                      xcu(ic, itr) = xcu(ic, itr) + (2.*(xht(itr)-xoi(ic&
&                       , itr)*(eta(ic)-eta(ic+1)))-(xoi(ic, itr)+xoi(ic&
&                       +1, itr))*eta(ic+1))*tem
                    END DO
                    DO itr=1,itrcr
                      DO l=ic,k
                        xoi(l, itr) = xoi(l, itr) + xcu(l, itr)
                      END DO
                    END DO
                  ELSE
!*FRICFAC*0.5
                    wfn = wfn*0.5*1.0
                  END IF
                  lambdsv(ic) = 1.000
!   CUMULUS FRICTION
                  IF (fricfac .LE. 0.0) THEN
!            RC(IC) = 0
!  NO CUMULUS FRICTION =========>>
                    RETURN
                  ELSE
                    wfn = wfn*fricfac*EXP(-(alm/friclambda))
                    tem = wfn*pri(k)
                    ucu(k) = ucu(k) + tem*(uoi(k-1)-uoi(k))
                    vcu(k) = vcu(k) + tem*(voi(k-1)-voi(k))
                    DO l=k-1,ic+1,-1
                      tem = wfn*pri(l)
                      ucu(l) = ucu(l) + tem*((uoi(l-1)-uoi(l))*eta(l)+(&
&                       uoi(l)-uoi(l+1))*eta(l+1))
                      vcu(l) = vcu(l) + tem*((voi(l-1)-voi(l))*eta(l)+(&
&                       voi(l)-voi(l+1))*eta(l+1))
                    END DO
                    tem = wfn*pri(ic)
                    ucu(ic) = ucu(ic) + (2.*(uht-uoi(ic)*(eta(ic)-eta(ic&
&                     +1)))-(uoi(ic)+uoi(ic+1))*eta(ic+1))*tem
                    vcu(ic) = vcu(ic) + (2.*(vht-voi(ic)*(eta(ic)-eta(ic&
&                     +1)))-(voi(ic)+voi(ic+1))*eta(ic+1))*tem
                    dissk0(ic) = eta(ic)*grav*wfnog*pri(ic)*0.5*((uht/&
&                     eta(ic)-uoi(ic))**2+(vht/eta(ic)-voi(ic))**2)
                    DO l=ic,k
                      uoi(l) = uoi(l) + ucu(l)
                      voi(l) = voi(l) + vcu(l)
                    END DO
!         RC(IC) = 0
                    RETURN
                  END IF
                END IF
              END IF
            END IF
          END IF
        END IF
      END IF
    END SUBROUTINE CLOUDE
    SUBROUTINE ACRITN(pl, plb, acr)
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: pl, plb
      REAL*8, INTENT(OUT) :: acr
      INTEGER :: iwk
!!REAL(8), PARAMETER :: FACM=0.5
      REAL*8, PARAMETER :: ph(15)=(/150.0, 200.0, 250.0, 300.0, 350.0, &
&       400.0, 450.0, 500.0, 550.0, 600.0, 650.0, 700.0, 750.0, 800.0, &
&       850.0/)
!!*FACM
      REAL*8, PARAMETER :: a(15)=(/1.6851, 1.1686, 0.7663, 0.5255, &
&       0.4100, 0.3677, 0.3151, 0.2216, 0.1521, 0.1082, 0.0750, 0.0664, &
&       0.0553, 0.0445, 0.0633/)
      INTRINSIC INT
      iwk = INT(pl*0.02 - 0.999999999)
      IF (iwk .GT. 1 .AND. iwk .LE. 15) THEN
        acr = a(iwk-1) + (pl-ph(iwk-1))*.02*(a(iwk)-a(iwk-1))
      ELSE IF (iwk .GT. 15) THEN
        acr = a(15)
      ELSE
        acr = a(1)
      END IF
      acr = acritfac*acr*(plb-pl)
      RETURN
    END SUBROUTINE ACRITN
    SUBROUTINE RNEVP()
      IMPLICIT NONE
      zet(k+1) = 0
      DO l=k,icmin,-1
        tem = poi(l)*(prj(l+1)-prj(l))*cpbg
        zet(l) = zet(l+1) + tem
      END DO
      DO l=icmin,k
        tem = pri(l)*grav
        cnv_prc3(i, l) = rns(l)*tem
      END DO
!! If hst is smoothed then adjusted precips may be negative
      IF (smooth_hst) THEN
        DO l=icmin,k
          IF (cnv_prc3(i, l) .LT. 0.) THEN
            qoi(l) = qoi(l) + cnv_prc3(i, l)
            poi(l) = poi(l) - cnv_prc3(i, l)*(alhl/cp)/prj(l+1)
            cnv_prc3(i, l) = 0.
          END IF
        END DO
      END IF
      RETURN
    END SUBROUTINE RNEVP
    SUBROUTINE HTEST()
      IMPLICIT NONE
      REAL*8, DIMENSION(k0) :: hol1
      INTEGER :: lminhol
      REAL*8 :: minhol
      INTRINSIC MIN
      INTRINSIC MAX
! HOL initialized here in order not to confuse Valgrind debugger
      hol = 0.
      lminhol = k + 1
      minhol = -999999.
      zet(k+1) = 0
      sht(k+1) = cp*poi(k)*prj(k+1)
      DO l=k,icmin,-1
        IF (qst(l)*rhmax .GT. qoi(l)) THEN
          qol(l) = qoi(l)
        ELSE
          qol(l) = qst(l)*rhmax
        END IF
        IF (0.000 .LT. qol(l)) THEN
          qol(l) = qol(l)
        ELSE
          qol(l) = 0.000
        END IF
        ssl(l) = cp*prj(l+1)*poi(l) + grav*zet(l+1)
        hol(l) = ssl(l) + qol(l)*alhl
        hst(l) = ssl(l) + qst(l)*alhl
        tem = poi(l)*(prj(l+1)-prj(l))*cpbg
        zet(l) = zet(l+1) + tem
        zol(l) = zet(l+1) + (prj(l+1)-prh(l))*poi(l)*cpbg
      END DO
      hol1 = hol
      DO l=k-1,icmin+1,-1
        hol1(l) = 0.25*hol(l+1) + 0.50*hol(l) + 0.25*hol(l-1)
        IF (minhol .GE. hol1(l) .OR. minhol .LT. 0.) THEN
          minhol = hol1(l)
          lminhol = l
        END IF
      END DO
      sige_minhol = sige(lminhol)
    END SUBROUTINE HTEST
    SUBROUTINE FINDDTLS()
      IMPLICIT NONE
!! NO SHALLOW CONV           ICL_V(L) = 56 - L
!         else
!            do L=1,N_DTL
!               call random_number ( SIGDT0 )
!               SIGDT0 = 1.00 - ( SIGDT0**RDTLEXPON )
!               SIGDT0 = SIGMIN+SIGDT0*(SIGMAX-SIGMIN)
!
!               do LL=ICMIN,K
!                  if ( (SIGE(LL+1)>=SIGDT0) .and. (SIGE(LL)<SIGDT0 ) ) ICL_V(L)=LL
!               enddo
!            end do
!         endif
      REAL*8 :: sigdt0, sigmax, sigmin
      INTEGER :: ll
      INTRINSIC ALLOCATED
!#ifndef __GFORTRAN__
!         integer :: THE_SEED(2)
!#else
!         integer :: THE_SEED(12)
!#endif
!         THE_SEED(1)=SEEDRAS(I,1)*IRAS(I) + SEEDRAS(I,2)*JRAS(I)
!         THE_SEED(2)=SEEDRAS(I,1)*JRAS(I) + SEEDRAS(I,2)*IRAS(I)
!         THE_SEED(1)=THE_SEED(1)*SEEDRAS(I,1)/( SEEDRAS(I,2) + 10)
!         THE_SEED(2)=THE_SEED(2)*SEEDRAS(I,1)/( SEEDRAS(I,2) + 10)
!         if(THE_SEED(1) == 0) THE_SEED(1) =  5
!         if(THE_SEED(2) == 0) THE_SEED(2) = -5
!#ifdef __GFORTRAN__
!         THE_SEED(3:12) = 0
!#endif
!
!         call random_seed(PUT=THE_SEED)
!dh RASNCL = -300
      sigmax = sige(k)
      sigmin = sige(icmin)
      IF (rasncl .LT. 0.0) n_dtl = k - icmin
!! NO SHALLOW CONV   N_DTL = 56 - ICMIN 
!            N_DTL = min( int( RASNCL ) , K-ICMIN )
      IF (ALLOCATED(icl_v)) THEN
        DEALLOCATE(icl_v)
      END IF
      ALLOCATE(icl_v(n_dtl))
!         if ( ( RASNCL < 0.0 ) .and. ( RASNCL >=-100.) ) then 
!            do L=1,N_DTL
!               ICL_V(L) = ICMIN + L - 1
!            enddo
!         else if ( RASNCL < -100.0 ) then 
      DO l=1,n_dtl
        icl_v(l) = k - l
      END DO
    END SUBROUTINE FINDDTLS
    SUBROUTINE STRAP(final)
      IMPLICIT NONE
!            CNV_CVW   (I,ICMIN:K-1)   =     CVW(ICMIN:K-1)
!            CNV_QC(I,ICMIN:K-1)       =  CLLB(ICMIN:K-1)
      INTEGER :: final
      REAL*8, DIMENSION(k0) :: wght, massf
      REAL*8 :: wght0, prcbl
      INTEGER, PARAMETER :: nrands=1
      REAL*8 :: rndu(nrands)
      INTEGER :: seedcbl(nrands), jj
! !DESCRIPTION: 
!   {\tt STRAP} is called: FINAL=0, to compute cloud base layer CBL properties
!   given a value K for the index of the upper {\em EDGE} of the CBL; FINAL=1
!   to redistribute convective tendencies within CBL
      INTEGER :: kk
      INTRINSIC SQRT
      INTRINSIC MAX
      INTRINSIC ABS
      INTRINSIC PRESENT
      INTRINSIC ALLOCATED
      REAL*8 :: abs1
      REAL*8 :: abs0
!  LOCAL VARIABLES FOR USE IN CLOUDE
!!IF (.NOT. PRESENT(FINAL)) THEN
      IF (final .EQ. 0) THEN
!!PRJ(ICMIN:K+1) = PKE(I,ICMIN:K+1)
        DO kk=icmin,k+1
          prj(kk) = pke(i, kk)
        END DO
! These initialized here in order not to confuse Valgrind debugger
        poi = 0.
! Do not believe it actually makes any difference.
        qoi = 0.
        uoi = 0.
        voi = 0.
        prs(icmin:k0+1) = ple(i, icmin:k0+1)
        poi(icmin:k) = tho(i, icmin:k)
        qoi(icmin:k) = qho(i, icmin:k)
        uoi(icmin:k) = uho(i, icmin:k)
        voi(icmin:k) = vho(i, icmin:k)
        wsp(icmin:k) = SQRT((uoi(icmin:k)-uoi(k))**2 + (voi(icmin:k)-voi&
&         (k))**2)
        qst(icmin:k) = qss(i, icmin:k)
        dqq(icmin:k) = dqs(i, icmin:k)
        IF (do_tracers) THEN
          DO itr=1,itrcr
            xoi(icmin:k, itr) = xho(i, icmin:k, itr)
          END DO
        END IF
!!! Mass fraction of each layer below cloud base
!!! contributed to aggregate cloudbase layer (CBL) 
        massf(:) = wgt0(i, :)
!!! RESET PRESSURE at bottom edge of CBL 
        prcbl = prs(k)
        DO l=k,k0
          prcbl = prcbl + massf(l)*(prs(l+1)-prs(l))
        END DO
        prs(k+1) = prcbl
        prj(k+1) = (prs(k+1)/1000.)**(mapl_rgas/mapl_cp)
        DO l=k,icmin,-1
          pol(l) = 0.5*(prs(l)+prs(l+1))
          prh(l) = (prs(l+1)*prj(l+1)-prs(l)*prj(l))/(onepkap*(prs(l+1)-&
&           prs(l)))
          pki(l) = 1.0/prh(l)
          dpt(l) = prh(l) - prj(l)
          dpb(l) = prj(l+1) - prh(l)
          pri(l) = .01/(prs(l+1)-prs(l))
        END DO
!!!!! RECALCULATE PROFILE QUAN. IN LOWEST STRAPPED LAYER
        IF (k .LE. k0) THEN
          poi(k) = 0.
          qoi(k) = 0.
          uoi(k) = 0.
          voi(k) = 0.
!! SPECIFY WEIGHTS GIVEN TO EACH LAYER WITHIN SUBCLOUD "SUPERLAYER"
          wght = 0.
          DO l=k,k0
            wght(l) = massf(l)*(ple(i, l+1)-ple(i, l))/(prs(k+1)-prs(k))
          END DO
          DO l=k,k0
            poi(k) = poi(k) + wght(l)*tho(i, l)
            qoi(k) = qoi(k) + wght(l)*qho(i, l)
            uoi(k) = uoi(k) + wght(l)*uho(i, l)
            voi(k) = voi(k) + wght(l)*vho(i, l)
          END DO
          IF (do_tracers) THEN
            xoi(k, :) = 0.
            DO itr=1,itrcr
              DO l=k,k0
                xoi(k, itr) = xoi(k, itr) + wght(l)*xho(i, l, itr)
              END DO
            END DO
          END IF
          tmpt(k) = poi(k)*prh(k)
          CALL DQSATPERT(dqq(k), qst(k), tmpt(k), pol(k), 1)
        END IF
!!DPTH_BL = CPBG*POI(K)*( PRJ(K+1)-PRJ(K) )
! seedras(1,2) are both integers passed from
! from GEOS_Moist w/ values 0 - 1000000
! rndu(:) = 1.0*( seedras(1)+seedras(2) )/2000000.
!dh tap produces code with error
        DO jj=1,nrands
          IF (seedras(i, 1)/1000000. .LT. 1e-6) THEN
            rndu(jj) = 1e-6
          ELSE
            rndu(jj) = seedras(i, 1)/1000000.
          END IF
        END DO
        IF (maxdallowed_d .GE. 0.) THEN
          abs0 = maxdallowed_d
        ELSE
          abs0 = -maxdallowed_d
        END IF
        IF (maxdallowed_s .GE. 0.) THEN
          abs1 = maxdallowed_s
        ELSE
          abs1 = -maxdallowed_s
        END IF
!!call congvec( npoints , seedcbl , rndu )
!            DPTH_BL   = ZCBL(I)
        mxdiam(i) = cnv_fraction(i)*abs0 + (1-cnv_fraction(i))*abs1
        IF (maxdallowed_d .GT. 0) mxdiam(i) = mxdiam(i)*rndu(1)**(-(1./&
&           2.))
! Make MXDIAM stochastic
        DO l=k,icmin,-1
!*
          bet(l) = dqq(l)*pki(l)
!*
          gam(l) = pki(l)/(1.0+lbcp*dqq(l))
          IF (l .LT. k) THEN
            ght(l+1) = gam(l)*dpb(l) + gam(l+1)*dpt(l+1)
            gm1(l+1) = 0.5*lbcp*(dqq(l)/(alhl*(1.0+lbcp*dqq(l)))+dqq(l+1&
&             )/(alhl*(1.0+lbcp*dqq(l+1))))
          END IF
        END DO
        tcu(icmin:k) = -(poi(icmin:k)*prh(icmin:k))
        qcu(icmin:k) = -qoi(icmin:k)
        rns = 0.
        cll = 0.
        rmf = 0.
        rmfd = 0.
        rmfc = 0.
        rmfp = 0.
        cll0 = 0.
        dll0 = 0.
        cllx = 0.
        dllx = 0.
        clli = 0.
        cllb = 0.
        poi_sv = poi
        qoi_sv = qoi
        uoi_sv = uoi
        voi_sv = voi
        IF (do_tracers) xoi_sv = xoi
        lambdsv = 0.0
        cvw = 0.0
        updfrc = 0.0
        updfrp = 0.0
        dissk0 = 0.0
      END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!    IF (PRESENT(FINAL)) THEN
      IF (final .EQ. 1) THEN
        tho(i, icmin:k-1) = poi(icmin:k-1)
        qho(i, icmin:k-1) = qoi(icmin:k-1)
        uho(i, icmin:k-1) = uoi(icmin:k-1)
        vho(i, icmin:k-1) = voi(icmin:k-1)
        cnv_updfrc(i, icmin:k-1) = updfrc(icmin:k-1)
!======================AER_CLOUD=============
!               CNV_NDROP   (I,ICMIN:K-1)  =    CNVNDROP(ICMIN:K-1) !DONIF
!               CNV_NICE   (I,ICMIN:K-1)   =     CNVNICE(ICMIN:K-1) !DONIF
!               CNV_FICE   (I,ICMIN:K-1)   =     CNVFICE(ICMIN:K-1) !DONIF
!! De-strap tendencies from RAS
!! specify weighting "SHAPE"
        wght = wgt1(i, :)
!! Scale properly by layer masses
        wght0 = 0.
        DO l=k,k0
          wght0 = wght0 + wght(l)*(ple(i, l+1)-ple(i, l))
        END DO
        wght0 = (prs(k+1)-prs(k))/wght0
        wght = wght0*wght
        DO l=k,k0
          tho(i, l) = tho(i, l) + wght(l)*(poi(k)-poi_sv(k))
          qho(i, l) = qho(i, l) + wght(l)*(qoi(k)-qoi_sv(k))
          uho(i, l) = uho(i, l) + wght(l)*(uoi(k)-uoi_sv(k))
          vho(i, l) = vho(i, l) + wght(l)*(voi(k)-voi_sv(k))
        END DO
        IF (do_tracers) THEN
          xho(i, icmin:k-1, :) = xoi(icmin:k-1, :)
          DO itr=1,itrcr
            DO l=k,k0
              xho(i, l, itr) = xho(i, l, itr) + wght(l)*(xoi(k, itr)-&
&               xoi_sv(k, itr))
            END DO
          END DO
        END IF
!            FLX (I,ICMIN:K) = RMF (ICMIN:K) * DDT/DAYLEN  !  (KG/m^2/s @ CLOUD BASE)
!  (KG/m^2/s @ CLOUD TOP)
        flxd(i, icmin:k) = rmfd(icmin:k)*ddt/daylen
!            FLXC(I,ICMIN:K) = RMFC(ICMIN:K) * DDT/DAYLEN  !  (KG/m^2/s @ CLOUD TOP)
!  (KG/m^2/s )
        clw(i, icmin:k) = cll(icmin:k)*ddt/daylen
        IF (PRESENT(disske)) disske(i, icmin:k-1) = dissk0(icmin:k-1)*&
&           ddt/daylen
!            FLX (I,1:ICMIN-1) = 0.
        flxd(i, 1:icmin-1) = 0.
!            FLXC(I,1:ICMIN-1) = 0.
        clw(i, 1:icmin-1) = 0.
        IF (k .LT. k0) THEN
!               FLX (I,K:K0) = 0.
          flxd(i, k:k0) = 0.
!               FLXC(I,K:K0) = 0.
          clw(i, k:k0) = 0.
        END IF
!            IRC (I,ICMIN:K-1) = RC(ICMIN:K-1)
        IF (ALLOCATED(icl_v)) THEN
          DEALLOCATE(icl_v)
        END IF
      END IF
      IF (final .EQ. 2) THEN
!            FLX (I,:) = 0.
        flxd(i, :) = 0.
!            FLXC(I,:) = 0.
        clw(i, :) = 0.
!            IRC (I,ICMIN:K-1) = RC(ICMIN:K-1)
      END IF
      RETURN
    END SUBROUTINE STRAP
  END SUBROUTINE RASE
!  Differentiation of rase in reverse (adjoint) mode, forward sweep (with options r8 split(GEOS_MoistGridComp.PRE_RASE GEOS_Moist
!GridComp.PRE_PROGNO_CLOUD ras.RASE rase.CLOUDE rase.ACRITN rase.RNEVP rase.HTEST rase.FINDDTLS rase.STRAP ras.SUNDQ3_ICE cloudne
!w.pdf_spread cloudnew.fix_up_clouds cloudnew.meltfrz cloudnew.hystpdf cloudnew.pdffrac cloudnew.pdfcondensate cloudnew.cnvsrc cl
!oudnew.evap3 cloudnew.subl3 cloudnew.autocon3 cloudnew.PRECIP3 cloudnew.ICEFALL cloudnew.SETTLE_VEL cloudnew.MARSHPALMQ2 cloudne
!w.MICRO_AA_BB_3 cloudnew.LDRADIUS3 cloudnew.ICE_FRACTION cloudnew.GET_ALHX3 cloudnew.ICEFRAC cloudnew.SUNDQ3_ICE3 qsat_util.DQSA
!TPERT)):
!   gradient     of useful results: clw cnv_prc3 tho qho vho cnv_updfrc
!                uho flxd
!   with respect to varying inputs: tho dqs qho vho tpert uho qss
!         pko, plo, phio, phie, qlo, qio,                  &
!FLX, FLXC,                  &
!         CNV_CVW,                                         &
!         CNV_QC,                                          &
!         ENTLAM,                                          &
!         CLAN,                                            &
!         HHO, HSO,PRECU,                                  &
!         AEROPROPS,                                   & !!!!!AER_CLOUD
!         CNV_FICE,                                        & 
!         CNV_NICE,                                        & 
!         CNV_NDROP,                                    &   !DONIF
  SUBROUTINE RASE_FWD(idim, irun, k0, icmin, dt, doras, cpo, alhlo, &
&   alhl1, tice, gravo, seedras, sige, kcbl, wgt0, wgt1, tpert, qpert, &
&   tho, qho, uho, vho, qss, dqs, cnv_fraction, rasal2_2d, co_auto, ple&
&   , pke, clw, flxd, cnv_prc3, cnv_updfrc, rasparams, itrcr, xho, &
&   triedlev_diag, fscav, disske)
    IMPLICIT NONE
!!!!!!!!!!!!!!!======================================     
!!!!!!!!!!   Subroutine ARGact: finds the activated droplet number from Abdul_Razzak and Ghan 2000.
!      ! Tailored for GOCART AEROSOL and works with AERO input 
!      !Written by Donifan Barahona
!      !!donifan.o.barahona@nasa.gov
!!!!!!!!!!!!!!!====================================
!
!      SUBROUTINE ARGact (TEMP, WX, NCPL_ACT, NCPL_AMB,  CDUST, CSOOT, LEV, ISBASE, DDUSTAMB, DSOOTAMB, ENT_PARAM)
!         !
!         Integer, intent(in)     ::  LEV    
!         LOGICAL,  intent(in)     ::   ISBASE
!
!         REAL(8), intent(inout)     ::   TEMP, WX, ENT_PARAM
!         REAL(8), intent(out)     ::   NCPL_ACT, NCPL_AMB, CSOOT, DSOOTAMB
!         REAL(8), DIMENSION(NDUSTMAX), INTENT(OUT) :: CDUST, DDUSTAMB
!         integer                 :: INDEX, NMODES, naux       
!
!         type(AerProps)  :: AER, auxaer      
!
!
!         REAL(8)     ::      kappa, alfa, beta, Akoh, G, T, smax, fi, gi, nui, &
!                       citai, ui, aux1, PACT,  Ntot, auxx, aux, auxconc, W, alph, aseasalt_aux, f_seasalt1
!         REAL(8), dimension (30) ::  SMI, TPI, SIGI 
!
!
!         SMI=0.0      
!         TPI = 0.0
!         SIGI =2.0
!         NCPL_ACT=0.0
!         NCPL_AMB=0.0
!         CDUST=0.0
!         CSOOT=0.0
!         DDUSTAMB =1.0e-9
!         DSOOTAMB= 1.0e-9
!         W=MIN(WX*(1.0-ENT_PARAM), 20.0)    
!         call init_Aer(AER)  
!         AER  =   AERO(LEV) 
!
!         PACT=0.0 !activation probability of entrained aerosol 
!   auxconc =0.0
!    aseasalt_aux  = 0.0
!!!!!!!!!!!activate aerosol transported from cloud base 
!             NMODES =  AER_BASE%nmods
!             TPI(1:nmodes) = AER_BASE%num(1:nmodes)
!             SIGI(1:nmodes) = AER_BASE%sig(1:nmodes)                          
!
!             
!             Ntot= 0.0
!              do index = 1, nmodes 
!	              if (AER_BASE%kap(index) .gt. 0.1) Ntot =  Ntot + TPI(index)  
!              end do
!         
!
!         if ((Ntot .lt. 1.0e4) .or. (TEMP .lt. 245.0) .or. (W .lt. 0.01)) then !no activation if aerosol < 1e-4 1/cm3           
! 
!            NCPL_ACT  = 0.0
!         else
!
!            ! Calculate constants. These fits were obtained from detailed correlations of physical properties. G is actually 1/G
!            T = min(max(TEMP, 243.0), 323.0)     
!            alfa=2.8915E-08*(T**2) - 2.1328E-05*T + 4.2523E-03
!            beta=exp(3.49996E-04*T**2 - 2.27938E-01*T + 4.20901E+01)
!            G=exp(-2.94362E-06*T**3 + 2.77941E-03*T**2 - 8.92889E-01*T + 1.18787E+02)
!            Akoh= 0.66e-6/T  !from Seinfeld and Pandis (1998)
!     
!            !=======================================================
!            !Activate droplets   
!            !=======================================================
!            !Calculate maximum supersaturation according to ARG2002
!
!            auxx=0.0 
!            
!            
!            DO INDEX = 1, NMODES            
!                
!               kappa=  max(AER_BASE%kap(INDEX), 0.001)
!             
!                  SMI (INDEX) = ((0.667*Akoh/AER_BASE%dpg(INDEX))**1.5)/SQRT(2.0*kappa)   ! Critical supersat for mode I   
!                  SMI=MAX(SMI, 1.0e-5)   
!                   
!              if ((TPI(INDEX) .gt. 1e4) .and.  (kappa .gt. 0.1)) then                       
!                  fi=0.5*exp(2.5*SIGI(INDEX)) !sigi is now log(sigi)
!                  gi=1.0+0.25*SIGI(INDEX)
!                  nui=((alfa*W*G)**1.5)/(2.0*MAPL_PI*980.0*beta*TPI(INDEX))
!                  citai = 0.667*Akoh*SQRT(alfa*W*G)
!                  aux1=fi*((citai/nui)**1.5) + gi*(SMI(INDEX)*SMI(INDEX)/(nui+(3.0*citai)))**0.75
!                  aux1=aux1/(SMI(INDEX)*SMI(INDEX))      
!                  auxx=auxx+aux1                  
!                end if
!            end do
!
!  !Calculate number of activated droplets
!            if (auxx .gt. 0.0) then
!               smax = 1/sqrt(auxx)
!               auxx=0.0
!
!                   DO INDEX = 1, NMODES
!                        if ((TPI(INDEX) .gt. 1e4) .and. (AER_BASE%kap(index) .gt. 0.1)) then
!                           ui=sqrt(2.0)*log(SMI(INDEX)/smax)/3.0
!                           aux1=0.5*TPI(INDEX)*(1.0-ERFAPP(ui))
!                           auxx=auxx+aux1
!                           AER_BASE%num(index) = max(TPI(INDEX) -aux1, 0.0) !remove already activated aerosol
!                        end if                    
!                   END DO
!                  NCPL_ACT=auxx             
!            else
!                  NCPL_ACT = 0.0             
!            end if
!
!     
! 
!       end if
!
!         !now filllup dust and soot number
!         NMODES =  AER%nmods
!
!         call getINsubset(1, AER,  auxaer)
!         CDUST(1:auxaer%nmods)= auxaer%num(1:auxaer%nmods)
!         DDUSTAMB(1:auxaer%nmods)= auxaer%dpg(1:auxaer%nmods)
!         call getINsubset(2, AER,  auxaer)
!         naux = max(auxaer%nmods, 1)
!         CSOOT= sum(auxaer%num) 
!         DSOOTAMB= sum(auxaer%dpg)/naux
!
!
!         PACT=1.0 ! fraction of entrained aerosol that is activated
!         auxconc =0.0
!         aseasalt_aux  = 0.0
!
!         do index = 1, nmodes 
!	       if (AER%kap(index) .gt. 0.8)  auxconc = AER%num(index) + auxconc
!           if (AER_BASE%kap(index) .gt. 0.8)    aseasalt_aux  = aseasalt_aux  + AER_BASE%num(index)*AER_BASE%dpg(index)*AER_BASE
!%dpg(index)*1.61*MAPL_PI !assumes a fixed sigma = 2.0
!
!         end do
!       aseasalt = max(aseasalt, aseasalt_aux)
!     
!	  NCPL_AMB=auxconc !Activate  entrained aerosol with kappa>0.8   
!
!
!
!      END SUBROUTINE ARGact
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      !=====Subroutine INfreezing=========
!      ! Freeze droplets in immersion and contact ice nucleation modes, according to Barahona et al. GMD (2014)
!!!!!!!!!!!!!!!!!!!!!!
!
!      subroutine INfreezing(QL, NL, NIN, NDUST, NSOOT, INDUST, INSOOT, TEMP, PRE, WSUB_, DDUST, DSOOT)  !NIN freezing tendency
!         REAL(8), INTENT( IN) :: TEMP, NSOOT, WSUB_, PRE, QL, NL, DSOOT
!         REAL(8), DIMENSION(NDUSTMAX), INTENT (IN) ::  NDUST, DDUST
!         REAL(8), DIMENSION(NDUSTMAX), INTENT (OUT) ::  INDUST
!         REAL(8), INTENT(OUT) ::  NIN, INSOOT 
!         REAL(8), DIMENSION(NDUSTMAX) ::  INDUST_C
!
!         REAL(8) :: a, b, c , d, Tx, n05, ui, aux, nssoot, nsdust, ninbkg, SI, acorr, &
!               dnsd, dnss, coolr, WSUB, ahet, INSOOT_C
!
!         REAL(8) :: nssoot_c, nsdust_c, mfp, nslip, ndfaer, viscosity, lam_r, taux, rho_a, &
!                fdust_drop, fsoot_drop, min_ns_dust, min_ns_soot, nsss, INsea, dnsss, min_ns_seasalt 
!
!         logical :: demott, Drop_mediated
!         integer :: ix
!
!         min_ns_dust= 3.75e6 !limits ice nucleation to -12 !new 02/10/14 
!         min_ns_soot= 3.75e9 !limits ice nucleation to -18
!         min_ns_seasalt = 4.0e2 !limits ice nucleation to -5
!         
!         demott=.false.
!         INDUST=0.0
!         INSOOT=0.0
!         INDUST_C=0.0
!         INSOOT_C=0.0   
!         NIN=0.0
!         Drop_mediated = .false.
!       INsea = 0.0 ! sea salt only in immersion
!   
!   
!   ! note for sea salt we just assume that it is equal to the current droplet concentration and take the area from the calculati
!on at cloud base
!
!         ! fraction of dust and soot within droplets
!         fdust_drop= FDROP_DUST
!         fsoot_drop = FDROP_SOOT
!
!
!         WSUB=MAX(MIN(WSUB_, 10.0), 0.8)
!         coolr=5.0e-3*WSUB  !Approximation to saturated cooling rate 
!         n05=sum(NDUST)+NSOOT
!
!         if (TEMP .gt. T_ICE_MAX) then 
!            return
!         end if
!
!         if (TEMP .lt. T_ICE_ALL   ) then 
!            return
!         end if
!
!
!         if ((QL .le. 1.0e-10) .or. (NL .le. 1.0)) then 
!            return
!         end if
!
!         !Background IN
!         ! SI at water saturation
!
!         rho_a = PRE*100.0/MAPL_RGAS/TEMP 
!         ninbkg=0.0
!         SI = -1.2379e-2+3.3595 !Ice supersat at water sat. Derived from Murphy and Koop 2005
!         !if (TEMP .lt. 260.0)  ninbkg=coolr*42.8*exp(3.88*si)*0.1 !tendency in IN from background IN. Derived from Phillips et 
!al. 2007
!
!
!         Tx = max(TEMP-273.16, -38.0 )
!
!         lam_r=min((MAPL_PI*950.0*NL/rho_a/QL)**(1./3.), 1.0e8)
!         viscosity=1.8e-5*(TEMP/298.0)**0.85    ! Viscosity (kg/m/s)
!         mfp=2.0*viscosity/(PRE  &                   ! Mean free path (m)
!               *sqrt(8.0*28.96e-3/(MAPL_PI*8.314409*TEMP)))        
!
!         if ((n05 .gt.1.0) .and. (TEMP .lt. 272.0)) then
!
!            nsdust=  max(exp(-0.517*Tx + 8.934)-min_ns_dust, 0.0) !From Niemand 2012
!            nssoot= max(1.0e4*exp(-0.0101*Tx*Tx - 0.8525*Tx + 0.7667)-min_ns_soot, 0.0) !Murray (review_ 2012)
!            dnsd  = 0.517*nsdust
!            dnss  = max(-(-2.0*0.0101*Tx -0.8525)*nssoot, 0.0)
!
!            !ns in  contact. It is assumed that in contact is T-3 immersion
!            taux=max(Tx-3.0, -35.0)
!            nsdust_c= max(exp(-0.517*taux + 8.934)-min_ns_dust, 0.0) !From Niemand 2012
!            nssoot_c= max(1.0e4*exp(-0.0101*taux*taux - 0.8525*taux + 0.7667)-min_ns_soot, 0.0) !Murray (review_ 2012)
!
!            aux=0.0        
!            acorr=2.7e7 !m2/m3 correction to the area due to non sphericity and aggregation. Assumes 10 m2/g (Murray 2011)
!
!
!            DO ix=1, NDUSTMAX
!               !Immersion
!               ahet=0.52*DDUST(ix)*DDUST(ix)*DDUST(ix)*acorr*exp(4.5*log(2.0)*log(2.0)*log(2.0))    !this needs to be improved
!
!               INDUST(ix) = NDUST(ix)*exp(-nsdust*ahet)* &
!                     dnsd*coolr*ahet*fdust_drop
!               !Contact   
!               nslip =1.0+(2.0*mfp/DDUST(ix))*(1.257+(0.4*exp(-(1.1*DDUST(ix)*0.5/mfp))))! Slip correction factor               
!     
!               ndfaer =1.381e-23*TEMP*nslip*(1.0-exp(-nsdust_c*ahet)) /(12.*MAPL_PI*viscosity*DDUST(ix))             
!               INDUST_C(ix) = 2.0*MAPL_PI*ndfaer*NDUST(ix)*NL/lam_r
!
!            END DO
!
!
!            acorr=8.0e7 !m2/m3 correction to the area due to non sphericity and aggregation  Assumes 50 m2/g (Popovicheva 2003) 
!      
!            ahet =0.52*DSOOT*DSOOT*DSOOT*acorr*exp(4.5*log(2.0)*log(2.0)*log(2.0))             
!            INSOOT=fsoot_drop*NSOOT*exp(-nssoot*ahet)*dnss*ahet*coolr !
!
!            nslip =1.0+(2.0*mfp/DSOOT)*(1.257+(0.4*exp(-(1.1*DSOOT*0.5/mfp))))! Slip correction factor                    
!            ndfaer =1.381e-23*TEMP*nslip*(1.0-exp(-nssoot_c*ahet)) /(12.*MAPL_PI*viscosity*DSOOT)             
!            INSOOT_c= 2.0*MAPL_PI*ndfaer*NSOOT*NL/lam_r
!
!         ! sea salt
!         nsss =  -0.459*TEMP +128.6235 ! from Demott et al. PNAS, 2015  
!         nsss=  max(exp(nsss)-min_ns_seasalt, 0.0)           
!    	 dnsss=  max(0.459*nsss, 0.0)
!         INsea= aseasalt*dnsss*coolr 
!
! 
!         end if
!
!	    NIN =ninbkg+ INSOOT + SUM(INDUST) + INSOOT_C + SUM(INDUST_C) + INsea!
!         INSOOT=INSOOT +INSOOT_C
!         INDUST =INDUST + INDUST_C
!
!	    
!      end subroutine INfreezing
!
!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      !===================Subroutine Qgrowth======================
!      !Partitions water and ice according to Korolev and Mazin 2005. Assume spheres by now.
!!!!!!!!!!!!!!!!!!!!!
!      subroutine Qgrowth(TEMP, PRE, QICE0, NIPRE, QT, NINUC, DQIG, RIM, FNDRIM) !freezing of IN according to Demott et al 2010 (
!everything SI)
!         REAL(8), INTENT(IN) :: TEMP, QICE0, NIPRE, NINUC,  PRE, QT
!         REAL(8), INTENT(INOUT) :: DQIG, RIM, FNDRIM
!
!         !REAL(8) :: A, denI, aux, Dco, SI, denA, DQold, DQnew 
!         REAL(8) :: DIFF, DENICE, DENAIR, K1, K2, K3, SI, AUX, DC, TEFF, PLo, TEo, TC, &
!               DQnew, DQold, rho_a, LWC, IWC, qmin_rim
!
!
!         if (TEMP .gt. 272.15) then
!            DQIG =0.0
!            return
!         end if
!
!         TC=TEMP-273.0 
!         PLo = max(PRE, 10.0) !limits  of the correlations 
!         TEo = max(190.0, TEMP)
!
!         qmin_rim = 1.0e-12
!
!
!         DENICE= 1000.0*(0.9167 - 1.75e-4*TC -5.0e-7*TC*TC) !From PK 97
!         DENAIR= PLo*100.0/MAPL_RGAS/TEMP
!         DIFF=(0.211*1013.25/(PLo+0.1))*(((TEo+0.1)/273.0)**1.94)*1e-4  !From Seinfeld and Pandis 2006
!
!         K1 = EXP(7.1170e-4*TEo*TEo-0.43563*TEo+78.744) 
!         K2 = EXP(-9.933e-3*TEo+25.26)
!         K3 = EXP(7.1772e-4*TEo*TEo-0.44055*TEo+73.996)
!
!
!         AUX= 210368.0 + 131.438*TEMP - (3.32373E6/TEMP)- (41729.1*LOG(TEMP)) !From Murphy and Koop 2005
!         SI=exp(-aux/8.314/TEMP)-1.0 !ratio of pw/pi-1
!         rho_a = PRE*100.0/MAPL_RGAS/TEMP 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!         if  ((NIPRE .gt. 1.0) .and. (QICE0 .gt. 1.0e-10)) then 
!            DC=max((QICE0/(NIPRE*500.0*MAPL_PI))**(0.333), 40.0e-6) !Assumme monodisperse size distribution about size distribut
!ion is made. 
!         else   
!            DC = 40.0e-6
!         end if
!
!         AUX=  NIPRE*DENICE*MAPL_PI*DC*DC
!         TEFF = DENAIR*2.0*((K1*DIFF+K2)*DC+(K3/0.1))
!
!         if  (AUX .gt. 1.0e-12) then 
!            TEFF=min(TEFF/AUX, 1.0e20)   
!            DQold= SI/TEFF  
!         else
!            DQold=0.0
!         end if
!
!         ! Calculate rimming fraction
!         IWC =  QICE0*rho_a
!         LWC =  max(QT-QICE0, 0.0)*rho_a
!         aux = DQold ! only due to deposition
!
!
!         !Account for rimming
!
!         if ((LWC .gt. qmin_rim)  .and. (IWC .gt. qmin_rim))  then 
!            RIM = 6.0e-5/(LWC*(IWC**0.17)) !Fom Lin and Colle, NRW, 2011
!            RIM = 1.0/(1.0+RIM)
!            RIM  = min (0.95, RIM)
!            DQold =  DQold*(1 + RIM/(1.0-RIM))
!            FNDrim =  max(min(rho_a*(DQold -aux)/LWC, 1.0), 0.0) !Fraction of liquid condensate removed due to riming   
!         END if
!
!
!
!
!!!!!!!!!!!!!!!!!!!!!!recently nucleated!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!         AUX=  NINUC*DENICE*MAPL_PI*20.0e-6*20.0e-6
!         TEFF = DENAIR*2.0*((K1*DIFF+K2)*DC+(K3/0.1))
!
!         if  (AUX .gt. 1.0e-12)  then 
!            TEFF=min(TEFF/AUX, 1.0e10)
!            DQnew= SI/TEFF
!         else
!            DQnew = 0.0
!         end if
!
!         DQIG = DQold+DQnew     
!
!      end subroutine Qgrowth
!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!
!
!      !*************************************************************
!      ! Function PDG07 (simplified background ice nucleation 
!      !                     spectra according to Phillips et. al. 2007).  
!      ! si is supersaturation wrt ice and T is in K 
!      !************************************************************  
!
!      subroutine PDG07_ice(si, Tx, N)     
!
!         REAL(8), intent(in) :: si, Tx
!         REAL(8), intent(out)  :: N 
!         N=0.0
!
!         !if (Tx .le. 243.0)then
!
!
!         N=1000.0*exp(-0.388)*(exp(3.88*si)-1.0)/0.76
!         !elseif (Tx .le. 260.0) then
!         !  N=60.0*exp(-0.639)*(exp(12.96*si)-1.0)/0.76   
!         !end if      
!      end subroutine PDG07_ice
!
!
!      !*********************************************************************
!*********************************************************************
!*********************************************************************
!******************** Relaxed Arakawa-Schubert ***********************
!************************ Parameterization ***************************
!********************** SCALAR RAS-1 VERSION  ************************
!************************* 31 DECEMBER 1999 **************************
!*********************************************************************
!************************** Developed By *****************************
!*********************************************************************
!************************ Shrinivas Moorthi **************************
!******************************* and *********************************
!************************** Max J. Suarez ****************************
!*********************************************************************
!******************** Laboratory for Atmospheres *********************
!****************** NASA/GSFC, Greenbelt, MD 20771 *******************
!*********************************************************************
!*********************************************************************
!  Input:
!  ------
! 
!     K0      : Number of vertical levels (increasing downwards)
!
!     DT      : Time step in seconds
!
!     RASAL   : Array of dimension K-1 containing relaxation parameters
!               for cloud-types detraining at those levels
!
!     CPO     : Specific heat at constant pressure (J/kg/K)
!
!     ALHLO   : Latent Heat of condensation (J/kg)
!
!     ALHL1   : Latent Heat of condensation + fusion (J/kg)
!
!     GRAVO   : Acceleration due to gravity (m/s^2)
!
!     PLE     : 2D array of dimension (IDIM,K0+1) containing pressure
!               in hPa at the interfaces of K-layers from top of the 
!               atmosphere to the bottom  (mb)
!
!     PKE     : 2D array of dimension (IDIM,K0+1) containing (PRS/P00) **
!               RKAP.  i.e. Exner function at layer edges.
!
!     PKL     : 2D array of dimension (IDIM,K0) ) containing the
!               Exner function at the layers.
!
!     QSS     : 2D array of dimension (IDIM,K0  ) containing the
!               saturation specific humidity at the layers. (kg/kg)
!
!     DQS     : 2D array of dimension (IDIM,K0  ) containing
!               d(qss)/dt at the layers.  (1/K)
!   
!     CNV_FRACTION    : 1D array of dimension (IDIM) containing
!               fraction of grid cell considered to be convective
!   
!  Update:
!  -------
!
!     THO     : 2D array of dimension (IDIM,K0) containing potential
!               temperature (K)
!
!     QHO     : 2D array of dimension (IDIM,K0) containing specific
!               humidity (kg/kg)
!
!     UHO     : 2D array of dimension (IDIM,K0) containing u-wind (m/s)
!
!     VHO     : 2D array of dimension (IDIM,K0) containing v-wind (m/s)
!
!  Output:
!  -------
!!
!     CLW     : 2D array of dimension (IDIM,K0) containing the
!               detrained cloud liquid water.  (kg/m^2/s) 
!
!     FLX     : 2D array of dimension (IDIM,K0) containing the
!               cloud-base mass flux for each cloud type ordered by
!               detrainment level.   (kg/m^2/s) 
!
!     FLXD    : 2D array of dimension (IDIM,K0) containing the
!               detrained  mass flux for each cloud type ordered by
!               detrainment level.   (kg/m^2/s) 
!
!     FLXC    : 2D array of dimension (IDIM,K0+1) containing the
!               total cloud mass flux for all cloud types through
!               the top of each level. (e.g., FLXC(K)=SUM(FLX(ICMIN:K))
!               and  FLXD(L) = FLXC(L+1)-FLSD(L) )
!                (kg/m^2/s) 
!
!     PRECU   : 1D (IDIM) Locally-handled convective precip
!               Zero if older version of RAS-1 is used. Nonzero if
!               RAS-2 is used.
!
!   AEROPROPS, Structure containing aerosol propoerties (in)
!   CNV_NICE, CNV_DROP.  Flux of ice crystals and droplet number at det level (1/m^2/s) (out)
!   CNV_FICE: Ice fraction in the detrained condensate (out)
!************************************************************************
!  ARGUMENTS
    INTEGER, INTENT(IN) :: idim, irun, k0, icmin
!, QLO, QIO, CLAN
    REAL*8, DIMENSION(idim, k0), INTENT(INOUT) :: tho, qho, uho, vho
!,PHIE
    REAL*8, DIMENSION(idim, k0 + 1), INTENT(IN) :: ple, pke
!,PLO,PKO,PHIO
    REAL*8, DIMENSION(idim, k0), INTENT(IN) :: qss, dqs
    REAL*8, DIMENSION(idim), INTENT(IN) :: cnv_fraction, rasal2_2d
    REAL*8, DIMENSION(k0 + 1), INTENT(IN) :: sige
!, FLX
    REAL*8, DIMENSION(idim, k0) :: clw, flxd
!      REAL(8), DIMENSION (IDIM,K0+1), INTENT(  OUT) ::  FLXC
    REAL*8, DIMENSION(idim, k0) :: cnv_prc3
!, CNV_CVW, CNV_QC
    REAL*8, DIMENSION(idim, k0) :: cnv_updfrc
!      REAL(8), DIMENSION (IDIM,K0  ), INTENT(  OUT) ::  ENTLAM
!      REAL(8), DIMENSION (IDIM,K0  ), INTENT(  OUT) ::  HHO, HSO
    REAL*8, INTENT(IN) :: dt, cpo, alhlo, gravo
    REAL*8, INTENT(IN) :: alhl1, tice
    INTEGER, INTENT(IN) :: itrcr
    INTEGER, DIMENSION(idim, 2), INTENT(IN) :: seedras
    INTEGER, DIMENSION(idim), INTENT(IN) :: kcbl
    INTEGER, DIMENSION(idim), INTENT(IN) :: doras
    REAL*8, DIMENSION(idim), INTENT(IN) :: tpert, qpert
    REAL*8, DIMENSION(idim), INTENT(IN) :: co_auto
    REAL*8, DIMENSION(idim) :: mxdiam
!      REAL(8), DIMENSION (IDIM     ), INTENT(  OUT) ::  PRECU
    REAL*8, DIMENSION(idim, k0), INTENT(IN) :: wgt0, wgt1
!     REAL(8), DIMENSION(:),          INTENT(IN   ) ::  RASPARAMS
    TYPE(RASPARAM_TYPE), INTENT(IN) :: rasparams
!      INTEGER, DIMENSION(IDIM,K0), INTENT(  OUT) ::  IRC
    REAL*8, OPTIONAL, INTENT(INOUT) :: xho(idim, k0, itrcr)
    REAL*8, OPTIONAL, INTENT(OUT) :: triedlev_diag(idim, k0)
    REAL*8, OPTIONAL, INTENT(OUT) :: disske(idim, k0)
! Fraction scavenged per km
    REAL*8, OPTIONAL, INTENT(IN) :: fscav(itrcr)
! = 0 (no scav), = 1 (full scav)
!  LOCALS
    REAL*8, DIMENSION(k0) :: poi_sv, qoi_sv, uoi_sv, voi_sv
    REAL*8, DIMENSION(k0) :: poi, qoi, uoi, voi, dqq, bet, gam, cll, &
&   tmpt
    REAL*8, DIMENSION(k0) :: poi_c, qoi_c
    REAL*8, DIMENSION(k0) :: prh, pri, ght, dpt, dpb, pki, dissk0, &
&   dissk1, clantnd
    REAL*8, DIMENSION(k0) :: tcu, qcu, ucu, vcu, cln, rns, pol, dm
    REAL*8, DIMENSION(k0) :: qst, ssl, rmf, rnn, rn1, rmfc, rmfp
    REAL*8, DIMENSION(k0) :: gms, eta, gmh, eht, gm1, hcc, rmfd
    REAL*8, DIMENSION(k0) :: hol, hst, qol, zol, hcld, cll0, cllx, clli&
&   , cllb
    REAL*8, DIMENSION(k0) :: wsp, lambdsv, bke, cvw, updfrc
    REAL*8, DIMENSION(k0) :: rasal, mtkwi, updfrp, bk2, bk3, dll0, dllx
    REAL*8, DIMENSION(itrcr) :: xht
    REAL*8, DIMENSION(k0, itrcr) :: xoi, xcu, xoi_sv
    REAL*8, DIMENSION(k0 + 1) :: prj, prs, qht, sht, zet, xyd, xyd0
!      INTEGER,  DIMENSION(K0-1) :: RC
    INTEGER :: k, my_pe
    REAL*8, DIMENSION(idim, k0) :: lambdsv2
    REAL*8 :: tx2, tx3, uht, vht, akm, acr, alm, tth, qqh, shtrg, wspbl&
&   , dqx
!, BKE
    REAL*8 :: wfn, tem, trg, trgexp, evp, wlq, qcc, mtkw_max
    REAL*8 :: shtrg_fac, sige_minhol, wfnog
    INTEGER :: i, ic, l, icl, itr, icl_c, n_dtl
    INTEGER :: ndtlexpon
    INTEGER, DIMENSION(:), ALLOCATABLE :: icl_v
!  RASE GLOBAL CONSTANTS
    REAL*8 :: grav, cp, alhl, cpbg, alhi, cpi, gravi, ddt, lbcp, obg, &
&   afc
!!!!!!!!!
    REAL*8 :: fricfac, dpth_bl, wupdrft, pblfrac, autorampb, co_zdep
    REAL*8 :: rasal1, rasal2, rasal2i, co_t, rasncl, friclambda, sdqvt1&
&   , sdqv2
    REAL*8 :: lambda_fac, strapping, acritfac, hmintrigger, lldisaggxp
    REAL*8 :: lambmx_fac, diammn_min, rdtlexpon, cli_crit, sdqv3, &
&   maxdallowed_d, maxdallowed_s
    REAL*8 :: rhmn, rhmx, cldmicro, fdrop_dust, fdrop_soot
    INTEGER :: kstrap, rasal_exp
    REAL*8 :: cld_radius, areal_frac, spect_mflx, cvw_cbase
!!!!!!!!!
    REAL*8, PARAMETER :: onepkap=1.+2./7., daylen=86400.0
!      REAL(8), PARAMETER :: PBLFRAC = 0.5
    REAL*8, PARAMETER :: rhmax=0.9999
!  LAMBDA LIMITS
    REAL*8 :: lambda_min
    REAL*8 :: lambda_max
!  TRIGGER PARAMETERS
! Density of liquid water in kg/m^3
    REAL*8, PARAMETER :: rho_w=1.0e3
    LOGICAL :: dyna_strapping, do_tracers, smooth_hst
    CHARACTER(len=esmf_maxstr) :: cbl_style
    REAL*8, DIMENSION(k0) :: tcu8, qcu8, pcu, flx8
    REAL*8, DIMENSION(k0, itrcr + 2) :: rcu
!, dpd, tla
    REAL*8 :: cup
    LOGICAL :: revap, wrkfun, calkpb, crtfun, lprnt, dndrft
    REAL*8, DIMENSION(k0) :: toi8, qoi8, prsm8, phil8, qli8, qii8, &
&   trcfac
    REAL*8, DIMENSION(k0) :: alfind, alfint, alfinq, rhc_ls
    REAL*8, DIMENSION(k0 + 1) :: prs8, phih8
    REAL*8, DIMENSION(k0, itrcr + 2) :: roi8
    REAL*8 :: fracbl, dt8, rasalf
    INTEGER :: kpbl
! no inhibition for =1.0
    REAL*8, SAVE :: max_neg_bouy=1.0
!!real*8 :: ALFINT = 0.5
!!real*8 :: ALFINQ = 0.5
! not used
    REAL*8, SAVE :: rhfacl=0.0
! no inhibition
    REAL*8, SAVE :: rhfacs=0.0
! 1 degree resolution
    REAL*8, SAVE :: garea=1.e10
!!real*8 :: ALFIND = 1.0
!!real*8 :: RHC_LS = 0.80
    REAL*8, SAVE :: dsfc=0.001
    REAL*8, SAVE :: cd=1.e-3
    REAL*8, SAVE :: wfnc=0.0
    REAL*8, SAVE :: tla=-1.0
    REAL*8, SAVE :: dpd=300.
!  SCAVANGING RELATED PARAMETERS
! layer thickness in km
    REAL*8 :: delzkm
! fraction of tracer *not* scavenged
    REAL*8 :: fnoscav
! Fraction scavenged per km
    REAL*8 :: fscav_(itrcr)
    INTRINSIC PRESENT
    INTRINSIC INT
    INTRINSIC SQRT
    INTRINSIC SUM
    INTRINSIC ALLOCATED
! ************************AER_CLOUD *********************************************
!      TYPE(AerProps),    DIMENSION(IDIM, K0),    INTENT(IN)    :: AEROPROPS !DONIF
!      REAL(8), DIMENSION (IDIM,K0), INTENT( OUT) ::  CNV_NDROP, CNV_FICE, CNV_NICE
!      REAL(8),  DIMENSION(K0) :: CNVNDROP, CNVNICE, CNVFICE !DONIF
!
!      TYPE(AerProps), DIMENSION(K0) ::  AERO   !AEROSOL VERTICAl ARRAY FOR ALL SPECIES    
!      REAL(8) ::  T_ICE_ALL, T_ICE_MAX, ASEASALT, f_seasalt 
!      INTEGER, PARAMETER :: NDUSTMAX = 10
!
!      INTEGER :: INDEX        
!      TYPE(AerProps) :: AERAUX, AER_BASE
!
!      CNV_FICE  =0.0
!      CNV_NDROP =0.0
!      CNV_NICE  =0.0
!      CNVFICE  =0.0
!      CNVNDROP =0.0
!      CNVNICE  =0.0
!      T_ICE_ALL= 238.0
!      T_ICE_MAX= MAPL_TICE
!      CLDMICRO = 0.0
!      FDROP_DUST = 0.1
!      FDROP_SOOT = 0.01
! *********************************************************************
    IF (PRESENT(fscav)) THEN
      fscav_ = fscav
    ELSE
! NO SCAVENGING BY DEFAULT
      fscav_ = 0.0
    END IF
    IF (irun .LE. 0) THEN
      CALL PUSHCONTROL1B(0)
    ELSE
      cnv_prc3 = 0.
      cnv_updfrc = 0.
!      CNV_CVW   =0. 
!      CNV_QC    =0.
!      ENTLAM    =0.
!!LAMBDSV2 = 0.
!      IRC = -2
!      SMOOTH_HST   = .TRUE. 
      smooth_hst = .false.
      fricfac = rasparams%cufricfac
! MAT CO_AUTO is now passed in from outside
!     in order to allow this code to run over im*jm 
!     columns
!CO_AUTO      = RASPARAMS(3)     !  ---  3
      cli_crit = rasparams%qc_crit_cn
      rasal1 = rasparams%rasal1
      rasncl = rasparams%rasncl
      friclambda = rasparams%cufriclambda
      sdqv2 = rasparams%sdqv2
      sdqv3 = rasparams%sdqv3
      sdqvt1 = rasparams%sdqvt1
      acritfac = rasparams%acritfac
      pblfrac = rasparams%pblfrac
      autorampb = rasparams%rasautorampb
      co_zdep = rasparams%autoc_cn_zdep
      maxdallowed_s = rasparams%maxdallowed_s
      maxdallowed_d = rasparams%maxdallowed_d
      rasal_exp = rasparams%rasal_exp
      rhmn = rasparams%ras_rhmin
      rhmx = rasparams%ras_rhfull
      cldmicro = rasparams%cldmicro
      do_tracers = PRESENT(xho) .AND. itrcr .GT. 0
      grav = gravo
      alhl = alhlo
      cp = cpo
      cpi = 1.0/cp
      alhi = 1.0/alhl
      gravi = 1.0/grav
      cpbg = cp*gravi
      ddt = daylen/dt
      lbcp = alhl*cpi
!      HHO = MAPL_UNDEF
!      HSO = MAPL_UNDEF
!#ifdef RAS2
!      dt8    =  dt
!
!      call set_ras_afc(dt8)
!      call ras_init(k0,1)
!
!      kpbl  = k0
!      lprnt = .false.
!      trcfac(:) = 1.0
!      revap  = .true.
!      wrkfun = .false.
!      calkpb = .true.
!      crtfun = .true.
!      dndrft = .true.
!      wfnc   = 0.0
!
!#endif
      DO i=1,irun
        IF (doras(i) .EQ. 0) THEN
          CALL PUSHCONTROL2B(0)
        ELSE
!===================AER_CLOUD
!            AERO =  AEROPROPS(I, :)
!            CNVFICE  =0.0
!            CNVNDROP =0.0
!            CNVNICE  =0.0
!#ifndef RAS2
!!CALL FINDBASE
          k = kcbl(i)
!         rc(icmin) = 0
          CALL FINDDTLS()
          IF (k .GT. 0) THEN
            CALL STRAP_FWD(0)
            CALL HTEST_FWD()
!            HHO(I,:) = HOL
!            HSO(I,:) = HST
            DO icl_c=1,n_dtl
              CALL PUSHINTEGER4(icl)
              icl = icl_v(icl_c)
              IF (do_tracers) xcu(icmin:, :) = 0.
! This change makes cumulus friction
              ucu(icmin:) = 0.
! correct.
              vcu(icmin:) = 0.
              IF (icl .GT. icmin) THEN
                CALL CLOUDE_FWD(icl)
                CALL PUSHCONTROL1B(1)
              ELSE
                CALL PUSHCONTROL1B(0)
              END IF
            END DO
            CALL PUSHINTEGER4(icl_c - 1)
!               ENTLAM(I,ICL) = ALM
            IF (SUM(rmf(icmin:k)) .GT. 0.0) THEN
              CALL RNEVP_FWD()
              CALL STRAP_FWD(1)
              CALL PUSHCONTROL2B(3)
            ELSE
              CALL STRAP_FWD(2)
              CALL PUSHCONTROL2B(2)
            END IF
          ELSE
            CALL STRAP_FWD(2)
            CALL PUSHCONTROL2B(1)
          END IF
        END IF
      END DO
!         PRECU = 0.0  ! Zero out precip - TBD w/ in progno_cloud
!#else
!
!
!         K = K0
!         CALL FINDDTLS
!
!         !===>  K      INPUT   THE RISE & THE INDEX OF THE SUBCLOUD LAYER
!         !===>  KD     INPUT   DETRAINMENT LEVEL ( 1<= KD < K )
!         !===>  M      INPUT   NUMBER OF TRACERS. MAY BE ZERO.
!         !===>  RASALF INPUT
!         !===>  FRACBL INPUT   MASS FLUX LIMIT
!         !===>  MAX_NEG_BOUY  INPUT  FRACTION OF NEG ALLOWED
!         !===>  ALFINT(K) INPUT   INTERFACE UPWIND (=1) PARAMETER FOR T
!         !===>  ALFINQ(K) INPUT   INTERFACE UPWIND (=1) PARAMETER FOR Q
!         !===>  RHFACL INPUT   CRITICAL RH AT PBL TOP (LAND)
!         !===>  RHFACS INPUT   CRITICAL RH AT PBL TOP (SEA) ???
!         !===>  GAREA  INPUT   AREA OF GRID BOX
!         !===>  ALFIND(K) INPUT   INTERFACE UPWIND (=1) PARAMETER FOR DOWNDRAFT
!         !===>  RHC_LS(K) INPUT   CRITICAL RH FOR RE-EVAP
!
!         ALFINT = 0.5
!         ALFINQ = 0.5
!         ALFIND = 0.5
!         RHC_LS = 0.8
!
!         tcu8  = 0.0
!         qcu8  = 0.0
!         rcu  = 0.0
!         pcu  = 0.0
!         flx8 = 0.0
!         cup  = 0.0
!
!         cloudloop: DO ICL = ICMIN,K0-1
!            IF (DO_TRACERS) XCU(ICMIN:,:) = 0.
!
!
!            rasalf = rasal(ICL)
!            FRACBL = pblfrac
!
!            toi8 = pko(i,:)*tho(i,:) ! tbdone
!            qoi8 = qho(i,:)
!            roi8(:,1) = uho(i,:)
!            roi8(:,2) = vho(i,:)
!            if(present(xho)) roi8(:,3:size(xho,3)+2) = xho(i,:,:) 
!
!            prs8  = ple(i,:)
!            prsm8 = plo(i,:) ! tbdone
!            phil8 = phio(i,:) ! tbdone
!            phih8 = phie(i,:) ! tbdone
!            qli8  = qlo(i,:) ! tbdone
!            qii8  = qio(i,:) ! tbdone
!
!            !===>  TOI(K)     INOUT   TEMPERATURE             KELVIN
!            !===>  QOI(K)     INOUT   SPECIFIC HUMIDITY       NON-DIMENSIONAL
!            !===>  ROI(K,M)   INOUT   TRACER                  ARBITRARY
!            !===>  PRS(K+1)   INPUT   PRESSURE @ EDGES        MB
!            !===>  PRSM(K)    INPUT   PRESSURE @ LAYERS       MB
!            !===>  PHIL(K)    INPUT   GEOPOTENTIAL @ LAYERS IN MKS units
!            !===>  PHIH(K+1)  INPUT   GEOPOTENTIAL @ EDGES  IN MKS units
!            !===>  QLI(K)     INOUT   LIQUID WATER            NON-DIMENSIONAL
!            !===>  QII(K)     INOUT   ICE                     NON-DIMENSIONAL
!            !===>  KPBL       INOUT   PBL TOP
!            !===>  DSFC       INOUT   GUSTINESS (m/s)
!            !===>  CD         INOUT   DRAG COEF. (NON-DIMENSIONAL)
!            !===>  LPRNT      INOUT
!            !===>  TRCFAC     INOUT   tracer pg factor for u AND V 
!            !===>  TCU(K  )   UPDATE  TEMPERATURE TENDENCY       DEG
!            !===>  QCU(K  )   UPDATE  WATER VAPOR TENDENCY       (G/G)
!            !===>  RCU(K,M)   UPDATE  TRACER TENDENCIES          ND
!            !===>  PCU(K-1)   UPDATE  PRECIP @ BASE OF LAYER     KG/M^2
!            !===>  FLX(K  )   UPDATE  UPDRAFT MASS FLUX @ TOP OF LAYER   KG/M^2
!            !===>  CUP        UPDATE  PRECIPITATION AT THE SURFACE KG/M^2
!            !===>  REVAP      UPDATE  LOGICAL TO CONTROL REEVAP IN ENVIRONMENT
!            !===>  DT         UPDATE  TIME STEP USED FOR GUSTINESS
!            !===>  WFNC       UPDATE  WORK FUNCTION (in if CRTFUN = F, out if WRKFUN = T ) 
!            !===>  WRKFUN     UPDATE  LOGICAL TO DO WRKFUNC ONLY
!            !===>  CALKBP     UPDATE  LOGICAL TO COMPUTE KPBL INTERNALLY
!            !===>  CRTFUN     UPDATE  LOGICAL TO CONTROL CRITICAL WORK FUNCTION
!            !===>  TLA        UPDATE  TILTING ANGLE FOR DD (NEG USE TABLE)
!            !===>  DNDRFT     INPUT   LOGICAL TO CONTROL DOWNDRAFT
!            !===>  DPD        INPUT   Minumum Cloud Depth for DOWNDRFAT Computation hPa
!
!            call CLOUD(                                                       &
!                  &                  K0, ICL, size(xcu,2)                            &
!                  &,                 RASALF, FRACBL, MAX_NEG_BOUY                    &
!                  &,                 ALFINT, ALFINQ, RHFACL, RHFACS, garea           &
!                  &,                 ALFIND, RHC_LS                                  &
!                  
!                  &,                 TOI8, QOI8, ROI8, PRS8, PRSM8, phil8, phih8     &
!                  &,                 QLI8, QII8, KPBL, DSFC                            &
!                  &,                 CD,lprnt, trcfac                                &
!                  
!                  &,                 TCU8, QCU8, RCU, PCU, FLX8                         &
!                  &,                 CUP, REVAP, DT8                                 &
!                  &,                 WFNC, WRKFUN, CALKPB, CRTFUN, TLA, DnDRFT, DPD)  
!
!            tho(i,:)  = toi8/pko(i,:)
!            qho(i,:)  = qoi8
!            uho(i,:)  = roi8(:,1) 
!            vho(i,:)  = roi8(:,2)
!            qlo(i,:)  = qli8 ! tbdone
!            qio(i,:)  = qii8  ! tbdone
!
!
!            if(present(xho)) xho(i,:,:) =roi8(:,3:size(xho,3)+2)
!         ENDDO cloudloop
!
!         ! temporary zeroing of outputs needed by rest of moist
!         !-----------------------------------------------------
!         cnv_prc3(i,:)   = 0.0
!         clw(i,:)        = 0.0
!         HHO(i,:)        = 0.0
!         HSO(i,:)        = 0.0
!         cnv_updfrc(i,:) = 0.0001
!         cnv_qc(i,:)     = 0.0
!         flxc(i,1:k0)    = flx8 / DT
!         flx(i,:)        = 0.0
!
!         flxd(i,1:k0-1)  = MAX( flxc(i,2:k0)- flxc(i,1:k0-1) , 0. )
!         flxd(i,k0)      = 0.0
!
!         precu(i)        = cup / DT   ! 
!
!         !
!         ! Tiedtke-style prognostic anvil fraction. 
!         ! Normally done in progno_cloud. Moved here
!         ! since RAS-2 manages its own detraining 
!         ! ice and liquid condensate.
!         !---------------------------------------------
!         dm           = 100.*(prs8(2:k0+1)-prs8(1:k0) )/grav
!         CLANTND      = FLXD(i,:) / DM
!         CLAN(i,:)    = CLAN(i,:) + CLANTND*DT
!         CLAN(i,:)    = MIN( CLAN(i,:) , 0.99 )
!
!         ! Now detraining mass flux has to be zeroed
!         ! or clouds will be crated twice
!         !--------------------------------
!         flxd(i,:)    = 0.
!
!#endif
      IF (ALLOCATED(icl_v)) THEN
        CALL PUSHREAL8(acr)
        CALL PUSHREAL8ARRAY(eht, k0)
        CALL PUSHREAL8ARRAY(gm1, k0)
        CALL PUSHREAL8ARRAY(hcc, k0)
        CALL PUSHREAL8(trg)
        CALL PUSHREAL8(tem)
        CALL PUSHREAL8ARRAY(cvw, k0)
        CALL PUSHREAL8ARRAY(dqq, k0)
        CALL PUSHREAL8(wlq)
        CALL PUSHREAL8ARRAY(ght, k0)
        CALL PUSHREAL8ARRAY(hst, k0)
        CALL PUSHREAL8ARRAY(bet, k0)
        CALL PUSHREAL8ARRAY(qht, k0 + 1)
        CALL PUSHREAL8(vht)
        CALL PUSHREAL8ARRAY(qoi, k0)
        CALL PUSHREAL8ARRAY(qol, k0)
        CALL PUSHREAL8(akm)
        CALL PUSHREAL8ARRAY(voi, k0)
        CALL PUSHREAL8(tx2)
        CALL PUSHREAL8(tx3)
        CALL PUSHREAL8ARRAY(rasal, k0)
        CALL PUSHREAL8ARRAY(hcld, k0)
        CALL PUSHREAL8ARRAY(eta, k0)
        CALL PUSHREAL8ARRAY(pki, k0)
        CALL PUSHREAL8ARRAY(sht, k0 + 1)
        CALL PUSHINTEGER4(n_dtl)
        CALL PUSHREAL8(rasal2i)
        CALL PUSHREAL8ARRAY(gmh, k0)
        CALL PUSHREAL8ARRAY(dpb, k0)
        CALL PUSHREAL8ARRAY(prh, k0)
        CALL PUSHREAL8ARRAY(pri, k0)
        CALL PUSHREAL8ARRAY(prj, k0 + 1)
        CALL PUSHREAL8(lambda_min)
        CALL PUSHREAL8(alm)
        CALL PUSHREAL8ARRAY(rnn, k0)
        CALL PUSHREAL8ARRAY(qst, k0)
        CALL PUSHREAL8(uht)
        CALL PUSHREAL8ARRAY(gms, k0)
        CALL PUSHREAL8ARRAY(poi, k0)
        CALL PUSHREAL8ARRAY(prs, k0 + 1)
        CALL PUSHINTEGER4(icl)
        CALL PUSHREAL8ARRAY(uoi, k0)
        CALL PUSHREAL8ARRAY(ssl, k0)
        CALL PUSHREAL8ARRAY(dpt, k0)
        CALL PUSHREAL8ARRAY(zet, k0 + 1)
        CALL PUSHREAL8ARRAY(zol, k0)
        CALL PUSHREAL8ARRAY(hol, k0)
        CALL PUSHREAL8ARRAY(rmfp, k0)
        CALL PUSHREAL8ARRAY(gam, k0)
        CALL PUSHCONTROL1B(1)
      ELSE
        CALL PUSHREAL8(acr)
        CALL PUSHREAL8ARRAY(eht, k0)
        CALL PUSHREAL8ARRAY(gm1, k0)
        CALL PUSHREAL8ARRAY(hcc, k0)
        CALL PUSHREAL8(trg)
        CALL PUSHREAL8(tem)
        CALL PUSHREAL8ARRAY(cvw, k0)
        CALL PUSHREAL8ARRAY(dqq, k0)
        CALL PUSHREAL8(wlq)
        CALL PUSHREAL8ARRAY(ght, k0)
        CALL PUSHREAL8ARRAY(hst, k0)
        CALL PUSHREAL8ARRAY(bet, k0)
        CALL PUSHREAL8ARRAY(qht, k0 + 1)
        CALL PUSHREAL8(vht)
        CALL PUSHREAL8ARRAY(qoi, k0)
        CALL PUSHREAL8ARRAY(qol, k0)
        CALL PUSHREAL8(akm)
        CALL PUSHREAL8ARRAY(voi, k0)
        CALL PUSHREAL8(tx2)
        CALL PUSHREAL8(tx3)
        CALL PUSHREAL8ARRAY(rasal, k0)
        CALL PUSHREAL8ARRAY(hcld, k0)
        CALL PUSHREAL8ARRAY(eta, k0)
        CALL PUSHREAL8ARRAY(pki, k0)
        CALL PUSHREAL8ARRAY(sht, k0 + 1)
        CALL PUSHINTEGER4(n_dtl)
        CALL PUSHREAL8(rasal2i)
        CALL PUSHREAL8ARRAY(gmh, k0)
        CALL PUSHREAL8ARRAY(dpb, k0)
        CALL PUSHREAL8ARRAY(prh, k0)
        CALL PUSHREAL8ARRAY(pri, k0)
        CALL PUSHREAL8ARRAY(prj, k0 + 1)
        CALL PUSHREAL8(lambda_min)
        CALL PUSHREAL8(alm)
        CALL PUSHREAL8ARRAY(rnn, k0)
        CALL PUSHREAL8ARRAY(qst, k0)
        CALL PUSHREAL8(uht)
        CALL PUSHREAL8ARRAY(gms, k0)
        CALL PUSHREAL8ARRAY(poi, k0)
        CALL PUSHREAL8ARRAY(prs, k0 + 1)
        CALL PUSHINTEGER4(icl)
        CALL PUSHREAL8ARRAY(uoi, k0)
        CALL PUSHREAL8ARRAY(ssl, k0)
        CALL PUSHREAL8ARRAY(dpt, k0)
        CALL PUSHREAL8ARRAY(zet, k0 + 1)
        CALL PUSHREAL8ARRAY(zol, k0)
        CALL PUSHREAL8ARRAY(hol, k0)
        CALL PUSHREAL8ARRAY(rmfp, k0)
        CALL PUSHREAL8ARRAY(gam, k0)
        CALL PUSHCONTROL1B(1)
      END IF
    END IF

  CONTAINS
!*********************************************************************
    SUBROUTINE CLOUDE(ic)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: ic
      REAL*8 :: deep_fact, cu_diam, wscale
!, dQx
      REAL*8 :: cli, te_a, c00_x, cli_crit_x, pete, toki, gmhx, hstx
      REAL*8 :: dt_lyr, rate, cvw_x, closs, f2, f3, f4, f5
      INTEGER :: k700
      INTRINSIC MIN
      INTRINSIC MAX
      INTRINSIC SQRT
      INTRINSIC EXP
      REAL*8 :: min2
      REAL*8 :: min1
      REAL*8 :: x5
      REAL*8 :: x4
      REAL*8 :: x3
      REAL*8 :: x2
      REAL*8 :: x1
      LOGICAL :: mask(k-ic+1)
      REAL*8 :: max2
      REAL*8 :: max1
      REAL*8 :: y1
!         !=============================AER_CLOUD local variables ====================
!         REAL(8) :: WBASE, NDROP, NICE, FP_D, FF_A, FP_I, FICE, &
!               NDROP_AMB, NSOOT_AMB, NSOOT, NIN, INSOOT, dCVW2, QICE, &
!               dQICE, dQIG, FPICE, dNICE, dNDROP, DSOOT_AMB, DSOOT, QLIQ, dQLIQ, FPRECIP, AUX, QT, &
!               MAXNICE, MAXNDROP, MINNICE, MINNDROP, NDROP_ACT, RIMM, FNDRIM, TminusTa, Tparcel , &
!	      alph_e, beta_e, RH_AMB, ECRIT
!
!         REAL(8), DIMENSION (NDUSTMAX) :: NDUST, NDUST_AMB, INDUST, DDUST_AMB, DDUST 
!         INTEGER :: INX
!
!
!         T_ICE_ALL = 238.0 !A little higher so there is ice at the freezing level
!         WBASE=1.0 
!         FICE=0.0
!         NICE=0.0
!         NDROP=0.0
!         NDROP_AMB=0.0
!         NDROP_ACT=0.0
!         NIN = 0.0
!         NDUST=0.0
!         NSOOT=0.0
!         NDUST_AMB=0.0
!         NSOOT_AMB = 0.0
!         dCVW2=0.0
!         QICE=0.0
!         QLIQ=0.0
!         FPICE = 0.0
!         INDUST=0.0
!         INSOOT= 0.0
!         QT= 0.0
!         FPRECIP=0.0   
!         FNDRIM = 0.0
!         RIMM = 0.0   
!         TminusTa = 0.0
!
!         f_seasalt = 0.0
!         aseasalt = 0.0
!
!         call init_Aer(AER_BASE)
!
!         !AER_CLOUD=============================
      alm = 0.
      IF (1. .GT. (qoi(k)/qst(k)-rhmn)/(rhmx-rhmn)) THEN
        trg = (qoi(k)/qst(k)-rhmn)/(rhmx-rhmn)
      ELSE
        trg = 1.
      END IF
      IF (0.0 .LT. (autorampb-sige(ic))/0.2) THEN
        y1 = (autorampb-sige(ic))/0.2
      ELSE
        y1 = 0.0
      END IF
      IF (1.0 .GT. y1) THEN
        f4 = y1
      ELSE
        f4 = 1.0
      END IF
! to 1 at SIG=AUTORAMPB-0.2
      IF (sige(ic) .GE. 0.5) THEN
        f5 = 1.0
      ELSE
        f5 = 1.0 - 2.*co_zdep*(0.5-sige(ic))
        IF (f5 .LT. 0.0) THEN
          f5 = 0.0
        ELSE
          f5 = f5
        END IF
      END IF
      IF (trg .LE. 1.0e-5) THEN
! TRIGGER  =========>>
!            RC(IC) = 7
        RETURN
      ELSE
!  RECOMPUTE SOUNDING UP TO DETRAINMENT LEVEL
        poi_c = poi
        qoi_c = qoi
        poi_c(k) = poi_c(k) + tpert(i)
        qoi_c(k) = qoi_c(k) + qpert(i)
        zet(k+1) = 0.
        sht(k+1) = cp*poi_c(k)*prj(k+1)
        DO l=k,ic,-1
          IF (qst(l)*rhmax .GT. qoi_c(l)) THEN
            qol(l) = qoi_c(l)
          ELSE
            qol(l) = qst(l)*rhmax
          END IF
          IF (0.000 .LT. qol(l)) THEN
            qol(l) = qol(l)
          ELSE
            qol(l) = 0.000
          END IF
          ssl(l) = cp*prj(l+1)*poi_c(l) + grav*zet(l+1)
          hol(l) = ssl(l) + qol(l)*alhl
          hst(l) = ssl(l) + qst(l)*alhl
          tem = poi_c(l)*(prj(l+1)-prj(l))*cpbg
          zet(l) = zet(l+1) + tem
          zol(l) = zet(l+1) + (prj(l+1)-prh(l))*poi_c(l)*cpbg
        END DO
        DO l=ic+1,k
          tem = (prj(l)-prh(l-1))/(prh(l)-prh(l-1))
          sht(l) = ssl(l-1) + tem*(ssl(l)-ssl(l-1))
          qht(l) = .5*(qol(l)+qol(l-1))
        END DO
! SMOOTH HSTAR W/ 1-2-1 Filter
        IF (smooth_hst) THEN
! save for later
          hstx = hst(ic)
          DO l=k-1,ic+1,-1
            hst(l) = 0.25*(hst(l+1)+hst(l-1)) + 0.5*hst(l)
          END DO
          DO l=ic,ic
            hst(l) = 0.5*hst(l+1) + 0.5*hst(l)
          END DO
        END IF
!  CALCULATE LAMBDA, ETA, AND WORKFUNCTION
        lambda_min = .2/mxdiam(i)
        lambda_max = .2/200.
!     LAMBDA_MIN = .2/(LAMBDA_FAC*DPTH_BL)
!     LAMBDA_MAX = .2/( MAX( LAMBMX_FAC*DPTH_BL , DIAMMN_MIN ) )
        IF (hol(k) .LE. hst(ic)) THEN
! CANNOT REACH IC LEVEL  ======>>
!            RC(IC) = 1
          RETURN
        ELSE
!  LAMBDA CALCULATION: MS-A18
          tem = (hst(ic)-hol(ic))*(zol(ic)-zet(ic+1))
          DO l=ic+1,k-1
            tem = tem + (hst(ic)-hol(l))*(zet(l)-zet(l+1))
          END DO
          IF (tem .LE. 0.0) THEN
! NO VALID LAMBDA  ============>>
!            RC(IC) = 2
            RETURN
          ELSE
            alm = (hol(k)-hst(ic))/tem
            IF (alm .GT. lambda_max) THEN
!            RC(IC) = 3
              RETURN
            ELSE
              toki = 1.0
              IF (alm .LT. lambda_min) toki = (alm/lambda_min)**2
!we can probably replace this by a actual distribution based on grid cell size
!RC(IC) = 6
!RETURN
!LAMBDSV(IC) = ALM
!  ETA CALCULATION: MS-A2
              DO l=ic+1,k
                eta(l) = 1.0 + alm*(zet(l)-zet(k))
              END DO
              eta(ic) = 1.0 + alm*(zol(ic)-zet(k))
!  WORKFUNCTION CALCULATION:  MS-A22
              wfn = 0.0
              hcc(k) = hol(k)
              DO l=k-1,ic+1,-1
                hcc(l) = hcc(l+1) + (eta(l)-eta(l+1))*hol(l)
                tem = hcc(l+1)*dpb(l) + hcc(l)*dpt(l)
                eht(l) = eta(l+1)*dpb(l) + eta(l)*dpt(l)
                wfn = wfn + (tem-eht(l)*hst(l))*gam(l)
              END DO
              hcc(ic) = hst(ic)*eta(ic)
              wfn = wfn + (hcc(ic+1)-hst(ic)*eta(ic+1))*gam(ic)*dpb(ic)
!  VERTICAL VELOCITY/KE CALCULATION (ADDED 12/2001 JTB)
              bk3(k) = 0.0
              bk2(k) = 0.0
              bke(k) = 0.0
              hcld(k) = hol(k)
              DO l=k-1,ic,-1
                hcld(l) = (eta(l+1)*hcld(l+1)+(eta(l)-eta(l+1))*hol(l))/&
&                 eta(l)
                tem = (hcld(l)-hst(l))*(zet(l)-zet(l+1))/(1.0+lbcp*dqq(l&
&                 ))
                IF (cldmicro .LE. 0.0) THEN
                  bke(l) = bke(l+1) + grav*tem/(cp*prj(l+1)*poi(l))
                  IF (tem .LT. 0.0) THEN
                    max1 = 0.0
                  ELSE
                    max1 = tem
                  END IF
                  bk2(l) = bk2(l+1) + grav*max1/(cp*prj(l+1)*poi(l))
                  IF (tem .GT. 0.0) THEN
                    min1 = 0.0
                  ELSE
                    min1 = tem
                  END IF
                  bk3(l) = bk3(l+1) + grav*min1/(cp*prj(l+1)*poi(l))
                  IF (bk2(l) .LT. 0.0) THEN
                    max2 = 0.0
                  ELSE
                    max2 = bk2(l)
                  END IF
                  cvw(l) = SQRT(2.0*max2)
                END IF
              END DO
! 1.0 / ( 5.0*ALM )
              cu_diam = 1000.
!   ALPHA CALCULATION 
              rasal2i = rasal2_2d(i)
              IF (zet(ic) .LT. 2000.) rasal(ic) = rasal1
              IF (zet(ic) .GE. 2000.) THEN
                IF (1.0 .GT. (zet(ic)-2000.)/8000.) THEN
                  min2 = (zet(ic)-2000.)/8000.
                ELSE
                  min2 = 1.0
                END IF
!WMP     RASAL(IC) = RASAL1 + (RASAL2i-RASAL1)*(ZET(IC) - 2000.)/8000.
                rasal(ic) = rasal1 + (rasal2i-rasal1)*min2**rasal_exp
              END IF
!WMP  RASAL(IC) = MIN( RASAL(IC) , 1.0e5 )
              rasal(ic) = dt/rasal(ic)
              mask(1:k-ic+1) = cvw(ic:k) .LT. 1.00
              WHERE (mask(1:k-ic+1)) 
                WHERE (mask(1:k-ic+1)) 
                  cvw(ic:k) = 1.00
                ELSEWHERE
                  cvw(ic:k) = cvw(ic:k)
                END WHERE
              END WHERE
!  NOTE THIS "CENTRALIZES" A KLUGE PRESENT IN OTHER LOCATIONS.
!  CLEAN UP SOME TIME.      -JTB 12/04/03
!  TEST FOR CRITICAL WORK FUNCTION
              CALL ACRITN(pol(ic), prs(k), acr)
              IF (wfn .LE. acr) THEN
! SUB-CRITICAL WORK FUNCTION ======>>
!            RC(IC) = 4
                RETURN
              ELSE
!  CLOUD TOP WATER AND MOMENTUM (TIMES ETA(IC)) MS-A16
! Tracer scavenging
! RAS loops over a series of plumes all having common cloud base level K 
! and different detrainment levels IC.  The plumes operate sequentially
! on the grid box mean quantities (wind, moisture, tracer) and so each
! subsequent plume is seeing the effects of previous plumes.  We parameterize
! scavenging following Liu et al. [JGR, 2001], their equation 1:
!  AEROSOL FRACTION SCAVENGED = 1 - exp(-FSCAV*DZ)
! where FSCAV is a specified scavenging efficiency [km-1] and DZ is the
! distance [km] the tracer traverses in the plume from it's entrainment
! level to its detrainment level.  We write the aerosol fraction surviving as:
!  FNOSCAV = exp(- FSCAV_(ITR) * DZ)
! The total scavenging is proportional to the convective mass flux, which
! is not explicitly solved for at this point.
                IF (do_tracers) THEN
                  DO itr=1,itrcr
!           Scavenging of the below cloud tracer
                    delzkm = (zet(ic)-zet(k))/1000.
                    x4 = EXP(-(fscav_(itr)*delzkm))
                    IF (x4 .GT. 1.) THEN
                      x1 = 1.
                    ELSE
                      x1 = x4
                    END IF
                    IF (x1 .LT. 0.) THEN
                      fnoscav = 0.
                    ELSE
                      fnoscav = x1
                    END IF
                    xht(itr) = xoi(k, itr)*fnoscav
                  END DO
                END IF
                wlq = qol(k)
                uht = uoi(k)
                vht = voi(k)
                rnn(k) = 0.
                cll0(k) = 0.
!print *, '========================================='
                DO l=k-1,ic,-1
                  tem = eta(l) - eta(l+1)
                  wlq = wlq + tem*qol(l)
                  uht = uht + tem*uoi(l)
                  vht = vht + tem*voi(l)
                  IF (do_tracers) THEN
                    DO itr=1,itrcr
!         Scavenging of the entrained tracer.  Updates transported tracer mass.
                      delzkm = (zet(ic)-zet(l+1))/1000.
                      x5 = EXP(-(fscav_(itr)*delzkm))
                      IF (x5 .GT. 1.) THEN
                        x2 = 1.
                      ELSE
                        x2 = x5
                      END IF
                      IF (x2 .LT. 0.) THEN
                        fnoscav = 0.
                      ELSE
                        fnoscav = x2
                      END IF
                      xht(itr) = xht(itr) + tem*xoi(l, itr)*fnoscav
                    END DO
                  END IF
!!!! How much condensate (CLI) is present here? 
                  IF (l .GT. ic) THEN
                    tx2 = 0.5*(qst(l)+qst(l-1))*eta(l)
                    tx3 = 0.5*(hst(l)+hst(l-1))*eta(l)
                    qcc = tx2 + gm1(l)*(hcc(l)-tx3)
                    cll0(l) = wlq - qcc
                  ELSE
                    cll0(l) = wlq - qst(ic)*eta(ic)
                  END IF
                  IF (cll0(l) .LT. 0.00) THEN
                    cll0(l) = 0.00
                  ELSE
                    cll0(l) = cll0(l)
                  END IF
! condensate (kg/kg)
                  cli = cll0(l)/eta(l)
! Temperature (K)
                  te_a = poi(l)*prh(l)
!=====================================================================
                  IF (cldmicro .LE. 0.0) THEN
!AER_CLOUD MIcrophysics considering activation and nucleation 
!               !recompute vertical velocity
!
!               Tparcel = TE_A
!               CVW(K) = 0.8  ! Assume a below cloud base  W of 0.8 m s-1            
!               BK2(K)   = 0.0
!
!     
!               TEM     = (HCLD(L)-HST(L) )/ (1.0+LBCP*DQQ(L))  
!               TminusTa = max(min(TEM/CP, 5.0), 0.0) !limit DT to 5 K. According to Wei, JAS, 1998   
!	     TEM =0.33*TminusTa*CO_AUTO(I)/TE_A !Bouyancy term, effciency =0.5 mwr Roode et al    	     
!
!               BK2(L)  = BK2(L+1) + GRAV * TEM*(ZET(L)-ZET(L+1)) 
!               BK2(L) = BK2(L) - (ZET(L)-ZET(L+1))*(BK2(L+1)*ALM + CLI*GRAV)  !Account for drag from entrainment of stagnat air 
!and condesate loading
!               CVW(L) = max(SQRT(  2.0* MAX( BK2(L) , 0.0 )  ), 1.0) 
!
!
!	    CVW_X = MIN(CVW(L), 50.0)
!               DT_LYR  =  max(( ZET(L)-ZET(L+1) )/CVW_X, 1.0) !Sanity check 
!               TEM   = ETA(L) - ETA(L+1)
!
!               Tparcel  =  TE_A + TminusTa
!
!
!
!!!!!!!!!!account for entrainment effects on activation !!!!!!!!!!!
!               ! Barahona and Nenes, JGR, 2007
!               alph_e = 2.8915e-8*Tparcel*Tparcel -2.1328e-5*Tparcel+4.2523e-3
!               beta_e = MAPL_ALHL*TminusTa/MAPL_RVAP/Tparcel/Tparcel
!               RH_AMB=QOI(L)/QST(L)
!               ECRIT  = max(RH_AMB -beta_e, 1.0e-6) 
!               ECRIT =  alph_e/ECRIT
!               ! print *, L, Tparcel, RH_AMB, ECRIT, ALM
!	           ECRIT =  ALM/ECRIT
!       ! ECRIT = 0.0 ! do not use this for now
!               !Print *, ECRIT
!
!
!               if (L .eq. K-1) then
!
!                  FICE=0.0
!                  NICE=0.0
!                  NDROP=0.0
!                  NIN =0.0
!                  NDUST_AMB =0.0
!                  NSOOT_AMB = 0.0
!                  NSOOT=0.0
!                  NDUST= 0.0
!
!        
!                  AER_BASE =  AERO(L)
!                  RATE=0.0
!                  FPRECIP=0.0
!                  !initial conditions     
!                  call ARGact(Tparcel, CVW_X, NDROP_ACT, NDROP_AMB, NDUST_AMB, NSOOT_AMB, L,  .true., DDUST_AMB, DSOOT_AMB, ECRI
!T) !cloud droplet number and INsource at cloud base 
!                  NDUST=NDUST_AMB
!                  NSOOT=NSOOT_AMB
!                  DDUST=DDUST_AMB
!                  DSOOT=DSOOT_AMB                                     
!
!               else 
!                  call ARGact(Tparcel, CVW_X, NDROP_ACT, NDROP_AMB, NDUST_AMB, NSOOT_AMB, L, .false., DDUST_AMB, DSOOT_AMB, ECRI
!T) !cloud droplet number above cloud base  
!
!               end if
!
!               QT = CLI
!               RATE = 0.0
!               FPRECIP = 0.0
!
!               if (QT .gt. 0.0) then
!
!                  ! if FICE is already >= 1.0 then the cloud is glaciated and there is no need to do anymore partitioning
!
!                  if (FICE .ge. 1.0) then
!
!
!                     CALL  Qremoval(RATE, FICE, FP_D, FP_I, Tparcel,  & 
!                           POL(L), QT,  NICE, NDROP, CVW_X, FPICE, &
!                           DT_LYR, RIMM, CO_AUTO(I)) 
!
!
!
!                     dNICE = -NICE*FP_I 
!                     NICE  =  (NICE +dNICE*DT_LYR)*ETA(L+1)/ETA(L) !ice
!
!                     MINNICE = max(QICE*(1.0-RATE*DT_LYR), 0.0)/4.0e-8!assuming maximum vol radius 250 microns
!                     MAXNICE = max(QICE*(1.0-RATE*DT_LYR), 0.0)/2.51e-12 !assuming minimum vol radius 10 microns
!
!                     NICE=MIN(max(NICE, MINNICE), MAXNICE)
!
!                     FICE = 1.0
!
!                  else 
!
!                     ! Cloud is not completely glaciated do the whole thing
!                     ! ALL this subroutines return tendencies
!
!
!                     CALL  INfreezing(QLIQ, NDROP, NIN, NDUST, NSOOT, INDUST, INSOOT, Tparcel, POL(L), CVW_X, DDUST, DSOOT)  !ca
!lculate the freezing fraction of the aerosol at this level
!
!                     NIN = min(NIN, NDROP/DT_LYR)
!
!                     call Qgrowth(Tparcel, POL(L), QICE, NICE, QT, NIN, dQIG, RIMM, FNDRIM)
!
!                     CALL  Qremoval(RATE, FICE, FP_D, FP_I, Tparcel,  & 
!                           POL(L), QT,  NICE, NDROP, CVW_X, FPICE, &
!                           DT_LYR,  RIMM, CO_AUTO(I)) 
!
!
!
!                     !ice number tendency: -precip + freezin
!                     dNICE = -NICE*FP_I  + NIN     
!                     NICE  =  (NICE +dNICE*DT_LYR)*ETA(L+1)/ETA(L) !ice
!                     NICE =max(NICE, 0.0)
!
!
!                     !ice mass tendency: growth - precip
!                     dQICE = -QICE*FPICE + dQIG
!                     QICE  =  min((QICE + dQICE*DT_LYR)*ETA(L+1)/ETA(L), QT) !ice
!                     QICE=max(min(QT, QICE), 0.0)
!
!
!                     ! Liquid Tendency: source/evap -  precip 
!                     !dQLIQ = max((CLI-QICE), -QLIQ)/DT_LYR -QLIQ*max(RATE-FPICE, 0.0) 
!                     ! dQLIQ = CLI*(1.0-RATE*DT_LYR)/DT_LYR -dQICE - QLIQ*max(RATE-FPICE, 0.0)           
!                     !QLIQ  =  max((QLIQ + dQLIQ*DT_LYR)*ETA(L+1)/ETA(L), 0.0) !liquid. This is actually diagnostic
!                     QLIQ=max((QT-QICE), 0.0)
!
!
!                     !droplet number tendency: -precip - freezin + activation + activated entrained aerosol 
!
!
!                     dNDROP =-NDROP*FP_D - NIN -  FNDRIM*NDROP/DT_LYR + max(NDROP_ACT-NDROP, 0.0)/DT_LYR          
!
!                     !dNDROP =-NDROP*FP_D - NIN -  FNDRIM*NDROP/DT_LYR + NDROP_ACT/DT_LYR         
!
!                     NDROP =  (NDROP + dNDROP*DT_LYR)*ETA(L+1)/ETA(L) + &
!                           (ZET(L) - ZET(L+1))*ALM*MAX((NDROP_AMB-NDROP), 0.0)
!
!                     !Aerosol tendency: Entrainment - freezing 
!
!                     NDUST = (NDUST - INDUST*DT_LYR)*ETA(L+1)/ETA(L) + &
!                           (ZET(L) - ZET(L+1))*ALM*MAX(NDUST_AMB-NDUST, 0.0) 
!
!                     NSOOT =  (NSOOT - INSOOT*DT_LYR)*ETA(L+1)/ETA(L)  + &     
!                           (ZET(L) - ZET(L+1))*ALM*MAX(NSOOT_AMB-NSOOT, 0.0)  
!
!
!                           
!                     !Update FICE and perform Sanity checks
!
!
!                     MINNDROP = (1.0-FICE)*QLIQ*max(1.0-RATE*DT_LYR, 0.0)/2.e-10    !assuming maximum vol radius 36 microns
!                     MAXNDROP = (1.0-FICE)*QLIQ*max(1.0-RATE*DT_LYR, 0.0)/3.35e-14 !assuming minimum vol radius 2 microns
!                     MINNICE = QICE/4.0e-8!assuming maximum vol radius 250 microns
!                     MAXNICE = QICE/2.51e-12 !assuming minimum vol radius 10 microns
!
!                     IF ((NICE .gt. MAXNICE) .or. (NICE .lt. MINNICE))   then    
!                        !print *, 'nilim', NICE*1e-6, MINNICE*1e-6, MAXNICE*1e-6
!                     END IF
!
!                     IF ((NDROP .gt. MAXNDROP) .or. (NDROP .lt. MINNDROP))      then 
!                        !print *, 'ndroplim', NDROP*1e-6, MINNDROP*1e-6, MAXNDROP*1e-6
!                     end if
!
!
!                     NSOOT=MAX(NSOOT, 0.0)
!                     NDUST=MAX(NDUST, 0.0)              
!
!                     NDROP=MIN(max(NDROP, MINNDROP), MAXNDROP)
!                     NICE=MIN(max(NICE, MINNICE), MAXNICE)
!
!                     FICE=max(min(QICE/QT, 1.0), 0.0)
!
!                     IF (FICE .ge. 1.0) THEN !Complete glaciation 
!                        NICE=NICE+NDROP 
!                        NDROP = 0.0
!                        QICE  = QT
!                        QLIQ= 0.0
!                     END IF
!
!                     IF (Tparcel .LT. T_ICE_ALL) THEN !instantaneous freezing
!                        NICE=NICE+NDROP 
!                        NDROP = 0.0
!                        FICE  = 1.0
!                        QICE  = QT
!                        QLIQ=0.0
!                     END IF
!
!                     IF (Tparcel .GT. T_ICE_MAX) THEN !instantaneous melting
!                        NDROP=NICE+NDROP 
!                        NICE = 0.0
!                        FICE  = 0.0
!                        QICE  = 0.0
!                        QLIQ=QT
!                     END IF
!
!                  END IF
!
!               else 
!
!                  FICE =0.0 
!                  QICE = 0.0
!                  QLIQ = 0.0
!                  NICE= 0.0 
!                  NDROP = 0.0
!                  RATE =0.0
!               end if
!
!               FPRECIP= RATE*DT_LYR
!
!               !RATE=RATE*F4
!               ! NDROP=NDROP*F4
!               !NICE=NICE*(1.0-F4)
!
!               !print *, TE_A, FICE, 'NICE', NICE*1e-6, 'NDROP', NDROP*1e-6, L 
!               !print *, 'FPI', FP_I*DT_LYR, 'FPD', FP_D*DT_LYR, 'FPICE', FPICE, 'FPRE', FPRECIP, QT, QLIQ
!
!Bacmeister 2006 microphysics
                    CALL SUNDQ3_ICE(te_a, sdqv2, sdqv3, sdqvt1, f2, f3)
! * F5  ! F4 reduces AUTO for shallow clouds, F5 modifies auto for deep clouds
                    c00_x = co_auto(i)*f2*f3*f4
                    cli_crit_x = cli_crit/(f2*f3)
                    rate = c00_x*(1.0-EXP(-(cli**2/cli_crit_x**2)))
                  END IF
                  IF (cvw(l) .LT. 1.00) THEN
                    cvw_x = 1.00
                  ELSE
                    cvw_x = cvw(l)
                  END IF
! really trust it at low values
! l.h.s. DT_LYR => time in layer (L,L+1)
                  dt_lyr = (zet(l)-zet(l+1))/cvw_x
                  closs = cll0(l)*rate*dt_lyr
                  IF (closs .GT. cll0(l)) THEN
                    closs = cll0(l)
                  ELSE
                    closs = closs
                  END IF
                  cll0(l) = cll0(l) - closs
                  dll0(l) = closs
                  IF (closs .GT. 0.) THEN
                    wlq = wlq - closs
                    rnn(l) = closs
                  ELSE
                    rnn(l) = 0.
                  END IF
                END DO
!AER_CLOUD=======================================
!            CNVNDROP(IC)=NDROP
!            CNVNICE(IC)=NICE
!            CNVFICE(IC)=FICE
                wlq = wlq - qst(ic)*eta(ic)
!     CALCULATE GAMMAS AND KERNEL
! MS-A30 (W/O GRAV)
                gms(k) = (sht(k)-ssl(k))*pri(k)
! MS-A31 (W/O GRAV)
                gmh(k) = gms(k) + (qht(k)-qol(k))*pri(k)*alhl
! MS-A37 (W/O GRAV)
                akm = gmh(k)*gam(k-1)*dpb(k-1)
                tx2 = gmh(k)
                DO l=k-1,ic+1,-1
                  gms(l) = (eta(l)*(sht(l)-ssl(l))+eta(l+1)*(ssl(l)-sht(&
&                   l+1)))*pri(l)
                  gmh(l) = gms(l) + (eta(l)*(qht(l)-qol(l))+eta(l+1)*(&
&                   qol(l)-qht(l+1)))*alhl*pri(l)
                  tx2 = tx2 + (eta(l)-eta(l+1))*gmh(l)
                  akm = akm - gms(l)*eht(l)*pki(l) + tx2*ght(l)
                END DO
                gms(ic) = eta(ic+1)*(ssl(ic)-sht(ic+1))*pri(ic)
                akm = akm - gms(ic)*eta(ic+1)*dpb(ic)*pki(ic)
                gmh(ic) = gms(ic) + (eta(ic+1)*(qol(ic)-qht(ic+1))*alhl+&
&                 eta(ic)*(hst(ic)-hol(ic)))*pri(ic)
                IF (smooth_hst) gmhx = gms(ic) + (eta(ic+1)*(qol(ic)-qht&
&                   (ic+1))*alhl+eta(ic)*(hstx-hol(ic)))*pri(ic)
!    CLOUD BASE MASS FLUX
                IF (akm .GE. 0.0 .OR. wlq .LT. 0.0) THEN
!  =========>
!            RC(IC) = 5
                  RETURN
                ELSE
! MS-A39 MASS-FLUX IN Pa/step
                  wfn = -((wfn-acr)/akm)
                  x3 = rasal(ic)*trg*toki*wfn
                  IF (x3 .GT. (prs(k+1)-prs(k))*(100.*pblfrac)) THEN
                    wfn = (prs(k+1)-prs(k))*(100.*pblfrac)
                  ELSE
                    wfn = x3
                  END IF
!    CUMULATIVE PRECIP AND CLOUD-BASE MASS FLUX FOR OUTPUT
                  wfnog = wfn*gravi
                  tem = wfn*gravi
! (kg/m^2/step)
                  cll(ic) = cll(ic) + wlq*tem
! (kg/m^2/step)
                  rmf(ic) = rmf(ic) + tem
! (kg/m^2/step)
                  rmfd(ic) = rmfd(ic) + tem*eta(ic)
                  DO l=ic+1,k
! (kg/m^2/step)
                    rmfp(l) = tem*eta(l)
! (kg/m^2/step)
                    rmfc(l) = rmfc(l) + rmfp(l)
                    dllx(l) = dllx(l) + tem*dll0(l)
                    IF (cvw(l) .GT. 0.0) THEN
                      updfrp(l) = rmfp(l)*(ddt/daylen)*1000./(cvw(l)*prs&
&                       (l))
                    ELSE
                      updfrp(l) = 0.0
                    END IF
! current cloud; incloud condensate        
                    clli(l) = cll0(l)/eta(l)
!  cumulative grid mean convective condensate        
                    cllb(l) = cllb(l) + updfrp(l)*clli(l)
                    updfrc(l) = updfrc(l) + updfrp(l)
                  END DO
!    THETA AND Q CHANGE DUE TO CLOUD TYPE IC
                  DO l=ic,k
! (kg/m^2/step)
                    rns(l) = rns(l) + rnn(l)*tem
                    gmh(l) = gmh(l)*wfn
                    gms(l) = gms(l)*wfn
                    qoi(l) = qoi(l) + (gmh(l)-gms(l))*alhi
                    poi(l) = poi(l) + gms(l)*pki(l)*cpi
                    qst(l) = qst(l) + gms(l)*bet(l)*cpi
                  END DO
                  IF (smooth_hst) THEN
                    gmhx = gmhx*wfn
                    dqx = (gmhx-gmh(ic))*alhi
                    rns(ic) = rns(ic) + dqx/(pri(ic)*grav)
                  END IF
                  IF (do_tracers) THEN
!*FRICFAC*0.5
                    wfn = wfn*0.5*1.0
                    tem = wfn*pri(k)
                    DO itr=1,itrcr
                      xcu(k, itr) = xcu(k, itr) + tem*(xoi(k-1, itr)-xoi&
&                       (k, itr))
                    END DO
                    DO itr=1,itrcr
                      DO l=k-1,ic+1,-1
                        tem = wfn*pri(l)
                        xcu(l, itr) = xcu(l, itr) + tem*((xoi(l-1, itr)-&
&                         xoi(l, itr))*eta(l)+(xoi(l, itr)-xoi(l+1, itr)&
&                         )*eta(l+1))
                      END DO
                    END DO
                    tem = wfn*pri(ic)
                    DO itr=1,itrcr
                      xcu(ic, itr) = xcu(ic, itr) + (2.*(xht(itr)-xoi(ic&
&                       , itr)*(eta(ic)-eta(ic+1)))-(xoi(ic, itr)+xoi(ic&
&                       +1, itr))*eta(ic+1))*tem
                    END DO
                    DO itr=1,itrcr
                      DO l=ic,k
                        xoi(l, itr) = xoi(l, itr) + xcu(l, itr)
                      END DO
                    END DO
                  ELSE
!*FRICFAC*0.5
                    wfn = wfn*0.5*1.0
                  END IF
                  lambdsv(ic) = 1.000
!   CUMULUS FRICTION
                  IF (fricfac .LE. 0.0) THEN
!            RC(IC) = 0
!  NO CUMULUS FRICTION =========>>
                    RETURN
                  ELSE
                    wfn = wfn*fricfac*EXP(-(alm/friclambda))
                    tem = wfn*pri(k)
                    ucu(k) = ucu(k) + tem*(uoi(k-1)-uoi(k))
                    vcu(k) = vcu(k) + tem*(voi(k-1)-voi(k))
                    DO l=k-1,ic+1,-1
                      tem = wfn*pri(l)
                      ucu(l) = ucu(l) + tem*((uoi(l-1)-uoi(l))*eta(l)+(&
&                       uoi(l)-uoi(l+1))*eta(l+1))
                      vcu(l) = vcu(l) + tem*((voi(l-1)-voi(l))*eta(l)+(&
&                       voi(l)-voi(l+1))*eta(l+1))
                    END DO
                    tem = wfn*pri(ic)
                    ucu(ic) = ucu(ic) + (2.*(uht-uoi(ic)*(eta(ic)-eta(ic&
&                     +1)))-(uoi(ic)+uoi(ic+1))*eta(ic+1))*tem
                    vcu(ic) = vcu(ic) + (2.*(vht-voi(ic)*(eta(ic)-eta(ic&
&                     +1)))-(voi(ic)+voi(ic+1))*eta(ic+1))*tem
                    dissk0(ic) = eta(ic)*grav*wfnog*pri(ic)*0.5*((uht/&
&                     eta(ic)-uoi(ic))**2+(vht/eta(ic)-voi(ic))**2)
                    DO l=ic,k
                      uoi(l) = uoi(l) + ucu(l)
                      voi(l) = voi(l) + vcu(l)
                    END DO
!         RC(IC) = 0
                    RETURN
                  END IF
                END IF
              END IF
            END IF
          END IF
        END IF
      END IF
    END SUBROUTINE CLOUDE
    SUBROUTINE ACRITN(pl, plb, acr)
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: pl, plb
      REAL*8, INTENT(OUT) :: acr
      INTEGER :: iwk
!!REAL(8), PARAMETER :: FACM=0.5
      REAL*8, PARAMETER :: ph(15)=(/150.0, 200.0, 250.0, 300.0, 350.0, &
&       400.0, 450.0, 500.0, 550.0, 600.0, 650.0, 700.0, 750.0, 800.0, &
&       850.0/)
!!*FACM
      REAL*8, PARAMETER :: a(15)=(/1.6851, 1.1686, 0.7663, 0.5255, &
&       0.4100, 0.3677, 0.3151, 0.2216, 0.1521, 0.1082, 0.0750, 0.0664, &
&       0.0553, 0.0445, 0.0633/)
      INTRINSIC INT
      iwk = INT(pl*0.02 - 0.999999999)
      IF (iwk .GT. 1 .AND. iwk .LE. 15) THEN
        acr = a(iwk-1) + (pl-ph(iwk-1))*.02*(a(iwk)-a(iwk-1))
      ELSE IF (iwk .GT. 15) THEN
        acr = a(15)
      ELSE
        acr = a(1)
      END IF
      acr = acritfac*acr*(plb-pl)
      RETURN
    END SUBROUTINE ACRITN
    SUBROUTINE RNEVP()
      IMPLICIT NONE
      zet(k+1) = 0
      DO l=k,icmin,-1
        tem = poi(l)*(prj(l+1)-prj(l))*cpbg
        zet(l) = zet(l+1) + tem
      END DO
      DO l=icmin,k
        tem = pri(l)*grav
        cnv_prc3(i, l) = rns(l)*tem
      END DO
!! If hst is smoothed then adjusted precips may be negative
      IF (smooth_hst) THEN
        DO l=icmin,k
          IF (cnv_prc3(i, l) .LT. 0.) THEN
            qoi(l) = qoi(l) + cnv_prc3(i, l)
            poi(l) = poi(l) - cnv_prc3(i, l)*(alhl/cp)/prj(l+1)
            cnv_prc3(i, l) = 0.
          END IF
        END DO
      END IF
      RETURN
    END SUBROUTINE RNEVP
    SUBROUTINE HTEST()
      IMPLICIT NONE
      REAL*8, DIMENSION(k0) :: hol1
      INTEGER :: lminhol
      REAL*8 :: minhol
      INTRINSIC MIN
      INTRINSIC MAX
! HOL initialized here in order not to confuse Valgrind debugger
      hol = 0.
      lminhol = k + 1
      minhol = -999999.
      zet(k+1) = 0
      sht(k+1) = cp*poi(k)*prj(k+1)
      DO l=k,icmin,-1
        IF (qst(l)*rhmax .GT. qoi(l)) THEN
          qol(l) = qoi(l)
        ELSE
          qol(l) = qst(l)*rhmax
        END IF
        IF (0.000 .LT. qol(l)) THEN
          qol(l) = qol(l)
        ELSE
          qol(l) = 0.000
        END IF
        ssl(l) = cp*prj(l+1)*poi(l) + grav*zet(l+1)
        hol(l) = ssl(l) + qol(l)*alhl
        hst(l) = ssl(l) + qst(l)*alhl
        tem = poi(l)*(prj(l+1)-prj(l))*cpbg
        zet(l) = zet(l+1) + tem
        zol(l) = zet(l+1) + (prj(l+1)-prh(l))*poi(l)*cpbg
      END DO
      hol1 = hol
      DO l=k-1,icmin+1,-1
        hol1(l) = 0.25*hol(l+1) + 0.50*hol(l) + 0.25*hol(l-1)
        IF (minhol .GE. hol1(l) .OR. minhol .LT. 0.) THEN
          minhol = hol1(l)
          lminhol = l
        END IF
      END DO
      sige_minhol = sige(lminhol)
    END SUBROUTINE HTEST
    SUBROUTINE FINDDTLS()
      IMPLICIT NONE
      REAL*8 :: sigdt0, sigmax, sigmin
      INTEGER :: ll
      INTRINSIC ALLOCATED
!#ifndef __GFORTRAN__
!         integer :: THE_SEED(2)
!#else
!         integer :: THE_SEED(12)
!#endif
!         THE_SEED(1)=SEEDRAS(I,1)*IRAS(I) + SEEDRAS(I,2)*JRAS(I)
!         THE_SEED(2)=SEEDRAS(I,1)*JRAS(I) + SEEDRAS(I,2)*IRAS(I)
!         THE_SEED(1)=THE_SEED(1)*SEEDRAS(I,1)/( SEEDRAS(I,2) + 10)
!         THE_SEED(2)=THE_SEED(2)*SEEDRAS(I,1)/( SEEDRAS(I,2) + 10)
!         if(THE_SEED(1) == 0) THE_SEED(1) =  5
!         if(THE_SEED(2) == 0) THE_SEED(2) = -5
!#ifdef __GFORTRAN__
!         THE_SEED(3:12) = 0
!#endif
!
!         call random_seed(PUT=THE_SEED)
!dh RASNCL = -300
      sigmax = sige(k)
      sigmin = sige(icmin)
      IF (rasncl .LT. 0.0) n_dtl = k - icmin
!! NO SHALLOW CONV   N_DTL = 56 - ICMIN 
!            N_DTL = min( int( RASNCL ) , K-ICMIN )
      IF (ALLOCATED(icl_v)) THEN
        DEALLOCATE(icl_v)
      END IF
      ALLOCATE(icl_v(n_dtl))
!         if ( ( RASNCL < 0.0 ) .and. ( RASNCL >=-100.) ) then 
!            do L=1,N_DTL
!               ICL_V(L) = ICMIN + L - 1
!            enddo
!         else if ( RASNCL < -100.0 ) then 
      DO l=1,n_dtl
        icl_v(l) = k - l
      END DO
    END SUBROUTINE FINDDTLS
    SUBROUTINE STRAP(final)
      IMPLICIT NONE
      INTEGER :: final
      REAL*8, DIMENSION(k0) :: wght, massf
      REAL*8 :: wght0, prcbl
      INTEGER, PARAMETER :: nrands=1
      REAL*8 :: rndu(nrands)
      INTEGER :: seedcbl(nrands), jj
! !DESCRIPTION: 
!   {\tt STRAP} is called: FINAL=0, to compute cloud base layer CBL properties
!   given a value K for the index of the upper {\em EDGE} of the CBL; FINAL=1
!   to redistribute convective tendencies within CBL
      INTEGER :: kk
      INTRINSIC SQRT
      INTRINSIC MAX
      INTRINSIC ABS
      INTRINSIC PRESENT
      INTRINSIC ALLOCATED
      REAL*8 :: abs1
      REAL*8 :: abs0
!  LOCAL VARIABLES FOR USE IN CLOUDE
!!IF (.NOT. PRESENT(FINAL)) THEN
      IF (final .EQ. 0) THEN
!!PRJ(ICMIN:K+1) = PKE(I,ICMIN:K+1)
        DO kk=icmin,k+1
          prj(kk) = pke(i, kk)
        END DO
! These initialized here in order not to confuse Valgrind debugger
        poi = 0.
! Do not believe it actually makes any difference.
        qoi = 0.
        uoi = 0.
        voi = 0.
        prs(icmin:k0+1) = ple(i, icmin:k0+1)
        poi(icmin:k) = tho(i, icmin:k)
        qoi(icmin:k) = qho(i, icmin:k)
        uoi(icmin:k) = uho(i, icmin:k)
        voi(icmin:k) = vho(i, icmin:k)
        wsp(icmin:k) = SQRT((uoi(icmin:k)-uoi(k))**2 + (voi(icmin:k)-voi&
&         (k))**2)
        qst(icmin:k) = qss(i, icmin:k)
        dqq(icmin:k) = dqs(i, icmin:k)
        IF (do_tracers) THEN
          DO itr=1,itrcr
            xoi(icmin:k, itr) = xho(i, icmin:k, itr)
          END DO
        END IF
!!! Mass fraction of each layer below cloud base
!!! contributed to aggregate cloudbase layer (CBL) 
        massf(:) = wgt0(i, :)
!!! RESET PRESSURE at bottom edge of CBL 
        prcbl = prs(k)
        DO l=k,k0
          prcbl = prcbl + massf(l)*(prs(l+1)-prs(l))
        END DO
        prs(k+1) = prcbl
        prj(k+1) = (prs(k+1)/1000.)**(mapl_rgas/mapl_cp)
        DO l=k,icmin,-1
          pol(l) = 0.5*(prs(l)+prs(l+1))
          prh(l) = (prs(l+1)*prj(l+1)-prs(l)*prj(l))/(onepkap*(prs(l+1)-&
&           prs(l)))
          pki(l) = 1.0/prh(l)
          dpt(l) = prh(l) - prj(l)
          dpb(l) = prj(l+1) - prh(l)
          pri(l) = .01/(prs(l+1)-prs(l))
        END DO
!!!!! RECALCULATE PROFILE QUAN. IN LOWEST STRAPPED LAYER
        IF (k .LE. k0) THEN
          poi(k) = 0.
          qoi(k) = 0.
          uoi(k) = 0.
          voi(k) = 0.
!! SPECIFY WEIGHTS GIVEN TO EACH LAYER WITHIN SUBCLOUD "SUPERLAYER"
          wght = 0.
          DO l=k,k0
            wght(l) = massf(l)*(ple(i, l+1)-ple(i, l))/(prs(k+1)-prs(k))
          END DO
          DO l=k,k0
            poi(k) = poi(k) + wght(l)*tho(i, l)
            qoi(k) = qoi(k) + wght(l)*qho(i, l)
            uoi(k) = uoi(k) + wght(l)*uho(i, l)
            voi(k) = voi(k) + wght(l)*vho(i, l)
          END DO
          IF (do_tracers) THEN
            xoi(k, :) = 0.
            DO itr=1,itrcr
              DO l=k,k0
                xoi(k, itr) = xoi(k, itr) + wght(l)*xho(i, l, itr)
              END DO
            END DO
          END IF
          tmpt(k) = poi(k)*prh(k)
          CALL DQSATPERT(dqq(k), qst(k), tmpt(k), pol(k), 1)
        END IF
!!DPTH_BL = CPBG*POI(K)*( PRJ(K+1)-PRJ(K) )
! seedras(1,2) are both integers passed from
! from GEOS_Moist w/ values 0 - 1000000
! rndu(:) = 1.0*( seedras(1)+seedras(2) )/2000000.
!dh tap produces code with error
        DO jj=1,nrands
          IF (seedras(i, 1)/1000000. .LT. 1e-6) THEN
            rndu(jj) = 1e-6
          ELSE
            rndu(jj) = seedras(i, 1)/1000000.
          END IF
        END DO
        IF (maxdallowed_d .GE. 0.) THEN
          abs0 = maxdallowed_d
        ELSE
          abs0 = -maxdallowed_d
        END IF
        IF (maxdallowed_s .GE. 0.) THEN
          abs1 = maxdallowed_s
        ELSE
          abs1 = -maxdallowed_s
        END IF
!!call congvec( npoints , seedcbl , rndu )
!            DPTH_BL   = ZCBL(I)
        mxdiam(i) = cnv_fraction(i)*abs0 + (1-cnv_fraction(i))*abs1
        IF (maxdallowed_d .GT. 0) mxdiam(i) = mxdiam(i)*rndu(1)**(-(1./&
&           2.))
! Make MXDIAM stochastic
        DO l=k,icmin,-1
!*
          bet(l) = dqq(l)*pki(l)
!*
          gam(l) = pki(l)/(1.0+lbcp*dqq(l))
          IF (l .LT. k) THEN
            ght(l+1) = gam(l)*dpb(l) + gam(l+1)*dpt(l+1)
            gm1(l+1) = 0.5*lbcp*(dqq(l)/(alhl*(1.0+lbcp*dqq(l)))+dqq(l+1&
&             )/(alhl*(1.0+lbcp*dqq(l+1))))
          END IF
        END DO
        tcu(icmin:k) = -(poi(icmin:k)*prh(icmin:k))
        qcu(icmin:k) = -qoi(icmin:k)
        rns = 0.
        cll = 0.
        rmf = 0.
        rmfd = 0.
        rmfc = 0.
        rmfp = 0.
        cll0 = 0.
        dll0 = 0.
        cllx = 0.
        dllx = 0.
        clli = 0.
        cllb = 0.
        poi_sv = poi
        qoi_sv = qoi
        uoi_sv = uoi
        voi_sv = voi
        IF (do_tracers) xoi_sv = xoi
        lambdsv = 0.0
        cvw = 0.0
        updfrc = 0.0
        updfrp = 0.0
        dissk0 = 0.0
      END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!    IF (PRESENT(FINAL)) THEN
      IF (final .EQ. 1) THEN
        tho(i, icmin:k-1) = poi(icmin:k-1)
        qho(i, icmin:k-1) = qoi(icmin:k-1)
        uho(i, icmin:k-1) = uoi(icmin:k-1)
        vho(i, icmin:k-1) = voi(icmin:k-1)
        cnv_updfrc(i, icmin:k-1) = updfrc(icmin:k-1)
!======================AER_CLOUD=============
!               CNV_NDROP   (I,ICMIN:K-1)  =    CNVNDROP(ICMIN:K-1) !DONIF
!               CNV_NICE   (I,ICMIN:K-1)   =     CNVNICE(ICMIN:K-1) !DONIF
!               CNV_FICE   (I,ICMIN:K-1)   =     CNVFICE(ICMIN:K-1) !DONIF
!! De-strap tendencies from RAS
!! specify weighting "SHAPE"
        wght = wgt1(i, :)
!! Scale properly by layer masses
        wght0 = 0.
        DO l=k,k0
          wght0 = wght0 + wght(l)*(ple(i, l+1)-ple(i, l))
        END DO
        wght0 = (prs(k+1)-prs(k))/wght0
        wght = wght0*wght
        DO l=k,k0
          tho(i, l) = tho(i, l) + wght(l)*(poi(k)-poi_sv(k))
          qho(i, l) = qho(i, l) + wght(l)*(qoi(k)-qoi_sv(k))
          uho(i, l) = uho(i, l) + wght(l)*(uoi(k)-uoi_sv(k))
          vho(i, l) = vho(i, l) + wght(l)*(voi(k)-voi_sv(k))
        END DO
        IF (do_tracers) THEN
          xho(i, icmin:k-1, :) = xoi(icmin:k-1, :)
          DO itr=1,itrcr
            DO l=k,k0
              xho(i, l, itr) = xho(i, l, itr) + wght(l)*(xoi(k, itr)-&
&               xoi_sv(k, itr))
            END DO
          END DO
        END IF
!            FLX (I,ICMIN:K) = RMF (ICMIN:K) * DDT/DAYLEN  !  (KG/m^2/s @ CLOUD BASE)
!  (KG/m^2/s @ CLOUD TOP)
        flxd(i, icmin:k) = rmfd(icmin:k)*ddt/daylen
!            FLXC(I,ICMIN:K) = RMFC(ICMIN:K) * DDT/DAYLEN  !  (KG/m^2/s @ CLOUD TOP)
!  (KG/m^2/s )
        clw(i, icmin:k) = cll(icmin:k)*ddt/daylen
        IF (PRESENT(disske)) disske(i, icmin:k-1) = dissk0(icmin:k-1)*&
&           ddt/daylen
!            FLX (I,1:ICMIN-1) = 0.
        flxd(i, 1:icmin-1) = 0.
!            FLXC(I,1:ICMIN-1) = 0.
        clw(i, 1:icmin-1) = 0.
        IF (k .LT. k0) THEN
!               FLX (I,K:K0) = 0.
          flxd(i, k:k0) = 0.
!               FLXC(I,K:K0) = 0.
          clw(i, k:k0) = 0.
        END IF
!            IRC (I,ICMIN:K-1) = RC(ICMIN:K-1)
        IF (ALLOCATED(icl_v)) THEN
          DEALLOCATE(icl_v)
        END IF
      END IF
      IF (final .EQ. 2) THEN
!            FLX (I,:) = 0.
        flxd(i, :) = 0.
!            FLXC(I,:) = 0.
        clw(i, :) = 0.
!            IRC (I,ICMIN:K-1) = RC(ICMIN:K-1)
      END IF
      RETURN
    END SUBROUTINE STRAP
    SUBROUTINE STRAP_FWD(final)
      IMPLICIT NONE
!            CNV_CVW   (I,ICMIN:K-1)   =     CVW(ICMIN:K-1)
!            CNV_QC(I,ICMIN:K-1)       =  CLLB(ICMIN:K-1)
      INTEGER :: final
      REAL*8, DIMENSION(k0) :: wght, massf
      REAL*8 :: wght0, prcbl
      INTEGER, PARAMETER :: nrands=1
      REAL*8 :: rndu(nrands)
      INTEGER :: seedcbl(nrands), jj
! !DESCRIPTION: 
!   {\tt STRAP} is called: FINAL=0, to compute cloud base layer CBL properties
!   given a value K for the index of the upper {\em EDGE} of the CBL; FINAL=1
!   to redistribute convective tendencies within CBL
      INTEGER :: kk
      INTRINSIC SQRT
      INTRINSIC MAX
      INTRINSIC ABS
      INTRINSIC PRESENT
      INTRINSIC ALLOCATED
      REAL*8 :: abs1
      REAL*8 :: abs0
!  LOCAL VARIABLES FOR USE IN CLOUDE
!!IF (.NOT. PRESENT(FINAL)) THEN
      IF (final .EQ. 0) THEN
!!PRJ(ICMIN:K+1) = PKE(I,ICMIN:K+1)
        DO kk=icmin,k+1
          CALL PUSHREAL8(prj(kk))
          prj(kk) = pke(i, kk)
        END DO
! These initialized here in order not to confuse Valgrind debugger
        CALL PUSHREAL8ARRAY(poi, k0)
        poi = 0.
! Do not believe it actually makes any difference.
        CALL PUSHREAL8ARRAY(qoi, k0)
        qoi = 0.
        CALL PUSHREAL8ARRAY(uoi, k0)
        uoi = 0.
        CALL PUSHREAL8ARRAY(voi, k0)
        voi = 0.
        CALL PUSHREAL8ARRAY(prs(icmin:k0+1), k0 - icmin + 2)
        prs(icmin:k0+1) = ple(i, icmin:k0+1)
        poi(icmin:k) = tho(i, icmin:k)
        qoi(icmin:k) = qho(i, icmin:k)
        uoi(icmin:k) = uho(i, icmin:k)
        voi(icmin:k) = vho(i, icmin:k)
        CALL PUSHREAL8ARRAY(qst(icmin:k), k - icmin + 1)
        qst(icmin:k) = qss(i, icmin:k)
        CALL PUSHREAL8ARRAY(dqq(icmin:k), k - icmin + 1)
        dqq(icmin:k) = dqs(i, icmin:k)
        IF (do_tracers) THEN
          DO itr=1,itrcr
            xoi(icmin:k, itr) = xho(i, icmin:k, itr)
          END DO
        END IF
!!! Mass fraction of each layer below cloud base
!!! contributed to aggregate cloudbase layer (CBL) 
        massf(:) = wgt0(i, :)
!!! RESET PRESSURE at bottom edge of CBL 
        prcbl = prs(k)
        DO l=k,k0
          prcbl = prcbl + massf(l)*(prs(l+1)-prs(l))
        END DO
        CALL PUSHREAL8(prs(k+1))
        prs(k+1) = prcbl
        CALL PUSHREAL8(prj(k+1))
        prj(k+1) = (prs(k+1)/1000.)**(mapl_rgas/mapl_cp)
        DO l=k,icmin,-1
          pol(l) = 0.5*(prs(l)+prs(l+1))
          CALL PUSHREAL8(prh(l))
          prh(l) = (prs(l+1)*prj(l+1)-prs(l)*prj(l))/(onepkap*(prs(l+1)-&
&           prs(l)))
          CALL PUSHREAL8(pki(l))
          pki(l) = 1.0/prh(l)
          CALL PUSHREAL8(dpt(l))
          dpt(l) = prh(l) - prj(l)
          CALL PUSHREAL8(dpb(l))
          dpb(l) = prj(l+1) - prh(l)
          CALL PUSHREAL8(pri(l))
          pri(l) = .01/(prs(l+1)-prs(l))
        END DO
!!!!! RECALCULATE PROFILE QUAN. IN LOWEST STRAPPED LAYER
        IF (k .LE. k0) THEN
          poi(k) = 0.
          qoi(k) = 0.
          uoi(k) = 0.
          voi(k) = 0.
!! SPECIFY WEIGHTS GIVEN TO EACH LAYER WITHIN SUBCLOUD "SUPERLAYER"
          wght = 0.
          DO l=k,k0
            wght(l) = massf(l)*(ple(i, l+1)-ple(i, l))/(prs(k+1)-prs(k))
          END DO
          DO l=k,k0
            poi(k) = poi(k) + wght(l)*tho(i, l)
            qoi(k) = qoi(k) + wght(l)*qho(i, l)
            uoi(k) = uoi(k) + wght(l)*uho(i, l)
            voi(k) = voi(k) + wght(l)*vho(i, l)
          END DO
          IF (do_tracers) THEN
            xoi(k, :) = 0.
            DO itr=1,itrcr
              DO l=k,k0
                xoi(k, itr) = xoi(k, itr) + wght(l)*xho(i, l, itr)
              END DO
            END DO
          END IF
          tmpt(k) = poi(k)*prh(k)
          CALL DQSATPERT_FWD(dqq(k), qst(k), tmpt(k), pol(k), 1)
          CALL PUSHCONTROL1B(1)
        ELSE
          CALL PUSHCONTROL1B(0)
        END IF
!!DPTH_BL = CPBG*POI(K)*( PRJ(K+1)-PRJ(K) )
! seedras(1,2) are both integers passed from
! from GEOS_Moist w/ values 0 - 1000000
! rndu(:) = 1.0*( seedras(1)+seedras(2) )/2000000.
!dh tap produces code with error
        DO jj=1,nrands
          IF (seedras(i, 1)/1000000. .LT. 1e-6) THEN
            CALL PUSHCONTROL1B(0)
            rndu(jj) = 1e-6
          ELSE
            CALL PUSHCONTROL1B(1)
            rndu(jj) = seedras(i, 1)/1000000.
          END IF
        END DO
        IF (maxdallowed_d .GE. 0.) THEN
          abs0 = maxdallowed_d
        ELSE
          abs0 = -maxdallowed_d
        END IF
        IF (maxdallowed_s .GE. 0.) THEN
          abs1 = maxdallowed_s
        ELSE
          abs1 = -maxdallowed_s
        END IF
!!call congvec( npoints , seedcbl , rndu )
!            DPTH_BL   = ZCBL(I)
        mxdiam(i) = cnv_fraction(i)*abs0 + (1-cnv_fraction(i))*abs1
        IF (maxdallowed_d .GT. 0) mxdiam(i) = mxdiam(i)*rndu(1)**(-(1./&
&           2.))
! Make MXDIAM stochastic
        DO l=k,icmin,-1
!*
          CALL PUSHREAL8(bet(l))
          bet(l) = dqq(l)*pki(l)
!*
          CALL PUSHREAL8(gam(l))
          gam(l) = pki(l)/(1.0+lbcp*dqq(l))
          IF (l .LT. k) THEN
            CALL PUSHREAL8(ght(l+1))
            ght(l+1) = gam(l)*dpb(l) + gam(l+1)*dpt(l+1)
            CALL PUSHREAL8(gm1(l+1))
            gm1(l+1) = 0.5*lbcp*(dqq(l)/(alhl*(1.0+lbcp*dqq(l)))+dqq(l+1&
&             )/(alhl*(1.0+lbcp*dqq(l+1))))
            CALL PUSHCONTROL1B(1)
          ELSE
            CALL PUSHCONTROL1B(0)
          END IF
        END DO
        rns = 0.
        cll = 0.
        rmf = 0.
        rmfd = 0.
        rmfc = 0.
        CALL PUSHREAL8ARRAY(rmfp, k0)
        rmfp = 0.
        cll0 = 0.
        dll0 = 0.
        dllx = 0.
        clli = 0.
        cllb = 0.
        poi_sv = poi
        qoi_sv = qoi
        uoi_sv = uoi
        voi_sv = voi
        IF (do_tracers) xoi_sv = xoi
        CALL PUSHREAL8ARRAY(cvw, k0)
        cvw = 0.0
        updfrc = 0.0
        updfrp = 0.0
        dissk0 = 0.0
        CALL PUSHCONTROL1B(0)
      ELSE
        CALL PUSHCONTROL1B(1)
      END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!    IF (PRESENT(FINAL)) THEN
      IF (final .EQ. 1) THEN
        CALL PUSHREAL8ARRAY(tho(i, icmin:k-1), k - icmin)
        tho(i, icmin:k-1) = poi(icmin:k-1)
        CALL PUSHREAL8ARRAY(qho(i, icmin:k-1), k - icmin)
        qho(i, icmin:k-1) = qoi(icmin:k-1)
        uho(i, icmin:k-1) = uoi(icmin:k-1)
        vho(i, icmin:k-1) = voi(icmin:k-1)
        cnv_updfrc(i, icmin:k-1) = updfrc(icmin:k-1)
!======================AER_CLOUD=============
!               CNV_NDROP   (I,ICMIN:K-1)  =    CNVNDROP(ICMIN:K-1) !DONIF
!               CNV_NICE   (I,ICMIN:K-1)   =     CNVNICE(ICMIN:K-1) !DONIF
!               CNV_FICE   (I,ICMIN:K-1)   =     CNVFICE(ICMIN:K-1) !DONIF
!! De-strap tendencies from RAS
!! specify weighting "SHAPE"
        CALL PUSHREAL8ARRAY(wght, k0)
        wght = wgt1(i, :)
!! Scale properly by layer masses
        wght0 = 0.
        DO l=k,k0
          wght0 = wght0 + wght(l)*(ple(i, l+1)-ple(i, l))
        END DO
        wght0 = (prs(k+1)-prs(k))/wght0
        wght = wght0*wght
        DO l=k,k0
          CALL PUSHREAL8(tho(i, l))
          tho(i, l) = tho(i, l) + wght(l)*(poi(k)-poi_sv(k))
          CALL PUSHREAL8(qho(i, l))
          qho(i, l) = qho(i, l) + wght(l)*(qoi(k)-qoi_sv(k))
          uho(i, l) = uho(i, l) + wght(l)*(uoi(k)-uoi_sv(k))
          vho(i, l) = vho(i, l) + wght(l)*(voi(k)-voi_sv(k))
        END DO
        IF (do_tracers) THEN
          xho(i, icmin:k-1, :) = xoi(icmin:k-1, :)
          DO itr=1,itrcr
            DO l=k,k0
              xho(i, l, itr) = xho(i, l, itr) + wght(l)*(xoi(k, itr)-&
&               xoi_sv(k, itr))
            END DO
          END DO
        END IF
!            FLX (I,ICMIN:K) = RMF (ICMIN:K) * DDT/DAYLEN  !  (KG/m^2/s @ CLOUD BASE)
!  (KG/m^2/s @ CLOUD TOP)
        flxd(i, icmin:k) = rmfd(icmin:k)*ddt/daylen
!            FLXC(I,ICMIN:K) = RMFC(ICMIN:K) * DDT/DAYLEN  !  (KG/m^2/s @ CLOUD TOP)
!  (KG/m^2/s )
        clw(i, icmin:k) = cll(icmin:k)*ddt/daylen
!            FLX (I,1:ICMIN-1) = 0.
        flxd(i, 1:icmin-1) = 0.
!            FLXC(I,1:ICMIN-1) = 0.
        clw(i, 1:icmin-1) = 0.
        IF (k .LT. k0) THEN
!               FLX (I,K:K0) = 0.
          flxd(i, k:k0) = 0.
!               FLXC(I,K:K0) = 0.
          clw(i, k:k0) = 0.
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
!            IRC (I,ICMIN:K-1) = RC(ICMIN:K-1)
        IF (ALLOCATED(icl_v)) THEN
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(0)
        END IF
      ELSE
        CALL PUSHCONTROL1B(1)
      END IF
      IF (final .EQ. 2) THEN
!            FLX (I,:) = 0.
        flxd(i, :) = 0.
!            FLXC(I,:) = 0.
        clw(i, :) = 0.
!            IRC (I,ICMIN:K-1) = RC(ICMIN:K-1)
        CALL PUSHREAL8ARRAY(wght, k0)
        CALL PUSHCONTROL1B(0)
      ELSE
        CALL PUSHREAL8ARRAY(wght, k0)
        CALL PUSHCONTROL1B(1)
      END IF
    END SUBROUTINE STRAP_FWD
    SUBROUTINE HTEST_FWD()
      IMPLICIT NONE
      REAL*8, DIMENSION(k0) :: hol1
      INTEGER :: lminhol
      REAL*8 :: minhol
      INTRINSIC MIN
      INTRINSIC MAX
! HOL initialized here in order not to confuse Valgrind debugger
      CALL PUSHREAL8ARRAY(hol, k0)
      hol = 0.
      CALL PUSHREAL8(zet(k+1))
      zet(k+1) = 0
      CALL PUSHREAL8(sht(k+1))
      sht(k+1) = cp*poi(k)*prj(k+1)
      DO l=k,icmin,-1
        IF (qst(l)*rhmax .GT. qoi(l)) THEN
          CALL PUSHREAL8(qol(l))
          qol(l) = qoi(l)
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHREAL8(qol(l))
          qol(l) = qst(l)*rhmax
          CALL PUSHCONTROL1B(1)
        END IF
        IF (0.000 .LT. qol(l)) THEN
          CALL PUSHCONTROL1B(0)
          qol(l) = qol(l)
        ELSE
          qol(l) = 0.000
          CALL PUSHCONTROL1B(1)
        END IF
        CALL PUSHREAL8(ssl(l))
        ssl(l) = cp*prj(l+1)*poi(l) + grav*zet(l+1)
        hol(l) = ssl(l) + qol(l)*alhl
        CALL PUSHREAL8(hst(l))
        hst(l) = ssl(l) + qst(l)*alhl
        CALL PUSHREAL8(tem)
        tem = poi(l)*(prj(l+1)-prj(l))*cpbg
        CALL PUSHREAL8(zet(l))
        zet(l) = zet(l+1) + tem
        CALL PUSHREAL8(zol(l))
        zol(l) = zet(l+1) + (prj(l+1)-prh(l))*poi(l)*cpbg
      END DO
    END SUBROUTINE HTEST_FWD
    SUBROUTINE CLOUDE_FWD(ic)
      IMPLICIT NONE
!=======================================
      INTEGER, INTENT(IN) :: ic
      REAL*8 :: deep_fact, cu_diam, wscale
!, dQx
      REAL*8 :: cli, te_a, c00_x, cli_crit_x, pete, toki, gmhx, hstx
      REAL*8 :: dt_lyr, rate, cvw_x, closs, f2, f3, f4, f5
      INTEGER :: k700
      INTRINSIC MIN
      INTRINSIC MAX
      INTRINSIC SQRT
      INTRINSIC EXP
      REAL*8 :: min2
      REAL*8 :: min1
      REAL*8 :: x5
      REAL*8 :: x4
      REAL*8 :: x3
      REAL*8 :: x2
      REAL*8 :: x1
      LOGICAL :: mask(k-ic+1)
      REAL*8 :: max2
      REAL*8 :: max1
      REAL*8 :: y1
!         !=============================AER_CLOUD local variables ====================
!         REAL(8) :: WBASE, NDROP, NICE, FP_D, FF_A, FP_I, FICE, &
!               NDROP_AMB, NSOOT_AMB, NSOOT, NIN, INSOOT, dCVW2, QICE, &
!               dQICE, dQIG, FPICE, dNICE, dNDROP, DSOOT_AMB, DSOOT, QLIQ, dQLIQ, FPRECIP, AUX, QT, &
!               MAXNICE, MAXNDROP, MINNICE, MINNDROP, NDROP_ACT, RIMM, FNDRIM, TminusTa, Tparcel , &
!	      alph_e, beta_e, RH_AMB, ECRIT
!
!         REAL(8), DIMENSION (NDUSTMAX) :: NDUST, NDUST_AMB, INDUST, DDUST_AMB, DDUST 
!         INTEGER :: INX
!
!
!         T_ICE_ALL = 238.0 !A little higher so there is ice at the freezing level
!         WBASE=1.0 
!         FICE=0.0
!         NICE=0.0
!         NDROP=0.0
!         NDROP_AMB=0.0
!         NDROP_ACT=0.0
!         NIN = 0.0
!         NDUST=0.0
!         NSOOT=0.0
!         NDUST_AMB=0.0
!         NSOOT_AMB = 0.0
!         dCVW2=0.0
!         QICE=0.0
!         QLIQ=0.0
!         FPICE = 0.0
!         INDUST=0.0
!         INSOOT= 0.0
!         QT= 0.0
!         FPRECIP=0.0   
!         FNDRIM = 0.0
!         RIMM = 0.0   
!         TminusTa = 0.0
!
!         f_seasalt = 0.0
!         aseasalt = 0.0
!
!         call init_Aer(AER_BASE)
!
!         !AER_CLOUD=============================
      IF (1. .GT. (qoi(k)/qst(k)-rhmn)/(rhmx-rhmn)) THEN
        CALL PUSHREAL8(trg)
        trg = (qoi(k)/qst(k)-rhmn)/(rhmx-rhmn)
        CALL PUSHCONTROL1B(0)
      ELSE
        CALL PUSHREAL8(trg)
        trg = 1.
        CALL PUSHCONTROL1B(1)
      END IF
      IF (0.0 .LT. (autorampb-sige(ic))/0.2) THEN
        y1 = (autorampb-sige(ic))/0.2
      ELSE
        y1 = 0.0
      END IF
      IF (1.0 .GT. y1) THEN
        f4 = y1
      ELSE
        f4 = 1.0
      END IF
      IF (trg .LE. 1.0e-5) THEN
        CALL PUSHCONTROL3B(0)
      ELSE
!  RECOMPUTE SOUNDING UP TO DETRAINMENT LEVEL
        poi_c = poi
        qoi_c = qoi
        poi_c(k) = poi_c(k) + tpert(i)
        qoi_c(k) = qoi_c(k) + qpert(i)
        CALL PUSHREAL8(zet(k+1))
        zet(k+1) = 0.
        CALL PUSHREAL8(sht(k+1))
        sht(k+1) = cp*poi_c(k)*prj(k+1)
        DO l=k,ic,-1
          IF (qst(l)*rhmax .GT. qoi_c(l)) THEN
            CALL PUSHREAL8(qol(l))
            qol(l) = qoi_c(l)
            CALL PUSHCONTROL1B(0)
          ELSE
            CALL PUSHREAL8(qol(l))
            qol(l) = qst(l)*rhmax
            CALL PUSHCONTROL1B(1)
          END IF
          IF (0.000 .LT. qol(l)) THEN
            CALL PUSHCONTROL1B(0)
            qol(l) = qol(l)
          ELSE
            qol(l) = 0.000
            CALL PUSHCONTROL1B(1)
          END IF
          CALL PUSHREAL8(ssl(l))
          ssl(l) = cp*prj(l+1)*poi_c(l) + grav*zet(l+1)
          CALL PUSHREAL8(hol(l))
          hol(l) = ssl(l) + qol(l)*alhl
          CALL PUSHREAL8(hst(l))
          hst(l) = ssl(l) + qst(l)*alhl
          CALL PUSHREAL8(tem)
          tem = poi_c(l)*(prj(l+1)-prj(l))*cpbg
          CALL PUSHREAL8(zet(l))
          zet(l) = zet(l+1) + tem
          CALL PUSHREAL8(zol(l))
          zol(l) = zet(l+1) + (prj(l+1)-prh(l))*poi_c(l)*cpbg
        END DO
        DO l=ic+1,k
          CALL PUSHREAL8(tem)
          tem = (prj(l)-prh(l-1))/(prh(l)-prh(l-1))
          CALL PUSHREAL8(sht(l))
          sht(l) = ssl(l-1) + tem*(ssl(l)-ssl(l-1))
          CALL PUSHREAL8(qht(l))
          qht(l) = .5*(qol(l)+qol(l-1))
        END DO
! SMOOTH HSTAR W/ 1-2-1 Filter
        IF (smooth_hst) THEN
! save for later
          hstx = hst(ic)
          DO l=k-1,ic+1,-1
            CALL PUSHREAL8(hst(l))
            hst(l) = 0.25*(hst(l+1)+hst(l-1)) + 0.5*hst(l)
          END DO
          DO l=ic,ic
            CALL PUSHREAL8(hst(l))
            hst(l) = 0.5*hst(l+1) + 0.5*hst(l)
          END DO
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
!  CALCULATE LAMBDA, ETA, AND WORKFUNCTION
        CALL PUSHREAL8(lambda_min)
        lambda_min = .2/mxdiam(i)
        lambda_max = .2/200.
!     LAMBDA_MIN = .2/(LAMBDA_FAC*DPTH_BL)
!     LAMBDA_MAX = .2/( MAX( LAMBMX_FAC*DPTH_BL , DIAMMN_MIN ) )
        IF (hol(k) .LE. hst(ic)) THEN
          CALL PUSHCONTROL3B(1)
        ELSE
!  LAMBDA CALCULATION: MS-A18
          CALL PUSHREAL8(tem)
          tem = (hst(ic)-hol(ic))*(zol(ic)-zet(ic+1))
          DO l=ic+1,k-1
            tem = tem + (hst(ic)-hol(l))*(zet(l)-zet(l+1))
          END DO
          IF (tem .LE. 0.0) THEN
            CALL PUSHCONTROL3B(2)
          ELSE
            CALL PUSHREAL8(alm)
            alm = (hol(k)-hst(ic))/tem
            IF (alm .GT. lambda_max) THEN
              CALL PUSHCONTROL3B(3)
            ELSE
              toki = 1.0
              IF (alm .LT. lambda_min) THEN
!we can probably replace this by a actual distribution based on grid cell size
                toki = (alm/lambda_min)**2
!RC(IC) = 6
!RETURN
                CALL PUSHCONTROL1B(1)
              ELSE
                CALL PUSHCONTROL1B(0)
              END IF
!LAMBDSV(IC) = ALM
!  ETA CALCULATION: MS-A2
              DO l=ic+1,k
                CALL PUSHREAL8(eta(l))
                eta(l) = 1.0 + alm*(zet(l)-zet(k))
              END DO
              CALL PUSHREAL8(eta(ic))
              eta(ic) = 1.0 + alm*(zol(ic)-zet(k))
!  WORKFUNCTION CALCULATION:  MS-A22
              wfn = 0.0
              CALL PUSHREAL8(hcc(k))
              hcc(k) = hol(k)
              DO l=k-1,ic+1,-1
                CALL PUSHREAL8(hcc(l))
                hcc(l) = hcc(l+1) + (eta(l)-eta(l+1))*hol(l)
                CALL PUSHREAL8(tem)
                tem = hcc(l+1)*dpb(l) + hcc(l)*dpt(l)
                CALL PUSHREAL8(eht(l))
                eht(l) = eta(l+1)*dpb(l) + eta(l)*dpt(l)
                wfn = wfn + (tem-eht(l)*hst(l))*gam(l)
              END DO
              CALL PUSHREAL8(hcc(ic))
              hcc(ic) = hst(ic)*eta(ic)
              wfn = wfn + (hcc(ic+1)-hst(ic)*eta(ic+1))*gam(ic)*dpb(ic)
!  VERTICAL VELOCITY/KE CALCULATION (ADDED 12/2001 JTB)
              bk3(k) = 0.0
              bk2(k) = 0.0
              bke(k) = 0.0
              CALL PUSHREAL8(hcld(k))
              hcld(k) = hol(k)
              DO l=k-1,ic,-1
                CALL PUSHREAL8(hcld(l))
                hcld(l) = (eta(l+1)*hcld(l+1)+(eta(l)-eta(l+1))*hol(l))/&
&                 eta(l)
                CALL PUSHREAL8(tem)
                tem = (hcld(l)-hst(l))*(zet(l)-zet(l+1))/(1.0+lbcp*dqq(l&
&                 ))
                IF (cldmicro .LE. 0.0) THEN
                  bke(l) = bke(l+1) + grav*tem/(cp*prj(l+1)*poi(l))
                  IF (tem .LT. 0.0) THEN
                    CALL PUSHREAL8(max1)
                    max1 = 0.0
                    CALL PUSHCONTROL1B(0)
                  ELSE
                    CALL PUSHREAL8(max1)
                    max1 = tem
                    CALL PUSHCONTROL1B(1)
                  END IF
                  bk2(l) = bk2(l+1) + grav*max1/(cp*prj(l+1)*poi(l))
                  IF (tem .GT. 0.0) THEN
                    min1 = 0.0
                  ELSE
                    min1 = tem
                  END IF
                  bk3(l) = bk3(l+1) + grav*min1/(cp*prj(l+1)*poi(l))
                  IF (bk2(l) .LT. 0.0) THEN
                    CALL PUSHREAL8(max2)
                    max2 = 0.0
                    CALL PUSHCONTROL1B(0)
                  ELSE
                    CALL PUSHREAL8(max2)
                    max2 = bk2(l)
                    CALL PUSHCONTROL1B(1)
                  END IF
                  CALL PUSHREAL8(cvw(l))
                  cvw(l) = SQRT(2.0*max2)
                  CALL PUSHCONTROL1B(1)
                ELSE
                  CALL PUSHCONTROL1B(0)
                END IF
              END DO
! 1.0 / ( 5.0*ALM )
!   ALPHA CALCULATION 
              CALL PUSHREAL8(rasal2i)
              rasal2i = rasal2_2d(i)
              IF (zet(ic) .LT. 2000.) THEN
                CALL PUSHREAL8(rasal(ic))
                rasal(ic) = rasal1
                CALL PUSHCONTROL1B(0)
              ELSE
                CALL PUSHCONTROL1B(1)
              END IF
              IF (zet(ic) .GE. 2000.) THEN
                IF (1.0 .GT. (zet(ic)-2000.)/8000.) THEN
                  min2 = (zet(ic)-2000.)/8000.
                  CALL PUSHCONTROL1B(0)
                ELSE
                  CALL PUSHCONTROL1B(1)
                  min2 = 1.0
                END IF
!WMP     RASAL(IC) = RASAL1 + (RASAL2i-RASAL1)*(ZET(IC) - 2000.)/8000.
                CALL PUSHREAL8(rasal(ic))
                rasal(ic) = rasal1 + (rasal2i-rasal1)*min2**rasal_exp
                CALL PUSHCONTROL1B(1)
              ELSE
                CALL PUSHCONTROL1B(0)
              END IF
!WMP  RASAL(IC) = MIN( RASAL(IC) , 1.0e5 )
              CALL PUSHREAL8(rasal(ic))
              rasal(ic) = dt/rasal(ic)
              mask(1:k-ic+1) = cvw(ic:k) .LT. 1.00
              CALL PUSHREAL8ARRAY(cvw(ic:k), k - ic + 1)
              WHERE (mask(1:k-ic+1)) cvw(ic:k) = 1.00
              CALL PUSHREAL8ARRAY(cvw(ic:k), k - ic + 1)
              WHERE (.NOT.mask(1:k-ic+1)) cvw(ic:k) = cvw(ic:k)
!  NOTE THIS "CENTRALIZES" A KLUGE PRESENT IN OTHER LOCATIONS.
!  CLEAN UP SOME TIME.      -JTB 12/04/03
!  TEST FOR CRITICAL WORK FUNCTION
              CALL ACRITN(pol(ic), prs(k), acr)
              IF (wfn .LE. acr) THEN
                CALL PUSHREAL8(max1)
                CALL PUSHREAL8(max2)
                CALL PUSHBOOLEANARRAY(mask, k - ic + 1)
                CALL PUSHREAL8(min2)
                CALL PUSHCONTROL3B(4)
              ELSE
!  CLOUD TOP WATER AND MOMENTUM (TIMES ETA(IC)) MS-A16
! Tracer scavenging
! RAS loops over a series of plumes all having common cloud base level K 
! and different detrainment levels IC.  The plumes operate sequentially
! on the grid box mean quantities (wind, moisture, tracer) and so each
! subsequent plume is seeing the effects of previous plumes.  We parameterize
! scavenging following Liu et al. [JGR, 2001], their equation 1:
!  AEROSOL FRACTION SCAVENGED = 1 - exp(-FSCAV*DZ)
! where FSCAV is a specified scavenging efficiency [km-1] and DZ is the
! distance [km] the tracer traverses in the plume from it's entrainment
! level to its detrainment level.  We write the aerosol fraction surviving as:
!  FNOSCAV = exp(- FSCAV_(ITR) * DZ)
! The total scavenging is proportional to the convective mass flux, which
! is not explicitly solved for at this point.
                IF (do_tracers) THEN
                  DO itr=1,itrcr
!           Scavenging of the below cloud tracer
                    delzkm = (zet(ic)-zet(k))/1000.
                    x4 = EXP(-(fscav_(itr)*delzkm))
                    IF (x4 .GT. 1.) THEN
                      x1 = 1.
                    ELSE
                      x1 = x4
                    END IF
                    IF (x1 .LT. 0.) THEN
                      fnoscav = 0.
                    ELSE
                      fnoscav = x1
                    END IF
                    xht(itr) = xoi(k, itr)*fnoscav
                  END DO
                  CALL PUSHCONTROL1B(0)
                ELSE
                  CALL PUSHCONTROL1B(1)
                END IF
                CALL PUSHREAL8(wlq)
                wlq = qol(k)
                CALL PUSHREAL8(uht)
                uht = uoi(k)
                CALL PUSHREAL8(vht)
                vht = voi(k)
                CALL PUSHREAL8(rnn(k))
                rnn(k) = 0.
                cll0(k) = 0.
!print *, '========================================='
                DO l=k-1,ic,-1
                  CALL PUSHREAL8(tem)
                  tem = eta(l) - eta(l+1)
                  wlq = wlq + tem*qol(l)
                  uht = uht + tem*uoi(l)
                  vht = vht + tem*voi(l)
                  IF (do_tracers) THEN
                    DO itr=1,itrcr
!         Scavenging of the entrained tracer.  Updates transported tracer mass.
                      delzkm = (zet(ic)-zet(l+1))/1000.
                      x5 = EXP(-(fscav_(itr)*delzkm))
                      IF (x5 .GT. 1.) THEN
                        x2 = 1.
                      ELSE
                        x2 = x5
                      END IF
                      IF (x2 .LT. 0.) THEN
                        fnoscav = 0.
                      ELSE
                        fnoscav = x2
                      END IF
                      xht(itr) = xht(itr) + tem*xoi(l, itr)*fnoscav
                    END DO
                    CALL PUSHCONTROL1B(0)
                  ELSE
                    CALL PUSHCONTROL1B(1)
                  END IF
!!!! How much condensate (CLI) is present here? 
                  IF (l .GT. ic) THEN
                    CALL PUSHREAL8(tx2)
                    tx2 = 0.5*(qst(l)+qst(l-1))*eta(l)
                    CALL PUSHREAL8(tx3)
                    tx3 = 0.5*(hst(l)+hst(l-1))*eta(l)
                    qcc = tx2 + gm1(l)*(hcc(l)-tx3)
                    cll0(l) = wlq - qcc
                    CALL PUSHCONTROL1B(1)
                  ELSE
                    cll0(l) = wlq - qst(ic)*eta(ic)
                    CALL PUSHCONTROL1B(0)
                  END IF
                  IF (cll0(l) .LT. 0.00) THEN
                    cll0(l) = 0.00
                    CALL PUSHCONTROL1B(0)
                  ELSE
                    CALL PUSHCONTROL1B(1)
                    cll0(l) = cll0(l)
                  END IF
! condensate (kg/kg)
                  cli = cll0(l)/eta(l)
! Temperature (K)
                  te_a = poi(l)*prh(l)
!=====================================================================
                  IF (cldmicro .GT. 0.0) THEN
                    CALL PUSHCONTROL1B(1)
                  ELSE
!AER_CLOUD MIcrophysics considering activation and nucleation 
!               !recompute vertical velocity
!
!               Tparcel = TE_A
!               CVW(K) = 0.8  ! Assume a below cloud base  W of 0.8 m s-1            
!               BK2(K)   = 0.0
!
!     
!               TEM     = (HCLD(L)-HST(L) )/ (1.0+LBCP*DQQ(L))  
!               TminusTa = max(min(TEM/CP, 5.0), 0.0) !limit DT to 5 K. According to Wei, JAS, 1998   
!	     TEM =0.33*TminusTa*CO_AUTO(I)/TE_A !Bouyancy term, effciency =0.5 mwr Roode et al    	     
!
!               BK2(L)  = BK2(L+1) + GRAV * TEM*(ZET(L)-ZET(L+1)) 
!               BK2(L) = BK2(L) - (ZET(L)-ZET(L+1))*(BK2(L+1)*ALM + CLI*GRAV)  !Account for drag from entrainment of stagnat air 
!and condesate loading
!               CVW(L) = max(SQRT(  2.0* MAX( BK2(L) , 0.0 )  ), 1.0) 
!
!
!	    CVW_X = MIN(CVW(L), 50.0)
!               DT_LYR  =  max(( ZET(L)-ZET(L+1) )/CVW_X, 1.0) !Sanity check 
!               TEM   = ETA(L) - ETA(L+1)
!
!               Tparcel  =  TE_A + TminusTa
!
!
!
!!!!!!!!!!account for entrainment effects on activation !!!!!!!!!!!
!               ! Barahona and Nenes, JGR, 2007
!               alph_e = 2.8915e-8*Tparcel*Tparcel -2.1328e-5*Tparcel+4.2523e-3
!               beta_e = MAPL_ALHL*TminusTa/MAPL_RVAP/Tparcel/Tparcel
!               RH_AMB=QOI(L)/QST(L)
!               ECRIT  = max(RH_AMB -beta_e, 1.0e-6) 
!               ECRIT =  alph_e/ECRIT
!               ! print *, L, Tparcel, RH_AMB, ECRIT, ALM
!	           ECRIT =  ALM/ECRIT
!       ! ECRIT = 0.0 ! do not use this for now
!               !Print *, ECRIT
!
!
!               if (L .eq. K-1) then
!
!                  FICE=0.0
!                  NICE=0.0
!                  NDROP=0.0
!                  NIN =0.0
!                  NDUST_AMB =0.0
!                  NSOOT_AMB = 0.0
!                  NSOOT=0.0
!                  NDUST= 0.0
!
!        
!                  AER_BASE =  AERO(L)
!                  RATE=0.0
!                  FPRECIP=0.0
!                  !initial conditions     
!                  call ARGact(Tparcel, CVW_X, NDROP_ACT, NDROP_AMB, NDUST_AMB, NSOOT_AMB, L,  .true., DDUST_AMB, DSOOT_AMB, ECRI
!T) !cloud droplet number and INsource at cloud base 
!                  NDUST=NDUST_AMB
!                  NSOOT=NSOOT_AMB
!                  DDUST=DDUST_AMB
!                  DSOOT=DSOOT_AMB                                     
!
!               else 
!                  call ARGact(Tparcel, CVW_X, NDROP_ACT, NDROP_AMB, NDUST_AMB, NSOOT_AMB, L, .false., DDUST_AMB, DSOOT_AMB, ECRI
!T) !cloud droplet number above cloud base  
!
!               end if
!
!               QT = CLI
!               RATE = 0.0
!               FPRECIP = 0.0
!
!               if (QT .gt. 0.0) then
!
!                  ! if FICE is already >= 1.0 then the cloud is glaciated and there is no need to do anymore partitioning
!
!                  if (FICE .ge. 1.0) then
!
!
!                     CALL  Qremoval(RATE, FICE, FP_D, FP_I, Tparcel,  & 
!                           POL(L), QT,  NICE, NDROP, CVW_X, FPICE, &
!                           DT_LYR, RIMM, CO_AUTO(I)) 
!
!
!
!                     dNICE = -NICE*FP_I 
!                     NICE  =  (NICE +dNICE*DT_LYR)*ETA(L+1)/ETA(L) !ice
!
!                     MINNICE = max(QICE*(1.0-RATE*DT_LYR), 0.0)/4.0e-8!assuming maximum vol radius 250 microns
!                     MAXNICE = max(QICE*(1.0-RATE*DT_LYR), 0.0)/2.51e-12 !assuming minimum vol radius 10 microns
!
!                     NICE=MIN(max(NICE, MINNICE), MAXNICE)
!
!                     FICE = 1.0
!
!                  else 
!
!                     ! Cloud is not completely glaciated do the whole thing
!                     ! ALL this subroutines return tendencies
!
!
!                     CALL  INfreezing(QLIQ, NDROP, NIN, NDUST, NSOOT, INDUST, INSOOT, Tparcel, POL(L), CVW_X, DDUST, DSOOT)  !ca
!lculate the freezing fraction of the aerosol at this level
!
!                     NIN = min(NIN, NDROP/DT_LYR)
!
!                     call Qgrowth(Tparcel, POL(L), QICE, NICE, QT, NIN, dQIG, RIMM, FNDRIM)
!
!                     CALL  Qremoval(RATE, FICE, FP_D, FP_I, Tparcel,  & 
!                           POL(L), QT,  NICE, NDROP, CVW_X, FPICE, &
!                           DT_LYR,  RIMM, CO_AUTO(I)) 
!
!
!
!                     !ice number tendency: -precip + freezin
!                     dNICE = -NICE*FP_I  + NIN     
!                     NICE  =  (NICE +dNICE*DT_LYR)*ETA(L+1)/ETA(L) !ice
!                     NICE =max(NICE, 0.0)
!
!
!                     !ice mass tendency: growth - precip
!                     dQICE = -QICE*FPICE + dQIG
!                     QICE  =  min((QICE + dQICE*DT_LYR)*ETA(L+1)/ETA(L), QT) !ice
!                     QICE=max(min(QT, QICE), 0.0)
!
!
!                     ! Liquid Tendency: source/evap -  precip 
!                     !dQLIQ = max((CLI-QICE), -QLIQ)/DT_LYR -QLIQ*max(RATE-FPICE, 0.0) 
!                     ! dQLIQ = CLI*(1.0-RATE*DT_LYR)/DT_LYR -dQICE - QLIQ*max(RATE-FPICE, 0.0)           
!                     !QLIQ  =  max((QLIQ + dQLIQ*DT_LYR)*ETA(L+1)/ETA(L), 0.0) !liquid. This is actually diagnostic
!                     QLIQ=max((QT-QICE), 0.0)
!
!
!                     !droplet number tendency: -precip - freezin + activation + activated entrained aerosol 
!
!
!                     dNDROP =-NDROP*FP_D - NIN -  FNDRIM*NDROP/DT_LYR + max(NDROP_ACT-NDROP, 0.0)/DT_LYR          
!
!                     !dNDROP =-NDROP*FP_D - NIN -  FNDRIM*NDROP/DT_LYR + NDROP_ACT/DT_LYR         
!
!                     NDROP =  (NDROP + dNDROP*DT_LYR)*ETA(L+1)/ETA(L) + &
!                           (ZET(L) - ZET(L+1))*ALM*MAX((NDROP_AMB-NDROP), 0.0)
!
!                     !Aerosol tendency: Entrainment - freezing 
!
!                     NDUST = (NDUST - INDUST*DT_LYR)*ETA(L+1)/ETA(L) + &
!                           (ZET(L) - ZET(L+1))*ALM*MAX(NDUST_AMB-NDUST, 0.0) 
!
!                     NSOOT =  (NSOOT - INSOOT*DT_LYR)*ETA(L+1)/ETA(L)  + &     
!                           (ZET(L) - ZET(L+1))*ALM*MAX(NSOOT_AMB-NSOOT, 0.0)  
!
!
!                           
!                     !Update FICE and perform Sanity checks
!
!
!                     MINNDROP = (1.0-FICE)*QLIQ*max(1.0-RATE*DT_LYR, 0.0)/2.e-10    !assuming maximum vol radius 36 microns
!                     MAXNDROP = (1.0-FICE)*QLIQ*max(1.0-RATE*DT_LYR, 0.0)/3.35e-14 !assuming minimum vol radius 2 microns
!                     MINNICE = QICE/4.0e-8!assuming maximum vol radius 250 microns
!                     MAXNICE = QICE/2.51e-12 !assuming minimum vol radius 10 microns
!
!                     IF ((NICE .gt. MAXNICE) .or. (NICE .lt. MINNICE))   then    
!                        !print *, 'nilim', NICE*1e-6, MINNICE*1e-6, MAXNICE*1e-6
!                     END IF
!
!                     IF ((NDROP .gt. MAXNDROP) .or. (NDROP .lt. MINNDROP))      then 
!                        !print *, 'ndroplim', NDROP*1e-6, MINNDROP*1e-6, MAXNDROP*1e-6
!                     end if
!
!
!                     NSOOT=MAX(NSOOT, 0.0)
!                     NDUST=MAX(NDUST, 0.0)              
!
!                     NDROP=MIN(max(NDROP, MINNDROP), MAXNDROP)
!                     NICE=MIN(max(NICE, MINNICE), MAXNICE)
!
!                     FICE=max(min(QICE/QT, 1.0), 0.0)
!
!                     IF (FICE .ge. 1.0) THEN !Complete glaciation 
!                        NICE=NICE+NDROP 
!                        NDROP = 0.0
!                        QICE  = QT
!                        QLIQ= 0.0
!                     END IF
!
!                     IF (Tparcel .LT. T_ICE_ALL) THEN !instantaneous freezing
!                        NICE=NICE+NDROP 
!                        NDROP = 0.0
!                        FICE  = 1.0
!                        QICE  = QT
!                        QLIQ=0.0
!                     END IF
!
!                     IF (Tparcel .GT. T_ICE_MAX) THEN !instantaneous melting
!                        NDROP=NICE+NDROP 
!                        NICE = 0.0
!                        FICE  = 0.0
!                        QICE  = 0.0
!                        QLIQ=QT
!                     END IF
!
!                  END IF
!
!               else 
!
!                  FICE =0.0 
!                  QICE = 0.0
!                  QLIQ = 0.0
!                  NICE= 0.0 
!                  NDROP = 0.0
!                  RATE =0.0
!               end if
!
!               FPRECIP= RATE*DT_LYR
!
!               !RATE=RATE*F4
!               ! NDROP=NDROP*F4
!               !NICE=NICE*(1.0-F4)
!
!               !print *, TE_A, FICE, 'NICE', NICE*1e-6, 'NDROP', NDROP*1e-6, L 
!               !print *, 'FPI', FP_I*DT_LYR, 'FPD', FP_D*DT_LYR, 'FPICE', FPICE, 'FPRE', FPRECIP, QT, QLIQ
!
!Bacmeister 2006 microphysics
                    CALL SUNDQ3_ICE_FWD(te_a, sdqv2, sdqv3, sdqvt1, f2, &
&                                 f3)
! * F5  ! F4 reduces AUTO for shallow clouds, F5 modifies auto for deep clouds
                    CALL PUSHREAL8(c00_x)
                    c00_x = co_auto(i)*f2*f3*f4
                    cli_crit_x = cli_crit/(f2*f3)
                    CALL PUSHREAL8(rate)
                    rate = c00_x*(1.0-EXP(-(cli**2/cli_crit_x**2)))
                    CALL PUSHCONTROL1B(0)
                  END IF
                  IF (cvw(l) .LT. 1.00) THEN
                    CALL PUSHREAL8(cvw_x)
                    cvw_x = 1.00
                    CALL PUSHCONTROL1B(0)
                  ELSE
                    CALL PUSHREAL8(cvw_x)
                    cvw_x = cvw(l)
                    CALL PUSHCONTROL1B(1)
                  END IF
! really trust it at low values
! l.h.s. DT_LYR => time in layer (L,L+1)
                  dt_lyr = (zet(l)-zet(l+1))/cvw_x
                  closs = cll0(l)*rate*dt_lyr
                  IF (closs .GT. cll0(l)) THEN
                    closs = cll0(l)
                    CALL PUSHCONTROL1B(0)
                  ELSE
                    closs = closs
                    CALL PUSHCONTROL1B(1)
                  END IF
                  CALL PUSHREAL8(cll0(l))
                  cll0(l) = cll0(l) - closs
                  dll0(l) = closs
                  IF (closs .GT. 0.) THEN
                    wlq = wlq - closs
                    CALL PUSHREAL8(rnn(l))
                    rnn(l) = closs
                    CALL PUSHCONTROL1B(1)
                  ELSE
                    CALL PUSHREAL8(rnn(l))
                    rnn(l) = 0.
                    CALL PUSHCONTROL1B(0)
                  END IF
                END DO
!AER_CLOUD=======================================
!            CNVNDROP(IC)=NDROP
!            CNVNICE(IC)=NICE
!            CNVFICE(IC)=FICE
                wlq = wlq - qst(ic)*eta(ic)
!     CALCULATE GAMMAS AND KERNEL
! MS-A30 (W/O GRAV)
                CALL PUSHREAL8(gms(k))
                gms(k) = (sht(k)-ssl(k))*pri(k)
! MS-A31 (W/O GRAV)
                CALL PUSHREAL8(gmh(k))
                gmh(k) = gms(k) + (qht(k)-qol(k))*pri(k)*alhl
! MS-A37 (W/O GRAV)
                CALL PUSHREAL8(akm)
                akm = gmh(k)*gam(k-1)*dpb(k-1)
                CALL PUSHREAL8(tx2)
                tx2 = gmh(k)
                DO l=k-1,ic+1,-1
                  CALL PUSHREAL8(gms(l))
                  gms(l) = (eta(l)*(sht(l)-ssl(l))+eta(l+1)*(ssl(l)-sht(&
&                   l+1)))*pri(l)
                  CALL PUSHREAL8(gmh(l))
                  gmh(l) = gms(l) + (eta(l)*(qht(l)-qol(l))+eta(l+1)*(&
&                   qol(l)-qht(l+1)))*alhl*pri(l)
                  CALL PUSHREAL8(tx2)
                  tx2 = tx2 + (eta(l)-eta(l+1))*gmh(l)
                  akm = akm - gms(l)*eht(l)*pki(l) + tx2*ght(l)
                END DO
                CALL PUSHREAL8(gms(ic))
                gms(ic) = eta(ic+1)*(ssl(ic)-sht(ic+1))*pri(ic)
                akm = akm - gms(ic)*eta(ic+1)*dpb(ic)*pki(ic)
                CALL PUSHREAL8(gmh(ic))
                gmh(ic) = gms(ic) + (eta(ic+1)*(qol(ic)-qht(ic+1))*alhl+&
&                 eta(ic)*(hst(ic)-hol(ic)))*pri(ic)
                IF (smooth_hst) THEN
                  gmhx = gms(ic) + (eta(ic+1)*(qol(ic)-qht(ic+1))*alhl+&
&                   eta(ic)*(hstx-hol(ic)))*pri(ic)
                  CALL PUSHCONTROL1B(0)
                ELSE
                  CALL PUSHCONTROL1B(1)
                END IF
!    CLOUD BASE MASS FLUX
                IF (akm .GE. 0.0 .OR. wlq .LT. 0.0) THEN
                  CALL PUSHREAL8(max1)
                  CALL PUSHREAL8(max2)
                  CALL PUSHBOOLEANARRAY(mask, k - ic + 1)
                  CALL PUSHREAL8(f2)
                  CALL PUSHREAL8(f3)
                  CALL PUSHREAL8(f4)
                  CALL PUSHREAL8(c00_x)
                  CALL PUSHREAL8(cvw_x)
                  CALL PUSHREAL8(min2)
                  CALL PUSHREAL8(rate)
                  CALL PUSHREAL8(hstx)
                  CALL PUSHCONTROL3B(5)
                ELSE
! MS-A39 MASS-FLUX IN Pa/step
                  CALL PUSHREAL8(wfn)
                  wfn = -((wfn-acr)/akm)
                  x3 = rasal(ic)*trg*toki*wfn
                  IF (x3 .GT. (prs(k+1)-prs(k))*(100.*pblfrac)) THEN
                    CALL PUSHREAL8(wfn)
                    wfn = (prs(k+1)-prs(k))*(100.*pblfrac)
                    CALL PUSHCONTROL1B(0)
                  ELSE
                    CALL PUSHREAL8(wfn)
                    wfn = x3
                    CALL PUSHCONTROL1B(1)
                  END IF
!    CUMULATIVE PRECIP AND CLOUD-BASE MASS FLUX FOR OUTPUT
                  wfnog = wfn*gravi
                  CALL PUSHREAL8(tem)
                  tem = wfn*gravi
! (kg/m^2/step)
                  cll(ic) = cll(ic) + wlq*tem
! (kg/m^2/step)
                  rmf(ic) = rmf(ic) + tem
! (kg/m^2/step)
                  rmfd(ic) = rmfd(ic) + tem*eta(ic)
                  DO l=ic+1,k
! (kg/m^2/step)
                    CALL PUSHREAL8(rmfp(l))
                    rmfp(l) = tem*eta(l)
! (kg/m^2/step)
                    rmfc(l) = rmfc(l) + rmfp(l)
                    dllx(l) = dllx(l) + tem*dll0(l)
                    IF (cvw(l) .GT. 0.0) THEN
                      updfrp(l) = rmfp(l)*(ddt/daylen)*1000./(cvw(l)*prs&
&                       (l))
                      CALL PUSHCONTROL1B(0)
                    ELSE
                      updfrp(l) = 0.0
                      CALL PUSHCONTROL1B(1)
                    END IF
! current cloud; incloud condensate        
                    clli(l) = cll0(l)/eta(l)
!  cumulative grid mean convective condensate        
                    cllb(l) = cllb(l) + updfrp(l)*clli(l)
                    updfrc(l) = updfrc(l) + updfrp(l)
                  END DO
!    THETA AND Q CHANGE DUE TO CLOUD TYPE IC
                  DO l=ic,k
! (kg/m^2/step)
                    rns(l) = rns(l) + rnn(l)*tem
                    CALL PUSHREAL8(gmh(l))
                    gmh(l) = gmh(l)*wfn
                    CALL PUSHREAL8(gms(l))
                    gms(l) = gms(l)*wfn
                    CALL PUSHREAL8(qoi(l))
                    qoi(l) = qoi(l) + (gmh(l)-gms(l))*alhi
                    CALL PUSHREAL8(poi(l))
                    poi(l) = poi(l) + gms(l)*pki(l)*cpi
                    CALL PUSHREAL8(qst(l))
                    qst(l) = qst(l) + gms(l)*bet(l)*cpi
                  END DO
                  IF (smooth_hst) THEN
                    CALL PUSHREAL8(gmhx)
                    gmhx = gmhx*wfn
                    dqx = (gmhx-gmh(ic))*alhi
                    rns(ic) = rns(ic) + dqx/(pri(ic)*grav)
                    CALL PUSHCONTROL1B(0)
                  ELSE
                    CALL PUSHCONTROL1B(1)
                  END IF
                  IF (do_tracers) THEN
!*FRICFAC*0.5
                    CALL PUSHREAL8(wfn)
                    wfn = wfn*0.5*1.0
                    CALL PUSHREAL8(tem)
                    tem = wfn*pri(k)
                    CALL PUSHCONTROL1B(0)
                    DO itr=1,itrcr
                      xcu(k, itr) = xcu(k, itr) + tem*(xoi(k-1, itr)-xoi&
&                       (k, itr))
                    END DO
                    DO itr=1,itrcr
                      DO l=k-1,ic+1,-1
                        tem = wfn*pri(l)
                        xcu(l, itr) = xcu(l, itr) + tem*((xoi(l-1, itr)-&
&                         xoi(l, itr))*eta(l)+(xoi(l, itr)-xoi(l+1, itr)&
&                         )*eta(l+1))
                      END DO
                    END DO
                    tem = wfn*pri(ic)
                    DO itr=1,itrcr
                      xcu(ic, itr) = xcu(ic, itr) + (2.*(xht(itr)-xoi(ic&
&                       , itr)*(eta(ic)-eta(ic+1)))-(xoi(ic, itr)+xoi(ic&
&                       +1, itr))*eta(ic+1))*tem
                    END DO
                    DO itr=1,itrcr
                      DO l=ic,k
                        xoi(l, itr) = xoi(l, itr) + xcu(l, itr)
                      END DO
                    END DO
                  ELSE
!*FRICFAC*0.5
                    CALL PUSHREAL8(wfn)
                    wfn = wfn*0.5*1.0
                    CALL PUSHCONTROL1B(1)
                  END IF
!   CUMULUS FRICTION
                  IF (fricfac .LE. 0.0) THEN
                    CALL PUSHREAL8(max1)
                    CALL PUSHREAL8(max2)
                    CALL PUSHBOOLEANARRAY(mask, k - ic + 1)
                    CALL PUSHREAL8(f2)
                    CALL PUSHREAL8(f3)
                    CALL PUSHREAL8(f4)
                    CALL PUSHREAL8(c00_x)
                    CALL PUSHREAL8(cvw_x)
                    CALL PUSHREAL8(toki)
                    CALL PUSHREAL8(min2)
                    CALL PUSHREAL8(rate)
                    CALL PUSHREAL8(hstx)
                    CALL PUSHCONTROL3B(6)
                  ELSE
                    CALL PUSHREAL8(wfn)
                    wfn = wfn*fricfac*EXP(-(alm/friclambda))
                    CALL PUSHREAL8(tem)
                    tem = wfn*pri(k)
                    ucu(k) = ucu(k) + tem*(uoi(k-1)-uoi(k))
                    vcu(k) = vcu(k) + tem*(voi(k-1)-voi(k))
                    DO l=k-1,ic+1,-1
                      tem = wfn*pri(l)
                      ucu(l) = ucu(l) + tem*((uoi(l-1)-uoi(l))*eta(l)+(&
&                       uoi(l)-uoi(l+1))*eta(l+1))
                      vcu(l) = vcu(l) + tem*((voi(l-1)-voi(l))*eta(l)+(&
&                       voi(l)-voi(l+1))*eta(l+1))
                    END DO
                    tem = wfn*pri(ic)
                    ucu(ic) = ucu(ic) + (2.*(uht-uoi(ic)*(eta(ic)-eta(ic&
&                     +1)))-(uoi(ic)+uoi(ic+1))*eta(ic+1))*tem
                    vcu(ic) = vcu(ic) + (2.*(vht-voi(ic)*(eta(ic)-eta(ic&
&                     +1)))-(voi(ic)+voi(ic+1))*eta(ic+1))*tem
                    dissk0(ic) = eta(ic)*grav*wfnog*pri(ic)*0.5*((uht/&
&                     eta(ic)-uoi(ic))**2+(vht/eta(ic)-voi(ic))**2)
                    DO l=ic,k
                      CALL PUSHREAL8(uoi(l))
                      uoi(l) = uoi(l) + ucu(l)
                      CALL PUSHREAL8(voi(l))
                      voi(l) = voi(l) + vcu(l)
                    END DO
                    CALL PUSHREAL8(max1)
                    CALL PUSHREAL8(max2)
                    CALL PUSHBOOLEANARRAY(mask, k - ic + 1)
                    CALL PUSHREAL8(f2)
                    CALL PUSHREAL8(f3)
                    CALL PUSHREAL8(f4)
                    CALL PUSHREAL8(c00_x)
                    CALL PUSHREAL8(cvw_x)
                    CALL PUSHREAL8(toki)
                    CALL PUSHREAL8(min2)
                    CALL PUSHREAL8(rate)
                    CALL PUSHREAL8(hstx)
                    CALL PUSHCONTROL3B(7)
                  END IF
                END IF
              END IF
            END IF
          END IF
        END IF
      END IF
    END SUBROUTINE CLOUDE_FWD
    SUBROUTINE RNEVP_FWD()
      IMPLICIT NONE
      CALL PUSHREAL8(zet(k+1))
      zet(k+1) = 0
      DO l=k,icmin,-1
        CALL PUSHREAL8(tem)
        tem = poi(l)*(prj(l+1)-prj(l))*cpbg
        CALL PUSHREAL8(zet(l))
        zet(l) = zet(l+1) + tem
      END DO
      DO l=icmin,k
        CALL PUSHREAL8(tem)
        tem = pri(l)*grav
        cnv_prc3(i, l) = rns(l)*tem
      END DO
!! If hst is smoothed then adjusted precips may be negative
      IF (smooth_hst) THEN
        DO l=icmin,k
          IF (cnv_prc3(i, l) .LT. 0.) THEN
            CALL PUSHREAL8(qoi(l))
            qoi(l) = qoi(l) + cnv_prc3(i, l)
            CALL PUSHREAL8(poi(l))
            poi(l) = poi(l) - cnv_prc3(i, l)*(alhl/cp)/prj(l+1)
            cnv_prc3(i, l) = 0.
            CALL PUSHCONTROL1B(1)
          ELSE
            CALL PUSHCONTROL1B(0)
          END IF
        END DO
        CALL PUSHCONTROL1B(0)
      ELSE
        CALL PUSHCONTROL1B(1)
      END IF
    END SUBROUTINE RNEVP_FWD
  END SUBROUTINE RASE_FWD
!  Differentiation of rase in reverse (adjoint) mode, backward sweep (with options r8 split(GEOS_MoistGridComp.PRE_RASE GEOS_Mois
!tGridComp.PRE_PROGNO_CLOUD ras.RASE rase.CLOUDE rase.ACRITN rase.RNEVP rase.HTEST rase.FINDDTLS rase.STRAP ras.SUNDQ3_ICE cloudn
!ew.pdf_spread cloudnew.fix_up_clouds cloudnew.meltfrz cloudnew.hystpdf cloudnew.pdffrac cloudnew.pdfcondensate cloudnew.cnvsrc c
!loudnew.evap3 cloudnew.subl3 cloudnew.autocon3 cloudnew.PRECIP3 cloudnew.ICEFALL cloudnew.SETTLE_VEL cloudnew.MARSHPALMQ2 cloudn
!ew.MICRO_AA_BB_3 cloudnew.LDRADIUS3 cloudnew.ICE_FRACTION cloudnew.GET_ALHX3 cloudnew.ICEFRAC cloudnew.SUNDQ3_ICE3 qsat_util.DQS
!ATPERT)):
!   gradient     of useful results: clw cnv_prc3 tho qho vho cnv_updfrc
!                uho flxd
!   with respect to varying inputs: tho dqs qho vho tpert uho qss
!         pko, plo, phio, phie, qlo, qio,                  &
!FLX, FLXC,                  &
!         CNV_CVW,                                         &
!         CNV_QC,                                          &
!         ENTLAM,                                          &
!         CLAN,                                            &
!         HHO, HSO,PRECU,                                  &
!         AEROPROPS,                                   & !!!!!AER_CLOUD
!         CNV_FICE,                                        & 
!         CNV_NICE,                                        & 
!         CNV_NDROP,                                    &   !DONIF
  SUBROUTINE RASE_BWD(idim, irun, k0, icmin, dt, doras, cpo, alhlo, &
&   alhl1, tice, gravo, seedras, sige, kcbl, wgt0, wgt1, tpert, tpert_ad&
&   , qpert, tho, tho_ad, qho, qho_ad, uho, uho_ad, vho, vho_ad, qss, &
&   qss_ad, dqs, dqs_ad, cnv_fraction, rasal2_2d, co_auto, ple, pke, clw&
&   , clw_ad, flxd, flxd_ad, cnv_prc3, cnv_prc3_ad, cnv_updfrc, &
&   cnv_updfrc_ad, rasparams, itrcr, xho, triedlev_diag, fscav, disske)
    IMPLICIT NONE
!!!!!!!!!!!!!!!======================================     
!!!!!!!!!!   Subroutine ARGact: finds the activated droplet number from Abdul_Razzak and Ghan 2000.
!      ! Tailored for GOCART AEROSOL and works with AERO input 
!      !Written by Donifan Barahona
!      !!donifan.o.barahona@nasa.gov
!!!!!!!!!!!!!!!====================================
!
!      SUBROUTINE ARGact (TEMP, WX, NCPL_ACT, NCPL_AMB,  CDUST, CSOOT, LEV, ISBASE, DDUSTAMB, DSOOTAMB, ENT_PARAM)
!         !
!         Integer, intent(in)     ::  LEV    
!         LOGICAL,  intent(in)     ::   ISBASE
!
!         REAL(8), intent(inout)     ::   TEMP, WX, ENT_PARAM
!         REAL(8), intent(out)     ::   NCPL_ACT, NCPL_AMB, CSOOT, DSOOTAMB
!         REAL(8), DIMENSION(NDUSTMAX), INTENT(OUT) :: CDUST, DDUSTAMB
!         integer                 :: INDEX, NMODES, naux       
!
!         type(AerProps)  :: AER, auxaer      
!
!
!         REAL(8)     ::      kappa, alfa, beta, Akoh, G, T, smax, fi, gi, nui, &
!                       citai, ui, aux1, PACT,  Ntot, auxx, aux, auxconc, W, alph, aseasalt_aux, f_seasalt1
!         REAL(8), dimension (30) ::  SMI, TPI, SIGI 
!
!
!         SMI=0.0      
!         TPI = 0.0
!         SIGI =2.0
!         NCPL_ACT=0.0
!         NCPL_AMB=0.0
!         CDUST=0.0
!         CSOOT=0.0
!         DDUSTAMB =1.0e-9
!         DSOOTAMB= 1.0e-9
!         W=MIN(WX*(1.0-ENT_PARAM), 20.0)    
!         call init_Aer(AER)  
!         AER  =   AERO(LEV) 
!
!         PACT=0.0 !activation probability of entrained aerosol 
!   auxconc =0.0
!    aseasalt_aux  = 0.0
!!!!!!!!!!!activate aerosol transported from cloud base 
!             NMODES =  AER_BASE%nmods
!             TPI(1:nmodes) = AER_BASE%num(1:nmodes)
!             SIGI(1:nmodes) = AER_BASE%sig(1:nmodes)                          
!
!             
!             Ntot= 0.0
!              do index = 1, nmodes 
!	              if (AER_BASE%kap(index) .gt. 0.1) Ntot =  Ntot + TPI(index)  
!              end do
!         
!
!         if ((Ntot .lt. 1.0e4) .or. (TEMP .lt. 245.0) .or. (W .lt. 0.01)) then !no activation if aerosol < 1e-4 1/cm3           
! 
!            NCPL_ACT  = 0.0
!         else
!
!            ! Calculate constants. These fits were obtained from detailed correlations of physical properties. G is actually 1/G
!            T = min(max(TEMP, 243.0), 323.0)     
!            alfa=2.8915E-08*(T**2) - 2.1328E-05*T + 4.2523E-03
!            beta=exp(3.49996E-04*T**2 - 2.27938E-01*T + 4.20901E+01)
!            G=exp(-2.94362E-06*T**3 + 2.77941E-03*T**2 - 8.92889E-01*T + 1.18787E+02)
!            Akoh= 0.66e-6/T  !from Seinfeld and Pandis (1998)
!     
!            !=======================================================
!            !Activate droplets   
!            !=======================================================
!            !Calculate maximum supersaturation according to ARG2002
!
!            auxx=0.0 
!            
!            
!            DO INDEX = 1, NMODES            
!                
!               kappa=  max(AER_BASE%kap(INDEX), 0.001)
!             
!                  SMI (INDEX) = ((0.667*Akoh/AER_BASE%dpg(INDEX))**1.5)/SQRT(2.0*kappa)   ! Critical supersat for mode I   
!                  SMI=MAX(SMI, 1.0e-5)   
!                   
!              if ((TPI(INDEX) .gt. 1e4) .and.  (kappa .gt. 0.1)) then                       
!                  fi=0.5*exp(2.5*SIGI(INDEX)) !sigi is now log(sigi)
!                  gi=1.0+0.25*SIGI(INDEX)
!                  nui=((alfa*W*G)**1.5)/(2.0*MAPL_PI*980.0*beta*TPI(INDEX))
!                  citai = 0.667*Akoh*SQRT(alfa*W*G)
!                  aux1=fi*((citai/nui)**1.5) + gi*(SMI(INDEX)*SMI(INDEX)/(nui+(3.0*citai)))**0.75
!                  aux1=aux1/(SMI(INDEX)*SMI(INDEX))      
!                  auxx=auxx+aux1                  
!                end if
!            end do
!
!  !Calculate number of activated droplets
!            if (auxx .gt. 0.0) then
!               smax = 1/sqrt(auxx)
!               auxx=0.0
!
!                   DO INDEX = 1, NMODES
!                        if ((TPI(INDEX) .gt. 1e4) .and. (AER_BASE%kap(index) .gt. 0.1)) then
!                           ui=sqrt(2.0)*log(SMI(INDEX)/smax)/3.0
!                           aux1=0.5*TPI(INDEX)*(1.0-ERFAPP(ui))
!                           auxx=auxx+aux1
!                           AER_BASE%num(index) = max(TPI(INDEX) -aux1, 0.0) !remove already activated aerosol
!                        end if                    
!                   END DO
!                  NCPL_ACT=auxx             
!            else
!                  NCPL_ACT = 0.0             
!            end if
!
!     
! 
!       end if
!
!         !now filllup dust and soot number
!         NMODES =  AER%nmods
!
!         call getINsubset(1, AER,  auxaer)
!         CDUST(1:auxaer%nmods)= auxaer%num(1:auxaer%nmods)
!         DDUSTAMB(1:auxaer%nmods)= auxaer%dpg(1:auxaer%nmods)
!         call getINsubset(2, AER,  auxaer)
!         naux = max(auxaer%nmods, 1)
!         CSOOT= sum(auxaer%num) 
!         DSOOTAMB= sum(auxaer%dpg)/naux
!
!
!         PACT=1.0 ! fraction of entrained aerosol that is activated
!         auxconc =0.0
!         aseasalt_aux  = 0.0
!
!         do index = 1, nmodes 
!	       if (AER%kap(index) .gt. 0.8)  auxconc = AER%num(index) + auxconc
!           if (AER_BASE%kap(index) .gt. 0.8)    aseasalt_aux  = aseasalt_aux  + AER_BASE%num(index)*AER_BASE%dpg(index)*AER_BASE
!%dpg(index)*1.61*MAPL_PI !assumes a fixed sigma = 2.0
!
!         end do
!       aseasalt = max(aseasalt, aseasalt_aux)
!     
!	  NCPL_AMB=auxconc !Activate  entrained aerosol with kappa>0.8   
!
!
!
!      END SUBROUTINE ARGact
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      !=====Subroutine INfreezing=========
!      ! Freeze droplets in immersion and contact ice nucleation modes, according to Barahona et al. GMD (2014)
!!!!!!!!!!!!!!!!!!!!!!
!
!      subroutine INfreezing(QL, NL, NIN, NDUST, NSOOT, INDUST, INSOOT, TEMP, PRE, WSUB_, DDUST, DSOOT)  !NIN freezing tendency
!         REAL(8), INTENT( IN) :: TEMP, NSOOT, WSUB_, PRE, QL, NL, DSOOT
!         REAL(8), DIMENSION(NDUSTMAX), INTENT (IN) ::  NDUST, DDUST
!         REAL(8), DIMENSION(NDUSTMAX), INTENT (OUT) ::  INDUST
!         REAL(8), INTENT(OUT) ::  NIN, INSOOT 
!         REAL(8), DIMENSION(NDUSTMAX) ::  INDUST_C
!
!         REAL(8) :: a, b, c , d, Tx, n05, ui, aux, nssoot, nsdust, ninbkg, SI, acorr, &
!               dnsd, dnss, coolr, WSUB, ahet, INSOOT_C
!
!         REAL(8) :: nssoot_c, nsdust_c, mfp, nslip, ndfaer, viscosity, lam_r, taux, rho_a, &
!                fdust_drop, fsoot_drop, min_ns_dust, min_ns_soot, nsss, INsea, dnsss, min_ns_seasalt 
!
!         logical :: demott, Drop_mediated
!         integer :: ix
!
!         min_ns_dust= 3.75e6 !limits ice nucleation to -12 !new 02/10/14 
!         min_ns_soot= 3.75e9 !limits ice nucleation to -18
!         min_ns_seasalt = 4.0e2 !limits ice nucleation to -5
!         
!         demott=.false.
!         INDUST=0.0
!         INSOOT=0.0
!         INDUST_C=0.0
!         INSOOT_C=0.0   
!         NIN=0.0
!         Drop_mediated = .false.
!       INsea = 0.0 ! sea salt only in immersion
!   
!   
!   ! note for sea salt we just assume that it is equal to the current droplet concentration and take the area from the calculati
!on at cloud base
!
!         ! fraction of dust and soot within droplets
!         fdust_drop= FDROP_DUST
!         fsoot_drop = FDROP_SOOT
!
!
!         WSUB=MAX(MIN(WSUB_, 10.0), 0.8)
!         coolr=5.0e-3*WSUB  !Approximation to saturated cooling rate 
!         n05=sum(NDUST)+NSOOT
!
!         if (TEMP .gt. T_ICE_MAX) then 
!            return
!         end if
!
!         if (TEMP .lt. T_ICE_ALL   ) then 
!            return
!         end if
!
!
!         if ((QL .le. 1.0e-10) .or. (NL .le. 1.0)) then 
!            return
!         end if
!
!         !Background IN
!         ! SI at water saturation
!
!         rho_a = PRE*100.0/MAPL_RGAS/TEMP 
!         ninbkg=0.0
!         SI = -1.2379e-2+3.3595 !Ice supersat at water sat. Derived from Murphy and Koop 2005
!         !if (TEMP .lt. 260.0)  ninbkg=coolr*42.8*exp(3.88*si)*0.1 !tendency in IN from background IN. Derived from Phillips et 
!al. 2007
!
!
!         Tx = max(TEMP-273.16, -38.0 )
!
!         lam_r=min((MAPL_PI*950.0*NL/rho_a/QL)**(1./3.), 1.0e8)
!         viscosity=1.8e-5*(TEMP/298.0)**0.85    ! Viscosity (kg/m/s)
!         mfp=2.0*viscosity/(PRE  &                   ! Mean free path (m)
!               *sqrt(8.0*28.96e-3/(MAPL_PI*8.314409*TEMP)))        
!
!         if ((n05 .gt.1.0) .and. (TEMP .lt. 272.0)) then
!
!            nsdust=  max(exp(-0.517*Tx + 8.934)-min_ns_dust, 0.0) !From Niemand 2012
!            nssoot= max(1.0e4*exp(-0.0101*Tx*Tx - 0.8525*Tx + 0.7667)-min_ns_soot, 0.0) !Murray (review_ 2012)
!            dnsd  = 0.517*nsdust
!            dnss  = max(-(-2.0*0.0101*Tx -0.8525)*nssoot, 0.0)
!
!            !ns in  contact. It is assumed that in contact is T-3 immersion
!            taux=max(Tx-3.0, -35.0)
!            nsdust_c= max(exp(-0.517*taux + 8.934)-min_ns_dust, 0.0) !From Niemand 2012
!            nssoot_c= max(1.0e4*exp(-0.0101*taux*taux - 0.8525*taux + 0.7667)-min_ns_soot, 0.0) !Murray (review_ 2012)
!
!            aux=0.0        
!            acorr=2.7e7 !m2/m3 correction to the area due to non sphericity and aggregation. Assumes 10 m2/g (Murray 2011)
!
!
!            DO ix=1, NDUSTMAX
!               !Immersion
!               ahet=0.52*DDUST(ix)*DDUST(ix)*DDUST(ix)*acorr*exp(4.5*log(2.0)*log(2.0)*log(2.0))    !this needs to be improved
!
!               INDUST(ix) = NDUST(ix)*exp(-nsdust*ahet)* &
!                     dnsd*coolr*ahet*fdust_drop
!               !Contact   
!               nslip =1.0+(2.0*mfp/DDUST(ix))*(1.257+(0.4*exp(-(1.1*DDUST(ix)*0.5/mfp))))! Slip correction factor               
!     
!               ndfaer =1.381e-23*TEMP*nslip*(1.0-exp(-nsdust_c*ahet)) /(12.*MAPL_PI*viscosity*DDUST(ix))             
!               INDUST_C(ix) = 2.0*MAPL_PI*ndfaer*NDUST(ix)*NL/lam_r
!
!            END DO
!
!
!            acorr=8.0e7 !m2/m3 correction to the area due to non sphericity and aggregation  Assumes 50 m2/g (Popovicheva 2003) 
!      
!            ahet =0.52*DSOOT*DSOOT*DSOOT*acorr*exp(4.5*log(2.0)*log(2.0)*log(2.0))             
!            INSOOT=fsoot_drop*NSOOT*exp(-nssoot*ahet)*dnss*ahet*coolr !
!
!            nslip =1.0+(2.0*mfp/DSOOT)*(1.257+(0.4*exp(-(1.1*DSOOT*0.5/mfp))))! Slip correction factor                    
!            ndfaer =1.381e-23*TEMP*nslip*(1.0-exp(-nssoot_c*ahet)) /(12.*MAPL_PI*viscosity*DSOOT)             
!            INSOOT_c= 2.0*MAPL_PI*ndfaer*NSOOT*NL/lam_r
!
!         ! sea salt
!         nsss =  -0.459*TEMP +128.6235 ! from Demott et al. PNAS, 2015  
!         nsss=  max(exp(nsss)-min_ns_seasalt, 0.0)           
!    	 dnsss=  max(0.459*nsss, 0.0)
!         INsea= aseasalt*dnsss*coolr 
!
! 
!         end if
!
!	    NIN =ninbkg+ INSOOT + SUM(INDUST) + INSOOT_C + SUM(INDUST_C) + INsea!
!         INSOOT=INSOOT +INSOOT_C
!         INDUST =INDUST + INDUST_C
!
!	    
!      end subroutine INfreezing
!
!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      !===================Subroutine Qgrowth======================
!      !Partitions water and ice according to Korolev and Mazin 2005. Assume spheres by now.
!!!!!!!!!!!!!!!!!!!!!
!      subroutine Qgrowth(TEMP, PRE, QICE0, NIPRE, QT, NINUC, DQIG, RIM, FNDRIM) !freezing of IN according to Demott et al 2010 (
!everything SI)
!         REAL(8), INTENT(IN) :: TEMP, QICE0, NIPRE, NINUC,  PRE, QT
!         REAL(8), INTENT(INOUT) :: DQIG, RIM, FNDRIM
!
!         !REAL(8) :: A, denI, aux, Dco, SI, denA, DQold, DQnew 
!         REAL(8) :: DIFF, DENICE, DENAIR, K1, K2, K3, SI, AUX, DC, TEFF, PLo, TEo, TC, &
!               DQnew, DQold, rho_a, LWC, IWC, qmin_rim
!
!
!         if (TEMP .gt. 272.15) then
!            DQIG =0.0
!            return
!         end if
!
!         TC=TEMP-273.0 
!         PLo = max(PRE, 10.0) !limits  of the correlations 
!         TEo = max(190.0, TEMP)
!
!         qmin_rim = 1.0e-12
!
!
!         DENICE= 1000.0*(0.9167 - 1.75e-4*TC -5.0e-7*TC*TC) !From PK 97
!         DENAIR= PLo*100.0/MAPL_RGAS/TEMP
!         DIFF=(0.211*1013.25/(PLo+0.1))*(((TEo+0.1)/273.0)**1.94)*1e-4  !From Seinfeld and Pandis 2006
!
!         K1 = EXP(7.1170e-4*TEo*TEo-0.43563*TEo+78.744) 
!         K2 = EXP(-9.933e-3*TEo+25.26)
!         K3 = EXP(7.1772e-4*TEo*TEo-0.44055*TEo+73.996)
!
!
!         AUX= 210368.0 + 131.438*TEMP - (3.32373E6/TEMP)- (41729.1*LOG(TEMP)) !From Murphy and Koop 2005
!         SI=exp(-aux/8.314/TEMP)-1.0 !ratio of pw/pi-1
!         rho_a = PRE*100.0/MAPL_RGAS/TEMP 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!         if  ((NIPRE .gt. 1.0) .and. (QICE0 .gt. 1.0e-10)) then 
!            DC=max((QICE0/(NIPRE*500.0*MAPL_PI))**(0.333), 40.0e-6) !Assumme monodisperse size distribution about size distribut
!ion is made. 
!         else   
!            DC = 40.0e-6
!         end if
!
!         AUX=  NIPRE*DENICE*MAPL_PI*DC*DC
!         TEFF = DENAIR*2.0*((K1*DIFF+K2)*DC+(K3/0.1))
!
!         if  (AUX .gt. 1.0e-12) then 
!            TEFF=min(TEFF/AUX, 1.0e20)   
!            DQold= SI/TEFF  
!         else
!            DQold=0.0
!         end if
!
!         ! Calculate rimming fraction
!         IWC =  QICE0*rho_a
!         LWC =  max(QT-QICE0, 0.0)*rho_a
!         aux = DQold ! only due to deposition
!
!
!         !Account for rimming
!
!         if ((LWC .gt. qmin_rim)  .and. (IWC .gt. qmin_rim))  then 
!            RIM = 6.0e-5/(LWC*(IWC**0.17)) !Fom Lin and Colle, NRW, 2011
!            RIM = 1.0/(1.0+RIM)
!            RIM  = min (0.95, RIM)
!            DQold =  DQold*(1 + RIM/(1.0-RIM))
!            FNDrim =  max(min(rho_a*(DQold -aux)/LWC, 1.0), 0.0) !Fraction of liquid condensate removed due to riming   
!         END if
!
!
!
!
!!!!!!!!!!!!!!!!!!!!!!recently nucleated!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!         AUX=  NINUC*DENICE*MAPL_PI*20.0e-6*20.0e-6
!         TEFF = DENAIR*2.0*((K1*DIFF+K2)*DC+(K3/0.1))
!
!         if  (AUX .gt. 1.0e-12)  then 
!            TEFF=min(TEFF/AUX, 1.0e10)
!            DQnew= SI/TEFF
!         else
!            DQnew = 0.0
!         end if
!
!         DQIG = DQold+DQnew     
!
!      end subroutine Qgrowth
!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!
!
!      !*************************************************************
!      ! Function PDG07 (simplified background ice nucleation 
!      !                     spectra according to Phillips et. al. 2007).  
!      ! si is supersaturation wrt ice and T is in K 
!      !************************************************************  
!
!      subroutine PDG07_ice(si, Tx, N)     
!
!         REAL(8), intent(in) :: si, Tx
!         REAL(8), intent(out)  :: N 
!         N=0.0
!
!         !if (Tx .le. 243.0)then
!
!
!         N=1000.0*exp(-0.388)*(exp(3.88*si)-1.0)/0.76
!         !elseif (Tx .le. 260.0) then
!         !  N=60.0*exp(-0.639)*(exp(12.96*si)-1.0)/0.76   
!         !end if      
!      end subroutine PDG07_ice
!
!
!      !*********************************************************************
    INTEGER, INTENT(IN) :: idim, irun, k0, icmin
    REAL*8, DIMENSION(idim, k0), INTENT(INOUT) :: tho, qho, uho, vho
    REAL*8, DIMENSION(idim, k0), INTENT(INOUT) :: tho_ad
    REAL*8, DIMENSION(idim, k0 + 1), INTENT(IN) :: ple, pke
    REAL*8, DIMENSION(idim, k0), INTENT(IN) :: qss, dqs
    REAL*8, DIMENSION(idim, k0) :: qss_ad, dqs_ad
    REAL*8, DIMENSION(idim), INTENT(IN) :: cnv_fraction, rasal2_2d
    REAL*8, DIMENSION(k0 + 1), INTENT(IN) :: sige
    REAL*8, DIMENSION(idim, k0) :: clw, flxd
    REAL*8, DIMENSION(idim, k0) :: clw_ad, flxd_ad
    REAL*8, DIMENSION(idim, k0) :: cnv_prc3
    REAL*8, DIMENSION(idim, k0) :: cnv_prc3_ad
    REAL*8, DIMENSION(idim, k0) :: cnv_updfrc
    REAL*8, DIMENSION(idim, k0) :: cnv_updfrc_ad
    REAL*8, INTENT(IN) :: dt, cpo, alhlo, gravo
    REAL*8, INTENT(IN) :: alhl1, tice
    INTEGER, INTENT(IN) :: itrcr
    INTEGER, DIMENSION(idim, 2), INTENT(IN) :: seedras
    INTEGER, DIMENSION(idim), INTENT(IN) :: kcbl
    INTEGER, DIMENSION(idim), INTENT(IN) :: doras
    REAL*8, DIMENSION(idim), INTENT(IN) :: tpert, qpert
    REAL*8, DIMENSION(idim) :: tpert_ad
    REAL*8, DIMENSION(idim), INTENT(IN) :: co_auto
    REAL*8, DIMENSION(idim) :: mxdiam
    REAL*8, DIMENSION(idim, k0), INTENT(IN) :: wgt0, wgt1
    TYPE(RASPARAM_TYPE), INTENT(IN) :: rasparams
    REAL*8, OPTIONAL, INTENT(INOUT) :: xho(idim, k0, itrcr)
    REAL*8, OPTIONAL, INTENT(OUT) :: triedlev_diag(idim, k0)
    REAL*8, OPTIONAL, INTENT(OUT) :: disske(idim, k0)
    REAL*8, OPTIONAL, INTENT(IN) :: fscav(itrcr)
    REAL*8, DIMENSION(k0) :: poi_sv, qoi_sv, uoi_sv, voi_sv
    REAL*8, DIMENSION(k0) :: poi_sv_ad, qoi_sv_ad, uoi_sv_ad, voi_sv_ad
    REAL*8, DIMENSION(k0) :: poi, qoi, uoi, voi, dqq, bet, gam, cll, &
&   tmpt
    REAL*8, DIMENSION(k0) :: poi_ad, qoi_ad, uoi_ad, voi_ad, dqq_ad, &
&   bet_ad, gam_ad, cll_ad, tmpt_ad
    REAL*8, DIMENSION(k0) :: poi_c, qoi_c
    REAL*8, DIMENSION(k0) :: poi_c_ad, qoi_c_ad
    REAL*8, DIMENSION(k0) :: prh, pri, ght, dpt, dpb, pki, dissk0, &
&   dissk1, clantnd
    REAL*8, DIMENSION(k0) :: ght_ad
    REAL*8, DIMENSION(k0) :: tcu, qcu, ucu, vcu, cln, rns, pol, dm
    REAL*8, DIMENSION(k0) :: ucu_ad, vcu_ad, rns_ad
    REAL*8, DIMENSION(k0) :: qst, ssl, rmf, rnn, rn1, rmfc, rmfp
    REAL*8, DIMENSION(k0) :: qst_ad, ssl_ad, rnn_ad, rmfp_ad
    REAL*8, DIMENSION(k0) :: gms, eta, gmh, eht, gm1, hcc, rmfd
    REAL*8, DIMENSION(k0) :: gms_ad, eta_ad, gmh_ad, eht_ad, gm1_ad, &
&   hcc_ad, rmfd_ad
    REAL*8, DIMENSION(k0) :: hol, hst, qol, zol, hcld, cll0, cllx, clli&
&   , cllb
    REAL*8, DIMENSION(k0) :: hol_ad, hst_ad, qol_ad, zol_ad, hcld_ad, &
&   cll0_ad
    REAL*8, DIMENSION(k0) :: wsp, lambdsv, bke, cvw, updfrc
    REAL*8, DIMENSION(k0) :: cvw_ad, updfrc_ad
    REAL*8, DIMENSION(k0) :: rasal, mtkwi, updfrp, bk2, bk3, dll0, dllx
    REAL*8, DIMENSION(k0) :: rasal_ad, updfrp_ad, bk2_ad
    REAL*8, DIMENSION(itrcr) :: xht
    REAL*8, DIMENSION(k0, itrcr) :: xoi, xcu, xoi_sv
    REAL*8, DIMENSION(k0 + 1) :: prj, prs, qht, sht, zet, xyd, xyd0
    REAL*8, DIMENSION(k0+1) :: qht_ad, sht_ad, zet_ad
    INTEGER :: k, my_pe
    REAL*8, DIMENSION(idim, k0) :: lambdsv2
    REAL*8 :: tx2, tx3, uht, vht, akm, acr, alm, tth, qqh, shtrg, wspbl&
&   , dqx
    REAL*8 :: tx2_ad, tx3_ad, uht_ad, vht_ad, akm_ad, alm_ad, dqx_ad
    REAL*8 :: wfn, tem, trg, trgexp, evp, wlq, qcc, mtkw_max
    REAL*8 :: wfn_ad, tem_ad, trg_ad, wlq_ad, qcc_ad
    REAL*8 :: shtrg_fac, sige_minhol, wfnog
    INTEGER :: i, ic, l, icl, itr, icl_c, n_dtl
    INTEGER :: ndtlexpon
    INTEGER, DIMENSION(:), ALLOCATABLE :: icl_v
    REAL*8 :: grav, cp, alhl, cpbg, alhi, cpi, gravi, ddt, lbcp, obg, &
&   afc
    REAL*8 :: fricfac, dpth_bl, wupdrft, pblfrac, autorampb, co_zdep
    REAL*8 :: rasal1, rasal2, rasal2i, co_t, rasncl, friclambda, sdqvt1&
&   , sdqv2
    REAL*8 :: lambda_fac, strapping, acritfac, hmintrigger, lldisaggxp
    REAL*8 :: lambmx_fac, diammn_min, rdtlexpon, cli_crit, sdqv3, &
&   maxdallowed_d, maxdallowed_s
    REAL*8 :: rhmn, rhmx, cldmicro, fdrop_dust, fdrop_soot
    INTEGER :: kstrap, rasal_exp
    REAL*8 :: cld_radius, areal_frac, spect_mflx, cvw_cbase
    REAL*8, PARAMETER :: onepkap=1.+2./7., daylen=86400.0
    REAL*8, PARAMETER :: rhmax=0.9999
    REAL*8 :: lambda_min
    REAL*8 :: lambda_max
    REAL*8, PARAMETER :: rho_w=1.0e3
    LOGICAL :: dyna_strapping, do_tracers, smooth_hst
    CHARACTER(len=esmf_maxstr) :: cbl_style
    REAL*8, DIMENSION(k0) :: tcu8, qcu8, pcu, flx8
    REAL*8, DIMENSION(k0, itrcr + 2) :: rcu
    REAL*8 :: cup
    LOGICAL :: revap, wrkfun, calkpb, crtfun, lprnt, dndrft
    REAL*8, DIMENSION(k0) :: toi8, qoi8, prsm8, phil8, qli8, qii8, &
&   trcfac
    REAL*8, DIMENSION(k0) :: alfind, alfint, alfinq, rhc_ls
    REAL*8, DIMENSION(k0 + 1) :: prs8, phih8
    REAL*8, DIMENSION(k0, itrcr + 2) :: roi8
    REAL*8 :: fracbl, dt8, rasalf
    INTEGER :: kpbl
    REAL*8, SAVE :: max_neg_bouy=1.0
    REAL*8, SAVE :: rhfacl=0.0
    REAL*8, SAVE :: rhfacs=0.0
    REAL*8, SAVE :: garea=1.e10
    REAL*8, SAVE :: dsfc=0.001
    REAL*8, SAVE :: cd=1.e-3
    REAL*8, SAVE :: wfnc=0.0
    REAL*8, SAVE :: tla=-1.0
    REAL*8, SAVE :: dpd=300.
    REAL*8 :: delzkm
    REAL*8 :: fnoscav
    REAL*8 :: fscav_(itrcr)
    INTRINSIC PRESENT
    INTRINSIC INT
    INTRINSIC SQRT
    INTRINSIC SUM
    INTRINSIC ALLOCATED
    INTEGER :: ad_to
    INTEGER :: branch
    REAL*8, DIMENSION(idim, k0), INTENT(INOUT) :: vho_ad
    REAL*8, DIMENSION(idim, k0), INTENT(INOUT) :: qho_ad
    REAL*8, DIMENSION(idim, k0), INTENT(INOUT) :: uho_ad
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      dqs_ad = 0.0_8
      tpert_ad = 0.0_8
      qss_ad = 0.0_8
    ELSE
      CALL POPREAL8ARRAY(gam, k0)
      CALL POPREAL8ARRAY(rmfp, k0)
      CALL POPREAL8ARRAY(hol, k0)
      CALL POPREAL8ARRAY(zol, k0)
      CALL POPREAL8ARRAY(zet, k0 + 1)
      CALL POPREAL8ARRAY(dpt, k0)
      CALL POPREAL8ARRAY(ssl, k0)
      CALL POPREAL8ARRAY(uoi, k0)
      CALL POPINTEGER4(icl)
      CALL POPREAL8ARRAY(prs, k0 + 1)
      CALL POPREAL8ARRAY(poi, k0)
      CALL POPREAL8ARRAY(gms, k0)
      CALL POPREAL8(uht)
      CALL POPREAL8ARRAY(qst, k0)
      CALL POPREAL8ARRAY(rnn, k0)
      CALL POPREAL8(alm)
      CALL POPREAL8(lambda_min)
      CALL POPREAL8ARRAY(prj, k0 + 1)
      CALL POPREAL8ARRAY(pri, k0)
      CALL POPREAL8ARRAY(prh, k0)
      CALL POPREAL8ARRAY(dpb, k0)
      CALL POPREAL8ARRAY(gmh, k0)
      CALL POPREAL8(rasal2i)
      CALL POPINTEGER4(n_dtl)
      CALL POPREAL8ARRAY(sht, k0 + 1)
      CALL POPREAL8ARRAY(pki, k0)
      CALL POPREAL8ARRAY(eta, k0)
      CALL POPREAL8ARRAY(hcld, k0)
      CALL POPREAL8ARRAY(rasal, k0)
      CALL POPREAL8(tx3)
      CALL POPREAL8(tx2)
      CALL POPREAL8ARRAY(voi, k0)
      CALL POPREAL8(akm)
      CALL POPREAL8ARRAY(qol, k0)
      CALL POPREAL8ARRAY(qoi, k0)
      CALL POPREAL8(vht)
      CALL POPREAL8ARRAY(qht, k0 + 1)
      CALL POPREAL8ARRAY(bet, k0)
      CALL POPREAL8ARRAY(hst, k0)
      CALL POPREAL8ARRAY(ght, k0)
      CALL POPREAL8(wlq)
      CALL POPREAL8ARRAY(dqq, k0)
      CALL POPREAL8ARRAY(cvw, k0)
      CALL POPREAL8(tem)
      CALL POPREAL8(trg)
      CALL POPREAL8ARRAY(hcc, k0)
      CALL POPREAL8ARRAY(gm1, k0)
      CALL POPREAL8ARRAY(eht, k0)
      CALL POPREAL8(acr)
      fricfac = rasparams%cufricfac
      rhmn = rasparams%ras_rhmin
      sdqv2 = rasparams%sdqv2
      ddt = daylen/dt
      cp = cpo
      cpi = 1.0/cp
      rhmx = rasparams%ras_rhfull
      sdqv3 = rasparams%sdqv3
      grav = gravo
      gravi = 1.0/grav
      rasal1 = rasparams%rasal1
      alhl = alhlo
      alhi = 1.0/alhl
      rasal_exp = rasparams%rasal_exp
      cli_crit = rasparams%qc_crit_cn
      cpbg = cp*gravi
      friclambda = rasparams%cufriclambda
      lbcp = alhl*cpi
      sdqvt1 = rasparams%sdqvt1
      dqs_ad = 0.0_8
      tpert_ad = 0.0_8
      qss_ad = 0.0_8
      eht_ad = 0.0_8
      gm1_ad = 0.0_8
      hcc_ad = 0.0_8
      cvw_ad = 0.0_8
      dqq_ad = 0.0_8
      ucu_ad = 0.0_8
      ght_ad = 0.0_8
      hst_ad = 0.0_8
      bet_ad = 0.0_8
      qoi_sv_ad = 0.0_8
      qht_ad = 0.0_8
      uoi_sv_ad = 0.0_8
      qoi_ad = 0.0_8
      qol_ad = 0.0_8
      voi_ad = 0.0_8
      bk2_ad = 0.0_8
      rasal_ad = 0.0_8
      hcld_ad = 0.0_8
      eta_ad = 0.0_8
      sht_ad = 0.0_8
      gmh_ad = 0.0_8
      rnn_ad = 0.0_8
      qst_ad = 0.0_8
      gms_ad = 0.0_8
      updfrc_ad = 0.0_8
      poi_ad = 0.0_8
      rns_ad = 0.0_8
      poi_sv_ad = 0.0_8
      uoi_ad = 0.0_8
      vcu_ad = 0.0_8
      rmfd_ad = 0.0_8
      voi_sv_ad = 0.0_8
      ssl_ad = 0.0_8
      zet_ad = 0.0_8
      updfrp_ad = 0.0_8
      zol_ad = 0.0_8
      cll0_ad = 0.0_8
      tmpt_ad = 0.0_8
      cll_ad = 0.0_8
      rmfp_ad = 0.0_8
      gam_ad = 0.0_8
      DO i=irun,1,-1
        CALL POPCONTROL2B(branch)
        IF (branch .LT. 2) THEN
          IF (branch .NE. 0) THEN
            k = kcbl(i)
            CALL STRAP_BWD(final=2)
          END IF
        ELSE
          IF (branch .EQ. 2) THEN
            k = kcbl(i)
            CALL STRAP_BWD(final=2)
          ELSE
            k = kcbl(i)
            CALL STRAP_BWD(final=1)
            CALL RNEVP_BWD()
          END IF
          hol_ad = 0.0_8
          CALL POPINTEGER4(ad_to)
          DO icl_c=ad_to,1,-1
            CALL POPCONTROL1B(branch)
            IF (branch .NE. 0) CALL CLOUDE_BWD(icl)
            vcu_ad(icmin:k0) = 0.0_8
            ucu_ad(icmin:k0) = 0.0_8
            CALL POPINTEGER4(icl)
          END DO
          CALL HTEST_BWD()
          CALL STRAP_BWD(final=0)
        END IF
      END DO
    END IF

  CONTAINS
!  Differentiation of cloude in reverse (adjoint) mode, forward sweep (with options r8 split(GEOS_MoistGridComp.PRE_RASE GEOS_Moi
!stGridComp.PRE_PROGNO_CLOUD ras.RASE rase.CLOUDE rase.ACRITN rase.RNEVP rase.HTEST rase.FINDDTLS rase.STRAP ras.SUNDQ3_ICE cloud
!new.pdf_spread cloudnew.fix_up_clouds cloudnew.meltfrz cloudnew.hystpdf cloudnew.pdffrac cloudnew.pdfcondensate cloudnew.cnvsrc 
!cloudnew.evap3 cloudnew.subl3 cloudnew.autocon3 cloudnew.PRECIP3 cloudnew.ICEFALL cloudnew.SETTLE_VEL cloudnew.MARSHPALMQ2 cloud
!new.MICRO_AA_BB_3 cloudnew.LDRADIUS3 cloudnew.ICE_FRACTION cloudnew.GET_ALHX3 cloudnew.ICEFRAC cloudnew.SUNDQ3_ICE3 qsat_util.DQ
!SATPERT)):
!   gradient     of useful results: tpert eht gm1 hcc cvw dqq ucu
!                ght hst bet qht qoi qol voi bk2 rasal hcld eta
!                sht gmh rnn qst gms updfrc poi rns uoi vcu rmfd
!                ssl zet updfrp zol cll0 hol cll rmfp gam
!   with respect to varying inputs: tpert eht gm1 hcc cvw dqq ucu
!                ght hst bet qht qoi qol voi bk2 rasal hcld eta
!                sht gmh rnn qst gms updfrc poi rns uoi vcu rmfd
!                ssl zet updfrp zol cll0 hol cll rmfp gam
!*********************************************************************
    SUBROUTINE CLOUDE_FWD(ic)
      IMPLICIT NONE
!=======================================
      INTEGER, INTENT(IN) :: ic
      REAL*8 :: deep_fact, cu_diam, wscale
!, dQx
      REAL*8 :: cli, te_a, c00_x, cli_crit_x, pete, toki, gmhx, hstx
      REAL*8 :: dt_lyr, rate, cvw_x, closs, f2, f3, f4, f5
      INTEGER :: k700
      INTRINSIC MIN
      INTRINSIC MAX
      INTRINSIC SQRT
      INTRINSIC EXP
      REAL*8 :: min2
      REAL*8 :: min1
      REAL*8 :: x5
      REAL*8 :: x4
      REAL*8 :: x3
      REAL*8 :: x2
      REAL*8 :: x1
      LOGICAL :: mask(k-ic+1)
      REAL*8 :: max2
      REAL*8 :: max1
      REAL*8 :: y1
!         !=============================AER_CLOUD local variables ====================
!         REAL(8) :: WBASE, NDROP, NICE, FP_D, FF_A, FP_I, FICE, &
!               NDROP_AMB, NSOOT_AMB, NSOOT, NIN, INSOOT, dCVW2, QICE, &
!               dQICE, dQIG, FPICE, dNICE, dNDROP, DSOOT_AMB, DSOOT, QLIQ, dQLIQ, FPRECIP, AUX, QT, &
!               MAXNICE, MAXNDROP, MINNICE, MINNDROP, NDROP_ACT, RIMM, FNDRIM, TminusTa, Tparcel , &
!	      alph_e, beta_e, RH_AMB, ECRIT
!
!         REAL(8), DIMENSION (NDUSTMAX) :: NDUST, NDUST_AMB, INDUST, DDUST_AMB, DDUST 
!         INTEGER :: INX
!
!
!         T_ICE_ALL = 238.0 !A little higher so there is ice at the freezing level
!         WBASE=1.0 
!         FICE=0.0
!         NICE=0.0
!         NDROP=0.0
!         NDROP_AMB=0.0
!         NDROP_ACT=0.0
!         NIN = 0.0
!         NDUST=0.0
!         NSOOT=0.0
!         NDUST_AMB=0.0
!         NSOOT_AMB = 0.0
!         dCVW2=0.0
!         QICE=0.0
!         QLIQ=0.0
!         FPICE = 0.0
!         INDUST=0.0
!         INSOOT= 0.0
!         QT= 0.0
!         FPRECIP=0.0   
!         FNDRIM = 0.0
!         RIMM = 0.0   
!         TminusTa = 0.0
!
!         f_seasalt = 0.0
!         aseasalt = 0.0
!
!         call init_Aer(AER_BASE)
!
!         !AER_CLOUD=============================
      IF (1. .GT. (qoi(k)/qst(k)-rhmn)/(rhmx-rhmn)) THEN
        CALL PUSHREAL8(trg)
        trg = (qoi(k)/qst(k)-rhmn)/(rhmx-rhmn)
        CALL PUSHCONTROL1B(0)
      ELSE
        CALL PUSHREAL8(trg)
        trg = 1.
        CALL PUSHCONTROL1B(1)
      END IF
      IF (0.0 .LT. (autorampb-sige(ic))/0.2) THEN
        y1 = (autorampb-sige(ic))/0.2
      ELSE
        y1 = 0.0
      END IF
      IF (1.0 .GT. y1) THEN
        f4 = y1
      ELSE
        f4 = 1.0
      END IF
      IF (trg .LE. 1.0e-5) THEN
        CALL PUSHCONTROL3B(0)
      ELSE
!  RECOMPUTE SOUNDING UP TO DETRAINMENT LEVEL
        poi_c = poi
        qoi_c = qoi
        poi_c(k) = poi_c(k) + tpert(i)
        qoi_c(k) = qoi_c(k) + qpert(i)
        CALL PUSHREAL8(zet(k+1))
        zet(k+1) = 0.
        CALL PUSHREAL8(sht(k+1))
        sht(k+1) = cp*poi_c(k)*prj(k+1)
        DO l=k,ic,-1
          IF (qst(l)*rhmax .GT. qoi_c(l)) THEN
            CALL PUSHREAL8(qol(l))
            qol(l) = qoi_c(l)
            CALL PUSHCONTROL1B(0)
          ELSE
            CALL PUSHREAL8(qol(l))
            qol(l) = qst(l)*rhmax
            CALL PUSHCONTROL1B(1)
          END IF
          IF (0.000 .LT. qol(l)) THEN
            CALL PUSHCONTROL1B(0)
            qol(l) = qol(l)
          ELSE
            qol(l) = 0.000
            CALL PUSHCONTROL1B(1)
          END IF
          CALL PUSHREAL8(ssl(l))
          ssl(l) = cp*prj(l+1)*poi_c(l) + grav*zet(l+1)
          CALL PUSHREAL8(hol(l))
          hol(l) = ssl(l) + qol(l)*alhl
          CALL PUSHREAL8(hst(l))
          hst(l) = ssl(l) + qst(l)*alhl
          CALL PUSHREAL8(tem)
          tem = poi_c(l)*(prj(l+1)-prj(l))*cpbg
          CALL PUSHREAL8(zet(l))
          zet(l) = zet(l+1) + tem
          CALL PUSHREAL8(zol(l))
          zol(l) = zet(l+1) + (prj(l+1)-prh(l))*poi_c(l)*cpbg
        END DO
        DO l=ic+1,k
          CALL PUSHREAL8(tem)
          tem = (prj(l)-prh(l-1))/(prh(l)-prh(l-1))
          CALL PUSHREAL8(sht(l))
          sht(l) = ssl(l-1) + tem*(ssl(l)-ssl(l-1))
          CALL PUSHREAL8(qht(l))
          qht(l) = .5*(qol(l)+qol(l-1))
        END DO
! SMOOTH HSTAR W/ 1-2-1 Filter
        IF (smooth_hst) THEN
! save for later
          hstx = hst(ic)
          DO l=k-1,ic+1,-1
            CALL PUSHREAL8(hst(l))
            hst(l) = 0.25*(hst(l+1)+hst(l-1)) + 0.5*hst(l)
          END DO
          DO l=ic,ic
            CALL PUSHREAL8(hst(l))
            hst(l) = 0.5*hst(l+1) + 0.5*hst(l)
          END DO
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
!  CALCULATE LAMBDA, ETA, AND WORKFUNCTION
        CALL PUSHREAL8(lambda_min)
        lambda_min = .2/mxdiam(i)
        lambda_max = .2/200.
!     LAMBDA_MIN = .2/(LAMBDA_FAC*DPTH_BL)
!     LAMBDA_MAX = .2/( MAX( LAMBMX_FAC*DPTH_BL , DIAMMN_MIN ) )
        IF (hol(k) .LE. hst(ic)) THEN
          CALL PUSHCONTROL3B(1)
        ELSE
!  LAMBDA CALCULATION: MS-A18
          CALL PUSHREAL8(tem)
          tem = (hst(ic)-hol(ic))*(zol(ic)-zet(ic+1))
          DO l=ic+1,k-1
            tem = tem + (hst(ic)-hol(l))*(zet(l)-zet(l+1))
          END DO
          IF (tem .LE. 0.0) THEN
            CALL PUSHCONTROL3B(2)
          ELSE
            CALL PUSHREAL8(alm)
            alm = (hol(k)-hst(ic))/tem
            IF (alm .GT. lambda_max) THEN
              CALL PUSHCONTROL3B(3)
            ELSE
              toki = 1.0
              IF (alm .LT. lambda_min) THEN
!we can probably replace this by a actual distribution based on grid cell size
                toki = (alm/lambda_min)**2
!RC(IC) = 6
!RETURN
                CALL PUSHCONTROL1B(1)
              ELSE
                CALL PUSHCONTROL1B(0)
              END IF
!LAMBDSV(IC) = ALM
!  ETA CALCULATION: MS-A2
              DO l=ic+1,k
                CALL PUSHREAL8(eta(l))
                eta(l) = 1.0 + alm*(zet(l)-zet(k))
              END DO
              CALL PUSHREAL8(eta(ic))
              eta(ic) = 1.0 + alm*(zol(ic)-zet(k))
!  WORKFUNCTION CALCULATION:  MS-A22
              wfn = 0.0
              CALL PUSHREAL8(hcc(k))
              hcc(k) = hol(k)
              DO l=k-1,ic+1,-1
                CALL PUSHREAL8(hcc(l))
                hcc(l) = hcc(l+1) + (eta(l)-eta(l+1))*hol(l)
                CALL PUSHREAL8(tem)
                tem = hcc(l+1)*dpb(l) + hcc(l)*dpt(l)
                CALL PUSHREAL8(eht(l))
                eht(l) = eta(l+1)*dpb(l) + eta(l)*dpt(l)
                wfn = wfn + (tem-eht(l)*hst(l))*gam(l)
              END DO
              CALL PUSHREAL8(hcc(ic))
              hcc(ic) = hst(ic)*eta(ic)
              wfn = wfn + (hcc(ic+1)-hst(ic)*eta(ic+1))*gam(ic)*dpb(ic)
!  VERTICAL VELOCITY/KE CALCULATION (ADDED 12/2001 JTB)
              bk3(k) = 0.0
              bk2(k) = 0.0
              bke(k) = 0.0
              CALL PUSHREAL8(hcld(k))
              hcld(k) = hol(k)
              DO l=k-1,ic,-1
                CALL PUSHREAL8(hcld(l))
                hcld(l) = (eta(l+1)*hcld(l+1)+(eta(l)-eta(l+1))*hol(l))/&
&                 eta(l)
                CALL PUSHREAL8(tem)
                tem = (hcld(l)-hst(l))*(zet(l)-zet(l+1))/(1.0+lbcp*dqq(l&
&                 ))
                IF (cldmicro .LE. 0.0) THEN
                  bke(l) = bke(l+1) + grav*tem/(cp*prj(l+1)*poi(l))
                  IF (tem .LT. 0.0) THEN
                    CALL PUSHREAL8(max1)
                    max1 = 0.0
                    CALL PUSHCONTROL1B(0)
                  ELSE
                    CALL PUSHREAL8(max1)
                    max1 = tem
                    CALL PUSHCONTROL1B(1)
                  END IF
                  bk2(l) = bk2(l+1) + grav*max1/(cp*prj(l+1)*poi(l))
                  IF (tem .GT. 0.0) THEN
                    min1 = 0.0
                  ELSE
                    min1 = tem
                  END IF
                  bk3(l) = bk3(l+1) + grav*min1/(cp*prj(l+1)*poi(l))
                  IF (bk2(l) .LT. 0.0) THEN
                    CALL PUSHREAL8(max2)
                    max2 = 0.0
                    CALL PUSHCONTROL1B(0)
                  ELSE
                    CALL PUSHREAL8(max2)
                    max2 = bk2(l)
                    CALL PUSHCONTROL1B(1)
                  END IF
                  CALL PUSHREAL8(cvw(l))
                  cvw(l) = SQRT(2.0*max2)
                  CALL PUSHCONTROL1B(1)
                ELSE
                  CALL PUSHCONTROL1B(0)
                END IF
              END DO
! 1.0 / ( 5.0*ALM )
!   ALPHA CALCULATION 
              CALL PUSHREAL8(rasal2i)
              rasal2i = rasal2_2d(i)
              IF (zet(ic) .LT. 2000.) THEN
                CALL PUSHREAL8(rasal(ic))
                rasal(ic) = rasal1
                CALL PUSHCONTROL1B(0)
              ELSE
                CALL PUSHCONTROL1B(1)
              END IF
              IF (zet(ic) .GE. 2000.) THEN
                IF (1.0 .GT. (zet(ic)-2000.)/8000.) THEN
                  min2 = (zet(ic)-2000.)/8000.
                  CALL PUSHCONTROL1B(0)
                ELSE
                  CALL PUSHCONTROL1B(1)
                  min2 = 1.0
                END IF
!WMP     RASAL(IC) = RASAL1 + (RASAL2i-RASAL1)*(ZET(IC) - 2000.)/8000.
                CALL PUSHREAL8(rasal(ic))
                rasal(ic) = rasal1 + (rasal2i-rasal1)*min2**rasal_exp
                CALL PUSHCONTROL1B(1)
              ELSE
                CALL PUSHCONTROL1B(0)
              END IF
!WMP  RASAL(IC) = MIN( RASAL(IC) , 1.0e5 )
              CALL PUSHREAL8(rasal(ic))
              rasal(ic) = dt/rasal(ic)
              mask(1:k-ic+1) = cvw(ic:k) .LT. 1.00
              CALL PUSHREAL8ARRAY(cvw(ic:k), k - ic + 1)
              WHERE (mask(1:k-ic+1)) cvw(ic:k) = 1.00
              CALL PUSHREAL8ARRAY(cvw(ic:k), k - ic + 1)
              WHERE (.NOT.mask(1:k-ic+1)) cvw(ic:k) = cvw(ic:k)
!  NOTE THIS "CENTRALIZES" A KLUGE PRESENT IN OTHER LOCATIONS.
!  CLEAN UP SOME TIME.      -JTB 12/04/03
!  TEST FOR CRITICAL WORK FUNCTION
              CALL ACRITN(pol(ic), prs(k), acr)
              IF (wfn .LE. acr) THEN
                CALL PUSHREAL8(max1)
                CALL PUSHREAL8(max2)
                CALL PUSHBOOLEANARRAY(mask, k - ic + 1)
                CALL PUSHREAL8(min2)
                CALL PUSHCONTROL3B(4)
              ELSE
!  CLOUD TOP WATER AND MOMENTUM (TIMES ETA(IC)) MS-A16
! Tracer scavenging
! RAS loops over a series of plumes all having common cloud base level K 
! and different detrainment levels IC.  The plumes operate sequentially
! on the grid box mean quantities (wind, moisture, tracer) and so each
! subsequent plume is seeing the effects of previous plumes.  We parameterize
! scavenging following Liu et al. [JGR, 2001], their equation 1:
!  AEROSOL FRACTION SCAVENGED = 1 - exp(-FSCAV*DZ)
! where FSCAV is a specified scavenging efficiency [km-1] and DZ is the
! distance [km] the tracer traverses in the plume from it's entrainment
! level to its detrainment level.  We write the aerosol fraction surviving as:
!  FNOSCAV = exp(- FSCAV_(ITR) * DZ)
! The total scavenging is proportional to the convective mass flux, which
! is not explicitly solved for at this point.
                IF (do_tracers) THEN
                  DO itr=1,itrcr
!           Scavenging of the below cloud tracer
                    delzkm = (zet(ic)-zet(k))/1000.
                    x4 = EXP(-(fscav_(itr)*delzkm))
                    IF (x4 .GT. 1.) THEN
                      x1 = 1.
                    ELSE
                      x1 = x4
                    END IF
                    IF (x1 .LT. 0.) THEN
                      fnoscav = 0.
                    ELSE
                      fnoscav = x1
                    END IF
                    xht(itr) = xoi(k, itr)*fnoscav
                  END DO
                  CALL PUSHCONTROL1B(0)
                ELSE
                  CALL PUSHCONTROL1B(1)
                END IF
                CALL PUSHREAL8(wlq)
                wlq = qol(k)
                CALL PUSHREAL8(uht)
                uht = uoi(k)
                CALL PUSHREAL8(vht)
                vht = voi(k)
                CALL PUSHREAL8(rnn(k))
                rnn(k) = 0.
                cll0(k) = 0.
!print *, '========================================='
                DO l=k-1,ic,-1
                  CALL PUSHREAL8(tem)
                  tem = eta(l) - eta(l+1)
                  wlq = wlq + tem*qol(l)
                  uht = uht + tem*uoi(l)
                  vht = vht + tem*voi(l)
                  IF (do_tracers) THEN
                    DO itr=1,itrcr
!         Scavenging of the entrained tracer.  Updates transported tracer mass.
                      delzkm = (zet(ic)-zet(l+1))/1000.
                      x5 = EXP(-(fscav_(itr)*delzkm))
                      IF (x5 .GT. 1.) THEN
                        x2 = 1.
                      ELSE
                        x2 = x5
                      END IF
                      IF (x2 .LT. 0.) THEN
                        fnoscav = 0.
                      ELSE
                        fnoscav = x2
                      END IF
                      xht(itr) = xht(itr) + tem*xoi(l, itr)*fnoscav
                    END DO
                    CALL PUSHCONTROL1B(0)
                  ELSE
                    CALL PUSHCONTROL1B(1)
                  END IF
!!!! How much condensate (CLI) is present here? 
                  IF (l .GT. ic) THEN
                    CALL PUSHREAL8(tx2)
                    tx2 = 0.5*(qst(l)+qst(l-1))*eta(l)
                    CALL PUSHREAL8(tx3)
                    tx3 = 0.5*(hst(l)+hst(l-1))*eta(l)
                    qcc = tx2 + gm1(l)*(hcc(l)-tx3)
                    cll0(l) = wlq - qcc
                    CALL PUSHCONTROL1B(1)
                  ELSE
                    cll0(l) = wlq - qst(ic)*eta(ic)
                    CALL PUSHCONTROL1B(0)
                  END IF
                  IF (cll0(l) .LT. 0.00) THEN
                    cll0(l) = 0.00
                    CALL PUSHCONTROL1B(0)
                  ELSE
                    CALL PUSHCONTROL1B(1)
                    cll0(l) = cll0(l)
                  END IF
! condensate (kg/kg)
                  cli = cll0(l)/eta(l)
! Temperature (K)
                  te_a = poi(l)*prh(l)
!=====================================================================
                  IF (cldmicro .GT. 0.0) THEN
                    CALL PUSHCONTROL1B(1)
                  ELSE
!AER_CLOUD MIcrophysics considering activation and nucleation 
!               !recompute vertical velocity
!
!               Tparcel = TE_A
!               CVW(K) = 0.8  ! Assume a below cloud base  W of 0.8 m s-1            
!               BK2(K)   = 0.0
!
!     
!               TEM     = (HCLD(L)-HST(L) )/ (1.0+LBCP*DQQ(L))  
!               TminusTa = max(min(TEM/CP, 5.0), 0.0) !limit DT to 5 K. According to Wei, JAS, 1998   
!	     TEM =0.33*TminusTa*CO_AUTO(I)/TE_A !Bouyancy term, effciency =0.5 mwr Roode et al    	     
!
!               BK2(L)  = BK2(L+1) + GRAV * TEM*(ZET(L)-ZET(L+1)) 
!               BK2(L) = BK2(L) - (ZET(L)-ZET(L+1))*(BK2(L+1)*ALM + CLI*GRAV)  !Account for drag from entrainment of stagnat air 
!and condesate loading
!               CVW(L) = max(SQRT(  2.0* MAX( BK2(L) , 0.0 )  ), 1.0) 
!
!
!	    CVW_X = MIN(CVW(L), 50.0)
!               DT_LYR  =  max(( ZET(L)-ZET(L+1) )/CVW_X, 1.0) !Sanity check 
!               TEM   = ETA(L) - ETA(L+1)
!
!               Tparcel  =  TE_A + TminusTa
!
!
!
!!!!!!!!!!account for entrainment effects on activation !!!!!!!!!!!
!               ! Barahona and Nenes, JGR, 2007
!               alph_e = 2.8915e-8*Tparcel*Tparcel -2.1328e-5*Tparcel+4.2523e-3
!               beta_e = MAPL_ALHL*TminusTa/MAPL_RVAP/Tparcel/Tparcel
!               RH_AMB=QOI(L)/QST(L)
!               ECRIT  = max(RH_AMB -beta_e, 1.0e-6) 
!               ECRIT =  alph_e/ECRIT
!               ! print *, L, Tparcel, RH_AMB, ECRIT, ALM
!	           ECRIT =  ALM/ECRIT
!       ! ECRIT = 0.0 ! do not use this for now
!               !Print *, ECRIT
!
!
!               if (L .eq. K-1) then
!
!                  FICE=0.0
!                  NICE=0.0
!                  NDROP=0.0
!                  NIN =0.0
!                  NDUST_AMB =0.0
!                  NSOOT_AMB = 0.0
!                  NSOOT=0.0
!                  NDUST= 0.0
!
!        
!                  AER_BASE =  AERO(L)
!                  RATE=0.0
!                  FPRECIP=0.0
!                  !initial conditions     
!                  call ARGact(Tparcel, CVW_X, NDROP_ACT, NDROP_AMB, NDUST_AMB, NSOOT_AMB, L,  .true., DDUST_AMB, DSOOT_AMB, ECRI
!T) !cloud droplet number and INsource at cloud base 
!                  NDUST=NDUST_AMB
!                  NSOOT=NSOOT_AMB
!                  DDUST=DDUST_AMB
!                  DSOOT=DSOOT_AMB                                     
!
!               else 
!                  call ARGact(Tparcel, CVW_X, NDROP_ACT, NDROP_AMB, NDUST_AMB, NSOOT_AMB, L, .false., DDUST_AMB, DSOOT_AMB, ECRI
!T) !cloud droplet number above cloud base  
!
!               end if
!
!               QT = CLI
!               RATE = 0.0
!               FPRECIP = 0.0
!
!               if (QT .gt. 0.0) then
!
!                  ! if FICE is already >= 1.0 then the cloud is glaciated and there is no need to do anymore partitioning
!
!                  if (FICE .ge. 1.0) then
!
!
!                     CALL  Qremoval(RATE, FICE, FP_D, FP_I, Tparcel,  & 
!                           POL(L), QT,  NICE, NDROP, CVW_X, FPICE, &
!                           DT_LYR, RIMM, CO_AUTO(I)) 
!
!
!
!                     dNICE = -NICE*FP_I 
!                     NICE  =  (NICE +dNICE*DT_LYR)*ETA(L+1)/ETA(L) !ice
!
!                     MINNICE = max(QICE*(1.0-RATE*DT_LYR), 0.0)/4.0e-8!assuming maximum vol radius 250 microns
!                     MAXNICE = max(QICE*(1.0-RATE*DT_LYR), 0.0)/2.51e-12 !assuming minimum vol radius 10 microns
!
!                     NICE=MIN(max(NICE, MINNICE), MAXNICE)
!
!                     FICE = 1.0
!
!                  else 
!
!                     ! Cloud is not completely glaciated do the whole thing
!                     ! ALL this subroutines return tendencies
!
!
!                     CALL  INfreezing(QLIQ, NDROP, NIN, NDUST, NSOOT, INDUST, INSOOT, Tparcel, POL(L), CVW_X, DDUST, DSOOT)  !ca
!lculate the freezing fraction of the aerosol at this level
!
!                     NIN = min(NIN, NDROP/DT_LYR)
!
!                     call Qgrowth(Tparcel, POL(L), QICE, NICE, QT, NIN, dQIG, RIMM, FNDRIM)
!
!                     CALL  Qremoval(RATE, FICE, FP_D, FP_I, Tparcel,  & 
!                           POL(L), QT,  NICE, NDROP, CVW_X, FPICE, &
!                           DT_LYR,  RIMM, CO_AUTO(I)) 
!
!
!
!                     !ice number tendency: -precip + freezin
!                     dNICE = -NICE*FP_I  + NIN     
!                     NICE  =  (NICE +dNICE*DT_LYR)*ETA(L+1)/ETA(L) !ice
!                     NICE =max(NICE, 0.0)
!
!
!                     !ice mass tendency: growth - precip
!                     dQICE = -QICE*FPICE + dQIG
!                     QICE  =  min((QICE + dQICE*DT_LYR)*ETA(L+1)/ETA(L), QT) !ice
!                     QICE=max(min(QT, QICE), 0.0)
!
!
!                     ! Liquid Tendency: source/evap -  precip 
!                     !dQLIQ = max((CLI-QICE), -QLIQ)/DT_LYR -QLIQ*max(RATE-FPICE, 0.0) 
!                     ! dQLIQ = CLI*(1.0-RATE*DT_LYR)/DT_LYR -dQICE - QLIQ*max(RATE-FPICE, 0.0)           
!                     !QLIQ  =  max((QLIQ + dQLIQ*DT_LYR)*ETA(L+1)/ETA(L), 0.0) !liquid. This is actually diagnostic
!                     QLIQ=max((QT-QICE), 0.0)
!
!
!                     !droplet number tendency: -precip - freezin + activation + activated entrained aerosol 
!
!
!                     dNDROP =-NDROP*FP_D - NIN -  FNDRIM*NDROP/DT_LYR + max(NDROP_ACT-NDROP, 0.0)/DT_LYR          
!
!                     !dNDROP =-NDROP*FP_D - NIN -  FNDRIM*NDROP/DT_LYR + NDROP_ACT/DT_LYR         
!
!                     NDROP =  (NDROP + dNDROP*DT_LYR)*ETA(L+1)/ETA(L) + &
!                           (ZET(L) - ZET(L+1))*ALM*MAX((NDROP_AMB-NDROP), 0.0)
!
!                     !Aerosol tendency: Entrainment - freezing 
!
!                     NDUST = (NDUST - INDUST*DT_LYR)*ETA(L+1)/ETA(L) + &
!                           (ZET(L) - ZET(L+1))*ALM*MAX(NDUST_AMB-NDUST, 0.0) 
!
!                     NSOOT =  (NSOOT - INSOOT*DT_LYR)*ETA(L+1)/ETA(L)  + &     
!                           (ZET(L) - ZET(L+1))*ALM*MAX(NSOOT_AMB-NSOOT, 0.0)  
!
!
!                           
!                     !Update FICE and perform Sanity checks
!
!
!                     MINNDROP = (1.0-FICE)*QLIQ*max(1.0-RATE*DT_LYR, 0.0)/2.e-10    !assuming maximum vol radius 36 microns
!                     MAXNDROP = (1.0-FICE)*QLIQ*max(1.0-RATE*DT_LYR, 0.0)/3.35e-14 !assuming minimum vol radius 2 microns
!                     MINNICE = QICE/4.0e-8!assuming maximum vol radius 250 microns
!                     MAXNICE = QICE/2.51e-12 !assuming minimum vol radius 10 microns
!
!                     IF ((NICE .gt. MAXNICE) .or. (NICE .lt. MINNICE))   then    
!                        !print *, 'nilim', NICE*1e-6, MINNICE*1e-6, MAXNICE*1e-6
!                     END IF
!
!                     IF ((NDROP .gt. MAXNDROP) .or. (NDROP .lt. MINNDROP))      then 
!                        !print *, 'ndroplim', NDROP*1e-6, MINNDROP*1e-6, MAXNDROP*1e-6
!                     end if
!
!
!                     NSOOT=MAX(NSOOT, 0.0)
!                     NDUST=MAX(NDUST, 0.0)              
!
!                     NDROP=MIN(max(NDROP, MINNDROP), MAXNDROP)
!                     NICE=MIN(max(NICE, MINNICE), MAXNICE)
!
!                     FICE=max(min(QICE/QT, 1.0), 0.0)
!
!                     IF (FICE .ge. 1.0) THEN !Complete glaciation 
!                        NICE=NICE+NDROP 
!                        NDROP = 0.0
!                        QICE  = QT
!                        QLIQ= 0.0
!                     END IF
!
!                     IF (Tparcel .LT. T_ICE_ALL) THEN !instantaneous freezing
!                        NICE=NICE+NDROP 
!                        NDROP = 0.0
!                        FICE  = 1.0
!                        QICE  = QT
!                        QLIQ=0.0
!                     END IF
!
!                     IF (Tparcel .GT. T_ICE_MAX) THEN !instantaneous melting
!                        NDROP=NICE+NDROP 
!                        NICE = 0.0
!                        FICE  = 0.0
!                        QICE  = 0.0
!                        QLIQ=QT
!                     END IF
!
!                  END IF
!
!               else 
!
!                  FICE =0.0 
!                  QICE = 0.0
!                  QLIQ = 0.0
!                  NICE= 0.0 
!                  NDROP = 0.0
!                  RATE =0.0
!               end if
!
!               FPRECIP= RATE*DT_LYR
!
!               !RATE=RATE*F4
!               ! NDROP=NDROP*F4
!               !NICE=NICE*(1.0-F4)
!
!               !print *, TE_A, FICE, 'NICE', NICE*1e-6, 'NDROP', NDROP*1e-6, L 
!               !print *, 'FPI', FP_I*DT_LYR, 'FPD', FP_D*DT_LYR, 'FPICE', FPICE, 'FPRE', FPRECIP, QT, QLIQ
!
!Bacmeister 2006 microphysics
                    CALL SUNDQ3_ICE_FWD(te_a, sdqv2, sdqv3, sdqvt1, f2, &
&                                 f3)
! * F5  ! F4 reduces AUTO for shallow clouds, F5 modifies auto for deep clouds
                    CALL PUSHREAL8(c00_x)
                    c00_x = co_auto(i)*f2*f3*f4
                    cli_crit_x = cli_crit/(f2*f3)
                    CALL PUSHREAL8(rate)
                    rate = c00_x*(1.0-EXP(-(cli**2/cli_crit_x**2)))
                    CALL PUSHCONTROL1B(0)
                  END IF
                  IF (cvw(l) .LT. 1.00) THEN
                    CALL PUSHREAL8(cvw_x)
                    cvw_x = 1.00
                    CALL PUSHCONTROL1B(0)
                  ELSE
                    CALL PUSHREAL8(cvw_x)
                    cvw_x = cvw(l)
                    CALL PUSHCONTROL1B(1)
                  END IF
! really trust it at low values
! l.h.s. DT_LYR => time in layer (L,L+1)
                  dt_lyr = (zet(l)-zet(l+1))/cvw_x
                  closs = cll0(l)*rate*dt_lyr
                  IF (closs .GT. cll0(l)) THEN
                    closs = cll0(l)
                    CALL PUSHCONTROL1B(0)
                  ELSE
                    closs = closs
                    CALL PUSHCONTROL1B(1)
                  END IF
                  CALL PUSHREAL8(cll0(l))
                  cll0(l) = cll0(l) - closs
                  dll0(l) = closs
                  IF (closs .GT. 0.) THEN
                    wlq = wlq - closs
                    CALL PUSHREAL8(rnn(l))
                    rnn(l) = closs
                    CALL PUSHCONTROL1B(1)
                  ELSE
                    CALL PUSHREAL8(rnn(l))
                    rnn(l) = 0.
                    CALL PUSHCONTROL1B(0)
                  END IF
                END DO
!AER_CLOUD=======================================
!            CNVNDROP(IC)=NDROP
!            CNVNICE(IC)=NICE
!            CNVFICE(IC)=FICE
                wlq = wlq - qst(ic)*eta(ic)
!     CALCULATE GAMMAS AND KERNEL
! MS-A30 (W/O GRAV)
                CALL PUSHREAL8(gms(k))
                gms(k) = (sht(k)-ssl(k))*pri(k)
! MS-A31 (W/O GRAV)
                CALL PUSHREAL8(gmh(k))
                gmh(k) = gms(k) + (qht(k)-qol(k))*pri(k)*alhl
! MS-A37 (W/O GRAV)
                CALL PUSHREAL8(akm)
                akm = gmh(k)*gam(k-1)*dpb(k-1)
                CALL PUSHREAL8(tx2)
                tx2 = gmh(k)
                DO l=k-1,ic+1,-1
                  CALL PUSHREAL8(gms(l))
                  gms(l) = (eta(l)*(sht(l)-ssl(l))+eta(l+1)*(ssl(l)-sht(&
&                   l+1)))*pri(l)
                  CALL PUSHREAL8(gmh(l))
                  gmh(l) = gms(l) + (eta(l)*(qht(l)-qol(l))+eta(l+1)*(&
&                   qol(l)-qht(l+1)))*alhl*pri(l)
                  CALL PUSHREAL8(tx2)
                  tx2 = tx2 + (eta(l)-eta(l+1))*gmh(l)
                  akm = akm - gms(l)*eht(l)*pki(l) + tx2*ght(l)
                END DO
                CALL PUSHREAL8(gms(ic))
                gms(ic) = eta(ic+1)*(ssl(ic)-sht(ic+1))*pri(ic)
                akm = akm - gms(ic)*eta(ic+1)*dpb(ic)*pki(ic)
                CALL PUSHREAL8(gmh(ic))
                gmh(ic) = gms(ic) + (eta(ic+1)*(qol(ic)-qht(ic+1))*alhl+&
&                 eta(ic)*(hst(ic)-hol(ic)))*pri(ic)
                IF (smooth_hst) THEN
                  gmhx = gms(ic) + (eta(ic+1)*(qol(ic)-qht(ic+1))*alhl+&
&                   eta(ic)*(hstx-hol(ic)))*pri(ic)
                  CALL PUSHCONTROL1B(0)
                ELSE
                  CALL PUSHCONTROL1B(1)
                END IF
!    CLOUD BASE MASS FLUX
                IF (akm .GE. 0.0 .OR. wlq .LT. 0.0) THEN
                  CALL PUSHREAL8(max1)
                  CALL PUSHREAL8(max2)
                  CALL PUSHBOOLEANARRAY(mask, k - ic + 1)
                  CALL PUSHREAL8(f2)
                  CALL PUSHREAL8(f3)
                  CALL PUSHREAL8(f4)
                  CALL PUSHREAL8(c00_x)
                  CALL PUSHREAL8(cvw_x)
                  CALL PUSHREAL8(min2)
                  CALL PUSHREAL8(rate)
                  CALL PUSHREAL8(hstx)
                  CALL PUSHCONTROL3B(5)
                ELSE
! MS-A39 MASS-FLUX IN Pa/step
                  CALL PUSHREAL8(wfn)
                  wfn = -((wfn-acr)/akm)
                  x3 = rasal(ic)*trg*toki*wfn
                  IF (x3 .GT. (prs(k+1)-prs(k))*(100.*pblfrac)) THEN
                    CALL PUSHREAL8(wfn)
                    wfn = (prs(k+1)-prs(k))*(100.*pblfrac)
                    CALL PUSHCONTROL1B(0)
                  ELSE
                    CALL PUSHREAL8(wfn)
                    wfn = x3
                    CALL PUSHCONTROL1B(1)
                  END IF
!    CUMULATIVE PRECIP AND CLOUD-BASE MASS FLUX FOR OUTPUT
                  wfnog = wfn*gravi
                  CALL PUSHREAL8(tem)
                  tem = wfn*gravi
! (kg/m^2/step)
                  cll(ic) = cll(ic) + wlq*tem
! (kg/m^2/step)
                  rmf(ic) = rmf(ic) + tem
! (kg/m^2/step)
                  rmfd(ic) = rmfd(ic) + tem*eta(ic)
                  DO l=ic+1,k
! (kg/m^2/step)
                    CALL PUSHREAL8(rmfp(l))
                    rmfp(l) = tem*eta(l)
! (kg/m^2/step)
                    rmfc(l) = rmfc(l) + rmfp(l)
                    dllx(l) = dllx(l) + tem*dll0(l)
                    IF (cvw(l) .GT. 0.0) THEN
                      updfrp(l) = rmfp(l)*(ddt/daylen)*1000./(cvw(l)*prs&
&                       (l))
                      CALL PUSHCONTROL1B(0)
                    ELSE
                      updfrp(l) = 0.0
                      CALL PUSHCONTROL1B(1)
                    END IF
! current cloud; incloud condensate        
                    clli(l) = cll0(l)/eta(l)
!  cumulative grid mean convective condensate        
                    cllb(l) = cllb(l) + updfrp(l)*clli(l)
                    updfrc(l) = updfrc(l) + updfrp(l)
                  END DO
!    THETA AND Q CHANGE DUE TO CLOUD TYPE IC
                  DO l=ic,k
! (kg/m^2/step)
                    rns(l) = rns(l) + rnn(l)*tem
                    CALL PUSHREAL8(gmh(l))
                    gmh(l) = gmh(l)*wfn
                    CALL PUSHREAL8(gms(l))
                    gms(l) = gms(l)*wfn
                    CALL PUSHREAL8(qoi(l))
                    qoi(l) = qoi(l) + (gmh(l)-gms(l))*alhi
                    CALL PUSHREAL8(poi(l))
                    poi(l) = poi(l) + gms(l)*pki(l)*cpi
                    CALL PUSHREAL8(qst(l))
                    qst(l) = qst(l) + gms(l)*bet(l)*cpi
                  END DO
                  IF (smooth_hst) THEN
                    CALL PUSHREAL8(gmhx)
                    gmhx = gmhx*wfn
                    dqx = (gmhx-gmh(ic))*alhi
                    rns(ic) = rns(ic) + dqx/(pri(ic)*grav)
                    CALL PUSHCONTROL1B(0)
                  ELSE
                    CALL PUSHCONTROL1B(1)
                  END IF
                  IF (do_tracers) THEN
!*FRICFAC*0.5
                    CALL PUSHREAL8(wfn)
                    wfn = wfn*0.5*1.0
                    CALL PUSHREAL8(tem)
                    tem = wfn*pri(k)
                    CALL PUSHCONTROL1B(0)
                    DO itr=1,itrcr
                      xcu(k, itr) = xcu(k, itr) + tem*(xoi(k-1, itr)-xoi&
&                       (k, itr))
                    END DO
                    DO itr=1,itrcr
                      DO l=k-1,ic+1,-1
                        tem = wfn*pri(l)
                        xcu(l, itr) = xcu(l, itr) + tem*((xoi(l-1, itr)-&
&                         xoi(l, itr))*eta(l)+(xoi(l, itr)-xoi(l+1, itr)&
&                         )*eta(l+1))
                      END DO
                    END DO
                    tem = wfn*pri(ic)
                    DO itr=1,itrcr
                      xcu(ic, itr) = xcu(ic, itr) + (2.*(xht(itr)-xoi(ic&
&                       , itr)*(eta(ic)-eta(ic+1)))-(xoi(ic, itr)+xoi(ic&
&                       +1, itr))*eta(ic+1))*tem
                    END DO
                    DO itr=1,itrcr
                      DO l=ic,k
                        xoi(l, itr) = xoi(l, itr) + xcu(l, itr)
                      END DO
                    END DO
                  ELSE
!*FRICFAC*0.5
                    CALL PUSHREAL8(wfn)
                    wfn = wfn*0.5*1.0
                    CALL PUSHCONTROL1B(1)
                  END IF
!   CUMULUS FRICTION
                  IF (fricfac .LE. 0.0) THEN
                    CALL PUSHREAL8(max1)
                    CALL PUSHREAL8(max2)
                    CALL PUSHBOOLEANARRAY(mask, k - ic + 1)
                    CALL PUSHREAL8(f2)
                    CALL PUSHREAL8(f3)
                    CALL PUSHREAL8(f4)
                    CALL PUSHREAL8(c00_x)
                    CALL PUSHREAL8(cvw_x)
                    CALL PUSHREAL8(toki)
                    CALL PUSHREAL8(min2)
                    CALL PUSHREAL8(rate)
                    CALL PUSHREAL8(hstx)
                    CALL PUSHCONTROL3B(6)
                  ELSE
                    CALL PUSHREAL8(wfn)
                    wfn = wfn*fricfac*EXP(-(alm/friclambda))
                    CALL PUSHREAL8(tem)
                    tem = wfn*pri(k)
                    ucu(k) = ucu(k) + tem*(uoi(k-1)-uoi(k))
                    vcu(k) = vcu(k) + tem*(voi(k-1)-voi(k))
                    DO l=k-1,ic+1,-1
                      tem = wfn*pri(l)
                      ucu(l) = ucu(l) + tem*((uoi(l-1)-uoi(l))*eta(l)+(&
&                       uoi(l)-uoi(l+1))*eta(l+1))
                      vcu(l) = vcu(l) + tem*((voi(l-1)-voi(l))*eta(l)+(&
&                       voi(l)-voi(l+1))*eta(l+1))
                    END DO
                    tem = wfn*pri(ic)
                    ucu(ic) = ucu(ic) + (2.*(uht-uoi(ic)*(eta(ic)-eta(ic&
&                     +1)))-(uoi(ic)+uoi(ic+1))*eta(ic+1))*tem
                    vcu(ic) = vcu(ic) + (2.*(vht-voi(ic)*(eta(ic)-eta(ic&
&                     +1)))-(voi(ic)+voi(ic+1))*eta(ic+1))*tem
                    dissk0(ic) = eta(ic)*grav*wfnog*pri(ic)*0.5*((uht/&
&                     eta(ic)-uoi(ic))**2+(vht/eta(ic)-voi(ic))**2)
                    DO l=ic,k
                      CALL PUSHREAL8(uoi(l))
                      uoi(l) = uoi(l) + ucu(l)
                      CALL PUSHREAL8(voi(l))
                      voi(l) = voi(l) + vcu(l)
                    END DO
                    CALL PUSHREAL8(max1)
                    CALL PUSHREAL8(max2)
                    CALL PUSHBOOLEANARRAY(mask, k - ic + 1)
                    CALL PUSHREAL8(f2)
                    CALL PUSHREAL8(f3)
                    CALL PUSHREAL8(f4)
                    CALL PUSHREAL8(c00_x)
                    CALL PUSHREAL8(cvw_x)
                    CALL PUSHREAL8(toki)
                    CALL PUSHREAL8(min2)
                    CALL PUSHREAL8(rate)
                    CALL PUSHREAL8(hstx)
                    CALL PUSHCONTROL3B(7)
                  END IF
                END IF
              END IF
            END IF
          END IF
        END IF
      END IF
    END SUBROUTINE CLOUDE_FWD
!  Differentiation of cloude in reverse (adjoint) mode, backward sweep (with options r8 split(GEOS_MoistGridComp.PRE_RASE GEOS_Mo
!istGridComp.PRE_PROGNO_CLOUD ras.RASE rase.CLOUDE rase.ACRITN rase.RNEVP rase.HTEST rase.FINDDTLS rase.STRAP ras.SUNDQ3_ICE clou
!dnew.pdf_spread cloudnew.fix_up_clouds cloudnew.meltfrz cloudnew.hystpdf cloudnew.pdffrac cloudnew.pdfcondensate cloudnew.cnvsrc
! cloudnew.evap3 cloudnew.subl3 cloudnew.autocon3 cloudnew.PRECIP3 cloudnew.ICEFALL cloudnew.SETTLE_VEL cloudnew.MARSHPALMQ2 clou
!dnew.MICRO_AA_BB_3 cloudnew.LDRADIUS3 cloudnew.ICE_FRACTION cloudnew.GET_ALHX3 cloudnew.ICEFRAC cloudnew.SUNDQ3_ICE3 qsat_util.D
!QSATPERT)):
!   gradient     of useful results: tpert eht gm1 hcc cvw dqq ucu
!                ght hst bet qht qoi qol voi bk2 rasal hcld eta
!                sht gmh rnn qst gms updfrc poi rns uoi vcu rmfd
!                ssl zet updfrp zol cll0 hol cll rmfp gam
!   with respect to varying inputs: tpert eht gm1 hcc cvw dqq ucu
!                ght hst bet qht qoi qol voi bk2 rasal hcld eta
!                sht gmh rnn qst gms updfrc poi rns uoi vcu rmfd
!                ssl zet updfrp zol cll0 hol cll rmfp gam
!*********************************************************************
    SUBROUTINE CLOUDE_BWD(ic)
      IMPLICIT NONE
!=======================================
      INTEGER, INTENT(IN) :: ic
      REAL*8 :: deep_fact, cu_diam, wscale
      REAL*8 :: cli, te_a, c00_x, cli_crit_x, pete, toki, gmhx, hstx
      REAL*8 :: cli_ad, te_a_ad, c00_x_ad, cli_crit_x_ad, toki_ad, &
&     gmhx_ad, hstx_ad
      REAL*8 :: dt_lyr, rate, cvw_x, closs, f2, f3, f4, f5
      REAL*8 :: dt_lyr_ad, rate_ad, cvw_x_ad, closs_ad, f2_ad
      INTEGER :: k700
      INTRINSIC MIN
      INTRINSIC MAX
      INTRINSIC SQRT
      INTRINSIC EXP
      INTEGER :: branch
      REAL*8 :: temp3
      REAL*8 :: temp2
      REAL*8 :: temp1
      REAL*8 :: temp0
      REAL*8 :: temp_ad39
      REAL*8 :: temp_ad38
      REAL*8 :: temp_ad37
      REAL*8 :: temp_ad36
      REAL*8 :: min2
      REAL*8 :: temp_ad35
      REAL*8 :: min1
      REAL*8 :: temp_ad34
      REAL*8 :: temp_ad33
      REAL*8 :: temp_ad32
      REAL*8 :: temp_ad31
      REAL*8 :: temp_ad30
      REAL*8 :: x5
      REAL*8 :: x4
      REAL*8 :: x3
      REAL*8 :: x2
      REAL*8 :: x1
      REAL*8 :: x3_ad
      REAL*8 :: temp_ad29
      REAL*8 :: temp_ad28
      REAL*8 :: temp_ad27
      REAL*8 :: temp_ad26
      REAL*8 :: temp_ad25
      REAL*8 :: temp_ad24
      REAL*8 :: temp_ad23
      REAL*8 :: temp_ad22
      REAL*8 :: max1_ad
      REAL*8 :: temp_ad21
      REAL*8 :: temp_ad20
      REAL*8 :: temp_ad9
      REAL*8 :: temp_ad8
      REAL*8 :: temp_ad7
      REAL*8 :: temp_ad6
      REAL*8 :: temp_ad5
      REAL*8 :: temp_ad4
      REAL*8 :: temp_ad3
      REAL*8 :: temp_ad2
      REAL*8 :: max2_ad
      REAL*8 :: temp_ad1
      REAL*8 :: temp_ad0
      REAL*8 :: min2_ad
      REAL*8 :: temp_ad
      REAL*8 :: temp_ad19
      REAL*8 :: temp_ad18
      REAL*8 :: temp_ad17
      REAL*8 :: temp_ad16
      REAL*8 :: temp_ad15
      REAL*8 :: temp_ad14
      REAL*8 :: temp_ad13
      REAL*8 :: temp_ad12
      REAL*8 :: temp_ad11
      REAL*8 :: temp_ad10
      LOGICAL :: mask(k-ic+1)
      REAL*8 :: temp
      REAL*8 :: max2
      REAL*8 :: max1
      REAL*8 :: temp6
      REAL*8 :: temp5
      REAL*8 :: y1
      REAL*8 :: temp4
      CALL POPCONTROL3B(branch)
      IF (branch .LT. 4) THEN
        IF (branch .LT. 2) THEN
          IF (branch .EQ. 0) THEN
            trg_ad = 0.0_8
            GOTO 130
          ELSE
            trg_ad = 0.0_8
            hstx_ad = 0.0_8
            GOTO 120
          END IF
        ELSE IF (branch .EQ. 2) THEN
          trg_ad = 0.0_8
          tem_ad = 0.0_8
          hstx_ad = 0.0_8
          GOTO 110
        ELSE
          trg_ad = 0.0_8
          alm_ad = 0.0_8
          hstx_ad = 0.0_8
        END IF
      ELSE
        IF (branch .LT. 6) THEN
          IF (branch .EQ. 4) THEN
            CALL POPREAL8(min2)
            CALL POPBOOLEANARRAY(mask, k - ic + 1)
            CALL POPREAL8(max2)
            CALL POPREAL8(max1)
            trg_ad = 0.0_8
            wfn_ad = 0.0_8
            alm_ad = 0.0_8
            toki_ad = 0.0_8
            hstx_ad = 0.0_8
            GOTO 100
          ELSE
            CALL POPREAL8(hstx)
            CALL POPREAL8(rate)
            CALL POPREAL8(min2)
            CALL POPREAL8(cvw_x)
            CALL POPREAL8(c00_x)
            CALL POPREAL8(f4)
            CALL POPREAL8(f3)
            CALL POPREAL8(f2)
            CALL POPBOOLEANARRAY(mask, k - ic + 1)
            CALL POPREAL8(max2)
            CALL POPREAL8(max1)
            trg_ad = 0.0_8
            wlq_ad = 0.0_8
            wfn_ad = 0.0_8
            vht_ad = 0.0_8
            akm_ad = 0.0_8
            alm_ad = 0.0_8
            uht_ad = 0.0_8
            toki_ad = 0.0_8
            gmhx_ad = 0.0_8
          END IF
        ELSE
          IF (branch .EQ. 6) THEN
            CALL POPREAL8(hstx)
            CALL POPREAL8(rate)
            CALL POPREAL8(min2)
            CALL POPREAL8(toki)
            CALL POPREAL8(cvw_x)
            CALL POPREAL8(c00_x)
            CALL POPREAL8(f4)
            CALL POPREAL8(f3)
            CALL POPREAL8(f2)
            CALL POPBOOLEANARRAY(mask, k - ic + 1)
            CALL POPREAL8(max2)
            CALL POPREAL8(max1)
            wfn_ad = 0.0_8
            vht_ad = 0.0_8
            alm_ad = 0.0_8
            uht_ad = 0.0_8
          ELSE
            CALL POPREAL8(hstx)
            CALL POPREAL8(rate)
            CALL POPREAL8(min2)
            CALL POPREAL8(toki)
            CALL POPREAL8(cvw_x)
            CALL POPREAL8(c00_x)
            CALL POPREAL8(f4)
            CALL POPREAL8(f3)
            CALL POPREAL8(f2)
            CALL POPBOOLEANARRAY(mask, k - ic + 1)
            CALL POPREAL8(max2)
            CALL POPREAL8(max1)
            DO l=k,ic,-1
              CALL POPREAL8(voi(l))
              vcu_ad(l) = vcu_ad(l) + voi_ad(l)
              CALL POPREAL8(uoi(l))
              ucu_ad(l) = ucu_ad(l) + uoi_ad(l)
            END DO
            temp_ad34 = tem*vcu_ad(ic)
            temp_ad35 = 2.*temp_ad34
            temp_ad36 = -(eta(ic+1)*temp_ad34)
            vht_ad = temp_ad35
            voi_ad(ic) = voi_ad(ic) + temp_ad36 - (eta(ic)-eta(ic+1))*&
&             temp_ad35
            eta_ad(ic) = eta_ad(ic) - voi(ic)*temp_ad35
            eta_ad(ic+1) = eta_ad(ic+1) + voi(ic)*temp_ad35 - (voi(ic)+&
&             voi(ic+1))*temp_ad34
            voi_ad(ic+1) = voi_ad(ic+1) + temp_ad36
            tem_ad = (2.*(uht-uoi(ic)*(eta(ic)-eta(ic+1)))-(uoi(ic)+uoi(&
&             ic+1))*eta(ic+1))*ucu_ad(ic) + (2.*(vht-voi(ic)*(eta(ic)-&
&             eta(ic+1)))-(voi(ic)+voi(ic+1))*eta(ic+1))*vcu_ad(ic)
            temp_ad37 = tem*ucu_ad(ic)
            temp_ad38 = 2.*temp_ad37
            temp_ad39 = -(eta(ic+1)*temp_ad37)
            uht_ad = temp_ad38
            uoi_ad(ic) = uoi_ad(ic) + temp_ad39 - (eta(ic)-eta(ic+1))*&
&             temp_ad38
            eta_ad(ic) = eta_ad(ic) - uoi(ic)*temp_ad38
            eta_ad(ic+1) = eta_ad(ic+1) + uoi(ic)*temp_ad38 - (uoi(ic)+&
&             uoi(ic+1))*temp_ad37
            uoi_ad(ic+1) = uoi_ad(ic+1) + temp_ad39
            wfn_ad = pri(ic)*tem_ad
            DO l=ic+1,k-1,1
              tem = wfn*pri(l)
              temp_ad30 = tem*vcu_ad(l)
              temp_ad31 = eta(l+1)*temp_ad30
              tem_ad = ((uoi(l-1)-uoi(l))*eta(l)+(uoi(l)-uoi(l+1))*eta(l&
&               +1))*ucu_ad(l) + ((voi(l-1)-voi(l))*eta(l)+(voi(l)-voi(l&
&               +1))*eta(l+1))*vcu_ad(l)
              voi_ad(l-1) = voi_ad(l-1) + eta(l)*temp_ad30
              voi_ad(l) = voi_ad(l) + temp_ad31 - eta(l)*temp_ad30
              eta_ad(l) = eta_ad(l) + (voi(l-1)-voi(l))*temp_ad30
              voi_ad(l+1) = voi_ad(l+1) - temp_ad31
              eta_ad(l+1) = eta_ad(l+1) + (voi(l)-voi(l+1))*temp_ad30
              temp_ad32 = tem*ucu_ad(l)
              temp_ad33 = eta(l+1)*temp_ad32
              uoi_ad(l-1) = uoi_ad(l-1) + eta(l)*temp_ad32
              uoi_ad(l) = uoi_ad(l) + temp_ad33 - eta(l)*temp_ad32
              eta_ad(l) = eta_ad(l) + (uoi(l-1)-uoi(l))*temp_ad32
              uoi_ad(l+1) = uoi_ad(l+1) - temp_ad33
              eta_ad(l+1) = eta_ad(l+1) + (uoi(l)-uoi(l+1))*temp_ad32
              wfn_ad = wfn_ad + pri(l)*tem_ad
            END DO
            tem = wfn*pri(k)
            tem_ad = (uoi(k-1)-uoi(k))*ucu_ad(k) + (voi(k-1)-voi(k))*&
&             vcu_ad(k)
            voi_ad(k-1) = voi_ad(k-1) + tem*vcu_ad(k)
            voi_ad(k) = voi_ad(k) - tem*vcu_ad(k)
            uoi_ad(k-1) = uoi_ad(k-1) + tem*ucu_ad(k)
            uoi_ad(k) = uoi_ad(k) - tem*ucu_ad(k)
            CALL POPREAL8(tem)
            wfn_ad = wfn_ad + pri(k)*tem_ad
            CALL POPREAL8(wfn)
            alm_ad = -(EXP(-(alm/friclambda))*wfn*fricfac*wfn_ad/&
&             friclambda)
            wfn_ad = fricfac*EXP(-(alm/friclambda))*wfn_ad
          END IF
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            CALL POPREAL8(tem)
            CALL POPREAL8(wfn)
            wfn_ad = 0.5*wfn_ad
          ELSE
            CALL POPREAL8(wfn)
            wfn_ad = 0.5*wfn_ad
          END IF
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            dqx_ad = rns_ad(ic)/(pri(ic)*grav)
            gmhx_ad = alhi*dqx_ad
            gmh_ad(ic) = gmh_ad(ic) - alhi*dqx_ad
            CALL POPREAL8(gmhx)
            wfn_ad = wfn_ad + gmhx*gmhx_ad
            gmhx_ad = wfn*gmhx_ad
          ELSE
            gmhx_ad = 0.0_8
          END IF
          tem_ad = 0.0_8
          DO l=k,ic,-1
            CALL POPREAL8(qst(l))
            gms_ad(l) = gms_ad(l) + pki(l)*cpi*poi_ad(l) - alhi*qoi_ad(l&
&             ) + cpi*bet(l)*qst_ad(l)
            bet_ad(l) = bet_ad(l) + cpi*gms(l)*qst_ad(l)
            CALL POPREAL8(poi(l))
            CALL POPREAL8(qoi(l))
            gmh_ad(l) = gmh_ad(l) + alhi*qoi_ad(l)
            CALL POPREAL8(gms(l))
            CALL POPREAL8(gmh(l))
            wfn_ad = wfn_ad + gmh(l)*gmh_ad(l) + gms(l)*gms_ad(l)
            gms_ad(l) = wfn*gms_ad(l)
            gmh_ad(l) = wfn*gmh_ad(l)
            rnn_ad(l) = rnn_ad(l) + tem*rns_ad(l)
            tem_ad = tem_ad + rnn(l)*rns_ad(l)
          END DO
          DO l=k,ic+1,-1
            updfrp_ad(l) = updfrp_ad(l) + updfrc_ad(l)
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              temp6 = daylen*prs(l)*cvw(l)
              temp_ad29 = ddt*1000.*updfrp_ad(l)/temp6
              rmfp_ad(l) = rmfp_ad(l) + temp_ad29
              cvw_ad(l) = cvw_ad(l) - rmfp(l)*daylen*prs(l)*temp_ad29/&
&               temp6
              updfrp_ad(l) = 0.0_8
            ELSE
              updfrp_ad(l) = 0.0_8
            END IF
            CALL POPREAL8(rmfp(l))
            tem_ad = tem_ad + eta(l)*rmfp_ad(l)
            eta_ad(l) = eta_ad(l) + tem*rmfp_ad(l)
            rmfp_ad(l) = 0.0_8
          END DO
          tem_ad = tem_ad + wlq*cll_ad(ic) + eta(ic)*rmfd_ad(ic)
          eta_ad(ic) = eta_ad(ic) + tem*rmfd_ad(ic)
          wlq_ad = tem*cll_ad(ic)
          CALL POPREAL8(tem)
          wfn_ad = wfn_ad + gravi*tem_ad
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            CALL POPREAL8(wfn)
            x3_ad = 0.0_8
          ELSE
            CALL POPREAL8(wfn)
            x3_ad = wfn_ad
          END IF
          temp_ad27 = trg*toki*x3_ad
          temp_ad28 = rasal(ic)*wfn*x3_ad
          rasal_ad(ic) = rasal_ad(ic) + wfn*temp_ad27
          wfn_ad = rasal(ic)*temp_ad27
          trg_ad = toki*temp_ad28
          toki_ad = trg*temp_ad28
          CALL POPREAL8(wfn)
          akm_ad = (wfn-acr)*wfn_ad/akm**2
          wfn_ad = -(wfn_ad/akm)
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          temp_ad25 = pri(ic)*gmhx_ad
          temp_ad26 = alhl*eta(ic+1)*temp_ad25
          gms_ad(ic) = gms_ad(ic) + gmhx_ad
          eta_ad(ic+1) = eta_ad(ic+1) + alhl*(qol(ic)-qht(ic+1))*&
&           temp_ad25
          qol_ad(ic) = qol_ad(ic) + temp_ad26
          qht_ad(ic+1) = qht_ad(ic+1) - temp_ad26
          eta_ad(ic) = eta_ad(ic) + (hstx-hol(ic))*temp_ad25
          hstx_ad = eta(ic)*temp_ad25
          hol_ad(ic) = hol_ad(ic) - eta(ic)*temp_ad25
        ELSE
          hstx_ad = 0.0_8
        END IF
        temp_ad23 = -(dpb(ic)*pki(ic)*akm_ad)
        CALL POPREAL8(gmh(ic))
        temp_ad21 = pri(ic)*gmh_ad(ic)
        temp_ad22 = alhl*eta(ic+1)*temp_ad21
        gms_ad(ic) = gms_ad(ic) + eta(ic+1)*temp_ad23 + gmh_ad(ic)
        eta_ad(ic+1) = eta_ad(ic+1) + alhl*(qol(ic)-qht(ic+1))*temp_ad21
        qol_ad(ic) = qol_ad(ic) + temp_ad22
        qht_ad(ic+1) = qht_ad(ic+1) - temp_ad22
        eta_ad(ic) = eta_ad(ic) + (hst(ic)-hol(ic))*temp_ad21
        hst_ad(ic) = hst_ad(ic) + eta(ic)*temp_ad21
        hol_ad(ic) = hol_ad(ic) - eta(ic)*temp_ad21
        gmh_ad(ic) = 0.0_8
        eta_ad(ic+1) = eta_ad(ic+1) + pri(ic)*(ssl(ic)-sht(ic+1))*gms_ad&
&         (ic) + gms(ic)*temp_ad23
        CALL POPREAL8(gms(ic))
        temp_ad24 = pri(ic)*eta(ic+1)*gms_ad(ic)
        ssl_ad(ic) = ssl_ad(ic) + temp_ad24
        sht_ad(ic+1) = sht_ad(ic+1) - temp_ad24
        gms_ad(ic) = 0.0_8
        tx2_ad = 0.0_8
        DO l=ic+1,k-1,1
          tx2_ad = tx2_ad + ght(l)*akm_ad
          gmh_ad(l) = gmh_ad(l) + (eta(l)-eta(l+1))*tx2_ad
          gms_ad(l) = gms_ad(l) + gmh_ad(l) - pki(l)*eht(l)*akm_ad
          eht_ad(l) = eht_ad(l) - pki(l)*gms(l)*akm_ad
          ght_ad(l) = ght_ad(l) + tx2*akm_ad
          CALL POPREAL8(tx2)
          eta_ad(l) = eta_ad(l) + gmh(l)*tx2_ad
          eta_ad(l+1) = eta_ad(l+1) - gmh(l)*tx2_ad
          CALL POPREAL8(gmh(l))
          temp_ad17 = alhl*pri(l)*gmh_ad(l)
          temp_ad18 = eta(l+1)*temp_ad17
          eta_ad(l) = eta_ad(l) + (qht(l)-qol(l))*temp_ad17
          qht_ad(l) = qht_ad(l) + eta(l)*temp_ad17
          qol_ad(l) = qol_ad(l) + temp_ad18 - eta(l)*temp_ad17
          eta_ad(l+1) = eta_ad(l+1) + (qol(l)-qht(l+1))*temp_ad17
          qht_ad(l+1) = qht_ad(l+1) - temp_ad18
          gmh_ad(l) = 0.0_8
          CALL POPREAL8(gms(l))
          temp_ad19 = pri(l)*gms_ad(l)
          temp_ad20 = eta(l+1)*temp_ad19
          eta_ad(l) = eta_ad(l) + (sht(l)-ssl(l))*temp_ad19
          sht_ad(l) = sht_ad(l) + eta(l)*temp_ad19
          ssl_ad(l) = ssl_ad(l) + temp_ad20 - eta(l)*temp_ad19
          eta_ad(l+1) = eta_ad(l+1) + (ssl(l)-sht(l+1))*temp_ad19
          sht_ad(l+1) = sht_ad(l+1) - temp_ad20
          gms_ad(l) = 0.0_8
        END DO
        temp_ad15 = dpb(k-1)*akm_ad
        CALL POPREAL8(tx2)
        gmh_ad(k) = gmh_ad(k) + gam(k-1)*temp_ad15 + tx2_ad
        CALL POPREAL8(akm)
        gam_ad(k-1) = gam_ad(k-1) + gmh(k)*temp_ad15
        CALL POPREAL8(gmh(k))
        temp_ad16 = pri(k)*alhl*gmh_ad(k)
        gms_ad(k) = gms_ad(k) + gmh_ad(k)
        qht_ad(k) = qht_ad(k) + temp_ad16
        qol_ad(k) = qol_ad(k) - temp_ad16
        gmh_ad(k) = 0.0_8
        CALL POPREAL8(gms(k))
        sht_ad(k) = sht_ad(k) + pri(k)*gms_ad(k)
        ssl_ad(k) = ssl_ad(k) - pri(k)*gms_ad(k)
        gms_ad(k) = 0.0_8
        qst_ad(ic) = qst_ad(ic) - eta(ic)*wlq_ad
        eta_ad(ic) = eta_ad(ic) - qst(ic)*wlq_ad
        f2_ad = 0.0_8
        rate_ad = 0.0_8
        DO l=ic,k-1,1
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            CALL POPREAL8(rnn(l))
            rnn_ad(l) = 0.0_8
            closs_ad = 0.0_8
          ELSE
            CALL POPREAL8(rnn(l))
            closs_ad = rnn_ad(l) - wlq_ad
            rnn_ad(l) = 0.0_8
          END IF
          CALL POPREAL8(cll0(l))
          closs_ad = closs_ad - cll0_ad(l)
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            dt_lyr = (zet(l)-zet(l+1))/cvw_x
            cll0_ad(l) = cll0_ad(l) + closs_ad
            closs_ad = 0.0_8
          ELSE
            dt_lyr = (zet(l)-zet(l+1))/cvw_x
          END IF
          cll0_ad(l) = cll0_ad(l) + rate*dt_lyr*closs_ad
          rate_ad = rate_ad + cll0(l)*dt_lyr*closs_ad
          dt_lyr_ad = cll0(l)*rate*closs_ad
          temp_ad14 = dt_lyr_ad/cvw_x
          zet_ad(l) = zet_ad(l) + temp_ad14
          zet_ad(l+1) = zet_ad(l+1) - temp_ad14
          cvw_x_ad = -((zet(l)-zet(l+1))*temp_ad14/cvw_x)
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            CALL POPREAL8(cvw_x)
          ELSE
            CALL POPREAL8(cvw_x)
            cvw_ad(l) = cvw_ad(l) + cvw_x_ad
          END IF
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            c00_x = co_auto(i)*f2*f3*f4
            cli_crit_x = cli_crit/(f2*f3)
            cli = cll0(l)/eta(l)
            CALL POPREAL8(rate)
            temp5 = cli_crit_x**2
            temp4 = cli**2/temp5
            temp_ad13 = EXP(-temp4)*c00_x*rate_ad/temp5
            c00_x_ad = (1.0-EXP(-temp4))*rate_ad
            cli_ad = 2*cli*temp_ad13
            cli_crit_x_ad = -(temp4*2*cli_crit_x*temp_ad13)
            f2_ad = f2_ad + f3*f4*co_auto(i)*c00_x_ad - cli_crit*&
&             cli_crit_x_ad/(f3*f2**2)
            CALL POPREAL8(c00_x)
            te_a = poi(l)*prh(l)
            CALL SUNDQ3_ICE_BWD(te_a, te_a_ad, sdqv2, sdqv3, sdqvt1, f2&
&                         , f2_ad, f3)
            rate_ad = 0.0_8
          ELSE
            te_a_ad = 0.0_8
            cli_ad = 0.0_8
          END IF
          poi_ad(l) = poi_ad(l) + prh(l)*te_a_ad
          temp_ad12 = cli_ad/eta(l)
          cll0_ad(l) = cll0_ad(l) + temp_ad12
          eta_ad(l) = eta_ad(l) - cll0(l)*temp_ad12/eta(l)
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) cll0_ad(l) = 0.0_8
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            wlq_ad = wlq_ad + cll0_ad(l)
            qst_ad(ic) = qst_ad(ic) - eta(ic)*cll0_ad(l)
            eta_ad(ic) = eta_ad(ic) - qst(ic)*cll0_ad(l)
            cll0_ad(l) = 0.0_8
          ELSE
            wlq_ad = wlq_ad + cll0_ad(l)
            qcc_ad = -cll0_ad(l)
            cll0_ad(l) = 0.0_8
            tx2_ad = qcc_ad
            gm1_ad(l) = gm1_ad(l) + (hcc(l)-tx3)*qcc_ad
            hcc_ad(l) = hcc_ad(l) + gm1(l)*qcc_ad
            tx3_ad = -(gm1(l)*qcc_ad)
            CALL POPREAL8(tx3)
            temp_ad10 = 0.5*eta(l)*tx3_ad
            hst_ad(l) = hst_ad(l) + temp_ad10
            hst_ad(l-1) = hst_ad(l-1) + temp_ad10
            eta_ad(l) = eta_ad(l) + 0.5*(qst(l)+qst(l-1))*tx2_ad + 0.5*(&
&             hst(l)+hst(l-1))*tx3_ad
            CALL POPREAL8(tx2)
            temp_ad11 = 0.5*eta(l)*tx2_ad
            qst_ad(l) = qst_ad(l) + temp_ad11
            qst_ad(l-1) = qst_ad(l-1) + temp_ad11
          END IF
          tem = eta(l) - eta(l+1)
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) itr = 0
          tem_ad = uoi(l)*uht_ad + qol(l)*wlq_ad + voi(l)*vht_ad
          voi_ad(l) = voi_ad(l) + tem*vht_ad
          uoi_ad(l) = uoi_ad(l) + tem*uht_ad
          qol_ad(l) = qol_ad(l) + tem*wlq_ad
          CALL POPREAL8(tem)
          eta_ad(l) = eta_ad(l) + tem_ad
          eta_ad(l+1) = eta_ad(l+1) - tem_ad
        END DO
        cll0_ad(k) = 0.0_8
        CALL POPREAL8(rnn(k))
        rnn_ad(k) = 0.0_8
        CALL POPREAL8(vht)
        voi_ad(k) = voi_ad(k) + vht_ad
        CALL POPREAL8(uht)
        uoi_ad(k) = uoi_ad(k) + uht_ad
        CALL POPREAL8(wlq)
        qol_ad(k) = qol_ad(k) + wlq_ad
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) itr = 0
 100    CALL POPREAL8ARRAY(cvw(ic:k), k - ic + 1)
        CALL POPREAL8ARRAY(cvw(ic:k), k - ic + 1)
        CALL POPREAL8(rasal(ic))
        rasal_ad(ic) = -(dt*rasal_ad(ic)/rasal(ic)**2)
        WHERE (mask(1:k-ic+1)) cvw_ad(ic:k) = 0.0_8
        CALL POPCONTROL1B(branch)
        IF (branch .NE. 0) THEN
          CALL POPREAL8(rasal(ic))
          IF (min2 .LE. 0.0_8 .AND. (rasal_exp .EQ. 0.0_8 .OR. rasal_exp&
&             .NE. INT(rasal_exp))) THEN
            min2_ad = 0.0
          ELSE
            min2_ad = (rasal2i-rasal1)*rasal_exp*min2**(rasal_exp-1)*&
&             rasal_ad(ic)
          END IF
          rasal_ad(ic) = 0.0_8
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) zet_ad(ic) = zet_ad(ic) + min2_ad/8000.
        END IF
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          CALL POPREAL8(rasal(ic))
          rasal_ad(ic) = 0.0_8
        END IF
        CALL POPREAL8(rasal2i)
        DO l=ic,k-1,1
          CALL POPCONTROL1B(branch)
          IF (branch .EQ. 0) THEN
            tem_ad = 0.0_8
          ELSE
            CALL POPREAL8(cvw(l))
            IF (2.0*max2 .EQ. 0.0_8) THEN
              max2_ad = 0.0
            ELSE
              max2_ad = cvw_ad(l)/SQRT(2.0*max2)
            END IF
            cvw_ad(l) = 0.0_8
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              CALL POPREAL8(max2)
            ELSE
              CALL POPREAL8(max2)
              bk2_ad(l) = bk2_ad(l) + max2_ad
            END IF
            temp3 = cp*prj(l+1)
            temp2 = temp3*poi(l)
            temp_ad9 = grav*bk2_ad(l)/temp2
            bk2_ad(l+1) = bk2_ad(l+1) + bk2_ad(l)
            max1_ad = temp_ad9
            poi_ad(l) = poi_ad(l) - max1*temp3*temp_ad9/temp2
            bk2_ad(l) = 0.0_8
            CALL POPCONTROL1B(branch)
            IF (branch .EQ. 0) THEN
              CALL POPREAL8(max1)
              tem_ad = 0.0_8
            ELSE
              CALL POPREAL8(max1)
              tem_ad = max1_ad
            END IF
          END IF
          CALL POPREAL8(tem)
          temp1 = lbcp*dqq(l) + 1.0
          temp_ad7 = tem_ad/temp1
          temp0 = zet(l) - zet(l+1)
          temp = hcld(l) - hst(l)
          hcld_ad(l) = hcld_ad(l) + temp0*temp_ad7
          hst_ad(l) = hst_ad(l) - temp0*temp_ad7
          zet_ad(l) = zet_ad(l) + temp*temp_ad7
          zet_ad(l+1) = zet_ad(l+1) - temp*temp_ad7
          dqq_ad(l) = dqq_ad(l) - temp*temp0*lbcp*temp_ad7/temp1
          CALL POPREAL8(hcld(l))
          temp_ad8 = hcld_ad(l)/eta(l)
          eta_ad(l+1) = eta_ad(l+1) + (hcld(l+1)-hol(l))*temp_ad8
          hcld_ad(l+1) = hcld_ad(l+1) + eta(l+1)*temp_ad8
          eta_ad(l) = eta_ad(l) + (hol(l)-(eta(l+1)*hcld(l+1)+(eta(l)-&
&           eta(l+1))*hol(l))/eta(l))*temp_ad8
          hol_ad(l) = hol_ad(l) + (eta(l)-eta(l+1))*temp_ad8
          hcld_ad(l) = 0.0_8
        END DO
        CALL POPREAL8(hcld(k))
        hol_ad(k) = hol_ad(k) + hcld_ad(k)
        hcld_ad(k) = 0.0_8
        bk2_ad(k) = 0.0_8
        temp_ad6 = dpb(ic)*gam(ic)*wfn_ad
        hcc_ad(ic+1) = hcc_ad(ic+1) + temp_ad6
        hst_ad(ic) = hst_ad(ic) + eta(ic)*hcc_ad(ic) - eta(ic+1)*&
&         temp_ad6
        eta_ad(ic+1) = eta_ad(ic+1) - hst(ic)*temp_ad6
        gam_ad(ic) = gam_ad(ic) + dpb(ic)*(hcc(ic+1)-hst(ic)*eta(ic+1))*&
&         wfn_ad
        CALL POPREAL8(hcc(ic))
        eta_ad(ic) = eta_ad(ic) + hst(ic)*hcc_ad(ic)
        hcc_ad(ic) = 0.0_8
        DO l=ic+1,k-1,1
          tem = hcc(l+1)*dpb(l) + hcc(l)*dpt(l)
          temp_ad5 = gam(l)*wfn_ad
          tem_ad = temp_ad5
          eht_ad(l) = eht_ad(l) - hst(l)*temp_ad5
          hst_ad(l) = hst_ad(l) - eht(l)*temp_ad5
          gam_ad(l) = gam_ad(l) + (tem-eht(l)*hst(l))*wfn_ad
          CALL POPREAL8(eht(l))
          eta_ad(l+1) = eta_ad(l+1) + dpb(l)*eht_ad(l)
          CALL POPREAL8(tem)
          hcc_ad(l+1) = hcc_ad(l+1) + dpb(l)*tem_ad
          hcc_ad(l) = hcc_ad(l) + dpt(l)*tem_ad
          CALL POPREAL8(hcc(l))
          hcc_ad(l+1) = hcc_ad(l+1) + hcc_ad(l)
          eta_ad(l) = eta_ad(l) + hol(l)*hcc_ad(l) + dpt(l)*eht_ad(l)
          eht_ad(l) = 0.0_8
          eta_ad(l+1) = eta_ad(l+1) - hol(l)*hcc_ad(l)
          hol_ad(l) = hol_ad(l) + (eta(l)-eta(l+1))*hcc_ad(l)
          hcc_ad(l) = 0.0_8
        END DO
        CALL POPREAL8(hcc(k))
        hol_ad(k) = hol_ad(k) + hcc_ad(k)
        hcc_ad(k) = 0.0_8
        CALL POPREAL8(eta(ic))
        alm_ad = alm_ad + (zol(ic)-zet(k))*eta_ad(ic)
        zol_ad(ic) = zol_ad(ic) + alm*eta_ad(ic)
        zet_ad(k) = zet_ad(k) - alm*eta_ad(ic)
        eta_ad(ic) = 0.0_8
        DO l=k,ic+1,-1
          CALL POPREAL8(eta(l))
          alm_ad = alm_ad + (zet(l)-zet(k))*eta_ad(l)
          zet_ad(l) = zet_ad(l) + alm*eta_ad(l)
          zet_ad(k) = zet_ad(k) - alm*eta_ad(l)
          eta_ad(l) = 0.0_8
        END DO
        CALL POPCONTROL1B(branch)
        IF (branch .NE. 0) alm_ad = alm_ad + 2*alm*toki_ad/lambda_min**2
      END IF
      CALL POPREAL8(alm)
      temp_ad4 = alm_ad/tem
      hol_ad(k) = hol_ad(k) + temp_ad4
      hst_ad(ic) = hst_ad(ic) - temp_ad4
      tem_ad = -((hol(k)-hst(ic))*temp_ad4/tem)
 110  DO l=k-1,ic+1,-1
        temp_ad2 = (zet(l)-zet(l+1))*tem_ad
        temp_ad3 = (hst(ic)-hol(l))*tem_ad
        hst_ad(ic) = hst_ad(ic) + temp_ad2
        hol_ad(l) = hol_ad(l) - temp_ad2
        zet_ad(l) = zet_ad(l) + temp_ad3
        zet_ad(l+1) = zet_ad(l+1) - temp_ad3
      END DO
      CALL POPREAL8(tem)
      temp_ad0 = (zol(ic)-zet(ic+1))*tem_ad
      temp_ad1 = (hst(ic)-hol(ic))*tem_ad
      hst_ad(ic) = hst_ad(ic) + temp_ad0
      hol_ad(ic) = hol_ad(ic) - temp_ad0
      zol_ad(ic) = zol_ad(ic) + temp_ad1
      zet_ad(ic+1) = zet_ad(ic+1) - temp_ad1
 120  CALL POPREAL8(lambda_min)
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        DO l=ic,ic,-1
          CALL POPREAL8(hst(l))
          hst_ad(l+1) = hst_ad(l+1) + 0.5*hst_ad(l)
          hst_ad(l) = 0.5*hst_ad(l)
        END DO
        DO l=ic+1,k-1,1
          CALL POPREAL8(hst(l))
          hst_ad(l+1) = hst_ad(l+1) + 0.25*hst_ad(l)
          hst_ad(l-1) = hst_ad(l-1) + 0.25*hst_ad(l)
          hst_ad(l) = 0.5*hst_ad(l)
        END DO
        hst_ad(ic) = hst_ad(ic) + hstx_ad
      END IF
      DO l=k,ic+1,-1
        CALL POPREAL8(qht(l))
        qol_ad(l) = qol_ad(l) + .5*qht_ad(l)
        qol_ad(l-1) = qol_ad(l-1) + .5*qht_ad(l)
        qht_ad(l) = 0.0_8
        tem = (prj(l)-prh(l-1))/(prh(l)-prh(l-1))
        CALL POPREAL8(sht(l))
        ssl_ad(l-1) = ssl_ad(l-1) + (1.0_8-tem)*sht_ad(l)
        ssl_ad(l) = ssl_ad(l) + tem*sht_ad(l)
        sht_ad(l) = 0.0_8
        CALL POPREAL8(tem)
      END DO
      poi_c_ad = 0.0_8
      qoi_c_ad = 0.0_8
      DO l=ic,k,1
        ssl_ad(l) = ssl_ad(l) + hol_ad(l) + hst_ad(l)
        CALL POPREAL8(zol(l))
        zet_ad(l+1) = zet_ad(l+1) + zet_ad(l) + zol_ad(l)
        CALL POPREAL8(zet(l))
        tem_ad = zet_ad(l)
        poi_c_ad(l) = poi_c_ad(l) + (prj(l+1)-prj(l))*cpbg*tem_ad + prj(&
&         l+1)*cp*ssl_ad(l) + (prj(l+1)-prh(l))*cpbg*zol_ad(l)
        zol_ad(l) = 0.0_8
        zet_ad(l) = 0.0_8
        CALL POPREAL8(tem)
        CALL POPREAL8(hst(l))
        qst_ad(l) = qst_ad(l) + alhl*hst_ad(l)
        hst_ad(l) = 0.0_8
        CALL POPREAL8(hol(l))
        qol_ad(l) = qol_ad(l) + alhl*hol_ad(l)
        hol_ad(l) = 0.0_8
        CALL POPREAL8(ssl(l))
        zet_ad(l+1) = zet_ad(l+1) + grav*ssl_ad(l)
        ssl_ad(l) = 0.0_8
        CALL POPCONTROL1B(branch)
        IF (branch .NE. 0) qol_ad(l) = 0.0_8
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          CALL POPREAL8(qol(l))
          qoi_c_ad(l) = qoi_c_ad(l) + qol_ad(l)
          qol_ad(l) = 0.0_8
        ELSE
          CALL POPREAL8(qol(l))
          qst_ad(l) = qst_ad(l) + rhmax*qol_ad(l)
          qol_ad(l) = 0.0_8
        END IF
      END DO
      CALL POPREAL8(sht(k+1))
      poi_c_ad(k) = poi_c_ad(k) + prj(k+1)*cp*sht_ad(k+1)
      sht_ad(k+1) = 0.0_8
      CALL POPREAL8(zet(k+1))
      zet_ad(k+1) = 0.0_8
      tpert_ad(i) = tpert_ad(i) + poi_c_ad(k)
      qoi_ad = qoi_ad + qoi_c_ad
      poi_ad = poi_ad + poi_c_ad
 130  CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        CALL POPREAL8(trg)
        temp_ad = trg_ad/((rhmx-rhmn)*qst(k))
        qoi_ad(k) = qoi_ad(k) + temp_ad
        qst_ad(k) = qst_ad(k) - qoi(k)*temp_ad/qst(k)
      ELSE
        CALL POPREAL8(trg)
      END IF
    END SUBROUTINE CLOUDE_BWD
!*********************************************************************
    SUBROUTINE CLOUDE(ic)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: ic
      REAL*8 :: deep_fact, cu_diam, wscale
!, dQx
      REAL*8 :: cli, te_a, c00_x, cli_crit_x, pete, toki, gmhx, hstx
      REAL*8 :: dt_lyr, rate, cvw_x, closs, f2, f3, f4, f5
      INTEGER :: k700
      INTRINSIC MIN
      INTRINSIC MAX
      INTRINSIC SQRT
      INTRINSIC EXP
      REAL*8 :: min2
      REAL*8 :: min1
      REAL*8 :: x5
      REAL*8 :: x4
      REAL*8 :: x3
      REAL*8 :: x2
      REAL*8 :: x1
      LOGICAL :: mask(k-ic+1)
      REAL*8 :: max2
      REAL*8 :: max1
      REAL*8 :: y1
!         !=============================AER_CLOUD local variables ====================
!         REAL(8) :: WBASE, NDROP, NICE, FP_D, FF_A, FP_I, FICE, &
!               NDROP_AMB, NSOOT_AMB, NSOOT, NIN, INSOOT, dCVW2, QICE, &
!               dQICE, dQIG, FPICE, dNICE, dNDROP, DSOOT_AMB, DSOOT, QLIQ, dQLIQ, FPRECIP, AUX, QT, &
!               MAXNICE, MAXNDROP, MINNICE, MINNDROP, NDROP_ACT, RIMM, FNDRIM, TminusTa, Tparcel , &
!	      alph_e, beta_e, RH_AMB, ECRIT
!
!         REAL(8), DIMENSION (NDUSTMAX) :: NDUST, NDUST_AMB, INDUST, DDUST_AMB, DDUST 
!         INTEGER :: INX
!
!
!         T_ICE_ALL = 238.0 !A little higher so there is ice at the freezing level
!         WBASE=1.0 
!         FICE=0.0
!         NICE=0.0
!         NDROP=0.0
!         NDROP_AMB=0.0
!         NDROP_ACT=0.0
!         NIN = 0.0
!         NDUST=0.0
!         NSOOT=0.0
!         NDUST_AMB=0.0
!         NSOOT_AMB = 0.0
!         dCVW2=0.0
!         QICE=0.0
!         QLIQ=0.0
!         FPICE = 0.0
!         INDUST=0.0
!         INSOOT= 0.0
!         QT= 0.0
!         FPRECIP=0.0   
!         FNDRIM = 0.0
!         RIMM = 0.0   
!         TminusTa = 0.0
!
!         f_seasalt = 0.0
!         aseasalt = 0.0
!
!         call init_Aer(AER_BASE)
!
!         !AER_CLOUD=============================
      alm = 0.
      IF (1. .GT. (qoi(k)/qst(k)-rhmn)/(rhmx-rhmn)) THEN
        trg = (qoi(k)/qst(k)-rhmn)/(rhmx-rhmn)
      ELSE
        trg = 1.
      END IF
      IF (0.0 .LT. (autorampb-sige(ic))/0.2) THEN
        y1 = (autorampb-sige(ic))/0.2
      ELSE
        y1 = 0.0
      END IF
      IF (1.0 .GT. y1) THEN
        f4 = y1
      ELSE
        f4 = 1.0
      END IF
! to 1 at SIG=AUTORAMPB-0.2
      IF (sige(ic) .GE. 0.5) THEN
        f5 = 1.0
      ELSE
        f5 = 1.0 - 2.*co_zdep*(0.5-sige(ic))
        IF (f5 .LT. 0.0) THEN
          f5 = 0.0
        ELSE
          f5 = f5
        END IF
      END IF
      IF (trg .LE. 1.0e-5) THEN
! TRIGGER  =========>>
!            RC(IC) = 7
        RETURN
      ELSE
!  RECOMPUTE SOUNDING UP TO DETRAINMENT LEVEL
        poi_c = poi
        qoi_c = qoi
        poi_c(k) = poi_c(k) + tpert(i)
        qoi_c(k) = qoi_c(k) + qpert(i)
        zet(k+1) = 0.
        sht(k+1) = cp*poi_c(k)*prj(k+1)
        DO l=k,ic,-1
          IF (qst(l)*rhmax .GT. qoi_c(l)) THEN
            qol(l) = qoi_c(l)
          ELSE
            qol(l) = qst(l)*rhmax
          END IF
          IF (0.000 .LT. qol(l)) THEN
            qol(l) = qol(l)
          ELSE
            qol(l) = 0.000
          END IF
          ssl(l) = cp*prj(l+1)*poi_c(l) + grav*zet(l+1)
          hol(l) = ssl(l) + qol(l)*alhl
          hst(l) = ssl(l) + qst(l)*alhl
          tem = poi_c(l)*(prj(l+1)-prj(l))*cpbg
          zet(l) = zet(l+1) + tem
          zol(l) = zet(l+1) + (prj(l+1)-prh(l))*poi_c(l)*cpbg
        END DO
        DO l=ic+1,k
          tem = (prj(l)-prh(l-1))/(prh(l)-prh(l-1))
          sht(l) = ssl(l-1) + tem*(ssl(l)-ssl(l-1))
          qht(l) = .5*(qol(l)+qol(l-1))
        END DO
! SMOOTH HSTAR W/ 1-2-1 Filter
        IF (smooth_hst) THEN
! save for later
          hstx = hst(ic)
          DO l=k-1,ic+1,-1
            hst(l) = 0.25*(hst(l+1)+hst(l-1)) + 0.5*hst(l)
          END DO
          DO l=ic,ic
            hst(l) = 0.5*hst(l+1) + 0.5*hst(l)
          END DO
        END IF
!  CALCULATE LAMBDA, ETA, AND WORKFUNCTION
        lambda_min = .2/mxdiam(i)
        lambda_max = .2/200.
!     LAMBDA_MIN = .2/(LAMBDA_FAC*DPTH_BL)
!     LAMBDA_MAX = .2/( MAX( LAMBMX_FAC*DPTH_BL , DIAMMN_MIN ) )
        IF (hol(k) .LE. hst(ic)) THEN
! CANNOT REACH IC LEVEL  ======>>
!            RC(IC) = 1
          RETURN
        ELSE
!  LAMBDA CALCULATION: MS-A18
          tem = (hst(ic)-hol(ic))*(zol(ic)-zet(ic+1))
          DO l=ic+1,k-1
            tem = tem + (hst(ic)-hol(l))*(zet(l)-zet(l+1))
          END DO
          IF (tem .LE. 0.0) THEN
! NO VALID LAMBDA  ============>>
!            RC(IC) = 2
            RETURN
          ELSE
            alm = (hol(k)-hst(ic))/tem
            IF (alm .GT. lambda_max) THEN
!            RC(IC) = 3
              RETURN
            ELSE
              toki = 1.0
              IF (alm .LT. lambda_min) toki = (alm/lambda_min)**2
!we can probably replace this by a actual distribution based on grid cell size
!RC(IC) = 6
!RETURN
!LAMBDSV(IC) = ALM
!  ETA CALCULATION: MS-A2
              DO l=ic+1,k
                eta(l) = 1.0 + alm*(zet(l)-zet(k))
              END DO
              eta(ic) = 1.0 + alm*(zol(ic)-zet(k))
!  WORKFUNCTION CALCULATION:  MS-A22
              wfn = 0.0
              hcc(k) = hol(k)
              DO l=k-1,ic+1,-1
                hcc(l) = hcc(l+1) + (eta(l)-eta(l+1))*hol(l)
                tem = hcc(l+1)*dpb(l) + hcc(l)*dpt(l)
                eht(l) = eta(l+1)*dpb(l) + eta(l)*dpt(l)
                wfn = wfn + (tem-eht(l)*hst(l))*gam(l)
              END DO
              hcc(ic) = hst(ic)*eta(ic)
              wfn = wfn + (hcc(ic+1)-hst(ic)*eta(ic+1))*gam(ic)*dpb(ic)
!  VERTICAL VELOCITY/KE CALCULATION (ADDED 12/2001 JTB)
              bk3(k) = 0.0
              bk2(k) = 0.0
              bke(k) = 0.0
              hcld(k) = hol(k)
              DO l=k-1,ic,-1
                hcld(l) = (eta(l+1)*hcld(l+1)+(eta(l)-eta(l+1))*hol(l))/&
&                 eta(l)
                tem = (hcld(l)-hst(l))*(zet(l)-zet(l+1))/(1.0+lbcp*dqq(l&
&                 ))
                IF (cldmicro .LE. 0.0) THEN
                  bke(l) = bke(l+1) + grav*tem/(cp*prj(l+1)*poi(l))
                  IF (tem .LT. 0.0) THEN
                    max1 = 0.0
                  ELSE
                    max1 = tem
                  END IF
                  bk2(l) = bk2(l+1) + grav*max1/(cp*prj(l+1)*poi(l))
                  IF (tem .GT. 0.0) THEN
                    min1 = 0.0
                  ELSE
                    min1 = tem
                  END IF
                  bk3(l) = bk3(l+1) + grav*min1/(cp*prj(l+1)*poi(l))
                  IF (bk2(l) .LT. 0.0) THEN
                    max2 = 0.0
                  ELSE
                    max2 = bk2(l)
                  END IF
                  cvw(l) = SQRT(2.0*max2)
                END IF
              END DO
! 1.0 / ( 5.0*ALM )
              cu_diam = 1000.
!   ALPHA CALCULATION 
              rasal2i = rasal2_2d(i)
              IF (zet(ic) .LT. 2000.) rasal(ic) = rasal1
              IF (zet(ic) .GE. 2000.) THEN
                IF (1.0 .GT. (zet(ic)-2000.)/8000.) THEN
                  min2 = (zet(ic)-2000.)/8000.
                ELSE
                  min2 = 1.0
                END IF
!WMP     RASAL(IC) = RASAL1 + (RASAL2i-RASAL1)*(ZET(IC) - 2000.)/8000.
                rasal(ic) = rasal1 + (rasal2i-rasal1)*min2**rasal_exp
              END IF
!WMP  RASAL(IC) = MIN( RASAL(IC) , 1.0e5 )
              rasal(ic) = dt/rasal(ic)
              mask(1:k-ic+1) = cvw(ic:k) .LT. 1.00
              WHERE (mask(1:k-ic+1)) 
                WHERE (mask(1:k-ic+1)) 
                  WHERE (mask(1:k-ic+1)) 
                    cvw(ic:k) = 1.00
                  ELSEWHERE
                    cvw(ic:k) = cvw(ic:k)
                  END WHERE
                END WHERE
              END WHERE
!  NOTE THIS "CENTRALIZES" A KLUGE PRESENT IN OTHER LOCATIONS.
!  CLEAN UP SOME TIME.      -JTB 12/04/03
!  TEST FOR CRITICAL WORK FUNCTION
              CALL ACRITN(pol(ic), prs(k), acr)
              IF (wfn .LE. acr) THEN
! SUB-CRITICAL WORK FUNCTION ======>>
!            RC(IC) = 4
                RETURN
              ELSE
!  CLOUD TOP WATER AND MOMENTUM (TIMES ETA(IC)) MS-A16
! Tracer scavenging
! RAS loops over a series of plumes all having common cloud base level K 
! and different detrainment levels IC.  The plumes operate sequentially
! on the grid box mean quantities (wind, moisture, tracer) and so each
! subsequent plume is seeing the effects of previous plumes.  We parameterize
! scavenging following Liu et al. [JGR, 2001], their equation 1:
!  AEROSOL FRACTION SCAVENGED = 1 - exp(-FSCAV*DZ)
! where FSCAV is a specified scavenging efficiency [km-1] and DZ is the
! distance [km] the tracer traverses in the plume from it's entrainment
! level to its detrainment level.  We write the aerosol fraction surviving as:
!  FNOSCAV = exp(- FSCAV_(ITR) * DZ)
! The total scavenging is proportional to the convective mass flux, which
! is not explicitly solved for at this point.
                IF (do_tracers) THEN
                  DO itr=1,itrcr
!           Scavenging of the below cloud tracer
                    delzkm = (zet(ic)-zet(k))/1000.
                    x4 = EXP(-(fscav_(itr)*delzkm))
                    IF (x4 .GT. 1.) THEN
                      x1 = 1.
                    ELSE
                      x1 = x4
                    END IF
                    IF (x1 .LT. 0.) THEN
                      fnoscav = 0.
                    ELSE
                      fnoscav = x1
                    END IF
                    xht(itr) = xoi(k, itr)*fnoscav
                  END DO
                END IF
                wlq = qol(k)
                uht = uoi(k)
                vht = voi(k)
                rnn(k) = 0.
                cll0(k) = 0.
!print *, '========================================='
                DO l=k-1,ic,-1
                  tem = eta(l) - eta(l+1)
                  wlq = wlq + tem*qol(l)
                  uht = uht + tem*uoi(l)
                  vht = vht + tem*voi(l)
                  IF (do_tracers) THEN
                    DO itr=1,itrcr
!         Scavenging of the entrained tracer.  Updates transported tracer mass.
                      delzkm = (zet(ic)-zet(l+1))/1000.
                      x5 = EXP(-(fscav_(itr)*delzkm))
                      IF (x5 .GT. 1.) THEN
                        x2 = 1.
                      ELSE
                        x2 = x5
                      END IF
                      IF (x2 .LT. 0.) THEN
                        fnoscav = 0.
                      ELSE
                        fnoscav = x2
                      END IF
                      xht(itr) = xht(itr) + tem*xoi(l, itr)*fnoscav
                    END DO
                  END IF
!!!! How much condensate (CLI) is present here? 
                  IF (l .GT. ic) THEN
                    tx2 = 0.5*(qst(l)+qst(l-1))*eta(l)
                    tx3 = 0.5*(hst(l)+hst(l-1))*eta(l)
                    qcc = tx2 + gm1(l)*(hcc(l)-tx3)
                    cll0(l) = wlq - qcc
                  ELSE
                    cll0(l) = wlq - qst(ic)*eta(ic)
                  END IF
                  IF (cll0(l) .LT. 0.00) THEN
                    cll0(l) = 0.00
                  ELSE
                    cll0(l) = cll0(l)
                  END IF
! condensate (kg/kg)
                  cli = cll0(l)/eta(l)
! Temperature (K)
                  te_a = poi(l)*prh(l)
!=====================================================================
                  IF (cldmicro .LE. 0.0) THEN
!AER_CLOUD MIcrophysics considering activation and nucleation 
!               !recompute vertical velocity
!
!               Tparcel = TE_A
!               CVW(K) = 0.8  ! Assume a below cloud base  W of 0.8 m s-1            
!               BK2(K)   = 0.0
!
!     
!               TEM     = (HCLD(L)-HST(L) )/ (1.0+LBCP*DQQ(L))  
!               TminusTa = max(min(TEM/CP, 5.0), 0.0) !limit DT to 5 K. According to Wei, JAS, 1998   
!	     TEM =0.33*TminusTa*CO_AUTO(I)/TE_A !Bouyancy term, effciency =0.5 mwr Roode et al    	     
!
!               BK2(L)  = BK2(L+1) + GRAV * TEM*(ZET(L)-ZET(L+1)) 
!               BK2(L) = BK2(L) - (ZET(L)-ZET(L+1))*(BK2(L+1)*ALM + CLI*GRAV)  !Account for drag from entrainment of stagnat air 
!and condesate loading
!               CVW(L) = max(SQRT(  2.0* MAX( BK2(L) , 0.0 )  ), 1.0) 
!
!
!	    CVW_X = MIN(CVW(L), 50.0)
!               DT_LYR  =  max(( ZET(L)-ZET(L+1) )/CVW_X, 1.0) !Sanity check 
!               TEM   = ETA(L) - ETA(L+1)
!
!               Tparcel  =  TE_A + TminusTa
!
!
!
!!!!!!!!!!account for entrainment effects on activation !!!!!!!!!!!
!               ! Barahona and Nenes, JGR, 2007
!               alph_e = 2.8915e-8*Tparcel*Tparcel -2.1328e-5*Tparcel+4.2523e-3
!               beta_e = MAPL_ALHL*TminusTa/MAPL_RVAP/Tparcel/Tparcel
!               RH_AMB=QOI(L)/QST(L)
!               ECRIT  = max(RH_AMB -beta_e, 1.0e-6) 
!               ECRIT =  alph_e/ECRIT
!               ! print *, L, Tparcel, RH_AMB, ECRIT, ALM
!	           ECRIT =  ALM/ECRIT
!       ! ECRIT = 0.0 ! do not use this for now
!               !Print *, ECRIT
!
!
!               if (L .eq. K-1) then
!
!                  FICE=0.0
!                  NICE=0.0
!                  NDROP=0.0
!                  NIN =0.0
!                  NDUST_AMB =0.0
!                  NSOOT_AMB = 0.0
!                  NSOOT=0.0
!                  NDUST= 0.0
!
!        
!                  AER_BASE =  AERO(L)
!                  RATE=0.0
!                  FPRECIP=0.0
!                  !initial conditions     
!                  call ARGact(Tparcel, CVW_X, NDROP_ACT, NDROP_AMB, NDUST_AMB, NSOOT_AMB, L,  .true., DDUST_AMB, DSOOT_AMB, ECRI
!T) !cloud droplet number and INsource at cloud base 
!                  NDUST=NDUST_AMB
!                  NSOOT=NSOOT_AMB
!                  DDUST=DDUST_AMB
!                  DSOOT=DSOOT_AMB                                     
!
!               else 
!                  call ARGact(Tparcel, CVW_X, NDROP_ACT, NDROP_AMB, NDUST_AMB, NSOOT_AMB, L, .false., DDUST_AMB, DSOOT_AMB, ECRI
!T) !cloud droplet number above cloud base  
!
!               end if
!
!               QT = CLI
!               RATE = 0.0
!               FPRECIP = 0.0
!
!               if (QT .gt. 0.0) then
!
!                  ! if FICE is already >= 1.0 then the cloud is glaciated and there is no need to do anymore partitioning
!
!                  if (FICE .ge. 1.0) then
!
!
!                     CALL  Qremoval(RATE, FICE, FP_D, FP_I, Tparcel,  & 
!                           POL(L), QT,  NICE, NDROP, CVW_X, FPICE, &
!                           DT_LYR, RIMM, CO_AUTO(I)) 
!
!
!
!                     dNICE = -NICE*FP_I 
!                     NICE  =  (NICE +dNICE*DT_LYR)*ETA(L+1)/ETA(L) !ice
!
!                     MINNICE = max(QICE*(1.0-RATE*DT_LYR), 0.0)/4.0e-8!assuming maximum vol radius 250 microns
!                     MAXNICE = max(QICE*(1.0-RATE*DT_LYR), 0.0)/2.51e-12 !assuming minimum vol radius 10 microns
!
!                     NICE=MIN(max(NICE, MINNICE), MAXNICE)
!
!                     FICE = 1.0
!
!                  else 
!
!                     ! Cloud is not completely glaciated do the whole thing
!                     ! ALL this subroutines return tendencies
!
!
!                     CALL  INfreezing(QLIQ, NDROP, NIN, NDUST, NSOOT, INDUST, INSOOT, Tparcel, POL(L), CVW_X, DDUST, DSOOT)  !ca
!lculate the freezing fraction of the aerosol at this level
!
!                     NIN = min(NIN, NDROP/DT_LYR)
!
!                     call Qgrowth(Tparcel, POL(L), QICE, NICE, QT, NIN, dQIG, RIMM, FNDRIM)
!
!                     CALL  Qremoval(RATE, FICE, FP_D, FP_I, Tparcel,  & 
!                           POL(L), QT,  NICE, NDROP, CVW_X, FPICE, &
!                           DT_LYR,  RIMM, CO_AUTO(I)) 
!
!
!
!                     !ice number tendency: -precip + freezin
!                     dNICE = -NICE*FP_I  + NIN     
!                     NICE  =  (NICE +dNICE*DT_LYR)*ETA(L+1)/ETA(L) !ice
!                     NICE =max(NICE, 0.0)
!
!
!                     !ice mass tendency: growth - precip
!                     dQICE = -QICE*FPICE + dQIG
!                     QICE  =  min((QICE + dQICE*DT_LYR)*ETA(L+1)/ETA(L), QT) !ice
!                     QICE=max(min(QT, QICE), 0.0)
!
!
!                     ! Liquid Tendency: source/evap -  precip 
!                     !dQLIQ = max((CLI-QICE), -QLIQ)/DT_LYR -QLIQ*max(RATE-FPICE, 0.0) 
!                     ! dQLIQ = CLI*(1.0-RATE*DT_LYR)/DT_LYR -dQICE - QLIQ*max(RATE-FPICE, 0.0)           
!                     !QLIQ  =  max((QLIQ + dQLIQ*DT_LYR)*ETA(L+1)/ETA(L), 0.0) !liquid. This is actually diagnostic
!                     QLIQ=max((QT-QICE), 0.0)
!
!
!                     !droplet number tendency: -precip - freezin + activation + activated entrained aerosol 
!
!
!                     dNDROP =-NDROP*FP_D - NIN -  FNDRIM*NDROP/DT_LYR + max(NDROP_ACT-NDROP, 0.0)/DT_LYR          
!
!                     !dNDROP =-NDROP*FP_D - NIN -  FNDRIM*NDROP/DT_LYR + NDROP_ACT/DT_LYR         
!
!                     NDROP =  (NDROP + dNDROP*DT_LYR)*ETA(L+1)/ETA(L) + &
!                           (ZET(L) - ZET(L+1))*ALM*MAX((NDROP_AMB-NDROP), 0.0)
!
!                     !Aerosol tendency: Entrainment - freezing 
!
!                     NDUST = (NDUST - INDUST*DT_LYR)*ETA(L+1)/ETA(L) + &
!                           (ZET(L) - ZET(L+1))*ALM*MAX(NDUST_AMB-NDUST, 0.0) 
!
!                     NSOOT =  (NSOOT - INSOOT*DT_LYR)*ETA(L+1)/ETA(L)  + &     
!                           (ZET(L) - ZET(L+1))*ALM*MAX(NSOOT_AMB-NSOOT, 0.0)  
!
!
!                           
!                     !Update FICE and perform Sanity checks
!
!
!                     MINNDROP = (1.0-FICE)*QLIQ*max(1.0-RATE*DT_LYR, 0.0)/2.e-10    !assuming maximum vol radius 36 microns
!                     MAXNDROP = (1.0-FICE)*QLIQ*max(1.0-RATE*DT_LYR, 0.0)/3.35e-14 !assuming minimum vol radius 2 microns
!                     MINNICE = QICE/4.0e-8!assuming maximum vol radius 250 microns
!                     MAXNICE = QICE/2.51e-12 !assuming minimum vol radius 10 microns
!
!                     IF ((NICE .gt. MAXNICE) .or. (NICE .lt. MINNICE))   then    
!                        !print *, 'nilim', NICE*1e-6, MINNICE*1e-6, MAXNICE*1e-6
!                     END IF
!
!                     IF ((NDROP .gt. MAXNDROP) .or. (NDROP .lt. MINNDROP))      then 
!                        !print *, 'ndroplim', NDROP*1e-6, MINNDROP*1e-6, MAXNDROP*1e-6
!                     end if
!
!
!                     NSOOT=MAX(NSOOT, 0.0)
!                     NDUST=MAX(NDUST, 0.0)              
!
!                     NDROP=MIN(max(NDROP, MINNDROP), MAXNDROP)
!                     NICE=MIN(max(NICE, MINNICE), MAXNICE)
!
!                     FICE=max(min(QICE/QT, 1.0), 0.0)
!
!                     IF (FICE .ge. 1.0) THEN !Complete glaciation 
!                        NICE=NICE+NDROP 
!                        NDROP = 0.0
!                        QICE  = QT
!                        QLIQ= 0.0
!                     END IF
!
!                     IF (Tparcel .LT. T_ICE_ALL) THEN !instantaneous freezing
!                        NICE=NICE+NDROP 
!                        NDROP = 0.0
!                        FICE  = 1.0
!                        QICE  = QT
!                        QLIQ=0.0
!                     END IF
!
!                     IF (Tparcel .GT. T_ICE_MAX) THEN !instantaneous melting
!                        NDROP=NICE+NDROP 
!                        NICE = 0.0
!                        FICE  = 0.0
!                        QICE  = 0.0
!                        QLIQ=QT
!                     END IF
!
!                  END IF
!
!               else 
!
!                  FICE =0.0 
!                  QICE = 0.0
!                  QLIQ = 0.0
!                  NICE= 0.0 
!                  NDROP = 0.0
!                  RATE =0.0
!               end if
!
!               FPRECIP= RATE*DT_LYR
!
!               !RATE=RATE*F4
!               ! NDROP=NDROP*F4
!               !NICE=NICE*(1.0-F4)
!
!               !print *, TE_A, FICE, 'NICE', NICE*1e-6, 'NDROP', NDROP*1e-6, L 
!               !print *, 'FPI', FP_I*DT_LYR, 'FPD', FP_D*DT_LYR, 'FPICE', FPICE, 'FPRE', FPRECIP, QT, QLIQ
!
!Bacmeister 2006 microphysics
                    CALL SUNDQ3_ICE(te_a, sdqv2, sdqv3, sdqvt1, f2, f3)
! * F5  ! F4 reduces AUTO for shallow clouds, F5 modifies auto for deep clouds
                    c00_x = co_auto(i)*f2*f3*f4
                    cli_crit_x = cli_crit/(f2*f3)
                    rate = c00_x*(1.0-EXP(-(cli**2/cli_crit_x**2)))
                  END IF
                  IF (cvw(l) .LT. 1.00) THEN
                    cvw_x = 1.00
                  ELSE
                    cvw_x = cvw(l)
                  END IF
! really trust it at low values
! l.h.s. DT_LYR => time in layer (L,L+1)
                  dt_lyr = (zet(l)-zet(l+1))/cvw_x
                  closs = cll0(l)*rate*dt_lyr
                  IF (closs .GT. cll0(l)) THEN
                    closs = cll0(l)
                  ELSE
                    closs = closs
                  END IF
                  cll0(l) = cll0(l) - closs
                  dll0(l) = closs
                  IF (closs .GT. 0.) THEN
                    wlq = wlq - closs
                    rnn(l) = closs
                  ELSE
                    rnn(l) = 0.
                  END IF
                END DO
!AER_CLOUD=======================================
!            CNVNDROP(IC)=NDROP
!            CNVNICE(IC)=NICE
!            CNVFICE(IC)=FICE
                wlq = wlq - qst(ic)*eta(ic)
!     CALCULATE GAMMAS AND KERNEL
! MS-A30 (W/O GRAV)
                gms(k) = (sht(k)-ssl(k))*pri(k)
! MS-A31 (W/O GRAV)
                gmh(k) = gms(k) + (qht(k)-qol(k))*pri(k)*alhl
! MS-A37 (W/O GRAV)
                akm = gmh(k)*gam(k-1)*dpb(k-1)
                tx2 = gmh(k)
                DO l=k-1,ic+1,-1
                  gms(l) = (eta(l)*(sht(l)-ssl(l))+eta(l+1)*(ssl(l)-sht(&
&                   l+1)))*pri(l)
                  gmh(l) = gms(l) + (eta(l)*(qht(l)-qol(l))+eta(l+1)*(&
&                   qol(l)-qht(l+1)))*alhl*pri(l)
                  tx2 = tx2 + (eta(l)-eta(l+1))*gmh(l)
                  akm = akm - gms(l)*eht(l)*pki(l) + tx2*ght(l)
                END DO
                gms(ic) = eta(ic+1)*(ssl(ic)-sht(ic+1))*pri(ic)
                akm = akm - gms(ic)*eta(ic+1)*dpb(ic)*pki(ic)
                gmh(ic) = gms(ic) + (eta(ic+1)*(qol(ic)-qht(ic+1))*alhl+&
&                 eta(ic)*(hst(ic)-hol(ic)))*pri(ic)
                IF (smooth_hst) gmhx = gms(ic) + (eta(ic+1)*(qol(ic)-qht&
&                   (ic+1))*alhl+eta(ic)*(hstx-hol(ic)))*pri(ic)
!    CLOUD BASE MASS FLUX
                IF (akm .GE. 0.0 .OR. wlq .LT. 0.0) THEN
!  =========>
!            RC(IC) = 5
                  RETURN
                ELSE
! MS-A39 MASS-FLUX IN Pa/step
                  wfn = -((wfn-acr)/akm)
                  x3 = rasal(ic)*trg*toki*wfn
                  IF (x3 .GT. (prs(k+1)-prs(k))*(100.*pblfrac)) THEN
                    wfn = (prs(k+1)-prs(k))*(100.*pblfrac)
                  ELSE
                    wfn = x3
                  END IF
!    CUMULATIVE PRECIP AND CLOUD-BASE MASS FLUX FOR OUTPUT
                  wfnog = wfn*gravi
                  tem = wfn*gravi
! (kg/m^2/step)
                  cll(ic) = cll(ic) + wlq*tem
! (kg/m^2/step)
                  rmf(ic) = rmf(ic) + tem
! (kg/m^2/step)
                  rmfd(ic) = rmfd(ic) + tem*eta(ic)
                  DO l=ic+1,k
! (kg/m^2/step)
                    rmfp(l) = tem*eta(l)
! (kg/m^2/step)
                    rmfc(l) = rmfc(l) + rmfp(l)
                    dllx(l) = dllx(l) + tem*dll0(l)
                    IF (cvw(l) .GT. 0.0) THEN
                      updfrp(l) = rmfp(l)*(ddt/daylen)*1000./(cvw(l)*prs&
&                       (l))
                    ELSE
                      updfrp(l) = 0.0
                    END IF
! current cloud; incloud condensate        
                    clli(l) = cll0(l)/eta(l)
!  cumulative grid mean convective condensate        
                    cllb(l) = cllb(l) + updfrp(l)*clli(l)
                    updfrc(l) = updfrc(l) + updfrp(l)
                  END DO
!    THETA AND Q CHANGE DUE TO CLOUD TYPE IC
                  DO l=ic,k
! (kg/m^2/step)
                    rns(l) = rns(l) + rnn(l)*tem
                    gmh(l) = gmh(l)*wfn
                    gms(l) = gms(l)*wfn
                    qoi(l) = qoi(l) + (gmh(l)-gms(l))*alhi
                    poi(l) = poi(l) + gms(l)*pki(l)*cpi
                    qst(l) = qst(l) + gms(l)*bet(l)*cpi
                  END DO
                  IF (smooth_hst) THEN
                    gmhx = gmhx*wfn
                    dqx = (gmhx-gmh(ic))*alhi
                    rns(ic) = rns(ic) + dqx/(pri(ic)*grav)
                  END IF
                  IF (do_tracers) THEN
!*FRICFAC*0.5
                    wfn = wfn*0.5*1.0
                    tem = wfn*pri(k)
                    DO itr=1,itrcr
                      xcu(k, itr) = xcu(k, itr) + tem*(xoi(k-1, itr)-xoi&
&                       (k, itr))
                    END DO
                    DO itr=1,itrcr
                      DO l=k-1,ic+1,-1
                        tem = wfn*pri(l)
                        xcu(l, itr) = xcu(l, itr) + tem*((xoi(l-1, itr)-&
&                         xoi(l, itr))*eta(l)+(xoi(l, itr)-xoi(l+1, itr)&
&                         )*eta(l+1))
                      END DO
                    END DO
                    tem = wfn*pri(ic)
                    DO itr=1,itrcr
                      xcu(ic, itr) = xcu(ic, itr) + (2.*(xht(itr)-xoi(ic&
&                       , itr)*(eta(ic)-eta(ic+1)))-(xoi(ic, itr)+xoi(ic&
&                       +1, itr))*eta(ic+1))*tem
                    END DO
                    DO itr=1,itrcr
                      DO l=ic,k
                        xoi(l, itr) = xoi(l, itr) + xcu(l, itr)
                      END DO
                    END DO
                  ELSE
!*FRICFAC*0.5
                    wfn = wfn*0.5*1.0
                  END IF
                  lambdsv(ic) = 1.000
!   CUMULUS FRICTION
                  IF (fricfac .LE. 0.0) THEN
!            RC(IC) = 0
!  NO CUMULUS FRICTION =========>>
                    RETURN
                  ELSE
                    wfn = wfn*fricfac*EXP(-(alm/friclambda))
                    tem = wfn*pri(k)
                    ucu(k) = ucu(k) + tem*(uoi(k-1)-uoi(k))
                    vcu(k) = vcu(k) + tem*(voi(k-1)-voi(k))
                    DO l=k-1,ic+1,-1
                      tem = wfn*pri(l)
                      ucu(l) = ucu(l) + tem*((uoi(l-1)-uoi(l))*eta(l)+(&
&                       uoi(l)-uoi(l+1))*eta(l+1))
                      vcu(l) = vcu(l) + tem*((voi(l-1)-voi(l))*eta(l)+(&
&                       voi(l)-voi(l+1))*eta(l+1))
                    END DO
                    tem = wfn*pri(ic)
                    ucu(ic) = ucu(ic) + (2.*(uht-uoi(ic)*(eta(ic)-eta(ic&
&                     +1)))-(uoi(ic)+uoi(ic+1))*eta(ic+1))*tem
                    vcu(ic) = vcu(ic) + (2.*(vht-voi(ic)*(eta(ic)-eta(ic&
&                     +1)))-(voi(ic)+voi(ic+1))*eta(ic+1))*tem
                    dissk0(ic) = eta(ic)*grav*wfnog*pri(ic)*0.5*((uht/&
&                     eta(ic)-uoi(ic))**2+(vht/eta(ic)-voi(ic))**2)
                    DO l=ic,k
                      uoi(l) = uoi(l) + ucu(l)
                      voi(l) = voi(l) + vcu(l)
                    END DO
!         RC(IC) = 0
                    RETURN
                  END IF
                END IF
              END IF
            END IF
          END IF
        END IF
      END IF
    END SUBROUTINE CLOUDE
    SUBROUTINE ACRITN(pl, plb, acr)
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: pl, plb
      REAL*8, INTENT(OUT) :: acr
      INTEGER :: iwk
!!REAL(8), PARAMETER :: FACM=0.5
      REAL*8, PARAMETER :: ph(15)=(/150.0, 200.0, 250.0, 300.0, 350.0, &
&       400.0, 450.0, 500.0, 550.0, 600.0, 650.0, 700.0, 750.0, 800.0, &
&       850.0/)
!!*FACM
      REAL*8, PARAMETER :: a(15)=(/1.6851, 1.1686, 0.7663, 0.5255, &
&       0.4100, 0.3677, 0.3151, 0.2216, 0.1521, 0.1082, 0.0750, 0.0664, &
&       0.0553, 0.0445, 0.0633/)
      INTRINSIC INT
      iwk = INT(pl*0.02 - 0.999999999)
      IF (iwk .GT. 1 .AND. iwk .LE. 15) THEN
        acr = a(iwk-1) + (pl-ph(iwk-1))*.02*(a(iwk)-a(iwk-1))
      ELSE IF (iwk .GT. 15) THEN
        acr = a(15)
      ELSE
        acr = a(1)
      END IF
      acr = acritfac*acr*(plb-pl)
      RETURN
    END SUBROUTINE ACRITN
!  Differentiation of rnevp in reverse (adjoint) mode, forward sweep (with options r8 split(GEOS_MoistGridComp.PRE_RASE GEOS_Mois
!tGridComp.PRE_PROGNO_CLOUD ras.RASE rase.CLOUDE rase.ACRITN rase.RNEVP rase.HTEST rase.FINDDTLS rase.STRAP ras.SUNDQ3_ICE cloudn
!ew.pdf_spread cloudnew.fix_up_clouds cloudnew.meltfrz cloudnew.hystpdf cloudnew.pdffrac cloudnew.pdfcondensate cloudnew.cnvsrc c
!loudnew.evap3 cloudnew.subl3 cloudnew.autocon3 cloudnew.PRECIP3 cloudnew.ICEFALL cloudnew.SETTLE_VEL cloudnew.MARSHPALMQ2 cloudn
!ew.MICRO_AA_BB_3 cloudnew.LDRADIUS3 cloudnew.ICE_FRACTION cloudnew.GET_ALHX3 cloudnew.ICEFRAC cloudnew.SUNDQ3_ICE3 qsat_util.DQS
!ATPERT)):
!   gradient     of useful results: cnv_prc3 qoi poi rns zet
!   with respect to varying inputs: cnv_prc3 qoi poi rns zet
    SUBROUTINE RNEVP_FWD()
      IMPLICIT NONE
      CALL PUSHREAL8(zet(k+1))
      zet(k+1) = 0
      DO l=k,icmin,-1
        CALL PUSHREAL8(tem)
        tem = poi(l)*(prj(l+1)-prj(l))*cpbg
        CALL PUSHREAL8(zet(l))
        zet(l) = zet(l+1) + tem
      END DO
      DO l=icmin,k
        CALL PUSHREAL8(tem)
        tem = pri(l)*grav
        cnv_prc3(i, l) = rns(l)*tem
      END DO
!! If hst is smoothed then adjusted precips may be negative
      IF (smooth_hst) THEN
        DO l=icmin,k
          IF (cnv_prc3(i, l) .LT. 0.) THEN
            CALL PUSHREAL8(qoi(l))
            qoi(l) = qoi(l) + cnv_prc3(i, l)
            CALL PUSHREAL8(poi(l))
            poi(l) = poi(l) - cnv_prc3(i, l)*(alhl/cp)/prj(l+1)
            cnv_prc3(i, l) = 0.
            CALL PUSHCONTROL1B(1)
          ELSE
            CALL PUSHCONTROL1B(0)
          END IF
        END DO
        CALL PUSHCONTROL1B(0)
      ELSE
        CALL PUSHCONTROL1B(1)
      END IF
    END SUBROUTINE RNEVP_FWD
!  Differentiation of rnevp in reverse (adjoint) mode, backward sweep (with options r8 split(GEOS_MoistGridComp.PRE_RASE GEOS_Moi
!stGridComp.PRE_PROGNO_CLOUD ras.RASE rase.CLOUDE rase.ACRITN rase.RNEVP rase.HTEST rase.FINDDTLS rase.STRAP ras.SUNDQ3_ICE cloud
!new.pdf_spread cloudnew.fix_up_clouds cloudnew.meltfrz cloudnew.hystpdf cloudnew.pdffrac cloudnew.pdfcondensate cloudnew.cnvsrc 
!cloudnew.evap3 cloudnew.subl3 cloudnew.autocon3 cloudnew.PRECIP3 cloudnew.ICEFALL cloudnew.SETTLE_VEL cloudnew.MARSHPALMQ2 cloud
!new.MICRO_AA_BB_3 cloudnew.LDRADIUS3 cloudnew.ICE_FRACTION cloudnew.GET_ALHX3 cloudnew.ICEFRAC cloudnew.SUNDQ3_ICE3 qsat_util.DQ
!SATPERT)):
!   gradient     of useful results: cnv_prc3 qoi poi rns zet
!   with respect to varying inputs: cnv_prc3 qoi poi rns zet
    SUBROUTINE RNEVP_BWD()
      IMPLICIT NONE
      INTEGER :: branch
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        DO l=k,icmin,-1
          CALL POPCONTROL1B(branch)
          IF (branch .NE. 0) THEN
            cnv_prc3_ad(i, l) = qoi_ad(l) - alhl*poi_ad(l)/(cp*prj(l+1))
            CALL POPREAL8(poi(l))
            CALL POPREAL8(qoi(l))
          END IF
        END DO
      END IF
      DO l=k,icmin,-1
        tem = pri(l)*grav
        rns_ad(l) = rns_ad(l) + tem*cnv_prc3_ad(i, l)
        cnv_prc3_ad(i, l) = 0.0_8
        CALL POPREAL8(tem)
      END DO
      DO l=icmin,k,1
        CALL POPREAL8(zet(l))
        zet_ad(l+1) = zet_ad(l+1) + zet_ad(l)
        tem_ad = zet_ad(l)
        zet_ad(l) = 0.0_8
        CALL POPREAL8(tem)
        poi_ad(l) = poi_ad(l) + (prj(l+1)-prj(l))*cpbg*tem_ad
      END DO
      CALL POPREAL8(zet(k+1))
      zet_ad(k+1) = 0.0_8
    END SUBROUTINE RNEVP_BWD
    SUBROUTINE RNEVP()
      IMPLICIT NONE
      zet(k+1) = 0
      DO l=k,icmin,-1
        tem = poi(l)*(prj(l+1)-prj(l))*cpbg
        zet(l) = zet(l+1) + tem
      END DO
      DO l=icmin,k
        tem = pri(l)*grav
        cnv_prc3(i, l) = rns(l)*tem
      END DO
!! If hst is smoothed then adjusted precips may be negative
      IF (smooth_hst) THEN
        DO l=icmin,k
          IF (cnv_prc3(i, l) .LT. 0.) THEN
            qoi(l) = qoi(l) + cnv_prc3(i, l)
            poi(l) = poi(l) - cnv_prc3(i, l)*(alhl/cp)/prj(l+1)
            cnv_prc3(i, l) = 0.
          END IF
        END DO
      END IF
      RETURN
    END SUBROUTINE RNEVP
!  Differentiation of htest in reverse (adjoint) mode, forward sweep (with options r8 split(GEOS_MoistGridComp.PRE_RASE GEOS_Mois
!tGridComp.PRE_PROGNO_CLOUD ras.RASE rase.CLOUDE rase.ACRITN rase.RNEVP rase.HTEST rase.FINDDTLS rase.STRAP ras.SUNDQ3_ICE cloudn
!ew.pdf_spread cloudnew.fix_up_clouds cloudnew.meltfrz cloudnew.hystpdf cloudnew.pdffrac cloudnew.pdfcondensate cloudnew.cnvsrc c
!loudnew.evap3 cloudnew.subl3 cloudnew.autocon3 cloudnew.PRECIP3 cloudnew.ICEFALL cloudnew.SETTLE_VEL cloudnew.MARSHPALMQ2 cloudn
!ew.MICRO_AA_BB_3 cloudnew.LDRADIUS3 cloudnew.ICE_FRACTION cloudnew.GET_ALHX3 cloudnew.ICEFRAC cloudnew.SUNDQ3_ICE3 qsat_util.DQS
!ATPERT)):
!   gradient     of useful results: hst qoi qol sht qst poi ssl
!                zet zol hol
!   with respect to varying inputs: hst qoi qol sht qst poi ssl
!                zet zol
    SUBROUTINE HTEST_FWD()
      IMPLICIT NONE
      REAL*8, DIMENSION(k0) :: hol1
      INTEGER :: lminhol
      REAL*8 :: minhol
      INTRINSIC MIN
      INTRINSIC MAX
! HOL initialized here in order not to confuse Valgrind debugger
      CALL PUSHREAL8ARRAY(hol, k0)
      hol = 0.
      CALL PUSHREAL8(zet(k+1))
      zet(k+1) = 0
      CALL PUSHREAL8(sht(k+1))
      sht(k+1) = cp*poi(k)*prj(k+1)
      DO l=k,icmin,-1
        IF (qst(l)*rhmax .GT. qoi(l)) THEN
          CALL PUSHREAL8(qol(l))
          qol(l) = qoi(l)
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHREAL8(qol(l))
          qol(l) = qst(l)*rhmax
          CALL PUSHCONTROL1B(1)
        END IF
        IF (0.000 .LT. qol(l)) THEN
          CALL PUSHCONTROL1B(0)
          qol(l) = qol(l)
        ELSE
          qol(l) = 0.000
          CALL PUSHCONTROL1B(1)
        END IF
        CALL PUSHREAL8(ssl(l))
        ssl(l) = cp*prj(l+1)*poi(l) + grav*zet(l+1)
        hol(l) = ssl(l) + qol(l)*alhl
        CALL PUSHREAL8(hst(l))
        hst(l) = ssl(l) + qst(l)*alhl
        CALL PUSHREAL8(tem)
        tem = poi(l)*(prj(l+1)-prj(l))*cpbg
        CALL PUSHREAL8(zet(l))
        zet(l) = zet(l+1) + tem
        CALL PUSHREAL8(zol(l))
        zol(l) = zet(l+1) + (prj(l+1)-prh(l))*poi(l)*cpbg
      END DO
    END SUBROUTINE HTEST_FWD
!  Differentiation of htest in reverse (adjoint) mode, backward sweep (with options r8 split(GEOS_MoistGridComp.PRE_RASE GEOS_Moi
!stGridComp.PRE_PROGNO_CLOUD ras.RASE rase.CLOUDE rase.ACRITN rase.RNEVP rase.HTEST rase.FINDDTLS rase.STRAP ras.SUNDQ3_ICE cloud
!new.pdf_spread cloudnew.fix_up_clouds cloudnew.meltfrz cloudnew.hystpdf cloudnew.pdffrac cloudnew.pdfcondensate cloudnew.cnvsrc 
!cloudnew.evap3 cloudnew.subl3 cloudnew.autocon3 cloudnew.PRECIP3 cloudnew.ICEFALL cloudnew.SETTLE_VEL cloudnew.MARSHPALMQ2 cloud
!new.MICRO_AA_BB_3 cloudnew.LDRADIUS3 cloudnew.ICE_FRACTION cloudnew.GET_ALHX3 cloudnew.ICEFRAC cloudnew.SUNDQ3_ICE3 qsat_util.DQ
!SATPERT)):
!   gradient     of useful results: hst qoi qol sht qst poi ssl
!                zet zol hol
!   with respect to varying inputs: hst qoi qol sht qst poi ssl
!                zet zol
    SUBROUTINE HTEST_BWD()
      IMPLICIT NONE
      REAL*8, DIMENSION(k0) :: hol1
      INTEGER :: lminhol
      REAL*8 :: minhol
      INTRINSIC MIN
      INTRINSIC MAX
      INTEGER :: branch
      DO l=icmin,k,1
        ssl_ad(l) = ssl_ad(l) + hol_ad(l) + hst_ad(l)
        CALL POPREAL8(zol(l))
        zet_ad(l+1) = zet_ad(l+1) + zet_ad(l) + zol_ad(l)
        CALL POPREAL8(zet(l))
        tem_ad = zet_ad(l)
        poi_ad(l) = poi_ad(l) + (prj(l+1)-prj(l))*cpbg*tem_ad + prj(l+1)&
&         *cp*ssl_ad(l) + (prj(l+1)-prh(l))*cpbg*zol_ad(l)
        zol_ad(l) = 0.0_8
        zet_ad(l) = 0.0_8
        CALL POPREAL8(tem)
        CALL POPREAL8(hst(l))
        qst_ad(l) = qst_ad(l) + alhl*hst_ad(l)
        hst_ad(l) = 0.0_8
        qol_ad(l) = qol_ad(l) + alhl*hol_ad(l)
        hol_ad(l) = 0.0_8
        CALL POPREAL8(ssl(l))
        zet_ad(l+1) = zet_ad(l+1) + grav*ssl_ad(l)
        ssl_ad(l) = 0.0_8
        CALL POPCONTROL1B(branch)
        IF (branch .NE. 0) qol_ad(l) = 0.0_8
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          CALL POPREAL8(qol(l))
          qoi_ad(l) = qoi_ad(l) + qol_ad(l)
          qol_ad(l) = 0.0_8
        ELSE
          CALL POPREAL8(qol(l))
          qst_ad(l) = qst_ad(l) + rhmax*qol_ad(l)
          qol_ad(l) = 0.0_8
        END IF
      END DO
      CALL POPREAL8(sht(k+1))
      poi_ad(k) = poi_ad(k) + prj(k+1)*cp*sht_ad(k+1)
      sht_ad(k+1) = 0.0_8
      CALL POPREAL8(zet(k+1))
      zet_ad(k+1) = 0.0_8
      CALL POPREAL8ARRAY(hol, k0)
    END SUBROUTINE HTEST_BWD
    SUBROUTINE HTEST()
      IMPLICIT NONE
      REAL*8, DIMENSION(k0) :: hol1
      INTEGER :: lminhol
      REAL*8 :: minhol
      INTRINSIC MIN
      INTRINSIC MAX
! HOL initialized here in order not to confuse Valgrind debugger
      hol = 0.
      lminhol = k + 1
      minhol = -999999.
      zet(k+1) = 0
      sht(k+1) = cp*poi(k)*prj(k+1)
      DO l=k,icmin,-1
        IF (qst(l)*rhmax .GT. qoi(l)) THEN
          qol(l) = qoi(l)
        ELSE
          qol(l) = qst(l)*rhmax
        END IF
        IF (0.000 .LT. qol(l)) THEN
          qol(l) = qol(l)
        ELSE
          qol(l) = 0.000
        END IF
        ssl(l) = cp*prj(l+1)*poi(l) + grav*zet(l+1)
        hol(l) = ssl(l) + qol(l)*alhl
        hst(l) = ssl(l) + qst(l)*alhl
        tem = poi(l)*(prj(l+1)-prj(l))*cpbg
        zet(l) = zet(l+1) + tem
        zol(l) = zet(l+1) + (prj(l+1)-prh(l))*poi(l)*cpbg
      END DO
      hol1 = hol
      DO l=k-1,icmin+1,-1
        hol1(l) = 0.25*hol(l+1) + 0.50*hol(l) + 0.25*hol(l-1)
        IF (minhol .GE. hol1(l) .OR. minhol .LT. 0.) THEN
          minhol = hol1(l)
          lminhol = l
        END IF
      END DO
      sige_minhol = sige(lminhol)
    END SUBROUTINE HTEST
    SUBROUTINE FINDDTLS()
      IMPLICIT NONE
      REAL*8 :: sigdt0, sigmax, sigmin
      INTEGER :: ll
      INTRINSIC ALLOCATED
!#ifndef __GFORTRAN__
!         integer :: THE_SEED(2)
!#else
!         integer :: THE_SEED(12)
!#endif
!         THE_SEED(1)=SEEDRAS(I,1)*IRAS(I) + SEEDRAS(I,2)*JRAS(I)
!         THE_SEED(2)=SEEDRAS(I,1)*JRAS(I) + SEEDRAS(I,2)*IRAS(I)
!         THE_SEED(1)=THE_SEED(1)*SEEDRAS(I,1)/( SEEDRAS(I,2) + 10)
!         THE_SEED(2)=THE_SEED(2)*SEEDRAS(I,1)/( SEEDRAS(I,2) + 10)
!         if(THE_SEED(1) == 0) THE_SEED(1) =  5
!         if(THE_SEED(2) == 0) THE_SEED(2) = -5
!#ifdef __GFORTRAN__
!         THE_SEED(3:12) = 0
!#endif
!
!         call random_seed(PUT=THE_SEED)
!dh RASNCL = -300
      sigmax = sige(k)
      sigmin = sige(icmin)
      IF (rasncl .LT. 0.0) n_dtl = k - icmin
!! NO SHALLOW CONV   N_DTL = 56 - ICMIN 
!            N_DTL = min( int( RASNCL ) , K-ICMIN )
      IF (ALLOCATED(icl_v)) THEN
        DEALLOCATE(icl_v)
      END IF
      ALLOCATE(icl_v(n_dtl))
!         if ( ( RASNCL < 0.0 ) .and. ( RASNCL >=-100.) ) then 
!            do L=1,N_DTL
!               ICL_V(L) = ICMIN + L - 1
!            enddo
!         else if ( RASNCL < -100.0 ) then 
      DO l=1,n_dtl
        icl_v(l) = k - l
      END DO
    END SUBROUTINE FINDDTLS
!  Differentiation of strap in reverse (adjoint) mode, forward sweep (with options r8 split(GEOS_MoistGridComp.PRE_RASE GEOS_Mois
!tGridComp.PRE_PROGNO_CLOUD ras.RASE rase.CLOUDE rase.ACRITN rase.RNEVP rase.HTEST rase.FINDDTLS rase.STRAP ras.SUNDQ3_ICE cloudn
!ew.pdf_spread cloudnew.fix_up_clouds cloudnew.meltfrz cloudnew.hystpdf cloudnew.pdffrac cloudnew.pdfcondensate cloudnew.cnvsrc c
!loudnew.evap3 cloudnew.subl3 cloudnew.autocon3 cloudnew.PRECIP3 cloudnew.ICEFALL cloudnew.SETTLE_VEL cloudnew.MARSHPALMQ2 cloudn
!ew.MICRO_AA_BB_3 cloudnew.LDRADIUS3 cloudnew.ICE_FRACTION cloudnew.GET_ALHX3 cloudnew.ICEFRAC cloudnew.SUNDQ3_ICE3 qsat_util.DQS
!ATPERT)):
!   gradient     of useful results: clw tho dqs qho vho cnv_updfrc
!                uho qss flxd gm1 cvw dqq ght bet qoi_sv uoi_sv
!                qoi voi qst updfrc poi rns poi_sv uoi rmfd voi_sv
!                updfrp cll0 tmpt cll rmfp gam
!   with respect to varying inputs: clw tho dqs qho vho cnv_updfrc
!                uho qss flxd gm1 cvw dqq ght bet qoi_sv uoi_sv
!                qoi voi qst updfrc poi rns poi_sv uoi rmfd voi_sv
!                updfrp cll0 tmpt cll rmfp gam
    SUBROUTINE STRAP_FWD(final)
      IMPLICIT NONE
!            CNV_CVW   (I,ICMIN:K-1)   =     CVW(ICMIN:K-1)
!            CNV_QC(I,ICMIN:K-1)       =  CLLB(ICMIN:K-1)
      INTEGER :: final
      REAL*8, DIMENSION(k0) :: wght, massf
      REAL*8 :: wght0, prcbl
      INTEGER, PARAMETER :: nrands=1
      REAL*8 :: rndu(nrands)
      INTEGER :: seedcbl(nrands), jj
! !DESCRIPTION: 
!   {\tt STRAP} is called: FINAL=0, to compute cloud base layer CBL properties
!   given a value K for the index of the upper {\em EDGE} of the CBL; FINAL=1
!   to redistribute convective tendencies within CBL
      INTEGER :: kk
      INTRINSIC SQRT
      INTRINSIC MAX
      INTRINSIC ABS
      INTRINSIC PRESENT
      INTRINSIC ALLOCATED
      REAL*8 :: abs1
      REAL*8 :: abs0
!  LOCAL VARIABLES FOR USE IN CLOUDE
!!IF (.NOT. PRESENT(FINAL)) THEN
      IF (final .EQ. 0) THEN
!!PRJ(ICMIN:K+1) = PKE(I,ICMIN:K+1)
        DO kk=icmin,k+1
          CALL PUSHREAL8(prj(kk))
          prj(kk) = pke(i, kk)
        END DO
! These initialized here in order not to confuse Valgrind debugger
        CALL PUSHREAL8ARRAY(poi, k0)
        poi = 0.
! Do not believe it actually makes any difference.
        CALL PUSHREAL8ARRAY(qoi, k0)
        qoi = 0.
        CALL PUSHREAL8ARRAY(uoi, k0)
        uoi = 0.
        CALL PUSHREAL8ARRAY(voi, k0)
        voi = 0.
        CALL PUSHREAL8ARRAY(prs(icmin:k0+1), k0 - icmin + 2)
        prs(icmin:k0+1) = ple(i, icmin:k0+1)
        poi(icmin:k) = tho(i, icmin:k)
        qoi(icmin:k) = qho(i, icmin:k)
        uoi(icmin:k) = uho(i, icmin:k)
        voi(icmin:k) = vho(i, icmin:k)
        CALL PUSHREAL8ARRAY(qst(icmin:k), k - icmin + 1)
        qst(icmin:k) = qss(i, icmin:k)
        CALL PUSHREAL8ARRAY(dqq(icmin:k), k - icmin + 1)
        dqq(icmin:k) = dqs(i, icmin:k)
        IF (do_tracers) THEN
          DO itr=1,itrcr
            xoi(icmin:k, itr) = xho(i, icmin:k, itr)
          END DO
        END IF
!!! Mass fraction of each layer below cloud base
!!! contributed to aggregate cloudbase layer (CBL) 
        massf(:) = wgt0(i, :)
!!! RESET PRESSURE at bottom edge of CBL 
        prcbl = prs(k)
        DO l=k,k0
          prcbl = prcbl + massf(l)*(prs(l+1)-prs(l))
        END DO
        CALL PUSHREAL8(prs(k+1))
        prs(k+1) = prcbl
        CALL PUSHREAL8(prj(k+1))
        prj(k+1) = (prs(k+1)/1000.)**(mapl_rgas/mapl_cp)
        DO l=k,icmin,-1
          pol(l) = 0.5*(prs(l)+prs(l+1))
          CALL PUSHREAL8(prh(l))
          prh(l) = (prs(l+1)*prj(l+1)-prs(l)*prj(l))/(onepkap*(prs(l+1)-&
&           prs(l)))
          CALL PUSHREAL8(pki(l))
          pki(l) = 1.0/prh(l)
          CALL PUSHREAL8(dpt(l))
          dpt(l) = prh(l) - prj(l)
          CALL PUSHREAL8(dpb(l))
          dpb(l) = prj(l+1) - prh(l)
          CALL PUSHREAL8(pri(l))
          pri(l) = .01/(prs(l+1)-prs(l))
        END DO
!!!!! RECALCULATE PROFILE QUAN. IN LOWEST STRAPPED LAYER
        IF (k .LE. k0) THEN
          poi(k) = 0.
          qoi(k) = 0.
          uoi(k) = 0.
          voi(k) = 0.
!! SPECIFY WEIGHTS GIVEN TO EACH LAYER WITHIN SUBCLOUD "SUPERLAYER"
          wght = 0.
          DO l=k,k0
            wght(l) = massf(l)*(ple(i, l+1)-ple(i, l))/(prs(k+1)-prs(k))
          END DO
          DO l=k,k0
            poi(k) = poi(k) + wght(l)*tho(i, l)
            qoi(k) = qoi(k) + wght(l)*qho(i, l)
            uoi(k) = uoi(k) + wght(l)*uho(i, l)
            voi(k) = voi(k) + wght(l)*vho(i, l)
          END DO
          IF (do_tracers) THEN
            xoi(k, :) = 0.
            DO itr=1,itrcr
              DO l=k,k0
                xoi(k, itr) = xoi(k, itr) + wght(l)*xho(i, l, itr)
              END DO
            END DO
          END IF
          tmpt(k) = poi(k)*prh(k)
          CALL DQSATPERT_FWD(dqq(k), qst(k), tmpt(k), pol(k), 1)
          CALL PUSHCONTROL1B(1)
        ELSE
          CALL PUSHCONTROL1B(0)
        END IF
!!DPTH_BL = CPBG*POI(K)*( PRJ(K+1)-PRJ(K) )
! seedras(1,2) are both integers passed from
! from GEOS_Moist w/ values 0 - 1000000
! rndu(:) = 1.0*( seedras(1)+seedras(2) )/2000000.
!dh tap produces code with error
        DO jj=1,nrands
          IF (seedras(i, 1)/1000000. .LT. 1e-6) THEN
            CALL PUSHCONTROL1B(0)
            rndu(jj) = 1e-6
          ELSE
            CALL PUSHCONTROL1B(1)
            rndu(jj) = seedras(i, 1)/1000000.
          END IF
        END DO
        IF (maxdallowed_d .GE. 0.) THEN
          abs0 = maxdallowed_d
        ELSE
          abs0 = -maxdallowed_d
        END IF
        IF (maxdallowed_s .GE. 0.) THEN
          abs1 = maxdallowed_s
        ELSE
          abs1 = -maxdallowed_s
        END IF
!!call congvec( npoints , seedcbl , rndu )
!            DPTH_BL   = ZCBL(I)
        mxdiam(i) = cnv_fraction(i)*abs0 + (1-cnv_fraction(i))*abs1
        IF (maxdallowed_d .GT. 0) mxdiam(i) = mxdiam(i)*rndu(1)**(-(1./&
&           2.))
! Make MXDIAM stochastic
        DO l=k,icmin,-1
!*
          CALL PUSHREAL8(bet(l))
          bet(l) = dqq(l)*pki(l)
!*
          CALL PUSHREAL8(gam(l))
          gam(l) = pki(l)/(1.0+lbcp*dqq(l))
          IF (l .LT. k) THEN
            CALL PUSHREAL8(ght(l+1))
            ght(l+1) = gam(l)*dpb(l) + gam(l+1)*dpt(l+1)
            CALL PUSHREAL8(gm1(l+1))
            gm1(l+1) = 0.5*lbcp*(dqq(l)/(alhl*(1.0+lbcp*dqq(l)))+dqq(l+1&
&             )/(alhl*(1.0+lbcp*dqq(l+1))))
            CALL PUSHCONTROL1B(1)
          ELSE
            CALL PUSHCONTROL1B(0)
          END IF
        END DO
        rns = 0.
        cll = 0.
        rmf = 0.
        rmfd = 0.
        rmfc = 0.
        CALL PUSHREAL8ARRAY(rmfp, k0)
        rmfp = 0.
        cll0 = 0.
        dll0 = 0.
        dllx = 0.
        clli = 0.
        cllb = 0.
        poi_sv = poi
        qoi_sv = qoi
        uoi_sv = uoi
        voi_sv = voi
        IF (do_tracers) xoi_sv = xoi
        CALL PUSHREAL8ARRAY(cvw, k0)
        cvw = 0.0
        updfrc = 0.0
        updfrp = 0.0
        dissk0 = 0.0
        CALL PUSHCONTROL1B(0)
      ELSE
        CALL PUSHCONTROL1B(1)
      END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!    IF (PRESENT(FINAL)) THEN
      IF (final .EQ. 1) THEN
        CALL PUSHREAL8ARRAY(tho(i, icmin:k-1), k - icmin)
        tho(i, icmin:k-1) = poi(icmin:k-1)
        CALL PUSHREAL8ARRAY(qho(i, icmin:k-1), k - icmin)
        qho(i, icmin:k-1) = qoi(icmin:k-1)
        uho(i, icmin:k-1) = uoi(icmin:k-1)
        vho(i, icmin:k-1) = voi(icmin:k-1)
        cnv_updfrc(i, icmin:k-1) = updfrc(icmin:k-1)
!======================AER_CLOUD=============
!               CNV_NDROP   (I,ICMIN:K-1)  =    CNVNDROP(ICMIN:K-1) !DONIF
!               CNV_NICE   (I,ICMIN:K-1)   =     CNVNICE(ICMIN:K-1) !DONIF
!               CNV_FICE   (I,ICMIN:K-1)   =     CNVFICE(ICMIN:K-1) !DONIF
!! De-strap tendencies from RAS
!! specify weighting "SHAPE"
        CALL PUSHREAL8ARRAY(wght, k0)
        wght = wgt1(i, :)
!! Scale properly by layer masses
        wght0 = 0.
        DO l=k,k0
          wght0 = wght0 + wght(l)*(ple(i, l+1)-ple(i, l))
        END DO
        wght0 = (prs(k+1)-prs(k))/wght0
        wght = wght0*wght
        DO l=k,k0
          CALL PUSHREAL8(tho(i, l))
          tho(i, l) = tho(i, l) + wght(l)*(poi(k)-poi_sv(k))
          CALL PUSHREAL8(qho(i, l))
          qho(i, l) = qho(i, l) + wght(l)*(qoi(k)-qoi_sv(k))
          uho(i, l) = uho(i, l) + wght(l)*(uoi(k)-uoi_sv(k))
          vho(i, l) = vho(i, l) + wght(l)*(voi(k)-voi_sv(k))
        END DO
        IF (do_tracers) THEN
          xho(i, icmin:k-1, :) = xoi(icmin:k-1, :)
          DO itr=1,itrcr
            DO l=k,k0
              xho(i, l, itr) = xho(i, l, itr) + wght(l)*(xoi(k, itr)-&
&               xoi_sv(k, itr))
            END DO
          END DO
        END IF
!            FLX (I,ICMIN:K) = RMF (ICMIN:K) * DDT/DAYLEN  !  (KG/m^2/s @ CLOUD BASE)
!  (KG/m^2/s @ CLOUD TOP)
        flxd(i, icmin:k) = rmfd(icmin:k)*ddt/daylen
!            FLXC(I,ICMIN:K) = RMFC(ICMIN:K) * DDT/DAYLEN  !  (KG/m^2/s @ CLOUD TOP)
!  (KG/m^2/s )
        clw(i, icmin:k) = cll(icmin:k)*ddt/daylen
!            FLX (I,1:ICMIN-1) = 0.
        flxd(i, 1:icmin-1) = 0.
!            FLXC(I,1:ICMIN-1) = 0.
        clw(i, 1:icmin-1) = 0.
        IF (k .LT. k0) THEN
!               FLX (I,K:K0) = 0.
          flxd(i, k:k0) = 0.
!               FLXC(I,K:K0) = 0.
          clw(i, k:k0) = 0.
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(1)
        END IF
!            IRC (I,ICMIN:K-1) = RC(ICMIN:K-1)
        IF (ALLOCATED(icl_v)) THEN
          CALL PUSHCONTROL1B(0)
        ELSE
          CALL PUSHCONTROL1B(0)
        END IF
      ELSE
        CALL PUSHCONTROL1B(1)
      END IF
      IF (final .EQ. 2) THEN
!            FLX (I,:) = 0.
        flxd(i, :) = 0.
!            FLXC(I,:) = 0.
        clw(i, :) = 0.
!            IRC (I,ICMIN:K-1) = RC(ICMIN:K-1)
        CALL PUSHREAL8ARRAY(wght, k0)
        CALL PUSHCONTROL1B(0)
      ELSE
        CALL PUSHREAL8ARRAY(wght, k0)
        CALL PUSHCONTROL1B(1)
      END IF
    END SUBROUTINE STRAP_FWD
!  Differentiation of strap in reverse (adjoint) mode, backward sweep (with options r8 split(GEOS_MoistGridComp.PRE_RASE GEOS_Moi
!stGridComp.PRE_PROGNO_CLOUD ras.RASE rase.CLOUDE rase.ACRITN rase.RNEVP rase.HTEST rase.FINDDTLS rase.STRAP ras.SUNDQ3_ICE cloud
!new.pdf_spread cloudnew.fix_up_clouds cloudnew.meltfrz cloudnew.hystpdf cloudnew.pdffrac cloudnew.pdfcondensate cloudnew.cnvsrc 
!cloudnew.evap3 cloudnew.subl3 cloudnew.autocon3 cloudnew.PRECIP3 cloudnew.ICEFALL cloudnew.SETTLE_VEL cloudnew.MARSHPALMQ2 cloud
!new.MICRO_AA_BB_3 cloudnew.LDRADIUS3 cloudnew.ICE_FRACTION cloudnew.GET_ALHX3 cloudnew.ICEFRAC cloudnew.SUNDQ3_ICE3 qsat_util.DQ
!SATPERT)):
!   gradient     of useful results: clw tho dqs qho vho cnv_updfrc
!                uho qss flxd gm1 cvw dqq ght bet qoi_sv uoi_sv
!                qoi voi qst updfrc poi rns poi_sv uoi rmfd voi_sv
!                updfrp cll0 tmpt cll rmfp gam
!   with respect to varying inputs: clw tho dqs qho vho cnv_updfrc
!                uho qss flxd gm1 cvw dqq ght bet qoi_sv uoi_sv
!                qoi voi qst updfrc poi rns poi_sv uoi rmfd voi_sv
!                updfrp cll0 tmpt cll rmfp gam
    SUBROUTINE STRAP_BWD(final)
      IMPLICIT NONE
!            CNV_CVW   (I,ICMIN:K-1)   =     CVW(ICMIN:K-1)
!            CNV_QC(I,ICMIN:K-1)       =  CLLB(ICMIN:K-1)
      INTEGER :: final
      REAL*8, DIMENSION(k0) :: wght, massf
      REAL*8 :: wght0, prcbl
      INTEGER, PARAMETER :: nrands=1
      REAL*8 :: rndu(nrands)
      INTEGER :: seedcbl(nrands), jj
      INTEGER :: kk
      INTRINSIC SQRT
      INTRINSIC MAX
      INTRINSIC ABS
      INTRINSIC PRESENT
      INTRINSIC ALLOCATED
      INTEGER :: branch
      REAL*8 :: temp1
      REAL*8 :: temp0
      REAL*8 :: temp_ad1
      REAL*8 :: temp_ad0
      REAL*8 :: temp_ad
      REAL*8 :: abs1
      REAL*8 :: abs0
      REAL*8 :: temp
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        CALL POPREAL8ARRAY(wght, k0)
        clw_ad(i, :) = 0.0_8
        flxd_ad(i, :) = 0.0_8
      ELSE
        CALL POPREAL8ARRAY(wght, k0)
      END IF
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        CALL POPCONTROL1B(branch)
        IF (branch .EQ. 0) THEN
          clw_ad(i, k:k0) = 0.0_8
          flxd_ad(i, k:k0) = 0.0_8
        END IF
        clw_ad(i, 1:icmin-1) = 0.0_8
        flxd_ad(i, 1:icmin-1) = 0.0_8
        cll_ad(icmin:k) = cll_ad(icmin:k) + ddt*clw_ad(i, icmin:k)/&
&         daylen
        clw_ad(i, icmin:k) = 0.0_8
        rmfd_ad(icmin:k) = rmfd_ad(icmin:k) + ddt*flxd_ad(i, icmin:k)/&
&         daylen
        flxd_ad(i, icmin:k) = 0.0_8
        DO l=k0,k,-1
          voi_ad(k) = voi_ad(k) + wght(l)*vho_ad(i, l)
          voi_sv_ad(k) = voi_sv_ad(k) - wght(l)*vho_ad(i, l)
          uoi_ad(k) = uoi_ad(k) + wght(l)*uho_ad(i, l)
          uoi_sv_ad(k) = uoi_sv_ad(k) - wght(l)*uho_ad(i, l)
          CALL POPREAL8(qho(i, l))
          qoi_ad(k) = qoi_ad(k) + wght(l)*qho_ad(i, l)
          qoi_sv_ad(k) = qoi_sv_ad(k) - wght(l)*qho_ad(i, l)
          CALL POPREAL8(tho(i, l))
          poi_ad(k) = poi_ad(k) + wght(l)*tho_ad(i, l)
          poi_sv_ad(k) = poi_sv_ad(k) - wght(l)*tho_ad(i, l)
        END DO
        CALL POPREAL8ARRAY(wght, k0)
        updfrc_ad(icmin:k-1) = updfrc_ad(icmin:k-1) + cnv_updfrc_ad(i, &
&         icmin:k-1)
        cnv_updfrc_ad(i, icmin:k-1) = 0.0_8
        voi_ad(icmin:k-1) = voi_ad(icmin:k-1) + vho_ad(i, icmin:k-1)
        vho_ad(i, icmin:k-1) = 0.0_8
        uoi_ad(icmin:k-1) = uoi_ad(icmin:k-1) + uho_ad(i, icmin:k-1)
        uho_ad(i, icmin:k-1) = 0.0_8
        CALL POPREAL8ARRAY(qho(i, icmin:k-1), k - icmin)
        qoi_ad(icmin:k-1) = qoi_ad(icmin:k-1) + qho_ad(i, icmin:k-1)
        qho_ad(i, icmin:k-1) = 0.0_8
        CALL POPREAL8ARRAY(tho(i, icmin:k-1), k - icmin)
        poi_ad(icmin:k-1) = poi_ad(icmin:k-1) + tho_ad(i, icmin:k-1)
        tho_ad(i, icmin:k-1) = 0.0_8
      END IF
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        CALL POPREAL8ARRAY(cvw, k0)
        voi_ad = voi_ad + voi_sv_ad
        uoi_ad = uoi_ad + uoi_sv_ad
        qoi_ad = qoi_ad + qoi_sv_ad
        poi_ad = poi_ad + poi_sv_ad
        CALL POPREAL8ARRAY(rmfp, k0)
        DO l=icmin,k,1
          CALL POPCONTROL1B(branch)
          IF (branch .NE. 0) THEN
            CALL POPREAL8(gm1(l+1))
            temp1 = alhl*(lbcp*dqq(l+1)+1.0)
            temp0 = alhl*(lbcp*dqq(l)+1.0)
            temp_ad = lbcp*0.5*gm1_ad(l+1)
            temp_ad0 = temp_ad/temp0
            temp_ad1 = temp_ad/temp1
            dqq_ad(l) = dqq_ad(l) + (1.0_8-alhl*dqq(l)*lbcp/temp0)*&
&             temp_ad0
            dqq_ad(l+1) = dqq_ad(l+1) + (1.0_8-alhl*dqq(l+1)*lbcp/temp1)&
&             *temp_ad1
            gm1_ad(l+1) = 0.0_8
            CALL POPREAL8(ght(l+1))
            gam_ad(l) = gam_ad(l) + dpb(l)*ght_ad(l+1)
            gam_ad(l+1) = gam_ad(l+1) + dpt(l+1)*ght_ad(l+1)
            ght_ad(l+1) = 0.0_8
          END IF
          CALL POPREAL8(gam(l))
          temp = lbcp*dqq(l) + 1.0
          dqq_ad(l) = dqq_ad(l) + pki(l)*bet_ad(l) - pki(l)*lbcp*gam_ad(&
&           l)/temp**2
          gam_ad(l) = 0.0_8
          CALL POPREAL8(bet(l))
          bet_ad(l) = 0.0_8
        END DO
        DO jj=nrands,1,-1
          CALL POPCONTROL1B(branch)
        END DO
        CALL POPCONTROL1B(branch)
        IF (branch .NE. 0) THEN
          CALL DQSATPERT_BWD(dqq(k), dqq_ad(k), qst(k), qst_ad(k), tmpt(&
&                      k), tmpt_ad(k), pol(k), 1)
          poi_ad(k) = poi_ad(k) + prh(k)*tmpt_ad(k)
          tmpt_ad(k) = 0.0_8
          DO l=k0,k,-1
            vho_ad(i, l) = vho_ad(i, l) + wght(l)*voi_ad(k)
            uho_ad(i, l) = uho_ad(i, l) + wght(l)*uoi_ad(k)
            qho_ad(i, l) = qho_ad(i, l) + wght(l)*qoi_ad(k)
            tho_ad(i, l) = tho_ad(i, l) + wght(l)*poi_ad(k)
          END DO
          voi_ad(k) = 0.0_8
          uoi_ad(k) = 0.0_8
          qoi_ad(k) = 0.0_8
          poi_ad(k) = 0.0_8
        END IF
        DO l=icmin,k,1
          CALL POPREAL8(pri(l))
          CALL POPREAL8(dpb(l))
          CALL POPREAL8(dpt(l))
          CALL POPREAL8(pki(l))
          CALL POPREAL8(prh(l))
        END DO
        CALL POPREAL8(prj(k+1))
        CALL POPREAL8(prs(k+1))
        CALL POPREAL8ARRAY(dqq(icmin:k), k - icmin + 1)
        dqs_ad(i, icmin:k) = dqs_ad(i, icmin:k) + dqq_ad(icmin:k)
        dqq_ad(icmin:k) = 0.0_8
        CALL POPREAL8ARRAY(qst(icmin:k), k - icmin + 1)
        qss_ad(i, icmin:k) = qss_ad(i, icmin:k) + qst_ad(icmin:k)
        qst_ad(icmin:k) = 0.0_8
        vho_ad(i, icmin:k) = vho_ad(i, icmin:k) + voi_ad(icmin:k)
        uho_ad(i, icmin:k) = uho_ad(i, icmin:k) + uoi_ad(icmin:k)
        qho_ad(i, icmin:k) = qho_ad(i, icmin:k) + qoi_ad(icmin:k)
        tho_ad(i, icmin:k) = tho_ad(i, icmin:k) + poi_ad(icmin:k)
        CALL POPREAL8ARRAY(prs(icmin:k0+1), k0 - icmin + 2)
        CALL POPREAL8ARRAY(voi, k0)
        CALL POPREAL8ARRAY(uoi, k0)
        CALL POPREAL8ARRAY(qoi, k0)
        CALL POPREAL8ARRAY(poi, k0)
        DO kk=k+1,icmin,-1
          CALL POPREAL8(prj(kk))
        END DO
        cvw_ad = 0.0_8
        qoi_sv_ad = 0.0_8
        uoi_sv_ad = 0.0_8
        qoi_ad = 0.0_8
        voi_ad = 0.0_8
        updfrc_ad = 0.0_8
        poi_ad = 0.0_8
        rns_ad = 0.0_8
        poi_sv_ad = 0.0_8
        uoi_ad = 0.0_8
        rmfd_ad = 0.0_8
        voi_sv_ad = 0.0_8
        updfrp_ad = 0.0_8
        cll0_ad = 0.0_8
        cll_ad = 0.0_8
        rmfp_ad = 0.0_8
      END IF
    END SUBROUTINE STRAP_BWD
    SUBROUTINE STRAP(final)
      IMPLICIT NONE
      INTEGER :: final
      REAL*8, DIMENSION(k0) :: wght, massf
      REAL*8 :: wght0, prcbl
      INTEGER, PARAMETER :: nrands=1
      REAL*8 :: rndu(nrands)
      INTEGER :: seedcbl(nrands), jj
! !DESCRIPTION: 
!   {\tt STRAP} is called: FINAL=0, to compute cloud base layer CBL properties
!   given a value K for the index of the upper {\em EDGE} of the CBL; FINAL=1
!   to redistribute convective tendencies within CBL
      INTEGER :: kk
      INTRINSIC SQRT
      INTRINSIC MAX
      INTRINSIC ABS
      INTRINSIC PRESENT
      INTRINSIC ALLOCATED
      REAL*8 :: abs1
      REAL*8 :: abs0
!  LOCAL VARIABLES FOR USE IN CLOUDE
!!IF (.NOT. PRESENT(FINAL)) THEN
      IF (final .EQ. 0) THEN
!!PRJ(ICMIN:K+1) = PKE(I,ICMIN:K+1)
        DO kk=icmin,k+1
          prj(kk) = pke(i, kk)
        END DO
! These initialized here in order not to confuse Valgrind debugger
        poi = 0.
! Do not believe it actually makes any difference.
        qoi = 0.
        uoi = 0.
        voi = 0.
        prs(icmin:k0+1) = ple(i, icmin:k0+1)
        poi(icmin:k) = tho(i, icmin:k)
        qoi(icmin:k) = qho(i, icmin:k)
        uoi(icmin:k) = uho(i, icmin:k)
        voi(icmin:k) = vho(i, icmin:k)
        wsp(icmin:k) = SQRT((uoi(icmin:k)-uoi(k))**2 + (voi(icmin:k)-voi&
&         (k))**2)
        qst(icmin:k) = qss(i, icmin:k)
        dqq(icmin:k) = dqs(i, icmin:k)
        IF (do_tracers) THEN
          DO itr=1,itrcr
            xoi(icmin:k, itr) = xho(i, icmin:k, itr)
          END DO
        END IF
!!! Mass fraction of each layer below cloud base
!!! contributed to aggregate cloudbase layer (CBL) 
        massf(:) = wgt0(i, :)
!!! RESET PRESSURE at bottom edge of CBL 
        prcbl = prs(k)
        DO l=k,k0
          prcbl = prcbl + massf(l)*(prs(l+1)-prs(l))
        END DO
        prs(k+1) = prcbl
        prj(k+1) = (prs(k+1)/1000.)**(mapl_rgas/mapl_cp)
        DO l=k,icmin,-1
          pol(l) = 0.5*(prs(l)+prs(l+1))
          prh(l) = (prs(l+1)*prj(l+1)-prs(l)*prj(l))/(onepkap*(prs(l+1)-&
&           prs(l)))
          pki(l) = 1.0/prh(l)
          dpt(l) = prh(l) - prj(l)
          dpb(l) = prj(l+1) - prh(l)
          pri(l) = .01/(prs(l+1)-prs(l))
        END DO
!!!!! RECALCULATE PROFILE QUAN. IN LOWEST STRAPPED LAYER
        IF (k .LE. k0) THEN
          poi(k) = 0.
          qoi(k) = 0.
          uoi(k) = 0.
          voi(k) = 0.
!! SPECIFY WEIGHTS GIVEN TO EACH LAYER WITHIN SUBCLOUD "SUPERLAYER"
          wght = 0.
          DO l=k,k0
            wght(l) = massf(l)*(ple(i, l+1)-ple(i, l))/(prs(k+1)-prs(k))
          END DO
          DO l=k,k0
            poi(k) = poi(k) + wght(l)*tho(i, l)
            qoi(k) = qoi(k) + wght(l)*qho(i, l)
            uoi(k) = uoi(k) + wght(l)*uho(i, l)
            voi(k) = voi(k) + wght(l)*vho(i, l)
          END DO
          IF (do_tracers) THEN
            xoi(k, :) = 0.
            DO itr=1,itrcr
              DO l=k,k0
                xoi(k, itr) = xoi(k, itr) + wght(l)*xho(i, l, itr)
              END DO
            END DO
          END IF
          tmpt(k) = poi(k)*prh(k)
          CALL DQSATPERT(dqq(k), qst(k), tmpt(k), pol(k), 1)
        END IF
!!DPTH_BL = CPBG*POI(K)*( PRJ(K+1)-PRJ(K) )
! seedras(1,2) are both integers passed from
! from GEOS_Moist w/ values 0 - 1000000
! rndu(:) = 1.0*( seedras(1)+seedras(2) )/2000000.
!dh tap produces code with error
        DO jj=1,nrands
          IF (seedras(i, 1)/1000000. .LT. 1e-6) THEN
            rndu(jj) = 1e-6
          ELSE
            rndu(jj) = seedras(i, 1)/1000000.
          END IF
        END DO
        IF (maxdallowed_d .GE. 0.) THEN
          abs0 = maxdallowed_d
        ELSE
          abs0 = -maxdallowed_d
        END IF
        IF (maxdallowed_s .GE. 0.) THEN
          abs1 = maxdallowed_s
        ELSE
          abs1 = -maxdallowed_s
        END IF
!!call congvec( npoints , seedcbl , rndu )
!            DPTH_BL   = ZCBL(I)
        mxdiam(i) = cnv_fraction(i)*abs0 + (1-cnv_fraction(i))*abs1
        IF (maxdallowed_d .GT. 0) mxdiam(i) = mxdiam(i)*rndu(1)**(-(1./&
&           2.))
! Make MXDIAM stochastic
        DO l=k,icmin,-1
!*
          bet(l) = dqq(l)*pki(l)
!*
          gam(l) = pki(l)/(1.0+lbcp*dqq(l))
          IF (l .LT. k) THEN
            ght(l+1) = gam(l)*dpb(l) + gam(l+1)*dpt(l+1)
            gm1(l+1) = 0.5*lbcp*(dqq(l)/(alhl*(1.0+lbcp*dqq(l)))+dqq(l+1&
&             )/(alhl*(1.0+lbcp*dqq(l+1))))
          END IF
        END DO
        tcu(icmin:k) = -(poi(icmin:k)*prh(icmin:k))
        qcu(icmin:k) = -qoi(icmin:k)
        rns = 0.
        cll = 0.
        rmf = 0.
        rmfd = 0.
        rmfc = 0.
        rmfp = 0.
        cll0 = 0.
        dll0 = 0.
        cllx = 0.
        dllx = 0.
        clli = 0.
        cllb = 0.
        poi_sv = poi
        qoi_sv = qoi
        uoi_sv = uoi
        voi_sv = voi
        IF (do_tracers) xoi_sv = xoi
        lambdsv = 0.0
        cvw = 0.0
        updfrc = 0.0
        updfrp = 0.0
        dissk0 = 0.0
      END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!    IF (PRESENT(FINAL)) THEN
      IF (final .EQ. 1) THEN
        tho(i, icmin:k-1) = poi(icmin:k-1)
        qho(i, icmin:k-1) = qoi(icmin:k-1)
        uho(i, icmin:k-1) = uoi(icmin:k-1)
        vho(i, icmin:k-1) = voi(icmin:k-1)
        cnv_updfrc(i, icmin:k-1) = updfrc(icmin:k-1)
!======================AER_CLOUD=============
!               CNV_NDROP   (I,ICMIN:K-1)  =    CNVNDROP(ICMIN:K-1) !DONIF
!               CNV_NICE   (I,ICMIN:K-1)   =     CNVNICE(ICMIN:K-1) !DONIF
!               CNV_FICE   (I,ICMIN:K-1)   =     CNVFICE(ICMIN:K-1) !DONIF
!! De-strap tendencies from RAS
!! specify weighting "SHAPE"
        wght = wgt1(i, :)
!! Scale properly by layer masses
        wght0 = 0.
        DO l=k,k0
          wght0 = wght0 + wght(l)*(ple(i, l+1)-ple(i, l))
        END DO
        wght0 = (prs(k+1)-prs(k))/wght0
        wght = wght0*wght
        DO l=k,k0
          tho(i, l) = tho(i, l) + wght(l)*(poi(k)-poi_sv(k))
          qho(i, l) = qho(i, l) + wght(l)*(qoi(k)-qoi_sv(k))
          uho(i, l) = uho(i, l) + wght(l)*(uoi(k)-uoi_sv(k))
          vho(i, l) = vho(i, l) + wght(l)*(voi(k)-voi_sv(k))
        END DO
        IF (do_tracers) THEN
          xho(i, icmin:k-1, :) = xoi(icmin:k-1, :)
          DO itr=1,itrcr
            DO l=k,k0
              xho(i, l, itr) = xho(i, l, itr) + wght(l)*(xoi(k, itr)-&
&               xoi_sv(k, itr))
            END DO
          END DO
        END IF
!            FLX (I,ICMIN:K) = RMF (ICMIN:K) * DDT/DAYLEN  !  (KG/m^2/s @ CLOUD BASE)
!  (KG/m^2/s @ CLOUD TOP)
        flxd(i, icmin:k) = rmfd(icmin:k)*ddt/daylen
!            FLXC(I,ICMIN:K) = RMFC(ICMIN:K) * DDT/DAYLEN  !  (KG/m^2/s @ CLOUD TOP)
!  (KG/m^2/s )
        clw(i, icmin:k) = cll(icmin:k)*ddt/daylen
        IF (PRESENT(disske)) disske(i, icmin:k-1) = dissk0(icmin:k-1)*&
&           ddt/daylen
!            FLX (I,1:ICMIN-1) = 0.
        flxd(i, 1:icmin-1) = 0.
!            FLXC(I,1:ICMIN-1) = 0.
        clw(i, 1:icmin-1) = 0.
        IF (k .LT. k0) THEN
!               FLX (I,K:K0) = 0.
          flxd(i, k:k0) = 0.
!               FLXC(I,K:K0) = 0.
          clw(i, k:k0) = 0.
        END IF
!            IRC (I,ICMIN:K-1) = RC(ICMIN:K-1)
        IF (ALLOCATED(icl_v)) THEN
          DEALLOCATE(icl_v)
        END IF
      END IF
      IF (final .EQ. 2) THEN
!            FLX (I,:) = 0.
        flxd(i, :) = 0.
!            FLXC(I,:) = 0.
        clw(i, :) = 0.
!            IRC (I,ICMIN:K-1) = RC(ICMIN:K-1)
      END IF
      RETURN
    END SUBROUTINE STRAP
  END SUBROUTINE RASE_BWD
!  Differentiation of sundq3_ice in reverse (adjoint) mode, forward sweep (with options r8 split(GEOS_MoistGridComp.PRE_RASE GEOS
!_MoistGridComp.PRE_PROGNO_CLOUD ras.RASE rase.CLOUDE rase.ACRITN rase.RNEVP rase.HTEST rase.FINDDTLS rase.STRAP ras.SUNDQ3_ICE c
!loudnew.pdf_spread cloudnew.fix_up_clouds cloudnew.meltfrz cloudnew.hystpdf cloudnew.pdffrac cloudnew.pdfcondensate cloudnew.cnv
!src cloudnew.evap3 cloudnew.subl3 cloudnew.autocon3 cloudnew.PRECIP3 cloudnew.ICEFALL cloudnew.SETTLE_VEL cloudnew.MARSHPALMQ2 c
!loudnew.MICRO_AA_BB_3 cloudnew.LDRADIUS3 cloudnew.ICE_FRACTION cloudnew.GET_ALHX3 cloudnew.ICEFRAC cloudnew.SUNDQ3_ICE3 qsat_uti
!l.DQSATPERT)):
!   gradient     of useful results: f2
!   with respect to varying inputs: temp f2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE SUNDQ3_ICE_FWD(temp, rate2, rate3, te1, f2, f3)
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: temp, rate2, rate3, te1
    REAL*8 :: f2, f3
!,RATE2,RATE3,TE1
    REAL*8 :: xx, yy, te0, te2, jump1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Ice - phase treatment totally invented
!!  Sharp increase in autoconversion in range
!!  ~~TE1 K ~< T < TE0 K .
!!  (JTB, 3/25/2003)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    te0 = 273.
    te2 = 200.
    jump1 = (rate2-1.0)/(te0-te1)**0.333
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Ice - phase treatment  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    IF (temp .GE. te0) THEN
      CALL PUSHREAL8(f2)
      f2 = 1.0
      CALL PUSHREAL8(f3)
      f3 = 1.0
      CALL PUSHCONTROL1B(0)
    ELSE
      CALL PUSHCONTROL1B(1)
    END IF
    IF (temp .GE. te1 .AND. temp .LT. te0) THEN
      CALL PUSHREAL8(f2)
      f2 = 1.0 + jump1*(te0-temp)**0.3333
      CALL PUSHREAL8(f3)
      f3 = 1.0
      CALL PUSHCONTROL1B(0)
    ELSE
      CALL PUSHCONTROL1B(1)
    END IF
    IF (temp .LT. te1) THEN
      CALL PUSHREAL8(f2)
      f2 = rate2 + (rate3-rate2)*(te1-temp)/(te1-te2)
      CALL PUSHREAL8(f3)
      f3 = 1.0
      CALL PUSHCONTROL1B(0)
    ELSE
      CALL PUSHCONTROL1B(1)
    END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IF (f2 .GT. 27.0) THEN
      CALL PUSHREAL8(f2)
      f2 = 27.0
      CALL PUSHCONTROL1B(1)
    ELSE
      CALL PUSHCONTROL1B(0)
    END IF
  END SUBROUTINE SUNDQ3_ICE_FWD
!  Differentiation of sundq3_ice in reverse (adjoint) mode, backward sweep (with options r8 split(GEOS_MoistGridComp.PRE_RASE GEO
!S_MoistGridComp.PRE_PROGNO_CLOUD ras.RASE rase.CLOUDE rase.ACRITN rase.RNEVP rase.HTEST rase.FINDDTLS rase.STRAP ras.SUNDQ3_ICE 
!cloudnew.pdf_spread cloudnew.fix_up_clouds cloudnew.meltfrz cloudnew.hystpdf cloudnew.pdffrac cloudnew.pdfcondensate cloudnew.cn
!vsrc cloudnew.evap3 cloudnew.subl3 cloudnew.autocon3 cloudnew.PRECIP3 cloudnew.ICEFALL cloudnew.SETTLE_VEL cloudnew.MARSHPALMQ2 
!cloudnew.MICRO_AA_BB_3 cloudnew.LDRADIUS3 cloudnew.ICE_FRACTION cloudnew.GET_ALHX3 cloudnew.ICEFRAC cloudnew.SUNDQ3_ICE3 qsat_ut
!il.DQSATPERT)):
!   gradient     of useful results: f2
!   with respect to varying inputs: temp f2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE SUNDQ3_ICE_BWD(temp, temp_ad, rate2, rate3, te1, f2, f2_ad&
&   , f3)
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: temp, rate2, rate3, te1
    REAL*8 :: temp_ad
    REAL*8 :: f2, f3
    REAL*8 :: f2_ad
    REAL*8 :: xx, yy, te0, te2, jump1
    INTEGER :: branch
    CALL POPCONTROL1B(branch)
    IF (branch .NE. 0) THEN
      CALL POPREAL8(f2)
      f2_ad = 0.0_8
    END IF
    te2 = 200.
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      CALL POPREAL8(f3)
      CALL POPREAL8(f2)
      temp_ad = -((rate3-rate2)*f2_ad/(te1-te2))
      f2_ad = 0.0_8
    ELSE
      temp_ad = 0.0_8
    END IF
    te0 = 273.
    jump1 = (rate2-1.0)/(te0-te1)**0.333
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      CALL POPREAL8(f3)
      CALL POPREAL8(f2)
      temp_ad = temp_ad - 0.3333*(te0-temp)**(-0.6667)*jump1*f2_ad
      f2_ad = 0.0_8
    END IF
    CALL POPCONTROL1B(branch)
    IF (branch .EQ. 0) THEN
      CALL POPREAL8(f3)
      CALL POPREAL8(f2)
      f2_ad = 0.0_8
    END IF
  END SUBROUTINE SUNDQ3_ICE_BWD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE SUNDQ3_ICE(temp, rate2, rate3, te1, f2, f3)
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: temp, rate2, rate3, te1
    REAL*8, INTENT(OUT) :: f2, f3
!,RATE2,RATE3,TE1
    REAL*8 :: xx, yy, te0, te2, jump1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Ice - phase treatment totally invented
!!  Sharp increase in autoconversion in range
!!  ~~TE1 K ~< T < TE0 K .
!!  (JTB, 3/25/2003)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    te0 = 273.
    te2 = 200.
    jump1 = (rate2-1.0)/(te0-te1)**0.333
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Ice - phase treatment  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    IF (temp .GE. te0) THEN
      f2 = 1.0
      f3 = 1.0
    END IF
    IF (temp .GE. te1 .AND. temp .LT. te0) THEN
      f2 = 1.0 + jump1*(te0-temp)**0.3333
      f3 = 1.0
    END IF
    IF (temp .LT. te1) THEN
      f2 = rate2 + (rate3-rate2)*(te1-temp)/(te1-te2)
      f3 = 1.0
    END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IF (f2 .GT. 27.0) f2 = 27.0
  END SUBROUTINE SUNDQ3_ICE

END MODULE RAS_ADM
