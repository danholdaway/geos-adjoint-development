#!/bin/csh

reset

$TAPENADE_HOME/bin/tapenade -r8 -d -tgtvarname _tl -tgtfuncname _tlm -adjvarname _ad -adjfuncname _adm -head "GEOS_MoistGridComp.moist_run (TH1 Q1 U1 V1 QLLS QLCN QILS QICN CLCN CLLS)/(TH1 Q1 U1 V1 QLLS QLCN QILS QICN CLCN CLLS)" GEOS_MoistGridComp.F90 MAPL_Constants.F90 qsat_util.F90 CLDPARAMS.F90 RASPARAMS.F90 ras.F90 cloudnew.F90

rename _d.f90 _tlm.F90 *_d.f90

#$TAPENADE_HOME/bin/tapenade -r8 -b -tgtvarname _tl -tgtfuncname _tlm -adjvarname _ad -adjfuncname _adm -head "GEOS_MoistGridComp.moist_run (TH1 Q1 U1 V1 QLLS QLCN QILS QICN CLCN CLLS)/(TH1 Q1 U1 V1 QLLS QLCN QILS QICN CLCN CLLS)" GEOS_MoistGridComp.F90 MAPL_Constants.F90 qsat_util.F90 CLDPARAMS.F90 RASPARAMS.F90 ras.F90 cloudnew.F90 -nocheckpoint "GEOS_MoistGridComp.PRE_RASE GEOS_MoistGridComp.PRE_PROGNO_CLOUD ras.RASE rase.CLOUDE rase.ACRITN rase.RNEVP rase.HTEST rase.FINDDTLS rase.STRAP ras.SUNDQ3_ICE cloudnew.pdf_spread cloudnew.fix_up_clouds cloudnew.meltfrz cloudnew.hystpdf cloudnew.pdffrac cloudnew.pdfcondensate cloudnew.cnvsrc cloudnew.evap3 cloudnew.subl3 cloudnew.autocon3 cloudnew.PRECIP3 cloudnew.ICEFALL cloudnew.SETTLE_VEL cloudnew.MARSHPALMQ2 cloudnew.MICRO_AA_BB_3 cloudnew.LDRADIUS3 cloudnew.ICE_FRACTION cloudnew.GET_ALHX3 cloudnew.ICEFRAC cloudnew.SUNDQ3_ICE3"

$TAPENADE_HOME/bin/tapenade -r8 -b -tgtvarname _tl -tgtfuncname _tlm -adjvarname _ad -adjfuncname _adm -head "GEOS_MoistGridComp.moist_run (TH1 Q1 U1 V1 QLLS QLCN QILS QICN CLCN CLLS)/(TH1 Q1 U1 V1 QLLS QLCN QILS QICN CLCN CLLS)" GEOS_MoistGridComp.F90 MAPL_Constants.F90 qsat_util.F90 CLDPARAMS.F90 RASPARAMS.F90 ras.F90 cloudnew.F90 -nocheckpoint GEOS_MoistGridComp.PRE_PROGNO_CLOUD

rename _b.f90 _adm.F90 *_b.f90

