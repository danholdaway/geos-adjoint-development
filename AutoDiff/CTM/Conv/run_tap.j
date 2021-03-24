#!/bin/csh

reset

source /discover/home/drholdaw/.cshrc

#Tapenade options
set opts = "-html -tgtvarname _tl -tgtfuncname _tlm -adjvarname _ad -adjfuncname _adm"


#Generate tangent linear code
#----------------------------
rm -f tlm/*.F90

$TAPENADE_HOME/bin/tapenade ${opts} -d -O tlm/ -head "GenericConvectionMethod_mod.doGenericConvectiveTransport (concentration%pArray3D)/(concentration%pArray3D)" GenericConvectionMethod_mod.F90 convectiveTransport_mod.F90 ESMF.F90 GmiArrayBundlePointer_mod.F90

cd tlm/
rename mod_d.f90 tlm.F90 *mod_d.f90
cd ../


#Generate adjoint code
#---------------------
rm -f adm/*.F90

$TAPENADE_HOME/bin/tapenade ${opts} -b -O adm/ -head "GenericConvectionMethod_mod.doGenericConvectiveTransport (concentration%pArray3D)/(concentration%pArray3D)" GenericConvectionMethod_mod.F90 convectiveTransport_mod.F90 ESMF.F90 GmiArrayBundlePointer_mod.F90

cd adm/
rename mod_b.f90 adm.F90 *mod_b.f90
cd ../



