#!/bin/csh

reset

source /discover/home/drholdaw/.cshrc

#Tapenade options
set opts = "-r8 -tgtvarname _tl -tgtfuncname _tlm -adjvarname _ad -adjfuncname _adm"


$TAPENADE_HOME/bin/tapenade ${opts} -b -head "dyn_core_mod.dyn_core (u v)/(u v)" adto.F90 array.F90 -nocheckpoint "dyn_core_mod.dyn_core dyn_core_mod.again"
