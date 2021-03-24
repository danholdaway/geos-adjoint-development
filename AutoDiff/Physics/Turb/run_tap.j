#!/bin/csh

reset

$TAPENADE_HOME/bin/tapenade -d -tgtvarname _tl -tgtfuncname _tlm -adjvarname _ad -adjfuncname _adm -head "G5TDSOLVER (Y)/(Y)" solver.f90 

$TAPENADE_HOME/bin/tapenade -b -tgtvarname _tl -tgtfuncname _tlm -adjvarname _ad -adjfuncname _adm -head "G5TDSOLVER (Y)/(Y)" solver.f90 


