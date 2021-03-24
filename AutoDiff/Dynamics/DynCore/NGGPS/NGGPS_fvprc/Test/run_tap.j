#!/bin/csh

reset

source /discover/home/drholdaw/.cshrc

rm -rf tlm/
mkdir tlm/

$TAPENADE_HOME/bin/tapenade -O tlm/ -d -tgtvarname _tl -tgtfuncname _tlm -adjvarname _ad -adjfuncname _adm -head "group_call_mod.group_call_test (u2 v2 u3 v3)/(u2 v2 u3 v3)" group_call.F90 group_subs.F90




