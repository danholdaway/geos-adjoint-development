#!/bin/csh

reset

$TAPENADE_HOME/bin/tapenade -d -html -tgtvarname _tl -tgtfuncname _tlm -adjvarname _ad -adjfuncname _adm -head "ropp_fm_state2state_gsi_1d (temp shum pres geop)/(state)" ropp_s2s.f90 ropp_fm_constants.F90

