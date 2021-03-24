#!/bin/csh

reset

$TAPENADE_HOME/bin/tapenade -tangent -head v2_mapz_q -outvars "q2" -vars "pe1 pe2 dp1 q4 q2" v2_mapz.f90 -output v2_mapz_q -html


#$TAPENADE_HOME/bin/tapenade -reverse -head v2_cs_profile_linear -outvars "a4" -vars "a4 delp" cs_profile.f90 -output v2_cs_profile_linear -html
