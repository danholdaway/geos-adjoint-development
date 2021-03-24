#!/bin/csh

reset

$TAPENADE_HOME/bin/tapenade -tangent -head fyppmlinear -outvars "flux" -vars "c q" fyppmlinear.f90 -output fyppmlinear -html

#$TAPENADE_HOME/bin/tapenade -tangent -head fxppm -outvars "flux" -vars "c q" fxppm.f90 -output fxppm -html


