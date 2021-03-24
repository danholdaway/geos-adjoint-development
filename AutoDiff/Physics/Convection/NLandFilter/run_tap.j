#!/bin/csh

reset

$TAPENADE_HOME/bin/tapenade -forward -inputlanguage fortran95 -outputlanguage fortran95 -head rase0 -outvars "THO QHO" -vars "THO QHO" convection.f95 -output convection -html

