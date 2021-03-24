#!/bin/csh

reset

$TAPENADE_HOME/bin/tapenade -tangent -head FILL_Friendly -outvars "Q" -vars "Q" fill_friendly.f90 -output FILL_Friendly -html

$TAPENADE_HOME/bin/tapenade -reverse -head FILL_Friendly -outvars "Q" -vars "Q" fill_friendly.f90 -output FILL_Friendly -html

