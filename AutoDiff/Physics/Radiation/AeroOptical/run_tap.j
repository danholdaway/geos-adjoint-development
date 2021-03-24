#!/bin/csh

reset

$TAPENADE_HOME/bin/tapenade -tangent -head aeroopt -outvars "TAUA_IRT SSAA_IRT ASYA_IRT TAUA_SOT SSAA_SOT ASYA_SOT" -vars "AEROT" aerooptical.f90 -output aeroopt -html

$TAPENADE_HOME/bin/tapenade -reverse -head aeroopt -outvars "TAUA_IRT SSAA_IRT ASYA_IRT TAUA_SOT SSAA_SOT ASYA_SOT" -vars "AEROT" aerooptical.f90 -output aeroopt -html


