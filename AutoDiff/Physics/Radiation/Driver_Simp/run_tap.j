#!/bin/csh

reset


$TAPENADE_HOME/bin/tapenade -tangent -head radiation_driver -outvars "PTT" -vars "PTT" driver.f95 -output radiation_driver -html

$TAPENADE_HOME/bin/tapenade -reverse -head radiation_driver -outvars "PTT" -vars "PTT" driver.f95 -output radiation_driver -html


rm -rf tapenadehtml
