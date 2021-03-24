#!/bin/csh

reset


$TAPENADE_HOME/bin/tapenade -tangent -head radiation_driver_dummy -outvars "PTT QVT O3T CFLST CFCNT QIT QLT DU001T DU002T DU003T DU004T DU005T" -vars "PTT QVT O3T CFLST CFCNT QIT QLT DU001T DU002T DU003T DU004T DU005T" driver.f95 -output radiation_driver -html

$TAPENADE_HOME/bin/tapenade -reverse -head radiation_driver_dummy -outvars "PTT QVT O3T CFLST CFCNT QIT QLT DU001T DU002T DU003T DU004T DU005T" -vars "PTT QVT O3T CFLST CFCNT QIT QLT DU001T DU002T DU003T DU004T DU005T" driver.f95 -output radiation_driver -html



