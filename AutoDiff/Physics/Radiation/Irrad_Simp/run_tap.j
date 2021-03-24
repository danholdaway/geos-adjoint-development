#!/bin/csh

reset

rm -f irrad_tl.F90 irrad_ad.F90 irrad_ad.F90

#Simplified
$TAPENADE_HOME/bin/tapenade -tangent -head irrad -outvars "flxu_dev flxd_dev dfdts_dev" -vars "ta_dev" irrad.f90 -output irrad -html

$TAPENADE_HOME/bin/tapenade -reverse -head irrad -outvars "flxu_dev flxd_dev dfdts_dev" -vars "ta_dev" irrad.f90 -output irrad -html


cat headfoot/tl_top.F90 irrad_d.f90 headfoot/tl_bot.F90 > irrad_tl.F90
cat headfoot/ad_top.F90 irrad_b.f90 headfoot/ad_bot.F90 > irrad_ad.F90

rm -f irrad_b.f90 irrad_d.f90

rm -r tapenadehtml

#cp irrad_tl.F90 $radpertdirdust
#cp irrad_ad.F90 $radpertdirdust

