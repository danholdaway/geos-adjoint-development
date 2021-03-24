#!/bin/csh

reset

rm -f sorad_tl.F90 sorad_ad.F90 sorad_ad.F90

#Simplified
$TAPENADE_HOME/bin/tapenade -tangent -head sorad -outvars "flx_dev" -vars "ta_dev wa_dev oa_dev cwc_dev fcld_dev reff_dev flx_dev taua_dev ssaa_dev asya_dev" sorad.f90 -output sorad -html

$TAPENADE_HOME/bin/tapenade -reverse -head sorad -outvars "flx_dev" -vars "ta_dev wa_dev oa_dev cwc_dev fcld_dev reff_dev flx_dev taua_dev ssaa_dev asya_dev" sorad.f90 -output sorad -html


cat headfoot/tl_top.F90 sorad_d.f90 headfoot/tl_bot.F90 > sorad_tl.F90
cat headfoot/ad_top.F90 sorad_b.f90 headfoot/ad_bot.F90 > sorad_ad.F90
cat headfoot/main_top.F90 sorad.f90 headfoot/main_bot.F90 > sorad.F90

rm -f sorad_b.f90 sorad_d.f90

cp -r sorad_tl.F90 $radpertdirdust
cp -r sorad_ad.F90 $radpertdirdust
cp -r sorad.F90    $radpertdirdust

