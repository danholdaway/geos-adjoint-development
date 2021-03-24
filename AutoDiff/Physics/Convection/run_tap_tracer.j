#!/bin/csh

reset

rm -f convection_tracer_tl.F90 convection_tracer_ad.F90

#Simplified
$TAPENADE_HOME/bin/tapenade -forward -inputlanguage fortran95 -outputlanguage fortran95 -head rase_tracer -outvars "XHO" -vars "XHO" convection_tracer.f95 -output convection_tracer -html

$TAPENADE_HOME/bin/tapenade -reverse -inputlanguage fortran95 -outputlanguage fortran95 -head rase_tracer -outvars "XHO" -vars "XHO" convection_tracer.f95 -output convection_tracer -html

#cat headfoot/tl_top.F90 convection_tracer_d.f95 headfoot/tl_bot.F90 > convection_tracer_tl.F90
#cat headfoot/ad_top.F90 convection_tracer_b.f95 headfoot/ad_bot.F90 > convection_tracer_ad.F90
#cat headfoot/main_top.F90 convection_tracer.f95 headfoot/main_bot.F90 > convection_tracer.F90

#rm -f convection_tracer_b.f95 convection_tracer_d.f95

#sed -i 's/RASE/RASE_TLAD/g' convection.F90

#cp -r convection_tl.F90 $moistpertdir
#cp -r convection_ad.F90 $moistpertdir
#cp -r convection.F90    $moistpertdir
