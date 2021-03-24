#!/bin/csh

reset

rm -f convection_tl.F90 convection_ad.F90 convection_ad.F90

#Simplified
$TAPENADE_HOME/bin/tapenade -forward -inputlanguage fortran95 -outputlanguage fortran95 -head rase -outvars "THO QHO UHO VHO CLW FLXD CNV_PRC3 CNV_UPDFRC" -vars "THO QHO UHO VHO" convection.f95 -output convection -html

$TAPENADE_HOME/bin/tapenade -reverse -inputlanguage fortran95 -outputlanguage fortran95 -head rase -outvars "THO QHO UHO VHO CLW FLXD CNV_PRC3 CNV_UPDFRC" -vars "THO QHO UHO VHO CLW FLXD CNV_PRC3 CNV_UPDFRC" convection.f95 -output convection -html

cat headfoot/tl_top.F90 convection_d.f95 headfoot/tl_bot.F90 > convection_tl.F90
cat headfoot/ad_top.F90 convection_b.f95 headfoot/ad_bot.F90 > convection_ad.F90
cat headfoot/main_top.F90 convection.f95 headfoot/main_bot.F90 > convection.F90

rm -f convection_b.f95 convection_d.f95

#sed -i 's/RASE/RASE_TLAD/g' convection.F90

#cp -r convection_tl.F90 $moistpertdir
#cp -r convection_ad.F90 $moistpertdir
#cp -r convection.F90    $moistpertdir
