#!/bin/csh

reset

$TAPENADE_HOME/bin/tapenade -tangent -head convection_driver -outvars "U V TH Q CNV_DQLDT CNV_MFD CNV_PRC3 CNV_UPDF" -vars "U V TH Q" convection_driver.f90 -output convection -html

#sed -i '1,/END SUBROUTINE DUMMY_D/d' convection_d.f90

cat headfoot/tl_top.f90 convection_d.f90 headfoot/tl_bot.f90 > convection_tl.F90
rm convection_d.f90

cp -r tapenadehtml log_tl


#firefox log_tl/tapenade.html &

$TAPENADE_HOME/bin/tapenade -reverse -head convection_driver -outvars "U V TH Q CNV_DQLDT CNV_MFD CNV_PRC3 CNV_UPDF" -vars "U V TH Q CNV_DQLDT CNV_MFD CNV_PRC3 CNV_UPDF" convection_driver.f90 -output convection -html

#sed -i '1,/END SUBROUTINE DUMMY_B/d' convection_b.f90

cat headfoot/ad_top.f90 convection_b.f90 headfoot/ad_bot.f90 > convection_ad.F90
rm convection_b.f90

cp -r tapenadehtml log_ad

#rm -r tapenadehtml
#rm convection_d.msg convection_b.msg

#firefox log_ad/tapenade.html &

rm *~


