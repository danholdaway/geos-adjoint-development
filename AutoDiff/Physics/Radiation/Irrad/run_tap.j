#!/bin/csh

reset

rm -f irrad_tl.F90 irrad_ad.F90 irrad_ad.F90

#Full
#$TAPENADE_HOME/bin/tapenade -tangent -head irrad -outvars "taua_dev ssaa_dev asya_dev flxu_dev flcu_dev flau_dev flxd_dev flcd_dev flad_dev dfdts_dev sfcem_dev taudiag_dev" -vars "ta_dev wa_dev oa_dev tb_dev n2o_dev ch4_dev cfc11_dev cfc12_dev cfc22_dev cwc_dev fcld_dev reff_dev fs_dev tg_dev eg_dev tv_dev ev_dev rv_dev taua_dev ssaa_dev asya_dev" irrad.F90 -output irrad -html

#$TAPENADE_HOME/bin/tapenade -reverse -head irrad -outvars "taua_dev ssaa_dev asya_dev flxu_dev flcu_dev flau_dev flxd_dev flcd_dev flad_dev dfdts_dev sfcem_dev taudiag_dev" -vars "ta_dev wa_dev oa_dev tb_dev n2o_dev ch4_dev cfc11_dev cfc12_dev cfc22_dev cwc_dev fcld_dev reff_dev fs_dev tg_dev eg_dev tv_dev ev_dev rv_dev taua_dev ssaa_dev asya_dev" irrad.F90 -output irrad -html


#Simplified
$TAPENADE_HOME/bin/tapenade -tangent -head irrad -outvars "flxu_dev flxd_dev dfdts_dev" -vars "ta_dev wa_dev oa_dev tb_dev cwc_dev fcld_dev reff_dev flxu_dev flxd_dev dfdts_dev taua_dev ssaa_dev asya_dev" irrad.f90 -output irrad -html

$TAPENADE_HOME/bin/tapenade -reverse -head irrad -outvars "flxu_dev flxd_dev dfdts_dev" -vars "ta_dev wa_dev oa_dev tb_dev cwc_dev fcld_dev reff_dev flxu_dev flxd_dev dfdts_dev taua_dev ssaa_dev asya_dev" irrad.f90 -output irrad -html


cat headfoot/tl_top.F90 irrad_d.f90 headfoot/tl_bot.F90 > irrad_tl.F90
cat headfoot/ad_top.F90 irrad_b.f90 headfoot/ad_bot.F90 > irrad_ad.F90
cat headfoot/main_top.F90 irrad.f90 headfoot/main_bot.F90 > irrad.F90

rm -f irrad_b.f90 irrad_d.f90

cp irrad_tl.F90 $radpertdirdust
cp irrad_ad.F90 $radpertdirdust
cp irrad.F90    $radpertdirdust

