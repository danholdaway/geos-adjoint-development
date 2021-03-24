#!/bin/csh

reset

$TAPENADE_HOME/bin/tapenade -forward -inputlanguage fortran95 -outputlanguage fortran95 -head cloud_driver -outvars "th q QI_LS QL_LS QI_con QL_con CF_con CF_ls" -vars "th q QI_LS QL_LS QI_con QL_con CF_con CF_ls CNV_DQLDT CNV_MFD CNV_PRC3 CNV_UPDF" cloud.f95 -output cloud -html

# COMBINE WITH THE FILTERING PARTS
# --------------------------------
#sed -e '/FILTERINGINSERT1/{r filt/tl1' -e '}' cloud_d.f95 > cloud_d1.f95
#sed -e '/also update convection clouds/{r filt/tl2' -e '}' cloud_d1.f95 > cloud_d2.f95
#sed -e '/&                 , t_ice_all, t_ice_max, icefrpwr, estblx)/{r filt/tl3' -e '}' cloud_d2.f95 > cloud_d3.f95
#sed -e '/FILTERINGINSERT2/{r filt/tl4' -e '}' cloud_d3.f95 > cloud_d4.f95
#rm cloud_d1.f95 cloud_d2.f95 cloud_d3.f95
#mv cloud_d4.f95 cloud_d.f95

# Combine with the header and footers
# -----------------------------------
cat headfoot/tl_top.f95 cloud_d.f95 headfoot/tl_bot.f95 > cloud_tl.F90

# Remove the original
# -------------------
rm cloud_d.f95

# Copy the html logs somewhere else
# ---------------------------------
#cp -r tapenadehtml log_tl


#firefox log_tl/tapenade.html &

$TAPENADE_HOME/bin/tapenade -reverse -inputlanguage fortran95 -outputlanguage fortran95 -head cloud_driver -outvars "th q QI_LS QL_LS QI_con QL_con CF_con CF_ls" -vars "th q QI_LS QL_LS QI_con QL_con CF_con CF_ls CNV_DQLDT CNV_MFD CNV_PRC3 CNV_UPDF" cloud.f95 -output cloud -html

# COMBINE WITH THE FILTERING PARTS
# --------------------------------
#sed -e '/FILTERINGINSERT1/{r filt/tl1' -e '}' cloud_b.f95 > cloud_b1.f95
#sed -e '/CALL LS_CLOUD_B/{r filt/ad2' -e 'b}' cloud_b1.f95 > cloud_b2.f95
#sed -e '/&                 , t_ice_all, t_ice_max, icefrpwr, estblx)/{r filt/ad3' -e '}' cloud_b2.f95 > cloud_b3.f95
#sed -i '/CALL LS_CLOUD_B/d' cloud_b3.f95
#sed -e '/cf_conb(i,j,k) = filt_lsc/{r filt/ad4' -e '}' cloud_b3.f95 > cloud_b4.f95
#sed -e '/thb = 0.0_8/{r filt/ad5' -e '}' cloud_b4.f95 > cloud_b5.f95
#sed -e '/tb = thb/{r filt/ad6' -e '}' cloud_b5.f95 > cloud_b6.f95

#rm cloud_b1.f95 cloud_b2.f95 cloud_b3.f95 cloud_b4.f95 cloud_b5.f95
#mv cloud_b6.f95 cloud_b.f95

# Combine with the header and footers
# -----------------------------------
cat headfoot/ad_top.f95 cloud_b.f95 headfoot/ad_bot.f95 > cloud_ad.F90

# Remove the original
# -------------------
rm cloud_b.f95

# Combine with the header and footers
# -----------------------------------
#cat headfoot/main_top.f95 cloud.f95 headfoot/main_bot.f95 > cloud.F90

# Copy the html logs somewhere else
# ---------------------------------
#cp -r tapenadehtml log_ad


# Remove the html logs
# --------------------
rm -r tapenadehtml

# Remove the messages
# -------------------
#rm cloud_d.msg cloud_b.msg

# Remove tmp files
# ----------------
rm *~

# Copy to the moistpert directory
# -------------------------------

#cp cloud_ad.F90 $moistpertdir
#cp cloud_tl.F90 $moistpertdir
#cp cloud.F90 $moistpertdir
