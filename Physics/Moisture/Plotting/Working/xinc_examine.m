close all
clear
clc

cd /discover/nobackup/drholdaw/tmp.8823/

U_bkg  = ncread('iau590.bkg.eta.20120317_00z.nc4','u');
V_bkg  = ncread('iau590.bkg.eta.20120317_00z.nc4','v');   
TV_bkg = ncread('iau590.bkg.eta.20120317_00z.nc4','tv');
DP_bkg = ncread('iau590.bkg.eta.20120317_00z.nc4','delp');
QV_bkg = ncread('iau590.bkg.eta.20120317_00z.nc4','sphu');
QL_bkg = ncread('iau590.bkg.eta.20120317_00z.nc4','qltot');
QI_bkg = ncread('iau590.bkg.eta.20120317_00z.nc4','qitot');
O3_bkg = ncread('iau590.bkg.eta.20120317_00z.nc4','ozone');

U_ana  = ncread('iau590.ana.eta.20120317_00z.nc4','u');
V_ana  = ncread('iau590.ana.eta.20120317_00z.nc4','v');   
TV_ana = ncread('iau590.ana.eta.20120317_00z.nc4','tv');
DP_ana = ncread('iau590.ana.eta.20120317_00z.nc4','delp');
QV_ana = ncread('iau590.ana.eta.20120317_00z.nc4','sphu');
QL_ana = ncread('iau590.ana.eta.20120317_00z.nc4','qltot');
QI_ana = ncread('iau590.ana.eta.20120317_00z.nc4','qitot');
O3_ana = ncread('iau590.ana.eta.20120317_00z.nc4','ozone');

U  = ncread('xinc.eta.nc4','u');
V  = ncread('xinc.eta.nc4','v');   
TV = ncread('xinc.eta.nc4','tv');
DP = ncread('xinc.eta.nc4','delp');
QV = ncread('xinc.eta.nc4','sphu');
QL = ncread('xinc.eta.nc4','qltot');
QI = ncread('xinc.eta.nc4','qitot');
O3 = ncread('xinc.eta.nc4','ozone');

U_0p01  = ncread('xinc_0.01.eta.nc4','u');
V_0p01  = ncread('xinc_0.01.eta.nc4','v');   
TV_0p01 = ncread('xinc_0.01.eta.nc4','tv');
DP_0p01 = ncread('xinc_0.01.eta.nc4','delp');
QV_0p01 = ncread('xinc_0.01.eta.nc4','sphu');
QL_0p01 = ncread('xinc_0.01.eta.nc4','qltot');
QI_0p01 = ncread('xinc_0.01.eta.nc4','qitot');
O3_0p01 = ncread('xinc_0.01.eta.nc4','ozone');

U_0p1  = ncread('xinc_0.1.eta.nc4','u');
V_0p1  = ncread('xinc_0.1.eta.nc4','v');   
TV_0p1 = ncread('xinc_0.1.eta.nc4','tv');
DP_0p1 = ncread('xinc_0.1.eta.nc4','delp');
QV_0p1 = ncread('xinc_0.1.eta.nc4','sphu');
QL_0p1 = ncread('xinc_0.1.eta.nc4','qltot');
QI_0p1 = ncread('xinc_0.1.eta.nc4','qitot');
O3_0p1 = ncread('xinc_0.1.eta.nc4','ozone');

U_0p2  = ncread('xinc_0.2.eta.nc4','u');
V_0p2  = ncread('xinc_0.2.eta.nc4','v');   
TV_0p2 = ncread('xinc_0.2.eta.nc4','tv');
DP_0p2 = ncread('xinc_0.2.eta.nc4','delp');
QV_0p2 = ncread('xinc_0.2.eta.nc4','sphu');
QL_0p2 = ncread('xinc_0.2.eta.nc4','qltot');
QI_0p2 = ncread('xinc_0.2.eta.nc4','qitot');
O3_0p2 = ncread('xinc_0.2.eta.nc4','ozone');

U_0p5  = ncread('xinc_0.5.eta.nc4','u');
V_0p5  = ncread('xinc_0.5.eta.nc4','v');   
TV_0p5 = ncread('xinc_0.5.eta.nc4','tv');
DP_0p5 = ncread('xinc_0.5.eta.nc4','delp');
QV_0p5 = ncread('xinc_0.5.eta.nc4','sphu');
QL_0p5 = ncread('xinc_0.5.eta.nc4','qltot');
QI_0p5 = ncread('xinc_0.5.eta.nc4','qitot');
O3_0p5 = ncread('xinc_0.5.eta.nc4','ozone');

U_1p0  = ncread('xinc_1.0.eta.nc4','u');
V_1p0  = ncread('xinc_1.0.eta.nc4','v');   
TV_1p0 = ncread('xinc_1.0.eta.nc4','tv');
DP_1p0 = ncread('xinc_1.0.eta.nc4','delp');
QV_1p0 = ncread('xinc_1.0.eta.nc4','sphu');
QL_1p0 = ncread('xinc_1.0.eta.nc4','qltot');
QI_1p0 = ncread('xinc_1.0.eta.nc4','qitot');
O3_1p0 = ncread('xinc_1.0.eta.nc4','ozone');

plot_level = 50;

contour(TV_ana(:,:,plot_level)'-TV_bkg(:,:,plot_level)')
colorbar

figure
contour(TV(:,:,plot_level)')
colorbar

figure
contour(TV_0p01(:,:,plot_level)')
colorbar

figure
contour(TV_0p1(:,:,plot_level)')
colorbar

figure
contour(TV_0p2(:,:,plot_level)')
colorbar

figure
contour(TV_0p5(:,:,plot_level)')
colorbar

figure
contour(TV_1p0(:,:,plot_level)')
colorbar


cd /home/drholdaw/Lin_Moist_Physics/DAS_plotting/







