close all
clear
clc

cd /discover/nobackup/drholdaw/tmp.8469/sens.20111115.000000/

time = ncread('fsens.eta.111114_06z.nc4','time');
lon = ncread('fsens.eta.111114_06z.nc4','lon');
lat = ncread('fsens.eta.111114_06z.nc4','lat');

U_nm = ncread('fsens.eta.111114_06z.nc4','U');
V_nm = ncread('fsens.eta.111114_06z.nc4','V');
TV_nm = ncread('fsens.eta.111114_06z.nc4','TV');
DP_nm = ncread('fsens.eta.111114_06z.nc4','DP');
QV_nm = ncread('fsens.eta.111114_06z.nc4','QV');
QL_nm = ncread('fsens.eta.111114_06z.nc4','QL');
QI_nm = ncread('fsens.eta.111114_06z.nc4','QI');
O3_nm = ncread('fsens.eta.111114_06z.nc4','O3');

U_wm = ncread('fsens.eta.111114_06z.nc4','U');
V_wm = ncread('fsens.eta.111114_06z.nc4','V');
TV_wm = ncread('fsens.eta.111114_06z.nc4','TV');
DP_wm = ncread('fsens.eta.111114_06z.nc4','DP');
QV_wm = ncread('fsens.eta.111114_06z.nc4','QV');
QL_wm = ncread('fsens.eta.111114_06z.nc4','QL');
QI_wm = ncread('fsens.eta.111114_06z.nc4','QI');
O3_wm = ncread('fsens.eta.111114_06z.nc4','O3');

cd /home/drholdaw/moist_adjoint/

U_diff = U_nm - U_wm;
V_diff = V_nm - V_wm;
TV_diff = TV_nm - TV_wm;
DP_diff = DP_nm - DP_wm;
QV_diff = QV_nm - QV_wm;
QL_diff = QL_nm - QL_wm;
QI_diff = QI_nm - QI_wm;
O3_diff = O3_nm - O3_wm;


plotlevel = 70;


Umax = max(max([U_nm(:,:,plotlevel); U_wm(:,:,plotlevel)]));
Umin = min(min([U_nm(:,:,plotlevel); U_wm(:,:,plotlevel)]));
Vmax = max(max([V_nm(:,:,plotlevel); V_wm(:,:,plotlevel)]));
Vmin = min(min([V_nm(:,:,plotlevel); V_wm(:,:,plotlevel)]));
TVmax = max(max([TV_nm(:,:,plotlevel); TV_wm(:,:,plotlevel)]));
TVmin = min(min([TV_nm(:,:,plotlevel); TV_wm(:,:,plotlevel)]));
QVmax = max(max([QV_nm(:,:,plotlevel); QV_wm(:,:,plotlevel)]));
QVmin = min(min([QV_nm(:,:,plotlevel); QV_wm(:,:,plotlevel)]));
DPmax = max(max([DP_nm(:,:,plotlevel); DP_wm(:,:,plotlevel)]));
DPmin = min(min([DP_nm(:,:,plotlevel); DP_wm(:,:,plotlevel)]));

% zerocheck = max(max(max([QL_diff(:,:,plotlevel); QI_diff(:,:,plotlevel); O3_diff(:,:,plotlevel)])))

figure

subplot(2,3,1)
contour(lon,lat,(U_nm(:,:,plotlevel))')
% caxis([Umin Umax])
title('U \prime')
colorbar

subplot(2,3,2)
contour(lon,lat,V_nm(:,:,plotlevel)')
% caxis([Vmin Vmax])
title('V \prime')
colorbar

subplot(2,3,3)
contour(lon,lat,TV_nm(:,:,plotlevel)')
% caxis([TVmin TVmax])
title('T \prime')
colorbar

subplot(2,3,4)
contour(lon,lat,QV_nm(:,:,plotlevel)')
title('Q \prime')
colorbar
% caxis([QVmin QVmax])

subplot(2,3,5)
contour(lon,lat,DP_nm(:,:,plotlevel)')
title('P \prime')
colorbar
% caxis([DPmin DPmax])


figure

subplot(2,3,1)
contour(lon,lat,(U_wm(:,:,plotlevel))')
title('U \prime')
% caxis([Umin Umax])
colorbar

subplot(2,3,2)
contour(lon,lat,V_wm(:,:,plotlevel)')
title('V \prime')
% caxis([Vmin Vmax])
colorbar

subplot(2,3,3)
contour(lon,lat,TV_wm(:,:,plotlevel)')
title('T \prime')
% caxis([TVmin TVmax])
colorbar

subplot(2,3,4)
contour(lon,lat,QV_wm(:,:,plotlevel)')
title('Q \prime')
colorbar
% caxis([QVmin QVmax])

subplot(2,3,5)
contour(lon,lat,DP_wm(:,:,plotlevel)')
title('P \prime')
colorbar
% caxis([DPmin DPmax])


figure
subplot(2,3,1)
contour(lon,lat,U_diff(:,:,plotlevel)')
colorbar

subplot(2,3,2)
contour(lon,lat,V_diff(:,:,plotlevel)')
colorbar

subplot(2,3,3)
contour(lon,lat,TV_diff(:,:,plotlevel)')
colorbar

subplot(2,3,4)
contour(lon,lat,QV_diff(:,:,plotlevel)')
colorbar

subplot(2,3,5)
contour(lon,lat,DP_diff(:,:,plotlevel)')
colorbar