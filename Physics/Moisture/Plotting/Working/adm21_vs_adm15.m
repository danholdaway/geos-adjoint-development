close all
clear
clc

cd /discover/nobackup/drholdaw/tmp.16412/sens.20120318.000000/

lon  = ncread('fsens.eta.20120318_00z.nc4','lon');
lat  = ncread('fsens.eta.20120318_00z.nc4','lat');

U_15z  = ncread('fsens.eta.20120318_00z.nc4','u');
V_15z  = ncread('fsens.eta.20120318_00z.nc4','v');
TV_15z = ncread('fsens.eta.20120318_00z.nc4','tv');
DP_15z = ncread('fsens.eta.20120318_00z.nc4','delp');
QV_15z = ncread('fsens.eta.20120318_00z.nc4','sphu');
QL_15z = ncread('fsens.eta.20120318_00z.nc4','qltot');
QI_15z = ncread('fsens.eta.20120318_00z.nc4','qitot');
O3_15z = ncread('fsens.eta.20120318_00z.nc4','ozone');

cd /discover/nobackup/drholdaw/tmp.30782/sens.20120318.000000/

U_21z  = ncread('fsens.eta.20120318_00z.nc4','u');
V_21z  = ncread('fsens.eta.20120318_00z.nc4','v');
TV_21z = ncread('fsens.eta.20120318_00z.nc4','tv');
DP_21z = ncread('fsens.eta.20120318_00z.nc4','delp');
QV_21z = ncread('fsens.eta.20120318_00z.nc4','sphu');
QL_21z = ncread('fsens.eta.20120318_00z.nc4','qltot');
QI_21z = ncread('fsens.eta.20120318_00z.nc4','qitot');
O3_21z = ncread('fsens.eta.20120318_00z.nc4','ozone');

cd /home/drholdaw/Lin_Moist_Physics/DAS_plotting/

plot_level = 50;

%GET [C,h] = contour LIMITS
Umax = max(max([U_15z(:,:,plot_level); U_21z(:,:,plot_level)]));
Umin = min(min([U_15z(:,:,plot_level); U_21z(:,:,plot_level)]));
Vmax = max(max([V_15z(:,:,plot_level); V_21z(:,:,plot_level)]));
Vmin = min(min([V_15z(:,:,plot_level); V_21z(:,:,plot_level)]));
Tmax = max(max([TV_15z(:,:,plot_level); TV_21z(:,:,plot_level)]));
Tmin = min(min([TV_15z(:,:,plot_level); TV_21z(:,:,plot_level)]));
Qmax = max(max([QV_15z(:,:,plot_level); QV_21z(:,:,plot_level)]));
Qmin = min(min([QV_15z(:,:,plot_level); QV_21z(:,:,plot_level)]));
Pmax = max(max([DP_15z(:,:,plot_level); DP_21z(:,:,plot_level)]));
Pmin = min(min([DP_15z(:,:,plot_level); DP_21z(:,:,plot_level)]));



figure
subplot(3,2,1)
[C,h] = contourf(lon,lat,U_15z(:,:,plot_level)','LineStyle','none');
caxis([Umin Umax])
cintu = get(h,'LevelStep');
colorbar

subplot(3,2,2)
[C,h] = contourf(lon,lat,V_15z(:,:,plot_level)','LineStyle','none');
caxis([Vmin Vmax])
cintv = get(h,'LevelStep');
colorbar

subplot(3,2,3)
[C,h] = contourf(lon,lat,TV_15z(:,:,plot_level)','LineStyle','none');
caxis([Tmin Tmax])
cintt = get(h,'LevelStep');
colorbar

subplot(3,2,4)
[C,h] = contourf(lon,lat,QV_15z(:,:,plot_level)','LineStyle','none');
caxis([Qmin Qmax])
cintq = get(h,'LevelStep');
colorbar

subplot(3,2,5)
[C,h] = contourf(lon,lat,DP_15z(:,:,plot_level)','LineStyle','none');
caxis([Pmin Pmax])
cintp = get(h,'LevelStep');
colorbar

figure
subplot(3,2,1)
[C,h] = contourf(lon,lat,U_21z(:,:,plot_level)','LineStyle','none','LevelStep',cintu);
caxis([Umin Umax])
colorbar

subplot(3,2,2)
[C,h] = contourf(lon,lat,V_21z(:,:,plot_level)','LineStyle','none','LevelStep',cintv);
caxis([Vmin Vmax])
colorbar

subplot(3,2,3)
[C,h] = contourf(lon,lat,TV_21z(:,:,plot_level)','LineStyle','none','LevelStep',cintt);
caxis([Tmin Tmax])
colorbar

subplot(3,2,4)
[C,h] = contourf(lon,lat,QV_21z(:,:,plot_level)','LineStyle','none','LevelStep',cintq);
caxis([Qmin Qmax])
colorbar

subplot(3,2,5)
[C,h] = contourf(lon,lat,DP_21z(:,:,plot_level)','LineStyle','none','LevelStep',cintp);
caxis([Pmin Pmax])
colorbar