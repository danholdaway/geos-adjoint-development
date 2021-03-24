close all
clear
clc

% DOT PRODUCT TEST RESULT
% (x,x)    =   1.971796897610801E-005
% (Tx,Tx)  =   1.480459420116413E-003
% (T'Tx,x) =   1.480459703412876E-003
% rel error =   1.913571283877516E-007

load coast
coast_lat = lat; clear lat
coast_lon = long; clear long
grey = 0.0;

cd /discover/nobackup/drholdaw/tmp.11265/sens.20120411.000000/

lat = ncread('JgradfX.eta.nc4','lat');
lon = ncread('JgradfX.eta.nc4','lon');

tlm_20120410_1200z_u = ncread('fsens.tlm.20120410_12z.nc4','u');
tlm_20120410_1200z_v = ncread('fsens.tlm.20120410_12z.nc4','v');
tlm_20120410_1200z_t = ncread('fsens.tlm.20120410_12z.nc4','tv');
tlm_20120410_1200z_q = ncread('fsens.tlm.20120410_12z.nc4','sphu');
tlm_20120410_1200z_p = ncread('fsens.tlm.20120410_12z.nc4','delp');

tlm_20120410_1400z_u = ncread('fsens.tlm.20120410_14z.nc4','u');
tlm_20120410_1400z_v = ncread('fsens.tlm.20120410_14z.nc4','v');
tlm_20120410_1400z_t = ncread('fsens.tlm.20120410_14z.nc4','tv');
tlm_20120410_1400z_q = ncread('fsens.tlm.20120410_14z.nc4','sphu');
tlm_20120410_1400z_p = ncread('fsens.tlm.20120410_14z.nc4','delp');

tlm_20120410_1500z_u = ncread('fsens.tlm.20120410_15z.nc4','u');
tlm_20120410_1500z_v = ncread('fsens.tlm.20120410_15z.nc4','v');
tlm_20120410_1500z_t = ncread('fsens.tlm.20120410_15z.nc4','tv');
tlm_20120410_1500z_q = ncread('fsens.tlm.20120410_15z.nc4','sphu');
tlm_20120410_1500z_p = ncread('fsens.tlm.20120410_15z.nc4','delp');

tlm_20120410_1520z_u = ncread('fsens.tlm.20120410_1520z.nc4','u');
tlm_20120410_1520z_v = ncread('fsens.tlm.20120410_1520z.nc4','v');
tlm_20120410_1520z_t = ncread('fsens.tlm.20120410_1520z.nc4','tv');
tlm_20120410_1520z_q = ncread('fsens.tlm.20120410_1520z.nc4','sphu');
tlm_20120410_1520z_p = ncread('fsens.tlm.20120410_1520z.nc4','delp');

tlm_20120410_1540z_u = ncread('fsens.tlm.20120410_1540z.nc4','u');
tlm_20120410_1540z_v = ncread('fsens.tlm.20120410_1540z.nc4','v');
tlm_20120410_1540z_t = ncread('fsens.tlm.20120410_1540z.nc4','tv');
tlm_20120410_1540z_q = ncread('fsens.tlm.20120410_1540z.nc4','sphu');
tlm_20120410_1540z_p = ncread('fsens.tlm.20120410_1540z.nc4','delp');

tlm_20120410_1600z_u = ncread('fsens.tlm.20120410_16z.nc4','u');
tlm_20120410_1600z_v = ncread('fsens.tlm.20120410_16z.nc4','v');
tlm_20120410_1600z_t = ncread('fsens.tlm.20120410_16z.nc4','tv');
tlm_20120410_1600z_q = ncread('fsens.tlm.20120410_16z.nc4','sphu');
tlm_20120410_1600z_p = ncread('fsens.tlm.20120410_16z.nc4','delp');

tlm_20120410_1800z_u = ncread('fsens.tlm.20120410_18z.nc4','u');
tlm_20120410_1800z_v = ncread('fsens.tlm.20120410_18z.nc4','v');
tlm_20120410_1800z_t = ncread('fsens.tlm.20120410_18z.nc4','tv');
tlm_20120410_1800z_q = ncread('fsens.tlm.20120410_18z.nc4','sphu');
tlm_20120410_1800z_p = ncread('fsens.tlm.20120410_18z.nc4','delp');

% tlm_20120411_0000z_u = ncread('fsens.tlm.20120411_00z.nc4','u');
% tlm_20120411_0000z_v = ncread('fsens.tlm.20120411_00z.nc4','v');
% tlm_20120411_0000z_t = ncread('fsens.tlm.20120411_00z.nc4','tv');
% tlm_20120411_0000z_q = ncread('fsens.tlm.20120411_00z.nc4','sphu');
% tlm_20120411_0000z_p = ncread('fsens.tlm.20120411_00z.nc4','delp');

cd /home/drholdaw/Lin_Moist_Physics/DAS_plotting/

plotlevel = 50;

figure('position',[151 236 1128 683])
% subplot(3,2,1)
contourf(lon,lat,tlm_20120410_1200z_u(:,:,plotlevel)','LineStyle','none')
hold on
plot(coast_lon,coast_lat,'Color',[grey grey grey])
colorbar
title('u')

figure('position',[151 236 1128 683])
% subplot(3,2,1)
contourf(lon,lat,tlm_20120410_1400z_u(:,:,plotlevel)','LineStyle','none')
hold on
plot(coast_lon,coast_lat,'Color',[grey grey grey])
colorbar
title('u')

figure('position',[151 236 1128 683])
% subplot(3,2,1)
contourf(lon,lat,tlm_20120410_1500z_u(:,:,plotlevel)','LineStyle','none')
hold on
plot(coast_lon,coast_lat,'Color',[grey grey grey])
colorbar
title('u')

figure('position',[151 236 1128 683])
% subplot(3,2,1)
contourf(lon,lat,tlm_20120410_1520z_u(:,:,plotlevel)','LineStyle','none')
hold on
plot(coast_lon,coast_lat,'Color',[grey grey grey])
colorbar
title('u')

figure('position',[151 236 1128 683])
% subplot(3,2,1)
contourf(lon,lat,tlm_20120410_1540z_u(:,:,plotlevel)','LineStyle','none')
hold on
plot(coast_lon,coast_lat,'Color',[grey grey grey])
colorbar
title('u')

figure('position',[151 236 1128 683])
% subplot(3,2,1)
contourf(lon,lat,tlm_20120410_1600z_u(:,:,plotlevel)','LineStyle','none')
hold on
plot(coast_lon,coast_lat,'Color',[grey grey grey])
colorbar
title('u')

figure('position',[151 236 1128 683])
% subplot(3,2,1)
contourf(lon,lat,tlm_20120410_1800z_u(:,:,plotlevel)','LineStyle','none')
hold on
plot(coast_lon,coast_lat,'Color',[grey grey grey])
colorbar
title('u')

% figure('position',[151 236 1128 683])
% % subplot(3,2,1)
% contourf(lon,lat,tlm_20120411_0000z_u(:,:,plotlevel)','LineStyle','none')
% hold on
% plot(coast_lon,coast_lat,'Color',[grey grey grey])
% colorbar
% title('u')
