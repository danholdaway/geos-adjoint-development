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
Jgradf_u = ncread('JgradfX.eta.nc4','U');
Jgradf_v = ncread('JgradfX.eta.nc4','V');
Jgradf_t = ncread('JgradfX.eta.nc4','TV');
Jgradf_q = ncread('JgradfX.eta.nc4','QV');
Jgradf_p = ncread('JgradfX.eta.nc4','DP');

adjoint24_u = ncread('fsens.adjoint.20120411_00z.nc4','u');
adjoint24_v = ncread('fsens.adjoint.20120411_00z.nc4','v');
adjoint24_t = ncread('fsens.adjoint.20120411_00z.nc4','tv');
adjoint24_q = ncread('fsens.adjoint.20120411_00z.nc4','sphu');
adjoint24_p = ncread('fsens.adjoint.20120411_00z.nc4','delp');

adjoint1_u = ncread('fsens.adjoint.20120410_0020z.nc4','u');
adjoint1_v = ncread('fsens.adjoint.20120410_0020z.nc4','v');
adjoint1_t = ncread('fsens.adjoint.20120410_0020z.nc4','tv');
adjoint1_q = ncread('fsens.adjoint.20120410_0020z.nc4','sphu');
adjoint1_p = ncread('fsens.adjoint.20120410_0020z.nc4','delp');

tlm24_u = ncread('fsens.tlm.20120410_12z.nc4','u');
tlm24_v = ncread('fsens.tlm.20120410_12z.nc4','v');
tlm24_t = ncread('fsens.tlm.20120410_12z.nc4','tv');
tlm24_q = ncread('fsens.tlm.20120410_12z.nc4','sphu');
tlm24_p = ncread('fsens.tlm.20120410_12z.nc4','delp');

tlm48_u = ncread('fsens.tlm.20120412_00z_reset.nc4','u');
tlm48_v = ncread('fsens.tlm.20120412_00z_reset.nc4','v');
tlm48_t = ncread('fsens.tlm.20120412_00z_reset.nc4','tv');
tlm48_q = ncread('fsens.tlm.20120412_00z_reset.nc4','sphu');
tlm48_p = ncread('fsens.tlm.20120412_00z_reset.nc4','delp');

tlm24_ai_u = ncread('fsens.tlm.20120411_00z_anainc.nc4','u');
tlm24_ai_v = ncread('fsens.tlm.20120411_00z_anainc.nc4','v');
tlm24_ai_t = ncread('fsens.tlm.20120411_00z_anainc.nc4','tv');
tlm24_ai_q = ncread('fsens.tlm.20120411_00z_anainc.nc4','sphu');
tlm24_ai_p = ncread('fsens.tlm.20120411_00z_anainc.nc4','delp');


cd /discover/nobackup/drholdaw/iau590a/asens/
adjoint24_asens_u = ncread('iau590a.fsens_txe.eta.20120409_15z+20120411_00z-20120410_00z.nc4','u');
adjoint24_asens_v = ncread('iau590a.fsens_txe.eta.20120409_15z+20120411_00z-20120410_00z.nc4','v');
adjoint24_asens_t = ncread('iau590a.fsens_txe.eta.20120409_15z+20120411_00z-20120410_00z.nc4','tv');
adjoint24_asens_q = ncread('iau590a.fsens_txe.eta.20120409_15z+20120411_00z-20120410_00z.nc4','sphu');
adjoint24_asens_p = ncread('iau590a.fsens_txe.eta.20120409_15z+20120411_00z-20120410_00z.nc4','delp');

cd /home/drholdaw/Lin_Moist_Physics/DAS_plotting/

plotlevel = 50;

% for plotlevel = 40:72;

% figure('position',[151 236 1128 683])
% subplot(3,2,1)
% contourf(lon,lat,Jgradf_u(:,:,plotlevel)','LineStyle','none')
% hold on
% plot(coast_lon,coast_lat,'Color',[grey grey grey])
% colorbar
% title('u')
% 
% subplot(3,2,2)
% contourf(lon,lat,Jgradf_v(:,:,plotlevel)','LineStyle','none')
% hold on
% plot(coast_lon,coast_lat,'Color',[grey grey grey])
% colorbar
% title('v')
% 
% subplot(3,2,3)
% contourf(lon,lat,Jgradf_t(:,:,plotlevel)','LineStyle','none')
% hold on
% plot(coast_lon,coast_lat,'Color',[grey grey grey])
% colorbar
% title('t')
% 
% subplot(3,2,4)
% contourf(lon,lat,Jgradf_q(:,:,plotlevel)','LineStyle','none')
% hold on
% plot(coast_lon,coast_lat,'Color',[grey grey grey])
% colorbar
% title('q')
% 
% subplot(3,2,5)
% contourf(lon,lat,Jgradf_p(:,:,plotlevel)','LineStyle','none')
% hold on
% plot(coast_lon,coast_lat,'Color',[grey grey grey])
% colorbar
% title('p')
% 
% % pause
% % close
% % 
% % end
% 
% % asd
% 
% figure('position',[151 236 1128 683])
% subplot(3,2,1)
% contourf(lon,lat,adjoint24_u(:,:,plotlevel)','LineStyle','none')
% hold on
% plot(coast_lon,coast_lat,'Color',[grey grey grey])
% colorbar
% title('u')
% 
% subplot(3,2,2)
% contourf(lon,lat,adjoint24_v(:,:,plotlevel)','LineStyle','none')
% hold on
% plot(coast_lon,coast_lat,'Color',[grey grey grey])
% colorbar
% title('v')
% 
% subplot(3,2,3)
% contourf(lon,lat,adjoint24_t(:,:,plotlevel)','LineStyle','none')
% hold on
% plot(coast_lon,coast_lat,'Color',[grey grey grey])
% colorbar
% title('t')
% 
% subplot(3,2,4)
% contourf(lon,lat,adjoint24_q(:,:,plotlevel)','LineStyle','none')
% hold on
% plot(coast_lon,coast_lat,'Color',[grey grey grey])
% colorbar
% title('q')
% 
% subplot(3,2,5)
% contourf(lon,lat,adjoint24_p(:,:,plotlevel)','LineStyle','none')
% hold on
% plot(coast_lon,coast_lat,'Color',[grey grey grey])
% colorbar
% title('p')
% 
% 
% figure('position',[151 236 1128 683])
% subplot(3,2,1)
% contourf(lon,lat,adjoint1_u(:,:,plotlevel)','LineStyle','none')
% hold on
% plot(coast_lon,coast_lat,'Color',[grey grey grey])
% colorbar
% title('u')
% 
% subplot(3,2,2)
% contourf(lon,lat,adjoint1_v(:,:,plotlevel)','LineStyle','none')
% hold on
% plot(coast_lon,coast_lat,'Color',[grey grey grey])
% colorbar
% title('v')
% 
% subplot(3,2,3)
% contourf(lon,lat,adjoint1_t(:,:,plotlevel)','LineStyle','none')
% hold on
% plot(coast_lon,coast_lat,'Color',[grey grey grey])
% colorbar
% title('t')
% 
% subplot(3,2,4)
% contourf(lon,lat,adjoint1_q(:,:,plotlevel)','LineStyle','none')
% hold on
% plot(coast_lon,coast_lat,'Color',[grey grey grey])
% colorbar
% title('q')
% 
% subplot(3,2,5)
% contourf(lon,lat,adjoint1_p(:,:,plotlevel)','LineStyle','none')
% hold on
% plot(coast_lon,coast_lat,'Color',[grey grey grey])
% colorbar
% title('p')


figure('position',[151 236 1128 683])
subplot(3,2,1)
contourf(lon,lat,tlm24_u(:,:,plotlevel)','LineStyle','none')
hold on
plot(coast_lon,coast_lat,'Color',[grey grey grey])
colorbar
title('u')

subplot(3,2,2)
contourf(lon,lat,tlm24_v(:,:,plotlevel)','LineStyle','none')
hold on
plot(coast_lon,coast_lat,'Color',[grey grey grey])
colorbar
title('v')

subplot(3,2,3)
contourf(lon,lat,tlm24_t(:,:,plotlevel)','LineStyle','none')
hold on
plot(coast_lon,coast_lat,'Color',[grey grey grey])
colorbar
title('t')

subplot(3,2,4)
contourf(lon,lat,tlm24_q(:,:,plotlevel)','LineStyle','none')
hold on
plot(coast_lon,coast_lat,'Color',[grey grey grey])
colorbar
title('q')

subplot(3,2,5)
contourf(lon,lat,tlm24_p(:,:,plotlevel)','LineStyle','none')
hold on
plot(coast_lon,coast_lat,'Color',[grey grey grey])
colorbar
title('p')


% figure('position',[151 236 1128 683])
% subplot(3,2,1)
% contourf(lon,lat,tlm24_ai_u(:,:,plotlevel)','LineStyle','none')
% hold on
% plot(coast_lon,coast_lat,'Color',[grey grey grey])
% colorbar
% title('u')
% 
% subplot(3,2,2)
% contourf(lon,lat,tlm24_ai_v(:,:,plotlevel)','LineStyle','none')
% hold on
% plot(coast_lon,coast_lat,'Color',[grey grey grey])
% colorbar
% title('v')
% 
% subplot(3,2,3)
% contourf(lon,lat,tlm24_ai_t(:,:,plotlevel)','LineStyle','none')
% hold on
% plot(coast_lon,coast_lat,'Color',[grey grey grey])
% colorbar
% title('t')
% 
% subplot(3,2,4)
% contourf(lon,lat,tlm24_ai_q(:,:,plotlevel)','LineStyle','none')
% hold on
% plot(coast_lon,coast_lat,'Color',[grey grey grey])
% colorbar
% title('q')
% 
% subplot(3,2,5)
% contourf(lon,lat,tlm24_ai_p(:,:,plotlevel)','LineStyle','none')
% hold on
% plot(coast_lon,coast_lat,'Color',[grey grey grey])
% colorbar
% title('p')
% 
% 
% figure('position',[151 236 1128 683])
% subplot(3,2,1)
% contourf(lon,lat,tlm48_u(:,:,plotlevel)','LineStyle','none')
% hold on
% plot(coast_lon,coast_lat,'Color',[grey grey grey])
% colorbar
% title('u')
% 
% subplot(3,2,2)
% contourf(lon,lat,tlm48_v(:,:,plotlevel)','LineStyle','none')
% hold on
% plot(coast_lon,coast_lat,'Color',[grey grey grey])
% colorbar
% title('v')
% 
% subplot(3,2,3)
% contourf(lon,lat,tlm48_t(:,:,plotlevel)','LineStyle','none')
% hold on
% plot(coast_lon,coast_lat,'Color',[grey grey grey])
% colorbar
% title('t')
% 
% subplot(3,2,4)
% contourf(lon,lat,tlm48_q(:,:,plotlevel)','LineStyle','none')
% hold on
% plot(coast_lon,coast_lat,'Color',[grey grey grey])
% colorbar
% title('q')
% 
% subplot(3,2,5)
% contourf(lon,lat,tlm48_p(:,:,plotlevel)','LineStyle','none')
% hold on
% plot(coast_lon,coast_lat,'Color',[grey grey grey])
% colorbar
% title('p')

