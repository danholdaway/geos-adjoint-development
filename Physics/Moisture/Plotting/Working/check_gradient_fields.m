close all
clear
clc

load coast
coast_lat = lat; clear lat
coast_lon = long; clear long

cd /discover/nobackup/drholdaw/ExperimentData/GradientFields/dryphy_patch17/

u_dry = ncread('iau590a_mp.fsens_txe.eta.20120316_21z+20120318_00z-20120317_00z.nc4','u');
v_dry = ncread('iau590a_mp.fsens_txe.eta.20120316_21z+20120318_00z-20120317_00z.nc4','v');
t_dry = ncread('iau590a_mp.fsens_txe.eta.20120316_21z+20120318_00z-20120317_00z.nc4','tv');
q_dry = ncread('iau590a_mp.fsens_txe.eta.20120316_21z+20120318_00z-20120317_00z.nc4','sphu');
p_dry = ncread('iau590a_mp.fsens_txe.eta.20120316_21z+20120318_00z-20120317_00z.nc4','delp');




% cd /discover/nobackup/drholdaw/ExperimentData/GradientFields/dryphys/
% cd /discover/nobackup/drholdaw/ExperimentData/GradientFields/filter_2_patch17/
% cd /discover/nobackup/drholdaw/ExperimentData/GradientFields/filter_5_patch17/
cd /discover/nobackup/drholdaw/ExperimentData/GradientFields/filter_7.5_patch17/
% cd /discover/nobackup/drholdaw/ExperimentData/GradientFields/filter_10_patch17/

lon = ncread('iau590a_mp.fsens_txe.eta.20120316_21z+20120318_00z-20120317_00z.nc4','lon');
lat = ncread('iau590a_mp.fsens_txe.eta.20120316_21z+20120318_00z-20120317_00z.nc4','lat');

u_txe = ncread('iau590a_mp.fsens_txe.eta.20120316_21z+20120318_00z-20120317_00z.nc4','u');
v_txe = ncread('iau590a_mp.fsens_txe.eta.20120316_21z+20120318_00z-20120317_00z.nc4','v');
t_txe = ncread('iau590a_mp.fsens_txe.eta.20120316_21z+20120318_00z-20120317_00z.nc4','tv');
q_txe = ncread('iau590a_mp.fsens_txe.eta.20120316_21z+20120318_00z-20120317_00z.nc4','sphu');
p_txe = ncread('iau590a_mp.fsens_txe.eta.20120316_21z+20120318_00z-20120317_00z.nc4','delp');

% cd /discover/nobackup/drholdaw/ExperimentData/GradientFields/dryphy_patch17/
% cd /discover/nobackup/drholdaw/ExperimentData/GradientFields/filter_5_patch17/
% cd /discover/nobackup/drholdaw/ExperimentData/GradientFields/filter_7.5_patch17/
cd /discover/nobackup/drholdaw/ExperimentData/GradientFields/filter_10_patch17/

u_twe = ncread('iau590a_mp.fsens_txe.eta.20120316_21z+20120318_00z-20120317_00z.nc4','u');
v_twe = ncread('iau590a_mp.fsens_txe.eta.20120316_21z+20120318_00z-20120317_00z.nc4','v');
t_twe = ncread('iau590a_mp.fsens_txe.eta.20120316_21z+20120318_00z-20120317_00z.nc4','tv');
q_twe = ncread('iau590a_mp.fsens_txe.eta.20120316_21z+20120318_00z-20120317_00z.nc4','sphu');
p_twe = ncread('iau590a_mp.fsens_txe.eta.20120316_21z+20120318_00z-20120317_00z.nc4','delp');

cd /home/drholdaw/Lin_Moist_Physics/DAS_plotting/

plotlevel = 63;
grey = 0.5;

up_dry = u_dry(:,:,plotlevel);
vp_dry = v_dry(:,:,plotlevel);
tp_dry = t_dry(:,:,plotlevel);
qp_dry = q_dry(:,:,plotlevel);
pp_dry = p_dry(:,:,plotlevel);

up_txe = u_txe(:,:,plotlevel);
vp_txe = v_txe(:,:,plotlevel);
tp_txe = t_txe(:,:,plotlevel);
qp_txe = q_txe(:,:,plotlevel);
pp_txe = p_txe(:,:,plotlevel);

up_twe = u_twe(:,:,plotlevel);
vp_twe = v_twe(:,:,plotlevel);
tp_twe = t_twe(:,:,plotlevel);
qp_twe = q_twe(:,:,plotlevel);
pp_twe = p_twe(:,:,plotlevel);


figure('position',[4 30 1152 889])
subplot(3,2,1)
[C,h] = contourf(lon,lat,up_dry','LineStyle','none');
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])
colorbar
caxis([-max(abs(up_dry(:))) max(abs(up_dry(:)))])
title('u')

subplot(3,2,2)
[C,h] = contourf(lon,lat,vp_dry','LineStyle','none');
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])
colorbar
caxis([-max(abs(vp_dry(:))) max(abs(vp_dry(:)))])
title('v')

subplot(3,2,3)
[C,h] = contourf(lon,lat,tp_dry','LineStyle','none');
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])
colorbar
title('t')
caxis([-max(abs(tp_dry(:))) max(abs(tp_dry(:)))])

subplot(3,2,4)
[C,h] = contourf(lon,lat,qp_dry','LineStyle','none');
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])
colorbar
caxis([-max(abs(qp_dry(:))) max(abs(qp_dry(:)))])
title('q')

subplot(3,2,5)
[C,h] = contourf(lon,lat,pp_dry','LineStyle','none');
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])
colorbar
caxis([-max(abs(pp_dry(:))) max(abs(pp_dry(:)))])
title('p')


figure('position',[707 30 1152 889])
subplot(3,2,1)
[C,h] = contourf(lon,lat,up_txe','LineStyle','none');
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])
colorbar
caxis([-max(abs(up_txe(:))) max(abs(up_txe(:)))])
title('u')

subplot(3,2,2)
[C,h] = contourf(lon,lat,vp_txe','LineStyle','none');
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])
colorbar
caxis([-max(abs(vp_txe(:))) max(abs(vp_txe(:)))])
title('v')

subplot(3,2,3)
[C,h] = contourf(lon,lat,tp_txe','LineStyle','none');
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])
colorbar
caxis([-max(abs(tp_txe(:))) max(abs(tp_txe(:)))])
title('t')

subplot(3,2,4)
[C,h] = contourf(lon,lat,qp_txe','LineStyle','none');
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])
colorbar
caxis([-max(abs(qp_txe(:))) max(abs(qp_txe(:)))])
title('q')

subplot(3,2,5)
[C,h] = contourf(lon,lat,pp_txe','LineStyle','none');
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])
colorbar
caxis([-max(abs(pp_txe(:))) max(abs(pp_txe(:)))])
title('p')


figure('position',[1407 47 1152 889])
subplot(3,2,1)
[C,h] = contourf(lon,lat,up_twe','LineStyle','none');
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])
colorbar
caxis([-max(abs(up_twe(:))) max(abs(up_twe(:)))])
title('u')

subplot(3,2,2)
[C,h] = contourf(lon,lat,vp_twe','LineStyle','none');
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])
colorbar
caxis([-max(abs(vp_twe(:))) max(abs(vp_twe(:)))])
title('v')

subplot(3,2,3)
[C,h] = contourf(lon,lat,tp_twe','LineStyle','none');
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])
colorbar
caxis([-max(abs(tp_twe(:))) max(abs(tp_twe(:)))])
title('t')

subplot(3,2,4)
[C,h] = contourf(lon,lat,qp_twe','LineStyle','none');
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])
colorbar
caxis([-max(abs(qp_twe(:))) max(abs(qp_twe(:)))])
title('q')

subplot(3,2,5)
[C,h] = contourf(lon,lat,pp_twe','LineStyle','none');
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])
colorbar
caxis([-max(abs(pp_twe(:))) max(abs(pp_twe(:)))])
title('p')