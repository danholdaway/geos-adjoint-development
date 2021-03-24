close all
clear
clc

load coast
coast_lat = lat; clear lat
coast_lon = long; clear long

cd /discover/nobackup/drholdaw/ExperimentData/Journal_Articles/moist_dev_mwr/

lon = ncread('fsens.eta.moi.nofilt.20120317_06z.nc4','lon');
lat = ncread('fsens.eta.moi.nofilt.20120317_06z.nc4','lat');

U_moi  = ncread('fsens.eta.moi.nofilt.20120317_06z.nc4','u');
V_moi  = ncread('fsens.eta.moi.nofilt.20120317_06z.nc4','v');
TV_moi = ncread('fsens.eta.moi.nofilt.20120317_06z.nc4','tv');
DP_moi = ncread('fsens.eta.moi.nofilt.20120317_06z.nc4','delp');
QV_moi = ncread('fsens.eta.moi.nofilt.20120317_06z.nc4','sphu');
QL_moi = ncread('fsens.eta.moi.nofilt.20120317_06z.nc4','qltot');
QI_moi = ncread('fsens.eta.moi.nofilt.20120317_06z.nc4','qitot');
O3_moi = ncread('fsens.eta.moi.nofilt.20120317_06z.nc4','ozone');

cd /home/drholdaw/Lin_Moist_Physics/Moist_Dev_MWR/

QV_moi_log = log(abs(QV_moi));

plotlevel = 50;
line_wid_cont = 0.8;
line_wid_det = 0.6;
fontsize = 8;
grey = 0.75;

maxabsQ = max(max(abs(QV_moi(:,:,plotlevel))));

%Make center of colormap white
cmap = colormap;
cmap(32,:) = [1 1 1];
cmap(33,:) = [1 1 1];

[C,h] = contour(lon,lat,QV_moi(:,:,plotlevel)','LineWidth',line_wid_cont);
caxis([-maxabsQ maxabsQ])
colormap(cmap)
colbar = colorbar;
title(colbar,'(Jkg^{-1}kg^{-1}kg)','FontSize',fontsize,'FontName','TimesNewRoman') 
title('\partial J/\partial q\prime','FontSize',fontsize,'FontName','TimesNewRoman')
hold on; load topo.mat;
% xlabel('Longitude (^\circ)','FontSize',fontsize,'FontName','TimesNewRoman')
% ylabel('Latitude','FontSize',fontsize,'FontName','TimesNewRoman')
plot(coast_lon,coast_lat,'Color',[grey grey grey])
axis equal; box on; xlim([-180 180]); ylim([-90 90])
set(gca,'XTick',[-180 -120 -60 0 60 120 180],'Ytick',[-90 -60 -30 0 30 60 90]);
set(gca,'XTickLabel',{'180W', '120W', '60W', '0', '60E', '120E', '180E'}...
       ,'YtickLabel',{'90S' '60S' '30S' '0' '30N' '60N' '90N'});
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

levstep = get(h,'LevelStep');
set(h,'LevelStep',levstep/2);

% levlist = get(h,'LevelList');
% indzero = find(levlist==0);
% levlist(indzero) = [];
% set(h,'LevelList',levlist)

pos = get(gcf,'Position');
set(gcf,'Position',[pos(1) pos(2) pos(3) 0.6*pos(4)])


saveas(gcf,'nofiltering.eps', 'psc2')


% figure
% contourf(lon,lat,QV_moi_log(:,:,plotlevel)','LineStyle','none')
% colormap(cmap)
% colbar = colorbar;
% hold on
% plot(coast_lon,coast_lat,'Color',[grey grey grey])
% levstep = get(h,'LevelStep');
% set(h,'LevelStep',levstep/2);