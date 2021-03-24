close all
clear
clc

load topo.mat;
topo_grenwich(1:180,181:360) = topo(1:180,1:180);
topo_grenwich(1:180,1:180) = topo(1:180,181:360);

% cd /discover/nobackup/drholdaw/iau590/asens/fsens_data/fsens_data_dryphys_storm/
cd /archive/u/drholdaw/fsens_data/iau590/fsens_data_dryphys_storm/


lon = ncread('iau590.fsens_txe.eta.20120320_21z+20120322_00z-20120321_00z.nc4','lon');
lat = ncread('iau590.fsens_txe.eta.20120320_21z+20120322_00z-20120321_00z.nc4','lat');

u_txe_dryp = ncread('iau590.fsens_txe.eta.20120320_21z+20120322_00z-20120321_00z.nc4','u');
v_txe_dryp = ncread('iau590.fsens_txe.eta.20120320_21z+20120322_00z-20120321_00z.nc4','v');
tv_txe_dryp = ncread('iau590.fsens_txe.eta.20120320_21z+20120322_00z-20120321_00z.nc4','tv');
q_txe_dryp = ncread('iau590.fsens_txe.eta.20120320_21z+20120322_00z-20120321_00z.nc4','sphu');

u_twe_dryp = ncread('iau590.fsens_twe.eta.20120320_21z+20120322_00z-20120321_00z.nc4','u');
v_twe_dryp = ncread('iau590.fsens_twe.eta.20120320_21z+20120322_00z-20120321_00z.nc4','v');
tv_twe_dryp = ncread('iau590.fsens_twe.eta.20120320_21z+20120322_00z-20120321_00z.nc4','tv');
q_twe_dryp = ncread('iau590.fsens_twe.eta.20120320_21z+20120322_00z-20120321_00z.nc4','sphu');

cd /archive/u/drholdaw/fsens_data/iau590/fsens_data_moiphys_storm/

u_txe_moip = ncread('iau590.fsens_txe.eta.20120320_21z+20120322_00z-20120321_00z.nc4','u');
v_txe_moip = ncread('iau590.fsens_txe.eta.20120320_21z+20120322_00z-20120321_00z.nc4','v');
tv_txe_moip = ncread('iau590.fsens_txe.eta.20120320_21z+20120322_00z-20120321_00z.nc4','tv');
q_txe_moip = ncread('iau590.fsens_txe.eta.20120320_21z+20120322_00z-20120321_00z.nc4','sphu');

u_twe_moip = ncread('iau590.fsens_twe.eta.20120320_21z+20120322_00z-20120321_00z.nc4','u');
v_twe_moip = ncread('iau590.fsens_twe.eta.20120320_21z+20120322_00z-20120321_00z.nc4','v');
tv_twe_moip = ncread('iau590.fsens_twe.eta.20120320_21z+20120322_00z-20120321_00z.nc4','tv');
q_twe_moip = ncread('iau590.fsens_twe.eta.20120320_21z+20120322_00z-20120321_00z.nc4','sphu');

cd /home/drholdaw/Lin_Moist_Physics/DAS_plotting/

boxhor = -37:-34;
boxver = 45:48;

plot_level = 54;

%Dry Physics
max_u_txe_dryp = max(max(abs(u_txe_dryp(:,:,plot_level))));
max_v_txe_dryp = max(max(abs(v_txe_dryp(:,:,plot_level))));
max_tv_txe_dryp = max(max(abs(tv_txe_dryp(:,:,plot_level))));
max_q_txe_dryp = max(max(abs(q_txe_dryp(:,:,plot_level))));

max_u_twe_dryp = max(max(abs(u_twe_dryp(:,:,plot_level))));
max_v_twe_dryp = max(max(abs(v_twe_dryp(:,:,plot_level))));
max_tv_twe_dryp = max(max(abs(tv_twe_dryp(:,:,plot_level))));
max_q_twe_dryp = max(max(abs(q_twe_dryp(:,:,plot_level))));

%Moist Physics
max_u_txe_moip = max(max(abs(u_txe_moip(:,:,plot_level))));
max_v_txe_moip = max(max(abs(v_txe_moip(:,:,plot_level))));
max_tv_txe_moip = max(max(abs(tv_txe_moip(:,:,plot_level))));
max_q_txe_moip = max(max(abs(q_txe_moip(:,:,plot_level))));

max_u_twe_moip = max(max(abs(u_twe_moip(:,:,plot_level))));
max_v_twe_moip = max(max(abs(v_twe_moip(:,:,plot_level))));
max_tv_twe_moip = max(max(abs(tv_twe_moip(:,:,plot_level))));
max_q_twe_moip = max(max(abs(q_twe_moip(:,:,plot_level))));

max_u_moip = max(max(abs([u_txe_moip(:,:,plot_level); u_twe_moip(:,:,plot_level)])));
max_v_moip = max(max(abs([v_txe_moip(:,:,plot_level); v_twe_moip(:,:,plot_level)])));
max_tv_moip = max(max(abs([tv_txe_moip(:,:,plot_level); tv_twe_moip(:,:,plot_level)])));
max_q_moip = max(max(abs([q_txe_moip(:,:,plot_level); q_twe_moip(:,:,plot_level)])));

max_u = max(max(abs([u_txe_dryp(:,:,plot_level); u_twe_dryp(:,:,plot_level); u_txe_moip(:,:,plot_level); u_twe_moip(:,:,plot_level)])));
max_v = max(max(abs([v_txe_dryp(:,:,plot_level); v_twe_dryp(:,:,plot_level); v_txe_moip(:,:,plot_level); v_twe_moip(:,:,plot_level)])));
max_tv = max(max(abs([tv_txe_dryp(:,:,plot_level); tv_twe_dryp(:,:,plot_level); tv_txe_moip(:,:,plot_level); tv_twe_moip(:,:,plot_level)])));
max_q = max(max(abs([q_txe_dryp(:,:,plot_level); q_twe_dryp(:,:,plot_level); q_txe_moip(:,:,plot_level); q_twe_moip(:,:,plot_level)])));


lon1 = -100;
lon2 = 0;

lat1 = 10;
lat2 = 80;

line_wid = 1;
fontsize = 14;

%%%%%%%%%%%%%%%%%%
%  Dry Physics   %
%%%%%%%%%%%%%%%%%%
figure
subplot(2,2,1)
[C,h] = contour(lon,lat,u_txe_dryp(:,:,plot_level)');
hold on; 
plot(boxhor,boxver(1)*ones(1,4),'k')
plot(boxhor,boxver(end)*ones(1,4),'k')
plot(boxhor(1)*ones(1,4),boxver,'k')
plot(boxhor(end)*ones(1,4),boxver,'k')
contour(-180:179,-89:90,topo_grenwich,[0 0],'k','LineWidth',1.0)

colorbar
xlim([lon1 lon2])
ylim([lat1 lat2])
caxis([-max_u max_u])
% caxis([-max_u_txe_dryp max_u_txe_dryp])
title('u','FontSize',fontsize,'FontName','TimesNewRoman')
step1 = get(h,'LevelStep');
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

pos = get(gca,'Position');
set(gca,'Position',[0.70*pos(1), pos(2), 1.4*pos(3), pos(4)])

subplot(2,2,2)
[C,h] = contour(lon,lat,v_txe_dryp(:,:,plot_level)');
hold on; 
plot(boxhor,boxver(1)*ones(1,4),'k')
plot(boxhor,boxver(end)*ones(1,4),'k')
plot(boxhor(1)*ones(1,4),boxver,'k')
plot(boxhor(end)*ones(1,4),boxver,'k')
contour(-180:179,-89:90,topo_grenwich,[0 0],'k','LineWidth',1.0)

colorbar
xlim([lon1 lon2])
ylim([lat1 lat2])
caxis([-max_v max_v])
% caxis([-max_v_txe_dryp max_v_txe_dryp])
title('v','FontSize',fontsize,'FontName','TimesNewRoman')
step2 = get(h,'LevelStep');
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

set(gca,'YtickLabel',[])

pos = get(gca,'Position');
set(gca,'Position',[0.95*pos(1), pos(2), 1.4*pos(3), pos(4)])


subplot(2,2,3)
[C,h] = contour(lon,lat,tv_txe_dryp(:,:,plot_level)');
hold on; 
plot(boxhor,boxver(1)*ones(1,4),'k')
plot(boxhor,boxver(end)*ones(1,4),'k')
plot(boxhor(1)*ones(1,4),boxver,'k')
plot(boxhor(end)*ones(1,4),boxver,'k')
contour(-180:179,-89:90,topo_grenwich,[0 0],'k','LineWidth',1.0)

colorbar
xlim([lon1 lon2])
ylim([lat1 lat2])
caxis([-max_tv max_tv])
% caxis([-max_tv_txe_dryp max_tv_txe_dryp])
title('T_v','FontSize',fontsize,'FontName','TimesNewRoman')
step3 = get(h,'LevelStep');
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

pos = get(gca,'Position');
set(gca,'Position',[0.70*pos(1), pos(2), 1.4*pos(3), pos(4)])

subplot(2,2,4)
[C,h] = contour(lon,lat,q_txe_dryp(:,:,plot_level)');
hold on; 
plot(boxhor,boxver(1)*ones(1,4),'k')
plot(boxhor,boxver(end)*ones(1,4),'k')
plot(boxhor(1)*ones(1,4),boxver,'k')
plot(boxhor(end)*ones(1,4),boxver,'k')
contour(-180:179,-89:90,topo_grenwich,[0 0],'k','LineWidth',1.0)

colorbar
xlim([lon1 lon2])
ylim([lat1 lat2])
% caxis([-max_q max_q])
caxis([-max_q_txe_dryp max_q_txe_dryp])
title('q','FontSize',fontsize,'FontName','TimesNewRoman')
step4 = get(h,'LevelStep');
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

set(gca,'YtickLabel',[])

pos = get(gca,'Position');
set(gca,'Position',[0.95*pos(1), pos(2), 1.4*pos(3), pos(4)])

% figure
% subplot(2,2,1)
% [C,h] = contour(lon,lat,u_twe_dryp(:,:,plot_level)','LevelStep',step1);
% hold on; 
% plot(boxhor,boxver(1)*ones(1,4),'k')
% plot(boxhor,boxver(end)*ones(1,4),'k')
% plot(boxhor(1)*ones(1,4),boxver,'k')
% plot(boxhor(end)*ones(1,4),boxver,'k')
% contour(-180:179,-89:90,topo_grenwich,[0 0],'k','LineWidth',1.0)
% 
% colorbar
% xlim([lon1 lon2])
% ylim([lat1 lat2])
% caxis([-max_u max_u])
% % caxis([-max_u_twe_dryp max_u_twe_dryp])
% title('u')
% 
% subplot(2,2,2)
% [C,h] = contour(lon,lat,v_twe_dryp(:,:,plot_level)','LevelStep',step2);
% hold on; 
% plot(boxhor,boxver(1)*ones(1,4),'k')
% plot(boxhor,boxver(end)*ones(1,4),'k')
% plot(boxhor(1)*ones(1,4),boxver,'k')
% plot(boxhor(end)*ones(1,4),boxver,'k')
% contour(-180:179,-89:90,topo_grenwich,[0 0],'k','LineWidth',1.0)
% 
% colorbar
% xlim([lon1 lon2])
% ylim([lat1 lat2])
% caxis([-max_v max_v])
% % caxis([-max_v_twe_dryp max_v_twe_dryp])
% title('v')
% 
% subplot(2,2,3)
% [C,h] = contour(lon,lat,tv_twe_dryp(:,:,plot_level)','LevelStep',step3);
% hold on; 
% plot(boxhor,boxver(1)*ones(1,4),'k')
% plot(boxhor,boxver(end)*ones(1,4),'k')
% plot(boxhor(1)*ones(1,4),boxver,'k')
% plot(boxhor(end)*ones(1,4),boxver,'k')
% contour(-180:179,-89:90,topo_grenwich,[0 0],'k','LineWidth',1.0)
% 
% colorbar
% xlim([lon1 lon2])
% ylim([lat1 lat2])
% caxis([-max_tv max_tv])
% % caxis([-max_tv_twe_dryp max_tv_twe_dryp])
% title('T_v')
% 
% subplot(2,2,4)
% [C,h] = contour(lon,lat,q_twe_dryp(:,:,plot_level)');
% hold on; 
% plot(boxhor,boxver(1)*ones(1,4),'k')
% plot(boxhor,boxver(end)*ones(1,4),'k')
% plot(boxhor(1)*ones(1,4),boxver,'k')
% plot(boxhor(end)*ones(1,4),boxver,'k')
% contour(-180:179,-89:90,topo_grenwich,[0 0],'k','LineWidth',1.0)
% 
% colorbar
% xlim([lon1 lon2])
% ylim([lat1 lat2])
% % caxis([-max_q max_q])
% caxis([-max_q_twe_dryp max_q_twe_dryp])
% title('q')


%%%%%%%%%%%%%%%%%%%%%
%   Moist Physics   %
%%%%%%%%%%%%%%%%%%%%%
figure
subplot(2,2,1)
[C,h] = contour(lon,lat,u_txe_moip(:,:,plot_level)','LevelStep',step1);
hold on; 
plot(boxhor,boxver(1)*ones(1,4),'k')
plot(boxhor,boxver(end)*ones(1,4),'k')
plot(boxhor(1)*ones(1,4),boxver,'k')
plot(boxhor(end)*ones(1,4),boxver,'k')
contour(-180:179,-89:90,topo_grenwich,[0 0],'k','LineWidth',1.0)

colorbar
xlim([lon1 lon2])
ylim([lat1 lat2])
caxis([-max_u max_u])
% caxis([-max_u_moip max_u_moip])
title('u','FontSize',fontsize,'FontName','TimesNewRoman')
% step1 = get(h,'LevelStep');
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

pos = get(gca,'Position');
set(gca,'Position',[0.70*pos(1), pos(2), 1.4*pos(3), pos(4)])

subplot(2,2,2)
[C,h] = contour(lon,lat,v_txe_moip(:,:,plot_level)','LevelStep',step2);
hold on; 
plot(boxhor,boxver(1)*ones(1,4),'k')
plot(boxhor,boxver(end)*ones(1,4),'k')
plot(boxhor(1)*ones(1,4),boxver,'k')
plot(boxhor(end)*ones(1,4),boxver,'k')
contour(-180:179,-89:90,topo_grenwich,[0 0],'k','LineWidth',1.0)

colorbar
xlim([lon1 lon2])
ylim([lat1 lat2])
caxis([-max_v max_v])
% caxis([-max_v_moip max_v_moip])
title('v','FontSize',fontsize,'FontName','TimesNewRoman')
% step2 = get(h,'LevelStep');
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

set(gca,'YtickLabel',[])

pos = get(gca,'Position');
set(gca,'Position',[0.95*pos(1), pos(2), 1.4*pos(3), pos(4)])

subplot(2,2,3)
[C,h] = contour(lon,lat,tv_txe_moip(:,:,plot_level)','LevelStep',step3);
hold on; 
plot(boxhor,boxver(1)*ones(1,4),'k')
plot(boxhor,boxver(end)*ones(1,4),'k')
plot(boxhor(1)*ones(1,4),boxver,'k')
plot(boxhor(end)*ones(1,4),boxver,'k')
contour(-180:179,-89:90,topo_grenwich,[0 0],'k','LineWidth',1.0)

colorbar
xlim([lon1 lon2])
ylim([lat1 lat2])
caxis([-max_tv max_tv])
% caxis([-max_tv_moip max_tv_moip])
title('T_v','FontSize',fontsize,'FontName','TimesNewRoman')
% step3 = get(h,'LevelStep');
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

pos = get(gca,'Position');
set(gca,'Position',[0.70*pos(1), pos(2), 1.4*pos(3), pos(4)])

subplot(2,2,4)
[C,h] = contour(lon,lat,q_txe_moip(:,:,plot_level)');
hold on; 
plot(boxhor,boxver(1)*ones(1,4),'k')
plot(boxhor,boxver(end)*ones(1,4),'k')
plot(boxhor(1)*ones(1,4),boxver,'k')
plot(boxhor(end)*ones(1,4),boxver,'k')
contour(-180:179,-89:90,topo_grenwich,[0 0],'k','LineWidth',1.0)

colorbar
xlim([lon1 lon2])
ylim([lat1 lat2])
% caxis([-max_q max_q])
caxis([-max_q_moip max_q_moip])
title('q','FontSize',fontsize,'FontName','TimesNewRoman')
step4 = get(h,'LevelStep');
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

set(gca,'YtickLabel',[])

pos = get(gca,'Position');
set(gca,'Position',[0.95*pos(1), pos(2), 1.4*pos(3), pos(4)])

% figure
% subplot(2,2,1)
% [C,h] = contour(lon,lat,u_twe_moip(:,:,plot_level)','LevelStep',step1);
% hold on; 
% plot(boxhor,boxver(1)*ones(1,4),'k')
% plot(boxhor,boxver(end)*ones(1,4),'k')
% plot(boxhor(1)*ones(1,4),boxver,'k')
% plot(boxhor(end)*ones(1,4),boxver,'k')
% contour(-180:179,-89:90,topo_grenwich,[0 0],'k','LineWidth',1.0)
% 
% colorbar
% xlim([lon1 lon2])
% ylim([lat1 lat2])
% caxis([-max_u max_u])
% % caxis([-max_u_moip max_u_moip])
% title('u')
% 
% subplot(2,2,2)
% [C,h] = contour(lon,lat,v_twe_moip(:,:,plot_level)','LevelStep',step2);
% hold on; 
% plot(boxhor,boxver(1)*ones(1,4),'k')
% plot(boxhor,boxver(end)*ones(1,4),'k')
% plot(boxhor(1)*ones(1,4),boxver,'k')
% plot(boxhor(end)*ones(1,4),boxver,'k')
% contour(-180:179,-89:90,topo_grenwich,[0 0],'k','LineWidth',1.0)
% 
% colorbar
% xlim([lon1 lon2])
% ylim([lat1 lat2])
% caxis([-max_v max_v])
% % caxis([-max_v_moip max_v_moip])
% title('v')
% 
% subplot(2,2,3)
% [C,h] = contour(lon,lat,tv_twe_moip(:,:,plot_level)','LevelStep',step3);
% hold on; 
% plot(boxhor,boxver(1)*ones(1,4),'k')
% plot(boxhor,boxver(end)*ones(1,4),'k')
% plot(boxhor(1)*ones(1,4),boxver,'k')
% plot(boxhor(end)*ones(1,4),boxver,'k')
% contour(-180:179,-89:90,topo_grenwich,[0 0],'k','LineWidth',1.0)
% 
% colorbar
% xlim([lon1 lon2])
% ylim([lat1 lat2])
% caxis([-max_tv max_tv])
% % caxis([-max_tv_moip max_tv_moip])
% title('T_v')
% 
% subplot(2,2,4)
% [C,h] = contour(lon,lat,q_twe_moip(:,:,plot_level)','LevelStep',step4);
% hold on; 
% plot(boxhor,boxver(1)*ones(1,4),'k')
% plot(boxhor,boxver(end)*ones(1,4),'k')
% plot(boxhor(1)*ones(1,4),boxver,'k')
% plot(boxhor(end)*ones(1,4),boxver,'k')
% contour(-180:179,-89:90,topo_grenwich,[0 0],'k','LineWidth',1.0)
% 
% colorbar
% xlim([lon1 lon2])
% ylim([lat1 lat2])
% % caxis([-max_q max_q])
% caxis([-max_q_moip max_q_moip])
% title('q')


