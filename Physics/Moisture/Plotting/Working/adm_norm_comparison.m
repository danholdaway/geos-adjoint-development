close all
clear
clc

%LOAD DRY NORM RSULTS
cd /discover/nobackup/drholdaw/tmp.31275/sens.20120318.000000/

lon = ncread('JgradfX.eta.nc4','lon');
lat = ncread('JgradfX.eta.nc4','lat');

U_dry_normd  = ncread('fsens.eta.dry.normd.20120317_06z.nc4','u');
V_dry_normd  = ncread('fsens.eta.dry.normd.20120317_06z.nc4','v');
TV_dry_normd = ncread('fsens.eta.dry.normd.20120317_06z.nc4','tv');
DP_dry_normd = ncread('fsens.eta.dry.normd.20120317_06z.nc4','delp');
QV_dry_normd = ncread('fsens.eta.dry.normd.20120317_06z.nc4','sphu');
QL_dry_normd = ncread('fsens.eta.dry.normd.20120317_06z.nc4','qltot');
QI_dry_normd = ncread('fsens.eta.dry.normd.20120317_06z.nc4','qitot');
O3_dry_normd = ncread('fsens.eta.dry.normd.20120317_06z.nc4','ozone');

U_moi_normd  = ncread('fsens.eta.moi.normd.20120317_06z.nc4','u');
V_moi_normd  = ncread('fsens.eta.moi.normd.20120317_06z.nc4','v');
TV_moi_normd = ncread('fsens.eta.moi.normd.20120317_06z.nc4','tv');
DP_moi_normd = ncread('fsens.eta.moi.normd.20120317_06z.nc4','delp');
QV_moi_normd = ncread('fsens.eta.moi.normd.20120317_06z.nc4','sphu');
QL_moi_normd = ncread('fsens.eta.moi.normd.20120317_06z.nc4','qltot');
QI_moi_normd = ncread('fsens.eta.moi.normd.20120317_06z.nc4','qitot');
O3_moi_normd = ncread('fsens.eta.moi.normd.20120317_06z.nc4','ozone');

%LOAD MOIST NORM RESULTS
% cd /discover/nobackup/drholdaw/tmp.21867/sens.20120318.000000/
% 
% U_dry_normm  = ncread('fsens.eta.dry.normm.20120317_06z.nc4','u');
% V_dry_normm  = ncread('fsens.eta.dry.normm.20120317_06z.nc4','v');
% TV_dry_normm = ncread('fsens.eta.dry.normm.20120317_06z.nc4','tv');
% DP_dry_normm = ncread('fsens.eta.dry.normm.20120317_06z.nc4','delp');
% QV_dry_normm = ncread('fsens.eta.dry.normm.20120317_06z.nc4','sphu');
% QL_dry_normm = ncread('fsens.eta.dry.normm.20120317_06z.nc4','qltot');
% QI_dry_normm = ncread('fsens.eta.dry.normm.20120317_06z.nc4','qitot');
% O3_dry_normm = ncread('fsens.eta.dry.normm.20120317_06z.nc4','ozone');
% 
% U_moi_normm  = ncread('fsens.eta.moi.normm.20120317_06z.nc4','u');
% V_moi_normm  = ncread('fsens.eta.moi.normm.20120317_06z.nc4','v');
% TV_moi_normm = ncread('fsens.eta.moi.normm.20120317_06z.nc4','tv');
% DP_moi_normm = ncread('fsens.eta.moi.normm.20120317_06z.nc4','delp');
% QV_moi_normm = ncread('fsens.eta.moi.normm.20120317_06z.nc4','sphu');
% QL_moi_normm = ncread('fsens.eta.moi.normm.20120317_06z.nc4','qltot');
% QI_moi_normm = ncread('fsens.eta.moi.normm.20120317_06z.nc4','qitot');
% O3_moi_normm = ncread('fsens.eta.moi.normm.20120317_06z.nc4','ozone');


plotlevel = 50;

scrsz = get(0,'ScreenSize');
figure('visible','on','Position',[0 scrsz(4) scrsz(3)/3 scrsz(4)/2])

subplot(4,2,1)
contourf(lon,lat,U_dry_normd(:,:,plotlevel)')
% caxis([Umin Umax])
colorbar
title('U \prime')

subplot(4,2,2)
contourf(lon,lat,V_dry_normd(:,:,plotlevel)')
% caxis([Vmin Vmax])
colorbar
title('V \prime')

subplot(4,2,3)
contourf(lon,lat,TV_dry_normd(:,:,plotlevel)')
% caxis([TVmin TVmax])
colorbar
title('\theta \prime')

subplot(4,2,4)
contourf(lon,lat,DP_dry_normd(:,:,plotlevel)')
% caxis([DPmin DPmax])
colorbar
title('P \prime')

subplot(4,2,5)
contourf(lon,lat,QV_dry_normd(:,:,plotlevel)')
% caxis([QVmin QVmax])
colorbar
title('Q_{tot} \prime')

subplot(4,2,6)
contourf(lon,lat,QL_dry_normd(:,:,plotlevel)')
% caxis([QLmin QLmax])
colorbar
title('Q_L \prime')

subplot(4,2,7)
contourf(lon,lat,QI_dry_normd(:,:,plotlevel)')
% caxis([QImin QImax])
colorbar
title('Q_I \prime')

subplot(4,2,8)
contourf(lon,lat,O3_dry_normd(:,:,plotlevel)')
% caxis([O3min O3max])
colorbar
title('O^3 \prime')


figure('visible','on','Position',[scrsz(3)/2 scrsz(4) scrsz(3)/3 scrsz(4)/2])
subplot(4,2,1)
contourf(lon,lat,U_moi_normd(:,:,plotlevel)')
% caxis([Umin Umax])
colorbar
title('U \prime')

subplot(4,2,2)
contourf(lon,lat,V_moi_normd(:,:,plotlevel)')
% caxis([Vmin Vmax])
colorbar
title('V \prime')

subplot(4,2,3)
contourf(lon,lat,TV_moi_normd(:,:,plotlevel)')
% caxis([TVmin TVmax])
colorbar
title('\theta \prime')

subplot(4,2,4)
contourf(lon,lat,DP_moi_normd(:,:,plotlevel)')
% caxis([DPmin DPmax])
colorbar
title('P \prime')

subplot(4,2,5)
contourf(lon,lat,QV_moi_normd(:,:,plotlevel)')
% caxis([QVmin QVmax])
colorbar
title('Q_{tot} \prime')

subplot(4,2,6)
contourf(lon,lat,QL_moi_normd(:,:,plotlevel)')
% caxis([QLmin QLmax])
colorbar
title('Q_L \prime')

subplot(4,2,7)
contourf(lon,lat,QI_moi_normd(:,:,plotlevel)')
% caxis([QImin QImax])
colorbar
title('Q_I \prime')

subplot(4,2,8)
contourf(lon,lat,O3_moi_normd(:,:,plotlevel)')
% caxis([O3min O3max])
colorbar
title('O^3 \prime')






% scrsz = get(0,'ScreenSize');
% figure('visible','on','Position',[0 scrsz(4) scrsz(3)/3 scrsz(4)/2])
% 
% subplot(4,2,1)
% contourf(lon,lat,U_dry_normm(:,:,plotlevel)')
% % caxis([Umin Umax])
% colorbar
% title('U \prime')
% 
% subplot(4,2,2)
% contourf(lon,lat,V_dry_normm(:,:,plotlevel)')
% % caxis([Vmin Vmax])
% colorbar
% title('V \prime')
% 
% subplot(4,2,3)
% contourf(lon,lat,TV_dry_normm(:,:,plotlevel)')
% % caxis([TVmin TVmax])
% colorbar
% title('\theta \prime')
% 
% subplot(4,2,4)
% contourf(lon,lat,DP_dry_normm(:,:,plotlevel)')
% % caxis([DPmin DPmax])
% colorbar
% title('P \prime')
% 
% subplot(4,2,5)
% contourf(lon,lat,QV_dry_normm(:,:,plotlevel)')
% % caxis([QVmin QVmax])
% colorbar
% title('Q_{tot} \prime')
% 
% subplot(4,2,6)
% contourf(lon,lat,QL_dry_normm(:,:,plotlevel)')
% % caxis([QLmin QLmax])
% colorbar
% title('Q_L \prime')
% 
% subplot(4,2,7)
% contourf(lon,lat,QI_dry_normm(:,:,plotlevel)')
% % caxis([QImin QImax])
% colorbar
% title('Q_I \prime')
% 
% subplot(4,2,8)
% contourf(lon,lat,O3_dry_normm(:,:,plotlevel)')
% % caxis([O3min O3max])
% colorbar
% title('O^3 \prime')
% 
% 
% figure('visible','on','Position',[scrsz(3)/2 scrsz(4) scrsz(3)/3 scrsz(4)/2])
% subplot(4,2,1)
% contourf(lon,lat,U_moi_normm(:,:,plotlevel)')
% % caxis([Umin Umax])
% colorbar
% title('U \prime')
% 
% subplot(4,2,2)
% contourf(lon,lat,V_moi_normm(:,:,plotlevel)')
% % caxis([Vmin Vmax])
% colorbar
% title('V \prime')
% 
% subplot(4,2,3)
% contourf(lon,lat,TV_moi_normm(:,:,plotlevel)')
% % caxis([TVmin TVmax])
% colorbar
% title('\theta \prime')
% 
% subplot(4,2,4)
% contourf(lon,lat,DP_moi_normm(:,:,plotlevel)')
% % caxis([DPmin DPmax])
% colorbar
% title('P \prime')
% 
% subplot(4,2,5)
% contourf(lon,lat,QV_moi_normm(:,:,plotlevel)')
% % caxis([QVmin QVmax])
% colorbar
% title('Q_{tot} \prime')
% 
% subplot(4,2,6)
% contourf(lon,lat,QL_moi_normm(:,:,plotlevel)')
% % caxis([QLmin QLmax])
% colorbar
% title('Q_L \prime')
% 
% subplot(4,2,7)
% contourf(lon,lat,QI_moi_normm(:,:,plotlevel)')
% % caxis([QImin QImax])
% colorbar
% title('Q_I \prime')
% 
% subplot(4,2,8)
% contourf(lon,lat,O3_moi_normm(:,:,plotlevel)')
% % caxis([O3min O3max])
% colorbar
% title('O^3 \prime')
