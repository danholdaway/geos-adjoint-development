close all
clear
clc

cd /discover/nobackup/drholdaw/tmp.31275/sens.20120318.000000/

lon = ncread('fvpert.eta.nc4','lon');
lat = ncread('fvpert.eta.nc4','lat');

U_a  = ncread('fsens.eta.moi.uon.20120317_06z.nc4','u');
V_a  = ncread('fsens.eta.moi.uon.20120317_06z.nc4','v');
TV_a = ncread('fsens.eta.moi.uon.20120317_06z.nc4','tv');
DP_a = ncread('fsens.eta.moi.uon.20120317_06z.nc4','delp');
QV_a = ncread('fsens.eta.moi.uon.20120317_06z.nc4','sphu');

U_b  = ncread('fsens.eta.moi.uoffnew.20120317_06z.nc4','u');
V_b  = ncread('fsens.eta.moi.uoffnew.20120317_06z.nc4','v');
TV_b = ncread('fsens.eta.moi.uoffnew.20120317_06z.nc4','tv');
DP_b = ncread('fsens.eta.moi.uoffnew.20120317_06z.nc4','delp');
QV_b = ncread('fsens.eta.moi.uoffnew.20120317_06z.nc4','sphu');

U_c  = ncread('fsens.eta.moi.alloff.20120317_06z.nc4','u');
V_c  = ncread('fsens.eta.moi.alloff.20120317_06z.nc4','v');
TV_c = ncread('fsens.eta.moi.alloff.20120317_06z.nc4','tv');
DP_c = ncread('fsens.eta.moi.alloff.20120317_06z.nc4','delp');
QV_c = ncread('fsens.eta.moi.alloff.20120317_06z.nc4','sphu');

cd /home/drholdaw/Lin_Moist_Physics/DAS_plotting/

plotlevel = 50;

Umax = max(max([U_a(:,:,plotlevel); U_b(:,:,plotlevel); U_c(:,:,plotlevel)]));
Umin = min(min([U_a(:,:,plotlevel); U_b(:,:,plotlevel); U_c(:,:,plotlevel)]));
Vmax = max(max([V_a(:,:,plotlevel); V_b(:,:,plotlevel); V_c(:,:,plotlevel)]));
Vmin = min(min([V_a(:,:,plotlevel); V_b(:,:,plotlevel); V_c(:,:,plotlevel)]));
TVmax = max(max([TV_a(:,:,plotlevel); TV_b(:,:,plotlevel); TV_c(:,:,plotlevel)]));
TVmin = min(min([TV_a(:,:,plotlevel); TV_b(:,:,plotlevel); TV_c(:,:,plotlevel)]));
DPmax = max(max([DP_a(:,:,plotlevel); DP_b(:,:,plotlevel); DP_c(:,:,plotlevel)]));
DPmin = min(min([DP_a(:,:,plotlevel); DP_b(:,:,plotlevel); DP_c(:,:,plotlevel)]));
QVmax = max(max([QV_a(:,:,plotlevel); QV_b(:,:,plotlevel); QV_c(:,:,plotlevel)]));
QVmin = min(min([QV_a(:,:,plotlevel); QV_b(:,:,plotlevel); QV_c(:,:,plotlevel)]));


U_e = U_a - U_b;
V_e = V_a - V_b;
TV_e = TV_a - TV_b;
DP_e = DP_a - DP_b;
QV_e = QV_a - QV_b;


scrsz = get(0,'ScreenSize');

    
figure('visible','on','Position',[0 scrsz(4) scrsz(3)/3 scrsz(4)/2])

subplot(2,3,1)
contour(lon,lat,U_a(:,:,plotlevel)')
caxis([Umin Umax])
colorbar
title('U \prime')

subplot(2,3,2)
contour(lon,lat,V_a(:,:,plotlevel)')
caxis([Vmin Vmax])
colorbar
title('V \prime')

subplot(2,3,3)
contour(lon,lat,TV_a(:,:,plotlevel)')
caxis([TVmin TVmax])
colorbar
title('\theta \prime')

subplot(2,3,4)
contour(lon,lat,DP_a(:,:,plotlevel)')
caxis([DPmin DPmax])
colorbar
title('P \prime')

subplot(2,3,5)
contour(lon,lat,QV_a(:,:,plotlevel)')
caxis([QVmin QVmax])
colorbar
title('Q_{tot} \prime')



figure('visible','on','Position',[scrsz(3)/2 scrsz(4) scrsz(3)/3 scrsz(4)/2])

subplot(2,3,1)
contour(lon,lat,U_b(:,:,plotlevel)')
caxis([Umin Umax])
colorbar
title('U \prime')

subplot(2,3,2)
contour(lon,lat,V_b(:,:,plotlevel)')
caxis([Vmin Vmax])
colorbar
title('V \prime')

subplot(2,3,3)
contour(lon,lat,TV_b(:,:,plotlevel)')
caxis([TVmin TVmax])
colorbar
title('\theta \prime')

subplot(2,3,4)
contour(lon,lat,DP_b(:,:,plotlevel)')
caxis([DPmin DPmax])
colorbar
title('P \prime')

subplot(2,3,5)
contour(lon,lat,QV_b(:,:,plotlevel)')
caxis([QVmin QVmax])
colorbar
title('Q_{tot} \prime')

figure('visible','on','Position',[scrsz(3)/2 scrsz(4) scrsz(3)/3 scrsz(4)/2])

subplot(2,3,1)
contour(lon,lat,U_c(:,:,plotlevel)')
caxis([Umin Umax])
colorbar
title('U \prime')

subplot(2,3,2)
contour(lon,lat,V_c(:,:,plotlevel)')
caxis([Vmin Vmax])
colorbar
title('V \prime')

subplot(2,3,3)
contour(lon,lat,TV_c(:,:,plotlevel)')
caxis([TVmin TVmax])
colorbar
title('\theta \prime')

subplot(2,3,4)
contour(lon,lat,DP_c(:,:,plotlevel)')
caxis([DPmin DPmax])
colorbar
title('P \prime')

subplot(2,3,5)
contour(lon,lat,QV_c(:,:,plotlevel)')
caxis([QVmin QVmax])
colorbar
title('Q_{tot} \prime')


figure('visible','on','Position',[scrsz(3)/2 scrsz(4) scrsz(3)/3 scrsz(4)/2])

subplot(2,3,1)
contour(lon,lat,U_e(:,:,plotlevel)')
% caxis([Umin Umax])
colorbar
title('U Error')

subplot(2,3,2)
contour(lon,lat,V_e(:,:,plotlevel)')
% caxis([Vmin Vmax])
colorbar
title('V Error')

subplot(2,3,3)
contour(lon,lat,TV_e(:,:,plotlevel)')
% caxis([TVmin TVmax])
colorbar
title('\theta Error')

subplot(2,3,4)
contour(lon,lat,DP_e(:,:,plotlevel)')
% caxis([DPmin DPmax])
colorbar
title('P Error')

subplot(2,3,5)
contour(lon,lat,QV_e(:,:,plotlevel)')
% caxis([QVmin QVmax])
colorbar
title('Q_{tot} Error')



