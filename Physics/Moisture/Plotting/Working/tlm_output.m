close all
clear
clc

cd /discover/nobackup/drholdaw/tmp.31275/sens.20120318.000000/

lon = ncread('fvpert.eta.dry.20120317_06z.nc4','lon');
lat = ncread('fvpert.eta.dry.20120317_06z.nc4','lat');

U_dry  = ncread('fvpert.eta.dry.20120317_06z.nc4','u');
V_dry  = ncread('fvpert.eta.dry.20120317_06z.nc4','v');
TV_dry = ncread('fvpert.eta.dry.20120317_06z.nc4','tv');
DP_dry = ncread('fvpert.eta.dry.20120317_06z.nc4','delp');
QV_dry = ncread('fvpert.eta.dry.20120317_06z.nc4','sphu');
% QL_dry = ncread('fvpert.eta.dry.20120317_06z.nc4','QL');
% QI_dry = ncread('fvpert.eta.dry.20120317_06z.nc4','QI');
% O3_dry = ncread('fvpert.eta.dry.20120317_06z.nc4','O3');

U_moi  = ncread('fvpert.eta.moi.20120317_06z.nc4','u');
V_moi  = ncread('fvpert.eta.moi.20120317_06z.nc4','v');
TV_moi = ncread('fvpert.eta.moi.20120317_06z.nc4','tv');
DP_moi = ncread('fvpert.eta.moi.20120317_06z.nc4','delp');
QV_moi = ncread('fvpert.eta.moi.20120317_06z.nc4','sphu');
% QL_moi = ncread('fvpert.eta.moi.20120317_06z.nc4','QL');
% QI_moi = ncread('fvpert.eta.moi.20120317_06z.nc4','QI');
% O3_moi = ncread('fvpert.eta.moi.20120317_06z.nc4','O3');


cd /home/drholdaw/Lin_Moist_Physics/DAS_plotting/

plotlevel = 72;

Umax = max(max([U_dry(:,:,plotlevel); U_moi(:,:,plotlevel)]));
Umin = min(min([U_dry(:,:,plotlevel); U_moi(:,:,plotlevel)]));
Vmax = max(max([V_dry(:,:,plotlevel); V_moi(:,:,plotlevel)]));
Vmin = min(min([V_dry(:,:,plotlevel); V_moi(:,:,plotlevel)]));
TVmax = max(max([TV_dry(:,:,plotlevel); TV_moi(:,:,plotlevel)]));
TVmin = min(min([TV_dry(:,:,plotlevel); TV_moi(:,:,plotlevel)]));
DPmax = max(max([DP_dry(:,:,plotlevel); DP_moi(:,:,plotlevel)]));
DPmin = min(min([DP_dry(:,:,plotlevel); DP_moi(:,:,plotlevel)]));
QVmax = max(max([QV_dry(:,:,plotlevel); QV_moi(:,:,plotlevel)]));
QVmin = min(min([QV_dry(:,:,plotlevel); QV_moi(:,:,plotlevel)]));
% QLmax = max(max([QL_dry(:,:,plotlevel); QL_moi(:,:,plotlevel)]));
% QLmin = min(min([QL_dry(:,:,plotlevel); QL_moi(:,:,plotlevel)]));
% QImax = max(max([QI_dry(:,:,plotlevel); QI_moi(:,:,plotlevel)]));
% QImin = min(min([QI_dry(:,:,plotlevel); QI_moi(:,:,plotlevel)]));
% O3max = max(max([O3_dry(:,:,plotlevel); O3_moi(:,:,plotlevel)]));
% O3min = min(min([O3_dry(:,:,plotlevel); O3_moi(:,:,plotlevel)]));


figure

subplot(3,2,1)
contourf(lon,lat,U_dry(:,:,plotlevel)')
caxis([Umin Umax])
colorbar

subplot(3,2,2)
contourf(lon,lat,V_dry(:,:,plotlevel)')
caxis([Vmin Vmax])
colorbar

subplot(3,2,3)
contourf(lon,lat,TV_dry(:,:,plotlevel)')
caxis([TVmin TVmax])
colorbar

subplot(3,2,4)
contourf(lon,lat,QV_dry(:,:,plotlevel)')
caxis([QVmin QVmax])
colorbar

subplot(3,2,5)
contourf(lon,lat,DP_dry(:,:,plotlevel)')
caxis([DPmin DPmax])
colorbar


figure

subplot(3,2,1)
contourf(lon,lat,U_moi(:,:,plotlevel)')
caxis([Umin Umax])
colorbar

subplot(3,2,2)
contourf(lon,lat,V_moi(:,:,plotlevel)')
caxis([Vmin Vmax])
colorbar

subplot(3,2,3)
contourf(lon,lat,TV_moi(:,:,plotlevel)')
caxis([TVmin TVmax])
colorbar

subplot(3,2,4)
contourf(lon,lat,QV_moi(:,:,plotlevel)')
caxis([QVmin QVmax])
colorbar

subplot(3,2,5)
contourf(lon,lat,DP_moi(:,:,plotlevel)')
caxis([DPmin DPmax])
colorbar
