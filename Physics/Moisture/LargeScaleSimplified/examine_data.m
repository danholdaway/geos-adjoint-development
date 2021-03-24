close all
clear
clc

cd ../ras_offline
grid_data

lev = 1:Nk;

cd /discover/nobackup/drholdaw/inputs_forRAS/high_vs_low/

U_start = ncread('high_vs_low_14491.geosgcm_rasinputs3d.20100702_1200z.nc4','U_preras');
V_start = ncread('high_vs_low_14491.geosgcm_rasinputs3d.20100702_1200z.nc4','V_preras');
TH_start = ncread('high_vs_low_14491.geosgcm_rasinputs3d.20100702_1200z.nc4','TH_preras');
Q_start = ncread('high_vs_low_14491.geosgcm_rasinputs3d.20100702_1200z.nc4','Q_preras');

Q_ras_GEOScheck = ncread('high_vs_low_14491.geosgcm_rasinputs3d.20100702_1200z.nc4','Q_total');

PRC2_ras_GEOScheck =  ncread('high_vs_low_14491.geosgcm_surf.20100702_1330z.nc4','PCU');
PRC2_lsc_GEOScheck =  ncread('high_vs_low_14491.geosgcm_surf.20100702_1330z.nc4','PLS');

cd /discover/nobackup/drholdaw/ExperimentData/zerodifftest/

U_ras = ncread('3DFields_bothmoist_indepsubs.nc4','U_moist');
V_ras = ncread('3DFields_bothmoist_indepsubs.nc4','V_moist');
TH_ras = ncread('3DFields_bothmoist_indepsubs.nc4','TH_moist');
Q_ras = ncread('3DFields_bothmoist_indepsubs.nc4','Q_moist');

PRC3_ras = ncread('3DFields_bothmoist_indepsubs.nc4','PRC3_ras');
PRC3_lsc = ncread('3DFields_bothmoist_indepsubs.nc4','PRC3_lsc');

cd /home/drholdaw/moist_adjoint/largecloud


TH_diff = U_start - U_ras;
%FIND HIGHEST LEVEL WHERE MOIST ACTS
for i = 1:72
    
    maxatlev(i) = max(max(TH_diff(:,:,i)));
    
end
maxatlev

PRC2_ras = zeros(Ni,Nj);
PRC2_lsc = zeros(Ni,Nj);

for i = 1:Ni
    for j = 1:Nj
        PRC2_ras(i,j) = sum(PRC3_ras(i,j,:));
        PRC2_lsc(i,j) = sum(PRC3_lsc(i,j,:));
    end
end

% figure
% contour(PRC2_ras')
% colorbar
% 
% figure
% contour(PRC2_lsc')
% colorbar
% 
% figure
% contour(PRC2_ras_GEOScheck')
% colorbar
% 
% figure
% contour(PRC2_lsc_GEOScheck')
% colorbar

level = 39;

figure
contour(TH_diff(:,:,level)')
colorbar
% 
% figure
% contour(Q_lsc(:,:,level)')
% colorbar