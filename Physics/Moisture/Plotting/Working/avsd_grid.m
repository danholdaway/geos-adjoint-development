close all
clear
clc

cd /discover/nobackup/drholdaw/NormCreator/

lat = ncread('verify_fcst.nc4','lat');
lon = ncread('verify_fcst.nc4','lon');
u = ncread('x0011dh_a.prog.eta.20130117_00z.nc4','u');
t = ncread('verify_fcst.nc4','tv');

cd /home/drholdaw/Lin_Moist_Physics/DAS_plotting/

ul = u(:,:,50)';
tl = t(:,:,50)';

contourf(lon,lat,ul,'LineStyle','none')
colorbar

figure
contourf(lon,lat,tl,'LineStyle','none')
colorbar