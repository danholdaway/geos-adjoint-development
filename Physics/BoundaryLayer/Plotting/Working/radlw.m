close all
clear
clc

cd /discover/nobackup/drholdaw/tmp.22292/prog_free/
file = 'x0011dh_a.prog.eta.20130116_00z.nc4';

lon = ncread(file,'lon');
lat = ncread(file,'lat');
lev = ncread(file,'lev');

radlongwave = ncread(file,'DTDTLWR');

QILS = ncread(file,'qils');
QLLS = ncread(file,'qlls');
QICN = ncread(file,'qicn');
QLCN = ncread(file,'qlcn');

cd /home/drholdaw/Lin_Moist_Physics/DAS_plotting/BoundaryLayer/

CWI = QILS + QICN; 
CWL = QLLS + QLCN;

CW = CWI + CWL;

pl = 30;

figure
contour(lon,lat,radlongwave(:,:,pl)')
colorbar


figure
contour(lon,lat,CW(:,:,pl)')
colorbar
