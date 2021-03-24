close all
clear
clc
mydir = pwd;

%Load MATLAB topography
load coast
coast_lat = lat; clear lat
coast_lon = long; clear long
grey = 0.75;

dir1 = '/discover/nobackup/drholdaw/btmp.25956/tlmrestarts_03z';

file1 = 'fvpert.eta.nc4';

cd(dir1)
file = file1;

lon = ncread(file,'lon');
lat = ncread(file,'lat');
lev = ncread(file,'lev');

DU001_1 = ncread(file,'DU001');
DU002_1 = ncread(file,'DU002');
DU003_1 = ncread(file,'DU003');
DU004_1 = ncread(file,'DU004');
DU005_1 = ncread(file,'DU005');

plot_lev = 65;

figure

contour(lon,lat,DU001_1(:,:,plot_lev)')
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])
colorbar