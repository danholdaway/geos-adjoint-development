close all
clear
clc
mydir = pwd;

%Load MATLAB topography
load coast
coast_lat = lat; clear lat
coast_lon = long; clear long
grey = 0.75;

cd /archive/u/drholdaw/C360_513_snamma/prog/Y2006/M08/D26/H21/

file = 'C360_513_snamma.prog.inst3d_met_p.20060826_21z+20060828_00z.nc4';

lon = ncread(file,'lon');
lat = ncread(file,'lat');
SLP = ncread(file,'SLP');

cd(mydir)

imin = find(lon == -80);
imax = find(lon ==  10);
jmin = find(lat ==   0);
jmax = find(lat ==  40);

plot(coast_lon,coast_lat,'Color',[grey grey grey])
hold on
contourf(lon(imin:imax),lat(jmin:jmax),SLP(imin:imax,jmin:jmax)');
plot(coast_lon,coast_lat,'Color',[grey grey grey])
xlim([-80 10])
ylim([0 40])

colorbar
