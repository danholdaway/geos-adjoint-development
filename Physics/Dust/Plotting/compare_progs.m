close all
clear
clc
mydir = pwd;

%Load MATLAB topography
load coast
coast_lat = lat; clear lat
coast_lon = long; clear long
grey = 0.75;

dir1 = '/discover/nobackup/drholdaw/btmp.25956/prog/prog_free';
dir2 = '/discover/nobackup/drholdaw/btmp.25956/prog/prog_replay';

file1 = 'v000_C180.prog.eta.20140201_0000z.nc4';
file2 = 'v000_C180.prog.eta.20140201_0000z.nc4';

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

cd(dir2)
file = file2;

DU001_2 = ncread(file,'DU001');
DU002_2 = ncread(file,'DU002');
DU003_2 = ncread(file,'DU003');
DU004_2 = ncread(file,'DU004');
DU005_2 = ncread(file,'DU005');

cd(mydir)

plot_lev = 65;

figure
set(gcf,'position',[719 74 560 845])

subplot(3,1,1)
contour(lon,lat,DU001_1(:,:,plot_lev)')
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])
colorbar
subplot(3,1,2)
contour(lon,lat,DU001_2(:,:,plot_lev)')
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])
colorbar
subplot(3,1,3)
contour(lon,lat,DU001_1(:,:,plot_lev)' - DU001_2(:,:,plot_lev)')
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])
colorbar

figure
set(gcf,'position',[719 74 560 845])

subplot(3,1,1)
contour(lon,lat,DU002_1(:,:,plot_lev)')
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])
colorbar
subplot(3,1,2)
contour(lon,lat,DU002_2(:,:,plot_lev)')
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])
colorbar
subplot(3,1,3)
contour(lon,lat,DU002_1(:,:,plot_lev)' - DU002_2(:,:,plot_lev)')
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])
colorbar

figure
set(gcf,'position',[719 74 560 845])

subplot(3,1,1)
contour(lon,lat,DU003_1(:,:,plot_lev)')
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])
colorbar
subplot(3,1,2)
contour(lon,lat,DU003_2(:,:,plot_lev)')
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])
colorbar
subplot(3,1,3)
contour(lon,lat,DU003_1(:,:,plot_lev)' - DU003_2(:,:,plot_lev)')
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])
colorbar

figure
set(gcf,'position',[719 74 560 845])

subplot(3,1,1)
contour(lon,lat,DU004_1(:,:,plot_lev)')
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])
colorbar
subplot(3,1,2)
contour(lon,lat,DU004_2(:,:,plot_lev)')
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])
colorbar
subplot(3,1,3)
contour(lon,lat,DU004_1(:,:,plot_lev)' - DU004_2(:,:,plot_lev)')
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])
colorbar

figure
set(gcf,'position',[719 74 560 845])

subplot(3,1,1)
contour(lon,lat,DU005_1(:,:,plot_lev)')
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])
colorbar
subplot(3,1,2)
contour(lon,lat,DU005_2(:,:,plot_lev)')
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])
colorbar
subplot(3,1,3)
contour(lon,lat,DU005_1(:,:,plot_lev)' - DU005_2(:,:,plot_lev)')
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])
colorbar
