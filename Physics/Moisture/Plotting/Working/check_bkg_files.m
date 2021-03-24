close all
clear
clc

load coast
coast_lat = lat; clear lat
coast_lon = long; clear long

cd /discover/nobackup/drholdaw/tmp.32203/

lon = ncread('x0009a.bkg.eta.20120702_21z.nc4','lat');
lat = ncread('x0009a.bkg.eta.20120702_21z.nc4','lon');

q = ncread('x0009a.bkg.eta.20120702_21z.nc4','u');
qitot = ncread('x0009a.bkg.eta.20120702_21z.nc4','qitot');
qltot = ncread('x0009a.bkg.eta.20120702_21z.nc4','qltot');

cd /home/drholdaw/Lin_Moist_Physics/DAS_plotting/


plot_level = 45;
grey = 0.5;

ql = q(:,:,plot_level);
for i = 1:576
    for j = 1:361
        qltotl(i,j) = sum(qltot(i,j,1:72));
        qitotl(i,j) = sum(qitot(i,j,1:72));
    end
end

figure
contourf(lat,lon,ql','LineStyle','none')
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])
colorbar

figure
contourf(lat,lon,qltotl','LineStyle','none')
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])
colorbar

figure
contourf(lat,lon,qitotl','LineStyle','none')
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])
colorbar