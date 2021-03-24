% close all
clear
clc
mydir = pwd;

cd /discover/nobackup/drholdaw/wrk.e5130/ADJRestart/
file = 'Jgradf.eta.nc4';

cd /discover/nobackup/drholdaw/wrk.e5130/sens/
file = 'e5130_dh.fsens.eta.20150412_0000z_P20000.nc4';

lon = ncread(file,'lon');
lat = ncread(file,'lat');
u = ncread(file,'u');
v = ncread(file,'v');
t = ncread(file,'tv');
q = ncread(file,'sphu');
p = ncread(file,'delp');
ql = ncread(file,'qltot');
qi = ncread(file,'qitot');
o3 = ncread(file,'ozone');

plot_level = 50;


figure
set(gcf,'position',[496 499 783 420])
contourf(lon,lat,u(:,:,plot_level)','LineStyle','none')
caxis([-7e-6 7e-6])
colorbar

figure
set(gcf,'position',[496 499 783 420])
contourf(lon,lat,v(:,:,plot_level)','LineStyle','none')
caxis([-9e-6 9e-6])
colorbar

figure
set(gcf,'position',[496 499 783 420])
contourf(lon,lat,t(:,:,plot_level)','LineStyle','none')
caxis([-8e-6 8e-6])
colorbar

asd

figure
set(gcf,'position',[496 499 783 420])
contourf(lon,lat,q(:,:,plot_level)','LineStyle','none')
colorbar

figure
set(gcf,'position',[496 499 783 420])
contourf(lon,lat,p(:,:,plot_level)','LineStyle','none')
colorbar

figure
set(gcf,'position',[496 499 783 420])
contourf(lon,lat,ql(:,:,plot_level)','LineStyle','none')
colorbar

figure
set(gcf,'position',[496 499 783 420])
contourf(lon,lat,qi(:,:,plot_level)','LineStyle','none')
colorbar




cd(mydir)