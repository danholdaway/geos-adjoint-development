close all
clear
clc
mydir = pwd;


dir1 = '/discover/nobackup/drholdaw/btmp.25956/sens.20140202.000000';
dir2 = '/discover/nobackup/drholdaw/btmp.25956/sens.20140202.000000';

file1 = 'fvpert.eta.20140201_04z.nc4';
file2 = 'fvpert.eta.20140201_04z_1.nc4';

cd(dir1)
file = file1;

lon = ncread(file,'lon');
lat = ncread(file,'lat');
lev = ncread(file,'lev');

u_1 = ncread(file,'U');
v_1 = ncread(file,'V');
t_1 = ncread(file,'TV');
q_1 = ncread(file,'QV');
p_1 = ncread(file,'DP');
qi_1 = ncread(file,'QI');
ql_1 = ncread(file,'QL');
o3_1 = ncread(file,'O3');

cd(dir2)
file = file2;

u_2 = ncread(file,'U');
v_2 = ncread(file,'V');
t_2 = ncread(file,'TV');
q_2 = ncread(file,'QV');
p_2 = ncread(file,'DP');
qi_2 = ncread(file,'QI');
ql_2 = ncread(file,'QL');
o3_2 = ncread(file,'O3');

cd(mydir)


x_all = [u_1-u_2; v_1-v_2; t_1-t_2; q_1-q_2; p_1-p_2; o3_1-o3_2; ql_1-ql_2; qi_1-qi_2];

max(abs(x_all(:)))

plot_lev = 65;
figure
contour(lon,lat,ql_1(:,:,plot_lev)')
% caxis([-3e-5 3e-5])
colorbar

figure
contour(lon,lat,ql_2(:,:,plot_lev)')
% caxis([-3e-5 3e-5])
colorbar

figure
contour(lon,lat,ql_1(:,:,plot_lev)'-ql_2(:,:,plot_lev)')
% caxis([-3e-5 3e-5])
colorbar


plot_lev = 45;
figure
contour(lon,lat,qi_1(:,:,plot_lev)')
% caxis([-15e-6 15e-6])
colorbar

figure
contour(lon,lat,qi_2(:,:,plot_lev)')
% caxis([-15e-6 15e-6])
colorbar

figure
contour(lon,lat,qi_1(:,:,plot_lev)'-qi_2(:,:,plot_lev)')
% caxis([-15e-6 15e-6])
colorbar




