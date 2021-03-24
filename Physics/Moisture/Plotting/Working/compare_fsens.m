close all
clear
clc

cd /archive/u/drholdaw/iau590a_mp/prog/Y2012/M03/D16/H15/

lon = ncread('iau590a_mp.fsens_txe.eta.20120316_15z+20120318_00z-20120317_00z.nc4','lon');
lat = ncread('iau590a_mp.fsens_txe.eta.20120316_15z+20120318_00z-20120317_00z.nc4','lat');

u_a = ncread('iau590a_mp.fsens_txe.eta.20120316_15z+20120318_00z-20120317_00z.nc4','u');
v_a = ncread('iau590a_mp.fsens_txe.eta.20120316_15z+20120318_00z-20120317_00z.nc4','u');
t_a = ncread('iau590a_mp.fsens_txe.eta.20120316_15z+20120318_00z-20120317_00z.nc4','tv');
q_a = ncread('iau590a_mp.fsens_txe.eta.20120316_15z+20120318_00z-20120317_00z.nc4','sphu');
p_a = ncread('iau590a_mp.fsens_txe.eta.20120316_15z+20120318_00z-20120317_00z.nc4','delp');

cd /archive/u/drholdaw/iau590a/prog_txe/Y2012/M03/D16/H15/

u_b = ncread('iau590a.fsens_txe.eta.20120316_15z+20120318_00z-20120317_00z.nc4','u');
v_b = ncread('iau590a.fsens_txe.eta.20120316_15z+20120318_00z-20120317_00z.nc4','u');
t_b = ncread('iau590a.fsens_txe.eta.20120316_15z+20120318_00z-20120317_00z.nc4','tv');
q_b = ncread('iau590a.fsens_txe.eta.20120316_15z+20120318_00z-20120317_00z.nc4','sphu');
p_b = ncread('iau590a.fsens_txe.eta.20120316_15z+20120318_00z-20120317_00z.nc4','delp');

cd /home/drholdaw/Lin_Moist_Physics/DAS_plotting/

u_e = u_a - u_b;
v_e = v_a - v_b;
t_e = t_a - t_b;
q_e = q_a - q_b;
p_e = p_a - p_b;


plot_level = 65;

figure

subplot(2,3,1)
contourf(lon,lat,u_a(:,:,plot_level)');
colorbar

subplot(2,3,2)
contourf(lon,lat,v_a(:,:,plot_level)');
colorbar

subplot(2,3,3)
contourf(lon,lat,t_a(:,:,plot_level)');
colorbar

subplot(2,3,4)
contourf(lon,lat,q_a(:,:,plot_level)');
colorbar

subplot(2,3,5)
contourf(lon,lat,p_a(:,:,plot_level)');
colorbar

figure

subplot(2,3,1)
contourf(lon,lat,u_b(:,:,plot_level)');
colorbar

subplot(2,3,2)
contourf(lon,lat,v_b(:,:,plot_level)');
colorbar

subplot(2,3,3)
contourf(lon,lat,t_b(:,:,plot_level)');
colorbar

subplot(2,3,4)
contourf(lon,lat,q_b(:,:,plot_level)');
colorbar

subplot(2,3,5)
contourf(lon,lat,p_b(:,:,plot_level)');
colorbar

figure

subplot(2,3,1)
contourf(lon,lat,u_e(:,:,plot_level)');
colorbar

subplot(2,3,2)
contourf(lon,lat,v_e(:,:,plot_level)');
colorbar

subplot(2,3,3)
contourf(lon,lat,t_e(:,:,plot_level)');
colorbar

subplot(2,3,4)
contourf(lon,lat,q_e(:,:,plot_level)');
colorbar

subplot(2,3,5)
contourf(lon,lat,p_e(:,:,plot_level)');
colorbar