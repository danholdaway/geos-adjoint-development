close all
clear
clc

cd /discover/nobackup/drholdaw/u001_C180_PH0/asens/stage/
file = 'u001_C180_PH0.fsens_txe.eta.20121231_15z+20130102_00z-20130101_00z.nc4';

lon = ncread(file,'lon');
lat = ncread(file,'lat');
lev = ncread(file,'lev');

u_ph15 = ncread(file,'u');
v_ph15 = ncread(file,'v');
t_ph15 = ncread(file,'tv');
q_ph15 = ncread(file,'sphu');

cd /discover/nobackup/drholdaw/u001_C180_PH0/asens/stage/
file = 'u001_C180_PH0.fsens_txe.eta.20121231_21z+20130102_00z-20130101_00z.nc4';

lon = ncread(file,'lon');
lat = ncread(file,'lat');
lev = ncread(file,'lev');

u_ph21 = ncread(file,'u');
v_ph21 = ncread(file,'v');
t_ph21 = ncread(file,'tv');
q_ph21 = ncread(file,'sphu');

cd /home/drholdaw/Lin_Moist_Physics/DAS_plotting/Cloud' Fraction'/


plot_level = 50;


u_15 = u_ph15(:,:,plot_level);
v_15 = v_ph15(:,:,plot_level);
t_15 = t_ph15(:,:,plot_level);
q_15 = q_ph15(:,:,plot_level);

u_21 = u_ph21(:,:,plot_level);
v_21 = v_ph21(:,:,plot_level);
t_21 = t_ph21(:,:,plot_level);
q_21 = q_ph21(:,:,plot_level);


figure
set(gcf,'position',[71 58 1173 837])

subplot(2,1,1)
contourf(lon,lat,q_15','LineStyle','none')
colorbar
caxis([-max(abs(u_15(:))) max(abs(u_15(:)))])

subplot(2,1,2)
contourf(lon,lat,q_21','LineStyle','none')
colorbar
caxis([-max(abs(u_15(:))) max(abs(u_15(:)))])