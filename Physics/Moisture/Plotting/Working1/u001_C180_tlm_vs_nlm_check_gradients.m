close all
clear
clc

cd /archive/u/drholdaw/u001_C180_PH1/prog/Y2013/M01/D02/H21
file = 'u001_C180_PH1.fsens_txe.eta.20130102_21z+20130104_00z-20130103_00z.nc4';

lon = ncread(file,'lon');
lat = ncread(file,'lat');
lev = ncread(file,'lev');

u_ph15 = ncread(file,'u');
v_ph15 = ncread(file,'v');
t_ph15 = ncread(file,'tv');
q_ph15 = ncread(file,'sphu');
p_ph15 = ncread(file,'delp');
qi_ph15 = ncread(file,'qitot');
ql_ph15 = ncread(file,'qltot');
o3_ph15 = ncread(file,'ozone');

% cd /archive/u/drholdaw/u001_C180_PH0/prog/Y2013/M01/D02/H15
cd /discover/nobackup/drholdaw/atmp.22292/sens.20130117.000000/
% file = 'u001_C180_PH0.fsens_txe.eta.20130102_15z+20130104_00z-20130103_00z.nc4';
file = 'u001_C180_PH1.fsens_txe.eta.20130102_21z+20130104_00z-20130103_00z.nc4'

u_ph21 = ncread(file,'u');
v_ph21 = ncread(file,'v');
t_ph21 = ncread(file,'tv');
q_ph21 = ncread(file,'sphu');
p_ph21 = ncread(file,'delp');
qi_ph21 = ncread(file,'qitot');
ql_ph21 = ncread(file,'qltot');
o3_ph21 = ncread(file,'ozone');

cd /home/drholdaw/Lin_Moist_Physics/DAS_plotting/Cloud' Fraction'/

u_er = u_ph21 - u_ph15;
v_er = v_ph21 - v_ph15;
t_er = t_ph21 - t_ph15;
q_er = q_ph21 - q_ph15;
p_er = p_ph21 - p_ph15;
qi_er = qi_ph21 - qi_ph15;
ql_er = ql_ph21 - ql_ph15;
o3_er = o3_ph21 - o3_ph15;

er = [u_er; v_er; t_er; q_er; p_er; qi_er; ql_er; o3_er];

max(abs(er(:)))

for plot_level = 40%:72;


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
    contourf(lon,lat,u_15','LineStyle','none')
    colorbar
    caxis([-max(abs(u_15(:))) max(abs(u_15(:)))])

    subplot(2,1,2)
    contourf(lon,lat,u_21','LineStyle','none')
    colorbar
    caxis([-max(abs(u_15(:))) max(abs(u_15(:)))])
    
%     pause
%     close
    
end