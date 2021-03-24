close all
clear
clc
mydir = pwd;

%Load MATLAB topography
load coast
coast_lat = lat; clear lat
coast_lon = long; clear long
grey = 0.75;
%Make center of colormap white
cmap = colormap;
cmap(32,:) = [1 1 1];
cmap(33,:) = [1 1 1];
close


im = 576;
jm = 361;
lm = 72;

cd /home/drholdaw/LinearisedPhysics/Inputs/
pref

%Load Free (background) State.
cd /discover/nobackup/drholdaw/btmp.25956/sens.20140202.000000/
% file = 'v000_C180.fvsens.eta.20140201_0000z_ST03z_DUST_P20000_S0.nc4';
file = 'v000_C180.fsens.eta.20140202_0000z-20140201_0000z_DUST_P220NODUST0_S0.nc4';

lon = ncread(file,'lon');
lat = ncread(file,'lat');
lev = ncread(file,'lev');

u_adm = ncread(file,'u');
v_adm = ncread(file,'v');
t_adm = ncread(file,'tv');
q_adm = ncread(file,'sphu');
p_adm = ncread(file,'delp');
qi_adm = ncread(file,'qitot');
ql_adm = ncread(file,'qltot');
o3_adm = ncread(file,'ozone');
d1_adm = ncread(file,'DU001');
d2_adm = ncread(file,'DU002');
d3_adm = ncread(file,'DU003');
d4_adm = ncread(file,'DU004');
d5_adm = ncread(file,'DU005');

cd /discover/nobackup/drholdaw/btmp.25956/prog/prog_free/
file = 'v000_C180.prog.eta.20140201_0000z.nc4';

u_nlm = ncread(file,'u');
v_nlm = ncread(file,'v');
t_nlm = ncread(file,'tv');
q_nlm = ncread(file,'sphu');
p_nlm = ncread(file,'delp');
qi_nlm = ncread(file,'qitot');
ql_nlm = ncread(file,'qltot');
o3_nlm = ncread(file,'ozone');
d1_nlm = ncread(file,'DU001');
d2_nlm = ncread(file,'DU002');
d3_nlm = ncread(file,'DU003');
d4_nlm = ncread(file,'DU004');
d5_nlm = ncread(file,'DU005');

cd(mydir)



plot_level = 65;

u_adm_l = u_adm(:,:,plot_level);
v_adm_l = v_adm(:,:,plot_level);
t_adm_l = t_adm(:,:,plot_level);
q_adm_l = q_adm(:,:,plot_level);
ql_adm_l = ql_adm(:,:,plot_level);
qi_adm_l = qi_adm(:,:,plot_level);
d1_adm_l = d1_adm(:,:,plot_level);
d2_adm_l = d2_adm(:,:,plot_level);
d3_adm_l = d3_adm(:,:,plot_level);
d4_adm_l = d4_adm(:,:,plot_level);
d5_adm_l = d5_adm(:,:,plot_level);
d_adm_l  = d1_adm(:,:,plot_level) + d2_adm(:,:,plot_level) + d3_adm(:,:,plot_level) + d4_adm(:,:,plot_level) + d5_adm(:,:,plot_level);

u_nlm_l = u_nlm(:,:,plot_level);
v_nlm_l = v_nlm(:,:,plot_level);
t_nlm_l = t_nlm(:,:,plot_level);
q_nlm_l = q_nlm(:,:,plot_level);
ql_nlm_l = ql_nlm(:,:,plot_level);
qi_nlm_l = qi_nlm(:,:,plot_level);
d1_nlm_l = d1_nlm(:,:,plot_level);
d2_nlm_l = d2_nlm(:,:,plot_level);
d3_nlm_l = d3_nlm(:,:,plot_level);
d4_nlm_l = d4_nlm(:,:,plot_level);
d5_nlm_l = d5_nlm(:,:,plot_level);
d_nlm_l  = d1_nlm(:,:,plot_level) + d2_nlm(:,:,plot_level) + d3_nlm(:,:,plot_level) + d4_nlm(:,:,plot_level) + d5_nlm(:,:,plot_level);

u_adm_scale_l = u_adm_l.*(abs(u_nlm_l)./max(abs(u_nlm_l(:))));
v_adm_scale_l = v_adm_l.*(abs(v_nlm_l)./max(abs(v_nlm_l(:))));
t_adm_scale_l = t_adm_l.*(abs(t_nlm_l)./max(abs(t_nlm_l(:))));
q_adm_scale_l = q_adm_l.*(abs(q_nlm_l)./max(abs(q_nlm_l(:))));
ql_adm_scale_l = ql_adm_l.*(abs(ql_nlm_l)./max(abs(ql_nlm_l(:))));
qi_adm_scale_l = qi_adm_l.*(abs(qi_nlm_l)./max(abs(qi_nlm_l(:))));
d1_adm_scale_l = d1_adm_l.*(abs(d1_nlm_l)./max(abs(d1_nlm_l(:))));
d2_adm_scale_l = d2_adm_l.*(abs(d2_nlm_l)./max(abs(d2_nlm_l(:))));
d3_adm_scale_l = d3_adm_l.*(abs(d3_nlm_l)./max(abs(d3_nlm_l(:))));
d4_adm_scale_l = d4_adm_l.*(abs(d4_nlm_l)./max(abs(d4_nlm_l(:))));
d5_adm_scale_l = d5_adm_l.*(abs(d5_nlm_l)./max(abs(d5_nlm_l(:))));
d_adm_scale_l = d_adm_l.*(abs(d_nlm_l)./max(abs(d_nlm_l(:))));

figure
set(gcf,'position',[97 86 1131 828])

subplot(3,2,1)
contour(lon,lat,d_adm_l')
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])
caxis([-max(abs(d_adm_l(:))) max(abs(d_adm_l(:)))])
colorbar
title('Dust Sum - Sensitivity')
colormap(cmap)

subplot(3,2,2)
contour(lon,lat,d1_adm_l')
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])
caxis([-max(abs(d1_adm_l(:))) max(abs(d1_adm_l(:)))])
colorbar
title('Dust bin 1 - Sensitivity')
colormap(cmap)

subplot(3,2,3)
contour(lon,lat,d2_adm_l')
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])
caxis([-max(abs(d2_adm_l(:))) max(abs(d2_adm_l(:)))])
colorbar
title('Dust bin 2 - Sensitivity')
colormap(cmap)

subplot(3,2,4)
contour(lon,lat,d3_adm_l')
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])
caxis([-max(abs(d3_adm_l(:))) max(abs(d3_adm_l(:)))])
colorbar
title('Dust bin 3 - Sensitivity')
colormap(cmap)

subplot(3,2,5)
contour(lon,lat,d4_adm_l')
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])
caxis([-max(abs(d4_adm_l(:))) max(abs(d4_adm_l(:)))])
colorbar
title('Dust bin 4 - Sensitivity')
colormap(cmap)

subplot(3,2,6)
contour(lon,lat,d5_adm_l')
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])
caxis([-max(abs(d5_adm_l(:))) max(abs(d5_adm_l(:)))])
colorbar
title('Dust bin 5 - Sensitivity')
colormap(cmap)

figure
set(gcf,'position',[97 86 1131 828])

subplot(3,2,1)
contour(lon,lat,d_adm_scale_l')
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])
caxis([-max(abs(d_adm_scale_l(:))) max(abs(d_adm_scale_l(:)))])
colorbar
title('Dust Sum - Sensitivity')
colormap(cmap)

subplot(3,2,2)
contour(lon,lat,d1_adm_scale_l')
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])
caxis([-max(abs(d1_adm_scale_l(:))) max(abs(d1_adm_scale_l(:)))])
colorbar
title('Dust bin 1 - Sensitivity')
colormap(cmap)

subplot(3,2,3)
contour(lon,lat,d2_adm_scale_l')
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])
caxis([-max(abs(d2_adm_scale_l(:))) max(abs(d2_adm_scale_l(:)))])
colorbar
title('Dust bin 2 - Sensitivity')
colormap(cmap)

subplot(3,2,4)
contour(lon,lat,d3_adm_scale_l')
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])
caxis([-max(abs(d3_adm_scale_l(:))) max(abs(d3_adm_scale_l(:)))])
colorbar
title('Dust bin 3 - Sensitivity')
colormap(cmap)

subplot(3,2,5)
contour(lon,lat,d4_adm_scale_l')
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])
caxis([-max(abs(d4_adm_scale_l(:))) max(abs(d4_adm_scale_l(:)))])
colorbar
title('Dust bin 4 - Sensitivity')
colormap(cmap)

subplot(3,2,6)
contour(lon,lat,d5_adm_scale_l')
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])
caxis([-max(abs(d5_adm_scale_l(:))) max(abs(d5_adm_scale_l(:)))])
colorbar
title('Dust bin 5 - Sensitivity')
colormap(cmap)

figure
set(gcf,'position',[97 86 1131 828])

subplot(3,2,1)
contour(lon,lat,u_adm_l')
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])
caxis([-max(abs(u_adm_l(:))) max(abs(u_adm_l(:)))])
colorbar
title('u - Sensitivity')
colormap(cmap)

subplot(3,2,2)
contour(lon,lat,v_adm_l')
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])
caxis([-max(abs(v_adm_l(:))) max(abs(v_adm_l(:)))])
colorbar
title('v - Sensitivity')
colormap(cmap)

subplot(3,2,3)
contour(lon,lat,t_adm_l')
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])
caxis([-max(abs(t_adm_l(:))) max(abs(t_adm_l(:)))])
colorbar
title('T_v - Sensitivity')
colormap(cmap)

subplot(3,2,4)
contour(lon,lat,q_adm_l')
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])
caxis([-max(abs(q_adm_l(:))) max(abs(q_adm_l(:)))])
colorbar
title('q - Sensitivity')
colormap(cmap)

subplot(3,2,5)
contour(lon,lat,ql_adm_l')
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])
caxis([-max(abs(ql_adm_l(:))) max(abs(ql_adm_l(:)))])
colorbar
title('ql - Sensitivity')
colormap(cmap)

subplot(3,2,6)
contour(lon,lat,qi_adm_l')
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])
caxis([-max(abs(d5_adm_l(:))) max(abs(d5_adm_l(:)))])
colorbar
title('qi - Sensitivity')
colormap(cmap)





figure
set(gcf,'position',[1380 78 1131 828])

subplot(3,2,1)
contour(lon,lat,d_nlm_l')
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])
colorbar
title('Dust Sum - Forecast')

subplot(3,2,2)
contour(lon,lat,d1_nlm_l')
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])
colorbar
title('Dust bin 1 - Forecast')

subplot(3,2,3)
contour(lon,lat,d2_nlm_l')
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])
colorbar
title('Dust bin 2 - Forecast')

subplot(3,2,4)
contour(lon,lat,d3_nlm_l')
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])
colorbar
title('Dust bin 3 - Forecast')

subplot(3,2,5)
contour(lon,lat,d4_nlm_l')
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])
colorbar
title('Dust bin 4 - Forecast')

subplot(3,2,6)
contour(lon,lat,d5_nlm_l')
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])
colorbar
title('Dust bin 5 - Forecast')

figure
set(gcf,'position',[1380 78 1131 828])

subplot(3,2,1)
contour(lon,lat,u_nlm_l')
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])
caxis([-max(abs(u_nlm_l(:))) max(abs(u_nlm_l(:)))])
colorbar
title('u - Forecast')

subplot(3,2,2)
contour(lon,lat,v_nlm_l')
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])
caxis([-max(abs(v_nlm_l(:))) max(abs(v_nlm_l(:)))])
colorbar
title('v - Forecast')

subplot(3,2,3)
contour(lon,lat,t_nlm_l')
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])
colorbar
title('T_v - Forecast')

subplot(3,2,4)
contour(lon,lat,q_nlm_l')
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])
colorbar
title('q - Forecast')

subplot(3,2,5)
contour(lon,lat,ql_nlm_l')
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])
colorbar
title('ql - Forecast')

subplot(3,2,6)
contour(lon,lat,qi_nlm_l')
hold on; 
plot(coast_lon,coast_lat,'Color',[grey grey grey])
colorbar
title('qi - Forecast')


