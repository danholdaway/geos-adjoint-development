close all
clc
clear

cd /discover/nobackup/drholdaw/tmp.27169/sens.20130117.000000
file1 = 'fsens.eta.20130117_00z_0.nc4';

% cd /discover/nobackup/drholdaw/x0011dh_a/TuningExperiments/filt5_rmd8_adjd
% file1 = 'x0011dh_a.fsens_twe.eta.20130115_21z+20130117_00z-20130116_00z_old.nc4';

lon = ncread(file1,'lon');
lat = ncread(file1,'lat');

u_a = ncread(file1,'u');
v_a = ncread(file1,'v');
t_a = ncread(file1,'tv');
q_a = ncread(file1,'sphu');
p_a = ncread(file1,'delp');
o3_a = ncread(file1,'ozone');

cd /discover/nobackup/drholdaw/tmp.27169/sens.20130117.000000
file2 = 'fsens.eta.20130117_00z_sbac.nc4';

% cd /discover/nobackup/drholdaw/x0011dh_a/TuningExperiments/filt5_rmd5_adjd
% file2 = 'x0011dh_a.fsens_twe.eta.20130109_21z+20130111_00z-20130110_00z.nc4';

u_b = ncread(file2,'u');
v_b = ncread(file2,'v');
t_b = ncread(file2,'tv');
q_b = ncread(file2,'sphu');
p_b = ncread(file2,'delp');
o3_b = ncread(file2,'ozone');

cd /discover/nobackup/drholdaw/tmp.27169/sens.20130117.000000
file3 = 'fsens.eta.20130117_00z_2.nc4';

u_c = ncread(file3,'u');
v_c = ncread(file3,'v');
t_c = ncread(file3,'tv');
q_c = ncread(file3,'sphu');
p_c = ncread(file3,'delp');
o3_c = ncread(file3,'ozone');

% cd /discover/nobackup/drholdaw/tmp.26094/sens.20130103.000000
% file4 = 'fsens.eta.20130106_0000z_ssl.nc4';
% 
% u_d = ncread(file4,'u');
% v_d = ncread(file4,'v');
% t_d = ncread(file4,'tv');
% q_d = ncread(file4,'sphu');
% p_d = ncread(file4,'delp');
% o3_d = ncread(file4,'ozone');

cd /home/drholdaw/Lin_Moist_Physics/DAS_plotting/

a = q_a;
b = q_b;
c = q_c;
% d = q_d;

plot_level = 63;

a_lev = a(:,:,plot_level);
b_lev = b(:,:,plot_level);
c_lev = c(:,:,plot_level);

max(max(abs(c_lev - a_lev)))

% d_lev = d(:,:,plot_level);

maxi_a = max(abs(a_lev(:)));
maxi_b = max(abs(b_lev(:)));
maxi_c = max(abs(c_lev(:)));
% maxi_d = max(abs(d_lev(:)));

cint = 2*maxi_c /200;

cint1 = 2*maxi_b /200;

contourf(lon,lat,((a_lev))','LineStyle','none')%,'LevelStep',cint)
colorbar
caxis([-maxi_a maxi_a])
% ylim([56 61])
% xlim([-64 -58])

figure
contourf(lon,lat,((b_lev))','LineStyle','none')%,'LevelStep',cint1)
colorbar
% colormap(cmap)
caxis([-maxi_b maxi_b])
% ylim([56 61])
% xlim([-64 -58])

figure
contourf(lon,lat,((c_lev))','LineStyle','none')%,'LevelStep',cint)
colorbar
% colormap(cmap)
caxis([-maxi_c maxi_c])
% ylim([56 61])
% xlim([-64 -58])

% figure
% contourf(lon,lat,d(:,:,plot_level)','LineStyle','none')%,'LevelStep',cint)
% colorbar
% % colormap(cmap)
% caxis([-maxi_d maxi_d])
% % ylim([56 61])
% % xlim([-64 -58])close all
clc
clear

cd /discover/nobackup/drholdaw/tmp.27169/sens.20130117.000000
file1 = 'fsens.eta.20130117_00z_0.nc4';

% cd /discover/nobackup/drholdaw/x0011dh_a/TuningExperiments/filt5_rmd8_adjd
% file1 = 'x0011dh_a.fsens_twe.eta.20130115_21z+20130117_00z-20130116_00z_old.nc4';

lon = ncread(file1,'lon');
lat = ncread(file1,'lat');

u_a = ncread(file1,'u');
v_a = ncread(file1,'v');
t_a = ncread(file1,'tv');
q_a = ncread(file1,'sphu');
p_a = ncread(file1,'delp');
o3_a = ncread(file1,'ozone');

cd /discover/nobackup/drholdaw/tmp.27169/sens.20130117.000000
file2 = 'fsens.eta.20130117_00z_sbac.nc4';

% cd /discover/nobackup/drholdaw/x0011dh_a/TuningExperiments/filt5_rmd5_adjd
% file2 = 'x0011dh_a.fsens_twe.eta.20130109_21z+20130111_00z-20130110_00z.nc4';

u_b = ncread(file2,'u');
v_b = ncread(file2,'v');
t_b = ncread(file2,'tv');
q_b = ncread(file2,'sphu');
p_b = ncread(file2,'delp');
o3_b = ncread(file2,'ozone');

cd /discover/nobackup/drholdaw/tmp.27169/sens.20130117.000000
file3 = 'fsens.eta.20130117_00z_2.nc4';

u_c = ncread(file3,'u');
v_c = ncread(file3,'v');
t_c = ncread(file3,'tv');
q_c = ncread(file3,'sphu');
p_c = ncread(file3,'delp');
o3_c = ncread(file3,'ozone');

% cd /discover/nobackup/drholdaw/tmp.26094/sens.20130103.000000
% file4 = 'fsens.eta.20130106_0000z_ssl.nc4';
% 
% u_d = ncread(file4,'u');
% v_d = ncread(file4,'v');
% t_d = ncread(file4,'tv');
% q_d = ncread(file4,'sphu');
% p_d = ncread(file4,'delp');
% o3_d = ncread(file4,'ozone');

cd /home/drholdaw/Lin_Moist_Physics/DAS_plotting/

a = q_a;
b = q_b;
c = q_c;
% d = q_d;

plot_level = 63;

a_lev = a(:,:,plot_level);
b_lev = b(:,:,plot_level);
c_lev = c(:,:,plot_level);

max(max(abs(c_lev - a_lev)))

% d_lev = d(:,:,plot_level);

maxi_a = max(abs(a_lev(:)));
maxi_b = max(abs(b_lev(:)));
maxi_c = max(abs(c_lev(:)));
% maxi_d = max(abs(d_lev(:)));

cint = 2*maxi_c /200;

cint1 = 2*maxi_b /200;

contourf(lon,lat,((a_lev))','LineStyle','none')%,'LevelStep',cint)
colorbar
caxis([-maxi_a maxi_a])
% ylim([56 61])
% xlim([-64 -58])

figure
contourf(lon,lat,((b_lev))','LineStyle','none')%,'LevelStep',cint1)
colorbar
% colormap(cmap)
caxis([-maxi_b maxi_b])
% ylim([56 61])
% xlim([-64 -58])

figure
contourf(lon,lat,((c_lev))','LineStyle','none')%,'LevelStep',cint)
colorbar
% colormap(cmap)
caxis([-maxi_c maxi_c])
% ylim([56 61])
% xlim([-64 -58])

% figure
% contourf(lon,lat,d(:,:,plot_level)','LineStyle','none')%,'LevelStep',cint)
% colorbar
% % colormap(cmap)
% caxis([-maxi_d maxi_d])
% % ylim([56 61])
% % xlim([-64 -58])