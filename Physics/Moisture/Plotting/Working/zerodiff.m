close all
clc
clear

dir1 = '/discover/nobackup/drholdaw/tmp.22292/traj1TMP/';
file1 = 'x0011dh_a.traj.lcv.20130117_0300z.nc4';

dir2 = '/discover/nobackup/drholdaw/tmp.22292/';
file2 = 'x0011dh_a.traj.lcv.20130117_0300z.nc4';

cd(dir1)

u_a = ncread(file1,'U');
v_a = ncread(file1,'V');
t_a = ncread(file1,'PT');
q_a = ncread(file1,'QV');
p_a = ncread(file1,'DP');
o3_a = ncread(file1,'QILS');

cd(dir2)

u_b = ncread(file2,'U');
v_b = ncread(file2,'V');
t_b = ncread(file2,'PT');
q_b = ncread(file2,'QV');
p_b = ncread(file2,'DP');
o3_b = ncread(file2,'QILS');


cd /home/drholdaw/Lin_Moist_Physics/DAS_plotting/


udiff = u_a(:,:,:)-u_b(:,:,:);
vdiff = v_a(:,:,:)-v_b(:,:,:);
tdiff = t_a(:,:,:)-t_b(:,:,:);
qdiff = q_a(:,:,:)-q_b(:,:,:);
pdiff = p_a(:,:,:)-p_b(:,:,:);
odiff = o3_a(:,:,:)-o3_b(:,:,:);

max(abs(udiff(:)))
max(abs(vdiff(:)))
max(abs(tdiff(:)))
max(abs(qdiff(:)))
max(abs(pdiff(:)))
max(abs(odiff(:)))


% contourf(q_a(:,:,63)','LineStyle','none')
% colorbar
% 
% figure
% contourf(q_b(:,:,63)','LineStyle','none')
% colorbar