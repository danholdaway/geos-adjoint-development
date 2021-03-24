close all
clc
clear

cd /discover/nobackup/drholdaw/tmp.22292/sens.20130117.000000/
file1 = 'Jgradf_twe.eta.nc4';

lon = ncread(file1,'lon');
lat = ncread(file1,'lat');

u_a = ncread(file1,'u');
v_a = ncread(file1,'v');
t_a = ncread(file1,'tv');
q_a = ncread(file1,'sphu');
p_a = ncread(file1,'delp');
qi_a = ncread(file1,'qitot');
ql_a = ncread(file1,'qltot');
o3_a = ncread(file1,'ozone');

[a,b,c] = size(u_a);

q_a = 10e-4*q_a;
qi_a = 10e-4*qi_a;
ql_a = 10e-4*ql_a;

ncwrite('Jgradf_twe.eta.nc4','sphu',q_a);
ncwrite('Jgradf_twe.eta.nc4','qitot',qi_a);
ncwrite('Jgradf_twe.eta.nc4','qltot',ql_a);


asd

cd /discover/nobackup/drholdaw/tmp.22292/sens.20130117.000000/
file2 = 'fvpertX.eta.nc4';

u_b = ncread(file2,'U');
v_b = ncread(file2,'V');
t_b = ncread(file2,'TV');
q_b = ncread(file2,'QV');
p_b = ncread(file2,'DP');
qi_b = ncread(file2,'QI');
ql_b = ncread(file2,'QL');
o3_b = ncread(file2,'O3');

qi_b = zeros(a,b,c);
ql_b = zeros(a,b,c);
o3_b = zeros(a,b,c);

ncwrite('fvpertX.eta.nc4','QI',qi_b);
ncwrite('fvpertX.eta.nc4','QL',ql_b);
ncwrite('fvpertX.eta.nc4','O3',o3_b);

cd /home/drholdaw/Lin_Moist_Physics/DAS_plotting/

% asd

clear


cd /discover/nobackup/drholdaw/tmp.22292/sens.20130117.000000/
file1 = 'fvpert.eta.nc4';

lon = ncread(file1,'lon');
lat = ncread(file1,'lat');

u_a = ncread(file1,'u');
v_a = ncread(file1,'v');
t_a = ncread(file1,'tv');
q_a = ncread(file1,'sphu');
p_a = ncread(file1,'delp');
qi_a = ncread(file1,'qitot');
ql_a = ncread(file1,'qltot');
o3_a = ncread(file1,'ozone');

cd /discover/nobackup/drholdaw/tmp.22292/sens.20130117.000000/
file2 = 'fvpertX.eta.nc4';

u_b = ncread(file2,'U');
v_b = ncread(file2,'V');
t_b = ncread(file2,'TV');
q_b = ncread(file2,'QV');
p_b = ncread(file2,'DP');
qi_b = ncread(file2,'QI');
ql_b = ncread(file2,'QL');
o3_b = ncread(file2,'O3');

cd /home/drholdaw/Lin_Moist_Physics/DAS_plotting/





plot_level = 63;

u_a_lev = u_a(:,:,plot_level);
v_a_lev = v_a(:,:,plot_level);
t_a_lev = t_a(:,:,plot_level);
q_a_lev = q_a(:,:,plot_level);
p_a_lev = p_a(:,:,plot_level);
qi_a_lev = qi_a(:,:,plot_level);
ql_a_lev = ql_a(:,:,plot_level);
o3_a_lev = o3_a(:,:,plot_level);

u_b_lev = u_b(:,:,plot_level);
v_b_lev = v_b(:,:,plot_level);
t_b_lev = t_b(:,:,plot_level);
q_b_lev = q_b(:,:,plot_level);
p_b_lev = p_b(:,:,plot_level);
qi_b_lev = qi_b(:,:,plot_level);
ql_b_lev = ql_b(:,:,plot_level);
o3_b_lev = o3_b(:,:,plot_level);

%Zerodiff
all_a = [u_a v_a t_a q_a p_a qi_a ql_a o3_a];
all_b = [u_b v_b t_b q_b p_b qi_b ql_b o3_b];

all_diff = all_a - all_b;
zd = max(abs(all_diff(:)));

if zd == 0 
    fprintf('Zero Diff Passed \n\n');
else
    fprintf('Not Zero Diff \n\n');
end

%Get Maximums for Levels
u_lev = [u_a_lev u_b_lev];
maxu = max(abs(u_lev(:)));
v_lev = [v_a_lev v_b_lev];
maxv = max(abs(v_lev(:)));
t_lev = [t_a_lev t_b_lev];
maxt = max(abs(t_lev(:)));
q_lev = [q_a_lev q_b_lev];
maxq = max(abs(q_lev(:)));
p_lev = [p_a_lev p_b_lev];
maxp = max(abs(p_lev(:)));
qi_lev = [qi_a_lev qi_b_lev];
maxqi = max(abs(qi_lev(:)));
ql_lev = [ql_a_lev ql_b_lev];
maxql = max(abs(ql_lev(:)));
o3_lev = [o3_a_lev o3_b_lev];
maxo3 = max(abs(o3_lev(:)));



figure
set(gcf,'position',[349 30 930 889])

subplot(4,2,1)
contourf(lon,lat,((u_a_lev))','LineStyle','none')
caxis([-maxu maxu])
title('u')
colorbar

subplot(4,2,2)
contourf(lon,lat,((v_a_lev))','LineStyle','none')
caxis([-maxv maxv])
title('v')
colorbar

subplot(4,2,3)
contourf(lon,lat,((t_a_lev))','LineStyle','none')
caxis([-maxt maxt])
title('\theta')
colorbar

subplot(4,2,4)
contourf(lon,lat,((q_a_lev))','LineStyle','none')
caxis([-maxq maxq])
title('q')
colorbar

subplot(4,2,5)
contourf(lon,lat,((p_a_lev))','LineStyle','none')
caxis([-maxp maxp])
title('dp')
colorbar

subplot(4,2,6)
contourf(lon,lat,((qi_a_lev))','LineStyle','none')
caxis([-maxqi maxqi])
title('Qi')
colorbar

subplot(4,2,7)
contourf(lon,lat,((ql_a_lev))','LineStyle','none')
caxis([-maxql maxql])
title('Ql')
colorbar

subplot(4,2,8)
contourf(lon,lat,((o3_a_lev))','LineStyle','none')
caxis([-maxo3 maxo3])
title('O3')
colorbar


figure
set(gcf,'position',[349 30 930 889])

subplot(4,2,1)
contourf(lon,lat,((u_b_lev))','LineStyle','none')
caxis([-maxu maxu])
title('u')
colorbar

subplot(4,2,2)
contourf(lon,lat,((v_b_lev))','LineStyle','none')
caxis([-maxv maxv])
title('v')
colorbar

subplot(4,2,3)
contourf(lon,lat,((t_b_lev))','LineStyle','none')
caxis([-maxt maxt])
title('\theta')
colorbar

subplot(4,2,4)
contourf(lon,lat,((q_b_lev))','LineStyle','none')
caxis([-maxq maxq])
title('q')
colorbar

subplot(4,2,5)
contourf(lon,lat,((p_b_lev))','LineStyle','none')
caxis([-maxp maxp])
title('dp')
colorbar

subplot(4,2,6)
contourf(lon,lat,((qi_b_lev))','LineStyle','none')
caxis([-maxqi maxqi])
title('Qi')
colorbar

subplot(4,2,7)
contourf(lon,lat,((ql_b_lev))','LineStyle','none')
caxis([-maxql maxql])
title('Ql')
colorbar

subplot(4,2,8)
contourf(lon,lat,((o3_b_lev))','LineStyle','none')
caxis([-maxo3 maxo3])
title('O3')
colorbar