close all
clear
clc

% cd /gpfsm/dnb31/drholdaw/tmp.22292/nlm_runs/free_qlscn/
% file1 = 'x0011dh_a.prog.eta.20130116_00z.nc4';

cd /gpfsm/dnb31/drholdaw/tmp.22292/
file1 = 'x0011dh_a.traj.lcv.20130116_0000z.nc4';

lon = ncread(file1,'lon');
lat = ncread(file1,'lat');

% qi_a = ncread(file1,'qitot');
% ql_a = ncread(file1,'qltot');
% qils_a = ncread(file1,'qils');
% qlls_a = ncread(file1,'qlls');
% qicn_a = ncread(file1,'qicn');
% qlcn_a = ncread(file1,'qlcn');
% cfls_a = ncread(file1,'cfls');
% cfcn_a = ncread(file1,'cfcn');

qi_a = ncread(file1,'QI');
ql_a = ncread(file1,'QL');
qils_a = ncread(file1,'QILS');
qlls_a = ncread(file1,'QLLS');
qicn_a = ncread(file1,'QICN');
qlcn_a = ncread(file1,'QLCN');
cfls_a = ncread(file1,'CFLS');
cfcn_a = ncread(file1,'CFCN');

cd /home/drholdaw/Lin_Moist_Physics/DAS_plotting/Cloud' Fraction'/

plot_level = 40;

[im,jm,lm] = size(qi_a);

qi_a_lev = qi_a(:,:,plot_level);
ql_a_lev = ql_a(:,:,plot_level);
qils_a_lev = qils_a(:,:,plot_level);
qlls_a_lev = qlls_a(:,:,plot_level);
qicn_a_lev = qicn_a(:,:,plot_level);
qlcn_a_lev = qlcn_a(:,:,plot_level);
cfls_a_lev = cfls_a(:,:,plot_level);
cfcn_a_lev = cfcn_a(:,:,plot_level);

qit = qils_a_lev + qicn_a_lev;
qlt = qlls_a_lev + qlcn_a_lev;

ils_frac = qils_a_lev./(qils_a_lev + qicn_a_lev);

qit_ana = rand(im,jm).*qit./10;

qils = qit.*ils_frac;

qils_ana = qit_ana.*ils_frac;


contour(lon,lat,qit'-qi_a_lev')
colorbar

figure
contour(lon,lat,qils_ana')
colorbar

figure
contour(lon,lat,qit')
colorbar

figure
contour(lon,lat,qils')
colorbar
