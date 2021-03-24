close all
clear
clc

%Choose model level to plot.
plot_level = 63;

fontsize = 13;

line_wid_cont = 1.0;
line_wid_det = 0.6;

linewid = 1.5;

grey = 0.75;

%Choose whether to compute the correlations.
calc_corr = 1;

%Choose Region,
region = 1;

%Choose whether you want plots of pert traj (0/1) and which model level to plot 
makeplots = 0;
plot_level = 50;

plot_u = 0;
plot_v = 0;
plot_t = 1;
plot_q = 0;
plot_p = 0;
plot_qi = 0;
plot_ql = 0;
plot_o3 = 0;

%Choose whether to compute the correlations.
calc_corr = 1;

%Choose Region,
% 1 = Global
% 2 = Tropics 30S to 30N, below 100hPa
% 3 = Northern hemisphere
% 4 = Southern hemisphere
region = 1;

%Choose the time being considered: 16th at 12z would be 6_12
% file_end = '6_0300';
% file_end = '6_0600';
file_end = '6_0900';
% file_end = '6_1200';
% file_end = '6_1500';
% file_end = '6_1800';
% file_end = '6_2100';
% file_end = '7_0000';

cd /home/drholdaw/Lin_Moist_Physics/Inputs/
pref
 
%Load Free (background) State.
cd /discover/nobackup/drholdaw/tmp.22292/prog/prog_free/
file = ['x0011dh_a.prog.eta.2013011',file_end(1:4),'z.nc4'];

lon = ncread(file,'lon');
lat = ncread(file,'lat');
lev = ncread(file,'lev');

q_free = ncread(file,'sphu');
qi_free = ncread(file,'qitot');
ql_free = ncread(file,'qltot');

%Load Perturbed (analysis) state.
cd /discover/nobackup/drholdaw/tmp.22292/prog/prog_replay10/
file = ['x0011dh_a.prog.eta.2013011',file_end(1:4),'z.nc4'];

q_replay = ncread(file,'sphu');
qi_replay = ncread(file,'qitot');
ql_replay = ncread(file,'qltot');

%Load first TLM state to compare.
cd /discover/nobackup/drholdaw/tmp.22292/sens.20130117.000000%/fvpert_old/
file = ['x0011dh_a.fvpert.eta.2013011',file_end,'z_newBL_03rst.nc4'];

q_tlmdry = ncread(file,'sphu');
qi_tlmdry = ncread(file,'qitot');
ql_tlmdry = ncread(file,'qltot');

%Load second TLM state to compare.
cd /discover/nobackup/drholdaw/tmp.22292/sens.20130117.000000%/fvpert_old/
file = ['x0011dh_a.fvpert.eta.2013011',file_end,'z_newMOInewBL_03rst.nc4'];

q_tlmmoi = ncread(file,'sphu');
qi_tlmmoi = ncread(file,'qitot');
ql_tlmmoi = ncread(file,'qltot');

cd /home/drholdaw/Lin_Moist_Physics/DAS_plotting/GSI' Talk'/

%Compute NL perturbation trajectory.
q_nlm = q_replay - q_free;
qi_nlm = qi_replay - qi_free;
ql_nlm = ql_replay - ql_free;


q_error_dry = abs(q_nlm - q_tlmdry);
q_error_moi = abs(q_nlm - q_tlmmoi);

q_error_dry = abs(q_nlm - q_tlmdry);
q_error_moi = abs(q_nlm - q_tlmmoi);



[IM,JM,LM] = size(q_nlm);

q_meanerror_dry = zeros(JM,LM);
q_meanerror_moi = zeros(JM,LM);

q_meanerror_dry(:,:) = mean(q_error_dry,1);
q_meanerror_moi(:,:) = mean(q_error_moi,1);

maxerror = max([max(q_meanerror_dry(:)) max(q_meanerror_moi(:))]);

figure
subplot(1,2,1)
contourf(q_meanerror_dry')
caxis([0 maxerror])
set(gca,'YDir','reverse')
ylim([40 72])

subplot(1,2,2)
contourf(q_meanerror_moi')
caxis([0 maxerror])
set(gca,'YDir','reverse')
ylim([40 72])

