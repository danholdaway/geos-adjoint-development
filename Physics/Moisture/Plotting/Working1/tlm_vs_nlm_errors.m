close all
clear
clc

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
 file_end = '6_0600';
% file_end = '6_1200';
% file_end = '6_1500';
% file_end = '7_0000';

cd /home/drholdaw/Lin_Moist_Physics/Inputs/
pref
 
%Load Free (background) State.
cd /discover/nobackup/drholdaw/tmp.22292/nlm_runs/free/
file = ['x0011dh_a.prog.eta.2013011',file_end(1:4),'z.nc4'];

lon = ncread(file,'lon');
lat = ncread(file,'lat');
lev = ncread(file,'lev');

u_free = ncread(file,'u');
v_free = ncread(file,'v');
t_free = ncread(file,'tv');
q_free = ncread(file,'sphu');
p_free = ncread(file,'delp');
qi_free = ncread(file,'qitot');
ql_free = ncread(file,'qltot');
o3_free = ncread(file,'ozone');

[im, jm, lm] = size(u_free);

%Load Perturbed (analysis) state.
cd /discover/nobackup/drholdaw/tmp.22292/nlm_runs/replay10/
file = ['x0011dh_a.prog.eta.2013011',file_end(1:4),'z.nc4'];

u_replay = ncread(file,'u');
v_replay = ncread(file,'v');
t_replay = ncread(file,'tv');
q_replay = ncread(file,'sphu');
p_replay = ncread(file,'delp');
qi_replay = ncread(file,'qitot');
ql_replay = ncread(file,'qltot');
o3_replay = ncread(file,'ozone');

%Load first TLM state to compare.
cd /discover/nobackup/drholdaw/tmp.22292/sens.20130117.000000%/fvpert_old/
file = ['x0011dh_a.fvpert.eta.2013011',file_end,'z_dry.nc4'];

u_tlmdry = ncread(file,'u');
v_tlmdry = ncread(file,'v');
t_tlmdry = ncread(file,'tv');
q_tlmdry = ncread(file,'sphu');
p_tlmdry = ncread(file,'delp');
qi_tlmdry = ncread(file,'qitot');
ql_tlmdry = ncread(file,'qltot');
o3_tlmdry = ncread(file,'ozone');

%Load second TLM state to compare.
cd /discover/nobackup/drholdaw/tmp.22292/sens.20130117.000000%/fvpert_old/
file = ['x0011dh_a.fvpert.eta.2013011',file_end,'z_moi.nc4'];

u_tlmmoi = ncread(file,'u');
v_tlmmoi = ncread(file,'v');
t_tlmmoi = ncread(file,'tv');
q_tlmmoi = ncread(file,'sphu');
p_tlmmoi = ncread(file,'delp');
qi_tlmmoi = ncread(file,'qitot');
ql_tlmmoi = ncread(file,'qltot');
o3_tlmmoi = ncread(file,'ozone');

cd /home/drholdaw/Lin_Moist_Physics/DAS_plotting/Cloud' Fraction'/

%Compute NL perturbation trajectory.
u_nlm = u_replay - u_free;
v_nlm = v_replay - v_free;
t_nlm = t_replay - t_free;
q_nlm = q_replay - q_free;
p_nlm = p_replay - p_free;
qi_nlm = qi_replay - qi_free;
ql_nlm = ql_replay - ql_free;
o3_nlm = o3_replay - o3_free;


fprintf('Computing Errors \n')

u_reler_dry = zeros(lm,1);
v_reler_dry = zeros(lm,1);
t_reler_dry = zeros(lm,1);
q_reler_dry = zeros(lm,1);
p_reler_dry = zeros(lm,1);
qi_reler_dry = zeros(lm,1);
ql_reler_dry = zeros(lm,1);
o3_reler_dry = zeros(lm,1);

u_reler_moi = zeros(lm,1);
v_reler_moi = zeros(lm,1);
t_reler_moi = zeros(lm,1);
q_reler_moi = zeros(lm,1);
p_reler_moi = zeros(lm,1);
qi_reler_moi = zeros(lm,1);
ql_reler_moi = zeros(lm,1);
o3_reler_moi = zeros(lm,1);

u_rmse_dry = zeros(lm,1);
v_rmse_dry = zeros(lm,1);
t_rmse_dry = zeros(lm,1);
q_rmse_dry = zeros(lm,1);
p_rmse_dry = zeros(lm,1);
qi_rmse_dry = zeros(lm,1);
ql_rmse_dry = zeros(lm,1);
o3_rmse_dry = zeros(lm,1);

u_rmse_moi = zeros(lm,1);
v_rmse_moi = zeros(lm,1);
t_rmse_moi = zeros(lm,1);
q_rmse_moi = zeros(lm,1);
p_rmse_moi = zeros(lm,1);
qi_rmse_moi = zeros(lm,1);
ql_rmse_moi = zeros(lm,1);
o3_rmse_moi = zeros(lm,1);

for k = 1:lm
        
    %Pick out just the level to be plotted.
    u_a_lev = u_nlm(:,:,k);
    v_a_lev = v_nlm(:,:,k);
    t_a_lev = t_nlm(:,:,k);
    q_a_lev = q_nlm(:,:,k);
    p_a_lev = p_nlm(:,:,k);
    qi_a_lev = qi_nlm(:,:,k);
    ql_a_lev = ql_nlm(:,:,k);
    o3_a_lev = o3_nlm(:,:,k);
    u_a_lev = u_a_lev(:);
    v_a_lev = v_a_lev(:);
    t_a_lev = t_a_lev(:);
    q_a_lev = q_a_lev(:);
    p_a_lev = p_a_lev(:);
    qi_a_lev = qi_a_lev(:);
    ql_a_lev = ql_a_lev(:);
    o3_a_lev = o3_a_lev(:);

    u_b_lev = u_tlmdry(:,:,k);
    v_b_lev = v_tlmdry(:,:,k);
    t_b_lev = t_tlmdry(:,:,k);
    q_b_lev = q_tlmdry(:,:,k);
    p_b_lev = p_tlmdry(:,:,k);
    qi_b_lev = qi_tlmdry(:,:,k);
    ql_b_lev = ql_tlmdry(:,:,k);
    o3_b_lev = o3_tlmdry(:,:,k);
    u_b_lev = u_b_lev(:);
    v_b_lev = v_b_lev(:);
    t_b_lev = t_b_lev(:);
    q_b_lev = q_b_lev(:);
    p_b_lev = p_b_lev(:);
    qi_b_lev = qi_b_lev(:);
    ql_b_lev = ql_b_lev(:);
    o3_b_lev = o3_b_lev(:);

    u_c_lev = u_tlmmoi(:,:,k);
    v_c_lev = v_tlmmoi(:,:,k);
    t_c_lev = t_tlmmoi(:,:,k);
    q_c_lev = q_tlmmoi(:,:,k);
    p_c_lev = p_tlmmoi(:,:,k);
    qi_c_lev = qi_tlmmoi(:,:,k);
    ql_c_lev = ql_tlmmoi(:,:,k);
    o3_c_lev = o3_tlmmoi(:,:,k);
    u_c_lev = u_c_lev(:);
    v_c_lev = v_c_lev(:);
    t_c_lev = t_c_lev(:);
    q_c_lev = q_c_lev(:);
    p_c_lev = p_c_lev(:);
    qi_c_lev = qi_c_lev(:);
    ql_c_lev = ql_c_lev(:);
    o3_c_lev = o3_c_lev(:);
    
    
    u_reler_dry(k) = mean(abs(u_b_lev - u_a_lev)./abs(u_a_lev));
    v_reler_dry(k) = mean(abs(v_b_lev - v_a_lev)./abs(v_a_lev));
    t_reler_dry(k) = mean(abs(t_b_lev - t_a_lev)./abs(t_a_lev));
    q_reler_dry(k) = mean(abs(q_b_lev - q_a_lev)./abs(q_a_lev));
    p_reler_dry(k) = mean(abs(p_b_lev - p_a_lev)./abs(p_a_lev));
    qi_reler_dry(k) = mean(abs(qi_b_lev - qi_a_lev)./abs(qi_a_lev));
    ql_reler_dry(k) = mean(abs(ql_b_lev - ql_a_lev)./abs(ql_a_lev));
    o3_reler_dry(k) = mean(abs(o3_b_lev - o3_a_lev)./abs(o3_a_lev));
    
    u_reler_moi(k) = mean(abs(u_c_lev - u_a_lev)./abs(u_a_lev));
    v_reler_moi(k) = mean(abs(v_c_lev - v_a_lev)./abs(v_a_lev));
    t_reler_moi(k) = mean(abs(t_c_lev - t_a_lev)./abs(t_a_lev));
    q_reler_moi(k) = mean(abs(q_c_lev - q_a_lev)./abs(q_a_lev));
    p_reler_moi(k) = mean(abs(p_c_lev - p_a_lev)./abs(p_a_lev));
    qi_reler_moi(k) = mean(abs(qi_c_lev - qi_a_lev)./abs(qi_a_lev));
    ql_reler_moi(k) = mean(abs(ql_c_lev - ql_a_lev)./abs(ql_a_lev));
    o3_reler_moi(k) = mean(abs(o3_c_lev - o3_a_lev)./abs(o3_a_lev));
    
    u_rmse_dry(k) = sqrt(mean((u_b_lev - u_a_lev).^2));
    v_rmse_dry(k) = sqrt(mean((v_b_lev - v_a_lev).^2));
    t_rmse_dry(k) = sqrt(mean((t_b_lev - t_a_lev).^2));
    q_rmse_dry(k) = sqrt(mean((q_b_lev - q_a_lev).^2));
    p_rmse_dry(k) = sqrt(mean((p_b_lev - p_a_lev).^2));
    qi_rmse_dry(k) = sqrt(mean((qi_b_lev - qi_a_lev).^2));
    ql_rmse_dry(k) = sqrt(mean((ql_b_lev - ql_a_lev).^2));
    o3_rmse_dry(k) = sqrt(mean((o3_b_lev - o3_a_lev).^2));
    
    u_rmse_moi(k) = sqrt(mean((u_c_lev - u_a_lev).^2));
    v_rmse_moi(k) = sqrt(mean((v_c_lev - v_a_lev).^2));
    t_rmse_moi(k) = sqrt(mean((t_c_lev - t_a_lev).^2));
    q_rmse_moi(k) = sqrt(mean((q_c_lev - q_a_lev).^2));
    p_rmse_moi(k) = sqrt(mean((p_c_lev - p_a_lev).^2));
    qi_rmse_moi(k) = sqrt(mean((qi_c_lev - qi_a_lev).^2));
    ql_rmse_moi(k) = sqrt(mean((ql_c_lev - ql_a_lev).^2));
    o3_rmse_moi(k) = sqrt(mean((o3_c_lev - o3_a_lev).^2));
    
end

lev_min = 40;
lev_max = 72;

mean_u_rmse = [ mean(u_rmse_dry(lev_min:lev_max,1)) mean(u_rmse_moi(lev_min:lev_max,1)) ];
fprintf('Mean u  RMSE. Dry: %f  Moist: %f  Rel-reduction: %f \n\n', mean_u_rmse(1), mean_u_rmse(2), (mean_u_rmse(1)-mean_u_rmse(2))/mean_u_rmse(2));
mean_v_rmse = [ mean(v_rmse_dry(lev_min:lev_max,1)) mean(v_rmse_moi(lev_min:lev_max,1)) ];
fprintf('Mean v  RMSE. Dry: %f  Moist: %f  Rel-reduction: %f \n\n', mean_v_rmse(1), mean_v_rmse(2), (mean_v_rmse(1)-mean_v_rmse(2))/mean_v_rmse(2));
mean_t_rmse = [ mean(t_rmse_dry(lev_min:lev_max,1)) mean(t_rmse_moi(lev_min:lev_max,1)) ];
fprintf('Mean t  RMSE. Dry: %f  Moist: %f  Rel-reduction: %f \n\n', mean_t_rmse(1), mean_t_rmse(2), (mean_t_rmse(1)-mean_t_rmse(2))/mean_t_rmse(2));
mean_q_rmse = [ mean(q_rmse_dry(lev_min:lev_max,1)) mean(q_rmse_moi(lev_min:lev_max,1)) ];
fprintf('Mean q  RMSE. Dry: %f  Moist: %f  Rel-reduction: %f \n\n', mean_q_rmse(1), mean_q_rmse(2), (mean_q_rmse(1)-mean_q_rmse(2))/mean_q_rmse(2));
mean_p_rmse = [ mean(p_rmse_dry(lev_min:lev_max,1)) mean(p_rmse_moi(lev_min:lev_max,1)) ];
fprintf('Mean p  RMSE. Dry: %f  Moist: %f  Rel-reduction: %f \n\n', mean_p_rmse(1), mean_p_rmse(2), (mean_p_rmse(1)-mean_p_rmse(2))/mean_p_rmse(2));
mean_qi_rmse = [ mean(qi_rmse_dry(lev_min:lev_max,1)) mean(qi_rmse_moi(lev_min:lev_max,1)) ];
fprintf('Mean qi  RMSE. Dry: %f  Moist: %f  Rel-reduction: %f \n\n', mean_qi_rmse(1), mean_qi_rmse(2), (mean_qi_rmse(1)-mean_qi_rmse(2))/mean_qi_rmse(2));
mean_ql_rmse = [ mean(ql_rmse_dry(lev_min:lev_max,1)) mean(ql_rmse_moi(lev_min:lev_max,1)) ];
fprintf('Mean ql  RMSE. Dry: %f  Moist: %f  Rel-reduction: %f \n\n', mean_ql_rmse(1), mean_ql_rmse(2), (mean_ql_rmse(1)-mean_ql_rmse(2))/mean_ql_rmse(2));
mean_o3_rmse = [ mean(o3_rmse_dry(lev_min:lev_max,1)) mean(o3_rmse_moi(lev_min:lev_max,1)) ];
fprintf('Mean o3  RMSE. Dry: %f  Moist: %f  Rel-reduction: %f \n\n', mean_o3_rmse(1), mean_o3_rmse(2), (mean_o3_rmse(1)-mean_o3_rmse(2))/mean_o3_rmse(2));


figure
set(gcf,'position',[3 343 1276 576])

subplot(1,6,1)
plot(u_rmse_dry(lev_min:lev_max),p_ref(lev_min+1:lev_max+1))
hold on
plot(mean(u_rmse_dry(lev_min:lev_max))*ones(1,lev_max-lev_min+1),p_ref(lev_min+1:lev_max+1),'--')
plot(u_rmse_moi(lev_min:lev_max),p_ref(lev_min+1:lev_max+1),'r')
plot(mean(u_rmse_moi(lev_min:lev_max))*ones(1,lev_max-lev_min+1),p_ref(lev_min+1:lev_max+1),'r--')
set(gca,'YDir','reverse')
ylim([p_ref(lev_min+1) p_ref(lev_max+1)])

subplot(1,6,2)
plot(v_rmse_dry(lev_min:lev_max),p_ref(lev_min+1:lev_max+1))
hold on
plot(mean(v_rmse_dry(lev_min:lev_max))*ones(1,lev_max-lev_min+1),p_ref(lev_min+1:lev_max+1),'--')
plot(v_rmse_moi(lev_min:lev_max),p_ref(lev_min+1:lev_max+1),'r')
plot(mean(v_rmse_moi(lev_min:lev_max))*ones(1,lev_max-lev_min+1),p_ref(lev_min+1:lev_max+1),'r--')
set(gca,'YDir','reverse')
ylim([p_ref(lev_min+1) p_ref(lev_max+1)])

subplot(1,6,3)
plot(t_rmse_dry(lev_min:lev_max),p_ref(lev_min+1:lev_max+1))
hold on
plot(mean(t_rmse_dry(lev_min:lev_max))*ones(1,lev_max-lev_min+1),p_ref(lev_min+1:lev_max+1),'--')
plot(t_rmse_moi(lev_min:lev_max),p_ref(lev_min+1:lev_max+1),'r')
plot(mean(t_rmse_moi(lev_min:lev_max))*ones(1,lev_max-lev_min+1),p_ref(lev_min+1:lev_max+1),'r--')
set(gca,'YDir','reverse')
ylim([p_ref(lev_min+1) p_ref(lev_max+1)])

subplot(1,6,4)
plot(q_rmse_dry(lev_min:lev_max),p_ref(lev_min+1:lev_max+1))
hold on
plot(mean(q_rmse_dry(lev_min:lev_max))*ones(1,lev_max-lev_min+1),p_ref(lev_min+1:lev_max+1),'--')
plot(q_rmse_moi(lev_min:lev_max),p_ref(lev_min+1:lev_max+1),'r')
plot(mean(q_rmse_moi(lev_min:lev_max))*ones(1,lev_max-lev_min+1),p_ref(lev_min+1:lev_max+1),'r--')
set(gca,'YDir','reverse')
ylim([p_ref(lev_min+1) p_ref(lev_max+1)])

subplot(1,6,5)
plot(qi_rmse_dry(lev_min:lev_max),p_ref(lev_min+1:lev_max+1))
hold on
plot(mean(qi_rmse_dry(lev_min:lev_max))*ones(1,lev_max-lev_min+1),p_ref(lev_min+1:lev_max+1),'--')
plot(qi_rmse_moi(lev_min:lev_max),p_ref(lev_min+1:lev_max+1),'r')
plot(mean(qi_rmse_moi(lev_min:lev_max))*ones(1,lev_max-lev_min+1),p_ref(lev_min+1:lev_max+1),'r--')
set(gca,'YDir','reverse')
ylim([p_ref(lev_min+1) p_ref(lev_max+1)])

subplot(1,6,6)
plot(ql_rmse_dry(lev_min:lev_max),p_ref(lev_min+1:lev_max+1))
hold on
plot(mean(ql_rmse_dry(lev_min:lev_max))*ones(1,lev_max-lev_min+1),p_ref(lev_min+1:lev_max+1),'--')
plot(ql_rmse_moi(lev_min:lev_max),p_ref(lev_min+1:lev_max+1),'r')
plot(mean(ql_rmse_moi(lev_min:lev_max))*ones(1,lev_max-lev_min+1),p_ref(lev_min+1:lev_max+1),'r--')
set(gca,'YDir','reverse')
ylim([p_ref(lev_min+1) p_ref(lev_max+1)])
