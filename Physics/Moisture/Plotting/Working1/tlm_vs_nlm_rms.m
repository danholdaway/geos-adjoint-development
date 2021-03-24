close all
clear
clc
mydir = pwd;


%Choose whether to compute the correlations.
calc_rms = 1;

%Choose Region,
% 1 = Global
% 2 = Tropics 23S to 23N, below 100hPa
% 3 = Northern hemisphere
% 4 = Southern hemisphere
region = 1;

if region == 1 
    fprintf('Computing RMSEs globally \n')
elseif region == 2
    fprintf('Computing RMSEs for the tropics \n')  
elseif region == 3
    fprintf('Computing RMSEs for the northern hemisphere \n')
elseif region == 4
    fprintf('Computing RMSEs for the southern hemisphere \n')
else
    fprintf('Not a valid region selection \n')
end

%Choose the time being considered: 16th at 12z would be 6_12
% file_end = '1_0100';
% file_end = '1_0300';
% file_end = '1_0600';
% file_end = '1_0900';
% file_end = '1_1200';
% file_end = '1_1500';
% file_end = '1_1800';
% file_end = '1_2100';
% file_end = '2_0000';
file_end = '2_0300';

% file_end = '2_0300'; % 24 hours


im = 576;
jm = 361;
lm = 72;

LevMin = 1;
LevMax = lm;

if region == 1
    LonMin = 1;
    LonMax = im;
    LatMin = 1;
    LatMax = jm;
elseif region == 2
    LonMin = 1;
    LonMax = im;
    LatMin = 135;
    LatMax = 227;
elseif region == 3
    LonMin = 1;
    LonMax = im;
    LatMin = ceil(jm/2);
    LatMax = jm;
elseif region == 4
    LonMin = 1;
    LonMax = im;
    LatMin = 1;
    LatMax = floor(jm/2);
end

cd /home/drholdaw/LinearisedPhysics/Inputs/
pref

%Load Free (background) State.
cd /discover/nobackup/drholdaw/ExperimentData/TracerDryRemap/onlytracer_hord1/
file = ['v000_C180.prog_free.eta.2014020',file_end,'z.nc4'];

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

%Load Perturbed (analysis) st
cd /discover/nobackup/drholdaw/ExperimentData/TracerDryRemap/onlytracer_hord1/
file = ['v000_C180.prog_replay.eta.2014020',file_end,'z.nc4'];

u_replay = ncread(file,'u');
v_replay = ncread(file,'v');
t_replay = ncread(file,'tv');
q_replay = ncread(file,'sphu');
p_replay = ncread(file,'delp');
qi_replay = ncread(file,'qitot');
ql_replay = ncread(file,'qltot');
o3_replay = ncread(file,'ozone');

%Load first TLM state to compare.
cd /discover/nobackup/drholdaw/btmp.25956/sens.20140202.000000/
file = ['v000_C180.fvpert.eta.2014020',file_end,'z_ST03z_P20000_S0_onlytracer_hord1.nc4'];

u_tlma = ncread(file,'U');
v_tlma = ncread(file,'V');
t_tlma = ncread(file,'TV');
q_tlma = ncread(file,'QV');
p_tlma = ncread(file,'DP');
qi_tlma = ncread(file,'QI');
ql_tlma = ncread(file,'QL');
o3_tlma = ncread(file,'O3');

%Load second TLM state to compare.
cd /discover/nobackup/drholdaw/btmp.25956/sens.20140202.000000/
file = ['v000_C180.fvpert.eta.2014020',file_end,'z_ST03z_P20000_S1_onlytracer_hord1.nc4'];

u_tlmb = ncread(file,'U');
v_tlmb = ncread(file,'V');
t_tlmb = ncread(file,'TV');
q_tlmb = ncread(file,'QV');
p_tlmb = ncread(file,'DP');
qi_tlmb = ncread(file,'QI');
ql_tlmb = ncread(file,'QL');
o3_tlmb = ncread(file,'O3');

cd(mydir)

%Compute NL perturbation trajectory.
u_nlm = u_replay - u_free;
v_nlm = v_replay - v_free;
t_nlm = t_replay - t_free;
q_nlm = q_replay - q_free;
p_nlm = p_replay - p_free;
qi_nlm = qi_replay - qi_free;
ql_nlm = ql_replay - ql_free;
o3_nlm = o3_replay - o3_free;


im = length(lon);
jm = length(lat);
lm = length(lev);

if region == 1
    lon_min = 1;
    lon_max = im;
    lat_min = 1;
    lat_max = jm;
elseif region == 2
    lon_min = 1;
    lon_max = im;
    lat_min = 135;
    lat_max = 227;
elseif region == 3
    lon_min = 1;
    lon_max = im;
    lat_min = ceil(jm/2);
    lat_max = jm;
elseif region == 4
    lon_min = 1;
    lon_max = im;
    lat_min = 1;
    lat_max = floor(jm/2);
   
end

lev_min = 1;
lev_max = lm;    

rms_u_a = zeros(lm,1);
rms_u_b = zeros(lm,1);
rms_u_c = zeros(lm,1);
rms_v_a = zeros(lm,1);
rms_v_b = zeros(lm,1);
rms_v_c = zeros(lm,1);
rms_t_a = zeros(lm,1);
rms_t_b = zeros(lm,1);
rms_t_c = zeros(lm,1);
rms_q_a = zeros(lm,1);
rms_q_b = zeros(lm,1);
rms_q_c = zeros(lm,1);
rms_p_a = zeros(lm,1);
rms_p_b = zeros(lm,1);
rms_p_c = zeros(lm,1);
rms_qi_a = zeros(lm,1);
rms_qi_b = zeros(lm,1);
rms_qi_c = zeros(lm,1);
rms_ql_a = zeros(lm,1);
rms_ql_b = zeros(lm,1);
rms_ql_c = zeros(lm,1);
rms_o3_a = zeros(lm,1);
rms_o3_b = zeros(lm,1);
rms_o3_c = zeros(lm,1);

n = (lon_max-lon_min+1)*(lat_max-lat_min+1);

%Loop over model all levels.
for k = 1:lm

    rms_u_a(k) = sqrt(  sum(sum(u_nlm(lon_min:lon_max,lat_min:lat_max,k).^2)) / n );
    rms_u_b(k) = sqrt(  sum(sum(u_tlma(lon_min:lon_max,lat_min:lat_max,k).^2)) / n );
    rms_u_c(k) = sqrt(  sum(sum(u_tlmb(lon_min:lon_max,lat_min:lat_max,k).^2)) / n );

    rms_v_a(k) = sqrt(  sum(sum(v_nlm(lon_min:lon_max,lat_min:lat_max,k).^2)) / n );
    rms_v_b(k) = sqrt(  sum(sum(v_tlma(lon_min:lon_max,lat_min:lat_max,k).^2)) / n );
    rms_v_c(k) = sqrt(  sum(sum(v_tlmb(lon_min:lon_max,lat_min:lat_max,k).^2)) / n );

    rms_t_a(k) = sqrt(  sum(sum(t_nlm(lon_min:lon_max,lat_min:lat_max,k).^2)) / n );
    rms_t_b(k) = sqrt(  sum(sum(t_tlma(lon_min:lon_max,lat_min:lat_max,k).^2)) / n );
    rms_t_c(k) = sqrt(  sum(sum(t_tlmb(lon_min:lon_max,lat_min:lat_max,k).^2)) / n );

    rms_q_a(k) = sqrt(  sum(sum(q_nlm(lon_min:lon_max,lat_min:lat_max,k).^2)) / n );
    rms_q_b(k) = sqrt(  sum(sum(q_tlma(lon_min:lon_max,lat_min:lat_max,k).^2)) / n );
    rms_q_c(k) = sqrt(  sum(sum(q_tlmb(lon_min:lon_max,lat_min:lat_max,k).^2)) / n );

    rms_p_a(k) = sqrt(  sum(sum(p_nlm(lon_min:lon_max,lat_min:lat_max,k).^2)) / n );
    rms_p_b(k) = sqrt(  sum(sum(p_tlma(lon_min:lon_max,lat_min:lat_max,k).^2)) / n );
    rms_p_c(k) = sqrt(  sum(sum(p_tlmb(lon_min:lon_max,lat_min:lat_max,k).^2)) / n );

    rms_qi_a(k) = sqrt(  sum(sum(qi_nlm(lon_min:lon_max,lat_min:lat_max,k).^2)) / n );
    rms_qi_b(k) = sqrt(  sum(sum(qi_tlma(lon_min:lon_max,lat_min:lat_max,k).^2)) / n );
    rms_qi_c(k) = sqrt(  sum(sum(qi_tlmb(lon_min:lon_max,lat_min:lat_max,k).^2)) / n );

    rms_ql_a(k) = sqrt(  sum(sum(ql_nlm(lon_min:lon_max,lat_min:lat_max,k).^2)) / n );
    rms_ql_b(k) = sqrt(  sum(sum(ql_tlma(lon_min:lon_max,lat_min:lat_max,k).^2)) / n );
    rms_ql_c(k) = sqrt(  sum(sum(ql_tlmb(lon_min:lon_max,lat_min:lat_max,k).^2)) / n );

    rms_o3_a(k) = sqrt(  sum(sum(o3_nlm(lon_min:lon_max,lat_min:lat_max,k).^2)) / n );
    rms_o3_b(k) = sqrt(  sum(sum(o3_tlma(lon_min:lon_max,lat_min:lat_max,k).^2)) / n );
    rms_o3_c(k) = sqrt(  sum(sum(o3_tlmb(lon_min:lon_max,lat_min:lat_max,k).^2)) / n );
end

format short


figure
set(gcf,'position',[3 343 1276 576])

subplot(1,8,1)
plot(rms_u_a(lev_min:lev_max),p_ref(lev_min+1:lev_max+1),'k')
hold on
plot(rms_u_b(lev_min:lev_max),p_ref(lev_min+1:lev_max+1),'g')
plot(rms_u_c(lev_min:lev_max),p_ref(lev_min+1:lev_max+1),'r')
set(gca,'YDir','reverse')
ylim([p_ref(30+1) p_ref(lev_max+1)])

subplot(1,8,2)
plot(rms_v_a(lev_min:lev_max),p_ref(lev_min+1:lev_max+1),'k')
hold on
plot(rms_v_b(lev_min:lev_max),p_ref(lev_min+1:lev_max+1),'g')
plot(rms_v_c(lev_min:lev_max),p_ref(lev_min+1:lev_max+1),'r')
set(gca,'YDir','reverse')
ylim([p_ref(30+1) p_ref(lev_max+1)])

subplot(1,8,3)
plot(rms_t_a(lev_min:lev_max),p_ref(lev_min+1:lev_max+1),'k')
hold on
plot(rms_t_b(lev_min:lev_max),p_ref(lev_min+1:lev_max+1),'g')
plot(rms_t_c(lev_min:lev_max),p_ref(lev_min+1:lev_max+1),'r')
set(gca,'YDir','reverse')
ylim([p_ref(30+1) p_ref(lev_max+1)])

subplot(1,8,4)
plot(rms_q_a(lev_min:lev_max),p_ref(lev_min+1:lev_max+1),'k')
hold on
plot(rms_q_b(lev_min:lev_max),p_ref(lev_min+1:lev_max+1),'g')
plot(rms_q_c(lev_min:lev_max),p_ref(lev_min+1:lev_max+1),'r')
set(gca,'YDir','reverse')
ylim([p_ref(30+1) p_ref(lev_max+1)])

subplot(1,8,5)
plot(rms_p_a(lev_min:lev_max),p_ref(lev_min+1:lev_max+1),'k')
hold on
plot(rms_p_b(lev_min:lev_max),p_ref(lev_min+1:lev_max+1),'g')
plot(rms_p_c(lev_min:lev_max),p_ref(lev_min+1:lev_max+1),'r')
set(gca,'YDir','reverse')
ylim([p_ref(30+1) p_ref(lev_max+1)])

subplot(1,8,6)
plot(rms_qi_a(lev_min:lev_max),p_ref(lev_min+1:lev_max+1),'k')
hold on
plot(rms_qi_b(lev_min:lev_max),p_ref(lev_min+1:lev_max+1),'g')
plot(rms_qi_c(lev_min:lev_max),p_ref(lev_min+1:lev_max+1),'r')
set(gca,'YDir','reverse')
ylim([p_ref(30+1) p_ref(lev_max+1)])

subplot(1,8,7)
plot(rms_ql_a(lev_min:lev_max),p_ref(lev_min+1:lev_max+1),'k')
hold on
plot(rms_ql_b(lev_min:lev_max),p_ref(lev_min+1:lev_max+1),'g')
plot(rms_ql_c(lev_min:lev_max),p_ref(lev_min+1:lev_max+1),'r')
set(gca,'YDir','reverse')
ylim([p_ref(30+1) p_ref(lev_max+1)])

subplot(1,8,8)
plot(rms_o3_a(lev_min:lev_max),p_ref(lev_min+1:lev_max+1),'k')
hold on
plot(rms_o3_b(lev_min:lev_max),p_ref(lev_min+1:lev_max+1),'g')
plot(rms_o3_c(lev_min:lev_max),p_ref(lev_min+1:lev_max+1),'r')
set(gca,'YDir','reverse')
ylim([p_ref(30+1) p_ref(lev_max+1)])




figure
set(gcf,'position',[3 343 1276 576])

subplot(1,8,1)
plot(rms_u_a(lev_min:lev_max),lev_min:lev_max,'k')
hold on
plot(rms_u_b(lev_min:lev_max),lev_min:lev_max,'g')
plot(rms_u_c(lev_min:lev_max),lev_min:lev_max,'r')
set(gca,'YDir','reverse')
ylim([lev_min lev_max])

subplot(1,8,2)
plot(rms_v_a(lev_min:lev_max),lev_min:lev_max,'k')
hold on
plot(rms_v_b(lev_min:lev_max),lev_min:lev_max,'g')
plot(rms_v_c(lev_min:lev_max),lev_min:lev_max,'r')
set(gca,'YDir','reverse')
ylim([lev_min lev_max])

subplot(1,8,3)
plot(rms_t_a(lev_min:lev_max),lev_min:lev_max,'k')
hold on
plot(rms_t_b(lev_min:lev_max),lev_min:lev_max,'g')
plot(rms_t_c(lev_min:lev_max),lev_min:lev_max,'r')
set(gca,'YDir','reverse')
ylim([lev_min lev_max])

subplot(1,8,4)
plot(rms_q_a(lev_min:lev_max),lev_min:lev_max,'k')
hold on
plot(rms_q_b(lev_min:lev_max),lev_min:lev_max,'g')
plot(rms_q_c(lev_min:lev_max),lev_min:lev_max,'r')
set(gca,'YDir','reverse')
ylim([lev_min lev_max])

subplot(1,8,5)
plot(rms_p_a(lev_min:lev_max),lev_min:lev_max,'k')
hold on
plot(rms_p_b(lev_min:lev_max),lev_min:lev_max,'g')
plot(rms_p_c(lev_min:lev_max),lev_min:lev_max,'r')
set(gca,'YDir','reverse')
ylim([lev_min lev_max])

subplot(1,8,6)
plot(rms_qi_a(lev_min:lev_max),lev_min:lev_max,'k')
hold on
plot(rms_qi_b(lev_min:lev_max),lev_min:lev_max,'g')
plot(rms_qi_c(lev_min:lev_max),lev_min:lev_max,'r')
set(gca,'YDir','reverse')
ylim([lev_min lev_max])

subplot(1,8,7)
plot(rms_ql_a(lev_min:lev_max),lev_min:lev_max,'k')
hold on
plot(rms_ql_b(lev_min:lev_max),lev_min:lev_max,'g')
plot(rms_ql_c(lev_min:lev_max),lev_min:lev_max,'r')
set(gca,'YDir','reverse')
ylim([lev_min lev_max])

subplot(1,8,8)
plot(rms_o3_a(lev_min:lev_max),lev_min:lev_max,'k')
hold on
plot(rms_o3_b(lev_min:lev_max),lev_min:lev_max,'g')
plot(rms_o3_c(lev_min:lev_max),lev_min:lev_max,'r')
set(gca,'YDir','reverse')
ylim([lev_min lev_max])
