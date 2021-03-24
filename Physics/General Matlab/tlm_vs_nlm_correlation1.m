close all
clear
clc
mydir = pwd;

% dir = 'wrk.sandy';
% dir = 'wrk.dust';
dir = 'btmp.25956';

% exp = 'sens_sandy';
% exp = 'C360_513_snamma';
exp = 'v000_C180';

%Set directories of experiment
% dir_free = '/ExperimentData/TracerDryRemap/onlytracer_hord1/';
% dir_repl = '/ExperimentData/TracerDryRemap/onlytracer_hord1/';
dir_free = dir;
dir_repl = dir;
dir_tlm1 = dir;
dir_tlm2 = dir;
dir_tlm3 = dir;

% Set experiment names
exp_free = exp;
exp_repl = exp;
exp_tlm1 = exp;
exp_tlm2 = exp;
exp_tlm3 = exp;

%Choose physics options for TLM
phy_tlm1 = '20000';
phy_tlm2 = '20000';
phy_tlm3 = '20000';

%Choose shapiro filter options for TLM
shp_tlm1 = '0';
shp_tlm2 = '0';
shp_tlm3 = '0';

%Choose 'other' options for TLM
opt_tlm1 = '';
opt_tlm2 = '';
opt_tlm3 = '';

%Choose highest model level
LevMin = 1;

%Choose start date
% datea = '20060827';
% datea = '20060827';
datea = '20140201';

%Choose start time
timea = '030000';

%Choose length of forecasts
lead = 24;

%Plotting options
lin_wid = 1.75;
fontsize = 11;
fontsize1 = 12;

%Choose Region,
% 1 = Global
% 2 = Tropics 23S to 23N, below 100hPa
% 3 = Northern hemisphere
% 4 = Southern hemisphere
region = 1;

if region == 1 
    fprintf(' \nWill compute correlations for the whole globe \n\n')
elseif region == 2
    fprintf('Will compute correlations for the tropics \n\n')  
elseif region == 3
    fprintf('Will compute correlations for the northern hemisphere \n\n')
elseif region == 4
    fprintf('Will compute correlations for the southern hemisphere \n\n')
else
    fprintf('Not a valid region selection \n\n')
end

cd /home/drholdaw/LinearisedPhysics/Inputs/
pref
p = 0.5*(p_ref(1:end-1) + p_ref(2:end));

date_start = datenum(str2double(datea(1:4)), str2double(datea(5:6)), str2double(datea(7:8)), ...
                     str2double(timea(1:2)), str2double(timea(3:4)), str2double(timea(5:6)));

datetimef = datestr(date_start + lead/24     , 'yyyymmddHHMMSS');
datetimev = datestr(date_start + 1           , 'yyyymmddHHMMSS');          

fprintf(' Reading in the states \n')

%Load Free (background) State.
dir = ['/discover/nobackup/drholdaw/',dir_free,'/prog/prog_free/']; cd(dir)
file = [exp_free,'.prog.eta.',datetimef(1:8),'_',datetimef(9:10),'z.nc4'];

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
dir = ['/discover/nobackup/drholdaw/',dir_repl,'/prog/prog_replay/']; cd(dir)
file = [exp_repl,'.prog.eta.',datetimef(1:8),'_',datetimef(9:10),'z.nc4'];

u_replay = ncread(file,'u');
v_replay = ncread(file,'v');
t_replay = ncread(file,'tv');
q_replay = ncread(file,'sphu');
p_replay = ncread(file,'delp');
qi_replay = ncread(file,'qitot');
ql_replay = ncread(file,'qltot');
o3_replay = ncread(file,'ozone');

%Load first Tlm state to compare.
dir = ['/discover/nobackup/drholdaw/',dir_tlm1,'/sens.',datetimev(1:8),'.000000']; cd(dir)
file = [exp_tlm1,'.fvpert.eta.',datetimef(1:8),'_',datetimef(9:12),'z_ST03z_P',phy_tlm1,'_S',shp_tlm1,opt_tlm1,'.nc4'];

fprintf(' TLM1 file: \n')
disp(file)

u_tl1 = ncread(file,'U');
v_tl1 = ncread(file,'V');
t_tl1 = ncread(file,'TV');
q_tl1 = ncread(file,'QV');
p_tl1 = ncread(file,'DP');
qi_tl1 = ncread(file,'QI');
ql_tl1 = ncread(file,'QL');
o3_tl1 = ncread(file,'O3');

%Load second Tlm state to compare.
dir = ['/discover/nobackup/drholdaw/',dir_tlm2,'/sens.',datetimev(1:8),'.000000']; cd(dir)
file = [exp_tlm2,'.fvpert.eta.',datetimef(1:8),'_',datetimef(9:12),'z_ST03z_P',phy_tlm2,'_S',shp_tlm2,opt_tlm2,'.nc4'];

fprintf(' TLM2 file: \n')
disp(file)

u_tl2 = ncread(file,'U');
v_tl2 = ncread(file,'V');
t_tl2 = ncread(file,'TV');
q_tl2 = ncread(file,'QV');
p_tl2 = ncread(file,'DP');
qi_tl2 = ncread(file,'QI');
ql_tl2 = ncread(file,'QL');
o3_tl2 = ncread(file,'O3');

%Load third Tlm state to compare.
dir = ['/discover/nobackup/drholdaw/',dir_tlm3,'/sens.',datetimev(1:8),'.000000']; cd(dir)
file = [exp_tlm3,'.fvpert.eta.',datetimef(1:8),'_',datetimef(9:12),'z_ST03z_P',phy_tlm3,'_S',shp_tlm3,opt_tlm3,'.nc4'];

fprintf(' TLM3 file: \n')
disp(file)

u_tl3 = ncread(file,'U');
v_tl3 = ncread(file,'V');
t_tl3 = ncread(file,'TV');
q_tl3 = ncread(file,'QV');
p_tl3 = ncread(file,'DP');
qi_tl3 = ncread(file,'QI');
ql_tl3 = ncread(file,'QL');
o3_tl3 = ncread(file,'O3');

cd(mydir)

fprintf(' Done reading in the states \n\n')

fprintf('  Correlation computation \n')

im = size(u_free,1);
jm = size(u_free,2);
lm = size(u_free,3);

if region == 1
    LonMin = 1;
    LonMax = im;
    LatMin = 1;
    LatMax = jm;
elseif region == 2
    LonMin = 1;
    LonMax = im;
    LatMin = find(lat == -24);
    LatMax = find(lat ==  24);
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

%Compute NL perturbation trajectory.
u_nlm = u_replay - u_free;
v_nlm = v_replay - v_free;
t_nlm = t_replay - t_free;
q_nlm = q_replay - q_free;
p_nlm = p_replay - p_free;
qi_nlm = qi_replay - qi_free;
ql_nlm = ql_replay - ql_free;
o3_nlm = o3_replay - o3_free;

%Weighting 
A = zeros(im,jm,lm);
A(1,:,1) = cosd(lat);
for i = 1:im
    for k = 1:lm
        A(i,:,k) = A(1,:,1);
    end
end
dA = sum(sum(A(1,:,1)));

%Variance, covariance and correlation
%Second index is number of variables.
var_nlm = zeros(lm-LevMin+1,8);
var_tl1 = zeros(lm-LevMin+1,8);
var_tl2 = zeros(lm-LevMin+1,8);

cov_tl1 = zeros(lm-LevMin+1,8);
cov_tl2 = zeros(lm-LevMin+1,8);

%Initialize correlations
corr_tl1 = zeros(lm-LevMin+1,8);
corr_tl2 = zeros(lm-LevMin+1,8);
corr_tl3 = zeros(lm-LevMin+1,8);

mean_nlm = zeros(im,jm,lm,8);
mean_tl1 = zeros(im,jm,lm,8);
mean_tl2 = zeros(im,jm,lm,8);
mean_tl3 = zeros(im,jm,lm,8);


for k = 1:lm
    mean_nlm(1,1,k,1) = mean(mean(u_nlm(LonMin:LonMax,LatMin:LatMax,k)));
    mean_tl1(1,1,k,1) = mean(mean(u_tl1(LonMin:LonMax,LatMin:LatMax,k)));
    mean_tl2(1,1,k,1) = mean(mean(u_tl2(LonMin:LonMax,LatMin:LatMax,k)));
    mean_tl3(1,1,k,1) = mean(mean(u_tl3(LonMin:LonMax,LatMin:LatMax,k)));
    
    mean_nlm(1,1,k,2) = mean(mean(v_nlm(LonMin:LonMax,LatMin:LatMax,k)));
    mean_tl1(1,1,k,2) = mean(mean(v_tl1(LonMin:LonMax,LatMin:LatMax,k)));
    mean_tl2(1,1,k,2) = mean(mean(v_tl2(LonMin:LonMax,LatMin:LatMax,k)));
    mean_tl3(1,1,k,2) = mean(mean(v_tl3(LonMin:LonMax,LatMin:LatMax,k)));
    
    mean_nlm(1,1,k,3) = mean(mean(t_nlm(LonMin:LonMax,LatMin:LatMax,k)));
    mean_tl1(1,1,k,3) = mean(mean(t_tl1(LonMin:LonMax,LatMin:LatMax,k)));
    mean_tl2(1,1,k,3) = mean(mean(t_tl2(LonMin:LonMax,LatMin:LatMax,k)));
    mean_tl3(1,1,k,3) = mean(mean(t_tl3(LonMin:LonMax,LatMin:LatMax,k)));
    
    mean_nlm(1,1,k,4) = mean(mean(q_nlm(LonMin:LonMax,LatMin:LatMax,k)));
    mean_tl1(1,1,k,4) = mean(mean(q_tl1(LonMin:LonMax,LatMin:LatMax,k)));
    mean_tl2(1,1,k,4) = mean(mean(q_tl2(LonMin:LonMax,LatMin:LatMax,k)));
    mean_tl3(1,1,k,4) = mean(mean(q_tl3(LonMin:LonMax,LatMin:LatMax,k)));
    
    mean_nlm(1,1,k,5) = mean(mean(p_nlm(LonMin:LonMax,LatMin:LatMax,k)));
    mean_tl1(1,1,k,5) = mean(mean(p_tl1(LonMin:LonMax,LatMin:LatMax,k)));
    mean_tl2(1,1,k,5) = mean(mean(p_tl2(LonMin:LonMax,LatMin:LatMax,k)));
    mean_tl3(1,1,k,5) = mean(mean(p_tl3(LonMin:LonMax,LatMin:LatMax,k)));
    
    mean_nlm(1,1,k,6) = mean(mean(qi_nlm(LonMin:LonMax,LatMin:LatMax,k)));
    mean_tl1(1,1,k,6) = mean(mean(qi_tl1(LonMin:LonMax,LatMin:LatMax,k)));
    mean_tl2(1,1,k,6) = mean(mean(qi_tl2(LonMin:LonMax,LatMin:LatMax,k)));
    mean_tl3(1,1,k,6) = mean(mean(qi_tl3(LonMin:LonMax,LatMin:LatMax,k)));
    
    mean_nlm(1,1,k,7) = mean(mean(ql_nlm(LonMin:LonMax,LatMin:LatMax,k)));
    mean_tl1(1,1,k,7) = mean(mean(ql_tl1(LonMin:LonMax,LatMin:LatMax,k)));
    mean_tl2(1,1,k,7) = mean(mean(ql_tl2(LonMin:LonMax,LatMin:LatMax,k)));
    mean_tl3(1,1,k,7) = mean(mean(ql_tl3(LonMin:LonMax,LatMin:LatMax,k)));
    
    mean_nlm(1,1,k,8) = mean(mean(o3_nlm(LonMin:LonMax,LatMin:LatMax,k)));
    mean_tl1(1,1,k,8) = mean(mean(o3_tl1(LonMin:LonMax,LatMin:LatMax,k)));
    mean_tl2(1,1,k,8) = mean(mean(o3_tl2(LonMin:LonMax,LatMin:LatMax,k)));
    mean_tl3(1,1,k,8) = mean(mean(o3_tl3(LonMin:LonMax,LatMin:LatMax,k)));
end

%Fill in array
for l = 1:8
    for i = LonMin+1:LonMax
        for j = LatMin+1:LatMax
            mean_nlm(i,j,:,l) = mean_nlm(1,1,:,l);
            mean_tl1(i,j,:,l) = mean_tl1(1,1,:,l);
            mean_tl2(i,j,:,l) = mean_tl2(1,1,:,l);
            mean_tl3(i,j,:,l) = mean_tl3(1,1,:,l);
        end
    end
end

% U winds
var_nlm(1:lm-LevMin+1,1) = sum(sum((u_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1) - mean_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1,1)).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
var_tl1(1:lm-LevMin+1,1) = sum(sum((u_tl1(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1) - mean_tl1(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1,1)).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
var_tl2(1:lm-LevMin+1,1) = sum(sum((u_tl2(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1) - mean_tl2(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1,1)).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
var_tl3(1:lm-LevMin+1,1) = sum(sum((u_tl3(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1) - mean_tl3(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1,1)).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;

cov_tl1(1:lm-LevMin+1,1) = sum(sum((u_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1) - mean_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1,1)).*(u_tl1(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1) - mean_tl1(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1,1)).*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
cov_tl2(1:lm-LevMin+1,1) = sum(sum((u_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1) - mean_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1,1)).*(u_tl2(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1) - mean_tl2(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1,1)).*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
cov_tl3(1:lm-LevMin+1,1) = sum(sum((u_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1) - mean_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1,1)).*(u_tl3(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1) - mean_tl3(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1,1)).*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;

corr_tl1(1:lm-LevMin+1,1) = cov_tl1(1:lm-LevMin+1,1)./sqrt(var_nlm(1:lm-LevMin+1,1).*var_tl1(1:lm-LevMin+1,1));
corr_tl2(1:lm-LevMin+1,1) = cov_tl2(1:lm-LevMin+1,1)./sqrt(var_nlm(1:lm-LevMin+1,1).*var_tl2(1:lm-LevMin+1,1));
corr_tl3(1:lm-LevMin+1,1) = cov_tl3(1:lm-LevMin+1,1)./sqrt(var_nlm(1:lm-LevMin+1,1).*var_tl3(1:lm-LevMin+1,1));

% V winds
var_nlm(1:lm-LevMin+1,2) = sum(sum((v_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1) - mean_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1,2)).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
var_tl1(1:lm-LevMin+1,2) = sum(sum((v_tl1(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1) - mean_tl1(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1,2)).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
var_tl2(1:lm-LevMin+1,2) = sum(sum((v_tl2(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1) - mean_tl2(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1,2)).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
var_tl3(1:lm-LevMin+1,2) = sum(sum((v_tl3(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1) - mean_tl3(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1,2)).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;

cov_tl1(1:lm-LevMin+1,2) = sum(sum((v_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1) - mean_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1,2)).*(v_tl1(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1) - mean_tl1(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1,2)).*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
cov_tl2(1:lm-LevMin+1,2) = sum(sum((v_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1) - mean_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1,2)).*(v_tl2(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1) - mean_tl2(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1,2)).*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
cov_tl3(1:lm-LevMin+1,2) = sum(sum((v_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1) - mean_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1,2)).*(v_tl3(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1) - mean_tl3(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1,2)).*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;

corr_tl1(1:lm-LevMin+1,2) = cov_tl1(1:lm-LevMin+1,2)./sqrt(var_nlm(1:lm-LevMin+1,2).*var_tl1(1:lm-LevMin+1,2));
corr_tl2(1:lm-LevMin+1,2) = cov_tl2(1:lm-LevMin+1,2)./sqrt(var_nlm(1:lm-LevMin+1,2).*var_tl2(1:lm-LevMin+1,2));
corr_tl3(1:lm-LevMin+1,2) = cov_tl3(1:lm-LevMin+1,2)./sqrt(var_nlm(1:lm-LevMin+1,2).*var_tl3(1:lm-LevMin+1,2));

% Temperature
var_nlm(1:lm-LevMin+1,3) = sum(sum((t_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1) - mean_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1,3)).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
var_tl1(1:lm-LevMin+1,3) = sum(sum((t_tl1(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1) - mean_tl1(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1,3)).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
var_tl2(1:lm-LevMin+1,3) = sum(sum((t_tl2(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1) - mean_tl2(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1,3)).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
var_tl3(1:lm-LevMin+1,3) = sum(sum((t_tl3(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1) - mean_tl3(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1,3)).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;

cov_tl1(1:lm-LevMin+1,3) = sum(sum((t_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1) - mean_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1,3)).*(t_tl1(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1) - mean_tl1(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1,3)).*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
cov_tl2(1:lm-LevMin+1,3) = sum(sum((t_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1) - mean_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1,3)).*(t_tl2(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1) - mean_tl2(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1,3)).*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
cov_tl3(1:lm-LevMin+1,3) = sum(sum((t_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1) - mean_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1,3)).*(t_tl3(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1) - mean_tl3(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1,3)).*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;

corr_tl1(1:lm-LevMin+1,3) = cov_tl1(1:lm-LevMin+1,3)./sqrt(var_nlm(1:lm-LevMin+1,3).*var_tl1(1:lm-LevMin+1,3));
corr_tl2(1:lm-LevMin+1,3) = cov_tl2(1:lm-LevMin+1,3)./sqrt(var_nlm(1:lm-LevMin+1,3).*var_tl2(1:lm-LevMin+1,3));
corr_tl3(1:lm-LevMin+1,3) = cov_tl3(1:lm-LevMin+1,3)./sqrt(var_nlm(1:lm-LevMin+1,3).*var_tl3(1:lm-LevMin+1,3));

% Specific Humidity
var_nlm(1:lm-LevMin+1,4) = sum(sum((q_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1) - mean_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1,4)).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
var_tl1(1:lm-LevMin+1,4) = sum(sum((q_tl1(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1) - mean_tl1(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1,4)).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
var_tl2(1:lm-LevMin+1,4) = sum(sum((q_tl2(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1) - mean_tl2(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1,4)).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
var_tl3(1:lm-LevMin+1,4) = sum(sum((q_tl3(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1) - mean_tl3(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1,4)).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;

cov_tl1(1:lm-LevMin+1,4) = sum(sum((q_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1) - mean_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1,4)).*(q_tl1(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1) - mean_tl1(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1,4)).*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
cov_tl2(1:lm-LevMin+1,4) = sum(sum((q_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1) - mean_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1,4)).*(q_tl2(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1) - mean_tl2(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1,4)).*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
cov_tl3(1:lm-LevMin+1,4) = sum(sum((q_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1) - mean_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1,4)).*(q_tl3(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1) - mean_tl3(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1,4)).*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;

corr_tl1(1:lm-LevMin+1,4) = cov_tl1(1:lm-LevMin+1,4)./sqrt(var_nlm(1:lm-LevMin+1,4).*var_tl1(1:lm-LevMin+1,4));
corr_tl2(1:lm-LevMin+1,4) = cov_tl2(1:lm-LevMin+1,4)./sqrt(var_nlm(1:lm-LevMin+1,4).*var_tl2(1:lm-LevMin+1,4));
corr_tl3(1:lm-LevMin+1,4) = cov_tl3(1:lm-LevMin+1,4)./sqrt(var_nlm(1:lm-LevMin+1,4).*var_tl3(1:lm-LevMin+1,4));

% Pressure
var_nlm(1:lm-LevMin+1,5) = sum(sum((p_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1) - mean_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1,5)).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
var_tl1(1:lm-LevMin+1,5) = sum(sum((p_tl1(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1) - mean_tl1(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1,5)).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
var_tl2(1:lm-LevMin+1,5) = sum(sum((p_tl2(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1) - mean_tl2(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1,5)).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
var_tl3(1:lm-LevMin+1,5) = sum(sum((p_tl3(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1) - mean_tl3(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1,5)).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;

cov_tl1(1:lm-LevMin+1,5) = sum(sum((p_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1) - mean_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1,5)).*(p_tl1(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1) - mean_tl1(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1,5)).*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
cov_tl2(1:lm-LevMin+1,5) = sum(sum((p_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1) - mean_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1,5)).*(p_tl2(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1) - mean_tl2(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1,5)).*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
cov_tl3(1:lm-LevMin+1,5) = sum(sum((p_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1) - mean_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1,5)).*(p_tl3(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1) - mean_tl3(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1,5)).*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;

corr_tl1(1:lm-LevMin+1,5) = cov_tl1(1:lm-LevMin+1,5)./sqrt(var_nlm(1:lm-LevMin+1,5).*var_tl1(1:lm-LevMin+1,5));
corr_tl2(1:lm-LevMin+1,5) = cov_tl2(1:lm-LevMin+1,5)./sqrt(var_nlm(1:lm-LevMin+1,5).*var_tl2(1:lm-LevMin+1,5));
corr_tl3(1:lm-LevMin+1,5) = cov_tl3(1:lm-LevMin+1,5)./sqrt(var_nlm(1:lm-LevMin+1,5).*var_tl3(1:lm-LevMin+1,5));

% Cloud Liquid Ice
var_nlm(1:lm-LevMin+1,6) = sum(sum((qi_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1) - mean_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1,6)).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
var_tl1(1:lm-LevMin+1,6) = sum(sum((qi_tl1(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1) - mean_tl1(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1,6)).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
var_tl2(1:lm-LevMin+1,6) = sum(sum((qi_tl2(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1) - mean_tl2(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1,6)).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
var_tl3(1:lm-LevMin+1,6) = sum(sum((qi_tl3(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1) - mean_tl3(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1,6)).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;

cov_tl1(1:lm-LevMin+1,6) = sum(sum((qi_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1) - mean_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1,6)).*(qi_tl1(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1) - mean_tl1(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1,6)).*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
cov_tl2(1:lm-LevMin+1,6) = sum(sum((qi_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1) - mean_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1,6)).*(qi_tl2(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1) - mean_tl2(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1,6)).*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
cov_tl3(1:lm-LevMin+1,6) = sum(sum((qi_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1) - mean_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1,6)).*(qi_tl3(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1) - mean_tl3(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1,6)).*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;

corr_tl1(1:lm-LevMin+1,6) = cov_tl1(1:lm-LevMin+1,6)./sqrt(var_nlm(1:lm-LevMin+1,6).*var_tl1(1:lm-LevMin+1,6));
corr_tl2(1:lm-LevMin+1,6) = cov_tl2(1:lm-LevMin+1,6)./sqrt(var_nlm(1:lm-LevMin+1,6).*var_tl2(1:lm-LevMin+1,6));
corr_tl3(1:lm-LevMin+1,6) = cov_tl3(1:lm-LevMin+1,6)./sqrt(var_nlm(1:lm-LevMin+1,6).*var_tl3(1:lm-LevMin+1,6));

% Cloud Liquid Water
var_nlm(1:lm-LevMin+1,7) = sum(sum((ql_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1) - mean_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1,7)).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
var_tl1(1:lm-LevMin+1,7) = sum(sum((ql_tl1(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1) - mean_tl1(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1,7)).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
var_tl2(1:lm-LevMin+1,7) = sum(sum((ql_tl2(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1) - mean_tl2(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1,7)).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
var_tl3(1:lm-LevMin+1,7) = sum(sum((ql_tl3(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1) - mean_tl3(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1,7)).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;

cov_tl1(1:lm-LevMin+1,7) = sum(sum((ql_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1) - mean_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1,7)).*(ql_tl1(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1) - mean_tl1(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1,7)).*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
cov_tl2(1:lm-LevMin+1,7) = sum(sum((ql_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1) - mean_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1,7)).*(ql_tl2(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1) - mean_tl2(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1,7)).*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
cov_tl3(1:lm-LevMin+1,7) = sum(sum((ql_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1) - mean_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1,7)).*(ql_tl3(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1) - mean_tl3(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1,7)).*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;

corr_tl1(1:lm-LevMin+1,7) = cov_tl1(1:lm-LevMin+1,7)./sqrt(var_nlm(1:lm-LevMin+1,7).*var_tl1(1:lm-LevMin+1,7));
corr_tl2(1:lm-LevMin+1,7) = cov_tl2(1:lm-LevMin+1,7)./sqrt(var_nlm(1:lm-LevMin+1,7).*var_tl2(1:lm-LevMin+1,7));
corr_tl3(1:lm-LevMin+1,7) = cov_tl3(1:lm-LevMin+1,7)./sqrt(var_nlm(1:lm-LevMin+1,7).*var_tl3(1:lm-LevMin+1,7));

% Ozone
var_nlm(1:lm-LevMin+1,8) = sum(sum((o3_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1) - mean_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1,8)).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
var_tl1(1:lm-LevMin+1,8) = sum(sum((o3_tl1(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1) - mean_tl1(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1,8)).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
var_tl2(1:lm-LevMin+1,8) = sum(sum((o3_tl2(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1) - mean_tl2(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1,8)).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
var_tl3(1:lm-LevMin+1,8) = sum(sum((o3_tl3(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1) - mean_tl3(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1,8)).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;

cov_tl1(1:lm-LevMin+1,8) = sum(sum((o3_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1) - mean_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1,8)).*(o3_tl1(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1) - mean_tl1(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1,8)).*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
cov_tl2(1:lm-LevMin+1,8) = sum(sum((o3_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1) - mean_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1,8)).*(o3_tl2(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1) - mean_tl2(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1,8)).*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
cov_tl3(1:lm-LevMin+1,8) = sum(sum((o3_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1) - mean_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1,8)).*(o3_tl3(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1) - mean_tl3(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1,8)).*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;

corr_tl1(1:lm-LevMin+1,8) = cov_tl1(1:lm-LevMin+1,8)./sqrt(var_nlm(1:lm-LevMin+1,8).*var_tl1(1:lm-LevMin+1,8));
corr_tl2(1:lm-LevMin+1,8) = cov_tl2(1:lm-LevMin+1,8)./sqrt(var_nlm(1:lm-LevMin+1,8).*var_tl2(1:lm-LevMin+1,8));
corr_tl3(1:lm-LevMin+1,8) = cov_tl3(1:lm-LevMin+1,8)./sqrt(var_nlm(1:lm-LevMin+1,8).*var_tl3(1:lm-LevMin+1,8));

fprintf('  Done correlation computation \n\n')



figure
set(gcf,'position',[3 343 1276 576])

subplot(1,8,1)
plot(corr_tl1(1:lm-LevMin+1,1),p(LevMin:lm),'b','LineWidth',lin_wid)
hold on
plot(corr_tl2(1:lm-LevMin+1,1),p(LevMin:lm),'r','LineWidth',lin_wid)
plot(corr_tl3(1:lm-LevMin+1,1),p(LevMin:lm),'g','LineWidth',lin_wid)
set(gca,'YDir','reverse')
xlim([0 1])
ylim([p(LevMin) p(lm)])
title('u (ms^{-1})','FontSize',fontsize1,'FontName','TimesNewRoman')
ylabel('Height (hPa)','FontSize',fontsize1,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

subplot(1,8,2)
plot(corr_tl1(1:lm-LevMin+1,2),p(LevMin:lm),'b','LineWidth',lin_wid)
hold on
plot(corr_tl2(1:lm-LevMin+1,2),p(LevMin:lm),'r','LineWidth',lin_wid)
plot(corr_tl3(1:lm-LevMin+1,2),p(LevMin:lm),'g','LineWidth',lin_wid)
set(gca,'YDir','reverse')
xlim([0 1])
ylim([p(LevMin) p(lm)])
title('v (ms^{-1})','FontSize',fontsize1,'FontName','TimesNewRoman')
set(gca,'YTickLabel',[])
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

subplot(1,8,3)
plot(corr_tl1(1:lm-LevMin+1,3),p(LevMin:lm),'b','LineWidth',lin_wid)
hold on
plot(corr_tl2(1:lm-LevMin+1,3),p(LevMin:lm),'r','LineWidth',lin_wid)
plot(corr_tl3(1:lm-LevMin+1,3),p(LevMin:lm),'g','LineWidth',lin_wid)
set(gca,'YDir','reverse')
xlim([0 1])
ylim([p(LevMin) p(lm)])
title('T_v (K)','FontSize',fontsize1,'FontName','TimesNewRoman')
set(gca,'YTickLabel',[])
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

subplot(1,8,4)
plot(corr_tl1(1:lm-LevMin+1,4),p(LevMin:lm),'b','LineWidth',lin_wid)
hold on
plot(corr_tl2(1:lm-LevMin+1,4),p(LevMin:lm),'r','LineWidth',lin_wid)
plot(corr_tl3(1:lm-LevMin+1,4),p(LevMin:lm),'g','LineWidth',lin_wid)
set(gca,'YDir','reverse')
xlim([0 1])
ylim([p(LevMin) p(lm)])
title('q (kgkg^{-1})','FontSize',fontsize1,'FontName','TimesNewRoman')
set(gca,'YTickLabel',[])
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

subplot(1,8,5)
plot(corr_tl1(1:lm-LevMin+1,5),p(LevMin:lm),'b','LineWidth',lin_wid)
hold on
plot(corr_tl2(1:lm-LevMin+1,5),p(LevMin:lm),'r','LineWidth',lin_wid)
plot(corr_tl3(1:lm-LevMin+1,5),p(LevMin:lm),'g','LineWidth',lin_wid)
set(gca,'YDir','reverse')
xlim([0 1])
ylim([p(LevMin) p(lm)])
title('p (Pa)','FontSize',fontsize1,'FontName','TimesNewRoman')
set(gca,'YTickLabel',[])
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

subplot(1,8,6)
plot(corr_tl1(1:lm-LevMin+1,6),p(LevMin:lm),'b','LineWidth',lin_wid)
hold on
plot(corr_tl2(1:lm-LevMin+1,6),p(LevMin:lm),'r','LineWidth',lin_wid)
plot(corr_tl3(1:lm-LevMin+1,6),p(LevMin:lm),'g','LineWidth',lin_wid)
set(gca,'YDir','reverse')
xlim([0 1])
ylim([p(LevMin) p(lm)])
title('q_i (kgkg^{-1})','FontSize',fontsize1,'FontName','TimesNewRoman')
set(gca,'YTickLabel',[])
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
% xlim([0 0.5])

subplot(1,8,7)
plot(corr_tl1(1:lm-LevMin+1,7),p(LevMin:lm),'b','LineWidth',lin_wid)
hold on
plot(corr_tl2(1:lm-LevMin+1,7),p(LevMin:lm),'r','LineWidth',lin_wid)
plot(corr_tl3(1:lm-LevMin+1,7),p(LevMin:lm),'g','LineWidth',lin_wid)
set(gca,'YDir','reverse')
xlim([0 1])
ylim([p(LevMin) p(lm)])
title('q_l (kgkg^{-1})','FontSize',fontsize1,'FontName','TimesNewRoman')
set(gca,'YTickLabel',[])
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
% xlim([0 0.5])

subplot(1,8,8)
plot(corr_tl1(1:lm-LevMin+1,8),p(LevMin:lm),'b','LineWidth',lin_wid)
hold on
plot(corr_tl2(1:lm-LevMin+1,8),p(LevMin:lm),'r','LineWidth',lin_wid)
plot(corr_tl3(1:lm-LevMin+1,8),p(LevMin:lm),'g','LineWidth',lin_wid)
set(gca,'YDir','reverse')
xlim([0 1])
ylim([p(0+1) p(lm)])
title('o3 (ppm)','FontSize',fontsize1,'FontName','TimesNewRoman')
ylabel('Height (hPa)','FontSize',fontsize1,'FontName','TimesNewRoman')

set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'YAxisLocation','right')
% xlim([0 0.5])





figure
set(gcf,'position',[3 343 1276 576])

subplot(1,8,1)
plot(corr_tl1(1:lm-LevMin+1,1),LevMin:lm,'b','LineWidth',lin_wid)
hold on
plot(corr_tl2(1:lm-LevMin+1,1),LevMin:lm,'r','LineWidth',lin_wid)
plot(corr_tl3(1:lm-LevMin+1,1),LevMin:lm,'g','LineWidth',lin_wid)
set(gca,'YDir','reverse')
xlim([0 1])
ylim([LevMin lm])
title('u (ms^{-1})','FontSize',fontsize1,'FontName','TimesNewRoman')
ylabel('Model level','FontSize',fontsize1,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

subplot(1,8,2)
plot(corr_tl1(1:lm-LevMin+1,2),LevMin:lm,'b','LineWidth',lin_wid)
hold on
plot(corr_tl2(1:lm-LevMin+1,2),LevMin:lm,'r','LineWidth',lin_wid)
plot(corr_tl3(1:lm-LevMin+1,2),LevMin:lm,'g','LineWidth',lin_wid)
set(gca,'YDir','reverse')
xlim([0 1])
ylim([LevMin lm])
title('v (ms^{-1})','FontSize',fontsize1,'FontName','TimesNewRoman')
set(gca,'YTickLabel',[])
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

subplot(1,8,3)
plot(corr_tl1(1:lm-LevMin+1,3),LevMin:lm,'b','LineWidth',lin_wid)
hold on
plot(corr_tl2(1:lm-LevMin+1,3),LevMin:lm,'r','LineWidth',lin_wid)
plot(corr_tl3(1:lm-LevMin+1,3),LevMin:lm,'g','LineWidth',lin_wid)
set(gca,'YDir','reverse')
xlim([0 1])
ylim([LevMin lm])
title('T_v (K)','FontSize',fontsize1,'FontName','TimesNewRoman')
set(gca,'YTickLabel',[])
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

subplot(1,8,4)
plot(corr_tl1(1:lm-LevMin+1,4),LevMin:lm,'b','LineWidth',lin_wid)
hold on
plot(corr_tl2(1:lm-LevMin+1,4),LevMin:lm,'r','LineWidth',lin_wid)
plot(corr_tl3(1:lm-LevMin+1,4),LevMin:lm,'g','LineWidth',lin_wid)
set(gca,'YDir','reverse')
xlim([0 1])
ylim([LevMin lm])
title('q (kgkg^{-1})','FontSize',fontsize1,'FontName','TimesNewRoman')
set(gca,'YTickLabel',[])
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

subplot(1,8,5)
plot(corr_tl1(1:lm-LevMin+1,5),LevMin:lm,'b','LineWidth',lin_wid)
hold on
plot(corr_tl2(1:lm-LevMin+1,5),LevMin:lm,'r','LineWidth',lin_wid)
plot(corr_tl3(1:lm-LevMin+1,5),LevMin:lm,'g','LineWidth',lin_wid)
set(gca,'YDir','reverse')
xlim([0 1])
ylim([LevMin lm])
title('p (Pa)','FontSize',fontsize1,'FontName','TimesNewRoman')
set(gca,'YTickLabel',[])
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

subplot(1,8,6)
plot(corr_tl1(1:lm-LevMin+1,6),LevMin:lm,'b','LineWidth',lin_wid)
hold on
plot(corr_tl2(1:lm-LevMin+1,6),LevMin:lm,'r','LineWidth',lin_wid)
plot(corr_tl3(1:lm-LevMin+1,6),LevMin:lm,'g','LineWidth',lin_wid)
set(gca,'YDir','reverse')
xlim([0 1])
ylim([LevMin lm])
title('q_i (kgkg^{-1})','FontSize',fontsize1,'FontName','TimesNewRoman')
set(gca,'YTickLabel',[])
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
% xlim([0 0.5])

subplot(1,8,7)
plot(corr_tl1(1:lm-LevMin+1,7),LevMin:lm,'b','LineWidth',lin_wid)
hold on
plot(corr_tl2(1:lm-LevMin+1,7),LevMin:lm,'r','LineWidth',lin_wid)
plot(corr_tl3(1:lm-LevMin+1,7),LevMin:lm,'g','LineWidth',lin_wid)
set(gca,'YDir','reverse')
xlim([0 1])
ylim([LevMin lm])
title('q_l (kgkg^{-1})','FontSize',fontsize1,'FontName','TimesNewRoman')
set(gca,'YTickLabel',[])
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
% xlim([0 0.5])

subplot(1,8,8)
plot(corr_tl1(1:lm-LevMin+1,8),LevMin:lm,'b','LineWidth',lin_wid)
hold on
plot(corr_tl2(1:lm-LevMin+1,8),LevMin:lm,'r','LineWidth',lin_wid)
plot(corr_tl3(1:lm-LevMin+1,8),LevMin:lm,'g','LineWidth',lin_wid)
set(gca,'YDir','reverse')
xlim([0 1])
ylim([LevMin lm])
title('o3 (ppm)','FontSize',fontsize1,'FontName','TimesNewRoman')
ylabel('Model level','FontSize',fontsize1,'FontName','TimesNewRoman')

set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'YAxisLocation','right')
% xlim([0 0.5])
