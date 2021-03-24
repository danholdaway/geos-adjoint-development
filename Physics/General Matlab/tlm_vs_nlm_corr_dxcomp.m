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

%Legend
leg1 = '\alpha = 1';
leg2 = '\alpha = 1/2';
leg3 = '\alpha = 1/10';

%Choose shapiro filter options for TLM
shp_tlm1 = '1';
shp_tlm2 = '1';
shp_tlm3 = '1';

%Choose 'other' options for TLM
opt = '_dry';
opt_tlm1 = '_traj_dry_a10';
opt_tlm2 = '_traj_dry_a05';
opt_tlm3 = '_traj_dry_a01';

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
dir = ['/discover/nobackup/drholdaw/',dir_free,'/prog/prog_free',opt,'/']; cd(dir)
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
dir = ['/discover/nobackup/drholdaw/',dir_repl,'/prog/prog_replay',opt,'/']; cd(dir)
file = [exp_repl,'.prog.eta.',datetimef(1:8),'_',datetimef(9:10),'z.nc4'];

u_replay1 = ncread(file,'u');
v_replay1 = ncread(file,'v');
t_replay1 = ncread(file,'tv');
q_replay1 = ncread(file,'sphu');
p_replay1 = ncread(file,'delp');
qi_replay1 = ncread(file,'qitot');
ql_replay1 = ncread(file,'qltot');
o3_replay1 = ncread(file,'ozone');

%Load Perturbed (analysis) st
dir = ['/discover/nobackup/drholdaw/',dir_repl,'/prog/prog_replay05',opt,'/']; cd(dir)
file = [exp_repl,'.prog.eta.',datetimef(1:8),'_',datetimef(9:10),'z.nc4'];

u_replay2 = ncread(file,'u');
v_replay2 = ncread(file,'v');
t_replay2 = ncread(file,'tv');
q_replay2 = ncread(file,'sphu');
p_replay2 = ncread(file,'delp');
qi_replay2 = ncread(file,'qitot');
ql_replay2 = ncread(file,'qltot');
o3_replay2 = ncread(file,'ozone');

%Load Perturbed (analysis) st
dir = ['/discover/nobackup/drholdaw/',dir_repl,'/prog/prog_replay01',opt,'/']; cd(dir)
file = [exp_repl,'.prog.eta.',datetimef(1:8),'_',datetimef(9:10),'z.nc4'];

u_replay3 = ncread(file,'u');
v_replay3 = ncread(file,'v');
t_replay3 = ncread(file,'tv');
q_replay3 = ncread(file,'sphu');
p_replay3 = ncread(file,'delp');
qi_replay3 = ncread(file,'qitot');
ql_replay3 = ncread(file,'qltot');
o3_replay3 = ncread(file,'ozone');

%Load first Tlm state to compare.
dir = ['/discover/nobackup/drholdaw/',dir_tlm1,'/sens.',datetimev(1:8),'.000000']; cd(dir)
file = [exp_tlm1,'.fvpert.eta.',datetimef(1:8),'_',datetimef(9:12),'z_ST03z_P',phy_tlm1,'_S',shp_tlm1,opt_tlm1,'.nc4'];

fprintf(' TLM1 file: \n')
disp(file)

u_tlm1 = ncread(file,'U');
v_tlm1 = ncread(file,'V');
t_tlm1 = ncread(file,'TV');
q_tlm1 = ncread(file,'QV');
p_tlm1 = ncread(file,'DP');
qi_tlm1 = ncread(file,'QI');
ql_tlm1 = ncread(file,'QL');
o3_tlm1 = ncread(file,'O3');

%Load second Tlm state to compare.
dir = ['/discover/nobackup/drholdaw/',dir_tlm2,'/sens.',datetimev(1:8),'.000000']; cd(dir)
file = [exp_tlm2,'.fvpert.eta.',datetimef(1:8),'_',datetimef(9:12),'z_ST03z_P',phy_tlm2,'_S',shp_tlm2,opt_tlm2,'.nc4'];

fprintf(' TLM2 file: \n')
disp(file)

u_tlm2 = ncread(file,'U');
v_tlm2 = ncread(file,'V');
t_tlm2 = ncread(file,'TV');
q_tlm2 = ncread(file,'QV');
p_tlm2 = ncread(file,'DP');
qi_tlm2 = ncread(file,'QI');
ql_tlm2 = ncread(file,'QL');
o3_tlm2 = ncread(file,'O3');

%Load third Tlm state to compare.
dir = ['/discover/nobackup/drholdaw/',dir_tlm3,'/sens.',datetimev(1:8),'.000000']; cd(dir)
file = [exp_tlm3,'.fvpert.eta.',datetimef(1:8),'_',datetimef(9:12),'z_ST03z_P',phy_tlm3,'_S',shp_tlm3,opt_tlm3,'.nc4'];

fprintf(' TLM3 file: \n')
disp(file)

u_tlm3 = ncread(file,'U');
v_tlm3 = ncread(file,'V');
t_tlm3 = ncread(file,'TV');
q_tlm3 = ncread(file,'QV');
p_tlm3 = ncread(file,'DP');
qi_tlm3 = ncread(file,'QI');
ql_tlm3 = ncread(file,'QL');
o3_tlm3 = ncread(file,'O3');

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
u_nlm1 = u_replay1 - u_free;
v_nlm1 = v_replay1 - v_free;
t_nlm1 = t_replay1 - t_free;
q_nlm1 = q_replay1 - q_free;
p_nlm1 = p_replay1 - p_free;
qi_nlm1 = qi_replay1 - qi_free;
ql_nlm1 = ql_replay1 - ql_free;
o3_nlm1 = o3_replay1 - o3_free;

u_nlm2 = u_replay2 - u_free;
v_nlm2 = v_replay2 - v_free;
t_nlm2 = t_replay2 - t_free;
q_nlm2 = q_replay2 - q_free;
p_nlm2 = p_replay2 - p_free;
qi_nlm2 = qi_replay2 - qi_free;
ql_nlm2 = ql_replay2 - ql_free;
o3_nlm2 = o3_replay2 - o3_free;

u_nlm3 = u_replay3 - u_free;
v_nlm3 = v_replay3 - v_free;
t_nlm3 = t_replay3 - t_free;
q_nlm3 = q_replay3 - q_free;
p_nlm3 = p_replay3 - p_free;
qi_nlm3 = qi_replay3 - qi_free;
ql_nlm3 = ql_replay3 - ql_free;
o3_nlm3 = o3_replay3 - o3_free;

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
var_nlm1 = zeros(lm-LevMin+1,8);
var_nlm2 = zeros(lm-LevMin+1,8);
var_nlm3 = zeros(lm-LevMin+1,8);

var_tlm1 = zeros(lm-LevMin+1,8);
var_tlm2 = zeros(lm-LevMin+1,8);
var_tlm3 = zeros(lm-LevMin+1,8);

cov_tlm1 = zeros(lm-LevMin+1,8);
cov_tlm2 = zeros(lm-LevMin+1,8);
cov_tlm3 = zeros(lm-LevMin+1,8);

%Initialize correlations
corr_tlm1 = zeros(lm-LevMin+1,8);
corr_tlm2 = zeros(lm-LevMin+1,8);
corr_tlm3 = zeros(lm-LevMin+1,8);

% U winds
var_nlm1(1:lm-LevMin+1,1) = sum(sum(u_nlm1(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
var_nlm2(1:lm-LevMin+1,1) = sum(sum(u_nlm2(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
var_nlm3(1:lm-LevMin+1,1) = sum(sum(u_nlm3(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;

var_tlm1(1:lm-LevMin+1,1) = sum(sum(u_tlm1(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
var_tlm2(1:lm-LevMin+1,1) = sum(sum(u_tlm2(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
var_tlm3(1:lm-LevMin+1,1) = sum(sum(u_tlm3(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;

cov_tlm1(1:lm-LevMin+1,1) = sum(sum(u_nlm1(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).*u_tlm1(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
cov_tlm2(1:lm-LevMin+1,1) = sum(sum(u_nlm2(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).*u_tlm2(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
cov_tlm3(1:lm-LevMin+1,1) = sum(sum(u_nlm3(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).*u_tlm3(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;

corr_tlm1(1:lm-LevMin+1,1) = cov_tlm1(1:lm-LevMin+1,1)./sqrt(var_nlm1(1:lm-LevMin+1,1).*var_tlm1(1:lm-LevMin+1,1));
corr_tlm2(1:lm-LevMin+1,1) = cov_tlm2(1:lm-LevMin+1,1)./sqrt(var_nlm2(1:lm-LevMin+1,1).*var_tlm2(1:lm-LevMin+1,1));
corr_tlm3(1:lm-LevMin+1,1) = cov_tlm3(1:lm-LevMin+1,1)./sqrt(var_nlm3(1:lm-LevMin+1,1).*var_tlm3(1:lm-LevMin+1,1));

% V winds
var_nlm1(1:lm-LevMin+1,2) = sum(sum(v_nlm1(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
var_nlm2(1:lm-LevMin+1,2) = sum(sum(v_nlm2(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
var_nlm3(1:lm-LevMin+1,2) = sum(sum(v_nlm3(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;

var_tlm1(1:lm-LevMin+1,2) = sum(sum(v_tlm1(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
var_tlm2(1:lm-LevMin+1,2) = sum(sum(v_tlm2(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
var_tlm3(1:lm-LevMin+1,2) = sum(sum(v_tlm3(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;

cov_tlm1(1:lm-LevMin+1,2) = sum(sum(v_nlm1(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).*v_tlm1(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
cov_tlm2(1:lm-LevMin+1,2) = sum(sum(v_nlm2(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).*v_tlm2(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
cov_tlm3(1:lm-LevMin+1,2) = sum(sum(v_nlm3(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).*v_tlm3(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;

corr_tlm1(1:lm-LevMin+1,2) = cov_tlm1(1:lm-LevMin+1,2)./sqrt(var_nlm1(1:lm-LevMin+1,2).*var_tlm1(1:lm-LevMin+1,2));
corr_tlm2(1:lm-LevMin+1,2) = cov_tlm2(1:lm-LevMin+1,2)./sqrt(var_nlm2(1:lm-LevMin+1,2).*var_tlm2(1:lm-LevMin+1,2));
corr_tlm3(1:lm-LevMin+1,2) = cov_tlm3(1:lm-LevMin+1,2)./sqrt(var_nlm3(1:lm-LevMin+1,2).*var_tlm3(1:lm-LevMin+1,2));

% Temperature
var_nlm1(1:lm-LevMin+1,3) = sum(sum(t_nlm1(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
var_nlm2(1:lm-LevMin+1,3) = sum(sum(t_nlm2(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
var_nlm3(1:lm-LevMin+1,3) = sum(sum(t_nlm3(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;

var_tlm1(1:lm-LevMin+1,3) = sum(sum(t_tlm1(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
var_tlm2(1:lm-LevMin+1,3) = sum(sum(t_tlm2(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
var_tlm3(1:lm-LevMin+1,3) = sum(sum(t_tlm3(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;

cov_tlm1(1:lm-LevMin+1,3) = sum(sum(t_nlm1(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).*t_tlm1(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
cov_tlm2(1:lm-LevMin+1,3) = sum(sum(t_nlm2(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).*t_tlm2(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
cov_tlm3(1:lm-LevMin+1,3) = sum(sum(t_nlm3(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).*t_tlm3(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;

corr_tlm1(1:lm-LevMin+1,3) = cov_tlm1(1:lm-LevMin+1,3)./sqrt(var_nlm1(1:lm-LevMin+1,3).*var_tlm1(1:lm-LevMin+1,3));
corr_tlm2(1:lm-LevMin+1,3) = cov_tlm2(1:lm-LevMin+1,3)./sqrt(var_nlm2(1:lm-LevMin+1,3).*var_tlm2(1:lm-LevMin+1,3));
corr_tlm3(1:lm-LevMin+1,3) = cov_tlm3(1:lm-LevMin+1,3)./sqrt(var_nlm3(1:lm-LevMin+1,3).*var_tlm3(1:lm-LevMin+1,3));

% Specific Humidity
var_nlm1(1:lm-LevMin+1,4) = sum(sum(q_nlm1(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
var_nlm2(1:lm-LevMin+1,4) = sum(sum(q_nlm2(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
var_nlm3(1:lm-LevMin+1,4) = sum(sum(q_nlm3(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;

var_tlm1(1:lm-LevMin+1,4) = sum(sum(q_tlm1(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
var_tlm2(1:lm-LevMin+1,4) = sum(sum(q_tlm2(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
var_tlm3(1:lm-LevMin+1,4) = sum(sum(q_tlm3(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;

cov_tlm1(1:lm-LevMin+1,4) = sum(sum(q_nlm1(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).*q_tlm1(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
cov_tlm2(1:lm-LevMin+1,4) = sum(sum(q_nlm2(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).*q_tlm2(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
cov_tlm3(1:lm-LevMin+1,4) = sum(sum(q_nlm3(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).*q_tlm3(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;

corr_tlm1(1:lm-LevMin+1,4) = cov_tlm1(1:lm-LevMin+1,4)./sqrt(var_nlm1(1:lm-LevMin+1,4).*var_tlm1(1:lm-LevMin+1,4));
corr_tlm2(1:lm-LevMin+1,4) = cov_tlm2(1:lm-LevMin+1,4)./sqrt(var_nlm2(1:lm-LevMin+1,4).*var_tlm2(1:lm-LevMin+1,4));
corr_tlm3(1:lm-LevMin+1,4) = cov_tlm3(1:lm-LevMin+1,4)./sqrt(var_nlm3(1:lm-LevMin+1,4).*var_tlm3(1:lm-LevMin+1,4));

% Pressure
var_nlm1(1:lm-LevMin+1,5) = sum(sum(p_nlm1(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
var_nlm2(1:lm-LevMin+1,5) = sum(sum(p_nlm2(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
var_nlm3(1:lm-LevMin+1,5) = sum(sum(p_nlm3(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;

var_tlm1(1:lm-LevMin+1,5) = sum(sum(p_tlm1(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
var_tlm2(1:lm-LevMin+1,5) = sum(sum(p_tlm2(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
var_tlm3(1:lm-LevMin+1,5) = sum(sum(p_tlm3(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;

cov_tlm1(1:lm-LevMin+1,5) = sum(sum(p_nlm1(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).*p_tlm1(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
cov_tlm2(1:lm-LevMin+1,5) = sum(sum(p_nlm2(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).*p_tlm2(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
cov_tlm3(1:lm-LevMin+1,5) = sum(sum(p_nlm3(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).*p_tlm3(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;

corr_tlm1(1:lm-LevMin+1,5) = cov_tlm1(1:lm-LevMin+1,5)./sqrt(var_nlm1(1:lm-LevMin+1,5).*var_tlm1(1:lm-LevMin+1,5));
corr_tlm2(1:lm-LevMin+1,5) = cov_tlm2(1:lm-LevMin+1,5)./sqrt(var_nlm2(1:lm-LevMin+1,5).*var_tlm2(1:lm-LevMin+1,5));
corr_tlm3(1:lm-LevMin+1,5) = cov_tlm3(1:lm-LevMin+1,5)./sqrt(var_nlm3(1:lm-LevMin+1,5).*var_tlm3(1:lm-LevMin+1,5));

% Cloud Liquid Water
var_nlm1(1:lm-LevMin+1,6) = sum(sum(qi_nlm1(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
var_nlm2(1:lm-LevMin+1,6) = sum(sum(qi_nlm2(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
var_nlm3(1:lm-LevMin+1,6) = sum(sum(qi_nlm3(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;

var_tlm1(1:lm-LevMin+1,6) = sum(sum(qi_tlm1(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
var_tlm2(1:lm-LevMin+1,6) = sum(sum(qi_tlm2(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
var_tlm3(1:lm-LevMin+1,6) = sum(sum(qi_tlm3(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;

cov_tlm1(1:lm-LevMin+1,6) = sum(sum(qi_nlm1(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).*qi_tlm1(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
cov_tlm2(1:lm-LevMin+1,6) = sum(sum(qi_nlm2(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).*qi_tlm2(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
cov_tlm3(1:lm-LevMin+1,6) = sum(sum(qi_nlm3(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).*qi_tlm3(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;

corr_tlm1(1:lm-LevMin+1,6) = cov_tlm1(1:lm-LevMin+1,6)./sqrt(var_nlm1(1:lm-LevMin+1,6).*var_tlm1(1:lm-LevMin+1,6));
corr_tlm2(1:lm-LevMin+1,6) = cov_tlm2(1:lm-LevMin+1,6)./sqrt(var_nlm2(1:lm-LevMin+1,6).*var_tlm2(1:lm-LevMin+1,6));
corr_tlm3(1:lm-LevMin+1,6) = cov_tlm3(1:lm-LevMin+1,6)./sqrt(var_nlm3(1:lm-LevMin+1,6).*var_tlm3(1:lm-LevMin+1,6));

% Cloud Liquid Ice
var_nlm1(1:lm-LevMin+1,7) = sum(sum(ql_nlm1(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
var_nlm2(1:lm-LevMin+1,7) = sum(sum(ql_nlm2(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
var_nlm3(1:lm-LevMin+1,7) = sum(sum(ql_nlm3(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;

var_tlm1(1:lm-LevMin+1,7) = sum(sum(ql_tlm1(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
var_tlm2(1:lm-LevMin+1,7) = sum(sum(ql_tlm2(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
var_tlm3(1:lm-LevMin+1,7) = sum(sum(ql_tlm3(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;

cov_tlm1(1:lm-LevMin+1,7) = sum(sum(ql_nlm1(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).*ql_tlm1(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
cov_tlm2(1:lm-LevMin+1,7) = sum(sum(ql_nlm2(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).*ql_tlm2(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
cov_tlm3(1:lm-LevMin+1,7) = sum(sum(ql_nlm3(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).*ql_tlm3(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;

corr_tlm1(1:lm-LevMin+1,7) = cov_tlm1(1:lm-LevMin+1,7)./sqrt(var_nlm1(1:lm-LevMin+1,7).*var_tlm1(1:lm-LevMin+1,7));
corr_tlm2(1:lm-LevMin+1,7) = cov_tlm2(1:lm-LevMin+1,7)./sqrt(var_nlm2(1:lm-LevMin+1,7).*var_tlm2(1:lm-LevMin+1,7));
corr_tlm3(1:lm-LevMin+1,7) = cov_tlm3(1:lm-LevMin+1,7)./sqrt(var_nlm3(1:lm-LevMin+1,7).*var_tlm3(1:lm-LevMin+1,7));

% Ozone
var_nlm1(1:lm-LevMin+1,8) = sum(sum(o3_nlm1(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
var_nlm2(1:lm-LevMin+1,8) = sum(sum(o3_nlm2(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
var_nlm3(1:lm-LevMin+1,8) = sum(sum(o3_nlm3(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;

var_tlm1(1:lm-LevMin+1,8) = sum(sum(o3_tlm1(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
var_tlm2(1:lm-LevMin+1,8) = sum(sum(o3_tlm2(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
var_tlm3(1:lm-LevMin+1,8) = sum(sum(o3_tlm3(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;

cov_tlm1(1:lm-LevMin+1,8) = sum(sum(o3_nlm1(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).*o3_tlm1(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
cov_tlm2(1:lm-LevMin+1,8) = sum(sum(o3_nlm2(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).*o3_tlm2(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
cov_tlm3(1:lm-LevMin+1,8) = sum(sum(o3_nlm3(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).*o3_tlm3(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;

corr_tlm1(1:lm-LevMin+1,8) = cov_tlm1(1:lm-LevMin+1,8)./sqrt(var_nlm1(1:lm-LevMin+1,8).*var_tlm1(1:lm-LevMin+1,8));
corr_tlm2(1:lm-LevMin+1,8) = cov_tlm2(1:lm-LevMin+1,8)./sqrt(var_nlm2(1:lm-LevMin+1,8).*var_tlm2(1:lm-LevMin+1,8));
corr_tlm3(1:lm-LevMin+1,8) = cov_tlm3(1:lm-LevMin+1,8)./sqrt(var_nlm3(1:lm-LevMin+1,8).*var_tlm3(1:lm-LevMin+1,8));

fprintf('  Done correlation computation \n\n')



figure
set(gcf,'position',[3 343 1276 576])

subplot(1,8,1)
plot(corr_tlm1(1:lm-LevMin+1,1),p(LevMin:lm),'b','LineWidth',lin_wid)
hold on
plot(corr_tlm2(1:lm-LevMin+1,1),p(LevMin:lm),'r','LineWidth',lin_wid)
plot(corr_tlm3(1:lm-LevMin+1,1),p(LevMin:lm),'g','LineWidth',lin_wid)
set(gca,'YDir','reverse')
xlim([0 1])
ylim([p(LevMin) p(lm)])
title('u (ms^{-1})','FontSize',fontsize1,'FontName','TimesNewRoman')
ylabel('Height (hPa)','FontSize',fontsize1,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

subplot(1,8,2)
plot(corr_tlm1(1:lm-LevMin+1,2),p(LevMin:lm),'b','LineWidth',lin_wid)
hold on
plot(corr_tlm2(1:lm-LevMin+1,2),p(LevMin:lm),'r','LineWidth',lin_wid)
plot(corr_tlm3(1:lm-LevMin+1,2),p(LevMin:lm),'g','LineWidth',lin_wid)
set(gca,'YDir','reverse')
xlim([0 1])
ylim([p(LevMin) p(lm)])
title('v (ms^{-1})','FontSize',fontsize1,'FontName','TimesNewRoman')
set(gca,'YTickLabel',[])
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

subplot(1,8,3)
plot(corr_tlm1(1:lm-LevMin+1,3),p(LevMin:lm),'b','LineWidth',lin_wid)
hold on
plot(corr_tlm2(1:lm-LevMin+1,3),p(LevMin:lm),'r','LineWidth',lin_wid)
plot(corr_tlm3(1:lm-LevMin+1,3),p(LevMin:lm),'g','LineWidth',lin_wid)
set(gca,'YDir','reverse')
xlim([0 1])
ylim([p(LevMin) p(lm)])
title('T_v (K)','FontSize',fontsize1,'FontName','TimesNewRoman')
set(gca,'YTickLabel',[])
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

subplot(1,8,4)
plot(corr_tlm1(1:lm-LevMin+1,4),p(LevMin:lm),'b','LineWidth',lin_wid)
hold on
plot(corr_tlm2(1:lm-LevMin+1,4),p(LevMin:lm),'r','LineWidth',lin_wid)
plot(corr_tlm3(1:lm-LevMin+1,4),p(LevMin:lm),'g','LineWidth',lin_wid)
set(gca,'YDir','reverse')
xlim([0 1])
ylim([p(LevMin) p(lm)])
title('q (kgkg^{-1})','FontSize',fontsize1,'FontName','TimesNewRoman')
set(gca,'YTickLabel',[])
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

subplot(1,8,5)
plot(corr_tlm1(1:lm-LevMin+1,5),p(LevMin:lm),'b','LineWidth',lin_wid)
hold on
plot(corr_tlm2(1:lm-LevMin+1,5),p(LevMin:lm),'r','LineWidth',lin_wid)
plot(corr_tlm3(1:lm-LevMin+1,5),p(LevMin:lm),'g','LineWidth',lin_wid)
set(gca,'YDir','reverse')
xlim([0 1])
ylim([p(LevMin) p(lm)])
title('p (Pa)','FontSize',fontsize1,'FontName','TimesNewRoman')
set(gca,'YTickLabel',[])
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

subplot(1,8,6)
plot(corr_tlm1(1:lm-LevMin+1,6),p(LevMin:lm),'b','LineWidth',lin_wid)
hold on
plot(corr_tlm2(1:lm-LevMin+1,6),p(LevMin:lm),'r','LineWidth',lin_wid)
plot(corr_tlm3(1:lm-LevMin+1,6),p(LevMin:lm),'g','LineWidth',lin_wid)
set(gca,'YDir','reverse')
xlim([0 1])
ylim([p(LevMin) p(lm)])
title('q_i (kgkg^{-1})','FontSize',fontsize1,'FontName','TimesNewRoman')
set(gca,'YTickLabel',[])
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

subplot(1,8,7)
plot(corr_tlm1(1:lm-LevMin+1,7),p(LevMin:lm),'b','LineWidth',lin_wid)
hold on
plot(corr_tlm2(1:lm-LevMin+1,7),p(LevMin:lm),'r','LineWidth',lin_wid)
plot(corr_tlm3(1:lm-LevMin+1,7),p(LevMin:lm),'g','LineWidth',lin_wid)
set(gca,'YDir','reverse')
xlim([0 1])
ylim([p(LevMin) p(lm)])
title('q_l (kgkg^{-1})','FontSize',fontsize1,'FontName','TimesNewRoman')
set(gca,'YTickLabel',[])
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

subplot(1,8,8)
plot(corr_tlm1(1:lm-LevMin+1,8),p(LevMin:lm),'b','LineWidth',lin_wid)
hold on
plot(corr_tlm2(1:lm-LevMin+1,8),p(LevMin:lm),'r','LineWidth',lin_wid)
plot(corr_tlm3(1:lm-LevMin+1,8),p(LevMin:lm),'g','LineWidth',lin_wid)
set(gca,'YDir','reverse')
xlim([0 1])
ylim([p(0+1) p(lm)])
title('o3 (ppm)','FontSize',fontsize1,'FontName','TimesNewRoman')
ylabel('Height (hPa)','FontSize',fontsize1,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'YAxisLocation','right')
legend(leg1,leg2,leg3)


figure
set(gcf,'position',[3 343 1276 576])

subplot(1,8,1)
plot(corr_tlm1(1:lm-LevMin+1,1),LevMin:lm,'b','LineWidth',lin_wid)
hold on
plot(corr_tlm2(1:lm-LevMin+1,1),LevMin:lm,'r','LineWidth',lin_wid)
plot(corr_tlm3(1:lm-LevMin+1,1),LevMin:lm,'g','LineWidth',lin_wid)
set(gca,'YDir','reverse')
xlim([0 1])
ylim([LevMin lm])
title('u (ms^{-1})','FontSize',fontsize1,'FontName','TimesNewRoman')
ylabel('Model level','FontSize',fontsize1,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

subplot(1,8,2)
plot(corr_tlm1(1:lm-LevMin+1,2),LevMin:lm,'b','LineWidth',lin_wid)
hold on
plot(corr_tlm2(1:lm-LevMin+1,2),LevMin:lm,'r','LineWidth',lin_wid)
plot(corr_tlm3(1:lm-LevMin+1,2),LevMin:lm,'g','LineWidth',lin_wid)
set(gca,'YDir','reverse')
xlim([0 1])
ylim([LevMin lm])
title('v (ms^{-1})','FontSize',fontsize1,'FontName','TimesNewRoman')
set(gca,'YTickLabel',[])
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

subplot(1,8,3)
plot(corr_tlm1(1:lm-LevMin+1,3),LevMin:lm,'b','LineWidth',lin_wid)
hold on
plot(corr_tlm2(1:lm-LevMin+1,3),LevMin:lm,'r','LineWidth',lin_wid)
plot(corr_tlm3(1:lm-LevMin+1,3),LevMin:lm,'g','LineWidth',lin_wid)
set(gca,'YDir','reverse')
xlim([0 1])
ylim([LevMin lm])
title('T_v (K)','FontSize',fontsize1,'FontName','TimesNewRoman')
set(gca,'YTickLabel',[])
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

subplot(1,8,4)
plot(corr_tlm1(1:lm-LevMin+1,4),LevMin:lm,'b','LineWidth',lin_wid)
hold on
plot(corr_tlm2(1:lm-LevMin+1,4),LevMin:lm,'r','LineWidth',lin_wid)
plot(corr_tlm3(1:lm-LevMin+1,4),LevMin:lm,'g','LineWidth',lin_wid)
set(gca,'YDir','reverse')
xlim([0 1])
ylim([LevMin lm])
title('q (kgkg^{-1})','FontSize',fontsize1,'FontName','TimesNewRoman')
set(gca,'YTickLabel',[])
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

subplot(1,8,5)
plot(corr_tlm1(1:lm-LevMin+1,5),LevMin:lm,'b','LineWidth',lin_wid)
hold on
plot(corr_tlm2(1:lm-LevMin+1,5),LevMin:lm,'r','LineWidth',lin_wid)
plot(corr_tlm3(1:lm-LevMin+1,5),LevMin:lm,'g','LineWidth',lin_wid)
set(gca,'YDir','reverse')
xlim([0 1])
ylim([LevMin lm])
title('p (Pa)','FontSize',fontsize1,'FontName','TimesNewRoman')
set(gca,'YTickLabel',[])
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

subplot(1,8,6)
plot(corr_tlm1(1:lm-LevMin+1,6),LevMin:lm,'b','LineWidth',lin_wid)
hold on
plot(corr_tlm2(1:lm-LevMin+1,6),LevMin:lm,'r','LineWidth',lin_wid)
plot(corr_tlm3(1:lm-LevMin+1,6),LevMin:lm,'g','LineWidth',lin_wid)
set(gca,'YDir','reverse')
xlim([0 1])
ylim([LevMin lm])
title('q_i (kgkg^{-1})','FontSize',fontsize1,'FontName','TimesNewRoman')
set(gca,'YTickLabel',[])
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

subplot(1,8,7)
plot(corr_tlm1(1:lm-LevMin+1,7),LevMin:lm,'b','LineWidth',lin_wid)
hold on
plot(corr_tlm2(1:lm-LevMin+1,7),LevMin:lm,'r','LineWidth',lin_wid)
plot(corr_tlm3(1:lm-LevMin+1,7),LevMin:lm,'g','LineWidth',lin_wid)
set(gca,'YDir','reverse')
xlim([0 1])
ylim([LevMin lm])
title('q_l (kgkg^{-1})','FontSize',fontsize1,'FontName','TimesNewRoman')
set(gca,'YTickLabel',[])
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

subplot(1,8,8)
plot(corr_tlm1(1:lm-LevMin+1,8),LevMin:lm,'b','LineWidth',lin_wid)
hold on
plot(corr_tlm2(1:lm-LevMin+1,8),LevMin:lm,'r','LineWidth',lin_wid)
plot(corr_tlm3(1:lm-LevMin+1,8),LevMin:lm,'g','LineWidth',lin_wid)
set(gca,'YDir','reverse')
xlim([0 1])
ylim([LevMin lm])
title('o3 (ppm)','FontSize',fontsize1,'FontName','TimesNewRoman')
ylabel('Model level','FontSize',fontsize1,'FontName','TimesNewRoman')
legend(leg1,leg2,leg3)
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'YAxisLocation','right')