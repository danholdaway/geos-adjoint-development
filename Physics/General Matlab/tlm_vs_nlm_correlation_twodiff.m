close all
clear
clc
mydir = pwd;

%Choose Region,
% 1 = Global
% 2 = Tropics 23S to 23N, below 100hPa
% 3 = Northern hemisphere
% 4 = Southern hemisphere
% 5 = Dust over sahara (30W to 40E and 10S to 50N)
region = 1;

if region == 1 
    fprintf('Computing correlations globally \n')
elseif region == 2
    fprintf('Computing correlations for the tropics \n')  
elseif region == 3
    fprintf('Computing correlations for the northern hemisphere \n')
elseif region == 4
    fprintf('Computing correlations for the southern hemisphere \n')
elseif region == 5
    fprintf('Computing correlations for the Saharan dust region \n')
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

LevMin = 30;
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
elseif region == 5
    LonMin = 241;
    LonMax = 353;
    LatMin = 161;
    LatMax = 282;
end


cd /home/drholdaw/LinearisedPhysics/Inputs/
pref

%Load Free (background) State.
cd /discover/nobackup/drholdaw/btmp.25956/prog/prog_free
file = ['v000_C180.prog.eta.2014020',file_end,'z.nc4'];

lon = ncread(file,'lon');
lat = ncread(file,'lat');
lev = ncread(file,'lev');

u_free1 = ncread(file,'u',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
v_free1 = ncread(file,'v',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
t_free1 = ncread(file,'tv',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
q_free1 = ncread(file,'sphu',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
p_free1 = ncread(file,'delp',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
qi_free1 = ncread(file,'qitot',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
ql_free1 = ncread(file,'qltot',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
o3_free1 = ncread(file,'ozone',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
% d1_free1 = ncread(file,'DU001',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
% d2_free1 = ncread(file,'DU002',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
% d3_free1 = ncread(file,'DU003',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
% d4_free1 = ncread(file,'DU004',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
% d5_free1 = ncread(file,'DU005',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);

%Load Perturbed (analysis) st
cd /discover/nobackup/drholdaw/btmp.25956/prog/prog_replay/
file = ['v000_C180.prog.eta.2014020',file_end,'z.nc4'];

u_replay1 = ncread(file,'u',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
v_replay1 = ncread(file,'v',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
t_replay1 = ncread(file,'tv',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
q_replay1 = ncread(file,'sphu',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
p_replay1 = ncread(file,'delp',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
qi_replay1 = ncread(file,'qitot',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
ql_replay1 = ncread(file,'qltot',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
o3_replay1 = ncread(file,'ozone',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
% d1_replay1 = ncread(file,'DU001',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
% d2_replay1 = ncread(file,'DU002',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
% d3_replay1 = ncread(file,'DU003',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
% d4_replay1 = ncread(file,'DU004',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
% d5_replay1 = ncread(file,'DU005',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);

%Load Free (background) State.
cd /discover/nobackup/drholdaw/btmp.25956/prog/prog_free
file = ['v000_C180.prog.eta.2014020',file_end,'z.nc4'];

u_free2 = ncread(file,'u',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
v_free2 = ncread(file,'v',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
t_free2 = ncread(file,'tv',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
q_free2 = ncread(file,'sphu',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
p_free2 = ncread(file,'delp',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
qi_free2 = ncread(file,'qitot',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
ql_free2 = ncread(file,'qltot',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
o3_free2 = ncread(file,'ozone',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
% d1_free2 = ncread(file,'DU001',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
% d2_free2 = ncread(file,'DU002',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
% d3_free2 = ncread(file,'DU003',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
% d4_free2 = ncread(file,'DU004',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
% d5_free2 = ncread(file,'DU005',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);

%Load Perturbed (analysis) st
cd /discover/nobackup/drholdaw/btmp.25956/prog/prog_replay/
file = ['v000_C180.prog.eta.2014020',file_end,'z.nc4'];

u_replay2 = ncread(file,'u',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
v_replay2 = ncread(file,'v',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
t_replay2 = ncread(file,'tv',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
q_replay2 = ncread(file,'sphu',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
p_replay2 = ncread(file,'delp',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
qi_replay2 = ncread(file,'qitot',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
ql_replay2 = ncread(file,'qltot',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
o3_replay2 = ncread(file,'ozone',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
% d1_replay2 = ncread(file,'DU001',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
% d2_replay2 = ncread(file,'DU002',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
% d3_replay2 = ncread(file,'DU003',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
% d4_replay2 = ncread(file,'DU004',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
% d5_replay2 = ncread(file,'DU005',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);

%Load first TLM state to compare.
cd /discover/nobackup/drholdaw/btmp.25956/sens.20140202.000000/
file = ['v000_C180.fvpert.eta.2014020',file_end,'z_ST03z_DUST_P20010_S0.nc4'];
% file = 'fvpertX.eta.nc4';

u_tlm1 = ncread(file,'U',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
v_tlm1 = ncread(file,'V',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
t_tlm1 = ncread(file,'TV',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
q_tlm1 = ncread(file,'QV',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
p_tlm1 = ncread(file,'DP',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
qi_tlm1 = ncread(file,'QI',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
ql_tlm1 = ncread(file,'QL',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
o3_tlm1 = ncread(file,'O3',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
% d1_tlm1 = ncread(file,'DU001',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
% d2_tlm1 = ncread(file,'DU002',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
% d3_tlm1 = ncread(file,'DU003',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
% d4_tlm1 = ncread(file,'DU004',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
% d5_tlm1 = ncread(file,'DU005',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);

%Load second TLM state to compare.
cd /discover/nobackup/drholdaw/btmp.25956/sens.20140202.000000/
file = ['v000_C180.fvpert.eta.2014020',file_end,'z_ST03z_DUST_P22010_S0.nc4'];
% file = 'fvpertX.eta.nc4';

u_tlm2 = ncread(file,'U',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
v_tlm2 = ncread(file,'V',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
t_tlm2 = ncread(file,'TV',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
q_tlm2 = ncread(file,'QV',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
p_tlm2 = ncread(file,'DP',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
qi_tlm2 = ncread(file,'QI',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
ql_tlm2 = ncread(file,'QL',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
o3_tlm2 = ncread(file,'O3',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
% d1_tlm2 = ncread(file,'DU001',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
% d2_tlm2 = ncread(file,'DU002',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
% d3_tlm2 = ncread(file,'DU003',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
% d4_tlm2 = ncread(file,'DU004',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
% d5_tlm2 = ncread(file,'DU005',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);

cd(mydir)

%Compute NL perturbation trajectory.
u_nlm1 = u_replay1 - u_free1;
v_nlm1 = v_replay1 - v_free1;
t_nlm1 = t_replay1 - t_free1;
q_nlm1 = q_replay1 - q_free1;
p_nlm1 = p_replay1 - p_free1;
qi_nlm1 = qi_replay1 - qi_free1;
ql_nlm1 = ql_replay1 - ql_free1;
o3_nlm1 = o3_replay1 - o3_free1;
% d1_nlm1 = d1_replay1 - d1_free1;
% d2_nlm1 = d2_replay1 - d2_free1;
% d3_nlm1 = d3_replay1 - d3_free1;
% d4_nlm1 = d4_replay1 - d4_free1;
% d5_nlm1 = d5_replay1 - d5_free1;

u_nlm2 = u_replay2 - u_free2;
v_nlm2 = v_replay2 - v_free2;
t_nlm2 = t_replay2 - t_free2;
q_nlm2 = q_replay2 - q_free2;
p_nlm2 = p_replay2 - p_free2;
qi_nlm2 = qi_replay2 - qi_free2;
ql_nlm2 = ql_replay2 - ql_free2;
o3_nlm2 = o3_replay2 - o3_free2;
% d1_nlm2 = d1_replay2 - d1_free2;
% d2_nlm2 = d2_replay2 - d2_free2;
% d3_nlm2 = d3_replay2 - d3_free2;
% d4_nlm2 = d4_replay2 - d4_free2;
% d5_nlm2 = d5_replay2 - d5_free2;

fprintf('Beginning correlation calculation \n')

%Weighting 
r_earth = 6378.1;
lat_pole = lat+lat(1)*-1;
A = zeros(im,jm,LevMax-LevMin+1);
for k = 1:LevMax-LevMin+1
    A(1,:,k) = (2*pi^2*r_earth^2/(length(lon)*length(lat))) * sind(lat_pole);
    for j = 2:length(lon)
        A(j,:,k) = A(1,:,k); 
    end
end
dA = sum(sum(A(1:LevMax-LevMin+1,:,1)));

%Variance, covariance and correlation
%Second index is number of variables.
var_nlm1 = zeros(LevMax-LevMin+1,13);
var_nlm2 = zeros(LevMax-LevMin+1,13);
var_tlm1 = zeros(LevMax-LevMin+1,13);
var_tlm2 = zeros(LevMax-LevMin+1,13);

cov_tlm1_nlm1 = zeros(LevMax-LevMin+1,13);
cov_tlm2_nlm2 = zeros(LevMax-LevMin+1,13);

%Initialize correlations
corr_tlm1_nlm1 = zeros(LevMax-LevMin+1,13);
corr_tlm2_nlm2 = zeros(LevMax-LevMin+1,13);
corr_tlm3 = zeros(LevMax-LevMin+1,13);
corr_tlm4 = zeros(LevMax-LevMin+1,13);
corr_tlm5 = zeros(LevMax-LevMin+1,13);
corr_tlm6 = zeros(LevMax-LevMin+1,13);

% U winds
var_nlm1(1:LevMax-LevMin+1,1) = sum(sum(u_nlm1(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;
var_nlm2(1:LevMax-LevMin+1,1) = sum(sum(u_nlm2(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;

var_tlm1(1:LevMax-LevMin+1,1) = sum(sum(u_tlm1(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;
var_tlm2(1:LevMax-LevMin+1,1) = sum(sum(u_tlm2(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;

cov_tlm1_nlm1(1:LevMax-LevMin+1,1) = sum(sum(u_nlm1(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*u_tlm1(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;
cov_tlm2_nlm2(1:LevMax-LevMin+1,1) = sum(sum(u_nlm2(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*u_tlm2(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;

corr_tlm1_nlm1(1:LevMax-LevMin+1,1) = cov_tlm1_nlm1(1:LevMax-LevMin+1,1)./sqrt(var_nlm1(1:LevMax-LevMin+1,1).*var_tlm1(1:LevMax-LevMin+1,1));
corr_tlm2_nlm2(1:LevMax-LevMin+1,1) = cov_tlm2_nlm2(1:LevMax-LevMin+1,1)./sqrt(var_nlm2(1:LevMax-LevMin+1,1).*var_tlm2(1:LevMax-LevMin+1,1));

% V winds
var_nlm1(1:LevMax-LevMin+1,2) = sum(sum(v_nlm1(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;
var_nlm2(1:LevMax-LevMin+1,2) = sum(sum(v_nlm2(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;

var_tlm1(1:LevMax-LevMin+1,2) = sum(sum(v_tlm1(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;
var_tlm2(1:LevMax-LevMin+1,2) = sum(sum(v_tlm2(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;

cov_tlm1_nlm1(1:LevMax-LevMin+1,2) = sum(sum(v_nlm1(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*v_tlm1(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;
cov_tlm2_nlm2(1:LevMax-LevMin+1,2) = sum(sum(v_nlm2(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*v_tlm2(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;

corr_tlm1_nlm1(1:LevMax-LevMin+1,2) = cov_tlm1_nlm1(1:LevMax-LevMin+1,2)./sqrt(var_nlm1(1:LevMax-LevMin+1,2).*var_tlm1(1:LevMax-LevMin+1,2));
corr_tlm2_nlm2(1:LevMax-LevMin+1,2) = cov_tlm2_nlm2(1:LevMax-LevMin+1,2)./sqrt(var_nlm2(1:LevMax-LevMin+1,2).*var_tlm2(1:LevMax-LevMin+1,2));

% Temperature
var_nlm1(1:LevMax-LevMin+1,3) = sum(sum(t_nlm1(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;
var_nlm2(1:LevMax-LevMin+1,3) = sum(sum(t_nlm2(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;

var_tlm1(1:LevMax-LevMin+1,3) = sum(sum(t_tlm1(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;
var_tlm2(1:LevMax-LevMin+1,3) = sum(sum(t_tlm2(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;

cov_tlm1_nlm1(1:LevMax-LevMin+1,3) = sum(sum(t_nlm1(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*t_tlm1(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;
cov_tlm2_nlm2(1:LevMax-LevMin+1,3) = sum(sum(t_nlm2(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*t_tlm2(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;

corr_tlm1_nlm1(1:LevMax-LevMin+1,3) = cov_tlm1_nlm1(1:LevMax-LevMin+1,3)./sqrt(var_nlm1(1:LevMax-LevMin+1,3).*var_tlm1(1:LevMax-LevMin+1,3));
corr_tlm2_nlm2(1:LevMax-LevMin+1,3) = cov_tlm2_nlm2(1:LevMax-LevMin+1,3)./sqrt(var_nlm2(1:LevMax-LevMin+1,3).*var_tlm2(1:LevMax-LevMin+1,3));

% Specific Humidity
var_nlm1(1:LevMax-LevMin+1,4) = sum(sum(q_nlm1(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;
var_nlm2(1:LevMax-LevMin+1,4) = sum(sum(q_nlm2(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;

var_tlm1(1:LevMax-LevMin+1,4) = sum(sum(q_tlm1(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;
var_tlm2(1:LevMax-LevMin+1,4) = sum(sum(q_tlm2(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;

cov_tlm1_nlm1(1:LevMax-LevMin+1,4) = sum(sum(q_nlm1(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*q_tlm1(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;
cov_tlm2_nlm2(1:LevMax-LevMin+1,4) = sum(sum(q_nlm2(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*q_tlm2(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;

corr_tlm1_nlm1(1:LevMax-LevMin+1,4) = cov_tlm1_nlm1(1:LevMax-LevMin+1,4)./sqrt(var_nlm1(1:LevMax-LevMin+1,4).*var_tlm1(1:LevMax-LevMin+1,4));
corr_tlm2_nlm2(1:LevMax-LevMin+1,4) = cov_tlm2_nlm2(1:LevMax-LevMin+1,4)./sqrt(var_nlm2(1:LevMax-LevMin+1,4).*var_tlm2(1:LevMax-LevMin+1,4));

% Pressure
var_nlm1(1:LevMax-LevMin+1,5) = sum(sum(p_nlm1(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;
var_nlm2(1:LevMax-LevMin+1,5) = sum(sum(p_nlm2(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;

var_tlm1(1:LevMax-LevMin+1,5) = sum(sum(p_tlm1(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;
var_tlm2(1:LevMax-LevMin+1,5) = sum(sum(p_tlm2(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;

cov_tlm1_nlm1(1:LevMax-LevMin+1,5) = sum(sum(p_nlm1(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*p_tlm1(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;
cov_tlm2_nlm2(1:LevMax-LevMin+1,5) = sum(sum(p_nlm2(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*p_tlm2(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;

corr_tlm1_nlm1(1:LevMax-LevMin+1,5) = cov_tlm1_nlm1(1:LevMax-LevMin+1,5)./sqrt(var_nlm1(1:LevMax-LevMin+1,5).*var_tlm1(1:LevMax-LevMin+1,5));
corr_tlm2_nlm2(1:LevMax-LevMin+1,5) = cov_tlm2_nlm2(1:LevMax-LevMin+1,5)./sqrt(var_nlm2(1:LevMax-LevMin+1,5).*var_tlm2(1:LevMax-LevMin+1,5));

% Cloud Liquid Water
var_nlm1(1:LevMax-LevMin+1,6) = sum(sum(ql_nlm1(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;
var_nlm2(1:LevMax-LevMin+1,6) = sum(sum(ql_nlm2(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;

var_tlm1(1:LevMax-LevMin+1,6) = sum(sum(ql_tlm1(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;
var_tlm2(1:LevMax-LevMin+1,6) = sum(sum(ql_tlm2(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;

cov_tlm1_nlm1(1:LevMax-LevMin+1,6) = sum(sum(ql_nlm1(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*ql_tlm1(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;
cov_tlm2_nlm2(1:LevMax-LevMin+1,6) = sum(sum(ql_nlm2(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*ql_tlm2(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;

corr_tlm1_nlm1(1:LevMax-LevMin+1,6) = cov_tlm1_nlm1(1:LevMax-LevMin+1,6)./sqrt(var_nlm1(1:LevMax-LevMin+1,6).*var_tlm1(1:LevMax-LevMin+1,6));
corr_tlm2_nlm2(1:LevMax-LevMin+1,6) = cov_tlm2_nlm2(1:LevMax-LevMin+1,6)./sqrt(var_nlm2(1:LevMax-LevMin+1,6).*var_tlm2(1:LevMax-LevMin+1,6));

% Cloud Liquid Ice
var_nlm1(1:LevMax-LevMin+1,7) = sum(sum(qi_nlm1(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;
var_nlm2(1:LevMax-LevMin+1,7) = sum(sum(qi_nlm2(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;

var_tlm1(1:LevMax-LevMin+1,7) = sum(sum(qi_tlm1(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;
var_tlm2(1:LevMax-LevMin+1,7) = sum(sum(qi_tlm2(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;

cov_tlm1_nlm1(1:LevMax-LevMin+1,7) = sum(sum(qi_nlm1(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*qi_tlm1(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;
cov_tlm2_nlm2(1:LevMax-LevMin+1,7) = sum(sum(qi_nlm2(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*qi_tlm2(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;

corr_tlm1_nlm1(1:LevMax-LevMin+1,7) = cov_tlm1_nlm1(1:LevMax-LevMin+1,7)./sqrt(var_nlm1(1:LevMax-LevMin+1,7).*var_tlm1(1:LevMax-LevMin+1,7));
corr_tlm2_nlm2(1:LevMax-LevMin+1,7) = cov_tlm2_nlm2(1:LevMax-LevMin+1,7)./sqrt(var_nlm2(1:LevMax-LevMin+1,7).*var_tlm2(1:LevMax-LevMin+1,7));

% Ozone
var_nlm1(1:LevMax-LevMin+1,8) = sum(sum(o3_nlm1(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;
var_nlm2(1:LevMax-LevMin+1,8) = sum(sum(o3_nlm2(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;

var_tlm1(1:LevMax-LevMin+1,8) = sum(sum(o3_tlm1(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;
var_tlm2(1:LevMax-LevMin+1,8) = sum(sum(o3_tlm2(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;

cov_tlm1_nlm1(1:LevMax-LevMin+1,8) = sum(sum(o3_nlm1(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*o3_tlm1(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;
cov_tlm2_nlm2(1:LevMax-LevMin+1,8) = sum(sum(o3_nlm2(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*o3_tlm2(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;

corr_tlm1_nlm1(1:LevMax-LevMin+1,8) = cov_tlm1_nlm1(1:LevMax-LevMin+1,8)./sqrt(var_nlm1(1:LevMax-LevMin+1,8).*var_tlm1(1:LevMax-LevMin+1,8));
corr_tlm2_nlm2(1:LevMax-LevMin+1,8) = cov_tlm2_nlm2(1:LevMax-LevMin+1,8)./sqrt(var_nlm2(1:LevMax-LevMin+1,8).*var_tlm2(1:LevMax-LevMin+1,8));


fprintf('  Done correlation calculation \n\n\n')


lin_wid = 1.75;
fontsize = 11;
fontsize1 = 12;

figure
set(gcf,'position',[3 343 1276 1100])

subplot(1,8,1)
plot(corr_tlm1_nlm1(1:LevMax-LevMin+1,1),p_ref(LevMin+1:LevMax+1),'b','LineWidth',lin_wid)
hold on
plot(corr_tlm2_nlm2(1:LevMax-LevMin+1,1),p_ref(LevMin+1:LevMax+1),'r-','LineWidth',lin_wid)
set(gca,'YDir','reverse')
xlim([0 1])
ylim([p_ref(LevMin+1) p_ref(LevMax)])
title('u (ms^{-1})','FontSize',fontsize1,'FontName','TimesNewRoman')
ylabel('Height (hPa)','FontSize',fontsize1,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

subplot(1,8,2)
plot(corr_tlm1_nlm1(1:LevMax-LevMin+1,2),p_ref(LevMin+1:LevMax+1),'b','LineWidth',lin_wid)
hold on
plot(corr_tlm2_nlm2(1:LevMax-LevMin+1,2),p_ref(LevMin+1:LevMax+1),'r-','LineWidth',lin_wid)
set(gca,'YDir','reverse')
xlim([0 1])
ylim([p_ref(LevMin+1) p_ref(LevMax)])
title('v (ms^{-1})','FontSize',fontsize1,'FontName','TimesNewRoman')
set(gca,'YTickLabel',[])
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

subplot(1,8,3)
plot(corr_tlm1_nlm1(1:LevMax-LevMin+1,3),p_ref(LevMin+1:LevMax+1),'b','LineWidth',lin_wid)
hold on
plot(corr_tlm2_nlm2(1:LevMax-LevMin+1,3),p_ref(LevMin+1:LevMax+1),'r-','LineWidth',lin_wid)
set(gca,'YDir','reverse')
xlim([0 1])
ylim([p_ref(LevMin+1) p_ref(LevMax)])
title('T_v (K)','FontSize',fontsize1,'FontName','TimesNewRoman')
set(gca,'YTickLabel',[])
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

subplot(1,8,4)
plot(corr_tlm1_nlm1(1:LevMax-LevMin+1,4),p_ref(LevMin+1:LevMax+1),'b','LineWidth',lin_wid)
hold on
plot(corr_tlm2_nlm2(1:LevMax-LevMin+1,4),p_ref(LevMin+1:LevMax+1),'r-','LineWidth',lin_wid)
set(gca,'YDir','reverse')
xlim([0 1])
ylim([p_ref(LevMin+1) p_ref(LevMax)])
title('q (kgkg^{-1})','FontSize',fontsize1,'FontName','TimesNewRoman')
set(gca,'YTickLabel',[])
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

subplot(1,8,5)
plot(corr_tlm1_nlm1(1:LevMax-LevMin+1,5),p_ref(LevMin+1:LevMax+1),'b','LineWidth',lin_wid)
hold on
plot(corr_tlm2_nlm2(1:LevMax-LevMin+1,5),p_ref(LevMin+1:LevMax+1),'r-','LineWidth',lin_wid)
set(gca,'YDir','reverse')
xlim([0 1])
ylim([p_ref(LevMin+1) p_ref(LevMax)])
title('p (Pa)','FontSize',fontsize1,'FontName','TimesNewRoman')
set(gca,'YTickLabel',[])
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

subplot(1,8,6)
plot(corr_tlm1_nlm1(1:LevMax-LevMin+1,6),p_ref(LevMin+1:LevMax+1),'b','LineWidth',lin_wid)
hold on
plot(corr_tlm2_nlm2(1:LevMax-LevMin+1,6),p_ref(LevMin+1:LevMax+1),'r-','LineWidth',lin_wid)
set(gca,'YDir','reverse')
xlim([0 1])
ylim([p_ref(LevMin+1) p_ref(LevMax)])
title('q_l (kgkg^{-1})','FontSize',fontsize1,'FontName','TimesNewRoman')
set(gca,'YTickLabel',[])
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
% xlim([0 0.5])

subplot(1,8,7)
plot(corr_tlm1_nlm1(1:LevMax-LevMin+1,7),p_ref(LevMin+1:LevMax+1),'b','LineWidth',lin_wid)
hold on
plot(corr_tlm2_nlm2(1:LevMax-LevMin+1,7),p_ref(LevMin+1:LevMax+1),'r-','LineWidth',lin_wid)
set(gca,'YDir','reverse')
xlim([0 1])
ylim([p_ref(LevMin+1) p_ref(LevMax)])
title('q_i (kgkg^{-1})','FontSize',fontsize1,'FontName','TimesNewRoman')
set(gca,'YTickLabel',[])
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
% xlim([0 0.5])

subplot(1,8,8)
plot(corr_tlm1_nlm1(1:LevMax-LevMin+1,8),p_ref(LevMin+1:LevMax+1),'b','LineWidth',lin_wid)
hold on
plot(corr_tlm2_nlm2(1:LevMax-LevMin+1,8),p_ref(LevMin+1:LevMax+1),'r-','LineWidth',lin_wid)
set(gca,'YDir','reverse')
xlim([0 1])
ylim([p_ref(LevMin+1) p_ref(LevMax)])
title('o3 (ppm)','FontSize',fontsize1,'FontName','TimesNewRoman')
ylabel('Height (hPa)','FontSize',fontsize1,'FontName','TimesNewRoman')

set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'YAxisLocation','right')
% xlim([0 0.5])
