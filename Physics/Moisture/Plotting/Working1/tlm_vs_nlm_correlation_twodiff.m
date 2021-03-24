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
file_end = '2_0000';
% file_end = '2_0300';

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
d1_free1 = ncread(file,'DU001',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
d2_free1 = ncread(file,'DU002',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
d3_free1 = ncread(file,'DU003',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
d4_free1 = ncread(file,'DU004',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
d5_free1 = ncread(file,'DU005',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);

%Load Perturbed (analysis) st
cd /discover/nobackup/drholdaw/btmp.25956/prog/prog_replay1/
file = ['v000_C180.prog.eta.2014020',file_end,'z.nc4'];

u_replay1 = ncread(file,'u',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
v_replay1 = ncread(file,'v',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
t_replay1 = ncread(file,'tv',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
q_replay1 = ncread(file,'sphu',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
p_replay1 = ncread(file,'delp',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
qi_replay1 = ncread(file,'qitot',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
ql_replay1 = ncread(file,'qltot',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
o3_replay1 = ncread(file,'ozone',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
d1_replay1 = ncread(file,'DU001',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
d2_replay1 = ncread(file,'DU002',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
d3_replay1 = ncread(file,'DU003',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
d4_replay1 = ncread(file,'DU004',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
d5_replay1 = ncread(file,'DU005',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);

%Load first TLM state to compare.
cd /discover/nobackup/drholdaw/btmp.25956/sens.20140202.000000/
file = ['v000_C180.fvpert.eta.2014020',file_end,'z_ST03z_DUST_P20000_S0.nc4'];
% file = 'fvpertX.eta.nc4';

u_tl1 = ncread(file,'U',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
v_tl1 = ncread(file,'V',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
t_tl1 = ncread(file,'TV',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
q_tl1 = ncread(file,'QV',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
p_tl1 = ncread(file,'DP',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
qi_tl1 = ncread(file,'QI',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
ql_tl1 = ncread(file,'QL',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
o3_tl1 = ncread(file,'O3',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
d1_tl1 = ncread(file,'DU001',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
d2_tl1 = ncread(file,'DU002',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
d3_tl1 = ncread(file,'DU003',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
d4_tl1 = ncread(file,'DU004',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
d5_tl1 = ncread(file,'DU005',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);

%Load second TLM state to compare.
cd /discover/nobackup/drholdaw/btmp.25956/sens.20140202.000000/
file = ['v000_C180.fvpert.eta.2014020',file_end,'z_ST03z_DUST_P22000_S0.nc4'];
% file = 'fvpertX.eta.nc4';

u_tl2 = ncread(file,'U',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
v_tl2 = ncread(file,'V',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
t_tl2 = ncread(file,'TV',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
q_tl2 = ncread(file,'QV',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
p_tl2 = ncread(file,'DP',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
qi_tl2 = ncread(file,'QI',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
ql_tl2 = ncread(file,'QL',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
o3_tl2 = ncread(file,'O3',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
d1_tl2 = ncread(file,'DU001',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
d2_tl2 = ncread(file,'DU002',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
d3_tl2 = ncread(file,'DU003',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
d4_tl2 = ncread(file,'DU004',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);
d5_tl2 = ncread(file,'DU005',[1 1 LevMin 1],[Inf Inf LevMax-LevMin+1 Inf]);

cd(mydir)

%Compute NL perturbation trajectory.
u_nlm = u_replay1 - u_free1;
v_nlm = v_replay1 - v_free1;
t_nlm = t_replay1 - t_free1;
q_nlm = q_replay1 - q_free1;
p_nlm = p_replay1 - p_free1;
qi_nlm = qi_replay1 - qi_free1;
ql_nlm = ql_replay1 - ql_free1;
o3_nlm = o3_replay1 - o3_free1;
d1_nlm = d1_replay1 - d1_free1;
d2_nlm = d2_replay1 - d2_free1;
d3_nlm = d3_replay1 - d3_free1;
d4_nlm = d4_replay1 - d4_free1;
d5_nlm = d5_replay1 - d5_free1;

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
var_nlm = zeros(LevMax-LevMin+1,13);
var_tl1 = zeros(LevMax-LevMin+1,13);
var_tl2 = zeros(LevMax-LevMin+1,13);

cov_tl1 = zeros(LevMax-LevMin+1,13);
cov_tl2 = zeros(LevMax-LevMin+1,13);

%Initialize correlations
corr_tl1 = zeros(LevMax-LevMin+1,13);
corr_tl2 = zeros(LevMax-LevMin+1,13);
corr_tl3 = zeros(LevMax-LevMin+1,13);
corr_tl4 = zeros(LevMax-LevMin+1,13);
corr_tl5 = zeros(LevMax-LevMin+1,13);
corr_tl6 = zeros(LevMax-LevMin+1,13);

% U winds
var_nlm(1:LevMax-LevMin+1,1) = sum(sum(u_nlm(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;
var_tl1(1:LevMax-LevMin+1,1) = sum(sum(u_tl1(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;
var_tl2(1:LevMax-LevMin+1,1) = sum(sum(u_tl2(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;

cov_tl1(1:LevMax-LevMin+1,1) = sum(sum(u_nlm(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*u_tl1(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;
cov_tl2(1:LevMax-LevMin+1,1) = sum(sum(u_nlm(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*u_tl2(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;

corr_tl1(1:LevMax-LevMin+1,1) = cov_tl1(1:LevMax-LevMin+1,1)./sqrt(var_nlm(1:LevMax-LevMin+1,1).*var_tl1(1:LevMax-LevMin+1,1));
corr_tl2(1:LevMax-LevMin+1,1) = cov_tl2(1:LevMax-LevMin+1,1)./sqrt(var_nlm(1:LevMax-LevMin+1,1).*var_tl2(1:LevMax-LevMin+1,1));

% V winds
var_nlm(1:LevMax-LevMin+1,2) = sum(sum(v_nlm(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;
var_tl1(1:LevMax-LevMin+1,2) = sum(sum(v_tl1(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;
var_tl2(1:LevMax-LevMin+1,2) = sum(sum(v_tl2(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;

cov_tl1(1:LevMax-LevMin+1,2) = sum(sum(v_nlm(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*v_tl1(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;
cov_tl2(1:LevMax-LevMin+1,2) = sum(sum(v_nlm(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*v_tl2(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;

corr_tl1(1:LevMax-LevMin+1,2) = cov_tl1(1:LevMax-LevMin+1,2)./sqrt(var_nlm(1:LevMax-LevMin+1,2).*var_tl1(1:LevMax-LevMin+1,2));
corr_tl2(1:LevMax-LevMin+1,2) = cov_tl2(1:LevMax-LevMin+1,2)./sqrt(var_nlm(1:LevMax-LevMin+1,2).*var_tl2(1:LevMax-LevMin+1,2));

% Temperature
var_nlm(1:LevMax-LevMin+1,3) = sum(sum(t_nlm(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;
var_tl1(1:LevMax-LevMin+1,3) = sum(sum(t_tl1(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;
var_tl2(1:LevMax-LevMin+1,3) = sum(sum(t_tl2(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;

cov_tl1(1:LevMax-LevMin+1,3) = sum(sum(t_nlm(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*t_tl1(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;
cov_tl2(1:LevMax-LevMin+1,3) = sum(sum(t_nlm(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*t_tl2(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;

corr_tl1(1:LevMax-LevMin+1,3) = cov_tl1(1:LevMax-LevMin+1,3)./sqrt(var_nlm(1:LevMax-LevMin+1,3).*var_tl1(1:LevMax-LevMin+1,3));
corr_tl2(1:LevMax-LevMin+1,3) = cov_tl2(1:LevMax-LevMin+1,3)./sqrt(var_nlm(1:LevMax-LevMin+1,3).*var_tl2(1:LevMax-LevMin+1,3));

% Specific Humidity
var_nlm(1:LevMax-LevMin+1,4) = sum(sum(q_nlm(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;
var_tl1(1:LevMax-LevMin+1,4) = sum(sum(q_tl1(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;
var_tl2(1:LevMax-LevMin+1,4) = sum(sum(q_tl2(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;

cov_tl1(1:LevMax-LevMin+1,4) = sum(sum(q_nlm(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*q_tl1(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;
cov_tl2(1:LevMax-LevMin+1,4) = sum(sum(q_nlm(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*q_tl2(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;

corr_tl1(1:LevMax-LevMin+1,4) = cov_tl1(1:LevMax-LevMin+1,4)./sqrt(var_nlm(1:LevMax-LevMin+1,4).*var_tl1(1:LevMax-LevMin+1,4));
corr_tl2(1:LevMax-LevMin+1,4) = cov_tl2(1:LevMax-LevMin+1,4)./sqrt(var_nlm(1:LevMax-LevMin+1,4).*var_tl2(1:LevMax-LevMin+1,4));

% Pressure
var_nlm(1:LevMax-LevMin+1,5) = sum(sum(p_nlm(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;
var_tl1(1:LevMax-LevMin+1,5) = sum(sum(p_tl1(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;
var_tl2(1:LevMax-LevMin+1,5) = sum(sum(p_tl2(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;

cov_tl1(1:LevMax-LevMin+1,5) = sum(sum(p_nlm(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*p_tl1(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;
cov_tl2(1:LevMax-LevMin+1,5) = sum(sum(p_nlm(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*p_tl2(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;

corr_tl1(1:LevMax-LevMin+1,5) = cov_tl1(1:LevMax-LevMin+1,5)./sqrt(var_nlm(1:LevMax-LevMin+1,5).*var_tl1(1:LevMax-LevMin+1,5));
corr_tl2(1:LevMax-LevMin+1,5) = cov_tl2(1:LevMax-LevMin+1,5)./sqrt(var_nlm(1:LevMax-LevMin+1,5).*var_tl2(1:LevMax-LevMin+1,5));

% Cloud Liquid Water
var_nlm(1:LevMax-LevMin+1,6) = sum(sum(ql_nlm(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;
var_tl1(1:LevMax-LevMin+1,6) = sum(sum(ql_tl1(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;
var_tl2(1:LevMax-LevMin+1,6) = sum(sum(ql_tl2(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;

cov_tl1(1:LevMax-LevMin+1,6) = sum(sum(ql_nlm(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*ql_tl1(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;
cov_tl2(1:LevMax-LevMin+1,6) = sum(sum(ql_nlm(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*ql_tl2(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;

corr_tl1(1:LevMax-LevMin+1,6) = cov_tl1(1:LevMax-LevMin+1,6)./sqrt(var_nlm(1:LevMax-LevMin+1,6).*var_tl1(1:LevMax-LevMin+1,6));
corr_tl2(1:LevMax-LevMin+1,6) = cov_tl2(1:LevMax-LevMin+1,6)./sqrt(var_nlm(1:LevMax-LevMin+1,6).*var_tl2(1:LevMax-LevMin+1,6));

% Cloud Liquid Ice
var_nlm(1:LevMax-LevMin+1,7) = sum(sum(qi_nlm(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;
var_tl1(1:LevMax-LevMin+1,7) = sum(sum(qi_tl1(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;
var_tl2(1:LevMax-LevMin+1,7) = sum(sum(qi_tl2(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;

cov_tl1(1:LevMax-LevMin+1,7) = sum(sum(qi_nlm(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*qi_tl1(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;
cov_tl2(1:LevMax-LevMin+1,7) = sum(sum(qi_nlm(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*qi_tl2(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;

corr_tl1(1:LevMax-LevMin+1,7) = cov_tl1(1:LevMax-LevMin+1,7)./sqrt(var_nlm(1:LevMax-LevMin+1,7).*var_tl1(1:LevMax-LevMin+1,7));
corr_tl2(1:LevMax-LevMin+1,7) = cov_tl2(1:LevMax-LevMin+1,7)./sqrt(var_nlm(1:LevMax-LevMin+1,7).*var_tl2(1:LevMax-LevMin+1,7));

% Ozone
var_nlm(1:LevMax-LevMin+1,8) = sum(sum(o3_nlm(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;
var_tl1(1:LevMax-LevMin+1,8) = sum(sum(o3_tl1(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;
var_tl2(1:LevMax-LevMin+1,8) = sum(sum(o3_tl2(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;

cov_tl1(1:LevMax-LevMin+1,8) = sum(sum(o3_nlm(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*o3_tl1(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;
cov_tl2(1:LevMax-LevMin+1,8) = sum(sum(o3_nlm(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*o3_tl2(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;

corr_tl1(1:LevMax-LevMin+1,8) = cov_tl1(1:LevMax-LevMin+1,8)./sqrt(var_nlm(1:LevMax-LevMin+1,8).*var_tl1(1:LevMax-LevMin+1,8));
corr_tl2(1:LevMax-LevMin+1,8) = cov_tl2(1:LevMax-LevMin+1,8)./sqrt(var_nlm(1:LevMax-LevMin+1,8).*var_tl2(1:LevMax-LevMin+1,8));

% Dust bin 1
var_nlm(1:LevMax-LevMin+1,9) = sum(sum(d1_nlm(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;
var_tl1(1:LevMax-LevMin+1,9) = sum(sum(d1_tl1(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;
var_tl2(1:LevMax-LevMin+1,9) = sum(sum(d1_tl2(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;

cov_tl1(1:LevMax-LevMin+1,9) = sum(sum(d1_nlm(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*d1_tl1(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;
cov_tl2(1:LevMax-LevMin+1,9) = sum(sum(d1_nlm(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*d1_tl2(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;

corr_tl1(1:LevMax-LevMin+1,9) = cov_tl1(1:LevMax-LevMin+1,9)./sqrt(var_nlm(1:LevMax-LevMin+1,9).*var_tl1(1:LevMax-LevMin+1,9));
corr_tl2(1:LevMax-LevMin+1,9) = cov_tl2(1:LevMax-LevMin+1,9)./sqrt(var_nlm(1:LevMax-LevMin+1,9).*var_tl2(1:LevMax-LevMin+1,9));

% Dust bin 2
var_nlm(1:LevMax-LevMin+1,10) = sum(sum(d2_nlm(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;
var_tl1(1:LevMax-LevMin+1,10) = sum(sum(d2_tl1(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;
var_tl2(1:LevMax-LevMin+1,10) = sum(sum(d2_tl2(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;

cov_tl1(1:LevMax-LevMin+1,10) = sum(sum(d2_nlm(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*d2_tl1(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;
cov_tl2(1:LevMax-LevMin+1,10) = sum(sum(d2_nlm(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*d2_tl2(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;

corr_tl1(1:LevMax-LevMin+1,10) = cov_tl1(1:LevMax-LevMin+1,10)./sqrt(var_nlm(1:LevMax-LevMin+1,10).*var_tl1(1:LevMax-LevMin+1,10));
corr_tl2(1:LevMax-LevMin+1,10) = cov_tl2(1:LevMax-LevMin+1,10)./sqrt(var_nlm(1:LevMax-LevMin+1,10).*var_tl2(1:LevMax-LevMin+1,10));

% Dust bin 3
var_nlm(1:LevMax-LevMin+1,11) = sum(sum(d3_nlm(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;
var_tl1(1:LevMax-LevMin+1,11) = sum(sum(d3_tl1(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;
var_tl2(1:LevMax-LevMin+1,11) = sum(sum(d3_tl2(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;

cov_tl1(1:LevMax-LevMin+1,11) = sum(sum(d3_nlm(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*d3_tl1(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;
cov_tl2(1:LevMax-LevMin+1,11) = sum(sum(d3_nlm(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*d3_tl2(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;

corr_tl1(1:LevMax-LevMin+1,11) = cov_tl1(1:LevMax-LevMin+1,11)./sqrt(var_nlm(1:LevMax-LevMin+1,11).*var_tl1(1:LevMax-LevMin+1,11));
corr_tl2(1:LevMax-LevMin+1,11) = cov_tl2(1:LevMax-LevMin+1,11)./sqrt(var_nlm(1:LevMax-LevMin+1,11).*var_tl2(1:LevMax-LevMin+1,11));

% Dust bin 4
var_nlm(1:LevMax-LevMin+1,12) = sum(sum(d4_nlm(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;
var_tl1(1:LevMax-LevMin+1,12) = sum(sum(d4_tl1(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;
var_tl2(1:LevMax-LevMin+1,12) = sum(sum(d4_tl2(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;

cov_tl1(1:LevMax-LevMin+1,12) = sum(sum(d4_nlm(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*d4_tl1(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;
cov_tl2(1:LevMax-LevMin+1,12) = sum(sum(d4_nlm(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*d4_tl2(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;

corr_tl1(1:LevMax-LevMin+1,12) = cov_tl1(1:LevMax-LevMin+1,12)./sqrt(var_nlm(1:LevMax-LevMin+1,12).*var_tl1(1:LevMax-LevMin+1,12));
corr_tl2(1:LevMax-LevMin+1,12) = cov_tl2(1:LevMax-LevMin+1,12)./sqrt(var_nlm(1:LevMax-LevMin+1,12).*var_tl2(1:LevMax-LevMin+1,12));

% Dust bin 5
var_nlm(1:LevMax-LevMin+1,13) = sum(sum(d5_nlm(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;
var_tl1(1:LevMax-LevMin+1,13) = sum(sum(d5_tl1(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;
var_tl2(1:LevMax-LevMin+1,13) = sum(sum(d5_tl2(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;

cov_tl1(1:LevMax-LevMin+1,13) = sum(sum(d5_nlm(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*d5_tl1(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;
cov_tl2(1:LevMax-LevMin+1,13) = sum(sum(d5_nlm(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*d5_tl2(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:LevMax-LevMin+1))) / dA;

corr_tl1(1:LevMax-LevMin+1,13) = cov_tl1(1:LevMax-LevMin+1,13)./sqrt(var_nlm(1:LevMax-LevMin+1,13).*var_tl1(1:LevMax-LevMin+1,13));
corr_tl2(1:LevMax-LevMin+1,13) = cov_tl2(1:LevMax-LevMin+1,13)./sqrt(var_nlm(1:LevMax-LevMin+1,13).*var_tl2(1:LevMax-LevMin+1,13));




fprintf('  Done correlation calculation \n\n\n')


lin_wid = 1.75;
fontsize = 11;
fontsize1 = 12;

figure
set(gcf,'position',[3 343 1276 1100])

subplot(2,8,1)
plot(corr_tl1(1:LevMax-LevMin+1,1),p_ref(LevMin+1:LevMax+1),'b','LineWidth',lin_wid)
hold on
plot(corr_tl2(1:LevMax-LevMin+1,1),p_ref(LevMin+1:LevMax+1),'r-','LineWidth',lin_wid)
set(gca,'YDir','reverse')
xlim([0 1])
ylim([p_ref(LevMin+1) p_ref(LevMax)])
title('u (ms^{-1})','FontSize',fontsize1,'FontName','TimesNewRoman')
ylabel('Height (hPa)','FontSize',fontsize1,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

subplot(2,8,2)
plot(corr_tl1(1:LevMax-LevMin+1,2),p_ref(LevMin+1:LevMax+1),'b','LineWidth',lin_wid)
hold on
plot(corr_tl2(1:LevMax-LevMin+1,2),p_ref(LevMin+1:LevMax+1),'r-','LineWidth',lin_wid)
set(gca,'YDir','reverse')
xlim([0 1])
ylim([p_ref(LevMin+1) p_ref(LevMax)])
title('v (ms^{-1})','FontSize',fontsize1,'FontName','TimesNewRoman')
set(gca,'YTickLabel',[])
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

subplot(2,8,3)
plot(corr_tl1(1:LevMax-LevMin+1,3),p_ref(LevMin+1:LevMax+1),'b','LineWidth',lin_wid)
hold on
plot(corr_tl2(1:LevMax-LevMin+1,3),p_ref(LevMin+1:LevMax+1),'r-','LineWidth',lin_wid)
set(gca,'YDir','reverse')
xlim([0 1])
ylim([p_ref(LevMin+1) p_ref(LevMax)])
title('T_v (K)','FontSize',fontsize1,'FontName','TimesNewRoman')
set(gca,'YTickLabel',[])
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

subplot(2,8,4)
plot(corr_tl1(1:LevMax-LevMin+1,4),p_ref(LevMin+1:LevMax+1),'b','LineWidth',lin_wid)
hold on
plot(corr_tl2(1:LevMax-LevMin+1,4),p_ref(LevMin+1:LevMax+1),'r-','LineWidth',lin_wid)
set(gca,'YDir','reverse')
xlim([0 1])
ylim([p_ref(LevMin+1) p_ref(LevMax)])
title('q (kgkg^{-1})','FontSize',fontsize1,'FontName','TimesNewRoman')
set(gca,'YTickLabel',[])
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

subplot(2,8,5)
plot(corr_tl1(1:LevMax-LevMin+1,5),p_ref(LevMin+1:LevMax+1),'b','LineWidth',lin_wid)
hold on
plot(corr_tl2(1:LevMax-LevMin+1,5),p_ref(LevMin+1:LevMax+1),'r-','LineWidth',lin_wid)
set(gca,'YDir','reverse')
xlim([0 1])
ylim([p_ref(LevMin+1) p_ref(LevMax)])
title('p (Pa)','FontSize',fontsize1,'FontName','TimesNewRoman')
set(gca,'YTickLabel',[])
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

subplot(2,8,6)
plot(corr_tl1(1:LevMax-LevMin+1,6),p_ref(LevMin+1:LevMax+1),'b','LineWidth',lin_wid)
hold on
plot(corr_tl2(1:LevMax-LevMin+1,6),p_ref(LevMin+1:LevMax+1),'r-','LineWidth',lin_wid)
set(gca,'YDir','reverse')
xlim([0 1])
ylim([p_ref(LevMin+1) p_ref(LevMax)])
title('q_l (kgkg^{-1})','FontSize',fontsize1,'FontName','TimesNewRoman')
set(gca,'YTickLabel',[])
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
% xlim([0 0.5])

subplot(2,8,7)
plot(corr_tl1(1:LevMax-LevMin+1,7),p_ref(LevMin+1:LevMax+1),'b','LineWidth',lin_wid)
hold on
plot(corr_tl2(1:LevMax-LevMin+1,7),p_ref(LevMin+1:LevMax+1),'r-','LineWidth',lin_wid)
set(gca,'YDir','reverse')
xlim([0 1])
ylim([p_ref(LevMin+1) p_ref(LevMax)])
title('q_i (kgkg^{-1})','FontSize',fontsize1,'FontName','TimesNewRoman')
set(gca,'YTickLabel',[])
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
% xlim([0 0.5])

subplot(2,8,8)
plot(corr_tl1(1:LevMax-LevMin+1,8),p_ref(LevMin+1:LevMax+1),'b','LineWidth',lin_wid)
hold on
plot(corr_tl2(1:LevMax-LevMin+1,8),p_ref(LevMin+1:LevMax+1),'r-','LineWidth',lin_wid)
set(gca,'YDir','reverse')
xlim([0 1])
ylim([p_ref(LevMin+1) p_ref(LevMax)])
title('o3 (ppm)','FontSize',fontsize1,'FontName','TimesNewRoman')
ylabel('Height (hPa)','FontSize',fontsize1,'FontName','TimesNewRoman')

set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'YAxisLocation','right')
% xlim([0 0.5])

subplot(2,8,9)
plot(corr_tl1(1:LevMax-LevMin+1,9),p_ref(LevMin+1:LevMax+1),'b','LineWidth',lin_wid)
hold on
plot(corr_tl2(1:LevMax-LevMin+1,9),p_ref(LevMin+1:LevMax+1),'r-','LineWidth',lin_wid)
set(gca,'YDir','reverse')
% xlim([0 1])
ylim([p_ref(LevMin+1) p_ref(LevMax)])
title('DU B1 (kgkg^{-1})','FontSize',fontsize1,'FontName','TimesNewRoman')
ylabel('Height (hPa)','FontSize',fontsize1,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

subplot(2,8,10)
plot(corr_tl1(1:LevMax-LevMin+1,10),p_ref(LevMin+1:LevMax+1),'b','LineWidth',lin_wid)
hold on
plot(corr_tl2(1:LevMax-LevMin+1,10),p_ref(LevMin+1:LevMax+1),'r-','LineWidth',lin_wid)
set(gca,'YDir','reverse')
% xlim([0 1])
ylim([p_ref(LevMin+1) p_ref(LevMax)])
title('DU B2 (kgkg^{-1})','FontSize',fontsize1,'FontName','TimesNewRoman')
set(gca,'YTickLabel',[])
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

subplot(2,8,11)
plot(corr_tl1(1:LevMax-LevMin+1,11),p_ref(LevMin+1:LevMax+1),'b','LineWidth',lin_wid)
hold on
plot(corr_tl2(1:LevMax-LevMin+1,11),p_ref(LevMin+1:LevMax+1),'r-','LineWidth',lin_wid)
set(gca,'YDir','reverse')
% xlim([0 1])
ylim([p_ref(LevMin+1) p_ref(LevMax)])
title('DU B3 (kgkg^{-1})','FontSize',fontsize1,'FontName','TimesNewRoman')
set(gca,'YTickLabel',[])
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

subplot(2,8,12)
plot(corr_tl1(1:LevMax-LevMin+1,12),p_ref(LevMin+1:LevMax+1),'b','LineWidth',lin_wid)
hold on
plot(corr_tl2(1:LevMax-LevMin+1,12),p_ref(LevMin+1:LevMax+1),'r-','LineWidth',lin_wid)
set(gca,'YDir','reverse')
% xlim([0 1])
ylim([p_ref(LevMin+1) p_ref(LevMax)])
title('DU B4 (kgkg^{-1})','FontSize',fontsize1,'FontName','TimesNewRoman')
set(gca,'YTickLabel',[])
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

subplot(2,8,13)
plot(corr_tl1(1:LevMax-LevMin+1,13),p_ref(LevMin+1:LevMax+1),'b','LineWidth',lin_wid)
hold on
plot(corr_tl2(1:LevMax-LevMin+1,13),p_ref(LevMin+1:LevMax+1),'r-','LineWidth',lin_wid)
set(gca,'YDir','reverse')
% xlim([0 1])
ylim([p_ref(LevMin+1) p_ref(LevMax)])
title('DU B5 (kgkg^{-1})','FontSize',fontsize1,'FontName','TimesNewRoman')
set(gca,'YTickLabel',[])
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

