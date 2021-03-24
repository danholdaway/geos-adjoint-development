close all
clear
clc
mydir = pwd;

%Working dir
edir = 'btmp.25956';

%Experiment name
exp = 'v000_C180';

%Choose physics options for TLM
phys = '20000';

%Choose HORD
hord = '4';

%Choose the time being considered: 16th at 12z would be 6_12
file_end = '1_1200';
% file_end = '2_0000';


%Plotting options
lin_wid = 1.5;
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
       

fprintf(' Reading in the states \n')

%Load Free (background) State.
cd /discover/nobackup/drholdaw/wrk.advec/Model_Output/hord_tr1/prog_free
file = ['v000_C180.prog.eta.2014020',file_end(1:6),'z.nc4'];

lon = ncread(file,'lon');
lat = ncread(file,'lat');
lev = ncread(file,'lev');

ql_free_hord1 = ncread(file,'qltot');

%Load Perturbed (analysis) st
cd /discover/nobackup/drholdaw/wrk.advec/Model_Output/hord_tr1/prog_replay_b_1p0001
file = ['v000_C180.prog.eta.2014020',file_end(1:6),'z.nc4'];

ql_replay_hord1 = ncread(file,'qltot');

%Load TLM state to compare.
cd /discover/nobackup/drholdaw/wrk.advec/Model_Output/hord_tr1/fvpert
file = ['v000_C180.fvpert.eta.2014020',file_end(1:6),'z_P20000.nc4'];

ql_tl_hord1 = ncread(file,'QL');



%Load Free (background) State.
cd /discover/nobackup/drholdaw/wrk.advec/Model_Output/hord_tr12/prog_free
file = ['v000_C180.prog.eta.2014020',file_end(1:6),'z.nc4'];
ql_free_hord12 = ncread(file,'qltot');

%Load Perturbed (analysis) st
cd /discover/nobackup/drholdaw/wrk.advec/Model_Output/hord_tr12/prog_replay_b_1p0001
file = ['v000_C180.prog.eta.2014020',file_end(1:6),'z.nc4'];
ql_replay_hord12 = ncread(file,'qltot');

%Load TLM state to compare.
cd /discover/nobackup/drholdaw/wrk.advec/Model_Output/hord_tr4/fvpert
file = ['v000_C180.fvpert.eta.2014020',file_end(1:6),'z_P20000.nc4'];
ql_tl_hord12 = ncread(file,'QL');

























fprintf(' Done reading in the states \n\n')

fprintf('  Correlation computation \n')

im = size(ql_free_hord1,1);
jm = size(ql_free_hord1,2);
lm = size(ql_free_hord1,3);

LonMin = 1;
LonMax = im;
LatMin = 1;
LatMax = jm;
LevMin = 1;
LevMax = lm;

cd /home/drholdaw/LinearisedPhysics/Inputs/
pref
p = 0.5*(p_ref(1:end-1)+p_ref(2:end));

%Compute NL perturbation trajectory.
ql_nlm_hord1 = ql_replay_hord1 - ql_free_hord1;
ql_nlm_hord12 = ql_replay_hord12 - ql_free_hord12;

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
var_nlm = zeros(lm-LevMin+1,2);
var_tl1 = zeros(lm-LevMin+1,2);

cov_tl1 = zeros(lm-LevMin+1,2);

%Initialize correlations
corr_tl1 = zeros(lm-LevMin+1,2);

ql_nlm = ql_nlm_hord1;
ql_tl1 = ql_tl_hord1;

% Cloud Liquid Water
var_nlm(1:lm-LevMin+1,2) = sum(sum(ql_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
var_tl1(1:lm-LevMin+1,2) = sum(sum(ql_tl1(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;

cov_tl1(1:lm-LevMin+1,2) = sum(sum(ql_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).*ql_tl1(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;

corr_tl1(1:lm-LevMin+1,2) = cov_tl1(1:lm-LevMin+1,2)./sqrt(var_nlm(1:lm-LevMin+1,2).*var_tl1(1:lm-LevMin+1,2));



figure
set(gcf,'position',[1481         363         812         514])

subplot(1,3,1)
plot(corr_tl1(1:lm-LevMin+1,2),p(LevMin:lm),'b','LineWidth',lin_wid)

xlim([0 1.01])
ylim([p(LevMin) p(lm)])
set(gca,'YDir','reverse')
title('q_l (kgkg^{-1})','FontSize',fontsize1,'FontName','TimesNewRoman')
xlabel('Correlation','FontSize',fontsize1,'FontName','TimesNewRoman')
ylabel('Pressure (hPa)','FontSize',fontsize1,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')


ql_nlm = ql_nlm_hord12;
ql_tl1 = ql_tl_hord12;

% Cloud Liquid Water
var_nlm(1:lm-LevMin+1,2) = sum(sum(ql_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
var_tl1(1:lm-LevMin+1,2) = sum(sum(ql_tl1(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;

cov_tl1(1:lm-LevMin+1,2) = sum(sum(ql_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).*ql_tl1(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;

corr_tl1(1:lm-LevMin+1,2) = cov_tl1(1:lm-LevMin+1,2)./sqrt(var_nlm(1:lm-LevMin+1,2).*var_tl1(1:lm-LevMin+1,2));




subplot(1,3,2)
plot(corr_tl1(1:lm-LevMin+1,2),p(LevMin:lm),'b','LineWidth',lin_wid)

xlim([0 1])
ylim([p(LevMin) p(lm)])
set(gca,'YDir','reverse')
title('q_l (kgkg^{-1})','FontSize',fontsize1,'FontName','TimesNewRoman')
xlabel('Correlation','FontSize',fontsize1,'FontName','TimesNewRoman')
ylabel('Pressure (hPa)','FontSize',fontsize1,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')






ql_nlm = ql_nlm_hord12;
ql_tl1 = ql_tl_hord1;

% Cloud Liquid Water
var_nlm(1:lm-LevMin+1,2) = sum(sum(ql_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
var_tl1(1:lm-LevMin+1,2) = sum(sum(ql_tl1(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;

cov_tl1(1:lm-LevMin+1,2) = sum(sum(ql_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).*ql_tl1(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;

corr_tl1(1:lm-LevMin+1,2) = cov_tl1(1:lm-LevMin+1,2)./sqrt(var_nlm(1:lm-LevMin+1,2).*var_tl1(1:lm-LevMin+1,2));



subplot(1,3,3)
plot(corr_tl1(1:lm-LevMin+1,2),p(LevMin:lm),'b','LineWidth',lin_wid)

xlim([0 1])
ylim([p(LevMin) p(lm)])
set(gca,'YDir','reverse')
title('q_l (kgkg^{-1})','FontSize',fontsize1,'FontName','TimesNewRoman')
xlabel('Correlation','FontSize',fontsize1,'FontName','TimesNewRoman')
ylabel('Pressure (hPa)','FontSize',fontsize1,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')


cd(mydir)