close all
clear
clc
mydir = pwd;

cd /home/drholdaw/LinearisedPhysics/Inputs/
pref
p = 0.5*(p_ref(1:end-1) + p_ref(2:end));

%Load Free (background) State.
cd /discover/nobackup/drholdaw/ExperimentData/BacmeisterPaper/ModelOutput/prog/free/
file = 'v000_C180.prog.eta.20140202_0300z.nc4';

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

%Load Perturbed (analysis) st
cd /discover/nobackup/drholdaw/ExperimentData/BacmeisterPaper/ModelOutput/prog/replay/
file = 'v000_C180.prog.eta.20140202_0300z.nc4';

u_replay = ncread(file,'u');
v_replay = ncread(file,'v');
t_replay = ncread(file,'tv');
q_replay = ncread(file,'sphu');
p_replay = ncread(file,'delp');
qi_replay = ncread(file,'qitot');
ql_replay = ncread(file,'qltot');

%Load first Tlm state to compare.
cd /discover/nobackup/drholdaw/wrk.bac/sens.20140202.000000/
file = 'v000_C180.fvpert.eta.20140202_0300z_ST03z_P20000_S1.nc4';

u_tl1 = ncread(file,'U');
v_tl1 = ncread(file,'V');
t_tl1 = ncread(file,'TV');
q_tl1 = ncread(file,'QV');
p_tl1 = ncread(file,'DP');
qi_tl1 = ncread(file,'QI');
ql_tl1 = ncread(file,'QL');

%Load second Tlm state to compare.
cd /discover/nobackup/drholdaw/wrk.bac/sens.20140202.000000/
file = 'v000_C180.fvpert.eta.20140202_0300z_ST03z_P22000_S1.nc4';

u_tl2 = ncread(file,'U');
v_tl2 = ncread(file,'V');
t_tl2 = ncread(file,'TV');
q_tl2 = ncread(file,'QV');
p_tl2 = ncread(file,'DP');
qi_tl2 = ncread(file,'QI');
ql_tl2 = ncread(file,'QL');

cd(mydir)

fprintf(' Done reading in the states \n\n')

fprintf('  Correlation computation \n')

im = size(u_free,1);
jm = size(u_free,2);
lm = size(u_free,3);


% LonMin = 1;
% LonMax = im;
% LatMin = 1;
% LatMax = jm;
% LevMin = 1;
% LevMax = lm;

LonMin = 1;
LonMax = im;
LatMin = find(lat == -24);
LatMax = find(lat ==  24);
LevMin = 1;
LevMax = lm;

%Compute NL perturbation trajectory.
u_nlm = u_replay - u_free;
v_nlm = v_replay - v_free;
t_nlm = t_replay - t_free;
q_nlm = q_replay - q_free;
p_nlm = p_replay - p_free;
qi_nlm = qi_replay - qi_free;
ql_nlm = ql_replay - ql_free;

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

% U winds
var_nlm(1:lm-LevMin+1,1) = sum(sum(u_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
var_tl1(1:lm-LevMin+1,1) = sum(sum(u_tl1(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
var_tl2(1:lm-LevMin+1,1) = sum(sum(u_tl2(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;

cov_tl1(1:lm-LevMin+1,1) = sum(sum(u_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).*u_tl1(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
cov_tl2(1:lm-LevMin+1,1) = sum(sum(u_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).*u_tl2(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;

corr_tl1(1:lm-LevMin+1,1) = cov_tl1(1:lm-LevMin+1,1)./sqrt(var_nlm(1:lm-LevMin+1,1).*var_tl1(1:lm-LevMin+1,1));
corr_tl2(1:lm-LevMin+1,1) = cov_tl2(1:lm-LevMin+1,1)./sqrt(var_nlm(1:lm-LevMin+1,1).*var_tl2(1:lm-LevMin+1,1));

% V winds
var_nlm(1:lm-LevMin+1,2) = sum(sum(v_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
var_tl1(1:lm-LevMin+1,2) = sum(sum(v_tl1(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
var_tl2(1:lm-LevMin+1,2) = sum(sum(v_tl2(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;

cov_tl1(1:lm-LevMin+1,2) = sum(sum(v_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).*v_tl1(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
cov_tl2(1:lm-LevMin+1,2) = sum(sum(v_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).*v_tl2(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;

corr_tl1(1:lm-LevMin+1,2) = cov_tl1(1:lm-LevMin+1,2)./sqrt(var_nlm(1:lm-LevMin+1,2).*var_tl1(1:lm-LevMin+1,2));
corr_tl2(1:lm-LevMin+1,2) = cov_tl2(1:lm-LevMin+1,2)./sqrt(var_nlm(1:lm-LevMin+1,2).*var_tl2(1:lm-LevMin+1,2));

% Temperature
var_nlm(1:lm-LevMin+1,3) = sum(sum(t_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
var_tl1(1:lm-LevMin+1,3) = sum(sum(t_tl1(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
var_tl2(1:lm-LevMin+1,3) = sum(sum(t_tl2(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;

cov_tl1(1:lm-LevMin+1,3) = sum(sum(t_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).*t_tl1(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
cov_tl2(1:lm-LevMin+1,3) = sum(sum(t_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).*t_tl2(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;

corr_tl1(1:lm-LevMin+1,3) = cov_tl1(1:lm-LevMin+1,3)./sqrt(var_nlm(1:lm-LevMin+1,3).*var_tl1(1:lm-LevMin+1,3));
corr_tl2(1:lm-LevMin+1,3) = cov_tl2(1:lm-LevMin+1,3)./sqrt(var_nlm(1:lm-LevMin+1,3).*var_tl2(1:lm-LevMin+1,3));

% Specific Humidity
var_nlm(1:lm-LevMin+1,4) = sum(sum(q_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
var_tl1(1:lm-LevMin+1,4) = sum(sum(q_tl1(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
var_tl2(1:lm-LevMin+1,4) = sum(sum(q_tl2(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;

cov_tl1(1:lm-LevMin+1,4) = sum(sum(q_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).*q_tl1(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
cov_tl2(1:lm-LevMin+1,4) = sum(sum(q_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).*q_tl2(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;

corr_tl1(1:lm-LevMin+1,4) = cov_tl1(1:lm-LevMin+1,4)./sqrt(var_nlm(1:lm-LevMin+1,4).*var_tl1(1:lm-LevMin+1,4));
corr_tl2(1:lm-LevMin+1,4) = cov_tl2(1:lm-LevMin+1,4)./sqrt(var_nlm(1:lm-LevMin+1,4).*var_tl2(1:lm-LevMin+1,4));

% Pressure
var_nlm(1:lm-LevMin+1,5) = sum(sum(p_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
var_tl1(1:lm-LevMin+1,5) = sum(sum(p_tl1(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
var_tl2(1:lm-LevMin+1,5) = sum(sum(p_tl2(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;

cov_tl1(1:lm-LevMin+1,5) = sum(sum(p_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).*p_tl1(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
cov_tl2(1:lm-LevMin+1,5) = sum(sum(p_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).*p_tl2(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;

corr_tl1(1:lm-LevMin+1,5) = cov_tl1(1:lm-LevMin+1,5)./sqrt(var_nlm(1:lm-LevMin+1,5).*var_tl1(1:lm-LevMin+1,5));
corr_tl2(1:lm-LevMin+1,5) = cov_tl2(1:lm-LevMin+1,5)./sqrt(var_nlm(1:lm-LevMin+1,5).*var_tl2(1:lm-LevMin+1,5));

% Cloud Liquid Water
var_nlm(1:lm-LevMin+1,6) = sum(sum(qi_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
var_tl1(1:lm-LevMin+1,6) = sum(sum(qi_tl1(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
var_tl2(1:lm-LevMin+1,6) = sum(sum(qi_tl2(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;

cov_tl1(1:lm-LevMin+1,6) = sum(sum(qi_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).*qi_tl1(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
cov_tl2(1:lm-LevMin+1,6) = sum(sum(qi_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).*qi_tl2(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;

corr_tl1(1:lm-LevMin+1,6) = cov_tl1(1:lm-LevMin+1,6)./sqrt(var_nlm(1:lm-LevMin+1,6).*var_tl1(1:lm-LevMin+1,6));
corr_tl2(1:lm-LevMin+1,6) = cov_tl2(1:lm-LevMin+1,6)./sqrt(var_nlm(1:lm-LevMin+1,6).*var_tl2(1:lm-LevMin+1,6));

% Cloud Liquid Ice
var_nlm(1:lm-LevMin+1,7) = sum(sum(ql_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
var_tl1(1:lm-LevMin+1,7) = sum(sum(ql_tl1(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
var_tl2(1:lm-LevMin+1,7) = sum(sum(ql_tl2(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;

cov_tl1(1:lm-LevMin+1,7) = sum(sum(ql_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).*ql_tl1(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
cov_tl2(1:lm-LevMin+1,7) = sum(sum(ql_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).*ql_tl2(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;

corr_tl1(1:lm-LevMin+1,7) = cov_tl1(1:lm-LevMin+1,7)./sqrt(var_nlm(1:lm-LevMin+1,7).*var_tl1(1:lm-LevMin+1,7));
corr_tl2(1:lm-LevMin+1,7) = cov_tl2(1:lm-LevMin+1,7)./sqrt(var_nlm(1:lm-LevMin+1,7).*var_tl2(1:lm-LevMin+1,7));

fprintf('  Done correlation computation \n\n')


fontsize = 14;
lin_wid = 1.5;


figure
set(gcf,'position',[3 343 1276 576])

subplot(1,7,1)
plot(corr_tl1(1:lm-LevMin+1,1),p(LevMin:lm),'b','LineWidth',lin_wid)
hold on
plot(corr_tl2(1:lm-LevMin+1,1),p(LevMin:lm),'r','LineWidth',lin_wid)
set(gca,'YDir','reverse')
xlim([0 1])
ylim([p(LevMin) p(lm)])
title('(a) u\prime','FontSize',fontsize,'FontName','TimesNewRoman')
ylabel('Height (hPa)','FontSize',fontsize,'FontName','TimesNewRoman')
xlabel('Correlation','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

subplot(1,7,2)
plot(corr_tl1(1:lm-LevMin+1,2),p(LevMin:lm),'b','LineWidth',lin_wid)
hold on
plot(corr_tl2(1:lm-LevMin+1,2),p(LevMin:lm),'r','LineWidth',lin_wid)
set(gca,'YDir','reverse')
xlim([0 1])
ylim([p(LevMin) p(lm)])
title('(b) v\prime','FontSize',fontsize,'FontName','TimesNewRoman')
xlabel('Correlation','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'YTickLabel',[])
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

subplot(1,7,3)
plot(corr_tl1(1:lm-LevMin+1,3),p(LevMin:lm),'b','LineWidth',lin_wid)
hold on
plot(corr_tl2(1:lm-LevMin+1,3),p(LevMin:lm),'r','LineWidth',lin_wid)
set(gca,'YDir','reverse')
xlim([0 1])
ylim([p(LevMin) p(lm)])
title('(c) T_v\prime','FontSize',fontsize,'FontName','TimesNewRoman')
xlabel('Correlation','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'YTickLabel',[])
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

subplot(1,7,4)
plot(corr_tl1(1:lm-LevMin+1,4),p(LevMin:lm),'b','LineWidth',lin_wid)
hold on
plot(corr_tl2(1:lm-LevMin+1,4),p(LevMin:lm),'r','LineWidth',lin_wid)
set(gca,'YDir','reverse')
xlim([0 1])
ylim([p(LevMin) p(lm)])
title('(d) q\prime','FontSize',fontsize,'FontName','TimesNewRoman')
xlabel('Correlation','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'YTickLabel',[])
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
xlim([0 1.0])

subplot(1,7,5)
plot(corr_tl1(1:lm-LevMin+1,7),p(LevMin:lm),'b','LineWidth',lin_wid)
hold on
plot(corr_tl2(1:lm-LevMin+1,7),p(LevMin:lm),'r','LineWidth',lin_wid)
set(gca,'YDir','reverse')
xlim([0 1])
ylim([p(LevMin) p(lm)])
title('(e) q_l\prime','FontSize',fontsize,'FontName','TimesNewRoman')
xlabel('Correlation','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'YTickLabel',[])
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
xlim([0 0.5])

subplot(1,7,6)
plot(corr_tl1(1:lm-LevMin+1,6),p(LevMin:lm),'b','LineWidth',lin_wid)
hold on
plot(corr_tl2(1:lm-LevMin+1,6),p(LevMin:lm),'r','LineWidth',lin_wid)
set(gca,'YDir','reverse')
xlim([0 1])
ylim([p(LevMin) p(lm)])
title('(f) q_i\prime','FontSize',fontsize,'FontName','TimesNewRoman')
xlabel('Correlation','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'YTickLabel',[])
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
xlim([0 0.5])

subplot(1,7,7)
plot(corr_tl1(1:lm-LevMin+1,5),p(LevMin:lm),'b','LineWidth',lin_wid)
hold on
plot(corr_tl2(1:lm-LevMin+1,5),p(LevMin:lm),'r','LineWidth',lin_wid)
set(gca,'YDir','reverse')
xlim([0 1])
ylim([p(LevMin) p(lm)])
title('(g) p\prime','FontSize',fontsize,'FontName','TimesNewRoman')
xlabel('Correlation','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'YAxisLocation','right')
ylabel('Height (hPa)','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
legend('Dry Model','With Moist','Location','NorthEast')




