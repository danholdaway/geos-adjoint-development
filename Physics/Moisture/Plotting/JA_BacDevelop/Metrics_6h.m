close all
clear
clc
mydir = pwd;

cd /home/drholdaw/LinearisedPhysics/Inputs/
pref
p = 0.5*(p_ref(1:end-1) + p_ref(2:end));

%Load Free (background) State.
cd /discover/nobackup/drholdaw/ExperimentData/BacmeisterPaper/ModelOutput/prog/free/
file = 'v000_C180.prog.eta.20140201_0900z.nc4';

lon = ncread(file,'lon');
lat = ncread(file,'lat');
lev = ncread(file,'lev');

qi_free = ncread(file,'qitot');
ql_free = ncread(file,'qltot');

%Load Perturbed (analysis) st
cd /discover/nobackup/drholdaw/ExperimentData/BacmeisterPaper/ModelOutput/prog/replay/
file = 'v000_C180.prog.eta.20140201_0900z.nc4';

qi_replay = ncread(file,'qitot');
ql_replay = ncread(file,'qltot');

%Load first Tlm state to compare.
cd /discover/nobackup/drholdaw/ExperimentData/BacmeisterPaper/ModelOutput/4DVAR
cd /discover/nobackup/drholdaw/wrk.bac/sens.20140202.000000/
% file = 'v000_C180.fvpert.eta.20140201_0900z_ST03z_P21000_S1.nc4';
file = 'v000_C180.fvpert.eta.20140201_0900z_ST03z_P21000_S1.nc4';


qi_tl1 = ncread(file,'QI');
ql_tl1 = ncread(file,'QL');

cd(mydir)

fprintf(' Done reading in the states \n\n')

fprintf('  Correlation computation \n')

im = size(qi_free,1);
jm = size(qi_free,2);
lm = size(qi_free,3);

%Compute NL perturbation trajectory.
qi_nlm = qi_replay - qi_free;
ql_nlm = ql_replay - ql_free;

LevMin = 1;

LonMin = 1;
LonMax = im;
LatMin = 1;
LatMax = jm;

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

cov_tl1 = zeros(lm-LevMin+1,8);

%Initialize correlations
corr_tl1 = zeros(lm-LevMin+1,8);

% Cloud Liquid Water
var_nlm(1:lm-LevMin+1,7) = sum(sum(qi_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
var_tl1(1:lm-LevMin+1,7) = sum(sum(qi_tl1(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
cov_tl1(1:lm-LevMin+1,7) = sum(sum(qi_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).*qi_tl1(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
corr_tl1(1:lm-LevMin+1,7) = cov_tl1(1:lm-LevMin+1,7)./sqrt(var_nlm(1:lm-LevMin+1,7).*var_tl1(1:lm-LevMin+1,7));

% Cloud Liquid Ice
var_nlm(1:lm-LevMin+1,6) = sum(sum(ql_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
var_tl1(1:lm-LevMin+1,6) = sum(sum(ql_tl1(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
cov_tl1(1:lm-LevMin+1,6) = sum(sum(ql_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).*ql_tl1(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
corr_tl1(1:lm-LevMin+1,6) = cov_tl1(1:lm-LevMin+1,6)./sqrt(var_nlm(1:lm-LevMin+1,6).*var_tl1(1:lm-LevMin+1,6));

corr_global_ql = corr_tl1(1:lm-LevMin+1,6);
corr_global_qi = corr_tl1(1:lm-LevMin+1,7);

rmsd_tl1 = zeros(lm,8);
n = (LonMax-LonMin+1)*(LatMax-LatMin+1);
%Loop over model all levels.
for k = 1:lm
    rmsd_tl1(k,6) = sqrt(  sum(sum( (ql_nlm(LonMin:LonMax,LatMin:LatMax,k) - ql_tl1(LonMin:LonMax,LatMin:LatMax,k)).^2 )) / n );
    rmsd_tl1(k,7) = sqrt(  sum(sum( (qi_nlm(LonMin:LonMax,LatMin:LatMax,k) - qi_tl1(LonMin:LonMax,LatMin:LatMax,k)).^2 )) / n );
end

rmsd_global_ql = rmsd_tl1(1:lm-LevMin+1,6);
rmsd_global_qi = rmsd_tl1(1:lm-LevMin+1,7);



rms_nlm = zeros(lm,8);
rms_tl1 = zeros(lm,8);

n = (LonMax-LonMin+1)*(LatMax-LatMin+1);

%Loop over model all levels.
for k = 1:lm
    rms_nlm(k,6) = sqrt(  sum(sum( ql_nlm(LonMin:LonMax,LatMin:LatMax,k).^2 )) / n );
    rms_tl1(k,6) = sqrt(  sum(sum( ql_tl1(LonMin:LonMax,LatMin:LatMax,k).^2 )) / n );
    
    rms_nlm(k,7) = sqrt(  sum(sum( qi_nlm(LonMin:LonMax,LatMin:LatMax,k).^2 )) / n );
    rms_tl1(k,7) = sqrt(  sum(sum( qi_tl1(LonMin:LonMax,LatMin:LatMax,k).^2 )) / n );
end

rms_global_qlnl = rms_nlm(1:lm-LevMin+1,6);
rms_global_qinl = rms_nlm(1:lm-LevMin+1,7);

rms_global_ql = rms_tl1(1:lm-LevMin+1,6);
rms_global_qi = rms_tl1(1:lm-LevMin+1,7);






LonMin = 1;
LonMax = im;
LatMin = find(lat == -24);
LatMax = find(lat ==  24);


%Variance, covariance and correlation
%Second index is number of variables.
var_nlm = zeros(lm-LevMin+1,8);
var_tl1 = zeros(lm-LevMin+1,8);

cov_tl1 = zeros(lm-LevMin+1,8);

%Initialize correlations
corr_tl1 = zeros(lm-LevMin+1,8);

% Cloud Liquid Water
var_nlm(1:lm-LevMin+1,7) = sum(sum(qi_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
var_tl1(1:lm-LevMin+1,7) = sum(sum(qi_tl1(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
cov_tl1(1:lm-LevMin+1,7) = sum(sum(qi_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).*qi_tl1(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
corr_tl1(1:lm-LevMin+1,7) = cov_tl1(1:lm-LevMin+1,7)./sqrt(var_nlm(1:lm-LevMin+1,7).*var_tl1(1:lm-LevMin+1,7));

% Cloud Liquid Ice
var_nlm(1:lm-LevMin+1,6) = sum(sum(ql_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
var_tl1(1:lm-LevMin+1,6) = sum(sum(ql_tl1(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).^2.*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
cov_tl1(1:lm-LevMin+1,6) = sum(sum(ql_nlm(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).*ql_tl1(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1).*A(LonMin:LonMax,LatMin:LatMax,1:lm-LevMin+1))) / dA;
corr_tl1(1:lm-LevMin+1,6) = cov_tl1(1:lm-LevMin+1,6)./sqrt(var_nlm(1:lm-LevMin+1,6).*var_tl1(1:lm-LevMin+1,6));

corr_tropics_ql = corr_tl1(1:lm-LevMin+1,6);
corr_tropics_qi = corr_tl1(1:lm-LevMin+1,7);

rmsd_tl1 = zeros(lm,8);
n = (LonMax-LonMin+1)*(LatMax-LatMin+1);
%Loop over model all levels.
for k = 1:lm
    rmsd_tl1(k,6) = sqrt(  sum(sum( (ql_nlm(LonMin:LonMax,LatMin:LatMax,k) - ql_tl1(LonMin:LonMax,LatMin:LatMax,k)).^2 )) / n );
    rmsd_tl1(k,7) = sqrt(  sum(sum( (qi_nlm(LonMin:LonMax,LatMin:LatMax,k) - qi_tl1(LonMin:LonMax,LatMin:LatMax,k)).^2 )) / n );
end

rmsd_tropics_ql = rmsd_tl1(1:lm-LevMin+1,6);
rmsd_tropics_qi = rmsd_tl1(1:lm-LevMin+1,7);



rms_nlm = zeros(lm,8);
rms_tl1 = zeros(lm,8);

n = (LonMax-LonMin+1)*(LatMax-LatMin+1);

%Loop over model all levels.
for k = 1:lm
    rms_nlm(k,6) = sqrt(  sum(sum( ql_nlm(LonMin:LonMax,LatMin:LatMax,k).^2 )) / n );
    rms_tl1(k,6) = sqrt(  sum(sum( ql_tl1(LonMin:LonMax,LatMin:LatMax,k).^2 )) / n );
    
    rms_nlm(k,7) = sqrt(  sum(sum( qi_nlm(LonMin:LonMax,LatMin:LatMax,k).^2 )) / n );
    rms_tl1(k,7) = sqrt(  sum(sum( qi_tl1(LonMin:LonMax,LatMin:LatMax,k).^2 )) / n );
end

rms_tropics_qlnl = rms_nlm(1:lm-LevMin+1,6);
rms_tropics_qinl = rms_nlm(1:lm-LevMin+1,7);

rms_tropics_ql = rms_tl1(1:lm-LevMin+1,6);
rms_tropics_qi = rms_tl1(1:lm-LevMin+1,7);

























grey = 0.75;
fontsize = 11;
lin_wid = 1.25;


figure
set(gcf,'position',[788   486   491   433])

subplot(1,2,1)
plot(corr_global_ql,p(LevMin:lm),'k','LineWidth',lin_wid)
hold on
plot(corr_tropics_ql,p(LevMin:lm),'Color',[grey grey grey],'LineWidth',lin_wid)

set(gca,'YDir','reverse')
ylim([p(LevMin) p(lm)])
title('(a) q_l','FontSize',fontsize,'FontName','TimesNewRoman')
ylabel('Height (hPa)','FontSize',fontsize,'FontName','TimesNewRoman')
xlabel('Correlation','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
xlim([0 0.7])
box on
legend('Global','Tropics','Location','NorthEast');
legend boxoff


subplot(1,2,2)
plot(corr_global_qi,p(LevMin:lm),'k','LineWidth',lin_wid)
hold on
plot(corr_tropics_qi,p(LevMin:lm),'Color',[grey grey grey],'LineWidth',lin_wid)

set(gca,'YDir','reverse')
ylim([p(LevMin) p(lm)])
title('(b) q_i','FontSize',fontsize,'FontName','TimesNewRoman')
ylabel('Height (hPa)','FontSize',fontsize,'FontName','TimesNewRoman')
xlabel('Correlation','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'YAxisLocation','right')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
xlim([0 0.7])
box on
pos = get(gca,'position');
set(gca,'position',[0.9*pos(1) pos(2) pos(3) pos(4)])



figure
set(gcf,'position',[788   486   491   433])

subplot(1,2,1)
plot(rmsd_global_ql,p(LevMin:lm),'k','LineWidth',lin_wid)
hold on
plot(rmsd_tropics_ql,p(LevMin:lm),'Color',[grey grey grey],'LineWidth',lin_wid)

set(gca,'YDir','reverse')
ylim([p(LevMin) p(lm)])
title('(a) q_l','FontSize',fontsize,'FontName','TimesNewRoman')
ylabel('Height (hPa)','FontSize',fontsize,'FontName','TimesNewRoman')
xlabel('RMSE (kgkg^{-1})','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
box on
legend('Global','Tropics','Location','NorthEast');
legend boxoff
xa = get(gca,'Xlim');
xlim([-xa(2)/20 xa(2)])

subplot(1,2,2)
plot(rmsd_global_qi,p(LevMin:lm),'k','LineWidth',lin_wid)
hold on
plot(rmsd_tropics_qi,p(LevMin:lm),'Color',[grey grey grey],'LineWidth',lin_wid)

set(gca,'YDir','reverse')
ylim([p(LevMin) p(lm)])
title('(b) q_i','FontSize',fontsize,'FontName','TimesNewRoman')
ylabel('Height (hPa)','FontSize',fontsize,'FontName','TimesNewRoman')
xlabel('RMSE (kgkg^{-1})','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'YAxisLocation','right')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
box on
pos = get(gca,'position');
set(gca,'position',[0.9*pos(1) pos(2) pos(3) pos(4)])
xa = get(gca,'Xlim');
xlim([-xa(2)/20 xa(2)])



figure
set(gcf,'position',[788   486   491   433])

subplot(1,2,1)
plot(rms_global_ql,p(LevMin:lm),'k','LineWidth',lin_wid)
hold on
plot(rms_global_qlnl,p(LevMin:lm),'Color',[grey grey grey],'LineWidth',lin_wid)

set(gca,'YDir','reverse')
ylim([p(LevMin) p(lm)])
title('(a) q_l','FontSize',fontsize,'FontName','TimesNewRoman')
ylabel('Height (hPa)','FontSize',fontsize,'FontName','TimesNewRoman')
xlabel('RMS (kgkg^{-1})','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
box on
legend('Linear Model','Nonlinear Model','Location','NorthEast');
legend boxoff
xa = get(gca,'Xlim');
xlim([-xa(2)/20 xa(2)])

subplot(1,2,2)
plot(rms_global_qi,p(LevMin:lm),'k','LineWidth',lin_wid)
hold on
plot(rms_global_qinl,p(LevMin:lm),'Color',[grey grey grey],'LineWidth',lin_wid)

set(gca,'YDir','reverse')
ylim([p(LevMin) p(lm)])
title('(b) q_i','FontSize',fontsize,'FontName','TimesNewRoman')
ylabel('Height (hPa)','FontSize',fontsize,'FontName','TimesNewRoman')
xlabel('RMS (kgkg^{-1})','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'YAxisLocation','right')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
box on
pos = get(gca,'position');
set(gca,'position',[0.9*pos(1) pos(2) pos(3) pos(4)])
xa = get(gca,'Xlim');
xlim([-xa(2)/20 xa(2)])
