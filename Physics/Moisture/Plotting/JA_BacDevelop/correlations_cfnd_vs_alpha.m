close all
clear
clc
mydir = pwd;
fontsize = 16;

cd /discover/nobackup/drholdaw/ExperimentData/BacmeisterPaper/Correlations/

cfnd_vs_alpha


figure
subplot(2,1,1)
set(gcf,'position',[494   498   785   840])
scatter(cf_lsd,alpha,'kx')
box on
xlabel('C_{LS}\prime','FontSize',fontsize,'FontName','TimesNewRoman')
ylabel('\alpha','FontSize',fontsize,'FontName','TimesNewRoman')
title('(a)','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
ylim([0 0.075])


qllsd_vs_alpha


subplot(2,1,2)
set(gcf,'position',[494   498   785   840])
scatter(ql_lsd,alpha,'kx')
box on
xlabel('\partial q_{l,LS}\prime (kgkg^{-1})','FontSize',fontsize,'FontName','TimesNewRoman')
ylabel('\alpha','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
title('(b)','FontSize',fontsize,'FontName','TimesNewRoman')
ylim([0 0.075])


qcnd_vs_denom

figure
set(gcf,'position',[494   498   785   355])
scatter(qcnd,denom,'kx')
box on
xlabel('q_{LS}\prime','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')



cd(mydir)