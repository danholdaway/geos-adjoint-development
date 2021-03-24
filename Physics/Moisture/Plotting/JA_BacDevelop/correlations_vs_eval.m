close all
clear
clc
mydir = pwd;
fontsize = 16;

cd /discover/nobackup/drholdaw/ExperimentData/BacmeisterPaper/Correlations/

cf_vs_eval


figure
subplot(2,1,1)
set(gcf,'position',[494   498   785   840])
scatter(cf_lsd,Me,'kx')
box on
xlabel('C_{LS}\prime','FontSize',fontsize,'FontName','TimesNewRoman')
ylabel('\lambda','FontSize',fontsize,'FontName','TimesNewRoman')
title('(a)','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')


clear Me

qlls_vs_eval0


subplot(2,1,2)
set(gcf,'position',[494   498   785   840])
scatter(ql_lsd,Me,'kx')
box on
xlabel('\partial q_{l,LS}\prime (kgkg^{-1})','FontSize',fontsize,'FontName','TimesNewRoman')
ylabel('\lambda','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
title('(b)','FontSize',fontsize,'FontName','TimesNewRoman')






cd(mydir)