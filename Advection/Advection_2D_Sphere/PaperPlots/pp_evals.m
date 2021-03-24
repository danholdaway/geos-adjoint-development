close all
clear
clc
mydir = pwd;

figvis = 'on';

rehash
mo_t1_maxevallwlim_smalldt
maxevallwlim = maxevals; clear maxevals
mo_t1_evalslw_T1604
e_i_lw = e_i;
e_r_lw = e_r;
mo_t1_eveclwlim_T1604
e_i_lwlim = e_i;
e_r_lwlim = e_r;
cd(mydir)

lat = -90:181/90:90;
lon = -179:2:180;

posmult = 1.0;
loweroffset = 0.00;
rightoffset = 0.00;
fontsize = 12;


C1 = figure('visible',figvis);

set(gcf,'position',[317 302 600 716])

subplot(2,4,1:3)
scatter(e_r_lw,e_i_lw,'k')
hold on
scatter(e_r_lwlim,e_i_lwlim,'kx')
box on
pos = get(gca,'position');
set(gca,'position',[pos(1)-0.02 pos(2)-0.02 posmult*pos(3) posmult*pos(4)])
ylabel('Im(\lambda)','FontSize',fontsize,'FontName','TimesNewRoman')
xlabel('Re(\lambda)','FontSize',fontsize,'FontName','TimesNewRoman')
title('(a) Eigenvalue spectrum','FontSize',fontsize,'FontName','TimesNewRoman')
legend('Lax-Wendroff','Limited Lax-Wendroff','Location','NorthEast','Orientation','horizontal')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
% xlim([-2100 50])
% ylim([-1400 1400])

subplot(2,4,4)
scatter(e_r_lw,e_i_lw,'k')
hold on
scatter(e_r_lwlim,e_i_lwlim,'kx')
box on
xlim([-50 50])
pos = get(gca,'position');
set(gca,'position',[pos(1)-0.02 pos(2)-0.02 posmult*pos(3) posmult*pos(4)])
set(gca,'YAxisLocation','right')
ylabel('Im(\lambda)','FontSize',fontsize,'FontName','TimesNewRoman')
xlabel('Re(\lambda)','FontSize',fontsize,'FontName','TimesNewRoman')
title('(b) Unstable region','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

subplot(2,4,5:8)
plot(maxevallwlim,'k')
pos = get(gca,'position');
set(gca,'position',[pos(1)-0.02 pos(2)-0.02 posmult*pos(3) posmult*pos(4)])
ylabel('Re(\lambda)','FontSize',fontsize,'FontName','TimesNewRoman')
xlabel('Timestep','FontSize',fontsize,'FontName','TimesNewRoman')
title('(c) Time series of maximum Re(\lambda)','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
ylim([-0.5 48])

pos = get(gca,'position');
set(gca,'position',[pos(1) 1.5*pos(2) pos(3) pos(4)])

cd(mydir)

asd

pause(0.1)
 
set(C1,'Units','Inches');
pos = get(C1,'Position');
set(C1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(C1,'evals_gauss.eps','-painters','-depsc2','-r300')



rehash
mo_t1_maxevallwlim_slotcyl
maxevallwlim = maxevals; clear maxevals
mo_t1_evalslw_T746_slotcyl
e_i_lw = e_i;
e_r_lw = e_r;
mo_t1_evalslwlim_T746_slotcyl
e_i_lwlim = e_i;
e_r_lwlim = e_r;
cd(mydir)


C1 = figure('visible',figvis);

set(gcf,'position',[600 302 600 716])

subplot(2,4,1:3)
scatter(e_r_lw,e_i_lw,'k')
hold on
scatter(e_r_lwlim,e_i_lwlim,'kx')
box on
pos = get(gca,'position');
set(gca,'position',[pos(1)-0.02 pos(2)-0.02 posmult*pos(3) posmult*pos(4)])
ylabel('Im(\lambda)','FontSize',fontsize,'FontName','TimesNewRoman')
xlabel('Re(\lambda)','FontSize',fontsize,'FontName','TimesNewRoman')
title('(a) Eigenvalue spectrum','FontSize',fontsize,'FontName','TimesNewRoman')
legend('Lax-Wendroff','Limited Lax-Wendroff','Location','NorthEast','Orientation','horizontal')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
xlim([-2100 50])
ylim([-1400 1400])

subplot(2,4,4)
scatter(e_r_lw,e_i_lw,'k')
hold on
scatter(e_r_lwlim,e_i_lwlim,'kx')
box on
xlim([-50 50])
pos = get(gca,'position');
set(gca,'position',[pos(1)-0.02 pos(2)-0.02 posmult*pos(3) posmult*pos(4)])
set(gca,'YAxisLocation','right')
ylabel('Im(\lambda)','FontSize',fontsize,'FontName','TimesNewRoman')
xlabel('Re(\lambda)','FontSize',fontsize,'FontName','TimesNewRoman')
title('(b) Unstable region','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

subplot(2,4,5:8)
plot(maxevallwlim,'k')
pos = get(gca,'position');
set(gca,'position',[pos(1)-0.02 pos(2)-0.02 posmult*pos(3) posmult*pos(4)])
ylabel('Re(\lambda)','FontSize',fontsize,'FontName','TimesNewRoman')
xlabel('Timestep','FontSize',fontsize,'FontName','TimesNewRoman')
title('(c) Time series of maximum Re(\lambda)','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
ylim([-0.1 12])

pos = get(gca,'position');
set(gca,'position',[pos(1) 1.5*pos(2) pos(3) pos(4)])

cd(mydir)

clearvars -except C1
pause(0.1)
 
set(C1,'Units','Inches');
pos = get(C1,'Position');
set(C1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(C1,'evals_slotc.eps','-painters','-depsc2','-r300')















