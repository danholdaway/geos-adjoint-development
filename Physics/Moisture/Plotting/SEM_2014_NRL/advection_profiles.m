close all
clear
clc
mydir = pwd;

cd /home/drholdaw/Advection/Advection_1D_Offline/
outputabcd

N = 64;
x = 0:1/(N-1):1;
ymin = -0.01;
ymax =  1.01;
fontsize = 12;
lin_wid = 1.5;


figure
set(gcf,'position',[1312 218 1220 680])

subplot(2,4,1)
plot(x,q_free,'b','LineWidth',lin_wid)
hold on
plot(x,q3rd_free,'r','LineWidth',lin_wid)
xlim([0 1])
title('3rd order upwind','FontSize',fontsize,'FontName','TimesNewRoman')
ylabel('q','FontSize',fontsize,'FontName','TimesNewRoman')
legend('t = 0','t = end','Location','West')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

subplot(2,4,5)
plot(x,q_tlm,'k','LineWidth',lin_wid)
hold on
plot(x,q3rd_repl-q3rd_free,'b','LineWidth',lin_wid)
plot(x,q3rd_tlm,'r--','LineWidth',lin_wid)
xlim([0 1])
title('TLM','FontSize',fontsize,'FontName','TimesNewRoman')
xlabel('x','FontSize',fontsize,'FontName','TimesNewRoman')
ylabel('q\prime','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
legend('t=0','NLM','TLM')


subplot(2,4,2)
plot(x,q_free,'b','LineWidth',lin_wid)
hold on
plot(x,qppm_free,'r','LineWidth',lin_wid)
ylim([ymin ymax])
xlim([0 1])
title('PPM With Lin','FontSize',fontsize,'FontName','TimesNewRoman')
ylabel('q','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')


subplot(2,4,6)
plot(x,q_tlm,'k','LineWidth',lin_wid)
hold on
plot(x,qppm_repl-qppm_free,'b','LineWidth',lin_wid)
plot(x,qppm_tlm,'r--','LineWidth',lin_wid)
xlim([0 1])
title('TLM','FontSize',fontsize,'FontName','TimesNewRoman')
xlabel('x','FontSize',fontsize,'FontName','TimesNewRoman')
ylabel('q\prime','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')




subplot(2,4,3)
plot(x,q_free,'b','LineWidth',lin_wid)
hold on
plot(x,qppm4_free,'r','LineWidth',lin_wid)
ylim([ymin ymax])
xlim([0 1])
title('PPM with Colella-Woddward lim','FontSize',fontsize,'FontName','TimesNewRoman')
ylabel('q','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')


subplot(2,4,7)
plot(x,q_tlm,'k','LineWidth',lin_wid)
hold on
plot(x,qppm4_repl-qppm4_free,'b','LineWidth',lin_wid)
plot(x,qppm4_tlm,'r--','LineWidth',lin_wid)
xlim([0 1])
title('TLM','FontSize',fontsize,'FontName','TimesNewRoman')
xlabel('x','FontSize',fontsize,'FontName','TimesNewRoman')
ylabel('q\prime','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')





subplot(2,4,4)
plot(x,q_free,'b','LineWidth',lin_wid)
hold on
plot(x,qslicebs_free,'r','LineWidth',lin_wid)
ylim([ymin ymax])
xlim([0 1])
title('SLICE with Ber-Stan Limiter','FontSize',fontsize,'FontName','TimesNewRoman')
ylabel('q','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')


subplot(2,4,8)
plot(x,q_tlm,'k','LineWidth',lin_wid)
hold on
plot(x,qslicebs_repl-qslicebs_free,'b','LineWidth',lin_wid)
plot(x,qslicebs_tlm,'r--','LineWidth',lin_wid)
xlim([0 1])
title('TLM','FontSize',fontsize,'FontName','TimesNewRoman')
xlabel('x','FontSize',fontsize,'FontName','TimesNewRoman')
ylabel('q\prime','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')










cd(mydir)