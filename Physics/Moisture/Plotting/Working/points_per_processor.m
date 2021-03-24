close all
clear
clc

i = 1;

proc(i) =   0; i = i +1;
proc(i) =   0; i = i +1;
proc(i) =   0; i = i +1;
proc(i) =   0; i = i +1;
proc(i) =   0; i = i +1;
proc(i) =   0; i = i +1;
proc(i) =   0; i = i +1;
proc(i) =   5; i = i +1;
proc(i) =  28; i = i +1;
proc(i) = 200; i = i +1;
proc(i) = 252; i = i +1;
proc(i) = 239; i = i +1;
proc(i) = 294; i = i +1;
proc(i) = 253; i = i +1;
proc(i) = 350; i = i +1;
proc(i) = 398; i = i +1;
proc(i) = 495; i = i +1;
proc(i) = 511; i = i +1;
proc(i) = 492; i = i +1;
proc(i) = 549; i = i +1;
proc(i) = 563; i = i +1;
proc(i) = 629; i = i +1;
proc(i) = 654; i = i +1;
proc(i) = 830; 


mean_per_proc = sum(proc)/(length(proc)-0);

mean_per_proc_vec = mean_per_proc*ones(1,26);

bar(proc)

hold on
plot(0:25,mean_per_proc_vec,'r','LineWidth',2.5)
xlim([0 24.6])

fontsize = 14;

xlabel('Processor Number','FontSize',fontsize,'FontName','TimesNewRoman')
ylabel('Numer of Convective Points','FontSize',fontsize,'FontName','TimesNewRoman')
title('Total Number of Convective Points = 6742','FontSize',fontsize,'FontName','TimesNewRoman')
legend('Number Per Processor','Mean','Location','NorthWest')

set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

saveas(gcf,'processor_dist.eps','psc2')