close all
clear
clc



i = 1;
x(i) = datenum(2014,04,20,12,00,00);   a(i) = 177.2611   ; b(i) = 89.1642; i = i+1;
x(i) = datenum(2014,04,20,18,00,00);   a(i) = 178.7100   ; b(i) = 92.5360; i = i+1;
x(i) = datenum(2014,04,21,00,00,00);   a(i) = 180.9398   ; b(i) = 94.1558; i = i+1;
x(i) = datenum(2014,04,21,06,00,00);   a(i) = 193.2705   ; b(i) = 92.4198; i = i+1;
x(i) = datenum(2014,04,21,12,00,00);   a(i) = 178.3461   ; b(i) = 91.1055; i = i+1;
x(i) = datenum(2014,04,21,18,00,00);   a(i) = 186.7132   ; b(i) = 94.4299; i = i+1;
x(i) = datenum(2014,04,22,00,00,00);   a(i) = 201.0412   ; b(i) = 93.2707; i = i+1;
x(i) = datenum(2014,04,22,06,00,00);   a(i) = 193.6542   ; b(i) = 91.8489; i = i+1;
x(i) = datenum(2014,04,22,12,00,00);   a(i) = 198.7870   ; b(i) = 93.9335; i = i+1;
x(i) = datenum(2014,04,22,18,00,00);   a(i) = 201.1398   ; b(i) = 92.4305; i = i+1;
x(i) = datenum(2014,04,23,00,00,00);   a(i) = 179.4767   ; b(i) = 92.4594; i = i+1;
x(i) = datenum(2014,04,23,06,00,00);   a(i) = 181.8802   ; b(i) = 92.0806; i = i+1;
x(i) = datenum(2014,04,23,12,00,00);   a(i) = 187.5875   ; b(i) = 92.9769; i = i+1;
x(i) = datenum(2014,04,23,18,00,00);   a(i) = 205.7498   ; b(i) = 93.4152; i = i+1;
x(i) = datenum(2014,04,24,00,00,00);   a(i) = 205.1510   ; b(i) = 95.1507; i = i+1;
x(i) = datenum(2014,04,24,06,00,00);   a(i) = 188.5005   ; b(i) = 92.3858; i = i+1;
x(i) = datenum(2014,04,24,12,00,00);   a(i) = 186.2454   ; b(i) = 92.2839; i = i+1;
x(i) = datenum(2014,04,24,18,00,00);   a(i) = 183.4356   ; b(i) = 92.2485; i = i+1;
x(i) = datenum(2014,04,25,00,00,00);   a(i) = 183.6512   ; b(i) = 91.8220; i = i+1;
x(i) = datenum(2014,04,25,06,00,00);   a(i) = 182.1001   ; b(i) = 90.7297; i = i+1;
x(i) = datenum(2014,04,25,12,00,00);   a(i) = 186.4037   ; b(i) = 91.9179; i = i+1;
% x(i) = datenum(2014,04,25,18,00,00);   a(i) = 195.9719   ; b(i) = ; i = i+1;
% x(i) = datenum(2014,04,26,00,00,00);   a(i) = 176.1509   ; b(i) = ; i = i+1;



fontsize = 13;
lin_wid = 1.5;
j = i-1;

figure
set(gcf,'position',[105 327 1097 542])
plot(x(1:j),a(1:j),'b-x','LineWidth',lin_wid)
hold on
plot(x(1:j),b(1:j),'r-x','LineWidth',lin_wid)
box on
datetick('x','dd-mm-yy')
ylim([0 220])
xlim([x(1) x(end)])
xlabel('Assimilation Time (2014)','FontSize',fontsize,'FontName','TimesNewRoman')
ylabel('Moist physics timing (seconds)','FontSize',fontsize,'FontName','TimesNewRoman')
title('Total time for linearized moist physics during minimization','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
legend('Old Moist','New 4DVAR Moist','Location','NorthWest')

legend boxoff
