close all
clear
clc

N = 64;
dx = 1/64;
C = 0.1;
dt = C*dx;

x = 0:dx:1;


q_1 = 0.5*(1.0+sin(2.0*pi*x));

q_2 = zeros(1,length(x));
q_3 = zeros(1,length(x));

for i = 1:length(x)
   
    if x(i) > 0.25 && x(i) < 0.75
        q_2(i) = 1;
    end
    
end

q_3(floor(N/2))=1.0;

fontsize = 13;
lin_wid = 1.5;
grey1 = 0.75;

plot(x,q_1,'b','LineWidth',lin_wid)
hold on
plot(x,q_2,'r','LineWidth',lin_wid)
plot(x,q_3,'g','LineWidth',lin_wid)
plot(x,q_2,'r','LineWidth',lin_wid)

ylim([-0.01 1.01])
xlabel('x','FontSize',fontsize,'FontName','TimesNewRoman')
ylabel('q','FontSize',fontsize,'FontName','TimesNewRoman')
title('Initial profiles','FontSize',fontsize,'FontName','TimesNewRoman')
legend('Test Case 1','Test Case 2','Test Case 3','location','West')
% legend boxoff
box on
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')