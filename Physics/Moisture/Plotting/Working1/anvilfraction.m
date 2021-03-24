close all
clear
clc

AF = 0:0.0001:1;


tempAF = 1./(1-AF);

dTdAF_new = (2./-((-1)./(1.-AF).^1.75)) .* -((-1)./(1.-AF).^2);
dTdAF = -((-1)./(1.-AF).^2);

dTdAF1 = exp(AF.^100000).^1000;

dTdAF_new = (1./dTdAF1).* -((-1)./(1.-AF).^2);

plot(AF,dTdAF)


figure
plot(AF,dTdAF_new,'r--')


% ylim([0 max(dTdAF_new)])


figure
plot(AF,dTdAF1)