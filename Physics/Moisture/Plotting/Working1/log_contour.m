close all
clear
clc


x = linspace(0,3,10000);
y = exp(-x).*(sin(10*x)-0.5);
plot(x,y)
figure
Y = sign(y).*log10(abs(y));
plot(x,Y)
yl = get(gca,'ytick')
% set(gca,'yticklabel',sign(yl).*10.^abs(yl))