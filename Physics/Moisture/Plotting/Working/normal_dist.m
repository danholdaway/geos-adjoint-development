close all
clear
clc

x = 0:0.01:1;

sigma = 0.01;
mu = 0.5;
a = 4;

y = exp(-(x-mu).^a/(2*sigma^2));

dydx = (-a*(x-mu).^(a-1).*exp(-(x-mu).^a/(2*sigma^2)))/(2*sigma^2);


plot(x,y)

figure
plot(x,dydx)