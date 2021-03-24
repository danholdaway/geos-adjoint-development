close all
clear
clc

qmn = 0.1;
qmx = 0.9;
qs = qmn:0.01:qmx;

cf = (qmx - qs)./(qmx-qmn);

plot(qs,cf)


qs = qmn:0.001:qmx;
cf1 = 0.5*(1+cos(pi*(qs-qmn)/(qmx-qmn)));

hold on
plot(qs,cf1,'r')

