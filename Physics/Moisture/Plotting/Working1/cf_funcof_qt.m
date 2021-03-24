close all
clear
clc


qmin = 0.1e-3;
qmax = 0.9e-3;

alpha = 0.06;

qt = 0.8e-3;

qs = 0:1e-5:1e-3;

dq = 2*alpha*qs;

Cf = (qt + alpha*qs - qs)./dq;

plot(qs,Cf)