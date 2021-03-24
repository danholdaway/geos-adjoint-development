close all
clear
clc

load qmnqmx.mat

q_mx = qmx(1); 
q_mn = qmn(1);

% qsmqmx = q_mn:(q_mx-q_mn)/100:q_mx;
qsmqmx = 0:0.01:5;


% f = ones(1,length(qsmqmx));

a = 5

% f = 0.5*(1-tanh(a*qsmqmx-3));

f = sech(100*qsmqmx);

plot(qsmqmx,f)