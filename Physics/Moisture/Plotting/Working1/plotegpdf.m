close all
clear
clc
mydir = pwd;


qt1    = 0.001610473116459;
qstar1 = 0.001612508437677;
alpha1 = 0.010557271308417;
sigma1 = alpha1*qstar1;
a1 = qt1+sigma1-qstar1;
cf = (qt1+sigma1-qstar1)/(2*sigma1)

y = 0:0.1:1;

i = length(y);

qt    = qt1*ones(1,i);
qstar = qstar1*ones(1,i);
alpha = alpha1*ones(1,i);
sigma = sigma1*ones(1,i);
a = a1*ones(1,i);


plot(qt,y,'k')
hold on
% plot(qt+alpha,y,'r')
% plot(qt-alpha,y,'r')
plot(qt+sigma,y,'r')
plot(qt-sigma,y,'r')
% plot(a,y,'g--')
