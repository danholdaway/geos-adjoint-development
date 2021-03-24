close all
clear
clc
% addpath(cd)

cd /discover/nobackup/drholdaw/tmp.22292/sens.20130117.000000/

addpath(cd)

i = 1;
j = 1;
cfd_vs_qmxqsqmn


cd /home/drholdaw/Lin_Moist_Physics/DAS_plotting/Cloud' Fraction'/

qsqmn = qs - qmn;
qmxqmn = qmx - qmn;
qmxqs = qmx - qs;
qratio = 100*qsqmn./qmxqmn;


figure
plot(qs,abs(cfod),'x')

figure
plot(qmn,abs(cfod),'x')

figure
plot(qmx,abs(cfod),'x')

figure
plot(qsqmn,abs(cfod),'x')

figure
plot(qmxqs,abs(cfod),'x')

figure
plot(qmxqmn,abs(cfod),'x')

figure
plot(qratio,abs(cfod),'x')



% figure
% plot(qmxmqmn,abs(cfod),'x')
% 
% figure
% plot(qmxqs,abs(cfod),'x')
% 
% figure
% semilogy(qsqmn,abs(cfod),'x')
% 
% figure
% plot(qratio,abs(cfod),'x')
% 
% % figure
% % plot(dq,abs(cfod),'x')
% 
% figure
% plot(qmxdqsd,abs(cfod),'x')
% 
% figure
% plot(qmxd,abs(cfod),'x')
% 
% figure
% plot(qsd,abs(cfod),'x')
% 
% figure
% plot(dqd,abs(cfod),'x')
% 
% figure
% plot(qtd,abs(cfod),'x')
% 
% figure
% plot(qsd,abs(cfod),'x')

