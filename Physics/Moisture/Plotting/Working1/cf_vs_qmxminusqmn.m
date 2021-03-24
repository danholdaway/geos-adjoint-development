close all
clear
clc
% addpath(cd)

cd /discover/nobackup/drholdaw/tmp.22292/sens.20130117.000000/

addpath(cd)

i = 1;
cfbt_vs_qmxmqmn
cfoba = cfob;
% qmxmqmna = func;
qmxmqmna = qmxmqmn;


i = 1;
cfbt_vs_qsmqmn1
cfobb = cfob;
% qsmqmnb = func;
qsmqmnb = qsmqmn;

i = 1;
cfbt_vs_qmxmqs
cfobc = cfob;
% qmxmqsc = func;
qmxmqsc = qmxmqs;

% i = 1;
% cfbt_vs_qcob
% cfobd = cfob;
% qcob = func;

clear qmxmqs qsmqmn qmxmqmn cfob

ii = 1;
jj = 1;
kk = 1;
cfbt_vs_q_problem1


cd /home/drholdaw/Lin_Moist_Physics/DAS_plotting/Cloud' Fraction'/

figure
scatter(qmxmqmna,cfoba,'x')
hold on
scatter(qmxmqmn,cfob,'rx')

figure
scatter(qmxmqsc,cfobc,'x')
hold on
scatter(qmxmqs,cfob,'rx')

figure
scatter(qsmqmnb,cfobb,'x')
hold on
scatter(qsmqmn,cfob,'rx')

qratiob = 100*qsmqmnb./qmxmqmna;
qratio = 100*qsmqmn./qmxmqmn;
figure
scatter(qratiob,cfobb,'x')
hold on
scatter(qratio,cfob,'rx')
scatter(qratio(end),cfob(end),'gx')

asd


[cfob1_sort, j] = sort(cfoba,'ascend');
qmxmqmn_sort = qmxmqmn(j);

[cfob2_sort, j] = sort(cfobb,'ascend');
qmxmqs_sort = qmxmqs(j);

[cfob3_sort, j] = sort(cfobc,'ascend');
qsmqmn_sort = qsmqmn(j);

[cfob4_sort, j] = sort(cfobd,'ascend');
qcob_sort = qcob(j);

figure
plot((qmxmqmn_sort),abs(cfob1_sort),'b')

figure
plot((qmxmqs_sort),abs(cfob2_sort),'b')

figure
plot(qsmqmn_sort,abs(cfob3_sort),'b')

figure
plot(100*qsmqmn./qmxmqmn,abs(cfob3_sort),'b')

figure
plot(abs(qcob_sort),abs(cfob4_sort),'b')

figure
plot(abs(qcob_sort.*qmxmqs_sort),abs(cfob4_sort),'-x')



%SHOW PERCENTAGES

qmxmqmn_sort1 = sort(qmxmqmn,'ascend');
qmxmqs_sort1 = sort(qmxmqs,'ascend');
qsmqmn_sort1 = sort(qsmqmn,'ascend');


i90 = round(90*i/100)
i95 = round(95*i/100)
i975 = round(97.5*i/100)

qmxmqmn_90 = qmxmqmn_sort1(i90);
qmxmqmn_95 = qmxmqmn_sort1(i95);
qmxmqmn_975 = qmxmqmn_sort1(i975);

qmxmqs_90 = qmxmqs_sort1(i90);
qmxmqs_95 = qmxmqs_sort1(i95);
qmxmqs_975 = qmxmqs_sort1(i975);

qsmqmn_90 = qsmqmn_sort1(i90);
qsmqmn_95 = qsmqmn_sort1(i95);
qsmqmn_975 = qsmqmn_sort1(i975);

k = 100;
% cf_vec = [min(cfob1_sort) : (max(cfob1_sort) - min(cfob1_sort))/(k-1) : max(cfob1_sort)];
cf_vec = [0 : (max(abs(cfob1_sort)) - 0)/(k-1) : max(abs(cfob1_sort))];


figure(1)
hold on
plot(qmxmqmn_90*ones(1,k),cf_vec,'r')
plot(qmxmqmn_95*ones(1,k),cf_vec,'g')
plot(qmxmqmn_975*ones(1,k),cf_vec,'k')

figure(2)
hold on
plot(qmxmqs_90*ones(1,k),cf_vec,'r')
plot(qmxmqs_95*ones(1,k),cf_vec,'g')
plot(qmxmqs_975*ones(1,k),cf_vec,'k')

figure(3)
hold on
plot(qsmqmn_90*ones(1,k),cf_vec,'r')
plot(qsmqmn_95*ones(1,k),cf_vec,'g')
plot(qsmqmn_975*ones(1,k),cf_vec,'k')


